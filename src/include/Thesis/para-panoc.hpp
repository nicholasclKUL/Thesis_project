#pragma once

#include <alpaqa/accelerators/lbfgs.hpp>
#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/inner/directions/panoc/lbfgs.hpp>
#include <alpaqa/implementation/inner/panoc.tpp>
#include <alpaqa/implementation/inner/panoc-ocp.tpp>
#include <alpaqa/implementation/inner/panoc-helpers.tpp>
#include <Kokkos_Core.hpp>

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace alpaqa {


template <Config Conf>
class ParaPANOCSolver {
  public:
    USING_ALPAQA_CONFIG(Conf);

    using Problem      = alpaqa::TypeErasedControlProblem<config_t>;
    using Params       = PANOCOCPParams<config_t>;
    using Stats        = PANOCOCPStats<config_t>;
    using ProgressInfo = PANOCOCPProgressInfo<config_t>;
    using SolveOptions = InnerSolveOptions<config_t>;

    ParaPANOCSolver(const Params &params) : params(params) {}

    Stats operator()(const Problem &problem,    // in
                     const SolveOptions &opts,  // in
                     rvec xu_init,              // inout
                     rvec y,                    // in
                     crvec μ,                   // in
                     rvec err_z,                // out
                     rvec g,                    // out
                     index_t nthrds);
    template <class P>
    Stats operator()(const P &problem, const SolveOptions &opts, rvec xu_init, rvec y,
                     crvec μ, rvec e, rvec g, index_t nthrds) {
        return operator()(Problem::template make<P>(problem), opts, xu_init, y, μ, e, g, nthrds);
    }

    /// Specify a callable that is invoked with some intermediate results on
    /// each iteration of the algorithm.
    /// @see @ref ProgressInfo
    ParaPANOCSolver &
    set_progress_callback(std::function<void(const ProgressInfo &)> cb) {
        this->progress_cb = cb;
        return *this;
    }

    std::string get_name() const;

    void stop() { stop_signal.stop(); }

    const Params &get_params() const { return params; }

  private:
    Params params;
    AtomicStopSignal stop_signal;
    std::function<void(const ProgressInfo &)> progress_cb;
    using Helpers = detail::PANOCHelpers<config_t>;

  public:
    std::ostream *os = &std::cout;
    
};



    // Implementation --------------------------------------------------------



template <Config Conf>
auto ParaPANOCSolver<Conf>::operator()(
    /// [in]    Problem description
    const Problem &problem,
    /// [in]    Solve options
    const SolveOptions &opts,
    /// [inout] Initial vector input for MS-OCP
    rvec xu,
    /// [inout] Lagrange multipliers @f$ y @f$
    rvec y,
    /// [in]    Penalty factors @f$ \mu @f$
    crvec μ,
    /// [out]   Slack variable error @f$ c(x) - \Pi_D(c(x) + \mu^{-1} y) @f$
    rvec err_z, 
    /// [out]   Continuity, f(x,u) - x+
    rvec g,
    /// [in]    Number of threads
    index_t nthrds) -> Stats {

    using std::chrono::nanoseconds;
    auto os         = opts.os ? opts.os : this->os;
    auto start_time = std::chrono::steady_clock::now();
    Stats s;
    LBFGS<config_t> lbfgs{params.lbfgs_params};

    const auto N    = problem.get_N()+1;      // horizon N
    const auto nu   = problem.get_nu();     // number of inputs
    const auto nx   = problem.get_nx();     // number of states
    const auto nc   = problem.get_nc();     // number of box constraints
    const auto nc_N = problem.get_nc_N();   // box constraints * horizon
    const auto nxu  = nx + nu;              // number of decision variables per stage
    const auto n    = (nxu*(N-1)) + nx;     // total number of decision variables
    const auto m    = nx*(N-1);             // total number of dynamic constraints

    // Allocate storage --------------------------------------------------------

    vec q(n); // Newton step, including states
    Box<config_t> U   = Box<config_t>::NaN(nu);   
    Box<config_t> D   = Box<config_t>::NaN(nx);
    Box<config_t> D_N = Box<config_t>::NaN(nx);
    Box<config_t> F   = Box<config_t>::NaN(nx); 
    Box<config_t> Dfd = Box<config_t>::NaN(m); 

    F.lowerbound.setConstant(0);
    F.upperbound.setConstant(0);

    // Iterates ----------------------------------------------------------------
 
    // Represents an iterate in the algorithm, keeping track of some
    // intermediate values and function evaluations.
    struct Iterate {
        vec xû;       //< Inputs u interleaved with states x after prox grad
        vec xu;       //< Inputs u interleaved with states x -> x₀, u₀, x₁, u₁,...,uₙ₋₁, xₙ 
        vec grad_ψ;   //< Gradient of cost w.r.t {x₀, u₀, x₁, u₁,...,uₙ₋₁, xₙ}
        vec fxu;      //< Dynamics f(x,u)
        vec hxu;      //< nonconvex mapping h(x,u)
        vec fxû;      //< Dynamics f(x,û)
        vec hxû;      //< nonconvex mapping h(x,û)
        vec lxu;      //< cost function at iterate -> l(x,u)
        vec lxû;      //< cost function after T(x,u) -> l(x,û)
        mat Jfxu;     //< Jacobian of dynamics Jf(x,u)
        vec qr;       //< Cost function gradient Jh(x,u)dl(x,u)
        vec g;        //< Dynamic constraints, f(x,u) - x+
        vec gz;       //< g(x) + Σ⁻¹y
        vec gz_hat;   //< g(x_hat) + Σ⁻¹y
        vec Π_D;      //< Projection of states and variables into their box constraints
        vec Π_D_hat;  //< Projection of states and variables into their box constraints
        vec Z;        //< (Σ*(f(xₖ,uₖ)-xₖ₊₁-Π_D(f(xₖ,uₖ)-xₖ₊₁+Σ⁻¹y))-y)
        vec p;        //< Proximal gradient step
        real_t ψxu      = NaN<config_t>;            //< Cost in x
        real_t ψxû      = NaN<config_t>;            //< Cost in x̂
        real_t γ        = NaN<config_t>;            //< Step size γ
        real_t L        = NaN<config_t>;            //< Lipschitz estimate L
        real_t pᵀp      = NaN<config_t>;            //< Norm squared of p
        real_t grad_ψᵀp = NaN<config_t>;            //< Dot product of gradient and p

        // @pre    @ref ψxu, @ref pᵀp, @pre grad_ψᵀp
        // @return φγ
        real_t fbe() const { return ψxu + pᵀp / (2 * γ) + grad_ψᵀp; }

        Iterate(length_t n, length_t m, length_t nxu, length_t N) :
            xu(n), xû(n), grad_ψ(n), fxu(m), fxû(m), hxu(n), hxû(n), 
            qr(n), g(m), gz(m), gz_hat(m), Π_D(m), Π_D_hat(m), Z(m), p(n), Jfxu(m, nxu), 
            lxu(N), lxû(N) {}
    
    } 
    
    iterates[2]{{n, m, nxu, N}, {n, m, nxu, N}};     
    Iterate *curr = &iterates[0];
    Iterate *next = &iterates[1];
    
    // Helper functions --------------------------------------------------------

    auto eval_proj_set = [&](const Box<config_t> &box, crvec x) {
        using binary_real_f = real_t (*)(real_t, real_t);
        return x.binaryExpr(box.lowerbound, binary_real_f(std::fmax))
                .binaryExpr(box.upperbound, binary_real_f(std::fmin));
    };

    auto eval_iterate = [&](int k, Iterate &It, crvec μ, crvec y){       

        Eigen::Index k_ = k;

        if (k_ == N-1){
            problem.eval_h_N(It.xu.segment(k*nxu,nx),
                            It.hxu.segment(k*nxu,nx));
            It.lxu(k) = problem.eval_l_N(It.hxu.segment(k*nxu,nx));
            problem.eval_q_N(It.xu.segment(k*nxu,nx),
                            It.hxu.segment(k*nxu,nx), 
                            It.qr.segment(k*nxu,nx));
        } else {
            problem.eval_f(k, It.xu.segment(k*nxu,nx),
                            It.xu.segment((k*nxu)+nx,nu), 
                            It.fxu.segment(k*nx,nx));
            problem.eval_jac_f(k, It.xu.segment(k*nxu,nx),
                            It.xu.segment((k*nxu)+nx,nu),
                            It.Jfxu.block(k*nx,0,nx,nxu));
            problem.eval_h(k, It.xu.segment(k*nxu,nx), 
                            It.xu.segment((k*nxu)+nx,nu),
                            It.hxu.segment(k*nxu,nxu));            
            It.lxu(k) = problem.eval_l(k, It.hxu.segment(k*nxu,nxu));
            problem.eval_qr(k, It.xu.segment(k*nxu,nxu), 
                            It.hxu.segment(k*nxu,nxu),
                            It.qr.segment(k*nxu,nxu));
            It.g.segment(k*nx,nx) = It.fxu.segment(k*nx,nx) - It.xu.segment((k+1)*nxu,nx);           
            It.gz.segment(k*nx,nx) = It.g.segment(k*nx,nx) + 
                                    ((μ).segment(k*nx,nx).asDiagonal().inverse()*(y.segment(k*nx,nx)));            
            It.Π_D.segment(k*nx,nx).noalias() = eval_proj_set(F, It.gz.segment(k*nx,nx));
            It.Z.segment(k*nx,nx)   = (μ.segment(k*nx,nx).cwiseProduct(It.fxu.segment(k*nx,nx)
                                    -It.xu.segment((k+1)*nxu,nx)-It.Π_D.segment(k*nx,nx)) 
                                    + y.segment(k*nx,nx));          
        }
    };

    auto eval_grad_ψ_k = [&](int k, Iterate &It, crvec μ, crvec y){

        Eigen::Index k_ = k;

        mat I; I.setIdentity(nxu,nx);

        if (k_ == N-1){
            It.grad_ψ.segment(k*nxu,nx) = It.qr.segment(k*nxu,nx) - It.Z.segment((k-1)*nx,nx);                            
        } 
        else if(k_ == 0){ //this part need optimization, unnecessary calculations done!
            It.grad_ψ.segment(k*nxu,nxu) = It.qr.segment(k*nxu,nxu) +
                        It.Jfxu.block(k*nx,0,nx,nxu).transpose() * It.Z.segment(k*nx,nx);
            It.grad_ψ.segment(0,nx).setZero();
        }
        else {
            It.grad_ψ.segment(k*nxu,nxu) = It.qr.segment(k*nxu,nxu) - I*It.Z.segment((k-1)*nx,nx) +
                        It.Jfxu.block(k*nx,0,nx,nxu).transpose() * It.Z.segment(k*nx,nx);
        }
    };

    auto eval_ψ_k = [&](int k, Iterate &It, crvec μ){

        Eigen::Index k_ = k;
        real_t ψxu_k = 0;

        if (k_ == N-1){
            return ψxu_k = It.lxu(k);
        }
        else {
            auto d = (It.gz.segment(k*nx,nx)-
                      eval_proj_set(F,It.gz.segment(k*nx,nx)));
            return ψxu_k = It.lxu(k) + 
                         0.5*d.transpose()*(μ).segment(k*nx,nx).asDiagonal()*d;
        }
    };

    auto eval_ψ_hat_k = [&](int k, Iterate &It, crvec μ, crvec y){
        
        real_t ψxû_k = 0;
        Eigen::Index i_ = k;

        if (i_ == N-1){
            problem.eval_h_N(It.xû.segment(k*nxu,nx),
                             It.hxû.segment(k*nxu,nx));
            It.lxû(k) = problem.eval_l_N(It.hxû.segment(k*nxu,nx));
            return ψxû_k = It.lxû(k);
        } 
        else {
            problem.eval_f(k, It.xû.segment(k*nxu,nx),
                          It.xû.segment(k*(nxu)+nx,nu), 
                          It.fxû.segment(k*nx,nx));
            problem.eval_h(k, It.xû.segment(k*nxu,nx), 
                            It.xû.segment((k*nxu)+nx,nu),
                            It.hxû.segment(k*nxu,nxu));
            It.lxû(k) = problem.eval_l(k, It.hxû.segment(k*nxu,nxu));                       
            It.gz_hat.segment(k*nx,nx) = It.fxû.segment(k*nx,nx) - 
                            It.xû.segment((k+1)*nxu,nx)  
                            + ((μ).segment(k*nx,nx).asDiagonal().inverse()*(y.segment(k*nx,nx)));
            It.Π_D_hat.segment(k*nx,nx) = eval_proj_set(F, It.gz_hat.segment(k*nx,nx));  
            auto d = (It.gz_hat.segment(k*nx,nx)-
                      eval_proj_set(F,It.gz_hat.segment(k*nx,nx))); 
            return ψxû_k = It.lxû(k) + 
                         0.5*d.transpose()*(μ).segment(k*nx,nx).asDiagonal()*d;            
        }
    };
    
    // Τγₖ(xₖ)
    auto eval_proj_grad_step_box = [&](const Box<config_t> &box, real_t γ, crvec x, crvec grad_ψ, rvec x̂,
                                        rvec p) {
        using binary_real_f = real_t (*)(real_t, real_t);
        p                   = (-γ * grad_ψ)
                .binaryExpr(box.lowerbound - x, binary_real_f(std::fmax))
                .binaryExpr(box.upperbound - x, binary_real_f(std::fmin));
        x̂ = x + p;
    };
    
    auto eval_prox = [&](int k, Iterate &It) {
        
        Eigen::Index k_ = k;

        if (k_==N-1){
            eval_proj_grad_step_box(D_N, It.γ, It.xu.segment(k*nxu,nx), 
                                    It.grad_ψ.segment(k*nxu,nx), 
                                    It.xû.segment(k*nxu,nx), It.p.segment(k*nxu,nx));
        }
        else{
            eval_proj_grad_step_box(D, It.γ, It.xu.segment(k*nxu,nx), 
                                    It.grad_ψ.segment(k*nxu,nx), 
                                    It.xû.segment(k*nxu,nx), It.p.segment(k*nxu,nx));
            eval_proj_grad_step_box(U, It.γ, It.xu.segment((k*nxu)+nx,nu), 
                                    It.grad_ψ.segment((k*nxu)+nx,nu), 
                                    It.xû.segment((k*nxu)+nx,nu), It.p.segment((k*nxu)+nx,nu));
        }
    };

    auto eval_pᵀp_grad_ψᵀp = [&](Iterate &It) {
        It.pᵀp      = It.p.squaredNorm();
        It.grad_ψᵀp = It.grad_ψ.dot(It.p);
    };

    auto fd_grad_ψ = [&](Iterate &it, crvec μ, crvec y) {
        
        alpaqa::ScopedMallocAllower ma;
        
        real_t h = 0.001;
        real_t ψxu_ = 0;
        index_t N   = problem.get_N(); 
        index_t nu  = problem.get_nu();
        index_t nx  = problem.get_nx();
        index_t nxu = nx + nu;
        index_t n   = (N-1)*(nxu) + nx;
        index_t m   = (N-1)*nx;
        vec lxu_    = vec::Zero(N);
        vec hxu_    = vec::Zero(n);
        vec xu_     = vec::Zero(n); 
        vec fxu_    = vec::Zero(m);
        vec d       = vec::Zero(m);
        vec v       = vec::Zero(m);
        vec g       = vec::Zero(m);
        Box D       = alpaqa::Box<config_t>::NaN(m);
        problem.get_D(D);

        for (index_t i = 0; i < n; ++i) {
            xu_ = it.xu;
            xu_(i) += h;
            for (index_t k = 0; k < N; ++k) {
                if (k == N-1) {
                    problem.eval_h_N(xu_.segment(k*nxu,nx),
                                        hxu_.segment(k*nxu,nx));
                    lxu_(k) = problem.eval_l_N(hxu_.segment(k*nxu,nx));
                } else {
                    problem.eval_f(k, xu_.segment(k*nxu,nx),
                                    xu_.segment(k*(nxu)+nx,nu), 
                                    fxu_.segment(k*nx,nx));
                    problem.eval_h(k, xu_.segment(k*nxu,nx), 
                                    xu_.segment((k*nxu)+nx,nu),
                                    hxu_.segment(k*nxu,nxu));
                    lxu_(k) = problem.eval_l(k, hxu_.segment(k*nxu,nxu)); 
                }
            }
            for (index_t k = 0; k < N-1; ++k) {
                g.segment(k*nx,nx) = fxu_.segment(k*nx,nx) - 
                                        xu_.segment((k+1)*nxu,nx);
            }
            d.noalias() = g + μ.asDiagonal().inverse()*y; 
            v.noalias() = d - eval_proj_set(F, d); 
            ψxu_ = lxu_.sum() + real_t(.5)*v.transpose()*(μ).asDiagonal()*v;
            it.grad_ψ(i) = (ψxu_ - it.ψxu)/h;
        }
    };

    auto calc_error_stop_crit = [this](
                                    real_t γ, crvec xuₖ, crvec grad_ψₖ,
                                    crvec pₖ, real_t pₖᵀpₖ, rvec work_xu,
                                    rvec work_p) {
        switch (params.stop_crit) {
            case PANOCStopCrit::ProjGradNorm: {
                return vec_util::norm_inf(pₖ);
            }
            case PANOCStopCrit::ProjGradNorm2: {
                return std::sqrt(pₖᵀpₖ);
            }
            case PANOCStopCrit::FPRNorm: {
                return vec_util::norm_inf(pₖ) / γ;
            }
            case PANOCStopCrit::FPRNorm2: {
                return std::sqrt(pₖᵀpₖ) / γ;
            }
            case PANOCStopCrit::ApproxKKT: [[fallthrough]];
            case PANOCStopCrit::ApproxKKT2: [[fallthrough]];
            case PANOCStopCrit::Ipopt: [[fallthrough]];
            case PANOCStopCrit::LBFGSBpp: [[fallthrough]];
            default:
                throw std::invalid_argument("Unsupported stopping criterion");
        }
    };

    auto check_all_stop_conditions =
        [this, &opts](
            /// [in]    Time elapsed since the start of the algorithm
            auto time_elapsed,
            /// [in]    The current iteration number
            unsigned iteration,
            /// [in]    Tolerance of the current iterate
            real_t εₖ,
            /// [in]    The number of successive iterations no progress was made
            unsigned no_progress) {
            auto max_time = params.max_time;
            if (opts.max_time)
                max_time = std::min(max_time, *opts.max_time);
            auto tolerance = opts.tolerance > 0 ? opts.tolerance : real_t(1e-8);
            bool out_of_time     = time_elapsed > max_time;
            bool out_of_iter     = iteration == params.max_iter;
            bool interrupted     = stop_signal.stop_requested();
            bool not_finite      = not std::isfinite(εₖ);
            bool conv            = εₖ <= tolerance;
            bool max_no_progress = no_progress > params.max_no_progress;
            return conv              ? SolverStatus::Converged
                   : out_of_time     ? SolverStatus::MaxTime
                   : out_of_iter     ? SolverStatus::MaxIter
                   : not_finite      ? SolverStatus::NotFinite
                   : max_no_progress ? SolverStatus::NoProgress
                   : interrupted     ? SolverStatus::Interrupted
                                     : SolverStatus::Busy;
        };

    auto qub_violated = [this](const Iterate &i) {
        real_t margin =
            (1 + std::abs(i.ψxu)) * params.quadratic_upperbound_tolerance_factor;
        return i.ψxû > i.ψxu + i.grad_ψᵀp + real_t(0.5) * i.L * i.pᵀp + margin;
    };

    auto linesearch_violated = [this](const Iterate &curr,
                                      const Iterate &next) {
        real_t σ  = params.β * (1 - curr.γ * curr.L) / (2 * curr.γ);
        real_t φγ = curr.fbe();
        real_t margin = (1 + std::abs(φγ)) * params.linesearch_tolerance_factor;
        return next.fbe() > φγ - σ * curr.pᵀp + margin;
    };

    auto initial_lipschitz_estimate =
        [&](
            /// Iterate, updates xu, ψ, grad_ψ, have_jacobians, L
            Iterate *cur,
            /// Iterate
            Iterate *nex,
            /// [in]    Finite difference step size relative to x
            real_t ε,
            /// [in]    Minimum absolute finite difference step size
            real_t δ,
            /// [in]    Minimum allowed Lipschitz estimate.
            real_t L_min,
            /// [in]    Maximum allowed Lipschitz estimate.
            real_t L_max) {

            // Calculate ψ(x₀), ∇ψ(x₀)
            Kokkos::parallel_for("iterate x₀, lipschitz", nthrds, [&](const int i){
                eval_iterate(i, *cur, μ, y);
            });
            Kokkos::fence();
            Kokkos::parallel_for("∇ψ(x₀), lipschitz", nthrds, [&](const int i){
                eval_grad_ψ_k(i, *cur, μ, y);
            });
            Kokkos::parallel_reduce("ψ(x₀), lipschitz", nthrds, [&](const int i, real_t &ψ_){
                ψ_ += eval_ψ_k(i, *cur, μ);
            },cur->ψxu);
            // Select a small step h for finite differences
            auto h        = cur->grad_ψ.unaryExpr([&](real_t g) {
                return g > 0 ? std::max(g * ε, δ) : std::min(g * ε, -δ);
            });
            real_t norm_h = h.norm();
            // work_xu = xu - h
            nex->xu = cur->xu - h;
            // Calculate ψ(x₀ - h)
            Kokkos::parallel_for("iterate (xu₀-h), lipschitz", nthrds, [&](const int i){
                eval_iterate(i, *nex, μ, y);
            });
            Kokkos::fence();
            Kokkos::parallel_for("ψ(xu₀-h), lipschitz", nthrds, [&](const int i){
                eval_grad_ψ_k(i, *nex, μ, y);
            });
            // Calculate ∇ψ(x₀ - h)
            Kokkos::parallel_reduce("∇ψ(xu₀-h), lipschitz", nthrds, [&](const int i, real_t &ψ_){
                ψ_ += eval_ψ_k(i, *nex, μ);
            },nex->ψxu);
            // Estimate Lipschitz constant using finite differences
            cur->L = (nex->grad_ψ - cur->grad_ψ).norm() / norm_h;
            cur->L = std::clamp(cur->L, L_min, L_max);

        };


    // Printing ----------------------------------------------------------------

    std::array<char, 64> print_buf;
    auto print_real = [&](real_t x) {
        return float_to_str_vw(print_buf, x, params.print_precision);
    };
    auto print_real3 = [&](real_t x) {
        return float_to_str_vw(print_buf, x, 3);
    };
    auto print_progress_1 = [&](unsigned k, real_t φₖ, real_t ψₖ, crvec grad_ψₖ,
                                real_t pₖᵀpₖ, real_t γₖ, real_t εₖ) {
        if (k == 0){
            *os << "┌─[ParaPANOC]\n";
        } else {
            *os << "├─ " << std::setw(6) << k << '\n';
            *os << "│   φγ = " << print_real(φₖ)               //
                << ",    ψ = " << print_real(ψₖ)               //
                << ", ‖∇ψ‖ = " << print_real(grad_ψₖ.norm())   //
                << ",  ‖p‖ = " << print_real(std::sqrt(pₖᵀpₖ)) //
                << ",  ‖q‖ = " << print_real(q.norm())         //
                << ",    γ = " << print_real(γₖ)               //
                << ",    ε = " << print_real(εₖ) << '\n';      //
        }
    };
    auto print_progress_2 = [&](crvec qₖ, real_t τₖ, bool did_gn, length_t nJ,
                                real_t min_rcond) {
        *os << "│  ‖q‖ = " << print_real(qₖ.norm())                       //
            << ",   #J = " << std::setw(7 + params.print_precision) << nJ //
            << ", cond = " << print_real3(real_t(1) / min_rcond)          //
            << ",    τ = " << print_real3(τₖ)                             //
            << ",    " << (did_gn ? "GN" : "L-BFGS")                      //
            << std::endl; // Flush for Python buffering
    };
    auto print_progress_3 = [&](Iterate &It, real_t &εₖ, real_t τ, unsigned k){
        if (k == 0){
            *os << "┌─[ParaPANOC]\n";
        } else {
            *os << "├───── " << k << " ──────" <<'\n';
            *os << "│    ψ = " << print_real(It.ψxu)               
                << ", ‖∇ψ‖ = " << print_real(It.grad_ψ.norm())   
                << ",  ‖p‖ = " << print_real(It.pᵀp)
                << ",  ‖q‖ = " << print_real(q.norm())
                << ",    τ = " << print_real3(τ) 
                << ",    γ = " << print_real(It.γ)               
                << ",    ε = " << print_real(εₖ) << '\n';
        }
    };
    auto print_progress_n = [&](SolverStatus status) {
        *os << "└─ " << status << " ──"
            << std::endl; // Flush for Python buffering
    };

    // Initialize inputs and initial state ----------------------------

    curr->xu = xu;   // initial state
    curr->xû = curr->xu;
    next->xu = curr->xu;
    next->xû = curr->xu;

    problem.get_U(U);               // input box constraints
    problem.get_D(D);               // general constraints
    problem.get_D_N(D_N);           // general terminal constraints

    // Make sure that we don't allocate any memory in the inner loop
    //ScopedMallocBlocker mb;

    // Estimate Lipschitz constant ---------------------------------------------

    // Finite difference approximation of ∇²ψ in starting point
    if (params.Lipschitz.L_0 <= 0) {
        initial_lipschitz_estimate(curr, next, params.Lipschitz.ε, params.Lipschitz.δ,
                                   params.L_min, params.L_max);
    }

    // Initial Lipschitz constant provided by the user
    else {
        curr->L = params.Lipschitz.L_0;
        // Calculate ψ(x₀), ∇ψ(x₀)
        Kokkos::parallel_for("iterate xu₀", nthrds, [&] (const int i){
            eval_iterate(i, *curr, μ, y);
        });
        Kokkos::fence();
        Kokkos::parallel_for("∇ψ(xu₀)", nthrds, [&] (const int i){
            eval_grad_ψ_k(i, *curr, μ, y);
        });
        Kokkos::parallel_reduce("ψ(xu₀)", nthrds, [&](const int i, real_t &ψ_){
            ψ_ += eval_ψ_k(i, *curr, μ);
        },curr->ψxu);
    }

    if (not std::isfinite(curr->L)) {
        s.status = SolverStatus::NotFinite;
        return s;
    }
    curr->γ = params.Lipschitz.Lγ_factor / curr->L;

    // First proximal gradient step:
    Kokkos::parallel_for("xû = T(xu₀)", nthrds, [&](const int i){
        eval_prox(i, *curr);
    });
    Kokkos::fence();
    eval_pᵀp_grad_ψᵀp(*curr);
    curr->ψxû = 0;
    Kokkos::parallel_reduce("ψ(xû) = T(xu₀)", nthrds, [&](const int i, real_t &ψ_){
        ψ_ += eval_ψ_hat_k(i, *curr, μ, y);
    },curr->ψxû);
    Kokkos::fence();
    // Initialize steps 
    unsigned k  = 0;
    real_t τ    = 1;

    // Keep track of how many successive iterations didn't update the iterate
    unsigned no_progress = 0;

    // Main PANOC loop
    // =========================================================================
    while (true) {

        // Check stop condition ------------------------------------------------

        real_t εₖ = calc_error_stop_crit(curr->γ, curr->xu, curr->grad_ψ,
                                         curr->p, curr->pᵀp, next->xû, next->p);

        // Print progress ------------------------------------------------------
        bool do_print =
            params.print_interval != 0 && k % params.print_interval == 0;
        if (do_print)
            print_progress_3(*curr, εₖ, τ ,k);
            //print_progress_1(k, curr->fbe(), curr->ψxu, curr->grad_ψ, curr->pᵀp,
            //                 curr->γ, εₖ);
        if (progress_cb) {
            //ScopedMallocAllower ma;
            alpaqa::detail::Timed t{s.time_progress_callback};
            progress_cb({.k             = k,
                         .xu            = curr->xu,
                         .p             = curr->p,
                         .norm_sq_p     = curr->pᵀp,
                         .x̂u            = curr->xû,
                         .φγ            = curr->fbe(),
                         .ψ             = curr->ψxu,
                         .grad_ψ        = curr->grad_ψ,
                         .ψ_hat         = curr->ψxû,
                         .q             = q,
                         .L             = curr->L,
                         .γ             = curr->γ,
                         .τ             = τ,
                         .ε             = εₖ,
                         .problem       = problem,
                         .params        = params});
        }

        // Return solution -----------------------------------------------------

        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        auto stop_status =
            check_all_stop_conditions(time_elapsed, k, εₖ, no_progress);
        if (stop_status != SolverStatus::Busy) {
            bool do_final_print = params.print_interval != 0;
            if (!do_print && do_final_print)
                print_progress_3(*curr, εₖ, τ, k);
                // print_progress_1(k, curr->fbe(), curr->ψxu, curr->grad_ψ, 
                //                  curr->pᵀp, curr->γ, εₖ);
            if (do_print || do_final_print)
                print_progress_n(stop_status);
            if (stop_status == SolverStatus::Converged ||
                stop_status == SolverStatus::Interrupted ||
                opts.always_overwrite_results) {
                    err_z = curr->g - curr->Π_D;
            }
            s.iterations   = k;
            s.ε            = εₖ;
            s.elapsed_time = duration_cast<nanoseconds>(time_elapsed);
            s.time_lqr_factor -= s.time_hessians;
            s.status   = stop_status;
            s.final_γ  = curr->γ;
            s.final_ψ  = curr->ψxû;
            s.final_h  = 0; // only box constraints
            s.final_φγ = curr->fbe();
            return s;
        }

        // Calculate quasi-Newton step -----------------------------------------

        real_t τ_init = NaN<config_t>;
        if (k == 0) { // Initialize L-BFGS
            lbfgs.resize(n);
            q = curr->p;
            τ_init = 0;
        }
        if (k > 0) {
            τ_init = 1;
            lbfgs.apply(q, curr->γ);
        }
        // Make sure quasi-Newton step is valid
        if (not q.allFinite()) {
            if (τ_init == 1) { // If we computed a quasi-Newton step
                ++s.lbfgs_failures;
                lbfgs.reset(); // Is there anything else we can do?
            }
            τ_init = 0;
        }
        // Line search ---------------------------------------------------------

        next->γ       = curr->γ;
        next->L       = curr->L;
        τ             = τ_init;
        real_t τ_prev = -1;

        // xₖ₊₁ = xₖ + pₖ
        auto take_safe_step = [&] {
            next->xu = curr->xû; // makes xₖ₊₁ = xₖ + pₖ = curr->xû, 
            next->ψxu = curr->ψxû; // and consequentely, ψ(xₖ₊₁) = ψ(xₖ) = curr->ψxû
            // Calculate ∇ψ(xₖ₊₁) .: ∇ψ(xₖ₊₁) = ∇ψ(xₖ); shouldn't we reuse ∇ψ(xₖ) from previous step??  
            Kokkos::parallel_for("xₖ₊₁, safe step", nthrds, [&](const int i){
                eval_iterate(i, *next, μ, y);
            });
            Kokkos::fence();
            Kokkos::parallel_for("∇ψ(xₖ₊₁), safe step", nthrds, [&] (const int i){
                eval_grad_ψ_k(i, *next, μ, y);
            });
        };

        // xₖ₊₁ = xₖ + (1-τ) pₖ + τ qₖ
        auto take_accelerated_step = [&](real_t τ) {
            if (τ == 1) { // → faster quasi-Newton step
                next->xu = curr->xu + q;
            }
            else {
                next->xu = curr->xu + (1-τ)*curr->p + τ*q;
            }
            // Calculate ψ(xₖ₊₁), ∇ψ(xₖ₊₁)
            Kokkos::parallel_for("xₖ₊₁, accelerated step", nthrds, [&](const int i){
                eval_iterate(i, *next, μ, y);
            });
            Kokkos::fence();
            Kokkos::parallel_for("∇ψ(xₖ₊₁), accelerated step", nthrds, [&] (const int i){
                eval_grad_ψ_k(i, *next, μ, y);
            });
            Kokkos::fence();
            next->ψxu = 0;
            Kokkos::parallel_reduce("ψ(xₖ₊₁), accelerated step", nthrds, [&](const int i, real_t &ψ_){
                ψ_ += eval_ψ_k(i, *next, μ);
            },next->ψxu); 
            Kokkos::fence();
        };

        // Backtracking line search loop
        while (!stop_signal.stop_requested()) {

            // Recompute step only if τ changed
            if (τ != τ_prev) {
                (τ != 0) && (params.disable_acceleration == false)? take_accelerated_step(τ) : take_safe_step();
                τ_prev = τ;
            } 

            Kokkos::fence();
            // Calculate x̂ₖ₊₁, pₖ₊₁ 
            Kokkos::parallel_for("xûₖ₊₁ = T(xuₖ₊₁)", nthrds, [&](const int i){
                eval_prox(i, *next);
            });
            Kokkos::fence();
            eval_pᵀp_grad_ψᵀp(*next);
            // Calculate ψ(x̂ₖ₊₁)
            next->ψxû = 0;
            Kokkos::parallel_reduce("ψ(xûₖ₊₁)", nthrds, [&](const int i, real_t &ψ_){
                ψ_ += eval_ψ_hat_k(i, *next, μ, y);
            },next->ψxû);

            // Quadratic upper bound
            if (next->L < params.L_max && qub_violated(*next)) { 
                next->γ /= 2;
                next->L *= 2;
                τ = τ_init;
                ++s.stepsize_backtracks;
                continue;
            }

            // Line search condition
            if (τ > 0 && linesearch_violated(*curr, *next)) {
                τ /= 2;
                if (τ < params.τ_min)
                    τ = 0;
                ++s.linesearch_backtracks;
                continue;
            }

            // QUB and line search satisfied
            break;
        }

        // If τ < τ_min the line search failed and we accepted the prox step
        s.linesearch_failures += (τ == 0 && τ_init > 0);
        s.τ_1_accepted += τ == 1;
        s.count_τ += 1;
        s.sum_τ += τ;

        // Check if we made any progress
        if (no_progress > 0 || k % params.max_no_progress == 0)
            no_progress = curr->xu == next->xu ? no_progress + 1 : 0;

        // Update L-BFGS -------------------------------------------------------
        if (curr->γ != next->γ){ // Flush L-BFGS if γ changed
            lbfgs.reset();
        } 
        s.lbfgs_rejected += not lbfgs.update(
                curr->xu, next->xu, curr->grad_ψ, next->grad_ψ,
                LBFGS<config_t>::Sign::Positive, true);

        // Advance step --------------------------------------------------------
        //std::cout<<curr->grad_ψ<<'\n'<<std::endl;
        std::swap(curr, next);
        ++k; 
        xu = curr->xu;
        g = curr->g;
    }
    throw std::logic_error("[PANOC] loop error");
}


} // namespace alpaqa