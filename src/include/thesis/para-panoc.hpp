#pragma once

#include <alpaqa/accelerators/lbfgs.hpp>
#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/inner/directions/panoc/lbfgs.hpp>
#include <alpaqa/implementation/inner/panoc.tpp>
#include <alpaqa/implementation/inner/panoc-ocp.tpp>
#include <alpaqa/implementation/inner/panoc-helpers.tpp>

#include <Kokkos_Core.hpp>

#include <thesis/para-panoc-helpers.hpp>

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <chrono>

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
                     index_t nthrds);
    template <class P>
    Stats operator()(const P &problem, const SolveOptions &opts, rvec xu_init, rvec y,
                     crvec μ, rvec e, index_t nthrds) {
        return operator()(Problem::template make<P>(problem), opts, xu_init, y, μ, e, nthrds);
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
    /// [in]    Number of threads
    index_t nthrds) -> Stats {

    using std::chrono::nanoseconds;
    auto os         = opts.os ? opts.os : this->os;
    auto start_time = std::chrono::steady_clock::now();
    Stats s;
    LBFGS<config_t> lbfgs{params.lbfgs_params};

    const auto N    = problem.get_N()+1;    // horizon N
    const auto nu   = problem.get_nu();     // number of inputs
    const auto nx   = problem.get_nx();     // number of states
    const auto nc   = problem.get_nc();     // number of box constraints
    const auto nc_N = problem.get_nc_N();   // box constraints * horizon
    const auto nxu  = nx + nu;              // number of decision variables per stage
    const auto n    = (nxu*(N-1)) + nx;     // total number of decision variables
    const auto m    = nx*(N);               // total number of dynamic constraints

    // Allocate storage --------------------------------------------------------

    vec q(n); // Newton step, including states
    Box<config_t> U   = Box<config_t>::NaN(nu);     // Input box constraints
    Box<config_t> D   = Box<config_t>::NaN(nx);     // State box constraints     
    Box<config_t> D_N = Box<config_t>::NaN(nx);     // Final state box constraints
    Box<config_t> F   = Box<config_t>::NaN(nx);     // Dynamic box constraints (set equal 0)

    F.lowerbound.setConstant(0);
    F.upperbound.setConstant(0);

    // Iterates ----------------------------------------------------------------

    // Represents an iterate in the algorithm, keeping track of some
    // intermediate values and function evaluations.
    struct Iterate {
        
        vec xû;         //< Inputs u interleaved with states x after prox grad
        vec xu;         //< Inputs u interleaved with states x -> x₀, u₀, x₁, u₁,...,uₙ₋₁, xₙ 
        vec grad_ψ;     //< Gradient of cost w.r.t {x₀, u₀, x₁, u₁,...,uₙ₋₁, xₙ}
        vec fxu;        //< Dynamics f(x,u)
        vec hxu;        //< nonconvex mapping h(x,u)
        vec fxû;        //< Dynamics f(x,û)
        vec hxû;        //< nonconvex mapping h(x,û)
        vec lxu;        //< cost function at iterate -> l(x,u)
        vec lxû;        //< cost function after T(x,u) -> l(x,û)
        mat Jfxu;       //< Jacobian of dynamics Jf(x,u)
        vec qr;         //< Cost function gradient Jh(x,u)dl(x,u)
        vec g;          //< Dynamic constraints, f(x,u) - x+
        vec gz;         //< g(x) + Σ⁻¹y
        vec gz_hat;     //< g(x_hat) + Σ⁻¹y
        vec Π_D;        //< Projection of states and variables into their box constraints
        vec Π_D_hat;    //< Projection of states and variables into their box constraints
        vec Z;          //< (Σ*(f(xₖ,uₖ)-xₖ₊₁-Π_D(f(xₖ,uₖ)-xₖ₊₁+Σ⁻¹y))-y)
        vec p;          //< Proximal gradient step
        vec v_GN;
        vec i_GN;
        vec j_GN;

        spmat spGN;         //< GN approximation of ∇²ψ
        mat GN;
        
        real_t ψxu      = alpaqa::NaN<config_t>;        //< Cost in x
        real_t ψxû      = alpaqa::NaN<config_t>;        //< Cost in x̂
        real_t γ        = alpaqa::NaN<config_t>;        //< Step size γ
        real_t L        = alpaqa::NaN<config_t>;        //< Lipschitz estimate L
        real_t pᵀp      = alpaqa::NaN<config_t>;        //< Norm squared of p
        real_t grad_ψᵀp = alpaqa::NaN<config_t>;        //< Dot product of gradient and p

        // @pre    @ref ψxu, @ref pᵀp, @pre grad_ψᵀp
        // @return φγ
        real_t fbe() const { return ψxu + pᵀp / (2 * γ) + grad_ψᵀp; }

        Iterate(length_t n, length_t m, length_t nx, length_t nu, length_t N) :
            xu(n), xû(n), grad_ψ(n), fxu(m-nx), fxû(m-nx), hxu(n), hxû(n), 
            qr(n), g(m), gz(m), gz_hat(m), Π_D(m), Π_D_hat(m), Z(m), p(n), 
            Jfxu(m-nx, nx+nu), spGN(n,n), GN(n,n), v_GN((2*nx+nu)*(2*nx+nu)*N), lxu(N),
            i_GN((2*nx+nu)*(2*nx+nu)*N), j_GN((2*nx+nu)*(2*nx+nu)*N), lxû(N) 
            {
                v_GN.setZero();
                i_GN.setZero();
                j_GN.setZero(); 
                                }

    } iterates[2]{{n, m, nx, nu, N}, {n, m, nx, nu, N}}; 

    Iterate *curr = &iterates[0];
    Iterate *next = &iterates[1];
    
    // Helper functions --------------------------------------------------------

    //create lambdas for Kokkos
    auto eval_iterate = [&](int k, Iterate &It, crvec μ, crvec y){       
        eval_iterate_fun(k, It, problem, μ, y, F);
    };
    auto eval_grad_ψ_k = [&](int k, Iterate &It, crvec μ, crvec y){
        eval_grad_ψ_k_fun(k, It, problem, μ, y);
    };
    auto eval_ψ_k = [&](int k, Iterate &It, crvec μ){
        return eval_ψ_k_fun(k, It, problem, μ, F);
    };
    auto eval_ψ_hat_k = [&](int k, Iterate &It, crvec μ, crvec y){
        return eval_ψ_hat_k_fun(k, It, problem, μ, y, F);        
    };
    auto eval_prox = [&](int k, Iterate &It) {
        eval_prox_fun(k, It, problem, μ, y, D_N, D, U);
    };
    auto eval_prox_stop_crit = [&problem, &D_N, &D, &U](int k, crvec grad_ψ, crvec xu, rvec work_xu, rvec work_p) {
        eval_prox_stop_crit_fun(k, problem, grad_ψ, xu, work_xu, work_p, D_N, D, U);
    };
    auto eval_GN_accelerator = [&](int k, Iterate &It, crvec μ){
        eval_GN_accelerator_fun(k, It, problem, μ);
    };

    //finite difference gradient approximation
    auto fd_grad_ψ = [&](int k, Iterate &It, crvec μ, crvec y) {
        fd_grad_ψ_fun(It, problem, μ, y, F, 0.0001);       
    };

    //vector products for stopping criteria
    auto eval_pᵀp_grad_ψᵀp = [&](Iterate &It) {
        It.pᵀp      = It.p.squaredNorm();
        It.grad_ψᵀp = It.grad_ψ.dot(It.p);
    };

    auto calc_error_stop_crit = [this, &eval_prox_stop_crit](
                                    real_t γ, crvec xuₖ, crvec grad_ψₖ,
                                    crvec pₖ, real_t pₖᵀpₖ, rvec work_xu,
                                    rvec work_p, index_t nthrds) {
        
        switch (params.stop_crit) {
            case PANOCStopCrit::ProjGradNorm: {
                return vec_util::norm_inf(pₖ);
            }
            case PANOCStopCrit::ProjGradNorm2: {
                return std::sqrt(pₖᵀpₖ);
            }
            case PANOCStopCrit::ProjGradUnitNorm: {
                Kokkos::parallel_for("xû = T(xu₀)", nthrds, [&](const int i){
                    eval_prox_stop_crit(i, grad_ψₖ, xuₖ, work_xu, work_p);
                });
                Kokkos::fence();
                return vec_util::norm_inf(xuₖ-work_xu);
            }
            case PANOCStopCrit::ProjGradUnitNorm2: {
                Kokkos::parallel_for("xû = T(xu₀)", nthrds, [&](const int i){
                    eval_prox_stop_crit(i, grad_ψₖ, xuₖ, work_xu, work_p);
                });
                Kokkos::fence();
                work_xu -= xuₖ;
                return work_xu.squaredNorm();
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
        real_t β  = params.linesearch_strictness_factor;
        real_t σ  = β * (1 - curr.γ * curr.L) / (2 * curr.γ);
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
    auto print_progress_1 = [&](Iterate &It, real_t &εₖ, real_t τ, unsigned k) {
        if (k == 0){
            *os << "┌─[ParaPANOC]\n";
        } else {
            *os << "├───── " << k << " ──────" <<'\n';
            *os << "│    ψ = " << print_real(It.ψxu)               
                << ", ‖∇ψ‖ = " << print_real(It.grad_ψ.norm())   
                << ",  ‖p‖ = " << print_real(It.pᵀp)
                << ",    γ = " << print_real(It.γ)               
                << ",    ε = " << print_real(εₖ) << '\n';
        }
    };
    auto print_progress_2 = [&](crvec qₖ, real_t τₖ, bool did_gn) {
        *os << ",  ‖q‖ = " << print_real(qₖ.norm())                       //
            << ",    τ = " << print_real3(τₖ)                             //
            << ",    " << (did_gn ? "GN" : "L-BFGS")                      //
            << std::endl; // Flush for Python buffering
    };
    auto print_progress_n = [&](SolverStatus status) {
        *os << "└─ " << status << " ──"
            << std::endl; // Flush for Python buffering
    };

    // Initialize inputs and initial state ----------------------------

    curr->xu = xu;   
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

    // Initialize lbfgs
    lbfgs.resize(n);
    
    // GN initial status
    unsigned k_gn = 0; 
    bool enable_gn;
    bool enable_gn_global;
    (params.gn_interval > 0) && (params.disable_acceleration == false) ? 
                            enable_gn_global = true : enable_gn_global = false;
    enable_gn = enable_gn_global;
    bool did_gn = enable_gn_global;

    // Keep track of how many successive iterations didn't update the iterate
    unsigned no_progress = 0;

    // Main PANOC loop
    // =========================================================================
    while (true) {

        // Check stop condition ------------------------------------------------

        real_t εₖ = calc_error_stop_crit(curr->γ, curr->xu, curr->grad_ψ,
                                         curr->p, curr->pᵀp, next->xû, next->p, nthrds);

        // Print progress ------------------------------------------------------
        bool do_print =
            params.print_interval != 0 && k % params.print_interval == 0;
        if (do_print){
            print_progress_1(*curr, εₖ, τ ,k);
            print_progress_2(q, τ, did_gn);
        }
        // Return solution -----------------------------------------------------

        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        auto stop_status =
            check_all_stop_conditions(time_elapsed, k, εₖ, no_progress);
        if (stop_status != SolverStatus::Busy) {
            bool do_final_print = params.print_interval != 0;
            if (!do_print && do_final_print){
                print_progress_1(*curr, εₖ, τ, k);
                print_progress_2(q, τ, did_gn);
            }
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

        real_t τ_init = NaN<config_t>;

        // Calculate Gauss-Newton step -----------------------------------------

        if (enable_gn == true){
            curr->v_GN.setZero();
            curr->i_GN.setZero();
            curr->j_GN.setZero();
            curr->GN.setZero();
            curr->spGN.setZero();
            Eigen::SparseLU<spmat> sp_lu;
            Kokkos::parallel_for("GN_hessian", nthrds, [&](const int i){
                eval_GN_accelerator(i, *curr, μ);
            });
            Kokkos::fence();
            fill_sp(curr->spGN,curr->GN);
            sp_lu.analyzePattern(curr->spGN);
            sp_lu.factorize(curr->spGN);
            q = - sp_lu.solve(curr->grad_ψ);
            τ_init = 1;
            k_gn = k + params.gn_interval;
        }

        // Calculate quasi-Newton step -----------------------------------------

        if (enable_gn == false){
            if (k == 0) { // Initialize L-BFGS
                q = curr->p;
                τ_init = 0;
            }
            if (k > 0) {
                q = curr->p;
                τ_init = lbfgs.apply(q, curr->γ)
                        ? 1
                        : 0; 
            }
            // Make sure quasi-Newton step is valid
            if (not q.allFinite()) {
                if (τ_init == 1) { // If we computed a quasi-Newton step
                    ++s.lbfgs_failures;
                    lbfgs.reset(); // Is there anything else we can do?
                }
                τ_init = 0;
            }
        }

        // Line search ---------------------------------------------------------

        next->γ       = curr->γ;
        next->L       = curr->L;
        τ             = τ_init;
        real_t τ_prev = -1;

        // xₖ₊₁ = xₖ + pₖ
        auto take_safe_step = [&] {
            next->xu = curr->xû;  
            next->ψxu = curr->ψxû; 
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
                (τ != 0) && (params.disable_acceleration == false) ? 
                                    take_accelerated_step(τ) : take_safe_step();
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
            Kokkos::fence();

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
                if (τ < params.min_linesearch_coefficient)
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
            no_progress = alpaqa::vec_util::norm_inf(curr->xu - next->xu) <= 0 ? no_progress + 1 : 0;

        // Check whether we used gn step in this PANOC step
        (enable_gn == true) ? did_gn = true : did_gn = false;

        ++k;         

        // Check if the solver should use GN step in the next step
        ((τ == 1) && (enable_gn == true)) || 
        ((k == k_gn) && (enable_gn_global == true)) ? 
                        enable_gn = true : enable_gn = false;

        // Update L-BFGS -------------------------------------------------------
        if ((curr->γ != next->γ) || enable_gn == true){ // Flush L-BFGS if γ changed
            lbfgs.reset();
        } 
        s.lbfgs_rejected += not lbfgs.update(
                curr->xu, next->xu, curr->p, next->p,
                LBFGS<config_t>::Sign::Negative);

        // Advance step --------------------------------------------------------
        std::swap(curr, next);
        xu = curr->xu;
    }
    throw std::logic_error("[PANOC] loop error");
}


} // namespace alpaqa