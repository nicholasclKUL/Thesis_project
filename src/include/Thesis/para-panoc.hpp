#pragma once

#include <alpaqa/inner/panoc.hpp>
#include <Kokkos_Core.hpp>
#include <Thesis/para-type-erased-problem.hpp>

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include <alpaqa/config/config.hpp>
#include <alpaqa/implementation/inner/panoc-helpers.tpp>
#include <alpaqa/implementation/util/print.tpp>
#include <alpaqa/util/alloc-check.hpp>
#include <alpaqa/util/quadmath/quadmath-print.hpp>

namespace alpaqa {


/* Create struct to carry parallel parameters and data structures
 to be passed to PANOCSolver */


/// PANOC solver for ALM.
/// @ingroup    grp_InnerSolvers
template <class DirectionT>
class ParaPANOCSolver {
  public:
    USING_ALPAQA_CONFIG_TEMPLATE(DirectionT::config_t);

    using Problem      = ParaTypeErasedProblem<config_t>;
    using Params       = PANOCParams<config_t>;
    using Direction    = DirectionT;
    using Stats        = PANOCStats<config_t>;
    using ProgressInfo = PANOCProgressInfo<config_t>;
    using SolveOptions = InnerSolveOptions<config_t>;

    ParaPANOCSolver(const Params &params)
        requires std::default_initializable<Direction>
        : params(params) {}
    ParaPANOCSolver(const Params &params, Direction &&direction)
        : params(params), direction(std::move(direction)) {}
    ParaPANOCSolver(const Params &params, const Direction &direction)
        : params(params), direction(direction) {}

    Stats operator()(const Problem &problem,   // in
                     const SolveOptions &opts, // in
                     rvec x,                   // inout
                     rvec y,                   // inout
                     crvec Σ,                  // in
                     rvec err_z,               // out
                     unsigned int &nprocs);

    template <class P>
    Stats operator()(const P &problem, const SolveOptions &opts, rvec u, rvec y,
                     crvec Σ, rvec e, unsigned int &nprocs) {
        return operator()(Problem::template make<P>(problem), opts, u, y, Σ, e, nprocs);
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

  public:
    Params params;
    AtomicStopSignal stop_signal;
    std::function<void(const ProgressInfo &)> progress_cb;
    using Helpers = detail::PANOCHelpers<config_t>;
    Direction direction;
    std::ostream *os = &std::cout;
};



    // Implementation --------------------------------------------------------



template <class DirectionProviderT>
std::string ParaPANOCSolver<DirectionProviderT>::get_name() const {
    return "ParaPANOCSolver<" + std::string(direction.get_name()) + ">";
}

template <class DirectionProviderT>
auto ParaPANOCSolver<DirectionProviderT>::operator()(
    /// [in]    Problem description
    const Problem &problem,
    /// [in]    Solve options
    const SolveOptions &opts,
    /// [inout] Decision variable @f$ x @f$
    rvec x,
    /// [inout] Lagrange multipliers @f$ y @f$
    rvec y,
    /// [in]    Constraint weights @f$ \Sigma @f$
    crvec Σ,
    /// [out]   Slack variable error @f$ g(x) - \Pi_D(g(x) + \Sigma^{-1} y) @f$
    rvec err_z,
    /// [in]    Number of processes
    unsigned int &nprocs) -> Stats {
 
    using std::chrono::nanoseconds;
    auto os         = opts.os ? opts.os : this->os;
    auto start_time = std::chrono::steady_clock::now();
    Stats s;
 
    const auto n = problem.get_n();
    const auto m = problem.get_m();
 
    // Represents an iterate in the algorithm, keeping track of some
    // intermediate values and function evaluations.
    struct Iterate {
        vec x;      //< Decision variables
        vec x̂;      //< Decision variables after proximal gradient step
        vec grad_ψ; //< Gradient of cost in x
        vec p;      //< Proximal gradient step in x
        vec ŷx̂;     //< Candidate Lagrange multipliers in x̂
        real_t ψx       = NaN<config_t>; //< Cost in x
        real_t ψx̂       = NaN<config_t>; //< Cost in x̂
        real_t γ        = NaN<config_t>; //< Step size γ
        real_t L        = NaN<config_t>; //< Lipschitz estimate L
        real_t pᵀp      = NaN<config_t>; //< Norm squared of p
        real_t grad_ψᵀp = NaN<config_t>; //< Dot product of gradient and p
        real_t hx̂       = NaN<config_t>; //< Non-smooth function value in x̂
 
        // @pre    @ref ψx, @ref hx̂ @ref pᵀp, @ref grad_ψᵀp
        // @return φγ
        real_t fbe() const { return ψx + hx̂ + pᵀp / (2 * γ) + grad_ψᵀp; }
 
        Iterate(length_t n, length_t m) : x(n), x̂(n), grad_ψ(n), p(n), ŷx̂(m) {}
    } iterates[2]{{n, m}, {n, m}};
    Iterate *curr = &iterates[0];
    Iterate *next = &iterates[1];
 
    bool need_grad_ψx̂ = Helpers::stop_crit_requires_grad_ψx̂(params.stop_crit);
    vec grad_ψx̂(n);
    vec work_n(n), work_m(m);
    vec q(n); // (quasi-)Newton step Hₖ pₖ
 
    // Helper functions --------------------------------------------------------
 
    auto qub_violated = [this](const Iterate &i) { // quadratic upperbound 
        real_t margin =
            (1 + std::abs(i.ψx)) * params.quadratic_upperbound_tolerance_factor;
        return i.ψx̂ > i.ψx + i.grad_ψᵀp + real_t(0.5) * i.L * i.pᵀp + margin;
    };
 
    auto linesearch_violated = [this](const Iterate &curr,
                                      const Iterate &next) {
        if (params.force_linesearch)
            return false;
        real_t σ  = params.β * (1 - curr.γ * curr.L) / (2 * curr.γ);
        real_t φγ = curr.fbe();
        real_t margin = (1 + std::abs(φγ)) * params.linesearch_tolerance_factor;
        return next.fbe() > φγ - σ * curr.pᵀp + margin;
    };
 
    // Problem functions -------------------------------------------------------
 
    auto eval_ψ_grad_ψ = [&problem, &y, &Σ, &work_n, &work_m](Iterate &i, unsigned int &nprocs) {
        i.ψx = problem.eval_ψ_grad_ψ(i.x, y, Σ, i.grad_ψ, work_n, work_m, nprocs);
    };
    auto eval_prox_grad_step = [&problem](Iterate &i, unsigned int &nprocs) {
        i.hx̂  = problem.eval_prox_grad_step(i.γ, i.x, i.grad_ψ, i.x̂, i.p, nprocs);
        i.pᵀp = i.p.squaredNorm();
        i.grad_ψᵀp = i.p.dot(i.grad_ψ);
    };
    auto eval_ψx̂ = [&problem, &y, &Σ](Iterate &i, unsigned int &nprocs) {
        i.ψx̂ = problem.eval_ψ(i.x̂, y, Σ, i.ŷx̂, nprocs);
    };
    auto eval_grad_ψx̂ = [&problem, &work_n](Iterate &i, rvec grad_ψx̂, unsigned int &nprocs) {
        problem.eval_grad_ψ_from_ŷ(i.x̂, i.ŷx̂, grad_ψx̂, work_n, nprocs);
    };
 
    // Printing ----------------------------------------------------------------
 
    std::array<char, 64> print_buf;
    auto print_real = [this, &print_buf](real_t x) {
        return float_to_str_vw(print_buf, x, params.print_precision);
    };
    auto print_real3 = [&print_buf](real_t x) {
        return float_to_str_vw(print_buf, x, 3);
    };
    auto print_progress = [&](unsigned k, real_t φγ, real_t ψ, crvec grad_ψ,
                              real_t pᵀp, crvec q, real_t γ, real_t τ,
                              real_t ε) {
        *os << "[PANOC] " << std::setw(6) << k << ": φγ = " << print_real(φγ)
            << ", ψ = " << print_real(ψ)
            << ", ‖∇ψ‖ = " << print_real(grad_ψ.norm())
            << ", ‖p‖ = " << print_real(std::sqrt(pᵀp))
            << ", γ = " << print_real(γ) << ", ε = " << print_real(ε);
        if (k > 0)
            *os << ", τ = " << print_real3(τ)
                << ", ‖q‖ = " << print_real(q.norm());
        *os << std::endl; // Flush for Python buffering
    };

    // Initialization ----------------------------------------------------------
    
    curr->x = x;
 
    // Estimate Lipschitz constant ---------------------------------------------
 
    // Finite difference approximation of ∇²ψ in starting point 
    if (params.Lipschitz.L_0 <= 0) {
            curr->L = Helpers::initial_lipschitz_estimate(
                problem, curr->x, y, Σ, params.Lipschitz.ε, params.Lipschitz.δ,
                params.L_min, params.L_max,
                /* in ⟹ out */ curr->ψx, curr->grad_ψ, curr->x̂, next->grad_ψ,
                work_n, work_m, nprocs);
        }

    // Initial Lipschitz constant provided by the user
    else {
        curr->L = params.Lipschitz.L_0;
        // Calculate ψ(xₖ), ∇ψ(x₀)
        eval_ψ_grad_ψ(*curr, nprocs);
    }
    if (not std::isfinite(curr->L)) {
        s.status = SolverStatus::NotFinite;
        return s;
    }
    curr->γ = params.Lipschitz.Lγ_factor / curr->L;
 
    // First proximal gradient step --------------------------------------------
 
    eval_prox_grad_step(*curr, nprocs);
    eval_ψx̂(*curr, nprocs);
 
    // Quadratic upper bound
    while (curr->L < params.L_max && qub_violated(*curr)) {
        curr->γ /= 2;
        curr->L *= 2;
        eval_prox_grad_step(*curr, nprocs);
        eval_ψx̂(*curr, nprocs);
    }

    // Loop data ---------------------------------------------------------------
 
    real_t τ_init = NaN<config_t>;
    unsigned k = 0;             // iteration
    real_t τ   = NaN<config_t>; // line search parameter
    // Keep track of how many successive iterations didn't update the iterate
    unsigned no_progress = 0;


    // Main PANOC loop
    // =========================================================================
 
    ScopedMallocBlocker mb; // Don't allocate in the inner loop
    while (true) {

        // Calculate ∇ψ(x̂ₖ)
        if (need_grad_ψx̂)
            eval_grad_ψx̂(*curr, grad_ψx̂, nprocs);
        bool have_grad_ψx̂ = need_grad_ψx̂;

        // Check stopping criteria ---------------------------------------------
 
        real_t εₖ = Helpers::calc_error_stop_crit(
            problem, params.stop_crit, curr->p, curr->γ, curr->x, curr->x̂,
            curr->ŷx̂, curr->grad_ψ, grad_ψx̂, work_n, next->p);
 
        // Print progress ------------------------------------------------------
 
        if (params.print_interval != 0 && k % params.print_interval == 0)
            print_progress(k, curr->fbe(), curr->ψx, curr->grad_ψ, curr->pᵀp, q,
                           curr->γ, τ, εₖ);
        if (progress_cb) {
            ScopedMallocAllower ma;
            progress_cb({.k          = k,
                         .x          = curr->x,
                         .p          = curr->p,
                         .norm_sq_p  = curr->pᵀp,
                         .x̂          = curr->x̂,
                         .φγ         = curr->fbe(),
                         .ψ          = curr->ψx,
                         .grad_ψ     = curr->grad_ψ,
                         .ψ_hat      = curr->ψx̂,
                         .grad_ψ_hat = grad_ψx̂,
                         .q          = q,
                         .L          = curr->L,
                         .γ          = curr->γ,
                         .τ          = τ,
                         .ε          = εₖ,
                         .Σ          = Σ,
                         .y          = y,
                         .problem    = problem,
                         .params     = params});
        }
 
        // Return solution -----------------------------------------------------
 
        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        auto stop_status  = Helpers::check_all_stop_conditions(
            params, opts, time_elapsed, k, stop_signal, εₖ, no_progress);
        if (stop_status != SolverStatus::Busy) {
            // TODO: move the computation of ẑ and g(x) to ALM?
            if (stop_status == SolverStatus::Converged ||
                stop_status == SolverStatus::Interrupted ||
                opts.always_overwrite_results) {
                auto &ŷ = curr->ŷx̂;
                if (err_z.size() > 0)
                    err_z = Σ.asDiagonal().inverse() * (ŷ - y);
                x = std::move(curr->x̂);
                y = std::move(curr->ŷx̂);
            }
            s.iterations   = k;
            s.ε            = εₖ;
            s.elapsed_time = duration_cast<nanoseconds>(time_elapsed);
            s.status       = stop_status;
            s.final_γ      = curr->γ;
            s.final_ψ      = curr->ψx̂;
            s.final_h      = curr->hx̂;
            s.final_φγ     = curr->fbe();
            return s;
        }
 
        // Calculate quasi-Newton step -----------------------------------------
 
        if (k == 0) { // Initialize L-BFGS
            ScopedMallocAllower ma;
            direction.initialize(problem, y, Σ, curr->γ, curr->x, curr->x̂,
                                 curr->p, curr->grad_ψ);
            τ_init = 0;
        }
        if (k > 0 || direction.has_initial_direction()) {
            τ_init = 1;
            direction.apply(curr->γ, curr->x, curr->x̂, curr->p, curr->grad_ψ,
                            q);
        }
        // Make sure quasi-Newton step is valid
        if (not q.allFinite()) {
            if (τ_init == 1) { // If we computed a quasi-Newton step
                ++s.lbfgs_failures;
                direction.reset(); // Is there anything else we can do?
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
            // Calculate ∇ψ(xₖ₊₁)
            if (not have_grad_ψx̂)
            Kokkos::parallel_for(nprocs, [=] (const unsigned int i) {
                eval_grad_ψx̂(*curr, grad_ψx̂, nprocs);
            });
            have_grad_ψx̂ = true;
            next->x      = curr->x̂; // → safe prox step
            next->ψx     = curr->ψx̂;
            next->grad_ψ.swap(grad_ψx̂);
        };
 
        // xₖ₊₁ = xₖ + (1-τ) pₖ + τ qₖ
        auto take_accelerated_step = [&](real_t τ) {
            if (τ == 1) // → faster quasi-Newton step
                next->x = curr->x + q;
            else
                next->x = curr->x + (1 - τ) * curr->p + τ * q;
            // Calculate ψ(xₖ₊₁), ∇ψ(xₖ₊₁)
            eval_ψ_grad_ψ(*next, nprocs);
        };

        while (!stop_signal.stop_requested()) {
 

            // Recompute step only if τ changed
            if (τ != τ_prev) {
                τ != 0 ? take_accelerated_step(τ) : take_safe_step();
                τ_prev = τ;
            }
 
            // If the cost is not finite, abandon the direction entirely, don't
            // even bother backtracking.
            if (τ > 0 && !std::isfinite(next->ψx)) {
                τ = 0;
                direction.reset();
                continue;
            }

            // Calculate x̂ₖ₊₁, ψ(x̂ₖ₊₁)
            eval_prox_grad_step(*next, nprocs);
            eval_ψx̂(*next, nprocs);

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
            no_progress = curr->x == next->x ? no_progress + 1 : 0;
 
        // Update L-BFGS -------------------------------------------------------
 
        if (curr->γ != next->γ) // Flush L-BFGS if γ changed
            direction.changed_γ(next->γ, curr->γ);
 
        s.lbfgs_rejected +=
            not direction.update(curr->γ, next->γ, curr->x, next->x, curr->p,
                                 next->p, curr->grad_ψ, next->grad_ψ);
 
        // Advance step --------------------------------------------------------
        std::swap(curr, next);
        ++k;
    }
    throw std::logic_error("[PANOC] loop error");
}


} // namespace alpaqa