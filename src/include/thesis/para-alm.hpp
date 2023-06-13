#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/export.hpp>
#include <alpaqa/inner/internal/solverstatus.hpp>
#include <alpaqa/outer/internal/alm-helpers.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#include <alpaqa/outer/alm.hpp>
#include <alpaqa/implementation/outer/alm.tpp>
#include <alpaqa/implementation/outer/internal/alm-helpers.tpp>
#include <alpaqa/implementation/util/print.tpp>
#include <alpaqa/inner/inner-solve-options.hpp>
#include <alpaqa/util/quadmath/quadmath-print.hpp>

#include <chrono>
#include <iostream>
#include <string>
#include <algorithm>
#include <iomanip>
#include <utility>

namespace alpaqa {

/// Augmented Lagrangian Method solver
///
/// @ingroup    grp_ParaALMSolver
template <class InnerSolverT>
class ParaALMSolver {
  public:
    USING_ALPAQA_CONFIG_TEMPLATE(InnerSolverT::config_t);

    using Problem     = TypeErasedControlProblem<config_t>;
    using Params      = ALMParams<config_t>;
    using InnerSolver = InnerSolverT;

    struct Stats {
        unsigned outer_iterations = 0;
        std::chrono::nanoseconds elapsed_time{};
        unsigned initial_penalty_reduced = 0;
        unsigned penalty_reduced = 0;
        unsigned inner_convergence_failures = 0;
        real_t ε = inf<config_t>;
        real_t δ = inf<config_t>;
        real_t norm_penalty = 0;
        SolverStatus status = SolverStatus::Busy;
        InnerStatsAccumulator<typename InnerSolver::Stats> inner;
    };

    ParaALMSolver(Params params, InnerSolver &&inner_solver)
        : params(params),
          inner_solver(std::forward<InnerSolver>(inner_solver)) {}
    ParaALMSolver(Params params, const InnerSolver &inner_solver)
        : params(params), inner_solver(inner_solver) {}

    Stats operator()(const Problem &problem, rvec x, rvec y, real_t ϵ, index_t nthrds);
    template <class P>
    Stats operator()(const P &problem, rvec x, rvec y, real_t ϵ, index_t nthrds) {
        return operator()(Problem::template make<P>(problem), x, y, nthrds);
    }

    std::string get_name() const {
        return "ParaALMSolver<" + inner_solver.get_name() + ">";
    }

    /// Abort the computation and return the result so far.
    /// Can be called from other threads or signal handlers.
    void stop() { inner_solver.stop(); }

    const Params &get_params() const { return params; }

  private:
    Params params;
    using Helpers = detail::ALMHelpers<config_t>;

  public:
    InnerSolver inner_solver;
    std::ostream *os = &std::cout;
};



    // Implementation --------------------------------------------------------


template <class InnerSolverT>
typename ParaALMSolver<InnerSolverT>::Stats
ParaALMSolver<InnerSolverT>::operator()(const Problem &p, rvec x, rvec y, real_t ϵ, index_t nt) {
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    auto start_time = std::chrono::steady_clock::now();

    // Check the problem dimensions etc.
    p.check();

    auto m = (p.get_nx())*(p.get_N()+1);
    if (m == 0) { // No general constraints, only box constraints
        Stats s;
        vec Σ(0), error(0);
        vec g = vec::Constant(m, NaN<config_t>);
        InnerSolveOptions<config_t> opts{
            .always_overwrite_results = true,
            .max_time                 = params.max_time,
            .tolerance                = params.tolerance,
            .os                       = os,
            .check                    = false,
        };
        auto ps              = inner_solver(p, opts, x, y, Σ, error, g, nt);
        bool inner_converged = ps.status == SolverStatus::Converged;
        auto time_elapsed    = std::chrono::steady_clock::now() - start_time;
        s.inner_convergence_failures = not inner_converged;
        s.inner += ps;
        s.ε                = ps.ε;
        s.δ                = 0;
        s.norm_penalty     = 0;
        s.outer_iterations = 1;
        s.elapsed_time     = duration_cast<nanoseconds>(time_elapsed);
        s.status           = ps.status;
        return s;
    }

    constexpr auto NaN               = alpaqa::NaN<config_t>;
    vec Σ                            = vec::Constant(m, NaN);
    vec Σ_old                        = vec::Constant(m, NaN);
    vec error_1                      = vec::Constant(m, NaN);
    vec error_2                      = vec::Constant(m, NaN);
    vec g                            = vec::Constant(m, NaN);
    vec zeros                        = vec::Constant(m, 0.);
    [[maybe_unused]] real_t norm_e_1 = NaN;
    [[maybe_unused]] real_t norm_e_2 = NaN;

    std::array<char, 64> printbuf;
    auto print_real = [&](real_t x) {
        return float_to_str_vw(printbuf, x, params.print_precision);
    };

    Stats s;

    // Initialize the penalty weights
    if (params.initial_penalty > 0) {
        Σ.fill(params.initial_penalty);
    }
    // Initial penalty weights from problem
    else {
        Helpers::initialize_penalty(p, params, x, Σ);
    }
    
    real_t ε                      = params.initial_tolerance;
    [[maybe_unused]] real_t ε_old = NaN;
    real_t Δ                      = params.penalty_update_factor;
    real_t ρ                      = params.tolerance_update_factor;

    unsigned num_successful_iters = 0;

    for (unsigned int i = 0; i < params.max_iter; ++i) {

        // Check if we're allowed to lower the penalty factor even further.
        bool out_of_penalty_factor_updates =
            (num_successful_iters == 0
                 ? s.initial_penalty_reduced == params.max_num_initial_retries
                 : s.penalty_reduced == params.max_num_retries) ||
            (s.initial_penalty_reduced + s.penalty_reduced ==
             params.max_total_num_retries);
        bool out_of_iter = i + 1 == params.max_iter;
        // If this is the final iteration, or the final chance to reduce the
        // penalty update factor, the inner solver can just return its results,
        // even if it doesn't converge.
        bool overwrite_results = out_of_iter || out_of_penalty_factor_updates;

        // Inner solver
        // ------------

        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        InnerSolveOptions<config_t> opts{
            .always_overwrite_results = overwrite_results,
            .max_time                 = params.max_time - time_elapsed,
            .tolerance                = ε,
            .os                       = os,
        };
        // Call the inner solver to minimize the augmented lagrangian for fixed
        // Lagrange multipliers y.
        auto ps = inner_solver(p, opts, x, y, Σ, error_2, g, nt);
        // Update lagrange multipliers
        y += Σ.asDiagonal() * g;
        // Check if the inner solver converged
        bool inner_converged = ps.status == SolverStatus::Converged;
        // Accumulate the inner solver statistics
        s.inner_convergence_failures += not inner_converged;
        s.inner += ps;

        time_elapsed     = std::chrono::steady_clock::now() - start_time;
        bool out_of_time = time_elapsed > params.max_time;
        //bool backtrack = true;
        
        bool backtrack =
            not inner_converged && not overwrite_results && not out_of_time;

        // Print statistics of current iteration
        if (params.print_interval != 0 && i % params.print_interval == 0) {
            real_t δ       = backtrack ? NaN : vec_util::norm_inf(error_2);
            auto color     = inner_converged ? "\x1b[0;32m" : "\x1b[0;31m";
            auto color_end = "\x1b[0m";
            *os << "[\x1b[0;34mALM\x1b[0m]   " << std::setw(5) << i
                << ": ‖Σ‖ = " << print_real(Σ.norm())
                << ", ‖y‖ = " << print_real(y.norm())
                << ", ‖g‖ = " << print_real(vec_util::norm_inf(g))
                << ", δ = " << print_real(δ) << ", ε = " << print_real(ps.ε)
                << ", Δ = " << print_real(Δ) << ", status = " << color
                << std::setw(13) << ps.status << color_end
                << ", iter = " << ps.iterations << std::endl; // Flush for Python buffering
        }

        if (ps.status == SolverStatus::Interrupted) {
            s.ε                = ps.ε;
            s.δ                = vec_util::norm_inf(error_2);
            s.norm_penalty     = Σ.norm();
            s.outer_iterations = i + 1;
            s.elapsed_time     = duration_cast<nanoseconds>(time_elapsed);
            s.status           = ps.status;
            return s;
        }

        // Backtrack and lower penalty if inner solver did not converge
        if (backtrack) {
            if (num_successful_iters > 0) {
                // We have a previous Σ and error
                // Recompute penalty with smaller Δ
                Δ = std::fmax(params.min_penalty_update_factor,
                              Δ * params.penalty_update_factor_lower);
                Helpers::update_penalty_weights(params, Δ, false, error_1,
                                                error_2, norm_e_1, norm_e_2,
                                                Σ_old, Σ, true);
                // Recompute the primal tolerance with larger ρ
                ρ = std::fmin(params.ρ_max, ρ * params.ρ_increase);
                //ε = std::fmax(ρ * ε_old, params.tolerance);
                ++s.penalty_reduced;
            } else {
                // We don't have a previous Σ, simply lower the current Σ and
                // increase ε
                Σ *= params.initial_penalty_lower;
                //ε *= params.initial_tolerance_increase;
                ε_old = std::exchange(ε, std::fmax((1e-3)*std::pow(alpaqa::vec_util::norm_inf(Σ),-1.01), params.tolerance));
                ++s.initial_penalty_reduced;
            }
        }

        // If the inner solver did converge, increase penalty
        else {
            error_2.swap(error_1);
            norm_e_2 = std::exchange(norm_e_1, vec_util::norm_inf(error_1));

            // Check the termination criteria
            real_t continuity = vec_util::norm_inf(g); // continuity of PDEs between stages
            bool alm_converged =
                ((ps.ε <= params.tolerance && inner_converged && norm_e_1 <= params.dual_tolerance) || (continuity <= params.tolerance*ϵ));
            bool exit = alm_converged || out_of_iter || out_of_time;
            if (exit) {
                s.ε                = ps.ε;
                s.δ                = norm_e_1;
                s.norm_penalty     = Σ.norm();
                s.outer_iterations = i + 1;
                s.elapsed_time     = duration_cast<nanoseconds>(time_elapsed);
                s.status           = alm_converged ? SolverStatus::Converged
                                     : out_of_time ? SolverStatus::MaxTime
                                     : out_of_iter ? SolverStatus::MaxIter
                                                   : SolverStatus::Busy;
                return s;
            }
            // After this line, Σ_old contains the penalty used in the current
            // (successful) iteration.
            Σ_old.swap(Σ);
            // Update Σ to contain the penalty to use on the next iteration.
            Helpers::update_penalty_weights(
                params, Δ, num_successful_iters == 0, error_1, error_2,
                norm_e_1, norm_e_2, Σ_old, Σ, true);
            // Lower the primal tolerance for the inner solver.
            ε_old = std::exchange(ε, std::fmax((1e-3)*std::pow(alpaqa::vec_util::norm_inf(Σ),-1.01), params.tolerance));
            // ε_old = std::exchange(ε, std::fmax(ρ * ε, params.tolerance));
            ++num_successful_iters;
        }
    }
    throw std::logic_error("[ALM]   loop error");
}

} // namespace alpaqa