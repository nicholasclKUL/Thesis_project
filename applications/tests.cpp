#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/print.hpp>

#include <Kokkos_Core.hpp>

#include <thesis/para-panoc.hpp>
#include <thesis/printing.hpp>
#include <nonlinear_example1.hpp>
#include <linear_dynamics.hpp>

#include <iomanip>
#include <iostream>

int main() {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

    alpaqa::ControlProblemWithCounters<LinearOCP> problem;

    // Problem dimensions
    //SS
    const auto n = problem.get_N() * problem.get_nu(),
               m = problem.get_N() * problem.get_nc() + problem.get_nc_N();
    //MS-Parallel
    const auto n_ms = (problem.get_nx()+problem.get_nu())*(problem.get_N())+problem.get_nx(),
               m_ms = problem.get_nx()*(problem.get_N()),
               nt = problem.get_N()+1;

    // Initial guess and other solver inputs
    //SS
    vec u = vec::Zero(n); // Inputs
    vec y = vec::Zero(m); // Lagrange multipliers
    vec μ = vec::Ones(m); // Penalty factors
    vec e(m);             // Constraint violation
    problem.get_x_init(u);

    //MS-Parallel
    vec xu = vec::Zero(n_ms);   // Inputs
    vec y_ms = vec::Zero(m_ms); // Lagrange multipliers
    vec μ_ms = vec::Ones(m_ms); // Penalty factors
    vec e_ms(m_ms);             // Constraint violation
    vec g_ms(m_ms);
    problem.get_x_init(xu);

    // Solver
    alpaqa::PANOCOCPParams<config_t> params;
    params.stop_crit      = alpaqa::PANOCStopCrit::ProjGradNorm2;
    params.gn_interval    = 1;
    params.print_interval = 1;
    params.max_iter = 100;
    params.disable_acceleration = false;
    auto tol = 1e-4;

    // Solve SS
    // alpaqa::PANOCOCPSolver<config_t> solver_ss{params};
    // auto stats_ss = solver_ss(problem, {.tolerance = tol}, u, y, μ, e);
    // //printing statistics:
    // printing::print_stats_inner(stats_ss, e);
    
    //Solve MS-Parallel
    Kokkos::initialize(Kokkos::InitializationSettings());
    alpaqa::ParaPANOCSolver<config_t> solver_ms{params};
    auto stats_ms = solver_ms(problem, {.tolerance = tol}, xu, y_ms, μ_ms, e_ms, g_ms, nt);
    Kokkos::finalize();
    // MS statistics
    printing::print_stats_inner(stats_ms, e_ms);

}