#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/print.hpp>

#include <Kokkos_Core.hpp>

#include <thesis/para-panoc.hpp>

// #include <nonlinear_dynamics.hpp>
#include <multi-RTAC.hpp>

#include <iomanip>
#include <iostream>

int main() {

    Kokkos::initialize(Kokkos::InitializationSettings());

    {

    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

    auto problem = alpaqa::TypeErasedControlProblem<config_t>::make<MultiRTAC>();

    // Problem dimensions
    //SS
    const auto n = problem.get_N() * problem.get_nu(),
               m = problem.get_N() * problem.get_nc() + problem.get_nc_N();
    //MS-Parallel
    const auto n_ms = (problem.get_nx()+problem.get_nu())*(problem.get_N())+problem.get_nx(),
               m_ms = problem.get_nx()*(problem.get_N()+1),
               nt = problem.get_N()+1;

    // Initial guess and other solver inputs
    //SS
    vec u = vec::Zero(n); // Inputs
    vec y = vec::Zero(m); // Lagrange multipliers
    vec μ = vec::Ones(m); // Penalty factor
    vec e(m);             // Constraint violation
    problem.get_x_init(u);

    //MS-Parallel
    vec xu = vec::Zero(n_ms);   // Inputs
    vec y_ms = vec::Zero(m_ms); // Lagrange multipliers
    vec μ_ms = vec::Ones(m_ms); // Penalty factors
    // μ_ms.bottomRows(problem.get_nx()) *= 3;
    vec e_ms(m_ms);             // Constraint violation
    vec g_ms(m_ms);
    problem.get_x_init(xu);


    // Solver
    alpaqa::PANOCOCPParams<config_t> params;
    params.stop_crit      = alpaqa::PANOCStopCrit::ProjGradNorm2;
    params.print_interval = 1;
    params.max_iter = 100;
    params.disable_acceleration = false;
    params.linesearch_tolerance_factor = 1e-04;
    params.quadratic_upperbound_tolerance_factor = 1e-03;
    auto tol = 1e-4;
    
    // Solve SS
    params.gn_interval = 1;
    alpaqa::ParaPANOCSolver<config_t> solver_ms_gn{params};
    auto stats_ms_gn = solver_ms_gn(problem, {.tolerance = tol}, xu, y_ms, μ_ms, e_ms, g_ms, nt);

    // Solve MS-Parallel
    problem.get_x_init(xu);
    params.gn_interval = 0;
    alpaqa::ParaPANOCSolver<config_t> solver_ms_lbfgs{params};
    auto stats_ms_lbfgs = solver_ms_lbfgs(problem, {.tolerance = tol}, xu, y_ms, μ_ms, e_ms, g_ms, nt);

    }

    Kokkos::finalize();

}