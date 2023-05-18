#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/print.hpp>

#include <Kokkos_Core.hpp>

#include <thesis/para-panoc.hpp>
#include <thesis/printing.hpp>
#include <nonlinear_example1.hpp>
#include <linear_dynamics.hpp>
#include <quadcopter_AD.hpp>
#include <quadcopter.hpp>

#include <iomanip>
#include <iostream>

int main() {

    Kokkos::initialize(Kokkos::InitializationSettings());

    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

    auto problem = alpaqa::TypeErasedControlProblem<config_t>::make<Quadcopter>();
    auto problemAD = alpaqa::TypeErasedControlProblem<config_t>::make<QuadcopterAD>();

    // Problem dimensions
    const auto n_ms = (problem.get_nx()+problem.get_nu())*(problem.get_N())+problem.get_nx(),
               m_ms = problem.get_nx()*(problem.get_N()),
               nt = problem.get_N()+1;

    vec xu = vec::Zero(n_ms);   // Inputs
    vec y_ms = vec::Zero(m_ms); // Lagrange multipliers
    vec μ_ms = vec::Ones(m_ms); // Penalty factors
    vec e_ms(m_ms); e_ms.setConstant(alpaqa::NaN<config_t>);
    vec g_ms(m_ms); g_ms.setConstant(alpaqa::NaN<config_t>);
    problem.get_x_init(xu);

    // Solver
    alpaqa::PANOCOCPParams<config_t> params;
    params.stop_crit      = alpaqa::PANOCStopCrit::ProjGradNorm2;
    params.gn_interval    = 0;
    params.print_interval = 1;
    params.max_iter = 100;
    params.disable_acceleration = false;
    auto tol = 1e-4;
    
    //Solve 
    alpaqa::ParaPANOCSolver<config_t> solver_ms{params};
    auto stats_ms = solver_ms(problem, {.tolerance = tol}, xu, y_ms, μ_ms, e_ms, g_ms, nt);
    printing::print_stats_inner(stats_ms, e_ms);

    // AD
    //Reset Values
    xu = vec::Zero(n_ms);   // Inputs
    y_ms = vec::Zero(m_ms); // Lagrange multipliers
    μ_ms = vec::Ones(m_ms); // Penalty factors
    e_ms.setConstant(alpaqa::NaN<config_t>);
    g_ms.setConstant(alpaqa::NaN<config_t>);
    problem.get_x_init(xu);
    //Solve
    alpaqa::ParaPANOCSolver<config_t> solver_ms_AD{params};
    auto stats_ms_AD = solver_ms_AD(problemAD, {.tolerance = tol}, xu, y_ms, μ_ms, e_ms, g_ms, nt);
    printing::print_stats_inner(stats_ms_AD, e_ms);

    Kokkos::finalize();
}

    