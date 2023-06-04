#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/outer/alm.hpp>
#include <alpaqa/panoc-alm.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/print.hpp>
#include <alpaqa/util/float.hpp>
#include <alpaqa/problem/box.hpp>

#include <Kokkos_Core.hpp>

#include <thesis/para-alm.hpp>
#include <thesis/para-panoc.hpp>
#include <thesis/printing.hpp>

#include <nonlinear_example1.hpp>
#include <linear_dynamics.hpp>
//#include <hanging_chain.hpp>
#include <quadcopter.hpp>
//#include <quadcopter_AD.hpp>

#include <iomanip>
#include <iostream>
#include <cassert>
#include <cmath>


int main() {

    Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(4));

    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

    auto problem = alpaqa::TypeErasedControlProblem<config_t>::make<Quadcopter>();

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
    vec x = vec::Zero(n); // Inputs
    vec y = vec::Zero(m); // Lagrange multipliers
    vec μ = vec::Ones(m); // Penalty factors
    vec e(m);             // Constraint violation
    problem.get_x_init(x);

    //MS-Parallel
    vec xu = vec::Zero(n_ms);   // Inputs
    vec y_ms = vec::Zero(m_ms); // Lagrange multipliers
    vec μ_ms = vec::Ones(m_ms); // Penalty factors
    vec e_ms(m_ms);             // Constraint violation
    vec g_ms(m_ms);             // Continuity violation
    problem.get_x_init(xu);

    // Solver
    // Inner:
    alpaqa::PANOCOCPParams<config_t> params;
    params.stop_crit = alpaqa::PANOCStopCrit::ProjGradNorm2;
    params.gn_interval = 0;
    params.print_interval = 2;
    params.max_iter = 10000;
    params.disable_acceleration = false;
    // Outer:
    alpaqa::ALMParams almparams; 
    almparams.tolerance = 1e-4;
    almparams.max_iter = 1000;
    almparams.print_interval = 2;

    //MS
    alpaqa::ParaALMSolver<alpaqa::ParaPANOCSolver<config_t>> almsolver_ms{almparams,{params}};
    auto stats_ms = almsolver_ms(problem, xu, y_ms, nt);
    Kokkos::finalize();
    printing::print_stats_outer(stats_ms);
    printing::print_solution(problem, xu);

    //SS
    alpaqa::ALMSolver<alpaqa::PANOCOCPSolver<config_t>> almsolver_ss{almparams,{params}};
    auto stats_ss = almsolver_ss(problem, x, y);
    printing::print_stats_outer(stats_ss);
    printing::print_solution_ss(problem,x);

}