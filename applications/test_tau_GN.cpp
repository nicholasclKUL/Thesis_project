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

#include <quadcopter.hpp>

#include <iomanip>
#include <iostream>
#include <cassert>
#include <cmath>


int main() {

    const int np = 16;
    Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(np));

    {

    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

    // Create Problem
    auto problem = alpaqa::TypeErasedControlProblem<config_t>::make<QuadcopterAD>();

    // Problem dimensions
    const auto n_ms = (problem.get_nx()+problem.get_nu())*(problem.get_N())+problem.get_nx(),
               m_ms = problem.get_nx()*(problem.get_N()+1),
               nt = problem.get_N()+1;

    // Containers
    vec xu = vec::Zero(n_ms);   // Inputs
    vec y_ms = vec::Zero(m_ms); // Lagrange multipliers
    vec μ_ms = vec::Ones(m_ms); // Penalty factors
    vec e_ms(m_ms);             // Constraint violation
    vec g_ms(m_ms);             // Continuity violation
    problem.get_x_init(xu);

    // Solver Options
    // Inner:
    alpaqa::PANOCOCPParams<config_t> params;
    params.stop_crit = alpaqa::PANOCStopCrit::ProjGradUnitNorm2;
    params.gn_interval = 0;
    params.print_interval = 0;
    params.max_iter = 5000;
    params.disable_acceleration = false;
    params.linesearch_tolerance_factor = 1e-02;
    params.quadratic_upperbound_tolerance_factor = 1e-01;
    params.max_time = std::chrono::minutes(10);
    // Outer:
    alpaqa::ALMParams almparams; 
    almparams.max_time = std::chrono::minutes(60);
    almparams.tolerance = 1e-4;
    almparams.max_iter = 500;
    almparams.print_interval = 1;
    real_t tol = 1e-0; 

    // Solve
    alpaqa::ParaALMSolver<alpaqa::ParaPANOCSolver<config_t>> almsolver_ms{almparams,{params}};
    auto stats_ms = almsolver_ms(problem, xu, y_ms, tol, nt);
    index_t sim_index = 4;
    std::string problem_name = "NonlinearOCP1b";
    printing::print_solution(problem, problem_name, xu, np, sim_index);
    printing::print_stats_outer(stats_ms, problem, problem_name, np, sim_index);

    }

    Kokkos::finalize();

}