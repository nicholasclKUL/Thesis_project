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
#include <hanging_chain.hpp>
#include <quadcopter.hpp>

#include <iomanip>
#include <iostream>
#include <cassert>
#include <cmath>


int main() {

    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

    // Create Problem
    auto problem = alpaqa::TypeErasedControlProblem<config_t>::make<NonlinearOCP1>();

    // Problem dimensions
    const auto n_ms = (problem.get_nx()+problem.get_nu())*(problem.get_N())+problem.get_nx(),
               m_ms = problem.get_nx()*(problem.get_N()),
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
    params.stop_crit = alpaqa::PANOCStopCrit::ProjGradNorm2;
    params.gn_interval = 0;
    params.print_interval = 0;
    params.max_iter = 10000;
    params.disable_acceleration = false;
    // Outer:
    alpaqa::ALMParams almparams; 
    almparams.tolerance = 1e-4;
    almparams.max_iter = 1000;
    almparams.print_interval = 0;

    // Solve
    alpaqa::ParaALMSolver<alpaqa::ParaPANOCSolver<config_t>> almsolver_ms{almparams,{params}};
    Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(1));
    auto stats_ms = almsolver_ms(problem, xu, y_ms, nt);
    Kokkos::finalize();
    printing::print_stats_outer(stats_ms);

}