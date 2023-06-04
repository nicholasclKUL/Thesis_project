#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/print.hpp>

#include <Kokkos_Core.hpp>

#include <thesis/para-panoc.hpp>
#include <thesis/printing.hpp>
#include <thesis/ocp-kkt-error.hpp>
#include <nonlinear_example1.hpp>
#include <linear_dynamics.hpp>
//#include <quadcopter_AD.hpp>
#include <hanging_chain.hpp>
#include <quadcopter.hpp>

#include <iomanip>
#include <iostream>

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

void progress_callback(const alpaqa::PANOCOCPProgressInfo<config_t>& info){
   auto states = info.x();
   std::cout<<"["<<states.transpose()<<"]"<<'\n\n'<<std::endl;
}

int main() {

    Kokkos::initialize(Kokkos::InitializationSettings());
    
    {

    auto problem = alpaqa::TypeErasedControlProblem<config_t>::make<HangingChain>();

    // Problem dimensions
    
    //SS:
    const auto n = problem.get_N() * problem.get_nu(),
               m = problem.get_N() * problem.get_nc() + problem.get_nc_N();
    
    //MS-Parallel:
    const auto n_ms = (problem.get_nx()+problem.get_nu())*(problem.get_N())+problem.get_nx(),
               m_ms = problem.get_nx()*(problem.get_N()),
               nt = problem.get_N()+1;

    // Initial guess and other solver inputs
    
    //SS:
    vec u = vec::Zero(n), // Inputs
        y = vec::Zero(m), // Lagrange multipliers
        μ = vec::Ones(m), // Penalty factors
        e(m);             // Constraint violation

    //MS-Parallel:
    vec xu = vec::Zero(n_ms),   // Inputs
        y_ms = vec::Zero(m_ms), // Lagrange multipliers
        μ_ms = vec::Ones(m_ms), // Penalty factors
        e_ms(m_ms),             // Constraint violation
        g_ms(m_ms);             // Continuity violation
    problem.get_x_init(xu);

    // Solver
    
    //Inner:
    alpaqa::PANOCOCPParams<config_t> params;
    params.stop_crit = alpaqa::PANOCStopCrit::ProjGradNorm2;
    params.gn_interval = 0; //GN disabled
    params.print_interval = 0;
    params.max_iter = 10000;
    params.disable_acceleration = false;
    
    //Outer:
    alpaqa::ALMParams almparams; 
    almparams.tolerance = 1e-6;
    almparams.max_iter = 1000;
    almparams.print_interval = 0;

    // Solving
    
    //MS:
    alpaqa::ParaALMSolver<alpaqa::ParaPANOCSolver<config_t>> almsolver_ms{almparams,{params}};
    auto stats_ms = almsolver_ms(problem, xu, y_ms, nt);
    printing::print_stats_outer(stats_ms);
    printing::print_solution(problem, xu);
    
    //SS:
    
    alpaqa::ALMSolver<alpaqa::PANOCOCPSolver<config_t>> almsolver_ss{almparams,{params}};
    almsolver_ss.inner_solver.set_progress_callback(progress_callback);
    auto stats_ss = almsolver_ss(problem, u, y);
    printing::print_stats_outer(stats_ss);
    printing::print_solution_ss(problem, u);
        
    }

    Kokkos::finalize();

}

    