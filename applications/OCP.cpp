#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/print.hpp>

#include <Kokkos_Core.hpp>

#include <thesis/para-panoc.hpp>
#include <thesis/para-alm.hpp>
#include <thesis/printing.hpp>
#include <nagumo_schlogl_ws.hpp>

#include <iomanip>
#include <iostream>

int main() {

    const int np = 16;
    Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(np));
    
    {

        USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

        // Create Problem
        auto problem = alpaqa::TypeErasedControlProblem<config_t>::make<Nagumo>();

        // Problem dimensions
        const auto n_ms = (problem.get_nx()+problem.get_nu())*(problem.get_N())+problem.get_nx(),
                    m_ms = problem.get_nx()*(problem.get_N()+1),
                    nt = problem.get_N()+1;

        // Containers
        vec xu = vec::Zero(n_ms);   // Inputs
        vec y_ms = vec::Zero(m_ms); // Lagrange multipliers
        problem.get_x_init(xu);

        // Solver Options
        //Inner:
        alpaqa::PANOCOCPParams<config_t> params;
        params.stop_crit = alpaqa::PANOCStopCrit::ProjGradUnitNorm2;
        params.gn_interval = 0; //GN disabled
        params.print_interval = 0;
        params.max_iter = 10000;
        params.disable_acceleration = false;
        params.linesearch_tolerance_factor = 1e-02;
        params.quadratic_upperbound_tolerance_factor = 1e-02;
        params.max_time = std::chrono::minutes(2);   
        //Outer:
        alpaqa::ALMParams almparams; 
        almparams.tolerance = 1e-5;
        almparams.max_iter = 300;
        almparams.print_interval = 1;
        almparams.max_time = std::chrono::minutes(10);
        real_t tol = 1;
        
        // Solve 
        alpaqa::ParaALMSolver<alpaqa::ParaPANOCSolver<config_t>> almsolver_ms{almparams,{params}};
        auto stats_ms = almsolver_ms(problem, xu, y_ms, tol, nt);
        printing::print_stats_outer(stats_ms);
        printing::print_solution(problem, xu);
    
    }

    Kokkos::finalize();
}

    