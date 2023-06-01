#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/print.hpp>

#include <Kokkos_Core.hpp>

#include <thesis/para-panoc.hpp>
#include <thesis/printing.hpp>
#include <nonlinear_example1.hpp>
#include <linear_dynamics.hpp>
// #include <quadcopter_AD.hpp>
#include <hanging_chain.hpp>
// #include <small_hanging_chain.hpp>

#include <iomanip>
#include <iostream>

int main() {

    Kokkos::initialize(Kokkos::InitializationSettings());
    
    {

        USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

        HangingChain hg;
            hg.a = hg.a/8.46;

        auto problem = alpaqa::TypeErasedControlProblem<config_t>::make<HangingChain>();
        
        // Optimization problem dimensions
        const auto n_ms = (problem.get_nx()+problem.get_nu())*(problem.get_N())+problem.get_nx(),
                m_ms = problem.get_nx()*(problem.get_N()),
                nt = problem.get_N()+1;

        vec xu = vec::Zero(n_ms);   // Inputs
        vec y_ms = vec::Zero(m_ms); // Lagrange multipliers
        vec μ_ms = vec::Ones(m_ms); // Penalty factors
        vec e_ms(m_ms); e_ms.setConstant(alpaqa::NaN<config_t>);
        vec g_ms(m_ms); g_ms.setConstant(alpaqa::NaN<config_t>);

        // Get steady-state for given initial values
        vec xu_ss(problem.get_nx()+problem.get_nu());
        vec fxu(problem.get_nx());
        hg.get_x_init(xu_ss);
        index_t k = 0;
        size_t N = 2e6;
        for (size_t i = 0; i < N; ++i){
            problem.eval_f(k, xu_ss.segment(0,problem.get_nx()), xu_ss.segment(problem.get_nx(),problem.get_nu()), fxu);
            xu_ss.topRows(problem.get_nx()) = fxu;
            xu_ss.bottomRows(problem.get_nu()).setConstant(0.);       
        }

        std::cout<<"States at T = "<<N*hg.Ts<<'\n'<<xu_ss.transpose()<<std::endl;

        // Assign steady-state to initial values of xu
        for (size_t i = 0; i < problem.get_N(); ++i){
            if (i == N-1){
                xu.segment(i*(problem.get_nx()+problem.get_nu()),problem.get_nx()) = fxu;
            }
            else {
                xu.segment(i*(problem.get_nx()+problem.get_nu()),problem.get_nx()+problem.get_nu()) = xu_ss;
            }
        }

        
        // Solver configuration
        // Inner:
        alpaqa::PANOCOCPParams<config_t> params;
        params.stop_crit = alpaqa::PANOCStopCrit::ProjGradNorm2;
        params.gn_interval = 0;
        params.print_interval = 0;
        params.max_iter = 10000;
        params.disable_acceleration = false;
        // Outer:
        alpaqa::ALMParams almparams; 
        almparams.tolerance = 1e-7;
        almparams.max_iter = 1000;
        almparams.print_interval = 0;
        almparams.max_time = std::chrono::minutes(20);

        // Solve 
        alpaqa::ParaALMSolver<alpaqa::ParaPANOCSolver<config_t>> almsolver_ms{almparams,{params}};
        auto stats_ms = almsolver_ms(problem, xu, y_ms, nt);
        printing::print_stats_outer(stats_ms);

        printing::print_solution(problem, xu);

        // index_t stride  = ((nx+dim)/2) - dim,
        //         nxu     = problem.get_nx() + problem.get_nu();      

        // for (size_t i = 0; i < problem.get_N()-1; ++i){
        //     std::cout << std::scientific << xu(i*nxu+((p_Nb-1)*p_dim-3)) 
        //             << ", " << std::setw(24) << std::scientific << xu(i*nxu+((p_Nb-1)*p_dim-2))
        //             << ", " << std::setw(24) << std::scientific << xu(i*nxu+((p_Nb-1)*p_dim-1))<<'\n';
        // }
        // std::cout<<std::endl;

        // std::string problem_name = "HangingChain"; 

        // printing::output_file(problem, problem_name, xu, p_Ts);
    
    }

    Kokkos::finalize();
}

    