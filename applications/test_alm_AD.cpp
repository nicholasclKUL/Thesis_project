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
#include <quadcopter.hpp>
#include <thesis/ocp-kkt-error.hpp>

#include <iomanip>
#include <iostream>

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

void progress_callback(const alpaqa::PANOCOCPProgressInfo<config_t>& info){
   auto states = info.x();
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
    real_t ϵ = 1e-0;            // Continuity tolerance factor
    vec   xu = vec::Zero(n_ms), // Inputs
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
    almparams.max_iter = 10000;
    almparams.print_interval = 0;

    // Solving
    
    //MS:
    alpaqa::ParaALMSolver<alpaqa::ParaPANOCSolver<config_t>> almsolver_ms{almparams,{params}};
    auto stats_ms = almsolver_ms(problem, xu, y_ms, ϵ, nt);
    printing::print_stats_outer(stats_ms);
    printing::print_solution(problem, xu);

    //SS:
    alpaqa::ALMSolver<alpaqa::PANOCOCPSolver<config_t>> almsolver_ss{almparams,{params}};
    almsolver_ss.inner_solver.set_progress_callback(progress_callback);
    auto stats_ss = almsolver_ss(problem, u, y);
    printing::print_stats_outer(stats_ss);
    printing::print_solution_ss(problem, u);

    // Check Optimality Conditions

    //MS:
    alpaqa::KKTiterate<config_t> it(n_ms,m_ms,problem.get_nu()+problem.get_nx());
    auto kkt = alpaqa::compute_kkt_error(problem, it, xu, y_ms, nt);
    printing::kkt_error(kkt);

    // //simulate states for ss control action:
    // vec x(problem.get_nx()), fxu(problem.get_nx()); 
    // problem.get_x_init(x);
    // for (size_t i = 0; i < problem.get_N(); ++i){
    //     problem.eval_f(i, x, u.segment(i*problem.get_nu(),problem.get_nu()),fxu);
    //     x = fxu;
    //     std::cout<<std::scientific<<"["<<fxu.transpose()<<"]"<<std::endl;
    // }
    

    // std::cout<<alpaqa::vec_util::norm_inf(y_ms)<<", "<<alpaqa::vec_util::norm_inf(it.qr)<<", "
    //         <<it.Jfxu.lpNorm<Eigen::Infinity>()<<'\n'<<std::endl;

    // std::cout<<"it.hxu"<<std::setw(16)<<"it.qr"<<std::setw(16)<<"it.grad_L"<<'\n';
    // index_t k = 0;
    // for (size_t i = 0; i < problem.get_N(); ++i){
    //     for (size_t j = 0; j < problem.get_nx()+problem.get_nu(); j++){
    //         std::cout<<it.hxu(i*(problem.get_nx()+problem.get_nu())+j)<<std::setw(16)
    //                  <<it.qr(i*(problem.get_nx()+problem.get_nu())+j)<<std::setw(16)
    //                  <<it.grad_L(i*(problem.get_nx()+problem.get_nu())+j)<<'\n';
    //     }
    //     std::cout<<'\n'<<std::endl;
    // }
    // std::cout<<std::endl;

    }

    Kokkos::finalize();

}

    