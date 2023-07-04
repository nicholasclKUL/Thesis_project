#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/print.hpp>

#include <Kokkos_Core.hpp>

#include <thesis/para-panoc.hpp>
#include <thesis/printing.hpp>
#include <thesis/ocp-kkt-error.hpp>

#include <nagumo_schlogl.hpp>

#include <iomanip>
#include <iostream>

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

void progress_callback(const alpaqa::PANOCOCPProgressInfo<config_t>& info){
   auto x = info.x();
}

int main() {

    Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(4));
    
    {

    auto problem = alpaqa::TypeErasedControlProblem<config_t>::make<Nagumo>();

    // Problem dimensions
    
    //SS:
    const auto n = problem.get_N() * problem.get_nu(),
               m = problem.get_N() * problem.get_nc() + problem.get_nc_N();
    
    //MS-Parallel:
    const auto n_ms = (problem.get_nx()+problem.get_nu())*(problem.get_N())+problem.get_nx(),
               m_ms = problem.get_nx()*(problem.get_N()+1),
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
    params.stop_crit = alpaqa::PANOCStopCrit::ProjGradUnitNorm2;
    params.gn_interval = 0; //GN disabled
    params.print_interval = 0;
    params.max_iter = 5000;
    params.disable_acceleration = false;
    params.linesearch_tolerance_factor = 1e-02;
    params.quadratic_upperbound_tolerance_factor = 1e-01;
    params.max_time = std::chrono::minutes(90);
    
    //Outer:
    alpaqa::ALMParams almparams; 
    almparams.tolerance = 1e-5;
    almparams.max_iter = 300;
    almparams.print_interval = 1;
    // almparams.initial_penalty = 1e6;
    almparams.max_time = std::chrono::minutes(60);

    // Solving
    
    //MS:
    alpaqa::ParaALMSolver<alpaqa::ParaPANOCSolver<config_t>> almsolver_ms{almparams,{params}};
    auto stats_ms = almsolver_ms(problem, xu, y_ms, ϵ, nt);
    printing::print_stats_outer(stats_ms);
    alpaqa::KKTiterate<config_t> it(n_ms,m_ms,problem.get_nx(),problem.get_nu()+problem.get_nx());
    auto kkt = alpaqa::compute_kkt_error(problem, it, xu, y_ms, nt);
    printing::kkt_error(kkt);
    printing::print_solution(problem, xu);
    index_t Ts = 1;
    std::string problem_name = "Nagumo";
    printing::output_file(problem, problem_name, xu, Ts);

    //SS:
    params.max_iter = 10000;
    params.gn_interval = 0;
    params.print_interval = 100;
    alpaqa::ALMSolver<alpaqa::PANOCOCPSolver<config_t>> almsolver_ss{almparams,{params}};
    // almsolver_ss.inner_solver.set_progress_callback(progress_callback);
    auto stats_ss = almsolver_ss(problem, u, y);
    printing::print_stats_outer(stats_ss);
    printing::print_solution_ss(problem, u);
    
    Nagumo hg;
    auto kkt_ss = alpaqa::compute_kkt_error(problem, u, hg.Q, hg.R);
    printing::kkt_error(kkt_ss);

    }

    Kokkos::finalize();

}

    