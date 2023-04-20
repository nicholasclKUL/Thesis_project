#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/print.hpp>

#include <Kokkos_Core.hpp>

#include <Thesis/para-panoc.hpp>
#include <nonlinear_example1.hpp>

#include <iomanip>
#include <iostream>

int main() {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

    auto problem = alpaqa::ocproblem_with_counters(alpaqa::NonlinearOCP1());

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
    vec u = vec::Zero(n); // Inputs
    vec y = vec::Zero(m); // Lagrange multipliers
    vec μ = vec::Ones(m); // Penalty factors
    vec e(m);             // Constraint violation
    problem.get_x_init(u);

    //MS-Parallel
    vec xu = vec::Zero(n_ms);   // Inputs
    vec y_ms = vec::Zero(m_ms); // Lagrange multipliers
    vec μ_ms = vec::Ones(m_ms); // Penalty factors
    vec e_ms(m_ms);             // Constraint violation
    vec g_ms(m_ms);
    problem.get_x_init(xu);

    // Solver
    alpaqa::PANOCOCPParams<config_t> params;
    params.stop_crit      = alpaqa::PANOCStopCrit::ProjGradNorm2;
    params.gn_interval    = 0;
    params.print_interval = 0;
    params.max_iter = 4000;
    params.disable_acceleration = true;
    auto tol = 1e-4;

    // Solve SS
    alpaqa::PANOCOCPSolver<config_t> solver_ss{params};
    auto stats_ss = solver_ss(problem, {.tolerance = tol}, u, y, μ, e);

    // SS statistics
    auto δ      = e.lpNorm<Eigen::Infinity>();
    auto time_s = std::chrono::duration<double>(stats_ss.elapsed_time).count();
    std::cout << '\n'
              << "solver:  " << solver_ss.get_name() << '\n'
              << "status:  " << stats_ss.status << '\n'
              << "ψ = " << alpaqa::float_to_str(stats_ss.final_ψ) << '\n'
              << "ε = " << alpaqa::float_to_str(stats_ss.ε) << '\n'
              << "δ = " << alpaqa::float_to_str(δ) << '\n'
              << "time: " << alpaqa::float_to_str(time_s, 3) << " s\n"
              << "iter:      " << std::setw(6) << stats_ss.iterations << '\n'
              << "line search backtrack: " << std::setw(6)
              << stats_ss.linesearch_backtracks << '\n'
              << "step size backtrack:   " << std::setw(6)
              << stats_ss.stepsize_backtracks << '\n'
              << "solution: ";
    alpaqa::print_python(std::cout, u) << std::endl;          

    //Solve MS-Parallel
    Kokkos::initialize(Kokkos::InitializationSettings());
    alpaqa::ParaPANOCSolver<config_t> solver_ms{params};
    auto stats_ms = solver_ms(problem, {.tolerance = tol}, xu, y_ms, μ_ms, e_ms, g_ms, nt);
    Kokkos::finalize();

    // MS statistics
    auto δ_ms    = e.lpNorm<Eigen::Infinity>();
    auto time_ms = std::chrono::duration<double>(stats_ms.elapsed_time).count();
    std::cout << '\n'
              << "status:  " << stats_ms.status << '\n'
              << "ψ = " << alpaqa::float_to_str(stats_ms.final_ψ) << '\n'
              << "ε = " << alpaqa::float_to_str(stats_ms.ε) << '\n'
              << "δ = " << alpaqa::float_to_str(δ_ms) << '\n'
              << "time: " << alpaqa::float_to_str(time_ms, 3) << " s\n"
              << "iter:      " << std::setw(6) << stats_ms.iterations << '\n'
              << "line search backtrack: " << std::setw(6)
              << stats_ms.linesearch_backtracks << '\n'
              << "step size backtrack:   " << std::setw(6)
              << stats_ms.stepsize_backtracks << '\n'
              << "solution: "<<'\n';

    for (size_t i = 0; i < problem.get_N()+1; ++i){
        std::cout<<"--Stage "<<i<<":"<<'\n';
        if (i < problem.get_N()){
            std::cout<<"x = ["<<xu.segment(i*(problem.get_nx()+problem.get_nu()),problem.get_nx()).transpose()<<"]"<<'\n'
                     <<"u = ["<<xu.segment(i*(problem.get_nx()+problem.get_nu())+problem.get_nx(),problem.get_nu()).transpose()<<"]"<<'\n';
        }
        else{
            std::cout<<"x = ["<<xu.segment(i*(problem.get_nx()+problem.get_nu()),problem.get_nx()).transpose()<<"]"<<std::endl;
        }      
    }
}