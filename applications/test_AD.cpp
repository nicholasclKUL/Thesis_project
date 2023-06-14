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
#include <thesis/ocp-kkt-error.hpp>

#include <iomanip>
#include <iostream>

int main() {

    Kokkos::initialize(Kokkos::InitializationSettings());

    {

    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

    auto problem = alpaqa::TypeErasedControlProblem<config_t>::make<Quadcopter>();
    auto problemAD = alpaqa::TypeErasedControlProblem<config_t>::make<QuadcopterAD>();

    // Problem dimensions
    const auto n_ms = (problem.get_nx()+problem.get_nu())*(problem.get_N())+problem.get_nx(),
               m_ms = problem.get_nx()*(problem.get_N()+1),
               nt = problem.get_N()+1;

    //MS-Parallel:
    real_t ϵ = 1e-0;            // Continuity tolerance factor
    vec   xu = vec::Zero(n_ms), // Inputs
        y_ms = vec::Zero(m_ms), // Lagrange multipliers
        μ_ms = vec::Ones(m_ms), // Penalty factors
        e_ms(m_ms),             // Constraint violation
        g_ms(m_ms);             // Continuity violation
    problem.get_x_init(xu);

    // Solver
    alpaqa::PANOCOCPParams<config_t> params;
    params.stop_crit      = alpaqa::PANOCStopCrit::FPRNorm;
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

    //Comparing Jacobians
    vec x(12), u(4); 
    mat Jfxu(12,16), JfxuAD(12,16), e(12,16);
    x << -2, 0.1, 3, -.5, .8, -1, .7, 2, 3, .4, -.8, 3; 
    u << 1, -0.5, 2, 0.73;
    Jfxu.setConstant(0.); JfxuAD.setConstant(0.);

    problemAD.eval_jac_f(0, x, u ,Jfxu);
    problemAD.eval_jac_f(0, x, u ,JfxuAD);
    e = Jfxu - Jfxu;

    std::cout<<'\n'<<Jfxu<<'\n'<<std::endl;
    std::cout<<'\n'<<JfxuAD<<'\n'<<std::endl;
    std::cout<<'\n'<<e<<'\n'<<std::endl;

    }

    Kokkos::finalize();
}

    