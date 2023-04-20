#include <alpaqa/config/config.hpp>
#include <Thesis/para-panoc.hpp>
#include <Thesis/para-alm.hpp>
#include <Thesis/nonlinear_example1.hpp>
#include <Kokkos_Core.hpp>
#include <iostream>

#include <alpaqa/config/config.hpp>
#include <alpaqa/util/float.hpp>
#include <alpaqa/problem/box.hpp>

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>



int main(){

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

auto problem = alpaqa::TypeErasedControlProblem<config_t>::make<OCProblem>();

Kokkos::initialize(Kokkos::InitializationSettings());

const auto n = problem.get_N() * problem.get_nu(),

m = problem.get_nx()*(problem.get_N()-1),

m2 = problem.get_N()*problem.get_nc() + problem.get_nc_N(),

nt = problem.get_N(),

nxu = (problem.get_nx()+problem.get_nu())*(problem.get_N()-1)+problem.get_nx();

// Initial guess and other solver inputs

vec u  = vec::Zero(n);      // Inputs (single shooting)
vec xu = vec::Ones(nxu);    // Inputs (multiple shooting)
vec g  = vec::Ones(m);      // constraints g(x,u)=0
vec y  = vec::Ones(m);      // Lagrange multipliers
vec μ  = vec::Ones(m);      // Penalty factors
vec e(m);                   // Constraint violation
problem.get_x_init(xu); 

// Solver Configurations 
// Inner:
alpaqa::PANOCOCPParams<config_t> params;
params.stop_crit = alpaqa::PANOCStopCrit::ProjGradNorm2;
params.gn_interval = 1;
params.print_interval = 0;
params.max_iter = 20;
//params.disable_acceleration;
// Outer:
alpaqa::ALMParams almparams; 
//almparams.ε = 1e-7; // tolerance;
almparams.δ = 1e-7;
almparams.ε_0 = 1e-7;
almparams.max_iter = 1000;
almparams.print_interval = 999;

// Solve
alpaqa::ParaALMSolver<alpaqa::ParaPANOCSolver<config_t>> almsolver{almparams,{params}};
auto stats = almsolver(problem, xu, y, nt);

// Print Solution
std::cout<<'\n'<<"--[FINAL REPORT]--"<<'\n'<<std::endl;

std::cout   << "status: " << stats.status << '\n'
            << "inner iterations: " << stats.inner.iterations << '\n'
            << "outer iterations: " << stats.outer_iterations << '\n'
            << "ε = " << stats.ε << '\n'
            << "δ = " << stats.δ << '\n'
            << "elapsed time:     "
            << std::chrono::duration<double>{stats.elapsed_time}.count()
            << " s" << '\n'
            << "avg τ = " << (stats.inner.sum_τ / stats.inner.count_τ) << '\n'
            << "L-BFGS rejected = " << stats.inner.lbfgs_rejected << '\n'
            << "L-BFGS failures = " << stats.inner.lbfgs_failures << '\n'
            << "Line search failures = " << stats.inner.linesearch_failures
            << '\n'
            << std::endl;

    for (size_t i = 0; i < problem.get_N(); ++i){
        std::cout<<"Stage "<<i<<":"<<'\n';
        if (i < problem.get_N()-1){
            std::cout<<"x = "<<xu.segment(i*(problem.get_nx()+problem.get_nu()),problem.get_nx()).transpose()<<'\n'
                     <<"u = "<<xu.segment(i*(problem.get_nx()+problem.get_nu())+problem.get_nx(),problem.get_nu()).transpose()<<'\n'<<std::endl;
        }
        else{
            std::cout<<"x = "<<xu.segment(i*(problem.get_nx()+problem.get_nu()),problem.get_nx()).transpose()<<std::endl;
        }      
    }
    
Kokkos::finalize();

}