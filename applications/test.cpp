#include <alpaqa/config/config.hpp>
#include <Thesis/para-panoc.hpp>
#include <Thesis/para-alm.hpp>
#include <Thesis/ocp-funcs.hpp>
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

auto problem = alpaqa::TypeErasedControlProblem<config_t>::make<alpaqa::OCProblem>();

Kokkos::initialize(Kokkos::InitializationSettings());

const auto n = problem.get_N() * problem.get_nu(),

m = problem.get_nx()*(problem.get_N()-1),

m2 = problem.get_N() * problem.get_nc() + problem.get_nc_N(),

nt = problem.get_N(),

nxu = (problem.get_nx()+problem.get_nu())*(problem.get_N()-1)+problem.get_nx();

// Initial guess and other solver inputs

vec u = vec::Zero(n); // Inputs (single shooting)

vec xu = vec::Ones(nxu); // Inputs (multiple shooting)

problem.get_x_init(xu); 

vec g = vec::Ones(m); // constraints g(x,u)=0

vec y = vec::Zero(m); // Lagrange multipliers

vec y2 = vec::Zero(m2);

vec μ = vec::Ones(m); // Penalty factors

vec μ2 = vec::Ones(m2);

vec e(m); // Constraint violation

vec e2(m2);

// Solver

alpaqa::PANOCOCPParams<config_t> params;

params.stop_crit = alpaqa::PANOCStopCrit::ProjGradNorm2;

params.gn_interval = 1;

params.print_interval = 1;

std::cout<<"initial guess x₀: "<<xu.transpose()<<'\n'<<std::endl;

// Solve

alpaqa::ParaPANOCSolver<config_t> solver{params};

auto stats = solver(problem, {.tolerance = 1e-8}, xu, y, μ, e, g, nt);

alpaqa::ALMParams almparams; almparams.ε = 1e-5; // tolerance

alpaqa::ParaALMSolver<alpaqa::ParaPANOCSolver<config_t>> almsolver{almparams,{params}};

auto almstats = almsolver(problem, xu, y, nt);

Kokkos::finalize();

}