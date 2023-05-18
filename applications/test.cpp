#include <quadcopter_AD.hpp>
#include <quadcopter.hpp>

#include <iostream>

int main() {

Kokkos::initialize(Kokkos::InitializationSettings());

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

auto problem_ad = alpaqa::TypeErasedControlProblem<config_t>::make<QuadcopterAD>();
auto problem = alpaqa::TypeErasedControlProblem<config_t>::make<Quadcopter>();

vec x, u;
mat J_ad(12,16);
mat J(12,16);

x = vec::Ones(12);
u = vec::Ones(4);

problem.eval_jac_f(0, x, u, J);
problem_ad.eval_jac_f(0, x, u, J_ad);

for (size_t i = 0; i < problem.get_nx(); i++){
    for (size_t j = 0; j < problem.get_nx()+problem.get_nu(); j++){
        std::cout<<std::scientific<<J(i,j)-J_ad(i,j)<<",  ";
    }
    std::cout<<'\n';
}

Kokkos::finalize();

}

