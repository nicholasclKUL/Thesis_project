#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/print.hpp>

namespace printing {

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

using Problem = alpaqa::TypeErasedControlProblem<config_t>;
using Stats   = alpaqa::PANOCOCPStats<config_t>;

void print_solution (Problem &problem, crvec xu) {
    for (length_t i = 0; i < problem.get_N()+1; ++i){
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

void print_stats (Stats &stats, crvec e){
    auto δ    = e.lpNorm<Eigen::Infinity>();
    auto time = std::chrono::duration<double>(stats.elapsed_time).count();
    std::cout << '\n'
              << "status:  " << stats.status << '\n'
              << "ψ = " << alpaqa::float_to_str(stats.final_ψ) << '\n'
              << "ε = " << alpaqa::float_to_str(stats.ε) << '\n'
              << "δ = " << alpaqa::float_to_str(δ) << '\n'
              << "time: " << alpaqa::float_to_str(time, 3) << " s\n"
              << "iter:      " << std::setw(6) << stats.iterations << '\n'
              << "line search backtrack: " << std::setw(6)
              << stats.linesearch_backtracks << '\n'
              << "step size backtrack:   " << std::setw(6)
              << stats.stepsize_backtracks << '\n'
              << "solution: "<<'\n';
}


}