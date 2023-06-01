#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/outer/alm.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/print.hpp>

#include <thesis/para-alm.hpp>
#include <thesis/para-panoc.hpp>

#include <iomanip>
#include <iostream>
#include <fstream>

namespace printing {

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

using Stats_inner       = alpaqa::PANOCOCPStats<config_t>;
using Stats_outer_ms    = alpaqa::ParaALMSolver<alpaqa::ParaPANOCSolver<alpaqa::DefaultConfig>>::Stats;
using Stats_outer_ss    = alpaqa::ALMSolver<alpaqa::PANOCOCPSolver<alpaqa::DefaultConfig>>::Stats;

template <typename P>
void print_solution (P &problem, crvec xu) {
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

template <typename P>
void print_solution_ss (P &problem, crvec u) {
        for (length_t i = 0; i < problem.get_N()+1; ++i){
        std::cout<<"--Stage "<<i<<":"<<'\n';
        if (i < problem.get_N()){
                     std::cout<<"u = ["<<u.segment(i*(problem.get_nu()),problem.get_nu()).transpose()<<"]"<<'\n';
        }      
    }
}

void print_stats_inner (Stats_inner &stats, crvec e){
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

void print_stats_outer (Stats_outer_ms stats, crvec g){
    std::cout<<'\n'<<"--[FINAL REPORT]--"<<'\n'<<std::endl;
    std::cout   << "status: " << stats.status << '\n'
                << "inner iterations: " << stats.inner.iterations << '\n'
                << "outer iterations: " << stats.outer_iterations << '\n'
                << "ε = " << stats.ε << '\n'
                << "δ = " << stats.δ << '\n'
                << "Δg = " << g.norm() << '\n'
                << "elapsed time:     "
                << std::chrono::duration<double>{stats.elapsed_time}.count()
                << " s" << '\n'
                << "avg τ = " << (stats.inner.sum_τ / stats.inner.count_τ) << '\n'
                << "L-BFGS rejected = " << stats.inner.lbfgs_rejected << '\n'
                << "L-BFGS failures = " << stats.inner.lbfgs_failures << '\n'
                << "Line search failures = " << stats.inner.linesearch_failures
                << '\n'
                << std::endl;
}

void print_stats_outer (Stats_outer_ms stats){
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
}

void print_stats_outer (Stats_outer_ss stats){
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
}

template <typename P>
void output_file (P &problem, std::string &problem_name, crvec xu, index_t Ts){
    auto N = problem.get_N(),
        nu = problem.get_nu(),
        nx = problem.get_nx();
    // create 1st option of possible file name
    std::string filename = problem_name + "_" + std::to_string(problem.get_N())
                        + "_" + std::to_string(Ts) + "_" + std::to_string(0) + ".csv";
    // search for available names
    FILE *file_check;
    unsigned int k = 1;  
    const char* filename_ = filename.c_str();
    file_check = std::fopen(filename_,"r");
    while ((file_check != NULL) && (k <= 1000)){
        filename = problem_name + "_" + std::to_string(problem.get_N())
                + "_" + std::to_string(Ts) + "_" + std::to_string(k) + ".csv";
        filename_ = filename.c_str();
        file_check = std::fopen(filename_,"r");
        ++k;
    }
    // open/create file
    std::ofstream myfile(filename_);
    // fill in header
    for (size_t i = 0; i < nx; ++i){
        myfile<<"x_"<<i<<","<<std::setw(14);
    }
    for (size_t i = 0; i < nu; ++i){
        myfile<<"u_"<<i<<","<<std::setw(14);
    }
    myfile<<std::endl;
    // fill in data for 
    for (size_t i = 0; i < N; ++i){
        for (size_t j = 0; j < nx+nu; ++j){
            myfile<<std::scientific<<xu((i*(nx+nu))+j)<<","<<std::setw(14);   
        }  
        myfile<<'\n';
    }
    myfile<<std::endl;   
    // close file
    myfile.close();
}

}