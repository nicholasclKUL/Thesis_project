#pragma once
#include <alpaqa/casadi/CasADiProblem.hpp>
#include <alpaqa/casadi/CasADiControlProblem.hpp>
#include <problem-loader-export.h>

struct Problem {
    USING_ALPAQA_CONFIG(alpaqa::EigenConfigd);
    alpaqa::CasADiProblem<config_t> problem;
    vec initial_guess = vec(problem.get_n());
};

struct ProblemOCP {
    USING_ALPAQA_CONFIG(alpaqa::EigenConfigd);
    alpaqa::CasADiControlProblem<config_t> problem;
    vec initial_guess = vec(problem.get_nx());
};

PROBLEM_LOADER_EXPORT Problem load_problem(const char *directory);
PROBLEM_LOADER_EXPORT ProblemOCP load_problem_ocp(const char *directory);