#pragma once
#include <alpaqa/casadi/CasADiControlProblem.hpp>
#include <problem-loader-export.h>

struct Problem {
    USING_ALPAQA_CONFIG(alpaqa::EigenConfigd);
    alpaqa::CasADiProblem<config_t> problem;
    vec initial_guess = vec(problem.get_n());
};

PROBLEM_LOADER_EXPORT Problem load_problem(const char *directory);