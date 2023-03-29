/**
 * @file
 * Loads a compiled CasADi problem and problem values from a .tsv file.
 */

#include <alpaqa/casadi/CasADiProblem.hpp>
#include <alpaqa/casadi/CasADiControlProblem.hpp>
#include <problem-loader.hpp>

#include <charconv>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <mutex>

namespace fs = std::filesystem;

namespace {

/// Hacky way to read a vector from a .tsv file
void read_vector(std::istream &is, Problem::vec &v) {
    std::string s;
    s.reserve(32);
    for (auto &vv : v)
        if (!(is >> s))
            throw std::runtime_error("read_vector extraction failed");
        else if (std::from_chars(&*s.begin(), &*s.end(), vv).ec != std::errc{})
            throw std::runtime_error("read_vector conversion failed");
    if (is.get() != '\n' && is)
        throw std::runtime_error("read_vector line not finished");
}

Problem load_problem(const fs::path &full_path, bool second_order) {
    static std::mutex mtx;
    std::unique_lock lck{mtx};
    // Load CasADi problem and allocate workspaces
    Problem problem{.problem = {full_path, 0, 0, 0, second_order}};
    lck.unlock();

    // Load numeric data
    std::ifstream bounds{fs::path{full_path}.replace_extension("tsv")};
    if (!bounds)
        throw std::runtime_error("Failed to open bounds file");
    read_vector(bounds, problem.problem.C.lowerbound);
    read_vector(bounds, problem.problem.C.upperbound);
    read_vector(bounds, problem.problem.D.lowerbound);
    read_vector(bounds, problem.problem.D.upperbound);
    read_vector(bounds, problem.problem.param);
    read_vector(bounds, problem.initial_guess);
    return problem;
}

ProblemOCP load_problem_ocp(const fs::path &full_path) {
    static std::mutex mtx;
    std::unique_lock lck{mtx};
    // Load CasADi problem and allocate workspaces
    ProblemOCP problem{.problem = {full_path, 0}};
    lck.unlock();

    // Load numeric data
    std::ifstream bounds{fs::path{full_path}.replace_extension("tsv")};
    if (!bounds)
        throw std::runtime_error("Failed to open bounds file");
    read_vector(bounds, problem.problem.D.lowerbound);
    read_vector(bounds, problem.problem.D.upperbound);
    read_vector(bounds, problem.problem.param);
    read_vector(bounds, problem.initial_guess);
    return problem;
}

}   // namespace

Problem load_problem(const char *directory) {
    return load_problem(fs::path{directory} / PROBLEM_DLL,
                        static_cast<bool>(PROBLEM_SECOND_ORDER));
}

ProblemOCP load_problem_ocp(const char *directory) {
    return load_problem_ocp(fs::path{directory}/ PROBLEM_DLL);
}

