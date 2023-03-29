#include <alpaqa/config/config.hpp>
#include <alpaqa/structured-panoc-alm.hpp>
#include <alpaqa/util/print.hpp>
#include <problem-loader.hpp>
USING_ALPAQA_CONFIG(Problem::config_t);

#include <cassert>
#include <filesystem>
#include <iostream>
namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    assert(argc >= 1);
    fs::path problem_path = fs::canonical(fs::path(argv[0])).parent_path();
    auto problem          = load_problem(problem_path.c_str());


    // Define the solvers to use
    using Accelerator = alpaqa::StructuredLBFGSDirection<config_t>;
    using InnerSolver = alpaqa::PANOCSolver<Accelerator>;
    using OuterSolver = alpaqa::ALMSolver<InnerSolver>;

    // Settings for the outer augmented Lagrangian method
    OuterSolver::Params almparam;
    almparam.ε              = 1e-8; // tolerance
    almparam.δ              = 1e-8;
    almparam.print_interval = 1;

    // Settings for the inner PANOC solver
    InnerSolver::Params panocparam;
    panocparam.max_iter       = 500;
    panocparam.print_interval = 25;
    panocparam.stop_crit      = alpaqa::PANOCStopCrit::ProjGradUnitNorm2;
    // Settings for the L-BFGS algorithm used by PANOC
    Accelerator::LBFGSParams lbfgsparam;
    lbfgsparam.memory = 40;

    // Create an ALM solver using PANOC as inner solver
    OuterSolver solver{
        almparam,                 // params for outer solver
        {panocparam, lbfgsparam}, // inner solver
    };

    // Initial guess
    vec x = problem.initial_guess;
    vec y = vec::Zero(problem.problem.get_m());

    // Wrap the problem to count the function evaluations
    auto counted_problem = alpaqa::problem_with_counters_ref(problem.problem);

    // Solve the problem
    auto stats = solver(counted_problem, x, y);

    // Print the results
    alpaqa::print_python(std::cout << "x = ", x);
    alpaqa::print_python(std::cout << "y = ", y);
    std::cout << "\nEvaluations:\n\n" << *counted_problem.evaluations << '\n';
    std::cout << "status: " << stats.status << '\n'
              << "f = " << problem.problem.eval_f(x) << '\n'
              << "ε = " << stats.ε << '\n'
              << "δ = " << stats.δ << '\n'
              << "inner: " << stats.inner.iterations << '\n'
              << "outer: " << stats.outer_iterations << '\n'
              << "time:  "
              << 1e3 * std::chrono::duration<double>(stats.elapsed_time).count()
              << " ms\n"
              << std::endl;
}
