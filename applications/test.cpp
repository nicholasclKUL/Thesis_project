#include <Thesis/para-panoc.hpp>
#include <alpaqa/inner/directions/panoc/lbfgs.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#include <alpaqa/inner/inner-solve-options.hpp>

#include <iostream>

int main(int argc, char* argv[]){

    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
 
    // Problem specification
    // minimize  ½ xᵀQx
    //  s.t.     Ax ≤ b
    struct Problem : alpaqa::BoxConstrProblem<config_t> {
        mat Q{n, n};
        mat A{m, n};
        vec b{m};
        mutable vec Qx{n};
 
        Problem() : alpaqa::BoxConstrProblem<config_t>{2, 1} {
            // Initialize problem matrices
            Q << 3, -1, -1, 3;
            A << 2, 1;
            b << -1;
 
            // Specify the bounds
            C.lowerbound = vec::Constant(n, -alpaqa::inf<config_t>);
            C.upperbound = vec::Constant(n, +alpaqa::inf<config_t>);
            D.lowerbound = vec::Constant(m, -alpaqa::inf<config_t>);
            D.upperbound = b;
        }
 
        // Evaluate the cost
        real_t eval_f(crvec x) const {
            Qx.noalias() = Q * x;
            return 0.5 * x.dot(Qx);
        }
        // Evaluat the gradient of the cost
        void eval_grad_f(crvec x, rvec gr) const { gr.noalias() = Q * x; }
        // Evaluate the constraints
        void eval_g(crvec x, rvec g) const { g.noalias() = A * x; }
        // Evaluate a matrix-vector product with the gradient of the constraints
        void eval_grad_g_prod(crvec x, crvec y, rvec gr) const {
            (void)x;
            gr.noalias() = A.transpose() * y;
        }
    };
 
    Problem problem;
 
    // Wrap the problem to count the function evaluations
    auto counted_problem = alpaqa::problem_with_counters_ref(problem);
 
    // Define the solvers to use
    using Accelerator = alpaqa::LBFGSDirection<config_t>;
    using InnerSolver = alpaqa::ParaPANOCSolver<Accelerator>;

    // Initialize parameters for ParaPANOC
    InnerSolver::Params panocparam;

    // Inatialize parameters for L-BFGS
    Accelerator::LBFGSParams lbfgsparam;   

    InnerSolver solver{
        panocparam, lbfgsparam
    };

        // Initial guess
    vec x(2);
    x << 2, 2; // decision variables
    vec y(1);
    y << 1; // Lagrange multipliers
    vec Σ(1);
    Σ << 0.1; // Penalty 
    vec err(1);
    err << 0; //error in z
    alpaqa::InnerSolveOptions<config_t> sol_opts;
 
    // Solve the problem
    auto stats = solver(counted_problem, sol_opts , x, y, Σ, err);
    // y and x have been overwritten by the solution

    // Print the results
    std::cout << '\n' << *counted_problem.evaluations << '\n';
    std::cout << "status: " << stats.status << '\n'
              << "f = " << problem.eval_f(x) << '\n'
              << "inner iterations: " << stats.iterations << '\n'
              << "ε = " << stats.ε << '\n'
              << "elapsed time:     "
              << std::chrono::duration<double>{stats.elapsed_time}.count()
              << " s" << '\n'
              << "x = " << x.transpose() << '\n'
              << "y = " << y.transpose() << '\n'
              << "avg τ = " << (stats.sum_τ / stats.count_τ) << '\n'
              << "L-BFGS rejected = " << stats.lbfgs_rejected << '\n'
              << "L-BFGS failures = " << stats.lbfgs_failures << '\n'
              << "Line search failures = " << stats.linesearch_failures
              << '\n'
              << std::endl;

}