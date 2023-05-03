#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/print.hpp>

#include <iomanip>
#include <iostream>

struct MyProblem {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    using Box = alpaqa::Box<config_t>;

    length_t N = 5,      ///< Horizon length
        nu     = 1,       ///< Number of inputs
        nx     = 5,       ///< Number of states
        nh     = nu + nx, ///< Number of stage outputs
        nh_N   = nx,      ///< Number of terminal outputs
        nc     = 0,       ///< Number of stage constraints
        nc_N   = 0;       ///< Number of terminal constraints

    mat A, B;

    MyProblem() : A(nx, nx), B(nx, nu) {
        A.setIdentity();
        B.setOnes();
    }

    [[nodiscard]] length_t get_N() const { return N; }
    [[nodiscard]] length_t get_nu() const { return nu; }
    [[nodiscard]] length_t get_nx() const { return nx; }
    [[nodiscard]] length_t get_nh() const { return nh; }
    [[nodiscard]] length_t get_nh_N() const { return nh_N; }
    [[nodiscard]] length_t get_nc() const { return nc; }
    [[nodiscard]] length_t get_nc_N() const { return nc_N; }

    void get_U(Box &U) const {
        alpaqa::ScopedMallocAllower ma;
        U.lowerbound.setConstant(-2);
        U.upperbound.setConstant(+2);
    }
    void get_D(Box &D) const {       
        alpaqa::ScopedMallocAllower ma;
        D.lowerbound.setConstant(0);
        D.upperbound.setConstant(0);}
    void get_D_N(Box &D) const {}

    void get_x_init(rvec x_init) const { x_init.setConstant(1.); }

    void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const {
        alpaqa::ScopedMallocAllower ma;
        fxu.noalias() = A * x + B * u;
    }
    void eval_jac_f(index_t timestep, crvec x, crvec u, rmat J_fxu) const {
        alpaqa::ScopedMallocAllower ma;
        J_fxu.leftCols(nx).noalias()  = A;
        J_fxu.rightCols(nu).noalias() = B;
    }
    void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p,
                          rvec grad_fxu_p) const {
        alpaqa::ScopedMallocAllower ma;
        grad_fxu_p.topRows(nx).noalias()    = A.transpose() * p;
        grad_fxu_p.bottomRows(nu).noalias() = B.transpose() * p;
    }
    void eval_h(index_t timestep, crvec x, crvec u, rvec h) const {
        alpaqa::ScopedMallocAllower ma;
        h.topRows(nx)    = x;
        h.bottomRows(nu) = u;
    }
    void eval_h_N(crvec x, rvec h) const { h = x; }

    [[nodiscard]] real_t eval_l(index_t timestep, crvec h) const {
        alpaqa::ScopedMallocAllower ma;
        return 0.5 * h.squaredNorm();
    }
    [[nodiscard]] real_t eval_l_N(crvec h) const {
        alpaqa::ScopedMallocAllower ma;
        return 5. * h.squaredNorm();
    }
    void eval_qr(index_t timestep, crvec xu, crvec h, rvec qr) const {
        alpaqa::ScopedMallocAllower ma;
        auto Jh_xu    = mat::Identity(nx + nu, nx + nu);
        auto &&grad_l = h;
        qr            = Jh_xu.transpose() * grad_l;
    }
    void eval_q_N(crvec x, crvec h, rvec q) const {
        alpaqa::ScopedMallocAllower ma;
        auto Jh_x     = mat::Identity(nx, nx);
        auto &&grad_l = 10 * h;
        q             = Jh_x.transpose() * grad_l;
    }
    void eval_add_Q(index_t timestep, crvec xu, crvec h, rmat Q) const {
        Q += mat::Identity(nx, nx);
    }
    void eval_add_Q_N(crvec x, crvec h, rmat Q) const {
        alpaqa::ScopedMallocAllower ma;
        Q += 10 * mat::Identity(nx, nx);
    }
    void eval_add_R_masked(index_t timestep, crvec xu, crvec h, crindexvec mask,
                           rmat R, rvec work) const {
        alpaqa::ScopedMallocAllower ma;
        const auto n = mask.size();
        R.noalias() += mat::Identity(n, n);
    }
    void eval_add_S_masked(index_t timestep, crvec xu, crvec h, crindexvec mask,
                           rmat S, rvec work) const {
        // Mixed derivatives are zero
        // S.noalias() += (...);
    }
    void eval_add_R_prod_masked(index_t timestep, crvec xu, crvec h,
                                crindexvec mask_J, crindexvec mask_K, crvec v,
                                rvec out, rvec work) const {
        // The following has no effect because R is diagonal, and J ∩ K = ∅
        alpaqa::ScopedMallocAllower ma;
        auto R = mat::Identity(nu, nu);
        out.noalias() += R(mask_J, mask_K) * v(mask_K);
    }
    void eval_add_S_prod_masked(index_t timestep, crvec xu, crvec h,
                                crindexvec mask_K, crvec v, rvec out,
                                rvec work) const {
        // Mixed derivatives are zero
        // using Eigen::indexing::all;
        // auto Sᵀ = (...);
        // out.noalias() += Sᵀ(all, mask_K) * v(mask_K);
    }
    [[nodiscard]] length_t get_R_work_size() const {
        // No workspace needed to share data between eval_add_R_masked and
        // eval_add_R_prod_masked
        return 0;
    }
    [[nodiscard]] length_t get_S_work_size() const { return 0; }

    void eval_proj_multipliers(rvec y, real_t M,
                               index_t penalty_alm_split) const {}

    void eval_proj_diff_g(crvec z, rvec p) const { p.setZero(); }

    void check() const {
        // You could do some sanity checks here
    }
};

int main() {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

    // Problem
#if 0
    // Option 1: simply use your own problem class directly, and let the solver
    //           convert it to a TypeErasedControlProblem for you.
    Problem problem;
#elif 0
    // Option 2: initialize a TypeErasedControlProblem using an instance of your
    // own problem class, so the solver doesn't have to create a copy of it.
    alpaqa::TypeErasedControlProblem<config_t> problem{Problem()};
#elif 0
    // Option 3: create a TypeErasedControlProblem from your own problem class
    //           directly, using in-place construction (most efficient).
    auto problem = alpaqa::TypeErasedControlProblem<config_t>::make<Problem>();
#else
    // Option 4: use your own problem class, and wrap it with evaluation
    //           counters (handy for benchmarking)
    auto problem = alpaqa::ocproblem_with_counters(MyProblem());
#endif

    // Problem dimensions
    const auto n = problem.get_N() * problem.get_nu(),
               m = problem.get_N() * problem.get_nc() + problem.get_nc_N();

    // Initial guess and other solver inputs
    vec u = vec::Zero(n); // Inputs (single shooting)
    vec y = vec::Zero(m); // Lagrange multipliers
    vec μ = vec::Ones(m); // Penalty factors
    vec e(m);             // Constraint violation

    // Solver
    alpaqa::PANOCOCPParams<config_t> params;
    params.stop_crit      = alpaqa::PANOCStopCrit::ProjGradUnitNorm;
    params.gn_interval    = 1; //for gn_interval != 1, program throws an exception
    params.print_interval = 1;
    alpaqa::PANOCOCPSolver<config_t> solver{params};

    // Solve
    auto stats = solver(problem, {.tolerance = 1e-8}, u, y, μ, e);

    // Print evaluation counters (if available)
    [](const auto &problem) {
        if constexpr (requires { problem.evaluations; })
            std::cout << '\n' << *problem.evaluations;
    }(problem);

    // Print statistics
    auto δ      = e.lpNorm<Eigen::Infinity>();
    auto time_s = std::chrono::duration<double>(stats.elapsed_time).count();
    std::cout << '\n'
              << "solver:  " << solver.get_name() << '\n'
              << "status:  " << stats.status << '\n'
              << "ψ = " << alpaqa::float_to_str(stats.final_ψ) << '\n'
              << "ε = " << alpaqa::float_to_str(stats.ε) << '\n'
              << "δ = " << alpaqa::float_to_str(δ) << '\n'
              << "time: " << alpaqa::float_to_str(time_s, 3) << " s\n"
              << "iter:      " << std::setw(6) << stats.iterations << '\n'
              << "line search backtrack: " << std::setw(6)
              << stats.linesearch_backtracks << '\n'
              << "step size backtrack:   " << std::setw(6)
              << stats.stepsize_backtracks << '\n'
              << "solution: ";
    // alpaqa::print_python(std::cout, u) << std::endl;
}