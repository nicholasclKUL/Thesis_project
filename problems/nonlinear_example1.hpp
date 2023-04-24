#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/print.hpp>
#include <iomanip>
#include <iostream>

namespace alpaqa{

struct NonlinearOCP1 {

    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    using Box = alpaqa::Box<config_t>;

    length_t N = 5,      ///< Horizon length
        nu     = 1,       ///< Number of inputs
        nx     = 100,       ///< Number of states
        nh     = nu + nx, ///< Number of stage outputs
        nh_N   = nx,      ///< Number of terminal outputs
        nc     = 0,       ///< Number of stage constraints
        nc_N   = 0;       ///< Number of terminal constraints

    mat A, B;

    NonlinearOCP1() : A(nx, nx), B(nx, nu) {
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
        U.lowerbound.setConstant(-alpaqa::inf<config_t>);
        U.upperbound.setConstant(+alpaqa::inf<config_t>);
    }
    void get_D(Box &D) const {       
        alpaqa::ScopedMallocAllower ma;
        D.lowerbound.setConstant(-alpaqa::inf<config_t>);
        D.upperbound.setConstant(+alpaqa::inf<config_t>);}

    void get_D_N(Box &D) const {}

    void get_x_init(rvec x_init) const { x_init.setConstant(1.); }

    void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const {
        alpaqa::ScopedMallocAllower ma;
        fxu.noalias() = x.asDiagonal() * x + B * u;
    }
    void eval_jac_f(index_t timestep, crvec x, crvec u, rmat J_fxu) const {
        alpaqa::ScopedMallocAllower ma;
        J_fxu.leftCols(nx).noalias()  = 2 * (A*x.asDiagonal());
        J_fxu.rightCols(nu).noalias() = B;
    }
    void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p,
                          rvec grad_fxu_p) const {
        alpaqa::ScopedMallocAllower ma;
        grad_fxu_p.topRows(nx).noalias()    = 2 * (A*x.asDiagonal()) * p;
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

    void eval_proj_multipliers(rvec y, real_t M) const {}

    void eval_proj_diff_g(crvec z, rvec p) const { p.setZero(); }

    void check() const {
        // You could do some sanity checks here
    }
};

}