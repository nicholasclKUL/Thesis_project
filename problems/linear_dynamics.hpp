#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/util/print.hpp>
#include <iomanip>
#include <iostream>

struct LinearOCP {

 USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    using Box = alpaqa::Box<config_t>;

    length_t N = 5,      ///< Horizon length
        nu     = 2,       ///< Number of inputs
        nx     = 3,       ///< Number of states
        nh     = nu + nx, ///< Number of stage outputs
        nh_N   = nx,      ///< Number of terminal outputs
        nc     = 0,       ///< Number of stage constraints
        nc_N   = 0;       ///< Number of terminal constraints

    mat A, B;

    LinearOCP() : A(nx, nx), B(nx, nu) {
        A.setIdentity();
        B.setOnes();
        A *= 2;
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
    void get_D([[maybe_unused]] Box &D) const {       
        alpaqa::ScopedMallocAllower ma;
        D.lowerbound.setConstant(-alpaqa::inf<config_t>);
        D.upperbound.setConstant(+alpaqa::inf<config_t>);}

    void get_D_N([[maybe_unused]] Box &D) const {}

    void get_x_init(rvec x_init) const { x_init.setConstant(1.); }

    void eval_f([[maybe_unused]] index_t timestep, 
                crvec x, crvec u, rvec fxu) const {
        alpaqa::ScopedMallocAllower ma;
        fxu.noalias() = A * x + B * u;
    }
    void eval_jac_f([[maybe_unused]] index_t timestep,
                    [[maybe_unused]] crvec x,
                    [[maybe_unused]] crvec u, rmat J_fxu) const {
        alpaqa::ScopedMallocAllower ma;
        J_fxu.leftCols(nx).noalias()  = A;
        J_fxu.rightCols(nu).noalias() = B;
    }
    void eval_grad_f_prod([[maybe_unused]] index_t timestep, 
                          [[maybe_unused]] crvec x,
                          [[maybe_unused]] crvec u, crvec p,
                          rvec grad_fxu_p) const {
        alpaqa::ScopedMallocAllower ma;
        grad_fxu_p.topRows(nx).noalias()    = A.transpose() * p;
        grad_fxu_p.bottomRows(nu).noalias() = B.transpose() * p;
    }
    void eval_h([[maybe_unused]] index_t timestep, crvec x, crvec u, rvec h) const {
        alpaqa::ScopedMallocAllower ma;
        h.topRows(nx)    = x;
        h.bottomRows(nu) = u;
    }
    void eval_h_N(crvec x, rvec h) const { h = x; }

    [[nodiscard]] real_t eval_l([[maybe_unused]] index_t timestep, crvec h) const {
        alpaqa::ScopedMallocAllower ma;
        return 0.5 * h.squaredNorm();
    }
    [[nodiscard]] real_t eval_l_N(crvec h) const {
        alpaqa::ScopedMallocAllower ma;
        return 5. * h.squaredNorm();
    }
    void eval_qr([[maybe_unused]] index_t timestep, 
                 [[maybe_unused]] crvec xu, crvec h, rvec qr) const {
        alpaqa::ScopedMallocAllower ma;
        auto Jh_xu    = mat::Identity(nx + nu, nx + nu);
        auto &&grad_l = h;
        qr            = Jh_xu.transpose() * grad_l;
    }
    void eval_q_N([[maybe_unused]] crvec x, crvec h, rvec q) const {
        alpaqa::ScopedMallocAllower ma;
        auto Jh_x     = mat::Identity(nx, nx);
        auto &&grad_l = 10 * h;
        q             = Jh_x.transpose() * grad_l;
    }
    void eval_add_Q([[maybe_unused]] index_t timestep, 
                    [[maybe_unused]] crvec xu, 
                    [[maybe_unused]] crvec h, rmat Q) const {
        alpaqa::ScopedMallocAllower ma;
        auto Jh_xu    = mat::Identity(nx + nu, nx + nu);
        Q.noalias()   = Jh_xu.transpose() * Jh_xu;
        // Q += mat::Identity(nx,nx);
    }
    void eval_add_Q_N([[maybe_unused]] crvec x,
                      [[maybe_unused]] crvec h, rmat Q) const {
        alpaqa::ScopedMallocAllower ma;
        auto Jh_x     = mat::Identity(nx, nx);
        Q.noalias()   = 10 * (Jh_x.transpose() * Jh_x);
        // Q += 10 * mat::Identity(nx,nx);
    }
    void eval_add_R_masked([[maybe_unused]] index_t timestep,
                           [[maybe_unused]] crvec xu, 
                           [[maybe_unused]] crvec h, crindexvec mask,
                           rmat R, 
                           [[maybe_unused]] rvec work) const {
        alpaqa::ScopedMallocAllower ma;
        const auto n = mask.size();
        R.noalias() += mat::Identity(n, n);
    }
    void eval_add_S_masked([[maybe_unused]] index_t timestep, 
                           [[maybe_unused]] crvec xu, 
                           [[maybe_unused]] crvec h, crindexvec mask,
                           rmat S, 
                           [[maybe_unused]] rvec work) const {
        // Mixed derivatives are zero
        // S.noalias() += (...);
    }
    void eval_add_R_prod_masked([[maybe_unused]] index_t timestep,
                                [[maybe_unused]] crvec xu, 
                                [[maybe_unused]] crvec h,
                                crindexvec mask_J, crindexvec mask_K, crvec v,
                                rvec out, 
                                [[maybe_unused]] rvec work) const {
        // The following has no effect because R is diagonal, and J ∩ K = ∅
        alpaqa::ScopedMallocAllower ma;
        auto R = mat::Identity(nu, nu);
        out.noalias() += R(mask_J, mask_K) * v(mask_K);
    }
    void eval_add_S_prod_masked([[maybe_unused]] index_t timestep,
                                [[maybe_unused]] crvec xu, 
                                [[maybe_unused]] crvec h,
                                [[maybe_unused]] crindexvec mask_K, 
                                [[maybe_unused]] crvec v, 
                                [[maybe_unused]] rvec out,
                                [[maybe_unused]] rvec work) const {
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

    void eval_proj_multipliers([[maybe_unused]] rvec y, 
                               [[maybe_unused]] real_t M) const {}

    void eval_proj_diff_g([[maybe_unused]] crvec z, 
                          [[maybe_unused]] rvec p) const { p.setZero(); }

    void check() const {
        // You could do some sanity checks here
    }
};

// ___________________________________________________________________________________________________________________________________________________________________________ //
// Linear Dynamics problem supported by AD

#include <thesis/auto-diff-sacado-helpers.hpp>

#include <time-discretization.hpp>

template <typename X, typename F, typename Params>
void dynamics(index_t k, X &xu, F &fxu, const Params &params){
    auto nx = fxu.size();
    auto nu = xu.size()-nx;
    for (size_t i = 0; i < nx; ++i){
        if (i < nu){
            fxu(i) = params.A(i,i)*xu(i) + params.B(i,i)*xu(i+nx);
        }
        else{
            fxu(i) = params.A(i,i)*xu(i);
        }
    }
};

const int p_N = 5, p_nx = 3, p_nu = 2, p_nh = 5;

using ADobj = AD<p_nx+p_nu,p_nx,p_nh>;

struct LinearOCPAD {

    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    using Box = alpaqa::Box<config_t>;

    struct Params {
        length_t N = p_N,       ///< Horizon length
            nu     = p_nu,      ///< Number of inputs
            nx     = p_nx,      ///< Number of states
            nh     = p_nh,      ///< Number of stage outputs
            nh_N   = p_nx,      ///< Number of terminal outputs
            nc     = 0,         ///< Number of stage constraints
            nc_N   = 0;         ///< Number of terminal constraints

        real_t T = 1.,
              Ts = T/real_t(N); 

        mat A, B;
    };

    Params params;

    mat A, B;

    ADobj ad_obj[p_N]; 

    LinearOCPAD() : A(params.nx, params.nx), B(params.nx, params.nu) {
        A.setIdentity();
        B.setIdentity();
        A *= -2;

        params.A = A; params.B = B;
    }

    [[nodiscard]] length_t get_N() const { return params.N; }
    [[nodiscard]] length_t get_nu() const { return params.nu; }
    [[nodiscard]] length_t get_nx() const { return params.nx; }
    [[nodiscard]] length_t get_nh() const { return params.nh; }
    [[nodiscard]] length_t get_nh_N() const { return params.nh_N; }
    [[nodiscard]] length_t get_nc() const { return params.nc; }
    [[nodiscard]] length_t get_nc_N() const { return params.nc_N; }

    void get_U(Box &U) const {
        alpaqa::ScopedMallocAllower ma;
        U.lowerbound.setConstant(-alpaqa::inf<config_t>);
        U.upperbound.setConstant(+alpaqa::inf<config_t>);
    }
    void get_D([[maybe_unused]] Box &D) const {       
        alpaqa::ScopedMallocAllower ma;
        D.lowerbound.setConstant(-alpaqa::inf<config_t>);
        D.upperbound.setConstant(+alpaqa::inf<config_t>);}

    void get_D_N([[maybe_unused]] Box &D) const {}

    void get_x_init(rvec x_init) const { x_init.setConstant(1.); }

    void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const { 
        alpaqa::ScopedMallocAllower ma;
        vec xu(params.nx+params.nu); xu << x, u;
        // fe(timestep, xu, fxu, params);
        rk4(timestep, xu, fxu, params);
        // std::cout<<fxu.transpose()<<'\n'<<xu.transpose()<<std::endl;
    } 
    void eval_jac_f(index_t timestep, crvec x, crvec u, rmat Jfxu) const {
        alpaqa::ScopedMallocAllower ma;
        assign_values_xu<p_nx+p_nu,p_nx>(x, u, ad_obj[timestep]);
        // fe(timestep, xu, fxu, params);
        rk4(timestep, ad_obj[timestep].xu_fad, ad_obj[timestep].fxu_fad, params, p_nx+p_nu, p_nx);
        assign_values<p_nx+p_nu,p_nx>(Jfxu, ad_obj[timestep]);
    }
    void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p,
                          rvec grad_fxu_p) const {
        alpaqa::ScopedMallocAllower ma;
        assign_values_xu<p_nx+p_nu,p_nx>(x, u, ad_obj[timestep]);
        // fe(timestep, xu, fxu, params);
        rk4(timestep, ad_obj[timestep].xu_fad, ad_obj[timestep].fxu_fad, params, p_nx+p_nu, p_nx);
        assign_values<p_nx+p_nu,p_nx> (grad_fxu_p, p, ad_obj[timestep]);
    }
    void eval_h([[maybe_unused]] index_t timestep, crvec x, crvec u, rvec h) const {
        alpaqa::ScopedMallocAllower ma;
        h.topRows(params.nx)    = x;
        h.bottomRows(params.nu) = u;
    }
    void eval_h_N(crvec x, rvec h) const { h = x; }

    [[nodiscard]] real_t eval_l([[maybe_unused]] index_t timestep, crvec h) const {
        alpaqa::ScopedMallocAllower ma;
        // std::cout<<h.transpose()<<std::endl;
        return 0.5 * h.squaredNorm();
    }
    [[nodiscard]] real_t eval_l_N(crvec h) const {
        alpaqa::ScopedMallocAllower ma;
        // std::cout<<h.transpose()<<std::endl;
        return 5. * h.squaredNorm();
    }
    void eval_qr([[maybe_unused]] index_t timestep, 
                 [[maybe_unused]] crvec xu, crvec h, rvec qr) const {
        alpaqa::ScopedMallocAllower ma;
        auto Jh_xu    = mat::Identity(params.nx + params.nu, params.nx + params.nu);
        auto &&grad_l = h;
        qr            = Jh_xu.transpose() * grad_l;
    }
    void eval_q_N([[maybe_unused]] crvec x, crvec h, rvec q) const {
        alpaqa::ScopedMallocAllower ma;
        auto Jh_x     = mat::Identity(params.nx, params.nx);
        auto &&grad_l = 10 * h;
        q             = Jh_x.transpose() * grad_l;
    }
    void eval_add_Q([[maybe_unused]] index_t timestep, 
                    [[maybe_unused]] crvec xu, 
                    [[maybe_unused]] crvec h, rmat Q) const {
        alpaqa::ScopedMallocAllower ma;
        auto Jh_xu    = mat::Identity(params.nx + params.nu, params.nx + params.nu);
        Q.noalias()   = Jh_xu.transpose() * Jh_xu;
        // Q += mat::Identity(nx,nx);
    }
    void eval_add_Q_N([[maybe_unused]] crvec x,
                      [[maybe_unused]] crvec h, rmat Q) const {
        alpaqa::ScopedMallocAllower ma;
        auto Jh_x     = mat::Identity(params.nx, params.nx);
        Q.noalias()   = 10 * (Jh_x.transpose() * Jh_x);
        // Q += 10 * mat::Identity(nx,nx);
    }
    void eval_add_R_masked([[maybe_unused]] index_t timestep,
                           [[maybe_unused]] crvec xu, 
                           [[maybe_unused]] crvec h, crindexvec mask,
                           rmat R, 
                           [[maybe_unused]] rvec work) const {
        alpaqa::ScopedMallocAllower ma;
        const auto n = mask.size();
        R.noalias() += mat::Identity(n, n);
    }
    void eval_add_S_masked([[maybe_unused]] index_t timestep, 
                           [[maybe_unused]] crvec xu, 
                           [[maybe_unused]] crvec h, crindexvec mask,
                           rmat S, 
                           [[maybe_unused]] rvec work) const {
        // Mixed derivatives are zero
        // S.noalias() += (...);
    }
    void eval_add_R_prod_masked([[maybe_unused]] index_t timestep,
                                [[maybe_unused]] crvec xu, 
                                [[maybe_unused]] crvec h,
                                crindexvec mask_J, crindexvec mask_K, crvec v,
                                rvec out, 
                                [[maybe_unused]] rvec work) const {
        // The following has no effect because R is diagonal, and J ∩ K = ∅
        alpaqa::ScopedMallocAllower ma;
        auto R = mat::Identity(params.nu, params.nu);
        out.noalias() += R(mask_J, mask_K) * v(mask_K);
    }
    void eval_add_S_prod_masked([[maybe_unused]] index_t timestep,
                                [[maybe_unused]] crvec xu, 
                                [[maybe_unused]] crvec h,
                                [[maybe_unused]] crindexvec mask_K, 
                                [[maybe_unused]] crvec v, 
                                [[maybe_unused]] rvec out,
                                [[maybe_unused]] rvec work) const {
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

    void eval_proj_multipliers([[maybe_unused]] rvec y, 
                               [[maybe_unused]] real_t M) const {}

    void eval_proj_diff_g([[maybe_unused]] crvec z, 
                          [[maybe_unused]] rvec p) const { p.setZero(); }

    void check() const {
        // You could do some sanity checks here
    }
};
