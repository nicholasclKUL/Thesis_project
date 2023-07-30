#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/print.hpp>
#include <iomanip>
#include <iostream>
#include <random>
#include <chrono>

#include <thesis/auto-diff-sacado-helpers.hpp>

#include <time-discretization.hpp>

struct NonlinearOCP1 {

    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    using Box = alpaqa::Box<config_t>;

    length_t N = 50,      ///< Horizon length
        nu     = 1,       ///< Number of inputs
        nx     = 50,       ///< Number of states
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

    void get_x_init(rvec x_init) const { x_init.setConstant(0.1); }

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

// ______________________________________________________________________________________________________________________________________________________________________ //
// Nonlinear Dynamics System no. 2 - inspired in the nonlinear system of "Scientific Software" assignment 4

template <typename X, typename F, typename Params>
void dynamics(index_t k, X &xu, F &fxu, const Params &params){
  for (size_t i = 0; i < params.nx; ++i) {
      if (i == 0){
          fxu(i) = std::pow(xu(params.nx),2);
      }
      else {
          fxu(i) = - 10 * std::pow(xu(i-1)/(real_t(i)),3) + xu(params.nx)*xu(i);
      } 
  }
};

const int p_N = 16, p_nx = 20, p_nu = 1, p_nh = 21;

using ADobj = AD<p_nx+p_nu,p_nx,p_nh>;

struct NonlinearOCP1b {

    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    using Box = alpaqa::Box<config_t>;

    struct Params {
        
        length_t N = p_N,             ///< Horizon length
            nu     = p_nu,            ///< Number of inputs
            nx     = p_nx,            ///< Number of states
            nh     = p_nu + p_nx,     ///< Number of stage outputs
            nh_N   = p_nx,            ///< Number of terminal outputs
            nc     = 0,               ///< Number of stage constraints
            nc_N   = 0;               ///< Number of terminal constraints

        real_t   T = 1.,
                Ts = T/real_t(N);
    
    };

unsigned long int n_seed = 11;

    vec Q, R, QR; 

    Params params;

    mat A, B;

    ADobj adobj[p_N]; 

    NonlinearOCP1b() : A(params.nx, params.nx), B(params.nx, params.nu) {
        A.setIdentity();
        B.setOnes();
        Q.setConstant(1.);
        R.setConstant(1.);
        QR << Q, R;
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
        U.lowerbound.setConstant(-1);
        U.upperbound.setConstant(+1);
    }
    void get_D(Box &D) const {       
        alpaqa::ScopedMallocAllower ma;
        D.lowerbound.setConstant(-1);
        D.upperbound.setConstant(+1);}

    void get_D_N(Box &D) const {}

    void get_x_init(rvec x_init) const { 
        
        real_t lower_bound = -1;
        real_t upper_bound = 1;
 
        std::uniform_real_distribution<real_t> unif(lower_bound,
                                           upper_bound);

        std::default_random_engine re;
        re.seed(n_seed);

        for (size_t i = 0; i < x_init.size(); ++i){
            x_init(i) = unif(re);
        }

    }

    void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const { 
        alpaqa::ScopedMallocAllower ma;
        vec xu(params.nx+params.nu); xu << x, u;
        fe(timestep, xu, fxu, params);
        // rk4(timestep, xu, fxu, params);
    } 
    void eval_jac_f(index_t timestep, crvec x, crvec u, rmat Jfxu) const {
        alpaqa::ScopedMallocAllower ma;
        assign_values_xu<p_nx+p_nu,p_nx>(x, u, adobj[timestep]);
        fe(timestep, adobj[timestep].xu_fad, adobj[timestep].fxu_fad, params);
        // rk4(timestep, adobj[timestep].xu_fad, adobj[timestep].fxu_fad, params, p_nx+p_nu, p_nx);
        assign_values<p_nx+p_nu,p_nx>(Jfxu, adobj[timestep]);
    }
    void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p,
                          rvec grad_fxu_p) const {
        alpaqa::ScopedMallocAllower ma;
        assign_values_xu<p_nx+p_nu,p_nx>(x, u, adobj[timestep]);
        fe(timestep, adobj[timestep].xu_fad, adobj[timestep].fxu_fad, params);
        // rk4(timestep, adobj[timestep].xu_fad, adobj[timestep].fxu_fad, params, p_nx+p_nu, p_nx);
        assign_values<p_nx+p_nu,p_nx> (grad_fxu_p, p, adobj[timestep]);
    }
    void eval_h(index_t timestep, crvec x, crvec u, rvec h) const {
        alpaqa::ScopedMallocAllower ma;
        h.topRows(params.nx)    = x;
        h.bottomRows(params.nu) = u;
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
        auto Jh_xu    = mat::Identity(params.nx + params.nu, params.nx + params.nu);
        auto &&grad_l = h;
        qr            = Jh_xu.transpose() * grad_l;
    }
    void eval_q_N(crvec x, crvec h, rvec q) const {
        alpaqa::ScopedMallocAllower ma;
        auto Jh_x     = mat::Identity(params.nx, params.nx);
        auto &&grad_l = 10 * h;
        q             = Jh_x.transpose() * grad_l;
    }
    void eval_add_Q(index_t timestep, crvec xu, crvec h, rmat Q) const {    
        alpaqa::ScopedMallocAllower ma;
        auto Jh_xu    = mat::Identity(params.nx + params.nu, params.nx + params.nu);
        Q.noalias()   = Jh_xu.transpose() * Jh_xu;
        // Q += mat::Identity(params.nx, params.nx);
    }
    void eval_add_Q_N(crvec x, crvec h, rmat Q) const {
        alpaqa::ScopedMallocAllower ma;
        auto Jh_x     = mat::Identity(params.nx, params.nx);
        Q.noalias()   = 10 * (Jh_x.transpose() * Jh_x);
        // Q += 10 * mat::Identity(params.nx, params.nx);
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
        auto R = mat::Identity(params.nu, params.nu);
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

// ______________________________________________________________________________________________________________________________________________________________________ //
// Nonlinear Dynamics System no. 3 - nuclear batch reactor problem from https://apmonitor.com/do/index.php/Main/DynamicOptimizationBenchmarks

struct NonlinearOCP1c{ 

  USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

  using Box = alpaqa::Box<config_t>;

  // Dynamics and discretization parameters:
  real_t Ts = 0.05;                 ///< Discretization step length

  // OCP parameters:
  length_t N = 20,                    ///< Horizon length
          nu = 1,                     ///< Number of inputs
          nx = 2,                     ///< Number of states  
          nh = nu + nx,               ///< Number of stage outputs
        nh_N = nx,                    ///< Number of terminal outputs
          nc = 0,                     ///< Number of stage constraints
        nc_N = 0,                     ///< Number of terminal constraints
           n = ((nx + nu) * N) - nu;  ///< Total number of decision variables 

  vec Q, R, QR;

  NonlinearOCP1c() : Q(nx), R(nu), QR(nx+nu) {
    Q << 10, 10;
    R << 0.005;
    QR << Q, R;
  }

  [[nodiscard]] length_t get_N() const { return N; }
  [[nodiscard]] length_t get_nu() const { return nu; }
  [[nodiscard]] length_t get_nx() const { return nx; }
  [[nodiscard]] length_t get_nh() const { return nh; }
  [[nodiscard]] length_t get_nh_N() const { return nh_N; }
  [[nodiscard]] length_t get_nc() const { return nc; }
  [[nodiscard]] length_t get_nc_N() const { return nc_N; }

  void get_U(Box &U) const {
    U.lowerbound.setConstant(-alpaqa::inf<config_t>);
    U.upperbound.setConstant(+alpaqa::inf<config_t>);
  }

  void get_D(Box &D) const {    
    D.lowerbound.setConstant(-alpaqa::inf<config_t>);
    D.upperbound.setConstant(+alpaqa::inf<config_t>);
  }

  void get_D_N(Box &D) const {}

  void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const { 
    fxu(0) = x(0) + Ts * (-x(0)*(u(0)+0.5*std::pow(u(0),2)));
    fxu(1) = x(1) + Ts * (u(0)*x(0));
  } // *discretized using Explicit Euler

  void eval_jac_f(index_t timestep, crvec x, crvec u, rmat Jfxu) const {
    Jfxu(0,0) = 1 + Ts*(u(0)+0.5*std::pow(u(0),2));   Jfxu(0,1) = 0;  Jfxu(0,2) = Ts * (x(0)+2*x(0)*u(0));     
    Jfxu(1,0) = Ts * u(0);                            Jfxu(1,1) = 1;  Jfxu(1,2) = Ts * x(0);
  }

  void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p,
                        rvec grad_fxu_p) const {
    grad_fxu_p(0) = (1 + Ts*(u(0)+0.5*std::pow(u(0),2)))*p(0) + (Ts * u(0))*p(1);
    grad_fxu_p(1) = 0*p(0) + 1*p(1);
    grad_fxu_p(2) = (Ts * (x(0)+2*x(0)*u(0)))*p(0) + (Ts * x(0))*p(1);
  }

  void eval_h(index_t timestep, crvec x, crvec u, rvec h) const {
    h.topRows(nx) = x;
    h.bottomRows(nu) = u;
  }

  void eval_h_N(crvec x, rvec h) const { 
    h.topRows(nx) = x;
  }

  [[nodiscard]] real_t eval_l([[maybe_unused]] index_t timestep, crvec h) const {
      alpaqa::ScopedMallocAllower ma;
      return h.transpose() * QR.asDiagonal() * h;
  }
  [[nodiscard]] real_t eval_l_N(crvec h) const {
      alpaqa::ScopedMallocAllower ma;
      return h.transpose() * Q.asDiagonal() * h;
  }
  void eval_qr([[maybe_unused]] index_t timestep, 
                [[maybe_unused]] crvec xu, crvec h, rvec qr) const {
      alpaqa::ScopedMallocAllower ma;
      auto Jh_xu    = mat::Identity(nx + nu, nx + nu);
      auto &&grad_l = (QR.asDiagonal() * h);
      qr            = Jh_xu.transpose() * grad_l;
  }
  void eval_q_N([[maybe_unused]] crvec x, crvec h, rvec q) const {
      alpaqa::ScopedMallocAllower ma;
      auto Jh_x     = mat::Identity(nx, nx);
      auto &&grad_l = 2 * Q.asDiagonal() * h;
      q             = Jh_x.transpose() * grad_l;
  }

  void eval_add_Q(index_t timestep, crvec xu, crvec h, rmat Q) const {
    Q.noalias() += mat::Identity(nx, nx);
  }

  void eval_add_Q_N(crvec x, crvec h, rmat Q) const {
    Q.noalias() += 10 * mat::Identity(nx, nx);
  }

  void eval_add_R_masked(index_t timestep, crvec xu, crvec h, crindexvec mask,
                         rmat R, rvec work) const {
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

  void eval_proj_multipliers(rvec y, real_t M) const {}

  void eval_proj_diff_g(crvec z, rvec p) const { p.setZero(); }

  void check() const {
      // You could do some sanity checks here
  }

  void get_x_init(rvec x_init) const { //need to be implemented! Check ExSession06-MPC 
    x_init.setConstant(0.5);    
  }

};
