#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/print.hpp>
#include <iomanip>
#include <iostream>

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

struct NonlinearOCP1b {

    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    using Box = alpaqa::Box<config_t>;

    length_t N = 4,      ///< Horizon length
        nu     = 1,       ///< Number of inputs
        nx     = 100,       ///< Number of states
        nh     = nu + nx, ///< Number of stage outputs
        nh_N   = nx,      ///< Number of terminal outputs
        nc     = 0,       ///< Number of stage constraints
        nc_N   = 0;       ///< Number of terminal constraints

    mat A, B;

    NonlinearOCP1b() : A(nx, nx), B(nx, nu) {
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
        for (size_t i = 0; i < nx; ++i) {
            if (i == 0){
                fxu(i) = u(0);
            }
            else {
                fxu(i) = std::pow(x(i-1),2);
            }
        }
        
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

struct NonlinearOCP1c{ // UNFINISHED!

  USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

  using Box = alpaqa::Box<config_t>;

  // Dynamics and discretization parameters:
  real_t Ts = 0.05;                 ///< Discretization step length

  // OCP parameters:
  length_t N = 20,                    ///< Horizon length
          nu = 1,                     ///< Number of inputs
          nx = 3,                     ///< Number of states  
          nh = nu + nx,               ///< Number of stage outputs
        nh_N = nx,                    ///< Number of terminal outputs
          nc = 0,                     ///< Number of stage constraints
        nc_N = 0,                     ///< Number of terminal constraints
           n = ((nx + nu) * N) - nu;  ///< Total number of decision variables 

  mat A, B;

  NonlinearOCP1c() : A(nx, nx), B(nx, nu) {
    A.setZero();
    B.setZero();
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
    fxu(0) = x(0) + Ts * (- 0.877*x(0) + x(2) - 0.088*x(0)*x(2) 
                          + 0.47*std::pow(x(0),2) - 0.019*std::pow(x(1),2) - x(2)*std::pow(x(0),2)
                          + 3.846*std::pow(x(0),3) 
                          - 0.215*u(0) + 0.28*std::pow(x(0),2)*u(0) + 0.47*std::pow(u(0),2)*x(0) + 0.63*std::pow(u(0),3) );
    fxu(1) = x(1) + Ts * (x(2));
    fxu(2) = x(2) + Ts * (- 4.208*x(0) - 0.396*x(2) - 0.47*std::pow(x(0),2) - 3.564*std::pow(x(0),3)
                          - 20.967*u(0) + 6.265*std::pow(x(0),2)*u(0) + 46*std::pow(u(0),2)*x(0) + 61.4*std::pow(u(0),3));
  } // *discretized using Explicit Euler

  void eval_jac_f(index_t timestep, crvec x, crvec u, rmat Jfxu) const {
   //[d/dx1]                     [d/dx2]                      [d/dx3]                     [d/du]
    Jfxu(0,0) = 1 + (Ts * 10);   Jfxu(0,1) = - (Ts * 10);     Jfxu(0,2) = 0;              Jfxu(0,3) = 0;
    Jfxu(1,0) = 0;               Jfxu(1,1) = 1 - Ts*2*x(1);   Jfxu(1,2) = 0;              Jfxu(1,3) = Ts;
    Jfxu(2,0) = Ts * x(1);       Jfxu(2,1) = Ts * x(1);       Jfxu(2,2) = 1 - (Ts * 3);   Jfxu(2,3) = 0;
  }

  void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p,
                        rvec grad_fxu_p) const {
    grad_fxu_p.topRows(nx).noalias() = A.transpose() * p;
    grad_fxu_p.bottomRows(nu).noalias() = B.transpose() * p;
  }

  void eval_h(index_t timestep, crvec x, crvec u, rvec h) const {
    h.topRows(nx) = x;
    h.bottomRows(nu) = u;
  }

  void eval_h_N(crvec x, rvec h) const { 
    h.topRows(nx) = x;
  }

  [[nodiscard]] real_t eval_l(index_t timestep, crvec h) const {
    return 10*h.topRows(nx).squaredNorm() + 0.005*h.bottomRows(nu).squaredNorm();
  }

  [[nodiscard]] real_t eval_l_N(crvec h) const {
    return 50*h.squaredNorm();
  }

  void eval_qr(index_t timestep, crvec xu, crvec h, rvec qr) const {
    auto Jh_xu = mat::Identity(nx + nu, nx + nu);
    vec grad_l;
    grad_l = vec::Zero(xu.size());
    grad_l.topRows(nx) = 20 * h.segment(0,nx);
    grad_l.bottomRows(nu) = 0.01 * h.segment(nx,nu); 
    qr = Jh_xu.transpose() * grad_l; 
  }

  void eval_q_N(crvec x, crvec h, rvec q) const {
    auto Jh_x = mat::Identity(nx, nx);
    auto &&grad_l = 100 * h;
    q = Jh_x.transpose() * grad_l;
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

  [[nodiscard]] length_t get_S_work_size() const { return 0; }

  void check() const {

  // You could do some sanity checks here

  }

  void get_x_init(rvec x_init) const { //need to be implemented! Check ExSession06-MPC 
    x_init(0) = 0.4655;
    x_init(1) = 0.;
    x_init(2) = 0.;
    x_init(3) = 0.;
    for (index_t i = 1; i < N; ++i)
    {
      if (i == N-1){
        eval_f(i, x_init.segment((i-1)*(nx+nu),nx), 
            x_init.segment((i-1)*(nx+nu)+nx,nu),
            x_init.segment(i*(nx+nu),nx));  
      } else {
        eval_f(i, x_init.segment((i-1)*(nx+nu),nx), 
            x_init.segment((i-1)*(nx+nu)+nx,nu),
            x_init.segment(i*(nx+nu),nx));
        x_init.segment(i*(nx+nu)+nx,nu).setConstant(0.);
      }
    }    
  }

};
