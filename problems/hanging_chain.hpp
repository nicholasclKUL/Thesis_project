#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/util/float.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <iomanip>
#include <iostream>
#include <functional>

struct HangingChain{

  USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

  using Box = alpaqa::Box<config_t>;

  // Hanging chain parameters:
  real_t Ts             = 0.01;
  length_t N_balls      = 5;      ///< Number of balls
  real_t β              = 25;
  real_t γ              = 1;
  real_t δ              = 0.01;
  real_t m              = 0.03;   ///< mass
  real_t D              = 0.1;    ///< spring constant
  real_t L              = 0.033;  ///< spring length
  real_t v_max          = 1;      ///< maximum actuator velocity
  real_t g_grav         = 9.81;   ///< Gravitational acceleration       [m/s²]
  real_t x_ref          = 0.3;    ///< desired position of x_N

  length_t N = 50,                     ///< Horizon length
          nu = 1,                     ///< Number of inputs
          nx = (N_balls * 2)-2,       ///< Number of states : x_0 to x_N and v_0 to v_N,
                                        // *remember f(x_0) = f(v_0) = f(v_N) = 0    
          nh = nu + nx,               ///< Number of stage outputs
        nh_N = nx,                    ///< Number of terminal outputs
          nc = 0,                     ///< Number of stage constraints
        nc_N = 0,                     ///< Number of terminal constraints
           n = ((nx + nu) * N) - nu;  ///< Total number of decision variables 

  mat A, B;

  HangingChain() : A(nx, nx), B(nx, nu) {
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
    U.lowerbound.setConstant(-1);
    U.upperbound.setConstant(+1);
  }

  void get_D(Box &D) const {    
    D.lowerbound.setConstant(-alpaqa::inf<config_t>);
    D.upperbound.setConstant(+alpaqa::inf<config_t>);
  }

  void get_D_N(Box &D) const {}

  void get_x_init(rvec x_init) const { 
    real_t Δ = .01;
    real_t k = 0;
    for (index_t i = 0; i < n; ++i){
      x_init(i) = 0.01 + k * Δ;
      ++k;
      if (k == real_t(nh)){
        k = 0;
      }
    } 
    std::cout<<x_init.transpose()<<std::endl;
  }

  void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const { 
    for (index_t i = 0; i < nx; ++i){
      if (i < nx/2){
        if (i == nx/2 - 1){
          fxu(i) = x(i) + Ts*u(0);
        }
        else{
          fxu(i) = x(i) + Ts*x((nx/2)+i);
        }   
      }
      else{
        if (i == nx/2){
          auto F1 = D * (1-(L/std::abs(x(i-nx/2)))) * (x(i-nx/2));
          auto F2 = D * (1-(L/std::abs(x(i-nx/2)-x(i+1-nx/2)))) * (x(i+1-nx/2)-x(i-nx/2));
          fxu(i) = x(i) + Ts * (((1/m) * (F2-F1)) + g_grav);
        }
        else if (i == nx-1){
          auto F1 = D * (1-(L/std::abs(x(i-nx/2)-x(i-1-nx/2)))) * (x(i-nx/2)-x(i-1-nx/2));
          fxu(i) = x(i) + Ts * (((1/m) * (-F1)) + g_grav);
        }
        else{
          auto F1 = D * (1-(L/std::abs(x(i-nx/2)-x(i-1-nx/2)))) * (x(i-nx/2)-x(i-1-nx/2));
          auto F2 = D * (1-(L/std::abs(x(i-nx/2)-x(i+1-nx/2)))) * (x(i+1-nx/2)-x(i-nx/2));
          fxu(i) = x(i) + Ts * (((1/m) * (F2-F1)) + g_grav);
        }
      }
      fxu(0) = 0;
    };
  } // *discretized using Explicit Euler

  void eval_jac_f(index_t timestep, crvec x, crvec u, rmat Jfxu) const {
    Jfxu.setConstant(0.);

    //Jfxu : rows = (2*N_balls) ; col = (2*N_balls + 1)
    
    //filling lower left corner
    //dx_{k-1}
    Jfxu.block(nx/2+1,0,nx/2-1,nx/2-1).setIdentity();     
    //dx_{k},dx_{k+1} 
    for (index_t i = 1; i < nx/2-1; ++i) {
      Jfxu(i+nx/2,i)   = 2 * ((std::abs(x(i)-x(i-1))) - 1);
      Jfxu(i+nx/2,i+1) = (1 - (2*L/std::abs(x(i+1)-x(i))));
    }
    Jfxu *= (D/m);
        
    //filling upper right corner
    Jfxu.block(0,nx/2,nx/2,nx/2).setIdentity();
    Jfxu(nx/2-1,nx) = 1;

    //setting the derivatives wrt to x_0 and v_0 equals to 0  
    Jfxu.row(0).setConstant(0.);
    Jfxu.col(0).setConstant(0.);
    Jfxu.row(nx/2).setConstant(0.);
    Jfxu.col(nx/2).setConstant(0.);
    Jfxu.row(nx-1).setConstant(0.);
    Jfxu.col(nx-1).setConstant(0.);

    //std::cout<<Jfxu.size()<<std::endl;
  }

  void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p,
                        rvec grad_fxu_p) const {
    grad_fxu_p.topRows(nx).noalias() = A.transpose() * p;
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
        Q += mat::Identity(nx, nx);
    }
    void eval_add_Q_N([[maybe_unused]] crvec x,
                      [[maybe_unused]] crvec h, rmat Q) const {
        alpaqa::ScopedMallocAllower ma;
        Q += 10 * mat::Identity(nx, nx);
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

