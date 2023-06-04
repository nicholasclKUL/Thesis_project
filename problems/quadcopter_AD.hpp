#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/util/float.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <Eigen/Core>
#include <Sacado.hpp>
#include <Kokkos_Core.hpp>

#include <thesis/auto-diff-sacado-helpers.hpp>

#include <iomanip>
#include <iostream>
#include <functional>
#include <vector>

//dynamic model described in https://www.codeproject.com/Articles/5335232/Optimal-Control-of-a-Quadcopter#l2

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

// Dynamics
template <typename X, typename J>
void dynamics(index_t &k, X &xu, J &fxu,
             real_t Ts, real_t a1, real_t a2, real_t a3,
             real_t b1, real_t b2, real_t M, real_t g){
    fxu(0) = xu(0) + Ts * xu(1);
    fxu(1) = xu(1) + Ts * (a1*xu(3)*xu(5) + b1*xu(13));
    fxu(2) = xu(2) + Ts * xu(3);    
    fxu(3) = xu(3) + Ts * (a2*xu(1)*xu(5) + b2*xu(14));
    fxu(4) = xu(4) + Ts * xu(5);
    fxu(5) = xu(5) + Ts * (a3*xu(1)*xu(3) + b1*xu(15));
    fxu(6) = xu(6) + Ts * xu(7);
    fxu(7) = xu(7) + Ts * (cos(xu(0))*sin(xu(2))*cos(xu(4)) + sin(xu(0))*sin(xu(4))) * xu(15) / M;
    fxu(8) = xu(8) + Ts * xu(9);
    fxu(9) = xu(9) + Ts * (cos(xu(0))*sin(xu(2))*sin(xu(4)) - sin(xu(0))*cos(xu(4))) * xu(15) / M;
    fxu(10) = xu(10) + Ts * xu(11);
    fxu(11) = xu(11) + Ts * (cos(xu(0))*cos(xu(2))) * xu(15) / (M-g);
}

// Horizon length, number of states and input for compile time allocation of AD views
const int p_N = 30, p_nx = 12, p_nu = 4, p_nh = 16; 

using AD_obj = AD<p_nx+p_nu,p_nx,p_nh>;

// Create the Control Problem object to pass on to the Solver
struct QuadcopterAD{

  using Box = alpaqa::Box<config_t>;

  length_t  T = 3;                      ///< Time horizon (s) 

  // OCP parameters:
  length_t  N = p_N,                   ///< Total horizon length
           nu = 4,                     ///< Number of inputs
           nx = 12,                    ///< Number of states  
           nh = nu + nx,               ///< Number of stage outputs
         nh_N = nx,                    ///< Number of terminal outputs
           nc = 0,                     ///< Number of stage constraints
         nc_N = 0,                     ///< Number of terminal constraints
            n = ((nx + nu) * N) - nu;  ///< Total number of decision variables 

  // Dynamics and discretization parameters:
  real_t Ts = real_t(T)/real_t(N),          ///< Discretization step length
         g  = 9.81,                     
         M  = 0.65,                 
         I  = 0.23,
         Jx = 7.5e-3,
         Jy = 7.5e-3,
         Jz = 1.3e-2,

         a1 = (Jy-Jz)/Jx,
         a2 = (Jz-Jx)/Jy,
         a3 = (Jx-Jy)/Jz,
         b1 = I/Jx,
         b2 = I/Jy,
         b3 = I/Jz;
  
  // QR matrix diagonal:
  vec Q, R, QR; 

  // Automatic Differentiation object:
  AD_obj ad_obj[p_N]; 

  QuadcopterAD() : Q(nx), R(nu), QR(nx+nu) {
    Q << 200, 200, 500, 50, 50, 50, 100, 100, 100, 90, 90, 90;
    R << 0.005, 0.005, 0.005, 0.005;
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

  void get_x_init(rvec x_init) const {
    if (x_init.size() <= nu*N){
      x_init.setConstant(0);
      x_init(0) = -0.25;
      x_init(1) = 0.51;
      x_init(2) = 0.32;
    }
    else{
      x_init.setConstant(0.);
      for (size_t i = 0; i < N-1; ++i){
        x_init(0+(i*(nx+nu))) = -0.25; //x
        x_init(1+(i*(nx+nu))) = 0.51;  //y
        x_init(2+(i*(nx+nu))) = 0.32;  //z
      }
    }
  }

  void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const { 
    alpaqa::ScopedMallocAllower ma;

    fxu(0) = x(0) + Ts * x(1);
    fxu(1) = x(1) + Ts * (a1*x(3)*x(5) + b1*u(1));
    fxu(2) = x(2) + Ts * x(3);    
    fxu(3) = x(3) + Ts * (a2*x(1)*x(5) + b2*u(2));
    fxu(4) = x(4) + Ts * x(5);
    fxu(5) = x(5) + Ts * (a3*x(1)*x(3) + b1*u(3));
    fxu(6) = x(6) + Ts * x(7);
    fxu(7) = x(7) + Ts * (cos(x(0))*sin(x(2))*cos(x(4)) + sin(x(0))*sin(x(4))) * u(0) / M;
    fxu(8) = x(8) + Ts * x(9);
    fxu(9) = x(9) + Ts * (cos(x(0))*sin(x(2))*sin(x(4)) - sin(x(0))*cos(x(4))) * u(0) / M;
    fxu(10) = x(10) + Ts * x(11);
    fxu(11) = x(11) + Ts * (cos(x(0))*cos(x(2))) * u(0) / (M-g);

  } // *discretized using Explicit Euler

  void eval_jac_f(index_t timestep, crvec x, crvec u, rmat Jfxu) const {
    alpaqa::ScopedMallocAllower ma;

    // Assign values of x and u to their respective views
    assign_values_xu<p_nx+p_nu,p_nx,p_nh>(x, u, ad_obj[timestep]);

    // Calculate Jacobian of system dynamics using AD
    dynamics(timestep, ad_obj[timestep].xu_fad, ad_obj[timestep].fxu_fad, 
              Ts, a1, a2, a3, b1, b2, M, g);
   
    // Assigning Jfxu AD view to input container
    assign_values<p_nx+p_nu,p_nx,p_nh>(Jfxu, ad_obj[timestep]);

  }

  void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p,
                        rvec grad_fxu_p) const {
    alpaqa::ScopedMallocAllower ma;

    assign_values_xu<p_nx+p_nu,p_nx,p_nh>(x, u, ad_obj[timestep]);

    dynamics(timestep, ad_obj[timestep].xu_fad, ad_obj[timestep].fxu_fad, 
              Ts, a1, a2, a3, b1, b2, M, g);

    assign_values<p_nx+p_nu,p_nx,p_nh>(grad_fxu_p, p, ad_obj[timestep]);

  }

  void eval_h([[maybe_unused]] index_t timestep, crvec x, crvec u, rvec h) const {        
    alpaqa::ScopedMallocAllower ma;
    h.topRows(nx)    = x;
    h.bottomRows(nu) = u;
  }

  void eval_h_N(crvec x, rvec h) const { h = x; }

  [[nodiscard]] real_t eval_l([[maybe_unused]] index_t timestep, crvec h) const {      
    alpaqa::ScopedMallocAllower ma;
    return h.transpose() * QR.asDiagonal() * h;
  }

  [[nodiscard]] real_t eval_l_N(crvec h) const {      
    alpaqa::ScopedMallocAllower ma;
    return 5 * h.transpose()* Q.asDiagonal() * h;
  }

  void eval_qr([[maybe_unused]] index_t timestep, 
              [[maybe_unused]] crvec xu, crvec h, rvec qr) const {
    alpaqa::ScopedMallocAllower ma;
    auto Jh_xu    = mat::Identity(nx + nu, nx + nu);
    auto &&grad_l = QR.asDiagonal() * h;
    qr            = Jh_xu.transpose() * grad_l;
  }

  void eval_q_N([[maybe_unused]] crvec x, crvec h, rvec q) const {
      alpaqa::ScopedMallocAllower ma;
      auto Jh_x     = mat::Identity(nx, nx);
      auto &&grad_l = 10 * Q.asDiagonal() * h;
      q             = Jh_x.transpose() * grad_l;
  }

  void eval_add_Q([[maybe_unused]] index_t timestep, 
                  [[maybe_unused]] crvec xu, 
                  [[maybe_unused]] crvec h, rmat Q) const {
      alpaqa::ScopedMallocAllower ma;
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