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
             real_t b1, real_t b2, real_t b3,
             real_t M, real_t g){

    fxu(0) = xu(0) + Ts * xu(1);
    fxu(1) = xu(1) + Ts * (a1*xu(3)*xu(5) + b1*xu(13));
    fxu(2) = xu(2) + Ts * xu(3);    
    fxu(3) = xu(3) + Ts * (a2*xu(1)*xu(5) + b2*xu(14));
    fxu(4) = xu(4) + Ts * xu(5);
    fxu(5) = xu(5) + Ts * (a3*xu(1)*xu(3) + b3*xu(15));
    fxu(6) = xu(6) + Ts * xu(7);
    fxu(7) = xu(7) + Ts * (cos(xu(0))*sin(xu(2))*cos(xu(4)) + sin(xu(0))*sin(xu(4))) * xu(12) / M;
    fxu(8) = xu(8) + Ts * xu(9);
    fxu(9) = xu(9) + Ts * (cos(xu(0))*sin(xu(2))*sin(xu(4)) - sin(xu(0))*cos(xu(4))) * xu(12) / M;
    fxu(10) = xu(10) + Ts * xu(11);
    fxu(11) = xu(11) + Ts * (cos(xu(0))*cos(xu(2))) * xu(12) / (M-g);

}

// Output Mapping
template <typename T>
void h(index_t &k, T &xu, T &h){
    h = xu;
}

// Cost function
template <typename T>
real_t l(index_t &k, T &h){
  0.5 * h.squaredNorm();
}

// Terminal Cost function
template <typename T>
real_t l(T &h){
  5 * h.squaredNorm();
}

// Horizon length, number of states and input for compile time allocation
const int p_N = 60, p_nx = 12, p_nu = 4;


// Create the Control Problem object to pass on to the Solver
struct QuadcopterAD{

  using Box = alpaqa::Box<config_t>;

  // OCP parameters:
  length_t N = p_N,                   ///< Horizon length
          nu = p_nu,                  ///< Number of inputs
          nx = p_nx,                  ///< Number of states  
          nh = nu + nx,               ///< Number of stage outputs
        nh_N = nx,                    ///< Number of terminal outputs
          nc = 0,                     ///< Number of stage constraints
        nc_N = 0,                     ///< Number of terminal constraints
           n = ((nx + nu) * N) - nu;  ///< Total number of decision variables 

  // Dynamics and discretization parameters:
  real_t Ts = 1/real_t(N),          ///< Discretization step length
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

  // Automatic Differentiation object:
  AD<p_nu+p_nx,p_nx> jac_ad[p_N]; 

  mat A, B;

  QuadcopterAD() : A(nx, nx), B(nx, nu) {}

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
    x_init.setConstant(1.);
  }

  void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const { 
    
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

    // Assign values of x and u to their respective views
    assign_values_xu<16,12>(x, u, jac_ad[timestep]);

    // Calculate Jacobian of system dynamics using AD
    dynamics(timestep, jac_ad[timestep].xu_fad, jac_ad[timestep].fxu_fad, Ts, a1, a2, a3, b1, b2, b3, M, g);
   
    // Assigning Jfxu AD view to input container
    for (size_t i = 0; i < nx; ++i){
      for (size_t j = 0; j < nx+nu; j++){
        Jfxu(i,j) = jac_ad[timestep].fxu_fad(i).fastAccessDx(j);
      }
    }

  }

  void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p,
                        rvec grad_fxu_p) const {

    // Assign values of x and u to their respective views
    assign_values_xu<16,12>(x, u, jac_ad[timestep]);

    // Calculate Jacobian of system dynamics using AD
    dynamics(timestep, jac_ad[timestep].xu_fad, jac_ad[timestep].fxu_fad, Ts, a1, a2, a3, b1, b2, b3, M, g);
    
    // Assigning grad-vector product to input container
    for (size_t i = 0; i < nx+nu; ++i){
      for (size_t j = 0; j < nx; j++){
        grad_fxu_p(i) += jac_ad[timestep].fxu_fad(j).fastAccessDx(i)*p(j);
      }
    }

  }

  void eval_h([[maybe_unused]] index_t timestep, crvec x, crvec u, rvec h) const {        
    h.topRows(nx)    = x;
    h.bottomRows(nu) = u;
  }

  void eval_h_N(crvec x, rvec h) const { h = x; }

  [[nodiscard]] real_t eval_l([[maybe_unused]] index_t timestep, crvec h) const {      
    return 0.5 * h.squaredNorm();
  }

  [[nodiscard]] real_t eval_l_N(crvec h) const {      
    return 5. * h.squaredNorm();
  }

  void eval_qr([[maybe_unused]] index_t timestep, 
              [[maybe_unused]] crvec xu, crvec h, rvec qr) const {
      
      auto Jh_xu    = mat::Identity(nx + nu, nx + nu);
      auto &&grad_l = h;
      qr            = Jh_xu.transpose() * grad_l;
  }

  void eval_q_N([[maybe_unused]] crvec x, crvec h, rvec q) const {
      
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
      
      Q += 10 * mat::Identity(nx, nx);
  }

  void eval_add_R_masked([[maybe_unused]] index_t timestep,
                        [[maybe_unused]] crvec xu, 
                        [[maybe_unused]] crvec h, crindexvec mask,
                        rmat R, 
                        [[maybe_unused]] rvec work) const {
      
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