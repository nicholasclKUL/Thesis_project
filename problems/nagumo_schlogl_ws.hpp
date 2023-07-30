#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/util/float.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <Eigen/Core>
#include <Sacado.hpp>
#include <Kokkos_Core.hpp>

#include <thesis/auto-diff-sacado-helpers.hpp>

#include <time-discretization.hpp>

#include <iomanip>
#include <iostream>
#include <functional>
#include <vector>
#include <ctime>
#include <random>
#include <chrono>

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

template <typename X, typename J, typename Params>
void dynamics(index_t &k, X &xu, J &fxu, const Params &params){    

auto Δx2 = params.Δx*params.Δx;

for (size_t i = 0; i < fxu.size(); ++i){
  if (i == 0){
    fxu(i) = (2/Δx2)*(-xu(i)+xu(i+1)) - (1/3)*(xu(i)*xu(i)*xu(i)) + xu(i) + xu(i+params.nx);
  }
  else if (i == fxu.size()-1){
    fxu(i) = (2/Δx2)*(xu(i-1)-xu(i)) - (1/3)*(xu(i)*xu(i)*xu(i)) + xu(i) + xu(i+params.nx);
  }
  else{
    fxu(i) = (1/Δx2)*(xu(i-1)-2*xu(i)+xu(i+1)) - (1/3)*(xu(i)*xu(i)*xu(i)) + xu(i) + xu(i+params.nx); 
  }
}

}; 

// Horizon length, number of states and input for compile time allocation of AD views
const int p_N = 16,
          p_nx = 10, 
          p_nu = 10,
          p_nh = p_nx + p_nu; 

// Define the AD object to use in the computations
using AD_obj = AD<p_nx+p_nu,p_nx,p_nh>;

// Create the struct of the Hanging Chain problem
struct Nagumo{

  USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

  using Box = alpaqa::Box<config_t>;

  // Hanging chain parameters:
  struct Params{
    real_t β          = .5;           ///< Q-value horizon-N 
    real_t γ          = .5;           ///< Q-value
    real_t δ          = 0.01;         ///< R-value
    real_t L          = 1.00;         ///< length
    real_t v_max      = 50.;          ///< maximum actuator velocity
    real_t dx1        = 0.35;
    real_t dx2        = 0.65;
  
    length_t   N = p_N,           ///< Horizon length
              nu = p_nu,          ///< Number of inputs
              nx = p_nx,          ///< Number of states : x_0 to x_N and v_0 to v_N, 
                                  // *remember f(x_0) = f(v_0) = f(v_N) = 0    
              nh = p_nh,          ///< Number of stage outputs
            nh_N = p_nx,          ///< Number of terminal outputs
              nc = 0,             ///< Number of stage constraints
            nc_N = 0,             ///< Number of terminal constraints
              n = ((nx+nu)*N)     ///< Total number of decision variables
                  - nu;  

    real_t    T = 0.01,
              Ts = T/real_t(N),   ///< time step size
              Δx = L/real_t(nx);  ///< spatial step size

  };

  Params params;

unsigned long int n_seed = 1;

  AD_obj ad_obj[p_N];       ///< AD objects for each stage (to avoid data race)

  vec Q,                    ///< Diagonal of matrix Q
      R,                    ///< Diagonal of matrix R
      QR;                   ///< Diagonals Q and R merged

  Nagumo() : Q(params.nx), R(params.nu), QR(params.nx+params.nu) {

    Q.setConstant(1.); 
    
    R.setConstant(0.);
    
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
    U.lowerbound.setConstant(-alpaqa::inf<config_t>);
    U.upperbound.setConstant(+alpaqa::inf<config_t>);
  }

  void get_D(Box &D) const {    
    D.lowerbound.setConstant(-alpaqa::inf<config_t>);
    D.upperbound.setConstant(+alpaqa::inf<config_t>);
  }

  void get_D_N(Box &D) const {}

  void get_x_init(rvec x_init) const {
    x_init.setConstant(0);
    // if (x_init.size() > params.nx+params.nu){
    //   for (size_t i = 0; i < params.N-1; ++i){
    //     for (size_t j = 0; j < params.nx; ++j){
    //       if (real_t(j)*params.Δx < params.dx1*params.L){
    //         x_init(i*(params.nx+params.nu)+j) = 0;
    //       }
    //       else if(params.dx2*params.L < real_t(j)*params.Δx){
    //         x_init(i*(params.nx+params.nu)+j) = 0;
    //       }
    //       else {
    //         x_init(i*(params.nx+params.nu)+j) = 2*(1+params.Δx*i);
    //       }
    //     }
    //   }
    // }
    // else {
      for (size_t j = 0; j < params.nx; ++j){
        if (real_t(j)*params.Δx < params.dx1*params.L){
          x_init(j) = 0;
        }
        else if(params.dx2*params.L < real_t(j)*params.Δx){
          x_init(j) = 0;
        }
        else {
          x_init(j) = 5;
        }
      }
    // }
  }

  void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const { 
    alpaqa::ScopedMallocAllower ma;
    vec xu(params.nx+params.nu); xu << x, u;
    // fe(timestep, xu, fxu, params);
    rk4(timestep, xu, fxu, params);
  } 
  void eval_jac_f(index_t timestep, crvec x, crvec u, rmat Jfxu) const {
    alpaqa::ScopedMallocAllower ma;
    assign_values_xu<p_nx+p_nu,p_nx>(x, u, ad_obj[timestep]);
    // fe(timestep, ad_obj[timestep].xu_fad, ad_obj[timestep].fxu_fad, params);
    rk4(timestep, ad_obj[timestep].xu_fad, ad_obj[timestep].fxu_fad, params, p_nx+p_nu, p_nx);
    assign_values<p_nx+p_nu,p_nx>(Jfxu, ad_obj[timestep]);
  }

  void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p,
                        rvec grad_fxu_p) const {
    alpaqa::ScopedMallocAllower ma;
    assign_values_xu<p_nx+p_nu,p_nx>(x, u, ad_obj[timestep]);
    // fe(timestep, ad_obj[timestep].xu_fad, ad_obj[timestep].fxu_fad, params);
    rk4(timestep, ad_obj[timestep].xu_fad, ad_obj[timestep].fxu_fad, params, p_nx+p_nu, p_nx);
    assign_values<p_nx+p_nu,p_nx> (grad_fxu_p, p, ad_obj[timestep]);
  }

  void eval_h([[maybe_unused]] index_t timestep, crvec x, crvec u, rvec h) const {
      alpaqa::ScopedMallocAllower ma;
      h.topRows(params.nx)    = x;
      h.bottomRows(params.nu) = u;
      // set target (set) points    
      for (size_t j = 0; j < params.nx; ++j){
        if (real_t(j)*params.Δx < params.dx1*params.L){
          h(j) = h(j);
        }
        else if(params.dx2*params.L < real_t(j)*params.Δx){
          h(j) = h(j);
        }
        else {
          h(j) = h(j) - 2;
        }
      }
  }
  void eval_h_N(crvec x, rvec h) const {       
      alpaqa::ScopedMallocAllower ma;
      h = x;
      // set target (set) points    
      for (size_t j = 0; j < params.nx; ++j){
        if (real_t(j)*params.Δx < params.dx1*params.L){
          h(j) = h(j);
        }
        else if(params.dx2*params.L < real_t(j)*params.Δx){
          h(j) = h(j);
        }
        else {
          h(j) = h(j) - 2;
        }
      }
  }
  [[nodiscard]] real_t eval_l([[maybe_unused]] index_t timestep, crvec h) const {
      alpaqa::ScopedMallocAllower ma;
      return params.γ * h.transpose() * QR.asDiagonal() * h;
  }
  [[nodiscard]] real_t eval_l_N(crvec h) const {
      alpaqa::ScopedMallocAllower ma;
      return params.β * h.transpose() * Q.asDiagonal() * h;
  }
  void eval_qr([[maybe_unused]] index_t timestep, 
                [[maybe_unused]] crvec xu, crvec h, rvec qr) const {
      alpaqa::ScopedMallocAllower ma;
      auto Jh_xu    = mat::Identity(params.nx + params.nu, params.nx + params.nu);
      auto &&grad_l = 2 * params.γ * (QR.asDiagonal() * h);
      qr            = Jh_xu.transpose() * grad_l;
  }
  void eval_q_N([[maybe_unused]] crvec x, crvec h, rvec q) const {
      alpaqa::ScopedMallocAllower ma;
      auto Jh_x     = mat::Identity(params.nx, params.nx);
      auto &&grad_l = 2 * params.β * Q.asDiagonal() * h;
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

