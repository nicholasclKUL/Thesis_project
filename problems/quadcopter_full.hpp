#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/util/float.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <iomanip>
#include <iostream>
#include <functional>

// __________________________________________________________________________________________________________________________________________________________________________ //
// Quadcopter using AD


#include <Eigen/Core>
#include <Sacado.hpp>
#include <Kokkos_Core.hpp>

#include <thesis/auto-diff-sacado-helpers.hpp>
#include <time-discretization.hpp>

#include <vector>
#include <climits>
#include <cstdlib>
#include <ctime>
#include <random>
#include <chrono>

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

// Dynamics
template <typename X, typename J, typename Params>
void dynamics(index_t &k, X &xu, J &fxu, const Params &params){    
    fxu(0) = xu(1) + xu(3)*(sin(xu(0))*tan(xu(2))) + xu(5)*(cos(xu(0))*tan(xu(2)));
    fxu(1) = - params.a1*xu(3)*xu(5) + params.b1*(xu(12)-xu(14));
    fxu(2) = xu(3)*cos(xu(0)) - xu(5)*sin(xu(0));    
    fxu(3) = - params.a2*xu(1)*xu(5) + params.b2*(xu(13)-xu(15));
    fxu(4) = (cos(xu(0))/cos(xu(2)))*xu(5) + (sin(xu(0))/cos(xu(2)))*xu(3);
    fxu(5) = - params.a3*xu(1)*xu(3) + params.b3*(xu(12)-xu(13)+xu(14)-xu(15));
    fxu(6) = xu(7);
    fxu(7) = params.c1*xu(7) + params.c2*(cos(xu(0))*sin(xu(2))*cos(xu(4)) + sin(xu(0))*sin(xu(4)))*(xu(12)+xu(13)+xu(14)+xu(15));
    fxu(8) = xu(9);
    fxu(9) = params.c1*xu(9) + params.c2*(cos(xu(0))*sin(xu(2))*sin(xu(4)) - sin(xu(0))*cos(xu(4)))*(xu(12)+xu(13)+xu(14)+xu(15));
    fxu(10) = xu(11);
    fxu(11) = - params.g + params.c1*xu(11) + params.c2*(cos(xu(0))*cos(xu(2)))*(xu(12)+xu(13)+xu(14)+xu(15));
};

// Horizon length, number of states and input for compile time allocation of AD views
const int p_N = 60, p_nx = 12, p_nu = 4, p_nh = p_nx + p_nu; 

using AD_obj = AD<p_nx+p_nu,p_nx,p_nh>;

// Create the Control Problem object to pass on to the Solver
struct QuadcopterFull{

  using Box = alpaqa::Box<config_t>;

  struct Params{
    length_t  T = 1.;                     ///< Time horizon (s) 

    // OCP parameters:
    length_t  N = p_N,                   ///< Total horizon length
            nu = p_nu,                     ///< Number of inputs
            nx = p_nx,                    ///< Number of states  
            nh = p_nh,               ///< Number of stage outputs
          nh_N = p_nx,                    ///< Number of terminal outputs
            nc = 0,                     ///< Number of stage constraints
          nc_N = 0,                     ///< Number of terminal constraints
              n = ((nx + nu) * N) - nu;  ///< Total number of decision variables 

    // Dynamics and discretization parameters:
    real_t Ts = real_t(T)/real_t(N),     ///< Discretization step length
          g  = 9.81,                     
          M  = 0.5,                 
          L  = 0.25,
          k  = 3e-6,
          b  = 1e-7,
          kd = 0.25,
          c  = 1e4,
          Jx = 5e-3,
          Jy = 5e-3,
          Jz = 1e-2,

          a1 = (Jy-Jz)/Jx,
          a2 = (Jz-Jx)/Jy,
          a3 = (Jx-Jy)/Jz,
          
          b1 = L*k*c/Jx,
          b2 = L*k*c/Jy,
          b3 = b*c/Jz,
          
          c1 = -(kd/M),
          c2 = k*c/M; 

  };

  Params params;

unsigned long int n_seed = 1;
   
  // QR matrix diagonal:
  vec Q, R, QR; 

  // Automatic Differentiation object:
  AD_obj ad_obj[p_N]; 

  QuadcopterFull() : Q(params.nx), R(params.nu), QR(params.nx+params.nu) {
    Q << 100, 50, 100, 50, 100, 50, 200, 90, 200, 90, 500, 90;
    R << 0.005, 0.005, 0.005, 0.005;
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
    if (x_init.size() <= params.nu*params.N){
      x_init.setConstant(0.0);
      x_init(6) = 0.3;
      x_init(8) = 0.17;
      x_init(10) = 0.2;
    }
    else{
      x_init.setConstant(0.0);
      for (size_t i = 0; i < params.N+1; ++i){
        std::default_random_engine re;
        re.seed(n_seed);

        std::uniform_real_distribution<real_t> unif1(-2, 2);
        x_init(6+(i*(params.nx+params.nu))) = unif1(re);

        std::uniform_real_distribution<real_t> unif2(-2, 2);
        x_init(8+(i*(params.nx+params.nu))) = unif2(re);

        std::uniform_real_distribution<real_t> unif3(-2, 2);
        x_init(10+(i*(params.nx+params.nu))) = unif3(re);   
      }
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
    assign_values_xu<p_nx+p_nu,p_nx>(x, u, ad_obj[timestep]);
    fe(timestep, ad_obj[timestep].xu_fad, ad_obj[timestep].fxu_fad, params);
    // rk4(timestep, ad_obj[timestep].xu_fad, ad_obj[timestep].fxu_fad, params, p_nx+p_nu, p_nx);
    assign_values<p_nx+p_nu,p_nx>(Jfxu, ad_obj[timestep]);
  }

  void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p,
                        rvec grad_fxu_p) const {
    alpaqa::ScopedMallocAllower ma;
    assign_values_xu<p_nx+p_nu,p_nx>(x, u, ad_obj[timestep]);
    fe(timestep, ad_obj[timestep].xu_fad, ad_obj[timestep].fxu_fad, params);
    // rk4(timestep, ad_obj[timestep].xu_fad, ad_obj[timestep].fxu_fad, params, p_nx+p_nu, p_nx);
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
    return h.transpose() * QR.asDiagonal() * h;
  }

  [[nodiscard]] real_t eval_l_N(crvec h) const {      
    alpaqa::ScopedMallocAllower ma;
    return 5 * h.transpose() * Q.asDiagonal() * h;
  }

  void eval_qr([[maybe_unused]] index_t timestep, 
              [[maybe_unused]] crvec xu, crvec h, rvec qr) const {
    alpaqa::ScopedMallocAllower ma;
    auto Jh_xu    = mat::Identity(params.nx + params.nu, params.nx + params.nu);
    auto &&grad_l = QR.asDiagonal() * h;
    qr            = Jh_xu.transpose() * grad_l;
  }

  void eval_q_N([[maybe_unused]] crvec x, crvec h, rvec q) const {
      alpaqa::ScopedMallocAllower ma;
      auto Jh_x     = mat::Identity(params.nx, params.nx);
      auto &&grad_l = 10 * Q.asDiagonal() * h;
      q             = Jh_x.transpose() * grad_l;
  }

  void eval_add_Q([[maybe_unused]] index_t timestep, 
                  [[maybe_unused]] crvec xu, 
                  [[maybe_unused]] crvec h, rmat Q) const {
      alpaqa::ScopedMallocAllower ma;
      Q += mat::Identity(params.nx, params.nx);   
  }

  void eval_add_Q_N([[maybe_unused]] crvec x,
                    [[maybe_unused]] crvec h, rmat Q) const {
      alpaqa::ScopedMallocAllower ma;
      Q += 10 * mat::Identity(params.nx, params.nx);
  }

  void eval_add_R_masked([[maybe_unused]] index_t timestep,
                        [[maybe_unused]] crvec xu, 
                        [[maybe_unused]] crvec h, crindexvec mask,
                        rmat R, 
                        [[maybe_unused]] rvec work) const {
      alpaqa::ScopedMallocAllower ma;
      const auto n = mask.size();
      R.noalias() += mat::Identity(params.n, params.n);
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