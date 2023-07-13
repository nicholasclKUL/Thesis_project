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
#include <climits>
#include <cstdlib>
#include <ctime>
#include <random>
#include <chrono>

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

template <typename X, typename J, typename Params>
void dynamics(index_t &k, X &xu, J &fxu, const Params &params){

  size_t j = 0;

  for (size_t i = 0; i < params.nx; i += 4){
    {
      fxu(i)    = xu(1+i);
      fxu(1+i)  = (-xu(i)+(params.ϵ*std::pow(xu(3+i),2)*sin(xu(2+i))))/(1-(std::pow(params.ϵ,2)*std::pow(cos(xu(2+i)),2)))
                - xu(params.nx+j)*(params.ϵ*cos(xu(2+i)))/(1-(std::pow(params.ϵ,2)*std::pow(cos(xu(2+i)),2)));
      fxu(2+i)  = xu(3+1);
      fxu(3+i)  = params.ϵ*cos(xu(2+i))*(xu(i)-(params.ϵ*std::pow(xu(3+i),2)*sin(xu(2+i))))/(1-(std::pow(params.ϵ,2)*std::pow(cos(xu(2+i)),2)))
                + xu(params.nx+j)*(1-(std::pow(params.ϵ,2)*std::pow(cos(xu(2+i)),2)));
    }
    j++;
  }

};

// Horizon length, number of states and input for compile time allocation of AD views
const int p_N = 16,
p_Nc = 2,
          p_nx = 4*p_Nc,
          p_nu = p_Nc,
          p_nh = p_nx + p_nu;

using AD_obj = AD<p_nx+p_nu,p_nx,p_nh>;

// Create the Control Problem object to pass on to the Solver
struct MultiRTAC{

  using Box = alpaqa::Box<config_t>;

  struct Params{

    real_t   M = 1.7428,
             m = 0.2739,
             e = 0.0537,
             I = 0.000884,
             k = 200,
            mc = 0.0017,
            kc = 0.166,
             ϵ = m*e/std::sqrt((I+m*e*e)*(M+m));

    length_t N_carts  = p_Nc;     ///< Number of balls + fixed and free-end

    length_t   N = p_N,           ///< Horizon length
              nu = p_nu,          ///< Number of inputs
              nx = p_nx,          ///< Number of states : x_0 to x_N and v_0 to v_N,
                                  // *remember f(x_0) = f(v_0) = f(v_N) = 0
              nh = p_nh,          ///< Number of stage outputs
            nh_N = p_nx,            ///< Number of terminal outputs
              nc = 0,             ///< Number of stage constraints
            nc_N = 0,             ///< Number of terminal constraints
              n = ((nx+nu)*N)    ///< Total number of decision variables
                  - nu;

    real_t    T = 1.,
              Ts = T/real_t(N);          ///< Sampling time

  };

  Params params;

unsigned long int n_seed = 1;

  // QR matrix diagonal:
  vec Q, R, QR;

  // Automatic Differentiation object:
  AD_obj ad_obj[p_N];

  MultiRTAC() : Q(params.nx), R(params.nu), QR(params.nx+params.nu) {
    Q.setConstant(1.);
    R.setConstant(.005);
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

    std::default_random_engine re;
    re.seed(n_seed);
    std::uniform_real_distribution<real_t> unif(-0.2, 0.2);

    x_init.setConstant(0.05);
    for (size_t i = 0; i < params.nx; i += 4*params.N_carts){
      for (size_t j = 0; j < params.N_carts; ++j){
        x_init(i+(4*j)) = real_t(j) + 1. + unif(re);
      }
    }

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