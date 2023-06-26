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

    // ẋₖ = vₖ
    for (size_t i = 0; i < (params.N_balls-1)*params.dim-params.dim; ++i){
      fxu(i) = xu(i+(params.N_balls-1)*params.dim);
    } 
    // ẋ₀ = 0 OR x₀ = p₀
    fxu((params.N_balls-1)*params.dim-3) = xu(xu.size()-3); 
    fxu((params.N_balls-1)*params.dim-2) = xu(xu.size()-2);
    fxu((params.N_balls-1)*params.dim-1) = xu(xu.size()-1);

    // v̇ₖ = (1/m) * (Fₖ₊₁ - Fₖ₋₁) + g
    for (size_t i = 0; i < (params.N_balls-2)*params.dim; i += size_t(params.dim)){
      if (i == 0){

        fxu(i+(params.N_balls-1)*params.dim) = ((params.D/params.m)*((1-(params.L/(std::sqrt(std::pow(xu(i+params.dim)-xu(i),2)+
            std::pow(xu(1+i+params.dim)-xu(1+i),2)+std::pow(xu(2+i+params.dim)-xu(2+i),2)))))*(xu(i+params.dim)-xu(i))
            -(1-(params.L/(std::sqrt(std::pow(xu(i)-params.xo,2)+std::pow(xu(1+i)-params.yo,2)+std::pow(xu(2+i)-params.zo,2)))))*(xu(i)-params.xo)));
        
        fxu(i+1+(params.N_balls-1)*params.dim) = ((params.D/params.m)*((1-(params.L/(std::sqrt(std::pow(xu(i+params.dim)-xu(i),2)+
            std::pow(xu(1+i+params.dim)-xu(1+i),2)+std::pow(xu(2+i+params.dim)-xu(2+i),2)))))*(xu(1+i+params.dim)-xu(1+i))
            -(1-(params.L/(std::sqrt(std::pow(xu(i)-params.xo,2)+std::pow(xu(1+i)-params.yo,2)+std::pow(xu(2+i)-params.zo,2)))))*(xu(1+i)-params.yo)));
        
        fxu(i+2+(params.N_balls-1)*params.dim) = ((params.D/params.m)*((1-(params.L/(std::sqrt(std::pow(xu(i+params.dim)-xu(i),2)+
            std::pow(xu(1+i+params.dim)-xu(1+i),2)+std::pow(xu(2+i+params.dim)-xu(2+i),2)))))*(xu(2+i+params.dim)-xu(2+i))
            -(1-(params.L/(std::sqrt(std::pow(xu(i)-params.xo,2)+std::pow(xu(1+i)-params.yo,2)+std::pow(xu(2+i)-params.zo,2)))))*(xu(2+i)-params.zo))+params.g);
      
      }
      else {

        fxu(i+(params.N_balls-1)*params.dim) = ((params.D/params.m)*((1-(params.L/(std::sqrt(std::pow(xu(i+params.dim)-xu(i),2)+
            std::pow(xu(1+i+params.dim)-xu(1+i),2)+std::pow(xu(2+i+params.dim)-xu(2+i),2)))))*(xu(i+params.dim)-xu(i))
            -(1-(params.L/(std::sqrt(std::pow(xu(i)-xu(i-params.dim),2)+std::pow(xu(1+i)-xu(1+i-params.dim),2)+std::pow(xu(2+i)-xu(2+i-params.dim),2)))))*(xu(i)-xu(i-params.dim))));
        
        fxu(i+1+(params.N_balls-1)*params.dim) = ((params.D/params.m)*((1-(params.L/(std::sqrt(std::pow(xu(i+params.dim)-xu(i),2)+
            std::pow(xu(1+i+params.dim)-xu(1+i),2)+std::pow(xu(2+i+params.dim)-xu(2+i),2)))))*(xu(1+i+params.dim)-xu(1+i))
            -(1-(params.L/(std::sqrt(std::pow(xu(i)-xu(i-params.dim),2)+std::pow(xu(1+i)-xu(1+i-params.dim),2)+std::pow(xu(2+i)-xu(2+i-params.dim),2)))))*(xu(1+i)-xu(1+i-params.dim))));
        
        fxu(i+2+(params.N_balls-1)*params.dim) = ((params.D/params.m)*((1-(params.L/(std::sqrt(std::pow(xu(i+params.dim)-xu(i),2)+
            std::pow(xu(1+i+params.dim)-xu(1+i),2)+std::pow(xu(2+i+params.dim)-xu(2+i),2)))))*(xu(2+i+params.dim)-xu(2+i))
            -(1-(params.L/(std::sqrt(std::pow(xu(i)-xu(i-params.dim),2)+std::pow(xu(1+i)-xu(1+i-params.dim),2)+std::pow(xu(2+i)-xu(2+i-params.dim),2)))))*(xu(2+i)-xu(2+i-params.dim)))+params.g);
      
      }
    } 
}; 

// Horizon length, number of states and input for compile time allocation of AD views
const int p_N = 30,
          p_Nb = 3,
          p_dim = 3,
          p_nx = ((p_Nb - 1) * p_dim * 2) - p_dim, //< ((N_balls - 1) * dim * 2) - dim
          p_nu = 3,
          p_nh = p_nx + p_nu; 

// Define the AD object to use in the computations
using AD_obj = AD<p_nx+p_nu,p_nx,p_nh>;

// Create the struct of the Hanging Chain problem
struct HangingChain{

  USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

  using Box = alpaqa::Box<config_t>;

  // Hanging chain parameters:
  struct Params{
    real_t β          = 25;       ///< Q-value for free-end 
    real_t γ          = 1;        ///< Q-value
    real_t δ          = 0.01;     ///< R-value
    real_t m          = 0.03;     ///< mass
    real_t D          = 1.0;      ///< spring constant
    real_t L          = 0.033;    ///< spring length
    real_t v_max      = 0.5;      ///< maximum actuator velocity
    real_t g          = -9.81;    ///< Gravitational acceleration, [m/s²]
    length_t N_balls  = p_Nb;     ///< Number of balls + fixed and free-end
    length_t dim      = p_dim;    ///< Dimensions
    real_t xo         = 0.0;      ///< Fixed-end position in x
    real_t yo         = 0.0;      ///< Fixed-end position in y
    real_t zo         = 0.0;      ///< Fixed-end position in z
    real_t a          = -m*g;
    real_t b          = (N_balls+1)*L/2;
  
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

  AD_obj ad_obj[p_N];           ///< AD objects for each stage (to avoid data race)

  vec p_ref,                    ///< desired position of x_N, y_N and z_N
         po,                    ///< fixed-end position
          Q,                    ///< Diagonal of matrix Q
          R,                    ///< Diagonal of matrix R
         QR;                    ///< Diagonals Q and R merged

  HangingChain() : p_ref(params.dim), po(params.dim), Q(params.nx), R(params.nu), QR(params.nx+params.nu) {

    po << params.xo, params.yo, params.zo;
    
    Q.segment(0,(params.N_balls-1)*params.dim).setConstant(0); 
    Q((params.N_balls-1)*params.dim-1) = params.β; 
    Q((params.N_balls-1)*params.dim-2) = params.β; 
    Q((params.N_balls-1)*params.dim-3) = params.β;
    Q.segment((params.N_balls-1)*params.dim,params.nx-(params.N_balls-1)*params.dim).setConstant(params.γ);
    
    R.setConstant(params.δ);
    
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
    U.lowerbound.setConstant(-1);
    U.upperbound.setConstant(+1);
  }

  void get_D(Box &D) const {    
    D.lowerbound.setConstant(-alpaqa::inf<config_t>);
    D.upperbound.setConstant(+alpaqa::inf<config_t>);
  }

  void get_D_N(Box &D) const {}

  void get_x_init(rvec x_init) const {
    if (x_init.size() <= 0){
      x_init.setConstant(0);
    }
    else {
      real_t k = 0;
      for (index_t i = 0; i < x_init.size(); i += 3){
        if (k*3 >= params.nx/2){
          x_init(i) = 0.0;
          x_init(i+1) = 0.0;
          x_init(i+2) = 0.0;
        } 
        else{
          x_init(i) = (k+1)*params.L;
          x_init(i+1) = 0;
          x_init(i+2) = params.a*std::cosh((x_init(i)-params.b)/params.a) - params.a*std::cosh(-params.b/params.a);
        }
        ++k;
        if (k*3 == real_t(p_nh)){
          k = 0;
        }
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
      // set target (set) points    
      std::default_random_engine re;
      re.seed(n_seed);
      std::uniform_real_distribution<real_t> unif(-0.2, 0.2);
      h((params.N_balls-1)*params.dim-3) = x((params.N_balls-1)*params.dim-3) - unif(re);
      h((params.N_balls-1)*params.dim-2) = x((params.N_balls-1)*params.dim-2) - unif(re);
      h((params.N_balls-1)*params.dim-1) = x((params.N_balls-1)*params.dim-1) - unif(re);
  }
  void eval_h_N(crvec x, rvec h) const {       
      h = x;
      // set target (set) points   
      std::default_random_engine re; 
      re.seed(n_seed);
      std::uniform_real_distribution<real_t> unif(-0.2, 0.2);
      h((params.N_balls-1)*params.dim-3) = x((params.N_balls-1)*params.dim-3) - unif(re);
      h((params.N_balls-1)*params.dim-2) = x((params.N_balls-1)*params.dim-2) - unif(re);
      h((params.N_balls-1)*params.dim-1) = x((params.N_balls-1)*params.dim-1) - unif(re); 
  }
  [[nodiscard]] real_t eval_l([[maybe_unused]] index_t timestep, crvec h) const {
      alpaqa::ScopedMallocAllower ma;
      return h.transpose() * QR.asDiagonal() * h;
  }
  [[nodiscard]] real_t eval_l_N(crvec h) const {
      alpaqa::ScopedMallocAllower ma;
      return params.β * h.transpose() * Q.asDiagonal() * h;
  }
  void eval_qr([[maybe_unused]] index_t timestep, 
                [[maybe_unused]] crvec xu, crvec h, rvec qr) const {
      alpaqa::ScopedMallocAllower ma;
      auto Jh_xu    = mat::Identity(params.nx + params.nu, params.nx + params.nu);
      auto &&grad_l = (QR.asDiagonal() * h) * params.β;
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

