#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/util/float.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <iomanip>
#include <iostream>
#include <functional>

//dynamic model described in https://www.codeproject.com/Articles/5335232/Optimal-Control-of-a-Quadcopter#l2

struct Quadcopter{

  USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

  using Box = alpaqa::Box<config_t>;

  length_t  T = 1;                      ///< Time horizon (s) 

  // OCP parameters:
  length_t N = 10,                     ///< Total horizon length
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

  Quadcopter() : Q(nx), R(nu), QR(nx+nu) {
    Q << 100, 50, 100, 50, 100, 50, 200, 90, 200, 90, 500, 90;
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
      x_init.setConstant(0.0);
      x_init(6) = 0.5;
      x_init(8) = -0.33;
      x_init(10) = 0.2;
    }
    else{
      x_init.setConstant(0.0);
      for (size_t i = 0; i < N+1; ++i){
        x_init(6+(i*(nx+nu))) = 0.5;
        x_init(8+(i*(nx+nu))) = -0.33;
        x_init(10+(i*(nx+nu))) = 0.2;  
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

    // "self"-derivatives
    Jfxu.setIdentity();
    
    // dx1/dx
    Jfxu(0,1) = Ts;
    // dx2/dx
    Jfxu(1,3) = Ts * (a1*x(5));   
    Jfxu(1,5) = Ts * (a1*x(3));
    // dx3/dx
    Jfxu(2,3) = Ts;
    // dx4/dx
    Jfxu(3,1) = Ts * (a2*x(5));    
    Jfxu(3,5) = Ts * (a2*x(1));
    // dx5/dx
    Jfxu(4,5) = Ts;
    // dx6/dx
    Jfxu(5,1) = Ts * (a3*x(3));    
    Jfxu(5,3) = Ts * (a3*x(1));
    // dx7/dx
    Jfxu(6,7) = Ts;
    // dx8/dx
    Jfxu(7,0) = Ts * (-sin(x(0))*sin(x(2))*cos(x(4)) + cos(x(0))*sin(x(4))) * u(0) / M;    
    Jfxu(7,2) = Ts * (cos(x(0))*cos(x(2))*cos(x(4))) * u(0) / M;    
    Jfxu(7,4) = Ts * (-cos(x(0))*sin(x(4))*sin(x(2)) + sin(x(0))*cos(x(4))) * u(0) / M;
    // dx9/dx
    Jfxu(8,9) = Ts;
    // dx10/dx
    Jfxu(9,0) = Ts * (-sin(x(0))*sin(x(2))*sin(x(4)) - cos(x(0))*cos(x(4))) * u(0) / M;    
    Jfxu(9,2) = Ts * (cos(x(2))*cos(x(0))*sin(x(4))) * u(0) / M;    
    Jfxu(9,4) = Ts * (cos(x(4))*cos(x(0))*sin(x(2)) + sin(x(4))*sin(x(0))) * u(0) / M;
    // dx11/dx
    Jfxu(10,11) = Ts;
    // dx12/dx
    Jfxu(11,0) = Ts * (-sin(x(0))*cos(x(2))) * u(0) / (M-g);    
    Jfxu(11,2) = Ts * (-sin(x(2))*cos(x(0))) * u(0) / (M-g);   
    
    // dx2/du
    Jfxu(1,13) = Ts * b1;
    // dx4/du
    Jfxu(3,14) = Ts * b2;
    // dx6/du
    Jfxu(5,15) = Ts * b3;
    // dx8/du
    Jfxu(7,12) = Ts * (cos(x(0))*sin(x(2))*cos(x(4)) + sin(x(0))*sin(x(4))) / M;
    // dx10/du
    Jfxu(9,12) = Ts * (cos(x(0))*sin(x(2))*sin(x(4)) - sin(x(0))*cos(x(4))) / M;
    // dx12/du
    Jfxu(11,12) = Ts * (cos(x(0))*cos(x(2))) / (M-g);

  }

  void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p,
                        rvec grad_fxu_p) const {
    alpaqa::ScopedMallocAllower ma;

    mat Jfxu(nx,nx+nu); Jfxu.setZero();
    
    // "self"-derivatives
    Jfxu.setIdentity();
    
    // dx1/dx
    Jfxu(0,1) = Ts;
    // dx2/dx
    Jfxu(1,3) = Ts * (a1*x(5));   
    Jfxu(1,5) = Ts * (a1*x(3));
    // dx3/dx
    Jfxu(2,3) = Ts;
    // dx4/dx
    Jfxu(3,1) = Ts * (a2*x(5));    
    Jfxu(3,5) = Ts * (a2*x(1));
    // dx5/dx
    Jfxu(4,5) = Ts;
    // dx6/dx
    Jfxu(5,1) = Ts * (a3*x(3));    
    Jfxu(5,3) = Ts * (a3*x(1));
    // dx7/dx
    Jfxu(6,7) = Ts;
    // dx8/dx
    Jfxu(7,0) = Ts * (-sin(x(0))*sin(x(2))*cos(x(4)) + cos(x(0))*sin(x(4))) * u(0) / M;    
    Jfxu(7,2) = Ts * (cos(x(0))*cos(x(2))*cos(x(4))) * u(0) / M;    
    Jfxu(7,4) = Ts * (-cos(x(0))*sin(x(4))*sin(x(2)) + sin(x(0))*cos(x(4))) * u(0) / M;
    // dx9/dx
    Jfxu(8,9) = Ts;
    // dx10/dx
    Jfxu(9,0) = Ts * (-sin(x(0))*sin(x(2))*sin(x(4)) - cos(x(0))*cos(x(4))) * u(0) / M;    
    Jfxu(9,2) = Ts * (cos(x(2))*cos(x(0))*sin(x(4))) * u(0) / M;    
    Jfxu(9,4) = Ts * (cos(x(4))*cos(x(0))*sin(x(2)) + sin(x(4))*sin(x(0))) * u(0) / M;
    // dx11/dx
    Jfxu(10,11) = Ts;
    // dx12/dx
    Jfxu(11,0) = Ts * (-sin(x(0))*cos(x(2))) * u(0) / (M-g);    
    Jfxu(11,2) = Ts * (-sin(x(2))*cos(x(0))) * u(0) / (M-g);   
    
    // dx2/du
    Jfxu(1,13) = Ts * b1;
    // dx4/du
    Jfxu(3,14) = Ts * b2;
    // dx6/du
    Jfxu(5,15) = Ts * b3;
    // dx8/du
    Jfxu(7,12) = Ts * (cos(x(0))*sin(x(2))*cos(x(4)) + sin(x(0))*sin(x(4))) / M;
    // dx10/du
    Jfxu(9,12) = Ts * (cos(x(0))*sin(x(2))*sin(x(4)) - sin(x(0))*cos(x(4))) / M;
    // dx12/du
    Jfxu(11,12) = Ts * (cos(x(0))*cos(x(2))) / (M-g);

    grad_fxu_p.noalias() = Jfxu.transpose() * p;
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
        Q += 10 * mat::Identity(nx, nx);;
    }
    void eval_add_R_masked([[maybe_unused]] index_t timestep,
                           [[maybe_unused]] crvec xu, 
                           [[maybe_unused]] crvec h, crindexvec mask,
                           rmat R, 
                           [[maybe_unused]] rvec work) const {
        alpaqa::ScopedMallocAllower ma;
        const auto n = mask.size();
        R.noalias() += 0.005 * mat::Identity(n, n);
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


//dynamic model described in https://www.codeproject.com/Articles/5335232/Optimal-Control-of-a-Quadcopter#l2

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

// Dynamics
template <typename X, typename J, typename Params>
void dynamics(index_t &k, X &xu, J &fxu, const Params &params){
    fxu(0) = xu(1);
    fxu(1) = (params.a1*xu(3)*xu(5) + params.b1*xu(13));
    fxu(2) = xu(3);    
    fxu(3) = (params.a2*xu(1)*xu(5) + params.b2*xu(14));
    fxu(4) = xu(5);
    fxu(5) = (params.a3*xu(1)*xu(3) + params.b1*xu(15));
    fxu(6) = xu(7);
    fxu(7) = (cos(xu(0))*sin(xu(2))*cos(xu(4)) + sin(xu(0))*sin(xu(4))) * xu(12) / params.M;
    fxu(8) = xu(9);
    fxu(9) = (cos(xu(0))*sin(xu(2))*sin(xu(4)) - sin(xu(0))*cos(xu(4))) * xu(12) / params.M;
    fxu(10) = xu(11);
    fxu(11) = (cos(xu(0))*cos(xu(2))) * xu(12) / (params.M-params.g);
};

// Horizon length, number of states and input for compile time allocation of AD views
const int p_N = 30, p_nx = 12, p_nu = 4, p_nh = p_nx + p_nu; 

using AD_obj = AD<p_nx+p_nu,p_nx,p_nh>;

// Create the Control Problem object to pass on to the Solver
struct QuadcopterAD{

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

  };

  Params params;

unsigned long int n_seed = 3;
   
  // QR matrix diagonal:
  vec Q, R, QR; 

  // Automatic Differentiation object:
  AD_obj ad_obj[p_N]; 

  QuadcopterAD() : Q(params.nx), R(params.nu), QR(params.nx+params.nu) {
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
    U.lowerbound.setConstant(-1);
    U.upperbound.setConstant(+1);
  }

  void get_D(Box &D) const {    
    D.lowerbound.setConstant(-alpaqa::inf<config_t>);
    D.upperbound.setConstant(+alpaqa::inf<config_t>);
  }

  void get_D_N(Box &D) const {}

  void get_x_init(rvec x_init) const {
    // Single-Shooting
    if (x_init.size() <= params.nu*params.N){
      x_init.setConstant(0.0);
      x_init(6) = 0.3;
      x_init(8) = 0.17;
      x_init(10) = 0.2;
    }
    // Multiple-Shooting
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