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

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

// discrete (Fwd Euler) dynamics function which accepts the AD objects as arguments [in/out]
template <typename X, typename J>
void dynamics(index_t &k, X &xu, J &fxu, const length_t dim, const length_t N_balls,
              const length_t &nx, const real_t &Ts, const real_t &D, const real_t &L, 
              const real_t &m, const real_t &g, const real_t &x_ref, const real_t &y_ref,
              const real_t &z_ref){    

    // ẋₖ = vₖ
    for (size_t i = 0; i < (N_balls-1)*dim-dim; ++i){
      fxu(i) = xu(i) + Ts * xu(i+(N_balls-1)*dim);
    } 
    // ẋ₀ = 0 OR x₀ = p₀
    fxu((N_balls-1)*dim-3) = xu((N_balls-1)*dim-3) + Ts * xu(xu.size()-3); 
    fxu((N_balls-1)*dim-2) = xu((N_balls-1)*dim-2) + Ts * xu(xu.size()-2);
    fxu((N_balls-1)*dim-1) = xu((N_balls-1)*dim-1) + Ts * xu(xu.size()-1);

    // v̇ₖ = (1/m) * (Fₖ₊₁ - Fₖ₋₁) + g
    for (size_t i = 0; i < (N_balls-2)*dim; i += size_t(dim)){
      if (i == 0){
        fxu(i+(N_balls-1)*dim) = xu(i+(N_balls-1)*dim) + Ts * ((D/m)*((1-(L/(std::sqrt(std::pow(xu(i+dim)-xu(i),2)+std::pow(xu(1+i+dim)-xu(1+i),2)+std::pow(xu(2+i+dim)-xu(2+i),2)))))*(xu(i+dim)-xu(i))
                                                              -(1-(L/(std::sqrt(std::pow(xu(i)-x_ref,2)+std::pow(xu(1+i)-y_ref,2)+std::pow(xu(2+i)-z_ref,2)))))*(xu(i)-x_ref)));
        fxu(i+1+(N_balls-1)*dim) = xu(i+1+(N_balls-1)*dim) + Ts * ((D/m)*((1-(L/(std::sqrt(std::pow(xu(i+dim)-xu(i),2)+std::pow(xu(1+i+dim)-xu(1+i),2)+std::pow(xu(2+i+dim)-xu(2+i),2)))))*(xu(1+i+dim)-xu(1+i))
                                                              -(1-(L/(std::sqrt(std::pow(xu(i)-x_ref,2)+std::pow(xu(1+i)-y_ref,2)+std::pow(xu(2+i)-z_ref,2)))))*(xu(1+i)-y_ref)));
        fxu(i+2+(N_balls-1)*dim) = xu(i+2+(N_balls-1)*dim) + Ts * ((D/m)*((1-(L/(std::sqrt(std::pow(xu(i+dim)-xu(i),2)+std::pow(xu(1+i+dim)-xu(1+i),2)+std::pow(xu(2+i+dim)-xu(2+i),2)))))*(xu(2+i+dim)-xu(2+i))
                                                              -(1-(L/(std::sqrt(std::pow(xu(i)-x_ref,2)+std::pow(xu(1+i)-y_ref,2)+std::pow(xu(2+i)-z_ref,2)))))*(xu(2+i)-z_ref))+g);
      }
      else {
        fxu(i+(N_balls-1)*dim) = xu(i+(N_balls-1)*dim) + Ts * ((D/m)*((1-(L/(std::sqrt(std::pow(xu(i+dim)-xu(i),2)+std::pow(xu(1+i+dim)-xu(1+i),2)+std::pow(xu(2+i+dim)-xu(2+i),2)))))*(xu(i+dim)-xu(i))
                                                              -(1-(L/(std::sqrt(std::pow(xu(i)-xu(i-dim),2)+std::pow(xu(1+i)-xu(1+i-dim),2)+std::pow(xu(2+i)-xu(2+i-dim),2)))))*(xu(i)-xu(i-dim))));
        fxu(i+1+(N_balls-1)*dim) = xu(i+1+(N_balls-1)*dim) + Ts * ((D/m)*((1-(L/(std::sqrt(std::pow(xu(i+dim)-xu(i),2)+std::pow(xu(1+i+dim)-xu(1+i),2)+std::pow(xu(2+i+dim)-xu(2+i),2)))))*(xu(1+i+dim)-xu(1+i))
                                                              -(1-(L/(std::sqrt(std::pow(xu(i)-xu(i-dim),2)+std::pow(xu(1+i)-xu(1+i-dim),2)+std::pow(xu(2+i)-xu(2+i-dim),2)))))*(xu(1+i)-xu(1+i-dim))));
        fxu(i+2+(N_balls-1)*dim) = xu(i+2+(N_balls-1)*dim) + Ts * ((D/m)*((1-(L/(std::sqrt(std::pow(xu(i+dim)-xu(i),2)+std::pow(xu(1+i+dim)-xu(1+i),2)+std::pow(xu(2+i+dim)-xu(2+i),2)))))*(xu(2+i+dim)-xu(2+i))
                                                              -(1-(L/(std::sqrt(std::pow(xu(i)-xu(i-dim),2)+std::pow(xu(1+i)-xu(1+i-dim),2)+std::pow(xu(2+i)-xu(2+i-dim),2)))))*(xu(2+i)-xu(2+i-dim)))+g);
      }
    } 
}; 

// Horizon length, number of states and input for compile time allocation of AD views
const int p_N = 10,
          p_Nb = 4,
          p_dim = 3,
          p_nx = ((p_Nb - 1) * p_dim * 2) - p_dim, //< ((N_balls - 1) * dim * 2) - dim
          p_nu = 3,
          p_nh = p_nx + p_nu; 

const real_t p_Ts = 0.0001;

// Define the AD object to use in the computations
using AD_obj = AD<p_nx+p_nu,p_nx,p_nh>;

// Create the struct of the Hanging Chain problem
struct HangingChain{

  USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

  using Box = alpaqa::Box<config_t>;

  // Hanging chain parameters:
  real_t β          = 25;       ///< Q-value for free-end 
  real_t γ          = 1;        ///< Q-value
  real_t δ          = 0.01;     ///< R-value
  real_t m          = 0.03;     ///< mass
  real_t D          = 1.0;      ///< spring constant
  real_t L          = 0.033;    ///< spring length
  real_t v_max      = 0.5;      ///< maximum actuator velocity
  real_t g          = -9.81;    ///< Gravitational acceleration       [m/s²]
  length_t N_balls  = p_Nb;     ///< Number of balls + fixed and free-end
  length_t dim      = p_dim;    ///< Dimensions
  real_t x_ref      =(N_balls)*L;      ///< Desired position of x_N
  real_t y_ref      = 0.0;      ///< Desired position of y_N
  real_t z_ref      = 0.0;      ///< Desired position of z_N
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
  real_t    Ts = p_Ts;          ///< Sampling time
  
  AD_obj ad_obj[p_N];           ///< AD objects for each stage (to avoid data race)

  vec p_ref,                    ///< desired position of x_N, y_N and z_N
         po,                    ///< fixed-end position
          Q,                    ///< Diagonal of matrix Q
          R,                    ///< Diagonal of matrix R
         QR;                    ///< Diagonals Q and R merged

  HangingChain() : p_ref(dim), po(dim), Q(nx), R(nu), QR(nx+nu) {
    
    p_ref << x_ref, y_ref, z_ref;
    po << xo, yo, zo;

    std::cout<<"The desired position is: "<<p_ref.transpose()<<std::endl;
    
    Q.segment(0,(N_balls-1)*dim).setConstant(0); 
    Q((N_balls-1)*dim-1) = β; Q((N_balls-1)*dim-2) = β; Q((N_balls-1)*dim-3) = β;
    Q.segment((N_balls-1)*dim,nx-(N_balls-1)*dim).setConstant(γ);
    
    R.setConstant(δ);
    
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
    U.lowerbound.setConstant(-v_max);
    U.upperbound.setConstant(+v_max);
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
        if (k*3 >= nx/2){
          x_init(i) = 0.0;
          x_init(i+1) = 0.0;
          x_init(i+2) = 0.0;
        } 
        else{
          x_init(i) = (k+1)*L;
          x_init(i+1) = 0;
          //x_init(i+2) = a*std::cosh((x_init(i)-b)/a) - a*std::cosh(-b/a);
          x_init(i+2) = 0;
        }
        ++k;
        if (k*3 == real_t(p_nh)){
          k = 0;
        }
      }
      length_t p_index = (N_balls-1)*dim-dim;
      std::cout<<"The starting point of the free-end is: "<<x_init.segment(p_index,3).transpose()<<std::endl;
    }
  }

  void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const { 
    // alpaqa::ScopedMallocAllower ma;

    auto &&p = x.segment(0, (N_balls-1)*dim);
    auto &&v = x.segment(p.size(),nx-p.size());

    // ẋₖ = vₖ
    for (size_t i = 0; i < (N_balls-1)*dim-dim; ++i){
      fxu(i) = p(i) + Ts * v(i);
    } 
    // ẋ₀ = 0 OR x₀ = p₀
    fxu((N_balls-1)*dim-3) = p((N_balls-1)*dim-3) + Ts * u(0); 
    fxu((N_balls-1)*dim-2) = p((N_balls-1)*dim-2) + Ts * u(1);
    fxu((N_balls-1)*dim-1) = p((N_balls-1)*dim-1) + Ts * u(2);

    // v̇ₖ = (1/m) * (Fₖ₊₁ - Fₖ₋₁) + g
    for (size_t i = 0; i < (N_balls-2)*dim; i += size_t(dim)){
      if (i == 0){
        fxu(i+(N_balls-1)*dim) = v(i) + Ts * ((D/m)*((1-(L/(std::sqrt(std::pow(p(i+dim)-p(i),2)+std::pow(p(1+i+dim)-p(1+i),2)+std::pow(p(2+i+dim)-p(2+i),2)))))*(p(i+dim)-p(i))
                                            -(1-(L/(std::sqrt(std::pow(p(i)-xo,2)+std::pow(p(1+i)-yo,2)+std::pow(p(2+i)-zo,2)))))*(p(i)-xo)));
        fxu(i+1+(N_balls-1)*dim) = v(i+1) + Ts * ((D/m)*((1-(L/(std::sqrt(std::pow(p(i+dim)-p(i),2)+std::pow(p(1+i+dim)-p(1+i),2)+std::pow(p(2+i+dim)-p(2+i),2)))))*(p(1+i+dim)-p(1+i))
                                            -(1-(L/(std::sqrt(std::pow(p(i)-xo,2)+std::pow(p(1+i)-yo,2)+std::pow(p(2+i)-zo,2)))))*(p(1+i)-yo)));
        fxu(i+2+(N_balls-1)*dim) = v(i+2) + Ts * ((D/m)*((1-(L/(std::sqrt(std::pow(p(i+dim)-p(i),2)+std::pow(p(1+i+dim)-p(1+i),2)+std::pow(p(2+i+dim)-p(2+i),2)))))*(p(2+i+dim)-p(2+i))
                                            -(1-(L/(std::sqrt(std::pow(p(i)-xo,2)+std::pow(p(1+i)-yo,2)+std::pow(p(2+i)-zo,2)))))*(p(2+i)-zo))+g);
      }
      else {
        fxu(i+(N_balls-1)*dim) = v(i) + Ts * ((D/m)*((1-(L/(std::sqrt(std::pow(p(i+dim)-p(i),2)+std::pow(p(1+i+dim)-p(1+i),2)+std::pow(p(2+i+dim)-p(2+i),2)))))*(p(i+dim)-p(i))
                                            -(1-(L/(std::sqrt(std::pow(p(i)-p(i-dim),2)+std::pow(p(1+i)-p(1+i-dim),2)+std::pow(p(2+i)-p(2+i-dim),2)))))*(p(i)-p(i-dim))));
        fxu(i+1+(N_balls-1)*dim) = v(i+1) + Ts * ((D/m)*((1-(L/(std::sqrt(std::pow(p(i+dim)-p(i),2)+std::pow(p(1+i+dim)-p(1+i),2)+std::pow(p(2+i+dim)-p(2+i),2)))))*(p(1+i+dim)-p(1+i))
                                            -(1-(L/(std::sqrt(std::pow(p(i)-p(i-dim),2)+std::pow(p(1+i)-p(1+i-dim),2)+std::pow(p(2+i)-p(2+i-dim),2)))))*(p(1+i)-p(1+i-dim))));
        fxu(i+2+(N_balls-1)*dim) = v(i+2) + Ts * ((D/m)*((1-(L/(std::sqrt(std::pow(p(i+dim)-p(i),2)+std::pow(p(1+i+dim)-p(1+i),2)+std::pow(p(2+i+dim)-p(2+i),2)))))*(p(2+i+dim)-p(2+i))
                                            -(1-(L/(std::sqrt(std::pow(p(i)-p(i-dim),2)+std::pow(p(1+i)-p(1+i-dim),2)+std::pow(p(2+i)-p(2+i-dim),2)))))*(p(2+i)-p(2+i-dim)))+g);
      }
    } 

  } // *discretized using Explicit Euler

  void eval_jac_f(index_t timestep, crvec x, crvec u, rmat Jfxu) const {
    // alpaqa::ScopedMallocAllower ma;
    assign_values_xu<p_nx+p_nu,p_nx>(x, u, ad_obj[timestep]);
    dynamics(timestep, ad_obj[timestep].xu_fad, ad_obj[timestep].fxu_fad, dim, N_balls, 
              nx, Ts, D, L, m, g, xo, yo, zo);
    assign_values<p_nx+p_nu,p_nx>(Jfxu, ad_obj[timestep]);
  }

  void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p,
                        rvec grad_fxu_p) const {
    // alpaqa::ScopedMallocAllower ma;
    assign_values_xu<p_nx+p_nu,p_nx>(x, u, ad_obj[timestep]);
    dynamics(timestep, ad_obj[timestep].xu_fad, ad_obj[timestep].fxu_fad, dim, N_balls,
              nx, Ts, D, L, m, g, xo, yo, zo);
    assign_values<p_nx+p_nu,p_nx> (grad_fxu_p, p, ad_obj[timestep]);
  }

  void eval_h([[maybe_unused]] index_t timestep, crvec x, crvec u, rvec h) const {
      // alpaqa::ScopedMallocAllower ma;
      h.topRows(nx)    = x;
      h.bottomRows(nu) = u;
      h((N_balls-1)*dim-3) = x((N_balls-1)*dim-3) - x_ref;
      h((N_balls-1)*dim-2) = x((N_balls-1)*dim-2) - y_ref;
      h((N_balls-1)*dim-1) = x((N_balls-1)*dim-1) - z_ref;
  }
  void eval_h_N(crvec x, rvec h) const {       
      h = x;
      h((N_balls-1)*dim-3) = x((N_balls-1)*dim-3) - x_ref;
      h((N_balls-1)*dim-2) = x((N_balls-1)*dim-2) - y_ref;
      h((N_balls-1)*dim-1) = x((N_balls-1)*dim-1) - z_ref; 
  }
  [[nodiscard]] real_t eval_l([[maybe_unused]] index_t timestep, crvec h) const {
      // alpaqa::ScopedMallocAllower ma;
      return h.transpose() * QR.asDiagonal() * h;
  }
  [[nodiscard]] real_t eval_l_N(crvec h) const {
      // alpaqa::ScopedMallocAllower ma;
      return  h.transpose() * Q.asDiagonal() * h;
  }
  void eval_qr([[maybe_unused]] index_t timestep, 
                [[maybe_unused]] crvec xu, crvec h, rvec qr) const {
      // alpaqa::ScopedMallocAllower ma;
      auto Jh_xu    = mat::Identity(nx + nu, nx + nu);
      auto &&grad_l = QR.asDiagonal() * h;
      qr            = Jh_xu.transpose() * grad_l;
  }
  void eval_q_N([[maybe_unused]] crvec x, crvec h, rvec q) const {
      // alpaqa::ScopedMallocAllower ma;
      auto Jh_x     = mat::Identity(nx, nx);
      auto &&grad_l = 2 * Q.asDiagonal() * h;
      q             = Jh_x.transpose() * grad_l;
  }
  void eval_add_Q([[maybe_unused]] index_t timestep, 
                  [[maybe_unused]] crvec xu, 
                  [[maybe_unused]] crvec h, rmat Q) const {
      Q += mat::Identity(nx, nx);
  }
  void eval_add_Q_N([[maybe_unused]] crvec x,
                    [[maybe_unused]] crvec h, rmat Q) const {
      // alpaqa::ScopedMallocAllower ma;
      Q += 10 * mat::Identity(nx, nx);
  }
  void eval_add_R_masked([[maybe_unused]] index_t timestep,
                          [[maybe_unused]] crvec xu, 
                          [[maybe_unused]] crvec h, crindexvec mask,
                          rmat R, 
                          [[maybe_unused]] rvec work) const {
      // alpaqa::ScopedMallocAllower ma;
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
      // alpaqa::ScopedMallocAllower ma;
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

