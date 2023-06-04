#pragma once

#include <alpaqa/problem/kkt-error.hpp>
#include <alpaqa/problem/ocproblem.hpp>

#include <Kokkos_Core.hpp>

namespace alpaqa{

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

struct Iterate2 {
    
    vec xû;         //< Inputs u interleaved with states x after prox grad
    vec xu;         //< Inputs u interleaved with states x -> x₀, u₀, x₁, u₁,...,uₙ₋₁, xₙ 
    vec grad_ψ;     //< Gradient of cost w.r.t {x₀, u₀, x₁, u₁,...,uₙ₋₁, xₙ}
    vec fxu;        //< Dynamics f(x,u)
    vec hxu;        //< nonconvex mapping h(x,u)
    vec fxû;        //< Dynamics f(x,û)
    vec hxû;        //< nonconvex mapping h(x,û)
    vec lxu;        //< cost function at iterate -> l(x,u)
    vec lxû;        //< cost function after T(x,u) -> l(x,û)
    mat Jfxu;       //< Jacobian of dynamics Jf(x,u)
    mat GN;         //< GN approximation of ∇²ψ
    vec qr;         //< Cost function gradient Jh(x,u)dl(x,u)
    vec g;          //< Dynamic constraints, f(x,u) - x+
    vec gz;         //< g(x) + Σ⁻¹y
    vec gz_hat;     //< g(x_hat) + Σ⁻¹y
    vec Π_D;        //< Projection of states and variables into their box constraints
    vec Π_D_hat;    //< Projection of states and variables into their box constraints
    vec Z;          //< (Σ*(f(xₖ,uₖ)-xₖ₊₁-Π_D(f(xₖ,uₖ)-xₖ₊₁+Σ⁻¹y))-y)
    vec p;          //< Proximal gradient step
    
    real_t ψxu      = alpaqa::NaN<config_t>;        //< Cost in x
    real_t ψxû      = alpaqa::NaN<config_t>;        //< Cost in x̂
    real_t γ        = alpaqa::NaN<config_t>;        //< Step size γ
    real_t L        = alpaqa::NaN<config_t>;        //< Lipschitz estimate L
    real_t pᵀp      = alpaqa::NaN<config_t>;        //< Norm squared of p
    real_t grad_ψᵀp = alpaqa::NaN<config_t>;        //< Dot product of gradient and p

    Iterate2(length_t n, length_t m, length_t nxu, length_t N) :
        xu(n), xû(n), grad_ψ(n), fxu(m), fxû(m), hxu(n), hxû(n), 
        qr(n), g(m), gz(m), gz_hat(m), Π_D(m), Π_D_hat(m), Z(m), p(n), 
        Jfxu(m, nxu), GN(n,n), lxu(N), lxû(N) {}

};

template <Config Conf, typename timestep>
void eval_grad_L(const alpaqa::TypeErasedControlProblem<Conf> &problem,
                alpaqa::Iterate2 It, alpaqa::crvec x, alpaqa::crvec y, timestep k){
    
    auto nx = problem.get_nx();
    auto nu = problem.get_nu();
    auto nxu = problem.get_nx() + problem.get_nu();
    auto N = problem.get_N() + 1;

    Eigen::Index k_ = k;

    if (k_ == N-1){
        It.grad_ψ.segment(k*nxu,nx) = It.qr.segment(k*nxu,nx) - It.y.segment((k-1)*nx,nx);                            
    } 
    else if(k_ == 0){
        It.grad_ψ.segment(nx,nu) = It.qr.segment(nx,nu) +
                    It.Jfxu.block(0,nx,nx,nu).transpose() * It.y.segment(k*nx,nx);
        It.grad_ψ.segment(0,nx).setZero();
    }
    else {
        It.grad_ψ.segment(k*nxu,nxu) = It.qr.segment(k*nxu,nxu) - mat::Identity(nxu, nx)*It.y.segment((k-1)*nx,nx) +
                    It.Jfxu.block(k*nx,0,nx,nxu).transpose() * It.y.segment(k*nx,nx);
    }
}

template <Config Conf>
KKTError <Conf> compute_kkt_error (const alpaqa::TypeErasedControlProblem<Conf> &problem,
                                    alpaqa::crvec xu, alpaqa::crvec y) {

    USING_ALPAQA_CONFIG(Conf);
    const auto n = xu.size(), m = y.size();
    vec z(n), grad_Lx(n), work(n), g(m), e(m);
    // Gradient of Lagrangian, ∇ℒ(x,y) = ∇f(x) + ∇g(x) y
    problem.eval_grad_L(xu, y, grad_Lx, work);
    // Eliminate normal cone of bound constraints, z = Π(x - ∇ℒ(x,y)) - x
    problem.eval_prox_grad_step(1, xu, grad_Lx, work, z);
    // Stationarity, ‖Π(x - ∇ℒ(x,y)) - x‖
    auto stationarity = alpaqa::vec_util::norm_inf(z);
    // Constraints, g(x)
    problem.eval_g(xu, g);
    // Distance to feasible set, e = g(x) - Π(g(x))
    problem.eval_proj_diff_g(g, e);
    // Constraint violation, ‖g(x) - Π(g(x))‖
    auto constr_violation = alpaqa::vec_util::norm_inf(e);
    // Complementary slackness
    real_t complementarity = std::inner_product(
        y.begin(), y.end(), e.begin(), real_t(0),
        [](real_t acc, real_t ye) { return std::fmax(acc, std::abs(ye)); }, std::multiplies<>{});
    return {.stationarity     = stationarity,
            .constr_violation = constr_violation,
            .complementarity  = complementarity};

}

}