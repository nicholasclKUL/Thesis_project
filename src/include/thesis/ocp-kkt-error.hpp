#pragma once

#include <alpaqa/problem/kkt-error.hpp>
#include <alpaqa/problem/ocproblem.hpp>

#include <Kokkos_Core.hpp>
#include <Sacado.hpp>

#include <thesis/para-panoc-helpers.hpp>

namespace alpaqa{

template <Config Conf>
struct KKTiterate {

    USING_ALPAQA_CONFIG(Conf);
    
    vec xu;         //< Inputs u interleaved with states x -> x₀, u₀, x₁, u₁,...,uₙ₋₁, xₙ 
    vec grad_L;     //< Gradient of cost w.r.t {x₀, u₀, x₁, u₁,...,uₙ₋₁, xₙ}
    vec fxu;        //< Dynamics f(x,u)
    vec hxu;        //< nonconvex mapping h(x,u)
    mat Jfxu;       //< Jacobian of dynamics Jf(x,u)
    vec qr;         //< Cost function gradient Jh(x,u)dl(x,u)
    vec g;          //< Dynamic constraints, f(x,u) - x+
    vec Π_xu;        //< Projection of states into their box constraints

    KKTiterate(length_t n, length_t m, length_t nx,length_t nxu) :
        xu(n), grad_L(n), fxu(m-nx), hxu(n), qr(n), g(m), Π_xu(n), 
        Jfxu(m-nx, nxu) {}

};

template <Config Conf, typename timestep>
void eval_KKTiterate_fun (timestep k, alpaqa::KKTiterate<Conf> &It, 
                        const alpaqa::TypeErasedControlProblem<Conf> &problem){
    
    auto nx = problem.get_nx();
    auto nu = problem.get_nu();
    auto nxu = problem.get_nx() + problem.get_nu(); 
    auto N = problem.get_N() + 1;

    alpaqa::Box<Conf> U(nu), D(nx);
    problem.get_U(U);
    problem.get_D(D);

    Eigen::Index k_ = k;

    vec<Conf> xo(nx);
    problem.get_x_init(xo);

    if (k_ == N-1){
        problem.eval_h_N(It.xu.segment(k*nxu,nx),
                        It.hxu.segment(k*nxu,nx));
        problem.eval_q_N(It.xu.segment(k*nxu,nx),
                        It.hxu.segment(k*nxu,nx), 
                        It.qr.segment(k*nxu,nx));
    } else {
        problem.eval_f(k, It.xu.segment(k*nxu,nx),
                        It.xu.segment((k*nxu)+nx,nu), 
                        It.fxu.segment(k*nx,nx));
        problem.eval_jac_f(k, It.xu.segment(k*nxu,nx),
                        It.xu.segment((k*nxu)+nx,nu),
                        It.Jfxu.block(k*nx,0,nx,nxu));
        problem.eval_h(k, It.xu.segment(k*nxu,nx), 
                        It.xu.segment((k*nxu)+nx,nu),
                        It.hxu.segment(k*nxu,nxu));
        problem.eval_qr(k, It.xu.segment(k*nxu,nxu), 
                        It.hxu.segment(k*nxu,nxu),
                        It.qr.segment(k*nxu,nxu));
        It.g.segment(k*nx,nx) = It.fxu.segment(k*nx,nx) - It.xu.segment((k+1)*nxu,nx);           
        It.Π_xu.segment(k*(nxu),nx).noalias() = eval_proj_set(D, It.xu.segment(k*(nxu),nx));
        It.Π_xu.segment(k*(nxu)+nx,nu).noalias() = eval_proj_set(U, It.xu.segment(k*(nxu)+nx,nu));          
    }
    if (k_ == 0){
        It.g.bottomRows(nx) = xo - It.xu.segment(0,nx);
    }
} 

template <Config Conf, typename timestep>
void eval_grad_L_fun(timestep k, alpaqa::KKTiterate<Conf> &It,
                     const alpaqa::TypeErasedControlProblem<Conf> &problem, 
                     alpaqa::crvec<Conf> y){
    
    USING_ALPAQA_CONFIG(Conf);

    auto nx = problem.get_nx();
    auto nu = problem.get_nu();
    auto nxu = problem.get_nx() + problem.get_nu();
    auto N = problem.get_N() + 1;

    Eigen::Index k_ = k;

    if (k_ == N-1){
        It.grad_L.segment(k*nxu,nx) = It.qr.segment(k*nxu,nx) - y.segment((k-1)*nx,nx);                       
    } 
    else if(k_ == 0){
        It.grad_L.segment(nx,nu) = It.qr.segment(nx,nu) +
                    It.Jfxu.block(0,nx,nx,nu).transpose() * y.segment(k*nx,nx);
        It.grad_L.segment(0,nx) = It.qr.segment(0,nx) - y.bottomRows(nx) +
                    It.Jfxu.block(0,0,nx,nx).transpose() * y.segment(0,nx);
    }
    else {
        It.grad_L.segment(k*nxu,nxu) = It.qr.segment(k*nxu,nxu) 
                                    - alpaqa::mat<Conf>::Identity(nxu, nx)*y.segment((k-1)*nx,nx) +
                                    It.Jfxu.block(k*nx,0,nx,nxu).transpose() * y.segment(k*nx,nx);
    }
}

template <Config Conf>
KKTError<Conf> compute_kkt_error (const alpaqa::TypeErasedControlProblem<Conf> &problem,
                                    alpaqa::KKTiterate<Conf> &It, alpaqa::crvec<Conf> xu, alpaqa::crvec<Conf> y, 
                                    alpaqa::index_t<Conf> nthrds) {
    
    USING_ALPAQA_CONFIG(Conf);

    const auto n = xu.size(), m = y.size(), nx = problem.get_nx(), nu = problem.get_nu(),
               N = problem.get_N();
    vec e(m+n), μ = vec::Zero(m);

    It.xu = xu;

    alpaqa::Box<Conf> F = Box<config_t>::NaN(nx*N),
                      U = Box<config_t>::NaN(nu*N), 
                      D = Box<config_t>::NaN(nx*N+nx);
    F.lowerbound(0.); F.upperbound(0.);
    problem.get_D(D);
    problem.get_U(U);

    // Create the lambdas to use Kokkos
    auto eval_iterate = [&](int k, alpaqa::KKTiterate<Conf> &It){       
        eval_KKTiterate_fun(k, It, problem);
    }; 
    auto eval_grad_L = [&](int k, alpaqa::KKTiterate<Conf> &It, alpaqa::crvec<Conf> y){       
        eval_grad_L_fun(k, It, problem, y);
    };
    // Evaluate iterate
    Kokkos::parallel_for("iterate xu₀", nthrds, [&] (const int k){
        eval_iterate(k, It);
    });
    Kokkos::fence(); 
    // Evaluate gradient of Lagrangian, ∇ℒ(x,y) = ∇f(x) + ∇g(x)y
    Kokkos::parallel_for("iterate xu₀", nthrds, [&] (const int k){
        eval_grad_L(k, It, y);
    });
    auto stationarity = alpaqa::vec_util::norm_inf(It.grad_L);
    // Distance to feasible set, e = g(x) - Π(g(x)),
    e << It.g , (xu - It.Π_xu);
    // Constraint violation, ‖g(x) - Π(g(x))‖
    auto constr_violation = alpaqa::vec_util::norm_inf(e);
    return {.stationarity     = stationarity,
            .constr_violation = constr_violation,
            .complementarity  = alpaqa::NaN<Conf> //no lagrange multipliers for states and inputs
            };

}

template <Config Conf>
KKTError<Conf> compute_kkt_error (const alpaqa::TypeErasedControlProblem<Conf> &problem,
                                    alpaqa::crvec<Conf> u, alpaqa::crvec<Conf> x) {
    
    USING_ALPAQA_CONFIG(Conf);

    const auto nx = problem.get_nx(), nu = problem.get_nu(),
               N = problem.get_N();
    vec e(nx*(N-1)), z(nx+nu), grad_L((nu+nx)*N); 

    for (size_t k = 0; k < N; ++k){
        if (k == N-1){
            problem.eval_q_N(x.bottomRows(nx),x.bottomRows(nx),grad_L.bottomRows(nx));
        }
        else {
            z << x.segment(k*nx,nx),u.segment(k*nu,nu); 
            problem.eval_qr(k,z,z,grad_L.segment(k*(nx+nu),nx+nu));
            problem.eval_f(k,x.segment(k*nx,nx),u.segment(k*nu,nu),e.segment(k*nx,nx));
            e.segment(k*nx,nx) -= x.segment((k+1)*nx,nx);
        }
    }

    auto stationarity = alpaqa::vec_util::norm_inf(grad_L);
    auto constr_violation = alpaqa::vec_util::norm_inf(e);
    return {.stationarity     = stationarity,
            .constr_violation = constr_violation,
            .complementarity  = alpaqa::NaN<Conf> //no lagrange multipliers for states and inputs
            };

}

template <Config Conf>
KKTError<Conf> compute_kkt_error (const alpaqa::TypeErasedControlProblem<Conf> &problem, alpaqa::crvec<Conf> u, 
                                  alpaqa::crvec<Conf> Qk, alpaqa::crvec<Conf> Rk) {
    
    USING_ALPAQA_CONFIG(Conf);

    const auto nx = problem.get_nx(), nu = problem.get_nu(),
               N = problem.get_N()+1;
    real_t fz = 0;
    vec grad_L(nu*(N-1)), x(nx*N), qr((nx+nu)*N-nu), 
        qn(nx), p(nx), pk(nx), xu((nx+nu)*N-nu);
    mat Jfxu(nx,nx+nu); 

    problem.get_x_init(x.segment(0,nx));
    for (size_t i = 0; i < N-1; ++i){
        problem.eval_qr(i, xu.segment(i*(nx+nu),nx+nu), xu.segment(i*(nx+nu),nx+nu), qr.segment(i*(nx+nu),nx+nu));
        fz += (.5) * x.segment(i*nx,nx).transpose() * Qk.asDiagonal() * x.segment(i*nx,nx); 
        fz += (.5) * u.segment(i*nu,nu).transpose() * Rk.asDiagonal() * u.segment(i*nu,nu);
        problem.eval_f(i, x.segment(i*nx,nx), u.segment(i*nu,nu), x.segment((i+1)*nx,nx));
    }
    fz += (.5) * x.bottomRows(nx).transpose() * Qk.asDiagonal() * x.bottomRows(nx);
    p = Qk.asDiagonal() * x.bottomRows(nx);
    for (size_t i = N-2; i > 0; --i){
        problem.eval_jac_f(i, x.segment(i*nx,nx), u.segment(i*nu,nu), Jfxu);
        grad_L.segment(i*nu,nu) = Rk.asDiagonal()*u.segment(i*nu,nu) +
                                  Jfxu.block(0,nx,nx,nu).transpose()*p;
        p = Qk.asDiagonal()*x.segment(i*nx,nx) + Jfxu.block(0,0,nx,nx).transpose()*p;
    }
    grad_L.topRows(nu) = Rk.asDiagonal()*u.topRows(nu) + Jfxu.block(0,nx,nx,nu).transpose()*p;

    std::cout<<grad_L.transpose()<<std::endl;

    auto stationarity = alpaqa::vec_util::norm_inf(grad_L);
    return {.stationarity     = stationarity,
            .constr_violation = alpaqa::NaN<Conf>,
            .complementarity  = alpaqa::NaN<Conf> //no lagrange multipliers for states and inputs
            };
            
}

}