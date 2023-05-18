#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/problem/box.hpp>

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

using Box = alpaqa::Box<config_t>;

auto eval_proj_set (const Box &box, crvec x) {
    using binary_real_f = real_t (*)(real_t, real_t);
    return x.binaryExpr(box.lowerbound, binary_real_f(std::fmax))
            .binaryExpr(box.upperbound, binary_real_f(std::fmin));
};

template <typename Iterate, typename Problem>
void eval_iterate_fun (int k, Iterate &It, Problem &problem, crvec μ, crvec y, Box &F){
    
    auto nx = problem.get_nx();
    auto nu = problem.get_nu();
    auto nxu = problem.get_nx() + problem.get_nu(); 
    auto N = problem.get_N() + 1;

    Eigen::Index k_ = k;

    if (k_ == N-1){
        problem.eval_h_N(It.xu.segment(k*nxu,nx),
                        It.hxu.segment(k*nxu,nx));
        It.lxu(k) = problem.eval_l_N(It.hxu.segment(k*nxu,nx));
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
        It.lxu(k) = problem.eval_l(k, It.hxu.segment(k*nxu,nxu));
        problem.eval_qr(k, It.xu.segment(k*nxu,nxu), 
                        It.hxu.segment(k*nxu,nxu),
                        It.qr.segment(k*nxu,nxu));
        It.g.segment(k*nx,nx) = It.fxu.segment(k*nx,nx) - It.xu.segment((k+1)*nxu,nx);           
        It.gz.segment(k*nx,nx) = It.g.segment(k*nx,nx) + 
                                ((μ).segment(k*nx,nx).asDiagonal().inverse()*(y.segment(k*nx,nx)));            
        It.Π_D.segment(k*nx,nx).noalias() = eval_proj_set(F, It.gz.segment(k*nx,nx));
        It.Z.segment(k*nx,nx)   = (μ.segment(k*nx,nx).cwiseProduct(It.fxu.segment(k*nx,nx)
                                -It.xu.segment((k+1)*nxu,nx)-It.Π_D.segment(k*nx,nx)) 
                                + y.segment(k*nx,nx));          
    }
} 

template <typename Iterate, typename Problem>
void eval_grad_ψ_k_fun (int k, Iterate &It, Problem &problem, crvec μ, crvec y){

    auto nx = problem.get_nx();
    auto nu = problem.get_nu();
    auto nxu = problem.get_nx() + problem.get_nu(); 
    auto N = problem.get_N() + 1;

    Eigen::Index k_ = k;

    if (k_ == N-1){
        It.grad_ψ.segment(k*nxu,nx) = It.qr.segment(k*nxu,nx) - It.Z.segment((k-1)*nx,nx);                            
    } 
    else if(k_ == 0){
        It.grad_ψ.segment(nx,nu) = It.qr.segment(nx,nu) +
                    It.Jfxu.block(0,nx,nx,nu).transpose() * It.Z.segment(k*nx,nx);
        It.grad_ψ.segment(0,nx).setZero();
    }
    else {
        It.grad_ψ.segment(k*nxu,nxu) = It.qr.segment(k*nxu,nxu) - mat::Identity(nxu, nx)*It.Z.segment((k-1)*nx,nx) +
                    It.Jfxu.block(k*nx,0,nx,nxu).transpose() * It.Z.segment(k*nx,nx);
    }
}

template <typename Iterate, typename Problem>
real_t eval_ψ_k_fun (int k, Iterate &It, Problem &problem, crvec μ, Box &F){
    
    auto nx = problem.get_nx();
    auto nxu = problem.get_nx() + problem.get_nu(); 
    auto N = problem.get_N() + 1;

    Eigen::Index k_ = k;
    real_t ψxu_k = 0;

    if (k_ == N-1){
        return ψxu_k = It.lxu(k);
    }
    else {
        auto d = (It.gz.segment(k*nx,nx)-
                eval_proj_set(F,It.gz.segment(k*nx,nx)));
        return ψxu_k = It.lxu(k) + 
                    0.5*d.transpose()*(μ).segment(k*nx,nx).asDiagonal()*d;
    }
} 

template <typename Iterate, typename Problem>
real_t eval_ψ_hat_k_fun (int k, Iterate &It, Problem &problem, crvec μ, crvec y, Box &F){
    
    auto nx = problem.get_nx();
    auto nu = problem.get_nu();
    auto nxu = problem.get_nx() + problem.get_nu(); 
    auto N = problem.get_N() + 1;

    real_t ψxû_k = 0;
    Eigen::Index i_ = k;

    if (i_ == N-1){
        problem.eval_h_N(It.xû.segment(k*nxu,nx),
                        It.hxû.segment(k*nxu,nx));
        It.lxû(k) = problem.eval_l_N(It.hxû.segment(k*nxu,nx));
        return ψxû_k = It.lxû(k);
    } 
    else {
        problem.eval_f(k, It.xû.segment(k*nxu,nx),
                    It.xû.segment(k*(nxu)+nx,nu), 
                    It.fxû.segment(k*nx,nx));
        problem.eval_h(k, It.xû.segment(k*nxu,nx), 
                        It.xû.segment((k*nxu)+nx,nu),
                        It.hxû.segment(k*nxu,nxu));
        It.lxû(k) = problem.eval_l(k, It.hxû.segment(k*nxu,nxu));                       
        It.gz_hat.segment(k*nx,nx) = It.fxû.segment(k*nx,nx) - 
                        It.xû.segment((k+1)*nxu,nx)  
                        + ((μ).segment(k*nx,nx).asDiagonal().inverse()*(y.segment(k*nx,nx)));
        It.Π_D_hat.segment(k*nx,nx) = eval_proj_set(F, It.gz_hat.segment(k*nx,nx));  
        auto d = (It.gz_hat.segment(k*nx,nx)-
                eval_proj_set(F,It.gz_hat.segment(k*nx,nx))); 
        return ψxû_k = It.lxû(k) + 
                    0.5*d.transpose()*(μ).segment(k*nx,nx).asDiagonal()*d;            
    }
}

auto eval_proj_grad_step_box (const Box &box, real_t γ, crvec x, crvec grad_ψ, rvec x̂,
                                    rvec p) {
    using binary_real_f = real_t (*)(real_t, real_t);
    p                   = (-γ * grad_ψ)
            .binaryExpr(box.lowerbound - x, binary_real_f(std::fmax))
            .binaryExpr(box.upperbound - x, binary_real_f(std::fmin));
    x̂ = x + p;
};

template <typename Iterate, typename Problem>
void eval_prox_fun (int k, Iterate &It, Problem &problem, crvec μ, crvec y, Box &D_N, Box &D, Box &U){
    
    auto nx = problem.get_nx();
    auto nu = problem.get_nu();
    auto nxu = problem.get_nx() + problem.get_nu(); 
    auto N = problem.get_N() + 1;

    Eigen::Index k_ = k;

    if (k_ == N-1){
        eval_proj_grad_step_box(D_N, It.γ, It.xu.segment(k*nxu,nx), 
                                It.grad_ψ.segment(k*nxu,nx), 
                                It.xû.segment(k*nxu,nx), It.p.segment(k*nxu,nx));
    } 
    else {
        eval_proj_grad_step_box(D, It.γ, It.xu.segment(k*nxu,nx), 
                                It.grad_ψ.segment(k*nxu,nx), 
                                It.xû.segment(k*nxu,nx), It.p.segment(k*nxu,nx));
        eval_proj_grad_step_box(U, It.γ, It.xu.segment((k*nxu)+nx,nu), 
                                It.grad_ψ.segment((k*nxu)+nx,nu), 
                                It.xû.segment((k*nxu)+nx,nu), It.p.segment((k*nxu)+nx,nu));
    }
}

template <typename Iterate, typename Problem>
void eval_GN_accelerator_fun (int k, Iterate &It, Problem &problem, crvec μ){

    auto nx = problem.get_nx();
    auto nxu = problem.get_nx() + problem.get_nu(); 
    auto N = problem.get_N() + 1;

    Eigen::Index k_ = k;

    mat work (nx,nxu+nx);

    if (k_ == N-1){
        mat work_eval_Q(nx,nx);
        problem.eval_add_Q_N(It.xu.segment(k*nxu,nx),
                        It.hxu.segment(k*nxu,nx), 
                        work_eval_Q);
        It.GN.block(nxu*k,nxu*k,nx,nx) += work_eval_Q;
    } 
    else if (k_ == 0){
        mat work_eval_Q(nxu,nxu);
        work.leftCols(nxu) = It.Jfxu.block(k*nx,0,nx,nxu);
        work.rightCols(nx).setIdentity();
        work.rightCols(nx) *= -1; 
        problem.eval_add_Q(k, It.xu.segment(k*nxu,nxu), 
                            It.hxu.segment(k*nxu,nxu),
                            work_eval_Q); // should find another way of evaluating it, other than using eval_add_Q and eval_add_Q_N
        It.GN.block(k*nxu,k*nxu,nxu,nxu) += work_eval_Q;
        It.GN.block(k*nxu,k*nxu,nxu+nx,nxu+nx) += (work.transpose() * 
                                        (μ).segment(k*nx,nx).asDiagonal()) *
                                        (work);
    } 
    else {
        mat work_eval_Q(nxu,nxu);
        work.leftCols(nxu) = It.Jfxu.block(k*nx,0,nx,nxu);
        work.rightCols(nx).setIdentity();
        work.rightCols(nx) *= -1;
        problem.eval_add_Q(k, It.xu.segment(k*nxu,nxu), 
                            It.hxu.segment(k*nxu,nxu),
                            work_eval_Q);
        It.GN.block(k*nxu,k*nxu,nxu,nxu) += work_eval_Q;
        It.GN.block(k*nxu,k*nxu,nxu+nx,nxu+nx) += (work.transpose() * 
                                        (μ).segment(k*nx,nx).asDiagonal()) *
                                        (work);
    }
}

template <typename Iterate, typename Problem>
void fd_grad_ψ_fun (int k, Iterate &It, Problem &problem, crvec μ, crvec y, Box &F){
    alpaqa::ScopedMallocAllower ma;
        
    real_t h = 0.001;
    real_t ψxu_ = 0;
    index_t N   = problem.get_N(); 
    index_t nu  = problem.get_nu();
    index_t nx  = problem.get_nx();
    index_t nxu = nx + nu;
    index_t n   = (N-1)*(nxu) + nx;
    index_t m   = (N-1)*nx;
    vec lxu_    = vec::Zero(N);
    vec hxu_    = vec::Zero(n);
    vec xu_     = vec::Zero(n); 
    vec fxu_    = vec::Zero(m);
    vec d       = vec::Zero(m);
    vec v       = vec::Zero(m);
    vec g       = vec::Zero(m);
    Box D       = alpaqa::Box<config_t>::NaN(m);
    problem.get_D(D);

    for (index_t i = 0; i < n; ++i) {
        xu_ = It.xu;
        xu_(i) += h;
        for (index_t k = 0; k < N; ++k) {
            if (k == N-1) {
                problem.eval_h_N(xu_.segment(k*nxu,nx),
                                    hxu_.segment(k*nxu,nx));
                lxu_(k) = problem.eval_l_N(hxu_.segment(k*nxu,nx));
            } else {
                problem.eval_f(k, xu_.segment(k*nxu,nx),
                                xu_.segment(k*(nxu)+nx,nu), 
                                fxu_.segment(k*nx,nx));
                problem.eval_h(k, xu_.segment(k*nxu,nx), 
                                xu_.segment((k*nxu)+nx,nu),
                                hxu_.segment(k*nxu,nxu));
                lxu_(k) = problem.eval_l(k, hxu_.segment(k*nxu,nxu)); 
            }
        }
        for (index_t k = 0; k < N-1; ++k) {
            g.segment(k*nx,nx) = fxu_.segment(k*nx,nx) - 
                                    xu_.segment((k+1)*nxu,nx);
        }
        d.noalias() = g + μ.asDiagonal().inverse()*y; 
        v.noalias() = d - eval_proj_set(F, d); 
        ψxu_ = lxu_.sum() + real_t(.5)*v.transpose()*(μ).asDiagonal()*v;
        It.grad_ψ(i) = (ψxu_ - It.ψxu)/h;
    }
}