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
    auto nh = problem.get_nh();
    auto nh_N = problem.get_nh_N();

    vec xo(nx);
    problem.get_x_init(xo);

    Eigen::Index k_ = k;

    if (k_ == N-1){
        problem.eval_h_N(It.xu.segment(k*nxu,nx),
                        It.hxu.segment(k*nxu,nx));
        It.lxu(k) = problem.eval_l_N(It.hxu.segment(k*nh,nh_N));
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
        It.lxu(k) = problem.eval_l(k, It.hxu.segment(k*nh,nh));
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
    if (k_ == 0){
        It.g.bottomRows(nx) = xo - It.xu.segment(0,nx);           
        It.gz.bottomRows(nx) = It.g.bottomRows(nx) + 
                                ((μ).bottomRows(nx).asDiagonal().inverse()*(y.bottomRows(nx)));            
        It.Π_D.bottomRows(nx).noalias() = eval_proj_set(F, It.gz.bottomRows(nx));
        It.Z.bottomRows(nx)   = (μ.bottomRows(nx).cwiseProduct(It.g.bottomRows(nx)) 
                                + y.bottomRows(nx));  
    }
} 

// template <typename Problem>
// void eval_grad_f_plus_prod_grad_g_y(int k, Problem &problem, rvec grad_L_or_ψ, mat Jfxu, crvec qr, crvec y){ 
    
//     auto nx = problem.get_nx();
//     auto nu = problem.get_nu();
//     auto nxu = problem.get_nx() + problem.get_nu(); 
//     auto N = problem.get_N() + 1;
    
//     Eigen::Index k_ = k;

//     if (k_ == N-1){
//         grad_L_or_ψ.segment(k*nxu,nx) = qr.segment(k*nxu,nx) - y.segment((k-1)*nx,nx);                            
//     } 
//     else if(k_ == 0){
//         grad_L_or_ψ.segment(nx,nu) = qr.segment(nx,nu) +
//                     Jfxu.block(0,nx,nx,nu).transpose() * y.segment(k*nx,nx);
//         grad_L_or_ψ.segment(0,nx).setZero();
//     }
//     else {
//         grad_L_or_ψ.segment(k*nxu,nxu) = qr.segment(k*nxu,nxu) - mat::Identity(nxu, nx)*y.segment((k-1)*nx,nx) +
//                     Jfxu.block(k*nx,0,nx,nxu).transpose() * y.segment(k*nx,nx);
//     }

// } 

template <typename Iterate, typename Problem>
void eval_grad_ψ_k_fun (int k, Iterate &It, Problem &problem, crvec μ, crvec y){

    // eval_grad_f_plus_prod_grad_g_y(k, problem, It.grad_ψ, It.Jfxu, It.qr, It.Z);

    auto nx = problem.get_nx();
    auto nu = problem.get_nu();
    auto nxu = problem.get_nx() + problem.get_nu(); 
    auto N = problem.get_N() + 1;

    Eigen::Index k_ = k;

    vec xo(nx);
    problem.get_x_init(xo);

    if (k_ == N-1){
        It.grad_ψ.segment(k*nxu,nx) = It.qr.segment(k*nxu,nx) - It.Z.segment((k-1)*nx,nx);                            
    } 
    else if(k_ == 0){
        It.grad_ψ.segment(nx,nu) = It.qr.segment(nx,nu) +
                    It.Jfxu.block(0,nx,nx,nu).transpose() * It.Z.segment(0,nx);
        It.grad_ψ.segment(0,nx) = It.qr.segment(0,nx) - It.Z.bottomRows(nx) +
                    It.Jfxu.block(0,0,nx,nx).transpose() * It.Z.segment(0,nx);
    }
    else {
        It.grad_ψ.segment(k*nxu,nxu) = It.qr.segment(k*nxu,nxu) - mat::Identity(nxu, nx)*It.Z.segment((k-1)*nx,nx) +
                    It.Jfxu.block(k*nx,0,nx,nxu).transpose() * It.Z.segment(k*nx,nx);
    }
}

template <typename Iterate, typename Problem>
real_t eval_ψ_k_fun (int k, Iterate &It, Problem &problem, crvec μ, Box &F){
    
    auto nx = problem.get_nx();
    auto N = problem.get_N() + 1;

    Eigen::Index k_ = k;
    real_t ψxu_k = 0;

    if (k_ == N-1){
        return ψxu_k = It.lxu(k);
    }
    else if (k_ == 0){
        auto d1 = (It.gz.segment(k*nx,nx)-
                eval_proj_set(F,It.gz.segment(k*nx,nx)));
        auto d2 = (It.gz.bottomRows(nx)-
                eval_proj_set(F,It.gz.bottomRows(nx)));
        return ψxu_k = It.lxu(k) + 
                    0.5*d1.transpose()*(μ).segment(k*nx,nx).asDiagonal()*d1  
                    + 0.5*d2.transpose()*(μ).bottomRows(nx).asDiagonal()*d2;
    }
    else {
        auto d1 = (It.gz.segment(k*nx,nx)-
                eval_proj_set(F,It.gz.segment(k*nx,nx)));
        return ψxu_k = It.lxu(k) + 
                    0.5*d1.transpose()*(μ).segment(k*nx,nx).asDiagonal()*d1;
    }
} 

template <typename Iterate, typename Problem>
real_t eval_ψ_hat_k_fun (int k, Iterate &It, Problem &problem, crvec μ, crvec y, Box &F){
    
    auto nx = problem.get_nx();
    auto nu = problem.get_nu();
    auto nxu = problem.get_nx() + problem.get_nu(); 
    auto N = problem.get_N() + 1;
    auto nh = problem.get_nh();
    auto nh_N = problem.get_nh_N();

    vec xo(nx);
    problem.get_x_init(xo);

    real_t ψxû_k = 0;
    Eigen::Index i_ = k;

    if (i_ == N-1){
        problem.eval_h_N(It.xû.segment(k*nxu,nx),
                        It.hxû.segment(k*nxu,nx));
        It.lxû(k) = problem.eval_l_N(It.hxû.segment(k*nh,nh_N));
        return ψxû_k = It.lxû(k);
    } 
    else if (i_ == 0){
        problem.eval_f(k, It.xû.segment(k*nxu,nx),
                    It.xû.segment(k*(nxu)+nx,nu), 
                    It.fxû.segment(k*nx,nx));
        problem.eval_h(k, It.xû.segment(k*nxu,nx), 
                        It.xû.segment((k*nxu)+nx,nu),
                        It.hxû.segment(k*nxu,nxu));
        It.lxû(k) = problem.eval_l(k, It.hxû.segment(k*nh,nh)); 
        // f(x_{0},u_{0}) = x_{1}                      
        It.gz_hat.segment(k*nx,nx) = It.fxû.segment(k*nx,nx) - 
                        It.xû.segment((k+1)*nxu,nx)  
                        + ((μ).segment(k*nx,nx).asDiagonal().inverse()*(y.segment(k*nx,nx)));
        It.Π_D_hat.segment(k*nx,nx) = eval_proj_set(F, It.gz_hat.segment(k*nx,nx));  
        auto d = (It.gz_hat.segment(k*nx,nx)-
                eval_proj_set(F,It.gz_hat.segment(k*nx,nx))); 
        // x_{0} = x            
        It.gz_hat.bottomRows(nx) = xo - It.xû.segment(k*nxu,nx) + 
                                ((μ).bottomRows(nx).asDiagonal().inverse()*(y.bottomRows(nx)));            
        It.Π_D_hat.bottomRows(nx).noalias() = eval_proj_set(F, It.gz_hat.bottomRows(nx)); 
        auto d1 = (It.gz_hat.bottomRows(nx)-
                eval_proj_set(F,It.gz_hat.bottomRows(nx)));
        //
        return ψxû_k = It.lxû(k) + 
                    0.5*d.transpose()*(μ).segment(k*nx,nx).asDiagonal()*d + 
                    0.5*d1.transpose()*(μ).bottomRows(nx).asDiagonal()*d1; 
    }
    else {
        problem.eval_f(k, It.xû.segment(k*nxu,nx),
                    It.xû.segment(k*(nxu)+nx,nu), 
                    It.fxû.segment(k*nx,nx));
        problem.eval_h(k, It.xû.segment(k*nxu,nx), 
                        It.xû.segment((k*nxu)+nx,nu),
                        It.hxû.segment(k*nxu,nxu));
        It.lxû(k) = problem.eval_l(k, It.hxû.segment(k*nh,nh));                       
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

template <typename Problem>
void eval_prox_stop_crit_fun (int k, Problem &problem, crvec grad_ψ, crvec xu, rvec work_xu, rvec work_p, Box &D_N, Box &D, Box &U){
    
    auto nx = problem.get_nx();
    auto nu = problem.get_nu();
    auto nxu = problem.get_nx() + problem.get_nu(); 
    auto N = problem.get_N() + 1;

    Eigen::Index k_ = k;

    if (k_ == N-1){
        eval_proj_grad_step_box(D_N, real_t(1.), xu.segment(k*nxu,nx), 
                                grad_ψ.segment(k*nxu,nx), 
                                work_xu.segment(k*nxu,nx), work_p.segment(k*nxu,nx));
    } 
    else {
        eval_proj_grad_step_box(D, real_t(1.), xu.segment(k*nxu,nx), 
                                grad_ψ.segment(k*nxu,nx), 
                                work_xu.segment(k*nxu,nx), work_p.segment(k*nxu,nx));
        eval_proj_grad_step_box(U, real_t(1.), xu.segment((k*nxu)+nx,nu), 
                                grad_ψ.segment((k*nxu)+nx,nu), 
                                work_xu.segment((k*nxu)+nx,nu), work_p.segment((k*nxu)+nx,nu));
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
        work.rightCols(nx)*((μ).bottomRows(nx).asDiagonal());
        problem.eval_add_Q(k, It.xu.segment(k*nxu,nxu), 
                            It.hxu.segment(k*nxu,nxu),
                            work_eval_Q);
        It.GN.block(k*nxu,k*nxu,nxu,nxu) += work_eval_Q;
        It.GN.block(k*nxu,k*nxu,nxu+nx,nxu+nx) += (work.transpose() * 
                                        (μ).segment(k*nx,nx).asDiagonal()) *
                                        (work);
    }
    else {
        mat work_eval_Q(nxu,nxu);
        work.leftCols(nxu) = It.Jfxu.block(k*nx,0,nx,nxu);
        work.rightCols(nx).setIdentity();
        work.rightCols(nx)*((μ).segment((k-1)*nx,nx).asDiagonal());
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
void fd_grad_ψ_fun (int k, Iterate &It, Problem &problem, crvec μ, crvec y, Box &F, real_t h){
    alpaqa::ScopedMallocAllower ma;
        
    real_t ψxu_     = 0;
    index_t N       = problem.get_N() + 1; 
    index_t nu      = problem.get_nu();
    index_t nx      = problem.get_nx();
    index_t nh      = problem.get_nh();
    index_t nh_N    = problem.get_nh_N();
    index_t nxu     = nx + nu;
    index_t n       = (N * nxu) - nu;
    index_t m       = N * nx;
    vec lxu_        = vec::Zero(N);
    vec hxu_        = vec::Zero(n);
    vec xu_         = vec::Zero(n); 
    vec fxu_        = vec::Zero(m);
    vec xo          = vec::Zero(nx);
    vec d           = vec::Zero(m);
    vec v           = vec::Zero(m);
    vec g           = vec::Zero(m);
    Box D           = alpaqa::Box<config_t>::NaN(m);
    problem.get_D(D);
    problem.get_x_init(xo);

    for (index_t i = 0; i < n; ++i) {
        xu_ = It.xu;
        xu_(i) += h;
        for (index_t k = 0; k < N; ++k) {
            if (k == N-1) {
                problem.eval_h_N(xu_.segment(k*nxu,nx),
                                    hxu_.segment(k*nxu,nx));
                lxu_(k) = problem.eval_l_N(hxu_.segment(k*nh,nh_N));
            } else {
                problem.eval_f(k, xu_.segment(k*nxu,nx),
                                xu_.segment(k*(nxu)+nx,nu), 
                                fxu_.segment(k*nx,nx));
                problem.eval_h(k, xu_.segment(k*nxu,nx), 
                                xu_.segment((k*nxu)+nx,nu),
                                hxu_.segment(k*nxu,nxu));
                lxu_(k) = problem.eval_l(k, hxu_.segment(k*nh,nh)); 
            }
        }
        for (index_t k = 0; k < N-1; ++k) {
            g.segment(k*nx,nx) = fxu_.segment(k*nx,nx) - 
                                    xu_.segment((k+1)*nxu,nx);
        }
        g.bottomRows(nx) = xu_.segment(0,nx) - xo;
        d.noalias() = g + μ.asDiagonal().inverse()*y; 
        for (index_t i = 0; i < N-1; ++i){
            v.segment(k*nx,nx) = d.segment(k*nx,nx) - eval_proj_set(F, d.segment(k*nx,nx));
        } 
        v.bottomRows(nx) = d.bottomRows(nx) - eval_proj_set(F, d.bottomRows(nx));
        ψxu_ = lxu_.sum() + real_t(.5)*v.transpose()*(μ).asDiagonal()*v;
        It.grad_ψ(i) = (ψxu_ - It.ψxu)/h;
    }
}

