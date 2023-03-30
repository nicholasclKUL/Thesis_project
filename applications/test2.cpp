#include <alpaqa/config/config.hpp>
#include <Thesis/ocp-funcs.hpp>
#include <Kokkos_Core.hpp>
#include <iostream>

#include <alpaqa/config/config.hpp>
#include <alpaqa/util/float.hpp>
#include <alpaqa/problem/box.hpp>

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

using Box = alpaqa::Box<config_t>;
using Problem = alpaqa::TypeErasedControlProblem<config_t>;

struct Iterate {

    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);  
    
    vec xu;             //< Inputs u interleaved with states x -> x₀, u₀, x₁, u₁,...,uₙ₋₁, xₙ 
    vec xû;             //< Inputs u interleaved with states x after prox grad
    vec grad_ψ;         //< Gradient of cost w.r.t {x₀, u₀, x₁, u₁,...,uₙ₋₁, xₙ}
    vec z;              //< f(xᵏ,uᵏ)-xᵏ⁺¹
    vec p;              //< Proximal gradient step w.r.t {x₀, u₀, x₁, u₁,...,uₙ₋₁, xₙ}
    vec u;              //< Inputs u (used for L-BFGS only) Do I need it???
    vec fxu;            //< Dynamics f(x,u)
    vec hxu;            //< nonconvex mapping h(x,u)
    mat Jfxu;           //< Jacobian of dynamics Jf(x,u)
    vec qr;             //< Cost function gradient Jh(x,u)dl(x,u)
    vec Π_D;            //< Projection of states and variables into their box constraints
    vec Z;              //< (Σ*(f(xₖ,uₖ)-xₖ₊₁-Π_D(f(xₖ,uₖ)-xₖ₊₁+Σ⁻¹y))-y)
    vec l;              //< cost function l(x,u) in each stage in the horizon N
    real_t ψxu = 0;

    Iterate(length_t n, length_t m, length_t nxu, length_t N) :
    xu(n), xû(n), grad_ψ(n), z(n), p(n), fxu(m), hxu(n), Jfxu(m,nxu),
    qr(n), Π_D(n), Z(n), l(N) {}
};


auto eval_proj_set = [](const Box BoxSet, crvec x) {
    using binary_real_f = real_t (*)(real_t, real_t);
    return x.binaryExpr(BoxSet.lowerbound, binary_real_f(std::fmax))
            .binaryExpr(BoxSet.upperbound, binary_real_f(std::fmin));
};


auto eval_iterate = [](int k, Problem &problem, Iterate &It, crvec μ, crvec y){       

    auto N          = problem.get_N();
    auto nx         = problem.get_nx();
    auto nu         = problem.get_nu();
    auto nxu        = nx+nu;
    Eigen::Index i_ = k;
    Box D           = alpaqa::Box<config_t>::NaN(nx);

    problem.get_D(D);

    if (i_ == N-1){
        
        problem.eval_h_N(It.xu.segment(k*nxu,nx),
                         It.hxu.segment(k*nxu,nx));

        It.l(k) = problem.eval_l_N(It.hxu.segment(k*nxu,nx));

        problem.eval_q_N(It.xu.segment(k*nxu,nx),
                         It.hxu.segment(k*nxu,nx), 
                         It.qr.segment(k*nxu,nx));

    } else {

        problem.eval_f(k, It.xu.segment(k*nxu,nx),
                          It.xu.segment(k*(nxu)+nx,nu), 
                          It.fxu.segment(k*nx,nx));

        problem.eval_jac_f(k, It.xu.segment(k*nxu,nx),
                          It.xu.segment(k*(nxu)+nx,nu),
                          It.Jfxu.block(k*nx,0,nx,nxu));
        
        problem.eval_h(k, It.xu.segment(k*nxu,nx), 
                          It.xu.segment((k*nxu)+nx,nu),
                          It.hxu.segment(k*nxu,nxu));
        
        It.l(k) = problem.eval_l(k, It.hxu.segment(k*nxu,nxu));

        problem.eval_qr(k, It.xu.segment(k*nxu,nxu), 
                           It.hxu.segment(k*nxu,nxu),
                           It.qr.segment(k*nxu,nxu));
        
        It.Π_D.segment(k*nx,nx) = It.fxu.segment(k*nx,nx) - 
                                It.xu.segment((k+1)*nxu,nx)  
                                + ((μ).segment(k*nx,nx).asDiagonal().inverse()*(y.segment(k*nx,nx)));
        
        It.Π_D.segment(k*nx,nx) = eval_proj_set(D, It.Π_D.segment(k*nx,nx));

        It.Z.segment(k*nx,nx)   = (μ.segment(k*nx,nx).cwiseProduct(It.fxu.segment(k*nx,nx)
                                -It.xu.segment((k+1)*nxu,nx)-It.Π_D.segment(k*nx,nx)) 
                                + y.segment(k*nx,nx));   
    
    }
};

auto eval_grad_ψ_k = [](int k, Problem &problem, Iterate &It, crvec μ, crvec y){

    auto N          = problem.get_N();
    auto nx         = problem.get_nx();
    auto nu         = problem.get_nu();
    auto nxu        = nx+nu;
    Eigen::Index k_ = k;

    mat I; I.setIdentity(nxu,nx);

    if (k_ == N-1) {
        It.grad_ψ.segment(k*nxu,nx) = It.qr.segment(k*nxu,nx) - It.Z.segment((k-1)*nx,nx);            
    } 
    else if (k_ == 0) {
        It.grad_ψ.segment(k*nxu,nxu) = It.qr.segment(k*nxu,nxu) +
                    It.Jfxu.block(k*nx,0,nx,nxu).transpose() * It.Z.segment(k*nx,nx);
    }
    else {
        It.grad_ψ.segment(k*nxu,nxu) = It.qr.segment(k*nxu,nxu) - I*It.Z.segment((k-1)*nx,nx) +
                    It.Jfxu.block(k*nx,0,nx,nxu).transpose() * It.Z.segment(k*nx,nx);                       
    }
};

auto eval_ψ_k = [](int k, Problem &problem, Iterate &It){

    auto N          = problem.get_N();
    auto nx         = problem.get_nx();
    Eigen::Index k_ = k;
    real_t ψ_k;

    if (k_ == N-1){
        return ψ_k = It.l(k);
    }
    else{
        return ψ_k = It.l(k) + (It.Π_D.segment(k*nx,nx).squaredNorm())/2;
    }
    
};

auto FD = [](Problem &problem, real_t &h, Iterate &it, crvec μ, crvec y){
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
        xu_ = it.xu;
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
        d = g + μ.asDiagonal().inverse()*y; 
        v = d - eval_proj_set(D, d); 
        ψxu_ = lxu_.sum() + real_t(.5)*v.transpose()*(μ).asDiagonal()*v;
        it.grad_ψ(i) = (ψxu_ - it.ψxu)/h;
    }
};

int main(){

    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

    auto problem = alpaqa::TypeErasedControlProblem<config_t>::make<alpaqa::OCProblem>();

    auto nx     = problem.get_nx();
    auto nu     = problem.get_nu();
    auto N      = problem.get_N();
    auto nxu    = nx + nu;
    auto n      = (nxu*(N-1))+nx;
    auto m      = nx*(N-1);

    real_t ψ    = 0;
    real_t ψ2   = 0;
    vec grad_ψ  = vec::Zero(n);

    Iterate ite(n, m, nxu, N);
    Iterate ite2(n, m, nxu, N);
    problem.get_x_init(ite.xu);
    problem.get_x_init(ite2.xu);

    srand(10);

    vec μ = vec::Ones(ite.fxu.size());
    vec y = vec::Zero(ite.fxu.size());

    // -------------------- Parallelization testing --------------------------- //

    Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(4));

    Kokkos::parallel_for(N, [&] (const int i){
            eval_iterate(i, problem, ite, μ, y);
            eval_iterate(i, problem, ite2, μ, y); 
        }
    );

    Kokkos::fence();

    Kokkos::parallel_for(N, [&] (const int i){ 
            eval_grad_ψ_k(i, problem, ite, μ, y);
        }
    );

    Kokkos::fence();

    Kokkos::parallel_reduce(N, [&] (const int i, real_t &ψ_){
            ψ_ += eval_ψ_k(i, problem, ite);
        },
    ψ);
    Kokkos::parallel_reduce(N, [&] (const int i, real_t &ψ_){
            ψ_ += eval_ψ_k(i, problem, ite2);
        },
    ψ2);
    ite2.ψxu = ψ2;

    Kokkos::fence();

    Kokkos::finalize();

    // ----------------------------- FD testing --------------------------------- //

    real_t h = 0.01;

    FD(problem, h, ite2, μ, y);

    // ------------------------------ Printing ---------------------------------- //

    Eigen::IOFormat CleanFmt(4, 0, ", ", "      \n", "[", "]");

    for (index_t i = 0; i < N; ++i){
        if (i==N-1){
            std::cout << "* Stage : " << i+1 << '\n'
                << "(x,u)   = " << ite.xu.segment(i*nxu,nx).transpose() << '\n'
                << "    Z   = " << ite.Z.segment((i-1)*nx,nx).transpose() << '\n'
                << "    qr  = " << ite.qr.segment(i*nxu,nx).transpose() << '\n'
                << "∇ψ(x,u) = " << ite.grad_ψ.segment(i*nxu,nx).transpose() << '\n'
                << "    l   = " << ite.l(i) << '\n'
                << std::endl;
        }
        else if (i==0){
            std::cout << "* Stage : " << i+1 << '\n'
                << "f(x,u)  = " << ite.fxu.segment(i*nx,nx).transpose() << '\n'
                << "(x,u)   = " << ite.xu.segment(i*nxu,nxu).transpose() << '\n'
                << "    Z   = " << ite.Z.segment(i*nx,nx).transpose() << '\n'
                << "    qr  = " << ite.qr.segment(i*nxu,nxu).transpose() << '\n'
                << "Jf(x,u) = " << ite.Jfxu.block(i*nx,0,nx,nxu).transpose().format(CleanFmt) << '\n'
                << "∇ψ(x,u) = " << ite.grad_ψ.segment(i*nxu,nxu).transpose() << '\n'
                << "    l   = " << ite.l(i) << '\n'
                << "Π_D     = " << ite.Π_D.segment(i*nx,nx).squaredNorm() << '\n'
                << std::endl;
        }
        else{
            std::cout << "* Stage : " << i+1 << '\n'
                << "f(x,u)  = " << ite.fxu.segment(i*nx,nx).transpose() << '\n'
                << "(x,u)   = " << ite.xu.segment(i*nxu,nxu).transpose() << '\n'
                << "    Z   = " << ite.Z.segment(i*nx,nx).transpose() << '\n'
                << "Jf(x,u) = " << ite.Jfxu.block(i*nx,0,nx,nxu).transpose().format(CleanFmt) << '\n'
                << "    qr  = " << ite.qr.segment(i*nxu,nxu).transpose() << '\n'
                << "∇ψ(x,u) = " << ite.grad_ψ.segment(i*nxu,nxu).transpose() << '\n'
                << "    l   = " << ite.l(i) << '\n'
                << "Π_D     = " << ite.Π_D.segment(i*nx,nx).squaredNorm() << '\n'
                << std::endl;
        }
    }
    std::cout<< "ψ(x,u)         = " << ψ <<std::endl;
    std::cout<< "∇ψ(x,u)        = " << ite.grad_ψ.transpose() <<std::endl;
    std::cout<< "∇ψ(x,u) (FD)   = " << ite2.grad_ψ.transpose() <<std::endl;  
}