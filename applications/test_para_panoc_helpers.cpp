#include <thesis/para-panoc-helpers.hpp>

struct Iterate {

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
    vec qr;         //< Cost function gradient Jh(x,u)dl(x,u)
    vec g;          //< Dynamic constraints, f(x,u) - x+
    vec gz;         //< g(x) + Σ⁻¹y
    vec gz_hat;     //< g(x_hat) + Σ⁻¹y
    vec Π_D;        //< Projection of states and variables into their box constraints
    vec Π_D_hat;    //< Projection of states and variables into their box constraints
    vec Z;          //< (Σ*(f(xₖ,uₖ)-xₖ₊₁-Π_D(f(xₖ,uₖ)-xₖ₊₁+Σ⁻¹y))-y)
    vec p;          //< Proximal gradient step
    vec v_GN;
    vec i_GN;
    vec j_GN;

    spmat GN;         //< GN approximation of ∇²ψ
    
    real_t ψxu      = alpaqa::NaN<config_t>;        //< Cost in x
    real_t ψxû      = alpaqa::NaN<config_t>;        //< Cost in x̂
    real_t γ        = alpaqa::NaN<config_t>;        //< Step size γ
    real_t L        = alpaqa::NaN<config_t>;        //< Lipschitz estimate L
    real_t pᵀp      = alpaqa::NaN<config_t>;        //< Norm squared of p
    real_t grad_ψᵀp = alpaqa::NaN<config_t>;        //< Dot product of gradient and p

    // @pre    @ref ψxu, @ref pᵀp, @pre grad_ψᵀp
    // @return φγ
    real_t fbe() const { return ψxu + pᵀp / (2 * γ) + grad_ψᵀp; }

    Iterate(length_t n, length_t m, length_t nx, length_t nu, length_t N) :
        xu(n), xû(n), grad_ψ(n), fxu(m-nx), fxû(m-nx), hxu(n), hxû(n), 
        qr(n), g(m), gz(m), gz_hat(m), Π_D(m), Π_D_hat(m), Z(m), p(n), 
        Jfxu(m-nx, nx+nu), GN(n,n), lxu(N), lxû(N), i_GN(n*n), j_GN(n*n), v_GN(n*n)
        {
            v_GN.setZero();
            i_GN.setZero();
            j_GN.setZero(); }

};

// void (spmat sp_A, rmat A){
//     for (size_t i = 0; i < A.rows(); ++i){
//         for (size_t j = 0; i < A.cols(); ++j){
//             sp_A.coeffRef(i,j) = A(i,j);
//         }
//     }
// }

int main(){

    Iterate it(4,4,4,4,4);
    mat A(4,4);
    spmat sp_A(4,4);

    A << 1, 0, 0, 0, 
        5, 0, 1, 0,
        0, 3, 0, 0,
        0, 0, 0, 2;

    fill_vecs_GN(it.v_GN,it.i_GN,it.j_GN,A,0,0);

    std::cout<<it.v_GN.transpose()<<'\n'
            <<it.i_GN.transpose()<<'\n'
            <<it.j_GN.transpose()<<'\n'<<std::endl;

    fill_sp_GN(it);

    std::cout<<it.GN<<std::endl;

    for (size_t i = 0; i < A.rows(); ++i){
        for (size_t j = 0; j < A.cols(); ++j){
            if (A(i,j) != 0){
                sp_A.coeffRef(i,j) = A(i,j);
            }
        }
    }

    std::cout<<sp_A<<std::endl;

}