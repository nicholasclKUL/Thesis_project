#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/print.hpp>

#include <Kokkos_Core.hpp>

#include <thesis/para-panoc.hpp>
#include <thesis/printing.hpp>
// #include <nonlinear_dynamics.hpp>
#include <linear_dynamics.hpp>
// #include <quadcopter.hpp>
// #include <hanging_chain.hpp>
#include <thesis/ocp-kkt-error.hpp>

#include <iomanip>
#include <iostream>

int main() {

    Kokkos::initialize(Kokkos::InitializationSettings());

    {

        USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

        auto problemAD = alpaqa::TypeErasedControlProblem<config_t>::make<LinearOCPAD>();

        LinearOCPAD locp;

        auto h = locp.params.Ts;
        auto A = locp.params.A;
        auto B = locp.params.B;
        auto nx = locp.params.nx;
        auto nu = locp.params.nu;

        vec xu(nx+nu),
                fxu(nx);

        mat Jfxu(nx,nx+nu),
                JfxuAD(nx,nx+nu);

        xu << 0.7, 1.0, 0.3, 0.5, 0.9;

        // Simulation using RK4

        problemAD.eval_f(0,xu.segment(0,nx),xu.segment(nx,nu),fxu);
        
        auto k1 = A * xu.segment(0,nx) + 
                B* xu.segment(nx,nu);

        auto k2 = A * (xu.segment(0,nx) + (h/2) * k1) + 
                B* xu.segment(nx,nu);

        auto k3 = A * (xu.segment(0,nx) + (h/2) * k2) + 
                B* xu.segment(nx,nu); 

        auto k4 = A * (xu.segment(0,nx) + (h) * k3) + 
                B* xu.segment(nx,nu);

        auto x1 = xu.segment(0,nx) + (h/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);

        std::cout<<x1-fxu<<std::endl;

        // Jacobian of RK4

        auto dk1_dx = h*A;
        auto dk1_du = B;

        auto dk2_dx = (std::pow(h,2)/2)*A*A + h*A;
        auto dk2_du = (std::pow(h,1)/2)*A*B + B;

        auto dk3_dx = (std::pow(h,3)/4)*A*A*A + (std::pow(h,2)/2)*A*A + h*A;
        auto dk3_du = (std::pow(h,2)/4)*A*A*B + (h/2)*A*B + B;

        auto dk4_dx = (std::pow(h,4)/4)*A*A*A*A + (std::pow(h,3)/2)*A*A*A + (std::pow(h,2))*A*A + h*A;
        auto dk4_du = (std::pow(h,3)/4)*A*A*A*B + (std::pow(h,2)/2)*A*A*B + h*A*B + B;

        Jfxu.leftCols(nx) = mat::Identity(nx,nx) + 
                (1./6.) * (dk1_dx + 2*dk2_dx + 2*dk3_dx + dk4_dx);
        Jfxu.rightCols(nu) = (h/6.) * (dk1_du + 2*dk2_du + 2*dk3_du + dk4_du);

        problemAD.eval_jac_f(0, xu.segment(0,nx), xu.segment(nx,nu), JfxuAD);

        std::cout<<Jfxu<<'\n'<<JfxuAD<<std::endl;

        std::cout<<Jfxu-JfxuAD<<std::endl;

    }

    Kokkos::finalize();
}

    