#pragma once

#include <Eigen/Core>
#include <Sacado.hpp>
#include <Kokkos_Core.hpp>
#include <alpaqa/config/config.hpp>
#include <thesis/auto-diff-sacado-helpers.hpp>
#include <functional>

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

template <typename X, typename F, typename Params>
void dynamics(index_t k, X &x, F &f, Params &params);

// Forward (Explicit) Euler

template <typename X, typename F, typename Params>
void fe(index_t &k, X &xu, F &fxu, const Params &params){
    dynamics(k, xu, fxu, params);
    for (index_t i = 0; i < fxu.size(); ++i){
      fxu(i) = xu(i) + params.Ts * fxu(i);
    }
}

// Explicit 4th-order Runge-Kutta 

template <typename X, typename K, typename Temp>
void assign_rk4_temp(X &xu, K &k, Temp &temp, real_t c){
    for (index_t i = 0; i < xu.size(); ++i){
      if (i < k.size()){
        temp(i) = xu(i) + c * k(i);
      }
      else {
        temp(i) = xu(i);
      }
    }    
}
template <typename X, typename F, typename Params>
void rk4(index_t &k, X &xu, F &fxu, const Params &params){

    // declare the intermediate containers
    X k1(fxu), k2(fxu), k3(fxu), k4(fxu), temp(xu);

    // k1
    dynamics(k, xu, k1, params);
    assign_rk4_temp(xu, k1, temp, real_t(params.Ts/2.));
    // k2
    dynamics(k, temp, k2, params);
    assign_rk4_temp(xu, k2, temp, real_t(params.Ts/2.));
    // k3
    dynamics(k, temp, k3, params);
    assign_rk4_temp(xu, k3, temp, real_t(params.Ts));
    // k4
    dynamics(k, temp, k4, params);

    // dx
    for (index_t i = 0; i < fxu.size(); ++i) {
        fxu(i) = xu(i) + params.Ts / 6.0 * (k1(i) + 2.0 * k2(i) + 2.0 * k3(i) + k4(i));
    }
}
template <typename X, typename F, typename Params>
void rk4(index_t &k, X &xu, F &fxu, const Params &params, const int p_nxu, const int p_nx){

    // declare the intermediate containers
    X k1("k1", p_nx, p_nxu), k2("k2", p_nx, p_nxu), k3("k3", p_nx, p_nxu), 
        k4("k4", p_nx, p_nxu), temp("temp", p_nxu, p_nxu);

    // k1
    dynamics(k, xu, k1, params);
    assign_rk4_temp(xu, k1, temp, real_t(params.Ts/2.));
    // k2
    dynamics(k, temp, k2, params);
    assign_rk4_temp(xu, k2, temp, real_t(params.Ts/2.));
    // k3
    dynamics(k, temp, k3, params);
    assign_rk4_temp(xu, k3, temp, real_t(params.Ts));
    // k4
    dynamics(k, temp, k4, params);

    // dx
    for (index_t i = 0; i < fxu.size(); ++i) {
        fxu(i) = xu(i) + (params.Ts / 6.0) * (k1(i) + 2.0 * k2(i) + 2.0 * k3(i) + k4(i));
    }
}