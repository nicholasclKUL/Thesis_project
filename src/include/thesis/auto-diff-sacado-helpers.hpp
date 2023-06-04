#pragma once

#include <Eigen/Core>
#include <Sacado.hpp>
#include <Kokkos_Core.hpp>
#include <alpaqa/config/config.hpp>

USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

// Struct for the AD views
template <const int p_nxu, const int p_nx, const int p_nh>
struct AD{
  
  typedef Sacado::Fad::SFad<real_t,p_nxu> FadType1;   ///< Sacado AD type for nxu derivatives
  typedef Sacado::Fad::SFad<real_t,p_nxu> FadType2;   ///< Sacado AD type for nxu derivatives
  
  FadType2 Qr;

  Kokkos::View<FadType1*> xu_fad,    ///< view of states and inputs
                          fxu_fad,   ///< view of dynamics
                          h_fad;     ///< view of output mapping
  Kokkos::View<FadType2*> l_fad;     ///< view of stages cost function

  const int nx,   ///< total number of states  
            nxu,  ///< total number of states and input
            nh;

  AD() : xu_fad("xu_fad", p_nxu, p_nxu), fxu_fad("fxu_fad", p_nx, p_nxu),
          h_fad("h_fad", p_nh, p_nxu), l_fad("l_fad", 1, p_nh), Qr(p_nh, 0, 0),
          nx(p_nx), nxu(p_nxu), nh(p_nh) {};

};

// Assign x and u values to their respectives AD views
template <const int nxu, const int nx, const int nh>
void assign_values_xu (crvec x, crvec u, const AD<nxu,nx,nh> &ad){
  index_t k = 0;
  for (size_t i = 0; i < x.size()+u.size(); ++i){
    if (i < x.size()){
      ad.xu_fad(i).val() = x(i);
    }
    else {
      ad.xu_fad(i).val() = u(k);
      k++;
    }
    ad.xu_fad(i).diff(i, x.size()+u.size()); //generate vector of duals
  }
}

template <const int nxu, const int nx, const int nh>
void assign_values_xu (crvec xu, const AD<nxu,nx,nh> &ad){
  for (size_t i = 0; i < xu.size(); ++i){
    ad.xu_fad(i).val() = xu(i);
    ad.xu_fad(i).diff(i, xu.size()); //generate vector of duals
  }
}

// template <const int nxu, const int nx, const int nh>
// void assign_values_Qr (crvec QR, const AD<nxu,nx,nh> &ad){
//   for (size_t i = 0; i < QR.size(); ++i){
//     ad.Qr(i).val() = QR(i)
//   }
// }

// Assign the Jacobian AD calculated view to its original container
template <const int nxu, const int nx, const int nh>
void assign_values (rmat Jfxu, const AD<nxu,nx,nh> &ad){
  for (size_t i = 0; i < nx; ++i){
    for (size_t j = 0; j < nxu; ++j){
      Jfxu(i,j) = ad.fxu_fad(i).fastAccessDx(j);
    }
  }
}

// Assign grad-vector product to input container
template <const int nxu, const int nx, const int nh>
void assign_values (rvec grad_fxu_p, crvec p, const AD<nxu,nx,nh> &ad){
  grad_fxu_p.setConstant(0.);
  for (size_t i = 0; i < nxu; ++i){
    for (size_t j = 0; j < nx; ++j){
      grad_fxu_p(i) += ad.fxu_fad(j).fastAccessDx(i)*p(j);
    }
  }
}

// template <const int nxu, const int nx, const int nh>
// void assign_values (rvec qr, const AD<nxu,nx,nh> &ad){
//   for (size_t i = 0; i < nxu; ++i){
//     for (size_t j = 0; j < nx; ++j){
//       qr(i) += ad.h_fad(j).fastAccessDx(i)*ad.l_fad(j).val();
//     }
//   }
// }