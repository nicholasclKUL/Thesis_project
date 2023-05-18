#pragma once

#include <Eigen/Core>
#include <Sacado.hpp>
#include <Kokkos_Core.hpp>

// Struct for the AD views
template <const int p_nxu, const int p_nx>
struct AD{
  
  typedef Sacado::Fad::SFad<real_t,p_nxu> FadType;
  
  Kokkos::View<FadType*> xu_fad, fxu_fad;
  
  const int nx, nxu;

  AD() : xu_fad("xu_fad", p_nxu, p_nxu), fxu_fad("fxu_fad", p_nx, p_nxu),
          nx(p_nx), nxu(p_nxu) {};

};

// Function to assign x and u values to their respectives AD views
template <const int nxu, const int nx>
void assign_values_xu (crvec x, crvec u, const AD<nxu,nx> &ad){
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