#ifndef _screen_H_
#define _screen_H_

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include <network_properties.H>

#include <cmath>

using namespace amrex;

struct plasma_state_t {

  Real qlam0z;
  Real qlam0zdt;
  //Real qlam0zdd;

  Real taufac;
  Real taufacdt;

  Real aa;
  Real daadt;
  //Real daadd;
};


namespace scrn {
  const Real fact = 1.25992104989487e0_rt;
  const Real co2        = (1.0_rt/3.0_rt) * 4.248719e3_rt;
  const Real gamefx     = 0.3e0_rt;
  const Real gamefs     = 0.8e0_rt;
  const Real h12_max    = 300.e0_rt;
};


AMREX_FORCE_INLINE
void
fill_plasma_state(plasma_state_t state, const Real temp, const Real dens, Real* y) {

  Real sum = 0.0_rt;
  for (int n = 0; n < NumSpec; n++) {
    sum += y[n];
  }
  Real abar = 1.0_rt / sum;
  Real ytot = sum;

  sum = 0.0_rt;
  Real sum2 = 0.0_rt;
  for (int n = 0; n < NumSpec; n++) {
    sum += zion[n]*y[n];
    sum2 += zion[n]*zion[n]*y[n];
  }

  Real zbar = sum * abar;
  Real z2bar = sum2 * abar;

  Real rr = dens * ytot;
  Real tempi = 1.0_rt / temp;
  Real dtempi = -tempi * tempi;
  Real deni = 1.0_rt / dens;

  Real pp = std::sqrt(rr*tempi*(z2bar + zbar));
  Real qq = 0.5_rt/pp *(z2bar + zbar);
  Real dppdt = qq*rr*dtempi;
  //Real dppdd = qq * ytot * tempi;

  state.qlam0z = 1.88e8_rt * tempi * pp;
  state.qlam0zdt = 1.88e8_rt * (dtempi*pp + tempi*dppdt);
  //state.qlam0zdd = 1.88e8_rt * tempi * dppdd;

  state.taufac = scrn::co2 * std::pow(tempi, 1.0_rt/3.0_rt);
  state.taufacdt = -(1.0_rt/3.0_rt) * state.taufac * tempi;

  qq = rr * zbar;
  Real xni = std::pow(qq, 1.0_rt/3.0_rt);
  //dxnidd = (1.0_rt/3.0_rt) * xni * deni;

  state.aa = 2.27493e5_rt * tempi * xni;
  state.daadt = 2.27493e5_rt * dtempi * xni;
  //state.daadd = 2.27493e5_rt * tempi * dxnidd;
}

#endif
