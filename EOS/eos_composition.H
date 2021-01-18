
#ifndef _eos_composition_H_
#define _eos_composition_H_

#include <AMReX_BLFort.H>
#include <network.H>
#include <eos_type.H>

using namespace amrex;

#ifdef NSE_THERMO
using namespace AuxZero;
#endif

struct eos_xderivs_t {
  Real dedX[NumSpec];
  Real dpdX[NumSpec];
  Real dhdX[NumSpec];
};

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void composition (T& state)
{
  // Calculate abar, the mean nucleon number,
  // zbar, the mean proton number,
  // mu, the mean molecular weight,
  // mu_e, the mean number of nucleons per electron, and
  // y_e, the electron fraction.

#ifdef NSE_THERMO

  state.mu_e = 1.0 / state.aux[iye];
  state.y_e = state.aux[iye];
  state.abar = state.aux[iabar];
  state.zbar = state.abar * state.y_e;

#else

  Real sum = 0;
  for (int n = 0; n < NumSpec; n++) {
    sum += state.xn[n] * zion[n] * aion_inv[n];
  }

  state.mu_e = 1.0 / sum;
  state.y_e = sum;

  sum = 0;

  for (int n = 0; n < NumSpec; n++) {
    sum += state.xn[n] * aion_inv[n];
  }

  state.abar = 1.0 / sum;
  state.zbar = state.abar * state.y_e;
#endif
}

#ifdef EXTRA_THERMO
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
eos_xderivs_t composition_derivatives (const T& state)
{
  eos_xderivs_t state_xderivs;

  for (int n = 0; n < NumSpec; n++) {
    state_xderivs.dpdX[n] =
      state.dpdA * (state.abar * aion_inv[n]) * (aion[n] - state.abar) +
      state.dpdZ * (state.abar * aion_inv[n]) * (zion[n] - state.zbar);

    state_xderivs.dedX[n] =
      state.dedA * (state.abar * aion_inv[n]) * (aion[n] - state.abar) +
      state.dedZ * (state.abar * aion_inv[n]) * (zion[n] - state.zbar);

    if (state.dpdr != 0.0) {

      state_xderivs.dhdX[n] = state_xderivs.dedX[n]
            + (state.p / (state.rho*state.rho) - state.dedr)
            * state_xderivs.dpdX[n] / state.dpdr;
    }
  }

  return state_xderivs;
}
#endif


#endif
