
#ifndef _eos_composition_H_
#define _eos_composition_H_

#include <AMReX_BLFort.H>
#include <network.H>
#include <eos_type.H>

using namespace amrex;

struct eos_xderivs_t {
  Real dedX[NumSpec];
  Real dpdX[NumSpec];
  Real dhdX[NumSpec];
};

AMREX_GPU_HOST_DEVICE inline
void composition(eos_t& state) {

  // Calculate abar, the mean nucleon number,
  // zbar, the mean proton number,
  // mu, the mean molecular weight,
  // mu_e, the mean number of nucleons per electron, and
  // y_e, the electron fraction.

  Real sum = 0;
  for (int n = 0; n < NumSpec; n++) {
    sum += state.xn[n] * zion[n] * aion_inv[n];
  }

  state.mu_e = 1.0 / sum;
  state.y_e = 1.0 / state.mu_e;

  sum = 0;

  for (int n = 0; n < NumSpec; n++) {
    sum += state.xn[n] * aion_inv[n];
  }

  state.abar = 1.0 / sum;
  state.zbar = state.abar / state.mu_e;

}

AMREX_GPU_HOST_DEVICE inline
eos_xderivs_t composition_derivatives(const eos_t& state) {

  eos_xderivs_t state_xderivs;

#ifdef EXTRA_THERMO
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

#endif

  return state_xderivs;
}

#endif
