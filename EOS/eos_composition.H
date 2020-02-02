#include <actual_network.H>
#include <eos.H>

#ifndef _eos_composition_H_
#define _eos_composition_H_

#include <AMReX_BLFort.H>

using namespace amrex;

struct eos_xderivs_t {
  Real dedX[NumSpec];
  Real dpdX[NumSpec];
  Real dhdX[NumSpec];
};


inline void composition(eos_t& state) {

  // Calculate abar, the mean nucleon number,
  // zbar, the mean proton number,
  // mu, the mean molecular weight,
  // mu_e, the mean number of nucleons per electron, and
  // y_e, the electron fraction.

  Real sum = 0;
  for (int n = 0; n < NumSpec; n++) {
    sum += state.xn[n] * zion[n] / aion[n];
  }

  state.mu_e = 1.0 / sum;
  state.y_e = 1.0 / state.mu_e;

  sum = 0;

  for (int n = 0; n < NumSpec; n++) {
    sum += state.xn[n] / aion[n];
  }

  state.abar = 1.0 / sum;
  state.zbar = state.abar / state.mu_e;

}

eos_xderivs_t compsition_derivative(eos_t& state);

#endif
