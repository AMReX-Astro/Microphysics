
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

// Given a set of mass fractions, calculate quantities that depend
// on the composition like abar and zbar.

AMREX_GPU_HOST_DEVICE inline
void subroutine (eos_t& state)
{

    // Calculate abar, the mean nucleon number,
    // zbar, the mean proton number,
    // mu, the mean molecular weight,
    // mu_e, the mean number of nucleons per electron, and
    // y_e, the electron fraction.

    state.y_e = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        state.y_e = state.y_e + state.xn[n] * zion[n] * aion_inv[n];
    }
    state.mu_e = 1.0_rt / state.y_e;

    Real sum = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        sum = sum + state.xn[n] * aion_inv[n];
    }
    state.abar = 1.0_rt / sum;
    state.zbar = state.abar / state.mu_e;

}



// Compute thermodynamic derivatives with respect to xn(:)

AMREX_GPU_HOST_DEVICE inline
void composition_derivatives (const eos_t& state, eos_derive_t& state_xderivs)
{

    // Get the mass of a nucleon from Avogadro's number.
    const Real m_nucleon = 1.0_rt / n_A;

    // Composition derivatives

#ifdef EXTRA_THERMO
    for (int n = 0; n < NumSpec; ++n) {
        state_xderivs.dpdX[n] = state.rho * (k_B / m_nucleon) * state.T * aion_inv[n];
        state_xderivs.dedX[n] = (k_B / m_nucleon) * state.T * aion_inv[n] / (gammas[n] - 1.0_rt);
        state_xderivs.dhdX[n] = state_xderivs.dedX[n] + (state.p / (state.rho * state.rho) - state.dedr) *
                                                        state_xderivs.dpdX[n] / state.dpdr;
    }
#endif

}

#endif
