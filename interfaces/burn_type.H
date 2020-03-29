#ifndef _burn_type_H_
#define _burn_type_H_

#include <network.H>

// A generic structure holding data necessary to do a nuclear burn.

// Set the number of independent variables -- this should be
// temperature, enuc + the number of species which participate
// in the evolution equations.

const int neqs = 2 + NumSpec;

// Indices of the temperature and energy variables in the work arrays.

const int net_itemp = NumSpec + 1;
const int net_ienuc = NumSpec + 2;

struct burn_t
{
    Real rho;
    Real T;
    Real e;
    Real xn[NumSpec];
#if NumAux > 0
    Real aux[NumAux];
#endif

    Real cv;
    Real cp;
    Real y_e;
    Real eta;
    Real cs;
    Real dx;
    Real abar;
    Real zbar;

    // Last temperature we evaluated the EOS at
    Real T_old;
    Real cv_old;
    Real cp_old;

    // Temperature derivatives of specific heat
    Real dcvdT;
    Real dcpdT;

    // Zone index information.
    int i, j, k;

    // diagnostics
    int n_rhs, n_jac;

    // Was the burn successful?
    bool success;
};

// Given an eos type, copy the data relevant to the burn type.

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eos_to_burn (const eos_t& eos_state, burn_t& burn_state)
{

    burn_state.rho  = eos_state.rho;
    burn_state.T    = eos_state.T;
    burn_state.e    = eos_state.e;
    for (int n = 0; n < NumSpec; ++n) {
        burn_state.xn[n] = eos_state.xn[n];
    }
#if NumAux > 0
    for (int n = 0; n < NumAux; ++n) {
        burn_state.aux[n] = eos_state.aux[n];
    }
#endif
    burn_state.cv   = eos_state.cv;
    burn_state.cp   = eos_state.cp;
    burn_state.y_e  = eos_state.y_e;
    burn_state.eta  = eos_state.eta;
    burn_state.cs   = eos_state.cs;
    burn_state.abar = eos_state.abar;
    burn_state.zbar = eos_state.zbar;

}



// Given a burn type, copy the data relevant to the eos type.

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void burn_to_eos (const burn_t& burn_state, eos_t& eos_state)
{

    eos_state.rho  = burn_state.rho;
    eos_state.T    = burn_state.T;
    eos_state.e    = burn_state.e;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = burn_state.xn[n];
    }
#if NumAux > 0
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.aux[n] = burn_state.aux[n];
    }
#endif
    eos_state.cv   = burn_state.cv;
    eos_state.cp   = burn_state.cp;
    eos_state.y_e  = burn_state.y_e;
    eos_state.eta  = burn_state.eta;
    eos_state.cs   = burn_state.cs;
    eos_state.abar = burn_state.abar;
    eos_state.zbar = burn_state.zbar;

}



AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void normalize_abundances_burn (burn_t& state)
{

    Real sum = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        state.xn[n] = amrex::max(small_x, amrex::min(1.0_rt, state.xn[n]));
        sum += state.xn[n];
    }
    for (int n = 0; n < NumSpec; ++n) {
        state.xn[n] /= sum;
    }

}

#endif
