#ifndef VODE_TYPE_STRANG_H
#define VODE_TYPE_STRANG_H

#include <eos.H>

template <typename I>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void integrator_to_burn (const I& int_state, burn_t& state)
{
    // Copy the integration data to the burn state.

    for (int n = 1; n <= NumSpec; ++n) {
        state.xn[n-1] = int_state.y(n);
    }

    state.e = int_state.y(net_ienuc);

    state.time = int_state.t;
}


template <typename I>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void burn_to_integrator (const burn_t& state, I& int_state)
{
    // Copy the integration data from the burn state.

    for (int n = 1; n <= NumSpec; ++n) {
        int_state.y(n) = state.xn[n-1];
    }

    int_state.y(net_ienuc) = state.e;
}


template <typename I>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void renormalize_species (I& int_state)
{
    Real sum = 0.0_rt;

    for (int n = 1; n <= NumSpec; ++n) {
        sum += int_state.y(n);
    }

    for (int n = 1; n <= NumSpec; ++n) {
        int_state.y(n) /= sum;
    }
}


template <typename I>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void clean_state (I& int_state)
{

    // Ensure that mass fractions always stay positive and less than or
    // equal to 1.

    if (do_species_clip) {
        for (int n = 1; n <= NumSpec; ++n) {
            int_state.y(n) = amrex::max(amrex::min(int_state.y(n), 1.0_rt), SMALL_X_SAFE);
        }
    }

    // for primordial chem, just enforce positivity
    if (use_number_densities) {
        for (int n = 1; n <= NumSpec; ++n) {
            int_state.y(n) = amrex::max(int_state.y(n), SMALL_X_SAFE);
        }
    }

    // Renormalize the abundances as necessary.

    if (renormalize_abundances) {
        renormalize_species(int_state);
    }

}

///
/// update the thermodynamics in the burn_t state -- this may involve an EOS call.
/// we also pass in the int_state so we get the latest values of the mass fractions,
/// temperature, and internal energy
///
template <typename I>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void update_thermodynamics (burn_t& state, const I& int_state)
{
#ifdef AUX_THERMO
    // we are handling the thermodynamics via the aux quantities, which
    // are stored in the burn_t here, so we need to update those based on
    // the current state.

    state.abar = 0.0_rt;
    state.y_e = 0.0_rt;
    for (int n = 1; n <= NumSpec; ++n) {
      state.abar += int_state.y(n) * aion_inv[n-1];
      state.y_e += int_state.y(n) * zion[n-1] * aion_inv[n-1];
    }
    state.abar = 1.0_rt / state.abar;
    state.zbar = state.abar * state.y_e;

    state.aux[AuxZero::iabar] = state.abar;
    state.aux[AuxZero::iye] = state.y_e;
    state.aux[AuxZero::ibea] = 0.0_rt; // this is not needed when we integrate

#endif

    for (int n = 1; n <= NumSpec; ++n) {
        state.xn[n-1] = int_state.y(n);
    }
    state.e = int_state.y(net_ienuc);

    // Get T from e (also updates composition quantities).

    if (call_eos_in_rhs) {
        eos(eos_input_re, state);
    }

    // override T if we are fixing it (e.g. due to
    // drive_initial_convection)
    if (state.T_fixed > 0.0_rt) {
        state.T = state.T_fixed;
    }

}

#endif
