#ifndef VODE_TYPE_STRANG_H
#define VODE_TYPE_STRANG_H


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void vode_to_burn (const dvode_t& vode_state, burn_t& state)
{
    // Copy the integration data to the burn state.

    for (int n = 1; n <= NumSpec; ++n) {
        state.xn[n-1] = vode_state.y(n);
    }

    state.e = vode_state.y(net_ienuc);
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void burn_to_vode (const burn_t& state, dvode_t& vode_state)
{
    // Copy the integration data from the burn state.

    for (int n = 1; n <= NumSpec; ++n) {
        vode_state.y(n) = state.xn[n-1];
    }

    vode_state.y(net_ienuc) = state.e;
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void renormalize_species (dvode_t& vode_state)
{
    Real sum = 0.0_rt;

    for (int n = 1; n <= NumSpec; ++n) {
        sum += vode_state.y(n);
    }

    for (int n = 1; n <= NumSpec; ++n) {
        vode_state.y(n) /= sum;
    }
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void clean_state (dvode_t& vode_state)
{

    // Ensure that mass fractions always stay positive and less than or
    // equal to 1.

    for (int n = 1; n <= NumSpec; ++n) {
        vode_state.y(n) = amrex::max(amrex::min(vode_state.y(n), 1.0_rt), SMALL_X_SAFE);
    }

    // Renormalize the abundances as necessary.

    if (renormalize_abundances) {
        renormalize_species(vode_state);
    }

}

///
/// update the thermodynamics in the burn_t state -- this may involve an EOS call.
/// we also pass in the vode_state so we get the latest values of the mass fractions,
/// temperature, and internal energy
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void update_thermodynamics (burn_t& state, const dvode_t& vode_state)
{
#ifdef NSE_THERMO
    // we are handling the thermodynamics via the aux quantities, which
    // are stored in the burn_t here, so we need to update those based on
    // the current state.

    state.abar = 0.0_rt;
    state.y_e = 0.0_rt;
    for (int n = 1; n <= NumSpec; ++n) {
      state.abar += vode_state.y(n) * aion_inv[n-1];
      state.y_e += vode_state.y(n) * zion[n-1] * aion_inv[n-1];
    }
    state.abar = 1.0_rt / state.abar;
    state.zbar = state.abar * state.y_e;

    state.aux[AuxZero::iabar] = state.abar;
    state.aux[AuxZero::iye] = state.y_e;
    state.aux[AuxZero::ibea] = 0.0_rt; // this is not needed when we integrate

#endif

    for (int n = 1; n <= NumSpec; ++n) {
        state.xn[n-1] = vode_state.y(n);
    }
    state.e = vode_state.y(net_ienuc);

    // Get T from e (also updates composition quantities).

    if (call_eos_in_rhs) {
        eos(eos_input_re, state);
    }
}

#endif
