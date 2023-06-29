#ifndef INTEGRATOR_TYPE_STRANG_H
#define INTEGRATOR_TYPE_STRANG_H

#include <eos.H>
#include <integrator_type.H>

using namespace integrator_rp;

template <typename I, typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void integrator_to_burn (const I& int_state, BurnT& state)
{
    // Copy the integration data to the burn state.

    for (int n = 1; n <= NumSpec; ++n) {
        state.xn[n-1] = int_state.y(n);
    }

    state.e = int_state.y(net_ienuc);

    if (scale_system) {
        state.e *= state.e_scale;
    }

    state.time = int_state.t;
}


template <typename BurnT, typename I>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void burn_to_integrator (const BurnT& state, I& int_state)
{
    // Copy the integration data from the burn state.

    for (int n = 1; n <= NumSpec; ++n) {
        int_state.y(n) = state.xn[n-1];
    }

    int_state.y(net_ienuc) = state.e;

    if (scale_system) {
        int_state.y(net_ienuc) /= state.e_scale;
    }
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

#endif
