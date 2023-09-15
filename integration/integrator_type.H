#ifndef INTEGRATOR_TYPE_H
#define INTEGRATOR_TYPE_H

#include <eos.H>

using namespace integrator_rp;

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void update_density_in_time(const Real time, BurnT& state)
{

    // some quantities are only advected.  Here we update those state
    // variables in burn_t's y[] array due only to advection.
    //
    // Note: we don't update the aux data here based on the X's -- that
    // needs to be done separately.
    //
    // we are always integrating from t = 0, so there is no offset
    // time needed here.  The indexing of ydot_a is based on
    // the indices in burn_t and is 0-based
    state.y[SRHO] = amrex::max(state.rho_orig + state.ydot_a[SRHO] * time, EOSData::mindens);

    // for consistency
    state.rho = state.y[SRHO];

}


///
/// update the thermodynamics in the burn_t state -- this may involve an EOS call.
/// we also pass in the int_state so we get the latest values of the mass fractions,
/// temperature, and internal energy
///
template <typename BurnT, typename I>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void update_thermodynamics (BurnT& state, const I& int_state)
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

    if (scale_system) {
        state.e *= state.e_scale;
    }

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
