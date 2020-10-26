#ifndef VODE_TYPE_STRANG_H
#define VODE_TYPE_STRANG_H


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void vode_to_burn (const dvode_t& vode_state, burn_t& state)
{
    // Copy the integration data to the burn state.

    for (int n = 1; n <= NumSpec; ++n) {
        state.xn[n-1] = vode_state.y(n);
    }

    state.T = vode_state.y(net_itemp);
    state.e = vode_state.y(net_ienuc);
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void burn_to_vode (const burn_t& state, dvode_t& vode_state)
{
    // Copy the integration data from the burn state.

    for (int n = 1; n <= NumSpec; ++n) {
        vode_state.y(n) = state.xn[n-1];
    }

    vode_state.y(net_itemp) = state.T;
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

    // Ensure that the temperature always stays within reasonable limits.

    vode_state.y(net_itemp) = amrex::min(MAX_TEMP, amrex::max(vode_state.y(net_itemp), EOSData::mintemp));

}

///
/// update the thermodynamics in the burn_t state -- this may involve an EOS call.
/// we also pass in the vode_state so we get the latest values of the mass fractions,
/// temperature, and internal energy
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void update_thermodynamics (burn_t& state, const dvode_t& vode_state)
{

    // Fill an EOS state using a combination of the burn data
    // (which has supplementary data such as density) and the
    // VODE state (which has the current xn, T, and e).

    eos_t eos_state;

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

    burn_to_eos(state, eos_state);

    for (int n = 1; n <= NumSpec; ++n) {
        eos_state.xn[n-1] = vode_state.y(n);
    }
    eos_state.T = vode_state.y(net_itemp);
    eos_state.e = vode_state.y(net_ienuc);

    // Evaluate the thermodynamics -- if desired. Note that we can
    // either do this as an unconditional update (call_eos_in_rhs)
    // or we optionally set a fraction dT_crit such that we don't
    // call the EOS if the last temperature we evaluated the EOS
    // at is relatively close to the current temperature. Otherwise
    // just do the composition calculations since that's needed to
    // construct dX/dt.

    if (call_eos_in_rhs && state.self_heat) {

        eos(eos_input_rt, eos_state);

    }
    else if (std::abs(eos_state.T - state.T_old) > dT_crit * eos_state.T && state.self_heat)
    {

        eos(eos_input_rt, eos_state);

        state.dcvdT = (eos_state.cv - state.cv_old) / (eos_state.T - state.T_old);
        state.dcpdT = (eos_state.cp - state.cp_old) / (eos_state.T - state.T_old);

        state.T_old = eos_state.T;
        state.cv_old = eos_state.cv;
        state.cp_old = eos_state.cp;

    }
    else {

        composition(eos_state);

    }

    eos_to_burn(eos_state, state);

}

#endif
