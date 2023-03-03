#ifndef INTEGRATOR_TYPE_SIMPLIFIED_SDC_H
#define INTEGRATOR_TYPE_SIMPLIFIED_SDC_H

#include <eos.H>
#include <eos_composition.H>
#include <burn_type.H>
#include <actual_rhs.H>
#if defined(SDC_EVOLVE_ENTHALPY)
#include <maestro_params.H>
#endif


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void update_density_in_time(const Real time, burn_t& state, [[maybe_unused]] const T& int_state)
{

    // some quantities are only advected.  Here we update those state
    // variables in burn_t's y[] array due only to advection.
    //
    // Note: we don't update the aux data here based on the X's -- that
    // needs to be done separately.
    //
    // Also note: int_state is only used as input for the case where
    // we need to construct rho from (rho X).

#if defined(SDC_EVOLVE_ENERGY)

    // we are always integrating from t = 0, so there is no offset
    // time needed here.  The indexing of ydot_a is based on
    // the indices in burn_t and is 0-based    
    state.y[SRHO] = amrex::max(state.rho_orig + state.ydot_a[SRHO] * time, EOSData::mindens);

    // for consistency
    state.rho = state.y[SRHO];

#elif defined(SDC_EVOLVE_ENTHALPY)

    // Keep density consistent with the partial densities.
    state.rho = 0.0_rt;
    for (int n = 1; n <= NumSpec; n++) {
        // use 1-based indexing into int_state.y
        state.rho += int_state.y(SFS+n);
    }
#endif

}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void renormalize_species(const Real time, burn_t& state, T& int_state)
{

    // Constrain the partial densities in int_state to sum to the
    // density.  This does not update the burn_t state.
    //
    // We only renormalize species when evolving energy because
    // when we evolve enthalpy, we define the density as
    // the sum of the partial densities rho*X for each species.

#ifdef SDC_EVOLVE_ENERGY

    // update rho, rho*u, etc.

    update_density_in_time(time, state, int_state);

    Real nspec_sum = 0.0_rt;
    for (int n = 1; n <= NumSpec; n++) {
        // use 1-based indexing
        nspec_sum += int_state.y(SFS+n);
    }
    nspec_sum /= state.y[SRHO];

    for (int n = 1; n <= NumSpec; n++) {
        int_state.y(SFS+n) /= nspec_sum;
    }

    // note: we don't carry the auxiliary data in int_state, so there
    // is no update to them done here

#endif

}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void clean_state(const Real time, burn_t& state, T& int_state)
{

    // Ensure that mass fractions always stay positive.

    if (do_species_clip) {
        for (int n = 1; n <= NumSpec; ++n) {
            // we use 1-based indexing, so we need to offset SFS
            int_state.y(SFS+n) = amrex::max(amrex::min(int_state.y(SFS+n), state.rho),
                                            state.rho * SMALL_X_SAFE);
        }
    }

    // renormalize abundances as necessary

    if (renormalize_abundances) {
        renormalize_species(time, state, int_state);
    }

}


// Given a burn state, fill the rpar and integration state data.
// this is intended only to be called once -- at the very start of the integration

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void burn_to_int(burn_t& state, T& int_state)
{

    for (int n = 0; n < SVAR_EVOLVE; n++) {
        // int_state uses 1-based indexing
        int_state.y(n+1) = state.y[n];
    }

    // store the original rho and rho e
#if defined(SDC_EVOLVE_ENERGY)

    state.rho_orig = state.y[SRHO];

    state.rhoe_orig = state.y[SEINT];

#else

    state.rho_orig = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        state.rho_orig += state.y[SFS+n];
    }

#endif

}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void int_to_burn(const Real time, const T& int_state, burn_t& state)
{

    // this makes burn_t represent the current integration state.  The
    // main thing this really does is compute the current temperature

    state.time = time;

    for (int n = 0; n < SVAR_EVOLVE; n++) {
        // note int_state uses 1-based indexing
        state.y[n] = int_state.y(n+1);
    }

    // update rho in the burn_t state
    // this may be redundant, but better to be safe

    update_density_in_time(time, state, int_state);

    Real rhoInv = 1.0_rt / state.rho;

#if defined(SDC_EVOLVE_ENERGY)
    for (int n = 0; n < NumSpec; n++) {
        // int_state uses 1-based indexing
        state.xn[n] = int_state.y(SFS+1+n) * rhoInv;
    }

#ifdef AUX_THERMO
    // make the aux data consistent with the int_state X's
    set_aux_comp_from_X(state);

    // also store it in the burn_t state solution
    for (int n = 0; n < NumAux; n++) {
        state.y[SFX+n] = state.rho * state.aux[n];
    }
#endif

    // set internal energy for EOS call

    state.e = int_state.y(SEINT+1) * rhoInv;

    state.T = state.T;  // initial guess

    eos(eos_input_re, state);

#elif defined(SDC_EVOLVE_ENTHALPY)
    eos_t eos_state;
    eos_state.rho = state.rho;
    for (int n = 0; n < NumSpec; n++) {
        // int_state uses 1-based indexing
        eos_state.xn[n] = int_state.y(SFS+1+n) * rhoInv;
    }

#ifdef AUX_THERMO
    // make the aux data consistent with the int_state X's
    set_aux_comp_from_X(eos_state);

    // also store it in the burn_t state
    for (int n = 0; n < NumAux; n++) {
        state.y[SFX+n] = state.rho * eos_state.aux[n];
    }
#endif

    if (maestro::use_tfromp) {
        // under this assumption, we are assuming that there is no base state evolution
        eos_state.p = state.p0;
    } else {
        eos_state.h = int_state.y(SENTH+1) * rhoInv;
    }

    eos_state.T = state.T;  // initial guess


    if (maestro::use_tfromp) {
        // NOT SURE IF THIS IS VALID
        // used to be an Abort statement
        eos(eos_input_rp, eos_state);
    } else {
        eos(eos_input_rh, eos_state);
    }

    // fill the rest of the burn_t state

    eos_to_burn(eos_state, state);

#endif

    // override T if we are fixing it (e.g. due to
    // drive_initial_convection)
    if (state.T_fixed > 0.0_rt) {
        state.T = state.T_fixed;
    }

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rhs_to_int([[maybe_unused]] const Real time,
                burn_t& state,
                RArray1D& ydot)
{


    // on input, ydot has the contributions to the RHS from the
    // reaction network.  Note that these are in terms of dY/dt
    // we will convert these in place to the form needed for SDC

    // convert from dY/dt to dX/dt. The species derivatives are the
    // first NumSpec components of ydot coming from the reaction network

    BL_ASSERT(SFS == 0);

    for (int n = 1; n <= NumSpec; n++) {
        ydot(n) = state.rho * aion[n-1] * ydot(n);
    }

    // rescale the energy to be d(rho e)/dt

#if defined(SDC_EVOLVE_ENERGY)

    BL_ASSERT(SEINT+1 == net_ienuc);

    ydot(SEINT+1) *= state.rho;

#elif defined(SDC_EVOLVE_ENTHALPY)

    ydot(SENTH+1) *= state.rho;

#endif

    // now add the contribution from the non-reacting sources --
    // including advection

    // Note: both ydot is 1-based

    for (int n = 0; n < SVAR_EVOLVE; n++) {
        ydot(n+1) += state.ydot_a[n];
    }

}

#endif
