#ifndef VODE_TYPE_SIMPLIFIED_SDC_H
#define VODE_TYPE_SIMPLIFIED_SDC_H

#include <eos.H>
#include <eos_composition.H>
#include <burn_type.H>
#ifdef NETWORK_HAS_CXX_IMPLEMENTATION
#include <actual_rhs.H>
#else
#include <fortran_to_cxx_actual_rhs.H>
#endif
#ifdef NSE_THERMO
#include <nse.H>
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void fill_unevolved_variables(const Real time, burn_t& state, dvode_t& vode_state)
{

    // some quantities are only advected.  Here we update those state
    // variables in burn_t's y[] array due only to advection.
    //
    // Note: we don't update the aux data here based on the X's -- that
    // needs to be done separately.
    //
    // Also note: vode_state is only used as input for the case where
    // we need to construct rho from (rho X).

#if defined(SDC_EVOLVE_ENERGY)

    // we are always integrating from t = 0, so there is no offset
    // time needed here.  The indexing of ydot_a is based on
    // the indices in burn_t and is 0-based
    state.y[SRHO] = state.rho_orig + state.ydot_a[SRHO] * time;

    // for consistency
    state.rho = state.y[SRHO];

    state.y[SMX] = state.umx_orig + state.ydot_a[SMX] * time;
    state.y[SMY] = state.umy_orig + state.ydot_a[SMY] * time;
    state.y[SMZ] = state.umz_orig + state.ydot_a[SMZ] * time;

#elif defined(SDC_EVOLVE_ENTHALPY)

    // Keep density consistent with the partial densities.
    state.rho = 0.0_rt;
    for (int n = 1; n <= NumSpec; n++) {
        // use 1-based indexing into vode_state.y
        state.rho += vode_state.y[SFS+n];
    }
#endif

}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void sdc_burn_to_eos (const burn_t& state, T& eos_state)
{
#if defined(SDC_EVOLVE_ENERGY)
    eos_state.rho = state.y[SRHO];

#elif defined(SDC_EVOLVE_ENTHALPY)
    eos_state.rho = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.rho += state.y[SFS+n];
    }
#endif

    eos_state.T = state.T;

    // species

    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state.y[SFS+n] / eos_state.rho;
    }

#if NAUX_NET > 0
    // aux
    for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = state.y[SFX+n] / eos_state.rho;
    }
    //set_nse_aux_from_X(eos_state);
#endif

    // we don't bother filling the other fields, since the EOS call
    // will do that for us
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void renormalize_species(const Real time, burn_t& state, dvode_t& vode_state)
{

    // Constrain the partial densities in vode_state to sum to the
    // density.  This does not update the burn_t state.
    //
    // We only renormalize species when evolving energy because
    // when we evolve enthalpy, we define the density as
    // the sum of the partial densities rho*X for each species.

#ifdef SDC_EVOLVE_ENERGY

    // update rho, rho*u, etc.

    fill_unevolved_variables(time, state, vode_state);

    Real nspec_sum = 0.0_rt;
    for (int n = 1; n <= NumSpec; n++) {
        // use 1-based indexing
        nspec_sum += vode_state.y(SFS+n);
    }
    nspec_sum /= state.y[SRHO];

    for (int n = 1; n <= NumSpec; n++) {
        vode_state.y(SFS+n) /= nspec_sum;
    }

    // note: we don't carry the auxiliary data in vode_state, so there
    // is no update to them done here

#endif

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void clean_state(const Real time, burn_t& state, dvode_t& vode_state)
{

    // This routine makes the state we are integrating in vode_state
    // consistent We will need to use some of the unevolved data from
    // state, but only vode_state is updated here

    // update rho, rho*u, etc.

    fill_unevolved_variables(time, state, vode_state);

    // Ensure that mass fractions always stay positive.

    for (int n = 1; n <= NumSpec; ++n) {
        // we use 1-based indexing, so we need to offset SFS
        vode_state.y(SFS+n) = amrex::max(amrex::min(vode_state.y(SFS+n), state.rho),
                                         state.rho * 1.e-200_rt);
    }


    // renormalize abundances as necessary

    if (renormalize_abundances) {
        renormalize_species(time, state, vode_state);
    }

#ifdef SDC_EVOLVE_ENERGY

    // Ensure that internal energy never goes above the maximum limit
    // provided by the EOS. Same for the internal energy implied by the
    // total energy (which we get by subtracting kinetic energy).

    eos_t eos_state;
    eos_state.rho = state.y[SRHO];
    eos_state.T = MAX_TEMP;
    for (int n = 0; n < NumSpec; n++) {
        // vode_state.y uses 1-based indexing
        eos_state.xn[n] = vode_state.y(SFS+1+n) / state.rho;
    }
#ifdef NSE_THERMO
    // make the aux data here to be consistent with the renormalized
    // abundances
    set_nse_aux_from_X(eos_state);
#endif

    eos(eos_input_rt, eos_state);

    Real max_e = eos_state.e;

    vode_state.y(SEINT+1) = amrex::min(state.rho * max_e, vode_state.y(SEINT+1));

    Real ke = vode_state.y(SEDEN+1) -
        0.5_rt * (state.y[SMX] * state.y[SMX] +
                  state.y[SMY] * state.y[SMY] +
                  state.y[SMZ] * state.y[SMZ]) / state.rho;

    vode_state.y(SEDEN+1) = amrex::min(state.rho * max_e + ke, vode_state.y(SEDEN+1));

#endif

}


// Given a burn state, fill the rpar and integration state data.
// this is intended only to be called once -- at the very start of the integration

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void burn_to_vode(burn_t& state, dvode_t& vode_state)
{

    for (int n = 0; n < SVAR_EVOLVE; n++) {
        // vode_state uses 1-based indexing
        vode_state.y(n+1) = state.y[n];
    }

    // store the original rho and rho U
#if defined(SDC_EVOLVE_ENERGY)

    state.rho_orig = state.y[SRHO];

    state.umx_orig = state.y[SMX];
    state.umy_orig = state.y[SMY];
    state.umz_orig = state.y[SMZ];

#else

    state.rho_orig = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        state.rho_orig += state.y[SFS+n];
    }

#endif

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void vode_to_burn(const Real time, dvode_t& vode_state, burn_t& state)
{

    // this makes burn_t represent the current integration state.  The
    // main thing this really does is compute the current temperature

    for (int n = 0; n < SVAR_EVOLVE; n++) {
        // note vode_state uses 1-based indexing
        state.y[n] = vode_state.y(n+1);
    }

    // update rho, rho*u, ... in the burn_t state

    fill_unevolved_variables(time, state, vode_state);


    Real rhoInv = 1.0_rt / state.rho;

    eos_t eos_state;
    eos_state.rho = state.rho;
    for (int n = 0; n < NumSpec; n++) {
        // vode_state uses 1-based indexing
        eos_state.xn[n] = vode_state.y(SFS+1+n) * rhoInv;
    }
#ifdef NSE_THERMO
    // make the aux data consistent with the vode_state X's
    set_nse_aux_from_X(eos_state);

    // also store it in the burn_t state
    for (int n = 0; n < NumAux; n++) {
        state.y[SFX+n] = state.rho * eos_state.aux[n];
    }
#endif


#if defined(SDC_EVOLVE_ENERGY)

    if (state.T_from_eden) {
        eos_state.e = (vode_state.y(SEDEN+1) -
                       0.5_rt * rhoInv * (state.y[SMX] * state.y[SMX] +
                                          state.y[SMY] * state.y[SMY] +
                                          state.y[SMZ] * state.y[SMZ])) * rhoInv;
    } else {
        eos_state.e = vode_state.y(SEINT+1) * rhoInv;
    }

#elif defined(SDC_EVOLVE_ENTHALPY)

    if (use_tfromp) {
        // under this assumption, we are assuming that there is no base state evolution
        eos_state.p = state.p0;
    } else {
        eos_state.h = vode_state.y(SENTH+1) * rhoInv;
    }

#endif

    // Give the temperature an initial guess -- use the geometric mean
    // of the minimum and maximum temperatures.

    // TODO: we should just change this to use state.T

    Real min_temp = EOSData::mintemp;
    Real max_temp = EOSData::maxtemp;

    eos_state.T = std::sqrt(min_temp * max_temp);

#if defined(SDC_EVOLVE_ENERGY)
    eos(eos_input_re, eos_state);

#elif defined(SDC_EVOLVE_ENTHALPY)

    if (use_tfromp) {
        // NOT SURE IF THIS IS VALID
        // used to be an Abort statement
        eos(eos_input_rp, eos_state);
    } else {
        eos(eos_input_rh, eos_state);
    }

#endif

    // fill the rest of the burn_t state

    eos_to_burn(eos_state, state);

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rhs_to_vode(const Real time, burn_t& state,
                 YdotNetArray1D& ydot_react, dvode_t& vode_state,
                 RArray1D& ydot)
{


    fill_unevolved_variables(time, state, vode_state);

    // ydot_react has just the contribution to the RHS from the
    // reaction network.  Note that these are in terms of dY/dt

    // start with the contribution from the non-reacting sources

    // Note: both ydot and ydot_react are 1-based
    for (int n = 0; n < SVAR_EVOLVE; n++) {
        ydot(n+1) = state.ydot_a[n];
    }

    // add in the reacting terms -- here we convert from dY/dt to dX/dt
    // The species derivatives are the first NumSpec components of
    // ydot_react

    for (int n = 1; n <= NumSpec; n++) {
        ydot(SFS+n) += state.rho * aion[n-1] * ydot_react(n);
    }

#if defined(SDC_EVOLVE_ENERGY)

    ydot(SEINT+1) += state.rho * ydot_react(net_ienuc);
    ydot(SEDEN+1) += state.rho * ydot_react(net_ienuc);

#elif defined(SDC_EVOLVE_ENTHALPY)

    ydot(SENTH+1) += state.rho * ydot_react(net_ienuc);

#endif

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void jac_to_vode(const Real time, burn_t& state,
                 JacNetArray2D& jac_react, dvode_t& vode_state,
                 RArray2D& jac)
{

    // this is only used with an analytic Jacobian.  At the moment, we
    // only support a dense Jacobian.


    // we come in with jac_react being the Jacobian of the reacting system
    // but we need to convert it to the SDC system


    // SVAR_EVOLVE doesn't include rho, but we will include it here in
    // the intermediate this affects both the Castro
    // (SDC_EVOLVE_ENERGY) and MAESTROeX (SDC_EVOLVE_ENTHALPY) systems.

#if defined(SDC_EVOLVE_ENERGY)
    constexpr int iwrho = 1;
    constexpr int iwfs = 2;
    constexpr int iwK = iwfs+NumSpec;
    constexpr int iwT = iwK+1;
    constexpr int iwvar = 3+NumSpec;
#else
    constexpr int iwrho = 1;
    constexpr int iwfs=2;
    constexpr int iwT = iwfs+NumSpec;
    constexpr int iwvar = 2+NumSpec;
#endif

    amrex::Array2D<Real, 1, SVAR_EVOLVE+1, 1, iwvar> dRdw = {0.0_rt};
    amrex::Array2D<Real, 1, iwvar, 1, SVAR_EVOLVE+1> dwdU = {0.0_rt};

    constexpr Real eps = 1.e-8_rt;

    // this is 0-based to be consistent with SFS, SEDEN, ...
    constexpr int SRHO_EXTRA = SVAR_EVOLVE;


    // jac_react has the derivatives with respect to the native
    // network variables, X, T. e.  It does not have derivatives with
    // respect to density, so we'll have to compute those ourselves.

    // The Jacobian from the nets is in terms of dYdot/dY, but we want
    // it was dXdot/dX, so convert here.
    for (int n = 1; n <= NumSpec; n++) {
        for (int m = 1; m <= neqs; m++) {
            jac_react(n,m) = jac_react(n,m) * aion[n-1];
        }
    }

    for (int m = 1; m <= neqs; m++) {
        for (int n = 1; n <= NumSpec; n++) {
            jac_react(m,n) = jac_react(m,n) * aion_inv[n-1];
        }
    }

    // also fill the ydot
    YdotNetArray1D ydot;
    vode_to_burn(time, vode_state, state);
    actual_rhs(state, ydot);

    // at this point, our Jacobian should be entirely in terms of X,
    // not Y.  Let's now fix the rhs terms themselves to be in terms of
    // dX/dt and not dY/dt.
    for (int n = 1; n <= NumSpec; n++) {
        ydot(n) = ydot(n) * aion[n-1];
    }

    // now perturb density and call the RHS to compute the derivative wrt rho
    // species rates come back in terms of molar fractions
    burn_t state_pert = state;
    state_pert.rho = state.rho * (1.0_rt + eps);

    YdotNetArray1D ydot_pert;
    actual_rhs(state_pert, ydot_pert);

    // make the rates dX/dt and not dY/dt
    for (int n = 1; n <= NumSpec; n++) {
        ydot_pert(n) = ydot_pert(n) * aion[n-1];
    }

#if defined(SDC_EVOLVE_ENERGY)

    // The system we integrate has the form (rho X_k, rho E, rho e), but we will temporarily augment
    // This with rho, giving U = (rho, rho X_k, rho E, rho e).
    //
    // The intermediate state, w, has the form w = (rho, X_k, K, T), where K is 1/2 |U|^2

    // First compute dR/dw using the Jacobian that comes from the
    // network.  Note: this doesn't include density derivatives, so we
    // compute those via differencing.

    // dR/dw has the form:
    //
    //  SFS         / d(rho X1dot)/drho  d(rho X1dot)/dX1   d(rho X1dit)/dX2   ...  0   d(rho X1dot)/dT \
    //              | d(rho X2dot)/drho  d(rho X2dot)/dX1   d(rho X2dot)/dX2   ...  0   d(rho X2dot)/dT |
    //  SFS-1+nspec |   ...                                                         0                   |
    //  SEINT       | d(rho Edot)/drho   d(rho Edot)/dX1    d(rho Edot)/dX2    ...  0   d(rho Edot)/dT  |
    //  SEDEN       | d(rho Edot)/drho   d(rho Edot)/dX1    d(rho Edot)/dX2    ...  0   d(rho Edot)/dT  |
    //  SRHO_EXTRA  \        0                  0                  0                0          0       /
    //
    //                                                                              ^
    //                                                                              K derivatives


    // fill the column of dRdw corresponding to the derivative
    // with respect to rho

    // keep in mind here that that we are using 1-based indexing but SFS, ... are 0-based.

    for (int m = 1; m <= NumSpec; m++) {
        // d( d(rho X_m)/dt)/drho
        dRdw(SFS+m, iwrho) = ydot(m) + state.rho * (ydot_pert(m) - ydot(m))/(eps * state.rho);
    }

    // d( d(rho e)/dt)/drho
    dRdw(SEINT+1, iwrho) = ydot(net_ienuc) +
        state.rho * (ydot_pert(net_ienuc) - ydot(net_ienuc))/(eps * state.rho);

    // d( d(rho E)/dt)/drho
    dRdw(SEDEN+1, iwrho) = dRdw(SEINT+1, iwrho);

    // fill the columns of dRdw corresponding to each derivative
    // with respect to species mass fraction
    for (int n = 1; n <= NumSpec; n++) {
        for (int m = 1; m <= NumSpec; m++) {
            // d( d(rho X_m)/dt)/dX_n
            dRdw(SFS+m, iwfs-1+n) = state.rho * jac_react(m, n);
        }

        // d( d(rho e)/dt)/dX_n
        dRdw(SEINT+1, iwfs-1+n) = state.rho * jac_react(net_ienuc, n);

        // d( d(rho E)/dt)/dX_n
        dRdw(SEDEN+1, iwfs-1+n) = state.rho * jac_react(net_ienuc, n);
    }

    // now fill the column corresponding to derivatives with respect to
    // temperature -- this column is iwT

    // d( d(rho X_m)/dt)/dT
    for (int m = 1; m <= NumSpec; m++) {
        dRdw(SFS+m, iwT) = state.rho * jac_react(m, net_itemp);
    }

    // d( d(rho e)/dt)/dT
    dRdw(SEINT+1, iwT) = state.rho * jac_react(net_ienuc, net_itemp);

    // d( d(rho E)/dt)/dT
    dRdw(SEDEN+1, iwT) = dRdw(SEINT+1, iwT);

    // for the K derivatives, dRdw(:, iwK), and the rho sources,
    // dRdw(SRHO_EXTRA, :), we don't need to do anything, because these
    // are already zeroed out

    // that completes dRdw


    // construct dwdU

    // kinetic energy, K = 1/2 |U|^2
    Real kineng = 0.5_rt * (state.y[SMX] * state.y[SMX] +
                            state.y[SMY] * state.y[SMY] +
                            state.y[SMZ] * state.y[SMZ]) /
        (state.rho * state.rho);

    // density row (iwrho)
    dwdU(iwrho, SRHO_EXTRA+1) = 1.0_rt;

    // species rows
    for (int m = 1; m <= NumSpec; m++) {
        dwdU(iwfs-1+m, SFS+m) = 1.0_rt / state.rho;
        dwdU(iwfs-1+m, SRHO_EXTRA+1) = -state.xn[m-1] / state.rho;
    }

    // K row
    dwdU(iwK, SRHO_EXTRA+1) = -kineng / state.rho;
    dwdU(iwK, SEINT+1) = -1.0_rt / state.rho;
    dwdU(iwK, SEDEN+1) = 1.0_rt / state.rho;

    // T row
    eos_t eos_state;
    eos_state.rho = state.rho;
    eos_state.T = 1.e4_rt;   // initial guess
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = vode_state.y(SFS+1+n) / state.rho;
    }
#ifdef NSE_THERMO
    set_nse_aux_from_X(eos_state);
#endif

    eos_state.e = vode_state.y(SEINT+1) / state.rho;

    eos(eos_input_re, eos_state);

    const eos_xderivs_t eos_xderivs = composition_derivatives(eos_state);

    // temperature row
    for (int n = 1; n <= NumSpec; n++) {
        dwdU(iwT, SFS+n) = -eos_xderivs.dedX[n-1]/  (state.rho * eos_state.dedT);
    }
    dwdU(iwT, SEINT+1) = 1.0_rt / (state.rho * eos_state.dedT);
    dwdU(iwT, SEDEN+1) = 0.0_rt;
    Real X_dedX_sum = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        X_dedX_sum += eos_state.xn[n] * eos_xderivs.dedX[n];
    }
    dwdU(iwT, SRHO_EXTRA+1) =
        (X_dedX_sum - state.rho * eos_state.dedr - eos_state.e) /
        (state.rho * eos_state.dedT);

#elif defined(SDC_EVOLVE_ENTHALPY)

    // Our R source has components for species and enthalpy only.  But
    // we will extend it here to include the mass density too to ensure
    // that we have a square matrix in dU/dw that we can take the
    // inverse of to use below.  When we compute the final Jacobian, we will
    // discard the density row.

    // Our jacobian, dR/dw has the form:
    //
    //  SFS         / d(rho X1dot)/drho  d(rho X1dot)/dX1   d(rho X1dit)/dX2   ...  d(rho X1dot)/dT \
    //              | d(rho X2dot)/drho  d(rho X2dot)/dX1   d(rho X2dot)/dX2   ...  d(rho X2dot)/dT |
    //  SFS-1+nspec |   ...                                                                         |
    //  SENTH       | d(rho h)/drho      d(rho h)/dX1       d(rho h)/dX2       ...  d(rho h)/dT     |
    //  SRHO_EXTRA  \ 0                  0                  0                       0               /


    // fill the column of dRdw corresponding to the derivative
    // with respect to rho
    for (int m = 1; m <= NumSpec; m++) {
        // d( d(rho X_m)/dt)/drho
        dRdw(SFS+m, iwrho) = ydot(m) + state.rho * (ydot_pert(m) - ydot(m))/(eps * state.rho);
    }

    // d( d(rho h)/dt)/drho
    dRdw(SENTH+1, iwrho) = ydot(net_ienuc) +
        state.rho * (ydot_pert(net_ienuc) - ydot(net_ienuc))/(eps * state.rho);

    // d( d(rho)/dt)/drho
    dRdw(SRHO_EXTRA+1, iwrho) = 0.0_rt;

    // fill the columns of dRdw corresponding to each derivative
    // with respect to species mass fraction
    for (int n = 1; n <= NumSpec; n++) {
        for (int m = 1; m <= NumSpec; m++) {
            // d( d(rho X_m)/dt)/dX_n
            dRdw(SFS+m, iwfs-1+n) = state.rho * jac_react(m, n);
        }

        // d( d(rho h)/dt)/dX_n
        dRdw(SENTH+1, iwfs-1+n) = state.rho * jac_react(net_ienuc, n);

        // d( d(rho)/dt)/dX_n
        dRdw(SRHO_EXTRA+1, iwfs-1+n) = 0.0_rt;
    }

    // now fill the column corresponding to derivatives with respect to
    // temperature -- this column is iwT

    // d( d(rho X_m)/dt)/dT
    for (int m = 1; m <= NumSpec; m++) {
        dRdw(SFS+m, iwT) = state.rho * jac_react(m, net_itemp);
    }

    // d( d(rho h)/dt)/dT
    dRdw(SENTH+1, iwT) = state.rho * jac_react(net_ienuc, net_itemp);

    // d( d(rho)/dt)/dT
    dRdw(SRHO_EXTRA+1, iwT) = 0.0_rt;

    // that completes dRdw

    // construct dwdU.  Here we take U = (rho X, rho h, rho)^T

    // density row (iwrho)
    dwdU(iwrho, SRHO_EXTRA+1) = 1.0_rt;

    // species rows
    for (int m = 1; m <= NumSpec; m++) {
        dwdU(iwfs-1+m, SFS+m) = 1.0_rt / state.rho;
        dwdU(iwfs-1+m, SRHO_EXTRA+1) = -state.xn[m-1] / state.rho;
    }

    // T row
    eos_t eos_state;
    eos_state.rho = state.rho;
    eos_state.T = 1.e4_rt;   // initial guess
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = vode_state.y(SFS+1+n) / state.rho;
    }
#ifdef NSE_THERMO
    set_nse_aux_from_X(eos_state);
#endif

    eos_state.h = vode_state.y(SENTH+1) / state.rho;

    eos(eos_input_rh, eos_state);

    const eos_xderivs_t eos_xderivs = composition_derivatives(eos_state);

    // temperature row
    for (int n = 1; n <= NumSpec; n++) {  
        dwdU(iwT, SFS+n) = -eos_xderivs.dhdX[n-1] / (state.rho * eos_state.dedT);
    }
    dwdU(iwT, SENTH+1) = 1.0_rt / (state.rho * eos_state.dhdT);
    Real X_dhdX_sum = 0.0;
    for (int n = 0; n < NumSpec; n++) {
        X_dhdX_sum += eos_state.xn[n] * eos_xderivs.dhdX[n];
    }
    dwdU(iwT, SRHO_EXTRA+1) =
        (X_dhdX_sum - state.rho * eos_state.dhdr - eos_state.h) /
        (state.rho * eos_state.dhdT);

#endif


    // compute J = dR/dw dw/dU

    // J is SVAR_EVOLVE x SVAR_EVOLVE, which will call m x n
    //
    // J = dR/dw dw/dU
    //
    //   dR/dw is SVAR_EVOLVE+1 x iwvar, which we call m x k
    //   dw/dU is iwvar x SVAR_EVOLVE+1, which we call k x n
    //

    // we need to cut out the density (SRHO_EXTRA) row and column of
    // the Jacobian, since that is not in our full SVAR_EVOLVE state
    for (int n = 1; n <= SVAR_EVOLVE; n++) {
        if (n == SRHO_EXTRA+1) continue;
        for (int m = 1; m <= SVAR_EVOLVE; m++) {
            if (m == SRHO_EXTRA+1) continue;

            jac(m, n) = 0.0_rt;
            for (int k = 1; k <= iwvar; k++) {
                jac(m, n) += dRdw(m, k) * dwdU(k, n);
            }
        }
    }

}

#endif
