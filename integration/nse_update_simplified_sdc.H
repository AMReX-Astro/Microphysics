#ifndef NSE_UPDATE_H
#define NSE_UPDATE_H

#include <iostream>
#include <fstream>
#include <actual_network.H>

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include <extern_parameters.H>

#include <burn_type.H>
#include <eos.H>

#ifdef NSE_TABLE
#include <nse_table.H>
#endif
#ifdef NSE_NET
#include <nse_solver.H>
#endif

using namespace amrex;

#if defined(NSE_TABLE)

///
/// this acts as an explicit Euler step for the system (rho e, rho aux)
/// on input, *_source are the reactive sources at time t0 and on output
/// they are the sources at time t0+dt
///
template <bool skip_X_fill>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void evolve(const Real rho0, const Real T0, const Real rhoe0, const Real *rhoaux0,
            const Real dt, const Real *ydot_a,
            Real& rhoe_new, Real* rhoaux_new, Real* rhoX_new,
            Real& rhoe_source, Real* rhoaux_source, const Real T_fixed) {

    // integrate rho to the new time

    Real rho_new = rho0 + dt * ydot_a[SRHO];

    // integrate rhoe and rhoaux using the input sources

    rhoe_new = rhoe0 + dt * ydot_a[SEINT] + dt * rhoe_source;
    for (int n = 0; n < NumAux; n++) {
        rhoaux_new[n] = rhoaux0[n] + dt * ydot_a[SFX+n] + dt * rhoaux_source[n];
    }

    // find the temperature at t0 + dt

    eos_t eos_state;
    eos_state.rho = rho_new;
    eos_state.e = rhoe_new / rho_new;
    eos_state.T = T0;
    for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = rhoaux_new[n] / rho_new;
    }

    Real T_new;
    if (T_fixed > 0) {
        T_new = T_fixed;
    } else {
        eos(eos_input_re, eos_state);
        T_new = eos_state.T;
    }

    // compute the plasma neutrino losses at t0 + dt

    Real snu{0.0};
    Real dsnudt{0.0};
    Real dsnudd{0.0};
    Real dsnuda{0.0};
    Real dsnudz{0.0};

    Real abar = eos_state.aux[iabar];
    Real zbar = abar * eos_state.aux[iye];

    constexpr int do_derivatives = 0;
    sneut5<do_derivatives>(T_new, rho_new, abar, zbar,
                           snu, dsnudt, dsnudd, dsnuda, dsnudz);

    // call the NSE table

    Real abar_out;
    Real bea_out;
    Real dyedt;
    Real dabardt;
    Real dbeadt;
    Real enu;

    nse_interp(T_new, rho_new, eos_state.aux[iye],
               abar_out, bea_out, dyedt, dabardt, dbeadt, enu, rhoX_new,
               skip_X_fill);

    if constexpr (!skip_X_fill) {
        for (int n = 0; n < NumSpec; ++n) {
            rhoX_new[n] *= rho_new;
        }
    }

    // construct the sources at time t0+dt

    // note that abar, (B/A) here have already seen advection
    // implicitly

    Real rho_bea_tilde = rho_new * bea_out - dt * ydot_a[SFX+ibea];
    Real rho_dBEA = rho_bea_tilde - rhoaux0[ibea]; // this is MeV / nucleon * g / cm**3

    rhoe_source = rho_dBEA * C::MeV2eV * C::ev2erg * C::n_A / dt;
    rhoe_source -= rho_new * (enu + snu);

    rhoaux_source[iabar] = 0.0;
    rhoaux_source[iye] = rho_new * dyedt;
    rhoaux_source[ibea] = rho_dBEA / dt;

    // reset the new <B/A> -- we want to ensure that if we use the
    // updated T(rho, e, Y_e, Abar) that we get exactly the same B/A,
    // so there is no instantaneous change in the binding energy

    rhoaux_new[ibea] = rho_new * bea_out;

}

///
/// update the state due to NSE changes for the simplified-SDC case
/// this version works with the tabulated NSE and requires AUX_THERMO
///
template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void sdc_nse_burn(BurnT& state, const Real dt) {

    using namespace AuxZero;

    state.success = true;
    state.n_rhs = 0;
    state.n_jac = 0;

    // we need the initial Ye to get the NSE state.  We also
    // want the initial rho <B/A> since we may not be entering
    // this routine already in NSE.  This means that there will
    // be an energy release from us instantaneously adjusting
    // into NSE (the first call to nse_interp) + an energy
    // release from the evolution of NSE over the timestep

    // we will write the energy update as:
    //
    // e^new = e^old + dbea_instantantous + dbea_over_dt
    //
    // where dbea_over_dt will be an integral over dt using
    // a simple RK2 method

    Real ye_in = state.y[SFX+iye] / state.rho;
    Real rho_bea_orig = state.y[SFX+ibea];
    Real rhoe_orig = state.y[SEINT];

    // density and momentum have no reactive sources
    Real rho_old = state.y[SRHO];

    state.y[SRHO] += dt * state.ydot_a[SRHO];
    state.y[SMX] += dt * state.ydot_a[SMX];
    state.y[SMY] += dt * state.ydot_a[SMY];
    state.y[SMZ] += dt * state.ydot_a[SMZ];

    // if we are doing drive_initial_convection, we want to use
    // the temperature that comes in through T_fixed

    Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;

    // get the instantaneous change in B/A, as well as the source
    // terms for the new evolution.

    Real abar_out{};
    Real bea_out{};
    Real dyedt{};
    Real dabardt{};
    Real dbeadt{};
    Real enu{};
    Real X[NumSpec] = {};

    {
        constexpr bool skip_X_fill{true};

        nse_interp(T_in, state.rho, ye_in,
                   abar_out, bea_out, dyedt, dabardt, dbeadt, enu, X, skip_X_fill);
    }

    // the starting state for integration uses the state that we just
    // adjusted instantaneously into NSE

    Real rhoaux_old[NumAux];
    rhoaux_old[iye] = rho_old * ye_in;
    rhoaux_old[iabar] = rho_old * abar_out;
    rhoaux_old[ibea] = rho_old * bea_out;

    Real rhoenuc = (rho_old * bea_out - rho_bea_orig) * C::MeV2eV * C::ev2erg * C::n_A;
    Real rhoe_old = state.y[SEINT] + rhoenuc;


    // get the starting temperature

    eos_t eos_state;
    eos_state.rho = rho_old;
    eos_state.e = rhoe_old / rho_old;
    for (int n = 0; n < NumAux; ++n) {
        eos_state.aux[n] = rhoaux_old[n] / rho_old;
    }
    eos_state.T = T_in;

    Real T_old;
    if (state.T_fixed > 0) {
        T_old = state.T_fixed;
    } else {
        eos(eos_input_re, eos_state);
        T_old = eos_state.T;
    }

    // start of integration

    // we want to integrate from t^n to t^{n+1/2}
    // start by computing the sources

    // the aux sources are just due to weak rates
    Real rhoaux_source[NumAux];
    rhoaux_source[iye] = rho_old * dyedt;
    rhoaux_source[ibea] = rho_old * dbeadt;
    rhoaux_source[iabar] = 0.0;

    // the energy source is due to plasma neutrinos + weak rate neutrinos

    Real snu{};
    Real dsnudt{};
    Real dsnudd{};
    Real dsnuda{};
    Real dsnudz{};
    constexpr int do_derivatives = 0;

    Real zbar = abar_out * ye_in;
    sneut5<do_derivatives>(T_old, rho_old, abar_out, zbar,
                           snu, dsnudt, dsnudd, dsnuda, dsnudz);

    Real rhoe_source = -rho_old * (snu + enu);

    Real rhoe_new;
    Real rhoaux_new[NumAux];
    Real rhoX_new[NumSpec];

    // evolve for dt/2 -- the sources will be overwritten with the
    // sources at the end point of this integration

    {
        constexpr bool skip_X_fill{true};

        evolve<skip_X_fill>(rho_old, T_old, rhoe_old, rhoaux_old,
                            0.5_rt*dt, state.ydot_a,
                            rhoe_new, rhoaux_new, rhoX_new,
                            rhoe_source, rhoaux_source, state.T_fixed);
    }

    // this gives us the sources at n+1/2.  Now we do the full dt step,
    // starting again from the beginning, but using the time-centered
    // sources

    {
        constexpr bool skip_X_fill{false};

        evolve<skip_X_fill>(rho_old, T_old, rhoe_old, rhoaux_old,
                            dt, state.ydot_a,
                            rhoe_new, rhoaux_new, rhoX_new,
                            rhoe_source, rhoaux_source, state.T_fixed);
    }

    // done with the integration.  This gave us a second-order
    // accurate update to the rhoe and Ye, and the aux data and mass
    // fractions are consistent with the NSE table at the end of the
    // integration

    // now update the aux quantities

    // the new mass fractions are just those that come from the table
    // make sure they are normalized
    Real sum_rhoX{0.0_rt};
    for (auto & rx : rhoX_new) {
        rx = amrex::max(state.y[SRHO] * small_x, amrex::min(state.y[SRHO], rx));
        sum_rhoX += rx;
    }

    for (auto & rx : rhoX_new) {
        rx /= sum_rhoX;
    }

    for (int n = 0; n < NumSpec; ++n) {
        state.y[SFS+n] = rhoX_new[n];
    }

    // aux data comes from the integration

    for (int n = 0; n < NumAux; ++n) {
        state.y[SFX+n] = rhoaux_new[n];
    }


    // density and momenta have already been updated

    // get the energy release from the change in rhoe over the timestep
    // excluding any advection, and use that to update the total energy

    Real rho_enucdot = (rhoe_new - rhoe_orig) / dt - state.ydot_a[SEINT];

    state.y[SEINT] = rhoe_new;
    state.y[SEDEN] += dt * state.ydot_a[SEDEN] + dt * rho_enucdot;

}

#else
///
/// update the state due to NSE changes for the simplified-SDC case
/// this version works with the self-consistent NSE
///
template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void sdc_nse_burn(BurnT& state, const Real dt) {

    state.success = true;
    state.n_rhs = 0;
    state.n_jac = 0;

    // store the initial mass fractions -- we will need these
    // to compute the energy release.

    Real X_old[NumSpec];

    for (int n = 0; n < NumSpec; ++n) {
        X_old[n] = state.y[SFS+n] / state.y[SRHO];
    }

    // density and momentum have no reactive sources
    Real rho_old = state.y[SRHO];

    state.y[SRHO] += dt * state.ydot_a[SRHO];
    state.y[SMX] += dt * state.ydot_a[SMX];
    state.y[SMY] += dt * state.ydot_a[SMY];
    state.y[SMZ] += dt * state.ydot_a[SMZ];

    // if we are doing drive_initial_convection, we want to use
    // the temperature that comes in through T_fixed

    Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;

    // get the neutrino loss term -- we want to use the state that we
    // came in here with, so the original Abar and Zbar
    Real snu{0.0};
    Real dsnudt{0.0};
    Real dsnudd{0.0};
    Real dsnuda{0.0};
    Real dsnudz{0.0};

    Real abar{0.0};
    Real zbar{0.0};
    for (int n = 0; n < NumSpec; ++n) {
        abar += X_old[n] * aion_inv[n];
        zbar += X_old[n] * zion[n] * aion_inv[n];
    }
    abar = 1.0 / abar;
    zbar *= abar;

    constexpr int do_derivatives = 0;
    sneut5<do_derivatives>(T_in, rho_old, abar, zbar,
                           snu, dsnudt, dsnudd, dsnuda, dsnudz);

    Real snu_old = snu;


    // if our network could return the evolution of Ye due to the
    // weak interactions, we would evaluate the NSE state here and
    // get dYe/dt.

    Real dyedt_old = 0.0;


    // predict the U^{n+1,*} state with only estimates of the X
    // updates with advection to dt and the neutrino loss term in
    // energy

    BurnT burn_state;
    burn_state.T = T_in; // initial guess
    burn_state.mu_p = state.mu_p;
    burn_state.mu_n = state.mu_n;

    Real rhoX_source[NumSpec] = {0.0_rt};

    Real rhoX_new[NumSpec];

    Real rhoe_new;
    Real rho_enucdot = -rho_old * snu;

    Real rho_half = 0.5_rt * (rho_old + state.y[SRHO]);

    // predict the new aux for the first iteration -- this is really
    // just including the advection bits

    for (int n = 0; n < NumSpec; n++) {
        rhoX_new[n] = state.y[SFS+n] + dt * state.ydot_a[SFS+n] + dt * rhoX_source[n];
    }

    burn_t nse_state;

    for (int iter = 0; iter < integrator_rp::nse_iters; iter++) {

        // update (rho e)^{n+1} based on the new energy generation rate
        rhoe_new = state.y[SEINT] + dt * state.ydot_a[SEINT] + dt * rho_enucdot;

        // call the EOS to get the updated T*

        Real T_new;
        burn_state.rho = state.y[SRHO];
        burn_state.e = rhoe_new / state.y[SRHO];
        for (int n = 0; n < NumSpec; n++) {
            burn_state.xn[n] = rhoX_new[n] / state.y[SRHO];
            burn_state.y[SFS+n] = rhoX_new[n];
        }

        if (state.T_fixed > 0) {
            T_new = state.T_fixed;
        } else {
            eos(eos_input_re, burn_state);
            T_new = burn_state.T;
        }
        burn_state.T = T_new;

        // solve for the NSE state for this network we need to call
        // get_nse_state with a burn_t.  We will have it compute
        // Ye from the input X's

        nse_state = get_actual_nse_state(burn_state);

        // compute the energy release.  The mass fractions in nse_state.xn[]
        // include the advective parts, so first we need to remove that.

        Real rhoX_tilde[NumSpec];
        for (int n = 0; n < NumSpec; ++n) {
            rhoX_tilde[n] = state.y[SRHO] * nse_state.xn[n] - dt * state.ydot_a[SFS+n];
        }

        Real dyedt = 0.0_rt;   // we can update this in the future by calling actual_rhs()

        // we want to compute (rho eps) = - N_A c^2 sum{m_i (rhoX_tilde - rhoX_old) / A_i}
        rho_enucdot = 0.0;
        for (int n = 0; n < NumSpec; ++n) {
            rho_enucdot += (rhoX_tilde[n] - rho_old * X_old[n]) *
                network::mion(n+1) * aion_inv[n];
        }
        rho_enucdot *= C::Legacy::enuc_conv2;

        // now get the updated neutrino term
        abar = 0.0;
        zbar = 0.0;
        for (int n = 0; n < NumSpec; ++n) {
            abar += nse_state.xn[n] * aion_inv[n];
            zbar += nse_state.xn[n] * zion[n] * aion_inv[n];
        }
        abar = 1.0 / abar;
        zbar *= abar;

        constexpr int do_derivatives = 0;
        sneut5<do_derivatives>(T_new, state.y[SRHO], abar, zbar,
                               snu, dsnudt, dsnudd, dsnuda, dsnudz);

        rho_enucdot -= 0.5_rt * rho_half * (snu_old + snu);

        // update the new state for the next pass -- this should
        // already implicitly have the advective portion included,
        // since it was there when we called the NSE state

        for (int n = 0; n < NumSpec; n++) {
            rhoX_new[n] = state.y[SRHO] * nse_state.xn[n];
        }

    }

    // now update the aux quantities

    // the new mass fractions are just those that come from the table
    // make sure they are normalized
    Real sum_X{0.0_rt};
    Real X[NumSpec] = {0.0_rt};
    for (int n = 0; n < NumSpec; ++n) {
        X[n] = amrex::max(small_x, amrex::min(1.0_rt, nse_state.xn[n]));
        sum_X += X[n];
    }

    for (int n = 0; n < NumSpec; ++n) {
        X[n] /= sum_X;
    }

    for (int n = 0; n < NumSpec; ++n) {
        state.y[SFS+n] = state.y[SRHO] * X[n];
    }


    // density and momenta have already been updated

    // update the total and internal energy now

    state.y[SEINT] += dt * state.ydot_a[SEINT] + dt * rho_enucdot;
    state.y[SEDEN] += dt * state.ydot_a[SEDEN] + dt * rho_enucdot;

    // store the chemical potentials
    state.mu_p = nse_state.mu_p;
    state.mu_n = nse_state.mu_n;

}
#endif

#endif
