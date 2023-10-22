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
    Real ye_in = state.y[SFX+iye] / state.rho;
    Real rho_bea_old = state.y[SFX+ibea];

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

    Real abar = state.y[SFX+iabar] / rho_old;
    Real zbar = abar * ye_in;

    constexpr int do_derivatives = 0;
    sneut5<do_derivatives>(T_in, rho_old, abar, zbar,
                           snu, dsnudt, dsnudd, dsnuda, dsnudz);

    Real snu_old = snu;

    // get the current NSE state from the table -- this will be used
    // to compute the NSE evolution sources

    Real abar_out;
    Real dq_out;
    Real dyedt;
    Real X[NumSpec];

    nse_interp(T_in, rho_old, ye_in,
               abar_out, dq_out, dyedt, X);

    Real dyedt_old = dyedt;


    // predict the U^{n+1,*} state with only estimates of the aux
    // reaction terms and advection to dt

    eos_re_t eos_state;
    eos_state.T = T_in; // initial guess

    // initial aux_sources
    Real aux_source[NumAux] = {0.0_rt};
    aux_source[iye] = rho_old * dyedt;

    Real rho_aux_new[NumAux];

    Real rhoe_new;
    Real rho_enucdot = -rho_old * snu;

    Real rho_half = 0.5_rt * (rho_old + state.y[SRHO]);

    // predict the new aux for the first iteration -- this is really
    // just including the advection bits

    for (int n = 0; n < NumAux; n++) {
        rho_aux_new[n] = state.y[SFX+n] + dt * state.ydot_a[SFX+n] + dt * aux_source[n];
    }

    for (int iter = 0; iter < 3; iter++) {

        // update (rho e)^{n+1} based on the new energy generation rate
        rhoe_new = state.y[SEINT] + dt * state.ydot_a[SEINT] + dt * rho_enucdot;

        // call the EOS to get the updated T*

        Real T_new;
        eos_state.rho = state.y[SRHO];
        eos_state.e = rhoe_new / state.y[SRHO];
        for (int n = 0; n < NumAux; n++) {
            eos_state.aux[n] = rho_aux_new[n] / state.y[SRHO];
        }

        if (state.T_fixed > 0) {
            T_new = state.T_fixed;
        } else {
            eos(eos_input_re, eos_state);
            T_new = eos_state.T;
        }

        // call the NSE table using the * state to get the t^{n+1}
        // source estimates

        nse_interp(T_new, eos_state.rho, eos_state.aux[iye],
                   abar_out, dq_out, dyedt, X);

        // note that the abar, dq (B/A), and X here have all already
        // seen advection implicitly


        // compute the energy release -- we need to remove the advective part

        Real rho_bea_tilde = state.y[SRHO] * dq_out - dt * state.ydot_a[SFX+ibea];
        Real rho_dBEA = rho_bea_tilde - rho_bea_old;  // this is MeV / nucleon * g / cm**3

        // convert the energy to erg / cm**3
        rho_enucdot = rho_dBEA * C::MeV2eV * C::ev2erg * C::n_A / dt;

        // now get the updated neutrino term
        zbar = abar_out * eos_state.aux[iye];
        sneut5<do_derivatives>(T_new, eos_state.rho, abar_out, zbar,
                               snu, dsnudt, dsnudd, dsnuda, dsnudz);

        rho_enucdot -= 0.5_rt * rho_half * (snu_old + snu);

        // update the new state for the next pass

        aux_source[iye] = 0.5_rt * rho_half * (dyedt_old + dyedt);
        rho_aux_new[iye] = state.y[SFX+iye] + dt * state.ydot_a[SFX+iye] + dt * aux_source[iye];

        rho_aux_new[iabar] = state.y[SRHO] * abar_out;
        rho_aux_new[ibea] = state.y[SRHO] * dq_out;

    }

    // now update the aux quantities

    // the new mass fractions are just those that come from the table
    // make sure they are normalized
    Real sum_X{0.0_rt};
    for (int n = 0; n < NumSpec; ++n) {
        X[n] = amrex::max(small_x, amrex::min(1.0_rt, X[n]));
        sum_X += X[n];
    }

    for (int n = 0; n < NumSpec; ++n) {
        X[n] /= sum_X;
    }

    for (int n = 0; n < NumSpec; ++n) {
        state.y[SFS+n] = state.y[SRHO] * X[n];
    }

    // aux data comes from the table (except Ye, which we predicted)

    state.y[SFX+iye] = eos_state.rho * eos_state.aux[iye];
    state.y[SFX+iabar] = eos_state.rho * abar_out;
    state.y[SFX+ibea] = eos_state.rho * dq_out;


    // density and momenta have already been updated

    // update the total and internal energy now

    state.y[SEINT] += dt * state.ydot_a[SEINT] + dt * rho_enucdot;
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

    for (int iter = 0; iter < 3; iter++) {

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
