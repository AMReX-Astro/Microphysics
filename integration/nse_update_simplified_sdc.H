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

///
/// update the state due to NSE changes for the simplified-SDC case
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void sdc_nse_burn(burn_t& state, const Real dt) {

    using namespace AuxZero;

    state.success = true;
    state.n_rhs = 0;
    state.n_jac = 0;

    // call the NSE table to get (dYe/dt)^n
    Real abar_out;
    Real dq_out;
    Real dyedt;
    Real X[NumSpec];

    Real ye_in = state.y[SFX+iye] / state.rho;

    // if we are doing drive_initial_convection, we want to use
    // the temperature that comes in through T_fixed

    Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;

#if defined (NSE_TABLE)
    // get the current NSE state from the table

    nse_interp(T_in, state.rho, ye_in,
               abar_out, dq_out, dyedt, X);

#elif defined (NSE_NET)
    // solve for the NSE state directly

    // compute abar, dq_out (B/A), and dyedt

#endif

    Real dyedt_old = dyedt;

    // density and momentum have no reactive sources
    Real rho_old = state.y[SRHO];
    Real rho_bea_old = state.y[SFX+ibea];

    state.y[SRHO] += dt * state.ydot_a[SRHO];
    state.y[SMX] += dt * state.ydot_a[SMX];
    state.y[SMY] += dt * state.ydot_a[SMY];
    state.y[SMZ] += dt * state.ydot_a[SMZ];

    // predict the U^{n+1,*} state with only estimates of the aux
    // reaction terms and advection to dt

    eos_re_t eos_state;
    eos_state.T = 1.e6_rt; // initial guess

    // initial aux_sources
    Real aux_source[NumAux] = {0.0_rt};

    aux_source[iye] = rho_old * dyedt;

    Real rhoe_new;
    Real rho_aux_new[NumAux];
    Real dq;
    Real enucdot = 0.0_rt;

    Real rho_half = 0.5_rt * (rho_old + state.y[SRHO]);

    for (int iter = 0; iter < 3; iter++) {

        rhoe_new = state.y[SEINT] + dt * state.ydot_a[SEINT] + dt * rho_half * enucdot;

        for (int n = 0; n < NumAux; n++) {
            rho_aux_new[n] = state.y[SFX+n] + dt * state.ydot_a[SFX+n] + dt * aux_source[n];
        }

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

#if defined(NSE_TABLE)
        // call the NSE table using the * state to get the t^{n+1}
        // source estimates

        nse_interp(T_new, eos_state.rho, eos_state.aux[iye],
                   abar_out, dq_out, dyedt, X);

#elif defined(NSE_NET)
        // solve for the NSE state for this network we need to call
        // get_nse_state with a burn_t, so we need to ensure that the
        // xn are in sync

        // compute abar, dq_out, and dyedt

#endif

        // compute the energy release

        // this is MeV / nucleon * g / cm**3
        dq = dq_out - rho_bea_old/rho_old;

        // convert the energy to erg / cm**3
        enucdot = dq * C::MeV2eV * C::ev2erg * C::n_A / dt;

        // update the aux source terms

        aux_source[iye] = 0.5_rt * rho_half * (dyedt_old + dyedt);
        aux_source[iabar] = rho_half * (abar_out - state.y[SFX+iabar] / rho_old) / dt;
        aux_source[ibea] = 0.0;   // this is actually not used in the EOS

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

    state.y[SEINT] += dt * state.ydot_a[SEINT] + dt * rho_half * enucdot;
    state.y[SEDEN] += dt * state.ydot_a[SEDEN] + dt * rho_half * enucdot;

}

#endif
