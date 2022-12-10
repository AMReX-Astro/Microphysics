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

using namespace amrex;

#ifdef SIMPLIFIED_SDC
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
    nse_interp(state.T, state.rho, ye_in,
               abar_out, dq_out, dyedt, X);

    Real dyedt_old = dyedt;

    // density and momentum have no reactive sources
    Real rho_old = state.y[SRHO];
    Real rho_bea_old = state.y[SFX+ibea];
    Real rho_abar_old = state.y[SFX+iabar];

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
    Real rho_enucdot = 0.0_rt;

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

        eos_state.rho = state.y[SRHO];
        eos_state.e = rhoe_new / state.y[SRHO];
        for (int n = 0; n < NumAux; n++) {
            eos_state.aux[n] = rho_aux_new[n] / state.y[SRHO];
        }

        eos(eos_input_re, eos_state);

        // call the NSE table using the * state to get the t^{n+1} aux / composition

        nse_interp(eos_state.T, eos_state.rho, eos_state.aux[iye],
                   abar_out, dq_out, dyedt, X);

        // note that the abar, dq (B/A), and X here have all already
        // seen advection implicitly

        // compute the energy release -- we need to remove the advective part

        Real rho_bea_tilde = state.y[SRHO] * dq_out - dt * state.ydot_a[SFX+ibea];
        Real rho_dBEA = rho_bea_tilde - rho_bea_old;  // this is MeV / nucleon * g / cm**3

        // convert the energy to erg / cm**3
        rho_enucdot  = rho_dBEA * C::MeV2eV * C::ev2erg * C::n_A / dt;

        // update the new state for the next pass

        aux_source[iye] = 0.5_rt * rho_half * (dyedt_old + dyedt);
        rho_aux_new[iye] = state.y[SFX+iye] + dt * state.ydot_a[SFX+iye] + dt * aux_source[iye];

        rho_aux_new[iabar] = state.y[SRHO] * abar_out;
        rho_aux_new[ibea] = state.y[SRHO] * dq_out;

    }

    // now update the aux quantities

    // the new mass fractions are just those that come from the table

    for (int n = 0; n < NumSpec; n++) {
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
/// update the state due to NSE changes for the Strang-split case
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_burn(burn_t& state, const Real dt) {

  using namespace AuxZero;

  // use the NSE table
  Real abar_out;
  Real dq_out;
  Real dyedt;
  Real X[NumSpec];

  // use the NSE table to get an estimate of dYe/dt

  nse_interp(state.T, state.rho, state.aux[iye],
             abar_out, dq_out, dyedt, X);

  // update Ye

  state.aux[iye] += dt * dyedt;

  // now get the composition from the table using the upated Ye

  nse_interp(state.T, state.rho, state.aux[iye],
             abar_out, dq_out, dyedt, X);

  state.success = true;
  state.n_rhs = 0;
  state.n_jac = 0;

  // this is MeV / nucleon -- here aux has not yet been updated, so we
  // access the old binding energy
  Real deltaq = dq_out - state.aux[ibea];

  // under-relaxation / inertia (see Ma et el. 2013)
  deltaq = eta * deltaq;

  state.aux[ibea] += deltaq;

  // convert the energy to erg / g
  Real enuc = deltaq * C::MeV2eV * C::ev2erg * C::n_A;

  state.e = enuc + state.e;

  state.aux[iabar] = abar_out;

  // store the new composition
  for (int n = 0; n < NumSpec; n++) {
    state.xn[n] = X[n];
  }

}
#endif

#endif
