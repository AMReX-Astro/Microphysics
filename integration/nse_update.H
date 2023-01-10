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

#if defined (NSE_TABLE)
    // get the current NSE state from the table

    nse_interp(state.T, state.rho, ye_in,
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

        eos_state.rho = state.y[SRHO];
        eos_state.e = rhoe_new / state.y[SRHO];
        for (int n = 0; n < NumAux; n++) {
            eos_state.aux[n] = rho_aux_new[n] / state.y[SRHO];
        }

        eos(eos_input_re, eos_state);

#if defined(NSE_TABLE)
        // call the NSE table using the * state to get the t^{n+1}
        // source estimates

        nse_interp(eos_state.T, eos_state.rho, eos_state.aux[iye],
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

#else
///
/// update the state due to NSE changes for the Strang-split case
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_burn(burn_t& state, const Real dt) {

  using namespace AuxZero;

#if defined(NSE_TABLE)

  // use the NSE table
  Real abar_out;
  Real dq_out;
  Real dyedt;
  Real X[NumSpec];

  // get the energy -- we are assuming that rho, T are valid on input

  eos(eos_input_rt, state);

  // call the NSE table using the * state to get the t^{n+1}
  // source estimates.  The thermodynamnics here is specified
  // in terms of the auxillary composition, Ye, abar, and B/A

  nse_interp(state.T, state.rho, state.aux[iye],
             abar_out, dq_out, dyedt, X);

  // update Ye

  state.aux[iye] += dt * dyedt;

  // now get the composition from the table using the upated Ye

  nse_interp(state.T, state.rho, state.aux[iye],
             abar_out, dq_out, dyedt, X);


  // this is MeV / nucleon -- here aux has not yet been updated, so we
  // access the old binding energy
  Real deltaq = dq_out - state.aux[ibea];

  // under-relaxation / inertia (see Ma et el. 2013)
  deltaq = eta * deltaq;

  state.aux[ibea] += deltaq;
  state.aux[iabar] = abar_out;

#elif defined(NSE_NET)

  // solve for the NSE state -- first compute Ye
  state.ye = 0.0_rt;
  for (int n = 0; n < NumSpec; ++n) {
      state.ye += zion[n] * state.xn[n] * aion_inv[n];
  }

  auto nse_state = get_actual_nse_state(state);

  // for now, we'll take dyedt = 0.  Later we can evaluate this by
  // calling the RHS with the NSE state.
  Real dyedt = 0.0_rt;

  // compute the change in binding energy -- this is the energy release
  // we want to do sum {B dY}, where Y is the molar fraction.
  // this will be in MeV / nucleon
  Real deltaq = 0.0_rt;
  for (int n = 0; n < NumSpec; ++n) {
      deltaq += network::bion(n+1) * aion_inv[n] * (nse_state.xn[n] - state.xn[n]);
  }

  // under-relaxation / inertia (see Ma et el. 2013)
  deltaq = eta * deltaq;

#endif

  state.success = true;
  state.n_rhs = 0;
  state.n_jac = 0;

  // convert the energy to erg / g
  Real enuc = deltaq * C::MeV2eV * C::ev2erg * C::n_A;

  state.e = enuc + state.e;


  // store the new composition

#if defined(NSE_TABLE)
  for (int n = 0; n < NumSpec; ++n) {
    state.xn[n] = X[n];
  }
#elif defined(NSE_NET)
  for (int n = 0; n < NumSpec; ++n) {
    state.xn[n] = nse_state.xn[n];
  }
#endif


}
#endif

#endif
