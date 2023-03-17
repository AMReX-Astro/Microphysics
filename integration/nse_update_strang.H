#ifndef NSE_UPDATE_STRANG_H
#define NSE_UPDATE_STRANG_H

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
/// update the state due to NSE changes for the Strang-split case
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_burn(burn_t& state, const Real dt) {

#if defined(NSE_TABLE)

  using namespace AuxZero;

  // use the NSE table
  Real abar_out;
  Real dq_out;
  Real dyedt;
  Real X[NumSpec];
#endif

  // get the energy -- we are assuming that rho, T are valid on input

  eos(eos_input_rt, state);

  Real e_in = state.e;

#if defined(NSE_TABLE)
  Real bea_in = state.aux[ibea];
  Real ye_in = state.aux[iye];
#endif

  // if we are doing drive_initial_convection, we want to use
  // the temperature that comes in through T_fixed

  Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;

#if defined(NSE_TABLE)
  // call the NSE table using the * state to get the t^{n+1}
  // source estimates.  The thermodynamnics here is specified
  // in terms of the auxillary composition, Ye, abar, and B/A

  nse_interp(T_in, state.rho, state.aux[iye],
             abar_out, dq_out, dyedt, X);

  // break the update into multiple substeps of length dt / nse_substeps

  Real dt_sub = dt / static_cast<Real>(nse_substeps);

  for (int iter = 0; iter < nse_substeps; ++iter) {

      // update Ye

      state.aux[iye] += dt_sub * dyedt;

      // now get the composition from the table using the upated Ye

      Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;

      nse_interp(T_in, state.rho, state.aux[iye],
                 abar_out, dq_out, dyedt, X);

      // this is MeV / nucleon -- here aux has not yet been updated, so we
      // access the old binding energy
      Real deltaq = dq_out - state.aux[ibea];

      state.aux[ibea] += deltaq;
      state.aux[iabar] = abar_out;

      // convert the energy to erg / g
      Real enuc = deltaq * C::MeV2eV * C::ev2erg * C::n_A;

      state.e += enuc;

      // update the temperature for the next pass

      eos(eos_input_re, state);

  }

#elif defined(NSE_NET)

  // solve for the NSE state -- first compute Ye
  state.y_e = 0.0_rt;
  for (int n = 0; n < NumSpec; ++n) {
      state.y_e += zion[n] * state.xn[n] * aion_inv[n];
  }

  // TODO: need to make this respect T_fixed

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

  // convert the energy to erg / g
  Real enuc = deltaq * C::MeV2eV * C::ev2erg * C::n_A;

  state.e = e_in + enuc;

#endif

  state.success = true;
  state.n_rhs = 0;
  state.n_jac = 0;


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
