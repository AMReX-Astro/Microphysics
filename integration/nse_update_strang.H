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
#include <nse_table_type.H>
#include <nse_table.H>
#endif
#ifdef NSE_NET
#include <nse_solver.H>
#endif

using namespace amrex;

///
/// update the state due to NSE changes for the Strang-split case
///
template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_burn(BurnT& state, const Real dt) {

#if defined(NSE_TABLE)

  using namespace AuxZero;

  // use the NSE table

  // get the energy -- we are assuming that rho, T are valid on input

  eos(eos_input_rt, state);

  // if we are doing drive_initial_convection, we want to use
  // the temperature that comes in through T_fixed

  Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;

  // call the NSE table using the * state to get the t^{n+1}
  // source estimates.  The thermodynamnics here is specified
  // in terms of the auxiliary composition, Ye, abar, and B/A

  nse_table_t nse_state;
  nse_state.T = T_in;
  nse_state.rho = state.rho;
  nse_state.Ye = state.aux[iye];

  nse_interp(nse_state);

  // update Ye

  state.aux[iye] += dt * nse_state.dyedt;

  // now get the composition from the table using the updated Ye

  nse_state.Ye = state.aux[iye];
  nse_interp(nse_state);


  // this is MeV / nucleon -- here aux has not yet been updated, so we
  // access the old binding energy
  Real deltaq = nse_state.bea - state.aux[ibea];

  state.aux[ibea] += deltaq;
  state.aux[iabar] = nse_state.abar;

#elif defined(NSE_NET)

  // solve for the NSE state -- first compute Ye
  state.y_e = 0.0_rt;
  for (int n = 0; n < NumSpec; ++n) {
      state.y_e += zion[n] * state.xn[n] * aion_inv[n];
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
    state.xn[n] = nse_state.X[n];
  }
#elif defined(NSE_NET)
  for (int n = 0; n < NumSpec; ++n) {
    state.xn[n] = nse_state.xn[n];
  }
#endif


}

#endif
