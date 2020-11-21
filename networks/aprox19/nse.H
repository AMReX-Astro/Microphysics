#ifndef NSE_H
#define NSE_H

#include <iostream>
#include <fstream>
#include <actual_network.H>

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include "extern_parameters.H"

#include "burn_type.H"

using namespace amrex;

AMREX_INLINE
void init_nse() {

  // set table parameters

  // read in table
  std::ifstream nse_table;

  std::cout << "reading the NSE table (C++) ..." << std::endl;

  nse_table.open("nse19.tbl", std::ios::in);

  Real the, tsi, tfe;

  for (int irho = 1; irho <= table::nden; irho++) {
    for (int it9 = 1; it9 <= table::ntemp; it9++) {
      for (int iye = 1; iye <= table::nye; iye++) {
        int j = (irho-1)*table::ntemp*table::nye + (it9-1)*table::nye + iye;

        nse_table >> table::ttlog(j) >> table::ddlog(j) >> table::yetab(j);
        nse_table >> the >> tsi >> tfe;
        nse_table >> table::abartab(j) >> table::ebtab(j) >> table::wratetab(j);
        for (int n = 1; n <= NumSpec; n++) {
          nse_table >> table::massfractab(n, j);
        }
      }
    }
  }

}


template <class state_t> 
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool in_nse(const state_t& state) {

  using namespace Species;

  bool nse_check = false;

  if (state.rho > rho_nse && state.T > T_nse) {

    // Ma checks on Fe-group (for our composition, this means Cr48, Fe52, Fe54, Ni56)
    // and C-group (for us, that is C12, N14)
    // and He-group (for us, that is H1, He3, He4)

    Real Fe_group = state.xn[Cr48-1] + state.xn[Fe52-1] +
                    state.xn[Fe54-1] + state.xn[Ni56-1];
    Real C_group = state.xn[C12-1] + state.xn[N14-1];
    Real He_group = state.xn[H1-1] + state.xn[He3-1] + state.xn[He4-1];

    if (Fe_group + He_group > He_Fe_nse && C_group < C_nse) {
      nse_check = true;
    }
  }

  return nse_check;

}

///
/// After doing a burn via integrating the network, update the composition
/// aux variables to be consistent with the new X's resulting from the burn.
///
template <class state_t> 
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void set_nse_aux_from_X(state_t& state) {

  using namespace AuxZero;

  state.aux[iye] = 0.0_rt;
  state.aux[iabar] = 0.0_rt;
  state.aux[ibea] = 0.0_rt;

  for (int n = 0; n < NumSpec; n++) {
    state.aux[iye] += state.xn[n] * zion[n] * aion_inv[n];
    state.aux[iabar] += state.xn[n] * aion_inv[n];
    state.aux[ibea] += state.xn[n] * aprox19::bion(n+1) * aion_inv[n];
  }
  state.aux[iabar] = 1.0_rt/state.aux[iabar];

}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_interp(const Real T, const Real rho, const Real ye,
                Real& abar, Real& dq, Real& dyedt, Real* X) {

  using namespace table;
  using namespace AuxZero;

  Real tlog = std::log10(T);
  Real rholog = std::log10(rho);
  Real yet = ye;

  if (tlog < 9.0_rt) {
    tlog = 9.0_rt;
  }

  if (tlog > 10.4_rt) {
    tlog = 10.4_rt;
  }

  int it1 = static_cast<int>((tlog - 9.0_rt)*50.0_rt - 1.e-6_rt);
  it1 += 1;
  int it2 = it1 + 1;

  if (rholog < 7.0_rt) {
    rholog = 7.0_rt;
  }

  if (rholog > 10.0_rt) {
    rholog = 10.0_rt;
  }

  int ir1 = static_cast<int>((rholog - 7.0_rt)*10.0_rt - 1.e-6_rt);
  ir1 += 1;
  int ir2 = ir1 + 1;

  if (yet < 0.40_rt) {
    yet = 0.40_rt;
  }

  if (yet > 0.50_rt) {
    yet = 0.50_rt;
  }

  int ic1 = static_cast<int>((0.50_rt - yet)/0.005_rt - 1.0e-6_rt);
  ic1 += 1;
  int ic2 = ic1 + 1;

  // find the eight interpolation points in the 1D arrays

  int it1r1c1 = (ir1-1)*71*21 + (it1-1)*21 + ic1;
  int it1r1c2 = (ir1-1)*71*21 + (it1-1)*21 + ic2;
  int it1r2c1 = (ir2-1)*71*21 + (it1-1)*21 + ic1;
  int it1r2c2 = (ir2-1)*71*21 + (it1-1)*21 + ic2;
  int it2r1c1 = (ir1-1)*71*21 + (it2-1)*21 + ic1;
  int it2r1c2 = (ir1-1)*71*21 + (it2-1)*21 + ic2;
  int it2r2c1 = (ir2-1)*71*21 + (it2-1)*21 + ic1;
  int it2r2c2 = (ir2-1)*71*21 + (it2-1)*21 + ic2;

  Real t0 = 9.0_rt + static_cast<Real>(it1-1)*0.02_rt;
  Real r0 = 7.0_rt + static_cast<Real>(ir1-1)*0.10_rt;
  Real x0 = 0.50_rt - static_cast<Real>(ic1-1)*0.005_rt;

  Real td = (tlog - t0)/0.02_rt;
  Real rd = (rholog - r0)/0.10_rt;
  Real xd = (x0-yet)/0.005_rt;
  xd = amrex::max(0.0_rt, xd);

  Real omtd = 1.0_rt - td;
  Real omrd = 1.0_rt - rd;
  Real omxd = 1.0_rt - xd;

  abar =
    abartab(it1r1c1)*omtd*omrd*omxd +
    abartab(it1r1c2)*omtd*omrd*xd +
    abartab(it1r2c1)*omtd*rd*omxd +
    abartab(it1r2c2)*omtd*rd*xd +
    abartab(it2r1c1)*td*omrd*omxd +
    abartab(it2r1c2)*td*omrd*xd +
    abartab(it2r2c1)*td*rd*omxd +
    abartab(it2r2c2)*td*rd*xd;

  dq =
    ebtab(it1r1c1)*omtd*omrd*omxd +
    ebtab(it1r1c2)*omtd*omrd*xd +
    ebtab(it1r2c1)*omtd*rd*omxd +
    ebtab(it1r2c2)*omtd*rd*xd +
    ebtab(it2r1c1)*td*omrd*omxd +
    ebtab(it2r1c2)*td*omrd*xd +
    ebtab(it2r2c1)*td*rd*omxd +
    ebtab(it2r2c2)*td*rd*xd;

  dyedt =
    wratetab(it1r1c1)*omtd*omrd*omxd +
    wratetab(it1r1c2)*omtd*omrd*xd +
    wratetab(it1r2c1)*omtd*rd*omxd +
    wratetab(it1r2c2)*omtd*rd*xd +
    wratetab(it2r1c1)*td*omrd*omxd +
    wratetab(it2r1c2)*td*omrd*xd +
    wratetab(it2r2c1)*td*rd*omxd +
    wratetab(it2r2c2)*td*rd*xd;

  // this is actually the sum of all e- capture and e+ decay, so if
  // e- capture dominates, this quantity is positive, but Ye should
  // decrease, so we swap the sign here.
  dyedt = -dyedt;

  for (int n = 1; n <= NumSpec; n++) {
    X[n-1] =
      massfractab(n, it1r1c1)*omtd*omrd*omxd +
      massfractab(n, it1r1c2)*omtd*omrd*xd +
      massfractab(n, it1r2c1)*omtd*rd*omxd +
      massfractab(n, it1r2c2)*omtd*rd*xd +
      massfractab(n, it2r1c1)*td*omrd*omxd +
      massfractab(n, it2r1c2)*td*omrd*xd +
      massfractab(n, it2r2c1)*td*rd*omxd +
      massfractab(n, it2r2c2)*td*rd*xd;
  }

}


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
