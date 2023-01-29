#ifndef NSE_TABLE_H
#define NSE_TABLE_H

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

AMREX_INLINE
void init_nse() {

  // set table parameters

  // read in table
  std::ifstream nse_table;

  amrex::Print() << "reading the NSE table (C++) ..." << std::endl;

  nse_table.open("nse19.tbl", std::ios::in);

  Real the, tsi, tfe;
  Real ttemp, tdens, tye;

  for (int irho = 1; irho <= table::nden; irho++) {
    for (int it9 = 1; it9 <= table::ntemp; it9++) {
      for (int iye = 1; iye <= table::nye; iye++) {
        int j = (irho-1)*table::ntemp*table::nye + (it9-1)*table::nye + iye;

        nse_table >> ttemp >> tdens >> tye;
        nse_table >> the >> tsi >> tfe;
        nse_table >> table::abartab(j) >> table::ebtab(j) >> table::wratetab(j);
        for (int n = 1; n <= NumSpec; n++) {
          nse_table >> table::massfractab(n, j);
        }
      }
    }
  }

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
int nse_idx(const int ir, const int it, const int ic) {
    // this uses a 1-based indexing
    return (ir-1) * table::ntemp * table::nye + (it-1) * table::nye + ic;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real nse_table_logT(const int it) {
    return 9.0_rt + static_cast<Real>(it-1) * table::dlogT;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real nse_table_logrho(const int ir) {
    return 7.0_rt + static_cast<Real>(ir-1) * table::dlogrho;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real nse_table_ye(const int ic) {
    return 0.50_rt - static_cast<Real>(ic-1) * table::dye;
}

// return the index in the table such that logrho[irho] < input density
// note: this is a 1-based index
AMREX_GPU_HOST_DEVICE AMREX_INLINE
int nse_get_logrho_index(const Real logrho) {

    int ir0 = static_cast<int>((logrho - 7.0_rt) / table::dlogrho - 1.e-6_rt);
    return ir0 + 1;
}

// return the index in the table such that logT[it] < input temperature
// note: this is a 1-based index
AMREX_GPU_HOST_DEVICE AMREX_INLINE
int nse_get_logT_index(const Real logT) {

    int it0 = static_cast<int>((logT - 9.0_rt) / table::dlogT - 1.e-6_rt);
    return it0 + 1;
}

// return the index in the table such that ye[ic] < input Ye
// note: this is a 1-based index
AMREX_GPU_HOST_DEVICE AMREX_INLINE
int nse_get_ye_index(const Real ye) {

    int ic0 = static_cast<int>((0.50_rt - ye) / table::dye - 1.0e-6_rt);
    return ic0 + 1;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_interp(const Real T, const Real rho, const Real ye,
                Real& abar, Real& dq, Real& dyedt, Real* X) {

  using namespace table;
  using namespace AuxZero;

  Real tlog = std::log10(T);
  Real rholog = std::log10(rho);
  Real yet = ye;

  Real rholog = amrex::max(7.0_rt, amrex::min(10.0_rt, std::log10(rho)));
  Real tlog = amrex::max(9.0_rt, amrex::min(10.4_rt, std::log10(T)));
  Real yet = amrex::max(0.40_rt, amrex::min(0.50_rt, ye));

  int ir1 = nse_get_logrho_index(rholog);
  int ir2 = ir1 + 1;

  int it1 = nse_get_logT_index(tlog);
  int it2 = it1 + 1;

  int ic1 = nse_get_ye_index(yet);
  int ic2 = ic1 + 1;

  // find the eight interpolation points in the 1D arrays

  int it1r1c1 = nse_idx(ir1, it1, ic1);
  int it1r1c2 = nse_idx(ir1, it1, ic2);
  int it1r2c1 = nse_idx(ir2, it1, ic1);
  int it1r2c2 = nse_idx(ir2, it1, ic2);
  int it2r1c1 = nse_idx(ir1, it2, ic1);
  int it2r1c2 = nse_idx(ir1, it2, ic2);
  int it2r2c1 = nse_idx(ir2, it2, ic1);
  int it2r2c2 = nse_idx(ir2, it2, ic2);

  Real t0 = nse_table_logT(it1);
  Real r0 = nse_table_logrho(ir1);
  Real x0 = nse_table_ye(ic1);

  Real td = (tlog - t0) / table::dlogT;
  Real rd = (rholog - r0) / table::dlogrho;
  Real xd = (x0-yet) / table::dye;
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

#endif
