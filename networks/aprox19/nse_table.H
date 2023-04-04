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
  if (nse_table.fail()) {
      amrex::Error("unable to open nse19.tbl");
  }

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
Real cubic(const Real* xs, const Real* fs, const Real dx, const Real x) {

    // fit a cubic of the form
    // f(x) = a (x - x_i)**3 + b (x - x_i)**2 + c (x - x_i) + d
    // to the data (xs, fs)
    // we take x_i to be x[1]

    Real a = (3 * fs[1] - 3 * fs[2] + fs[3] - fs[0]) / (6 * std::pow(dx, 3));
    Real b = (-2 * fs[1] + fs[2] + fs[0]) / (2 * dx * dx);
    Real c = (-3 * fs[1] + 6 * fs[2] - fs[3] - 2 * fs[0]) / (6 * dx);
    Real d = fs[1];

    return a * std::pow(x - xs[1], 3) + b * std::pow(x - xs[1], 2) + c * (x - xs[1]) + d;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real trilinear(const int ir1, const int it1, const int ic1,
               const Real rho, const Real temp, const Real ye, const T& data) {

    // find the eight interpolation points in the 1D arrays

    int it1r1c1 = nse_idx(ir1,   it1,   ic1);
    int it1r1c2 = nse_idx(ir1,   it1,   ic1+1);
    int it1r2c1 = nse_idx(ir1+1, it1,   ic1);
    int it1r2c2 = nse_idx(ir1+1, it1,   ic1+1);
    int it2r1c1 = nse_idx(ir1,   it1+1, ic1);
    int it2r1c2 = nse_idx(ir1,   it1+1, ic1+1);
    int it2r2c1 = nse_idx(ir1+1, it1+1, ic1);
    int it2r2c2 = nse_idx(ir1+1, it1+1, ic1+1);

    Real t0 = nse_table_logT(it1);
    Real r0 = nse_table_logrho(ir1);
    Real x0 = nse_table_ye(ic1);

    Real td = (temp - t0) / table::dlogT;
    Real rd = (rho - r0) / table::dlogrho;
    Real xd = (x0 - ye) / table::dye;
    xd = amrex::max(0.0_rt, xd);

    Real omtd = 1.0_rt - td;
    Real omrd = 1.0_rt - rd;
    Real omxd = 1.0_rt - xd;

    Real val = data(it1r1c1) * omtd * omrd * omxd +
               data(it1r1c2) * omtd * omrd * xd +
               data(it1r2c1) * omtd * rd * omxd +
               data(it1r2c2) * omtd * rd * xd +
               data(it2r1c1) * td * omrd * omxd +
               data(it2r1c2) * td * omrd * xd +
               data(it2r2c1) * td * rd * omxd +
               data(it2r2c2) * td * rd * xd;

    return val;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real tricubic(const int ir0, const int it0, const int ic0,
              const Real rho, const Real temp, const Real ye, const T& data) {

    Real yes[] = {nse_table_ye(ic0),
                  nse_table_ye(ic0+1),
                  nse_table_ye(ic0+2),
                  nse_table_ye(ic0+3)};

    Real Ts[] = {nse_table_logT(it0),
                 nse_table_logT(it0+1),
                 nse_table_logT(it0+2),
                 nse_table_logT(it0+3)};

    Real rhos[] = {nse_table_logrho(ir0),
                   nse_table_logrho(ir0+1),
                   nse_table_logrho(ir0+2),
                   nse_table_logrho(ir0+3)};

    // first do the 16 ye interpolations

    // the first index will be rho and the second will be T
    Real d1[4][4];

    for (int ii = 0; ii < 4; ++ii) {
        for (int jj = 0; jj < 4; ++jj) {

            Real _d[] = {data(nse_idx(ir0+ii, it0+jj, ic0)),
                         data(nse_idx(ir0+ii, it0+jj, ic0+1)),
                         data(nse_idx(ir0+ii, it0+jj, ic0+2)),
                         data(nse_idx(ir0+ii, it0+jj, ic0+3))};

            // note that the ye values are monotonically decreasing,
            // so the "dx" needs to be negative
            d1[ii][jj] = cubic(yes, _d, -table::dye, ye);
        }
    }

    // now do the 4 T interpolations (one in each rho plane)

    Real d2[4];

    for (int ii = 0; ii < 4; ++ii) {

        Real _d[] = {d1[ii][0], d1[ii][1], d1[ii][2], d1[ii][3]};
        d2[ii] = cubic(Ts, _d, table::dlogT, temp);
    }

    // finally do the remaining interpolation over rho

    Real val = cubic(rhos, d2, table::dlogrho, rho);

    return val;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_interp(const Real T, const Real rho, const Real ye,
                Real& abar, Real& dq, Real& dyedt, Real* X, bool skip_X_fill=false) {

    // if skip_X_fill = true then we don't fill X[] with the mass fractions.

    using namespace table;
    using namespace AuxZero;

    Real rholog = amrex::max(7.0_rt, amrex::min(10.0_rt, std::log10(rho)));
    Real tlog = amrex::max(9.0_rt, amrex::min(10.4_rt, std::log10(T)));
    Real yet = amrex::max(0.40_rt, amrex::min(0.50_rt, ye));

    if (nse_table_interp_linear) {

        int ir1 = nse_get_logrho_index(rholog);
        int it1 = nse_get_logT_index(tlog);
        int ic1 = nse_get_ye_index(yet);

        abar = trilinear(ir1, it1, ic1, rholog, tlog, yet, abartab);
        dq = trilinear(ir1, it1, ic1, rholog, tlog, yet, ebtab);
        dyedt = trilinear(ir1, it1, ic1, rholog, tlog, yet, wratetab);

        // this is actually the sum of all e- capture and e+ decay, so if
        // e- capture dominates, this quantity is positive, but Ye should
        // decrease, so we swap the sign here.

        dyedt = -dyedt;

        // massfractab is 2-d, so we wrap the access in a lambda already
        // indexing the component

        if (! skip_X_fill) {
            for (int n = 1; n <= NumSpec; n++) {
                Real _X = trilinear(ir1, it1, ic1, rholog, tlog, yet,
                                    [=] (const int i) {return massfractab(n, i);});
                X[n-1] = amrex::max(0.0_rt, amrex::min(1.0_rt, _X));
            }
        }

    } else {

        // for a cubic interpolant, we need 4 points that span the data value
        // for temperature, these will be it0, it0+1, it0+2, it0+3
        // with the idea that the temperature we want is between it0+1 and it0+2
        // so we offset one to the left and also ensure that we don't go off the table

        int ir0 = nse_get_logrho_index(rholog) - 1;
        ir0 = amrex::max(1, amrex::min(table::nden-3, ir0));

        int it0 = nse_get_logT_index(tlog) - 1;
        it0 = amrex::max(1, amrex::min(table::ntemp-3, it0));

        int ic0 = nse_get_ye_index(yet) - 1;
        ic0 = amrex::max(1, amrex::min(table::nye-3, ic0));

        abar = tricubic(ir0, it0, ic0, rholog, tlog, yet, abartab);
        dq = tricubic(ir0, it0, ic0, rholog, tlog, yet, ebtab);
        dyedt = tricubic(ir0, it0, ic0, rholog, tlog, yet, wratetab);

        // this is actually the sum of all e- capture and e+ decay, so if
        // e- capture dominates, this quantity is positive, but Ye should
        // decrease, so we swap the sign here.

        dyedt = -dyedt;

        // massfractab is 2-d, so we wrap the access in a lambda already
        // indexing the component

        if (! skip_X_fill) {
            for (int n = 1; n <= NumSpec; n++) {
                Real _X = tricubic(ir0, it0, ic0, rholog, tlog, yet,
                                   [=] (const int i) {return massfractab(n, i);});
                X[n-1] = amrex::max(0.0_rt, amrex::min(1.0_rt, _X));
            }
        }
    }

}

#endif
