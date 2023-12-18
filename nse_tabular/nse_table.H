#ifndef NSE_TABLE_H
#define NSE_TABLE_H

#include <iostream>
#include <fstream>

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include <extern_parameters.H>
#include <nse_table_data.H>
#include <nse_table_size.H>
#include <nse_table_type.H>


using namespace amrex;
using namespace network_rp;

AMREX_INLINE
void init_nse() {

  // set table parameters

  // read in table
  std::ifstream nse_table_file;

  amrex::Print() << "reading the NSE table (C++) ..." << std::endl;

  nse_table_file.open(nse_table_size::table_name, std::ios::in);
  if (nse_table_file.fail()) {
      amrex::Error("unable to open NSE table: " + nse_table_size::table_name);
  }

  Real ttemp, tdens, tye;

  // skip the header -- it is 4 lines
  std::string line;
  std::getline(nse_table_file, line);
  std::getline(nse_table_file, line);
  std::getline(nse_table_file, line);
  std::getline(nse_table_file, line);

  for (int irho = 1; irho <= nse_table_size::nden; irho++) {
      for (int it = 1; it <= nse_table_size::ntemp; it++) {
          for (int iye = 1; iye <= nse_table_size::nye; iye++) {
              int j = (irho-1) * nse_table_size::ntemp * nse_table_size::nye +
                  (it-1) * nse_table_size::nye + iye;

              std::getline(nse_table_file, line);
              if (line.empty()) {
                  amrex::Error("Error reading from the NSE table");
              }
              std::istringstream data(line);
              data >> ttemp >> tdens >> tye;
              data >> nse_table::abartab(j)
                   >> nse_table::beatab(j)
                   >> nse_table::dyedttab(j)
                   >> nse_table::dabardttab(j)
                   >> nse_table::dbeadttab(j)
                   >> nse_table::enutab(j);
              for (int n = 1; n <= NumSpec; n++) {
                  data >> nse_table::massfractab(n, j);
              }
          }
      }
  }

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
int nse_idx(const int ir, const int it, const int ic) {
    // this uses a 1-based indexing
    return (ir-1) * nse_table_size::ntemp * nse_table_size::nye + (it-1) * nse_table_size::nye + ic;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real nse_table_logT(const int it) {
    return nse_table_size::logT_min + static_cast<Real>(it-1) * nse_table_size::dlogT;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real nse_table_logrho(const int ir) {
    return nse_table_size::logrho_min + static_cast<Real>(ir-1) * nse_table_size::dlogrho;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real nse_table_ye(const int ic) {
    return nse_table_size::ye_max - static_cast<Real>(ic-1) * nse_table_size::dye;
}

// return the index in the table such that logrho[irho] < input density
// note: this is a 1-based index
AMREX_GPU_HOST_DEVICE AMREX_INLINE
int nse_get_logrho_index(const Real logrho) {

    int ir0 = static_cast<int>((logrho - nse_table_size::logrho_min) /
                               nse_table_size::dlogrho - 1.e-6_rt);
    return ir0 + 1;
}

// return the index in the table such that logT[it] < input temperature
// note: this is a 1-based index
AMREX_GPU_HOST_DEVICE AMREX_INLINE
int nse_get_logT_index(const Real logT) {

    int it0 = static_cast<int>((logT - nse_table_size::logT_min) /
                               nse_table_size::dlogT - 1.e-6_rt);
    return it0 + 1;
}

// return the index in the table such that ye[ic] < input Ye
// note: this is a 1-based index
AMREX_GPU_HOST_DEVICE AMREX_INLINE
int nse_get_ye_index(const Real ye) {

    int ic0 = static_cast<int>((nse_table_size::ye_max - ye) /
                               nse_table_size::dye - 1.0e-6_rt);
    return ic0 + 1;
}

///
/// given 4 points (xs, fs), with spacing dx, return the interplated
/// value of f at point x by fitting a cubic to the points
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real cubic(const Real* xs, const Real* fs, const Real dx, const Real x) {

    // fit a cubic of the form
    // f(x) = a (x - x_i)**3 + b (x - x_i)**2 + c (x - x_i) + d
    // to the data (xs, fs)
    // we take x_i to be x[1]

    Real a = (3 * fs[1] - 3 * fs[2] + fs[3] - fs[0]) / (6 * amrex::Math::powi<3>(dx));
    Real b = (-2 * fs[1] + fs[2] + fs[0]) / (2 * dx * dx);
    Real c = (-3 * fs[1] + 6 * fs[2] - fs[3] - 2 * fs[0]) / (6 * dx);
    Real d = fs[1];

    return a * amrex::Math::powi<3>(x - xs[1]) +
           b * amrex::Math::powi<2>(x - xs[1]) + c * (x - xs[1]) + d;

}

///
/// given 4 points (xs, fs), with spacing dx between the xs, return
/// the derivative of f at point x by fitting a cubic to the
/// points and differentiating the interpolant
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real cubic_deriv(const Real* xs, const Real* fs, const Real dx, const Real x) {

    // fit a cubic of the form
    // f(x) = a (x - x_i)**3 + b (x - x_i)**2 + c (x - x_i) + d
    // to the data (xs, fs)
    // we take x_i to be x[1]
    // then return dfdx = 3 a (x - x_i)**2 + 2 b (x - x_i) + c

    Real a = (3 * fs[1] - 3 * fs[2] + fs[3] - fs[0]) / (6 * amrex::Math::powi<3>(dx));
    Real b = (-2 * fs[1] + fs[2] + fs[0]) / (2 * dx * dx);
    Real c = (-3 * fs[1] + 6 * fs[2] - fs[3] - 2 * fs[0]) / (6 * dx);
    //Real d = fs[1];

    return 3.0_rt * a * amrex::Math::powi<2>(x - xs[1]) + 2.0_rt * b * (x - xs[1]) + c;

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

    Real td = (temp - t0) / nse_table_size::dlogT;
    Real rd = (rho - r0) / nse_table_size::dlogrho;
    Real xd = (x0 - ye) / nse_table_size::dye;
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
            d1[ii][jj] = cubic(yes, _d, -nse_table_size::dye, ye);
        }
    }

    // now do the 4 T interpolations (one in each rho plane)

    Real d2[4];

    for (int ii = 0; ii < 4; ++ii) {

        Real _d[] = {d1[ii][0], d1[ii][1], d1[ii][2], d1[ii][3]};
        d2[ii] = cubic(Ts, _d, nse_table_size::dlogT, temp);
    }

    // finally do the remaining interpolation over rho

    Real val = cubic(rhos, d2, nse_table_size::dlogrho, rho);

    return val;

}

///
/// take the temperature derivative of a table quantity by differentiating
/// the cubic interpolant
///
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real tricubic_dT(const int ir0, const int it0, const int ic0,
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
            d1[ii][jj] = cubic(yes, _d, -nse_table_size::dye, ye);
        }
    }

    // now do the 4 rho interpolations (one in each T plane)

    Real d2[4];

    for (int jj = 0; jj < 4; ++jj) {

        Real _d[] = {d1[0][jj], d1[1][jj], d1[2][jj], d1[3][jj]};
        d2[jj] = cubic(rhos, _d, nse_table_size::dlogrho, rho);
    }

    // finally do the remaining interpolation over T, but return
    // the derivative of the interpolant

    Real val = cubic_deriv(Ts, d2, nse_table_size::dlogT, temp);

    return val;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_interp(nse_table_t& nse_state, bool skip_X_fill=false) {

    // if skip_X_fill = true then we don't fill X[] with the mass fractions.

    using namespace nse_table;
    using namespace AuxZero;

    Real rholog = std::log10(nse_state.rho);
    {
        Real rmin = nse_table_size::logrho_min;
        Real rmax = nse_table_size::logrho_max;

        rholog = std::clamp(rholog, rmin, rmax);
    }

    Real tlog = std::log10(nse_state.T);
    {
        Real tmin = nse_table_size::logT_min;
        Real tmax = nse_table_size::logT_max;

        tlog = std::clamp(tlog, tmin, tmax);
    }

    Real yet = nse_state.Ye;
    {
        Real yemin = nse_table_size::ye_min;
        Real yemax = nse_table_size::ye_max;

        yet = std::clamp(yet, yemin, yemax);
    }

    if (nse_table_interp_linear) {

        int ir1 = nse_get_logrho_index(rholog);
        int it1 = nse_get_logT_index(tlog);
        int ic1 = nse_get_ye_index(yet);

        nse_state.abar = trilinear(ir1, it1, ic1, rholog, tlog, yet, abartab);
        nse_state.bea = trilinear(ir1, it1, ic1, rholog, tlog, yet, beatab);
        nse_state.dyedt = trilinear(ir1, it1, ic1, rholog, tlog, yet, dyedttab);
        nse_state.dabardt = trilinear(ir1, it1, ic1, rholog, tlog, yet, dabardttab);
        nse_state.dbeadt = trilinear(ir1, it1, ic1, rholog, tlog, yet, dbeadttab);
        nse_state.e_nu = trilinear(ir1, it1, ic1, rholog, tlog, yet, enutab);

        // massfractab is 2-d, so we wrap the access in a lambda already
        // indexing the component

        if (! skip_X_fill) {
            for (int n = 1; n <= NumSpec; n++) {
                Real _X = trilinear(ir1, it1, ic1, rholog, tlog, yet,
                                    [=] (const int i) {return massfractab(n, i);});
                nse_state.X[n-1] = std::clamp(_X, 0.0_rt, 1.0_rt);
            }
        }

    } else {

        // for a cubic interpolant, we need 4 points that span the data value
        // for temperature, these will be it0, it0+1, it0+2, it0+3
        // with the idea that the temperature we want is between it0+1 and it0+2
        // so we offset one to the left and also ensure that we don't go off the table

        int ir0 = nse_get_logrho_index(rholog) - 1;
        ir0 = std::clamp(ir0, 1, nse_table_size::nden-3);

        int it0 = nse_get_logT_index(tlog) - 1;
        it0 = std::clamp(it0, 1, nse_table_size::ntemp-3);

        int ic0 = nse_get_ye_index(yet) - 1;
        ic0 = std::clamp(ic0, 1, nse_table_size::nye-3);

        nse_state.abar = tricubic(ir0, it0, ic0, rholog, tlog, yet, abartab);
        nse_state.bea = tricubic(ir0, it0, ic0, rholog, tlog, yet, beatab);
        nse_state.dyedt = tricubic(ir0, it0, ic0, rholog, tlog, yet, dyedttab);
        nse_state.dabardt = tricubic(ir0, it0, ic0, rholog, tlog, yet, dabardttab);
        nse_state.dbeadt = tricubic(ir0, it0, ic0, rholog, tlog, yet, dbeadttab);
        nse_state.e_nu = tricubic(ir0, it0, ic0, rholog, tlog, yet, enutab);

        // massfractab is 2-d, so we wrap the access in a lambda already
        // indexing the component

        if (! skip_X_fill) {
            for (int n = 1; n <= NumSpec; n++) {
                Real _X = tricubic(ir0, it0, ic0, rholog, tlog, yet,
                                   [=] (const int i) {return massfractab(n, i);});
                nse_state.X[n-1] = std::clamp(_X, 0.0_rt, 1.0_rt);
            }
        }
    }

}


///
/// compute the temperature derivative of the table quantity data
/// at the point T, rho, ye by using cubic interpolation
///
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real
nse_interp_dT(const Real temp, const Real rho, const Real ye, const T& data) {

    Real rholog = std::log10(rho);
    {
        Real rmin = nse_table_size::logrho_min;
        Real rmax = nse_table_size::logrho_max;

        rholog = std::clamp(rholog, rmin, rmax);
    }

    Real tlog = std::log10(temp);
    {
        Real tmin = nse_table_size::logT_min;
        Real tmax = nse_table_size::logT_max;

        tlog = std::clamp(tlog, tmin, tmax);
    }

    Real yet = ye;
    {
        Real yemin = nse_table_size::ye_min;
        Real yemax = nse_table_size::ye_max;

        yet = std::clamp(yet, yemin, yemax);
    }

    int ir0 = nse_get_logrho_index(rholog) - 1;
    ir0 = std::clamp(ir0, 1, nse_table_size::nden-3);

    int it0 = nse_get_logT_index(tlog) - 1;
    it0 = std::clamp(it0, 1, nse_table_size::ntemp-3);

    int ic0 = nse_get_ye_index(yet) - 1;
    ic0 = std::clamp(ic0, 1, nse_table_size::nye-3);

    // note: this is returning the derivative wrt log10(T), so we need to
    // convert to d/dT

    Real ddatadT = tricubic_dT(ir0, it0, ic0, rholog, tlog, yet, data) / (std::log(10.0_rt) * temp);

    return ddatadT;

}

#endif
