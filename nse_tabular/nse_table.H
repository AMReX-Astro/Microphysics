#ifndef NSE_TABLE_H
#define NSE_TABLE_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <sstream>

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include <extern_parameters.H>
#include <nse_table_data.H>
#include <nse_table_size.H>
#include <nse_table_type.H>


using namespace amrex::literals;
using namespace network_rp;

///
/// given a rho, T, and Ye index, return the 1-d index into the NSE table arrays
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int nse_idx(const int ir, const int it, const int ic) {
    // this uses a 1-based indexing
    int idx = (ir-1) * nse_table_size::ntemp * nse_table_size::nye + (it-1) * nse_table_size::nye + ic;
    AMREX_ASSERT(idx >= 1 && idx <= nse_table::npts);
    return idx;
}

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

  // skip the header -- it is 4 lines
  std::string line;
  std::getline(nse_table_file, line);
  std::getline(nse_table_file, line);
  std::getline(nse_table_file, line);
  std::getline(nse_table_file, line);

  for (int irho = 1; irho <= nse_table_size::nden; irho++) {
      for (int it = 1; it <= nse_table_size::ntemp; it++) {
          for (int iye = 1; iye <= nse_table_size::nye; iye++) {
              const int j = nse_idx(irho, it, iye);

              const auto row_error = [=] (const std::string& reason) {
                  amrex::Error("NSE table " + nse_table_size::table_name +
                               ", line " + std::to_string(j + 4) +
                               ", index " + std::to_string(j) +
                               " (irho=" + std::to_string(irho) +
                               ", it=" + std::to_string(it) +
                               ", iye=" + std::to_string(iye) + "): " + reason);
              };
              if (!std::getline(nse_table_file, line) || line.empty()) {
                  row_error("missing or empty row");
              }
              std::istringstream data(line);
              amrex::Real ttemp, tdens, tye;
              data >> tdens >> ttemp >> tye;
              data >> nse_table::abartab(j)
                   >> nse_table::beatab(j)
                   >> nse_table::dyedttab(j)
                   >> nse_table::dabardttab(j)
                   >> nse_table::dbeadttab(j)
                   >> nse_table::enutab(j);
              for (int n = 1; n <= NumSpec; n++) {
                  data >> nse_table::massfractab(n, j);
              }
              if (data.fail()) {
                  row_error("incomplete or invalid numeric data");
              }
              data >> std::ws;
              if (!data.eof()) {
                  row_error("unexpected trailing data");
              }

              for (const auto value : {tdens, ttemp, tye,
                                       nse_table::abartab(j), nse_table::beatab(j),
                                       nse_table::dyedttab(j), nse_table::dabardttab(j),
                                       nse_table::dbeadttab(j), nse_table::enutab(j)}) {
                  if (!std::isfinite(value)) {
                      row_error("non-finite coordinate or table value");
                  }
              }
              for (int n = 1; n <= NumSpec; n++) {
                  if (!std::isfinite(nse_table::massfractab(n, j))) {
                      row_error("non-finite mass fraction");
                  }
              }

              // Coordinates are printed to ten decimal places. Also allow
              // roundoff when amrex::Real uses single precision.
              const auto matches_grid = [] (amrex::Real value, amrex::Real expected) {
                  const amrex::Real tol = 1.e-10_rt +
                      4 * std::numeric_limits<amrex::Real>::epsilon() * std::abs(expected);
                  return std::abs(value - expected) <= tol;
              };
              const amrex::Real logrho = nse_table_size::logrho_min +
                  (irho - 1) * nse_table_size::dlogrho;
              const amrex::Real logT = nse_table_size::logT_min +
                  (it - 1) * nse_table_size::dlogT;
              const amrex::Real ye = nse_table_size::ye_max -
                  (iye - 1) * nse_table_size::dye;
              if (!matches_grid(tdens, logrho) || !matches_grid(ttemp, logT) ||
                  !matches_grid(tye, ye)) {
                  row_error("coordinates do not match the expected (logrho, logT, Ye) grid: (" +
                            std::to_string(logrho) + ", " + std::to_string(logT) +
                            ", " + std::to_string(ye) + ")");
              }
          }
      }
  }

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real nse_table_logT(const int it) {
    return nse_table_size::logT_min + static_cast<amrex::Real>(it-1) * nse_table_size::dlogT;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real nse_table_logrho(const int ir) {
    return nse_table_size::logrho_min + static_cast<amrex::Real>(ir-1) * nse_table_size::dlogrho;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real nse_table_ye(const int ic) {
    return nse_table_size::ye_max - static_cast<amrex::Real>(ic-1) * nse_table_size::dye;
}

// return the index in the table such that logrho[irho] < input density
// note: this is a 1-based index
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int nse_get_logrho_index(const amrex::Real logrho) {

    int ir0 = static_cast<int>((logrho - nse_table_size::logrho_min) *
                               nse_table_size::inv_dlogrho - 1.e-6_rt);
    return ir0 + 1;
}

// return the index in the table such that logT[it] < input temperature
// note: this is a 1-based index
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int nse_get_logT_index(const amrex::Real logT) {

    int it0 = static_cast<int>((logT - nse_table_size::logT_min) *
                               nse_table_size::inv_dlogT - 1.e-6_rt);
    return it0 + 1;
}

// return the index in the table such that ye[ic] < input Ye
// note: this is a 1-based index
AMREX_GPU_HOST_DEVICE AMREX_INLINE
int nse_get_ye_index(const amrex::Real ye) {

    int ic0 = static_cast<int>((nse_table_size::ye_max - ye) *
                               nse_table_size::inv_dye - 1.0e-6_rt);
    return ic0 + 1;
}


// Slope Limiter to preserve monotonicity
// See C. Moler, Numerical Computing with Matlab, 2004.
// :doi:`10.1137/1.9780898717952`, Chapter 3.4
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real
limit_slope(const amrex::Real dm, const amrex::Real dp) {

    // note: we assume that dx is constant in our interval

    if (dm * dp <= 0.0_rt) {
        return 0.0_rt;
    }

    return 2.0_rt / (1.0_rt / dm + 1.0_rt / dp);
}


///
/// given 4 points (xs, fs), with spacing dx, return the interplated
/// value of f at point x by fitting a monotone cubic Hermite
/// polynomial to the points:
///
/// H(x) = f0 H_{1,0}(x) + f1 H_{1,1}(x) + f'0 Ĥ_{1,0}(x) + f'1 Ĥ_{1,0}(x)
///
/// The slope f' is limited to preserve monotonicity.
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real cubic(const amrex::Real* xs, const amrex::Real* fs,
                  const amrex::Real dx, const amrex::Real x,
                  const bool lo_bnd, const bool hi_bnd) {

    // given the 4 points xs, the point x is in [xs[1], xs[2]]
    // we'll use i0 = 1 as the reference for the interpolation

    // we have already guaranteed that all of the xs are in the
    // domain, so we don't need to worry about boundaries.

    // Compute the slope at i0 and i0 + 1
    const auto i0 = 1;

    const auto dx_inv = 1.0_rt / dx;

    amrex::Real m0, m1;
    m1 = m0 = (fs[i0+1] - fs[i0]) * dx_inv;

    // Compute the slope from the other side and limit the slope to
    // preserve monotonicity
    if (! lo_bnd) {
        amrex::Real d_im1 = (fs[i0] - fs[i0-1]) * dx_inv;
        m0 = limit_slope(d_im1, m0);
    }

    if (! hi_bnd) {
        amrex::Real d_ip1 = (fs[i0+2] - fs[i0+1]) * dx_inv;
        m1 = limit_slope(m1, d_ip1);
    }

    // Compute cubic Hermite polynomial basis
    const amrex::Real t = (x - xs[i0]) * dx_inv;
    const amrex::Real t2 = t * t;
    const amrex::Real t3 = t2 * t;

    const amrex::Real H1 = -2.0_rt * t3 + 3.0_rt * t2;
    const amrex::Real H0 = -H1 + 1.0_rt;
    const amrex::Real Hhat1 = dx * (t3 - t2);
    const amrex::Real Hhat0 = dx * (t3 - 2.0_rt * t2 + t);

    // Construct the interpolant
    // H_3(x) = f0 H_{1,0}(x) + f1 H_{1,1}(x) + f'0 Ĥ_{1,0}(x) + f'1 Ĥ_{1,0}(x)

    const amrex::Real fv = fs[i0] * H0 + fs[i0+1] * H1 + m0 * Hhat0 + m1 * Hhat1;
    return fv;
}

///
/// given 4 points (xs, fs), with spacing dx between the xs, return
/// the derivative of f at point x by fitting a cubic to the
/// points and differentiating the interpolant
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real cubic_deriv(const amrex::Real* xs, const amrex::Real* fs,
                        const amrex::Real dx, const amrex::Real x,
                        const bool lo_bnd, const bool hi_bnd) {

    // given the 4 points xs, the point x is in [xs[1], xs[2]]
    // we'll use i0 = 1 as the reference for the interpolation

    // we have already guaranteed that all of the xs are in the
    // domain, so we don't need to worry about boundaries.

    // Compute the slope at i0 and i0 + 1
    const auto i0 = 1;

    const auto dx_inv = 1.0_rt / dx;

    amrex::Real m0, m1;
    m1 = m0 = (fs[i0+1] - fs[i0]) * dx_inv;

    // Compute the slope from the other side and limit the slope to
    // preserve monotonicity
    if (! lo_bnd) {
        amrex::Real d_im1 = (fs[i0] - fs[i0-1]) * dx_inv;
        m0 = limit_slope(d_im1, m0);
    }

    if (! hi_bnd) {
        amrex::Real d_ip1 = (fs[i0+2] - fs[i0+1]) * dx_inv;
        m1 = limit_slope(m1, d_ip1);
    }

    // Compute cubic Hermite polynomial basis
    const amrex::Real t = (x - xs[i0]) * dx_inv;
    const amrex::Real t2 = t * t;

    const amrex::Real dH1dx = -6.0_rt * (t2 - t) * dx_inv;
    const amrex::Real dH0dx = -dH1dx;
    const amrex::Real dHhat1dx = 3.0_rt * t2 - 2.0_rt * t;
    const amrex::Real dHhat0dx = 3.0_rt * t2 - 4.0_rt * t + 1.0_rt;

    const amrex::Real dfdx = fs[i0] * dH0dx + fs[i0+1] * dH1dx +
                             m0 * dHhat0dx + m1 * dHhat1dx;

    return dfdx;

}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real trilinear(const int ir1, const int it1, const int ic1,
               const amrex::Real rho, const amrex::Real temp, const amrex::Real ye, const T& data) {

    // find the eight interpolation points in the 1D arrays

    int it1r1c1 = nse_idx(ir1,   it1,   ic1);
    int it1r1c2 = nse_idx(ir1,   it1,   ic1+1);
    int it1r2c1 = nse_idx(ir1+1, it1,   ic1);
    int it1r2c2 = nse_idx(ir1+1, it1,   ic1+1);
    int it2r1c1 = nse_idx(ir1,   it1+1, ic1);
    int it2r1c2 = nse_idx(ir1,   it1+1, ic1+1);
    int it2r2c1 = nse_idx(ir1+1, it1+1, ic1);
    int it2r2c2 = nse_idx(ir1+1, it1+1, ic1+1);

    amrex::Real t0 = nse_table_logT(it1);
    amrex::Real r0 = nse_table_logrho(ir1);
    amrex::Real x0 = nse_table_ye(ic1);

    amrex::Real td = (temp - t0) * nse_table_size::inv_dlogT;
    amrex::Real rd = (rho - r0) * nse_table_size::inv_dlogrho;
    amrex::Real xd = (x0 - ye) * nse_table_size::inv_dye;
    xd = amrex::max(0.0_rt, xd);

    amrex::Real omtd = 1.0_rt - td;
    amrex::Real omrd = 1.0_rt - rd;
    amrex::Real omxd = 1.0_rt - xd;

    amrex::Real val = data(it1r1c1) * omtd * omrd * omxd +
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
amrex::Real tricubic(const int ir0, const int it0, const int ic0,
              const amrex::Real rho, const amrex::Real temp, const amrex::Real ye, const T& data) {

    // as we come in, we expect temp to be in it0:it0+1, and we also
    // are assured that it0 and it0+1 don't go out of bounds.  Likewise
    // for rho and ye.  Our 1D cubic interpolation routine will handle
    // the boundary conditions for us.

    const amrex::Real yes[] = {nse_table_ye(std::max(1, ic0-1)),
                               nse_table_ye(ic0),
                               nse_table_ye(ic0+1),
                               nse_table_ye(std::min(nse_table_size::nye, ic0+2))};

    const amrex::Real Ts[] = {nse_table_logT(std::max(1, it0-1)),
                              nse_table_logT(it0),
                              nse_table_logT(it0+1),
                              nse_table_logT(std::min(nse_table_size::ntemp, it0+2))};

    const amrex::Real rhos[] = {nse_table_logrho(std::max(1, ir0-1)),
                                nse_table_logrho(ir0),
                                nse_table_logrho(ir0+1),
                                nse_table_logrho(std::min(nse_table_size::nden, ir0+2))};

    // first do the 16 ye interpolations

    // the first index will be rho and the second will be T
    amrex::Real d1[4][4];

    bool lo_bnd = (ic0 == 1);
    bool hi_bnd = (ic0 + 1 == nse_table_size::nye);

    for (int ii = 0; ii < 4; ++ii) {
        int idx_r = std::max(1, std::min(ir0 - 1 + ii, nse_table_size::nden));

        for (int jj = 0; jj < 4; ++jj) {
            int idx_t = std::max(1, std::min(it0 - 1 + jj, nse_table_size::ntemp));

            const amrex::Real _d[] = {data(nse_idx(idx_r, idx_t, std::max(1, ic0-1))),
                                      data(nse_idx(idx_r, idx_t, ic0)),
                                      data(nse_idx(idx_r, idx_t, ic0+1)),
                                      data(nse_idx(idx_r, idx_t, std::min(ic0+2, nse_table_size::nye)))};

            // note that the ye values are monotonically decreasing,
            // so the "dx" needs to be negative
            d1[ii][jj] = cubic(yes, _d, -nse_table_size::dye, ye, lo_bnd, hi_bnd);
        }
    }

    // now do the 4 T interpolations (one in each rho plane)

    amrex::Real d2[4];

    lo_bnd = (it0 == 1);
    hi_bnd = (it0 + 1 == nse_table_size::ntemp);

    for (int ii = 0; ii < 4; ++ii) {
        const amrex::Real _d[] = {d1[ii][0], d1[ii][1], d1[ii][2], d1[ii][3]};
        d2[ii] = cubic(Ts, _d, nse_table_size::dlogT, temp, lo_bnd, hi_bnd);
    }

    // finally do the remaining interpolation over rho

    lo_bnd = (ir0 == 1);
    hi_bnd = (ir0 + 1 == nse_table_size::nden);

    amrex::Real val = cubic(rhos, d2, nse_table_size::dlogrho, rho, lo_bnd, hi_bnd);

    return val;

}

///
/// take the temperature derivative of a table quantity by differentiating
/// the cubic interpolant
///
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real tricubic_dT(const int ir0, const int it0, const int ic0,
                 const amrex::Real rho, const amrex::Real temp, const amrex::Real ye, const T& data) {

    // as we come in, we expect temp to be in it0:it0+1, and we also
    // are assured that it0 and it0+1 don't go out of bounds.  Likewise
    // for rho and ye.  Our 1D cubic interpolation routine will handle
    // the boundary conditions for us.

    const amrex::Real yes[] = {nse_table_ye(std::max(1, ic0-1)),
                               nse_table_ye(ic0),
                               nse_table_ye(ic0+1),
                               nse_table_ye(std::min(nse_table_size::nye, ic0+2))};

    const amrex::Real Ts[] = {nse_table_logT(std::max(1, it0-1)),
                              nse_table_logT(it0),
                              nse_table_logT(it0+1),
                              nse_table_logT(std::min(nse_table_size::ntemp, it0+2))};

    const amrex::Real rhos[] = {nse_table_logrho(std::max(1, ir0-1)),
                                nse_table_logrho(ir0),
                                nse_table_logrho(ir0+1),
                                nse_table_logrho(std::min(nse_table_size::nden, ir0+2))};

    // first do the 16 ye interpolations

    // the first index will be rho and the second will be T
    amrex::Real d1[4][4];

    bool lo_bnd = (ic0 == 1);
    bool hi_bnd = (ic0 + 1 == nse_table_size::nye);

    for (int ii = 0; ii < 4; ++ii) {
        int idx_r = std::max(1, std::min(ir0 - 1 + ii, nse_table_size::nden));

        for (int jj = 0; jj < 4; ++jj) {
            int idx_t = std::max(1, std::min(it0 - 1 + jj, nse_table_size::ntemp));

            const amrex::Real _d[] = {data(nse_idx(idx_r, idx_t, std::max(1, ic0-1))),
                                      data(nse_idx(idx_r, idx_t, ic0)),
                                      data(nse_idx(idx_r, idx_t, ic0+1)),
                                      data(nse_idx(idx_r, idx_t, std::min(ic0+2, nse_table_size::nye)))};

            // note that the ye values are monotonically decreasing,
            // so the "dx" needs to be negative
            d1[ii][jj] = cubic(yes, _d, -nse_table_size::dye, ye, lo_bnd, hi_bnd);
        }
    }

    // now do the 4 rho interpolations (one in each T plane)

    amrex::Real d2[4];

    lo_bnd = (ir0 == 1);
    hi_bnd = (ir0 + 1 == nse_table_size::nden);

    for (int jj = 0; jj < 4; ++jj) {
        const amrex::Real _d[] = {d1[0][jj], d1[1][jj], d1[2][jj], d1[3][jj]};
        d2[jj] = cubic(rhos, _d, nse_table_size::dlogrho, rho, lo_bnd, hi_bnd);
    }

    // finally do the remaining interpolation over T, but return
    // the derivative of the interpolant

    lo_bnd = (it0 == 1);
    hi_bnd = (it0 + 1 == nse_table_size::ntemp);

    amrex::Real val = cubic_deriv(Ts, d2, nse_table_size::dlogT, temp, lo_bnd, hi_bnd);

    return val;

}


///
/// take the density derivative of a table quantity by differentiating
/// the cubic interpolant
///
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real tricubic_drho(const int ir0, const int it0, const int ic0,
                   const amrex::Real rho, const amrex::Real temp, const amrex::Real ye, const T& data) {

    // as we come in, we expect temp to be in it0:it0+1, and we also
    // are assured that it0 and it0+1 don't go out of bounds.  Likewise
    // for rho and ye.  Our 1D cubic interpolation routine will handle
    // the boundary conditions for us.

    const amrex::Real yes[] = {nse_table_ye(std::max(1, ic0-1)),
                               nse_table_ye(ic0),
                               nse_table_ye(ic0+1),
                               nse_table_ye(std::min(nse_table_size::nye, ic0+2))};

    const amrex::Real Ts[] = {nse_table_logT(std::max(1, it0-1)),
                              nse_table_logT(it0),
                              nse_table_logT(it0+1),
                              nse_table_logT(std::min(nse_table_size::ntemp, it0+2))};

    const amrex::Real rhos[] = {nse_table_logrho(std::max(1, ir0-1)),
                                nse_table_logrho(ir0),
                                nse_table_logrho(ir0+1),
                                nse_table_logrho(std::min(nse_table_size::nden, ir0+2))};

    // first do the 16 ye interpolations

    // the first index will be rho and the second will be T
    amrex::Real d1[4][4];

    bool lo_bnd = (ic0 == 1);
    bool hi_bnd = (ic0 + 1 == nse_table_size::nye);

    for (int ii = 0; ii < 4; ++ii) {
        int idx_r = std::max(1, std::min(ir0 - 1 + ii, nse_table_size::nden));

        for (int jj = 0; jj < 4; ++jj) {
            int idx_t = std::max(1, std::min(it0 - 1 + jj, nse_table_size::ntemp));

            const amrex::Real _d[] = {data(nse_idx(idx_r, idx_t, std::max(1, ic0-1))),
                                      data(nse_idx(idx_r, idx_t, ic0)),
                                      data(nse_idx(idx_r, idx_t, ic0+1)),
                                      data(nse_idx(idx_r, idx_t, std::min(ic0+2, nse_table_size::nye)))};

            // note that the ye values are monotonically decreasing,
            // so the "dx" needs to be negative
            d1[ii][jj] = cubic(yes, _d, -nse_table_size::dye, ye, lo_bnd, hi_bnd);
        }
    }

    // now do the 4 T interpolations (one in each rho plane)

    amrex::Real d2[4];

    lo_bnd = (it0 == 1);
    hi_bnd = (it0 + 1 == nse_table_size::ntemp);

    for (int ii = 0; ii < 4; ++ii) {

        const amrex::Real _d[] = {d1[ii][0], d1[ii][1], d1[ii][2], d1[ii][3]};
        d2[ii] = cubic(Ts, _d, nse_table_size::dlogT, temp, lo_bnd, hi_bnd);
    }

    // finally do the remaining interpolation over rho, but return
    // the derivative of the interpolant


    lo_bnd = (ir0 == 1);
    hi_bnd = (ir0 + 1 == nse_table_size::nden);

    amrex::Real val = cubic_deriv(rhos, d2, nse_table_size::dlogrho, rho, lo_bnd, hi_bnd);

    return val;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_interp(nse_table_t& nse_state, bool skip_X_fill=false) {

    // if skip_X_fill = true then we don't fill X[] with the mass fractions.

    using namespace nse_table;
    using namespace AuxZero;

    amrex::Real rholog = std::log10(nse_state.rho);
    {
        amrex::Real rmin = nse_table_size::logrho_min;
        amrex::Real rmax = nse_table_size::logrho_max;

        rholog = amrex::Clamp(rholog, rmin, rmax);
    }

    amrex::Real tlog = std::log10(nse_state.T);
    {
        amrex::Real tmin = nse_table_size::logT_min;
        amrex::Real tmax = nse_table_size::logT_max;

        tlog = amrex::Clamp(tlog, tmin, tmax);
    }

    amrex::Real yet = nse_state.Ye;
    {
        amrex::Real yemin = nse_table_size::ye_min;
        amrex::Real yemax = nse_table_size::ye_max;

        yet = amrex::Clamp(yet, yemin, yemax);
    }

    if (nse_table_interp_linear) {

        int ir1 = nse_get_logrho_index(rholog);
        int it1 = nse_get_logT_index(tlog);
        int ic1 = nse_get_ye_index(yet);

        nse_state.abar = trilinear(ir1, it1, ic1, rholog, tlog, yet, abartab);
        nse_state.bea = trilinear(ir1, it1, ic1, rholog, tlog, yet, beatab);
        nse_state.dyedt = trilinear(ir1, it1, ic1, rholog, tlog, yet, dyedttab);
        nse_state.dbeadt = trilinear(ir1, it1, ic1, rholog, tlog, yet, dbeadttab);
        nse_state.e_nu = trilinear(ir1, it1, ic1, rholog, tlog, yet, enutab);

        // massfractab is 2-d, so we wrap the access in a lambda already
        // indexing the component

        if (! skip_X_fill) {
            for (int n = 1; n <= NumSpec; n++) {
                amrex::Real _X = trilinear(ir1, it1, ic1, rholog, tlog, yet,
                                    [=] (const int i) {return massfractab(n, i);});
                nse_state.X[n-1] = amrex::Clamp(_X, 0.0_rt, 1.0_rt);
            }
        }

    } else {

        // for a cubic interpolant, we need 4 points that span the data value
        // for temperature, these will be it0-1, it0, it0+1, it0+2
        // with the idea that the temperature we want is between it0 and it0+1

        // we will handle the case where it0 or it0+1 is at the bounds of
        // the table in the caller.

        int ir0 = nse_get_logrho_index(rholog);
        int it0 = nse_get_logT_index(tlog);
        int ic0 = nse_get_ye_index(yet);

        nse_state.abar = tricubic(ir0, it0, ic0, rholog, tlog, yet, abartab);
        nse_state.bea = tricubic(ir0, it0, ic0, rholog, tlog, yet, beatab);
        nse_state.dyedt = tricubic(ir0, it0, ic0, rholog, tlog, yet, dyedttab);
        nse_state.dbeadt = tricubic(ir0, it0, ic0, rholog, tlog, yet, dbeadttab);
        nse_state.e_nu = tricubic(ir0, it0, ic0, rholog, tlog, yet, enutab);

        // massfractab is 2-d, so we wrap the access in a lambda already
        // indexing the component

        if (! skip_X_fill) {
            for (int n = 1; n <= NumSpec; n++) {
                amrex::Real _X = tricubic(ir0, it0, ic0, rholog, tlog, yet,
                                   [=] (const int i) {return massfractab(n, i);});
                nse_state.X[n-1] = amrex::Clamp(_X, 0.0_rt, 1.0_rt);
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
amrex::Real
nse_interp_dT(const amrex::Real temp, const amrex::Real rho, const amrex::Real ye, const T& data) {

    amrex::Real rholog = std::log10(rho);
    {
        amrex::Real rmin = nse_table_size::logrho_min;
        amrex::Real rmax = nse_table_size::logrho_max;

        rholog = amrex::Clamp(rholog, rmin, rmax);
    }

    amrex::Real tlog = std::log10(temp);
    {
        amrex::Real tmin = nse_table_size::logT_min;
        amrex::Real tmax = nse_table_size::logT_max;

        tlog = amrex::Clamp(tlog, tmin, tmax);
    }

    amrex::Real yet = ye;
    {
        amrex::Real yemin = nse_table_size::ye_min;
        amrex::Real yemax = nse_table_size::ye_max;

        yet = amrex::Clamp(yet, yemin, yemax);
    }

    int ir0 = nse_get_logrho_index(rholog);
    int it0 = nse_get_logT_index(tlog);
    int ic0 = nse_get_ye_index(yet);

    // note: this is returning the derivative wrt log10(T), so we need to
    // convert to d/dT

    amrex::Real ddatadT = tricubic_dT(ir0, it0, ic0, rholog, tlog, yet, data) / (M_LN10 * temp);

    return ddatadT;


}


///
/// compute the density derivative of the table quantity data
/// at the point T, rho, ye by using cubic interpolation
///
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real
nse_interp_drho(const amrex::Real temp, const amrex::Real rho, const amrex::Real ye, const T& data) {

    amrex::Real rholog = std::log10(rho);
    {
        amrex::Real rmin = nse_table_size::logrho_min;
        amrex::Real rmax = nse_table_size::logrho_max;

        rholog = amrex::Clamp(rholog, rmin, rmax);
    }

    amrex::Real tlog = std::log10(temp);
    {
        amrex::Real tmin = nse_table_size::logT_min;
        amrex::Real tmax = nse_table_size::logT_max;

        tlog = amrex::Clamp(tlog, tmin, tmax);
    }

    amrex::Real yet = ye;
    {
        amrex::Real yemin = nse_table_size::ye_min;
        amrex::Real yemax = nse_table_size::ye_max;

        yet = amrex::Clamp(yet, yemin, yemax);
    }

    int ir0 = nse_get_logrho_index(rholog);
    int it0 = nse_get_logT_index(tlog);
    int ic0 = nse_get_ye_index(yet);

    // note: this is returning the derivative wrt log10(rho), so we need to
    // convert to d/drho

    amrex::Real ddatadrho = tricubic_drho(ir0, it0, ic0, rholog, tlog, yet, data) / (M_LN10 * rho);

    return ddatadrho;

}

#endif
