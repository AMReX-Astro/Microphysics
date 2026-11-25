#ifndef TABLE_RATES_H
#define TABLE_RATES_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <AMReX_Array.H>
#include <AMReX_Math.H>
#include <interp_tools.H>

using namespace amrex::literals;


// Table is expected to be in terms of dens*ye and temp (logarithmic, cgs units)
// Table energy units are expected in terms of ergs

// all tables are expected to have columns:
// Log(rhoY)     Log(T)   mu    dQ    Vs    Log(e-cap-rate)   Log(nu-energy-loss)  Log(gamma-energy)
// Log(g/cm^3)   Log(K)   erg   erg   erg   Log(1/s)          Log(erg/s)           Log(erg/s)
//

constexpr int num_tables = 1;

enum TableVars : std::uint8_t
{
    jtab_rate    = 1,
    jtab_nuloss  = 2,
    jtab_gamma   = 3,
    num_vars = jtab_gamma
};


struct table_t
{
    int ntemp;
    int nrhoy;
    int nvars;
    int nheader;
};


struct weak_rate_t
{
    amrex::Real rate;
    amrex::Real drate_dT;
    amrex::Real drate_drhoye;
    amrex::Real enu;
    amrex::Real denu_dT;
    amrex::Real denu_drhoye;
    amrex::Real gamma;
};


namespace rate_tables
{
    // n --> p
    inline AMREX_GPU_DEVICE table_t j_n_p_meta{.ntemp=13, .nrhoy=11, .nvars=3, .nheader=5};
    inline AMREX_GPU_DEVICE amrex::Array1D<amrex::Real, 1, 11> j_n_p_rhoy{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0};
    inline AMREX_GPU_DEVICE amrex::Array1D<amrex::Real, 1, 13> j_n_p_temp{7.0, 8.0, 8.30103, 8.60206, 8.845098, 9.0, 9.176091, 9.30103, 9.477121, 9.69897, 10.0, 10.477121, 11.0};
    // Array3D is column-major (Fortran-ordering).  T varies fastest, then rho Ye, then the component
    inline AMREX_GPU_DEVICE amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, num_vars>
         j_n_p_data{
-2.959,          -2.959,          -2.959,          -2.959,       -2.958959,       -2.957841,       -2.941107,       -2.874417,       -2.540055,       -1.720786,        -0.41411,        1.802004,           4.347,  // RATE
-2.959,          -2.959,          -2.959,          -2.959,       -2.958984,       -2.957898,       -2.941192,       -2.874417,       -2.540055,       -1.720786,        -0.41411,        1.802004,           4.347,
-2.959,          -2.959,          -2.959,          -2.959,       -2.958998,       -2.958314,       -2.941776,       -2.874975,       -2.540055,       -1.721736,        -0.41411,        1.802004,           4.347,
-2.961,          -2.961,          -2.961,          -2.961,           -2.96,       -2.959896,       -2.947722,       -2.881381,       -2.544233,       -1.721736,        -0.41411,        1.802004,           4.347,
-2.98,          -2.979,          -2.979,          -2.978,          -2.975,        -2.97299,       -2.967462,       -2.925218,       -2.577269,       -1.729383,       -0.415108,        1.802004,           4.347,
-3.14,           -3.14,          -3.138,          -3.131,          -3.116,       -3.097999,       -3.069766,       -3.041782,       -2.813017,       -1.806664,       -0.423096,        1.801004,           4.347,
-5.499,          -5.442,          -5.306,          -4.989,          -4.592,          -4.297,       -3.954978,       -3.726186,       -3.408689,       -2.366458,       -0.508967,        1.798004,           4.347,
-99.69797,          -65.57,         -35.616,         -20.209,         -13.299,         -10.394,       -8.019984,       -6.772438,       -5.436532,       -3.791028,       -1.150873,        1.768005,           4.346,
-99.69797,       -99.69797,       -99.69797,         -54.663,           -33.0,         -24.198,      -17.245984,       -13.71544,      -10.109646,       -6.677964,       -2.765936,        1.478006,           4.338,
-99.69797,       -99.69797,       -99.69797,       -99.69797,         -75.759,         -54.136,      -37.214984,       -28.70344,      -20.122703,      -12.726964,       -5.885901,         0.27701,           4.255,
-99.69797,       -99.69797,       -99.69797,       -99.69797,       -99.69797,       -99.69797,      -80.281984,       -61.00844,      -41.669646,      -25.673964,      -12.403901,        -2.04999,           3.625,
-9.07229,        -9.07229,        -9.07229,        -9.07229,        -9.07229,        -9.06729,        -8.99729,        -8.76529,        -8.08629,        -7.02729,        -5.49029,        -2.85229,       0.1987104,  // NU
-9.07229,        -9.07229,        -9.07229,        -9.07229,        -9.07229,        -9.06729,        -8.99729,        -8.76529,        -8.08629,        -7.02729,        -5.49029,        -2.85229,       0.1987104,
-9.07229,        -9.07229,        -9.07229,        -9.07229,        -9.07229,        -9.06929,        -9.00029,        -8.76729,        -8.08729,        -7.02729,        -5.49029,        -2.85229,       0.1987104,
-9.07629,        -9.07529,        -9.07529,        -9.07529,        -9.07429,        -9.07329,        -9.02029,        -8.78329,        -8.09129,        -7.02829,        -5.49129,        -2.85229,       0.1987104,
-9.10529,        -9.10429,        -9.10329,        -9.10029,        -9.09629,        -9.09229,        -9.07629,        -8.90629,        -8.13929,        -7.03629,        -5.49129,        -2.85229,       0.1987104,
-9.34529,        -9.34329,        -9.34029,        -9.32529,        -9.29529,        -9.26229,        -9.21529,        -9.15329,        -8.51229,        -7.11629,        -5.50029,        -2.85229,       0.1987104,
-12.57129,       -12.46129,       -12.21829,       -11.70729,       -11.12229,       -10.71129,       -10.25629,        -9.96229,        -9.42029,        -7.71429,        -5.58729,        -2.85529,       0.1987104,
-105.7943,       -72.95929,       -42.71029,       -27.01829,       -19.88929,       -16.85929,       -14.36829,       -13.05429,       -11.53729,        -9.20529,        -6.23429,        -2.88629,       0.1977104,
-105.7943,       -105.7943,       -105.7943,       -61.47129,       -39.59029,       -30.66429,       -23.59429,       -19.99829,       -16.21129,       -12.09729,        -7.85229,        -3.17929,       0.1897104,
-105.7943,       -105.7943,       -105.7943,       -105.7943,       -82.34929,       -60.60229,       -43.56329,       -34.98629,       -26.22429,       -18.14629,       -10.97129,        -4.38329,       0.1057104,
-105.7943,       -105.7943,       -105.7943,       -105.7943,       -105.7943,       -105.7943,       -86.63029,       -67.29129,       -47.77129,       -31.09329,       -17.48929,        -6.71029,      -0.5272896,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,  // GAMMA
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
         };

}


AMREX_INLINE AMREX_GPU_DEVICE
amrex::Real
evaluate_linear_1d(const amrex::Real fhi, const amrex::Real flo, const amrex::Real xhi, const amrex::Real xlo, const amrex::Real x)
{
    // This function is a 1-D linear interpolator, that keeps x constant to xlo or xhi, based
    // on the side, if x is outside [xlo, xhi] respectively.

    amrex::Real xx = amrex::Clamp(x, xlo, xhi);
    amrex::Real f = flo + (fhi - flo) * (xx - xlo) / (xhi - xlo);

    return f;
}

AMREX_INLINE AMREX_GPU_DEVICE
amrex::Real
evaluate_linear_2d(const amrex::Real fip1jp1, const amrex::Real fip1j, const amrex::Real fijp1, const amrex::Real fij,
          const amrex::Real xhi, const amrex::Real xlo, const amrex::Real yhi, const amrex::Real ylo,
          const amrex::Real x, const amrex::Real y)
{
    // This is the 2-D linear interpolator, as an extension of evaluate_linear_1d.

    const amrex::Real u = (amrex::Clamp(x, xlo, xhi) - xlo) / (xhi - xlo);
    const amrex::Real v = (amrex::Clamp(y, ylo, yhi) - ylo) / (yhi - ylo);

    const amrex::Real E = fij;
    const amrex::Real C = fijp1 - fij;
    const amrex::Real B = fip1j - fij;
    const amrex::Real A = fip1jp1 - fip1j - fijp1 + fij;

    return A * u * v + B * u + C * v + E;
}


template<typename R, typename T, typename D>
AMREX_INLINE AMREX_GPU_DEVICE
amrex::Real
evaluate_vars(const int irhoy_lo, const int jtemp_lo,
              const R& log_rhoy_table, const T& log_temp_table, const D& data,
              const amrex::Real log_rhoy, const amrex::Real log_temp, const int component)
{
    // This function evaluates the 2-D interpolator, for several pairs of rho_ye and temperature.

    int irhoy_hi = irhoy_lo + 1;
    int jtemp_hi = jtemp_lo + 1;

    amrex::Real rhoy_lo = log_rhoy_table(irhoy_lo);
    amrex::Real rhoy_hi = log_rhoy_table(irhoy_hi);

    amrex::Real t_lo = log_temp_table(jtemp_lo);
    amrex::Real t_hi = log_temp_table(jtemp_hi);

    amrex::Real fij     = data(jtemp_lo, irhoy_lo, component);
    amrex::Real fip1j   = data(jtemp_lo, irhoy_hi, component);
    amrex::Real fijp1   = data(jtemp_hi, irhoy_lo, component);
    amrex::Real fip1jp1 = data(jtemp_hi, irhoy_hi, component);

    amrex::Real r = evaluate_linear_2d(fip1jp1, fip1j, fijp1, fij,
                                       rhoy_hi, rhoy_lo, t_hi, t_lo, log_rhoy, log_temp);

    return r;
}


template<typename R, typename T, typename D>
AMREX_INLINE AMREX_GPU_DEVICE
amrex::Real
evaluate_dr_dtemp(const int irhoy_lo, const int jtemp_lo,
                  const table_t& table_meta, const R& log_rhoy_table, const T& log_temp_table, const D& data,
                  const amrex::Real log_rhoy, const amrex::Real log_temp, const int component)
{

    // Differentiate the same bilinear log-rate interpolant as evaluate_vars.
    // At temperature knots, use the cell selected by find_index_extrap:
    // the cell to the right, except at the upper endpoint (the last cell).

    const int irhoy_hi = irhoy_lo + 1;
    const int jtemp_hi = jtemp_lo + 1;

    const amrex::Real log_temp_lo = log_temp_table(jtemp_lo);
    const amrex::Real log_temp_hi = log_temp_table(jtemp_hi);

    // The interpolated rate is constant outside the temperature table.
    if (log_temp < log_temp_lo || log_temp > log_temp_hi) {
        return 0.0_rt;
    }

    const amrex::Real log_rhoy_lo = log_rhoy_table(irhoy_lo);
    const amrex::Real log_rhoy_hi = log_rhoy_table(irhoy_hi);

    const amrex::Real fij     = data(jtemp_lo, irhoy_lo, component);
    const amrex::Real fip1j   = data(jtemp_lo, irhoy_hi, component);
    const amrex::Real fijp1   = data(jtemp_hi, irhoy_lo, component);
    const amrex::Real fip1jp1 = data(jtemp_hi, irhoy_hi, component);

    // compute the gradients in the T direction: dlog(r)/dlog(T) |_i and
    // dlog(r)/dlog(T) |_{i+1}

    const amrex::Real dlogT = log_temp_hi - log_temp_lo;
    const amrex::Real dlogr_dlogT_i = (fijp1 - fij) / dlogT;
    const amrex::Real dlogr_dlogT_ip1 = (fip1jp1 - fip1j) / dlogT;

    return evaluate_linear_1d(dlogr_dlogT_ip1, dlogr_dlogT_i,
                              log_rhoy_hi, log_rhoy_lo, log_rhoy);
}


template<typename R, typename T, typename D>
AMREX_INLINE AMREX_GPU_DEVICE
amrex::Real
evaluate_dr_drhoYe(const int irhoy_lo, const int jtemp_lo,
                   [[maybe_unused]] const table_t& table_meta,
                   const R& log_rhoy_table, const T& log_temp_table, const D& data,
                   const amrex::Real log_rhoy, const amrex::Real log_temp, const int component)
{

    // Differentiate the same bilinear log-rate interpolant as evaluate_vars.
    // At (rho Ye) knots, use the cell selected by find_index_extrap:
    // the cell to the right, except at the upper endpoint (the last cell).

    const int irhoy_hi = irhoy_lo + 1;
    const int jtemp_hi = jtemp_lo + 1;

    const amrex::Real log_rhoy_lo = log_rhoy_table(irhoy_lo);
    const amrex::Real log_rhoy_hi = log_rhoy_table(irhoy_hi);

    // The interpolated rate is constant outside the density table.
    if (log_rhoy < log_rhoy_lo || log_rhoy > log_rhoy_hi) {
        return 0.0_rt;
    }

    const amrex::Real log_temp_lo = log_temp_table(jtemp_lo);
    const amrex::Real log_temp_hi = log_temp_table(jtemp_hi);

    const amrex::Real fij     = data(jtemp_lo, irhoy_lo, component);
    const amrex::Real fip1j   = data(jtemp_lo, irhoy_hi, component);
    const amrex::Real fijp1   = data(jtemp_hi, irhoy_lo, component);
    const amrex::Real fip1jp1 = data(jtemp_hi, irhoy_hi, component);

    // compute the gradients in the rhoY direction: dlog(r)/dlog(ρY_e) |_j and
    // dlog(r)/dlog(ρY_e) |_{j+1}

    const amrex::Real dlogrhoy = log_rhoy_hi - log_rhoy_lo;
    const amrex::Real dlogr_dlogrhoy_j = (fip1j - fij) / dlogrhoy;
    const amrex::Real dlogr_dlogrhoy_jp1 = (fip1jp1 - fijp1) / dlogrhoy;

    return evaluate_linear_1d(dlogr_dlogrhoy_jp1, dlogr_dlogrhoy_j,
                              log_temp_hi, log_temp_lo, log_temp);
}


template <int do_T_derivatives, typename R, typename T, typename D>
AMREX_INLINE AMREX_GPU_DEVICE
void
get_entries(const int irhoy_lo, const int jtemp_lo,
            const table_t& table_meta, const R& log_rhoy_table, const T& log_temp_table, const D& data,
            const amrex::Real log_rhoy, const amrex::Real log_temp,
            weak_rate_t& entries)
{

    // only interpolate the quantities we currently use

    entries.rate = evaluate_vars(irhoy_lo, jtemp_lo,
                                 log_rhoy_table, log_temp_table, data,
                                 log_rhoy, log_temp, jtab_rate);

    entries.enu = evaluate_vars(irhoy_lo, jtemp_lo,
                                log_rhoy_table, log_temp_table, data,
                                log_rhoy, log_temp, jtab_nuloss);

    entries.gamma = evaluate_vars(irhoy_lo, jtemp_lo,
                                  log_rhoy_table, log_temp_table, data,
                                  log_rhoy, log_temp, jtab_gamma);

    if constexpr (do_T_derivatives) {
        entries.drate_dT =
            evaluate_dr_dtemp(irhoy_lo, jtemp_lo,
                              table_meta, log_rhoy_table, log_temp_table, data,
                              log_rhoy, log_temp, jtab_rate);

        entries.drate_drhoye =
            evaluate_dr_drhoYe(irhoy_lo, jtemp_lo,
                               table_meta, log_rhoy_table, log_temp_table, data,
                               log_rhoy, log_temp, jtab_rate);

        entries.denu_dT =
            evaluate_dr_dtemp(irhoy_lo, jtemp_lo,
                              table_meta, log_rhoy_table, log_temp_table, data,
                              log_rhoy, log_temp, jtab_nuloss);

        entries.denu_drhoye =
            evaluate_dr_drhoYe(irhoy_lo, jtemp_lo,
                               table_meta, log_rhoy_table, log_temp_table, data,
                               log_rhoy, log_temp, jtab_nuloss);

    }
}

template <int do_T_derivatives, typename R, typename T, typename D>
AMREX_INLINE AMREX_GPU_DEVICE
void
tabular_evaluate(const table_t& table_meta,
                 const R& log_rhoy_table, const T& log_temp_table, const D& data,
                 const amrex::Real log_rhoy, const amrex::Real rhoy,
                 const amrex::Real log_temp, const amrex::Real temp,
                 weak_rate_t& entries)
{

    // Get the table entries at this rhoy, temp

    // get the index into the table -- all components are
    // stored with the same indices

    const int irhoy_lo = interp_net::find_index_extrap(log_rhoy, log_rhoy_table);
    const int jtemp_lo = interp_net::find_index_extrap(log_temp, log_temp_table);

    get_entries<do_T_derivatives>(irhoy_lo, jtemp_lo,
                                  table_meta, log_rhoy_table, log_temp_table, data,
                                  log_rhoy, log_temp, entries);

    // Fill outputs: rate, d(rate)/dT, and (negative) neutrino loss
    // contribution to energy generation

    entries.rate = amrex::Math::exp10(entries.rate);
    if constexpr (do_T_derivatives) {
        entries.drate_dT = entries.rate * entries.drate_dT / temp;
        entries.drate_drhoye = entries.rate * entries.drate_drhoye / rhoy;
    } else {
        entries.drate_dT = 0.0;
        entries.drate_drhoye = 0.0;
    }
    entries.enu = -amrex::Math::exp10(entries.enu);
    if constexpr (do_T_derivatives) {
        entries.denu_dT = entries.enu * entries.denu_dT / temp;
        entries.denu_drhoye = entries.enu * entries.denu_drhoye / rhoy;
    } else {
        entries.denu_dT = 0.0;
        entries.denu_drhoye = 0.0;
    }
    entries.gamma = amrex::Math::exp10(entries.gamma);
}

#endif
