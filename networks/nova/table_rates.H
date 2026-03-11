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

const int num_tables = 2;

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

// we add an additional index, k_index_dlogr_dlogt used for computing
// the derivative of Log(rate) with respect of Log(temperature) by
// using the table values. It isn't an index into the table but into
// the 'entries' array. Is important to mention that although we
// compute dlogr/dlogT as the computed quantity in 'entries', we
// pursue ultimately dr/dt as the final desired quantity to be
// computed for this index.

constexpr int k_index_dlogr_dlogt = num_vars + 1;
constexpr int add_vars = 1;  // 1 Additional Var in entries


namespace rate_tables
{
    // F17 --> O17
    inline AMREX_GPU_MANAGED table_t j_F17_O17_meta{.ntemp=12, .nrhoy=11, .nvars=3, .nheader=5};
    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_F17_O17_rhoy{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0};
    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_F17_O17_temp{7.0, 8.0, 8.30103, 8.60206, 8.845098, 9.0, 9.176091, 9.30103, 9.477121, 9.69897, 10.0, 10.477121};
    // Array3D is column-major (Fortran-ordering).  T varies fastest, then rho Ye, then the component
    inline AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, num_vars>
         j_F17_O17_data{
-1.98599,       -1.985996,       -1.985997,       -1.985998,       -1.985991,       -1.985811,       -1.984179,       -1.977601,       -1.918501,       -1.559899,       -0.462121,        1.718247,  // RATE
-1.985903,       -1.985964,       -1.985971,       -1.985976,       -1.985976,       -1.985801,       -1.984173,       -1.977571,       -1.918501,       -1.559899,       -0.462121,        1.718247,
-1.985273,       -1.985641,       -1.985715,       -1.985762,       -1.985784,       -1.985681,       -1.984074,       -1.977479,       -1.918324,       -1.559899,       -0.462121,        1.718247,
-1.982066,       -1.982755,        -1.98325,       -1.983647,       -1.983854,       -1.983902,       -1.982956,       -1.976495,       -1.917434,       -1.559232,       -0.462121,        1.718247,
-1.962293,       -1.962505,       -1.963077,       -1.964279,       -1.965373,       -1.965787,        -1.96628,       -1.962645,       -1.906863,       -1.553549,       -0.461157,        1.718247,
-1.822514,       -1.822514,       -1.822514,       -1.821886,       -1.820939,       -1.819671,       -1.816492,       -1.812226,       -1.783616,       -1.497551,       -0.452441,        1.719242,
-1.201147,       -1.201147,       -1.200311,       -1.199474,       -1.196963,       -1.193611,       -1.183682,       -1.172312,       -1.140726,       -1.040908,       -0.370812,        1.722232,
-0.075628,       -0.075628,       -0.075628,       -0.075628,       -0.073653,       -0.071677,        -0.06675,       -0.060857,       -0.041154,        0.012185,        0.215715,        1.752107,
1.273239,        1.273239,        1.273239,        1.273239,        1.273239,        1.274239,        1.276237,        1.280233,        1.289225,        1.314208,        1.405246,        2.033211,
2.796007,        2.796007,        2.796007,        2.796007,        2.796007,        2.796007,        2.797007,        2.799007,        2.802007,        2.812007,        2.843009,        3.076122,
4.404,           4.404,           4.404,           4.404,           4.404,           4.404,           4.404,           4.405,           4.406,           4.408,           4.416,        4.486005,
-7.78129,        -7.78129,        -7.78129,        -7.78129,        -7.78129,        -7.78129,        -7.77529,        -7.74729,        -7.56429,        -6.87729,        -5.49729,        -2.92229,  // NU
-7.78129,        -7.78129,        -7.78129,        -7.78129,        -7.78129,        -7.78129,        -7.77529,        -7.74729,        -7.56429,        -6.87729,        -5.49729,        -2.92229,
-7.77929,        -7.78029,        -7.78029,        -7.78029,        -7.78029,        -7.78029,        -7.77429,        -7.74729,        -7.56329,        -6.87729,        -5.49729,        -2.92229,
-7.77029,        -7.77229,        -7.77329,        -7.77429,        -7.77529,        -7.77529,        -7.77129,        -7.74429,        -7.56129,        -6.87629,        -5.49729,        -2.92229,
-7.71829,        -7.71829,        -7.71929,        -7.72229,        -7.72429,        -7.72429,        -7.72329,        -7.70629,        -7.53729,        -6.86929,        -5.49629,        -2.92229,
-7.41529,        -7.41529,        -7.41429,        -7.41329,        -7.40829,        -7.40229,        -7.38929,        -7.37129,        -7.29229,        -6.79729,        -5.48829,        -2.92129,
-6.53729,        -6.53629,        -6.53629,        -6.53429,        -6.53029,        -6.52229,        -6.50629,        -6.48329,        -6.42629,        -6.25829,        -5.40329,        -2.91829,
-5.25129,        -5.25129,        -5.25129,        -5.25029,        -5.24829,        -5.24529,        -5.23729,        -5.22629,        -5.19729,        -5.11129,        -4.79429,        -2.88729,
-3.71829,        -3.71829,        -3.71829,        -3.71829,        -3.71729,        -3.71629,        -3.71329,        -3.70929,        -3.69629,        -3.66029,        -3.51729,        -2.60029,
-1.96929,        -1.96929,        -1.96929,        -1.96929,        -1.96929,        -1.96829,        -1.96729,        -1.96529,        -1.95929,        -1.94629,        -1.89529,        -1.50629,
-0.07329,        -0.07329,        -0.07329,        -0.07329,        -0.07329,        -0.07329,        -0.07329,        -0.07229,        -0.07029,        -0.06629,        -0.05229,         0.07271,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,  // GAMMA
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
         };

    // O17 --> F17
    inline AMREX_GPU_MANAGED table_t j_O17_F17_meta{.ntemp=12, .nrhoy=11, .nvars=3, .nheader=5};
    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_O17_F17_rhoy{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0};
    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_O17_F17_temp{7.0, 8.0, 8.30103, 8.60206, 8.845098, 9.0, 9.176091, 9.30103, 9.477121, 9.69897, 10.0, 10.477121};
    // Array3D is column-major (Fortran-ordering).  T varies fastest, then rho Ye, then the component
    inline AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, num_vars>
         j_O17_F17_data{
-99.69897,       -99.69897,         -72.639,         -36.441,         -21.007,         -15.597,         -11.237,          -8.922,       -6.409967,       -4.069196,       -1.723001,        1.202117,  // RATE
-99.69897,       -99.69897,         -73.639,         -37.441,          -21.42,         -15.619,         -11.238,          -8.923,       -6.409967,       -4.069196,       -1.723001,        1.202117,
-99.69897,       -99.69897,         -74.645,         -38.443,         -22.368,         -15.825,         -11.253,          -8.926,       -6.409967,       -4.069196,       -1.723001,        1.202117,
-99.69897,       -99.69897,           -75.7,          -39.46,         -23.374,         -16.647,         -11.395,          -8.957,       -6.415967,       -4.070194,       -1.723997,        1.202117,
-99.69897,       -99.69897,         -77.203,         -40.631,         -24.441,         -17.678,         -12.128,          -9.246,       -6.468963,       -4.078179,       -1.723997,        1.202117,
-99.69897,       -99.69897,         -81.385,         -42.905,         -25.991,         -18.982,         -13.265,      -10.212996,       -6.911898,       -4.161015,       -1.732959,        1.202117,
-99.69897,       -99.69897,          -93.93,         -49.222,         -29.674,         -21.645,      -15.175994,      -11.777875,        -8.12459,       -4.770369,        -1.82157,        1.199117,
-99.69897,       -99.69897,       -99.69897,          -64.67,       -38.52992,      -27.872783,      -19.361149,      -14.944239,      -10.306988,       -6.238059,       -2.468687,        1.168123,
-99.69897,       -99.69897,       -99.69897,      -98.710944,      -57.848526,      -41.325098,      -28.297405,      -21.662021,      -14.845418,       -9.070587,       -4.075741,         0.87518,
-99.69897,       -99.69897,       -99.69897,       -99.69897,       -99.69897,      -71.263098,      -48.267405,      -36.650021,       -24.85878,      -15.119379,        -7.19457,        -0.32867,
-99.69897,       -99.69897,       -99.69897,       -99.69897,       -99.69897,       -99.69897,      -91.333914,       -68.95548,      -46.405418,      -28.066379,       -13.71257,        -2.65565,
-105.79529,      -105.79529,       -79.71329,       -43.20829,       -27.52229,       -21.94929,       -17.40029,       -14.94929,       -12.24129,        -9.65029,        -6.95829,        -3.50429,  // NU
-105.79529,      -105.79529,       -80.71429,       -44.20829,       -27.93529,       -21.97129,       -17.40229,       -14.95029,       -12.24129,        -9.65029,        -6.95829,        -3.50429,
-105.79529,      -105.79529,       -81.71929,       -45.21029,       -28.88429,       -22.17729,       -17.41629,       -14.95329,       -12.24229,        -9.65029,        -6.95829,        -3.50429,
-105.79529,      -105.79529,       -82.77429,       -46.22829,       -29.89029,       -22.99929,       -17.55829,       -14.98429,       -12.24729,        -9.65029,        -6.95829,        -3.50429,
-105.79529,      -105.79529,       -84.27729,       -47.39929,       -30.95629,       -24.03029,       -18.29129,       -15.27329,       -12.30129,        -9.65929,        -6.95929,        -3.50429,
-105.79529,      -105.79529,       -88.46029,       -49.67329,       -32.50629,       -25.33429,       -19.42829,       -16.23929,       -12.74329,        -9.74129,        -6.96829,        -3.50429,
-105.79529,      -105.79529,      -101.00429,       -55.98929,       -36.19029,       -27.99729,       -21.33929,       -17.80529,       -13.95529,       -10.35029,        -7.05629,        -3.50729,
-105.79529,      -105.79529,      -105.79529,       -71.43729,       -45.04529,       -34.22129,       -25.51129,       -20.95629,       -16.13629,       -11.83029,        -7.70729,        -3.53929,
-105.79529,      -105.79529,      -105.79529,      -105.71329,       -64.38629,       -47.70729,       -34.50029,       -27.73529,       -20.72829,       -14.69429,        -9.32329,        -3.83329,
-105.79529,      -105.79529,      -105.79529,      -105.79529,      -105.79529,       -77.64529,       -54.47029,       -42.72329,       -30.74129,       -20.74329,       -12.44229,        -5.03929,
-105.79529,      -105.79529,      -105.79529,      -105.79529,      -105.79529,      -105.79529,       -97.53629,       -75.02829,       -52.28829,       -33.69029,       -18.96029,        -7.36629,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,  // GAMMA
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
-100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,          -100.0,
         };

}


AMREX_INLINE AMREX_GPU_HOST_DEVICE
amrex::Real
evaluate_linear_1d(const amrex::Real fhi, const amrex::Real flo, const amrex::Real xhi, const amrex::Real xlo, const amrex::Real x)
{
    // This function is a 1-D linear interpolator, that keeps x constant to xlo or xhi, based
    // on the side, if x is outside [xlo, xhi] respectively.

    amrex::Real xx = amrex::Clamp(x, xlo, xhi);
    amrex::Real f = flo + (fhi - flo) * (xx - xlo) / (xhi - xlo);

    return f;
}

AMREX_INLINE AMREX_GPU_HOST_DEVICE
amrex::Real
evaluate_linear_2d(const amrex::Real fip1jp1, const amrex::Real fip1j, const amrex::Real fijp1, const amrex::Real fij,
          const amrex::Real xhi, const amrex::Real xlo, const amrex::Real yhi, const amrex::Real ylo,
          const amrex::Real x, const amrex::Real y)
{
    // This is the 2-D linear interpolator, as an extension of evaluate_linear_1d.

    amrex::Real f;
    amrex::Real dx = xhi - xlo;
    amrex::Real dxinv = 1.0_rt / dx;
    amrex::Real dy = yhi - ylo;
    amrex::Real dyinv = 1.0_rt / dy;

    amrex::Real E =  fij;
    amrex::Real C = (fijp1 - fij) * dyinv;
    amrex::Real B = (fip1j - fij) * dxinv;
    amrex::Real A = (fip1jp1 - B * dx - C * dy - E) * (dxinv * dyinv);

    amrex::Real xx = amrex::Clamp(x, xlo, xhi);
    amrex::Real yy = amrex::Clamp(y, ylo, yhi);

    f =  A * (xx - xlo) * (yy - ylo) +
         B * (xx - xlo) +
         C * (yy - ylo) +
         E;

    return f;
}


template<typename R, typename T, typename D>
AMREX_INLINE AMREX_GPU_HOST_DEVICE
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
AMREX_INLINE AMREX_GPU_HOST_DEVICE
amrex::Real
evaluate_dr_dtemp(const int irhoy_lo, const int jtemp_lo,
                  const table_t& table_meta, const R& log_rhoy_table, const T& log_temp_table, const D& data,
                  const amrex::Real log_rhoy, const amrex::Real log_temp)
{
    // The main objective of this function is compute dlogr_dlogt.

    const int irhoy_hi = irhoy_lo + 1;
    const int jtemp_hi = jtemp_lo + 1;

    amrex::Real dlogr_dlogt;

    //Now we compute the forward finite difference on the boundary

    if ((jtemp_lo - 1 < 1) || (jtemp_hi + 1 > table_meta.ntemp)) {

        // In this case we are in the boundaries of the table.
        // At the boundary, we compute the forward-j finite difference
        // to compute dlogr_dlogt_i and dlogr_dlogt_ip1, using the
        // following stencil:
        //
        //
        //             fijp1-----------fip1jp1
        //              |                 |
        //              |                 |
        //              |                 |
        //              |                 |
        //              |                 |
        //              |                 |
        //              |                 |
        //             fij-------------fip1j
        //
        // with the following result:
        //
        //            dlogr_dlogt_i --------dlogr_dlogt_ip1
        //
        // Finally, we perform a 1d-linear interpolation between dlogr_dlogt_ip1
        // and dlogr_dlogt_i to compute dlogr_dlogt

        const amrex::Real log_rhoy_lo = log_rhoy_table(irhoy_lo);
        const amrex::Real log_rhoy_hi = log_rhoy_table(irhoy_hi);

        const amrex::Real log_temp_lo = log_temp_table(jtemp_lo);
        const amrex::Real log_temp_hi = log_temp_table(jtemp_hi);

        const amrex::Real fij     = data(jtemp_lo, irhoy_lo, jtab_rate);
        const amrex::Real fip1j   = data(jtemp_lo, irhoy_hi, jtab_rate);
        const amrex::Real fijp1   = data(jtemp_hi, irhoy_lo, jtab_rate);
        const amrex::Real fip1jp1 = data(jtemp_hi, irhoy_hi, jtab_rate);

        const amrex::Real dlogr_dlogt_i   = (fijp1 - fij) / (log_temp_hi - log_temp_lo);
        const amrex::Real dlogr_dlogt_ip1 = (fip1jp1 - fip1j) / (log_temp_hi - log_temp_lo);

        if ((log_temp < log_temp_lo) || (log_temp > log_temp_hi)) {
            dlogr_dlogt = 0.0_rt;
        } else {
            dlogr_dlogt = evaluate_linear_1d(dlogr_dlogt_ip1, dlogr_dlogt_i, log_rhoy_hi, log_rhoy_lo, log_rhoy);
        }

    } else {

        // In this case, we use a bigger stencil to reconstruct the
        // temperature derivatives in the j and j+1 temperature positions,
        // using the cetral-j finite differences:
        //
        //              fijp2 ------------fip1jp2
        //                |                 |
        //                |                 |
        //                |                 |
        //                |                 |
        //                |                 |
        //                |                 |
        //                |                 |
        //               fijp1------------fip1jp1
        //                |                 |
        //                |                 |
        //                |                 |
        //                |                 |
        //                |                 |
        //                |                 |
        //                |                 |
        //               fij------------- fip1j
        //                |                 |
        //                |                 |
        //                |                 |
        //                |                 |
        //                |                 |
        //                |                 |
        //                |                 |
        //               fijm1------------fip1jm1
        //
        // with the following result:
        //
        //
        //            dr_dt_ijp1 --------dr_dt_ip1jp1
        //                |                 |
        //                |                 |
        //                |                 |
        //                |                 |
        //                |                 |
        //                |                 |
        //                |                 |
        //            dr_dt_ij-----------dr_dt_ip1j
        //
        // Finally, we perform a 2d-linear interpolation to
        // compute dlogr_dlogt.

        const amrex::Real log_temp_jm1  = log_temp_table(jtemp_lo-1);
        const amrex::Real log_temp_j    = log_temp_table(jtemp_lo);
        const amrex::Real log_temp_jp1  = log_temp_table(jtemp_hi);
        const amrex::Real log_temp_jp2  = log_temp_table(jtemp_hi+1);

        const amrex::Real log_rhoy_lo = log_rhoy_table(irhoy_lo);
        const amrex::Real log_rhoy_hi = log_rhoy_table(irhoy_hi);

        const amrex::Real fijm1   = data(jtemp_lo-1, irhoy_lo, jtab_rate);
        const amrex::Real fij     = data(jtemp_lo, irhoy_lo, jtab_rate);
        const amrex::Real fijp1   = data(jtemp_hi, irhoy_lo, jtab_rate);
        const amrex::Real fijp2   = data(jtemp_hi+1, irhoy_lo, jtab_rate);

        const amrex::Real fip1jm1 = data(jtemp_lo-1, irhoy_hi, jtab_rate);
        const amrex::Real fip1j   = data(jtemp_lo, irhoy_hi, jtab_rate);
        const amrex::Real fip1jp1 = data(jtemp_hi, irhoy_hi, jtab_rate);
        const amrex::Real fip1jp2 = data(jtemp_hi+1, irhoy_hi, jtab_rate);

        const amrex::Real dlogTinv1 = 1.0_rt / (log_temp_jp1 - log_temp_jm1);
        const amrex::Real dlogTinv2 = 1.0_rt / (log_temp_jp2 - log_temp_j);
        const amrex::Real dlogr_dlogt_ij     = (fijp1 - fijm1) * dlogTinv1;
        const amrex::Real dlogr_dlogt_ijp1   = (fijp2 - fij) * dlogTinv2;
        const amrex::Real dlogr_dlogt_ip1j   = (fip1jp1 - fip1jm1) * dlogTinv1;
        const amrex::Real dlogr_dlogt_ip1jp1 = (fip1jp2 - fip1j) * dlogTinv2;

        dlogr_dlogt = evaluate_linear_2d(dlogr_dlogt_ip1jp1,  dlogr_dlogt_ip1j, dlogr_dlogt_ijp1, dlogr_dlogt_ij,
                                         log_rhoy_hi, log_rhoy_lo, log_temp_jp1, log_temp_j,
                                         log_rhoy, log_temp);

    }
    return dlogr_dlogt;
}


template <typename R, typename T, typename D>
AMREX_INLINE AMREX_GPU_HOST_DEVICE
void
get_entries(const int irhoy_lo, const int jtemp_lo,
            const table_t& table_meta, const R& log_rhoy_table, const T& log_temp_table, const D& data,
            const amrex::Real log_rhoy, const amrex::Real log_temp, amrex::Array1D<amrex::Real, 1, num_vars+add_vars>& entries)
{

    // only interpolate the quantities we currently use

    entries(jtab_rate) = evaluate_vars(irhoy_lo, jtemp_lo,
                                       log_rhoy_table, log_temp_table, data,
                                       log_rhoy, log_temp, jtab_rate);

    entries(jtab_nuloss) = evaluate_vars(irhoy_lo, jtemp_lo,
                                         log_rhoy_table, log_temp_table, data,
                                         log_rhoy, log_temp, jtab_nuloss);

    entries(jtab_gamma) = evaluate_vars(irhoy_lo, jtemp_lo,
                                        log_rhoy_table, log_temp_table, data,
                                        log_rhoy, log_temp, jtab_gamma);

    entries(k_index_dlogr_dlogt) =
        evaluate_dr_dtemp(irhoy_lo, jtemp_lo,
                          table_meta, log_rhoy_table, log_temp_table, data,
                          log_rhoy, log_temp);
}

template <typename R, typename T, typename D>
AMREX_INLINE AMREX_GPU_HOST_DEVICE
void
tabular_evaluate(const table_t& table_meta,
                 const R& log_rhoy_table, const T& log_temp_table, const D& data,
                 const amrex::Real log_rhoy, const amrex::Real log_temp, const amrex::Real temp,
                 amrex::Real& rate, amrex::Real& drate_dt, amrex::Real& edot_nu, amrex::Real& edot_gamma)
{
    amrex::Array1D<amrex::Real, 1, num_vars+add_vars> entries;

    // Get the table entries at this rhoy, temp

    // get the index into the table -- all components are
    // stored with the same indices

    const int irhoy_lo = interp_net::find_index_extrap(log_rhoy, log_rhoy_table);
    const int jtemp_lo = interp_net::find_index_extrap(log_temp, log_temp_table);

    get_entries(irhoy_lo, jtemp_lo,
                table_meta, log_rhoy_table, log_temp_table, data,
                log_rhoy, log_temp, entries);

    // Fill outputs: rate, d(rate)/d(temperature), and
    // (negative) neutrino loss contribution to energy generation

    rate       = amrex::Math::exp10(entries(jtab_rate));
    drate_dt   = rate * entries(k_index_dlogr_dlogt) / temp;
    edot_nu    = -amrex::Math::exp10(entries(jtab_nuloss));
    edot_gamma = amrex::Math::exp10(entries(jtab_gamma));
}

#endif
