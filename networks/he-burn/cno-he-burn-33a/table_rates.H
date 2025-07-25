#ifndef TABLE_RATES_H
#define TABLE_RATES_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <AMReX_Array.H>

using namespace amrex;

void init_tabular();

// Table is expected to be in terms of dens*ye and temp (logarithmic, cgs units)
// Table energy units are expected in terms of ergs

// all tables are expected to have columns:
// Log(rhoY)     Log(T)   mu    dQ    Vs    Log(e-cap-rate)   Log(nu-energy-loss)  Log(gamma-energy)
// Log(g/cm^3)   Log(K)   erg   erg   erg   Log(1/s)          Log(erg/s)           Log(erg/s)
//

const int num_tables = 10;

enum TableVars : std::uint8_t
{
    jtab_mu      = 1,
    jtab_dq      = 2,
    jtab_vs      = 3,
    jtab_rate    = 4,
    jtab_nuloss  = 5,
    jtab_gamma   = 6,
    num_vars = jtab_gamma
};


struct table_t
{
    int ntemp;
    int nrhoy;
    int nvars;
    int nheader;
};

// we add a 7th index, k_index_dlogr_dlogt used for computing the derivative
// of Log(rate) with respect of Log(temperature) by using the table
// values. It isn't an index into the table but into the 'entries'
// array. Is important to mention that although we compute dlogr/dlogT is
// the computed quantity in 'entries', we pursue ultimately
// dr/dt as the final desired quantity to be computed for this index.

const int  k_index_dlogr_dlogt  = 7;
const int add_vars              = 1;  // 1 Additional Var in entries


namespace rate_tables
{
    extern AMREX_GPU_MANAGED table_t j_F17_O17_meta;
    extern AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_F17_O17_data;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_F17_O17_rhoy;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_F17_O17_temp;

    extern AMREX_GPU_MANAGED table_t j_O17_F17_meta;
    extern AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_O17_F17_data;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_O17_F17_rhoy;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_O17_F17_temp;

    extern AMREX_GPU_MANAGED table_t j_F18_Ne18_meta;
    extern AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_F18_Ne18_data;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_F18_Ne18_rhoy;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_F18_Ne18_temp;

    extern AMREX_GPU_MANAGED table_t j_F18_O18_meta;
    extern AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_F18_O18_data;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_F18_O18_rhoy;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_F18_O18_temp;

    extern AMREX_GPU_MANAGED table_t j_Ne18_F18_meta;
    extern AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Ne18_F18_data;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ne18_F18_rhoy;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Ne18_F18_temp;

    extern AMREX_GPU_MANAGED table_t j_O18_F18_meta;
    extern AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_O18_F18_data;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_O18_F18_rhoy;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_O18_F18_temp;

    extern AMREX_GPU_MANAGED table_t j_F19_Ne19_meta;
    extern AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_F19_Ne19_data;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_F19_Ne19_rhoy;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_F19_Ne19_temp;

    extern AMREX_GPU_MANAGED table_t j_Ne19_F19_meta;
    extern AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Ne19_F19_data;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ne19_F19_rhoy;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Ne19_F19_temp;

    extern AMREX_GPU_MANAGED table_t j_Mg22_Na22_meta;
    extern AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Mg22_Na22_data;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mg22_Na22_rhoy;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Mg22_Na22_temp;

    extern AMREX_GPU_MANAGED table_t j_Na22_Mg22_meta;
    extern AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Na22_Mg22_data;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Na22_Mg22_rhoy;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Na22_Mg22_temp;

}

template <typename R, typename T, typename D>
void init_tab_info(const table_t& tf, const std::string& file, R& log_rhoy_table, T& log_temp_table, D& data)
{
    // This function initializes the selected tabular-rate tables. From the tables we are interested
    // on the rate, neutrino-energy-loss and the gamma-energy entries.

    std::ifstream table;
    table.open(file);

    if (!table.is_open()) {
        // the table was not present or we could not open it; abort
        amrex::Error("table could not be opened");
    }

    std::string line;

    // read and skip over the header

    for (int i = 0; i < tf.nheader; ++i) {
        std::getline(table, line);
    }

    // now the data -- there are 2 extra columns, for log_temp and log_rhoy

    for (int j = 1; j <= tf.nrhoy; ++j) {
        for (int i = 1; i <= tf.ntemp; ++i) {
            std::getline(table, line);
            if (line.empty()) {
                amrex::Error("Error reading table data");
            }

            std::istringstream sdata(line);

            sdata >> log_rhoy_table(j) >> log_temp_table(i);

            for (int n = 1; n <= tf.nvars; ++n) {
                sdata >> data(i, j, n);
            }
        }
    }
    table.close();
}


template <typename V>
AMREX_INLINE AMREX_GPU_HOST_DEVICE
int vector_index_lu(const int vlen, const V& vector, const amrex::Real fvar)
{

    // Returns the greatest index of vector for which vector(index) < fvar.
    // Return 1 if fvar < vector(1)
    // Return size(vector)-1 if fvar > vector(size(vector))
    // The interval [index, index+1] brackets fvar for fvar within the range of vector.

    int index;

    if (fvar < vector(1)) {
        index = 1;
    } else if (fvar > vector(vlen)) {
        index = vlen - 1;
    } else {
        int nup = vlen;
        int ndn = 1;
        for (int i = 1; i <= vlen; ++i) {
            int j = ndn + (nup - ndn)/2;
            if (fvar < vector(j)) {
                nup = j;
            } else {
                ndn = j;
            }
            if ((nup - ndn) == 1) {
                break;
            }
        }
        index = ndn;
    }
    return index;
}


AMREX_INLINE AMREX_GPU_HOST_DEVICE
amrex::Real
evaluate_linear_1d(const amrex::Real fhi, const amrex::Real flo, const amrex::Real xhi, const amrex::Real xlo, const amrex::Real x)
{
    // This function is a 1-D linear interpolator, that keeps x constant to xlo or xhi, based
    // on the side, if x is outside [xlo, xhi] respectively.

    amrex::Real xx = Clamp(x, xlo, xhi);
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
    amrex::Real dy = yhi - ylo;

    amrex::Real E =  fij;
    amrex::Real C = (fijp1 - fij) / dy;
    amrex::Real B = (fip1j - fij) / dx;
    amrex::Real A = (fip1jp1 - B * dx - C * dy - E) / (dx * dy);

    amrex::Real xx = Clamp(x, xlo, xhi);
    amrex::Real yy = Clamp(y, ylo, yhi);

    f =  A * (xx - xlo) * (yy - ylo) +
         B * (xx - xlo) +
         C * (yy - ylo) +
         E;

    return f;
}


template<typename R, typename T, typename D>
AMREX_INLINE AMREX_GPU_HOST_DEVICE
amrex::Real
evaluate_vars(const table_t& table_meta, const R& log_rhoy_table, const T& log_temp_table, const D& data,
                    const amrex::Real log_rhoy, const amrex::Real log_temp, const int component)
{
    // This function evaluates the 2-D interpolator, for several pairs of rho_ye and temperature.

    int jtemp_lo = vector_index_lu(table_meta.ntemp, log_temp_table, log_temp);
    int jtemp_hi = jtemp_lo + 1;

    int irhoy_lo = vector_index_lu(table_meta.nrhoy, log_rhoy_table, log_rhoy);
    int irhoy_hi = irhoy_lo + 1;

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
evaluate_dr_dtemp(const table_t& table_meta, const R& log_rhoy_table, const T& log_temp_table, const D& data,
                  const amrex::Real log_rhoy, const amrex::Real log_temp)
{
    // The main objective of this function is compute dlogr_dlogt.

    int irhoy_lo = vector_index_lu(table_meta.nrhoy, log_rhoy_table, log_rhoy);
    int irhoy_hi = irhoy_lo + 1;

    int jtemp_lo = vector_index_lu(table_meta.ntemp, log_temp_table, log_temp);
    int jtemp_hi = jtemp_lo + 1;

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

        amrex::Real log_rhoy_lo = log_rhoy_table(irhoy_lo);
        amrex::Real log_rhoy_hi = log_rhoy_table(irhoy_hi);

        amrex::Real log_temp_lo = log_temp_table(jtemp_lo);
        amrex::Real log_temp_hi = log_temp_table(jtemp_hi);

        amrex::Real fij     = data(jtemp_lo, irhoy_lo, jtab_rate);
        amrex::Real fip1j   = data(jtemp_lo, irhoy_hi, jtab_rate);
        amrex::Real fijp1   = data(jtemp_hi, irhoy_lo, jtab_rate);
        amrex::Real fip1jp1 = data(jtemp_hi, irhoy_hi, jtab_rate);

        amrex::Real dlogr_dlogt_i   = (fijp1 - fij) / (log_temp_hi - log_temp_lo);
        amrex::Real dlogr_dlogt_ip1 = (fip1jp1 - fip1j) / (log_temp_hi - log_temp_lo);

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

        amrex::Real log_temp_jm1  = log_temp_table(jtemp_lo-1);
        amrex::Real log_temp_j    = log_temp_table(jtemp_lo);
        amrex::Real log_temp_jp1  = log_temp_table(jtemp_hi);
        amrex::Real log_temp_jp2  = log_temp_table(jtemp_hi+1);

        amrex::Real log_rhoy_lo = log_rhoy_table(irhoy_lo);
        amrex::Real log_rhoy_hi = log_rhoy_table(irhoy_hi);

        amrex::Real fijm1   = data(jtemp_lo-1, irhoy_lo, jtab_rate);
        amrex::Real fij     = data(jtemp_lo, irhoy_lo, jtab_rate);
        amrex::Real fijp1   = data(jtemp_hi, irhoy_lo, jtab_rate);
        amrex::Real fijp2   = data(jtemp_hi+1, irhoy_lo, jtab_rate);

        amrex::Real fip1jm1 = data(jtemp_lo-1, irhoy_hi, jtab_rate);
        amrex::Real fip1j   = data(jtemp_lo, irhoy_hi, jtab_rate);
        amrex::Real fip1jp1 = data(jtemp_hi, irhoy_hi, jtab_rate);
        amrex::Real fip1jp2 = data(jtemp_hi+1, irhoy_hi, jtab_rate);

        amrex::Real dlogr_dlogt_ij     = (fijp1 - fijm1)/(log_temp_jp1 - log_temp_jm1);
        amrex::Real dlogr_dlogt_ijp1   = (fijp2 - fij)/(log_temp_jp2 - log_temp_j);
        amrex::Real dlogr_dlogt_ip1j   = (fip1jp1 - fip1jm1)/(log_temp_jp1 - log_temp_jm1);
        amrex::Real dlogr_dlogt_ip1jp1 = (fip1jp2 - fip1j)/(log_temp_jp2 - log_temp_j);

        dlogr_dlogt = evaluate_linear_2d(dlogr_dlogt_ip1jp1,  dlogr_dlogt_ip1j, dlogr_dlogt_ijp1, dlogr_dlogt_ij,
                                         log_rhoy_hi, log_rhoy_lo, log_temp_jp1, log_temp_j,
                                         log_rhoy, log_temp);

    }
    return dlogr_dlogt;
}


template <typename R, typename T, typename D>
AMREX_INLINE AMREX_GPU_HOST_DEVICE
void
get_entries(const table_t& table_meta, const R& log_rhoy_table, const T& log_temp_table, const D& data,
            const amrex::Real log_rhoy, const amrex::Real log_temp, amrex::Array1D<amrex::Real, 1, num_vars+1>& entries)
{
    for (int ivar = 1; ivar <= num_vars; ivar++) {
        entries(ivar) = evaluate_vars(table_meta, log_rhoy_table, log_temp_table, data,
                                         log_rhoy, log_temp, ivar);
    }

    entries(k_index_dlogr_dlogt)  = evaluate_dr_dtemp(table_meta, log_rhoy_table, log_temp_table, data,
                                             log_rhoy, log_temp);
}

template <typename R, typename T, typename D>
AMREX_INLINE AMREX_GPU_HOST_DEVICE
void
tabular_evaluate(const table_t& table_meta,
                 const R& log_rhoy_table, const T& log_temp_table, const D& data,
                 const amrex::Real rhoy, const amrex::Real temp,
                 amrex::Real& rate, amrex::Real& drate_dt, amrex::Real& edot_nu, amrex::Real& edot_gamma)
{
    amrex::Array1D<amrex::Real, 1, num_vars+1> entries;

    // Get the table entries at this rhoy, temp

    amrex::Real log_rhoy = std::log10(rhoy);
    amrex::Real log_temp = std::log10(temp);

    get_entries(table_meta, log_rhoy_table, log_temp_table, data,
                log_rhoy, log_temp, entries);

    // Fill outputs: rate, d(rate)/d(temperature), and
    // (negative) neutrino loss contribution to energy generation

    rate       = std::pow(10.0_rt, entries(jtab_rate));
    drate_dt   = rate * entries(k_index_dlogr_dlogt) / temp;
    edot_nu    = -std::pow(10.0_rt, entries(jtab_nuloss));
    edot_gamma = std::pow(10.0_rt, entries(jtab_gamma));
}

#endif
