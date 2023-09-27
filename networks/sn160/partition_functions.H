#ifndef PARTITION_FUNCTIONS_H
#define PARTITION_FUNCTIONS_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <tfactors.H>
#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex;
using namespace Species;

namespace part_fun {



    // interpolation routine

    template <int npts>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void interpolate_pf(const Real t9, const Real (&temp_array)[npts], const Real (&pf_array)[npts],
                        Real& pf, Real& dpf_dT) {

        if (t9 >= temp_array[0] && t9 < temp_array[npts-1]) {

            // find the largest temperature element <= t9 using a binary search

            int left = 0;
            int right = npts;

            while (left < right) {
                int mid = (left + right) / 2;
                if (temp_array[mid] > t9) {
                    right = mid;
                } else {
                    left = mid + 1;
                }
            }

            const int idx = right - 1;

            // now we have temp_array[idx] <= t9 < temp_array[idx+1]

            // construct the slope -- this is (log10(pf_{i+1}) - log10(pf_i)) / (T_{i+1} - T_i)

            Real slope = (pf_array[idx+1] - pf_array[idx]) / (temp_array[idx+1] - temp_array[idx]);

            // find the PF

            Real log10_pf = pf_array[idx] + slope * (t9 - temp_array[idx]);
            pf = std::pow(10.0_rt, log10_pf);

            // find the derivative (with respect to T, not T9)

            Real dpf_dT9 = pf * M_LN10 * slope;
            dpf_dT = dpf_dT9 / 1.e9_rt;

        } else {

            // T < the smallest T or >= the largest T in the partition function table
            pf = 1.0;
            dpf_dT = 0.0;

        }

    }

}

// main interface

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_partition_function(const int inuc, [[maybe_unused]] const tf_t& tfactors,
                            Real& pf, Real& dpf_dT) {

    // inuc is the 1-based index for the species

    switch (inuc) {


    default:

        pf = 1.0_rt;
        dpf_dT = 0.0_rt;

    }

}

// spins

AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Real get_spin_state(const int inuc) {

    Real spin = -1.0;

    switch (inuc) {

    case He4:
    case C12:
    case C14:
    case O14:
    case O16:
    case O18:
    case Ne18:
    case Ne20:
    case Ne22:
    case Mg24:
    case Mg26:
    case Si28:
    case Si30:
    case Si32:
    case S32:
    case S34:
    case S36:
    case Cl34:
    case Ar36:
    case Ar38:
    case Ar40:
    case Ca40:
    case Ca42:
    case Ca44:
    case Ca46:
    case Ca48:
    case Ti44:
    case Ti46:
    case Ti48:
    case Ti50:
    case V46:
    case Cr48:
    case Cr50:
    case Cr52:
    case Cr54:
    case Mn50:
    case Fe52:
    case Fe54:
    case Fe56:
    case Fe58:
    case Co54:
    case Ni56:
    case Ni58:
    case Ni60:
    case Ni62:
    case Ni64:
    case Zn60:
    case Zn62:
    case Zn64:
    case Zn66:
    case Ga62:
    case Ge64:
        spin = 1;
        break;

    case N:
    case H1:
    case He3:
    case C13:
    case N13:
    case N15:
    case O15:
    case F19:
    case Ne19:
    case Si29:
    case P29:
    case P31:
    case P33:
    case Fe57:
    case Ni63:
        spin = 2;
        break;

    case H2:
    case Li6:
    case N14:
    case F18:
    case P30:
    case P32:
    case Cu58:
    case Cu62:
    case Cu64:
        spin = 3;
        break;

    case Li7:
    case Be7:
    case Be9:
    case B11:
    case Ne21:
    case Na21:
    case Na23:
    case Mg23:
    case Si31:
    case S33:
    case S35:
    case Cl33:
    case Cl35:
    case Cl37:
    case Ar37:
    case K37:
    case K39:
    case K41:
    case Ti51:
    case V47:
    case Cr53:
    case Fe55:
    case Ni57:
    case Ni59:
    case Ni61:
    case Cu57:
    case Cu59:
    case Cu61:
    case Cu63:
    case Cu65:
    case Zn59:
    case Zn61:
    case Zn63:
    case Ga63:
        spin = 4;
        break;

    case B8:
    case Cl36:
    case Sc44:
    case Co58:
    case Cu60:
        spin = 5;
        break;

    case O17:
    case F17:
    case Mg25:
    case Al25:
    case Al27:
    case Ti47:
    case Cr49:
    case Mn51:
    case Mn55:
    case Zn65:
        spin = 6;
        break;

    case B10:
    case Na22:
    case K38:
    case V52:
    case Mn54:
        spin = 7;
        break;

    case Ar39:
    case Ca41:
    case Ca43:
    case Ca45:
    case Ca47:
    case Sc43:
    case Sc45:
    case Sc47:
    case Sc49:
    case Ti45:
    case Ti49:
    case V49:
    case V51:
    case Cr51:
    case Mn53:
    case Fe53:
    case Co55:
    case Co57:
    case Co59:
        spin = 8;
        break;

    case K40:
    case Sc46:
    case V48:
    case Co56:
        spin = 9;
        break;

    case Al26:
        spin = 11;
        break;

    case Sc48:
    case V50:
    case Mn52:
        spin = 13;
        break;


    }

    return spin;

}


#endif
