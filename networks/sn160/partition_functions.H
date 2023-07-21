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

        // find the first temperature element larger than t9
        auto upper = std::upper_bound(std::begin(temp_array), std::end(temp_array), t9);

        int idx = -1;
        if (upper != std::end(temp_array)) {
            // base index is one before the upper element
            // this will be -1 if t9 < smallest entry in the table
            idx = std::distance(temp_array, upper) - 1;
        }

        if (idx >= 0) {

            // construct the slope -- this is (log10(pf_{i+1}) - log10(pf_i)) / (T_{i+1} - T_i)

            Real slope = (pf_array[idx+1] - pf_array[idx]) / (temp_array[idx+1] - temp_array[idx]);

            // find the PF

            Real log10_pf = pf_array[idx] + slope * (t9 - temp_array[idx]);
            pf = std::pow(10.0_rt, log10_pf);

            // find the derivative (with respect to T, not T9)

            Real dpf_dT9 = pf * M_LN10 * slope;
            dpf_dT = dpf_dT9 / 1.e9_rt;

        } else {

            // T < the smallest T in the partition function table
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

    case N:
        spin = 2;
        break;

    case H1:
        spin = 2;
        break;

    case H2:
        spin = 3;
        break;

    case He3:
        spin = 2;
        break;

    case He4:
        spin = 1;
        break;

    case Li6:
        spin = 3;
        break;

    case Li7:
        spin = 4;
        break;

    case Be7:
        spin = 4;
        break;

    case Be9:
        spin = 4;
        break;

    case B8:
        spin = 5;
        break;

    case B10:
        spin = 7;
        break;

    case B11:
        spin = 4;
        break;

    case C12:
        spin = 1;
        break;

    case C13:
        spin = 2;
        break;

    case C14:
        spin = 1;
        break;

    case N13:
        spin = 2;
        break;

    case N14:
        spin = 3;
        break;

    case N15:
        spin = 2;
        break;

    case O14:
        spin = 1;
        break;

    case O15:
        spin = 2;
        break;

    case O16:
        spin = 1;
        break;

    case O17:
        spin = 6;
        break;

    case O18:
        spin = 1;
        break;

    case F17:
        spin = 6;
        break;

    case F18:
        spin = 3;
        break;

    case F19:
        spin = 2;
        break;

    case Ne18:
        spin = 1;
        break;

    case Ne19:
        spin = 2;
        break;

    case Ne20:
        spin = 1;
        break;

    case Ne21:
        spin = 4;
        break;

    case Ne22:
        spin = 1;
        break;

    case Na21:
        spin = 4;
        break;

    case Na22:
        spin = 7;
        break;

    case Na23:
        spin = 4;
        break;

    case Mg23:
        spin = 4;
        break;

    case Mg24:
        spin = 1;
        break;

    case Mg25:
        spin = 6;
        break;

    case Mg26:
        spin = 1;
        break;

    case Al25:
        spin = 6;
        break;

    case Al26:
        spin = 11;
        break;

    case Al27:
        spin = 6;
        break;

    case Si28:
        spin = 1;
        break;

    case Si29:
        spin = 2;
        break;

    case Si30:
        spin = 1;
        break;

    case Si31:
        spin = 4;
        break;

    case Si32:
        spin = 1;
        break;

    case P29:
        spin = 2;
        break;

    case P30:
        spin = 3;
        break;

    case P31:
        spin = 2;
        break;

    case P32:
        spin = 3;
        break;

    case P33:
        spin = 2;
        break;

    case S32:
        spin = 1;
        break;

    case S33:
        spin = 4;
        break;

    case S34:
        spin = 1;
        break;

    case S35:
        spin = 4;
        break;

    case S36:
        spin = 1;
        break;

    case Cl33:
        spin = 4;
        break;

    case Cl34:
        spin = 1;
        break;

    case Cl35:
        spin = 4;
        break;

    case Cl36:
        spin = 5;
        break;

    case Cl37:
        spin = 4;
        break;

    case Ar36:
        spin = 1;
        break;

    case Ar37:
        spin = 4;
        break;

    case Ar38:
        spin = 1;
        break;

    case Ar39:
        spin = 8;
        break;

    case Ar40:
        spin = 1;
        break;

    case K37:
        spin = 4;
        break;

    case K38:
        spin = 7;
        break;

    case K39:
        spin = 4;
        break;

    case K40:
        spin = 9;
        break;

    case K41:
        spin = 4;
        break;

    case Ca40:
        spin = 1;
        break;

    case Ca41:
        spin = 8;
        break;

    case Ca42:
        spin = 1;
        break;

    case Ca43:
        spin = 8;
        break;

    case Ca44:
        spin = 1;
        break;

    case Ca45:
        spin = 8;
        break;

    case Ca46:
        spin = 1;
        break;

    case Ca47:
        spin = 8;
        break;

    case Ca48:
        spin = 1;
        break;

    case Sc43:
        spin = 8;
        break;

    case Sc44:
        spin = 5;
        break;

    case Sc45:
        spin = 8;
        break;

    case Sc46:
        spin = 9;
        break;

    case Sc47:
        spin = 8;
        break;

    case Sc48:
        spin = 13;
        break;

    case Sc49:
        spin = 8;
        break;

    case Ti44:
        spin = 1;
        break;

    case Ti45:
        spin = 8;
        break;

    case Ti46:
        spin = 1;
        break;

    case Ti47:
        spin = 6;
        break;

    case Ti48:
        spin = 1;
        break;

    case Ti49:
        spin = 8;
        break;

    case Ti50:
        spin = 1;
        break;

    case Ti51:
        spin = 4;
        break;

    case V46:
        spin = 1;
        break;

    case V47:
        spin = 4;
        break;

    case V48:
        spin = 9;
        break;

    case V49:
        spin = 8;
        break;

    case V50:
        spin = 13;
        break;

    case V51:
        spin = 8;
        break;

    case V52:
        spin = 7;
        break;

    case Cr48:
        spin = 1;
        break;

    case Cr49:
        spin = 6;
        break;

    case Cr50:
        spin = 1;
        break;

    case Cr51:
        spin = 8;
        break;

    case Cr52:
        spin = 1;
        break;

    case Cr53:
        spin = 4;
        break;

    case Cr54:
        spin = 1;
        break;

    case Mn50:
        spin = 1;
        break;

    case Mn51:
        spin = 6;
        break;

    case Mn52:
        spin = 13;
        break;

    case Mn53:
        spin = 8;
        break;

    case Mn54:
        spin = 7;
        break;

    case Mn55:
        spin = 6;
        break;

    case Fe52:
        spin = 1;
        break;

    case Fe53:
        spin = 8;
        break;

    case Fe54:
        spin = 1;
        break;

    case Fe55:
        spin = 4;
        break;

    case Fe56:
        spin = 1;
        break;

    case Fe57:
        spin = 2;
        break;

    case Fe58:
        spin = 1;
        break;

    case Co54:
        spin = 1;
        break;

    case Co55:
        spin = 8;
        break;

    case Co56:
        spin = 9;
        break;

    case Co57:
        spin = 8;
        break;

    case Co58:
        spin = 5;
        break;

    case Co59:
        spin = 8;
        break;

    case Ni56:
        spin = 1;
        break;

    case Ni57:
        spin = 4;
        break;

    case Ni58:
        spin = 1;
        break;

    case Ni59:
        spin = 4;
        break;

    case Ni60:
        spin = 1;
        break;

    case Ni61:
        spin = 4;
        break;

    case Ni62:
        spin = 1;
        break;

    case Ni63:
        spin = 2;
        break;

    case Ni64:
        spin = 1;
        break;

    case Cu57:
        spin = 4;
        break;

    case Cu58:
        spin = 3;
        break;

    case Cu59:
        spin = 4;
        break;

    case Cu60:
        spin = 5;
        break;

    case Cu61:
        spin = 4;
        break;

    case Cu62:
        spin = 3;
        break;

    case Cu63:
        spin = 4;
        break;

    case Cu64:
        spin = 3;
        break;

    case Cu65:
        spin = 4;
        break;

    case Zn59:
        spin = 4;
        break;

    case Zn60:
        spin = 1;
        break;

    case Zn61:
        spin = 4;
        break;

    case Zn62:
        spin = 1;
        break;

    case Zn63:
        spin = 4;
        break;

    case Zn64:
        spin = 1;
        break;

    case Zn65:
        spin = 6;
        break;

    case Zn66:
        spin = 1;
        break;

    case Ga62:
        spin = 1;
        break;

    case Ga63:
        spin = 4;
        break;

    case Ge64:
        spin = 1;
        break;


    }

    return spin;

}


#endif
