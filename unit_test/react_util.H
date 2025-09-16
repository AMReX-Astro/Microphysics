#ifndef REACT_UTIL_H
#define REACT_UTIL_H

#include <network.H>
#include <extern_parameters.H>

using namespace amrex::literals;

struct init_t {
    int nprim;
    int is1;
    int is2;
    int is3;

    int n1;
    int n2;
    int n3;

    amrex::Real Xp_min;
    amrex::Real Xp_max;
};

AMREX_INLINE
init_t setup_composition(const int nz) {

    // get the primary species indices
    init_t comp_data;

    comp_data.nprim = 0;

    // we absolutely require primary_species_1 to be defined
    comp_data.is1 = network_spec_index(unit_test_rp::primary_species_1);
    if (comp_data.is1 >= 0) {
        comp_data.nprim++;
    } else {
        amrex::Error("Error: invalid primary_species_1");
    }

    // we'll check the next species, but if it is not valid,
    // we'll just move on
    comp_data.is2 = network_spec_index(unit_test_rp::primary_species_2);
    if (comp_data.is2 >= 0) {
        comp_data.nprim++;

        // only consider primary_species_3 if primary_species_2
        // was defined
        comp_data.is3 = network_spec_index(unit_test_rp::primary_species_3);
        if (comp_data.is3 >= 0) {
            comp_data.nprim++;
        }
    }

    if (comp_data.nprim == 0) {
        amrex::Error("ERROR: no primary species set");
    }

    // figure out how many zones to allocate to the each of the primary
    // species and the extrema for the primary species

    if (comp_data.nprim == 1) {
        comp_data.n1 = nz;
        comp_data.n2 = 0;
        comp_data.n3 = 0;

        comp_data.Xp_min = 0.2_rt;
        comp_data.Xp_max = 0.9_rt;

    } else if (comp_data.nprim == 2) {
        comp_data.n1 = nz/2;
        comp_data.n2 = nz - comp_data.n1;
        comp_data.n3 = 0;

        comp_data.Xp_min = 0.2_rt;
        comp_data.Xp_max = 0.8_rt;

    } else if (comp_data.nprim == 3) {
        comp_data.n1 = nz/3;
        comp_data.n2 = nz/3;
        comp_data.n3 = nz - comp_data.n1 - comp_data.n2;

        comp_data.Xp_min = 0.2_rt;
        comp_data.Xp_max = 0.7_rt;

    }

    return comp_data;
}


AMREX_INLINE AMREX_GPU_HOST_DEVICE
void get_xn(const int k, const init_t cd, amrex::Real *xn_zone, bool uniform_composition=false) {

    for (int n = 0; n < NumSpec; n++) {
        xn_zone[n] = 0.0_rt;
    }

    if (uniform_composition) {

        for (int n = 0; n < NumSpec; ++n) {
            xn_zone[n] = 1.0_rt / NumSpec;
        }

    }
    else {

        if (k < cd.n1) {
            if (cd.nprim >= 2) xn_zone[cd.is2] = cd.Xp_min/2;
            if (cd.nprim >= 3) xn_zone[cd.is3] = cd.Xp_min/2;

            amrex::Real dX = (cd.Xp_max - cd.Xp_min) /
                amrex::max((cd.n1 - 1), 1);
            xn_zone[cd.is1] = cd.Xp_min + k * dX;

        } else if (cd.nprim >= 2 && k < cd.n1 + cd.n2) {
            xn_zone[cd.is1] = cd.Xp_min/2;
            if (cd.nprim >= 3) {
                xn_zone[cd.is3] = cd.Xp_min/2;
            }
            amrex::Real dX = (cd.Xp_max - cd.Xp_min) /
                amrex::max((cd.n2 - 1), 1);
            xn_zone[cd.is2] = cd.Xp_min + (k - cd.n1) * dX;

        } else {
            xn_zone[cd.is1] = cd.Xp_min/2;
            xn_zone[cd.is2] = cd.Xp_min/2;

            amrex::Real dX = (cd.Xp_max - cd.Xp_min) /
                amrex::max((cd.n3 - 1), 1);
            xn_zone[cd.is3] = cd.Xp_min + (k - (cd.n1 + cd.n2)) * dX;

        }

        amrex::Real excess = 0.0_rt;
        for (int n = 0; n < NumSpec; n++) {
            excess += xn_zone[n];
        }
        excess = 1.0_rt - excess;

        for (int n = 0; n < NumSpec; n++) {
            if (n == cd.is1 ||
                (cd.nprim >= 2 && n == cd.is2) ||
                (cd.nprim >= 3 && n == cd.is3)) {
                continue;
            }

            xn_zone[n] = excess / (NumSpec - cd.nprim);
        }

    }


    // normalize -- just in case

    amrex::Real sum_X = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        sum_X += xn_zone[n];
    }
    for (int n = 0; n < NumSpec; n++) {
        xn_zone[n] = xn_zone[n] / sum_X;
    }

}

///
/// return the value of the runtime parameter XN, where N is
/// the integer index passed in
///

AMREX_INLINE AMREX_GPU_HOST_DEVICE
amrex::Real get_xn(const int index, bool uniform_composition=false) {

    amrex::Real mass_fraction{};

    if (uniform_composition) {
        mass_fraction = 1.0_rt / NumSpec;
        return mass_fraction;
    }

    switch (index) {

    case 1:
        mass_fraction = unit_test_rp::X1;
        break;
    case 2:
        mass_fraction = unit_test_rp::X2;
        break;
    case 3:
        mass_fraction = unit_test_rp::X3;
        break;
    case 4:
        mass_fraction = unit_test_rp::X4;
        break;
    case 5:
        mass_fraction = unit_test_rp::X5;
        break;
    case 6:
        mass_fraction = unit_test_rp::X6;
        break;
    case 7:
        mass_fraction = unit_test_rp::X7;
        break;
    case 8:
        mass_fraction = unit_test_rp::X8;
        break;
    case 9:
        mass_fraction = unit_test_rp::X9;
        break;
    case 10:
        mass_fraction = unit_test_rp::X10;
        break;
    case 11:
        mass_fraction = unit_test_rp::X11;
        break;
    case 12:
        mass_fraction = unit_test_rp::X12;
        break;
    case 13:
        mass_fraction = unit_test_rp::X13;
        break;
    case 14:
        mass_fraction = unit_test_rp::X14;
        break;
    case 15:
        mass_fraction = unit_test_rp::X15;
        break;
    case 16:
        mass_fraction = unit_test_rp::X16;
        break;
    case 17:
        mass_fraction = unit_test_rp::X17;
        break;
    case 18:
        mass_fraction = unit_test_rp::X18;
        break;
    case 19:
        mass_fraction = unit_test_rp::X19;
        break;
    case 20:
        mass_fraction = unit_test_rp::X20;
        break;
    case 21:
        mass_fraction = unit_test_rp::X21;
        break;
    case 22:
        mass_fraction = unit_test_rp::X22;
        break;
    case 23:
        mass_fraction = unit_test_rp::X23;
        break;
    case 24:
        mass_fraction = unit_test_rp::X24;
        break;
    case 25:
        mass_fraction = unit_test_rp::X25;
        break;
    case 26:
        mass_fraction = unit_test_rp::X26;
        break;
    case 27:
        mass_fraction = unit_test_rp::X27;
        break;
    case 28:
        mass_fraction = unit_test_rp::X28;
        break;
    case 29:
        mass_fraction = unit_test_rp::X29;
        break;
    case 30:
        mass_fraction = unit_test_rp::X30;
        break;
    case 31:
        mass_fraction = unit_test_rp::X31;
        break;
    case 32:
        mass_fraction = unit_test_rp::X32;
        break;
    case 33:
        mass_fraction = unit_test_rp::X33;
        break;
    case 34:
        mass_fraction = unit_test_rp::X34;
        break;
    case 35:
        mass_fraction = unit_test_rp::X35;
        break;

    }

    return mass_fraction;

}

#endif
