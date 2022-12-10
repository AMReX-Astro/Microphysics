#ifndef REACT_UTIL_H
#define REACT_UTIL_H

using namespace amrex;

struct init_t {
    int nprim;
    int is1;
    int is2;
    int is3;

    int n1;
    int n2;
    int n3;

    Real Xp_min;
    Real Xp_max;
};

AMREX_INLINE
init_t setup_composition(const int nz) {

    // get the primary species indices
    init_t comp_data;

    comp_data.nprim = 0;

    // we absolutely require primary_species_1 to be defined
    comp_data.is1 = network_spec_index(primary_species_1);
    if (comp_data.is1 >= 0) {
        comp_data.nprim++;
    } else {
        amrex::Error("Error: invalid primary_species_1");
    }

    // we'll check the next species, but if it is not valid,
    // we'll just move on
    comp_data.is2 = network_spec_index(primary_species_2);
    if (comp_data.is2 >= 0) {
        comp_data.nprim++;

        // only consider primary_species_3 if primary_species_2
        // was defined
        comp_data.is3 = network_spec_index(primary_species_3);
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
void get_xn(const int k, const init_t cd, Real *xn_zone, int uniform_composition = 0) {

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

            Real dX = (cd.Xp_max - cd.Xp_min) /
                amrex::max((cd.n1 - 1), 1);
            xn_zone[cd.is1] = cd.Xp_min + k * dX;

        } else if (cd.nprim >= 2 && k < cd.n1 + cd.n2) {
            xn_zone[cd.is1] = cd.Xp_min/2;
            if (cd.nprim >= 3) {
                xn_zone[cd.is3] = cd.Xp_min/2;
            }
            Real dX = (cd.Xp_max - cd.Xp_min) /
                amrex::max((cd.n2 - 1), 1);
            xn_zone[cd.is2] = cd.Xp_min + (k - cd.n1) * dX;

        } else {
            xn_zone[cd.is1] = cd.Xp_min/2;
            xn_zone[cd.is2] = cd.Xp_min/2;

            Real dX = (cd.Xp_max - cd.Xp_min) /
                amrex::max((cd.n3 - 1), 1);
            xn_zone[cd.is3] = cd.Xp_min + (k - (cd.n1 + cd.n2)) * dX;

        }

        Real excess = 0.0_rt;
        for (int n = 0; n < NumSpec; n++) {
            excess += xn_zone[n];
        }
        excess = 1.0_rt - excess;

        for (int n = 0; n < NumSpec; n++) {
            if (n == cd.is1 ||
                (n == cd.is2 && cd.nprim >= 2) ||
                (n == cd.is3 && cd.nprim >= 3)) {
                continue;
            }

            xn_zone[n] = excess / (NumSpec - cd.nprim);
        }

    }


    // normalize -- just in case

    Real sum_X = 0.0_rt;
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
Real get_xn(const int index) {

    Real mass_fraction{};

    switch (index) {

    case 1:
        mass_fraction = X1;
        break;
    case 2:
        mass_fraction = X2;
        break;
    case 3:
        mass_fraction = X3;
        break;
    case 4:
        mass_fraction = X4;
        break;
    case 5:
        mass_fraction = X5;
        break;
    case 6:
        mass_fraction = X6;
        break;
    case 7:
        mass_fraction = X7;
        break;
    case 8:
        mass_fraction = X8;
        break;
    case 9:
        mass_fraction = X9;
        break;
    case 10:
        mass_fraction = X10;
        break;
    case 11:
        mass_fraction = X11;
        break;
    case 12:
        mass_fraction = X12;
        break;
    case 13:
        mass_fraction = X13;
        break;
    case 14:
        mass_fraction = X14;
        break;
    case 15:
        mass_fraction = X15;
        break;
    case 16:
        mass_fraction = X16;
        break;
    case 17:
        mass_fraction = X17;
        break;
    case 18:
        mass_fraction = X18;
        break;
    case 19:
        mass_fraction = X19;
        break;
    case 20:
        mass_fraction = X20;
        break;
    case 21:
        mass_fraction = X21;
        break;
    case 22:
        mass_fraction = X22;
        break;
    case 23:
        mass_fraction = X23;
        break;
    case 24:
        mass_fraction = X24;
        break;
    case 25:
        mass_fraction = X25;
        break;
    case 26:
        mass_fraction = X26;
        break;
    case 27:
        mass_fraction = X27;
        break;
    case 28:
        mass_fraction = X28;
        break;
    case 29:
        mass_fraction = X29;
        break;
    case 30:
        mass_fraction = X30;
        break;
    case 31:
        mass_fraction = X31;
        break;
    case 32:
        mass_fraction = X32;
        break;
    case 33:
        mass_fraction = X33;
        break;
    case 34:
        mass_fraction = X34;
        break;
    case 35:
        mass_fraction = X35;
        break;

    }

    return mass_fraction;

}

#endif

