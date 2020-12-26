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

    comp_data.is1 = network_spec_index(primary_species_1);
    if (comp_data.is1 >= 0) {
        comp_data.nprim++;
    }

    comp_data.is2 = network_spec_index(primary_species_2);
    if (comp_data.is2 >= 0) {
        comp_data.nprim++;
    }

    comp_data.is3 = network_spec_index(primary_species_3);
    if (comp_data.is3 >= 0) {
        comp_data.nprim++;
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
void get_xn(const int k, const init_t cd, Real *xn_zone) {

    for (int n = 0; n < NumSpec; n++) {
        xn_zone[n] = 0.0_rt;
    }

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


    // normalize -- just in case

    Real sum_X = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        sum_X += xn_zone[n];
    }
    for (int n = 0; n < NumSpec; n++) {
        xn_zone[n] = xn_zone[n] / sum_X;
    }

}
#endif

