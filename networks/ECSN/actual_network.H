#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex;

void actual_network_init();

const std::string network_name = "pynucastro-cxx";

namespace network
{
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
}

namespace Rates
{

    enum NetworkRates
    {
        k_ne20__he4_o16 = 1,
        k_he4_o16__ne20 = 2,
        k_he4_ne20__mg24 = 3,
        k_he4_mg24__si28 = 4,
        k_p_al27__si28 = 5,
        k_he4_al27__p31 = 6,
        k_he4_si28__s32 = 7,
        k_p_p31__s32 = 8,
        k_o16_o16__p_p31 = 9,
        k_o16_o16__he4_si28 = 10,
        k_he4_mg24__p_al27 = 11,
        k_p_al27__he4_mg24 = 12,
        k_he4_si28__p_p31 = 13,
        k_p_p31__he4_si28 = 14,
        k_f20__o20 = 15,
        k_ne20__f20 = 16,
        k_o20__f20 = 17,
        k_f20__ne20 = 18,
        NumRates = k_f20__ne20
    };

    // For each rate, we need: rate, drate/dT, screening, dscreening/dT

    const int NumGroups = 4;

    // number of reaclib rates

    const int NrateReaclib = 14;
    const int NumReaclibSets = 34;

    // number of tabular rates

    const int NrateTabular = 4;

}

namespace reaclib_rates
{

    // Temperature coefficient arrays (numbers correspond to reaction
    // numbers in net_info)

    extern AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 1, 7, 1, Rates::NumReaclibSets> ctemp_rate;

    // Index into ctemp_rate, dimension 2, where each rate's
    // coefficients start

    extern AMREX_GPU_MANAGED amrex::Array1D<int, 1, Rates::NrateReaclib> rate_start_idx;

    // Reaction multiplicities-1 (how many rates contribute - 1)

    extern AMREX_GPU_MANAGED amrex::Array1D<int, 1, Rates::NrateReaclib> rate_extra_mult;

}

#endif
