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
        k_n13__c13__weak__wc12 = 1,
        k_o14__n14__weak__wc12 = 2,
        k_o15__n15__weak__wc12 = 3,
        k_f17__o17__weak__wc12 = 4,
        k_p_c12__n13 = 5,
        k_he4_c12__o16 = 6,
        k_p_c13__n14 = 7,
        k_p_n13__o14 = 8,
        k_p_n14__o15 = 9,
        k_he4_n14__f18 = 10,
        k_p_n15__o16 = 11,
        k_p_o16__f17 = 12,
        k_p_o17__f18 = 13,
        k_he4_n13__p_o16 = 14,
        k_p_n15__he4_c12 = 15,
        k_he4_o14__p_f17 = 16,
        k_p_o17__he4_n14 = 17,
        k_p_f18__he4_o15 = 18,
        k_he4_he4_he4__c12 = 19,
        NumRates = k_he4_he4_he4__c12
    };

    // For each rate, we need: rate, drate/dT, screening, dscreening/dT

    const int NumGroups = 4;

    // number of reaclib rates

    const int NrateReaclib = 19;
    const int NumReaclibSets = 48;

    // number of tabular rates

    const int NrateTabular = 0;

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
