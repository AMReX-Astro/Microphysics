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
        k_b8__he4_he4__weak__wc12 = 5,
        k_p_p__d__weak__bet_pos_ = 6,
        k_p_p__d__weak__electron_capture = 7,
        k_p_d__he3 = 8,
        k_d_d__he4 = 9,
        k_p_he3__he4__weak__bet_pos_ = 10,
        k_he4_he3__be7 = 11,
        k_p_be7__b8 = 12,
        k_p_c12__n13 = 13,
        k_he4_c12__o16 = 14,
        k_p_c13__n14 = 15,
        k_p_n13__o14 = 16,
        k_p_n14__o15 = 17,
        k_he4_n14__f18 = 18,
        k_p_n15__o16 = 19,
        k_p_o16__f17 = 20,
        k_p_o17__f18 = 21,
        k_d_he3__p_he4 = 22,
        k_he4_n13__p_o16 = 23,
        k_p_n15__he4_c12 = 24,
        k_he4_o14__p_f17 = 25,
        k_p_o17__he4_n14 = 26,
        k_p_f18__he4_o15 = 27,
        k_he3_he3__p_p_he4 = 28,
        k_d_be7__p_he4_he4 = 29,
        k_he3_be7__p_p_he4_he4 = 30,
        k_he4_he4_he4__c12 = 31,
        NumRates = k_he4_he4_he4__c12
    };

    // number of reaclib rates

    const int NrateReaclib = 31;
    const int NumReaclibSets = 64;

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
