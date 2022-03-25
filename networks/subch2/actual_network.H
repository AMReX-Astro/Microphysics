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
        k_c14__n14__weak__wc12 = 1,
        k_f18__o18__weak__wc12 = 2,
        k_n13__p_c12 = 3,
        k_o16__he4_c12 = 4,
        k_o18__he4_c14 = 5,
        k_f18__he4_n14 = 6,
        k_ne20__he4_o16 = 7,
        k_mg24__p_na23 = 8,
        k_mg24__he4_ne20 = 9,
        k_al27__he4_na23 = 10,
        k_si28__p_al27 = 11,
        k_si28__he4_mg24 = 12,
        k_p31__he4_al27 = 13,
        k_s32__p_p31 = 14,
        k_s32__he4_si28 = 15,
        k_cl35__he4_p31 = 16,
        k_ar36__p_cl35 = 17,
        k_ar36__he4_s32 = 18,
        k_k39__he4_cl35 = 19,
        k_ca40__p_k39 = 20,
        k_ca40__he4_ar36 = 21,
        k_sc43__he4_k39 = 22,
        k_ti44__p_sc43 = 23,
        k_ti44__he4_ca40 = 24,
        k_v47__he4_sc43 = 25,
        k_cr48__p_v47 = 26,
        k_cr48__he4_ti44 = 27,
        k_mn51__he4_v47 = 28,
        k_fe52__p_mn51 = 29,
        k_fe52__he4_cr48 = 30,
        k_co55__he4_mn51 = 31,
        k_ni56__p_co55 = 32,
        k_ni56__he4_fe52 = 33,
        k_c12__he4_he4_he4 = 34,
        k_p_c12__n13 = 35,
        k_he4_c12__o16 = 36,
        k_he4_c14__o18 = 37,
        k_he4_n14__f18 = 38,
        k_he4_o16__ne20 = 39,
        k_he4_ne20__mg24 = 40,
        k_p_na23__mg24 = 41,
        k_he4_na23__al27 = 42,
        k_he4_mg24__si28 = 43,
        k_p_al27__si28 = 44,
        k_he4_al27__p31 = 45,
        k_he4_si28__s32 = 46,
        k_p_p31__s32 = 47,
        k_he4_p31__cl35 = 48,
        k_he4_s32__ar36 = 49,
        k_p_cl35__ar36 = 50,
        k_he4_cl35__k39 = 51,
        k_he4_ar36__ca40 = 52,
        k_p_k39__ca40 = 53,
        k_he4_k39__sc43 = 54,
        k_he4_ca40__ti44 = 55,
        k_p_sc43__ti44 = 56,
        k_he4_sc43__v47 = 57,
        k_he4_ti44__cr48 = 58,
        k_p_v47__cr48 = 59,
        k_he4_v47__mn51 = 60,
        k_he4_cr48__fe52 = 61,
        k_p_mn51__fe52 = 62,
        k_he4_mn51__co55 = 63,
        k_he4_fe52__ni56 = 64,
        k_p_co55__ni56 = 65,
        k_c12_c12__p_na23 = 66,
        k_c12_c12__he4_ne20 = 67,
        k_he4_n13__p_o16 = 68,
        k_p_o16__he4_n13 = 69,
        k_c12_o16__p_al27 = 70,
        k_c12_o16__he4_mg24 = 71,
        k_o16_o16__p_p31 = 72,
        k_o16_o16__he4_si28 = 73,
        k_he4_f18__p_ne21 = 74,
        k_he4_ne20__p_na23 = 75,
        k_he4_ne20__c12_c12 = 76,
        k_c12_ne20__p_p31 = 77,
        k_c12_ne20__he4_si28 = 78,
        k_p_ne21__he4_f18 = 79,
        k_p_na23__he4_ne20 = 80,
        k_p_na23__c12_c12 = 81,
        k_he4_mg24__p_al27 = 82,
        k_he4_mg24__c12_o16 = 83,
        k_p_al27__he4_mg24 = 84,
        k_p_al27__c12_o16 = 85,
        k_he4_si28__p_p31 = 86,
        k_he4_si28__c12_ne20 = 87,
        k_he4_si28__o16_o16 = 88,
        k_p_p31__he4_si28 = 89,
        k_p_p31__c12_ne20 = 90,
        k_p_p31__o16_o16 = 91,
        k_he4_s32__p_cl35 = 92,
        k_p_cl35__he4_s32 = 93,
        k_he4_ar36__p_k39 = 94,
        k_p_k39__he4_ar36 = 95,
        k_he4_ca40__p_sc43 = 96,
        k_p_sc43__he4_ca40 = 97,
        k_he4_ti44__p_v47 = 98,
        k_p_v47__he4_ti44 = 99,
        k_he4_cr48__p_mn51 = 100,
        k_p_mn51__he4_cr48 = 101,
        k_he4_fe52__p_co55 = 102,
        k_p_co55__he4_fe52 = 103,
        k_he4_he4_he4__c12 = 104,
        k_c12_c12__n_mg23 = 105,
        k_n_mg23__mg24 = 106,
        k_o16_o16__n_s31 = 107,
        k_n_s31__s32 = 108,
        k_c12_o16__n_si27 = 109,
        k_n_si27__si28 = 110,
        NumRates = k_n_si27__si28
    };

    // For each rate, we need: rate, drate/dT, screening, dscreening/dT

    const int NumGroups = 4;

    // number of reaclib rates

    const int NrateReaclib = 110;
    const int NumReaclibSets = 180;

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
