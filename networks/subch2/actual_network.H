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
        k_na22__ne22__weak__wc12 = 3,
        k_n13__p_c12 = 4,
        k_o16__he4_c12 = 5,
        k_o18__he4_c14 = 6,
        k_f18__he4_n14 = 7,
        k_ne20__he4_o16 = 8,
        k_ne22__he4_o18 = 9,
        k_na22__p_ne21 = 10,
        k_na22__he4_f18 = 11,
        k_na23__p_ne22 = 12,
        k_mg24__p_na23 = 13,
        k_mg24__he4_ne20 = 14,
        k_al27__he4_na23 = 15,
        k_si28__p_al27 = 16,
        k_si28__he4_mg24 = 17,
        k_p31__he4_al27 = 18,
        k_s32__p_p31 = 19,
        k_s32__he4_si28 = 20,
        k_cl35__he4_p31 = 21,
        k_ar36__p_cl35 = 22,
        k_ar36__he4_s32 = 23,
        k_k39__he4_cl35 = 24,
        k_ca40__p_k39 = 25,
        k_ca40__he4_ar36 = 26,
        k_sc43__he4_k39 = 27,
        k_ti44__p_sc43 = 28,
        k_ti44__he4_ca40 = 29,
        k_v47__he4_sc43 = 30,
        k_cr48__p_v47 = 31,
        k_cr48__he4_ti44 = 32,
        k_mn51__he4_v47 = 33,
        k_fe52__p_mn51 = 34,
        k_fe52__he4_cr48 = 35,
        k_co55__he4_mn51 = 36,
        k_ni56__p_co55 = 37,
        k_ni56__he4_fe52 = 38,
        k_c12__he4_he4_he4 = 39,
        k_p_c12__n13 = 40,
        k_he4_c12__o16 = 41,
        k_he4_c14__o18 = 42,
        k_he4_n14__f18 = 43,
        k_he4_o16__ne20 = 44,
        k_he4_o18__ne22 = 45,
        k_he4_f18__na22 = 46,
        k_he4_ne20__mg24 = 47,
        k_p_ne21__na22 = 48,
        k_p_ne22__na23 = 49,
        k_p_na23__mg24 = 50,
        k_he4_na23__al27 = 51,
        k_he4_mg24__si28 = 52,
        k_p_al27__si28 = 53,
        k_he4_al27__p31 = 54,
        k_he4_si28__s32 = 55,
        k_p_p31__s32 = 56,
        k_he4_p31__cl35 = 57,
        k_he4_s32__ar36 = 58,
        k_p_cl35__ar36 = 59,
        k_he4_cl35__k39 = 60,
        k_he4_ar36__ca40 = 61,
        k_p_k39__ca40 = 62,
        k_he4_k39__sc43 = 63,
        k_he4_ca40__ti44 = 64,
        k_p_sc43__ti44 = 65,
        k_he4_sc43__v47 = 66,
        k_he4_ti44__cr48 = 67,
        k_p_v47__cr48 = 68,
        k_he4_v47__mn51 = 69,
        k_he4_cr48__fe52 = 70,
        k_p_mn51__fe52 = 71,
        k_he4_mn51__co55 = 72,
        k_he4_fe52__ni56 = 73,
        k_p_co55__ni56 = 74,
        k_c12_c12__p_na23 = 75,
        k_c12_c12__he4_ne20 = 76,
        k_he4_n13__p_o16 = 77,
        k_p_o16__he4_n13 = 78,
        k_c12_o16__p_al27 = 79,
        k_c12_o16__he4_mg24 = 80,
        k_o16_o16__p_p31 = 81,
        k_o16_o16__he4_si28 = 82,
        k_he4_f18__p_ne21 = 83,
        k_he4_ne20__p_na23 = 84,
        k_he4_ne20__c12_c12 = 85,
        k_c12_ne20__p_p31 = 86,
        k_c12_ne20__he4_si28 = 87,
        k_p_ne21__he4_f18 = 88,
        k_p_na23__he4_ne20 = 89,
        k_p_na23__c12_c12 = 90,
        k_he4_mg24__p_al27 = 91,
        k_he4_mg24__c12_o16 = 92,
        k_p_al27__he4_mg24 = 93,
        k_p_al27__c12_o16 = 94,
        k_he4_si28__p_p31 = 95,
        k_he4_si28__c12_ne20 = 96,
        k_he4_si28__o16_o16 = 97,
        k_p_p31__he4_si28 = 98,
        k_p_p31__c12_ne20 = 99,
        k_p_p31__o16_o16 = 100,
        k_he4_s32__p_cl35 = 101,
        k_p_cl35__he4_s32 = 102,
        k_he4_ar36__p_k39 = 103,
        k_p_k39__he4_ar36 = 104,
        k_he4_ca40__p_sc43 = 105,
        k_p_sc43__he4_ca40 = 106,
        k_he4_ti44__p_v47 = 107,
        k_p_v47__he4_ti44 = 108,
        k_he4_cr48__p_mn51 = 109,
        k_p_mn51__he4_cr48 = 110,
        k_he4_fe52__p_co55 = 111,
        k_p_co55__he4_fe52 = 112,
        k_he4_he4_he4__c12 = 113,
        k_c12_c12__n_mg23 = 114,
        k_n_mg23__mg24 = 115,
        k_o16_o16__n_s31 = 116,
        k_n_s31__s32 = 117,
        k_c12_o16__n_si27 = 118,
        k_n_si27__si28 = 119,
        NumRates = k_n_si27__si28
    };

    // For each rate, we need: rate, drate/dT, screening, dscreening/dT

    const int NumGroups = 4;

    // number of reaclib rates

    const int NrateReaclib = 119;
    const int NumReaclibSets = 207;

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
