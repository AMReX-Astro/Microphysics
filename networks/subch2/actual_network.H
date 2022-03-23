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
        k_n__p__weak__wc12 = 1,
        k_c14__n14__weak__wc12 = 2,
        k_f18__o18__weak__wc12 = 3,
        k_mg23__na23__weak__wc12 = 4,
        k_si27__al27__weak__wc12 = 5,
        k_s31__p31__weak__wc12 = 6,
        k_n13__p_c12 = 7,
        k_n14__n_n13 = 8,
        k_o16__he4_c12 = 9,
        k_o18__he4_c14 = 10,
        k_f18__he4_n14 = 11,
        k_ne20__he4_o16 = 12,
        k_ne21__n_ne20 = 13,
        k_mg24__n_mg23 = 14,
        k_mg24__p_na23 = 15,
        k_mg24__he4_ne20 = 16,
        k_al27__he4_na23 = 17,
        k_si27__he4_mg23 = 18,
        k_si28__n_si27 = 19,
        k_si28__p_al27 = 20,
        k_si28__he4_mg24 = 21,
        k_p31__he4_al27 = 22,
        k_s31__he4_si27 = 23,
        k_s32__n_s31 = 24,
        k_s32__p_p31 = 25,
        k_s32__he4_si28 = 26,
        k_cl35__he4_p31 = 27,
        k_ar36__p_cl35 = 28,
        k_ar36__he4_s32 = 29,
        k_k39__he4_cl35 = 30,
        k_ca40__p_k39 = 31,
        k_ca40__he4_ar36 = 32,
        k_sc43__he4_k39 = 33,
        k_ti44__p_sc43 = 34,
        k_ti44__he4_ca40 = 35,
        k_v47__he4_sc43 = 36,
        k_cr48__p_v47 = 37,
        k_cr48__he4_ti44 = 38,
        k_mn51__he4_v47 = 39,
        k_fe52__p_mn51 = 40,
        k_fe52__he4_cr48 = 41,
        k_co55__he4_mn51 = 42,
        k_ni56__p_co55 = 43,
        k_ni56__he4_fe52 = 44,
        k_c12__he4_he4_he4 = 45,
        k_p_c12__n13 = 46,
        k_he4_c12__o16 = 47,
        k_he4_c14__o18 = 48,
        k_n_n13__n14 = 49,
        k_he4_n14__f18 = 50,
        k_he4_o16__ne20 = 51,
        k_n_ne20__ne21 = 52,
        k_he4_ne20__mg24 = 53,
        k_p_na23__mg24 = 54,
        k_he4_na23__al27 = 55,
        k_n_mg23__mg24 = 56,
        k_he4_mg23__si27 = 57,
        k_he4_mg24__si28 = 58,
        k_p_al27__si28 = 59,
        k_he4_al27__p31 = 60,
        k_n_si27__si28 = 61,
        k_he4_si27__s31 = 62,
        k_he4_si28__s32 = 63,
        k_p_p31__s32 = 64,
        k_he4_p31__cl35 = 65,
        k_n_s31__s32 = 66,
        k_he4_s32__ar36 = 67,
        k_p_cl35__ar36 = 68,
        k_he4_cl35__k39 = 69,
        k_he4_ar36__ca40 = 70,
        k_p_k39__ca40 = 71,
        k_he4_k39__sc43 = 72,
        k_he4_ca40__ti44 = 73,
        k_p_sc43__ti44 = 74,
        k_he4_sc43__v47 = 75,
        k_he4_ti44__cr48 = 76,
        k_p_v47__cr48 = 77,
        k_he4_v47__mn51 = 78,
        k_he4_cr48__fe52 = 79,
        k_p_mn51__fe52 = 80,
        k_he4_mn51__co55 = 81,
        k_he4_fe52__ni56 = 82,
        k_p_co55__ni56 = 83,
        k_c12_c12__n_mg23 = 84,
        k_c12_c12__p_na23 = 85,
        k_c12_c12__he4_ne20 = 86,
        k_p_c14__n_n14 = 87,
        k_he4_n13__p_o16 = 88,
        k_n_n14__p_c14 = 89,
        k_p_o16__he4_n13 = 90,
        k_c12_o16__n_si27 = 91,
        k_c12_o16__p_al27 = 92,
        k_c12_o16__he4_mg24 = 93,
        k_o16_o16__n_s31 = 94,
        k_o16_o16__p_p31 = 95,
        k_o16_o16__he4_si28 = 96,
        k_p_o18__n_f18 = 97,
        k_he4_o18__n_ne21 = 98,
        k_n_f18__p_o18 = 99,
        k_he4_f18__p_ne21 = 100,
        k_he4_ne20__n_mg23 = 101,
        k_he4_ne20__p_na23 = 102,
        k_he4_ne20__c12_c12 = 103,
        k_c12_ne20__n_s31 = 104,
        k_c12_ne20__p_p31 = 105,
        k_c12_ne20__he4_si28 = 106,
        k_n_ne21__he4_o18 = 107,
        k_p_ne21__he4_f18 = 108,
        k_he4_ne21__n_mg24 = 109,
        k_p_na23__n_mg23 = 110,
        k_p_na23__he4_ne20 = 111,
        k_p_na23__c12_c12 = 112,
        k_n_mg23__p_na23 = 113,
        k_n_mg23__he4_ne20 = 114,
        k_n_mg23__c12_c12 = 115,
        k_n_mg24__he4_ne21 = 116,
        k_he4_mg24__n_si27 = 117,
        k_he4_mg24__p_al27 = 118,
        k_he4_mg24__c12_o16 = 119,
        k_p_al27__n_si27 = 120,
        k_p_al27__he4_mg24 = 121,
        k_p_al27__c12_o16 = 122,
        k_n_si27__p_al27 = 123,
        k_n_si27__he4_mg24 = 124,
        k_n_si27__c12_o16 = 125,
        k_he4_si28__n_s31 = 126,
        k_he4_si28__p_p31 = 127,
        k_he4_si28__c12_ne20 = 128,
        k_he4_si28__o16_o16 = 129,
        k_p_p31__n_s31 = 130,
        k_p_p31__he4_si28 = 131,
        k_p_p31__c12_ne20 = 132,
        k_p_p31__o16_o16 = 133,
        k_n_s31__p_p31 = 134,
        k_n_s31__he4_si28 = 135,
        k_n_s31__c12_ne20 = 136,
        k_n_s31__o16_o16 = 137,
        k_he4_s32__p_cl35 = 138,
        k_p_cl35__he4_s32 = 139,
        k_he4_ar36__p_k39 = 140,
        k_p_k39__he4_ar36 = 141,
        k_he4_ca40__p_sc43 = 142,
        k_p_sc43__he4_ca40 = 143,
        k_he4_ti44__p_v47 = 144,
        k_p_v47__he4_ti44 = 145,
        k_he4_cr48__p_mn51 = 146,
        k_p_mn51__he4_cr48 = 147,
        k_he4_fe52__p_co55 = 148,
        k_p_co55__he4_fe52 = 149,
        k_he4_he4_he4__c12 = 150,
        NumRates = k_he4_he4_he4__c12
    };

    // For each rate, we need: rate, drate/dT, screening, dscreening/dT

    const int NumGroups = 4;

    // number of reaclib rates

    const int NrateReaclib = 150;
    const int NumReaclibSets = 232;

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
