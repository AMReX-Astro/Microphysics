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
        k_mg24__he4_ne20 = 8,
        k_si28__p_al27 = 9,
        k_si28__he4_mg24 = 10,
        k_p31__he4_al27 = 11,
        k_s32__p_p31 = 12,
        k_s32__he4_si28 = 13,
        k_cl35__he4_p31 = 14,
        k_ar36__p_cl35 = 15,
        k_ar36__he4_s32 = 16,
        k_k39__he4_cl35 = 17,
        k_ca40__p_k39 = 18,
        k_ca40__he4_ar36 = 19,
        k_sc43__he4_k39 = 20,
        k_ti44__p_sc43 = 21,
        k_ti44__he4_ca40 = 22,
        k_v47__he4_sc43 = 23,
        k_cr48__p_v47 = 24,
        k_cr48__he4_ti44 = 25,
        k_mn51__he4_v47 = 26,
        k_fe52__p_mn51 = 27,
        k_fe52__he4_cr48 = 28,
        k_co55__he4_mn51 = 29,
        k_ni56__p_co55 = 30,
        k_ni56__he4_fe52 = 31,
        k_c12__he4_he4_he4 = 32,
        k_p_c12__n13 = 33,
        k_he4_c12__o16 = 34,
        k_he4_c14__o18 = 35,
        k_he4_n14__f18 = 36,
        k_he4_o16__ne20 = 37,
        k_he4_ne20__mg24 = 38,
        k_he4_mg24__si28 = 39,
        k_p_al27__si28 = 40,
        k_he4_al27__p31 = 41,
        k_he4_si28__s32 = 42,
        k_p_p31__s32 = 43,
        k_he4_p31__cl35 = 44,
        k_he4_s32__ar36 = 45,
        k_p_cl35__ar36 = 46,
        k_he4_cl35__k39 = 47,
        k_he4_ar36__ca40 = 48,
        k_p_k39__ca40 = 49,
        k_he4_k39__sc43 = 50,
        k_he4_ca40__ti44 = 51,
        k_p_sc43__ti44 = 52,
        k_he4_sc43__v47 = 53,
        k_he4_ti44__cr48 = 54,
        k_p_v47__cr48 = 55,
        k_he4_v47__mn51 = 56,
        k_he4_cr48__fe52 = 57,
        k_p_mn51__fe52 = 58,
        k_he4_mn51__co55 = 59,
        k_he4_fe52__ni56 = 60,
        k_p_co55__ni56 = 61,
        k_c12_c12__he4_ne20 = 62,
        k_he4_n13__p_o16 = 63,
        k_p_o16__he4_n13 = 64,
        k_c12_o16__p_al27 = 65,
        k_c12_o16__he4_mg24 = 66,
        k_o16_o16__p_p31 = 67,
        k_o16_o16__he4_si28 = 68,
        k_he4_f18__p_ne21 = 69,
        k_he4_ne20__c12_c12 = 70,
        k_c12_ne20__p_p31 = 71,
        k_c12_ne20__he4_si28 = 72,
        k_p_ne21__he4_f18 = 73,
        k_he4_mg24__p_al27 = 74,
        k_he4_mg24__c12_o16 = 75,
        k_p_al27__he4_mg24 = 76,
        k_p_al27__c12_o16 = 77,
        k_he4_si28__p_p31 = 78,
        k_he4_si28__c12_ne20 = 79,
        k_he4_si28__o16_o16 = 80,
        k_p_p31__he4_si28 = 81,
        k_p_p31__c12_ne20 = 82,
        k_p_p31__o16_o16 = 83,
        k_he4_s32__p_cl35 = 84,
        k_p_cl35__he4_s32 = 85,
        k_he4_ar36__p_k39 = 86,
        k_p_k39__he4_ar36 = 87,
        k_he4_ca40__p_sc43 = 88,
        k_p_sc43__he4_ca40 = 89,
        k_he4_ti44__p_v47 = 90,
        k_p_v47__he4_ti44 = 91,
        k_he4_cr48__p_mn51 = 92,
        k_p_mn51__he4_cr48 = 93,
        k_he4_fe52__p_co55 = 94,
        k_p_co55__he4_fe52 = 95,
        k_he4_he4_he4__c12 = 96,
        NumRates = k_he4_he4_he4__c12
    };

    // For each rate, we need: rate, drate/dT, screening, dscreening/dT

    const int NumGroups = 4;

    // number of reaclib rates

    const int NrateReaclib = 96;
    const int NumReaclibSets = 158;

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
