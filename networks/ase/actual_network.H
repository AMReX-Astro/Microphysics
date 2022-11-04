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
        k_n13__p_c12 = 2,
        k_n14__n_n13 = 3,
        k_o16__he4_c12 = 4,
        k_f18__he4_n14 = 5,
        k_ne20__he4_o16 = 6,
        k_ne21__n_ne20 = 7,
        k_na22__p_ne21 = 8,
        k_na22__he4_f18 = 9,
        k_na23__n_na22 = 10,
        k_mg24__p_na23 = 11,
        k_mg24__he4_ne20 = 12,
        k_al27__he4_na23 = 13,
        k_si28__p_al27 = 14,
        k_si28__he4_mg24 = 15,
        k_p31__he4_al27 = 16,
        k_s32__p_p31 = 17,
        k_s32__he4_si28 = 18,
        k_ar36__he4_s32 = 19,
        k_ca40__he4_ar36 = 20,
        k_ti44__he4_ca40 = 21,
        k_cr48__he4_ti44 = 22,
        k_fe52__he4_cr48 = 23,
        k_ni56__he4_fe52 = 24,
        k_c12__he4_he4_he4 = 25,
        k_p_c12__n13 = 26,
        k_he4_c12__o16 = 27,
        k_n_n13__n14 = 28,
        k_he4_n14__f18 = 29,
        k_he4_o16__ne20 = 30,
        k_he4_f18__na22 = 31,
        k_n_ne20__ne21 = 32,
        k_he4_ne20__mg24 = 33,
        k_p_ne21__na22 = 34,
        k_n_na22__na23 = 35,
        k_p_na23__mg24 = 36,
        k_he4_na23__al27 = 37,
        k_he4_mg24__si28 = 38,
        k_p_al27__si28 = 39,
        k_he4_al27__p31 = 40,
        k_he4_si28__s32 = 41,
        k_p_p31__s32 = 42,
        k_he4_s32__ar36 = 43,
        k_he4_ar36__ca40 = 44,
        k_he4_ca40__ti44 = 45,
        k_he4_ti44__cr48 = 46,
        k_he4_cr48__fe52 = 47,
        k_he4_fe52__ni56 = 48,
        k_c12_c12__p_na23 = 49,
        k_c12_c12__he4_ne20 = 50,
        k_he4_n13__p_o16 = 51,
        k_p_o16__he4_n13 = 52,
        k_c12_o16__p_al27 = 53,
        k_c12_o16__he4_mg24 = 54,
        k_o16_o16__p_p31 = 55,
        k_o16_o16__he4_si28 = 56,
        k_he4_f18__p_ne21 = 57,
        k_he4_ne20__p_na23 = 58,
        k_he4_ne20__c12_c12 = 59,
        k_c12_ne20__p_p31 = 60,
        k_c12_ne20__he4_si28 = 61,
        k_p_ne21__he4_f18 = 62,
        k_he4_ne21__n_mg24 = 63,
        k_p_na23__he4_ne20 = 64,
        k_p_na23__c12_c12 = 65,
        k_n_mg24__he4_ne21 = 66,
        k_he4_mg24__p_al27 = 67,
        k_he4_mg24__c12_o16 = 68,
        k_p_al27__he4_mg24 = 69,
        k_p_al27__c12_o16 = 70,
        k_he4_si28__p_p31 = 71,
        k_he4_si28__c12_ne20 = 72,
        k_he4_si28__o16_o16 = 73,
        k_p_p31__he4_si28 = 74,
        k_p_p31__c12_ne20 = 75,
        k_p_p31__o16_o16 = 76,
        k_he4_he4_he4__c12 = 77,
        NumRates = k_he4_he4_he4__c12
    };

    // number of reaclib rates

    const int NrateReaclib = 77;

    // number of tabular rates

    const int NrateTabular = 0;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "n__p__weak__wc12",  // 1,
        "n13__p_c12",  // 2,
        "n14__n_n13",  // 3,
        "o16__he4_c12",  // 4,
        "f18__he4_n14",  // 5,
        "ne20__he4_o16",  // 6,
        "ne21__n_ne20",  // 7,
        "na22__p_ne21",  // 8,
        "na22__he4_f18",  // 9,
        "na23__n_na22",  // 10,
        "mg24__p_na23",  // 11,
        "mg24__he4_ne20",  // 12,
        "al27__he4_na23",  // 13,
        "si28__p_al27",  // 14,
        "si28__he4_mg24",  // 15,
        "p31__he4_al27",  // 16,
        "s32__p_p31",  // 17,
        "s32__he4_si28",  // 18,
        "ar36__he4_s32",  // 19,
        "ca40__he4_ar36",  // 20,
        "ti44__he4_ca40",  // 21,
        "cr48__he4_ti44",  // 22,
        "fe52__he4_cr48",  // 23,
        "ni56__he4_fe52",  // 24,
        "c12__he4_he4_he4",  // 25,
        "p_c12__n13",  // 26,
        "he4_c12__o16",  // 27,
        "n_n13__n14",  // 28,
        "he4_n14__f18",  // 29,
        "he4_o16__ne20",  // 30,
        "he4_f18__na22",  // 31,
        "n_ne20__ne21",  // 32,
        "he4_ne20__mg24",  // 33,
        "p_ne21__na22",  // 34,
        "n_na22__na23",  // 35,
        "p_na23__mg24",  // 36,
        "he4_na23__al27",  // 37,
        "he4_mg24__si28",  // 38,
        "p_al27__si28",  // 39,
        "he4_al27__p31",  // 40,
        "he4_si28__s32",  // 41,
        "p_p31__s32",  // 42,
        "he4_s32__ar36",  // 43,
        "he4_ar36__ca40",  // 44,
        "he4_ca40__ti44",  // 45,
        "he4_ti44__cr48",  // 46,
        "he4_cr48__fe52",  // 47,
        "he4_fe52__ni56",  // 48,
        "c12_c12__p_na23",  // 49,
        "c12_c12__he4_ne20",  // 50,
        "he4_n13__p_o16",  // 51,
        "p_o16__he4_n13",  // 52,
        "c12_o16__p_al27",  // 53,
        "c12_o16__he4_mg24",  // 54,
        "o16_o16__p_p31",  // 55,
        "o16_o16__he4_si28",  // 56,
        "he4_f18__p_ne21",  // 57,
        "he4_ne20__p_na23",  // 58,
        "he4_ne20__c12_c12",  // 59,
        "c12_ne20__p_p31",  // 60,
        "c12_ne20__he4_si28",  // 61,
        "p_ne21__he4_f18",  // 62,
        "he4_ne21__n_mg24",  // 63,
        "p_na23__he4_ne20",  // 64,
        "p_na23__c12_c12",  // 65,
        "n_mg24__he4_ne21",  // 66,
        "he4_mg24__p_al27",  // 67,
        "he4_mg24__c12_o16",  // 68,
        "p_al27__he4_mg24",  // 69,
        "p_al27__c12_o16",  // 70,
        "he4_si28__p_p31",  // 71,
        "he4_si28__c12_ne20",  // 72,
        "he4_si28__o16_o16",  // 73,
        "p_p31__he4_si28",  // 74,
        "p_p31__c12_ne20",  // 75,
        "p_p31__o16_o16",  // 76,
        "he4_he4_he4__c12"  // 77,
    };

}

#endif
