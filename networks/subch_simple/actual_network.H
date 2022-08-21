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
        k_n13__p_c12 = 1,
        k_o16__he4_c12 = 2,
        k_f18__he4_n14 = 3,
        k_ne20__he4_o16 = 4,
        k_na22__p_ne21 = 5,
        k_na22__he4_f18 = 6,
        k_mg24__p_na23 = 7,
        k_mg24__he4_ne20 = 8,
        k_al27__he4_na23 = 9,
        k_si28__p_al27 = 10,
        k_si28__he4_mg24 = 11,
        k_p31__he4_al27 = 12,
        k_s32__p_p31 = 13,
        k_s32__he4_si28 = 14,
        k_c12__he4_he4_he4 = 15,
        k_p_c12__n13 = 16,
        k_he4_c12__o16 = 17,
        k_he4_n14__f18 = 18,
        k_he4_o16__ne20 = 19,
        k_he4_f18__na22 = 20,
        k_he4_ne20__mg24 = 21,
        k_p_ne21__na22 = 22,
        k_p_na23__mg24 = 23,
        k_he4_na23__al27 = 24,
        k_he4_mg24__si28 = 25,
        k_p_al27__si28 = 26,
        k_he4_al27__p31 = 27,
        k_he4_si28__s32 = 28,
        k_p_p31__s32 = 29,
        k_c12_c12__p_na23 = 30,
        k_c12_c12__he4_ne20 = 31,
        k_he4_n13__p_o16 = 32,
        k_p_o16__he4_n13 = 33,
        k_c12_o16__p_al27 = 34,
        k_c12_o16__he4_mg24 = 35,
        k_o16_o16__p_p31 = 36,
        k_o16_o16__he4_si28 = 37,
        k_he4_f18__p_ne21 = 38,
        k_he4_ne20__p_na23 = 39,
        k_p_ne21__he4_f18 = 40,
        k_p_na23__he4_ne20 = 41,
        k_he4_mg24__p_al27 = 42,
        k_p_al27__he4_mg24 = 43,
        k_he4_si28__p_p31 = 44,
        k_p_p31__he4_si28 = 45,
        k_he4_he4_he4__c12 = 46,
        k_c12_c12__mg24__modified = 47,
        k_o16_o16__s32__modified = 48,
        k_c12_o16__si28__modified = 49,
        k_he4_s32__ar36 = 50,
        k_he4_s32__p_cl35 = 51,
        k_p_cl35__ar36 = 52,
        k_ar36__he4_s32 = 53,
        k_ar36__p_cl35 = 54,
        k_p_cl35__he4_s32 = 55,
        k_he4_ar36__ca40 = 56,
        k_he4_ar36__p_k39 = 57,
        k_p_k39__ca40 = 58,
        k_ca40__he4_ar36 = 59,
        k_ca40__p_k39 = 60,
        k_p_k39__he4_ar36 = 61,
        k_he4_ca40__ti44 = 62,
        k_he4_ca40__p_sc43 = 63,
        k_p_sc43__ti44 = 64,
        k_ti44__he4_ca40 = 65,
        k_ti44__p_sc43 = 66,
        k_p_sc43__he4_ca40 = 67,
        k_he4_ti44__cr48 = 68,
        k_he4_ti44__p_v47 = 69,
        k_p_v47__cr48 = 70,
        k_cr48__he4_ti44 = 71,
        k_cr48__p_v47 = 72,
        k_p_v47__he4_ti44 = 73,
        k_he4_cr48__fe52 = 74,
        k_he4_cr48__p_mn51 = 75,
        k_p_mn51__fe52 = 76,
        k_fe52__he4_cr48 = 77,
        k_fe52__p_mn51 = 78,
        k_p_mn51__he4_cr48 = 79,
        k_he4_fe52__ni56 = 80,
        k_he4_fe52__p_co55 = 81,
        k_p_co55__ni56 = 82,
        k_ni56__he4_fe52 = 83,
        k_ni56__p_co55 = 84,
        k_p_co55__he4_fe52 = 85,
        k_s32_he4__ar36__approx = 86,
        k_ar36__s32_he4__approx = 87,
        k_ar36_he4__ca40__approx = 88,
        k_ca40__ar36_he4__approx = 89,
        k_ca40_he4__ti44__approx = 90,
        k_ti44__ca40_he4__approx = 91,
        k_ti44_he4__cr48__approx = 92,
        k_cr48__ti44_he4__approx = 93,
        k_cr48_he4__fe52__approx = 94,
        k_fe52__cr48_he4__approx = 95,
        k_fe52_he4__ni56__approx = 96,
        k_ni56__fe52_he4__approx = 97,
        NumRates = k_ni56__fe52_he4__approx
    };

    // number of reaclib rates

    const int NrateReaclib = 85;
    const int NumReaclibSets = 161;

    // number of tabular rates

    const int NrateTabular = 0;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "n13__p_c12",  // 1,
        "o16__he4_c12",  // 2,
        "f18__he4_n14",  // 3,
        "ne20__he4_o16",  // 4,
        "na22__p_ne21",  // 5,
        "na22__he4_f18",  // 6,
        "mg24__p_na23",  // 7,
        "mg24__he4_ne20",  // 8,
        "al27__he4_na23",  // 9,
        "si28__p_al27",  // 10,
        "si28__he4_mg24",  // 11,
        "p31__he4_al27",  // 12,
        "s32__p_p31",  // 13,
        "s32__he4_si28",  // 14,
        "c12__he4_he4_he4",  // 15,
        "p_c12__n13",  // 16,
        "he4_c12__o16",  // 17,
        "he4_n14__f18",  // 18,
        "he4_o16__ne20",  // 19,
        "he4_f18__na22",  // 20,
        "he4_ne20__mg24",  // 21,
        "p_ne21__na22",  // 22,
        "p_na23__mg24",  // 23,
        "he4_na23__al27",  // 24,
        "he4_mg24__si28",  // 25,
        "p_al27__si28",  // 26,
        "he4_al27__p31",  // 27,
        "he4_si28__s32",  // 28,
        "p_p31__s32",  // 29,
        "c12_c12__p_na23",  // 30,
        "c12_c12__he4_ne20",  // 31,
        "he4_n13__p_o16",  // 32,
        "p_o16__he4_n13",  // 33,
        "c12_o16__p_al27",  // 34,
        "c12_o16__he4_mg24",  // 35,
        "o16_o16__p_p31",  // 36,
        "o16_o16__he4_si28",  // 37,
        "he4_f18__p_ne21",  // 38,
        "he4_ne20__p_na23",  // 39,
        "p_ne21__he4_f18",  // 40,
        "p_na23__he4_ne20",  // 41,
        "he4_mg24__p_al27",  // 42,
        "p_al27__he4_mg24",  // 43,
        "he4_si28__p_p31",  // 44,
        "p_p31__he4_si28",  // 45,
        "he4_he4_he4__c12",  // 46,
        "c12_c12__mg24__modified",  // 47,
        "o16_o16__s32__modified",  // 48,
        "c12_o16__si28__modified",  // 49,
        "he4_s32__ar36",  // 50,
        "he4_s32__p_cl35",  // 51,
        "p_cl35__ar36",  // 52,
        "ar36__he4_s32",  // 53,
        "ar36__p_cl35",  // 54,
        "p_cl35__he4_s32",  // 55,
        "he4_ar36__ca40",  // 56,
        "he4_ar36__p_k39",  // 57,
        "p_k39__ca40",  // 58,
        "ca40__he4_ar36",  // 59,
        "ca40__p_k39",  // 60,
        "p_k39__he4_ar36",  // 61,
        "he4_ca40__ti44",  // 62,
        "he4_ca40__p_sc43",  // 63,
        "p_sc43__ti44",  // 64,
        "ti44__he4_ca40",  // 65,
        "ti44__p_sc43",  // 66,
        "p_sc43__he4_ca40",  // 67,
        "he4_ti44__cr48",  // 68,
        "he4_ti44__p_v47",  // 69,
        "p_v47__cr48",  // 70,
        "cr48__he4_ti44",  // 71,
        "cr48__p_v47",  // 72,
        "p_v47__he4_ti44",  // 73,
        "he4_cr48__fe52",  // 74,
        "he4_cr48__p_mn51",  // 75,
        "p_mn51__fe52",  // 76,
        "fe52__he4_cr48",  // 77,
        "fe52__p_mn51",  // 78,
        "p_mn51__he4_cr48",  // 79,
        "he4_fe52__ni56",  // 80,
        "he4_fe52__p_co55",  // 81,
        "p_co55__ni56",  // 82,
        "ni56__he4_fe52",  // 83,
        "ni56__p_co55",  // 84,
        "p_co55__he4_fe52",  // 85,
        "s32_he4__ar36__approx",  // 86,
        "ar36__s32_he4__approx",  // 87,
        "ar36_he4__ca40__approx",  // 88,
        "ca40__ar36_he4__approx",  // 89,
        "ca40_he4__ti44__approx",  // 90,
        "ti44__ca40_he4__approx",  // 91,
        "ti44_he4__cr48__approx",  // 92,
        "cr48__ti44_he4__approx",  // 93,
        "cr48_he4__fe52__approx",  // 94,
        "fe52__cr48_he4__approx",  // 95,
        "fe52_he4__ni56__approx",  // 96,
        "ni56__fe52_he4__approx"  // 97,
    };

}

#endif
