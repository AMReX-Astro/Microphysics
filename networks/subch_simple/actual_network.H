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
        k_si28__p_al27 = 9,
        k_si28__he4_mg24 = 10,
        k_s32__p_p31 = 11,
        k_s32__he4_si28 = 12,
        k_c12__he4_he4_he4 = 13,
        k_p_c12__n13 = 14,
        k_he4_c12__o16 = 15,
        k_he4_n14__f18 = 16,
        k_he4_o16__ne20 = 17,
        k_he4_f18__na22 = 18,
        k_he4_ne20__mg24 = 19,
        k_p_ne21__na22 = 20,
        k_p_na23__mg24 = 21,
        k_he4_mg24__si28 = 22,
        k_p_al27__si28 = 23,
        k_he4_si28__s32 = 24,
        k_p_p31__s32 = 25,
        k_c12_c12__p_na23 = 26,
        k_c12_c12__he4_ne20 = 27,
        k_he4_n13__p_o16 = 28,
        k_p_o16__he4_n13 = 29,
        k_c12_o16__p_al27 = 30,
        k_c12_o16__he4_mg24 = 31,
        k_o16_o16__p_p31 = 32,
        k_o16_o16__he4_si28 = 33,
        k_he4_f18__p_ne21 = 34,
        k_he4_ne20__p_na23 = 35,
        k_p_ne21__he4_f18 = 36,
        k_p_na23__he4_ne20 = 37,
        k_he4_mg24__p_al27 = 38,
        k_p_al27__he4_mg24 = 39,
        k_he4_si28__p_p31 = 40,
        k_p_p31__he4_si28 = 41,
        k_he4_he4_he4__c12 = 42,
        k_c12_c12__mg24__modified = 43,
        k_o16_o16__s32__modified = 44,
        k_c12_o16__si28__modified = 45,
        k_he4_s32__ar36__removed = 46,
        k_he4_s32__p_cl35__removed = 47,
        k_p_cl35__ar36__removed = 48,
        k_ar36__he4_s32__removed = 49,
        k_ar36__p_cl35__removed = 50,
        k_p_cl35__he4_s32__removed = 51,
        k_he4_ar36__ca40__removed = 52,
        k_he4_ar36__p_k39__removed = 53,
        k_p_k39__ca40__removed = 54,
        k_ca40__he4_ar36__removed = 55,
        k_ca40__p_k39__removed = 56,
        k_p_k39__he4_ar36__removed = 57,
        k_he4_ca40__ti44__removed = 58,
        k_he4_ca40__p_sc43__removed = 59,
        k_p_sc43__ti44__removed = 60,
        k_ti44__he4_ca40__removed = 61,
        k_ti44__p_sc43__removed = 62,
        k_p_sc43__he4_ca40__removed = 63,
        k_he4_ti44__cr48__removed = 64,
        k_he4_ti44__p_v47__removed = 65,
        k_p_v47__cr48__removed = 66,
        k_cr48__he4_ti44__removed = 67,
        k_cr48__p_v47__removed = 68,
        k_p_v47__he4_ti44__removed = 69,
        k_he4_cr48__fe52__removed = 70,
        k_he4_cr48__p_mn51__removed = 71,
        k_p_mn51__fe52__removed = 72,
        k_fe52__he4_cr48__removed = 73,
        k_fe52__p_mn51__removed = 74,
        k_p_mn51__he4_cr48__removed = 75,
        k_he4_fe52__ni56__removed = 76,
        k_he4_fe52__p_co55__removed = 77,
        k_p_co55__ni56__removed = 78,
        k_ni56__he4_fe52__removed = 79,
        k_ni56__p_co55__removed = 80,
        k_p_co55__he4_fe52__removed = 81,
        k_s32_he4__ar36__approx = 82,
        k_ar36__s32_he4__approx = 83,
        k_ar36_he4__ca40__approx = 84,
        k_ca40__ar36_he4__approx = 85,
        k_ca40_he4__ti44__approx = 86,
        k_ti44__ca40_he4__approx = 87,
        k_ti44_he4__cr48__approx = 88,
        k_cr48__ti44_he4__approx = 89,
        k_cr48_he4__fe52__approx = 90,
        k_fe52__cr48_he4__approx = 91,
        k_fe52_he4__ni56__approx = 92,
        k_ni56__fe52_he4__approx = 93,
        NumRates = k_ni56__fe52_he4__approx
    };

    // number of reaclib rates

    const int NrateReaclib = 81;

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
        "si28__p_al27",  // 9,
        "si28__he4_mg24",  // 10,
        "s32__p_p31",  // 11,
        "s32__he4_si28",  // 12,
        "c12__he4_he4_he4",  // 13,
        "p_c12__n13",  // 14,
        "he4_c12__o16",  // 15,
        "he4_n14__f18",  // 16,
        "he4_o16__ne20",  // 17,
        "he4_f18__na22",  // 18,
        "he4_ne20__mg24",  // 19,
        "p_ne21__na22",  // 20,
        "p_na23__mg24",  // 21,
        "he4_mg24__si28",  // 22,
        "p_al27__si28",  // 23,
        "he4_si28__s32",  // 24,
        "p_p31__s32",  // 25,
        "c12_c12__p_na23",  // 26,
        "c12_c12__he4_ne20",  // 27,
        "he4_n13__p_o16",  // 28,
        "p_o16__he4_n13",  // 29,
        "c12_o16__p_al27",  // 30,
        "c12_o16__he4_mg24",  // 31,
        "o16_o16__p_p31",  // 32,
        "o16_o16__he4_si28",  // 33,
        "he4_f18__p_ne21",  // 34,
        "he4_ne20__p_na23",  // 35,
        "p_ne21__he4_f18",  // 36,
        "p_na23__he4_ne20",  // 37,
        "he4_mg24__p_al27",  // 38,
        "p_al27__he4_mg24",  // 39,
        "he4_si28__p_p31",  // 40,
        "p_p31__he4_si28",  // 41,
        "he4_he4_he4__c12",  // 42,
        "c12_c12__mg24__modified",  // 43,
        "o16_o16__s32__modified",  // 44,
        "c12_o16__si28__modified",  // 45,
        "he4_s32__ar36__removed",  // 46,
        "he4_s32__p_cl35__removed",  // 47,
        "p_cl35__ar36__removed",  // 48,
        "ar36__he4_s32__removed",  // 49,
        "ar36__p_cl35__removed",  // 50,
        "p_cl35__he4_s32__removed",  // 51,
        "he4_ar36__ca40__removed",  // 52,
        "he4_ar36__p_k39__removed",  // 53,
        "p_k39__ca40__removed",  // 54,
        "ca40__he4_ar36__removed",  // 55,
        "ca40__p_k39__removed",  // 56,
        "p_k39__he4_ar36__removed",  // 57,
        "he4_ca40__ti44__removed",  // 58,
        "he4_ca40__p_sc43__removed",  // 59,
        "p_sc43__ti44__removed",  // 60,
        "ti44__he4_ca40__removed",  // 61,
        "ti44__p_sc43__removed",  // 62,
        "p_sc43__he4_ca40__removed",  // 63,
        "he4_ti44__cr48__removed",  // 64,
        "he4_ti44__p_v47__removed",  // 65,
        "p_v47__cr48__removed",  // 66,
        "cr48__he4_ti44__removed",  // 67,
        "cr48__p_v47__removed",  // 68,
        "p_v47__he4_ti44__removed",  // 69,
        "he4_cr48__fe52__removed",  // 70,
        "he4_cr48__p_mn51__removed",  // 71,
        "p_mn51__fe52__removed",  // 72,
        "fe52__he4_cr48__removed",  // 73,
        "fe52__p_mn51__removed",  // 74,
        "p_mn51__he4_cr48__removed",  // 75,
        "he4_fe52__ni56__removed",  // 76,
        "he4_fe52__p_co55__removed",  // 77,
        "p_co55__ni56__removed",  // 78,
        "ni56__he4_fe52__removed",  // 79,
        "ni56__p_co55__removed",  // 80,
        "p_co55__he4_fe52__removed",  // 81,
        "s32_he4__ar36__approx",  // 82,
        "ar36__s32_he4__approx",  // 83,
        "ar36_he4__ca40__approx",  // 84,
        "ca40__ar36_he4__approx",  // 85,
        "ca40_he4__ti44__approx",  // 86,
        "ti44__ca40_he4__approx",  // 87,
        "ti44_he4__cr48__approx",  // 88,
        "cr48__ti44_he4__approx",  // 89,
        "cr48_he4__fe52__approx",  // 90,
        "fe52__cr48_he4__approx",  // 91,
        "fe52_he4__ni56__approx",  // 92,
        "ni56__fe52_he4__approx"  // 93,
    };

}

#endif
