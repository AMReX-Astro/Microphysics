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
        k_p_c12__n13 = 1,
        k_he4_c12__o16 = 2,
        k_he4_n14__f18 = 3,
        k_he4_o16__ne20 = 4,
        k_he4_f18__na22 = 5,
        k_he4_ne20__mg24 = 6,
        k_p_ne21__na22 = 7,
        k_p_na23__mg24 = 8,
        k_he4_mg24__si28 = 9,
        k_p_al27__si28 = 10,
        k_he4_si28__s32 = 11,
        k_p_p31__s32 = 12,
        k_c12_c12__p_na23 = 13,
        k_c12_c12__he4_ne20 = 14,
        k_he4_n13__p_o16 = 15,
        k_c12_o16__p_al27 = 16,
        k_c12_o16__he4_mg24 = 17,
        k_o16_o16__p_p31 = 18,
        k_o16_o16__he4_si28 = 19,
        k_he4_f18__p_ne21 = 20,
        k_p_na23__he4_ne20 = 21,
        k_p_al27__he4_mg24 = 22,
        k_p_p31__he4_si28 = 23,
        k_he4_he4_he4__c12 = 24,
        k_c12_c12__mg24__modified = 25,
        k_o16_o16__s32__modified = 26,
        k_c12_o16__si28__modified = 27,
        k_he4_s32__ar36 = 28,
        k_he4_s32__p_cl35__derived = 29,
        k_p_cl35__ar36 = 30,
        k_ar36__he4_s32__derived = 31,
        k_ar36__p_cl35__derived = 32,
        k_p_cl35__he4_s32 = 33,
        k_he4_ar36__ca40 = 34,
        k_he4_ar36__p_k39__derived = 35,
        k_p_k39__ca40 = 36,
        k_ca40__he4_ar36__derived = 37,
        k_ca40__p_k39__derived = 38,
        k_p_k39__he4_ar36 = 39,
        k_he4_ca40__ti44 = 40,
        k_he4_ca40__p_sc43__derived = 41,
        k_p_sc43__ti44 = 42,
        k_ti44__he4_ca40__derived = 43,
        k_ti44__p_sc43__derived = 44,
        k_p_sc43__he4_ca40 = 45,
        k_he4_ti44__cr48 = 46,
        k_he4_ti44__p_v47 = 47,
        k_p_v47__cr48 = 48,
        k_cr48__he4_ti44__derived = 49,
        k_cr48__p_v47__derived = 50,
        k_p_v47__he4_ti44__derived = 51,
        k_he4_cr48__fe52 = 52,
        k_he4_cr48__p_mn51 = 53,
        k_p_mn51__fe52 = 54,
        k_fe52__he4_cr48__derived = 55,
        k_fe52__p_mn51__derived = 56,
        k_p_mn51__he4_cr48__derived = 57,
        k_he4_fe52__ni56 = 58,
        k_he4_fe52__p_co55 = 59,
        k_p_co55__ni56 = 60,
        k_ni56__he4_fe52__derived = 61,
        k_ni56__p_co55__derived = 62,
        k_p_co55__he4_fe52__derived = 63,
        k_s32_he4__ar36__approx = 64,
        k_ar36__s32_he4__approx = 65,
        k_ar36_he4__ca40__approx = 66,
        k_ca40__ar36_he4__approx = 67,
        k_ca40_he4__ti44__approx = 68,
        k_ti44__ca40_he4__approx = 69,
        k_ti44_he4__cr48__approx = 70,
        k_cr48__ti44_he4__approx = 71,
        k_cr48_he4__fe52__approx = 72,
        k_fe52__cr48_he4__approx = 73,
        k_fe52_he4__ni56__approx = 74,
        k_ni56__fe52_he4__approx = 75,
        k_n13__p_c12__derived = 76,
        k_o16__he4_c12__derived = 77,
        k_f18__he4_n14__derived = 78,
        k_ne20__he4_o16__derived = 79,
        k_na22__p_ne21__derived = 80,
        k_na22__he4_f18__derived = 81,
        k_mg24__p_na23__derived = 82,
        k_mg24__he4_ne20__derived = 83,
        k_si28__p_al27__derived = 84,
        k_si28__he4_mg24__derived = 85,
        k_s32__p_p31__derived = 86,
        k_s32__he4_si28__derived = 87,
        k_c12__he4_he4_he4__derived = 88,
        k_p_o16__he4_n13__derived = 89,
        k_he4_ne20__p_na23__derived = 90,
        k_p_ne21__he4_f18__derived = 91,
        k_he4_mg24__p_al27__derived = 92,
        k_he4_si28__p_p31__derived = 93,
        NumRates = k_he4_si28__p_p31__derived
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
        "p_c12__n13",  // 1,
        "he4_c12__o16",  // 2,
        "he4_n14__f18",  // 3,
        "he4_o16__ne20",  // 4,
        "he4_f18__na22",  // 5,
        "he4_ne20__mg24",  // 6,
        "p_ne21__na22",  // 7,
        "p_na23__mg24",  // 8,
        "he4_mg24__si28",  // 9,
        "p_al27__si28",  // 10,
        "he4_si28__s32",  // 11,
        "p_p31__s32",  // 12,
        "c12_c12__p_na23",  // 13,
        "c12_c12__he4_ne20",  // 14,
        "he4_n13__p_o16",  // 15,
        "c12_o16__p_al27",  // 16,
        "c12_o16__he4_mg24",  // 17,
        "o16_o16__p_p31",  // 18,
        "o16_o16__he4_si28",  // 19,
        "he4_f18__p_ne21",  // 20,
        "p_na23__he4_ne20",  // 21,
        "p_al27__he4_mg24",  // 22,
        "p_p31__he4_si28",  // 23,
        "he4_he4_he4__c12",  // 24,
        "c12_c12__mg24__modified",  // 25,
        "o16_o16__s32__modified",  // 26,
        "c12_o16__si28__modified",  // 27,
        "he4_s32__ar36",  // 28,
        "he4_s32__p_cl35__derived",  // 29,
        "p_cl35__ar36",  // 30,
        "ar36__he4_s32__derived",  // 31,
        "ar36__p_cl35__derived",  // 32,
        "p_cl35__he4_s32",  // 33,
        "he4_ar36__ca40",  // 34,
        "he4_ar36__p_k39__derived",  // 35,
        "p_k39__ca40",  // 36,
        "ca40__he4_ar36__derived",  // 37,
        "ca40__p_k39__derived",  // 38,
        "p_k39__he4_ar36",  // 39,
        "he4_ca40__ti44",  // 40,
        "he4_ca40__p_sc43__derived",  // 41,
        "p_sc43__ti44",  // 42,
        "ti44__he4_ca40__derived",  // 43,
        "ti44__p_sc43__derived",  // 44,
        "p_sc43__he4_ca40",  // 45,
        "he4_ti44__cr48",  // 46,
        "he4_ti44__p_v47",  // 47,
        "p_v47__cr48",  // 48,
        "cr48__he4_ti44__derived",  // 49,
        "cr48__p_v47__derived",  // 50,
        "p_v47__he4_ti44__derived",  // 51,
        "he4_cr48__fe52",  // 52,
        "he4_cr48__p_mn51",  // 53,
        "p_mn51__fe52",  // 54,
        "fe52__he4_cr48__derived",  // 55,
        "fe52__p_mn51__derived",  // 56,
        "p_mn51__he4_cr48__derived",  // 57,
        "he4_fe52__ni56",  // 58,
        "he4_fe52__p_co55",  // 59,
        "p_co55__ni56",  // 60,
        "ni56__he4_fe52__derived",  // 61,
        "ni56__p_co55__derived",  // 62,
        "p_co55__he4_fe52__derived",  // 63,
        "s32_he4__ar36__approx",  // 64,
        "ar36__s32_he4__approx",  // 65,
        "ar36_he4__ca40__approx",  // 66,
        "ca40__ar36_he4__approx",  // 67,
        "ca40_he4__ti44__approx",  // 68,
        "ti44__ca40_he4__approx",  // 69,
        "ti44_he4__cr48__approx",  // 70,
        "cr48__ti44_he4__approx",  // 71,
        "cr48_he4__fe52__approx",  // 72,
        "fe52__cr48_he4__approx",  // 73,
        "fe52_he4__ni56__approx",  // 74,
        "ni56__fe52_he4__approx",  // 75,
        "n13__p_c12__derived",  // 76,
        "o16__he4_c12__derived",  // 77,
        "f18__he4_n14__derived",  // 78,
        "ne20__he4_o16__derived",  // 79,
        "na22__p_ne21__derived",  // 80,
        "na22__he4_f18__derived",  // 81,
        "mg24__p_na23__derived",  // 82,
        "mg24__he4_ne20__derived",  // 83,
        "si28__p_al27__derived",  // 84,
        "si28__he4_mg24__derived",  // 85,
        "s32__p_p31__derived",  // 86,
        "s32__he4_si28__derived",  // 87,
        "c12__he4_he4_he4__derived",  // 88,
        "p_o16__he4_n13__derived",  // 89,
        "he4_ne20__p_na23__derived",  // 90,
        "p_ne21__he4_f18__derived",  // 91,
        "he4_mg24__p_al27__derived",  // 92,
        "he4_si28__p_p31__derived"  // 93,
    };

}

#endif
