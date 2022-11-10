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
        k_p_cl35__ar36 = 29,
        k_p_cl35__he4_s32 = 30,
        k_he4_ar36__ca40 = 31,
        k_p_k39__ca40 = 32,
        k_p_k39__he4_ar36 = 33,
        k_he4_ca40__ti44 = 34,
        k_p_sc43__ti44 = 35,
        k_p_sc43__he4_ca40 = 36,
        k_he4_ti44__cr48 = 37,
        k_he4_ti44__p_v47 = 38,
        k_p_v47__cr48 = 39,
        k_he4_cr48__fe52 = 40,
        k_he4_cr48__p_mn51 = 41,
        k_p_mn51__fe52 = 42,
        k_he4_fe52__ni56 = 43,
        k_he4_fe52__p_co55 = 44,
        k_p_co55__ni56 = 45,
        k_s32_he4__ar36__approx = 46,
        k_ar36__s32_he4__approx = 47,
        k_ar36_he4__ca40__approx = 48,
        k_ca40__ar36_he4__approx = 49,
        k_ca40_he4__ti44__approx = 50,
        k_ti44__ca40_he4__approx = 51,
        k_ti44_he4__cr48__approx = 52,
        k_cr48__ti44_he4__approx = 53,
        k_cr48_he4__fe52__approx = 54,
        k_fe52__cr48_he4__approx = 55,
        k_fe52_he4__ni56__approx = 56,
        k_ni56__fe52_he4__approx = 57,
        k_n13__p_c12__derived = 58,
        k_o16__he4_c12__derived = 59,
        k_f18__he4_n14__derived = 60,
        k_ne20__he4_o16__derived = 61,
        k_na22__p_ne21__derived = 62,
        k_na22__he4_f18__derived = 63,
        k_mg24__p_na23__derived = 64,
        k_mg24__he4_ne20__derived = 65,
        k_si28__p_al27__derived = 66,
        k_si28__he4_mg24__derived = 67,
        k_s32__p_p31__derived = 68,
        k_s32__he4_si28__derived = 69,
        k_c12__he4_he4_he4__derived = 70,
        k_p_o16__he4_n13__derived = 71,
        k_he4_ne20__p_na23__derived = 72,
        k_p_ne21__he4_f18__derived = 73,
        k_he4_mg24__p_al27__derived = 74,
        k_he4_si28__p_p31__derived = 75,
        k_he4_s32__p_cl35__derived = 76,
        k_ar36__he4_s32__derived = 77,
        k_ar36__p_cl35__derived = 78,
        k_he4_ar36__p_k39__derived = 79,
        k_ca40__he4_ar36__derived = 80,
        k_ca40__p_k39__derived = 81,
        k_he4_ca40__p_sc43__derived = 82,
        k_ti44__he4_ca40__derived = 83,
        k_ti44__p_sc43__derived = 84,
        k_cr48__he4_ti44__derived = 85,
        k_cr48__p_v47__derived = 86,
        k_p_v47__he4_ti44__derived = 87,
        k_fe52__he4_cr48__derived = 88,
        k_fe52__p_mn51__derived = 89,
        k_p_mn51__he4_cr48__derived = 90,
        k_ni56__he4_fe52__derived = 91,
        k_ni56__p_co55__derived = 92,
        k_p_co55__he4_fe52__derived = 93,
        NumRates = k_p_co55__he4_fe52__derived
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
        "p_cl35__ar36",  // 29,
        "p_cl35__he4_s32",  // 30,
        "he4_ar36__ca40",  // 31,
        "p_k39__ca40",  // 32,
        "p_k39__he4_ar36",  // 33,
        "he4_ca40__ti44",  // 34,
        "p_sc43__ti44",  // 35,
        "p_sc43__he4_ca40",  // 36,
        "he4_ti44__cr48",  // 37,
        "he4_ti44__p_v47",  // 38,
        "p_v47__cr48",  // 39,
        "he4_cr48__fe52",  // 40,
        "he4_cr48__p_mn51",  // 41,
        "p_mn51__fe52",  // 42,
        "he4_fe52__ni56",  // 43,
        "he4_fe52__p_co55",  // 44,
        "p_co55__ni56",  // 45,
        "s32_he4__ar36__approx",  // 46,
        "ar36__s32_he4__approx",  // 47,
        "ar36_he4__ca40__approx",  // 48,
        "ca40__ar36_he4__approx",  // 49,
        "ca40_he4__ti44__approx",  // 50,
        "ti44__ca40_he4__approx",  // 51,
        "ti44_he4__cr48__approx",  // 52,
        "cr48__ti44_he4__approx",  // 53,
        "cr48_he4__fe52__approx",  // 54,
        "fe52__cr48_he4__approx",  // 55,
        "fe52_he4__ni56__approx",  // 56,
        "ni56__fe52_he4__approx",  // 57,
        "n13__p_c12__derived",  // 58,
        "o16__he4_c12__derived",  // 59,
        "f18__he4_n14__derived",  // 60,
        "ne20__he4_o16__derived",  // 61,
        "na22__p_ne21__derived",  // 62,
        "na22__he4_f18__derived",  // 63,
        "mg24__p_na23__derived",  // 64,
        "mg24__he4_ne20__derived",  // 65,
        "si28__p_al27__derived",  // 66,
        "si28__he4_mg24__derived",  // 67,
        "s32__p_p31__derived",  // 68,
        "s32__he4_si28__derived",  // 69,
        "c12__he4_he4_he4__derived",  // 70,
        "p_o16__he4_n13__derived",  // 71,
        "he4_ne20__p_na23__derived",  // 72,
        "p_ne21__he4_f18__derived",  // 73,
        "he4_mg24__p_al27__derived",  // 74,
        "he4_si28__p_p31__derived",  // 75,
        "he4_s32__p_cl35__derived",  // 76,
        "ar36__he4_s32__derived",  // 77,
        "ar36__p_cl35__derived",  // 78,
        "he4_ar36__p_k39__derived",  // 79,
        "ca40__he4_ar36__derived",  // 80,
        "ca40__p_k39__derived",  // 81,
        "he4_ca40__p_sc43__derived",  // 82,
        "ti44__he4_ca40__derived",  // 83,
        "ti44__p_sc43__derived",  // 84,
        "cr48__he4_ti44__derived",  // 85,
        "cr48__p_v47__derived",  // 86,
        "p_v47__he4_ti44__derived",  // 87,
        "fe52__he4_cr48__derived",  // 88,
        "fe52__p_mn51__derived",  // 89,
        "p_mn51__he4_cr48__derived",  // 90,
        "ni56__he4_fe52__derived",  // 91,
        "ni56__p_co55__derived",  // 92,
        "p_co55__he4_fe52__derived"  // 93,
    };

}

#endif
