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
        k_n13_to_p_c12 = 1,
        k_o16_to_he4_c12 = 2,
        k_f18_to_he4_n14 = 3,
        k_ne20_to_he4_o16 = 4,
        k_na22_to_p_ne21 = 5,
        k_na22_to_he4_f18 = 6,
        k_mg24_to_p_na23 = 7,
        k_mg24_to_he4_ne20 = 8,
        k_al27_to_he4_na23 = 9,
        k_si28_to_p_al27 = 10,
        k_si28_to_he4_mg24 = 11,
        k_p31_to_he4_al27 = 12,
        k_s32_to_p_p31 = 13,
        k_s32_to_he4_si28 = 14,
        k_cl35_to_he4_p31 = 15,
        k_ar36_to_p_cl35 = 16,
        k_ar36_to_he4_s32 = 17,
        k_k39_to_he4_cl35 = 18,
        k_ca40_to_p_k39 = 19,
        k_ca40_to_he4_ar36 = 20,
        k_sc43_to_he4_k39 = 21,
        k_ti44_to_p_sc43 = 22,
        k_ti44_to_he4_ca40 = 23,
        k_v47_to_he4_sc43 = 24,
        k_cr48_to_p_v47 = 25,
        k_cr48_to_he4_ti44 = 26,
        k_mn51_to_he4_v47 = 27,
        k_fe52_to_p_mn51 = 28,
        k_fe52_to_he4_cr48 = 29,
        k_co55_to_he4_mn51 = 30,
        k_ni56_to_p_co55 = 31,
        k_ni56_to_he4_fe52 = 32,
        k_c12_to_he4_he4_he4 = 33,
        k_p_c12_to_n13 = 34,
        k_he4_c12_to_o16 = 35,
        k_he4_n14_to_f18 = 36,
        k_he4_o16_to_ne20 = 37,
        k_he4_f18_to_na22 = 38,
        k_he4_ne20_to_mg24 = 39,
        k_p_ne21_to_na22 = 40,
        k_p_na23_to_mg24 = 41,
        k_he4_na23_to_al27 = 42,
        k_he4_mg24_to_si28 = 43,
        k_p_al27_to_si28 = 44,
        k_he4_al27_to_p31 = 45,
        k_he4_si28_to_s32 = 46,
        k_p_p31_to_s32 = 47,
        k_he4_p31_to_cl35 = 48,
        k_he4_s32_to_ar36 = 49,
        k_p_cl35_to_ar36 = 50,
        k_he4_cl35_to_k39 = 51,
        k_he4_ar36_to_ca40 = 52,
        k_p_k39_to_ca40 = 53,
        k_he4_k39_to_sc43 = 54,
        k_he4_ca40_to_ti44 = 55,
        k_p_sc43_to_ti44 = 56,
        k_he4_sc43_to_v47 = 57,
        k_he4_ti44_to_cr48 = 58,
        k_p_v47_to_cr48 = 59,
        k_he4_v47_to_mn51 = 60,
        k_he4_cr48_to_fe52 = 61,
        k_p_mn51_to_fe52 = 62,
        k_he4_mn51_to_co55 = 63,
        k_he4_fe52_to_ni56 = 64,
        k_p_co55_to_ni56 = 65,
        k_c12_c12_to_p_na23 = 66,
        k_c12_c12_to_he4_ne20 = 67,
        k_he4_n13_to_p_o16 = 68,
        k_p_o16_to_he4_n13 = 69,
        k_c12_o16_to_p_al27 = 70,
        k_c12_o16_to_he4_mg24 = 71,
        k_o16_o16_to_p_p31 = 72,
        k_o16_o16_to_he4_si28 = 73,
        k_he4_f18_to_p_ne21 = 74,
        k_he4_ne20_to_p_na23 = 75,
        k_he4_ne20_to_c12_c12 = 76,
        k_c12_ne20_to_p_p31 = 77,
        k_c12_ne20_to_he4_si28 = 78,
        k_p_ne21_to_he4_f18 = 79,
        k_p_na23_to_he4_ne20 = 80,
        k_p_na23_to_c12_c12 = 81,
        k_he4_mg24_to_p_al27 = 82,
        k_he4_mg24_to_c12_o16 = 83,
        k_p_al27_to_he4_mg24 = 84,
        k_p_al27_to_c12_o16 = 85,
        k_he4_si28_to_p_p31 = 86,
        k_he4_si28_to_c12_ne20 = 87,
        k_he4_si28_to_o16_o16 = 88,
        k_p_p31_to_he4_si28 = 89,
        k_p_p31_to_c12_ne20 = 90,
        k_p_p31_to_o16_o16 = 91,
        k_he4_s32_to_p_cl35 = 92,
        k_p_cl35_to_he4_s32 = 93,
        k_he4_ar36_to_p_k39 = 94,
        k_p_k39_to_he4_ar36 = 95,
        k_he4_ca40_to_p_sc43 = 96,
        k_p_sc43_to_he4_ca40 = 97,
        k_he4_ti44_to_p_v47 = 98,
        k_p_v47_to_he4_ti44 = 99,
        k_he4_cr48_to_p_mn51 = 100,
        k_p_mn51_to_he4_cr48 = 101,
        k_he4_fe52_to_p_co55 = 102,
        k_p_co55_to_he4_fe52 = 103,
        k_he4_he4_he4_to_c12 = 104,
        k_c12_c12_to_mg24_modified = 105,
        k_o16_o16_to_s32_modified = 106,
        k_c12_o16_to_si28_modified = 107,
        NumRates = k_c12_o16_to_si28_modified
    };

    // number of reaclib rates

    const int NrateReaclib = 107;

    // number of tabular rates

    const int NrateTabular = 0;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "n13_to_p_c12",  // 1,
        "o16_to_he4_c12",  // 2,
        "f18_to_he4_n14",  // 3,
        "ne20_to_he4_o16",  // 4,
        "na22_to_p_ne21",  // 5,
        "na22_to_he4_f18",  // 6,
        "mg24_to_p_na23",  // 7,
        "mg24_to_he4_ne20",  // 8,
        "al27_to_he4_na23",  // 9,
        "si28_to_p_al27",  // 10,
        "si28_to_he4_mg24",  // 11,
        "p31_to_he4_al27",  // 12,
        "s32_to_p_p31",  // 13,
        "s32_to_he4_si28",  // 14,
        "cl35_to_he4_p31",  // 15,
        "ar36_to_p_cl35",  // 16,
        "ar36_to_he4_s32",  // 17,
        "k39_to_he4_cl35",  // 18,
        "ca40_to_p_k39",  // 19,
        "ca40_to_he4_ar36",  // 20,
        "sc43_to_he4_k39",  // 21,
        "ti44_to_p_sc43",  // 22,
        "ti44_to_he4_ca40",  // 23,
        "v47_to_he4_sc43",  // 24,
        "cr48_to_p_v47",  // 25,
        "cr48_to_he4_ti44",  // 26,
        "mn51_to_he4_v47",  // 27,
        "fe52_to_p_mn51",  // 28,
        "fe52_to_he4_cr48",  // 29,
        "co55_to_he4_mn51",  // 30,
        "ni56_to_p_co55",  // 31,
        "ni56_to_he4_fe52",  // 32,
        "c12_to_he4_he4_he4",  // 33,
        "p_c12_to_n13",  // 34,
        "he4_c12_to_o16",  // 35,
        "he4_n14_to_f18",  // 36,
        "he4_o16_to_ne20",  // 37,
        "he4_f18_to_na22",  // 38,
        "he4_ne20_to_mg24",  // 39,
        "p_ne21_to_na22",  // 40,
        "p_na23_to_mg24",  // 41,
        "he4_na23_to_al27",  // 42,
        "he4_mg24_to_si28",  // 43,
        "p_al27_to_si28",  // 44,
        "he4_al27_to_p31",  // 45,
        "he4_si28_to_s32",  // 46,
        "p_p31_to_s32",  // 47,
        "he4_p31_to_cl35",  // 48,
        "he4_s32_to_ar36",  // 49,
        "p_cl35_to_ar36",  // 50,
        "he4_cl35_to_k39",  // 51,
        "he4_ar36_to_ca40",  // 52,
        "p_k39_to_ca40",  // 53,
        "he4_k39_to_sc43",  // 54,
        "he4_ca40_to_ti44",  // 55,
        "p_sc43_to_ti44",  // 56,
        "he4_sc43_to_v47",  // 57,
        "he4_ti44_to_cr48",  // 58,
        "p_v47_to_cr48",  // 59,
        "he4_v47_to_mn51",  // 60,
        "he4_cr48_to_fe52",  // 61,
        "p_mn51_to_fe52",  // 62,
        "he4_mn51_to_co55",  // 63,
        "he4_fe52_to_ni56",  // 64,
        "p_co55_to_ni56",  // 65,
        "c12_c12_to_p_na23",  // 66,
        "c12_c12_to_he4_ne20",  // 67,
        "he4_n13_to_p_o16",  // 68,
        "p_o16_to_he4_n13",  // 69,
        "c12_o16_to_p_al27",  // 70,
        "c12_o16_to_he4_mg24",  // 71,
        "o16_o16_to_p_p31",  // 72,
        "o16_o16_to_he4_si28",  // 73,
        "he4_f18_to_p_ne21",  // 74,
        "he4_ne20_to_p_na23",  // 75,
        "he4_ne20_to_c12_c12",  // 76,
        "c12_ne20_to_p_p31",  // 77,
        "c12_ne20_to_he4_si28",  // 78,
        "p_ne21_to_he4_f18",  // 79,
        "p_na23_to_he4_ne20",  // 80,
        "p_na23_to_c12_c12",  // 81,
        "he4_mg24_to_p_al27",  // 82,
        "he4_mg24_to_c12_o16",  // 83,
        "p_al27_to_he4_mg24",  // 84,
        "p_al27_to_c12_o16",  // 85,
        "he4_si28_to_p_p31",  // 86,
        "he4_si28_to_c12_ne20",  // 87,
        "he4_si28_to_o16_o16",  // 88,
        "p_p31_to_he4_si28",  // 89,
        "p_p31_to_c12_ne20",  // 90,
        "p_p31_to_o16_o16",  // 91,
        "he4_s32_to_p_cl35",  // 92,
        "p_cl35_to_he4_s32",  // 93,
        "he4_ar36_to_p_k39",  // 94,
        "p_k39_to_he4_ar36",  // 95,
        "he4_ca40_to_p_sc43",  // 96,
        "p_sc43_to_he4_ca40",  // 97,
        "he4_ti44_to_p_v47",  // 98,
        "p_v47_to_he4_ti44",  // 99,
        "he4_cr48_to_p_mn51",  // 100,
        "p_mn51_to_he4_cr48",  // 101,
        "he4_fe52_to_p_co55",  // 102,
        "p_co55_to_he4_fe52",  // 103,
        "he4_he4_he4_to_c12",  // 104,
        "c12_c12_to_mg24_modified",  // 105,
        "o16_o16_to_s32_modified",  // 106,
        "c12_o16_to_si28_modified"  // 107,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int h1_index = 0;
    constexpr int n_index = -1;
    constexpr int he4_index = 1;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
