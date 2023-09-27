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
        k_c12_to_he4_he4_he4 = 15,
        k_p_c12_to_n13 = 16,
        k_he4_c12_to_o16 = 17,
        k_he4_n14_to_f18 = 18,
        k_he4_o16_to_ne20 = 19,
        k_he4_f18_to_na22 = 20,
        k_he4_ne20_to_mg24 = 21,
        k_p_ne21_to_na22 = 22,
        k_p_na23_to_mg24 = 23,
        k_he4_na23_to_al27 = 24,
        k_he4_mg24_to_si28 = 25,
        k_p_al27_to_si28 = 26,
        k_he4_al27_to_p31 = 27,
        k_he4_si28_to_s32 = 28,
        k_p_p31_to_s32 = 29,
        k_c12_c12_to_p_na23 = 30,
        k_c12_c12_to_he4_ne20 = 31,
        k_he4_n13_to_p_o16 = 32,
        k_p_o16_to_he4_n13 = 33,
        k_c12_o16_to_p_al27 = 34,
        k_c12_o16_to_he4_mg24 = 35,
        k_o16_o16_to_p_p31 = 36,
        k_o16_o16_to_he4_si28 = 37,
        k_he4_f18_to_p_ne21 = 38,
        k_he4_ne20_to_p_na23 = 39,
        k_he4_ne20_to_c12_c12 = 40,
        k_c12_ne20_to_p_p31 = 41,
        k_c12_ne20_to_he4_si28 = 42,
        k_p_ne21_to_he4_f18 = 43,
        k_p_na23_to_he4_ne20 = 44,
        k_p_na23_to_c12_c12 = 45,
        k_he4_mg24_to_p_al27 = 46,
        k_he4_mg24_to_c12_o16 = 47,
        k_p_al27_to_he4_mg24 = 48,
        k_p_al27_to_c12_o16 = 49,
        k_he4_si28_to_p_p31 = 50,
        k_he4_si28_to_c12_ne20 = 51,
        k_he4_si28_to_o16_o16 = 52,
        k_p_p31_to_he4_si28 = 53,
        k_p_p31_to_c12_ne20 = 54,
        k_p_p31_to_o16_o16 = 55,
        k_he4_he4_he4_to_c12 = 56,
        k_c12_c12_to_mg24_modified = 57,
        k_o16_o16_to_s32_modified = 58,
        k_c12_o16_to_si28_modified = 59,
        k_he4_s32_to_ar36_removed = 60,
        k_he4_s32_to_p_cl35_removed = 61,
        k_p_cl35_to_ar36_removed = 62,
        k_ar36_to_he4_s32_removed = 63,
        k_ar36_to_p_cl35_removed = 64,
        k_p_cl35_to_he4_s32_removed = 65,
        k_he4_ar36_to_ca40_removed = 66,
        k_he4_ar36_to_p_k39_removed = 67,
        k_p_k39_to_ca40_removed = 68,
        k_ca40_to_he4_ar36_removed = 69,
        k_ca40_to_p_k39_removed = 70,
        k_p_k39_to_he4_ar36_removed = 71,
        k_he4_ca40_to_ti44_removed = 72,
        k_he4_ca40_to_p_sc43_removed = 73,
        k_p_sc43_to_ti44_removed = 74,
        k_ti44_to_he4_ca40_removed = 75,
        k_ti44_to_p_sc43_removed = 76,
        k_p_sc43_to_he4_ca40_removed = 77,
        k_he4_ti44_to_cr48_removed = 78,
        k_he4_ti44_to_p_v47_removed = 79,
        k_p_v47_to_cr48_removed = 80,
        k_cr48_to_he4_ti44_removed = 81,
        k_cr48_to_p_v47_removed = 82,
        k_p_v47_to_he4_ti44_removed = 83,
        k_he4_cr48_to_fe52_removed = 84,
        k_he4_cr48_to_p_mn51_removed = 85,
        k_p_mn51_to_fe52_removed = 86,
        k_fe52_to_he4_cr48_removed = 87,
        k_fe52_to_p_mn51_removed = 88,
        k_p_mn51_to_he4_cr48_removed = 89,
        k_he4_fe52_to_ni56_removed = 90,
        k_he4_fe52_to_p_co55_removed = 91,
        k_p_co55_to_ni56_removed = 92,
        k_ni56_to_he4_fe52_removed = 93,
        k_ni56_to_p_co55_removed = 94,
        k_p_co55_to_he4_fe52_removed = 95,
        k_s32_he4_to_ar36_approx = 96,
        k_ar36_to_s32_he4_approx = 97,
        k_ar36_he4_to_ca40_approx = 98,
        k_ca40_to_ar36_he4_approx = 99,
        k_ca40_he4_to_ti44_approx = 100,
        k_ti44_to_ca40_he4_approx = 101,
        k_ti44_he4_to_cr48_approx = 102,
        k_cr48_to_ti44_he4_approx = 103,
        k_cr48_he4_to_fe52_approx = 104,
        k_fe52_to_cr48_he4_approx = 105,
        k_fe52_he4_to_ni56_approx = 106,
        k_ni56_to_fe52_he4_approx = 107,
        NumRates = k_ni56_to_fe52_he4_approx
    };

    // number of reaclib rates

    const int NrateReaclib = 95;

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
        "c12_to_he4_he4_he4",  // 15,
        "p_c12_to_n13",  // 16,
        "he4_c12_to_o16",  // 17,
        "he4_n14_to_f18",  // 18,
        "he4_o16_to_ne20",  // 19,
        "he4_f18_to_na22",  // 20,
        "he4_ne20_to_mg24",  // 21,
        "p_ne21_to_na22",  // 22,
        "p_na23_to_mg24",  // 23,
        "he4_na23_to_al27",  // 24,
        "he4_mg24_to_si28",  // 25,
        "p_al27_to_si28",  // 26,
        "he4_al27_to_p31",  // 27,
        "he4_si28_to_s32",  // 28,
        "p_p31_to_s32",  // 29,
        "c12_c12_to_p_na23",  // 30,
        "c12_c12_to_he4_ne20",  // 31,
        "he4_n13_to_p_o16",  // 32,
        "p_o16_to_he4_n13",  // 33,
        "c12_o16_to_p_al27",  // 34,
        "c12_o16_to_he4_mg24",  // 35,
        "o16_o16_to_p_p31",  // 36,
        "o16_o16_to_he4_si28",  // 37,
        "he4_f18_to_p_ne21",  // 38,
        "he4_ne20_to_p_na23",  // 39,
        "he4_ne20_to_c12_c12",  // 40,
        "c12_ne20_to_p_p31",  // 41,
        "c12_ne20_to_he4_si28",  // 42,
        "p_ne21_to_he4_f18",  // 43,
        "p_na23_to_he4_ne20",  // 44,
        "p_na23_to_c12_c12",  // 45,
        "he4_mg24_to_p_al27",  // 46,
        "he4_mg24_to_c12_o16",  // 47,
        "p_al27_to_he4_mg24",  // 48,
        "p_al27_to_c12_o16",  // 49,
        "he4_si28_to_p_p31",  // 50,
        "he4_si28_to_c12_ne20",  // 51,
        "he4_si28_to_o16_o16",  // 52,
        "p_p31_to_he4_si28",  // 53,
        "p_p31_to_c12_ne20",  // 54,
        "p_p31_to_o16_o16",  // 55,
        "he4_he4_he4_to_c12",  // 56,
        "c12_c12_to_mg24_modified",  // 57,
        "o16_o16_to_s32_modified",  // 58,
        "c12_o16_to_si28_modified",  // 59,
        "he4_s32_to_ar36_removed",  // 60,
        "he4_s32_to_p_cl35_removed",  // 61,
        "p_cl35_to_ar36_removed",  // 62,
        "ar36_to_he4_s32_removed",  // 63,
        "ar36_to_p_cl35_removed",  // 64,
        "p_cl35_to_he4_s32_removed",  // 65,
        "he4_ar36_to_ca40_removed",  // 66,
        "he4_ar36_to_p_k39_removed",  // 67,
        "p_k39_to_ca40_removed",  // 68,
        "ca40_to_he4_ar36_removed",  // 69,
        "ca40_to_p_k39_removed",  // 70,
        "p_k39_to_he4_ar36_removed",  // 71,
        "he4_ca40_to_ti44_removed",  // 72,
        "he4_ca40_to_p_sc43_removed",  // 73,
        "p_sc43_to_ti44_removed",  // 74,
        "ti44_to_he4_ca40_removed",  // 75,
        "ti44_to_p_sc43_removed",  // 76,
        "p_sc43_to_he4_ca40_removed",  // 77,
        "he4_ti44_to_cr48_removed",  // 78,
        "he4_ti44_to_p_v47_removed",  // 79,
        "p_v47_to_cr48_removed",  // 80,
        "cr48_to_he4_ti44_removed",  // 81,
        "cr48_to_p_v47_removed",  // 82,
        "p_v47_to_he4_ti44_removed",  // 83,
        "he4_cr48_to_fe52_removed",  // 84,
        "he4_cr48_to_p_mn51_removed",  // 85,
        "p_mn51_to_fe52_removed",  // 86,
        "fe52_to_he4_cr48_removed",  // 87,
        "fe52_to_p_mn51_removed",  // 88,
        "p_mn51_to_he4_cr48_removed",  // 89,
        "he4_fe52_to_ni56_removed",  // 90,
        "he4_fe52_to_p_co55_removed",  // 91,
        "p_co55_to_ni56_removed",  // 92,
        "ni56_to_he4_fe52_removed",  // 93,
        "ni56_to_p_co55_removed",  // 94,
        "p_co55_to_he4_fe52_removed",  // 95,
        "s32_he4_to_ar36_approx",  // 96,
        "ar36_to_s32_he4_approx",  // 97,
        "ar36_he4_to_ca40_approx",  // 98,
        "ca40_to_ar36_he4_approx",  // 99,
        "ca40_he4_to_ti44_approx",  // 100,
        "ti44_to_ca40_he4_approx",  // 101,
        "ti44_he4_to_cr48_approx",  // 102,
        "cr48_to_ti44_he4_approx",  // 103,
        "cr48_he4_to_fe52_approx",  // 104,
        "fe52_to_cr48_he4_approx",  // 105,
        "fe52_he4_to_ni56_approx",  // 106,
        "ni56_to_fe52_he4_approx"  // 107,
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
