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
        k_p_c12_to_n13 = 1,
        k_he4_c12_to_o16 = 2,
        k_he4_n14_to_f18 = 3,
        k_he4_o16_to_ne20 = 4,
        k_he4_f18_to_na22 = 5,
        k_he4_ne20_to_mg24 = 6,
        k_p_ne21_to_na22 = 7,
        k_n_na22_to_na23 = 8,
        k_p_na23_to_mg24 = 9,
        k_he4_mg24_to_si28 = 10,
        k_p_al27_to_si28 = 11,
        k_he4_si28_to_s32 = 12,
        k_p_p31_to_s32 = 13,
        k_c12_c12_to_p_na23 = 14,
        k_c12_c12_to_he4_ne20 = 15,
        k_he4_n13_to_p_o16 = 16,
        k_c12_o16_to_p_al27 = 17,
        k_c12_o16_to_he4_mg24 = 18,
        k_o16_o16_to_p_p31 = 19,
        k_o16_o16_to_he4_si28 = 20,
        k_he4_f18_to_p_ne21 = 21,
        k_he4_ne21_to_n_mg24 = 22,
        k_p_na23_to_he4_ne20 = 23,
        k_p_al27_to_he4_mg24 = 24,
        k_p_p31_to_he4_si28 = 25,
        k_he4_he4_he4_to_c12 = 26,
        k_c12_c12_to_mg24_modified = 27,
        k_o16_o16_to_s32_modified = 28,
        k_c12_o16_to_si28_modified = 29,
        k_he4_s32_to_ar36_removed = 30,
        k_p_cl35_to_ar36_removed = 31,
        k_p_cl35_to_he4_s32_removed = 32,
        k_he4_ar36_to_ca40_removed = 33,
        k_p_k39_to_ca40_removed = 34,
        k_p_k39_to_he4_ar36_removed = 35,
        k_he4_ca40_to_ti44_removed = 36,
        k_p_sc43_to_ti44_removed = 37,
        k_p_sc43_to_he4_ca40_removed = 38,
        k_he4_ti44_to_cr48_removed = 39,
        k_he4_ti44_to_p_v47_removed = 40,
        k_p_v47_to_cr48_removed = 41,
        k_he4_cr48_to_fe52_removed = 42,
        k_he4_cr48_to_p_mn51_removed = 43,
        k_p_mn51_to_fe52_removed = 44,
        k_he4_fe52_to_ni56_removed = 45,
        k_he4_fe52_to_p_co55_removed = 46,
        k_p_co55_to_ni56_removed = 47,
        k_s32_he4_to_ar36_approx = 48,
        k_ar36_to_s32_he4_approx = 49,
        k_ar36_he4_to_ca40_approx = 50,
        k_ca40_to_ar36_he4_approx = 51,
        k_ca40_he4_to_ti44_approx = 52,
        k_ti44_to_ca40_he4_approx = 53,
        k_ti44_he4_to_cr48_approx = 54,
        k_cr48_to_ti44_he4_approx = 55,
        k_cr48_he4_to_fe52_approx = 56,
        k_fe52_to_cr48_he4_approx = 57,
        k_fe52_he4_to_ni56_approx = 58,
        k_ni56_to_fe52_he4_approx = 59,
        k_n13_to_p_c12_derived = 60,
        k_o16_to_he4_c12_derived = 61,
        k_f18_to_he4_n14_derived = 62,
        k_ne20_to_he4_o16_derived = 63,
        k_na22_to_p_ne21_derived = 64,
        k_na22_to_he4_f18_derived = 65,
        k_na23_to_n_na22_derived = 66,
        k_mg24_to_p_na23_derived = 67,
        k_mg24_to_he4_ne20_derived = 68,
        k_si28_to_p_al27_derived = 69,
        k_si28_to_he4_mg24_derived = 70,
        k_s32_to_p_p31_derived = 71,
        k_s32_to_he4_si28_derived = 72,
        k_c12_to_he4_he4_he4_derived = 73,
        k_p_o16_to_he4_n13_derived = 74,
        k_he4_ne20_to_p_na23_derived = 75,
        k_p_ne21_to_he4_f18_derived = 76,
        k_n_mg24_to_he4_ne21_derived = 77,
        k_he4_mg24_to_p_al27_derived = 78,
        k_he4_si28_to_p_p31_derived = 79,
        k_he4_s32_to_p_cl35_derived_removed = 80,
        k_ar36_to_he4_s32_derived_removed = 81,
        k_ar36_to_p_cl35_derived_removed = 82,
        k_he4_ar36_to_p_k39_derived_removed = 83,
        k_ca40_to_he4_ar36_derived_removed = 84,
        k_ca40_to_p_k39_derived_removed = 85,
        k_he4_ca40_to_p_sc43_derived_removed = 86,
        k_ti44_to_he4_ca40_derived_removed = 87,
        k_ti44_to_p_sc43_derived_removed = 88,
        k_cr48_to_he4_ti44_derived_removed = 89,
        k_cr48_to_p_v47_derived_removed = 90,
        k_p_v47_to_he4_ti44_derived_removed = 91,
        k_fe52_to_he4_cr48_derived_removed = 92,
        k_fe52_to_p_mn51_derived_removed = 93,
        k_p_mn51_to_he4_cr48_derived_removed = 94,
        k_ni56_to_he4_fe52_derived_removed = 95,
        k_ni56_to_p_co55_derived_removed = 96,
        k_p_co55_to_he4_fe52_derived_removed = 97,
        NumRates = k_p_co55_to_he4_fe52_derived_removed
    };

    // number of reaclib rates

    const int NrateReaclib = 85;

    // number of tabular rates

    const int NrateTabular = 0;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "p_c12_to_n13",  // 1,
        "he4_c12_to_o16",  // 2,
        "he4_n14_to_f18",  // 3,
        "he4_o16_to_ne20",  // 4,
        "he4_f18_to_na22",  // 5,
        "he4_ne20_to_mg24",  // 6,
        "p_ne21_to_na22",  // 7,
        "n_na22_to_na23",  // 8,
        "p_na23_to_mg24",  // 9,
        "he4_mg24_to_si28",  // 10,
        "p_al27_to_si28",  // 11,
        "he4_si28_to_s32",  // 12,
        "p_p31_to_s32",  // 13,
        "c12_c12_to_p_na23",  // 14,
        "c12_c12_to_he4_ne20",  // 15,
        "he4_n13_to_p_o16",  // 16,
        "c12_o16_to_p_al27",  // 17,
        "c12_o16_to_he4_mg24",  // 18,
        "o16_o16_to_p_p31",  // 19,
        "o16_o16_to_he4_si28",  // 20,
        "he4_f18_to_p_ne21",  // 21,
        "he4_ne21_to_n_mg24",  // 22,
        "p_na23_to_he4_ne20",  // 23,
        "p_al27_to_he4_mg24",  // 24,
        "p_p31_to_he4_si28",  // 25,
        "he4_he4_he4_to_c12",  // 26,
        "c12_c12_to_mg24_modified",  // 27,
        "o16_o16_to_s32_modified",  // 28,
        "c12_o16_to_si28_modified",  // 29,
        "he4_s32_to_ar36_removed",  // 30,
        "p_cl35_to_ar36_removed",  // 31,
        "p_cl35_to_he4_s32_removed",  // 32,
        "he4_ar36_to_ca40_removed",  // 33,
        "p_k39_to_ca40_removed",  // 34,
        "p_k39_to_he4_ar36_removed",  // 35,
        "he4_ca40_to_ti44_removed",  // 36,
        "p_sc43_to_ti44_removed",  // 37,
        "p_sc43_to_he4_ca40_removed",  // 38,
        "he4_ti44_to_cr48_removed",  // 39,
        "he4_ti44_to_p_v47_removed",  // 40,
        "p_v47_to_cr48_removed",  // 41,
        "he4_cr48_to_fe52_removed",  // 42,
        "he4_cr48_to_p_mn51_removed",  // 43,
        "p_mn51_to_fe52_removed",  // 44,
        "he4_fe52_to_ni56_removed",  // 45,
        "he4_fe52_to_p_co55_removed",  // 46,
        "p_co55_to_ni56_removed",  // 47,
        "s32_he4_to_ar36_approx",  // 48,
        "ar36_to_s32_he4_approx",  // 49,
        "ar36_he4_to_ca40_approx",  // 50,
        "ca40_to_ar36_he4_approx",  // 51,
        "ca40_he4_to_ti44_approx",  // 52,
        "ti44_to_ca40_he4_approx",  // 53,
        "ti44_he4_to_cr48_approx",  // 54,
        "cr48_to_ti44_he4_approx",  // 55,
        "cr48_he4_to_fe52_approx",  // 56,
        "fe52_to_cr48_he4_approx",  // 57,
        "fe52_he4_to_ni56_approx",  // 58,
        "ni56_to_fe52_he4_approx",  // 59,
        "n13_to_p_c12_derived",  // 60,
        "o16_to_he4_c12_derived",  // 61,
        "f18_to_he4_n14_derived",  // 62,
        "ne20_to_he4_o16_derived",  // 63,
        "na22_to_p_ne21_derived",  // 64,
        "na22_to_he4_f18_derived",  // 65,
        "na23_to_n_na22_derived",  // 66,
        "mg24_to_p_na23_derived",  // 67,
        "mg24_to_he4_ne20_derived",  // 68,
        "si28_to_p_al27_derived",  // 69,
        "si28_to_he4_mg24_derived",  // 70,
        "s32_to_p_p31_derived",  // 71,
        "s32_to_he4_si28_derived",  // 72,
        "c12_to_he4_he4_he4_derived",  // 73,
        "p_o16_to_he4_n13_derived",  // 74,
        "he4_ne20_to_p_na23_derived",  // 75,
        "p_ne21_to_he4_f18_derived",  // 76,
        "n_mg24_to_he4_ne21_derived",  // 77,
        "he4_mg24_to_p_al27_derived",  // 78,
        "he4_si28_to_p_p31_derived",  // 79,
        "he4_s32_to_p_cl35_derived_removed",  // 80,
        "ar36_to_he4_s32_derived_removed",  // 81,
        "ar36_to_p_cl35_derived_removed",  // 82,
        "he4_ar36_to_p_k39_derived_removed",  // 83,
        "ca40_to_he4_ar36_derived_removed",  // 84,
        "ca40_to_p_k39_derived_removed",  // 85,
        "he4_ca40_to_p_sc43_derived_removed",  // 86,
        "ti44_to_he4_ca40_derived_removed",  // 87,
        "ti44_to_p_sc43_derived_removed",  // 88,
        "cr48_to_he4_ti44_derived_removed",  // 89,
        "cr48_to_p_v47_derived_removed",  // 90,
        "p_v47_to_he4_ti44_derived_removed",  // 91,
        "fe52_to_he4_cr48_derived_removed",  // 92,
        "fe52_to_p_mn51_derived_removed",  // 93,
        "p_mn51_to_he4_cr48_derived_removed",  // 94,
        "ni56_to_he4_fe52_derived_removed",  // 95,
        "ni56_to_p_co55_derived_removed",  // 96,
        "p_co55_to_he4_fe52_derived_removed"  // 97,
    };

}

#endif
