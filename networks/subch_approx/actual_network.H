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
        k_N13_to_p_C12 = 1,
        k_O16_to_He4_C12 = 2,
        k_F18_to_He4_N14 = 3,
        k_Ne20_to_He4_O16 = 4,
        k_Na22_to_p_Ne21 = 5,
        k_Na22_to_He4_F18 = 6,
        k_Mg24_to_p_Na23 = 7,
        k_Mg24_to_He4_Ne20 = 8,
        k_Al27_to_He4_Na23 = 9,
        k_Si28_to_p_Al27 = 10,
        k_Si28_to_He4_Mg24 = 11,
        k_P31_to_He4_Al27 = 12,
        k_S32_to_p_P31 = 13,
        k_S32_to_He4_Si28 = 14,
        k_C12_to_He4_He4_He4 = 15,
        k_p_C12_to_N13 = 16,
        k_He4_C12_to_O16 = 17,
        k_He4_N14_to_F18 = 18,
        k_He4_O16_to_Ne20 = 19,
        k_He4_F18_to_Na22 = 20,
        k_He4_Ne20_to_Mg24 = 21,
        k_p_Ne21_to_Na22 = 22,
        k_p_Na23_to_Mg24 = 23,
        k_He4_Na23_to_Al27 = 24,
        k_He4_Mg24_to_Si28 = 25,
        k_p_Al27_to_Si28 = 26,
        k_He4_Al27_to_P31 = 27,
        k_He4_Si28_to_S32 = 28,
        k_p_P31_to_S32 = 29,
        k_C12_C12_to_p_Na23 = 30,
        k_C12_C12_to_He4_Ne20 = 31,
        k_He4_N13_to_p_O16 = 32,
        k_p_O16_to_He4_N13 = 33,
        k_C12_O16_to_p_Al27 = 34,
        k_C12_O16_to_He4_Mg24 = 35,
        k_O16_O16_to_p_P31 = 36,
        k_O16_O16_to_He4_Si28 = 37,
        k_He4_F18_to_p_Ne21 = 38,
        k_He4_Ne20_to_p_Na23 = 39,
        k_He4_Ne20_to_C12_C12 = 40,
        k_C12_Ne20_to_p_P31 = 41,
        k_C12_Ne20_to_He4_Si28 = 42,
        k_p_Ne21_to_He4_F18 = 43,
        k_p_Na23_to_He4_Ne20 = 44,
        k_p_Na23_to_C12_C12 = 45,
        k_He4_Mg24_to_p_Al27 = 46,
        k_He4_Mg24_to_C12_O16 = 47,
        k_p_Al27_to_He4_Mg24 = 48,
        k_p_Al27_to_C12_O16 = 49,
        k_He4_Si28_to_p_P31 = 50,
        k_He4_Si28_to_C12_Ne20 = 51,
        k_He4_Si28_to_O16_O16 = 52,
        k_p_P31_to_He4_Si28 = 53,
        k_p_P31_to_C12_Ne20 = 54,
        k_p_P31_to_O16_O16 = 55,
        k_He4_He4_He4_to_C12 = 56,
        k_C12_C12_to_Mg24_modified = 57,
        k_O16_O16_to_S32_modified = 58,
        k_C12_O16_to_Si28_modified = 59,
        k_He4_S32_to_Ar36_removed = 60,
        k_He4_S32_to_p_Cl35_removed = 61,
        k_p_Cl35_to_Ar36_removed = 62,
        k_Ar36_to_He4_S32_removed = 63,
        k_Ar36_to_p_Cl35_removed = 64,
        k_p_Cl35_to_He4_S32_removed = 65,
        k_He4_Ar36_to_Ca40_removed = 66,
        k_He4_Ar36_to_p_K39_removed = 67,
        k_p_K39_to_Ca40_removed = 68,
        k_Ca40_to_He4_Ar36_removed = 69,
        k_Ca40_to_p_K39_removed = 70,
        k_p_K39_to_He4_Ar36_removed = 71,
        k_He4_Ca40_to_Ti44_removed = 72,
        k_He4_Ca40_to_p_Sc43_removed = 73,
        k_p_Sc43_to_Ti44_removed = 74,
        k_Ti44_to_He4_Ca40_removed = 75,
        k_Ti44_to_p_Sc43_removed = 76,
        k_p_Sc43_to_He4_Ca40_removed = 77,
        k_He4_Ti44_to_Cr48_removed = 78,
        k_He4_Ti44_to_p_V47_removed = 79,
        k_p_V47_to_Cr48_removed = 80,
        k_Cr48_to_He4_Ti44_removed = 81,
        k_Cr48_to_p_V47_removed = 82,
        k_p_V47_to_He4_Ti44_removed = 83,
        k_He4_Cr48_to_Fe52_removed = 84,
        k_He4_Cr48_to_p_Mn51_removed = 85,
        k_p_Mn51_to_Fe52_removed = 86,
        k_Fe52_to_He4_Cr48_removed = 87,
        k_Fe52_to_p_Mn51_removed = 88,
        k_p_Mn51_to_He4_Cr48_removed = 89,
        k_He4_Fe52_to_Ni56_removed = 90,
        k_He4_Fe52_to_p_Co55_removed = 91,
        k_p_Co55_to_Ni56_removed = 92,
        k_Ni56_to_He4_Fe52_removed = 93,
        k_Ni56_to_p_Co55_removed = 94,
        k_p_Co55_to_He4_Fe52_removed = 95,
        k_S32_He4_to_Ar36_approx = 96,
        k_Ar36_to_S32_He4_approx = 97,
        k_Ar36_He4_to_Ca40_approx = 98,
        k_Ca40_to_Ar36_He4_approx = 99,
        k_Ca40_He4_to_Ti44_approx = 100,
        k_Ti44_to_Ca40_He4_approx = 101,
        k_Ti44_He4_to_Cr48_approx = 102,
        k_Cr48_to_Ti44_He4_approx = 103,
        k_Cr48_He4_to_Fe52_approx = 104,
        k_Fe52_to_Cr48_He4_approx = 105,
        k_Fe52_He4_to_Ni56_approx = 106,
        k_Ni56_to_Fe52_He4_approx = 107,
        NumRates = k_Ni56_to_Fe52_He4_approx
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
        "N13_to_p_C12",  // 1,
        "O16_to_He4_C12",  // 2,
        "F18_to_He4_N14",  // 3,
        "Ne20_to_He4_O16",  // 4,
        "Na22_to_p_Ne21",  // 5,
        "Na22_to_He4_F18",  // 6,
        "Mg24_to_p_Na23",  // 7,
        "Mg24_to_He4_Ne20",  // 8,
        "Al27_to_He4_Na23",  // 9,
        "Si28_to_p_Al27",  // 10,
        "Si28_to_He4_Mg24",  // 11,
        "P31_to_He4_Al27",  // 12,
        "S32_to_p_P31",  // 13,
        "S32_to_He4_Si28",  // 14,
        "C12_to_He4_He4_He4",  // 15,
        "p_C12_to_N13",  // 16,
        "He4_C12_to_O16",  // 17,
        "He4_N14_to_F18",  // 18,
        "He4_O16_to_Ne20",  // 19,
        "He4_F18_to_Na22",  // 20,
        "He4_Ne20_to_Mg24",  // 21,
        "p_Ne21_to_Na22",  // 22,
        "p_Na23_to_Mg24",  // 23,
        "He4_Na23_to_Al27",  // 24,
        "He4_Mg24_to_Si28",  // 25,
        "p_Al27_to_Si28",  // 26,
        "He4_Al27_to_P31",  // 27,
        "He4_Si28_to_S32",  // 28,
        "p_P31_to_S32",  // 29,
        "C12_C12_to_p_Na23",  // 30,
        "C12_C12_to_He4_Ne20",  // 31,
        "He4_N13_to_p_O16",  // 32,
        "p_O16_to_He4_N13",  // 33,
        "C12_O16_to_p_Al27",  // 34,
        "C12_O16_to_He4_Mg24",  // 35,
        "O16_O16_to_p_P31",  // 36,
        "O16_O16_to_He4_Si28",  // 37,
        "He4_F18_to_p_Ne21",  // 38,
        "He4_Ne20_to_p_Na23",  // 39,
        "He4_Ne20_to_C12_C12",  // 40,
        "C12_Ne20_to_p_P31",  // 41,
        "C12_Ne20_to_He4_Si28",  // 42,
        "p_Ne21_to_He4_F18",  // 43,
        "p_Na23_to_He4_Ne20",  // 44,
        "p_Na23_to_C12_C12",  // 45,
        "He4_Mg24_to_p_Al27",  // 46,
        "He4_Mg24_to_C12_O16",  // 47,
        "p_Al27_to_He4_Mg24",  // 48,
        "p_Al27_to_C12_O16",  // 49,
        "He4_Si28_to_p_P31",  // 50,
        "He4_Si28_to_C12_Ne20",  // 51,
        "He4_Si28_to_O16_O16",  // 52,
        "p_P31_to_He4_Si28",  // 53,
        "p_P31_to_C12_Ne20",  // 54,
        "p_P31_to_O16_O16",  // 55,
        "He4_He4_He4_to_C12",  // 56,
        "C12_C12_to_Mg24_modified",  // 57,
        "O16_O16_to_S32_modified",  // 58,
        "C12_O16_to_Si28_modified",  // 59,
        "He4_S32_to_Ar36_removed",  // 60,
        "He4_S32_to_p_Cl35_removed",  // 61,
        "p_Cl35_to_Ar36_removed",  // 62,
        "Ar36_to_He4_S32_removed",  // 63,
        "Ar36_to_p_Cl35_removed",  // 64,
        "p_Cl35_to_He4_S32_removed",  // 65,
        "He4_Ar36_to_Ca40_removed",  // 66,
        "He4_Ar36_to_p_K39_removed",  // 67,
        "p_K39_to_Ca40_removed",  // 68,
        "Ca40_to_He4_Ar36_removed",  // 69,
        "Ca40_to_p_K39_removed",  // 70,
        "p_K39_to_He4_Ar36_removed",  // 71,
        "He4_Ca40_to_Ti44_removed",  // 72,
        "He4_Ca40_to_p_Sc43_removed",  // 73,
        "p_Sc43_to_Ti44_removed",  // 74,
        "Ti44_to_He4_Ca40_removed",  // 75,
        "Ti44_to_p_Sc43_removed",  // 76,
        "p_Sc43_to_He4_Ca40_removed",  // 77,
        "He4_Ti44_to_Cr48_removed",  // 78,
        "He4_Ti44_to_p_V47_removed",  // 79,
        "p_V47_to_Cr48_removed",  // 80,
        "Cr48_to_He4_Ti44_removed",  // 81,
        "Cr48_to_p_V47_removed",  // 82,
        "p_V47_to_He4_Ti44_removed",  // 83,
        "He4_Cr48_to_Fe52_removed",  // 84,
        "He4_Cr48_to_p_Mn51_removed",  // 85,
        "p_Mn51_to_Fe52_removed",  // 86,
        "Fe52_to_He4_Cr48_removed",  // 87,
        "Fe52_to_p_Mn51_removed",  // 88,
        "p_Mn51_to_He4_Cr48_removed",  // 89,
        "He4_Fe52_to_Ni56_removed",  // 90,
        "He4_Fe52_to_p_Co55_removed",  // 91,
        "p_Co55_to_Ni56_removed",  // 92,
        "Ni56_to_He4_Fe52_removed",  // 93,
        "Ni56_to_p_Co55_removed",  // 94,
        "p_Co55_to_He4_Fe52_removed",  // 95,
        "S32_He4_to_Ar36_approx",  // 96,
        "Ar36_to_S32_He4_approx",  // 97,
        "Ar36_He4_to_Ca40_approx",  // 98,
        "Ca40_to_Ar36_He4_approx",  // 99,
        "Ca40_He4_to_Ti44_approx",  // 100,
        "Ti44_to_Ca40_He4_approx",  // 101,
        "Ti44_He4_to_Cr48_approx",  // 102,
        "Cr48_to_Ti44_He4_approx",  // 103,
        "Cr48_He4_to_Fe52_approx",  // 104,
        "Fe52_to_Cr48_He4_approx",  // 105,
        "Fe52_He4_to_Ni56_approx",  // 106,
        "Ni56_to_Fe52_He4_approx"  // 107,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int H1_index = 0;
    constexpr int N_index = -1;
    constexpr int He4_index = 1;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
