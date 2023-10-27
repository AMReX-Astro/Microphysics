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
        k_N13_to_C13_weak_wc12 = 1,
        k_O14_to_N14_weak_wc12 = 2,
        k_O15_to_N15_weak_wc12 = 3,
        k_F17_to_O17_weak_wc12 = 4,
        k_F18_to_O18_weak_wc12 = 5,
        k_Ne18_to_F18_weak_wc12 = 6,
        k_Ne19_to_F19_weak_wc12 = 7,
        k_p_C12_to_N13 = 8,
        k_He4_C12_to_O16 = 9,
        k_p_C13_to_N14 = 10,
        k_p_N13_to_O14 = 11,
        k_p_N14_to_O15 = 12,
        k_He4_N14_to_F18 = 13,
        k_p_N15_to_O16 = 14,
        k_He4_N15_to_F19 = 15,
        k_He4_O14_to_Ne18 = 16,
        k_He4_O15_to_Ne19 = 17,
        k_p_O16_to_F17 = 18,
        k_He4_O16_to_Ne20 = 19,
        k_p_O17_to_F18 = 20,
        k_p_O18_to_F19 = 21,
        k_p_F17_to_Ne18 = 22,
        k_p_F18_to_Ne19 = 23,
        k_p_F19_to_Ne20 = 24,
        k_He4_Ne18_to_Mg22 = 25,
        k_He4_Ne20_to_Mg24 = 26,
        k_C12_C12_to_He4_Ne20 = 27,
        k_He4_N13_to_p_O16 = 28,
        k_p_N15_to_He4_C12 = 29,
        k_He4_O14_to_p_F17 = 30,
        k_C12_O16_to_He4_Mg24 = 31,
        k_p_O17_to_He4_N14 = 32,
        k_p_O18_to_He4_N15 = 33,
        k_p_F18_to_He4_O15 = 34,
        k_p_F19_to_He4_O16 = 35,
        k_p_Ne20_to_He4_F17 = 36,
        k_He4_He4_He4_to_C12 = 37,
        NumRates = k_He4_He4_He4_to_C12
    };

    // number of reaclib rates

    const int NrateReaclib = 37;

    // number of tabular rates

    const int NrateTabular = 0;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "N13_to_C13_weak_wc12",  // 1,
        "O14_to_N14_weak_wc12",  // 2,
        "O15_to_N15_weak_wc12",  // 3,
        "F17_to_O17_weak_wc12",  // 4,
        "F18_to_O18_weak_wc12",  // 5,
        "Ne18_to_F18_weak_wc12",  // 6,
        "Ne19_to_F19_weak_wc12",  // 7,
        "p_C12_to_N13",  // 8,
        "He4_C12_to_O16",  // 9,
        "p_C13_to_N14",  // 10,
        "p_N13_to_O14",  // 11,
        "p_N14_to_O15",  // 12,
        "He4_N14_to_F18",  // 13,
        "p_N15_to_O16",  // 14,
        "He4_N15_to_F19",  // 15,
        "He4_O14_to_Ne18",  // 16,
        "He4_O15_to_Ne19",  // 17,
        "p_O16_to_F17",  // 18,
        "He4_O16_to_Ne20",  // 19,
        "p_O17_to_F18",  // 20,
        "p_O18_to_F19",  // 21,
        "p_F17_to_Ne18",  // 22,
        "p_F18_to_Ne19",  // 23,
        "p_F19_to_Ne20",  // 24,
        "He4_Ne18_to_Mg22",  // 25,
        "He4_Ne20_to_Mg24",  // 26,
        "C12_C12_to_He4_Ne20",  // 27,
        "He4_N13_to_p_O16",  // 28,
        "p_N15_to_He4_C12",  // 29,
        "He4_O14_to_p_F17",  // 30,
        "C12_O16_to_He4_Mg24",  // 31,
        "p_O17_to_He4_N14",  // 32,
        "p_O18_to_He4_N15",  // 33,
        "p_F18_to_He4_O15",  // 34,
        "p_F19_to_He4_O16",  // 35,
        "p_Ne20_to_He4_F17",  // 36,
        "He4_He4_He4_to_C12"  // 37,
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
