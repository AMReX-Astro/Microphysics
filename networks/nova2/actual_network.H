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
        k_B8_to_He4_He4_weak_wc12 = 5,
        k_p_p_to_d_weak_bet_pos_ = 6,
        k_p_p_to_d_weak_electron_capture = 7,
        k_p_d_to_He3 = 8,
        k_d_d_to_He4 = 9,
        k_p_He3_to_He4_weak_bet_pos_ = 10,
        k_He4_He3_to_Be7 = 11,
        k_p_Be7_to_B8 = 12,
        k_p_C12_to_N13 = 13,
        k_He4_C12_to_O16 = 14,
        k_p_C13_to_N14 = 15,
        k_p_N13_to_O14 = 16,
        k_p_N14_to_O15 = 17,
        k_He4_N14_to_F18 = 18,
        k_p_N15_to_O16 = 19,
        k_p_O16_to_F17 = 20,
        k_p_O17_to_F18 = 21,
        k_d_He3_to_p_He4 = 22,
        k_He4_N13_to_p_O16 = 23,
        k_p_N15_to_He4_C12 = 24,
        k_He4_O14_to_p_F17 = 25,
        k_p_O17_to_He4_N14 = 26,
        k_p_F18_to_He4_O15 = 27,
        k_He3_He3_to_p_p_He4 = 28,
        k_d_Be7_to_p_He4_He4 = 29,
        k_He3_Be7_to_p_p_He4_He4 = 30,
        k_He4_He4_He4_to_C12 = 31,
        NumRates = k_He4_He4_He4_to_C12
    };

    // number of reaclib rates

    const int NrateReaclib = 31;

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
        "B8_to_He4_He4_weak_wc12",  // 5,
        "p_p_to_d_weak_bet_pos_",  // 6,
        "p_p_to_d_weak_electron_capture",  // 7,
        "p_d_to_He3",  // 8,
        "d_d_to_He4",  // 9,
        "p_He3_to_He4_weak_bet_pos_",  // 10,
        "He4_He3_to_Be7",  // 11,
        "p_Be7_to_B8",  // 12,
        "p_C12_to_N13",  // 13,
        "He4_C12_to_O16",  // 14,
        "p_C13_to_N14",  // 15,
        "p_N13_to_O14",  // 16,
        "p_N14_to_O15",  // 17,
        "He4_N14_to_F18",  // 18,
        "p_N15_to_O16",  // 19,
        "p_O16_to_F17",  // 20,
        "p_O17_to_F18",  // 21,
        "d_He3_to_p_He4",  // 22,
        "He4_N13_to_p_O16",  // 23,
        "p_N15_to_He4_C12",  // 24,
        "He4_O14_to_p_F17",  // 25,
        "p_O17_to_He4_N14",  // 26,
        "p_F18_to_He4_O15",  // 27,
        "He3_He3_to_p_p_He4",  // 28,
        "d_Be7_to_p_He4_He4",  // 29,
        "He3_Be7_to_p_p_He4_He4",  // 30,
        "He4_He4_He4_to_C12"  // 31,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int H1_index = 0;
    constexpr int N_index = -1;
    constexpr int He4_index = 3;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
