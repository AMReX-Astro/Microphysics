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
        k_n13_to_c13_weak_wc12 = 1,
        k_o14_to_n14_weak_wc12 = 2,
        k_o15_to_n15_weak_wc12 = 3,
        k_f17_to_o17_weak_wc12 = 4,
        k_b8_to_he4_he4_weak_wc12 = 5,
        k_p_p_to_d_weak_bet_pos_ = 6,
        k_p_p_to_d_weak_electron_capture = 7,
        k_p_d_to_he3 = 8,
        k_d_d_to_he4 = 9,
        k_p_he3_to_he4_weak_bet_pos_ = 10,
        k_he4_he3_to_be7 = 11,
        k_p_be7_to_b8 = 12,
        k_p_c12_to_n13 = 13,
        k_he4_c12_to_o16 = 14,
        k_p_c13_to_n14 = 15,
        k_p_n13_to_o14 = 16,
        k_p_n14_to_o15 = 17,
        k_he4_n14_to_f18 = 18,
        k_p_n15_to_o16 = 19,
        k_p_o16_to_f17 = 20,
        k_p_o17_to_f18 = 21,
        k_d_he3_to_p_he4 = 22,
        k_he4_n13_to_p_o16 = 23,
        k_p_n15_to_he4_c12 = 24,
        k_he4_o14_to_p_f17 = 25,
        k_p_o17_to_he4_n14 = 26,
        k_p_f18_to_he4_o15 = 27,
        k_he3_he3_to_p_p_he4 = 28,
        k_d_be7_to_p_he4_he4 = 29,
        k_he3_be7_to_p_p_he4_he4 = 30,
        k_he4_he4_he4_to_c12 = 31,
        NumRates = k_he4_he4_he4_to_c12
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
        "n13_to_c13_weak_wc12",  // 1,
        "o14_to_n14_weak_wc12",  // 2,
        "o15_to_n15_weak_wc12",  // 3,
        "f17_to_o17_weak_wc12",  // 4,
        "b8_to_he4_he4_weak_wc12",  // 5,
        "p_p_to_d_weak_bet_pos_",  // 6,
        "p_p_to_d_weak_electron_capture",  // 7,
        "p_d_to_he3",  // 8,
        "d_d_to_he4",  // 9,
        "p_he3_to_he4_weak_bet_pos_",  // 10,
        "he4_he3_to_be7",  // 11,
        "p_be7_to_b8",  // 12,
        "p_c12_to_n13",  // 13,
        "he4_c12_to_o16",  // 14,
        "p_c13_to_n14",  // 15,
        "p_n13_to_o14",  // 16,
        "p_n14_to_o15",  // 17,
        "he4_n14_to_f18",  // 18,
        "p_n15_to_o16",  // 19,
        "p_o16_to_f17",  // 20,
        "p_o17_to_f18",  // 21,
        "d_he3_to_p_he4",  // 22,
        "he4_n13_to_p_o16",  // 23,
        "p_n15_to_he4_c12",  // 24,
        "he4_o14_to_p_f17",  // 25,
        "p_o17_to_he4_n14",  // 26,
        "p_f18_to_he4_o15",  // 27,
        "he3_he3_to_p_p_he4",  // 28,
        "d_be7_to_p_he4_he4",  // 29,
        "he3_be7_to_p_p_he4_he4",  // 30,
        "he4_he4_he4_to_c12"  // 31,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int h1_index = 0;
    constexpr int n_index = -1;
    constexpr int he4_index = 3;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
