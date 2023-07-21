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
        k_f18_to_o18_weak_wc12 = 5,
        k_ne18_to_f18_weak_wc12 = 6,
        k_ne19_to_f19_weak_wc12 = 7,
        k_p_c12_to_n13 = 8,
        k_he4_c12_to_o16 = 9,
        k_p_c13_to_n14 = 10,
        k_p_n13_to_o14 = 11,
        k_p_n14_to_o15 = 12,
        k_he4_n14_to_f18 = 13,
        k_p_n15_to_o16 = 14,
        k_he4_n15_to_f19 = 15,
        k_he4_o14_to_ne18 = 16,
        k_he4_o15_to_ne19 = 17,
        k_p_o16_to_f17 = 18,
        k_he4_o16_to_ne20 = 19,
        k_p_o17_to_f18 = 20,
        k_p_o18_to_f19 = 21,
        k_p_f17_to_ne18 = 22,
        k_p_f18_to_ne19 = 23,
        k_p_f19_to_ne20 = 24,
        k_he4_ne18_to_mg22 = 25,
        k_he4_ne20_to_mg24 = 26,
        k_c12_c12_to_he4_ne20 = 27,
        k_he4_n13_to_p_o16 = 28,
        k_p_n15_to_he4_c12 = 29,
        k_he4_o14_to_p_f17 = 30,
        k_c12_o16_to_he4_mg24 = 31,
        k_p_o17_to_he4_n14 = 32,
        k_p_o18_to_he4_n15 = 33,
        k_p_f18_to_he4_o15 = 34,
        k_p_f19_to_he4_o16 = 35,
        k_p_ne20_to_he4_f17 = 36,
        k_he4_he4_he4_to_c12 = 37,
        NumRates = k_he4_he4_he4_to_c12
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
        "n13_to_c13_weak_wc12",  // 1,
        "o14_to_n14_weak_wc12",  // 2,
        "o15_to_n15_weak_wc12",  // 3,
        "f17_to_o17_weak_wc12",  // 4,
        "f18_to_o18_weak_wc12",  // 5,
        "ne18_to_f18_weak_wc12",  // 6,
        "ne19_to_f19_weak_wc12",  // 7,
        "p_c12_to_n13",  // 8,
        "he4_c12_to_o16",  // 9,
        "p_c13_to_n14",  // 10,
        "p_n13_to_o14",  // 11,
        "p_n14_to_o15",  // 12,
        "he4_n14_to_f18",  // 13,
        "p_n15_to_o16",  // 14,
        "he4_n15_to_f19",  // 15,
        "he4_o14_to_ne18",  // 16,
        "he4_o15_to_ne19",  // 17,
        "p_o16_to_f17",  // 18,
        "he4_o16_to_ne20",  // 19,
        "p_o17_to_f18",  // 20,
        "p_o18_to_f19",  // 21,
        "p_f17_to_ne18",  // 22,
        "p_f18_to_ne19",  // 23,
        "p_f19_to_ne20",  // 24,
        "he4_ne18_to_mg22",  // 25,
        "he4_ne20_to_mg24",  // 26,
        "c12_c12_to_he4_ne20",  // 27,
        "he4_n13_to_p_o16",  // 28,
        "p_n15_to_he4_c12",  // 29,
        "he4_o14_to_p_f17",  // 30,
        "c12_o16_to_he4_mg24",  // 31,
        "p_o17_to_he4_n14",  // 32,
        "p_o18_to_he4_n15",  // 33,
        "p_f18_to_he4_o15",  // 34,
        "p_f19_to_he4_o16",  // 35,
        "p_ne20_to_he4_f17",  // 36,
        "he4_he4_he4_to_c12"  // 37,
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
