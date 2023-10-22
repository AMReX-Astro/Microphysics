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
        k_he4_o16_to_ne20 = 3,
        k_he4_ne20_to_mg24 = 4,
        k_p_na23_to_mg24 = 5,
        k_n_mg23_to_mg24 = 6,
        k_c12_c12_to_n_mg23 = 7,
        k_c12_c12_to_p_na23 = 8,
        k_c12_c12_to_he4_ne20 = 9,
        k_he4_n13_to_p_o16 = 10,
        k_c12_o16_to_he4_mg24 = 11,
        k_he4_ne20_to_p_na23 = 12,
        k_p_na23_to_he4_ne20 = 13,
        k_n_mg23_to_p_na23 = 14,
        k_n_mg23_to_he4_ne20 = 15,
        k_n_mg23_to_c12_c12 = 16,
        k_he4_he4_he4_to_c12 = 17,
        k_na23_to_ne23 = 18,
        k_ne23_to_na23 = 19,
        k_mg23_to_na23 = 20,
        k_n_to_p = 21,
        k_p_to_n = 22,
        NumRates = k_p_to_n
    };

    // number of reaclib rates

    const int NrateReaclib = 17;

    // number of tabular rates

    const int NrateTabular = 5;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "p_c12_to_n13",  // 1,
        "he4_c12_to_o16",  // 2,
        "he4_o16_to_ne20",  // 3,
        "he4_ne20_to_mg24",  // 4,
        "p_na23_to_mg24",  // 5,
        "n_mg23_to_mg24",  // 6,
        "c12_c12_to_n_mg23",  // 7,
        "c12_c12_to_p_na23",  // 8,
        "c12_c12_to_he4_ne20",  // 9,
        "he4_n13_to_p_o16",  // 10,
        "c12_o16_to_he4_mg24",  // 11,
        "he4_ne20_to_p_na23",  // 12,
        "p_na23_to_he4_ne20",  // 13,
        "n_mg23_to_p_na23",  // 14,
        "n_mg23_to_he4_ne20",  // 15,
        "n_mg23_to_c12_c12",  // 16,
        "he4_he4_he4_to_c12",  // 17,
        "na23_to_ne23",  // 18,
        "ne23_to_na23",  // 19,
        "mg23_to_na23",  // 20,
        "n_to_p",  // 21,
        "p_to_n"  // 22,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int h1_index = 1;
    constexpr int n_index = 0;
    constexpr int he4_index = 2;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
