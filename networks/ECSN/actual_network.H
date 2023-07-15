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
        k_ne20_to_he4_o16 = 1,
        k_he4_o16_to_ne20 = 2,
        k_he4_ne20_to_mg24 = 3,
        k_he4_mg24_to_si28 = 4,
        k_p_al27_to_si28 = 5,
        k_he4_al27_to_p31 = 6,
        k_he4_si28_to_s32 = 7,
        k_p_p31_to_s32 = 8,
        k_o16_o16_to_p_p31 = 9,
        k_o16_o16_to_he4_si28 = 10,
        k_he4_mg24_to_p_al27 = 11,
        k_p_al27_to_he4_mg24 = 12,
        k_he4_si28_to_p_p31 = 13,
        k_p_p31_to_he4_si28 = 14,
        k_f20_to_o20 = 15,
        k_ne20_to_f20 = 16,
        k_o20_to_f20 = 17,
        k_f20_to_ne20 = 18,
        NumRates = k_f20_to_ne20
    };

    // number of reaclib rates

    const int NrateReaclib = 14;

    // number of tabular rates

    const int NrateTabular = 4;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "ne20_to_he4_o16",  // 1,
        "he4_o16_to_ne20",  // 2,
        "he4_ne20_to_mg24",  // 3,
        "he4_mg24_to_si28",  // 4,
        "p_al27_to_si28",  // 5,
        "he4_al27_to_p31",  // 6,
        "he4_si28_to_s32",  // 7,
        "p_p31_to_s32",  // 8,
        "o16_o16_to_p_p31",  // 9,
        "o16_o16_to_he4_si28",  // 10,
        "he4_mg24_to_p_al27",  // 11,
        "p_al27_to_he4_mg24",  // 12,
        "he4_si28_to_p_p31",  // 13,
        "p_p31_to_he4_si28",  // 14,
        "f20_to_o20",  // 15,
        "ne20_to_f20",  // 16,
        "o20_to_f20",  // 17,
        "f20_to_ne20"  // 18,
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
    extern AMREX_GPU_MANAGED bool initialized;
}
#endif

#endif
