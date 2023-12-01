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
        k_Ne20_to_He4_O16 = 1,
        k_He4_O16_to_Ne20 = 2,
        k_He4_Ne20_to_Mg24 = 3,
        k_He4_Mg24_to_Si28 = 4,
        k_p_Al27_to_Si28 = 5,
        k_He4_Al27_to_P31 = 6,
        k_He4_Si28_to_S32 = 7,
        k_p_P31_to_S32 = 8,
        k_O16_O16_to_p_P31 = 9,
        k_O16_O16_to_He4_Si28 = 10,
        k_He4_Mg24_to_p_Al27 = 11,
        k_p_Al27_to_He4_Mg24 = 12,
        k_He4_Si28_to_p_P31 = 13,
        k_p_P31_to_He4_Si28 = 14,
        k_F20_to_O20 = 15,
        k_Ne20_to_F20 = 16,
        k_O20_to_F20 = 17,
        k_F20_to_Ne20 = 18,
        NumRates = k_F20_to_Ne20
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
        "Ne20_to_He4_O16",  // 1,
        "He4_O16_to_Ne20",  // 2,
        "He4_Ne20_to_Mg24",  // 3,
        "He4_Mg24_to_Si28",  // 4,
        "p_Al27_to_Si28",  // 5,
        "He4_Al27_to_P31",  // 6,
        "He4_Si28_to_S32",  // 7,
        "p_P31_to_S32",  // 8,
        "O16_O16_to_p_P31",  // 9,
        "O16_O16_to_He4_Si28",  // 10,
        "He4_Mg24_to_p_Al27",  // 11,
        "p_Al27_to_He4_Mg24",  // 12,
        "He4_Si28_to_p_P31",  // 13,
        "p_P31_to_He4_Si28",  // 14,
        "F20_to_O20",  // 15,
        "Ne20_to_F20",  // 16,
        "O20_to_F20",  // 17,
        "F20_to_Ne20"  // 18,
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
