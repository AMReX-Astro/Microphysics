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
        k_p_C12_to_N13 = 1,
        k_He4_C12_to_O16 = 2,
        k_He4_O16_to_Ne20 = 3,
        k_He4_Ne20_to_Mg24 = 4,
        k_p_Na23_to_Mg24 = 5,
        k_n_Mg23_to_Mg24 = 6,
        k_C12_C12_to_n_Mg23 = 7,
        k_C12_C12_to_p_Na23 = 8,
        k_C12_C12_to_He4_Ne20 = 9,
        k_He4_N13_to_p_O16 = 10,
        k_C12_O16_to_He4_Mg24 = 11,
        k_He4_Ne20_to_p_Na23 = 12,
        k_p_Na23_to_He4_Ne20 = 13,
        k_n_Mg23_to_p_Na23 = 14,
        k_n_Mg23_to_He4_Ne20 = 15,
        k_n_Mg23_to_C12_C12 = 16,
        k_He4_He4_He4_to_C12 = 17,
        k_Na23_to_Ne23 = 18,
        k_Ne23_to_Na23 = 19,
        k_Mg23_to_Na23 = 20,
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
        "p_C12_to_N13",  // 1,
        "He4_C12_to_O16",  // 2,
        "He4_O16_to_Ne20",  // 3,
        "He4_Ne20_to_Mg24",  // 4,
        "p_Na23_to_Mg24",  // 5,
        "n_Mg23_to_Mg24",  // 6,
        "C12_C12_to_n_Mg23",  // 7,
        "C12_C12_to_p_Na23",  // 8,
        "C12_C12_to_He4_Ne20",  // 9,
        "He4_N13_to_p_O16",  // 10,
        "C12_O16_to_He4_Mg24",  // 11,
        "He4_Ne20_to_p_Na23",  // 12,
        "p_Na23_to_He4_Ne20",  // 13,
        "n_Mg23_to_p_Na23",  // 14,
        "n_Mg23_to_He4_Ne20",  // 15,
        "n_Mg23_to_C12_C12",  // 16,
        "He4_He4_He4_to_C12",  // 17,
        "Na23_to_Ne23",  // 18,
        "Ne23_to_Na23",  // 19,
        "Mg23_to_Na23",  // 20,
        "n_to_p",  // 21,
        "p_to_n"  // 22,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int H1_index = 1;
    constexpr int N_index = 0;
    constexpr int He4_index = 2;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
