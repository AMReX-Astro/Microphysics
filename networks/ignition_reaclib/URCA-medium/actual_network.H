#ifndef actual_network_H
#define actual_network_H

#include <string_view>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Loop.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex;

void actual_network_init();

constexpr std::string_view network_name = "pynucastro-cxx";

namespace network
{

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr amrex::Real bion () {
        using namespace Species;

        static_assert(spec >= 1 && spec <= NumSpec);

        // Set the binding energy of the element

        if constexpr (spec == N) {
            return 0.0_rt;
        }
        else if constexpr (spec == H1) {
            return 0.0_rt;
        }
        else if constexpr (spec == He4) {
            return 28.295662457999697_rt;
        }
        else if constexpr (spec == C12) {
            return 92.16173498399803_rt;
        }
        else if constexpr (spec == N13) {
            return 94.10522604799917_rt;
        }
        else if constexpr (spec == O16) {
            return 127.6193154119992_rt;
        }
        else if constexpr (spec == Ne20) {
            return 160.64482384000075_rt;
        }
        else if constexpr (spec == Ne23) {
            return 182.97089593999772_rt;
        }
        else if constexpr (spec == Na23) {
            return 186.56435240400242_rt;
        }
        else if constexpr (spec == Mg23) {
            return 181.7258218679999_rt;
        }
        else if constexpr (spec == Mg24) {
            return 198.2570479679962_rt;
        }


        // Return zero if we don't recognize the species.
        return 0.0_rt;
    }

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr amrex::Real mion () {
        using namespace Species;

        static_assert(spec >= 1 && spec <= NumSpec);

        // Set the mass of the element

        if constexpr (spec == N) {
            return 1.674927498034172e-24_rt;
        }
        else if constexpr (spec == H1) {
            return 1.6735328377636005e-24_rt;
        }
        else if constexpr (spec == He4) {
            return 6.646479071584587e-24_rt;
        }
        else if constexpr (spec == C12) {
            return 1.99264687992e-23_rt;
        }
        else if constexpr (spec == N13) {
            return 2.1596537049448796e-23_rt;
        }
        else if constexpr (spec == O16) {
            return 2.6560180592333686e-23_rt;
        }
        else if constexpr (spec == Ne20) {
            return 3.3198227947612416e-23_rt;
        }
        else if constexpr (spec == Ne23) {
            return 3.818321060312286e-23_rt;
        }
        else if constexpr (spec == Na23) {
            return 3.817541002484691e-23_rt;
        }
        else if constexpr (spec == Mg23) {
            return 3.8182640828719474e-23_rt;
        }
        else if constexpr (spec == Mg24) {
            return 3.9828098739467446e-23_rt;
        }


        // Return zero if we don't recognize the species.
        return 0.0_rt;
    }

    // Legacy (non-templated) interfaces

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real bion (int spec) {
        using namespace Species;

        amrex::Real b = 0.0_rt;

        // Set the binding energy of the element
        constexpr_for<1, NumSpec+1>([&] (auto n) {
            if (n == spec) {
                b = bion<n>();
            }
        });

        return b;
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real mion (int spec) {
        using namespace Species;

        amrex::Real m = 0.0_rt;

        constexpr_for<1, NumSpec+1>([&] (auto n) {
            if (n == spec) {
                m = mion<n>();
            }
        });

        return m;
    }
}

namespace Rates
{

    enum NetworkRates :
    std::uint8_t
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
