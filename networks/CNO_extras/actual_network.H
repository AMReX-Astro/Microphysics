#ifndef actual_network_H
#define actual_network_H

#include <string_view>

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Loop.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex::literals;

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

        if constexpr (spec == H1) {
            return 0.0_rt;
        }
        else if constexpr (spec == He4) {
            return 28.295662457999697_rt;
        }
        else if constexpr (spec == C12) {
            return 92.16173498399803_rt;
        }
        else if constexpr (spec == C13) {
            return 97.10804378399916_rt;
        }
        else if constexpr (spec == N13) {
            return 94.10522604799917_rt;
        }
        else if constexpr (spec == N14) {
            return 104.65860734799753_rt;
        }
        else if constexpr (spec == N15) {
            return 115.49190414799887_rt;
        }
        else if constexpr (spec == O14) {
            return 98.73189611199996_rt;
        }
        else if constexpr (spec == O15) {
            return 111.95539521199862_rt;
        }
        else if constexpr (spec == O16) {
            return 127.6193154119992_rt;
        }
        else if constexpr (spec == O17) {
            return 131.76239561199873_rt;
        }
        else if constexpr (spec == O18) {
            return 139.8077658120019_rt;
        }
        else if constexpr (spec == F17) {
            return 128.21958437599824_rt;
        }
        else if constexpr (spec == F18) {
            return 137.36950247599816_rt;
        }
        else if constexpr (spec == F19) {
            return 147.80136567599766_rt;
        }
        else if constexpr (spec == Ne18) {
            return 132.14265544000227_rt;
        }
        else if constexpr (spec == Ne19) {
            return 143.7795235400008_rt;
        }
        else if constexpr (spec == Ne20) {
            return 160.64482384000075_rt;
        }
        else if constexpr (spec == Mg22) {
            return 168.58082376800303_rt;
        }
        else if constexpr (spec == Mg24) {
            return 198.2570479679962_rt;
        }
        else if constexpr (spec == Fe56) {
            return 492.2599506639962_rt;
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

        if constexpr (spec == H1) {
            return 1.6735328377636005e-24_rt;
        }
        else if constexpr (spec == He4) {
            return 6.646479071584587e-24_rt;
        }
        else if constexpr (spec == C12) {
            return 1.99264687992e-23_rt;
        }
        else if constexpr (spec == C13) {
            return 2.1592578700883826e-23_rt;
        }
        else if constexpr (spec == N13) {
            return 2.1596537049448796e-23_rt;
        }
        else if constexpr (spec == N14) {
            return 2.3252651436495096e-23_rt;
        }
        else if constexpr (spec == N15) {
            return 2.490826682883827e-23_rt;
        }
        else if constexpr (spec == O14) {
            return 2.3261822098665393e-23_rt;
        }
        else if constexpr (spec == O15) {
            return 2.4913176568383383e-23_rt;
        }
        else if constexpr (spec == O16) {
            return 2.6560180592333686e-23_rt;
        }
        else if constexpr (spec == O17) {
            return 2.822772237905709e-23_rt;
        }
        else if constexpr (spec == O18) {
            return 2.9888307701990364e-23_rt;
        }
        else if constexpr (spec == F17) {
            return 2.8232643353472436e-23_rt;
        }
        else if constexpr (spec == F18) {
            return 2.989125964092377e-23_rt;
        }
        else if constexpr (spec == F19) {
            return 3.1547590653659727e-23_rt;
        }
        else if constexpr (spec == Ne18) {
            return 2.9899182681834417e-23_rt;
        }
        else if constexpr (spec == Ne19) {
            return 3.155336557821978e-23_rt;
        }
        else if constexpr (spec == Ne20) {
            return 3.3198227947612416e-23_rt;
        }
        else if constexpr (spec == Mg22) {
            return 3.653114641825797e-23_rt;
        }
        else if constexpr (spec == Mg24) {
            return 3.9828098739467446e-23_rt;
        }
        else if constexpr (spec == Fe56) {
            return 9.288214565329033e-23_rt;
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
        amrex::constexpr_for<1, NumSpec+1>([&] (auto n) {
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

        amrex::constexpr_for<1, NumSpec+1>([&] (auto n) {
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
        k_N13_to_C13_weak_wc12 = 1,
        k_O14_to_N14_weak_wc12 = 2,
        k_O15_to_N15_weak_wc12 = 3,
        k_F17_to_O17_weak_wc12 = 4,
        k_F18_to_O18_weak_wc12 = 5,
        k_Ne18_to_F18_weak_wc12 = 6,
        k_Ne19_to_F19_weak_wc12 = 7,
        k_N13_to_p_C12 = 8,
        k_N14_to_p_C13 = 9,
        k_O14_to_p_N13 = 10,
        k_O15_to_p_N14 = 11,
        k_O16_to_p_N15 = 12,
        k_O16_to_He4_C12 = 13,
        k_F17_to_p_O16 = 14,
        k_F18_to_p_O17 = 15,
        k_F18_to_He4_N14 = 16,
        k_F19_to_p_O18 = 17,
        k_F19_to_He4_N15 = 18,
        k_Ne18_to_p_F17 = 19,
        k_Ne18_to_He4_O14 = 20,
        k_Ne19_to_p_F18 = 21,
        k_Ne19_to_He4_O15 = 22,
        k_Ne20_to_p_F19 = 23,
        k_Ne20_to_He4_O16 = 24,
        k_Mg22_to_He4_Ne18 = 25,
        k_Mg24_to_He4_Ne20 = 26,
        k_C12_to_He4_He4_He4 = 27,
        k_p_C12_to_N13 = 28,
        k_He4_C12_to_O16 = 29,
        k_p_C13_to_N14 = 30,
        k_p_N13_to_O14 = 31,
        k_p_N14_to_O15 = 32,
        k_He4_N14_to_F18 = 33,
        k_p_N15_to_O16 = 34,
        k_He4_N15_to_F19 = 35,
        k_He4_O14_to_Ne18 = 36,
        k_He4_O15_to_Ne19 = 37,
        k_p_O16_to_F17 = 38,
        k_He4_O16_to_Ne20 = 39,
        k_p_O17_to_F18 = 40,
        k_p_O18_to_F19 = 41,
        k_p_F17_to_Ne18 = 42,
        k_p_F18_to_Ne19 = 43,
        k_p_F19_to_Ne20 = 44,
        k_He4_Ne18_to_Mg22 = 45,
        k_He4_Ne20_to_Mg24 = 46,
        k_He4_C12_to_p_N15 = 47,
        k_C12_C12_to_He4_Ne20 = 48,
        k_He4_N13_to_p_O16 = 49,
        k_He4_N14_to_p_O17 = 50,
        k_p_N15_to_He4_C12 = 51,
        k_He4_N15_to_p_O18 = 52,
        k_He4_O14_to_p_F17 = 53,
        k_He4_O15_to_p_F18 = 54,
        k_p_O16_to_He4_N13 = 55,
        k_He4_O16_to_p_F19 = 56,
        k_C12_O16_to_He4_Mg24 = 57,
        k_p_O17_to_He4_N14 = 58,
        k_p_O18_to_He4_N15 = 59,
        k_p_F17_to_He4_O14 = 60,
        k_He4_F17_to_p_Ne20 = 61,
        k_p_F18_to_He4_O15 = 62,
        k_p_F19_to_He4_O16 = 63,
        k_p_Ne20_to_He4_F17 = 64,
        k_He4_Ne20_to_C12_C12 = 65,
        k_He4_Mg24_to_C12_O16 = 66,
        k_He4_He4_He4_to_C12 = 67,
        NumRates = k_He4_He4_He4_to_C12
    };

    // number of reaclib rates

    const int NrateReaclib = 67;

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
        "N13_to_p_C12",  // 8,
        "N14_to_p_C13",  // 9,
        "O14_to_p_N13",  // 10,
        "O15_to_p_N14",  // 11,
        "O16_to_p_N15",  // 12,
        "O16_to_He4_C12",  // 13,
        "F17_to_p_O16",  // 14,
        "F18_to_p_O17",  // 15,
        "F18_to_He4_N14",  // 16,
        "F19_to_p_O18",  // 17,
        "F19_to_He4_N15",  // 18,
        "Ne18_to_p_F17",  // 19,
        "Ne18_to_He4_O14",  // 20,
        "Ne19_to_p_F18",  // 21,
        "Ne19_to_He4_O15",  // 22,
        "Ne20_to_p_F19",  // 23,
        "Ne20_to_He4_O16",  // 24,
        "Mg22_to_He4_Ne18",  // 25,
        "Mg24_to_He4_Ne20",  // 26,
        "C12_to_He4_He4_He4",  // 27,
        "p_C12_to_N13",  // 28,
        "He4_C12_to_O16",  // 29,
        "p_C13_to_N14",  // 30,
        "p_N13_to_O14",  // 31,
        "p_N14_to_O15",  // 32,
        "He4_N14_to_F18",  // 33,
        "p_N15_to_O16",  // 34,
        "He4_N15_to_F19",  // 35,
        "He4_O14_to_Ne18",  // 36,
        "He4_O15_to_Ne19",  // 37,
        "p_O16_to_F17",  // 38,
        "He4_O16_to_Ne20",  // 39,
        "p_O17_to_F18",  // 40,
        "p_O18_to_F19",  // 41,
        "p_F17_to_Ne18",  // 42,
        "p_F18_to_Ne19",  // 43,
        "p_F19_to_Ne20",  // 44,
        "He4_Ne18_to_Mg22",  // 45,
        "He4_Ne20_to_Mg24",  // 46,
        "He4_C12_to_p_N15",  // 47,
        "C12_C12_to_He4_Ne20",  // 48,
        "He4_N13_to_p_O16",  // 49,
        "He4_N14_to_p_O17",  // 50,
        "p_N15_to_He4_C12",  // 51,
        "He4_N15_to_p_O18",  // 52,
        "He4_O14_to_p_F17",  // 53,
        "He4_O15_to_p_F18",  // 54,
        "p_O16_to_He4_N13",  // 55,
        "He4_O16_to_p_F19",  // 56,
        "C12_O16_to_He4_Mg24",  // 57,
        "p_O17_to_He4_N14",  // 58,
        "p_O18_to_He4_N15",  // 59,
        "p_F17_to_He4_O14",  // 60,
        "He4_F17_to_p_Ne20",  // 61,
        "p_F18_to_He4_O15",  // 62,
        "p_F19_to_He4_O16",  // 63,
        "p_Ne20_to_He4_F17",  // 64,
        "He4_Ne20_to_C12_C12",  // 65,
        "He4_Mg24_to_C12_O16",  // 66,
        "He4_He4_He4_to_C12"  // 67,
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

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, amrex::Order::C> rate_indices;
}
#endif

#endif
