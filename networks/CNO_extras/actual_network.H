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
        k_N13_to_C13_reaclib = 1,
        k_O14_to_N14_reaclib = 2,
        k_O15_to_N15_reaclib = 3,
        k_p_C12_to_N13_reaclib = 4,
        k_He4_C12_to_O16_reaclib = 5,
        k_p_C13_to_N14_reaclib = 6,
        k_p_N13_to_O14_reaclib = 7,
        k_p_N14_to_O15_reaclib = 8,
        k_He4_N14_to_F18_reaclib = 9,
        k_p_N15_to_O16_reaclib = 10,
        k_He4_N15_to_F19_reaclib = 11,
        k_He4_O14_to_Ne18_reaclib = 12,
        k_He4_O15_to_Ne19_reaclib = 13,
        k_p_O16_to_F17_reaclib = 14,
        k_He4_O16_to_Ne20_reaclib = 15,
        k_p_O17_to_F18_reaclib = 16,
        k_p_O18_to_F19_reaclib = 17,
        k_p_F17_to_Ne18_reaclib = 18,
        k_p_F18_to_Ne19_reaclib = 19,
        k_p_F19_to_Ne20_reaclib = 20,
        k_He4_Ne18_to_Mg22_reaclib = 21,
        k_He4_Ne20_to_Mg24_reaclib = 22,
        k_C12_C12_to_He4_Ne20_reaclib = 23,
        k_He4_N13_to_p_O16_reaclib = 24,
        k_p_N15_to_He4_C12_reaclib = 25,
        k_He4_O14_to_p_F17_reaclib = 26,
        k_C12_O16_to_He4_Mg24_reaclib = 27,
        k_p_O17_to_He4_N14_reaclib = 28,
        k_p_O18_to_He4_N15_reaclib = 29,
        k_p_F18_to_He4_O15_reaclib = 30,
        k_p_F19_to_He4_O16_reaclib = 31,
        k_p_Ne20_to_He4_F17_reaclib = 32,
        k_He4_He4_He4_to_C12_reaclib = 33,
        k_F17_to_O17_weaktab = 34,
        k_O17_to_F17_weaktab = 35,
        k_F18_to_Ne18_weaktab = 36,
        k_F18_to_O18_weaktab = 37,
        k_Ne18_to_F18_weaktab = 38,
        k_O18_to_F18_weaktab = 39,
        k_F19_to_Ne19_weaktab = 40,
        k_Ne19_to_F19_weaktab = 41,
        k_N13_to_p_C12_derived = 42,
        k_N14_to_p_C13_derived = 43,
        k_O14_to_p_N13_derived = 44,
        k_O15_to_p_N14_derived = 45,
        k_O16_to_p_N15_derived = 46,
        k_O16_to_He4_C12_derived = 47,
        k_F17_to_p_O16_derived = 48,
        k_F18_to_p_O17_derived = 49,
        k_F18_to_He4_N14_derived = 50,
        k_F19_to_p_O18_derived = 51,
        k_F19_to_He4_N15_derived = 52,
        k_Ne18_to_p_F17_derived = 53,
        k_Ne18_to_He4_O14_derived = 54,
        k_Ne19_to_p_F18_derived = 55,
        k_Ne19_to_He4_O15_derived = 56,
        k_Ne20_to_p_F19_derived = 57,
        k_Ne20_to_He4_O16_derived = 58,
        k_Mg22_to_He4_Ne18_derived = 59,
        k_Mg24_to_He4_Ne20_derived = 60,
        k_C12_to_He4_He4_He4_derived = 61,
        k_He4_C12_to_p_N15_derived = 62,
        k_He4_N14_to_p_O17_derived = 63,
        k_He4_N15_to_p_O18_derived = 64,
        k_He4_O15_to_p_F18_derived = 65,
        k_p_O16_to_He4_N13_derived = 66,
        k_He4_O16_to_p_F19_derived = 67,
        k_p_F17_to_He4_O14_derived = 68,
        k_He4_F17_to_p_Ne20_derived = 69,
        k_He4_Ne20_to_C12_C12_derived = 70,
        k_He4_Mg24_to_C12_O16_derived = 71,
        NumRates = k_He4_Mg24_to_C12_O16_derived
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "N13_to_C13_reaclib",  // 1,
        "O14_to_N14_reaclib",  // 2,
        "O15_to_N15_reaclib",  // 3,
        "p_C12_to_N13_reaclib",  // 4,
        "He4_C12_to_O16_reaclib",  // 5,
        "p_C13_to_N14_reaclib",  // 6,
        "p_N13_to_O14_reaclib",  // 7,
        "p_N14_to_O15_reaclib",  // 8,
        "He4_N14_to_F18_reaclib",  // 9,
        "p_N15_to_O16_reaclib",  // 10,
        "He4_N15_to_F19_reaclib",  // 11,
        "He4_O14_to_Ne18_reaclib",  // 12,
        "He4_O15_to_Ne19_reaclib",  // 13,
        "p_O16_to_F17_reaclib",  // 14,
        "He4_O16_to_Ne20_reaclib",  // 15,
        "p_O17_to_F18_reaclib",  // 16,
        "p_O18_to_F19_reaclib",  // 17,
        "p_F17_to_Ne18_reaclib",  // 18,
        "p_F18_to_Ne19_reaclib",  // 19,
        "p_F19_to_Ne20_reaclib",  // 20,
        "He4_Ne18_to_Mg22_reaclib",  // 21,
        "He4_Ne20_to_Mg24_reaclib",  // 22,
        "C12_C12_to_He4_Ne20_reaclib",  // 23,
        "He4_N13_to_p_O16_reaclib",  // 24,
        "p_N15_to_He4_C12_reaclib",  // 25,
        "He4_O14_to_p_F17_reaclib",  // 26,
        "C12_O16_to_He4_Mg24_reaclib",  // 27,
        "p_O17_to_He4_N14_reaclib",  // 28,
        "p_O18_to_He4_N15_reaclib",  // 29,
        "p_F18_to_He4_O15_reaclib",  // 30,
        "p_F19_to_He4_O16_reaclib",  // 31,
        "p_Ne20_to_He4_F17_reaclib",  // 32,
        "He4_He4_He4_to_C12_reaclib",  // 33,
        "F17_to_O17_weaktab",  // 34,
        "O17_to_F17_weaktab",  // 35,
        "F18_to_Ne18_weaktab",  // 36,
        "F18_to_O18_weaktab",  // 37,
        "Ne18_to_F18_weaktab",  // 38,
        "O18_to_F18_weaktab",  // 39,
        "F19_to_Ne19_weaktab",  // 40,
        "Ne19_to_F19_weaktab",  // 41,
        "N13_to_p_C12_derived",  // 42,
        "N14_to_p_C13_derived",  // 43,
        "O14_to_p_N13_derived",  // 44,
        "O15_to_p_N14_derived",  // 45,
        "O16_to_p_N15_derived",  // 46,
        "O16_to_He4_C12_derived",  // 47,
        "F17_to_p_O16_derived",  // 48,
        "F18_to_p_O17_derived",  // 49,
        "F18_to_He4_N14_derived",  // 50,
        "F19_to_p_O18_derived",  // 51,
        "F19_to_He4_N15_derived",  // 52,
        "Ne18_to_p_F17_derived",  // 53,
        "Ne18_to_He4_O14_derived",  // 54,
        "Ne19_to_p_F18_derived",  // 55,
        "Ne19_to_He4_O15_derived",  // 56,
        "Ne20_to_p_F19_derived",  // 57,
        "Ne20_to_He4_O16_derived",  // 58,
        "Mg22_to_He4_Ne18_derived",  // 59,
        "Mg24_to_He4_Ne20_derived",  // 60,
        "C12_to_He4_He4_He4_derived",  // 61,
        "He4_C12_to_p_N15_derived",  // 62,
        "He4_N14_to_p_O17_derived",  // 63,
        "He4_N15_to_p_O18_derived",  // 64,
        "He4_O15_to_p_F18_derived",  // 65,
        "p_O16_to_He4_N13_derived",  // 66,
        "He4_O16_to_p_F19_derived",  // 67,
        "p_F17_to_He4_O14_derived",  // 68,
        "He4_F17_to_p_Ne20_derived",  // 69,
        "He4_Ne20_to_C12_C12_derived",  // 70,
        "He4_Mg24_to_C12_O16_derived"  // 71,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    // A RatePair data table, 2D array of shape (NumNSERatePairs, 10).
    // Each row corresponds to a RatePair returned by _get_nse_rate_pairs().
    // Use enum names following NetworkSpecies and NetworkRates.

    //  - Columns 1–3: 1-based index of reactant nuclei (forward rate)
    //  - Columns 4–6: 1-based index of product nuclei (forward rate)
    //  - Columns 7-8: 1-based index of non (n, H1, He4) nuclei.
    //  - Column  9:   1-based index of the forward rate
    //  - Column  10:  1-based index of the reverse rate

    using namespace Species;
    using namespace Rates;

    constexpr int H1_index = 1;
    constexpr int N_index = -1;
    constexpr int He4_index = 2;
    constexpr int NumNSERatePairs = 28;

    inline AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, NumNSERatePairs, 1, 10, amrex::Order::C> rate_pair_data {
        H1, C12, -1, N13, -1, -1, C12, N13, k_p_C12_to_N13_reaclib, k_N13_to_p_C12_derived,
        He4, C12, -1, O16, -1, -1, C12, O16, k_He4_C12_to_O16_reaclib, k_O16_to_He4_C12_derived,
        H1, C13, -1, N14, -1, -1, C13, N14, k_p_C13_to_N14_reaclib, k_N14_to_p_C13_derived,
        H1, N13, -1, O14, -1, -1, N13, O14, k_p_N13_to_O14_reaclib, k_O14_to_p_N13_derived,
        H1, N14, -1, O15, -1, -1, N14, O15, k_p_N14_to_O15_reaclib, k_O15_to_p_N14_derived,
        He4, N14, -1, F18, -1, -1, N14, F18, k_He4_N14_to_F18_reaclib, k_F18_to_He4_N14_derived,
        H1, N15, -1, O16, -1, -1, N15, O16, k_p_N15_to_O16_reaclib, k_O16_to_p_N15_derived,
        He4, N15, -1, F19, -1, -1, N15, F19, k_He4_N15_to_F19_reaclib, k_F19_to_He4_N15_derived,
        He4, O14, -1, Ne18, -1, -1, O14, Ne18, k_He4_O14_to_Ne18_reaclib, k_Ne18_to_He4_O14_derived,
        He4, O15, -1, Ne19, -1, -1, O15, Ne19, k_He4_O15_to_Ne19_reaclib, k_Ne19_to_He4_O15_derived,
        H1, O16, -1, F17, -1, -1, O16, F17, k_p_O16_to_F17_reaclib, k_F17_to_p_O16_derived,
        He4, O16, -1, Ne20, -1, -1, O16, Ne20, k_He4_O16_to_Ne20_reaclib, k_Ne20_to_He4_O16_derived,
        H1, O17, -1, F18, -1, -1, O17, F18, k_p_O17_to_F18_reaclib, k_F18_to_p_O17_derived,
        H1, O18, -1, F19, -1, -1, O18, F19, k_p_O18_to_F19_reaclib, k_F19_to_p_O18_derived,
        H1, F17, -1, Ne18, -1, -1, F17, Ne18, k_p_F17_to_Ne18_reaclib, k_Ne18_to_p_F17_derived,
        H1, F18, -1, Ne19, -1, -1, F18, Ne19, k_p_F18_to_Ne19_reaclib, k_Ne19_to_p_F18_derived,
        H1, F19, -1, Ne20, -1, -1, F19, Ne20, k_p_F19_to_Ne20_reaclib, k_Ne20_to_p_F19_derived,
        He4, Ne18, -1, Mg22, -1, -1, Ne18, Mg22, k_He4_Ne18_to_Mg22_reaclib, k_Mg22_to_He4_Ne18_derived,
        He4, Ne20, -1, Mg24, -1, -1, Ne20, Mg24, k_He4_Ne20_to_Mg24_reaclib, k_Mg24_to_He4_Ne20_derived,
        He4, N13, -1, H1, O16, -1, N13, O16, k_He4_N13_to_p_O16_reaclib, k_p_O16_to_He4_N13_derived,
        H1, N15, -1, He4, C12, -1, N15, C12, k_p_N15_to_He4_C12_reaclib, k_He4_C12_to_p_N15_derived,
        He4, O14, -1, H1, F17, -1, O14, F17, k_He4_O14_to_p_F17_reaclib, k_p_F17_to_He4_O14_derived,
        H1, O17, -1, He4, N14, -1, O17, N14, k_p_O17_to_He4_N14_reaclib, k_He4_N14_to_p_O17_derived,
        H1, O18, -1, He4, N15, -1, O18, N15, k_p_O18_to_He4_N15_reaclib, k_He4_N15_to_p_O18_derived,
        H1, F18, -1, He4, O15, -1, F18, O15, k_p_F18_to_He4_O15_reaclib, k_He4_O15_to_p_F18_derived,
        H1, F19, -1, He4, O16, -1, F19, O16, k_p_F19_to_He4_O16_reaclib, k_He4_O16_to_p_F19_derived,
        He4, He4, He4, C12, -1, -1, C12, -1, k_He4_He4_He4_to_C12_reaclib, k_C12_to_He4_He4_He4_derived,
        He4, F17, -1, H1, Ne20, -1, F17, Ne20, k_He4_F17_to_p_Ne20_derived, k_p_Ne20_to_He4_F17_reaclib
    };
}
#endif

#endif
