#ifndef ACTUAL_NETWORK_H
#define ACTUAL_NETWORK_H

#include <random>
#include <string_view>

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Loop.H>

#include <extern_parameters.H>
#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex::literals;


inline
void actual_network_init() {


}

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
        else if constexpr (spec == He3) {
            return 7.718041347999588_rt;
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
        else if constexpr (spec == He3) {
            return 5.008234515140786e-24_rt;
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
        k_N13_to_C13_beta_pos_reaclib = 1,
        k_O14_to_N14_beta_pos_reaclib = 2,
        k_O15_to_N15_beta_pos_reaclib = 3,
        k_C12_p_to_N13_reaclib = 4,
        k_C12_He4_to_O16_reaclib = 5,
        k_C13_p_to_N14_reaclib = 6,
        k_N13_p_to_O14_reaclib = 7,
        k_N14_p_to_O15_reaclib = 8,
        k_N14_He4_to_F18_reaclib = 9,
        k_N15_p_to_O16_reaclib = 10,
        k_N15_He4_to_F19_reaclib = 11,
        k_O14_He4_to_Ne18_reaclib = 12,
        k_O15_He4_to_Ne19_reaclib = 13,
        k_O16_p_to_F17_reaclib = 14,
        k_O16_He4_to_Ne20_reaclib = 15,
        k_O17_p_to_F18_reaclib = 16,
        k_O18_p_to_F19_reaclib = 17,
        k_F17_p_to_Ne18_reaclib = 18,
        k_F18_p_to_Ne19_reaclib = 19,
        k_F19_p_to_Ne20_reaclib = 20,
        k_Ne20_He4_to_Mg24_reaclib = 21,
        k_N13_He4_to_p_O16_reaclib = 22,
        k_N15_p_to_He4_C12_reaclib = 23,
        k_O14_He4_to_p_F17_reaclib = 24,
        k_O17_p_to_He4_N14_reaclib = 25,
        k_O18_p_to_He4_N15_reaclib = 26,
        k_F18_p_to_He4_O15_reaclib = 27,
        k_F19_p_to_He4_O16_reaclib = 28,
        k_Ne20_p_to_He4_F17_reaclib = 29,
        k_He4_He4_He4_to_C12_reaclib = 30,
        k_Ne18_He4_to_p_Na21_reaclib = 31,
        k_Ne18_He4_to_Mg22_reaclib = 32,
        k_Na21_p_to_Mg22_reaclib = 33,
        k_Ne19_p_to_Na20_reaclib = 34,
        k_Ne20_p_to_Na21_reaclib = 35,
        k_p_p_to_d_beta_pos_reaclib = 36,
        k_p_p_to_d_electron_capture_reaclib = 37,
        k_He3_He3_to_p_p_He4_reaclib = 38,
        k_He3_p_to_He4_beta_pos_reaclib = 39,
        k_He3_He4_to_Be7_reaclib = 40,
        k_F17_to_O17_electron_capture_weaktab = 41,
        k_O17_to_F17_beta_neg_weaktab = 42,
        k_F18_to_Ne18_beta_neg_weaktab = 43,
        k_F18_to_O18_electron_capture_weaktab = 44,
        k_Ne18_to_F18_electron_capture_weaktab = 45,
        k_O18_to_F18_beta_neg_weaktab = 46,
        k_F19_to_Ne19_beta_neg_weaktab = 47,
        k_Ne19_to_F19_electron_capture_weaktab = 48,
        k_Ne18_He4_to_Mg22_approx = 49,
        k_Ne19_3p_to_Mg22_beta_pos_modified = 50,
        k_Ne20_p_p_to_Mg22_modified = 51,
        k_p_p_p_to_He3_beta_pos_modified = 52,
        k_p_p_p_to_He3_electron_capture_modified = 53,
        k_He3_He4_p_to_He4_He4_beta_pos_modified = 54,
        k_N13_to_p_C12_derived = 55,
        k_N14_to_p_C13_derived = 56,
        k_O14_to_p_N13_derived = 57,
        k_O15_to_p_N14_derived = 58,
        k_O16_to_p_N15_derived = 59,
        k_O16_to_He4_C12_derived = 60,
        k_F17_to_p_O16_derived = 61,
        k_F18_to_p_O17_derived = 62,
        k_F18_to_He4_N14_derived = 63,
        k_F19_to_p_O18_derived = 64,
        k_F19_to_He4_N15_derived = 65,
        k_Ne18_to_p_F17_derived = 66,
        k_Ne18_to_He4_O14_derived = 67,
        k_Ne19_to_p_F18_derived = 68,
        k_Ne19_to_He4_O15_derived = 69,
        k_Ne20_to_p_F19_derived = 70,
        k_Ne20_to_He4_O16_derived = 71,
        k_Mg22_to_He4_Ne18_derived = 72,
        k_Mg24_to_He4_Ne20_derived = 73,
        k_C12_to_3He4_derived = 74,
        k_C12_He4_to_p_N15_derived = 75,
        k_N14_He4_to_p_O17_derived = 76,
        k_N15_He4_to_p_O18_derived = 77,
        k_O15_He4_to_p_F18_derived = 78,
        k_O16_p_to_He4_N13_derived = 79,
        k_O16_He4_to_p_F19_derived = 80,
        k_F17_p_to_He4_O14_derived = 81,
        k_F17_He4_to_p_Ne20_derived = 82,
        k_Na21_p_to_He4_Ne18_derived = 83,
        k_Mg22_to_p_Na21_derived = 84,
        NumRates = k_Mg22_to_p_Na21_derived
    };

    enum ScreenPairs :
    std::uint8_t
    {
        k_He4_Ne18 = 1,
        k_p_F18 = 2,
        k_p_C12 = 3,
        k_p_Na21 = 4,
        k_He4_O16 = 5,
        k_He3_He3 = 6,
        k_He4_F17 = 7,
        k_He4_N14 = 8,
        k_He3_He4 = 9,
        k_p_F19 = 10,
        k_p_C13 = 11,
        k_He4_C12 = 12,
        k_p_p = 13,
        k_p_O17 = 14,
        k_p_N15 = 15,
        k_p_Ne20 = 16,
        k_p_O18 = 17,
        k_p_He3 = 18,
        k_He4_N15 = 19,
        k_He4_Ne20 = 20,
        k_p_N13 = 21,
        k_He4_O14 = 22,
        k_He4_He4 = 23,
        k_He4_Be8 = 24,
        k_p_O16 = 25,
        k_p_N14 = 26,
        k_p_F17 = 27,
        k_p_Ne19 = 28,
        k_He4_O15 = 29,
        k_He4_N13 = 30,
        NumScreenPairs = k_He4_N13
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "N13_to_C13_beta_pos_reaclib",  // 1,
        "O14_to_N14_beta_pos_reaclib",  // 2,
        "O15_to_N15_beta_pos_reaclib",  // 3,
        "C12_p_to_N13_reaclib",  // 4,
        "C12_He4_to_O16_reaclib",  // 5,
        "C13_p_to_N14_reaclib",  // 6,
        "N13_p_to_O14_reaclib",  // 7,
        "N14_p_to_O15_reaclib",  // 8,
        "N14_He4_to_F18_reaclib",  // 9,
        "N15_p_to_O16_reaclib",  // 10,
        "N15_He4_to_F19_reaclib",  // 11,
        "O14_He4_to_Ne18_reaclib",  // 12,
        "O15_He4_to_Ne19_reaclib",  // 13,
        "O16_p_to_F17_reaclib",  // 14,
        "O16_He4_to_Ne20_reaclib",  // 15,
        "O17_p_to_F18_reaclib",  // 16,
        "O18_p_to_F19_reaclib",  // 17,
        "F17_p_to_Ne18_reaclib",  // 18,
        "F18_p_to_Ne19_reaclib",  // 19,
        "F19_p_to_Ne20_reaclib",  // 20,
        "Ne20_He4_to_Mg24_reaclib",  // 21,
        "N13_He4_to_p_O16_reaclib",  // 22,
        "N15_p_to_He4_C12_reaclib",  // 23,
        "O14_He4_to_p_F17_reaclib",  // 24,
        "O17_p_to_He4_N14_reaclib",  // 25,
        "O18_p_to_He4_N15_reaclib",  // 26,
        "F18_p_to_He4_O15_reaclib",  // 27,
        "F19_p_to_He4_O16_reaclib",  // 28,
        "Ne20_p_to_He4_F17_reaclib",  // 29,
        "He4_He4_He4_to_C12_reaclib",  // 30,
        "Ne18_He4_to_p_Na21_reaclib",  // 31,
        "Ne18_He4_to_Mg22_reaclib",  // 32,
        "Na21_p_to_Mg22_reaclib",  // 33,
        "Ne19_p_to_Na20_reaclib",  // 34,
        "Ne20_p_to_Na21_reaclib",  // 35,
        "p_p_to_d_beta_pos_reaclib",  // 36,
        "p_p_to_d_electron_capture_reaclib",  // 37,
        "He3_He3_to_p_p_He4_reaclib",  // 38,
        "He3_p_to_He4_beta_pos_reaclib",  // 39,
        "He3_He4_to_Be7_reaclib",  // 40,
        "F17_to_O17_electron_capture_weaktab",  // 41,
        "O17_to_F17_beta_neg_weaktab",  // 42,
        "F18_to_Ne18_beta_neg_weaktab",  // 43,
        "F18_to_O18_electron_capture_weaktab",  // 44,
        "Ne18_to_F18_electron_capture_weaktab",  // 45,
        "O18_to_F18_beta_neg_weaktab",  // 46,
        "F19_to_Ne19_beta_neg_weaktab",  // 47,
        "Ne19_to_F19_electron_capture_weaktab",  // 48,
        "Ne18_He4_to_Mg22_approx",  // 49,
        "Ne19_3p_to_Mg22_beta_pos_modified",  // 50,
        "Ne20_p_p_to_Mg22_modified",  // 51,
        "p_p_p_to_He3_beta_pos_modified",  // 52,
        "p_p_p_to_He3_electron_capture_modified",  // 53,
        "He3_He4_p_to_He4_He4_beta_pos_modified",  // 54,
        "N13_to_p_C12_derived",  // 55,
        "N14_to_p_C13_derived",  // 56,
        "O14_to_p_N13_derived",  // 57,
        "O15_to_p_N14_derived",  // 58,
        "O16_to_p_N15_derived",  // 59,
        "O16_to_He4_C12_derived",  // 60,
        "F17_to_p_O16_derived",  // 61,
        "F18_to_p_O17_derived",  // 62,
        "F18_to_He4_N14_derived",  // 63,
        "F19_to_p_O18_derived",  // 64,
        "F19_to_He4_N15_derived",  // 65,
        "Ne18_to_p_F17_derived",  // 66,
        "Ne18_to_He4_O14_derived",  // 67,
        "Ne19_to_p_F18_derived",  // 68,
        "Ne19_to_He4_O15_derived",  // 69,
        "Ne20_to_p_F19_derived",  // 70,
        "Ne20_to_He4_O16_derived",  // 71,
        "Mg22_to_He4_Ne18_derived",  // 72,
        "Mg24_to_He4_Ne20_derived",  // 73,
        "C12_to_3He4_derived",  // 74,
        "C12_He4_to_p_N15_derived",  // 75,
        "N14_He4_to_p_O17_derived",  // 76,
        "N15_He4_to_p_O18_derived",  // 77,
        "O15_He4_to_p_F18_derived",  // 78,
        "O16_p_to_He4_N13_derived",  // 79,
        "O16_He4_to_p_F19_derived",  // 80,
        "F17_p_to_He4_O14_derived",  // 81,
        "F17_He4_to_p_Ne20_derived",  // 82,
        "Na21_p_to_He4_Ne18_derived",  // 83,
        "Mg22_to_p_Na21_derived"  // 84,
    };

    template<std::uint8_t rate>
    AMREX_GPU_DEVICE AMREX_INLINE
    constexpr amrex::Real get_p_random() {
        return 0.0_rt;
    }

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
    constexpr int He4_index = 3;
    constexpr int NumNSERatePairs = 28;

    inline AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, NumNSERatePairs, 1, 10, amrex::Order::C> rate_pair_data {
        H1, C12, -1, N13, -1, -1, C12, N13, k_C12_p_to_N13_reaclib, k_N13_to_p_C12_derived,
        He4, C12, -1, O16, -1, -1, C12, O16, k_C12_He4_to_O16_reaclib, k_O16_to_He4_C12_derived,
        H1, C13, -1, N14, -1, -1, C13, N14, k_C13_p_to_N14_reaclib, k_N14_to_p_C13_derived,
        H1, N13, -1, O14, -1, -1, N13, O14, k_N13_p_to_O14_reaclib, k_O14_to_p_N13_derived,
        H1, N14, -1, O15, -1, -1, N14, O15, k_N14_p_to_O15_reaclib, k_O15_to_p_N14_derived,
        He4, N14, -1, F18, -1, -1, N14, F18, k_N14_He4_to_F18_reaclib, k_F18_to_He4_N14_derived,
        H1, N15, -1, O16, -1, -1, N15, O16, k_N15_p_to_O16_reaclib, k_O16_to_p_N15_derived,
        He4, N15, -1, F19, -1, -1, N15, F19, k_N15_He4_to_F19_reaclib, k_F19_to_He4_N15_derived,
        He4, O14, -1, Ne18, -1, -1, O14, Ne18, k_O14_He4_to_Ne18_reaclib, k_Ne18_to_He4_O14_derived,
        He4, O15, -1, Ne19, -1, -1, O15, Ne19, k_O15_He4_to_Ne19_reaclib, k_Ne19_to_He4_O15_derived,
        H1, O16, -1, F17, -1, -1, O16, F17, k_O16_p_to_F17_reaclib, k_F17_to_p_O16_derived,
        He4, O16, -1, Ne20, -1, -1, O16, Ne20, k_O16_He4_to_Ne20_reaclib, k_Ne20_to_He4_O16_derived,
        H1, O17, -1, F18, -1, -1, O17, F18, k_O17_p_to_F18_reaclib, k_F18_to_p_O17_derived,
        H1, O18, -1, F19, -1, -1, O18, F19, k_O18_p_to_F19_reaclib, k_F19_to_p_O18_derived,
        H1, F17, -1, Ne18, -1, -1, F17, Ne18, k_F17_p_to_Ne18_reaclib, k_Ne18_to_p_F17_derived,
        H1, F18, -1, Ne19, -1, -1, F18, Ne19, k_F18_p_to_Ne19_reaclib, k_Ne19_to_p_F18_derived,
        H1, F19, -1, Ne20, -1, -1, F19, Ne20, k_F19_p_to_Ne20_reaclib, k_Ne20_to_p_F19_derived,
        He4, Ne20, -1, Mg24, -1, -1, Ne20, Mg24, k_Ne20_He4_to_Mg24_reaclib, k_Mg24_to_He4_Ne20_derived,
        He4, N13, -1, H1, O16, -1, N13, O16, k_N13_He4_to_p_O16_reaclib, k_O16_p_to_He4_N13_derived,
        H1, N15, -1, He4, C12, -1, N15, C12, k_N15_p_to_He4_C12_reaclib, k_C12_He4_to_p_N15_derived,
        He4, O14, -1, H1, F17, -1, O14, F17, k_O14_He4_to_p_F17_reaclib, k_F17_p_to_He4_O14_derived,
        H1, O17, -1, He4, N14, -1, O17, N14, k_O17_p_to_He4_N14_reaclib, k_N14_He4_to_p_O17_derived,
        H1, O18, -1, He4, N15, -1, O18, N15, k_O18_p_to_He4_N15_reaclib, k_N15_He4_to_p_O18_derived,
        H1, F18, -1, He4, O15, -1, F18, O15, k_F18_p_to_He4_O15_reaclib, k_O15_He4_to_p_F18_derived,
        H1, F19, -1, He4, O16, -1, F19, O16, k_F19_p_to_He4_O16_reaclib, k_O16_He4_to_p_F19_derived,
        He4, He4, He4, C12, -1, -1, C12, -1, k_He4_He4_He4_to_C12_reaclib, k_C12_to_3He4_derived,
        He4, F17, -1, H1, Ne20, -1, F17, Ne20, k_F17_He4_to_p_Ne20_derived, k_Ne20_p_to_He4_F17_reaclib,
        Ne18, He4, -1, Mg22, -1, -1, Ne18, Mg22, k_Ne18_He4_to_Mg22_approx, k_Mg22_to_He4_Ne18_derived
    };
}
#endif

#endif
