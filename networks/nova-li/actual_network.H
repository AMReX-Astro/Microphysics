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
        else if constexpr (spec == H2) {
            return 2.2245662689997516_rt;
        }
        else if constexpr (spec == He3) {
            return 7.718041347999588_rt;
        }
        else if constexpr (spec == He4) {
            return 28.295662457999697_rt;
        }
        else if constexpr (spec == Li7) {
            return 39.245080591999795_rt;
        }
        else if constexpr (spec == Be7) {
            return 37.60083855599987_rt;
        }
        else if constexpr (spec == B8) {
            return 37.7372096199997_rt;
        }
        else if constexpr (spec == B11) {
            return 76.20505591999972_rt;
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
        else if constexpr (spec == Ne21) {
            return 167.40598973999658_rt;
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
        else if constexpr (spec == H2) {
            return 3.3444946862179e-24_rt;
        }
        else if constexpr (spec == He3) {
            return 5.008234515140786e-24_rt;
        }
        else if constexpr (spec == He4) {
            return 6.646479071584587e-24_rt;
        }
        else if constexpr (spec == Li7) {
            return 1.1650347794645127e-23_rt;
        }
        else if constexpr (spec == Be7) {
            return 1.1651884262042064e-23_rt;
        }
        else if constexpr (spec == B8) {
            return 1.332517399630266e-23_rt;
        }
        else if constexpr (spec == B11) {
            return 1.828138132559934e-23_rt;
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
        else if constexpr (spec == Ne21) {
            return 3.4861102572650884e-23_rt;
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
        k_Be7_to_Li7_electron_capture_reaclib = 1,
        k_N13_to_C13_beta_pos_reaclib = 2,
        k_O14_to_N14_beta_pos_reaclib = 3,
        k_O15_to_N15_beta_pos_reaclib = 4,
        k_B8_to_p_Be7_reaclib = 5,
        k_B8_to_He4_He4_beta_pos_reaclib = 6,
        k_F17_to_p_O16_reaclib = 7,
        k_p_p_to_d_beta_pos_reaclib = 8,
        k_p_p_to_d_electron_capture_reaclib = 9,
        k_d_p_to_He3_reaclib = 10,
        k_d_d_to_He4_reaclib = 11,
        k_He3_p_to_He4_beta_pos_reaclib = 12,
        k_He3_He4_to_Be7_reaclib = 13,
        k_Li7_He4_to_B11_reaclib = 14,
        k_Be7_p_to_B8_reaclib = 15,
        k_B11_p_to_C12_reaclib = 16,
        k_C12_p_to_N13_reaclib = 17,
        k_C12_He4_to_O16_reaclib = 18,
        k_C13_p_to_N14_reaclib = 19,
        k_N13_p_to_O14_reaclib = 20,
        k_N14_p_to_O15_reaclib = 21,
        k_N14_He4_to_F18_reaclib = 22,
        k_N15_p_to_O16_reaclib = 23,
        k_N15_He4_to_F19_reaclib = 24,
        k_O14_He4_to_Ne18_reaclib = 25,
        k_O15_He4_to_Ne19_reaclib = 26,
        k_O16_p_to_F17_reaclib = 27,
        k_O16_He4_to_Ne20_reaclib = 28,
        k_O17_p_to_F18_reaclib = 29,
        k_O17_He4_to_Ne21_reaclib = 30,
        k_O18_p_to_F19_reaclib = 31,
        k_F17_p_to_Ne18_reaclib = 32,
        k_F18_p_to_Ne19_reaclib = 33,
        k_F19_p_to_Ne20_reaclib = 34,
        k_He3_d_to_p_He4_reaclib = 35,
        k_Li7_p_to_He4_He4_reaclib = 36,
        k_N13_He4_to_p_O16_reaclib = 37,
        k_N15_p_to_He4_C12_reaclib = 38,
        k_O14_He4_to_p_F17_reaclib = 39,
        k_O17_p_to_He4_N14_reaclib = 40,
        k_O18_p_to_He4_N15_reaclib = 41,
        k_F17_He4_to_p_Ne20_reaclib = 42,
        k_F18_p_to_He4_O15_reaclib = 43,
        k_F18_He4_to_p_Ne21_reaclib = 44,
        k_F19_p_to_He4_O16_reaclib = 45,
        k_He3_He3_to_p_p_He4_reaclib = 46,
        k_Be7_d_to_p_He4_He4_reaclib = 47,
        k_B11_p_to_3He4_reaclib = 48,
        k_Be7_He3_to_p_p_He4_He4_reaclib = 49,
        k_He4_He4_He4_to_C12_reaclib = 50,
        k_F17_to_O17_electron_capture_weaktab = 51,
        k_F18_to_O18_electron_capture_weaktab = 52,
        k_Ne18_to_F18_electron_capture_weaktab = 53,
        k_Ne19_to_F19_electron_capture_weaktab = 54,
        NumRates = k_Ne19_to_F19_electron_capture_weaktab
    };

    enum ScreenPairs :
    std::uint8_t
    {
        k_p_Li7 = 1,
        k_p_F18 = 2,
        k_p_C12 = 3,
        k_He3_He3 = 4,
        k_He4_O16 = 5,
        k_He4_F17 = 6,
        k_He4_N14 = 7,
        k_He3_He4 = 8,
        k_p_F19 = 9,
        k_p_C13 = 10,
        k_p_B11 = 11,
        k_He3_Be7 = 12,
        k_He4_Li7 = 13,
        k_He4_F18 = 14,
        k_He4_C12 = 15,
        k_p_p = 16,
        k_p_O17 = 17,
        k_p_N15 = 18,
        k_p_O18 = 19,
        k_p_d = 20,
        k_p_He3 = 21,
        k_He4_O17 = 22,
        k_d_d = 23,
        k_He4_N15 = 24,
        k_d_He3 = 25,
        k_p_Be7 = 26,
        k_d_Be7 = 27,
        k_p_N13 = 28,
        k_He4_O14 = 29,
        k_He4_He4 = 30,
        k_He4_Be8 = 31,
        k_p_O16 = 32,
        k_p_N14 = 33,
        k_p_F17 = 34,
        k_He4_O15 = 35,
        k_He4_N13 = 36,
        NumScreenPairs = k_He4_N13
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "Be7_to_Li7_electron_capture_reaclib",  // 1,
        "N13_to_C13_beta_pos_reaclib",  // 2,
        "O14_to_N14_beta_pos_reaclib",  // 3,
        "O15_to_N15_beta_pos_reaclib",  // 4,
        "B8_to_p_Be7_reaclib",  // 5,
        "B8_to_He4_He4_beta_pos_reaclib",  // 6,
        "F17_to_p_O16_reaclib",  // 7,
        "p_p_to_d_beta_pos_reaclib",  // 8,
        "p_p_to_d_electron_capture_reaclib",  // 9,
        "d_p_to_He3_reaclib",  // 10,
        "d_d_to_He4_reaclib",  // 11,
        "He3_p_to_He4_beta_pos_reaclib",  // 12,
        "He3_He4_to_Be7_reaclib",  // 13,
        "Li7_He4_to_B11_reaclib",  // 14,
        "Be7_p_to_B8_reaclib",  // 15,
        "B11_p_to_C12_reaclib",  // 16,
        "C12_p_to_N13_reaclib",  // 17,
        "C12_He4_to_O16_reaclib",  // 18,
        "C13_p_to_N14_reaclib",  // 19,
        "N13_p_to_O14_reaclib",  // 20,
        "N14_p_to_O15_reaclib",  // 21,
        "N14_He4_to_F18_reaclib",  // 22,
        "N15_p_to_O16_reaclib",  // 23,
        "N15_He4_to_F19_reaclib",  // 24,
        "O14_He4_to_Ne18_reaclib",  // 25,
        "O15_He4_to_Ne19_reaclib",  // 26,
        "O16_p_to_F17_reaclib",  // 27,
        "O16_He4_to_Ne20_reaclib",  // 28,
        "O17_p_to_F18_reaclib",  // 29,
        "O17_He4_to_Ne21_reaclib",  // 30,
        "O18_p_to_F19_reaclib",  // 31,
        "F17_p_to_Ne18_reaclib",  // 32,
        "F18_p_to_Ne19_reaclib",  // 33,
        "F19_p_to_Ne20_reaclib",  // 34,
        "He3_d_to_p_He4_reaclib",  // 35,
        "Li7_p_to_He4_He4_reaclib",  // 36,
        "N13_He4_to_p_O16_reaclib",  // 37,
        "N15_p_to_He4_C12_reaclib",  // 38,
        "O14_He4_to_p_F17_reaclib",  // 39,
        "O17_p_to_He4_N14_reaclib",  // 40,
        "O18_p_to_He4_N15_reaclib",  // 41,
        "F17_He4_to_p_Ne20_reaclib",  // 42,
        "F18_p_to_He4_O15_reaclib",  // 43,
        "F18_He4_to_p_Ne21_reaclib",  // 44,
        "F19_p_to_He4_O16_reaclib",  // 45,
        "He3_He3_to_p_p_He4_reaclib",  // 46,
        "Be7_d_to_p_He4_He4_reaclib",  // 47,
        "B11_p_to_3He4_reaclib",  // 48,
        "Be7_He3_to_p_p_He4_He4_reaclib",  // 49,
        "He4_He4_He4_to_C12_reaclib",  // 50,
        "F17_to_O17_electron_capture_weaktab",  // 51,
        "F18_to_O18_electron_capture_weaktab",  // 52,
        "Ne18_to_F18_electron_capture_weaktab",  // 53,
        "Ne19_to_F19_electron_capture_weaktab"  // 54,
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
    constexpr int He4_index = 4;
    constexpr int NumNSERatePairs = 2;

    inline AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, NumNSERatePairs, 1, 10, amrex::Order::C> rate_pair_data {
        H1, Be7, -1, B8, -1, -1, Be7, B8, k_Be7_p_to_B8_reaclib, k_B8_to_p_Be7_reaclib,
        H1, O16, -1, F17, -1, -1, O16, F17, k_O16_p_to_F17_reaclib, k_F17_to_p_O16_reaclib
    };
}
#endif

#endif
