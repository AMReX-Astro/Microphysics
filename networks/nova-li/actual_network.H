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
        k_Be7_to_Li7_weak_electron_capture = 1,
        k_N13_to_C13_weak_wc12 = 2,
        k_O14_to_N14_weak_wc12 = 3,
        k_O15_to_N15_weak_wc12 = 4,
        k_B8_to_p_Be7 = 5,
        k_B8_to_He4_He4_weak_wc12 = 6,
        k_F17_to_p_O16 = 7,
        k_p_p_to_d_weak_bet_pos_ = 8,
        k_p_p_to_d_weak_electron_capture = 9,
        k_p_d_to_He3 = 10,
        k_d_d_to_He4 = 11,
        k_p_He3_to_He4_weak_bet_pos_ = 12,
        k_He4_He3_to_Be7 = 13,
        k_He4_Li7_to_B11 = 14,
        k_p_Be7_to_B8 = 15,
        k_p_B11_to_C12 = 16,
        k_p_C12_to_N13 = 17,
        k_He4_C12_to_O16 = 18,
        k_p_C13_to_N14 = 19,
        k_p_N13_to_O14 = 20,
        k_p_N14_to_O15 = 21,
        k_He4_N14_to_F18 = 22,
        k_p_N15_to_O16 = 23,
        k_He4_N15_to_F19 = 24,
        k_He4_O14_to_Ne18 = 25,
        k_He4_O15_to_Ne19 = 26,
        k_p_O16_to_F17 = 27,
        k_He4_O16_to_Ne20 = 28,
        k_p_O17_to_F18 = 29,
        k_He4_O17_to_Ne21 = 30,
        k_p_O18_to_F19 = 31,
        k_p_F17_to_Ne18 = 32,
        k_p_F18_to_Ne19 = 33,
        k_p_F19_to_Ne20 = 34,
        k_d_He3_to_p_He4 = 35,
        k_p_Li7_to_He4_He4 = 36,
        k_He4_N13_to_p_O16 = 37,
        k_p_N15_to_He4_C12 = 38,
        k_He4_O14_to_p_F17 = 39,
        k_p_O17_to_He4_N14 = 40,
        k_p_O18_to_He4_N15 = 41,
        k_He4_F17_to_p_Ne20 = 42,
        k_p_F18_to_He4_O15 = 43,
        k_He4_F18_to_p_Ne21 = 44,
        k_p_F19_to_He4_O16 = 45,
        k_He3_He3_to_p_p_He4 = 46,
        k_d_Be7_to_p_He4_He4 = 47,
        k_p_B11_to_He4_He4_He4 = 48,
        k_He3_Be7_to_p_p_He4_He4 = 49,
        k_He4_He4_He4_to_C12 = 50,
        k_F17_to_O17 = 51,
        k_F18_to_O18 = 52,
        k_Ne18_to_F18 = 53,
        k_Ne19_to_F19 = 54,
        NumRates = k_Ne19_to_F19
    };

    // number of reaclib rates

    const int NrateReaclib = 50;

    // number of tabular rates

    const int NrateTabular = 4;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "Be7_to_Li7_weak_electron_capture",  // 1,
        "N13_to_C13_weak_wc12",  // 2,
        "O14_to_N14_weak_wc12",  // 3,
        "O15_to_N15_weak_wc12",  // 4,
        "B8_to_p_Be7",  // 5,
        "B8_to_He4_He4_weak_wc12",  // 6,
        "F17_to_p_O16",  // 7,
        "p_p_to_d_weak_bet_pos_",  // 8,
        "p_p_to_d_weak_electron_capture",  // 9,
        "p_d_to_He3",  // 10,
        "d_d_to_He4",  // 11,
        "p_He3_to_He4_weak_bet_pos_",  // 12,
        "He4_He3_to_Be7",  // 13,
        "He4_Li7_to_B11",  // 14,
        "p_Be7_to_B8",  // 15,
        "p_B11_to_C12",  // 16,
        "p_C12_to_N13",  // 17,
        "He4_C12_to_O16",  // 18,
        "p_C13_to_N14",  // 19,
        "p_N13_to_O14",  // 20,
        "p_N14_to_O15",  // 21,
        "He4_N14_to_F18",  // 22,
        "p_N15_to_O16",  // 23,
        "He4_N15_to_F19",  // 24,
        "He4_O14_to_Ne18",  // 25,
        "He4_O15_to_Ne19",  // 26,
        "p_O16_to_F17",  // 27,
        "He4_O16_to_Ne20",  // 28,
        "p_O17_to_F18",  // 29,
        "He4_O17_to_Ne21",  // 30,
        "p_O18_to_F19",  // 31,
        "p_F17_to_Ne18",  // 32,
        "p_F18_to_Ne19",  // 33,
        "p_F19_to_Ne20",  // 34,
        "d_He3_to_p_He4",  // 35,
        "p_Li7_to_He4_He4",  // 36,
        "He4_N13_to_p_O16",  // 37,
        "p_N15_to_He4_C12",  // 38,
        "He4_O14_to_p_F17",  // 39,
        "p_O17_to_He4_N14",  // 40,
        "p_O18_to_He4_N15",  // 41,
        "He4_F17_to_p_Ne20",  // 42,
        "p_F18_to_He4_O15",  // 43,
        "He4_F18_to_p_Ne21",  // 44,
        "p_F19_to_He4_O16",  // 45,
        "He3_He3_to_p_p_He4",  // 46,
        "d_Be7_to_p_He4_He4",  // 47,
        "p_B11_to_He4_He4_He4",  // 48,
        "He3_Be7_to_p_p_He4_He4",  // 49,
        "He4_He4_He4_to_C12",  // 50,
        "F17_to_O17",  // 51,
        "F18_to_O18",  // 52,
        "Ne18_to_F18",  // 53,
        "Ne19_to_F19"  // 54,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int H1_index = 0;
    constexpr int N_index = -1;
    constexpr int He4_index = 3;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, amrex::Order::C> rate_indices;
}
#endif

#endif
