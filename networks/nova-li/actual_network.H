#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Loop.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex;

void actual_network_init();

const std::string network_name = "pynucastro-cxx";

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
        else if constexpr (spec == Li6) {
            return 31.99398709199886_rt;
        }
        else if constexpr (spec == Li7) {
            return 39.245080591999795_rt;
        }
        else if constexpr (spec == Be7) {
            return 37.60083855599987_rt;
        }
        else if constexpr (spec == Be8) {
            return 56.499486655999135_rt;
        }
        else if constexpr (spec == B8) {
            return 37.7372096199997_rt;
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
        static_assert(spec >= 1 && spec <= NumSpec);

        constexpr amrex::Real A = NetworkProperties::aion(spec);
        constexpr amrex::Real Z = NetworkProperties::zion(spec);

        return (A - Z) * C::Legacy::m_n + Z * (C::Legacy::m_p + C::Legacy::m_e) - bion<spec>() * C::Legacy::MeV2gr;
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
        k_Be7_to_Li7_weak_electron_capture = 1,
        k_B8_to_Be8_weak_wc17 = 2,
        k_N13_to_C13_weak_wc12 = 3,
        k_O14_to_N14_weak_wc12 = 4,
        k_O15_to_N15_weak_wc12 = 5,
        k_F17_to_O17_weak_wc12 = 6,
        k_F18_to_O18_weak_wc12 = 7,
        k_Ne18_to_F18_weak_wc12 = 8,
        k_Ne19_to_F19_weak_wc12 = 9,
        k_He3_to_p_d = 10,
        k_He4_to_d_d = 11,
        k_Li6_to_He4_d = 12,
        k_Be7_to_p_Li6 = 13,
        k_Be7_to_He4_He3 = 14,
        k_B8_to_p_Be7 = 15,
        k_B8_to_He4_He4_weak_wc12 = 16,
        k_N13_to_p_C12 = 17,
        k_N14_to_p_C13 = 18,
        k_O14_to_p_N13 = 19,
        k_O15_to_p_N14 = 20,
        k_O16_to_p_N15 = 21,
        k_O16_to_He4_C12 = 22,
        k_F17_to_p_O16 = 23,
        k_F18_to_p_O17 = 24,
        k_F18_to_He4_N14 = 25,
        k_F19_to_p_O18 = 26,
        k_F19_to_He4_N15 = 27,
        k_Ne18_to_p_F17 = 28,
        k_Ne18_to_He4_O14 = 29,
        k_Ne19_to_p_F18 = 30,
        k_Ne19_to_He4_O15 = 31,
        k_Ne20_to_p_F19 = 32,
        k_Ne20_to_He4_O16 = 33,
        k_Ne21_to_He4_O17 = 34,
        k_C12_to_He4_He4_He4 = 35,
        k_p_p_to_d_weak_bet_pos_ = 36,
        k_p_p_to_d_weak_electron_capture = 37,
        k_p_d_to_He3 = 38,
        k_d_d_to_He4 = 39,
        k_He4_d_to_Li6 = 40,
        k_p_He3_to_He4_weak_bet_pos_ = 41,
        k_He4_He3_to_Be7 = 42,
        k_p_Li6_to_Be7 = 43,
        k_p_Be7_to_B8 = 44,
        k_p_C12_to_N13 = 45,
        k_He4_C12_to_O16 = 46,
        k_p_C13_to_N14 = 47,
        k_p_N13_to_O14 = 48,
        k_p_N14_to_O15 = 49,
        k_He4_N14_to_F18 = 50,
        k_p_N15_to_O16 = 51,
        k_He4_N15_to_F19 = 52,
        k_He4_O14_to_Ne18 = 53,
        k_He4_O15_to_Ne19 = 54,
        k_p_O16_to_F17 = 55,
        k_He4_O16_to_Ne20 = 56,
        k_p_O17_to_F18 = 57,
        k_He4_O17_to_Ne21 = 58,
        k_p_O18_to_F19 = 59,
        k_p_F17_to_Ne18 = 60,
        k_p_F18_to_Ne19 = 61,
        k_p_F19_to_Ne20 = 62,
        k_d_He3_to_p_He4 = 63,
        k_He4_He3_to_p_Li6 = 64,
        k_p_He4_to_d_He3 = 65,
        k_He4_He4_to_p_Li7 = 66,
        k_p_Li6_to_He4_He3 = 67,
        k_d_Li6_to_p_Li7 = 68,
        k_p_Li7_to_d_Li6 = 69,
        k_p_Li7_to_He4_He4 = 70,
        k_He4_C12_to_p_N15 = 71,
        k_C12_C12_to_He4_Ne20 = 72,
        k_He4_N13_to_p_O16 = 73,
        k_He4_N14_to_p_O17 = 74,
        k_p_N15_to_He4_C12 = 75,
        k_He4_N15_to_p_O18 = 76,
        k_He4_O14_to_p_F17 = 77,
        k_He4_O15_to_p_F18 = 78,
        k_p_O16_to_He4_N13 = 79,
        k_He4_O16_to_p_F19 = 80,
        k_p_O17_to_He4_N14 = 81,
        k_p_O18_to_He4_N15 = 82,
        k_p_F17_to_He4_O14 = 83,
        k_He4_F17_to_p_Ne20 = 84,
        k_p_F18_to_He4_O15 = 85,
        k_He4_F18_to_p_Ne21 = 86,
        k_p_F19_to_He4_O16 = 87,
        k_p_Ne20_to_He4_F17 = 88,
        k_He4_Ne20_to_C12_C12 = 89,
        k_p_Ne21_to_He4_F18 = 90,
        k_He3_He3_to_p_p_He4 = 91,
        k_d_Be7_to_p_He4_He4 = 92,
        k_He3_Be7_to_p_p_He4_He4 = 93,
        k_He4_He4_He4_to_C12 = 94,
        k_p_p_He4_to_He3_He3 = 95,
        k_p_He4_He4_to_d_Be7 = 96,
        k_p_p_He4_He4_to_He3_Be7 = 97,
        NumRates = k_p_p_He4_He4_to_He3_Be7
    };

    // number of reaclib rates

    const int NrateReaclib = 97;

    // number of tabular rates

    const int NrateTabular = 0;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "Be7_to_Li7_weak_electron_capture",  // 1,
        "B8_to_Be8_weak_wc17",  // 2,
        "N13_to_C13_weak_wc12",  // 3,
        "O14_to_N14_weak_wc12",  // 4,
        "O15_to_N15_weak_wc12",  // 5,
        "F17_to_O17_weak_wc12",  // 6,
        "F18_to_O18_weak_wc12",  // 7,
        "Ne18_to_F18_weak_wc12",  // 8,
        "Ne19_to_F19_weak_wc12",  // 9,
        "He3_to_p_d",  // 10,
        "He4_to_d_d",  // 11,
        "Li6_to_He4_d",  // 12,
        "Be7_to_p_Li6",  // 13,
        "Be7_to_He4_He3",  // 14,
        "B8_to_p_Be7",  // 15,
        "B8_to_He4_He4_weak_wc12",  // 16,
        "N13_to_p_C12",  // 17,
        "N14_to_p_C13",  // 18,
        "O14_to_p_N13",  // 19,
        "O15_to_p_N14",  // 20,
        "O16_to_p_N15",  // 21,
        "O16_to_He4_C12",  // 22,
        "F17_to_p_O16",  // 23,
        "F18_to_p_O17",  // 24,
        "F18_to_He4_N14",  // 25,
        "F19_to_p_O18",  // 26,
        "F19_to_He4_N15",  // 27,
        "Ne18_to_p_F17",  // 28,
        "Ne18_to_He4_O14",  // 29,
        "Ne19_to_p_F18",  // 30,
        "Ne19_to_He4_O15",  // 31,
        "Ne20_to_p_F19",  // 32,
        "Ne20_to_He4_O16",  // 33,
        "Ne21_to_He4_O17",  // 34,
        "C12_to_He4_He4_He4",  // 35,
        "p_p_to_d_weak_bet_pos_",  // 36,
        "p_p_to_d_weak_electron_capture",  // 37,
        "p_d_to_He3",  // 38,
        "d_d_to_He4",  // 39,
        "He4_d_to_Li6",  // 40,
        "p_He3_to_He4_weak_bet_pos_",  // 41,
        "He4_He3_to_Be7",  // 42,
        "p_Li6_to_Be7",  // 43,
        "p_Be7_to_B8",  // 44,
        "p_C12_to_N13",  // 45,
        "He4_C12_to_O16",  // 46,
        "p_C13_to_N14",  // 47,
        "p_N13_to_O14",  // 48,
        "p_N14_to_O15",  // 49,
        "He4_N14_to_F18",  // 50,
        "p_N15_to_O16",  // 51,
        "He4_N15_to_F19",  // 52,
        "He4_O14_to_Ne18",  // 53,
        "He4_O15_to_Ne19",  // 54,
        "p_O16_to_F17",  // 55,
        "He4_O16_to_Ne20",  // 56,
        "p_O17_to_F18",  // 57,
        "He4_O17_to_Ne21",  // 58,
        "p_O18_to_F19",  // 59,
        "p_F17_to_Ne18",  // 60,
        "p_F18_to_Ne19",  // 61,
        "p_F19_to_Ne20",  // 62,
        "d_He3_to_p_He4",  // 63,
        "He4_He3_to_p_Li6",  // 64,
        "p_He4_to_d_He3",  // 65,
        "He4_He4_to_p_Li7",  // 66,
        "p_Li6_to_He4_He3",  // 67,
        "d_Li6_to_p_Li7",  // 68,
        "p_Li7_to_d_Li6",  // 69,
        "p_Li7_to_He4_He4",  // 70,
        "He4_C12_to_p_N15",  // 71,
        "C12_C12_to_He4_Ne20",  // 72,
        "He4_N13_to_p_O16",  // 73,
        "He4_N14_to_p_O17",  // 74,
        "p_N15_to_He4_C12",  // 75,
        "He4_N15_to_p_O18",  // 76,
        "He4_O14_to_p_F17",  // 77,
        "He4_O15_to_p_F18",  // 78,
        "p_O16_to_He4_N13",  // 79,
        "He4_O16_to_p_F19",  // 80,
        "p_O17_to_He4_N14",  // 81,
        "p_O18_to_He4_N15",  // 82,
        "p_F17_to_He4_O14",  // 83,
        "He4_F17_to_p_Ne20",  // 84,
        "p_F18_to_He4_O15",  // 85,
        "He4_F18_to_p_Ne21",  // 86,
        "p_F19_to_He4_O16",  // 87,
        "p_Ne20_to_He4_F17",  // 88,
        "He4_Ne20_to_C12_C12",  // 89,
        "p_Ne21_to_He4_F18",  // 90,
        "He3_He3_to_p_p_He4",  // 91,
        "d_Be7_to_p_He4_He4",  // 92,
        "He3_Be7_to_p_p_He4_He4",  // 93,
        "He4_He4_He4_to_C12",  // 94,
        "p_p_He4_to_He3_He3",  // 95,
        "p_He4_He4_to_d_Be7",  // 96,
        "p_p_He4_He4_to_He3_Be7"  // 97,
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

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
