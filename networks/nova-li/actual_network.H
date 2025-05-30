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
        k_B8_to_He4_He4_weak_wc12 = 10,
        k_p_p_to_d_weak_bet_pos_ = 11,
        k_p_p_to_d_weak_electron_capture = 12,
        k_p_d_to_He3 = 13,
        k_d_d_to_He4 = 14,
        k_He4_d_to_Li6 = 15,
        k_p_He3_to_He4_weak_bet_pos_ = 16,
        k_He4_He3_to_Be7 = 17,
        k_p_Li6_to_Be7 = 18,
        k_p_Be7_to_B8 = 19,
        k_p_C12_to_N13 = 20,
        k_He4_C12_to_O16 = 21,
        k_p_C13_to_N14 = 22,
        k_p_N13_to_O14 = 23,
        k_p_N14_to_O15 = 24,
        k_He4_N14_to_F18 = 25,
        k_p_N15_to_O16 = 26,
        k_He4_N15_to_F19 = 27,
        k_He4_O14_to_Ne18 = 28,
        k_He4_O15_to_Ne19 = 29,
        k_p_O16_to_F17 = 30,
        k_He4_O16_to_Ne20 = 31,
        k_p_O17_to_F18 = 32,
        k_He4_O17_to_Ne21 = 33,
        k_p_O18_to_F19 = 34,
        k_p_F17_to_Ne18 = 35,
        k_p_F18_to_Ne19 = 36,
        k_p_F19_to_Ne20 = 37,
        k_d_He3_to_p_He4 = 38,
        k_p_Li6_to_He4_He3 = 39,
        k_d_Li6_to_p_Li7 = 40,
        k_p_Li7_to_He4_He4 = 41,
        k_C12_C12_to_He4_Ne20 = 42,
        k_He4_N13_to_p_O16 = 43,
        k_p_N15_to_He4_C12 = 44,
        k_He4_O14_to_p_F17 = 45,
        k_p_O17_to_He4_N14 = 46,
        k_p_O18_to_He4_N15 = 47,
        k_p_F18_to_He4_O15 = 48,
        k_He4_F18_to_p_Ne21 = 49,
        k_p_F19_to_He4_O16 = 50,
        k_p_Ne20_to_He4_F17 = 51,
        k_He3_He3_to_p_p_He4 = 52,
        k_d_Be7_to_p_He4_He4 = 53,
        k_He3_Be7_to_p_p_He4_He4 = 54,
        k_He4_He4_He4_to_C12 = 55,
        NumRates = k_He4_He4_He4_to_C12
    };

    // number of reaclib rates

    const int NrateReaclib = 55;

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
        "B8_to_He4_He4_weak_wc12",  // 10,
        "p_p_to_d_weak_bet_pos_",  // 11,
        "p_p_to_d_weak_electron_capture",  // 12,
        "p_d_to_He3",  // 13,
        "d_d_to_He4",  // 14,
        "He4_d_to_Li6",  // 15,
        "p_He3_to_He4_weak_bet_pos_",  // 16,
        "He4_He3_to_Be7",  // 17,
        "p_Li6_to_Be7",  // 18,
        "p_Be7_to_B8",  // 19,
        "p_C12_to_N13",  // 20,
        "He4_C12_to_O16",  // 21,
        "p_C13_to_N14",  // 22,
        "p_N13_to_O14",  // 23,
        "p_N14_to_O15",  // 24,
        "He4_N14_to_F18",  // 25,
        "p_N15_to_O16",  // 26,
        "He4_N15_to_F19",  // 27,
        "He4_O14_to_Ne18",  // 28,
        "He4_O15_to_Ne19",  // 29,
        "p_O16_to_F17",  // 30,
        "He4_O16_to_Ne20",  // 31,
        "p_O17_to_F18",  // 32,
        "He4_O17_to_Ne21",  // 33,
        "p_O18_to_F19",  // 34,
        "p_F17_to_Ne18",  // 35,
        "p_F18_to_Ne19",  // 36,
        "p_F19_to_Ne20",  // 37,
        "d_He3_to_p_He4",  // 38,
        "p_Li6_to_He4_He3",  // 39,
        "d_Li6_to_p_Li7",  // 40,
        "p_Li7_to_He4_He4",  // 41,
        "C12_C12_to_He4_Ne20",  // 42,
        "He4_N13_to_p_O16",  // 43,
        "p_N15_to_He4_C12",  // 44,
        "He4_O14_to_p_F17",  // 45,
        "p_O17_to_He4_N14",  // 46,
        "p_O18_to_He4_N15",  // 47,
        "p_F18_to_He4_O15",  // 48,
        "He4_F18_to_p_Ne21",  // 49,
        "p_F19_to_He4_O16",  // 50,
        "p_Ne20_to_He4_F17",  // 51,
        "He3_He3_to_p_p_He4",  // 52,
        "d_Be7_to_p_He4_He4",  // 53,
        "He3_Be7_to_p_p_He4_He4",  // 54,
        "He4_He4_He4_to_C12"  // 55,
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
