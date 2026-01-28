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
        else if constexpr (spec == Be7) {
            return 37.60083855599987_rt;
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
        else if constexpr (spec == F17) {
            return 128.21958437599824_rt;
        }
        else if constexpr (spec == F18) {
            return 137.36950247599816_rt;
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
        else if constexpr (spec == Be7) {
            return 1.1651884262042064e-23_rt;
        }
        else if constexpr (spec == B8) {
            return 1.332517399630266e-23_rt;
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
        else if constexpr (spec == F17) {
            return 2.8232643353472436e-23_rt;
        }
        else if constexpr (spec == F18) {
            return 2.989125964092377e-23_rt;
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
        k_F17_to_O17_reaclib = 4,
        k_B8_to_He4_He4_reaclib = 5,
        k_p_p_to_d_reaclib_bet_pos = 6,
        k_p_p_to_d_reaclib_electron_capture = 7,
        k_p_d_to_He3_reaclib = 8,
        k_d_d_to_He4_reaclib = 9,
        k_p_He3_to_He4_reaclib = 10,
        k_He4_He3_to_Be7_reaclib = 11,
        k_p_Be7_to_B8_reaclib = 12,
        k_p_C12_to_N13_reaclib = 13,
        k_He4_C12_to_O16_reaclib = 14,
        k_p_C13_to_N14_reaclib = 15,
        k_p_N13_to_O14_reaclib = 16,
        k_p_N14_to_O15_reaclib = 17,
        k_He4_N14_to_F18_reaclib = 18,
        k_p_N15_to_O16_reaclib = 19,
        k_p_O16_to_F17_reaclib = 20,
        k_p_O17_to_F18_reaclib = 21,
        k_d_He3_to_p_He4_reaclib = 22,
        k_He4_N13_to_p_O16_reaclib = 23,
        k_p_N15_to_He4_C12_reaclib = 24,
        k_He4_O14_to_p_F17_reaclib = 25,
        k_p_O17_to_He4_N14_reaclib = 26,
        k_p_F18_to_He4_O15_reaclib = 27,
        k_He3_He3_to_p_p_He4_reaclib = 28,
        k_d_Be7_to_p_He4_He4_reaclib = 29,
        k_He3_Be7_to_p_p_He4_He4_reaclib = 30,
        k_He4_He4_He4_to_C12_reaclib = 31,
        NumRates = k_He4_He4_He4_to_C12_reaclib
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "N13_to_C13_reaclib",  // 1,
        "O14_to_N14_reaclib",  // 2,
        "O15_to_N15_reaclib",  // 3,
        "F17_to_O17_reaclib",  // 4,
        "B8_to_He4_He4_reaclib",  // 5,
        "p_p_to_d_reaclib_bet_pos",  // 6,
        "p_p_to_d_reaclib_electron_capture",  // 7,
        "p_d_to_He3_reaclib",  // 8,
        "d_d_to_He4_reaclib",  // 9,
        "p_He3_to_He4_reaclib",  // 10,
        "He4_He3_to_Be7_reaclib",  // 11,
        "p_Be7_to_B8_reaclib",  // 12,
        "p_C12_to_N13_reaclib",  // 13,
        "He4_C12_to_O16_reaclib",  // 14,
        "p_C13_to_N14_reaclib",  // 15,
        "p_N13_to_O14_reaclib",  // 16,
        "p_N14_to_O15_reaclib",  // 17,
        "He4_N14_to_F18_reaclib",  // 18,
        "p_N15_to_O16_reaclib",  // 19,
        "p_O16_to_F17_reaclib",  // 20,
        "p_O17_to_F18_reaclib",  // 21,
        "d_He3_to_p_He4_reaclib",  // 22,
        "He4_N13_to_p_O16_reaclib",  // 23,
        "p_N15_to_He4_C12_reaclib",  // 24,
        "He4_O14_to_p_F17_reaclib",  // 25,
        "p_O17_to_He4_N14_reaclib",  // 26,
        "p_F18_to_He4_O15_reaclib",  // 27,
        "He3_He3_to_p_p_He4_reaclib",  // 28,
        "d_Be7_to_p_He4_He4_reaclib",  // 29,
        "He3_Be7_to_p_p_He4_He4_reaclib",  // 30,
        "He4_He4_He4_to_C12_reaclib"  // 31,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    // A RatePair data table, 2D array of shape (NumNSERatePairs, 8).
    // Each row corresponds to a RatePair returned by _get_nse_rate_pairs().
    // Nuclei indices follow NetworkSpecies and rate indices follow NetworkRates.

    //  - Columns 1–3: 1-based index of reactant nuclei (forward rate)
    //  - Columns 4–6: 1-based index of product nuclei (forward rate)
    //  - Columns 7-8: 1-based index of non (n, H1, He4) nuclei.
    //  - Column  9:   1-based index of the forward rate
    //  - Column  10:  1-based index of the reverse rate

    constexpr int H1_index = 1;
    constexpr int N_index = -1;
    constexpr int He4_index = 4;
    constexpr int NumNSERatePairs = 0;

    inline AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, NumNSERatePairs, 1, 10, amrex::Order::C> rate_pair_data {
    };
}
#endif

#endif
