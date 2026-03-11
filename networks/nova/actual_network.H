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
        k_B8_to_He4_He4_reaclib = 4,
        k_p_p_to_d_reaclib_beta_pos = 5,
        k_p_p_to_d_reaclib_electron_capture = 6,
        k_p_d_to_He3_reaclib = 7,
        k_d_d_to_He4_reaclib = 8,
        k_p_He3_to_He4_reaclib = 9,
        k_He4_He3_to_Be7_reaclib = 10,
        k_p_Be7_to_B8_reaclib = 11,
        k_p_C12_to_N13_reaclib = 12,
        k_He4_C12_to_O16_reaclib = 13,
        k_p_C13_to_N14_reaclib = 14,
        k_p_N13_to_O14_reaclib = 15,
        k_p_N14_to_O15_reaclib = 16,
        k_He4_N14_to_F18_reaclib = 17,
        k_p_N15_to_O16_reaclib = 18,
        k_p_O16_to_F17_reaclib = 19,
        k_p_O17_to_F18_reaclib = 20,
        k_d_He3_to_p_He4_reaclib = 21,
        k_He4_N13_to_p_O16_reaclib = 22,
        k_p_N15_to_He4_C12_reaclib = 23,
        k_He4_O14_to_p_F17_reaclib = 24,
        k_p_O17_to_He4_N14_reaclib = 25,
        k_p_F18_to_He4_O15_reaclib = 26,
        k_He3_He3_to_p_p_He4_reaclib = 27,
        k_d_Be7_to_p_He4_He4_reaclib = 28,
        k_He3_Be7_to_p_p_He4_He4_reaclib = 29,
        k_He4_He4_He4_to_C12_reaclib = 30,
        k_F17_to_O17_weaktab = 31,
        k_O17_to_F17_weaktab = 32,
        k_He3_to_p_d_derived = 33,
        k_He4_to_d_d_derived = 34,
        k_Be7_to_He4_He3_derived = 35,
        k_B8_to_p_Be7_derived = 36,
        k_N13_to_p_C12_derived = 37,
        k_N14_to_p_C13_derived = 38,
        k_O14_to_p_N13_derived = 39,
        k_O15_to_p_N14_derived = 40,
        k_O16_to_p_N15_derived = 41,
        k_O16_to_He4_C12_derived = 42,
        k_F17_to_p_O16_derived = 43,
        k_F18_to_p_O17_derived = 44,
        k_F18_to_He4_N14_derived = 45,
        k_C12_to_He4_He4_He4_derived = 46,
        k_p_He4_to_d_He3_derived = 47,
        k_He4_C12_to_p_N15_derived = 48,
        k_He4_N14_to_p_O17_derived = 49,
        k_He4_O15_to_p_F18_derived = 50,
        k_p_O16_to_He4_N13_derived = 51,
        k_p_F17_to_He4_O14_derived = 52,
        k_p_p_He4_to_He3_He3_derived = 53,
        k_p_He4_He4_to_d_Be7_derived = 54,
        k_p_p_He4_He4_to_He3_Be7_derived = 55,
        NumRates = k_p_p_He4_He4_to_He3_Be7_derived
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "N13_to_C13_reaclib",  // 1,
        "O14_to_N14_reaclib",  // 2,
        "O15_to_N15_reaclib",  // 3,
        "B8_to_He4_He4_reaclib",  // 4,
        "p_p_to_d_reaclib_beta_pos",  // 5,
        "p_p_to_d_reaclib_electron_capture",  // 6,
        "p_d_to_He3_reaclib",  // 7,
        "d_d_to_He4_reaclib",  // 8,
        "p_He3_to_He4_reaclib",  // 9,
        "He4_He3_to_Be7_reaclib",  // 10,
        "p_Be7_to_B8_reaclib",  // 11,
        "p_C12_to_N13_reaclib",  // 12,
        "He4_C12_to_O16_reaclib",  // 13,
        "p_C13_to_N14_reaclib",  // 14,
        "p_N13_to_O14_reaclib",  // 15,
        "p_N14_to_O15_reaclib",  // 16,
        "He4_N14_to_F18_reaclib",  // 17,
        "p_N15_to_O16_reaclib",  // 18,
        "p_O16_to_F17_reaclib",  // 19,
        "p_O17_to_F18_reaclib",  // 20,
        "d_He3_to_p_He4_reaclib",  // 21,
        "He4_N13_to_p_O16_reaclib",  // 22,
        "p_N15_to_He4_C12_reaclib",  // 23,
        "He4_O14_to_p_F17_reaclib",  // 24,
        "p_O17_to_He4_N14_reaclib",  // 25,
        "p_F18_to_He4_O15_reaclib",  // 26,
        "He3_He3_to_p_p_He4_reaclib",  // 27,
        "d_Be7_to_p_He4_He4_reaclib",  // 28,
        "He3_Be7_to_p_p_He4_He4_reaclib",  // 29,
        "He4_He4_He4_to_C12_reaclib",  // 30,
        "F17_to_O17_weaktab",  // 31,
        "O17_to_F17_weaktab",  // 32,
        "He3_to_p_d_derived",  // 33,
        "He4_to_d_d_derived",  // 34,
        "Be7_to_He4_He3_derived",  // 35,
        "B8_to_p_Be7_derived",  // 36,
        "N13_to_p_C12_derived",  // 37,
        "N14_to_p_C13_derived",  // 38,
        "O14_to_p_N13_derived",  // 39,
        "O15_to_p_N14_derived",  // 40,
        "O16_to_p_N15_derived",  // 41,
        "O16_to_He4_C12_derived",  // 42,
        "F17_to_p_O16_derived",  // 43,
        "F18_to_p_O17_derived",  // 44,
        "F18_to_He4_N14_derived",  // 45,
        "C12_to_He4_He4_He4_derived",  // 46,
        "p_He4_to_d_He3_derived",  // 47,
        "He4_C12_to_p_N15_derived",  // 48,
        "He4_N14_to_p_O17_derived",  // 49,
        "He4_O15_to_p_F18_derived",  // 50,
        "p_O16_to_He4_N13_derived",  // 51,
        "p_F17_to_He4_O14_derived",  // 52,
        "p_p_He4_to_He3_He3_derived",  // 53,
        "p_He4_He4_to_d_Be7_derived",  // 54,
        "p_p_He4_He4_to_He3_Be7_derived"  // 55,
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
    constexpr int He4_index = 4;
    constexpr int NumNSERatePairs = 20;

    inline AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, NumNSERatePairs, 1, 10, amrex::Order::C> rate_pair_data {
        H1, H2, -1, He3, -1, -1, H2, He3, k_p_d_to_He3_reaclib, k_He3_to_p_d_derived,
        H2, H2, -1, He4, -1, -1, H2, H2, k_d_d_to_He4_reaclib, k_He4_to_d_d_derived,
        He4, He3, -1, Be7, -1, -1, He3, Be7, k_He4_He3_to_Be7_reaclib, k_Be7_to_He4_He3_derived,
        H1, Be7, -1, B8, -1, -1, Be7, B8, k_p_Be7_to_B8_reaclib, k_B8_to_p_Be7_derived,
        H1, C12, -1, N13, -1, -1, C12, N13, k_p_C12_to_N13_reaclib, k_N13_to_p_C12_derived,
        He4, C12, -1, O16, -1, -1, C12, O16, k_He4_C12_to_O16_reaclib, k_O16_to_He4_C12_derived,
        H1, C13, -1, N14, -1, -1, C13, N14, k_p_C13_to_N14_reaclib, k_N14_to_p_C13_derived,
        H1, N13, -1, O14, -1, -1, N13, O14, k_p_N13_to_O14_reaclib, k_O14_to_p_N13_derived,
        H1, N14, -1, O15, -1, -1, N14, O15, k_p_N14_to_O15_reaclib, k_O15_to_p_N14_derived,
        He4, N14, -1, F18, -1, -1, N14, F18, k_He4_N14_to_F18_reaclib, k_F18_to_He4_N14_derived,
        H1, N15, -1, O16, -1, -1, N15, O16, k_p_N15_to_O16_reaclib, k_O16_to_p_N15_derived,
        H1, O16, -1, F17, -1, -1, O16, F17, k_p_O16_to_F17_reaclib, k_F17_to_p_O16_derived,
        H1, O17, -1, F18, -1, -1, O17, F18, k_p_O17_to_F18_reaclib, k_F18_to_p_O17_derived,
        H2, He3, -1, H1, He4, -1, H2, He3, k_d_He3_to_p_He4_reaclib, k_p_He4_to_d_He3_derived,
        He4, N13, -1, H1, O16, -1, N13, O16, k_He4_N13_to_p_O16_reaclib, k_p_O16_to_He4_N13_derived,
        H1, N15, -1, He4, C12, -1, N15, C12, k_p_N15_to_He4_C12_reaclib, k_He4_C12_to_p_N15_derived,
        He4, O14, -1, H1, F17, -1, O14, F17, k_He4_O14_to_p_F17_reaclib, k_p_F17_to_He4_O14_derived,
        H1, O17, -1, He4, N14, -1, O17, N14, k_p_O17_to_He4_N14_reaclib, k_He4_N14_to_p_O17_derived,
        H1, F18, -1, He4, O15, -1, F18, O15, k_p_F18_to_He4_O15_reaclib, k_He4_O15_to_p_F18_derived,
        He4, He4, He4, C12, -1, -1, C12, -1, k_He4_He4_He4_to_C12_reaclib, k_C12_to_He4_He4_He4_derived
    };
}
#endif

#endif
