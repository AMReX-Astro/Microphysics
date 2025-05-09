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
        else if constexpr (spec == Ne21) {
            return 167.40598973999658_rt;
        }
        else if constexpr (spec == Na22) {
            return 174.14457080400098_rt;
        }
        else if constexpr (spec == Na23) {
            return 186.56435240400242_rt;
        }
        else if constexpr (spec == Mg22) {
            return 168.58082376800303_rt;
        }
        else if constexpr (spec == Mg24) {
            return 198.2570479679962_rt;
        }
        else if constexpr (spec == Al27) {
            return 224.95193723199915_rt;
        }
        else if constexpr (spec == Si28) {
            return 236.53684539599638_rt;
        }
        else if constexpr (spec == P31) {
            return 262.9161999600037_rt;
        }
        else if constexpr (spec == S32) {
            return 271.78016372399725_rt;
        }
        else if constexpr (spec == Ar36) {
            return 306.7167469519991_rt;
        }
        else if constexpr (spec == Ca40) {
            return 342.05218528000114_rt;
        }
        else if constexpr (spec == Ti44) {
            return 375.47496160800074_rt;
        }
        else if constexpr (spec == Cr48) {
            return 411.4679399359957_rt;
        }
        else if constexpr (spec == Fe52) {
            return 447.6996182639923_rt;
        }
        else if constexpr (spec == Ni56) {
            return 483.9956965919919_rt;
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
        k_Ne21_to_He4_O17 = 25,
        k_Na22_to_p_Ne21 = 26,
        k_Na22_to_He4_F18 = 27,
        k_Na23_to_He4_F19 = 28,
        k_Mg22_to_He4_Ne18 = 29,
        k_Mg24_to_p_Na23 = 30,
        k_Mg24_to_He4_Ne20 = 31,
        k_Si28_to_p_Al27 = 32,
        k_Si28_to_He4_Mg24 = 33,
        k_S32_to_p_P31 = 34,
        k_S32_to_He4_Si28 = 35,
        k_C12_to_He4_He4_He4 = 36,
        k_p_C12_to_N13 = 37,
        k_He4_C12_to_O16 = 38,
        k_p_C13_to_N14 = 39,
        k_p_N13_to_O14 = 40,
        k_p_N14_to_O15 = 41,
        k_He4_N14_to_F18 = 42,
        k_p_N15_to_O16 = 43,
        k_He4_N15_to_F19 = 44,
        k_He4_O14_to_Ne18 = 45,
        k_He4_O15_to_Ne19 = 46,
        k_p_O16_to_F17 = 47,
        k_He4_O16_to_Ne20 = 48,
        k_p_O17_to_F18 = 49,
        k_He4_O17_to_Ne21 = 50,
        k_p_O18_to_F19 = 51,
        k_p_F17_to_Ne18 = 52,
        k_p_F18_to_Ne19 = 53,
        k_He4_F18_to_Na22 = 54,
        k_p_F19_to_Ne20 = 55,
        k_He4_F19_to_Na23 = 56,
        k_He4_Ne18_to_Mg22 = 57,
        k_He4_Ne20_to_Mg24 = 58,
        k_p_Ne21_to_Na22 = 59,
        k_p_Na23_to_Mg24 = 60,
        k_He4_Mg24_to_Si28 = 61,
        k_p_Al27_to_Si28 = 62,
        k_He4_Si28_to_S32 = 63,
        k_p_P31_to_S32 = 64,
        k_He4_C12_to_p_N15 = 65,
        k_C12_C12_to_p_Na23 = 66,
        k_C12_C12_to_He4_Ne20 = 67,
        k_He4_N13_to_p_O16 = 68,
        k_He4_N14_to_p_O17 = 69,
        k_p_N15_to_He4_C12 = 70,
        k_He4_N15_to_p_O18 = 71,
        k_He4_O14_to_p_F17 = 72,
        k_He4_O15_to_p_F18 = 73,
        k_p_O16_to_He4_N13 = 74,
        k_He4_O16_to_p_F19 = 75,
        k_C12_O16_to_p_Al27 = 76,
        k_C12_O16_to_He4_Mg24 = 77,
        k_O16_O16_to_p_P31 = 78,
        k_O16_O16_to_He4_Si28 = 79,
        k_p_O17_to_He4_N14 = 80,
        k_p_O18_to_He4_N15 = 81,
        k_p_F17_to_He4_O14 = 82,
        k_He4_F17_to_p_Ne20 = 83,
        k_p_F18_to_He4_O15 = 84,
        k_He4_F18_to_p_Ne21 = 85,
        k_p_F19_to_He4_O16 = 86,
        k_He4_Ne19_to_p_Na22 = 87,
        k_p_Ne20_to_He4_F17 = 88,
        k_He4_Ne20_to_p_Na23 = 89,
        k_p_Ne21_to_He4_F18 = 90,
        k_p_Na22_to_He4_Ne19 = 91,
        k_p_Na23_to_He4_Ne20 = 92,
        k_He4_Mg24_to_p_Al27 = 93,
        k_p_Al27_to_He4_Mg24 = 94,
        k_He4_Si28_to_p_P31 = 95,
        k_p_P31_to_He4_Si28 = 96,
        k_He4_He4_He4_to_C12 = 97,
        k_C12_C12_to_Mg24_modified = 98,
        k_O16_O16_to_S32_modified = 99,
        k_C12_O16_to_Si28_modified = 100,
        k_He4_S32_to_Ar36_removed = 101,
        k_He4_S32_to_p_Cl35_removed = 102,
        k_p_Cl35_to_Ar36_removed = 103,
        k_Ar36_to_He4_S32_removed = 104,
        k_Ar36_to_p_Cl35_removed = 105,
        k_p_Cl35_to_He4_S32_removed = 106,
        k_He4_Ar36_to_Ca40_removed = 107,
        k_He4_Ar36_to_p_K39_removed = 108,
        k_p_K39_to_Ca40_removed = 109,
        k_Ca40_to_He4_Ar36_removed = 110,
        k_Ca40_to_p_K39_removed = 111,
        k_p_K39_to_He4_Ar36_removed = 112,
        k_He4_Ca40_to_Ti44_removed = 113,
        k_He4_Ca40_to_p_Sc43_removed = 114,
        k_p_Sc43_to_Ti44_removed = 115,
        k_Ti44_to_He4_Ca40_removed = 116,
        k_Ti44_to_p_Sc43_removed = 117,
        k_p_Sc43_to_He4_Ca40_removed = 118,
        k_He4_Ti44_to_Cr48_removed = 119,
        k_He4_Ti44_to_p_V47_removed = 120,
        k_p_V47_to_Cr48_removed = 121,
        k_Cr48_to_He4_Ti44_removed = 122,
        k_Cr48_to_p_V47_removed = 123,
        k_p_V47_to_He4_Ti44_removed = 124,
        k_He4_Cr48_to_Fe52_removed = 125,
        k_He4_Cr48_to_p_Mn51_removed = 126,
        k_p_Mn51_to_Fe52_removed = 127,
        k_Fe52_to_He4_Cr48_removed = 128,
        k_Fe52_to_p_Mn51_removed = 129,
        k_p_Mn51_to_He4_Cr48_removed = 130,
        k_He4_Fe52_to_Ni56_removed = 131,
        k_He4_Fe52_to_p_Co55_removed = 132,
        k_p_Co55_to_Ni56_removed = 133,
        k_Ni56_to_He4_Fe52_removed = 134,
        k_Ni56_to_p_Co55_removed = 135,
        k_p_Co55_to_He4_Fe52_removed = 136,
        k_Mg22_to_Na22 = 137,
        k_Na22_to_Mg22 = 138,
        k_S32_He4_to_Ar36_approx = 139,
        k_Ar36_to_S32_He4_approx = 140,
        k_Ar36_He4_to_Ca40_approx = 141,
        k_Ca40_to_Ar36_He4_approx = 142,
        k_Ca40_He4_to_Ti44_approx = 143,
        k_Ti44_to_Ca40_He4_approx = 144,
        k_Ti44_He4_to_Cr48_approx = 145,
        k_Cr48_to_Ti44_He4_approx = 146,
        k_Cr48_He4_to_Fe52_approx = 147,
        k_Fe52_to_Cr48_He4_approx = 148,
        k_Fe52_He4_to_Ni56_approx = 149,
        k_Ni56_to_Fe52_He4_approx = 150,
        NumRates = k_Ni56_to_Fe52_He4_approx
    };

    // number of reaclib rates

    const int NrateReaclib = 136;

    // number of tabular rates

    const int NrateTabular = 2;

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
        "Ne21_to_He4_O17",  // 25,
        "Na22_to_p_Ne21",  // 26,
        "Na22_to_He4_F18",  // 27,
        "Na23_to_He4_F19",  // 28,
        "Mg22_to_He4_Ne18",  // 29,
        "Mg24_to_p_Na23",  // 30,
        "Mg24_to_He4_Ne20",  // 31,
        "Si28_to_p_Al27",  // 32,
        "Si28_to_He4_Mg24",  // 33,
        "S32_to_p_P31",  // 34,
        "S32_to_He4_Si28",  // 35,
        "C12_to_He4_He4_He4",  // 36,
        "p_C12_to_N13",  // 37,
        "He4_C12_to_O16",  // 38,
        "p_C13_to_N14",  // 39,
        "p_N13_to_O14",  // 40,
        "p_N14_to_O15",  // 41,
        "He4_N14_to_F18",  // 42,
        "p_N15_to_O16",  // 43,
        "He4_N15_to_F19",  // 44,
        "He4_O14_to_Ne18",  // 45,
        "He4_O15_to_Ne19",  // 46,
        "p_O16_to_F17",  // 47,
        "He4_O16_to_Ne20",  // 48,
        "p_O17_to_F18",  // 49,
        "He4_O17_to_Ne21",  // 50,
        "p_O18_to_F19",  // 51,
        "p_F17_to_Ne18",  // 52,
        "p_F18_to_Ne19",  // 53,
        "He4_F18_to_Na22",  // 54,
        "p_F19_to_Ne20",  // 55,
        "He4_F19_to_Na23",  // 56,
        "He4_Ne18_to_Mg22",  // 57,
        "He4_Ne20_to_Mg24",  // 58,
        "p_Ne21_to_Na22",  // 59,
        "p_Na23_to_Mg24",  // 60,
        "He4_Mg24_to_Si28",  // 61,
        "p_Al27_to_Si28",  // 62,
        "He4_Si28_to_S32",  // 63,
        "p_P31_to_S32",  // 64,
        "He4_C12_to_p_N15",  // 65,
        "C12_C12_to_p_Na23",  // 66,
        "C12_C12_to_He4_Ne20",  // 67,
        "He4_N13_to_p_O16",  // 68,
        "He4_N14_to_p_O17",  // 69,
        "p_N15_to_He4_C12",  // 70,
        "He4_N15_to_p_O18",  // 71,
        "He4_O14_to_p_F17",  // 72,
        "He4_O15_to_p_F18",  // 73,
        "p_O16_to_He4_N13",  // 74,
        "He4_O16_to_p_F19",  // 75,
        "C12_O16_to_p_Al27",  // 76,
        "C12_O16_to_He4_Mg24",  // 77,
        "O16_O16_to_p_P31",  // 78,
        "O16_O16_to_He4_Si28",  // 79,
        "p_O17_to_He4_N14",  // 80,
        "p_O18_to_He4_N15",  // 81,
        "p_F17_to_He4_O14",  // 82,
        "He4_F17_to_p_Ne20",  // 83,
        "p_F18_to_He4_O15",  // 84,
        "He4_F18_to_p_Ne21",  // 85,
        "p_F19_to_He4_O16",  // 86,
        "He4_Ne19_to_p_Na22",  // 87,
        "p_Ne20_to_He4_F17",  // 88,
        "He4_Ne20_to_p_Na23",  // 89,
        "p_Ne21_to_He4_F18",  // 90,
        "p_Na22_to_He4_Ne19",  // 91,
        "p_Na23_to_He4_Ne20",  // 92,
        "He4_Mg24_to_p_Al27",  // 93,
        "p_Al27_to_He4_Mg24",  // 94,
        "He4_Si28_to_p_P31",  // 95,
        "p_P31_to_He4_Si28",  // 96,
        "He4_He4_He4_to_C12",  // 97,
        "C12_C12_to_Mg24_modified",  // 98,
        "O16_O16_to_S32_modified",  // 99,
        "C12_O16_to_Si28_modified",  // 100,
        "He4_S32_to_Ar36_removed",  // 101,
        "He4_S32_to_p_Cl35_removed",  // 102,
        "p_Cl35_to_Ar36_removed",  // 103,
        "Ar36_to_He4_S32_removed",  // 104,
        "Ar36_to_p_Cl35_removed",  // 105,
        "p_Cl35_to_He4_S32_removed",  // 106,
        "He4_Ar36_to_Ca40_removed",  // 107,
        "He4_Ar36_to_p_K39_removed",  // 108,
        "p_K39_to_Ca40_removed",  // 109,
        "Ca40_to_He4_Ar36_removed",  // 110,
        "Ca40_to_p_K39_removed",  // 111,
        "p_K39_to_He4_Ar36_removed",  // 112,
        "He4_Ca40_to_Ti44_removed",  // 113,
        "He4_Ca40_to_p_Sc43_removed",  // 114,
        "p_Sc43_to_Ti44_removed",  // 115,
        "Ti44_to_He4_Ca40_removed",  // 116,
        "Ti44_to_p_Sc43_removed",  // 117,
        "p_Sc43_to_He4_Ca40_removed",  // 118,
        "He4_Ti44_to_Cr48_removed",  // 119,
        "He4_Ti44_to_p_V47_removed",  // 120,
        "p_V47_to_Cr48_removed",  // 121,
        "Cr48_to_He4_Ti44_removed",  // 122,
        "Cr48_to_p_V47_removed",  // 123,
        "p_V47_to_He4_Ti44_removed",  // 124,
        "He4_Cr48_to_Fe52_removed",  // 125,
        "He4_Cr48_to_p_Mn51_removed",  // 126,
        "p_Mn51_to_Fe52_removed",  // 127,
        "Fe52_to_He4_Cr48_removed",  // 128,
        "Fe52_to_p_Mn51_removed",  // 129,
        "p_Mn51_to_He4_Cr48_removed",  // 130,
        "He4_Fe52_to_Ni56_removed",  // 131,
        "He4_Fe52_to_p_Co55_removed",  // 132,
        "p_Co55_to_Ni56_removed",  // 133,
        "Ni56_to_He4_Fe52_removed",  // 134,
        "Ni56_to_p_Co55_removed",  // 135,
        "p_Co55_to_He4_Fe52_removed",  // 136,
        "Mg22_to_Na22",  // 137,
        "Na22_to_Mg22",  // 138,
        "S32_He4_to_Ar36_approx",  // 139,
        "Ar36_to_S32_He4_approx",  // 140,
        "Ar36_He4_to_Ca40_approx",  // 141,
        "Ca40_to_Ar36_He4_approx",  // 142,
        "Ca40_He4_to_Ti44_approx",  // 143,
        "Ti44_to_Ca40_He4_approx",  // 144,
        "Ti44_He4_to_Cr48_approx",  // 145,
        "Cr48_to_Ti44_He4_approx",  // 146,
        "Cr48_He4_to_Fe52_approx",  // 147,
        "Fe52_to_Cr48_He4_approx",  // 148,
        "Fe52_He4_to_Ni56_approx",  // 149,
        "Ni56_to_Fe52_He4_approx"  // 150,
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

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
