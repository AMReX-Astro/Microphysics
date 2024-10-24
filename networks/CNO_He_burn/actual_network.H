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

    enum NetworkRates
    {
        k_N13_to_C13_weak_wc12 = 1,
        k_O14_to_N14_weak_wc12 = 2,
        k_O15_to_N15_weak_wc12 = 3,
        k_F17_to_O17_weak_wc12 = 4,
        k_F18_to_O18_weak_wc12 = 5,
        k_Ne18_to_F18_weak_wc12 = 6,
        k_Ne19_to_F19_weak_wc12 = 7,
        k_Mg22_to_Na22_weak_wc12 = 8,
        k_p_C12_to_N13 = 9,
        k_He4_C12_to_O16 = 10,
        k_p_C13_to_N14 = 11,
        k_p_N13_to_O14 = 12,
        k_p_N14_to_O15 = 13,
        k_He4_N14_to_F18 = 14,
        k_p_N15_to_O16 = 15,
        k_He4_N15_to_F19 = 16,
        k_He4_O14_to_Ne18 = 17,
        k_He4_O15_to_Ne19 = 18,
        k_p_O16_to_F17 = 19,
        k_He4_O16_to_Ne20 = 20,
        k_p_O17_to_F18 = 21,
        k_He4_O17_to_Ne21 = 22,
        k_p_O18_to_F19 = 23,
        k_p_F17_to_Ne18 = 24,
        k_p_F18_to_Ne19 = 25,
        k_He4_F18_to_Na22 = 26,
        k_p_F19_to_Ne20 = 27,
        k_He4_F19_to_Na23 = 28,
        k_He4_Ne18_to_Mg22 = 29,
        k_He4_Ne20_to_Mg24 = 30,
        k_p_Ne21_to_Na22 = 31,
        k_p_Na23_to_Mg24 = 32,
        k_He4_Mg24_to_Si28 = 33,
        k_p_Al27_to_Si28 = 34,
        k_He4_Si28_to_S32 = 35,
        k_p_P31_to_S32 = 36,
        k_C12_C12_to_p_Na23 = 37,
        k_C12_C12_to_He4_Ne20 = 38,
        k_He4_N13_to_p_O16 = 39,
        k_p_N15_to_He4_C12 = 40,
        k_He4_O14_to_p_F17 = 41,
        k_C12_O16_to_p_Al27 = 42,
        k_C12_O16_to_He4_Mg24 = 43,
        k_O16_O16_to_p_P31 = 44,
        k_O16_O16_to_He4_Si28 = 45,
        k_p_O17_to_He4_N14 = 46,
        k_p_O18_to_He4_N15 = 47,
        k_p_F18_to_He4_O15 = 48,
        k_He4_F18_to_p_Ne21 = 49,
        k_p_F19_to_He4_O16 = 50,
        k_He4_Ne19_to_p_Na22 = 51,
        k_p_Ne20_to_He4_F17 = 52,
        k_p_Na23_to_He4_Ne20 = 53,
        k_p_Al27_to_He4_Mg24 = 54,
        k_p_P31_to_He4_Si28 = 55,
        k_He4_He4_He4_to_C12 = 56,
        k_C12_C12_to_Mg24_modified = 57,
        k_O16_O16_to_S32_modified = 58,
        k_C12_O16_to_Si28_modified = 59,
        k_He4_S32_to_Ar36_removed = 60,
        k_p_Cl35_to_Ar36_removed = 61,
        k_p_Cl35_to_He4_S32_removed = 62,
        k_He4_Ar36_to_Ca40_removed = 63,
        k_p_K39_to_Ca40_removed = 64,
        k_p_K39_to_He4_Ar36_removed = 65,
        k_He4_Ca40_to_Ti44_removed = 66,
        k_p_Sc43_to_Ti44_removed = 67,
        k_p_Sc43_to_He4_Ca40_removed = 68,
        k_He4_Ti44_to_Cr48_removed = 69,
        k_He4_Ti44_to_p_V47_removed = 70,
        k_p_V47_to_Cr48_removed = 71,
        k_He4_Cr48_to_Fe52_removed = 72,
        k_He4_Cr48_to_p_Mn51_removed = 73,
        k_p_Mn51_to_Fe52_removed = 74,
        k_He4_Fe52_to_Ni56_removed = 75,
        k_He4_Fe52_to_p_Co55_removed = 76,
        k_p_Co55_to_Ni56_removed = 77,
        k_S32_He4_to_Ar36_approx = 78,
        k_Ar36_to_S32_He4_approx = 79,
        k_Ar36_He4_to_Ca40_approx = 80,
        k_Ca40_to_Ar36_He4_approx = 81,
        k_Ca40_He4_to_Ti44_approx = 82,
        k_Ti44_to_Ca40_He4_approx = 83,
        k_Ti44_He4_to_Cr48_approx = 84,
        k_Cr48_to_Ti44_He4_approx = 85,
        k_Cr48_He4_to_Fe52_approx = 86,
        k_Fe52_to_Cr48_He4_approx = 87,
        k_Fe52_He4_to_Ni56_approx = 88,
        k_Ni56_to_Fe52_He4_approx = 89,
        k_N13_to_p_C12_derived = 90,
        k_N14_to_p_C13_derived = 91,
        k_O14_to_p_N13_derived = 92,
        k_O15_to_p_N14_derived = 93,
        k_O16_to_p_N15_derived = 94,
        k_O16_to_He4_C12_derived = 95,
        k_F17_to_p_O16_derived = 96,
        k_F18_to_p_O17_derived = 97,
        k_F18_to_He4_N14_derived = 98,
        k_F19_to_p_O18_derived = 99,
        k_F19_to_He4_N15_derived = 100,
        k_Ne18_to_p_F17_derived = 101,
        k_Ne18_to_He4_O14_derived = 102,
        k_Ne19_to_p_F18_derived = 103,
        k_Ne19_to_He4_O15_derived = 104,
        k_Ne20_to_p_F19_derived = 105,
        k_Ne20_to_He4_O16_derived = 106,
        k_Ne21_to_He4_O17_derived = 107,
        k_Na22_to_p_Ne21_derived = 108,
        k_Na22_to_He4_F18_derived = 109,
        k_Na23_to_He4_F19_derived = 110,
        k_Mg22_to_He4_Ne18_derived = 111,
        k_Mg24_to_p_Na23_derived = 112,
        k_Mg24_to_He4_Ne20_derived = 113,
        k_Si28_to_p_Al27_derived = 114,
        k_Si28_to_He4_Mg24_derived = 115,
        k_S32_to_p_P31_derived = 116,
        k_S32_to_He4_Si28_derived = 117,
        k_C12_to_He4_He4_He4_derived = 118,
        k_He4_C12_to_p_N15_derived = 119,
        k_He4_N14_to_p_O17_derived = 120,
        k_He4_N15_to_p_O18_derived = 121,
        k_He4_O15_to_p_F18_derived = 122,
        k_p_O16_to_He4_N13_derived = 123,
        k_He4_O16_to_p_F19_derived = 124,
        k_p_F17_to_He4_O14_derived = 125,
        k_He4_F17_to_p_Ne20_derived = 126,
        k_He4_Ne20_to_p_Na23_derived = 127,
        k_p_Ne21_to_He4_F18_derived = 128,
        k_p_Na22_to_He4_Ne19_derived = 129,
        k_He4_Mg24_to_p_Al27_derived = 130,
        k_He4_Si28_to_p_P31_derived = 131,
        k_He4_S32_to_p_Cl35_derived_removed = 132,
        k_Ar36_to_He4_S32_derived_removed = 133,
        k_Ar36_to_p_Cl35_derived_removed = 134,
        k_He4_Ar36_to_p_K39_derived_removed = 135,
        k_Ca40_to_He4_Ar36_derived_removed = 136,
        k_Ca40_to_p_K39_derived_removed = 137,
        k_He4_Ca40_to_p_Sc43_derived_removed = 138,
        k_Ti44_to_He4_Ca40_derived_removed = 139,
        k_Ti44_to_p_Sc43_derived_removed = 140,
        k_Cr48_to_He4_Ti44_derived_removed = 141,
        k_Cr48_to_p_V47_derived_removed = 142,
        k_p_V47_to_He4_Ti44_derived_removed = 143,
        k_Fe52_to_He4_Cr48_derived_removed = 144,
        k_Fe52_to_p_Mn51_derived_removed = 145,
        k_p_Mn51_to_He4_Cr48_derived_removed = 146,
        k_Ni56_to_He4_Fe52_derived_removed = 147,
        k_Ni56_to_p_Co55_derived_removed = 148,
        k_p_Co55_to_He4_Fe52_derived_removed = 149,
        NumRates = k_p_Co55_to_He4_Fe52_derived_removed
    };

    // number of reaclib rates

    const int NrateReaclib = 137;

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
        "Mg22_to_Na22_weak_wc12",  // 8,
        "p_C12_to_N13",  // 9,
        "He4_C12_to_O16",  // 10,
        "p_C13_to_N14",  // 11,
        "p_N13_to_O14",  // 12,
        "p_N14_to_O15",  // 13,
        "He4_N14_to_F18",  // 14,
        "p_N15_to_O16",  // 15,
        "He4_N15_to_F19",  // 16,
        "He4_O14_to_Ne18",  // 17,
        "He4_O15_to_Ne19",  // 18,
        "p_O16_to_F17",  // 19,
        "He4_O16_to_Ne20",  // 20,
        "p_O17_to_F18",  // 21,
        "He4_O17_to_Ne21",  // 22,
        "p_O18_to_F19",  // 23,
        "p_F17_to_Ne18",  // 24,
        "p_F18_to_Ne19",  // 25,
        "He4_F18_to_Na22",  // 26,
        "p_F19_to_Ne20",  // 27,
        "He4_F19_to_Na23",  // 28,
        "He4_Ne18_to_Mg22",  // 29,
        "He4_Ne20_to_Mg24",  // 30,
        "p_Ne21_to_Na22",  // 31,
        "p_Na23_to_Mg24",  // 32,
        "He4_Mg24_to_Si28",  // 33,
        "p_Al27_to_Si28",  // 34,
        "He4_Si28_to_S32",  // 35,
        "p_P31_to_S32",  // 36,
        "C12_C12_to_p_Na23",  // 37,
        "C12_C12_to_He4_Ne20",  // 38,
        "He4_N13_to_p_O16",  // 39,
        "p_N15_to_He4_C12",  // 40,
        "He4_O14_to_p_F17",  // 41,
        "C12_O16_to_p_Al27",  // 42,
        "C12_O16_to_He4_Mg24",  // 43,
        "O16_O16_to_p_P31",  // 44,
        "O16_O16_to_He4_Si28",  // 45,
        "p_O17_to_He4_N14",  // 46,
        "p_O18_to_He4_N15",  // 47,
        "p_F18_to_He4_O15",  // 48,
        "He4_F18_to_p_Ne21",  // 49,
        "p_F19_to_He4_O16",  // 50,
        "He4_Ne19_to_p_Na22",  // 51,
        "p_Ne20_to_He4_F17",  // 52,
        "p_Na23_to_He4_Ne20",  // 53,
        "p_Al27_to_He4_Mg24",  // 54,
        "p_P31_to_He4_Si28",  // 55,
        "He4_He4_He4_to_C12",  // 56,
        "C12_C12_to_Mg24_modified",  // 57,
        "O16_O16_to_S32_modified",  // 58,
        "C12_O16_to_Si28_modified",  // 59,
        "He4_S32_to_Ar36_removed",  // 60,
        "p_Cl35_to_Ar36_removed",  // 61,
        "p_Cl35_to_He4_S32_removed",  // 62,
        "He4_Ar36_to_Ca40_removed",  // 63,
        "p_K39_to_Ca40_removed",  // 64,
        "p_K39_to_He4_Ar36_removed",  // 65,
        "He4_Ca40_to_Ti44_removed",  // 66,
        "p_Sc43_to_Ti44_removed",  // 67,
        "p_Sc43_to_He4_Ca40_removed",  // 68,
        "He4_Ti44_to_Cr48_removed",  // 69,
        "He4_Ti44_to_p_V47_removed",  // 70,
        "p_V47_to_Cr48_removed",  // 71,
        "He4_Cr48_to_Fe52_removed",  // 72,
        "He4_Cr48_to_p_Mn51_removed",  // 73,
        "p_Mn51_to_Fe52_removed",  // 74,
        "He4_Fe52_to_Ni56_removed",  // 75,
        "He4_Fe52_to_p_Co55_removed",  // 76,
        "p_Co55_to_Ni56_removed",  // 77,
        "S32_He4_to_Ar36_approx",  // 78,
        "Ar36_to_S32_He4_approx",  // 79,
        "Ar36_He4_to_Ca40_approx",  // 80,
        "Ca40_to_Ar36_He4_approx",  // 81,
        "Ca40_He4_to_Ti44_approx",  // 82,
        "Ti44_to_Ca40_He4_approx",  // 83,
        "Ti44_He4_to_Cr48_approx",  // 84,
        "Cr48_to_Ti44_He4_approx",  // 85,
        "Cr48_He4_to_Fe52_approx",  // 86,
        "Fe52_to_Cr48_He4_approx",  // 87,
        "Fe52_He4_to_Ni56_approx",  // 88,
        "Ni56_to_Fe52_He4_approx",  // 89,
        "N13_to_p_C12_derived",  // 90,
        "N14_to_p_C13_derived",  // 91,
        "O14_to_p_N13_derived",  // 92,
        "O15_to_p_N14_derived",  // 93,
        "O16_to_p_N15_derived",  // 94,
        "O16_to_He4_C12_derived",  // 95,
        "F17_to_p_O16_derived",  // 96,
        "F18_to_p_O17_derived",  // 97,
        "F18_to_He4_N14_derived",  // 98,
        "F19_to_p_O18_derived",  // 99,
        "F19_to_He4_N15_derived",  // 100,
        "Ne18_to_p_F17_derived",  // 101,
        "Ne18_to_He4_O14_derived",  // 102,
        "Ne19_to_p_F18_derived",  // 103,
        "Ne19_to_He4_O15_derived",  // 104,
        "Ne20_to_p_F19_derived",  // 105,
        "Ne20_to_He4_O16_derived",  // 106,
        "Ne21_to_He4_O17_derived",  // 107,
        "Na22_to_p_Ne21_derived",  // 108,
        "Na22_to_He4_F18_derived",  // 109,
        "Na23_to_He4_F19_derived",  // 110,
        "Mg22_to_He4_Ne18_derived",  // 111,
        "Mg24_to_p_Na23_derived",  // 112,
        "Mg24_to_He4_Ne20_derived",  // 113,
        "Si28_to_p_Al27_derived",  // 114,
        "Si28_to_He4_Mg24_derived",  // 115,
        "S32_to_p_P31_derived",  // 116,
        "S32_to_He4_Si28_derived",  // 117,
        "C12_to_He4_He4_He4_derived",  // 118,
        "He4_C12_to_p_N15_derived",  // 119,
        "He4_N14_to_p_O17_derived",  // 120,
        "He4_N15_to_p_O18_derived",  // 121,
        "He4_O15_to_p_F18_derived",  // 122,
        "p_O16_to_He4_N13_derived",  // 123,
        "He4_O16_to_p_F19_derived",  // 124,
        "p_F17_to_He4_O14_derived",  // 125,
        "He4_F17_to_p_Ne20_derived",  // 126,
        "He4_Ne20_to_p_Na23_derived",  // 127,
        "p_Ne21_to_He4_F18_derived",  // 128,
        "p_Na22_to_He4_Ne19_derived",  // 129,
        "He4_Mg24_to_p_Al27_derived",  // 130,
        "He4_Si28_to_p_P31_derived",  // 131,
        "He4_S32_to_p_Cl35_derived_removed",  // 132,
        "Ar36_to_He4_S32_derived_removed",  // 133,
        "Ar36_to_p_Cl35_derived_removed",  // 134,
        "He4_Ar36_to_p_K39_derived_removed",  // 135,
        "Ca40_to_He4_Ar36_derived_removed",  // 136,
        "Ca40_to_p_K39_derived_removed",  // 137,
        "He4_Ca40_to_p_Sc43_derived_removed",  // 138,
        "Ti44_to_He4_Ca40_derived_removed",  // 139,
        "Ti44_to_p_Sc43_derived_removed",  // 140,
        "Cr48_to_He4_Ti44_derived_removed",  // 141,
        "Cr48_to_p_V47_derived_removed",  // 142,
        "p_V47_to_He4_Ti44_derived_removed",  // 143,
        "Fe52_to_He4_Cr48_derived_removed",  // 144,
        "Fe52_to_p_Mn51_derived_removed",  // 145,
        "p_Mn51_to_He4_Cr48_derived_removed",  // 146,
        "Ni56_to_He4_Fe52_derived_removed",  // 147,
        "Ni56_to_p_Co55_derived_removed",  // 148,
        "p_Co55_to_He4_Fe52_derived_removed"  // 149,
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
