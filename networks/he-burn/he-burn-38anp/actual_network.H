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

        if constexpr (spec == N) {
            return 0.0_rt;
        }
        else if constexpr (spec == H1) {
            return 0.0_rt;
        }
        else if constexpr (spec == P_nse) {
            return 0.0_rt;
        }
        else if constexpr (spec == He4) {
            return 28.295662457999697_rt;
        }
        else if constexpr (spec == C12) {
            return 92.16173498399803_rt;
        }
        else if constexpr (spec == N13) {
            return 94.10522604799917_rt;
        }
        else if constexpr (spec == N14) {
            return 104.65860734799753_rt;
        }
        else if constexpr (spec == O16) {
            return 127.6193154119992_rt;
        }
        else if constexpr (spec == O17) {
            return 131.76239561199873_rt;
        }
        else if constexpr (spec == F18) {
            return 137.36950247599816_rt;
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
        else if constexpr (spec == Mg24) {
            return 198.2570479679962_rt;
        }
        else if constexpr (spec == Mg25) {
            return 205.5875680680001_rt;
        }
        else if constexpr (spec == Al27) {
            return 224.95193723199915_rt;
        }
        else if constexpr (spec == Si28) {
            return 236.53684539599638_rt;
        }
        else if constexpr (spec == Si29) {
            return 245.01044789599834_rt;
        }
        else if constexpr (spec == P30) {
            return 250.60519745999775_rt;
        }
        else if constexpr (spec == P31) {
            return 262.9161999600037_rt;
        }
        else if constexpr (spec == S31) {
            return 256.7358385240004_rt;
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
        else if constexpr (spec == Mn51) {
            return 440.321747199996_rt;
        }
        else if constexpr (spec == Fe52) {
            return 447.6996182639923_rt;
        }
        else if constexpr (spec == Fe54) {
            return 471.76475446399854_rt;
        }
        else if constexpr (spec == Fe56) {
            return 492.2599506639962_rt;
        }
        else if constexpr (spec == Co55) {
            return 476.82912552799826_rt;
        }
        else if constexpr (spec == Co56) {
            return 486.91094362799777_rt;
        }
        else if constexpr (spec == Co57) {
            return 498.28746172798856_rt;
        }
        else if constexpr (spec == Ni56) {
            return 483.9956965919919_rt;
        }
        else if constexpr (spec == Ni58) {
            return 506.4596327920008_rt;
        }
        else if constexpr (spec == Cu59) {
            return 509.8782038560021_rt;
        }
        else if constexpr (spec == Zn60) {
            return 514.9831749199948_rt;
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
        k_S31_to_P31_weak_wc12 = 1,
        k_p_C12_to_N13 = 2,
        k_He4_C12_to_O16 = 3,
        k_He4_N14_to_F18 = 4,
        k_He4_O16_to_Ne20 = 5,
        k_p_O17_to_F18 = 6,
        k_He4_O17_to_Ne21 = 7,
        k_He4_F18_to_Na22 = 8,
        k_He4_Ne20_to_Mg24 = 9,
        k_p_Ne21_to_Na22 = 10,
        k_He4_Ne21_to_Mg25 = 11,
        k_p_Na23_to_Mg24 = 12,
        k_He4_Na23_to_Al27 = 13,
        k_He4_Mg24_to_Si28 = 14,
        k_He4_Mg25_to_Si29 = 15,
        k_p_Al27_to_Si28 = 16,
        k_He4_Al27_to_P31 = 17,
        k_He4_Si28_to_S32 = 18,
        k_p_Si29_to_P30 = 19,
        k_p_P30_to_S31 = 20,
        k_p_P31_to_S32 = 21,
        k_He4_Cr48_to_Fe52 = 22,
        k_p_nse_Mn51_to_Fe52 = 23,
        k_He4_Mn51_to_Co55 = 24,
        k_He4_Fe52_to_Ni56 = 25,
        k_p_nse_Co55_to_Ni56 = 26,
        k_He4_Co55_to_Cu59 = 27,
        k_He4_Ni56_to_Zn60 = 28,
        k_p_nse_Cu59_to_Zn60 = 29,
        k_C12_C12_to_p_Na23 = 30,
        k_C12_C12_to_He4_Ne20 = 31,
        k_He4_N13_to_p_O16 = 32,
        k_C12_O16_to_p_Al27 = 33,
        k_C12_O16_to_He4_Mg24 = 34,
        k_O16_O16_to_p_P31 = 35,
        k_O16_O16_to_He4_Si28 = 36,
        k_p_O17_to_He4_N14 = 37,
        k_He4_F18_to_p_Ne21 = 38,
        k_He4_Na22_to_p_Mg25 = 39,
        k_p_Na23_to_He4_Ne20 = 40,
        k_p_Al27_to_He4_Mg24 = 41,
        k_p_P31_to_He4_Si28 = 42,
        k_He4_Cr48_to_p_nse_Mn51 = 43,
        k_He4_Fe52_to_p_nse_Co55 = 44,
        k_p_nse_Cu59_to_He4_Ni56 = 45,
        k_He4_He4_He4_to_C12 = 46,
        k_C12_C12_to_Mg24_modified = 47,
        k_O16_O16_to_S32_modified = 48,
        k_C12_O16_to_Si28_modified = 49,
        k_p_nse_Fe54_to_Co55 = 50,
        k_He4_Fe54_to_Ni58 = 51,
        k_p_nse_Fe56_to_Co57 = 52,
        k_n_Co55_to_Co56 = 53,
        k_n_Co56_to_Co57 = 54,
        k_p_nse_Co57_to_Ni58 = 55,
        k_p_nse_Ni58_to_Cu59 = 56,
        k_He4_Mn51_to_p_nse_Fe54 = 57,
        k_He4_Co55_to_p_nse_Ni58 = 58,
        k_n_Co56_to_p_nse_Fe56 = 59,
        k_p_nse_Co57_to_He4_Fe54 = 60,
        k_n_Ni56_to_p_nse_Co56 = 61,
        k_n_Cu59_to_He4_Co56 = 62,
        k_He4_S32_to_Ar36_removed = 63,
        k_p_Cl35_to_Ar36_removed = 64,
        k_p_Cl35_to_He4_S32_removed = 65,
        k_He4_Ar36_to_Ca40_removed = 66,
        k_p_K39_to_Ca40_removed = 67,
        k_p_K39_to_He4_Ar36_removed = 68,
        k_He4_Ca40_to_Ti44_removed = 69,
        k_p_Sc43_to_Ti44_removed = 70,
        k_p_Sc43_to_He4_Ca40_removed = 71,
        k_He4_Ti44_to_Cr48_removed = 72,
        k_He4_Ti44_to_p_V47_removed = 73,
        k_p_V47_to_Cr48_removed = 74,
        k_n_Fe52_to_Fe53_removed = 75,
        k_n_Fe53_to_Fe54_removed = 76,
        k_n_Fe54_to_Fe55_removed = 77,
        k_n_Fe55_to_Fe56_removed = 78,
        k_n_Ni56_to_Ni57_removed = 79,
        k_n_Ni57_to_Ni58_removed = 80,
        k_Co56_to_Fe56 = 81,
        k_Co56_to_Ni56 = 82,
        k_Fe56_to_Co56 = 83,
        k_n_to_p = 84,
        k_Ni56_to_Co56 = 85,
        k_p_to_n = 86,
        k_S32_He4_to_Ar36_approx = 87,
        k_Ar36_to_S32_He4_approx = 88,
        k_Ar36_He4_to_Ca40_approx = 89,
        k_Ca40_to_Ar36_He4_approx = 90,
        k_Ca40_He4_to_Ti44_approx = 91,
        k_Ti44_to_Ca40_He4_approx = 92,
        k_Ti44_He4_to_Cr48_approx = 93,
        k_Cr48_to_Ti44_He4_approx = 94,
        k_Fe52_n_n_to_Fe54_approx = 95,
        k_Fe54_to_Fe52_n_n_approx = 96,
        k_Fe54_n_n_to_Fe56_approx = 97,
        k_Fe56_to_Fe54_n_n_approx = 98,
        k_Ni56_n_n_to_Ni58_approx = 99,
        k_Ni58_to_Ni56_n_n_approx = 100,
        k_N13_to_p_C12_derived = 101,
        k_O16_to_He4_C12_derived = 102,
        k_F18_to_p_O17_derived = 103,
        k_F18_to_He4_N14_derived = 104,
        k_Ne20_to_He4_O16_derived = 105,
        k_Ne21_to_He4_O17_derived = 106,
        k_Na22_to_p_Ne21_derived = 107,
        k_Na22_to_He4_F18_derived = 108,
        k_Mg24_to_p_Na23_derived = 109,
        k_Mg24_to_He4_Ne20_derived = 110,
        k_Mg25_to_He4_Ne21_derived = 111,
        k_Al27_to_He4_Na23_derived = 112,
        k_Si28_to_p_Al27_derived = 113,
        k_Si28_to_He4_Mg24_derived = 114,
        k_Si29_to_He4_Mg25_derived = 115,
        k_P30_to_p_Si29_derived = 116,
        k_P31_to_He4_Al27_derived = 117,
        k_S31_to_p_P30_derived = 118,
        k_S32_to_p_P31_derived = 119,
        k_S32_to_He4_Si28_derived = 120,
        k_Fe52_to_p_nse_Mn51_derived = 121,
        k_Fe52_to_He4_Cr48_derived = 122,
        k_Co55_to_He4_Mn51_derived = 123,
        k_Ni56_to_p_nse_Co55_derived = 124,
        k_Ni56_to_He4_Fe52_derived = 125,
        k_Cu59_to_He4_Co55_derived = 126,
        k_Zn60_to_p_nse_Cu59_derived = 127,
        k_Zn60_to_He4_Ni56_derived = 128,
        k_C12_to_He4_He4_He4_derived = 129,
        k_He4_N14_to_p_O17_derived = 130,
        k_p_O16_to_He4_N13_derived = 131,
        k_He4_Ne20_to_p_Na23_derived = 132,
        k_p_Ne21_to_He4_F18_derived = 133,
        k_He4_Mg24_to_p_Al27_derived = 134,
        k_p_Mg25_to_He4_Na22_derived = 135,
        k_He4_Si28_to_p_P31_derived = 136,
        k_p_nse_Mn51_to_He4_Cr48_derived = 137,
        k_p_nse_Co55_to_He4_Fe52_derived = 138,
        k_He4_Ni56_to_p_nse_Cu59_derived = 139,
        k_Co55_to_p_nse_Fe54_derived = 140,
        k_Co56_to_n_Co55_derived = 141,
        k_Co57_to_n_Co56_derived = 142,
        k_Co57_to_p_nse_Fe56_derived = 143,
        k_Ni58_to_p_nse_Co57_derived = 144,
        k_Ni58_to_He4_Fe54_derived = 145,
        k_Cu59_to_p_nse_Ni58_derived = 146,
        k_p_nse_Fe54_to_He4_Mn51_derived = 147,
        k_He4_Fe54_to_p_nse_Co57_derived = 148,
        k_p_nse_Fe56_to_n_Co56_derived = 149,
        k_p_nse_Co56_to_n_Ni56_derived = 150,
        k_He4_Co56_to_n_Cu59_derived = 151,
        k_p_nse_Ni58_to_He4_Co55_derived = 152,
        k_He4_S32_to_p_Cl35_derived_removed = 153,
        k_Ar36_to_He4_S32_derived_removed = 154,
        k_Ar36_to_p_Cl35_derived_removed = 155,
        k_He4_Ar36_to_p_K39_derived_removed = 156,
        k_Ca40_to_He4_Ar36_derived_removed = 157,
        k_Ca40_to_p_K39_derived_removed = 158,
        k_He4_Ca40_to_p_Sc43_derived_removed = 159,
        k_Ti44_to_He4_Ca40_derived_removed = 160,
        k_Ti44_to_p_Sc43_derived_removed = 161,
        k_Cr48_to_He4_Ti44_derived_removed = 162,
        k_Cr48_to_p_V47_derived_removed = 163,
        k_p_V47_to_He4_Ti44_derived_removed = 164,
        k_Fe54_to_n_Fe53_derived_removed = 165,
        k_Fe53_to_n_Fe52_derived_removed = 166,
        k_Fe56_to_n_Fe55_derived_removed = 167,
        k_Fe55_to_n_Fe54_derived_removed = 168,
        k_Ni58_to_n_Ni57_derived_removed = 169,
        k_Ni57_to_n_Ni56_derived_removed = 170,
        NumRates = k_Ni57_to_n_Ni56_derived_removed
    };

    // number of reaclib rates

    const int NrateReaclib = 150;

    // number of tabular rates

    const int NrateTabular = 6;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "S31_to_P31_weak_wc12",  // 1,
        "p_C12_to_N13",  // 2,
        "He4_C12_to_O16",  // 3,
        "He4_N14_to_F18",  // 4,
        "He4_O16_to_Ne20",  // 5,
        "p_O17_to_F18",  // 6,
        "He4_O17_to_Ne21",  // 7,
        "He4_F18_to_Na22",  // 8,
        "He4_Ne20_to_Mg24",  // 9,
        "p_Ne21_to_Na22",  // 10,
        "He4_Ne21_to_Mg25",  // 11,
        "p_Na23_to_Mg24",  // 12,
        "He4_Na23_to_Al27",  // 13,
        "He4_Mg24_to_Si28",  // 14,
        "He4_Mg25_to_Si29",  // 15,
        "p_Al27_to_Si28",  // 16,
        "He4_Al27_to_P31",  // 17,
        "He4_Si28_to_S32",  // 18,
        "p_Si29_to_P30",  // 19,
        "p_P30_to_S31",  // 20,
        "p_P31_to_S32",  // 21,
        "He4_Cr48_to_Fe52",  // 22,
        "p_nse_Mn51_to_Fe52",  // 23,
        "He4_Mn51_to_Co55",  // 24,
        "He4_Fe52_to_Ni56",  // 25,
        "p_nse_Co55_to_Ni56",  // 26,
        "He4_Co55_to_Cu59",  // 27,
        "He4_Ni56_to_Zn60",  // 28,
        "p_nse_Cu59_to_Zn60",  // 29,
        "C12_C12_to_p_Na23",  // 30,
        "C12_C12_to_He4_Ne20",  // 31,
        "He4_N13_to_p_O16",  // 32,
        "C12_O16_to_p_Al27",  // 33,
        "C12_O16_to_He4_Mg24",  // 34,
        "O16_O16_to_p_P31",  // 35,
        "O16_O16_to_He4_Si28",  // 36,
        "p_O17_to_He4_N14",  // 37,
        "He4_F18_to_p_Ne21",  // 38,
        "He4_Na22_to_p_Mg25",  // 39,
        "p_Na23_to_He4_Ne20",  // 40,
        "p_Al27_to_He4_Mg24",  // 41,
        "p_P31_to_He4_Si28",  // 42,
        "He4_Cr48_to_p_nse_Mn51",  // 43,
        "He4_Fe52_to_p_nse_Co55",  // 44,
        "p_nse_Cu59_to_He4_Ni56",  // 45,
        "He4_He4_He4_to_C12",  // 46,
        "C12_C12_to_Mg24_modified",  // 47,
        "O16_O16_to_S32_modified",  // 48,
        "C12_O16_to_Si28_modified",  // 49,
        "p_nse_Fe54_to_Co55",  // 50,
        "He4_Fe54_to_Ni58",  // 51,
        "p_nse_Fe56_to_Co57",  // 52,
        "n_Co55_to_Co56",  // 53,
        "n_Co56_to_Co57",  // 54,
        "p_nse_Co57_to_Ni58",  // 55,
        "p_nse_Ni58_to_Cu59",  // 56,
        "He4_Mn51_to_p_nse_Fe54",  // 57,
        "He4_Co55_to_p_nse_Ni58",  // 58,
        "n_Co56_to_p_nse_Fe56",  // 59,
        "p_nse_Co57_to_He4_Fe54",  // 60,
        "n_Ni56_to_p_nse_Co56",  // 61,
        "n_Cu59_to_He4_Co56",  // 62,
        "He4_S32_to_Ar36_removed",  // 63,
        "p_Cl35_to_Ar36_removed",  // 64,
        "p_Cl35_to_He4_S32_removed",  // 65,
        "He4_Ar36_to_Ca40_removed",  // 66,
        "p_K39_to_Ca40_removed",  // 67,
        "p_K39_to_He4_Ar36_removed",  // 68,
        "He4_Ca40_to_Ti44_removed",  // 69,
        "p_Sc43_to_Ti44_removed",  // 70,
        "p_Sc43_to_He4_Ca40_removed",  // 71,
        "He4_Ti44_to_Cr48_removed",  // 72,
        "He4_Ti44_to_p_V47_removed",  // 73,
        "p_V47_to_Cr48_removed",  // 74,
        "n_Fe52_to_Fe53_removed",  // 75,
        "n_Fe53_to_Fe54_removed",  // 76,
        "n_Fe54_to_Fe55_removed",  // 77,
        "n_Fe55_to_Fe56_removed",  // 78,
        "n_Ni56_to_Ni57_removed",  // 79,
        "n_Ni57_to_Ni58_removed",  // 80,
        "Co56_to_Fe56",  // 81,
        "Co56_to_Ni56",  // 82,
        "Fe56_to_Co56",  // 83,
        "n_to_p",  // 84,
        "Ni56_to_Co56",  // 85,
        "p_to_n",  // 86,
        "S32_He4_to_Ar36_approx",  // 87,
        "Ar36_to_S32_He4_approx",  // 88,
        "Ar36_He4_to_Ca40_approx",  // 89,
        "Ca40_to_Ar36_He4_approx",  // 90,
        "Ca40_He4_to_Ti44_approx",  // 91,
        "Ti44_to_Ca40_He4_approx",  // 92,
        "Ti44_He4_to_Cr48_approx",  // 93,
        "Cr48_to_Ti44_He4_approx",  // 94,
        "Fe52_n_n_to_Fe54_approx",  // 95,
        "Fe54_to_Fe52_n_n_approx",  // 96,
        "Fe54_n_n_to_Fe56_approx",  // 97,
        "Fe56_to_Fe54_n_n_approx",  // 98,
        "Ni56_n_n_to_Ni58_approx",  // 99,
        "Ni58_to_Ni56_n_n_approx",  // 100,
        "N13_to_p_C12_derived",  // 101,
        "O16_to_He4_C12_derived",  // 102,
        "F18_to_p_O17_derived",  // 103,
        "F18_to_He4_N14_derived",  // 104,
        "Ne20_to_He4_O16_derived",  // 105,
        "Ne21_to_He4_O17_derived",  // 106,
        "Na22_to_p_Ne21_derived",  // 107,
        "Na22_to_He4_F18_derived",  // 108,
        "Mg24_to_p_Na23_derived",  // 109,
        "Mg24_to_He4_Ne20_derived",  // 110,
        "Mg25_to_He4_Ne21_derived",  // 111,
        "Al27_to_He4_Na23_derived",  // 112,
        "Si28_to_p_Al27_derived",  // 113,
        "Si28_to_He4_Mg24_derived",  // 114,
        "Si29_to_He4_Mg25_derived",  // 115,
        "P30_to_p_Si29_derived",  // 116,
        "P31_to_He4_Al27_derived",  // 117,
        "S31_to_p_P30_derived",  // 118,
        "S32_to_p_P31_derived",  // 119,
        "S32_to_He4_Si28_derived",  // 120,
        "Fe52_to_p_nse_Mn51_derived",  // 121,
        "Fe52_to_He4_Cr48_derived",  // 122,
        "Co55_to_He4_Mn51_derived",  // 123,
        "Ni56_to_p_nse_Co55_derived",  // 124,
        "Ni56_to_He4_Fe52_derived",  // 125,
        "Cu59_to_He4_Co55_derived",  // 126,
        "Zn60_to_p_nse_Cu59_derived",  // 127,
        "Zn60_to_He4_Ni56_derived",  // 128,
        "C12_to_He4_He4_He4_derived",  // 129,
        "He4_N14_to_p_O17_derived",  // 130,
        "p_O16_to_He4_N13_derived",  // 131,
        "He4_Ne20_to_p_Na23_derived",  // 132,
        "p_Ne21_to_He4_F18_derived",  // 133,
        "He4_Mg24_to_p_Al27_derived",  // 134,
        "p_Mg25_to_He4_Na22_derived",  // 135,
        "He4_Si28_to_p_P31_derived",  // 136,
        "p_nse_Mn51_to_He4_Cr48_derived",  // 137,
        "p_nse_Co55_to_He4_Fe52_derived",  // 138,
        "He4_Ni56_to_p_nse_Cu59_derived",  // 139,
        "Co55_to_p_nse_Fe54_derived",  // 140,
        "Co56_to_n_Co55_derived",  // 141,
        "Co57_to_n_Co56_derived",  // 142,
        "Co57_to_p_nse_Fe56_derived",  // 143,
        "Ni58_to_p_nse_Co57_derived",  // 144,
        "Ni58_to_He4_Fe54_derived",  // 145,
        "Cu59_to_p_nse_Ni58_derived",  // 146,
        "p_nse_Fe54_to_He4_Mn51_derived",  // 147,
        "He4_Fe54_to_p_nse_Co57_derived",  // 148,
        "p_nse_Fe56_to_n_Co56_derived",  // 149,
        "p_nse_Co56_to_n_Ni56_derived",  // 150,
        "He4_Co56_to_n_Cu59_derived",  // 151,
        "p_nse_Ni58_to_He4_Co55_derived",  // 152,
        "He4_S32_to_p_Cl35_derived_removed",  // 153,
        "Ar36_to_He4_S32_derived_removed",  // 154,
        "Ar36_to_p_Cl35_derived_removed",  // 155,
        "He4_Ar36_to_p_K39_derived_removed",  // 156,
        "Ca40_to_He4_Ar36_derived_removed",  // 157,
        "Ca40_to_p_K39_derived_removed",  // 158,
        "He4_Ca40_to_p_Sc43_derived_removed",  // 159,
        "Ti44_to_He4_Ca40_derived_removed",  // 160,
        "Ti44_to_p_Sc43_derived_removed",  // 161,
        "Cr48_to_He4_Ti44_derived_removed",  // 162,
        "Cr48_to_p_V47_derived_removed",  // 163,
        "p_V47_to_He4_Ti44_derived_removed",  // 164,
        "Fe54_to_n_Fe53_derived_removed",  // 165,
        "Fe53_to_n_Fe52_derived_removed",  // 166,
        "Fe56_to_n_Fe55_derived_removed",  // 167,
        "Fe55_to_n_Fe54_derived_removed",  // 168,
        "Ni58_to_n_Ni57_derived_removed",  // 169,
        "Ni57_to_n_Ni56_derived_removed"  // 170,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int H1_index = 1;
    constexpr int N_index = 0;
    constexpr int He4_index = 3;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
