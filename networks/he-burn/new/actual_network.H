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
        else if constexpr (spec == F19) {
            return 147.80136567599766_rt;
        }
        else if constexpr (spec == Ne20) {
            return 160.64482384000075_rt;
        }
        else if constexpr (spec == Ne21) {
            return 167.40598973999658_rt;
        }
        else if constexpr (spec == Ne22) {
            return 177.77024384000106_rt;
        }
        else if constexpr (spec == Na22) {
            return 174.14457080400098_rt;
        }
        else if constexpr (spec == Na23) {
            return 186.56435240400242_rt;
        }
        else if constexpr (spec == Mg23) {
            return 181.7258218679999_rt;
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

    enum NetworkRates :
    std::uint8_t
    {
        k_p_C12_to_N13 = 1,
        k_He4_C12_to_O16 = 2,
        k_He4_N14_to_F18 = 3,
        k_He4_O16_to_Ne20 = 4,
        k_p_O17_to_F18 = 5,
        k_He4_O17_to_Ne21 = 6,
        k_He4_F18_to_Na22 = 7,
        k_p_F19_to_Ne20 = 8,
        k_He4_F19_to_Na23 = 9,
        k_He4_Ne20_to_Mg24 = 10,
        k_p_Ne21_to_Na22 = 11,
        k_p_Ne22_to_Na23 = 12,
        k_p_Na22_to_Mg23 = 13,
        k_p_Na23_to_Mg24 = 14,
        k_He4_Mg24_to_Si28 = 15,
        k_p_Al27_to_Si28 = 16,
        k_He4_Si28_to_S32 = 17,
        k_p_P31_to_S32 = 18,
        k_He4_Cr48_to_Fe52 = 19,
        k_p_nse_Mn51_to_Fe52 = 20,
        k_He4_Mn51_to_Co55 = 21,
        k_He4_Fe52_to_Ni56 = 22,
        k_p_nse_Co55_to_Ni56 = 23,
        k_He4_Co55_to_Cu59 = 24,
        k_He4_Ni56_to_Zn60 = 25,
        k_p_nse_Cu59_to_Zn60 = 26,
        k_C12_C12_to_p_Na23 = 27,
        k_C12_C12_to_He4_Ne20 = 28,
        k_He4_N13_to_p_O16 = 29,
        k_C12_O16_to_p_Al27 = 30,
        k_C12_O16_to_He4_Mg24 = 31,
        k_O16_O16_to_p_P31 = 32,
        k_O16_O16_to_He4_Si28 = 33,
        k_p_O17_to_He4_N14 = 34,
        k_He4_F18_to_p_Ne21 = 35,
        k_p_F19_to_He4_O16 = 36,
        k_He4_F19_to_p_Ne22 = 37,
        k_p_Na23_to_He4_Ne20 = 38,
        k_p_Al27_to_He4_Mg24 = 39,
        k_p_P31_to_He4_Si28 = 40,
        k_He4_Cr48_to_p_nse_Mn51 = 41,
        k_He4_Fe52_to_p_nse_Co55 = 42,
        k_p_nse_Cu59_to_He4_Ni56 = 43,
        k_He4_He4_He4_to_C12 = 44,
        k_C12_C12_to_Mg24_modified = 45,
        k_O16_O16_to_S32_modified = 46,
        k_C12_O16_to_Si28_modified = 47,
        k_p_nse_Fe54_to_Co55 = 48,
        k_He4_Fe54_to_Ni58 = 49,
        k_p_nse_Fe56_to_Co57 = 50,
        k_n_Co55_to_Co56 = 51,
        k_n_Co56_to_Co57 = 52,
        k_p_nse_Co57_to_Ni58 = 53,
        k_p_nse_Ni58_to_Cu59 = 54,
        k_He4_Mn51_to_p_nse_Fe54 = 55,
        k_He4_Co55_to_p_nse_Ni58 = 56,
        k_n_Co56_to_p_nse_Fe56 = 57,
        k_p_nse_Co57_to_He4_Fe54 = 58,
        k_n_Ni56_to_p_nse_Co56 = 59,
        k_n_Cu59_to_He4_Co56 = 60,
        k_He4_S32_to_Ar36_removed = 61,
        k_p_Cl35_to_Ar36_removed = 62,
        k_p_Cl35_to_He4_S32_removed = 63,
        k_He4_Ar36_to_Ca40_removed = 64,
        k_p_K39_to_Ca40_removed = 65,
        k_p_K39_to_He4_Ar36_removed = 66,
        k_He4_Ca40_to_Ti44_removed = 67,
        k_p_Sc43_to_Ti44_removed = 68,
        k_p_Sc43_to_He4_Ca40_removed = 69,
        k_He4_Ti44_to_Cr48_removed = 70,
        k_He4_Ti44_to_p_V47_removed = 71,
        k_p_V47_to_Cr48_removed = 72,
        k_n_Fe52_to_Fe53_removed = 73,
        k_n_Fe53_to_Fe54_removed = 74,
        k_n_Fe54_to_Fe55_removed = 75,
        k_n_Fe55_to_Fe56_removed = 76,
        k_n_Ni56_to_Ni57_removed = 77,
        k_n_Ni57_to_Ni58_removed = 78,
        k_Mg23_to_Na23 = 79,
        k_Na22_to_Ne22 = 80,
        k_Na23_to_Mg23 = 81,
        k_Ne22_to_Na22 = 82,
        k_Co56_to_Fe56 = 83,
        k_Co56_to_Ni56 = 84,
        k_Fe56_to_Co56 = 85,
        k_n_to_p = 86,
        k_Ni56_to_Co56 = 87,
        k_p_to_n = 88,
        k_S32_He4_to_Ar36_approx = 89,
        k_Ar36_to_S32_He4_approx = 90,
        k_Ar36_He4_to_Ca40_approx = 91,
        k_Ca40_to_Ar36_He4_approx = 92,
        k_Ca40_He4_to_Ti44_approx = 93,
        k_Ti44_to_Ca40_He4_approx = 94,
        k_Ti44_He4_to_Cr48_approx = 95,
        k_Cr48_to_Ti44_He4_approx = 96,
        k_Fe52_n_n_to_Fe54_approx = 97,
        k_Fe54_to_Fe52_n_n_approx = 98,
        k_Fe54_n_n_to_Fe56_approx = 99,
        k_Fe56_to_Fe54_n_n_approx = 100,
        k_Ni56_n_n_to_Ni58_approx = 101,
        k_Ni58_to_Ni56_n_n_approx = 102,
        k_N13_to_p_C12_derived = 103,
        k_O16_to_He4_C12_derived = 104,
        k_F18_to_p_O17_derived = 105,
        k_F18_to_He4_N14_derived = 106,
        k_Ne20_to_p_F19_derived = 107,
        k_Ne20_to_He4_O16_derived = 108,
        k_Ne21_to_He4_O17_derived = 109,
        k_Na22_to_p_Ne21_derived = 110,
        k_Na22_to_He4_F18_derived = 111,
        k_Na23_to_p_Ne22_derived = 112,
        k_Na23_to_He4_F19_derived = 113,
        k_Mg23_to_p_Na22_derived = 114,
        k_Mg24_to_p_Na23_derived = 115,
        k_Mg24_to_He4_Ne20_derived = 116,
        k_Si28_to_p_Al27_derived = 117,
        k_Si28_to_He4_Mg24_derived = 118,
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
        k_He4_O16_to_p_F19_derived = 132,
        k_He4_Ne20_to_p_Na23_derived = 133,
        k_p_Ne21_to_He4_F18_derived = 134,
        k_p_Ne22_to_He4_F19_derived = 135,
        k_He4_Mg24_to_p_Al27_derived = 136,
        k_He4_Si28_to_p_P31_derived = 137,
        k_p_nse_Mn51_to_He4_Cr48_derived = 138,
        k_p_nse_Co55_to_He4_Fe52_derived = 139,
        k_He4_Ni56_to_p_nse_Cu59_derived = 140,
        k_Co55_to_p_nse_Fe54_derived = 141,
        k_Co56_to_n_Co55_derived = 142,
        k_Co57_to_n_Co56_derived = 143,
        k_Co57_to_p_nse_Fe56_derived = 144,
        k_Ni58_to_p_nse_Co57_derived = 145,
        k_Ni58_to_He4_Fe54_derived = 146,
        k_Cu59_to_p_nse_Ni58_derived = 147,
        k_p_nse_Fe54_to_He4_Mn51_derived = 148,
        k_He4_Fe54_to_p_nse_Co57_derived = 149,
        k_p_nse_Fe56_to_n_Co56_derived = 150,
        k_p_nse_Co56_to_n_Ni56_derived = 151,
        k_He4_Co56_to_n_Cu59_derived = 152,
        k_p_nse_Ni58_to_He4_Co55_derived = 153,
        k_He4_S32_to_p_Cl35_derived_removed = 154,
        k_Ar36_to_He4_S32_derived_removed = 155,
        k_Ar36_to_p_Cl35_derived_removed = 156,
        k_He4_Ar36_to_p_K39_derived_removed = 157,
        k_Ca40_to_He4_Ar36_derived_removed = 158,
        k_Ca40_to_p_K39_derived_removed = 159,
        k_He4_Ca40_to_p_Sc43_derived_removed = 160,
        k_Ti44_to_He4_Ca40_derived_removed = 161,
        k_Ti44_to_p_Sc43_derived_removed = 162,
        k_Cr48_to_He4_Ti44_derived_removed = 163,
        k_Cr48_to_p_V47_derived_removed = 164,
        k_p_V47_to_He4_Ti44_derived_removed = 165,
        k_Fe54_to_n_Fe53_derived_removed = 166,
        k_Fe53_to_n_Fe52_derived_removed = 167,
        k_Fe56_to_n_Fe55_derived_removed = 168,
        k_Fe55_to_n_Fe54_derived_removed = 169,
        k_Ni58_to_n_Ni57_derived_removed = 170,
        k_Ni57_to_n_Ni56_derived_removed = 171,
        NumRates = k_Ni57_to_n_Ni56_derived_removed
    };

    // number of reaclib rates

    const int NrateReaclib = 147;

    // number of tabular rates

    const int NrateTabular = 10;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "p_C12_to_N13",  // 1,
        "He4_C12_to_O16",  // 2,
        "He4_N14_to_F18",  // 3,
        "He4_O16_to_Ne20",  // 4,
        "p_O17_to_F18",  // 5,
        "He4_O17_to_Ne21",  // 6,
        "He4_F18_to_Na22",  // 7,
        "p_F19_to_Ne20",  // 8,
        "He4_F19_to_Na23",  // 9,
        "He4_Ne20_to_Mg24",  // 10,
        "p_Ne21_to_Na22",  // 11,
        "p_Ne22_to_Na23",  // 12,
        "p_Na22_to_Mg23",  // 13,
        "p_Na23_to_Mg24",  // 14,
        "He4_Mg24_to_Si28",  // 15,
        "p_Al27_to_Si28",  // 16,
        "He4_Si28_to_S32",  // 17,
        "p_P31_to_S32",  // 18,
        "He4_Cr48_to_Fe52",  // 19,
        "p_nse_Mn51_to_Fe52",  // 20,
        "He4_Mn51_to_Co55",  // 21,
        "He4_Fe52_to_Ni56",  // 22,
        "p_nse_Co55_to_Ni56",  // 23,
        "He4_Co55_to_Cu59",  // 24,
        "He4_Ni56_to_Zn60",  // 25,
        "p_nse_Cu59_to_Zn60",  // 26,
        "C12_C12_to_p_Na23",  // 27,
        "C12_C12_to_He4_Ne20",  // 28,
        "He4_N13_to_p_O16",  // 29,
        "C12_O16_to_p_Al27",  // 30,
        "C12_O16_to_He4_Mg24",  // 31,
        "O16_O16_to_p_P31",  // 32,
        "O16_O16_to_He4_Si28",  // 33,
        "p_O17_to_He4_N14",  // 34,
        "He4_F18_to_p_Ne21",  // 35,
        "p_F19_to_He4_O16",  // 36,
        "He4_F19_to_p_Ne22",  // 37,
        "p_Na23_to_He4_Ne20",  // 38,
        "p_Al27_to_He4_Mg24",  // 39,
        "p_P31_to_He4_Si28",  // 40,
        "He4_Cr48_to_p_nse_Mn51",  // 41,
        "He4_Fe52_to_p_nse_Co55",  // 42,
        "p_nse_Cu59_to_He4_Ni56",  // 43,
        "He4_He4_He4_to_C12",  // 44,
        "C12_C12_to_Mg24_modified",  // 45,
        "O16_O16_to_S32_modified",  // 46,
        "C12_O16_to_Si28_modified",  // 47,
        "p_nse_Fe54_to_Co55",  // 48,
        "He4_Fe54_to_Ni58",  // 49,
        "p_nse_Fe56_to_Co57",  // 50,
        "n_Co55_to_Co56",  // 51,
        "n_Co56_to_Co57",  // 52,
        "p_nse_Co57_to_Ni58",  // 53,
        "p_nse_Ni58_to_Cu59",  // 54,
        "He4_Mn51_to_p_nse_Fe54",  // 55,
        "He4_Co55_to_p_nse_Ni58",  // 56,
        "n_Co56_to_p_nse_Fe56",  // 57,
        "p_nse_Co57_to_He4_Fe54",  // 58,
        "n_Ni56_to_p_nse_Co56",  // 59,
        "n_Cu59_to_He4_Co56",  // 60,
        "He4_S32_to_Ar36_removed",  // 61,
        "p_Cl35_to_Ar36_removed",  // 62,
        "p_Cl35_to_He4_S32_removed",  // 63,
        "He4_Ar36_to_Ca40_removed",  // 64,
        "p_K39_to_Ca40_removed",  // 65,
        "p_K39_to_He4_Ar36_removed",  // 66,
        "He4_Ca40_to_Ti44_removed",  // 67,
        "p_Sc43_to_Ti44_removed",  // 68,
        "p_Sc43_to_He4_Ca40_removed",  // 69,
        "He4_Ti44_to_Cr48_removed",  // 70,
        "He4_Ti44_to_p_V47_removed",  // 71,
        "p_V47_to_Cr48_removed",  // 72,
        "n_Fe52_to_Fe53_removed",  // 73,
        "n_Fe53_to_Fe54_removed",  // 74,
        "n_Fe54_to_Fe55_removed",  // 75,
        "n_Fe55_to_Fe56_removed",  // 76,
        "n_Ni56_to_Ni57_removed",  // 77,
        "n_Ni57_to_Ni58_removed",  // 78,
        "Mg23_to_Na23",  // 79,
        "Na22_to_Ne22",  // 80,
        "Na23_to_Mg23",  // 81,
        "Ne22_to_Na22",  // 82,
        "Co56_to_Fe56",  // 83,
        "Co56_to_Ni56",  // 84,
        "Fe56_to_Co56",  // 85,
        "n_to_p",  // 86,
        "Ni56_to_Co56",  // 87,
        "p_to_n",  // 88,
        "S32_He4_to_Ar36_approx",  // 89,
        "Ar36_to_S32_He4_approx",  // 90,
        "Ar36_He4_to_Ca40_approx",  // 91,
        "Ca40_to_Ar36_He4_approx",  // 92,
        "Ca40_He4_to_Ti44_approx",  // 93,
        "Ti44_to_Ca40_He4_approx",  // 94,
        "Ti44_He4_to_Cr48_approx",  // 95,
        "Cr48_to_Ti44_He4_approx",  // 96,
        "Fe52_n_n_to_Fe54_approx",  // 97,
        "Fe54_to_Fe52_n_n_approx",  // 98,
        "Fe54_n_n_to_Fe56_approx",  // 99,
        "Fe56_to_Fe54_n_n_approx",  // 100,
        "Ni56_n_n_to_Ni58_approx",  // 101,
        "Ni58_to_Ni56_n_n_approx",  // 102,
        "N13_to_p_C12_derived",  // 103,
        "O16_to_He4_C12_derived",  // 104,
        "F18_to_p_O17_derived",  // 105,
        "F18_to_He4_N14_derived",  // 106,
        "Ne20_to_p_F19_derived",  // 107,
        "Ne20_to_He4_O16_derived",  // 108,
        "Ne21_to_He4_O17_derived",  // 109,
        "Na22_to_p_Ne21_derived",  // 110,
        "Na22_to_He4_F18_derived",  // 111,
        "Na23_to_p_Ne22_derived",  // 112,
        "Na23_to_He4_F19_derived",  // 113,
        "Mg23_to_p_Na22_derived",  // 114,
        "Mg24_to_p_Na23_derived",  // 115,
        "Mg24_to_He4_Ne20_derived",  // 116,
        "Si28_to_p_Al27_derived",  // 117,
        "Si28_to_He4_Mg24_derived",  // 118,
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
        "He4_O16_to_p_F19_derived",  // 132,
        "He4_Ne20_to_p_Na23_derived",  // 133,
        "p_Ne21_to_He4_F18_derived",  // 134,
        "p_Ne22_to_He4_F19_derived",  // 135,
        "He4_Mg24_to_p_Al27_derived",  // 136,
        "He4_Si28_to_p_P31_derived",  // 137,
        "p_nse_Mn51_to_He4_Cr48_derived",  // 138,
        "p_nse_Co55_to_He4_Fe52_derived",  // 139,
        "He4_Ni56_to_p_nse_Cu59_derived",  // 140,
        "Co55_to_p_nse_Fe54_derived",  // 141,
        "Co56_to_n_Co55_derived",  // 142,
        "Co57_to_n_Co56_derived",  // 143,
        "Co57_to_p_nse_Fe56_derived",  // 144,
        "Ni58_to_p_nse_Co57_derived",  // 145,
        "Ni58_to_He4_Fe54_derived",  // 146,
        "Cu59_to_p_nse_Ni58_derived",  // 147,
        "p_nse_Fe54_to_He4_Mn51_derived",  // 148,
        "He4_Fe54_to_p_nse_Co57_derived",  // 149,
        "p_nse_Fe56_to_n_Co56_derived",  // 150,
        "p_nse_Co56_to_n_Ni56_derived",  // 151,
        "He4_Co56_to_n_Cu59_derived",  // 152,
        "p_nse_Ni58_to_He4_Co55_derived",  // 153,
        "He4_S32_to_p_Cl35_derived_removed",  // 154,
        "Ar36_to_He4_S32_derived_removed",  // 155,
        "Ar36_to_p_Cl35_derived_removed",  // 156,
        "He4_Ar36_to_p_K39_derived_removed",  // 157,
        "Ca40_to_He4_Ar36_derived_removed",  // 158,
        "Ca40_to_p_K39_derived_removed",  // 159,
        "He4_Ca40_to_p_Sc43_derived_removed",  // 160,
        "Ti44_to_He4_Ca40_derived_removed",  // 161,
        "Ti44_to_p_Sc43_derived_removed",  // 162,
        "Cr48_to_He4_Ti44_derived_removed",  // 163,
        "Cr48_to_p_V47_derived_removed",  // 164,
        "p_V47_to_He4_Ti44_derived_removed",  // 165,
        "Fe54_to_n_Fe53_derived_removed",  // 166,
        "Fe53_to_n_Fe52_derived_removed",  // 167,
        "Fe56_to_n_Fe55_derived_removed",  // 168,
        "Fe55_to_n_Fe54_derived_removed",  // 169,
        "Ni58_to_n_Ni57_derived_removed",  // 170,
        "Ni57_to_n_Ni56_derived_removed"  // 171,
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
