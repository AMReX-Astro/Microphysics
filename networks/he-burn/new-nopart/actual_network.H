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
        k_N13_to_p_C12 = 1,
        k_O16_to_He4_C12 = 2,
        k_F18_to_p_O17 = 3,
        k_F18_to_He4_N14 = 4,
        k_Ne20_to_p_F19 = 5,
        k_Ne20_to_He4_O16 = 6,
        k_Ne21_to_He4_O17 = 7,
        k_Na22_to_p_Ne21 = 8,
        k_Na22_to_He4_F18 = 9,
        k_Na23_to_p_Ne22 = 10,
        k_Na23_to_He4_F19 = 11,
        k_Mg23_to_p_Na22 = 12,
        k_Mg24_to_p_Na23 = 13,
        k_Mg24_to_He4_Ne20 = 14,
        k_Si28_to_p_Al27 = 15,
        k_Si28_to_He4_Mg24 = 16,
        k_S32_to_p_P31 = 17,
        k_S32_to_He4_Si28 = 18,
        k_Fe52_to_p_nse_Mn51 = 19,
        k_Fe52_to_He4_Cr48 = 20,
        k_Co55_to_He4_Mn51 = 21,
        k_Ni56_to_p_nse_Co55 = 22,
        k_Ni56_to_He4_Fe52 = 23,
        k_Cu59_to_He4_Co55 = 24,
        k_Zn60_to_p_nse_Cu59 = 25,
        k_Zn60_to_He4_Ni56 = 26,
        k_C12_to_He4_He4_He4 = 27,
        k_p_C12_to_N13 = 28,
        k_He4_C12_to_O16 = 29,
        k_He4_N14_to_F18 = 30,
        k_He4_O16_to_Ne20 = 31,
        k_p_O17_to_F18 = 32,
        k_He4_O17_to_Ne21 = 33,
        k_He4_F18_to_Na22 = 34,
        k_p_F19_to_Ne20 = 35,
        k_He4_F19_to_Na23 = 36,
        k_He4_Ne20_to_Mg24 = 37,
        k_p_Ne21_to_Na22 = 38,
        k_p_Ne22_to_Na23 = 39,
        k_p_Na22_to_Mg23 = 40,
        k_p_Na23_to_Mg24 = 41,
        k_He4_Mg24_to_Si28 = 42,
        k_p_Al27_to_Si28 = 43,
        k_He4_Si28_to_S32 = 44,
        k_p_P31_to_S32 = 45,
        k_He4_Cr48_to_Fe52 = 46,
        k_p_nse_Mn51_to_Fe52 = 47,
        k_He4_Mn51_to_Co55 = 48,
        k_He4_Fe52_to_Ni56 = 49,
        k_p_nse_Co55_to_Ni56 = 50,
        k_He4_Co55_to_Cu59 = 51,
        k_He4_Ni56_to_Zn60 = 52,
        k_p_nse_Cu59_to_Zn60 = 53,
        k_C12_C12_to_p_Na23 = 54,
        k_C12_C12_to_He4_Ne20 = 55,
        k_He4_N13_to_p_O16 = 56,
        k_He4_N14_to_p_O17 = 57,
        k_p_O16_to_He4_N13 = 58,
        k_He4_O16_to_p_F19 = 59,
        k_C12_O16_to_p_Al27 = 60,
        k_C12_O16_to_He4_Mg24 = 61,
        k_O16_O16_to_p_P31 = 62,
        k_O16_O16_to_He4_Si28 = 63,
        k_p_O17_to_He4_N14 = 64,
        k_He4_F18_to_p_Ne21 = 65,
        k_p_F19_to_He4_O16 = 66,
        k_He4_F19_to_p_Ne22 = 67,
        k_He4_Ne20_to_p_Na23 = 68,
        k_p_Ne21_to_He4_F18 = 69,
        k_p_Ne22_to_He4_F19 = 70,
        k_p_Na23_to_He4_Ne20 = 71,
        k_He4_Mg24_to_p_Al27 = 72,
        k_p_Al27_to_He4_Mg24 = 73,
        k_He4_Si28_to_p_P31 = 74,
        k_p_P31_to_He4_Si28 = 75,
        k_He4_Cr48_to_p_nse_Mn51 = 76,
        k_p_nse_Mn51_to_He4_Cr48 = 77,
        k_He4_Fe52_to_p_nse_Co55 = 78,
        k_p_nse_Co55_to_He4_Fe52 = 79,
        k_He4_Ni56_to_p_nse_Cu59 = 80,
        k_p_nse_Cu59_to_He4_Ni56 = 81,
        k_He4_He4_He4_to_C12 = 82,
        k_C12_C12_to_Mg24_modified = 83,
        k_O16_O16_to_S32_modified = 84,
        k_C12_O16_to_Si28_modified = 85,
        k_Co55_to_p_nse_Fe54 = 86,
        k_Co56_to_n_Co55 = 87,
        k_Co57_to_n_Co56 = 88,
        k_Co57_to_p_nse_Fe56 = 89,
        k_Ni58_to_p_nse_Co57 = 90,
        k_Ni58_to_He4_Fe54 = 91,
        k_Cu59_to_p_nse_Ni58 = 92,
        k_p_nse_Fe54_to_Co55 = 93,
        k_He4_Fe54_to_Ni58 = 94,
        k_p_nse_Fe56_to_Co57 = 95,
        k_n_Co55_to_Co56 = 96,
        k_n_Co56_to_Co57 = 97,
        k_p_nse_Co57_to_Ni58 = 98,
        k_p_nse_Ni58_to_Cu59 = 99,
        k_He4_Mn51_to_p_nse_Fe54 = 100,
        k_p_nse_Fe54_to_He4_Mn51 = 101,
        k_He4_Fe54_to_p_nse_Co57 = 102,
        k_p_nse_Fe56_to_n_Co56 = 103,
        k_He4_Co55_to_p_nse_Ni58 = 104,
        k_n_Co56_to_p_nse_Fe56 = 105,
        k_p_nse_Co56_to_n_Ni56 = 106,
        k_He4_Co56_to_n_Cu59 = 107,
        k_p_nse_Co57_to_He4_Fe54 = 108,
        k_n_Ni56_to_p_nse_Co56 = 109,
        k_p_nse_Ni58_to_He4_Co55 = 110,
        k_n_Cu59_to_He4_Co56 = 111,
        k_He4_S32_to_Ar36_removed = 112,
        k_He4_S32_to_p_Cl35_removed = 113,
        k_p_Cl35_to_Ar36_removed = 114,
        k_Ar36_to_He4_S32_removed = 115,
        k_Ar36_to_p_Cl35_removed = 116,
        k_p_Cl35_to_He4_S32_removed = 117,
        k_He4_Ar36_to_Ca40_removed = 118,
        k_He4_Ar36_to_p_K39_removed = 119,
        k_p_K39_to_Ca40_removed = 120,
        k_Ca40_to_He4_Ar36_removed = 121,
        k_Ca40_to_p_K39_removed = 122,
        k_p_K39_to_He4_Ar36_removed = 123,
        k_He4_Ca40_to_Ti44_removed = 124,
        k_He4_Ca40_to_p_Sc43_removed = 125,
        k_p_Sc43_to_Ti44_removed = 126,
        k_Ti44_to_He4_Ca40_removed = 127,
        k_Ti44_to_p_Sc43_removed = 128,
        k_p_Sc43_to_He4_Ca40_removed = 129,
        k_He4_Ti44_to_Cr48_removed = 130,
        k_He4_Ti44_to_p_V47_removed = 131,
        k_p_V47_to_Cr48_removed = 132,
        k_Cr48_to_He4_Ti44_removed = 133,
        k_Cr48_to_p_V47_removed = 134,
        k_p_V47_to_He4_Ti44_removed = 135,
        k_n_Fe52_to_Fe53_removed = 136,
        k_n_Fe53_to_Fe54_removed = 137,
        k_Fe54_to_n_Fe53_removed = 138,
        k_Fe53_to_n_Fe52_removed = 139,
        k_n_Fe54_to_Fe55_removed = 140,
        k_n_Fe55_to_Fe56_removed = 141,
        k_Fe56_to_n_Fe55_removed = 142,
        k_Fe55_to_n_Fe54_removed = 143,
        k_n_Ni56_to_Ni57_removed = 144,
        k_n_Ni57_to_Ni58_removed = 145,
        k_Ni58_to_n_Ni57_removed = 146,
        k_Ni57_to_n_Ni56_removed = 147,
        k_Mg23_to_Na23 = 148,
        k_Na22_to_Ne22 = 149,
        k_Na23_to_Mg23 = 150,
        k_Ne22_to_Na22 = 151,
        k_Co56_to_Fe56 = 152,
        k_Co56_to_Ni56 = 153,
        k_Fe56_to_Co56 = 154,
        k_n_to_p = 155,
        k_Ni56_to_Co56 = 156,
        k_p_to_n = 157,
        k_S32_He4_to_Ar36_approx = 158,
        k_Ar36_to_S32_He4_approx = 159,
        k_Ar36_He4_to_Ca40_approx = 160,
        k_Ca40_to_Ar36_He4_approx = 161,
        k_Ca40_He4_to_Ti44_approx = 162,
        k_Ti44_to_Ca40_He4_approx = 163,
        k_Ti44_He4_to_Cr48_approx = 164,
        k_Cr48_to_Ti44_He4_approx = 165,
        k_Fe52_n_n_to_Fe54_approx = 166,
        k_Fe54_to_Fe52_n_n_approx = 167,
        k_Fe54_n_n_to_Fe56_approx = 168,
        k_Fe56_to_Fe54_n_n_approx = 169,
        k_Ni56_n_n_to_Ni58_approx = 170,
        k_Ni58_to_Ni56_n_n_approx = 171,
        NumRates = k_Ni58_to_Ni56_n_n_approx
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
        "N13_to_p_C12",  // 1,
        "O16_to_He4_C12",  // 2,
        "F18_to_p_O17",  // 3,
        "F18_to_He4_N14",  // 4,
        "Ne20_to_p_F19",  // 5,
        "Ne20_to_He4_O16",  // 6,
        "Ne21_to_He4_O17",  // 7,
        "Na22_to_p_Ne21",  // 8,
        "Na22_to_He4_F18",  // 9,
        "Na23_to_p_Ne22",  // 10,
        "Na23_to_He4_F19",  // 11,
        "Mg23_to_p_Na22",  // 12,
        "Mg24_to_p_Na23",  // 13,
        "Mg24_to_He4_Ne20",  // 14,
        "Si28_to_p_Al27",  // 15,
        "Si28_to_He4_Mg24",  // 16,
        "S32_to_p_P31",  // 17,
        "S32_to_He4_Si28",  // 18,
        "Fe52_to_p_nse_Mn51",  // 19,
        "Fe52_to_He4_Cr48",  // 20,
        "Co55_to_He4_Mn51",  // 21,
        "Ni56_to_p_nse_Co55",  // 22,
        "Ni56_to_He4_Fe52",  // 23,
        "Cu59_to_He4_Co55",  // 24,
        "Zn60_to_p_nse_Cu59",  // 25,
        "Zn60_to_He4_Ni56",  // 26,
        "C12_to_He4_He4_He4",  // 27,
        "p_C12_to_N13",  // 28,
        "He4_C12_to_O16",  // 29,
        "He4_N14_to_F18",  // 30,
        "He4_O16_to_Ne20",  // 31,
        "p_O17_to_F18",  // 32,
        "He4_O17_to_Ne21",  // 33,
        "He4_F18_to_Na22",  // 34,
        "p_F19_to_Ne20",  // 35,
        "He4_F19_to_Na23",  // 36,
        "He4_Ne20_to_Mg24",  // 37,
        "p_Ne21_to_Na22",  // 38,
        "p_Ne22_to_Na23",  // 39,
        "p_Na22_to_Mg23",  // 40,
        "p_Na23_to_Mg24",  // 41,
        "He4_Mg24_to_Si28",  // 42,
        "p_Al27_to_Si28",  // 43,
        "He4_Si28_to_S32",  // 44,
        "p_P31_to_S32",  // 45,
        "He4_Cr48_to_Fe52",  // 46,
        "p_nse_Mn51_to_Fe52",  // 47,
        "He4_Mn51_to_Co55",  // 48,
        "He4_Fe52_to_Ni56",  // 49,
        "p_nse_Co55_to_Ni56",  // 50,
        "He4_Co55_to_Cu59",  // 51,
        "He4_Ni56_to_Zn60",  // 52,
        "p_nse_Cu59_to_Zn60",  // 53,
        "C12_C12_to_p_Na23",  // 54,
        "C12_C12_to_He4_Ne20",  // 55,
        "He4_N13_to_p_O16",  // 56,
        "He4_N14_to_p_O17",  // 57,
        "p_O16_to_He4_N13",  // 58,
        "He4_O16_to_p_F19",  // 59,
        "C12_O16_to_p_Al27",  // 60,
        "C12_O16_to_He4_Mg24",  // 61,
        "O16_O16_to_p_P31",  // 62,
        "O16_O16_to_He4_Si28",  // 63,
        "p_O17_to_He4_N14",  // 64,
        "He4_F18_to_p_Ne21",  // 65,
        "p_F19_to_He4_O16",  // 66,
        "He4_F19_to_p_Ne22",  // 67,
        "He4_Ne20_to_p_Na23",  // 68,
        "p_Ne21_to_He4_F18",  // 69,
        "p_Ne22_to_He4_F19",  // 70,
        "p_Na23_to_He4_Ne20",  // 71,
        "He4_Mg24_to_p_Al27",  // 72,
        "p_Al27_to_He4_Mg24",  // 73,
        "He4_Si28_to_p_P31",  // 74,
        "p_P31_to_He4_Si28",  // 75,
        "He4_Cr48_to_p_nse_Mn51",  // 76,
        "p_nse_Mn51_to_He4_Cr48",  // 77,
        "He4_Fe52_to_p_nse_Co55",  // 78,
        "p_nse_Co55_to_He4_Fe52",  // 79,
        "He4_Ni56_to_p_nse_Cu59",  // 80,
        "p_nse_Cu59_to_He4_Ni56",  // 81,
        "He4_He4_He4_to_C12",  // 82,
        "C12_C12_to_Mg24_modified",  // 83,
        "O16_O16_to_S32_modified",  // 84,
        "C12_O16_to_Si28_modified",  // 85,
        "Co55_to_p_nse_Fe54",  // 86,
        "Co56_to_n_Co55",  // 87,
        "Co57_to_n_Co56",  // 88,
        "Co57_to_p_nse_Fe56",  // 89,
        "Ni58_to_p_nse_Co57",  // 90,
        "Ni58_to_He4_Fe54",  // 91,
        "Cu59_to_p_nse_Ni58",  // 92,
        "p_nse_Fe54_to_Co55",  // 93,
        "He4_Fe54_to_Ni58",  // 94,
        "p_nse_Fe56_to_Co57",  // 95,
        "n_Co55_to_Co56",  // 96,
        "n_Co56_to_Co57",  // 97,
        "p_nse_Co57_to_Ni58",  // 98,
        "p_nse_Ni58_to_Cu59",  // 99,
        "He4_Mn51_to_p_nse_Fe54",  // 100,
        "p_nse_Fe54_to_He4_Mn51",  // 101,
        "He4_Fe54_to_p_nse_Co57",  // 102,
        "p_nse_Fe56_to_n_Co56",  // 103,
        "He4_Co55_to_p_nse_Ni58",  // 104,
        "n_Co56_to_p_nse_Fe56",  // 105,
        "p_nse_Co56_to_n_Ni56",  // 106,
        "He4_Co56_to_n_Cu59",  // 107,
        "p_nse_Co57_to_He4_Fe54",  // 108,
        "n_Ni56_to_p_nse_Co56",  // 109,
        "p_nse_Ni58_to_He4_Co55",  // 110,
        "n_Cu59_to_He4_Co56",  // 111,
        "He4_S32_to_Ar36_removed",  // 112,
        "He4_S32_to_p_Cl35_removed",  // 113,
        "p_Cl35_to_Ar36_removed",  // 114,
        "Ar36_to_He4_S32_removed",  // 115,
        "Ar36_to_p_Cl35_removed",  // 116,
        "p_Cl35_to_He4_S32_removed",  // 117,
        "He4_Ar36_to_Ca40_removed",  // 118,
        "He4_Ar36_to_p_K39_removed",  // 119,
        "p_K39_to_Ca40_removed",  // 120,
        "Ca40_to_He4_Ar36_removed",  // 121,
        "Ca40_to_p_K39_removed",  // 122,
        "p_K39_to_He4_Ar36_removed",  // 123,
        "He4_Ca40_to_Ti44_removed",  // 124,
        "He4_Ca40_to_p_Sc43_removed",  // 125,
        "p_Sc43_to_Ti44_removed",  // 126,
        "Ti44_to_He4_Ca40_removed",  // 127,
        "Ti44_to_p_Sc43_removed",  // 128,
        "p_Sc43_to_He4_Ca40_removed",  // 129,
        "He4_Ti44_to_Cr48_removed",  // 130,
        "He4_Ti44_to_p_V47_removed",  // 131,
        "p_V47_to_Cr48_removed",  // 132,
        "Cr48_to_He4_Ti44_removed",  // 133,
        "Cr48_to_p_V47_removed",  // 134,
        "p_V47_to_He4_Ti44_removed",  // 135,
        "n_Fe52_to_Fe53_removed",  // 136,
        "n_Fe53_to_Fe54_removed",  // 137,
        "Fe54_to_n_Fe53_removed",  // 138,
        "Fe53_to_n_Fe52_removed",  // 139,
        "n_Fe54_to_Fe55_removed",  // 140,
        "n_Fe55_to_Fe56_removed",  // 141,
        "Fe56_to_n_Fe55_removed",  // 142,
        "Fe55_to_n_Fe54_removed",  // 143,
        "n_Ni56_to_Ni57_removed",  // 144,
        "n_Ni57_to_Ni58_removed",  // 145,
        "Ni58_to_n_Ni57_removed",  // 146,
        "Ni57_to_n_Ni56_removed",  // 147,
        "Mg23_to_Na23",  // 148,
        "Na22_to_Ne22",  // 149,
        "Na23_to_Mg23",  // 150,
        "Ne22_to_Na22",  // 151,
        "Co56_to_Fe56",  // 152,
        "Co56_to_Ni56",  // 153,
        "Fe56_to_Co56",  // 154,
        "n_to_p",  // 155,
        "Ni56_to_Co56",  // 156,
        "p_to_n",  // 157,
        "S32_He4_to_Ar36_approx",  // 158,
        "Ar36_to_S32_He4_approx",  // 159,
        "Ar36_He4_to_Ca40_approx",  // 160,
        "Ca40_to_Ar36_He4_approx",  // 161,
        "Ca40_He4_to_Ti44_approx",  // 162,
        "Ti44_to_Ca40_He4_approx",  // 163,
        "Ti44_He4_to_Cr48_approx",  // 164,
        "Cr48_to_Ti44_He4_approx",  // 165,
        "Fe52_n_n_to_Fe54_approx",  // 166,
        "Fe54_to_Fe52_n_n_approx",  // 167,
        "Fe54_n_n_to_Fe56_approx",  // 168,
        "Fe56_to_Fe54_n_n_approx",  // 169,
        "Ni56_n_n_to_Ni58_approx",  // 170,
        "Ni58_to_Ni56_n_n_approx"  // 171,
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
