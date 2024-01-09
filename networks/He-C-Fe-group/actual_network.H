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
    constexpr Real bion () {
        using namespace Species;

        static_assert(spec >= 1 && spec <= NumSpec);

        // Set the binding energy of the element

        if constexpr (spec == N) {
            return 0.0_rt;
        }
        else if constexpr (spec == H1) {
            return 0.0_rt;
        }
        else if constexpr (spec == He4) {
            return 28.29566_rt;
        }
        else if constexpr (spec == C12) {
            return 92.16172800000001_rt;
        }
        else if constexpr (spec == N13) {
            return 94.105219_rt;
        }
        else if constexpr (spec == N14) {
            return 104.65859599999999_rt;
        }
        else if constexpr (spec == O16) {
            return 127.619296_rt;
        }
        else if constexpr (spec == F18) {
            return 137.369484_rt;
        }
        else if constexpr (spec == Ne20) {
            return 160.6448_rt;
        }
        else if constexpr (spec == Ne21) {
            return 167.405973_rt;
        }
        else if constexpr (spec == Na22) {
            return 174.144674_rt;
        }
        else if constexpr (spec == Na23) {
            return 186.56433900000002_rt;
        }
        else if constexpr (spec == Mg24) {
            return 198.25701600000002_rt;
        }
        else if constexpr (spec == Al27) {
            return 224.951931_rt;
        }
        else if constexpr (spec == Si28) {
            return 236.536832_rt;
        }
        else if constexpr (spec == P31) {
            return 262.91617699999995_rt;
        }
        else if constexpr (spec == S32) {
            return 271.78012800000005_rt;
        }
        else if constexpr (spec == Ar36) {
            return 306.716724_rt;
        }
        else if constexpr (spec == Ca40) {
            return 342.05212000000006_rt;
        }
        else if constexpr (spec == Ti44) {
            return 375.47488000000004_rt;
        }
        else if constexpr (spec == Cr48) {
            return 411.46891200000005_rt;
        }
        else if constexpr (spec == Mn51) {
            return 440.32237200000003_rt;
        }
        else if constexpr (spec == Mn55) {
            return 482.07621_rt;
        }
        else if constexpr (spec == Fe52) {
            return 447.697848_rt;
        }
        else if constexpr (spec == Fe53) {
            return 458.386347_rt;
        }
        else if constexpr (spec == Fe54) {
            return 471.76462799999996_rt;
        }
        else if constexpr (spec == Fe55) {
            return 481.06272499999994_rt;
        }
        else if constexpr (spec == Fe56) {
            return 492.2598239999999_rt;
        }
        else if constexpr (spec == Co55) {
            return 476.82899_rt;
        }
        else if constexpr (spec == Co56) {
            return 486.91081599999995_rt;
        }
        else if constexpr (spec == Co57) {
            return 498.287274_rt;
        }
        else if constexpr (spec == Ni56) {
            return 483.995624_rt;
        }
        else if constexpr (spec == Ni57) {
            return 494.2431810000001_rt;
        }
        else if constexpr (spec == Ni58) {
            return 506.45942199999996_rt;
        }
        else if constexpr (spec == Cu59) {
            return 509.878_rt;
        }
        else if constexpr (spec == Zn60) {
            return 514.983_rt;
        }


        // Return zero if we don't recognize the species.
        return 0.0_rt;
    }

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr Real mion () {
        static_assert(spec >= 1 && spec <= NumSpec);

        constexpr Real A = NetworkProperties::aion(spec);
        constexpr Real Z = NetworkProperties::zion(spec);

        return (A - Z) * C::Legacy::m_n + Z * (C::Legacy::m_p + C::Legacy::m_e) - bion<spec>() * C::Legacy::MeV2gr;
    }

    // Legacy (non-templated) interfaces

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real bion (int spec) {
        using namespace Species;

        Real b = 0.0_rt;

        // Set the binding energy of the element
        constexpr_for<1, NumSpec+1>([&] (auto n) {
            if (n == spec) {
                b = bion<n>();
            }
        });

        return b;
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real mion (int spec) {
        using namespace Species;

        Real m = 0.0_rt;

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
        k_p_C12_to_N13 = 1,
        k_He4_C12_to_O16 = 2,
        k_He4_N14_to_F18 = 3,
        k_He4_O16_to_Ne20 = 4,
        k_He4_F18_to_Na22 = 5,
        k_He4_Ne20_to_Mg24 = 6,
        k_p_Ne21_to_Na22 = 7,
        k_p_Na23_to_Mg24 = 8,
        k_He4_Mg24_to_Si28 = 9,
        k_p_Al27_to_Si28 = 10,
        k_He4_Si28_to_S32 = 11,
        k_p_P31_to_S32 = 12,
        k_He4_Cr48_to_Fe52 = 13,
        k_p_Mn51_to_Fe52 = 14,
        k_He4_Mn51_to_Co55 = 15,
        k_He4_Fe52_to_Ni56 = 16,
        k_p_Co55_to_Ni56 = 17,
        k_He4_Co55_to_Cu59 = 18,
        k_He4_Ni56_to_Zn60 = 19,
        k_p_Cu59_to_Zn60 = 20,
        k_C12_C12_to_p_Na23 = 21,
        k_C12_C12_to_He4_Ne20 = 22,
        k_He4_N13_to_p_O16 = 23,
        k_C12_O16_to_p_Al27 = 24,
        k_C12_O16_to_He4_Mg24 = 25,
        k_O16_O16_to_p_P31 = 26,
        k_O16_O16_to_He4_Si28 = 27,
        k_He4_F18_to_p_Ne21 = 28,
        k_p_Na23_to_He4_Ne20 = 29,
        k_p_Al27_to_He4_Mg24 = 30,
        k_p_P31_to_He4_Si28 = 31,
        k_He4_Cr48_to_p_Mn51 = 32,
        k_He4_Fe52_to_p_Co55 = 33,
        k_p_Cu59_to_He4_Ni56 = 34,
        k_He4_He4_He4_to_C12 = 35,
        k_C12_C12_to_Mg24_modified = 36,
        k_O16_O16_to_S32_modified = 37,
        k_C12_O16_to_Si28_modified = 38,
        k_p_Mn55_to_Fe56 = 39,
        k_n_Fe52_to_Fe53 = 40,
        k_n_Fe53_to_Fe54 = 41,
        k_He4_Fe53_to_Ni57 = 42,
        k_n_Fe54_to_Fe55 = 43,
        k_p_Fe54_to_Co55 = 44,
        k_He4_Fe54_to_Ni58 = 45,
        k_n_Fe55_to_Fe56 = 46,
        k_p_Fe55_to_Co56 = 47,
        k_p_Fe56_to_Co57 = 48,
        k_n_Co55_to_Co56 = 49,
        k_n_Co56_to_Co57 = 50,
        k_p_Co56_to_Ni57 = 51,
        k_p_Co57_to_Ni58 = 52,
        k_n_Ni56_to_Ni57 = 53,
        k_n_Ni57_to_Ni58 = 54,
        k_p_Ni58_to_Cu59 = 55,
        k_He4_Mn51_to_p_Fe54 = 56,
        k_He4_Fe53_to_p_Co56 = 57,
        k_n_Fe55_to_p_Mn55 = 58,
        k_n_Co55_to_p_Fe55 = 59,
        k_He4_Co55_to_p_Ni58 = 60,
        k_n_Co56_to_p_Fe56 = 61,
        k_p_Co57_to_He4_Fe54 = 62,
        k_n_Ni56_to_p_Co56 = 63,
        k_n_Ni56_to_He4_Fe53 = 64,
        k_n_Ni57_to_p_Co57 = 65,
        k_n_Ni57_to_He4_Fe54 = 66,
        k_n_Ni58_to_He4_Fe55 = 67,
        k_n_Cu59_to_He4_Co56 = 68,
        k_n_Zn60_to_He4_Ni57 = 69,
        k_He4_S32_to_Ar36_removed = 70,
        k_p_Cl35_to_Ar36_removed = 71,
        k_p_Cl35_to_He4_S32_removed = 72,
        k_He4_Ar36_to_Ca40_removed = 73,
        k_p_K39_to_Ca40_removed = 74,
        k_p_K39_to_He4_Ar36_removed = 75,
        k_He4_Ca40_to_Ti44_removed = 76,
        k_p_Sc43_to_Ti44_removed = 77,
        k_p_Sc43_to_He4_Ca40_removed = 78,
        k_He4_Ti44_to_Cr48_removed = 79,
        k_He4_Ti44_to_p_V47_removed = 80,
        k_p_V47_to_Cr48_removed = 81,
        k_Co55_to_Fe55 = 82,
        k_Co56_to_Fe56 = 83,
        k_Co56_to_Ni56 = 84,
        k_Co57_to_Ni57 = 85,
        k_Fe55_to_Co55 = 86,
        k_Fe55_to_Mn55 = 87,
        k_Fe56_to_Co56 = 88,
        k_Mn55_to_Fe55 = 89,
        k_n_to_p = 90,
        k_Ni56_to_Co56 = 91,
        k_Ni57_to_Co57 = 92,
        k_p_to_n = 93,
        k_S32_He4_to_Ar36_approx = 94,
        k_Ar36_to_S32_He4_approx = 95,
        k_Ar36_He4_to_Ca40_approx = 96,
        k_Ca40_to_Ar36_He4_approx = 97,
        k_Ca40_He4_to_Ti44_approx = 98,
        k_Ti44_to_Ca40_He4_approx = 99,
        k_Ti44_He4_to_Cr48_approx = 100,
        k_Cr48_to_Ti44_He4_approx = 101,
        k_N13_to_p_C12_derived = 102,
        k_O16_to_He4_C12_derived = 103,
        k_F18_to_He4_N14_derived = 104,
        k_Ne20_to_He4_O16_derived = 105,
        k_Na22_to_p_Ne21_derived = 106,
        k_Na22_to_He4_F18_derived = 107,
        k_Mg24_to_p_Na23_derived = 108,
        k_Mg24_to_He4_Ne20_derived = 109,
        k_Si28_to_p_Al27_derived = 110,
        k_Si28_to_He4_Mg24_derived = 111,
        k_S32_to_p_P31_derived = 112,
        k_S32_to_He4_Si28_derived = 113,
        k_Fe52_to_p_Mn51_derived = 114,
        k_Fe52_to_He4_Cr48_derived = 115,
        k_Co55_to_He4_Mn51_derived = 116,
        k_Ni56_to_p_Co55_derived = 117,
        k_Ni56_to_He4_Fe52_derived = 118,
        k_Cu59_to_He4_Co55_derived = 119,
        k_Zn60_to_p_Cu59_derived = 120,
        k_Zn60_to_He4_Ni56_derived = 121,
        k_C12_to_He4_He4_He4_derived = 122,
        k_p_O16_to_He4_N13_derived = 123,
        k_He4_Ne20_to_p_Na23_derived = 124,
        k_p_Ne21_to_He4_F18_derived = 125,
        k_He4_Mg24_to_p_Al27_derived = 126,
        k_He4_Si28_to_p_P31_derived = 127,
        k_p_Mn51_to_He4_Cr48_derived = 128,
        k_p_Co55_to_He4_Fe52_derived = 129,
        k_He4_Ni56_to_p_Cu59_derived = 130,
        k_Fe53_to_n_Fe52_derived = 131,
        k_Fe54_to_n_Fe53_derived = 132,
        k_Fe55_to_n_Fe54_derived = 133,
        k_Fe56_to_n_Fe55_derived = 134,
        k_Fe56_to_p_Mn55_derived = 135,
        k_Co55_to_p_Fe54_derived = 136,
        k_Co56_to_n_Co55_derived = 137,
        k_Co56_to_p_Fe55_derived = 138,
        k_Co57_to_n_Co56_derived = 139,
        k_Co57_to_p_Fe56_derived = 140,
        k_Ni57_to_n_Ni56_derived = 141,
        k_Ni57_to_p_Co56_derived = 142,
        k_Ni57_to_He4_Fe53_derived = 143,
        k_Ni58_to_n_Ni57_derived = 144,
        k_Ni58_to_p_Co57_derived = 145,
        k_Ni58_to_He4_Fe54_derived = 146,
        k_Cu59_to_p_Ni58_derived = 147,
        k_p_Mn55_to_n_Fe55_derived = 148,
        k_He4_Fe53_to_n_Ni56_derived = 149,
        k_p_Fe54_to_He4_Mn51_derived = 150,
        k_He4_Fe54_to_n_Ni57_derived = 151,
        k_He4_Fe54_to_p_Co57_derived = 152,
        k_p_Fe55_to_n_Co55_derived = 153,
        k_He4_Fe55_to_n_Ni58_derived = 154,
        k_p_Fe56_to_n_Co56_derived = 155,
        k_p_Co56_to_n_Ni56_derived = 156,
        k_p_Co56_to_He4_Fe53_derived = 157,
        k_He4_Co56_to_n_Cu59_derived = 158,
        k_p_Co57_to_n_Ni57_derived = 159,
        k_He4_Ni57_to_n_Zn60_derived = 160,
        k_p_Ni58_to_He4_Co55_derived = 161,
        k_He4_S32_to_p_Cl35_derived_removed = 162,
        k_Ar36_to_He4_S32_derived_removed = 163,
        k_Ar36_to_p_Cl35_derived_removed = 164,
        k_He4_Ar36_to_p_K39_derived_removed = 165,
        k_Ca40_to_He4_Ar36_derived_removed = 166,
        k_Ca40_to_p_K39_derived_removed = 167,
        k_He4_Ca40_to_p_Sc43_derived_removed = 168,
        k_Ti44_to_He4_Ca40_derived_removed = 169,
        k_Ti44_to_p_Sc43_derived_removed = 170,
        k_Cr48_to_He4_Ti44_derived_removed = 171,
        k_Cr48_to_p_V47_derived_removed = 172,
        k_p_V47_to_He4_Ti44_derived_removed = 173,
        NumRates = k_p_V47_to_He4_Ti44_derived_removed
    };

    // number of reaclib rates

    const int NrateReaclib = 153;

    // number of tabular rates

    const int NrateTabular = 12;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "p_C12_to_N13",  // 1,
        "He4_C12_to_O16",  // 2,
        "He4_N14_to_F18",  // 3,
        "He4_O16_to_Ne20",  // 4,
        "He4_F18_to_Na22",  // 5,
        "He4_Ne20_to_Mg24",  // 6,
        "p_Ne21_to_Na22",  // 7,
        "p_Na23_to_Mg24",  // 8,
        "He4_Mg24_to_Si28",  // 9,
        "p_Al27_to_Si28",  // 10,
        "He4_Si28_to_S32",  // 11,
        "p_P31_to_S32",  // 12,
        "He4_Cr48_to_Fe52",  // 13,
        "p_Mn51_to_Fe52",  // 14,
        "He4_Mn51_to_Co55",  // 15,
        "He4_Fe52_to_Ni56",  // 16,
        "p_Co55_to_Ni56",  // 17,
        "He4_Co55_to_Cu59",  // 18,
        "He4_Ni56_to_Zn60",  // 19,
        "p_Cu59_to_Zn60",  // 20,
        "C12_C12_to_p_Na23",  // 21,
        "C12_C12_to_He4_Ne20",  // 22,
        "He4_N13_to_p_O16",  // 23,
        "C12_O16_to_p_Al27",  // 24,
        "C12_O16_to_He4_Mg24",  // 25,
        "O16_O16_to_p_P31",  // 26,
        "O16_O16_to_He4_Si28",  // 27,
        "He4_F18_to_p_Ne21",  // 28,
        "p_Na23_to_He4_Ne20",  // 29,
        "p_Al27_to_He4_Mg24",  // 30,
        "p_P31_to_He4_Si28",  // 31,
        "He4_Cr48_to_p_Mn51",  // 32,
        "He4_Fe52_to_p_Co55",  // 33,
        "p_Cu59_to_He4_Ni56",  // 34,
        "He4_He4_He4_to_C12",  // 35,
        "C12_C12_to_Mg24_modified",  // 36,
        "O16_O16_to_S32_modified",  // 37,
        "C12_O16_to_Si28_modified",  // 38,
        "p_Mn55_to_Fe56",  // 39,
        "n_Fe52_to_Fe53",  // 40,
        "n_Fe53_to_Fe54",  // 41,
        "He4_Fe53_to_Ni57",  // 42,
        "n_Fe54_to_Fe55",  // 43,
        "p_Fe54_to_Co55",  // 44,
        "He4_Fe54_to_Ni58",  // 45,
        "n_Fe55_to_Fe56",  // 46,
        "p_Fe55_to_Co56",  // 47,
        "p_Fe56_to_Co57",  // 48,
        "n_Co55_to_Co56",  // 49,
        "n_Co56_to_Co57",  // 50,
        "p_Co56_to_Ni57",  // 51,
        "p_Co57_to_Ni58",  // 52,
        "n_Ni56_to_Ni57",  // 53,
        "n_Ni57_to_Ni58",  // 54,
        "p_Ni58_to_Cu59",  // 55,
        "He4_Mn51_to_p_Fe54",  // 56,
        "He4_Fe53_to_p_Co56",  // 57,
        "n_Fe55_to_p_Mn55",  // 58,
        "n_Co55_to_p_Fe55",  // 59,
        "He4_Co55_to_p_Ni58",  // 60,
        "n_Co56_to_p_Fe56",  // 61,
        "p_Co57_to_He4_Fe54",  // 62,
        "n_Ni56_to_p_Co56",  // 63,
        "n_Ni56_to_He4_Fe53",  // 64,
        "n_Ni57_to_p_Co57",  // 65,
        "n_Ni57_to_He4_Fe54",  // 66,
        "n_Ni58_to_He4_Fe55",  // 67,
        "n_Cu59_to_He4_Co56",  // 68,
        "n_Zn60_to_He4_Ni57",  // 69,
        "He4_S32_to_Ar36_removed",  // 70,
        "p_Cl35_to_Ar36_removed",  // 71,
        "p_Cl35_to_He4_S32_removed",  // 72,
        "He4_Ar36_to_Ca40_removed",  // 73,
        "p_K39_to_Ca40_removed",  // 74,
        "p_K39_to_He4_Ar36_removed",  // 75,
        "He4_Ca40_to_Ti44_removed",  // 76,
        "p_Sc43_to_Ti44_removed",  // 77,
        "p_Sc43_to_He4_Ca40_removed",  // 78,
        "He4_Ti44_to_Cr48_removed",  // 79,
        "He4_Ti44_to_p_V47_removed",  // 80,
        "p_V47_to_Cr48_removed",  // 81,
        "Co55_to_Fe55",  // 82,
        "Co56_to_Fe56",  // 83,
        "Co56_to_Ni56",  // 84,
        "Co57_to_Ni57",  // 85,
        "Fe55_to_Co55",  // 86,
        "Fe55_to_Mn55",  // 87,
        "Fe56_to_Co56",  // 88,
        "Mn55_to_Fe55",  // 89,
        "n_to_p",  // 90,
        "Ni56_to_Co56",  // 91,
        "Ni57_to_Co57",  // 92,
        "p_to_n",  // 93,
        "S32_He4_to_Ar36_approx",  // 94,
        "Ar36_to_S32_He4_approx",  // 95,
        "Ar36_He4_to_Ca40_approx",  // 96,
        "Ca40_to_Ar36_He4_approx",  // 97,
        "Ca40_He4_to_Ti44_approx",  // 98,
        "Ti44_to_Ca40_He4_approx",  // 99,
        "Ti44_He4_to_Cr48_approx",  // 100,
        "Cr48_to_Ti44_He4_approx",  // 101,
        "N13_to_p_C12_derived",  // 102,
        "O16_to_He4_C12_derived",  // 103,
        "F18_to_He4_N14_derived",  // 104,
        "Ne20_to_He4_O16_derived",  // 105,
        "Na22_to_p_Ne21_derived",  // 106,
        "Na22_to_He4_F18_derived",  // 107,
        "Mg24_to_p_Na23_derived",  // 108,
        "Mg24_to_He4_Ne20_derived",  // 109,
        "Si28_to_p_Al27_derived",  // 110,
        "Si28_to_He4_Mg24_derived",  // 111,
        "S32_to_p_P31_derived",  // 112,
        "S32_to_He4_Si28_derived",  // 113,
        "Fe52_to_p_Mn51_derived",  // 114,
        "Fe52_to_He4_Cr48_derived",  // 115,
        "Co55_to_He4_Mn51_derived",  // 116,
        "Ni56_to_p_Co55_derived",  // 117,
        "Ni56_to_He4_Fe52_derived",  // 118,
        "Cu59_to_He4_Co55_derived",  // 119,
        "Zn60_to_p_Cu59_derived",  // 120,
        "Zn60_to_He4_Ni56_derived",  // 121,
        "C12_to_He4_He4_He4_derived",  // 122,
        "p_O16_to_He4_N13_derived",  // 123,
        "He4_Ne20_to_p_Na23_derived",  // 124,
        "p_Ne21_to_He4_F18_derived",  // 125,
        "He4_Mg24_to_p_Al27_derived",  // 126,
        "He4_Si28_to_p_P31_derived",  // 127,
        "p_Mn51_to_He4_Cr48_derived",  // 128,
        "p_Co55_to_He4_Fe52_derived",  // 129,
        "He4_Ni56_to_p_Cu59_derived",  // 130,
        "Fe53_to_n_Fe52_derived",  // 131,
        "Fe54_to_n_Fe53_derived",  // 132,
        "Fe55_to_n_Fe54_derived",  // 133,
        "Fe56_to_n_Fe55_derived",  // 134,
        "Fe56_to_p_Mn55_derived",  // 135,
        "Co55_to_p_Fe54_derived",  // 136,
        "Co56_to_n_Co55_derived",  // 137,
        "Co56_to_p_Fe55_derived",  // 138,
        "Co57_to_n_Co56_derived",  // 139,
        "Co57_to_p_Fe56_derived",  // 140,
        "Ni57_to_n_Ni56_derived",  // 141,
        "Ni57_to_p_Co56_derived",  // 142,
        "Ni57_to_He4_Fe53_derived",  // 143,
        "Ni58_to_n_Ni57_derived",  // 144,
        "Ni58_to_p_Co57_derived",  // 145,
        "Ni58_to_He4_Fe54_derived",  // 146,
        "Cu59_to_p_Ni58_derived",  // 147,
        "p_Mn55_to_n_Fe55_derived",  // 148,
        "He4_Fe53_to_n_Ni56_derived",  // 149,
        "p_Fe54_to_He4_Mn51_derived",  // 150,
        "He4_Fe54_to_n_Ni57_derived",  // 151,
        "He4_Fe54_to_p_Co57_derived",  // 152,
        "p_Fe55_to_n_Co55_derived",  // 153,
        "He4_Fe55_to_n_Ni58_derived",  // 154,
        "p_Fe56_to_n_Co56_derived",  // 155,
        "p_Co56_to_n_Ni56_derived",  // 156,
        "p_Co56_to_He4_Fe53_derived",  // 157,
        "He4_Co56_to_n_Cu59_derived",  // 158,
        "p_Co57_to_n_Ni57_derived",  // 159,
        "He4_Ni57_to_n_Zn60_derived",  // 160,
        "p_Ni58_to_He4_Co55_derived",  // 161,
        "He4_S32_to_p_Cl35_derived_removed",  // 162,
        "Ar36_to_He4_S32_derived_removed",  // 163,
        "Ar36_to_p_Cl35_derived_removed",  // 164,
        "He4_Ar36_to_p_K39_derived_removed",  // 165,
        "Ca40_to_He4_Ar36_derived_removed",  // 166,
        "Ca40_to_p_K39_derived_removed",  // 167,
        "He4_Ca40_to_p_Sc43_derived_removed",  // 168,
        "Ti44_to_He4_Ca40_derived_removed",  // 169,
        "Ti44_to_p_Sc43_derived_removed",  // 170,
        "Cr48_to_He4_Ti44_derived_removed",  // 171,
        "Cr48_to_p_V47_derived_removed",  // 172,
        "p_V47_to_He4_Ti44_derived_removed"  // 173,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int H1_index = 1;
    constexpr int N_index = 0;
    constexpr int He4_index = 2;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
