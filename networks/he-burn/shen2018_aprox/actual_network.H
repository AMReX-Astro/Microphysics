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
        else if constexpr (spec == B11) {
            return 76.20505591999972_rt;
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
        else if constexpr (spec == Ne20) {
            return 160.64482384000075_rt;
        }
        else if constexpr (spec == Ne22) {
            return 177.77024384000106_rt;
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
        else if constexpr (spec == Mg26) {
            return 216.68065016799665_rt;
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
        else if constexpr (spec == Si30) {
            return 255.61964649599756_rt;
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
        else if constexpr (spec == Cl35) {
            return 298.20976388800045_rt;
        }
        else if constexpr (spec == Ar36) {
            return 306.7167469519991_rt;
        }
        else if constexpr (spec == Ar37) {
            return 315.50420305199077_rt;
        }
        else if constexpr (spec == Ar38) {
            return 327.3426711519933_rt;
        }
        else if constexpr (spec == Ar39) {
            return 333.9411592519973_rt;
        }
        else if constexpr (spec == K39) {
            return 333.7240072160057_rt;
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
        else if constexpr (spec == Fe56) {
            return 492.2599506639962_rt;
        }
        else if constexpr (spec == Co55) {
            return 476.82912552799826_rt;
        }
        else if constexpr (spec == Ni56) {
            return 483.9956965919919_rt;
        }
        else if constexpr (spec == Ni58) {
            return 506.4596327920008_rt;
        }
        else if constexpr (spec == Ni59) {
            return 515.4588508919987_rt;
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
        k_p_B11_to_C12 = 1,
        k_p_C12_to_N13 = 2,
        k_He4_C12_to_O16 = 3,
        k_n_N13_to_N14 = 4,
        k_n_O16_to_O17 = 5,
        k_He4_O16_to_Ne20 = 6,
        k_He4_Ne20_to_Mg24 = 7,
        k_p_Ne22_to_Na23 = 8,
        k_He4_Ne22_to_Mg26 = 9,
        k_p_Na23_to_Mg24 = 10,
        k_He4_Na23_to_Al27 = 11,
        k_n_Mg24_to_Mg25 = 12,
        k_He4_Mg24_to_Si28 = 13,
        k_n_Mg25_to_Mg26 = 14,
        k_He4_Mg25_to_Si29 = 15,
        k_p_Mg26_to_Al27 = 16,
        k_He4_Mg26_to_Si30 = 17,
        k_p_Al27_to_Si28 = 18,
        k_He4_Al27_to_P31 = 19,
        k_n_Si28_to_Si29 = 20,
        k_He4_Si28_to_S32 = 21,
        k_n_Si29_to_Si30 = 22,
        k_p_Si29_to_P30 = 23,
        k_p_Si30_to_P31 = 24,
        k_n_P30_to_P31 = 25,
        k_p_P30_to_S31 = 26,
        k_p_P31_to_S32 = 27,
        k_He4_P31_to_Cl35 = 28,
        k_n_S31_to_S32 = 29,
        k_He4_S32_to_Ar36 = 30,
        k_p_Cl35_to_Ar36 = 31,
        k_He4_Cl35_to_K39 = 32,
        k_n_Ar36_to_Ar37 = 33,
        k_He4_Ar36_to_Ca40 = 34,
        k_n_Ar37_to_Ar38 = 35,
        k_n_Ar38_to_Ar39 = 36,
        k_p_Ar38_to_K39 = 37,
        k_p_K39_to_Ca40 = 38,
        k_He4_Fe52_to_Ni56 = 39,
        k_p_nse_Co55_to_Ni56 = 40,
        k_n_Ni58_to_Ni59 = 41,
        k_He4_B11_to_n_N14 = 42,
        k_C12_C12_to_p_Na23 = 43,
        k_C12_C12_to_He4_Ne20 = 44,
        k_He4_N13_to_p_O16 = 45,
        k_C12_O16_to_p_Al27 = 46,
        k_C12_O16_to_He4_Mg24 = 47,
        k_O16_O16_to_n_S31 = 48,
        k_O16_O16_to_p_P31 = 49,
        k_O16_O16_to_He4_Si28 = 50,
        k_p_O17_to_He4_N14 = 51,
        k_He4_O17_to_n_Ne20 = 52,
        k_C12_Ne20_to_n_S31 = 53,
        k_C12_Ne20_to_p_P31 = 54,
        k_C12_Ne20_to_He4_Si28 = 55,
        k_He4_Ne22_to_n_Mg25 = 56,
        k_p_Na23_to_He4_Ne20 = 57,
        k_He4_Na23_to_p_Mg26 = 58,
        k_He4_Mg25_to_n_Si28 = 59,
        k_He4_Mg26_to_n_Si29 = 60,
        k_p_Al27_to_He4_Mg24 = 61,
        k_He4_Al27_to_n_P30 = 62,
        k_He4_Al27_to_p_Si30 = 63,
        k_n_P30_to_p_Si30 = 64,
        k_p_P31_to_He4_Si28 = 65,
        k_n_S31_to_p_P31 = 66,
        k_n_S31_to_He4_Si28 = 67,
        k_n_S32_to_He4_Si29 = 68,
        k_p_Cl35_to_He4_S32 = 69,
        k_He4_Cl35_to_p_Ar38 = 70,
        k_n_K39_to_p_Ar39 = 71,
        k_p_K39_to_He4_Ar36 = 72,
        k_n_Ca40_to_He4_Ar37 = 73,
        k_He4_Fe52_to_p_nse_Co55 = 74,
        k_He4_Co55_to_p_nse_Ni58 = 75,
        k_n_Ni59_to_He4_Fe56 = 76,
        k_p_B11_to_He4_He4_He4 = 77,
        k_He4_He4_He4_to_C12 = 78,
        k_He4_Ca40_to_Ti44_removed = 79,
        k_p_Sc43_to_Ti44_removed = 80,
        k_p_Sc43_to_He4_Ca40_removed = 81,
        k_He4_Ti44_to_Cr48_removed = 82,
        k_He4_Ti44_to_p_V47_removed = 83,
        k_p_V47_to_Cr48_removed = 84,
        k_He4_Cr48_to_Fe52_removed = 85,
        k_He4_Cr48_to_p_Mn51_removed = 86,
        k_p_Mn51_to_Fe52_removed = 87,
        k_Ar39_to_K39 = 88,
        k_K39_to_Ar39 = 89,
        k_P30_to_Si30 = 90,
        k_P31_to_S31 = 91,
        k_S31_to_P31 = 92,
        k_Si30_to_P30 = 93,
        k_n_to_p = 94,
        k_p_to_n = 95,
        k_Ca40_He4_to_Ti44_approx = 96,
        k_Ti44_to_Ca40_He4_approx = 97,
        k_Ti44_He4_to_Cr48_approx = 98,
        k_Cr48_to_Ti44_He4_approx = 99,
        k_Cr48_He4_to_Fe52_approx = 100,
        k_Fe52_to_Cr48_He4_approx = 101,
        k_C12_to_p_B11_derived = 102,
        k_N13_to_p_C12_derived = 103,
        k_N14_to_n_N13_derived = 104,
        k_O16_to_He4_C12_derived = 105,
        k_O17_to_n_O16_derived = 106,
        k_Ne20_to_He4_O16_derived = 107,
        k_Na23_to_p_Ne22_derived = 108,
        k_Mg24_to_p_Na23_derived = 109,
        k_Mg24_to_He4_Ne20_derived = 110,
        k_Mg25_to_n_Mg24_derived = 111,
        k_Mg26_to_n_Mg25_derived = 112,
        k_Mg26_to_He4_Ne22_derived = 113,
        k_Al27_to_p_Mg26_derived = 114,
        k_Al27_to_He4_Na23_derived = 115,
        k_Si28_to_p_Al27_derived = 116,
        k_Si28_to_He4_Mg24_derived = 117,
        k_Si29_to_n_Si28_derived = 118,
        k_Si29_to_He4_Mg25_derived = 119,
        k_Si30_to_n_Si29_derived = 120,
        k_Si30_to_He4_Mg26_derived = 121,
        k_P30_to_p_Si29_derived = 122,
        k_P31_to_n_P30_derived = 123,
        k_P31_to_p_Si30_derived = 124,
        k_P31_to_He4_Al27_derived = 125,
        k_S31_to_p_P30_derived = 126,
        k_S32_to_n_S31_derived = 127,
        k_S32_to_p_P31_derived = 128,
        k_S32_to_He4_Si28_derived = 129,
        k_Cl35_to_He4_P31_derived = 130,
        k_Ar36_to_p_Cl35_derived = 131,
        k_Ar36_to_He4_S32_derived = 132,
        k_Ar37_to_n_Ar36_derived = 133,
        k_Ar38_to_n_Ar37_derived = 134,
        k_Ar39_to_n_Ar38_derived = 135,
        k_K39_to_p_Ar38_derived = 136,
        k_K39_to_He4_Cl35_derived = 137,
        k_Ca40_to_p_K39_derived = 138,
        k_Ca40_to_He4_Ar36_derived = 139,
        k_Ni56_to_p_nse_Co55_derived = 140,
        k_Ni56_to_He4_Fe52_derived = 141,
        k_Ni59_to_n_Ni58_derived = 142,
        k_C12_to_He4_He4_He4_derived = 143,
        k_n_N14_to_He4_B11_derived = 144,
        k_He4_N14_to_p_O17_derived = 145,
        k_p_O16_to_He4_N13_derived = 146,
        k_n_Ne20_to_He4_O17_derived = 147,
        k_He4_Ne20_to_p_Na23_derived = 148,
        k_He4_Ne20_to_C12_C12_derived = 149,
        k_p_Na23_to_C12_C12_derived = 150,
        k_He4_Mg24_to_p_Al27_derived = 151,
        k_He4_Mg24_to_C12_O16_derived = 152,
        k_n_Mg25_to_He4_Ne22_derived = 153,
        k_p_Mg26_to_He4_Na23_derived = 154,
        k_p_Al27_to_C12_O16_derived = 155,
        k_n_Si28_to_He4_Mg25_derived = 156,
        k_He4_Si28_to_n_S31_derived = 157,
        k_He4_Si28_to_p_P31_derived = 158,
        k_He4_Si28_to_C12_Ne20_derived = 159,
        k_He4_Si28_to_O16_O16_derived = 160,
        k_n_Si29_to_He4_Mg26_derived = 161,
        k_He4_Si29_to_n_S32_derived = 162,
        k_p_Si30_to_n_P30_derived = 163,
        k_p_Si30_to_He4_Al27_derived = 164,
        k_n_P30_to_He4_Al27_derived = 165,
        k_p_P31_to_n_S31_derived = 166,
        k_p_P31_to_C12_Ne20_derived = 167,
        k_p_P31_to_O16_O16_derived = 168,
        k_n_S31_to_C12_Ne20_derived = 169,
        k_n_S31_to_O16_O16_derived = 170,
        k_He4_S32_to_p_Cl35_derived = 171,
        k_He4_Ar36_to_p_K39_derived = 172,
        k_He4_Ar37_to_n_Ca40_derived = 173,
        k_p_Ar38_to_He4_Cl35_derived = 174,
        k_p_Ar39_to_n_K39_derived = 175,
        k_He4_Fe56_to_n_Ni59_derived = 176,
        k_p_nse_Co55_to_He4_Fe52_derived = 177,
        k_p_nse_Ni58_to_He4_Co55_derived = 178,
        k_He4_He4_He4_to_p_B11_derived = 179,
        k_He4_Ca40_to_p_Sc43_derived_removed = 180,
        k_Ti44_to_He4_Ca40_derived_removed = 181,
        k_Ti44_to_p_Sc43_derived_removed = 182,
        k_Cr48_to_He4_Ti44_derived_removed = 183,
        k_Cr48_to_p_V47_derived_removed = 184,
        k_p_V47_to_He4_Ti44_derived_removed = 185,
        k_Fe52_to_He4_Cr48_derived_removed = 186,
        k_Fe52_to_p_Mn51_derived_removed = 187,
        k_p_Mn51_to_He4_Cr48_derived_removed = 188,
        NumRates = k_p_Mn51_to_He4_Cr48_derived_removed
    };

    // number of reaclib rates

    const int NrateReaclib = 174;

    // number of tabular rates

    const int NrateTabular = 8;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "p_B11_to_C12",  // 1,
        "p_C12_to_N13",  // 2,
        "He4_C12_to_O16",  // 3,
        "n_N13_to_N14",  // 4,
        "n_O16_to_O17",  // 5,
        "He4_O16_to_Ne20",  // 6,
        "He4_Ne20_to_Mg24",  // 7,
        "p_Ne22_to_Na23",  // 8,
        "He4_Ne22_to_Mg26",  // 9,
        "p_Na23_to_Mg24",  // 10,
        "He4_Na23_to_Al27",  // 11,
        "n_Mg24_to_Mg25",  // 12,
        "He4_Mg24_to_Si28",  // 13,
        "n_Mg25_to_Mg26",  // 14,
        "He4_Mg25_to_Si29",  // 15,
        "p_Mg26_to_Al27",  // 16,
        "He4_Mg26_to_Si30",  // 17,
        "p_Al27_to_Si28",  // 18,
        "He4_Al27_to_P31",  // 19,
        "n_Si28_to_Si29",  // 20,
        "He4_Si28_to_S32",  // 21,
        "n_Si29_to_Si30",  // 22,
        "p_Si29_to_P30",  // 23,
        "p_Si30_to_P31",  // 24,
        "n_P30_to_P31",  // 25,
        "p_P30_to_S31",  // 26,
        "p_P31_to_S32",  // 27,
        "He4_P31_to_Cl35",  // 28,
        "n_S31_to_S32",  // 29,
        "He4_S32_to_Ar36",  // 30,
        "p_Cl35_to_Ar36",  // 31,
        "He4_Cl35_to_K39",  // 32,
        "n_Ar36_to_Ar37",  // 33,
        "He4_Ar36_to_Ca40",  // 34,
        "n_Ar37_to_Ar38",  // 35,
        "n_Ar38_to_Ar39",  // 36,
        "p_Ar38_to_K39",  // 37,
        "p_K39_to_Ca40",  // 38,
        "He4_Fe52_to_Ni56",  // 39,
        "p_nse_Co55_to_Ni56",  // 40,
        "n_Ni58_to_Ni59",  // 41,
        "He4_B11_to_n_N14",  // 42,
        "C12_C12_to_p_Na23",  // 43,
        "C12_C12_to_He4_Ne20",  // 44,
        "He4_N13_to_p_O16",  // 45,
        "C12_O16_to_p_Al27",  // 46,
        "C12_O16_to_He4_Mg24",  // 47,
        "O16_O16_to_n_S31",  // 48,
        "O16_O16_to_p_P31",  // 49,
        "O16_O16_to_He4_Si28",  // 50,
        "p_O17_to_He4_N14",  // 51,
        "He4_O17_to_n_Ne20",  // 52,
        "C12_Ne20_to_n_S31",  // 53,
        "C12_Ne20_to_p_P31",  // 54,
        "C12_Ne20_to_He4_Si28",  // 55,
        "He4_Ne22_to_n_Mg25",  // 56,
        "p_Na23_to_He4_Ne20",  // 57,
        "He4_Na23_to_p_Mg26",  // 58,
        "He4_Mg25_to_n_Si28",  // 59,
        "He4_Mg26_to_n_Si29",  // 60,
        "p_Al27_to_He4_Mg24",  // 61,
        "He4_Al27_to_n_P30",  // 62,
        "He4_Al27_to_p_Si30",  // 63,
        "n_P30_to_p_Si30",  // 64,
        "p_P31_to_He4_Si28",  // 65,
        "n_S31_to_p_P31",  // 66,
        "n_S31_to_He4_Si28",  // 67,
        "n_S32_to_He4_Si29",  // 68,
        "p_Cl35_to_He4_S32",  // 69,
        "He4_Cl35_to_p_Ar38",  // 70,
        "n_K39_to_p_Ar39",  // 71,
        "p_K39_to_He4_Ar36",  // 72,
        "n_Ca40_to_He4_Ar37",  // 73,
        "He4_Fe52_to_p_nse_Co55",  // 74,
        "He4_Co55_to_p_nse_Ni58",  // 75,
        "n_Ni59_to_He4_Fe56",  // 76,
        "p_B11_to_He4_He4_He4",  // 77,
        "He4_He4_He4_to_C12",  // 78,
        "He4_Ca40_to_Ti44_removed",  // 79,
        "p_Sc43_to_Ti44_removed",  // 80,
        "p_Sc43_to_He4_Ca40_removed",  // 81,
        "He4_Ti44_to_Cr48_removed",  // 82,
        "He4_Ti44_to_p_V47_removed",  // 83,
        "p_V47_to_Cr48_removed",  // 84,
        "He4_Cr48_to_Fe52_removed",  // 85,
        "He4_Cr48_to_p_Mn51_removed",  // 86,
        "p_Mn51_to_Fe52_removed",  // 87,
        "Ar39_to_K39",  // 88,
        "K39_to_Ar39",  // 89,
        "P30_to_Si30",  // 90,
        "P31_to_S31",  // 91,
        "S31_to_P31",  // 92,
        "Si30_to_P30",  // 93,
        "n_to_p",  // 94,
        "p_to_n",  // 95,
        "Ca40_He4_to_Ti44_approx",  // 96,
        "Ti44_to_Ca40_He4_approx",  // 97,
        "Ti44_He4_to_Cr48_approx",  // 98,
        "Cr48_to_Ti44_He4_approx",  // 99,
        "Cr48_He4_to_Fe52_approx",  // 100,
        "Fe52_to_Cr48_He4_approx",  // 101,
        "C12_to_p_B11_derived",  // 102,
        "N13_to_p_C12_derived",  // 103,
        "N14_to_n_N13_derived",  // 104,
        "O16_to_He4_C12_derived",  // 105,
        "O17_to_n_O16_derived",  // 106,
        "Ne20_to_He4_O16_derived",  // 107,
        "Na23_to_p_Ne22_derived",  // 108,
        "Mg24_to_p_Na23_derived",  // 109,
        "Mg24_to_He4_Ne20_derived",  // 110,
        "Mg25_to_n_Mg24_derived",  // 111,
        "Mg26_to_n_Mg25_derived",  // 112,
        "Mg26_to_He4_Ne22_derived",  // 113,
        "Al27_to_p_Mg26_derived",  // 114,
        "Al27_to_He4_Na23_derived",  // 115,
        "Si28_to_p_Al27_derived",  // 116,
        "Si28_to_He4_Mg24_derived",  // 117,
        "Si29_to_n_Si28_derived",  // 118,
        "Si29_to_He4_Mg25_derived",  // 119,
        "Si30_to_n_Si29_derived",  // 120,
        "Si30_to_He4_Mg26_derived",  // 121,
        "P30_to_p_Si29_derived",  // 122,
        "P31_to_n_P30_derived",  // 123,
        "P31_to_p_Si30_derived",  // 124,
        "P31_to_He4_Al27_derived",  // 125,
        "S31_to_p_P30_derived",  // 126,
        "S32_to_n_S31_derived",  // 127,
        "S32_to_p_P31_derived",  // 128,
        "S32_to_He4_Si28_derived",  // 129,
        "Cl35_to_He4_P31_derived",  // 130,
        "Ar36_to_p_Cl35_derived",  // 131,
        "Ar36_to_He4_S32_derived",  // 132,
        "Ar37_to_n_Ar36_derived",  // 133,
        "Ar38_to_n_Ar37_derived",  // 134,
        "Ar39_to_n_Ar38_derived",  // 135,
        "K39_to_p_Ar38_derived",  // 136,
        "K39_to_He4_Cl35_derived",  // 137,
        "Ca40_to_p_K39_derived",  // 138,
        "Ca40_to_He4_Ar36_derived",  // 139,
        "Ni56_to_p_nse_Co55_derived",  // 140,
        "Ni56_to_He4_Fe52_derived",  // 141,
        "Ni59_to_n_Ni58_derived",  // 142,
        "C12_to_He4_He4_He4_derived",  // 143,
        "n_N14_to_He4_B11_derived",  // 144,
        "He4_N14_to_p_O17_derived",  // 145,
        "p_O16_to_He4_N13_derived",  // 146,
        "n_Ne20_to_He4_O17_derived",  // 147,
        "He4_Ne20_to_p_Na23_derived",  // 148,
        "He4_Ne20_to_C12_C12_derived",  // 149,
        "p_Na23_to_C12_C12_derived",  // 150,
        "He4_Mg24_to_p_Al27_derived",  // 151,
        "He4_Mg24_to_C12_O16_derived",  // 152,
        "n_Mg25_to_He4_Ne22_derived",  // 153,
        "p_Mg26_to_He4_Na23_derived",  // 154,
        "p_Al27_to_C12_O16_derived",  // 155,
        "n_Si28_to_He4_Mg25_derived",  // 156,
        "He4_Si28_to_n_S31_derived",  // 157,
        "He4_Si28_to_p_P31_derived",  // 158,
        "He4_Si28_to_C12_Ne20_derived",  // 159,
        "He4_Si28_to_O16_O16_derived",  // 160,
        "n_Si29_to_He4_Mg26_derived",  // 161,
        "He4_Si29_to_n_S32_derived",  // 162,
        "p_Si30_to_n_P30_derived",  // 163,
        "p_Si30_to_He4_Al27_derived",  // 164,
        "n_P30_to_He4_Al27_derived",  // 165,
        "p_P31_to_n_S31_derived",  // 166,
        "p_P31_to_C12_Ne20_derived",  // 167,
        "p_P31_to_O16_O16_derived",  // 168,
        "n_S31_to_C12_Ne20_derived",  // 169,
        "n_S31_to_O16_O16_derived",  // 170,
        "He4_S32_to_p_Cl35_derived",  // 171,
        "He4_Ar36_to_p_K39_derived",  // 172,
        "He4_Ar37_to_n_Ca40_derived",  // 173,
        "p_Ar38_to_He4_Cl35_derived",  // 174,
        "p_Ar39_to_n_K39_derived",  // 175,
        "He4_Fe56_to_n_Ni59_derived",  // 176,
        "p_nse_Co55_to_He4_Fe52_derived",  // 177,
        "p_nse_Ni58_to_He4_Co55_derived",  // 178,
        "He4_He4_He4_to_p_B11_derived",  // 179,
        "He4_Ca40_to_p_Sc43_derived_removed",  // 180,
        "Ti44_to_He4_Ca40_derived_removed",  // 181,
        "Ti44_to_p_Sc43_derived_removed",  // 182,
        "Cr48_to_He4_Ti44_derived_removed",  // 183,
        "Cr48_to_p_V47_derived_removed",  // 184,
        "p_V47_to_He4_Ti44_derived_removed",  // 185,
        "Fe52_to_He4_Cr48_derived_removed",  // 186,
        "Fe52_to_p_Mn51_derived_removed",  // 187,
        "p_Mn51_to_He4_Cr48_derived_removed"  // 188,
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
