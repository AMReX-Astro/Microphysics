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
        else if constexpr (spec == Sc43) {
            return 366.82549054399715_rt;
        }
        else if constexpr (spec == Ti44) {
            return 375.47496160800074_rt;
        }
        else if constexpr (spec == V47) {
            return 403.3650388719979_rt;
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
        k_He4_K39_to_Sc43 = 39,
        k_He4_Ca40_to_Ti44 = 40,
        k_p_Sc43_to_Ti44 = 41,
        k_He4_Sc43_to_V47 = 42,
        k_He4_Ti44_to_Cr48 = 43,
        k_p_V47_to_Cr48 = 44,
        k_He4_V47_to_Mn51 = 45,
        k_He4_Cr48_to_Fe52 = 46,
        k_p_Mn51_to_Fe52 = 47,
        k_He4_Mn51_to_Co55 = 48,
        k_He4_Fe52_to_Ni56 = 49,
        k_p_Co55_to_Ni56 = 50,
        k_n_Ni58_to_Ni59 = 51,
        k_He4_B11_to_n_N14 = 52,
        k_C12_C12_to_p_Na23 = 53,
        k_C12_C12_to_He4_Ne20 = 54,
        k_He4_N13_to_p_O16 = 55,
        k_C12_O16_to_p_Al27 = 56,
        k_C12_O16_to_He4_Mg24 = 57,
        k_O16_O16_to_n_S31 = 58,
        k_O16_O16_to_p_P31 = 59,
        k_O16_O16_to_He4_Si28 = 60,
        k_p_O17_to_He4_N14 = 61,
        k_He4_O17_to_n_Ne20 = 62,
        k_C12_Ne20_to_n_S31 = 63,
        k_C12_Ne20_to_p_P31 = 64,
        k_C12_Ne20_to_He4_Si28 = 65,
        k_He4_Ne22_to_n_Mg25 = 66,
        k_p_Na23_to_He4_Ne20 = 67,
        k_He4_Na23_to_p_Mg26 = 68,
        k_He4_Mg25_to_n_Si28 = 69,
        k_He4_Mg26_to_n_Si29 = 70,
        k_p_Al27_to_He4_Mg24 = 71,
        k_He4_Al27_to_n_P30 = 72,
        k_He4_Al27_to_p_Si30 = 73,
        k_n_P30_to_p_Si30 = 74,
        k_p_P31_to_He4_Si28 = 75,
        k_n_S31_to_p_P31 = 76,
        k_n_S31_to_He4_Si28 = 77,
        k_n_S32_to_He4_Si29 = 78,
        k_p_Cl35_to_He4_S32 = 79,
        k_He4_Cl35_to_p_Ar38 = 80,
        k_n_K39_to_p_Ar39 = 81,
        k_p_K39_to_He4_Ar36 = 82,
        k_n_Ca40_to_He4_Ar37 = 83,
        k_p_Sc43_to_He4_Ca40 = 84,
        k_He4_Ti44_to_p_V47 = 85,
        k_He4_Cr48_to_p_Mn51 = 86,
        k_He4_Fe52_to_p_Co55 = 87,
        k_He4_Co55_to_p_Ni58 = 88,
        k_n_Ni59_to_He4_Fe56 = 89,
        k_p_B11_to_He4_He4_He4 = 90,
        k_He4_He4_He4_to_C12 = 91,
        k_Ar39_to_K39 = 92,
        k_K39_to_Ar39 = 93,
        k_P30_to_Si30 = 94,
        k_P31_to_S31 = 95,
        k_S31_to_P31 = 96,
        k_Si30_to_P30 = 97,
        k_n_to_p = 98,
        k_p_to_n = 99,
        k_C12_to_p_B11_derived = 100,
        k_N13_to_p_C12_derived = 101,
        k_N14_to_n_N13_derived = 102,
        k_O16_to_He4_C12_derived = 103,
        k_O17_to_n_O16_derived = 104,
        k_Ne20_to_He4_O16_derived = 105,
        k_Na23_to_p_Ne22_derived = 106,
        k_Mg24_to_p_Na23_derived = 107,
        k_Mg24_to_He4_Ne20_derived = 108,
        k_Mg25_to_n_Mg24_derived = 109,
        k_Mg26_to_n_Mg25_derived = 110,
        k_Mg26_to_He4_Ne22_derived = 111,
        k_Al27_to_p_Mg26_derived = 112,
        k_Al27_to_He4_Na23_derived = 113,
        k_Si28_to_p_Al27_derived = 114,
        k_Si28_to_He4_Mg24_derived = 115,
        k_Si29_to_n_Si28_derived = 116,
        k_Si29_to_He4_Mg25_derived = 117,
        k_Si30_to_n_Si29_derived = 118,
        k_Si30_to_He4_Mg26_derived = 119,
        k_P30_to_p_Si29_derived = 120,
        k_P31_to_n_P30_derived = 121,
        k_P31_to_p_Si30_derived = 122,
        k_P31_to_He4_Al27_derived = 123,
        k_S31_to_p_P30_derived = 124,
        k_S32_to_n_S31_derived = 125,
        k_S32_to_p_P31_derived = 126,
        k_S32_to_He4_Si28_derived = 127,
        k_Cl35_to_He4_P31_derived = 128,
        k_Ar36_to_p_Cl35_derived = 129,
        k_Ar36_to_He4_S32_derived = 130,
        k_Ar37_to_n_Ar36_derived = 131,
        k_Ar38_to_n_Ar37_derived = 132,
        k_Ar39_to_n_Ar38_derived = 133,
        k_K39_to_p_Ar38_derived = 134,
        k_K39_to_He4_Cl35_derived = 135,
        k_Ca40_to_p_K39_derived = 136,
        k_Ca40_to_He4_Ar36_derived = 137,
        k_Sc43_to_He4_K39_derived = 138,
        k_Ti44_to_p_Sc43_derived = 139,
        k_Ti44_to_He4_Ca40_derived = 140,
        k_V47_to_He4_Sc43_derived = 141,
        k_Cr48_to_p_V47_derived = 142,
        k_Cr48_to_He4_Ti44_derived = 143,
        k_Mn51_to_He4_V47_derived = 144,
        k_Fe52_to_p_Mn51_derived = 145,
        k_Fe52_to_He4_Cr48_derived = 146,
        k_Co55_to_He4_Mn51_derived = 147,
        k_Ni56_to_p_Co55_derived = 148,
        k_Ni56_to_He4_Fe52_derived = 149,
        k_Ni59_to_n_Ni58_derived = 150,
        k_C12_to_He4_He4_He4_derived = 151,
        k_n_N14_to_He4_B11_derived = 152,
        k_He4_N14_to_p_O17_derived = 153,
        k_p_O16_to_He4_N13_derived = 154,
        k_n_Ne20_to_He4_O17_derived = 155,
        k_He4_Ne20_to_p_Na23_derived = 156,
        k_He4_Ne20_to_C12_C12_derived = 157,
        k_p_Na23_to_C12_C12_derived = 158,
        k_He4_Mg24_to_p_Al27_derived = 159,
        k_He4_Mg24_to_C12_O16_derived = 160,
        k_n_Mg25_to_He4_Ne22_derived = 161,
        k_p_Mg26_to_He4_Na23_derived = 162,
        k_p_Al27_to_C12_O16_derived = 163,
        k_n_Si28_to_He4_Mg25_derived = 164,
        k_He4_Si28_to_n_S31_derived = 165,
        k_He4_Si28_to_p_P31_derived = 166,
        k_He4_Si28_to_C12_Ne20_derived = 167,
        k_He4_Si28_to_O16_O16_derived = 168,
        k_n_Si29_to_He4_Mg26_derived = 169,
        k_He4_Si29_to_n_S32_derived = 170,
        k_p_Si30_to_n_P30_derived = 171,
        k_p_Si30_to_He4_Al27_derived = 172,
        k_n_P30_to_He4_Al27_derived = 173,
        k_p_P31_to_n_S31_derived = 174,
        k_p_P31_to_C12_Ne20_derived = 175,
        k_p_P31_to_O16_O16_derived = 176,
        k_n_S31_to_C12_Ne20_derived = 177,
        k_n_S31_to_O16_O16_derived = 178,
        k_He4_S32_to_p_Cl35_derived = 179,
        k_He4_Ar36_to_p_K39_derived = 180,
        k_He4_Ar37_to_n_Ca40_derived = 181,
        k_p_Ar38_to_He4_Cl35_derived = 182,
        k_p_Ar39_to_n_K39_derived = 183,
        k_He4_Ca40_to_p_Sc43_derived = 184,
        k_p_V47_to_He4_Ti44_derived = 185,
        k_p_Mn51_to_He4_Cr48_derived = 186,
        k_He4_Fe56_to_n_Ni59_derived = 187,
        k_p_Co55_to_He4_Fe52_derived = 188,
        k_p_Ni58_to_He4_Co55_derived = 189,
        k_He4_He4_He4_to_p_B11_derived = 190,
        NumRates = k_He4_He4_He4_to_p_B11_derived
    };

    // number of reaclib rates

    const int NrateReaclib = 182;

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
        "He4_K39_to_Sc43",  // 39,
        "He4_Ca40_to_Ti44",  // 40,
        "p_Sc43_to_Ti44",  // 41,
        "He4_Sc43_to_V47",  // 42,
        "He4_Ti44_to_Cr48",  // 43,
        "p_V47_to_Cr48",  // 44,
        "He4_V47_to_Mn51",  // 45,
        "He4_Cr48_to_Fe52",  // 46,
        "p_Mn51_to_Fe52",  // 47,
        "He4_Mn51_to_Co55",  // 48,
        "He4_Fe52_to_Ni56",  // 49,
        "p_Co55_to_Ni56",  // 50,
        "n_Ni58_to_Ni59",  // 51,
        "He4_B11_to_n_N14",  // 52,
        "C12_C12_to_p_Na23",  // 53,
        "C12_C12_to_He4_Ne20",  // 54,
        "He4_N13_to_p_O16",  // 55,
        "C12_O16_to_p_Al27",  // 56,
        "C12_O16_to_He4_Mg24",  // 57,
        "O16_O16_to_n_S31",  // 58,
        "O16_O16_to_p_P31",  // 59,
        "O16_O16_to_He4_Si28",  // 60,
        "p_O17_to_He4_N14",  // 61,
        "He4_O17_to_n_Ne20",  // 62,
        "C12_Ne20_to_n_S31",  // 63,
        "C12_Ne20_to_p_P31",  // 64,
        "C12_Ne20_to_He4_Si28",  // 65,
        "He4_Ne22_to_n_Mg25",  // 66,
        "p_Na23_to_He4_Ne20",  // 67,
        "He4_Na23_to_p_Mg26",  // 68,
        "He4_Mg25_to_n_Si28",  // 69,
        "He4_Mg26_to_n_Si29",  // 70,
        "p_Al27_to_He4_Mg24",  // 71,
        "He4_Al27_to_n_P30",  // 72,
        "He4_Al27_to_p_Si30",  // 73,
        "n_P30_to_p_Si30",  // 74,
        "p_P31_to_He4_Si28",  // 75,
        "n_S31_to_p_P31",  // 76,
        "n_S31_to_He4_Si28",  // 77,
        "n_S32_to_He4_Si29",  // 78,
        "p_Cl35_to_He4_S32",  // 79,
        "He4_Cl35_to_p_Ar38",  // 80,
        "n_K39_to_p_Ar39",  // 81,
        "p_K39_to_He4_Ar36",  // 82,
        "n_Ca40_to_He4_Ar37",  // 83,
        "p_Sc43_to_He4_Ca40",  // 84,
        "He4_Ti44_to_p_V47",  // 85,
        "He4_Cr48_to_p_Mn51",  // 86,
        "He4_Fe52_to_p_Co55",  // 87,
        "He4_Co55_to_p_Ni58",  // 88,
        "n_Ni59_to_He4_Fe56",  // 89,
        "p_B11_to_He4_He4_He4",  // 90,
        "He4_He4_He4_to_C12",  // 91,
        "Ar39_to_K39",  // 92,
        "K39_to_Ar39",  // 93,
        "P30_to_Si30",  // 94,
        "P31_to_S31",  // 95,
        "S31_to_P31",  // 96,
        "Si30_to_P30",  // 97,
        "n_to_p",  // 98,
        "p_to_n",  // 99,
        "C12_to_p_B11_derived",  // 100,
        "N13_to_p_C12_derived",  // 101,
        "N14_to_n_N13_derived",  // 102,
        "O16_to_He4_C12_derived",  // 103,
        "O17_to_n_O16_derived",  // 104,
        "Ne20_to_He4_O16_derived",  // 105,
        "Na23_to_p_Ne22_derived",  // 106,
        "Mg24_to_p_Na23_derived",  // 107,
        "Mg24_to_He4_Ne20_derived",  // 108,
        "Mg25_to_n_Mg24_derived",  // 109,
        "Mg26_to_n_Mg25_derived",  // 110,
        "Mg26_to_He4_Ne22_derived",  // 111,
        "Al27_to_p_Mg26_derived",  // 112,
        "Al27_to_He4_Na23_derived",  // 113,
        "Si28_to_p_Al27_derived",  // 114,
        "Si28_to_He4_Mg24_derived",  // 115,
        "Si29_to_n_Si28_derived",  // 116,
        "Si29_to_He4_Mg25_derived",  // 117,
        "Si30_to_n_Si29_derived",  // 118,
        "Si30_to_He4_Mg26_derived",  // 119,
        "P30_to_p_Si29_derived",  // 120,
        "P31_to_n_P30_derived",  // 121,
        "P31_to_p_Si30_derived",  // 122,
        "P31_to_He4_Al27_derived",  // 123,
        "S31_to_p_P30_derived",  // 124,
        "S32_to_n_S31_derived",  // 125,
        "S32_to_p_P31_derived",  // 126,
        "S32_to_He4_Si28_derived",  // 127,
        "Cl35_to_He4_P31_derived",  // 128,
        "Ar36_to_p_Cl35_derived",  // 129,
        "Ar36_to_He4_S32_derived",  // 130,
        "Ar37_to_n_Ar36_derived",  // 131,
        "Ar38_to_n_Ar37_derived",  // 132,
        "Ar39_to_n_Ar38_derived",  // 133,
        "K39_to_p_Ar38_derived",  // 134,
        "K39_to_He4_Cl35_derived",  // 135,
        "Ca40_to_p_K39_derived",  // 136,
        "Ca40_to_He4_Ar36_derived",  // 137,
        "Sc43_to_He4_K39_derived",  // 138,
        "Ti44_to_p_Sc43_derived",  // 139,
        "Ti44_to_He4_Ca40_derived",  // 140,
        "V47_to_He4_Sc43_derived",  // 141,
        "Cr48_to_p_V47_derived",  // 142,
        "Cr48_to_He4_Ti44_derived",  // 143,
        "Mn51_to_He4_V47_derived",  // 144,
        "Fe52_to_p_Mn51_derived",  // 145,
        "Fe52_to_He4_Cr48_derived",  // 146,
        "Co55_to_He4_Mn51_derived",  // 147,
        "Ni56_to_p_Co55_derived",  // 148,
        "Ni56_to_He4_Fe52_derived",  // 149,
        "Ni59_to_n_Ni58_derived",  // 150,
        "C12_to_He4_He4_He4_derived",  // 151,
        "n_N14_to_He4_B11_derived",  // 152,
        "He4_N14_to_p_O17_derived",  // 153,
        "p_O16_to_He4_N13_derived",  // 154,
        "n_Ne20_to_He4_O17_derived",  // 155,
        "He4_Ne20_to_p_Na23_derived",  // 156,
        "He4_Ne20_to_C12_C12_derived",  // 157,
        "p_Na23_to_C12_C12_derived",  // 158,
        "He4_Mg24_to_p_Al27_derived",  // 159,
        "He4_Mg24_to_C12_O16_derived",  // 160,
        "n_Mg25_to_He4_Ne22_derived",  // 161,
        "p_Mg26_to_He4_Na23_derived",  // 162,
        "p_Al27_to_C12_O16_derived",  // 163,
        "n_Si28_to_He4_Mg25_derived",  // 164,
        "He4_Si28_to_n_S31_derived",  // 165,
        "He4_Si28_to_p_P31_derived",  // 166,
        "He4_Si28_to_C12_Ne20_derived",  // 167,
        "He4_Si28_to_O16_O16_derived",  // 168,
        "n_Si29_to_He4_Mg26_derived",  // 169,
        "He4_Si29_to_n_S32_derived",  // 170,
        "p_Si30_to_n_P30_derived",  // 171,
        "p_Si30_to_He4_Al27_derived",  // 172,
        "n_P30_to_He4_Al27_derived",  // 173,
        "p_P31_to_n_S31_derived",  // 174,
        "p_P31_to_C12_Ne20_derived",  // 175,
        "p_P31_to_O16_O16_derived",  // 176,
        "n_S31_to_C12_Ne20_derived",  // 177,
        "n_S31_to_O16_O16_derived",  // 178,
        "He4_S32_to_p_Cl35_derived",  // 179,
        "He4_Ar36_to_p_K39_derived",  // 180,
        "He4_Ar37_to_n_Ca40_derived",  // 181,
        "p_Ar38_to_He4_Cl35_derived",  // 182,
        "p_Ar39_to_n_K39_derived",  // 183,
        "He4_Ca40_to_p_Sc43_derived",  // 184,
        "p_V47_to_He4_Ti44_derived",  // 185,
        "p_Mn51_to_He4_Cr48_derived",  // 186,
        "He4_Fe56_to_n_Ni59_derived",  // 187,
        "p_Co55_to_He4_Fe52_derived",  // 188,
        "p_Ni58_to_He4_Co55_derived",  // 189,
        "He4_He4_He4_to_p_B11_derived"  // 190,
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
