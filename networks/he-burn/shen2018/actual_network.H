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
        k_C12_to_p_B11 = 1,
        k_N13_to_p_C12 = 2,
        k_N14_to_n_N13 = 3,
        k_O16_to_He4_C12 = 4,
        k_O17_to_n_O16 = 5,
        k_Ne20_to_He4_O16 = 6,
        k_Na23_to_p_Ne22 = 7,
        k_Mg24_to_p_Na23 = 8,
        k_Mg24_to_He4_Ne20 = 9,
        k_Mg25_to_n_Mg24 = 10,
        k_Mg26_to_n_Mg25 = 11,
        k_Mg26_to_He4_Ne22 = 12,
        k_Al27_to_p_Mg26 = 13,
        k_Al27_to_He4_Na23 = 14,
        k_Si28_to_p_Al27 = 15,
        k_Si28_to_He4_Mg24 = 16,
        k_Si29_to_n_Si28 = 17,
        k_Si29_to_He4_Mg25 = 18,
        k_Si30_to_n_Si29 = 19,
        k_Si30_to_He4_Mg26 = 20,
        k_P30_to_p_Si29 = 21,
        k_P31_to_n_P30 = 22,
        k_P31_to_p_Si30 = 23,
        k_P31_to_He4_Al27 = 24,
        k_S31_to_p_P30 = 25,
        k_S32_to_n_S31 = 26,
        k_S32_to_p_P31 = 27,
        k_S32_to_He4_Si28 = 28,
        k_Cl35_to_He4_P31 = 29,
        k_Ar36_to_p_Cl35 = 30,
        k_Ar36_to_He4_S32 = 31,
        k_Ar37_to_n_Ar36 = 32,
        k_Ar38_to_n_Ar37 = 33,
        k_Ar39_to_n_Ar38 = 34,
        k_K39_to_p_Ar38 = 35,
        k_K39_to_He4_Cl35 = 36,
        k_Ca40_to_p_K39 = 37,
        k_Ca40_to_He4_Ar36 = 38,
        k_Sc43_to_He4_K39 = 39,
        k_Ti44_to_p_Sc43 = 40,
        k_Ti44_to_He4_Ca40 = 41,
        k_V47_to_He4_Sc43 = 42,
        k_Cr48_to_p_V47 = 43,
        k_Cr48_to_He4_Ti44 = 44,
        k_Mn51_to_He4_V47 = 45,
        k_Fe52_to_p_Mn51 = 46,
        k_Fe52_to_He4_Cr48 = 47,
        k_Co55_to_He4_Mn51 = 48,
        k_Ni56_to_p_Co55 = 49,
        k_Ni56_to_He4_Fe52 = 50,
        k_Ni59_to_n_Ni58 = 51,
        k_C12_to_He4_He4_He4 = 52,
        k_p_B11_to_C12 = 53,
        k_p_C12_to_N13 = 54,
        k_He4_C12_to_O16 = 55,
        k_n_N13_to_N14 = 56,
        k_n_O16_to_O17 = 57,
        k_He4_O16_to_Ne20 = 58,
        k_He4_Ne20_to_Mg24 = 59,
        k_p_Ne22_to_Na23 = 60,
        k_He4_Ne22_to_Mg26 = 61,
        k_p_Na23_to_Mg24 = 62,
        k_He4_Na23_to_Al27 = 63,
        k_n_Mg24_to_Mg25 = 64,
        k_He4_Mg24_to_Si28 = 65,
        k_n_Mg25_to_Mg26 = 66,
        k_He4_Mg25_to_Si29 = 67,
        k_p_Mg26_to_Al27 = 68,
        k_He4_Mg26_to_Si30 = 69,
        k_p_Al27_to_Si28 = 70,
        k_He4_Al27_to_P31 = 71,
        k_n_Si28_to_Si29 = 72,
        k_He4_Si28_to_S32 = 73,
        k_n_Si29_to_Si30 = 74,
        k_p_Si29_to_P30 = 75,
        k_p_Si30_to_P31 = 76,
        k_n_P30_to_P31 = 77,
        k_p_P30_to_S31 = 78,
        k_p_P31_to_S32 = 79,
        k_He4_P31_to_Cl35 = 80,
        k_n_S31_to_S32 = 81,
        k_He4_S32_to_Ar36 = 82,
        k_p_Cl35_to_Ar36 = 83,
        k_He4_Cl35_to_K39 = 84,
        k_n_Ar36_to_Ar37 = 85,
        k_He4_Ar36_to_Ca40 = 86,
        k_n_Ar37_to_Ar38 = 87,
        k_n_Ar38_to_Ar39 = 88,
        k_p_Ar38_to_K39 = 89,
        k_p_K39_to_Ca40 = 90,
        k_He4_K39_to_Sc43 = 91,
        k_He4_Ca40_to_Ti44 = 92,
        k_p_Sc43_to_Ti44 = 93,
        k_He4_Sc43_to_V47 = 94,
        k_He4_Ti44_to_Cr48 = 95,
        k_p_V47_to_Cr48 = 96,
        k_He4_V47_to_Mn51 = 97,
        k_He4_Cr48_to_Fe52 = 98,
        k_p_Mn51_to_Fe52 = 99,
        k_He4_Mn51_to_Co55 = 100,
        k_He4_Fe52_to_Ni56 = 101,
        k_p_Co55_to_Ni56 = 102,
        k_n_Ni58_to_Ni59 = 103,
        k_He4_B11_to_n_N14 = 104,
        k_C12_C12_to_p_Na23 = 105,
        k_C12_C12_to_He4_Ne20 = 106,
        k_He4_N13_to_p_O16 = 107,
        k_n_N14_to_He4_B11 = 108,
        k_He4_N14_to_p_O17 = 109,
        k_p_O16_to_He4_N13 = 110,
        k_C12_O16_to_p_Al27 = 111,
        k_C12_O16_to_He4_Mg24 = 112,
        k_O16_O16_to_n_S31 = 113,
        k_O16_O16_to_p_P31 = 114,
        k_O16_O16_to_He4_Si28 = 115,
        k_p_O17_to_He4_N14 = 116,
        k_He4_O17_to_n_Ne20 = 117,
        k_n_Ne20_to_He4_O17 = 118,
        k_He4_Ne20_to_p_Na23 = 119,
        k_He4_Ne20_to_C12_C12 = 120,
        k_C12_Ne20_to_n_S31 = 121,
        k_C12_Ne20_to_p_P31 = 122,
        k_C12_Ne20_to_He4_Si28 = 123,
        k_He4_Ne22_to_n_Mg25 = 124,
        k_p_Na23_to_He4_Ne20 = 125,
        k_p_Na23_to_C12_C12 = 126,
        k_He4_Na23_to_p_Mg26 = 127,
        k_He4_Mg24_to_p_Al27 = 128,
        k_He4_Mg24_to_C12_O16 = 129,
        k_n_Mg25_to_He4_Ne22 = 130,
        k_He4_Mg25_to_n_Si28 = 131,
        k_p_Mg26_to_He4_Na23 = 132,
        k_He4_Mg26_to_n_Si29 = 133,
        k_p_Al27_to_He4_Mg24 = 134,
        k_p_Al27_to_C12_O16 = 135,
        k_He4_Al27_to_n_P30 = 136,
        k_He4_Al27_to_p_Si30 = 137,
        k_n_Si28_to_He4_Mg25 = 138,
        k_He4_Si28_to_n_S31 = 139,
        k_He4_Si28_to_p_P31 = 140,
        k_He4_Si28_to_C12_Ne20 = 141,
        k_He4_Si28_to_O16_O16 = 142,
        k_n_Si29_to_He4_Mg26 = 143,
        k_He4_Si29_to_n_S32 = 144,
        k_p_Si30_to_n_P30 = 145,
        k_p_Si30_to_He4_Al27 = 146,
        k_n_P30_to_p_Si30 = 147,
        k_n_P30_to_He4_Al27 = 148,
        k_p_P31_to_n_S31 = 149,
        k_p_P31_to_He4_Si28 = 150,
        k_p_P31_to_C12_Ne20 = 151,
        k_p_P31_to_O16_O16 = 152,
        k_n_S31_to_p_P31 = 153,
        k_n_S31_to_He4_Si28 = 154,
        k_n_S31_to_C12_Ne20 = 155,
        k_n_S31_to_O16_O16 = 156,
        k_n_S32_to_He4_Si29 = 157,
        k_He4_S32_to_p_Cl35 = 158,
        k_p_Cl35_to_He4_S32 = 159,
        k_He4_Cl35_to_p_Ar38 = 160,
        k_He4_Ar36_to_p_K39 = 161,
        k_He4_Ar37_to_n_Ca40 = 162,
        k_p_Ar38_to_He4_Cl35 = 163,
        k_p_Ar39_to_n_K39 = 164,
        k_n_K39_to_p_Ar39 = 165,
        k_p_K39_to_He4_Ar36 = 166,
        k_n_Ca40_to_He4_Ar37 = 167,
        k_He4_Ca40_to_p_Sc43 = 168,
        k_p_Sc43_to_He4_Ca40 = 169,
        k_He4_Ti44_to_p_V47 = 170,
        k_p_V47_to_He4_Ti44 = 171,
        k_He4_Cr48_to_p_Mn51 = 172,
        k_p_Mn51_to_He4_Cr48 = 173,
        k_He4_Fe52_to_p_Co55 = 174,
        k_He4_Fe56_to_n_Ni59 = 175,
        k_p_Co55_to_He4_Fe52 = 176,
        k_He4_Co55_to_p_Ni58 = 177,
        k_p_Ni58_to_He4_Co55 = 178,
        k_n_Ni59_to_He4_Fe56 = 179,
        k_p_B11_to_He4_He4_He4 = 180,
        k_He4_He4_He4_to_C12 = 181,
        k_He4_He4_He4_to_p_B11 = 182,
        k_Ar39_to_K39 = 183,
        k_K39_to_Ar39 = 184,
        k_P30_to_Si30 = 185,
        k_P31_to_S31 = 186,
        k_S31_to_P31 = 187,
        k_Si30_to_P30 = 188,
        k_n_to_p = 189,
        k_p_to_n = 190,
        NumRates = k_p_to_n
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
        "C12_to_p_B11",  // 1,
        "N13_to_p_C12",  // 2,
        "N14_to_n_N13",  // 3,
        "O16_to_He4_C12",  // 4,
        "O17_to_n_O16",  // 5,
        "Ne20_to_He4_O16",  // 6,
        "Na23_to_p_Ne22",  // 7,
        "Mg24_to_p_Na23",  // 8,
        "Mg24_to_He4_Ne20",  // 9,
        "Mg25_to_n_Mg24",  // 10,
        "Mg26_to_n_Mg25",  // 11,
        "Mg26_to_He4_Ne22",  // 12,
        "Al27_to_p_Mg26",  // 13,
        "Al27_to_He4_Na23",  // 14,
        "Si28_to_p_Al27",  // 15,
        "Si28_to_He4_Mg24",  // 16,
        "Si29_to_n_Si28",  // 17,
        "Si29_to_He4_Mg25",  // 18,
        "Si30_to_n_Si29",  // 19,
        "Si30_to_He4_Mg26",  // 20,
        "P30_to_p_Si29",  // 21,
        "P31_to_n_P30",  // 22,
        "P31_to_p_Si30",  // 23,
        "P31_to_He4_Al27",  // 24,
        "S31_to_p_P30",  // 25,
        "S32_to_n_S31",  // 26,
        "S32_to_p_P31",  // 27,
        "S32_to_He4_Si28",  // 28,
        "Cl35_to_He4_P31",  // 29,
        "Ar36_to_p_Cl35",  // 30,
        "Ar36_to_He4_S32",  // 31,
        "Ar37_to_n_Ar36",  // 32,
        "Ar38_to_n_Ar37",  // 33,
        "Ar39_to_n_Ar38",  // 34,
        "K39_to_p_Ar38",  // 35,
        "K39_to_He4_Cl35",  // 36,
        "Ca40_to_p_K39",  // 37,
        "Ca40_to_He4_Ar36",  // 38,
        "Sc43_to_He4_K39",  // 39,
        "Ti44_to_p_Sc43",  // 40,
        "Ti44_to_He4_Ca40",  // 41,
        "V47_to_He4_Sc43",  // 42,
        "Cr48_to_p_V47",  // 43,
        "Cr48_to_He4_Ti44",  // 44,
        "Mn51_to_He4_V47",  // 45,
        "Fe52_to_p_Mn51",  // 46,
        "Fe52_to_He4_Cr48",  // 47,
        "Co55_to_He4_Mn51",  // 48,
        "Ni56_to_p_Co55",  // 49,
        "Ni56_to_He4_Fe52",  // 50,
        "Ni59_to_n_Ni58",  // 51,
        "C12_to_He4_He4_He4",  // 52,
        "p_B11_to_C12",  // 53,
        "p_C12_to_N13",  // 54,
        "He4_C12_to_O16",  // 55,
        "n_N13_to_N14",  // 56,
        "n_O16_to_O17",  // 57,
        "He4_O16_to_Ne20",  // 58,
        "He4_Ne20_to_Mg24",  // 59,
        "p_Ne22_to_Na23",  // 60,
        "He4_Ne22_to_Mg26",  // 61,
        "p_Na23_to_Mg24",  // 62,
        "He4_Na23_to_Al27",  // 63,
        "n_Mg24_to_Mg25",  // 64,
        "He4_Mg24_to_Si28",  // 65,
        "n_Mg25_to_Mg26",  // 66,
        "He4_Mg25_to_Si29",  // 67,
        "p_Mg26_to_Al27",  // 68,
        "He4_Mg26_to_Si30",  // 69,
        "p_Al27_to_Si28",  // 70,
        "He4_Al27_to_P31",  // 71,
        "n_Si28_to_Si29",  // 72,
        "He4_Si28_to_S32",  // 73,
        "n_Si29_to_Si30",  // 74,
        "p_Si29_to_P30",  // 75,
        "p_Si30_to_P31",  // 76,
        "n_P30_to_P31",  // 77,
        "p_P30_to_S31",  // 78,
        "p_P31_to_S32",  // 79,
        "He4_P31_to_Cl35",  // 80,
        "n_S31_to_S32",  // 81,
        "He4_S32_to_Ar36",  // 82,
        "p_Cl35_to_Ar36",  // 83,
        "He4_Cl35_to_K39",  // 84,
        "n_Ar36_to_Ar37",  // 85,
        "He4_Ar36_to_Ca40",  // 86,
        "n_Ar37_to_Ar38",  // 87,
        "n_Ar38_to_Ar39",  // 88,
        "p_Ar38_to_K39",  // 89,
        "p_K39_to_Ca40",  // 90,
        "He4_K39_to_Sc43",  // 91,
        "He4_Ca40_to_Ti44",  // 92,
        "p_Sc43_to_Ti44",  // 93,
        "He4_Sc43_to_V47",  // 94,
        "He4_Ti44_to_Cr48",  // 95,
        "p_V47_to_Cr48",  // 96,
        "He4_V47_to_Mn51",  // 97,
        "He4_Cr48_to_Fe52",  // 98,
        "p_Mn51_to_Fe52",  // 99,
        "He4_Mn51_to_Co55",  // 100,
        "He4_Fe52_to_Ni56",  // 101,
        "p_Co55_to_Ni56",  // 102,
        "n_Ni58_to_Ni59",  // 103,
        "He4_B11_to_n_N14",  // 104,
        "C12_C12_to_p_Na23",  // 105,
        "C12_C12_to_He4_Ne20",  // 106,
        "He4_N13_to_p_O16",  // 107,
        "n_N14_to_He4_B11",  // 108,
        "He4_N14_to_p_O17",  // 109,
        "p_O16_to_He4_N13",  // 110,
        "C12_O16_to_p_Al27",  // 111,
        "C12_O16_to_He4_Mg24",  // 112,
        "O16_O16_to_n_S31",  // 113,
        "O16_O16_to_p_P31",  // 114,
        "O16_O16_to_He4_Si28",  // 115,
        "p_O17_to_He4_N14",  // 116,
        "He4_O17_to_n_Ne20",  // 117,
        "n_Ne20_to_He4_O17",  // 118,
        "He4_Ne20_to_p_Na23",  // 119,
        "He4_Ne20_to_C12_C12",  // 120,
        "C12_Ne20_to_n_S31",  // 121,
        "C12_Ne20_to_p_P31",  // 122,
        "C12_Ne20_to_He4_Si28",  // 123,
        "He4_Ne22_to_n_Mg25",  // 124,
        "p_Na23_to_He4_Ne20",  // 125,
        "p_Na23_to_C12_C12",  // 126,
        "He4_Na23_to_p_Mg26",  // 127,
        "He4_Mg24_to_p_Al27",  // 128,
        "He4_Mg24_to_C12_O16",  // 129,
        "n_Mg25_to_He4_Ne22",  // 130,
        "He4_Mg25_to_n_Si28",  // 131,
        "p_Mg26_to_He4_Na23",  // 132,
        "He4_Mg26_to_n_Si29",  // 133,
        "p_Al27_to_He4_Mg24",  // 134,
        "p_Al27_to_C12_O16",  // 135,
        "He4_Al27_to_n_P30",  // 136,
        "He4_Al27_to_p_Si30",  // 137,
        "n_Si28_to_He4_Mg25",  // 138,
        "He4_Si28_to_n_S31",  // 139,
        "He4_Si28_to_p_P31",  // 140,
        "He4_Si28_to_C12_Ne20",  // 141,
        "He4_Si28_to_O16_O16",  // 142,
        "n_Si29_to_He4_Mg26",  // 143,
        "He4_Si29_to_n_S32",  // 144,
        "p_Si30_to_n_P30",  // 145,
        "p_Si30_to_He4_Al27",  // 146,
        "n_P30_to_p_Si30",  // 147,
        "n_P30_to_He4_Al27",  // 148,
        "p_P31_to_n_S31",  // 149,
        "p_P31_to_He4_Si28",  // 150,
        "p_P31_to_C12_Ne20",  // 151,
        "p_P31_to_O16_O16",  // 152,
        "n_S31_to_p_P31",  // 153,
        "n_S31_to_He4_Si28",  // 154,
        "n_S31_to_C12_Ne20",  // 155,
        "n_S31_to_O16_O16",  // 156,
        "n_S32_to_He4_Si29",  // 157,
        "He4_S32_to_p_Cl35",  // 158,
        "p_Cl35_to_He4_S32",  // 159,
        "He4_Cl35_to_p_Ar38",  // 160,
        "He4_Ar36_to_p_K39",  // 161,
        "He4_Ar37_to_n_Ca40",  // 162,
        "p_Ar38_to_He4_Cl35",  // 163,
        "p_Ar39_to_n_K39",  // 164,
        "n_K39_to_p_Ar39",  // 165,
        "p_K39_to_He4_Ar36",  // 166,
        "n_Ca40_to_He4_Ar37",  // 167,
        "He4_Ca40_to_p_Sc43",  // 168,
        "p_Sc43_to_He4_Ca40",  // 169,
        "He4_Ti44_to_p_V47",  // 170,
        "p_V47_to_He4_Ti44",  // 171,
        "He4_Cr48_to_p_Mn51",  // 172,
        "p_Mn51_to_He4_Cr48",  // 173,
        "He4_Fe52_to_p_Co55",  // 174,
        "He4_Fe56_to_n_Ni59",  // 175,
        "p_Co55_to_He4_Fe52",  // 176,
        "He4_Co55_to_p_Ni58",  // 177,
        "p_Ni58_to_He4_Co55",  // 178,
        "n_Ni59_to_He4_Fe56",  // 179,
        "p_B11_to_He4_He4_He4",  // 180,
        "He4_He4_He4_to_C12",  // 181,
        "He4_He4_He4_to_p_B11",  // 182,
        "Ar39_to_K39",  // 183,
        "K39_to_Ar39",  // 184,
        "P30_to_Si30",  // 185,
        "P31_to_S31",  // 186,
        "S31_to_P31",  // 187,
        "Si30_to_P30",  // 188,
        "n_to_p",  // 189,
        "p_to_n"  // 190,
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
