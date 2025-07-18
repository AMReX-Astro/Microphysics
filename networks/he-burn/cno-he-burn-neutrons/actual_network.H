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

        if constexpr (spec == N) {
            return 0.0_rt;
        }
        else if constexpr (spec == H1) {
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
        else if constexpr (spec == Ne22) {
            return 177.77024384000106_rt;
        }
        else if constexpr (spec == Na21) {
            return 163.07672270400144_rt;
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
        using namespace Species;

        static_assert(spec >= 1 && spec <= NumSpec);

        // Set the mass of the element

        if constexpr (spec == N) {
            return 1.674927498034172e-24_rt;
        }
        else if constexpr (spec == H1) {
            return 1.6735328377636005e-24_rt;
        }
        else if constexpr (spec == He4) {
            return 6.646479071584587e-24_rt;
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
        else if constexpr (spec == O18) {
            return 2.9888307701990364e-23_rt;
        }
        else if constexpr (spec == F17) {
            return 2.8232643353472436e-23_rt;
        }
        else if constexpr (spec == F18) {
            return 2.989125964092377e-23_rt;
        }
        else if constexpr (spec == F19) {
            return 3.1547590653659727e-23_rt;
        }
        else if constexpr (spec == Ne18) {
            return 2.9899182681834417e-23_rt;
        }
        else if constexpr (spec == Ne19) {
            return 3.155336557821978e-23_rt;
        }
        else if constexpr (spec == Ne20) {
            return 3.3198227947612416e-23_rt;
        }
        else if constexpr (spec == Ne21) {
            return 3.4861102572650884e-23_rt;
        }
        else if constexpr (spec == Ne22) {
            return 3.651755410955497e-23_rt;
        }
        else if constexpr (spec == Na21) {
            return 3.486742553187392e-23_rt;
        }
        else if constexpr (spec == Na22) {
            return 3.652262279854593e-23_rt;
        }
        else if constexpr (spec == Na23) {
            return 3.817541002484691e-23_rt;
        }
        else if constexpr (spec == Mg24) {
            return 3.9828098739467446e-23_rt;
        }
        else if constexpr (spec == Al27) {
            return 4.480389861070653e-23_rt;
        }
        else if constexpr (spec == Si28) {
            return 4.6456779473820677e-23_rt;
        }
        else if constexpr (spec == P31) {
            return 5.14331418367544e-23_rt;
        }
        else if constexpr (spec == S32) {
            return 5.309087322384128e-23_rt;
        }
        else if constexpr (spec == Ar36) {
            return 5.972551377884467e-23_rt;
        }
        else if constexpr (spec == Ca40) {
            return 6.635944331004904e-23_rt;
        }
        else if constexpr (spec == Ti44) {
            return 7.299678247096977e-23_rt;
        }
        else if constexpr (spec == Cr48) {
            return 7.962953983065421e-23_rt;
        }
        else if constexpr (spec == Fe52) {
            return 8.626187166893794e-23_rt;
        }
        else if constexpr (spec == Ni56) {
            return 9.289408870379396e-23_rt;
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
        k_N13_to_C13_weak_wc12 = 1,
        k_O14_to_N14_weak_wc12 = 2,
        k_O15_to_N15_weak_wc12 = 3,
        k_n_C12_to_C13 = 4,
        k_p_C12_to_N13 = 5,
        k_He4_C12_to_O16 = 6,
        k_p_C13_to_N14 = 7,
        k_n_N13_to_N14 = 8,
        k_p_N13_to_O14 = 9,
        k_n_N14_to_N15 = 10,
        k_p_N14_to_O15 = 11,
        k_He4_N14_to_F18 = 12,
        k_p_N15_to_O16 = 13,
        k_He4_N15_to_F19 = 14,
        k_n_O14_to_O15 = 15,
        k_He4_O14_to_Ne18 = 16,
        k_n_O15_to_O16 = 17,
        k_He4_O15_to_Ne19 = 18,
        k_n_O16_to_O17 = 19,
        k_p_O16_to_F17 = 20,
        k_He4_O16_to_Ne20 = 21,
        k_n_O17_to_O18 = 22,
        k_p_O17_to_F18 = 23,
        k_He4_O17_to_Ne21 = 24,
        k_p_O18_to_F19 = 25,
        k_He4_O18_to_Ne22 = 26,
        k_n_F17_to_F18 = 27,
        k_p_F17_to_Ne18 = 28,
        k_He4_F17_to_Na21 = 29,
        k_n_F18_to_F19 = 30,
        k_p_F18_to_Ne19 = 31,
        k_He4_F18_to_Na22 = 32,
        k_p_F19_to_Ne20 = 33,
        k_He4_F19_to_Na23 = 34,
        k_n_Ne18_to_Ne19 = 35,
        k_n_Ne19_to_Ne20 = 36,
        k_n_Ne20_to_Ne21 = 37,
        k_p_Ne20_to_Na21 = 38,
        k_He4_Ne20_to_Mg24 = 39,
        k_n_Ne21_to_Ne22 = 40,
        k_p_Ne21_to_Na22 = 41,
        k_p_Ne22_to_Na23 = 42,
        k_n_Na21_to_Na22 = 43,
        k_n_Na22_to_Na23 = 44,
        k_p_Na23_to_Mg24 = 45,
        k_He4_Mg24_to_Si28 = 46,
        k_p_Al27_to_Si28 = 47,
        k_He4_Si28_to_S32 = 48,
        k_p_P31_to_S32 = 49,
        k_C12_C12_to_p_Na23 = 50,
        k_C12_C12_to_He4_Ne20 = 51,
        k_p_C13_to_n_N13 = 52,
        k_He4_C13_to_n_O16 = 53,
        k_He4_N13_to_p_O16 = 54,
        k_p_N14_to_n_O14 = 55,
        k_He4_N14_to_n_F17 = 56,
        k_p_N15_to_n_O15 = 57,
        k_p_N15_to_He4_C12 = 58,
        k_He4_O14_to_p_F17 = 59,
        k_n_O15_to_He4_C12 = 60,
        k_C12_O16_to_p_Al27 = 61,
        k_C12_O16_to_He4_Mg24 = 62,
        k_O16_O16_to_p_P31 = 63,
        k_O16_O16_to_He4_Si28 = 64,
        k_p_O17_to_He4_N14 = 65,
        k_He4_O17_to_n_Ne20 = 66,
        k_p_O18_to_He4_N15 = 67,
        k_He4_O18_to_n_Ne21 = 68,
        k_n_F17_to_p_O17 = 69,
        k_n_F18_to_p_O18 = 70,
        k_n_F18_to_He4_N15 = 71,
        k_p_F18_to_He4_O15 = 72,
        k_He4_F18_to_p_Ne21 = 73,
        k_p_F19_to_n_Ne19 = 74,
        k_p_F19_to_He4_O16 = 75,
        k_He4_F19_to_p_Ne22 = 76,
        k_n_Ne18_to_p_F18 = 77,
        k_n_Ne18_to_He4_O15 = 78,
        k_He4_Ne18_to_p_Na21 = 79,
        k_n_Ne19_to_He4_O16 = 80,
        k_He4_Ne19_to_p_Na22 = 81,
        k_p_Ne20_to_He4_F17 = 82,
        k_He4_Ne21_to_n_Mg24 = 83,
        k_n_Na21_to_p_Ne21 = 84,
        k_n_Na21_to_He4_F18 = 85,
        k_n_Na22_to_p_Ne22 = 86,
        k_n_Na22_to_He4_F19 = 87,
        k_p_Na23_to_He4_Ne20 = 88,
        k_p_Mg24_to_He4_Na21 = 89,
        k_p_Al27_to_He4_Mg24 = 90,
        k_p_P31_to_He4_Si28 = 91,
        k_He4_He4_He4_to_C12 = 92,
        k_C12_C12_to_Mg24_modified = 93,
        k_O16_O16_to_S32_modified = 94,
        k_C12_O16_to_Si28_modified = 95,
        k_He4_S32_to_Ar36_removed = 96,
        k_p_Cl35_to_Ar36_removed = 97,
        k_p_Cl35_to_He4_S32_removed = 98,
        k_He4_Ar36_to_Ca40_removed = 99,
        k_p_K39_to_Ca40_removed = 100,
        k_p_K39_to_He4_Ar36_removed = 101,
        k_He4_Ca40_to_Ti44_removed = 102,
        k_p_Sc43_to_Ti44_removed = 103,
        k_p_Sc43_to_He4_Ca40_removed = 104,
        k_He4_Ti44_to_Cr48_removed = 105,
        k_He4_Ti44_to_p_V47_removed = 106,
        k_p_V47_to_Cr48_removed = 107,
        k_He4_Cr48_to_Fe52_removed = 108,
        k_He4_Cr48_to_p_Mn51_removed = 109,
        k_p_Mn51_to_Fe52_removed = 110,
        k_He4_Fe52_to_Ni56_removed = 111,
        k_He4_Fe52_to_p_Co55_removed = 112,
        k_p_Co55_to_Ni56_removed = 113,
        k_n_to_p = 114,
        k_p_to_n = 115,
        k_F17_to_O17 = 116,
        k_O17_to_F17 = 117,
        k_F18_to_Ne18 = 118,
        k_F18_to_O18 = 119,
        k_Ne18_to_F18 = 120,
        k_O18_to_F18 = 121,
        k_F19_to_Ne19 = 122,
        k_Ne19_to_F19 = 123,
        k_Na21_to_Ne21 = 124,
        k_Ne21_to_Na21 = 125,
        k_Na22_to_Ne22 = 126,
        k_Ne22_to_Na22 = 127,
        k_S32_He4_to_Ar36_approx = 128,
        k_Ar36_to_S32_He4_approx = 129,
        k_Ar36_He4_to_Ca40_approx = 130,
        k_Ca40_to_Ar36_He4_approx = 131,
        k_Ca40_He4_to_Ti44_approx = 132,
        k_Ti44_to_Ca40_He4_approx = 133,
        k_Ti44_He4_to_Cr48_approx = 134,
        k_Cr48_to_Ti44_He4_approx = 135,
        k_Cr48_He4_to_Fe52_approx = 136,
        k_Fe52_to_Cr48_He4_approx = 137,
        k_Fe52_He4_to_Ni56_approx = 138,
        k_Ni56_to_Fe52_He4_approx = 139,
        k_C13_to_n_C12_derived = 140,
        k_N13_to_p_C12_derived = 141,
        k_N14_to_n_N13_derived = 142,
        k_N14_to_p_C13_derived = 143,
        k_N15_to_n_N14_derived = 144,
        k_O14_to_p_N13_derived = 145,
        k_O15_to_n_O14_derived = 146,
        k_O15_to_p_N14_derived = 147,
        k_O16_to_n_O15_derived = 148,
        k_O16_to_p_N15_derived = 149,
        k_O16_to_He4_C12_derived = 150,
        k_O17_to_n_O16_derived = 151,
        k_O18_to_n_O17_derived = 152,
        k_F17_to_p_O16_derived = 153,
        k_F18_to_n_F17_derived = 154,
        k_F18_to_p_O17_derived = 155,
        k_F18_to_He4_N14_derived = 156,
        k_F19_to_n_F18_derived = 157,
        k_F19_to_p_O18_derived = 158,
        k_F19_to_He4_N15_derived = 159,
        k_Ne18_to_p_F17_derived = 160,
        k_Ne18_to_He4_O14_derived = 161,
        k_Ne19_to_n_Ne18_derived = 162,
        k_Ne19_to_p_F18_derived = 163,
        k_Ne19_to_He4_O15_derived = 164,
        k_Ne20_to_n_Ne19_derived = 165,
        k_Ne20_to_p_F19_derived = 166,
        k_Ne20_to_He4_O16_derived = 167,
        k_Ne21_to_n_Ne20_derived = 168,
        k_Ne21_to_He4_O17_derived = 169,
        k_Ne22_to_n_Ne21_derived = 170,
        k_Ne22_to_He4_O18_derived = 171,
        k_Na21_to_p_Ne20_derived = 172,
        k_Na21_to_He4_F17_derived = 173,
        k_Na22_to_n_Na21_derived = 174,
        k_Na22_to_p_Ne21_derived = 175,
        k_Na22_to_He4_F18_derived = 176,
        k_Na23_to_n_Na22_derived = 177,
        k_Na23_to_p_Ne22_derived = 178,
        k_Na23_to_He4_F19_derived = 179,
        k_Mg24_to_p_Na23_derived = 180,
        k_Mg24_to_He4_Ne20_derived = 181,
        k_Si28_to_p_Al27_derived = 182,
        k_Si28_to_He4_Mg24_derived = 183,
        k_S32_to_p_P31_derived = 184,
        k_S32_to_He4_Si28_derived = 185,
        k_C12_to_He4_He4_He4_derived = 186,
        k_He4_C12_to_n_O15_derived = 187,
        k_He4_C12_to_p_N15_derived = 188,
        k_n_N13_to_p_C13_derived = 189,
        k_He4_N14_to_p_O17_derived = 190,
        k_He4_N15_to_n_F18_derived = 191,
        k_He4_N15_to_p_O18_derived = 192,
        k_n_O14_to_p_N14_derived = 193,
        k_n_O15_to_p_N15_derived = 194,
        k_He4_O15_to_n_Ne18_derived = 195,
        k_He4_O15_to_p_F18_derived = 196,
        k_n_O16_to_He4_C13_derived = 197,
        k_p_O16_to_He4_N13_derived = 198,
        k_He4_O16_to_n_Ne19_derived = 199,
        k_He4_O16_to_p_F19_derived = 200,
        k_p_O17_to_n_F17_derived = 201,
        k_p_O18_to_n_F18_derived = 202,
        k_n_F17_to_He4_N14_derived = 203,
        k_p_F17_to_He4_O14_derived = 204,
        k_He4_F17_to_p_Ne20_derived = 205,
        k_p_F18_to_n_Ne18_derived = 206,
        k_He4_F18_to_n_Na21_derived = 207,
        k_He4_F19_to_n_Na22_derived = 208,
        k_n_Ne19_to_p_F19_derived = 209,
        k_n_Ne20_to_He4_O17_derived = 210,
        k_He4_Ne20_to_p_Na23_derived = 211,
        k_n_Ne21_to_He4_O18_derived = 212,
        k_p_Ne21_to_n_Na21_derived = 213,
        k_p_Ne21_to_He4_F18_derived = 214,
        k_p_Ne22_to_n_Na22_derived = 215,
        k_p_Ne22_to_He4_F19_derived = 216,
        k_p_Na21_to_He4_Ne18_derived = 217,
        k_He4_Na21_to_p_Mg24_derived = 218,
        k_p_Na22_to_He4_Ne19_derived = 219,
        k_n_Mg24_to_He4_Ne21_derived = 220,
        k_He4_Mg24_to_p_Al27_derived = 221,
        k_He4_Si28_to_p_P31_derived = 222,
        k_He4_S32_to_p_Cl35_derived_removed = 223,
        k_Ar36_to_He4_S32_derived_removed = 224,
        k_Ar36_to_p_Cl35_derived_removed = 225,
        k_He4_Ar36_to_p_K39_derived_removed = 226,
        k_Ca40_to_He4_Ar36_derived_removed = 227,
        k_Ca40_to_p_K39_derived_removed = 228,
        k_He4_Ca40_to_p_Sc43_derived_removed = 229,
        k_Ti44_to_He4_Ca40_derived_removed = 230,
        k_Ti44_to_p_Sc43_derived_removed = 231,
        k_Cr48_to_He4_Ti44_derived_removed = 232,
        k_Cr48_to_p_V47_derived_removed = 233,
        k_p_V47_to_He4_Ti44_derived_removed = 234,
        k_Fe52_to_He4_Cr48_derived_removed = 235,
        k_Fe52_to_p_Mn51_derived_removed = 236,
        k_p_Mn51_to_He4_Cr48_derived_removed = 237,
        k_Ni56_to_He4_Fe52_derived_removed = 238,
        k_Ni56_to_p_Co55_derived_removed = 239,
        k_p_Co55_to_He4_Fe52_derived_removed = 240,
        NumRates = k_p_Co55_to_He4_Fe52_derived_removed
    };

    // number of reaclib rates

    const int NrateReaclib = 214;

    // number of tabular rates

    const int NrateTabular = 14;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "N13_to_C13_weak_wc12",  // 1,
        "O14_to_N14_weak_wc12",  // 2,
        "O15_to_N15_weak_wc12",  // 3,
        "n_C12_to_C13",  // 4,
        "p_C12_to_N13",  // 5,
        "He4_C12_to_O16",  // 6,
        "p_C13_to_N14",  // 7,
        "n_N13_to_N14",  // 8,
        "p_N13_to_O14",  // 9,
        "n_N14_to_N15",  // 10,
        "p_N14_to_O15",  // 11,
        "He4_N14_to_F18",  // 12,
        "p_N15_to_O16",  // 13,
        "He4_N15_to_F19",  // 14,
        "n_O14_to_O15",  // 15,
        "He4_O14_to_Ne18",  // 16,
        "n_O15_to_O16",  // 17,
        "He4_O15_to_Ne19",  // 18,
        "n_O16_to_O17",  // 19,
        "p_O16_to_F17",  // 20,
        "He4_O16_to_Ne20",  // 21,
        "n_O17_to_O18",  // 22,
        "p_O17_to_F18",  // 23,
        "He4_O17_to_Ne21",  // 24,
        "p_O18_to_F19",  // 25,
        "He4_O18_to_Ne22",  // 26,
        "n_F17_to_F18",  // 27,
        "p_F17_to_Ne18",  // 28,
        "He4_F17_to_Na21",  // 29,
        "n_F18_to_F19",  // 30,
        "p_F18_to_Ne19",  // 31,
        "He4_F18_to_Na22",  // 32,
        "p_F19_to_Ne20",  // 33,
        "He4_F19_to_Na23",  // 34,
        "n_Ne18_to_Ne19",  // 35,
        "n_Ne19_to_Ne20",  // 36,
        "n_Ne20_to_Ne21",  // 37,
        "p_Ne20_to_Na21",  // 38,
        "He4_Ne20_to_Mg24",  // 39,
        "n_Ne21_to_Ne22",  // 40,
        "p_Ne21_to_Na22",  // 41,
        "p_Ne22_to_Na23",  // 42,
        "n_Na21_to_Na22",  // 43,
        "n_Na22_to_Na23",  // 44,
        "p_Na23_to_Mg24",  // 45,
        "He4_Mg24_to_Si28",  // 46,
        "p_Al27_to_Si28",  // 47,
        "He4_Si28_to_S32",  // 48,
        "p_P31_to_S32",  // 49,
        "C12_C12_to_p_Na23",  // 50,
        "C12_C12_to_He4_Ne20",  // 51,
        "p_C13_to_n_N13",  // 52,
        "He4_C13_to_n_O16",  // 53,
        "He4_N13_to_p_O16",  // 54,
        "p_N14_to_n_O14",  // 55,
        "He4_N14_to_n_F17",  // 56,
        "p_N15_to_n_O15",  // 57,
        "p_N15_to_He4_C12",  // 58,
        "He4_O14_to_p_F17",  // 59,
        "n_O15_to_He4_C12",  // 60,
        "C12_O16_to_p_Al27",  // 61,
        "C12_O16_to_He4_Mg24",  // 62,
        "O16_O16_to_p_P31",  // 63,
        "O16_O16_to_He4_Si28",  // 64,
        "p_O17_to_He4_N14",  // 65,
        "He4_O17_to_n_Ne20",  // 66,
        "p_O18_to_He4_N15",  // 67,
        "He4_O18_to_n_Ne21",  // 68,
        "n_F17_to_p_O17",  // 69,
        "n_F18_to_p_O18",  // 70,
        "n_F18_to_He4_N15",  // 71,
        "p_F18_to_He4_O15",  // 72,
        "He4_F18_to_p_Ne21",  // 73,
        "p_F19_to_n_Ne19",  // 74,
        "p_F19_to_He4_O16",  // 75,
        "He4_F19_to_p_Ne22",  // 76,
        "n_Ne18_to_p_F18",  // 77,
        "n_Ne18_to_He4_O15",  // 78,
        "He4_Ne18_to_p_Na21",  // 79,
        "n_Ne19_to_He4_O16",  // 80,
        "He4_Ne19_to_p_Na22",  // 81,
        "p_Ne20_to_He4_F17",  // 82,
        "He4_Ne21_to_n_Mg24",  // 83,
        "n_Na21_to_p_Ne21",  // 84,
        "n_Na21_to_He4_F18",  // 85,
        "n_Na22_to_p_Ne22",  // 86,
        "n_Na22_to_He4_F19",  // 87,
        "p_Na23_to_He4_Ne20",  // 88,
        "p_Mg24_to_He4_Na21",  // 89,
        "p_Al27_to_He4_Mg24",  // 90,
        "p_P31_to_He4_Si28",  // 91,
        "He4_He4_He4_to_C12",  // 92,
        "C12_C12_to_Mg24_modified",  // 93,
        "O16_O16_to_S32_modified",  // 94,
        "C12_O16_to_Si28_modified",  // 95,
        "He4_S32_to_Ar36_removed",  // 96,
        "p_Cl35_to_Ar36_removed",  // 97,
        "p_Cl35_to_He4_S32_removed",  // 98,
        "He4_Ar36_to_Ca40_removed",  // 99,
        "p_K39_to_Ca40_removed",  // 100,
        "p_K39_to_He4_Ar36_removed",  // 101,
        "He4_Ca40_to_Ti44_removed",  // 102,
        "p_Sc43_to_Ti44_removed",  // 103,
        "p_Sc43_to_He4_Ca40_removed",  // 104,
        "He4_Ti44_to_Cr48_removed",  // 105,
        "He4_Ti44_to_p_V47_removed",  // 106,
        "p_V47_to_Cr48_removed",  // 107,
        "He4_Cr48_to_Fe52_removed",  // 108,
        "He4_Cr48_to_p_Mn51_removed",  // 109,
        "p_Mn51_to_Fe52_removed",  // 110,
        "He4_Fe52_to_Ni56_removed",  // 111,
        "He4_Fe52_to_p_Co55_removed",  // 112,
        "p_Co55_to_Ni56_removed",  // 113,
        "n_to_p",  // 114,
        "p_to_n",  // 115,
        "F17_to_O17",  // 116,
        "O17_to_F17",  // 117,
        "F18_to_Ne18",  // 118,
        "F18_to_O18",  // 119,
        "Ne18_to_F18",  // 120,
        "O18_to_F18",  // 121,
        "F19_to_Ne19",  // 122,
        "Ne19_to_F19",  // 123,
        "Na21_to_Ne21",  // 124,
        "Ne21_to_Na21",  // 125,
        "Na22_to_Ne22",  // 126,
        "Ne22_to_Na22",  // 127,
        "S32_He4_to_Ar36_approx",  // 128,
        "Ar36_to_S32_He4_approx",  // 129,
        "Ar36_He4_to_Ca40_approx",  // 130,
        "Ca40_to_Ar36_He4_approx",  // 131,
        "Ca40_He4_to_Ti44_approx",  // 132,
        "Ti44_to_Ca40_He4_approx",  // 133,
        "Ti44_He4_to_Cr48_approx",  // 134,
        "Cr48_to_Ti44_He4_approx",  // 135,
        "Cr48_He4_to_Fe52_approx",  // 136,
        "Fe52_to_Cr48_He4_approx",  // 137,
        "Fe52_He4_to_Ni56_approx",  // 138,
        "Ni56_to_Fe52_He4_approx",  // 139,
        "C13_to_n_C12_derived",  // 140,
        "N13_to_p_C12_derived",  // 141,
        "N14_to_n_N13_derived",  // 142,
        "N14_to_p_C13_derived",  // 143,
        "N15_to_n_N14_derived",  // 144,
        "O14_to_p_N13_derived",  // 145,
        "O15_to_n_O14_derived",  // 146,
        "O15_to_p_N14_derived",  // 147,
        "O16_to_n_O15_derived",  // 148,
        "O16_to_p_N15_derived",  // 149,
        "O16_to_He4_C12_derived",  // 150,
        "O17_to_n_O16_derived",  // 151,
        "O18_to_n_O17_derived",  // 152,
        "F17_to_p_O16_derived",  // 153,
        "F18_to_n_F17_derived",  // 154,
        "F18_to_p_O17_derived",  // 155,
        "F18_to_He4_N14_derived",  // 156,
        "F19_to_n_F18_derived",  // 157,
        "F19_to_p_O18_derived",  // 158,
        "F19_to_He4_N15_derived",  // 159,
        "Ne18_to_p_F17_derived",  // 160,
        "Ne18_to_He4_O14_derived",  // 161,
        "Ne19_to_n_Ne18_derived",  // 162,
        "Ne19_to_p_F18_derived",  // 163,
        "Ne19_to_He4_O15_derived",  // 164,
        "Ne20_to_n_Ne19_derived",  // 165,
        "Ne20_to_p_F19_derived",  // 166,
        "Ne20_to_He4_O16_derived",  // 167,
        "Ne21_to_n_Ne20_derived",  // 168,
        "Ne21_to_He4_O17_derived",  // 169,
        "Ne22_to_n_Ne21_derived",  // 170,
        "Ne22_to_He4_O18_derived",  // 171,
        "Na21_to_p_Ne20_derived",  // 172,
        "Na21_to_He4_F17_derived",  // 173,
        "Na22_to_n_Na21_derived",  // 174,
        "Na22_to_p_Ne21_derived",  // 175,
        "Na22_to_He4_F18_derived",  // 176,
        "Na23_to_n_Na22_derived",  // 177,
        "Na23_to_p_Ne22_derived",  // 178,
        "Na23_to_He4_F19_derived",  // 179,
        "Mg24_to_p_Na23_derived",  // 180,
        "Mg24_to_He4_Ne20_derived",  // 181,
        "Si28_to_p_Al27_derived",  // 182,
        "Si28_to_He4_Mg24_derived",  // 183,
        "S32_to_p_P31_derived",  // 184,
        "S32_to_He4_Si28_derived",  // 185,
        "C12_to_He4_He4_He4_derived",  // 186,
        "He4_C12_to_n_O15_derived",  // 187,
        "He4_C12_to_p_N15_derived",  // 188,
        "n_N13_to_p_C13_derived",  // 189,
        "He4_N14_to_p_O17_derived",  // 190,
        "He4_N15_to_n_F18_derived",  // 191,
        "He4_N15_to_p_O18_derived",  // 192,
        "n_O14_to_p_N14_derived",  // 193,
        "n_O15_to_p_N15_derived",  // 194,
        "He4_O15_to_n_Ne18_derived",  // 195,
        "He4_O15_to_p_F18_derived",  // 196,
        "n_O16_to_He4_C13_derived",  // 197,
        "p_O16_to_He4_N13_derived",  // 198,
        "He4_O16_to_n_Ne19_derived",  // 199,
        "He4_O16_to_p_F19_derived",  // 200,
        "p_O17_to_n_F17_derived",  // 201,
        "p_O18_to_n_F18_derived",  // 202,
        "n_F17_to_He4_N14_derived",  // 203,
        "p_F17_to_He4_O14_derived",  // 204,
        "He4_F17_to_p_Ne20_derived",  // 205,
        "p_F18_to_n_Ne18_derived",  // 206,
        "He4_F18_to_n_Na21_derived",  // 207,
        "He4_F19_to_n_Na22_derived",  // 208,
        "n_Ne19_to_p_F19_derived",  // 209,
        "n_Ne20_to_He4_O17_derived",  // 210,
        "He4_Ne20_to_p_Na23_derived",  // 211,
        "n_Ne21_to_He4_O18_derived",  // 212,
        "p_Ne21_to_n_Na21_derived",  // 213,
        "p_Ne21_to_He4_F18_derived",  // 214,
        "p_Ne22_to_n_Na22_derived",  // 215,
        "p_Ne22_to_He4_F19_derived",  // 216,
        "p_Na21_to_He4_Ne18_derived",  // 217,
        "He4_Na21_to_p_Mg24_derived",  // 218,
        "p_Na22_to_He4_Ne19_derived",  // 219,
        "n_Mg24_to_He4_Ne21_derived",  // 220,
        "He4_Mg24_to_p_Al27_derived",  // 221,
        "He4_Si28_to_p_P31_derived",  // 222,
        "He4_S32_to_p_Cl35_derived_removed",  // 223,
        "Ar36_to_He4_S32_derived_removed",  // 224,
        "Ar36_to_p_Cl35_derived_removed",  // 225,
        "He4_Ar36_to_p_K39_derived_removed",  // 226,
        "Ca40_to_He4_Ar36_derived_removed",  // 227,
        "Ca40_to_p_K39_derived_removed",  // 228,
        "He4_Ca40_to_p_Sc43_derived_removed",  // 229,
        "Ti44_to_He4_Ca40_derived_removed",  // 230,
        "Ti44_to_p_Sc43_derived_removed",  // 231,
        "Cr48_to_He4_Ti44_derived_removed",  // 232,
        "Cr48_to_p_V47_derived_removed",  // 233,
        "p_V47_to_He4_Ti44_derived_removed",  // 234,
        "Fe52_to_He4_Cr48_derived_removed",  // 235,
        "Fe52_to_p_Mn51_derived_removed",  // 236,
        "p_Mn51_to_He4_Cr48_derived_removed",  // 237,
        "Ni56_to_He4_Fe52_derived_removed",  // 238,
        "Ni56_to_p_Co55_derived_removed",  // 239,
        "p_Co55_to_He4_Fe52_derived_removed"  // 240,
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

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, amrex::Order::C> rate_indices;
}
#endif

#endif
