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

        if constexpr (spec == H1) {
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
        k_N13_to_C13_reaclib = 1,
        k_O14_to_N14_reaclib = 2,
        k_O15_to_N15_reaclib = 3,
        k_p_C12_to_N13_reaclib = 4,
        k_He4_C12_to_O16_reaclib = 5,
        k_p_C13_to_N14_reaclib = 6,
        k_p_N13_to_O14_reaclib = 7,
        k_p_N14_to_O15_reaclib = 8,
        k_He4_N14_to_F18_reaclib = 9,
        k_p_N15_to_O16_reaclib = 10,
        k_He4_N15_to_F19_reaclib = 11,
        k_He4_O14_to_Ne18_reaclib = 12,
        k_He4_O15_to_Ne19_reaclib = 13,
        k_p_O16_to_F17_reaclib = 14,
        k_He4_O16_to_Ne20_reaclib = 15,
        k_p_O17_to_F18_reaclib = 16,
        k_He4_O17_to_Ne21_reaclib = 17,
        k_p_O18_to_F19_reaclib = 18,
        k_He4_O18_to_Ne22_reaclib = 19,
        k_p_F17_to_Ne18_reaclib = 20,
        k_He4_F17_to_Na21_reaclib = 21,
        k_p_F18_to_Ne19_reaclib = 22,
        k_He4_F18_to_Na22_reaclib = 23,
        k_p_F19_to_Ne20_reaclib = 24,
        k_He4_F19_to_Na23_reaclib = 25,
        k_p_Ne20_to_Na21_reaclib = 26,
        k_He4_Ne20_to_Mg24_reaclib = 27,
        k_p_Ne21_to_Na22_reaclib = 28,
        k_p_Ne22_to_Na23_reaclib = 29,
        k_p_Na23_to_Mg24_reaclib = 30,
        k_He4_Na23_to_Al27_reaclib = 31,
        k_He4_Mg24_to_Si28_reaclib = 32,
        k_p_Al27_to_Si28_reaclib = 33,
        k_He4_Al27_to_P31_reaclib = 34,
        k_He4_Si28_to_S32_reaclib = 35,
        k_p_P31_to_S32_reaclib = 36,
        k_C12_C12_to_p_Na23_reaclib = 37,
        k_C12_C12_to_He4_Ne20_reaclib = 38,
        k_He4_N13_to_p_O16_reaclib = 39,
        k_p_N15_to_He4_C12_reaclib = 40,
        k_He4_O14_to_p_F17_reaclib = 41,
        k_C12_O16_to_p_Al27_reaclib = 42,
        k_C12_O16_to_He4_Mg24_reaclib = 43,
        k_O16_O16_to_p_P31_reaclib = 44,
        k_O16_O16_to_He4_Si28_reaclib = 45,
        k_p_O17_to_He4_N14_reaclib = 46,
        k_p_O18_to_He4_N15_reaclib = 47,
        k_p_F18_to_He4_O15_reaclib = 48,
        k_He4_F18_to_p_Ne21_reaclib = 49,
        k_p_F19_to_He4_O16_reaclib = 50,
        k_He4_F19_to_p_Ne22_reaclib = 51,
        k_He4_Ne18_to_p_Na21_reaclib = 52,
        k_He4_Ne19_to_p_Na22_reaclib = 53,
        k_p_Ne20_to_He4_F17_reaclib = 54,
        k_p_Na23_to_He4_Ne20_reaclib = 55,
        k_p_Mg24_to_He4_Na21_reaclib = 56,
        k_p_Al27_to_He4_Mg24_reaclib = 57,
        k_p_P31_to_He4_Si28_reaclib = 58,
        k_He4_He4_He4_to_C12_reaclib = 59,
        k_C12_C12_to_n_Mg23_removed = 60,
        k_O16_O16_to_n_S31_removed = 61,
        k_C12_O16_to_n_Si27_removed = 62,
        k_He4_Na22_to_Al26_removed = 63,
        k_p_Mg24_to_Al25_removed = 64,
        k_He4_S32_to_Ar36_removed = 65,
        k_p_Cl35_to_Ar36_removed = 66,
        k_p_Cl35_to_He4_S32_removed = 67,
        k_He4_Ar36_to_Ca40_removed = 68,
        k_p_K39_to_Ca40_removed = 69,
        k_p_K39_to_He4_Ar36_removed = 70,
        k_He4_Ca40_to_Ti44_removed = 71,
        k_p_Sc43_to_Ti44_removed = 72,
        k_p_Sc43_to_He4_Ca40_removed = 73,
        k_He4_Ti44_to_Cr48_removed = 74,
        k_He4_Ti44_to_p_V47_removed = 75,
        k_p_V47_to_Cr48_removed = 76,
        k_He4_Cr48_to_Fe52_removed = 77,
        k_He4_Cr48_to_p_Mn51_removed = 78,
        k_p_Mn51_to_Fe52_removed = 79,
        k_He4_Fe52_to_Ni56_removed = 80,
        k_He4_Fe52_to_p_Co55_removed = 81,
        k_p_Co55_to_Ni56_removed = 82,
        k_F17_to_O17_weaktab = 83,
        k_O17_to_F17_weaktab = 84,
        k_F18_to_Ne18_weaktab = 85,
        k_F18_to_O18_weaktab = 86,
        k_Ne18_to_F18_weaktab = 87,
        k_O18_to_F18_weaktab = 88,
        k_F19_to_Ne19_weaktab = 89,
        k_Ne19_to_F19_weaktab = 90,
        k_Na21_to_Ne21_weaktab = 91,
        k_Ne21_to_Na21_weaktab = 92,
        k_Na22_to_Ne22_weaktab = 93,
        k_Ne22_to_Na22_weaktab = 94,
        k_S32_He4_to_Ar36_approx = 95,
        k_Ar36_to_S32_He4_approx = 96,
        k_Ar36_He4_to_Ca40_approx = 97,
        k_Ca40_to_Ar36_He4_approx = 98,
        k_Ca40_He4_to_Ti44_approx = 99,
        k_Ti44_to_Ca40_He4_approx = 100,
        k_Ti44_He4_to_Cr48_approx = 101,
        k_Cr48_to_Ti44_He4_approx = 102,
        k_Cr48_He4_to_Fe52_approx = 103,
        k_Fe52_to_Cr48_He4_approx = 104,
        k_Fe52_He4_to_Ni56_approx = 105,
        k_Ni56_to_Fe52_He4_approx = 106,
        k_C12_C12_to_Mg24_modified = 107,
        k_O16_O16_to_S32_modified = 108,
        k_C12_O16_to_Si28_modified = 109,
        k_He4_Na22_to_Si28_modified = 110,
        k_p_Mg24_to_Na22_He4_modified = 111,
        k_N13_to_p_C12_derived = 112,
        k_N14_to_p_C13_derived = 113,
        k_O14_to_p_N13_derived = 114,
        k_O15_to_p_N14_derived = 115,
        k_O16_to_p_N15_derived = 116,
        k_O16_to_He4_C12_derived = 117,
        k_F17_to_p_O16_derived = 118,
        k_F18_to_p_O17_derived = 119,
        k_F18_to_He4_N14_derived = 120,
        k_F19_to_p_O18_derived = 121,
        k_F19_to_He4_N15_derived = 122,
        k_Ne18_to_p_F17_derived = 123,
        k_Ne18_to_He4_O14_derived = 124,
        k_Ne19_to_p_F18_derived = 125,
        k_Ne19_to_He4_O15_derived = 126,
        k_Ne20_to_p_F19_derived = 127,
        k_Ne20_to_He4_O16_derived = 128,
        k_Ne21_to_He4_O17_derived = 129,
        k_Ne22_to_He4_O18_derived = 130,
        k_Na21_to_p_Ne20_derived = 131,
        k_Na21_to_He4_F17_derived = 132,
        k_Na22_to_p_Ne21_derived = 133,
        k_Na22_to_He4_F18_derived = 134,
        k_Na23_to_p_Ne22_derived = 135,
        k_Na23_to_He4_F19_derived = 136,
        k_Mg24_to_p_Na23_derived = 137,
        k_Mg24_to_He4_Ne20_derived = 138,
        k_Al27_to_He4_Na23_derived = 139,
        k_Si28_to_p_Al27_derived = 140,
        k_Si28_to_He4_Mg24_derived = 141,
        k_P31_to_He4_Al27_derived = 142,
        k_S32_to_p_P31_derived = 143,
        k_S32_to_He4_Si28_derived = 144,
        k_C12_to_He4_He4_He4_derived = 145,
        k_He4_C12_to_p_N15_derived = 146,
        k_He4_N14_to_p_O17_derived = 147,
        k_He4_N15_to_p_O18_derived = 148,
        k_He4_O15_to_p_F18_derived = 149,
        k_p_O16_to_He4_N13_derived = 150,
        k_He4_O16_to_p_F19_derived = 151,
        k_p_F17_to_He4_O14_derived = 152,
        k_He4_F17_to_p_Ne20_derived = 153,
        k_He4_Ne20_to_p_Na23_derived = 154,
        k_p_Ne21_to_He4_F18_derived = 155,
        k_p_Ne22_to_He4_F19_derived = 156,
        k_p_Na21_to_He4_Ne18_derived = 157,
        k_He4_Na21_to_p_Mg24_derived = 158,
        k_p_Na22_to_He4_Ne19_derived = 159,
        k_He4_Mg24_to_p_Al27_derived = 160,
        k_He4_Si28_to_p_P31_derived = 161,
        k_He4_S32_to_p_Cl35_removed = 162,
        k_Ar36_to_He4_S32_removed = 163,
        k_Ar36_to_p_Cl35_removed = 164,
        k_He4_Ar36_to_p_K39_removed = 165,
        k_Ca40_to_He4_Ar36_removed = 166,
        k_Ca40_to_p_K39_removed = 167,
        k_He4_Ca40_to_p_Sc43_removed = 168,
        k_Ti44_to_He4_Ca40_removed = 169,
        k_Ti44_to_p_Sc43_removed = 170,
        k_Cr48_to_He4_Ti44_removed = 171,
        k_Cr48_to_p_V47_removed = 172,
        k_p_V47_to_He4_Ti44_removed = 173,
        k_Fe52_to_He4_Cr48_removed = 174,
        k_Fe52_to_p_Mn51_removed = 175,
        k_p_Mn51_to_He4_Cr48_removed = 176,
        k_Ni56_to_He4_Fe52_removed = 177,
        k_Ni56_to_p_Co55_removed = 178,
        k_p_Co55_to_He4_Fe52_removed = 179,
        NumRates = k_p_Co55_to_He4_Fe52_removed
    };

    enum ScreenPairs :
    std::uint8_t
    {
        k_p_C12 = 1,
        k_He4_C12 = 2,
        k_p_C13 = 3,
        k_p_N13 = 4,
        k_p_N14 = 5,
        k_He4_N14 = 6,
        k_p_N15 = 7,
        k_He4_N15 = 8,
        k_He4_O14 = 9,
        k_He4_O15 = 10,
        k_p_O16 = 11,
        k_He4_O16 = 12,
        k_p_O17 = 13,
        k_He4_O17 = 14,
        k_p_O18 = 15,
        k_He4_O18 = 16,
        k_p_F17 = 17,
        k_He4_F17 = 18,
        k_p_F18 = 19,
        k_He4_F18 = 20,
        k_p_F19 = 21,
        k_He4_F19 = 22,
        k_p_Ne20 = 23,
        k_He4_Ne20 = 24,
        k_p_Ne21 = 25,
        k_p_Ne22 = 26,
        k_p_Na23 = 27,
        k_He4_Na23 = 28,
        k_He4_Mg24 = 29,
        k_p_Al27 = 30,
        k_He4_Al27 = 31,
        k_He4_Si28 = 32,
        k_p_P31 = 33,
        k_C12_C12 = 34,
        k_He4_N13 = 35,
        k_C12_O16 = 36,
        k_O16_O16 = 37,
        k_He4_Ne18 = 38,
        k_He4_Ne19 = 39,
        k_p_Mg24 = 40,
        k_He4_He4 = 41,
        k_He4_Be8 = 42,
        k_p_Na21 = 43,
        k_He4_Na21 = 44,
        k_p_Na22 = 45,
        k_He4_Na22 = 46,
        k_He4_S32 = 47,
        k_p_Cl35 = 48,
        k_He4_Ar36 = 49,
        k_p_K39 = 50,
        k_He4_Ca40 = 51,
        k_p_Sc43 = 52,
        k_He4_Ti44 = 53,
        k_p_V47 = 54,
        k_He4_Cr48 = 55,
        k_p_Mn51 = 56,
        k_He4_Fe52 = 57,
        k_p_Co55 = 58,
        NumScreenPairs = k_p_Co55
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "N13_to_C13_reaclib",  // 1,
        "O14_to_N14_reaclib",  // 2,
        "O15_to_N15_reaclib",  // 3,
        "p_C12_to_N13_reaclib",  // 4,
        "He4_C12_to_O16_reaclib",  // 5,
        "p_C13_to_N14_reaclib",  // 6,
        "p_N13_to_O14_reaclib",  // 7,
        "p_N14_to_O15_reaclib",  // 8,
        "He4_N14_to_F18_reaclib",  // 9,
        "p_N15_to_O16_reaclib",  // 10,
        "He4_N15_to_F19_reaclib",  // 11,
        "He4_O14_to_Ne18_reaclib",  // 12,
        "He4_O15_to_Ne19_reaclib",  // 13,
        "p_O16_to_F17_reaclib",  // 14,
        "He4_O16_to_Ne20_reaclib",  // 15,
        "p_O17_to_F18_reaclib",  // 16,
        "He4_O17_to_Ne21_reaclib",  // 17,
        "p_O18_to_F19_reaclib",  // 18,
        "He4_O18_to_Ne22_reaclib",  // 19,
        "p_F17_to_Ne18_reaclib",  // 20,
        "He4_F17_to_Na21_reaclib",  // 21,
        "p_F18_to_Ne19_reaclib",  // 22,
        "He4_F18_to_Na22_reaclib",  // 23,
        "p_F19_to_Ne20_reaclib",  // 24,
        "He4_F19_to_Na23_reaclib",  // 25,
        "p_Ne20_to_Na21_reaclib",  // 26,
        "He4_Ne20_to_Mg24_reaclib",  // 27,
        "p_Ne21_to_Na22_reaclib",  // 28,
        "p_Ne22_to_Na23_reaclib",  // 29,
        "p_Na23_to_Mg24_reaclib",  // 30,
        "He4_Na23_to_Al27_reaclib",  // 31,
        "He4_Mg24_to_Si28_reaclib",  // 32,
        "p_Al27_to_Si28_reaclib",  // 33,
        "He4_Al27_to_P31_reaclib",  // 34,
        "He4_Si28_to_S32_reaclib",  // 35,
        "p_P31_to_S32_reaclib",  // 36,
        "C12_C12_to_p_Na23_reaclib",  // 37,
        "C12_C12_to_He4_Ne20_reaclib",  // 38,
        "He4_N13_to_p_O16_reaclib",  // 39,
        "p_N15_to_He4_C12_reaclib",  // 40,
        "He4_O14_to_p_F17_reaclib",  // 41,
        "C12_O16_to_p_Al27_reaclib",  // 42,
        "C12_O16_to_He4_Mg24_reaclib",  // 43,
        "O16_O16_to_p_P31_reaclib",  // 44,
        "O16_O16_to_He4_Si28_reaclib",  // 45,
        "p_O17_to_He4_N14_reaclib",  // 46,
        "p_O18_to_He4_N15_reaclib",  // 47,
        "p_F18_to_He4_O15_reaclib",  // 48,
        "He4_F18_to_p_Ne21_reaclib",  // 49,
        "p_F19_to_He4_O16_reaclib",  // 50,
        "He4_F19_to_p_Ne22_reaclib",  // 51,
        "He4_Ne18_to_p_Na21_reaclib",  // 52,
        "He4_Ne19_to_p_Na22_reaclib",  // 53,
        "p_Ne20_to_He4_F17_reaclib",  // 54,
        "p_Na23_to_He4_Ne20_reaclib",  // 55,
        "p_Mg24_to_He4_Na21_reaclib",  // 56,
        "p_Al27_to_He4_Mg24_reaclib",  // 57,
        "p_P31_to_He4_Si28_reaclib",  // 58,
        "He4_He4_He4_to_C12_reaclib",  // 59,
        "C12_C12_to_n_Mg23_removed",  // 60,
        "O16_O16_to_n_S31_removed",  // 61,
        "C12_O16_to_n_Si27_removed",  // 62,
        "He4_Na22_to_Al26_removed",  // 63,
        "p_Mg24_to_Al25_removed",  // 64,
        "He4_S32_to_Ar36_removed",  // 65,
        "p_Cl35_to_Ar36_removed",  // 66,
        "p_Cl35_to_He4_S32_removed",  // 67,
        "He4_Ar36_to_Ca40_removed",  // 68,
        "p_K39_to_Ca40_removed",  // 69,
        "p_K39_to_He4_Ar36_removed",  // 70,
        "He4_Ca40_to_Ti44_removed",  // 71,
        "p_Sc43_to_Ti44_removed",  // 72,
        "p_Sc43_to_He4_Ca40_removed",  // 73,
        "He4_Ti44_to_Cr48_removed",  // 74,
        "He4_Ti44_to_p_V47_removed",  // 75,
        "p_V47_to_Cr48_removed",  // 76,
        "He4_Cr48_to_Fe52_removed",  // 77,
        "He4_Cr48_to_p_Mn51_removed",  // 78,
        "p_Mn51_to_Fe52_removed",  // 79,
        "He4_Fe52_to_Ni56_removed",  // 80,
        "He4_Fe52_to_p_Co55_removed",  // 81,
        "p_Co55_to_Ni56_removed",  // 82,
        "F17_to_O17_weaktab",  // 83,
        "O17_to_F17_weaktab",  // 84,
        "F18_to_Ne18_weaktab",  // 85,
        "F18_to_O18_weaktab",  // 86,
        "Ne18_to_F18_weaktab",  // 87,
        "O18_to_F18_weaktab",  // 88,
        "F19_to_Ne19_weaktab",  // 89,
        "Ne19_to_F19_weaktab",  // 90,
        "Na21_to_Ne21_weaktab",  // 91,
        "Ne21_to_Na21_weaktab",  // 92,
        "Na22_to_Ne22_weaktab",  // 93,
        "Ne22_to_Na22_weaktab",  // 94,
        "S32_He4_to_Ar36_approx",  // 95,
        "Ar36_to_S32_He4_approx",  // 96,
        "Ar36_He4_to_Ca40_approx",  // 97,
        "Ca40_to_Ar36_He4_approx",  // 98,
        "Ca40_He4_to_Ti44_approx",  // 99,
        "Ti44_to_Ca40_He4_approx",  // 100,
        "Ti44_He4_to_Cr48_approx",  // 101,
        "Cr48_to_Ti44_He4_approx",  // 102,
        "Cr48_He4_to_Fe52_approx",  // 103,
        "Fe52_to_Cr48_He4_approx",  // 104,
        "Fe52_He4_to_Ni56_approx",  // 105,
        "Ni56_to_Fe52_He4_approx",  // 106,
        "C12_C12_to_Mg24_modified",  // 107,
        "O16_O16_to_S32_modified",  // 108,
        "C12_O16_to_Si28_modified",  // 109,
        "He4_Na22_to_Si28_modified",  // 110,
        "p_Mg24_to_Na22_He4_modified",  // 111,
        "N13_to_p_C12_derived",  // 112,
        "N14_to_p_C13_derived",  // 113,
        "O14_to_p_N13_derived",  // 114,
        "O15_to_p_N14_derived",  // 115,
        "O16_to_p_N15_derived",  // 116,
        "O16_to_He4_C12_derived",  // 117,
        "F17_to_p_O16_derived",  // 118,
        "F18_to_p_O17_derived",  // 119,
        "F18_to_He4_N14_derived",  // 120,
        "F19_to_p_O18_derived",  // 121,
        "F19_to_He4_N15_derived",  // 122,
        "Ne18_to_p_F17_derived",  // 123,
        "Ne18_to_He4_O14_derived",  // 124,
        "Ne19_to_p_F18_derived",  // 125,
        "Ne19_to_He4_O15_derived",  // 126,
        "Ne20_to_p_F19_derived",  // 127,
        "Ne20_to_He4_O16_derived",  // 128,
        "Ne21_to_He4_O17_derived",  // 129,
        "Ne22_to_He4_O18_derived",  // 130,
        "Na21_to_p_Ne20_derived",  // 131,
        "Na21_to_He4_F17_derived",  // 132,
        "Na22_to_p_Ne21_derived",  // 133,
        "Na22_to_He4_F18_derived",  // 134,
        "Na23_to_p_Ne22_derived",  // 135,
        "Na23_to_He4_F19_derived",  // 136,
        "Mg24_to_p_Na23_derived",  // 137,
        "Mg24_to_He4_Ne20_derived",  // 138,
        "Al27_to_He4_Na23_derived",  // 139,
        "Si28_to_p_Al27_derived",  // 140,
        "Si28_to_He4_Mg24_derived",  // 141,
        "P31_to_He4_Al27_derived",  // 142,
        "S32_to_p_P31_derived",  // 143,
        "S32_to_He4_Si28_derived",  // 144,
        "C12_to_He4_He4_He4_derived",  // 145,
        "He4_C12_to_p_N15_derived",  // 146,
        "He4_N14_to_p_O17_derived",  // 147,
        "He4_N15_to_p_O18_derived",  // 148,
        "He4_O15_to_p_F18_derived",  // 149,
        "p_O16_to_He4_N13_derived",  // 150,
        "He4_O16_to_p_F19_derived",  // 151,
        "p_F17_to_He4_O14_derived",  // 152,
        "He4_F17_to_p_Ne20_derived",  // 153,
        "He4_Ne20_to_p_Na23_derived",  // 154,
        "p_Ne21_to_He4_F18_derived",  // 155,
        "p_Ne22_to_He4_F19_derived",  // 156,
        "p_Na21_to_He4_Ne18_derived",  // 157,
        "He4_Na21_to_p_Mg24_derived",  // 158,
        "p_Na22_to_He4_Ne19_derived",  // 159,
        "He4_Mg24_to_p_Al27_derived",  // 160,
        "He4_Si28_to_p_P31_derived",  // 161,
        "He4_S32_to_p_Cl35_removed",  // 162,
        "Ar36_to_He4_S32_removed",  // 163,
        "Ar36_to_p_Cl35_removed",  // 164,
        "He4_Ar36_to_p_K39_removed",  // 165,
        "Ca40_to_He4_Ar36_removed",  // 166,
        "Ca40_to_p_K39_removed",  // 167,
        "He4_Ca40_to_p_Sc43_removed",  // 168,
        "Ti44_to_He4_Ca40_removed",  // 169,
        "Ti44_to_p_Sc43_removed",  // 170,
        "Cr48_to_He4_Ti44_removed",  // 171,
        "Cr48_to_p_V47_removed",  // 172,
        "p_V47_to_He4_Ti44_removed",  // 173,
        "Fe52_to_He4_Cr48_removed",  // 174,
        "Fe52_to_p_Mn51_removed",  // 175,
        "p_Mn51_to_He4_Cr48_removed",  // 176,
        "Ni56_to_He4_Fe52_removed",  // 177,
        "Ni56_to_p_Co55_removed",  // 178,
        "p_Co55_to_He4_Fe52_removed"  // 179,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    // A RatePair data table, 2D array of shape (NumNSERatePairs, 10).
    // Each row corresponds to a RatePair returned by _get_nse_rate_pairs().
    // Use enum names following NetworkSpecies and NetworkRates.

    //  - Columns 1–3: 1-based index of reactant nuclei (forward rate)
    //  - Columns 4–6: 1-based index of product nuclei (forward rate)
    //  - Columns 7-8: 1-based index of non (n, H1, He4) nuclei.
    //  - Column  9:   1-based index of the forward rate
    //  - Column  10:  1-based index of the reverse rate

    using namespace Species;
    using namespace Rates;

    constexpr int H1_index = 1;
    constexpr int N_index = -1;
    constexpr int He4_index = 2;
    constexpr int NumNSERatePairs = 56;

    inline AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, NumNSERatePairs, 1, 10, amrex::Order::C> rate_pair_data {
        H1, C12, -1, N13, -1, -1, C12, N13, k_p_C12_to_N13_reaclib, k_N13_to_p_C12_derived,
        He4, C12, -1, O16, -1, -1, C12, O16, k_He4_C12_to_O16_reaclib, k_O16_to_He4_C12_derived,
        H1, C13, -1, N14, -1, -1, C13, N14, k_p_C13_to_N14_reaclib, k_N14_to_p_C13_derived,
        H1, N13, -1, O14, -1, -1, N13, O14, k_p_N13_to_O14_reaclib, k_O14_to_p_N13_derived,
        H1, N14, -1, O15, -1, -1, N14, O15, k_p_N14_to_O15_reaclib, k_O15_to_p_N14_derived,
        He4, N14, -1, F18, -1, -1, N14, F18, k_He4_N14_to_F18_reaclib, k_F18_to_He4_N14_derived,
        H1, N15, -1, O16, -1, -1, N15, O16, k_p_N15_to_O16_reaclib, k_O16_to_p_N15_derived,
        He4, N15, -1, F19, -1, -1, N15, F19, k_He4_N15_to_F19_reaclib, k_F19_to_He4_N15_derived,
        He4, O14, -1, Ne18, -1, -1, O14, Ne18, k_He4_O14_to_Ne18_reaclib, k_Ne18_to_He4_O14_derived,
        He4, O15, -1, Ne19, -1, -1, O15, Ne19, k_He4_O15_to_Ne19_reaclib, k_Ne19_to_He4_O15_derived,
        H1, O16, -1, F17, -1, -1, O16, F17, k_p_O16_to_F17_reaclib, k_F17_to_p_O16_derived,
        He4, O16, -1, Ne20, -1, -1, O16, Ne20, k_He4_O16_to_Ne20_reaclib, k_Ne20_to_He4_O16_derived,
        H1, O17, -1, F18, -1, -1, O17, F18, k_p_O17_to_F18_reaclib, k_F18_to_p_O17_derived,
        He4, O17, -1, Ne21, -1, -1, O17, Ne21, k_He4_O17_to_Ne21_reaclib, k_Ne21_to_He4_O17_derived,
        H1, O18, -1, F19, -1, -1, O18, F19, k_p_O18_to_F19_reaclib, k_F19_to_p_O18_derived,
        He4, O18, -1, Ne22, -1, -1, O18, Ne22, k_He4_O18_to_Ne22_reaclib, k_Ne22_to_He4_O18_derived,
        H1, F17, -1, Ne18, -1, -1, F17, Ne18, k_p_F17_to_Ne18_reaclib, k_Ne18_to_p_F17_derived,
        He4, F17, -1, Na21, -1, -1, F17, Na21, k_He4_F17_to_Na21_reaclib, k_Na21_to_He4_F17_derived,
        H1, F18, -1, Ne19, -1, -1, F18, Ne19, k_p_F18_to_Ne19_reaclib, k_Ne19_to_p_F18_derived,
        He4, F18, -1, Na22, -1, -1, F18, Na22, k_He4_F18_to_Na22_reaclib, k_Na22_to_He4_F18_derived,
        H1, F19, -1, Ne20, -1, -1, F19, Ne20, k_p_F19_to_Ne20_reaclib, k_Ne20_to_p_F19_derived,
        He4, F19, -1, Na23, -1, -1, F19, Na23, k_He4_F19_to_Na23_reaclib, k_Na23_to_He4_F19_derived,
        H1, Ne20, -1, Na21, -1, -1, Ne20, Na21, k_p_Ne20_to_Na21_reaclib, k_Na21_to_p_Ne20_derived,
        He4, Ne20, -1, Mg24, -1, -1, Ne20, Mg24, k_He4_Ne20_to_Mg24_reaclib, k_Mg24_to_He4_Ne20_derived,
        H1, Ne21, -1, Na22, -1, -1, Ne21, Na22, k_p_Ne21_to_Na22_reaclib, k_Na22_to_p_Ne21_derived,
        H1, Ne22, -1, Na23, -1, -1, Ne22, Na23, k_p_Ne22_to_Na23_reaclib, k_Na23_to_p_Ne22_derived,
        H1, Na23, -1, Mg24, -1, -1, Na23, Mg24, k_p_Na23_to_Mg24_reaclib, k_Mg24_to_p_Na23_derived,
        He4, Na23, -1, Al27, -1, -1, Na23, Al27, k_He4_Na23_to_Al27_reaclib, k_Al27_to_He4_Na23_derived,
        He4, Mg24, -1, Si28, -1, -1, Mg24, Si28, k_He4_Mg24_to_Si28_reaclib, k_Si28_to_He4_Mg24_derived,
        H1, Al27, -1, Si28, -1, -1, Al27, Si28, k_p_Al27_to_Si28_reaclib, k_Si28_to_p_Al27_derived,
        He4, Al27, -1, P31, -1, -1, Al27, P31, k_He4_Al27_to_P31_reaclib, k_P31_to_He4_Al27_derived,
        He4, Si28, -1, S32, -1, -1, Si28, S32, k_He4_Si28_to_S32_reaclib, k_S32_to_He4_Si28_derived,
        H1, P31, -1, S32, -1, -1, P31, S32, k_p_P31_to_S32_reaclib, k_S32_to_p_P31_derived,
        He4, N13, -1, H1, O16, -1, N13, O16, k_He4_N13_to_p_O16_reaclib, k_p_O16_to_He4_N13_derived,
        H1, N15, -1, He4, C12, -1, N15, C12, k_p_N15_to_He4_C12_reaclib, k_He4_C12_to_p_N15_derived,
        He4, O14, -1, H1, F17, -1, O14, F17, k_He4_O14_to_p_F17_reaclib, k_p_F17_to_He4_O14_derived,
        H1, O17, -1, He4, N14, -1, O17, N14, k_p_O17_to_He4_N14_reaclib, k_He4_N14_to_p_O17_derived,
        H1, O18, -1, He4, N15, -1, O18, N15, k_p_O18_to_He4_N15_reaclib, k_He4_N15_to_p_O18_derived,
        H1, F18, -1, He4, O15, -1, F18, O15, k_p_F18_to_He4_O15_reaclib, k_He4_O15_to_p_F18_derived,
        He4, F18, -1, H1, Ne21, -1, F18, Ne21, k_He4_F18_to_p_Ne21_reaclib, k_p_Ne21_to_He4_F18_derived,
        H1, F19, -1, He4, O16, -1, F19, O16, k_p_F19_to_He4_O16_reaclib, k_He4_O16_to_p_F19_derived,
        He4, F19, -1, H1, Ne22, -1, F19, Ne22, k_He4_F19_to_p_Ne22_reaclib, k_p_Ne22_to_He4_F19_derived,
        He4, Ne18, -1, H1, Na21, -1, Ne18, Na21, k_He4_Ne18_to_p_Na21_reaclib, k_p_Na21_to_He4_Ne18_derived,
        He4, Ne19, -1, H1, Na22, -1, Ne19, Na22, k_He4_Ne19_to_p_Na22_reaclib, k_p_Na22_to_He4_Ne19_derived,
        H1, Na23, -1, He4, Ne20, -1, Na23, Ne20, k_p_Na23_to_He4_Ne20_reaclib, k_He4_Ne20_to_p_Na23_derived,
        H1, Al27, -1, He4, Mg24, -1, Al27, Mg24, k_p_Al27_to_He4_Mg24_reaclib, k_He4_Mg24_to_p_Al27_derived,
        H1, P31, -1, He4, Si28, -1, P31, Si28, k_p_P31_to_He4_Si28_reaclib, k_He4_Si28_to_p_P31_derived,
        He4, He4, He4, C12, -1, -1, C12, -1, k_He4_He4_He4_to_C12_reaclib, k_C12_to_He4_He4_He4_derived,
        He4, F17, -1, H1, Ne20, -1, F17, Ne20, k_He4_F17_to_p_Ne20_derived, k_p_Ne20_to_He4_F17_reaclib,
        He4, Na21, -1, H1, Mg24, -1, Na21, Mg24, k_He4_Na21_to_p_Mg24_derived, k_p_Mg24_to_He4_Na21_reaclib,
        S32, He4, -1, Ar36, -1, -1, S32, Ar36, k_S32_He4_to_Ar36_approx, k_Ar36_to_S32_He4_approx,
        Ar36, He4, -1, Ca40, -1, -1, Ar36, Ca40, k_Ar36_He4_to_Ca40_approx, k_Ca40_to_Ar36_He4_approx,
        Ca40, He4, -1, Ti44, -1, -1, Ca40, Ti44, k_Ca40_He4_to_Ti44_approx, k_Ti44_to_Ca40_He4_approx,
        Ti44, He4, -1, Cr48, -1, -1, Ti44, Cr48, k_Ti44_He4_to_Cr48_approx, k_Cr48_to_Ti44_He4_approx,
        Cr48, He4, -1, Fe52, -1, -1, Cr48, Fe52, k_Cr48_He4_to_Fe52_approx, k_Fe52_to_Cr48_He4_approx,
        Fe52, He4, -1, Ni56, -1, -1, Fe52, Ni56, k_Fe52_He4_to_Ni56_approx, k_Ni56_to_Fe52_He4_approx
    };
}
#endif

#endif
