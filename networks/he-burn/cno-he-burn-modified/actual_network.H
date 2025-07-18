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
        k_N13_to_C13_weak_wc12 = 1,
        k_O14_to_N14_weak_wc12 = 2,
        k_O15_to_N15_weak_wc12 = 3,
        k_p_C12_to_N13 = 4,
        k_He4_C12_to_O16 = 5,
        k_p_C13_to_N14 = 6,
        k_p_N13_to_O14 = 7,
        k_p_N14_to_O15 = 8,
        k_He4_N14_to_F18 = 9,
        k_p_N15_to_O16 = 10,
        k_He4_N15_to_F19 = 11,
        k_He4_O14_to_Ne18 = 12,
        k_He4_O15_to_Ne19 = 13,
        k_p_O16_to_F17 = 14,
        k_He4_O16_to_Ne20 = 15,
        k_p_O17_to_F18 = 16,
        k_He4_O17_to_Ne21 = 17,
        k_p_O18_to_F19 = 18,
        k_He4_O18_to_Ne22 = 19,
        k_p_F17_to_Ne18 = 20,
        k_He4_F17_to_Na21 = 21,
        k_p_F18_to_Ne19 = 22,
        k_He4_F18_to_Na22 = 23,
        k_p_F19_to_Ne20 = 24,
        k_He4_F19_to_Na23 = 25,
        k_p_Ne20_to_Na21 = 26,
        k_He4_Ne20_to_Mg24 = 27,
        k_p_Ne21_to_Na22 = 28,
        k_p_Ne22_to_Na23 = 29,
        k_p_Na23_to_Mg24 = 30,
        k_He4_Mg24_to_Si28 = 31,
        k_p_Al27_to_Si28 = 32,
        k_He4_Si28_to_S32 = 33,
        k_p_P31_to_S32 = 34,
        k_C12_C12_to_p_Na23 = 35,
        k_C12_C12_to_He4_Ne20 = 36,
        k_He4_N13_to_p_O16 = 37,
        k_p_N15_to_He4_C12 = 38,
        k_He4_O14_to_p_F17 = 39,
        k_C12_O16_to_p_Al27 = 40,
        k_C12_O16_to_He4_Mg24 = 41,
        k_O16_O16_to_p_P31 = 42,
        k_O16_O16_to_He4_Si28 = 43,
        k_p_O17_to_He4_N14 = 44,
        k_p_O18_to_He4_N15 = 45,
        k_p_F18_to_He4_O15 = 46,
        k_He4_F18_to_p_Ne21 = 47,
        k_p_F19_to_He4_O16 = 48,
        k_He4_F19_to_p_Ne22 = 49,
        k_He4_Ne18_to_p_Na21 = 50,
        k_He4_Ne19_to_p_Na22 = 51,
        k_p_Ne20_to_He4_F17 = 52,
        k_p_Na23_to_He4_Ne20 = 53,
        k_p_Mg24_to_He4_Na21 = 54,
        k_p_Al27_to_He4_Mg24 = 55,
        k_p_P31_to_He4_Si28 = 56,
        k_He4_He4_He4_to_C12 = 57,
        k_C12_C12_to_Mg24_modified = 58,
        k_O16_O16_to_S32_modified = 59,
        k_C12_O16_to_Si28_modified = 60,
        k_He4_Na22_to_Al26_removed = 61,
        k_p_Mg24_to_Al25_removed = 62,
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
        k_He4_Cr48_to_Fe52_removed = 75,
        k_He4_Cr48_to_p_Mn51_removed = 76,
        k_p_Mn51_to_Fe52_removed = 77,
        k_He4_Fe52_to_Ni56_removed = 78,
        k_He4_Fe52_to_p_Co55_removed = 79,
        k_p_Co55_to_Ni56_removed = 80,
        k_F17_to_O17 = 81,
        k_O17_to_F17 = 82,
        k_F18_to_Ne18 = 83,
        k_F18_to_O18 = 84,
        k_Ne18_to_F18 = 85,
        k_O18_to_F18 = 86,
        k_F19_to_Ne19 = 87,
        k_Ne19_to_F19 = 88,
        k_Na21_to_Ne21 = 89,
        k_Ne21_to_Na21 = 90,
        k_Na22_to_Ne22 = 91,
        k_Ne22_to_Na22 = 92,
        k_S32_He4_to_Ar36_approx = 93,
        k_Ar36_to_S32_He4_approx = 94,
        k_Ar36_He4_to_Ca40_approx = 95,
        k_Ca40_to_Ar36_He4_approx = 96,
        k_Ca40_He4_to_Ti44_approx = 97,
        k_Ti44_to_Ca40_He4_approx = 98,
        k_Ti44_He4_to_Cr48_approx = 99,
        k_Cr48_to_Ti44_He4_approx = 100,
        k_Cr48_He4_to_Fe52_approx = 101,
        k_Fe52_to_Cr48_He4_approx = 102,
        k_Fe52_He4_to_Ni56_approx = 103,
        k_Ni56_to_Fe52_He4_approx = 104,
        k_He4_Na22_to_Si28_modified = 105,
        k_p_Mg24_to_Na22_He4_modified = 106,
        k_N13_to_p_C12_derived = 107,
        k_N14_to_p_C13_derived = 108,
        k_O14_to_p_N13_derived = 109,
        k_O15_to_p_N14_derived = 110,
        k_O16_to_p_N15_derived = 111,
        k_O16_to_He4_C12_derived = 112,
        k_F17_to_p_O16_derived = 113,
        k_F18_to_p_O17_derived = 114,
        k_F18_to_He4_N14_derived = 115,
        k_F19_to_p_O18_derived = 116,
        k_F19_to_He4_N15_derived = 117,
        k_Ne18_to_p_F17_derived = 118,
        k_Ne18_to_He4_O14_derived = 119,
        k_Ne19_to_p_F18_derived = 120,
        k_Ne19_to_He4_O15_derived = 121,
        k_Ne20_to_p_F19_derived = 122,
        k_Ne20_to_He4_O16_derived = 123,
        k_Ne21_to_He4_O17_derived = 124,
        k_Ne22_to_He4_O18_derived = 125,
        k_Na21_to_p_Ne20_derived = 126,
        k_Na21_to_He4_F17_derived = 127,
        k_Na22_to_p_Ne21_derived = 128,
        k_Na22_to_He4_F18_derived = 129,
        k_Na23_to_p_Ne22_derived = 130,
        k_Na23_to_He4_F19_derived = 131,
        k_Mg24_to_p_Na23_derived = 132,
        k_Mg24_to_He4_Ne20_derived = 133,
        k_Si28_to_p_Al27_derived = 134,
        k_Si28_to_He4_Mg24_derived = 135,
        k_S32_to_p_P31_derived = 136,
        k_S32_to_He4_Si28_derived = 137,
        k_C12_to_He4_He4_He4_derived = 138,
        k_He4_C12_to_p_N15_derived = 139,
        k_He4_N14_to_p_O17_derived = 140,
        k_He4_N15_to_p_O18_derived = 141,
        k_He4_O15_to_p_F18_derived = 142,
        k_p_O16_to_He4_N13_derived = 143,
        k_He4_O16_to_p_F19_derived = 144,
        k_p_F17_to_He4_O14_derived = 145,
        k_He4_F17_to_p_Ne20_derived = 146,
        k_He4_Ne20_to_p_Na23_derived = 147,
        k_p_Ne21_to_He4_F18_derived = 148,
        k_p_Ne22_to_He4_F19_derived = 149,
        k_p_Na21_to_He4_Ne18_derived = 150,
        k_He4_Na21_to_p_Mg24_derived = 151,
        k_p_Na22_to_He4_Ne19_derived = 152,
        k_He4_Mg24_to_p_Al27_derived = 153,
        k_He4_Si28_to_p_P31_derived = 154,
        k_He4_S32_to_p_Cl35_derived_removed = 155,
        k_Ar36_to_He4_S32_derived_removed = 156,
        k_Ar36_to_p_Cl35_derived_removed = 157,
        k_He4_Ar36_to_p_K39_derived_removed = 158,
        k_Ca40_to_He4_Ar36_derived_removed = 159,
        k_Ca40_to_p_K39_derived_removed = 160,
        k_He4_Ca40_to_p_Sc43_derived_removed = 161,
        k_Ti44_to_He4_Ca40_derived_removed = 162,
        k_Ti44_to_p_Sc43_derived_removed = 163,
        k_Cr48_to_He4_Ti44_derived_removed = 164,
        k_Cr48_to_p_V47_derived_removed = 165,
        k_p_V47_to_He4_Ti44_derived_removed = 166,
        k_Fe52_to_He4_Cr48_derived_removed = 167,
        k_Fe52_to_p_Mn51_derived_removed = 168,
        k_p_Mn51_to_He4_Cr48_derived_removed = 169,
        k_Ni56_to_He4_Fe52_derived_removed = 170,
        k_Ni56_to_p_Co55_derived_removed = 171,
        k_p_Co55_to_He4_Fe52_derived_removed = 172,
        NumRates = k_p_Co55_to_He4_Fe52_derived_removed
    };

    // number of reaclib rates

    const int NrateReaclib = 146;

    // number of tabular rates

    const int NrateTabular = 12;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "N13_to_C13_weak_wc12",  // 1,
        "O14_to_N14_weak_wc12",  // 2,
        "O15_to_N15_weak_wc12",  // 3,
        "p_C12_to_N13",  // 4,
        "He4_C12_to_O16",  // 5,
        "p_C13_to_N14",  // 6,
        "p_N13_to_O14",  // 7,
        "p_N14_to_O15",  // 8,
        "He4_N14_to_F18",  // 9,
        "p_N15_to_O16",  // 10,
        "He4_N15_to_F19",  // 11,
        "He4_O14_to_Ne18",  // 12,
        "He4_O15_to_Ne19",  // 13,
        "p_O16_to_F17",  // 14,
        "He4_O16_to_Ne20",  // 15,
        "p_O17_to_F18",  // 16,
        "He4_O17_to_Ne21",  // 17,
        "p_O18_to_F19",  // 18,
        "He4_O18_to_Ne22",  // 19,
        "p_F17_to_Ne18",  // 20,
        "He4_F17_to_Na21",  // 21,
        "p_F18_to_Ne19",  // 22,
        "He4_F18_to_Na22",  // 23,
        "p_F19_to_Ne20",  // 24,
        "He4_F19_to_Na23",  // 25,
        "p_Ne20_to_Na21",  // 26,
        "He4_Ne20_to_Mg24",  // 27,
        "p_Ne21_to_Na22",  // 28,
        "p_Ne22_to_Na23",  // 29,
        "p_Na23_to_Mg24",  // 30,
        "He4_Mg24_to_Si28",  // 31,
        "p_Al27_to_Si28",  // 32,
        "He4_Si28_to_S32",  // 33,
        "p_P31_to_S32",  // 34,
        "C12_C12_to_p_Na23",  // 35,
        "C12_C12_to_He4_Ne20",  // 36,
        "He4_N13_to_p_O16",  // 37,
        "p_N15_to_He4_C12",  // 38,
        "He4_O14_to_p_F17",  // 39,
        "C12_O16_to_p_Al27",  // 40,
        "C12_O16_to_He4_Mg24",  // 41,
        "O16_O16_to_p_P31",  // 42,
        "O16_O16_to_He4_Si28",  // 43,
        "p_O17_to_He4_N14",  // 44,
        "p_O18_to_He4_N15",  // 45,
        "p_F18_to_He4_O15",  // 46,
        "He4_F18_to_p_Ne21",  // 47,
        "p_F19_to_He4_O16",  // 48,
        "He4_F19_to_p_Ne22",  // 49,
        "He4_Ne18_to_p_Na21",  // 50,
        "He4_Ne19_to_p_Na22",  // 51,
        "p_Ne20_to_He4_F17",  // 52,
        "p_Na23_to_He4_Ne20",  // 53,
        "p_Mg24_to_He4_Na21",  // 54,
        "p_Al27_to_He4_Mg24",  // 55,
        "p_P31_to_He4_Si28",  // 56,
        "He4_He4_He4_to_C12",  // 57,
        "C12_C12_to_Mg24_modified",  // 58,
        "O16_O16_to_S32_modified",  // 59,
        "C12_O16_to_Si28_modified",  // 60,
        "He4_Na22_to_Al26_removed",  // 61,
        "p_Mg24_to_Al25_removed",  // 62,
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
        "He4_Cr48_to_Fe52_removed",  // 75,
        "He4_Cr48_to_p_Mn51_removed",  // 76,
        "p_Mn51_to_Fe52_removed",  // 77,
        "He4_Fe52_to_Ni56_removed",  // 78,
        "He4_Fe52_to_p_Co55_removed",  // 79,
        "p_Co55_to_Ni56_removed",  // 80,
        "F17_to_O17",  // 81,
        "O17_to_F17",  // 82,
        "F18_to_Ne18",  // 83,
        "F18_to_O18",  // 84,
        "Ne18_to_F18",  // 85,
        "O18_to_F18",  // 86,
        "F19_to_Ne19",  // 87,
        "Ne19_to_F19",  // 88,
        "Na21_to_Ne21",  // 89,
        "Ne21_to_Na21",  // 90,
        "Na22_to_Ne22",  // 91,
        "Ne22_to_Na22",  // 92,
        "S32_He4_to_Ar36_approx",  // 93,
        "Ar36_to_S32_He4_approx",  // 94,
        "Ar36_He4_to_Ca40_approx",  // 95,
        "Ca40_to_Ar36_He4_approx",  // 96,
        "Ca40_He4_to_Ti44_approx",  // 97,
        "Ti44_to_Ca40_He4_approx",  // 98,
        "Ti44_He4_to_Cr48_approx",  // 99,
        "Cr48_to_Ti44_He4_approx",  // 100,
        "Cr48_He4_to_Fe52_approx",  // 101,
        "Fe52_to_Cr48_He4_approx",  // 102,
        "Fe52_He4_to_Ni56_approx",  // 103,
        "Ni56_to_Fe52_He4_approx",  // 104,
        "He4_Na22_to_Si28_modified",  // 105,
        "p_Mg24_to_Na22_He4_modified",  // 106,
        "N13_to_p_C12_derived",  // 107,
        "N14_to_p_C13_derived",  // 108,
        "O14_to_p_N13_derived",  // 109,
        "O15_to_p_N14_derived",  // 110,
        "O16_to_p_N15_derived",  // 111,
        "O16_to_He4_C12_derived",  // 112,
        "F17_to_p_O16_derived",  // 113,
        "F18_to_p_O17_derived",  // 114,
        "F18_to_He4_N14_derived",  // 115,
        "F19_to_p_O18_derived",  // 116,
        "F19_to_He4_N15_derived",  // 117,
        "Ne18_to_p_F17_derived",  // 118,
        "Ne18_to_He4_O14_derived",  // 119,
        "Ne19_to_p_F18_derived",  // 120,
        "Ne19_to_He4_O15_derived",  // 121,
        "Ne20_to_p_F19_derived",  // 122,
        "Ne20_to_He4_O16_derived",  // 123,
        "Ne21_to_He4_O17_derived",  // 124,
        "Ne22_to_He4_O18_derived",  // 125,
        "Na21_to_p_Ne20_derived",  // 126,
        "Na21_to_He4_F17_derived",  // 127,
        "Na22_to_p_Ne21_derived",  // 128,
        "Na22_to_He4_F18_derived",  // 129,
        "Na23_to_p_Ne22_derived",  // 130,
        "Na23_to_He4_F19_derived",  // 131,
        "Mg24_to_p_Na23_derived",  // 132,
        "Mg24_to_He4_Ne20_derived",  // 133,
        "Si28_to_p_Al27_derived",  // 134,
        "Si28_to_He4_Mg24_derived",  // 135,
        "S32_to_p_P31_derived",  // 136,
        "S32_to_He4_Si28_derived",  // 137,
        "C12_to_He4_He4_He4_derived",  // 138,
        "He4_C12_to_p_N15_derived",  // 139,
        "He4_N14_to_p_O17_derived",  // 140,
        "He4_N15_to_p_O18_derived",  // 141,
        "He4_O15_to_p_F18_derived",  // 142,
        "p_O16_to_He4_N13_derived",  // 143,
        "He4_O16_to_p_F19_derived",  // 144,
        "p_F17_to_He4_O14_derived",  // 145,
        "He4_F17_to_p_Ne20_derived",  // 146,
        "He4_Ne20_to_p_Na23_derived",  // 147,
        "p_Ne21_to_He4_F18_derived",  // 148,
        "p_Ne22_to_He4_F19_derived",  // 149,
        "p_Na21_to_He4_Ne18_derived",  // 150,
        "He4_Na21_to_p_Mg24_derived",  // 151,
        "p_Na22_to_He4_Ne19_derived",  // 152,
        "He4_Mg24_to_p_Al27_derived",  // 153,
        "He4_Si28_to_p_P31_derived",  // 154,
        "He4_S32_to_p_Cl35_derived_removed",  // 155,
        "Ar36_to_He4_S32_derived_removed",  // 156,
        "Ar36_to_p_Cl35_derived_removed",  // 157,
        "He4_Ar36_to_p_K39_derived_removed",  // 158,
        "Ca40_to_He4_Ar36_derived_removed",  // 159,
        "Ca40_to_p_K39_derived_removed",  // 160,
        "He4_Ca40_to_p_Sc43_derived_removed",  // 161,
        "Ti44_to_He4_Ca40_derived_removed",  // 162,
        "Ti44_to_p_Sc43_derived_removed",  // 163,
        "Cr48_to_He4_Ti44_derived_removed",  // 164,
        "Cr48_to_p_V47_derived_removed",  // 165,
        "p_V47_to_He4_Ti44_derived_removed",  // 166,
        "Fe52_to_He4_Cr48_derived_removed",  // 167,
        "Fe52_to_p_Mn51_derived_removed",  // 168,
        "p_Mn51_to_He4_Cr48_derived_removed",  // 169,
        "Ni56_to_He4_Fe52_derived_removed",  // 170,
        "Ni56_to_p_Co55_derived_removed",  // 171,
        "p_Co55_to_He4_Fe52_derived_removed"  // 172,
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

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, amrex::Order::C> rate_indices;
}
#endif

#endif
