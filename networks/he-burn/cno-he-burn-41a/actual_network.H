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
        else if constexpr (spec == Mg22) {
            return 168.58082376800303_rt;
        }
        else if constexpr (spec == Mg23) {
            return 181.7258218679999_rt;
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
        else if constexpr (spec == Al25) {
            return 200.52841103199896_rt;
        }
        else if constexpr (spec == Al26) {
            return 211.89389913199557_rt;
        }
        else if constexpr (spec == Al27) {
            return 224.95193723199915_rt;
        }
        else if constexpr (spec == Si27) {
            return 219.3572401960009_rt;
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
        else if constexpr (spec == Mg22) {
            return 3.653114641825797e-23_rt;
        }
        else if constexpr (spec == Mg23) {
            return 3.8182640828719474e-23_rt;
        }
        else if constexpr (spec == Mg24) {
            return 3.9828098739467446e-23_rt;
        }
        else if constexpr (spec == Mg25) {
            return 4.148995839845366e-23_rt;
        }
        else if constexpr (spec == Mg26) {
            return 4.314511068143473e-23_rt;
        }
        else if constexpr (spec == Al25) {
            return 4.149758250478667e-23_rt;
        }
        else if constexpr (spec == Al26) {
            return 4.315224917996432e-23_rt;
        }
        else if constexpr (spec == Al27) {
            return 4.480389861070653e-23_rt;
        }
        else if constexpr (spec == Si27) {
            return 4.4812477403805044e-23_rt;
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
        k_He4_Ne18_to_Mg22 = 26,
        k_He4_Ne19_to_Mg23 = 27,
        k_p_Ne20_to_Na21 = 28,
        k_He4_Ne20_to_Mg24 = 29,
        k_p_Ne21_to_Na22 = 30,
        k_He4_Ne21_to_Mg25 = 31,
        k_p_Ne22_to_Na23 = 32,
        k_He4_Ne22_to_Mg26 = 33,
        k_p_Na21_to_Mg22 = 34,
        k_He4_Na21_to_Al25 = 35,
        k_p_Na22_to_Mg23 = 36,
        k_He4_Na22_to_Al26 = 37,
        k_p_Na23_to_Mg24 = 38,
        k_He4_Mg23_to_Si27 = 39,
        k_p_Mg24_to_Al25 = 40,
        k_He4_Mg24_to_Si28 = 41,
        k_p_Mg25_to_Al26 = 42,
        k_p_Mg26_to_Al27 = 43,
        k_p_Al26_to_Si27 = 44,
        k_p_Al27_to_Si28 = 45,
        k_He4_Si28_to_S32 = 46,
        k_p_P31_to_S32 = 47,
        k_C12_C12_to_p_Na23 = 48,
        k_C12_C12_to_He4_Ne20 = 49,
        k_He4_N13_to_p_O16 = 50,
        k_p_N15_to_He4_C12 = 51,
        k_He4_O14_to_p_F17 = 52,
        k_C12_O16_to_p_Al27 = 53,
        k_C12_O16_to_He4_Mg24 = 54,
        k_O16_O16_to_p_P31 = 55,
        k_O16_O16_to_He4_Si28 = 56,
        k_p_O17_to_He4_N14 = 57,
        k_p_O18_to_He4_N15 = 58,
        k_p_F18_to_He4_O15 = 59,
        k_He4_F18_to_p_Ne21 = 60,
        k_p_F19_to_He4_O16 = 61,
        k_He4_F19_to_p_Ne22 = 62,
        k_He4_Ne18_to_p_Na21 = 63,
        k_He4_Ne19_to_p_Na22 = 64,
        k_p_Ne20_to_He4_F17 = 65,
        k_He4_Na22_to_p_Mg25 = 66,
        k_p_Na23_to_He4_Ne20 = 67,
        k_He4_Na23_to_p_Mg26 = 68,
        k_He4_Mg22_to_p_Al25 = 69,
        k_He4_Mg23_to_p_Al26 = 70,
        k_p_Mg24_to_He4_Na21 = 71,
        k_He4_Al25_to_p_Si28 = 72,
        k_p_Al27_to_He4_Mg24 = 73,
        k_p_P31_to_He4_Si28 = 74,
        k_He4_He4_He4_to_C12 = 75,
        k_C12_C12_to_Mg24_modified = 76,
        k_O16_O16_to_S32_modified = 77,
        k_C12_O16_to_Si28_modified = 78,
        k_He4_S32_to_Ar36_removed = 79,
        k_p_Cl35_to_Ar36_removed = 80,
        k_p_Cl35_to_He4_S32_removed = 81,
        k_He4_Ar36_to_Ca40_removed = 82,
        k_p_K39_to_Ca40_removed = 83,
        k_p_K39_to_He4_Ar36_removed = 84,
        k_He4_Ca40_to_Ti44_removed = 85,
        k_p_Sc43_to_Ti44_removed = 86,
        k_p_Sc43_to_He4_Ca40_removed = 87,
        k_He4_Ti44_to_Cr48_removed = 88,
        k_He4_Ti44_to_p_V47_removed = 89,
        k_p_V47_to_Cr48_removed = 90,
        k_He4_Cr48_to_Fe52_removed = 91,
        k_He4_Cr48_to_p_Mn51_removed = 92,
        k_p_Mn51_to_Fe52_removed = 93,
        k_He4_Fe52_to_Ni56_removed = 94,
        k_He4_Fe52_to_p_Co55_removed = 95,
        k_p_Co55_to_Ni56_removed = 96,
        k_F17_to_O17 = 97,
        k_O17_to_F17 = 98,
        k_F18_to_Ne18 = 99,
        k_F18_to_O18 = 100,
        k_Ne18_to_F18 = 101,
        k_O18_to_F18 = 102,
        k_F19_to_Ne19 = 103,
        k_Ne19_to_F19 = 104,
        k_Na21_to_Ne21 = 105,
        k_Ne21_to_Na21 = 106,
        k_Mg22_to_Na22 = 107,
        k_Na22_to_Mg22 = 108,
        k_Na22_to_Ne22 = 109,
        k_Ne22_to_Na22 = 110,
        k_Mg23_to_Na23 = 111,
        k_Na23_to_Mg23 = 112,
        k_Al25_to_Mg25 = 113,
        k_Mg25_to_Al25 = 114,
        k_Al26_to_Mg26 = 115,
        k_Mg26_to_Al26 = 116,
        k_Al27_to_Si27 = 117,
        k_Si27_to_Al27 = 118,
        k_S32_He4_to_Ar36_approx = 119,
        k_Ar36_to_S32_He4_approx = 120,
        k_Ar36_He4_to_Ca40_approx = 121,
        k_Ca40_to_Ar36_He4_approx = 122,
        k_Ca40_He4_to_Ti44_approx = 123,
        k_Ti44_to_Ca40_He4_approx = 124,
        k_Ti44_He4_to_Cr48_approx = 125,
        k_Cr48_to_Ti44_He4_approx = 126,
        k_Cr48_He4_to_Fe52_approx = 127,
        k_Fe52_to_Cr48_He4_approx = 128,
        k_Fe52_He4_to_Ni56_approx = 129,
        k_Ni56_to_Fe52_He4_approx = 130,
        k_N13_to_p_C12_derived = 131,
        k_N14_to_p_C13_derived = 132,
        k_O14_to_p_N13_derived = 133,
        k_O15_to_p_N14_derived = 134,
        k_O16_to_p_N15_derived = 135,
        k_O16_to_He4_C12_derived = 136,
        k_F17_to_p_O16_derived = 137,
        k_F18_to_p_O17_derived = 138,
        k_F18_to_He4_N14_derived = 139,
        k_F19_to_p_O18_derived = 140,
        k_F19_to_He4_N15_derived = 141,
        k_Ne18_to_p_F17_derived = 142,
        k_Ne18_to_He4_O14_derived = 143,
        k_Ne19_to_p_F18_derived = 144,
        k_Ne19_to_He4_O15_derived = 145,
        k_Ne20_to_p_F19_derived = 146,
        k_Ne20_to_He4_O16_derived = 147,
        k_Ne21_to_He4_O17_derived = 148,
        k_Ne22_to_He4_O18_derived = 149,
        k_Na21_to_p_Ne20_derived = 150,
        k_Na21_to_He4_F17_derived = 151,
        k_Na22_to_p_Ne21_derived = 152,
        k_Na22_to_He4_F18_derived = 153,
        k_Na23_to_p_Ne22_derived = 154,
        k_Na23_to_He4_F19_derived = 155,
        k_Mg22_to_p_Na21_derived = 156,
        k_Mg22_to_He4_Ne18_derived = 157,
        k_Mg23_to_p_Na22_derived = 158,
        k_Mg23_to_He4_Ne19_derived = 159,
        k_Mg24_to_p_Na23_derived = 160,
        k_Mg24_to_He4_Ne20_derived = 161,
        k_Mg25_to_He4_Ne21_derived = 162,
        k_Mg26_to_He4_Ne22_derived = 163,
        k_Al25_to_p_Mg24_derived = 164,
        k_Al25_to_He4_Na21_derived = 165,
        k_Al26_to_p_Mg25_derived = 166,
        k_Al26_to_He4_Na22_derived = 167,
        k_Al27_to_p_Mg26_derived = 168,
        k_Si27_to_p_Al26_derived = 169,
        k_Si27_to_He4_Mg23_derived = 170,
        k_Si28_to_p_Al27_derived = 171,
        k_Si28_to_He4_Mg24_derived = 172,
        k_S32_to_p_P31_derived = 173,
        k_S32_to_He4_Si28_derived = 174,
        k_C12_to_He4_He4_He4_derived = 175,
        k_He4_C12_to_p_N15_derived = 176,
        k_He4_N14_to_p_O17_derived = 177,
        k_He4_N15_to_p_O18_derived = 178,
        k_He4_O15_to_p_F18_derived = 179,
        k_p_O16_to_He4_N13_derived = 180,
        k_He4_O16_to_p_F19_derived = 181,
        k_p_F17_to_He4_O14_derived = 182,
        k_He4_F17_to_p_Ne20_derived = 183,
        k_He4_Ne20_to_p_Na23_derived = 184,
        k_p_Ne21_to_He4_F18_derived = 185,
        k_p_Ne22_to_He4_F19_derived = 186,
        k_p_Na21_to_He4_Ne18_derived = 187,
        k_He4_Na21_to_p_Mg24_derived = 188,
        k_p_Na22_to_He4_Ne19_derived = 189,
        k_He4_Mg24_to_p_Al27_derived = 190,
        k_p_Mg25_to_He4_Na22_derived = 191,
        k_p_Mg26_to_He4_Na23_derived = 192,
        k_p_Al25_to_He4_Mg22_derived = 193,
        k_p_Al26_to_He4_Mg23_derived = 194,
        k_p_Si28_to_He4_Al25_derived = 195,
        k_He4_Si28_to_p_P31_derived = 196,
        k_He4_S32_to_p_Cl35_derived_removed = 197,
        k_Ar36_to_He4_S32_derived_removed = 198,
        k_Ar36_to_p_Cl35_derived_removed = 199,
        k_He4_Ar36_to_p_K39_derived_removed = 200,
        k_Ca40_to_He4_Ar36_derived_removed = 201,
        k_Ca40_to_p_K39_derived_removed = 202,
        k_He4_Ca40_to_p_Sc43_derived_removed = 203,
        k_Ti44_to_He4_Ca40_derived_removed = 204,
        k_Ti44_to_p_Sc43_derived_removed = 205,
        k_Cr48_to_He4_Ti44_derived_removed = 206,
        k_Cr48_to_p_V47_derived_removed = 207,
        k_p_V47_to_He4_Ti44_derived_removed = 208,
        k_Fe52_to_He4_Cr48_derived_removed = 209,
        k_Fe52_to_p_Mn51_derived_removed = 210,
        k_p_Mn51_to_He4_Cr48_derived_removed = 211,
        k_Ni56_to_He4_Fe52_derived_removed = 212,
        k_Ni56_to_p_Co55_derived_removed = 213,
        k_p_Co55_to_He4_Fe52_derived_removed = 214,
        NumRates = k_p_Co55_to_He4_Fe52_derived_removed
    };

    // number of reaclib rates

    const int NrateReaclib = 180;

    // number of tabular rates

    const int NrateTabular = 22;

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
        "He4_Ne18_to_Mg22",  // 26,
        "He4_Ne19_to_Mg23",  // 27,
        "p_Ne20_to_Na21",  // 28,
        "He4_Ne20_to_Mg24",  // 29,
        "p_Ne21_to_Na22",  // 30,
        "He4_Ne21_to_Mg25",  // 31,
        "p_Ne22_to_Na23",  // 32,
        "He4_Ne22_to_Mg26",  // 33,
        "p_Na21_to_Mg22",  // 34,
        "He4_Na21_to_Al25",  // 35,
        "p_Na22_to_Mg23",  // 36,
        "He4_Na22_to_Al26",  // 37,
        "p_Na23_to_Mg24",  // 38,
        "He4_Mg23_to_Si27",  // 39,
        "p_Mg24_to_Al25",  // 40,
        "He4_Mg24_to_Si28",  // 41,
        "p_Mg25_to_Al26",  // 42,
        "p_Mg26_to_Al27",  // 43,
        "p_Al26_to_Si27",  // 44,
        "p_Al27_to_Si28",  // 45,
        "He4_Si28_to_S32",  // 46,
        "p_P31_to_S32",  // 47,
        "C12_C12_to_p_Na23",  // 48,
        "C12_C12_to_He4_Ne20",  // 49,
        "He4_N13_to_p_O16",  // 50,
        "p_N15_to_He4_C12",  // 51,
        "He4_O14_to_p_F17",  // 52,
        "C12_O16_to_p_Al27",  // 53,
        "C12_O16_to_He4_Mg24",  // 54,
        "O16_O16_to_p_P31",  // 55,
        "O16_O16_to_He4_Si28",  // 56,
        "p_O17_to_He4_N14",  // 57,
        "p_O18_to_He4_N15",  // 58,
        "p_F18_to_He4_O15",  // 59,
        "He4_F18_to_p_Ne21",  // 60,
        "p_F19_to_He4_O16",  // 61,
        "He4_F19_to_p_Ne22",  // 62,
        "He4_Ne18_to_p_Na21",  // 63,
        "He4_Ne19_to_p_Na22",  // 64,
        "p_Ne20_to_He4_F17",  // 65,
        "He4_Na22_to_p_Mg25",  // 66,
        "p_Na23_to_He4_Ne20",  // 67,
        "He4_Na23_to_p_Mg26",  // 68,
        "He4_Mg22_to_p_Al25",  // 69,
        "He4_Mg23_to_p_Al26",  // 70,
        "p_Mg24_to_He4_Na21",  // 71,
        "He4_Al25_to_p_Si28",  // 72,
        "p_Al27_to_He4_Mg24",  // 73,
        "p_P31_to_He4_Si28",  // 74,
        "He4_He4_He4_to_C12",  // 75,
        "C12_C12_to_Mg24_modified",  // 76,
        "O16_O16_to_S32_modified",  // 77,
        "C12_O16_to_Si28_modified",  // 78,
        "He4_S32_to_Ar36_removed",  // 79,
        "p_Cl35_to_Ar36_removed",  // 80,
        "p_Cl35_to_He4_S32_removed",  // 81,
        "He4_Ar36_to_Ca40_removed",  // 82,
        "p_K39_to_Ca40_removed",  // 83,
        "p_K39_to_He4_Ar36_removed",  // 84,
        "He4_Ca40_to_Ti44_removed",  // 85,
        "p_Sc43_to_Ti44_removed",  // 86,
        "p_Sc43_to_He4_Ca40_removed",  // 87,
        "He4_Ti44_to_Cr48_removed",  // 88,
        "He4_Ti44_to_p_V47_removed",  // 89,
        "p_V47_to_Cr48_removed",  // 90,
        "He4_Cr48_to_Fe52_removed",  // 91,
        "He4_Cr48_to_p_Mn51_removed",  // 92,
        "p_Mn51_to_Fe52_removed",  // 93,
        "He4_Fe52_to_Ni56_removed",  // 94,
        "He4_Fe52_to_p_Co55_removed",  // 95,
        "p_Co55_to_Ni56_removed",  // 96,
        "F17_to_O17",  // 97,
        "O17_to_F17",  // 98,
        "F18_to_Ne18",  // 99,
        "F18_to_O18",  // 100,
        "Ne18_to_F18",  // 101,
        "O18_to_F18",  // 102,
        "F19_to_Ne19",  // 103,
        "Ne19_to_F19",  // 104,
        "Na21_to_Ne21",  // 105,
        "Ne21_to_Na21",  // 106,
        "Mg22_to_Na22",  // 107,
        "Na22_to_Mg22",  // 108,
        "Na22_to_Ne22",  // 109,
        "Ne22_to_Na22",  // 110,
        "Mg23_to_Na23",  // 111,
        "Na23_to_Mg23",  // 112,
        "Al25_to_Mg25",  // 113,
        "Mg25_to_Al25",  // 114,
        "Al26_to_Mg26",  // 115,
        "Mg26_to_Al26",  // 116,
        "Al27_to_Si27",  // 117,
        "Si27_to_Al27",  // 118,
        "S32_He4_to_Ar36_approx",  // 119,
        "Ar36_to_S32_He4_approx",  // 120,
        "Ar36_He4_to_Ca40_approx",  // 121,
        "Ca40_to_Ar36_He4_approx",  // 122,
        "Ca40_He4_to_Ti44_approx",  // 123,
        "Ti44_to_Ca40_He4_approx",  // 124,
        "Ti44_He4_to_Cr48_approx",  // 125,
        "Cr48_to_Ti44_He4_approx",  // 126,
        "Cr48_He4_to_Fe52_approx",  // 127,
        "Fe52_to_Cr48_He4_approx",  // 128,
        "Fe52_He4_to_Ni56_approx",  // 129,
        "Ni56_to_Fe52_He4_approx",  // 130,
        "N13_to_p_C12_derived",  // 131,
        "N14_to_p_C13_derived",  // 132,
        "O14_to_p_N13_derived",  // 133,
        "O15_to_p_N14_derived",  // 134,
        "O16_to_p_N15_derived",  // 135,
        "O16_to_He4_C12_derived",  // 136,
        "F17_to_p_O16_derived",  // 137,
        "F18_to_p_O17_derived",  // 138,
        "F18_to_He4_N14_derived",  // 139,
        "F19_to_p_O18_derived",  // 140,
        "F19_to_He4_N15_derived",  // 141,
        "Ne18_to_p_F17_derived",  // 142,
        "Ne18_to_He4_O14_derived",  // 143,
        "Ne19_to_p_F18_derived",  // 144,
        "Ne19_to_He4_O15_derived",  // 145,
        "Ne20_to_p_F19_derived",  // 146,
        "Ne20_to_He4_O16_derived",  // 147,
        "Ne21_to_He4_O17_derived",  // 148,
        "Ne22_to_He4_O18_derived",  // 149,
        "Na21_to_p_Ne20_derived",  // 150,
        "Na21_to_He4_F17_derived",  // 151,
        "Na22_to_p_Ne21_derived",  // 152,
        "Na22_to_He4_F18_derived",  // 153,
        "Na23_to_p_Ne22_derived",  // 154,
        "Na23_to_He4_F19_derived",  // 155,
        "Mg22_to_p_Na21_derived",  // 156,
        "Mg22_to_He4_Ne18_derived",  // 157,
        "Mg23_to_p_Na22_derived",  // 158,
        "Mg23_to_He4_Ne19_derived",  // 159,
        "Mg24_to_p_Na23_derived",  // 160,
        "Mg24_to_He4_Ne20_derived",  // 161,
        "Mg25_to_He4_Ne21_derived",  // 162,
        "Mg26_to_He4_Ne22_derived",  // 163,
        "Al25_to_p_Mg24_derived",  // 164,
        "Al25_to_He4_Na21_derived",  // 165,
        "Al26_to_p_Mg25_derived",  // 166,
        "Al26_to_He4_Na22_derived",  // 167,
        "Al27_to_p_Mg26_derived",  // 168,
        "Si27_to_p_Al26_derived",  // 169,
        "Si27_to_He4_Mg23_derived",  // 170,
        "Si28_to_p_Al27_derived",  // 171,
        "Si28_to_He4_Mg24_derived",  // 172,
        "S32_to_p_P31_derived",  // 173,
        "S32_to_He4_Si28_derived",  // 174,
        "C12_to_He4_He4_He4_derived",  // 175,
        "He4_C12_to_p_N15_derived",  // 176,
        "He4_N14_to_p_O17_derived",  // 177,
        "He4_N15_to_p_O18_derived",  // 178,
        "He4_O15_to_p_F18_derived",  // 179,
        "p_O16_to_He4_N13_derived",  // 180,
        "He4_O16_to_p_F19_derived",  // 181,
        "p_F17_to_He4_O14_derived",  // 182,
        "He4_F17_to_p_Ne20_derived",  // 183,
        "He4_Ne20_to_p_Na23_derived",  // 184,
        "p_Ne21_to_He4_F18_derived",  // 185,
        "p_Ne22_to_He4_F19_derived",  // 186,
        "p_Na21_to_He4_Ne18_derived",  // 187,
        "He4_Na21_to_p_Mg24_derived",  // 188,
        "p_Na22_to_He4_Ne19_derived",  // 189,
        "He4_Mg24_to_p_Al27_derived",  // 190,
        "p_Mg25_to_He4_Na22_derived",  // 191,
        "p_Mg26_to_He4_Na23_derived",  // 192,
        "p_Al25_to_He4_Mg22_derived",  // 193,
        "p_Al26_to_He4_Mg23_derived",  // 194,
        "p_Si28_to_He4_Al25_derived",  // 195,
        "He4_Si28_to_p_P31_derived",  // 196,
        "He4_S32_to_p_Cl35_derived_removed",  // 197,
        "Ar36_to_He4_S32_derived_removed",  // 198,
        "Ar36_to_p_Cl35_derived_removed",  // 199,
        "He4_Ar36_to_p_K39_derived_removed",  // 200,
        "Ca40_to_He4_Ar36_derived_removed",  // 201,
        "Ca40_to_p_K39_derived_removed",  // 202,
        "He4_Ca40_to_p_Sc43_derived_removed",  // 203,
        "Ti44_to_He4_Ca40_derived_removed",  // 204,
        "Ti44_to_p_Sc43_derived_removed",  // 205,
        "Cr48_to_He4_Ti44_derived_removed",  // 206,
        "Cr48_to_p_V47_derived_removed",  // 207,
        "p_V47_to_He4_Ti44_derived_removed",  // 208,
        "Fe52_to_He4_Cr48_derived_removed",  // 209,
        "Fe52_to_p_Mn51_derived_removed",  // 210,
        "p_Mn51_to_He4_Cr48_derived_removed",  // 211,
        "Ni56_to_He4_Fe52_derived_removed",  // 212,
        "Ni56_to_p_Co55_derived_removed",  // 213,
        "p_Co55_to_He4_Fe52_derived_removed"  // 214,
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
