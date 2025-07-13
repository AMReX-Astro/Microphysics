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
        else if constexpr (spec == Na22) {
            return 174.14457080400098_rt;
        }
        else if constexpr (spec == Na23) {
            return 186.56435240400242_rt;
        }
        else if constexpr (spec == Mg22) {
            return 168.58082376800303_rt;
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
        else if constexpr (spec == Na22) {
            return 3.652262279854593e-23_rt;
        }
        else if constexpr (spec == Na23) {
            return 3.817541002484691e-23_rt;
        }
        else if constexpr (spec == Mg22) {
            return 3.653114641825797e-23_rt;
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
        k_p_F17_to_Ne18 = 19,
        k_p_F18_to_Ne19 = 20,
        k_He4_F18_to_Na22 = 21,
        k_p_F19_to_Ne20 = 22,
        k_He4_F19_to_Na23 = 23,
        k_He4_Ne18_to_Mg22 = 24,
        k_He4_Ne20_to_Mg24 = 25,
        k_p_Ne21_to_Na22 = 26,
        k_p_Na23_to_Mg24 = 27,
        k_He4_Mg24_to_Si28 = 28,
        k_p_Al27_to_Si28 = 29,
        k_He4_Si28_to_S32 = 30,
        k_p_P31_to_S32 = 31,
        k_C12_C12_to_p_Na23 = 32,
        k_C12_C12_to_He4_Ne20 = 33,
        k_He4_N13_to_p_O16 = 34,
        k_p_N15_to_He4_C12 = 35,
        k_He4_O14_to_p_F17 = 36,
        k_C12_O16_to_p_Al27 = 37,
        k_C12_O16_to_He4_Mg24 = 38,
        k_O16_O16_to_p_P31 = 39,
        k_O16_O16_to_He4_Si28 = 40,
        k_p_O17_to_He4_N14 = 41,
        k_p_O18_to_He4_N15 = 42,
        k_p_F18_to_He4_O15 = 43,
        k_He4_F18_to_p_Ne21 = 44,
        k_p_F19_to_He4_O16 = 45,
        k_He4_Ne19_to_p_Na22 = 46,
        k_p_Ne20_to_He4_F17 = 47,
        k_p_Na23_to_He4_Ne20 = 48,
        k_p_Al27_to_He4_Mg24 = 49,
        k_p_P31_to_He4_Si28 = 50,
        k_He4_He4_He4_to_C12 = 51,
        k_C12_C12_to_n_Mg23_removed = 52,
        k_O16_O16_to_n_S31_removed = 53,
        k_C12_O16_to_n_Si27_removed = 54,
        k_He4_S32_to_Ar36_removed = 55,
        k_p_Cl35_to_Ar36_removed = 56,
        k_p_Cl35_to_He4_S32_removed = 57,
        k_He4_Ar36_to_Ca40_removed = 58,
        k_p_K39_to_Ca40_removed = 59,
        k_p_K39_to_He4_Ar36_removed = 60,
        k_He4_Ca40_to_Ti44_removed = 61,
        k_p_Sc43_to_Ti44_removed = 62,
        k_p_Sc43_to_He4_Ca40_removed = 63,
        k_He4_Ti44_to_Cr48_removed = 64,
        k_He4_Ti44_to_p_V47_removed = 65,
        k_p_V47_to_Cr48_removed = 66,
        k_He4_Cr48_to_Fe52_removed = 67,
        k_He4_Cr48_to_p_Mn51_removed = 68,
        k_p_Mn51_to_Fe52_removed = 69,
        k_He4_Fe52_to_Ni56_removed = 70,
        k_He4_Fe52_to_p_Co55_removed = 71,
        k_p_Co55_to_Ni56_removed = 72,
        k_F17_to_O17 = 73,
        k_O17_to_F17 = 74,
        k_F18_to_Ne18 = 75,
        k_F18_to_O18 = 76,
        k_Ne18_to_F18 = 77,
        k_O18_to_F18 = 78,
        k_F19_to_Ne19 = 79,
        k_Ne19_to_F19 = 80,
        k_Mg22_to_Na22 = 81,
        k_Na22_to_Mg22 = 82,
        k_S32_He4_to_Ar36_approx = 83,
        k_Ar36_to_S32_He4_approx = 84,
        k_Ar36_He4_to_Ca40_approx = 85,
        k_Ca40_to_Ar36_He4_approx = 86,
        k_Ca40_He4_to_Ti44_approx = 87,
        k_Ti44_to_Ca40_He4_approx = 88,
        k_Ti44_He4_to_Cr48_approx = 89,
        k_Cr48_to_Ti44_He4_approx = 90,
        k_Cr48_He4_to_Fe52_approx = 91,
        k_Fe52_to_Cr48_He4_approx = 92,
        k_Fe52_He4_to_Ni56_approx = 93,
        k_Ni56_to_Fe52_He4_approx = 94,
        k_C12_C12_to_Mg24_modified = 95,
        k_O16_O16_to_S32_modified = 96,
        k_C12_O16_to_Si28_modified = 97,
        k_N13_to_p_C12_derived = 98,
        k_N14_to_p_C13_derived = 99,
        k_O14_to_p_N13_derived = 100,
        k_O15_to_p_N14_derived = 101,
        k_O16_to_p_N15_derived = 102,
        k_O16_to_He4_C12_derived = 103,
        k_F17_to_p_O16_derived = 104,
        k_F18_to_p_O17_derived = 105,
        k_F18_to_He4_N14_derived = 106,
        k_F19_to_p_O18_derived = 107,
        k_F19_to_He4_N15_derived = 108,
        k_Ne18_to_p_F17_derived = 109,
        k_Ne18_to_He4_O14_derived = 110,
        k_Ne19_to_p_F18_derived = 111,
        k_Ne19_to_He4_O15_derived = 112,
        k_Ne20_to_p_F19_derived = 113,
        k_Ne20_to_He4_O16_derived = 114,
        k_Ne21_to_He4_O17_derived = 115,
        k_Na22_to_p_Ne21_derived = 116,
        k_Na22_to_He4_F18_derived = 117,
        k_Na23_to_He4_F19_derived = 118,
        k_Mg22_to_He4_Ne18_derived = 119,
        k_Mg24_to_p_Na23_derived = 120,
        k_Mg24_to_He4_Ne20_derived = 121,
        k_Si28_to_p_Al27_derived = 122,
        k_Si28_to_He4_Mg24_derived = 123,
        k_S32_to_p_P31_derived = 124,
        k_S32_to_He4_Si28_derived = 125,
        k_C12_to_He4_He4_He4_derived = 126,
        k_He4_C12_to_p_N15_derived = 127,
        k_He4_N14_to_p_O17_derived = 128,
        k_He4_N15_to_p_O18_derived = 129,
        k_He4_O15_to_p_F18_derived = 130,
        k_p_O16_to_He4_N13_derived = 131,
        k_He4_O16_to_p_F19_derived = 132,
        k_p_F17_to_He4_O14_derived = 133,
        k_He4_F17_to_p_Ne20_derived = 134,
        k_He4_Ne20_to_p_Na23_derived = 135,
        k_p_Ne21_to_He4_F18_derived = 136,
        k_p_Na22_to_He4_Ne19_derived = 137,
        k_He4_Mg24_to_p_Al27_derived = 138,
        k_He4_Si28_to_p_P31_derived = 139,
        k_He4_S32_to_p_Cl35_derived_removed = 140,
        k_Ar36_to_He4_S32_derived_removed = 141,
        k_Ar36_to_p_Cl35_derived_removed = 142,
        k_He4_Ar36_to_p_K39_derived_removed = 143,
        k_Ca40_to_He4_Ar36_derived_removed = 144,
        k_Ca40_to_p_K39_derived_removed = 145,
        k_He4_Ca40_to_p_Sc43_derived_removed = 146,
        k_Ti44_to_He4_Ca40_derived_removed = 147,
        k_Ti44_to_p_Sc43_derived_removed = 148,
        k_Cr48_to_He4_Ti44_derived_removed = 149,
        k_Cr48_to_p_V47_derived_removed = 150,
        k_p_V47_to_He4_Ti44_derived_removed = 151,
        k_Fe52_to_He4_Cr48_derived_removed = 152,
        k_Fe52_to_p_Mn51_derived_removed = 153,
        k_p_Mn51_to_He4_Cr48_derived_removed = 154,
        k_Ni56_to_He4_Fe52_derived_removed = 155,
        k_Ni56_to_p_Co55_derived_removed = 156,
        k_p_Co55_to_He4_Fe52_derived_removed = 157,
        NumRates = k_p_Co55_to_He4_Fe52_derived_removed
    };

    // number of reaclib rates

    const int NrateReaclib = 132;

    // number of tabular rates

    const int NrateTabular = 10;

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
        "p_F17_to_Ne18",  // 19,
        "p_F18_to_Ne19",  // 20,
        "He4_F18_to_Na22",  // 21,
        "p_F19_to_Ne20",  // 22,
        "He4_F19_to_Na23",  // 23,
        "He4_Ne18_to_Mg22",  // 24,
        "He4_Ne20_to_Mg24",  // 25,
        "p_Ne21_to_Na22",  // 26,
        "p_Na23_to_Mg24",  // 27,
        "He4_Mg24_to_Si28",  // 28,
        "p_Al27_to_Si28",  // 29,
        "He4_Si28_to_S32",  // 30,
        "p_P31_to_S32",  // 31,
        "C12_C12_to_p_Na23",  // 32,
        "C12_C12_to_He4_Ne20",  // 33,
        "He4_N13_to_p_O16",  // 34,
        "p_N15_to_He4_C12",  // 35,
        "He4_O14_to_p_F17",  // 36,
        "C12_O16_to_p_Al27",  // 37,
        "C12_O16_to_He4_Mg24",  // 38,
        "O16_O16_to_p_P31",  // 39,
        "O16_O16_to_He4_Si28",  // 40,
        "p_O17_to_He4_N14",  // 41,
        "p_O18_to_He4_N15",  // 42,
        "p_F18_to_He4_O15",  // 43,
        "He4_F18_to_p_Ne21",  // 44,
        "p_F19_to_He4_O16",  // 45,
        "He4_Ne19_to_p_Na22",  // 46,
        "p_Ne20_to_He4_F17",  // 47,
        "p_Na23_to_He4_Ne20",  // 48,
        "p_Al27_to_He4_Mg24",  // 49,
        "p_P31_to_He4_Si28",  // 50,
        "He4_He4_He4_to_C12",  // 51,
        "C12_C12_to_n_Mg23_removed",  // 52,
        "O16_O16_to_n_S31_removed",  // 53,
        "C12_O16_to_n_Si27_removed",  // 54,
        "He4_S32_to_Ar36_removed",  // 55,
        "p_Cl35_to_Ar36_removed",  // 56,
        "p_Cl35_to_He4_S32_removed",  // 57,
        "He4_Ar36_to_Ca40_removed",  // 58,
        "p_K39_to_Ca40_removed",  // 59,
        "p_K39_to_He4_Ar36_removed",  // 60,
        "He4_Ca40_to_Ti44_removed",  // 61,
        "p_Sc43_to_Ti44_removed",  // 62,
        "p_Sc43_to_He4_Ca40_removed",  // 63,
        "He4_Ti44_to_Cr48_removed",  // 64,
        "He4_Ti44_to_p_V47_removed",  // 65,
        "p_V47_to_Cr48_removed",  // 66,
        "He4_Cr48_to_Fe52_removed",  // 67,
        "He4_Cr48_to_p_Mn51_removed",  // 68,
        "p_Mn51_to_Fe52_removed",  // 69,
        "He4_Fe52_to_Ni56_removed",  // 70,
        "He4_Fe52_to_p_Co55_removed",  // 71,
        "p_Co55_to_Ni56_removed",  // 72,
        "F17_to_O17",  // 73,
        "O17_to_F17",  // 74,
        "F18_to_Ne18",  // 75,
        "F18_to_O18",  // 76,
        "Ne18_to_F18",  // 77,
        "O18_to_F18",  // 78,
        "F19_to_Ne19",  // 79,
        "Ne19_to_F19",  // 80,
        "Mg22_to_Na22",  // 81,
        "Na22_to_Mg22",  // 82,
        "S32_He4_to_Ar36_approx",  // 83,
        "Ar36_to_S32_He4_approx",  // 84,
        "Ar36_He4_to_Ca40_approx",  // 85,
        "Ca40_to_Ar36_He4_approx",  // 86,
        "Ca40_He4_to_Ti44_approx",  // 87,
        "Ti44_to_Ca40_He4_approx",  // 88,
        "Ti44_He4_to_Cr48_approx",  // 89,
        "Cr48_to_Ti44_He4_approx",  // 90,
        "Cr48_He4_to_Fe52_approx",  // 91,
        "Fe52_to_Cr48_He4_approx",  // 92,
        "Fe52_He4_to_Ni56_approx",  // 93,
        "Ni56_to_Fe52_He4_approx",  // 94,
        "C12_C12_to_Mg24_modified",  // 95,
        "O16_O16_to_S32_modified",  // 96,
        "C12_O16_to_Si28_modified",  // 97,
        "N13_to_p_C12_derived",  // 98,
        "N14_to_p_C13_derived",  // 99,
        "O14_to_p_N13_derived",  // 100,
        "O15_to_p_N14_derived",  // 101,
        "O16_to_p_N15_derived",  // 102,
        "O16_to_He4_C12_derived",  // 103,
        "F17_to_p_O16_derived",  // 104,
        "F18_to_p_O17_derived",  // 105,
        "F18_to_He4_N14_derived",  // 106,
        "F19_to_p_O18_derived",  // 107,
        "F19_to_He4_N15_derived",  // 108,
        "Ne18_to_p_F17_derived",  // 109,
        "Ne18_to_He4_O14_derived",  // 110,
        "Ne19_to_p_F18_derived",  // 111,
        "Ne19_to_He4_O15_derived",  // 112,
        "Ne20_to_p_F19_derived",  // 113,
        "Ne20_to_He4_O16_derived",  // 114,
        "Ne21_to_He4_O17_derived",  // 115,
        "Na22_to_p_Ne21_derived",  // 116,
        "Na22_to_He4_F18_derived",  // 117,
        "Na23_to_He4_F19_derived",  // 118,
        "Mg22_to_He4_Ne18_derived",  // 119,
        "Mg24_to_p_Na23_derived",  // 120,
        "Mg24_to_He4_Ne20_derived",  // 121,
        "Si28_to_p_Al27_derived",  // 122,
        "Si28_to_He4_Mg24_derived",  // 123,
        "S32_to_p_P31_derived",  // 124,
        "S32_to_He4_Si28_derived",  // 125,
        "C12_to_He4_He4_He4_derived",  // 126,
        "He4_C12_to_p_N15_derived",  // 127,
        "He4_N14_to_p_O17_derived",  // 128,
        "He4_N15_to_p_O18_derived",  // 129,
        "He4_O15_to_p_F18_derived",  // 130,
        "p_O16_to_He4_N13_derived",  // 131,
        "He4_O16_to_p_F19_derived",  // 132,
        "p_F17_to_He4_O14_derived",  // 133,
        "He4_F17_to_p_Ne20_derived",  // 134,
        "He4_Ne20_to_p_Na23_derived",  // 135,
        "p_Ne21_to_He4_F18_derived",  // 136,
        "p_Na22_to_He4_Ne19_derived",  // 137,
        "He4_Mg24_to_p_Al27_derived",  // 138,
        "He4_Si28_to_p_P31_derived",  // 139,
        "He4_S32_to_p_Cl35_derived_removed",  // 140,
        "Ar36_to_He4_S32_derived_removed",  // 141,
        "Ar36_to_p_Cl35_derived_removed",  // 142,
        "He4_Ar36_to_p_K39_derived_removed",  // 143,
        "Ca40_to_He4_Ar36_derived_removed",  // 144,
        "Ca40_to_p_K39_derived_removed",  // 145,
        "He4_Ca40_to_p_Sc43_derived_removed",  // 146,
        "Ti44_to_He4_Ca40_derived_removed",  // 147,
        "Ti44_to_p_Sc43_derived_removed",  // 148,
        "Cr48_to_He4_Ti44_derived_removed",  // 149,
        "Cr48_to_p_V47_derived_removed",  // 150,
        "p_V47_to_He4_Ti44_derived_removed",  // 151,
        "Fe52_to_He4_Cr48_derived_removed",  // 152,
        "Fe52_to_p_Mn51_derived_removed",  // 153,
        "p_Mn51_to_He4_Cr48_derived_removed",  // 154,
        "Ni56_to_He4_Fe52_derived_removed",  // 155,
        "Ni56_to_p_Co55_derived_removed",  // 156,
        "p_Co55_to_He4_Fe52_derived_removed"  // 157,
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
