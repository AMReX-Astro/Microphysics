#ifndef actual_network_H
#define actual_network_H

#include <string_view>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Loop.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex;

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
        k_N13_to_C13_weak_wc12 = 1,
        k_O14_to_N14_weak_wc12 = 2,
        k_O15_to_N15_weak_wc12 = 3,
        k_N13_to_p_C12 = 4,
        k_N14_to_p_C13 = 5,
        k_O14_to_p_N13 = 6,
        k_O15_to_p_N14 = 7,
        k_O16_to_p_N15 = 8,
        k_O16_to_He4_C12 = 9,
        k_F17_to_p_O16 = 10,
        k_F18_to_p_O17 = 11,
        k_F18_to_He4_N14 = 12,
        k_F19_to_p_O18 = 13,
        k_F19_to_He4_N15 = 14,
        k_Ne18_to_p_F17 = 15,
        k_Ne18_to_He4_O14 = 16,
        k_Ne19_to_p_F18 = 17,
        k_Ne19_to_He4_O15 = 18,
        k_Ne20_to_p_F19 = 19,
        k_Ne20_to_He4_O16 = 20,
        k_Ne21_to_He4_O17 = 21,
        k_Na22_to_p_Ne21 = 22,
        k_Na22_to_He4_F18 = 23,
        k_Na23_to_He4_F19 = 24,
        k_Mg22_to_He4_Ne18 = 25,
        k_Mg24_to_p_Na23 = 26,
        k_Mg24_to_He4_Ne20 = 27,
        k_Si28_to_p_Al27 = 28,
        k_Si28_to_He4_Mg24 = 29,
        k_S32_to_p_P31 = 30,
        k_S32_to_He4_Si28 = 31,
        k_C12_to_He4_He4_He4 = 32,
        k_p_C12_to_N13 = 33,
        k_He4_C12_to_O16 = 34,
        k_p_C13_to_N14 = 35,
        k_p_N13_to_O14 = 36,
        k_p_N14_to_O15 = 37,
        k_He4_N14_to_F18 = 38,
        k_p_N15_to_O16 = 39,
        k_He4_N15_to_F19 = 40,
        k_He4_O14_to_Ne18 = 41,
        k_He4_O15_to_Ne19 = 42,
        k_p_O16_to_F17 = 43,
        k_He4_O16_to_Ne20 = 44,
        k_p_O17_to_F18 = 45,
        k_He4_O17_to_Ne21 = 46,
        k_p_O18_to_F19 = 47,
        k_p_F17_to_Ne18 = 48,
        k_p_F18_to_Ne19 = 49,
        k_He4_F18_to_Na22 = 50,
        k_p_F19_to_Ne20 = 51,
        k_He4_F19_to_Na23 = 52,
        k_He4_Ne18_to_Mg22 = 53,
        k_He4_Ne20_to_Mg24 = 54,
        k_p_Ne21_to_Na22 = 55,
        k_p_Na23_to_Mg24 = 56,
        k_He4_Mg24_to_Si28 = 57,
        k_p_Al27_to_Si28 = 58,
        k_He4_Si28_to_S32 = 59,
        k_p_P31_to_S32 = 60,
        k_He4_C12_to_p_N15 = 61,
        k_C12_C12_to_p_Na23 = 62,
        k_C12_C12_to_He4_Ne20 = 63,
        k_He4_N13_to_p_O16 = 64,
        k_He4_N14_to_p_O17 = 65,
        k_p_N15_to_He4_C12 = 66,
        k_He4_N15_to_p_O18 = 67,
        k_He4_O14_to_p_F17 = 68,
        k_He4_O15_to_p_F18 = 69,
        k_p_O16_to_He4_N13 = 70,
        k_He4_O16_to_p_F19 = 71,
        k_C12_O16_to_p_Al27 = 72,
        k_C12_O16_to_He4_Mg24 = 73,
        k_O16_O16_to_p_P31 = 74,
        k_O16_O16_to_He4_Si28 = 75,
        k_p_O17_to_He4_N14 = 76,
        k_p_O18_to_He4_N15 = 77,
        k_p_F17_to_He4_O14 = 78,
        k_He4_F17_to_p_Ne20 = 79,
        k_p_F18_to_He4_O15 = 80,
        k_He4_F18_to_p_Ne21 = 81,
        k_p_F19_to_He4_O16 = 82,
        k_He4_Ne19_to_p_Na22 = 83,
        k_p_Ne20_to_He4_F17 = 84,
        k_He4_Ne20_to_p_Na23 = 85,
        k_p_Ne21_to_He4_F18 = 86,
        k_p_Na22_to_He4_Ne19 = 87,
        k_p_Na23_to_He4_Ne20 = 88,
        k_He4_Mg24_to_p_Al27 = 89,
        k_p_Al27_to_He4_Mg24 = 90,
        k_He4_Si28_to_p_P31 = 91,
        k_p_P31_to_He4_Si28 = 92,
        k_He4_He4_He4_to_C12 = 93,
        k_C12_C12_to_Mg24_modified = 94,
        k_O16_O16_to_S32_modified = 95,
        k_C12_O16_to_Si28_modified = 96,
        k_He4_S32_to_Ar36_removed = 97,
        k_He4_S32_to_p_Cl35_removed = 98,
        k_p_Cl35_to_Ar36_removed = 99,
        k_Ar36_to_He4_S32_removed = 100,
        k_Ar36_to_p_Cl35_removed = 101,
        k_p_Cl35_to_He4_S32_removed = 102,
        k_He4_Ar36_to_Ca40_removed = 103,
        k_He4_Ar36_to_p_K39_removed = 104,
        k_p_K39_to_Ca40_removed = 105,
        k_Ca40_to_He4_Ar36_removed = 106,
        k_Ca40_to_p_K39_removed = 107,
        k_p_K39_to_He4_Ar36_removed = 108,
        k_He4_Ca40_to_Ti44_removed = 109,
        k_He4_Ca40_to_p_Sc43_removed = 110,
        k_p_Sc43_to_Ti44_removed = 111,
        k_Ti44_to_He4_Ca40_removed = 112,
        k_Ti44_to_p_Sc43_removed = 113,
        k_p_Sc43_to_He4_Ca40_removed = 114,
        k_He4_Ti44_to_Cr48_removed = 115,
        k_He4_Ti44_to_p_V47_removed = 116,
        k_p_V47_to_Cr48_removed = 117,
        k_Cr48_to_He4_Ti44_removed = 118,
        k_Cr48_to_p_V47_removed = 119,
        k_p_V47_to_He4_Ti44_removed = 120,
        k_He4_Cr48_to_Fe52_removed = 121,
        k_He4_Cr48_to_p_Mn51_removed = 122,
        k_p_Mn51_to_Fe52_removed = 123,
        k_Fe52_to_He4_Cr48_removed = 124,
        k_Fe52_to_p_Mn51_removed = 125,
        k_p_Mn51_to_He4_Cr48_removed = 126,
        k_He4_Fe52_to_Ni56_removed = 127,
        k_He4_Fe52_to_p_Co55_removed = 128,
        k_p_Co55_to_Ni56_removed = 129,
        k_Ni56_to_He4_Fe52_removed = 130,
        k_Ni56_to_p_Co55_removed = 131,
        k_p_Co55_to_He4_Fe52_removed = 132,
        k_F17_to_O17 = 133,
        k_O17_to_F17 = 134,
        k_F18_to_Ne18 = 135,
        k_F18_to_O18 = 136,
        k_Ne18_to_F18 = 137,
        k_O18_to_F18 = 138,
        k_F19_to_Ne19 = 139,
        k_Ne19_to_F19 = 140,
        k_Mg22_to_Na22 = 141,
        k_Na22_to_Mg22 = 142,
        k_S32_He4_to_Ar36_approx = 143,
        k_Ar36_to_S32_He4_approx = 144,
        k_Ar36_He4_to_Ca40_approx = 145,
        k_Ca40_to_Ar36_He4_approx = 146,
        k_Ca40_He4_to_Ti44_approx = 147,
        k_Ti44_to_Ca40_He4_approx = 148,
        k_Ti44_He4_to_Cr48_approx = 149,
        k_Cr48_to_Ti44_He4_approx = 150,
        k_Cr48_He4_to_Fe52_approx = 151,
        k_Fe52_to_Cr48_He4_approx = 152,
        k_Fe52_He4_to_Ni56_approx = 153,
        k_Ni56_to_Fe52_He4_approx = 154,
        NumRates = k_Ni56_to_Fe52_He4_approx
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
        "N13_to_p_C12",  // 4,
        "N14_to_p_C13",  // 5,
        "O14_to_p_N13",  // 6,
        "O15_to_p_N14",  // 7,
        "O16_to_p_N15",  // 8,
        "O16_to_He4_C12",  // 9,
        "F17_to_p_O16",  // 10,
        "F18_to_p_O17",  // 11,
        "F18_to_He4_N14",  // 12,
        "F19_to_p_O18",  // 13,
        "F19_to_He4_N15",  // 14,
        "Ne18_to_p_F17",  // 15,
        "Ne18_to_He4_O14",  // 16,
        "Ne19_to_p_F18",  // 17,
        "Ne19_to_He4_O15",  // 18,
        "Ne20_to_p_F19",  // 19,
        "Ne20_to_He4_O16",  // 20,
        "Ne21_to_He4_O17",  // 21,
        "Na22_to_p_Ne21",  // 22,
        "Na22_to_He4_F18",  // 23,
        "Na23_to_He4_F19",  // 24,
        "Mg22_to_He4_Ne18",  // 25,
        "Mg24_to_p_Na23",  // 26,
        "Mg24_to_He4_Ne20",  // 27,
        "Si28_to_p_Al27",  // 28,
        "Si28_to_He4_Mg24",  // 29,
        "S32_to_p_P31",  // 30,
        "S32_to_He4_Si28",  // 31,
        "C12_to_He4_He4_He4",  // 32,
        "p_C12_to_N13",  // 33,
        "He4_C12_to_O16",  // 34,
        "p_C13_to_N14",  // 35,
        "p_N13_to_O14",  // 36,
        "p_N14_to_O15",  // 37,
        "He4_N14_to_F18",  // 38,
        "p_N15_to_O16",  // 39,
        "He4_N15_to_F19",  // 40,
        "He4_O14_to_Ne18",  // 41,
        "He4_O15_to_Ne19",  // 42,
        "p_O16_to_F17",  // 43,
        "He4_O16_to_Ne20",  // 44,
        "p_O17_to_F18",  // 45,
        "He4_O17_to_Ne21",  // 46,
        "p_O18_to_F19",  // 47,
        "p_F17_to_Ne18",  // 48,
        "p_F18_to_Ne19",  // 49,
        "He4_F18_to_Na22",  // 50,
        "p_F19_to_Ne20",  // 51,
        "He4_F19_to_Na23",  // 52,
        "He4_Ne18_to_Mg22",  // 53,
        "He4_Ne20_to_Mg24",  // 54,
        "p_Ne21_to_Na22",  // 55,
        "p_Na23_to_Mg24",  // 56,
        "He4_Mg24_to_Si28",  // 57,
        "p_Al27_to_Si28",  // 58,
        "He4_Si28_to_S32",  // 59,
        "p_P31_to_S32",  // 60,
        "He4_C12_to_p_N15",  // 61,
        "C12_C12_to_p_Na23",  // 62,
        "C12_C12_to_He4_Ne20",  // 63,
        "He4_N13_to_p_O16",  // 64,
        "He4_N14_to_p_O17",  // 65,
        "p_N15_to_He4_C12",  // 66,
        "He4_N15_to_p_O18",  // 67,
        "He4_O14_to_p_F17",  // 68,
        "He4_O15_to_p_F18",  // 69,
        "p_O16_to_He4_N13",  // 70,
        "He4_O16_to_p_F19",  // 71,
        "C12_O16_to_p_Al27",  // 72,
        "C12_O16_to_He4_Mg24",  // 73,
        "O16_O16_to_p_P31",  // 74,
        "O16_O16_to_He4_Si28",  // 75,
        "p_O17_to_He4_N14",  // 76,
        "p_O18_to_He4_N15",  // 77,
        "p_F17_to_He4_O14",  // 78,
        "He4_F17_to_p_Ne20",  // 79,
        "p_F18_to_He4_O15",  // 80,
        "He4_F18_to_p_Ne21",  // 81,
        "p_F19_to_He4_O16",  // 82,
        "He4_Ne19_to_p_Na22",  // 83,
        "p_Ne20_to_He4_F17",  // 84,
        "He4_Ne20_to_p_Na23",  // 85,
        "p_Ne21_to_He4_F18",  // 86,
        "p_Na22_to_He4_Ne19",  // 87,
        "p_Na23_to_He4_Ne20",  // 88,
        "He4_Mg24_to_p_Al27",  // 89,
        "p_Al27_to_He4_Mg24",  // 90,
        "He4_Si28_to_p_P31",  // 91,
        "p_P31_to_He4_Si28",  // 92,
        "He4_He4_He4_to_C12",  // 93,
        "C12_C12_to_Mg24_modified",  // 94,
        "O16_O16_to_S32_modified",  // 95,
        "C12_O16_to_Si28_modified",  // 96,
        "He4_S32_to_Ar36_removed",  // 97,
        "He4_S32_to_p_Cl35_removed",  // 98,
        "p_Cl35_to_Ar36_removed",  // 99,
        "Ar36_to_He4_S32_removed",  // 100,
        "Ar36_to_p_Cl35_removed",  // 101,
        "p_Cl35_to_He4_S32_removed",  // 102,
        "He4_Ar36_to_Ca40_removed",  // 103,
        "He4_Ar36_to_p_K39_removed",  // 104,
        "p_K39_to_Ca40_removed",  // 105,
        "Ca40_to_He4_Ar36_removed",  // 106,
        "Ca40_to_p_K39_removed",  // 107,
        "p_K39_to_He4_Ar36_removed",  // 108,
        "He4_Ca40_to_Ti44_removed",  // 109,
        "He4_Ca40_to_p_Sc43_removed",  // 110,
        "p_Sc43_to_Ti44_removed",  // 111,
        "Ti44_to_He4_Ca40_removed",  // 112,
        "Ti44_to_p_Sc43_removed",  // 113,
        "p_Sc43_to_He4_Ca40_removed",  // 114,
        "He4_Ti44_to_Cr48_removed",  // 115,
        "He4_Ti44_to_p_V47_removed",  // 116,
        "p_V47_to_Cr48_removed",  // 117,
        "Cr48_to_He4_Ti44_removed",  // 118,
        "Cr48_to_p_V47_removed",  // 119,
        "p_V47_to_He4_Ti44_removed",  // 120,
        "He4_Cr48_to_Fe52_removed",  // 121,
        "He4_Cr48_to_p_Mn51_removed",  // 122,
        "p_Mn51_to_Fe52_removed",  // 123,
        "Fe52_to_He4_Cr48_removed",  // 124,
        "Fe52_to_p_Mn51_removed",  // 125,
        "p_Mn51_to_He4_Cr48_removed",  // 126,
        "He4_Fe52_to_Ni56_removed",  // 127,
        "He4_Fe52_to_p_Co55_removed",  // 128,
        "p_Co55_to_Ni56_removed",  // 129,
        "Ni56_to_He4_Fe52_removed",  // 130,
        "Ni56_to_p_Co55_removed",  // 131,
        "p_Co55_to_He4_Fe52_removed",  // 132,
        "F17_to_O17",  // 133,
        "O17_to_F17",  // 134,
        "F18_to_Ne18",  // 135,
        "F18_to_O18",  // 136,
        "Ne18_to_F18",  // 137,
        "O18_to_F18",  // 138,
        "F19_to_Ne19",  // 139,
        "Ne19_to_F19",  // 140,
        "Mg22_to_Na22",  // 141,
        "Na22_to_Mg22",  // 142,
        "S32_He4_to_Ar36_approx",  // 143,
        "Ar36_to_S32_He4_approx",  // 144,
        "Ar36_He4_to_Ca40_approx",  // 145,
        "Ca40_to_Ar36_He4_approx",  // 146,
        "Ca40_He4_to_Ti44_approx",  // 147,
        "Ti44_to_Ca40_He4_approx",  // 148,
        "Ti44_He4_to_Cr48_approx",  // 149,
        "Cr48_to_Ti44_He4_approx",  // 150,
        "Cr48_He4_to_Fe52_approx",  // 151,
        "Fe52_to_Cr48_He4_approx",  // 152,
        "Fe52_He4_to_Ni56_approx",  // 153,
        "Ni56_to_Fe52_He4_approx"  // 154,
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

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
