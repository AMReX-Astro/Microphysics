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
        else if constexpr (spec == Al25) {
            return 200.52841103199896_rt;
        }
        else if constexpr (spec == Al26) {
            return 211.89389913199557_rt;
        }
        else if constexpr (spec == Al27) {
            return 224.95193723199915_rt;
        }
        else if constexpr (spec == Si28) {
            return 236.53684539599638_rt;
        }
        else if constexpr (spec == P29) {
            return 239.28581936000046_rt;
        }
        else if constexpr (spec == P30) {
            return 250.60519745999775_rt;
        }
        else if constexpr (spec == P31) {
            return 262.9161999600037_rt;
        }
        else if constexpr (spec == S32) {
            return 271.78016372399725_rt;
        }
        else if constexpr (spec == Cl33) {
            return 274.05689768799857_rt;
        }
        else if constexpr (spec == Cl34) {
            return 285.5650057879975_rt;
        }
        else if constexpr (spec == Ar36) {
            return 306.7167469519991_rt;
        }
        else if constexpr (spec == K37) {
            return 308.5743760159967_rt;
        }
        else if constexpr (spec == K38) {
            return 320.64625411599985_rt;
        }
        else if constexpr (spec == Ca40) {
            return 342.05218528000114_rt;
        }
        else if constexpr (spec == Sc41) {
            return 343.1371143439974_rt;
        }
        else if constexpr (spec == Sc42) {
            return 354.68707244399405_rt;
        }
        else if constexpr (spec == Ti44) {
            return 375.47496160800074_rt;
        }
        else if constexpr (spec == V45) {
            return 377.1017326720029_rt;
        }
        else if constexpr (spec == V46) {
            return 390.3625507720062_rt;
        }
        else if constexpr (spec == Cr48) {
            return 411.4679399359957_rt;
        }
        else if constexpr (spec == Mn49) {
            return 413.55581099999836_rt;
        }
        else if constexpr (spec == Mn50) {
            return 426.63411910000286_rt;
        }
        else if constexpr (spec == Fe52) {
            return 447.6996182639923_rt;
        }
        else if constexpr (spec == Fe53) {
            return 458.38633636399754_rt;
        }
        else if constexpr (spec == Fe54) {
            return 471.76475446399854_rt;
        }
        else if constexpr (spec == Fe55) {
            return 481.06287256399816_rt;
        }
        else if constexpr (spec == Fe56) {
            return 492.2599506639962_rt;
        }
        else if constexpr (spec == Co54) {
            return 462.7379074280034_rt;
        }
        else if constexpr (spec == Co55) {
            return 476.82912552799826_rt;
        }
        else if constexpr (spec == Co56) {
            return 486.91094362799777_rt;
        }
        else if constexpr (spec == Co57) {
            return 498.28746172798856_rt;
        }
        else if constexpr (spec == Ni56) {
            return 483.9956965919919_rt;
        }
        else if constexpr (spec == Ni57) {
            return 494.2434146919986_rt;
        }
        else if constexpr (spec == Ni58) {
            return 506.4596327920008_rt;
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
        else if constexpr (spec == Al25) {
            return 4.149758250478667e-23_rt;
        }
        else if constexpr (spec == Al26) {
            return 4.315224917996432e-23_rt;
        }
        else if constexpr (spec == Al27) {
            return 4.480389861070653e-23_rt;
        }
        else if constexpr (spec == Si28) {
            return 4.6456779473820677e-23_rt;
        }
        else if constexpr (spec == P29) {
            return 4.812541182037512e-23_rt;
        }
        else if constexpr (spec == P30) {
            return 4.978016069409398e-23_rt;
        }
        else if constexpr (spec == P31) {
            return 5.14331418367544e-23_rt;
        }
        else if constexpr (spec == S32) {
            return 5.309087322384128e-23_rt;
        }
        else if constexpr (spec == Cl33) {
            return 5.476034741466158e-23_rt;
        }
        else if constexpr (spec == Cl34) {
            return 5.641475984659597e-23_rt;
        }
        else if constexpr (spec == Ar36) {
            return 5.972551377884467e-23_rt;
        }
        else if constexpr (spec == K37) {
            return 6.139573509201139e-23_rt;
        }
        else if constexpr (spec == K38) {
            return 6.304914251263421e-23_rt;
        }
        else if constexpr (spec == Ca40) {
            return 6.635944331004904e-23_rt;
        }
        else if constexpr (spec == Sc41) {
            return 6.80310420860826e-23_rt;
        }
        else if constexpr (spec == Sc42) {
            return 6.968537991361557e-23_rt;
        }
        else if constexpr (spec == Ti44) {
            return 7.299678247096977e-23_rt;
        }
        else if constexpr (spec == V45) {
            return 7.466741532590238e-23_rt;
        }
        else if constexpr (spec == V46) {
            return 7.631870326846012e-23_rt;
        }
        else if constexpr (spec == Cr48) {
            return 7.962953983065421e-23_rt;
        }
        else if constexpr (spec == Mn49) {
            return 8.129935070017476e-23_rt;
        }
        else if constexpr (spec == Mn50) {
            return 8.295096399635981e-23_rt;
        }
        else if constexpr (spec == Fe52) {
            return 8.626187166893794e-23_rt;
        }
        else if constexpr (spec == Fe53) {
            return 8.791774836154814e-23_rt;
        }
        else if constexpr (spec == Fe54) {
            return 8.95688266630639e-23_rt;
        }
        else if constexpr (spec == Fe55) {
            return 9.122717876001846e-23_rt;
        }
        else if constexpr (spec == Fe56) {
            return 9.288214565329033e-23_rt;
        }
        else if constexpr (spec == Co54) {
            return 8.958352381927671e-23_rt;
        }
        else if constexpr (spec == Co55) {
            return 9.123333143937473e-23_rt;
        }
        else if constexpr (spec == Co56) {
            return 9.289028646418131e-23_rt;
        }
        else if constexpr (spec == Co57) {
            return 9.454493347659797e-23_rt;
        }
        else if constexpr (spec == Ni56) {
            return 9.289408870379396e-23_rt;
        }
        else if constexpr (spec == Ni57) {
            return 9.455074798498772e-23_rt;
        }
        else if constexpr (spec == Ni58) {
            return 9.620389809618879e-23_rt;
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
    std::uint16_t
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
        k_He4_Na21_to_Al25 = 30,
        k_He4_Na22_to_Al26 = 31,
        k_p_Na23_to_Mg24 = 32,
        k_p_Mg24_to_Al25 = 33,
        k_He4_Mg24_to_Si28 = 34,
        k_He4_Al25_to_P29 = 35,
        k_He4_Al26_to_P30 = 36,
        k_p_Al27_to_Si28 = 37,
        k_p_Si28_to_P29 = 38,
        k_He4_Si28_to_S32 = 39,
        k_He4_P29_to_Cl33 = 40,
        k_He4_P30_to_Cl34 = 41,
        k_p_P31_to_S32 = 42,
        k_p_S32_to_Cl33 = 43,
        k_He4_Cl33_to_K37 = 44,
        k_He4_Cl34_to_K38 = 45,
        k_p_Ar36_to_K37 = 46,
        k_He4_K37_to_Sc41 = 47,
        k_He4_K38_to_Sc42 = 48,
        k_p_Ca40_to_Sc41 = 49,
        k_He4_Sc41_to_V45 = 50,
        k_He4_Sc42_to_V46 = 51,
        k_p_Ti44_to_V45 = 52,
        k_He4_V45_to_Mn49 = 53,
        k_He4_V46_to_Mn50 = 54,
        k_p_Cr48_to_Mn49 = 55,
        k_He4_Mn50_to_Co54 = 56,
        k_He4_Fe52_to_Ni56 = 57,
        k_p_Co55_to_Ni56 = 58,
        k_C12_C12_to_p_Na23 = 59,
        k_C12_C12_to_He4_Ne20 = 60,
        k_He4_N13_to_p_O16 = 61,
        k_p_N15_to_He4_C12 = 62,
        k_He4_O14_to_p_F17 = 63,
        k_C12_O16_to_p_Al27 = 64,
        k_C12_O16_to_He4_Mg24 = 65,
        k_O16_O16_to_p_P31 = 66,
        k_O16_O16_to_He4_Si28 = 67,
        k_p_O17_to_He4_N14 = 68,
        k_p_O18_to_He4_N15 = 69,
        k_p_F18_to_He4_O15 = 70,
        k_He4_F18_to_p_Ne21 = 71,
        k_p_F19_to_He4_O16 = 72,
        k_He4_F19_to_p_Ne22 = 73,
        k_He4_Ne18_to_p_Na21 = 74,
        k_He4_Ne19_to_p_Na22 = 75,
        k_p_Ne20_to_He4_F17 = 76,
        k_p_Na23_to_He4_Ne20 = 77,
        k_p_Mg24_to_He4_Na21 = 78,
        k_He4_Al25_to_p_Si28 = 79,
        k_p_Al27_to_He4_Mg24 = 80,
        k_He4_P29_to_p_S32 = 81,
        k_p_P31_to_He4_Si28 = 82,
        k_He4_Cl33_to_p_Ar36 = 83,
        k_He4_K37_to_p_Ca40 = 84,
        k_He4_Sc41_to_p_Ti44 = 85,
        k_He4_V45_to_p_Cr48 = 86,
        k_He4_Mn49_to_p_Fe52 = 87,
        k_He4_Fe52_to_p_Co55 = 88,
        k_He4_He4_He4_to_C12 = 89,
        k_C12_C12_to_Mg24_modified = 90,
        k_O16_O16_to_S32_modified = 91,
        k_C12_O16_to_Si28_modified = 92,
        k_n_Fe52_to_Fe53 = 93,
        k_n_Fe53_to_Fe54 = 94,
        k_He4_Fe53_to_Ni57 = 95,
        k_n_Fe54_to_Fe55 = 96,
        k_p_Fe54_to_Co55 = 97,
        k_He4_Fe54_to_Ni58 = 98,
        k_n_Fe55_to_Fe56 = 99,
        k_p_Fe55_to_Co56 = 100,
        k_p_Fe56_to_Co57 = 101,
        k_n_Co55_to_Co56 = 102,
        k_n_Co56_to_Co57 = 103,
        k_p_Co56_to_Ni57 = 104,
        k_p_Co57_to_Ni58 = 105,
        k_n_Ni56_to_Ni57 = 106,
        k_n_Ni57_to_Ni58 = 107,
        k_He4_Fe53_to_p_Co56 = 108,
        k_n_Co55_to_p_Fe55 = 109,
        k_He4_Co55_to_p_Ni58 = 110,
        k_n_Co56_to_p_Fe56 = 111,
        k_p_Co57_to_He4_Fe54 = 112,
        k_n_Ni56_to_p_Co56 = 113,
        k_n_Ni56_to_He4_Fe53 = 114,
        k_n_Ni57_to_p_Co57 = 115,
        k_n_Ni57_to_He4_Fe54 = 116,
        k_n_Ni58_to_He4_Fe55 = 117,
        k_He4_S32_to_Ar36_removed = 118,
        k_p_Cl35_to_Ar36_removed = 119,
        k_p_Cl35_to_He4_S32_removed = 120,
        k_He4_Ar36_to_Ca40_removed = 121,
        k_p_K39_to_Ca40_removed = 122,
        k_p_K39_to_He4_Ar36_removed = 123,
        k_He4_Ca40_to_Ti44_removed = 124,
        k_p_Sc43_to_Ti44_removed = 125,
        k_p_Sc43_to_He4_Ca40_removed = 126,
        k_He4_Ti44_to_Cr48_removed = 127,
        k_He4_Ti44_to_p_V47_removed = 128,
        k_p_V47_to_Cr48_removed = 129,
        k_He4_Cr48_to_Fe52_removed = 130,
        k_He4_Cr48_to_p_Mn51_removed = 131,
        k_p_Mn51_to_Fe52_removed = 132,
        k_Co54_to_Fe54 = 133,
        k_Fe54_to_Co54 = 134,
        k_Co55_to_Fe55 = 135,
        k_Fe55_to_Co55 = 136,
        k_Co56_to_Fe56 = 137,
        k_Co56_to_Ni56 = 138,
        k_Fe56_to_Co56 = 139,
        k_Ni56_to_Co56 = 140,
        k_Co57_to_Ni57 = 141,
        k_Ni57_to_Co57 = 142,
        k_n_to_p = 143,
        k_p_to_n = 144,
        k_F17_to_O17 = 145,
        k_O17_to_F17 = 146,
        k_F18_to_Ne18 = 147,
        k_F18_to_O18 = 148,
        k_Ne18_to_F18 = 149,
        k_O18_to_F18 = 150,
        k_F19_to_Ne19 = 151,
        k_Ne19_to_F19 = 152,
        k_Na21_to_Ne21 = 153,
        k_Ne21_to_Na21 = 154,
        k_Na22_to_Ne22 = 155,
        k_Ne22_to_Na22 = 156,
        k_S32_He4_to_Ar36_approx = 157,
        k_Ar36_to_S32_He4_approx = 158,
        k_Ar36_He4_to_Ca40_approx = 159,
        k_Ca40_to_Ar36_He4_approx = 160,
        k_Ca40_He4_to_Ti44_approx = 161,
        k_Ti44_to_Ca40_He4_approx = 162,
        k_Ti44_He4_to_Cr48_approx = 163,
        k_Cr48_to_Ti44_He4_approx = 164,
        k_Cr48_He4_to_Fe52_approx = 165,
        k_Fe52_to_Cr48_He4_approx = 166,
        k_N13_to_p_C12_derived = 167,
        k_N14_to_p_C13_derived = 168,
        k_O14_to_p_N13_derived = 169,
        k_O15_to_p_N14_derived = 170,
        k_O16_to_p_N15_derived = 171,
        k_O16_to_He4_C12_derived = 172,
        k_F17_to_p_O16_derived = 173,
        k_F18_to_p_O17_derived = 174,
        k_F18_to_He4_N14_derived = 175,
        k_F19_to_p_O18_derived = 176,
        k_F19_to_He4_N15_derived = 177,
        k_Ne18_to_p_F17_derived = 178,
        k_Ne18_to_He4_O14_derived = 179,
        k_Ne19_to_p_F18_derived = 180,
        k_Ne19_to_He4_O15_derived = 181,
        k_Ne20_to_p_F19_derived = 182,
        k_Ne20_to_He4_O16_derived = 183,
        k_Ne21_to_He4_O17_derived = 184,
        k_Ne22_to_He4_O18_derived = 185,
        k_Na21_to_p_Ne20_derived = 186,
        k_Na21_to_He4_F17_derived = 187,
        k_Na22_to_p_Ne21_derived = 188,
        k_Na22_to_He4_F18_derived = 189,
        k_Na23_to_p_Ne22_derived = 190,
        k_Na23_to_He4_F19_derived = 191,
        k_Mg24_to_p_Na23_derived = 192,
        k_Mg24_to_He4_Ne20_derived = 193,
        k_Al25_to_p_Mg24_derived = 194,
        k_Al25_to_He4_Na21_derived = 195,
        k_Al26_to_He4_Na22_derived = 196,
        k_Si28_to_p_Al27_derived = 197,
        k_Si28_to_He4_Mg24_derived = 198,
        k_P29_to_p_Si28_derived = 199,
        k_P29_to_He4_Al25_derived = 200,
        k_P30_to_He4_Al26_derived = 201,
        k_S32_to_p_P31_derived = 202,
        k_S32_to_He4_Si28_derived = 203,
        k_Cl33_to_p_S32_derived = 204,
        k_Cl33_to_He4_P29_derived = 205,
        k_Cl34_to_He4_P30_derived = 206,
        k_K37_to_p_Ar36_derived = 207,
        k_K37_to_He4_Cl33_derived = 208,
        k_K38_to_He4_Cl34_derived = 209,
        k_Sc41_to_p_Ca40_derived = 210,
        k_Sc41_to_He4_K37_derived = 211,
        k_Sc42_to_He4_K38_derived = 212,
        k_V45_to_p_Ti44_derived = 213,
        k_V45_to_He4_Sc41_derived = 214,
        k_V46_to_He4_Sc42_derived = 215,
        k_Mn49_to_p_Cr48_derived = 216,
        k_Mn49_to_He4_V45_derived = 217,
        k_Mn50_to_He4_V46_derived = 218,
        k_Co54_to_He4_Mn50_derived = 219,
        k_Ni56_to_p_Co55_derived = 220,
        k_Ni56_to_He4_Fe52_derived = 221,
        k_C12_to_He4_He4_He4_derived = 222,
        k_He4_C12_to_p_N15_derived = 223,
        k_He4_N14_to_p_O17_derived = 224,
        k_He4_N15_to_p_O18_derived = 225,
        k_He4_O15_to_p_F18_derived = 226,
        k_p_O16_to_He4_N13_derived = 227,
        k_He4_O16_to_p_F19_derived = 228,
        k_p_F17_to_He4_O14_derived = 229,
        k_He4_F17_to_p_Ne20_derived = 230,
        k_He4_Ne20_to_p_Na23_derived = 231,
        k_p_Ne21_to_He4_F18_derived = 232,
        k_p_Ne22_to_He4_F19_derived = 233,
        k_p_Na21_to_He4_Ne18_derived = 234,
        k_He4_Na21_to_p_Mg24_derived = 235,
        k_p_Na22_to_He4_Ne19_derived = 236,
        k_He4_Mg24_to_p_Al27_derived = 237,
        k_p_Si28_to_He4_Al25_derived = 238,
        k_He4_Si28_to_p_P31_derived = 239,
        k_p_S32_to_He4_P29_derived = 240,
        k_p_Ar36_to_He4_Cl33_derived = 241,
        k_p_Ca40_to_He4_K37_derived = 242,
        k_p_Ti44_to_He4_Sc41_derived = 243,
        k_p_Cr48_to_He4_V45_derived = 244,
        k_p_Fe52_to_He4_Mn49_derived = 245,
        k_p_Co55_to_He4_Fe52_derived = 246,
        k_Fe53_to_n_Fe52_derived = 247,
        k_Fe54_to_n_Fe53_derived = 248,
        k_Fe55_to_n_Fe54_derived = 249,
        k_Fe56_to_n_Fe55_derived = 250,
        k_Co55_to_p_Fe54_derived = 251,
        k_Co56_to_n_Co55_derived = 252,
        k_Co56_to_p_Fe55_derived = 253,
        k_Co57_to_n_Co56_derived = 254,
        k_Co57_to_p_Fe56_derived = 255,
        k_Ni57_to_n_Ni56_derived = 256,
        k_Ni57_to_p_Co56_derived = 257,
        k_Ni57_to_He4_Fe53_derived = 258,
        k_Ni58_to_n_Ni57_derived = 259,
        k_Ni58_to_p_Co57_derived = 260,
        k_Ni58_to_He4_Fe54_derived = 261,
        k_He4_Fe53_to_n_Ni56_derived = 262,
        k_He4_Fe54_to_n_Ni57_derived = 263,
        k_He4_Fe54_to_p_Co57_derived = 264,
        k_p_Fe55_to_n_Co55_derived = 265,
        k_He4_Fe55_to_n_Ni58_derived = 266,
        k_p_Fe56_to_n_Co56_derived = 267,
        k_p_Co56_to_n_Ni56_derived = 268,
        k_p_Co56_to_He4_Fe53_derived = 269,
        k_p_Co57_to_n_Ni57_derived = 270,
        k_p_Ni58_to_He4_Co55_derived = 271,
        k_He4_S32_to_p_Cl35_derived_removed = 272,
        k_Ar36_to_He4_S32_derived_removed = 273,
        k_Ar36_to_p_Cl35_derived_removed = 274,
        k_He4_Ar36_to_p_K39_derived_removed = 275,
        k_Ca40_to_He4_Ar36_derived_removed = 276,
        k_Ca40_to_p_K39_derived_removed = 277,
        k_He4_Ca40_to_p_Sc43_derived_removed = 278,
        k_Ti44_to_He4_Ca40_derived_removed = 279,
        k_Ti44_to_p_Sc43_derived_removed = 280,
        k_Cr48_to_He4_Ti44_derived_removed = 281,
        k_Cr48_to_p_V47_derived_removed = 282,
        k_p_V47_to_He4_Ti44_derived_removed = 283,
        k_Fe52_to_He4_Cr48_derived_removed = 284,
        k_Fe52_to_p_Mn51_derived_removed = 285,
        k_p_Mn51_to_He4_Cr48_derived_removed = 286,
        NumRates = k_p_Mn51_to_He4_Cr48_derived_removed
    };

    // number of reaclib rates

    const int NrateReaclib = 252;

    // number of tabular rates

    const int NrateTabular = 24;

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
        "He4_Na21_to_Al25",  // 30,
        "He4_Na22_to_Al26",  // 31,
        "p_Na23_to_Mg24",  // 32,
        "p_Mg24_to_Al25",  // 33,
        "He4_Mg24_to_Si28",  // 34,
        "He4_Al25_to_P29",  // 35,
        "He4_Al26_to_P30",  // 36,
        "p_Al27_to_Si28",  // 37,
        "p_Si28_to_P29",  // 38,
        "He4_Si28_to_S32",  // 39,
        "He4_P29_to_Cl33",  // 40,
        "He4_P30_to_Cl34",  // 41,
        "p_P31_to_S32",  // 42,
        "p_S32_to_Cl33",  // 43,
        "He4_Cl33_to_K37",  // 44,
        "He4_Cl34_to_K38",  // 45,
        "p_Ar36_to_K37",  // 46,
        "He4_K37_to_Sc41",  // 47,
        "He4_K38_to_Sc42",  // 48,
        "p_Ca40_to_Sc41",  // 49,
        "He4_Sc41_to_V45",  // 50,
        "He4_Sc42_to_V46",  // 51,
        "p_Ti44_to_V45",  // 52,
        "He4_V45_to_Mn49",  // 53,
        "He4_V46_to_Mn50",  // 54,
        "p_Cr48_to_Mn49",  // 55,
        "He4_Mn50_to_Co54",  // 56,
        "He4_Fe52_to_Ni56",  // 57,
        "p_Co55_to_Ni56",  // 58,
        "C12_C12_to_p_Na23",  // 59,
        "C12_C12_to_He4_Ne20",  // 60,
        "He4_N13_to_p_O16",  // 61,
        "p_N15_to_He4_C12",  // 62,
        "He4_O14_to_p_F17",  // 63,
        "C12_O16_to_p_Al27",  // 64,
        "C12_O16_to_He4_Mg24",  // 65,
        "O16_O16_to_p_P31",  // 66,
        "O16_O16_to_He4_Si28",  // 67,
        "p_O17_to_He4_N14",  // 68,
        "p_O18_to_He4_N15",  // 69,
        "p_F18_to_He4_O15",  // 70,
        "He4_F18_to_p_Ne21",  // 71,
        "p_F19_to_He4_O16",  // 72,
        "He4_F19_to_p_Ne22",  // 73,
        "He4_Ne18_to_p_Na21",  // 74,
        "He4_Ne19_to_p_Na22",  // 75,
        "p_Ne20_to_He4_F17",  // 76,
        "p_Na23_to_He4_Ne20",  // 77,
        "p_Mg24_to_He4_Na21",  // 78,
        "He4_Al25_to_p_Si28",  // 79,
        "p_Al27_to_He4_Mg24",  // 80,
        "He4_P29_to_p_S32",  // 81,
        "p_P31_to_He4_Si28",  // 82,
        "He4_Cl33_to_p_Ar36",  // 83,
        "He4_K37_to_p_Ca40",  // 84,
        "He4_Sc41_to_p_Ti44",  // 85,
        "He4_V45_to_p_Cr48",  // 86,
        "He4_Mn49_to_p_Fe52",  // 87,
        "He4_Fe52_to_p_Co55",  // 88,
        "He4_He4_He4_to_C12",  // 89,
        "C12_C12_to_Mg24_modified",  // 90,
        "O16_O16_to_S32_modified",  // 91,
        "C12_O16_to_Si28_modified",  // 92,
        "n_Fe52_to_Fe53",  // 93,
        "n_Fe53_to_Fe54",  // 94,
        "He4_Fe53_to_Ni57",  // 95,
        "n_Fe54_to_Fe55",  // 96,
        "p_Fe54_to_Co55",  // 97,
        "He4_Fe54_to_Ni58",  // 98,
        "n_Fe55_to_Fe56",  // 99,
        "p_Fe55_to_Co56",  // 100,
        "p_Fe56_to_Co57",  // 101,
        "n_Co55_to_Co56",  // 102,
        "n_Co56_to_Co57",  // 103,
        "p_Co56_to_Ni57",  // 104,
        "p_Co57_to_Ni58",  // 105,
        "n_Ni56_to_Ni57",  // 106,
        "n_Ni57_to_Ni58",  // 107,
        "He4_Fe53_to_p_Co56",  // 108,
        "n_Co55_to_p_Fe55",  // 109,
        "He4_Co55_to_p_Ni58",  // 110,
        "n_Co56_to_p_Fe56",  // 111,
        "p_Co57_to_He4_Fe54",  // 112,
        "n_Ni56_to_p_Co56",  // 113,
        "n_Ni56_to_He4_Fe53",  // 114,
        "n_Ni57_to_p_Co57",  // 115,
        "n_Ni57_to_He4_Fe54",  // 116,
        "n_Ni58_to_He4_Fe55",  // 117,
        "He4_S32_to_Ar36_removed",  // 118,
        "p_Cl35_to_Ar36_removed",  // 119,
        "p_Cl35_to_He4_S32_removed",  // 120,
        "He4_Ar36_to_Ca40_removed",  // 121,
        "p_K39_to_Ca40_removed",  // 122,
        "p_K39_to_He4_Ar36_removed",  // 123,
        "He4_Ca40_to_Ti44_removed",  // 124,
        "p_Sc43_to_Ti44_removed",  // 125,
        "p_Sc43_to_He4_Ca40_removed",  // 126,
        "He4_Ti44_to_Cr48_removed",  // 127,
        "He4_Ti44_to_p_V47_removed",  // 128,
        "p_V47_to_Cr48_removed",  // 129,
        "He4_Cr48_to_Fe52_removed",  // 130,
        "He4_Cr48_to_p_Mn51_removed",  // 131,
        "p_Mn51_to_Fe52_removed",  // 132,
        "Co54_to_Fe54",  // 133,
        "Fe54_to_Co54",  // 134,
        "Co55_to_Fe55",  // 135,
        "Fe55_to_Co55",  // 136,
        "Co56_to_Fe56",  // 137,
        "Co56_to_Ni56",  // 138,
        "Fe56_to_Co56",  // 139,
        "Ni56_to_Co56",  // 140,
        "Co57_to_Ni57",  // 141,
        "Ni57_to_Co57",  // 142,
        "n_to_p",  // 143,
        "p_to_n",  // 144,
        "F17_to_O17",  // 145,
        "O17_to_F17",  // 146,
        "F18_to_Ne18",  // 147,
        "F18_to_O18",  // 148,
        "Ne18_to_F18",  // 149,
        "O18_to_F18",  // 150,
        "F19_to_Ne19",  // 151,
        "Ne19_to_F19",  // 152,
        "Na21_to_Ne21",  // 153,
        "Ne21_to_Na21",  // 154,
        "Na22_to_Ne22",  // 155,
        "Ne22_to_Na22",  // 156,
        "S32_He4_to_Ar36_approx",  // 157,
        "Ar36_to_S32_He4_approx",  // 158,
        "Ar36_He4_to_Ca40_approx",  // 159,
        "Ca40_to_Ar36_He4_approx",  // 160,
        "Ca40_He4_to_Ti44_approx",  // 161,
        "Ti44_to_Ca40_He4_approx",  // 162,
        "Ti44_He4_to_Cr48_approx",  // 163,
        "Cr48_to_Ti44_He4_approx",  // 164,
        "Cr48_He4_to_Fe52_approx",  // 165,
        "Fe52_to_Cr48_He4_approx",  // 166,
        "N13_to_p_C12_derived",  // 167,
        "N14_to_p_C13_derived",  // 168,
        "O14_to_p_N13_derived",  // 169,
        "O15_to_p_N14_derived",  // 170,
        "O16_to_p_N15_derived",  // 171,
        "O16_to_He4_C12_derived",  // 172,
        "F17_to_p_O16_derived",  // 173,
        "F18_to_p_O17_derived",  // 174,
        "F18_to_He4_N14_derived",  // 175,
        "F19_to_p_O18_derived",  // 176,
        "F19_to_He4_N15_derived",  // 177,
        "Ne18_to_p_F17_derived",  // 178,
        "Ne18_to_He4_O14_derived",  // 179,
        "Ne19_to_p_F18_derived",  // 180,
        "Ne19_to_He4_O15_derived",  // 181,
        "Ne20_to_p_F19_derived",  // 182,
        "Ne20_to_He4_O16_derived",  // 183,
        "Ne21_to_He4_O17_derived",  // 184,
        "Ne22_to_He4_O18_derived",  // 185,
        "Na21_to_p_Ne20_derived",  // 186,
        "Na21_to_He4_F17_derived",  // 187,
        "Na22_to_p_Ne21_derived",  // 188,
        "Na22_to_He4_F18_derived",  // 189,
        "Na23_to_p_Ne22_derived",  // 190,
        "Na23_to_He4_F19_derived",  // 191,
        "Mg24_to_p_Na23_derived",  // 192,
        "Mg24_to_He4_Ne20_derived",  // 193,
        "Al25_to_p_Mg24_derived",  // 194,
        "Al25_to_He4_Na21_derived",  // 195,
        "Al26_to_He4_Na22_derived",  // 196,
        "Si28_to_p_Al27_derived",  // 197,
        "Si28_to_He4_Mg24_derived",  // 198,
        "P29_to_p_Si28_derived",  // 199,
        "P29_to_He4_Al25_derived",  // 200,
        "P30_to_He4_Al26_derived",  // 201,
        "S32_to_p_P31_derived",  // 202,
        "S32_to_He4_Si28_derived",  // 203,
        "Cl33_to_p_S32_derived",  // 204,
        "Cl33_to_He4_P29_derived",  // 205,
        "Cl34_to_He4_P30_derived",  // 206,
        "K37_to_p_Ar36_derived",  // 207,
        "K37_to_He4_Cl33_derived",  // 208,
        "K38_to_He4_Cl34_derived",  // 209,
        "Sc41_to_p_Ca40_derived",  // 210,
        "Sc41_to_He4_K37_derived",  // 211,
        "Sc42_to_He4_K38_derived",  // 212,
        "V45_to_p_Ti44_derived",  // 213,
        "V45_to_He4_Sc41_derived",  // 214,
        "V46_to_He4_Sc42_derived",  // 215,
        "Mn49_to_p_Cr48_derived",  // 216,
        "Mn49_to_He4_V45_derived",  // 217,
        "Mn50_to_He4_V46_derived",  // 218,
        "Co54_to_He4_Mn50_derived",  // 219,
        "Ni56_to_p_Co55_derived",  // 220,
        "Ni56_to_He4_Fe52_derived",  // 221,
        "C12_to_He4_He4_He4_derived",  // 222,
        "He4_C12_to_p_N15_derived",  // 223,
        "He4_N14_to_p_O17_derived",  // 224,
        "He4_N15_to_p_O18_derived",  // 225,
        "He4_O15_to_p_F18_derived",  // 226,
        "p_O16_to_He4_N13_derived",  // 227,
        "He4_O16_to_p_F19_derived",  // 228,
        "p_F17_to_He4_O14_derived",  // 229,
        "He4_F17_to_p_Ne20_derived",  // 230,
        "He4_Ne20_to_p_Na23_derived",  // 231,
        "p_Ne21_to_He4_F18_derived",  // 232,
        "p_Ne22_to_He4_F19_derived",  // 233,
        "p_Na21_to_He4_Ne18_derived",  // 234,
        "He4_Na21_to_p_Mg24_derived",  // 235,
        "p_Na22_to_He4_Ne19_derived",  // 236,
        "He4_Mg24_to_p_Al27_derived",  // 237,
        "p_Si28_to_He4_Al25_derived",  // 238,
        "He4_Si28_to_p_P31_derived",  // 239,
        "p_S32_to_He4_P29_derived",  // 240,
        "p_Ar36_to_He4_Cl33_derived",  // 241,
        "p_Ca40_to_He4_K37_derived",  // 242,
        "p_Ti44_to_He4_Sc41_derived",  // 243,
        "p_Cr48_to_He4_V45_derived",  // 244,
        "p_Fe52_to_He4_Mn49_derived",  // 245,
        "p_Co55_to_He4_Fe52_derived",  // 246,
        "Fe53_to_n_Fe52_derived",  // 247,
        "Fe54_to_n_Fe53_derived",  // 248,
        "Fe55_to_n_Fe54_derived",  // 249,
        "Fe56_to_n_Fe55_derived",  // 250,
        "Co55_to_p_Fe54_derived",  // 251,
        "Co56_to_n_Co55_derived",  // 252,
        "Co56_to_p_Fe55_derived",  // 253,
        "Co57_to_n_Co56_derived",  // 254,
        "Co57_to_p_Fe56_derived",  // 255,
        "Ni57_to_n_Ni56_derived",  // 256,
        "Ni57_to_p_Co56_derived",  // 257,
        "Ni57_to_He4_Fe53_derived",  // 258,
        "Ni58_to_n_Ni57_derived",  // 259,
        "Ni58_to_p_Co57_derived",  // 260,
        "Ni58_to_He4_Fe54_derived",  // 261,
        "He4_Fe53_to_n_Ni56_derived",  // 262,
        "He4_Fe54_to_n_Ni57_derived",  // 263,
        "He4_Fe54_to_p_Co57_derived",  // 264,
        "p_Fe55_to_n_Co55_derived",  // 265,
        "He4_Fe55_to_n_Ni58_derived",  // 266,
        "p_Fe56_to_n_Co56_derived",  // 267,
        "p_Co56_to_n_Ni56_derived",  // 268,
        "p_Co56_to_He4_Fe53_derived",  // 269,
        "p_Co57_to_n_Ni57_derived",  // 270,
        "p_Ni58_to_He4_Co55_derived",  // 271,
        "He4_S32_to_p_Cl35_derived_removed",  // 272,
        "Ar36_to_He4_S32_derived_removed",  // 273,
        "Ar36_to_p_Cl35_derived_removed",  // 274,
        "He4_Ar36_to_p_K39_derived_removed",  // 275,
        "Ca40_to_He4_Ar36_derived_removed",  // 276,
        "Ca40_to_p_K39_derived_removed",  // 277,
        "He4_Ca40_to_p_Sc43_derived_removed",  // 278,
        "Ti44_to_He4_Ca40_derived_removed",  // 279,
        "Ti44_to_p_Sc43_derived_removed",  // 280,
        "Cr48_to_He4_Ti44_derived_removed",  // 281,
        "Cr48_to_p_V47_derived_removed",  // 282,
        "p_V47_to_He4_Ti44_derived_removed",  // 283,
        "Fe52_to_He4_Cr48_derived_removed",  // 284,
        "Fe52_to_p_Mn51_derived_removed",  // 285,
        "p_Mn51_to_He4_Cr48_derived_removed"  // 286,
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
