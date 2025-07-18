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
        else if constexpr (spec == Mn51) {
            return 440.321747199996_rt;
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
        else if constexpr (spec == Mn51) {
            return 8.460149108098279e-23_rt;
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
        k_He4_Na23_to_Al27 = 33,
        k_p_Mg24_to_Al25 = 34,
        k_He4_Mg24_to_Si28 = 35,
        k_He4_Al25_to_P29 = 36,
        k_He4_Al26_to_P30 = 37,
        k_p_Al27_to_Si28 = 38,
        k_He4_Al27_to_P31 = 39,
        k_p_Si28_to_P29 = 40,
        k_He4_Si28_to_S32 = 41,
        k_He4_P29_to_Cl33 = 42,
        k_He4_P30_to_Cl34 = 43,
        k_p_P31_to_S32 = 44,
        k_p_S32_to_Cl33 = 45,
        k_He4_Cl33_to_K37 = 46,
        k_He4_Cl34_to_K38 = 47,
        k_p_Ar36_to_K37 = 48,
        k_He4_K37_to_Sc41 = 49,
        k_He4_K38_to_Sc42 = 50,
        k_p_Ca40_to_Sc41 = 51,
        k_He4_Sc41_to_V45 = 52,
        k_He4_Sc42_to_V46 = 53,
        k_p_Ti44_to_V45 = 54,
        k_He4_V45_to_Mn49 = 55,
        k_He4_V46_to_Mn50 = 56,
        k_p_Cr48_to_Mn49 = 57,
        k_He4_Cr48_to_Fe52 = 58,
        k_He4_Mn50_to_Co54 = 59,
        k_p_Mn51_to_Fe52 = 60,
        k_He4_Mn51_to_Co55 = 61,
        k_He4_Fe52_to_Ni56 = 62,
        k_p_Co55_to_Ni56 = 63,
        k_C12_C12_to_p_Na23 = 64,
        k_C12_C12_to_He4_Ne20 = 65,
        k_He4_N13_to_p_O16 = 66,
        k_p_N15_to_He4_C12 = 67,
        k_He4_O14_to_p_F17 = 68,
        k_C12_O16_to_p_Al27 = 69,
        k_C12_O16_to_He4_Mg24 = 70,
        k_O16_O16_to_p_P31 = 71,
        k_O16_O16_to_He4_Si28 = 72,
        k_p_O17_to_He4_N14 = 73,
        k_p_O18_to_He4_N15 = 74,
        k_p_F18_to_He4_O15 = 75,
        k_He4_F18_to_p_Ne21 = 76,
        k_p_F19_to_He4_O16 = 77,
        k_He4_F19_to_p_Ne22 = 78,
        k_He4_Ne18_to_p_Na21 = 79,
        k_He4_Ne19_to_p_Na22 = 80,
        k_p_Ne20_to_He4_F17 = 81,
        k_p_Na23_to_He4_Ne20 = 82,
        k_p_Mg24_to_He4_Na21 = 83,
        k_He4_Al25_to_p_Si28 = 84,
        k_p_Al27_to_He4_Mg24 = 85,
        k_He4_P29_to_p_S32 = 86,
        k_p_P31_to_He4_Si28 = 87,
        k_He4_Cl33_to_p_Ar36 = 88,
        k_He4_K37_to_p_Ca40 = 89,
        k_He4_Sc41_to_p_Ti44 = 90,
        k_He4_V45_to_p_Cr48 = 91,
        k_He4_Cr48_to_p_Mn51 = 92,
        k_He4_Mn49_to_p_Fe52 = 93,
        k_He4_Fe52_to_p_Co55 = 94,
        k_He4_He4_He4_to_C12 = 95,
        k_C12_C12_to_Mg24_modified = 96,
        k_O16_O16_to_S32_modified = 97,
        k_C12_O16_to_Si28_modified = 98,
        k_n_Mn49_to_Mn50 = 99,
        k_n_Mn50_to_Mn51 = 100,
        k_n_Fe52_to_Fe53 = 101,
        k_n_Fe53_to_Fe54 = 102,
        k_p_Fe53_to_Co54 = 103,
        k_He4_Fe53_to_Ni57 = 104,
        k_n_Fe54_to_Fe55 = 105,
        k_p_Fe54_to_Co55 = 106,
        k_He4_Fe54_to_Ni58 = 107,
        k_n_Fe55_to_Fe56 = 108,
        k_p_Fe55_to_Co56 = 109,
        k_p_Fe56_to_Co57 = 110,
        k_n_Co54_to_Co55 = 111,
        k_n_Co55_to_Co56 = 112,
        k_n_Co56_to_Co57 = 113,
        k_p_Co56_to_Ni57 = 114,
        k_p_Co57_to_Ni58 = 115,
        k_n_Ni56_to_Ni57 = 116,
        k_n_Ni57_to_Ni58 = 117,
        k_He4_Mn50_to_p_Fe53 = 118,
        k_He4_Mn51_to_p_Fe54 = 119,
        k_He4_Fe53_to_p_Co56 = 120,
        k_n_Co54_to_p_Fe54 = 121,
        k_n_Co54_to_He4_Mn51 = 122,
        k_He4_Co54_to_p_Ni57 = 123,
        k_n_Co55_to_p_Fe55 = 124,
        k_He4_Co55_to_p_Ni58 = 125,
        k_n_Co56_to_p_Fe56 = 126,
        k_p_Co57_to_He4_Fe54 = 127,
        k_n_Ni56_to_p_Co56 = 128,
        k_n_Ni56_to_He4_Fe53 = 129,
        k_n_Ni57_to_p_Co57 = 130,
        k_n_Ni57_to_He4_Fe54 = 131,
        k_n_Ni58_to_He4_Fe55 = 132,
        k_He4_S32_to_Ar36_removed = 133,
        k_p_Cl35_to_Ar36_removed = 134,
        k_p_Cl35_to_He4_S32_removed = 135,
        k_He4_Ar36_to_Ca40_removed = 136,
        k_p_K39_to_Ca40_removed = 137,
        k_p_K39_to_He4_Ar36_removed = 138,
        k_He4_Ca40_to_Ti44_removed = 139,
        k_p_Sc43_to_Ti44_removed = 140,
        k_p_Sc43_to_He4_Ca40_removed = 141,
        k_He4_Ti44_to_Cr48_removed = 142,
        k_He4_Ti44_to_p_V47_removed = 143,
        k_p_V47_to_Cr48_removed = 144,
        k_Co54_to_Fe54 = 145,
        k_Fe54_to_Co54 = 146,
        k_Co55_to_Fe55 = 147,
        k_Fe55_to_Co55 = 148,
        k_Co56_to_Fe56 = 149,
        k_Co56_to_Ni56 = 150,
        k_Fe56_to_Co56 = 151,
        k_Ni56_to_Co56 = 152,
        k_Co57_to_Ni57 = 153,
        k_Ni57_to_Co57 = 154,
        k_n_to_p = 155,
        k_p_to_n = 156,
        k_F17_to_O17 = 157,
        k_O17_to_F17 = 158,
        k_F18_to_Ne18 = 159,
        k_F18_to_O18 = 160,
        k_Ne18_to_F18 = 161,
        k_O18_to_F18 = 162,
        k_F19_to_Ne19 = 163,
        k_Ne19_to_F19 = 164,
        k_Na21_to_Ne21 = 165,
        k_Ne21_to_Na21 = 166,
        k_Na22_to_Ne22 = 167,
        k_Ne22_to_Na22 = 168,
        k_S32_He4_to_Ar36_approx = 169,
        k_Ar36_to_S32_He4_approx = 170,
        k_Ar36_He4_to_Ca40_approx = 171,
        k_Ca40_to_Ar36_He4_approx = 172,
        k_Ca40_He4_to_Ti44_approx = 173,
        k_Ti44_to_Ca40_He4_approx = 174,
        k_Ti44_He4_to_Cr48_approx = 175,
        k_Cr48_to_Ti44_He4_approx = 176,
        k_N13_to_p_C12_derived = 177,
        k_N14_to_p_C13_derived = 178,
        k_O14_to_p_N13_derived = 179,
        k_O15_to_p_N14_derived = 180,
        k_O16_to_p_N15_derived = 181,
        k_O16_to_He4_C12_derived = 182,
        k_F17_to_p_O16_derived = 183,
        k_F18_to_p_O17_derived = 184,
        k_F18_to_He4_N14_derived = 185,
        k_F19_to_p_O18_derived = 186,
        k_F19_to_He4_N15_derived = 187,
        k_Ne18_to_p_F17_derived = 188,
        k_Ne18_to_He4_O14_derived = 189,
        k_Ne19_to_p_F18_derived = 190,
        k_Ne19_to_He4_O15_derived = 191,
        k_Ne20_to_p_F19_derived = 192,
        k_Ne20_to_He4_O16_derived = 193,
        k_Ne21_to_He4_O17_derived = 194,
        k_Ne22_to_He4_O18_derived = 195,
        k_Na21_to_p_Ne20_derived = 196,
        k_Na21_to_He4_F17_derived = 197,
        k_Na22_to_p_Ne21_derived = 198,
        k_Na22_to_He4_F18_derived = 199,
        k_Na23_to_p_Ne22_derived = 200,
        k_Na23_to_He4_F19_derived = 201,
        k_Mg24_to_p_Na23_derived = 202,
        k_Mg24_to_He4_Ne20_derived = 203,
        k_Al25_to_p_Mg24_derived = 204,
        k_Al25_to_He4_Na21_derived = 205,
        k_Al26_to_He4_Na22_derived = 206,
        k_Al27_to_He4_Na23_derived = 207,
        k_Si28_to_p_Al27_derived = 208,
        k_Si28_to_He4_Mg24_derived = 209,
        k_P29_to_p_Si28_derived = 210,
        k_P29_to_He4_Al25_derived = 211,
        k_P30_to_He4_Al26_derived = 212,
        k_P31_to_He4_Al27_derived = 213,
        k_S32_to_p_P31_derived = 214,
        k_S32_to_He4_Si28_derived = 215,
        k_Cl33_to_p_S32_derived = 216,
        k_Cl33_to_He4_P29_derived = 217,
        k_Cl34_to_He4_P30_derived = 218,
        k_K37_to_p_Ar36_derived = 219,
        k_K37_to_He4_Cl33_derived = 220,
        k_K38_to_He4_Cl34_derived = 221,
        k_Sc41_to_p_Ca40_derived = 222,
        k_Sc41_to_He4_K37_derived = 223,
        k_Sc42_to_He4_K38_derived = 224,
        k_V45_to_p_Ti44_derived = 225,
        k_V45_to_He4_Sc41_derived = 226,
        k_V46_to_He4_Sc42_derived = 227,
        k_Mn49_to_p_Cr48_derived = 228,
        k_Mn49_to_He4_V45_derived = 229,
        k_Mn50_to_He4_V46_derived = 230,
        k_Fe52_to_p_Mn51_derived = 231,
        k_Fe52_to_He4_Cr48_derived = 232,
        k_Co54_to_He4_Mn50_derived = 233,
        k_Co55_to_He4_Mn51_derived = 234,
        k_Ni56_to_p_Co55_derived = 235,
        k_Ni56_to_He4_Fe52_derived = 236,
        k_C12_to_He4_He4_He4_derived = 237,
        k_He4_C12_to_p_N15_derived = 238,
        k_He4_N14_to_p_O17_derived = 239,
        k_He4_N15_to_p_O18_derived = 240,
        k_He4_O15_to_p_F18_derived = 241,
        k_p_O16_to_He4_N13_derived = 242,
        k_He4_O16_to_p_F19_derived = 243,
        k_p_F17_to_He4_O14_derived = 244,
        k_He4_F17_to_p_Ne20_derived = 245,
        k_He4_Ne20_to_p_Na23_derived = 246,
        k_p_Ne21_to_He4_F18_derived = 247,
        k_p_Ne22_to_He4_F19_derived = 248,
        k_p_Na21_to_He4_Ne18_derived = 249,
        k_He4_Na21_to_p_Mg24_derived = 250,
        k_p_Na22_to_He4_Ne19_derived = 251,
        k_He4_Mg24_to_p_Al27_derived = 252,
        k_p_Si28_to_He4_Al25_derived = 253,
        k_He4_Si28_to_p_P31_derived = 254,
        k_p_S32_to_He4_P29_derived = 255,
        k_p_Ar36_to_He4_Cl33_derived = 256,
        k_p_Ca40_to_He4_K37_derived = 257,
        k_p_Ti44_to_He4_Sc41_derived = 258,
        k_p_Cr48_to_He4_V45_derived = 259,
        k_p_Mn51_to_He4_Cr48_derived = 260,
        k_p_Fe52_to_He4_Mn49_derived = 261,
        k_p_Co55_to_He4_Fe52_derived = 262,
        k_Mn50_to_n_Mn49_derived = 263,
        k_Mn51_to_n_Mn50_derived = 264,
        k_Fe53_to_n_Fe52_derived = 265,
        k_Fe54_to_n_Fe53_derived = 266,
        k_Fe55_to_n_Fe54_derived = 267,
        k_Fe56_to_n_Fe55_derived = 268,
        k_Co54_to_p_Fe53_derived = 269,
        k_Co55_to_n_Co54_derived = 270,
        k_Co55_to_p_Fe54_derived = 271,
        k_Co56_to_n_Co55_derived = 272,
        k_Co56_to_p_Fe55_derived = 273,
        k_Co57_to_n_Co56_derived = 274,
        k_Co57_to_p_Fe56_derived = 275,
        k_Ni57_to_n_Ni56_derived = 276,
        k_Ni57_to_p_Co56_derived = 277,
        k_Ni57_to_He4_Fe53_derived = 278,
        k_Ni58_to_n_Ni57_derived = 279,
        k_Ni58_to_p_Co57_derived = 280,
        k_Ni58_to_He4_Fe54_derived = 281,
        k_He4_Mn51_to_n_Co54_derived = 282,
        k_p_Fe53_to_He4_Mn50_derived = 283,
        k_He4_Fe53_to_n_Ni56_derived = 284,
        k_p_Fe54_to_n_Co54_derived = 285,
        k_p_Fe54_to_He4_Mn51_derived = 286,
        k_He4_Fe54_to_n_Ni57_derived = 287,
        k_He4_Fe54_to_p_Co57_derived = 288,
        k_p_Fe55_to_n_Co55_derived = 289,
        k_He4_Fe55_to_n_Ni58_derived = 290,
        k_p_Fe56_to_n_Co56_derived = 291,
        k_p_Co56_to_n_Ni56_derived = 292,
        k_p_Co56_to_He4_Fe53_derived = 293,
        k_p_Co57_to_n_Ni57_derived = 294,
        k_p_Ni57_to_He4_Co54_derived = 295,
        k_p_Ni58_to_He4_Co55_derived = 296,
        k_He4_S32_to_p_Cl35_derived_removed = 297,
        k_Ar36_to_He4_S32_derived_removed = 298,
        k_Ar36_to_p_Cl35_derived_removed = 299,
        k_He4_Ar36_to_p_K39_derived_removed = 300,
        k_Ca40_to_He4_Ar36_derived_removed = 301,
        k_Ca40_to_p_K39_derived_removed = 302,
        k_He4_Ca40_to_p_Sc43_derived_removed = 303,
        k_Ti44_to_He4_Ca40_derived_removed = 304,
        k_Ti44_to_p_Sc43_derived_removed = 305,
        k_Cr48_to_He4_Ti44_derived_removed = 306,
        k_Cr48_to_p_V47_derived_removed = 307,
        k_p_V47_to_He4_Ti44_derived_removed = 308,
        NumRates = k_p_V47_to_He4_Ti44_derived_removed
    };

    // number of reaclib rates

    const int NrateReaclib = 276;

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
        "He4_Na23_to_Al27",  // 33,
        "p_Mg24_to_Al25",  // 34,
        "He4_Mg24_to_Si28",  // 35,
        "He4_Al25_to_P29",  // 36,
        "He4_Al26_to_P30",  // 37,
        "p_Al27_to_Si28",  // 38,
        "He4_Al27_to_P31",  // 39,
        "p_Si28_to_P29",  // 40,
        "He4_Si28_to_S32",  // 41,
        "He4_P29_to_Cl33",  // 42,
        "He4_P30_to_Cl34",  // 43,
        "p_P31_to_S32",  // 44,
        "p_S32_to_Cl33",  // 45,
        "He4_Cl33_to_K37",  // 46,
        "He4_Cl34_to_K38",  // 47,
        "p_Ar36_to_K37",  // 48,
        "He4_K37_to_Sc41",  // 49,
        "He4_K38_to_Sc42",  // 50,
        "p_Ca40_to_Sc41",  // 51,
        "He4_Sc41_to_V45",  // 52,
        "He4_Sc42_to_V46",  // 53,
        "p_Ti44_to_V45",  // 54,
        "He4_V45_to_Mn49",  // 55,
        "He4_V46_to_Mn50",  // 56,
        "p_Cr48_to_Mn49",  // 57,
        "He4_Cr48_to_Fe52",  // 58,
        "He4_Mn50_to_Co54",  // 59,
        "p_Mn51_to_Fe52",  // 60,
        "He4_Mn51_to_Co55",  // 61,
        "He4_Fe52_to_Ni56",  // 62,
        "p_Co55_to_Ni56",  // 63,
        "C12_C12_to_p_Na23",  // 64,
        "C12_C12_to_He4_Ne20",  // 65,
        "He4_N13_to_p_O16",  // 66,
        "p_N15_to_He4_C12",  // 67,
        "He4_O14_to_p_F17",  // 68,
        "C12_O16_to_p_Al27",  // 69,
        "C12_O16_to_He4_Mg24",  // 70,
        "O16_O16_to_p_P31",  // 71,
        "O16_O16_to_He4_Si28",  // 72,
        "p_O17_to_He4_N14",  // 73,
        "p_O18_to_He4_N15",  // 74,
        "p_F18_to_He4_O15",  // 75,
        "He4_F18_to_p_Ne21",  // 76,
        "p_F19_to_He4_O16",  // 77,
        "He4_F19_to_p_Ne22",  // 78,
        "He4_Ne18_to_p_Na21",  // 79,
        "He4_Ne19_to_p_Na22",  // 80,
        "p_Ne20_to_He4_F17",  // 81,
        "p_Na23_to_He4_Ne20",  // 82,
        "p_Mg24_to_He4_Na21",  // 83,
        "He4_Al25_to_p_Si28",  // 84,
        "p_Al27_to_He4_Mg24",  // 85,
        "He4_P29_to_p_S32",  // 86,
        "p_P31_to_He4_Si28",  // 87,
        "He4_Cl33_to_p_Ar36",  // 88,
        "He4_K37_to_p_Ca40",  // 89,
        "He4_Sc41_to_p_Ti44",  // 90,
        "He4_V45_to_p_Cr48",  // 91,
        "He4_Cr48_to_p_Mn51",  // 92,
        "He4_Mn49_to_p_Fe52",  // 93,
        "He4_Fe52_to_p_Co55",  // 94,
        "He4_He4_He4_to_C12",  // 95,
        "C12_C12_to_Mg24_modified",  // 96,
        "O16_O16_to_S32_modified",  // 97,
        "C12_O16_to_Si28_modified",  // 98,
        "n_Mn49_to_Mn50",  // 99,
        "n_Mn50_to_Mn51",  // 100,
        "n_Fe52_to_Fe53",  // 101,
        "n_Fe53_to_Fe54",  // 102,
        "p_Fe53_to_Co54",  // 103,
        "He4_Fe53_to_Ni57",  // 104,
        "n_Fe54_to_Fe55",  // 105,
        "p_Fe54_to_Co55",  // 106,
        "He4_Fe54_to_Ni58",  // 107,
        "n_Fe55_to_Fe56",  // 108,
        "p_Fe55_to_Co56",  // 109,
        "p_Fe56_to_Co57",  // 110,
        "n_Co54_to_Co55",  // 111,
        "n_Co55_to_Co56",  // 112,
        "n_Co56_to_Co57",  // 113,
        "p_Co56_to_Ni57",  // 114,
        "p_Co57_to_Ni58",  // 115,
        "n_Ni56_to_Ni57",  // 116,
        "n_Ni57_to_Ni58",  // 117,
        "He4_Mn50_to_p_Fe53",  // 118,
        "He4_Mn51_to_p_Fe54",  // 119,
        "He4_Fe53_to_p_Co56",  // 120,
        "n_Co54_to_p_Fe54",  // 121,
        "n_Co54_to_He4_Mn51",  // 122,
        "He4_Co54_to_p_Ni57",  // 123,
        "n_Co55_to_p_Fe55",  // 124,
        "He4_Co55_to_p_Ni58",  // 125,
        "n_Co56_to_p_Fe56",  // 126,
        "p_Co57_to_He4_Fe54",  // 127,
        "n_Ni56_to_p_Co56",  // 128,
        "n_Ni56_to_He4_Fe53",  // 129,
        "n_Ni57_to_p_Co57",  // 130,
        "n_Ni57_to_He4_Fe54",  // 131,
        "n_Ni58_to_He4_Fe55",  // 132,
        "He4_S32_to_Ar36_removed",  // 133,
        "p_Cl35_to_Ar36_removed",  // 134,
        "p_Cl35_to_He4_S32_removed",  // 135,
        "He4_Ar36_to_Ca40_removed",  // 136,
        "p_K39_to_Ca40_removed",  // 137,
        "p_K39_to_He4_Ar36_removed",  // 138,
        "He4_Ca40_to_Ti44_removed",  // 139,
        "p_Sc43_to_Ti44_removed",  // 140,
        "p_Sc43_to_He4_Ca40_removed",  // 141,
        "He4_Ti44_to_Cr48_removed",  // 142,
        "He4_Ti44_to_p_V47_removed",  // 143,
        "p_V47_to_Cr48_removed",  // 144,
        "Co54_to_Fe54",  // 145,
        "Fe54_to_Co54",  // 146,
        "Co55_to_Fe55",  // 147,
        "Fe55_to_Co55",  // 148,
        "Co56_to_Fe56",  // 149,
        "Co56_to_Ni56",  // 150,
        "Fe56_to_Co56",  // 151,
        "Ni56_to_Co56",  // 152,
        "Co57_to_Ni57",  // 153,
        "Ni57_to_Co57",  // 154,
        "n_to_p",  // 155,
        "p_to_n",  // 156,
        "F17_to_O17",  // 157,
        "O17_to_F17",  // 158,
        "F18_to_Ne18",  // 159,
        "F18_to_O18",  // 160,
        "Ne18_to_F18",  // 161,
        "O18_to_F18",  // 162,
        "F19_to_Ne19",  // 163,
        "Ne19_to_F19",  // 164,
        "Na21_to_Ne21",  // 165,
        "Ne21_to_Na21",  // 166,
        "Na22_to_Ne22",  // 167,
        "Ne22_to_Na22",  // 168,
        "S32_He4_to_Ar36_approx",  // 169,
        "Ar36_to_S32_He4_approx",  // 170,
        "Ar36_He4_to_Ca40_approx",  // 171,
        "Ca40_to_Ar36_He4_approx",  // 172,
        "Ca40_He4_to_Ti44_approx",  // 173,
        "Ti44_to_Ca40_He4_approx",  // 174,
        "Ti44_He4_to_Cr48_approx",  // 175,
        "Cr48_to_Ti44_He4_approx",  // 176,
        "N13_to_p_C12_derived",  // 177,
        "N14_to_p_C13_derived",  // 178,
        "O14_to_p_N13_derived",  // 179,
        "O15_to_p_N14_derived",  // 180,
        "O16_to_p_N15_derived",  // 181,
        "O16_to_He4_C12_derived",  // 182,
        "F17_to_p_O16_derived",  // 183,
        "F18_to_p_O17_derived",  // 184,
        "F18_to_He4_N14_derived",  // 185,
        "F19_to_p_O18_derived",  // 186,
        "F19_to_He4_N15_derived",  // 187,
        "Ne18_to_p_F17_derived",  // 188,
        "Ne18_to_He4_O14_derived",  // 189,
        "Ne19_to_p_F18_derived",  // 190,
        "Ne19_to_He4_O15_derived",  // 191,
        "Ne20_to_p_F19_derived",  // 192,
        "Ne20_to_He4_O16_derived",  // 193,
        "Ne21_to_He4_O17_derived",  // 194,
        "Ne22_to_He4_O18_derived",  // 195,
        "Na21_to_p_Ne20_derived",  // 196,
        "Na21_to_He4_F17_derived",  // 197,
        "Na22_to_p_Ne21_derived",  // 198,
        "Na22_to_He4_F18_derived",  // 199,
        "Na23_to_p_Ne22_derived",  // 200,
        "Na23_to_He4_F19_derived",  // 201,
        "Mg24_to_p_Na23_derived",  // 202,
        "Mg24_to_He4_Ne20_derived",  // 203,
        "Al25_to_p_Mg24_derived",  // 204,
        "Al25_to_He4_Na21_derived",  // 205,
        "Al26_to_He4_Na22_derived",  // 206,
        "Al27_to_He4_Na23_derived",  // 207,
        "Si28_to_p_Al27_derived",  // 208,
        "Si28_to_He4_Mg24_derived",  // 209,
        "P29_to_p_Si28_derived",  // 210,
        "P29_to_He4_Al25_derived",  // 211,
        "P30_to_He4_Al26_derived",  // 212,
        "P31_to_He4_Al27_derived",  // 213,
        "S32_to_p_P31_derived",  // 214,
        "S32_to_He4_Si28_derived",  // 215,
        "Cl33_to_p_S32_derived",  // 216,
        "Cl33_to_He4_P29_derived",  // 217,
        "Cl34_to_He4_P30_derived",  // 218,
        "K37_to_p_Ar36_derived",  // 219,
        "K37_to_He4_Cl33_derived",  // 220,
        "K38_to_He4_Cl34_derived",  // 221,
        "Sc41_to_p_Ca40_derived",  // 222,
        "Sc41_to_He4_K37_derived",  // 223,
        "Sc42_to_He4_K38_derived",  // 224,
        "V45_to_p_Ti44_derived",  // 225,
        "V45_to_He4_Sc41_derived",  // 226,
        "V46_to_He4_Sc42_derived",  // 227,
        "Mn49_to_p_Cr48_derived",  // 228,
        "Mn49_to_He4_V45_derived",  // 229,
        "Mn50_to_He4_V46_derived",  // 230,
        "Fe52_to_p_Mn51_derived",  // 231,
        "Fe52_to_He4_Cr48_derived",  // 232,
        "Co54_to_He4_Mn50_derived",  // 233,
        "Co55_to_He4_Mn51_derived",  // 234,
        "Ni56_to_p_Co55_derived",  // 235,
        "Ni56_to_He4_Fe52_derived",  // 236,
        "C12_to_He4_He4_He4_derived",  // 237,
        "He4_C12_to_p_N15_derived",  // 238,
        "He4_N14_to_p_O17_derived",  // 239,
        "He4_N15_to_p_O18_derived",  // 240,
        "He4_O15_to_p_F18_derived",  // 241,
        "p_O16_to_He4_N13_derived",  // 242,
        "He4_O16_to_p_F19_derived",  // 243,
        "p_F17_to_He4_O14_derived",  // 244,
        "He4_F17_to_p_Ne20_derived",  // 245,
        "He4_Ne20_to_p_Na23_derived",  // 246,
        "p_Ne21_to_He4_F18_derived",  // 247,
        "p_Ne22_to_He4_F19_derived",  // 248,
        "p_Na21_to_He4_Ne18_derived",  // 249,
        "He4_Na21_to_p_Mg24_derived",  // 250,
        "p_Na22_to_He4_Ne19_derived",  // 251,
        "He4_Mg24_to_p_Al27_derived",  // 252,
        "p_Si28_to_He4_Al25_derived",  // 253,
        "He4_Si28_to_p_P31_derived",  // 254,
        "p_S32_to_He4_P29_derived",  // 255,
        "p_Ar36_to_He4_Cl33_derived",  // 256,
        "p_Ca40_to_He4_K37_derived",  // 257,
        "p_Ti44_to_He4_Sc41_derived",  // 258,
        "p_Cr48_to_He4_V45_derived",  // 259,
        "p_Mn51_to_He4_Cr48_derived",  // 260,
        "p_Fe52_to_He4_Mn49_derived",  // 261,
        "p_Co55_to_He4_Fe52_derived",  // 262,
        "Mn50_to_n_Mn49_derived",  // 263,
        "Mn51_to_n_Mn50_derived",  // 264,
        "Fe53_to_n_Fe52_derived",  // 265,
        "Fe54_to_n_Fe53_derived",  // 266,
        "Fe55_to_n_Fe54_derived",  // 267,
        "Fe56_to_n_Fe55_derived",  // 268,
        "Co54_to_p_Fe53_derived",  // 269,
        "Co55_to_n_Co54_derived",  // 270,
        "Co55_to_p_Fe54_derived",  // 271,
        "Co56_to_n_Co55_derived",  // 272,
        "Co56_to_p_Fe55_derived",  // 273,
        "Co57_to_n_Co56_derived",  // 274,
        "Co57_to_p_Fe56_derived",  // 275,
        "Ni57_to_n_Ni56_derived",  // 276,
        "Ni57_to_p_Co56_derived",  // 277,
        "Ni57_to_He4_Fe53_derived",  // 278,
        "Ni58_to_n_Ni57_derived",  // 279,
        "Ni58_to_p_Co57_derived",  // 280,
        "Ni58_to_He4_Fe54_derived",  // 281,
        "He4_Mn51_to_n_Co54_derived",  // 282,
        "p_Fe53_to_He4_Mn50_derived",  // 283,
        "He4_Fe53_to_n_Ni56_derived",  // 284,
        "p_Fe54_to_n_Co54_derived",  // 285,
        "p_Fe54_to_He4_Mn51_derived",  // 286,
        "He4_Fe54_to_n_Ni57_derived",  // 287,
        "He4_Fe54_to_p_Co57_derived",  // 288,
        "p_Fe55_to_n_Co55_derived",  // 289,
        "He4_Fe55_to_n_Ni58_derived",  // 290,
        "p_Fe56_to_n_Co56_derived",  // 291,
        "p_Co56_to_n_Ni56_derived",  // 292,
        "p_Co56_to_He4_Fe53_derived",  // 293,
        "p_Co57_to_n_Ni57_derived",  // 294,
        "p_Ni57_to_He4_Co54_derived",  // 295,
        "p_Ni58_to_He4_Co55_derived",  // 296,
        "He4_S32_to_p_Cl35_derived_removed",  // 297,
        "Ar36_to_He4_S32_derived_removed",  // 298,
        "Ar36_to_p_Cl35_derived_removed",  // 299,
        "He4_Ar36_to_p_K39_derived_removed",  // 300,
        "Ca40_to_He4_Ar36_derived_removed",  // 301,
        "Ca40_to_p_K39_derived_removed",  // 302,
        "He4_Ca40_to_p_Sc43_derived_removed",  // 303,
        "Ti44_to_He4_Ca40_derived_removed",  // 304,
        "Ti44_to_p_Sc43_derived_removed",  // 305,
        "Cr48_to_He4_Ti44_derived_removed",  // 306,
        "Cr48_to_p_V47_derived_removed",  // 307,
        "p_V47_to_He4_Ti44_derived_removed"  // 308,
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
