#ifndef actual_network_H
#define actual_network_H

#include <random>
#include <string_view>

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Loop.H>

#include <extern_parameters.H>
#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex::literals;


inline
void actual_network_init() {


}

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
        else if constexpr (spec == C14) {
            return 105.28447818399945_rt;
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
        else if constexpr (spec == O16) {
            return 127.6193154119992_rt;
        }
        else if constexpr (spec == O17) {
            return 131.76239561199873_rt;
        }
        else if constexpr (spec == O18) {
            return 139.8077658120019_rt;
        }
        else if constexpr (spec == F18) {
            return 137.36950247599816_rt;
        }
        else if constexpr (spec == F19) {
            return 147.80136567599766_rt;
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
        else if constexpr (spec == Na22) {
            return 174.14457080400098_rt;
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
        else if constexpr (spec == Al26) {
            return 211.89389913199557_rt;
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
        else if constexpr (spec == P30) {
            return 250.60519745999775_rt;
        }
        else if constexpr (spec == P31) {
            return 262.9161999600037_rt;
        }
        else if constexpr (spec == S32) {
            return 271.78016372399725_rt;
        }
        else if constexpr (spec == S33) {
            return 280.42180302400084_rt;
        }
        else if constexpr (spec == Cl34) {
            return 285.5650057879975_rt;
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
        else if constexpr (spec == K38) {
            return 320.64625411599985_rt;
        }
        else if constexpr (spec == K39) {
            return 333.7240072160057_rt;
        }
        else if constexpr (spec == Ca40) {
            return 342.05218528000114_rt;
        }
        else if constexpr (spec == Ca41) {
            return 350.4150113799915_rt;
        }
        else if constexpr (spec == Sc42) {
            return 354.68707244399405_rt;
        }
        else if constexpr (spec == Sc43) {
            return 366.82549054399715_rt;
        }
        else if constexpr (spec == Ti44) {
            return 375.47496160800074_rt;
        }
        else if constexpr (spec == Ti45) {
            return 385.0079797080034_rt;
        }
        else if constexpr (spec == V46) {
            return 390.3625507720062_rt;
        }
        else if constexpr (spec == V47) {
            return 403.3650388719979_rt;
        }
        else if constexpr (spec == Cr48) {
            return 411.4679399359957_rt;
        }
        else if constexpr (spec == Cr49) {
            return 422.05065803600155_rt;
        }
        else if constexpr (spec == Mn50) {
            return 426.63411910000286_rt;
        }
        else if constexpr (spec == Mn51) {
            return 440.321747199996_rt;
        }
        else if constexpr (spec == Mn52) {
            return 450.8612552999912_rt;
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
        else if constexpr (spec == Cu59) {
            return 509.8782038560021_rt;
        }
        else if constexpr (spec == Zn60) {
            return 514.9831749199948_rt;
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
        else if constexpr (spec == C14) {
            return 2.3252930380658474e-23_rt;
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
        else if constexpr (spec == O16) {
            return 2.6560180592333686e-23_rt;
        }
        else if constexpr (spec == O17) {
            return 2.822772237905709e-23_rt;
        }
        else if constexpr (spec == O18) {
            return 2.9888307701990364e-23_rt;
        }
        else if constexpr (spec == F18) {
            return 2.989125964092377e-23_rt;
        }
        else if constexpr (spec == F19) {
            return 3.1547590653659727e-23_rt;
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
        else if constexpr (spec == Na22) {
            return 3.652262279854593e-23_rt;
        }
        else if constexpr (spec == Na23) {
            return 3.817541002484691e-23_rt;
        }
        else if constexpr (spec == Mg24) {
            return 3.9828098739467446e-23_rt;
        }
        else if constexpr (spec == Mg25) {
            return 4.148995839845366e-23_rt;
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
        else if constexpr (spec == Si29) {
            return 4.811660140333913e-23_rt;
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
        else if constexpr (spec == S33) {
            return 5.475039560055494e-23_rt;
        }
        else if constexpr (spec == Cl34) {
            return 5.641475984659597e-23_rt;
        }
        else if constexpr (spec == Cl35) {
            return 5.806714601585715e-23_rt;
        }
        else if constexpr (spec == Ar36) {
            return 5.972551377884467e-23_rt;
        }
        else if constexpr (spec == Ar37) {
            return 6.138477621350144e-23_rt;
        }
        else if constexpr (spec == K38) {
            return 6.304914251263421e-23_rt;
        }
        else if constexpr (spec == K39) {
            return 6.470075679819664e-23_rt;
        }
        else if constexpr (spec == Ca40) {
            return 6.635944331004904e-23_rt;
        }
        else if constexpr (spec == Ca41) {
            return 6.801946271643761e-23_rt;
        }
        else if constexpr (spec == Sc42) {
            return 6.968537991361557e-23_rt;
        }
        else if constexpr (spec == Sc43) {
            return 7.133866871591414e-23_rt;
        }
        else if constexpr (spec == Ti44) {
            return 7.299678247096977e-23_rt;
        }
        else if constexpr (spec == Ti45) {
            return 7.465471582063894e-23_rt;
        }
        else if constexpr (spec == V46) {
            return 7.631870326846012e-23_rt;
        }
        else if constexpr (spec == V47) {
            return 7.797045172607207e-23_rt;
        }
        else if constexpr (spec == Cr48) {
            return 7.962953983065421e-23_rt;
        }
        else if constexpr (spec == Cr49) {
            return 8.128560192010424e-23_rt;
        }
        else if constexpr (spec == Mn50) {
            return 8.295096399635981e-23_rt;
        }
        else if constexpr (spec == Mn51) {
            return 8.460149108098279e-23_rt;
        }
        else if constexpr (spec == Mn52) {
            return 8.625763019925445e-23_rt;
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
        else if constexpr (spec == Cu59) {
            return 9.787133677749025e-23_rt;
        }
        else if constexpr (spec == Zn60) {
            return 9.953576917772707e-23_rt;
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
        k_C14_to_N14_reaclib = 1,
        k_N13_to_C13_reaclib = 2,
        k_p_C12_to_N13_reaclib = 3,
        k_He4_C12_to_O16_reaclib = 4,
        k_p_C13_to_N14_reaclib = 5,
        k_p_C14_to_N15_reaclib = 6,
        k_He4_C14_to_O18_reaclib = 7,
        k_He4_N14_to_F18_reaclib = 8,
        k_p_N15_to_O16_reaclib = 9,
        k_He4_N15_to_F19_reaclib = 10,
        k_He4_O16_to_Ne20_reaclib = 11,
        k_p_O17_to_F18_reaclib = 12,
        k_He4_O17_to_Ne21_reaclib = 13,
        k_p_O18_to_F19_reaclib = 14,
        k_He4_O18_to_Ne22_reaclib = 15,
        k_He4_F18_to_Na22_reaclib = 16,
        k_p_F19_to_Ne20_reaclib = 17,
        k_He4_F19_to_Na23_reaclib = 18,
        k_He4_Ne20_to_Mg24_reaclib = 19,
        k_p_Ne21_to_Na22_reaclib = 20,
        k_He4_Ne21_to_Mg25_reaclib = 21,
        k_p_Ne22_to_Na23_reaclib = 22,
        k_He4_Na22_to_Al26_reaclib = 23,
        k_p_Na23_to_Mg24_reaclib = 24,
        k_He4_Na23_to_Al27_reaclib = 25,
        k_He4_Mg24_to_Si28_reaclib = 26,
        k_p_Mg25_to_Al26_reaclib = 27,
        k_He4_Mg25_to_Si29_reaclib = 28,
        k_He4_Al26_to_P30_reaclib = 29,
        k_p_Al27_to_Si28_reaclib = 30,
        k_He4_Al27_to_P31_reaclib = 31,
        k_He4_Si28_to_S32_reaclib = 32,
        k_p_Si29_to_P30_reaclib = 33,
        k_He4_Si29_to_S33_reaclib = 34,
        k_He4_P30_to_Cl34_reaclib = 35,
        k_p_P31_to_S32_reaclib = 36,
        k_He4_P31_to_Cl35_reaclib = 37,
        k_He4_S32_to_Ar36_reaclib = 38,
        k_p_S33_to_Cl34_reaclib = 39,
        k_He4_S33_to_Ar37_reaclib = 40,
        k_He4_Cl34_to_K38_reaclib = 41,
        k_p_Cl35_to_Ar36_reaclib = 42,
        k_He4_Cl35_to_K39_reaclib = 43,
        k_He4_Ar36_to_Ca40_reaclib = 44,
        k_p_Ar37_to_K38_reaclib = 45,
        k_He4_Ar37_to_Ca41_reaclib = 46,
        k_He4_K38_to_Sc42_reaclib = 47,
        k_p_K39_to_Ca40_reaclib = 48,
        k_He4_K39_to_Sc43_reaclib = 49,
        k_He4_Ca40_to_Ti44_reaclib = 50,
        k_p_Ca41_to_Sc42_reaclib = 51,
        k_He4_Ca41_to_Ti45_reaclib = 52,
        k_He4_Sc42_to_V46_reaclib = 53,
        k_p_Sc43_to_Ti44_reaclib = 54,
        k_He4_Sc43_to_V47_reaclib = 55,
        k_He4_Ti44_to_Cr48_reaclib = 56,
        k_p_Ti45_to_V46_reaclib = 57,
        k_He4_Ti45_to_Cr49_reaclib = 58,
        k_He4_V46_to_Mn50_reaclib = 59,
        k_p_V47_to_Cr48_reaclib = 60,
        k_He4_V47_to_Mn51_reaclib = 61,
        k_He4_Cr48_to_Fe52_reaclib = 62,
        k_p_Cr49_to_Mn50_reaclib = 63,
        k_p_Mn51_to_Fe52_reaclib = 64,
        k_He4_Mn51_to_Co55_reaclib = 65,
        k_He4_Fe52_to_Ni56_reaclib = 66,
        k_p_Co55_to_Ni56_reaclib = 67,
        k_C12_C12_to_p_Na23_reaclib = 68,
        k_C12_C12_to_He4_Ne20_reaclib = 69,
        k_He4_N13_to_p_O16_reaclib = 70,
        k_p_N15_to_He4_C12_reaclib = 71,
        k_C12_O16_to_p_Al27_reaclib = 72,
        k_C12_O16_to_He4_Mg24_reaclib = 73,
        k_O16_O16_to_p_P31_reaclib = 74,
        k_O16_O16_to_He4_Si28_reaclib = 75,
        k_p_O17_to_He4_N14_reaclib = 76,
        k_p_O18_to_He4_N15_reaclib = 77,
        k_He4_F18_to_p_Ne21_reaclib = 78,
        k_p_F19_to_He4_O16_reaclib = 79,
        k_He4_F19_to_p_Ne22_reaclib = 80,
        k_C12_Ne20_to_p_P31_reaclib = 81,
        k_C12_Ne20_to_He4_Si28_reaclib = 82,
        k_He4_Na22_to_p_Mg25_reaclib = 83,
        k_p_Na23_to_He4_Ne20_reaclib = 84,
        k_He4_Al26_to_p_Si29_reaclib = 85,
        k_p_Al27_to_He4_Mg24_reaclib = 86,
        k_He4_P30_to_p_S33_reaclib = 87,
        k_p_P31_to_He4_Si28_reaclib = 88,
        k_He4_Cl34_to_p_Ar37_reaclib = 89,
        k_p_Cl35_to_He4_S32_reaclib = 90,
        k_He4_K38_to_p_Ca41_reaclib = 91,
        k_p_K39_to_He4_Ar36_reaclib = 92,
        k_He4_Sc42_to_p_Ti45_reaclib = 93,
        k_p_Sc43_to_He4_Ca40_reaclib = 94,
        k_He4_Ti44_to_p_V47_reaclib = 95,
        k_He4_V46_to_p_Cr49_reaclib = 96,
        k_He4_Cr48_to_p_Mn51_reaclib = 97,
        k_He4_Cr49_to_p_Mn52_reaclib = 98,
        k_He4_Fe52_to_p_Co55_reaclib = 99,
        k_He4_He4_He4_to_C12_reaclib = 100,
        k_n_Mn50_to_Mn51_reaclib = 101,
        k_n_Mn51_to_Mn52_reaclib = 102,
        k_p_Mn52_to_Fe53_reaclib = 103,
        k_He4_Mn52_to_Co56_reaclib = 104,
        k_n_Fe52_to_Fe53_reaclib = 105,
        k_n_Fe53_to_Fe54_reaclib = 106,
        k_He4_Fe53_to_Ni57_reaclib = 107,
        k_n_Fe54_to_Fe55_reaclib = 108,
        k_p_Fe54_to_Co55_reaclib = 109,
        k_He4_Fe54_to_Ni58_reaclib = 110,
        k_n_Fe55_to_Fe56_reaclib = 111,
        k_p_Fe55_to_Co56_reaclib = 112,
        k_p_Fe56_to_Co57_reaclib = 113,
        k_n_Co55_to_Co56_reaclib = 114,
        k_He4_Co55_to_Cu59_reaclib = 115,
        k_n_Co56_to_Co57_reaclib = 116,
        k_p_Co56_to_Ni57_reaclib = 117,
        k_p_Co57_to_Ni58_reaclib = 118,
        k_n_Ni56_to_Ni57_reaclib = 119,
        k_He4_Ni56_to_Zn60_reaclib = 120,
        k_n_Ni57_to_Ni58_reaclib = 121,
        k_p_Ni58_to_Cu59_reaclib = 122,
        k_p_Cu59_to_Zn60_reaclib = 123,
        k_He4_Mn50_to_p_Fe53_reaclib = 124,
        k_He4_Mn51_to_p_Fe54_reaclib = 125,
        k_He4_Mn52_to_p_Fe55_reaclib = 126,
        k_n_Fe52_to_p_Mn52_reaclib = 127,
        k_He4_Fe53_to_p_Co56_reaclib = 128,
        k_n_Co55_to_p_Fe55_reaclib = 129,
        k_n_Co55_to_He4_Mn52_reaclib = 130,
        k_He4_Co55_to_p_Ni58_reaclib = 131,
        k_n_Co56_to_p_Fe56_reaclib = 132,
        k_p_Co57_to_He4_Fe54_reaclib = 133,
        k_n_Ni56_to_p_Co56_reaclib = 134,
        k_n_Ni56_to_He4_Fe53_reaclib = 135,
        k_n_Ni57_to_p_Co57_reaclib = 136,
        k_n_Ni57_to_He4_Fe54_reaclib = 137,
        k_n_Ni58_to_He4_Fe55_reaclib = 138,
        k_n_Cu59_to_He4_Co56_reaclib = 139,
        k_p_Cu59_to_He4_Ni56_reaclib = 140,
        k_n_Zn60_to_He4_Ni57_reaclib = 141,
        k_Fe52_to_Mn52_weaktab = 142,
        k_Mn52_to_Fe52_weaktab = 143,
        k_Co55_to_Fe55_weaktab = 144,
        k_Fe55_to_Co55_weaktab = 145,
        k_Co56_to_Fe56_weaktab = 146,
        k_Co56_to_Ni56_weaktab = 147,
        k_Fe56_to_Co56_weaktab = 148,
        k_Ni56_to_Co56_weaktab = 149,
        k_Co57_to_Ni57_weaktab = 150,
        k_Ni57_to_Co57_weaktab = 151,
        k_n_to_p_weaktab = 152,
        k_p_to_n_weaktab = 153,
        k_F18_to_O18_weaktab = 154,
        k_O18_to_F18_weaktab = 155,
        k_Na22_to_Ne22_weaktab = 156,
        k_Ne22_to_Na22_weaktab = 157,
        k_N13_to_p_C12_derived = 158,
        k_N14_to_p_C13_derived = 159,
        k_N15_to_p_C14_derived = 160,
        k_O16_to_p_N15_derived = 161,
        k_O16_to_He4_C12_derived = 162,
        k_O18_to_He4_C14_derived = 163,
        k_F18_to_p_O17_derived = 164,
        k_F18_to_He4_N14_derived = 165,
        k_F19_to_p_O18_derived = 166,
        k_F19_to_He4_N15_derived = 167,
        k_Ne20_to_p_F19_derived = 168,
        k_Ne20_to_He4_O16_derived = 169,
        k_Ne21_to_He4_O17_derived = 170,
        k_Ne22_to_He4_O18_derived = 171,
        k_Na22_to_p_Ne21_derived = 172,
        k_Na22_to_He4_F18_derived = 173,
        k_Na23_to_p_Ne22_derived = 174,
        k_Na23_to_He4_F19_derived = 175,
        k_Mg24_to_p_Na23_derived = 176,
        k_Mg24_to_He4_Ne20_derived = 177,
        k_Mg25_to_He4_Ne21_derived = 178,
        k_Al26_to_p_Mg25_derived = 179,
        k_Al26_to_He4_Na22_derived = 180,
        k_Al27_to_He4_Na23_derived = 181,
        k_Si28_to_p_Al27_derived = 182,
        k_Si28_to_He4_Mg24_derived = 183,
        k_Si29_to_He4_Mg25_derived = 184,
        k_P30_to_p_Si29_derived = 185,
        k_P30_to_He4_Al26_derived = 186,
        k_P31_to_He4_Al27_derived = 187,
        k_S32_to_p_P31_derived = 188,
        k_S32_to_He4_Si28_derived = 189,
        k_S33_to_He4_Si29_derived = 190,
        k_Cl34_to_p_S33_derived = 191,
        k_Cl34_to_He4_P30_derived = 192,
        k_Cl35_to_He4_P31_derived = 193,
        k_Ar36_to_p_Cl35_derived = 194,
        k_Ar36_to_He4_S32_derived = 195,
        k_Ar37_to_He4_S33_derived = 196,
        k_K38_to_p_Ar37_derived = 197,
        k_K38_to_He4_Cl34_derived = 198,
        k_K39_to_He4_Cl35_derived = 199,
        k_Ca40_to_p_K39_derived = 200,
        k_Ca40_to_He4_Ar36_derived = 201,
        k_Ca41_to_He4_Ar37_derived = 202,
        k_Sc42_to_p_Ca41_derived = 203,
        k_Sc42_to_He4_K38_derived = 204,
        k_Sc43_to_He4_K39_derived = 205,
        k_Ti44_to_p_Sc43_derived = 206,
        k_Ti44_to_He4_Ca40_derived = 207,
        k_Ti45_to_He4_Ca41_derived = 208,
        k_V46_to_p_Ti45_derived = 209,
        k_V46_to_He4_Sc42_derived = 210,
        k_V47_to_He4_Sc43_derived = 211,
        k_Cr48_to_p_V47_derived = 212,
        k_Cr48_to_He4_Ti44_derived = 213,
        k_Cr49_to_He4_Ti45_derived = 214,
        k_Mn50_to_p_Cr49_derived = 215,
        k_Mn50_to_He4_V46_derived = 216,
        k_Mn51_to_He4_V47_derived = 217,
        k_Fe52_to_p_Mn51_derived = 218,
        k_Fe52_to_He4_Cr48_derived = 219,
        k_Co55_to_He4_Mn51_derived = 220,
        k_Ni56_to_p_Co55_derived = 221,
        k_Ni56_to_He4_Fe52_derived = 222,
        k_C12_to_He4_He4_He4_derived = 223,
        k_He4_C12_to_p_N15_derived = 224,
        k_He4_N14_to_p_O17_derived = 225,
        k_He4_N15_to_p_O18_derived = 226,
        k_p_O16_to_He4_N13_derived = 227,
        k_He4_O16_to_p_F19_derived = 228,
        k_He4_Ne20_to_p_Na23_derived = 229,
        k_He4_Ne20_to_C12_C12_derived = 230,
        k_p_Ne21_to_He4_F18_derived = 231,
        k_p_Ne22_to_He4_F19_derived = 232,
        k_p_Na23_to_C12_C12_derived = 233,
        k_He4_Mg24_to_p_Al27_derived = 234,
        k_He4_Mg24_to_C12_O16_derived = 235,
        k_p_Mg25_to_He4_Na22_derived = 236,
        k_p_Al27_to_C12_O16_derived = 237,
        k_He4_Si28_to_p_P31_derived = 238,
        k_He4_Si28_to_C12_Ne20_derived = 239,
        k_He4_Si28_to_O16_O16_derived = 240,
        k_p_Si29_to_He4_Al26_derived = 241,
        k_p_P31_to_C12_Ne20_derived = 242,
        k_p_P31_to_O16_O16_derived = 243,
        k_He4_S32_to_p_Cl35_derived = 244,
        k_p_S33_to_He4_P30_derived = 245,
        k_He4_Ar36_to_p_K39_derived = 246,
        k_p_Ar37_to_He4_Cl34_derived = 247,
        k_He4_Ca40_to_p_Sc43_derived = 248,
        k_p_Ca41_to_He4_K38_derived = 249,
        k_p_Ti45_to_He4_Sc42_derived = 250,
        k_p_V47_to_He4_Ti44_derived = 251,
        k_p_Cr49_to_He4_V46_derived = 252,
        k_p_Mn51_to_He4_Cr48_derived = 253,
        k_p_Mn52_to_He4_Cr49_derived = 254,
        k_p_Co55_to_He4_Fe52_derived = 255,
        k_Mn51_to_n_Mn50_derived = 256,
        k_Mn52_to_n_Mn51_derived = 257,
        k_Fe53_to_n_Fe52_derived = 258,
        k_Fe53_to_p_Mn52_derived = 259,
        k_Fe54_to_n_Fe53_derived = 260,
        k_Fe55_to_n_Fe54_derived = 261,
        k_Fe56_to_n_Fe55_derived = 262,
        k_Co55_to_p_Fe54_derived = 263,
        k_Co56_to_n_Co55_derived = 264,
        k_Co56_to_p_Fe55_derived = 265,
        k_Co56_to_He4_Mn52_derived = 266,
        k_Co57_to_n_Co56_derived = 267,
        k_Co57_to_p_Fe56_derived = 268,
        k_Ni57_to_n_Ni56_derived = 269,
        k_Ni57_to_p_Co56_derived = 270,
        k_Ni57_to_He4_Fe53_derived = 271,
        k_Ni58_to_n_Ni57_derived = 272,
        k_Ni58_to_p_Co57_derived = 273,
        k_Ni58_to_He4_Fe54_derived = 274,
        k_Cu59_to_p_Ni58_derived = 275,
        k_Cu59_to_He4_Co55_derived = 276,
        k_Zn60_to_p_Cu59_derived = 277,
        k_Zn60_to_He4_Ni56_derived = 278,
        k_p_Mn52_to_n_Fe52_derived = 279,
        k_He4_Mn52_to_n_Co55_derived = 280,
        k_p_Fe53_to_He4_Mn50_derived = 281,
        k_He4_Fe53_to_n_Ni56_derived = 282,
        k_p_Fe54_to_He4_Mn51_derived = 283,
        k_He4_Fe54_to_n_Ni57_derived = 284,
        k_He4_Fe54_to_p_Co57_derived = 285,
        k_p_Fe55_to_n_Co55_derived = 286,
        k_p_Fe55_to_He4_Mn52_derived = 287,
        k_He4_Fe55_to_n_Ni58_derived = 288,
        k_p_Fe56_to_n_Co56_derived = 289,
        k_p_Co56_to_n_Ni56_derived = 290,
        k_p_Co56_to_He4_Fe53_derived = 291,
        k_He4_Co56_to_n_Cu59_derived = 292,
        k_p_Co57_to_n_Ni57_derived = 293,
        k_He4_Ni56_to_p_Cu59_derived = 294,
        k_He4_Ni57_to_n_Zn60_derived = 295,
        k_p_Ni58_to_He4_Co55_derived = 296,
        NumRates = k_p_Ni58_to_He4_Co55_derived
    };

    enum ScreenPairs :
    std::uint8_t
    {
        k_p_Sc43 = 1,
        k_p_Fe56 = 2,
        k_He4_Ti44 = 3,
        k_He4_Fe55 = 4,
        k_He4_Mg24 = 5,
        k_C12_O16 = 6,
        k_p_Cr49 = 7,
        k_He4_Na22 = 8,
        k_p_V47 = 9,
        k_p_Co55 = 10,
        k_C12_C12 = 11,
        k_p_Fe53 = 12,
        k_He4_Cr49 = 13,
        k_p_Al27 = 14,
        k_He4_V47 = 15,
        k_p_Mg25 = 16,
        k_p_Ne22 = 17,
        k_He4_O18 = 18,
        k_p_Mn52 = 19,
        k_He4_Ni56 = 20,
        k_p_P31 = 21,
        k_p_Si29 = 22,
        k_p_Ni58 = 23,
        k_p_O16 = 24,
        k_He4_Mn52 = 25,
        k_O16_O16 = 26,
        k_He4_P31 = 27,
        k_He4_N13 = 28,
        k_He4_Si29 = 29,
        k_p_Ar37 = 30,
        k_p_F19 = 31,
        k_He4_Ca40 = 32,
        k_He4_K38 = 33,
        k_He4_Ar36 = 34,
        k_He4_F18 = 35,
        k_He4_Sc43 = 36,
        k_p_Ne21 = 37,
        k_He4_Ne20 = 38,
        k_p_Co56 = 39,
        k_He4_Mn50 = 40,
        k_p_Fe54 = 41,
        k_He4_Co55 = 42,
        k_He4_Fe53 = 43,
        k_He4_Al27 = 44,
        k_He4_Mg25 = 45,
        k_p_Cl35 = 46,
        k_p_S33 = 47,
        k_He4_Cl34 = 48,
        k_p_C12 = 49,
        k_He4_O16 = 50,
        k_He4_N14 = 51,
        k_p_Ti45 = 52,
        k_He4_Sc42 = 53,
        k_p_O17 = 54,
        k_p_C14 = 55,
        k_He4_Ar37 = 56,
        k_p_N15 = 57,
        k_He4_F19 = 58,
        k_He4_C14 = 59,
        k_He4_Al26 = 60,
        k_He4_Ne21 = 61,
        k_He4_S32 = 62,
        k_p_Cu59 = 63,
        k_He4_P30 = 64,
        k_He4_Co56 = 65,
        k_He4_Si28 = 66,
        k_He4_Fe54 = 67,
        k_He4_Ni57 = 68,
        k_He4_Be8 = 69,
        k_p_Ca41 = 70,
        k_p_Co57 = 71,
        k_p_K39 = 72,
        k_p_Fe55 = 73,
        k_He4_Cl35 = 74,
        k_p_C13 = 75,
        k_He4_S33 = 76,
        k_He4_C12 = 77,
        k_He4_Ca41 = 78,
        k_He4_K39 = 79,
        k_p_O18 = 80,
        k_He4_Cr48 = 81,
        k_He4_Ti45 = 82,
        k_He4_V46 = 83,
        k_He4_O17 = 84,
        k_He4_N15 = 85,
        k_p_Mn51 = 86,
        k_p_Na23 = 87,
        k_He4_Fe52 = 88,
        k_He4_Mn51 = 89,
        k_He4_He4 = 90,
        k_He4_Na23 = 91,
        k_C12_Ne20 = 92,
        NumScreenPairs = k_C12_Ne20
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "C14_to_N14_reaclib",  // 1,
        "N13_to_C13_reaclib",  // 2,
        "p_C12_to_N13_reaclib",  // 3,
        "He4_C12_to_O16_reaclib",  // 4,
        "p_C13_to_N14_reaclib",  // 5,
        "p_C14_to_N15_reaclib",  // 6,
        "He4_C14_to_O18_reaclib",  // 7,
        "He4_N14_to_F18_reaclib",  // 8,
        "p_N15_to_O16_reaclib",  // 9,
        "He4_N15_to_F19_reaclib",  // 10,
        "He4_O16_to_Ne20_reaclib",  // 11,
        "p_O17_to_F18_reaclib",  // 12,
        "He4_O17_to_Ne21_reaclib",  // 13,
        "p_O18_to_F19_reaclib",  // 14,
        "He4_O18_to_Ne22_reaclib",  // 15,
        "He4_F18_to_Na22_reaclib",  // 16,
        "p_F19_to_Ne20_reaclib",  // 17,
        "He4_F19_to_Na23_reaclib",  // 18,
        "He4_Ne20_to_Mg24_reaclib",  // 19,
        "p_Ne21_to_Na22_reaclib",  // 20,
        "He4_Ne21_to_Mg25_reaclib",  // 21,
        "p_Ne22_to_Na23_reaclib",  // 22,
        "He4_Na22_to_Al26_reaclib",  // 23,
        "p_Na23_to_Mg24_reaclib",  // 24,
        "He4_Na23_to_Al27_reaclib",  // 25,
        "He4_Mg24_to_Si28_reaclib",  // 26,
        "p_Mg25_to_Al26_reaclib",  // 27,
        "He4_Mg25_to_Si29_reaclib",  // 28,
        "He4_Al26_to_P30_reaclib",  // 29,
        "p_Al27_to_Si28_reaclib",  // 30,
        "He4_Al27_to_P31_reaclib",  // 31,
        "He4_Si28_to_S32_reaclib",  // 32,
        "p_Si29_to_P30_reaclib",  // 33,
        "He4_Si29_to_S33_reaclib",  // 34,
        "He4_P30_to_Cl34_reaclib",  // 35,
        "p_P31_to_S32_reaclib",  // 36,
        "He4_P31_to_Cl35_reaclib",  // 37,
        "He4_S32_to_Ar36_reaclib",  // 38,
        "p_S33_to_Cl34_reaclib",  // 39,
        "He4_S33_to_Ar37_reaclib",  // 40,
        "He4_Cl34_to_K38_reaclib",  // 41,
        "p_Cl35_to_Ar36_reaclib",  // 42,
        "He4_Cl35_to_K39_reaclib",  // 43,
        "He4_Ar36_to_Ca40_reaclib",  // 44,
        "p_Ar37_to_K38_reaclib",  // 45,
        "He4_Ar37_to_Ca41_reaclib",  // 46,
        "He4_K38_to_Sc42_reaclib",  // 47,
        "p_K39_to_Ca40_reaclib",  // 48,
        "He4_K39_to_Sc43_reaclib",  // 49,
        "He4_Ca40_to_Ti44_reaclib",  // 50,
        "p_Ca41_to_Sc42_reaclib",  // 51,
        "He4_Ca41_to_Ti45_reaclib",  // 52,
        "He4_Sc42_to_V46_reaclib",  // 53,
        "p_Sc43_to_Ti44_reaclib",  // 54,
        "He4_Sc43_to_V47_reaclib",  // 55,
        "He4_Ti44_to_Cr48_reaclib",  // 56,
        "p_Ti45_to_V46_reaclib",  // 57,
        "He4_Ti45_to_Cr49_reaclib",  // 58,
        "He4_V46_to_Mn50_reaclib",  // 59,
        "p_V47_to_Cr48_reaclib",  // 60,
        "He4_V47_to_Mn51_reaclib",  // 61,
        "He4_Cr48_to_Fe52_reaclib",  // 62,
        "p_Cr49_to_Mn50_reaclib",  // 63,
        "p_Mn51_to_Fe52_reaclib",  // 64,
        "He4_Mn51_to_Co55_reaclib",  // 65,
        "He4_Fe52_to_Ni56_reaclib",  // 66,
        "p_Co55_to_Ni56_reaclib",  // 67,
        "C12_C12_to_p_Na23_reaclib",  // 68,
        "C12_C12_to_He4_Ne20_reaclib",  // 69,
        "He4_N13_to_p_O16_reaclib",  // 70,
        "p_N15_to_He4_C12_reaclib",  // 71,
        "C12_O16_to_p_Al27_reaclib",  // 72,
        "C12_O16_to_He4_Mg24_reaclib",  // 73,
        "O16_O16_to_p_P31_reaclib",  // 74,
        "O16_O16_to_He4_Si28_reaclib",  // 75,
        "p_O17_to_He4_N14_reaclib",  // 76,
        "p_O18_to_He4_N15_reaclib",  // 77,
        "He4_F18_to_p_Ne21_reaclib",  // 78,
        "p_F19_to_He4_O16_reaclib",  // 79,
        "He4_F19_to_p_Ne22_reaclib",  // 80,
        "C12_Ne20_to_p_P31_reaclib",  // 81,
        "C12_Ne20_to_He4_Si28_reaclib",  // 82,
        "He4_Na22_to_p_Mg25_reaclib",  // 83,
        "p_Na23_to_He4_Ne20_reaclib",  // 84,
        "He4_Al26_to_p_Si29_reaclib",  // 85,
        "p_Al27_to_He4_Mg24_reaclib",  // 86,
        "He4_P30_to_p_S33_reaclib",  // 87,
        "p_P31_to_He4_Si28_reaclib",  // 88,
        "He4_Cl34_to_p_Ar37_reaclib",  // 89,
        "p_Cl35_to_He4_S32_reaclib",  // 90,
        "He4_K38_to_p_Ca41_reaclib",  // 91,
        "p_K39_to_He4_Ar36_reaclib",  // 92,
        "He4_Sc42_to_p_Ti45_reaclib",  // 93,
        "p_Sc43_to_He4_Ca40_reaclib",  // 94,
        "He4_Ti44_to_p_V47_reaclib",  // 95,
        "He4_V46_to_p_Cr49_reaclib",  // 96,
        "He4_Cr48_to_p_Mn51_reaclib",  // 97,
        "He4_Cr49_to_p_Mn52_reaclib",  // 98,
        "He4_Fe52_to_p_Co55_reaclib",  // 99,
        "He4_He4_He4_to_C12_reaclib",  // 100,
        "n_Mn50_to_Mn51_reaclib",  // 101,
        "n_Mn51_to_Mn52_reaclib",  // 102,
        "p_Mn52_to_Fe53_reaclib",  // 103,
        "He4_Mn52_to_Co56_reaclib",  // 104,
        "n_Fe52_to_Fe53_reaclib",  // 105,
        "n_Fe53_to_Fe54_reaclib",  // 106,
        "He4_Fe53_to_Ni57_reaclib",  // 107,
        "n_Fe54_to_Fe55_reaclib",  // 108,
        "p_Fe54_to_Co55_reaclib",  // 109,
        "He4_Fe54_to_Ni58_reaclib",  // 110,
        "n_Fe55_to_Fe56_reaclib",  // 111,
        "p_Fe55_to_Co56_reaclib",  // 112,
        "p_Fe56_to_Co57_reaclib",  // 113,
        "n_Co55_to_Co56_reaclib",  // 114,
        "He4_Co55_to_Cu59_reaclib",  // 115,
        "n_Co56_to_Co57_reaclib",  // 116,
        "p_Co56_to_Ni57_reaclib",  // 117,
        "p_Co57_to_Ni58_reaclib",  // 118,
        "n_Ni56_to_Ni57_reaclib",  // 119,
        "He4_Ni56_to_Zn60_reaclib",  // 120,
        "n_Ni57_to_Ni58_reaclib",  // 121,
        "p_Ni58_to_Cu59_reaclib",  // 122,
        "p_Cu59_to_Zn60_reaclib",  // 123,
        "He4_Mn50_to_p_Fe53_reaclib",  // 124,
        "He4_Mn51_to_p_Fe54_reaclib",  // 125,
        "He4_Mn52_to_p_Fe55_reaclib",  // 126,
        "n_Fe52_to_p_Mn52_reaclib",  // 127,
        "He4_Fe53_to_p_Co56_reaclib",  // 128,
        "n_Co55_to_p_Fe55_reaclib",  // 129,
        "n_Co55_to_He4_Mn52_reaclib",  // 130,
        "He4_Co55_to_p_Ni58_reaclib",  // 131,
        "n_Co56_to_p_Fe56_reaclib",  // 132,
        "p_Co57_to_He4_Fe54_reaclib",  // 133,
        "n_Ni56_to_p_Co56_reaclib",  // 134,
        "n_Ni56_to_He4_Fe53_reaclib",  // 135,
        "n_Ni57_to_p_Co57_reaclib",  // 136,
        "n_Ni57_to_He4_Fe54_reaclib",  // 137,
        "n_Ni58_to_He4_Fe55_reaclib",  // 138,
        "n_Cu59_to_He4_Co56_reaclib",  // 139,
        "p_Cu59_to_He4_Ni56_reaclib",  // 140,
        "n_Zn60_to_He4_Ni57_reaclib",  // 141,
        "Fe52_to_Mn52_weaktab",  // 142,
        "Mn52_to_Fe52_weaktab",  // 143,
        "Co55_to_Fe55_weaktab",  // 144,
        "Fe55_to_Co55_weaktab",  // 145,
        "Co56_to_Fe56_weaktab",  // 146,
        "Co56_to_Ni56_weaktab",  // 147,
        "Fe56_to_Co56_weaktab",  // 148,
        "Ni56_to_Co56_weaktab",  // 149,
        "Co57_to_Ni57_weaktab",  // 150,
        "Ni57_to_Co57_weaktab",  // 151,
        "n_to_p_weaktab",  // 152,
        "p_to_n_weaktab",  // 153,
        "F18_to_O18_weaktab",  // 154,
        "O18_to_F18_weaktab",  // 155,
        "Na22_to_Ne22_weaktab",  // 156,
        "Ne22_to_Na22_weaktab",  // 157,
        "N13_to_p_C12_derived",  // 158,
        "N14_to_p_C13_derived",  // 159,
        "N15_to_p_C14_derived",  // 160,
        "O16_to_p_N15_derived",  // 161,
        "O16_to_He4_C12_derived",  // 162,
        "O18_to_He4_C14_derived",  // 163,
        "F18_to_p_O17_derived",  // 164,
        "F18_to_He4_N14_derived",  // 165,
        "F19_to_p_O18_derived",  // 166,
        "F19_to_He4_N15_derived",  // 167,
        "Ne20_to_p_F19_derived",  // 168,
        "Ne20_to_He4_O16_derived",  // 169,
        "Ne21_to_He4_O17_derived",  // 170,
        "Ne22_to_He4_O18_derived",  // 171,
        "Na22_to_p_Ne21_derived",  // 172,
        "Na22_to_He4_F18_derived",  // 173,
        "Na23_to_p_Ne22_derived",  // 174,
        "Na23_to_He4_F19_derived",  // 175,
        "Mg24_to_p_Na23_derived",  // 176,
        "Mg24_to_He4_Ne20_derived",  // 177,
        "Mg25_to_He4_Ne21_derived",  // 178,
        "Al26_to_p_Mg25_derived",  // 179,
        "Al26_to_He4_Na22_derived",  // 180,
        "Al27_to_He4_Na23_derived",  // 181,
        "Si28_to_p_Al27_derived",  // 182,
        "Si28_to_He4_Mg24_derived",  // 183,
        "Si29_to_He4_Mg25_derived",  // 184,
        "P30_to_p_Si29_derived",  // 185,
        "P30_to_He4_Al26_derived",  // 186,
        "P31_to_He4_Al27_derived",  // 187,
        "S32_to_p_P31_derived",  // 188,
        "S32_to_He4_Si28_derived",  // 189,
        "S33_to_He4_Si29_derived",  // 190,
        "Cl34_to_p_S33_derived",  // 191,
        "Cl34_to_He4_P30_derived",  // 192,
        "Cl35_to_He4_P31_derived",  // 193,
        "Ar36_to_p_Cl35_derived",  // 194,
        "Ar36_to_He4_S32_derived",  // 195,
        "Ar37_to_He4_S33_derived",  // 196,
        "K38_to_p_Ar37_derived",  // 197,
        "K38_to_He4_Cl34_derived",  // 198,
        "K39_to_He4_Cl35_derived",  // 199,
        "Ca40_to_p_K39_derived",  // 200,
        "Ca40_to_He4_Ar36_derived",  // 201,
        "Ca41_to_He4_Ar37_derived",  // 202,
        "Sc42_to_p_Ca41_derived",  // 203,
        "Sc42_to_He4_K38_derived",  // 204,
        "Sc43_to_He4_K39_derived",  // 205,
        "Ti44_to_p_Sc43_derived",  // 206,
        "Ti44_to_He4_Ca40_derived",  // 207,
        "Ti45_to_He4_Ca41_derived",  // 208,
        "V46_to_p_Ti45_derived",  // 209,
        "V46_to_He4_Sc42_derived",  // 210,
        "V47_to_He4_Sc43_derived",  // 211,
        "Cr48_to_p_V47_derived",  // 212,
        "Cr48_to_He4_Ti44_derived",  // 213,
        "Cr49_to_He4_Ti45_derived",  // 214,
        "Mn50_to_p_Cr49_derived",  // 215,
        "Mn50_to_He4_V46_derived",  // 216,
        "Mn51_to_He4_V47_derived",  // 217,
        "Fe52_to_p_Mn51_derived",  // 218,
        "Fe52_to_He4_Cr48_derived",  // 219,
        "Co55_to_He4_Mn51_derived",  // 220,
        "Ni56_to_p_Co55_derived",  // 221,
        "Ni56_to_He4_Fe52_derived",  // 222,
        "C12_to_He4_He4_He4_derived",  // 223,
        "He4_C12_to_p_N15_derived",  // 224,
        "He4_N14_to_p_O17_derived",  // 225,
        "He4_N15_to_p_O18_derived",  // 226,
        "p_O16_to_He4_N13_derived",  // 227,
        "He4_O16_to_p_F19_derived",  // 228,
        "He4_Ne20_to_p_Na23_derived",  // 229,
        "He4_Ne20_to_C12_C12_derived",  // 230,
        "p_Ne21_to_He4_F18_derived",  // 231,
        "p_Ne22_to_He4_F19_derived",  // 232,
        "p_Na23_to_C12_C12_derived",  // 233,
        "He4_Mg24_to_p_Al27_derived",  // 234,
        "He4_Mg24_to_C12_O16_derived",  // 235,
        "p_Mg25_to_He4_Na22_derived",  // 236,
        "p_Al27_to_C12_O16_derived",  // 237,
        "He4_Si28_to_p_P31_derived",  // 238,
        "He4_Si28_to_C12_Ne20_derived",  // 239,
        "He4_Si28_to_O16_O16_derived",  // 240,
        "p_Si29_to_He4_Al26_derived",  // 241,
        "p_P31_to_C12_Ne20_derived",  // 242,
        "p_P31_to_O16_O16_derived",  // 243,
        "He4_S32_to_p_Cl35_derived",  // 244,
        "p_S33_to_He4_P30_derived",  // 245,
        "He4_Ar36_to_p_K39_derived",  // 246,
        "p_Ar37_to_He4_Cl34_derived",  // 247,
        "He4_Ca40_to_p_Sc43_derived",  // 248,
        "p_Ca41_to_He4_K38_derived",  // 249,
        "p_Ti45_to_He4_Sc42_derived",  // 250,
        "p_V47_to_He4_Ti44_derived",  // 251,
        "p_Cr49_to_He4_V46_derived",  // 252,
        "p_Mn51_to_He4_Cr48_derived",  // 253,
        "p_Mn52_to_He4_Cr49_derived",  // 254,
        "p_Co55_to_He4_Fe52_derived",  // 255,
        "Mn51_to_n_Mn50_derived",  // 256,
        "Mn52_to_n_Mn51_derived",  // 257,
        "Fe53_to_n_Fe52_derived",  // 258,
        "Fe53_to_p_Mn52_derived",  // 259,
        "Fe54_to_n_Fe53_derived",  // 260,
        "Fe55_to_n_Fe54_derived",  // 261,
        "Fe56_to_n_Fe55_derived",  // 262,
        "Co55_to_p_Fe54_derived",  // 263,
        "Co56_to_n_Co55_derived",  // 264,
        "Co56_to_p_Fe55_derived",  // 265,
        "Co56_to_He4_Mn52_derived",  // 266,
        "Co57_to_n_Co56_derived",  // 267,
        "Co57_to_p_Fe56_derived",  // 268,
        "Ni57_to_n_Ni56_derived",  // 269,
        "Ni57_to_p_Co56_derived",  // 270,
        "Ni57_to_He4_Fe53_derived",  // 271,
        "Ni58_to_n_Ni57_derived",  // 272,
        "Ni58_to_p_Co57_derived",  // 273,
        "Ni58_to_He4_Fe54_derived",  // 274,
        "Cu59_to_p_Ni58_derived",  // 275,
        "Cu59_to_He4_Co55_derived",  // 276,
        "Zn60_to_p_Cu59_derived",  // 277,
        "Zn60_to_He4_Ni56_derived",  // 278,
        "p_Mn52_to_n_Fe52_derived",  // 279,
        "He4_Mn52_to_n_Co55_derived",  // 280,
        "p_Fe53_to_He4_Mn50_derived",  // 281,
        "He4_Fe53_to_n_Ni56_derived",  // 282,
        "p_Fe54_to_He4_Mn51_derived",  // 283,
        "He4_Fe54_to_n_Ni57_derived",  // 284,
        "He4_Fe54_to_p_Co57_derived",  // 285,
        "p_Fe55_to_n_Co55_derived",  // 286,
        "p_Fe55_to_He4_Mn52_derived",  // 287,
        "He4_Fe55_to_n_Ni58_derived",  // 288,
        "p_Fe56_to_n_Co56_derived",  // 289,
        "p_Co56_to_n_Ni56_derived",  // 290,
        "p_Co56_to_He4_Fe53_derived",  // 291,
        "He4_Co56_to_n_Cu59_derived",  // 292,
        "p_Co57_to_n_Ni57_derived",  // 293,
        "He4_Ni56_to_p_Cu59_derived",  // 294,
        "He4_Ni57_to_n_Zn60_derived",  // 295,
        "p_Ni58_to_He4_Co55_derived"  // 296,
    };

    template<std::uint8_t rate>
    AMREX_GPU_DEVICE AMREX_INLINE
    constexpr amrex::Real get_p_random() {
        return 0.0_rt;
    }

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

    constexpr int H1_index = 2;
    constexpr int N_index = 1;
    constexpr int He4_index = 3;
    constexpr int NumNSERatePairs = 131;

    inline AMREX_GPU_MANAGED amrex::Array2D<short, 1, NumNSERatePairs, 1, 10, amrex::Order::C> rate_pair_data {
        H1, C12, -1, N13, -1, -1, C12, N13, k_p_C12_to_N13_reaclib, k_N13_to_p_C12_derived,
        He4, C12, -1, O16, -1, -1, C12, O16, k_He4_C12_to_O16_reaclib, k_O16_to_He4_C12_derived,
        H1, C13, -1, N14, -1, -1, C13, N14, k_p_C13_to_N14_reaclib, k_N14_to_p_C13_derived,
        H1, C14, -1, N15, -1, -1, C14, N15, k_p_C14_to_N15_reaclib, k_N15_to_p_C14_derived,
        He4, C14, -1, O18, -1, -1, C14, O18, k_He4_C14_to_O18_reaclib, k_O18_to_He4_C14_derived,
        He4, N14, -1, F18, -1, -1, N14, F18, k_He4_N14_to_F18_reaclib, k_F18_to_He4_N14_derived,
        H1, N15, -1, O16, -1, -1, N15, O16, k_p_N15_to_O16_reaclib, k_O16_to_p_N15_derived,
        He4, N15, -1, F19, -1, -1, N15, F19, k_He4_N15_to_F19_reaclib, k_F19_to_He4_N15_derived,
        He4, O16, -1, Ne20, -1, -1, O16, Ne20, k_He4_O16_to_Ne20_reaclib, k_Ne20_to_He4_O16_derived,
        H1, O17, -1, F18, -1, -1, O17, F18, k_p_O17_to_F18_reaclib, k_F18_to_p_O17_derived,
        He4, O17, -1, Ne21, -1, -1, O17, Ne21, k_He4_O17_to_Ne21_reaclib, k_Ne21_to_He4_O17_derived,
        H1, O18, -1, F19, -1, -1, O18, F19, k_p_O18_to_F19_reaclib, k_F19_to_p_O18_derived,
        He4, O18, -1, Ne22, -1, -1, O18, Ne22, k_He4_O18_to_Ne22_reaclib, k_Ne22_to_He4_O18_derived,
        He4, F18, -1, Na22, -1, -1, F18, Na22, k_He4_F18_to_Na22_reaclib, k_Na22_to_He4_F18_derived,
        H1, F19, -1, Ne20, -1, -1, F19, Ne20, k_p_F19_to_Ne20_reaclib, k_Ne20_to_p_F19_derived,
        He4, F19, -1, Na23, -1, -1, F19, Na23, k_He4_F19_to_Na23_reaclib, k_Na23_to_He4_F19_derived,
        He4, Ne20, -1, Mg24, -1, -1, Ne20, Mg24, k_He4_Ne20_to_Mg24_reaclib, k_Mg24_to_He4_Ne20_derived,
        H1, Ne21, -1, Na22, -1, -1, Ne21, Na22, k_p_Ne21_to_Na22_reaclib, k_Na22_to_p_Ne21_derived,
        He4, Ne21, -1, Mg25, -1, -1, Ne21, Mg25, k_He4_Ne21_to_Mg25_reaclib, k_Mg25_to_He4_Ne21_derived,
        H1, Ne22, -1, Na23, -1, -1, Ne22, Na23, k_p_Ne22_to_Na23_reaclib, k_Na23_to_p_Ne22_derived,
        He4, Na22, -1, Al26, -1, -1, Na22, Al26, k_He4_Na22_to_Al26_reaclib, k_Al26_to_He4_Na22_derived,
        H1, Na23, -1, Mg24, -1, -1, Na23, Mg24, k_p_Na23_to_Mg24_reaclib, k_Mg24_to_p_Na23_derived,
        He4, Na23, -1, Al27, -1, -1, Na23, Al27, k_He4_Na23_to_Al27_reaclib, k_Al27_to_He4_Na23_derived,
        He4, Mg24, -1, Si28, -1, -1, Mg24, Si28, k_He4_Mg24_to_Si28_reaclib, k_Si28_to_He4_Mg24_derived,
        H1, Mg25, -1, Al26, -1, -1, Mg25, Al26, k_p_Mg25_to_Al26_reaclib, k_Al26_to_p_Mg25_derived,
        He4, Mg25, -1, Si29, -1, -1, Mg25, Si29, k_He4_Mg25_to_Si29_reaclib, k_Si29_to_He4_Mg25_derived,
        He4, Al26, -1, P30, -1, -1, Al26, P30, k_He4_Al26_to_P30_reaclib, k_P30_to_He4_Al26_derived,
        H1, Al27, -1, Si28, -1, -1, Al27, Si28, k_p_Al27_to_Si28_reaclib, k_Si28_to_p_Al27_derived,
        He4, Al27, -1, P31, -1, -1, Al27, P31, k_He4_Al27_to_P31_reaclib, k_P31_to_He4_Al27_derived,
        He4, Si28, -1, S32, -1, -1, Si28, S32, k_He4_Si28_to_S32_reaclib, k_S32_to_He4_Si28_derived,
        H1, Si29, -1, P30, -1, -1, Si29, P30, k_p_Si29_to_P30_reaclib, k_P30_to_p_Si29_derived,
        He4, Si29, -1, S33, -1, -1, Si29, S33, k_He4_Si29_to_S33_reaclib, k_S33_to_He4_Si29_derived,
        He4, P30, -1, Cl34, -1, -1, P30, Cl34, k_He4_P30_to_Cl34_reaclib, k_Cl34_to_He4_P30_derived,
        H1, P31, -1, S32, -1, -1, P31, S32, k_p_P31_to_S32_reaclib, k_S32_to_p_P31_derived,
        He4, P31, -1, Cl35, -1, -1, P31, Cl35, k_He4_P31_to_Cl35_reaclib, k_Cl35_to_He4_P31_derived,
        He4, S32, -1, Ar36, -1, -1, S32, Ar36, k_He4_S32_to_Ar36_reaclib, k_Ar36_to_He4_S32_derived,
        H1, S33, -1, Cl34, -1, -1, S33, Cl34, k_p_S33_to_Cl34_reaclib, k_Cl34_to_p_S33_derived,
        He4, S33, -1, Ar37, -1, -1, S33, Ar37, k_He4_S33_to_Ar37_reaclib, k_Ar37_to_He4_S33_derived,
        He4, Cl34, -1, K38, -1, -1, Cl34, K38, k_He4_Cl34_to_K38_reaclib, k_K38_to_He4_Cl34_derived,
        H1, Cl35, -1, Ar36, -1, -1, Cl35, Ar36, k_p_Cl35_to_Ar36_reaclib, k_Ar36_to_p_Cl35_derived,
        He4, Cl35, -1, K39, -1, -1, Cl35, K39, k_He4_Cl35_to_K39_reaclib, k_K39_to_He4_Cl35_derived,
        He4, Ar36, -1, Ca40, -1, -1, Ar36, Ca40, k_He4_Ar36_to_Ca40_reaclib, k_Ca40_to_He4_Ar36_derived,
        H1, Ar37, -1, K38, -1, -1, Ar37, K38, k_p_Ar37_to_K38_reaclib, k_K38_to_p_Ar37_derived,
        He4, Ar37, -1, Ca41, -1, -1, Ar37, Ca41, k_He4_Ar37_to_Ca41_reaclib, k_Ca41_to_He4_Ar37_derived,
        He4, K38, -1, Sc42, -1, -1, K38, Sc42, k_He4_K38_to_Sc42_reaclib, k_Sc42_to_He4_K38_derived,
        H1, K39, -1, Ca40, -1, -1, K39, Ca40, k_p_K39_to_Ca40_reaclib, k_Ca40_to_p_K39_derived,
        He4, K39, -1, Sc43, -1, -1, K39, Sc43, k_He4_K39_to_Sc43_reaclib, k_Sc43_to_He4_K39_derived,
        He4, Ca40, -1, Ti44, -1, -1, Ca40, Ti44, k_He4_Ca40_to_Ti44_reaclib, k_Ti44_to_He4_Ca40_derived,
        H1, Ca41, -1, Sc42, -1, -1, Ca41, Sc42, k_p_Ca41_to_Sc42_reaclib, k_Sc42_to_p_Ca41_derived,
        He4, Ca41, -1, Ti45, -1, -1, Ca41, Ti45, k_He4_Ca41_to_Ti45_reaclib, k_Ti45_to_He4_Ca41_derived,
        He4, Sc42, -1, V46, -1, -1, Sc42, V46, k_He4_Sc42_to_V46_reaclib, k_V46_to_He4_Sc42_derived,
        H1, Sc43, -1, Ti44, -1, -1, Sc43, Ti44, k_p_Sc43_to_Ti44_reaclib, k_Ti44_to_p_Sc43_derived,
        He4, Sc43, -1, V47, -1, -1, Sc43, V47, k_He4_Sc43_to_V47_reaclib, k_V47_to_He4_Sc43_derived,
        He4, Ti44, -1, Cr48, -1, -1, Ti44, Cr48, k_He4_Ti44_to_Cr48_reaclib, k_Cr48_to_He4_Ti44_derived,
        H1, Ti45, -1, V46, -1, -1, Ti45, V46, k_p_Ti45_to_V46_reaclib, k_V46_to_p_Ti45_derived,
        He4, Ti45, -1, Cr49, -1, -1, Ti45, Cr49, k_He4_Ti45_to_Cr49_reaclib, k_Cr49_to_He4_Ti45_derived,
        He4, V46, -1, Mn50, -1, -1, V46, Mn50, k_He4_V46_to_Mn50_reaclib, k_Mn50_to_He4_V46_derived,
        H1, V47, -1, Cr48, -1, -1, V47, Cr48, k_p_V47_to_Cr48_reaclib, k_Cr48_to_p_V47_derived,
        He4, V47, -1, Mn51, -1, -1, V47, Mn51, k_He4_V47_to_Mn51_reaclib, k_Mn51_to_He4_V47_derived,
        He4, Cr48, -1, Fe52, -1, -1, Cr48, Fe52, k_He4_Cr48_to_Fe52_reaclib, k_Fe52_to_He4_Cr48_derived,
        H1, Cr49, -1, Mn50, -1, -1, Cr49, Mn50, k_p_Cr49_to_Mn50_reaclib, k_Mn50_to_p_Cr49_derived,
        H1, Mn51, -1, Fe52, -1, -1, Mn51, Fe52, k_p_Mn51_to_Fe52_reaclib, k_Fe52_to_p_Mn51_derived,
        He4, Mn51, -1, Co55, -1, -1, Mn51, Co55, k_He4_Mn51_to_Co55_reaclib, k_Co55_to_He4_Mn51_derived,
        He4, Fe52, -1, Ni56, -1, -1, Fe52, Ni56, k_He4_Fe52_to_Ni56_reaclib, k_Ni56_to_He4_Fe52_derived,
        H1, Co55, -1, Ni56, -1, -1, Co55, Ni56, k_p_Co55_to_Ni56_reaclib, k_Ni56_to_p_Co55_derived,
        He4, N13, -1, H1, O16, -1, N13, O16, k_He4_N13_to_p_O16_reaclib, k_p_O16_to_He4_N13_derived,
        H1, N15, -1, He4, C12, -1, N15, C12, k_p_N15_to_He4_C12_reaclib, k_He4_C12_to_p_N15_derived,
        H1, O17, -1, He4, N14, -1, O17, N14, k_p_O17_to_He4_N14_reaclib, k_He4_N14_to_p_O17_derived,
        H1, O18, -1, He4, N15, -1, O18, N15, k_p_O18_to_He4_N15_reaclib, k_He4_N15_to_p_O18_derived,
        He4, F18, -1, H1, Ne21, -1, F18, Ne21, k_He4_F18_to_p_Ne21_reaclib, k_p_Ne21_to_He4_F18_derived,
        H1, F19, -1, He4, O16, -1, F19, O16, k_p_F19_to_He4_O16_reaclib, k_He4_O16_to_p_F19_derived,
        He4, F19, -1, H1, Ne22, -1, F19, Ne22, k_He4_F19_to_p_Ne22_reaclib, k_p_Ne22_to_He4_F19_derived,
        He4, Na22, -1, H1, Mg25, -1, Na22, Mg25, k_He4_Na22_to_p_Mg25_reaclib, k_p_Mg25_to_He4_Na22_derived,
        H1, Na23, -1, He4, Ne20, -1, Na23, Ne20, k_p_Na23_to_He4_Ne20_reaclib, k_He4_Ne20_to_p_Na23_derived,
        He4, Al26, -1, H1, Si29, -1, Al26, Si29, k_He4_Al26_to_p_Si29_reaclib, k_p_Si29_to_He4_Al26_derived,
        H1, Al27, -1, He4, Mg24, -1, Al27, Mg24, k_p_Al27_to_He4_Mg24_reaclib, k_He4_Mg24_to_p_Al27_derived,
        He4, P30, -1, H1, S33, -1, P30, S33, k_He4_P30_to_p_S33_reaclib, k_p_S33_to_He4_P30_derived,
        H1, P31, -1, He4, Si28, -1, P31, Si28, k_p_P31_to_He4_Si28_reaclib, k_He4_Si28_to_p_P31_derived,
        He4, Cl34, -1, H1, Ar37, -1, Cl34, Ar37, k_He4_Cl34_to_p_Ar37_reaclib, k_p_Ar37_to_He4_Cl34_derived,
        H1, Cl35, -1, He4, S32, -1, Cl35, S32, k_p_Cl35_to_He4_S32_reaclib, k_He4_S32_to_p_Cl35_derived,
        He4, K38, -1, H1, Ca41, -1, K38, Ca41, k_He4_K38_to_p_Ca41_reaclib, k_p_Ca41_to_He4_K38_derived,
        H1, K39, -1, He4, Ar36, -1, K39, Ar36, k_p_K39_to_He4_Ar36_reaclib, k_He4_Ar36_to_p_K39_derived,
        He4, Sc42, -1, H1, Ti45, -1, Sc42, Ti45, k_He4_Sc42_to_p_Ti45_reaclib, k_p_Ti45_to_He4_Sc42_derived,
        H1, Sc43, -1, He4, Ca40, -1, Sc43, Ca40, k_p_Sc43_to_He4_Ca40_reaclib, k_He4_Ca40_to_p_Sc43_derived,
        He4, V46, -1, H1, Cr49, -1, V46, Cr49, k_He4_V46_to_p_Cr49_reaclib, k_p_Cr49_to_He4_V46_derived,
        He4, Cr48, -1, H1, Mn51, -1, Cr48, Mn51, k_He4_Cr48_to_p_Mn51_reaclib, k_p_Mn51_to_He4_Cr48_derived,
        He4, Cr49, -1, H1, Mn52, -1, Cr49, Mn52, k_He4_Cr49_to_p_Mn52_reaclib, k_p_Mn52_to_He4_Cr49_derived,
        He4, Fe52, -1, H1, Co55, -1, Fe52, Co55, k_He4_Fe52_to_p_Co55_reaclib, k_p_Co55_to_He4_Fe52_derived,
        He4, He4, He4, C12, -1, -1, C12, -1, k_He4_He4_He4_to_C12_reaclib, k_C12_to_He4_He4_He4_derived,
        N, Mn50, -1, Mn51, -1, -1, Mn50, Mn51, k_n_Mn50_to_Mn51_reaclib, k_Mn51_to_n_Mn50_derived,
        N, Mn51, -1, Mn52, -1, -1, Mn51, Mn52, k_n_Mn51_to_Mn52_reaclib, k_Mn52_to_n_Mn51_derived,
        H1, Mn52, -1, Fe53, -1, -1, Mn52, Fe53, k_p_Mn52_to_Fe53_reaclib, k_Fe53_to_p_Mn52_derived,
        He4, Mn52, -1, Co56, -1, -1, Mn52, Co56, k_He4_Mn52_to_Co56_reaclib, k_Co56_to_He4_Mn52_derived,
        N, Fe52, -1, Fe53, -1, -1, Fe52, Fe53, k_n_Fe52_to_Fe53_reaclib, k_Fe53_to_n_Fe52_derived,
        N, Fe53, -1, Fe54, -1, -1, Fe53, Fe54, k_n_Fe53_to_Fe54_reaclib, k_Fe54_to_n_Fe53_derived,
        He4, Fe53, -1, Ni57, -1, -1, Fe53, Ni57, k_He4_Fe53_to_Ni57_reaclib, k_Ni57_to_He4_Fe53_derived,
        N, Fe54, -1, Fe55, -1, -1, Fe54, Fe55, k_n_Fe54_to_Fe55_reaclib, k_Fe55_to_n_Fe54_derived,
        H1, Fe54, -1, Co55, -1, -1, Fe54, Co55, k_p_Fe54_to_Co55_reaclib, k_Co55_to_p_Fe54_derived,
        He4, Fe54, -1, Ni58, -1, -1, Fe54, Ni58, k_He4_Fe54_to_Ni58_reaclib, k_Ni58_to_He4_Fe54_derived,
        N, Fe55, -1, Fe56, -1, -1, Fe55, Fe56, k_n_Fe55_to_Fe56_reaclib, k_Fe56_to_n_Fe55_derived,
        H1, Fe55, -1, Co56, -1, -1, Fe55, Co56, k_p_Fe55_to_Co56_reaclib, k_Co56_to_p_Fe55_derived,
        H1, Fe56, -1, Co57, -1, -1, Fe56, Co57, k_p_Fe56_to_Co57_reaclib, k_Co57_to_p_Fe56_derived,
        N, Co55, -1, Co56, -1, -1, Co55, Co56, k_n_Co55_to_Co56_reaclib, k_Co56_to_n_Co55_derived,
        He4, Co55, -1, Cu59, -1, -1, Co55, Cu59, k_He4_Co55_to_Cu59_reaclib, k_Cu59_to_He4_Co55_derived,
        N, Co56, -1, Co57, -1, -1, Co56, Co57, k_n_Co56_to_Co57_reaclib, k_Co57_to_n_Co56_derived,
        H1, Co56, -1, Ni57, -1, -1, Co56, Ni57, k_p_Co56_to_Ni57_reaclib, k_Ni57_to_p_Co56_derived,
        H1, Co57, -1, Ni58, -1, -1, Co57, Ni58, k_p_Co57_to_Ni58_reaclib, k_Ni58_to_p_Co57_derived,
        N, Ni56, -1, Ni57, -1, -1, Ni56, Ni57, k_n_Ni56_to_Ni57_reaclib, k_Ni57_to_n_Ni56_derived,
        He4, Ni56, -1, Zn60, -1, -1, Ni56, Zn60, k_He4_Ni56_to_Zn60_reaclib, k_Zn60_to_He4_Ni56_derived,
        N, Ni57, -1, Ni58, -1, -1, Ni57, Ni58, k_n_Ni57_to_Ni58_reaclib, k_Ni58_to_n_Ni57_derived,
        H1, Ni58, -1, Cu59, -1, -1, Ni58, Cu59, k_p_Ni58_to_Cu59_reaclib, k_Cu59_to_p_Ni58_derived,
        H1, Cu59, -1, Zn60, -1, -1, Cu59, Zn60, k_p_Cu59_to_Zn60_reaclib, k_Zn60_to_p_Cu59_derived,
        He4, Mn50, -1, H1, Fe53, -1, Mn50, Fe53, k_He4_Mn50_to_p_Fe53_reaclib, k_p_Fe53_to_He4_Mn50_derived,
        He4, Mn51, -1, H1, Fe54, -1, Mn51, Fe54, k_He4_Mn51_to_p_Fe54_reaclib, k_p_Fe54_to_He4_Mn51_derived,
        He4, Mn52, -1, H1, Fe55, -1, Mn52, Fe55, k_He4_Mn52_to_p_Fe55_reaclib, k_p_Fe55_to_He4_Mn52_derived,
        N, Fe52, -1, H1, Mn52, -1, Fe52, Mn52, k_n_Fe52_to_p_Mn52_reaclib, k_p_Mn52_to_n_Fe52_derived,
        He4, Fe53, -1, H1, Co56, -1, Fe53, Co56, k_He4_Fe53_to_p_Co56_reaclib, k_p_Co56_to_He4_Fe53_derived,
        N, Co55, -1, H1, Fe55, -1, Co55, Fe55, k_n_Co55_to_p_Fe55_reaclib, k_p_Fe55_to_n_Co55_derived,
        N, Co55, -1, He4, Mn52, -1, Co55, Mn52, k_n_Co55_to_He4_Mn52_reaclib, k_He4_Mn52_to_n_Co55_derived,
        He4, Co55, -1, H1, Ni58, -1, Co55, Ni58, k_He4_Co55_to_p_Ni58_reaclib, k_p_Ni58_to_He4_Co55_derived,
        N, Co56, -1, H1, Fe56, -1, Co56, Fe56, k_n_Co56_to_p_Fe56_reaclib, k_p_Fe56_to_n_Co56_derived,
        H1, Co57, -1, He4, Fe54, -1, Co57, Fe54, k_p_Co57_to_He4_Fe54_reaclib, k_He4_Fe54_to_p_Co57_derived,
        N, Ni56, -1, H1, Co56, -1, Ni56, Co56, k_n_Ni56_to_p_Co56_reaclib, k_p_Co56_to_n_Ni56_derived,
        N, Ni56, -1, He4, Fe53, -1, Ni56, Fe53, k_n_Ni56_to_He4_Fe53_reaclib, k_He4_Fe53_to_n_Ni56_derived,
        N, Ni57, -1, H1, Co57, -1, Ni57, Co57, k_n_Ni57_to_p_Co57_reaclib, k_p_Co57_to_n_Ni57_derived,
        N, Ni57, -1, He4, Fe54, -1, Ni57, Fe54, k_n_Ni57_to_He4_Fe54_reaclib, k_He4_Fe54_to_n_Ni57_derived,
        N, Ni58, -1, He4, Fe55, -1, Ni58, Fe55, k_n_Ni58_to_He4_Fe55_reaclib, k_He4_Fe55_to_n_Ni58_derived,
        N, Cu59, -1, He4, Co56, -1, Cu59, Co56, k_n_Cu59_to_He4_Co56_reaclib, k_He4_Co56_to_n_Cu59_derived,
        H1, Cu59, -1, He4, Ni56, -1, Cu59, Ni56, k_p_Cu59_to_He4_Ni56_reaclib, k_He4_Ni56_to_p_Cu59_derived,
        N, Zn60, -1, He4, Ni57, -1, Zn60, Ni57, k_n_Zn60_to_He4_Ni57_reaclib, k_He4_Ni57_to_n_Zn60_derived,
        H1, V47, -1, He4, Ti44, -1, V47, Ti44, k_p_V47_to_He4_Ti44_derived, k_He4_Ti44_to_p_V47_reaclib
    };
}
#endif

#endif
