#ifndef ACTUAL_NETWORK_H
#define ACTUAL_NETWORK_H

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
        else if constexpr (spec == O16) {
            return 127.6193154119992_rt;
        }
        else if constexpr (spec == Ne20) {
            return 160.64482384000075_rt;
        }
        else if constexpr (spec == Mg23) {
            return 181.7258218679999_rt;
        }
        else if constexpr (spec == Mg24) {
            return 198.2570479679962_rt;
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
        else if constexpr (spec == Si30) {
            return 255.61964649599756_rt;
        }
        else if constexpr (spec == P29) {
            return 239.28581936000046_rt;
        }
        else if constexpr (spec == P31) {
            return 262.9161999600037_rt;
        }
        else if constexpr (spec == S30) {
            return 243.68124042399725_rt;
        }
        else if constexpr (spec == S32) {
            return 271.78016372399725_rt;
        }
        else if constexpr (spec == S34) {
            return 291.83895282400044_rt;
        }
        else if constexpr (spec == Cl33) {
            return 274.05689768799857_rt;
        }
        else if constexpr (spec == Cl35) {
            return 298.20976388800045_rt;
        }
        else if constexpr (spec == Ar34) {
            return 278.7208587519963_rt;
        }
        else if constexpr (spec == Ar36) {
            return 306.7167469519991_rt;
        }
        else if constexpr (spec == Ar38) {
            return 327.3426711519933_rt;
        }
        else if constexpr (spec == K37) {
            return 308.5743760159967_rt;
        }
        else if constexpr (spec == K39) {
            return 333.7240072160057_rt;
        }
        else if constexpr (spec == Ca38) {
            return 313.1216470799991_rt;
        }
        else if constexpr (spec == Ca40) {
            return 342.05218528000114_rt;
        }
        else if constexpr (spec == Ca42) {
            return 361.8957094799989_rt;
        }
        else if constexpr (spec == Sc41) {
            return 343.1371143439974_rt;
        }
        else if constexpr (spec == Sc43) {
            return 366.82549054399715_rt;
        }
        else if constexpr (spec == Ti42) {
            return 346.88807540800184_rt;
        }
        else if constexpr (spec == Ti44) {
            return 375.47496160800074_rt;
        }
        else if constexpr (spec == Ti46) {
            return 398.19726780800556_rt;
        }
        else if constexpr (spec == V45) {
            return 377.1017326720029_rt;
        }
        else if constexpr (spec == V47) {
            return 403.3650388719979_rt;
        }
        else if constexpr (spec == V49) {
            return 425.4628050719984_rt;
        }
        else if constexpr (spec == Cr46) {
            return 381.97630373600254_rt;
        }
        else if constexpr (spec == Cr48) {
            return 411.4679399359957_rt;
        }
        else if constexpr (spec == Cr50) {
            return 435.0509361359945_rt;
        }
        else if constexpr (spec == Cr52) {
            return 456.35172233599593_rt;
        }
        else if constexpr (spec == Mn49) {
            return 413.55581099999836_rt;
        }
        else if constexpr (spec == Mn51) {
            return 440.321747199996_rt;
        }
        else if constexpr (spec == Mn53) {
            return 462.9114833999993_rt;
        }
        else if constexpr (spec == Fe50) {
            return 417.7008820640039_rt;
        }
        else if constexpr (spec == Fe52) {
            return 447.6996182639923_rt;
        }
        else if constexpr (spec == Fe54) {
            return 471.76475446399854_rt;
        }
        else if constexpr (spec == Fe56) {
            return 492.2599506639962_rt;
        }
        else if constexpr (spec == Co53) {
            return 449.3158893279979_rt;
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
        else if constexpr (spec == O16) {
            return 2.6560180592333686e-23_rt;
        }
        else if constexpr (spec == Ne20) {
            return 3.3198227947612416e-23_rt;
        }
        else if constexpr (spec == Mg23) {
            return 3.8182640828719474e-23_rt;
        }
        else if constexpr (spec == Mg24) {
            return 3.9828098739467446e-23_rt;
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
        else if constexpr (spec == Si30) {
            return 4.977261628701015e-23_rt;
        }
        else if constexpr (spec == P29) {
            return 4.812541182037512e-23_rt;
        }
        else if constexpr (spec == P31) {
            return 5.14331418367544e-23_rt;
        }
        else if constexpr (spec == S30) {
            return 4.979110910837843e-23_rt;
        }
        else if constexpr (spec == S32) {
            return 5.309087322384128e-23_rt;
        }
        else if constexpr (spec == S34) {
            return 5.64049701803872e-23_rt;
        }
        else if constexpr (spec == Cl33) {
            return 5.476034741466158e-23_rt;
        }
        else if constexpr (spec == Cl35) {
            return 5.806714601585715e-23_rt;
        }
        else if constexpr (spec == Ar34) {
            return 5.642556598663247e-23_rt;
        }
        else if constexpr (spec == Ar36) {
            return 5.972551377884467e-23_rt;
        }
        else if constexpr (spec == Ar38) {
            return 6.303859972524341e-23_rt;
        }
        else if constexpr (spec == K37) {
            return 6.139573509201139e-23_rt;
        }
        else if constexpr (spec == K39) {
            return 6.470075679819664e-23_rt;
        }
        else if constexpr (spec == Ca38) {
            return 6.30611616828019e-23_rt;
        }
        else if constexpr (spec == Ca40) {
            return 6.635944331004904e-23_rt;
        }
        else if constexpr (spec == Ca42) {
            return 6.967392401113526e-23_rt;
        }
        else if constexpr (spec == Sc41) {
            return 6.80310420860826e-23_rt;
        }
        else if constexpr (spec == Sc43) {
            return 7.133866871591414e-23_rt;
        }
        else if constexpr (spec == Ti42) {
            return 6.969788822838791e-23_rt;
        }
        else if constexpr (spec == Ti44) {
            return 7.299678247096977e-23_rt;
        }
        else if constexpr (spec == Ti46) {
            return 7.630613127700392e-23_rt;
        }
        else if constexpr (spec == V45) {
            return 7.466741532590238e-23_rt;
        }
        else if constexpr (spec == V47) {
            return 7.797045172607207e-23_rt;
        }
        else if constexpr (spec == V49) {
            return 8.128091387578276e-23_rt;
        }
        else if constexpr (spec == Cr46) {
            return 7.633225845144595e-23_rt;
        }
        else if constexpr (spec == Cr48) {
            return 7.962953983065421e-23_rt;
        }
        else if constexpr (spec == Cr50) {
            return 8.293735431739904e-23_rt;
        }
        else if constexpr (spec == Cr52) {
            return 8.624923721300802e-23_rt;
        }
        else if constexpr (spec == Mn49) {
            return 8.129935070017476e-23_rt;
        }
        else if constexpr (spec == Mn51) {
            return 8.460149108098279e-23_rt;
        }
        else if constexpr (spec == Mn53) {
            return 8.79110762145079e-23_rt;
        }
        else if constexpr (spec == Fe50) {
            return 8.296549427759014e-23_rt;
        }
        else if constexpr (spec == Fe52) {
            return 8.626187166893794e-23_rt;
        }
        else if constexpr (spec == Fe54) {
            return 8.95688266630639e-23_rt;
        }
        else if constexpr (spec == Fe56) {
            return 9.288214565329033e-23_rt;
        }
        else if constexpr (spec == Co53) {
            return 8.793252324182073e-23_rt;
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
        k_C12_He4_to_O16_reaclib = 1,
        k_O16_He4_to_Ne20_reaclib = 2,
        k_Mg23_n_to_Mg24_reaclib = 3,
        k_Mg23_He4_to_Si27_reaclib = 4,
        k_Mg24_He4_to_Si28_reaclib = 5,
        k_Al27_p_to_Si28_reaclib = 6,
        k_Al27_He4_to_P31_reaclib = 7,
        k_Si27_n_to_Si28_reaclib = 8,
        k_Si28_p_to_P29_reaclib = 9,
        k_Si28_He4_to_S32_reaclib = 10,
        k_Si30_p_to_P31_reaclib = 11,
        k_Si30_He4_to_S34_reaclib = 12,
        k_P29_p_to_S30_reaclib = 13,
        k_P29_He4_to_Cl33_reaclib = 14,
        k_P31_p_to_S32_reaclib = 15,
        k_P31_He4_to_Cl35_reaclib = 16,
        k_S30_He4_to_Ar34_reaclib = 17,
        k_S32_p_to_Cl33_reaclib = 18,
        k_S32_He4_to_Ar36_reaclib = 19,
        k_S34_p_to_Cl35_reaclib = 20,
        k_S34_He4_to_Ar38_reaclib = 21,
        k_Cl33_p_to_Ar34_reaclib = 22,
        k_Cl33_He4_to_K37_reaclib = 23,
        k_Cl35_p_to_Ar36_reaclib = 24,
        k_Cl35_He4_to_K39_reaclib = 25,
        k_Ar34_He4_to_Ca38_reaclib = 26,
        k_Ar36_p_to_K37_reaclib = 27,
        k_Ar36_He4_to_Ca40_reaclib = 28,
        k_Ar38_p_to_K39_reaclib = 29,
        k_Ar38_He4_to_Ca42_reaclib = 30,
        k_K37_p_to_Ca38_reaclib = 31,
        k_K37_He4_to_Sc41_reaclib = 32,
        k_K39_p_to_Ca40_reaclib = 33,
        k_K39_He4_to_Sc43_reaclib = 34,
        k_Ca38_He4_to_Ti42_reaclib = 35,
        k_Ca40_p_to_Sc41_reaclib = 36,
        k_Ca40_He4_to_Ti44_reaclib = 37,
        k_Ca42_p_to_Sc43_reaclib = 38,
        k_Ca42_He4_to_Ti46_reaclib = 39,
        k_Sc41_p_to_Ti42_reaclib = 40,
        k_Sc41_He4_to_V45_reaclib = 41,
        k_Sc43_p_to_Ti44_reaclib = 42,
        k_Sc43_He4_to_V47_reaclib = 43,
        k_Ti42_He4_to_Cr46_reaclib = 44,
        k_Ti44_p_to_V45_reaclib = 45,
        k_Ti44_He4_to_Cr48_reaclib = 46,
        k_Ti46_p_to_V47_reaclib = 47,
        k_Ti46_He4_to_Cr50_reaclib = 48,
        k_V45_p_to_Cr46_reaclib = 49,
        k_V45_He4_to_Mn49_reaclib = 50,
        k_V47_p_to_Cr48_reaclib = 51,
        k_V47_He4_to_Mn51_reaclib = 52,
        k_V49_p_to_Cr50_reaclib = 53,
        k_V49_He4_to_Mn53_reaclib = 54,
        k_Cr46_He4_to_Fe50_reaclib = 55,
        k_Cr48_p_to_Mn49_reaclib = 56,
        k_Cr48_He4_to_Fe52_reaclib = 57,
        k_Cr50_p_to_Mn51_reaclib = 58,
        k_Cr50_He4_to_Fe54_reaclib = 59,
        k_Cr52_p_to_Mn53_reaclib = 60,
        k_Cr52_He4_to_Fe56_reaclib = 61,
        k_Mn49_p_to_Fe50_reaclib = 62,
        k_Mn49_He4_to_Co53_reaclib = 63,
        k_Mn51_p_to_Fe52_reaclib = 64,
        k_Mn51_He4_to_Co55_reaclib = 65,
        k_Mn53_p_to_Fe54_reaclib = 66,
        k_Mn53_He4_to_Co57_reaclib = 67,
        k_Fe52_p_to_Co53_reaclib = 68,
        k_Fe52_He4_to_Ni56_reaclib = 69,
        k_Fe54_p_to_Co55_reaclib = 70,
        k_Fe54_He4_to_Ni58_reaclib = 71,
        k_Fe56_p_to_Co57_reaclib = 72,
        k_Co55_n_to_Co56_reaclib = 73,
        k_Co55_p_to_Ni56_reaclib = 74,
        k_Co55_He4_to_Cu59_reaclib = 75,
        k_Co56_n_to_Co57_reaclib = 76,
        k_Co57_p_to_Ni58_reaclib = 77,
        k_Ni56_He4_to_Zn60_reaclib = 78,
        k_Ni58_p_to_Cu59_reaclib = 79,
        k_Cu59_p_to_Zn60_reaclib = 80,
        k_O16_C12_to_p_Al27_reaclib = 81,
        k_O16_C12_to_He4_Mg24_reaclib = 82,
        k_O16_O16_to_p_P31_reaclib = 83,
        k_O16_O16_to_He4_Si28_reaclib = 84,
        k_Ne20_C12_to_p_P31_reaclib = 85,
        k_Ne20_C12_to_He4_Si28_reaclib = 86,
        k_Mg23_n_to_He4_Ne20_reaclib = 87,
        k_Mg23_n_to_C12_C12_reaclib = 88,
        k_Al27_p_to_He4_Mg24_reaclib = 89,
        k_Al27_He4_to_p_Si30_reaclib = 90,
        k_Si27_n_to_p_Al27_reaclib = 91,
        k_Si27_n_to_He4_Mg24_reaclib = 92,
        k_Si27_n_to_C12_O16_reaclib = 93,
        k_P29_He4_to_p_S32_reaclib = 94,
        k_P31_p_to_He4_Si28_reaclib = 95,
        k_P31_He4_to_p_S34_reaclib = 96,
        k_S30_n_to_He4_Si27_reaclib = 97,
        k_S30_He4_to_p_Cl33_reaclib = 98,
        k_Cl33_He4_to_p_Ar36_reaclib = 99,
        k_Cl35_p_to_He4_S32_reaclib = 100,
        k_Cl35_He4_to_p_Ar38_reaclib = 101,
        k_Ar34_He4_to_p_K37_reaclib = 102,
        k_K37_He4_to_p_Ca40_reaclib = 103,
        k_K39_p_to_He4_Ar36_reaclib = 104,
        k_Ca38_He4_to_p_Sc41_reaclib = 105,
        k_Ca42_p_to_He4_K39_reaclib = 106,
        k_Sc41_He4_to_p_Ti44_reaclib = 107,
        k_Sc43_p_to_He4_Ca40_reaclib = 108,
        k_Sc43_He4_to_p_Ti46_reaclib = 109,
        k_Ti42_He4_to_p_V45_reaclib = 110,
        k_Ti44_He4_to_p_V47_reaclib = 111,
        k_V45_He4_to_p_Cr48_reaclib = 112,
        k_V47_He4_to_p_Cr50_reaclib = 113,
        k_V49_p_to_He4_Ti46_reaclib = 114,
        k_V49_He4_to_p_Cr52_reaclib = 115,
        k_Cr46_He4_to_p_Mn49_reaclib = 116,
        k_Cr48_He4_to_p_Mn51_reaclib = 117,
        k_Mn49_He4_to_p_Fe52_reaclib = 118,
        k_Mn51_He4_to_p_Fe54_reaclib = 119,
        k_Mn53_p_to_He4_Cr50_reaclib = 120,
        k_Mn53_He4_to_p_Fe56_reaclib = 121,
        k_Fe50_He4_to_p_Co53_reaclib = 122,
        k_Fe52_He4_to_p_Co55_reaclib = 123,
        k_Co53_He4_to_p_Ni56_reaclib = 124,
        k_Co55_He4_to_p_Ni58_reaclib = 125,
        k_Co56_n_to_p_Fe56_reaclib = 126,
        k_Co56_n_to_He4_Mn53_reaclib = 127,
        k_Co57_p_to_He4_Fe54_reaclib = 128,
        k_Ni56_n_to_p_Co56_reaclib = 129,
        k_Cu59_n_to_He4_Co56_reaclib = 130,
        k_Cu59_p_to_He4_Ni56_reaclib = 131,
        k_He4_He4_He4_to_C12_reaclib = 132,
        k_Ni56_n_to_Ni57_reaclib = 133,
        k_Ni57_n_to_Ni58_reaclib = 134,
        k_Co53_n_to_Co54_reaclib = 135,
        k_Co54_n_to_Co55_reaclib = 136,
        k_Fe50_n_to_Fe51_reaclib = 137,
        k_Fe51_n_to_Fe52_reaclib = 138,
        k_Fe52_n_to_Fe53_reaclib = 139,
        k_Fe53_n_to_Fe54_reaclib = 140,
        k_Fe54_n_to_Fe55_reaclib = 141,
        k_Fe55_n_to_Fe56_reaclib = 142,
        k_Mn49_n_to_Mn50_reaclib = 143,
        k_Mn50_n_to_Mn51_reaclib = 144,
        k_Mn51_n_to_Mn52_reaclib = 145,
        k_Mn52_n_to_Mn53_reaclib = 146,
        k_Cr46_n_to_Cr47_reaclib = 147,
        k_Cr47_n_to_Cr48_reaclib = 148,
        k_Cr48_n_to_Cr49_reaclib = 149,
        k_Cr49_n_to_Cr50_reaclib = 150,
        k_Cr50_n_to_Cr51_reaclib = 151,
        k_Cr51_n_to_Cr52_reaclib = 152,
        k_V45_n_to_V46_reaclib = 153,
        k_V46_n_to_V47_reaclib = 154,
        k_V47_n_to_V48_reaclib = 155,
        k_V48_n_to_V49_reaclib = 156,
        k_Ti42_n_to_Ti43_reaclib = 157,
        k_Ti43_n_to_Ti44_reaclib = 158,
        k_Ti44_n_to_Ti45_reaclib = 159,
        k_Ti45_n_to_Ti46_reaclib = 160,
        k_Sc41_n_to_Sc42_reaclib = 161,
        k_Sc42_n_to_Sc43_reaclib = 162,
        k_Ca38_n_to_Ca39_reaclib = 163,
        k_Ca39_n_to_Ca40_reaclib = 164,
        k_Ca40_n_to_Ca41_reaclib = 165,
        k_Ca41_n_to_Ca42_reaclib = 166,
        k_Ar34_n_to_Ar35_reaclib = 167,
        k_Ar35_n_to_Ar36_reaclib = 168,
        k_Ar36_n_to_Ar37_reaclib = 169,
        k_Ar37_n_to_Ar38_reaclib = 170,
        k_Si28_n_to_Si29_reaclib = 171,
        k_Si29_n_to_Si30_reaclib = 172,
        k_S30_n_to_S31_reaclib = 173,
        k_S31_n_to_S32_reaclib = 174,
        k_S32_n_to_S33_reaclib = 175,
        k_S33_n_to_S34_reaclib = 176,
        k_C12_C12_to_p_Na23_reaclib = 177,
        k_Na23_p_to_He4_Ne20_reaclib = 178,
        k_Na23_p_to_Mg24_reaclib = 179,
        k_C12_C12_to_He4_Ne20_reaclib = 180,
        k_Ne20_He4_to_Mg24_reaclib = 181,
        k_Al27_to_Si27_beta_neg_weaktab = 182,
        k_Co56_to_Fe56_electron_capture_weaktab = 183,
        k_Co56_to_Ni56_beta_neg_weaktab = 184,
        k_Fe56_to_Co56_beta_neg_weaktab = 185,
        k_Ni56_to_Co56_electron_capture_weaktab = 186,
        k_n_to_p_beta_neg_weaktab = 187,
        k_p_to_n_electron_capture_weaktab = 188,
        k_Si27_to_Al27_electron_capture_weaktab = 189,
        k_Ni56_n_n_to_Ni58_approx = 190,
        k_Ni58_to_Ni56_n_n_approx = 191,
        k_Co53_n_n_to_Co55_approx = 192,
        k_Co55_to_Co53_n_n_approx = 193,
        k_Fe50_n_n_to_Fe52_approx = 194,
        k_Fe52_to_Fe50_n_n_approx = 195,
        k_Fe52_n_n_to_Fe54_approx = 196,
        k_Fe54_to_Fe52_n_n_approx = 197,
        k_Fe54_n_n_to_Fe56_approx = 198,
        k_Fe56_to_Fe54_n_n_approx = 199,
        k_Mn49_n_n_to_Mn51_approx = 200,
        k_Mn51_to_Mn49_n_n_approx = 201,
        k_Mn51_n_n_to_Mn53_approx = 202,
        k_Mn53_to_Mn51_n_n_approx = 203,
        k_Cr46_n_n_to_Cr48_approx = 204,
        k_Cr48_to_Cr46_n_n_approx = 205,
        k_Cr48_n_n_to_Cr50_approx = 206,
        k_Cr50_to_Cr48_n_n_approx = 207,
        k_Cr50_n_n_to_Cr52_approx = 208,
        k_Cr52_to_Cr50_n_n_approx = 209,
        k_V45_n_n_to_V47_approx = 210,
        k_V47_to_V45_n_n_approx = 211,
        k_V47_n_n_to_V49_approx = 212,
        k_V49_to_V47_n_n_approx = 213,
        k_Ti42_n_n_to_Ti44_approx = 214,
        k_Ti44_to_Ti42_n_n_approx = 215,
        k_Ti44_n_n_to_Ti46_approx = 216,
        k_Ti46_to_Ti44_n_n_approx = 217,
        k_Sc41_n_n_to_Sc43_approx = 218,
        k_Sc43_to_Sc41_n_n_approx = 219,
        k_Ca38_n_n_to_Ca40_approx = 220,
        k_Ca40_to_Ca38_n_n_approx = 221,
        k_Ca40_n_n_to_Ca42_approx = 222,
        k_Ca42_to_Ca40_n_n_approx = 223,
        k_Ar34_n_n_to_Ar36_approx = 224,
        k_Ar36_to_Ar34_n_n_approx = 225,
        k_Ar36_n_n_to_Ar38_approx = 226,
        k_Ar38_to_Ar36_n_n_approx = 227,
        k_Si28_n_n_to_Si30_approx = 228,
        k_Si30_to_Si28_n_n_approx = 229,
        k_S30_n_n_to_S32_approx = 230,
        k_S32_to_S30_n_n_approx = 231,
        k_S32_n_n_to_S34_approx = 232,
        k_S34_to_S32_n_n_approx = 233,
        k_C12_C12_to_Mg24_approx = 234,
        k_Mg24_to_C12_C12_approx = 235,
        k_C12_C12_to_Ne20_He4_approx = 236,
        k_Ne20_He4_to_C12_C12_approx = 237,
        k_Ne20_He4_to_Mg24_approx = 238,
        k_Mg24_to_Ne20_He4_approx = 239,
        k_O16_to_He4_C12_derived = 240,
        k_Ne20_to_He4_O16_derived = 241,
        k_Mg24_to_n_Mg23_derived = 242,
        k_Si27_to_He4_Mg23_derived = 243,
        k_Si28_to_n_Si27_derived = 244,
        k_Si28_to_p_Al27_derived = 245,
        k_Si28_to_He4_Mg24_derived = 246,
        k_P29_to_p_Si28_derived = 247,
        k_P31_to_p_Si30_derived = 248,
        k_P31_to_He4_Al27_derived = 249,
        k_S30_to_p_P29_derived = 250,
        k_S32_to_p_P31_derived = 251,
        k_S32_to_He4_Si28_derived = 252,
        k_S34_to_He4_Si30_derived = 253,
        k_Cl33_to_p_S32_derived = 254,
        k_Cl33_to_He4_P29_derived = 255,
        k_Cl35_to_p_S34_derived = 256,
        k_Cl35_to_He4_P31_derived = 257,
        k_Ar34_to_p_Cl33_derived = 258,
        k_Ar34_to_He4_S30_derived = 259,
        k_Ar36_to_p_Cl35_derived = 260,
        k_Ar36_to_He4_S32_derived = 261,
        k_Ar38_to_He4_S34_derived = 262,
        k_K37_to_p_Ar36_derived = 263,
        k_K37_to_He4_Cl33_derived = 264,
        k_K39_to_p_Ar38_derived = 265,
        k_K39_to_He4_Cl35_derived = 266,
        k_Ca38_to_p_K37_derived = 267,
        k_Ca38_to_He4_Ar34_derived = 268,
        k_Ca40_to_p_K39_derived = 269,
        k_Ca40_to_He4_Ar36_derived = 270,
        k_Ca42_to_He4_Ar38_derived = 271,
        k_Sc41_to_p_Ca40_derived = 272,
        k_Sc41_to_He4_K37_derived = 273,
        k_Sc43_to_p_Ca42_derived = 274,
        k_Sc43_to_He4_K39_derived = 275,
        k_Ti42_to_p_Sc41_derived = 276,
        k_Ti42_to_He4_Ca38_derived = 277,
        k_Ti44_to_p_Sc43_derived = 278,
        k_Ti44_to_He4_Ca40_derived = 279,
        k_Ti46_to_He4_Ca42_derived = 280,
        k_V45_to_p_Ti44_derived = 281,
        k_V45_to_He4_Sc41_derived = 282,
        k_V47_to_p_Ti46_derived = 283,
        k_V47_to_He4_Sc43_derived = 284,
        k_Cr46_to_p_V45_derived = 285,
        k_Cr46_to_He4_Ti42_derived = 286,
        k_Cr48_to_p_V47_derived = 287,
        k_Cr48_to_He4_Ti44_derived = 288,
        k_Cr50_to_p_V49_derived = 289,
        k_Cr50_to_He4_Ti46_derived = 290,
        k_Mn49_to_p_Cr48_derived = 291,
        k_Mn49_to_He4_V45_derived = 292,
        k_Mn51_to_p_Cr50_derived = 293,
        k_Mn51_to_He4_V47_derived = 294,
        k_Mn53_to_p_Cr52_derived = 295,
        k_Mn53_to_He4_V49_derived = 296,
        k_Fe50_to_p_Mn49_derived = 297,
        k_Fe50_to_He4_Cr46_derived = 298,
        k_Fe52_to_p_Mn51_derived = 299,
        k_Fe52_to_He4_Cr48_derived = 300,
        k_Fe54_to_p_Mn53_derived = 301,
        k_Fe54_to_He4_Cr50_derived = 302,
        k_Fe56_to_He4_Cr52_derived = 303,
        k_Co53_to_p_Fe52_derived = 304,
        k_Co53_to_He4_Mn49_derived = 305,
        k_Co55_to_p_Fe54_derived = 306,
        k_Co55_to_He4_Mn51_derived = 307,
        k_Co56_to_n_Co55_derived = 308,
        k_Co57_to_n_Co56_derived = 309,
        k_Co57_to_p_Fe56_derived = 310,
        k_Co57_to_He4_Mn53_derived = 311,
        k_Ni56_to_p_Co55_derived = 312,
        k_Ni56_to_He4_Fe52_derived = 313,
        k_Ni58_to_p_Co57_derived = 314,
        k_Ni58_to_He4_Fe54_derived = 315,
        k_Cu59_to_p_Ni58_derived = 316,
        k_Cu59_to_He4_Co55_derived = 317,
        k_Zn60_to_p_Cu59_derived = 318,
        k_Zn60_to_He4_Ni56_derived = 319,
        k_C12_to_3He4_derived = 320,
        k_C12_C12_to_n_Mg23_derived = 321,
        k_O16_C12_to_n_Si27_derived = 322,
        k_Ne20_He4_to_n_Mg23_derived = 323,
        k_Mg24_He4_to_n_Si27_derived = 324,
        k_Mg24_He4_to_p_Al27_derived = 325,
        k_Mg24_He4_to_C12_O16_derived = 326,
        k_Al27_p_to_n_Si27_derived = 327,
        k_Al27_p_to_C12_O16_derived = 328,
        k_Si27_He4_to_n_S30_derived = 329,
        k_Si28_He4_to_p_P31_derived = 330,
        k_Si28_He4_to_C12_Ne20_derived = 331,
        k_Si28_He4_to_O16_O16_derived = 332,
        k_Si30_p_to_He4_Al27_derived = 333,
        k_P31_p_to_C12_Ne20_derived = 334,
        k_P31_p_to_O16_O16_derived = 335,
        k_S32_p_to_He4_P29_derived = 336,
        k_S32_He4_to_p_Cl35_derived = 337,
        k_S34_p_to_He4_P31_derived = 338,
        k_Cl33_p_to_He4_S30_derived = 339,
        k_Ar36_p_to_He4_Cl33_derived = 340,
        k_Ar36_He4_to_p_K39_derived = 341,
        k_Ar38_p_to_He4_Cl35_derived = 342,
        k_K37_p_to_He4_Ar34_derived = 343,
        k_K39_He4_to_p_Ca42_derived = 344,
        k_Ca40_p_to_He4_K37_derived = 345,
        k_Ca40_He4_to_p_Sc43_derived = 346,
        k_Sc41_p_to_He4_Ca38_derived = 347,
        k_Ti44_p_to_He4_Sc41_derived = 348,
        k_Ti46_p_to_He4_Sc43_derived = 349,
        k_Ti46_He4_to_p_V49_derived = 350,
        k_V45_p_to_He4_Ti42_derived = 351,
        k_V47_p_to_He4_Ti44_derived = 352,
        k_Cr48_p_to_He4_V45_derived = 353,
        k_Cr50_p_to_He4_V47_derived = 354,
        k_Cr50_He4_to_p_Mn53_derived = 355,
        k_Cr52_p_to_He4_V49_derived = 356,
        k_Mn49_p_to_He4_Cr46_derived = 357,
        k_Mn51_p_to_He4_Cr48_derived = 358,
        k_Mn53_He4_to_n_Co56_derived = 359,
        k_Fe52_p_to_He4_Mn49_derived = 360,
        k_Fe54_p_to_He4_Mn51_derived = 361,
        k_Fe54_He4_to_p_Co57_derived = 362,
        k_Fe56_p_to_n_Co56_derived = 363,
        k_Fe56_p_to_He4_Mn53_derived = 364,
        k_Co53_p_to_He4_Fe50_derived = 365,
        k_Co55_p_to_He4_Fe52_derived = 366,
        k_Co56_p_to_n_Ni56_derived = 367,
        k_Co56_He4_to_n_Cu59_derived = 368,
        k_Ni56_p_to_He4_Co53_derived = 369,
        k_Ni56_He4_to_p_Cu59_derived = 370,
        k_Ni58_p_to_He4_Co55_derived = 371,
        k_Ni57_to_n_Ni56_derived = 372,
        k_Ni58_to_n_Ni57_derived = 373,
        k_Co54_to_n_Co53_derived = 374,
        k_Co55_to_n_Co54_derived = 375,
        k_Fe51_to_n_Fe50_derived = 376,
        k_Fe52_to_n_Fe51_derived = 377,
        k_Fe53_to_n_Fe52_derived = 378,
        k_Fe54_to_n_Fe53_derived = 379,
        k_Fe55_to_n_Fe54_derived = 380,
        k_Fe56_to_n_Fe55_derived = 381,
        k_Mn50_to_n_Mn49_derived = 382,
        k_Mn51_to_n_Mn50_derived = 383,
        k_Mn52_to_n_Mn51_derived = 384,
        k_Mn53_to_n_Mn52_derived = 385,
        k_Cr47_to_n_Cr46_derived = 386,
        k_Cr48_to_n_Cr47_derived = 387,
        k_Cr49_to_n_Cr48_derived = 388,
        k_Cr50_to_n_Cr49_derived = 389,
        k_Cr51_to_n_Cr50_derived = 390,
        k_Cr52_to_n_Cr51_derived = 391,
        k_V46_to_n_V45_derived = 392,
        k_V47_to_n_V46_derived = 393,
        k_V48_to_n_V47_derived = 394,
        k_V49_to_n_V48_derived = 395,
        k_Ti43_to_n_Ti42_derived = 396,
        k_Ti44_to_n_Ti43_derived = 397,
        k_Ti45_to_n_Ti44_derived = 398,
        k_Ti46_to_n_Ti45_derived = 399,
        k_Sc42_to_n_Sc41_derived = 400,
        k_Sc43_to_n_Sc42_derived = 401,
        k_Ca39_to_n_Ca38_derived = 402,
        k_Ca40_to_n_Ca39_derived = 403,
        k_Ca41_to_n_Ca40_derived = 404,
        k_Ca42_to_n_Ca41_derived = 405,
        k_Ar35_to_n_Ar34_derived = 406,
        k_Ar36_to_n_Ar35_derived = 407,
        k_Ar37_to_n_Ar36_derived = 408,
        k_Ar38_to_n_Ar37_derived = 409,
        k_Si29_to_n_Si28_derived = 410,
        k_Si30_to_n_Si29_derived = 411,
        k_S31_to_n_S30_derived = 412,
        k_S32_to_n_S31_derived = 413,
        k_S33_to_n_S32_derived = 414,
        k_S34_to_n_S33_derived = 415,
        k_Na23_p_to_C12_C12_derived = 416,
        k_Mg24_to_p_Na23_derived = 417,
        k_Ne20_He4_to_C12_C12_derived = 418,
        k_Ne20_He4_to_p_Na23_derived = 419,
        k_Mg24_to_He4_Ne20_derived = 420,
        NumRates = k_Mg24_to_He4_Ne20_derived
    };

    enum ScreenPairs :
    std::uint8_t
    {
        k_He4_Mn53 = 1,
        k_p_Ti44 = 2,
        k_p_Ca42 = 3,
        k_p_Sc43 = 4,
        k_p_Fe56 = 5,
        k_p_Ar38 = 6,
        k_He4_Ti44 = 7,
        k_He4_Ca42 = 8,
        k_He4_Mg24 = 9,
        k_C12_O16 = 10,
        k_He4_Ar38 = 11,
        k_p_V47 = 12,
        k_He4_Co53 = 13,
        k_p_Co55 = 14,
        k_C12_C12 = 15,
        k_He4_Cl33 = 16,
        k_p_Al27 = 17,
        k_He4_Ti46 = 18,
        k_He4_V47 = 19,
        k_He4_Ni56 = 20,
        k_p_P31 = 21,
        k_p_Ni58 = 22,
        k_He4_V49 = 23,
        k_He4_Cr46 = 24,
        k_O16_O16 = 25,
        k_He4_P31 = 26,
        k_p_Mn49 = 27,
        k_p_Cr52 = 28,
        k_p_S34 = 29,
        k_He4_Fe50 = 30,
        k_He4_Ca40 = 31,
        k_He4_Mn49 = 32,
        k_He4_Ar36 = 33,
        k_He4_Cr52 = 34,
        k_He4_S34 = 35,
        k_He4_S30 = 36,
        k_He4_Sc43 = 37,
        k_He4_Ne20 = 38,
        k_p_Co56 = 39,
        k_p_Si28 = 40,
        k_p_Fe54 = 41,
        k_He4_Co55 = 42,
        k_He4_Al27 = 43,
        k_p_V45 = 44,
        k_p_Cl35 = 45,
        k_p_Si30 = 46,
        k_He4_O16 = 47,
        k_He4_Si30 = 48,
        k_p_P29 = 49,
        k_He4_Mg23 = 50,
        k_He4_P29 = 51,
        k_p_Fe52 = 52,
        k_p_Cr50 = 53,
        k_p_S32 = 54,
        k_He4_Ca38 = 55,
        k_He4_Ar34 = 56,
        k_p_Cu59 = 57,
        k_He4_S32 = 58,
        k_p_Sc41 = 59,
        k_He4_Ti42 = 60,
        k_p_Mn53 = 61,
        k_He4_Sc41 = 62,
        k_He4_Co56 = 63,
        k_He4_Si28 = 64,
        k_He4_Fe54 = 65,
        k_He4_Be8 = 66,
        k_He4_V45 = 67,
        k_p_Co57 = 68,
        k_p_K39 = 69,
        k_He4_Cl35 = 70,
        k_p_Co53 = 71,
        k_He4_C12 = 72,
        k_He4_K39 = 73,
        k_p_Cr48 = 74,
        k_p_Cl33 = 75,
        k_p_Ti46 = 76,
        k_He4_Cr48 = 77,
        k_He4_Si27 = 78,
        k_p_K37 = 79,
        k_p_Ni56 = 80,
        k_p_Mn51 = 81,
        k_p_Na23 = 82,
        k_p_V49 = 83,
        k_He4_Fe52 = 84,
        k_He4_K37 = 85,
        k_He4_Cr50 = 86,
        k_He4_Mn51 = 87,
        k_He4_He4 = 88,
        k_C12_Ne20 = 89,
        k_p_Ca40 = 90,
        k_p_Ar36 = 91,
        NumScreenPairs = k_p_Ar36
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "C12_He4_to_O16_reaclib",  // 1,
        "O16_He4_to_Ne20_reaclib",  // 2,
        "Mg23_n_to_Mg24_reaclib",  // 3,
        "Mg23_He4_to_Si27_reaclib",  // 4,
        "Mg24_He4_to_Si28_reaclib",  // 5,
        "Al27_p_to_Si28_reaclib",  // 6,
        "Al27_He4_to_P31_reaclib",  // 7,
        "Si27_n_to_Si28_reaclib",  // 8,
        "Si28_p_to_P29_reaclib",  // 9,
        "Si28_He4_to_S32_reaclib",  // 10,
        "Si30_p_to_P31_reaclib",  // 11,
        "Si30_He4_to_S34_reaclib",  // 12,
        "P29_p_to_S30_reaclib",  // 13,
        "P29_He4_to_Cl33_reaclib",  // 14,
        "P31_p_to_S32_reaclib",  // 15,
        "P31_He4_to_Cl35_reaclib",  // 16,
        "S30_He4_to_Ar34_reaclib",  // 17,
        "S32_p_to_Cl33_reaclib",  // 18,
        "S32_He4_to_Ar36_reaclib",  // 19,
        "S34_p_to_Cl35_reaclib",  // 20,
        "S34_He4_to_Ar38_reaclib",  // 21,
        "Cl33_p_to_Ar34_reaclib",  // 22,
        "Cl33_He4_to_K37_reaclib",  // 23,
        "Cl35_p_to_Ar36_reaclib",  // 24,
        "Cl35_He4_to_K39_reaclib",  // 25,
        "Ar34_He4_to_Ca38_reaclib",  // 26,
        "Ar36_p_to_K37_reaclib",  // 27,
        "Ar36_He4_to_Ca40_reaclib",  // 28,
        "Ar38_p_to_K39_reaclib",  // 29,
        "Ar38_He4_to_Ca42_reaclib",  // 30,
        "K37_p_to_Ca38_reaclib",  // 31,
        "K37_He4_to_Sc41_reaclib",  // 32,
        "K39_p_to_Ca40_reaclib",  // 33,
        "K39_He4_to_Sc43_reaclib",  // 34,
        "Ca38_He4_to_Ti42_reaclib",  // 35,
        "Ca40_p_to_Sc41_reaclib",  // 36,
        "Ca40_He4_to_Ti44_reaclib",  // 37,
        "Ca42_p_to_Sc43_reaclib",  // 38,
        "Ca42_He4_to_Ti46_reaclib",  // 39,
        "Sc41_p_to_Ti42_reaclib",  // 40,
        "Sc41_He4_to_V45_reaclib",  // 41,
        "Sc43_p_to_Ti44_reaclib",  // 42,
        "Sc43_He4_to_V47_reaclib",  // 43,
        "Ti42_He4_to_Cr46_reaclib",  // 44,
        "Ti44_p_to_V45_reaclib",  // 45,
        "Ti44_He4_to_Cr48_reaclib",  // 46,
        "Ti46_p_to_V47_reaclib",  // 47,
        "Ti46_He4_to_Cr50_reaclib",  // 48,
        "V45_p_to_Cr46_reaclib",  // 49,
        "V45_He4_to_Mn49_reaclib",  // 50,
        "V47_p_to_Cr48_reaclib",  // 51,
        "V47_He4_to_Mn51_reaclib",  // 52,
        "V49_p_to_Cr50_reaclib",  // 53,
        "V49_He4_to_Mn53_reaclib",  // 54,
        "Cr46_He4_to_Fe50_reaclib",  // 55,
        "Cr48_p_to_Mn49_reaclib",  // 56,
        "Cr48_He4_to_Fe52_reaclib",  // 57,
        "Cr50_p_to_Mn51_reaclib",  // 58,
        "Cr50_He4_to_Fe54_reaclib",  // 59,
        "Cr52_p_to_Mn53_reaclib",  // 60,
        "Cr52_He4_to_Fe56_reaclib",  // 61,
        "Mn49_p_to_Fe50_reaclib",  // 62,
        "Mn49_He4_to_Co53_reaclib",  // 63,
        "Mn51_p_to_Fe52_reaclib",  // 64,
        "Mn51_He4_to_Co55_reaclib",  // 65,
        "Mn53_p_to_Fe54_reaclib",  // 66,
        "Mn53_He4_to_Co57_reaclib",  // 67,
        "Fe52_p_to_Co53_reaclib",  // 68,
        "Fe52_He4_to_Ni56_reaclib",  // 69,
        "Fe54_p_to_Co55_reaclib",  // 70,
        "Fe54_He4_to_Ni58_reaclib",  // 71,
        "Fe56_p_to_Co57_reaclib",  // 72,
        "Co55_n_to_Co56_reaclib",  // 73,
        "Co55_p_to_Ni56_reaclib",  // 74,
        "Co55_He4_to_Cu59_reaclib",  // 75,
        "Co56_n_to_Co57_reaclib",  // 76,
        "Co57_p_to_Ni58_reaclib",  // 77,
        "Ni56_He4_to_Zn60_reaclib",  // 78,
        "Ni58_p_to_Cu59_reaclib",  // 79,
        "Cu59_p_to_Zn60_reaclib",  // 80,
        "O16_C12_to_p_Al27_reaclib",  // 81,
        "O16_C12_to_He4_Mg24_reaclib",  // 82,
        "O16_O16_to_p_P31_reaclib",  // 83,
        "O16_O16_to_He4_Si28_reaclib",  // 84,
        "Ne20_C12_to_p_P31_reaclib",  // 85,
        "Ne20_C12_to_He4_Si28_reaclib",  // 86,
        "Mg23_n_to_He4_Ne20_reaclib",  // 87,
        "Mg23_n_to_C12_C12_reaclib",  // 88,
        "Al27_p_to_He4_Mg24_reaclib",  // 89,
        "Al27_He4_to_p_Si30_reaclib",  // 90,
        "Si27_n_to_p_Al27_reaclib",  // 91,
        "Si27_n_to_He4_Mg24_reaclib",  // 92,
        "Si27_n_to_C12_O16_reaclib",  // 93,
        "P29_He4_to_p_S32_reaclib",  // 94,
        "P31_p_to_He4_Si28_reaclib",  // 95,
        "P31_He4_to_p_S34_reaclib",  // 96,
        "S30_n_to_He4_Si27_reaclib",  // 97,
        "S30_He4_to_p_Cl33_reaclib",  // 98,
        "Cl33_He4_to_p_Ar36_reaclib",  // 99,
        "Cl35_p_to_He4_S32_reaclib",  // 100,
        "Cl35_He4_to_p_Ar38_reaclib",  // 101,
        "Ar34_He4_to_p_K37_reaclib",  // 102,
        "K37_He4_to_p_Ca40_reaclib",  // 103,
        "K39_p_to_He4_Ar36_reaclib",  // 104,
        "Ca38_He4_to_p_Sc41_reaclib",  // 105,
        "Ca42_p_to_He4_K39_reaclib",  // 106,
        "Sc41_He4_to_p_Ti44_reaclib",  // 107,
        "Sc43_p_to_He4_Ca40_reaclib",  // 108,
        "Sc43_He4_to_p_Ti46_reaclib",  // 109,
        "Ti42_He4_to_p_V45_reaclib",  // 110,
        "Ti44_He4_to_p_V47_reaclib",  // 111,
        "V45_He4_to_p_Cr48_reaclib",  // 112,
        "V47_He4_to_p_Cr50_reaclib",  // 113,
        "V49_p_to_He4_Ti46_reaclib",  // 114,
        "V49_He4_to_p_Cr52_reaclib",  // 115,
        "Cr46_He4_to_p_Mn49_reaclib",  // 116,
        "Cr48_He4_to_p_Mn51_reaclib",  // 117,
        "Mn49_He4_to_p_Fe52_reaclib",  // 118,
        "Mn51_He4_to_p_Fe54_reaclib",  // 119,
        "Mn53_p_to_He4_Cr50_reaclib",  // 120,
        "Mn53_He4_to_p_Fe56_reaclib",  // 121,
        "Fe50_He4_to_p_Co53_reaclib",  // 122,
        "Fe52_He4_to_p_Co55_reaclib",  // 123,
        "Co53_He4_to_p_Ni56_reaclib",  // 124,
        "Co55_He4_to_p_Ni58_reaclib",  // 125,
        "Co56_n_to_p_Fe56_reaclib",  // 126,
        "Co56_n_to_He4_Mn53_reaclib",  // 127,
        "Co57_p_to_He4_Fe54_reaclib",  // 128,
        "Ni56_n_to_p_Co56_reaclib",  // 129,
        "Cu59_n_to_He4_Co56_reaclib",  // 130,
        "Cu59_p_to_He4_Ni56_reaclib",  // 131,
        "He4_He4_He4_to_C12_reaclib",  // 132,
        "Ni56_n_to_Ni57_reaclib",  // 133,
        "Ni57_n_to_Ni58_reaclib",  // 134,
        "Co53_n_to_Co54_reaclib",  // 135,
        "Co54_n_to_Co55_reaclib",  // 136,
        "Fe50_n_to_Fe51_reaclib",  // 137,
        "Fe51_n_to_Fe52_reaclib",  // 138,
        "Fe52_n_to_Fe53_reaclib",  // 139,
        "Fe53_n_to_Fe54_reaclib",  // 140,
        "Fe54_n_to_Fe55_reaclib",  // 141,
        "Fe55_n_to_Fe56_reaclib",  // 142,
        "Mn49_n_to_Mn50_reaclib",  // 143,
        "Mn50_n_to_Mn51_reaclib",  // 144,
        "Mn51_n_to_Mn52_reaclib",  // 145,
        "Mn52_n_to_Mn53_reaclib",  // 146,
        "Cr46_n_to_Cr47_reaclib",  // 147,
        "Cr47_n_to_Cr48_reaclib",  // 148,
        "Cr48_n_to_Cr49_reaclib",  // 149,
        "Cr49_n_to_Cr50_reaclib",  // 150,
        "Cr50_n_to_Cr51_reaclib",  // 151,
        "Cr51_n_to_Cr52_reaclib",  // 152,
        "V45_n_to_V46_reaclib",  // 153,
        "V46_n_to_V47_reaclib",  // 154,
        "V47_n_to_V48_reaclib",  // 155,
        "V48_n_to_V49_reaclib",  // 156,
        "Ti42_n_to_Ti43_reaclib",  // 157,
        "Ti43_n_to_Ti44_reaclib",  // 158,
        "Ti44_n_to_Ti45_reaclib",  // 159,
        "Ti45_n_to_Ti46_reaclib",  // 160,
        "Sc41_n_to_Sc42_reaclib",  // 161,
        "Sc42_n_to_Sc43_reaclib",  // 162,
        "Ca38_n_to_Ca39_reaclib",  // 163,
        "Ca39_n_to_Ca40_reaclib",  // 164,
        "Ca40_n_to_Ca41_reaclib",  // 165,
        "Ca41_n_to_Ca42_reaclib",  // 166,
        "Ar34_n_to_Ar35_reaclib",  // 167,
        "Ar35_n_to_Ar36_reaclib",  // 168,
        "Ar36_n_to_Ar37_reaclib",  // 169,
        "Ar37_n_to_Ar38_reaclib",  // 170,
        "Si28_n_to_Si29_reaclib",  // 171,
        "Si29_n_to_Si30_reaclib",  // 172,
        "S30_n_to_S31_reaclib",  // 173,
        "S31_n_to_S32_reaclib",  // 174,
        "S32_n_to_S33_reaclib",  // 175,
        "S33_n_to_S34_reaclib",  // 176,
        "C12_C12_to_p_Na23_reaclib",  // 177,
        "Na23_p_to_He4_Ne20_reaclib",  // 178,
        "Na23_p_to_Mg24_reaclib",  // 179,
        "C12_C12_to_He4_Ne20_reaclib",  // 180,
        "Ne20_He4_to_Mg24_reaclib",  // 181,
        "Al27_to_Si27_beta_neg_weaktab",  // 182,
        "Co56_to_Fe56_electron_capture_weaktab",  // 183,
        "Co56_to_Ni56_beta_neg_weaktab",  // 184,
        "Fe56_to_Co56_beta_neg_weaktab",  // 185,
        "Ni56_to_Co56_electron_capture_weaktab",  // 186,
        "n_to_p_beta_neg_weaktab",  // 187,
        "p_to_n_electron_capture_weaktab",  // 188,
        "Si27_to_Al27_electron_capture_weaktab",  // 189,
        "Ni56_n_n_to_Ni58_approx",  // 190,
        "Ni58_to_Ni56_n_n_approx",  // 191,
        "Co53_n_n_to_Co55_approx",  // 192,
        "Co55_to_Co53_n_n_approx",  // 193,
        "Fe50_n_n_to_Fe52_approx",  // 194,
        "Fe52_to_Fe50_n_n_approx",  // 195,
        "Fe52_n_n_to_Fe54_approx",  // 196,
        "Fe54_to_Fe52_n_n_approx",  // 197,
        "Fe54_n_n_to_Fe56_approx",  // 198,
        "Fe56_to_Fe54_n_n_approx",  // 199,
        "Mn49_n_n_to_Mn51_approx",  // 200,
        "Mn51_to_Mn49_n_n_approx",  // 201,
        "Mn51_n_n_to_Mn53_approx",  // 202,
        "Mn53_to_Mn51_n_n_approx",  // 203,
        "Cr46_n_n_to_Cr48_approx",  // 204,
        "Cr48_to_Cr46_n_n_approx",  // 205,
        "Cr48_n_n_to_Cr50_approx",  // 206,
        "Cr50_to_Cr48_n_n_approx",  // 207,
        "Cr50_n_n_to_Cr52_approx",  // 208,
        "Cr52_to_Cr50_n_n_approx",  // 209,
        "V45_n_n_to_V47_approx",  // 210,
        "V47_to_V45_n_n_approx",  // 211,
        "V47_n_n_to_V49_approx",  // 212,
        "V49_to_V47_n_n_approx",  // 213,
        "Ti42_n_n_to_Ti44_approx",  // 214,
        "Ti44_to_Ti42_n_n_approx",  // 215,
        "Ti44_n_n_to_Ti46_approx",  // 216,
        "Ti46_to_Ti44_n_n_approx",  // 217,
        "Sc41_n_n_to_Sc43_approx",  // 218,
        "Sc43_to_Sc41_n_n_approx",  // 219,
        "Ca38_n_n_to_Ca40_approx",  // 220,
        "Ca40_to_Ca38_n_n_approx",  // 221,
        "Ca40_n_n_to_Ca42_approx",  // 222,
        "Ca42_to_Ca40_n_n_approx",  // 223,
        "Ar34_n_n_to_Ar36_approx",  // 224,
        "Ar36_to_Ar34_n_n_approx",  // 225,
        "Ar36_n_n_to_Ar38_approx",  // 226,
        "Ar38_to_Ar36_n_n_approx",  // 227,
        "Si28_n_n_to_Si30_approx",  // 228,
        "Si30_to_Si28_n_n_approx",  // 229,
        "S30_n_n_to_S32_approx",  // 230,
        "S32_to_S30_n_n_approx",  // 231,
        "S32_n_n_to_S34_approx",  // 232,
        "S34_to_S32_n_n_approx",  // 233,
        "C12_C12_to_Mg24_approx",  // 234,
        "Mg24_to_C12_C12_approx",  // 235,
        "C12_C12_to_Ne20_He4_approx",  // 236,
        "Ne20_He4_to_C12_C12_approx",  // 237,
        "Ne20_He4_to_Mg24_approx",  // 238,
        "Mg24_to_Ne20_He4_approx",  // 239,
        "O16_to_He4_C12_derived",  // 240,
        "Ne20_to_He4_O16_derived",  // 241,
        "Mg24_to_n_Mg23_derived",  // 242,
        "Si27_to_He4_Mg23_derived",  // 243,
        "Si28_to_n_Si27_derived",  // 244,
        "Si28_to_p_Al27_derived",  // 245,
        "Si28_to_He4_Mg24_derived",  // 246,
        "P29_to_p_Si28_derived",  // 247,
        "P31_to_p_Si30_derived",  // 248,
        "P31_to_He4_Al27_derived",  // 249,
        "S30_to_p_P29_derived",  // 250,
        "S32_to_p_P31_derived",  // 251,
        "S32_to_He4_Si28_derived",  // 252,
        "S34_to_He4_Si30_derived",  // 253,
        "Cl33_to_p_S32_derived",  // 254,
        "Cl33_to_He4_P29_derived",  // 255,
        "Cl35_to_p_S34_derived",  // 256,
        "Cl35_to_He4_P31_derived",  // 257,
        "Ar34_to_p_Cl33_derived",  // 258,
        "Ar34_to_He4_S30_derived",  // 259,
        "Ar36_to_p_Cl35_derived",  // 260,
        "Ar36_to_He4_S32_derived",  // 261,
        "Ar38_to_He4_S34_derived",  // 262,
        "K37_to_p_Ar36_derived",  // 263,
        "K37_to_He4_Cl33_derived",  // 264,
        "K39_to_p_Ar38_derived",  // 265,
        "K39_to_He4_Cl35_derived",  // 266,
        "Ca38_to_p_K37_derived",  // 267,
        "Ca38_to_He4_Ar34_derived",  // 268,
        "Ca40_to_p_K39_derived",  // 269,
        "Ca40_to_He4_Ar36_derived",  // 270,
        "Ca42_to_He4_Ar38_derived",  // 271,
        "Sc41_to_p_Ca40_derived",  // 272,
        "Sc41_to_He4_K37_derived",  // 273,
        "Sc43_to_p_Ca42_derived",  // 274,
        "Sc43_to_He4_K39_derived",  // 275,
        "Ti42_to_p_Sc41_derived",  // 276,
        "Ti42_to_He4_Ca38_derived",  // 277,
        "Ti44_to_p_Sc43_derived",  // 278,
        "Ti44_to_He4_Ca40_derived",  // 279,
        "Ti46_to_He4_Ca42_derived",  // 280,
        "V45_to_p_Ti44_derived",  // 281,
        "V45_to_He4_Sc41_derived",  // 282,
        "V47_to_p_Ti46_derived",  // 283,
        "V47_to_He4_Sc43_derived",  // 284,
        "Cr46_to_p_V45_derived",  // 285,
        "Cr46_to_He4_Ti42_derived",  // 286,
        "Cr48_to_p_V47_derived",  // 287,
        "Cr48_to_He4_Ti44_derived",  // 288,
        "Cr50_to_p_V49_derived",  // 289,
        "Cr50_to_He4_Ti46_derived",  // 290,
        "Mn49_to_p_Cr48_derived",  // 291,
        "Mn49_to_He4_V45_derived",  // 292,
        "Mn51_to_p_Cr50_derived",  // 293,
        "Mn51_to_He4_V47_derived",  // 294,
        "Mn53_to_p_Cr52_derived",  // 295,
        "Mn53_to_He4_V49_derived",  // 296,
        "Fe50_to_p_Mn49_derived",  // 297,
        "Fe50_to_He4_Cr46_derived",  // 298,
        "Fe52_to_p_Mn51_derived",  // 299,
        "Fe52_to_He4_Cr48_derived",  // 300,
        "Fe54_to_p_Mn53_derived",  // 301,
        "Fe54_to_He4_Cr50_derived",  // 302,
        "Fe56_to_He4_Cr52_derived",  // 303,
        "Co53_to_p_Fe52_derived",  // 304,
        "Co53_to_He4_Mn49_derived",  // 305,
        "Co55_to_p_Fe54_derived",  // 306,
        "Co55_to_He4_Mn51_derived",  // 307,
        "Co56_to_n_Co55_derived",  // 308,
        "Co57_to_n_Co56_derived",  // 309,
        "Co57_to_p_Fe56_derived",  // 310,
        "Co57_to_He4_Mn53_derived",  // 311,
        "Ni56_to_p_Co55_derived",  // 312,
        "Ni56_to_He4_Fe52_derived",  // 313,
        "Ni58_to_p_Co57_derived",  // 314,
        "Ni58_to_He4_Fe54_derived",  // 315,
        "Cu59_to_p_Ni58_derived",  // 316,
        "Cu59_to_He4_Co55_derived",  // 317,
        "Zn60_to_p_Cu59_derived",  // 318,
        "Zn60_to_He4_Ni56_derived",  // 319,
        "C12_to_3He4_derived",  // 320,
        "C12_C12_to_n_Mg23_derived",  // 321,
        "O16_C12_to_n_Si27_derived",  // 322,
        "Ne20_He4_to_n_Mg23_derived",  // 323,
        "Mg24_He4_to_n_Si27_derived",  // 324,
        "Mg24_He4_to_p_Al27_derived",  // 325,
        "Mg24_He4_to_C12_O16_derived",  // 326,
        "Al27_p_to_n_Si27_derived",  // 327,
        "Al27_p_to_C12_O16_derived",  // 328,
        "Si27_He4_to_n_S30_derived",  // 329,
        "Si28_He4_to_p_P31_derived",  // 330,
        "Si28_He4_to_C12_Ne20_derived",  // 331,
        "Si28_He4_to_O16_O16_derived",  // 332,
        "Si30_p_to_He4_Al27_derived",  // 333,
        "P31_p_to_C12_Ne20_derived",  // 334,
        "P31_p_to_O16_O16_derived",  // 335,
        "S32_p_to_He4_P29_derived",  // 336,
        "S32_He4_to_p_Cl35_derived",  // 337,
        "S34_p_to_He4_P31_derived",  // 338,
        "Cl33_p_to_He4_S30_derived",  // 339,
        "Ar36_p_to_He4_Cl33_derived",  // 340,
        "Ar36_He4_to_p_K39_derived",  // 341,
        "Ar38_p_to_He4_Cl35_derived",  // 342,
        "K37_p_to_He4_Ar34_derived",  // 343,
        "K39_He4_to_p_Ca42_derived",  // 344,
        "Ca40_p_to_He4_K37_derived",  // 345,
        "Ca40_He4_to_p_Sc43_derived",  // 346,
        "Sc41_p_to_He4_Ca38_derived",  // 347,
        "Ti44_p_to_He4_Sc41_derived",  // 348,
        "Ti46_p_to_He4_Sc43_derived",  // 349,
        "Ti46_He4_to_p_V49_derived",  // 350,
        "V45_p_to_He4_Ti42_derived",  // 351,
        "V47_p_to_He4_Ti44_derived",  // 352,
        "Cr48_p_to_He4_V45_derived",  // 353,
        "Cr50_p_to_He4_V47_derived",  // 354,
        "Cr50_He4_to_p_Mn53_derived",  // 355,
        "Cr52_p_to_He4_V49_derived",  // 356,
        "Mn49_p_to_He4_Cr46_derived",  // 357,
        "Mn51_p_to_He4_Cr48_derived",  // 358,
        "Mn53_He4_to_n_Co56_derived",  // 359,
        "Fe52_p_to_He4_Mn49_derived",  // 360,
        "Fe54_p_to_He4_Mn51_derived",  // 361,
        "Fe54_He4_to_p_Co57_derived",  // 362,
        "Fe56_p_to_n_Co56_derived",  // 363,
        "Fe56_p_to_He4_Mn53_derived",  // 364,
        "Co53_p_to_He4_Fe50_derived",  // 365,
        "Co55_p_to_He4_Fe52_derived",  // 366,
        "Co56_p_to_n_Ni56_derived",  // 367,
        "Co56_He4_to_n_Cu59_derived",  // 368,
        "Ni56_p_to_He4_Co53_derived",  // 369,
        "Ni56_He4_to_p_Cu59_derived",  // 370,
        "Ni58_p_to_He4_Co55_derived",  // 371,
        "Ni57_to_n_Ni56_derived",  // 372,
        "Ni58_to_n_Ni57_derived",  // 373,
        "Co54_to_n_Co53_derived",  // 374,
        "Co55_to_n_Co54_derived",  // 375,
        "Fe51_to_n_Fe50_derived",  // 376,
        "Fe52_to_n_Fe51_derived",  // 377,
        "Fe53_to_n_Fe52_derived",  // 378,
        "Fe54_to_n_Fe53_derived",  // 379,
        "Fe55_to_n_Fe54_derived",  // 380,
        "Fe56_to_n_Fe55_derived",  // 381,
        "Mn50_to_n_Mn49_derived",  // 382,
        "Mn51_to_n_Mn50_derived",  // 383,
        "Mn52_to_n_Mn51_derived",  // 384,
        "Mn53_to_n_Mn52_derived",  // 385,
        "Cr47_to_n_Cr46_derived",  // 386,
        "Cr48_to_n_Cr47_derived",  // 387,
        "Cr49_to_n_Cr48_derived",  // 388,
        "Cr50_to_n_Cr49_derived",  // 389,
        "Cr51_to_n_Cr50_derived",  // 390,
        "Cr52_to_n_Cr51_derived",  // 391,
        "V46_to_n_V45_derived",  // 392,
        "V47_to_n_V46_derived",  // 393,
        "V48_to_n_V47_derived",  // 394,
        "V49_to_n_V48_derived",  // 395,
        "Ti43_to_n_Ti42_derived",  // 396,
        "Ti44_to_n_Ti43_derived",  // 397,
        "Ti45_to_n_Ti44_derived",  // 398,
        "Ti46_to_n_Ti45_derived",  // 399,
        "Sc42_to_n_Sc41_derived",  // 400,
        "Sc43_to_n_Sc42_derived",  // 401,
        "Ca39_to_n_Ca38_derived",  // 402,
        "Ca40_to_n_Ca39_derived",  // 403,
        "Ca41_to_n_Ca40_derived",  // 404,
        "Ca42_to_n_Ca41_derived",  // 405,
        "Ar35_to_n_Ar34_derived",  // 406,
        "Ar36_to_n_Ar35_derived",  // 407,
        "Ar37_to_n_Ar36_derived",  // 408,
        "Ar38_to_n_Ar37_derived",  // 409,
        "Si29_to_n_Si28_derived",  // 410,
        "Si30_to_n_Si29_derived",  // 411,
        "S31_to_n_S30_derived",  // 412,
        "S32_to_n_S31_derived",  // 413,
        "S33_to_n_S32_derived",  // 414,
        "S34_to_n_S33_derived",  // 415,
        "Na23_p_to_C12_C12_derived",  // 416,
        "Mg24_to_p_Na23_derived",  // 417,
        "Ne20_He4_to_C12_C12_derived",  // 418,
        "Ne20_He4_to_p_Na23_derived",  // 419,
        "Mg24_to_He4_Ne20_derived"  // 420,
    };

    template<std::uint16_t rate>
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
    constexpr int NumNSERatePairs = 125;

    inline AMREX_GPU_MANAGED amrex::Array2D<short, 1, NumNSERatePairs, 1, 10, amrex::Order::C> rate_pair_data {
        He4, C12, -1, O16, -1, -1, C12, O16, k_C12_He4_to_O16_reaclib, k_O16_to_He4_C12_derived,
        He4, O16, -1, Ne20, -1, -1, O16, Ne20, k_O16_He4_to_Ne20_reaclib, k_Ne20_to_He4_O16_derived,
        N, Mg23, -1, Mg24, -1, -1, Mg23, Mg24, k_Mg23_n_to_Mg24_reaclib, k_Mg24_to_n_Mg23_derived,
        He4, Mg23, -1, Si27, -1, -1, Mg23, Si27, k_Mg23_He4_to_Si27_reaclib, k_Si27_to_He4_Mg23_derived,
        He4, Mg24, -1, Si28, -1, -1, Mg24, Si28, k_Mg24_He4_to_Si28_reaclib, k_Si28_to_He4_Mg24_derived,
        H1, Al27, -1, Si28, -1, -1, Al27, Si28, k_Al27_p_to_Si28_reaclib, k_Si28_to_p_Al27_derived,
        He4, Al27, -1, P31, -1, -1, Al27, P31, k_Al27_He4_to_P31_reaclib, k_P31_to_He4_Al27_derived,
        N, Si27, -1, Si28, -1, -1, Si27, Si28, k_Si27_n_to_Si28_reaclib, k_Si28_to_n_Si27_derived,
        H1, Si28, -1, P29, -1, -1, Si28, P29, k_Si28_p_to_P29_reaclib, k_P29_to_p_Si28_derived,
        He4, Si28, -1, S32, -1, -1, Si28, S32, k_Si28_He4_to_S32_reaclib, k_S32_to_He4_Si28_derived,
        H1, Si30, -1, P31, -1, -1, Si30, P31, k_Si30_p_to_P31_reaclib, k_P31_to_p_Si30_derived,
        He4, Si30, -1, S34, -1, -1, Si30, S34, k_Si30_He4_to_S34_reaclib, k_S34_to_He4_Si30_derived,
        H1, P29, -1, S30, -1, -1, P29, S30, k_P29_p_to_S30_reaclib, k_S30_to_p_P29_derived,
        He4, P29, -1, Cl33, -1, -1, P29, Cl33, k_P29_He4_to_Cl33_reaclib, k_Cl33_to_He4_P29_derived,
        H1, P31, -1, S32, -1, -1, P31, S32, k_P31_p_to_S32_reaclib, k_S32_to_p_P31_derived,
        He4, P31, -1, Cl35, -1, -1, P31, Cl35, k_P31_He4_to_Cl35_reaclib, k_Cl35_to_He4_P31_derived,
        He4, S30, -1, Ar34, -1, -1, S30, Ar34, k_S30_He4_to_Ar34_reaclib, k_Ar34_to_He4_S30_derived,
        H1, S32, -1, Cl33, -1, -1, S32, Cl33, k_S32_p_to_Cl33_reaclib, k_Cl33_to_p_S32_derived,
        He4, S32, -1, Ar36, -1, -1, S32, Ar36, k_S32_He4_to_Ar36_reaclib, k_Ar36_to_He4_S32_derived,
        H1, S34, -1, Cl35, -1, -1, S34, Cl35, k_S34_p_to_Cl35_reaclib, k_Cl35_to_p_S34_derived,
        He4, S34, -1, Ar38, -1, -1, S34, Ar38, k_S34_He4_to_Ar38_reaclib, k_Ar38_to_He4_S34_derived,
        H1, Cl33, -1, Ar34, -1, -1, Cl33, Ar34, k_Cl33_p_to_Ar34_reaclib, k_Ar34_to_p_Cl33_derived,
        He4, Cl33, -1, K37, -1, -1, Cl33, K37, k_Cl33_He4_to_K37_reaclib, k_K37_to_He4_Cl33_derived,
        H1, Cl35, -1, Ar36, -1, -1, Cl35, Ar36, k_Cl35_p_to_Ar36_reaclib, k_Ar36_to_p_Cl35_derived,
        He4, Cl35, -1, K39, -1, -1, Cl35, K39, k_Cl35_He4_to_K39_reaclib, k_K39_to_He4_Cl35_derived,
        He4, Ar34, -1, Ca38, -1, -1, Ar34, Ca38, k_Ar34_He4_to_Ca38_reaclib, k_Ca38_to_He4_Ar34_derived,
        H1, Ar36, -1, K37, -1, -1, Ar36, K37, k_Ar36_p_to_K37_reaclib, k_K37_to_p_Ar36_derived,
        He4, Ar36, -1, Ca40, -1, -1, Ar36, Ca40, k_Ar36_He4_to_Ca40_reaclib, k_Ca40_to_He4_Ar36_derived,
        H1, Ar38, -1, K39, -1, -1, Ar38, K39, k_Ar38_p_to_K39_reaclib, k_K39_to_p_Ar38_derived,
        He4, Ar38, -1, Ca42, -1, -1, Ar38, Ca42, k_Ar38_He4_to_Ca42_reaclib, k_Ca42_to_He4_Ar38_derived,
        H1, K37, -1, Ca38, -1, -1, K37, Ca38, k_K37_p_to_Ca38_reaclib, k_Ca38_to_p_K37_derived,
        He4, K37, -1, Sc41, -1, -1, K37, Sc41, k_K37_He4_to_Sc41_reaclib, k_Sc41_to_He4_K37_derived,
        H1, K39, -1, Ca40, -1, -1, K39, Ca40, k_K39_p_to_Ca40_reaclib, k_Ca40_to_p_K39_derived,
        He4, K39, -1, Sc43, -1, -1, K39, Sc43, k_K39_He4_to_Sc43_reaclib, k_Sc43_to_He4_K39_derived,
        He4, Ca38, -1, Ti42, -1, -1, Ca38, Ti42, k_Ca38_He4_to_Ti42_reaclib, k_Ti42_to_He4_Ca38_derived,
        H1, Ca40, -1, Sc41, -1, -1, Ca40, Sc41, k_Ca40_p_to_Sc41_reaclib, k_Sc41_to_p_Ca40_derived,
        He4, Ca40, -1, Ti44, -1, -1, Ca40, Ti44, k_Ca40_He4_to_Ti44_reaclib, k_Ti44_to_He4_Ca40_derived,
        H1, Ca42, -1, Sc43, -1, -1, Ca42, Sc43, k_Ca42_p_to_Sc43_reaclib, k_Sc43_to_p_Ca42_derived,
        He4, Ca42, -1, Ti46, -1, -1, Ca42, Ti46, k_Ca42_He4_to_Ti46_reaclib, k_Ti46_to_He4_Ca42_derived,
        H1, Sc41, -1, Ti42, -1, -1, Sc41, Ti42, k_Sc41_p_to_Ti42_reaclib, k_Ti42_to_p_Sc41_derived,
        He4, Sc41, -1, V45, -1, -1, Sc41, V45, k_Sc41_He4_to_V45_reaclib, k_V45_to_He4_Sc41_derived,
        H1, Sc43, -1, Ti44, -1, -1, Sc43, Ti44, k_Sc43_p_to_Ti44_reaclib, k_Ti44_to_p_Sc43_derived,
        He4, Sc43, -1, V47, -1, -1, Sc43, V47, k_Sc43_He4_to_V47_reaclib, k_V47_to_He4_Sc43_derived,
        He4, Ti42, -1, Cr46, -1, -1, Ti42, Cr46, k_Ti42_He4_to_Cr46_reaclib, k_Cr46_to_He4_Ti42_derived,
        H1, Ti44, -1, V45, -1, -1, Ti44, V45, k_Ti44_p_to_V45_reaclib, k_V45_to_p_Ti44_derived,
        He4, Ti44, -1, Cr48, -1, -1, Ti44, Cr48, k_Ti44_He4_to_Cr48_reaclib, k_Cr48_to_He4_Ti44_derived,
        H1, Ti46, -1, V47, -1, -1, Ti46, V47, k_Ti46_p_to_V47_reaclib, k_V47_to_p_Ti46_derived,
        He4, Ti46, -1, Cr50, -1, -1, Ti46, Cr50, k_Ti46_He4_to_Cr50_reaclib, k_Cr50_to_He4_Ti46_derived,
        H1, V45, -1, Cr46, -1, -1, V45, Cr46, k_V45_p_to_Cr46_reaclib, k_Cr46_to_p_V45_derived,
        He4, V45, -1, Mn49, -1, -1, V45, Mn49, k_V45_He4_to_Mn49_reaclib, k_Mn49_to_He4_V45_derived,
        H1, V47, -1, Cr48, -1, -1, V47, Cr48, k_V47_p_to_Cr48_reaclib, k_Cr48_to_p_V47_derived,
        He4, V47, -1, Mn51, -1, -1, V47, Mn51, k_V47_He4_to_Mn51_reaclib, k_Mn51_to_He4_V47_derived,
        H1, V49, -1, Cr50, -1, -1, V49, Cr50, k_V49_p_to_Cr50_reaclib, k_Cr50_to_p_V49_derived,
        He4, V49, -1, Mn53, -1, -1, V49, Mn53, k_V49_He4_to_Mn53_reaclib, k_Mn53_to_He4_V49_derived,
        He4, Cr46, -1, Fe50, -1, -1, Cr46, Fe50, k_Cr46_He4_to_Fe50_reaclib, k_Fe50_to_He4_Cr46_derived,
        H1, Cr48, -1, Mn49, -1, -1, Cr48, Mn49, k_Cr48_p_to_Mn49_reaclib, k_Mn49_to_p_Cr48_derived,
        He4, Cr48, -1, Fe52, -1, -1, Cr48, Fe52, k_Cr48_He4_to_Fe52_reaclib, k_Fe52_to_He4_Cr48_derived,
        H1, Cr50, -1, Mn51, -1, -1, Cr50, Mn51, k_Cr50_p_to_Mn51_reaclib, k_Mn51_to_p_Cr50_derived,
        He4, Cr50, -1, Fe54, -1, -1, Cr50, Fe54, k_Cr50_He4_to_Fe54_reaclib, k_Fe54_to_He4_Cr50_derived,
        H1, Cr52, -1, Mn53, -1, -1, Cr52, Mn53, k_Cr52_p_to_Mn53_reaclib, k_Mn53_to_p_Cr52_derived,
        He4, Cr52, -1, Fe56, -1, -1, Cr52, Fe56, k_Cr52_He4_to_Fe56_reaclib, k_Fe56_to_He4_Cr52_derived,
        H1, Mn49, -1, Fe50, -1, -1, Mn49, Fe50, k_Mn49_p_to_Fe50_reaclib, k_Fe50_to_p_Mn49_derived,
        He4, Mn49, -1, Co53, -1, -1, Mn49, Co53, k_Mn49_He4_to_Co53_reaclib, k_Co53_to_He4_Mn49_derived,
        H1, Mn51, -1, Fe52, -1, -1, Mn51, Fe52, k_Mn51_p_to_Fe52_reaclib, k_Fe52_to_p_Mn51_derived,
        He4, Mn51, -1, Co55, -1, -1, Mn51, Co55, k_Mn51_He4_to_Co55_reaclib, k_Co55_to_He4_Mn51_derived,
        H1, Mn53, -1, Fe54, -1, -1, Mn53, Fe54, k_Mn53_p_to_Fe54_reaclib, k_Fe54_to_p_Mn53_derived,
        He4, Mn53, -1, Co57, -1, -1, Mn53, Co57, k_Mn53_He4_to_Co57_reaclib, k_Co57_to_He4_Mn53_derived,
        H1, Fe52, -1, Co53, -1, -1, Fe52, Co53, k_Fe52_p_to_Co53_reaclib, k_Co53_to_p_Fe52_derived,
        He4, Fe52, -1, Ni56, -1, -1, Fe52, Ni56, k_Fe52_He4_to_Ni56_reaclib, k_Ni56_to_He4_Fe52_derived,
        H1, Fe54, -1, Co55, -1, -1, Fe54, Co55, k_Fe54_p_to_Co55_reaclib, k_Co55_to_p_Fe54_derived,
        He4, Fe54, -1, Ni58, -1, -1, Fe54, Ni58, k_Fe54_He4_to_Ni58_reaclib, k_Ni58_to_He4_Fe54_derived,
        H1, Fe56, -1, Co57, -1, -1, Fe56, Co57, k_Fe56_p_to_Co57_reaclib, k_Co57_to_p_Fe56_derived,
        N, Co55, -1, Co56, -1, -1, Co55, Co56, k_Co55_n_to_Co56_reaclib, k_Co56_to_n_Co55_derived,
        H1, Co55, -1, Ni56, -1, -1, Co55, Ni56, k_Co55_p_to_Ni56_reaclib, k_Ni56_to_p_Co55_derived,
        He4, Co55, -1, Cu59, -1, -1, Co55, Cu59, k_Co55_He4_to_Cu59_reaclib, k_Cu59_to_He4_Co55_derived,
        N, Co56, -1, Co57, -1, -1, Co56, Co57, k_Co56_n_to_Co57_reaclib, k_Co57_to_n_Co56_derived,
        H1, Co57, -1, Ni58, -1, -1, Co57, Ni58, k_Co57_p_to_Ni58_reaclib, k_Ni58_to_p_Co57_derived,
        He4, Ni56, -1, Zn60, -1, -1, Ni56, Zn60, k_Ni56_He4_to_Zn60_reaclib, k_Zn60_to_He4_Ni56_derived,
        H1, Ni58, -1, Cu59, -1, -1, Ni58, Cu59, k_Ni58_p_to_Cu59_reaclib, k_Cu59_to_p_Ni58_derived,
        H1, Cu59, -1, Zn60, -1, -1, Cu59, Zn60, k_Cu59_p_to_Zn60_reaclib, k_Zn60_to_p_Cu59_derived,
        N, Mg23, -1, He4, Ne20, -1, Mg23, Ne20, k_Mg23_n_to_He4_Ne20_reaclib, k_Ne20_He4_to_n_Mg23_derived,
        H1, Al27, -1, He4, Mg24, -1, Al27, Mg24, k_Al27_p_to_He4_Mg24_reaclib, k_Mg24_He4_to_p_Al27_derived,
        He4, Al27, -1, H1, Si30, -1, Al27, Si30, k_Al27_He4_to_p_Si30_reaclib, k_Si30_p_to_He4_Al27_derived,
        N, Si27, -1, H1, Al27, -1, Si27, Al27, k_Si27_n_to_p_Al27_reaclib, k_Al27_p_to_n_Si27_derived,
        N, Si27, -1, He4, Mg24, -1, Si27, Mg24, k_Si27_n_to_He4_Mg24_reaclib, k_Mg24_He4_to_n_Si27_derived,
        He4, P29, -1, H1, S32, -1, P29, S32, k_P29_He4_to_p_S32_reaclib, k_S32_p_to_He4_P29_derived,
        H1, P31, -1, He4, Si28, -1, P31, Si28, k_P31_p_to_He4_Si28_reaclib, k_Si28_He4_to_p_P31_derived,
        He4, P31, -1, H1, S34, -1, P31, S34, k_P31_He4_to_p_S34_reaclib, k_S34_p_to_He4_P31_derived,
        N, S30, -1, He4, Si27, -1, S30, Si27, k_S30_n_to_He4_Si27_reaclib, k_Si27_He4_to_n_S30_derived,
        He4, S30, -1, H1, Cl33, -1, S30, Cl33, k_S30_He4_to_p_Cl33_reaclib, k_Cl33_p_to_He4_S30_derived,
        He4, Cl33, -1, H1, Ar36, -1, Cl33, Ar36, k_Cl33_He4_to_p_Ar36_reaclib, k_Ar36_p_to_He4_Cl33_derived,
        H1, Cl35, -1, He4, S32, -1, Cl35, S32, k_Cl35_p_to_He4_S32_reaclib, k_S32_He4_to_p_Cl35_derived,
        He4, Cl35, -1, H1, Ar38, -1, Cl35, Ar38, k_Cl35_He4_to_p_Ar38_reaclib, k_Ar38_p_to_He4_Cl35_derived,
        He4, Ar34, -1, H1, K37, -1, Ar34, K37, k_Ar34_He4_to_p_K37_reaclib, k_K37_p_to_He4_Ar34_derived,
        He4, K37, -1, H1, Ca40, -1, K37, Ca40, k_K37_He4_to_p_Ca40_reaclib, k_Ca40_p_to_He4_K37_derived,
        H1, K39, -1, He4, Ar36, -1, K39, Ar36, k_K39_p_to_He4_Ar36_reaclib, k_Ar36_He4_to_p_K39_derived,
        He4, Ca38, -1, H1, Sc41, -1, Ca38, Sc41, k_Ca38_He4_to_p_Sc41_reaclib, k_Sc41_p_to_He4_Ca38_derived,
        H1, Ca42, -1, He4, K39, -1, Ca42, K39, k_Ca42_p_to_He4_K39_reaclib, k_K39_He4_to_p_Ca42_derived,
        He4, Sc41, -1, H1, Ti44, -1, Sc41, Ti44, k_Sc41_He4_to_p_Ti44_reaclib, k_Ti44_p_to_He4_Sc41_derived,
        H1, Sc43, -1, He4, Ca40, -1, Sc43, Ca40, k_Sc43_p_to_He4_Ca40_reaclib, k_Ca40_He4_to_p_Sc43_derived,
        He4, Sc43, -1, H1, Ti46, -1, Sc43, Ti46, k_Sc43_He4_to_p_Ti46_reaclib, k_Ti46_p_to_He4_Sc43_derived,
        He4, Ti42, -1, H1, V45, -1, Ti42, V45, k_Ti42_He4_to_p_V45_reaclib, k_V45_p_to_He4_Ti42_derived,
        He4, V45, -1, H1, Cr48, -1, V45, Cr48, k_V45_He4_to_p_Cr48_reaclib, k_Cr48_p_to_He4_V45_derived,
        He4, V47, -1, H1, Cr50, -1, V47, Cr50, k_V47_He4_to_p_Cr50_reaclib, k_Cr50_p_to_He4_V47_derived,
        H1, V49, -1, He4, Ti46, -1, V49, Ti46, k_V49_p_to_He4_Ti46_reaclib, k_Ti46_He4_to_p_V49_derived,
        He4, V49, -1, H1, Cr52, -1, V49, Cr52, k_V49_He4_to_p_Cr52_reaclib, k_Cr52_p_to_He4_V49_derived,
        He4, Cr46, -1, H1, Mn49, -1, Cr46, Mn49, k_Cr46_He4_to_p_Mn49_reaclib, k_Mn49_p_to_He4_Cr46_derived,
        He4, Cr48, -1, H1, Mn51, -1, Cr48, Mn51, k_Cr48_He4_to_p_Mn51_reaclib, k_Mn51_p_to_He4_Cr48_derived,
        He4, Mn49, -1, H1, Fe52, -1, Mn49, Fe52, k_Mn49_He4_to_p_Fe52_reaclib, k_Fe52_p_to_He4_Mn49_derived,
        He4, Mn51, -1, H1, Fe54, -1, Mn51, Fe54, k_Mn51_He4_to_p_Fe54_reaclib, k_Fe54_p_to_He4_Mn51_derived,
        H1, Mn53, -1, He4, Cr50, -1, Mn53, Cr50, k_Mn53_p_to_He4_Cr50_reaclib, k_Cr50_He4_to_p_Mn53_derived,
        He4, Mn53, -1, H1, Fe56, -1, Mn53, Fe56, k_Mn53_He4_to_p_Fe56_reaclib, k_Fe56_p_to_He4_Mn53_derived,
        He4, Fe50, -1, H1, Co53, -1, Fe50, Co53, k_Fe50_He4_to_p_Co53_reaclib, k_Co53_p_to_He4_Fe50_derived,
        He4, Fe52, -1, H1, Co55, -1, Fe52, Co55, k_Fe52_He4_to_p_Co55_reaclib, k_Co55_p_to_He4_Fe52_derived,
        He4, Co53, -1, H1, Ni56, -1, Co53, Ni56, k_Co53_He4_to_p_Ni56_reaclib, k_Ni56_p_to_He4_Co53_derived,
        He4, Co55, -1, H1, Ni58, -1, Co55, Ni58, k_Co55_He4_to_p_Ni58_reaclib, k_Ni58_p_to_He4_Co55_derived,
        N, Co56, -1, H1, Fe56, -1, Co56, Fe56, k_Co56_n_to_p_Fe56_reaclib, k_Fe56_p_to_n_Co56_derived,
        N, Co56, -1, He4, Mn53, -1, Co56, Mn53, k_Co56_n_to_He4_Mn53_reaclib, k_Mn53_He4_to_n_Co56_derived,
        H1, Co57, -1, He4, Fe54, -1, Co57, Fe54, k_Co57_p_to_He4_Fe54_reaclib, k_Fe54_He4_to_p_Co57_derived,
        N, Ni56, -1, H1, Co56, -1, Ni56, Co56, k_Ni56_n_to_p_Co56_reaclib, k_Co56_p_to_n_Ni56_derived,
        N, Cu59, -1, He4, Co56, -1, Cu59, Co56, k_Cu59_n_to_He4_Co56_reaclib, k_Co56_He4_to_n_Cu59_derived,
        H1, Cu59, -1, He4, Ni56, -1, Cu59, Ni56, k_Cu59_p_to_He4_Ni56_reaclib, k_Ni56_He4_to_p_Cu59_derived,
        He4, He4, He4, C12, -1, -1, C12, -1, k_He4_He4_He4_to_C12_reaclib, k_C12_to_3He4_derived,
        H1, V47, -1, He4, Ti44, -1, V47, Ti44, k_V47_p_to_He4_Ti44_derived, k_Ti44_He4_to_p_V47_reaclib,
        Ne20, He4, -1, Mg24, -1, -1, Ne20, Mg24, k_Ne20_He4_to_Mg24_approx, k_Mg24_to_Ne20_He4_approx
    };
}
#endif

#endif
