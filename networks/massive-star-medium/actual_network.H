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
        else if constexpr (spec == Si29) {
            return 245.01044789599834_rt;
        }
        else if constexpr (spec == Si30) {
            return 255.61964649599756_rt;
        }
        else if constexpr (spec == Si31) {
            return 262.207042596001_rt;
        }
        else if constexpr (spec == P29) {
            return 239.28581936000046_rt;
        }
        else if constexpr (spec == P31) {
            return 262.9161999600037_rt;
        }
        else if constexpr (spec == P33) {
            return 280.95569175999844_rt;
        }
        else if constexpr (spec == S32) {
            return 271.78016372399725_rt;
        }
        else if constexpr (spec == S33) {
            return 280.42180302400084_rt;
        }
        else if constexpr (spec == S34) {
            return 291.83895282400044_rt;
        }
        else if constexpr (spec == S35) {
            return 298.8247909240017_rt;
        }
        else if constexpr (spec == Cl33) {
            return 274.05689768799857_rt;
        }
        else if constexpr (spec == Cl35) {
            return 298.20976388800045_rt;
        }
        else if constexpr (spec == Cl37) {
            return 317.10042008799786_rt;
        }
        else if constexpr (spec == Ar36) {
            return 306.7167469519991_rt;
        }
        else if constexpr (spec == Ar37) {
            return 315.50420305199077_rt;
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
        else if constexpr (spec == Sc45) {
            return 387.8523267439959_rt;
        }
        else if constexpr (spec == Ti44) {
            return 375.47496160800074_rt;
        }
        else if constexpr (spec == Ti46) {
            return 398.19726780800556_rt;
        }
        else if constexpr (spec == Ti48) {
            return 418.70458400799544_rt;
        }
        else if constexpr (spec == V47) {
            return 403.3650388719979_rt;
        }
        else if constexpr (spec == V49) {
            return 425.4628050719984_rt;
        }
        else if constexpr (spec == V51) {
            return 445.84635127200454_rt;
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
        else if constexpr (spec == Cr54) {
            return 474.01022853599716_rt;
        }
        else if constexpr (spec == Mn51) {
            return 440.321747199996_rt;
        }
        else if constexpr (spec == Mn53) {
            return 462.9114833999993_rt;
        }
        else if constexpr (spec == Mn55) {
            return 482.07635959999607_rt;
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
        else if constexpr (spec == Fe58) {
            return 509.9507268639936_rt;
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
        else if constexpr (spec == Co59) {
            return 517.3141979280044_rt;
        }
        else if constexpr (spec == Ni56) {
            return 483.9956965919919_rt;
        }
        else if constexpr (spec == Ni58) {
            return 506.4596327920008_rt;
        }
        else if constexpr (spec == Ni60) {
            return 526.8465689919976_rt;
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
        else if constexpr (spec == Si29) {
            return 4.811660140333913e-23_rt;
        }
        else if constexpr (spec == Si30) {
            return 4.977261628701015e-23_rt;
        }
        else if constexpr (spec == Si31) {
            return 5.143580068485421e-23_rt;
        }
        else if constexpr (spec == P29) {
            return 4.812541182037512e-23_rt;
        }
        else if constexpr (spec == P31) {
            return 5.14331418367544e-23_rt;
        }
        else if constexpr (spec == P33) {
            return 5.475083851770546e-23_rt;
        }
        else if constexpr (spec == S32) {
            return 5.309087322384128e-23_rt;
        }
        else if constexpr (spec == S33) {
            return 5.475039560055494e-23_rt;
        }
        else if constexpr (spec == S34) {
            return 5.64049701803872e-23_rt;
        }
        else if constexpr (spec == S35) {
            return 5.806744429084988e-23_rt;
        }
        else if constexpr (spec == Cl33) {
            return 5.476034741466158e-23_rt;
        }
        else if constexpr (spec == Cl35) {
            return 5.806714601585715e-23_rt;
        }
        else if constexpr (spec == Cl37) {
            return 6.13833253584433e-23_rt;
        }
        else if constexpr (spec == Ar36) {
            return 5.972551377884467e-23_rt;
        }
        else if constexpr (spec == Ar37) {
            return 6.138477621350144e-23_rt;
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
        else if constexpr (spec == Sc45) {
            return 7.465103997175656e-23_rt;
        }
        else if constexpr (spec == Ti44) {
            return 7.299678247096977e-23_rt;
        }
        else if constexpr (spec == Ti46) {
            return 7.630613127700392e-23_rt;
        }
        else if constexpr (spec == Ti48) {
            return 7.961942866136786e-23_rt;
        }
        else if constexpr (spec == V47) {
            return 7.797045172607207e-23_rt;
        }
        else if constexpr (spec == V49) {
            return 8.128091387578276e-23_rt;
        }
        else if constexpr (spec == V51) {
            return 8.459443190021272e-23_rt;
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
        else if constexpr (spec == Cr54) {
            return 8.956761306248089e-23_rt;
        }
        else if constexpr (spec == Mn51) {
            return 8.460149108098279e-23_rt;
        }
        else if constexpr (spec == Mn53) {
            return 8.79110762145079e-23_rt;
        }
        else if constexpr (spec == Mn55) {
            return 9.122676671554189e-23_rt;
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
        else if constexpr (spec == Fe58) {
            return 9.620046397626297e-23_rt;
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
        else if constexpr (spec == Co59) {
            return 9.786087023454982e-23_rt;
        }
        else if constexpr (spec == Ni56) {
            return 9.289408870379396e-23_rt;
        }
        else if constexpr (spec == Ni58) {
            return 9.620389809618879e-23_rt;
        }
        else if constexpr (spec == Ni60) {
            return 9.951741007739486e-23_rt;
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
        k_O16_to_He4_C12_reaclib = 1,
        k_Ne20_to_He4_O16_reaclib = 2,
        k_Mg24_to_n_Mg23_reaclib = 3,
        k_Si27_to_He4_Mg23_reaclib = 4,
        k_Si28_to_n_Si27_reaclib = 5,
        k_Si28_to_p_Al27_reaclib = 6,
        k_Si28_to_He4_Mg24_reaclib = 7,
        k_Si29_to_n_Si28_reaclib = 8,
        k_Si30_to_n_Si29_reaclib = 9,
        k_Si31_to_n_Si30_reaclib = 10,
        k_P29_to_p_Si28_reaclib = 11,
        k_P31_to_p_Si30_reaclib = 12,
        k_P31_to_He4_Al27_reaclib = 13,
        k_S32_to_p_P31_reaclib = 14,
        k_S32_to_He4_Si28_reaclib = 15,
        k_S33_to_n_S32_reaclib = 16,
        k_S33_to_He4_Si29_reaclib = 17,
        k_S34_to_n_S33_reaclib = 18,
        k_S34_to_p_P33_reaclib = 19,
        k_S34_to_He4_Si30_reaclib = 20,
        k_S35_to_n_S34_reaclib = 21,
        k_S35_to_He4_Si31_reaclib = 22,
        k_Cl33_to_p_S32_reaclib = 23,
        k_Cl33_to_He4_P29_reaclib = 24,
        k_Cl35_to_p_S34_reaclib = 25,
        k_Cl35_to_He4_P31_reaclib = 26,
        k_Cl37_to_He4_P33_reaclib = 27,
        k_Ar36_to_p_Cl35_reaclib = 28,
        k_Ar36_to_He4_S32_reaclib = 29,
        k_Ar37_to_n_Ar36_reaclib = 30,
        k_Ar37_to_He4_S33_reaclib = 31,
        k_Ar38_to_n_Ar37_reaclib = 32,
        k_Ar38_to_p_Cl37_reaclib = 33,
        k_Ar38_to_He4_S34_reaclib = 34,
        k_K37_to_p_Ar36_reaclib = 35,
        k_K37_to_He4_Cl33_reaclib = 36,
        k_K39_to_p_Ar38_reaclib = 37,
        k_K39_to_He4_Cl35_reaclib = 38,
        k_Ca40_to_p_K39_reaclib = 39,
        k_Ca40_to_He4_Ar36_reaclib = 40,
        k_Ca42_to_He4_Ar38_reaclib = 41,
        k_Sc41_to_p_Ca40_reaclib = 42,
        k_Sc41_to_He4_K37_reaclib = 43,
        k_Sc43_to_p_Ca42_reaclib = 44,
        k_Sc43_to_He4_K39_reaclib = 45,
        k_Ti44_to_p_Sc43_reaclib = 46,
        k_Ti44_to_He4_Ca40_reaclib = 47,
        k_Ti46_to_p_Sc45_reaclib = 48,
        k_Ti46_to_He4_Ca42_reaclib = 49,
        k_V47_to_p_Ti46_reaclib = 50,
        k_V47_to_He4_Sc43_reaclib = 51,
        k_V49_to_p_Ti48_reaclib = 52,
        k_V49_to_He4_Sc45_reaclib = 53,
        k_Cr48_to_p_V47_reaclib = 54,
        k_Cr48_to_He4_Ti44_reaclib = 55,
        k_Cr50_to_p_V49_reaclib = 56,
        k_Cr50_to_He4_Ti46_reaclib = 57,
        k_Cr52_to_p_V51_reaclib = 58,
        k_Cr52_to_He4_Ti48_reaclib = 59,
        k_Mn51_to_p_Cr50_reaclib = 60,
        k_Mn51_to_He4_V47_reaclib = 61,
        k_Mn53_to_p_Cr52_reaclib = 62,
        k_Mn53_to_He4_V49_reaclib = 63,
        k_Mn55_to_p_Cr54_reaclib = 64,
        k_Mn55_to_He4_V51_reaclib = 65,
        k_Fe52_to_p_Mn51_reaclib = 66,
        k_Fe52_to_He4_Cr48_reaclib = 67,
        k_Fe54_to_p_Mn53_reaclib = 68,
        k_Fe54_to_He4_Cr50_reaclib = 69,
        k_Fe56_to_p_Mn55_reaclib = 70,
        k_Fe56_to_He4_Cr52_reaclib = 71,
        k_Fe58_to_He4_Cr54_reaclib = 72,
        k_Co55_to_p_Fe54_reaclib = 73,
        k_Co55_to_He4_Mn51_reaclib = 74,
        k_Co56_to_n_Co55_reaclib = 75,
        k_Co57_to_n_Co56_reaclib = 76,
        k_Co57_to_p_Fe56_reaclib = 77,
        k_Co57_to_He4_Mn53_reaclib = 78,
        k_Co59_to_p_Fe58_reaclib = 79,
        k_Co59_to_He4_Mn55_reaclib = 80,
        k_Ni56_to_p_Co55_reaclib = 81,
        k_Ni56_to_He4_Fe52_reaclib = 82,
        k_Ni58_to_p_Co57_reaclib = 83,
        k_Ni58_to_He4_Fe54_reaclib = 84,
        k_Ni60_to_p_Co59_reaclib = 85,
        k_Ni60_to_He4_Fe56_reaclib = 86,
        k_Cu59_to_p_Ni58_reaclib = 87,
        k_Cu59_to_He4_Co55_reaclib = 88,
        k_Zn60_to_p_Cu59_reaclib = 89,
        k_Zn60_to_He4_Ni56_reaclib = 90,
        k_C12_to_3He4_reaclib = 91,
        k_C12_He4_to_O16_reaclib = 92,
        k_O16_He4_to_Ne20_reaclib = 93,
        k_Mg23_n_to_Mg24_reaclib = 94,
        k_Mg23_He4_to_Si27_reaclib = 95,
        k_Mg24_He4_to_Si28_reaclib = 96,
        k_Al27_p_to_Si28_reaclib = 97,
        k_Al27_He4_to_P31_reaclib = 98,
        k_Si27_n_to_Si28_reaclib = 99,
        k_Si28_n_to_Si29_reaclib = 100,
        k_Si28_p_to_P29_reaclib = 101,
        k_Si28_He4_to_S32_reaclib = 102,
        k_Si29_n_to_Si30_reaclib = 103,
        k_Si29_He4_to_S33_reaclib = 104,
        k_Si30_n_to_Si31_reaclib = 105,
        k_Si30_p_to_P31_reaclib = 106,
        k_Si30_He4_to_S34_reaclib = 107,
        k_Si31_He4_to_S35_reaclib = 108,
        k_P29_He4_to_Cl33_reaclib = 109,
        k_P31_p_to_S32_reaclib = 110,
        k_P31_He4_to_Cl35_reaclib = 111,
        k_P33_p_to_S34_reaclib = 112,
        k_P33_He4_to_Cl37_reaclib = 113,
        k_S32_n_to_S33_reaclib = 114,
        k_S32_p_to_Cl33_reaclib = 115,
        k_S32_He4_to_Ar36_reaclib = 116,
        k_S33_n_to_S34_reaclib = 117,
        k_S33_He4_to_Ar37_reaclib = 118,
        k_S34_n_to_S35_reaclib = 119,
        k_S34_p_to_Cl35_reaclib = 120,
        k_S34_He4_to_Ar38_reaclib = 121,
        k_Cl33_He4_to_K37_reaclib = 122,
        k_Cl35_p_to_Ar36_reaclib = 123,
        k_Cl35_He4_to_K39_reaclib = 124,
        k_Cl37_p_to_Ar38_reaclib = 125,
        k_Ar36_n_to_Ar37_reaclib = 126,
        k_Ar36_p_to_K37_reaclib = 127,
        k_Ar36_He4_to_Ca40_reaclib = 128,
        k_Ar37_n_to_Ar38_reaclib = 129,
        k_Ar38_p_to_K39_reaclib = 130,
        k_Ar38_He4_to_Ca42_reaclib = 131,
        k_K37_He4_to_Sc41_reaclib = 132,
        k_K39_p_to_Ca40_reaclib = 133,
        k_K39_He4_to_Sc43_reaclib = 134,
        k_Ca40_p_to_Sc41_reaclib = 135,
        k_Ca40_He4_to_Ti44_reaclib = 136,
        k_Ca42_p_to_Sc43_reaclib = 137,
        k_Ca42_He4_to_Ti46_reaclib = 138,
        k_Sc43_p_to_Ti44_reaclib = 139,
        k_Sc43_He4_to_V47_reaclib = 140,
        k_Sc45_p_to_Ti46_reaclib = 141,
        k_Sc45_He4_to_V49_reaclib = 142,
        k_Ti44_He4_to_Cr48_reaclib = 143,
        k_Ti46_p_to_V47_reaclib = 144,
        k_Ti46_He4_to_Cr50_reaclib = 145,
        k_Ti48_p_to_V49_reaclib = 146,
        k_Ti48_He4_to_Cr52_reaclib = 147,
        k_V47_p_to_Cr48_reaclib = 148,
        k_V47_He4_to_Mn51_reaclib = 149,
        k_V49_p_to_Cr50_reaclib = 150,
        k_V49_He4_to_Mn53_reaclib = 151,
        k_V51_p_to_Cr52_reaclib = 152,
        k_V51_He4_to_Mn55_reaclib = 153,
        k_Cr48_He4_to_Fe52_reaclib = 154,
        k_Cr50_p_to_Mn51_reaclib = 155,
        k_Cr50_He4_to_Fe54_reaclib = 156,
        k_Cr52_p_to_Mn53_reaclib = 157,
        k_Cr52_He4_to_Fe56_reaclib = 158,
        k_Cr54_p_to_Mn55_reaclib = 159,
        k_Cr54_He4_to_Fe58_reaclib = 160,
        k_Mn51_p_to_Fe52_reaclib = 161,
        k_Mn51_He4_to_Co55_reaclib = 162,
        k_Mn53_p_to_Fe54_reaclib = 163,
        k_Mn53_He4_to_Co57_reaclib = 164,
        k_Mn55_p_to_Fe56_reaclib = 165,
        k_Mn55_He4_to_Co59_reaclib = 166,
        k_Fe52_He4_to_Ni56_reaclib = 167,
        k_Fe54_p_to_Co55_reaclib = 168,
        k_Fe54_He4_to_Ni58_reaclib = 169,
        k_Fe56_p_to_Co57_reaclib = 170,
        k_Fe56_He4_to_Ni60_reaclib = 171,
        k_Fe58_p_to_Co59_reaclib = 172,
        k_Co55_n_to_Co56_reaclib = 173,
        k_Co55_p_to_Ni56_reaclib = 174,
        k_Co55_He4_to_Cu59_reaclib = 175,
        k_Co56_n_to_Co57_reaclib = 176,
        k_Co57_p_to_Ni58_reaclib = 177,
        k_Co59_p_to_Ni60_reaclib = 178,
        k_Ni56_He4_to_Zn60_reaclib = 179,
        k_Ni58_p_to_Cu59_reaclib = 180,
        k_Cu59_p_to_Zn60_reaclib = 181,
        k_C12_C12_to_n_Mg23_reaclib = 182,
        k_O16_C12_to_n_Si27_reaclib = 183,
        k_O16_C12_to_p_Al27_reaclib = 184,
        k_O16_C12_to_He4_Mg24_reaclib = 185,
        k_O16_O16_to_p_P31_reaclib = 186,
        k_O16_O16_to_He4_Si28_reaclib = 187,
        k_Ne20_He4_to_n_Mg23_reaclib = 188,
        k_Ne20_C12_to_p_P31_reaclib = 189,
        k_Ne20_C12_to_He4_Si28_reaclib = 190,
        k_Mg23_n_to_He4_Ne20_reaclib = 191,
        k_Mg23_n_to_C12_C12_reaclib = 192,
        k_Mg24_He4_to_n_Si27_reaclib = 193,
        k_Mg24_He4_to_p_Al27_reaclib = 194,
        k_Mg24_He4_to_C12_O16_reaclib = 195,
        k_Al27_p_to_n_Si27_reaclib = 196,
        k_Al27_p_to_He4_Mg24_reaclib = 197,
        k_Al27_p_to_C12_O16_reaclib = 198,
        k_Al27_He4_to_p_Si30_reaclib = 199,
        k_Si27_n_to_p_Al27_reaclib = 200,
        k_Si27_n_to_He4_Mg24_reaclib = 201,
        k_Si27_n_to_C12_O16_reaclib = 202,
        k_Si28_He4_to_p_P31_reaclib = 203,
        k_Si28_He4_to_C12_Ne20_reaclib = 204,
        k_Si28_He4_to_O16_O16_reaclib = 205,
        k_Si29_p_to_n_P29_reaclib = 206,
        k_Si29_He4_to_n_S32_reaclib = 207,
        k_Si30_p_to_He4_Al27_reaclib = 208,
        k_Si30_He4_to_n_S33_reaclib = 209,
        k_Si30_He4_to_p_P33_reaclib = 210,
        k_Si31_p_to_n_P31_reaclib = 211,
        k_Si31_He4_to_n_S34_reaclib = 212,
        k_P29_n_to_p_Si29_reaclib = 213,
        k_P29_He4_to_p_S32_reaclib = 214,
        k_P31_n_to_p_Si31_reaclib = 215,
        k_P31_p_to_He4_Si28_reaclib = 216,
        k_P31_p_to_C12_Ne20_reaclib = 217,
        k_P31_p_to_O16_O16_reaclib = 218,
        k_P31_He4_to_p_S34_reaclib = 219,
        k_P33_p_to_n_S33_reaclib = 220,
        k_P33_p_to_He4_Si30_reaclib = 221,
        k_S32_n_to_He4_Si29_reaclib = 222,
        k_S32_p_to_He4_P29_reaclib = 223,
        k_S32_He4_to_p_Cl35_reaclib = 224,
        k_S33_n_to_p_P33_reaclib = 225,
        k_S33_n_to_He4_Si30_reaclib = 226,
        k_S33_p_to_n_Cl33_reaclib = 227,
        k_S33_He4_to_n_Ar36_reaclib = 228,
        k_S34_n_to_He4_Si31_reaclib = 229,
        k_S34_p_to_He4_P31_reaclib = 230,
        k_S34_He4_to_n_Ar37_reaclib = 231,
        k_S34_He4_to_p_Cl37_reaclib = 232,
        k_S35_p_to_n_Cl35_reaclib = 233,
        k_S35_He4_to_n_Ar38_reaclib = 234,
        k_Cl33_n_to_p_S33_reaclib = 235,
        k_Cl33_He4_to_p_Ar36_reaclib = 236,
        k_Cl35_n_to_p_S35_reaclib = 237,
        k_Cl35_p_to_He4_S32_reaclib = 238,
        k_Cl35_He4_to_p_Ar38_reaclib = 239,
        k_Cl37_p_to_n_Ar37_reaclib = 240,
        k_Cl37_p_to_He4_S34_reaclib = 241,
        k_Ar36_n_to_He4_S33_reaclib = 242,
        k_Ar36_p_to_He4_Cl33_reaclib = 243,
        k_Ar36_He4_to_p_K39_reaclib = 244,
        k_Ar37_n_to_p_Cl37_reaclib = 245,
        k_Ar37_n_to_He4_S34_reaclib = 246,
        k_Ar37_p_to_n_K37_reaclib = 247,
        k_Ar37_He4_to_n_Ca40_reaclib = 248,
        k_Ar38_n_to_He4_S35_reaclib = 249,
        k_Ar38_p_to_He4_Cl35_reaclib = 250,
        k_K37_n_to_p_Ar37_reaclib = 251,
        k_K37_He4_to_p_Ca40_reaclib = 252,
        k_K39_p_to_He4_Ar36_reaclib = 253,
        k_K39_He4_to_p_Ca42_reaclib = 254,
        k_Ca40_n_to_He4_Ar37_reaclib = 255,
        k_Ca40_p_to_He4_K37_reaclib = 256,
        k_Ca40_He4_to_p_Sc43_reaclib = 257,
        k_Ca42_p_to_He4_K39_reaclib = 258,
        k_Ca42_He4_to_p_Sc45_reaclib = 259,
        k_Sc41_He4_to_p_Ti44_reaclib = 260,
        k_Sc43_p_to_He4_Ca40_reaclib = 261,
        k_Sc43_He4_to_p_Ti46_reaclib = 262,
        k_Sc45_p_to_He4_Ca42_reaclib = 263,
        k_Sc45_He4_to_p_Ti48_reaclib = 264,
        k_Ti44_p_to_He4_Sc41_reaclib = 265,
        k_Ti44_He4_to_p_V47_reaclib = 266,
        k_Ti46_p_to_He4_Sc43_reaclib = 267,
        k_Ti46_He4_to_p_V49_reaclib = 268,
        k_Ti48_p_to_He4_Sc45_reaclib = 269,
        k_Ti48_He4_to_p_V51_reaclib = 270,
        k_V47_p_to_He4_Ti44_reaclib = 271,
        k_V47_He4_to_p_Cr50_reaclib = 272,
        k_V49_p_to_He4_Ti46_reaclib = 273,
        k_V49_He4_to_p_Cr52_reaclib = 274,
        k_V51_p_to_He4_Ti48_reaclib = 275,
        k_V51_He4_to_p_Cr54_reaclib = 276,
        k_Cr48_He4_to_p_Mn51_reaclib = 277,
        k_Cr50_p_to_He4_V47_reaclib = 278,
        k_Cr50_He4_to_p_Mn53_reaclib = 279,
        k_Cr52_p_to_He4_V49_reaclib = 280,
        k_Cr52_He4_to_p_Mn55_reaclib = 281,
        k_Cr54_p_to_He4_V51_reaclib = 282,
        k_Mn51_p_to_He4_Cr48_reaclib = 283,
        k_Mn51_He4_to_p_Fe54_reaclib = 284,
        k_Mn53_p_to_He4_Cr50_reaclib = 285,
        k_Mn53_He4_to_n_Co56_reaclib = 286,
        k_Mn53_He4_to_p_Fe56_reaclib = 287,
        k_Mn55_p_to_He4_Cr52_reaclib = 288,
        k_Mn55_He4_to_p_Fe58_reaclib = 289,
        k_Fe52_He4_to_p_Co55_reaclib = 290,
        k_Fe54_p_to_He4_Mn51_reaclib = 291,
        k_Fe54_He4_to_p_Co57_reaclib = 292,
        k_Fe56_p_to_n_Co56_reaclib = 293,
        k_Fe56_p_to_He4_Mn53_reaclib = 294,
        k_Fe56_He4_to_p_Co59_reaclib = 295,
        k_Fe58_p_to_He4_Mn55_reaclib = 296,
        k_Co55_p_to_He4_Fe52_reaclib = 297,
        k_Co55_He4_to_p_Ni58_reaclib = 298,
        k_Co56_n_to_p_Fe56_reaclib = 299,
        k_Co56_n_to_He4_Mn53_reaclib = 300,
        k_Co56_p_to_n_Ni56_reaclib = 301,
        k_Co56_He4_to_n_Cu59_reaclib = 302,
        k_Co57_p_to_He4_Fe54_reaclib = 303,
        k_Co57_He4_to_p_Ni60_reaclib = 304,
        k_Co59_p_to_He4_Fe56_reaclib = 305,
        k_Ni56_n_to_p_Co56_reaclib = 306,
        k_Ni56_He4_to_p_Cu59_reaclib = 307,
        k_Ni58_p_to_He4_Co55_reaclib = 308,
        k_Ni60_p_to_He4_Co57_reaclib = 309,
        k_Cu59_n_to_He4_Co56_reaclib = 310,
        k_Cu59_p_to_He4_Ni56_reaclib = 311,
        k_He4_He4_He4_to_C12_reaclib = 312,
        k_Ni56_n_to_Ni57_reaclib = 313,
        k_Ni57_to_n_Ni56_reaclib = 314,
        k_Ni57_n_to_Ni58_reaclib = 315,
        k_Ni58_to_n_Ni57_reaclib = 316,
        k_Ni58_n_to_Ni59_reaclib = 317,
        k_Ni59_to_n_Ni58_reaclib = 318,
        k_Ni59_n_to_Ni60_reaclib = 319,
        k_Ni60_to_n_Ni59_reaclib = 320,
        k_Co57_n_to_Co58_reaclib = 321,
        k_Co58_to_n_Co57_reaclib = 322,
        k_Co58_n_to_Co59_reaclib = 323,
        k_Co59_to_n_Co58_reaclib = 324,
        k_Fe52_n_to_Fe53_reaclib = 325,
        k_Fe53_to_n_Fe52_reaclib = 326,
        k_Fe53_n_to_Fe54_reaclib = 327,
        k_Fe54_to_n_Fe53_reaclib = 328,
        k_Fe54_n_to_Fe55_reaclib = 329,
        k_Fe55_to_n_Fe54_reaclib = 330,
        k_Fe55_n_to_Fe56_reaclib = 331,
        k_Fe56_to_n_Fe55_reaclib = 332,
        k_Fe56_n_to_Fe57_reaclib = 333,
        k_Fe57_to_n_Fe56_reaclib = 334,
        k_Fe57_n_to_Fe58_reaclib = 335,
        k_Fe58_to_n_Fe57_reaclib = 336,
        k_Mn51_n_to_Mn52_reaclib = 337,
        k_Mn52_to_n_Mn51_reaclib = 338,
        k_Mn52_n_to_Mn53_reaclib = 339,
        k_Mn53_to_n_Mn52_reaclib = 340,
        k_Mn53_n_to_Mn54_reaclib = 341,
        k_Mn54_to_n_Mn53_reaclib = 342,
        k_Mn54_n_to_Mn55_reaclib = 343,
        k_Mn55_to_n_Mn54_reaclib = 344,
        k_Cr48_n_to_Cr49_reaclib = 345,
        k_Cr49_to_n_Cr48_reaclib = 346,
        k_Cr49_n_to_Cr50_reaclib = 347,
        k_Cr50_to_n_Cr49_reaclib = 348,
        k_Cr50_n_to_Cr51_reaclib = 349,
        k_Cr51_to_n_Cr50_reaclib = 350,
        k_Cr51_n_to_Cr52_reaclib = 351,
        k_Cr52_to_n_Cr51_reaclib = 352,
        k_Cr52_n_to_Cr53_reaclib = 353,
        k_Cr53_to_n_Cr52_reaclib = 354,
        k_Cr53_n_to_Cr54_reaclib = 355,
        k_Cr54_to_n_Cr53_reaclib = 356,
        k_V47_n_to_V48_reaclib = 357,
        k_V48_to_n_V47_reaclib = 358,
        k_V48_n_to_V49_reaclib = 359,
        k_V49_to_n_V48_reaclib = 360,
        k_V49_n_to_V50_reaclib = 361,
        k_V50_to_n_V49_reaclib = 362,
        k_V50_n_to_V51_reaclib = 363,
        k_V51_to_n_V50_reaclib = 364,
        k_Ti44_n_to_Ti45_reaclib = 365,
        k_Ti45_to_n_Ti44_reaclib = 366,
        k_Ti45_n_to_Ti46_reaclib = 367,
        k_Ti46_to_n_Ti45_reaclib = 368,
        k_Ti46_n_to_Ti47_reaclib = 369,
        k_Ti47_to_n_Ti46_reaclib = 370,
        k_Ti47_n_to_Ti48_reaclib = 371,
        k_Ti48_to_n_Ti47_reaclib = 372,
        k_Sc41_n_to_Sc42_reaclib = 373,
        k_Sc42_to_n_Sc41_reaclib = 374,
        k_Sc42_n_to_Sc43_reaclib = 375,
        k_Sc43_to_n_Sc42_reaclib = 376,
        k_Sc43_n_to_Sc44_reaclib = 377,
        k_Sc44_to_n_Sc43_reaclib = 378,
        k_Sc44_n_to_Sc45_reaclib = 379,
        k_Sc45_to_n_Sc44_reaclib = 380,
        k_Ca40_n_to_Ca41_reaclib = 381,
        k_Ca41_to_n_Ca40_reaclib = 382,
        k_Ca41_n_to_Ca42_reaclib = 383,
        k_Ca42_to_n_Ca41_reaclib = 384,
        k_C12_C12_to_p_Na23_reaclib = 385,
        k_Na23_p_to_C12_C12_reaclib = 386,
        k_Na23_p_to_He4_Ne20_reaclib = 387,
        k_Na23_p_to_Mg24_reaclib = 388,
        k_Mg24_to_p_Na23_reaclib = 389,
        k_C12_C12_to_He4_Ne20_reaclib = 390,
        k_Ne20_He4_to_C12_C12_reaclib = 391,
        k_Ne20_He4_to_p_Na23_reaclib = 392,
        k_Ne20_He4_to_Mg24_reaclib = 393,
        k_Mg24_to_He4_Ne20_reaclib = 394,
        k_Al27_to_Si27_beta_neg_weaktab = 395,
        k_Si27_to_Al27_electron_capture_weaktab = 396,
        k_P29_to_Si29_electron_capture_weaktab = 397,
        k_Si29_to_P29_beta_neg_weaktab = 398,
        k_P31_to_Si31_electron_capture_weaktab = 399,
        k_Si31_to_P31_beta_neg_weaktab = 400,
        k_Cl33_to_S33_electron_capture_weaktab = 401,
        k_P33_to_S33_beta_neg_weaktab = 402,
        k_S33_to_Cl33_beta_neg_weaktab = 403,
        k_S33_to_P33_electron_capture_weaktab = 404,
        k_Cl35_to_S35_electron_capture_weaktab = 405,
        k_S35_to_Cl35_beta_neg_weaktab = 406,
        k_Ar37_to_Cl37_electron_capture_weaktab = 407,
        k_Ar37_to_K37_beta_neg_weaktab = 408,
        k_Cl37_to_Ar37_beta_neg_weaktab = 409,
        k_K37_to_Ar37_electron_capture_weaktab = 410,
        k_Co56_to_Fe56_electron_capture_weaktab = 411,
        k_Co56_to_Ni56_beta_neg_weaktab = 412,
        k_Fe56_to_Co56_beta_neg_weaktab = 413,
        k_Ni56_to_Co56_electron_capture_weaktab = 414,
        k_n_to_p_beta_neg_weaktab = 415,
        k_p_to_n_electron_capture_weaktab = 416,
        k_Ni56_n_n_to_Ni58_approx = 417,
        k_Ni58_to_Ni56_n_n_approx = 418,
        k_Ni58_n_n_to_Ni60_approx = 419,
        k_Ni60_to_Ni58_n_n_approx = 420,
        k_Co57_n_n_to_Co59_approx = 421,
        k_Co59_to_Co57_n_n_approx = 422,
        k_Fe52_n_n_to_Fe54_approx = 423,
        k_Fe54_to_Fe52_n_n_approx = 424,
        k_Fe54_n_n_to_Fe56_approx = 425,
        k_Fe56_to_Fe54_n_n_approx = 426,
        k_Fe56_n_n_to_Fe58_approx = 427,
        k_Fe58_to_Fe56_n_n_approx = 428,
        k_Mn51_n_n_to_Mn53_approx = 429,
        k_Mn53_to_Mn51_n_n_approx = 430,
        k_Mn53_n_n_to_Mn55_approx = 431,
        k_Mn55_to_Mn53_n_n_approx = 432,
        k_Cr48_n_n_to_Cr50_approx = 433,
        k_Cr50_to_Cr48_n_n_approx = 434,
        k_Cr50_n_n_to_Cr52_approx = 435,
        k_Cr52_to_Cr50_n_n_approx = 436,
        k_Cr52_n_n_to_Cr54_approx = 437,
        k_Cr54_to_Cr52_n_n_approx = 438,
        k_V47_n_n_to_V49_approx = 439,
        k_V49_to_V47_n_n_approx = 440,
        k_V49_n_n_to_V51_approx = 441,
        k_V51_to_V49_n_n_approx = 442,
        k_Ti44_n_n_to_Ti46_approx = 443,
        k_Ti46_to_Ti44_n_n_approx = 444,
        k_Ti46_n_n_to_Ti48_approx = 445,
        k_Ti48_to_Ti46_n_n_approx = 446,
        k_Sc41_n_n_to_Sc43_approx = 447,
        k_Sc43_to_Sc41_n_n_approx = 448,
        k_Sc43_n_n_to_Sc45_approx = 449,
        k_Sc45_to_Sc43_n_n_approx = 450,
        k_Ca40_n_n_to_Ca42_approx = 451,
        k_Ca42_to_Ca40_n_n_approx = 452,
        k_C12_C12_to_Mg24_approx = 453,
        k_Mg24_to_C12_C12_approx = 454,
        k_C12_C12_to_Ne20_He4_approx = 455,
        k_Ne20_He4_to_C12_C12_approx = 456,
        k_Ne20_He4_to_Mg24_approx = 457,
        k_Mg24_to_Ne20_He4_approx = 458,
        NumRates = k_Mg24_to_Ne20_He4_approx
    };

    enum ScreenPairs :
    std::uint8_t
    {
        k_p_P33 = 1,
        k_He4_Mn53 = 2,
        k_p_Si31 = 3,
        k_p_Ti44 = 4,
        k_He4_Ti48 = 5,
        k_p_Ca42 = 6,
        k_p_Sc43 = 7,
        k_p_Fe56 = 8,
        k_p_Ar38 = 9,
        k_p_Cr54 = 10,
        k_He4_Ti44 = 11,
        k_He4_Co57 = 12,
        k_He4_Ca42 = 13,
        k_He4_Mg24 = 14,
        k_C12_O16 = 15,
        k_He4_Ar38 = 16,
        k_p_V47 = 17,
        k_He4_V51 = 18,
        k_p_Sc45 = 19,
        k_p_Co55 = 20,
        k_C12_C12 = 21,
        k_He4_Cl33 = 22,
        k_p_Al27 = 23,
        k_He4_Ti46 = 24,
        k_He4_V47 = 25,
        k_He4_Sc45 = 26,
        k_He4_Ni56 = 27,
        k_p_P31 = 28,
        k_p_Ni58 = 29,
        k_p_Si29 = 30,
        k_He4_V49 = 31,
        k_O16_O16 = 32,
        k_He4_P31 = 33,
        k_He4_Si29 = 34,
        k_p_Cr52 = 35,
        k_p_S34 = 36,
        k_p_Ar37 = 37,
        k_p_Cl37 = 38,
        k_He4_Ca40 = 39,
        k_He4_Ar36 = 40,
        k_He4_P33 = 41,
        k_He4_Cr52 = 42,
        k_He4_S34 = 43,
        k_He4_Si31 = 44,
        k_p_Mn55 = 45,
        k_He4_Sc43 = 46,
        k_He4_Fe56 = 47,
        k_He4_Cr54 = 48,
        k_He4_Mn55 = 49,
        k_He4_Ne20 = 50,
        k_p_Co56 = 51,
        k_p_Si28 = 52,
        k_p_Fe54 = 53,
        k_He4_Co55 = 54,
        k_He4_Al27 = 55,
        k_p_Cl35 = 56,
        k_p_S33 = 57,
        k_p_Si30 = 58,
        k_He4_O16 = 59,
        k_He4_Si30 = 60,
        k_p_S35 = 61,
        k_He4_Mg23 = 62,
        k_He4_Ar37 = 63,
        k_p_Co59 = 64,
        k_He4_S35 = 65,
        k_He4_P29 = 66,
        k_p_Cr50 = 67,
        k_p_S32 = 68,
        k_p_Cu59 = 69,
        k_He4_S32 = 70,
        k_p_Mn53 = 71,
        k_He4_Sc41 = 72,
        k_He4_Co56 = 73,
        k_He4_Si28 = 74,
        k_He4_Fe54 = 75,
        k_p_Ti48 = 76,
        k_He4_Be8 = 77,
        k_p_Co57 = 78,
        k_p_K39 = 79,
        k_He4_Cl35 = 80,
        k_He4_S33 = 81,
        k_p_V51 = 82,
        k_He4_C12 = 83,
        k_He4_K39 = 84,
        k_p_Ti46 = 85,
        k_p_Ni60 = 86,
        k_p_Fe58 = 87,
        k_He4_Cr48 = 88,
        k_p_Mn51 = 89,
        k_p_Na23 = 90,
        k_p_V49 = 91,
        k_He4_Fe52 = 92,
        k_He4_K37 = 93,
        k_He4_Cr50 = 94,
        k_He4_Mn51 = 95,
        k_He4_He4 = 96,
        k_C12_Ne20 = 97,
        k_p_Ca40 = 98,
        k_p_Ar36 = 99,
        NumScreenPairs = k_p_Ar36
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "O16_to_He4_C12_reaclib",  // 1,
        "Ne20_to_He4_O16_reaclib",  // 2,
        "Mg24_to_n_Mg23_reaclib",  // 3,
        "Si27_to_He4_Mg23_reaclib",  // 4,
        "Si28_to_n_Si27_reaclib",  // 5,
        "Si28_to_p_Al27_reaclib",  // 6,
        "Si28_to_He4_Mg24_reaclib",  // 7,
        "Si29_to_n_Si28_reaclib",  // 8,
        "Si30_to_n_Si29_reaclib",  // 9,
        "Si31_to_n_Si30_reaclib",  // 10,
        "P29_to_p_Si28_reaclib",  // 11,
        "P31_to_p_Si30_reaclib",  // 12,
        "P31_to_He4_Al27_reaclib",  // 13,
        "S32_to_p_P31_reaclib",  // 14,
        "S32_to_He4_Si28_reaclib",  // 15,
        "S33_to_n_S32_reaclib",  // 16,
        "S33_to_He4_Si29_reaclib",  // 17,
        "S34_to_n_S33_reaclib",  // 18,
        "S34_to_p_P33_reaclib",  // 19,
        "S34_to_He4_Si30_reaclib",  // 20,
        "S35_to_n_S34_reaclib",  // 21,
        "S35_to_He4_Si31_reaclib",  // 22,
        "Cl33_to_p_S32_reaclib",  // 23,
        "Cl33_to_He4_P29_reaclib",  // 24,
        "Cl35_to_p_S34_reaclib",  // 25,
        "Cl35_to_He4_P31_reaclib",  // 26,
        "Cl37_to_He4_P33_reaclib",  // 27,
        "Ar36_to_p_Cl35_reaclib",  // 28,
        "Ar36_to_He4_S32_reaclib",  // 29,
        "Ar37_to_n_Ar36_reaclib",  // 30,
        "Ar37_to_He4_S33_reaclib",  // 31,
        "Ar38_to_n_Ar37_reaclib",  // 32,
        "Ar38_to_p_Cl37_reaclib",  // 33,
        "Ar38_to_He4_S34_reaclib",  // 34,
        "K37_to_p_Ar36_reaclib",  // 35,
        "K37_to_He4_Cl33_reaclib",  // 36,
        "K39_to_p_Ar38_reaclib",  // 37,
        "K39_to_He4_Cl35_reaclib",  // 38,
        "Ca40_to_p_K39_reaclib",  // 39,
        "Ca40_to_He4_Ar36_reaclib",  // 40,
        "Ca42_to_He4_Ar38_reaclib",  // 41,
        "Sc41_to_p_Ca40_reaclib",  // 42,
        "Sc41_to_He4_K37_reaclib",  // 43,
        "Sc43_to_p_Ca42_reaclib",  // 44,
        "Sc43_to_He4_K39_reaclib",  // 45,
        "Ti44_to_p_Sc43_reaclib",  // 46,
        "Ti44_to_He4_Ca40_reaclib",  // 47,
        "Ti46_to_p_Sc45_reaclib",  // 48,
        "Ti46_to_He4_Ca42_reaclib",  // 49,
        "V47_to_p_Ti46_reaclib",  // 50,
        "V47_to_He4_Sc43_reaclib",  // 51,
        "V49_to_p_Ti48_reaclib",  // 52,
        "V49_to_He4_Sc45_reaclib",  // 53,
        "Cr48_to_p_V47_reaclib",  // 54,
        "Cr48_to_He4_Ti44_reaclib",  // 55,
        "Cr50_to_p_V49_reaclib",  // 56,
        "Cr50_to_He4_Ti46_reaclib",  // 57,
        "Cr52_to_p_V51_reaclib",  // 58,
        "Cr52_to_He4_Ti48_reaclib",  // 59,
        "Mn51_to_p_Cr50_reaclib",  // 60,
        "Mn51_to_He4_V47_reaclib",  // 61,
        "Mn53_to_p_Cr52_reaclib",  // 62,
        "Mn53_to_He4_V49_reaclib",  // 63,
        "Mn55_to_p_Cr54_reaclib",  // 64,
        "Mn55_to_He4_V51_reaclib",  // 65,
        "Fe52_to_p_Mn51_reaclib",  // 66,
        "Fe52_to_He4_Cr48_reaclib",  // 67,
        "Fe54_to_p_Mn53_reaclib",  // 68,
        "Fe54_to_He4_Cr50_reaclib",  // 69,
        "Fe56_to_p_Mn55_reaclib",  // 70,
        "Fe56_to_He4_Cr52_reaclib",  // 71,
        "Fe58_to_He4_Cr54_reaclib",  // 72,
        "Co55_to_p_Fe54_reaclib",  // 73,
        "Co55_to_He4_Mn51_reaclib",  // 74,
        "Co56_to_n_Co55_reaclib",  // 75,
        "Co57_to_n_Co56_reaclib",  // 76,
        "Co57_to_p_Fe56_reaclib",  // 77,
        "Co57_to_He4_Mn53_reaclib",  // 78,
        "Co59_to_p_Fe58_reaclib",  // 79,
        "Co59_to_He4_Mn55_reaclib",  // 80,
        "Ni56_to_p_Co55_reaclib",  // 81,
        "Ni56_to_He4_Fe52_reaclib",  // 82,
        "Ni58_to_p_Co57_reaclib",  // 83,
        "Ni58_to_He4_Fe54_reaclib",  // 84,
        "Ni60_to_p_Co59_reaclib",  // 85,
        "Ni60_to_He4_Fe56_reaclib",  // 86,
        "Cu59_to_p_Ni58_reaclib",  // 87,
        "Cu59_to_He4_Co55_reaclib",  // 88,
        "Zn60_to_p_Cu59_reaclib",  // 89,
        "Zn60_to_He4_Ni56_reaclib",  // 90,
        "C12_to_3He4_reaclib",  // 91,
        "C12_He4_to_O16_reaclib",  // 92,
        "O16_He4_to_Ne20_reaclib",  // 93,
        "Mg23_n_to_Mg24_reaclib",  // 94,
        "Mg23_He4_to_Si27_reaclib",  // 95,
        "Mg24_He4_to_Si28_reaclib",  // 96,
        "Al27_p_to_Si28_reaclib",  // 97,
        "Al27_He4_to_P31_reaclib",  // 98,
        "Si27_n_to_Si28_reaclib",  // 99,
        "Si28_n_to_Si29_reaclib",  // 100,
        "Si28_p_to_P29_reaclib",  // 101,
        "Si28_He4_to_S32_reaclib",  // 102,
        "Si29_n_to_Si30_reaclib",  // 103,
        "Si29_He4_to_S33_reaclib",  // 104,
        "Si30_n_to_Si31_reaclib",  // 105,
        "Si30_p_to_P31_reaclib",  // 106,
        "Si30_He4_to_S34_reaclib",  // 107,
        "Si31_He4_to_S35_reaclib",  // 108,
        "P29_He4_to_Cl33_reaclib",  // 109,
        "P31_p_to_S32_reaclib",  // 110,
        "P31_He4_to_Cl35_reaclib",  // 111,
        "P33_p_to_S34_reaclib",  // 112,
        "P33_He4_to_Cl37_reaclib",  // 113,
        "S32_n_to_S33_reaclib",  // 114,
        "S32_p_to_Cl33_reaclib",  // 115,
        "S32_He4_to_Ar36_reaclib",  // 116,
        "S33_n_to_S34_reaclib",  // 117,
        "S33_He4_to_Ar37_reaclib",  // 118,
        "S34_n_to_S35_reaclib",  // 119,
        "S34_p_to_Cl35_reaclib",  // 120,
        "S34_He4_to_Ar38_reaclib",  // 121,
        "Cl33_He4_to_K37_reaclib",  // 122,
        "Cl35_p_to_Ar36_reaclib",  // 123,
        "Cl35_He4_to_K39_reaclib",  // 124,
        "Cl37_p_to_Ar38_reaclib",  // 125,
        "Ar36_n_to_Ar37_reaclib",  // 126,
        "Ar36_p_to_K37_reaclib",  // 127,
        "Ar36_He4_to_Ca40_reaclib",  // 128,
        "Ar37_n_to_Ar38_reaclib",  // 129,
        "Ar38_p_to_K39_reaclib",  // 130,
        "Ar38_He4_to_Ca42_reaclib",  // 131,
        "K37_He4_to_Sc41_reaclib",  // 132,
        "K39_p_to_Ca40_reaclib",  // 133,
        "K39_He4_to_Sc43_reaclib",  // 134,
        "Ca40_p_to_Sc41_reaclib",  // 135,
        "Ca40_He4_to_Ti44_reaclib",  // 136,
        "Ca42_p_to_Sc43_reaclib",  // 137,
        "Ca42_He4_to_Ti46_reaclib",  // 138,
        "Sc43_p_to_Ti44_reaclib",  // 139,
        "Sc43_He4_to_V47_reaclib",  // 140,
        "Sc45_p_to_Ti46_reaclib",  // 141,
        "Sc45_He4_to_V49_reaclib",  // 142,
        "Ti44_He4_to_Cr48_reaclib",  // 143,
        "Ti46_p_to_V47_reaclib",  // 144,
        "Ti46_He4_to_Cr50_reaclib",  // 145,
        "Ti48_p_to_V49_reaclib",  // 146,
        "Ti48_He4_to_Cr52_reaclib",  // 147,
        "V47_p_to_Cr48_reaclib",  // 148,
        "V47_He4_to_Mn51_reaclib",  // 149,
        "V49_p_to_Cr50_reaclib",  // 150,
        "V49_He4_to_Mn53_reaclib",  // 151,
        "V51_p_to_Cr52_reaclib",  // 152,
        "V51_He4_to_Mn55_reaclib",  // 153,
        "Cr48_He4_to_Fe52_reaclib",  // 154,
        "Cr50_p_to_Mn51_reaclib",  // 155,
        "Cr50_He4_to_Fe54_reaclib",  // 156,
        "Cr52_p_to_Mn53_reaclib",  // 157,
        "Cr52_He4_to_Fe56_reaclib",  // 158,
        "Cr54_p_to_Mn55_reaclib",  // 159,
        "Cr54_He4_to_Fe58_reaclib",  // 160,
        "Mn51_p_to_Fe52_reaclib",  // 161,
        "Mn51_He4_to_Co55_reaclib",  // 162,
        "Mn53_p_to_Fe54_reaclib",  // 163,
        "Mn53_He4_to_Co57_reaclib",  // 164,
        "Mn55_p_to_Fe56_reaclib",  // 165,
        "Mn55_He4_to_Co59_reaclib",  // 166,
        "Fe52_He4_to_Ni56_reaclib",  // 167,
        "Fe54_p_to_Co55_reaclib",  // 168,
        "Fe54_He4_to_Ni58_reaclib",  // 169,
        "Fe56_p_to_Co57_reaclib",  // 170,
        "Fe56_He4_to_Ni60_reaclib",  // 171,
        "Fe58_p_to_Co59_reaclib",  // 172,
        "Co55_n_to_Co56_reaclib",  // 173,
        "Co55_p_to_Ni56_reaclib",  // 174,
        "Co55_He4_to_Cu59_reaclib",  // 175,
        "Co56_n_to_Co57_reaclib",  // 176,
        "Co57_p_to_Ni58_reaclib",  // 177,
        "Co59_p_to_Ni60_reaclib",  // 178,
        "Ni56_He4_to_Zn60_reaclib",  // 179,
        "Ni58_p_to_Cu59_reaclib",  // 180,
        "Cu59_p_to_Zn60_reaclib",  // 181,
        "C12_C12_to_n_Mg23_reaclib",  // 182,
        "O16_C12_to_n_Si27_reaclib",  // 183,
        "O16_C12_to_p_Al27_reaclib",  // 184,
        "O16_C12_to_He4_Mg24_reaclib",  // 185,
        "O16_O16_to_p_P31_reaclib",  // 186,
        "O16_O16_to_He4_Si28_reaclib",  // 187,
        "Ne20_He4_to_n_Mg23_reaclib",  // 188,
        "Ne20_C12_to_p_P31_reaclib",  // 189,
        "Ne20_C12_to_He4_Si28_reaclib",  // 190,
        "Mg23_n_to_He4_Ne20_reaclib",  // 191,
        "Mg23_n_to_C12_C12_reaclib",  // 192,
        "Mg24_He4_to_n_Si27_reaclib",  // 193,
        "Mg24_He4_to_p_Al27_reaclib",  // 194,
        "Mg24_He4_to_C12_O16_reaclib",  // 195,
        "Al27_p_to_n_Si27_reaclib",  // 196,
        "Al27_p_to_He4_Mg24_reaclib",  // 197,
        "Al27_p_to_C12_O16_reaclib",  // 198,
        "Al27_He4_to_p_Si30_reaclib",  // 199,
        "Si27_n_to_p_Al27_reaclib",  // 200,
        "Si27_n_to_He4_Mg24_reaclib",  // 201,
        "Si27_n_to_C12_O16_reaclib",  // 202,
        "Si28_He4_to_p_P31_reaclib",  // 203,
        "Si28_He4_to_C12_Ne20_reaclib",  // 204,
        "Si28_He4_to_O16_O16_reaclib",  // 205,
        "Si29_p_to_n_P29_reaclib",  // 206,
        "Si29_He4_to_n_S32_reaclib",  // 207,
        "Si30_p_to_He4_Al27_reaclib",  // 208,
        "Si30_He4_to_n_S33_reaclib",  // 209,
        "Si30_He4_to_p_P33_reaclib",  // 210,
        "Si31_p_to_n_P31_reaclib",  // 211,
        "Si31_He4_to_n_S34_reaclib",  // 212,
        "P29_n_to_p_Si29_reaclib",  // 213,
        "P29_He4_to_p_S32_reaclib",  // 214,
        "P31_n_to_p_Si31_reaclib",  // 215,
        "P31_p_to_He4_Si28_reaclib",  // 216,
        "P31_p_to_C12_Ne20_reaclib",  // 217,
        "P31_p_to_O16_O16_reaclib",  // 218,
        "P31_He4_to_p_S34_reaclib",  // 219,
        "P33_p_to_n_S33_reaclib",  // 220,
        "P33_p_to_He4_Si30_reaclib",  // 221,
        "S32_n_to_He4_Si29_reaclib",  // 222,
        "S32_p_to_He4_P29_reaclib",  // 223,
        "S32_He4_to_p_Cl35_reaclib",  // 224,
        "S33_n_to_p_P33_reaclib",  // 225,
        "S33_n_to_He4_Si30_reaclib",  // 226,
        "S33_p_to_n_Cl33_reaclib",  // 227,
        "S33_He4_to_n_Ar36_reaclib",  // 228,
        "S34_n_to_He4_Si31_reaclib",  // 229,
        "S34_p_to_He4_P31_reaclib",  // 230,
        "S34_He4_to_n_Ar37_reaclib",  // 231,
        "S34_He4_to_p_Cl37_reaclib",  // 232,
        "S35_p_to_n_Cl35_reaclib",  // 233,
        "S35_He4_to_n_Ar38_reaclib",  // 234,
        "Cl33_n_to_p_S33_reaclib",  // 235,
        "Cl33_He4_to_p_Ar36_reaclib",  // 236,
        "Cl35_n_to_p_S35_reaclib",  // 237,
        "Cl35_p_to_He4_S32_reaclib",  // 238,
        "Cl35_He4_to_p_Ar38_reaclib",  // 239,
        "Cl37_p_to_n_Ar37_reaclib",  // 240,
        "Cl37_p_to_He4_S34_reaclib",  // 241,
        "Ar36_n_to_He4_S33_reaclib",  // 242,
        "Ar36_p_to_He4_Cl33_reaclib",  // 243,
        "Ar36_He4_to_p_K39_reaclib",  // 244,
        "Ar37_n_to_p_Cl37_reaclib",  // 245,
        "Ar37_n_to_He4_S34_reaclib",  // 246,
        "Ar37_p_to_n_K37_reaclib",  // 247,
        "Ar37_He4_to_n_Ca40_reaclib",  // 248,
        "Ar38_n_to_He4_S35_reaclib",  // 249,
        "Ar38_p_to_He4_Cl35_reaclib",  // 250,
        "K37_n_to_p_Ar37_reaclib",  // 251,
        "K37_He4_to_p_Ca40_reaclib",  // 252,
        "K39_p_to_He4_Ar36_reaclib",  // 253,
        "K39_He4_to_p_Ca42_reaclib",  // 254,
        "Ca40_n_to_He4_Ar37_reaclib",  // 255,
        "Ca40_p_to_He4_K37_reaclib",  // 256,
        "Ca40_He4_to_p_Sc43_reaclib",  // 257,
        "Ca42_p_to_He4_K39_reaclib",  // 258,
        "Ca42_He4_to_p_Sc45_reaclib",  // 259,
        "Sc41_He4_to_p_Ti44_reaclib",  // 260,
        "Sc43_p_to_He4_Ca40_reaclib",  // 261,
        "Sc43_He4_to_p_Ti46_reaclib",  // 262,
        "Sc45_p_to_He4_Ca42_reaclib",  // 263,
        "Sc45_He4_to_p_Ti48_reaclib",  // 264,
        "Ti44_p_to_He4_Sc41_reaclib",  // 265,
        "Ti44_He4_to_p_V47_reaclib",  // 266,
        "Ti46_p_to_He4_Sc43_reaclib",  // 267,
        "Ti46_He4_to_p_V49_reaclib",  // 268,
        "Ti48_p_to_He4_Sc45_reaclib",  // 269,
        "Ti48_He4_to_p_V51_reaclib",  // 270,
        "V47_p_to_He4_Ti44_reaclib",  // 271,
        "V47_He4_to_p_Cr50_reaclib",  // 272,
        "V49_p_to_He4_Ti46_reaclib",  // 273,
        "V49_He4_to_p_Cr52_reaclib",  // 274,
        "V51_p_to_He4_Ti48_reaclib",  // 275,
        "V51_He4_to_p_Cr54_reaclib",  // 276,
        "Cr48_He4_to_p_Mn51_reaclib",  // 277,
        "Cr50_p_to_He4_V47_reaclib",  // 278,
        "Cr50_He4_to_p_Mn53_reaclib",  // 279,
        "Cr52_p_to_He4_V49_reaclib",  // 280,
        "Cr52_He4_to_p_Mn55_reaclib",  // 281,
        "Cr54_p_to_He4_V51_reaclib",  // 282,
        "Mn51_p_to_He4_Cr48_reaclib",  // 283,
        "Mn51_He4_to_p_Fe54_reaclib",  // 284,
        "Mn53_p_to_He4_Cr50_reaclib",  // 285,
        "Mn53_He4_to_n_Co56_reaclib",  // 286,
        "Mn53_He4_to_p_Fe56_reaclib",  // 287,
        "Mn55_p_to_He4_Cr52_reaclib",  // 288,
        "Mn55_He4_to_p_Fe58_reaclib",  // 289,
        "Fe52_He4_to_p_Co55_reaclib",  // 290,
        "Fe54_p_to_He4_Mn51_reaclib",  // 291,
        "Fe54_He4_to_p_Co57_reaclib",  // 292,
        "Fe56_p_to_n_Co56_reaclib",  // 293,
        "Fe56_p_to_He4_Mn53_reaclib",  // 294,
        "Fe56_He4_to_p_Co59_reaclib",  // 295,
        "Fe58_p_to_He4_Mn55_reaclib",  // 296,
        "Co55_p_to_He4_Fe52_reaclib",  // 297,
        "Co55_He4_to_p_Ni58_reaclib",  // 298,
        "Co56_n_to_p_Fe56_reaclib",  // 299,
        "Co56_n_to_He4_Mn53_reaclib",  // 300,
        "Co56_p_to_n_Ni56_reaclib",  // 301,
        "Co56_He4_to_n_Cu59_reaclib",  // 302,
        "Co57_p_to_He4_Fe54_reaclib",  // 303,
        "Co57_He4_to_p_Ni60_reaclib",  // 304,
        "Co59_p_to_He4_Fe56_reaclib",  // 305,
        "Ni56_n_to_p_Co56_reaclib",  // 306,
        "Ni56_He4_to_p_Cu59_reaclib",  // 307,
        "Ni58_p_to_He4_Co55_reaclib",  // 308,
        "Ni60_p_to_He4_Co57_reaclib",  // 309,
        "Cu59_n_to_He4_Co56_reaclib",  // 310,
        "Cu59_p_to_He4_Ni56_reaclib",  // 311,
        "He4_He4_He4_to_C12_reaclib",  // 312,
        "Ni56_n_to_Ni57_reaclib",  // 313,
        "Ni57_to_n_Ni56_reaclib",  // 314,
        "Ni57_n_to_Ni58_reaclib",  // 315,
        "Ni58_to_n_Ni57_reaclib",  // 316,
        "Ni58_n_to_Ni59_reaclib",  // 317,
        "Ni59_to_n_Ni58_reaclib",  // 318,
        "Ni59_n_to_Ni60_reaclib",  // 319,
        "Ni60_to_n_Ni59_reaclib",  // 320,
        "Co57_n_to_Co58_reaclib",  // 321,
        "Co58_to_n_Co57_reaclib",  // 322,
        "Co58_n_to_Co59_reaclib",  // 323,
        "Co59_to_n_Co58_reaclib",  // 324,
        "Fe52_n_to_Fe53_reaclib",  // 325,
        "Fe53_to_n_Fe52_reaclib",  // 326,
        "Fe53_n_to_Fe54_reaclib",  // 327,
        "Fe54_to_n_Fe53_reaclib",  // 328,
        "Fe54_n_to_Fe55_reaclib",  // 329,
        "Fe55_to_n_Fe54_reaclib",  // 330,
        "Fe55_n_to_Fe56_reaclib",  // 331,
        "Fe56_to_n_Fe55_reaclib",  // 332,
        "Fe56_n_to_Fe57_reaclib",  // 333,
        "Fe57_to_n_Fe56_reaclib",  // 334,
        "Fe57_n_to_Fe58_reaclib",  // 335,
        "Fe58_to_n_Fe57_reaclib",  // 336,
        "Mn51_n_to_Mn52_reaclib",  // 337,
        "Mn52_to_n_Mn51_reaclib",  // 338,
        "Mn52_n_to_Mn53_reaclib",  // 339,
        "Mn53_to_n_Mn52_reaclib",  // 340,
        "Mn53_n_to_Mn54_reaclib",  // 341,
        "Mn54_to_n_Mn53_reaclib",  // 342,
        "Mn54_n_to_Mn55_reaclib",  // 343,
        "Mn55_to_n_Mn54_reaclib",  // 344,
        "Cr48_n_to_Cr49_reaclib",  // 345,
        "Cr49_to_n_Cr48_reaclib",  // 346,
        "Cr49_n_to_Cr50_reaclib",  // 347,
        "Cr50_to_n_Cr49_reaclib",  // 348,
        "Cr50_n_to_Cr51_reaclib",  // 349,
        "Cr51_to_n_Cr50_reaclib",  // 350,
        "Cr51_n_to_Cr52_reaclib",  // 351,
        "Cr52_to_n_Cr51_reaclib",  // 352,
        "Cr52_n_to_Cr53_reaclib",  // 353,
        "Cr53_to_n_Cr52_reaclib",  // 354,
        "Cr53_n_to_Cr54_reaclib",  // 355,
        "Cr54_to_n_Cr53_reaclib",  // 356,
        "V47_n_to_V48_reaclib",  // 357,
        "V48_to_n_V47_reaclib",  // 358,
        "V48_n_to_V49_reaclib",  // 359,
        "V49_to_n_V48_reaclib",  // 360,
        "V49_n_to_V50_reaclib",  // 361,
        "V50_to_n_V49_reaclib",  // 362,
        "V50_n_to_V51_reaclib",  // 363,
        "V51_to_n_V50_reaclib",  // 364,
        "Ti44_n_to_Ti45_reaclib",  // 365,
        "Ti45_to_n_Ti44_reaclib",  // 366,
        "Ti45_n_to_Ti46_reaclib",  // 367,
        "Ti46_to_n_Ti45_reaclib",  // 368,
        "Ti46_n_to_Ti47_reaclib",  // 369,
        "Ti47_to_n_Ti46_reaclib",  // 370,
        "Ti47_n_to_Ti48_reaclib",  // 371,
        "Ti48_to_n_Ti47_reaclib",  // 372,
        "Sc41_n_to_Sc42_reaclib",  // 373,
        "Sc42_to_n_Sc41_reaclib",  // 374,
        "Sc42_n_to_Sc43_reaclib",  // 375,
        "Sc43_to_n_Sc42_reaclib",  // 376,
        "Sc43_n_to_Sc44_reaclib",  // 377,
        "Sc44_to_n_Sc43_reaclib",  // 378,
        "Sc44_n_to_Sc45_reaclib",  // 379,
        "Sc45_to_n_Sc44_reaclib",  // 380,
        "Ca40_n_to_Ca41_reaclib",  // 381,
        "Ca41_to_n_Ca40_reaclib",  // 382,
        "Ca41_n_to_Ca42_reaclib",  // 383,
        "Ca42_to_n_Ca41_reaclib",  // 384,
        "C12_C12_to_p_Na23_reaclib",  // 385,
        "Na23_p_to_C12_C12_reaclib",  // 386,
        "Na23_p_to_He4_Ne20_reaclib",  // 387,
        "Na23_p_to_Mg24_reaclib",  // 388,
        "Mg24_to_p_Na23_reaclib",  // 389,
        "C12_C12_to_He4_Ne20_reaclib",  // 390,
        "Ne20_He4_to_C12_C12_reaclib",  // 391,
        "Ne20_He4_to_p_Na23_reaclib",  // 392,
        "Ne20_He4_to_Mg24_reaclib",  // 393,
        "Mg24_to_He4_Ne20_reaclib",  // 394,
        "Al27_to_Si27_beta_neg_weaktab",  // 395,
        "Si27_to_Al27_electron_capture_weaktab",  // 396,
        "P29_to_Si29_electron_capture_weaktab",  // 397,
        "Si29_to_P29_beta_neg_weaktab",  // 398,
        "P31_to_Si31_electron_capture_weaktab",  // 399,
        "Si31_to_P31_beta_neg_weaktab",  // 400,
        "Cl33_to_S33_electron_capture_weaktab",  // 401,
        "P33_to_S33_beta_neg_weaktab",  // 402,
        "S33_to_Cl33_beta_neg_weaktab",  // 403,
        "S33_to_P33_electron_capture_weaktab",  // 404,
        "Cl35_to_S35_electron_capture_weaktab",  // 405,
        "S35_to_Cl35_beta_neg_weaktab",  // 406,
        "Ar37_to_Cl37_electron_capture_weaktab",  // 407,
        "Ar37_to_K37_beta_neg_weaktab",  // 408,
        "Cl37_to_Ar37_beta_neg_weaktab",  // 409,
        "K37_to_Ar37_electron_capture_weaktab",  // 410,
        "Co56_to_Fe56_electron_capture_weaktab",  // 411,
        "Co56_to_Ni56_beta_neg_weaktab",  // 412,
        "Fe56_to_Co56_beta_neg_weaktab",  // 413,
        "Ni56_to_Co56_electron_capture_weaktab",  // 414,
        "n_to_p_beta_neg_weaktab",  // 415,
        "p_to_n_electron_capture_weaktab",  // 416,
        "Ni56_n_n_to_Ni58_approx",  // 417,
        "Ni58_to_Ni56_n_n_approx",  // 418,
        "Ni58_n_n_to_Ni60_approx",  // 419,
        "Ni60_to_Ni58_n_n_approx",  // 420,
        "Co57_n_n_to_Co59_approx",  // 421,
        "Co59_to_Co57_n_n_approx",  // 422,
        "Fe52_n_n_to_Fe54_approx",  // 423,
        "Fe54_to_Fe52_n_n_approx",  // 424,
        "Fe54_n_n_to_Fe56_approx",  // 425,
        "Fe56_to_Fe54_n_n_approx",  // 426,
        "Fe56_n_n_to_Fe58_approx",  // 427,
        "Fe58_to_Fe56_n_n_approx",  // 428,
        "Mn51_n_n_to_Mn53_approx",  // 429,
        "Mn53_to_Mn51_n_n_approx",  // 430,
        "Mn53_n_n_to_Mn55_approx",  // 431,
        "Mn55_to_Mn53_n_n_approx",  // 432,
        "Cr48_n_n_to_Cr50_approx",  // 433,
        "Cr50_to_Cr48_n_n_approx",  // 434,
        "Cr50_n_n_to_Cr52_approx",  // 435,
        "Cr52_to_Cr50_n_n_approx",  // 436,
        "Cr52_n_n_to_Cr54_approx",  // 437,
        "Cr54_to_Cr52_n_n_approx",  // 438,
        "V47_n_n_to_V49_approx",  // 439,
        "V49_to_V47_n_n_approx",  // 440,
        "V49_n_n_to_V51_approx",  // 441,
        "V51_to_V49_n_n_approx",  // 442,
        "Ti44_n_n_to_Ti46_approx",  // 443,
        "Ti46_to_Ti44_n_n_approx",  // 444,
        "Ti46_n_n_to_Ti48_approx",  // 445,
        "Ti48_to_Ti46_n_n_approx",  // 446,
        "Sc41_n_n_to_Sc43_approx",  // 447,
        "Sc43_to_Sc41_n_n_approx",  // 448,
        "Sc43_n_n_to_Sc45_approx",  // 449,
        "Sc45_to_Sc43_n_n_approx",  // 450,
        "Ca40_n_n_to_Ca42_approx",  // 451,
        "Ca42_to_Ca40_n_n_approx",  // 452,
        "C12_C12_to_Mg24_approx",  // 453,
        "Mg24_to_C12_C12_approx",  // 454,
        "C12_C12_to_Ne20_He4_approx",  // 455,
        "Ne20_He4_to_C12_C12_approx",  // 456,
        "Ne20_He4_to_Mg24_approx",  // 457,
        "Mg24_to_Ne20_He4_approx"  // 458,
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
    constexpr int NumNSERatePairs = 149;

    inline AMREX_GPU_MANAGED amrex::Array2D<short, 1, NumNSERatePairs, 1, 10, amrex::Order::C> rate_pair_data {
        He4, C12, -1, O16, -1, -1, C12, O16, k_C12_He4_to_O16_reaclib, k_O16_to_He4_C12_reaclib,
        He4, O16, -1, Ne20, -1, -1, O16, Ne20, k_O16_He4_to_Ne20_reaclib, k_Ne20_to_He4_O16_reaclib,
        N, Mg23, -1, Mg24, -1, -1, Mg23, Mg24, k_Mg23_n_to_Mg24_reaclib, k_Mg24_to_n_Mg23_reaclib,
        He4, Mg23, -1, Si27, -1, -1, Mg23, Si27, k_Mg23_He4_to_Si27_reaclib, k_Si27_to_He4_Mg23_reaclib,
        He4, Mg24, -1, Si28, -1, -1, Mg24, Si28, k_Mg24_He4_to_Si28_reaclib, k_Si28_to_He4_Mg24_reaclib,
        H1, Al27, -1, Si28, -1, -1, Al27, Si28, k_Al27_p_to_Si28_reaclib, k_Si28_to_p_Al27_reaclib,
        He4, Al27, -1, P31, -1, -1, Al27, P31, k_Al27_He4_to_P31_reaclib, k_P31_to_He4_Al27_reaclib,
        N, Si27, -1, Si28, -1, -1, Si27, Si28, k_Si27_n_to_Si28_reaclib, k_Si28_to_n_Si27_reaclib,
        N, Si28, -1, Si29, -1, -1, Si28, Si29, k_Si28_n_to_Si29_reaclib, k_Si29_to_n_Si28_reaclib,
        H1, Si28, -1, P29, -1, -1, Si28, P29, k_Si28_p_to_P29_reaclib, k_P29_to_p_Si28_reaclib,
        He4, Si28, -1, S32, -1, -1, Si28, S32, k_Si28_He4_to_S32_reaclib, k_S32_to_He4_Si28_reaclib,
        N, Si29, -1, Si30, -1, -1, Si29, Si30, k_Si29_n_to_Si30_reaclib, k_Si30_to_n_Si29_reaclib,
        He4, Si29, -1, S33, -1, -1, Si29, S33, k_Si29_He4_to_S33_reaclib, k_S33_to_He4_Si29_reaclib,
        N, Si30, -1, Si31, -1, -1, Si30, Si31, k_Si30_n_to_Si31_reaclib, k_Si31_to_n_Si30_reaclib,
        H1, Si30, -1, P31, -1, -1, Si30, P31, k_Si30_p_to_P31_reaclib, k_P31_to_p_Si30_reaclib,
        He4, Si30, -1, S34, -1, -1, Si30, S34, k_Si30_He4_to_S34_reaclib, k_S34_to_He4_Si30_reaclib,
        He4, Si31, -1, S35, -1, -1, Si31, S35, k_Si31_He4_to_S35_reaclib, k_S35_to_He4_Si31_reaclib,
        He4, P29, -1, Cl33, -1, -1, P29, Cl33, k_P29_He4_to_Cl33_reaclib, k_Cl33_to_He4_P29_reaclib,
        H1, P31, -1, S32, -1, -1, P31, S32, k_P31_p_to_S32_reaclib, k_S32_to_p_P31_reaclib,
        He4, P31, -1, Cl35, -1, -1, P31, Cl35, k_P31_He4_to_Cl35_reaclib, k_Cl35_to_He4_P31_reaclib,
        H1, P33, -1, S34, -1, -1, P33, S34, k_P33_p_to_S34_reaclib, k_S34_to_p_P33_reaclib,
        He4, P33, -1, Cl37, -1, -1, P33, Cl37, k_P33_He4_to_Cl37_reaclib, k_Cl37_to_He4_P33_reaclib,
        N, S32, -1, S33, -1, -1, S32, S33, k_S32_n_to_S33_reaclib, k_S33_to_n_S32_reaclib,
        H1, S32, -1, Cl33, -1, -1, S32, Cl33, k_S32_p_to_Cl33_reaclib, k_Cl33_to_p_S32_reaclib,
        He4, S32, -1, Ar36, -1, -1, S32, Ar36, k_S32_He4_to_Ar36_reaclib, k_Ar36_to_He4_S32_reaclib,
        N, S33, -1, S34, -1, -1, S33, S34, k_S33_n_to_S34_reaclib, k_S34_to_n_S33_reaclib,
        He4, S33, -1, Ar37, -1, -1, S33, Ar37, k_S33_He4_to_Ar37_reaclib, k_Ar37_to_He4_S33_reaclib,
        N, S34, -1, S35, -1, -1, S34, S35, k_S34_n_to_S35_reaclib, k_S35_to_n_S34_reaclib,
        H1, S34, -1, Cl35, -1, -1, S34, Cl35, k_S34_p_to_Cl35_reaclib, k_Cl35_to_p_S34_reaclib,
        He4, S34, -1, Ar38, -1, -1, S34, Ar38, k_S34_He4_to_Ar38_reaclib, k_Ar38_to_He4_S34_reaclib,
        He4, Cl33, -1, K37, -1, -1, Cl33, K37, k_Cl33_He4_to_K37_reaclib, k_K37_to_He4_Cl33_reaclib,
        H1, Cl35, -1, Ar36, -1, -1, Cl35, Ar36, k_Cl35_p_to_Ar36_reaclib, k_Ar36_to_p_Cl35_reaclib,
        He4, Cl35, -1, K39, -1, -1, Cl35, K39, k_Cl35_He4_to_K39_reaclib, k_K39_to_He4_Cl35_reaclib,
        H1, Cl37, -1, Ar38, -1, -1, Cl37, Ar38, k_Cl37_p_to_Ar38_reaclib, k_Ar38_to_p_Cl37_reaclib,
        N, Ar36, -1, Ar37, -1, -1, Ar36, Ar37, k_Ar36_n_to_Ar37_reaclib, k_Ar37_to_n_Ar36_reaclib,
        H1, Ar36, -1, K37, -1, -1, Ar36, K37, k_Ar36_p_to_K37_reaclib, k_K37_to_p_Ar36_reaclib,
        He4, Ar36, -1, Ca40, -1, -1, Ar36, Ca40, k_Ar36_He4_to_Ca40_reaclib, k_Ca40_to_He4_Ar36_reaclib,
        N, Ar37, -1, Ar38, -1, -1, Ar37, Ar38, k_Ar37_n_to_Ar38_reaclib, k_Ar38_to_n_Ar37_reaclib,
        H1, Ar38, -1, K39, -1, -1, Ar38, K39, k_Ar38_p_to_K39_reaclib, k_K39_to_p_Ar38_reaclib,
        He4, Ar38, -1, Ca42, -1, -1, Ar38, Ca42, k_Ar38_He4_to_Ca42_reaclib, k_Ca42_to_He4_Ar38_reaclib,
        He4, K37, -1, Sc41, -1, -1, K37, Sc41, k_K37_He4_to_Sc41_reaclib, k_Sc41_to_He4_K37_reaclib,
        H1, K39, -1, Ca40, -1, -1, K39, Ca40, k_K39_p_to_Ca40_reaclib, k_Ca40_to_p_K39_reaclib,
        He4, K39, -1, Sc43, -1, -1, K39, Sc43, k_K39_He4_to_Sc43_reaclib, k_Sc43_to_He4_K39_reaclib,
        H1, Ca40, -1, Sc41, -1, -1, Ca40, Sc41, k_Ca40_p_to_Sc41_reaclib, k_Sc41_to_p_Ca40_reaclib,
        He4, Ca40, -1, Ti44, -1, -1, Ca40, Ti44, k_Ca40_He4_to_Ti44_reaclib, k_Ti44_to_He4_Ca40_reaclib,
        H1, Ca42, -1, Sc43, -1, -1, Ca42, Sc43, k_Ca42_p_to_Sc43_reaclib, k_Sc43_to_p_Ca42_reaclib,
        He4, Ca42, -1, Ti46, -1, -1, Ca42, Ti46, k_Ca42_He4_to_Ti46_reaclib, k_Ti46_to_He4_Ca42_reaclib,
        H1, Sc43, -1, Ti44, -1, -1, Sc43, Ti44, k_Sc43_p_to_Ti44_reaclib, k_Ti44_to_p_Sc43_reaclib,
        He4, Sc43, -1, V47, -1, -1, Sc43, V47, k_Sc43_He4_to_V47_reaclib, k_V47_to_He4_Sc43_reaclib,
        H1, Sc45, -1, Ti46, -1, -1, Sc45, Ti46, k_Sc45_p_to_Ti46_reaclib, k_Ti46_to_p_Sc45_reaclib,
        He4, Sc45, -1, V49, -1, -1, Sc45, V49, k_Sc45_He4_to_V49_reaclib, k_V49_to_He4_Sc45_reaclib,
        He4, Ti44, -1, Cr48, -1, -1, Ti44, Cr48, k_Ti44_He4_to_Cr48_reaclib, k_Cr48_to_He4_Ti44_reaclib,
        H1, Ti46, -1, V47, -1, -1, Ti46, V47, k_Ti46_p_to_V47_reaclib, k_V47_to_p_Ti46_reaclib,
        He4, Ti46, -1, Cr50, -1, -1, Ti46, Cr50, k_Ti46_He4_to_Cr50_reaclib, k_Cr50_to_He4_Ti46_reaclib,
        H1, Ti48, -1, V49, -1, -1, Ti48, V49, k_Ti48_p_to_V49_reaclib, k_V49_to_p_Ti48_reaclib,
        He4, Ti48, -1, Cr52, -1, -1, Ti48, Cr52, k_Ti48_He4_to_Cr52_reaclib, k_Cr52_to_He4_Ti48_reaclib,
        H1, V47, -1, Cr48, -1, -1, V47, Cr48, k_V47_p_to_Cr48_reaclib, k_Cr48_to_p_V47_reaclib,
        He4, V47, -1, Mn51, -1, -1, V47, Mn51, k_V47_He4_to_Mn51_reaclib, k_Mn51_to_He4_V47_reaclib,
        H1, V49, -1, Cr50, -1, -1, V49, Cr50, k_V49_p_to_Cr50_reaclib, k_Cr50_to_p_V49_reaclib,
        He4, V49, -1, Mn53, -1, -1, V49, Mn53, k_V49_He4_to_Mn53_reaclib, k_Mn53_to_He4_V49_reaclib,
        H1, V51, -1, Cr52, -1, -1, V51, Cr52, k_V51_p_to_Cr52_reaclib, k_Cr52_to_p_V51_reaclib,
        He4, V51, -1, Mn55, -1, -1, V51, Mn55, k_V51_He4_to_Mn55_reaclib, k_Mn55_to_He4_V51_reaclib,
        He4, Cr48, -1, Fe52, -1, -1, Cr48, Fe52, k_Cr48_He4_to_Fe52_reaclib, k_Fe52_to_He4_Cr48_reaclib,
        H1, Cr50, -1, Mn51, -1, -1, Cr50, Mn51, k_Cr50_p_to_Mn51_reaclib, k_Mn51_to_p_Cr50_reaclib,
        He4, Cr50, -1, Fe54, -1, -1, Cr50, Fe54, k_Cr50_He4_to_Fe54_reaclib, k_Fe54_to_He4_Cr50_reaclib,
        H1, Cr52, -1, Mn53, -1, -1, Cr52, Mn53, k_Cr52_p_to_Mn53_reaclib, k_Mn53_to_p_Cr52_reaclib,
        He4, Cr52, -1, Fe56, -1, -1, Cr52, Fe56, k_Cr52_He4_to_Fe56_reaclib, k_Fe56_to_He4_Cr52_reaclib,
        H1, Cr54, -1, Mn55, -1, -1, Cr54, Mn55, k_Cr54_p_to_Mn55_reaclib, k_Mn55_to_p_Cr54_reaclib,
        He4, Cr54, -1, Fe58, -1, -1, Cr54, Fe58, k_Cr54_He4_to_Fe58_reaclib, k_Fe58_to_He4_Cr54_reaclib,
        H1, Mn51, -1, Fe52, -1, -1, Mn51, Fe52, k_Mn51_p_to_Fe52_reaclib, k_Fe52_to_p_Mn51_reaclib,
        He4, Mn51, -1, Co55, -1, -1, Mn51, Co55, k_Mn51_He4_to_Co55_reaclib, k_Co55_to_He4_Mn51_reaclib,
        H1, Mn53, -1, Fe54, -1, -1, Mn53, Fe54, k_Mn53_p_to_Fe54_reaclib, k_Fe54_to_p_Mn53_reaclib,
        He4, Mn53, -1, Co57, -1, -1, Mn53, Co57, k_Mn53_He4_to_Co57_reaclib, k_Co57_to_He4_Mn53_reaclib,
        H1, Mn55, -1, Fe56, -1, -1, Mn55, Fe56, k_Mn55_p_to_Fe56_reaclib, k_Fe56_to_p_Mn55_reaclib,
        He4, Mn55, -1, Co59, -1, -1, Mn55, Co59, k_Mn55_He4_to_Co59_reaclib, k_Co59_to_He4_Mn55_reaclib,
        He4, Fe52, -1, Ni56, -1, -1, Fe52, Ni56, k_Fe52_He4_to_Ni56_reaclib, k_Ni56_to_He4_Fe52_reaclib,
        H1, Fe54, -1, Co55, -1, -1, Fe54, Co55, k_Fe54_p_to_Co55_reaclib, k_Co55_to_p_Fe54_reaclib,
        He4, Fe54, -1, Ni58, -1, -1, Fe54, Ni58, k_Fe54_He4_to_Ni58_reaclib, k_Ni58_to_He4_Fe54_reaclib,
        H1, Fe56, -1, Co57, -1, -1, Fe56, Co57, k_Fe56_p_to_Co57_reaclib, k_Co57_to_p_Fe56_reaclib,
        He4, Fe56, -1, Ni60, -1, -1, Fe56, Ni60, k_Fe56_He4_to_Ni60_reaclib, k_Ni60_to_He4_Fe56_reaclib,
        H1, Fe58, -1, Co59, -1, -1, Fe58, Co59, k_Fe58_p_to_Co59_reaclib, k_Co59_to_p_Fe58_reaclib,
        N, Co55, -1, Co56, -1, -1, Co55, Co56, k_Co55_n_to_Co56_reaclib, k_Co56_to_n_Co55_reaclib,
        H1, Co55, -1, Ni56, -1, -1, Co55, Ni56, k_Co55_p_to_Ni56_reaclib, k_Ni56_to_p_Co55_reaclib,
        He4, Co55, -1, Cu59, -1, -1, Co55, Cu59, k_Co55_He4_to_Cu59_reaclib, k_Cu59_to_He4_Co55_reaclib,
        N, Co56, -1, Co57, -1, -1, Co56, Co57, k_Co56_n_to_Co57_reaclib, k_Co57_to_n_Co56_reaclib,
        H1, Co57, -1, Ni58, -1, -1, Co57, Ni58, k_Co57_p_to_Ni58_reaclib, k_Ni58_to_p_Co57_reaclib,
        H1, Co59, -1, Ni60, -1, -1, Co59, Ni60, k_Co59_p_to_Ni60_reaclib, k_Ni60_to_p_Co59_reaclib,
        He4, Ni56, -1, Zn60, -1, -1, Ni56, Zn60, k_Ni56_He4_to_Zn60_reaclib, k_Zn60_to_He4_Ni56_reaclib,
        H1, Ni58, -1, Cu59, -1, -1, Ni58, Cu59, k_Ni58_p_to_Cu59_reaclib, k_Cu59_to_p_Ni58_reaclib,
        H1, Cu59, -1, Zn60, -1, -1, Cu59, Zn60, k_Cu59_p_to_Zn60_reaclib, k_Zn60_to_p_Cu59_reaclib,
        N, Mg23, -1, He4, Ne20, -1, Mg23, Ne20, k_Mg23_n_to_He4_Ne20_reaclib, k_Ne20_He4_to_n_Mg23_reaclib,
        H1, Al27, -1, He4, Mg24, -1, Al27, Mg24, k_Al27_p_to_He4_Mg24_reaclib, k_Mg24_He4_to_p_Al27_reaclib,
        He4, Al27, -1, H1, Si30, -1, Al27, Si30, k_Al27_He4_to_p_Si30_reaclib, k_Si30_p_to_He4_Al27_reaclib,
        N, Si27, -1, H1, Al27, -1, Si27, Al27, k_Si27_n_to_p_Al27_reaclib, k_Al27_p_to_n_Si27_reaclib,
        N, Si27, -1, He4, Mg24, -1, Si27, Mg24, k_Si27_n_to_He4_Mg24_reaclib, k_Mg24_He4_to_n_Si27_reaclib,
        H1, Si31, -1, N, P31, -1, Si31, P31, k_Si31_p_to_n_P31_reaclib, k_P31_n_to_p_Si31_reaclib,
        He4, Si31, -1, N, S34, -1, Si31, S34, k_Si31_He4_to_n_S34_reaclib, k_S34_n_to_He4_Si31_reaclib,
        N, P29, -1, H1, Si29, -1, P29, Si29, k_P29_n_to_p_Si29_reaclib, k_Si29_p_to_n_P29_reaclib,
        He4, P29, -1, H1, S32, -1, P29, S32, k_P29_He4_to_p_S32_reaclib, k_S32_p_to_He4_P29_reaclib,
        H1, P31, -1, He4, Si28, -1, P31, Si28, k_P31_p_to_He4_Si28_reaclib, k_Si28_He4_to_p_P31_reaclib,
        He4, P31, -1, H1, S34, -1, P31, S34, k_P31_He4_to_p_S34_reaclib, k_S34_p_to_He4_P31_reaclib,
        H1, P33, -1, He4, Si30, -1, P33, Si30, k_P33_p_to_He4_Si30_reaclib, k_Si30_He4_to_p_P33_reaclib,
        N, S32, -1, He4, Si29, -1, S32, Si29, k_S32_n_to_He4_Si29_reaclib, k_Si29_He4_to_n_S32_reaclib,
        N, S33, -1, H1, P33, -1, S33, P33, k_S33_n_to_p_P33_reaclib, k_P33_p_to_n_S33_reaclib,
        N, S33, -1, He4, Si30, -1, S33, Si30, k_S33_n_to_He4_Si30_reaclib, k_Si30_He4_to_n_S33_reaclib,
        He4, S35, -1, N, Ar38, -1, S35, Ar38, k_S35_He4_to_n_Ar38_reaclib, k_Ar38_n_to_He4_S35_reaclib,
        N, Cl33, -1, H1, S33, -1, Cl33, S33, k_Cl33_n_to_p_S33_reaclib, k_S33_p_to_n_Cl33_reaclib,
        He4, Cl33, -1, H1, Ar36, -1, Cl33, Ar36, k_Cl33_He4_to_p_Ar36_reaclib, k_Ar36_p_to_He4_Cl33_reaclib,
        N, Cl35, -1, H1, S35, -1, Cl35, S35, k_Cl35_n_to_p_S35_reaclib, k_S35_p_to_n_Cl35_reaclib,
        H1, Cl35, -1, He4, S32, -1, Cl35, S32, k_Cl35_p_to_He4_S32_reaclib, k_S32_He4_to_p_Cl35_reaclib,
        He4, Cl35, -1, H1, Ar38, -1, Cl35, Ar38, k_Cl35_He4_to_p_Ar38_reaclib, k_Ar38_p_to_He4_Cl35_reaclib,
        H1, Cl37, -1, He4, S34, -1, Cl37, S34, k_Cl37_p_to_He4_S34_reaclib, k_S34_He4_to_p_Cl37_reaclib,
        N, Ar36, -1, He4, S33, -1, Ar36, S33, k_Ar36_n_to_He4_S33_reaclib, k_S33_He4_to_n_Ar36_reaclib,
        N, Ar37, -1, H1, Cl37, -1, Ar37, Cl37, k_Ar37_n_to_p_Cl37_reaclib, k_Cl37_p_to_n_Ar37_reaclib,
        N, Ar37, -1, He4, S34, -1, Ar37, S34, k_Ar37_n_to_He4_S34_reaclib, k_S34_He4_to_n_Ar37_reaclib,
        N, K37, -1, H1, Ar37, -1, K37, Ar37, k_K37_n_to_p_Ar37_reaclib, k_Ar37_p_to_n_K37_reaclib,
        He4, K37, -1, H1, Ca40, -1, K37, Ca40, k_K37_He4_to_p_Ca40_reaclib, k_Ca40_p_to_He4_K37_reaclib,
        H1, K39, -1, He4, Ar36, -1, K39, Ar36, k_K39_p_to_He4_Ar36_reaclib, k_Ar36_He4_to_p_K39_reaclib,
        N, Ca40, -1, He4, Ar37, -1, Ca40, Ar37, k_Ca40_n_to_He4_Ar37_reaclib, k_Ar37_He4_to_n_Ca40_reaclib,
        H1, Ca42, -1, He4, K39, -1, Ca42, K39, k_Ca42_p_to_He4_K39_reaclib, k_K39_He4_to_p_Ca42_reaclib,
        He4, Sc41, -1, H1, Ti44, -1, Sc41, Ti44, k_Sc41_He4_to_p_Ti44_reaclib, k_Ti44_p_to_He4_Sc41_reaclib,
        H1, Sc43, -1, He4, Ca40, -1, Sc43, Ca40, k_Sc43_p_to_He4_Ca40_reaclib, k_Ca40_He4_to_p_Sc43_reaclib,
        He4, Sc43, -1, H1, Ti46, -1, Sc43, Ti46, k_Sc43_He4_to_p_Ti46_reaclib, k_Ti46_p_to_He4_Sc43_reaclib,
        H1, Sc45, -1, He4, Ca42, -1, Sc45, Ca42, k_Sc45_p_to_He4_Ca42_reaclib, k_Ca42_He4_to_p_Sc45_reaclib,
        He4, Sc45, -1, H1, Ti48, -1, Sc45, Ti48, k_Sc45_He4_to_p_Ti48_reaclib, k_Ti48_p_to_He4_Sc45_reaclib,
        H1, V47, -1, He4, Ti44, -1, V47, Ti44, k_V47_p_to_He4_Ti44_reaclib, k_Ti44_He4_to_p_V47_reaclib,
        He4, V47, -1, H1, Cr50, -1, V47, Cr50, k_V47_He4_to_p_Cr50_reaclib, k_Cr50_p_to_He4_V47_reaclib,
        H1, V49, -1, He4, Ti46, -1, V49, Ti46, k_V49_p_to_He4_Ti46_reaclib, k_Ti46_He4_to_p_V49_reaclib,
        He4, V49, -1, H1, Cr52, -1, V49, Cr52, k_V49_He4_to_p_Cr52_reaclib, k_Cr52_p_to_He4_V49_reaclib,
        H1, V51, -1, He4, Ti48, -1, V51, Ti48, k_V51_p_to_He4_Ti48_reaclib, k_Ti48_He4_to_p_V51_reaclib,
        He4, Cr48, -1, H1, Mn51, -1, Cr48, Mn51, k_Cr48_He4_to_p_Mn51_reaclib, k_Mn51_p_to_He4_Cr48_reaclib,
        H1, Cr54, -1, He4, V51, -1, Cr54, V51, k_Cr54_p_to_He4_V51_reaclib, k_V51_He4_to_p_Cr54_reaclib,
        He4, Mn51, -1, H1, Fe54, -1, Mn51, Fe54, k_Mn51_He4_to_p_Fe54_reaclib, k_Fe54_p_to_He4_Mn51_reaclib,
        H1, Mn53, -1, He4, Cr50, -1, Mn53, Cr50, k_Mn53_p_to_He4_Cr50_reaclib, k_Cr50_He4_to_p_Mn53_reaclib,
        He4, Mn53, -1, H1, Fe56, -1, Mn53, Fe56, k_Mn53_He4_to_p_Fe56_reaclib, k_Fe56_p_to_He4_Mn53_reaclib,
        H1, Mn55, -1, He4, Cr52, -1, Mn55, Cr52, k_Mn55_p_to_He4_Cr52_reaclib, k_Cr52_He4_to_p_Mn55_reaclib,
        He4, Fe52, -1, H1, Co55, -1, Fe52, Co55, k_Fe52_He4_to_p_Co55_reaclib, k_Co55_p_to_He4_Fe52_reaclib,
        H1, Fe58, -1, He4, Mn55, -1, Fe58, Mn55, k_Fe58_p_to_He4_Mn55_reaclib, k_Mn55_He4_to_p_Fe58_reaclib,
        He4, Co55, -1, H1, Ni58, -1, Co55, Ni58, k_Co55_He4_to_p_Ni58_reaclib, k_Ni58_p_to_He4_Co55_reaclib,
        N, Co56, -1, H1, Fe56, -1, Co56, Fe56, k_Co56_n_to_p_Fe56_reaclib, k_Fe56_p_to_n_Co56_reaclib,
        N, Co56, -1, He4, Mn53, -1, Co56, Mn53, k_Co56_n_to_He4_Mn53_reaclib, k_Mn53_He4_to_n_Co56_reaclib,
        H1, Co57, -1, He4, Fe54, -1, Co57, Fe54, k_Co57_p_to_He4_Fe54_reaclib, k_Fe54_He4_to_p_Co57_reaclib,
        He4, Co57, -1, H1, Ni60, -1, Co57, Ni60, k_Co57_He4_to_p_Ni60_reaclib, k_Ni60_p_to_He4_Co57_reaclib,
        H1, Co59, -1, He4, Fe56, -1, Co59, Fe56, k_Co59_p_to_He4_Fe56_reaclib, k_Fe56_He4_to_p_Co59_reaclib,
        N, Ni56, -1, H1, Co56, -1, Ni56, Co56, k_Ni56_n_to_p_Co56_reaclib, k_Co56_p_to_n_Ni56_reaclib,
        N, Cu59, -1, He4, Co56, -1, Cu59, Co56, k_Cu59_n_to_He4_Co56_reaclib, k_Co56_He4_to_n_Cu59_reaclib,
        H1, Cu59, -1, He4, Ni56, -1, Cu59, Ni56, k_Cu59_p_to_He4_Ni56_reaclib, k_Ni56_He4_to_p_Cu59_reaclib,
        He4, He4, He4, C12, -1, -1, C12, -1, k_He4_He4_He4_to_C12_reaclib, k_C12_to_3He4_reaclib,
        Ne20, He4, -1, Mg24, -1, -1, Ne20, Mg24, k_Ne20_He4_to_Mg24_approx, k_Mg24_to_Ne20_He4_approx
    };
}
#endif

#endif
