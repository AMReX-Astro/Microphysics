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
        k_C12_He4_to_O16_reaclib = 1,
        k_O16_He4_to_Ne20_reaclib = 2,
        k_Mg23_n_to_Mg24_reaclib = 3,
        k_Mg23_He4_to_Si27_reaclib = 4,
        k_Mg24_He4_to_Si28_reaclib = 5,
        k_Al27_p_to_Si28_reaclib = 6,
        k_Al27_He4_to_P31_reaclib = 7,
        k_Si27_n_to_Si28_reaclib = 8,
        k_Si28_n_to_Si29_reaclib = 9,
        k_Si28_p_to_P29_reaclib = 10,
        k_Si28_He4_to_S32_reaclib = 11,
        k_Si29_n_to_Si30_reaclib = 12,
        k_Si29_He4_to_S33_reaclib = 13,
        k_Si30_n_to_Si31_reaclib = 14,
        k_Si30_p_to_P31_reaclib = 15,
        k_Si30_He4_to_S34_reaclib = 16,
        k_Si31_He4_to_S35_reaclib = 17,
        k_P29_He4_to_Cl33_reaclib = 18,
        k_P31_p_to_S32_reaclib = 19,
        k_P31_He4_to_Cl35_reaclib = 20,
        k_P33_p_to_S34_reaclib = 21,
        k_P33_He4_to_Cl37_reaclib = 22,
        k_S32_n_to_S33_reaclib = 23,
        k_S32_p_to_Cl33_reaclib = 24,
        k_S32_He4_to_Ar36_reaclib = 25,
        k_S33_n_to_S34_reaclib = 26,
        k_S33_He4_to_Ar37_reaclib = 27,
        k_S34_n_to_S35_reaclib = 28,
        k_S34_p_to_Cl35_reaclib = 29,
        k_S34_He4_to_Ar38_reaclib = 30,
        k_Cl33_He4_to_K37_reaclib = 31,
        k_Cl35_p_to_Ar36_reaclib = 32,
        k_Cl35_He4_to_K39_reaclib = 33,
        k_Cl37_p_to_Ar38_reaclib = 34,
        k_Ar36_n_to_Ar37_reaclib = 35,
        k_Ar36_p_to_K37_reaclib = 36,
        k_Ar36_He4_to_Ca40_reaclib = 37,
        k_Ar37_n_to_Ar38_reaclib = 38,
        k_Ar38_p_to_K39_reaclib = 39,
        k_Ar38_He4_to_Ca42_reaclib = 40,
        k_K37_He4_to_Sc41_reaclib = 41,
        k_K39_p_to_Ca40_reaclib = 42,
        k_K39_He4_to_Sc43_reaclib = 43,
        k_Ca40_p_to_Sc41_reaclib = 44,
        k_Ca40_He4_to_Ti44_reaclib = 45,
        k_Ca42_p_to_Sc43_reaclib = 46,
        k_Ca42_He4_to_Ti46_reaclib = 47,
        k_Sc43_p_to_Ti44_reaclib = 48,
        k_Sc43_He4_to_V47_reaclib = 49,
        k_Sc45_p_to_Ti46_reaclib = 50,
        k_Sc45_He4_to_V49_reaclib = 51,
        k_Ti44_He4_to_Cr48_reaclib = 52,
        k_Ti46_p_to_V47_reaclib = 53,
        k_Ti46_He4_to_Cr50_reaclib = 54,
        k_Ti48_p_to_V49_reaclib = 55,
        k_Ti48_He4_to_Cr52_reaclib = 56,
        k_V47_p_to_Cr48_reaclib = 57,
        k_V47_He4_to_Mn51_reaclib = 58,
        k_V49_p_to_Cr50_reaclib = 59,
        k_V49_He4_to_Mn53_reaclib = 60,
        k_V51_p_to_Cr52_reaclib = 61,
        k_V51_He4_to_Mn55_reaclib = 62,
        k_Cr48_He4_to_Fe52_reaclib = 63,
        k_Cr50_p_to_Mn51_reaclib = 64,
        k_Cr50_He4_to_Fe54_reaclib = 65,
        k_Cr52_p_to_Mn53_reaclib = 66,
        k_Cr52_He4_to_Fe56_reaclib = 67,
        k_Cr54_p_to_Mn55_reaclib = 68,
        k_Cr54_He4_to_Fe58_reaclib = 69,
        k_Mn51_p_to_Fe52_reaclib = 70,
        k_Mn51_He4_to_Co55_reaclib = 71,
        k_Mn53_p_to_Fe54_reaclib = 72,
        k_Mn53_He4_to_Co57_reaclib = 73,
        k_Mn55_p_to_Fe56_reaclib = 74,
        k_Mn55_He4_to_Co59_reaclib = 75,
        k_Fe52_He4_to_Ni56_reaclib = 76,
        k_Fe54_p_to_Co55_reaclib = 77,
        k_Fe54_He4_to_Ni58_reaclib = 78,
        k_Fe56_p_to_Co57_reaclib = 79,
        k_Fe56_He4_to_Ni60_reaclib = 80,
        k_Fe58_p_to_Co59_reaclib = 81,
        k_Co55_n_to_Co56_reaclib = 82,
        k_Co55_p_to_Ni56_reaclib = 83,
        k_Co55_He4_to_Cu59_reaclib = 84,
        k_Co56_n_to_Co57_reaclib = 85,
        k_Co57_p_to_Ni58_reaclib = 86,
        k_Co59_p_to_Ni60_reaclib = 87,
        k_Ni56_He4_to_Zn60_reaclib = 88,
        k_Ni58_p_to_Cu59_reaclib = 89,
        k_Cu59_p_to_Zn60_reaclib = 90,
        k_O16_C12_to_p_Al27_reaclib = 91,
        k_O16_C12_to_He4_Mg24_reaclib = 92,
        k_O16_O16_to_p_P31_reaclib = 93,
        k_O16_O16_to_He4_Si28_reaclib = 94,
        k_Ne20_C12_to_p_P31_reaclib = 95,
        k_Ne20_C12_to_He4_Si28_reaclib = 96,
        k_Mg23_n_to_He4_Ne20_reaclib = 97,
        k_Mg23_n_to_C12_C12_reaclib = 98,
        k_Al27_p_to_He4_Mg24_reaclib = 99,
        k_Al27_He4_to_p_Si30_reaclib = 100,
        k_Si27_n_to_p_Al27_reaclib = 101,
        k_Si27_n_to_He4_Mg24_reaclib = 102,
        k_Si27_n_to_C12_O16_reaclib = 103,
        k_Si31_p_to_n_P31_reaclib = 104,
        k_Si31_He4_to_n_S34_reaclib = 105,
        k_P29_n_to_p_Si29_reaclib = 106,
        k_P29_He4_to_p_S32_reaclib = 107,
        k_P31_p_to_He4_Si28_reaclib = 108,
        k_P31_He4_to_p_S34_reaclib = 109,
        k_P33_p_to_He4_Si30_reaclib = 110,
        k_S32_n_to_He4_Si29_reaclib = 111,
        k_S33_n_to_p_P33_reaclib = 112,
        k_S33_n_to_He4_Si30_reaclib = 113,
        k_S34_He4_to_n_Ar37_reaclib = 114,
        k_S34_He4_to_p_Cl37_reaclib = 115,
        k_S35_He4_to_n_Ar38_reaclib = 116,
        k_Cl33_n_to_p_S33_reaclib = 117,
        k_Cl33_He4_to_p_Ar36_reaclib = 118,
        k_Cl35_n_to_p_S35_reaclib = 119,
        k_Cl35_p_to_He4_S32_reaclib = 120,
        k_Cl35_He4_to_p_Ar38_reaclib = 121,
        k_Ar36_n_to_He4_S33_reaclib = 122,
        k_Ar37_n_to_p_Cl37_reaclib = 123,
        k_K37_n_to_p_Ar37_reaclib = 124,
        k_K37_He4_to_p_Ca40_reaclib = 125,
        k_K39_p_to_He4_Ar36_reaclib = 126,
        k_Ca40_n_to_He4_Ar37_reaclib = 127,
        k_Ca42_p_to_He4_K39_reaclib = 128,
        k_Sc41_He4_to_p_Ti44_reaclib = 129,
        k_Sc43_p_to_He4_Ca40_reaclib = 130,
        k_Sc43_He4_to_p_Ti46_reaclib = 131,
        k_Sc45_p_to_He4_Ca42_reaclib = 132,
        k_Sc45_He4_to_p_Ti48_reaclib = 133,
        k_Ti44_He4_to_p_V47_reaclib = 134,
        k_V47_He4_to_p_Cr50_reaclib = 135,
        k_V49_p_to_He4_Ti46_reaclib = 136,
        k_V49_He4_to_p_Cr52_reaclib = 137,
        k_V51_p_to_He4_Ti48_reaclib = 138,
        k_Cr48_He4_to_p_Mn51_reaclib = 139,
        k_Cr54_p_to_He4_V51_reaclib = 140,
        k_Mn51_He4_to_p_Fe54_reaclib = 141,
        k_Mn53_p_to_He4_Cr50_reaclib = 142,
        k_Mn53_He4_to_p_Fe56_reaclib = 143,
        k_Mn55_p_to_He4_Cr52_reaclib = 144,
        k_Fe52_He4_to_p_Co55_reaclib = 145,
        k_Fe58_p_to_He4_Mn55_reaclib = 146,
        k_Co55_He4_to_p_Ni58_reaclib = 147,
        k_Co56_n_to_p_Fe56_reaclib = 148,
        k_Co56_n_to_He4_Mn53_reaclib = 149,
        k_Co57_p_to_He4_Fe54_reaclib = 150,
        k_Co57_He4_to_p_Ni60_reaclib = 151,
        k_Co59_p_to_He4_Fe56_reaclib = 152,
        k_Ni56_n_to_p_Co56_reaclib = 153,
        k_Cu59_n_to_He4_Co56_reaclib = 154,
        k_Cu59_p_to_He4_Ni56_reaclib = 155,
        k_He4_He4_He4_to_C12_reaclib = 156,
        k_Ni56_n_to_Ni57_reaclib = 157,
        k_Ni57_n_to_Ni58_reaclib = 158,
        k_Ni58_n_to_Ni59_reaclib = 159,
        k_Ni59_n_to_Ni60_reaclib = 160,
        k_Co57_n_to_Co58_reaclib = 161,
        k_Co58_n_to_Co59_reaclib = 162,
        k_Fe52_n_to_Fe53_reaclib = 163,
        k_Fe53_n_to_Fe54_reaclib = 164,
        k_Fe54_n_to_Fe55_reaclib = 165,
        k_Fe55_n_to_Fe56_reaclib = 166,
        k_Fe56_n_to_Fe57_reaclib = 167,
        k_Fe57_n_to_Fe58_reaclib = 168,
        k_Mn51_n_to_Mn52_reaclib = 169,
        k_Mn52_n_to_Mn53_reaclib = 170,
        k_Mn53_n_to_Mn54_reaclib = 171,
        k_Mn54_n_to_Mn55_reaclib = 172,
        k_Cr48_n_to_Cr49_reaclib = 173,
        k_Cr49_n_to_Cr50_reaclib = 174,
        k_Cr50_n_to_Cr51_reaclib = 175,
        k_Cr51_n_to_Cr52_reaclib = 176,
        k_Cr52_n_to_Cr53_reaclib = 177,
        k_Cr53_n_to_Cr54_reaclib = 178,
        k_V47_n_to_V48_reaclib = 179,
        k_V48_n_to_V49_reaclib = 180,
        k_V49_n_to_V50_reaclib = 181,
        k_V50_n_to_V51_reaclib = 182,
        k_Ti44_n_to_Ti45_reaclib = 183,
        k_Ti45_n_to_Ti46_reaclib = 184,
        k_Ti46_n_to_Ti47_reaclib = 185,
        k_Ti47_n_to_Ti48_reaclib = 186,
        k_Sc41_n_to_Sc42_reaclib = 187,
        k_Sc42_n_to_Sc43_reaclib = 188,
        k_Sc43_n_to_Sc44_reaclib = 189,
        k_Sc44_n_to_Sc45_reaclib = 190,
        k_Ca40_n_to_Ca41_reaclib = 191,
        k_Ca41_n_to_Ca42_reaclib = 192,
        k_C12_C12_to_p_Na23_reaclib = 193,
        k_Na23_p_to_He4_Ne20_reaclib = 194,
        k_Na23_p_to_Mg24_reaclib = 195,
        k_C12_C12_to_He4_Ne20_reaclib = 196,
        k_Ne20_He4_to_Mg24_reaclib = 197,
        k_Al27_to_Si27_beta_neg_weaktab = 198,
        k_P29_to_Si29_electron_capture_weaktab = 199,
        k_Si29_to_P29_beta_neg_weaktab = 200,
        k_P31_to_Si31_electron_capture_weaktab = 201,
        k_Si31_to_P31_beta_neg_weaktab = 202,
        k_Cl33_to_S33_electron_capture_weaktab = 203,
        k_P33_to_S33_beta_neg_weaktab = 204,
        k_S33_to_Cl33_beta_neg_weaktab = 205,
        k_S33_to_P33_electron_capture_weaktab = 206,
        k_Cl35_to_S35_electron_capture_weaktab = 207,
        k_S35_to_Cl35_beta_neg_weaktab = 208,
        k_Ar37_to_Cl37_electron_capture_weaktab = 209,
        k_Ar37_to_K37_beta_neg_weaktab = 210,
        k_Cl37_to_Ar37_beta_neg_weaktab = 211,
        k_K37_to_Ar37_electron_capture_weaktab = 212,
        k_Co56_to_Fe56_electron_capture_weaktab = 213,
        k_Co56_to_Ni56_beta_neg_weaktab = 214,
        k_Fe56_to_Co56_beta_neg_weaktab = 215,
        k_Ni56_to_Co56_electron_capture_weaktab = 216,
        k_n_to_p_beta_neg_weaktab = 217,
        k_p_to_n_electron_capture_weaktab = 218,
        k_Si27_to_Al27_electron_capture_weaktab = 219,
        k_Ni56_n_n_to_Ni58_approx = 220,
        k_Ni58_to_Ni56_n_n_approx = 221,
        k_Ni58_n_n_to_Ni60_approx = 222,
        k_Ni60_to_Ni58_n_n_approx = 223,
        k_Co57_n_n_to_Co59_approx = 224,
        k_Co59_to_Co57_n_n_approx = 225,
        k_Fe52_n_n_to_Fe54_approx = 226,
        k_Fe54_to_Fe52_n_n_approx = 227,
        k_Fe54_n_n_to_Fe56_approx = 228,
        k_Fe56_to_Fe54_n_n_approx = 229,
        k_Fe56_n_n_to_Fe58_approx = 230,
        k_Fe58_to_Fe56_n_n_approx = 231,
        k_Mn51_n_n_to_Mn53_approx = 232,
        k_Mn53_to_Mn51_n_n_approx = 233,
        k_Mn53_n_n_to_Mn55_approx = 234,
        k_Mn55_to_Mn53_n_n_approx = 235,
        k_Cr48_n_n_to_Cr50_approx = 236,
        k_Cr50_to_Cr48_n_n_approx = 237,
        k_Cr50_n_n_to_Cr52_approx = 238,
        k_Cr52_to_Cr50_n_n_approx = 239,
        k_Cr52_n_n_to_Cr54_approx = 240,
        k_Cr54_to_Cr52_n_n_approx = 241,
        k_V47_n_n_to_V49_approx = 242,
        k_V49_to_V47_n_n_approx = 243,
        k_V49_n_n_to_V51_approx = 244,
        k_V51_to_V49_n_n_approx = 245,
        k_Ti44_n_n_to_Ti46_approx = 246,
        k_Ti46_to_Ti44_n_n_approx = 247,
        k_Ti46_n_n_to_Ti48_approx = 248,
        k_Ti48_to_Ti46_n_n_approx = 249,
        k_Sc41_n_n_to_Sc43_approx = 250,
        k_Sc43_to_Sc41_n_n_approx = 251,
        k_Sc43_n_n_to_Sc45_approx = 252,
        k_Sc45_to_Sc43_n_n_approx = 253,
        k_Ca40_n_n_to_Ca42_approx = 254,
        k_Ca42_to_Ca40_n_n_approx = 255,
        k_C12_C12_to_Mg24_approx = 256,
        k_Mg24_to_C12_C12_approx = 257,
        k_C12_C12_to_Ne20_He4_approx = 258,
        k_Ne20_He4_to_C12_C12_approx = 259,
        k_Ne20_He4_to_Mg24_approx = 260,
        k_Mg24_to_Ne20_He4_approx = 261,
        k_O16_to_He4_C12_derived = 262,
        k_Ne20_to_He4_O16_derived = 263,
        k_Mg24_to_n_Mg23_derived = 264,
        k_Si27_to_He4_Mg23_derived = 265,
        k_Si28_to_n_Si27_derived = 266,
        k_Si28_to_p_Al27_derived = 267,
        k_Si28_to_He4_Mg24_derived = 268,
        k_Si29_to_n_Si28_derived = 269,
        k_Si30_to_n_Si29_derived = 270,
        k_Si31_to_n_Si30_derived = 271,
        k_P29_to_p_Si28_derived = 272,
        k_P31_to_p_Si30_derived = 273,
        k_P31_to_He4_Al27_derived = 274,
        k_S32_to_p_P31_derived = 275,
        k_S32_to_He4_Si28_derived = 276,
        k_S33_to_n_S32_derived = 277,
        k_S33_to_He4_Si29_derived = 278,
        k_S34_to_n_S33_derived = 279,
        k_S34_to_p_P33_derived = 280,
        k_S34_to_He4_Si30_derived = 281,
        k_S35_to_n_S34_derived = 282,
        k_S35_to_He4_Si31_derived = 283,
        k_Cl33_to_p_S32_derived = 284,
        k_Cl33_to_He4_P29_derived = 285,
        k_Cl35_to_p_S34_derived = 286,
        k_Cl35_to_He4_P31_derived = 287,
        k_Cl37_to_He4_P33_derived = 288,
        k_Ar36_to_p_Cl35_derived = 289,
        k_Ar36_to_He4_S32_derived = 290,
        k_Ar37_to_n_Ar36_derived = 291,
        k_Ar37_to_He4_S33_derived = 292,
        k_Ar38_to_n_Ar37_derived = 293,
        k_Ar38_to_p_Cl37_derived = 294,
        k_Ar38_to_He4_S34_derived = 295,
        k_K37_to_p_Ar36_derived = 296,
        k_K37_to_He4_Cl33_derived = 297,
        k_K39_to_p_Ar38_derived = 298,
        k_K39_to_He4_Cl35_derived = 299,
        k_Ca40_to_p_K39_derived = 300,
        k_Ca40_to_He4_Ar36_derived = 301,
        k_Ca42_to_He4_Ar38_derived = 302,
        k_Sc41_to_p_Ca40_derived = 303,
        k_Sc41_to_He4_K37_derived = 304,
        k_Sc43_to_p_Ca42_derived = 305,
        k_Sc43_to_He4_K39_derived = 306,
        k_Ti44_to_p_Sc43_derived = 307,
        k_Ti44_to_He4_Ca40_derived = 308,
        k_Ti46_to_p_Sc45_derived = 309,
        k_Ti46_to_He4_Ca42_derived = 310,
        k_V47_to_p_Ti46_derived = 311,
        k_V47_to_He4_Sc43_derived = 312,
        k_V49_to_p_Ti48_derived = 313,
        k_V49_to_He4_Sc45_derived = 314,
        k_Cr48_to_p_V47_derived = 315,
        k_Cr48_to_He4_Ti44_derived = 316,
        k_Cr50_to_p_V49_derived = 317,
        k_Cr50_to_He4_Ti46_derived = 318,
        k_Cr52_to_p_V51_derived = 319,
        k_Cr52_to_He4_Ti48_derived = 320,
        k_Mn51_to_p_Cr50_derived = 321,
        k_Mn51_to_He4_V47_derived = 322,
        k_Mn53_to_p_Cr52_derived = 323,
        k_Mn53_to_He4_V49_derived = 324,
        k_Mn55_to_p_Cr54_derived = 325,
        k_Mn55_to_He4_V51_derived = 326,
        k_Fe52_to_p_Mn51_derived = 327,
        k_Fe52_to_He4_Cr48_derived = 328,
        k_Fe54_to_p_Mn53_derived = 329,
        k_Fe54_to_He4_Cr50_derived = 330,
        k_Fe56_to_p_Mn55_derived = 331,
        k_Fe56_to_He4_Cr52_derived = 332,
        k_Fe58_to_He4_Cr54_derived = 333,
        k_Co55_to_p_Fe54_derived = 334,
        k_Co55_to_He4_Mn51_derived = 335,
        k_Co56_to_n_Co55_derived = 336,
        k_Co57_to_n_Co56_derived = 337,
        k_Co57_to_p_Fe56_derived = 338,
        k_Co57_to_He4_Mn53_derived = 339,
        k_Co59_to_p_Fe58_derived = 340,
        k_Co59_to_He4_Mn55_derived = 341,
        k_Ni56_to_p_Co55_derived = 342,
        k_Ni56_to_He4_Fe52_derived = 343,
        k_Ni58_to_p_Co57_derived = 344,
        k_Ni58_to_He4_Fe54_derived = 345,
        k_Ni60_to_p_Co59_derived = 346,
        k_Ni60_to_He4_Fe56_derived = 347,
        k_Cu59_to_p_Ni58_derived = 348,
        k_Cu59_to_He4_Co55_derived = 349,
        k_Zn60_to_p_Cu59_derived = 350,
        k_Zn60_to_He4_Ni56_derived = 351,
        k_C12_to_3He4_derived = 352,
        k_C12_C12_to_n_Mg23_derived = 353,
        k_O16_C12_to_n_Si27_derived = 354,
        k_Ne20_He4_to_n_Mg23_derived = 355,
        k_Mg24_He4_to_n_Si27_derived = 356,
        k_Mg24_He4_to_p_Al27_derived = 357,
        k_Mg24_He4_to_C12_O16_derived = 358,
        k_Al27_p_to_n_Si27_derived = 359,
        k_Al27_p_to_C12_O16_derived = 360,
        k_Si28_He4_to_p_P31_derived = 361,
        k_Si28_He4_to_C12_Ne20_derived = 362,
        k_Si28_He4_to_O16_O16_derived = 363,
        k_Si29_p_to_n_P29_derived = 364,
        k_Si29_He4_to_n_S32_derived = 365,
        k_Si30_p_to_He4_Al27_derived = 366,
        k_Si30_He4_to_n_S33_derived = 367,
        k_Si30_He4_to_p_P33_derived = 368,
        k_P31_n_to_p_Si31_derived = 369,
        k_P31_p_to_C12_Ne20_derived = 370,
        k_P31_p_to_O16_O16_derived = 371,
        k_P33_p_to_n_S33_derived = 372,
        k_S32_p_to_He4_P29_derived = 373,
        k_S32_He4_to_p_Cl35_derived = 374,
        k_S33_p_to_n_Cl33_derived = 375,
        k_S33_He4_to_n_Ar36_derived = 376,
        k_S34_n_to_He4_Si31_derived = 377,
        k_S34_p_to_He4_P31_derived = 378,
        k_S35_p_to_n_Cl35_derived = 379,
        k_Cl37_p_to_n_Ar37_derived = 380,
        k_Cl37_p_to_He4_S34_derived = 381,
        k_Ar36_p_to_He4_Cl33_derived = 382,
        k_Ar36_He4_to_p_K39_derived = 383,
        k_Ar37_n_to_He4_S34_derived = 384,
        k_Ar37_p_to_n_K37_derived = 385,
        k_Ar37_He4_to_n_Ca40_derived = 386,
        k_Ar38_n_to_He4_S35_derived = 387,
        k_Ar38_p_to_He4_Cl35_derived = 388,
        k_K39_He4_to_p_Ca42_derived = 389,
        k_Ca40_p_to_He4_K37_derived = 390,
        k_Ca40_He4_to_p_Sc43_derived = 391,
        k_Ca42_He4_to_p_Sc45_derived = 392,
        k_Ti44_p_to_He4_Sc41_derived = 393,
        k_Ti46_p_to_He4_Sc43_derived = 394,
        k_Ti46_He4_to_p_V49_derived = 395,
        k_Ti48_p_to_He4_Sc45_derived = 396,
        k_Ti48_He4_to_p_V51_derived = 397,
        k_V47_p_to_He4_Ti44_derived = 398,
        k_V51_He4_to_p_Cr54_derived = 399,
        k_Cr50_p_to_He4_V47_derived = 400,
        k_Cr50_He4_to_p_Mn53_derived = 401,
        k_Cr52_p_to_He4_V49_derived = 402,
        k_Cr52_He4_to_p_Mn55_derived = 403,
        k_Mn51_p_to_He4_Cr48_derived = 404,
        k_Mn53_He4_to_n_Co56_derived = 405,
        k_Mn55_He4_to_p_Fe58_derived = 406,
        k_Fe54_p_to_He4_Mn51_derived = 407,
        k_Fe54_He4_to_p_Co57_derived = 408,
        k_Fe56_p_to_n_Co56_derived = 409,
        k_Fe56_p_to_He4_Mn53_derived = 410,
        k_Fe56_He4_to_p_Co59_derived = 411,
        k_Co55_p_to_He4_Fe52_derived = 412,
        k_Co56_p_to_n_Ni56_derived = 413,
        k_Co56_He4_to_n_Cu59_derived = 414,
        k_Ni56_He4_to_p_Cu59_derived = 415,
        k_Ni58_p_to_He4_Co55_derived = 416,
        k_Ni60_p_to_He4_Co57_derived = 417,
        k_Ni57_to_n_Ni56_derived = 418,
        k_Ni58_to_n_Ni57_derived = 419,
        k_Ni59_to_n_Ni58_derived = 420,
        k_Ni60_to_n_Ni59_derived = 421,
        k_Co58_to_n_Co57_derived = 422,
        k_Co59_to_n_Co58_derived = 423,
        k_Fe53_to_n_Fe52_derived = 424,
        k_Fe54_to_n_Fe53_derived = 425,
        k_Fe55_to_n_Fe54_derived = 426,
        k_Fe56_to_n_Fe55_derived = 427,
        k_Fe57_to_n_Fe56_derived = 428,
        k_Fe58_to_n_Fe57_derived = 429,
        k_Mn52_to_n_Mn51_derived = 430,
        k_Mn53_to_n_Mn52_derived = 431,
        k_Mn54_to_n_Mn53_derived = 432,
        k_Mn55_to_n_Mn54_derived = 433,
        k_Cr49_to_n_Cr48_derived = 434,
        k_Cr50_to_n_Cr49_derived = 435,
        k_Cr51_to_n_Cr50_derived = 436,
        k_Cr52_to_n_Cr51_derived = 437,
        k_Cr53_to_n_Cr52_derived = 438,
        k_Cr54_to_n_Cr53_derived = 439,
        k_V48_to_n_V47_derived = 440,
        k_V49_to_n_V48_derived = 441,
        k_V50_to_n_V49_derived = 442,
        k_V51_to_n_V50_derived = 443,
        k_Ti45_to_n_Ti44_derived = 444,
        k_Ti46_to_n_Ti45_derived = 445,
        k_Ti47_to_n_Ti46_derived = 446,
        k_Ti48_to_n_Ti47_derived = 447,
        k_Sc42_to_n_Sc41_derived = 448,
        k_Sc43_to_n_Sc42_derived = 449,
        k_Sc44_to_n_Sc43_derived = 450,
        k_Sc45_to_n_Sc44_derived = 451,
        k_Ca41_to_n_Ca40_derived = 452,
        k_Ca42_to_n_Ca41_derived = 453,
        k_Na23_p_to_C12_C12_derived = 454,
        k_Mg24_to_p_Na23_derived = 455,
        k_Ne20_He4_to_C12_C12_derived = 456,
        k_Ne20_He4_to_p_Na23_derived = 457,
        k_Mg24_to_He4_Ne20_derived = 458,
        NumRates = k_Mg24_to_He4_Ne20_derived
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
        "C12_He4_to_O16_reaclib",  // 1,
        "O16_He4_to_Ne20_reaclib",  // 2,
        "Mg23_n_to_Mg24_reaclib",  // 3,
        "Mg23_He4_to_Si27_reaclib",  // 4,
        "Mg24_He4_to_Si28_reaclib",  // 5,
        "Al27_p_to_Si28_reaclib",  // 6,
        "Al27_He4_to_P31_reaclib",  // 7,
        "Si27_n_to_Si28_reaclib",  // 8,
        "Si28_n_to_Si29_reaclib",  // 9,
        "Si28_p_to_P29_reaclib",  // 10,
        "Si28_He4_to_S32_reaclib",  // 11,
        "Si29_n_to_Si30_reaclib",  // 12,
        "Si29_He4_to_S33_reaclib",  // 13,
        "Si30_n_to_Si31_reaclib",  // 14,
        "Si30_p_to_P31_reaclib",  // 15,
        "Si30_He4_to_S34_reaclib",  // 16,
        "Si31_He4_to_S35_reaclib",  // 17,
        "P29_He4_to_Cl33_reaclib",  // 18,
        "P31_p_to_S32_reaclib",  // 19,
        "P31_He4_to_Cl35_reaclib",  // 20,
        "P33_p_to_S34_reaclib",  // 21,
        "P33_He4_to_Cl37_reaclib",  // 22,
        "S32_n_to_S33_reaclib",  // 23,
        "S32_p_to_Cl33_reaclib",  // 24,
        "S32_He4_to_Ar36_reaclib",  // 25,
        "S33_n_to_S34_reaclib",  // 26,
        "S33_He4_to_Ar37_reaclib",  // 27,
        "S34_n_to_S35_reaclib",  // 28,
        "S34_p_to_Cl35_reaclib",  // 29,
        "S34_He4_to_Ar38_reaclib",  // 30,
        "Cl33_He4_to_K37_reaclib",  // 31,
        "Cl35_p_to_Ar36_reaclib",  // 32,
        "Cl35_He4_to_K39_reaclib",  // 33,
        "Cl37_p_to_Ar38_reaclib",  // 34,
        "Ar36_n_to_Ar37_reaclib",  // 35,
        "Ar36_p_to_K37_reaclib",  // 36,
        "Ar36_He4_to_Ca40_reaclib",  // 37,
        "Ar37_n_to_Ar38_reaclib",  // 38,
        "Ar38_p_to_K39_reaclib",  // 39,
        "Ar38_He4_to_Ca42_reaclib",  // 40,
        "K37_He4_to_Sc41_reaclib",  // 41,
        "K39_p_to_Ca40_reaclib",  // 42,
        "K39_He4_to_Sc43_reaclib",  // 43,
        "Ca40_p_to_Sc41_reaclib",  // 44,
        "Ca40_He4_to_Ti44_reaclib",  // 45,
        "Ca42_p_to_Sc43_reaclib",  // 46,
        "Ca42_He4_to_Ti46_reaclib",  // 47,
        "Sc43_p_to_Ti44_reaclib",  // 48,
        "Sc43_He4_to_V47_reaclib",  // 49,
        "Sc45_p_to_Ti46_reaclib",  // 50,
        "Sc45_He4_to_V49_reaclib",  // 51,
        "Ti44_He4_to_Cr48_reaclib",  // 52,
        "Ti46_p_to_V47_reaclib",  // 53,
        "Ti46_He4_to_Cr50_reaclib",  // 54,
        "Ti48_p_to_V49_reaclib",  // 55,
        "Ti48_He4_to_Cr52_reaclib",  // 56,
        "V47_p_to_Cr48_reaclib",  // 57,
        "V47_He4_to_Mn51_reaclib",  // 58,
        "V49_p_to_Cr50_reaclib",  // 59,
        "V49_He4_to_Mn53_reaclib",  // 60,
        "V51_p_to_Cr52_reaclib",  // 61,
        "V51_He4_to_Mn55_reaclib",  // 62,
        "Cr48_He4_to_Fe52_reaclib",  // 63,
        "Cr50_p_to_Mn51_reaclib",  // 64,
        "Cr50_He4_to_Fe54_reaclib",  // 65,
        "Cr52_p_to_Mn53_reaclib",  // 66,
        "Cr52_He4_to_Fe56_reaclib",  // 67,
        "Cr54_p_to_Mn55_reaclib",  // 68,
        "Cr54_He4_to_Fe58_reaclib",  // 69,
        "Mn51_p_to_Fe52_reaclib",  // 70,
        "Mn51_He4_to_Co55_reaclib",  // 71,
        "Mn53_p_to_Fe54_reaclib",  // 72,
        "Mn53_He4_to_Co57_reaclib",  // 73,
        "Mn55_p_to_Fe56_reaclib",  // 74,
        "Mn55_He4_to_Co59_reaclib",  // 75,
        "Fe52_He4_to_Ni56_reaclib",  // 76,
        "Fe54_p_to_Co55_reaclib",  // 77,
        "Fe54_He4_to_Ni58_reaclib",  // 78,
        "Fe56_p_to_Co57_reaclib",  // 79,
        "Fe56_He4_to_Ni60_reaclib",  // 80,
        "Fe58_p_to_Co59_reaclib",  // 81,
        "Co55_n_to_Co56_reaclib",  // 82,
        "Co55_p_to_Ni56_reaclib",  // 83,
        "Co55_He4_to_Cu59_reaclib",  // 84,
        "Co56_n_to_Co57_reaclib",  // 85,
        "Co57_p_to_Ni58_reaclib",  // 86,
        "Co59_p_to_Ni60_reaclib",  // 87,
        "Ni56_He4_to_Zn60_reaclib",  // 88,
        "Ni58_p_to_Cu59_reaclib",  // 89,
        "Cu59_p_to_Zn60_reaclib",  // 90,
        "O16_C12_to_p_Al27_reaclib",  // 91,
        "O16_C12_to_He4_Mg24_reaclib",  // 92,
        "O16_O16_to_p_P31_reaclib",  // 93,
        "O16_O16_to_He4_Si28_reaclib",  // 94,
        "Ne20_C12_to_p_P31_reaclib",  // 95,
        "Ne20_C12_to_He4_Si28_reaclib",  // 96,
        "Mg23_n_to_He4_Ne20_reaclib",  // 97,
        "Mg23_n_to_C12_C12_reaclib",  // 98,
        "Al27_p_to_He4_Mg24_reaclib",  // 99,
        "Al27_He4_to_p_Si30_reaclib",  // 100,
        "Si27_n_to_p_Al27_reaclib",  // 101,
        "Si27_n_to_He4_Mg24_reaclib",  // 102,
        "Si27_n_to_C12_O16_reaclib",  // 103,
        "Si31_p_to_n_P31_reaclib",  // 104,
        "Si31_He4_to_n_S34_reaclib",  // 105,
        "P29_n_to_p_Si29_reaclib",  // 106,
        "P29_He4_to_p_S32_reaclib",  // 107,
        "P31_p_to_He4_Si28_reaclib",  // 108,
        "P31_He4_to_p_S34_reaclib",  // 109,
        "P33_p_to_He4_Si30_reaclib",  // 110,
        "S32_n_to_He4_Si29_reaclib",  // 111,
        "S33_n_to_p_P33_reaclib",  // 112,
        "S33_n_to_He4_Si30_reaclib",  // 113,
        "S34_He4_to_n_Ar37_reaclib",  // 114,
        "S34_He4_to_p_Cl37_reaclib",  // 115,
        "S35_He4_to_n_Ar38_reaclib",  // 116,
        "Cl33_n_to_p_S33_reaclib",  // 117,
        "Cl33_He4_to_p_Ar36_reaclib",  // 118,
        "Cl35_n_to_p_S35_reaclib",  // 119,
        "Cl35_p_to_He4_S32_reaclib",  // 120,
        "Cl35_He4_to_p_Ar38_reaclib",  // 121,
        "Ar36_n_to_He4_S33_reaclib",  // 122,
        "Ar37_n_to_p_Cl37_reaclib",  // 123,
        "K37_n_to_p_Ar37_reaclib",  // 124,
        "K37_He4_to_p_Ca40_reaclib",  // 125,
        "K39_p_to_He4_Ar36_reaclib",  // 126,
        "Ca40_n_to_He4_Ar37_reaclib",  // 127,
        "Ca42_p_to_He4_K39_reaclib",  // 128,
        "Sc41_He4_to_p_Ti44_reaclib",  // 129,
        "Sc43_p_to_He4_Ca40_reaclib",  // 130,
        "Sc43_He4_to_p_Ti46_reaclib",  // 131,
        "Sc45_p_to_He4_Ca42_reaclib",  // 132,
        "Sc45_He4_to_p_Ti48_reaclib",  // 133,
        "Ti44_He4_to_p_V47_reaclib",  // 134,
        "V47_He4_to_p_Cr50_reaclib",  // 135,
        "V49_p_to_He4_Ti46_reaclib",  // 136,
        "V49_He4_to_p_Cr52_reaclib",  // 137,
        "V51_p_to_He4_Ti48_reaclib",  // 138,
        "Cr48_He4_to_p_Mn51_reaclib",  // 139,
        "Cr54_p_to_He4_V51_reaclib",  // 140,
        "Mn51_He4_to_p_Fe54_reaclib",  // 141,
        "Mn53_p_to_He4_Cr50_reaclib",  // 142,
        "Mn53_He4_to_p_Fe56_reaclib",  // 143,
        "Mn55_p_to_He4_Cr52_reaclib",  // 144,
        "Fe52_He4_to_p_Co55_reaclib",  // 145,
        "Fe58_p_to_He4_Mn55_reaclib",  // 146,
        "Co55_He4_to_p_Ni58_reaclib",  // 147,
        "Co56_n_to_p_Fe56_reaclib",  // 148,
        "Co56_n_to_He4_Mn53_reaclib",  // 149,
        "Co57_p_to_He4_Fe54_reaclib",  // 150,
        "Co57_He4_to_p_Ni60_reaclib",  // 151,
        "Co59_p_to_He4_Fe56_reaclib",  // 152,
        "Ni56_n_to_p_Co56_reaclib",  // 153,
        "Cu59_n_to_He4_Co56_reaclib",  // 154,
        "Cu59_p_to_He4_Ni56_reaclib",  // 155,
        "He4_He4_He4_to_C12_reaclib",  // 156,
        "Ni56_n_to_Ni57_reaclib",  // 157,
        "Ni57_n_to_Ni58_reaclib",  // 158,
        "Ni58_n_to_Ni59_reaclib",  // 159,
        "Ni59_n_to_Ni60_reaclib",  // 160,
        "Co57_n_to_Co58_reaclib",  // 161,
        "Co58_n_to_Co59_reaclib",  // 162,
        "Fe52_n_to_Fe53_reaclib",  // 163,
        "Fe53_n_to_Fe54_reaclib",  // 164,
        "Fe54_n_to_Fe55_reaclib",  // 165,
        "Fe55_n_to_Fe56_reaclib",  // 166,
        "Fe56_n_to_Fe57_reaclib",  // 167,
        "Fe57_n_to_Fe58_reaclib",  // 168,
        "Mn51_n_to_Mn52_reaclib",  // 169,
        "Mn52_n_to_Mn53_reaclib",  // 170,
        "Mn53_n_to_Mn54_reaclib",  // 171,
        "Mn54_n_to_Mn55_reaclib",  // 172,
        "Cr48_n_to_Cr49_reaclib",  // 173,
        "Cr49_n_to_Cr50_reaclib",  // 174,
        "Cr50_n_to_Cr51_reaclib",  // 175,
        "Cr51_n_to_Cr52_reaclib",  // 176,
        "Cr52_n_to_Cr53_reaclib",  // 177,
        "Cr53_n_to_Cr54_reaclib",  // 178,
        "V47_n_to_V48_reaclib",  // 179,
        "V48_n_to_V49_reaclib",  // 180,
        "V49_n_to_V50_reaclib",  // 181,
        "V50_n_to_V51_reaclib",  // 182,
        "Ti44_n_to_Ti45_reaclib",  // 183,
        "Ti45_n_to_Ti46_reaclib",  // 184,
        "Ti46_n_to_Ti47_reaclib",  // 185,
        "Ti47_n_to_Ti48_reaclib",  // 186,
        "Sc41_n_to_Sc42_reaclib",  // 187,
        "Sc42_n_to_Sc43_reaclib",  // 188,
        "Sc43_n_to_Sc44_reaclib",  // 189,
        "Sc44_n_to_Sc45_reaclib",  // 190,
        "Ca40_n_to_Ca41_reaclib",  // 191,
        "Ca41_n_to_Ca42_reaclib",  // 192,
        "C12_C12_to_p_Na23_reaclib",  // 193,
        "Na23_p_to_He4_Ne20_reaclib",  // 194,
        "Na23_p_to_Mg24_reaclib",  // 195,
        "C12_C12_to_He4_Ne20_reaclib",  // 196,
        "Ne20_He4_to_Mg24_reaclib",  // 197,
        "Al27_to_Si27_beta_neg_weaktab",  // 198,
        "P29_to_Si29_electron_capture_weaktab",  // 199,
        "Si29_to_P29_beta_neg_weaktab",  // 200,
        "P31_to_Si31_electron_capture_weaktab",  // 201,
        "Si31_to_P31_beta_neg_weaktab",  // 202,
        "Cl33_to_S33_electron_capture_weaktab",  // 203,
        "P33_to_S33_beta_neg_weaktab",  // 204,
        "S33_to_Cl33_beta_neg_weaktab",  // 205,
        "S33_to_P33_electron_capture_weaktab",  // 206,
        "Cl35_to_S35_electron_capture_weaktab",  // 207,
        "S35_to_Cl35_beta_neg_weaktab",  // 208,
        "Ar37_to_Cl37_electron_capture_weaktab",  // 209,
        "Ar37_to_K37_beta_neg_weaktab",  // 210,
        "Cl37_to_Ar37_beta_neg_weaktab",  // 211,
        "K37_to_Ar37_electron_capture_weaktab",  // 212,
        "Co56_to_Fe56_electron_capture_weaktab",  // 213,
        "Co56_to_Ni56_beta_neg_weaktab",  // 214,
        "Fe56_to_Co56_beta_neg_weaktab",  // 215,
        "Ni56_to_Co56_electron_capture_weaktab",  // 216,
        "n_to_p_beta_neg_weaktab",  // 217,
        "p_to_n_electron_capture_weaktab",  // 218,
        "Si27_to_Al27_electron_capture_weaktab",  // 219,
        "Ni56_n_n_to_Ni58_approx",  // 220,
        "Ni58_to_Ni56_n_n_approx",  // 221,
        "Ni58_n_n_to_Ni60_approx",  // 222,
        "Ni60_to_Ni58_n_n_approx",  // 223,
        "Co57_n_n_to_Co59_approx",  // 224,
        "Co59_to_Co57_n_n_approx",  // 225,
        "Fe52_n_n_to_Fe54_approx",  // 226,
        "Fe54_to_Fe52_n_n_approx",  // 227,
        "Fe54_n_n_to_Fe56_approx",  // 228,
        "Fe56_to_Fe54_n_n_approx",  // 229,
        "Fe56_n_n_to_Fe58_approx",  // 230,
        "Fe58_to_Fe56_n_n_approx",  // 231,
        "Mn51_n_n_to_Mn53_approx",  // 232,
        "Mn53_to_Mn51_n_n_approx",  // 233,
        "Mn53_n_n_to_Mn55_approx",  // 234,
        "Mn55_to_Mn53_n_n_approx",  // 235,
        "Cr48_n_n_to_Cr50_approx",  // 236,
        "Cr50_to_Cr48_n_n_approx",  // 237,
        "Cr50_n_n_to_Cr52_approx",  // 238,
        "Cr52_to_Cr50_n_n_approx",  // 239,
        "Cr52_n_n_to_Cr54_approx",  // 240,
        "Cr54_to_Cr52_n_n_approx",  // 241,
        "V47_n_n_to_V49_approx",  // 242,
        "V49_to_V47_n_n_approx",  // 243,
        "V49_n_n_to_V51_approx",  // 244,
        "V51_to_V49_n_n_approx",  // 245,
        "Ti44_n_n_to_Ti46_approx",  // 246,
        "Ti46_to_Ti44_n_n_approx",  // 247,
        "Ti46_n_n_to_Ti48_approx",  // 248,
        "Ti48_to_Ti46_n_n_approx",  // 249,
        "Sc41_n_n_to_Sc43_approx",  // 250,
        "Sc43_to_Sc41_n_n_approx",  // 251,
        "Sc43_n_n_to_Sc45_approx",  // 252,
        "Sc45_to_Sc43_n_n_approx",  // 253,
        "Ca40_n_n_to_Ca42_approx",  // 254,
        "Ca42_to_Ca40_n_n_approx",  // 255,
        "C12_C12_to_Mg24_approx",  // 256,
        "Mg24_to_C12_C12_approx",  // 257,
        "C12_C12_to_Ne20_He4_approx",  // 258,
        "Ne20_He4_to_C12_C12_approx",  // 259,
        "Ne20_He4_to_Mg24_approx",  // 260,
        "Mg24_to_Ne20_He4_approx",  // 261,
        "O16_to_He4_C12_derived",  // 262,
        "Ne20_to_He4_O16_derived",  // 263,
        "Mg24_to_n_Mg23_derived",  // 264,
        "Si27_to_He4_Mg23_derived",  // 265,
        "Si28_to_n_Si27_derived",  // 266,
        "Si28_to_p_Al27_derived",  // 267,
        "Si28_to_He4_Mg24_derived",  // 268,
        "Si29_to_n_Si28_derived",  // 269,
        "Si30_to_n_Si29_derived",  // 270,
        "Si31_to_n_Si30_derived",  // 271,
        "P29_to_p_Si28_derived",  // 272,
        "P31_to_p_Si30_derived",  // 273,
        "P31_to_He4_Al27_derived",  // 274,
        "S32_to_p_P31_derived",  // 275,
        "S32_to_He4_Si28_derived",  // 276,
        "S33_to_n_S32_derived",  // 277,
        "S33_to_He4_Si29_derived",  // 278,
        "S34_to_n_S33_derived",  // 279,
        "S34_to_p_P33_derived",  // 280,
        "S34_to_He4_Si30_derived",  // 281,
        "S35_to_n_S34_derived",  // 282,
        "S35_to_He4_Si31_derived",  // 283,
        "Cl33_to_p_S32_derived",  // 284,
        "Cl33_to_He4_P29_derived",  // 285,
        "Cl35_to_p_S34_derived",  // 286,
        "Cl35_to_He4_P31_derived",  // 287,
        "Cl37_to_He4_P33_derived",  // 288,
        "Ar36_to_p_Cl35_derived",  // 289,
        "Ar36_to_He4_S32_derived",  // 290,
        "Ar37_to_n_Ar36_derived",  // 291,
        "Ar37_to_He4_S33_derived",  // 292,
        "Ar38_to_n_Ar37_derived",  // 293,
        "Ar38_to_p_Cl37_derived",  // 294,
        "Ar38_to_He4_S34_derived",  // 295,
        "K37_to_p_Ar36_derived",  // 296,
        "K37_to_He4_Cl33_derived",  // 297,
        "K39_to_p_Ar38_derived",  // 298,
        "K39_to_He4_Cl35_derived",  // 299,
        "Ca40_to_p_K39_derived",  // 300,
        "Ca40_to_He4_Ar36_derived",  // 301,
        "Ca42_to_He4_Ar38_derived",  // 302,
        "Sc41_to_p_Ca40_derived",  // 303,
        "Sc41_to_He4_K37_derived",  // 304,
        "Sc43_to_p_Ca42_derived",  // 305,
        "Sc43_to_He4_K39_derived",  // 306,
        "Ti44_to_p_Sc43_derived",  // 307,
        "Ti44_to_He4_Ca40_derived",  // 308,
        "Ti46_to_p_Sc45_derived",  // 309,
        "Ti46_to_He4_Ca42_derived",  // 310,
        "V47_to_p_Ti46_derived",  // 311,
        "V47_to_He4_Sc43_derived",  // 312,
        "V49_to_p_Ti48_derived",  // 313,
        "V49_to_He4_Sc45_derived",  // 314,
        "Cr48_to_p_V47_derived",  // 315,
        "Cr48_to_He4_Ti44_derived",  // 316,
        "Cr50_to_p_V49_derived",  // 317,
        "Cr50_to_He4_Ti46_derived",  // 318,
        "Cr52_to_p_V51_derived",  // 319,
        "Cr52_to_He4_Ti48_derived",  // 320,
        "Mn51_to_p_Cr50_derived",  // 321,
        "Mn51_to_He4_V47_derived",  // 322,
        "Mn53_to_p_Cr52_derived",  // 323,
        "Mn53_to_He4_V49_derived",  // 324,
        "Mn55_to_p_Cr54_derived",  // 325,
        "Mn55_to_He4_V51_derived",  // 326,
        "Fe52_to_p_Mn51_derived",  // 327,
        "Fe52_to_He4_Cr48_derived",  // 328,
        "Fe54_to_p_Mn53_derived",  // 329,
        "Fe54_to_He4_Cr50_derived",  // 330,
        "Fe56_to_p_Mn55_derived",  // 331,
        "Fe56_to_He4_Cr52_derived",  // 332,
        "Fe58_to_He4_Cr54_derived",  // 333,
        "Co55_to_p_Fe54_derived",  // 334,
        "Co55_to_He4_Mn51_derived",  // 335,
        "Co56_to_n_Co55_derived",  // 336,
        "Co57_to_n_Co56_derived",  // 337,
        "Co57_to_p_Fe56_derived",  // 338,
        "Co57_to_He4_Mn53_derived",  // 339,
        "Co59_to_p_Fe58_derived",  // 340,
        "Co59_to_He4_Mn55_derived",  // 341,
        "Ni56_to_p_Co55_derived",  // 342,
        "Ni56_to_He4_Fe52_derived",  // 343,
        "Ni58_to_p_Co57_derived",  // 344,
        "Ni58_to_He4_Fe54_derived",  // 345,
        "Ni60_to_p_Co59_derived",  // 346,
        "Ni60_to_He4_Fe56_derived",  // 347,
        "Cu59_to_p_Ni58_derived",  // 348,
        "Cu59_to_He4_Co55_derived",  // 349,
        "Zn60_to_p_Cu59_derived",  // 350,
        "Zn60_to_He4_Ni56_derived",  // 351,
        "C12_to_3He4_derived",  // 352,
        "C12_C12_to_n_Mg23_derived",  // 353,
        "O16_C12_to_n_Si27_derived",  // 354,
        "Ne20_He4_to_n_Mg23_derived",  // 355,
        "Mg24_He4_to_n_Si27_derived",  // 356,
        "Mg24_He4_to_p_Al27_derived",  // 357,
        "Mg24_He4_to_C12_O16_derived",  // 358,
        "Al27_p_to_n_Si27_derived",  // 359,
        "Al27_p_to_C12_O16_derived",  // 360,
        "Si28_He4_to_p_P31_derived",  // 361,
        "Si28_He4_to_C12_Ne20_derived",  // 362,
        "Si28_He4_to_O16_O16_derived",  // 363,
        "Si29_p_to_n_P29_derived",  // 364,
        "Si29_He4_to_n_S32_derived",  // 365,
        "Si30_p_to_He4_Al27_derived",  // 366,
        "Si30_He4_to_n_S33_derived",  // 367,
        "Si30_He4_to_p_P33_derived",  // 368,
        "P31_n_to_p_Si31_derived",  // 369,
        "P31_p_to_C12_Ne20_derived",  // 370,
        "P31_p_to_O16_O16_derived",  // 371,
        "P33_p_to_n_S33_derived",  // 372,
        "S32_p_to_He4_P29_derived",  // 373,
        "S32_He4_to_p_Cl35_derived",  // 374,
        "S33_p_to_n_Cl33_derived",  // 375,
        "S33_He4_to_n_Ar36_derived",  // 376,
        "S34_n_to_He4_Si31_derived",  // 377,
        "S34_p_to_He4_P31_derived",  // 378,
        "S35_p_to_n_Cl35_derived",  // 379,
        "Cl37_p_to_n_Ar37_derived",  // 380,
        "Cl37_p_to_He4_S34_derived",  // 381,
        "Ar36_p_to_He4_Cl33_derived",  // 382,
        "Ar36_He4_to_p_K39_derived",  // 383,
        "Ar37_n_to_He4_S34_derived",  // 384,
        "Ar37_p_to_n_K37_derived",  // 385,
        "Ar37_He4_to_n_Ca40_derived",  // 386,
        "Ar38_n_to_He4_S35_derived",  // 387,
        "Ar38_p_to_He4_Cl35_derived",  // 388,
        "K39_He4_to_p_Ca42_derived",  // 389,
        "Ca40_p_to_He4_K37_derived",  // 390,
        "Ca40_He4_to_p_Sc43_derived",  // 391,
        "Ca42_He4_to_p_Sc45_derived",  // 392,
        "Ti44_p_to_He4_Sc41_derived",  // 393,
        "Ti46_p_to_He4_Sc43_derived",  // 394,
        "Ti46_He4_to_p_V49_derived",  // 395,
        "Ti48_p_to_He4_Sc45_derived",  // 396,
        "Ti48_He4_to_p_V51_derived",  // 397,
        "V47_p_to_He4_Ti44_derived",  // 398,
        "V51_He4_to_p_Cr54_derived",  // 399,
        "Cr50_p_to_He4_V47_derived",  // 400,
        "Cr50_He4_to_p_Mn53_derived",  // 401,
        "Cr52_p_to_He4_V49_derived",  // 402,
        "Cr52_He4_to_p_Mn55_derived",  // 403,
        "Mn51_p_to_He4_Cr48_derived",  // 404,
        "Mn53_He4_to_n_Co56_derived",  // 405,
        "Mn55_He4_to_p_Fe58_derived",  // 406,
        "Fe54_p_to_He4_Mn51_derived",  // 407,
        "Fe54_He4_to_p_Co57_derived",  // 408,
        "Fe56_p_to_n_Co56_derived",  // 409,
        "Fe56_p_to_He4_Mn53_derived",  // 410,
        "Fe56_He4_to_p_Co59_derived",  // 411,
        "Co55_p_to_He4_Fe52_derived",  // 412,
        "Co56_p_to_n_Ni56_derived",  // 413,
        "Co56_He4_to_n_Cu59_derived",  // 414,
        "Ni56_He4_to_p_Cu59_derived",  // 415,
        "Ni58_p_to_He4_Co55_derived",  // 416,
        "Ni60_p_to_He4_Co57_derived",  // 417,
        "Ni57_to_n_Ni56_derived",  // 418,
        "Ni58_to_n_Ni57_derived",  // 419,
        "Ni59_to_n_Ni58_derived",  // 420,
        "Ni60_to_n_Ni59_derived",  // 421,
        "Co58_to_n_Co57_derived",  // 422,
        "Co59_to_n_Co58_derived",  // 423,
        "Fe53_to_n_Fe52_derived",  // 424,
        "Fe54_to_n_Fe53_derived",  // 425,
        "Fe55_to_n_Fe54_derived",  // 426,
        "Fe56_to_n_Fe55_derived",  // 427,
        "Fe57_to_n_Fe56_derived",  // 428,
        "Fe58_to_n_Fe57_derived",  // 429,
        "Mn52_to_n_Mn51_derived",  // 430,
        "Mn53_to_n_Mn52_derived",  // 431,
        "Mn54_to_n_Mn53_derived",  // 432,
        "Mn55_to_n_Mn54_derived",  // 433,
        "Cr49_to_n_Cr48_derived",  // 434,
        "Cr50_to_n_Cr49_derived",  // 435,
        "Cr51_to_n_Cr50_derived",  // 436,
        "Cr52_to_n_Cr51_derived",  // 437,
        "Cr53_to_n_Cr52_derived",  // 438,
        "Cr54_to_n_Cr53_derived",  // 439,
        "V48_to_n_V47_derived",  // 440,
        "V49_to_n_V48_derived",  // 441,
        "V50_to_n_V49_derived",  // 442,
        "V51_to_n_V50_derived",  // 443,
        "Ti45_to_n_Ti44_derived",  // 444,
        "Ti46_to_n_Ti45_derived",  // 445,
        "Ti47_to_n_Ti46_derived",  // 446,
        "Ti48_to_n_Ti47_derived",  // 447,
        "Sc42_to_n_Sc41_derived",  // 448,
        "Sc43_to_n_Sc42_derived",  // 449,
        "Sc44_to_n_Sc43_derived",  // 450,
        "Sc45_to_n_Sc44_derived",  // 451,
        "Ca41_to_n_Ca40_derived",  // 452,
        "Ca42_to_n_Ca41_derived",  // 453,
        "Na23_p_to_C12_C12_derived",  // 454,
        "Mg24_to_p_Na23_derived",  // 455,
        "Ne20_He4_to_C12_C12_derived",  // 456,
        "Ne20_He4_to_p_Na23_derived",  // 457,
        "Mg24_to_He4_Ne20_derived"  // 458,
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
        He4, C12, -1, O16, -1, -1, C12, O16, k_C12_He4_to_O16_reaclib, k_O16_to_He4_C12_derived,
        He4, O16, -1, Ne20, -1, -1, O16, Ne20, k_O16_He4_to_Ne20_reaclib, k_Ne20_to_He4_O16_derived,
        N, Mg23, -1, Mg24, -1, -1, Mg23, Mg24, k_Mg23_n_to_Mg24_reaclib, k_Mg24_to_n_Mg23_derived,
        He4, Mg23, -1, Si27, -1, -1, Mg23, Si27, k_Mg23_He4_to_Si27_reaclib, k_Si27_to_He4_Mg23_derived,
        He4, Mg24, -1, Si28, -1, -1, Mg24, Si28, k_Mg24_He4_to_Si28_reaclib, k_Si28_to_He4_Mg24_derived,
        H1, Al27, -1, Si28, -1, -1, Al27, Si28, k_Al27_p_to_Si28_reaclib, k_Si28_to_p_Al27_derived,
        He4, Al27, -1, P31, -1, -1, Al27, P31, k_Al27_He4_to_P31_reaclib, k_P31_to_He4_Al27_derived,
        N, Si27, -1, Si28, -1, -1, Si27, Si28, k_Si27_n_to_Si28_reaclib, k_Si28_to_n_Si27_derived,
        N, Si28, -1, Si29, -1, -1, Si28, Si29, k_Si28_n_to_Si29_reaclib, k_Si29_to_n_Si28_derived,
        H1, Si28, -1, P29, -1, -1, Si28, P29, k_Si28_p_to_P29_reaclib, k_P29_to_p_Si28_derived,
        He4, Si28, -1, S32, -1, -1, Si28, S32, k_Si28_He4_to_S32_reaclib, k_S32_to_He4_Si28_derived,
        N, Si29, -1, Si30, -1, -1, Si29, Si30, k_Si29_n_to_Si30_reaclib, k_Si30_to_n_Si29_derived,
        He4, Si29, -1, S33, -1, -1, Si29, S33, k_Si29_He4_to_S33_reaclib, k_S33_to_He4_Si29_derived,
        N, Si30, -1, Si31, -1, -1, Si30, Si31, k_Si30_n_to_Si31_reaclib, k_Si31_to_n_Si30_derived,
        H1, Si30, -1, P31, -1, -1, Si30, P31, k_Si30_p_to_P31_reaclib, k_P31_to_p_Si30_derived,
        He4, Si30, -1, S34, -1, -1, Si30, S34, k_Si30_He4_to_S34_reaclib, k_S34_to_He4_Si30_derived,
        He4, Si31, -1, S35, -1, -1, Si31, S35, k_Si31_He4_to_S35_reaclib, k_S35_to_He4_Si31_derived,
        He4, P29, -1, Cl33, -1, -1, P29, Cl33, k_P29_He4_to_Cl33_reaclib, k_Cl33_to_He4_P29_derived,
        H1, P31, -1, S32, -1, -1, P31, S32, k_P31_p_to_S32_reaclib, k_S32_to_p_P31_derived,
        He4, P31, -1, Cl35, -1, -1, P31, Cl35, k_P31_He4_to_Cl35_reaclib, k_Cl35_to_He4_P31_derived,
        H1, P33, -1, S34, -1, -1, P33, S34, k_P33_p_to_S34_reaclib, k_S34_to_p_P33_derived,
        He4, P33, -1, Cl37, -1, -1, P33, Cl37, k_P33_He4_to_Cl37_reaclib, k_Cl37_to_He4_P33_derived,
        N, S32, -1, S33, -1, -1, S32, S33, k_S32_n_to_S33_reaclib, k_S33_to_n_S32_derived,
        H1, S32, -1, Cl33, -1, -1, S32, Cl33, k_S32_p_to_Cl33_reaclib, k_Cl33_to_p_S32_derived,
        He4, S32, -1, Ar36, -1, -1, S32, Ar36, k_S32_He4_to_Ar36_reaclib, k_Ar36_to_He4_S32_derived,
        N, S33, -1, S34, -1, -1, S33, S34, k_S33_n_to_S34_reaclib, k_S34_to_n_S33_derived,
        He4, S33, -1, Ar37, -1, -1, S33, Ar37, k_S33_He4_to_Ar37_reaclib, k_Ar37_to_He4_S33_derived,
        N, S34, -1, S35, -1, -1, S34, S35, k_S34_n_to_S35_reaclib, k_S35_to_n_S34_derived,
        H1, S34, -1, Cl35, -1, -1, S34, Cl35, k_S34_p_to_Cl35_reaclib, k_Cl35_to_p_S34_derived,
        He4, S34, -1, Ar38, -1, -1, S34, Ar38, k_S34_He4_to_Ar38_reaclib, k_Ar38_to_He4_S34_derived,
        He4, Cl33, -1, K37, -1, -1, Cl33, K37, k_Cl33_He4_to_K37_reaclib, k_K37_to_He4_Cl33_derived,
        H1, Cl35, -1, Ar36, -1, -1, Cl35, Ar36, k_Cl35_p_to_Ar36_reaclib, k_Ar36_to_p_Cl35_derived,
        He4, Cl35, -1, K39, -1, -1, Cl35, K39, k_Cl35_He4_to_K39_reaclib, k_K39_to_He4_Cl35_derived,
        H1, Cl37, -1, Ar38, -1, -1, Cl37, Ar38, k_Cl37_p_to_Ar38_reaclib, k_Ar38_to_p_Cl37_derived,
        N, Ar36, -1, Ar37, -1, -1, Ar36, Ar37, k_Ar36_n_to_Ar37_reaclib, k_Ar37_to_n_Ar36_derived,
        H1, Ar36, -1, K37, -1, -1, Ar36, K37, k_Ar36_p_to_K37_reaclib, k_K37_to_p_Ar36_derived,
        He4, Ar36, -1, Ca40, -1, -1, Ar36, Ca40, k_Ar36_He4_to_Ca40_reaclib, k_Ca40_to_He4_Ar36_derived,
        N, Ar37, -1, Ar38, -1, -1, Ar37, Ar38, k_Ar37_n_to_Ar38_reaclib, k_Ar38_to_n_Ar37_derived,
        H1, Ar38, -1, K39, -1, -1, Ar38, K39, k_Ar38_p_to_K39_reaclib, k_K39_to_p_Ar38_derived,
        He4, Ar38, -1, Ca42, -1, -1, Ar38, Ca42, k_Ar38_He4_to_Ca42_reaclib, k_Ca42_to_He4_Ar38_derived,
        He4, K37, -1, Sc41, -1, -1, K37, Sc41, k_K37_He4_to_Sc41_reaclib, k_Sc41_to_He4_K37_derived,
        H1, K39, -1, Ca40, -1, -1, K39, Ca40, k_K39_p_to_Ca40_reaclib, k_Ca40_to_p_K39_derived,
        He4, K39, -1, Sc43, -1, -1, K39, Sc43, k_K39_He4_to_Sc43_reaclib, k_Sc43_to_He4_K39_derived,
        H1, Ca40, -1, Sc41, -1, -1, Ca40, Sc41, k_Ca40_p_to_Sc41_reaclib, k_Sc41_to_p_Ca40_derived,
        He4, Ca40, -1, Ti44, -1, -1, Ca40, Ti44, k_Ca40_He4_to_Ti44_reaclib, k_Ti44_to_He4_Ca40_derived,
        H1, Ca42, -1, Sc43, -1, -1, Ca42, Sc43, k_Ca42_p_to_Sc43_reaclib, k_Sc43_to_p_Ca42_derived,
        He4, Ca42, -1, Ti46, -1, -1, Ca42, Ti46, k_Ca42_He4_to_Ti46_reaclib, k_Ti46_to_He4_Ca42_derived,
        H1, Sc43, -1, Ti44, -1, -1, Sc43, Ti44, k_Sc43_p_to_Ti44_reaclib, k_Ti44_to_p_Sc43_derived,
        He4, Sc43, -1, V47, -1, -1, Sc43, V47, k_Sc43_He4_to_V47_reaclib, k_V47_to_He4_Sc43_derived,
        H1, Sc45, -1, Ti46, -1, -1, Sc45, Ti46, k_Sc45_p_to_Ti46_reaclib, k_Ti46_to_p_Sc45_derived,
        He4, Sc45, -1, V49, -1, -1, Sc45, V49, k_Sc45_He4_to_V49_reaclib, k_V49_to_He4_Sc45_derived,
        He4, Ti44, -1, Cr48, -1, -1, Ti44, Cr48, k_Ti44_He4_to_Cr48_reaclib, k_Cr48_to_He4_Ti44_derived,
        H1, Ti46, -1, V47, -1, -1, Ti46, V47, k_Ti46_p_to_V47_reaclib, k_V47_to_p_Ti46_derived,
        He4, Ti46, -1, Cr50, -1, -1, Ti46, Cr50, k_Ti46_He4_to_Cr50_reaclib, k_Cr50_to_He4_Ti46_derived,
        H1, Ti48, -1, V49, -1, -1, Ti48, V49, k_Ti48_p_to_V49_reaclib, k_V49_to_p_Ti48_derived,
        He4, Ti48, -1, Cr52, -1, -1, Ti48, Cr52, k_Ti48_He4_to_Cr52_reaclib, k_Cr52_to_He4_Ti48_derived,
        H1, V47, -1, Cr48, -1, -1, V47, Cr48, k_V47_p_to_Cr48_reaclib, k_Cr48_to_p_V47_derived,
        He4, V47, -1, Mn51, -1, -1, V47, Mn51, k_V47_He4_to_Mn51_reaclib, k_Mn51_to_He4_V47_derived,
        H1, V49, -1, Cr50, -1, -1, V49, Cr50, k_V49_p_to_Cr50_reaclib, k_Cr50_to_p_V49_derived,
        He4, V49, -1, Mn53, -1, -1, V49, Mn53, k_V49_He4_to_Mn53_reaclib, k_Mn53_to_He4_V49_derived,
        H1, V51, -1, Cr52, -1, -1, V51, Cr52, k_V51_p_to_Cr52_reaclib, k_Cr52_to_p_V51_derived,
        He4, V51, -1, Mn55, -1, -1, V51, Mn55, k_V51_He4_to_Mn55_reaclib, k_Mn55_to_He4_V51_derived,
        He4, Cr48, -1, Fe52, -1, -1, Cr48, Fe52, k_Cr48_He4_to_Fe52_reaclib, k_Fe52_to_He4_Cr48_derived,
        H1, Cr50, -1, Mn51, -1, -1, Cr50, Mn51, k_Cr50_p_to_Mn51_reaclib, k_Mn51_to_p_Cr50_derived,
        He4, Cr50, -1, Fe54, -1, -1, Cr50, Fe54, k_Cr50_He4_to_Fe54_reaclib, k_Fe54_to_He4_Cr50_derived,
        H1, Cr52, -1, Mn53, -1, -1, Cr52, Mn53, k_Cr52_p_to_Mn53_reaclib, k_Mn53_to_p_Cr52_derived,
        He4, Cr52, -1, Fe56, -1, -1, Cr52, Fe56, k_Cr52_He4_to_Fe56_reaclib, k_Fe56_to_He4_Cr52_derived,
        H1, Cr54, -1, Mn55, -1, -1, Cr54, Mn55, k_Cr54_p_to_Mn55_reaclib, k_Mn55_to_p_Cr54_derived,
        He4, Cr54, -1, Fe58, -1, -1, Cr54, Fe58, k_Cr54_He4_to_Fe58_reaclib, k_Fe58_to_He4_Cr54_derived,
        H1, Mn51, -1, Fe52, -1, -1, Mn51, Fe52, k_Mn51_p_to_Fe52_reaclib, k_Fe52_to_p_Mn51_derived,
        He4, Mn51, -1, Co55, -1, -1, Mn51, Co55, k_Mn51_He4_to_Co55_reaclib, k_Co55_to_He4_Mn51_derived,
        H1, Mn53, -1, Fe54, -1, -1, Mn53, Fe54, k_Mn53_p_to_Fe54_reaclib, k_Fe54_to_p_Mn53_derived,
        He4, Mn53, -1, Co57, -1, -1, Mn53, Co57, k_Mn53_He4_to_Co57_reaclib, k_Co57_to_He4_Mn53_derived,
        H1, Mn55, -1, Fe56, -1, -1, Mn55, Fe56, k_Mn55_p_to_Fe56_reaclib, k_Fe56_to_p_Mn55_derived,
        He4, Mn55, -1, Co59, -1, -1, Mn55, Co59, k_Mn55_He4_to_Co59_reaclib, k_Co59_to_He4_Mn55_derived,
        He4, Fe52, -1, Ni56, -1, -1, Fe52, Ni56, k_Fe52_He4_to_Ni56_reaclib, k_Ni56_to_He4_Fe52_derived,
        H1, Fe54, -1, Co55, -1, -1, Fe54, Co55, k_Fe54_p_to_Co55_reaclib, k_Co55_to_p_Fe54_derived,
        He4, Fe54, -1, Ni58, -1, -1, Fe54, Ni58, k_Fe54_He4_to_Ni58_reaclib, k_Ni58_to_He4_Fe54_derived,
        H1, Fe56, -1, Co57, -1, -1, Fe56, Co57, k_Fe56_p_to_Co57_reaclib, k_Co57_to_p_Fe56_derived,
        He4, Fe56, -1, Ni60, -1, -1, Fe56, Ni60, k_Fe56_He4_to_Ni60_reaclib, k_Ni60_to_He4_Fe56_derived,
        H1, Fe58, -1, Co59, -1, -1, Fe58, Co59, k_Fe58_p_to_Co59_reaclib, k_Co59_to_p_Fe58_derived,
        N, Co55, -1, Co56, -1, -1, Co55, Co56, k_Co55_n_to_Co56_reaclib, k_Co56_to_n_Co55_derived,
        H1, Co55, -1, Ni56, -1, -1, Co55, Ni56, k_Co55_p_to_Ni56_reaclib, k_Ni56_to_p_Co55_derived,
        He4, Co55, -1, Cu59, -1, -1, Co55, Cu59, k_Co55_He4_to_Cu59_reaclib, k_Cu59_to_He4_Co55_derived,
        N, Co56, -1, Co57, -1, -1, Co56, Co57, k_Co56_n_to_Co57_reaclib, k_Co57_to_n_Co56_derived,
        H1, Co57, -1, Ni58, -1, -1, Co57, Ni58, k_Co57_p_to_Ni58_reaclib, k_Ni58_to_p_Co57_derived,
        H1, Co59, -1, Ni60, -1, -1, Co59, Ni60, k_Co59_p_to_Ni60_reaclib, k_Ni60_to_p_Co59_derived,
        He4, Ni56, -1, Zn60, -1, -1, Ni56, Zn60, k_Ni56_He4_to_Zn60_reaclib, k_Zn60_to_He4_Ni56_derived,
        H1, Ni58, -1, Cu59, -1, -1, Ni58, Cu59, k_Ni58_p_to_Cu59_reaclib, k_Cu59_to_p_Ni58_derived,
        H1, Cu59, -1, Zn60, -1, -1, Cu59, Zn60, k_Cu59_p_to_Zn60_reaclib, k_Zn60_to_p_Cu59_derived,
        N, Mg23, -1, He4, Ne20, -1, Mg23, Ne20, k_Mg23_n_to_He4_Ne20_reaclib, k_Ne20_He4_to_n_Mg23_derived,
        H1, Al27, -1, He4, Mg24, -1, Al27, Mg24, k_Al27_p_to_He4_Mg24_reaclib, k_Mg24_He4_to_p_Al27_derived,
        He4, Al27, -1, H1, Si30, -1, Al27, Si30, k_Al27_He4_to_p_Si30_reaclib, k_Si30_p_to_He4_Al27_derived,
        N, Si27, -1, H1, Al27, -1, Si27, Al27, k_Si27_n_to_p_Al27_reaclib, k_Al27_p_to_n_Si27_derived,
        N, Si27, -1, He4, Mg24, -1, Si27, Mg24, k_Si27_n_to_He4_Mg24_reaclib, k_Mg24_He4_to_n_Si27_derived,
        H1, Si31, -1, N, P31, -1, Si31, P31, k_Si31_p_to_n_P31_reaclib, k_P31_n_to_p_Si31_derived,
        He4, Si31, -1, N, S34, -1, Si31, S34, k_Si31_He4_to_n_S34_reaclib, k_S34_n_to_He4_Si31_derived,
        N, P29, -1, H1, Si29, -1, P29, Si29, k_P29_n_to_p_Si29_reaclib, k_Si29_p_to_n_P29_derived,
        He4, P29, -1, H1, S32, -1, P29, S32, k_P29_He4_to_p_S32_reaclib, k_S32_p_to_He4_P29_derived,
        H1, P31, -1, He4, Si28, -1, P31, Si28, k_P31_p_to_He4_Si28_reaclib, k_Si28_He4_to_p_P31_derived,
        He4, P31, -1, H1, S34, -1, P31, S34, k_P31_He4_to_p_S34_reaclib, k_S34_p_to_He4_P31_derived,
        H1, P33, -1, He4, Si30, -1, P33, Si30, k_P33_p_to_He4_Si30_reaclib, k_Si30_He4_to_p_P33_derived,
        N, S32, -1, He4, Si29, -1, S32, Si29, k_S32_n_to_He4_Si29_reaclib, k_Si29_He4_to_n_S32_derived,
        N, S33, -1, H1, P33, -1, S33, P33, k_S33_n_to_p_P33_reaclib, k_P33_p_to_n_S33_derived,
        N, S33, -1, He4, Si30, -1, S33, Si30, k_S33_n_to_He4_Si30_reaclib, k_Si30_He4_to_n_S33_derived,
        He4, S35, -1, N, Ar38, -1, S35, Ar38, k_S35_He4_to_n_Ar38_reaclib, k_Ar38_n_to_He4_S35_derived,
        N, Cl33, -1, H1, S33, -1, Cl33, S33, k_Cl33_n_to_p_S33_reaclib, k_S33_p_to_n_Cl33_derived,
        He4, Cl33, -1, H1, Ar36, -1, Cl33, Ar36, k_Cl33_He4_to_p_Ar36_reaclib, k_Ar36_p_to_He4_Cl33_derived,
        N, Cl35, -1, H1, S35, -1, Cl35, S35, k_Cl35_n_to_p_S35_reaclib, k_S35_p_to_n_Cl35_derived,
        H1, Cl35, -1, He4, S32, -1, Cl35, S32, k_Cl35_p_to_He4_S32_reaclib, k_S32_He4_to_p_Cl35_derived,
        He4, Cl35, -1, H1, Ar38, -1, Cl35, Ar38, k_Cl35_He4_to_p_Ar38_reaclib, k_Ar38_p_to_He4_Cl35_derived,
        N, Ar36, -1, He4, S33, -1, Ar36, S33, k_Ar36_n_to_He4_S33_reaclib, k_S33_He4_to_n_Ar36_derived,
        N, Ar37, -1, H1, Cl37, -1, Ar37, Cl37, k_Ar37_n_to_p_Cl37_reaclib, k_Cl37_p_to_n_Ar37_derived,
        N, K37, -1, H1, Ar37, -1, K37, Ar37, k_K37_n_to_p_Ar37_reaclib, k_Ar37_p_to_n_K37_derived,
        He4, K37, -1, H1, Ca40, -1, K37, Ca40, k_K37_He4_to_p_Ca40_reaclib, k_Ca40_p_to_He4_K37_derived,
        H1, K39, -1, He4, Ar36, -1, K39, Ar36, k_K39_p_to_He4_Ar36_reaclib, k_Ar36_He4_to_p_K39_derived,
        N, Ca40, -1, He4, Ar37, -1, Ca40, Ar37, k_Ca40_n_to_He4_Ar37_reaclib, k_Ar37_He4_to_n_Ca40_derived,
        H1, Ca42, -1, He4, K39, -1, Ca42, K39, k_Ca42_p_to_He4_K39_reaclib, k_K39_He4_to_p_Ca42_derived,
        He4, Sc41, -1, H1, Ti44, -1, Sc41, Ti44, k_Sc41_He4_to_p_Ti44_reaclib, k_Ti44_p_to_He4_Sc41_derived,
        H1, Sc43, -1, He4, Ca40, -1, Sc43, Ca40, k_Sc43_p_to_He4_Ca40_reaclib, k_Ca40_He4_to_p_Sc43_derived,
        He4, Sc43, -1, H1, Ti46, -1, Sc43, Ti46, k_Sc43_He4_to_p_Ti46_reaclib, k_Ti46_p_to_He4_Sc43_derived,
        H1, Sc45, -1, He4, Ca42, -1, Sc45, Ca42, k_Sc45_p_to_He4_Ca42_reaclib, k_Ca42_He4_to_p_Sc45_derived,
        He4, Sc45, -1, H1, Ti48, -1, Sc45, Ti48, k_Sc45_He4_to_p_Ti48_reaclib, k_Ti48_p_to_He4_Sc45_derived,
        He4, V47, -1, H1, Cr50, -1, V47, Cr50, k_V47_He4_to_p_Cr50_reaclib, k_Cr50_p_to_He4_V47_derived,
        H1, V49, -1, He4, Ti46, -1, V49, Ti46, k_V49_p_to_He4_Ti46_reaclib, k_Ti46_He4_to_p_V49_derived,
        He4, V49, -1, H1, Cr52, -1, V49, Cr52, k_V49_He4_to_p_Cr52_reaclib, k_Cr52_p_to_He4_V49_derived,
        H1, V51, -1, He4, Ti48, -1, V51, Ti48, k_V51_p_to_He4_Ti48_reaclib, k_Ti48_He4_to_p_V51_derived,
        He4, Cr48, -1, H1, Mn51, -1, Cr48, Mn51, k_Cr48_He4_to_p_Mn51_reaclib, k_Mn51_p_to_He4_Cr48_derived,
        H1, Cr54, -1, He4, V51, -1, Cr54, V51, k_Cr54_p_to_He4_V51_reaclib, k_V51_He4_to_p_Cr54_derived,
        He4, Mn51, -1, H1, Fe54, -1, Mn51, Fe54, k_Mn51_He4_to_p_Fe54_reaclib, k_Fe54_p_to_He4_Mn51_derived,
        H1, Mn53, -1, He4, Cr50, -1, Mn53, Cr50, k_Mn53_p_to_He4_Cr50_reaclib, k_Cr50_He4_to_p_Mn53_derived,
        He4, Mn53, -1, H1, Fe56, -1, Mn53, Fe56, k_Mn53_He4_to_p_Fe56_reaclib, k_Fe56_p_to_He4_Mn53_derived,
        H1, Mn55, -1, He4, Cr52, -1, Mn55, Cr52, k_Mn55_p_to_He4_Cr52_reaclib, k_Cr52_He4_to_p_Mn55_derived,
        He4, Fe52, -1, H1, Co55, -1, Fe52, Co55, k_Fe52_He4_to_p_Co55_reaclib, k_Co55_p_to_He4_Fe52_derived,
        H1, Fe58, -1, He4, Mn55, -1, Fe58, Mn55, k_Fe58_p_to_He4_Mn55_reaclib, k_Mn55_He4_to_p_Fe58_derived,
        He4, Co55, -1, H1, Ni58, -1, Co55, Ni58, k_Co55_He4_to_p_Ni58_reaclib, k_Ni58_p_to_He4_Co55_derived,
        N, Co56, -1, H1, Fe56, -1, Co56, Fe56, k_Co56_n_to_p_Fe56_reaclib, k_Fe56_p_to_n_Co56_derived,
        N, Co56, -1, He4, Mn53, -1, Co56, Mn53, k_Co56_n_to_He4_Mn53_reaclib, k_Mn53_He4_to_n_Co56_derived,
        H1, Co57, -1, He4, Fe54, -1, Co57, Fe54, k_Co57_p_to_He4_Fe54_reaclib, k_Fe54_He4_to_p_Co57_derived,
        He4, Co57, -1, H1, Ni60, -1, Co57, Ni60, k_Co57_He4_to_p_Ni60_reaclib, k_Ni60_p_to_He4_Co57_derived,
        H1, Co59, -1, He4, Fe56, -1, Co59, Fe56, k_Co59_p_to_He4_Fe56_reaclib, k_Fe56_He4_to_p_Co59_derived,
        N, Ni56, -1, H1, Co56, -1, Ni56, Co56, k_Ni56_n_to_p_Co56_reaclib, k_Co56_p_to_n_Ni56_derived,
        N, Cu59, -1, He4, Co56, -1, Cu59, Co56, k_Cu59_n_to_He4_Co56_reaclib, k_Co56_He4_to_n_Cu59_derived,
        H1, Cu59, -1, He4, Ni56, -1, Cu59, Ni56, k_Cu59_p_to_He4_Ni56_reaclib, k_Ni56_He4_to_p_Cu59_derived,
        He4, He4, He4, C12, -1, -1, C12, -1, k_He4_He4_He4_to_C12_reaclib, k_C12_to_3He4_derived,
        H1, Cl37, -1, He4, S34, -1, Cl37, S34, k_Cl37_p_to_He4_S34_derived, k_S34_He4_to_p_Cl37_reaclib,
        N, Ar37, -1, He4, S34, -1, Ar37, S34, k_Ar37_n_to_He4_S34_derived, k_S34_He4_to_n_Ar37_reaclib,
        H1, V47, -1, He4, Ti44, -1, V47, Ti44, k_V47_p_to_He4_Ti44_derived, k_Ti44_He4_to_p_V47_reaclib,
        Ne20, He4, -1, Mg24, -1, -1, Ne20, Mg24, k_Ne20_He4_to_Mg24_approx, k_Mg24_to_Ne20_He4_approx
    };
}
#endif

#endif
