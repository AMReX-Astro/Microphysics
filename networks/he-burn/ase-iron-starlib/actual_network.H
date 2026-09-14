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


#define STARLIB

namespace starlib {

    constexpr std::uint8_t NumStarLibRates = 65;
    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumStarLibRates> prand{};
}
inline
void actual_network_init() {


    if (network_rp::starlib_seed > 0) {
        // generate Gaussian random numbers
        std::mt19937 generator(network_rp::starlib_seed);

        // normal distribution centered on 0 with a width of 1
        std::normal_distribution rn(0.0, 1.0);

        // generate the random numbers -- every process will
        // get the same numbers here, since we are using the
        // same seed.
        for (int n = 1; n <= starlib::NumStarLibRates; ++n) {
            starlib::prand(n) = rn(generator);
        }
    }
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
        else if constexpr (spec == N13) {
            return 94.10522604799917_rt;
        }
        else if constexpr (spec == O16) {
            return 127.6193154119992_rt;
        }
        else if constexpr (spec == Ne20) {
            return 160.64482384000075_rt;
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
        else if constexpr (spec == Mn51) {
            return 440.321747199996_rt;
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
        else if constexpr (spec == N13) {
            return 2.1596537049448796e-23_rt;
        }
        else if constexpr (spec == O16) {
            return 2.6560180592333686e-23_rt;
        }
        else if constexpr (spec == Ne20) {
            return 3.3198227947612416e-23_rt;
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
        else if constexpr (spec == Mn51) {
            return 8.460149108098279e-23_rt;
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
    std::uint8_t
    {
        k_Co56_to_Fe56_weaktab = 1,
        k_Co56_to_Ni56_weaktab = 2,
        k_Fe56_to_Co56_weaktab = 3,
        k_Ni56_to_Co56_weaktab = 4,
        k_n_to_p_weaktab = 5,
        k_p_to_n_weaktab = 6,
        k_C12_C12_to_p_Na23_starlib = 7,
        k_C12_C12_to_He4_Ne20_starlib = 8,
        k_C12_O16_to_p_Al27_starlib = 9,
        k_C12_O16_to_He4_Mg24_starlib = 10,
        k_O16_O16_to_p_P31_starlib = 11,
        k_O16_O16_to_He4_Si28_starlib = 12,
        k_He4_He4_He4_to_C12_starlib = 13,
        k_p_C12_to_N13_starlib = 14,
        k_He4_C12_to_O16_starlib = 15,
        k_He4_N13_to_p_O16_starlib = 16,
        k_He4_O16_to_Ne20_starlib = 17,
        k_He4_Ne20_to_Mg24_starlib = 18,
        k_p_Na23_to_He4_Ne20_starlib = 19,
        k_p_Na23_to_Mg24_starlib = 20,
        k_He4_Mg24_to_Si28_starlib = 21,
        k_p_Al27_to_He4_Mg24_starlib = 22,
        k_p_Al27_to_Si28_starlib = 23,
        k_He4_Si28_to_S32_starlib = 24,
        k_p_P31_to_He4_Si28_starlib = 25,
        k_p_P31_to_S32_starlib = 26,
        k_He4_Cr48_to_p_Mn51_starlib = 27,
        k_He4_Cr48_to_Fe52_starlib = 28,
        k_p_Mn51_to_Fe52_starlib = 29,
        k_He4_Mn51_to_Co55_starlib = 30,
        k_He4_Fe52_to_p_Co55_starlib = 31,
        k_He4_Fe52_to_Ni56_starlib = 32,
        k_p_Co55_to_Ni56_starlib = 33,
        k_He4_Co55_to_Cu59_starlib = 34,
        k_He4_Ni56_to_Zn60_starlib = 35,
        k_p_Cu59_to_He4_Ni56_starlib = 36,
        k_p_Cu59_to_Zn60_starlib = 37,
        k_C12_C12_to_n_Mg23_starlib = 38,
        k_O16_O16_to_n_S31_starlib = 39,
        k_C12_O16_to_n_Si27_starlib = 40,
        k_He4_Mn51_to_p_Fe54_starlib = 41,
        k_p_Fe54_to_Co55_starlib = 42,
        k_p_Fe56_to_Co57_starlib = 43,
        k_He4_Fe54_to_Ni58_starlib = 44,
        k_n_Co55_to_Co56_starlib = 45,
        k_n_Co56_to_p_Fe56_starlib = 46,
        k_n_Co56_to_Co57_starlib = 47,
        k_p_Co57_to_He4_Fe54_starlib = 48,
        k_p_Co57_to_Ni58_starlib = 49,
        k_He4_Co55_to_p_Ni58_starlib = 50,
        k_n_Ni56_to_p_Co56_starlib = 51,
        k_p_Ni58_to_Cu59_starlib = 52,
        k_n_Cu59_to_He4_Co56_starlib = 53,
        k_He4_S32_to_Ar36_starlib = 54,
        k_p_Cl35_to_Ar36_starlib = 55,
        k_p_Cl35_to_He4_S32_starlib = 56,
        k_He4_Ar36_to_Ca40_starlib = 57,
        k_p_K39_to_Ca40_starlib = 58,
        k_p_K39_to_He4_Ar36_starlib = 59,
        k_He4_Ca40_to_Ti44_starlib = 60,
        k_p_Sc43_to_Ti44_starlib = 61,
        k_p_Sc43_to_He4_Ca40_starlib = 62,
        k_He4_Ti44_to_Cr48_starlib = 63,
        k_p_V47_to_Cr48_starlib = 64,
        k_p_V47_to_He4_Ti44_starlib = 65,
        k_n_Fe52_to_Fe53_starlib = 66,
        k_n_Fe53_to_Fe54_starlib = 67,
        k_n_Fe54_to_Fe55_starlib = 68,
        k_n_Fe55_to_Fe56_starlib = 69,
        k_n_Ni56_to_Ni57_starlib = 70,
        k_n_Ni57_to_Ni58_starlib = 71,
        k_S32_He4_to_Ar36_approx = 72,
        k_Ar36_to_S32_He4_approx = 73,
        k_Ar36_He4_to_Ca40_approx = 74,
        k_Ca40_to_Ar36_He4_approx = 75,
        k_Ca40_He4_to_Ti44_approx = 76,
        k_Ti44_to_Ca40_He4_approx = 77,
        k_Ti44_He4_to_Cr48_approx = 78,
        k_Cr48_to_Ti44_He4_approx = 79,
        k_Fe52_n_n_to_Fe54_approx = 80,
        k_Fe54_to_Fe52_n_n_approx = 81,
        k_Fe54_n_n_to_Fe56_approx = 82,
        k_Fe56_to_Fe54_n_n_approx = 83,
        k_Ni56_n_n_to_Ni58_approx = 84,
        k_Ni58_to_Ni56_n_n_approx = 85,
        k_C12_C12_to_Mg24_modified = 86,
        k_O16_O16_to_S32_modified = 87,
        k_C12_O16_to_Si28_modified = 88,
        k_C12_to_He4_He4_He4_derived = 89,
        k_He4_Ne20_to_C12_C12_derived = 90,
        k_p_Na23_to_C12_C12_derived = 91,
        k_He4_Mg24_to_C12_O16_derived = 92,
        k_p_Al27_to_C12_O16_derived = 93,
        k_He4_Si28_to_O16_O16_derived = 94,
        k_p_P31_to_O16_O16_derived = 95,
        k_p_O16_to_He4_N13_derived = 96,
        k_N13_to_p_C12_derived = 97,
        k_O16_to_He4_C12_derived = 98,
        k_Ne20_to_He4_O16_derived = 99,
        k_He4_Ne20_to_p_Na23_derived = 100,
        k_Mg24_to_p_Na23_derived = 101,
        k_Mg24_to_He4_Ne20_derived = 102,
        k_He4_Mg24_to_p_Al27_derived = 103,
        k_Si28_to_p_Al27_derived = 104,
        k_Si28_to_He4_Mg24_derived = 105,
        k_He4_Si28_to_p_P31_derived = 106,
        k_S32_to_p_P31_derived = 107,
        k_S32_to_He4_Si28_derived = 108,
        k_p_Mn51_to_He4_Cr48_derived = 109,
        k_p_Co55_to_He4_Fe52_derived = 110,
        k_Fe52_to_p_Mn51_derived = 111,
        k_Fe52_to_He4_Cr48_derived = 112,
        k_Co55_to_He4_Mn51_derived = 113,
        k_Ni56_to_p_Co55_derived = 114,
        k_Ni56_to_He4_Fe52_derived = 115,
        k_He4_Ni56_to_p_Cu59_derived = 116,
        k_Cu59_to_He4_Co55_derived = 117,
        k_Zn60_to_p_Cu59_derived = 118,
        k_Zn60_to_He4_Ni56_derived = 119,
        k_Mg24_to_C12_C12_derived = 120,
        k_S32_to_O16_O16_derived = 121,
        k_Si28_to_C12_O16_derived = 122,
        k_p_Fe54_to_He4_Mn51_derived = 123,
        k_p_Fe56_to_n_Co56_derived = 124,
        k_He4_Fe54_to_p_Co57_derived = 125,
        k_p_Ni58_to_He4_Co55_derived = 126,
        k_Co55_to_p_Fe54_derived = 127,
        k_Co56_to_n_Co55_derived = 128,
        k_Co57_to_n_Co56_derived = 129,
        k_Co57_to_p_Fe56_derived = 130,
        k_p_Co56_to_n_Ni56_derived = 131,
        k_Ni58_to_p_Co57_derived = 132,
        k_Ni58_to_He4_Fe54_derived = 133,
        k_He4_Co56_to_n_Cu59_derived = 134,
        k_Cu59_to_p_Ni58_derived = 135,
        k_He4_S32_to_p_Cl35_derived = 136,
        k_Ar36_to_He4_S32_derived = 137,
        k_Ar36_to_p_Cl35_derived = 138,
        k_He4_Ar36_to_p_K39_derived = 139,
        k_Ca40_to_He4_Ar36_derived = 140,
        k_Ca40_to_p_K39_derived = 141,
        k_He4_Ca40_to_p_Sc43_derived = 142,
        k_Ti44_to_He4_Ca40_derived = 143,
        k_Ti44_to_p_Sc43_derived = 144,
        k_He4_Ti44_to_p_V47_derived = 145,
        k_Cr48_to_He4_Ti44_derived = 146,
        k_Cr48_to_p_V47_derived = 147,
        k_Fe54_to_n_Fe53_derived = 148,
        k_Fe53_to_n_Fe52_derived = 149,
        k_Fe56_to_n_Fe55_derived = 150,
        k_Fe55_to_n_Fe54_derived = 151,
        k_Ni58_to_n_Ni57_derived = 152,
        k_Ni57_to_n_Ni56_derived = 153,
        NumRates = k_Ni57_to_n_Ni56_derived
    };

    enum ScreenPairs :
    std::uint8_t
    {
        k_p_C12 = 1,
        k_He4_O16 = 2,
        k_p_Co57 = 3,
        k_p_K39 = 4,
        k_He4_Ca40 = 5,
        k_He4_Ar36 = 6,
        k_p_Sc43 = 7,
        k_p_Fe56 = 8,
        k_He4_C12 = 9,
        k_He4_Ti44 = 10,
        k_He4_Mg24 = 11,
        k_C12_O16 = 12,
        k_p_V47 = 13,
        k_He4_Cr48 = 14,
        k_C12_C12 = 15,
        k_p_Co55 = 16,
        k_He4_Ne20 = 17,
        k_p_Al27 = 18,
        k_p_Mn51 = 19,
        k_p_Na23 = 20,
        k_p_Co56 = 21,
        k_p_Fe54 = 22,
        k_He4_Fe52 = 23,
        k_p_Cu59 = 24,
        k_He4_S32 = 25,
        k_He4_Mn51 = 26,
        k_He4_Co55 = 27,
        k_He4_Ni56 = 28,
        k_p_P31 = 29,
        k_He4_He4 = 30,
        k_p_Ni58 = 31,
        k_p_O16 = 32,
        k_He4_Co56 = 33,
        k_He4_Si28 = 34,
        k_He4_Fe54 = 35,
        k_p_Cl35 = 36,
        k_O16_O16 = 37,
        k_He4_Be8 = 38,
        k_He4_N13 = 39,
        NumScreenPairs = k_He4_N13
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "Co56_to_Fe56_weaktab",  // 1,
        "Co56_to_Ni56_weaktab",  // 2,
        "Fe56_to_Co56_weaktab",  // 3,
        "Ni56_to_Co56_weaktab",  // 4,
        "n_to_p_weaktab",  // 5,
        "p_to_n_weaktab",  // 6,
        "C12_C12_to_p_Na23_starlib",  // 7,
        "C12_C12_to_He4_Ne20_starlib",  // 8,
        "C12_O16_to_p_Al27_starlib",  // 9,
        "C12_O16_to_He4_Mg24_starlib",  // 10,
        "O16_O16_to_p_P31_starlib",  // 11,
        "O16_O16_to_He4_Si28_starlib",  // 12,
        "He4_He4_He4_to_C12_starlib",  // 13,
        "p_C12_to_N13_starlib",  // 14,
        "He4_C12_to_O16_starlib",  // 15,
        "He4_N13_to_p_O16_starlib",  // 16,
        "He4_O16_to_Ne20_starlib",  // 17,
        "He4_Ne20_to_Mg24_starlib",  // 18,
        "p_Na23_to_He4_Ne20_starlib",  // 19,
        "p_Na23_to_Mg24_starlib",  // 20,
        "He4_Mg24_to_Si28_starlib",  // 21,
        "p_Al27_to_He4_Mg24_starlib",  // 22,
        "p_Al27_to_Si28_starlib",  // 23,
        "He4_Si28_to_S32_starlib",  // 24,
        "p_P31_to_He4_Si28_starlib",  // 25,
        "p_P31_to_S32_starlib",  // 26,
        "He4_Cr48_to_p_Mn51_starlib",  // 27,
        "He4_Cr48_to_Fe52_starlib",  // 28,
        "p_Mn51_to_Fe52_starlib",  // 29,
        "He4_Mn51_to_Co55_starlib",  // 30,
        "He4_Fe52_to_p_Co55_starlib",  // 31,
        "He4_Fe52_to_Ni56_starlib",  // 32,
        "p_Co55_to_Ni56_starlib",  // 33,
        "He4_Co55_to_Cu59_starlib",  // 34,
        "He4_Ni56_to_Zn60_starlib",  // 35,
        "p_Cu59_to_He4_Ni56_starlib",  // 36,
        "p_Cu59_to_Zn60_starlib",  // 37,
        "C12_C12_to_n_Mg23_starlib",  // 38,
        "O16_O16_to_n_S31_starlib",  // 39,
        "C12_O16_to_n_Si27_starlib",  // 40,
        "He4_Mn51_to_p_Fe54_starlib",  // 41,
        "p_Fe54_to_Co55_starlib",  // 42,
        "p_Fe56_to_Co57_starlib",  // 43,
        "He4_Fe54_to_Ni58_starlib",  // 44,
        "n_Co55_to_Co56_starlib",  // 45,
        "n_Co56_to_p_Fe56_starlib",  // 46,
        "n_Co56_to_Co57_starlib",  // 47,
        "p_Co57_to_He4_Fe54_starlib",  // 48,
        "p_Co57_to_Ni58_starlib",  // 49,
        "He4_Co55_to_p_Ni58_starlib",  // 50,
        "n_Ni56_to_p_Co56_starlib",  // 51,
        "p_Ni58_to_Cu59_starlib",  // 52,
        "n_Cu59_to_He4_Co56_starlib",  // 53,
        "He4_S32_to_Ar36_starlib",  // 54,
        "p_Cl35_to_Ar36_starlib",  // 55,
        "p_Cl35_to_He4_S32_starlib",  // 56,
        "He4_Ar36_to_Ca40_starlib",  // 57,
        "p_K39_to_Ca40_starlib",  // 58,
        "p_K39_to_He4_Ar36_starlib",  // 59,
        "He4_Ca40_to_Ti44_starlib",  // 60,
        "p_Sc43_to_Ti44_starlib",  // 61,
        "p_Sc43_to_He4_Ca40_starlib",  // 62,
        "He4_Ti44_to_Cr48_starlib",  // 63,
        "p_V47_to_Cr48_starlib",  // 64,
        "p_V47_to_He4_Ti44_starlib",  // 65,
        "n_Fe52_to_Fe53_starlib",  // 66,
        "n_Fe53_to_Fe54_starlib",  // 67,
        "n_Fe54_to_Fe55_starlib",  // 68,
        "n_Fe55_to_Fe56_starlib",  // 69,
        "n_Ni56_to_Ni57_starlib",  // 70,
        "n_Ni57_to_Ni58_starlib",  // 71,
        "S32_He4_to_Ar36_approx",  // 72,
        "Ar36_to_S32_He4_approx",  // 73,
        "Ar36_He4_to_Ca40_approx",  // 74,
        "Ca40_to_Ar36_He4_approx",  // 75,
        "Ca40_He4_to_Ti44_approx",  // 76,
        "Ti44_to_Ca40_He4_approx",  // 77,
        "Ti44_He4_to_Cr48_approx",  // 78,
        "Cr48_to_Ti44_He4_approx",  // 79,
        "Fe52_n_n_to_Fe54_approx",  // 80,
        "Fe54_to_Fe52_n_n_approx",  // 81,
        "Fe54_n_n_to_Fe56_approx",  // 82,
        "Fe56_to_Fe54_n_n_approx",  // 83,
        "Ni56_n_n_to_Ni58_approx",  // 84,
        "Ni58_to_Ni56_n_n_approx",  // 85,
        "C12_C12_to_Mg24_modified",  // 86,
        "O16_O16_to_S32_modified",  // 87,
        "C12_O16_to_Si28_modified",  // 88,
        "C12_to_He4_He4_He4_derived",  // 89,
        "He4_Ne20_to_C12_C12_derived",  // 90,
        "p_Na23_to_C12_C12_derived",  // 91,
        "He4_Mg24_to_C12_O16_derived",  // 92,
        "p_Al27_to_C12_O16_derived",  // 93,
        "He4_Si28_to_O16_O16_derived",  // 94,
        "p_P31_to_O16_O16_derived",  // 95,
        "p_O16_to_He4_N13_derived",  // 96,
        "N13_to_p_C12_derived",  // 97,
        "O16_to_He4_C12_derived",  // 98,
        "Ne20_to_He4_O16_derived",  // 99,
        "He4_Ne20_to_p_Na23_derived",  // 100,
        "Mg24_to_p_Na23_derived",  // 101,
        "Mg24_to_He4_Ne20_derived",  // 102,
        "He4_Mg24_to_p_Al27_derived",  // 103,
        "Si28_to_p_Al27_derived",  // 104,
        "Si28_to_He4_Mg24_derived",  // 105,
        "He4_Si28_to_p_P31_derived",  // 106,
        "S32_to_p_P31_derived",  // 107,
        "S32_to_He4_Si28_derived",  // 108,
        "p_Mn51_to_He4_Cr48_derived",  // 109,
        "p_Co55_to_He4_Fe52_derived",  // 110,
        "Fe52_to_p_Mn51_derived",  // 111,
        "Fe52_to_He4_Cr48_derived",  // 112,
        "Co55_to_He4_Mn51_derived",  // 113,
        "Ni56_to_p_Co55_derived",  // 114,
        "Ni56_to_He4_Fe52_derived",  // 115,
        "He4_Ni56_to_p_Cu59_derived",  // 116,
        "Cu59_to_He4_Co55_derived",  // 117,
        "Zn60_to_p_Cu59_derived",  // 118,
        "Zn60_to_He4_Ni56_derived",  // 119,
        "Mg24_to_C12_C12_derived",  // 120,
        "S32_to_O16_O16_derived",  // 121,
        "Si28_to_C12_O16_derived",  // 122,
        "p_Fe54_to_He4_Mn51_derived",  // 123,
        "p_Fe56_to_n_Co56_derived",  // 124,
        "He4_Fe54_to_p_Co57_derived",  // 125,
        "p_Ni58_to_He4_Co55_derived",  // 126,
        "Co55_to_p_Fe54_derived",  // 127,
        "Co56_to_n_Co55_derived",  // 128,
        "Co57_to_n_Co56_derived",  // 129,
        "Co57_to_p_Fe56_derived",  // 130,
        "p_Co56_to_n_Ni56_derived",  // 131,
        "Ni58_to_p_Co57_derived",  // 132,
        "Ni58_to_He4_Fe54_derived",  // 133,
        "He4_Co56_to_n_Cu59_derived",  // 134,
        "Cu59_to_p_Ni58_derived",  // 135,
        "He4_S32_to_p_Cl35_derived",  // 136,
        "Ar36_to_He4_S32_derived",  // 137,
        "Ar36_to_p_Cl35_derived",  // 138,
        "He4_Ar36_to_p_K39_derived",  // 139,
        "Ca40_to_He4_Ar36_derived",  // 140,
        "Ca40_to_p_K39_derived",  // 141,
        "He4_Ca40_to_p_Sc43_derived",  // 142,
        "Ti44_to_He4_Ca40_derived",  // 143,
        "Ti44_to_p_Sc43_derived",  // 144,
        "He4_Ti44_to_p_V47_derived",  // 145,
        "Cr48_to_He4_Ti44_derived",  // 146,
        "Cr48_to_p_V47_derived",  // 147,
        "Fe54_to_n_Fe53_derived",  // 148,
        "Fe53_to_n_Fe52_derived",  // 149,
        "Fe56_to_n_Fe55_derived",  // 150,
        "Fe55_to_n_Fe54_derived",  // 151,
        "Ni58_to_n_Ni57_derived",  // 152,
        "Ni57_to_n_Ni56_derived"  // 153,
    };

    template<std::uint8_t rate>
    AMREX_GPU_DEVICE AMREX_INLINE
    constexpr amrex::Real get_p_random() {
        if constexpr (rate == k_C12_C12_to_p_Na23_starlib) {
            return starlib::prand(1);
        }

        if constexpr (rate == k_C12_C12_to_He4_Ne20_starlib) {
            return starlib::prand(2);
        }

        if constexpr (rate == k_C12_O16_to_p_Al27_starlib) {
            return starlib::prand(3);
        }

        if constexpr (rate == k_C12_O16_to_He4_Mg24_starlib) {
            return starlib::prand(4);
        }

        if constexpr (rate == k_O16_O16_to_p_P31_starlib) {
            return starlib::prand(5);
        }

        if constexpr (rate == k_O16_O16_to_He4_Si28_starlib) {
            return starlib::prand(6);
        }

        if constexpr (rate == k_He4_He4_He4_to_C12_starlib) {
            return starlib::prand(7);
        }

        if constexpr (rate == k_p_C12_to_N13_starlib) {
            return starlib::prand(8);
        }

        if constexpr (rate == k_He4_C12_to_O16_starlib) {
            return starlib::prand(9);
        }

        if constexpr (rate == k_He4_N13_to_p_O16_starlib) {
            return starlib::prand(10);
        }

        if constexpr (rate == k_He4_O16_to_Ne20_starlib) {
            return starlib::prand(11);
        }

        if constexpr (rate == k_He4_Ne20_to_Mg24_starlib) {
            return starlib::prand(12);
        }

        if constexpr (rate == k_p_Na23_to_He4_Ne20_starlib) {
            return starlib::prand(13);
        }

        if constexpr (rate == k_p_Na23_to_Mg24_starlib) {
            return starlib::prand(14);
        }

        if constexpr (rate == k_He4_Mg24_to_Si28_starlib) {
            return starlib::prand(15);
        }

        if constexpr (rate == k_p_Al27_to_He4_Mg24_starlib) {
            return starlib::prand(16);
        }

        if constexpr (rate == k_p_Al27_to_Si28_starlib) {
            return starlib::prand(17);
        }

        if constexpr (rate == k_He4_Si28_to_S32_starlib) {
            return starlib::prand(18);
        }

        if constexpr (rate == k_p_P31_to_He4_Si28_starlib) {
            return starlib::prand(19);
        }

        if constexpr (rate == k_p_P31_to_S32_starlib) {
            return starlib::prand(20);
        }

        if constexpr (rate == k_He4_Cr48_to_p_Mn51_starlib) {
            return starlib::prand(21);
        }

        if constexpr (rate == k_He4_Cr48_to_Fe52_starlib) {
            return starlib::prand(22);
        }

        if constexpr (rate == k_p_Mn51_to_Fe52_starlib) {
            return starlib::prand(23);
        }

        if constexpr (rate == k_He4_Mn51_to_Co55_starlib) {
            return starlib::prand(24);
        }

        if constexpr (rate == k_He4_Fe52_to_p_Co55_starlib) {
            return starlib::prand(25);
        }

        if constexpr (rate == k_He4_Fe52_to_Ni56_starlib) {
            return starlib::prand(26);
        }

        if constexpr (rate == k_p_Co55_to_Ni56_starlib) {
            return starlib::prand(27);
        }

        if constexpr (rate == k_He4_Co55_to_Cu59_starlib) {
            return starlib::prand(28);
        }

        if constexpr (rate == k_He4_Ni56_to_Zn60_starlib) {
            return starlib::prand(29);
        }

        if constexpr (rate == k_p_Cu59_to_He4_Ni56_starlib) {
            return starlib::prand(30);
        }

        if constexpr (rate == k_p_Cu59_to_Zn60_starlib) {
            return starlib::prand(31);
        }

        if constexpr (rate == k_C12_C12_to_n_Mg23_starlib) {
            return starlib::prand(32);
        }

        if constexpr (rate == k_O16_O16_to_n_S31_starlib) {
            return starlib::prand(33);
        }

        if constexpr (rate == k_C12_O16_to_n_Si27_starlib) {
            return starlib::prand(34);
        }

        if constexpr (rate == k_He4_Mn51_to_p_Fe54_starlib) {
            return starlib::prand(35);
        }

        if constexpr (rate == k_p_Fe54_to_Co55_starlib) {
            return starlib::prand(36);
        }

        if constexpr (rate == k_p_Fe56_to_Co57_starlib) {
            return starlib::prand(37);
        }

        if constexpr (rate == k_He4_Fe54_to_Ni58_starlib) {
            return starlib::prand(38);
        }

        if constexpr (rate == k_n_Co55_to_Co56_starlib) {
            return starlib::prand(39);
        }

        if constexpr (rate == k_n_Co56_to_p_Fe56_starlib) {
            return starlib::prand(40);
        }

        if constexpr (rate == k_n_Co56_to_Co57_starlib) {
            return starlib::prand(41);
        }

        if constexpr (rate == k_p_Co57_to_He4_Fe54_starlib) {
            return starlib::prand(42);
        }

        if constexpr (rate == k_p_Co57_to_Ni58_starlib) {
            return starlib::prand(43);
        }

        if constexpr (rate == k_He4_Co55_to_p_Ni58_starlib) {
            return starlib::prand(44);
        }

        if constexpr (rate == k_n_Ni56_to_p_Co56_starlib) {
            return starlib::prand(45);
        }

        if constexpr (rate == k_p_Ni58_to_Cu59_starlib) {
            return starlib::prand(46);
        }

        if constexpr (rate == k_n_Cu59_to_He4_Co56_starlib) {
            return starlib::prand(47);
        }

        if constexpr (rate == k_He4_S32_to_Ar36_starlib) {
            return starlib::prand(48);
        }

        if constexpr (rate == k_p_Cl35_to_Ar36_starlib) {
            return starlib::prand(49);
        }

        if constexpr (rate == k_p_Cl35_to_He4_S32_starlib) {
            return starlib::prand(50);
        }

        if constexpr (rate == k_He4_Ar36_to_Ca40_starlib) {
            return starlib::prand(51);
        }

        if constexpr (rate == k_p_K39_to_Ca40_starlib) {
            return starlib::prand(52);
        }

        if constexpr (rate == k_p_K39_to_He4_Ar36_starlib) {
            return starlib::prand(53);
        }

        if constexpr (rate == k_He4_Ca40_to_Ti44_starlib) {
            return starlib::prand(54);
        }

        if constexpr (rate == k_p_Sc43_to_Ti44_starlib) {
            return starlib::prand(55);
        }

        if constexpr (rate == k_p_Sc43_to_He4_Ca40_starlib) {
            return starlib::prand(56);
        }

        if constexpr (rate == k_He4_Ti44_to_Cr48_starlib) {
            return starlib::prand(57);
        }

        if constexpr (rate == k_p_V47_to_Cr48_starlib) {
            return starlib::prand(58);
        }

        if constexpr (rate == k_p_V47_to_He4_Ti44_starlib) {
            return starlib::prand(59);
        }

        if constexpr (rate == k_n_Fe52_to_Fe53_starlib) {
            return starlib::prand(60);
        }

        if constexpr (rate == k_n_Fe53_to_Fe54_starlib) {
            return starlib::prand(61);
        }

        if constexpr (rate == k_n_Fe54_to_Fe55_starlib) {
            return starlib::prand(62);
        }

        if constexpr (rate == k_n_Fe55_to_Fe56_starlib) {
            return starlib::prand(63);
        }

        if constexpr (rate == k_n_Ni56_to_Ni57_starlib) {
            return starlib::prand(64);
        }

        if constexpr (rate == k_n_Ni57_to_Ni58_starlib) {
            return starlib::prand(65);
        }

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
    constexpr int NumNSERatePairs = 42;

    inline AMREX_GPU_MANAGED amrex::Array2D<short, 1, NumNSERatePairs, 1, 10, amrex::Order::C> rate_pair_data {
        He4, He4, He4, C12, -1, -1, C12, -1, k_He4_He4_He4_to_C12_starlib, k_C12_to_He4_He4_He4_derived,
        H1, C12, -1, N13, -1, -1, C12, N13, k_p_C12_to_N13_starlib, k_N13_to_p_C12_derived,
        He4, C12, -1, O16, -1, -1, C12, O16, k_He4_C12_to_O16_starlib, k_O16_to_He4_C12_derived,
        He4, N13, -1, H1, O16, -1, N13, O16, k_He4_N13_to_p_O16_starlib, k_p_O16_to_He4_N13_derived,
        He4, O16, -1, Ne20, -1, -1, O16, Ne20, k_He4_O16_to_Ne20_starlib, k_Ne20_to_He4_O16_derived,
        He4, Ne20, -1, Mg24, -1, -1, Ne20, Mg24, k_He4_Ne20_to_Mg24_starlib, k_Mg24_to_He4_Ne20_derived,
        H1, Na23, -1, He4, Ne20, -1, Na23, Ne20, k_p_Na23_to_He4_Ne20_starlib, k_He4_Ne20_to_p_Na23_derived,
        H1, Na23, -1, Mg24, -1, -1, Na23, Mg24, k_p_Na23_to_Mg24_starlib, k_Mg24_to_p_Na23_derived,
        He4, Mg24, -1, Si28, -1, -1, Mg24, Si28, k_He4_Mg24_to_Si28_starlib, k_Si28_to_He4_Mg24_derived,
        H1, Al27, -1, He4, Mg24, -1, Al27, Mg24, k_p_Al27_to_He4_Mg24_starlib, k_He4_Mg24_to_p_Al27_derived,
        H1, Al27, -1, Si28, -1, -1, Al27, Si28, k_p_Al27_to_Si28_starlib, k_Si28_to_p_Al27_derived,
        He4, Si28, -1, S32, -1, -1, Si28, S32, k_He4_Si28_to_S32_starlib, k_S32_to_He4_Si28_derived,
        H1, P31, -1, He4, Si28, -1, P31, Si28, k_p_P31_to_He4_Si28_starlib, k_He4_Si28_to_p_P31_derived,
        H1, P31, -1, S32, -1, -1, P31, S32, k_p_P31_to_S32_starlib, k_S32_to_p_P31_derived,
        He4, Cr48, -1, H1, Mn51, -1, Cr48, Mn51, k_He4_Cr48_to_p_Mn51_starlib, k_p_Mn51_to_He4_Cr48_derived,
        He4, Cr48, -1, Fe52, -1, -1, Cr48, Fe52, k_He4_Cr48_to_Fe52_starlib, k_Fe52_to_He4_Cr48_derived,
        H1, Mn51, -1, Fe52, -1, -1, Mn51, Fe52, k_p_Mn51_to_Fe52_starlib, k_Fe52_to_p_Mn51_derived,
        He4, Mn51, -1, Co55, -1, -1, Mn51, Co55, k_He4_Mn51_to_Co55_starlib, k_Co55_to_He4_Mn51_derived,
        He4, Fe52, -1, H1, Co55, -1, Fe52, Co55, k_He4_Fe52_to_p_Co55_starlib, k_p_Co55_to_He4_Fe52_derived,
        He4, Fe52, -1, Ni56, -1, -1, Fe52, Ni56, k_He4_Fe52_to_Ni56_starlib, k_Ni56_to_He4_Fe52_derived,
        H1, Co55, -1, Ni56, -1, -1, Co55, Ni56, k_p_Co55_to_Ni56_starlib, k_Ni56_to_p_Co55_derived,
        He4, Co55, -1, Cu59, -1, -1, Co55, Cu59, k_He4_Co55_to_Cu59_starlib, k_Cu59_to_He4_Co55_derived,
        He4, Ni56, -1, Zn60, -1, -1, Ni56, Zn60, k_He4_Ni56_to_Zn60_starlib, k_Zn60_to_He4_Ni56_derived,
        H1, Cu59, -1, He4, Ni56, -1, Cu59, Ni56, k_p_Cu59_to_He4_Ni56_starlib, k_He4_Ni56_to_p_Cu59_derived,
        H1, Cu59, -1, Zn60, -1, -1, Cu59, Zn60, k_p_Cu59_to_Zn60_starlib, k_Zn60_to_p_Cu59_derived,
        He4, Mn51, -1, H1, Fe54, -1, Mn51, Fe54, k_He4_Mn51_to_p_Fe54_starlib, k_p_Fe54_to_He4_Mn51_derived,
        H1, Fe54, -1, Co55, -1, -1, Fe54, Co55, k_p_Fe54_to_Co55_starlib, k_Co55_to_p_Fe54_derived,
        H1, Fe56, -1, Co57, -1, -1, Fe56, Co57, k_p_Fe56_to_Co57_starlib, k_Co57_to_p_Fe56_derived,
        He4, Fe54, -1, Ni58, -1, -1, Fe54, Ni58, k_He4_Fe54_to_Ni58_starlib, k_Ni58_to_He4_Fe54_derived,
        N, Co55, -1, Co56, -1, -1, Co55, Co56, k_n_Co55_to_Co56_starlib, k_Co56_to_n_Co55_derived,
        N, Co56, -1, H1, Fe56, -1, Co56, Fe56, k_n_Co56_to_p_Fe56_starlib, k_p_Fe56_to_n_Co56_derived,
        N, Co56, -1, Co57, -1, -1, Co56, Co57, k_n_Co56_to_Co57_starlib, k_Co57_to_n_Co56_derived,
        H1, Co57, -1, He4, Fe54, -1, Co57, Fe54, k_p_Co57_to_He4_Fe54_starlib, k_He4_Fe54_to_p_Co57_derived,
        H1, Co57, -1, Ni58, -1, -1, Co57, Ni58, k_p_Co57_to_Ni58_starlib, k_Ni58_to_p_Co57_derived,
        He4, Co55, -1, H1, Ni58, -1, Co55, Ni58, k_He4_Co55_to_p_Ni58_starlib, k_p_Ni58_to_He4_Co55_derived,
        N, Ni56, -1, H1, Co56, -1, Ni56, Co56, k_n_Ni56_to_p_Co56_starlib, k_p_Co56_to_n_Ni56_derived,
        H1, Ni58, -1, Cu59, -1, -1, Ni58, Cu59, k_p_Ni58_to_Cu59_starlib, k_Cu59_to_p_Ni58_derived,
        N, Cu59, -1, He4, Co56, -1, Cu59, Co56, k_n_Cu59_to_He4_Co56_starlib, k_He4_Co56_to_n_Cu59_derived,
        S32, He4, -1, Ar36, -1, -1, S32, Ar36, k_S32_He4_to_Ar36_approx, k_Ar36_to_S32_He4_approx,
        Ar36, He4, -1, Ca40, -1, -1, Ar36, Ca40, k_Ar36_He4_to_Ca40_approx, k_Ca40_to_Ar36_He4_approx,
        Ca40, He4, -1, Ti44, -1, -1, Ca40, Ti44, k_Ca40_He4_to_Ti44_approx, k_Ti44_to_Ca40_He4_approx,
        Ti44, He4, -1, Cr48, -1, -1, Ti44, Cr48, k_Ti44_He4_to_Cr48_approx, k_Cr48_to_Ti44_He4_approx
    };
}
#endif

#endif
