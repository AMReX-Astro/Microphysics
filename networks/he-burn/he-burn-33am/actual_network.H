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
        else if constexpr (spec == N13) {
            return 94.10522604799917_rt;
        }
        else if constexpr (spec == N14) {
            return 104.65860734799753_rt;
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
        else if constexpr (spec == Mn55) {
            return 482.07635959999607_rt;
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
        else if constexpr (spec == N13) {
            return 2.1596537049448796e-23_rt;
        }
        else if constexpr (spec == N14) {
            return 2.3252651436495096e-23_rt;
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
        else if constexpr (spec == Mn55) {
            return 9.122676671554189e-23_rt;
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
    std::uint8_t
    {
        k_C12_p_to_N13_reaclib = 1,
        k_C12_He4_to_O16_reaclib = 2,
        k_O16_He4_to_Ne20_reaclib = 3,
        k_Ne20_He4_to_Mg24_reaclib = 4,
        k_Na23_p_to_Mg24_reaclib = 5,
        k_Mg24_He4_to_Si28_reaclib = 6,
        k_Al27_p_to_Si28_reaclib = 7,
        k_Si28_He4_to_S32_reaclib = 8,
        k_P31_p_to_S32_reaclib = 9,
        k_Cr48_He4_to_Fe52_reaclib = 10,
        k_Mn51_p_to_Fe52_reaclib = 11,
        k_Mn51_He4_to_Co55_reaclib = 12,
        k_Fe52_He4_to_Ni56_reaclib = 13,
        k_Co55_p_to_Ni56_reaclib = 14,
        k_Co55_He4_to_Cu59_reaclib = 15,
        k_Ni56_He4_to_Zn60_reaclib = 16,
        k_Cu59_p_to_Zn60_reaclib = 17,
        k_C12_C12_to_p_Na23_reaclib = 18,
        k_C12_C12_to_He4_Ne20_reaclib = 19,
        k_N13_He4_to_p_O16_reaclib = 20,
        k_O16_C12_to_p_Al27_reaclib = 21,
        k_O16_C12_to_He4_Mg24_reaclib = 22,
        k_O16_O16_to_p_P31_reaclib = 23,
        k_O16_O16_to_He4_Si28_reaclib = 24,
        k_Na23_p_to_He4_Ne20_reaclib = 25,
        k_Al27_p_to_He4_Mg24_reaclib = 26,
        k_P31_p_to_He4_Si28_reaclib = 27,
        k_Cr48_He4_to_p_Mn51_reaclib = 28,
        k_Fe52_He4_to_p_Co55_reaclib = 29,
        k_Cu59_p_to_He4_Ni56_reaclib = 30,
        k_He4_He4_He4_to_C12_reaclib = 31,
        k_N14_He4_to_F18_reaclib = 32,
        k_O16_p_to_F17_reaclib = 33,
        k_C12_C12_to_n_Mg23_reaclib = 34,
        k_O16_O16_to_n_S31_reaclib = 35,
        k_O16_C12_to_n_Si27_reaclib = 36,
        k_Mn55_p_to_Fe56_reaclib = 37,
        k_Fe52_n_to_Fe53_reaclib = 38,
        k_Fe53_n_to_Fe54_reaclib = 39,
        k_Fe53_He4_to_Ni57_reaclib = 40,
        k_Fe54_n_to_Fe55_reaclib = 41,
        k_Fe54_p_to_Co55_reaclib = 42,
        k_Fe54_He4_to_Ni58_reaclib = 43,
        k_Fe55_n_to_Fe56_reaclib = 44,
        k_Fe55_p_to_Co56_reaclib = 45,
        k_Fe56_p_to_Co57_reaclib = 46,
        k_Co55_n_to_Co56_reaclib = 47,
        k_Co56_n_to_Co57_reaclib = 48,
        k_Co56_p_to_Ni57_reaclib = 49,
        k_Co57_p_to_Ni58_reaclib = 50,
        k_Ni56_n_to_Ni57_reaclib = 51,
        k_Ni57_n_to_Ni58_reaclib = 52,
        k_Ni58_p_to_Cu59_reaclib = 53,
        k_Mn51_He4_to_p_Fe54_reaclib = 54,
        k_Fe53_He4_to_p_Co56_reaclib = 55,
        k_Fe55_n_to_p_Mn55_reaclib = 56,
        k_Co55_n_to_p_Fe55_reaclib = 57,
        k_Co55_He4_to_p_Ni58_reaclib = 58,
        k_Co56_n_to_p_Fe56_reaclib = 59,
        k_Co57_p_to_He4_Fe54_reaclib = 60,
        k_Ni56_n_to_p_Co56_reaclib = 61,
        k_Ni56_n_to_He4_Fe53_reaclib = 62,
        k_Ni57_n_to_p_Co57_reaclib = 63,
        k_Ni57_n_to_He4_Fe54_reaclib = 64,
        k_Ni58_n_to_He4_Fe55_reaclib = 65,
        k_Cu59_n_to_He4_Co56_reaclib = 66,
        k_Zn60_n_to_He4_Ni57_reaclib = 67,
        k_S32_He4_to_Ar36_reaclib = 68,
        k_Cl35_p_to_He4_S32_reaclib = 69,
        k_Cl35_p_to_Ar36_reaclib = 70,
        k_Ar36_He4_to_Ca40_reaclib = 71,
        k_K39_p_to_He4_Ar36_reaclib = 72,
        k_K39_p_to_Ca40_reaclib = 73,
        k_Ca40_He4_to_Ti44_reaclib = 74,
        k_Sc43_p_to_He4_Ca40_reaclib = 75,
        k_Sc43_p_to_Ti44_reaclib = 76,
        k_Ti44_He4_to_p_V47_reaclib = 77,
        k_Ti44_He4_to_Cr48_reaclib = 78,
        k_V47_p_to_Cr48_reaclib = 79,
        k_Co55_to_Fe55_electron_capture_weaktab = 80,
        k_Fe55_to_Co55_beta_neg_weaktab = 81,
        k_Fe55_to_Mn55_electron_capture_weaktab = 82,
        k_Mn55_to_Fe55_beta_neg_weaktab = 83,
        k_Co56_to_Fe56_electron_capture_weaktab = 84,
        k_Co56_to_Ni56_beta_neg_weaktab = 85,
        k_Fe56_to_Co56_beta_neg_weaktab = 86,
        k_Ni56_to_Co56_electron_capture_weaktab = 87,
        k_Co57_to_Ni57_beta_neg_weaktab = 88,
        k_Ni57_to_Co57_electron_capture_weaktab = 89,
        k_n_to_p_beta_neg_weaktab = 90,
        k_p_to_n_electron_capture_weaktab = 91,
        k_S32_He4_to_Ar36_approx = 92,
        k_Ar36_to_S32_He4_approx = 93,
        k_Ar36_He4_to_Ca40_approx = 94,
        k_Ca40_to_Ar36_He4_approx = 95,
        k_Ca40_He4_to_Ti44_approx = 96,
        k_Ti44_to_Ca40_He4_approx = 97,
        k_Ti44_He4_to_Cr48_approx = 98,
        k_Cr48_to_Ti44_He4_approx = 99,
        k_N14_1p5He4_to_Ne20_modified = 100,
        k_O16_p_p_to_N14_He4_beta_pos_modified = 101,
        k_C12_C12_to_Mg24_modified = 102,
        k_O16_O16_to_S32_modified = 103,
        k_O16_C12_to_Si28_modified = 104,
        k_N13_to_p_C12_derived = 105,
        k_O16_to_He4_C12_derived = 106,
        k_Ne20_to_He4_O16_derived = 107,
        k_Mg24_to_p_Na23_derived = 108,
        k_Mg24_to_He4_Ne20_derived = 109,
        k_Si28_to_p_Al27_derived = 110,
        k_Si28_to_He4_Mg24_derived = 111,
        k_S32_to_p_P31_derived = 112,
        k_S32_to_He4_Si28_derived = 113,
        k_Fe52_to_p_Mn51_derived = 114,
        k_Fe52_to_He4_Cr48_derived = 115,
        k_Co55_to_He4_Mn51_derived = 116,
        k_Ni56_to_p_Co55_derived = 117,
        k_Ni56_to_He4_Fe52_derived = 118,
        k_Cu59_to_He4_Co55_derived = 119,
        k_Zn60_to_p_Cu59_derived = 120,
        k_Zn60_to_He4_Ni56_derived = 121,
        k_C12_to_3He4_derived = 122,
        k_O16_p_to_He4_N13_derived = 123,
        k_Ne20_He4_to_p_Na23_derived = 124,
        k_Mg24_He4_to_p_Al27_derived = 125,
        k_Si28_He4_to_p_P31_derived = 126,
        k_Mn51_p_to_He4_Cr48_derived = 127,
        k_Co55_p_to_He4_Fe52_derived = 128,
        k_Ni56_He4_to_p_Cu59_derived = 129,
        k_Fe53_to_n_Fe52_derived = 130,
        k_Fe54_to_n_Fe53_derived = 131,
        k_Fe55_to_n_Fe54_derived = 132,
        k_Fe56_to_n_Fe55_derived = 133,
        k_Fe56_to_p_Mn55_derived = 134,
        k_Co55_to_p_Fe54_derived = 135,
        k_Co56_to_n_Co55_derived = 136,
        k_Co56_to_p_Fe55_derived = 137,
        k_Co57_to_n_Co56_derived = 138,
        k_Co57_to_p_Fe56_derived = 139,
        k_Ni57_to_n_Ni56_derived = 140,
        k_Ni57_to_p_Co56_derived = 141,
        k_Ni57_to_He4_Fe53_derived = 142,
        k_Ni58_to_n_Ni57_derived = 143,
        k_Ni58_to_p_Co57_derived = 144,
        k_Ni58_to_He4_Fe54_derived = 145,
        k_Cu59_to_p_Ni58_derived = 146,
        k_Mn55_p_to_n_Fe55_derived = 147,
        k_Fe53_He4_to_n_Ni56_derived = 148,
        k_Fe54_p_to_He4_Mn51_derived = 149,
        k_Fe54_He4_to_n_Ni57_derived = 150,
        k_Fe54_He4_to_p_Co57_derived = 151,
        k_Fe55_p_to_n_Co55_derived = 152,
        k_Fe55_He4_to_n_Ni58_derived = 153,
        k_Fe56_p_to_n_Co56_derived = 154,
        k_Co56_p_to_n_Ni56_derived = 155,
        k_Co56_p_to_He4_Fe53_derived = 156,
        k_Co56_He4_to_n_Cu59_derived = 157,
        k_Co57_p_to_n_Ni57_derived = 158,
        k_Ni57_He4_to_n_Zn60_derived = 159,
        k_Ni58_p_to_He4_Co55_derived = 160,
        k_S32_He4_to_p_Cl35_derived = 161,
        k_Ar36_to_He4_S32_derived = 162,
        k_Ar36_to_p_Cl35_derived = 163,
        k_Ar36_He4_to_p_K39_derived = 164,
        k_Ca40_to_He4_Ar36_derived = 165,
        k_Ca40_to_p_K39_derived = 166,
        k_Ca40_He4_to_p_Sc43_derived = 167,
        k_Ti44_to_He4_Ca40_derived = 168,
        k_Ti44_to_p_Sc43_derived = 169,
        k_V47_p_to_He4_Ti44_derived = 170,
        k_Cr48_to_He4_Ti44_derived = 171,
        k_Cr48_to_p_V47_derived = 172,
        NumRates = k_Cr48_to_p_V47_derived
    };

    enum ScreenPairs :
    std::uint8_t
    {
        k_p_C12 = 1,
        k_He4_O16 = 2,
        k_He4_N14 = 3,
        k_p_Co57 = 4,
        k_p_K39 = 5,
        k_p_Fe55 = 6,
        k_He4_Ca40 = 7,
        k_He4_Ar36 = 8,
        k_p_Sc43 = 9,
        k_p_Fe56 = 10,
        k_He4_C12 = 11,
        k_He4_Ti44 = 12,
        k_He4_Fe55 = 13,
        k_p_Mn55 = 14,
        k_He4_Mg24 = 15,
        k_C12_O16 = 16,
        k_p_V47 = 17,
        k_He4_Cr48 = 18,
        k_p_Co55 = 19,
        k_C12_C12 = 20,
        k_He4_Ne20 = 21,
        k_p_Al27 = 22,
        k_p_Mn51 = 23,
        k_p_Na23 = 24,
        k_p_Co56 = 25,
        k_p_Fe54 = 26,
        k_He4_Fe52 = 27,
        k_p_Cu59 = 28,
        k_He4_S32 = 29,
        k_He4_Mn51 = 30,
        k_He4_Co55 = 31,
        k_He4_Ni56 = 32,
        k_He4_Fe53 = 33,
        k_p_P31 = 34,
        k_He4_He4 = 35,
        k_p_Ni58 = 36,
        k_p_O16 = 37,
        k_He4_Co56 = 38,
        k_He4_Si28 = 39,
        k_He4_Fe54 = 40,
        k_He4_Ni57 = 41,
        k_p_Cl35 = 42,
        k_O16_O16 = 43,
        k_He4_Be8 = 44,
        k_He4_N13 = 45,
        NumScreenPairs = k_He4_N13
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "C12_p_to_N13_reaclib",  // 1,
        "C12_He4_to_O16_reaclib",  // 2,
        "O16_He4_to_Ne20_reaclib",  // 3,
        "Ne20_He4_to_Mg24_reaclib",  // 4,
        "Na23_p_to_Mg24_reaclib",  // 5,
        "Mg24_He4_to_Si28_reaclib",  // 6,
        "Al27_p_to_Si28_reaclib",  // 7,
        "Si28_He4_to_S32_reaclib",  // 8,
        "P31_p_to_S32_reaclib",  // 9,
        "Cr48_He4_to_Fe52_reaclib",  // 10,
        "Mn51_p_to_Fe52_reaclib",  // 11,
        "Mn51_He4_to_Co55_reaclib",  // 12,
        "Fe52_He4_to_Ni56_reaclib",  // 13,
        "Co55_p_to_Ni56_reaclib",  // 14,
        "Co55_He4_to_Cu59_reaclib",  // 15,
        "Ni56_He4_to_Zn60_reaclib",  // 16,
        "Cu59_p_to_Zn60_reaclib",  // 17,
        "C12_C12_to_p_Na23_reaclib",  // 18,
        "C12_C12_to_He4_Ne20_reaclib",  // 19,
        "N13_He4_to_p_O16_reaclib",  // 20,
        "O16_C12_to_p_Al27_reaclib",  // 21,
        "O16_C12_to_He4_Mg24_reaclib",  // 22,
        "O16_O16_to_p_P31_reaclib",  // 23,
        "O16_O16_to_He4_Si28_reaclib",  // 24,
        "Na23_p_to_He4_Ne20_reaclib",  // 25,
        "Al27_p_to_He4_Mg24_reaclib",  // 26,
        "P31_p_to_He4_Si28_reaclib",  // 27,
        "Cr48_He4_to_p_Mn51_reaclib",  // 28,
        "Fe52_He4_to_p_Co55_reaclib",  // 29,
        "Cu59_p_to_He4_Ni56_reaclib",  // 30,
        "He4_He4_He4_to_C12_reaclib",  // 31,
        "N14_He4_to_F18_reaclib",  // 32,
        "O16_p_to_F17_reaclib",  // 33,
        "C12_C12_to_n_Mg23_reaclib",  // 34,
        "O16_O16_to_n_S31_reaclib",  // 35,
        "O16_C12_to_n_Si27_reaclib",  // 36,
        "Mn55_p_to_Fe56_reaclib",  // 37,
        "Fe52_n_to_Fe53_reaclib",  // 38,
        "Fe53_n_to_Fe54_reaclib",  // 39,
        "Fe53_He4_to_Ni57_reaclib",  // 40,
        "Fe54_n_to_Fe55_reaclib",  // 41,
        "Fe54_p_to_Co55_reaclib",  // 42,
        "Fe54_He4_to_Ni58_reaclib",  // 43,
        "Fe55_n_to_Fe56_reaclib",  // 44,
        "Fe55_p_to_Co56_reaclib",  // 45,
        "Fe56_p_to_Co57_reaclib",  // 46,
        "Co55_n_to_Co56_reaclib",  // 47,
        "Co56_n_to_Co57_reaclib",  // 48,
        "Co56_p_to_Ni57_reaclib",  // 49,
        "Co57_p_to_Ni58_reaclib",  // 50,
        "Ni56_n_to_Ni57_reaclib",  // 51,
        "Ni57_n_to_Ni58_reaclib",  // 52,
        "Ni58_p_to_Cu59_reaclib",  // 53,
        "Mn51_He4_to_p_Fe54_reaclib",  // 54,
        "Fe53_He4_to_p_Co56_reaclib",  // 55,
        "Fe55_n_to_p_Mn55_reaclib",  // 56,
        "Co55_n_to_p_Fe55_reaclib",  // 57,
        "Co55_He4_to_p_Ni58_reaclib",  // 58,
        "Co56_n_to_p_Fe56_reaclib",  // 59,
        "Co57_p_to_He4_Fe54_reaclib",  // 60,
        "Ni56_n_to_p_Co56_reaclib",  // 61,
        "Ni56_n_to_He4_Fe53_reaclib",  // 62,
        "Ni57_n_to_p_Co57_reaclib",  // 63,
        "Ni57_n_to_He4_Fe54_reaclib",  // 64,
        "Ni58_n_to_He4_Fe55_reaclib",  // 65,
        "Cu59_n_to_He4_Co56_reaclib",  // 66,
        "Zn60_n_to_He4_Ni57_reaclib",  // 67,
        "S32_He4_to_Ar36_reaclib",  // 68,
        "Cl35_p_to_He4_S32_reaclib",  // 69,
        "Cl35_p_to_Ar36_reaclib",  // 70,
        "Ar36_He4_to_Ca40_reaclib",  // 71,
        "K39_p_to_He4_Ar36_reaclib",  // 72,
        "K39_p_to_Ca40_reaclib",  // 73,
        "Ca40_He4_to_Ti44_reaclib",  // 74,
        "Sc43_p_to_He4_Ca40_reaclib",  // 75,
        "Sc43_p_to_Ti44_reaclib",  // 76,
        "Ti44_He4_to_p_V47_reaclib",  // 77,
        "Ti44_He4_to_Cr48_reaclib",  // 78,
        "V47_p_to_Cr48_reaclib",  // 79,
        "Co55_to_Fe55_electron_capture_weaktab",  // 80,
        "Fe55_to_Co55_beta_neg_weaktab",  // 81,
        "Fe55_to_Mn55_electron_capture_weaktab",  // 82,
        "Mn55_to_Fe55_beta_neg_weaktab",  // 83,
        "Co56_to_Fe56_electron_capture_weaktab",  // 84,
        "Co56_to_Ni56_beta_neg_weaktab",  // 85,
        "Fe56_to_Co56_beta_neg_weaktab",  // 86,
        "Ni56_to_Co56_electron_capture_weaktab",  // 87,
        "Co57_to_Ni57_beta_neg_weaktab",  // 88,
        "Ni57_to_Co57_electron_capture_weaktab",  // 89,
        "n_to_p_beta_neg_weaktab",  // 90,
        "p_to_n_electron_capture_weaktab",  // 91,
        "S32_He4_to_Ar36_approx",  // 92,
        "Ar36_to_S32_He4_approx",  // 93,
        "Ar36_He4_to_Ca40_approx",  // 94,
        "Ca40_to_Ar36_He4_approx",  // 95,
        "Ca40_He4_to_Ti44_approx",  // 96,
        "Ti44_to_Ca40_He4_approx",  // 97,
        "Ti44_He4_to_Cr48_approx",  // 98,
        "Cr48_to_Ti44_He4_approx",  // 99,
        "N14_1p5He4_to_Ne20_modified",  // 100,
        "O16_p_p_to_N14_He4_beta_pos_modified",  // 101,
        "C12_C12_to_Mg24_modified",  // 102,
        "O16_O16_to_S32_modified",  // 103,
        "O16_C12_to_Si28_modified",  // 104,
        "N13_to_p_C12_derived",  // 105,
        "O16_to_He4_C12_derived",  // 106,
        "Ne20_to_He4_O16_derived",  // 107,
        "Mg24_to_p_Na23_derived",  // 108,
        "Mg24_to_He4_Ne20_derived",  // 109,
        "Si28_to_p_Al27_derived",  // 110,
        "Si28_to_He4_Mg24_derived",  // 111,
        "S32_to_p_P31_derived",  // 112,
        "S32_to_He4_Si28_derived",  // 113,
        "Fe52_to_p_Mn51_derived",  // 114,
        "Fe52_to_He4_Cr48_derived",  // 115,
        "Co55_to_He4_Mn51_derived",  // 116,
        "Ni56_to_p_Co55_derived",  // 117,
        "Ni56_to_He4_Fe52_derived",  // 118,
        "Cu59_to_He4_Co55_derived",  // 119,
        "Zn60_to_p_Cu59_derived",  // 120,
        "Zn60_to_He4_Ni56_derived",  // 121,
        "C12_to_3He4_derived",  // 122,
        "O16_p_to_He4_N13_derived",  // 123,
        "Ne20_He4_to_p_Na23_derived",  // 124,
        "Mg24_He4_to_p_Al27_derived",  // 125,
        "Si28_He4_to_p_P31_derived",  // 126,
        "Mn51_p_to_He4_Cr48_derived",  // 127,
        "Co55_p_to_He4_Fe52_derived",  // 128,
        "Ni56_He4_to_p_Cu59_derived",  // 129,
        "Fe53_to_n_Fe52_derived",  // 130,
        "Fe54_to_n_Fe53_derived",  // 131,
        "Fe55_to_n_Fe54_derived",  // 132,
        "Fe56_to_n_Fe55_derived",  // 133,
        "Fe56_to_p_Mn55_derived",  // 134,
        "Co55_to_p_Fe54_derived",  // 135,
        "Co56_to_n_Co55_derived",  // 136,
        "Co56_to_p_Fe55_derived",  // 137,
        "Co57_to_n_Co56_derived",  // 138,
        "Co57_to_p_Fe56_derived",  // 139,
        "Ni57_to_n_Ni56_derived",  // 140,
        "Ni57_to_p_Co56_derived",  // 141,
        "Ni57_to_He4_Fe53_derived",  // 142,
        "Ni58_to_n_Ni57_derived",  // 143,
        "Ni58_to_p_Co57_derived",  // 144,
        "Ni58_to_He4_Fe54_derived",  // 145,
        "Cu59_to_p_Ni58_derived",  // 146,
        "Mn55_p_to_n_Fe55_derived",  // 147,
        "Fe53_He4_to_n_Ni56_derived",  // 148,
        "Fe54_p_to_He4_Mn51_derived",  // 149,
        "Fe54_He4_to_n_Ni57_derived",  // 150,
        "Fe54_He4_to_p_Co57_derived",  // 151,
        "Fe55_p_to_n_Co55_derived",  // 152,
        "Fe55_He4_to_n_Ni58_derived",  // 153,
        "Fe56_p_to_n_Co56_derived",  // 154,
        "Co56_p_to_n_Ni56_derived",  // 155,
        "Co56_p_to_He4_Fe53_derived",  // 156,
        "Co56_He4_to_n_Cu59_derived",  // 157,
        "Co57_p_to_n_Ni57_derived",  // 158,
        "Ni57_He4_to_n_Zn60_derived",  // 159,
        "Ni58_p_to_He4_Co55_derived",  // 160,
        "S32_He4_to_p_Cl35_derived",  // 161,
        "Ar36_to_He4_S32_derived",  // 162,
        "Ar36_to_p_Cl35_derived",  // 163,
        "Ar36_He4_to_p_K39_derived",  // 164,
        "Ca40_to_He4_Ar36_derived",  // 165,
        "Ca40_to_p_K39_derived",  // 166,
        "Ca40_He4_to_p_Sc43_derived",  // 167,
        "Ti44_to_He4_Ca40_derived",  // 168,
        "Ti44_to_p_Sc43_derived",  // 169,
        "V47_p_to_He4_Ti44_derived",  // 170,
        "Cr48_to_He4_Ti44_derived",  // 171,
        "Cr48_to_p_V47_derived"  // 172,
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
    constexpr int NumNSERatePairs = 60;

    inline AMREX_GPU_MANAGED amrex::Array2D<short, 1, NumNSERatePairs, 1, 10, amrex::Order::C> rate_pair_data {
        H1, C12, -1, N13, -1, -1, C12, N13, k_C12_p_to_N13_reaclib, k_N13_to_p_C12_derived,
        He4, C12, -1, O16, -1, -1, C12, O16, k_C12_He4_to_O16_reaclib, k_O16_to_He4_C12_derived,
        He4, O16, -1, Ne20, -1, -1, O16, Ne20, k_O16_He4_to_Ne20_reaclib, k_Ne20_to_He4_O16_derived,
        He4, Ne20, -1, Mg24, -1, -1, Ne20, Mg24, k_Ne20_He4_to_Mg24_reaclib, k_Mg24_to_He4_Ne20_derived,
        H1, Na23, -1, Mg24, -1, -1, Na23, Mg24, k_Na23_p_to_Mg24_reaclib, k_Mg24_to_p_Na23_derived,
        He4, Mg24, -1, Si28, -1, -1, Mg24, Si28, k_Mg24_He4_to_Si28_reaclib, k_Si28_to_He4_Mg24_derived,
        H1, Al27, -1, Si28, -1, -1, Al27, Si28, k_Al27_p_to_Si28_reaclib, k_Si28_to_p_Al27_derived,
        He4, Si28, -1, S32, -1, -1, Si28, S32, k_Si28_He4_to_S32_reaclib, k_S32_to_He4_Si28_derived,
        H1, P31, -1, S32, -1, -1, P31, S32, k_P31_p_to_S32_reaclib, k_S32_to_p_P31_derived,
        He4, Cr48, -1, Fe52, -1, -1, Cr48, Fe52, k_Cr48_He4_to_Fe52_reaclib, k_Fe52_to_He4_Cr48_derived,
        H1, Mn51, -1, Fe52, -1, -1, Mn51, Fe52, k_Mn51_p_to_Fe52_reaclib, k_Fe52_to_p_Mn51_derived,
        He4, Mn51, -1, Co55, -1, -1, Mn51, Co55, k_Mn51_He4_to_Co55_reaclib, k_Co55_to_He4_Mn51_derived,
        He4, Fe52, -1, Ni56, -1, -1, Fe52, Ni56, k_Fe52_He4_to_Ni56_reaclib, k_Ni56_to_He4_Fe52_derived,
        H1, Co55, -1, Ni56, -1, -1, Co55, Ni56, k_Co55_p_to_Ni56_reaclib, k_Ni56_to_p_Co55_derived,
        He4, Co55, -1, Cu59, -1, -1, Co55, Cu59, k_Co55_He4_to_Cu59_reaclib, k_Cu59_to_He4_Co55_derived,
        He4, Ni56, -1, Zn60, -1, -1, Ni56, Zn60, k_Ni56_He4_to_Zn60_reaclib, k_Zn60_to_He4_Ni56_derived,
        H1, Cu59, -1, Zn60, -1, -1, Cu59, Zn60, k_Cu59_p_to_Zn60_reaclib, k_Zn60_to_p_Cu59_derived,
        He4, N13, -1, H1, O16, -1, N13, O16, k_N13_He4_to_p_O16_reaclib, k_O16_p_to_He4_N13_derived,
        H1, Na23, -1, He4, Ne20, -1, Na23, Ne20, k_Na23_p_to_He4_Ne20_reaclib, k_Ne20_He4_to_p_Na23_derived,
        H1, Al27, -1, He4, Mg24, -1, Al27, Mg24, k_Al27_p_to_He4_Mg24_reaclib, k_Mg24_He4_to_p_Al27_derived,
        H1, P31, -1, He4, Si28, -1, P31, Si28, k_P31_p_to_He4_Si28_reaclib, k_Si28_He4_to_p_P31_derived,
        He4, Cr48, -1, H1, Mn51, -1, Cr48, Mn51, k_Cr48_He4_to_p_Mn51_reaclib, k_Mn51_p_to_He4_Cr48_derived,
        He4, Fe52, -1, H1, Co55, -1, Fe52, Co55, k_Fe52_He4_to_p_Co55_reaclib, k_Co55_p_to_He4_Fe52_derived,
        H1, Cu59, -1, He4, Ni56, -1, Cu59, Ni56, k_Cu59_p_to_He4_Ni56_reaclib, k_Ni56_He4_to_p_Cu59_derived,
        He4, He4, He4, C12, -1, -1, C12, -1, k_He4_He4_He4_to_C12_reaclib, k_C12_to_3He4_derived,
        H1, Mn55, -1, Fe56, -1, -1, Mn55, Fe56, k_Mn55_p_to_Fe56_reaclib, k_Fe56_to_p_Mn55_derived,
        N, Fe52, -1, Fe53, -1, -1, Fe52, Fe53, k_Fe52_n_to_Fe53_reaclib, k_Fe53_to_n_Fe52_derived,
        N, Fe53, -1, Fe54, -1, -1, Fe53, Fe54, k_Fe53_n_to_Fe54_reaclib, k_Fe54_to_n_Fe53_derived,
        He4, Fe53, -1, Ni57, -1, -1, Fe53, Ni57, k_Fe53_He4_to_Ni57_reaclib, k_Ni57_to_He4_Fe53_derived,
        N, Fe54, -1, Fe55, -1, -1, Fe54, Fe55, k_Fe54_n_to_Fe55_reaclib, k_Fe55_to_n_Fe54_derived,
        H1, Fe54, -1, Co55, -1, -1, Fe54, Co55, k_Fe54_p_to_Co55_reaclib, k_Co55_to_p_Fe54_derived,
        He4, Fe54, -1, Ni58, -1, -1, Fe54, Ni58, k_Fe54_He4_to_Ni58_reaclib, k_Ni58_to_He4_Fe54_derived,
        N, Fe55, -1, Fe56, -1, -1, Fe55, Fe56, k_Fe55_n_to_Fe56_reaclib, k_Fe56_to_n_Fe55_derived,
        H1, Fe55, -1, Co56, -1, -1, Fe55, Co56, k_Fe55_p_to_Co56_reaclib, k_Co56_to_p_Fe55_derived,
        H1, Fe56, -1, Co57, -1, -1, Fe56, Co57, k_Fe56_p_to_Co57_reaclib, k_Co57_to_p_Fe56_derived,
        N, Co55, -1, Co56, -1, -1, Co55, Co56, k_Co55_n_to_Co56_reaclib, k_Co56_to_n_Co55_derived,
        N, Co56, -1, Co57, -1, -1, Co56, Co57, k_Co56_n_to_Co57_reaclib, k_Co57_to_n_Co56_derived,
        H1, Co56, -1, Ni57, -1, -1, Co56, Ni57, k_Co56_p_to_Ni57_reaclib, k_Ni57_to_p_Co56_derived,
        H1, Co57, -1, Ni58, -1, -1, Co57, Ni58, k_Co57_p_to_Ni58_reaclib, k_Ni58_to_p_Co57_derived,
        N, Ni56, -1, Ni57, -1, -1, Ni56, Ni57, k_Ni56_n_to_Ni57_reaclib, k_Ni57_to_n_Ni56_derived,
        N, Ni57, -1, Ni58, -1, -1, Ni57, Ni58, k_Ni57_n_to_Ni58_reaclib, k_Ni58_to_n_Ni57_derived,
        H1, Ni58, -1, Cu59, -1, -1, Ni58, Cu59, k_Ni58_p_to_Cu59_reaclib, k_Cu59_to_p_Ni58_derived,
        He4, Mn51, -1, H1, Fe54, -1, Mn51, Fe54, k_Mn51_He4_to_p_Fe54_reaclib, k_Fe54_p_to_He4_Mn51_derived,
        He4, Fe53, -1, H1, Co56, -1, Fe53, Co56, k_Fe53_He4_to_p_Co56_reaclib, k_Co56_p_to_He4_Fe53_derived,
        N, Fe55, -1, H1, Mn55, -1, Fe55, Mn55, k_Fe55_n_to_p_Mn55_reaclib, k_Mn55_p_to_n_Fe55_derived,
        N, Co55, -1, H1, Fe55, -1, Co55, Fe55, k_Co55_n_to_p_Fe55_reaclib, k_Fe55_p_to_n_Co55_derived,
        He4, Co55, -1, H1, Ni58, -1, Co55, Ni58, k_Co55_He4_to_p_Ni58_reaclib, k_Ni58_p_to_He4_Co55_derived,
        N, Co56, -1, H1, Fe56, -1, Co56, Fe56, k_Co56_n_to_p_Fe56_reaclib, k_Fe56_p_to_n_Co56_derived,
        H1, Co57, -1, He4, Fe54, -1, Co57, Fe54, k_Co57_p_to_He4_Fe54_reaclib, k_Fe54_He4_to_p_Co57_derived,
        N, Ni56, -1, H1, Co56, -1, Ni56, Co56, k_Ni56_n_to_p_Co56_reaclib, k_Co56_p_to_n_Ni56_derived,
        N, Ni56, -1, He4, Fe53, -1, Ni56, Fe53, k_Ni56_n_to_He4_Fe53_reaclib, k_Fe53_He4_to_n_Ni56_derived,
        N, Ni57, -1, H1, Co57, -1, Ni57, Co57, k_Ni57_n_to_p_Co57_reaclib, k_Co57_p_to_n_Ni57_derived,
        N, Ni57, -1, He4, Fe54, -1, Ni57, Fe54, k_Ni57_n_to_He4_Fe54_reaclib, k_Fe54_He4_to_n_Ni57_derived,
        N, Ni58, -1, He4, Fe55, -1, Ni58, Fe55, k_Ni58_n_to_He4_Fe55_reaclib, k_Fe55_He4_to_n_Ni58_derived,
        N, Cu59, -1, He4, Co56, -1, Cu59, Co56, k_Cu59_n_to_He4_Co56_reaclib, k_Co56_He4_to_n_Cu59_derived,
        N, Zn60, -1, He4, Ni57, -1, Zn60, Ni57, k_Zn60_n_to_He4_Ni57_reaclib, k_Ni57_He4_to_n_Zn60_derived,
        S32, He4, -1, Ar36, -1, -1, S32, Ar36, k_S32_He4_to_Ar36_approx, k_Ar36_to_S32_He4_approx,
        Ar36, He4, -1, Ca40, -1, -1, Ar36, Ca40, k_Ar36_He4_to_Ca40_approx, k_Ca40_to_Ar36_He4_approx,
        Ca40, He4, -1, Ti44, -1, -1, Ca40, Ti44, k_Ca40_He4_to_Ti44_approx, k_Ti44_to_Ca40_He4_approx,
        Ti44, He4, -1, Cr48, -1, -1, Ti44, Cr48, k_Ti44_He4_to_Cr48_approx, k_Cr48_to_Ti44_He4_approx
    };
}
#endif

#endif
