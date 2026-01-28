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
        k_p_C12_to_N13_reaclib = 1,
        k_He4_C12_to_O16_reaclib = 2,
        k_He4_O16_to_Ne20_reaclib = 3,
        k_He4_Ne20_to_Mg24_reaclib = 4,
        k_p_Na23_to_Mg24_reaclib = 5,
        k_He4_Mg24_to_Si28_reaclib = 6,
        k_p_Al27_to_Si28_reaclib = 7,
        k_He4_Si28_to_S32_reaclib = 8,
        k_p_P31_to_S32_reaclib = 9,
        k_He4_Cr48_to_Fe52_reaclib = 10,
        k_p_Mn51_to_Fe52_reaclib = 11,
        k_He4_Mn51_to_Co55_reaclib = 12,
        k_He4_Fe52_to_Ni56_reaclib = 13,
        k_p_Co55_to_Ni56_reaclib = 14,
        k_C12_C12_to_p_Na23_reaclib = 15,
        k_C12_C12_to_He4_Ne20_reaclib = 16,
        k_He4_N13_to_p_O16_reaclib = 17,
        k_C12_O16_to_p_Al27_reaclib = 18,
        k_C12_O16_to_He4_Mg24_reaclib = 19,
        k_O16_O16_to_p_P31_reaclib = 20,
        k_O16_O16_to_He4_Si28_reaclib = 21,
        k_p_Na23_to_He4_Ne20_reaclib = 22,
        k_p_Al27_to_He4_Mg24_reaclib = 23,
        k_p_P31_to_He4_Si28_reaclib = 24,
        k_He4_Cr48_to_p_Mn51_reaclib = 25,
        k_He4_Fe52_to_p_Co55_reaclib = 26,
        k_He4_He4_He4_to_C12_reaclib = 27,
        k_C12_C12_to_n_Mg23_removed = 28,
        k_O16_O16_to_n_S31_removed = 29,
        k_C12_O16_to_n_Si27_removed = 30,
        k_p_Fe54_to_Co55_reaclib = 31,
        k_He4_Fe54_to_Ni58_reaclib = 32,
        k_p_Fe56_to_Co57_reaclib = 33,
        k_n_Co55_to_Co56_reaclib = 34,
        k_n_Co56_to_Co57_reaclib = 35,
        k_p_Co57_to_Ni58_reaclib = 36,
        k_He4_Mn51_to_p_Fe54_reaclib = 37,
        k_He4_Co55_to_p_Ni58_reaclib = 38,
        k_n_Co56_to_p_Fe56_reaclib = 39,
        k_p_Co57_to_He4_Fe54_reaclib = 40,
        k_n_Ni56_to_p_Co56_reaclib = 41,
        k_He4_S32_to_Ar36_removed = 42,
        k_p_Cl35_to_Ar36_removed = 43,
        k_p_Cl35_to_He4_S32_removed = 44,
        k_He4_Ar36_to_Ca40_removed = 45,
        k_p_K39_to_Ca40_removed = 46,
        k_p_K39_to_He4_Ar36_removed = 47,
        k_He4_Ca40_to_Ti44_removed = 48,
        k_p_Sc43_to_Ti44_removed = 49,
        k_p_Sc43_to_He4_Ca40_removed = 50,
        k_He4_Ti44_to_Cr48_removed = 51,
        k_He4_Ti44_to_p_V47_removed = 52,
        k_p_V47_to_Cr48_removed = 53,
        k_n_Fe52_to_Fe53_removed = 54,
        k_n_Fe53_to_Fe54_removed = 55,
        k_n_Fe54_to_Fe55_removed = 56,
        k_n_Fe55_to_Fe56_removed = 57,
        k_n_Ni56_to_Ni57_removed = 58,
        k_n_Ni57_to_Ni58_removed = 59,
        k_Co56_to_Fe56_weaktab = 60,
        k_Co56_to_Ni56_weaktab = 61,
        k_Fe56_to_Co56_weaktab = 62,
        k_Ni56_to_Co56_weaktab = 63,
        k_n_to_p_weaktab = 64,
        k_p_to_n_weaktab = 65,
        k_S32_He4_to_Ar36_approx = 66,
        k_Ar36_to_S32_He4_approx = 67,
        k_Ar36_He4_to_Ca40_approx = 68,
        k_Ca40_to_Ar36_He4_approx = 69,
        k_Ca40_He4_to_Ti44_approx = 70,
        k_Ti44_to_Ca40_He4_approx = 71,
        k_Ti44_He4_to_Cr48_approx = 72,
        k_Cr48_to_Ti44_He4_approx = 73,
        k_Fe52_n_n_to_Fe54_approx = 74,
        k_Fe54_to_Fe52_n_n_approx = 75,
        k_Fe54_n_n_to_Fe56_approx = 76,
        k_Fe56_to_Fe54_n_n_approx = 77,
        k_Ni56_n_n_to_Ni58_approx = 78,
        k_Ni58_to_Ni56_n_n_approx = 79,
        k_C12_C12_to_Mg24_modified = 80,
        k_O16_O16_to_S32_modified = 81,
        k_C12_O16_to_Si28_modified = 82,
        k_N13_to_p_C12_derived = 83,
        k_O16_to_He4_C12_derived = 84,
        k_Ne20_to_He4_O16_derived = 85,
        k_Mg24_to_p_Na23_derived = 86,
        k_Mg24_to_He4_Ne20_derived = 87,
        k_Si28_to_p_Al27_derived = 88,
        k_Si28_to_He4_Mg24_derived = 89,
        k_S32_to_p_P31_derived = 90,
        k_S32_to_He4_Si28_derived = 91,
        k_Fe52_to_p_Mn51_derived = 92,
        k_Fe52_to_He4_Cr48_derived = 93,
        k_Co55_to_He4_Mn51_derived = 94,
        k_Ni56_to_p_Co55_derived = 95,
        k_Ni56_to_He4_Fe52_derived = 96,
        k_C12_to_He4_He4_He4_derived = 97,
        k_p_O16_to_He4_N13_derived = 98,
        k_He4_Ne20_to_p_Na23_derived = 99,
        k_He4_Ne20_to_C12_C12_derived = 100,
        k_p_Na23_to_C12_C12_derived = 101,
        k_He4_Mg24_to_p_Al27_derived = 102,
        k_He4_Mg24_to_C12_O16_derived = 103,
        k_p_Al27_to_C12_O16_derived = 104,
        k_He4_Si28_to_p_P31_derived = 105,
        k_He4_Si28_to_O16_O16_derived = 106,
        k_p_P31_to_O16_O16_derived = 107,
        k_p_Mn51_to_He4_Cr48_derived = 108,
        k_p_Co55_to_He4_Fe52_derived = 109,
        k_Mg24_to_C12_C12_derived = 110,
        k_S32_to_O16_O16_derived = 111,
        k_Si28_to_C12_O16_derived = 112,
        k_Co55_to_p_Fe54_derived = 113,
        k_Co56_to_n_Co55_derived = 114,
        k_Co57_to_n_Co56_derived = 115,
        k_Co57_to_p_Fe56_derived = 116,
        k_Ni58_to_p_Co57_derived = 117,
        k_Ni58_to_He4_Fe54_derived = 118,
        k_p_Fe54_to_He4_Mn51_derived = 119,
        k_He4_Fe54_to_p_Co57_derived = 120,
        k_p_Fe56_to_n_Co56_derived = 121,
        k_p_Co56_to_n_Ni56_derived = 122,
        k_p_Ni58_to_He4_Co55_derived = 123,
        k_He4_S32_to_p_Cl35_removed = 124,
        k_Ar36_to_He4_S32_removed = 125,
        k_Ar36_to_p_Cl35_removed = 126,
        k_He4_Ar36_to_p_K39_removed = 127,
        k_Ca40_to_He4_Ar36_removed = 128,
        k_Ca40_to_p_K39_removed = 129,
        k_He4_Ca40_to_p_Sc43_removed = 130,
        k_Ti44_to_He4_Ca40_removed = 131,
        k_Ti44_to_p_Sc43_removed = 132,
        k_Cr48_to_He4_Ti44_removed = 133,
        k_Cr48_to_p_V47_removed = 134,
        k_p_V47_to_He4_Ti44_removed = 135,
        k_Fe54_to_n_Fe53_removed = 136,
        k_Fe53_to_n_Fe52_removed = 137,
        k_Fe56_to_n_Fe55_removed = 138,
        k_Fe55_to_n_Fe54_removed = 139,
        k_Ni58_to_n_Ni57_removed = 140,
        k_Ni57_to_n_Ni56_removed = 141,
        NumRates = k_Ni57_to_n_Ni56_removed
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "p_C12_to_N13_reaclib",  // 1,
        "He4_C12_to_O16_reaclib",  // 2,
        "He4_O16_to_Ne20_reaclib",  // 3,
        "He4_Ne20_to_Mg24_reaclib",  // 4,
        "p_Na23_to_Mg24_reaclib",  // 5,
        "He4_Mg24_to_Si28_reaclib",  // 6,
        "p_Al27_to_Si28_reaclib",  // 7,
        "He4_Si28_to_S32_reaclib",  // 8,
        "p_P31_to_S32_reaclib",  // 9,
        "He4_Cr48_to_Fe52_reaclib",  // 10,
        "p_Mn51_to_Fe52_reaclib",  // 11,
        "He4_Mn51_to_Co55_reaclib",  // 12,
        "He4_Fe52_to_Ni56_reaclib",  // 13,
        "p_Co55_to_Ni56_reaclib",  // 14,
        "C12_C12_to_p_Na23_reaclib",  // 15,
        "C12_C12_to_He4_Ne20_reaclib",  // 16,
        "He4_N13_to_p_O16_reaclib",  // 17,
        "C12_O16_to_p_Al27_reaclib",  // 18,
        "C12_O16_to_He4_Mg24_reaclib",  // 19,
        "O16_O16_to_p_P31_reaclib",  // 20,
        "O16_O16_to_He4_Si28_reaclib",  // 21,
        "p_Na23_to_He4_Ne20_reaclib",  // 22,
        "p_Al27_to_He4_Mg24_reaclib",  // 23,
        "p_P31_to_He4_Si28_reaclib",  // 24,
        "He4_Cr48_to_p_Mn51_reaclib",  // 25,
        "He4_Fe52_to_p_Co55_reaclib",  // 26,
        "He4_He4_He4_to_C12_reaclib",  // 27,
        "C12_C12_to_n_Mg23_removed",  // 28,
        "O16_O16_to_n_S31_removed",  // 29,
        "C12_O16_to_n_Si27_removed",  // 30,
        "p_Fe54_to_Co55_reaclib",  // 31,
        "He4_Fe54_to_Ni58_reaclib",  // 32,
        "p_Fe56_to_Co57_reaclib",  // 33,
        "n_Co55_to_Co56_reaclib",  // 34,
        "n_Co56_to_Co57_reaclib",  // 35,
        "p_Co57_to_Ni58_reaclib",  // 36,
        "He4_Mn51_to_p_Fe54_reaclib",  // 37,
        "He4_Co55_to_p_Ni58_reaclib",  // 38,
        "n_Co56_to_p_Fe56_reaclib",  // 39,
        "p_Co57_to_He4_Fe54_reaclib",  // 40,
        "n_Ni56_to_p_Co56_reaclib",  // 41,
        "He4_S32_to_Ar36_removed",  // 42,
        "p_Cl35_to_Ar36_removed",  // 43,
        "p_Cl35_to_He4_S32_removed",  // 44,
        "He4_Ar36_to_Ca40_removed",  // 45,
        "p_K39_to_Ca40_removed",  // 46,
        "p_K39_to_He4_Ar36_removed",  // 47,
        "He4_Ca40_to_Ti44_removed",  // 48,
        "p_Sc43_to_Ti44_removed",  // 49,
        "p_Sc43_to_He4_Ca40_removed",  // 50,
        "He4_Ti44_to_Cr48_removed",  // 51,
        "He4_Ti44_to_p_V47_removed",  // 52,
        "p_V47_to_Cr48_removed",  // 53,
        "n_Fe52_to_Fe53_removed",  // 54,
        "n_Fe53_to_Fe54_removed",  // 55,
        "n_Fe54_to_Fe55_removed",  // 56,
        "n_Fe55_to_Fe56_removed",  // 57,
        "n_Ni56_to_Ni57_removed",  // 58,
        "n_Ni57_to_Ni58_removed",  // 59,
        "Co56_to_Fe56_weaktab",  // 60,
        "Co56_to_Ni56_weaktab",  // 61,
        "Fe56_to_Co56_weaktab",  // 62,
        "Ni56_to_Co56_weaktab",  // 63,
        "n_to_p_weaktab",  // 64,
        "p_to_n_weaktab",  // 65,
        "S32_He4_to_Ar36_approx",  // 66,
        "Ar36_to_S32_He4_approx",  // 67,
        "Ar36_He4_to_Ca40_approx",  // 68,
        "Ca40_to_Ar36_He4_approx",  // 69,
        "Ca40_He4_to_Ti44_approx",  // 70,
        "Ti44_to_Ca40_He4_approx",  // 71,
        "Ti44_He4_to_Cr48_approx",  // 72,
        "Cr48_to_Ti44_He4_approx",  // 73,
        "Fe52_n_n_to_Fe54_approx",  // 74,
        "Fe54_to_Fe52_n_n_approx",  // 75,
        "Fe54_n_n_to_Fe56_approx",  // 76,
        "Fe56_to_Fe54_n_n_approx",  // 77,
        "Ni56_n_n_to_Ni58_approx",  // 78,
        "Ni58_to_Ni56_n_n_approx",  // 79,
        "C12_C12_to_Mg24_modified",  // 80,
        "O16_O16_to_S32_modified",  // 81,
        "C12_O16_to_Si28_modified",  // 82,
        "N13_to_p_C12_derived",  // 83,
        "O16_to_He4_C12_derived",  // 84,
        "Ne20_to_He4_O16_derived",  // 85,
        "Mg24_to_p_Na23_derived",  // 86,
        "Mg24_to_He4_Ne20_derived",  // 87,
        "Si28_to_p_Al27_derived",  // 88,
        "Si28_to_He4_Mg24_derived",  // 89,
        "S32_to_p_P31_derived",  // 90,
        "S32_to_He4_Si28_derived",  // 91,
        "Fe52_to_p_Mn51_derived",  // 92,
        "Fe52_to_He4_Cr48_derived",  // 93,
        "Co55_to_He4_Mn51_derived",  // 94,
        "Ni56_to_p_Co55_derived",  // 95,
        "Ni56_to_He4_Fe52_derived",  // 96,
        "C12_to_He4_He4_He4_derived",  // 97,
        "p_O16_to_He4_N13_derived",  // 98,
        "He4_Ne20_to_p_Na23_derived",  // 99,
        "He4_Ne20_to_C12_C12_derived",  // 100,
        "p_Na23_to_C12_C12_derived",  // 101,
        "He4_Mg24_to_p_Al27_derived",  // 102,
        "He4_Mg24_to_C12_O16_derived",  // 103,
        "p_Al27_to_C12_O16_derived",  // 104,
        "He4_Si28_to_p_P31_derived",  // 105,
        "He4_Si28_to_O16_O16_derived",  // 106,
        "p_P31_to_O16_O16_derived",  // 107,
        "p_Mn51_to_He4_Cr48_derived",  // 108,
        "p_Co55_to_He4_Fe52_derived",  // 109,
        "Mg24_to_C12_C12_derived",  // 110,
        "S32_to_O16_O16_derived",  // 111,
        "Si28_to_C12_O16_derived",  // 112,
        "Co55_to_p_Fe54_derived",  // 113,
        "Co56_to_n_Co55_derived",  // 114,
        "Co57_to_n_Co56_derived",  // 115,
        "Co57_to_p_Fe56_derived",  // 116,
        "Ni58_to_p_Co57_derived",  // 117,
        "Ni58_to_He4_Fe54_derived",  // 118,
        "p_Fe54_to_He4_Mn51_derived",  // 119,
        "He4_Fe54_to_p_Co57_derived",  // 120,
        "p_Fe56_to_n_Co56_derived",  // 121,
        "p_Co56_to_n_Ni56_derived",  // 122,
        "p_Ni58_to_He4_Co55_derived",  // 123,
        "He4_S32_to_p_Cl35_removed",  // 124,
        "Ar36_to_He4_S32_removed",  // 125,
        "Ar36_to_p_Cl35_removed",  // 126,
        "He4_Ar36_to_p_K39_removed",  // 127,
        "Ca40_to_He4_Ar36_removed",  // 128,
        "Ca40_to_p_K39_removed",  // 129,
        "He4_Ca40_to_p_Sc43_removed",  // 130,
        "Ti44_to_He4_Ca40_removed",  // 131,
        "Ti44_to_p_Sc43_removed",  // 132,
        "Cr48_to_He4_Ti44_removed",  // 133,
        "Cr48_to_p_V47_removed",  // 134,
        "p_V47_to_He4_Ti44_removed",  // 135,
        "Fe54_to_n_Fe53_removed",  // 136,
        "Fe53_to_n_Fe52_removed",  // 137,
        "Fe56_to_n_Fe55_removed",  // 138,
        "Fe55_to_n_Fe54_removed",  // 139,
        "Ni58_to_n_Ni57_removed",  // 140,
        "Ni57_to_n_Ni56_removed"  // 141,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    // A RatePair data table, 2D array of shape (NumNSERatePairs, 8).
    // Each row corresponds to a RatePair returned by _get_nse_rate_pairs().
    // Nuclei indices follow NetworkSpecies and rate indices follow NetworkRates.

    //  - Columns 1–3: 1-based index of reactant nuclei (forward rate)
    //  - Columns 4–6: 1-based index of product nuclei (forward rate)
    //  - Columns 7-8: 1-based index of non (n, H1, He4) nuclei.
    //  - Column  9:   1-based index of the forward rate
    //  - Column  10:  1-based index of the reverse rate

    constexpr int H1_index = 2;
    constexpr int N_index = 1;
    constexpr int He4_index = 3;
    constexpr int NumNSERatePairs = 36;

    inline AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, NumNSERatePairs, 1, 10, amrex::Order::C> rate_pair_data {
        2, 4, -1, 5, -1, -1, 4, 5, 1, 83,  // fr: p_C12_to_N13_reaclib, rr: N13_to_p_C12_derived
        3, 4, -1, 6, -1, -1, 4, 6, 2, 84,  // fr: He4_C12_to_O16_reaclib, rr: O16_to_He4_C12_derived
        3, 6, -1, 7, -1, -1, 6, 7, 3, 85,  // fr: He4_O16_to_Ne20_reaclib, rr: Ne20_to_He4_O16_derived
        3, 7, -1, 9, -1, -1, 7, 9, 4, 87,  // fr: He4_Ne20_to_Mg24_reaclib, rr: Mg24_to_He4_Ne20_derived
        2, 8, -1, 9, -1, -1, 8, 9, 5, 86,  // fr: p_Na23_to_Mg24_reaclib, rr: Mg24_to_p_Na23_derived
        3, 9, -1, 11, -1, -1, 9, 11, 6, 89,  // fr: He4_Mg24_to_Si28_reaclib, rr: Si28_to_He4_Mg24_derived
        2, 10, -1, 11, -1, -1, 10, 11, 7, 88,  // fr: p_Al27_to_Si28_reaclib, rr: Si28_to_p_Al27_derived
        3, 11, -1, 13, -1, -1, 11, 13, 8, 91,  // fr: He4_Si28_to_S32_reaclib, rr: S32_to_He4_Si28_derived
        2, 12, -1, 13, -1, -1, 12, 13, 9, 90,  // fr: p_P31_to_S32_reaclib, rr: S32_to_p_P31_derived
        3, 17, -1, 19, -1, -1, 17, 19, 10, 93,  // fr: He4_Cr48_to_Fe52_reaclib, rr: Fe52_to_He4_Cr48_derived
        2, 18, -1, 19, -1, -1, 18, 19, 11, 92,  // fr: p_Mn51_to_Fe52_reaclib, rr: Fe52_to_p_Mn51_derived
        3, 18, -1, 22, -1, -1, 18, 22, 12, 94,  // fr: He4_Mn51_to_Co55_reaclib, rr: Co55_to_He4_Mn51_derived
        3, 19, -1, 25, -1, -1, 19, 25, 13, 96,  // fr: He4_Fe52_to_Ni56_reaclib, rr: Ni56_to_He4_Fe52_derived
        2, 22, -1, 25, -1, -1, 22, 25, 14, 95,  // fr: p_Co55_to_Ni56_reaclib, rr: Ni56_to_p_Co55_derived
        3, 5, -1, 2, 6, -1, 5, 6, 17, 98,  // fr: He4_N13_to_p_O16_reaclib, rr: p_O16_to_He4_N13_derived
        2, 8, -1, 3, 7, -1, 8, 7, 22, 99,  // fr: p_Na23_to_He4_Ne20_reaclib, rr: He4_Ne20_to_p_Na23_derived
        2, 10, -1, 3, 9, -1, 10, 9, 23, 102,  // fr: p_Al27_to_He4_Mg24_reaclib, rr: He4_Mg24_to_p_Al27_derived
        2, 12, -1, 3, 11, -1, 12, 11, 24, 105,  // fr: p_P31_to_He4_Si28_reaclib, rr: He4_Si28_to_p_P31_derived
        3, 17, -1, 2, 18, -1, 17, 18, 25, 108,  // fr: He4_Cr48_to_p_Mn51_reaclib, rr: p_Mn51_to_He4_Cr48_derived
        3, 19, -1, 2, 22, -1, 19, 22, 26, 109,  // fr: He4_Fe52_to_p_Co55_reaclib, rr: p_Co55_to_He4_Fe52_derived
        3, 3, 3, 4, -1, -1, 4, -1, 27, 97,  // fr: He4_He4_He4_to_C12_reaclib, rr: C12_to_He4_He4_He4_derived
        2, 20, -1, 22, -1, -1, 20, 22, 31, 113,  // fr: p_Fe54_to_Co55_reaclib, rr: Co55_to_p_Fe54_derived
        3, 20, -1, 26, -1, -1, 20, 26, 32, 118,  // fr: He4_Fe54_to_Ni58_reaclib, rr: Ni58_to_He4_Fe54_derived
        2, 21, -1, 24, -1, -1, 21, 24, 33, 116,  // fr: p_Fe56_to_Co57_reaclib, rr: Co57_to_p_Fe56_derived
        1, 22, -1, 23, -1, -1, 22, 23, 34, 114,  // fr: n_Co55_to_Co56_reaclib, rr: Co56_to_n_Co55_derived
        1, 23, -1, 24, -1, -1, 23, 24, 35, 115,  // fr: n_Co56_to_Co57_reaclib, rr: Co57_to_n_Co56_derived
        2, 24, -1, 26, -1, -1, 24, 26, 36, 117,  // fr: p_Co57_to_Ni58_reaclib, rr: Ni58_to_p_Co57_derived
        3, 18, -1, 2, 20, -1, 18, 20, 37, 119,  // fr: He4_Mn51_to_p_Fe54_reaclib, rr: p_Fe54_to_He4_Mn51_derived
        3, 22, -1, 2, 26, -1, 22, 26, 38, 123,  // fr: He4_Co55_to_p_Ni58_reaclib, rr: p_Ni58_to_He4_Co55_derived
        1, 23, -1, 2, 21, -1, 23, 21, 39, 121,  // fr: n_Co56_to_p_Fe56_reaclib, rr: p_Fe56_to_n_Co56_derived
        2, 24, -1, 3, 20, -1, 24, 20, 40, 120,  // fr: p_Co57_to_He4_Fe54_reaclib, rr: He4_Fe54_to_p_Co57_derived
        1, 25, -1, 2, 23, -1, 25, 23, 41, 122,  // fr: n_Ni56_to_p_Co56_reaclib, rr: p_Co56_to_n_Ni56_derived
        13, 3, -1, 14, -1, -1, 13, 14, 66, 67,  // fr: S32_He4_to_Ar36_approx, rr: Ar36_to_S32_He4_approx
        14, 3, -1, 15, -1, -1, 14, 15, 68, 69,  // fr: Ar36_He4_to_Ca40_approx, rr: Ca40_to_Ar36_He4_approx
        15, 3, -1, 16, -1, -1, 15, 16, 70, 71,  // fr: Ca40_He4_to_Ti44_approx, rr: Ti44_to_Ca40_He4_approx
        16, 3, -1, 17, -1, -1, 16, 17, 72, 73  // fr: Ti44_He4_to_Cr48_approx, rr: Cr48_to_Ti44_He4_approx
    };
}
#endif

#endif
