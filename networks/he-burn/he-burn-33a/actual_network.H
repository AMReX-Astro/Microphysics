#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Loop.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex::literals;

void actual_network_init();

const std::string network_name = "pynucastro-cxx";

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
        k_p_C12_to_N13 = 1,
        k_He4_C12_to_O16 = 2,
        k_He4_O16_to_Ne20 = 3,
        k_He4_Ne20_to_Mg24 = 4,
        k_p_Na23_to_Mg24 = 5,
        k_He4_Mg24_to_Si28 = 6,
        k_p_Al27_to_Si28 = 7,
        k_He4_Si28_to_S32 = 8,
        k_p_P31_to_S32 = 9,
        k_He4_Cr48_to_Fe52 = 10,
        k_p_Mn51_to_Fe52 = 11,
        k_He4_Mn51_to_Co55 = 12,
        k_He4_Fe52_to_Ni56 = 13,
        k_p_Co55_to_Ni56 = 14,
        k_He4_Co55_to_Cu59 = 15,
        k_He4_Ni56_to_Zn60 = 16,
        k_p_Cu59_to_Zn60 = 17,
        k_C12_C12_to_p_Na23 = 18,
        k_C12_C12_to_He4_Ne20 = 19,
        k_He4_N13_to_p_O16 = 20,
        k_C12_O16_to_p_Al27 = 21,
        k_C12_O16_to_He4_Mg24 = 22,
        k_O16_O16_to_p_P31 = 23,
        k_O16_O16_to_He4_Si28 = 24,
        k_p_Na23_to_He4_Ne20 = 25,
        k_p_Al27_to_He4_Mg24 = 26,
        k_p_P31_to_He4_Si28 = 27,
        k_He4_Cr48_to_p_Mn51 = 28,
        k_He4_Fe52_to_p_Co55 = 29,
        k_p_Cu59_to_He4_Ni56 = 30,
        k_He4_He4_He4_to_C12 = 31,
        k_He4_N14_to_F18_removed = 32,
        k_p_O16_to_F17_removed = 33,
        k_C12_C12_to_n_Mg23_removed = 34,
        k_O16_O16_to_n_S31_removed = 35,
        k_C12_O16_to_n_Si27_removed = 36,
        k_p_Mn55_to_Fe56 = 37,
        k_n_Fe52_to_Fe53 = 38,
        k_n_Fe53_to_Fe54 = 39,
        k_He4_Fe53_to_Ni57 = 40,
        k_n_Fe54_to_Fe55 = 41,
        k_p_Fe54_to_Co55 = 42,
        k_He4_Fe54_to_Ni58 = 43,
        k_n_Fe55_to_Fe56 = 44,
        k_p_Fe55_to_Co56 = 45,
        k_p_Fe56_to_Co57 = 46,
        k_n_Co55_to_Co56 = 47,
        k_n_Co56_to_Co57 = 48,
        k_p_Co56_to_Ni57 = 49,
        k_p_Co57_to_Ni58 = 50,
        k_n_Ni56_to_Ni57 = 51,
        k_n_Ni57_to_Ni58 = 52,
        k_p_Ni58_to_Cu59 = 53,
        k_He4_Mn51_to_p_Fe54 = 54,
        k_He4_Fe53_to_p_Co56 = 55,
        k_n_Fe55_to_p_Mn55 = 56,
        k_n_Co55_to_p_Fe55 = 57,
        k_He4_Co55_to_p_Ni58 = 58,
        k_n_Co56_to_p_Fe56 = 59,
        k_p_Co57_to_He4_Fe54 = 60,
        k_n_Ni56_to_p_Co56 = 61,
        k_n_Ni56_to_He4_Fe53 = 62,
        k_n_Ni57_to_p_Co57 = 63,
        k_n_Ni57_to_He4_Fe54 = 64,
        k_n_Ni58_to_He4_Fe55 = 65,
        k_n_Cu59_to_He4_Co56 = 66,
        k_n_Zn60_to_He4_Ni57 = 67,
        k_He4_S32_to_Ar36_removed = 68,
        k_p_Cl35_to_Ar36_removed = 69,
        k_p_Cl35_to_He4_S32_removed = 70,
        k_He4_Ar36_to_Ca40_removed = 71,
        k_p_K39_to_Ca40_removed = 72,
        k_p_K39_to_He4_Ar36_removed = 73,
        k_He4_Ca40_to_Ti44_removed = 74,
        k_p_Sc43_to_Ti44_removed = 75,
        k_p_Sc43_to_He4_Ca40_removed = 76,
        k_He4_Ti44_to_Cr48_removed = 77,
        k_He4_Ti44_to_p_V47_removed = 78,
        k_p_V47_to_Cr48_removed = 79,
        k_Co55_to_Fe55 = 80,
        k_Fe55_to_Co55 = 81,
        k_Fe55_to_Mn55 = 82,
        k_Mn55_to_Fe55 = 83,
        k_Co56_to_Fe56 = 84,
        k_Co56_to_Ni56 = 85,
        k_Fe56_to_Co56 = 86,
        k_Ni56_to_Co56 = 87,
        k_Co57_to_Ni57 = 88,
        k_Ni57_to_Co57 = 89,
        k_n_to_p = 90,
        k_p_to_n = 91,
        k_S32_He4_to_Ar36_approx = 92,
        k_Ar36_to_S32_He4_approx = 93,
        k_Ar36_He4_to_Ca40_approx = 94,
        k_Ca40_to_Ar36_He4_approx = 95,
        k_Ca40_He4_to_Ti44_approx = 96,
        k_Ti44_to_Ca40_He4_approx = 97,
        k_Ti44_He4_to_Cr48_approx = 98,
        k_Cr48_to_Ti44_He4_approx = 99,
        k_He4_N14_to_Ne20_modified = 100,
        k_p_O16_to_N14_He4_modified = 101,
        k_C12_C12_to_Mg24_modified = 102,
        k_O16_O16_to_S32_modified = 103,
        k_C12_O16_to_Si28_modified = 104,
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
        k_C12_to_He4_He4_He4_derived = 122,
        k_p_O16_to_He4_N13_derived = 123,
        k_He4_Ne20_to_p_Na23_derived = 124,
        k_He4_Mg24_to_p_Al27_derived = 125,
        k_He4_Si28_to_p_P31_derived = 126,
        k_p_Mn51_to_He4_Cr48_derived = 127,
        k_p_Co55_to_He4_Fe52_derived = 128,
        k_He4_Ni56_to_p_Cu59_derived = 129,
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
        k_p_Mn55_to_n_Fe55_derived = 147,
        k_He4_Fe53_to_n_Ni56_derived = 148,
        k_p_Fe54_to_He4_Mn51_derived = 149,
        k_He4_Fe54_to_n_Ni57_derived = 150,
        k_He4_Fe54_to_p_Co57_derived = 151,
        k_p_Fe55_to_n_Co55_derived = 152,
        k_He4_Fe55_to_n_Ni58_derived = 153,
        k_p_Fe56_to_n_Co56_derived = 154,
        k_p_Co56_to_n_Ni56_derived = 155,
        k_p_Co56_to_He4_Fe53_derived = 156,
        k_He4_Co56_to_n_Cu59_derived = 157,
        k_p_Co57_to_n_Ni57_derived = 158,
        k_He4_Ni57_to_n_Zn60_derived = 159,
        k_p_Ni58_to_He4_Co55_derived = 160,
        k_He4_S32_to_p_Cl35_derived_removed = 161,
        k_Ar36_to_He4_S32_derived_removed = 162,
        k_Ar36_to_p_Cl35_derived_removed = 163,
        k_He4_Ar36_to_p_K39_derived_removed = 164,
        k_Ca40_to_He4_Ar36_derived_removed = 165,
        k_Ca40_to_p_K39_derived_removed = 166,
        k_He4_Ca40_to_p_Sc43_derived_removed = 167,
        k_Ti44_to_He4_Ca40_derived_removed = 168,
        k_Ti44_to_p_Sc43_derived_removed = 169,
        k_Cr48_to_He4_Ti44_derived_removed = 170,
        k_Cr48_to_p_V47_derived_removed = 171,
        k_p_V47_to_He4_Ti44_derived_removed = 172,
        NumRates = k_p_V47_to_He4_Ti44_derived_removed
    };

    // number of reaclib rates

    const int NrateReaclib = 147;

    // number of tabular rates

    const int NrateTabular = 12;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "p_C12_to_N13",  // 1,
        "He4_C12_to_O16",  // 2,
        "He4_O16_to_Ne20",  // 3,
        "He4_Ne20_to_Mg24",  // 4,
        "p_Na23_to_Mg24",  // 5,
        "He4_Mg24_to_Si28",  // 6,
        "p_Al27_to_Si28",  // 7,
        "He4_Si28_to_S32",  // 8,
        "p_P31_to_S32",  // 9,
        "He4_Cr48_to_Fe52",  // 10,
        "p_Mn51_to_Fe52",  // 11,
        "He4_Mn51_to_Co55",  // 12,
        "He4_Fe52_to_Ni56",  // 13,
        "p_Co55_to_Ni56",  // 14,
        "He4_Co55_to_Cu59",  // 15,
        "He4_Ni56_to_Zn60",  // 16,
        "p_Cu59_to_Zn60",  // 17,
        "C12_C12_to_p_Na23",  // 18,
        "C12_C12_to_He4_Ne20",  // 19,
        "He4_N13_to_p_O16",  // 20,
        "C12_O16_to_p_Al27",  // 21,
        "C12_O16_to_He4_Mg24",  // 22,
        "O16_O16_to_p_P31",  // 23,
        "O16_O16_to_He4_Si28",  // 24,
        "p_Na23_to_He4_Ne20",  // 25,
        "p_Al27_to_He4_Mg24",  // 26,
        "p_P31_to_He4_Si28",  // 27,
        "He4_Cr48_to_p_Mn51",  // 28,
        "He4_Fe52_to_p_Co55",  // 29,
        "p_Cu59_to_He4_Ni56",  // 30,
        "He4_He4_He4_to_C12",  // 31,
        "He4_N14_to_F18_removed",  // 32,
        "p_O16_to_F17_removed",  // 33,
        "C12_C12_to_n_Mg23_removed",  // 34,
        "O16_O16_to_n_S31_removed",  // 35,
        "C12_O16_to_n_Si27_removed",  // 36,
        "p_Mn55_to_Fe56",  // 37,
        "n_Fe52_to_Fe53",  // 38,
        "n_Fe53_to_Fe54",  // 39,
        "He4_Fe53_to_Ni57",  // 40,
        "n_Fe54_to_Fe55",  // 41,
        "p_Fe54_to_Co55",  // 42,
        "He4_Fe54_to_Ni58",  // 43,
        "n_Fe55_to_Fe56",  // 44,
        "p_Fe55_to_Co56",  // 45,
        "p_Fe56_to_Co57",  // 46,
        "n_Co55_to_Co56",  // 47,
        "n_Co56_to_Co57",  // 48,
        "p_Co56_to_Ni57",  // 49,
        "p_Co57_to_Ni58",  // 50,
        "n_Ni56_to_Ni57",  // 51,
        "n_Ni57_to_Ni58",  // 52,
        "p_Ni58_to_Cu59",  // 53,
        "He4_Mn51_to_p_Fe54",  // 54,
        "He4_Fe53_to_p_Co56",  // 55,
        "n_Fe55_to_p_Mn55",  // 56,
        "n_Co55_to_p_Fe55",  // 57,
        "He4_Co55_to_p_Ni58",  // 58,
        "n_Co56_to_p_Fe56",  // 59,
        "p_Co57_to_He4_Fe54",  // 60,
        "n_Ni56_to_p_Co56",  // 61,
        "n_Ni56_to_He4_Fe53",  // 62,
        "n_Ni57_to_p_Co57",  // 63,
        "n_Ni57_to_He4_Fe54",  // 64,
        "n_Ni58_to_He4_Fe55",  // 65,
        "n_Cu59_to_He4_Co56",  // 66,
        "n_Zn60_to_He4_Ni57",  // 67,
        "He4_S32_to_Ar36_removed",  // 68,
        "p_Cl35_to_Ar36_removed",  // 69,
        "p_Cl35_to_He4_S32_removed",  // 70,
        "He4_Ar36_to_Ca40_removed",  // 71,
        "p_K39_to_Ca40_removed",  // 72,
        "p_K39_to_He4_Ar36_removed",  // 73,
        "He4_Ca40_to_Ti44_removed",  // 74,
        "p_Sc43_to_Ti44_removed",  // 75,
        "p_Sc43_to_He4_Ca40_removed",  // 76,
        "He4_Ti44_to_Cr48_removed",  // 77,
        "He4_Ti44_to_p_V47_removed",  // 78,
        "p_V47_to_Cr48_removed",  // 79,
        "Co55_to_Fe55",  // 80,
        "Fe55_to_Co55",  // 81,
        "Fe55_to_Mn55",  // 82,
        "Mn55_to_Fe55",  // 83,
        "Co56_to_Fe56",  // 84,
        "Co56_to_Ni56",  // 85,
        "Fe56_to_Co56",  // 86,
        "Ni56_to_Co56",  // 87,
        "Co57_to_Ni57",  // 88,
        "Ni57_to_Co57",  // 89,
        "n_to_p",  // 90,
        "p_to_n",  // 91,
        "S32_He4_to_Ar36_approx",  // 92,
        "Ar36_to_S32_He4_approx",  // 93,
        "Ar36_He4_to_Ca40_approx",  // 94,
        "Ca40_to_Ar36_He4_approx",  // 95,
        "Ca40_He4_to_Ti44_approx",  // 96,
        "Ti44_to_Ca40_He4_approx",  // 97,
        "Ti44_He4_to_Cr48_approx",  // 98,
        "Cr48_to_Ti44_He4_approx",  // 99,
        "He4_N14_to_Ne20_modified",  // 100,
        "p_O16_to_N14_He4_modified",  // 101,
        "C12_C12_to_Mg24_modified",  // 102,
        "O16_O16_to_S32_modified",  // 103,
        "C12_O16_to_Si28_modified",  // 104,
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
        "C12_to_He4_He4_He4_derived",  // 122,
        "p_O16_to_He4_N13_derived",  // 123,
        "He4_Ne20_to_p_Na23_derived",  // 124,
        "He4_Mg24_to_p_Al27_derived",  // 125,
        "He4_Si28_to_p_P31_derived",  // 126,
        "p_Mn51_to_He4_Cr48_derived",  // 127,
        "p_Co55_to_He4_Fe52_derived",  // 128,
        "He4_Ni56_to_p_Cu59_derived",  // 129,
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
        "p_Mn55_to_n_Fe55_derived",  // 147,
        "He4_Fe53_to_n_Ni56_derived",  // 148,
        "p_Fe54_to_He4_Mn51_derived",  // 149,
        "He4_Fe54_to_n_Ni57_derived",  // 150,
        "He4_Fe54_to_p_Co57_derived",  // 151,
        "p_Fe55_to_n_Co55_derived",  // 152,
        "He4_Fe55_to_n_Ni58_derived",  // 153,
        "p_Fe56_to_n_Co56_derived",  // 154,
        "p_Co56_to_n_Ni56_derived",  // 155,
        "p_Co56_to_He4_Fe53_derived",  // 156,
        "He4_Co56_to_n_Cu59_derived",  // 157,
        "p_Co57_to_n_Ni57_derived",  // 158,
        "He4_Ni57_to_n_Zn60_derived",  // 159,
        "p_Ni58_to_He4_Co55_derived",  // 160,
        "He4_S32_to_p_Cl35_derived_removed",  // 161,
        "Ar36_to_He4_S32_derived_removed",  // 162,
        "Ar36_to_p_Cl35_derived_removed",  // 163,
        "He4_Ar36_to_p_K39_derived_removed",  // 164,
        "Ca40_to_He4_Ar36_derived_removed",  // 165,
        "Ca40_to_p_K39_derived_removed",  // 166,
        "He4_Ca40_to_p_Sc43_derived_removed",  // 167,
        "Ti44_to_He4_Ca40_derived_removed",  // 168,
        "Ti44_to_p_Sc43_derived_removed",  // 169,
        "Cr48_to_He4_Ti44_derived_removed",  // 170,
        "Cr48_to_p_V47_derived_removed",  // 171,
        "p_V47_to_He4_Ti44_derived_removed"  // 172,
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
