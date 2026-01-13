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

        if constexpr (spec == H1) {
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
        else if constexpr (spec == F17) {
            return 128.21958437599824_rt;
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
        else if constexpr (spec == Fe52) {
            return 447.6996182639923_rt;
        }
        else if constexpr (spec == Ni56) {
            return 483.9956965919919_rt;
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

        if constexpr (spec == H1) {
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
        else if constexpr (spec == F17) {
            return 2.8232643353472436e-23_rt;
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
        else if constexpr (spec == Fe52) {
            return 8.626187166893794e-23_rt;
        }
        else if constexpr (spec == Ni56) {
            return 9.289408870379396e-23_rt;
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
        k_C12_C12_to_p_Na23_reaclib = 10,
        k_C12_C12_to_He4_Ne20_reaclib = 11,
        k_He4_N13_to_p_O16_reaclib = 12,
        k_C12_O16_to_p_Al27_reaclib = 13,
        k_C12_O16_to_He4_Mg24_reaclib = 14,
        k_O16_O16_to_p_P31_reaclib = 15,
        k_O16_O16_to_He4_Si28_reaclib = 16,
        k_p_Ne20_to_He4_F17_reaclib = 17,
        k_p_Na23_to_He4_Ne20_reaclib = 18,
        k_p_Al27_to_He4_Mg24_reaclib = 19,
        k_p_P31_to_He4_Si28_reaclib = 20,
        k_He4_He4_He4_to_C12_reaclib = 21,
        k_C12_C12_to_n_Mg23_removed = 22,
        k_O16_O16_to_n_S31_removed = 23,
        k_C12_O16_to_n_Si27_removed = 24,
        k_He4_S32_to_Ar36_removed = 25,
        k_p_Cl35_to_Ar36_removed = 26,
        k_p_Cl35_to_He4_S32_removed = 27,
        k_He4_Ar36_to_Ca40_removed = 28,
        k_p_K39_to_Ca40_removed = 29,
        k_p_K39_to_He4_Ar36_removed = 30,
        k_He4_Ca40_to_Ti44_removed = 31,
        k_p_Sc43_to_Ti44_removed = 32,
        k_p_Sc43_to_He4_Ca40_removed = 33,
        k_He4_Ti44_to_Cr48_removed = 34,
        k_He4_Ti44_to_p_V47_removed = 35,
        k_p_V47_to_Cr48_removed = 36,
        k_He4_Cr48_to_Fe52_removed = 37,
        k_He4_Cr48_to_p_Mn51_removed = 38,
        k_p_Mn51_to_Fe52_removed = 39,
        k_He4_Fe52_to_Ni56_removed = 40,
        k_He4_Fe52_to_p_Co55_removed = 41,
        k_p_Co55_to_Ni56_removed = 42,
        k_O16_p_to_F17_iliadis = 43,
        k_S32_He4_to_Ar36_approx = 44,
        k_Ar36_to_S32_He4_approx = 45,
        k_Ar36_He4_to_Ca40_approx = 46,
        k_Ca40_to_Ar36_He4_approx = 47,
        k_Ca40_He4_to_Ti44_approx = 48,
        k_Ti44_to_Ca40_He4_approx = 49,
        k_Ti44_He4_to_Cr48_approx = 50,
        k_Cr48_to_Ti44_He4_approx = 51,
        k_Cr48_He4_to_Fe52_approx = 52,
        k_Fe52_to_Cr48_He4_approx = 53,
        k_Fe52_He4_to_Ni56_approx = 54,
        k_Ni56_to_Fe52_He4_approx = 55,
        k_C12_C12_to_Mg24_modified = 56,
        k_O16_O16_to_S32_modified = 57,
        k_C12_O16_to_Si28_modified = 58,
        k_N13_to_p_C12_derived = 59,
        k_O16_to_He4_C12_derived = 60,
        k_Ne20_to_He4_O16_derived = 61,
        k_Mg24_to_p_Na23_derived = 62,
        k_Mg24_to_He4_Ne20_derived = 63,
        k_Si28_to_p_Al27_derived = 64,
        k_Si28_to_He4_Mg24_derived = 65,
        k_S32_to_p_P31_derived = 66,
        k_S32_to_He4_Si28_derived = 67,
        k_C12_to_He4_He4_He4_derived = 68,
        k_p_O16_to_He4_N13_derived = 69,
        k_He4_F17_to_p_Ne20_derived = 70,
        k_He4_Ne20_to_p_Na23_derived = 71,
        k_He4_Ne20_to_C12_C12_derived = 72,
        k_p_Na23_to_C12_C12_derived = 73,
        k_He4_Mg24_to_p_Al27_derived = 74,
        k_He4_Mg24_to_C12_O16_derived = 75,
        k_p_Al27_to_C12_O16_derived = 76,
        k_He4_Si28_to_p_P31_derived = 77,
        k_He4_Si28_to_O16_O16_derived = 78,
        k_p_P31_to_O16_O16_derived = 79,
        k_Mg24_to_C12_C12_derived = 80,
        k_S32_to_O16_O16_derived = 81,
        k_Si28_to_C12_O16_derived = 82,
        k_F17_to_O16_p_derived = 83,
        k_He4_S32_to_p_Cl35_removed = 84,
        k_Ar36_to_He4_S32_removed = 85,
        k_Ar36_to_p_Cl35_removed = 86,
        k_He4_Ar36_to_p_K39_removed = 87,
        k_Ca40_to_He4_Ar36_removed = 88,
        k_Ca40_to_p_K39_removed = 89,
        k_He4_Ca40_to_p_Sc43_removed = 90,
        k_Ti44_to_He4_Ca40_removed = 91,
        k_Ti44_to_p_Sc43_removed = 92,
        k_Cr48_to_He4_Ti44_removed = 93,
        k_Cr48_to_p_V47_removed = 94,
        k_p_V47_to_He4_Ti44_removed = 95,
        k_Fe52_to_He4_Cr48_removed = 96,
        k_Fe52_to_p_Mn51_removed = 97,
        k_p_Mn51_to_He4_Cr48_removed = 98,
        k_Ni56_to_He4_Fe52_removed = 99,
        k_Ni56_to_p_Co55_removed = 100,
        k_p_Co55_to_He4_Fe52_removed = 101,
        NumRates = k_p_Co55_to_He4_Fe52_removed
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
        "C12_C12_to_p_Na23_reaclib",  // 10,
        "C12_C12_to_He4_Ne20_reaclib",  // 11,
        "He4_N13_to_p_O16_reaclib",  // 12,
        "C12_O16_to_p_Al27_reaclib",  // 13,
        "C12_O16_to_He4_Mg24_reaclib",  // 14,
        "O16_O16_to_p_P31_reaclib",  // 15,
        "O16_O16_to_He4_Si28_reaclib",  // 16,
        "p_Ne20_to_He4_F17_reaclib",  // 17,
        "p_Na23_to_He4_Ne20_reaclib",  // 18,
        "p_Al27_to_He4_Mg24_reaclib",  // 19,
        "p_P31_to_He4_Si28_reaclib",  // 20,
        "He4_He4_He4_to_C12_reaclib",  // 21,
        "C12_C12_to_n_Mg23_removed",  // 22,
        "O16_O16_to_n_S31_removed",  // 23,
        "C12_O16_to_n_Si27_removed",  // 24,
        "He4_S32_to_Ar36_removed",  // 25,
        "p_Cl35_to_Ar36_removed",  // 26,
        "p_Cl35_to_He4_S32_removed",  // 27,
        "He4_Ar36_to_Ca40_removed",  // 28,
        "p_K39_to_Ca40_removed",  // 29,
        "p_K39_to_He4_Ar36_removed",  // 30,
        "He4_Ca40_to_Ti44_removed",  // 31,
        "p_Sc43_to_Ti44_removed",  // 32,
        "p_Sc43_to_He4_Ca40_removed",  // 33,
        "He4_Ti44_to_Cr48_removed",  // 34,
        "He4_Ti44_to_p_V47_removed",  // 35,
        "p_V47_to_Cr48_removed",  // 36,
        "He4_Cr48_to_Fe52_removed",  // 37,
        "He4_Cr48_to_p_Mn51_removed",  // 38,
        "p_Mn51_to_Fe52_removed",  // 39,
        "He4_Fe52_to_Ni56_removed",  // 40,
        "He4_Fe52_to_p_Co55_removed",  // 41,
        "p_Co55_to_Ni56_removed",  // 42,
        "O16_p_to_F17_iliadis",  // 43,
        "S32_He4_to_Ar36_approx",  // 44,
        "Ar36_to_S32_He4_approx",  // 45,
        "Ar36_He4_to_Ca40_approx",  // 46,
        "Ca40_to_Ar36_He4_approx",  // 47,
        "Ca40_He4_to_Ti44_approx",  // 48,
        "Ti44_to_Ca40_He4_approx",  // 49,
        "Ti44_He4_to_Cr48_approx",  // 50,
        "Cr48_to_Ti44_He4_approx",  // 51,
        "Cr48_He4_to_Fe52_approx",  // 52,
        "Fe52_to_Cr48_He4_approx",  // 53,
        "Fe52_He4_to_Ni56_approx",  // 54,
        "Ni56_to_Fe52_He4_approx",  // 55,
        "C12_C12_to_Mg24_modified",  // 56,
        "O16_O16_to_S32_modified",  // 57,
        "C12_O16_to_Si28_modified",  // 58,
        "N13_to_p_C12_derived",  // 59,
        "O16_to_He4_C12_derived",  // 60,
        "Ne20_to_He4_O16_derived",  // 61,
        "Mg24_to_p_Na23_derived",  // 62,
        "Mg24_to_He4_Ne20_derived",  // 63,
        "Si28_to_p_Al27_derived",  // 64,
        "Si28_to_He4_Mg24_derived",  // 65,
        "S32_to_p_P31_derived",  // 66,
        "S32_to_He4_Si28_derived",  // 67,
        "C12_to_He4_He4_He4_derived",  // 68,
        "p_O16_to_He4_N13_derived",  // 69,
        "He4_F17_to_p_Ne20_derived",  // 70,
        "He4_Ne20_to_p_Na23_derived",  // 71,
        "He4_Ne20_to_C12_C12_derived",  // 72,
        "p_Na23_to_C12_C12_derived",  // 73,
        "He4_Mg24_to_p_Al27_derived",  // 74,
        "He4_Mg24_to_C12_O16_derived",  // 75,
        "p_Al27_to_C12_O16_derived",  // 76,
        "He4_Si28_to_p_P31_derived",  // 77,
        "He4_Si28_to_O16_O16_derived",  // 78,
        "p_P31_to_O16_O16_derived",  // 79,
        "Mg24_to_C12_C12_derived",  // 80,
        "S32_to_O16_O16_derived",  // 81,
        "Si28_to_C12_O16_derived",  // 82,
        "F17_to_O16_p_derived",  // 83,
        "He4_S32_to_p_Cl35_removed",  // 84,
        "Ar36_to_He4_S32_removed",  // 85,
        "Ar36_to_p_Cl35_removed",  // 86,
        "He4_Ar36_to_p_K39_removed",  // 87,
        "Ca40_to_He4_Ar36_removed",  // 88,
        "Ca40_to_p_K39_removed",  // 89,
        "He4_Ca40_to_p_Sc43_removed",  // 90,
        "Ti44_to_He4_Ca40_removed",  // 91,
        "Ti44_to_p_Sc43_removed",  // 92,
        "Cr48_to_He4_Ti44_removed",  // 93,
        "Cr48_to_p_V47_removed",  // 94,
        "p_V47_to_He4_Ti44_removed",  // 95,
        "Fe52_to_He4_Cr48_removed",  // 96,
        "Fe52_to_p_Mn51_removed",  // 97,
        "p_Mn51_to_He4_Cr48_removed",  // 98,
        "Ni56_to_He4_Fe52_removed",  // 99,
        "Ni56_to_p_Co55_removed",  // 100,
        "p_Co55_to_He4_Fe52_removed"  // 101,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int H1_index = 0;
    constexpr int N_index = -1;
    constexpr int He4_index = 1;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, Rates::NumRates, 1, 7, amrex::Order::C> rate_indices;
}
#endif

#endif
