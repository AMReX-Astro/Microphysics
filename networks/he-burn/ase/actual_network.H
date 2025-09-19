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
        k_p_C12_to_N13 = 1,
        k_He4_C12_to_O16 = 2,
        k_He4_O16_to_Ne20 = 3,
        k_He4_Ne20_to_Mg24 = 4,
        k_p_Na23_to_Mg24 = 5,
        k_He4_Mg24_to_Si28 = 6,
        k_p_Al27_to_Si28 = 7,
        k_He4_Si28_to_S32 = 8,
        k_p_P31_to_S32 = 9,
        k_C12_C12_to_p_Na23 = 10,
        k_C12_C12_to_He4_Ne20 = 11,
        k_He4_N13_to_p_O16 = 12,
        k_C12_O16_to_p_Al27 = 13,
        k_C12_O16_to_He4_Mg24 = 14,
        k_O16_O16_to_p_P31 = 15,
        k_O16_O16_to_He4_Si28 = 16,
        k_p_Na23_to_He4_Ne20 = 17,
        k_p_Al27_to_He4_Mg24 = 18,
        k_p_P31_to_He4_Si28 = 19,
        k_He4_He4_He4_to_C12 = 20,
        k_C12_C12_to_n_Mg23_removed = 21,
        k_O16_O16_to_n_S31_removed = 22,
        k_C12_O16_to_n_Si27_removed = 23,
        k_He4_S32_to_Ar36_removed = 24,
        k_p_Cl35_to_Ar36_removed = 25,
        k_p_Cl35_to_He4_S32_removed = 26,
        k_He4_Ar36_to_Ca40_removed = 27,
        k_p_K39_to_Ca40_removed = 28,
        k_p_K39_to_He4_Ar36_removed = 29,
        k_He4_Ca40_to_Ti44_removed = 30,
        k_p_Sc43_to_Ti44_removed = 31,
        k_p_Sc43_to_He4_Ca40_removed = 32,
        k_He4_Ti44_to_Cr48_removed = 33,
        k_He4_Ti44_to_p_V47_removed = 34,
        k_p_V47_to_Cr48_removed = 35,
        k_He4_Cr48_to_Fe52_removed = 36,
        k_He4_Cr48_to_p_Mn51_removed = 37,
        k_p_Mn51_to_Fe52_removed = 38,
        k_He4_Fe52_to_Ni56_removed = 39,
        k_He4_Fe52_to_p_Co55_removed = 40,
        k_p_Co55_to_Ni56_removed = 41,
        k_S32_He4_to_Ar36_approx = 42,
        k_Ar36_to_S32_He4_approx = 43,
        k_Ar36_He4_to_Ca40_approx = 44,
        k_Ca40_to_Ar36_He4_approx = 45,
        k_Ca40_He4_to_Ti44_approx = 46,
        k_Ti44_to_Ca40_He4_approx = 47,
        k_Ti44_He4_to_Cr48_approx = 48,
        k_Cr48_to_Ti44_He4_approx = 49,
        k_Cr48_He4_to_Fe52_approx = 50,
        k_Fe52_to_Cr48_He4_approx = 51,
        k_Fe52_He4_to_Ni56_approx = 52,
        k_Ni56_to_Fe52_He4_approx = 53,
        k_C12_C12_to_Mg24_modified = 54,
        k_O16_O16_to_S32_modified = 55,
        k_C12_O16_to_Si28_modified = 56,
        k_N13_to_p_C12_derived = 57,
        k_O16_to_He4_C12_derived = 58,
        k_Ne20_to_He4_O16_derived = 59,
        k_Mg24_to_p_Na23_derived = 60,
        k_Mg24_to_He4_Ne20_derived = 61,
        k_Si28_to_p_Al27_derived = 62,
        k_Si28_to_He4_Mg24_derived = 63,
        k_S32_to_p_P31_derived = 64,
        k_S32_to_He4_Si28_derived = 65,
        k_C12_to_He4_He4_He4_derived = 66,
        k_p_O16_to_He4_N13_derived = 67,
        k_He4_Ne20_to_p_Na23_derived = 68,
        k_He4_Ne20_to_C12_C12_derived = 69,
        k_p_Na23_to_C12_C12_derived = 70,
        k_He4_Mg24_to_p_Al27_derived = 71,
        k_He4_Mg24_to_C12_O16_derived = 72,
        k_p_Al27_to_C12_O16_derived = 73,
        k_He4_Si28_to_p_P31_derived = 74,
        k_He4_Si28_to_O16_O16_derived = 75,
        k_p_P31_to_O16_O16_derived = 76,
        k_Mg24_to_C12_C12_derived = 77,
        k_S32_to_O16_O16_derived = 78,
        k_Si28_to_C12_O16_derived = 79,
        k_He4_S32_to_p_Cl35_derived_removed = 80,
        k_Ar36_to_He4_S32_derived_removed = 81,
        k_Ar36_to_p_Cl35_derived_removed = 82,
        k_He4_Ar36_to_p_K39_derived_removed = 83,
        k_Ca40_to_He4_Ar36_derived_removed = 84,
        k_Ca40_to_p_K39_derived_removed = 85,
        k_He4_Ca40_to_p_Sc43_derived_removed = 86,
        k_Ti44_to_He4_Ca40_derived_removed = 87,
        k_Ti44_to_p_Sc43_derived_removed = 88,
        k_Cr48_to_He4_Ti44_derived_removed = 89,
        k_Cr48_to_p_V47_derived_removed = 90,
        k_p_V47_to_He4_Ti44_derived_removed = 91,
        k_Fe52_to_He4_Cr48_derived_removed = 92,
        k_Fe52_to_p_Mn51_derived_removed = 93,
        k_p_Mn51_to_He4_Cr48_derived_removed = 94,
        k_Ni56_to_He4_Fe52_derived_removed = 95,
        k_Ni56_to_p_Co55_derived_removed = 96,
        k_p_Co55_to_He4_Fe52_derived_removed = 97,
        NumRates = k_p_Co55_to_He4_Fe52_derived_removed
    };

    // number of reaclib rates

    const int NrateReaclib = 82;

    // number of tabular rates

    const int NrateTabular = 0;

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
        "C12_C12_to_p_Na23",  // 10,
        "C12_C12_to_He4_Ne20",  // 11,
        "He4_N13_to_p_O16",  // 12,
        "C12_O16_to_p_Al27",  // 13,
        "C12_O16_to_He4_Mg24",  // 14,
        "O16_O16_to_p_P31",  // 15,
        "O16_O16_to_He4_Si28",  // 16,
        "p_Na23_to_He4_Ne20",  // 17,
        "p_Al27_to_He4_Mg24",  // 18,
        "p_P31_to_He4_Si28",  // 19,
        "He4_He4_He4_to_C12",  // 20,
        "C12_C12_to_n_Mg23_removed",  // 21,
        "O16_O16_to_n_S31_removed",  // 22,
        "C12_O16_to_n_Si27_removed",  // 23,
        "He4_S32_to_Ar36_removed",  // 24,
        "p_Cl35_to_Ar36_removed",  // 25,
        "p_Cl35_to_He4_S32_removed",  // 26,
        "He4_Ar36_to_Ca40_removed",  // 27,
        "p_K39_to_Ca40_removed",  // 28,
        "p_K39_to_He4_Ar36_removed",  // 29,
        "He4_Ca40_to_Ti44_removed",  // 30,
        "p_Sc43_to_Ti44_removed",  // 31,
        "p_Sc43_to_He4_Ca40_removed",  // 32,
        "He4_Ti44_to_Cr48_removed",  // 33,
        "He4_Ti44_to_p_V47_removed",  // 34,
        "p_V47_to_Cr48_removed",  // 35,
        "He4_Cr48_to_Fe52_removed",  // 36,
        "He4_Cr48_to_p_Mn51_removed",  // 37,
        "p_Mn51_to_Fe52_removed",  // 38,
        "He4_Fe52_to_Ni56_removed",  // 39,
        "He4_Fe52_to_p_Co55_removed",  // 40,
        "p_Co55_to_Ni56_removed",  // 41,
        "S32_He4_to_Ar36_approx",  // 42,
        "Ar36_to_S32_He4_approx",  // 43,
        "Ar36_He4_to_Ca40_approx",  // 44,
        "Ca40_to_Ar36_He4_approx",  // 45,
        "Ca40_He4_to_Ti44_approx",  // 46,
        "Ti44_to_Ca40_He4_approx",  // 47,
        "Ti44_He4_to_Cr48_approx",  // 48,
        "Cr48_to_Ti44_He4_approx",  // 49,
        "Cr48_He4_to_Fe52_approx",  // 50,
        "Fe52_to_Cr48_He4_approx",  // 51,
        "Fe52_He4_to_Ni56_approx",  // 52,
        "Ni56_to_Fe52_He4_approx",  // 53,
        "C12_C12_to_Mg24_modified",  // 54,
        "O16_O16_to_S32_modified",  // 55,
        "C12_O16_to_Si28_modified",  // 56,
        "N13_to_p_C12_derived",  // 57,
        "O16_to_He4_C12_derived",  // 58,
        "Ne20_to_He4_O16_derived",  // 59,
        "Mg24_to_p_Na23_derived",  // 60,
        "Mg24_to_He4_Ne20_derived",  // 61,
        "Si28_to_p_Al27_derived",  // 62,
        "Si28_to_He4_Mg24_derived",  // 63,
        "S32_to_p_P31_derived",  // 64,
        "S32_to_He4_Si28_derived",  // 65,
        "C12_to_He4_He4_He4_derived",  // 66,
        "p_O16_to_He4_N13_derived",  // 67,
        "He4_Ne20_to_p_Na23_derived",  // 68,
        "He4_Ne20_to_C12_C12_derived",  // 69,
        "p_Na23_to_C12_C12_derived",  // 70,
        "He4_Mg24_to_p_Al27_derived",  // 71,
        "He4_Mg24_to_C12_O16_derived",  // 72,
        "p_Al27_to_C12_O16_derived",  // 73,
        "He4_Si28_to_p_P31_derived",  // 74,
        "He4_Si28_to_O16_O16_derived",  // 75,
        "p_P31_to_O16_O16_derived",  // 76,
        "Mg24_to_C12_C12_derived",  // 77,
        "S32_to_O16_O16_derived",  // 78,
        "Si28_to_C12_O16_derived",  // 79,
        "He4_S32_to_p_Cl35_derived_removed",  // 80,
        "Ar36_to_He4_S32_derived_removed",  // 81,
        "Ar36_to_p_Cl35_derived_removed",  // 82,
        "He4_Ar36_to_p_K39_derived_removed",  // 83,
        "Ca40_to_He4_Ar36_derived_removed",  // 84,
        "Ca40_to_p_K39_derived_removed",  // 85,
        "He4_Ca40_to_p_Sc43_derived_removed",  // 86,
        "Ti44_to_He4_Ca40_derived_removed",  // 87,
        "Ti44_to_p_Sc43_derived_removed",  // 88,
        "Cr48_to_He4_Ti44_derived_removed",  // 89,
        "Cr48_to_p_V47_derived_removed",  // 90,
        "p_V47_to_He4_Ti44_derived_removed",  // 91,
        "Fe52_to_He4_Cr48_derived_removed",  // 92,
        "Fe52_to_p_Mn51_derived_removed",  // 93,
        "p_Mn51_to_He4_Cr48_derived_removed",  // 94,
        "Ni56_to_He4_Fe52_derived_removed",  // 95,
        "Ni56_to_p_Co55_derived_removed",  // 96,
        "p_Co55_to_He4_Fe52_derived_removed"  // 97,
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

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, amrex::Order::C> rate_indices;
}
#endif

#endif
