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
        k_He4_N14_to_F18_removed = 21,
        k_p_O16_to_F17_removed = 22,
        k_C12_C12_to_n_Mg23_removed = 23,
        k_O16_O16_to_n_S31_removed = 24,
        k_C12_O16_to_n_Si27_removed = 25,
        k_He4_S32_to_Ar36_removed = 26,
        k_p_Cl35_to_Ar36_removed = 27,
        k_p_Cl35_to_He4_S32_removed = 28,
        k_He4_Ar36_to_Ca40_removed = 29,
        k_p_K39_to_Ca40_removed = 30,
        k_p_K39_to_He4_Ar36_removed = 31,
        k_He4_Ca40_to_Ti44_removed = 32,
        k_p_Sc43_to_Ti44_removed = 33,
        k_p_Sc43_to_He4_Ca40_removed = 34,
        k_He4_Ti44_to_Cr48_removed = 35,
        k_He4_Ti44_to_p_V47_removed = 36,
        k_p_V47_to_Cr48_removed = 37,
        k_He4_Cr48_to_Fe52_removed = 38,
        k_He4_Cr48_to_p_Mn51_removed = 39,
        k_p_Mn51_to_Fe52_removed = 40,
        k_He4_Fe52_to_Ni56_removed = 41,
        k_He4_Fe52_to_p_Co55_removed = 42,
        k_p_Co55_to_Ni56_removed = 43,
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
        k_He4_N14_to_Ne20_modified = 56,
        k_p_O16_to_N14_He4_modified = 57,
        k_C12_C12_to_Mg24_modified = 58,
        k_O16_O16_to_S32_modified = 59,
        k_C12_O16_to_Si28_modified = 60,
        k_N13_to_p_C12_derived = 61,
        k_O16_to_He4_C12_derived = 62,
        k_Ne20_to_He4_O16_derived = 63,
        k_Mg24_to_p_Na23_derived = 64,
        k_Mg24_to_He4_Ne20_derived = 65,
        k_Si28_to_p_Al27_derived = 66,
        k_Si28_to_He4_Mg24_derived = 67,
        k_S32_to_p_P31_derived = 68,
        k_S32_to_He4_Si28_derived = 69,
        k_C12_to_He4_He4_He4_derived = 70,
        k_p_O16_to_He4_N13_derived = 71,
        k_He4_Ne20_to_p_Na23_derived = 72,
        k_He4_Mg24_to_p_Al27_derived = 73,
        k_He4_Si28_to_p_P31_derived = 74,
        k_He4_S32_to_p_Cl35_derived_removed = 75,
        k_Ar36_to_He4_S32_derived_removed = 76,
        k_Ar36_to_p_Cl35_derived_removed = 77,
        k_He4_Ar36_to_p_K39_derived_removed = 78,
        k_Ca40_to_He4_Ar36_derived_removed = 79,
        k_Ca40_to_p_K39_derived_removed = 80,
        k_He4_Ca40_to_p_Sc43_derived_removed = 81,
        k_Ti44_to_He4_Ca40_derived_removed = 82,
        k_Ti44_to_p_Sc43_derived_removed = 83,
        k_Cr48_to_He4_Ti44_derived_removed = 84,
        k_Cr48_to_p_V47_derived_removed = 85,
        k_p_V47_to_He4_Ti44_derived_removed = 86,
        k_Fe52_to_He4_Cr48_derived_removed = 87,
        k_Fe52_to_p_Mn51_derived_removed = 88,
        k_p_Mn51_to_He4_Cr48_derived_removed = 89,
        k_Ni56_to_He4_Fe52_derived_removed = 90,
        k_Ni56_to_p_Co55_derived_removed = 91,
        k_p_Co55_to_He4_Fe52_derived_removed = 92,
        NumRates = k_p_Co55_to_He4_Fe52_derived_removed
    };

    // number of reaclib rates

    const int NrateReaclib = 75;

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
        "He4_N14_to_F18_removed",  // 21,
        "p_O16_to_F17_removed",  // 22,
        "C12_C12_to_n_Mg23_removed",  // 23,
        "O16_O16_to_n_S31_removed",  // 24,
        "C12_O16_to_n_Si27_removed",  // 25,
        "He4_S32_to_Ar36_removed",  // 26,
        "p_Cl35_to_Ar36_removed",  // 27,
        "p_Cl35_to_He4_S32_removed",  // 28,
        "He4_Ar36_to_Ca40_removed",  // 29,
        "p_K39_to_Ca40_removed",  // 30,
        "p_K39_to_He4_Ar36_removed",  // 31,
        "He4_Ca40_to_Ti44_removed",  // 32,
        "p_Sc43_to_Ti44_removed",  // 33,
        "p_Sc43_to_He4_Ca40_removed",  // 34,
        "He4_Ti44_to_Cr48_removed",  // 35,
        "He4_Ti44_to_p_V47_removed",  // 36,
        "p_V47_to_Cr48_removed",  // 37,
        "He4_Cr48_to_Fe52_removed",  // 38,
        "He4_Cr48_to_p_Mn51_removed",  // 39,
        "p_Mn51_to_Fe52_removed",  // 40,
        "He4_Fe52_to_Ni56_removed",  // 41,
        "He4_Fe52_to_p_Co55_removed",  // 42,
        "p_Co55_to_Ni56_removed",  // 43,
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
        "He4_N14_to_Ne20_modified",  // 56,
        "p_O16_to_N14_He4_modified",  // 57,
        "C12_C12_to_Mg24_modified",  // 58,
        "O16_O16_to_S32_modified",  // 59,
        "C12_O16_to_Si28_modified",  // 60,
        "N13_to_p_C12_derived",  // 61,
        "O16_to_He4_C12_derived",  // 62,
        "Ne20_to_He4_O16_derived",  // 63,
        "Mg24_to_p_Na23_derived",  // 64,
        "Mg24_to_He4_Ne20_derived",  // 65,
        "Si28_to_p_Al27_derived",  // 66,
        "Si28_to_He4_Mg24_derived",  // 67,
        "S32_to_p_P31_derived",  // 68,
        "S32_to_He4_Si28_derived",  // 69,
        "C12_to_He4_He4_He4_derived",  // 70,
        "p_O16_to_He4_N13_derived",  // 71,
        "He4_Ne20_to_p_Na23_derived",  // 72,
        "He4_Mg24_to_p_Al27_derived",  // 73,
        "He4_Si28_to_p_P31_derived",  // 74,
        "He4_S32_to_p_Cl35_derived_removed",  // 75,
        "Ar36_to_He4_S32_derived_removed",  // 76,
        "Ar36_to_p_Cl35_derived_removed",  // 77,
        "He4_Ar36_to_p_K39_derived_removed",  // 78,
        "Ca40_to_He4_Ar36_derived_removed",  // 79,
        "Ca40_to_p_K39_derived_removed",  // 80,
        "He4_Ca40_to_p_Sc43_derived_removed",  // 81,
        "Ti44_to_He4_Ca40_derived_removed",  // 82,
        "Ti44_to_p_Sc43_derived_removed",  // 83,
        "Cr48_to_He4_Ti44_derived_removed",  // 84,
        "Cr48_to_p_V47_derived_removed",  // 85,
        "p_V47_to_He4_Ti44_derived_removed",  // 86,
        "Fe52_to_He4_Cr48_derived_removed",  // 87,
        "Fe52_to_p_Mn51_derived_removed",  // 88,
        "p_Mn51_to_He4_Cr48_derived_removed",  // 89,
        "Ni56_to_He4_Fe52_derived_removed",  // 90,
        "Ni56_to_p_Co55_derived_removed",  // 91,
        "p_Co55_to_He4_Fe52_derived_removed"  // 92,
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
