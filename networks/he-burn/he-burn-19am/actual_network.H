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
        k_p_Na23_to_He4_Ne20_reaclib = 17,
        k_p_Al27_to_He4_Mg24_reaclib = 18,
        k_p_P31_to_He4_Si28_reaclib = 19,
        k_He4_He4_He4_to_C12_reaclib = 20,
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
        k_He4_S32_to_p_Cl35_removed = 75,
        k_Ar36_to_He4_S32_removed = 76,
        k_Ar36_to_p_Cl35_removed = 77,
        k_He4_Ar36_to_p_K39_removed = 78,
        k_Ca40_to_He4_Ar36_removed = 79,
        k_Ca40_to_p_K39_removed = 80,
        k_He4_Ca40_to_p_Sc43_removed = 81,
        k_Ti44_to_He4_Ca40_removed = 82,
        k_Ti44_to_p_Sc43_removed = 83,
        k_Cr48_to_He4_Ti44_removed = 84,
        k_Cr48_to_p_V47_removed = 85,
        k_p_V47_to_He4_Ti44_removed = 86,
        k_Fe52_to_He4_Cr48_removed = 87,
        k_Fe52_to_p_Mn51_removed = 88,
        k_p_Mn51_to_He4_Cr48_removed = 89,
        k_Ni56_to_He4_Fe52_removed = 90,
        k_Ni56_to_p_Co55_removed = 91,
        k_p_Co55_to_He4_Fe52_removed = 92,
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
        "p_Na23_to_He4_Ne20_reaclib",  // 17,
        "p_Al27_to_He4_Mg24_reaclib",  // 18,
        "p_P31_to_He4_Si28_reaclib",  // 19,
        "He4_He4_He4_to_C12_reaclib",  // 20,
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
        "He4_S32_to_p_Cl35_removed",  // 75,
        "Ar36_to_He4_S32_removed",  // 76,
        "Ar36_to_p_Cl35_removed",  // 77,
        "He4_Ar36_to_p_K39_removed",  // 78,
        "Ca40_to_He4_Ar36_removed",  // 79,
        "Ca40_to_p_K39_removed",  // 80,
        "He4_Ca40_to_p_Sc43_removed",  // 81,
        "Ti44_to_He4_Ca40_removed",  // 82,
        "Ti44_to_p_Sc43_removed",  // 83,
        "Cr48_to_He4_Ti44_removed",  // 84,
        "Cr48_to_p_V47_removed",  // 85,
        "p_V47_to_He4_Ti44_removed",  // 86,
        "Fe52_to_He4_Cr48_removed",  // 87,
        "Fe52_to_p_Mn51_removed",  // 88,
        "p_Mn51_to_He4_Cr48_removed",  // 89,
        "Ni56_to_He4_Fe52_removed",  // 90,
        "Ni56_to_p_Co55_removed",  // 91,
        "p_Co55_to_He4_Fe52_removed"  // 92,
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

    constexpr int H1_index = 1;
    constexpr int N_index = -1;
    constexpr int He4_index = 2;
    constexpr int NumNSERatePairs = 20;

    inline AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, NumNSERatePairs, 1, 10, amrex::Order::C> rate_pair_data {
        1, 3, -1, 4, -1, -1, 3, 4, 1, 61,  // fr: p_C12_to_N13_reaclib, rr: N13_to_p_C12_derived
        2, 3, -1, 6, -1, -1, 3, 6, 2, 62,  // fr: He4_C12_to_O16_reaclib, rr: O16_to_He4_C12_derived
        2, 6, -1, 7, -1, -1, 6, 7, 3, 63,  // fr: He4_O16_to_Ne20_reaclib, rr: Ne20_to_He4_O16_derived
        2, 7, -1, 9, -1, -1, 7, 9, 4, 65,  // fr: He4_Ne20_to_Mg24_reaclib, rr: Mg24_to_He4_Ne20_derived
        1, 8, -1, 9, -1, -1, 8, 9, 5, 64,  // fr: p_Na23_to_Mg24_reaclib, rr: Mg24_to_p_Na23_derived
        2, 9, -1, 11, -1, -1, 9, 11, 6, 67,  // fr: He4_Mg24_to_Si28_reaclib, rr: Si28_to_He4_Mg24_derived
        1, 10, -1, 11, -1, -1, 10, 11, 7, 66,  // fr: p_Al27_to_Si28_reaclib, rr: Si28_to_p_Al27_derived
        2, 11, -1, 13, -1, -1, 11, 13, 8, 69,  // fr: He4_Si28_to_S32_reaclib, rr: S32_to_He4_Si28_derived
        1, 12, -1, 13, -1, -1, 12, 13, 9, 68,  // fr: p_P31_to_S32_reaclib, rr: S32_to_p_P31_derived
        2, 4, -1, 1, 6, -1, 4, 6, 12, 71,  // fr: He4_N13_to_p_O16_reaclib, rr: p_O16_to_He4_N13_derived
        1, 8, -1, 2, 7, -1, 8, 7, 17, 72,  // fr: p_Na23_to_He4_Ne20_reaclib, rr: He4_Ne20_to_p_Na23_derived
        1, 10, -1, 2, 9, -1, 10, 9, 18, 73,  // fr: p_Al27_to_He4_Mg24_reaclib, rr: He4_Mg24_to_p_Al27_derived
        1, 12, -1, 2, 11, -1, 12, 11, 19, 74,  // fr: p_P31_to_He4_Si28_reaclib, rr: He4_Si28_to_p_P31_derived
        2, 2, 2, 3, -1, -1, 3, -1, 20, 70,  // fr: He4_He4_He4_to_C12_reaclib, rr: C12_to_He4_He4_He4_derived
        13, 2, -1, 14, -1, -1, 13, 14, 44, 45,  // fr: S32_He4_to_Ar36_approx, rr: Ar36_to_S32_He4_approx
        14, 2, -1, 15, -1, -1, 14, 15, 46, 47,  // fr: Ar36_He4_to_Ca40_approx, rr: Ca40_to_Ar36_He4_approx
        15, 2, -1, 16, -1, -1, 15, 16, 48, 49,  // fr: Ca40_He4_to_Ti44_approx, rr: Ti44_to_Ca40_He4_approx
        16, 2, -1, 17, -1, -1, 16, 17, 50, 51,  // fr: Ti44_He4_to_Cr48_approx, rr: Cr48_to_Ti44_He4_approx
        17, 2, -1, 18, -1, -1, 17, 18, 52, 53,  // fr: Cr48_He4_to_Fe52_approx, rr: Fe52_to_Cr48_He4_approx
        18, 2, -1, 19, -1, -1, 18, 19, 54, 55  // fr: Fe52_He4_to_Ni56_approx, rr: Ni56_to_Fe52_He4_approx
    };
}
#endif

#endif
