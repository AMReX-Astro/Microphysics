#ifndef actual_network_H
#define actual_network_H

#include <string_view>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Loop.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex;

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
        else if constexpr (spec == P_nse) {
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
        else if constexpr (spec == F18) {
            return 137.36950247599816_rt;
        }
        else if constexpr (spec == Ne20) {
            return 160.64482384000075_rt;
        }
        else if constexpr (spec == Ne21) {
            return 167.40598973999658_rt;
        }
        else if constexpr (spec == Na22) {
            return 174.14457080400098_rt;
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
        static_assert(spec >= 1 && spec <= NumSpec);

        constexpr amrex::Real A = NetworkProperties::aion(spec);
        constexpr amrex::Real Z = NetworkProperties::zion(spec);

        return (A - Z) * C::Legacy::m_n + Z * (C::Legacy::m_p + C::Legacy::m_e) - bion<spec>() * C::Legacy::MeV2gr;
    }

    // Legacy (non-templated) interfaces

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real bion (int spec) {
        using namespace Species;

        amrex::Real b = 0.0_rt;

        // Set the binding energy of the element
        constexpr_for<1, NumSpec+1>([&] (auto n) {
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

        constexpr_for<1, NumSpec+1>([&] (auto n) {
            if (n == spec) {
                m = mion<n>();
            }
        });

        return m;
    }
}

namespace Rates
{

    enum NetworkRates
    {
        k_p_C12_to_N13 = 1,
        k_He4_C12_to_O16 = 2,
        k_He4_N14_to_F18 = 3,
        k_He4_O16_to_Ne20 = 4,
        k_He4_F18_to_Na22 = 5,
        k_He4_Ne20_to_Mg24 = 6,
        k_p_Ne21_to_Na22 = 7,
        k_p_Na23_to_Mg24 = 8,
        k_He4_Mg24_to_Si28 = 9,
        k_p_Al27_to_Si28 = 10,
        k_He4_Si28_to_S32 = 11,
        k_p_P31_to_S32 = 12,
        k_He4_Cr48_to_Fe52 = 13,
        k_p_nse_Mn51_to_Fe52 = 14,
        k_He4_Mn51_to_Co55 = 15,
        k_He4_Fe52_to_Ni56 = 16,
        k_p_nse_Co55_to_Ni56 = 17,
        k_C12_C12_to_p_Na23 = 18,
        k_C12_C12_to_He4_Ne20 = 19,
        k_He4_N13_to_p_O16 = 20,
        k_C12_O16_to_p_Al27 = 21,
        k_C12_O16_to_He4_Mg24 = 22,
        k_O16_O16_to_p_P31 = 23,
        k_O16_O16_to_He4_Si28 = 24,
        k_He4_F18_to_p_Ne21 = 25,
        k_p_Na23_to_He4_Ne20 = 26,
        k_p_Al27_to_He4_Mg24 = 27,
        k_p_P31_to_He4_Si28 = 28,
        k_He4_Cr48_to_p_nse_Mn51 = 29,
        k_He4_Fe52_to_p_nse_Co55 = 30,
        k_He4_He4_He4_to_C12 = 31,
        k_C12_C12_to_Mg24_modified = 32,
        k_O16_O16_to_S32_modified = 33,
        k_C12_O16_to_Si28_modified = 34,
        k_p_nse_Fe54_to_Co55 = 35,
        k_He4_Fe54_to_Ni58 = 36,
        k_p_nse_Fe56_to_Co57 = 37,
        k_n_Co55_to_Co56 = 38,
        k_n_Co56_to_Co57 = 39,
        k_p_nse_Co57_to_Ni58 = 40,
        k_He4_Mn51_to_p_nse_Fe54 = 41,
        k_He4_Co55_to_p_nse_Ni58 = 42,
        k_n_Co56_to_p_nse_Fe56 = 43,
        k_p_nse_Co57_to_He4_Fe54 = 44,
        k_n_Ni56_to_p_nse_Co56 = 45,
        k_He4_S32_to_Ar36_removed = 46,
        k_p_Cl35_to_Ar36_removed = 47,
        k_p_Cl35_to_He4_S32_removed = 48,
        k_He4_Ar36_to_Ca40_removed = 49,
        k_p_K39_to_Ca40_removed = 50,
        k_p_K39_to_He4_Ar36_removed = 51,
        k_He4_Ca40_to_Ti44_removed = 52,
        k_p_Sc43_to_Ti44_removed = 53,
        k_p_Sc43_to_He4_Ca40_removed = 54,
        k_He4_Ti44_to_Cr48_removed = 55,
        k_He4_Ti44_to_p_V47_removed = 56,
        k_p_V47_to_Cr48_removed = 57,
        k_n_Fe52_to_Fe53_removed = 58,
        k_n_Fe53_to_Fe54_removed = 59,
        k_n_Fe54_to_Fe55_removed = 60,
        k_n_Fe55_to_Fe56_removed = 61,
        k_n_Ni56_to_Ni57_removed = 62,
        k_n_Ni57_to_Ni58_removed = 63,
        k_Co56_to_Fe56 = 64,
        k_Co56_to_Ni56 = 65,
        k_Fe56_to_Co56 = 66,
        k_n_to_p = 67,
        k_Ni56_to_Co56 = 68,
        k_p_to_n = 69,
        k_S32_He4_to_Ar36_approx = 70,
        k_Ar36_to_S32_He4_approx = 71,
        k_Ar36_He4_to_Ca40_approx = 72,
        k_Ca40_to_Ar36_He4_approx = 73,
        k_Ca40_He4_to_Ti44_approx = 74,
        k_Ti44_to_Ca40_He4_approx = 75,
        k_Ti44_He4_to_Cr48_approx = 76,
        k_Cr48_to_Ti44_He4_approx = 77,
        k_Fe52_n_n_to_Fe54_approx = 78,
        k_Fe54_to_Fe52_n_n_approx = 79,
        k_Fe54_n_n_to_Fe56_approx = 80,
        k_Fe56_to_Fe54_n_n_approx = 81,
        k_Ni56_n_n_to_Ni58_approx = 82,
        k_Ni58_to_Ni56_n_n_approx = 83,
        k_N13_to_p_C12_derived = 84,
        k_O16_to_He4_C12_derived = 85,
        k_F18_to_He4_N14_derived = 86,
        k_Ne20_to_He4_O16_derived = 87,
        k_Na22_to_p_Ne21_derived = 88,
        k_Na22_to_He4_F18_derived = 89,
        k_Mg24_to_p_Na23_derived = 90,
        k_Mg24_to_He4_Ne20_derived = 91,
        k_Si28_to_p_Al27_derived = 92,
        k_Si28_to_He4_Mg24_derived = 93,
        k_S32_to_p_P31_derived = 94,
        k_S32_to_He4_Si28_derived = 95,
        k_Fe52_to_p_nse_Mn51_derived = 96,
        k_Fe52_to_He4_Cr48_derived = 97,
        k_Co55_to_He4_Mn51_derived = 98,
        k_Ni56_to_p_nse_Co55_derived = 99,
        k_Ni56_to_He4_Fe52_derived = 100,
        k_C12_to_He4_He4_He4_derived = 101,
        k_p_O16_to_He4_N13_derived = 102,
        k_He4_Ne20_to_p_Na23_derived = 103,
        k_p_Ne21_to_He4_F18_derived = 104,
        k_He4_Mg24_to_p_Al27_derived = 105,
        k_He4_Si28_to_p_P31_derived = 106,
        k_p_nse_Mn51_to_He4_Cr48_derived = 107,
        k_p_nse_Co55_to_He4_Fe52_derived = 108,
        k_Co55_to_p_nse_Fe54_derived = 109,
        k_Co56_to_n_Co55_derived = 110,
        k_Co57_to_n_Co56_derived = 111,
        k_Co57_to_p_nse_Fe56_derived = 112,
        k_Ni58_to_p_nse_Co57_derived = 113,
        k_Ni58_to_He4_Fe54_derived = 114,
        k_p_nse_Fe54_to_He4_Mn51_derived = 115,
        k_He4_Fe54_to_p_nse_Co57_derived = 116,
        k_p_nse_Fe56_to_n_Co56_derived = 117,
        k_p_nse_Co56_to_n_Ni56_derived = 118,
        k_p_nse_Ni58_to_He4_Co55_derived = 119,
        k_He4_S32_to_p_Cl35_derived_removed = 120,
        k_Ar36_to_He4_S32_derived_removed = 121,
        k_Ar36_to_p_Cl35_derived_removed = 122,
        k_He4_Ar36_to_p_K39_derived_removed = 123,
        k_Ca40_to_He4_Ar36_derived_removed = 124,
        k_Ca40_to_p_K39_derived_removed = 125,
        k_He4_Ca40_to_p_Sc43_derived_removed = 126,
        k_Ti44_to_He4_Ca40_derived_removed = 127,
        k_Ti44_to_p_Sc43_derived_removed = 128,
        k_Cr48_to_He4_Ti44_derived_removed = 129,
        k_Cr48_to_p_V47_derived_removed = 130,
        k_p_V47_to_He4_Ti44_derived_removed = 131,
        k_Fe54_to_n_Fe53_derived_removed = 132,
        k_Fe53_to_n_Fe52_derived_removed = 133,
        k_Fe56_to_n_Fe55_derived_removed = 134,
        k_Fe55_to_n_Fe54_derived_removed = 135,
        k_Ni58_to_n_Ni57_derived_removed = 136,
        k_Ni57_to_n_Ni56_derived_removed = 137,
        NumRates = k_Ni57_to_n_Ni56_derived_removed
    };

    // number of reaclib rates

    const int NrateReaclib = 117;

    // number of tabular rates

    const int NrateTabular = 6;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "p_C12_to_N13",  // 1,
        "He4_C12_to_O16",  // 2,
        "He4_N14_to_F18",  // 3,
        "He4_O16_to_Ne20",  // 4,
        "He4_F18_to_Na22",  // 5,
        "He4_Ne20_to_Mg24",  // 6,
        "p_Ne21_to_Na22",  // 7,
        "p_Na23_to_Mg24",  // 8,
        "He4_Mg24_to_Si28",  // 9,
        "p_Al27_to_Si28",  // 10,
        "He4_Si28_to_S32",  // 11,
        "p_P31_to_S32",  // 12,
        "He4_Cr48_to_Fe52",  // 13,
        "p_nse_Mn51_to_Fe52",  // 14,
        "He4_Mn51_to_Co55",  // 15,
        "He4_Fe52_to_Ni56",  // 16,
        "p_nse_Co55_to_Ni56",  // 17,
        "C12_C12_to_p_Na23",  // 18,
        "C12_C12_to_He4_Ne20",  // 19,
        "He4_N13_to_p_O16",  // 20,
        "C12_O16_to_p_Al27",  // 21,
        "C12_O16_to_He4_Mg24",  // 22,
        "O16_O16_to_p_P31",  // 23,
        "O16_O16_to_He4_Si28",  // 24,
        "He4_F18_to_p_Ne21",  // 25,
        "p_Na23_to_He4_Ne20",  // 26,
        "p_Al27_to_He4_Mg24",  // 27,
        "p_P31_to_He4_Si28",  // 28,
        "He4_Cr48_to_p_nse_Mn51",  // 29,
        "He4_Fe52_to_p_nse_Co55",  // 30,
        "He4_He4_He4_to_C12",  // 31,
        "C12_C12_to_Mg24_modified",  // 32,
        "O16_O16_to_S32_modified",  // 33,
        "C12_O16_to_Si28_modified",  // 34,
        "p_nse_Fe54_to_Co55",  // 35,
        "He4_Fe54_to_Ni58",  // 36,
        "p_nse_Fe56_to_Co57",  // 37,
        "n_Co55_to_Co56",  // 38,
        "n_Co56_to_Co57",  // 39,
        "p_nse_Co57_to_Ni58",  // 40,
        "He4_Mn51_to_p_nse_Fe54",  // 41,
        "He4_Co55_to_p_nse_Ni58",  // 42,
        "n_Co56_to_p_nse_Fe56",  // 43,
        "p_nse_Co57_to_He4_Fe54",  // 44,
        "n_Ni56_to_p_nse_Co56",  // 45,
        "He4_S32_to_Ar36_removed",  // 46,
        "p_Cl35_to_Ar36_removed",  // 47,
        "p_Cl35_to_He4_S32_removed",  // 48,
        "He4_Ar36_to_Ca40_removed",  // 49,
        "p_K39_to_Ca40_removed",  // 50,
        "p_K39_to_He4_Ar36_removed",  // 51,
        "He4_Ca40_to_Ti44_removed",  // 52,
        "p_Sc43_to_Ti44_removed",  // 53,
        "p_Sc43_to_He4_Ca40_removed",  // 54,
        "He4_Ti44_to_Cr48_removed",  // 55,
        "He4_Ti44_to_p_V47_removed",  // 56,
        "p_V47_to_Cr48_removed",  // 57,
        "n_Fe52_to_Fe53_removed",  // 58,
        "n_Fe53_to_Fe54_removed",  // 59,
        "n_Fe54_to_Fe55_removed",  // 60,
        "n_Fe55_to_Fe56_removed",  // 61,
        "n_Ni56_to_Ni57_removed",  // 62,
        "n_Ni57_to_Ni58_removed",  // 63,
        "Co56_to_Fe56",  // 64,
        "Co56_to_Ni56",  // 65,
        "Fe56_to_Co56",  // 66,
        "n_to_p",  // 67,
        "Ni56_to_Co56",  // 68,
        "p_to_n",  // 69,
        "S32_He4_to_Ar36_approx",  // 70,
        "Ar36_to_S32_He4_approx",  // 71,
        "Ar36_He4_to_Ca40_approx",  // 72,
        "Ca40_to_Ar36_He4_approx",  // 73,
        "Ca40_He4_to_Ti44_approx",  // 74,
        "Ti44_to_Ca40_He4_approx",  // 75,
        "Ti44_He4_to_Cr48_approx",  // 76,
        "Cr48_to_Ti44_He4_approx",  // 77,
        "Fe52_n_n_to_Fe54_approx",  // 78,
        "Fe54_to_Fe52_n_n_approx",  // 79,
        "Fe54_n_n_to_Fe56_approx",  // 80,
        "Fe56_to_Fe54_n_n_approx",  // 81,
        "Ni56_n_n_to_Ni58_approx",  // 82,
        "Ni58_to_Ni56_n_n_approx",  // 83,
        "N13_to_p_C12_derived",  // 84,
        "O16_to_He4_C12_derived",  // 85,
        "F18_to_He4_N14_derived",  // 86,
        "Ne20_to_He4_O16_derived",  // 87,
        "Na22_to_p_Ne21_derived",  // 88,
        "Na22_to_He4_F18_derived",  // 89,
        "Mg24_to_p_Na23_derived",  // 90,
        "Mg24_to_He4_Ne20_derived",  // 91,
        "Si28_to_p_Al27_derived",  // 92,
        "Si28_to_He4_Mg24_derived",  // 93,
        "S32_to_p_P31_derived",  // 94,
        "S32_to_He4_Si28_derived",  // 95,
        "Fe52_to_p_nse_Mn51_derived",  // 96,
        "Fe52_to_He4_Cr48_derived",  // 97,
        "Co55_to_He4_Mn51_derived",  // 98,
        "Ni56_to_p_nse_Co55_derived",  // 99,
        "Ni56_to_He4_Fe52_derived",  // 100,
        "C12_to_He4_He4_He4_derived",  // 101,
        "p_O16_to_He4_N13_derived",  // 102,
        "He4_Ne20_to_p_Na23_derived",  // 103,
        "p_Ne21_to_He4_F18_derived",  // 104,
        "He4_Mg24_to_p_Al27_derived",  // 105,
        "He4_Si28_to_p_P31_derived",  // 106,
        "p_nse_Mn51_to_He4_Cr48_derived",  // 107,
        "p_nse_Co55_to_He4_Fe52_derived",  // 108,
        "Co55_to_p_nse_Fe54_derived",  // 109,
        "Co56_to_n_Co55_derived",  // 110,
        "Co57_to_n_Co56_derived",  // 111,
        "Co57_to_p_nse_Fe56_derived",  // 112,
        "Ni58_to_p_nse_Co57_derived",  // 113,
        "Ni58_to_He4_Fe54_derived",  // 114,
        "p_nse_Fe54_to_He4_Mn51_derived",  // 115,
        "He4_Fe54_to_p_nse_Co57_derived",  // 116,
        "p_nse_Fe56_to_n_Co56_derived",  // 117,
        "p_nse_Co56_to_n_Ni56_derived",  // 118,
        "p_nse_Ni58_to_He4_Co55_derived",  // 119,
        "He4_S32_to_p_Cl35_derived_removed",  // 120,
        "Ar36_to_He4_S32_derived_removed",  // 121,
        "Ar36_to_p_Cl35_derived_removed",  // 122,
        "He4_Ar36_to_p_K39_derived_removed",  // 123,
        "Ca40_to_He4_Ar36_derived_removed",  // 124,
        "Ca40_to_p_K39_derived_removed",  // 125,
        "He4_Ca40_to_p_Sc43_derived_removed",  // 126,
        "Ti44_to_He4_Ca40_derived_removed",  // 127,
        "Ti44_to_p_Sc43_derived_removed",  // 128,
        "Cr48_to_He4_Ti44_derived_removed",  // 129,
        "Cr48_to_p_V47_derived_removed",  // 130,
        "p_V47_to_He4_Ti44_derived_removed",  // 131,
        "Fe54_to_n_Fe53_derived_removed",  // 132,
        "Fe53_to_n_Fe52_derived_removed",  // 133,
        "Fe56_to_n_Fe55_derived_removed",  // 134,
        "Fe55_to_n_Fe54_derived_removed",  // 135,
        "Ni58_to_n_Ni57_derived_removed",  // 136,
        "Ni57_to_n_Ni56_derived_removed"  // 137,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int H1_index = 1;
    constexpr int N_index = 0;
    constexpr int He4_index = 3;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
