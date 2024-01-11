#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Loop.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex;

void actual_network_init();

const std::string network_name = "pynucastro-cxx";

namespace network
{

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr Real bion () {
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
            return 28.29566_rt;
        }
        else if constexpr (spec == C12) {
            return 92.16172800000001_rt;
        }
        else if constexpr (spec == N13) {
            return 94.105219_rt;
        }
        else if constexpr (spec == N14) {
            return 104.65859599999999_rt;
        }
        else if constexpr (spec == O16) {
            return 127.619296_rt;
        }
        else if constexpr (spec == F18) {
            return 137.369484_rt;
        }
        else if constexpr (spec == Ne20) {
            return 160.6448_rt;
        }
        else if constexpr (spec == Ne21) {
            return 167.405973_rt;
        }
        else if constexpr (spec == Na22) {
            return 174.144674_rt;
        }
        else if constexpr (spec == Na23) {
            return 186.56433900000002_rt;
        }
        else if constexpr (spec == Mg24) {
            return 198.25701600000002_rt;
        }
        else if constexpr (spec == Al27) {
            return 224.951931_rt;
        }
        else if constexpr (spec == Si28) {
            return 236.536832_rt;
        }
        else if constexpr (spec == P31) {
            return 262.91617699999995_rt;
        }
        else if constexpr (spec == S32) {
            return 271.78012800000005_rt;
        }
        else if constexpr (spec == Ar36) {
            return 306.716724_rt;
        }
        else if constexpr (spec == Ca40) {
            return 342.05212000000006_rt;
        }
        else if constexpr (spec == Ti44) {
            return 375.47488000000004_rt;
        }
        else if constexpr (spec == Cr48) {
            return 411.46891200000005_rt;
        }
        else if constexpr (spec == Fe52) {
            return 447.697848_rt;
        }
        else if constexpr (spec == Ni56) {
            return 483.995624_rt;
        }


        // Return zero if we don't recognize the species.
        return 0.0_rt;
    }

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr Real mion () {
        static_assert(spec >= 1 && spec <= NumSpec);

        constexpr Real A = NetworkProperties::aion(spec);
        constexpr Real Z = NetworkProperties::zion(spec);

        return (A - Z) * C::Legacy::m_n + Z * (C::Legacy::m_p + C::Legacy::m_e) - bion<spec>() * C::Legacy::MeV2gr;
    }

    // Legacy (non-templated) interfaces

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real bion (int spec) {
        using namespace Species;

        Real b = 0.0_rt;

        // Set the binding energy of the element
        constexpr_for<1, NumSpec+1>([&] (auto n) {
            if (n == spec) {
                b = bion<n>();
            }
        });

        return b;
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real mion (int spec) {
        using namespace Species;

        Real m = 0.0_rt;

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
        k_n_Ne20_to_Ne21 = 6,
        k_He4_Ne20_to_Mg24 = 7,
        k_p_Ne21_to_Na22 = 8,
        k_n_Na22_to_Na23 = 9,
        k_p_Na23_to_Mg24 = 10,
        k_He4_Mg24_to_Si28 = 11,
        k_p_Al27_to_Si28 = 12,
        k_He4_Si28_to_S32 = 13,
        k_p_P31_to_S32 = 14,
        k_C12_C12_to_p_Na23 = 15,
        k_C12_C12_to_He4_Ne20 = 16,
        k_He4_N13_to_p_O16 = 17,
        k_C12_O16_to_p_Al27 = 18,
        k_C12_O16_to_He4_Mg24 = 19,
        k_O16_O16_to_p_P31 = 20,
        k_O16_O16_to_He4_Si28 = 21,
        k_He4_F18_to_p_Ne21 = 22,
        k_He4_Ne21_to_n_Mg24 = 23,
        k_p_Na23_to_He4_Ne20 = 24,
        k_p_Al27_to_He4_Mg24 = 25,
        k_p_P31_to_He4_Si28 = 26,
        k_He4_He4_He4_to_C12 = 27,
        k_C12_C12_to_Mg24_modified = 28,
        k_O16_O16_to_S32_modified = 29,
        k_C12_O16_to_Si28_modified = 30,
        k_He4_S32_to_Ar36_removed = 31,
        k_p_Cl35_to_Ar36_removed = 32,
        k_p_Cl35_to_He4_S32_removed = 33,
        k_He4_Ar36_to_Ca40_removed = 34,
        k_p_K39_to_Ca40_removed = 35,
        k_p_K39_to_He4_Ar36_removed = 36,
        k_He4_Ca40_to_Ti44_removed = 37,
        k_p_Sc43_to_Ti44_removed = 38,
        k_p_Sc43_to_He4_Ca40_removed = 39,
        k_He4_Ti44_to_Cr48_removed = 40,
        k_He4_Ti44_to_p_V47_removed = 41,
        k_p_V47_to_Cr48_removed = 42,
        k_He4_Cr48_to_Fe52_removed = 43,
        k_He4_Cr48_to_p_Mn51_removed = 44,
        k_p_Mn51_to_Fe52_removed = 45,
        k_He4_Fe52_to_Ni56_removed = 46,
        k_He4_Fe52_to_p_Co55_removed = 47,
        k_p_Co55_to_Ni56_removed = 48,
        k_S32_He4_to_Ar36_approx = 49,
        k_Ar36_to_S32_He4_approx = 50,
        k_Ar36_He4_to_Ca40_approx = 51,
        k_Ca40_to_Ar36_He4_approx = 52,
        k_Ca40_He4_to_Ti44_approx = 53,
        k_Ti44_to_Ca40_He4_approx = 54,
        k_Ti44_He4_to_Cr48_approx = 55,
        k_Cr48_to_Ti44_He4_approx = 56,
        k_Cr48_He4_to_Fe52_approx = 57,
        k_Fe52_to_Cr48_He4_approx = 58,
        k_Fe52_He4_to_Ni56_approx = 59,
        k_Ni56_to_Fe52_He4_approx = 60,
        k_N13_to_p_C12_derived = 61,
        k_O16_to_He4_C12_derived = 62,
        k_F18_to_He4_N14_derived = 63,
        k_Ne20_to_He4_O16_derived = 64,
        k_Ne21_to_n_Ne20_derived = 65,
        k_Na22_to_p_Ne21_derived = 66,
        k_Na22_to_He4_F18_derived = 67,
        k_Na23_to_n_Na22_derived = 68,
        k_Mg24_to_p_Na23_derived = 69,
        k_Mg24_to_He4_Ne20_derived = 70,
        k_Si28_to_p_Al27_derived = 71,
        k_Si28_to_He4_Mg24_derived = 72,
        k_S32_to_p_P31_derived = 73,
        k_S32_to_He4_Si28_derived = 74,
        k_C12_to_He4_He4_He4_derived = 75,
        k_p_O16_to_He4_N13_derived = 76,
        k_He4_Ne20_to_p_Na23_derived = 77,
        k_p_Ne21_to_He4_F18_derived = 78,
        k_n_Mg24_to_He4_Ne21_derived = 79,
        k_He4_Mg24_to_p_Al27_derived = 80,
        k_He4_Si28_to_p_P31_derived = 81,
        k_He4_S32_to_p_Cl35_derived_removed = 82,
        k_Ar36_to_He4_S32_derived_removed = 83,
        k_Ar36_to_p_Cl35_derived_removed = 84,
        k_He4_Ar36_to_p_K39_derived_removed = 85,
        k_Ca40_to_He4_Ar36_derived_removed = 86,
        k_Ca40_to_p_K39_derived_removed = 87,
        k_He4_Ca40_to_p_Sc43_derived_removed = 88,
        k_Ti44_to_He4_Ca40_derived_removed = 89,
        k_Ti44_to_p_Sc43_derived_removed = 90,
        k_Cr48_to_He4_Ti44_derived_removed = 91,
        k_Cr48_to_p_V47_derived_removed = 92,
        k_p_V47_to_He4_Ti44_derived_removed = 93,
        k_Fe52_to_He4_Cr48_derived_removed = 94,
        k_Fe52_to_p_Mn51_derived_removed = 95,
        k_p_Mn51_to_He4_Cr48_derived_removed = 96,
        k_Ni56_to_He4_Fe52_derived_removed = 97,
        k_Ni56_to_p_Co55_derived_removed = 98,
        k_p_Co55_to_He4_Fe52_derived_removed = 99,
        NumRates = k_p_Co55_to_He4_Fe52_derived_removed
    };

    // number of reaclib rates

    const int NrateReaclib = 87;

    // number of tabular rates

    const int NrateTabular = 0;

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
        "n_Ne20_to_Ne21",  // 6,
        "He4_Ne20_to_Mg24",  // 7,
        "p_Ne21_to_Na22",  // 8,
        "n_Na22_to_Na23",  // 9,
        "p_Na23_to_Mg24",  // 10,
        "He4_Mg24_to_Si28",  // 11,
        "p_Al27_to_Si28",  // 12,
        "He4_Si28_to_S32",  // 13,
        "p_P31_to_S32",  // 14,
        "C12_C12_to_p_Na23",  // 15,
        "C12_C12_to_He4_Ne20",  // 16,
        "He4_N13_to_p_O16",  // 17,
        "C12_O16_to_p_Al27",  // 18,
        "C12_O16_to_He4_Mg24",  // 19,
        "O16_O16_to_p_P31",  // 20,
        "O16_O16_to_He4_Si28",  // 21,
        "He4_F18_to_p_Ne21",  // 22,
        "He4_Ne21_to_n_Mg24",  // 23,
        "p_Na23_to_He4_Ne20",  // 24,
        "p_Al27_to_He4_Mg24",  // 25,
        "p_P31_to_He4_Si28",  // 26,
        "He4_He4_He4_to_C12",  // 27,
        "C12_C12_to_Mg24_modified",  // 28,
        "O16_O16_to_S32_modified",  // 29,
        "C12_O16_to_Si28_modified",  // 30,
        "He4_S32_to_Ar36_removed",  // 31,
        "p_Cl35_to_Ar36_removed",  // 32,
        "p_Cl35_to_He4_S32_removed",  // 33,
        "He4_Ar36_to_Ca40_removed",  // 34,
        "p_K39_to_Ca40_removed",  // 35,
        "p_K39_to_He4_Ar36_removed",  // 36,
        "He4_Ca40_to_Ti44_removed",  // 37,
        "p_Sc43_to_Ti44_removed",  // 38,
        "p_Sc43_to_He4_Ca40_removed",  // 39,
        "He4_Ti44_to_Cr48_removed",  // 40,
        "He4_Ti44_to_p_V47_removed",  // 41,
        "p_V47_to_Cr48_removed",  // 42,
        "He4_Cr48_to_Fe52_removed",  // 43,
        "He4_Cr48_to_p_Mn51_removed",  // 44,
        "p_Mn51_to_Fe52_removed",  // 45,
        "He4_Fe52_to_Ni56_removed",  // 46,
        "He4_Fe52_to_p_Co55_removed",  // 47,
        "p_Co55_to_Ni56_removed",  // 48,
        "S32_He4_to_Ar36_approx",  // 49,
        "Ar36_to_S32_He4_approx",  // 50,
        "Ar36_He4_to_Ca40_approx",  // 51,
        "Ca40_to_Ar36_He4_approx",  // 52,
        "Ca40_He4_to_Ti44_approx",  // 53,
        "Ti44_to_Ca40_He4_approx",  // 54,
        "Ti44_He4_to_Cr48_approx",  // 55,
        "Cr48_to_Ti44_He4_approx",  // 56,
        "Cr48_He4_to_Fe52_approx",  // 57,
        "Fe52_to_Cr48_He4_approx",  // 58,
        "Fe52_He4_to_Ni56_approx",  // 59,
        "Ni56_to_Fe52_He4_approx",  // 60,
        "N13_to_p_C12_derived",  // 61,
        "O16_to_He4_C12_derived",  // 62,
        "F18_to_He4_N14_derived",  // 63,
        "Ne20_to_He4_O16_derived",  // 64,
        "Ne21_to_n_Ne20_derived",  // 65,
        "Na22_to_p_Ne21_derived",  // 66,
        "Na22_to_He4_F18_derived",  // 67,
        "Na23_to_n_Na22_derived",  // 68,
        "Mg24_to_p_Na23_derived",  // 69,
        "Mg24_to_He4_Ne20_derived",  // 70,
        "Si28_to_p_Al27_derived",  // 71,
        "Si28_to_He4_Mg24_derived",  // 72,
        "S32_to_p_P31_derived",  // 73,
        "S32_to_He4_Si28_derived",  // 74,
        "C12_to_He4_He4_He4_derived",  // 75,
        "p_O16_to_He4_N13_derived",  // 76,
        "He4_Ne20_to_p_Na23_derived",  // 77,
        "p_Ne21_to_He4_F18_derived",  // 78,
        "n_Mg24_to_He4_Ne21_derived",  // 79,
        "He4_Mg24_to_p_Al27_derived",  // 80,
        "He4_Si28_to_p_P31_derived",  // 81,
        "He4_S32_to_p_Cl35_derived_removed",  // 82,
        "Ar36_to_He4_S32_derived_removed",  // 83,
        "Ar36_to_p_Cl35_derived_removed",  // 84,
        "He4_Ar36_to_p_K39_derived_removed",  // 85,
        "Ca40_to_He4_Ar36_derived_removed",  // 86,
        "Ca40_to_p_K39_derived_removed",  // 87,
        "He4_Ca40_to_p_Sc43_derived_removed",  // 88,
        "Ti44_to_He4_Ca40_derived_removed",  // 89,
        "Ti44_to_p_Sc43_derived_removed",  // 90,
        "Cr48_to_He4_Ti44_derived_removed",  // 91,
        "Cr48_to_p_V47_derived_removed",  // 92,
        "p_V47_to_He4_Ti44_derived_removed",  // 93,
        "Fe52_to_He4_Cr48_derived_removed",  // 94,
        "Fe52_to_p_Mn51_derived_removed",  // 95,
        "p_Mn51_to_He4_Cr48_derived_removed",  // 96,
        "Ni56_to_He4_Fe52_derived_removed",  // 97,
        "Ni56_to_p_Co55_derived_removed",  // 98,
        "p_Co55_to_He4_Fe52_derived_removed"  // 99,
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

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
