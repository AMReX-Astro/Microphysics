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
    constexpr amrex::Real bion () {
        using namespace Species;

        static_assert(spec >= 1 && spec <= NumSpec);

        // Set the binding energy of the element

        if constexpr (spec == H1) {
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
        k_N13_to_p_C12 = 1,
        k_O16_to_He4_C12 = 2,
        k_F18_to_He4_N14 = 3,
        k_Ne20_to_He4_O16 = 4,
        k_Na22_to_p_Ne21 = 5,
        k_Na22_to_He4_F18 = 6,
        k_Mg24_to_p_Na23 = 7,
        k_Mg24_to_He4_Ne20 = 8,
        k_Si28_to_p_Al27 = 9,
        k_Si28_to_He4_Mg24 = 10,
        k_S32_to_p_P31 = 11,
        k_S32_to_He4_Si28 = 12,
        k_C12_to_He4_He4_He4 = 13,
        k_p_C12_to_N13 = 14,
        k_He4_C12_to_O16 = 15,
        k_He4_N14_to_F18 = 16,
        k_He4_O16_to_Ne20 = 17,
        k_He4_F18_to_Na22 = 18,
        k_He4_Ne20_to_Mg24 = 19,
        k_p_Ne21_to_Na22 = 20,
        k_p_Na23_to_Mg24 = 21,
        k_He4_Mg24_to_Si28 = 22,
        k_p_Al27_to_Si28 = 23,
        k_He4_Si28_to_S32 = 24,
        k_p_P31_to_S32 = 25,
        k_C12_C12_to_p_Na23 = 26,
        k_C12_C12_to_He4_Ne20 = 27,
        k_He4_N13_to_p_O16 = 28,
        k_p_O16_to_He4_N13 = 29,
        k_C12_O16_to_p_Al27 = 30,
        k_C12_O16_to_He4_Mg24 = 31,
        k_O16_O16_to_p_P31 = 32,
        k_O16_O16_to_He4_Si28 = 33,
        k_He4_F18_to_p_Ne21 = 34,
        k_He4_Ne20_to_p_Na23 = 35,
        k_p_Ne21_to_He4_F18 = 36,
        k_p_Na23_to_He4_Ne20 = 37,
        k_He4_Mg24_to_p_Al27 = 38,
        k_p_Al27_to_He4_Mg24 = 39,
        k_He4_Si28_to_p_P31 = 40,
        k_p_P31_to_He4_Si28 = 41,
        k_He4_He4_He4_to_C12 = 42,
        k_C12_C12_to_Mg24_modified = 43,
        k_O16_O16_to_S32_modified = 44,
        k_C12_O16_to_Si28_modified = 45,
        k_He4_S32_to_Ar36_removed = 46,
        k_He4_S32_to_p_Cl35_removed = 47,
        k_p_Cl35_to_Ar36_removed = 48,
        k_Ar36_to_He4_S32_removed = 49,
        k_Ar36_to_p_Cl35_removed = 50,
        k_p_Cl35_to_He4_S32_removed = 51,
        k_He4_Ar36_to_Ca40_removed = 52,
        k_He4_Ar36_to_p_K39_removed = 53,
        k_p_K39_to_Ca40_removed = 54,
        k_Ca40_to_He4_Ar36_removed = 55,
        k_Ca40_to_p_K39_removed = 56,
        k_p_K39_to_He4_Ar36_removed = 57,
        k_He4_Ca40_to_Ti44_removed = 58,
        k_He4_Ca40_to_p_Sc43_removed = 59,
        k_p_Sc43_to_Ti44_removed = 60,
        k_Ti44_to_He4_Ca40_removed = 61,
        k_Ti44_to_p_Sc43_removed = 62,
        k_p_Sc43_to_He4_Ca40_removed = 63,
        k_He4_Ti44_to_Cr48_removed = 64,
        k_He4_Ti44_to_p_V47_removed = 65,
        k_p_V47_to_Cr48_removed = 66,
        k_Cr48_to_He4_Ti44_removed = 67,
        k_Cr48_to_p_V47_removed = 68,
        k_p_V47_to_He4_Ti44_removed = 69,
        k_He4_Cr48_to_Fe52_removed = 70,
        k_He4_Cr48_to_p_Mn51_removed = 71,
        k_p_Mn51_to_Fe52_removed = 72,
        k_Fe52_to_He4_Cr48_removed = 73,
        k_Fe52_to_p_Mn51_removed = 74,
        k_p_Mn51_to_He4_Cr48_removed = 75,
        k_He4_Fe52_to_Ni56_removed = 76,
        k_He4_Fe52_to_p_Co55_removed = 77,
        k_p_Co55_to_Ni56_removed = 78,
        k_Ni56_to_He4_Fe52_removed = 79,
        k_Ni56_to_p_Co55_removed = 80,
        k_p_Co55_to_He4_Fe52_removed = 81,
        k_S32_He4_to_Ar36_approx = 82,
        k_Ar36_to_S32_He4_approx = 83,
        k_Ar36_He4_to_Ca40_approx = 84,
        k_Ca40_to_Ar36_He4_approx = 85,
        k_Ca40_He4_to_Ti44_approx = 86,
        k_Ti44_to_Ca40_He4_approx = 87,
        k_Ti44_He4_to_Cr48_approx = 88,
        k_Cr48_to_Ti44_He4_approx = 89,
        k_Cr48_He4_to_Fe52_approx = 90,
        k_Fe52_to_Cr48_He4_approx = 91,
        k_Fe52_He4_to_Ni56_approx = 92,
        k_Ni56_to_Fe52_He4_approx = 93,
        NumRates = k_Ni56_to_Fe52_He4_approx
    };

    // number of reaclib rates

    const int NrateReaclib = 81;

    // number of tabular rates

    const int NrateTabular = 0;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "N13_to_p_C12",  // 1,
        "O16_to_He4_C12",  // 2,
        "F18_to_He4_N14",  // 3,
        "Ne20_to_He4_O16",  // 4,
        "Na22_to_p_Ne21",  // 5,
        "Na22_to_He4_F18",  // 6,
        "Mg24_to_p_Na23",  // 7,
        "Mg24_to_He4_Ne20",  // 8,
        "Si28_to_p_Al27",  // 9,
        "Si28_to_He4_Mg24",  // 10,
        "S32_to_p_P31",  // 11,
        "S32_to_He4_Si28",  // 12,
        "C12_to_He4_He4_He4",  // 13,
        "p_C12_to_N13",  // 14,
        "He4_C12_to_O16",  // 15,
        "He4_N14_to_F18",  // 16,
        "He4_O16_to_Ne20",  // 17,
        "He4_F18_to_Na22",  // 18,
        "He4_Ne20_to_Mg24",  // 19,
        "p_Ne21_to_Na22",  // 20,
        "p_Na23_to_Mg24",  // 21,
        "He4_Mg24_to_Si28",  // 22,
        "p_Al27_to_Si28",  // 23,
        "He4_Si28_to_S32",  // 24,
        "p_P31_to_S32",  // 25,
        "C12_C12_to_p_Na23",  // 26,
        "C12_C12_to_He4_Ne20",  // 27,
        "He4_N13_to_p_O16",  // 28,
        "p_O16_to_He4_N13",  // 29,
        "C12_O16_to_p_Al27",  // 30,
        "C12_O16_to_He4_Mg24",  // 31,
        "O16_O16_to_p_P31",  // 32,
        "O16_O16_to_He4_Si28",  // 33,
        "He4_F18_to_p_Ne21",  // 34,
        "He4_Ne20_to_p_Na23",  // 35,
        "p_Ne21_to_He4_F18",  // 36,
        "p_Na23_to_He4_Ne20",  // 37,
        "He4_Mg24_to_p_Al27",  // 38,
        "p_Al27_to_He4_Mg24",  // 39,
        "He4_Si28_to_p_P31",  // 40,
        "p_P31_to_He4_Si28",  // 41,
        "He4_He4_He4_to_C12",  // 42,
        "C12_C12_to_Mg24_modified",  // 43,
        "O16_O16_to_S32_modified",  // 44,
        "C12_O16_to_Si28_modified",  // 45,
        "He4_S32_to_Ar36_removed",  // 46,
        "He4_S32_to_p_Cl35_removed",  // 47,
        "p_Cl35_to_Ar36_removed",  // 48,
        "Ar36_to_He4_S32_removed",  // 49,
        "Ar36_to_p_Cl35_removed",  // 50,
        "p_Cl35_to_He4_S32_removed",  // 51,
        "He4_Ar36_to_Ca40_removed",  // 52,
        "He4_Ar36_to_p_K39_removed",  // 53,
        "p_K39_to_Ca40_removed",  // 54,
        "Ca40_to_He4_Ar36_removed",  // 55,
        "Ca40_to_p_K39_removed",  // 56,
        "p_K39_to_He4_Ar36_removed",  // 57,
        "He4_Ca40_to_Ti44_removed",  // 58,
        "He4_Ca40_to_p_Sc43_removed",  // 59,
        "p_Sc43_to_Ti44_removed",  // 60,
        "Ti44_to_He4_Ca40_removed",  // 61,
        "Ti44_to_p_Sc43_removed",  // 62,
        "p_Sc43_to_He4_Ca40_removed",  // 63,
        "He4_Ti44_to_Cr48_removed",  // 64,
        "He4_Ti44_to_p_V47_removed",  // 65,
        "p_V47_to_Cr48_removed",  // 66,
        "Cr48_to_He4_Ti44_removed",  // 67,
        "Cr48_to_p_V47_removed",  // 68,
        "p_V47_to_He4_Ti44_removed",  // 69,
        "He4_Cr48_to_Fe52_removed",  // 70,
        "He4_Cr48_to_p_Mn51_removed",  // 71,
        "p_Mn51_to_Fe52_removed",  // 72,
        "Fe52_to_He4_Cr48_removed",  // 73,
        "Fe52_to_p_Mn51_removed",  // 74,
        "p_Mn51_to_He4_Cr48_removed",  // 75,
        "He4_Fe52_to_Ni56_removed",  // 76,
        "He4_Fe52_to_p_Co55_removed",  // 77,
        "p_Co55_to_Ni56_removed",  // 78,
        "Ni56_to_He4_Fe52_removed",  // 79,
        "Ni56_to_p_Co55_removed",  // 80,
        "p_Co55_to_He4_Fe52_removed",  // 81,
        "S32_He4_to_Ar36_approx",  // 82,
        "Ar36_to_S32_He4_approx",  // 83,
        "Ar36_He4_to_Ca40_approx",  // 84,
        "Ca40_to_Ar36_He4_approx",  // 85,
        "Ca40_He4_to_Ti44_approx",  // 86,
        "Ti44_to_Ca40_He4_approx",  // 87,
        "Ti44_He4_to_Cr48_approx",  // 88,
        "Cr48_to_Ti44_He4_approx",  // 89,
        "Cr48_He4_to_Fe52_approx",  // 90,
        "Fe52_to_Cr48_He4_approx",  // 91,
        "Fe52_He4_to_Ni56_approx",  // 92,
        "Ni56_to_Fe52_He4_approx"  // 93,
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

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
