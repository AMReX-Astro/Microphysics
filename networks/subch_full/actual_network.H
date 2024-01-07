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
        else if constexpr (spec == Cl35) {
            return 298.20973000000004_rt;
        }
        else if constexpr (spec == Ar36) {
            return 306.716724_rt;
        }
        else if constexpr (spec == K39) {
            return 333.723975_rt;
        }
        else if constexpr (spec == Ca40) {
            return 342.05212000000006_rt;
        }
        else if constexpr (spec == Sc43) {
            return 366.825475_rt;
        }
        else if constexpr (spec == Ti44) {
            return 375.47488000000004_rt;
        }
        else if constexpr (spec == V47) {
            return 403.36457500000006_rt;
        }
        else if constexpr (spec == Cr48) {
            return 411.46891200000005_rt;
        }
        else if constexpr (spec == Mn51) {
            return 440.32237200000003_rt;
        }
        else if constexpr (spec == Fe52) {
            return 447.697848_rt;
        }
        else if constexpr (spec == Co55) {
            return 476.82899_rt;
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
        k_N13_to_p_C12 = 1,
        k_O16_to_He4_C12 = 2,
        k_F18_to_He4_N14 = 3,
        k_Ne20_to_He4_O16 = 4,
        k_Na22_to_p_Ne21 = 5,
        k_Na22_to_He4_F18 = 6,
        k_Mg24_to_p_Na23 = 7,
        k_Mg24_to_He4_Ne20 = 8,
        k_Al27_to_He4_Na23 = 9,
        k_Si28_to_p_Al27 = 10,
        k_Si28_to_He4_Mg24 = 11,
        k_P31_to_He4_Al27 = 12,
        k_S32_to_p_P31 = 13,
        k_S32_to_He4_Si28 = 14,
        k_Cl35_to_He4_P31 = 15,
        k_Ar36_to_p_Cl35 = 16,
        k_Ar36_to_He4_S32 = 17,
        k_K39_to_He4_Cl35 = 18,
        k_Ca40_to_p_K39 = 19,
        k_Ca40_to_He4_Ar36 = 20,
        k_Sc43_to_He4_K39 = 21,
        k_Ti44_to_p_Sc43 = 22,
        k_Ti44_to_He4_Ca40 = 23,
        k_V47_to_He4_Sc43 = 24,
        k_Cr48_to_p_V47 = 25,
        k_Cr48_to_He4_Ti44 = 26,
        k_Mn51_to_He4_V47 = 27,
        k_Fe52_to_p_Mn51 = 28,
        k_Fe52_to_He4_Cr48 = 29,
        k_Co55_to_He4_Mn51 = 30,
        k_Ni56_to_p_Co55 = 31,
        k_Ni56_to_He4_Fe52 = 32,
        k_C12_to_He4_He4_He4 = 33,
        k_p_C12_to_N13 = 34,
        k_He4_C12_to_O16 = 35,
        k_He4_N14_to_F18 = 36,
        k_He4_O16_to_Ne20 = 37,
        k_He4_F18_to_Na22 = 38,
        k_He4_Ne20_to_Mg24 = 39,
        k_p_Ne21_to_Na22 = 40,
        k_p_Na23_to_Mg24 = 41,
        k_He4_Na23_to_Al27 = 42,
        k_He4_Mg24_to_Si28 = 43,
        k_p_Al27_to_Si28 = 44,
        k_He4_Al27_to_P31 = 45,
        k_He4_Si28_to_S32 = 46,
        k_p_P31_to_S32 = 47,
        k_He4_P31_to_Cl35 = 48,
        k_He4_S32_to_Ar36 = 49,
        k_p_Cl35_to_Ar36 = 50,
        k_He4_Cl35_to_K39 = 51,
        k_He4_Ar36_to_Ca40 = 52,
        k_p_K39_to_Ca40 = 53,
        k_He4_K39_to_Sc43 = 54,
        k_He4_Ca40_to_Ti44 = 55,
        k_p_Sc43_to_Ti44 = 56,
        k_He4_Sc43_to_V47 = 57,
        k_He4_Ti44_to_Cr48 = 58,
        k_p_V47_to_Cr48 = 59,
        k_He4_V47_to_Mn51 = 60,
        k_He4_Cr48_to_Fe52 = 61,
        k_p_Mn51_to_Fe52 = 62,
        k_He4_Mn51_to_Co55 = 63,
        k_He4_Fe52_to_Ni56 = 64,
        k_p_Co55_to_Ni56 = 65,
        k_C12_C12_to_p_Na23 = 66,
        k_C12_C12_to_He4_Ne20 = 67,
        k_He4_N13_to_p_O16 = 68,
        k_p_O16_to_He4_N13 = 69,
        k_C12_O16_to_p_Al27 = 70,
        k_C12_O16_to_He4_Mg24 = 71,
        k_O16_O16_to_p_P31 = 72,
        k_O16_O16_to_He4_Si28 = 73,
        k_He4_F18_to_p_Ne21 = 74,
        k_He4_Ne20_to_p_Na23 = 75,
        k_He4_Ne20_to_C12_C12 = 76,
        k_C12_Ne20_to_p_P31 = 77,
        k_C12_Ne20_to_He4_Si28 = 78,
        k_p_Ne21_to_He4_F18 = 79,
        k_p_Na23_to_He4_Ne20 = 80,
        k_p_Na23_to_C12_C12 = 81,
        k_He4_Mg24_to_p_Al27 = 82,
        k_He4_Mg24_to_C12_O16 = 83,
        k_p_Al27_to_He4_Mg24 = 84,
        k_p_Al27_to_C12_O16 = 85,
        k_He4_Si28_to_p_P31 = 86,
        k_He4_Si28_to_C12_Ne20 = 87,
        k_He4_Si28_to_O16_O16 = 88,
        k_p_P31_to_He4_Si28 = 89,
        k_p_P31_to_C12_Ne20 = 90,
        k_p_P31_to_O16_O16 = 91,
        k_He4_S32_to_p_Cl35 = 92,
        k_p_Cl35_to_He4_S32 = 93,
        k_He4_Ar36_to_p_K39 = 94,
        k_p_K39_to_He4_Ar36 = 95,
        k_He4_Ca40_to_p_Sc43 = 96,
        k_p_Sc43_to_He4_Ca40 = 97,
        k_He4_Ti44_to_p_V47 = 98,
        k_p_V47_to_He4_Ti44 = 99,
        k_He4_Cr48_to_p_Mn51 = 100,
        k_p_Mn51_to_He4_Cr48 = 101,
        k_He4_Fe52_to_p_Co55 = 102,
        k_p_Co55_to_He4_Fe52 = 103,
        k_He4_He4_He4_to_C12 = 104,
        k_C12_C12_to_Mg24_modified = 105,
        k_O16_O16_to_S32_modified = 106,
        k_C12_O16_to_Si28_modified = 107,
        NumRates = k_C12_O16_to_Si28_modified
    };

    // number of reaclib rates

    const int NrateReaclib = 107;

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
        "Al27_to_He4_Na23",  // 9,
        "Si28_to_p_Al27",  // 10,
        "Si28_to_He4_Mg24",  // 11,
        "P31_to_He4_Al27",  // 12,
        "S32_to_p_P31",  // 13,
        "S32_to_He4_Si28",  // 14,
        "Cl35_to_He4_P31",  // 15,
        "Ar36_to_p_Cl35",  // 16,
        "Ar36_to_He4_S32",  // 17,
        "K39_to_He4_Cl35",  // 18,
        "Ca40_to_p_K39",  // 19,
        "Ca40_to_He4_Ar36",  // 20,
        "Sc43_to_He4_K39",  // 21,
        "Ti44_to_p_Sc43",  // 22,
        "Ti44_to_He4_Ca40",  // 23,
        "V47_to_He4_Sc43",  // 24,
        "Cr48_to_p_V47",  // 25,
        "Cr48_to_He4_Ti44",  // 26,
        "Mn51_to_He4_V47",  // 27,
        "Fe52_to_p_Mn51",  // 28,
        "Fe52_to_He4_Cr48",  // 29,
        "Co55_to_He4_Mn51",  // 30,
        "Ni56_to_p_Co55",  // 31,
        "Ni56_to_He4_Fe52",  // 32,
        "C12_to_He4_He4_He4",  // 33,
        "p_C12_to_N13",  // 34,
        "He4_C12_to_O16",  // 35,
        "He4_N14_to_F18",  // 36,
        "He4_O16_to_Ne20",  // 37,
        "He4_F18_to_Na22",  // 38,
        "He4_Ne20_to_Mg24",  // 39,
        "p_Ne21_to_Na22",  // 40,
        "p_Na23_to_Mg24",  // 41,
        "He4_Na23_to_Al27",  // 42,
        "He4_Mg24_to_Si28",  // 43,
        "p_Al27_to_Si28",  // 44,
        "He4_Al27_to_P31",  // 45,
        "He4_Si28_to_S32",  // 46,
        "p_P31_to_S32",  // 47,
        "He4_P31_to_Cl35",  // 48,
        "He4_S32_to_Ar36",  // 49,
        "p_Cl35_to_Ar36",  // 50,
        "He4_Cl35_to_K39",  // 51,
        "He4_Ar36_to_Ca40",  // 52,
        "p_K39_to_Ca40",  // 53,
        "He4_K39_to_Sc43",  // 54,
        "He4_Ca40_to_Ti44",  // 55,
        "p_Sc43_to_Ti44",  // 56,
        "He4_Sc43_to_V47",  // 57,
        "He4_Ti44_to_Cr48",  // 58,
        "p_V47_to_Cr48",  // 59,
        "He4_V47_to_Mn51",  // 60,
        "He4_Cr48_to_Fe52",  // 61,
        "p_Mn51_to_Fe52",  // 62,
        "He4_Mn51_to_Co55",  // 63,
        "He4_Fe52_to_Ni56",  // 64,
        "p_Co55_to_Ni56",  // 65,
        "C12_C12_to_p_Na23",  // 66,
        "C12_C12_to_He4_Ne20",  // 67,
        "He4_N13_to_p_O16",  // 68,
        "p_O16_to_He4_N13",  // 69,
        "C12_O16_to_p_Al27",  // 70,
        "C12_O16_to_He4_Mg24",  // 71,
        "O16_O16_to_p_P31",  // 72,
        "O16_O16_to_He4_Si28",  // 73,
        "He4_F18_to_p_Ne21",  // 74,
        "He4_Ne20_to_p_Na23",  // 75,
        "He4_Ne20_to_C12_C12",  // 76,
        "C12_Ne20_to_p_P31",  // 77,
        "C12_Ne20_to_He4_Si28",  // 78,
        "p_Ne21_to_He4_F18",  // 79,
        "p_Na23_to_He4_Ne20",  // 80,
        "p_Na23_to_C12_C12",  // 81,
        "He4_Mg24_to_p_Al27",  // 82,
        "He4_Mg24_to_C12_O16",  // 83,
        "p_Al27_to_He4_Mg24",  // 84,
        "p_Al27_to_C12_O16",  // 85,
        "He4_Si28_to_p_P31",  // 86,
        "He4_Si28_to_C12_Ne20",  // 87,
        "He4_Si28_to_O16_O16",  // 88,
        "p_P31_to_He4_Si28",  // 89,
        "p_P31_to_C12_Ne20",  // 90,
        "p_P31_to_O16_O16",  // 91,
        "He4_S32_to_p_Cl35",  // 92,
        "p_Cl35_to_He4_S32",  // 93,
        "He4_Ar36_to_p_K39",  // 94,
        "p_K39_to_He4_Ar36",  // 95,
        "He4_Ca40_to_p_Sc43",  // 96,
        "p_Sc43_to_He4_Ca40",  // 97,
        "He4_Ti44_to_p_V47",  // 98,
        "p_V47_to_He4_Ti44",  // 99,
        "He4_Cr48_to_p_Mn51",  // 100,
        "p_Mn51_to_He4_Cr48",  // 101,
        "He4_Fe52_to_p_Co55",  // 102,
        "p_Co55_to_He4_Fe52",  // 103,
        "He4_He4_He4_to_C12",  // 104,
        "C12_C12_to_Mg24_modified",  // 105,
        "O16_O16_to_S32_modified",  // 106,
        "C12_O16_to_Si28_modified"  // 107,
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
