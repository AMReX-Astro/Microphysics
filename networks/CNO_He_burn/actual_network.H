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
        else if constexpr (spec == C13) {
            return 97.108037_rt;
        }
        else if constexpr (spec == N13) {
            return 94.105219_rt;
        }
        else if constexpr (spec == N14) {
            return 104.65859599999999_rt;
        }
        else if constexpr (spec == N15) {
            return 115.4919_rt;
        }
        else if constexpr (spec == O14) {
            return 98.731892_rt;
        }
        else if constexpr (spec == O15) {
            return 111.95538_rt;
        }
        else if constexpr (spec == O16) {
            return 127.619296_rt;
        }
        else if constexpr (spec == O17) {
            return 131.76237600000002_rt;
        }
        else if constexpr (spec == O18) {
            return 139.807746_rt;
        }
        else if constexpr (spec == F17) {
            return 128.21957600000002_rt;
        }
        else if constexpr (spec == F18) {
            return 137.369484_rt;
        }
        else if constexpr (spec == F19) {
            return 147.801342_rt;
        }
        else if constexpr (spec == Ne18) {
            return 132.142626_rt;
        }
        else if constexpr (spec == Ne19) {
            return 143.779517_rt;
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
        else if constexpr (spec == Mg22) {
            return 168.58074200000001_rt;
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
        k_N13_to_C13_weak_wc12 = 1,
        k_O14_to_N14_weak_wc12 = 2,
        k_O15_to_N15_weak_wc12 = 3,
        k_F17_to_O17_weak_wc12 = 4,
        k_F18_to_O18_weak_wc12 = 5,
        k_Ne18_to_F18_weak_wc12 = 6,
        k_Ne19_to_F19_weak_wc12 = 7,
        k_Mg22_to_Na22_weak_wc12 = 8,
        k_N13_to_p_C12 = 9,
        k_N14_to_p_C13 = 10,
        k_O14_to_p_N13 = 11,
        k_O15_to_p_N14 = 12,
        k_O16_to_p_N15 = 13,
        k_O16_to_He4_C12 = 14,
        k_F17_to_p_O16 = 15,
        k_F18_to_p_O17 = 16,
        k_F18_to_He4_N14 = 17,
        k_F19_to_p_O18 = 18,
        k_F19_to_He4_N15 = 19,
        k_Ne18_to_p_F17 = 20,
        k_Ne18_to_He4_O14 = 21,
        k_Ne19_to_p_F18 = 22,
        k_Ne19_to_He4_O15 = 23,
        k_Ne20_to_p_F19 = 24,
        k_Ne20_to_He4_O16 = 25,
        k_Ne21_to_He4_O17 = 26,
        k_Na22_to_p_Ne21 = 27,
        k_Na22_to_He4_F18 = 28,
        k_Na23_to_He4_F19 = 29,
        k_Mg22_to_He4_Ne18 = 30,
        k_Mg24_to_p_Na23 = 31,
        k_Mg24_to_He4_Ne20 = 32,
        k_Si28_to_p_Al27 = 33,
        k_Si28_to_He4_Mg24 = 34,
        k_S32_to_p_P31 = 35,
        k_S32_to_He4_Si28 = 36,
        k_C12_to_He4_He4_He4 = 37,
        k_p_C12_to_N13 = 38,
        k_He4_C12_to_O16 = 39,
        k_p_C13_to_N14 = 40,
        k_p_N13_to_O14 = 41,
        k_p_N14_to_O15 = 42,
        k_He4_N14_to_F18 = 43,
        k_p_N15_to_O16 = 44,
        k_He4_N15_to_F19 = 45,
        k_He4_O14_to_Ne18 = 46,
        k_He4_O15_to_Ne19 = 47,
        k_p_O16_to_F17 = 48,
        k_He4_O16_to_Ne20 = 49,
        k_p_O17_to_F18 = 50,
        k_He4_O17_to_Ne21 = 51,
        k_p_O18_to_F19 = 52,
        k_p_F17_to_Ne18 = 53,
        k_p_F18_to_Ne19 = 54,
        k_He4_F18_to_Na22 = 55,
        k_p_F19_to_Ne20 = 56,
        k_He4_F19_to_Na23 = 57,
        k_He4_Ne18_to_Mg22 = 58,
        k_He4_Ne20_to_Mg24 = 59,
        k_p_Ne21_to_Na22 = 60,
        k_p_Na23_to_Mg24 = 61,
        k_He4_Mg24_to_Si28 = 62,
        k_p_Al27_to_Si28 = 63,
        k_He4_Si28_to_S32 = 64,
        k_p_P31_to_S32 = 65,
        k_He4_C12_to_p_N15 = 66,
        k_C12_C12_to_p_Na23 = 67,
        k_C12_C12_to_He4_Ne20 = 68,
        k_He4_N13_to_p_O16 = 69,
        k_He4_N14_to_p_O17 = 70,
        k_p_N15_to_He4_C12 = 71,
        k_He4_N15_to_p_O18 = 72,
        k_He4_O14_to_p_F17 = 73,
        k_He4_O15_to_p_F18 = 74,
        k_p_O16_to_He4_N13 = 75,
        k_He4_O16_to_p_F19 = 76,
        k_C12_O16_to_p_Al27 = 77,
        k_C12_O16_to_He4_Mg24 = 78,
        k_O16_O16_to_p_P31 = 79,
        k_O16_O16_to_He4_Si28 = 80,
        k_p_O17_to_He4_N14 = 81,
        k_p_O18_to_He4_N15 = 82,
        k_p_F17_to_He4_O14 = 83,
        k_He4_F17_to_p_Ne20 = 84,
        k_p_F18_to_He4_O15 = 85,
        k_He4_F18_to_p_Ne21 = 86,
        k_p_F19_to_He4_O16 = 87,
        k_He4_Ne19_to_p_Na22 = 88,
        k_p_Ne20_to_He4_F17 = 89,
        k_He4_Ne20_to_p_Na23 = 90,
        k_p_Ne21_to_He4_F18 = 91,
        k_p_Na22_to_He4_Ne19 = 92,
        k_p_Na23_to_He4_Ne20 = 93,
        k_He4_Mg24_to_p_Al27 = 94,
        k_p_Al27_to_He4_Mg24 = 95,
        k_He4_Si28_to_p_P31 = 96,
        k_p_P31_to_He4_Si28 = 97,
        k_He4_He4_He4_to_C12 = 98,
        k_C12_C12_to_Mg24_modified = 99,
        k_O16_O16_to_S32_modified = 100,
        k_C12_O16_to_Si28_modified = 101,
        k_He4_S32_to_Ar36_removed = 102,
        k_He4_S32_to_p_Cl35_removed = 103,
        k_p_Cl35_to_Ar36_removed = 104,
        k_Ar36_to_He4_S32_removed = 105,
        k_Ar36_to_p_Cl35_removed = 106,
        k_p_Cl35_to_He4_S32_removed = 107,
        k_He4_Ar36_to_Ca40_removed = 108,
        k_He4_Ar36_to_p_K39_removed = 109,
        k_p_K39_to_Ca40_removed = 110,
        k_Ca40_to_He4_Ar36_removed = 111,
        k_Ca40_to_p_K39_removed = 112,
        k_p_K39_to_He4_Ar36_removed = 113,
        k_He4_Ca40_to_Ti44_removed = 114,
        k_He4_Ca40_to_p_Sc43_removed = 115,
        k_p_Sc43_to_Ti44_removed = 116,
        k_Ti44_to_He4_Ca40_removed = 117,
        k_Ti44_to_p_Sc43_removed = 118,
        k_p_Sc43_to_He4_Ca40_removed = 119,
        k_He4_Ti44_to_Cr48_removed = 120,
        k_He4_Ti44_to_p_V47_removed = 121,
        k_p_V47_to_Cr48_removed = 122,
        k_Cr48_to_He4_Ti44_removed = 123,
        k_Cr48_to_p_V47_removed = 124,
        k_p_V47_to_He4_Ti44_removed = 125,
        k_He4_Cr48_to_Fe52_removed = 126,
        k_He4_Cr48_to_p_Mn51_removed = 127,
        k_p_Mn51_to_Fe52_removed = 128,
        k_Fe52_to_He4_Cr48_removed = 129,
        k_Fe52_to_p_Mn51_removed = 130,
        k_p_Mn51_to_He4_Cr48_removed = 131,
        k_He4_Fe52_to_Ni56_removed = 132,
        k_He4_Fe52_to_p_Co55_removed = 133,
        k_p_Co55_to_Ni56_removed = 134,
        k_Ni56_to_He4_Fe52_removed = 135,
        k_Ni56_to_p_Co55_removed = 136,
        k_p_Co55_to_He4_Fe52_removed = 137,
        k_S32_He4_to_Ar36_approx = 138,
        k_Ar36_to_S32_He4_approx = 139,
        k_Ar36_He4_to_Ca40_approx = 140,
        k_Ca40_to_Ar36_He4_approx = 141,
        k_Ca40_He4_to_Ti44_approx = 142,
        k_Ti44_to_Ca40_He4_approx = 143,
        k_Ti44_He4_to_Cr48_approx = 144,
        k_Cr48_to_Ti44_He4_approx = 145,
        k_Cr48_He4_to_Fe52_approx = 146,
        k_Fe52_to_Cr48_He4_approx = 147,
        k_Fe52_He4_to_Ni56_approx = 148,
        k_Ni56_to_Fe52_He4_approx = 149,
        NumRates = k_Ni56_to_Fe52_He4_approx
    };

    // number of reaclib rates

    const int NrateReaclib = 137;

    // number of tabular rates

    const int NrateTabular = 0;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "N13_to_C13_weak_wc12",  // 1,
        "O14_to_N14_weak_wc12",  // 2,
        "O15_to_N15_weak_wc12",  // 3,
        "F17_to_O17_weak_wc12",  // 4,
        "F18_to_O18_weak_wc12",  // 5,
        "Ne18_to_F18_weak_wc12",  // 6,
        "Ne19_to_F19_weak_wc12",  // 7,
        "Mg22_to_Na22_weak_wc12",  // 8,
        "N13_to_p_C12",  // 9,
        "N14_to_p_C13",  // 10,
        "O14_to_p_N13",  // 11,
        "O15_to_p_N14",  // 12,
        "O16_to_p_N15",  // 13,
        "O16_to_He4_C12",  // 14,
        "F17_to_p_O16",  // 15,
        "F18_to_p_O17",  // 16,
        "F18_to_He4_N14",  // 17,
        "F19_to_p_O18",  // 18,
        "F19_to_He4_N15",  // 19,
        "Ne18_to_p_F17",  // 20,
        "Ne18_to_He4_O14",  // 21,
        "Ne19_to_p_F18",  // 22,
        "Ne19_to_He4_O15",  // 23,
        "Ne20_to_p_F19",  // 24,
        "Ne20_to_He4_O16",  // 25,
        "Ne21_to_He4_O17",  // 26,
        "Na22_to_p_Ne21",  // 27,
        "Na22_to_He4_F18",  // 28,
        "Na23_to_He4_F19",  // 29,
        "Mg22_to_He4_Ne18",  // 30,
        "Mg24_to_p_Na23",  // 31,
        "Mg24_to_He4_Ne20",  // 32,
        "Si28_to_p_Al27",  // 33,
        "Si28_to_He4_Mg24",  // 34,
        "S32_to_p_P31",  // 35,
        "S32_to_He4_Si28",  // 36,
        "C12_to_He4_He4_He4",  // 37,
        "p_C12_to_N13",  // 38,
        "He4_C12_to_O16",  // 39,
        "p_C13_to_N14",  // 40,
        "p_N13_to_O14",  // 41,
        "p_N14_to_O15",  // 42,
        "He4_N14_to_F18",  // 43,
        "p_N15_to_O16",  // 44,
        "He4_N15_to_F19",  // 45,
        "He4_O14_to_Ne18",  // 46,
        "He4_O15_to_Ne19",  // 47,
        "p_O16_to_F17",  // 48,
        "He4_O16_to_Ne20",  // 49,
        "p_O17_to_F18",  // 50,
        "He4_O17_to_Ne21",  // 51,
        "p_O18_to_F19",  // 52,
        "p_F17_to_Ne18",  // 53,
        "p_F18_to_Ne19",  // 54,
        "He4_F18_to_Na22",  // 55,
        "p_F19_to_Ne20",  // 56,
        "He4_F19_to_Na23",  // 57,
        "He4_Ne18_to_Mg22",  // 58,
        "He4_Ne20_to_Mg24",  // 59,
        "p_Ne21_to_Na22",  // 60,
        "p_Na23_to_Mg24",  // 61,
        "He4_Mg24_to_Si28",  // 62,
        "p_Al27_to_Si28",  // 63,
        "He4_Si28_to_S32",  // 64,
        "p_P31_to_S32",  // 65,
        "He4_C12_to_p_N15",  // 66,
        "C12_C12_to_p_Na23",  // 67,
        "C12_C12_to_He4_Ne20",  // 68,
        "He4_N13_to_p_O16",  // 69,
        "He4_N14_to_p_O17",  // 70,
        "p_N15_to_He4_C12",  // 71,
        "He4_N15_to_p_O18",  // 72,
        "He4_O14_to_p_F17",  // 73,
        "He4_O15_to_p_F18",  // 74,
        "p_O16_to_He4_N13",  // 75,
        "He4_O16_to_p_F19",  // 76,
        "C12_O16_to_p_Al27",  // 77,
        "C12_O16_to_He4_Mg24",  // 78,
        "O16_O16_to_p_P31",  // 79,
        "O16_O16_to_He4_Si28",  // 80,
        "p_O17_to_He4_N14",  // 81,
        "p_O18_to_He4_N15",  // 82,
        "p_F17_to_He4_O14",  // 83,
        "He4_F17_to_p_Ne20",  // 84,
        "p_F18_to_He4_O15",  // 85,
        "He4_F18_to_p_Ne21",  // 86,
        "p_F19_to_He4_O16",  // 87,
        "He4_Ne19_to_p_Na22",  // 88,
        "p_Ne20_to_He4_F17",  // 89,
        "He4_Ne20_to_p_Na23",  // 90,
        "p_Ne21_to_He4_F18",  // 91,
        "p_Na22_to_He4_Ne19",  // 92,
        "p_Na23_to_He4_Ne20",  // 93,
        "He4_Mg24_to_p_Al27",  // 94,
        "p_Al27_to_He4_Mg24",  // 95,
        "He4_Si28_to_p_P31",  // 96,
        "p_P31_to_He4_Si28",  // 97,
        "He4_He4_He4_to_C12",  // 98,
        "C12_C12_to_Mg24_modified",  // 99,
        "O16_O16_to_S32_modified",  // 100,
        "C12_O16_to_Si28_modified",  // 101,
        "He4_S32_to_Ar36_removed",  // 102,
        "He4_S32_to_p_Cl35_removed",  // 103,
        "p_Cl35_to_Ar36_removed",  // 104,
        "Ar36_to_He4_S32_removed",  // 105,
        "Ar36_to_p_Cl35_removed",  // 106,
        "p_Cl35_to_He4_S32_removed",  // 107,
        "He4_Ar36_to_Ca40_removed",  // 108,
        "He4_Ar36_to_p_K39_removed",  // 109,
        "p_K39_to_Ca40_removed",  // 110,
        "Ca40_to_He4_Ar36_removed",  // 111,
        "Ca40_to_p_K39_removed",  // 112,
        "p_K39_to_He4_Ar36_removed",  // 113,
        "He4_Ca40_to_Ti44_removed",  // 114,
        "He4_Ca40_to_p_Sc43_removed",  // 115,
        "p_Sc43_to_Ti44_removed",  // 116,
        "Ti44_to_He4_Ca40_removed",  // 117,
        "Ti44_to_p_Sc43_removed",  // 118,
        "p_Sc43_to_He4_Ca40_removed",  // 119,
        "He4_Ti44_to_Cr48_removed",  // 120,
        "He4_Ti44_to_p_V47_removed",  // 121,
        "p_V47_to_Cr48_removed",  // 122,
        "Cr48_to_He4_Ti44_removed",  // 123,
        "Cr48_to_p_V47_removed",  // 124,
        "p_V47_to_He4_Ti44_removed",  // 125,
        "He4_Cr48_to_Fe52_removed",  // 126,
        "He4_Cr48_to_p_Mn51_removed",  // 127,
        "p_Mn51_to_Fe52_removed",  // 128,
        "Fe52_to_He4_Cr48_removed",  // 129,
        "Fe52_to_p_Mn51_removed",  // 130,
        "p_Mn51_to_He4_Cr48_removed",  // 131,
        "He4_Fe52_to_Ni56_removed",  // 132,
        "He4_Fe52_to_p_Co55_removed",  // 133,
        "p_Co55_to_Ni56_removed",  // 134,
        "Ni56_to_He4_Fe52_removed",  // 135,
        "Ni56_to_p_Co55_removed",  // 136,
        "p_Co55_to_He4_Fe52_removed",  // 137,
        "S32_He4_to_Ar36_approx",  // 138,
        "Ar36_to_S32_He4_approx",  // 139,
        "Ar36_He4_to_Ca40_approx",  // 140,
        "Ca40_to_Ar36_He4_approx",  // 141,
        "Ca40_He4_to_Ti44_approx",  // 142,
        "Ti44_to_Ca40_He4_approx",  // 143,
        "Ti44_He4_to_Cr48_approx",  // 144,
        "Cr48_to_Ti44_He4_approx",  // 145,
        "Cr48_He4_to_Fe52_approx",  // 146,
        "Fe52_to_Cr48_He4_approx",  // 147,
        "Fe52_He4_to_Ni56_approx",  // 148,
        "Ni56_to_Fe52_He4_approx"  // 149,
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
