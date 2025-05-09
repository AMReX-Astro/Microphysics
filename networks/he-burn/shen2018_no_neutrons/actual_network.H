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

        if constexpr (spec == N) {
            return 0.0_rt;
        }
        else if constexpr (spec == H1) {
            return 0.0_rt;
        }
        else if constexpr (spec == He4) {
            return 28.295662457999697_rt;
        }
        else if constexpr (spec == B11) {
            return 76.20505591999972_rt;
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
        else if constexpr (spec == Ne22) {
            return 177.77024384000106_rt;
        }
        else if constexpr (spec == Na23) {
            return 186.56435240400242_rt;
        }
        else if constexpr (spec == Mg24) {
            return 198.2570479679962_rt;
        }
        else if constexpr (spec == Mg25) {
            return 205.5875680680001_rt;
        }
        else if constexpr (spec == Mg26) {
            return 216.68065016799665_rt;
        }
        else if constexpr (spec == Al27) {
            return 224.95193723199915_rt;
        }
        else if constexpr (spec == Si28) {
            return 236.53684539599638_rt;
        }
        else if constexpr (spec == Si29) {
            return 245.01044789599834_rt;
        }
        else if constexpr (spec == Si30) {
            return 255.61964649599756_rt;
        }
        else if constexpr (spec == P30) {
            return 250.60519745999775_rt;
        }
        else if constexpr (spec == P31) {
            return 262.9161999600037_rt;
        }
        else if constexpr (spec == S32) {
            return 271.78016372399725_rt;
        }
        else if constexpr (spec == Cl35) {
            return 298.20976388800045_rt;
        }
        else if constexpr (spec == Ar36) {
            return 306.7167469519991_rt;
        }
        else if constexpr (spec == Ar38) {
            return 327.3426711519933_rt;
        }
        else if constexpr (spec == K39) {
            return 333.7240072160057_rt;
        }
        else if constexpr (spec == Ca40) {
            return 342.05218528000114_rt;
        }
        else if constexpr (spec == Sc43) {
            return 366.82549054399715_rt;
        }
        else if constexpr (spec == Ti44) {
            return 375.47496160800074_rt;
        }
        else if constexpr (spec == V47) {
            return 403.3650388719979_rt;
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
        else if constexpr (spec == Fe56) {
            return 492.2599506639962_rt;
        }
        else if constexpr (spec == Co55) {
            return 476.82912552799826_rt;
        }
        else if constexpr (spec == Ni56) {
            return 483.9956965919919_rt;
        }
        else if constexpr (spec == Ni58) {
            return 506.4596327920008_rt;
        }
        else if constexpr (spec == Ni59) {
            return 515.4588508919987_rt;
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

    enum NetworkRates :
    std::uint8_t
    {
        k_C12_to_p_B11 = 1,
        k_N13_to_p_C12 = 2,
        k_O16_to_He4_C12 = 3,
        k_Ne20_to_He4_O16 = 4,
        k_Na23_to_p_Ne22 = 5,
        k_Mg24_to_p_Na23 = 6,
        k_Mg24_to_He4_Ne20 = 7,
        k_Mg26_to_He4_Ne22 = 8,
        k_Al27_to_p_Mg26 = 9,
        k_Al27_to_He4_Na23 = 10,
        k_Si28_to_p_Al27 = 11,
        k_Si28_to_He4_Mg24 = 12,
        k_Si29_to_He4_Mg25 = 13,
        k_Si30_to_He4_Mg26 = 14,
        k_P30_to_p_Si29 = 15,
        k_P31_to_p_Si30 = 16,
        k_P31_to_He4_Al27 = 17,
        k_S32_to_p_P31 = 18,
        k_S32_to_He4_Si28 = 19,
        k_Cl35_to_He4_P31 = 20,
        k_Ar36_to_p_Cl35 = 21,
        k_Ar36_to_He4_S32 = 22,
        k_K39_to_p_Ar38 = 23,
        k_K39_to_He4_Cl35 = 24,
        k_Ca40_to_p_K39 = 25,
        k_Ca40_to_He4_Ar36 = 26,
        k_Sc43_to_He4_K39 = 27,
        k_Ti44_to_p_Sc43 = 28,
        k_Ti44_to_He4_Ca40 = 29,
        k_V47_to_He4_Sc43 = 30,
        k_Cr48_to_p_V47 = 31,
        k_Cr48_to_He4_Ti44 = 32,
        k_Mn51_to_He4_V47 = 33,
        k_Fe52_to_p_Mn51 = 34,
        k_Fe52_to_He4_Cr48 = 35,
        k_Co55_to_He4_Mn51 = 36,
        k_Ni56_to_p_Co55 = 37,
        k_Ni56_to_He4_Fe52 = 38,
        k_C12_to_He4_He4_He4 = 39,
        k_p_B11_to_C12 = 40,
        k_p_C12_to_N13 = 41,
        k_He4_C12_to_O16 = 42,
        k_He4_O16_to_Ne20 = 43,
        k_He4_Ne20_to_Mg24 = 44,
        k_p_Ne22_to_Na23 = 45,
        k_He4_Ne22_to_Mg26 = 46,
        k_p_Na23_to_Mg24 = 47,
        k_He4_Na23_to_Al27 = 48,
        k_He4_Mg24_to_Si28 = 49,
        k_He4_Mg25_to_Si29 = 50,
        k_p_Mg26_to_Al27 = 51,
        k_He4_Mg26_to_Si30 = 52,
        k_p_Al27_to_Si28 = 53,
        k_He4_Al27_to_P31 = 54,
        k_He4_Si28_to_S32 = 55,
        k_p_Si29_to_P30 = 56,
        k_p_Si30_to_P31 = 57,
        k_p_P31_to_S32 = 58,
        k_He4_P31_to_Cl35 = 59,
        k_He4_S32_to_Ar36 = 60,
        k_p_Cl35_to_Ar36 = 61,
        k_He4_Cl35_to_K39 = 62,
        k_He4_Ar36_to_Ca40 = 63,
        k_p_Ar38_to_K39 = 64,
        k_p_K39_to_Ca40 = 65,
        k_He4_K39_to_Sc43 = 66,
        k_He4_Ca40_to_Ti44 = 67,
        k_p_Sc43_to_Ti44 = 68,
        k_He4_Sc43_to_V47 = 69,
        k_He4_Ti44_to_Cr48 = 70,
        k_p_V47_to_Cr48 = 71,
        k_He4_V47_to_Mn51 = 72,
        k_He4_Cr48_to_Fe52 = 73,
        k_p_Mn51_to_Fe52 = 74,
        k_He4_Mn51_to_Co55 = 75,
        k_He4_Fe52_to_Ni56 = 76,
        k_p_Co55_to_Ni56 = 77,
        k_C12_C12_to_p_Na23 = 78,
        k_C12_C12_to_He4_Ne20 = 79,
        k_He4_N13_to_p_O16 = 80,
        k_p_O16_to_He4_N13 = 81,
        k_C12_O16_to_p_Al27 = 82,
        k_C12_O16_to_He4_Mg24 = 83,
        k_O16_O16_to_p_P31 = 84,
        k_O16_O16_to_He4_Si28 = 85,
        k_He4_Ne20_to_p_Na23 = 86,
        k_He4_Ne20_to_C12_C12 = 87,
        k_C12_Ne20_to_p_P31 = 88,
        k_C12_Ne20_to_He4_Si28 = 89,
        k_p_Na23_to_He4_Ne20 = 90,
        k_p_Na23_to_C12_C12 = 91,
        k_He4_Na23_to_p_Mg26 = 92,
        k_He4_Mg24_to_p_Al27 = 93,
        k_He4_Mg24_to_C12_O16 = 94,
        k_p_Mg26_to_He4_Na23 = 95,
        k_p_Al27_to_He4_Mg24 = 96,
        k_p_Al27_to_C12_O16 = 97,
        k_He4_Al27_to_p_Si30 = 98,
        k_He4_Si28_to_p_P31 = 99,
        k_He4_Si28_to_C12_Ne20 = 100,
        k_He4_Si28_to_O16_O16 = 101,
        k_p_Si30_to_He4_Al27 = 102,
        k_p_P31_to_He4_Si28 = 103,
        k_p_P31_to_C12_Ne20 = 104,
        k_p_P31_to_O16_O16 = 105,
        k_He4_S32_to_p_Cl35 = 106,
        k_p_Cl35_to_He4_S32 = 107,
        k_He4_Cl35_to_p_Ar38 = 108,
        k_He4_Ar36_to_p_K39 = 109,
        k_p_Ar38_to_He4_Cl35 = 110,
        k_p_K39_to_He4_Ar36 = 111,
        k_He4_Ca40_to_p_Sc43 = 112,
        k_p_Sc43_to_He4_Ca40 = 113,
        k_He4_Ti44_to_p_V47 = 114,
        k_p_V47_to_He4_Ti44 = 115,
        k_He4_Cr48_to_p_Mn51 = 116,
        k_p_Mn51_to_He4_Cr48 = 117,
        k_He4_Fe52_to_p_Co55 = 118,
        k_p_Co55_to_He4_Fe52 = 119,
        k_He4_Co55_to_p_Ni58 = 120,
        k_p_Ni58_to_He4_Co55 = 121,
        k_p_B11_to_He4_He4_He4 = 122,
        k_He4_He4_He4_to_C12 = 123,
        k_He4_He4_He4_to_p_B11 = 124,
        k_Ni59_to_n_Ni58 = 125,
        k_n_Ni58_to_Ni59 = 126,
        k_He4_Fe56_to_n_Ni59 = 127,
        k_n_Ni59_to_He4_Fe56 = 128,
        k_P30_to_Si30 = 129,
        k_Si30_to_P30 = 130,
        k_n_to_p = 131,
        k_p_to_n = 132,
        NumRates = k_p_to_n
    };

    // number of reaclib rates

    const int NrateReaclib = 128;

    // number of tabular rates

    const int NrateTabular = 4;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "C12_to_p_B11",  // 1,
        "N13_to_p_C12",  // 2,
        "O16_to_He4_C12",  // 3,
        "Ne20_to_He4_O16",  // 4,
        "Na23_to_p_Ne22",  // 5,
        "Mg24_to_p_Na23",  // 6,
        "Mg24_to_He4_Ne20",  // 7,
        "Mg26_to_He4_Ne22",  // 8,
        "Al27_to_p_Mg26",  // 9,
        "Al27_to_He4_Na23",  // 10,
        "Si28_to_p_Al27",  // 11,
        "Si28_to_He4_Mg24",  // 12,
        "Si29_to_He4_Mg25",  // 13,
        "Si30_to_He4_Mg26",  // 14,
        "P30_to_p_Si29",  // 15,
        "P31_to_p_Si30",  // 16,
        "P31_to_He4_Al27",  // 17,
        "S32_to_p_P31",  // 18,
        "S32_to_He4_Si28",  // 19,
        "Cl35_to_He4_P31",  // 20,
        "Ar36_to_p_Cl35",  // 21,
        "Ar36_to_He4_S32",  // 22,
        "K39_to_p_Ar38",  // 23,
        "K39_to_He4_Cl35",  // 24,
        "Ca40_to_p_K39",  // 25,
        "Ca40_to_He4_Ar36",  // 26,
        "Sc43_to_He4_K39",  // 27,
        "Ti44_to_p_Sc43",  // 28,
        "Ti44_to_He4_Ca40",  // 29,
        "V47_to_He4_Sc43",  // 30,
        "Cr48_to_p_V47",  // 31,
        "Cr48_to_He4_Ti44",  // 32,
        "Mn51_to_He4_V47",  // 33,
        "Fe52_to_p_Mn51",  // 34,
        "Fe52_to_He4_Cr48",  // 35,
        "Co55_to_He4_Mn51",  // 36,
        "Ni56_to_p_Co55",  // 37,
        "Ni56_to_He4_Fe52",  // 38,
        "C12_to_He4_He4_He4",  // 39,
        "p_B11_to_C12",  // 40,
        "p_C12_to_N13",  // 41,
        "He4_C12_to_O16",  // 42,
        "He4_O16_to_Ne20",  // 43,
        "He4_Ne20_to_Mg24",  // 44,
        "p_Ne22_to_Na23",  // 45,
        "He4_Ne22_to_Mg26",  // 46,
        "p_Na23_to_Mg24",  // 47,
        "He4_Na23_to_Al27",  // 48,
        "He4_Mg24_to_Si28",  // 49,
        "He4_Mg25_to_Si29",  // 50,
        "p_Mg26_to_Al27",  // 51,
        "He4_Mg26_to_Si30",  // 52,
        "p_Al27_to_Si28",  // 53,
        "He4_Al27_to_P31",  // 54,
        "He4_Si28_to_S32",  // 55,
        "p_Si29_to_P30",  // 56,
        "p_Si30_to_P31",  // 57,
        "p_P31_to_S32",  // 58,
        "He4_P31_to_Cl35",  // 59,
        "He4_S32_to_Ar36",  // 60,
        "p_Cl35_to_Ar36",  // 61,
        "He4_Cl35_to_K39",  // 62,
        "He4_Ar36_to_Ca40",  // 63,
        "p_Ar38_to_K39",  // 64,
        "p_K39_to_Ca40",  // 65,
        "He4_K39_to_Sc43",  // 66,
        "He4_Ca40_to_Ti44",  // 67,
        "p_Sc43_to_Ti44",  // 68,
        "He4_Sc43_to_V47",  // 69,
        "He4_Ti44_to_Cr48",  // 70,
        "p_V47_to_Cr48",  // 71,
        "He4_V47_to_Mn51",  // 72,
        "He4_Cr48_to_Fe52",  // 73,
        "p_Mn51_to_Fe52",  // 74,
        "He4_Mn51_to_Co55",  // 75,
        "He4_Fe52_to_Ni56",  // 76,
        "p_Co55_to_Ni56",  // 77,
        "C12_C12_to_p_Na23",  // 78,
        "C12_C12_to_He4_Ne20",  // 79,
        "He4_N13_to_p_O16",  // 80,
        "p_O16_to_He4_N13",  // 81,
        "C12_O16_to_p_Al27",  // 82,
        "C12_O16_to_He4_Mg24",  // 83,
        "O16_O16_to_p_P31",  // 84,
        "O16_O16_to_He4_Si28",  // 85,
        "He4_Ne20_to_p_Na23",  // 86,
        "He4_Ne20_to_C12_C12",  // 87,
        "C12_Ne20_to_p_P31",  // 88,
        "C12_Ne20_to_He4_Si28",  // 89,
        "p_Na23_to_He4_Ne20",  // 90,
        "p_Na23_to_C12_C12",  // 91,
        "He4_Na23_to_p_Mg26",  // 92,
        "He4_Mg24_to_p_Al27",  // 93,
        "He4_Mg24_to_C12_O16",  // 94,
        "p_Mg26_to_He4_Na23",  // 95,
        "p_Al27_to_He4_Mg24",  // 96,
        "p_Al27_to_C12_O16",  // 97,
        "He4_Al27_to_p_Si30",  // 98,
        "He4_Si28_to_p_P31",  // 99,
        "He4_Si28_to_C12_Ne20",  // 100,
        "He4_Si28_to_O16_O16",  // 101,
        "p_Si30_to_He4_Al27",  // 102,
        "p_P31_to_He4_Si28",  // 103,
        "p_P31_to_C12_Ne20",  // 104,
        "p_P31_to_O16_O16",  // 105,
        "He4_S32_to_p_Cl35",  // 106,
        "p_Cl35_to_He4_S32",  // 107,
        "He4_Cl35_to_p_Ar38",  // 108,
        "He4_Ar36_to_p_K39",  // 109,
        "p_Ar38_to_He4_Cl35",  // 110,
        "p_K39_to_He4_Ar36",  // 111,
        "He4_Ca40_to_p_Sc43",  // 112,
        "p_Sc43_to_He4_Ca40",  // 113,
        "He4_Ti44_to_p_V47",  // 114,
        "p_V47_to_He4_Ti44",  // 115,
        "He4_Cr48_to_p_Mn51",  // 116,
        "p_Mn51_to_He4_Cr48",  // 117,
        "He4_Fe52_to_p_Co55",  // 118,
        "p_Co55_to_He4_Fe52",  // 119,
        "He4_Co55_to_p_Ni58",  // 120,
        "p_Ni58_to_He4_Co55",  // 121,
        "p_B11_to_He4_He4_He4",  // 122,
        "He4_He4_He4_to_C12",  // 123,
        "He4_He4_He4_to_p_B11",  // 124,
        "Ni59_to_n_Ni58",  // 125,
        "n_Ni58_to_Ni59",  // 126,
        "He4_Fe56_to_n_Ni59",  // 127,
        "n_Ni59_to_He4_Fe56",  // 128,
        "P30_to_Si30",  // 129,
        "Si30_to_P30",  // 130,
        "n_to_p",  // 131,
        "p_to_n"  // 132,
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
