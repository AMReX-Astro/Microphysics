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
        else if constexpr (spec == O17) {
            return 131.76239561199873_rt;
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
        else if constexpr (spec == Cu59) {
            return 509.8782038560021_rt;
        }
        else if constexpr (spec == Zn60) {
            return 514.9831749199948_rt;
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
        k_p_C12_to_N13 = 1,
        k_He4_C12_to_O16 = 2,
        k_n_N13_to_N14 = 3,
        k_n_O16_to_O17 = 4,
        k_He4_O16_to_Ne20 = 5,
        k_He4_Ne20_to_Mg24 = 6,
        k_p_Na23_to_Mg24 = 7,
        k_He4_Mg24_to_Si28 = 8,
        k_p_Al27_to_Si28 = 9,
        k_He4_Si28_to_S32 = 10,
        k_p_P31_to_S32 = 11,
        k_He4_Cr48_to_Fe52 = 12,
        k_p_nse_Mn51_to_Fe52 = 13,
        k_He4_Mn51_to_Co55 = 14,
        k_He4_Fe52_to_Ni56 = 15,
        k_p_nse_Co55_to_Ni56 = 16,
        k_He4_Co55_to_Cu59 = 17,
        k_He4_Ni56_to_Zn60 = 18,
        k_p_nse_Cu59_to_Zn60 = 19,
        k_C12_C12_to_p_Na23 = 20,
        k_C12_C12_to_He4_Ne20 = 21,
        k_He4_N13_to_p_O16 = 22,
        k_C12_O16_to_p_Al27 = 23,
        k_C12_O16_to_He4_Mg24 = 24,
        k_O16_O16_to_p_P31 = 25,
        k_O16_O16_to_He4_Si28 = 26,
        k_p_O17_to_He4_N14 = 27,
        k_He4_O17_to_n_Ne20 = 28,
        k_p_Na23_to_He4_Ne20 = 29,
        k_p_Al27_to_He4_Mg24 = 30,
        k_p_P31_to_He4_Si28 = 31,
        k_He4_Cr48_to_p_nse_Mn51 = 32,
        k_He4_Fe52_to_p_nse_Co55 = 33,
        k_p_nse_Cu59_to_He4_Ni56 = 34,
        k_He4_He4_He4_to_C12 = 35,
        k_C12_C12_to_Mg24_modified = 36,
        k_O16_O16_to_S32_modified = 37,
        k_C12_O16_to_Si28_modified = 38,
        k_p_nse_Fe54_to_Co55 = 39,
        k_He4_Fe54_to_Ni58 = 40,
        k_p_nse_Fe56_to_Co57 = 41,
        k_n_Co55_to_Co56 = 42,
        k_n_Co56_to_Co57 = 43,
        k_p_nse_Co57_to_Ni58 = 44,
        k_p_nse_Ni58_to_Cu59 = 45,
        k_He4_Mn51_to_p_nse_Fe54 = 46,
        k_He4_Co55_to_p_nse_Ni58 = 47,
        k_n_Co56_to_p_nse_Fe56 = 48,
        k_p_nse_Co57_to_He4_Fe54 = 49,
        k_n_Ni56_to_p_nse_Co56 = 50,
        k_n_Cu59_to_He4_Co56 = 51,
        k_He4_S32_to_Ar36_removed = 52,
        k_p_Cl35_to_Ar36_removed = 53,
        k_p_Cl35_to_He4_S32_removed = 54,
        k_He4_Ar36_to_Ca40_removed = 55,
        k_p_K39_to_Ca40_removed = 56,
        k_p_K39_to_He4_Ar36_removed = 57,
        k_He4_Ca40_to_Ti44_removed = 58,
        k_p_Sc43_to_Ti44_removed = 59,
        k_p_Sc43_to_He4_Ca40_removed = 60,
        k_He4_Ti44_to_Cr48_removed = 61,
        k_He4_Ti44_to_p_V47_removed = 62,
        k_p_V47_to_Cr48_removed = 63,
        k_n_Fe52_to_Fe53_removed = 64,
        k_n_Fe53_to_Fe54_removed = 65,
        k_n_Fe54_to_Fe55_removed = 66,
        k_n_Fe55_to_Fe56_removed = 67,
        k_n_Ni56_to_Ni57_removed = 68,
        k_n_Ni57_to_Ni58_removed = 69,
        k_Co56_to_Fe56 = 70,
        k_Co56_to_Ni56 = 71,
        k_Fe56_to_Co56 = 72,
        k_n_to_p = 73,
        k_Ni56_to_Co56 = 74,
        k_p_to_n = 75,
        k_S32_He4_to_Ar36_approx = 76,
        k_Ar36_to_S32_He4_approx = 77,
        k_Ar36_He4_to_Ca40_approx = 78,
        k_Ca40_to_Ar36_He4_approx = 79,
        k_Ca40_He4_to_Ti44_approx = 80,
        k_Ti44_to_Ca40_He4_approx = 81,
        k_Ti44_He4_to_Cr48_approx = 82,
        k_Cr48_to_Ti44_He4_approx = 83,
        k_Fe52_n_n_to_Fe54_approx = 84,
        k_Fe54_to_Fe52_n_n_approx = 85,
        k_Fe54_n_n_to_Fe56_approx = 86,
        k_Fe56_to_Fe54_n_n_approx = 87,
        k_Ni56_n_n_to_Ni58_approx = 88,
        k_Ni58_to_Ni56_n_n_approx = 89,
        k_N13_to_p_C12_derived = 90,
        k_N14_to_n_N13_derived = 91,
        k_O16_to_He4_C12_derived = 92,
        k_O17_to_n_O16_derived = 93,
        k_Ne20_to_He4_O16_derived = 94,
        k_Mg24_to_p_Na23_derived = 95,
        k_Mg24_to_He4_Ne20_derived = 96,
        k_Si28_to_p_Al27_derived = 97,
        k_Si28_to_He4_Mg24_derived = 98,
        k_S32_to_p_P31_derived = 99,
        k_S32_to_He4_Si28_derived = 100,
        k_Fe52_to_p_nse_Mn51_derived = 101,
        k_Fe52_to_He4_Cr48_derived = 102,
        k_Co55_to_He4_Mn51_derived = 103,
        k_Ni56_to_p_nse_Co55_derived = 104,
        k_Ni56_to_He4_Fe52_derived = 105,
        k_Cu59_to_He4_Co55_derived = 106,
        k_Zn60_to_p_nse_Cu59_derived = 107,
        k_Zn60_to_He4_Ni56_derived = 108,
        k_C12_to_He4_He4_He4_derived = 109,
        k_He4_N14_to_p_O17_derived = 110,
        k_p_O16_to_He4_N13_derived = 111,
        k_n_Ne20_to_He4_O17_derived = 112,
        k_He4_Ne20_to_p_Na23_derived = 113,
        k_He4_Mg24_to_p_Al27_derived = 114,
        k_He4_Si28_to_p_P31_derived = 115,
        k_p_nse_Mn51_to_He4_Cr48_derived = 116,
        k_p_nse_Co55_to_He4_Fe52_derived = 117,
        k_He4_Ni56_to_p_nse_Cu59_derived = 118,
        k_Co55_to_p_nse_Fe54_derived = 119,
        k_Co56_to_n_Co55_derived = 120,
        k_Co57_to_n_Co56_derived = 121,
        k_Co57_to_p_nse_Fe56_derived = 122,
        k_Ni58_to_p_nse_Co57_derived = 123,
        k_Ni58_to_He4_Fe54_derived = 124,
        k_Cu59_to_p_nse_Ni58_derived = 125,
        k_p_nse_Fe54_to_He4_Mn51_derived = 126,
        k_He4_Fe54_to_p_nse_Co57_derived = 127,
        k_p_nse_Fe56_to_n_Co56_derived = 128,
        k_p_nse_Co56_to_n_Ni56_derived = 129,
        k_He4_Co56_to_n_Cu59_derived = 130,
        k_p_nse_Ni58_to_He4_Co55_derived = 131,
        k_He4_S32_to_p_Cl35_derived_removed = 132,
        k_Ar36_to_He4_S32_derived_removed = 133,
        k_Ar36_to_p_Cl35_derived_removed = 134,
        k_He4_Ar36_to_p_K39_derived_removed = 135,
        k_Ca40_to_He4_Ar36_derived_removed = 136,
        k_Ca40_to_p_K39_derived_removed = 137,
        k_He4_Ca40_to_p_Sc43_derived_removed = 138,
        k_Ti44_to_He4_Ca40_derived_removed = 139,
        k_Ti44_to_p_Sc43_derived_removed = 140,
        k_Cr48_to_He4_Ti44_derived_removed = 141,
        k_Cr48_to_p_V47_derived_removed = 142,
        k_p_V47_to_He4_Ti44_derived_removed = 143,
        k_Fe54_to_n_Fe53_derived_removed = 144,
        k_Fe53_to_n_Fe52_derived_removed = 145,
        k_Fe56_to_n_Fe55_derived_removed = 146,
        k_Fe55_to_n_Fe54_derived_removed = 147,
        k_Ni58_to_n_Ni57_derived_removed = 148,
        k_Ni57_to_n_Ni56_derived_removed = 149,
        NumRates = k_Ni57_to_n_Ni56_derived_removed
    };

    // number of reaclib rates

    const int NrateReaclib = 129;

    // number of tabular rates

    const int NrateTabular = 6;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "p_C12_to_N13",  // 1,
        "He4_C12_to_O16",  // 2,
        "n_N13_to_N14",  // 3,
        "n_O16_to_O17",  // 4,
        "He4_O16_to_Ne20",  // 5,
        "He4_Ne20_to_Mg24",  // 6,
        "p_Na23_to_Mg24",  // 7,
        "He4_Mg24_to_Si28",  // 8,
        "p_Al27_to_Si28",  // 9,
        "He4_Si28_to_S32",  // 10,
        "p_P31_to_S32",  // 11,
        "He4_Cr48_to_Fe52",  // 12,
        "p_nse_Mn51_to_Fe52",  // 13,
        "He4_Mn51_to_Co55",  // 14,
        "He4_Fe52_to_Ni56",  // 15,
        "p_nse_Co55_to_Ni56",  // 16,
        "He4_Co55_to_Cu59",  // 17,
        "He4_Ni56_to_Zn60",  // 18,
        "p_nse_Cu59_to_Zn60",  // 19,
        "C12_C12_to_p_Na23",  // 20,
        "C12_C12_to_He4_Ne20",  // 21,
        "He4_N13_to_p_O16",  // 22,
        "C12_O16_to_p_Al27",  // 23,
        "C12_O16_to_He4_Mg24",  // 24,
        "O16_O16_to_p_P31",  // 25,
        "O16_O16_to_He4_Si28",  // 26,
        "p_O17_to_He4_N14",  // 27,
        "He4_O17_to_n_Ne20",  // 28,
        "p_Na23_to_He4_Ne20",  // 29,
        "p_Al27_to_He4_Mg24",  // 30,
        "p_P31_to_He4_Si28",  // 31,
        "He4_Cr48_to_p_nse_Mn51",  // 32,
        "He4_Fe52_to_p_nse_Co55",  // 33,
        "p_nse_Cu59_to_He4_Ni56",  // 34,
        "He4_He4_He4_to_C12",  // 35,
        "C12_C12_to_Mg24_modified",  // 36,
        "O16_O16_to_S32_modified",  // 37,
        "C12_O16_to_Si28_modified",  // 38,
        "p_nse_Fe54_to_Co55",  // 39,
        "He4_Fe54_to_Ni58",  // 40,
        "p_nse_Fe56_to_Co57",  // 41,
        "n_Co55_to_Co56",  // 42,
        "n_Co56_to_Co57",  // 43,
        "p_nse_Co57_to_Ni58",  // 44,
        "p_nse_Ni58_to_Cu59",  // 45,
        "He4_Mn51_to_p_nse_Fe54",  // 46,
        "He4_Co55_to_p_nse_Ni58",  // 47,
        "n_Co56_to_p_nse_Fe56",  // 48,
        "p_nse_Co57_to_He4_Fe54",  // 49,
        "n_Ni56_to_p_nse_Co56",  // 50,
        "n_Cu59_to_He4_Co56",  // 51,
        "He4_S32_to_Ar36_removed",  // 52,
        "p_Cl35_to_Ar36_removed",  // 53,
        "p_Cl35_to_He4_S32_removed",  // 54,
        "He4_Ar36_to_Ca40_removed",  // 55,
        "p_K39_to_Ca40_removed",  // 56,
        "p_K39_to_He4_Ar36_removed",  // 57,
        "He4_Ca40_to_Ti44_removed",  // 58,
        "p_Sc43_to_Ti44_removed",  // 59,
        "p_Sc43_to_He4_Ca40_removed",  // 60,
        "He4_Ti44_to_Cr48_removed",  // 61,
        "He4_Ti44_to_p_V47_removed",  // 62,
        "p_V47_to_Cr48_removed",  // 63,
        "n_Fe52_to_Fe53_removed",  // 64,
        "n_Fe53_to_Fe54_removed",  // 65,
        "n_Fe54_to_Fe55_removed",  // 66,
        "n_Fe55_to_Fe56_removed",  // 67,
        "n_Ni56_to_Ni57_removed",  // 68,
        "n_Ni57_to_Ni58_removed",  // 69,
        "Co56_to_Fe56",  // 70,
        "Co56_to_Ni56",  // 71,
        "Fe56_to_Co56",  // 72,
        "n_to_p",  // 73,
        "Ni56_to_Co56",  // 74,
        "p_to_n",  // 75,
        "S32_He4_to_Ar36_approx",  // 76,
        "Ar36_to_S32_He4_approx",  // 77,
        "Ar36_He4_to_Ca40_approx",  // 78,
        "Ca40_to_Ar36_He4_approx",  // 79,
        "Ca40_He4_to_Ti44_approx",  // 80,
        "Ti44_to_Ca40_He4_approx",  // 81,
        "Ti44_He4_to_Cr48_approx",  // 82,
        "Cr48_to_Ti44_He4_approx",  // 83,
        "Fe52_n_n_to_Fe54_approx",  // 84,
        "Fe54_to_Fe52_n_n_approx",  // 85,
        "Fe54_n_n_to_Fe56_approx",  // 86,
        "Fe56_to_Fe54_n_n_approx",  // 87,
        "Ni56_n_n_to_Ni58_approx",  // 88,
        "Ni58_to_Ni56_n_n_approx",  // 89,
        "N13_to_p_C12_derived",  // 90,
        "N14_to_n_N13_derived",  // 91,
        "O16_to_He4_C12_derived",  // 92,
        "O17_to_n_O16_derived",  // 93,
        "Ne20_to_He4_O16_derived",  // 94,
        "Mg24_to_p_Na23_derived",  // 95,
        "Mg24_to_He4_Ne20_derived",  // 96,
        "Si28_to_p_Al27_derived",  // 97,
        "Si28_to_He4_Mg24_derived",  // 98,
        "S32_to_p_P31_derived",  // 99,
        "S32_to_He4_Si28_derived",  // 100,
        "Fe52_to_p_nse_Mn51_derived",  // 101,
        "Fe52_to_He4_Cr48_derived",  // 102,
        "Co55_to_He4_Mn51_derived",  // 103,
        "Ni56_to_p_nse_Co55_derived",  // 104,
        "Ni56_to_He4_Fe52_derived",  // 105,
        "Cu59_to_He4_Co55_derived",  // 106,
        "Zn60_to_p_nse_Cu59_derived",  // 107,
        "Zn60_to_He4_Ni56_derived",  // 108,
        "C12_to_He4_He4_He4_derived",  // 109,
        "He4_N14_to_p_O17_derived",  // 110,
        "p_O16_to_He4_N13_derived",  // 111,
        "n_Ne20_to_He4_O17_derived",  // 112,
        "He4_Ne20_to_p_Na23_derived",  // 113,
        "He4_Mg24_to_p_Al27_derived",  // 114,
        "He4_Si28_to_p_P31_derived",  // 115,
        "p_nse_Mn51_to_He4_Cr48_derived",  // 116,
        "p_nse_Co55_to_He4_Fe52_derived",  // 117,
        "He4_Ni56_to_p_nse_Cu59_derived",  // 118,
        "Co55_to_p_nse_Fe54_derived",  // 119,
        "Co56_to_n_Co55_derived",  // 120,
        "Co57_to_n_Co56_derived",  // 121,
        "Co57_to_p_nse_Fe56_derived",  // 122,
        "Ni58_to_p_nse_Co57_derived",  // 123,
        "Ni58_to_He4_Fe54_derived",  // 124,
        "Cu59_to_p_nse_Ni58_derived",  // 125,
        "p_nse_Fe54_to_He4_Mn51_derived",  // 126,
        "He4_Fe54_to_p_nse_Co57_derived",  // 127,
        "p_nse_Fe56_to_n_Co56_derived",  // 128,
        "p_nse_Co56_to_n_Ni56_derived",  // 129,
        "He4_Co56_to_n_Cu59_derived",  // 130,
        "p_nse_Ni58_to_He4_Co55_derived",  // 131,
        "He4_S32_to_p_Cl35_derived_removed",  // 132,
        "Ar36_to_He4_S32_derived_removed",  // 133,
        "Ar36_to_p_Cl35_derived_removed",  // 134,
        "He4_Ar36_to_p_K39_derived_removed",  // 135,
        "Ca40_to_He4_Ar36_derived_removed",  // 136,
        "Ca40_to_p_K39_derived_removed",  // 137,
        "He4_Ca40_to_p_Sc43_derived_removed",  // 138,
        "Ti44_to_He4_Ca40_derived_removed",  // 139,
        "Ti44_to_p_Sc43_derived_removed",  // 140,
        "Cr48_to_He4_Ti44_derived_removed",  // 141,
        "Cr48_to_p_V47_derived_removed",  // 142,
        "p_V47_to_He4_Ti44_derived_removed",  // 143,
        "Fe54_to_n_Fe53_derived_removed",  // 144,
        "Fe53_to_n_Fe52_derived_removed",  // 145,
        "Fe56_to_n_Fe55_derived_removed",  // 146,
        "Fe55_to_n_Fe54_derived_removed",  // 147,
        "Ni58_to_n_Ni57_derived_removed",  // 148,
        "Ni57_to_n_Ni56_derived_removed"  // 149,
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
