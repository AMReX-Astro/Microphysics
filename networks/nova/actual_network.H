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
        else if constexpr (spec == F17) {
            return 128.21957600000002_rt;
        }
        else if constexpr (spec == F18) {
            return 137.369484_rt;
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
        k_N13_to_C13_weak_wc12 = 1,
        k_O14_to_N14_weak_wc12 = 2,
        k_O15_to_N15_weak_wc12 = 3,
        k_F17_to_O17_weak_wc12 = 4,
        k_p_C12_to_N13 = 5,
        k_He4_C12_to_O16 = 6,
        k_p_C13_to_N14 = 7,
        k_p_N13_to_O14 = 8,
        k_p_N14_to_O15 = 9,
        k_He4_N14_to_F18 = 10,
        k_p_N15_to_O16 = 11,
        k_p_O16_to_F17 = 12,
        k_p_O17_to_F18 = 13,
        k_He4_N13_to_p_O16 = 14,
        k_p_N15_to_He4_C12 = 15,
        k_He4_O14_to_p_F17 = 16,
        k_p_O17_to_He4_N14 = 17,
        k_p_F18_to_He4_O15 = 18,
        k_He4_He4_He4_to_C12 = 19,
        NumRates = k_He4_He4_He4_to_C12
    };

    // number of reaclib rates

    const int NrateReaclib = 19;

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
        "p_C12_to_N13",  // 5,
        "He4_C12_to_O16",  // 6,
        "p_C13_to_N14",  // 7,
        "p_N13_to_O14",  // 8,
        "p_N14_to_O15",  // 9,
        "He4_N14_to_F18",  // 10,
        "p_N15_to_O16",  // 11,
        "p_O16_to_F17",  // 12,
        "p_O17_to_F18",  // 13,
        "He4_N13_to_p_O16",  // 14,
        "p_N15_to_He4_C12",  // 15,
        "He4_O14_to_p_F17",  // 16,
        "p_O17_to_He4_N14",  // 17,
        "p_F18_to_He4_O15",  // 18,
        "He4_He4_He4_to_C12"  // 19,
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
