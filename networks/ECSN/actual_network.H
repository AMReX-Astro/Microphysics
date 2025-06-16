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
            return 28.295662457999697_rt;
        }
        else if constexpr (spec == O16) {
            return 127.6193154119992_rt;
        }
        else if constexpr (spec == O20) {
            return 151.37138571200194_rt;
        }
        else if constexpr (spec == F20) {
            return 154.40270167600102_rt;
        }
        else if constexpr (spec == Ne20) {
            return 160.64482384000075_rt;
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
        else if constexpr (spec == O16) {
            return 2.6560180592333686e-23_rt;
        }
        else if constexpr (spec == O20) {
            return 3.3217548673186867e-23_rt;
        }
        else if constexpr (spec == F20) {
            return 3.321075020137486e-23_rt;
        }
        else if constexpr (spec == Ne20) {
            return 3.3198227947612416e-23_rt;
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


        // Return zero if we don't recognize the species.
        return 0.0_rt;
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
        k_Ne20_to_He4_O16 = 1,
        k_He4_O16_to_Ne20 = 2,
        k_He4_Ne20_to_Mg24 = 3,
        k_He4_Mg24_to_Si28 = 4,
        k_p_Al27_to_Si28 = 5,
        k_He4_Al27_to_P31 = 6,
        k_He4_Si28_to_S32 = 7,
        k_p_P31_to_S32 = 8,
        k_O16_O16_to_p_P31 = 9,
        k_O16_O16_to_He4_Si28 = 10,
        k_He4_Mg24_to_p_Al27 = 11,
        k_p_Al27_to_He4_Mg24 = 12,
        k_He4_Si28_to_p_P31 = 13,
        k_p_P31_to_He4_Si28 = 14,
        k_F20_to_Ne20 = 15,
        k_F20_to_O20 = 16,
        k_Ne20_to_F20 = 17,
        k_O20_to_F20 = 18,
        NumRates = k_O20_to_F20
    };

    // number of reaclib rates

    const int NrateReaclib = 14;

    // number of tabular rates

    const int NrateTabular = 4;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "Ne20_to_He4_O16",  // 1,
        "He4_O16_to_Ne20",  // 2,
        "He4_Ne20_to_Mg24",  // 3,
        "He4_Mg24_to_Si28",  // 4,
        "p_Al27_to_Si28",  // 5,
        "He4_Al27_to_P31",  // 6,
        "He4_Si28_to_S32",  // 7,
        "p_P31_to_S32",  // 8,
        "O16_O16_to_p_P31",  // 9,
        "O16_O16_to_He4_Si28",  // 10,
        "He4_Mg24_to_p_Al27",  // 11,
        "p_Al27_to_He4_Mg24",  // 12,
        "He4_Si28_to_p_P31",  // 13,
        "p_P31_to_He4_Si28",  // 14,
        "F20_to_Ne20",  // 15,
        "F20_to_O20",  // 16,
        "Ne20_to_F20",  // 17,
        "O20_to_F20"  // 18,
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
