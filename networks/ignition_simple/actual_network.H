#ifndef actual_network_H
#define actual_network_H

#define NEW_NETWORK_IMPLEMENTATION

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

using namespace amrex::literals;

#include <fundamental_constants.H>
#include <network_properties.H>
#include <rhs_type.H>
#include <aprox_rates.H>

AMREX_INLINE
void actual_network_init() {}

const std::string network_name = "ignition_simple";

namespace network
{
    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr amrex::Real bion () {
        using namespace Species;

        static_assert(spec >= 1 && spec <= NumSpec);

        // Set the binding energy of the element
        if constexpr (spec == C12) {
            return 92.16294e0_rt;
        }
        else if constexpr (spec == O16) {
            return 127.62093e0_rt;
        }
        else if constexpr (spec == Mg24) {
            return 198.25790e0_rt;
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

#ifdef NSE_NET
namespace NSE_INDEX
{
    // P_index is for photoionization proton
    constexpr int P_index= -1;
    constexpr int H1_index = -1;
    constexpr int N_index = -1;
    constexpr int He4_index = -1;
}
#endif

namespace Rates
{
    enum NetworkRates
    {
        C12_C12_to_Mg24 = 1,
        NumRates = C12_C12_to_Mg24
    };
}

namespace RHS
{
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr rhs_t rhs_data (int rate)
    {
        using namespace Species;
        using namespace Rates;

        rhs_t data;

        switch (rate) {

        case C12_C12_to_Mg24:
            data.species_A = C12;
            data.species_D = Mg24;

            data.number_A = 2.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 2;
            data.exponent_D = 1;

            // This network doesn't allow the reverse reaction
            data.reverse_branching_ratio = 0.0_rt;
            break;
        }

        return data;
    }

    template<int rate>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void evaluate_analytical_rate (const rhs_state_t& state, rate_t& rates)
    {
        using namespace Species;
        using namespace Rates;

        if constexpr (rate == C12_C12_to_Mg24) {
            rate_c12c12(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
    }

    template<int rate>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void postprocess_rate ([[maybe_unused]] const rhs_state_t& state, [[maybe_unused]] rate_t& rates,
                           [[maybe_unused]] rate_t& rates1, [[maybe_unused]] rate_t& rates2, [[maybe_unused]] rate_t& rates3)
    {
        using namespace Species;
        using namespace Rates;

        // Nothing to do for this network.
    }

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real ener_gener_rate ([[maybe_unused]] const rhs_state_t& rhs_state, amrex::Real const& dydt)
    {
        return dydt * network::mion<spec>() * C::Legacy::enuc_conv2;
    }
} // namespace RHS

#endif
