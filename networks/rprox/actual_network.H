#ifndef actual_network_H
#define actual_network_H

#define NEW_NETWORK_IMPLEMENTATION

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>

#include <fundamental_constants.H>
#include <network_properties.H>
#include <microphysics_math.H>
#include <rhs_type.H>
#include <rprox_rates.H>

using namespace amrex;

AMREX_INLINE
void actual_network_init() {}

const std::string network_name = "rprox";

namespace network
{
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> ebin;
}

namespace Rates
{
    enum NetworkRates
    {
        He4_He4_He4_to_C12 = 1,
        C12_H1_to_N13,
        C12_2H1_to_O14,
        N14_H1_to_O15,
        O14_to_N14,
        O14_H1_to_O15,
        O14_He4_to_F17_H1,
        O15_H1_to_C12_He4,
        O15_to_N15,
        O15_He4_3H1_to_Mg22,
        O15_H1_to_O16,
        O15_He4_to_Ne19,
        O16_H1_to_F17,
        O16_He4_to_Ne20,
        O16_He4_2H1_to_Mg22,
        F17_H1_to_Ne18,
        F17_He4_H1_to_Mg22,
        Ne18_to_F18,
        Ne18_He4_to_Na21_H1,
        Ne19_to_F19,
        Ne19_H1_to_Na20,
        Si26_He4_to_P29_H1,
        Mg22_2He4_to_S30,
        Mg22_8H1_to_S30,
        Ti44_He4_to_V47_H1,
        S30_6_5He4_to_Ni56,
        S30_26H1_to_Ni56,
        NumRates = S30_26H1_to_Ni56
    };
};

namespace RHS {

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr rhs_t rhs_data (int rate)
    {
        using namespace Species;
        using namespace Rates;

        rhs_t data;

        switch (rate) {

        case He4_He4_He4_to_C12:
            data.species_A = He4;
            data.species_D = C12;

            data.number_A = 3.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 3;
            data.exponent_D = 1;
            break;

        case C12_H1_to_N13:
            data.species_A = C12;
            data.species_B = H1;
            data.species_D = N13;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case C12_2H1_to_O14:
            data.species_A = C12;
            data.species_B = H1;
            data.species_D = O14;

            data.number_A = 1.0_rt;
            data.number_B = 2.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = C12_H1_to_N13;
            break;

        // Beta decay of O14 -> N14
        case O14_to_N14:
            data.species_A = O14;
            data.species_D = N14;

            data.number_A = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_D = 1;
            break;

        case O14_H1_to_O15:
            data.species_A = O14;
            data.species_B = H1;
            data.species_D = O15;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 0;
            data.exponent_D = 1;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = O14_to_N14;
            break;

        // Proton capture on N14 -> O15
        case N14_H1_to_O15:
            data.species_A = N14;
            data.species_B = H1;
            data.species_D = O15;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case O14_He4_to_F17_H1:
            data.species_A = O14;
            data.species_B = He4;
            data.species_D = F17;
            data.species_E = H1;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        // Beta decay of O15 -> N15; corresponds to former irwk15o
        case O15_to_N15:
            data.species_A = O15;
            data.species_D = N15;

            data.number_A = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_D = 1;
            break;

        // Corresponds to lambda_CNO from WW81 (row 4 of Table 13)
        case O15_H1_to_C12_He4:
            data.species_A = O15;
            data.species_B = H1;
            data.species_D = C12;
            data.species_E = He4;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = N14_H1_to_O15;
            data.additional_reaction_2 = O15_to_N15;
            break;

        // This corresponds to what WW81 call O15 + p. But it is more
        // complicated than that: the rate is constructed constructed from
        // an alpha capture on O15, followed by beta decay from Ne19 to
        // F19, followed by proton capture to Na20, followed by alpha
        // decay to O16. The net result is indeed O15 + p -> O16 catalyzed
        // by an alpha, but in our reactions scheme the effective rate depends
        // only explicitly on y(He4), but not on y(H1) (which is hidden in the
        // reaction rate), so we need to involve helium here as a "reactant."
        case O15_H1_to_O16:
            data.species_A = O15;
            data.species_B = He4;
            data.species_C = H1;
            data.species_D = O16;
            data.species_E = He4;

            data.number_A = 1.0_rt;
            data.number_B = 0.0_rt;
            data.number_C = 1.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 0.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_C = 0;
            data.exponent_D = 1;
            data.exponent_E = 1;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = O15_He4_to_Ne19;
            data.additional_reaction_2 = Ne19_to_F19;
            data.additional_reaction_3 = Ne19_H1_to_Na20;
            break;

        case O15_He4_3H1_to_Mg22:
            data.species_A = O15;
            data.species_B = He4;
            data.species_C = H1;
            data.species_D = Mg22;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_C = 3.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = Ne19_to_F19;
            data.additional_reaction_2 = Ne19_H1_to_Na20;
            data.additional_reaction_3 = O15_He4_to_Ne19;
            break;

        // Alpha capture on O15
        case O15_He4_to_Ne19:
            data.species_A = O15;
            data.species_B = He4;
            data.species_D = Ne19;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        // Proton capture on O16
        case O16_H1_to_F17:
            data.species_A = O16;
            data.species_B = H1;
            data.species_D = F17;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.screen_reverse_reaction = 0; // For consistency with original rprox
            break;

        // Alpha capture on O16
        case O16_He4_to_Ne20:
            data.species_A = O16;
            data.species_B = He4;
            data.species_D = Ne20;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case O16_He4_2H1_to_Mg22:
            data.species_A = O16;
            data.species_B = He4;
            data.species_C = H1;
            data.species_D = Mg22;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_C = 2.0_rt;
            data.number_D = 1.0_rt;

            // It is assumed that the slow rate is the alpha capture,
            // and the proton captures happen instantaneously afterward.
            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_C = 0;
            data.exponent_D = 1;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = O16_He4_to_Ne20;
            break;

        // Proton capture on F17
        case F17_H1_to_Ne18:
            data.species_A = F17;
            data.species_B = H1;
            data.species_D = Ne18;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1.0_rt;
            data.exponent_B = 1.0_rt;
            data.exponent_D = 1.0_rt;
            break;

        case F17_He4_H1_to_Mg22:
            data.species_A = F17;
            data.species_B = He4;
            data.species_C = H1;
            data.species_D = Mg22;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_C = 1.0_rt;
            data.number_D = 1.0_rt;

            // A factor of y(He4) will be added into the rate if we are taking this path
            data.exponent_A = 1;
            data.exponent_B = 0;
            data.exponent_C = 1;
            data.exponent_D = 1;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = F17_H1_to_Ne18;
            data.additional_reaction_2 = Ne18_to_F18;
            data.additional_reaction_3 = Ne18_He4_to_Na21_H1;
            break;

        // Beta decay of Ne18 -> F18; contributes to former irs1
        case Ne18_to_F18:
            data.species_A = Ne18;
            data.species_D = F18;

            data.number_A = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_D = 1;
            break;

        // Alpha capture on Ne18; contributes to former irs1
        case Ne18_He4_to_Na21_H1:
            data.species_A = Ne18;
            data.species_B = He4;
            data.species_D = Na21;
            data.species_E = H1;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        // Beta decay of Ne19 -> F19; contributes to former irr1
        case Ne19_to_F19:
            data.species_A = Ne19;
            data.species_D = F19;

            data.number_A = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_D = 1;
            break;

        // Proton capture on Ne19 -> Na20; contributes to former irr1
        case Ne19_H1_to_Na20:
            data.species_A = Ne19;
            data.species_B = H1;
            data.species_D = Na20;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        // Alpha capture on Si27; contributes to former irlambda1
        case Si26_He4_to_P29_H1:
            data.species_A = Si26;
            data.species_B = He4;
            data.species_D = P29;
            data.species_E = H1;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        // Capture of 2 alphas onto Mg22; corresponds to one branch of the former irlambda1
        case Mg22_2He4_to_S30:
            data.species_A = Mg22;
            data.species_B = He4;
            data.species_D = S30;

            data.number_A = 1.0_rt;
            data.number_B = 2.0_rt;
            data.number_D = 1.0_rt;

            // A factor of y(He4) will be added into the rate if we are taking this path
            data.exponent_A = 1;
            data.exponent_B = 0;
            data.exponent_D = 1;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = Si26_He4_to_P29_H1;
            break;

        // Capture of 8 protons onto Mg22; corresponds to the other branch of the former irlambda1
        case Mg22_8H1_to_S30:
            data.species_A = Mg22;
            data.species_B = H1;
            data.species_D = S30;

            data.number_A = 1.0_rt;
            data.number_B = 8.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 0;
            data.exponent_D = 1;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = Si26_He4_to_P29_H1;
            break;

        // Alpha capture on Ti44; contributes to former irlambda2
        case Ti44_He4_to_V47_H1:
            data.species_A = Ti44;
            data.species_B = He4;
            data.species_D = V47;
            data.species_E = H1;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        // Capture of 6.5 alphas onto S30; corresponds to one branch of the former irlambda2
        case S30_6_5He4_to_Ni56:
            data.species_A = S30;
            data.species_B = He4;
            data.species_D = Ni56;

            data.number_A = 1.0_rt;
            data.number_B = 6.5_rt;
            data.number_D = 1.0_rt;

            // A factor of y(He4) will be added into the rate if we are taking this path
            data.exponent_A = 1;
            data.exponent_B = 0;
            data.exponent_D = 1;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = Ti44_He4_to_V47_H1;
            break;

        // Capture of 26 protons onto S30; corresponds to the other branch of the former irlambda2
        case S30_26H1_to_Ni56:
            data.species_A = S30;
            data.species_B = H1;
            data.species_D = Ni56;

            data.number_A = 1.0_rt;
            data.number_B = 26.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 0;
            data.exponent_D = 1;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = Ti44_He4_to_V47_H1;
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

        if constexpr (rate == He4_He4_He4_to_C12) {
            rate_he4_he4_he4_to_c12(state.tf, rates.fr, rates.frdt);
        }
        else if constexpr (rate == C12_H1_to_N13) {
            rate_p_c12_to_n13(state.tf, rates.fr, rates.frdt);
        }
        else if constexpr (rate == N14_H1_to_O15) {
            rate_p_n14_to_o15(state.tf, rates.fr, rates.frdt);
        }
        else if constexpr (rate == O14_to_N14) {
            rate_o14_to_n14(state.tf, rates.fr, rates.frdt);
        }
        else if constexpr (rate == O14_He4_to_F17_H1) {
            rate_he4_o14_to_p_f17(state.tf, rates.fr, rates.frdt);
        }
        else if constexpr (rate == O15_to_N15) {
            rate_o15_to_n15(state.tf, rates.fr, rates.frdt);
        }
        else if constexpr (rate == O15_He4_to_Ne19) {
            rate_he4_o15_to_ne19(state.tf, rates.fr, rates.frdt);
        }
        else if constexpr (rate == O16_H1_to_F17) {
            rate_p_o16_to_f17(state.tf, rates.fr, rates.frdt);
            rate_f17_to_p_o16(state.tf, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == O16_He4_to_Ne20) {
            rate_he4_o16_to_ne20(state.tf, rates.fr, rates.frdt);
        }
        else if constexpr (rate == F17_H1_to_Ne18) {
            rate_p_f17_to_ne18(state.tf, rates.fr, rates.frdt);
        }
        else if constexpr (rate == Ne18_He4_to_Na21_H1) {
            rate_he4_ne18_to_p_na21(state.tf, rates.fr, rates.frdt);
        }
        else if constexpr (rate == Ne19_to_F19) {
            rate_ne19_to_f19(state.tf, rates.fr, rates.frdt);
        }
        else if constexpr (rate == Ne19_H1_to_Na20) {
            rate_p_ne19_to_na20(state.tf, rates.fr, rates.frdt);
        }
        else if constexpr (rate == Si26_He4_to_P29_H1) {
            rate_he4_si26_to_p_p29(state.tf, rates.fr, rates.frdt);
        }
        else if constexpr (rate == Ti44_He4_to_V47_H1) {
            rate_he4_ti44_to_p_v47(state.tf, rates.fr, rates.frdt);
        }
    }

    template<int rate>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void postprocess_rate ([[maybe_unused]] const rhs_state_t& state, [[maybe_unused]] rate_t& rates,
                           [[maybe_unused]] rate_t& rates1, [[maybe_unused]] rate_t& rates2, [[maybe_unused]] rate_t& rates3)
    {
        using namespace Species;
        using namespace Rates;

        if constexpr (rate == C12_2H1_to_O14) {
            // Rate 1 == C12_H1_to_N13

            rates.fr = rates1.fr;
            rates.frdt = rates1.frdt;
        }
        else if constexpr (rate == O14_H1_to_O15) {
            // Rate 1 == O14_to_N14

            rates.fr = rates1.fr;
            rates.frdt = rates1.frdt;
        }
        else if constexpr (rate == O15_H1_to_C12_He4) {
            // Set as the minimum of the cold CNO and hot CNO cases respectively
            rates.fr = amrex::min(rates1.fr, rates2.fr);

            if (rates.fr == rates2.fr) {
                // For the hot CNO case, divide the rate by y(H1) because the proton capture
                // occurs on the "extra" N15, not on O15 which is represented in the network.
                // There is no T derivative to bring along since this is a weak rate.
                rates.fr /= state.y(H1);
            }
            else {
                // Cold CNO case
                rates.frdt = rates1.frdt;
            }
        }
        else if constexpr (rate == O15_H1_to_O16) {
            // Rate 1 == O15_He4_to_Ne19
            // Rate 2 == Ne19_to_F19
            // Rate 3 == Ne19_H1_to_Na20

            Real rate_irr1 = rates2.fr / (rates2.fr + state.y(H1) * rates3.fr);
            Real dratedt_irr1 = -rate_irr1 * state.y(H1) * rates3.frdt / (rates2.fr + state.y(H1) * rates3.fr);

            rates.fr = rates1.fr * rate_irr1;
            rates.frdt = rates1.frdt * rate_irr1 + rates1.fr * dratedt_irr1;
        }
        else if constexpr (rate == O15_He4_3H1_to_Mg22) {
            // Rate 1 == Ne19_to_F19 (wk19ne in the original rprox)
            // Rate 2 == Ne19_H1_to_Na20
            // Rate 3 == O15_He4_to_Ne19 (ir15ago in the original rprox)

            Real rate_irr1 = rates1.fr / (rates1.fr + state.y(H1) * rates2.fr);
            Real dratedt_irr1 = -rate_irr1 * state.y(H1) * rates2.frdt / (rates1.fr + state.y(H1) * rates2.fr);

            rates.fr = rates3.fr * (1.0_rt - rate_irr1);
            rates.frdt = rates3.frdt * (1.0_rt - rate_irr1) - rates3.fr * dratedt_irr1;
        }
        else if constexpr (rate == O16_He4_2H1_to_Mg22) {
            // Rate 1 == O16_He4_to_Ne20

            rates.fr = rates1.fr;
            rates.frdt = rates1.frdt;
        }
        else if constexpr (rate == F17_He4_H1_to_Mg22) {
            // Rate 1 == F17_H1_to_Ne18
            // Rate 2 == Ne18_to_F18
            // Rate 3 == Ne18_He4_to_Na21_H1

            Real rate_irs1 = rates2.fr / (rates2.fr + state.y(He4) * rates3.fr);
            Real dratedt_irs1 = -rate_irs1 * state.y(He4) * rates3.frdt / (rates2.fr + state.y(He4) * rates3.fr);

            rates.fr = rates1.fr * (1.0_rt - rate_irs1);
            rates.frdt = rates1.frdt * (1.0_rt - rate_irs1) - rates1.fr * dratedt_irs1;
        }
        else if constexpr (rate == Mg22_2He4_to_S30 || rate == Mg22_8H1_to_S30) {
            // Rate 1 == Si26_He4_to_P29_H1

            // Check if this proceeds via p-captures or (a,p) reactions.
            // The Lweak is from WW81, equation C15.

            // These should probably be updated with current rates.
            constexpr Real Lweak = 1.05e0_rt; // this is for NS
            // constexpr Real Lweak = 0.107e0_rt; // this is for lower densities

            rates.fr = amrex::max(Lweak, state.y(He4) * rates1.fr);

            if (rates.fr == Lweak) {
                // Proton capture case
                if constexpr (rate == Mg22_8H1_to_S30) {
                    rates.frdt = 0.0_rt;
                }
                else {
                    rates.fr = 0.0_rt;
                }
            }
            else {
                // Alpha capture case
                if constexpr (rate == Mg22_2He4_to_S30) {
                    rates.frdt = state.y(He4) * rates1.frdt;
                }
                else {
                    rates.fr = 0.0_rt;
                }
            }
        }
        else if constexpr (rate == S30_6_5He4_to_Ni56 || rate == S30_26H1_to_Ni56) {
            // Rate 1 == Ti44_He4_to_V47_H1

            // check if this proceeds via p-captures or (a,p) reactions
            // use 44ti(a,p)v47 as a typical limiting rate for the (a,p) process

            constexpr Real la2 = 1.0_rt / 15.0_rt; // mean rate from 30s to 56ni from p-capture and beta decays

            rates.fr = amrex::max(la2, state.y(He4) * rates1.fr);

            if (rates.fr == la2) {
                // Proton capture case
                if constexpr (rate == S30_26H1_to_Ni56) {
                    rates.frdt = 0.0_rt;
                }
                else {
                    rates.fr = 0.0_rt;
                }
            }
            else {
                // Alpha capture case
                if constexpr (rate == S30_6_5He4_to_Ni56) {
                    rates.frdt = state.y(He4) * rates1.frdt;
                }
                else {
                    rates.fr = 0.0_rt;
                }
            }
        }
    }

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ener_gener_rate ([[maybe_unused]] const rhs_state_t& rhs_state, Real const& dydt)
    {
        //return dydt * network::mion<spec>() * C::Legacy::enuc_conv2;
        return 0.0_rt;
    }

} // namespace RHS

#endif
