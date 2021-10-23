#ifndef _actual_network_H_
#define _actual_network_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

using namespace amrex;

#include <fundamental_constants.H>
#include <network_properties.H>
#include <network_utilities.H>
#include <extern_parameters.H>
#include <tfactors.H>
#include <aprox_rates.H>

void actual_network_init();

const std::string network_name = "iso7";

namespace network
{
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> wion;
}

namespace Rates {
    enum NetworkRatesFR {
        C12_He4_to_O16 = 1,
        He4_He4_He4_to_C12,
        C12_C12_to_Ne20_He4,
        C12_O16_to_Mg24_He4,
        C12_O16_to_Si28,
        O16_O16_to_Si28_He4,
        O16_He4_to_Ne20,
        Ne20_He4_to_Mg24,
        Mg24_He4_to_Si28,
        Ca40_He4_to_Ti44,
        Si28_7He4_to_Ni56,
        NumRatesFR = Si28_7He4_to_Ni56
    };

    enum NetworkRates {C12_He4_to_O16_forward = 1,
                       C12_He4_to_O16_reverse,
                       He4_He4_He4_to_C12_forward,
                       He4_He4_He4_to_C12_reverse,
                       C12_C12_to_Ne20_He4_forward,
                       C12_C12_to_Ne20_He4_reverse,
                       C12_O16_to_Mg24_He4_forward,
                       C12_O16_to_Mg24_He4_reverse,
                       C12_O16_to_Si28_forward,
                       C12_O16_to_Si28_reverse,
                       O16_O16_to_Si28_He4_forward,
                       O16_O16_to_Si28_He4_reverse,
                       O16_He4_to_Ne20_forward,
                       O16_He4_to_Ne20_reverse,
                       Ne20_He4_to_Mg24_forward,
                       Ne20_He4_to_Mg24_reverse,
                       Mg24_He4_to_Si28_forward,
                       Mg24_He4_to_Si28_reverse,
                       Ca40_He4_to_Ti44_forward,
                       Ca40_He4_to_Ti44_reverse,
                       Si28_7He4_to_Ni56_forward,
                       Si28_7He4_to_Ni56_reverse,
                       NumRates = Si28_7He4_to_Ni56_reverse};

    const int NumGroups = 2 + NumSpec;
}

namespace RHS {

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr rhs_t rhs_data (int rate)
    {
        using namespace Species;
        using namespace Rates;

        rhs_t data{};

        data.species_A = -1;
        data.species_B = -1;
        data.species_C = -1;
        data.species_D = -1;
        data.species_E = -1;
        data.species_F = -1;

        data.number_A = 0;
        data.number_B = 0;
        data.number_C = 0;
        data.number_D = 0;
        data.number_E = 0;
        data.number_F = 0;

        data.exponent_A = 0;
        data.exponent_B = 0;
        data.exponent_C = 0;
        data.exponent_D = 0;
        data.exponent_E = 0;
        data.exponent_F = 0;

        data.forward_branching_ratio = 1.0_rt;
        data.reverse_branching_ratio = 1.0_rt;

        data.additional_reaction_1 = -1;
        data.additional_reaction_2 = -1;

        switch (rate) {

        case He4_He4_He4_to_C12:
            data.species_A = He4;
            data.species_D = C12;

            data.number_A = 3;
            data.number_D = 1;

            data.exponent_A = 3;
            data.exponent_D = 1;
            break;

        case C12_He4_to_O16:
            data.species_A = C12;
            data.species_B = He4;
            data.species_D = O16;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case C12_C12_to_Ne20_He4:
            data.species_A = C12;
            data.species_D = Ne20;
            data.species_E = He4;

            data.number_A = 2;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 2;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case C12_O16_to_Mg24_He4:
            data.species_A = C12;
            data.species_B = O16;
            data.species_D = Mg24;
            data.species_E = He4;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;

            // Equal probability branching ratio between (Mg24 + He4) and (Si28) endpoints
            data.forward_branching_ratio = 0.5_rt;
            break;

        case C12_O16_to_Si28:
            data.species_A = C12;
            data.species_B = O16;
            data.species_D = Si28;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            // Equal probability branching ratio between (Mg24 + He4) and (Si28) endpoints
            data.forward_branching_ratio = 0.5_rt;
            break;

        case O16_O16_to_Si28_He4:
            data.species_A = O16;
            data.species_D = Si28;
            data.species_E = He4;

            data.number_A = 2;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 2;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case O16_He4_to_Ne20:
            data.species_A = O16;
            data.species_B = He4;
            data.species_D = Ne20;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Ne20_He4_to_Mg24:
            data.species_A = Ne20;
            data.species_B = He4;
            data.species_D = Mg24;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Mg24_He4_to_Si28:
            data.species_A = Mg24;
            data.species_B = He4;
            data.species_D = Si28;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Ca40_He4_to_Ti44:
            data.species_A = Ca40;
            data.species_B = He4;
            data.species_D = Ti44;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Si28_7He4_to_Ni56:
            data.species_A = Si28;
            data.species_B = He4;
            data.species_D = Ni56;

            data.number_A = 1;
            data.number_B = 7;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.additional_reaction_1 = Ca40_He4_to_Ti44;
            break;

        }

        return data;
    }

    template<int rate>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void evaluate_rate (const tf_t& tf, Real& fr, Real& frdt, Real& rr, Real& rrdt)
    {
        using namespace Rates;

        if constexpr (rate == C12_He4_to_O16) {
            rate_c12ag(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == He4_He4_He4_to_C12) {
            rate_triplealf(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == C12_C12_to_Ne20_He4) {
            rate_c12c12(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == C12_O16_to_Mg24_He4) {
            rate_c12o16(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == C12_O16_to_Si28) {
            rate_c12o16(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == O16_O16_to_Si28_He4) {
            rate_o16o16(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == O16_He4_to_Ne20) {
            rate_o16ag(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Ne20_He4_to_Mg24) {
            rate_ne20ag(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Mg24_He4_to_Si28) {
            rate_mg24ag(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Ca40_He4_to_Ti44) {
            rate_ca40ag(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
    }

    template<int rate>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void postprocess_rate (const Real& rho, const tf_t& tf,
                           const Array1D<Real, 1, NumSpec>& y,
                           Real& fr, Real& frdt, Real& rr, Real& rrdt,
                           Real& fr1, Real& fr1dt, Real& rr1, Real& rr1dt,
                           Real& fr2, Real& fr2dt, Real& rr2, Real& rr2dt)
    {
        using namespace Species;
        using namespace Rates;

        if constexpr (rate == Si28_7He4_to_Ni56) {
            // first rate corresponds to Ca40_He4_to_Ti44

            if (tf.t9 > 2.5_rt && y(C12) + y(O16) <= 4.0e-3_rt) {

                Real t992  = tf.t972 * tf.t9;
                Real t9i92 = 1.0e0_rt / t992;

                Real yeff_ca40   = t9i92 * std::exp(239.42_rt*tf.t9i - 74.741_rt);
                Real yeff_ca40dt = -yeff_ca40*(239.42_rt*tf.t9i2 + 4.5e0_rt*tf.t9i);

                Real yeff_ti44   = t992  * std::exp(-274.12_rt*tf.t9i + 74.914_rt);
                Real yeff_ti44dt = yeff_ti44*(274.12_rt*tf.t9i2 + 4.5e0_rt*tf.t9i);

                Real denom = std::pow(rho * y(He4), 3.0e0_rt);

                fr   = yeff_ca40 * denom * fr1;
                frdt = (yeff_ca40dt * fr1 + yeff_ca40 * fr1dt) * denom * 1.0e-9_rt;

                Real zz = 1.0e0_rt / denom;
                rr = amrex::min(1.0e10_rt, yeff_ti44 * rr1 * zz);

                if (rr == 1.0e10_rt) {
                    rrdt = 0.0e0_rt;
                } else {
                    rrdt  = (yeff_ti44dt * rr1 + yeff_ti44 * rr1dt) * zz * 1.0e-9_rt;
                }
            }
        }
    }
}

#endif
