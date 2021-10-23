#ifndef _actual_rhs_H_
#define _actual_rhs_H_

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include <extern_parameters.H>
#include <microphysics_math.H>
#include <actual_network.H>
#include <burn_type.H>
#include <tfactors.H>
#include <rate_type.H>
#include <screen.H>
#include <sneut5.H>
#include <aprox_rates.H>
#include <ArrayUtilities.H>
#include <rhs_utilities.H>
#include <jacobian_utilities.H>

using namespace amrex;
using namespace ArrayUtil;

void actual_rhs_init();

namespace RateTable
{
    constexpr Real tab_tlo = 6.0e0_rt;
    constexpr Real tab_thi = 10.0e0_rt;
    constexpr int tab_per_decade = 2000;
    constexpr int nrattab = static_cast<int>(tab_thi - tab_tlo) * tab_per_decade + 1;
    constexpr int tab_imax = static_cast<int>(tab_thi - tab_tlo) * tab_per_decade + 1;
    constexpr Real tab_tstp = (tab_thi - tab_tlo) / static_cast<Real>(tab_imax - 1);

    extern AMREX_GPU_MANAGED Array3D<Real, 1, Rates::NumRatesFR, 1, 2, 1, nrattab> rattab;
    extern AMREX_GPU_MANAGED Array3D<Real, 1, Rates::NumRatesFR, 1, 2, 1, nrattab> drattabdt;
    extern AMREX_GPU_MANAGED Array1D<Real, 1, nrattab> ttab;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void set_iso7rat()
{
    using namespace Rates;
    using namespace RateTable;

    for (int i = 1; i <= tab_imax; ++i) {
        Real temp = tab_tlo + static_cast<Real>(i-1) * tab_tstp;
        temp = std::pow(10.0e0_rt, temp);

        ttab(i) = temp;

        // Get the temperature factors
        tf_t tf = get_tfactors(temp);

        RHS::constexpr_for<1, NumRatesFR+1>([&] (auto n)
        {
            constexpr int rate = n;

            Real fr = 0.0_rt;
            Real frdt = 0.0_rt;

            Real rr = 0.0_rt;
            Real rrdt = 0.0_rt;

            RHS::evaluate_rate<rate>(tf, fr, frdt, rr, rrdt);

            rattab(rate, 1, i)    = fr;
            rattab(rate, 2, i)    = rr;
            drattabdt(rate, 1, i) = frdt;
            drattabdt(rate, 2, i) = rrdt;
       });
    }
}

template<int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void construct_rate (const Real& rho, const tf_t& tf, const Array1D<Real, 1, NumSpec>& y,
                     const plasma_state_t& pstate,
                     const Real& alfa, const Real& beta, const Real& gama, const Real& delt, const int& iat,
                     Real& fr, Real& frdt, Real& rr, Real& rrdt)
{
    using namespace Species;
    using namespace Rates;
    using namespace RateTable;

    if (use_tables) {
        fr = (alfa * rattab(rate, 1, iat  ) +
              beta * rattab(rate, 1, iat+1) +
              gama * rattab(rate, 1, iat+2) +
              delt * rattab(rate, 1, iat+3));

        rr = (alfa * rattab(rate, 2, iat  ) +
              beta * rattab(rate, 2, iat+1) +
              gama * rattab(rate, 2, iat+2) +
              delt * rattab(rate, 2, iat+3));

        frdt = (alfa * drattabdt(rate, 1, iat  ) +
                beta * drattabdt(rate, 1, iat+1) +
                gama * drattabdt(rate, 1, iat+2) +
                delt * drattabdt(rate, 1, iat+3));

        rrdt = (alfa * drattabdt(rate, 2, iat  ) +
                beta * drattabdt(rate, 2, iat+1) +
                gama * drattabdt(rate, 2, iat+2) +
                delt * drattabdt(rate, 2, iat+3));
    }
    else {
        RHS::evaluate_rate<rate>(tf, fr, frdt, rr, rrdt);
    }

    // Set the density dependence

    int forward_exponent = RHS::density_exponent_forward<rate>();
    int reverse_exponent = RHS::density_exponent_reverse<rate>();

    Real density_term_forward = std::pow(rho, forward_exponent);
    Real density_term_reverse = std::pow(rho, reverse_exponent);

    fr *= density_term_forward;
    frdt *= density_term_forward;
    rr *= density_term_reverse;
    rrdt *= density_term_reverse;

    // The screening behavior depends on the type of reaction. We provide screening
    // here for the reaction classes we know about, and any other reactions are unscreened.

    constexpr rhs_t rhs_data = RHS::rhs_data(rate);

    if constexpr (rhs_data.number_A == 1 && rhs_data.number_B == 1 && rhs_data.number_C == 0) {
        // Forward reaction is A + B, screen using these two species

        constexpr Real Z1 = NetworkProperties::zion(rhs_data.species_A);
        constexpr Real A1 = NetworkProperties::aion(rhs_data.species_A);

        constexpr Real Z2 = NetworkProperties::zion(rhs_data.species_B);
        constexpr Real A2 = NetworkProperties::aion(rhs_data.species_B);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        Real sc, scdt;
        actual_screen5(pstate, scn_fac, sc, scdt);

        frdt = frdt * sc + fr * scdt;
        fr   = fr * sc;
    }

    if constexpr (rhs_data.number_A == 2 && rhs_data.number_B == 0 && rhs_data.number_C == 0) {
        // Forward reaction is A + A, screen using just this species

        constexpr Real Z1 = NetworkProperties::zion(rhs_data.species_A);
        constexpr Real A1 = NetworkProperties::aion(rhs_data.species_A);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z1, A1);

        Real sc, scdt;
        actual_screen5(pstate, scn_fac, sc, scdt);

        frdt = frdt * sc + fr * scdt;
        fr   = fr * sc;
    }

    if constexpr (rhs_data.number_A == 3 && rhs_data.number_B == 0 && rhs_data.number_C == 0) {
        // Forward reaction is triple alpha or an equivalent, screen using A + A
        // and then A + X where X has twice the number of protons and neutrons.

        constexpr Real Z1 = NetworkProperties::zion(rhs_data.species_A);
        constexpr Real A1 = NetworkProperties::aion(rhs_data.species_A);

        constexpr auto scn_fac1 = scrn::calculate_screen_factor(Z1, A1, Z1, A1);

        Real sc1, sc1dt;
        actual_screen5(pstate, scn_fac1, sc1, sc1dt);

        constexpr Real Z2 = 2.0_rt * Z1;
        constexpr Real A2 = 2.0_rt * A1;

        constexpr auto scn_fac2 = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        Real sc2, sc2dt;
        actual_screen5(pstate, scn_fac2, sc2, sc2dt);

        // Compute combined screening factor

        Real sc = sc1 * sc2;
        Real scdt = sc1dt * sc2 + sc1 * sc2dt;

        frdt = frdt * sc + fr * scdt;
        fr   = fr * sc;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_rates(burn_t const& state, rate_fr_t& r, rate_fr_t& drdt)
{
    using namespace Species;
    using namespace Rates;
    using namespace RateTable;

    Real rho, temp;
    Array1D<Real, 1, NumSpec> y;

    // Get the data from the state
    rho  = state.rho;
    temp = state.T;

    for (int i = 1; i <= NumSpec; ++i)
        y(i) = state.xn[i-1] * aion_inv[i-1];


    // Set up the state data, which is the same for all screening factors.
    plasma_state_t pstate;
    fill_plasma_state(pstate, temp, rho, y);

    // Get the temperature factors
    tf_t tf = get_tfactors(temp);

    Real alfa, beta, gama, delt;
    int iat;

    if (use_tables) {
        // hash locate
        constexpr int mp = 4;
        iat = static_cast<int>((std::log10(temp) - tab_tlo) / tab_tstp) + 1;
        iat = amrex::max(1, amrex::min(iat - 1, tab_imax - mp + 1));

        // setup the lagrange interpolation coefficients for a cubic
        Real x  = temp;
        Real x1 = ttab(iat);
        Real x2 = ttab(iat+1);
        Real x3 = ttab(iat+2);
        Real x4 = ttab(iat+3);
        Real a  = x - x1;
        Real b  = x - x2;
        Real c  = x - x3;
        Real d  = x - x4;
        Real e  = x1 - x2;
        Real f  = x1 - x3;
        Real g  = x1 - x4;
        Real h  = x2 - x3;
        Real p  = x2 - x4;
        Real q  = x3 - x4;
        alfa =  b * c * d / (e * f * g);
        beta = -a * c * d / (e * h * p);
        gama =  a * b * d / (f * h * q);
        delt = -a * b * c / (g * p * q);
    }

    RHS::constexpr_for<1, Rates::NumRatesFR+1>([&] (auto n)
    {
        constexpr int rate = n;

        constexpr rhs_t rhs_data = RHS::rhs_data(rate);

        // If this is an intermediate reaction for another reaction,
        // skip it now. We recognize intermediate reactions because
        // at least one of the species will not be a "real" species.

        if constexpr (rhs_data.species_A > NumSpec ||
                      rhs_data.species_B > NumSpec ||
                      rhs_data.species_C > NumSpec ||
                      rhs_data.species_D > NumSpec ||
                      rhs_data.species_E > NumSpec ||
                      rhs_data.species_F > NumSpec) {
            return;
        }

        // Construct the (screened and density-scaled) rate.

        Real fr = 0.0, frdt = 0.0, rr = 0.0, rrdt = 0.0;

        construct_rate<rate>(rho, tf, y, pstate,
                             alfa, beta, gama, delt, iat,
                             fr, frdt, rr, rrdt);

        // Construct all intermediate rates needed to augment this reaction.
        // To keep the problem bounded we assume that there are no more than
        // two intermediate reactions needed, however it is straightforward
        // to add more as needed (or to generalize this into a loop over all
        // rates if really needed). If any are present, postprocess the rate
        // using the additional terms.

        Real fr1 = 0.0, fr1dt = 0.0, rr1 = 0.0, rr1dt = 0.0;
        Real fr2 = 0.0, fr2dt = 0.0, rr2 = 0.0, rr2dt = 0.0;

        constexpr int rate1 = rhs_data.additional_reaction_1;
        constexpr int rate2 = rhs_data.additional_reaction_2;

        if constexpr (rate1 >= 0 || rate2 >= 0) {
            if constexpr (rate1 >= 0) {
                construct_rate<rate1>(rho, tf, y, pstate,
                                      alfa, beta, gama, delt, iat,
                                      fr1, fr1dt, rr1, rr1dt);
            }

            if constexpr (rate2 >= 0) {
                construct_rate<rate1>(rho, tf, y, pstate,
                                      alfa, beta, gama, delt, iat,
                                      fr2, fr2dt, rr2, rr2dt);
            }

            RHS::postprocess_rate<rate>(rho, tf, y,
                                        fr, frdt, rr, rrdt,
                                        fr1, fr1dt, rr1, rr1dt,
                                        fr2, fr2dt, rr2, rr2dt);
        }

        // Handle any branching ratios.

        fr *= rhs_data.forward_branching_ratio;
        rr *= rhs_data.reverse_branching_ratio;

        frdt *= rhs_data.forward_branching_ratio;
        rrdt *= rhs_data.reverse_branching_ratio;

        // Store the final rate.

        r.forward(rate) = fr;
        drdt.forward(rate) = frdt;

        r.reverse(rate) = rr;
        drdt.reverse(rate) = rrdt;
    });
}

template<int spec>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real ener_gener_rate (Real const& dydt)
{
    using namespace network;

    // This is basically e = m c**2
    return dydt * mion(spec) * C::Legacy::enuc_conv2;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_rhs(burn_t& state, Array1D<Real, 1, neqs>& ydot)
{
    // Evaluate the rates

    rate_fr_t r, drdt;
    get_rates(state, r, drdt);

    // Call the RHS to get dydt

    Array1D<Real, 1, NumSpec> spec_rhs = RHS::species_rhs(state, r);

    for (int n = 1; n <= NumSpec; ++n) {
        ydot(n) = spec_rhs(n);
    }

    // Instantaneous energy generation rate

    ydot(net_ienuc) = RHS::ener_rhs(state, spec_rhs);
}


// Analytical Jacobian
template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_jac(burn_t& state, MatrixType& jac)
{
    rate_fr_t r, drdt;

    Real b1, sneut, dsneutdt, dsneutdd, snuda, snudz;
    Real rho, temp, abar, zbar;

    jac.zero();

    get_rates(state, r, drdt);

    // Get the data from the state

    rho  = state.rho;
    temp = state.T;
    abar = state.abar;
    zbar = state.zbar;

    // Species Jacobian elements with respect to other species

    RHS::constexpr_for<1, NumSpec+1>([&] (auto n1)
    {
        constexpr int spec1 = n1;
        RHS::constexpr_for<1, NumSpec+1>([&] (auto n2)
        {
            constexpr int spec2 = n2;
            RHS::constexpr_for<1, Rates::NumRatesFR+1>([&] (auto n3)
            {
                constexpr int rate = n3;
                if constexpr (RHS::is_rate_used<spec1, rate>()) {
                    jac(spec1, spec2) += RHS::jac_term<spec1, spec2, rate>(state, r);
                }
            });
        });
    });

    // Energy generation rate Jacobian elements with respect to species

    RHS::constexpr_for<1, NumSpec+1>([&] (auto j)
    {
        constexpr int spec2 = j;
        jac(net_ienuc, spec2) = 0.0;

        RHS::constexpr_for<1, NumSpec+1>([&] (auto i)
        {
            constexpr int spec1 = i;
            jac(net_ienuc, j) += ener_gener_rate<spec1>(jac(spec1, j));
        });
    });

    // Account for the thermal neutrino losses

    sneut5(temp, rho, abar, zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);

    for (int j = 1; j <= NumSpec; ++j) {
       b1 = (-abar * abar * snuda + (zion[j-1] - zbar) * abar * snudz);
       jac.add(net_ienuc, j, -b1);
    }

    // Evaluate the Jacobian elements with respect to energy by
    // calling the RHS using d(rate) / dT and then transform them
    // to our energy integration variable.

    Real jac_e_T = 0.0;

    RHS::constexpr_for<1, NumSpec+1>([&] (auto n1)
    {
        constexpr int species = n1;

        Real dydt = 0.0;

        RHS::constexpr_for<1, Rates::NumRatesFR+1>([&] (auto n2)
        {
            constexpr int rate = n2;
            if constexpr (RHS::is_rate_used<species, rate>()) {
                dydt += RHS::rhs_term<species, rate>(state, drdt);
            }
        });

        jac(species, net_ienuc) = temperature_to_energy_jacobian(state, dydt);

        jac_e_T += ener_gener_rate<species>(dydt);
    });

    jac_e_T -= dsneutdt;

    jac(net_ienuc, net_ienuc) = temperature_to_energy_jacobian(state, jac_e_T);
}

#endif
