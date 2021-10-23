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
void iso7tab(const Real btemp, const Real bden,
             rate_fr_t& r, rate_fr_t& drdt)
{
    using namespace RateTable;
    using namespace Rates;

    constexpr int mp = 4;

    int iat;
    Real x, x1, x2, x3, x4;
    Real a, b, c, d, e, f, g, h, p, q;
    Real alfa, beta, gama, delt;

    // hash locate
    iat = static_cast<int>((std::log10(btemp) - tab_tlo)/tab_tstp) + 1;
    iat = amrex::max(1, amrex::min(iat - 1, tab_imax - mp + 1));

    // setup the lagrange interpolation coefficients for a cubic
    x  = btemp;
    x1 = ttab(iat);
    x2 = ttab(iat+1);
    x3 = ttab(iat+2);
    x4 = ttab(iat+3);
    a  = x - x1;
    b  = x - x2;
    c  = x - x3;
    d  = x - x4;
    e  = x1 - x2;
    f  = x1 - x3;
    g  = x1 - x4;
    h  = x2 - x3;
    p  = x2 - x4;
    q  = x3 - x4;
    alfa =  b*c*d/(e*f*g);
    beta = -a*c*d/(e*h*p);
    gama =  a*b*d/(f*h*q);
    delt = -a*b*c/(g*p*q);

    RHS::constexpr_for<1, Rates::NumRatesFR+1>([&] (auto n)
    {
        constexpr int rate = n;

        // The rate array is flattened and for each reaction should contain the
        // forward rate and then the backward rate. So for any given rate index
        // we need to fill in 2 * index - 1 and 2 * index in the legacy flattened
        // scheme used by the rate array (whereas at every index the table has both
        // the forward and reverse terms).

        r.forward(rate) = (alfa * rattab(rate, 1, iat  ) +
                           beta * rattab(rate, 1, iat+1) +
                           gama * rattab(rate, 1, iat+2) +
                           delt * rattab(rate, 1, iat+3));

        r.reverse(rate) = (alfa * rattab(rate, 2, iat  ) +
                           beta * rattab(rate, 2, iat+1) +
                           gama * rattab(rate, 2, iat+2) +
                           delt * rattab(rate, 2, iat+3));

        drdt.forward(rate) = (alfa * drattabdt(rate, 1, iat  ) +
                              beta * drattabdt(rate, 1, iat+1) +
                              gama * drattabdt(rate, 1, iat+2) +
                              delt * drattabdt(rate, 1, iat+3));

        drdt.reverse(rate) = (alfa * drattabdt(rate, 2, iat  ) +
                              beta * drattabdt(rate, 2, iat+1) +
                              gama * drattabdt(rate, 2, iat+2) +
                              delt * drattabdt(rate, 2, iat+3));
    });
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void iso7rat(const Real btemp, const Real bden,
             rate_fr_t& r, rate_fr_t& drdt)
{
    using namespace Rates;

    // this routine generates unscreened
    // nuclear reaction rates for the iso7 network.

     for (int i = 1; i <= Rates::NumRatesFR; ++i) {
       r.forward(i) = 0.0_rt;
       r.reverse(i) = 0.0_rt;
       drdt.forward(i) = 0.0_rt;
       drdt.reverse(i) = 0.0_rt;
    }

    if (btemp < 1.0e6_rt) return;

    // get the temperature factors
    tf_t tf = get_tfactors(btemp);

    RHS::constexpr_for<1, NumRatesFR+1>([&] (auto n)
    {
        constexpr int rate = n;

        Real fr = 0.0_rt;
        Real frdt = 0.0_rt;

        Real rr = 0.0_rt;
        Real rrdt = 0.0_rt;

        RHS::evaluate_rate<rate>(tf, fr, frdt, rr, rrdt);

        r.forward(rate) = fr;
        drdt.forward(rate) = frdt,

        r.reverse(rate) = rr;
        drdt.reverse(rate) = rrdt;
    });
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void set_iso7rat()
{
    using namespace RateTable;

    Real btemp;
    Real bden = 1.0e0_rt;
    rate_fr_t rate;
    rate_fr_t dratedt;

    for (int i = 1; i <= tab_imax; ++i) {

       btemp = tab_tlo + static_cast<Real>(i-1) * tab_tstp;
       btemp = std::pow(10.0e0_rt, btemp);

       iso7rat(btemp, bden, rate, dratedt);

       ttab(i) = btemp;

       for (int j = 1; j <= Rates::NumRatesFR; ++j) {

           rattab(j, 1, i)    = rate.forward(j);
           rattab(j, 2, i)    = rate.reverse(j);
           drattabdt(j, 1, i) = dratedt.forward(j);
           drattabdt(j, 2, i) = dratedt.reverse(j);

       }
    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void screen_iso7(const Real btemp, const Real bden,
                 Array1D<Real, 1, NumSpec> const& y,
                 rate_fr_t& r, rate_fr_t& drdt)
{
    using namespace Species;
    using namespace Rates;

    /*
    this routine computes the screening factors
    and applies them to the raw reaction rates,
    producing the final reaction rates used by the
    right hand sides and jacobian matrix elements
    */

    int jscr;
    Real sc1a,sc1adt,sc1add,sc2a,sc2adt,sc2add;
    Real sc3a,sc3adt,sc3add,abar,zbar,ye,z2bar;
    Real t992,t9i92,yeff_ca40,yeff_ca40dt,yeff_ti44,yeff_ti44dt;
    Real denom,denomdt,denomdd,xx,zz;

    plasma_state_t pstate;

    // get the temperature factors
    tf_t tf = get_tfactors(btemp);

    // Set up the state data, which is the same for all screening factors.

    fill_plasma_state(pstate, btemp, bden, y);

    RHS::constexpr_for<1, NumRatesFR+1>([&] (auto n)
    {
        constexpr int rate = n;

        constexpr rhs_t rhs_data = RHS::rhs_data(rate);

        // The screening behavior depends on the type of reaction. We provide screening
        // here for the reaction classes we know about, and any other reactions are unscreened.

        if constexpr (rhs_data.number_A == 1 && rhs_data.number_B == 1 && rhs_data.number_C == 0) {
            // Forward reaction is A + B, screen using these two species

            constexpr Real Z1 = NetworkProperties::zion(rhs_data.species_A);
            constexpr Real A1 = NetworkProperties::aion(rhs_data.species_A);

            constexpr Real Z2 = NetworkProperties::zion(rhs_data.species_B);
            constexpr Real A2 = NetworkProperties::aion(rhs_data.species_B);

            constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

            Real sc, scdt;
            actual_screen5(pstate, scn_fac, sc, scdt);

            drdt.forward(rate) = drdt.forward(rate) * sc + r.forward(rate) * scdt;
            r.forward(rate)    = r.forward(rate) * sc;
        }

        if constexpr (rhs_data.number_A == 2 && rhs_data.number_B == 0 && rhs_data.number_C == 0) {
            // Forward reaction is A + A, screen using just this species

            constexpr Real Z1 = NetworkProperties::zion(rhs_data.species_A);
            constexpr Real A1 = NetworkProperties::aion(rhs_data.species_A);

            constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z1, A1);

            Real sc, scdt;
            actual_screen5(pstate, scn_fac, sc, scdt);

            drdt.forward(rate) = drdt.forward(rate) * sc + r.forward(rate) * scdt;
            r.forward(rate)    = r.forward(rate) * sc;
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

            drdt.forward(rate) = drdt.forward(rate) * sc + r.forward(rate) * scdt;
            r.forward(rate)    = r.forward(rate) * sc;
        }
    });

    /*
    the publication, timmes, woosley & hoffman apjs, 129, 377
    has a typo on page 393, where its says "y(C12)+y(O16) .gt. 0.004"
    it should be less than or equal to, since the idea is this piece
    gets activated during silicon buring, after all the c + o from
    oxygen burning is gone.
    */

    if (tf.t9 > 2.5_rt && y(C12) + y(O16) <= 4.0e-3_rt) {

       t992  = tf.t972 * tf.t9;
       t9i92 = 1.0e0_rt / t992;

       yeff_ca40   = t9i92 * std::exp(239.42_rt*tf.t9i - 74.741_rt);
       yeff_ca40dt = -yeff_ca40*(239.42_rt*tf.t9i2 + 4.5e0_rt*tf.t9i);

       yeff_ti44   = t992  * std::exp(-274.12_rt*tf.t9i + 74.914_rt);
       yeff_ti44dt = yeff_ti44*(274.12_rt*tf.t9i2 + 4.5e0_rt*tf.t9i);

       denom     = std::pow(bden * y(He4), 3.0e0_rt);

       r.forward(Si28_7He4_to_Ni56)     = yeff_ca40 * denom * r.forward(Ca40_He4_to_Ti44);
       drdt.forward(Si28_7He4_to_Ni56)  = (yeff_ca40dt * r.forward(Ca40_He4_to_Ti44) +
                                              yeff_ca40 * drdt.forward(Ca40_He4_to_Ti44)) * denom * 1.0e-9_rt;

       if (denom != 0.0_rt) {

          zz     = 1.0e0_rt / denom;
          r.reverse(Si28_7He4_to_Ni56) = amrex::min(1.0e10_rt, yeff_ti44 * r.reverse(Ca40_He4_to_Ti44) * zz);

          if (r.reverse(Si28_7He4_to_Ni56) == 1.0e10_rt) {
             drdt.reverse(Si28_7He4_to_Ni56)  = 0.0e0_rt;
          } else {
             drdt.reverse(Si28_7He4_to_Ni56)  = (yeff_ti44dt * r.reverse(Ca40_He4_to_Ti44) +
                                                    yeff_ti44 * drdt.reverse(Ca40_He4_to_Ti44)) * zz * 1.0e-9_rt;
          }

       }

    }

    RHS::constexpr_for<1, NumRatesFR+1>([&] (auto n)
    {
        constexpr int rate = n;

        constexpr rhs_t rhs_data = RHS::rhs_data(rate);

        r.forward(rate) *= rhs_data.forward_branching_ratio;
        r.reverse(rate) *= rhs_data.reverse_branching_ratio;

        drdt.forward(rate) *= rhs_data.forward_branching_ratio;
        drdt.reverse(rate) *= rhs_data.reverse_branching_ratio;
    });
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_rates(burn_t const& state, rate_fr_t& r, rate_fr_t& drdt)
{
    Real rho, temp;
    Array1D<Real, 1, NumSpec> y;

    // Get the data from the state
    rho  = state.rho;
    temp = state.T;

    for (int i = 1; i <= NumSpec; ++i)
        y(i) = state.xn[i-1] * aion_inv[i-1];

    // Get the raw reaction rates
    if (use_tables)
        iso7tab(temp, rho, r, drdt);
    else
        iso7rat(temp, rho, r, drdt);

    RHS::constexpr_for<1, Rates::NumRatesFR+1>([&] (auto n)
    {
        constexpr int rate = n;

        // Set the density dependence

        int forward_exponent = RHS::density_exponent_forward<rate>();
        int reverse_exponent = RHS::density_exponent_reverse<rate>();

        Real density_term_forward = std::pow(rho, forward_exponent);
        Real density_term_reverse = std::pow(rho, reverse_exponent);

        r.forward(rate) *= density_term_forward;
        drdt.forward(rate) *= density_term_forward;
        r.reverse(rate) *= density_term_reverse;
        drdt.reverse(rate) *= density_term_reverse;
    });

    // Do the screening here because the corrections depend on the composition
    screen_iso7(temp, rho, y, r, drdt);
}

template<class T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ener_gener_rate(T const& dydt, Real& enuc)
{

    using namespace network;

    // Computes the instantaneous energy generation rate

    Real Xdot = 0.0_rt;

    // Sum the mass fraction time derivatives
    for (int i = 1; i <= NumSpec; ++i) {
        Xdot += dydt(i) * mion(i);
    }

    // This is basically e = m c**2
    enuc = Xdot * C::Legacy::enuc_conv2;

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

    rate_fr_t rate, dratedt;

    bool deriva;

    Real b1, sneut, dsneutdt, dsneutdd, snuda, snudz;

    Real rho, temp, abar, zbar;
    Array1D<Real, 1, NumSpec> y;

    jac.zero();

    get_rates(state, rate, dratedt);

    // Get the data from the state

    rho  = state.rho;
    temp = state.T;
    abar = state.abar;
    zbar = state.zbar;

    for (int i = 1; i <= NumSpec; ++i)
        y(i) = state.xn[i-1] * aion_inv[i-1];

    // Species Jacobian elements with respect to other species

    auto spec_jac = RHS::species_jac(state, rate);

    for (int spec1 = 1; spec1 <= NumSpec; ++spec1) {
        for (int spec2 = 1; spec2 <= NumSpec; ++spec2) {
            jac(spec1, spec2) = spec_jac(spec1, spec2);
        }
    }

    // Energy generation rate Jacobian elements with respect to species

    for (int j = 1; j <= NumSpec; ++j) {
        auto jac_slice_2 = [&](int i) -> Real { return jac.get(i, j); };
        ener_gener_rate(jac_slice_2, jac(net_ienuc,j));
    }

    // Account for the thermal neutrino losses

    sneut5(temp, rho, abar, zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);

    for (int j = 1; j <= NumSpec; ++j) {
       b1 = (-abar * abar * snuda + (zion[j-1] - zbar) * abar * snudz);
       jac.add(net_ienuc, j, -b1);
    }

    // Evaluate the Jacobian elements with respect to energy by
    // calling the RHS using d(rate) / dT and then transform them
    // to our energy integration variable.

    Array1D<Real, 1, NumSpec> yderivs = RHS::species_rhs(state, dratedt);

    for (int i = 1; i <= NumSpec; ++i) {
        jac(i, net_ienuc) = temperature_to_energy_jacobian(state, yderivs(i));
    }

    Real jac_e_T;
    ener_gener_rate(yderivs, jac_e_T);
    jac_e_T -= dsneutdt;
    jac(net_ienuc, net_ienuc) = temperature_to_energy_jacobian(state, jac_e_T);

}


AMREX_INLINE
void set_up_screening_factors()
{
    // Compute and store the more expensive screening factors

    using namespace Species;

    // note: we need to set these up in the same order that we evaluate the
    // rates in actual_rhs.H (yes, it's ugly)
    int jscr = 0;
    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[He4-1], aion[He4-1]);
    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Be8-1], aion[Be8-1]);
    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[He4-1], aion[He4-1]);
    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[C12-1], aion[C12-1]);
    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[O16-1], aion[O16-1]);
    add_screening_factor(jscr++, zion[O16-1], aion[O16-1], zion[O16-1], aion[O16-1]);
    add_screening_factor(jscr++, zion[O16-1], aion[O16-1], zion[He4-1], aion[He4-1]);
    add_screening_factor(jscr++, zion[Ne20-1], aion[Ne20-1], zion[He4-1], aion[He4-1]);
    add_screening_factor(jscr++, zion[Mg24-1], aion[Mg24-1], zion[He4-1], aion[He4-1]);
    add_screening_factor(jscr++, zion[Ca40-1], aion[Ca40-1], zion[He4-1], aion[He4-1]);

}

#endif
