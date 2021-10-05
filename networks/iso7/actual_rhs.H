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
             rate_fr_t& rate, rate_fr_t& dratedt)
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

    // crank off the raw reaction rates
    for (int rateindex = 1; rateindex <= Rates::NumRatesFR; ++rateindex) {

        Real forward_dtab = 0.0_rt;
        Real reverse_dtab = 0.0_rt;

        int forward_exponent, reverse_exponent;

        // Set the density dependence

        switch (rateindex) {

        case (C12_He4_to_O16):
            forward_exponent = RHS::density_exponent_forward<C12_He4_to_O16>();
            reverse_exponent = RHS::density_exponent_reverse<C12_He4_to_O16>();

            forward_dtab = std::pow(bden, forward_exponent);
            reverse_dtab = std::pow(bden, reverse_exponent);
            break;

        case (He4_He4_He4_to_C12):
            forward_exponent = RHS::density_exponent_forward<He4_He4_He4_to_C12>();
            reverse_exponent = RHS::density_exponent_reverse<He4_He4_He4_to_C12>();

            forward_dtab = std::pow(bden, forward_exponent);
            reverse_dtab = std::pow(bden, reverse_exponent);
            break;

        case (C12_C12_to_Ne20_He4):
            forward_exponent = RHS::density_exponent_forward<C12_C12_to_Ne20_He4>();
            reverse_exponent = RHS::density_exponent_reverse<C12_C12_to_Ne20_He4>();

            forward_dtab = std::pow(bden, forward_exponent);
            reverse_dtab = std::pow(bden, reverse_exponent); // rate is zero in this net
            break;

        case (C12_O16_to_Mg24_He4):
            forward_exponent = RHS::density_exponent_forward<C12_O16_to_Mg24_He4>();
            reverse_exponent = RHS::density_exponent_reverse<C12_O16_to_Mg24_He4>();

            forward_dtab = std::pow(bden, forward_exponent);
            reverse_dtab = std::pow(bden, reverse_exponent); // rate is zero in this net
            break;

        case (C12_O16_to_Si28):
            forward_exponent = RHS::density_exponent_forward<C12_O16_to_Si28>();
            reverse_exponent = RHS::density_exponent_reverse<C12_O16_to_Si28>();

            forward_dtab = std::pow(bden, forward_exponent);
            reverse_dtab = std::pow(bden, reverse_exponent); // rate is zero in this net
            break;

        case (O16_O16_to_Si28_He4):
            forward_exponent = RHS::density_exponent_forward<O16_O16_to_Si28_He4>();
            reverse_exponent = RHS::density_exponent_reverse<O16_O16_to_Si28_He4>();

            forward_dtab = std::pow(bden, forward_exponent);
            reverse_dtab = std::pow(bden, reverse_exponent); // rate is zero in this net
            break;

        case (O16_He4_to_Ne20):
            forward_exponent = RHS::density_exponent_forward<O16_He4_to_Ne20>();
            reverse_exponent = RHS::density_exponent_reverse<O16_He4_to_Ne20>();

            forward_dtab = std::pow(bden, forward_exponent);
            reverse_dtab = std::pow(bden, reverse_exponent);
            break;

        case (Ne20_He4_to_Mg24):
            forward_exponent = RHS::density_exponent_forward<Ne20_He4_to_Mg24>();
            reverse_exponent = RHS::density_exponent_reverse<Ne20_He4_to_Mg24>();

            forward_dtab = std::pow(bden, forward_exponent);
            reverse_dtab = std::pow(bden, reverse_exponent);
            break;

        case (Mg24_He4_to_Si28):
            forward_exponent = RHS::density_exponent_forward<Mg24_He4_to_Si28>();
            reverse_exponent = RHS::density_exponent_reverse<Mg24_He4_to_Si28>();

            forward_dtab = std::pow(bden, forward_exponent);
            reverse_dtab = std::pow(bden, reverse_exponent);
            break;

        case (Ca40_He4_to_Ti44):
            forward_dtab = bden;
            reverse_dtab = 1.0e0_rt;
            break;

        case (Si28_7He4_to_Ni56):
            // This reaction is evaluated in the screening routine
            forward_dtab = 0.0e0_rt;
            reverse_dtab = 0.0e0_rt;
            break;

        }

        // The rate array is flattened and for each reaction should contain the
        // forward rate and then the backward rate. So for any given rate index
        // we need to fill in 2 * index - 1 and 2 * index in the legacy flattened
        // scheme used by the rate array (whereas at every index the table has both
        // the forward and reverse terms).

        rate.forward(rateindex) = (alfa * rattab(rateindex, 1, iat  ) +
                                   beta * rattab(rateindex, 1, iat+1) +
                                   gama * rattab(rateindex, 1, iat+2) +
                                   delt * rattab(rateindex, 1, iat+3)) * forward_dtab;

        rate.reverse(rateindex) = (alfa * rattab(rateindex, 2, iat  ) +
                                   beta * rattab(rateindex, 2, iat+1) +
                                   gama * rattab(rateindex, 2, iat+2) +
                                   delt * rattab(rateindex, 2, iat+3)) * reverse_dtab;

        dratedt.forward(rateindex) = (alfa * drattabdt(rateindex, 1, iat  ) +
                                      beta * drattabdt(rateindex, 1, iat+1) +
                                      gama * drattabdt(rateindex, 1, iat+2) +
                                      delt * drattabdt(rateindex, 1, iat+3)) * forward_dtab;

        dratedt.reverse(rateindex) = (alfa * drattabdt(rateindex, 2, iat  ) +
                                      beta * drattabdt(rateindex, 2, iat+1) +
                                      gama * drattabdt(rateindex, 2, iat+2) +
                                      delt * drattabdt(rateindex, 2, iat+3)) * reverse_dtab;

    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void iso7rat(const Real btemp, const Real bden,
             rate_fr_t& rate, rate_fr_t& dratedt)
{
    using namespace Rates;

    // this routine generates unscreened
    // nuclear reaction rates for the iso7 network.

     for (int i = 1; i <= Rates::NumRatesFR; ++i) {
       rate.forward(i) = 0.0_rt;
       rate.reverse(i) = 0.0_rt;
       dratedt.forward(i) = 0.0_rt;
       dratedt.reverse(i) = 0.0_rt;
    }

    if (btemp < 1.0e6_rt) return;

    // get the temperature factors
    tf_t tf = get_tfactors(btemp);

    Real density_term_forward, density_term_reverse;

    // c12(a,g)o16
    rate_c12ag(tf, 1.0_rt,
               rate.forward(C12_He4_to_O16), dratedt.forward(C12_He4_to_O16),
               rate.reverse(C12_He4_to_O16), dratedt.reverse(C12_He4_to_O16));

    density_term_forward = std::pow(bden, RHS::density_exponent_forward<C12_He4_to_O16>());
    density_term_reverse = std::pow(bden, RHS::density_exponent_reverse<C12_He4_to_O16>());

    rate.forward(C12_He4_to_O16) *= density_term_forward;
    dratedt.forward(C12_He4_to_O16) *= density_term_forward;
    rate.reverse(C12_He4_to_O16) *= density_term_reverse;
    dratedt.reverse(C12_He4_to_O16) *= density_term_reverse;

    // triple alpha to c12
    rate_triplealf(tf, 1.0_rt,
                   rate.forward(He4_He4_He4_to_C12), dratedt.forward(He4_He4_He4_to_C12),
                   rate.reverse(He4_He4_He4_to_C12), dratedt.reverse(He4_He4_He4_to_C12));

    density_term_forward = std::pow(bden, RHS::density_exponent_forward<He4_He4_He4_to_C12>());
    density_term_reverse = std::pow(bden, RHS::density_exponent_reverse<He4_He4_He4_to_C12>());

    rate.forward(He4_He4_He4_to_C12) *= density_term_forward;
    dratedt.forward(He4_He4_He4_to_C12) *= density_term_forward;
    rate.reverse(He4_He4_He4_to_C12) *= density_term_reverse;
    dratedt.reverse(He4_He4_He4_to_C12) *= density_term_reverse;

    // c12 + c12
    rate_c12c12(tf, 1.0_rt,
                rate.forward(C12_C12_to_Ne20_He4), dratedt.forward(C12_C12_to_Ne20_He4),
                rate.reverse(C12_C12_to_Ne20_He4), dratedt.reverse(C12_C12_to_Ne20_He4));

    density_term_forward = std::pow(bden, RHS::density_exponent_forward<C12_C12_to_Ne20_He4>());
    density_term_reverse = std::pow(bden, RHS::density_exponent_reverse<C12_C12_to_Ne20_He4>());

    rate.forward(C12_C12_to_Ne20_He4) *= density_term_forward;
    dratedt.forward(C12_C12_to_Ne20_He4) *= density_term_forward;
    rate.reverse(C12_C12_to_Ne20_He4) *= density_term_reverse;
    dratedt.reverse(C12_C12_to_Ne20_He4) *= density_term_reverse;

    // c12 + o16
    rate_c12o16(tf, 1.0_rt,
                rate.forward(C12_O16_to_Mg24_He4), dratedt.forward(C12_O16_to_Mg24_He4),
                rate.reverse(C12_O16_to_Mg24_He4), dratedt.reverse(C12_O16_to_Mg24_He4));

    density_term_forward = std::pow(bden, RHS::density_exponent_forward<C12_O16_to_Mg24_He4>());
    density_term_reverse = std::pow(bden, RHS::density_exponent_reverse<C12_O16_to_Mg24_He4>());

    rate.forward(C12_O16_to_Mg24_He4) *= density_term_forward;
    dratedt.forward(C12_O16_to_Mg24_He4) *= density_term_forward;
    rate.reverse(C12_O16_to_Mg24_He4) *= density_term_reverse;
    dratedt.reverse(C12_O16_to_Mg24_He4) *= density_term_reverse;

    // 16o + 16o
    rate_o16o16(tf, 1.0_rt,
                rate.forward(O16_O16_to_Si28_He4), dratedt.forward(O16_O16_to_Si28_He4),
                rate.reverse(O16_O16_to_Si28_He4), dratedt.reverse(O16_O16_to_Si28_He4));

    density_term_forward = std::pow(bden, RHS::density_exponent_forward<O16_O16_to_Si28_He4>());
    density_term_reverse = std::pow(bden, RHS::density_exponent_reverse<O16_O16_to_Si28_He4>());

    rate.forward(O16_O16_to_Si28_He4) *= density_term_forward;
    dratedt.forward(O16_O16_to_Si28_He4) *= density_term_forward;
    rate.reverse(O16_O16_to_Si28_He4) *= density_term_reverse;
    dratedt.reverse(O16_O16_to_Si28_He4) *= density_term_reverse;

    // o16(a,g)ne20
    rate_o16ag(tf, 1.0_rt,
               rate.forward(O16_He4_to_Ne20), dratedt.forward(O16_He4_to_Ne20),
               rate.reverse(O16_He4_to_Ne20), dratedt.reverse(O16_He4_to_Ne20));

    density_term_forward = std::pow(bden, RHS::density_exponent_forward<O16_He4_to_Ne20>());
    density_term_reverse = std::pow(bden, RHS::density_exponent_reverse<O16_He4_to_Ne20>());

    rate.forward(O16_He4_to_Ne20) *= density_term_forward;
    dratedt.forward(O16_He4_to_Ne20) *= density_term_forward;
    rate.reverse(O16_He4_to_Ne20) *= density_term_reverse;
    dratedt.reverse(O16_He4_to_Ne20) *= density_term_reverse;

    // ne20(a,g)mg24
    rate_ne20ag(tf, 1.0_rt,
                rate.forward(Ne20_He4_to_Mg24), dratedt.forward(Ne20_He4_to_Mg24),
                rate.reverse(Ne20_He4_to_Mg24), dratedt.reverse(Ne20_He4_to_Mg24));

    density_term_forward = std::pow(bden, RHS::density_exponent_forward<Ne20_He4_to_Mg24>());
    density_term_reverse = std::pow(bden, RHS::density_exponent_reverse<Ne20_He4_to_Mg24>());

    rate.forward(Ne20_He4_to_Mg24) *= density_term_forward;
    dratedt.forward(Ne20_He4_to_Mg24) *= density_term_forward;
    rate.reverse(Ne20_He4_to_Mg24) *= density_term_reverse;
    dratedt.reverse(Ne20_He4_to_Mg24) *= density_term_reverse;

    // mg24(a,g)si28
    rate_mg24ag(tf, 1.0_rt,
                rate.forward(Mg24_He4_to_Si28), dratedt.forward(Mg24_He4_to_Si28),
                rate.reverse(Mg24_He4_to_Si28), dratedt.reverse(Mg24_He4_to_Si28));

    density_term_forward = std::pow(bden, RHS::density_exponent_forward<Mg24_He4_to_Si28>());
    density_term_reverse = std::pow(bden, RHS::density_exponent_reverse<Mg24_He4_to_Si28>());

    rate.forward(Mg24_He4_to_Si28) *= density_term_forward;
    dratedt.forward(Mg24_He4_to_Si28) *= density_term_forward;
    rate.reverse(Mg24_He4_to_Si28) *= density_term_reverse;
    dratedt.reverse(Mg24_He4_to_Si28) *= density_term_reverse;

    // ca40(a,g)ti44
    rate_ca40ag(tf, bden,
                rate.forward(Ca40_He4_to_Ti44), dratedt.forward(Ca40_He4_to_Ti44),
                rate.reverse(Ca40_He4_to_Ti44), dratedt.reverse(Ca40_He4_to_Ti44));

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
                 rate_fr_t& rate, rate_fr_t& dratedt,
                 Array2D<Real, 1, NumSpec, 1, Rates::NumRates>& dratedY)
{
    using namespace Species;
    using namespace ExtraSpecies;
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

    // initialize
    for (int j = 1; j <= NumSpec; ++j) {
        for (int i = 1; i <= Rates::NumRates; ++i) {
            dratedY(j, i) = 0.0_rt;
        }
    }

    // get the temperature factors
    tf_t tf = get_tfactors(btemp);

    // Set up the state data, which is the same for all screening factors.

    fill_plasma_state(pstate, btemp, bden, y);

    // first the always fun triple alpha and its inverse
    jscr = 0;
    screen5(pstate,jscr,
            zion[He4-1], aion[He4-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);

    jscr++;
    screen5(pstate,jscr,
            zion[He4-1], aion[He4-1], extra_zion[Be8-1], extra_aion[Be8-1],
            sc2a,sc2adt,sc2add);

    sc3a   = sc1a * sc2a;
    sc3adt = sc1adt*sc2a + sc1a*sc2adt;

    dratedt.forward(He4_He4_He4_to_C12) = dratedt.forward(He4_He4_He4_to_C12) * sc3a + rate.forward(He4_He4_He4_to_C12) * sc3adt;
    rate.forward(He4_He4_He4_to_C12)    = rate.forward(He4_He4_He4_to_C12) * sc3a;

    // c12 to o16
    jscr++;
    screen5(pstate,jscr,
            zion[C12-1], aion[C12-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);

    dratedt.forward(C12_He4_to_O16)  = dratedt.forward(C12_He4_to_O16) * sc1a + rate.forward(C12_He4_to_O16) * sc1adt;
    rate.forward(C12_He4_to_O16)     = rate.forward(C12_He4_to_O16) * sc1a;

    // c12 + c12
    jscr++;
    screen5(pstate,jscr,
            zion[C12-1], aion[C12-1], zion[C12-1], aion[C12-1],
            sc1a,sc1adt,sc1add);

    dratedt.forward(C12_C12_to_Ne20_He4) = dratedt.forward(C12_C12_to_Ne20_He4) * sc1a + rate.forward(C12_C12_to_Ne20_He4) * sc1adt;
    rate.forward(C12_C12_to_Ne20_He4)    = rate.forward(C12_C12_to_Ne20_He4) * sc1a;

    // c12 + o16
    jscr++;
    screen5(pstate,jscr,
            zion[C12-1], aion[C12-1], zion[O16-1], aion[O16-1],
            sc1a,sc1adt,sc1add);

    dratedt.forward(C12_O16_to_Mg24_He4) = dratedt.forward(C12_O16_to_Mg24_He4) * sc1a + rate.forward(C12_O16_to_Mg24_He4) * sc1adt;
    rate.forward(C12_O16_to_Mg24_He4)    = rate.forward(C12_O16_to_Mg24_He4) * sc1a;

    // Equal probability branching ratio between (Mg24 + He4) and (Si28) endpoints

    rate.forward(C12_O16_to_Mg24_He4) *= 0.5;
    dratedt.forward(C12_O16_to_Mg24_He4) *= 0.5;

    rate.forward(C12_O16_to_Si28) = rate.forward(C12_O16_to_Mg24_He4);
    dratedt.forward(C12_O16_to_Si28) = dratedt.forward(C12_O16_to_Mg24_He4);

    // o16 + o16
    jscr++;
    screen5(pstate,jscr,
            zion[O16-1], aion[O16-1], zion[O16-1], aion[O16-1],
            sc1a,sc1adt,sc1add);

    dratedt.forward(O16_O16_to_Si28_He4) = dratedt.forward(O16_O16_to_Si28_He4) * sc1a + rate.forward(O16_O16_to_Si28_He4) * sc1adt;
    rate.forward(O16_O16_to_Si28_He4)    = rate.forward(O16_O16_to_Si28_He4) * sc1a;

    // o16 to ne20
    jscr++;
    screen5(pstate,jscr,
            zion[O16-1], aion[O16-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);

    dratedt.forward(O16_He4_to_Ne20) = dratedt.forward(O16_He4_to_Ne20) * sc1a + rate.forward(O16_He4_to_Ne20) * sc1adt;
    rate.forward(O16_He4_to_Ne20)    = rate.forward(O16_He4_to_Ne20) * sc1a;

    // ne20 to mg24
    jscr++;
    screen5(pstate,jscr,
            zion[Ne20-1], aion[Ne20-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);

    dratedt.forward(Ne20_He4_to_Mg24) = dratedt.forward(Ne20_He4_to_Mg24) * sc1a + rate.forward(Ne20_He4_to_Mg24) * sc1adt;
    rate.forward(Ne20_He4_to_Mg24)    = rate.forward(Ne20_He4_to_Mg24) * sc1a;

    // mg24 to si28
    jscr++;
    screen5(pstate,jscr,
            zion[Mg24-1], aion[Mg24-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);

    dratedt.forward(Mg24_He4_to_Si28) = dratedt.forward(Mg24_He4_to_Si28) * sc1a + rate.forward(Mg24_He4_to_Si28) * sc1adt;
    rate.forward(Mg24_He4_to_Si28)    = rate.forward(Mg24_He4_to_Si28) * sc1a;

    // ca40 to ti44
    jscr++;
    screen5(pstate,jscr,
            extra_zion[Ca40-1], extra_aion[Ca40-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);

    dratedt.forward(Ca40_He4_to_Ti44) = dratedt.forward(Ca40_He4_to_Ti44) * sc1a + rate.forward(Ca40_He4_to_Ti44) * sc1adt;
    rate.forward(Ca40_He4_to_Ti44)    = rate.forward(Ca40_He4_to_Ti44) * sc1a;

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

       rate.forward(Si28_7He4_to_Ni56)     = yeff_ca40 * denom * rate.forward(Ca40_He4_to_Ti44);
       dratedt.forward(Si28_7He4_to_Ni56)  = (yeff_ca40dt * rate.forward(Ca40_He4_to_Ti44) +
                                              yeff_ca40 * dratedt.forward(Ca40_He4_to_Ti44)) * denom * 1.0e-9_rt;

       if (denom != 0.0_rt) {

          zz     = 1.0e0_rt / denom;
          rate.reverse(Si28_7He4_to_Ni56) = amrex::min(1.0e10_rt, yeff_ti44 * rate.reverse(Ca40_He4_to_Ti44) * zz);

          if (rate.reverse(Si28_7He4_to_Ni56) == 1.0e10_rt) {
             dratedt.reverse(Si28_7He4_to_Ni56)  = 0.0e0_rt;
          } else {
             dratedt.reverse(Si28_7He4_to_Ni56)  = (yeff_ti44dt * rate.reverse(Ca40_He4_to_Ti44) +
                                                    yeff_ti44 * dratedt.reverse(Ca40_He4_to_Ti44)) * zz * 1.0e-9_rt;
          }

       }

    }
}



// RHS version (single rate group)
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void evaluate_rates(burn_t const& state, rate_fr_t& rate)
{
    Real rho, temp;
    Array1D<Real, 1, NumSpec> y;

    // Temporary arrays for storing unused derivatives.
    // A useful cleanup would be to rewrite iso7tab, iso7rat,
    // and screen_iso7 so that there are versions that do not
    // evaluate the derivatives.
    rate_fr_t dratedt;
    Array2D<Real, 1, NumSpec, 1, Rates::NumRates> dratedY;

    // Get the data from the state
    rho  = state.rho;
    temp = state.T;

    for (int i = 1; i <= NumSpec; ++i) {
        y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    // Get the raw reaction rates
    if (use_tables)
        iso7tab(temp, rho, rate, dratedt);
    else
        iso7rat(temp, rho, rate, dratedt);

    // Do the screening (depends on composition)
    screen_iso7(temp, rho, y, rate, dratedt, dratedY);
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_rates(burn_t const& state, rate_fr_t& rate, rate_fr_t& dratedt)
{
    Real rho, temp;
    Array1D<Real, 1, NumSpec> y;
    Array2D<Real, 1, NumSpec, 1, Rates::NumRates> dratedY;

    // Get the data from the state
    rho  = state.rho;
    temp = state.T;

    for (int i = 1; i <= NumSpec; ++i)
        y(i) = state.xn[i-1] * aion_inv[i-1];

    // Get the raw reaction rates
    if (use_tables)
        iso7tab(temp, rho, rate, dratedt);
    else
        iso7rat(temp, rho, rate, dratedt);

    // Do the screening here because the corrections depend on the composition
    screen_iso7(temp, rho, y, rate, dratedt, dratedY);
}


template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void dfdy_isotopes_iso7(Array1D<Real, 1, NumSpec> const& y,
                        burn_t const& state, rate_fr_t rate,
                        MatrixType& jac)
{

    // this routine sets up the dense iso7 jacobian for the isotopes

    using namespace Species;
    using namespace Rates;

    auto spec_jac = RHS::species_jac(state, rate);

    for (int spec1 = 1; spec1 <= NumSpec; ++spec1) {
        for (int spec2 = 1; spec2 <= NumSpec; ++spec2) {
            jac(spec1, spec2) = spec_jac(spec1, spec2);
        }
    }

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

    rate_fr_t rr;
    evaluate_rates(state, rr);

    // Call the RHS to get dydt

    Array1D<Real, 1, NumSpec> spec_rhs = RHS::species_rhs(state, rr);

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

    dfdy_isotopes_iso7(y, state, rate, jac);

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
    using namespace ExtraSpecies;

    // note: we need to set these up in the same order that we evaluate the
    // rates in actual_rhs.H (yes, it's ugly)
    int jscr = 0;
    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[He4-1], aion[He4-1]);
    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], extra_zion[Be8-1], extra_aion[Be8-1]);
    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[He4-1], aion[He4-1]);
    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[C12-1], aion[C12-1]);
    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[O16-1], aion[O16-1]);
    add_screening_factor(jscr++, zion[O16-1], aion[O16-1], zion[O16-1], aion[O16-1]);
    add_screening_factor(jscr++, zion[O16-1], aion[O16-1], zion[He4-1], aion[He4-1]);
    add_screening_factor(jscr++, zion[Ne20-1], aion[Ne20-1], zion[He4-1], aion[He4-1]);
    add_screening_factor(jscr++, zion[Mg24-1], aion[Mg24-1], zion[He4-1], aion[He4-1]);
    add_screening_factor(jscr++, extra_zion[Ca40-1], extra_aion[Ca40-1], zion[He4-1], aion[He4-1]);

}

#endif
