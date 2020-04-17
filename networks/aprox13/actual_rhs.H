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
#include <temperature_integration.H>

using namespace amrex;
using namespace ArrayUtil;

void actual_rhs_init();

namespace RateTable
{
    constexpr Real tab_tlo = 6.0e0_rt;
    constexpr Real tab_thi = 10.0e0_rt;
    constexpr int tab_per_decade = 500;
    constexpr int nrattab = static_cast<int>(tab_thi - tab_tlo) * tab_per_decade + 1;
    constexpr int tab_imax = static_cast<int>(tab_thi - tab_tlo) * tab_per_decade + 1;
    constexpr Real tab_tstp = (tab_thi - tab_tlo) / static_cast<Real>(tab_imax - 1);

    extern AMREX_GPU_MANAGED Array2D<Real, 1, Rates::NumRates, 1, nrattab> rattab;
    extern AMREX_GPU_MANAGED Array2D<Real, 1, Rates::NumRates, 1, nrattab> drattabdt;
    extern AMREX_GPU_MANAGED Array1D<Real, 1, nrattab> ttab;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void aprox13tab(const Real btemp, const Real bden, rate_t& rr)
{
    using namespace RateTable;

    constexpr int mp = 4;

    int iat;
    Real x, x1, x2, x3, x4;
    Real a, b, c, d, e, f, g, h, p, q;
    Real alfa, beta, gama, delt;
    Array1D<Real, 1, Rates::NumRates> dtab;

    // Set the density dependence array
    {
        using namespace Rates;
        dtab(ir3a)   = bden*bden;
        dtab(irg3a)  = 1.0e0_rt;
        dtab(ircag)  = bden;
        dtab(iroga)  = 1.0e0_rt;
        dtab(ir1212) = bden;
        dtab(ir1216) = bden;
        dtab(ir1616) = bden;
        dtab(iroag)  = bden;
        dtab(irnega) = 1.0e0_rt;
        dtab(irneag) = bden;
        dtab(irmgga) = 1.0e0_rt;
        dtab(irmgag) = bden;
        dtab(irsiga) = 1.0e0_rt;
        dtab(irmgap) = bden;
        dtab(iralpa) = bden;
        dtab(iralpg) = bden;
        dtab(irsigp) = 1.0e0_rt;
        dtab(irsiag) = bden;
        dtab(irsga)  = 1.0e0_rt;
        dtab(irppa)  = bden;
        dtab(irsiap) = bden;
        dtab(irppg)  = bden;
        dtab(irsgp)  = 1.0e0_rt;
        dtab(irsag)  = bden;
        dtab(irarga) = 1.0e0_rt;
        dtab(irsap)  = bden;
        dtab(irclpa) = bden;
        dtab(irclpg) = bden;
        dtab(irargp) = 1.0e0_rt;
        dtab(irarag) = bden;
        dtab(ircaga) = 1.0e0_rt;
        dtab(irarap) = bden;
        dtab(irkpa)  = bden;
        dtab(irkpg)  = bden;
        dtab(ircagp) = 1.0e0_rt;
        dtab(ircaag) = bden;
        dtab(irtiga) = 1.0e0_rt;
        dtab(ircaap) = bden;
        dtab(irscpa) = bden;
        dtab(irscpg) = bden;
        dtab(irtigp) = 1.0e0_rt;
        dtab(irtiag) = bden;
        dtab(ircrga) = 1.0e0_rt;
        dtab(irtiap) = bden;
        dtab(irvpa)  = bden;
        dtab(irvpg)  = bden;
        dtab(ircrgp) = 1.0e0_rt;
        dtab(ircrag) = bden;
        dtab(irfega) = 1.0e0_rt;
        dtab(ircrap) = bden;
        dtab(irmnpa) = bden;
        dtab(irmnpg) = bden;
        dtab(irfegp) = 1.0e0_rt;
        dtab(irfeag) = bden;
        dtab(irniga) = 1.0e0_rt;
        dtab(irfeap) = bden;
        dtab(ircopa) = bden;
        dtab(ircopg) = bden;
        dtab(irnigp) = 1.0e0_rt;
        dtab(irr1)   = 0.0e0_rt;
        dtab(irs1)   = 0.0e0_rt;
        dtab(irt1)   = 0.0e0_rt;
        dtab(iru1)   = 0.0e0_rt;
        dtab(irv1)   = 0.0e0_rt;
        dtab(irw1)   = 0.0e0_rt;
        dtab(irx1)   = 0.0e0_rt;
        dtab(iry1)   = 0.0e0_rt;
    }

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
    for (int j = 1; j <= Rates::NumRates; ++j) {

       rr.rates(1,j) = (  alfa * rattab(j,iat  )
                        + beta * rattab(j,iat+1)
                        + gama * rattab(j,iat+2)
                        + delt * rattab(j,iat+3) ) * dtab(j);

       rr.rates(2,j) = (  alfa * drattabdt(j,iat  )
                        + beta * drattabdt(j,iat+1)
                        + gama * drattabdt(j,iat+2)
                        + delt * drattabdt(j,iat+3) ) * dtab(j);

    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void aprox13rat(const Real btemp, const Real bden, rate_t& rr)
{
    using namespace Rates;

    // this routine generates unscreened
    // nuclear reaction rates for the aprox13 network.

    Real rrate,drratedt;

    for (int i = 1; i <= Rates::NumRates; ++i) {
       rr.rates(1,i) = 0.0_rt;
       rr.rates(2,i) = 0.0_rt;
    }

    if (btemp < 1.0e6_rt) return;


    // get the temperature factors
    tf_t tf = get_tfactors(btemp);


    // Determine which c12(a,g)o16 rate to use
    if (use_c12ag_deboer17) {
        // deboer + 2017 c12(a,g)o16 rate
        rate_c12ag_deboer17(tf,bden,
                            rr.rates(1,ircag),rr.rates(2,ircag),
                            rr.rates(1,iroga),rr.rates(2,iroga));
    } else {
        // 1.7 times cf88 c12(a,g)o16 rate
        rate_c12ag(tf,bden,
                   rr.rates(1,ircag),rr.rates(2,ircag),
                   rr.rates(1,iroga),rr.rates(2,iroga));
    }

    // triple alpha to c12
    rate_triplealf(tf,bden,
                   rr.rates(1,ir3a),rr.rates(2,ir3a),
                   rr.rates(1,irg3a),rr.rates(2,irg3a));

    // c12 + c12
    rate_c12c12(tf,bden,
                rr.rates(1,ir1212),rr.rates(2,ir1212),
                rrate,drratedt);

    // c12 + o16
    rate_c12o16(tf,bden,
                rr.rates(1,ir1216),rr.rates(2,ir1216),
                rrate,drratedt);

    // o16 + o16
    rate_o16o16(tf,bden,
                rr.rates(1,ir1616),rr.rates(2,ir1616),
                rrate,drratedt);

    // o16(a,g)ne20
    rate_o16ag(tf,bden,
               rr.rates(1,iroag),rr.rates(2,iroag),
               rr.rates(1,irnega),rr.rates(2,irnega));

    // ne20(a,g)mg24
    rate_ne20ag(tf,bden,
                rr.rates(1,irneag),rr.rates(2,irneag),
                rr.rates(1,irmgga),rr.rates(2,irmgga));

    // mg24(a,g)si28
    rate_mg24ag(tf,bden,
                rr.rates(1,irmgag),rr.rates(2,irmgag),
                rr.rates(1,irsiga),rr.rates(2,irsiga));

    // mg24(a,p)al27
    rate_mg24ap(tf,bden,
                rr.rates(1,irmgap),rr.rates(2,irmgap),
                rr.rates(1,iralpa),rr.rates(2,iralpa));

    // al27(p,g)si28
    rate_al27pg(tf,bden,
                rr.rates(1,iralpg),rr.rates(2,iralpg),
                rr.rates(1,irsigp),rr.rates(2,irsigp));

    // si28(a,g)s32
    rate_si28ag(tf,bden,
                rr.rates(1,irsiag),rr.rates(2,irsiag),
                rr.rates(1,irsga),rr.rates(2,irsga));

    // si28(a,p)p31
    rate_si28ap(tf,bden,
                rr.rates(1,irsiap),rr.rates(2,irsiap),
                rr.rates(1,irppa),rr.rates(2,irppa));

    // p31(p,g)s32
    rate_p31pg(tf,bden,
               rr.rates(1,irppg),rr.rates(2,irppg),
               rr.rates(1,irsgp),rr.rates(2,irsgp));

    // s32(a,g)ar36
    rate_s32ag(tf,bden,
               rr.rates(1,irsag),rr.rates(2,irsag),
               rr.rates(1,irarga),rr.rates(2,irarga));

    // s32(a,p)cl35
    rate_s32ap(tf,bden,
               rr.rates(1,irsap),rr.rates(2,irsap),
               rr.rates(1,irclpa),rr.rates(2,irclpa));

    // cl35(p,g)ar36
    rate_cl35pg(tf,bden,
                rr.rates(1,irclpg),rr.rates(2,irclpg),
                rr.rates(1,irargp),rr.rates(2,irargp));

    // ar36(a,g)ca40
    rate_ar36ag(tf,bden,
                rr.rates(1,irarag),rr.rates(2,irarag),
                rr.rates(1,ircaga),rr.rates(2,ircaga));

    // ar36(a,p)k39
    rate_ar36ap(tf,bden,
                rr.rates(1,irarap),rr.rates(2,irarap),
                rr.rates(1,irkpa),rr.rates(2,irkpa));

    // k39(p,g)ca40
    rate_k39pg(tf,bden,
               rr.rates(1,irkpg),rr.rates(2,irkpg),
               rr.rates(1,ircagp),rr.rates(2,ircagp));

    // ca40(a,g)ti44
    rate_ca40ag(tf,bden,
                rr.rates(1,ircaag),rr.rates(2,ircaag),
                rr.rates(1,irtiga),rr.rates(2,irtiga));

    // ca40(a,p)sc43
    rate_ca40ap(tf,bden,
                rr.rates(1,ircaap),rr.rates(2,ircaap),
                rr.rates(1,irscpa),rr.rates(2,irscpa));

    // sc43(p,g)ti44
    rate_sc43pg(tf,bden,
                rr.rates(1,irscpg),rr.rates(2,irscpg),
                rr.rates(1,irtigp),rr.rates(2,irtigp));

    // ti44(a,g)cr48
    rate_ti44ag(tf,bden,
                rr.rates(1,irtiag),rr.rates(2,irtiag),
                rr.rates(1,ircrga),rr.rates(2,ircrga));

    // ti44(a,p)v47
    rate_ti44ap(tf,bden,
                rr.rates(1,irtiap),rr.rates(2,irtiap), 
                rr.rates(1,irvpa),rr.rates(2,irvpa));

    // v47(p,g)cr48
    rate_v47pg(tf,bden,
               rr.rates(1,irvpg),rr.rates(2,irvpg), 
               rr.rates(1,ircrgp),rr.rates(2,ircrgp));

    // cr48(a,g)fe52
    rate_cr48ag(tf,bden,
                rr.rates(1,ircrag),rr.rates(2,ircrag),
                rr.rates(1,irfega),rr.rates(2,irfega));

    // cr48(a,p)mn51
    rate_cr48ap(tf,bden,
                rr.rates(1,ircrap),rr.rates(2,ircrap),
                rr.rates(1,irmnpa),rr.rates(2,irmnpa));

    // mn51(p,g)fe52
    rate_mn51pg(tf,bden,
                rr.rates(1,irmnpg),rr.rates(2,irmnpg),
                rr.rates(1,irfegp),rr.rates(2,irfegp));

    // fe52(a,g)ni56
    rate_fe52ag(tf,bden,
                rr.rates(1,irfeag),rr.rates(2,irfeag),
                rr.rates(1,irniga),rr.rates(2,irniga));

    // fe52(a,p)co55
    rate_fe52ap(tf,bden,
                rr.rates(1,irfeap),rr.rates(2,irfeap),
                rr.rates(1,ircopa),rr.rates(2,ircopa));

    // co55(p,g)ni56
    rate_co55pg(tf,bden,
                rr.rates(1,ircopg),rr.rates(2,ircopg),
                rr.rates(1,irnigp),rr.rates(2,irnigp));
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void set_aprox13rat()
{
    using namespace RateTable;

    Real btemp;
    Real bden = 1.0e0_rt;
    rate_t rr;

    for (int i = 1; i <= tab_imax; ++i) {

       btemp = tab_tlo + static_cast<Real>(i-1) * tab_tstp;
       btemp = std::pow(10.0e0_rt, btemp);

       aprox13rat(btemp, bden, rr);

       ttab(i) = btemp;

       for (int j = 1; j <= Rates::NumRates; ++j) {

          rattab(j,i)    = rr.rates(1,j);
          drattabdt(j,i) = rr.rates(2,j);

       }
    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void screen_aprox13(const Real btemp, const Real bden,
                    Array1D<Real, 1, NumSpec> const& y,
                    rate_t& rr)
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
    Real sc1a,sc1adt,sc2a,sc2adt,sc3a,sc3adt;
    Real sc1add,sc2add;
    Real denom,denomdt,zz;
    Real ratraw;
    plasma_state_t state;

    // Set up the state data, which is the same for all screening factors.

    fill_plasma_state(state, btemp, bden, y);

    // first the always fun triple alpha and its inverse
    jscr = 0;
    screen5(state,jscr,
            zion[He4-1], aion[He4-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);

    jscr++;
    screen5(state,jscr,
            zion[He4-1], aion[He4-1], 4.0_rt, 8.0_rt,
            sc2a,sc2adt,sc2add);

    sc3a   = sc1a * sc2a;
    sc3adt = sc1adt*sc2a + sc1a*sc2adt;

    ratraw = rr.rates(1,ir3a);
    rr.rates(1,ir3a) = ratraw * sc3a;
    rr.rates(2,ir3a) = rr.rates(2,ir3a)*sc3a + ratraw*sc3adt;

    ratraw = rr.rates(1,irg3a);
    rr.rates(1,irg3a) = ratraw * sc3a;
    rr.rates(2,irg3a) = rr.rates(2,irg3a)*sc3a + ratraw*sc3adt;

    // c12 to o16
    // c12(a,g)o16
    jscr++;
    screen5(state,jscr,
            zion[C12-1], aion[C12-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);

    ratraw = rr.rates(1,ircag);
    rr.rates(1,ircag)  = ratraw * sc1a;
    rr.rates(2,ircag)  = rr.rates(2,ircag)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,iroga);
    rr.rates(1,iroga)  = ratraw * sc1a;
    rr.rates(2,iroga)  = rr.rates(2,iroga)*sc1a + ratraw*sc1adt;

    // c12 + c12
    jscr++;
    screen5(state,jscr,
            zion[C12-1], aion[C12-1], zion[C12-1], aion[C12-1],
            sc1a,sc1adt,sc1add);

    ratraw = rr.rates(1,ir1212);
    rr.rates(1,ir1212) = ratraw * sc1a;
    rr.rates(2,ir1212) = rr.rates(2,ir1212)*sc1a + ratraw*sc1adt;

    // c12 + o16
    jscr++;
    screen5(state,jscr,
            zion[C12-1], aion[C12-1], zion[O16-1], aion[O16-1], 
            sc1a,sc1adt,sc1add);

    ratraw = rr.rates(1,ir1216);
    rr.rates(1,ir1216) = ratraw * sc1a;
    rr.rates(2,ir1216) = rr.rates(2,ir1216)*sc1a + ratraw*sc1adt;

    // o16 + o16
    jscr++;
    screen5(state,jscr,
            zion[O16-1], aion[O16-1], zion[O16-1], aion[O16-1],
            sc1a,sc1adt,sc1add);

    ratraw = rr.rates(1,ir1616);
    rr.rates(1,ir1616) = ratraw * sc1a;
    rr.rates(2,ir1616) = rr.rates(2,ir1616)*sc1a + ratraw*sc1adt;

    // o16 to ne20
    jscr++;
    screen5(state,jscr,
            zion[O16-1], aion[O16-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);


    // o16(a,g)ne20
    ratraw = rr.rates(1,iroag);
    rr.rates(1,iroag) = ratraw * sc1a;
    rr.rates(2,iroag) = rr.rates(2,iroag)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irnega);
    rr.rates(1,irnega) = ratraw * sc1a;
    rr.rates(2,irnega) = rr.rates(2,irnega)*sc1a + ratraw*sc1adt;

    // ne20 to mg24
    jscr++;
    screen5(state,jscr,
            zion[Ne20-1], aion[Ne20-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);


    // ne20(a,g)mg24
    ratraw = rr.rates(1,irneag);
    rr.rates(1,irneag) = ratraw * sc1a;
    rr.rates(2,irneag) = rr.rates(2,irneag)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irmgga);
    rr.rates(1,irmgga) = ratraw * sc1a;
    rr.rates(2,irmgga) = rr.rates(2,irmgga)*sc1a + ratraw*sc1adt;

    // mg24 to si28
    jscr++;
    screen5(state,jscr,
            zion[Mg24-1], aion[Mg24-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);


    // mg24(a,g)si28
    ratraw = rr.rates(1,irmgag);
    rr.rates(1,irmgag) = ratraw * sc1a;
    rr.rates(2,irmgag) = rr.rates(2,irmgag)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irsiga);
    rr.rates(1,irsiga) = ratraw * sc1a;
    rr.rates(2,irsiga) = rr.rates(2,irsiga)*sc1a + ratraw*sc1adt;


    // mg24(a,p)al27
    ratraw = rr.rates(1,irmgap);
    rr.rates(1,irmgap) = ratraw * sc1a;
    rr.rates(2,irmgap) = rr.rates(2,irmgap)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,iralpa);
    rr.rates(1,iralpa) = ratraw * sc1a;
    rr.rates(2,iralpa) = rr.rates(2,iralpa)*sc1a + ratraw*sc1adt;

    jscr++;
    screen5(state,jscr,
            13.0_rt, 27.0_rt, 1.0_rt, 1.0_rt,
            sc1a,sc1adt,sc1add);

    // al27(p,g)si28
    ratraw = rr.rates(1,iralpg);
    rr.rates(1,iralpg) = ratraw * sc1a;
    rr.rates(2,iralpg) = rr.rates(2,iralpg)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irsigp);
    rr.rates(1,irsigp) = ratraw * sc1a;
    rr.rates(2,irsigp) = rr.rates(2,irsigp)*sc1a + ratraw*sc1adt;


    // si28 to s32
    jscr++;
    screen5(state,jscr,
            zion[Si28-1], aion[Si28-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);


    // si28(a,g)s32
    ratraw = rr.rates(1,irsiag);
    rr.rates(1,irsiag) = ratraw * sc1a;
    rr.rates(2,irsiag) = rr.rates(2,irsiag)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irsga);
    rr.rates(1,irsga) = ratraw * sc1a;
    rr.rates(2,irsga) = rr.rates(2,irsga)*sc1a + ratraw*sc1adt;


    // si28(a,p)p31
    ratraw = rr.rates(1,irsiap);
    rr.rates(1,irsiap) = ratraw * sc1a;
    rr.rates(2,irsiap) = rr.rates(2,irsiap)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irppa);
    rr.rates(1,irppa)  = ratraw * sc1a;
    rr.rates(2,irppa)  = rr.rates(2,irppa)*sc1a  + ratraw*sc1adt;


    jscr++;
    screen5(state,jscr,
            15.0_rt, 31.0_rt, 1.0_rt, 1.0_rt,
            sc1a,sc1adt,sc1add);

    // p31(p,g)s32
    ratraw = rr.rates(1,irppg);
    rr.rates(1,irppg)  = ratraw * sc1a;
    rr.rates(2,irppg)  = rr.rates(2,irppg)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irsgp);
    rr.rates(1,irsgp)  = ratraw * sc1a;
    rr.rates(2,irsgp)  = rr.rates(2,irsgp)*sc1a + ratraw*sc1adt;


    // s32 to ar36
    jscr++;
    screen5(state,jscr,
            zion[S32-1], aion[S32-1], zion[He4-1], aion[He4-1], 
            sc1a,sc1adt,sc1add);


    // s32(a,g)ar36
    ratraw = rr.rates(1,irsag);
    rr.rates(1,irsag)  = ratraw * sc1a;
    rr.rates(2,irsag)  = rr.rates(2,irsag)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irarga);
    rr.rates(1,irarga)  = ratraw * sc1a;
    rr.rates(2,irarga)  = rr.rates(2,irarga)*sc1a + ratraw*sc1adt;

    // s32(a,p)cl35
    ratraw = rr.rates(1,irsap);
    rr.rates(1,irsap)  = ratraw * sc1a;
    rr.rates(2,irsap)  = rr.rates(2,irsap)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irclpa);
    rr.rates(1,irclpa) = ratraw * sc1a;
    rr.rates(2,irclpa) = rr.rates(2,irclpa)*sc1a + ratraw*sc1adt;


    jscr++;
    screen5(state,jscr,
            17.0_rt, 35.0_rt, 1.0_rt, 1.0_rt,
            sc1a,sc1adt,sc1add);

    // cl35(p,g)ar36
    ratraw = rr.rates(1,irclpg);
    rr.rates(1,irclpg) = ratraw * sc1a;
    rr.rates(2,irclpg) = rr.rates(2,irclpg)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irargp);
    rr.rates(1,irargp) = ratraw * sc1a;
    rr.rates(2,irargp) = rr.rates(2,irargp)*sc1a + ratraw*sc1adt;


    // ar36 to ca40
    jscr++;
    screen5(state,jscr,
            zion[Ar36-1], aion[Ar36-1], zion[He4-1], aion[He4-1], 
            sc1a,sc1adt,sc1add);


    // ar36(a,g)ca40
    ratraw = rr.rates(1,irarag);
    rr.rates(1,irarag) = ratraw * sc1a;
    rr.rates(2,irarag) = rr.rates(2,irarag)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,ircaga);
    rr.rates(1,ircaga) = ratraw * sc1a;
    rr.rates(2,ircaga) = rr.rates(2,ircaga)*sc1a + ratraw*sc1adt;


    // ar36(a,p)k39
    ratraw = rr.rates(1,irarap);
    rr.rates(1,irarap) = ratraw * sc1a;
    rr.rates(2,irarap) = rr.rates(2,irarap)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irkpa);
    rr.rates(1,irkpa) = ratraw * sc1a;
    rr.rates(2,irkpa) = rr.rates(2,irkpa)*sc1a  + ratraw*sc1adt;


    jscr++;
    screen5(state,jscr,
            19.0_rt, 39.0_rt, 1.0_rt, 1.0_rt,
            sc1a,sc1adt,sc1add);

    // k39(p,g)ca40
    ratraw = rr.rates(1,irkpg);
    rr.rates(1,irkpg) = ratraw * sc1a;
    rr.rates(2,irkpg) = rr.rates(2,irkpg)*sc1a  + ratraw*sc1adt;

    ratraw = rr.rates(1,ircagp);
    rr.rates(1,ircagp) = ratraw * sc1a;
    rr.rates(2,ircagp) = rr.rates(2,ircagp)*sc1a  + ratraw*sc1adt;


    // ca40 to ti44
    jscr++;
    screen5(state,jscr,
            zion[Ca40-1], aion[Ca40-1], zion[He4-1], aion[He4-1], 
            sc1a,sc1adt,sc1add);


    // ca40(a,g)ti44
    ratraw = rr.rates(1,ircaag);
    rr.rates(1,ircaag) = ratraw * sc1a;
    rr.rates(2,ircaag) = rr.rates(2,ircaag)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irtiga);
    rr.rates(1,irtiga) = ratraw * sc1a;
    rr.rates(2,irtiga) = rr.rates(2,irtiga)*sc1a + ratraw*sc1adt;


    // ca40(a,p)sc43
    ratraw = rr.rates(1,ircaap);
    rr.rates(1,ircaap) = ratraw * sc1a;
    rr.rates(2,ircaap) = rr.rates(2,ircaap)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irscpa);
    rr.rates(1,irscpa) = ratraw * sc1a;
    rr.rates(2,irscpa) = rr.rates(2,irscpa)*sc1a + ratraw*sc1adt;


    jscr++;
    screen5(state,jscr,
            21.0_rt, 43.0_rt, 1.0_rt, 1.0_rt,
            sc1a,sc1adt,sc1add);

    // sc43(p,g)ti44
    ratraw = rr.rates(1,irscpg);
    rr.rates(1,irscpg) = ratraw * sc1a;
    rr.rates(2,irscpg) = rr.rates(2,irscpg)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irtigp);
    rr.rates(1,irtigp) = ratraw * sc1a;
    rr.rates(2,irtigp) = rr.rates(2,irtigp)*sc1a + ratraw*sc1adt;


    // ti44 to cr48
    jscr++;
    screen5(state,jscr,
            zion[Ti44-1], aion[Ti44-1], zion[He4-1], aion[He4-1], 
            sc1a,sc1adt,sc1add);


    // ti44(a,g)cr48
    ratraw = rr.rates(1,irtiag);
    rr.rates(1,irtiag) = ratraw * sc1a;
    rr.rates(2,irtiag) = rr.rates(2,irtiag)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,ircrga);
    rr.rates(1,ircrga) = ratraw * sc1a;
    rr.rates(2,ircrga) = rr.rates(2,ircrga)*sc1a + ratraw*sc1adt;

    // ti44(a,p)v47
    ratraw = rr.rates(1,irtiap);
    rr.rates(1,irtiap) = ratraw * sc1a;
    rr.rates(2,irtiap) = rr.rates(2,irtiap)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irvpa);
    rr.rates(1,irvpa) = ratraw * sc1a;
    rr.rates(2,irvpa) = rr.rates(2,irvpa)*sc1a  + ratraw*sc1adt;


    jscr++;
    screen5(state,jscr,
            23.0_rt, 47.0_rt, 1.0_rt, 1.0_rt,
            sc1a,sc1adt,sc1add);

    // v47(p,g)cr48
    ratraw = rr.rates(1,irvpg);
    rr.rates(1,irvpg) = ratraw * sc1a;
    rr.rates(2,irvpg) = rr.rates(2,irvpg)*sc1a  + ratraw*sc1adt;

    ratraw = rr.rates(1,ircrgp);
    rr.rates(1,ircrgp)  = ratraw * sc1a;
    rr.rates(2,ircrgp)  = rr.rates(2,ircrgp)*sc1a  + ratraw*sc1adt;


    // cr48 to fe52
    jscr++;
    screen5(state,jscr,
            zion[Cr48-1], aion[Cr48-1], zion[He4-1], aion[He4-1], 
            sc1a,sc1adt,sc1add);


    // cr48(a,g)fe52
    ratraw = rr.rates(1,ircrag);
    rr.rates(1,ircrag) = ratraw * sc1a;
    rr.rates(2,ircrag) = rr.rates(2,ircrag)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irfega);
    rr.rates(1,irfega) = ratraw * sc1a;
    rr.rates(2,irfega) = rr.rates(2,irfega)*sc1a + ratraw*sc1adt;


    // cr48(a,p)mn51
    ratraw = rr.rates(1,ircrap);
    rr.rates(1,ircrap) = ratraw * sc1a;
    rr.rates(2,ircrap) = rr.rates(2,ircrap)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irmnpa);
    rr.rates(1,irmnpa) = ratraw * sc1a;
    rr.rates(2,irmnpa) = rr.rates(2,irmnpa)*sc1a + ratraw*sc1adt;


    jscr++;
    screen5(state,jscr,
            25.0_rt, 51.0_rt, 1.0_rt, 1.0_rt,
            sc1a,sc1adt,sc1add);

    // mn51(p,g)fe52
    ratraw = rr.rates(1,irmnpg);
    rr.rates(1,irmnpg) = ratraw * sc1a;
    rr.rates(2,irmnpg) = rr.rates(2,irmnpg)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irfegp);
    rr.rates(1,irfegp) = ratraw * sc1a;
    rr.rates(2,irfegp) = rr.rates(2,irfegp)*sc1a + ratraw*sc1adt;


    // fe52 to ni56
    jscr++;
    screen5(state,jscr,
            zion[Fe52-1], aion[Fe52-1], zion[He4-1], aion[He4-1], 
            sc1a,sc1adt,sc1add);


    // fe52(a,g)ni56
    ratraw = rr.rates(1,irfeag);
    rr.rates(1,irfeag) = ratraw * sc1a;
    rr.rates(2,irfeag) = rr.rates(2,irfeag)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irniga);
    rr.rates(1,irniga) = ratraw * sc1a;
    rr.rates(2,irniga) = rr.rates(2,irniga)*sc1a + ratraw*sc1adt;


    // fe52(a,p)co55
    ratraw = rr.rates(1,irfeap);
    rr.rates(1,irfeap) = ratraw * sc1a;
    rr.rates(2,irfeap) = rr.rates(2,irfeap)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,ircopa);
    rr.rates(1,ircopa) = ratraw * sc1a;
    rr.rates(2,ircopa) = rr.rates(2,ircopa)*sc1a + ratraw*sc1adt;


    jscr++;
    screen5(state,jscr,
            27.0_rt, 55.0_rt, 1.0_rt, 1.0_rt,
            sc1a,sc1adt,sc1add);


    // co55(p,g)ni56
    ratraw = rr.rates(1,ircopg);
    rr.rates(1,ircopg) = ratraw * sc1a;
    rr.rates(2,ircopg) = rr.rates(2,ircopg)*sc1a + ratraw*sc1adt;

    ratraw = rr.rates(1,irnigp);
    rr.rates(1,irnigp) = ratraw * sc1a;
    rr.rates(2,irnigp) = rr.rates(2,irnigp)*sc1a + ratraw*sc1adt;


    // now form those lovely dummy proton link rates

    // mg24(a,p)27al(p,g)28si
    rr.rates(1,irr1)  = 0.0e0_rt;
    rr.rates(2,irr1)  = 0.0e0_rt;
    denom    = rr.rates(1,iralpa) + rr.rates(1,iralpg);
    denomdt  = rr.rates(2,iralpa) + rr.rates(2,iralpg);

    if (denom > 1.0e-30_rt) {
       zz = 1.0e0_rt/denom;
       rr.rates(1,irr1) = rr.rates(1,iralpa)*zz;
       rr.rates(2,irr1) = (rr.rates(2,iralpa) - rr.rates(1,irr1)*denomdt)*zz;
    }

    // si28(a,p)p31(p,g)s32
    rr.rates(1,irs1)  = 0.0e0_rt;
    rr.rates(2,irs1)  = 0.0e0_rt;
    denom    = rr.rates(1,irppa) + rr.rates(1,irppg);
    denomdt  = rr.rates(2,irppa) + rr.rates(2,irppg);
    if (denom > 1.0e-30_rt) {
       zz = 1.0e0_rt/denom;
       rr.rates(1,irs1) = rr.rates(1,irppa)*zz;
       rr.rates(2,irs1) = (rr.rates(2,irppa) - rr.rates(1,irs1)*denomdt)*zz;
    }

    // s32(a,p)cl35(p,g)ar36
    rr.rates(1,irt1)  = 0.0e0_rt;
    rr.rates(2,irt1)  = 0.0e0_rt;
    denom    = rr.rates(1,irclpa) + rr.rates(1,irclpg);
    denomdt  = rr.rates(2,irclpa) + rr.rates(2,irclpg);
    if (denom > 1.0e-30_rt) {
       zz = 1.0e0_rt/denom;
       rr.rates(1,irt1) = rr.rates(1,irclpa)*zz;
       rr.rates(2,irt1) = (rr.rates(2,irclpa) - rr.rates(1,irt1)*denomdt)*zz;
    }

    // ar36(a,p)k39(p,g)ca40
    rr.rates(1,iru1)  = 0.0e0_rt;
    rr.rates(2,iru1)  = 0.0e0_rt;
    denom    = rr.rates(1,irkpa) + rr.rates(1,irkpg);
    denomdt  = rr.rates(2,irkpa) + rr.rates(2,irkpg);
    if (denom > 1.0e-30_rt) {
       zz   = 1.0e0_rt/denom;
       rr.rates(1,iru1)   = rr.rates(1,irkpa)*zz;
       rr.rates(2,iru1) = (rr.rates(2,irkpa) - rr.rates(1,iru1)*denomdt)*zz;
    }

    // ca40(a,p)sc43(p,g)ti44
    rr.rates(1,irv1)  = 0.0e0_rt;
    rr.rates(2,irv1)  = 0.0e0_rt;
    denom    = rr.rates(1,irscpa) + rr.rates(1,irscpg);
    denomdt  = rr.rates(2,irscpa) + rr.rates(2,irscpg);
    if (denom > 1.0e-30_rt) {
       zz  = 1.0e0_rt/denom;
       rr.rates(1,irv1) = rr.rates(1,irscpa)*zz;
       rr.rates(2,irv1) = (rr.rates(2,irscpa) - rr.rates(1,irv1)*denomdt)*zz;
    }

    // ti44(a,p)v47(p,g)cr48
    rr.rates(1,irw1) = 0.0e0_rt;
    rr.rates(2,irw1) = 0.0e0_rt;
    denom    = rr.rates(1,irvpa) + rr.rates(1,irvpg);
    denomdt  = rr.rates(2,irvpa) + rr.rates(2,irvpg);
    if (denom > 1.0e-30_rt) {
       zz = 1.0e0_rt/denom;
       rr.rates(1,irw1) = rr.rates(1,irvpa)*zz;
       rr.rates(2,irw1) = (rr.rates(2,irvpa) - rr.rates(1,irw1)*denomdt)*zz;
    }

    // cr48(a,p)mn51(p,g)fe52
    rr.rates(1,irx1) = 0.0e0_rt;
    rr.rates(2,irx1) = 0.0e0_rt;
    denom    = rr.rates(1,irmnpa) + rr.rates(1,irmnpg);
    denomdt  = rr.rates(2,irmnpa) + rr.rates(2,irmnpg);
    if (denom > 1.0e-30_rt) {
       zz = 1.0e0_rt/denom;
       rr.rates(1,irx1) = rr.rates(1,irmnpa)*zz;
       rr.rates(2,irx1) = (rr.rates(2,irmnpa) - rr.rates(1,irx1)*denomdt)*zz;
    }

    // fe52(a,p)co55(p,g)ni56
    rr.rates(1,iry1) = 0.0e0_rt;
    rr.rates(2,iry1) = 0.0e0_rt;
    denom    = rr.rates(1,ircopa) + rr.rates(1,ircopg);
    denomdt  = rr.rates(2,ircopa) + rr.rates(2,ircopg);
    if (denom > 1.0e-30_rt) {
       zz = 1.0e0_rt/denom;
       rr.rates(1,iry1) = rr.rates(1,ircopa)*zz;
       rr.rates(2,iry1) = (rr.rates(2,ircopa) - rr.rates(1,iry1)*denomdt)*zz;
    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void evaluate_rates(burn_t const& state, rate_t& rr)
{
    Real rho, temp;
    Array1D<Real, 1, NumSpec> y;

    // Get the data from the state
    rho  = state.rho;
    temp = state.T;

    for (int i = 1; i <= NumSpec; ++i) {
        y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    // Get the raw reaction rates
    if (use_tables) {
        aprox13tab(temp, rho, rr);
    } else {
        aprox13rat(temp, rho, rr);
    }

    // Do the screening here because the corrections depend on the composition
    screen_aprox13(temp, rho, y, rr);

    rr.T_eval = temp;
}


template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void dfdy_isotopes_aprox13(Array1D<Real, 1, NumSpec> const& y,
                           burn_t const& state, rate_t const& rr,
                           MatrixType& jac)
{

    // this routine sets up the dense aprox13 jacobian for the isotopes

    using namespace Species;
    using namespace Rates;

    {
        // he4 jacobian elements
        // d(he4)/d(he4)
        Array1D<Real, 1, 20> b {
            -1.5e0_rt * y(He4) * y(He4) * rr.rates(1,ir3a),
            -y(C12)  * rr.rates(1,ircag),
            -y(O16)  * rr.rates(1,iroag),
            -y(Ne20) * rr.rates(1,irneag),
            -y(Mg24) * rr.rates(1,irmgag),
            -y(Si28) * rr.rates(1,irsiag),
            -y(S32)  * rr.rates(1,irsag),
            -y(Ar36) * rr.rates(1,irarag),
            -y(Ca40) * rr.rates(1,ircaag),
            -y(Ti44) * rr.rates(1,irtiag),
            -y(Cr48) * rr.rates(1,ircrag),
            -y(Fe52) * rr.rates(1,irfeag),
            -y(Mg24) * rr.rates(1,irmgap) * (1.0e0_rt-rr.rates(1,irr1)),
            -y(Si28) * rr.rates(1,irsiap) * (1.0e0_rt-rr.rates(1,irs1)),
            -y(S32)  * rr.rates(1,irsap)  * (1.0e0_rt-rr.rates(1,irt1)),
            -y(Ar36) * rr.rates(1,irarap) * (1.0e0_rt-rr.rates(1,iru1)),
            -y(Ca40) * rr.rates(1,ircaap) * (1.0e0_rt-rr.rates(1,irv1)),
            -y(Ti44) * rr.rates(1,irtiap) * (1.0e0_rt-rr.rates(1,irw1)),
            -y(Cr48) * rr.rates(1,ircrap) * (1.0e0_rt-rr.rates(1,irx1)),
            -y(Fe52) * rr.rates(1,irfeap) * (1.0e0_rt-rr.rates(1,iry1))
        };

        jac(He4,He4) = esum20(b);
    }

    {
        // (he4)/d(o16)
        Array1D<Real, 1, 5> b {
            0.5e0_rt * y(C12) * rr.rates(1,ir1216),
            1.12e0_rt * 0.5e0_rt*y(O16) * rr.rates(1,ir1616),
            0.68e0_rt * rr.rates(1,irs1) * 0.5e0_rt*y(O16) * rr.rates(1,ir1616),
            rr.rates(1,iroga),
            -y(He4) * rr.rates(1,iroag),
        };

        jac(He4,O16) = esum5(b);
    }

    {
        Array1D<Real, 1, 4> b;

        // d(he4)/d(c12)
        b(1) =  y(C12) * rr.rates(1,ir1212);
        b(2) =  0.5e0_rt * y(O16) * rr.rates(1,ir1216);
        b(3) =  3.0e0_rt * rr.rates(1,irg3a);
        b(4) = -y(He4) * rr.rates(1,ircag);

        jac(He4,C12) = esum4(b);

        // d(he4)/d(si28)
        b(1) =  rr.rates(1,irsiga);
        b(2) = -y(He4) * rr.rates(1,irsiag);
        b(3) = -y(He4) * rr.rates(1,irsiap) * (1.0e0_rt-rr.rates(1,irs1));
        b(4) =  rr.rates(1,irr1) * rr.rates(1,irsigp);

        jac(He4,Si28) = esum4(b);

        // d(he4)/d(s32)
        b(1) =  rr.rates(1,irsga);
        b(2) = -y(He4) * rr.rates(1,irsag);
        b(3) = -y(He4) * rr.rates(1,irsap) * (1.0e0_rt-rr.rates(1,irt1));
        b(4) =  rr.rates(1,irs1) * rr.rates(1,irsgp);

        jac(He4,S32) = esum4(b);

        // d(he4)/d(ar36)
        b(1)  =  rr.rates(1,irarga);
        b(2)  = -y(He4) * rr.rates(1,irarag);
        b(3)  = -y(He4) * rr.rates(1,irarap) * (1.0e0_rt-rr.rates(1,iru1));
        b(4)  =  rr.rates(1,irt1) * rr.rates(1,irargp);

        jac(He4,Ar36) = esum4(b);

        // d(he4)/d(ca40)
        b(1) =  rr.rates(1,ircaga);
        b(2) = -y(He4) * rr.rates(1,ircaag);
        b(3) = -y(He4) * rr.rates(1,ircaap) * (1.0e0_rt-rr.rates(1,irv1));
        b(4) =  rr.rates(1,iru1) * rr.rates(1,ircagp);

        jac(He4,Ca40) = esum4(b);

        // d(he4)/d(ti44)
        b(1) =  rr.rates(1,irtiga);
        b(2) = -y(He4) * rr.rates(1,irtiag);
        b(3) = -y(He4) * rr.rates(1,irtiap) * (1.0e0_rt-rr.rates(1,irw1));
        b(4) =  rr.rates(1,irv1) * rr.rates(1,irtigp);

        jac(He4,Ti44) = esum4(b);

        // d(he4)/d(cr48)
        b(1) =  rr.rates(1,ircrga);
        b(2) = -y(He4) * rr.rates(1,ircrag);
        b(3) = -y(He4) * rr.rates(1,ircrap) * (1.0e0_rt-rr.rates(1,irx1));
        b(4) =  rr.rates(1,irw1) * rr.rates(1,ircrgp);

        jac(He4,Cr48) = esum4(b);

        // d(he4)/d(fe52)
        b(1) =  rr.rates(1,irfega);
        b(2) = -y(He4) * rr.rates(1,irfeag);
        b(3) = -y(He4) * rr.rates(1,irfeap) * (1.0e0_rt-rr.rates(1,iry1));
        b(4) =  rr.rates(1,irx1) * rr.rates(1,irfegp);

        jac(He4,Fe52) = esum4(b);

        // d(c12)/d(c12)
        b(1) = -2.0e0_rt * y(C12) * rr.rates(1,ir1212);
        b(2) = -y(O16) * rr.rates(1,ir1216);
        b(3) = -rr.rates(1,irg3a);
        b(4) = -y(He4) * rr.rates(1,ircag);

        jac(C12,C12) = esum4(b);

        // d(o16)/d(o16)
        b(1) = -y(C12) * rr.rates(1,ir1216);
        b(2) = -2.0e0_rt * y(O16) * rr.rates(1,ir1616);
        b(3) = -y(He4) * rr.rates(1,iroag);
        b(4) = -rr.rates(1,iroga);

        jac(O16,O16) = esum4(b);

        // si28 jacobian elements
        // d(si28)/d(he4)
        b(1) =  y(Mg24) * rr.rates(1,irmgag);
        b(2) = -y(Si28) * rr.rates(1,irsiag);
        b(3) =  y(Mg24) * rr.rates(1,irmgap) * (1.0e0_rt-rr.rates(1,irr1));
        b(4) = -y(Si28) * rr.rates(1,irsiap) * (1.0e0_rt-rr.rates(1,irs1));

        jac(Si28,He4) = esum4(b);

        // d(si28)/d(si28)
        b(1) =  -y(He4) * rr.rates(1,irsiag);
        b(2) = -rr.rates(1,irsiga);
        b(3) = -rr.rates(1,irr1) * rr.rates(1,irsigp);
        b(4) = -y(He4) * rr.rates(1,irsiap) * (1.0e0_rt-rr.rates(1,irs1));

        jac(Si28,Si28) = esum4(b);

        // s32 jacobian elements
        // d(s32)/d(he4)
        b(1) =  y(Si28) * rr.rates(1,irsiag);
        b(2) = -y(S32) * rr.rates(1,irsag);
        b(3) =  y(Si28) * rr.rates(1,irsiap) * (1.0e0_rt-rr.rates(1,irs1));
        b(4) = -y(S32) * rr.rates(1,irsap) * (1.0e0_rt-rr.rates(1,irt1));

        jac(S32,He4) = esum4(b);

        // d(s32)/d(s32)
        b(1) = -y(He4) * rr.rates(1,irsag);
        b(2) = -rr.rates(1,irsga);
        b(3) = -rr.rates(1,irs1) * rr.rates(1,irsgp);
        b(4) = -y(He4) * rr.rates(1,irsap) * (1.0e0_rt-rr.rates(1,irt1));

        jac(S32,S32) = esum4(b);

        // ar36 jacobian elements
        // d(ar36)/d(he4)
        b(1) =  y(S32)  * rr.rates(1,irsag);
        b(2) = -y(Ar36) * rr.rates(1,irarag);
        b(3) =  y(S32)  * rr.rates(1,irsap) * (1.0e0_rt-rr.rates(1,irt1));
        b(4) = -y(Ar36) * rr.rates(1,irarap) * (1.0e0_rt-rr.rates(1,iru1));

        jac(Ar36,He4) = esum4(b);

        // d(ar36)/d(ar36)
        b(1) = -y(He4) * rr.rates(1,irarag);
        b(2) = -rr.rates(1,irarga);
        b(3) = -rr.rates(1,irt1) * rr.rates(1,irargp);
        b(4) = -y(He4) * rr.rates(1,irarap) * (1.0e0_rt-rr.rates(1,iru1));

        jac(Ar36,Ar36) = esum4(b);

        // ca40 jacobian elements
        // d(ca40)/d(he4)
        b(1)  =  y(Ar36) * rr.rates(1,irarag);
        b(2)  = -y(Ca40) * rr.rates(1,ircaag);
        b(3)  =  y(Ar36) * rr.rates(1,irarap)*(1.0e0_rt-rr.rates(1,iru1));
        b(4)  = -y(Ca40) * rr.rates(1,ircaap)*(1.0e0_rt-rr.rates(1,irv1));

        jac(Ca40,He4) = esum4(b);

        // d(ca40)/d(ca40)
        b(1) = -y(He4) * rr.rates(1, ircaag);
        b(2) = -rr.rates(1, ircaga);
        b(3) = -rr.rates(1, ircagp) * rr.rates(1, iru1);
        b(4) = -y(He4) * rr.rates(1, ircaap)*(1.0e0_rt-rr.rates(1, irv1));

        jac(Ca40,Ca40) = esum4(b);

        // ti44 jacobian elements
        // d(ti44)/d(he4)
        b(1) =  y(Ca40) * rr.rates(1, ircaag);
        b(2) = -y(Ti44) * rr.rates(1, irtiag);
        b(3) =  y(Ca40) * rr.rates(1, ircaap)*(1.0e0_rt-rr.rates(1, irv1));
        b(4) = -y(Ti44) * rr.rates(1, irtiap)*(1.0e0_rt-rr.rates(1, irw1));

        jac(Ti44,He4) = esum4(b);

        // d(ti44)/d(ti44)
        b(1) = -y(He4) * rr.rates(1, irtiag);
        b(2) = -rr.rates(1, irtiga);
        b(3) = -rr.rates(1, irv1) * rr.rates(1, irtigp);
        b(4) = -y(He4) * rr.rates(1, irtiap)*(1.0e0_rt-rr.rates(1, irw1));

        jac(Ti44,Ti44) = esum4(b);

        // cr48 jacobian elements
        // d(cr48)/d(he4)
        b(1) =  y(Ti44) * rr.rates(1, irtiag);
        b(2) = -y(Cr48) * rr.rates(1, ircrag);
        b(3) =  y(Ti44) * rr.rates(1, irtiap)*(1.0e0_rt-rr.rates(1, irw1));
        b(4) = -y(Cr48) * rr.rates(1, ircrap)*(1.0e0_rt-rr.rates(1, irx1));

        jac(Cr48,He4) = esum4(b);

        // d(cr48)/d(cr48)
        b(1) = -y(He4) * rr.rates(1, ircrag);
        b(2) = -rr.rates(1, ircrga);
        b(3) = -rr.rates(1, irw1) * rr.rates(1, ircrgp);
        b(4) = -y(He4) * rr.rates(1, ircrap)*(1.0e0_rt-rr.rates(1, irx1));

        jac(Cr48,Cr48) = esum4(b);

        // fe52 jacobian elements
        // d(fe52)/d(he4)
        b(1) =  y(Cr48) * rr.rates(1, ircrag);
        b(2) = -y(Fe52) * rr.rates(1, irfeag);
        b(3) =  y(Cr48) * rr.rates(1, ircrap) * (1.0e0_rt-rr.rates(1, irx1));
        b(4) = -y(Fe52) * rr.rates(1, irfeap) * (1.0e0_rt-rr.rates(1, iry1));

        jac(Fe52,He4) = esum4(b);

        // d(fe52)/d(fe52)
        b(1) = -y(He4) * rr.rates(1, irfeag);
        b(2) = -rr.rates(1, irfega);
        b(3) = -rr.rates(1, irx1) * rr.rates(1, irfegp);
        b(4) = -y(He4) * rr.rates(1, irfeap) * (1.0e0_rt-rr.rates(1, iry1));

        jac(Fe52,Fe52) = esum4(b);
    }

    {
        Array1D<Real, 1, 3> b;

        // d(he4)/d(mg24)
        b(1) =  rr.rates(1,irmgga);
        b(2) = -y(He4) * rr.rates(1,irmgag);
        b(3) = -y(He4) * rr.rates(1,irmgap) * (1.0e0_rt-rr.rates(1,irr1));

        jac(He4,Mg24) = esum3(b);

        // mg24 jacobian elements
        // d(mg24)/d(he4)
        b(1) =  y(Ne20) * rr.rates(1,irneag);
        b(2) = -y(Mg24) * rr.rates(1,irmgag);
        b(3) = -y(Mg24) * rr.rates(1,irmgap) * (1.0e0_rt-rr.rates(1,irr1));

        jac(Mg24,He4) = esum3(b);

        // d(mg24)/d(mg24)
        b(1) = -y(He4) * rr.rates(1,irmgag);
        b(2) = -rr.rates(1,irmgga);
        b(3) = -y(He4) * rr.rates(1,irmgap) * (1.0e0_rt-rr.rates(1,irr1));

        jac(Mg24,Mg24) = esum3(b);

        // d(si28)/d(o16)
        b(1) = 0.5e0_rt * y(C12) * rr.rates(1,ir1216);
        b(2) = 1.12e0_rt * 0.5e0_rt*y(O16) * rr.rates(1,ir1616);
        b(3) = 0.68e0_rt * 0.5e0_rt*y(O16) * rr.rates(1,irs1) * rr.rates(1,ir1616);

        jac(Si28,O16) = esum3(b);
    }

    {
        Array1D<Real, 1, 2> b;

        // d(he4)/d(ne20)
        b(1) =  rr.rates(1,irnega);
        b(2) = -y(He4) * rr.rates(1,irneag);

        jac(He4,Ne20) = b(1) + b(2);

        // d(he4)/d(ni56)
        b(1) = rr.rates(1,irniga);
        b(2) = rr.rates(1,iry1) * rr.rates(1,irnigp);

        jac(He4,Ni56) = b(1) + b(2);

        // c12 jacobian elements
        // d(c12)/d(he4)
        b(1) =  0.5e0_rt * y(He4) * y(He4) * rr.rates(1,ir3a);
        b(2) = -y(C12) * rr.rates(1,ircag);

        jac(C12,He4) = b(1) + b(2);

        // d(c12)/d(o16)
        b(1) = -y(C12) * rr.rates(1,ir1216);
        b(2) =  rr.rates(1,iroga);

        jac(C12,O16) = b(1) + b(2);

        // o16 jacobian elements
        // d(o16)/d(he4)
        b(1) =  y(C12)*rr.rates(1,ircag);
        b(2) = -y(O16)*rr.rates(1,iroag);

        jac(O16,He4) = b(1) + b(2);

        // d(o16)/d(c12)
        b(1) = -y(O16)*rr.rates(1,ir1216);
        b(2) =  y(He4)*rr.rates(1,ircag);

        jac(O16,C12) = b(1) + b(2);

        // ne20 jacobian elements
        // d(ne20)/d(he4)
        b(1) =  y(O16) * rr.rates(1,iroag);
        b(2) = -y(Ne20) * rr.rates(1,irneag);

        jac(Ne20,He4) = b(1) + b(2);

        // d(ne20)/d(ne20)
        b(1) = -y(He4) * rr.rates(1,irneag);
        b(2) = -rr.rates(1,irnega);

        jac(Ne20,Ne20) = b(1) + b(2);

        // d(mg24)/d(si28)
        b(1) = rr.rates(1,irsiga);
        b(2) = rr.rates(1,irr1) * rr.rates(1,irsigp);

        jac(Mg24,Si28) = b(1) + b(2);

        // d(si28)/d(mg24)
        b(1) =  y(He4) * rr.rates(1,irmgag);
        b(2) =  y(He4) * rr.rates(1,irmgap) * (1.0e0_rt-rr.rates(1,irr1));

        jac(Si28,Mg24) = b(1) + b(2);

        // d(si28)/d(s32)
        b(1) = rr.rates(1,irsga);
        b(2) = rr.rates(1,irs1) * rr.rates(1,irsgp);

        jac(Si28,S32) = b(1) + b(2);

        // d(s32)/d(o16)
        b(1) = 0.68e0_rt*0.5e0_rt*y(O16)*rr.rates(1,ir1616)*(1.0e0_rt-rr.rates(1,irs1));
        b(2) = 0.2e0_rt * 0.5e0_rt*y(O16) * rr.rates(1,ir1616);

        jac(S32,O16) = b(1) + b(2);

        // d(s32)/d(si28)
        b(1)  =y(He4) * rr.rates(1,irsiag);
        b(2) = y(He4) * rr.rates(1,irsiap) * (1.0e0_rt-rr.rates(1,irs1));

        jac(S32,Si28) = b(1) + b(2);

        // d(s32)/d(ar36)
        b(1) = rr.rates(1,irarga);
        b(2) = rr.rates(1,irt1) * rr.rates(1,irargp);

        jac(S32,Ar36) = b(1) + b(2);

        // d(ar36)/d(s32)
        b(1) = y(He4) * rr.rates(1,irsag);
        b(2) = y(He4) * rr.rates(1,irsap) * (1.0e0_rt-rr.rates(1,irt1));

        jac(Ar36,S32) = b(1) + b(2);

        // d(ar36)/d(ca40)
        b(1) = rr.rates(1,ircaga);
        b(2) = rr.rates(1,ircagp) * rr.rates(1,iru1);

        jac(Ar36,Ca40) = b(1) + b(2);

        // d(ca40)/d(ar36)
        b(1) =  y(He4) * rr.rates(1,irarag);
        b(2) =  y(He4) * rr.rates(1,irarap)*(1.0e0_rt-rr.rates(1,iru1));

        jac(Ca40,Ar36) = b(1) + b(2);

        // d(ca40)/d(ti44)
        b(1) = rr.rates(1, irtiga);
        b(2) = rr.rates(1, irtigp) * rr.rates(1, irv1);

        jac(Ca40,Ti44) = b(1) + b(2);

        // d(ti44)/d(ca40)
        b(1) =  y(He4) * rr.rates(1, ircaag);
        b(2) =  y(He4) * rr.rates(1, ircaap)*(1.0e0_rt-rr.rates(1, irv1));

        jac(Ti44,Ca40) = b(1) + b(2);

        // d(ti44)/d(cr48)
        b(1) = rr.rates(1, ircrga);
        b(2) = rr.rates(1, irw1) * rr.rates(1, ircrgp);

        jac(Ti44,Cr48) = b(1) + b(2);

        // d(cr48)/d(ti44)
        b(1) =  y(He4) * rr.rates(1, irtiag);
        b(2) =  y(He4) * rr.rates(1, irtiap)*(1.0e0_rt-rr.rates(1, irw1));

        jac(Cr48,Ti44) = b(1) + b(2);

        // d(cr48)/d(fe52)
        b(1) = rr.rates(1, irfega);
        b(2) = rr.rates(1, irx1) * rr.rates(1, irfegp);

        jac(Cr48,Fe52) = b(1) + b(2);

        // d(fe52)/d(cr48)
        b(1) = y(He4) * rr.rates(1, ircrag);
        b(2) = y(He4) * rr.rates(1, ircrap) * (1.0e0_rt-rr.rates(1, irx1));

        jac(Fe52,Cr48) = b(1) + b(2);

        // d(fe52)/d(ni56)
        b(1) = rr.rates(1, irniga);
        b(2) = rr.rates(1, iry1) * rr.rates(1, irnigp);

        jac(Fe52,Ni56) = b(1) + b(2);

        // ni56 jacobian elements
        // d(ni56)/d(he4)
        b(1) =  y(Fe52) * rr.rates(1, irfeag);
        b(2) =  y(Fe52) * rr.rates(1, irfeap) * (1.0e0_rt-rr.rates(1, iry1));

        jac(Ni56,He4) = b(1) + b(2);

        // d(ni56)/d(fe52)
        b(1) = y(He4) * rr.rates(1, irfeag);
        b(2) = y(He4) * rr.rates(1, irfeap) * (1.0e0_rt-rr.rates(1, iry1));

        jac(Ni56,Fe52) = b(1) + b(2);

        // d(ni56)/d(ni56)
        b(1) = -rr.rates(1, irniga);
        b(2) = -rr.rates(1, iry1) * rr.rates(1, irnigp);

        jac(Ni56,Ni56) = b(1) + b(2);
    }

    // d(o16)/d(ne20)
    jac(O16,Ne20) = rr.rates(1,irnega);

    // d(ne20)/d(c12)
    jac(Ne20,C12) = y(C12) * rr.rates(1,ir1212);

    // d(ne20)/d(o16)
    jac(Ne20,O16) = y(He4) * rr.rates(1,iroag);

    // d(ne20)/d(mg24)
    jac(Ne20,Mg24) = rr.rates(1,irmgga);

    // d(mg24)/d(c12)
    jac(Mg24,C12) = 0.5e0_rt * y(O16) * rr.rates(1,ir1216);

    // d(mg24)/d(o16)
    jac(Mg24,O16) = 0.5e0_rt * y(C12) * rr.rates(1,ir1216);

    // d(mg24)/d(ne20)
    jac(Mg24,Ne20) = y(He4) * rr.rates(1,irneag);

    // d(si28)/d(c12)
    jac(Si28,C12) = 0.5e0_rt * y(O16) * rr.rates(1,ir1216);

}

template<class T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ener_gener_rate(T const& dydt, Real& enuc)
{

    using namespace aprox13;

    // Computes the instantaneous energy generation rate

    Real Xdot = 0.0_rt;

    // Sum the mass fraction time derivatives
    for (int i = 1; i <= NumSpec; ++i) {
        Xdot += dydt(i) * mion(i);
    }

    // This is basically e = m c**2
    enuc = Xdot * C::Legacy::enuc_conv2;

}


// Evaluates the right hand side of the aprox13 ODEs
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rhs(Array1D<Real, 1, NumSpec> const& y, rate_t const& rr,
         Array1D<Real, 1, neqs>& dydt,
         bool deriva, bool for_jacobian_tderiv)
{
    using namespace Species;
    using namespace Rates;

    // deriva is used in forming the analytic Jacobian to get
    // the derivative wrt A

    const int index_rate = for_jacobian_tderiv ? 2 : 1;

    for (int i = 1; i <= NumSpec; ++i) {
        dydt(i) = 0.0_rt;
    }

    Array1D<Real, 1, 17> a;

    // he4 reactions
    // heavy ion reactions
    a(1)  = 0.5e0_rt * y(C12) * y(C12) * rr.rates(index_rate, ir1212);
    a(2)  = 0.5e0_rt * y(C12) * y(O16) * rr.rates(index_rate, ir1216);
    a(3)  = 0.56e0_rt * 0.5e0_rt * y(O16) * y(O16) * rr.rates(index_rate, ir1616);

    dydt(He4) = dydt(He4) + esum3(a);

    // (a,g) and (g,a) reactions
    a(1)  = -0.5e0_rt * y(He4) * y(He4) * y(He4) * rr.rates(index_rate, ir3a);
    a(2)  =  3.0e0_rt * y(C12) * rr.rates(index_rate, irg3a);
    a(3)  = -y(He4)  * y(C12) * rr.rates(index_rate, ircag);
    a(4)  =  y(O16)  * rr.rates(index_rate, iroga);
    a(5)  = -y(He4)  * y(O16) * rr.rates(index_rate, iroag);
    a(6)  =  y(Ne20) * rr.rates(index_rate, irnega);
    a(7)  = -y(He4)  * y(Ne20) * rr.rates(index_rate, irneag);
    a(8)  =  y(Mg24) * rr.rates(index_rate, irmgga);
    a(9)  = -y(He4)  * y(Mg24)* rr.rates(index_rate, irmgag);
    a(10) =  y(Si28) * rr.rates(index_rate, irsiga);
    a(11) = -y(He4)  * y(Si28)*rr.rates(index_rate, irsiag);
    a(12) =  y(S32)  * rr.rates(index_rate, irsga);

    dydt(He4) = dydt(He4) + esum12(a);

    a(1)  = -y(He4)  * y(S32) * rr.rates(index_rate, irsag);
    a(2)  =  y(Ar36) * rr.rates(index_rate, irarga);
    a(3)  = -y(He4)  * y(Ar36)*rr.rates(index_rate, irarag);
    a(4)  =  y(Ca40) * rr.rates(index_rate, ircaga);
    a(5)  = -y(He4)  * y(Ca40)*rr.rates(index_rate, ircaag);
    a(6)  =  y(Ti44) * rr.rates(index_rate, irtiga);
    a(7)  = -y(He4)  * y(Ti44)*rr.rates(index_rate, irtiag);
    a(8)  =  y(Cr48) * rr.rates(index_rate, ircrga);
    a(9)  = -y(He4)  * y(Cr48)*rr.rates(index_rate, ircrag);
    a(10) =  y(Fe52) * rr.rates(index_rate, irfega);
    a(11) = -y(He4)  * y(Fe52) * rr.rates(index_rate, irfeag);
    a(12) =  y(Ni56) * rr.rates(index_rate, irniga);

    dydt(He4) = dydt(He4) + esum12(a);

    // (a,p)(p,g) and (g,p)(p,a) reactions

    if (!deriva) {

       a(1)  =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16)*rr.rates(index_rate, irs1)*rr.rates(index_rate, ir1616);
       a(2)  = -y(He4)  * y(Mg24) * rr.rates(index_rate, irmgap)*(1.0e0_rt-rr.rates(index_rate, irr1));
       a(3)  =  y(Si28) * rr.rates(index_rate, irsigp) * rr.rates(index_rate, irr1);
       a(4)  = -y(He4)  * y(Si28) * rr.rates(index_rate, irsiap)*(1.0e0_rt-rr.rates(index_rate, irs1));
       a(5)  =  y(S32)  * rr.rates(index_rate, irsgp) * rr.rates(index_rate, irs1);
       a(6)  = -y(He4)  * y(S32) * rr.rates(index_rate, irsap)*(1.0e0_rt-rr.rates(index_rate, irt1));
       a(7)  =  y(Ar36) * rr.rates(index_rate, irargp) * rr.rates(index_rate, irt1);
       a(8)  = -y(He4)  * y(Ar36) * rr.rates(index_rate, irarap)*(1.0e0_rt-rr.rates(index_rate, iru1));
       a(9)  =  y(Ca40) * rr.rates(index_rate, ircagp) * rr.rates(index_rate, iru1);
       a(10) = -y(He4)  * y(Ca40) * rr.rates(index_rate, ircaap)*(1.0e0_rt-rr.rates(index_rate, irv1));
       a(11) =  y(Ti44) * rr.rates(index_rate, irtigp) * rr.rates(index_rate, irv1);
       a(12) = -y(He4)  * y(Ti44) * rr.rates(index_rate, irtiap)*(1.0e0_rt-rr.rates(index_rate, irw1));
       a(13) =  y(Cr48) * rr.rates(index_rate, ircrgp) * rr.rates(index_rate, irw1);
       a(14) = -y(He4)  * y(Cr48) * rr.rates(index_rate, ircrap)*(1.0e0_rt-rr.rates(index_rate, irx1));
       a(15) =  y(Fe52) * rr.rates(index_rate, irfegp) * rr.rates(index_rate, irx1);
       a(16) = -y(He4)  * y(Fe52) * rr.rates(index_rate, irfeap)*(1.0e0_rt-rr.rates(index_rate, iry1));
       a(17) =  y(Ni56) * rr.rates(index_rate, irnigp) * rr.rates(index_rate, iry1);

       dydt(He4) = dydt(He4) + esum17(a);

    } else {

       a(1)  =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rr.rates(1, irs1) * rr.rates(index_rate, ir1616);
       a(2)  =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rr.rates(index_rate, irs1) * rr.rates(1, ir1616);
       a(3)  = -y(He4)*y(Mg24) * rr.rates(index_rate, irmgap)*(1.0e0_rt - rr.rates(1, irr1));
       a(4)  =  y(He4)*y(Mg24) * rr.rates(1, irmgap)*rr.rates(index_rate, irr1);
       a(5)  =  y(Si28) * rr.rates(1, irsigp) * rr.rates(index_rate, irr1);
       a(6)  =  y(Si28) * rr.rates(index_rate, irsigp) * rr.rates(1, irr1);
       a(7)  = -y(He4)*y(Si28) * rr.rates(index_rate, irsiap)*(1.0e0_rt - rr.rates(1, irs1));
       a(8)  =  y(He4)*y(Si28) * rr.rates(1, irsiap) * rr.rates(index_rate, irs1);
       a(9)  =  y(S32)  * rr.rates(1, irsgp) * rr.rates(index_rate, irs1);
       a(10) =  y(S32)  * rr.rates(index_rate, irsgp) * rr.rates(1, irs1);

       dydt(He4) = dydt(He4) + esum10(a);

       a(1)  = -y(He4)*y(S32) * rr.rates(index_rate, irsap)*(1.0e0_rt - rr.rates(1, irt1));
       a(2)  =  y(He4)*y(S32) * rr.rates(1, irsap)*rr.rates(index_rate, irt1);
       a(3)  =  y(Ar36) * rr.rates(1, irargp) * rr.rates(index_rate, irt1);
       a(4)  =  y(Ar36) * rr.rates(index_rate, irargp) * rr.rates(1, irt1);
       a(5)  = -y(He4)*y(Ar36) * rr.rates(index_rate, irarap)*(1.0e0_rt - rr.rates(1, iru1));
       a(6)  =  y(He4)*y(Ar36) * rr.rates(1, irarap)*rr.rates(index_rate, iru1);
       a(7)  =  y(Ca40) * rr.rates(1, ircagp) * rr.rates(index_rate, iru1);
       a(8)  =  y(Ca40) * rr.rates(index_rate, ircagp) * rr.rates(1, iru1);
       a(9)  = -y(He4)*y(Ca40) * rr.rates(index_rate, ircaap)*(1.0e0_rt-rr.rates(1, irv1));
       a(10) =  y(He4)*y(Ca40) * rr.rates(1, ircaap)*rr.rates(index_rate, irv1);
       a(11) =  y(Ti44) * rr.rates(1, irtigp) * rr.rates(index_rate, irv1);
       a(12) =  y(Ti44) * rr.rates(index_rate, irtigp) * rr.rates(1, irv1);

       dydt(He4) = dydt(He4) + esum12(a);

       a(1)  = -y(He4)*y(Ti44) * rr.rates(index_rate, irtiap)*(1.0e0_rt - rr.rates(1, irw1));
       a(2)  =  y(He4)*y(Ti44) * rr.rates(1, irtiap)*rr.rates(index_rate, irw1);
       a(3)  =  y(Cr48) * rr.rates(1, ircrgp) * rr.rates(index_rate, irw1);
       a(4)  =  y(Cr48) * rr.rates(index_rate, ircrgp) * rr.rates(1, irw1);
       a(5)  = -y(He4)*y(Cr48) * rr.rates(index_rate, ircrap)*(1.0e0_rt - rr.rates(1, irx1));
       a(6)  =  y(He4)*y(Cr48) * rr.rates(1, ircrap)*rr.rates(index_rate, irx1);
       a(7)  =  y(Fe52) * rr.rates(1, irfegp) * rr.rates(index_rate, irx1);
       a(8)  =  y(Fe52) * rr.rates(index_rate, irfegp) * rr.rates(1, irx1);
       a(9)  = -y(He4)*y(Fe52) * rr.rates(index_rate, irfeap)*(1.0e0_rt - rr.rates(1, iry1));
       a(10) =  y(He4)*y(Fe52) * rr.rates(1, irfeap)*rr.rates(index_rate, iry1);
       a(11) =  y(Ni56) * rr.rates(1, irnigp) * rr.rates(index_rate, iry1);
       a(12) =  y(Ni56) * rr.rates(index_rate, irnigp) * rr.rates(1, iry1);

       dydt(He4) = dydt(He4) + esum12(a);
    
    }

    // c12 reactions
    a(1) = -y(C12) * y(C12) * rr.rates(index_rate, ir1212);
    a(2) = -y(C12) * y(O16) * rr.rates(index_rate, ir1216);
    a(3) =  y(He4) * y(He4) * y(He4) * rr.rates(index_rate, ir3a) / 6.0_rt;;
    a(4) = -y(C12) * rr.rates(index_rate, irg3a);
    a(5) = -y(C12) * y(He4) * rr.rates(index_rate, ircag);
    a(6) =  y(O16) * rr.rates(index_rate, iroga);

    dydt(C12) = dydt(C12) + esum6(a);


    // o16 reactions
    a(1) = -y(C12) * y(O16) * rr.rates(index_rate, ir1216);
    a(2) = -y(O16) * y(O16) * rr.rates(index_rate, ir1616);
    a(3) =  y(C12) * y(He4) * rr.rates(index_rate, ircag);
    a(4) = -y(O16) * y(He4) * rr.rates(index_rate, iroag);
    a(5) = -y(O16) * rr.rates(index_rate, iroga);
    a(6) =  y(Ne20) * rr.rates(index_rate, irnega);

    dydt(O16) = dydt(O16) + esum6(a);


    // ne20 reactions
    a(1) =  0.5e0_rt * y(C12) * y(C12) * rr.rates(index_rate, ir1212);
    a(2) =  y(O16) * y(He4) * rr.rates(index_rate, iroag);
    a(3) = -y(Ne20) * y(He4) * rr.rates(index_rate, irneag);
    a(4) = -y(Ne20) * rr.rates(index_rate, irnega);
    a(5) =  y(Mg24) * rr.rates(index_rate, irmgga);

    dydt(Ne20) = dydt(Ne20) + esum5(a);


    // mg24 reactions
    a(1) =  0.5e0_rt * y(C12) * y(O16) * rr.rates(index_rate, ir1216);
    a(2) =  y(Ne20) * y(He4) * rr.rates(index_rate, irneag);
    a(3) = -y(Mg24) * y(He4) * rr.rates(index_rate, irmgag);
    a(4) = -y(Mg24) * rr.rates(index_rate, irmgga);
    a(5) =  y(Si28) * rr.rates(index_rate, irsiga);

    dydt(Mg24) = dydt(Mg24) + esum5(a);

    if (!deriva) {

       a(1) = -y(Mg24) * y(He4) * rr.rates(index_rate, irmgap)*(1.0e0_rt-rr.rates(index_rate, irr1));
       a(2) =  y(Si28) * rr.rates(index_rate, irr1) * rr.rates(index_rate, irsigp);

       dydt(Mg24) = dydt(Mg24) + a(1) + a(2);

    } else {

       a(1) = -y(Mg24)*y(He4) * rr.rates(index_rate, irmgap)*(1.0e0_rt - rr.rates(1, irr1));
       a(2) =  y(Mg24)*y(He4) * rr.rates(1, irmgap)*rr.rates(index_rate, irr1);
       a(3) =  y(Si28) * rr.rates(1, irr1) * rr.rates(index_rate, irsigp);
       a(4) =  y(Si28) * rr.rates(index_rate, irr1) * rr.rates(1, irsigp);

       dydt(Mg24) = dydt(Mg24) + esum4(a);

    }


    // si28 reactions
    a(1) =  0.5e0_rt * y(C12) * y(O16) * rr.rates(index_rate, ir1216);
    a(2) =  0.56e0_rt * 0.5e0_rt*y(O16) * y(O16) * rr.rates(index_rate, ir1616);
    a(3) =  y(Mg24) * y(He4) * rr.rates(index_rate, irmgag);
    a(4) = -y(Si28) * y(He4) * rr.rates(index_rate, irsiag);
    a(5) = -y(Si28) * rr.rates(index_rate, irsiga);
    a(6) =  y(S32)  * rr.rates(index_rate, irsga);

    dydt(Si28) = dydt(Si28) + esum6(a);

    if (!deriva) {

       a(1) =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16)*rr.rates(index_rate, irs1)*rr.rates(index_rate, ir1616);
       a(2) =  y(Mg24) * y(He4) * rr.rates(index_rate, irmgap)*(1.0e0_rt-rr.rates(index_rate, irr1));
       a(3) = -y(Si28) * rr.rates(index_rate, irr1) * rr.rates(index_rate, irsigp);
       a(4) = -y(Si28) * y(He4) * rr.rates(index_rate, irsiap)*(1.0e0_rt-rr.rates(index_rate, irs1));
       a(5) =  y(S32)  * rr.rates(index_rate, irs1) * rr.rates(index_rate, irsgp);

       dydt(Si28) = dydt(Si28) + esum5(a);

    } else {

       a(1)  =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rr.rates(1, irs1)*rr.rates(index_rate, ir1616);
       a(2)  =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rr.rates(index_rate, irs1)*rr.rates(1, ir1616);
       a(3)  =  y(Mg24)*y(He4) * rr.rates(index_rate, irmgap)*(1.0e0_rt - rr.rates(1, irr1));
       a(4)  = -y(Mg24)*y(He4) * rr.rates(1, irmgap)*rr.rates(index_rate, irr1);
       a(5)  = -y(Si28) * rr.rates(1, irr1) * rr.rates(index_rate, irsigp);
       a(6)  = -y(Si28) * rr.rates(index_rate, irr1) * rr.rates(1, irsigp);
       a(7)  = -y(Si28)*y(He4) * rr.rates(index_rate, irsiap)*(1.0e0_rt - rr.rates(1, irs1));
       a(8)  =  y(Si28)*y(He4) * rr.rates(1, irsiap)*rr.rates(index_rate, irs1);
       a(9)  = y(S32) * rr.rates(1, irs1) * rr.rates(index_rate, irsgp);
       a(10) = y(S32) * rr.rates(index_rate, irs1) * rr.rates(1, irsgp);

       dydt(Si28) = dydt(Si28) + esum10(a);

    }


    // s32 reactions
    a(1) =  0.1e0_rt * 0.5e0_rt*y(O16) * y(O16) * rr.rates(index_rate, ir1616);
    a(2) =  y(Si28) * y(He4) * rr.rates(index_rate, irsiag);
    a(3) = -y(S32) * y(He4) * rr.rates(index_rate, irsag);
    a(4) = -y(S32) * rr.rates(index_rate, irsga);
    a(5) =  y(Ar36) * rr.rates(index_rate, irarga);

    dydt(S32) = dydt(S32) + esum5(a);


    if (!deriva) {

       a(1) =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16)* rr.rates(index_rate, ir1616)*(1.0e0_rt-rr.rates(index_rate, irs1));
       a(2) =  y(Si28) * y(He4) * rr.rates(index_rate, irsiap)*(1.0e0_rt-rr.rates(index_rate, irs1));
       a(3) = -y(S32) * rr.rates(index_rate, irs1) * rr.rates(index_rate, irsgp);
       a(4) = -y(S32) * y(He4) * rr.rates(index_rate, irsap)*(1.0e0_rt-rr.rates(index_rate, irt1));
       a(5) =  y(Ar36) * rr.rates(index_rate, irt1) * rr.rates(index_rate, irargp);

       dydt(S32) = dydt(S32) + esum5(a);

    } else {

       a(1)  =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rr.rates(index_rate, ir1616)*(1.0e0_rt-rr.rates(1, irs1));
       a(2)  = -0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rr.rates(1, ir1616)*rr.rates(index_rate, irs1);
       a(3)  =  y(Si28)*y(He4) * rr.rates(index_rate, irsiap)*(1.0e0_rt-rr.rates(1, irs1));
       a(4)  = -y(Si28)*y(He4) * rr.rates(1, irsiap)*rr.rates(index_rate, irs1);
       a(5)  = -y(S32) * rr.rates(1, irs1) * rr.rates(index_rate, irsgp);
       a(6)  = -y(S32) * rr.rates(index_rate, irs1) * rr.rates(1, irsgp);
       a(7)  = -y(S32)*y(He4) * rr.rates(index_rate, irsap)*(1.0e0_rt-rr.rates(1, irt1));
       a(8)  =  y(S32)*y(He4) * rr.rates(1, irsap)*rr.rates(index_rate, irt1);
       a(9)  =  y(Ar36) * rr.rates(1, irt1) * rr.rates(index_rate, irargp);
       a(10) =  y(Ar36) * rr.rates(index_rate, irt1) * rr.rates(1, irargp);

       dydt(S32) = dydt(S32) + esum10(a);

    }


    // ar36 reactions
    a(1) =  y(S32)  * y(He4) * rr.rates(index_rate, irsag);
    a(2) = -y(Ar36) * y(He4) * rr.rates(index_rate, irarag);
    a(3) = -y(Ar36) * rr.rates(index_rate, irarga);
    a(4) =  y(Ca40) * rr.rates(index_rate, ircaga);

    dydt(Ar36) = dydt(Ar36) + esum4(a);

    if (!deriva) {

       a(1) = y(S32)  * y(He4) * rr.rates(index_rate, irsap)*(1.0e0_rt-rr.rates(index_rate, irt1));
       a(2) = -y(Ar36) * rr.rates(index_rate, irt1) * rr.rates(index_rate, irargp);
       a(3) = -y(Ar36) * y(He4) * rr.rates(index_rate, irarap)*(1.0e0_rt-rr.rates(index_rate, iru1));
       a(4) =  y(Ca40) * rr.rates(index_rate, ircagp) * rr.rates(index_rate, iru1);

       dydt(Ar36) = dydt(Ar36) + esum4(a);

    } else {

       a(1) =  y(S32)*y(He4) * rr.rates(index_rate, irsap)*(1.0e0_rt - rr.rates(1, irt1));
       a(2) = -y(S32)*y(He4) * rr.rates(1, irsap)*rr.rates(index_rate, irt1);
       a(3) = -y(Ar36) * rr.rates(1, irt1) * rr.rates(index_rate, irargp);
       a(4) = -y(Ar36) * rr.rates(index_rate, irt1) * rr.rates(1, irargp);
       a(5) = -y(Ar36)*y(He4) * rr.rates(index_rate, irarap)*(1.0e0_rt-rr.rates(1, iru1));
       a(6) =  y(Ar36)*y(He4) * rr.rates(1, irarap)*rr.rates(index_rate, iru1);
       a(7) =  y(Ca40) * rr.rates(1, ircagp) * rr.rates(index_rate, iru1);
       a(8) =  y(Ca40) * rr.rates(index_rate, ircagp) * rr.rates(1, iru1);

       dydt(Ar36) = dydt(Ar36) + esum8(a);

    }


    // ca40 reactions
    a(1) =  y(Ar36) * y(He4) * rr.rates(index_rate, irarag);
    a(2) = -y(Ca40) * y(He4) * rr.rates(index_rate, ircaag);
    a(3) = -y(Ca40) * rr.rates(index_rate, ircaga);
    a(4) =  y(Ti44) * rr.rates(index_rate, irtiga);

    dydt(Ca40) = dydt(Ca40) + esum4(a);

    if (!deriva) {

       a(1) =  y(Ar36) * y(He4) * rr.rates(index_rate, irarap)*(1.0e0_rt-rr.rates(index_rate, iru1));
       a(2) = -y(Ca40) * rr.rates(index_rate, ircagp) * rr.rates(index_rate, iru1);
       a(3) = -y(Ca40) * y(He4) * rr.rates(index_rate, ircaap)*(1.0e0_rt-rr.rates(index_rate, irv1));
       a(4) =  y(Ti44) * rr.rates(index_rate, irtigp) * rr.rates(index_rate, irv1);

       dydt(Ca40) = dydt(Ca40) + esum4(a);

    } else {

       a(1) =  y(Ar36)*y(He4) * rr.rates(index_rate, irarap)*(1.0e0_rt-rr.rates(1, iru1));
       a(2) = -y(Ar36)*y(He4) * rr.rates(1, irarap)*rr.rates(index_rate, iru1);
       a(3) = -y(Ca40) * rr.rates(1, ircagp) * rr.rates(index_rate, iru1);
       a(4) = -y(Ca40) * rr.rates(index_rate, ircagp) * rr.rates(1, iru1);
       a(5) = -y(Ca40)*y(He4) * rr.rates(index_rate, ircaap)*(1.0e0_rt-rr.rates(1, irv1));
       a(6) =  y(Ca40)*y(He4) * rr.rates(1, ircaap)*rr.rates(index_rate, irv1);
       a(7) =  y(Ti44) * rr.rates(1, irtigp) * rr.rates(index_rate, irv1);
       a(8) =  y(Ti44) * rr.rates(index_rate, irtigp) * rr.rates(1, irv1);

       dydt(Ca40) = dydt(Ca40) + esum8(a);

    }


    // ti44 reactions
    a(1) =  y(Ca40) * y(He4) * rr.rates(index_rate, ircaag);
    a(2) = -y(Ti44) * y(He4) * rr.rates(index_rate, irtiag);
    a(3) = -y(Ti44) * rr.rates(index_rate, irtiga);
    a(4) =  y(Cr48) * rr.rates(index_rate, ircrga);

    dydt(Ti44) = dydt(Ti44) + esum4(a);

    if (!deriva) {

       a(1) =  y(Ca40) * y(He4) * rr.rates(index_rate, ircaap)*(1.0e0_rt-rr.rates(index_rate, irv1));
       a(2) = -y(Ti44) * rr.rates(index_rate, irv1) * rr.rates(index_rate, irtigp);
       a(3) = -y(Ti44) * y(He4) * rr.rates(index_rate, irtiap)*(1.0e0_rt-rr.rates(index_rate, irw1));
       a(4) =  y(Cr48) * rr.rates(index_rate, irw1) * rr.rates(index_rate, ircrgp);

       dydt(Ti44) = dydt(Ti44) + esum4(a);

    } else {

       a(1) =  y(Ca40)*y(He4) * rr.rates(index_rate, ircaap)*(1.0e0_rt-rr.rates(1, irv1));
       a(2) = -y(Ca40)*y(He4) * rr.rates(1, ircaap)*rr.rates(index_rate, irv1);
       a(3) = -y(Ti44) * rr.rates(1, irv1) * rr.rates(index_rate, irtigp);
       a(4) = -y(Ti44) * rr.rates(index_rate, irv1) * rr.rates(1, irtigp);
       a(5) = -y(Ti44)*y(He4) * rr.rates(index_rate, irtiap)*(1.0e0_rt-rr.rates(1, irw1));
       a(6) =  y(Ti44)*y(He4) * rr.rates(1, irtiap)*rr.rates(index_rate, irw1);
       a(7) =  y(Cr48) * rr.rates(1, irw1) * rr.rates(index_rate, ircrgp);
       a(8) =  y(Cr48) * rr.rates(index_rate, irw1) * rr.rates(1, ircrgp);

       dydt(Ti44) = dydt(Ti44) + esum8(a);

    }


    // cr48 reactions
    a(1) =  y(Ti44) * y(He4) * rr.rates(index_rate, irtiag);
    a(2) = -y(Cr48) * y(He4) * rr.rates(index_rate, ircrag);
    a(3) = -y(Cr48) * rr.rates(index_rate, ircrga);
    a(4) =  y(Fe52) * rr.rates(index_rate, irfega);

    dydt(Cr48) = dydt(Cr48) + esum4(a);

    if (!deriva) {

       a(1) =  y(Ti44) * y(He4) * rr.rates(index_rate, irtiap)*(1.0e0_rt-rr.rates(index_rate, irw1));
       a(2) = -y(Cr48) * rr.rates(index_rate, irw1) * rr.rates(index_rate, ircrgp);
       a(3) = -y(Cr48) * y(He4) * rr.rates(index_rate, ircrap)*(1.0e0_rt-rr.rates(index_rate, irx1));
       a(4) =  y(Fe52) * rr.rates(index_rate, irx1) * rr.rates(index_rate, irfegp);

       dydt(Cr48) = dydt(Cr48) + esum4(a);

    } else {

       a(1) =  y(Ti44)*y(He4) * rr.rates(index_rate, irtiap)*(1.0e0_rt-rr.rates(1, irw1));
       a(2) = -y(Ti44)*y(He4) * rr.rates(1, irtiap)*rr.rates(index_rate, irw1);
       a(3) = -y(Cr48) * rr.rates(1, irw1) * rr.rates(index_rate, ircrgp);
       a(4) = -y(Cr48) * rr.rates(index_rate, irw1) * rr.rates(1, ircrgp);
       a(5) = -y(Cr48)*y(He4) * rr.rates(index_rate, ircrap)*(1.0e0_rt-rr.rates(1, irx1));
       a(6) =  y(Cr48)*y(He4) * rr.rates(1, ircrap)*rr.rates(index_rate, irx1);
       a(7) =  y(Fe52) * rr.rates(1, irx1) * rr.rates(index_rate, irfegp);
       a(8) =  y(Fe52) * rr.rates(index_rate, irx1) * rr.rates(1, irfegp);

       dydt(Cr48) = dydt(Cr48) + esum8(a);

    }


    // fe52 reactions
    a(1) =  y(Cr48) * y(He4) * rr.rates(index_rate, ircrag);
    a(2) = -y(Fe52) * y(He4) * rr.rates(index_rate, irfeag);
    a(3) = -y(Fe52) * rr.rates(index_rate, irfega);
    a(4) =  y(Ni56) * rr.rates(index_rate, irniga);

    dydt(Fe52) = dydt(Fe52) + esum4(a);

    if (!deriva) {

       a(1) =  y(Cr48) * y(He4) * rr.rates(index_rate, ircrap)*(1.0e0_rt-rr.rates(index_rate, irx1));
       a(2) = -y(Fe52) * rr.rates(index_rate, irx1) * rr.rates(index_rate, irfegp);
       a(3) = -y(Fe52) * y(He4) * rr.rates(index_rate, irfeap)*(1.0e0_rt-rr.rates(index_rate, iry1));
       a(4) =  y(Ni56) * rr.rates(index_rate, iry1) * rr.rates(index_rate, irnigp);

       dydt(Fe52) = dydt(Fe52) + esum4(a);

    } else {

       a(1) =  y(Cr48)*y(He4) * rr.rates(index_rate, ircrap)*(1.0e0_rt-rr.rates(1, irx1));
       a(2) = -y(Cr48)*y(He4) * rr.rates(1, ircrap)*rr.rates(index_rate, irx1);
       a(3) = -y(Fe52) * rr.rates(1, irx1) * rr.rates(index_rate, irfegp);
       a(4) = -y(Fe52) * rr.rates(index_rate, irx1) * rr.rates(1, irfegp);
       a(5) = -y(Fe52)*y(He4) * rr.rates(index_rate, irfeap)*(1.0e0_rt-rr.rates(1, iry1));
       a(6) =  y(Fe52)*y(He4) * rr.rates(1, irfeap)*rr.rates(index_rate, iry1);
       a(7) =  y(Ni56) * rr.rates(1, iry1) * rr.rates(index_rate, irnigp);
       a(8) =  y(Ni56) * rr.rates(index_rate, iry1) * rr.rates(1, irnigp);

       dydt(Fe52) = dydt(Fe52) + esum8(a);

    }


    // ni56 reactions
    a(1) =  y(Fe52) * y(He4) * rr.rates(index_rate, irfeag);
    a(2) = -y(Ni56) * rr.rates(index_rate, irniga);

    dydt(Ni56) = dydt(Ni56) + a(1) + a(2);

    if (!deriva) {

       a(1) =  y(Fe52) * y(He4) * rr.rates(index_rate, irfeap)*(1.0e0_rt-rr.rates(index_rate, iry1));
       a(2) = -y(Ni56) * rr.rates(index_rate, iry1) * rr.rates(index_rate, irnigp);

       dydt(Ni56) = dydt(Ni56) + a(1) + a(2);

    } else {

       a(1) =  y(Fe52)*y(He4) * rr.rates(index_rate, irfeap)*(1.0e0_rt-rr.rates(1, iry1));
       a(2) = -y(Fe52)*y(He4) * rr.rates(1, irfeap)*rr.rates(index_rate, iry1);
       a(3) = -y(Ni56) * rr.rates(1, iry1) * rr.rates(index_rate, irnigp);
       a(4) = -y(Ni56) * rr.rates(index_rate, iry1) * rr.rates(1, irnigp);

       dydt(Ni56) = dydt(Ni56) + esum4(a);

    }

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_rhs(burn_t& state, Array1D<Real, 1, neqs>& ydot)
{

    /*
     This routine sets up the system of ODE's for the aprox13
     nuclear reactions.  This is an alpha chain + heavy ion network
     with (a,p)(p,g) links.

     Isotopes: he4,  c12,  o16,  ne20, mg24, si28, s32,
               ar36, ca40, ti44, cr48, fe52, ni56
    */

    rate_t rr;

    Real sneut, dsneutdt, dsneutdd, snuda, snudz;
    Real enuc;

    Real rho, temp, abar, zbar;
    Array1D<Real, 1, NumSpec> y;

    // Initialize ydot to 0

    for (int i = 1; i <= neqs; ++i) {
        ydot(i) = 0.0_rt;
    }

    // Evaluate the rates

    evaluate_rates(state, rr);

    // Get the data from the state

    rho  = state.rho;
    temp = state.T;
    abar = state.abar;
    zbar = state.zbar;

    for (int i = 1; i <= NumSpec; ++i) {
        y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    // Call the RHS to actually get dydt.

    bool deriva = false;
    bool for_jacobian_tderiv = false;
    rhs(y, rr, ydot, deriva, for_jacobian_tderiv);

    // Instantaneous energy generation rate -- this needs molar fractions

    ener_gener_rate(ydot, enuc);

    // Get the neutrino losses

    sneut5(temp, rho, abar, zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);

    // Append the energy equation (this is erg/g/s)

    ydot(net_ienuc) = enuc - sneut;

    // Append the temperature equation

    temperature_rhs(state, ydot);
}


// Analytical Jacobian
template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_jac(burn_t& state, MatrixType& jac)
{

    rate_t rr;

    bool deriva;

    Real b1, sneut, dsneutdt, dsneutdd, snuda, snudz;

    Real rho, temp, abar, zbar;
    Array1D<Real, 1, NumSpec> y;
    Array1D<Real, 1, neqs> yderivs;

    // Initialize jac to 0

    jac.zero();

    // Evaluate the rates

    evaluate_rates(state, rr);

    // Get the data from the state

    rho  = state.rho;
    temp = state.T;
    abar = state.abar;
    zbar = state.zbar;

    for (int i = 1; i <= NumSpec; ++i) {
        y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    // Species Jacobian elements with respect to other species

    dfdy_isotopes_aprox13(y, state, rr, jac);

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

    // Evaluate the Jacobian elements with respect to temperature by
    // calling the RHS using d(rate) / dT

    deriva = true;
    bool for_jacobian_tderiv = true;
    rhs(y, rr, yderivs, deriva, for_jacobian_tderiv);

    for (int i = 1; i <= NumSpec; ++i) {
        jac(i,net_itemp) = yderivs(i);
    }

    ener_gener_rate(yderivs, jac(net_ienuc,net_itemp));

    jac(net_ienuc,net_itemp) -= dsneutdt;

    // Temperature Jacobian elements

    temperature_jac(state, jac);

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

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], 4.0e0_rt, 8.0e0_rt);

    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[C12-1], aion[C12-1]);

    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[O16-1], aion[O16-1]);

    add_screening_factor(jscr++, zion[O16-1], aion[O16-1], zion[O16-1], aion[O16-1]);

    add_screening_factor(jscr++, zion[O16-1], aion[O16-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, zion[Ne20-1], aion[Ne20-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, zion[Mg24-1], aion[Mg24-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, 13.0e0_rt, 27.0e0_rt, 1.0e0_rt, 1.0e0_rt);

    add_screening_factor(jscr++, zion[Si28-1], aion[Si28-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, 15.0e0_rt, 31.0e0_rt, 1.0e0_rt, 1.0e0_rt);

    add_screening_factor(jscr++, zion[S32-1], aion[S32-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, 17.0e0_rt, 35.0e0_rt, 1.0e0_rt, 1.0e0_rt);

    add_screening_factor(jscr++, zion[Ar36-1], aion[Ar36-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, 19.0e0_rt, 39.0e0_rt, 1.0e0_rt, 1.0e0_rt);

    add_screening_factor(jscr++, zion[Ca40-1], aion[Ca40-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, 21.0e0_rt, 43.0e0_rt, 1.0e0_rt, 1.0e0_rt);

    add_screening_factor(jscr++, zion[Ti44-1], aion[Ti44-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, 23.0e0_rt, 47.0e0_rt, 1.0e0_rt, 1.0e0_rt);

    add_screening_factor(jscr++, zion[Cr48-1], aion[Cr48-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, 25.0e0_rt, 51.0e0_rt, 1.0e0_rt, 1.0e0_rt);

    add_screening_factor(jscr++, zion[Fe52-1], aion[Fe52-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, 27.0e0_rt, 55.0e0_rt, 1.0e0_rt, 1.0e0_rt);

}

#endif
