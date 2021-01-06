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
void aprox13tab(const Real btemp, const Real bden, Array1D<rate_t, 1, Rates::NumGroups>& rr)
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

       rr(1).rates(j) = (  alfa * rattab(j,iat  )
                        + beta * rattab(j,iat+1)
                        + gama * rattab(j,iat+2)
                        + delt * rattab(j,iat+3) ) * dtab(j);

       rr(2).rates(j) = (  alfa * drattabdt(j,iat  )
                        + beta * drattabdt(j,iat+1)
                        + gama * drattabdt(j,iat+2)
                        + delt * drattabdt(j,iat+3) ) * dtab(j);

    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void aprox13rat(const Real btemp, const Real bden, Array1D<rate_t, 1, Rates::NumGroups>& rr)
{
    using namespace Rates;

    // this routine generates unscreened
    // nuclear reaction rates for the aprox13 network.

    Real rrate,drratedt;

    for (int i = 1; i <= Rates::NumRates; ++i) {
       rr(1).rates(i) = 0.0_rt;
       rr(2).rates(i) = 0.0_rt;
    }

    if (btemp < 1.0e6_rt) return;


    // get the temperature factors
    tf_t tf = get_tfactors(btemp);


    // Determine which c12(a,g)o16 rate to use
    if (use_c12ag_deboer17) {
        // deboer + 2017 c12(a,g)o16 rate
        rate_c12ag_deboer17(tf,bden,
                            rr(1).rates(ircag),rr(2).rates(ircag),
                            rr(1).rates(iroga),rr(2).rates(iroga));
    } else {
        // 1.7 times cf88 c12(a,g)o16 rate
        rate_c12ag(tf,bden,
                   rr(1).rates(ircag),rr(2).rates(ircag),
                   rr(1).rates(iroga),rr(2).rates(iroga));
    }

    // triple alpha to c12
    rate_triplealf(tf,bden,
                   rr(1).rates(ir3a),rr(2).rates(ir3a),
                   rr(1).rates(irg3a),rr(2).rates(irg3a));

    // c12 + c12
    rate_c12c12(tf,bden,
                rr(1).rates(ir1212),rr(2).rates(ir1212),
                rrate,drratedt);

    // c12 + o16
    rate_c12o16(tf,bden,
                rr(1).rates(ir1216),rr(2).rates(ir1216),
                rrate,drratedt);

    // o16 + o16
    rate_o16o16(tf,bden,
                rr(1).rates(ir1616),rr(2).rates(ir1616),
                rrate,drratedt);

    // o16(a,g)ne20
    rate_o16ag(tf,bden,
               rr(1).rates(iroag),rr(2).rates(iroag),
               rr(1).rates(irnega),rr(2).rates(irnega));

    // ne20(a,g)mg24
    rate_ne20ag(tf,bden,
                rr(1).rates(irneag),rr(2).rates(irneag),
                rr(1).rates(irmgga),rr(2).rates(irmgga));

    // mg24(a,g)si28
    rate_mg24ag(tf,bden,
                rr(1).rates(irmgag),rr(2).rates(irmgag),
                rr(1).rates(irsiga),rr(2).rates(irsiga));

    // mg24(a,p)al27
    rate_mg24ap(tf,bden,
                rr(1).rates(irmgap),rr(2).rates(irmgap),
                rr(1).rates(iralpa),rr(2).rates(iralpa));

    // al27(p,g)si28
    rate_al27pg(tf,bden,
                rr(1).rates(iralpg),rr(2).rates(iralpg),
                rr(1).rates(irsigp),rr(2).rates(irsigp));

    // si28(a,g)s32
    rate_si28ag(tf,bden,
                rr(1).rates(irsiag),rr(2).rates(irsiag),
                rr(1).rates(irsga),rr(2).rates(irsga));

    // si28(a,p)p31
    rate_si28ap(tf,bden,
                rr(1).rates(irsiap),rr(2).rates(irsiap),
                rr(1).rates(irppa),rr(2).rates(irppa));

    // p31(p,g)s32
    rate_p31pg(tf,bden,
               rr(1).rates(irppg),rr(2).rates(irppg),
               rr(1).rates(irsgp),rr(2).rates(irsgp));

    // s32(a,g)ar36
    rate_s32ag(tf,bden,
               rr(1).rates(irsag),rr(2).rates(irsag),
               rr(1).rates(irarga),rr(2).rates(irarga));

    // s32(a,p)cl35
    rate_s32ap(tf,bden,
               rr(1).rates(irsap),rr(2).rates(irsap),
               rr(1).rates(irclpa),rr(2).rates(irclpa));

    // cl35(p,g)ar36
    rate_cl35pg(tf,bden,
                rr(1).rates(irclpg),rr(2).rates(irclpg),
                rr(1).rates(irargp),rr(2).rates(irargp));

    // ar36(a,g)ca40
    rate_ar36ag(tf,bden,
                rr(1).rates(irarag),rr(2).rates(irarag),
                rr(1).rates(ircaga),rr(2).rates(ircaga));

    // ar36(a,p)k39
    rate_ar36ap(tf,bden,
                rr(1).rates(irarap),rr(2).rates(irarap),
                rr(1).rates(irkpa),rr(2).rates(irkpa));

    // k39(p,g)ca40
    rate_k39pg(tf,bden,
               rr(1).rates(irkpg),rr(2).rates(irkpg),
               rr(1).rates(ircagp),rr(2).rates(ircagp));

    // ca40(a,g)ti44
    rate_ca40ag(tf,bden,
                rr(1).rates(ircaag),rr(2).rates(ircaag),
                rr(1).rates(irtiga),rr(2).rates(irtiga));

    // ca40(a,p)sc43
    rate_ca40ap(tf,bden,
                rr(1).rates(ircaap),rr(2).rates(ircaap),
                rr(1).rates(irscpa),rr(2).rates(irscpa));

    // sc43(p,g)ti44
    rate_sc43pg(tf,bden,
                rr(1).rates(irscpg),rr(2).rates(irscpg),
                rr(1).rates(irtigp),rr(2).rates(irtigp));

    // ti44(a,g)cr48
    rate_ti44ag(tf,bden,
                rr(1).rates(irtiag),rr(2).rates(irtiag),
                rr(1).rates(ircrga),rr(2).rates(ircrga));

    // ti44(a,p)v47
    rate_ti44ap(tf,bden,
                rr(1).rates(irtiap),rr(2).rates(irtiap), 
                rr(1).rates(irvpa),rr(2).rates(irvpa));

    // v47(p,g)cr48
    rate_v47pg(tf,bden,
               rr(1).rates(irvpg),rr(2).rates(irvpg), 
               rr(1).rates(ircrgp),rr(2).rates(ircrgp));

    // cr48(a,g)fe52
    rate_cr48ag(tf,bden,
                rr(1).rates(ircrag),rr(2).rates(ircrag),
                rr(1).rates(irfega),rr(2).rates(irfega));

    // cr48(a,p)mn51
    rate_cr48ap(tf,bden,
                rr(1).rates(ircrap),rr(2).rates(ircrap),
                rr(1).rates(irmnpa),rr(2).rates(irmnpa));

    // mn51(p,g)fe52
    rate_mn51pg(tf,bden,
                rr(1).rates(irmnpg),rr(2).rates(irmnpg),
                rr(1).rates(irfegp),rr(2).rates(irfegp));

    // fe52(a,g)ni56
    rate_fe52ag(tf,bden,
                rr(1).rates(irfeag),rr(2).rates(irfeag),
                rr(1).rates(irniga),rr(2).rates(irniga));

    // fe52(a,p)co55
    rate_fe52ap(tf,bden,
                rr(1).rates(irfeap),rr(2).rates(irfeap),
                rr(1).rates(ircopa),rr(2).rates(ircopa));

    // co55(p,g)ni56
    rate_co55pg(tf,bden,
                rr(1).rates(ircopg),rr(2).rates(ircopg),
                rr(1).rates(irnigp),rr(2).rates(irnigp));
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void set_aprox13rat()
{
    using namespace RateTable;

    Real btemp;
    Real bden = 1.0e0_rt;
    Array1D<rate_t, 1, Rates::NumGroups> rr;

    for (int i = 1; i <= tab_imax; ++i) {

       btemp = tab_tlo + static_cast<Real>(i-1) * tab_tstp;
       btemp = std::pow(10.0e0_rt, btemp);

       aprox13rat(btemp, bden, rr);

       ttab(i) = btemp;

       for (int j = 1; j <= Rates::NumRates; ++j) {

          rattab(j,i)    = rr(1).rates(j);
          drattabdt(j,i) = rr(2).rates(j);

       }
    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void screen_aprox13(const Real btemp, const Real bden,
                    Array1D<Real, 1, NumSpec> const& y,
                    Array1D<rate_t, 1, Rates::NumGroups>& rr)
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

    ratraw = rr(1).rates(ir3a);
    rr(1).rates(ir3a) = ratraw * sc3a;
    rr(2).rates(ir3a) = rr(2).rates(ir3a)*sc3a + ratraw*sc3adt;

    ratraw = rr(1).rates(irg3a);
    rr(1).rates(irg3a) = ratraw * sc3a;
    rr(2).rates(irg3a) = rr(2).rates(irg3a)*sc3a + ratraw*sc3adt;

    // c12 to o16
    // c12(a,g)o16
    jscr++;
    screen5(state,jscr,
            zion[C12-1], aion[C12-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);

    ratraw = rr(1).rates(ircag);
    rr(1).rates(ircag)  = ratraw * sc1a;
    rr(2).rates(ircag)  = rr(2).rates(ircag)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(iroga);
    rr(1).rates(iroga)  = ratraw * sc1a;
    rr(2).rates(iroga)  = rr(2).rates(iroga)*sc1a + ratraw*sc1adt;

    // c12 + c12
    jscr++;
    screen5(state,jscr,
            zion[C12-1], aion[C12-1], zion[C12-1], aion[C12-1],
            sc1a,sc1adt,sc1add);

    ratraw = rr(1).rates(ir1212);
    rr(1).rates(ir1212) = ratraw * sc1a;
    rr(2).rates(ir1212) = rr(2).rates(ir1212)*sc1a + ratraw*sc1adt;

    // c12 + o16
    jscr++;
    screen5(state,jscr,
            zion[C12-1], aion[C12-1], zion[O16-1], aion[O16-1], 
            sc1a,sc1adt,sc1add);

    ratraw = rr(1).rates(ir1216);
    rr(1).rates(ir1216) = ratraw * sc1a;
    rr(2).rates(ir1216) = rr(2).rates(ir1216)*sc1a + ratraw*sc1adt;

    // o16 + o16
    jscr++;
    screen5(state,jscr,
            zion[O16-1], aion[O16-1], zion[O16-1], aion[O16-1],
            sc1a,sc1adt,sc1add);

    ratraw = rr(1).rates(ir1616);
    rr(1).rates(ir1616) = ratraw * sc1a;
    rr(2).rates(ir1616) = rr(2).rates(ir1616)*sc1a + ratraw*sc1adt;

    // o16 to ne20
    jscr++;
    screen5(state,jscr,
            zion[O16-1], aion[O16-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);


    // o16(a,g)ne20
    ratraw = rr(1).rates(iroag);
    rr(1).rates(iroag) = ratraw * sc1a;
    rr(2).rates(iroag) = rr(2).rates(iroag)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irnega);
    rr(1).rates(irnega) = ratraw * sc1a;
    rr(2).rates(irnega) = rr(2).rates(irnega)*sc1a + ratraw*sc1adt;

    // ne20 to mg24
    jscr++;
    screen5(state,jscr,
            zion[Ne20-1], aion[Ne20-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);


    // ne20(a,g)mg24
    ratraw = rr(1).rates(irneag);
    rr(1).rates(irneag) = ratraw * sc1a;
    rr(2).rates(irneag) = rr(2).rates(irneag)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irmgga);
    rr(1).rates(irmgga) = ratraw * sc1a;
    rr(2).rates(irmgga) = rr(2).rates(irmgga)*sc1a + ratraw*sc1adt;

    // mg24 to si28
    jscr++;
    screen5(state,jscr,
            zion[Mg24-1], aion[Mg24-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);


    // mg24(a,g)si28
    ratraw = rr(1).rates(irmgag);
    rr(1).rates(irmgag) = ratraw * sc1a;
    rr(2).rates(irmgag) = rr(2).rates(irmgag)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irsiga);
    rr(1).rates(irsiga) = ratraw * sc1a;
    rr(2).rates(irsiga) = rr(2).rates(irsiga)*sc1a + ratraw*sc1adt;


    // mg24(a,p)al27
    ratraw = rr(1).rates(irmgap);
    rr(1).rates(irmgap) = ratraw * sc1a;
    rr(2).rates(irmgap) = rr(2).rates(irmgap)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(iralpa);
    rr(1).rates(iralpa) = ratraw * sc1a;
    rr(2).rates(iralpa) = rr(2).rates(iralpa)*sc1a + ratraw*sc1adt;

    jscr++;
    screen5(state,jscr,
            13.0_rt, 27.0_rt, 1.0_rt, 1.0_rt,
            sc1a,sc1adt,sc1add);

    // al27(p,g)si28
    ratraw = rr(1).rates(iralpg);
    rr(1).rates(iralpg) = ratraw * sc1a;
    rr(2).rates(iralpg) = rr(2).rates(iralpg)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irsigp);
    rr(1).rates(irsigp) = ratraw * sc1a;
    rr(2).rates(irsigp) = rr(2).rates(irsigp)*sc1a + ratraw*sc1adt;


    // si28 to s32
    jscr++;
    screen5(state,jscr,
            zion[Si28-1], aion[Si28-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);


    // si28(a,g)s32
    ratraw = rr(1).rates(irsiag);
    rr(1).rates(irsiag) = ratraw * sc1a;
    rr(2).rates(irsiag) = rr(2).rates(irsiag)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irsga);
    rr(1).rates(irsga) = ratraw * sc1a;
    rr(2).rates(irsga) = rr(2).rates(irsga)*sc1a + ratraw*sc1adt;


    // si28(a,p)p31
    ratraw = rr(1).rates(irsiap);
    rr(1).rates(irsiap) = ratraw * sc1a;
    rr(2).rates(irsiap) = rr(2).rates(irsiap)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irppa);
    rr(1).rates(irppa)  = ratraw * sc1a;
    rr(2).rates(irppa)  = rr(2).rates(irppa)*sc1a  + ratraw*sc1adt;


    jscr++;
    screen5(state,jscr,
            15.0_rt, 31.0_rt, 1.0_rt, 1.0_rt,
            sc1a,sc1adt,sc1add);

    // p31(p,g)s32
    ratraw = rr(1).rates(irppg);
    rr(1).rates(irppg)  = ratraw * sc1a;
    rr(2).rates(irppg)  = rr(2).rates(irppg)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irsgp);
    rr(1).rates(irsgp)  = ratraw * sc1a;
    rr(2).rates(irsgp)  = rr(2).rates(irsgp)*sc1a + ratraw*sc1adt;


    // s32 to ar36
    jscr++;
    screen5(state,jscr,
            zion[S32-1], aion[S32-1], zion[He4-1], aion[He4-1], 
            sc1a,sc1adt,sc1add);


    // s32(a,g)ar36
    ratraw = rr(1).rates(irsag);
    rr(1).rates(irsag)  = ratraw * sc1a;
    rr(2).rates(irsag)  = rr(2).rates(irsag)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irarga);
    rr(1).rates(irarga)  = ratraw * sc1a;
    rr(2).rates(irarga)  = rr(2).rates(irarga)*sc1a + ratraw*sc1adt;

    // s32(a,p)cl35
    ratraw = rr(1).rates(irsap);
    rr(1).rates(irsap)  = ratraw * sc1a;
    rr(2).rates(irsap)  = rr(2).rates(irsap)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irclpa);
    rr(1).rates(irclpa) = ratraw * sc1a;
    rr(2).rates(irclpa) = rr(2).rates(irclpa)*sc1a + ratraw*sc1adt;


    jscr++;
    screen5(state,jscr,
            17.0_rt, 35.0_rt, 1.0_rt, 1.0_rt,
            sc1a,sc1adt,sc1add);

    // cl35(p,g)ar36
    ratraw = rr(1).rates(irclpg);
    rr(1).rates(irclpg) = ratraw * sc1a;
    rr(2).rates(irclpg) = rr(2).rates(irclpg)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irargp);
    rr(1).rates(irargp) = ratraw * sc1a;
    rr(2).rates(irargp) = rr(2).rates(irargp)*sc1a + ratraw*sc1adt;


    // ar36 to ca40
    jscr++;
    screen5(state,jscr,
            zion[Ar36-1], aion[Ar36-1], zion[He4-1], aion[He4-1], 
            sc1a,sc1adt,sc1add);


    // ar36(a,g)ca40
    ratraw = rr(1).rates(irarag);
    rr(1).rates(irarag) = ratraw * sc1a;
    rr(2).rates(irarag) = rr(2).rates(irarag)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(ircaga);
    rr(1).rates(ircaga) = ratraw * sc1a;
    rr(2).rates(ircaga) = rr(2).rates(ircaga)*sc1a + ratraw*sc1adt;


    // ar36(a,p)k39
    ratraw = rr(1).rates(irarap);
    rr(1).rates(irarap) = ratraw * sc1a;
    rr(2).rates(irarap) = rr(2).rates(irarap)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irkpa);
    rr(1).rates(irkpa) = ratraw * sc1a;
    rr(2).rates(irkpa) = rr(2).rates(irkpa)*sc1a  + ratraw*sc1adt;


    jscr++;
    screen5(state,jscr,
            19.0_rt, 39.0_rt, 1.0_rt, 1.0_rt,
            sc1a,sc1adt,sc1add);

    // k39(p,g)ca40
    ratraw = rr(1).rates(irkpg);
    rr(1).rates(irkpg) = ratraw * sc1a;
    rr(2).rates(irkpg) = rr(2).rates(irkpg)*sc1a  + ratraw*sc1adt;

    ratraw = rr(1).rates(ircagp);
    rr(1).rates(ircagp) = ratraw * sc1a;
    rr(2).rates(ircagp) = rr(2).rates(ircagp)*sc1a  + ratraw*sc1adt;


    // ca40 to ti44
    jscr++;
    screen5(state,jscr,
            zion[Ca40-1], aion[Ca40-1], zion[He4-1], aion[He4-1], 
            sc1a,sc1adt,sc1add);


    // ca40(a,g)ti44
    ratraw = rr(1).rates(ircaag);
    rr(1).rates(ircaag) = ratraw * sc1a;
    rr(2).rates(ircaag) = rr(2).rates(ircaag)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irtiga);
    rr(1).rates(irtiga) = ratraw * sc1a;
    rr(2).rates(irtiga) = rr(2).rates(irtiga)*sc1a + ratraw*sc1adt;


    // ca40(a,p)sc43
    ratraw = rr(1).rates(ircaap);
    rr(1).rates(ircaap) = ratraw * sc1a;
    rr(2).rates(ircaap) = rr(2).rates(ircaap)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irscpa);
    rr(1).rates(irscpa) = ratraw * sc1a;
    rr(2).rates(irscpa) = rr(2).rates(irscpa)*sc1a + ratraw*sc1adt;


    jscr++;
    screen5(state,jscr,
            21.0_rt, 43.0_rt, 1.0_rt, 1.0_rt,
            sc1a,sc1adt,sc1add);

    // sc43(p,g)ti44
    ratraw = rr(1).rates(irscpg);
    rr(1).rates(irscpg) = ratraw * sc1a;
    rr(2).rates(irscpg) = rr(2).rates(irscpg)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irtigp);
    rr(1).rates(irtigp) = ratraw * sc1a;
    rr(2).rates(irtigp) = rr(2).rates(irtigp)*sc1a + ratraw*sc1adt;


    // ti44 to cr48
    jscr++;
    screen5(state,jscr,
            zion[Ti44-1], aion[Ti44-1], zion[He4-1], aion[He4-1], 
            sc1a,sc1adt,sc1add);


    // ti44(a,g)cr48
    ratraw = rr(1).rates(irtiag);
    rr(1).rates(irtiag) = ratraw * sc1a;
    rr(2).rates(irtiag) = rr(2).rates(irtiag)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(ircrga);
    rr(1).rates(ircrga) = ratraw * sc1a;
    rr(2).rates(ircrga) = rr(2).rates(ircrga)*sc1a + ratraw*sc1adt;

    // ti44(a,p)v47
    ratraw = rr(1).rates(irtiap);
    rr(1).rates(irtiap) = ratraw * sc1a;
    rr(2).rates(irtiap) = rr(2).rates(irtiap)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irvpa);
    rr(1).rates(irvpa) = ratraw * sc1a;
    rr(2).rates(irvpa) = rr(2).rates(irvpa)*sc1a  + ratraw*sc1adt;


    jscr++;
    screen5(state,jscr,
            23.0_rt, 47.0_rt, 1.0_rt, 1.0_rt,
            sc1a,sc1adt,sc1add);

    // v47(p,g)cr48
    ratraw = rr(1).rates(irvpg);
    rr(1).rates(irvpg) = ratraw * sc1a;
    rr(2).rates(irvpg) = rr(2).rates(irvpg)*sc1a  + ratraw*sc1adt;

    ratraw = rr(1).rates(ircrgp);
    rr(1).rates(ircrgp)  = ratraw * sc1a;
    rr(2).rates(ircrgp)  = rr(2).rates(ircrgp)*sc1a  + ratraw*sc1adt;


    // cr48 to fe52
    jscr++;
    screen5(state,jscr,
            zion[Cr48-1], aion[Cr48-1], zion[He4-1], aion[He4-1], 
            sc1a,sc1adt,sc1add);


    // cr48(a,g)fe52
    ratraw = rr(1).rates(ircrag);
    rr(1).rates(ircrag) = ratraw * sc1a;
    rr(2).rates(ircrag) = rr(2).rates(ircrag)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irfega);
    rr(1).rates(irfega) = ratraw * sc1a;
    rr(2).rates(irfega) = rr(2).rates(irfega)*sc1a + ratraw*sc1adt;


    // cr48(a,p)mn51
    ratraw = rr(1).rates(ircrap);
    rr(1).rates(ircrap) = ratraw * sc1a;
    rr(2).rates(ircrap) = rr(2).rates(ircrap)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irmnpa);
    rr(1).rates(irmnpa) = ratraw * sc1a;
    rr(2).rates(irmnpa) = rr(2).rates(irmnpa)*sc1a + ratraw*sc1adt;


    jscr++;
    screen5(state,jscr,
            25.0_rt, 51.0_rt, 1.0_rt, 1.0_rt,
            sc1a,sc1adt,sc1add);

    // mn51(p,g)fe52
    ratraw = rr(1).rates(irmnpg);
    rr(1).rates(irmnpg) = ratraw * sc1a;
    rr(2).rates(irmnpg) = rr(2).rates(irmnpg)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irfegp);
    rr(1).rates(irfegp) = ratraw * sc1a;
    rr(2).rates(irfegp) = rr(2).rates(irfegp)*sc1a + ratraw*sc1adt;


    // fe52 to ni56
    jscr++;
    screen5(state,jscr,
            zion[Fe52-1], aion[Fe52-1], zion[He4-1], aion[He4-1], 
            sc1a,sc1adt,sc1add);


    // fe52(a,g)ni56
    ratraw = rr(1).rates(irfeag);
    rr(1).rates(irfeag) = ratraw * sc1a;
    rr(2).rates(irfeag) = rr(2).rates(irfeag)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irniga);
    rr(1).rates(irniga) = ratraw * sc1a;
    rr(2).rates(irniga) = rr(2).rates(irniga)*sc1a + ratraw*sc1adt;


    // fe52(a,p)co55
    ratraw = rr(1).rates(irfeap);
    rr(1).rates(irfeap) = ratraw * sc1a;
    rr(2).rates(irfeap) = rr(2).rates(irfeap)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(ircopa);
    rr(1).rates(ircopa) = ratraw * sc1a;
    rr(2).rates(ircopa) = rr(2).rates(ircopa)*sc1a + ratraw*sc1adt;


    jscr++;
    screen5(state,jscr,
            27.0_rt, 55.0_rt, 1.0_rt, 1.0_rt,
            sc1a,sc1adt,sc1add);


    // co55(p,g)ni56
    ratraw = rr(1).rates(ircopg);
    rr(1).rates(ircopg) = ratraw * sc1a;
    rr(2).rates(ircopg) = rr(2).rates(ircopg)*sc1a + ratraw*sc1adt;

    ratraw = rr(1).rates(irnigp);
    rr(1).rates(irnigp) = ratraw * sc1a;
    rr(2).rates(irnigp) = rr(2).rates(irnigp)*sc1a + ratraw*sc1adt;


    // now form those lovely dummy proton link rates

    // mg24(a,p)27al(p,g)28si
    rr(1).rates(irr1)  = 0.0e0_rt;
    rr(2).rates(irr1)  = 0.0e0_rt;
    denom    = rr(1).rates(iralpa) + rr(1).rates(iralpg);
    denomdt  = rr(2).rates(iralpa) + rr(2).rates(iralpg);

    if (denom > 1.0e-30_rt) {
       zz = 1.0e0_rt/denom;
       rr(1).rates(irr1) = rr(1).rates(iralpa)*zz;
       rr(2).rates(irr1) = (rr(2).rates(iralpa) - rr(1).rates(irr1)*denomdt)*zz;
    }

    // si28(a,p)p31(p,g)s32
    rr(1).rates(irs1)  = 0.0e0_rt;
    rr(2).rates(irs1)  = 0.0e0_rt;
    denom    = rr(1).rates(irppa) + rr(1).rates(irppg);
    denomdt  = rr(2).rates(irppa) + rr(2).rates(irppg);
    if (denom > 1.0e-30_rt) {
       zz = 1.0e0_rt/denom;
       rr(1).rates(irs1) = rr(1).rates(irppa)*zz;
       rr(2).rates(irs1) = (rr(2).rates(irppa) - rr(1).rates(irs1)*denomdt)*zz;
    }

    // s32(a,p)cl35(p,g)ar36
    rr(1).rates(irt1)  = 0.0e0_rt;
    rr(2).rates(irt1)  = 0.0e0_rt;
    denom    = rr(1).rates(irclpa) + rr(1).rates(irclpg);
    denomdt  = rr(2).rates(irclpa) + rr(2).rates(irclpg);
    if (denom > 1.0e-30_rt) {
       zz = 1.0e0_rt/denom;
       rr(1).rates(irt1) = rr(1).rates(irclpa)*zz;
       rr(2).rates(irt1) = (rr(2).rates(irclpa) - rr(1).rates(irt1)*denomdt)*zz;
    }

    // ar36(a,p)k39(p,g)ca40
    rr(1).rates(iru1)  = 0.0e0_rt;
    rr(2).rates(iru1)  = 0.0e0_rt;
    denom    = rr(1).rates(irkpa) + rr(1).rates(irkpg);
    denomdt  = rr(2).rates(irkpa) + rr(2).rates(irkpg);
    if (denom > 1.0e-30_rt) {
       zz   = 1.0e0_rt/denom;
       rr(1).rates(iru1)   = rr(1).rates(irkpa)*zz;
       rr(2).rates(iru1) = (rr(2).rates(irkpa) - rr(1).rates(iru1)*denomdt)*zz;
    }

    // ca40(a,p)sc43(p,g)ti44
    rr(1).rates(irv1)  = 0.0e0_rt;
    rr(2).rates(irv1)  = 0.0e0_rt;
    denom    = rr(1).rates(irscpa) + rr(1).rates(irscpg);
    denomdt  = rr(2).rates(irscpa) + rr(2).rates(irscpg);
    if (denom > 1.0e-30_rt) {
       zz  = 1.0e0_rt/denom;
       rr(1).rates(irv1) = rr(1).rates(irscpa)*zz;
       rr(2).rates(irv1) = (rr(2).rates(irscpa) - rr(1).rates(irv1)*denomdt)*zz;
    }

    // ti44(a,p)v47(p,g)cr48
    rr(1).rates(irw1) = 0.0e0_rt;
    rr(2).rates(irw1) = 0.0e0_rt;
    denom    = rr(1).rates(irvpa) + rr(1).rates(irvpg);
    denomdt  = rr(2).rates(irvpa) + rr(2).rates(irvpg);
    if (denom > 1.0e-30_rt) {
       zz = 1.0e0_rt/denom;
       rr(1).rates(irw1) = rr(1).rates(irvpa)*zz;
       rr(2).rates(irw1) = (rr(2).rates(irvpa) - rr(1).rates(irw1)*denomdt)*zz;
    }

    // cr48(a,p)mn51(p,g)fe52
    rr(1).rates(irx1) = 0.0e0_rt;
    rr(2).rates(irx1) = 0.0e0_rt;
    denom    = rr(1).rates(irmnpa) + rr(1).rates(irmnpg);
    denomdt  = rr(2).rates(irmnpa) + rr(2).rates(irmnpg);
    if (denom > 1.0e-30_rt) {
       zz = 1.0e0_rt/denom;
       rr(1).rates(irx1) = rr(1).rates(irmnpa)*zz;
       rr(2).rates(irx1) = (rr(2).rates(irmnpa) - rr(1).rates(irx1)*denomdt)*zz;
    }

    // fe52(a,p)co55(p,g)ni56
    rr(1).rates(iry1) = 0.0e0_rt;
    rr(2).rates(iry1) = 0.0e0_rt;
    denom    = rr(1).rates(ircopa) + rr(1).rates(ircopg);
    denomdt  = rr(2).rates(ircopa) + rr(2).rates(ircopg);
    if (denom > 1.0e-30_rt) {
       zz = 1.0e0_rt/denom;
       rr(1).rates(iry1) = rr(1).rates(ircopa)*zz;
       rr(2).rates(iry1) = (rr(2).rates(ircopa) - rr(1).rates(iry1)*denomdt)*zz;
    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void evaluate_rates(burn_t const& state, Array1D<rate_t, 1, Rates::NumGroups>& rr)
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
}


template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void dfdy_isotopes_aprox13(Array1D<Real, 1, NumSpec> const& y,
                           burn_t const& /*state*/, Array1D<rate_t, 1, Rates::NumGroups> const& rr,
                           MatrixType& jac)
{

    // this routine sets up the dense aprox13 jacobian for the isotopes

    using namespace Species;
    using namespace Rates;

    {
        // he4 jacobian elements
        // d(he4)/d(he4)
        Array1D<Real, 1, 20> b {
            -1.5e0_rt * y(He4) * y(He4) * rr(1).rates(ir3a),
            -y(C12)  * rr(1).rates(ircag),
            -y(O16)  * rr(1).rates(iroag),
            -y(Ne20) * rr(1).rates(irneag),
            -y(Mg24) * rr(1).rates(irmgag),
            -y(Si28) * rr(1).rates(irsiag),
            -y(S32)  * rr(1).rates(irsag),
            -y(Ar36) * rr(1).rates(irarag),
            -y(Ca40) * rr(1).rates(ircaag),
            -y(Ti44) * rr(1).rates(irtiag),
            -y(Cr48) * rr(1).rates(ircrag),
            -y(Fe52) * rr(1).rates(irfeag),
            -y(Mg24) * rr(1).rates(irmgap) * (1.0e0_rt-rr(1).rates(irr1)),
            -y(Si28) * rr(1).rates(irsiap) * (1.0e0_rt-rr(1).rates(irs1)),
            -y(S32)  * rr(1).rates(irsap)  * (1.0e0_rt-rr(1).rates(irt1)),
            -y(Ar36) * rr(1).rates(irarap) * (1.0e0_rt-rr(1).rates(iru1)),
            -y(Ca40) * rr(1).rates(ircaap) * (1.0e0_rt-rr(1).rates(irv1)),
            -y(Ti44) * rr(1).rates(irtiap) * (1.0e0_rt-rr(1).rates(irw1)),
            -y(Cr48) * rr(1).rates(ircrap) * (1.0e0_rt-rr(1).rates(irx1)),
            -y(Fe52) * rr(1).rates(irfeap) * (1.0e0_rt-rr(1).rates(iry1))
        };

        jac(He4,He4) = esum20(b);
    }

    {
        // (he4)/d(o16)
        Array1D<Real, 1, 5> b {
            0.5e0_rt * y(C12) * rr(1).rates(ir1216),
            1.12e0_rt * 0.5e0_rt*y(O16) * rr(1).rates(ir1616),
            0.68e0_rt * rr(1).rates(irs1) * 0.5e0_rt*y(O16) * rr(1).rates(ir1616),
            rr(1).rates(iroga),
            -y(He4) * rr(1).rates(iroag),
        };

        jac(He4,O16) = esum5(b);
    }

    {
        Array1D<Real, 1, 4> b;

        // d(he4)/d(c12)
        b(1) =  y(C12) * rr(1).rates(ir1212);
        b(2) =  0.5e0_rt * y(O16) * rr(1).rates(ir1216);
        b(3) =  3.0e0_rt * rr(1).rates(irg3a);
        b(4) = -y(He4) * rr(1).rates(ircag);

        jac(He4,C12) = esum4(b);

        // d(he4)/d(si28)
        b(1) =  rr(1).rates(irsiga);
        b(2) = -y(He4) * rr(1).rates(irsiag);
        b(3) = -y(He4) * rr(1).rates(irsiap) * (1.0e0_rt-rr(1).rates(irs1));
        b(4) =  rr(1).rates(irr1) * rr(1).rates(irsigp);

        jac(He4,Si28) = esum4(b);

        // d(he4)/d(s32)
        b(1) =  rr(1).rates(irsga);
        b(2) = -y(He4) * rr(1).rates(irsag);
        b(3) = -y(He4) * rr(1).rates(irsap) * (1.0e0_rt-rr(1).rates(irt1));
        b(4) =  rr(1).rates(irs1) * rr(1).rates(irsgp);

        jac(He4,S32) = esum4(b);

        // d(he4)/d(ar36)
        b(1)  =  rr(1).rates(irarga);
        b(2)  = -y(He4) * rr(1).rates(irarag);
        b(3)  = -y(He4) * rr(1).rates(irarap) * (1.0e0_rt-rr(1).rates(iru1));
        b(4)  =  rr(1).rates(irt1) * rr(1).rates(irargp);

        jac(He4,Ar36) = esum4(b);

        // d(he4)/d(ca40)
        b(1) =  rr(1).rates(ircaga);
        b(2) = -y(He4) * rr(1).rates(ircaag);
        b(3) = -y(He4) * rr(1).rates(ircaap) * (1.0e0_rt-rr(1).rates(irv1));
        b(4) =  rr(1).rates(iru1) * rr(1).rates(ircagp);

        jac(He4,Ca40) = esum4(b);

        // d(he4)/d(ti44)
        b(1) =  rr(1).rates(irtiga);
        b(2) = -y(He4) * rr(1).rates(irtiag);
        b(3) = -y(He4) * rr(1).rates(irtiap) * (1.0e0_rt-rr(1).rates(irw1));
        b(4) =  rr(1).rates(irv1) * rr(1).rates(irtigp);

        jac(He4,Ti44) = esum4(b);

        // d(he4)/d(cr48)
        b(1) =  rr(1).rates(ircrga);
        b(2) = -y(He4) * rr(1).rates(ircrag);
        b(3) = -y(He4) * rr(1).rates(ircrap) * (1.0e0_rt-rr(1).rates(irx1));
        b(4) =  rr(1).rates(irw1) * rr(1).rates(ircrgp);

        jac(He4,Cr48) = esum4(b);

        // d(he4)/d(fe52)
        b(1) =  rr(1).rates(irfega);
        b(2) = -y(He4) * rr(1).rates(irfeag);
        b(3) = -y(He4) * rr(1).rates(irfeap) * (1.0e0_rt-rr(1).rates(iry1));
        b(4) =  rr(1).rates(irx1) * rr(1).rates(irfegp);

        jac(He4,Fe52) = esum4(b);

        // d(c12)/d(c12)
        b(1) = -2.0e0_rt * y(C12) * rr(1).rates(ir1212);
        b(2) = -y(O16) * rr(1).rates(ir1216);
        b(3) = -rr(1).rates(irg3a);
        b(4) = -y(He4) * rr(1).rates(ircag);

        jac(C12,C12) = esum4(b);

        // d(o16)/d(o16)
        b(1) = -y(C12) * rr(1).rates(ir1216);
        b(2) = -2.0e0_rt * y(O16) * rr(1).rates(ir1616);
        b(3) = -y(He4) * rr(1).rates(iroag);
        b(4) = -rr(1).rates(iroga);

        jac(O16,O16) = esum4(b);

        // si28 jacobian elements
        // d(si28)/d(he4)
        b(1) =  y(Mg24) * rr(1).rates(irmgag);
        b(2) = -y(Si28) * rr(1).rates(irsiag);
        b(3) =  y(Mg24) * rr(1).rates(irmgap) * (1.0e0_rt-rr(1).rates(irr1));
        b(4) = -y(Si28) * rr(1).rates(irsiap) * (1.0e0_rt-rr(1).rates(irs1));

        jac(Si28,He4) = esum4(b);

        // d(si28)/d(si28)
        b(1) =  -y(He4) * rr(1).rates(irsiag);
        b(2) = -rr(1).rates(irsiga);
        b(3) = -rr(1).rates(irr1) * rr(1).rates(irsigp);
        b(4) = -y(He4) * rr(1).rates(irsiap) * (1.0e0_rt-rr(1).rates(irs1));

        jac(Si28,Si28) = esum4(b);

        // s32 jacobian elements
        // d(s32)/d(he4)
        b(1) =  y(Si28) * rr(1).rates(irsiag);
        b(2) = -y(S32) * rr(1).rates(irsag);
        b(3) =  y(Si28) * rr(1).rates(irsiap) * (1.0e0_rt-rr(1).rates(irs1));
        b(4) = -y(S32) * rr(1).rates(irsap) * (1.0e0_rt-rr(1).rates(irt1));

        jac(S32,He4) = esum4(b);

        // d(s32)/d(s32)
        b(1) = -y(He4) * rr(1).rates(irsag);
        b(2) = -rr(1).rates(irsga);
        b(3) = -rr(1).rates(irs1) * rr(1).rates(irsgp);
        b(4) = -y(He4) * rr(1).rates(irsap) * (1.0e0_rt-rr(1).rates(irt1));

        jac(S32,S32) = esum4(b);

        // ar36 jacobian elements
        // d(ar36)/d(he4)
        b(1) =  y(S32)  * rr(1).rates(irsag);
        b(2) = -y(Ar36) * rr(1).rates(irarag);
        b(3) =  y(S32)  * rr(1).rates(irsap) * (1.0e0_rt-rr(1).rates(irt1));
        b(4) = -y(Ar36) * rr(1).rates(irarap) * (1.0e0_rt-rr(1).rates(iru1));

        jac(Ar36,He4) = esum4(b);

        // d(ar36)/d(ar36)
        b(1) = -y(He4) * rr(1).rates(irarag);
        b(2) = -rr(1).rates(irarga);
        b(3) = -rr(1).rates(irt1) * rr(1).rates(irargp);
        b(4) = -y(He4) * rr(1).rates(irarap) * (1.0e0_rt-rr(1).rates(iru1));

        jac(Ar36,Ar36) = esum4(b);

        // ca40 jacobian elements
        // d(ca40)/d(he4)
        b(1)  =  y(Ar36) * rr(1).rates(irarag);
        b(2)  = -y(Ca40) * rr(1).rates(ircaag);
        b(3)  =  y(Ar36) * rr(1).rates(irarap)*(1.0e0_rt-rr(1).rates(iru1));
        b(4)  = -y(Ca40) * rr(1).rates(ircaap)*(1.0e0_rt-rr(1).rates(irv1));

        jac(Ca40,He4) = esum4(b);

        // d(ca40)/d(ca40)
        b(1) = -y(He4) * rr(1).rates(ircaag);
        b(2) = -rr(1).rates(ircaga);
        b(3) = -rr(1).rates(ircagp) * rr(1).rates(iru1);
        b(4) = -y(He4) * rr(1).rates(ircaap)*(1.0e0_rt-rr(1).rates(irv1));

        jac(Ca40,Ca40) = esum4(b);

        // ti44 jacobian elements
        // d(ti44)/d(he4)
        b(1) =  y(Ca40) * rr(1).rates(ircaag);
        b(2) = -y(Ti44) * rr(1).rates(irtiag);
        b(3) =  y(Ca40) * rr(1).rates(ircaap)*(1.0e0_rt-rr(1).rates(irv1));
        b(4) = -y(Ti44) * rr(1).rates(irtiap)*(1.0e0_rt-rr(1).rates(irw1));

        jac(Ti44,He4) = esum4(b);

        // d(ti44)/d(ti44)
        b(1) = -y(He4) * rr(1).rates(irtiag);
        b(2) = -rr(1).rates(irtiga);
        b(3) = -rr(1).rates(irv1) * rr(1).rates(irtigp);
        b(4) = -y(He4) * rr(1).rates(irtiap)*(1.0e0_rt-rr(1).rates(irw1));

        jac(Ti44,Ti44) = esum4(b);

        // cr48 jacobian elements
        // d(cr48)/d(he4)
        b(1) =  y(Ti44) * rr(1).rates(irtiag);
        b(2) = -y(Cr48) * rr(1).rates(ircrag);
        b(3) =  y(Ti44) * rr(1).rates(irtiap)*(1.0e0_rt-rr(1).rates(irw1));
        b(4) = -y(Cr48) * rr(1).rates(ircrap)*(1.0e0_rt-rr(1).rates(irx1));

        jac(Cr48,He4) = esum4(b);

        // d(cr48)/d(cr48)
        b(1) = -y(He4) * rr(1).rates(ircrag);
        b(2) = -rr(1).rates(ircrga);
        b(3) = -rr(1).rates(irw1) * rr(1).rates(ircrgp);
        b(4) = -y(He4) * rr(1).rates(ircrap)*(1.0e0_rt-rr(1).rates(irx1));

        jac(Cr48,Cr48) = esum4(b);

        // fe52 jacobian elements
        // d(fe52)/d(he4)
        b(1) =  y(Cr48) * rr(1).rates(ircrag);
        b(2) = -y(Fe52) * rr(1).rates(irfeag);
        b(3) =  y(Cr48) * rr(1).rates(ircrap) * (1.0e0_rt-rr(1).rates(irx1));
        b(4) = -y(Fe52) * rr(1).rates(irfeap) * (1.0e0_rt-rr(1).rates(iry1));

        jac(Fe52,He4) = esum4(b);

        // d(fe52)/d(fe52)
        b(1) = -y(He4) * rr(1).rates(irfeag);
        b(2) = -rr(1).rates(irfega);
        b(3) = -rr(1).rates(irx1) * rr(1).rates(irfegp);
        b(4) = -y(He4) * rr(1).rates(irfeap) * (1.0e0_rt-rr(1).rates(iry1));

        jac(Fe52,Fe52) = esum4(b);
    }

    {
        Array1D<Real, 1, 3> b;

        // d(he4)/d(mg24)
        b(1) =  rr(1).rates(irmgga);
        b(2) = -y(He4) * rr(1).rates(irmgag);
        b(3) = -y(He4) * rr(1).rates(irmgap) * (1.0e0_rt-rr(1).rates(irr1));

        jac(He4,Mg24) = esum3(b);

        // mg24 jacobian elements
        // d(mg24)/d(he4)
        b(1) =  y(Ne20) * rr(1).rates(irneag);
        b(2) = -y(Mg24) * rr(1).rates(irmgag);
        b(3) = -y(Mg24) * rr(1).rates(irmgap) * (1.0e0_rt-rr(1).rates(irr1));

        jac(Mg24,He4) = esum3(b);

        // d(mg24)/d(mg24)
        b(1) = -y(He4) * rr(1).rates(irmgag);
        b(2) = -rr(1).rates(irmgga);
        b(3) = -y(He4) * rr(1).rates(irmgap) * (1.0e0_rt-rr(1).rates(irr1));

        jac(Mg24,Mg24) = esum3(b);

        // d(si28)/d(o16)
        b(1) = 0.5e0_rt * y(C12) * rr(1).rates(ir1216);
        b(2) = 1.12e0_rt * 0.5e0_rt*y(O16) * rr(1).rates(ir1616);
        b(3) = 0.68e0_rt * 0.5e0_rt*y(O16) * rr(1).rates(irs1) * rr(1).rates(ir1616);

        jac(Si28,O16) = esum3(b);
    }

    {
        Array1D<Real, 1, 2> b;

        // d(he4)/d(ne20)
        b(1) =  rr(1).rates(irnega);
        b(2) = -y(He4) * rr(1).rates(irneag);

        jac(He4,Ne20) = b(1) + b(2);

        // d(he4)/d(ni56)
        b(1) = rr(1).rates(irniga);
        b(2) = rr(1).rates(iry1) * rr(1).rates(irnigp);

        jac(He4,Ni56) = b(1) + b(2);

        // c12 jacobian elements
        // d(c12)/d(he4)
        b(1) =  0.5e0_rt * y(He4) * y(He4) * rr(1).rates(ir3a);
        b(2) = -y(C12) * rr(1).rates(ircag);

        jac(C12,He4) = b(1) + b(2);

        // d(c12)/d(o16)
        b(1) = -y(C12) * rr(1).rates(ir1216);
        b(2) =  rr(1).rates(iroga);

        jac(C12,O16) = b(1) + b(2);

        // o16 jacobian elements
        // d(o16)/d(he4)
        b(1) =  y(C12)*rr(1).rates(ircag);
        b(2) = -y(O16)*rr(1).rates(iroag);

        jac(O16,He4) = b(1) + b(2);

        // d(o16)/d(c12)
        b(1) = -y(O16)*rr(1).rates(ir1216);
        b(2) =  y(He4)*rr(1).rates(ircag);

        jac(O16,C12) = b(1) + b(2);

        // ne20 jacobian elements
        // d(ne20)/d(he4)
        b(1) =  y(O16) * rr(1).rates(iroag);
        b(2) = -y(Ne20) * rr(1).rates(irneag);

        jac(Ne20,He4) = b(1) + b(2);

        // d(ne20)/d(ne20)
        b(1) = -y(He4) * rr(1).rates(irneag);
        b(2) = -rr(1).rates(irnega);

        jac(Ne20,Ne20) = b(1) + b(2);

        // d(mg24)/d(si28)
        b(1) = rr(1).rates(irsiga);
        b(2) = rr(1).rates(irr1) * rr(1).rates(irsigp);

        jac(Mg24,Si28) = b(1) + b(2);

        // d(si28)/d(mg24)
        b(1) =  y(He4) * rr(1).rates(irmgag);
        b(2) =  y(He4) * rr(1).rates(irmgap) * (1.0e0_rt-rr(1).rates(irr1));

        jac(Si28,Mg24) = b(1) + b(2);

        // d(si28)/d(s32)
        b(1) = rr(1).rates(irsga);
        b(2) = rr(1).rates(irs1) * rr(1).rates(irsgp);

        jac(Si28,S32) = b(1) + b(2);

        // d(s32)/d(o16)
        b(1) = 0.68e0_rt*0.5e0_rt*y(O16)*rr(1).rates(ir1616)*(1.0e0_rt-rr(1).rates(irs1));
        b(2) = 0.2e0_rt * 0.5e0_rt*y(O16) * rr(1).rates(ir1616);

        jac(S32,O16) = b(1) + b(2);

        // d(s32)/d(si28)
        b(1)  =y(He4) * rr(1).rates(irsiag);
        b(2) = y(He4) * rr(1).rates(irsiap) * (1.0e0_rt-rr(1).rates(irs1));

        jac(S32,Si28) = b(1) + b(2);

        // d(s32)/d(ar36)
        b(1) = rr(1).rates(irarga);
        b(2) = rr(1).rates(irt1) * rr(1).rates(irargp);

        jac(S32,Ar36) = b(1) + b(2);

        // d(ar36)/d(s32)
        b(1) = y(He4) * rr(1).rates(irsag);
        b(2) = y(He4) * rr(1).rates(irsap) * (1.0e0_rt-rr(1).rates(irt1));

        jac(Ar36,S32) = b(1) + b(2);

        // d(ar36)/d(ca40)
        b(1) = rr(1).rates(ircaga);
        b(2) = rr(1).rates(ircagp) * rr(1).rates(iru1);

        jac(Ar36,Ca40) = b(1) + b(2);

        // d(ca40)/d(ar36)
        b(1) =  y(He4) * rr(1).rates(irarag);
        b(2) =  y(He4) * rr(1).rates(irarap)*(1.0e0_rt-rr(1).rates(iru1));

        jac(Ca40,Ar36) = b(1) + b(2);

        // d(ca40)/d(ti44)
        b(1) = rr(1).rates(irtiga);
        b(2) = rr(1).rates(irtigp) * rr(1).rates(irv1);

        jac(Ca40,Ti44) = b(1) + b(2);

        // d(ti44)/d(ca40)
        b(1) =  y(He4) * rr(1).rates(ircaag);
        b(2) =  y(He4) * rr(1).rates(ircaap)*(1.0e0_rt-rr(1).rates(irv1));

        jac(Ti44,Ca40) = b(1) + b(2);

        // d(ti44)/d(cr48)
        b(1) = rr(1).rates(ircrga);
        b(2) = rr(1).rates(irw1) * rr(1).rates(ircrgp);

        jac(Ti44,Cr48) = b(1) + b(2);

        // d(cr48)/d(ti44)
        b(1) =  y(He4) * rr(1).rates(irtiag);
        b(2) =  y(He4) * rr(1).rates(irtiap)*(1.0e0_rt-rr(1).rates(irw1));

        jac(Cr48,Ti44) = b(1) + b(2);

        // d(cr48)/d(fe52)
        b(1) = rr(1).rates(irfega);
        b(2) = rr(1).rates(irx1) * rr(1).rates(irfegp);

        jac(Cr48,Fe52) = b(1) + b(2);

        // d(fe52)/d(cr48)
        b(1) = y(He4) * rr(1).rates(ircrag);
        b(2) = y(He4) * rr(1).rates(ircrap) * (1.0e0_rt-rr(1).rates(irx1));

        jac(Fe52,Cr48) = b(1) + b(2);

        // d(fe52)/d(ni56)
        b(1) = rr(1).rates(irniga);
        b(2) = rr(1).rates(iry1) * rr(1).rates(irnigp);

        jac(Fe52,Ni56) = b(1) + b(2);

        // ni56 jacobian elements
        // d(ni56)/d(he4)
        b(1) =  y(Fe52) * rr(1).rates(irfeag);
        b(2) =  y(Fe52) * rr(1).rates(irfeap) * (1.0e0_rt-rr(1).rates(iry1));

        jac(Ni56,He4) = b(1) + b(2);

        // d(ni56)/d(fe52)
        b(1) = y(He4) * rr(1).rates(irfeag);
        b(2) = y(He4) * rr(1).rates(irfeap) * (1.0e0_rt-rr(1).rates(iry1));

        jac(Ni56,Fe52) = b(1) + b(2);

        // d(ni56)/d(ni56)
        b(1) = -rr(1).rates(irniga);
        b(2) = -rr(1).rates(iry1) * rr(1).rates(irnigp);

        jac(Ni56,Ni56) = b(1) + b(2);
    }

    // d(o16)/d(ne20)
    jac(O16,Ne20) = rr(1).rates(irnega);

    // d(ne20)/d(c12)
    jac(Ne20,C12) = y(C12) * rr(1).rates(ir1212);

    // d(ne20)/d(o16)
    jac(Ne20,O16) = y(He4) * rr(1).rates(iroag);

    // d(ne20)/d(mg24)
    jac(Ne20,Mg24) = rr(1).rates(irmgga);

    // d(mg24)/d(c12)
    jac(Mg24,C12) = 0.5e0_rt * y(O16) * rr(1).rates(ir1216);

    // d(mg24)/d(o16)
    jac(Mg24,O16) = 0.5e0_rt * y(C12) * rr(1).rates(ir1216);

    // d(mg24)/d(ne20)
    jac(Mg24,Ne20) = y(He4) * rr(1).rates(irneag);

    // d(si28)/d(c12)
    jac(Si28,C12) = 0.5e0_rt * y(O16) * rr(1).rates(ir1216);

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
void rhs(Array1D<Real, 1, NumSpec> const& y, Array1D<rate_t, 1, Rates::NumGroups> const& rr,
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
    a(1)  = 0.5e0_rt * y(C12) * y(C12) * rr(index_rate).rates(ir1212);
    a(2)  = 0.5e0_rt * y(C12) * y(O16) * rr(index_rate).rates(ir1216);
    a(3)  = 0.56e0_rt * 0.5e0_rt * y(O16) * y(O16) * rr(index_rate).rates(ir1616);

    dydt(He4) = dydt(He4) + esum3(a);

    // (a,g) and (g,a) reactions
    a(1)  = -0.5e0_rt * y(He4) * y(He4) * y(He4) * rr(index_rate).rates(ir3a);
    a(2)  =  3.0e0_rt * y(C12) * rr(index_rate).rates(irg3a);
    a(3)  = -y(He4)  * y(C12) * rr(index_rate).rates(ircag);
    a(4)  =  y(O16)  * rr(index_rate).rates(iroga);
    a(5)  = -y(He4)  * y(O16) * rr(index_rate).rates(iroag);
    a(6)  =  y(Ne20) * rr(index_rate).rates(irnega);
    a(7)  = -y(He4)  * y(Ne20) * rr(index_rate).rates(irneag);
    a(8)  =  y(Mg24) * rr(index_rate).rates(irmgga);
    a(9)  = -y(He4)  * y(Mg24)* rr(index_rate).rates(irmgag);
    a(10) =  y(Si28) * rr(index_rate).rates(irsiga);
    a(11) = -y(He4)  * y(Si28)*rr(index_rate).rates(irsiag);
    a(12) =  y(S32)  * rr(index_rate).rates(irsga);

    dydt(He4) = dydt(He4) + esum12(a);

    a(1)  = -y(He4)  * y(S32) * rr(index_rate).rates(irsag);
    a(2)  =  y(Ar36) * rr(index_rate).rates(irarga);
    a(3)  = -y(He4)  * y(Ar36)*rr(index_rate).rates(irarag);
    a(4)  =  y(Ca40) * rr(index_rate).rates(ircaga);
    a(5)  = -y(He4)  * y(Ca40)*rr(index_rate).rates(ircaag);
    a(6)  =  y(Ti44) * rr(index_rate).rates(irtiga);
    a(7)  = -y(He4)  * y(Ti44)*rr(index_rate).rates(irtiag);
    a(8)  =  y(Cr48) * rr(index_rate).rates(ircrga);
    a(9)  = -y(He4)  * y(Cr48)*rr(index_rate).rates(ircrag);
    a(10) =  y(Fe52) * rr(index_rate).rates(irfega);
    a(11) = -y(He4)  * y(Fe52) * rr(index_rate).rates(irfeag);
    a(12) =  y(Ni56) * rr(index_rate).rates(irniga);

    dydt(He4) = dydt(He4) + esum12(a);

    // (a,p)(p,g) and (g,p)(p,a) reactions

    if (!deriva) {

       a(1)  =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16)*rr(index_rate).rates(irs1)*rr(index_rate).rates(ir1616);
       a(2)  = -y(He4)  * y(Mg24) * rr(index_rate).rates(irmgap)*(1.0e0_rt-rr(index_rate).rates(irr1));
       a(3)  =  y(Si28) * rr(index_rate).rates(irsigp) * rr(index_rate).rates(irr1);
       a(4)  = -y(He4)  * y(Si28) * rr(index_rate).rates(irsiap)*(1.0e0_rt-rr(index_rate).rates(irs1));
       a(5)  =  y(S32)  * rr(index_rate).rates(irsgp) * rr(index_rate).rates(irs1);
       a(6)  = -y(He4)  * y(S32) * rr(index_rate).rates(irsap)*(1.0e0_rt-rr(index_rate).rates(irt1));
       a(7)  =  y(Ar36) * rr(index_rate).rates(irargp) * rr(index_rate).rates(irt1);
       a(8)  = -y(He4)  * y(Ar36) * rr(index_rate).rates(irarap)*(1.0e0_rt-rr(index_rate).rates(iru1));
       a(9)  =  y(Ca40) * rr(index_rate).rates(ircagp) * rr(index_rate).rates(iru1);
       a(10) = -y(He4)  * y(Ca40) * rr(index_rate).rates(ircaap)*(1.0e0_rt-rr(index_rate).rates(irv1));
       a(11) =  y(Ti44) * rr(index_rate).rates(irtigp) * rr(index_rate).rates(irv1);
       a(12) = -y(He4)  * y(Ti44) * rr(index_rate).rates(irtiap)*(1.0e0_rt-rr(index_rate).rates(irw1));
       a(13) =  y(Cr48) * rr(index_rate).rates(ircrgp) * rr(index_rate).rates(irw1);
       a(14) = -y(He4)  * y(Cr48) * rr(index_rate).rates(ircrap)*(1.0e0_rt-rr(index_rate).rates(irx1));
       a(15) =  y(Fe52) * rr(index_rate).rates(irfegp) * rr(index_rate).rates(irx1);
       a(16) = -y(He4)  * y(Fe52) * rr(index_rate).rates(irfeap)*(1.0e0_rt-rr(index_rate).rates(iry1));
       a(17) =  y(Ni56) * rr(index_rate).rates(irnigp) * rr(index_rate).rates(iry1);

       dydt(He4) = dydt(He4) + esum17(a);

    } else {

       a(1)  =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rr(1).rates(irs1) * rr(index_rate).rates(ir1616);
       a(2)  =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rr(index_rate).rates(irs1) * rr(1).rates(ir1616);
       a(3)  = -y(He4)*y(Mg24) * rr(index_rate).rates(irmgap)*(1.0e0_rt - rr(1).rates(irr1));
       a(4)  =  y(He4)*y(Mg24) * rr(1).rates(irmgap)*rr(index_rate).rates(irr1);
       a(5)  =  y(Si28) * rr(1).rates(irsigp) * rr(index_rate).rates(irr1);
       a(6)  =  y(Si28) * rr(index_rate).rates(irsigp) * rr(1).rates(irr1);
       a(7)  = -y(He4)*y(Si28) * rr(index_rate).rates(irsiap)*(1.0e0_rt - rr(1).rates(irs1));
       a(8)  =  y(He4)*y(Si28) * rr(1).rates(irsiap) * rr(index_rate).rates(irs1);
       a(9)  =  y(S32)  * rr(1).rates(irsgp) * rr(index_rate).rates(irs1);
       a(10) =  y(S32)  * rr(index_rate).rates(irsgp) * rr(1).rates(irs1);

       dydt(He4) = dydt(He4) + esum10(a);

       a(1)  = -y(He4)*y(S32) * rr(index_rate).rates(irsap)*(1.0e0_rt - rr(1).rates(irt1));
       a(2)  =  y(He4)*y(S32) * rr(1).rates(irsap)*rr(index_rate).rates(irt1);
       a(3)  =  y(Ar36) * rr(1).rates(irargp) * rr(index_rate).rates(irt1);
       a(4)  =  y(Ar36) * rr(index_rate).rates(irargp) * rr(1).rates(irt1);
       a(5)  = -y(He4)*y(Ar36) * rr(index_rate).rates(irarap)*(1.0e0_rt - rr(1).rates(iru1));
       a(6)  =  y(He4)*y(Ar36) * rr(1).rates(irarap)*rr(index_rate).rates(iru1);
       a(7)  =  y(Ca40) * rr(1).rates(ircagp) * rr(index_rate).rates(iru1);
       a(8)  =  y(Ca40) * rr(index_rate).rates(ircagp) * rr(1).rates(iru1);
       a(9)  = -y(He4)*y(Ca40) * rr(index_rate).rates(ircaap)*(1.0e0_rt-rr(1).rates(irv1));
       a(10) =  y(He4)*y(Ca40) * rr(1).rates(ircaap)*rr(index_rate).rates(irv1);
       a(11) =  y(Ti44) * rr(1).rates(irtigp) * rr(index_rate).rates(irv1);
       a(12) =  y(Ti44) * rr(index_rate).rates(irtigp) * rr(1).rates(irv1);

       dydt(He4) = dydt(He4) + esum12(a);

       a(1)  = -y(He4)*y(Ti44) * rr(index_rate).rates(irtiap)*(1.0e0_rt - rr(1).rates(irw1));
       a(2)  =  y(He4)*y(Ti44) * rr(1).rates(irtiap)*rr(index_rate).rates(irw1);
       a(3)  =  y(Cr48) * rr(1).rates(ircrgp) * rr(index_rate).rates(irw1);
       a(4)  =  y(Cr48) * rr(index_rate).rates(ircrgp) * rr(1).rates(irw1);
       a(5)  = -y(He4)*y(Cr48) * rr(index_rate).rates(ircrap)*(1.0e0_rt - rr(1).rates(irx1));
       a(6)  =  y(He4)*y(Cr48) * rr(1).rates(ircrap)*rr(index_rate).rates(irx1);
       a(7)  =  y(Fe52) * rr(1).rates(irfegp) * rr(index_rate).rates(irx1);
       a(8)  =  y(Fe52) * rr(index_rate).rates(irfegp) * rr(1).rates(irx1);
       a(9)  = -y(He4)*y(Fe52) * rr(index_rate).rates(irfeap)*(1.0e0_rt - rr(1).rates(iry1));
       a(10) =  y(He4)*y(Fe52) * rr(1).rates(irfeap)*rr(index_rate).rates(iry1);
       a(11) =  y(Ni56) * rr(1).rates(irnigp) * rr(index_rate).rates(iry1);
       a(12) =  y(Ni56) * rr(index_rate).rates(irnigp) * rr(1).rates(iry1);

       dydt(He4) = dydt(He4) + esum12(a);
    
    }

    // c12 reactions
    a(1) = -y(C12) * y(C12) * rr(index_rate).rates(ir1212);
    a(2) = -y(C12) * y(O16) * rr(index_rate).rates(ir1216);
    a(3) =  y(He4) * y(He4) * y(He4) * rr(index_rate).rates(ir3a) / 6.0_rt;;
    a(4) = -y(C12) * rr(index_rate).rates(irg3a);
    a(5) = -y(C12) * y(He4) * rr(index_rate).rates(ircag);
    a(6) =  y(O16) * rr(index_rate).rates(iroga);

    dydt(C12) = dydt(C12) + esum6(a);


    // o16 reactions
    a(1) = -y(C12) * y(O16) * rr(index_rate).rates(ir1216);
    a(2) = -y(O16) * y(O16) * rr(index_rate).rates(ir1616);
    a(3) =  y(C12) * y(He4) * rr(index_rate).rates(ircag);
    a(4) = -y(O16) * y(He4) * rr(index_rate).rates(iroag);
    a(5) = -y(O16) * rr(index_rate).rates(iroga);
    a(6) =  y(Ne20) * rr(index_rate).rates(irnega);

    dydt(O16) = dydt(O16) + esum6(a);


    // ne20 reactions
    a(1) =  0.5e0_rt * y(C12) * y(C12) * rr(index_rate).rates(ir1212);
    a(2) =  y(O16) * y(He4) * rr(index_rate).rates(iroag);
    a(3) = -y(Ne20) * y(He4) * rr(index_rate).rates(irneag);
    a(4) = -y(Ne20) * rr(index_rate).rates(irnega);
    a(5) =  y(Mg24) * rr(index_rate).rates(irmgga);

    dydt(Ne20) = dydt(Ne20) + esum5(a);


    // mg24 reactions
    a(1) =  0.5e0_rt * y(C12) * y(O16) * rr(index_rate).rates(ir1216);
    a(2) =  y(Ne20) * y(He4) * rr(index_rate).rates(irneag);
    a(3) = -y(Mg24) * y(He4) * rr(index_rate).rates(irmgag);
    a(4) = -y(Mg24) * rr(index_rate).rates(irmgga);
    a(5) =  y(Si28) * rr(index_rate).rates(irsiga);

    dydt(Mg24) = dydt(Mg24) + esum5(a);

    if (!deriva) {

       a(1) = -y(Mg24) * y(He4) * rr(index_rate).rates(irmgap)*(1.0e0_rt-rr(index_rate).rates(irr1));
       a(2) =  y(Si28) * rr(index_rate).rates(irr1) * rr(index_rate).rates(irsigp);

       dydt(Mg24) = dydt(Mg24) + a(1) + a(2);

    } else {

       a(1) = -y(Mg24)*y(He4) * rr(index_rate).rates(irmgap)*(1.0e0_rt - rr(1).rates(irr1));
       a(2) =  y(Mg24)*y(He4) * rr(1).rates(irmgap)*rr(index_rate).rates(irr1);
       a(3) =  y(Si28) * rr(1).rates(irr1) * rr(index_rate).rates(irsigp);
       a(4) =  y(Si28) * rr(index_rate).rates(irr1) * rr(1).rates(irsigp);

       dydt(Mg24) = dydt(Mg24) + esum4(a);

    }


    // si28 reactions
    a(1) =  0.5e0_rt * y(C12) * y(O16) * rr(index_rate).rates(ir1216);
    a(2) =  0.56e0_rt * 0.5e0_rt*y(O16) * y(O16) * rr(index_rate).rates(ir1616);
    a(3) =  y(Mg24) * y(He4) * rr(index_rate).rates(irmgag);
    a(4) = -y(Si28) * y(He4) * rr(index_rate).rates(irsiag);
    a(5) = -y(Si28) * rr(index_rate).rates(irsiga);
    a(6) =  y(S32)  * rr(index_rate).rates(irsga);

    dydt(Si28) = dydt(Si28) + esum6(a);

    if (!deriva) {

       a(1) =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16)*rr(index_rate).rates(irs1)*rr(index_rate).rates(ir1616);
       a(2) =  y(Mg24) * y(He4) * rr(index_rate).rates(irmgap)*(1.0e0_rt-rr(index_rate).rates(irr1));
       a(3) = -y(Si28) * rr(index_rate).rates(irr1) * rr(index_rate).rates(irsigp);
       a(4) = -y(Si28) * y(He4) * rr(index_rate).rates(irsiap)*(1.0e0_rt-rr(index_rate).rates(irs1));
       a(5) =  y(S32)  * rr(index_rate).rates(irs1) * rr(index_rate).rates(irsgp);

       dydt(Si28) = dydt(Si28) + esum5(a);

    } else {

       a(1)  =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rr(1).rates(irs1)*rr(index_rate).rates(ir1616);
       a(2)  =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rr(index_rate).rates(irs1)*rr(1).rates(ir1616);
       a(3)  =  y(Mg24)*y(He4) * rr(index_rate).rates(irmgap)*(1.0e0_rt - rr(1).rates(irr1));
       a(4)  = -y(Mg24)*y(He4) * rr(1).rates(irmgap)*rr(index_rate).rates(irr1);
       a(5)  = -y(Si28) * rr(1).rates(irr1) * rr(index_rate).rates(irsigp);
       a(6)  = -y(Si28) * rr(index_rate).rates(irr1) * rr(1).rates(irsigp);
       a(7)  = -y(Si28)*y(He4) * rr(index_rate).rates(irsiap)*(1.0e0_rt - rr(1).rates(irs1));
       a(8)  =  y(Si28)*y(He4) * rr(1).rates(irsiap)*rr(index_rate).rates(irs1);
       a(9)  = y(S32) * rr(1).rates(irs1) * rr(index_rate).rates(irsgp);
       a(10) = y(S32) * rr(index_rate).rates(irs1) * rr(1).rates(irsgp);

       dydt(Si28) = dydt(Si28) + esum10(a);

    }


    // s32 reactions
    a(1) =  0.1e0_rt * 0.5e0_rt*y(O16) * y(O16) * rr(index_rate).rates(ir1616);
    a(2) =  y(Si28) * y(He4) * rr(index_rate).rates(irsiag);
    a(3) = -y(S32) * y(He4) * rr(index_rate).rates(irsag);
    a(4) = -y(S32) * rr(index_rate).rates(irsga);
    a(5) =  y(Ar36) * rr(index_rate).rates(irarga);

    dydt(S32) = dydt(S32) + esum5(a);


    if (!deriva) {

       a(1) =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16)* rr(index_rate).rates(ir1616)*(1.0e0_rt-rr(index_rate).rates(irs1));
       a(2) =  y(Si28) * y(He4) * rr(index_rate).rates(irsiap)*(1.0e0_rt-rr(index_rate).rates(irs1));
       a(3) = -y(S32) * rr(index_rate).rates(irs1) * rr(index_rate).rates(irsgp);
       a(4) = -y(S32) * y(He4) * rr(index_rate).rates(irsap)*(1.0e0_rt-rr(index_rate).rates(irt1));
       a(5) =  y(Ar36) * rr(index_rate).rates(irt1) * rr(index_rate).rates(irargp);

       dydt(S32) = dydt(S32) + esum5(a);

    } else {

       a(1)  =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rr(index_rate).rates(ir1616)*(1.0e0_rt-rr(1).rates(irs1));
       a(2)  = -0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rr(1).rates(ir1616)*rr(index_rate).rates(irs1);
       a(3)  =  y(Si28)*y(He4) * rr(index_rate).rates(irsiap)*(1.0e0_rt-rr(1).rates(irs1));
       a(4)  = -y(Si28)*y(He4) * rr(1).rates(irsiap)*rr(index_rate).rates(irs1);
       a(5)  = -y(S32) * rr(1).rates(irs1) * rr(index_rate).rates(irsgp);
       a(6)  = -y(S32) * rr(index_rate).rates(irs1) * rr(1).rates(irsgp);
       a(7)  = -y(S32)*y(He4) * rr(index_rate).rates(irsap)*(1.0e0_rt-rr(1).rates(irt1));
       a(8)  =  y(S32)*y(He4) * rr(1).rates(irsap)*rr(index_rate).rates(irt1);
       a(9)  =  y(Ar36) * rr(1).rates(irt1) * rr(index_rate).rates(irargp);
       a(10) =  y(Ar36) * rr(index_rate).rates(irt1) * rr(1).rates(irargp);

       dydt(S32) = dydt(S32) + esum10(a);

    }


    // ar36 reactions
    a(1) =  y(S32)  * y(He4) * rr(index_rate).rates(irsag);
    a(2) = -y(Ar36) * y(He4) * rr(index_rate).rates(irarag);
    a(3) = -y(Ar36) * rr(index_rate).rates(irarga);
    a(4) =  y(Ca40) * rr(index_rate).rates(ircaga);

    dydt(Ar36) = dydt(Ar36) + esum4(a);

    if (!deriva) {

       a(1) = y(S32)  * y(He4) * rr(index_rate).rates(irsap)*(1.0e0_rt-rr(index_rate).rates(irt1));
       a(2) = -y(Ar36) * rr(index_rate).rates(irt1) * rr(index_rate).rates(irargp);
       a(3) = -y(Ar36) * y(He4) * rr(index_rate).rates(irarap)*(1.0e0_rt-rr(index_rate).rates(iru1));
       a(4) =  y(Ca40) * rr(index_rate).rates(ircagp) * rr(index_rate).rates(iru1);

       dydt(Ar36) = dydt(Ar36) + esum4(a);

    } else {

       a(1) =  y(S32)*y(He4) * rr(index_rate).rates(irsap)*(1.0e0_rt - rr(1).rates(irt1));
       a(2) = -y(S32)*y(He4) * rr(1).rates(irsap)*rr(index_rate).rates(irt1);
       a(3) = -y(Ar36) * rr(1).rates(irt1) * rr(index_rate).rates(irargp);
       a(4) = -y(Ar36) * rr(index_rate).rates(irt1) * rr(1).rates(irargp);
       a(5) = -y(Ar36)*y(He4) * rr(index_rate).rates(irarap)*(1.0e0_rt-rr(1).rates(iru1));
       a(6) =  y(Ar36)*y(He4) * rr(1).rates(irarap)*rr(index_rate).rates(iru1);
       a(7) =  y(Ca40) * rr(1).rates(ircagp) * rr(index_rate).rates(iru1);
       a(8) =  y(Ca40) * rr(index_rate).rates(ircagp) * rr(1).rates(iru1);

       dydt(Ar36) = dydt(Ar36) + esum8(a);

    }


    // ca40 reactions
    a(1) =  y(Ar36) * y(He4) * rr(index_rate).rates(irarag);
    a(2) = -y(Ca40) * y(He4) * rr(index_rate).rates(ircaag);
    a(3) = -y(Ca40) * rr(index_rate).rates(ircaga);
    a(4) =  y(Ti44) * rr(index_rate).rates(irtiga);

    dydt(Ca40) = dydt(Ca40) + esum4(a);

    if (!deriva) {

       a(1) =  y(Ar36) * y(He4) * rr(index_rate).rates(irarap)*(1.0e0_rt-rr(index_rate).rates(iru1));
       a(2) = -y(Ca40) * rr(index_rate).rates(ircagp) * rr(index_rate).rates(iru1);
       a(3) = -y(Ca40) * y(He4) * rr(index_rate).rates(ircaap)*(1.0e0_rt-rr(index_rate).rates(irv1));
       a(4) =  y(Ti44) * rr(index_rate).rates(irtigp) * rr(index_rate).rates(irv1);

       dydt(Ca40) = dydt(Ca40) + esum4(a);

    } else {

       a(1) =  y(Ar36)*y(He4) * rr(index_rate).rates(irarap)*(1.0e0_rt-rr(1).rates(iru1));
       a(2) = -y(Ar36)*y(He4) * rr(1).rates(irarap)*rr(index_rate).rates(iru1);
       a(3) = -y(Ca40) * rr(1).rates(ircagp) * rr(index_rate).rates(iru1);
       a(4) = -y(Ca40) * rr(index_rate).rates(ircagp) * rr(1).rates(iru1);
       a(5) = -y(Ca40)*y(He4) * rr(index_rate).rates(ircaap)*(1.0e0_rt-rr(1).rates(irv1));
       a(6) =  y(Ca40)*y(He4) * rr(1).rates(ircaap)*rr(index_rate).rates(irv1);
       a(7) =  y(Ti44) * rr(1).rates(irtigp) * rr(index_rate).rates(irv1);
       a(8) =  y(Ti44) * rr(index_rate).rates(irtigp) * rr(1).rates(irv1);

       dydt(Ca40) = dydt(Ca40) + esum8(a);

    }


    // ti44 reactions
    a(1) =  y(Ca40) * y(He4) * rr(index_rate).rates(ircaag);
    a(2) = -y(Ti44) * y(He4) * rr(index_rate).rates(irtiag);
    a(3) = -y(Ti44) * rr(index_rate).rates(irtiga);
    a(4) =  y(Cr48) * rr(index_rate).rates(ircrga);

    dydt(Ti44) = dydt(Ti44) + esum4(a);

    if (!deriva) {

       a(1) =  y(Ca40) * y(He4) * rr(index_rate).rates(ircaap)*(1.0e0_rt-rr(index_rate).rates(irv1));
       a(2) = -y(Ti44) * rr(index_rate).rates(irv1) * rr(index_rate).rates(irtigp);
       a(3) = -y(Ti44) * y(He4) * rr(index_rate).rates(irtiap)*(1.0e0_rt-rr(index_rate).rates(irw1));
       a(4) =  y(Cr48) * rr(index_rate).rates(irw1) * rr(index_rate).rates(ircrgp);

       dydt(Ti44) = dydt(Ti44) + esum4(a);

    } else {

       a(1) =  y(Ca40)*y(He4) * rr(index_rate).rates(ircaap)*(1.0e0_rt-rr(1).rates(irv1));
       a(2) = -y(Ca40)*y(He4) * rr(1).rates(ircaap)*rr(index_rate).rates(irv1);
       a(3) = -y(Ti44) * rr(1).rates(irv1) * rr(index_rate).rates(irtigp);
       a(4) = -y(Ti44) * rr(index_rate).rates(irv1) * rr(1).rates(irtigp);
       a(5) = -y(Ti44)*y(He4) * rr(index_rate).rates(irtiap)*(1.0e0_rt-rr(1).rates(irw1));
       a(6) =  y(Ti44)*y(He4) * rr(1).rates(irtiap)*rr(index_rate).rates(irw1);
       a(7) =  y(Cr48) * rr(1).rates(irw1) * rr(index_rate).rates(ircrgp);
       a(8) =  y(Cr48) * rr(index_rate).rates(irw1) * rr(1).rates(ircrgp);

       dydt(Ti44) = dydt(Ti44) + esum8(a);

    }


    // cr48 reactions
    a(1) =  y(Ti44) * y(He4) * rr(index_rate).rates(irtiag);
    a(2) = -y(Cr48) * y(He4) * rr(index_rate).rates(ircrag);
    a(3) = -y(Cr48) * rr(index_rate).rates(ircrga);
    a(4) =  y(Fe52) * rr(index_rate).rates(irfega);

    dydt(Cr48) = dydt(Cr48) + esum4(a);

    if (!deriva) {

       a(1) =  y(Ti44) * y(He4) * rr(index_rate).rates(irtiap)*(1.0e0_rt-rr(index_rate).rates(irw1));
       a(2) = -y(Cr48) * rr(index_rate).rates(irw1) * rr(index_rate).rates(ircrgp);
       a(3) = -y(Cr48) * y(He4) * rr(index_rate).rates(ircrap)*(1.0e0_rt-rr(index_rate).rates(irx1));
       a(4) =  y(Fe52) * rr(index_rate).rates(irx1) * rr(index_rate).rates(irfegp);

       dydt(Cr48) = dydt(Cr48) + esum4(a);

    } else {

       a(1) =  y(Ti44)*y(He4) * rr(index_rate).rates(irtiap)*(1.0e0_rt-rr(1).rates(irw1));
       a(2) = -y(Ti44)*y(He4) * rr(1).rates(irtiap)*rr(index_rate).rates(irw1);
       a(3) = -y(Cr48) * rr(1).rates(irw1) * rr(index_rate).rates(ircrgp);
       a(4) = -y(Cr48) * rr(index_rate).rates(irw1) * rr(1).rates(ircrgp);
       a(5) = -y(Cr48)*y(He4) * rr(index_rate).rates(ircrap)*(1.0e0_rt-rr(1).rates(irx1));
       a(6) =  y(Cr48)*y(He4) * rr(1).rates(ircrap)*rr(index_rate).rates(irx1);
       a(7) =  y(Fe52) * rr(1).rates(irx1) * rr(index_rate).rates(irfegp);
       a(8) =  y(Fe52) * rr(index_rate).rates(irx1) * rr(1).rates(irfegp);

       dydt(Cr48) = dydt(Cr48) + esum8(a);

    }


    // fe52 reactions
    a(1) =  y(Cr48) * y(He4) * rr(index_rate).rates(ircrag);
    a(2) = -y(Fe52) * y(He4) * rr(index_rate).rates(irfeag);
    a(3) = -y(Fe52) * rr(index_rate).rates(irfega);
    a(4) =  y(Ni56) * rr(index_rate).rates(irniga);

    dydt(Fe52) = dydt(Fe52) + esum4(a);

    if (!deriva) {

       a(1) =  y(Cr48) * y(He4) * rr(index_rate).rates(ircrap)*(1.0e0_rt-rr(index_rate).rates(irx1));
       a(2) = -y(Fe52) * rr(index_rate).rates(irx1) * rr(index_rate).rates(irfegp);
       a(3) = -y(Fe52) * y(He4) * rr(index_rate).rates(irfeap)*(1.0e0_rt-rr(index_rate).rates(iry1));
       a(4) =  y(Ni56) * rr(index_rate).rates(iry1) * rr(index_rate).rates(irnigp);

       dydt(Fe52) = dydt(Fe52) + esum4(a);

    } else {

       a(1) =  y(Cr48)*y(He4) * rr(index_rate).rates(ircrap)*(1.0e0_rt-rr(1).rates(irx1));
       a(2) = -y(Cr48)*y(He4) * rr(1).rates(ircrap)*rr(index_rate).rates(irx1);
       a(3) = -y(Fe52) * rr(1).rates(irx1) * rr(index_rate).rates(irfegp);
       a(4) = -y(Fe52) * rr(index_rate).rates(irx1) * rr(1).rates(irfegp);
       a(5) = -y(Fe52)*y(He4) * rr(index_rate).rates(irfeap)*(1.0e0_rt-rr(1).rates(iry1));
       a(6) =  y(Fe52)*y(He4) * rr(1).rates(irfeap)*rr(index_rate).rates(iry1);
       a(7) =  y(Ni56) * rr(1).rates(iry1) * rr(index_rate).rates(irnigp);
       a(8) =  y(Ni56) * rr(index_rate).rates(iry1) * rr(1).rates(irnigp);

       dydt(Fe52) = dydt(Fe52) + esum8(a);

    }


    // ni56 reactions
    a(1) =  y(Fe52) * y(He4) * rr(index_rate).rates(irfeag);
    a(2) = -y(Ni56) * rr(index_rate).rates(irniga);

    dydt(Ni56) = dydt(Ni56) + a(1) + a(2);

    if (!deriva) {

       a(1) =  y(Fe52) * y(He4) * rr(index_rate).rates(irfeap)*(1.0e0_rt-rr(index_rate).rates(iry1));
       a(2) = -y(Ni56) * rr(index_rate).rates(iry1) * rr(index_rate).rates(irnigp);

       dydt(Ni56) = dydt(Ni56) + a(1) + a(2);

    } else {

       a(1) =  y(Fe52)*y(He4) * rr(index_rate).rates(irfeap)*(1.0e0_rt-rr(1).rates(iry1));
       a(2) = -y(Fe52)*y(He4) * rr(1).rates(irfeap)*rr(index_rate).rates(iry1);
       a(3) = -y(Ni56) * rr(1).rates(iry1) * rr(index_rate).rates(irnigp);
       a(4) = -y(Ni56) * rr(index_rate).rates(iry1) * rr(1).rates(irnigp);

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

    Array1D<rate_t, 1, Rates::NumGroups> rr;

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

#ifdef STRANG
    // Append the temperature equation

    ydot(net_itemp) = temperature_rhs(state, ydot(net_ienuc));
#endif
}


// Analytical Jacobian
template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_jac(burn_t& state, MatrixType& jac)
{

    Array1D<rate_t, 1, Rates::NumGroups> rr;

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
