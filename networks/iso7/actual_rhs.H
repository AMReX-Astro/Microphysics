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
#include <ArrayUtilities.H>

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

    extern AMREX_GPU_MANAGED Array2D<Real, 1, Rates::NumRates, 1, nrattab> rattab;
    extern AMREX_GPU_MANAGED Array2D<Real, 1, Rates::NumRates, 1, nrattab> drattabdt;
    extern AMREX_GPU_MANAGED Array1D<Real, 1, nrattab> ttab;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void iso7tab(const Real btemp, const Real bden,
             Array1D<Real, 1, Rates::NumRates>& rate,
             Array1D<Real, 1, Rates::NumRates>& dratedt)
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
        dtab(ircag)  = bden;
        dtab(iroga)  = 1.0e0_rt;
        dtab(ir3a)   = bden*bden;
        dtab(irg3a)  = 1.0e0_rt;
        dtab(ir1212) = bden;
        dtab(ir1216) = bden;
        dtab(ir1616) = bden;
        dtab(iroag)  = bden;
        dtab(irnega) = 1.0e0_rt;
        dtab(irneag) = bden;
        dtab(irmgga) = 1.0e0_rt;
        dtab(irmgag) = bden;
        dtab(irsiga) = 1.0e0_rt;
        dtab(ircaag) = bden;
        dtab(irtiga) = 1.0e0_rt;
        dtab(irsi2ni) = 0.0e0_rt;
        dtab(irni2si) = 0.0e0_rt;
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

       rate(j) = (  alfa * rattab(j,iat  )
                  + beta * rattab(j,iat+1)
                  + gama * rattab(j,iat+2)
                  + delt * rattab(j,iat+3) ) * dtab(j);

       dratedt(j) = (  alfa * drattabdt(j,iat  )
                     + beta * drattabdt(j,iat+1)
                     + gama * drattabdt(j,iat+2)
                     + delt * drattabdt(j,iat+3) ) * dtab(j);

    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void iso7rat(const Real btemp, const Real bden,
             Array1D<Real, 1, Rates::NumRates>& rate,
             Array1D<Real, 1, Rates::NumRates>& dratedt)
{
    using namespace Rates;

    // this routine generates unscreened
    // nuclear reaction rates for the iso7 network.

    Real rrate,drratedt;

    for (int i = 1; i <= Rates::NumRates; ++i) {
       rate(i)    = 0.0_rt;
       dratedt(i) = 0.0_rt;
    }

    if (btemp < 1.0e6_rt) return;


    // get the temperature factors
    tf_t tf = get_tfactors(btemp);

    // Determine which c12(a,g)o16 rate to use
    if (use_c12ag_deboer17) {
        // deboer + 2017 c12(a,g)o16 rate
        rate_c12ag_deboer17(tf,bden,
                            rate(ircag),dratedt(ircag),
                            rate(iroga),dratedt(iroga));
    } else {
        // 1.7 times cf88 c12(a,g)o16 rate
        rate_c12ag(tf,bden,
                   rate(ircag),dratedt(ircag),
                   rate(iroga),dratedt(iroga));
    }

    // triple alpha to c12
    rate_triplealf(tf,bden,
                   rate(ir3a),dratedt(ir3a),
                   rate(irg3a),dratedt(irg3a));

    // c12 + c12
    rate_c12c12(tf,bden,
                rate(ir1212),dratedt(ir1212),
                rrate,drratedt);

    // c12 + o16
    rate_c12o16(tf,bden,
                rate(ir1216),dratedt(ir1216),
                rrate,drratedt);

    // 16o + 16o
    rate_o16o16(tf,bden,
                rate(ir1616),dratedt(ir1616),
                rrate,drratedt);

    // o16(a,g)ne20
    rate_o16ag(tf,bden,
               rate(iroag),dratedt(iroag),
               rate(irnega),dratedt(irnega));

    // ne20(a,g)mg24
    rate_ne20ag(tf,bden,
                rate(irneag),dratedt(irneag),
                rate(irmgga),dratedt(irmgga));

    // mg24(a,g)si28
    rate_mg24ag(tf,bden,
                rate(irmgag),dratedt(irmgag),
                rate(irsiga),dratedt(irsiga));

    // ca40(a,g)ti44
    rate_ca40ag(tf,bden,
                rate(ircaag),dratedt(ircaag),
                rate(irtiga),dratedt(irtiga));

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void set_iso7rat()
{
    using namespace RateTable;

    Real btemp;
    Real bden = 1.0e0_rt;
    Array1D<Real, 1, Rates::NumRates> rate;
    Array1D<Real, 1, Rates::NumRates> dratedt;

    for (int i = 1; i <= tab_imax; ++i) {

       btemp = tab_tlo + static_cast<Real>(i-1) * tab_tstp;
       btemp = std::pow(10.0e0_rt, btemp);

       iso7rat(btemp, bden, rate, dratedt);

       ttab(i) = btemp;

       for (int j = 1; j <= Rates::NumRates; ++j) {

          rattab(j,i)    = rate(j);
          drattabdt(j,i) = dratedt(j);

       }
    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void screen_iso7(const Real btemp, const Real bden,
                 Array1D<Real, 1, NumSpec> const& y,
                 Array1D<Real, 1, Rates::NumRates>& rate,
                 Array1D<Real, 1, Rates::NumRates>& dratedt,
                 Array1D<Real, Rates::irsi2ni, Rates::irni2si>& dratedy1,
                 Array1D<Real, Rates::irsi2ni, Rates::irni2si>& dratedy2)
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

    // initialize
    for (int i = Rates::irsi2ni; i <= Rates::irni2si; ++i) {
        dratedy1(i) = 0.0_rt;
        dratedy2(i) = 0.0_rt;
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
            zion[He4-1], aion[He4-1], 4.0_rt, 8.0_rt,
            sc2a,sc2adt,sc2add);

    sc3a   = sc1a * sc2a;
    sc3adt = sc1adt*sc2a + sc1a*sc2adt;

    dratedt(ir3a) = dratedt(ir3a) * sc3a + rate(ir3a) * sc3adt;
    rate(ir3a)    = rate(ir3a) * sc3a;

    // c12 to o16
    jscr++;
    screen5(pstate,jscr,
            zion[C12-1], aion[C12-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);

    dratedt(ircag)  = dratedt(ircag) * sc1a + rate(ircag) * sc1adt;
    rate(ircag)     = rate(ircag) * sc1a;

    // c12 + c12
    jscr++;
    screen5(pstate,jscr,
            zion[C12-1], aion[C12-1], zion[C12-1], aion[C12-1],
            sc1a,sc1adt,sc1add);

    dratedt(ir1212) = dratedt(ir1212) * sc1a + rate(ir1212) * sc1adt;
    rate(ir1212)    = rate(ir1212) * sc1a;

    // c12 + o16
    jscr++;
    screen5(pstate,jscr,
            zion[C12-1], aion[C12-1], zion[O16-1], aion[O16-1],
            sc1a,sc1adt,sc1add);

    dratedt(ir1216) = dratedt(ir1216) * sc1a + rate(ir1216) * sc1adt;
    rate(ir1216)    = rate(ir1216) * sc1a;

    // o16 + o16
    jscr++;
    screen5(pstate,jscr,
            zion[O16-1], aion[O16-1], zion[O16-1], aion[O16-1],
            sc1a,sc1adt,sc1add);

    dratedt(ir1616) = dratedt(ir1616) * sc1a + rate(ir1616) * sc1adt;
    rate(ir1616)    = rate(ir1616) * sc1a;

    // o16 to ne20
    jscr++;
    screen5(pstate,jscr,
            zion[O16-1], aion[O16-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);

    dratedt(iroag) = dratedt(iroag) * sc1a + rate(iroag) * sc1adt;
    rate(iroag)    = rate(iroag) * sc1a;

    // ne20 to mg24
    jscr++;
    screen5(pstate,jscr,
            zion[Ne20-1], aion[Ne20-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);

    dratedt(irneag) = dratedt(irneag) * sc1a + rate(irneag) * sc1adt;
    rate(irneag)    = rate(irneag) * sc1a;

    // mg24 to si28
    jscr++;
    screen5(pstate,jscr,
            zion[Mg24-1], aion[Mg24-1], zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);

    dratedt(irmgag) = dratedt(irmgag) * sc1a + rate(irmgag) * sc1adt;
    rate(irmgag)    = rate(irmgag) * sc1a;

    // ca40 to ti44
    jscr++;
    screen5(pstate,jscr,
            20.0_rt, 40.0_rt, zion[He4-1], aion[He4-1],
            sc1a,sc1adt,sc1add);

    dratedt(ircaag) = dratedt(ircaag) * sc1a + rate(ircaag) * sc1adt;
    rate(ircaag)    = rate(ircaag) * sc1a;

    /*
    the publication, timmes, woosley & hoffman apjs, 129, 377
    has a typo on page 393, where its says "y(C12)+y(O16) .gt. 0.004"
    it should be less than or equal to, since the idea is this piece
    gets activated during silicon buring, after all the c + o from
    oxygen burning is gone.
    */

    if (tf.t9 > 2.5_rt && y(C12) + y(O16) <= 4.0e-3_rt) {

       t992  = tf.t972 * tf.t9;
       t9i92 = 1.0e0_rt/t992;

       yeff_ca40   = t9i92 * std::exp(239.42_rt*tf.t9i - 74.741_rt);
       yeff_ca40dt = -yeff_ca40*(239.42_rt*tf.t9i2 + 4.5e0_rt*tf.t9i);

       yeff_ti44   = t992  * std::exp(-274.12_rt*tf.t9i + 74.914_rt);
       yeff_ti44dt = yeff_ti44*(274.12_rt*tf.t9i2 + 4.5e0_rt*tf.t9i);

       denom     = std::pow(bden * y(He4), 3.0e0_rt);

       rate(irsi2ni)     = yeff_ca40*denom*rate(ircaag)*y(Si28);
       dratedy1(irsi2ni) = 3.0e0_rt * rate(irsi2ni)/y(He4);
       dratedy2(irsi2ni) = yeff_ca40*denom*rate(ircaag);
       dratedt(irsi2ni)  = (yeff_ca40dt*rate(ircaag)
                          + yeff_ca40*dratedt(ircaag))*denom*y(Si28)*1.0e-9_rt;

       if (denom != 0.0_rt) {

          zz     = 1.0e0_rt/denom;
          rate(irni2si) = amrex::min(1.0e10_rt,yeff_ti44*rate(irtiga)*zz);

          if (rate(irni2si) == 1.0e10_rt) {
             dratedy1(irni2si) = 0.0e0_rt;
             dratedt(irni2si)  = 0.0e0_rt;
          } else {
             dratedy1(irni2si) = -3.0e0_rt * rate(irni2si)/y(He4);
             dratedt(irni2si)  = (yeff_ti44dt*rate(irtiga)
                                + yeff_ti44*dratedt(irtiga))*zz*1.0e-9_rt;
          }

       }

    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_rates(burn_t const& state, rate_t& rr)
{
    Real rho, temp;
    Array1D<Real, 1, NumSpec> y;
    Array1D<Real, 1, Rates::NumRates> rate, dratedt;
    Array1D<Real, Rates::irsi2ni, Rates::irni2si> dratedy1;
    Array1D<Real, Rates::irsi2ni, Rates::irni2si> dratedy2;

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
    screen_iso7(temp, rho, y, rate, dratedt, dratedy1, dratedy2);

    // Save the rate data, for the Jacobian later if we need it.
    for (int i = 1; i <= Rates::NumRates; ++i) {
        rr.rates(1,i) = rate(i);
        rr.rates(2,i) = dratedt(i);
    }

    for (int i = Rates::irsi2ni; i <= Rates::irni2si; ++i) {
        rr.rates(3,i) = dratedy1(i);
        rr.rates(4,i) = dratedy2(i);
    }

    rr.T_eval = temp;
}


template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void dfdy_isotopes_iso7(Array1D<Real, 1, NumSpec> const& y,
                        burn_t const& state, rate_t const& rr,
                        MatrixType& jac)
{

    // this routine sets up the dense iso7 jacobian for the isotopes

    using namespace Species;
    using namespace Rates;

    // set up the jacobian
    // 4he jacobian elements
    // d(he4)/d(he4)
    jac(He4,He4) = esum(-1.5e0_rt * y(He4) * y(He4) * rr.rates(1,ir3a),
                        -y(C12) * rr.rates(1,ircag),
                        -y(O16) * rr.rates(1,iroag),
                        -y(Ne20) * rr.rates(1,irneag),
                        -y(Mg24) * rr.rates(1,irmgag),
                        -7.0e0_rt * rr.rates(1,irsi2ni),
                        -7.0e0_rt * rr.rates(3,irsi2ni) * y(He4),
                        7.0e0_rt * rr.rates(3,irni2si) * y(Ni56));

    // d(he4)/d(c12)
    jac(He4,C12) = esum(3.0e0_rt * rr.rates(1,irg3a),
                        -y(He4) * rr.rates(1,ircag),
                        y(C12) * rr.rates(1,ir1212),
                        0.5e0_rt * y(O16) * rr.rates(1,ir1216));


    // d(he4)/d(o16)
    jac(He4,O16) = esum(rr.rates(1,iroga),
                        0.5e0_rt * y(C12) * rr.rates(1,ir1216),
                        y(O16) * rr.rates(1,ir1616),
                        -y(He4) * rr.rates(1,iroag));


    // d(c12)/d(c12)
    jac(C12,C12) = esum(-rr.rates(1,irg3a),
                        -y(He4) * rr.rates(1,ircag),
                        -2.0e0_rt * y(C12) * rr.rates(1,ir1212),
                        -y(O16) * rr.rates(1,ir1216));


    // d(o16)/d(o16)
    jac(O16,O16) = esum(-rr.rates(1,iroga),
                        -y(C12) * rr.rates(1,ir1216),
                        -2.0e0_rt * y(O16) * rr.rates(1,ir1616),
                        -y(He4) * rr.rates(1,iroag));


    // 28si jacobian elements
    // d(si28)/d(he4)
    jac(Si28,He4) = esum(y(Mg24) * rr.rates(1,irmgag),
                            -rr.rates(1,irsi2ni),
                            -rr.rates(3,irsi2ni) * y(He4),
                            rr.rates(3,irni2si) * y(Ni56));

    // ni56 jacobian elements
    // d(ni56)/d(he4)
    jac(Ni56,He4) = esum(rr.rates(1,irsi2ni),
                            rr.rates(3,irsi2ni) * y(He4),
                            -rr.rates(3,irni2si) * y(Ni56));

    // d(he4)/d(ne20)
    jac(He4,Ne20) = rr.rates(1,irnega) +
                    -y(He4) * rr.rates(1,irneag);


    // d(he4)/d(mg24)
    jac(He4,Mg24) = rr.rates(1,irmgga) +
                    -y(He4) * rr.rates(1,irmgag);


    // d(he4)/d(si28)
    jac(He4,Si28) = rr.rates(1,irsiga) +
                    -7.0e0_rt * rr.rates(4,irsi2ni) * y(He4);


    // 12c jacobian elements
    // d(c12)/d(he4)
    jac(C12,He4) = 0.5e0_rt * y(He4) * y(He4) * rr.rates(1,ir3a) +
                    -y(C12) * rr.rates(1,ircag);


    // d(c12)/d(o16)
    jac(C12,O16) = rr.rates(1,iroga) +
                    -y(C12) * rr.rates(1,ir1216);


    // 16o jacobian elements
    // d(o16)/d(he4)
    jac(O16,He4) = y(C12) * rr.rates(1,ircag) +
                    -y(O16) * rr.rates(1,iroag);


    // d(o16)/d(c12)
    jac(O16,C12) = y(He4) * rr.rates(1,ircag) +
                    -y(O16) * rr.rates(1,ir1216);


    // 24mg jacobian elements
    // d(mg24)/d(he4)
    jac(Mg24,He4) = y(Ne20) * rr.rates(1,irneag) +
                    -y(Mg24) * rr.rates(1,irmgag);


    // d(mg24)/d(mg24)
    jac(Mg24,Mg24) = -rr.rates(1,irmgga) +
                        -y(He4) * rr.rates(1,irmgag);


    // d(si28)/d(o16)
    jac(Si28,O16) = y(O16) * rr.rates(1,ir1616) +
                    0.5e0_rt * y(C12) * rr.rates(1,ir1216);


    // d(si28)/d(si28)
    jac(Si28,Si28) = -rr.rates(1,irsiga) +
                        -rr.rates(4,irsi2ni) * y(He4);

    // d(he4)/d(ni56)
    jac(He4,Ni56) = 7.0e0_rt * rr.rates(1,irni2si);

    // d(o16)/d(ne20)
    jac(O16,Ne20) = rr.rates(1,irnega);

    // 20ne jacobian elements
    // d(ne20)/d(he4)
    jac(Ne20,He4) = y(O16) * rr.rates(1,iroag) - y(Ne20) * rr.rates(1,irneag);

    // d(ne20)/d(c12)
    jac(Ne20,C12) = y(C12) * rr.rates(1,ir1212);

    // d(ne20)/d(o16)
    jac(Ne20,O16) = y(He4) * rr.rates(1,iroag);

    // d(ne20)/d(ne20)
    jac(Ne20,Ne20) = -rr.rates(1,irnega) - y(He4) * rr.rates(1,irneag);

    // d(ne20)/d(mg24)
    jac(Ne20,Mg24) = rr.rates(1,irmgga);

    // d(mg24)/d(c12)
    jac(Mg24,C12) = 0.5e0_rt * y(O16) * rr.rates(1,ir1216);

    // d(mg24)/d(o16)
    jac(Mg24,O16) = 0.5e0_rt * y(C12) * rr.rates(1,ir1216);

    // d(mg24)/d(ne20)
    jac(Mg24,Ne20) = y(He4) * rr.rates(1,irneag);

    // d(mg24)/d(si28)
    jac(Mg24,Si28) = rr.rates(1,irsiga);

    // d(si28)/d(c12)
    jac(Si28,C12) = 0.5e0_rt * y(O16) * rr.rates(1,ir1216);

    // d(si28)/d(mg24)
    jac(Si28,Mg24) = y(He4) * rr.rates(1,irmgag);

    // d(si28)/d(ni56)
    jac(Si28,Ni56) = rr.rates(1,irni2si);

    // d(ni56)/d(si28)
    jac(Ni56,Si28) = rr.rates(4,irsi2ni) * y(He4);

    // d(ni56)/d(ni56)
    jac(Ni56,Ni56) = -rr.rates(1,irni2si);

}


template<class T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real ener_gener_rate(T const& dydt)
{

    using namespace iso7;

    // Computes the instantaneous energy generation rate

    Real Xdot = 0.0_rt;

    // Sum the mass fraction time derivatives
    for (int i = 1; i <= NumSpec; ++i) {
        Xdot += dydt(i) * mion(i);
    }

    // This is basically e = m c**2
    return Xdot * C::Legacy::enuc_conv2;

}


// Evaluates the right hand side of the iso7 ODEs
template<class dydt_T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rhs(Array1D<Real, 1, NumSpec> const& y, rate_t const& rr,
         dydt_T& dydt,
         bool deriva, bool for_jacobian_tderiv)
{
    using namespace Species;
    using namespace Rates;

    // deriva is used in forming the analytic Jacobian to get
    // the derivative wrt A

    const int index_rate = for_jacobian_tderiv ? 2 : 1;

    for (int i = 1; i <= neqs; ++i) {
        dydt(i) = 0.0_rt;
    }

    // set up the system of ode's :
    // 4he reactions
    dydt(He4) = esum(   3.0e0_rt * y(C12) * rr.rates(index_rate,irg3a),
                        -0.5e0_rt * y(He4) * y(He4) * y(He4) * rr.rates(index_rate,ir3a),
                        y(O16) * rr.rates(index_rate,iroga),
                        -y(C12) * y(He4) * rr.rates(index_rate,ircag),
                        0.5e0_rt * y(C12) * y(C12) * rr.rates(index_rate,ir1212),
                        0.5e0_rt * y(C12) * y(O16) * rr.rates(index_rate,ir1216),
                        0.5e0_rt * y(O16) * y(O16) * rr.rates(index_rate,ir1616),
                        -y(O16) * y(He4) * rr.rates(index_rate,iroag),
                        y(Ne20) * rr.rates(index_rate,irnega),
                        y(Mg24) * rr.rates(index_rate,irmgga),
                        -y(Ne20) * y(He4) * rr.rates(index_rate,irneag),
                        y(Si28) * rr.rates(index_rate,irsiga),
                        -y(Mg24) * y(He4) * rr.rates(index_rate,irmgag),
                        -7.0e0_rt * rr.rates(index_rate,irsi2ni) * y(He4),
                        7.0e0_rt * rr.rates(index_rate,irni2si) * y(Ni56)   );

    // 12c reactions
    dydt(C12) = esum(   y(He4) * y(He4) * y(He4) * rr.rates(index_rate,ir3a) / 6.0_rt,
                        -y(C12) * rr.rates(index_rate,irg3a),
                        y(O16) * rr.rates(index_rate,iroga),
                        -y(C12) * y(He4) * rr.rates(index_rate,ircag),
                        -y(C12) * y(C12) * rr.rates(index_rate,ir1212),
                        -y(C12) * y(O16) * rr.rates(index_rate,ir1216)   );


    // 16o reactions
    dydt(O16) = esum(   -y(O16) * rr.rates(index_rate,iroga);
                        y(C12) * y(He4) * rr.rates(index_rate,ircag);
                        -y(C12) * y(O16) * rr.rates(index_rate,ir1216);
                        -y(O16) * y(O16) * rr.rates(index_rate,ir1616);
                        -y(O16) * y(He4) * rr.rates(index_rate,iroag);
                        y(Ne20) * rr.rates(index_rate,irnega)   );


    // 20ne reactions
    dydt(Ne20) = esum(  0.5e0_rt * y(C12) * y(C12) * rr.rates(index_rate,ir1212),
                        y(O16) * y(He4) * rr.rates(index_rate,iroag),
                        -y(Ne20) * rr.rates(index_rate,irnega),
                        y(Mg24) * rr.rates(index_rate,irmgga),
                        -y(Ne20) * y(He4) * rr.rates(index_rate,irneag)  );


    // 24mg reactions
    dydt(Mg24) = esum(  0.5e0_rt * y(C12) * y(O16) * rr.rates(index_rate,ir1216),
                        -y(Mg24) * rr.rates(index_rate,irmgga),
                        y(Ne20) * y(He4) * rr.rates(index_rate,irneag),
                        y(Si28) * rr.rates(index_rate,irsiga),
                        -y(Mg24) * y(He4) * rr.rates(index_rate,irmgag)  );


    // 28si reactions
    dydt(Si28) = esum(  0.5e0_rt * y(C12) * y(O16) * rr.rates(index_rate,ir1216),
                        0.5e0_rt * y(O16) * y(O16) * rr.rates(index_rate,ir1616),
                        -y(Si28) * rr.rates(index_rate,irsiga),
                        y(Mg24) * y(He4) * rr.rates(index_rate,irmgag),
                        -rr.rates(index_rate,irsi2ni) * y(He4),
                        rr.rates(index_rate,irni2si) * y(Ni56)  );


    // ni56 reactions
    dydt(Ni56) = rr.rates(index_rate,irsi2ni) * y(He4) +
                    -rr.rates(index_rate,irni2si) * y(Ni56);

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_rhs(burn_t& state, Array1D<Real, 1, neqs>& ydot)
{

    /*
     This routine sets up the system of ODE's for the iso7
     nuclear reactions.  This is an alpha chain + heavy ion network
     with (a,p)(p,g) links up to silicon, as well as a Si <-> Ni link.

     Isotopes: he4,  c12,  o16,  ne20, mg24, si28, ni56
    */

    rate_t rr;

    Real sneut, dsneutdt, dsneutdd, snuda, snudz;
    Real enuc;

    Real rho, temp, abar, zbar;
    Array1D<Real, 1, NumSpec> y;

    // Get the data from the state

    rho  = state.rho;
    temp = state.T;
    abar = state.abar;
    zbar = state.zbar;

    for (int i = 1; i <= NumSpec; ++i) {
        y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    get_rates(state, rr);

    // Call the RHS to actually get dydt.

    bool deriva = false;
    bool for_jacobian_tderiv = false;
    rhs(y, rr, ydot, deriva, for_jacobian_tderiv);

    // Instantaneous energy generation rate -- this needs molar fractions

    enuc = ener_gener_rate(ydot);

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

    jac.zero();

    get_rates(state, rr);

    // Get the data from the state

    rho  = state.rho;
    temp = state.T;
    abar = state.abar;
    zbar = state.zbar;

    for (int i = 1; i <= NumSpec; ++i)
        y(i) = state.xn[i-1] * aion_inv[i-1];

    // Species Jacobian elements with respect to other species

    dfdy_isotopes_iso7(y, state, rr, jac);

    // Energy generation rate Jacobian elements with respect to species

    for (int j = 1; j <= NumSpec; ++j) {
        auto jac_slice_2 = [&](int i) -> Real { return jac.get(i, j); };
        jac(net_ienuc,j) = ener_gener_rate(jac_slice_2);
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

    jac(net_ienuc,net_itemp) = ener_gener_rate(yderivs);

    jac.add(net_ienuc, net_itemp, -dsneutdt);

    // Temperature Jacobian elements

    temperature_jac(state, jac);

}


// RHS & Analytical Jacobian
template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_rhs_jac(burn_t& state, Array1D<Real, 1, neqs>& ydot, MatrixType& jac)
{

    // Initialize jac to 0

    jac.zero();

    // Initialize ydot to 0

    for (int i = 1; i <= neqs; ++i) {
       ydot(i) = 0.0_rt;
    }

    // Get the data from the state

    const Real rho  = state.rho;
    const Real temp = state.T;
    const Real abar = state.abar;
    const Real zbar = state.zbar;

    Array1D<Real, 1, NumSpec> y;
    for (int i = 1; i <= NumSpec; ++i) {
        y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    rate_t rr;
    evaluate_rates(state, rr);

    // Get the species part of the RHS
    {
        const bool deriva = false;
        const bool for_jacobian_tderiv = false;
        rhs(y, rr, ydot, deriva, for_jacobian_tderiv);
    }

    // Species Jacobian elements with respect to other species

    dfdy_isotopes_iso7(y, state, rr, jac);

    // Energy generate rate Jacobian element with respect to temperature

    auto jac_species_temp = [&](int i) -> Real& {
        return jac(i, net_itemp);
    };

    {
        const bool deriva = true;
        const bool for_jacobian_tderiv = true;
        rhs(y, rr, jac_species_temp, deriva, for_jacobian_tderiv);
    }

    // Account for the thermal neutrino losses

    Real b1, sneut, dsneutdt, dsneutdd, snuda, snudz;

    sneut5(temp, rho, abar, zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);

    // Energy-temperature element of the Jacobian

    jac(net_ienuc,net_itemp) = ener_gener_rate(jac_species_temp) - dsneutdt;

    // Instantaneous energy generation rate (erg/g/s) -- this needs molar fractions

    ydot(net_ienuc) = ener_gener_rate(ydot) - sneut;

    // Energy generation rate Jacobian elements with respect to species

    for (int j = 1; j <= NumSpec; ++j) {
        auto jac_slice_2 = [&](int i) -> Real { return jac.get(i, j); };

        // thermal neutrino contribution is b1
        b1 = (-abar * abar * snuda + (zion[j-1] - zbar) * abar * snudz);

        jac(net_ienuc,j) = ener_gener_rate(jac_slice_2) - b1;
    }

    // Set the temperature RHS

    temperature_rhs(state, ydot);

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
    add_screening_factor(jscr++, zion[He4-1], aion[He4-1],  4.0e0_rt,  8.0e0_rt);
    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[He4-1], aion[He4-1]);
    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[C12-1], aion[C12-1]);
    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[O16-1], aion[O16-1]);
    add_screening_factor(jscr++, zion[O16-1], aion[O16-1], zion[O16-1], aion[O16-1]);
    add_screening_factor(jscr++, zion[O16-1], aion[O16-1], zion[He4-1], aion[He4-1]);
    add_screening_factor(jscr++, zion[Ne20-1], aion[Ne20-1], zion[He4-1], aion[He4-1]);
    add_screening_factor(jscr++, zion[Mg24-1], aion[Mg24-1], zion[He4-1], aion[He4-1]);
    add_screening_factor(jscr++,  20.0e0_rt,  40.0e0_rt, zion[He4-1], aion[He4-1]);

}

#endif
