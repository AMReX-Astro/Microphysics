#ifndef _actual_rhs_H_
#define _actual_rhs_H_

#include <AMReX_CONSTANTS.H>
#include <AMReX_REAL.H>

#include <extern_parameters.H>
#include <tfactors.H>
#include <actual_network.H>
#include <aprox_rates.H>

using namespace amrex;

namespace
{
    // does anything actually use this?
    constexpr Real c54 = 56.0e0_rt/54.0e0_rt;
}

namespace RateTable
{
    constexpr Real tab_tlo = 6.0e0_rt;
    constexpr Real tab_thi = 10.0e0_rt;
    constexpr int tab_per_decade = 2000;
    constexpr int nrattab = static_cast<int>(tab_thi - tab_tlo) * tab_per_decade + 1;
    constexpr int tab_imax = static_cast<int>(tab_thi - tab_tlo) * tab_per_decade + 1;
    constexpr Real tab_tstp = (tab_thi - tab_tlo) / static_cast<Real>(tab_imax - 1);

    extern AMREX_GPU_MANAGED Array2D<Real, 0, Rates::NRates-1, 0, nrattab-1> rattab;
    extern AMREX_GPU_MANAGED Array2D<Real, 0, Rates::NRates-1, 0, nrattab-1> drattabdt;
    extern AMREX_GPU_MANAGED Array1D<Real, 0, nrattab-1> ttab;
}

AMREX_INLINE
void actual_rhs_init()
{
    //use screening_module, only: screening_init

    rates_init()

    // set_up_screening_factors()

    // screening_init()

    if (use_tables)
    {
        Print() << "\nInitializing iso7 rate table\n";
        set_iso7rat();
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_rates(burn_t const& state, rate_t& rr)
{
    amrex::Real rho, temp;
    amrex::Array1D<amrex::Real, 0, Species::NSpec-1> y;
    amrex::Array1D<amrex::Real, 0, Rates::NRates-1> rate, dratedt;
    amrex::Array1D<amrex::Real, Rates::irsi2ni, Rates::irni2si> dratedy1;
    amrex::Array1D<amrex::Real, Rates::irsi2ni, Rates::irni2si> dratedy2;

    // Get the data from the state
    rho  = state.rho;
    temp = state.T;

    for (int i = 0; i < Species::NSpec;; ++i)
        y(i) = state.xn(i) * aion_inv[i];

    // Get the raw reaction rates
    if (use_tables)
        iso7tab(temp, rho, rate, dratedt);
    else
        iso7rat(temp, rho, rate, dratedt);

    // Do the screening here because the corrections depend on the composition
    screen_iso7(temp, rho, y, rate, dratedt, dratedy1, dratedy2);

    // Save the rate data, for the Jacobian later if we need it.
    for (int i = 0; i < Rates::NRates; ++i) {
        rr.rates(1,i) = rate(i);
        rr.rates(2,i) = dratedt(i);
    }

    for (int i = Rates::irsi2ni, i <= Rates::irni2si; ++i) {
        rr.rates(3,i) = dratedy1(i);
        rr.rates(4,i) = dratedy2(i);
    }

    rr.T_eval = temp;
}

void iso7tab(const Real btemp, const Real bden,
             Array1D<Real, 0, Rates::NRates-1>& rate,
             Array1D<Real, 0, Rates::NRates-1>& dratedt)
{
    constexpr int mp = 4;

    int j, iat;
    Real x, x1, x2, x3, x4;
    Real a, b, c, d, e, f, g, h, p, q;
    Real alfa, beta, gama, delt;
    Array1D<Real, 0, Rates::NRates-1> dtab;

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
    iat = amrex::max(1, amrex::min(iat - 1, tab_imax - mp + 1)) - 1;

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
    for (int j = 0; j < Rates::NRates; ++j) {

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


void set_iso7rat()
{

    Real btemp;
    Real bden = 1.0e0_rt;
    Array1D<Real, 0, Rates::NRates-1> rate;
    Array1D<Real, 0, Rates::NRates-1> dratedt;

    for (int i = 0; i < tab_imax; ++i) {

       btemp = tab_tlo + static_cast<Real>(i) * tab_tstp;
       btemp = std::pow(10.0e0_rt, btemp);

       iso7rat(btemp, bden, rate, dratedt);

       ttab(i) = btemp;

       for (int j = 0; j < Rates::NRates; ++j) {

          rattab(j,i)    = rate(j);
          drattabdt(j,i) = dratedt(j);

       }
    }
}


void actual_rhs(burn_t const& state, Array1D<Real, 0, Burn::neqs-1>& ydot)
{
    // need to port sneut5
    /* use sneut_module, only: sneut5 */
    // need to port temperature_rhs
    /* use temperature_integration_module, only: temperature_rhs */

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
    Array1D<Real, 0, Species::NSpec-1> y;

    // Get the data from the state

    rho  = state.rho;
    temp = state.T;
    abar = state.abar;
    zbar = state.zbar;

    for (int i = 0; i < Species::NSpec; ++i) {
        y(i) = state.xn(i) * aion_inv(i);
    }

    get_rates(state, rr);

    // Call the RHS to actually get dydt.

    bool deriva = false;
    bool for_jacobian_tderiv = false;
    rhs(y, rr, ydot, deriva, for_jacobian_tderiv);

    // Instantaneous energy generation rate -- this needs molar fractions

    ener_gener_rate(ydot, enuc);

    // Get the neutrino losses

    sneut5(temp, rho, abar, zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);

    // Append the energy equation (this is erg/g/s)

    ydot(Burn::ienuc) = enuc - sneut;

    // Append the temperature equation

    temperature_rhs(state, ydot);
}


  // Analytical Jacobian

void actual_jac(burn_t const& state, Array2D<Real, 0, Burn::njrows-1, 0, Burn::njcols-1>& jac)
{

    // use sneut_module, only: sneut5
    // use temperature_integration_module, only: temperature_jac

    rate_t rr;

    bool deriva;

    Real b1, sneut, dsneutdt, dsneutdd, snuda, snudz;

    int j;

    Real rho, temp, abar, zbar;
    Array1D<Real, 0, Species::NSpec-1> y;
    Array1D<Real, 0, Species::NSpec-1> yderivs;

    for (int j = 0; j < Burn::njcols; ++j) {
        for (int i = 0; i < Burn::njrows; ++i) {
            jac(i,j) = 0.0e0_rt;
        }
    }

    get_rates(state, rr);

    // Get the data from the state

    rho  = state.rho;
    temp = state.T;
    abar = state.abar;
    zbar = state.zbar;

    for (int i = 0; i < Species::NSpec; ++i)
        y(i) = state.xn(i) * aion_inv(i);

    // Species Jacobian elements with respect to other species

    dfdy_isotopes_iso7(y, state, rr, jac);

    // Energy generation rate Jacobian elements with respect to species

    for (int j = 0; j < Species::NSpec; ++j) {
        for (int i = 0; i < Species::NSpec; ++i) {
            ener_gener_rate(jac(i,j), jac(Burn::ienuc,j));
        }
    }

    // Account for the thermal neutrino losses

    sneut5(temp, rho, abar, zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);

    for (int j = 0; j < Species::NSpec; ++j) {
       b1 = (-abar * abar * snuda + (zion(j) - zbar) * abar * snudz);
       jac(Burn::ienuc,j) = jac(Burn::ienuc,j) - b1;
    }

    // Evaluate the Jacobian elements with respect to temperature by
    // calling the RHS using d(rate) / dT

    deriva = true;
    bool for_jacobian_tderiv = true;
    rhs(y, rr, yderivs, deriva, for_jacobian_tderiv);

    for (int i = 0; i < Species::NSpec; ++i) {
        jac(i,Burn::itemp) = yderivs(i);
    }

    for (int i = 0; i < Species::NSpec; ++i) {
        ener_gener_rate(jac(i,Burn::itemp), jac(Burn::ienuc,Burn::itemp));
    }

    jac(Burn::ienuc,Burn::itemp) = jac(Burn::ienuc,Burn::itemp) - dsneutdt;

    // Temperature Jacobian elements

    temperature_jac(state, jac);

}


// Evaluates the right hand side of the iso7 ODEs

void rhs(Array1D<Real, 0, Species::NSpec-1> const& y, rate_t const& rr,
         Array1D<Real, 0, Species::NSpec-1>& dydt,
         bool deriva, bool for_jacobian_tderiv)
{
    using namespace Species;
    using namespace Rates;

    // need to port esum functions ...
    /* use microphysics_math_module, only: esum5, esum6, esum15 ! function */

    // deriva is used in forming the analytic Jacobian to get
    // the derivative wrt A

    const int index_rate = for_jacobian_tderiv ? 2 : 1;
    Array1D<Real, 1, 15> a;

    for (int i = 0; i < Species::NSpec; ++i)
        dydt(i) = zero;

    // set up the system of ode's :
    // 4he reactions
    a(1)  =  3.0e0_rt * y(C12) * rr.rates(index_rate,irg3a);
    a(2)  = -0.5e0_rt * y(He4) * y(He4) * y(He4) * rr.rates(index_rate,ir3a);
    a(3)  =  y(O16) * rr.rates(index_rate,iroga);
    a(4)  = -y(C12) * y(He4) * rr.rates(index_rate,ircag);
    a(5)  =  0.5e0_rt * y(C12) * y(C12) * rr.rates(index_rate,ir1212);
    a(6)  =  0.5e0_rt * y(C12) * y(O16) * rr.rates(index_rate,ir1216);
    a(7)  =  0.5e0_rt * y(O16) * y(O16) * rr.rates(index_rate,ir1616);
    a(8)  = -y(O16) * y(He4) * rr.rates(index_rate,iroag);
    a(9)  =  y(Ne20) * rr.rates(index_rate,irnega);
    a(10) =  y(Mg24) * rr.rates(index_rate,irmgga);
    a(11) = -y(Ne20) * y(He4) * rr.rates(index_rate,irneag);
    a(12) =  y(Si28) * rr.rates(index_rate,irsiga);
    a(13) = -y(Mg24) * y(He4) * rr.rates(index_rate,irmgag);
    a(14) = -7.0e0_rt * rr.rates(index_rate,irsi2ni) * y(He4);
    a(15) =  7.0e0_rt * rr.rates(index_rate,irni2si) * y(Ni56);

    dydt(He4) = esum15(a);

    // 12c reactions
    a(1) =  sixth * y(He4) * y(He4) * y(He4) * rr.rates(index_rate,ir3a);
    a(2) = -y(C12) * rr.rates(index_rate,irg3a);
    a(3) =  y(O16) * rr.rates(index_rate,iroga);
    a(4) = -y(C12) * y(He4) * rr.rates(index_rate,ircag);
    a(5) = -y(C12) * y(C12) * rr.rates(index_rate,ir1212);
    a(6) = -y(C12) * y(O16) * rr.rates(index_rate,ir1216);

    dydt(C12) = esum6(a);

    // 16o reactions
    a(1) = -y(O16) * rr.rates(index_rate,iroga);
    a(2) =  y(C12) * y(He4) * rr.rates(index_rate,ircag);
    a(3) = -y(C12) * y(O16) * rr.rates(index_rate,ir1216);
    a(4) = -y(O16) * y(O16) * rr.rates(index_rate,ir1616);
    a(5) = -y(O16) * y(He4) * rr.rates(index_rate,iroag);
    a(6) =  y(Ne20) * rr.rates(index_rate,irnega);

    dydt(O16) = esum6(a);

    // 20ne reactions
    a(1) =  0.5e0_rt * y(C12) * y(C12) * rr.rates(index_rate,ir1212);
    a(2) =  y(O16) * y(He4) * rr.rates(index_rate,iroag);
    a(3) = -y(Ne20) * rr.rates(index_rate,irnega);
    a(4) =  y(Mg24) * rr.rates(index_rate,irmgga);
    a(5) = -y(Ne20) * y(He4) * rr.rates(index_rate,irneag);

    dydt(Ne20) = esum5(a);

    // 24mg reactions
    a(1) =  0.5e0_rt * y(C12) * y(O16) * rr.rates(index_rate,ir1216);
    a(2) = -y(Mg24) * rr.rates(index_rate,irmgga);
    a(3) =  y(Ne20) * y(He4) * rr.rates(index_rate,irneag);
    a(4) =  y(Si28) * rr.rates(index_rate,irsiga);
    a(5) = -y(Mg24) * y(He4) * rr.rates(index_rate,irmgag);

    dydt(Mg24) = esum5(a);

    // 28si reactions
    a(1) =  0.5e0_rt * y(C12) * y(O16) * rr.rates(index_rate,ir1216);
    a(2) =  0.5e0_rt * y(O16) * y(O16) * rr.rates(index_rate,ir1616);
    a(3) = -y(Si28) * rr.rates(index_rate,irsiga);
    a(4) =  y(Mg24) * y(He4) * rr.rates(index_rate,irmgag);
    a(5) = -rr.rates(index_rate,irsi2ni) * y(He4);
    a(6) =  rr.rates(index_rate,irni2si) * y(Ni56);

    dydt(Si28) = esum6(a);

    // ni56 reactions
    a(1) =  rr.rates(index_rate,irsi2ni) * y(He4);
    a(2) = -rr.rates(index_rate,irni2si) * y(Ni56);

    dydt(Ni56) = a(1) + a(2);

}


void iso7rat(const Real btemp, const Real bden,
             Array1D<Real, 0, Rates::NRates-1>& rate,
             Array1D<Real, 0, Rates::NRates-1>& dratedt)
{

    // this routine generates unscreened
    // nuclear reaction rates for the iso7 network.

    Real rrate,drratedt,drratedd,drratedd2;
    Real ff1,dff1dt,dff1dd,ff2,dff2dt,dff2dd,tot,dtotdt,dtotdd,invtot;
    tf_t tf;

    for (int i = 0; i < Rates::NRates; ++i) {
       rate(i)    = zero;
       dratedt(i) = zero;
    }

    if (btemp < 1.0e6_rt) return;


    // get the temperature factors
    get_tfactors(btemp, tf);

    // Determine which c12(a,g)o16 rate to use
    if (use_c12ag_deboer17) {
        // deboer + 2017 c12(a,g)o16 rate
        rate_c12ag_deboer17(tf,bden,
                            rate(ircag),dratedt(ircag),drratedd,
                            rate(iroga),dratedt(iroga),drratedd2);
    } else {
        // 1.7 times cf88 c12(a,g)o16 rate
        rate_c12ag(tf,bden,
                   rate(ircag),dratedt(ircag),drratedd,
                   rate(iroga),dratedt(iroga),drratedd2);
    }

    // triple alpha to c12
    rate_tripalf(tf,bden,
                 rate(ir3a),dratedt(ir3a),drratedd,
                 rate(irg3a),dratedt(irg3a),drratedd2);

    // c12 + c12
    rate_c12c12(tf,bden,
                rate(ir1212),dratedt(ir1212),drratedd2,
                rrate,drratedt,drratedd);

    // c12 + o16
    rate_c12o16(tf,bden,
                rate(ir1216),dratedt(ir1216),drratedd2,
                rrate,drratedt,drratedd);

    // 16o + 16o
    rate_o16o16(tf,bden,
                rate(ir1616),dratedt(ir1616),drratedd2,
                rrate,drratedt,drratedd);

    // o16(a,g)ne20
    rate_o16ag(tf,bden,
               rate(iroag),dratedt(iroag),drratedd,
               rate(irnega),dratedt(irnega),drratedd2);

    // ne20(a,g)mg24
    rate_ne20ag(tf,bden,
                rate(irneag),dratedt(irneag),drratedd,
                rate(irmgga),dratedt(irmgga),drratedd2);

    // mg24(a,g)si28
    rate_mg24ag(tf,bden,
                rate(irmgag),dratedt(irmgag),drratedd,
                rate(irsiga),dratedt(irsiga),drratedd2);

    // ca40(a,g)ti44
    rate_ca40ag(tf,bden,
                rate(ircaag),dratedt(ircaag),drratedd,
                rate(irtiga),dratedt(irtiga),drratedd2);

}

void screen_iso7(const Real btemp, const Real bden,
                 Array1D<Real, 0, Species::NSpec-1> const& y,
                 Array1D<Real, 0, Species::NRates-1>& rate,
                 Array1D<Real, 0, Species::NRates-1>& dratedt,
                 Array1D<Real, Rates::irsi2ni, Rates::irni2si>& dratedy1,
                 Array1D<Real, Rates::irsi2ni, Rates::irni2si>& dratedy2)
{
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

    plasma_state pstate;
    tf_t tf;

    // initialize
    for (int i = Rates::irsi2ni; i <= Rates::irni2si; ++i) {
        dratedy1(i) = zero;
        dratedy2(i) = zero;
    }

    // get the temperature factors
    get_tfactors(btemp, tf);

    // Set up the state data, which is the same for all screening factors.

    fill_plasma_state(pstate, btemp, bden, y);

    // first the always fun triple alpha and its inverse
    jscr = 1;
    screen5(pstate,jscr,sc1a,sc1adt,sc1add);

    jscr = jscr + 1;
    screen5(pstate,jscr,sc2a,sc2adt,sc2add);

    sc3a   = sc1a * sc2a;
    sc3adt = sc1adt*sc2a + sc1a*sc2adt;

    dratedt(ir3a) = dratedt(ir3a) * sc3a + rate(ir3a) * sc3adt;
    rate(ir3a)    = rate(ir3a) * sc3a;

    // c12 to o16
    jscr = jscr + 1;
    screen5(pstate,jscr,sc1a,sc1adt,sc1add);

    dratedt(ircag)  = dratedt(ircag) * sc1a + rate(ircag) * sc1adt;
    rate(ircag)     = rate(ircag) * sc1a;

    // c12 + c12
    jscr = jscr + 1;
    screen5(pstate,jscr,sc1a,sc1adt,sc1add);

    dratedt(ir1212) = dratedt(ir1212) * sc1a + rate(ir1212) * sc1adt;
    rate(ir1212)    = rate(ir1212) * sc1a;

    // c12 + o16
    jscr = jscr + 1;
    screen5(pstate,jscr,sc1a,sc1adt,sc1add);

    dratedt(ir1216) = dratedt(ir1216) * sc1a + rate(ir1216) * sc1adt;
    rate(ir1216)    = rate(ir1216) * sc1a;

    // o16 + o16
    jscr = jscr + 1;
    screen5(pstate,jscr,sc1a,sc1adt,sc1add);

    dratedt(ir1616) = dratedt(ir1616) * sc1a + rate(ir1616) * sc1adt;
    rate(ir1616)    = rate(ir1616) * sc1a;

    // o16 to ne20
    jscr = jscr + 1;
    screen5(pstate,jscr,sc1a,sc1adt,sc1add);

    dratedt(iroag) = dratedt(iroag) * sc1a + rate(iroag) * sc1adt;
    rate(iroag)    = rate(iroag) * sc1a;

    // o16 to mg24
    jscr = jscr + 1;
    screen5(pstate,jscr,sc1a,sc1adt,sc1add);

    dratedt(irneag) = dratedt(irneag) * sc1a + rate(irneag) * sc1adt;
    rate(irneag)    = rate(irneag) * sc1a;

    // mg24 to si28
    jscr = jscr + 1;
    screen5(pstate,jscr,sc1a,sc1adt,sc1add);

    dratedt(irmgag) = dratedt(irmgag) * sc1a + rate(irmgag) * sc1adt;
    rate(irmgag)    = rate(irmgag) * sc1a;

    // ca40 to ti44
    jscr = jscr + 1;
    screen5(pstate,jscr,sc1a,sc1adt,sc1add);

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


void dfdy_isotopes_iso7(Array1D<Real, 0, Species::NSpec-1> const& y,
                        burn_t const& state, rate_t const& rr,
                        Array2D<Real,0, njrows-1, 0, njcols-1>& jac)
{

    // this routine sets up the dense iso7 jacobian for the isotopes

    Array1D<Real, 1, 8> b;

    // set up the jacobian
    // 4he jacobian elements
    // d(he4)/d(he4)
    b(1) = -1.5e0_rt * y(He4) * y(He4) * rr.rates(1,ir3a);
    b(2) = -y(C12) * rr.rates(1,ircag);
    b(3) = -y(O16) * rr.rates(1,iroag);
    b(4) = -y(Ne20) * rr.rates(1,irneag);
    b(5) = -y(Mg24) * rr.rates(1,irmgag);
    b(6) = -7.0e0_rt * rr.rates(1,irsi2ni);
    b(7) = -7.0e0_rt * rr.rates(3,irsi2ni) * y(He4);
    b(8) =  7.0e0_rt * rr.rates(3,irni2si) * y(Ni56);

    jac(He4,He4) = esum8(b);

    // d(he4)/d(c12)
    b(1) =  3.0e0_rt * rr.rates(1,irg3a);
    b(2) = -y(He4) * rr.rates(1,ircag);
    b(3) =  y(C12) * rr.rates(1,ir1212);
    b(4) =  0.5e0_rt * y(O16) * rr.rates(1,ir1216);

    jac(He4,C12) = esum4(b);

    // d(he4)/d(o16)
    b(1) =  rr.rates(1,iroga);
    b(2) =  0.5e0_rt * y(C12) * rr.rates(1,ir1216);
    b(3) =  y(O16) * rr.rates(1,ir1616);
    b(4) = -y(He4) * rr.rates(1,iroag);

    jac(He4,O16) = esum4(b);

    // d(he4)/d(ne20)
    b(1) =  rr.rates(1,irnega);
    b(2) = -y(He4) * rr.rates(1,irneag);

    jac(He4,Ne20) = b(1) + b(2);

    // d(he4)/d(mg24)
    b(1) =  rr.rates(1,irmgga);
    b(2) = -y(He4) * rr.rates(1,irmgag);

    jac(He4,Mg24) = b(1) + b(2);

    // d(he4)/d(si28)
    b(1) =  rr.rates(1,irsiga);
    b(2) = -7.0e0_rt * rr.rates(4,irsi2ni) * y(He4);

    jac(He4,Si28) = b(1) + b(2);

    // d(he4)/d(ni56)
    b(1) =  7.0e0_rt * rr.rates(1,irni2si);

    jac(He4,Ni56) = b(1);



    // 12c jacobian elements
    // d(c12)/d(he4)
    b(1) =  0.5e0_rt * y(He4) * y(He4) * rr.rates(1,ir3a);
    b(2) = -y(C12) * rr.rates(1,ircag);

    jac(C12,He4) = b(1) + b(2);

    // d(c12)/d(c12)
    b(1) = -rr.rates(1,irg3a);
    b(2) = -y(He4) * rr.rates(1,ircag);
    b(3) = -2.0e0_rt * y(C12) * rr.rates(1,ir1212);
    b(4) = -y(O16) * rr.rates(1,ir1216);

    jac(C12,C12) = esum4(b);

    // d(c12)/d(o16)
    b(1) =  rr.rates(1,iroga);
    b(2) = -y(C12) * rr.rates(1,ir1216);

    jac(C12,O16) = b(1) + b(2);


    // 16o jacobian elements
    // d(o16)/d(he4)
    b(1) =  y(C12) * rr.rates(1,ircag);
    b(2) = -y(O16) * rr.rates(1,iroag);

    jac(O16,He4) = b(1) + b(2);

    // d(o16)/d(c12)
    b(1) =  y(He4) * rr.rates(1,ircag);
    b(2) = -y(O16) * rr.rates(1,ir1216);

    jac(O16,C12) = b(1) + b(2);

    // d(o16)/d(o16)
    b(1) = -rr.rates(1,iroga);
    b(2) = -y(C12) * rr.rates(1,ir1216);
    b(3) = -2.0e0_rt * y(O16) * rr.rates(1,ir1616);
    b(4) = -y(He4) * rr.rates(1,iroag);

    jac(O16,O16) = esum4(b);

    // d(o16)/d(ne20)
    b(1) =  rr.rates(1,irnega);

    jac(O16,Ne20) = b(1);



    // 20ne jacobian elements
    // d(ne20)/d(he4)
    b(1) =  y(O16) * rr.rates(1,iroag) - y(Ne20) * rr.rates(1,irneag);

    jac(Ne20,He4) = b(1);

    // d(ne20)/d(c12)
    b(1) =  y(C12) * rr.rates(1,ir1212);

    jac(Ne20,C12) = b(1);

    // d(ne20)/d(o16)
    b(1) =  y(He4) * rr.rates(1,iroag);

    jac(Ne20,O16) = b(1);

    // d(ne20)/d(ne20)
    b(1) = -rr.rates(1,irnega) - y(He4) * rr.rates(1,irneag);

    jac(Ne20,Ne20) = b(1);

    // d(ne20)/d(mg24)
    b(1) =  rr.rates(1,irmgga);

    jac(Ne20,Mg24) = b(1);



    // 24mg jacobian elements
    // d(mg24)/d(he4)
    b(1) =  y(Ne20) * rr.rates(1,irneag);
    b(2) = -y(Mg24) * rr.rates(1,irmgag);

    jac(Mg24,He4) = b(1) + b(2);

    // d(mg24)/d(c12)
    b(1) =  0.5e0_rt * y(O16) * rr.rates(1,ir1216);

    jac(Mg24,C12) = b(1);

    // d(mg24)/d(o16)
    b(1) =  0.5e0_rt * y(C12) * rr.rates(1,ir1216);

    jac(Mg24,O16) = b(1);

    // d(mg24)/d(ne20)
    b(1) =  y(He4) * rr.rates(1,irneag);

    jac(Mg24,Ne20) = b(1);

    // d(mg24)/d(mg24)
    b(1) = -rr.rates(1,irmgga);
    b(2) = -y(He4) * rr.rates(1,irmgag);

    jac(Mg24,Mg24) = b(1) + b(2);

    // d(mg24)/d(si28)
    b(1) =  rr.rates(1,irsiga);

    jac(Mg24,Si28) = b(1);

    // 28si jacobian elements
    // d(si28)/d(he4)
    b(1) =  y(Mg24) * rr.rates(1,irmgag);
    b(2) = -rr.rates(1,irsi2ni);
    b(3) = -rr.rates(3,irsi2ni) * y(He4);
    b(4) =  rr.rates(3,irni2si) * y(Ni56);

    jac(Si28,He4) = esum4(b);

    // d(si28)/d(c12)
    b(1) =  0.5e0_rt * y(O16) * rr.rates(1,ir1216);

    jac(Si28,C12) = b(1);

    // d(si28)/d(o16)
    b(1) =  y(O16) * rr.rates(1,ir1616);
    b(2) =  0.5e0_rt * y(C12) * rr.rates(1,ir1216);

    jac(Si28,O16) = b(1) + b(2);

    // d(si28)/d(mg24)
    b(1) =  y(He4) * rr.rates(1,irmgag);

    jac(Si28,Mg24) = b(1);

    // d(si28)/d(si28)
    b(1) = -rr.rates(1,irsiga);
    b(2) = -rr.rates(4,irsi2ni) * y(He4);

    jac(Si28,Si28) = b(1) + b(2);

    // d(si28)/d(ni56)
    b(1) =  rr.rates(1,irni2si);

    jac(Si28,Ni56) = b(1);

    // ni56 jacobian elements
    // d(ni56)/d(he4)
    b(1) =  rr.rates(1,irsi2ni);
    b(2) =  rr.rates(3,irsi2ni) * y(He4);
    b(3) = -rr.rates(3,irni2si) * y(Ni56);

    jac(Ni56,He4) = esum3(b);

    // d(ni56)/d(si28)
    b(1) = rr.rates(4,irsi2ni) * y(He4);

    jac(Ni56,Si28) = b(1);

    // d(ni56)/d(ni56)
    b(1) = -rr.rates(1,irni2si);

    jac(Ni56,Ni56) = b(1);

}


void ener_gener_rate(Array1D<Real, 0, Burn::neqs-1> const& dydt, Real& enuc)
{

    // Computes the instantaneous energy generation rate

    Real Xdot = zero;

    // Sum the mass fraction time derivatives
    for (int i = 0; i < Species::NSpec; ++i) {
        Xdot += dydt(i) * mion(i);
    }

    // This is basically e = m c**2
    enuc = Xdot * C::enuc_conv2;

}


void set_up_screening_factors()
{
    // Compute and store the more expensive screening factors

    // note: we need to set these up in the same order that we evaluate the
    // rates in actual_rhs.H (yes, it's ugly)
    add_screening_factor(0, zion( He4), aion( He4), zion(He4), aion(He4));
    add_screening_factor(1, zion( He4), aion( He4),  4.0e0_rt,  8.0e0_rt);
    add_screening_factor(2, zion( C12), aion( C12), zion(He4), aion(He4));
    add_screening_factor(3, zion( C12), aion( C12), zion(C12), aion(C12));
    add_screening_factor(4, zion( C12), aion( C12), zion(O16), aion(O16));
    add_screening_factor(5, zion( C12), aion( C12), zion(O16), aion(O16));
    add_screening_factor(6, zion( O16), aion( O16), zion(He4), aion(He4));
    add_screening_factor(7, zion(Ne20), aion(Ne20), zion(He4), aion(He4));
    add_screening_factor(8, zion(Mg24), aion(Mg24), zion(He4), aion(He4));
    add_screening_factor(9,  20.0e0_rt,  40.0e0_rt, zion(He4), aion(He4));

}

#endif
