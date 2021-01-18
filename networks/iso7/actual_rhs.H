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
#include <rhs_utilities.H>

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
                 Array2D<Real, 1, NumSpec, 1, Rates::NumRates>& dratedY)
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
       dratedY(He4, irsi2ni) = 3.0e0_rt * rate(irsi2ni)/y(He4);
       dratedY(Si28, irsi2ni) = yeff_ca40*denom*rate(ircaag);
       dratedt(irsi2ni)  = (yeff_ca40dt*rate(ircaag)
                          + yeff_ca40*dratedt(ircaag))*denom*y(Si28)*1.0e-9_rt;

       if (denom != 0.0_rt) {

          zz     = 1.0e0_rt/denom;
          rate(irni2si) = amrex::min(1.0e10_rt,yeff_ti44*rate(irtiga)*zz);

          if (rate(irni2si) == 1.0e10_rt) {
             dratedY(He4, irni2si) = 0.0e0_rt;
             dratedt(irni2si)  = 0.0e0_rt;
          } else {
             dratedY(He4, irni2si) = -3.0e0_rt * rate(irni2si)/y(He4);
             dratedt(irni2si)  = (yeff_ti44dt*rate(irtiga)
                                + yeff_ti44*dratedt(irtiga))*zz*1.0e-9_rt;
          }

       }

    }
}



// RHS version (single rate group)
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void evaluate_rates(burn_t const& state, rate_t& rr)
{
    Real rho, temp;
    Array1D<Real, 1, NumSpec> y;

    // Temporary arrays for storing unused derivatives.
    // A useful cleanup would be to rewrite iso7tab, iso7rat,
    // and screen_iso7 so that there are versions that do not
    // evaluate the derivatives.
    Array1D<Real, 1, Rates::NumRates> dratedt;
    Array2D<Real, 1, NumSpec, 1, Rates::NumRates> dratedY;

    // Get the data from the state
    rho  = state.rho;
    temp = state.T;

    for (int i = 1; i <= NumSpec; ++i) {
        y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    // Get the raw reaction rates
    if (use_tables)
        iso7tab(temp, rho, rr.rates, dratedt);
    else
        iso7rat(temp, rho, rr.rates, dratedt);

    // Do the screening (depends on composition)
    screen_iso7(temp, rho, y, rr.rates, dratedt, dratedY);
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_rates(burn_t const& state, Array1D<rate_t, 1, Rates::NumGroups>& rr)
{
    Real rho, temp;
    Array1D<Real, 1, NumSpec> y;
    Array1D<Real, 1, Rates::NumRates> rate, dratedt;
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

    // Save the rate data, for the Jacobian later if we need it.
    for (int i = 1; i <= Rates::NumRates; ++i) {
        rr(1).rates(i) = rate(i);
        rr(2).rates(i) = dratedt(i);
    }

    for (int j = 1; j <= NumSpec; ++j) {
        for (int i = 1; i <= Rates::NumRates; ++i) {
            rr(j+2).rates(i) = dratedY(j, i);
        }
    }
}


template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void dfdy_isotopes_iso7(Array1D<Real, 1, NumSpec> const& y,
                        burn_t const& state, Array1D<rate_t, 1, Rates::NumGroups> const& rr,
                        MatrixType& jac)
{

    // this routine sets up the dense iso7 jacobian for the isotopes

    using namespace Species;
    using namespace Rates;

    for (int spec1 = 1; spec1 <= NumSpec; ++spec1) {

        for (int spec2 = 1; spec2 <= NumSpec; ++spec2) {

            // d(spec1) / d(spec2)

            Array1D<Real, 1, Rates::NumRates> b = {0.0};

            int b_ctr = 0;

            for (int rate = 1; rate <= NumRates; ++rate) {
                jac_t jac_data = RHS::jac_data(spec1, spec2, rate);

                if (std::abs(jac_data.prefactor) > 0.0_rt) {
                    ++b_ctr;

                    Real term = jac_data.prefactor;

                    if (jac_data.specindex1 >= 0) {
                        term *= y(jac_data.specindex1);
                    }

                    if (jac_data.specindex2 >= 0) {
                        term *= y(jac_data.specindex2);
                    }

                    if (jac_data.specindex3 >= 0) {
                        term *= y(jac_data.specindex3);
                    }

                    term *= rr(1).rates(rate);

                    b(b_ctr) = term;
                }

                if (std::abs(jac_data.rate_prefactor) > 0.0_rt) {
                    ++b_ctr;

                    Real term = jac_data.rate_prefactor;

                    if (jac_data.rate_specindex1 >= 0) {
                        term *= y(jac_data.rate_specindex1);
                    }

                    if (jac_data.rate_specindex2 >= 0) {
                        term *= y(jac_data.rate_specindex2);
                    }

                    if (jac_data.rate_specindex3 >= 0) {
                        term *= y(jac_data.rate_specindex3);
                    }

                    term *= rr(2 + spec2).rates(rate);

                    b(b_ctr) = term;
                }
            }

            jac(spec1, spec2) = esum<Rates::NumRates>(b);

        }

    }

}

template<class T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ener_gener_rate(T const& dydt, Real& enuc)
{

    using namespace iso7;

    // Computes the instantaneous energy generation rate

    Real Xdot = 0.0_rt;

    // Sum the mass fraction time derivatives
    for (int i = 1; i <= NumSpec; ++i) {
        Xdot += dydt(i) * mion(i);
    }

    // This is basically e = m c**2
    enuc = Xdot * C::Legacy::enuc_conv2;

}


// This version adds the neutrino losses directly
// and is intended for use in the RHS only.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real ener_rhs(const burn_t& state, Array1D<Real, 1, NumSpec>& dydt)
{
    using namespace iso7;

    Real Xdot = 0.0_rt;

    // Sum the mass fraction time derivatives

    for (int i = 1; i <= NumSpec; ++i) {
        Xdot += dydt(i) * mion(i);
    }

    // This is basically e = m c**2

    Real dedt = Xdot * C::Legacy::enuc_conv2;

    // Get the neutrino losses

    Real sneut, dsneutdt, dsneutdd, snuda, snudz;
    sneut5(state.T, state.rho, state.abar, state.zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);

    dedt -= sneut;

    return dedt;
}



AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_rhs(burn_t& state, Array1D<Real, 1, neqs>& ydot)
{
    // Evaluate the rates

    rate_t rr;
    evaluate_rates(state, rr);

    // Call the RHS to get dydt

    Array1D<Real, 1, NumSpec> spec_rhs = species_rhs(state, rr);

    for (int n = 1; n <= NumSpec; ++n) {
        ydot(n) = spec_rhs(n);
    }

    // Instantaneous energy generation rate

    ydot(net_ienuc) = ener_rhs(state, spec_rhs);

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

    const int index_rate = 2;

    Array1D<Real, 1, NumSpec> yderivs = species_rhs(state, rr(index_rate));

    for (int i = 1; i <= NumSpec; ++i) {
        jac(i,net_itemp) = yderivs(i);
    }

    ener_gener_rate(yderivs, jac(net_ienuc,net_itemp));

    jac.add(net_ienuc, net_itemp, -dsneutdt);

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
