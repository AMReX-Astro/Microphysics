#ifndef actual_rhs_H
#define actual_rhs_H

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include <extern_parameters.H>
#include <microphysics_math.H>
#include <actual_network.H>
#include <burn_type.H>
#include <rate_type.H>
#include <screen.H>
#include <temperature_integration.H>
#include <fundamental_constants.H>

using namespace amrex;
using namespace ArrayUtil;
using namespace Rates;

void actual_rhs_init ();

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void screen (const Real& temp, const Real& dens, const Array1D<Real, 1, NumSpec>& ymol,
             Array1D<Real, 1, NumRates>& rates, Array1D<Real, 1, NumRates>& dratesdt)
{
    using namespace Species;
    using namespace Rates;

    Real scorr1, dscorr1dt, dscorr1dd;
    Real scorr2, dscorr2dt, dscorr2dd;
    Real scorr,  dscorrdt, dscorrdd;

    Array1D<Real, 1, NumRates> rates_in, dratesdt_in;

    for (int i = 1; i <= NumRates; ++i) {
        rates_in(i)    = rates(i);
        dratesdt_in(i) = dratesdt(i);
    }

    plasma_state_t state;
    fill_plasma_state(state, temp, dens, ymol);

    int jscr = 0;
    screen5(state, jscr,
            zion[He4-1], aion[He4-1], zion[He4-1], aion[He4-1],
            scorr1, dscorr1dt, dscorr1dd);

    jscr++;
    screen5(state, jscr,
            zion[He4-1], aion[He4-1], 4.0e0_rt, 8.0e0_rt,
            scorr2, dscorr2dt, dscorr2dd);

    scorr    = scorr1 * scorr2;
    dscorrdt = dscorr1dt * scorr2 + scorr1 * dscorr2dt;

    rates(ir3a)    = rates_in(ir3a) * scorr;
    dratesdt(ir3a) = dratesdt_in(ir3a) * scorr + rates_in(ir3a) * dscorrdt;

    // C12 + alpha --> O16
    jscr++;
    screen5(state, jscr,
            zion[C12-1], aion[C12-1], zion[He4-1], aion[He4-1],
            scorr, dscorrdt, dscorrdd);

    rates(ircago)    = rates_in(ircago) * scorr;
    dratesdt(ircago) = dratesdt_in(ircago) * scorr + rates_in(ircago) * dscorrdt;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void make_rates (const Real& temp, const Real& dens,
                 Array1D<Real, 1, NumRates>& rates,
                 Array1D<Real, 1, NumRates>& dratesdt)
{
    using namespace Rates;

    // rates given in terms of molar fractions

    Real t9r, t9r32, t9ri, t9ri2, t9, t9i, t913, t9i13, t923;
    Real t9i23, t943, t9i43, t932, t9i32, t953, t9i53, t92, t9i2;

    Real  term,  dtermdt;
    Real r2abe, dr2abedt;
    Real rbeac, drbeacdt;

    // some scratch terms
    Real a, b, b1, b2, c;
    Real dadt, dbdt, db1dt, db2dt, dcdt;

    constexpr Real ONE_TWELVTH  = 1.0_rt / 12.0_rt;
    constexpr Real FIVE_SIXTHS  = 5.0_rt / 6.0_rt;
    constexpr Real FIVE_THIRDS  = 5.0_rt / 3.0_rt;
    constexpr Real THREE_HALVES = 3.0_rt / 2.0_rt;
    constexpr Real T2T9         = 1.0e-9_rt;

    t9r   = temp * T2T9;
    t9r32 = std::pow(t9r, THREE_HALVES);
    t9ri  = 1.0_rt / t9r;
    t9ri2 = t9ri * t9ri;
    t9    = amrex::min(t9r, 10.0e0_rt);
    t9i   = 1.0_rt / t9;
    t913  = std::pow(t9, (1.0_rt / 3.0_rt));
    t9i13 = 1.0_rt / t913;
    t923  = std::pow(t9, (2.0_rt / 3.0_rt));
    t9i23 = 1.0_rt / t923;
    t943  = std::pow(t9, (4.0_rt / 3.0_rt));
    t9i43 = 1.0_rt / t943;
    t932  = std::pow(t9, THREE_HALVES);
    t9i32 = 1.0_rt / t932;
    t953  = std::pow(t9, FIVE_THIRDS);
    t9i53 = 1.0_rt / t953;
    t92   = t9 * t9;
    t9i2  = 1.0_rt / t92;

    for (int i = 1; i <= NumRates; ++i) {
        rates(i)    = 0.0_rt;
        dratesdt(i) = 0.0_rt;
    }

    // triple alpha to c12 (divided by 3! = 6)
    // from cf88

    // q = -0.092;     2 He4 --> Be8
    a    = (7.4e5_rt * t9i32) * std::exp(-1.0663_rt * t9i);
    dadt = -a * THREE_HALVES * t9i + a * t9i2 * 1.0663_rt;

    b    =  4.164e9_rt * t9i23 * std::exp(-13.49_rt * t9i13 - t92 / 9.604e-3_rt);
    dbdt = -b * (2.0_rt / 3.0_rt) * t9i + b * (13.49_rt * (1.0_rt / 3.0_rt) * t9i43 - 2.0_rt * t9 / 9.604e-3_rt);

    c    = 1.0_rt + 3.1e-2_rt * t913 + 8.009_rt * t923 + 1.732_rt * t9 + 49.883_rt * t943 + 27.426_rt * t953;
    dcdt = 3.1e-2_rt * (1.0_rt / 3.0_rt) * t9i23 + 8.009_rt * (2.0_rt / 3.0_rt) * t9i13 +
           1.732_rt + 49.883_rt * (4.0_rt / 3.0_rt) * t913 + 27.426_rt * FIVE_THIRDS * t923;

    r2abe    = a + b * c;
    dr2abedt = dadt + dbdt * c + b * dcdt;


    // q = 7.367;      He4 + Be8 --> C12
    a    = (130_rt * t9i32) * std::exp(-3.3364_rt * t9i);
    dadt = -a * THREE_HALVES * t9i + a * 3.3364_rt * t9i2;

    b    = 2.51e7_rt * t9i23 * std::exp(-23.57_rt * t9i13 - t92 / 0.055225_rt);
    dbdt = b * (2.0_rt / 3.0_rt) * t9i + b * (23.57_rt * (1.0_rt / 3.0_rt) * t9i43 - 2.0_rt * t9 / 0.055225_rt);

    c    = 1.0_rt + 0.018_rt * t913 + 5.249_rt * t923 + 0.65_rt * t9 + 19.176_rt * t943 + 6.034_rt * t953;
    dcdt = 0.018_rt * (1.0_rt / 3.0_rt) * t9i23 + 5.249_rt * (2.0_rt / 3.0_rt) * t9i13 +
           0.65_rt + 19.176_rt * (4.0_rt / 3.0_rt) * t913 + 6.034_rt * FIVE_THIRDS * t923;

    rbeac    = a + b * c;
    drbeacdt = dadt + dbdt * c + b * dcdt;

    // q = 7.275;      total reaction

    a    = 2.9e-16_rt * r2abe * rbeac;
    dadt = 2.9e-16_rt * (dr2abedt * rbeac + r2abe * drbeacdt);

    if (t9 > 8e-2_rt) {
       b    = 0.1_rt * 1.35e-7_rt * t9i32 * std::exp(-24.811_rt * t9i);
       dbdt = -b * (2.0_rt / 3.0_rt) * t9i + 24.811_rt * b * t9i2;

       term    = a + b;
       dtermdt = dadt + dbdt;
    } else {
       b1    = 1.0_rt + 4.0_rt * std::exp(-std::pow(2.5e-2_rt * t9i, 3.263_rt));
       db1dt = 4.0_rt * 3.263_rt * std::pow(2.5e-2_rt * t9i, 3.263_rt) *
               t9i * std::exp(-std::pow(2.5e-2_rt * t9i, 3.263_rt));

       b2    = 1.0_rt + 4.0_rt * std::exp(-std::pow(t9 / 2.5e-2_rt, 9.227_rt));
       db2dt = -4.0_rt * 9.227_rt * std::pow(t9 / 2.5e-2_rt, 9.227_rt) *
                t9i * std::exp(-std::pow(t9 / 2.5e-2_rt, 9.227_rt));

       b    = 1.e-2_rt + 0.2_rt * b1 / b2;
       dbdt = 0.2_rt * (db1dt / b2 - b1 * db2dt / (b2 * b2));

       c    = 0.1_rt * 1.35e-7_rt * t9i32 * std::exp(-24.811_rt * t9i);
       dcdt = -c * THREE_HALVES * t9i + 24.811_rt * c * t9i2;

       term    = a * b + c;
       dtermdt = dadt * b + a * dbdt + dcdt;
    }

    rates(ir3a)    = term * (dens * dens) / 6.0_rt;
    dratesdt(ir3a) = dtermdt * T2T9 * (dens * dens) / 6.0_rt;


    // c12 + he4 to o16
    // 1.7 time cf88 rate: see Weaver & Woosley Phy. Rep. 227 (1993)
    //                     and Garnett Nuc. Phys. A. 621  (1997)
    // q = 7.162
    b1    = 1.0_rt + 0.0489_rt * t9i23;
    db1dt = -0.0489_rt * (2.0_rt / 3.0_rt) * t9i53;

    b2    = -32.120_rt * t9i13 - t92/(3.496_rt * 3.496_rt);
    db2dt = 32.120_rt * (1.0_rt / 3.0_rt) * t9i43 - 2.0_rt * t9 / (3.496_rt * 3.496_rt);

    a    = t9i2 * b1 * b1 * std::exp(b2);
    dadt = a * (-2.0_rt * t9i + 2.0_rt * db1dt / b1 + db2dt);

    //------------------------

    b2    = 1.0_rt + 0.2654_rt * t9i23;
    db2dt = -0.2654_rt * (2.0_rt / 3.0_rt) * t9i53;

    c    = -32.120_rt * t9i13;
    dcdt = -(1.0_rt / 3.0_rt) * c * t9i;

    b1    = t9i2 * b2 * b2 * std::exp(c);
    db1dt = b1 * (-2.0_rt * t9i + 2.0_rt * db2dt / b2 + dcdt);

    //------------------------

    c    = -27.499_rt * t9i;
    dcdt = - c * t9i;

    b2    = t9i32 * std::exp(c);
    db2dt = b2 * (-THREE_HALVES * t9i + dcdt);

    //------------------------

    c    = t92 * t92 * t9 * std::exp(-15.541_rt * t9i);
    dcdt = c * (5.0_rt * t9i + 15.541_rt * t9i2);

    //------------------------

    term    = 1.04e8_rt * a + 1.76e8_rt * b1 + 1.25e3_rt * b2 + 1.43e-2_rt * c;
    dtermdt = 1.04e8_rt * dadt + 1.76e8_rt * db1dt + 1.25e3_rt * db2dt + 1.43e-2_rt * dcdt;

    term    = 1.7_rt * term;
    dtermdt = 1.7_rt * term;

    rates(ircago)    = term * dens;
    dratesdt(ircago) = dtermdt * T2T9 * dens;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_rates (burn_t& state, Array1D<rate_t, 1, Rates::NumGroups>& rr)
{
    using namespace Rates;

    Real temp = state.T;
    Real dens = state.rho;

    Array1D<Real, 1, NumSpec> ymol;
    for (int i = 1; i <= NumSpec; ++i) {
        ymol(i) = state.xn[i-1] * aion_inv[i-1];
    }

    Array1D<Real, 1, NumRates> rates, dratesdt;
    make_rates(temp, dens, rates, dratesdt);

    screen(temp, dens, ymol, rates, dratesdt);

    for (int i = 1; i <= NumRates; ++i) {
        rr(1).rates(i) = rates(i);
        rr(2).rates(i) = dratesdt(i);
    }
}


template<class T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ener_gener_rate (T& dydt, Real& enuc)
{
    using namespace triple_alpha_plus_cago;

    enuc = 0.0_rt;

    for (int i = 1; i <= NumSpec; ++i) {
        enuc += dydt(i) * bion(i);
    }
    enuc *= C::MeV2eV * C::ev2erg * C::n_A;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void dydt (const Array1D<Real, 1, NumSpec>& ymol, const Array1D<Real, 1, NumRates>& rates,
           Array1D<Real, 1, NumSpec>& ydot)
{
    using namespace Species;
    using namespace Rates;

    for (int i = 1; i <= NumSpec; ++i) {
        ydot(i) = 0.0_rt;
    }

    // He4 reactions
    ydot(He4) = -3.0_rt * ymol(He4) * ymol(He4) * ymol(He4) * rates(ir3a)
                -1.0_rt * ymol(C12) * ymol(He4) * rates(ircago);

    // C12 reactions
    ydot(C12) = 1.0_rt * ymol(He4) * ymol(He4) * ymol(He4) * rates(ir3a)
                -1.0_rt * ymol(C12) * ymol(He4) * rates(ircago);

    // O16 reactions
    ydot(O16) = 1.0_rt * ymol(C12) * ymol(He4) * rates(ircago);
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_rhs (burn_t& state, Array1D<Real, 1, neqs>& ydot)
{
    using namespace Rates;

    for (int i = 1; i <= neqs; ++i) {
        ydot(i) = 0.0_rt;
    }

    Array1D<Real, 1, NumSpec> ymol;
    for (int i = 1; i <= NumSpec; ++i) {
        ymol(i) = state.xn[i-1] * aion_inv[i-1];
    }

    // set up the species ODEs for the reaction network
    // species inputs are in molar fractions but come out in mass fractions

    Array1D<rate_t, 1, Rates::NumGroups> rr;
    get_rates(state, rr);

    Array1D<Real, 1, NumRates> rates;
    for (int i = 1; i <= NumRates; ++i) {
        rates(i) = rr(1).rates(i);
    }

    Array1D<Real, 1, NumSpec> yderivs;
    dydt(ymol, rates, yderivs);

    for (int i = 1; i <= NumSpec; ++i) {
        ydot(i) = yderivs(i);
    }

    // Energy generation rate

    ener_gener_rate(ydot, ydot(net_ienuc));

#ifdef STRANG
    ydot(net_itemp) = temperature_rhs(state, ydot(net_ienuc));
#endif
}


template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_jac (burn_t& state, MatrixType& jac)
{
    using namespace Species;
    using namespace Rates;

    Array1D<rate_t, 1, Rates::NumGroups> rr;
    get_rates(state, rr);

    Array1D<Real, 1, NumRates> rates, dratesdt;
    for (int i = 1; i <= NumRates; ++i) {
        rates(i)    = rr(1).rates(i);
        dratesdt(i) = rr(2).rates(i);
    }

    // initialize
    jac.zero();

    Array1D<Real, 1, NumSpec> ymol;
    for (int i = 1; i <= NumSpec; ++i) {
        ymol(i) = state.xn[i-1] * aion_inv[i-1];
    }

    // ======================================================================
    // THESE ARE IN TERMS OF MOLAR FRACTIONS

    // helium jacobian elements
    jac(He4, He4) = -9.0_rt * ymol(He4) * ymol(He4) * rates(ir3a) - 1.0_rt * ymol(C12) * rates(ircago);
    jac(He4, C12) = -1.0_rt * ymol(He4) * rates(ircago);

    // carbon jacobian elements
    jac(C12,He4) =  3.0_rt * ymol(He4) * ymol(He4) * rates(ir3a) - 1.0_rt * ymol(C12) * rates(ircago);
    jac(C12,C12) = -1.0_rt * ymol(He4) * rates(ircago);

    // oxygen jacobian elements
    jac(O16, He4) = 1.0_rt * ymol(C12) * rates(ircago);
    jac(O16, C12) = 1.0_rt * ymol(He4) * rates(ircago);

    // ======================================================================

    // Add the temperature derivatives: df(y_i) / dT

    Array1D<Real, 1, NumSpec> ydot;
    dydt(ymol, dratesdt, ydot);

    for (int j = 1; j <= NumSpec; ++j) {
        jac(j, net_itemp) = ydot(j);
    }

    // Energy generation rate Jacobian elements with respect to species

    for (int j = 1; j <= NumSpec; ++j) {
        auto jac_slice_2 = [&](int i) -> Real { return jac.get(i, j); };
        ener_gener_rate(jac_slice_2, jac(net_ienuc, j));
    }

    // Jacobian elements with respect to temperature

    ener_gener_rate(ydot, jac(net_ienuc, net_itemp));

    temperature_jac(state, jac);
}


// Compute and store the more expensive screening factors

AMREX_INLINE
void set_up_screening_factors ()
{
    using namespace Species;

    // note: it is critical that these are called in the exact order
    // that the screening calls are done in the RHS routine, since we
    // use that order in the screening

    int jscr = 0;
    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], 4.0e0_rt, 8.0e0_rt);

    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[He4-1], aion[He4-1]);
}

#endif
