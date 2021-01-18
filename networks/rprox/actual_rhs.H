#ifndef actual_rhs_H
#define actual_rhs_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <extern_parameters.H>
#include <actual_network.H>
#include <rate_type.H>
#include <temperature_integration.H>
#include <rprox_tfactors.H>
#include <rprox_rates.H>
#include <screen.H>

using namespace amrex;
using namespace ArrayUtil;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void make_rates (const Real& T9, const Real& dens, Array1D<Real, 1, NumSpec>& y,
                 burn_t& state, Array1D<rate_t, 1, Rates::NumGroups>& rr)
{
    using namespace Species;
    using namespace Rates;

    // some numbers from appendix C in WW81; these should probably be
    // updated with current rates
    constexpr Real Lweak = 1.05e0_rt; // this is for NS
    // constexpr Real Lweak = 0.107e0_rt; // this is for lower densities

    constexpr Real la2 = 1.0_rt / 15.0_rt; // mean rate from 30s to 56ni from p-capture and beta decays

    for (int j = 1; j <= NumGroups; ++j) {
        for (int i = 1; i <= NumRates; ++i) {
            rr(j).rates(i) = 0.0_rt; // Zero out rates
        }
    }

    tf_t tfactors = calc_tfactors(T9);

    // locally used rates
    Real rate, dratedt;

    // some common parameters
    rr(1).rates(irLweak) = Lweak;
    rr(1).rates(irla2)   = la2;

    // weak rates first
    //
    // 14o(beta nu)14n
    rate_o14_to_n14(tfactors, rate, dratedt);
    rr(1).rates(irwk14o) = rate;

    // 15o(beta nu)15n
    rate_o15_to_n15(tfactors, rate, dratedt);
    rr(1).rates(irwk15o) = rate;

    // 17f(beta nu)17o
    rate_f17_to_o17(tfactors, rate, dratedt);
    rr(1).rates(irwk17f) = rate;

    // these weak rates aren't needed outside of this routine
    // 18ne(beta nu)18f
    Real wk18ne;
    rate_ne18_to_f18(tfactors, wk18ne, dratedt);
    // 19ne(beta nu)19f
    Real wk19ne;
    rate_ne19_to_f19(tfactors, wk19ne, dratedt);

    // now the remaining rates -- we'll do screening inline here

    Real sc1a, sc1adt, sc1add;
    Real sc2a, sc2adt, sc2add;
    Real sc3a, sc3adt, sc3add;

    plasma_state_t pstate;
    fill_plasma_state(pstate, T9 * 1.e9_rt, dens, y);

    // 12c(p,g)13n
    rate_p_c12_to_n13(tfactors, rate, dratedt);

    int jscr = 0;
    screen5(pstate, jscr,
            zion[C12-1], aion[C12-1], zion[H1-1], aion[H1-1],
            sc1a, sc1adt, sc1add);

    rr(1).rates(irpg12c) = dens * rate * sc1a;
    rr(2).rates(irpg12c) = dens * (dratedt * sc1a + rate * sc1adt);

    // triple alpha
    rate_he4_he4_he4_to_c12(tfactors, rate, dratedt);

    jscr = 1;
    screen5(pstate, jscr,
            zion[He4-1], aion[He4-1], zion[He4-1], aion[He4-1],
            sc1a, sc1adt, sc1add);

    jscr = 2;
    screen5(pstate, jscr,
            zion[He4-1], aion[He4-1], 4.0_rt, 8.0_rt,
            sc2a, sc2adt, sc2add);

    sc3a = sc1a * sc2a;
    sc3adt = sc1adt*sc2a + sc1a*sc2adt;

    rr(1).rates(ir3a) = dens * dens * rate * sc3a;
    rr(2).rates(ir3a) = dens * dens * (dratedt * sc3a + rate * sc3adt);

    // 17f(p,g)18ne
    rate_p_f17_to_ne18(tfactors, rate, dratedt);

    jscr = 3;
    screen5(pstate, jscr,
            zion[F17-1], aion[F17-1], zion[H1-1], aion[H1-1],
            sc1a, sc1adt, sc1add);

    rr(1).rates(irpg17f) = dens * rate * sc1a;
    rr(2).rates(irpg17f) = dens * (dratedt * sc1a + rate * sc1adt);

    // 17f(g,p)16o -- not screened
    rate_f17_to_p_o16(tfactors, rate, dratedt);
    rr(1).rates(irgp17f) = rate;
    rr(2).rates(irgp17f) = dratedt;

    // 15o(a,g)19ne
    rate_he4_o15_to_ne19(tfactors, rate, dratedt);

    jscr = 4;
    screen5(pstate, jscr,
            zion[O15-1], aion[O15-1], zion[He4-1], aion[He4-1],
            sc1a, sc1adt, sc1add);

    rr(1).rates(irag15o) = dens * rate * sc1a;
    rr(2).rates(irag15o) = dens * (dratedt * sc1a + rate * sc1adt);

    // 16o(a,g)20ne
    rate_he4_o16_to_ne20(tfactors, rate, dratedt);

    jscr = 5;
    screen5(pstate, jscr,
            zion[O16-1], aion[O16-1], zion[He4-1], aion[He4-1],
            sc1a, sc1adt, sc1add);

    rr(1).rates(irag16o) = dens * rate * sc1a;
    rr(2).rates(irag16o) = dens * (dratedt * sc1a + rate * sc1adt);

    // 16o(p,g)17f
    rate_p_o16_to_f17(tfactors, rate, dratedt);

    jscr = 6;
    screen5(pstate, jscr,
            zion[O16-1], aion[O16-1], zion[H1-1], aion[H1-1],
            sc1a, sc1adt, sc1add);

    rr(1).rates(irpg16o) = dens * rate * sc1a;
    rr(2).rates(irpg16o) = dens * (dratedt * sc1a + rate * sc1adt);

    // 14o(a,p)17f
    rate_he4_o14_to_p_f17(tfactors, rate, dratedt);

    jscr = 7;
    screen5(pstate, jscr,
            zion[O14-1], aion[O14-1], zion[He4-1], aion[He4-1],
            sc1a, sc1adt, sc1add);

    rr(1).rates(irap14o) = dens * rate * sc1a;
    rr(2).rates(irap14o) = dens * (dratedt * sc1a + rate * sc1adt);

    // limit CNO as minimum between 14n(p,g)15o and 15o(beta nu)15n
    // we store the limited rate in irlambCNO; this is lambda_CNO in WW81
    rate_p_n14_to_o15(tfactors, rate, dratedt);

    jscr = 8;
    screen5(pstate, jscr,
            7.0_rt, 14.0_rt, zion[H1-1], aion[H1-1],
            sc1a, sc1adt, sc1add);

    rr(1).rates(irlambCNO) = amrex::min(rr(1).rates(irwk15o),
                                        rate * dens * sc1a * y(H1));
    if (rr(1).rates(irlambCNO) < rr(1).rates(irwk15o)) {
        rr(2).rates(irlambCNO) = dens * y(H1) * (dratedt * sc1a + rate * sc1adt);
        rr(3).rates(dlambCNOdh1) = rate * dens * sc1a;
    }

    // 22mg(...)30s
    // check if this proceeds via p-captures or (a,p) reactions
    // the Lweak is from WW81, eqn C15
    // we store the rate in irlambda1; this is the lambda1 in WW81
    rate_he4_si26_to_p_p29(tfactors, rate, dratedt);

    jscr = 9;
    screen5(pstate, jscr,
            14.0_rt, 26.0_rt, zion[He4-1], aion[He4-1],
            sc1a, sc1adt, sc1add);

    rr(1).rates(irlambda1) = amrex::max(rr(1).rates(irLweak),
                                        dens * y(He4) * rate * sc1a);
    if (rr(1).rates(irlambda1) > rr(1).rates(irLweak)) {
        rr(2).rates(irlambda1) = dens * y(He4) * (dratedt * sc1a + rate *  sc1adt);
        rr(3).rates(dlambda1dhe4) = dens * rate * sc1a;
        // use the sign of state.rates(1, irlambda1) to indicate the value of delta1 in WW81
        // if delta1 = 1, then we multiply the rate by -1
        rr(1).rates(irlambda1) = -1.0_rt * rr(1).rates(irlambda1);
    }

    // 30s(...) 56ni
    // check if this proceeds via p-captures or (a,p) reactions
    // use 44ti(a,p)v47 as a typical limiting rate for the (a,p) process
    // store this in irlambda2; this is lambda2 in WW81
    rate_he4_ti44_to_p_v47(tfactors, rate, dratedt);

    jscr = 10;
    screen5(pstate, jscr,
            22.0_rt, 44.0_rt, zion[He4-1], aion[He4-1],
            sc1a, sc1adt, sc1add);

    rr(1).rates(irlambda2) = amrex::max(rr(1).rates(irla2),
                                        dens * y(He4) * rate * sc1a);
    if (rr(1).rates(irlambda2) > rr(1).rates(irla2)) {
        rr(2).rates(irlambda2) = dens * y(He4) * (dratedt * sc1a + rate * sc1adt);
        rr(3).rates(dlambda2dhe4) = dens * rate * sc1a;
        // use the sign of rr(1).rates(irlambda2) to indicate the value of delta2
        // if delta2 = 1, then we multiply the rate by -1
        rr(1).rates(irlambda2) = -1.0_rt * rr(1).rates(irlambda2);
    }

    // form s1 from WW81; branching ratio for 18ne beta decay (wk18ne) vs (a,p)
    // store result in irs1
    // 18ne(a,p)21na
    rate_he4_ne18_to_p_na21(tfactors, rate, dratedt);

    jscr = 11;
    screen5(pstate, jscr,
            10.0_rt, 18.0_rt, zion[He4-1], aion[He4-1],
            sc1a, sc1adt, sc1add);

    rr(1).rates(irs1) = wk18ne / (wk18ne + dens * y(He4) * rate * sc1a);
    rr(2).rates(irs1) = -rr(1).rates(irs1) * dens * y(He4) * (dratedt * sc1a + rate * sc1adt)
                         / (wk18ne + dens * y(He4) * rate * sc1a);
    rr(3).rates(drs1dhe4) = -rr(1).rates(irs1) * dens * rate * sc1a
                             / (wk18ne + dens * y(He4) * rate * sc1a);

    // form r1 from WW81; ranching ratio for 19ne beta decay (wk19ne) vs (p,g)
    // store result in irr1
    // 19ne(p,g)20na
    rate_p_ne19_to_na20(tfactors, rate, dratedt);

    jscr = 12;
    screen5(pstate, jscr,
            10.0_rt, 19.0_rt, zion[H1-1], aion[H1-1],
            sc1a, sc1adt, sc1add);

    rr(1).rates(irr1) = wk19ne / (wk19ne + dens * y(H1) * rate * sc1a);
    rr(2).rates(irr1) = -rr(1).rates(irr1) * dens * y(H1) * (dratedt * sc1a + rate * sc1adt)
                         / (wk19ne + dens * y(H1) * rate * sc1a);
    rr(3).rates(drr1dh1) = -rr(1).rates(irr1) * dens * rate * sc1a
                            / (wk19ne + dens * y(H1) * rate * sc1a);


    //....
    //....  additional coding for proton capture on 56ni to heavier elements
    //....   kludge    56ni+56p -> 2 (56ni) at a rate given by min
    //....   of 56ni(pg) and 57cu decay rate
    //....
    //....  use 56ni rate from wallace and woosley 1981
    Real T9i32 = tfactors.T9i * std::sqrt(tfactors.T9i);
    Real r56pg = dens * (1.29e-02_rt * std::exp(-4.897_rt * tfactors.T9i)
                 + 7.065e+03_rt * std::exp(-20.33_rt * tfactors.T9i)) * T9i32;
    Real dr56pgdt = -(3.0_rt / 2.0_rt) * r56pg * tfactors.T9i +
                    dens * T9i32 * tfactors.T9i*tfactors.T9i *
                    (4.897_rt * 1.29e-2_rt * std::exp(-4.897_rt * tfactors.T9i)
                     + 20.33_rt * 7.065e3_rt * std::exp(-20.33_rt * tfactors.T9i));
    //....  use generic proton separation energy of 400 kev
    //....  8.02 -> 4.64
    //      cutoni=2.08e-10_rt*dens*exp(8.02*T9m1)/T932
    Real cutoni = 2.08e-10_rt * dens * std::exp(4.642_rt * tfactors.T9i) * T9i32;
    Real dcutonidt = cutoni * tfactors.T9i * (-(3.0_rt / 2.0_rt) - 4.642_rt * tfactors.T9i);
    Real r57decay = 3.54_rt;
    Real dr56eff = amrex::min(r56pg, cutoni * r57decay);
    //   rr(3).rates(r56eff) = d56eff;
    //   if (d56eff < r56pg) rr(3).rates(dr56effdt) = r57decay * dcutonidt;
    rr(3).rates(r56eff) = 0.0_rt;
    rr(3).rates(dr56effdt) = 0.0_rt;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void make_ydots (const Array1D<Real, 1, NumSpec>& ymol, const Real& T9,
                 burn_t& state, Array1D<rate_t, 1, Rates::NumGroups>& rr,
                 Array1D<Real, 1, neqs>& dydt, bool doing_dratesdt)
{
    using namespace Species;
    using namespace Rates;

    // initialize
    for (int i = 1; i <= NumSpec; ++i) {
        dydt(i) = 0.0_rt;
    }

    // check to see if we are doing this with the t-derivatives
    // if so, offset our starting index in the rate groups

    int rate_idx = 1;
    if (doing_dratesdt) {
       rate_idx = 2;
    }

    if (!doing_dratesdt) {
        Real ddelta1 = 0.0_rt;
        Real ddelta2 = 0.0_rt;
        // figure out the delta's; we used negative rates to indicate delta=1
        if (rr(1).rates(irlambda1) < 0.0_rt) {
            ddelta1 = 1.0_rt;
            rr(1).rates(irlambda1) = -1.0_rt * rr(1).rates(irlambda1);
        }
        if (rr(1).rates(irlambda2) < 0.0_rt) {
            ddelta2 = 1.0_rt;
            rr(1).rates(irlambda2) = -1.0_rt * rr(1).rates(irlambda2);
        }
        rr(3).rates(delta1) = ddelta1;
        rr(3).rates(delta2) = ddelta2;
    }

    // setup ODEs
    //
    //....
    //.... 12c = 1
    //....
    dydt(C12) = -ymol(C12) * ymol(H1) * rr(rate_idx).rates(irpg12c) +
                 ymol(He4) * ymol(He4) * ymol(He4) * rr(rate_idx).rates(ir3a) / 6.0_rt +
                 ymol(O15) * rr(rate_idx).rates(irlambCNO);
    //....
    //.... 14o = 2
    //....
    dydt(O14) = -ymol(O14) * ymol(He4) * rr(rate_idx).rates(irap14o) -
                 ymol(O14) * rr(rate_idx).rates(irwk14o) +
                 ymol(C12) * ymol(H1) * rr(rate_idx).rates(irpg12c);
    //....
    //.... 15o = 3
    //....
    dydt(O15) = ymol(O14) * rr(rate_idx).rates(irwk14o) -
                ymol(O15) * ymol(He4) * rr(rate_idx).rates(irag15o) -
                ymol(O15) * rr(rate_idx).rates(irlambCNO) +
                ymol(F17) * ymol(H1) * rr(rate_idx).rates(irpg17f) * rr(rate_idx).rates(irs1) +
                ymol(F17) * rr(rate_idx).rates(irwk17f);
    //....
    //.... 16o = 4
    //....
    dydt(O16) = ymol(F17) * rr(rate_idx).rates(irgp17f) -
                ymol(O16) * ymol(H1) * rr(rate_idx).rates(irpg16o) +
                ymol(O15) * ymol(He4) * rr(rate_idx).rates(irr1) * rr(rate_idx).rates(irag15o) -
                ymol(O16) * ymol(He4) * rr(rate_idx).rates(irag16o);
    //....
    //.... 17f = 5
    //....
    dydt(F17) = ymol(O14) * ymol(He4) * rr(rate_idx).rates(irap14o) +
                 ymol(O16) * ymol(H1) * rr(rate_idx).rates(irpg16o) -
                 ymol(F17) * rr(rate_idx).rates(irgp17f) -
                 ymol(F17) * ymol(H1) * rr(rate_idx).rates(irpg17f) -
                 ymol(F17) * rr(rate_idx).rates(irwk17f);
    //....
    //.... 22mg = 6
    //....
    dydt(Mg22) = ymol(O16) * ymol(He4) * rr(rate_idx).rates(irag16o) +
                  ymol(F17) * ymol(H1) * rr(rate_idx).rates(irpg17f) * (1.0_rt - rr(rate_idx).rates(irs1)) +
                  ymol(O15) * ymol(He4) * rr(rate_idx).rates(irag15o) * (1.0_rt - rr(rate_idx).rates(irr1)) -
                  ymol(Mg22) * rr(rate_idx).rates(irlambda1);
    //....
    //.... 30s = 7
    //....
    dydt(S30) = ymol(Mg22) * rr(rate_idx).rates(irlambda1) - ymol(S30) * rr(rate_idx).rates(irlambda2);
    //....
    //.... amax (56ni) = 8  (note that WW81 have a typo -- they write lambda1 here)
    //....
    dydt(Ni56) = ymol(S30) * rr(rate_idx).rates(irlambda2);
    //....
    //.... 4he (alpha) = 9
    //....
    dydt(He4) = -ymol(He4) * ymol(He4) * ymol(He4) * 0.5_rt * rr(rate_idx).rates(ir3a) +
                 ymol(O15) * rr(rate_idx).rates(irlambCNO) -
                 ymol(O14) * ymol(He4) * rr(rate_idx).rates(irap14o) +
                 ymol(F17) * ymol(H1) * rr(rate_idx).rates(irpg17f) * rr(rate_idx).rates(irs1) -
                 ymol(O15) * ymol(He4) * rr(rate_idx).rates(irag15o) * (1.0_rt - rr(rate_idx).rates(irr1)) -
                 ymol(O16) * ymol(He4) * rr(rate_idx).rates(irag16o) -
                 ymol(F17) * ymol(H1) * rr(rate_idx).rates(irpg17f) * (1.0_rt - rr(rate_idx).rates(irs1)) +
                 ymol(F17) * rr(rate_idx).rates(irwk17f) -
                 2.0_rt * ymol(Mg22) * rr(rate_idx).rates(irlambda1) * rr(3).rates(delta1) -
                 6.5e0_rt * ymol(S30) * rr(rate_idx).rates(irlambda2) * rr(3).rates(delta2);
    //....
    //.... 1h (p) = 10
    //....
    dydt(H1) = -ymol(O14) * rr(rate_idx).rates(irwk14o) -
                ymol(O15) * rr(rate_idx).rates(irlambCNO) -
                2.0_rt * ymol(C12) * ymol(H1) * rr(rate_idx).rates(irpg12c) +
                ymol(O14) * ymol(He4) * rr(rate_idx).rates(irap14o) -
                2.0_rt * ymol(F17) * ymol(H1) * rr(rate_idx).rates(irpg17f) * rr(rate_idx).rates(irs1) +
                ymol(F17) * rr(rate_idx).rates(irgp17f) -
                ymol(O16) * ymol(H1) * rr(rate_idx).rates(irpg16o) -
                ymol(O15) * ymol(He4) * rr(rate_idx).rates(irag15o) * rr(rate_idx).rates(irr1) -
                2.0_rt * ymol(O16) * ymol(He4) * rr(rate_idx).rates(irag16o) -
                3.0_rt * ymol(O15) * ymol(He4) * rr(rate_idx).rates(irag15o) * (1.0_rt - rr(rate_idx).rates(irr1)) -
                ymol(F17) * ymol(H1) * rr(rate_idx).rates(irpg17f) * (1.0_rt - rr(rate_idx).rates(irs1)) -
                2.0_rt * ymol(F17) * rr(rate_idx).rates(irwk17f) -
                8.0_rt * ymol(Mg22) * rr(rate_idx).rates(irlambda1) * (1.0_rt - rr(3).rates(delta1)) -
                26.e0_rt * ymol(S30) * rr(rate_idx).rates(irlambda2) * (1.0_rt - rr(3).rates(delta2));


    if (!doing_dratesdt) {
        dydt(Ni56) = dydt(Ni56) + ymol(Ni56) * ymol(H1) * rr(3).rates(r56eff);
        dydt(H1) = dydt(H1) - 56.0e0_rt * ymol(Ni56) * ymol(H1) * rr(3).rates(r56eff);
    } else {
        dydt(Ni56) = dydt(Ni56) + ymol(Ni56) * ymol(H1) * rr(3).rates(dr56effdt);
        dydt(H1) = dydt(H1) - 56.0e0_rt * ymol(Ni56) * ymol(H1 ) * rr(3).rates(dr56effdt);
    }
}


template<class T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ener_gener_rate(T& dydt, Real& enuc)
{
    enuc = 0.0_rt;

    for (int i = 1; i <= NumSpec; ++i) {
        enuc -= dydt(i) * aion[i-1] * rprox::ebin(i);
    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_rhs (burn_t& state, Array1D<Real, 1, neqs>& ydot)
{
    for (int i = 1; i <= neqs; ++i) {
        ydot(i) = 0.0_rt;
    }

    Real dens = state.rho;
    Real T9   = state.T * 1.e-9_rt;

    Array1D<Real, 1, NumSpec> y;
    for (int i = 1; i <= NumSpec; ++i) {
        y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    // build the rates; weak rates are the wk* variables
    Array1D<rate_t, 1, Rates::NumGroups> rr;
    make_rates(T9, dens, y, state, rr);

    // set up the ODEs for the species
    make_ydots(y, T9, state, rr, ydot, false);

    // Energy release
    ener_gener_rate(ydot, ydot(net_ienuc));

    // Temperature ODE
    temperature_rhs(state, ydot);
}


template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_jac(burn_t& state, MatrixType& jac)
{
    using namespace Species;
    using namespace Rates;

    // initialize
    jac.zero();

    Array1D<Real, 1, NumSpec> ymol;
    for (int i = 1; i <= NumSpec; ++i) {
        ymol(i) = state.xn[i-1] * aion_inv[i-1];
    }
    Real T9 = state.T * 1.e-9_rt;

    Real dens = state.rho;

    // build the rates; weak rates are the wk* variables
    Array1D<rate_t, 1, Rates::NumGroups> rr;
    make_rates(T9, dens, ymol, state, rr);


    // carbon-12
    jac(C12, C12) = -ymol(H1) * rr(1).rates(irpg12c);
    jac(C12, O15) = rr(1).rates(irlambCNO);
    jac(C12, He4) = 0.5_rt * ymol(He4) * ymol(He4) * rr(1).rates(ir3a);
    jac(C12, H1)  = -ymol(C12) * rr(1).rates(irpg12c) + ymol(O15) * rr(3).rates(dlambCNOdh1);

    // oxygen-14
    jac(O14, C12) = ymol(H1) * rr(1).rates(irpg12c);
    jac(O14, O14) = -ymol(He4) * rr(1).rates(irap14o) - rr(1).rates(irwk14o);
    jac(O14, He4) = -ymol(O14) * rr(1).rates(irap14o);
    jac(O14, H1)  = ymol(C12) * rr(1).rates(irpg12c);

    // oxygen-15
    jac(O15, O14) = rr(1).rates(irwk14o);
    jac(O15, O15) = -ymol(He4) * rr(1).rates(irag15o) - rr(1).rates(irlambCNO);
    jac(O15, F17) = ymol(H1) * rr(1).rates(irpg17f) * rr(1).rates(irs1) + rr(1).rates(irwk17f);
    jac(O15, He4) = -ymol(O15) * rr(1).rates(irag15o) +
                    ymol(F17) * ymol(H1) * rr(1).rates(irpg17f) * rr(3).rates(drs1dhe4);
    jac(O15, H1)  = ymol(F17) * rr(1).rates(irpg17f) * rr(1).rates(irs1) -
                    ymol(O15) * rr(3).rates(dlambCNOdh1);

    // oxygen-16
    jac(O16, O15) = ymol(He4) * rr(1).rates(irr1) * rr(1).rates(irag15o);
    jac(O16, O16) = -ymol(H1) * rr(1).rates(irpg16o) -
                    ymol(He4) * rr(1).rates(irag16o);
    jac(O16, F17) = rr(1).rates(irgp17f);
    jac(O16, He4) = ymol(O15) * rr(1).rates(irr1) * rr(1).rates(irag15o) -
                    ymol(O16) * rr(1).rates(irag16o);
    jac(O16, H1)  = -ymol(O16) * rr(1).rates(irpg16o) +
                    ymol(O15) * ymol(He4) * rr(3).rates(drr1dh1) * rr(1).rates(irag15o);

    // flourine-17
    jac(F17, O14) = ymol(He4) * rr(1).rates(irap14o);
    jac(F17, O16) = ymol(H1) * rr(1).rates(irpg16o);
    jac(F17, F17) = -rr(1).rates(irgp17f) -
                    ymol(H1) * rr(1).rates(irpg17f) -
                    rr(1).rates(irwk17f);
    jac(F17, He4) = ymol(O14) * rr(1).rates(irap14o);
    jac(F17, H1)  = ymol(O16) * rr(1).rates(irpg16o) -
                    ymol(F17) * rr(1).rates(irpg17f);

    // magnesium-22
    jac(Mg22, O15) = ymol(He4) * rr(1).rates(irag15o) * (1.0_rt - rr(1).rates(irr1));
    jac(Mg22, O16) = ymol(He4) * rr(1).rates(irag16o);
    jac(Mg22, F17) = ymol(H1) * rr(1).rates(irpg17f) * (1.0_rt - rr(1).rates(irs1));
    jac(Mg22, Mg22) = -rr(1).rates(irlambda1);
    jac(Mg22, He4) = ymol(O16) * rr(1).rates(irag16o) +
                     ymol(O15) * rr(1).rates(irag15o) * (1.0_rt - rr(1).rates(irr1)) -
                     ymol(F17) * ymol(H1) * rr(1).rates(irpg17f) * rr(3).rates(drs1dhe4) -
                     ymol(Mg22) * rr(3).rates(dlambda1dhe4);
    jac(Mg22, H1)  = ymol(F17) * rr(1).rates(irpg17f) * (1.0_rt - rr(1).rates(irs1)) -
                     ymol(O15) * ymol(He4) * rr(1).rates(irag15o) * rr(3).rates(drr1dh1);

    // sulfur-30
    jac(S30, Mg22) = rr(1).rates(irlambda1);
    jac(S30, S30)  = -rr(1).rates(irlambda2);
    jac(S30, He4)  = ymol(Mg22) * rr(3).rates(dlambda1dhe4) -
                     ymol(S30) * rr(3).rates(dlambda2dhe4);

    // nickel-56
    jac(Ni56, S30) = rr(1).rates(irlambda2);
    jac(Ni56, Ni56) = ymol(H1) * rr(3).rates(r56eff);
    jac(Ni56, He4) = ymol(S30) * rr(3).rates(dlambda2dhe4);
    jac(Ni56, H1) = ymol(Ni56) * rr(3).rates(r56eff);

    // helium-4
    jac(He4, O14) = -ymol(He4) * rr(1).rates(irap14o);
    jac(He4, O15) = rr(1).rates(irlambCNO) -
                    ymol(He4) * rr(1).rates(irag15o) * (1.0_rt - rr(1).rates(irr1));
    jac(He4, O16) = -ymol(He4) * rr(1).rates(irag16o);
    jac(He4, F17) = ymol(H1) * rr(1).rates(irpg17f) * rr(1).rates(irs1) -
                    ymol(H1) * rr(1).rates(irpg17f) * (1.0_rt - rr(1).rates(irs1)) +
                    rr(1).rates(irwk17f);
    jac(He4, Mg22) = -2.0_rt * rr(1).rates(irlambda1) * rr(3).rates(delta1);
    jac(He4, S30) = -6.5e0_rt * rr(1).rates(irlambda2) * rr(3).rates(delta2);
    jac(He4, He4) = -3.0_rt * ymol(He4) * ymol(He4) * 0.5_rt * rr(1).rates(ir3a) -
                    ymol(O14) * rr(1).rates(irap14o) -
                    ymol(O16) * rr(1).rates(irag16o) -
                    ymol(O15) * rr(1).rates(irag15o) * (1.0_rt - rr(1).rates(irr1)) +
                    ymol(F17) * ymol(H1) * rr(1).rates(irpg17f) * rr(3).rates(drs1dhe4) +
                    ymol(F17) * ymol(H1) * rr(1).rates(irpg17f) * rr(3).rates(drs1dhe4) -
                    2.0_rt * ymol(Mg22) * rr(3).rates(dlambda1dhe4) * rr(3).rates(delta1) -
                    6.5e0_rt * ymol(S30) * rr(3).rates(dlambda2dhe4) * rr(3).rates(delta2);
    jac(He4, H1)  = ymol(F17) * rr(1).rates(irpg17f) * rr(1).rates(irs1) -
                    ymol(F17) * rr(1).rates(irpg17f) * (1.0_rt - rr(1).rates(irs1)) +
                    ymol(O15) * rr(3).rates(dlambCNOdh1) +
                    ymol(O15) * ymol(He4) * rr(1).rates(irag15o) * rr(3).rates(drr1dh1);

    // hydrogen-1
    jac(H1, C12) = -2.0_rt * ymol(H1) * rr(1).rates(irpg12c);
    jac(H1, O14) = ymol(He4) * rr(1).rates(irap14o) - rr(1).rates(irwk14o);
    jac(H1, O15) = -rr(1).rates(irlambCNO) -
                   ymol(He4) * rr(1).rates(irag15o) * rr(1).rates(irr1) -
                   3.0_rt * ymol(He4) * rr(1).rates(irag15o) * (1.0_rt - rr(1).rates(irr1));
    jac(H1, O16) = -ymol(H1) * rr(1).rates(irpg16o) - 2.0_rt * ymol(He4) * rr(1).rates(irag16o);
    jac(H1, F17) = -2.0_rt * ymol(H1) * rr(1).rates(irpg17f) * rr(1).rates(irs1) +
                   rr(1).rates(irgp17f) -
                   ymol(H1) * rr(1).rates(irpg17f) * (1.0_rt - rr(1).rates(irs1)) -
                   2.0_rt * rr(1).rates(irwk17f);
    jac(H1, Mg22) = -8.0_rt * rr(1).rates(irlambda1) * (1.0_rt - rr(3).rates(delta1));
    jac(H1, S30)  = -26.e0_rt * rr(1).rates(irlambda2) * (1.0_rt - rr(3).rates(delta2));
    jac(H1, Ni56) = -56.0e0_rt * ymol(H1) * rr(3).rates(r56eff);
    jac(H1, He4) = ymol(O14) * rr(1).rates(irap14o) -
                   ymol(O15) * rr(1).rates(irag15o) * rr(1).rates(irr1) -
                   2.0_rt * ymol(O16) * rr(1).rates(irag16o) -
                   3.0_rt * ymol(O15) * rr(1).rates(irag15o) * (1.0_rt - rr(1).rates(irr1)) -
                   ymol(F17) * ymol(H1) * rr(1).rates(irpg17f) * rr(3).rates(drs1dhe4) -
                   8.0_rt * ymol(Mg22) * rr(3).rates(dlambda1dhe4) * (1.0_rt - rr(3).rates(delta1)) -
                   26.e0_rt * ymol(S30) * rr(3).rates(dlambda2dhe4) * (1.0_rt - rr(3).rates(delta2));
    jac(H1, H1)  = -2.0_rt * ymol(C12) * rr(1).rates(irpg12c) -
                   2.0_rt * ymol(F17) * rr(1).rates(irpg17f) * rr(1).rates(irs1) -
                   ymol(O16) * rr(1).rates(irpg16o) -
                   ymol(F17) * rr(1).rates(irpg17f) * (1.0_rt - rr(1).rates(irs1)) -
                   ymol(O15) * rr(3).rates(dlambCNOdh1) +
                   2.0_rt * ymol(O15) * ymol(He4) * rr(1).rates(irag15o) * rr(3).rates(drr1dh1) -
                   56.0e0_rt * ymol(Ni56) * rr(3).rates(r56eff);

    // temperature derivatives df(Y)/df(T)
    Array1D<Real, 1, neqs> ydot;
    make_ydots(ymol, T9, state, rr, ydot, true);

    for (int i = 1; i <= NumSpec; ++i) {
        jac(i, net_itemp) = ydot(i);
    }

    // Energy generation rate Jacobian elements with respect to species

    for (int j = 1; j <= NumSpec; ++j) {
        auto jac_slice_2 = [&](int i) -> Real { return jac.get(i, j); };
        ener_gener_rate(jac_slice_2, jac(net_ienuc, j));
    }

    // Jacobian elements with respect to temperature

    ener_gener_rate(ydot, jac(net_ienuc, net_itemp));

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


    // 0: C12(p,g)N13
    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[H1-1], aion[H1-1]);

    // 1, 2: 3-alpha
    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[He4-1], aion[He4-1]);
    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], 4.0_rt, 8.0_rt);

    // 3: F17(p,g)Ne18
    add_screening_factor(jscr++, zion[F17-1], aion[F17-1], zion[H1-1], aion[H1-1]);

    // F17(g,p)O16 is unscreened

    // 4: O15(a,g)Ne19
    add_screening_factor(jscr++, zion[O15-1], aion[O15-1], zion[He4-1], aion[He4-1]);

    // 5: O16(a,g)Ne20
    add_screening_factor(jscr++, zion[O16-1], aion[O16-1], zion[He4-1], aion[He4-1]);

    // 6: O16(p,g)F17
    add_screening_factor(jscr++, zion[O16-1], aion[O16-1], zion[H1-1], aion[H1-1]);

    // 7: O14(a,p)F17
    add_screening_factor(jscr++, zion[O14-1], aion[O14-1], zion[He4-1], aion[He4-1]);

    // 8: N14(p,g)O15
    add_screening_factor(jscr++, 7.0_rt, 14.0_rt, zion[H1-1], aion[H1-1]);

    // 9: Si26(a,p)P29
    add_screening_factor(jscr++, 14.0_rt, 26.0_rt, zion[He4-1], aion[He4-1]);

    // 10: Ti44(a,p)V47
    add_screening_factor(jscr++, 22.0_rt, 44.0_rt, zion[He4-1], aion[He4-1]);

    // 11: Ne18(a,p)Na21
    add_screening_factor(jscr++, 10.0_rt, 18.0e0_rt, zion[He4-1], aion[He4-1]);

    // 12: Ne19(p,g)Na20
    add_screening_factor(jscr++, 10.0_rt, 19.0_rt, zion[H1-1], aion[H1-1]);

}

AMREX_INLINE
void actual_rhs_init () {

    screening_init();

    set_up_screening_factors();

}



#endif
