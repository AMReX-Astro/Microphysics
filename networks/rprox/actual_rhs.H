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

using namespace amrex;
using namespace ArrayUtil;

AMREX_INLINE
void actual_rhs_init () {}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void make_rates (const Real& T9, const Real& dens, Array1D<Real, 1, NumSpec>& y,
                 burn_t& state, rate_t& rr)
{
    using namespace Species;
    using namespace Rates;

    // locally used rates
    Real rate, dratedt, wk18ne, wk19ne;
    Real r56pg, dr56pgdt, cutoni, dcutonidt, r57decay, dr56eff, ddr56effdt;
    Real T9i32;

    // some numbers from appendix C in WW81; these should probably be
    // updated with current rates
    constexpr Real Lweak = 1.05e0_rt; // this is for NS
    // constexpr Real Lweak = 0.107e0_rt; // this is for lower densities

    constexpr Real la2 = 1.0_rt / 15.0_rt; // mean rate from 30s to 56ni from p-capture and beta decays

    for (int j = 1; j <= NumGroups; ++j) {
        for (int i = 1; i <= NumRates; ++i) {
            rr.rates(j, i) = 0.0_rt; // Zero out rates
        }
    }

    rr.T_eval = T9 * 1.e9_rt;

    tf_t tfactors = calc_tfactors(T9);

    // some common parameters
    rr.rates(1, irLweak) = Lweak;
    rr.rates(1, irla2)   = la2;

    // weak rates first
    //
    // 14o(beta nu)14n
    rate_o14_to_n14(tfactors, rate, dratedt);
    rr.rates(1, irwk14o) = rate;

    // 15o(beta nu)15n
    rate_o15_to_n15(tfactors, rate, dratedt);
    rr.rates(1, irwk15o) = rate;

    // 17f(beta nu)17o
    rate_f17_to_o17(tfactors, rate, dratedt);
    rr.rates(1, irwk17f) = rate;

    // these weak rates aren't needed outside of this routine
    // 18ne(beta nu)18f
    rate_ne18_to_f18(tfactors, wk18ne, dratedt);
    // 19ne(beta nu)19f
    rate_ne19_to_f19(tfactors, wk19ne, dratedt);

    // 12c(p,g)13n
    rate_p_c12_to_n13(tfactors, rate, dratedt);
    rr.rates(1, irpg12c) = dens * rate;
    rr.rates(2, irpg12c) = dens * dratedt;

    // triple alpha
    rate_he4_he4_he4_to_c12(tfactors, rate, dratedt);
    rr.rates(1, ir3a) = dens * dens * rate;
    rr.rates(2, ir3a) = dens * dens * dratedt;

    // 17f(p,g)18ne
    rate_p_f17_to_ne18(tfactors, rate, dratedt);
    rr.rates(1, irpg17f) = dens * rate;
    rr.rates(2, irpg17f) = dens * dratedt;

    // 17f(g,p)16o
    rate_f17_to_p_o16(tfactors, rate, dratedt);
    rr.rates(1, irgp17f) = rate;
    rr.rates(2, irgp17f) = dratedt;

    // 15o(a,g)19ne
    rate_he4_o15_to_ne19(tfactors, rate, dratedt);
    rr.rates(1, irag15o) = dens * rate;
    rr.rates(2, irag15o) = dens * dratedt;

    // 16o(a,g)20ne
    rate_he4_o16_to_ne20(tfactors, rate, dratedt);
    rr.rates(1, irag16o) = dens * rate;
    rr.rates(2, irag16o) = dens * dratedt;

    // 16o(p,g)17f
    rate_p_o16_to_f17(tfactors, rate, dratedt);
    rr.rates(1, irpg16o) = dens * rate;
    rr.rates(2, irpg16o) = dens * dratedt;

    // 14o(a,p)17f
    rate_he4_o14_to_p_f17(tfactors, rate, dratedt);
    rr.rates(1, irap14o) = dens * rate;
    rr.rates(2, irap14o) = dens * dratedt;

    // limit CNO as minimum between 14n(p,g)15o and 15o(beta nu)15n
    // we store the limited rate in irlambCNO; this is lambda_CNO in WW81
    rate_p_n14_to_o15(tfactors, rate, dratedt);
    rr.rates(1, irlambCNO) = amrex::min(rr.rates(1, irwk15o), rate * dens * y(H1));
    if (rr.rates(1, irlambCNO) < rr.rates(1,irwk15o)) {
        rr.rates(2, irlambCNO) = dens * y(H1) * dratedt;
        rr.rates(3, dlambCNOdh1) = rate * dens;
    }

    // 22mg(...)30s
    // check if this proceeds via p-captures or (a,p) reactions
    // the Lweak is from WW81, eqn C15
    // we store the rate in irlambda1; this is the lambda1 in WW81
    rate_he4_si26_to_p_p29(tfactors, rate, dratedt);
    rr.rates(1, irlambda1) = amrex::max(rr.rates(1, irLweak), dens * y(He4) * rate);
    if (rr.rates(1, irlambda1) > rr.rates(1,irLweak)) {
        rr.rates(2, irlambda1) = dens * y(He4) * dratedt;
        rr.rates(3, dlambda1dhe4) = dens * rate;
        // use the sign of state.rates(1, irlambda1) to indicate the value of delta1 in WW81
        // if delta1 = 1, then we multiply the rate by -1
        rr.rates(1, irlambda1) = -1.0_rt * rr.rates(1, irlambda1);
    }

    // 30s(...) 56ni
    // check if this proceeds via p-captures or (a,p) reactions
    // use 44ti(a,p)v47 as a typical limiting rate for the (a,p) process
    // store this in irlambda2; this is lambda2 in WW81
    rate_he4_ti44_to_p_v47(tfactors, rate, dratedt);
    rr.rates(1, irlambda2) = amrex::max(rr.rates(1, irla2), dens * y(He4) * rate);
    if (rr.rates(1, irlambda2) > rr.rates(1, irla2)) {
        rr.rates(2, irlambda2) = dens * y(He4) * dratedt;
        rr.rates(3, dlambda2dhe4) = dens * rate;
        // use the sign of rr.rates(1, irlambda2) to indicate the value of delta2
        // if delta2 = 1, then we multiply the rate by -1
        rr.rates(1, irlambda2) = -1.0_rt * rr.rates(1, irlambda2);
    }

    // form s1 from WW81; branching ratio for 18ne beta decay (wk18ne) vs (a,p)
    // store result in irs1
    // 18ne(a,p)21na
    rate_he4_ne18_to_p_na21(tfactors, rate, dratedt);
    rr.rates(1, irs1) = wk18ne / (wk18ne + dens * y(He4) * rate);
    rr.rates(2, irs1) = -rr.rates(1, irs1) * dens * y(He4) * dratedt 
                         / (wk18ne + dens * y(He4) * rate);
    rr.rates(3, drs1dhe4) = -rr.rates(1, irs1) * dens * rate
                             / (wk18ne + dens * y(He4) * rate);

    // form r1 from WW81; ranching ratio for 19ne beta decay (wk19ne) vs (p,g)
    // store result in irr1
    // 19ne(p,g)20na
    rate_p_ne19_to_na20(tfactors, rate, dratedt);
    rr.rates(1, irr1) = wk19ne / (wk19ne + dens * y(H1) * rate);
    rr.rates(2, irr1) = -rr.rates(1,irr1) * dens * y(H1) * dratedt
                         / (wk19ne + dens * y(H1) * rate);
    rr.rates(3, drr1dh1) = -rr.rates(1, irr1) * dens * rate
                            / (wk19ne + dens * y(H1) * rate);


    //....
    //....  additional coding for proton capture on 56ni to heavier elements
    //....   kludge    56ni+56p -> 2 (56ni) at a rate given by min
    //....   of 56ni(pg) and 57cu decay rate
    //....
    //....  use 56ni rate from wallace and woosley 1981
    T9i32 = tfactors.T9i * std::sqrt(tfactors.T9i);
    r56pg = dens * (1.29e-02_rt * std::exp(-4.897_rt * tfactors.T9i)
            + 7.065e+03_rt * std::exp(-20.33_rt * tfactors.T9i)) * T9i32;
    dr56pgdt = -(3.0_rt / 2.0_rt) * r56pg * tfactors.T9i +
               dens * T9i32 * tfactors.T9i*tfactors.T9i *
               (4.897_rt * 1.29e-2_rt * std::exp(-4.897_rt * tfactors.T9i)
                + 20.33_rt * 7.065e3_rt * std::exp(-20.33_rt * tfactors.T9i));
    //....  use generic proton separation energy of 400 kev
    //....  8.02 -> 4.64
    //      cutoni=2.08e-10_rt*dens*exp(8.02*T9m1)/T932
    cutoni = 2.08e-10_rt * dens * std::exp(4.642_rt * tfactors.T9i) * T9i32;
    dcutonidt = cutoni * tfactors.T9i * (-(3.0_rt / 2.0_rt) - 4.642_rt * tfactors.T9i);
    r57decay = 3.54_rt;
    dr56eff = amrex::min(r56pg, cutoni * r57decay);
    //   rr.rates(3, r56eff) = d56eff;
    //   if (d56eff < r56pg) rr.rates(3, dr56effdt) = r57decay * dcutonidt;
    rr.rates(3, r56eff) = 0.0_rt;
    rr.rates(3, dr56effdt) = 0.0_rt;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void make_ydots (const Array1D<Real, 1, NumSpec>& ymol, const Real& T9,
                 burn_t& state, rate_t& rr,
                 Array1D<Real, 1, neqs>& dydt, bool doing_dratesdt)
{
    using namespace Species;
    using namespace Rates;

    // initialize
    for (int i = 1; i <= NumSpec; ++i) {
        dydt(i) = 0.0_rt;
    }

    Real dens = state.rho;

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
        if (rr.rates(1, irlambda1) < 0.0_rt) {
            ddelta1 = 1.0_rt;
            rr.rates(1, irlambda1) = -1.0_rt * rr.rates(1, irlambda1);
        }
        if (rr.rates(1, irlambda2) < 0.0_rt) {
            ddelta2 = 1.0_rt;
            rr.rates(1, irlambda2) = -1.0_rt * rr.rates(1, irlambda2);
        }
        rr.rates(3, delta1) = ddelta1;
        rr.rates(3, delta2) = ddelta2;
    }

    // setup ODEs
    //
    //....
    //.... 12c = 1
    //....
    dydt(C12) = -ymol(C12) * ymol(H1) * rr.rates(rate_idx,irpg12c) +
                 ymol(He4) * ymol(He4) * ymol(He4) * rr.rates(rate_idx, ir3a) / 6.0_rt +
                 ymol(O15) * rr.rates(rate_idx, irlambCNO);
    //....
    //.... 14o = 2
    //....
    dydt(O14) = -ymol(O14) * ymol(He4) * rr.rates(rate_idx, irap14o) -
                 ymol(O14) * rr.rates(rate_idx, irwk14o) +
                 ymol(C12) * ymol(H1) * rr.rates(rate_idx, irpg12c);
    //....
    //.... 15o = 3
    //....
    dydt(O15) = ymol(O14) * rr.rates(rate_idx, irwk14o) -
                ymol(O15) * ymol(He4) * rr.rates(rate_idx, irag15o) -
                ymol(O15) * rr.rates(rate_idx, irlambCNO) +
                ymol(F17) * ymol(H1) * rr.rates(rate_idx, irpg17f) * rr.rates(rate_idx, irs1) +
                ymol(F17) * rr.rates(rate_idx, irwk17f);
    //....
    //.... 16o = 4
    //....
    dydt(O16) = ymol(F17) * rr.rates(rate_idx, irgp17f) -
                ymol(O16) * ymol(H1) * rr.rates(rate_idx, irpg16o) +
                ymol(O15) * ymol(He4) * rr.rates(rate_idx, irr1) * rr.rates(rate_idx, irag15o) -
                ymol(O16) * ymol(He4) * rr.rates(rate_idx, irag16o);
    //....
    //.... 17f = 5
    //....
    dydt(F17) = ymol(O14) * ymol(He4) * rr.rates(rate_idx, irap14o) +
                 ymol(O16) * ymol(H1) * rr.rates(rate_idx, irpg16o) -
                 ymol(F17) * rr.rates(rate_idx, irgp17f) -
                 ymol(F17) * ymol(H1) * rr.rates(rate_idx, irpg17f) -
                 ymol(F17) * rr.rates(rate_idx, irwk17f);
    //....
    //.... 22mg = 6
    //....
    dydt(Mg22) = ymol(O16) * ymol(He4) * rr.rates(rate_idx, irag16o) +
                  ymol(F17) * ymol(H1) * rr.rates(rate_idx, irpg17f) * (1.0_rt - rr.rates(rate_idx, irs1)) +
                  ymol(O15) * ymol(He4) * rr.rates(rate_idx, irag15o) * (1.0_rt - rr.rates(rate_idx, irr1)) -
                  ymol(Mg22) * rr.rates(rate_idx, irlambda1);
    //....
    //.... 30s = 7
    //....
    dydt(S30) = ymol(Mg22) * rr.rates(rate_idx, irlambda1) - ymol(S30) * rr.rates(rate_idx, irlambda2);
    //....
    //.... amax (56ni) = 8  (note that WW81 have a typo -- they write lambda1 here)
    //....
    dydt(Ni56) = ymol(S30) * rr.rates(rate_idx, irlambda2);
    //....
    //.... 4he (alpha) = 9
    //....
    dydt(He4) = -ymol(He4) * ymol(He4) * ymol(He4) * 0.5_rt * rr.rates(rate_idx, ir3a) +
                 ymol(O15) * rr.rates(rate_idx, irlambCNO) -
                 ymol(O14) * ymol(He4) * rr.rates(rate_idx, irap14o) +
                 ymol(F17) * ymol(H1) * rr.rates(rate_idx, irpg17f) * rr.rates(rate_idx, irs1) -
                 ymol(O15) * ymol(He4) * rr.rates(rate_idx, irag15o) * (1.0_rt - rr.rates(rate_idx, irr1)) -
                 ymol(O16) * ymol(He4) * rr.rates(rate_idx, irag16o) -
                 ymol(F17) * ymol(H1) * rr.rates(rate_idx, irpg17f) * (1.0_rt - rr.rates(rate_idx, irs1)) +
                 ymol(F17) * rr.rates(rate_idx, irwk17f) -
                 2.0_rt * ymol(Mg22) * rr.rates(rate_idx, irlambda1) * rr.rates(3, delta1) -
                 6.5e0_rt * ymol(S30) * rr.rates(rate_idx, irlambda2) * rr.rates(3, delta2);
    //....
    //.... 1h (p) = 10
    //....
    dydt(H1) = -ymol(O14) * rr.rates(rate_idx, irwk14o) -
                ymol(O15) * rr.rates(rate_idx, irlambCNO) -
                2.0_rt * ymol(C12) * ymol(H1) * rr.rates(rate_idx, irpg12c) +
                ymol(O14) * ymol(He4) * rr.rates(rate_idx, irap14o) -
                2.0_rt * ymol(F17) * ymol(H1) * rr.rates(rate_idx, irpg17f) * rr.rates(rate_idx, irs1) +
                ymol(F17) * rr.rates(rate_idx, irgp17f) -
                ymol(O16) * ymol(H1) * rr.rates(rate_idx, irpg16o) -
                ymol(O15) * ymol(He4) * rr.rates(rate_idx, irag15o) * rr.rates(rate_idx, irr1) -
                2.0_rt * ymol(O16) * ymol(He4) * rr.rates(rate_idx, irag16o) -
                3.0_rt * ymol(O15) * ymol(He4) * rr.rates(rate_idx, irag15o) * (1.0_rt - rr.rates(rate_idx, irr1)) -
                ymol(F17) * ymol(H1) * rr.rates(rate_idx, irpg17f) * (1.0_rt - rr.rates(rate_idx, irs1)) -
                2.0_rt * ymol(F17) * rr.rates(rate_idx, irwk17f) -
                8.0_rt * ymol(Mg22) * rr.rates(rate_idx, irlambda1) * (1.0_rt - rr.rates(3, delta1)) -
                26.e0_rt * ymol(S30) * rr.rates(rate_idx, irlambda2) * (1.0_rt - rr.rates(3, delta2));


    if (!doing_dratesdt) {
        dydt(Ni56) = dydt(Ni56) + ymol(Ni56) * ymol(H1) * rr.rates(3, r56eff);
        dydt(H1) = dydt(H1) - 56.0e0_rt * ymol(Ni56) * ymol(H1) * rr.rates(3, r56eff);
    } else {
        dydt(Ni56) = dydt(Ni56) + ymol(Ni56) * ymol(H1) * rr.rates(3, dr56effdt);
        dydt(H1) = dydt(H1) - 56.0e0_rt * ymol(Ni56) * ymol(H1 ) * rr.rates(3, dr56effdt);
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
    rate_t rr;
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
    rate_t rr;
    make_rates(T9, dens, ymol, state, rr);


    // carbon-12
    jac(C12, C12) = -ymol(H1) * rr.rates(1, irpg12c);
    jac(C12, O15) = rr.rates(1, irlambCNO);
    jac(C12, He4) = 0.5_rt * ymol(He4) * ymol(He4) * rr.rates(1, ir3a);
    jac(C12, H1)  = -ymol(C12) * rr.rates(1, irpg12c) + ymol(O15) * rr.rates(3, dlambCNOdh1);

    // oxygen-14
    jac(O14, C12) = ymol(H1) * rr.rates(1, irpg12c);
    jac(O14, O14) = -ymol(He4) * rr.rates(1, irap14o) - rr.rates(1, irwk14o);
    jac(O14, He4) = -ymol(O14) * rr.rates(1, irap14o);
    jac(O14, H1)  = ymol(C12) * rr.rates(1, irpg12c);

    // oxygen-15
    jac(O15, O14) = rr.rates(1, irwk14o);
    jac(O15, O15) = -ymol(He4) * rr.rates(1, irag15o) - rr.rates(1, irlambCNO);
    jac(O15, F17) = ymol(H1) * rr.rates(1, irpg17f) * rr.rates(1, irs1) + rr.rates(1, irwk17f);
    jac(O15, He4) = -ymol(O15) * rr.rates(1, irag15o) +
                    ymol(F17) * ymol(H1) * rr.rates(1, irpg17f) * rr.rates(3, drs1dhe4);
    jac(O15, H1)  = ymol(F17) * rr.rates(1, irpg17f) * rr.rates(1, irs1) -
                    ymol(O15) * rr.rates(3, dlambCNOdh1);

    // oxygen-16
    jac(O16, O15) = ymol(He4) * rr.rates(1, irr1) * rr.rates(1, irag15o);
    jac(O16, O16) = -ymol(H1) * rr.rates(1, irpg16o) -
                    ymol(He4) * rr.rates(1, irag16o);
    jac(O16, F17) = rr.rates(1, irgp17f);
    jac(O16, He4) = ymol(O15) * rr.rates(1, irr1) * rr.rates(1, irag15o) -
                    ymol(O16) * rr.rates(1, irag16o);
    jac(O16, H1)  = -ymol(O16) * rr.rates(1, irpg16o) +
                    ymol(O15) * ymol(He4) * rr.rates(3, drr1dh1) * rr.rates(1, irag15o);

    // flourine-17
    jac(F17, O14) = ymol(He4) * rr.rates(1, irap14o);
    jac(F17, O16) = ymol(H1) * rr.rates(1, irpg16o);
    jac(F17, F17) = -rr.rates(1, irgp17f) -
                    ymol(H1) * rr.rates(1, irpg17f) -
                    rr.rates(1, irwk17f);
    jac(F17, He4) = ymol(O14) * rr.rates(1, irap14o);
    jac(F17, H1)  = ymol(O16) * rr.rates(1, irpg16o) -
                    ymol(F17) * rr.rates(1, irpg17f);

    // magnesium-22
    jac(Mg22, O15) = ymol(He4) * rr.rates(1, irag15o) * (1.0_rt - rr.rates(1, irr1));
    jac(Mg22, O16) = ymol(He4) * rr.rates(1, irag16o);
    jac(Mg22, F17) = ymol(H1) * rr.rates(1, irpg17f) * (1.0_rt - rr.rates(1, irs1));
    jac(Mg22, Mg22) = -rr.rates(1, irlambda1);
    jac(Mg22, He4) = ymol(O16) * rr.rates(1, irag16o) +
                     ymol(O15) * rr.rates(1, irag15o) * (1.0_rt - rr.rates(1, irr1)) -
                     ymol(F17) * ymol(H1) * rr.rates(1, irpg17f) * rr.rates(3, drs1dhe4) -
                     ymol(Mg22) * rr.rates(3, dlambda1dhe4);
    jac(Mg22, H1)  = ymol(F17) * rr.rates(1, irpg17f) * (1.0_rt - rr.rates(1, irs1)) -
                     ymol(O15) * ymol(He4) * rr.rates(1, irag15o) * rr.rates(3, drr1dh1);

    // sulfur-30
    jac(S30, Mg22) = rr.rates(1, irlambda1);
    jac(S30, S30)  = -rr.rates(1, irlambda2);
    jac(S30, He4)  = ymol(Mg22) * rr.rates(3, dlambda1dhe4) -
                     ymol(S30) * rr.rates(3, dlambda2dhe4);

    // nickel-56
    jac(Ni56, S30) = rr.rates(1, irlambda2);
    jac(Ni56, Ni56) = ymol(H1) * rr.rates(3, r56eff);
    jac(Ni56, He4) = ymol(S30) * rr.rates(3, dlambda2dhe4);
    jac(Ni56, H1) = ymol(Ni56) * rr.rates(3, r56eff);

    // helium-4
    jac(He4, O14) = -ymol(He4) * rr.rates(1, irap14o);
    jac(He4, O15) = rr.rates(1, irlambCNO) -
                    ymol(He4) * rr.rates(1, irag15o) * (1.0_rt - rr.rates(1, irr1));
    jac(He4, O16) = -ymol(He4) * rr.rates(1, irag16o);
    jac(He4, F17) = ymol(H1) * rr.rates(1, irpg17f) * rr.rates(1, irs1) -
                    ymol(H1) * rr.rates(1, irpg17f) * (1.0_rt - rr.rates(1, irs1)) +
                    rr.rates(1, irwk17f);
    jac(He4, Mg22) = -2.0_rt * rr.rates(1, irlambda1) * rr.rates(3, delta1);
    jac(He4, S30) = -6.5e0_rt * rr.rates(1, irlambda2) * rr.rates(3, delta2);
    jac(He4, He4) = -3.0_rt * ymol(He4) * ymol(He4) * 0.5_rt * rr.rates(1, ir3a) -
                    ymol(O14) * rr.rates(1, irap14o) -
                    ymol(O16) * rr.rates(1, irag16o) -
                    ymol(O15) * rr.rates(1, irag15o) * (1.0_rt - rr.rates(1, irr1)) +
                    ymol(F17) * ymol(H1) * rr.rates(1, irpg17f) * rr.rates(3, drs1dhe4) +
                    ymol(F17) * ymol(H1) * rr.rates(1, irpg17f) * rr.rates(3, drs1dhe4) -
                    2.0_rt * ymol(Mg22) * rr.rates(3, dlambda1dhe4) * rr.rates(3, delta1) -
                    6.5e0_rt * ymol(S30) * rr.rates(3, dlambda2dhe4) * rr.rates(3, delta2);
    jac(He4, H1)  = ymol(F17) * rr.rates(1, irpg17f) * rr.rates(1, irs1) -
                    ymol(F17) * rr.rates(1, irpg17f) * (1.0_rt - rr.rates(1, irs1)) +
                    ymol(O15) * rr.rates(3, dlambCNOdh1) +
                    ymol(O15) * ymol(He4) * rr.rates(1, irag15o) * rr.rates(3, drr1dh1);

    // hydrogen-1
    jac(H1, C12) = -2.0_rt * ymol(H1) * rr.rates(1, irpg12c);
    jac(H1, O14) = ymol(He4) * rr.rates(1, irap14o) - rr.rates(1, irwk14o);
    jac(H1, O15) = -rr.rates(1, irlambCNO) -
                   ymol(He4) * rr.rates(1, irag15o) * rr.rates(1, irr1) -
                   3.0_rt * ymol(He4) * rr.rates(1, irag15o) * (1.0_rt - rr.rates(1, irr1));
    jac(H1, O16) = -ymol(H1) * rr.rates(1, irpg16o) - 2.0_rt * ymol(He4) * rr.rates(1, irag16o);
    jac(H1, F17) = -2.0_rt * ymol(H1) * rr.rates(1, irpg17f) * rr.rates(1, irs1) +
                   rr.rates(1, irgp17f) -
                   ymol(H1) * rr.rates(1, irpg17f) * (1.0_rt - rr.rates(1, irs1)) -
                   2.0_rt * rr.rates(1, irwk17f);
    jac(H1, Mg22) = -8.0_rt * rr.rates(1, irlambda1) * (1.0_rt - rr.rates(3, delta1));
    jac(H1, S30)  = -26.e0_rt * rr.rates(1, irlambda2) * (1.0_rt - rr.rates(3, delta2));
    jac(H1, Ni56) = -56.0e0_rt * ymol(H1) * rr.rates(3, r56eff);
    jac(H1, He4) = ymol(O14) * rr.rates(1, irap14o) -
                   ymol(O15) * rr.rates(1, irag15o) * rr.rates(1, irr1) -
                   2.0_rt * ymol(O16) * rr.rates(1, irag16o) -
                   3.0_rt * ymol(O15) * rr.rates(1, irag15o) * (1.0_rt - rr.rates(1, irr1)) -
                   ymol(F17) * ymol(H1) * rr.rates(1, irpg17f) * rr.rates(3, drs1dhe4) -
                   8.0_rt * ymol(Mg22) * rr.rates(3, dlambda1dhe4) * (1.0_rt - rr.rates(3, delta1)) -
                   26.e0_rt * ymol(S30) * rr.rates(3, dlambda2dhe4) * (1.0_rt - rr.rates(3, delta2));
    jac(H1, H1)  = -2.0_rt * ymol(C12) * rr.rates(1, irpg12c) -
                   2.0_rt * ymol(F17) * rr.rates(1, irpg17f) * rr.rates(1, irs1) -
                   ymol(O16) * rr.rates(1, irpg16o) -
                   ymol(F17) * rr.rates(1, irpg17f) * (1.0_rt - rr.rates(1, irs1)) -
                   ymol(O15) * rr.rates(3, dlambCNOdh1) +
                   2.0_rt * ymol(O15) * ymol(He4) * rr.rates(1, irag15o) * rr.rates(3, drr1dh1) -
                   56.0e0_rt * ymol(Ni56) * rr.rates(3, r56eff);

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

#endif
