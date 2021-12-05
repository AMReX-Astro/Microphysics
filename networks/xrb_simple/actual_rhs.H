#ifndef actual_rhs_H
#define actual_rhs_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <extern_parameters.H>
#include <actual_network.H>
#include <rate_type.H>
#include <jacobian_utilities.H>
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


    for (int j = 1; j <= NumGroups; ++j) {
        for (int i = 1; i <= NumRates; ++i) {
            rr(j).rates(i) = 0.0_rt; // Zero out rates
        }
    }

    tf_t tfactors = calc_tfactors(T9);

    // locally used rates
    Real rate, dratedt;

    // TODO - temperature derivatives for use in analytic jacobian

    // helium burning - has been divided by 6
    rr(1).rates(ir3a) = 2.79e-8_rt * tfactors.t9m3 *
        std::exp(-4.403_rt * tfactors.t9m1) * dens * dens / 6.0_rt;

    // 15o(ag)19ne
    rr(1).rates(irag15) = (19.0_rt * std::pow(T9, 2.85_rt) *
                           std::exp(-7.356_rt * tfactors.t9m1) +
                           0.395_rt * tfactors.t9m32 *
                           std::exp(-5.849_rt * tfactors.t9m1)) * dens;

    ! 14o(ap)17f
    ! an approx good below t9=0.5  is to drop first term, then
    !  rap14 = 3.31e+04_rt*t9m32*exp(-11.733*t9m1)
    !   1  + 1.79e+07_rt*t9m32*exp(-22.609/t9)+9000.*t9113*exp(-12.517/t9)
    rr % rates(1,irap14) = 1.68e13_rt*t9m23 * &
         exp(-39.388e0_rt*t9m13-(t9/0.717e0_rt)**2) &
         *(ONE+0.011e0_rt*t913+13.117e0_rt*t923+0.971e0_rt*t9+85.295e0_rt*t943 &
         +16.06e0_rt*t943)+3.31e4_rt*t9m32*exp(-11.733e0_rt*t9m1) &
         +1.79e7_rt*t9m32*exp(-22.609e0_rt*t9m1) &
         +9.00e+03_rt*t9113*exp(-12.517e0_rt*t9m1)
    rr % rates(1,irap14) = dens* rr % rates(1,irap14)

    ! 18ne(ap)21na
    rr % rates(1,irap18) = exp(56.59577e0_rt-2.006856e0_rt*t9m1 &
         +26.05828e0_rt*t9m13 &
         -87.88732e0_rt*t913 + 3.718809e0_rt*t9 - 0.1767444e0_rt*t953 &
         + 46.971960e0_rt*t9log)
    rr % rates(1,irap18) = dens * rr % rates(1,irap18)

    ! weak rates
    rr % rates(1,irwk14o) = wk14o
    rr % rates(1,irwk15o) = wk15o

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
        enuc -= dydt(i) * aion[i-1] * network::ebin(i);
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

    // energy derivatives df(Y)/df(e)
    Array1D<Real, 1, neqs> ydot;
    make_ydots(ymol, T9, state, rr, ydot, true);

    for (int i = 1; i <= NumSpec; ++i) {
        jac(i, net_ienuc) = temperature_to_energy_jacobian(state, ydot(i));
    }

    // Energy generation rate Jacobian elements with respect to species

    for (int j = 1; j <= NumSpec; ++j) {
        auto jac_slice_2 = [&](int i) -> Real { return jac.get(i, j); };
        ener_gener_rate(jac_slice_2, jac(net_ienuc, j));
    }

    // Energy generation rate Jacobian element with respect to energy
    Real jac_e_T;
    ener_gener_rate(ydot, jac_e_T);
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
