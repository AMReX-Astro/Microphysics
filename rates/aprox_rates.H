#ifndef APROX_RATES_H
#define APROX_RATES_H

#include <aprox_rates_data.H>
#include <extern_parameters.H>
#include <microphysics_autodiff.H>
#include <tfactors.H>
#include <cmath>

using namespace network_rp;

inline
void rates_init()
{
    // store rates
    const amrex::Vector<amrex::Real> datn0 {-4.363_rt, -3.091_rt, -1.275_rt, 1.073_rt, 3.035_rt, 4.825_rt,
        -4.17_rt, -2.964_rt, -1.177_rt, 1.085_rt, 3.037_rt, 4.826_rt,
        -3.834_rt, -2.727_rt, -1.039_rt, 1.104_rt, 3.04_rt, 4.826_rt,
        -3.284_rt, -2.418_rt, -0.882_rt, 1.129_rt, 3.043_rt, 4.827_rt,
        -2.691_rt, -2.093_rt, -0.719_rt, 1.159_rt, 3.048_rt, 4.827_rt,
        -2.1675_rt, -1.7668_rt, -0.5573_rt, 1.1947_rt, 3.0527_rt, 4.8272_rt,
        -1.7095_rt, -1.4462_rt, -0.3991_rt, 1.2358_rt, 3.0577_rt, 4.8276_rt,
        -1.3119_rt, -1.1451_rt, -0.2495_rt, 1.2818_rt, 3.0648_rt, 4.8284_rt,
        -0.9812_rt, -0.8612_rt, -0.1084_rt, 1.3336_rt, 3.0738_rt, 4.8295_rt,
        -0.682_rt, -0.595_rt, 0.028_rt, 1.386_rt, 3.084_rt, 4.831_rt,
        -0.4046_rt, -0.3523_rt, 0.1605_rt, 1.4364_rt, 3.0957_rt, 4.8333_rt,
        -0.1636_rt, -0.1352_rt, 0.2879_rt, 1.4861_rt, 3.1092_rt, 4.8365_rt,
        0.0461_rt, 0.0595_rt, 0.4105_rt, 1.5354_rt, 3.1242_rt, 4.8405_rt,
        0.2295_rt, 0.235_rt, 0.5289_rt, 1.5842_rt, 3.1405_rt, 4.845_rt};
    const amrex::Vector<amrex::Real> datn1 {-4.539_rt, -3.097_rt, -1.134_rt, 1.525_rt, 3.907_rt, 6.078_rt,
        -4.199_rt, -2.905_rt, -1.024_rt, 1.545_rt, 3.91_rt, 6.079_rt,
        -3.736_rt, -2.602_rt, -0.851_rt, 1.578_rt, 3.916_rt, 6.08_rt,
        -3.052_rt, -2.206_rt, -0.636_rt, 1.623_rt, 3.923_rt, 6.081_rt,
        -2.31_rt, -1.766_rt, -0.396_rt, 1.678_rt, 3.931_rt, 6.082_rt,
        -1.6631_rt, -1.319_rt, -0.1438_rt, 1.7471_rt, 3.9409_rt, 6.0829_rt,
        -1.1064_rt, -0.8828_rt, 0.1094_rt, 1.8279_rt, 3.9534_rt, 6.0841_rt,
        -0.6344_rt, -0.496_rt, 0.3395_rt, 1.9168_rt, 3.9699_rt, 6.0862_rt,
        -0.2568_rt, -0.1555_rt, 0.5489_rt, 2.0163_rt, 3.9906_rt, 6.0893_rt,
        0.081_rt, 0.158_rt, 0.746_rt, 2.114_rt, 4.013_rt, 6.093_rt,
        0.3961_rt, 0.4448_rt, 0.9304_rt, 2.2026_rt, 4.0363_rt, 6.0976_rt,
        0.6673_rt, 0.6964_rt, 1.0985_rt, 2.2849_rt, 4.0614_rt, 6.1033_rt,
        0.9009_rt, 0.9175_rt, 1.2525_rt, 2.3619_rt, 4.0882_rt, 6.1099_rt,
        1.1032_rt, 1.113_rt, 1.3947_rt, 2.4345_rt, 4.1161_rt, 6.1171_rt};

    for (int i = 1; i <= 6; ++i) {
        for (int j = 1; j <= 14; ++j) {
            datn(1,i,j) = datn0[(i-1)+6*(j-1)];
            datn(2,i,j) = datn1[(i-1)+6*(j-1)];
        }
    }

    // Evaluate the cubic interp parameters for ni56 electron capture
    // which is used in the langanke subroutine.
    for (auto k = 2; k <= 4; ++k) {
        rfdm(k) = 1.0_rt / ((rv(k-1)-rv(k))*(rv(k-1)-rv(k+1))*(rv(k-1)-rv(k+2)));
        rfd0(k) = 1.0_rt / ((rv(k)-rv(k-1))*(rv(k)-rv(k+1))*(rv(k)-rv(k+2)));
        rfd1(k) = 1.0_rt / ((rv(k+1)-rv(k-1))*(rv(k+1)-rv(k))*(rv(k+1)-rv(k+2)));
        rfd2(k) = 1.0_rt / ((rv(k+2)-rv(k-1))*(rv(k+2)-rv(k))*(rv(k+2)-rv(k+1)));
    }

    for (auto j = 2; j <= 12; ++j) {
        tfdm(j) = 1.0_rt / ((tv(j-1)-tv(j))*(tv(j-1)-tv(j+1))*(tv(j-1)-tv(j+2)));
        tfd0(j) = 1.0_rt / ((tv(j)-tv(j-1))*(tv(j)-tv(j+1))*(tv(j)-tv(j+2)));
        tfd1(j) = 1.0_rt / ((tv(j+1)-tv(j-1))*(tv(j+1)-tv(j))*(tv(j+1)-tv(j+2)));
        tfd2(j) = 1.0_rt / ((tv(j+2)-tv(j-1))*(tv(j+2)-tv(j))*(tv(j+2)-tv(j+1)));
    }
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_c12ag(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    if (use_c12ag_deboer17) {
        // This version computes the nuclear reaction rate for 12C(a,g)16O and its inverse
        // using fit parameters from Deboer et al. 2017 (https://doi.org/10.1103/RevModPhys.89.035007).

        // from Table XXVI of deboer + 2017
        // non-resonant contributions to the reaction
        constexpr amrex::Real a0_nr = 24.1e0_rt;
        constexpr amrex::Real a1_nr = 0e0_rt;
        constexpr amrex::Real a2_nr = -32e0_rt;
        constexpr amrex::Real a3_nr = -5.9e0_rt;
        constexpr amrex::Real a4_nr = 1.8e0_rt;
        constexpr amrex::Real a5_nr = -0.17e0_rt;
        constexpr amrex::Real a6_nr = -2.0_rt/3.0_rt;

        const amrex::Real term_a0_nr = std::exp(a0_nr);
        const dual_t term_a1_nr = admath::exp(a1_nr*tf.t9i);
        const dual_t term_a2_nr = admath::exp(a2_nr*tf.t9i13);
        const dual_t term_a3_nr = admath::exp(a3_nr*tf.t913);
        const dual_t term_a4_nr = admath::exp(a4_nr*tf.t9);
        const dual_t term_a5_nr = admath::exp(a5_nr*tf.t953);
        const dual_t term_a6_nr = admath::pow(tf.t9,a6_nr);

        const dual_t term_nr = term_a0_nr * term_a1_nr * term_a2_nr *
                             term_a3_nr * term_a4_nr * term_a5_nr *
                             term_a6_nr;

        // resonant contributions to the reaction
        constexpr amrex::Real a0_r = 7.4e0_rt;
        constexpr amrex::Real a1_r = -30e0_rt;
        constexpr amrex::Real a2_r = 0e0_rt;
        constexpr amrex::Real a3_r = 0e0_rt;
        constexpr amrex::Real a4_r = 0e0_rt;
        constexpr amrex::Real a5_r = 0e0_rt;
        constexpr amrex::Real a6_r = -3.0e0_rt/2.0e0_rt;

        const amrex::Real term_a0_r = std::exp(a0_r);
        const dual_t term_a1_r = admath::exp(a1_r*tf.t9i);
        const dual_t term_a2_r = admath::exp(a2_r*tf.t9i13);
        const dual_t term_a3_r = admath::exp(a3_r*tf.t913);
        const dual_t term_a4_r = admath::exp(a4_r*tf.t9);
        const dual_t term_a5_r = admath::exp(a5_r*tf.t953);
        const dual_t term_a6_r = admath::pow(tf.t9,a6_r);

        const dual_t term_r = term_a0_r * term_a1_r * term_a2_r *
                            term_a3_r * term_a4_r * term_a5_r *
                            term_a6_r;

        // full rate is the sum of resonant and non-resonant contributions
        const dual_t term = term_nr + term_r;

        fr    = term * den;

        // first term is 9.8685e9_rt * T9**(2/3) * (M0*M1/M3)**(3/2)
        // see iliadis 2007 eqn. 3.44
        // ratio of partition functions are assumed to be unity
        const dual_t rev    = 5.1345573e10_rt * tf.t932 * admath::exp(-83.114082_rt*tf.t9i);

        rr     = rev * term;

    }
    else {

        constexpr amrex::Real q1 = 1.0e0_rt/12.222016e0_rt;

        // c12(a,g)o16
        const dual_t aa   = 1.0e0_rt + 0.0489e0_rt*tf.t9i23;

        const dual_t bb   = tf.t92*aa*aa;

        const dual_t cc   = admath::exp(-32.120e0_rt*tf.t9i13 - tf.t92*q1);

        const dual_t dd   = 1.0e0_rt + 0.2654e0_rt*tf.t9i23;

        const dual_t ee   = tf.t92*dd*dd;

        const dual_t ff   = admath::exp(-32.120e0_rt*tf.t9i13);

        const dual_t gg   = 1.25e3_rt * tf.t9i32 * admath::exp(-27.499_rt*tf.t9i);

        const dual_t hh   = 1.43e-2_rt * tf.t95 * admath::exp(-15.541_rt*tf.t9i);

        dual_t zz   = 1.0e0_rt/bb;
        const dual_t f1   = cc*zz;

        zz   = 1.0e0_rt/ee;
        const dual_t f2   = ff*zz;

        dual_t term = 1.04e8_rt*f1  + 1.76e8_rt*f2 + gg + hh;

        // 1.7 times cf88 value
        term     *= 1.7e0_rt;

        fr    = term * den;

        const dual_t rev = 5.13e10_rt * tf.t932 * admath::exp(-83.111_rt*tf.t9i);

        rr     = rev * term;

    }

}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_triplealf(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    constexpr amrex::Real rc28 = 0.1_rt;
    constexpr amrex::Real q1 = 1.0_rt / 0.009604e0_rt;
    constexpr amrex::Real q2 = 1.0e0_rt/0.055225e0_rt;

    // triple alfa to c12
    // this is a(a,g)be8
    const dual_t aa    = 7.40e5_rt * tf.t9i32 * admath::exp(-1.0663_rt*tf.t9i);
    const dual_t bb    = 4.164e9_rt * tf.t9i23 * admath::exp(-13.49_rt*tf.t9i13 - tf.t92*q1);
    const dual_t cc    = 1.0e0_rt + 0.031_rt*tf.t913 + 8.009_rt*tf.t923 + 1.732_rt*tf.t9
          + 49.883_rt*tf.t943 + 27.426_rt*tf.t953;
    const dual_t r2abe    = aa + bb * cc;

    // this is be8(a,g)c12
    const dual_t dd    = 130.0e0_rt * tf.t9i32 * admath::exp(-3.3364_rt*tf.t9i);
    const dual_t ee    = 2.510e7_rt * tf.t9i23 * admath::exp(-23.57_rt*tf.t9i13 - tf.t92*q2);
    const dual_t ff    = 1.0e0_rt + 0.018_rt*tf.t913 + 5.249_rt*tf.t923 + 0.650_rt*tf.t9 +
         19.176_rt*tf.t943 + 6.034_rt*tf.t953;
    const dual_t rbeac    = dd + ee * ff;

    // a factor
    const dual_t xx    = rc28 * 1.35e-07_rt * tf.t9i32 * admath::exp(-24.811_rt*tf.t9i);

    dual_t term;
    // high temperature rate
    if (tf.t9 > 0.08_rt) {
        term    = 2.90e-16_rt * r2abe * rbeac + xx;

    // low temperature rate
    } else {
        const dual_t uu   = 0.8e0_rt*admath::exp(-admath::pow(0.025_rt*tf.t9i,3.263_rt));
        const dual_t yy   = 0.2e0_rt + uu;
        //    ! fxt yy   = 0.01 + 0.2e0_rt + uu;
        const dual_t vv   = 4.0e0_rt*admath::exp(-admath::pow(tf.t9/0.025_rt, 9.227_rt));
        const dual_t zz   = 1.0e0_rt + vv;
        const dual_t zzinv   = 1.0e0_rt/zz;
        const dual_t f1   = 0.01e0_rt + yy * zzinv;
        //    ! fxt f1   = yy * 1/zz;
        term = 2.90e-16_rt * r2abe * rbeac * f1 +  xx;
    }

    // rates
    //      term    = 1.2e0_rt * term

    fr    = autodiff::eval(term * den) * den;

    const dual_t rev    = 2.00e20_rt*tf.t93*admath::exp(-84.424_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_c12c12(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // c12 + c12 reaction
    // this is the C12(C12,g)Mg24 rate from Caughlin & Fowler (1988)

    const dual_t aa      = 1.0e0_rt + 0.0396_rt*tf.t9;
    const dual_t zz      = 1.0e0_rt/aa;
    const dual_t t9a     = tf.t9*zz;
    const dual_t t9a13   = admath::cbrt(t9a);
    const dual_t t9a56   = admath::pow(t9a,5.0_rt/6.0_rt);
    const dual_t term    = 4.27e26_rt * t9a56 * tf.t9i32 *
         admath::exp(-84.165_rt/t9a13 - 2.12e-03_rt*tf.t93);

    // rates
    fr    = den * term;
    rr    = 0.0e0_rt;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_c12o16(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // c12 + o16 reaction; see cf88 references 47-4
    dual_t term;
    if (tf.t9 >= 0.5_rt) {
        dual_t aa     = 1.0e0_rt + 0.055_rt*tf.t9;
        dual_t zz     = 1.0e0_rt/aa;
        dual_t t9a    = tf.t9*zz;
        dual_t t9a13  = admath::cbrt(t9a);
        dual_t t9a23  = t9a13*t9a13;
        dual_t t9a56  = admath::pow(t9a,5.0_rt/6.0_rt);

        aa      = admath::exp(-0.18_rt*t9a*t9a);
        dual_t bb      = 1.06e-03_rt*admath::exp(2.562_rt*t9a23);
        dual_t cc      = aa + bb;
        zz      = 1.0e0_rt/cc;
        term    = 1.72e31_rt * t9a56 * tf.t9i32 * admath::exp(-106.594_rt/t9a13) * zz;

    } else {
        // term    = 2.6288035e-29_rt
        term    = 0.0e0_rt;
    }

    // the rates
    fr    = den * term;
    rr    = 0.0e0_rt;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_o16o16(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // o16 + o16
    // this is the O16(O16,g)S32 rate from Caughlin & Fowler (1988)

    const dual_t term  = 7.10e36_rt * tf.t9i23 *
         admath::exp(-135.93_rt * tf.t9i13 - 0.629_rt*tf.t923
         - 0.445_rt*tf.t943 + 0.0103_rt*tf.t92);

    // rates
    fr    = den * term;
    rr    = 0.0e0_rt;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_o16ag(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    constexpr amrex::Real q1 = 1.0e0_rt/2.515396e0_rt;

    // o16(a,g)ne20
    const dual_t term1   = 9.37e9_rt * tf.t9i23 * admath::exp(-39.757_rt*tf.t9i13 - tf.t92*q1);

    const dual_t aa      = 62.1_rt * tf.t9i32 * admath::exp(-10.297_rt*tf.t9i);
    const dual_t bb      = 538.0e0_rt * tf.t9i32 * admath::exp(-12.226_rt*tf.t9i);
    const dual_t cc      = 13.0e0_rt * tf.t92 * admath::exp(-20.093_rt*tf.t9i);
    const dual_t term2   = aa + bb + cc;

    const dual_t term    = term1 + term2;

    // the rates
    fr    = den * term;

    const dual_t rev      = 5.65e10_rt*tf.t932*admath::exp(-54.937_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_ne20ag(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    constexpr amrex::Real rc102 = 0.1_rt;
    constexpr amrex::Real q1 = 1.0e0_rt/4.923961e0_rt;

    // ne20(a,g)mg24
    dual_t aa   = 4.11e11_rt * tf.t9i23 * admath::exp(-46.766_rt*tf.t9i13 - tf.t92*q1);
    dual_t bb   = 1.0e0_rt + 0.009_rt*tf.t913 + 0.882_rt*tf.t923 + 0.055_rt*tf.t9
         + 0.749_rt*tf.t943 + 0.119_rt*tf.t953;
    const dual_t term1  = aa * bb;

    aa   = 5.27e3_rt * tf.t9i32 * admath::exp(-15.869_rt*tf.t9i);
    bb   = 6.51e3_rt * tf.t912 * admath::exp(-16.223_rt*tf.t9i);
    const dual_t term2  = aa + bb;

    aa   = 42.1_rt * tf.t9i32 * admath::exp(-9.115_rt*tf.t9i);
    bb   =  32.0_rt * tf.t9i23 * admath::exp(-9.383_rt*tf.t9i);
    const dual_t term3  = rc102 * (aa + bb);

    aa  = 5.0e0_rt*admath::exp(-18.960_rt*tf.t9i);
    bb  = 1.0e0_rt + aa;
    const dual_t zz      = 1.0e0_rt/bb;
    const dual_t term    = (term1 + term2 + term3)*zz;

    // the rates
    fr    = den * term;

    const dual_t rev      = 6.01e10_rt * tf.t932 * admath::exp(-108.059_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_mg24ag(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    constexpr amrex::Real rc121 = 0.1e0_rt;

    // 24mg(a,g)28si
    const dual_t aa    = 4.78e1_rt * tf.t9i32 * admath::exp(-13.506_rt*tf.t9i);
    const dual_t bb    =  2.38e3_rt * tf.t9i32 * admath::exp(-15.218_rt*tf.t9i);
    const dual_t cc    = 2.47e2_rt * tf.t932 * admath::exp(-15.147_rt*tf.t9i);
    const dual_t dd    = rc121 * 1.72e-09_rt * tf.t9i32 * admath::exp(-5.028_rt*tf.t9i);
    const dual_t ee    = rc121* 1.25e-03_rt * tf.t9i32 * admath::exp(-7.929_rt*tf.t9i);
    const dual_t ff    = rc121 * 2.43e1_rt * tf.t9i * admath::exp(-11.523_rt*tf.t9i);

    const dual_t gg    = 5.0e0_rt*admath::exp(-15.882_rt*tf.t9i);
    const dual_t hh    = 1.0e0_rt + gg;
    const dual_t hhi   = 1.0e0_rt/hh;
    const dual_t term    = (aa + bb + cc + dd + ee + ff) * hhi;

    // the rates
    fr    = den * term;

    const dual_t rev      = 6.27e10_rt * tf.t932 * admath::exp(-115.862_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_mg24ap(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    constexpr amrex::Real rc148 = 0.1_rt;
    constexpr amrex::Real q1 = 1.0_rt / 0.024649e0_rt;

    // 24mg(a,p)al27
    dual_t aa     = 1.10e8_rt * tf.t9i23 * admath::exp(-23.261_rt*tf.t9i13 - tf.t92*q1);
    dual_t bb     =  1.0e0_rt + 0.018_rt*tf.t913 + 12.85_rt*tf.t923 + 1.61_rt*tf.t9
         + 89.87_rt*tf.t943 + 28.66_rt*tf.t953;
    const dual_t term1  = aa * bb;

    aa     = 129.0e0_rt * tf.t9i32 * admath::exp(-2.517_rt*tf.t9i);
    bb     = 5660.0e0_rt * tf.t972 * admath::exp(-3.421_rt*tf.t9i);
    const dual_t cc     = rc148 * 3.89e-08_rt * tf.t9i32 * admath::exp(-0.853_rt*tf.t9i);
    const dual_t dd     = rc148 * 8.18e-09_rt * tf.t9i32 * admath::exp(-1.001_rt*tf.t9i);
    const dual_t term2  = aa + bb + cc + dd;

    const dual_t ee     = 1.0_rt/3.0_rt*admath::exp(-9.792_rt*tf.t9i);
    const dual_t ff     =  2.0_rt/3.0_rt * admath::exp(-11.773_rt*tf.t9i);
    const dual_t gg     = 1.0e0_rt + ee + ff;

    const dual_t term    = (term1 + term2)/gg;

    // the rates
    const dual_t rev      = 1.81_rt * admath::exp(-18.572_rt*tf.t9i);

    fr    = den * rev * term;
    rr    = den * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_al27pg(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // al27(p,g)si28
    // champagne 1996

    const dual_t aa  = 1.32e9_rt * tf.t9i23 * admath::exp(-23.26_rt*tf.t9i13);
    const dual_t bb  = autodiff::eval(3.22e-10_rt * tf.t9i32 * admath::exp(-0.836_rt*tf.t9i))*0.17_rt;
    const dual_t cc  = 1.74e0_rt * tf.t9i32 * admath::exp(-2.269_rt*tf.t9i);
    const dual_t dd  = 9.92e0_rt * tf.t9i32 * admath::exp(-2.492_rt*tf.t9i);
    const dual_t ee  = 4.29e1_rt * tf.t9i32 * admath::exp(-3.273_rt*tf.t9i);
    const dual_t ff  = 1.34e2_rt * tf.t9i32 * admath::exp(-3.654_rt*tf.t9i);
    const dual_t gg  = 1.77e4_rt * admath::pow(tf.t9,0.53_rt) * admath::exp(-4.588_rt*tf.t9i);
    const dual_t term    = aa + bb + cc + dd + ee + ff + gg;

    // rates
    fr    = den * term;

    const dual_t rev      = 1.13e11_rt * tf.t932 * admath::exp(-134.434_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_al27pg_old(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    constexpr amrex::Real rc147 = 0.1_rt;
    constexpr amrex::Real q1 = 1.0e0_rt/0.024025e0_rt;

    // 27al(p,g)si28  cf88
    const dual_t aa  = 1.67e8_rt * tf.t9i23 * admath::exp(-23.261_rt*tf.t9i13 - tf.t92*q1);
    const dual_t bb  = 1.0e0_rt + 0.018_rt*tf.t913 + 5.81_rt*tf.t923 + 0.728_rt*tf.t9
         + 27.31_rt*tf.t943 + 8.71_rt*tf.t953;
    const dual_t cc  = aa*bb;

    const dual_t dd  = 2.20e0_rt * tf.t9i32 * admath::exp(-2.269_rt*tf.t9i);
    const dual_t ee  = 1.22e1_rt * tf.t9i32 * admath::exp(-2.491_rt*tf.t9i);
    const dual_t ff  =  1.50e4_rt * tf.t9 * admath::exp(-4.112_rt*tf.t9i);
    const dual_t gg  = rc147 * 6.50e-10_rt * tf.t9i32 * admath::exp(-0.853_rt*tf.t9i);
    const dual_t hh  = rc147 * 1.63e-10_rt * tf.t9i32 * admath::exp(-1.001_rt*tf.t9i);

    const dual_t xx     = 1.0_rt/3.0_rt*admath::exp(-9.792_rt*tf.t9i);
    const dual_t yy     =  2.0_rt/3.0_rt * admath::exp(-11.773_rt*tf.t9i);
    const dual_t zz     = 1.0e0_rt + xx + yy;
    const dual_t pp      = 1.0e0_rt/zz;
    const dual_t term    = (cc + dd + ee + ff + gg + hh)*pp;

    // the rates
    fr    = den * term;

    const dual_t rev      = 1.13e11_rt*tf.t932*admath::exp(-134.434_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_si28ag(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // si28(a,g)s32
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 6.340e-2_rt*z + 2.541e-3_rt*z2 - 2.900e-4_rt*z3;

    const dual_t term    = 4.82e22_rt * tf.t9i23 * admath::exp(-61.015_rt * tf.t9i13 * aa);

    // the rates
    fr    = den * term;

    const dual_t rev      = 6.461e10_rt * tf.t932 * admath::exp(-80.643_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_si28ap(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // si28(a,p)p31
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 2.798e-3_rt*z + 2.763e-3_rt*z2 - 2.341e-4_rt*z3;

    const dual_t term    = 4.16e13_rt * tf.t9i23 * admath::exp(-25.631_rt * tf.t9i13 * aa);

    // the rates
    const dual_t rev      = 0.5825e0_rt * admath::exp(-22.224_rt*tf.t9i);

    fr    = den * rev * term;
    rr    = den * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_p31pg(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // p31(p,g)s32
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 1.928e-1_rt*z - 1.540e-2_rt*z2 + 6.444e-4_rt*z3;

    const dual_t term    = 1.08e16_rt * tf.t9i23 * admath::exp(-27.042_rt * tf.t9i13 * aa);

    // the rates
    fr    = den * term;

    const dual_t rev      = 3.764e10_rt * tf.t932 * admath::exp(-102.865_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_s32ag(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // s32(a,g)ar36
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 4.913e-2_rt*z + 4.637e-3_rt*z2 - 4.067e-4_rt*z3;

    const dual_t term    = 1.16e24_rt * tf.t9i23 * admath::exp(-66.690_rt * tf.t9i13 * aa);

    // the rates
    fr    = den * term;

    const dual_t rev      = 6.616e10_rt * tf.t932 * admath::exp(-77.080_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_s32ap(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // s32(a,p)cl35
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 1.041e-1_rt*z - 1.368e-2_rt*z2 + 6.969e-4_rt*z3;

    const dual_t term    = 1.27e16_rt * tf.t9i23 * admath::exp(-31.044_rt * tf.t9i13 * aa);

    // the rates
    const dual_t rev      = 1.144_rt * admath::exp(-21.643_rt*tf.t9i);

    fr    = den * rev * term;
    rr    = den * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_cl35pg(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // cl35(p,g)ar36
    const dual_t aa    = 1.0e0_rt + 1.761e-1_rt*tf.t9 - 1.322e-2_rt*tf.t92 + 5.245e-4_rt*tf.t93;
    const dual_t term    =  4.48e16_rt * tf.t9i23 * admath::exp(-29.483_rt * tf.t9i13 * aa);

    // the rates
    fr    = den * term;

    const dual_t rev      = 7.568e10_rt*tf.t932*admath::exp(-98.722_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_ar36ag(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // ar36(a,g)ca40
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 1.458e-1_rt*z - 1.069e-2_rt*z2 + 3.790e-4_rt*z3;

    const dual_t term    = 2.81e30_rt * tf.t9i23 * admath::exp(-78.271_rt * tf.t9i13 * aa);

    // the rates
    fr    = den * term;

    const dual_t rev      = 6.740e10_rt * tf.t932 * admath::exp(-81.711_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_ar36ap(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // ar36(a,p)k39
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 4.826e-3_rt*z - 5.534e-3_rt*z2 + 4.021e-4_rt*z3;

    const dual_t term    = 2.76e13_rt * tf.t9i23 * admath::exp(-34.922_rt * tf.t9i13 * aa);

    // the rates
    const dual_t rev      = 1.128_rt*admath::exp(-14.959_rt*tf.t9i);

    fr    = den * rev * term;

    rr    = den * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_k39pg(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // k39(p,g)ca40
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 1.622e-1_rt*z - 1.119e-2_rt*z2 + 3.910e-4_rt*z3;

    const dual_t term    = 4.09e16_rt * tf.t9i23 * admath::exp(-31.727_rt * tf.t9i13 * aa);

    // the rates
    fr    = den * term;

    const dual_t rev      = 7.600e10_rt * tf.t932 * admath::exp(-96.657_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_ca40ag(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // ca40(a,g)ti44
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 1.650e-2_rt*z + 5.973e-3_rt*z2 - 3.889e-04_rt*z3;

    const dual_t term    = 4.66e24_rt * tf.t9i23 * admath::exp(-76.435_rt * tf.t9i13 * aa);

    // the rates
    fr    = den * term;

    const dual_t rev      = 6.843e10_rt * tf.t932 * admath::exp(-59.510_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_ca40ap(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // ca40(a,p)sc43
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt - 1.206e-2_rt*z + 7.753e-3_rt*z2 - 5.071e-4_rt*z3;

    const dual_t term    = 4.54e14_rt * tf.t9i23 * admath::exp(-32.177_rt * tf.t9i13 * aa);

    // the rates
    const dual_t rev      = 2.229_rt * admath::exp(-40.966_rt*tf.t9i);

    fr    = den * rev * term;
    rr    = den * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_sc43pg(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // sc43(p,g)ca40
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 1.023e-1_rt*z - 2.242e-3_rt*z2 - 5.463e-5_rt*z3;

    const dual_t term    = 3.85e16_rt * tf.t9i23 * admath::exp(-33.234_rt * tf.t9i13 * aa);

    // the rates
    fr    = den * term;

    const dual_t rev      = 1.525e11_rt * tf.t932 * admath::exp(-100.475_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_ti44ag(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // ti44(a,g)cr48
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 1.066e-1_rt*z - 1.102e-2_rt*z2 + 5.324e-4_rt*z3;

    const dual_t term    = 1.37e26_rt * tf.t9i23 * admath::exp(-81.227_rt * tf.t9i13 * aa);

    // the rates
    fr    = den * term;

    const dual_t rev      = 6.928e10_rt*tf.t932*admath::exp(-89.289_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_ti44ap(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // ti44(a,p)v47
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 2.655e-2_rt*z - 3.947e-3_rt*z2 + 2.522e-4_rt*z3;

    const dual_t term    = 6.54e20_rt * tf.t9i23 * admath::exp(-66.678_rt * tf.t9i13 * aa);

    // the rates
    const dual_t rev      = 1.104_rt * admath::exp(-4.723_rt*tf.t9i);

    fr    = den * rev * term;
    rr    = den * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_v47pg(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // v47(p,g)cr48
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 9.979e-2_rt*z - 2.269e-3_rt*z2 - 6.662e-5_rt*z3;

    const dual_t term    = 2.05e17_rt * tf.t9i23 * admath::exp(-35.568_rt * tf.t9i13 * aa);

    // the rates
    fr    = den * term;

    const dual_t rev      = 7.649e10_rt*tf.t932*admath::exp(-93.999_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_cr48ag(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // cr48(a,g)fe52
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 6.325e-2_rt*z - 5.671e-3_rt*z2 + 2.848e-4_rt*z3;

    const dual_t term    = 1.04e23_rt * tf.t9i23 * admath::exp(-81.420_rt * tf.t9i13 * aa);

    // the rates
    fr    = den * term;

    const dual_t rev      = 7.001e10_rt * tf.t932 * admath::exp(-92.177_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_cr48ap(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // cr48(a,p)mn51
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 1.384e-2_rt*z + 1.081e-3_rt*z2 - 5.933e-5_rt*z3;

    const dual_t term    = 1.83e26_rt * tf.t9i23 * admath::exp(-86.741_rt * tf.t9i13 * aa);

    // the rates
    fr    = den * term;

    const dual_t rev      = 0.6087_rt*admath::exp(-6.510_rt*tf.t9i);

    rr    = den * rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_mn51pg(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // mn51(p,g)fe52
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 8.922e-2_rt*z - 1.256e-3_rt*z2 - 9.453e-5_rt*z3;

    const dual_t term    = 3.77e17_rt * tf.t9i23 * admath::exp(-37.516_rt * tf.t9i13 * aa);

    // the rates
    fr    = den * term;

    const dual_t rev      = 1.150e11_rt*tf.t932*admath::exp(-85.667_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_fe52ag(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // fe52(a,g)ni56
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 7.846e-2_rt*z - 7.430e-3_rt*z2 + 3.723e-4_rt*z3;

    const dual_t term    = 1.05e27_rt * tf.t9i23 * admath::exp(-91.674_rt * tf.t9i13 * aa);

    // the rates
    fr    = den * term;

    const dual_t rev      = 7.064e10_rt*tf.t932*admath::exp(-92.850_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_fe52ap(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // fe52(a,p)co55
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 1.367e-2_rt*z + 7.428e-4_rt*z2 - 3.050e-5_rt*z3;

    const dual_t term    = 1.30e27_rt * tf.t9i23 * admath::exp(-91.674_rt * tf.t9i13 * aa);

    // the rates
    fr    = den * term;

    const dual_t rev      = 0.4597_rt*admath::exp(-9.470_rt*tf.t9i);

    rr    = den * rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_co55pg(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // co55(p,g)ni56
    const dual_t z     = admath::min(10.0e0_rt, tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 9.894e-2_rt*z - 3.131e-3_rt*z2 - 2.160e-5_rt*z3;

    const dual_t term    = 1.21e18_rt * tf.t9i23 * admath::exp(-39.604_rt * tf.t9i13 * aa);

    // the rates
    fr    = den * term;

    const dual_t rev      = 1.537e11_rt*tf.t932*admath::exp(-83.382_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_pp(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // p(p,e+nu)d
    dual_t term;
    if (tf.t9 <= 3.0_rt) {
       const dual_t aa   = 4.01e-15_rt * tf.t9i23 * admath::exp(-3.380e0_rt*tf.t9i13);
       const dual_t bb   = 1.0e0_rt + 0.123e0_rt*tf.t913 + 1.09e0_rt*tf.t923 + 0.938e0_rt*tf.t9;
       term    = aa * bb;
    } else {
       term    = 1.1581136e-15_rt;
    }

    // rate
    fr    = den * term;
    rr    = 0.0e0_rt;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_png(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // p(n,g)d
    // smith,kawano,malany 1992

    const dual_t aa      = 1.0e0_rt - 0.8504_rt*tf.t912 + 0.4895_rt*tf.t9
         - 0.09623_rt*tf.t932 + 8.471e-3_rt*tf.t92
         - 2.80e-4_rt*tf.t952;
    const dual_t term    = 4.742e4_rt * aa;

    // wagoner,schramm 1977
    //      aa      = 1.0e0_rt - 0.86*tf.t912 + 0.429*tf.t9
    //      term    = 4.4e4_rt * aa

    // rates
    fr    = den * term;

    const dual_t rev      = 4.71e9_rt * tf.t932 * admath::exp(-25.82_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_dpg(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // d(p,g)he3
    const dual_t aa      = 2.24e3_rt * tf.t9i23 * admath::exp(-3.720_rt*tf.t9i13);
    const dual_t bb      = 1.0e0_rt + 0.112_rt*tf.t913 + 3.38_rt*tf.t923 + 2.65_rt*tf.t9;
    const dual_t term    = aa * bb;

    // rates
    fr    = den * term;

    const dual_t rev      = 1.63e10_rt * tf.t932 * admath::exp(-63.750_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_he3ng(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // he3(n,g)he4
    const dual_t term    = 6.62_rt * (1.0e0_rt + 905.0_rt*tf.t9);

    // rates
    fr    = den * term;

    const dual_t rev      = 2.61e10_rt * tf.t932 * admath::exp(-238.81_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_he3he3(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // he3(he3,2p)he4
    const dual_t aa   = 6.04e10_rt * tf.t9i23 * admath::exp(-12.276_rt*tf.t9i13);
    const dual_t bb   = 1.0e0_rt + 0.034_rt*tf.t913 - 0.522_rt*tf.t923 - 0.124_rt*tf.t9
         + 0.353_rt*tf.t943 + 0.213_rt*tf.t953;
    const dual_t term    = aa * bb;

    // rates
    fr    = den * term;

    const dual_t rev      = 3.39e-10_rt * tf.t9i32 * admath::exp(-149.230_rt*tf.t9i);

    rr    = den * den * rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_he3he4(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // he3(he4,g)be7
    const dual_t aa      = 1.0e0_rt + 0.0495_rt*tf.t9;

    dual_t zz      = 1.0e0_rt/aa;
    const dual_t t9a     = tf.t9*zz;
    const dual_t t9a13   = admath::cbrt(t9a);
    const dual_t t9a56   = admath::pow(t9a,5.0_rt/6.0_rt);

    const dual_t term    = 5.61e6_rt * t9a56 * tf.t9i32 * admath::exp(-12.826_rt/t9a13);

    // rates
    fr    = den * term;

    const dual_t rev      = 1.11e10_rt * tf.t932 * admath::exp(-18.423_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_c12pg(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    constexpr amrex::Real q1 = 1.0e0_rt/2.25e0_rt;

    // c12(p,g)13n
    const dual_t aa   = 2.04e7_rt * tf.t9i23 * admath::exp(-13.69_rt*tf.t9i13 - tf.t92*q1);
    const dual_t bb   = 1.0e0_rt + 0.03_rt*tf.t913 + 1.19_rt*tf.t923 + 0.254_rt*tf.t9
         + 2.06_rt*tf.t943 + 1.12_rt*tf.t953;
    const dual_t cc   = aa * bb;

    const dual_t dd   = 1.08e5_rt * tf.t9i32 * admath::exp(-4.925_rt*tf.t9i);

    const dual_t ee   = 2.15e5_rt * tf.t9i32 * admath::exp(-18.179_rt*tf.t9i);

    const dual_t term    = cc + dd + ee;

    // rates
    fr    = den * term;

    const dual_t rev      = 8.84e9_rt * tf.t932 * admath::exp(-22.553_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_n14pg(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    constexpr amrex::Real q1 = 1.0e0_rt/10.850436e0_rt;

    // n14(p,g)o15
    const dual_t aa  = 4.90e7_rt * tf.t9i23 * admath::exp(-15.228_rt*tf.t9i13 - tf.t92*q1);
    const dual_t bb   = 1.0e0_rt + 0.027_rt*tf.t913 - 0.778_rt*tf.t923 - 0.149_rt*tf.t9
         + 0.261_rt*tf.t943 + 0.127_rt*tf.t953;
    const dual_t cc   = aa * bb;

    const dual_t dd   = 2.37e3_rt * tf.t9i32 * admath::exp(-3.011_rt*tf.t9i);

    const dual_t ee   = 2.19e4_rt * admath::exp(-12.530_rt*tf.t9i);

    const dual_t term    = cc + dd + ee;

    // rates
    fr    = den * term;

    const dual_t rev    = 2.70e10_rt * tf.t932 * admath::exp(-84.678_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_n15pg(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    constexpr amrex::Real q1 = 1.0_rt / 0.2025_rt;

    // n15(p,g)o16
    const dual_t aa  = 9.78e8_rt * tf.t9i23 * admath::exp(-15.251_rt*tf.t9i13 - tf.t92*q1);
    const dual_t bb   = 1.0e0_rt  + 0.027_rt*tf.t913 + 0.219_rt*tf.t923 + 0.042_rt*tf.t9
         + 6.83_rt*tf.t943 + 3.32_rt*tf.t953;
    const dual_t cc   = aa * bb;

    const dual_t dd   = 1.11e4_rt*tf.t9i32*admath::exp(-3.328_rt*tf.t9i);

    const dual_t ee   = 1.49e4_rt*tf.t9i32*admath::exp(-4.665_rt*tf.t9i);

    const dual_t ff   = 3.8e6_rt*tf.t9i32*admath::exp(-11.048_rt*tf.t9i);

    const dual_t term    = cc + dd + ee + ff;

    // rates
    fr    = den * term;

    const dual_t rev      = 3.62e10_rt * tf.t932 * admath::exp(-140.734_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_n15pa(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    constexpr amrex::Real theta = 0.1_rt;
    constexpr amrex::Real q1 = 1.0_rt / 0.272484_rt;

    // n15(p,a)c12
    const dual_t aa  = 1.08e12_rt*tf.t9i23*admath::exp(-15.251_rt*tf.t9i13 - tf.t92*q1);
    const dual_t bb   = 1.0e0_rt + 0.027_rt*tf.t913 + 2.62_rt*tf.t923 + 0.501_rt*tf.t9
         + 5.36_rt*tf.t943 + 2.60_rt*tf.t953;
    const dual_t cc   = aa * bb;

    const dual_t dd   = 1.19e8_rt * tf.t9i32 * admath::exp(-3.676_rt*tf.t9i);

    const dual_t ee   = 5.41e8_rt * tf.t9i12 * admath::exp(-8.926_rt*tf.t9i);

    const dual_t ff   = theta * 4.72e8_rt * tf.t9i32 * admath::exp(-7.721_rt*tf.t9i);

    const dual_t gg   = theta * 2.20e9_rt * tf.t9i32 * admath::exp(-11.418_rt*tf.t9i);

    const dual_t term    = cc + dd + ee + ff + gg;

    // rates
    fr    = den * term;

    const dual_t rev      = 7.06e-01_rt*admath::exp(-57.625_rt*tf.t9i);

    rr    = den * rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_o16pg(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // o16(p,g)f17
    const dual_t aa  = admath::exp(-0.728_rt*tf.t923);
    const dual_t bb  = 1.0e0_rt + 2.13_rt * (1.0e0_rt - aa);
    const dual_t cc  = tf.t923 * bb;
    const dual_t dd   = admath::exp(-16.692_rt*tf.t9i13);
    const dual_t zz   = 1.0e0_rt/cc;
    const dual_t ee   = dd*zz;
    const dual_t term    = 1.50e8_rt * ee;

    // rates
    fr    = den * term;

    const dual_t rev      = 3.03e9_rt*tf.t932*admath::exp(-6.968_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_n14ag(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    constexpr amrex::Real q1 = 1.0e0_rt/0.776161e0_rt;

    // n14(a,g)f18
    const dual_t aa  = 7.78e9_rt * tf.t9i23 * admath::exp(-36.031_rt*tf.t9i13- tf.t92*q1);
    const dual_t bb   = 1.0e0_rt + 0.012_rt*tf.t913 + 1.45_rt*tf.t923 + 0.117_rt*tf.t9
         + 1.97_rt*tf.t943 + 0.406_rt*tf.t953;
    const dual_t cc   = aa * bb;

    const dual_t dd   = 2.36e-10_rt * tf.t9i32 * admath::exp(-2.798_rt*tf.t9i);

    const dual_t ee   = 2.03_rt * tf.t9i32 * admath::exp(-5.054_rt*tf.t9i);

    const dual_t ff   = 1.15e4_rt * tf.t9i23 * admath::exp(-12.310_rt*tf.t9i);

    const dual_t term    = cc + dd + ee + ff;

    // rates
    fr    = den * term;

    const dual_t rev      = 5.42e10_rt * tf.t932 * admath::exp(-51.236_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_fe52ng(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // fe52(n,g)fe53
    const dual_t tq2     = tf.t9 - 0.348e0_rt;
    const dual_t term    = 9.604e5_rt * admath::exp(-0.0626_rt*tq2);

    // the rates
    fr    = den * term;

    const dual_t rev      = 2.43e9_rt * tf.t932 * admath::exp(-123.951_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_fe53ng(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // fe53(n,g)fe54
    const dual_t tq1   = tf.t9/0.348_rt;
    const dual_t tq10  = admath::pow(tq1, 0.10_rt);
    const dual_t tq2   = tf.t9 - 0.348e0_rt;

    const dual_t term    = 1.817e6_rt * tq10 * admath::exp(-0.06319_rt*tq2);

    // the rates
    fr    = den * term;

    const dual_t rev      = 1.56e11_rt * tf.t932 * admath::exp(-155.284_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_fe54ng(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // fe54(n,g)fe55
    const dual_t aa   =  2.307390e1_rt - 7.931795e-02_rt * tf.t9i + 7.535681e0_rt * tf.t9i13
         - 1.595025e1_rt * tf.t913 + 1.377715e0_rt * tf.t9 - 1.291479e-01_rt * tf.t953
         + 6.707473e0_rt * admath::log(tf.t9);

    dual_t term;
    if (aa < 200.0_rt) {
       term    = admath::exp(aa);
    } else {
       term    = std::exp(200.0e0_rt);
    }

    const dual_t bb  = 4.800293e9_rt * tf.t932 * admath::exp(-1.078986e2_rt * tf.t9i);

    // reverse rate
    rr    = term*bb;

    // forward rate
    fr    = term*den;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_fe54pg(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // fe54(p,g)co55
    const dual_t z     = admath::min(10.0e0_rt,tf.t9);
    const dual_t z2    = z*z;
    const dual_t z3    = z2*z;
    const dual_t aa    = 1.0e0_rt + 9.593e-2_rt*z - 3.445e-3_rt*z2 + 8.594e-5_rt*z3;

    const dual_t term    = 4.51e17_rt * tf.t9i23 * admath::exp(-38.483_rt * tf.t9i13 * aa);

    // the rates
    fr    = den * term;

    const dual_t rev      = 2.400e9_rt * tf.t932 * admath::exp(-58.605_rt*tf.t9i);

    rr    = rev * term;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_fe54ap(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // fe54(a,p)co57
    const dual_t aa   =  3.97474900e1_rt - 6.06543100e0_rt * tf.t9i + 1.63239600e2_rt * tf.t9i13
         - 2.20457700e2_rt * tf.t913 + 8.63980400e0_rt * tf.t9 - 3.45841300e-01_rt * tf.t953
         + 1.31464200e2_rt * admath::log(tf.t9);

    dual_t term;
    if (aa < 200.0_rt) {
       term    = admath::exp(aa);
    } else {
       term    = std::exp(200.0e0_rt);
    }

    const dual_t bb  = 2.16896000e0_rt  * admath::exp(-2.05631700e1_rt * tf.t9i);

    // reverse rate
    rr    = term*den;

    // forward rate
    fr    = rr*bb;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_fe55ng(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // fe55(n,g)fe56
    const dual_t aa   =  1.954115e1_rt - 6.834029e-02_rt * tf.t9i + 5.379859e0_rt * tf.t9i13
         - 8.758150e0_rt * tf.t913 + 5.285107e-01_rt * tf.t9 - 4.973739e-02_rt  * tf.t953
         + 4.065564e0_rt  * admath::log(tf.t9);

    dual_t term;
    if (aa < 200.0_rt) {
       term    = admath::exp(aa);
    } else {
       term    = std::exp(200.0e0_rt);
    }

    const dual_t bb  = 7.684279e10_rt  * tf.t932 * admath::exp(-1.299472e2_rt  * tf.t9i);

    // reverse rate
    rr    = term*bb;

    // forward rate
    fr    = term*den;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FLATTEN inline
void rate_fe56pg(const tf_t<dual_t>& tf, const amrex::Real den, dual_t& fr, dual_t& rr)
{
    // fe56(p,g)co57

    const auto aa   =  1.755960e2_rt - 7.018872e0_rt * tf.t9i + 2.800131e2_rt * tf.t9i13
         - 4.749343e2_rt * tf.t913 + 2.683860e1_rt * tf.t9 - 1.542324e0_rt  * tf.t953
         + 2.315911e2_rt  * admath::log(tf.t9);

    dual_t term;
    if (aa < 200.0_rt) {
       term    = admath::exp(aa);
    } else {
       term    = std::exp(200.0e0_rt);
    }

    const dual_t bb  = 2.402486e9_rt * tf.t932 * admath::exp(-6.995192e1_rt * tf.t9i);

    // reverse rate
    rr    = term*bb;

    // forward rate
    fr    = term*den;
}


// this routine evaluates Langanke et al. 2000 fits for the ni56 electron
// capture rate rn56ec and neutrino loss rate sn56ec

// input:
// y56 = nickel56 molar abundance
// ye  = electron to baryon number, zbar/abar

// output:
// rn56ec = ni56 electron capture rate
// sn56ec = ni56 neutrino loss rate
AMREX_GPU_HOST_DEVICE inline
void langanke(const amrex::Real btemp, const amrex::Real bden,
              const amrex::Real y56, const amrex::Real ye,
              amrex::Real& rn56ec, amrex::Real& sn56ec)
{
    // calculate ni56 electron capture and neutrino loss rates
    rn56ec = 0.0_rt;
    sn56ec = 0.0_rt;

    if ( (btemp < 1.0e9_rt) || (bden*ye < 1.0e6_rt)) return;

    const amrex::Real t9    = amrex::min(btemp, 1.4e10_rt) * 1.0e-9_rt;
    const amrex::Real r     = amrex::max(6.0e0_rt, amrex::min(11.0e0_rt, std::log10(bden*ye)));
    const int jp    = amrex::min(amrex::max(2, int(t9)), 12);
    const int kp    = amrex::min(amrex::max(2, int(r)-5), 4);
    const amrex::Real rfm   = r - rv(kp-1);
    const amrex::Real rf0   = r - rv(kp);
    const amrex::Real rf1   = r - rv(kp+1);
    const amrex::Real rf2   = r - rv(kp+2);
    const amrex::Real dfacm = rf0*rf1*rf2*rfdm(kp);
    const amrex::Real dfac0 = rfm*rf1*rf2*rfd0(kp);
    const amrex::Real dfac1 = rfm*rf0*rf2*rfd1(kp);
    const amrex::Real dfac2 = rfm*rf0*rf1*rfd2(kp);
    const amrex::Real tfm   = t9 - tv(jp-1);
    const amrex::Real tf0   = t9 - tv(jp);
    const amrex::Real tf1   = t9 - tv(jp+1);
    const amrex::Real tf2   = t9 - tv(jp+2);
    const amrex::Real tfacm = tf0*tf1*tf2*tfdm(jp);
    const amrex::Real tfac0 = tfm*tf1*tf2*tfd0(jp);
    const amrex::Real tfac1 = tfm*tf0*tf2*tfd1(jp);
    const amrex::Real tfac2 = tfm*tf0*tf1*tfd2(jp);

    amrex::Array1D<amrex::Real, 1, 2> rnt;
    amrex::Array2D<amrex::Real, 1, 2, 1, 14> rne;

    // evaluate the spline fits
    for (auto jr = 1; jr <= 2; ++jr) {
        for (auto jd = jp-1; jd <= jp+2; ++jd) {
            rne(jr,jd) = dfacm*datn(jr,kp-1,jd) + dfac0*datn(jr,kp,jd)
                + dfac1*datn(jr,kp+1,jd) + dfac2*datn(jr,kp+2,jd);
        }
        rnt(jr) =  tfacm*rne(jr,jp-1) + tfac0*rne(jr,jp)
                 + tfac1*rne(jr,jp+1) + tfac2*rne(jr,jp+2);
    }

    // set the output
    rn56ec = std::pow(10.0e0_rt, rnt(1));
    sn56ec = 6.022548e23_rt * 1.60218e-6_rt * y56 * std::pow(10.0e0_rt, rnt(2));
}

// given the electron degeneracy parameter etakep (chemical potential
// without the electron's rest mass divided by kt) and the temperature
// temp, this routine calculates rates for
// electron capture on protons rpen (captures/sec/proton),
// positron capture on neutrons rnep (captures/sec/neutron),
// and their associated neutrino energy loss rates
// spenc (erg/sec/proton) and snepc (erg/sec/neutron)
AMREX_GPU_HOST_DEVICE inline
void ecapnuc(const amrex::Real etakep, const amrex::Real temp,
             amrex::Real& rpen, amrex::Real& rnep,
             amrex::Real& spenc, amrex::Real& snepc)
{
    constexpr amrex::Real qn1    = -2.0716446e-06_rt;
    constexpr amrex::Real ftinv  = 1.0e0_rt/1083.9269e0_rt;
    constexpr amrex::Real twoln  = 0.6931472e0_rt;
    constexpr amrex::Real cmk5   = 1.3635675e-49_rt;
    constexpr amrex::Real cmk6   = 2.2993864e-59_rt;
    constexpr amrex::Real bk     = 1.38062e-16_rt;
    constexpr amrex::Real qn2    = 2.0716446e-06_rt;
    constexpr amrex::Real c2me   = 8.1872665e-07_rt;
    constexpr amrex::Real pi     = 3.1415927e0_rt;
    constexpr amrex::Real pi2    = pi * pi;

    // c2me is the constant used to convert the neutrino energy
    // loss rate from mec2/s (as in the paper) to ergs/particle/sec.

    // initialize
    rpen   = 0.0e0_rt;
    rnep   = 0.0e0_rt;
    amrex::Real bktinv = 1.0e0_rt/(bk *temp);
    amrex::Real qn     = qn1;

    // chemical potential including the electron rest mass
    amrex::Real etaef = etakep + c2me*bktinv;
    amrex::Real eta, t5, rie1, rie2, rjv1, rjv2;

    // iflag=1 is for electrons,  iflag=2 is for positrons
    for (auto iflag = 1; iflag <= 2; ++iflag) {
        amrex::Real etael;
        if (iflag == 1) {
            etael = qn2*bktinv;
        } else {
            // iflag == 2
            etael = c2me*bktinv;
            etaef = -etaef;
        }

        t5    = temp*temp*temp*temp*temp;
        const amrex::Real zetan = qn*bktinv;
        eta   = etaef - etael;

        // protect from overflowing with large eta values
        const amrex::Real exeta = (eta <= 6.8e2_rt) ? std::exp(eta) : 0.0_rt;
        const amrex::Real etael2 = etael*etael;
        const amrex::Real etael3 = etael2*etael;
        const amrex::Real etael4 = etael3*etael;
        const amrex::Real etael5 = etael4*etael;
        const amrex::Real zetan2 = zetan*zetan;
        const amrex::Real f0 = (eta <= 6.8e2_rt) ? std::log(1.0e0_rt + exeta) : eta;

        // if eta <= 0., the following fermi integrals apply
        const amrex::Real f1l = exeta;
        const amrex::Real f2l = 2.0e0_rt   * f1l;
        const amrex::Real f3l = 6.0e0_rt   * f1l;
        const amrex::Real f4l = 24.0e0_rt  * f1l;
        const amrex::Real f5l = 120.0e0_rt * f1l;

        // if eta > 0., the following fermi integrals apply:
        amrex::Real f1g = 0.0e0_rt;
        amrex::Real f2g = 0.0e0_rt;
        amrex::Real f3g = 0.0e0_rt;
        amrex::Real f4g = 0.0e0_rt;
        amrex::Real f5g = 0.0e0_rt;
        if (eta > 0.0_rt) {
            const amrex::Real exmeta = std::exp(-eta);
            const amrex::Real eta2   = eta*eta;
            const amrex::Real eta3   = eta2*eta;
            const amrex::Real eta4   = eta3*eta;
            f1g = 0.5e0_rt*eta2 + 2.0e0_rt - exmeta;
            f2g = eta3*1.0_rt/3.0_rt + 4.0e0_rt*eta + 2.0e0_rt*exmeta;
            f3g = 0.25e0_rt*eta4 + 0.5e0_rt*pi2*eta2 + 12.0e0_rt - 6.0e0_rt*exmeta;
            f4g = 0.2e0_rt*eta4*eta + 2.0e0_rt*pi2*1.0_rt/3.0_rt*eta3 + 48.0e0_rt*eta
                    + 24.0e0_rt*exmeta;
            f5g = eta4*eta2*1.0_rt/6.0_rt + 5.0e0_rt*1.0_rt/6.0_rt*pi2*eta4
                    + 7.0e0_rt*1.0_rt/6.0_rt*pi2*eta2  + 240.0e0_rt -120.e0_rt*exmeta;
        }

        // factors which are multiplied by the fermi integrals
        const amrex::Real fac3 = 2.0e0_rt*zetan + 4.0e0_rt*etael;
        const amrex::Real fac2 = 6.0e0_rt*etael2 + 6.0e0_rt*etael*zetan + zetan2;
        const amrex::Real fac1 = 4.0e0_rt*etael3 + 6.0e0_rt*etael2*zetan + 2.0e0_rt*etael*zetan2;
        const amrex::Real fac0 = etael4 + 2.0e0_rt*zetan*etael3 + etael2*zetan2;

        // electron capture rates onto protons with no blocking
        rie1 = f4l + fac3*f3l + fac2*f2l + fac1*f1l + fac0*f0;
        rie2 = f4g + fac3*f3g + fac2*f2g + fac1*f1g + fac0*f0;

        // neutrino emission rate for electron capture:
        const amrex::Real facv4 = 5.0e0_rt*etael + 3.0e0_rt*zetan;
        const amrex::Real facv3 = 10.0e0_rt*etael2 + 12.0e0_rt*etael*zetan + 3.0e0_rt*zetan2;
        const amrex::Real facv2 = 10.0e0_rt*etael3 + 18.0e0_rt*etael2*zetan
            + 9.0e0_rt*etael*zetan2 + zetan2*zetan;
        const amrex::Real facv1 = 5.0e0_rt*etael4 + 12.0e0_rt*etael3*zetan
            + 9.0e0_rt*etael2*zetan2 + 2.0e0_rt*etael*zetan2*zetan;
        const amrex::Real facv0 = etael5 + 3.0e0_rt*etael4*zetan
            + 3.0e0_rt*etael3*zetan2 + etael2*zetan2*zetan;
        rjv1  = f5l + facv4*f4l + facv3*f3l
            + facv2*f2l + facv1*f1l + facv0*f0;
        rjv2  = f5g + facv4*f4g + facv3*f3g
            + facv2*f2g + facv1*f1g + facv0*f0;

        // for electrons capture onto protons
        if (iflag < 2) {
            if (eta <= 0.0_rt) {
                rpen  = twoln*cmk5*t5*rie1*ftinv;
                spenc = twoln*cmk6*t5*temp*rjv1*ftinv*c2me;
            } else {
                rpen = twoln*cmk5*t5*rie2*ftinv;
                spenc = twoln*cmk6*t5*temp*rjv2*ftinv*c2me;
            }
            qn = qn2;
        }
    }

    // for positrons capture onto neutrons
    if (eta <= 0.0_rt) {
        rnep  = twoln*cmk5*t5*rie1*ftinv;
        snepc = twoln*cmk6*t5*temp*rjv1*ftinv*c2me;
    } else {
        rnep  = twoln*cmk5*t5*rie2*ftinv;
        snepc = twoln*cmk6*t5*temp*rjv2*ftinv*c2me;
    }
}

#endif
