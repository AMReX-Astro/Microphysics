#ifndef _aprox_rates_H_
#define _aprox_rates_H_

#include <aprox_rates_data.H>
#include <tfactors.H>
#include <cmath>

using namespace amrex;

inline 
void rates_init()
{
    // allocate multi-dimensional array(s)
    datn.resize(2);
    for (auto i = 0; i < 2; ++i) {
        datn[i].resize(6);
    }

    // store rates 
    static const Real datn0[] = {-4.363_rt, -3.091_rt, -1.275_rt, 1.073_rt, 3.035_rt, 4.825_rt,
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
    static const Real datn1[] = {-4.539_rt, -3.097_rt, -1.134_rt, 1.525_rt, 3.907_rt, 6.078_rt,
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

    for (auto i = 0; i < 6; ++i) {
        datn[0][i].assign(datn0+i*14, datn0+(i+1)*14);
        datn[1][i].assign(datn1+i*14, datn1+(i+1)*14);
    }

    // Evaluate the cubic interp parameters for ni56 electron capture
    // which is used in the langanke subroutine.
    for (auto k = 1; k < 4; ++k) {
        rfdm[k] = 1.0_rt / ((rv[k-1]-rv[k])*(rv[k-1]-rv[k+1])*(rv[k-1]-rv[k+2]));
        rfd0[k] = 1.0_rt / ((rv[k]-rv[k-1])*(rv[k]-rv[k+1])*(rv[k]-rv[k+2]));
        rfd1[k] = 1.0_rt / ((rv[k+1]-rv[k-1])*(rv[k+1]-rv[k])*(rv[k+1]-rv[k+2]));
        rfd2[k] = 1.0_rt / ((rv[k+2]-rv[k-1])*(rv[k+2]-rv[k])*(rv[k+2]-rv[k+1]));
    }

    for (auto j = 1; j < 12; ++j) {
        tfdm[j] = 1.0_rt / ((tv[j-1]-tv[j])*(tv[j-1]-tv[j+1])*(tv[j-1]-tv[j+2]));
        tfd0[j] = 1.0_rt / ((tv[j]-tv[j-1])*(tv[j]-tv[j+1])*(tv[j]-tv[j+2]));
        tfd1[j] = 1.0_rt / ((tv[j+1]-tv[j-1])*(tv[j+1]-tv[j])*(tv[j+1]-tv[j+2]));
        tfd2[j] = 1.0_rt / ((tv[j+2]-tv[j-1])*(tv[j+2]-tv[j])*(tv[j+2]-tv[j+1]));
    }
}

AMREX_GPU_HOST_DEVICE inline 
void rate_c12ag(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    const Real q1 = 1.0e0_rt/12.222016e0_rt;

    // ! c12(a,g)o16
    const Real aa   = 1.0e0_rt + 0.0489e0_rt*tf.t9i23;
    const Real daa  = -twoth*0.0489e0_rt*tf.t9i53

    const Real bb   = tf.t92*aa*aa
    const Real dbb  = 2.0e0_rt*(bb*tf.t9i + tf.t92*aa*daa)

    const Real cc   = exp(-32.120e0_rt*tf.t9i13 - tf.t92*q1)
    const Real dcc  = cc * (oneth*32.120e0_rt*tf.t9i43 - 2.0e0_rt*tf.t9*q1)

    const Real dd   = 1.0e0_rt + 0.2654e0_rt*tf.t9i23
    const Real ddd  = -twoth*0.2654e0_rt*tf.t9i53

    const Real ee   = tf.t92*dd*dd
    const Real dee  = 2.0e0_rt*(ee*tf.t9i + tf.t92*dd*ddd)

    const Real ff   = exp(-32.120e0_rt*tf.t9i13)
    const Real dff  = ff * oneth*32.120e0_rt*tf.t9i43

    const Real gg   = 1.25e3_rt * tf.t9i32 * exp(-27.499_rt*tf.t9i)
    const Real dgg  = gg*(-1.5e0_rt*tf.t9i + 27.499_rt*tf.t9i2)

    const Real hh   = 1.43e-2_rt * tf.t95 * exp(-15.541_rt*tf.t9i)
    const Real dhh  = hh*(5.0e0_rt*tf.t9i + 15.541_rt*tf.t9i2)

    const Real zz   = 1.0e0_rt/bb
    const Real f1   = cc*zz
    const Real df1  = (dcc - f1*dbb)*zz

    const Real zz   = 1.0e0_rt/ee
    const Real f2   = ff*zz
    const Real df2  = (dff - f2*dee)*zz

    Real term = 1.04e8_rt*f1  + 1.76e8_rt*f2 + gg + hh
    Real dtermdt = 1.04e8_rt*df1 + 1.76e8_rt*df2 + dgg + dhh

    // ! 1.7 times cf88 value
    term     *= 1.7e0_rt;
    dtermdt  *= 1.7e0_rt;

    fr    = term * den;
    dfrdt = dtermdt * den * 1.0e-9_rt;
    // !dfrdd = term;

    const Real rev    = 5.13e10_rt * tf.t932 * exp(-83.111_rt*tf.t9i);
    drevdt = rev*(1.5e0_rt*tf.t9i + 83.111_rt*tf.t9i2);

    rr     = rev * term;
    drrdt  = (drevdt*term + rev*dtermdt) * 1.0e-9_rt;
    // !drrdd  = 0.0e0_rt;
}

// ! This routine computes the nuclear reaction rate for 12C(a,g)16O and its inverse 
// ! using fit parameters from Deboer et al. 2017 (https://doi.org/10.1103/RevModPhys.89.035007).
AMREX_GPU_HOST_DEVICE inline 
void rate_c12ag_deboer17(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! from Table XXVI of deboer + 2017
    // ! non-resonant contributions to the reaction
    Real a0_nr = 24.1e0_rt
    Real a1_nr = 0e0_rt 
    Real a2_nr = -32e0_rt
    Real a3_nr = -5.9e0_rt
    Real a4_nr = 1.8e0_rt
    Real a5_nr = -0.17e0_rt
    Real a6_nr = -twoth

    Real term_a0_nr = exp(a0_nr)
    Real term_a1_nr = exp(a1_nr*tf.t9i)
    Real term_a2_nr = exp(a2_nr*tf.t9i13)
    Real term_a3_nr = exp(a3_nr*tf.t913)
    Real term_a4_nr = exp(a4_nr*tf.t9)
    Real term_a5_nr = exp(a5_nr*tf.t953)
    Real term_a6_nr = tf.t9**a6_nr

    Real term_nr = term_a0_nr * term_a1_nr * term_a2_nr * &
              term_a3_nr * term_a4_nr * term_a5_nr * &
              term_a6_nr

    Real dterm_a0_nr = 0e0_rt
    Real dterm_a1_nr = 0e0_rt
    Real dterm_a2_nr = -a2_nr*tf.t9i43*term_a2_nr/3e0_rt
    Real dterm_a3_nr = a3_nr*tf.t9i23*term_a3_nr/3e0_rt
    Real dterm_a4_nr = a4_nr*term_a4_nr
    Real dterm_a5_nr = a5_nr*tf.t923*term_a5_nr*fiveth
    Real dterm_a6_nr = tf.t9i*a6_nr*tf.t9**a6_nr

    Real dterm_nr = (term_a0_nr * term_a1_nr * dterm_a2_nr * term_a3_nr * term_a4_nr * term_a5_nr * term_a6_nr) + &
               (term_a0_nr * term_a1_nr * term_a2_nr * dterm_a3_nr * term_a4_nr * term_a5_nr * term_a6_nr) + &
               (term_a0_nr * term_a1_nr * term_a2_nr * term_a3_nr * dterm_a4_nr * term_a5_nr * term_a6_nr) + &
               (term_a0_nr * term_a1_nr * term_a2_nr * term_a3_nr * term_a4_nr * dterm_a5_nr * term_a6_nr) + &
               (term_a0_nr * term_a1_nr * term_a2_nr * term_a3_nr * term_a4_nr * term_a5_nr * dterm_a6_nr)

    // ! resonant contributions to the reaction
    a0_r = 7.4e0_rt
    a1_r = -30e0_rt
    a2_r = 0e0_rt
    a3_r = 0e0_rt
    a4_r = 0e0_rt
    a5_r = 0e0_rt
    a6_r = -3.0e0_rt/2.0e0_rt

    term_a0_r = exp(a0_r)
    term_a1_r = exp(a1_r*tf.t9i)
    term_a2_r = exp(a2_r*tf.t9i13)
    term_a3_r = exp(a3_r*tf.t913)
    term_a4_r = exp(a4_r*tf.t9)
    term_a5_r = exp(a5_r*tf.t953)
    term_a6_r = tf.t9**a6_r

    term_r = term_a0_r * term_a1_r * term_a2_r * &
              term_a3_r * term_a4_r * term_a5_r * &
              term_a6_r

    dterm_a0_r = 0e0_rt
    dterm_a1_r = -a1_r*tf.t9i2*term_a1_r
    dterm_a2_r = 0e0_rt
    dterm_a3_r = 0e0_rt
    dterm_a4_r = 0e0_rt
    dterm_a5_r = 0e0_rt
    dterm_a6_r = tf.t9i*a6_r*tf.t9**a6_r
    
    dterm_r = (term_a0_r * dterm_a1_r * term_a6_r) + &
              (term_a0_r * term_a1_r * dterm_a6_r)
    

    // ! full rate is the sum of resonant and non-resonant contributions
    term = term_nr + term_r 
    dtermdt = dterm_nr + dterm_r

    fr    = term * den
    dfrdt = dtermdt * den * 1.0e-9_rt

    // ! first term is 9.8685e9_rt * T9**(2/3) * (M0*M1/M3)**(3/2) 
    // ! see iliadis 2007 eqn. 3.44
    // ! ratio of partition functions are assumed to be unity
    const Real rev    = 5.1345573e10_rt * tf.t932 * exp(-83.114082_rt*tf.t9i)
    drevdt = rev*(1.5e0_rt*tf.t9i + 83.114082_rt*tf.t9i2)

    rr     = rev * term
    drrdt  = (drevdt*term + rev*dtermdt) * 1.0e-9_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_tripalf(tf_t tf, Real& den, Real& fr, 
                  Real& dfrdt, Real& dfrdd, 
                  Real& rr, Real& drrdt, 
                  Real& drrdd) 
{
    const Real rc28 = 0.1_rt;
    const Real q1 = 1.0_rt / 0.009604e0_rt;
    const Real q2 = 1.0e0_rt/0.055225e0_rt;

    // ! triple alfa to c12
    // ! this is a(a,g)be8
    aa    = 7.40e+05_rt * tf.t9i32 * exp(-1.0663_rt*tf.t9i)
    daa   = aa*(-1.5e0_rt*tf.t9i  + 1.0663_rt*tf.t9i2)

    bb    = 4.164e+09_rt * tf.t9i23 * exp(-13.49_rt*tf.t9i13 - tf.t92*q1)
    dbb   = bb*(-twoth*tf.t9i + oneth*13.49_rt*tf.t9i43 - 2.0e0_rt*tf.t9*q1)

    cc    = 1.0e0_rt + 0.031_rt*tf.t913 + 8.009_rt*tf.t923 + 1.732_rt*tf.t9 
          + 49.883_rt*tf.t943 + 27.426_rt*tf.t953
    dcc   = oneth*0.031_rt*tf.t9i23 + twoth*8.009_rt*tf.t9i13 + 1.732_rt 
          + fourth*49.883_rt*tf.t913 + fiveth*27.426_rt*tf.t923

    r2abe    = aa + bb * cc
    dr2abedt = daa + dbb*cc + bb*dcc


    // ! this is be8(a,g)c12
    dd    = 130.0e0_rt * tf.t9i32 * exp(-3.3364_rt*tf.t9i)
    ddd   = dd*(-1.5e0_rt*tf.t9i + 3.3364_rt*tf.t9i2)

    ee    = 2.510e+07_rt * tf.t9i23 * exp(-23.57_rt*tf.t9i13 - tf.t92*q2)
    dee   = ee*(-twoth*tf.t9i + oneth*23.57_rt*tf.t9i43 - 2.0e0_rt*tf.t9*q2)

    ff    = 1.0e0_rt + 0.018_rt*tf.t913 + 5.249_rt*tf.t923 + 0.650_rt*tf.t9 + 
         19.176_rt*tf.t943 + 6.034_rt*tf.t953
    dff   = oneth*0.018_rt*tf.t9i23 + twoth*5.249_rt*tf.t9i13 + 0.650_rt 
          + fourth*19.176_rt*tf.t913 + fiveth*6.034_rt*tf.t923

    rbeac    = dd + ee * ff
    drbeacdt = ddd + dee * ff + ee * dff


    // ! a factor
    xx    = rc28 * 1.35e-07_rt * tf.t9i32 * exp(-24.811_rt*tf.t9i)
    dxx   = xx*(-1.5e0_rt*tf.t9i + 24.811_rt*tf.t9i2)


    // ! high temperature rate
    if (tf.t9.gt.0.08_rt) then
        term    = 2.90e-16_rt * r2abe * rbeac + xx
        dtermdt =   2.90e-16_rt * dr2abedt * rbeac 
                    + 2.90e-16_rt * r2abe * drbeacdt 
                    + dxx

    // ! low temperature rate
    else
        uu   = 0.8e0_rt*exp(-(0.025_rt*tf.t9i)**3.263_rt)
        yy   = 0.2e0_rt + uu
        //    ! fxt yy   = 0.01 + 0.2e0_rt + uu
        dyy  = uu * 3.263_rt*(0.025_rt*tf.t9i)**2.263_rt * (0.025_rt*tf.t9i2)
        vv   = 4.0e0_rt*exp(-(tf.t9/0.025_rt)**9.227_rt)
        zz   = 1.0e0_rt + vv
        dzz  = vv * 9.227_rt*(tf.t9/0.025_rt)**8.227_rt * 40.0e0_rt
        aa   = 1.0e0_rt/zz
        f1   = 0.01e0_rt + yy * aa
        //    ! fxt f1   = yy * aa
        df1  = (dyy - f1*dzz)*aa
        term = 2.90e-16_rt * r2abe * rbeac * f1 +  xx
        dtermdt =   2.90e-16_rt * dr2abedt * rbeac * f1 
                    + 2.90e-16_rt * r2abe * drbeacdt * f1 
                    + 2.90e-16_rt * r2abe * rbeac * df1 
                    + dxx
    end if


    // ! rates
    // !      term    = 1.2e0_rt * term
    // !      dtermdt = 1.2e0_rt * term

    fr    = term * den * den
    dfrdt = dtermdt * den * den * 1.0e-9_rt
    !dfrdd = 2.0e0_rt * term * den

    rev    = 2.00e+20_rt*tf.t93*exp(-84.424_rt*tf.t9i)
    drevdt = rev*(3.0e0_rt*tf.t9i + 84.424_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_c12c12(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! c12 + c12 reaction
    const Real aa      = 1.0e0_rt + 0.0396_rt*tf.t9
    const Real zz      = 1.0e0_rt/aa

    const Real t9a     = tf.t9*zz
    const Real dt9a    = (1.0e0_rt -  t9a*0.0396_rt)*zz

    const Real zz      = dt9a/t9a
    const Real t9a13   = t9a**oneth
    const Real dt9a13  = oneth*t9a13*zz

    const Real t9a56   = t9a**fivsix
    const Real dt9a56  = fivsix*t9a56*zz

    term    = 4.27e+26_rt * t9a56 * tf.t9i32 * &
         exp(-84.165_rt/t9a13 - 2.12e-03_rt*tf.t93)
    dtermdt = term*(dt9a56/t9a56 - 1.5e0_rt*tf.t9i &
            + 84.165_rt/t9a13**2*dt9a13 - 6.36e-3_rt*tf.t92)

    // ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    rr    = 0.0e0_rt
    drrdt = 0.0e0_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_c12o16(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! c12 + o16 reaction; see cf88 references 47-4
    if (tf.t9.ge.0.5_rt) then
        Real aa     = 1.0e0_rt + 0.055_rt*tf.t9
        Real zz     = 1.0e0_rt/aa

        Real t9a    = tf.t9*zz
        Real dt9a   = (1.0e0_rt - t9a*0.055_rt)*zz

        zz     = dt9a/t9a
        Real t9a13  = t9a**oneth
        Real dt9a13 = oneth*t9a13*zz

        Real t9a23  = t9a13*t9a13
        Real dt9a23 = 2.0e0_rt * t9a13 * dt9a13

        Real t9a56  = t9a**fivsix
        Real dt9a56 = fivsix*t9a56*zz

        aa      = exp(-0.18_rt*t9a*t9a)
        Real daa     = -aa * 0.36_rt * t9a * dt9a

        Real bb      = 1.06e-03_rt*exp(2.562_rt*t9a23)
        Real dbb     = bb * 2.562_rt * dt9a23

        Real cc      = aa + bb
        Real dcc     = daa + dbb

        zz      = 1.0e0_rt/cc
        term    = 1.72e+31_rt * t9a56 * tf.t9i32 * exp(-106.594_rt/t9a13) * zz
        dtermdt = term*(dt9a56/t9a56 - 1.5e0_rt*tf.t9i 
                        + 106.594_rt/t9a23*dt9a13 - zz*dcc)

    else
        // ! term    = 2.6288035e-29_rt
        term    = 0.0e0_rt
        dtermdt = 0.0e0_rt
    endif


    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    rr    = 0.0e0_rt
    drrdt = 0.0e0_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_o16o16(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! o16 + o16
    term  = 7.10e36_rt * tf.t9i23 * 
         exp(-135.93_rt * tf.t9i13 - 0.629_rt*tf.t923 
         - 0.445_rt*tf.t943 + 0.0103_rt*tf.t9*tf.t9)

    dtermdt = -twoth*term*tf.t9i 
         + term * (oneth*135.93_rt*tf.t9i43 - twoth*0.629_rt*tf.t9i13 
         - fourth*0.445_rt*tf.t913 + 0.0206_rt*tf.t9)


    // ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    rr    = 0.0e0_rt
    drrdt = 0.0e0_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_o16ag(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! o16(a,g)ne20
    const Real term1   = 9.37e9_rt * tf.t9i23 * exp(-39.757_rt*tf.t9i13 - tf.t92*q1)
    const Real dterm1  = term1*(-twoth*tf.t9i + oneth*39.757_rt*tf.t9i43 - 2.0e0_rt*tf.t9*q1)

    const Real aa      = 62.1_rt * tf.t9i32 * exp(-10.297_rt*tf.t9i)
    const Real daa     = aa*(-1.5e0_rt*tf.t9i + 10.297_rt*tf.t9i2)

    const Real bb      = 538.0e0_rt * tf.t9i32 * exp(-12.226_rt*tf.t9i)
    const Real dbb     = bb*(-1.5e0_rt*tf.t9i + 12.226_rt*tf.t9i2)

    const Real cc      = 13.0e0_rt * tf.t92 * exp(-20.093_rt*tf.t9i)
    const Real dcc     = cc*(2.0e0_rt*tf.t9i + 20.093_rt*tf.t9i2)

    const Real term2   = aa + bb + cc
    const Real dterm2  = daa + dbb + dcc

    term    = term1 + term2
    dtermdt = dterm1 + dterm2

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    rev      = 5.65e+10_rt*tf.t932*exp(-54.937_rt*tf.t9i)
    drevdt   = rev*(1.5e0_rt*tf.t9i + 54.937_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_ne20ag(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! ne20(a,g)mg24
    Real aa   = 4.11e+11_rt * tf.t9i23 * exp(-46.766_rt*tf.t9i13 - tf.t92*q1)
    Real daa  = aa*(-twoth*tf.t9i + oneth*46.766_rt*tf.t9i43 - 2.0e0_rt*tf.t9*q1)

    Real bb   = 1.0e0_rt + 0.009_rt*tf.t913 + 0.882_rt*tf.t923 + 0.055_rt*tf.t9 &
         + 0.749_rt*tf.t943 + 0.119_rt*tf.t953
    Real dbb  = oneth*0.009_rt*tf.t9i23 + twoth*0.882_rt*tf.t9i13 + 0.055_rt &
         + fourth*0.749_rt*tf.t913 + fiveth*0.119_rt*tf.t923

    Real term1  = aa * bb
    Real dterm1 = daa * bb + aa * dbb

    aa   = 5.27e+03_rt * tf.t9i32 * exp(-15.869_rt*tf.t9i)
    daa  = aa*(-1.5e0_rt*tf.t9i + 15.869_rt*tf.t9i2)

    bb   = 6.51e+03_rt * tf.t912 * exp(-16.223_rt*tf.t9i)
    dbb  = bb*(0.5e0_rt*tf.t9i + 16.223_rt*tf.t9i2)

    Real term2  = aa + bb
    Real dterm2 = daa + dbb

    aa   = 42.1_rt * tf.t9i32 * exp(-9.115_rt*tf.t9i)
    daa  = aa*(-1.5e0_rt*tf.t9i + 9.115_rt*tf.t9i2)

    bb   =  32.0_rt * tf.t9i23 * exp(-9.383_rt*tf.t9i)
    dbb  = bb*(-twoth*tf.t9i + 9.383_rt*tf.t9i2)

    Real term3  = rc102 * (aa + bb)
    Real dterm3 = rc102 * (daa + dbb)

    aa  = 5.0e0_rt*exp(-18.960_rt*tf.t9i)
    daa = aa*18.960_rt*tf.t9i2

    bb  = 1.0e0_rt + aa
    dbb = daa

    Real zz      = 1.0e0_rt/bb
    term    = (term1 + term2 + term3)*zz
    dtermdt = ((dterm1 + dterm2 + dterm3) - term*dbb)*zz

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    Real rev      = 6.01e+10_rt * tf.t932 * exp(-108.059_rt*tf.t9i)
    drevdt   = rev*(1.5e0_rt*tf.t9i + 108.059_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_mg24ag(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! 24mg(a,g)28si
    const Real aa    = 4.78e+01_rt * tf.t9i32 * exp(-13.506_rt*tf.t9i)
    const Real daa   = aa*(-1.5e0_rt*tf.t9i + 13.506_rt*tf.t9i2)

    const Real bb    =  2.38e+03_rt * tf.t9i32 * exp(-15.218_rt*tf.t9i)
    const Real dbb   = bb*(-1.5e0_rt*tf.t9i + 15.218_rt*tf.t9i2)

    const Real cc    = 2.47e+02_rt * tf.t932 * exp(-15.147_rt*tf.t9i)
    const Real dcc   = cc*(1.5e0_rt*tf.t9i + 15.147_rt*tf.t9i2)

    const Real dd    = rc121 * 1.72e-09_rt * tf.t9i32 * exp(-5.028_rt*tf.t9i)
    const Real ddd   = dd*(-1.5e0_rt*tf.t9i + 5.028_rt*tf.t9i2)

    const Real ee    = rc121* 1.25e-03_rt * tf.t9i32 * exp(-7.929_rt*tf.t9i)
    const Real dee   = ee*(-1.5e0_rt*tf.t9i + 7.929_rt*tf.t9i2)

    const Real ff    = rc121 * 2.43e+01_rt * tf.t9i * exp(-11.523_rt*tf.t9i)
    const Real dff   = ff*(-tf.t9i + 11.523_rt*tf.t9i2)

    const Real gg    = 5.0e0_rt*exp(-15.882_rt*tf.t9i)
    const Real dgg   = gg*15.882_rt*tf.t9i2

    const Real hh    = 1.0e0_rt + gg
    const Real hhi   = 1.0e0_rt/hh

    term    = (aa + bb + cc + dd + ee + ff) * hhi
    dtermdt = (daa + dbb + dcc + ddd + dee + dff - term*dgg) * hhi

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 6.27e+10_rt * tf.t932 * exp(-115.862_rt*tf.t9i)
    drevdt   = rev*(1.5e0_rt*tf.t9i + 115.862_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_mg24ap(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    const Real rc148 = 0.1_rt;
    const Real q1 = 1.0_rt / 0.024649e0_rt;

    // ! 24mg(a,p)al27
    Real aa     = 1.10e+08_rt * tf.t9i23 * exp(-23.261_rt*tf.t9i13 - tf.t92*q1)
    Real daa    = -twoth*aa*tf.t9i + aa*(23.261_rt*tf.t9i43 - 2.0e0_rt*tf.t9*q1)

    Real bb     =  1.0e0_rt + 0.018_rt*tf.t913 + 12.85_rt*tf.t923 + 1.61_rt*tf.t9 &
         + 89.87_rt*tf.t943 + 28.66_rt*tf.t953
    Real dbb    = oneth*0.018_rt*tf.t9i23 + twoth*12.85_rt*tf.t9i13 + 1.61_rt &
           + fourth*89.87_rt*tf.t913 + fiveth*28.66_rt*tf.t923

    Real term1  = aa * bb
    Real dterm1 = daa * bb + aa * dbb

    aa     = 129.0e0_rt * tf.t9i32 * exp(-2.517_rt*tf.t9i)
    daa    = -1.5e0_rt*aa*tf.t9i + aa*2.517_rt*tf.t9i2

    bb     = 5660.0e0_rt * tf.t972 * exp(-3.421_rt*tf.t9i)
    dbb    = 3.5e0_rt*bb*tf.t9i +  bb*3.421_rt*tf.t9i2

    Real cc     = rc148 * 3.89e-08_rt * tf.t9i32 * exp(-0.853_rt*tf.t9i)
    Real dcc    = -1.5e0_rt*cc*tf.t9i + cc*0.853_rt*tf.t9i2

    Real dd     = rc148 * 8.18e-09_rt * tf.t9i32 * exp(-1.001_rt*tf.t9i)
    Real ddd    = -1.5e0_rt*dd*tf.t9i + dd*1.001_rt*tf.t9i2

    Real term2  = aa + bb + cc + dd
    Real dterm2 = daa + dbb + dcc + ddd

    Real ee     = oneth*exp(-9.792_rt*tf.t9i)
    Real dee    = ee*9.792_rt*tf.t9i2

    Real ff     =  twoth * exp(-11.773_rt*tf.t9i)
    Real dff    = ff*11.773_rt*tf.t9i2

    Real gg     = 1.0e0_rt + ee + ff
    Real dgg    = dee + dff

    term    = (term1 + term2)/gg
    dtermdt = ((dterm1 + dterm2) - term*dgg)/gg

    // ! the rates
    Real rev      = 1.81_rt * exp(-18.572_rt*tf.t9i)
    Real drevdt   = rev*18.572_rt*tf.t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !dfrdd = rev * term

    rr    = den * term
    drrdt = den * dtermdt * 1.0e-9_rt
    // !drrdd = term
}

AMREX_GPU_HOST_DEVICE inline 
void rate_al27pg(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! al27(p,g)si28
    // ! champagne 1996

    const Real aa  = 1.32e+09_rt * tf.t9i23 * exp(-23.26_rt*tf.t9i13)
    const Real daa = aa*(-twoth*tf.t9i + oneth*23.26_rt*tf.t9i43)

    const Real bb  = 3.22e-10_rt * tf.t9i32 * exp(-0.836_rt*tf.t9i)*0.17_rt
    const Real dbb = bb*(-1.5e0_rt*tf.t9i + 0.836_rt*tf.t9i2)

    const Real cc  = 1.74e+00_rt * tf.t9i32 * exp(-2.269_rt*tf.t9i)
    const Real dcc = cc*(-1.5e0_rt*tf.t9i + 2.269_rt*tf.t9i2)

    const Real dd  = 9.92e+00_rt * tf.t9i32 * exp(-2.492_rt*tf.t9i)
    const Real ddd = dd*(-1.5e0_rt*tf.t9i + 2.492_rt*tf.t9i2)

    const Real ee  = 4.29e+01_rt * tf.t9i32 * exp(-3.273_rt*tf.t9i)
    const Real dee = ee*(-1.5e0_rt*tf.t9i + 3.273_rt*tf.t9i2)

    const Real ff  = 1.34e+02_rt * tf.t9i32 * exp(-3.654_rt*tf.t9i)
    const Real dff = ff*(-1.5e0_rt*tf.t9i + 3.654_rt*tf.t9i2)

    const Real gg  = 1.77e+04_rt * (tf.t9**0.53_rt) * exp(-4.588_rt*tf.t9i)
    const Real dgg = gg*(0.53_rt*tf.t9i + 4.588_rt*tf.t9i2)

    term    = aa + bb + cc + dd + ee + ff + gg
    dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg

    // ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 1.13e+11_rt * tf.t932 * exp(-134.434_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 134.434_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_al27pg_old(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    const Real rc147 = 0.1_rt;
    const Real q1 = 1.0e0_rt/0.024025e0_rt;

    // ! 27al(p,g)si28  cf88
    const Real aa  = 1.67e+08_rt * tf.t9i23 * exp(-23.261_rt*tf.t9i13 - tf.t92*q1)
    const Real daa = aa*(-twoth*tf.t9i + oneth*23.261_rt*tf.t9i43 - 2.0e0_rt*tf.t9*q1)

    const Real bb  = 1.0e0_rt + 0.018_rt*tf.t913 + 5.81_rt*tf.t923 + 0.728_rt*tf.t9 &
         + 27.31_rt*tf.t943 + 8.71_rt*tf.t953
    const Real dbb = oneth*0.018_rt*tf.t9i23 + twoth*5.81_rt*tf.t9i13 + 0.728_rt &
         + fourth*27.31_rt*tf.t913 + fiveth*8.71_rt*tf.t923

    const Real cc  = aa*bb
    const Real dcc = daa*bb + aa*dbb

    const Real dd  = 2.20e+00_rt * tf.t9i32 * exp(-2.269_rt*tf.t9i)
    const Real ddd = dd*(-1.5e0_rt*tf.t9i + 2.269_rt*tf.t9i2)

    const Real ee  = 1.22e+01_rt * tf.t9i32 * exp(-2.491_rt*tf.t9i)
    const Real dee = ee*(-1.5e0_rt*tf.t9i + 2.491_rt*tf.t9i2)

    const Real ff  =  1.50e+04_rt * tf.t9 * exp(-4.112_rt*tf.t9i)
    const Real dff = ff*(tf.t9i + 4.112_rt*tf.t9i2)

    const Real gg  = rc147 * 6.50e-10_rt * tf.t9i32 * exp(-0.853_rt*tf.t9i)
    const Real dgg = gg*(-1.5e0_rt*tf.t9i + 0.853_rt*tf.t9i2)

    const Real hh  = rc147 * 1.63e-10_rt * tf.t9i32 * exp(-1.001_rt*tf.t9i)
    const Real dhh = hh*(-1.5e0_rt*tf.t9i + 1.001_rt*tf.t9i2)

    const Real xx     = oneth*exp(-9.792_rt*tf.t9i)
    const Real dxx    = xx*9.792_rt*tf.t9i2

    const Real yy     =  twoth * exp(-11.773_rt*tf.t9i)
    const Real dyy    = yy*11.773_rt*tf.t9i2

    const Real zz     = 1.0e0_rt + xx + yy
    const Real dzz    = dxx + dyy

    const Real pp      = 1.0e0_rt/zz
    const Real term    = (cc + dd + ee + ff + gg + hh)*pp
    dtermdt = ((dcc + ddd + dee + dff + dgg + dhh) - term*dzz)*pp


    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 1.13e+11_rt*tf.t932*exp(-134.434_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 134.434_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_si28ag(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! si28(a,g)s32
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 6.340e-2_rt*z + 2.541e-3_rt*z2 - 2.900e-4_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 6.340e-2_rt + 2.0e0_rt*2.541e-3_rt*tf.t9 - 3.0e0_rt*2.900e-4_rt*tf.t92
    // if (z .eq. 10.0_rt) then
    //    daa = 0
    // else
    //    daa   = 6.340e-2_rt + 2.0e0_rt*2.541e-3_rt*tf.t9 - 3.0e0_rt*2.900e-4_rt*tf.t92
    // end if

    const Real term    = 4.82e+22_rt * tf.t9i23 * exp(-61.015_rt * tf.t9i13 * aa)
    const Real dtermdt = term*(-twoth*tf.t9i + 61.015_rt*tf.t9i13*(oneth*tf.t9i*aa - daa))
  
    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 6.461e+10_rt * tf.t932 * exp(-80.643_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 80.643_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_si28ap(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! si28(a,p)p31
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 2.798e-3_rt*z + 2.763e-3_rt*z2 - 2.341e-4_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 2.798e-3_rt + 2.0e0_rt*2.763e-3_rt*tf.t9 - 3.0e0_rt*2.341e-4_rt*tf.t92;
    // if (z .eq. 10.0_rt) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 2.798e-3_rt + 2.0e0_rt*2.763e-3_rt*tf.t9 - 3.0e0_rt*2.341e-4_rt*tf.t92
    // end if

    const Real term    = 4.16e+13_rt * tf.t9i23 * exp(-25.631_rt * tf.t9i13 * aa)
    const Real dtermdt = -twoth*term*tf.t9i + term*25.631_rt*tf.t9i13*(oneth*tf.t9i*aa - daa)


    // ! the rates
    const Real rev      = 0.5825e0_rt * exp(-22.224_rt*tf.t9i)
    const Real drevdt   = rev*22.224_rt*tf.t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !dfrdd = rev * term

    rr    = den * term
    drrdt = den * dtermdt * 1.0e-9_rt
    // !drrdd = term
}

AMREX_GPU_HOST_DEVICE inline 
void rate_p31pg(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! p31(p,g)s32
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 1.928e-1_rt*z - 1.540e-2_rt*z2 + 6.444e-4_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 1.928e-1_rt - 2.0e0_rt*1.540e-2_rt*tf.t9 + 3.0e0_rt*6.444e-4_rt*tf.t92;
    // if (z .eq. 10.0_rt) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 1.928e-1_rt - 2.0e0_rt*1.540e-2_rt*tf.t9 + 3.0e0_rt*6.444e-4_rt*tf.t92
    // end if

    const Real term    = 1.08e+16_rt * tf.t9i23 * exp(-27.042_rt * tf.t9i13 * aa)
    const Real dtermdt = term*(-twoth*tf.t9i + 27.042_rt*tf.t9i13*(oneth*tf.t9i*aa - daa))

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 3.764e+10_rt * tf.t932 * exp(-102.865_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 102.865_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_s32ag(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! s32(a,g)ar36
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 4.913e-2_rt*z + 4.637e-3_rt*z2 - 4.067e-4_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 4.913e-2_rt + 2.0e0_rt*4.637e-3_rt*tf.t9 - 3.0e0_rt*4.067e-4_rt*tf.t92;
    // if (z .eq. 10.0_rt) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 4.913e-2_rt + 2.0e0_rt*4.637e-3_rt*tf.t9 - 3.0e0_rt*4.067e-4_rt*tf.t92
    // end if

    const Real term    = 1.16e+24_rt * tf.t9i23 * exp(-66.690_rt * tf.t9i13 * aa)
    const Real dtermdt = term*(-twoth*tf.t9i + 66.690_rt*tf.t9i13*(oneth*tf.t9i*aa - daa))


    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 6.616e+10_rt * tf.t932 * exp(-77.080_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 77.080_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_s32ap(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! s32(a,p)cl35
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 1.041e-1_rt*z - 1.368e-2_rt*z2 + 6.969e-4_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 1.041e-1_rt - 2.0e0_rt*1.368e-2_rt*tf.t9 + 3.0e0_rt*6.969e-4_rt*tf.t92;
    // if (z .eq. 10) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 1.041e-1_rt - 2.0e0_rt*1.368e-2_rt*tf.t9 + 3.0e0_rt*6.969e-4_rt*tf.t92
    // end if

    const Real term    = 1.27e+16_rt * tf.t9i23 * exp(-31.044_rt * tf.t9i13 * aa)
    const Real dtermdt = -twoth*term*tf.t9i + term*31.044_rt*tf.t9i13*(oneth*tf.t9i*aa - daa)

    // ! the rates
    const Real rev      = 1.144_rt * exp(-21.643_rt*tf.t9i)
    const Real drevdt   = rev*21.643_rt*tf.t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    // !dfrdd = rev * term

    rr    = den * term
    drrdt = den * dtermdt * 1.0e-9_rt
    // !drrdd = term
}

AMREX_GPU_HOST_DEVICE inline 
void rate_cl35pg(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! cl35(p,g)ar36
    const Real aa    = 1.0e0_rt + 1.761e-1_rt*tf.t9 - 1.322e-2_rt*tf.t92 + 5.245e-4_rt*tf.t93
    const Real daa   = 1.761e-1_rt - 2.0e0_rt*1.322e-2_rt*tf.t9 + 3.0e0_rt*5.245e-4_rt*tf.t92
  
    const Real term    =  4.48e+16_rt * tf.t9i23 * exp(-29.483_rt * tf.t9i13 * aa)
    const Real dtermdt = term*(-twoth*tf.t9i + 29.483_rt*tf.t9i13*(oneth*tf.t9i*aa - daa))

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 7.568e+10_rt*tf.t932*exp(-98.722_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 98.722_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_ar36ag(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! ar36(a,g)ca40
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 1.458e-1_rt*z - 1.069e-2_rt*z2 + 3.790e-4_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 1.458e-1_rt - 2.0e0_rt*1.069e-2_rt*tf.t9 + 3.0e0_rt*3.790e-4_rt*tf.t92;
    // if (z .eq. 10.0_rt) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 1.458e-1_rt - 2.0e0_rt*1.069e-2_rt*tf.t9 + 3.0e0_rt*3.790e-4_rt*tf.t92
    // end if

    const Real term    = 2.81e+30_rt * tf.t9i23 * exp(-78.271_rt * tf.t9i13 * aa)
    const Real dtermdt = term*(-twoth*tf.t9i + 78.271_rt*tf.t9i13*(oneth*tf.t9i*aa - daa))

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 6.740e+10_rt * tf.t932 * exp(-81.711_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 81.711_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_ar36ap(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! ar36(a,p)k39
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 4.826e-3_rt*z - 5.534e-3_rt*z2 + 4.021e-4_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 4.826e-3_rt - 2.0e0_rt*5.534e-3_rt*tf.t9 + 3.0e0_rt*4.021e-4_rt*tf.t92;
    // if (z .eq. 10.0_rt) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 4.826e-3_rt - 2.0e0_rt*5.534e-3_rt*tf.t9 + 3.0e0_rt*4.021e-4_rt*tf.t92
    // end if

    const Real term    = 2.76e+13_rt * tf.t9i23 * exp(-34.922_rt * tf.t9i13 * aa)
    const Real dtermdt = -twoth*term*tf.t9i + term*34.922_rt*tf.t9i13*(oneth*tf.t9i*aa - daa)

    // ! the rates
    const Real rev      = 1.128_rt*exp(-14.959_rt*tf.t9i)
    const Real drevdt   = rev*14.959_rt*tf.t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    // !dfrdd = rev * term

    rr    = den * term
    drrdt = den * dtermdt * 1.0e-9_rt
    // !drrdd = term
}

AMREX_GPU_HOST_DEVICE inline 
void rate_k39pg(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! k39(p,g)ca40
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 1.622e-1_rt*z - 1.119e-2_rt*z2 + 3.910e-4_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 1.622e-1_rt - 2.0e0_rt*1.119e-2_rt*tf.t9 + 3.0e0_rt*3.910e-4_rt*tf.t92;
    // if (z .eq. 10) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 1.622e-1_rt - 2.0e0_rt*1.119e-2_rt*tf.t9 + 3.0e0_rt*3.910e-4_rt*tf.t92
    // end if

    const Real term    = 4.09e+16_rt * tf.t9i23 * exp(-31.727_rt * tf.t9i13 * aa)
    const Real dtermdt = term*(-twoth*tf.t9i + 31.727_rt*tf.t9i13*(oneth*tf.t9i*aa - daa))

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 7.600e+10_rt * tf.t932 * exp(-96.657_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 96.657_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_ca40ag(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! ca40(a,g)ti44
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 1.650e-2_rt*z + 5.973e-3_rt*z2 - 3.889e-04_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 1.650e-2_rt + 2.0e0_rt*5.973e-3_rt*tf.t9 - 3.0e0_rt*3.889e-4_rt*tf.t92;
    // if (z .eq. 10.0_rt) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 1.650e-2_rt + 2.0e0_rt*5.973e-3_rt*tf.t9 - 3.0e0_rt*3.889e-4_rt*tf.t92
    // end if

    const Real term    = 4.66e+24_rt * tf.t9i23 * exp(-76.435_rt * tf.t9i13 * aa)
    const Real dtermdt = term*(-twoth*tf.t9i + 76.435_rt*tf.t9i13*(oneth*tf.t9i*aa - daa))


    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 6.843e+10_rt * tf.t932 * exp(-59.510_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 59.510_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_ca40ap(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! ca40(a,p)sc43
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt - 1.206e-2_rt*z + 7.753e-3_rt*z2 - 5.071e-4_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : -1.206e-2_rt + 2.0e0_rt*7.753e-3_rt*tf.t9 - 3.0e0_rt*5.071e-4_rt*tf.t92;
    // if (z .eq. 10.0_rt) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = -1.206e-2_rt + 2.0e0_rt*7.753e-3_rt*tf.t9 - 3.0e0_rt*5.071e-4_rt*tf.t92
    // end if

    const Real term    = 4.54e+14_rt * tf.t9i23 * exp(-32.177_rt * tf.t9i13 * aa)
    const Real dtermdt = -twoth*term*tf.t9i + term*32.177_rt*tf.t9i13*(oneth*tf.t9i*aa - daa)


    // ! the rates
    const Real rev      = 2.229_rt * exp(-40.966_rt*tf.t9i)
    const Real drevdt   = rev*40.966_rt*tf.t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    // !dfrdd = rev * term

    rr    = den * term
    drrdt = den * dtermdt * 1.0e-9_rt
    // !drrdd = term
}

AMREX_GPU_HOST_DEVICE inline 
void rate_sc43pg(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! sc43(p,g)ca40
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 1.023e-1_rt*z - 2.242e-3_rt*z2 - 5.463e-5_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 1.023e-1_rt - 2.0e0_rt*2.242e-3_rt*tf.t9 - 3.0e0_rt*5.463e-5_rt*tf.t92;
    // if (z .eq. 10.0_rt) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 1.023e-1_rt - 2.0e0_rt*2.242e-3_rt*tf.t9 - 3.0e0_rt*5.463e-5_rt*tf.t92
    // end if

    const Real term    = 3.85e+16_rt * tf.t9i23 * exp(-33.234_rt * tf.t9i13 * aa)
    const Real dtermdt = term*(-twoth*tf.t9i + 33.234_rt*tf.t9i13*(oneth*tf.t9i*aa - daa))

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 1.525e+11_rt * tf.t932 * exp(-100.475_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 100.475_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_ti44ag(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! ti44(a,g)cr48
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 1.066e-1_rt*z - 1.102e-2_rt*z2 + 5.324e-4_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 1.066e-1_rt - 2.0e0_rt*1.102e-2_rt*tf.t9 + 3.0e0_rt*5.324e-4_rt*tf.t92
    // if (z .eq. 10.0_rt) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 1.066e-1_rt - 2.0e0_rt*1.102e-2_rt*tf.t9 + 3.0e0_rt*5.324e-4_rt*tf.t92
    // end if

    const Real term    = 1.37e+26_rt * tf.t9i23 * exp(-81.227_rt * tf.t9i13 * aa)
    const Real dtermdt = term*(-twoth*tf.t9i + 81.227_rt*tf.t9i13*(oneth*tf.t9i*aa - daa))

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 6.928e+10_rt*tf.t932*exp(-89.289_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 89.289_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_ti44ap(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! ti44(a,p)v47
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 2.655e-2_rt*z - 3.947e-3_rt*z2 + 2.522e-4_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 2.655e-2_rt - 2.0e0_rt*3.947e-3_rt*tf.t9 + 3.0e0_rt*2.522e-4_rt*tf.t92;
    // if (z .eq. 10.0_rt) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 2.655e-2_rt - 2.0e0_rt*3.947e-3_rt*tf.t9 + 3.0e0_rt*2.522e-4_rt*tf.t92
    // end if

    const Real term    = 6.54e+20_rt * tf.t9i23 * exp(-66.678_rt * tf.t9i13 * aa)
    const Real dtermdt = -twoth*term*tf.t9i + term*66.678_rt*tf.t9i13*(oneth*tf.t9i*aa - daa)

    // ! the rates
    const Real rev      = 1.104_rt * exp(-4.723_rt*tf.t9i)
    const Real drevdt   = rev*4.723_rt*tf.t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    // !dfrdd = rev * term

    rr    = den * term
    drrdt = den * dtermdt * 1.0e-9_rt
    // !drrdd = term
}

AMREX_GPU_HOST_DEVICE inline 
void rate_v47pg(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! v47(p,g)cr48
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 9.979e-2_rt*z - 2.269e-3_rt*z2 - 6.662e-5_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0 : 9.979e-2_rt - 2.0e0_rt*2.269e-3_rt*tf.t9 - 3.0e0_rt*6.662e-5_rt*tf.t92;
    // if (z .eq. 10.0_rt) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 9.979e-2_rt - 2.0e0_rt*2.269e-3_rt*tf.t9 - 3.0e0_rt*6.662e-5_rt*tf.t92
    // end if

    const Real term    = 2.05e+17_rt * tf.t9i23 * exp(-35.568_rt * tf.t9i13 * aa)
    const Real dtermdt = term*(-twoth*tf.t9i + 35.568_rt*tf.t9i13*(oneth*tf.t9i*aa - daa))

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 7.649e+10_rt*tf.t932*exp(-93.999_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 93.999_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_cr48ag(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! cr48(a,g)fe52
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 6.325e-2_rt*z - 5.671e-3_rt*z2 + 2.848e-4_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 6.325e-2_rt - 2.0e0_rt*5.671e-3_rt*tf.t9 + 3.0e0_rt*2.848e-4_rt*tf.t92
    // if (z .eq. 10.0_rt) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 6.325e-2_rt - 2.0e0_rt*5.671e-3_rt*tf.t9 + 3.0e0_rt*2.848e-4_rt*tf.t92
    // end if

    const Real term    = 1.04e+23_rt * tf.t9i23 * exp(-81.420_rt * tf.t9i13 * aa)
    const Real dtermdt = term*(-twoth*tf.t9i + 81.420_rt*tf.t9i13*(oneth*tf.t9i*aa - daa))

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 7.001e+10_rt * tf.t932 * exp(-92.177_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 92.177_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_cr48ap(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! cr48(a,p)mn51
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 1.384e-2_rt*z + 1.081e-3_rt*z2 - 5.933e-5_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 1.384e-2_rt + 2.0e0_rt*1.081e-3_rt*tf.t9 - 3.0e0_rt*5.933e-5_rt*tf.t92;
    // if (z .eq. 10.0_rt) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 1.384e-2_rt + 2.0e0_rt*1.081e-3_rt*tf.t9 - 3.0e0_rt*5.933e-5_rt*tf.t92
    // end if

    const Real term    = 1.83e+26_rt * tf.t9i23 * exp(-86.741_rt * tf.t9i13 * aa)
    const Real dtermdt = -twoth*term*tf.t9i + term*86.741_rt*tf.t9i13*(oneth*tf.t9i*aa - daa)

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 0.6087_rt*exp(-6.510_rt*tf.t9i)
    const Real drevdt   = rev*6.510_rt*tf.t9i2

    rr    = den * rev * term
    drrdt = den * (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    // !drrdd = rev * term
}

AMREX_GPU_HOST_DEVICE inline 
void rate_mn51pg(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! mn51(p,g)fe52
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 8.922e-2_rt*z - 1.256e-3_rt*z2 - 9.453e-5_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 8.922e-2_rt - 2.0e0_rt*1.256e-3_rt*tf.t9 - 3.0e0_rt*9.453e-5_rt*tf.t92;
    // if (z .eq. 10.0_rt) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 8.922e-2_rt - 2.0e0_rt*1.256e-3_rt*tf.t9 - 3.0e0_rt*9.453e-5_rt*tf.t92
    // end if

    const Real term    = 3.77e+17_rt * tf.t9i23 * exp(-37.516_rt * tf.t9i13 * aa)
    const Real dtermdt = term*(-twoth*tf.t9i + 37.516_rt*tf.t9i13*(oneth*tf.t9i*aa - daa))

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 1.150e+11_rt*tf.t932*exp(-85.667_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 85.667_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_fe52ag(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! fe52(a,g)ni56
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 7.846e-2_rt*z - 7.430e-3_rt*z2 + 3.723e-4_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 7.846e-2_rt - 2.0e0_rt*7.430e-3_rt*tf.t9 + 3.0e0_rt*3.723e-4_rt*tf.t92;
    // if (z .eq. 10.0_rt) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 7.846e-2_rt - 2.0e0_rt*7.430e-3_rt*tf.t9 + 3.0e0_rt*3.723e-4_rt*tf.t92
    // end if

    const Real term    = 1.05e+27_rt * tf.t9i23 * exp(-91.674_rt * tf.t9i13 * aa)
    const Real dtermdt = term*(-twoth*tf.t9i + 91.674_rt*tf.t9i13*(oneth*tf.t9i*aa - daa))

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 7.064e+10_rt*tf.t932*exp(-92.850_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 92.850_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_fe52ap(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! fe52(a,p)co55
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 1.367e-2_rt*z + 7.428e-4_rt*z2 - 3.050e-5_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 1.367e-2_rt + 2.0e0_rt*7.428e-4_rt*tf.t9 - 3.0e0_rt*3.050e-5_rt*tf.t92;
    // if (z .eq. 10.0_rt) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 1.367e-2_rt + 2.0e0_rt*7.428e-4_rt*tf.t9 - 3.0e0_rt*3.050e-5_rt*tf.t92
    // end if

    const Real term    = 1.30e+27_rt * tf.t9i23 * exp(-91.674_rt * tf.t9i13 * aa)
    const Real dtermdt = -twoth*term*tf.t9i + term*91.674_rt*tf.t9i13*(oneth*tf.t9i*aa - daa)

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 0.4597_rt*exp(-9.470_rt*tf.t9i)
    const Real drevdt   = rev*9.470_rt*tf.t9i2

    rr    = den * rev * term
    drrdt = den * (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    // !drrdd = rev * term
}

AMREX_GPU_HOST_DEVICE inline 
void rate_co55pg(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! co55(p,g)ni56
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 9.894e-2_rt*z - 3.131e-3_rt*z2 - 2.160e-5_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 9.894e-2_rt - 2.0e0_rt*3.131e-3_rt*tf.t9 - 3.0e0_rt*2.160e-5_rt*tf.t92;
    // if (z .eq. 10.0_rt) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 9.894e-2_rt - 2.0e0_rt*3.131e-3_rt*tf.t9 - 3.0e0_rt*2.160e-5_rt*tf.t92
    // end if

    const Real term    = 1.21e+18_rt * tf.t9i23 * exp(-39.604_rt * tf.t9i13 * aa)
    const Real dtermdt = term*(-twoth*tf.t9i + 39.604_rt*tf.t9i13*(oneth*tf.t9i*aa - daa))

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 1.537e+11_rt*tf.t932*exp(-83.382_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 83.382_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_pp(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! p(p,e+nu)d
    Real term, dtermdt;
    if (tf.t9 .le. 3.0_rt) then
       const Real aa   = 4.01e-15_rt * tf.t9i23 * exp(-3.380e0_rt*tf.t9i13)
       const Real daa  = aa*(-twoth*tf.t9i + oneth*3.380e0_rt*tf.t9i43)

       const Real bb   = 1.0e0_rt + 0.123e0_rt*tf.t913 + 1.09e0_rt*tf.t923 + 0.938e0_rt*tf.t9
       const Real dbb  = oneth*0.123e0_rt*tf.t9i23 + twoth*1.09e0_rt*tf.t9i13 + 0.938e0_rt

       term    = aa * bb
       dtermdt = daa * bb + aa * dbb

    else
       term    = 1.1581136e-15_rt
       dtermdt = 0.0e0_rt
    end if

    // ! rate
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    rr    = 0.0e0_rt
    drrdt = 0.0e0_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_png(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! p(n,g)d
    // ! smith,kawano,malany 1992

    const Real aa      = 1.0e0_rt - 0.8504_rt*tf.t912 + 0.4895_rt*tf.t9 \
         - 0.09623_rt*tf.t932 + 8.471e-3_rt*tf.t92 \
         - 2.80e-4_rt*tf.t952

    const Real daa     =  -0.5e0_rt*0.8504_rt*tf.t9i12 + 0.4895_rt \
         - 1.5e0_rt*0.09623_rt*tf.t912 + 2.0e0_rt*8.471e-3_rt*tf.t9 \
         - 2.5e0_rt*2.80e-4_rt*tf.t932

    const Real term    = 4.742e4_rt * aa
    const Real dtermdt = 4.742e4_rt * daa

    // ! wagoner,schramm 1977
    // !      aa      = 1.0e0_rt - 0.86*tf.t912 + 0.429*tf.t9
    // !      daa     =  -0.5e0_rt*0.86*tf.t9i12 + 0.429

    // !      term    = 4.4e4_rt * aa
    // !      dtermdt = 4.4e4_rt * daa

    // ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 4.71e+09_rt * tf.t932 * exp(-25.82_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 25.82_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_dpg(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! d(p,g)he3
    const Real aa      = 2.24e+03_rt * tf.t9i23 * exp(-3.720_rt*tf.t9i13)
    const Real daa     = aa*(-twoth*tf.t9i + oneth*3.720_rt*tf.t9i43)

    const Real bb      = 1.0e0_rt + 0.112_rt*tf.t913 + 3.38_rt*tf.t923 + 2.65_rt*tf.t9
    const Real dbb     = oneth*0.112_rt*tf.t9i23 + twoth*3.38_rt*tf.t9i13 + 2.65_rt

    const Real term    = aa * bb
    const Real dtermdt = daa * bb + aa * dbb

    // ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 1.63e+10_rt * tf.t932 * exp(-63.750_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 63.750_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_he3ng(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! he3(n,g)he4
    const Real term    = 6.62_rt * (1.0e0_rt + 905.0_rt*tf.t9)
    const Real dtermdt = 5.9911e3_rt

    // ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 2.61e+10_rt * tf.t932 * exp(-238.81_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 238.81_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_he3he3(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! he3(he3,2p)he4
    const Real aa   = 6.04e+10_rt * tf.t9i23 * exp(-12.276_rt*tf.t9i13)
    const Real daa  = aa*(-twoth*tf.t9i + oneth*12.276_rt*tf.t9i43)

    const Real bb   = 1.0e0_rt + 0.034_rt*tf.t913 - 0.522_rt*tf.t923 - 0.124_rt*tf.t9 &
         + 0.353_rt*tf.t943 + 0.213_rt*tf.t953
    const Real dbb  = oneth*0.034_rt*tf.t9i23 - twoth*0.522_rt*tf.t9i13 - 0.124_rt &
         + fourth*0.353_rt*tf.t913 + fiveth*0.213_rt*tf.t923

    const Real term    = aa * bb
    const Real dtermdt = daa*bb + aa*dbb

    // ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 3.39e-10_rt * tf.t9i32 * exp(-149.230_rt*tf.t9i)
    const Real drevdt   = rev*(-1.5e0_rt*tf.t9i + 149.230_rt*tf.t9i2)

    rr    = den * den * rev * term
    drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    // !drrdd = 2.0e0_rt * den * rev * term
}

AMREX_GPU_HOST_DEVICE inline 
void rate_he3he4(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! he3(he4,g)be7
    const Real aa      = 1.0e0_rt + 0.0495_rt*tf.t9
    const Real daa     = 0.0495_rt

    const Real zz      = 1.0e0_rt/aa
    const Real t9a     = tf.t9*zz
    const Real dt9a    = (1.0e0_rt - t9a*daa)*zz

    const Real zz      = dt9a/t9a
    const Real t9a13   = t9a**oneth
    const Real dt9a13  = oneth*t9a13*zz

    const Real t9a56   = t9a**fivsix
    const Real dt9a56  = fivsix*t9a56*zz

    const Real term    = 5.61e+6_rt * t9a56 * tf.t9i32 * exp(-12.826_rt/t9a13)
    const Real dtermdt = term*(dt9a56/t9a56 - 1.5e0_rt*tf.t9i &
         + 12.826_rt/t9a13**2 * dt9a13)

    // ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 1.11e+10_rt * tf.t932 * exp(-18.423_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 18.423_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_c12pg(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! c12(p,g)13n
    const Real aa   = 2.04e+07_rt * tf.t9i23 * exp(-13.69_rt*tf.t9i13 - tf.t92*q1)
    const Real daa  = aa*(-twoth*tf.t9i + oneth*13.69_rt*tf.t9i43 - 2.0e0_rt*tf.t9*q1)

    const Real bb   = 1.0e0_rt + 0.03_rt*tf.t913 + 1.19_rt*tf.t923 + 0.254_rt*tf.t9 \
         + 2.06_rt*tf.t943 + 1.12_rt*tf.t953
    const Real dbb  = oneth*0.03_rt*tf.t9i23 + twoth*1.19_rt*tf.t9i13 + 0.254_rt \
         + fourth*2.06_rt*tf.t913 + fiveth*1.12_rt*tf.t923

    const Real cc   = aa * bb
    const Real dcc  = daa*bb + aa*dbb

    const Real dd   = 1.08e+05_rt * tf.t9i32 * exp(-4.925_rt*tf.t9i)
    const Real ddd  = dd*(-1.5e0_rt*tf.t9i + 4.925_rt*tf.t9i2)

    const Real ee   = 2.15e+05_rt * tf.t9i32 * exp(-18.179_rt*tf.t9i)
    const Real dee  = ee*(-1.5e0_rt*tf.t9i + 18.179_rt*tf.t9i2)

    const Real term    = cc + dd + ee
    const Real dtermdt = dcc + ddd + dee

    // ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    dfrdd = term

    const Real rev      = 8.84e+09_rt * tf.t932 * exp(-22.553_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 22.553_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_n14pg(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! n14(p,g)o15
    const Real aa  = 4.90e+07_rt * tf.t9i23 * exp(-15.228_rt*tf.t9i13 - tf.t92*q1)
    const Real daa = aa*(-twoth*tf.t9i + oneth*15.228_rt*tf.t9i43 - 2.0e0_rt*tf.t9*q1)

    const Real bb   = 1.0e0_rt + 0.027_rt*tf.t913 - 0.778_rt*tf.t923 - 0.149_rt*tf.t9 &
         + 0.261_rt*tf.t943 + 0.127_rt*tf.t953
    const Real dbb  = oneth*0.027_rt*tf.t9i23 - twoth*0.778_rt*tf.t9i13 - 0.149_rt &
         + fourth*0.261_rt*tf.t913 + fiveth*0.127_rt*tf.t923

    const Real cc   = aa * bb
    const Real dcc  = daa*bb + aa*dbb

    const Real dd   = 2.37e+03_rt * tf.t9i32 * exp(-3.011_rt*tf.t9i)
    const Real ddd  = dd*(-1.5e0_rt*tf.t9i + 3.011_rt*tf.t9i2)

    const Real ee   = 2.19e+04_rt * exp(-12.530_rt*tf.t9i)
    const Real dee  = ee*12.530_rt*tf.t9i2

    const Real term    = cc + dd + ee
    const Real dtermdt = dcc + ddd + dee

    // ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev    = 2.70e+10_rt * tf.t932 * exp(-84.678_rt*tf.t9i)
    const Real drevdt = rev*(1.5e0_rt*tf.t9i + 84.678_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_n15pg(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    const Real q1 = 1.0_rt / 0.2025_rt;

    // ! n15(p,g)o16
    const Real aa  = 9.78e+08_rt * tf.t9i23 * exp(-15.251_rt*tf.t9i13 - tf.t92*q1)
    const Real daa = aa*(-twoth*tf.t9i + oneth*15.251_rt*tf.t9i43 - 2.0e0_rt*tf.t9*q1)

    const Real bb   = 1.0e0_rt  + 0.027_rt*tf.t913 + 0.219_rt*tf.t923 + 0.042_rt*tf.t9 \
         + 6.83_rt*tf.t943 + 3.32_rt*tf.t953
    const Real dbb  = oneth*0.027_rt*tf.t9i23 + twoth*0.219_rt*tf.t9i13 + 0.042_rt \
         + fourth*6.83_rt*tf.t913 + fiveth*3.32_rt*tf.t923

    const Real cc   = aa * bb
    const Real dcc  = daa*bb + aa*dbb

    const Real dd   = 1.11e+04_rt*tf.t9i32*exp(-3.328_rt*tf.t9i)
    const Real ddd  = dd*(-1.5e0_rt*tf.t9i + 3.328_rt*tf.t9i2)

    const Real ee   = 1.49e+04_rt*tf.t9i32*exp(-4.665_rt*tf.t9i)
    const Real dee  = ee*(-1.5e0_rt*tf.t9i + 4.665_rt*tf.t9i2)

    const Real ff   = 3.8e+06_rt*tf.t9i32*exp(-11.048_rt*tf.t9i)
    const Real dff  = ff*(-1.5e0_rt*tf.t9i + 11.048_rt*tf.t9i2)

    const Real term    = cc + dd + ee + ff
    const Real dtermdt = dcc + ddd + dee + dff

    // ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 3.62e+10_rt * tf.t932 * exp(-140.734_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 140.734_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_n15pa(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    const Real theta = 0.1_rt;
    const Real q1 = 1.0_rt / 0.272484_rt;

    // ! n15(p,a)c12
    const Real aa  = 1.08e+12_rt*tf.t9i23*exp(-15.251_rt*tf.t9i13 - tf.t92*q1)
    const Real daa = aa*(-twoth*tf.t9i + oneth*15.251_rt*tf.t9i43 - 2.0e0_rt*tf.t9*q1)

    const Real bb   = 1.0e0_rt + 0.027_rt*tf.t913 + 2.62_rt*tf.t923 + 0.501_rt*tf.t9 &
         + 5.36_rt*tf.t943 + 2.60_rt*tf.t953
    const Real dbb  = oneth*0.027_rt*tf.t9i23 + twoth*2.62_rt*tf.t9i13 + 0.501_rt &
         + fourth*5.36_rt*tf.t913 + fiveth*2.60_rt*tf.t923

    const Real cc   = aa * bb
    const Real dcc  = daa*bb + aa*dbb

    const Real dd   = 1.19e+08_rt * tf.t9i32 * exp(-3.676_rt*tf.t9i)
    const Real ddd  = dd*(-1.5e0_rt*tf.t9i + 3.676_rt*tf.t9i2)

    const Real ee   = 5.41e+08_rt * tf.t9i12 * exp(-8.926_rt*tf.t9i)
    const Real dee  = ee*(-0.5e0_rt*tf.t9i + 8.926_rt*tf.t9i2)

    const Real ff   = theta * 4.72e+08_rt * tf.t9i32 * exp(-7.721_rt*tf.t9i)
    const Real dff  = ff*(-1.5e0_rt*tf.t9i + 7.721_rt*tf.t9i2)

    const Real gg   = theta * 2.20e+09_rt * tf.t9i32 * exp(-11.418_rt*tf.t9i)
    const Real dgg  = gg*(-1.5e0_rt*tf.t9i + 11.418_rt*tf.t9i2)

    const Real term    = cc + dd + ee + ff + gg
    const Real dtermdt = dcc + ddd + dee + dff + dgg

    // ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 7.06e-01_rt*exp(-57.625_rt*tf.t9i)
    const Real drevdt   = rev*57.625_rt*tf.t9i2

    rr    = den * rev * term
    drrdt = den * (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    // !drrdd = rev * term
}

AMREX_GPU_HOST_DEVICE inline 
void rate_o16pg(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! o16(p,g)f17
    const Real aa  = exp(-0.728_rt*tf.t923)
    const Real daa = -twoth*aa*0.728_rt*tf.t9i13

    const Real bb  = 1.0e0_rt + 2.13_rt * (1.0e0_rt - aa)
    const Real dbb = -2.13_rt*daa

    const Real cc  = tf.t923 * bb
    const Real dcc = twoth*cc*tf.t9i + tf.t923*dbb

    const Real dd   = exp(-16.692_rt*tf.t9i13)
    const Real ddd  = oneth*dd*16.692_rt*tf.t9i43

    const Real zz   = 1.0e0_rt/cc
    const Real ee   = dd*zz
    const Real dee  = (ddd - ee*dcc)*zz

    const Real term    = 1.50e+08_rt * ee
    const Real dtermdt = 1.50e+08_rt * dee

    // ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 3.03e+09_rt*tf.t932*exp(-6.968_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 6.968_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_n14ag(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    const Real q1 = 1.0e0_rt/0.776161e0_rt;

    // ! n14(a,g)f18
    const Real aa  = 7.78e+09_rt * tf.t9i23 * exp(-36.031_rt*tf.t9i13- tf.t92*q1)
    const Real daa = aa*(-twoth*tf.t9i + oneth*36.031_rt*tf.t9i43 - 2.0e0_rt*tf.t9*q1)

    const Real bb   = 1.0e0_rt + 0.012_rt*tf.t913 + 1.45_rt*tf.t923 + 0.117_rt*tf.t9 \
         + 1.97_rt*tf.t943 + 0.406_rt*tf.t953
    const Real dbb  = oneth*0.012_rt*tf.t9i23 + twoth*1.45_rt*tf.t9i13 + 0.117_rt \
         + fourth*1.97_rt*tf.t913 + fiveth*0.406_rt*tf.t923

    const Real cc   = aa * bb
    const Real dcc  = daa*bb + aa*dbb

    const Real dd   = 2.36e-10_rt * tf.t9i32 * exp(-2.798_rt*tf.t9i)
    const Real ddd  = dd*(-1.5e0_rt*tf.t9i + 2.798_rt*tf.t9i2)

    const Real ee   = 2.03_rt * tf.t9i32 * exp(-5.054_rt*tf.t9i)
    const Real dee  = ee*(-1.5e0_rt*tf.t9i + 5.054_rt*tf.t9i2)

    const Real ff   = 1.15e+04_rt * tf.t9i23 * exp(-12.310_rt*tf.t9i)
    const Real dff  = ff*(-twoth*tf.t9i + 12.310_rt*tf.t9i2)

    const Real term    = cc + dd + ee + ff
    const Real dtermdt = dcc + ddd + dee + dff

    // ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 5.42e+10_rt * tf.t932 * exp(-51.236_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 51.236_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_fe52ng(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! fe52(n,g)fe53
    const Real tq2     = tf.t9 - 0.348e0_rt
    const Real term    = 9.604e+05_rt * exp(-0.0626_rt*tq2)
    const Real dtermdt = -term*0.0626_rt

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 2.43e+09_rt * tf.t932 * exp(-123.951_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 123.951_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_fe53ng(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! fe53(n,g)fe54
    const Real tq1   = tf.t9/0.348_rt
    const Real tq10  = tq1**0.10_rt
    const Real dtq10 = 0.1e0_rt*tq10/(0.348_rt*tq1)
    const Real tq2   = tf.t9 - 0.348e0_rt

    const Real term    = 1.817e+06_rt * tq10 * exp(-0.06319_rt*tq2)
    const Real dtermdt = term/tq10*dtq10 - term*0.06319_rt

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    const Real rev      = 1.56e+11_rt * tf.t932 * exp(-155.284_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 155.284_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_fe54ng(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! fe54(n,g)fe55
    const Real aa   =  2.307390e+01_rt - 7.931795e-02_rt * tf.t9i + 7.535681e+00_rt * tf.t9i13 \
         - 1.595025e+01_rt * tf.t913 + 1.377715e+00_rt * tf.t9 - 1.291479e-01_rt * tf.t953 \
         + 6.707473e+00_rt * log(tf.t9)

    const Real daa  =  7.931795e-02_rt * tf.t9i2 - oneth * 7.535681e+00_rt * tf.t9i43 \
         - oneth * 1.595025e+01_rt *tf.t9i23 + 1.377715e+00_rt - fiveth * 1.291479e-01_rt *tf.t923 \
         + 6.707473e+00_rt * tf.t9i

    Real term, dtermdt;
    if (aa .lt. 200.0_rt) then
       term    = exp(aa)
       dtermdt = term*daa*1.0e-9_rt
    else
       term    = exp(200.0e0_rt)
       dtermdt = 0.0e0_rt
    end if

    const Real bb  = 4.800293e+09_rt * tf.t932 * exp(-1.078986e+02_rt * tf.t9i)
    const Real dbb = bb*(1.5e0_rt*tf.t9i + 1.078986e+02_rt * tf.t9i2)

    // ! reverse rate
    rr    = term*bb
    drrdt = dtermdt*bb + term*dbb*1.0e-9_rt
    // !drrdd = 0.0e0_rt

    // ! forward rate
    // !dfrdd = term
    fr    = term*den
    dfrdt = dtermdt*den
}

AMREX_GPU_HOST_DEVICE inline 
void rate_fe54pg(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! fe54(p,g)co55
    const Real z     = min(tf.t9,10.0e0_rt)
    const Real z2    = z*z
    const Real z3    = z2*z
    const Real aa    = 1.0e0_rt + 9.593e-2_rt*z - 3.445e-3_rt*z2 + 8.594e-5_rt*z3
    const Real daa = (z == 10.0_rt) ? 0.0_rt : 9.593e-2_rt - 2.0e0_rt*3.445e-3_rt*tf.t9 + 3.0e0_rt*8.594e-5_rt*tf.t92;
    // if (z .eq. 10.0_rt) then
    //    daa = 0.0e0_rt
    // else
    //    daa   = 9.593e-2_rt - 2.0e0_rt*3.445e-3_rt*tf.t9 + 3.0e0_rt*8.594e-5_rt*tf.t92
    // end if

    const Real term    = 4.51e+17_rt * tf.t9i23 * exp(-38.483_rt * tf.t9i13 * aa)
    const Real dtermdt = term*(-twoth*tf.t9i + 38.483_rt*tf.t9i13*(oneth*tf.t9i*aa - daa))

    // ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    // !dfrdd = term

    rconst Real ev      = 2.400e+09_rt * tf.t932 * exp(-58.605_rt*tf.t9i)
    const Real drevdt   = rev*(1.5e0_rt*tf.t9i + 58.605_rt*tf.t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    // !drrdd = 0.0e0_rt
}

AMREX_GPU_HOST_DEVICE inline 
void rate_fe54ap(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! fe54(a,p)co57
    const Real aa   =  3.97474900e+01_rt - 6.06543100e+00_rt * tf.t9i + 1.63239600e+02_rt * tf.t9i13 \
         - 2.20457700e+02_rt * tf.t913 + 8.63980400e+00_rt * tf.t9 - 3.45841300e-01_rt * tf.t953 \
         + 1.31464200e+02_rt * log(tf.t9)

    const Real daa  =  6.06543100e+00_rt * tf.t9i2 - oneth * 1.63239600e+02_rt * tf.t9i43 \
         - oneth * 2.20457700e+02_rt * tf.t9i23 + 8.63980400e+00_rt - fiveth * 3.45841300e-01_rt * tf.t923 \
         + 1.31464200e+02_rt  * tf.t9i

    Real term, dtermdt;
    if (aa .lt. 200.0_rt) then
       term    = exp(aa)
       dtermdt = term*daa*1.0e-9_rt
    else
       term    = exp(200.0e0_rt)
       dtermdt = 0.0e0_rt
    end if

    const Real bb  = 2.16896000e+00_rt  * exp(-2.05631700e+01_rt * tf.t9i)
    const Real dbb = bb * 2.05631700e+01_rt * tf.t9i2

    // ! reverse rate
    // !drrdd = term
    rr    = term*den
    drrdt = dtermdt*den

    // ! forward rate
    fr    = rr*bb
    dfrdt = drrdt*bb + rr*dbb*1.0e-9_rt
    // !dfrdd = drrdd*bb
}

AMREX_GPU_HOST_DEVICE inline 
void rate_fe55ng(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! fe55(n,g)fe56
    const Real aa   =  1.954115e+01_rt - 6.834029e-02_rt * tf.t9i + 5.379859e+00_rt * tf.t9i13 \
         - 8.758150e+00_rt * tf.t913 + 5.285107e-01_rt * tf.t9 - 4.973739e-02_rt  * tf.t953 \
         + 4.065564e+00_rt  * log(tf.t9)

    const Real daa  =  6.834029e-02_rt * tf.t9i2 - oneth * 5.379859e+00_rt * tf.t9i43 \
         - oneth * 8.758150e+00_rt * tf.t9i23 + 5.285107e-01_rt - fiveth * 4.973739e-02_rt  *tf.t923 \
         + 4.065564e+00_rt  * tf.t9i

    Real term, dtermdt;
    if (aa .lt. 200.0_rt) then
       term    = exp(aa)
       dtermdt = term*daa*1.0e-9_rt
    else
       term    = exp(200.0e0_rt)
       dtermdt = 0.0e0_rt
    end if

    const Real bb  = 7.684279e+10_rt  * tf.t932 * exp(-1.299472e+02_rt  * tf.t9i)
    const Real dbb = bb*(1.5e0_rt*tf.t9i + 1.299472e+02_rt * tf.t9i2)

    // ! reverse rate
    rr    = term*bb
    drrdt = dtermdt*bb + term*dbb*1.0e-9_rt
    // !drrdd = 0.0e0_rt

    // ! forward rate
    // !dfrdd = term
    fr    = term*den
    dfrdt = dtermdt*den
}

AMREX_GPU_HOST_DEVICE inline 
void rate_fe56pg(tf_t tf, Real& den, Real& fr, 
                Real& dfrdt, Real& dfrdd, 
                Real& rr, Real& drrdt, 
                Real& drrdd) 
{
    // ! fe56(p,g)co57

    const Real aa   =  1.755960e+02_rt - 7.018872e+00_rt * tf.t9i + 2.800131e+02_rt * tf.t9i13 \
         - 4.749343e+02_rt * tf.t913 + 2.683860e+01_rt * tf.t9 - 1.542324e+00_rt  * tf.t953 \
         + 2.315911e+02_rt  * log(tf.t9)

    const Real daa  =  7.018872e+00_rt * tf.t9i2 - oneth * 2.800131e+02_rt * tf.t9i43 \
         - oneth * 4.749343e+02_rt * tf.t9i23 + 2.683860e+01_rt - fiveth * 1.542324e+00_rt  *tf.t923 \
         + 2.315911e+02_rt  * tf.t9i

    Real term, dtermdt;
    // if (aa .lt. 200.0_rt) then
    //    term    = exp(aa)
    //    dtermdt = term*daa*1.0e-9_rt
    // else
    //    term    = exp(200.0e0_rt)
    //    dtermdt = 0.0e0_rt
    // end if

    const Real bb  = 2.402486e+09_rt * tf.t932 * exp(-6.995192e+01_rt * tf.t9i)
    const Real dbb = bb*(1.5e0_rt*tf.t9i + 6.995192e+01_rt * tf.t9i2)

    // ! reverse rate
    rr    = term*bb
    drrdt = dtermdt*bb + term*dbb*1.0e-9_rt
    // !drrdd = 0.0e0_rt

    // ! forward rate
    // !dfrdd = term
    fr    = term*den
    dfrdt = dtermdt*den
}


// ! this routine evaluates Langanke et al. 2000 fits for the ni56 electron
// ! capture rate rn56ec and neutrino loss rate sn56ec

// ! input:
// ! y56 = nickel56 molar abundance
// ! ye  = electron to baryon number, zbar/abar

// ! output:
// ! rn56ec = ni56 electron capture rate
// ! sn56ec = ni56 neutrino loss rate
AMREX_GPU_HOST_DEVICE inline 
void langanke(const Real btemp, const Real bden, 
              const Real y56, const Real ye, 
              Real& rn56ec, Real& sn56ec) 
{
    // ! calculate ni56 electron capture and neutrino loss rates
    rn56ec = 0.0_rt
    sn56ec = 0.0_rt
    if ( (btemp .lt. 1.0e9_rt) .or. (bden*ye .lt. 1.0e6_rt)) return
    t9    = min(btemp,1.4e10_rt) * 1.0e-9_rt
    r     = max(6.0e0_rt,min(11.0e0_rt,log10(bden*ye)))
    jp    = min(max(2,int(t9)),12)
    kp    = min(max(2,int(r)-5),4)
    rfm   = r - rv(kp-1)
    rf0   = r - rv(kp)
    rf1   = r - rv(kp+1)
    rf2   = r - rv(kp+2)
    dfacm = rf0*rf1*rf2*rfdm(kp)
    dfac0 = rfm*rf1*rf2*rfd0(kp)
    dfac1 = rfm*rf0*rf2*rfd1(kp)
    dfac2 = rfm*rf0*rf1*rfd2(kp)
    tfm   = t9 - tv(jp-1)
    tf0   = t9 - tv(jp)
    tf1   = t9 - tv(jp+1)
    tf2   = t9 - tv(jp+2)
    tfacm = tf0*tf1*tf2*tfdm(jp)
    tfac0 = tfm*tf1*tf2*tfd0(jp)
    tfac1 = tfm*tf0*tf2*tfd1(jp)
    tfac2 = tfm*tf0*tf1*tfd2(jp)

    // ! evaluate the spline fits
    do jr = 1,2
       do jd = jp-1,jp+2
          rne(jr,jd) =   dfacm*datn(jr,kp-1,jd) + dfac0*datn(jr,kp,jd) &
               + dfac1*datn(jr,kp+1,jd) + dfac2*datn(jr,kp+2,jd)
       enddo
       rnt(jr) =  tfacm*rne(jr,jp-1) + tfac0*rne(jr,jp) &
            + tfac1*rne(jr,jp+1) + tfac2*rne(jr,jp+2)
    enddo

    // ! set the output
    rn56ec = 10.0e0_rt**rnt(1)
    sn56ec = 6.022548e+23_rt * 1.60218e-6_rt * y56 * 10.0e0_rt**rnt(2)
   
    // !write(*,*) "btemp",btemp, "bden", bden, "t9",t9,"r",r, "rn56ec",rn56ec, "sn56ec", sn56ec
}

// ! given the electron degeneracy parameter etakep (chemical potential
// ! without the electron's rest mass divided by kt) and the temperature
// ! temp, this routine calculates rates for
// ! electron capture on protons rpen (captures/sec/proton),
// ! positron capture on neutrons rnep (captures/sec/neutron),
// ! and their associated neutrino energy loss rates
// ! spenc (erg/sec/proton) and snepc (erg/sec/neutron)
AMREX_GPU_HOST_DEVICE inline 
void ecapnuc(Real& etakep, Real& temp, 
             Real& rpen, Real& rnep, 
             Real& spenc, Real& snepc) 
{
    const Real qn1    = -2.0716446e-06_rt
    const Real ftinv  = 1.0e0_rt/1083.9269e0_rt
    const Real twoln  = 0.6931472e0_rt
    const Real cmk5   = 1.3635675e-49_rt
    const Real cmk6   = 2.2993864e-59_rt
    const Real bk     = 1.38062e-16_rt
    const Real qn2    = 2.0716446e-06_rt
    const Real c2me   = 8.1872665e-07_rt
    const Real xmp    = 1.6726485e-24_rt
    const Real xmn    = 1.6749543e-24_rt
    const Real qndeca = 1.2533036e-06_rt
    const Real tmean  = 886.7e0_rt

    // ! tmean and qndeca are the mean lifetime and decay energy of the neutron
    // ! xmp,xnp are masses of the p and n in grams.
    // ! c2me is the constant used to convert the neutrino energy
    // ! loss rate from mec2/s (as in the paper) to ergs/particle/sec.

    // ! initialize
    rpen   = 0.0e0_rt
    rnep   = 0.0e0_rt
    Real spen   = 0.0e0_rt
    Real snep   = 0.0e0_rt
    Real t9     = temp * 1.0e-9_rt
    Real bktinv = 1.0e0_rt/(bk *temp)
    int iflag  = 0;
    Real qn     = qn1

    // ! chemical potential including the electron rest mass
    Real etaef = etakep + c2me*bktinv

    // ! iflag=1 is for electrons,  iflag=2 is for positrons
    for (auto iflag = 1; iflag <= 2; ++iflag) {
// 502 iflag = iflag + 1
        // iflag++;
        if (iflag == 1) {
            etael = qn2*bktinv
        } else if (iflag == 2) {
            etael = c2me*bktinv
            etaef = -etaef
        }

        t5    = temp*temp*temp*temp*temp
        zetan = qn*bktinv
        eta   = etaef - etael

        // ! protect from overflowing with large eta values
        if (eta .le. 6.8e2_rt) {
            exeta = exp(eta)
        } else {
            exeta = 0.0e0_rt
        }
        etael2 = etael*etael
        etael3 = etael2*etael
        etael4 = etael3*etael
        etael5 = etael4*etael
        zetan2 = zetan*zetan
        if (eta .le. 6.8e2_rt) {
            f0 = log(1.0e0_rt + exeta)
        } else {
            f0 = eta
        }

        // ! if eta le. 0., the following fermi integrals apply
        f1l = exeta
        f2l = 2.0e0_rt   * f1l
        f3l = 6.0e0_rt   * f1l
        f4l = 24.0e0_rt  * f1l
        f5l = 120.0e0_rt * f1l

        // ! if eta gt. 0., the following fermi integrals apply:
        f1g = 0.0e0_rt
        f2g = 0.0e0_rt
        f3g = 0.0e0_rt
        f4g = 0.0e0_rt
        f5g = 0.0e0_rt
        if (eta .gt. 0.0_rt) {
            exmeta = exp(-eta)
            eta2   = eta*eta
            eta3   = eta2*eta
            eta4   = eta3*eta
            f1g = 0.5e0_rt*eta2 + 2.0e0_rt - exmeta
            f2g = eta3*third + 4.0e0_rt*eta + 2.0e0_rt*exmeta
            f3g = 0.25e0_rt*eta4 + 0.5e0_rt*pi2*eta2 + 12.0e0_rt - 6.0e0_rt*exmeta
            f4g = 0.2e0_rt*eta4*eta + 2.0e0_rt*pi2*third*eta3 + 48.0e0_rt*eta \
                    + 24.0e0_rt*exmeta
            f5g = eta4*eta2*sixth + 5.0e0_rt*sixth*pi2*eta4 \
                    + 7.0e0_rt*sixth*pi2*eta2  + 240.0e0_rt -120.e0_rt*exmeta
        }

        // ! factors which are multiplied by the fermi integrals
        fac3 = 2.0e0_rt*zetan + 4.0e0_rt*etael
        fac2 = 6.0e0_rt*etael2 + 6.0e0_rt*etael*zetan + zetan2
        fac1 = 4.0e0_rt*etael3 + 6.0e0_rt*etael2*zetan + 2.0e0_rt*etael*zetan2
        fac0 = etael4 + 2.0e0_rt*zetan*etael3 + etael2*zetan2

        // ! electron capture rates onto protons with no blocking
        rie1 = f4l + fac3*f3l + fac2*f2l + fac1*f1l + fac0*f0
        rie2 = f4g + fac3*f3g + fac2*f2g + fac1*f1g + fac0*f0

        // ! neutrino emission rate for electron capture:
        facv4 = 5.0e0_rt*etael + 3.0e0_rt*zetan
        facv3 = 10.0e0_rt*etael2 + 12.0e0_rt*etael*zetan + 3.0e0_rt*zetan2
        facv2 = 10.0e0_rt*etael3 + 18.0e0_rt*etael2*zetan \
            + 9.0e0_rt*etael*zetan2 + zetan2*zetan
        facv1 = 5.0e0_rt*etael4 + 12.0e0_rt*etael3*zetan \
            + 9.0e0_rt*etael2*zetan2 + 2.0e0_rt*etael*zetan2*zetan
        facv0 = etael5 + 3.0e0_rt*etael4*zetan \
            + 3.0e0_rt*etael3*zetan2 + etael2*zetan2*zetan
        rjv1  = f5l + facv4*f4l + facv3*f3l \
            + facv2*f2l + facv1*f1l + facv0*f0
        rjv2  = f5g + facv4*f4g + facv3*f3g \
            + facv2*f2g + facv1*f1g + facv0*f0

        // ! for electrons capture onto protons
        // if (iflag.eq.2) go to 503
        if (iflag < 2) {
            // if (eta.gt.0.) go to 505
            if (eta <= 0.0_rt) {
            // rpen  = twoln*cmk5*t5*rie1*ftinv
            // spen  = twoln*cmk6*t5*temp*rjv1*ftinv
            // spenc = twoln*cmk6*t5*temp*rjv1*ftinv*c2me
            // // go to 504
            } else {
                // 505 rpen = twoln*cmk5*t5*rie2*ftinv
                rpen = twoln*cmk5*t5*rie2*ftinv
                spen = twoln*cmk6*t5*temp*rjv2*ftinv
                spenc = twoln*cmk6*t5*temp*rjv2*ftinv*c2me
            }
            // 504 continue
            qn = qn2
            // go to 502
        }
    }

    // ! for positrons capture onto neutrons
// 503 if (eta.gt.0.) go to 507
    if (eta <= 0.0_rt) {
    rnep  = twoln*cmk5*t5*rie1*ftinv
    snep  = twoln*cmk6*t5*temp*rjv1*ftinv
    snepc = twoln*cmk6*t5*temp*rjv1*ftinv*c2me
    // go to 506
    } else {
// 507 rnep  = twoln*cmk5*t5*rie2*ftinv
    rnep  = twoln*cmk5*t5*rie2*ftinv
    snep  = twoln*cmk6*t5*temp*rjv2*ftinv
    snepc = twoln*cmk6*t5*temp*rjv2*ftinv*c2me
    }
// 506 continue
    // return
}

#endif