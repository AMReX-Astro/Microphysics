#ifndef _aprox_rates_H_
#define _aprox_rates_H_

#include <aprox_rates_data.H>
#include <tfactors.H>
#include <cmath>

inline 
void rates_init()
{
    // allocate multi-dimensional array(s)
    datn.resize(2);
    for (auto i = 0; i < 2; ++i) {
        datn[i].resize(6);
    }

    // store rates 
    static const amrex::Real datn0[] = {-4.363_rt, -3.091_rt, -1.275_rt, 1.073_rt, 3.035_rt, 4.825_rt,
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
    static const amrex::Real datn1[] = {-4.539_rt, -3.097_rt, -1.134_rt, 1.525_rt, 3.907_rt, 6.078_rt,
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
void rate_c12ag(tf_t tf, amrex::Real& den, amrex::Real& fr, 
                amrex::Real& dfrdt, amrex::Real& dfrdd, 
                amrex::Real& rr, amrex::Real& drrdt, 
                amrex::Real& drrdd) 
{

}

#endif