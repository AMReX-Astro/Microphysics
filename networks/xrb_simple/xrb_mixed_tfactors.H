#ifndef rprox_tfactors_H
#define rprox_tfactors_H

struct tf_t {
    Real t9m1;
    Real t9m13;
    Real t9m23;
    Real t9m32;
    Real t9m3;
    Real t913;
    Real t923;
    Real t943;
    Real t953;
    Real t9113;
    Real t9log;
};

AMREX_GPU_HOST_DEVICE AMREX_INLINE
tf_t calc_tfactors (const Real& T9)
{
    tf_t tfactors;

    tfactors.t9m1 = 1.0_rt / t9;
    tfactors.t9m13 = std::pow(t9m1, 1.0_rt/3.0_rt);
    tfactors.t9m23 = t9m13 * t9m13;
    tfactors.t9m32 = t9m1 * std::sqrt(t9m1);
    tfactors.t9m3 = t9m1 * t9m1 * t9m1;
    tfactors.t913 = std::pow(t9, 1.0_rt/3.0_rt);
    tfactors.t923 = t913 * t913;
    tfactors.t943 = t923 * t923;
    tfactors.t953 = t943 * t913;
    tfactors.t9113 = t953 * t953 * t913;
    tfactors.t9log = std::log(t9);

    return tfactors;
}

#endif
