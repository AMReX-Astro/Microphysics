#ifndef TFACTORS_H
#define TFACTORS_H

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <microphysics_autodiff.H>
#include <cmath>

using namespace amrex::literals;

template <typename dual_t>
struct tf_t {
    dual_t temp;
    dual_t t9;
    dual_t t92;
    dual_t t93;
    // dual_t t94;
    dual_t t95;
    // dual_t t96;
    dual_t t912;
    dual_t t932;
    dual_t t952;
    dual_t t972;
    dual_t t913;
    dual_t t923;
    dual_t t943;
    dual_t t953;
    // dual_t t973;
    // dual_t t9113;
    // dual_t t914;
    // dual_t t934;
    // dual_t t954;
    // dual_t t974;
    // dual_t t915;
    // dual_t t935;
    // dual_t t945;
    // dual_t t965;
    // dual_t t917;
    // dual_t t927;
    // dual_t t947;
    // dual_t t918;
    // dual_t t938;
    // dual_t t958;
    dual_t t9i;
    dual_t t9i2;
    // dual_t t9i3;
    dual_t t9i12;
    dual_t t9i32;
    // dual_t t9i52;
    // dual_t t9i72;
    dual_t t9i13;
    dual_t t9i23;
    dual_t t9i43;
    dual_t t9i53;
    // dual_t t9i14;
    // dual_t t9i34;
    // dual_t t9i54;
    // dual_t t9i15;
    // dual_t t9i35;
    // dual_t t9i45;
    // dual_t t9i65;
    // dual_t t9i17;
    // dual_t t9i27;
    // dual_t t9i47;
    // dual_t t9i18;
    // dual_t t9i38;
    // dual_t t9i58;
    // dual_t t916;
    // dual_t t976;
    // dual_t t9i76;
    dual_t lnt9;
};

template <typename dual_t>
AMREX_GPU_HOST_DEVICE inline
tf_t<dual_t> get_tfactors(dual_t temp)
{
    tf_t<dual_t> tf;

    tf.temp = temp;

    tf.t9    = temp * 1.0e-9_rt;
    tf.t92   = tf.t9*tf.t9;
    tf.t93   = tf.t9*tf.t92;
    // tf.t94   = tf.t9*tf.t93;
    // tf.t95   = tf.t9*tf.t94;
    tf.t95   = tf.t92*tf.t93;
    // tf.t96   = tf.t9*tf.t95;

    tf.t912  = admath::sqrt(tf.t9);
    tf.t932  = tf.t9*tf.t912;
    tf.t952  = tf.t9*tf.t932;
    // tf.t972  = tf.t9*tf.t952;
    tf.t972  = tf.t92*tf.t932;

    tf.t913  = admath::cbrt(tf.t9);
    tf.t923  = tf.t913*tf.t913;
    tf.t943  = tf.t9*tf.t913;
    tf.t953  = tf.t9*tf.t923;
    // tf.t973  = tf.t953*tf.t923;
    // tf.t9113 = tf.t973*tf.t943;

    // tf.t914  = admath::pow(tf.t9, 0.25e0_rt);
    // tf.t934  = tf.t914*tf.t914*tf.t914;
    // tf.t954  = tf.t9*tf.t914;
    // tf.t974  = tf.t9*tf.t934;

    // tf.t915  = admath::pow(tf.t9, 0.2_rt);
    // tf.t935  = tf.t915*tf.t915*tf.t915;
    // tf.t945  = tf.t915 * tf.t935;
    // tf.t965  = tf.t9 * tf.t915;

    // tf.t916  = admath::pow(tf.t9, 1.0_rt/6.0_rt);
    // tf.t976  = tf.t9 * tf.t916;
    // tf.t9i76 = 1.0e0_rt/tf.t976;

    // tf.t917  = admath::pow(tf.t9, 1.0_rt/7.0_rt);
    // tf.t927  = tf.t917*tf.t917;
    // tf.t947  = tf.t927*tf.t927;

    // tf.t918  = admath::sqrt(tf.t914);
    // tf.t938  = tf.t918*tf.t918*tf.t918;
    // tf.t958  = tf.t938*tf.t918*tf.t918;

    tf.t9i   = 1.0e0_rt/tf.t9;
    tf.t9i2  = tf.t9i*tf.t9i;
    // tf.t9i3  = tf.t9i2*tf.t9i;

    tf.t9i12 = 1.0e0_rt/tf.t912;
    tf.t9i32 = tf.t9i*tf.t9i12;
    // tf.t9i52 = tf.t9i*tf.t9i32;
    // tf.t9i72 = tf.t9i*tf.t9i52;

    tf.t9i13 = 1.0e0_rt/tf.t913;
    tf.t9i23 = tf.t9i13*tf.t9i13;
    tf.t9i43 = tf.t9i*tf.t9i13;
    tf.t9i53 = tf.t9i*tf.t9i23;

    // tf.t9i14 = 1.0e0_rt/tf.t914;
    // tf.t9i34 = tf.t9i14*tf.t9i14*tf.t9i14;
    // tf.t9i54 = tf.t9i*tf.t9i14;

    // tf.t9i15 = 1.0e0_rt/tf.t915;
    // tf.t9i35 = tf.t9i15*tf.t9i15*tf.t9i15;
    // tf.t9i45 = tf.t9i15 * tf.t9i35;
    // tf.t9i65 = tf.t9i*tf.t9i15;

    // tf.t9i17 = 1.0e0_rt/tf.t917;
    // tf.t9i27 = tf.t9i17*tf.t9i17;
    // tf.t9i47 = tf.t9i27*tf.t9i27;

    // tf.t9i18 = 1.0e0_rt/tf.t918;
    // tf.t9i38 = tf.t9i18*tf.t9i18*tf.t9i18;
    // tf.t9i58 = tf.t9i38*tf.t9i18*tf.t9i18;

    tf.lnt9 = admath::log(tf.t9);

    return tf;
}

#endif
