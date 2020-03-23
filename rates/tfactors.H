#ifndef _tfactors_H_
#define _tfactors_H_

#include <AMReX.H>

struct tf_t {
    amrex::Real temp;
    amrex::Real t9;
    amrex::Real t92;
    amrex::Real t93;
    // amrex::Real t94;
    amrex::Real t95;
    // amrex::Real t96;
    amrex::Real t912;
    amrex::Real t932;
    amrex::Real t952;
    amrex::Real t972;
    amrex::Real t913;
    amrex::Real t923;
    amrex::Real t943;
    amrex::Real t953;
    // amrex::Real t973;
    // amrex::Real t9113;
    // amrex::Real t914;
    // amrex::Real t934;
    // amrex::Real t954;
    // amrex::Real t974;
    // amrex::Real t915;
    // amrex::Real t935;
    // amrex::Real t945;
    // amrex::Real t965;
    // amrex::Real t917;
    // amrex::Real t927;
    // amrex::Real t947;
    // amrex::Real t918;
    // amrex::Real t938;
    // amrex::Real t958;
    amrex::Real t9i;
    amrex::Real t9i2;
    // amrex::Real t9i3;
    amrex::Real t9i12;
    amrex::Real t9i32;
    // amrex::Real t9i52;
    // amrex::Real t9i72;
    amrex::Real t9i13;
    amrex::Real t9i23;
    amrex::Real t9i43;
    amrex::Real t9i53;
    // amrex::Real t9i14;
    // amrex::Real t9i34;
    // amrex::Real t9i54;
    // amrex::Real t9i15;
    // amrex::Real t9i35;
    // amrex::Real t9i45;
    // amrex::Real t9i65;
    // amrex::Real t9i17;
    // amrex::Real t9i27;
    // amrex::Real t9i47;
    // amrex::Real t9i18;
    // amrex::Real t9i38;
    // amrex::Real t9i58;
    // amrex::Real t916;
    // amrex::Real t976;
    // amrex::Real t9i76;
};

AMREX_GPU_HOST_DEVICE inline 
tf_t get_tfactors(amrex::Real temp)
{
    tf_t tf; 

    tf.temp = temp;

    tf.t9    = temp * 1.0e-9_rt;
    tf.t92   = tf.t9*tf.t9;
    tf.t93   = tf.t9*tf.t92;
    // tf.t94   = tf.t9*tf.t93;
    // tf.t95   = tf.t9*tf.t94;
    tf.t95   = tf.t92*tf.t93;
    // tf.t96   = tf.t9*tf.t95;
    
    tf.t912  = std::sqrt(tf.t9);
    tf.t932  = tf.t9*tf.t912;
    tf.t952  = tf.t9*tf.t932;
    // tf.t972  = tf.t9*tf.t952;
    tf.t972  = tf.t92*tf.t932;
    
    tf.t913  = std::pow(tf.t9, 1.0_rt/3.0_rt);
    tf.t923  = tf.t913*tf.t913;
    tf.t943  = tf.t9*tf.t913;
    tf.t953  = tf.t9*tf.t923;
    // tf.t973  = tf.t953*tf.t923;
    // tf.t9113 = tf.t973*tf.t943;

    // tf.t914  = std::pow(tf.t9, 0.25e0_rt);
    // tf.t934  = tf.t914*tf.t914*tf.t914;
    // tf.t954  = tf.t9*tf.t914;
    // tf.t974  = tf.t9*tf.t934;
    
    // tf.t915  = std::pow(tf.t9, 0.2_rt);
    // tf.t935  = tf.t915*tf.t915*tf.t915;
    // tf.t945  = tf.t915 * tf.t935;
    // tf.t965  = tf.t9 * tf.t915;
    
    // tf.t916  = std::pow(tf.t9, 1.0_rt/6.0_rt);
    // tf.t976  = tf.t9 * tf.t916;
    // tf.t9i76 = 1.0e0_rt/tf.t976;
    
    // tf.t917  = std::pow(tf.t9, 1.0_rt/7.0_rt);
    // tf.t927  = tf.t917*tf.t917;
    // tf.t947  = tf.t927*tf.t927;
    
    // tf.t918  = std::sqrt(tf.t914);
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

    return tf;
}

#endif