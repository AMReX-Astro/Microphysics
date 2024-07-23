#ifndef rprox_rates_H
#define rprox_rates_H

#include <AMReX_REAL.H>
#include <tfactors.H>
#include <microphysics_autodiff.H>

template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_c12_to_n13 (const tf_t<number_t>& tfactors, number_t& rate)
{
    rate = 0.0_rt;

    number_t ct9i = tfactors.t9i * (-3.77849e0_rt);
    number_t ct9i13 = tfactors.t9i13 * (-5.10735e0_rt);
    number_t ct913 = tfactors.t913 * (-2.24111e0_rt);
    number_t ct9 = tfactors.t9 * (0.148883e0_rt);
    number_t clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r0 = admath::exp(17.5428e0_rt +
                       ct9i +
                       ct9i13 +
                       ct913 +
                       ct9 +
                       clnt9);

    ct9i13 = tfactors.t9i13 * (-13.692e0_rt);
    ct913 = tfactors.t913 * (-0.230881e0_rt);
    ct9 = tfactors.t9 * (4.44362e0_rt);
    number_t ct953 = tfactors.t953 * (-3.15898e0_rt);
    clnt9 = tfactors.lnt9 * (-0.666667e0_rt);

    number_t r1 = admath::exp(17.1482e0_rt +
                       ct9i13 +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    rate = r0 + r1 ;
}


template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_f17_to_p_o16 (const tf_t<number_t>& tfactors, number_t& rate)
{
    rate = 0.0_rt;

    number_t ct9i = tfactors.t9i * (-6.96583e0_rt);
    number_t ct9i13 = tfactors.t9i13 * (-16.696e0_rt);
    number_t ct913 = tfactors.t913 * (-1.16252e0_rt);
    number_t ct9 = tfactors.t9 * (0.267703e0_rt);
    number_t ct953 = tfactors.t953 * (-0.0338411e0_rt);
    number_t clnt9 = tfactors.lnt9 * (0.833333e0_rt);

    number_t r0 = admath::exp(40.9135e0_rt +
                       ct9i +
                       ct9i13 +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    rate = r0;
}


template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_f17_to_o17 ([[maybe_unused]] const tf_t<number_t>& tfactors, number_t& rate)
{
    rate = std::exp(-4.53302e0_rt);
}


template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_f17_to_ne18 (const tf_t<number_t>& tfactors, number_t& rate)
{
    rate = 0.0_rt;

    number_t ct9i = tfactors.t9i * (-0.0323504e0_rt);
    number_t ct9i13 = tfactors.t9i13 * (-14.2191e0_rt);
    number_t ct913 = tfactors.t913 * (34.0647e0_rt);
    number_t ct9 = tfactors.t9 * (-16.5698e0_rt);
    number_t ct953 = tfactors.t953 * (2.48116e0_rt);
    number_t clnt9 = tfactors.lnt9 * (-2.13376e0_rt);

    number_t r0 = admath::exp(-7.84708e0_rt +
                       ct9i +
                       ct9i13 +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    ct9i = tfactors.t9i * (-4.95969e0_rt);
    ct9i13 = tfactors.t9i13 * (-21.3249e0_rt);
    ct913 = tfactors.t913 * (-0.230774e0_rt);
    ct9 = tfactors.t9 * (0.917931e0_rt);
    ct953 = tfactors.t953 * (-0.0440377e0_rt);
    clnt9 = tfactors.lnt9 * (-7.36014e0_rt);

    number_t r1 = admath::exp(27.5778e0_rt +
                       ct9i +
                       ct9i13 +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    rate = r0 + r1;
}


template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_he4_he4_to_c12 (const tf_t<number_t>& tfactors, number_t& rate)
{
    rate = 0.0_rt;

    number_t ct9i = tfactors.t9i * (-1.02446e0_rt);
    number_t ct9i13 = tfactors.t9i13 * (-23.57e0_rt);
    number_t ct913 = tfactors.t913 * (20.4886e0_rt);
    number_t ct9 = tfactors.t9 * (-12.9882e0_rt);
    number_t ct953 = tfactors.t953 * (-20.0e0_rt);
    number_t clnt9 = tfactors.lnt9 * (-2.16667e0_rt);

    number_t r0 = admath::exp(-11.7884e0_rt +
                       ct9i +
                       ct9i13 +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    ct9i13 = tfactors.t9i13 * (-37.06e0_rt);
    ct913 = tfactors.t913 * (29.3493e0_rt);
    ct9 = tfactors.t9 * (-115.507e0_rt);
    ct953 = tfactors.t953 * (-10.0e0_rt);
    clnt9 = tfactors.lnt9 * (-1.33333e0_rt);

    number_t r1 = admath::exp(-0.971052e0_rt +
                       ct9i13 +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    ct9i = tfactors.t9i * (-4.12656e0_rt);
    ct9i13 = tfactors.t9i13 * (-13.49e0_rt);
    ct913 = tfactors.t913 * (21.4259e0_rt);
    ct9 = tfactors.t9 * (-1.34769e0_rt);
    ct953 = tfactors.t953 * (0.0879816e0_rt);
    clnt9 = tfactors.lnt9 * (-13.1653e0_rt);

    number_t r2 = admath::exp(-24.3505e0_rt +
                       ct9i +
                       ct9i13 +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    rate = r0 + r1 + r2;
}


template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_n14_to_o15 (const tf_t<number_t>& tfactors, number_t& rate)
{
    rate = 0.0_rt;

    number_t ct9i13 = tfactors.t9i13 * (-15.193e0_rt);
    number_t ct913 = tfactors.t913 * (-4.63975e0_rt);
    number_t ct9 = tfactors.t9 * (9.73458e0_rt);
    number_t ct953 = tfactors.t953 * (-9.55051e0_rt);
    number_t clnt9 = tfactors.lnt9 * (0.333333e0_rt);

    number_t r0 = admath::exp(20.1169e0_rt +
                       ct9i13 +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    ct9i13 = tfactors.t9i13 * (-15.193e0_rt);
    ct913 = tfactors.t913 * (-0.161954e0_rt);
    ct9 = tfactors.t9 * (-7.52123e0_rt);
    ct953 = tfactors.t953 * (-0.987565e0_rt);
    clnt9 = tfactors.lnt9 * (-0.666667e0_rt);

    number_t r1 = admath::exp(17.01e0_rt +
                       ct9i13 +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    number_t ct9i = tfactors.t9i * (-4.891e0_rt);
    clnt9 = tfactors.lnt9 * (0.0682e0_rt);

    number_t r2 = admath::exp(6.73578e0_rt +
                       ct9i +
                       clnt9);

    ct9i = tfactors.t9i * (-2.998e0_rt);
    clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r3 = admath::exp(7.65444e0_rt +
                       ct9i +
                       clnt9);

    rate = r0 + r1 + r2 + r3;
}


template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_ne18_to_p_na21 (const tf_t<number_t>& tfactors, number_t& rate)
{
    rate = 0.0_rt;

    number_t ct9i = tfactors.t9i * (-24.7176e0_rt);
    number_t clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r0 = admath::exp(19.4058e0_rt +
                       ct9i +
                       clnt9);

    ct9i = tfactors.t9i * (-0.452576e0_rt);
    clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r1 = admath::exp(-137.358e0_rt +
                       ct9i +
                       clnt9);

    ct9i = tfactors.t9i * (-10.885e0_rt);
    clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r2 = admath::exp(6.39797e0_rt +
                       ct9i +
                       clnt9);

    ct9i = tfactors.t9i * (-7.45009e0_rt);
    clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r3 = admath::exp(-1.15641e0_rt +
                       ct9i +
                       +clnt9);

    ct9i = tfactors.t9i * (-5.97632e0_rt);
    clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r4 = admath::exp(-6.65137e0_rt +
                       ct9i +
                       clnt9);

    ct9i = tfactors.t9i * (-3.14078e0_rt);
    number_t ct9i13 = tfactors.t9i13 * (-17.7759e0_rt);
    number_t ct913 = tfactors.t913 * (36.0724e0_rt);
    number_t ct9 = tfactors.t9 * (-5.34039e0_rt);
    number_t ct953 = tfactors.t953 * (0.382679e0_rt);
    clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r5 = admath::exp(-13.3392e0_rt +
                       ct9i +
                       ct9i13 +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    ct9i = tfactors.t9i * (-2.81989e0_rt);
    clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r6 = admath::exp(-28.6929e0_rt +
                       ct9i +
                       clnt9);

    rate = r0 + r1 + r2 + r3 + r4 + r5 + r6;
}


template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ne18_to_f18 ([[maybe_unused]] const tf_t<number_t>& tfactors, number_t& rate)
{
    rate = std::exp(-0.880534e0_rt);
}


template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ne19_to_f19 ([[maybe_unused]] const tf_t<number_t>& tfactors, number_t& rate)
{
    rate = std::exp(-3.21258e0_rt);
}


template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_ne19_to_na20 (const tf_t<number_t>& tfactors, number_t& rate)
{
    rate = 0.0_rt;

    number_t ct9i = tfactors.t9i * (-5.07623e0_rt);
    number_t ct913 = tfactors.t913 * (1.23704e0_rt);
    number_t ct9 = tfactors.t9 * (0.337618e0_rt);
    number_t ct953 = tfactors.t953 * (-0.0562825e0_rt);
    number_t clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r0 = admath::exp(5.63289e0_rt +
                       ct9i +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    number_t ct9i13 = tfactors.t9i13 * (-19.5908e0_rt);
    ct913 = tfactors.t913 * (-2.37696e0_rt);
    ct9 = tfactors.t9 * (3.26815e0_rt);
    ct953 = tfactors.t953 * (-1.06524e0_rt);
    clnt9 = tfactors.lnt9 * (-0.666667e0_rt);

    number_t r1 = admath::exp(17.822e0_rt +
                       ct9i13 +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    rate = r0 + r1;
}


template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_o14_to_p_f17 (const tf_t<number_t>& tfactors, number_t& rate)
{
    rate = 0.0_rt;

    number_t ct9i = tfactors.t9i * (-12.0223e0_rt);
    number_t clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r0 = admath::exp(12.1289e0_rt +
                       ct9i +
                       clnt9);

    ct9i = tfactors.t9i * (-26.0e0_rt);
    clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r1 = admath::exp(18.6518e0_rt +
                       ct9i +
                       clnt9);

    number_t ct9i13 = tfactors.t9i13 * (-39.388e0_rt);
    number_t ct913 = tfactors.t913 * (-17.4673e0_rt);
    number_t ct9 = tfactors.t9 * (35.3029e0_rt);
    number_t ct953 = tfactors.t953 * (-24.8162e0_rt);
    clnt9 = tfactors.lnt9 * (-0.666667e0_rt);

    number_t r2 = admath::exp(40.8358e0_rt +
                       ct9i13 +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    ct9i = tfactors.t9i * (-22.51e0_rt);
    clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r3 = admath::exp(16.3087e0_rt +
                       ct9i +
                       clnt9);

    ct9i = tfactors.t9i * (-13.6e0_rt);
    clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r4 = admath::exp(11.1184e0_rt +
                       ct9i +
                       clnt9);

    ct9i = tfactors.t9i * (-0.453036e0_rt);
    clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r5 = admath::exp(-106.091e0_rt +
                       ct9i +
                       clnt9);

    rate = r0 + r1 + r2 + r3 + r4 + r5;
}


template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_o14_to_n14 ([[maybe_unused]] const tf_t<number_t>& tfactors, number_t& rate)
{
    rate = std::exp(-4.62412e0_rt);
}


template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_o15_to_ne19 (const tf_t<number_t>& tfactors, number_t& rate)
{
    rate = 0.0_rt;

    number_t ct9i = tfactors.t9i * (-5.88439e0_rt);
    number_t clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r0 = admath::exp(-0.0452465e0_rt +
                       ct9i +
                       clnt9);

    number_t ct9i13 = tfactors.t9i13 * (-39.578e0_rt);
    number_t ct953 = tfactors.t953 * (-3.0e0_rt);
    clnt9 = tfactors.lnt9 * (-0.666667e0_rt);

    number_t r1 = admath::exp(26.2914e0_rt +
                       ct9i13 +
                       ct953 +
                       clnt9);

    ct9i = tfactors.t9i * (-4.20439e0_rt);
    ct9i13 = tfactors.t9i13 * (-3.24609e0_rt);
    number_t ct913 = tfactors.t913 * (44.4647e0_rt);
    number_t ct9 = tfactors.t9 * (-9.79962e0_rt);
    ct953 = tfactors.t953 * (0.841782e0_rt);
    clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r2 = admath::exp(-32.2496e0_rt +
                       ct9i +
                       ct9i13 +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    rate = r0 + r1 + r2;
}


template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_o15_to_n15 ([[maybe_unused]] const tf_t<number_t>& tfactors, number_t& rate)
{
    rate = std::exp(-5.1725e0_rt);
}


template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_o16_to_ne20 (const tf_t<number_t>& tfactors, number_t& rate)
{
    rate = 0.0_rt;

    number_t ct9i = tfactors.t9i * (-10.3585e0_rt);
    number_t clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r0 = admath::exp(3.88571e0_rt +
                       ct9i +
                       clnt9);

    number_t ct9i13 = tfactors.t9i13 * (-39.7262e0_rt);
    number_t ct913 = tfactors.t913 * (-0.210799e0_rt);
    number_t ct9 = tfactors.t9 * (0.442879e0_rt);
    number_t ct953 = tfactors.t953 * (-0.0797753e0_rt);
    clnt9 = tfactors.lnt9 * (-0.666667e0_rt);

    number_t r1 = admath::exp(23.903e0_rt +
                       ct9i13 +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    ct9i = tfactors.t9i * (-12.7643e0_rt);
    ct913 = tfactors.t913 * (-3.65925e0_rt);
    ct9 = tfactors.t9 * (0.714224e0_rt);
    ct953 = tfactors.t953 * (-0.00107508e0_rt);
    clnt9 = tfactors.lnt9 * (-1.5e0_rt);

    number_t r2 = admath::exp(9.50848e0_rt +
                       ct9i +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    rate = r0 + r1 + r2;
}


template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_o16_to_f17 (const tf_t<number_t>& tfactors, number_t& rate)
{
    rate = 0.0_rt;

    number_t ct9i13 = tfactors.t9i13 * (-16.696e0_rt);
    number_t ct913 = tfactors.t913 * (-1.16252e0_rt);
    number_t ct9 = tfactors.t9 * (0.267703e0_rt);
    number_t ct953 = tfactors.t953 * (-0.0338411e0_rt);
    number_t clnt9 = tfactors.lnt9 * (-0.666667e0_rt);

    number_t r0 = admath::exp(19.0904e0_rt +
                       ct9i13 +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    rate = r0;
}


template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_si26_to_p_p29 (const tf_t<number_t>& tfactors, number_t& rate)
{
    rate = 0.0_rt;

    number_t ct9i13 = tfactors.t9i13 * (-59.3013e0_rt);
    number_t ct913 = tfactors.t913 * (0.480742e0_rt);
    number_t ct9 = tfactors.t9 * (-0.834505e0_rt);
    number_t ct953 = tfactors.t953 * (0.0621841e0_rt);
    number_t clnt9 = tfactors.lnt9 * (-0.666667e0_rt);

    number_t r0 = admath::exp(48.8732e0_rt +
                       ct9i13 +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    rate = r0;
}


template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_ti44_to_p_v47 (const tf_t<number_t>& tfactors, number_t& rate)
{
    rate = 0.0_rt;

    number_t ct9i = tfactors.t9i * (-9.07869e0_rt);
    number_t ct9i13 = tfactors.t9i13 * (5.56533e0_rt);
    number_t ct913 = tfactors.t913 * (18.4415e0_rt);
    number_t ct9 = tfactors.t9 * (-4.10095e0_rt);
    number_t ct953 = tfactors.t953 * (0.24244e0_rt);
    number_t clnt9 = tfactors.lnt9 * (16.0516e0_rt);

    number_t r0 = admath::exp(-34.2468e0_rt +
                       ct9i +
                       ct9i13 +
                       ct913 +
                       ct9 +
                       ct953 +
                       clnt9);

    rate = r0;
}

#endif
