#ifndef rprox_tfactors_H
#define rprox_tfactors_H

struct tf_t {
    Real T9;
    Real T9i;
    Real T9i13;
    Real T913;
    Real T953;
    Real lnT9;
};

AMREX_GPU_HOST_DEVICE AMREX_INLINE
tf_t calc_tfactors (const Real& T9)
{
    tf_t tfactors;

    tfactors.T9 = T9;
    tfactors.T9i = 1.e0_rt / tfactors.T9;
    tfactors.T9i13 = std::pow(tfactors.T9i, 1.0_rt / 3.0_rt);
    tfactors.T913 = std::pow(tfactors.T9, 1.0_rt / 3.0_rt);
    tfactors.T953 = tfactors.T9 * tfactors.T913 * tfactors.T913;
    tfactors.lnT9 = std::log(tfactors.T9);

    return tfactors;
}

#endif
