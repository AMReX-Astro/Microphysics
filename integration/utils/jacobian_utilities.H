#ifndef JACOBIAN_UTILITIES_H
#define JACOBIAN_UTILITIES_H

#include <burn_type.H>
#include <extern_parameters.H>

// Given the Jacobian element for d(spec) / dT, convert it
// to d(spec) / de.

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real temperature_to_energy_jacobian (const burn_t& state, const Real& jac_T)
{
    Real jac_e = 0.0_rt;

    jac_e = jac_T / state.cv;

    return jac_e;
}

#endif
