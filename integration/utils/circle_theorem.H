#ifndef CIRCLE_THEOREM_H
#define CIRCLE_THEOREM_H

#include <ArrayUtilities.H>
#ifdef STRANG
#include <integrator_type_strang.H>
#include <integrator_rhs_strang.H>
#endif
#ifdef SDC
#include <integrator_type_sdc.H>
#include <integrator_rhs_sdc.H>
#endif
#include <limits>
#include <numerical_jacobian.H>

template<typename BurnT, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void circle_theorem_sprad(const amrex::Real time, BurnT& state, T& int_state, amrex::Real& sprad)
{

    ArrayUtil::MathArray2D<1, INT_NEQS, 1, INT_NEQS> jac_array;

    if (integrator_rp::jacobian == 1) {
        jac(time, state, int_state, jac_array);
    } else {
        integrator_to_burn(int_state, state);
        jac_info_t jac_info;
        jac_info.h = 0.0_rt;
        numerical_jac(state, jac_info, jac_array);
    }

    // the Gershgorin circle theorem says that the spectral radius is <
    // max_i ( -a_{ii} + sum_{j,j!=i} |a_{ij}|)

    sprad = std::numeric_limits<amrex::Real>::lowest();

    for (int irow = 1; irow <= INT_NEQS; ++irow) {
        amrex::Real rho = -jac_array(irow, irow);
        for (int jcol = 1; jcol <= INT_NEQS; ++jcol) {
            if (jcol == irow) {
                continue;
            }
            rho += std::abs(jac_array(irow, jcol));
        }
        sprad = std::max(sprad, rho);
    }

}

#endif
