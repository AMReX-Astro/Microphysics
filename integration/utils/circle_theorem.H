#ifndef CIRCLE_THEOREM_H
#define CIRCLE_THEOREM_H

#include <ArrayUtilities.H>
#include <integrator_rhs_strang.H>
#include <limits>

template<typename BurnT, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void circle_theorem_sprad(const Real time, BurnT& state, T& int_state, Real& sprad)
{

    ArrayUtil::MathArray2D<1, INT_NEQS, 1, INT_NEQS> jac_array;

    jac(time, state, int_state, jac_array);

    // the Gershgorin circle theorem says that the spectral radius is <
    // max_i ( -a_{ii} + sum_{j,j!=i} |a_{ij}|)

    sprad = std::numeric_limits<Real>::lowest();

    for (int irow = 1; irow <= INT_NEQS; ++irow) {
        Real rho = -jac_array(irow, irow);
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
