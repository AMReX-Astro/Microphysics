#ifndef _vode_rhs_H_
#define _vode_rhs_H_

#include <AMReX_REAL.H>
#include <cuvode_parameters.H>
#include <vode_rpar.H>

// The f_rhs routine provides the right-hand-side for the DVODE solver.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void f_rhs (Real time, dvode_t& dvode_state, Real* ydot)
{

    ydot[0] = -.04e0_rt * dvode_state.y[0] + 1.e4_rt * dvode_state.y[1] * dvode_state.y[2];
    ydot[2] = 3.e7_rt * dvode_state.y[1] * dvode_state.y[1];
    ydot[1] = -ydot[0] - ydot[2];

}

// Analytical Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void jac(Real time, dvode_t& dvode_state, Real pd[][VODE_NEQS])
{

    pd[0][0] = -.04e0_rt;
    pd[1][0] = 1.e4_rt * dvode_state.y[2];
    pd[2][0] = 1.e4_rt * dvode_state.y[1];
    pd[0][1] = .04e0_rt;
    pd[2][1] = -pd[2][0];
    pd[1][2] = 6.e7_rt * dvode_state.y[1];
    pd[1][1] = -pd[1][0] - pd[1][2];

}

#endif
