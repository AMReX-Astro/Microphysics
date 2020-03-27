#ifndef _vode_rhs_H_
#define _vode_rhs_H_

#include <AMReX_REAL.H>
#include <cuvode_parameters.H>
#include <vode_rpar.H>

// The f_rhs routine provides the right-hand-side for the DVODE solver.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void f_rhs (Real time, dvode_t& dvode_state,
            Array1D<Real, 1, VODE_NEQS>& ydot)
{

    ydot(1) = -.04e0_rt * dvode_state.y(1) + 1.e4_rt * dvode_state.y(2) * dvode_state.y(3);
    ydot(3) = 3.e7_rt * dvode_state.y(2) * dvode_state.y(2);
    ydot(2) = -ydot(1) - ydot(3);

}

// Analytical Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void jac(Real time, dvode_t& dvode_state,
         Array2D<Real, 1, VODE_NEQS, 1, VODE_NEQS>& pd)
{

    pd(1,1) = -.04e0_rt;
    pd(1,2) = 1.e4_rt * dvode_state.y(3);
    pd(1,3) = 1.e4_rt * dvode_state.y(2);
    pd(2,1) = .04e0_rt;
    pd(2,3) = -pd(1,3);
    pd(3,1) = 0.0_rt;
    pd(3,2) = 6.e7_rt * dvode_state.y(2);
    pd(2,2) = -pd(1,2) - pd(3,2);
    pd(3,3) = 0.0_rt;

}

#endif
