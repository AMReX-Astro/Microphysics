#ifndef _react_zones_H_
#define _react_zones_H_

#include <cuvode.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void do_react(int i, int j, int k, int ncomp, const Array4<Real>& state, Real dt)
{

    dvode_t dvode_state;

    // Use an analytic Jacobian
    dvode_state.jacobian = 1;

    // Set the absolute tolerances
    dvode_state.atol[0] = 1.e-8_rt;
    dvode_state.atol[1] = 1.e-14_rt;
    dvode_state.atol[2] = 1.e-6_rt;

    // Set the relative tolerances
    dvode_state.rtol[0] = 1.e-4_rt;
    dvode_state.rtol[1] = 1.e-4_rt;
    dvode_state.rtol[2] = 1.e-4_rt;

    // We want VODE to re-initialize each time we call it.
    dvode_state.istate = 1;

    // Take no more than 500 steps.
    dvode_state.MXSTEP = 500;

    // Initialize the integration time and set the final time to dt
    dvode_state.t = 0.0_rt;
    dvode_state.tout = dt;

    // Initialize the initial conditions
    for (int n = 0; n < ncomp; ++n) {
        dvode_state.y[n] = state(i, j, k, n);
    }

    // Call the integration routine.
    dvode(dvode_state);

    // Check if the integration failed
    if (dvode_state.istate < 0) {
#ifndef AMREX_USE_CUDA       
        std::cout << "ERROR: integration failed" << std::endl;
        std::cout << "istate = " << dvode_state.istate << std::endl;
        std::cout << "time = " << dvode_state.t << std::endl;
        amrex::Abort("integration failed");
#endif
    }
             
    // Store the final result
    for (int n = 0; n < ncomp; ++n) {
        state(i, j, k, n) = dvode_state.y[n];
    }

}

#endif
