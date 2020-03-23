#ifndef _react_zones_H_
#define _react_zones_H_

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void init_state(int i, int j, int k, int nx, int ny, int nz, Array4<Real>& state)
{
    int n = i + j * nx + k * nx * ny;
    int npts = nx * ny * nz;
    state(i,j,k,1) = 1.0_rt - ((Real) n) / ((Real) npts);
    state(i,j,k,2) = 0.0_rt;
    state(i,j,k,3) = 0.0_rt;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void do_react(int i, int j, int k, int ncomp, Array4<Real>& state, Real dt)
{

    use cuvode_parameters_module, only: MF_ANALYTIC_JAC, MF_NUMERICAL_JAC
    use cuvode_types_module, only: dvode_t
    use cuvode_module, only: dvode

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), ncomp)
    integer, intent(in), value :: ncomp
    real(rt), intent(in), value :: dt

    integer         :: ii, jj, kk, n

    ! VODE variables
    type (dvode_t) :: dvode_state
   

    !$gpu

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
    dvode_state.T = ZERO;
    dvode_state.TOUT = dt;

    // Initialize the initial conditions
    for (int n = 0; n < ncomp; ++n) {
        dvode_state.y[n] = state(i, j, k, n);
    }

    // Call the integration routine.
    dvode(dvode_state);

    // Check if the integration failed
    if (dvode_state.istate < 0) {
#ifndef AMREX_USE_CUDA       
        std::cout << 'ERROR: integration failed' << std::endl;
        std::cout << 'istate = ' << dvode_state.istate << std::endl;
        std::cout << 'time = ' << dvode_state.T << std::endl;
        amrex::Abort("integration failed");
#endif
    }
             
    // Store the final result
    for (int n = 0; n < ncomp ++n) {
        state(i, j, k, n) = dvode_state.y[n];
    }

}

#endif
