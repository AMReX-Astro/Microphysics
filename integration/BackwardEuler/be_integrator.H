#ifndef BE_INTEGRATOR_H
#define BE_INTEGRATOR_H

#include <AMReX_Algorithm.H>

#include <be_type.H>
#include <network.H>
#include <actual_network.H>
#ifndef NEW_NETWORK_IMPLEMENTATION
#include <actual_rhs.H>
#endif
#include <burn_type.H>
#include <linpack.H>
#include <numerical_jacobian.H>
#ifdef STRANG
#include <integrator_rhs_strang.H>
#endif
#ifdef SDC
#include <integrator_rhs_sdc.H>
#endif
#include <integrator_data.H>
#include <initial_timestep.H>

///
/// update state.xn[] and state.e through a timestep dt
/// state is updated in place -- if we are unsuccessful, we reset it
/// we return an error code
///
template <typename BurnT, typename BeT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int single_step (BurnT& state, BeT& be, const amrex::Real dt)
{
    constexpr int int_neqs = integrator_neqs<BurnT>();

    int ierr = IERR_SUCCESS;
    bool converged = false;

    // create our current guess for the solution -- just as a first
    // order explicit prediction

    amrex::Array1D<amrex::Real, 1, int_neqs> ydot;

    rhs(be.t, state, be, ydot);

    be.n_rhs += 1;

    amrex::Array1D<amrex::Real, 1, int_neqs> y_old;

    for (int n = 1; n <= int_neqs; n++) {
        y_old(n) = be.y(n);
        be.y(n) = be.y(n) + dt * ydot(n);
    }

    // Newton loop

    for (int iter = 1; iter <= integrator_rp::max_iter; iter++) {

        // work with the current guess

        // get the ydots for our current guess of y

        rhs(be.t, state, be, ydot);
        be.n_rhs += 1;

        // construct the Jacobian

        if (be.jacobian_type == 1) {
            jac(be.t, state, be, be.jac);
        } else {
            jac_info_t jac_info;
            jac_info.h = dt;
            numerical_jac(state, jac_info, be.jac);
            be.n_rhs += (NumSpec+1);
        }

        be.n_jac++;

        // construct the matrix for the linear system
        // (I - dt J) dy^{n+1} = rhs

        for (int n = 1; n <= int_neqs; n++) {
            for (int m  = 1; m <= int_neqs; m++) {
                be.jac(m, n) *= -dt;
                if (m == n) {
                    be.jac(m, n) = 1.0_rt + be.jac(m, n);
                }
            }
        }

        // construct the RHS of our linear system

        amrex::Array1D<amrex::Real, 1, int_neqs> b;
        for (int n = 1; n <= int_neqs; n++) {
            b(n) = y_old(n) - be.y(n) + dt * ydot(n);
        }

        // solve the linear system

        int ierr_linpack;
        IArray1D pivot;

        if (integrator_rp::linalg_do_pivoting == 1) {
            constexpr bool allow_pivot{true};
            dgefa<int_neqs, allow_pivot>(be.jac, pivot, ierr_linpack);
        } else {
            constexpr bool allow_pivot{false};
            dgefa<int_neqs, allow_pivot>(be.jac, pivot, ierr_linpack);
        }

        if (ierr_linpack != 0) {
            ierr = IERR_LU_DECOMPOSITION_ERROR;
            break;
        }

        if (integrator_rp::linalg_do_pivoting == 1) {
            constexpr bool allow_pivot{true};
            dgesl<int_neqs, allow_pivot>(be.jac, pivot, b);
        } else {
            constexpr bool allow_pivot{false};
            dgesl<int_neqs, allow_pivot>(be.jac, pivot, b);
        }

        // update our current guess for the solution

        for (int n = 1; n <= int_neqs; n++) {
            be.y(n) += b(n);
        }

        // check to see if we converged
        // we compute the norms

        amrex::Real y_norm = 0.0;
        amrex::Real b_norm = 0.0;
        for (int n = 1; n <= int_neqs; n++) {
            y_norm += be.y(n) * be.y(n);
            b_norm += b(n) * b(n);
        }
        y_norm = std::sqrt(y_norm / int_neqs);
        b_norm = std::sqrt(b_norm / int_neqs);

        if (b_norm < integrator_rp::tol * y_norm) {
            converged = true;
            break;
        }

    }

    // we are done iterating -- did we converge?

    if (! converged) {

        if (ierr == IERR_SUCCESS) {

            // if we didn't set another error, then we probably ran
            // out of iterations, so set nonconvergence

            ierr = IERR_CORRECTOR_CONVERGENCE;

            // reset the solution to the original
            for (int n = 1; n <= int_neqs; n++) {
                be.y(n) = y_old(n);
            }

        }

    }

    return ierr;

}

template <typename BurnT, typename BeT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int be_integrator (BurnT& state, BeT& be)
{
    constexpr int int_neqs = integrator_neqs<BurnT>();

    be.n_rhs = 0;
    be.n_jac = 0;
    be.n_step = 0;

    int ierr;

    if (integrator_rp::do_single_step == 1) {
        // Do a single step to the final time
        const amrex::Real dt_single_step = be.tout - be.t;
        ierr = single_step(state, be, dt_single_step);
        be.t = be.tout;
        ++be.n_step;
        return ierr;
    }

    // estimate the timestep

    amrex::Array1D<amrex::Real, 1, int_neqs> ydot;
    rhs(be.t, state, be, ydot);

    be.n_rhs += 1;

    amrex::Real dt_sub = initial_react_dt(state, be, ydot);

    // main timestepping loop

    while (be.t < (1.0_rt - timestep_safety_factor) * be.tout &&
           be.n_step < integrator_rp::ode_max_steps) {

        // store the current solution -- we'll revert to this if a step fails
        amrex::Array1D<amrex::Real, 1, int_neqs> y_old;
        for (int n = 1; n <= int_neqs; ++n) {
            y_old(n) = be.y(n);
        }

        // don't go too far

        if (be.t + dt_sub > be.tout) {
            dt_sub = be.tout - be.t;
        }

        // our strategy is to take 2 steps at dt/2 and one at dt and
        // to compute the error from those


        // try to take a step dt

        // first do 2 (fine) dt/2 steps

        amrex::Array1D<amrex::Real, 1, int_neqs> y_fine;

        ierr = single_step(state, be, dt_sub/2);
        if (ierr == IERR_SUCCESS) {
            ierr = single_step(state, be, dt_sub/2);

            // store the fine dt solution

            for (int n = 1; n <= int_neqs; ++n) {
                y_fine(n) = be.y(n);
            }

            // now that single (coarse) dt step
            // first reset the solution
            for (int n = 1; n <= int_neqs; ++n) {
                be.y(n) = y_old(n);
            }
            ierr = single_step(state, be, dt_sub);
        }

        // define a weight for each variable to use in checking the error

        amrex::Array1D<amrex::Real, 1, int_neqs> w;
        for (int n = 1; n <= NumSpec; n++) {
            w(n) = 1.0_rt / (be.rtol_spec * std::abs(y_fine(n)) + be.atol_spec);
        }
        w(net_ienuc) = 1.0_rt / (be.rtol_enuc * std::abs(y_fine(net_ienuc)) + be.atol_enuc);

        // now look for w |y_fine - y_coarse| < 1

        amrex::Real rel_error = 0.0_rt;
        for (int n = 1; n <= NumSpec; n++) {
            rel_error = amrex::max(rel_error, w(n) * std::abs(y_fine(n) - be.y(n)));
        }
        rel_error = amrex::max(rel_error, w(net_ienuc) * std::abs(y_fine(net_ienuc) - be.y(net_ienuc)));

        bool step_success = false;
        if (rel_error < 1.0_rt) {
            step_success = true;
        }

        if (ierr == IERR_SUCCESS && step_success) {

            // there were no problems with the mechanics of the step
            // (linear algebra, etc.) and we met our error
            // goals.

            // y_fine has the current best solution

            be.t += dt_sub;

            for (int n = 1; n <= int_neqs; ++n) {
                be.y(n) = y_fine(n);
            }

            // can we potentially increase the timestep?
            // backward-Euler has a local truncation error of dt**2

            amrex::Real dt_new = dt_sub * std::sqrt(1.0_rt / rel_error);
            dt_sub = amrex::Clamp(dt_new, dt_sub / 2.0, 2.0 * dt_sub);

        } else {

            // roll back the solution
            for (int n = 1; n <= int_neqs; ++n) {
                be.y(n) = y_old(n);
            }

            // adjust the timestep and try again
            dt_sub /= 2;

        }

        ++be.n_step;
        ++be.n_step;

    }

    if (be.n_step >= integrator_rp::ode_max_steps) {
        ierr = IERR_TOO_MANY_STEPS;
    }

    return ierr;

}

#endif
