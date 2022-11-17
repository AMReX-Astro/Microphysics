#ifndef BE_INTEGRATOR_H
#define BE_INTEGRATOR_H

#include <be_type.H>
#include <network.H>
#include <actual_network.H>
#include <actual_rhs.H>
#include <burn_type.H>
#include <linpack.H>
#include <numerical_jacobian.H>
#include <be_type_strang.H>

///
/// update state.xn[] and state.e through a timestep dt
/// state is updated in place -- if we are unsuccessful, we reset it
/// we return an error code
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int single_step (burn_t& state, be_t& be, const Real dt)
{

    clean_state(state);

    // pack the initial state into a vector -- we need this for the
    // RHS of the linear system and to reset the state if we are
    // unsuccessful here

    for (int n = 1; n <= NumSpec; n++) {
        be.y(n) = state.xn[n-1];
    }
    be.y(net_ienuc) = state.e;

    int ierr = BE_SUCCESS;
    bool converged = false;

    // create our current guess for the solution -- just as a first
    // order explicit prediction

    Array1D<Real, 1, neqs> ydot;

    evaluate_rhs(state, ydot);

    be.n_rhs += 1;

    Array1D<Real, 1, BE_NEQS> y_np1;

    for (int n = 1; n <= BE_NEQS; n++) {
        y_np1(n) = be.y(n) + dt * ydot(n);
    }

    // Newton loop

    for (int iter = 1; iter <= max_iter; iter++) {

        // work with the current guess

        for (int n = 0; n < NumSpec; n++) {
            state.xn[n] = y_np1(n+1);
        }
        state.e = y_np1(net_ienuc);

        clean_state(state);

        // construct the Jacobian

        if (jacobian == 1) {
            actual_jac(state, be.jac);

            // We integrate X, not Y
            for (int j = 1; j <= NumSpec; ++j) {
                for (int i = 1; i <= neqs; ++i) {
                    be.jac(j, i) *= aion[j-1];
                    be.jac(i, j) *= aion_inv[j-1];
                }
            }

        } else {
            jac_info_t jac_info;
            jac_info.h = dt;
            numerical_jac(state, jac_info, be.jac);
        }

        if (!integrate_energy) {
            for (int j = 1; j <= BE_NEQS; j++) {
                be.jac(net_ienuc, j) = 0.0_rt;
            }
        }

        // construct the matrix for the linear system
        // (I - dt J) dy^{n+1} = rhs

        for (int m  = 1; m <= BE_NEQS; m++) {
            for (int n = 1; n <= BE_NEQS; n++) {
                be.jac(m, n) *= -dt;
                if (m == n) {
                    be.jac(m, n) = 1.0_rt + be.jac(m, n);
                }
            }
        }

        // get the ydots for our current guess of y

        evaluate_rhs(state, ydot);

        be.n_rhs += 1;

        // construct the RHS of our linear system

        Array1D<Real, 1, BE_NEQS> b;
        for (int n = 1; n <= BE_NEQS; n++) {
            b(n) = be.y(n) - y_np1(n) + dt * ydot(n);
        }

        // solve the linear system

        int ierr_linpack;
        IArray1D pivot;
        dgefa<BE_NEQS>(be.jac, pivot, ierr_linpack);

        if (ierr_linpack != 0) {
            ierr = BE_LU_DECOMPOSITION_ERROR;
            break;
        }

        dgesl<BE_NEQS>(be.jac, pivot, b);

        // update our current guess for the solution

        for (int n = 1; n <= BE_NEQS; n++) {
            y_np1(n) += b(n);
        }

        // check to see if we converged
        // we compute the norms

        Real y_norm = 0.0;
        Real b_norm = 0.0;
        for (int n = 1; n <= BE_NEQS; n++) {
            y_norm += y_np1(n) * y_np1(n);
            b_norm += b(n) * b(n);
        }
        y_norm = std::sqrt(y_norm / BE_NEQS);
        b_norm = std::sqrt(b_norm / BE_NEQS);

        if (b_norm < tol * y_norm) {
            converged = true;
            break;
        }

    }

    // we are done iterating -- did we converge?

    if (! converged) {

        if (ierr == BE_SUCCESS) {

            // if we didn't set another error, then we probably ran
            // out of iterations, so set nonconvergence

            ierr = BE_NONCONVERGENCE;
        }


    } else {

        // update the state to the latest values

        for (int n = 1; n <= BE_NEQS; n++) {
            be.y(n) = y_np1(n);
        }

    }

    // update the state -- if we didn't converge, then this just
    // restores the old state

    for (int n = 0; n < NumSpec; n++) {
        state.xn[n] = be.y(n+1);
    }
    state.e = be.y(net_ienuc);

    return ierr;

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int be_integrator (burn_t& state, be_t& be)
{

    be.n_rhs = 0;
    be.n_step = 0;

    int ierr;

    // estimate the timestep

    Array1D<Real, 1, neqs> ydot;
    evaluate_rhs(state, ydot);

    be.n_rhs += 1;

    Real dt_sub = initial_dt(state, be.tout, ydot);

    // main timestepping loop

    while (be.t < (1.0_rt - timestep_safety_factor) * be.tout &&
           be.n_step < ode_max_steps) {

        // don't go too far

        if (be.t + dt_sub > be.tout) {
            dt_sub = be.tout - be.t;
        }

        // our strategy is to take 2 steps at dt/2 and one at dt and
        // to compute the error from those


        // try to take a step dt

        burn_t state_coarse_dt = state;

        // first do 2 (fine) dt/2 steps

        ierr = single_step(state, be, dt_sub/2);
        if (ierr == BE_SUCCESS) {
            ierr = single_step(state, be, dt_sub/2);
        }

        // now that single (coarse) dt step

        if (ierr == BE_SUCCESS) {
            ierr = single_step(state_coarse_dt, be, dt_sub);
        }

#ifdef VODELIKE_ERROR
        // define a weight for each variable to use in checking the error

        Array1D<Real, 1, BE_NEQS> w;
        for (int n = 1; n <= NumSpec; n++) {
            w(n) = 1.0_rt / (rtol_spec * std::abs(state.xn[n-1]) + atol_spec);
        }
        w(net_ienuc) = 1.0_rt / (rtol_enuc * std::abs(state.e) + atol_enuc);

        // now look for w |y_fine - y_coarse| < 1

        Real rel_error = 0.0_rt;
        for (int n = 1; n <= NumSpec; n++) {
            rel_error = amrex::max(rel_error, w(n) * std::abs(state.xn[n-1] - state_coarse_dt.xn[n-1]));
        }
        rel_error = amrex::max(rel_error, w(net_ienuc) * std::abs(state.e - state_coarse_dt.e));

        bool step_success = false;
        if (rel_error < 1.0_rt) {
            step_success = true;
        }
#else

        // compute a relative error only for the species that are more
        // abundant that atol_spec

        Real rel_error = 0.0_rt;
        for (int n = 1; n <= NumSpec; n++) {
            if (state.xn[n-1] > atol_spec) {
                rel_error = amrex::max(rel_error,
                                       std::abs((state.xn[n-1] - state_coarse_dt.xn[n-1]) /
                                                (state.xn[n-1] + atol_spec)));
            }
        }
        rel_error = amrex::max(rel_error,
                               std::abs((state.e - state_coarse_dt.e) / state.e));

        bool step_success = false;
        if (rel_error < rtol_spec) {
            step_success = true;
        }

#endif


        if (ierr == BE_SUCCESS) {

            // there were no problems with the mechanics of the step
            // (linear algebra, etc.).  Now check if we met our error
            // goals.

            // state has the current best solution

            if (step_success) {
                be.t += dt_sub;

                // can we potentially increase the timestep?
                // backward-Euler has a local truncation error of dt**2

#ifdef VODELIKE_ERROR
                Real dt_new = dt_sub * std::pow(1.0_rt / rel_error, 0.5_rt);
#else
                Real dt_new = dt_sub * std::pow(rtol_spec / rel_error, 0.5_rt);
#endif
                dt_sub = amrex::min(amrex::max(dt_new, dt_sub / 2.0), 2.0 * dt_sub);

            } else {

                // we did not make our error goals, so just cut the timestep and leave the old solution in place

                dt_sub /= 2;

            }
        } else {

            // adjust the timestep and try again
            dt_sub /= 2;

        }

        ++be.n_step;
        ++be.n_step;

    }

    if (be.n_step >= ode_max_steps) {
        ierr = BE_TOO_MANY_STEPS;
    }

    if (ierr == BE_SUCCESS) {
        // do a final clean_state on the solution

        clean_state(state);
    }

    return ierr;

}

#endif
