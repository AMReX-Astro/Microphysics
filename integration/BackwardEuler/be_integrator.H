#ifndef BE_INTEGRATOR_H
#define BE_INTEGRATOR_H

#include <be_type.H>
#include <network.H>
#include <actual_network.H>
#include <actual_rhs.H>
#include <burn_type.H>
#include <linpack.H>
#include <numerical_jacobian.H>
#ifdef STRANG
#include <integrator_rhs_strang.H>
#endif
#ifdef SIMPLIFIED_SDC
#include <integrator_rhs_simplified_sdc.H>
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real initial_dt (burn_t& state, be_t& be,
                 Array1D<Real, 1, BE_NEQS>& ydot)
{
    // this is a version of the initial timestep estimator from VODE
    // we come in with be.y() storing the initial solution

    // initial lower and upper bounds on the timestep

    Real hL = 100.0_rt * std::numeric_limits<Real>::epsilon() * be.tout;
    Real hU = 0.1_rt * be.tout;

    // initial guess for the itteraction

    Real h = std::sqrt(hL * hU);
    Real h_old = 10.0_rt * h;

    // Iterate on ddydtt = (RHS(t + h, y + h * dydt) - dydt) / h

    Array1D<Real, 1, BE_NEQS> ewt;
    Array1D<Real, 1, BE_NEQS> ydot_temp;
    Array1D<Real, 1, BE_NEQS> ddydtt;

    // save the old state and update the stored state
    Array1D<Real, 1, BE_NEQS> y_old;
    for (int n = 1; n <= BE_NEQS; ++n) {
        y_old(n) = be.y(n);
    }

    for (int n = 1; n <= 4; n++) {

       h_old = h;

       // Get the error weighting -- this is similar to VODE's dewset
       // routine

       for (int ii = 1; ii <= NumSpec; ii++) {
           ewt(ii) = be.rtol_spec * std::abs(y_old(ii)) + be.atol_spec;
       }
       ewt(net_ienuc) = be.rtol_enuc * std::abs(y_old(net_ienuc)) + be.atol_enuc;

       // Construct the trial point.

       for (int ii = 1; ii <= BE_NEQS; ii++) {
           be.y(ii) = y_old(ii) + h * ydot(ii);
       }

       // Call the RHS, then estimate the finite difference.

       rhs(be.t, state, be, ydot_temp);

       for (int ii = 1; ii <= BE_NEQS; ii++) {
           ddydtt(ii) = (ydot_temp(ii) - ydot(ii)) / h;
       }

       Real yddnorm = 0.0_rt;
       for (int ii = 1; ii <= BE_NEQS; ii++) {
           yddnorm += std::pow(ddydtt(ii) * ewt(ii), 2);
       }
       yddnorm = std::sqrt(yddnorm / BE_NEQS);

       if (yddnorm*hU*hU > 2.0_rt) {
          h = std::sqrt(2.0_rt / yddnorm);
       } else {
          h = std::sqrt(h * hU);
       }

       if (h_old < 2.0_rt * h && h_old > 0.5_rt * h) {
           break;
       }

    }

    // Save the final timestep, with a bias factor.

    Real dt = h / 2.0_rt;
    dt = amrex::min(amrex::max(h, hL), hU);

    dt = amrex::min(dt, ode_max_dt);

    // restore the old time solution
    for (int n = 1; n <= BE_NEQS; ++n) {
        be.y(n) = y_old(n);
    }

    return dt;
}


///
/// update state.xn[] and state.e through a timestep dt
/// state is updated in place -- if we are unsuccessful, we reset it
/// we return an error code
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int single_step (burn_t& state, be_t& be, const Real dt)
{

    int ierr = BE_SUCCESS;
    bool converged = false;

    // create our current guess for the solution -- just as a first
    // order explicit prediction

    Array1D<Real, 1, neqs> ydot;

    rhs(be.t, state, be, ydot);

    be.n_rhs += 1;

    Array1D<Real, 1, BE_NEQS> y_old;

    for (int n = 1; n <= BE_NEQS; n++) {
        y_old(n) = be.y(n);
        be.y(n) = be.y(n) + dt * ydot(n);
    }

    // Newton loop

    for (int iter = 1; iter <= max_iter; iter++) {

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

        // construct the RHS of our linear system

        Array1D<Real, 1, BE_NEQS> b;
        for (int n = 1; n <= BE_NEQS; n++) {
            b(n) = y_old(n) - be.y(n) + dt * ydot(n);
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
            be.y(n) += b(n);
        }

        // check to see if we converged
        // we compute the norms

        Real y_norm = 0.0;
        Real b_norm = 0.0;
        for (int n = 1; n <= BE_NEQS; n++) {
            y_norm += be.y(n) * be.y(n);
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

            // reset the solution to the original
            for (int n = 1; n <= BE_NEQS; n++) {
                be.y(n) = y_old(n);
            }

        }

    }

    return ierr;

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int be_integrator (burn_t& state, be_t& be)
{

    be.n_rhs = 0;
    be.n_jac = 0;
    be.n_step = 0;

    int ierr;

    // estimate the timestep

    Array1D<Real, 1, neqs> ydot;
    rhs(be.t, state, be, ydot);

    be.n_rhs += 1;

    Real dt_sub = initial_dt(state, be, ydot);

    // main timestepping loop

    while (be.t < (1.0_rt - timestep_safety_factor) * be.tout &&
           be.n_step < ode_max_steps) {

        // store the current solution -- we'll revert to this if a step fails
        Array1D<Real, 1, BE_NEQS> y_old;
        for (int n = 1; n <= BE_NEQS; ++n) {
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

        Array1D<Real, 1, BE_NEQS> y_fine;

        ierr = single_step(state, be, dt_sub/2);
        if (ierr == BE_SUCCESS) {
            ierr = single_step(state, be, dt_sub/2);

            // store the fine dt solution

            for (int n = 1; n <= BE_NEQS; ++n) {
                y_fine(n) = be.y(n);
            }

            // now that single (coarse) dt step
            // first reset the solution
            for (int n = 1; n <= BE_NEQS; ++n) {
                be.y(n) = y_old(n);
            }
            ierr = single_step(state, be, dt_sub);
        }

        // define a weight for each variable to use in checking the error

        Array1D<Real, 1, BE_NEQS> w;
        for (int n = 1; n <= NumSpec; n++) {
            w(n) = 1.0_rt / (be.rtol_spec * std::abs(y_fine(n)) + be.atol_spec);
        }
        w(net_ienuc) = 1.0_rt / (be.rtol_enuc * std::abs(y_fine(net_ienuc)) + be.atol_enuc);

        // now look for w |y_fine - y_coarse| < 1

        Real rel_error = 0.0_rt;
        for (int n = 1; n <= NumSpec; n++) {
            rel_error = amrex::max(rel_error, w(n) * std::abs(y_fine(n) - be.y(n)));
        }
        rel_error = amrex::max(rel_error, w(net_ienuc) * std::abs(y_fine(net_ienuc) - be.y(net_ienuc)));

        bool step_success = false;
        if (rel_error < 1.0_rt) {
            step_success = true;
        }

        if (ierr == BE_SUCCESS && step_success) {

            // there were no problems with the mechanics of the step
            // (linear algebra, etc.) and we met our error
            // goals.

            // y_fine has the current best solution

            be.t += dt_sub;

            for (int n = 1; n <= BE_NEQS; ++n) {
                be.y(n) = y_fine(n);
            }

            // can we potentially increase the timestep?
            // backward-Euler has a local truncation error of dt**2

            Real dt_new = dt_sub * std::pow(1.0_rt / rel_error, 0.5_rt);
            dt_sub = amrex::min(amrex::max(dt_new, dt_sub / 2.0), 2.0 * dt_sub);

        } else {

            // roll back the solution
            for (int n = 1; n <= BE_NEQS; ++n) {
                be.y(n) = y_old(n);
            }

            // adjust the timestep and try again
            dt_sub /= 2;

        }

        ++be.n_step;
        ++be.n_step;

    }

    if (be.n_step >= ode_max_steps) {
        ierr = BE_TOO_MANY_STEPS;
    }

    return ierr;

}

#endif
