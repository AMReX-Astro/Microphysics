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
int single_step (BurnT& state, BeT& be, const amrex::Real time_start,
                 const amrex::Real dt)
{
    constexpr int int_neqs = integrator_neqs<BurnT>();

    int ierr = IERR_SUCCESS;
    bool converged = false;
    const amrex::Real time_end = time_start + dt;
    const amrex::Real time_save = be.t;

    // create our current guess for the solution -- just as a first
    // order explicit prediction

    amrex::Array1D<amrex::Real, 1, int_neqs> ydot;

    const bool do_clean_state = ! integrator_rp::use_number_densities;

    be.t = time_start;
    rhs(time_start, state, be, ydot, false, do_clean_state);

    be.n_rhs += 1;

    amrex::Array1D<amrex::Real, 1, int_neqs> y_old;

    for (int n = 1; n <= int_neqs; n++) {
        y_old(n) = be.y(n);
        be.y(n) = be.y(n) + dt * ydot(n);
    }

    // Newton loop

    for (int iter = 1; iter <= integrator_rp::max_iter; iter++) {

        // work with the current guess

        if (integrator_rp::do_corrector_validation && !integrator_rp::use_number_densities) {
#ifdef SDC
            const amrex::Real rho_current = state.rho_orig + time_end * state.ydot_a[SRHO];
            const amrex::Real thresh = integrator_rp::species_failure_tolerance * rho_current;
#else
            const amrex::Real thresh = integrator_rp::species_failure_tolerance;
#endif
            bool fail{};
            for (int i = 1; i <= NumSpec; ++i) {
                if (be.y(i) < -thresh) {
                    fail = true;
                    break;
                }
            }

            if (fail) {
                ierr = IERR_BAD_STATE_IN_CORRECTOR;
                break;
            }
        }

        // get the ydots for our current guess of y

        be.t = time_end;
        rhs(time_end, state, be, ydot, false, do_clean_state);
        be.n_rhs += 1;

        // construct the Jacobian

        if (be.jacobian_type == 1) {
            jac(time_end, state, be, be.jac, do_clean_state);
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
        }

        // reset the solution to the original on any failure
        for (int n = 1; n <= int_neqs; n++) {
            be.y(n) = y_old(n);
        }

    }

    be.t = time_save;

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
        ierr = single_step(state, be, be.t, dt_single_step);
        if (ierr == IERR_SUCCESS) {
            be.t = be.tout;
            ++be.n_step;
        }
        return ierr;
    }

    // estimate the timestep

    amrex::Array1D<amrex::Real, 1, int_neqs> ydot;
    rhs(be.t, state, be, ydot, false,
        ! integrator_rp::use_number_densities);

    be.n_rhs += 1;

    amrex::Real dt_sub = initial_react_dt(state, be, ydot);
    int n_attempt = 0;
    int n_reject = 0;
    int n_consecutive_reject = 0;
    bool previous_step_rejected = false;

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

        const amrex::Real t_old = be.t;

        // Estimate the local truncation error by comparing the backward
        // Euler candidate against a forward Euler candidate over the same dt.

        amrex::Array1D<amrex::Real, 1, int_neqs> y_forward;
        rhs(t_old, state, be, ydot, false,
            ! integrator_rp::use_number_densities);
        be.n_rhs += 1;

        for (int n = 1; n <= int_neqs; ++n) {
            y_forward(n) = y_old(n) + dt_sub * ydot(n);
        }

        ierr = single_step(state, be, t_old, dt_sub);

        amrex::Real rel_error = 0.0_rt;
        amrex::Real max_error = 0.0_rt;
        bool step_success = false;
        int err_idx = -1;
        const amrex::Real dt_attempt = dt_sub;
        const bool have_error_estimate = ierr == IERR_SUCCESS;
        if (have_error_estimate) {
            // define a weight for each variable to use in checking the error

            amrex::Array1D<amrex::Real, 1, int_neqs> w;
            for (int n = 1; n <= NumSpec; n++) {
                w(n) = 1.0_rt / (be.rtol_spec * std::abs(be.y(n)) + be.atol_spec);
            }
            w(net_ienuc) = 1.0_rt / (be.rtol_enuc * std::abs(be.y(net_ienuc)) + be.atol_enuc);

            // Now look for weighted |y_backward - y_forward| < 1.  Use the
            // same RMS norm style as VODE for the controller, while keeping
            // the maximum component around for diagnostics.

            for (int n = 1; n <= NumSpec; n++) {
                const amrex::Real component_error = w(n) * std::abs(be.y(n) - y_forward(n));
                rel_error += component_error * component_error;
                if (component_error > max_error) {
                    max_error = component_error;
                    err_idx = n;
                }
            }
            const amrex::Real enuc_error = w(net_ienuc) * std::abs(be.y(net_ienuc) - y_forward(net_ienuc));
            rel_error += enuc_error * enuc_error;
            if (enuc_error > max_error) {
                max_error = enuc_error;
                err_idx = net_ienuc;
            }
            constexpr amrex::Real TQ_ORDER1 = 2.0_rt;
            rel_error = std::sqrt(rel_error / int_neqs) / TQ_ORDER1;

            if (rel_error < 1.0_rt) {
                step_success = true;
            }

            if (step_success && integrator_rp::use_number_densities) {
                for (int n = 1; n <= NumSpec; ++n) {
                    if (be.y(n) < -integrator_rp::species_failure_tolerance) {
                        step_success = false;
                        break;
                    }
                }
            }
        }

        if (ierr == IERR_SUCCESS && step_success) {

            // there were no problems with the mechanics of the step
            // (linear algebra, etc.) and we met our error
            // goals.

            be.t += dt_sub;

            // Can we potentially increase the timestep?  Backward-Euler has
            // a local truncation error of O(dt**2).  Use a VODE-like biased
            // growth factor and threshold, but continue to estimate the error
            // from the forward Euler / backward Euler comparison.

            constexpr amrex::Real ADDON = 1.0e-6_rt;
            constexpr amrex::Real BIAS2 = 2.0_rt;
            constexpr amrex::Real ETAMX = 10.0_rt;
            constexpr amrex::Real THRESH = 1.5_rt;
            constexpr amrex::Real tiny_error = 1.0e-30_rt;

            amrex::Real dt_factor = 1.0_rt;
            const amrex::Real safe_error = amrex::max(rel_error, tiny_error);
            const amrex::Real eta = 1.0_rt / (std::sqrt(BIAS2 * safe_error) + ADDON);

            if (! previous_step_rejected && eta >= THRESH) {
                dt_factor = amrex::min(eta, ETAMX);
            }

            dt_sub *= dt_factor;
            previous_step_rejected = false;
            n_consecutive_reject = 0;

#ifndef AMREX_USE_GPU
            if (integrator_rp::burner_verbose &&
                (be.n_step < 40 || be.n_step % 1000 == 0)) {
                std::cout << "[BE_STEP] accept n=" << be.n_step
                          << " t_new=" << be.t
                          << " dt=" << dt_attempt
                          << " rel_error=" << rel_error
                          << " max_error=" << max_error
                          << " err_idx=" << err_idx
                          << " eta=" << eta
                          << " dt_factor=" << dt_factor
                          << " dt_next=" << dt_sub
                          << " y=";
                for (int n = 1; n <= NumSpec; ++n) {
                    std::cout << be.y(n) << " ";
                }
                std::cout << std::endl;
            }
#endif

            ++be.n_step;

        } else {

            // roll back the solution
            for (int n = 1; n <= int_neqs; ++n) {
                be.y(n) = y_old(n);
            }

            // Adjust the timestep and try again.  This mirrors VODE's
            // same-order error-failure step-size update for order 1, using
            // the forward Euler / backward Euler error estimate.

            constexpr amrex::Real ADDON = 1.0e-6_rt;
            constexpr amrex::Real BIAS2 = 2.0_rt;
            constexpr amrex::Real ETAMIN = 0.1_rt;
            constexpr amrex::Real ETAMXF = 0.2_rt;
            constexpr amrex::Real tiny_error = 1.0e-30_rt;

            amrex::Real dt_factor = ETAMIN;
            if (have_error_estimate) {
                const amrex::Real safe_error = amrex::max(rel_error, tiny_error);
                dt_factor = 1.0_rt / (std::sqrt(BIAS2 * safe_error) + ADDON);
                dt_factor = amrex::max(dt_factor, ETAMIN);
                if (n_consecutive_reject >= 1) {
                    dt_factor = amrex::min(dt_factor, ETAMXF);
                }
            }

            dt_sub *= dt_factor;
            previous_step_rejected = true;
            ++n_consecutive_reject;
            ++n_reject;

#ifndef AMREX_USE_GPU
            if (integrator_rp::burner_verbose &&
                (be.n_step < 40 || be.n_step % 1000 == 0)) {
                std::cout << "[BE_STEP] reject n=" << be.n_step
                          << " t=" << be.t
                          << " dt=" << dt_attempt
                          << " dt_next=" << dt_sub
                          << " ierr=" << ierr
                          << " have_error_estimate=" << have_error_estimate
                          << " rel_error=" << rel_error
                          << " max_error=" << max_error
                          << " err_idx=" << err_idx
                          << " dt_factor=" << dt_factor
                          << " consecutive_rejects=" << n_consecutive_reject
                          << std::endl;
            }
#endif

        }

        ++n_attempt;

    }

    if (be.n_step >= integrator_rp::ode_max_steps) {
        ierr = IERR_TOO_MANY_STEPS;
    }

#ifndef AMREX_USE_GPU
    if (integrator_rp::burner_verbose) {
        std::cout << "[BE_SUMMARY] attempts=" << n_attempt
                  << " accepts=" << be.n_step
                  << " rejects=" << n_reject
                  << std::endl;
    }
#endif

    return ierr;

}

#endif
