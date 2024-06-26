#ifndef VODE_DVJAC_H
#define VODE_DVJAC_H

#include <vode_type.H>
#ifndef NEW_NETWORK_IMPLEMENTATION
#include <linpack.H>
#endif
#ifdef STRANG
#include <integrator_rhs_strang.H>
#endif
#ifdef SDC
#include <integrator_rhs_sdc.H>
#endif

template <typename BurnT, typename DvodeT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dvjac (int& IERPJ, BurnT& state, DvodeT& vstate)
{
    // dvjac is called by dvnlsd to compute and process the matrix
    // P = I - h*rl1*J , where J is an approximation to the Jacobian
    // that we obtain either through direct evaluation or caching from
    // a previous evaluation. P is then subjected to LU decomposition
    // in preparation for later solution of linear systems with P as
    // coefficient matrix. This is done by DGEFA.

    constexpr int int_neqs = integrator_neqs<BurnT>();

    IERPJ = 0;

#ifdef ALLOW_JACOBIAN_CACHING
    // See whether the Jacobian should be evaluated. Start by basing
    // the decision on whether we're caching the Jacobian.

    int evaluate_jacobian = 1;

    if (use_jacobian_caching) {
        evaluate_jacobian = 0;
    }

    if (use_jacobian_caching) {
        // Now evaluate the cases where we're caching the Jacobian but aren't
        // going to be using the cached Jacobian.

        // On the first step we don't have a cached Jacobian. Also, after enough
        // steps, we consider the cached Jacobian too old and will want to re-evaluate
        // it, so we look at whether the step of the last Jacobian evaluation (NSLJ)
        // is more than max_steps_between_jacobian_evals steps in the past.
        if (vstate.n_step == 0 || vstate.n_step > vstate.NSLJ + max_steps_between_jacobian_evals) {
            evaluate_jacobian = 1;
        }

        // See the non-linear solver for details on these conditions.
        if (vstate.ICF == 1 && vstate.DRC < CCMXJ) {
            evaluate_jacobian = 1;
        }

        if (vstate.ICF == 2) {
            evaluate_jacobian = 1;
        }

    }

    if (evaluate_jacobian == 1) {
#endif

        // We want to evaluate the Jacobian -- now the path depends on
        // whether we're using the numerical or analytic Jacobian.

        if (vstate.jacobian_type == 1) {

            // For the analytic Jacobian, call the user-supplied function.

            // Increment the Jacobian evaluation counter.
            vstate.n_jac += 1;

            // Refresh the timestep marker for the last Jacobian evaluation.
            vstate.NSLJ = vstate.n_step;

            // Indicate that the Jacobian is current for this solve.
            vstate.JCUR = 1;

            // Initialize the Jacobian to zero
            vstate.jac.zero();

            jac(vstate.tn, state, vstate, vstate.jac);

#ifdef ALLOW_JACOBIAN_CACHING
            // Store the Jacobian if we're caching.
            if (use_jacobian_caching == 1) {
                vstate.jac_save = vstate.jac;
            }
#endif

        }
        else {

            // For the numerical Jacobian, make N calls to the RHS to approximate it.

            // Increment the Jacobian evaluation counter.
            vstate.n_jac += 1;

            // Refresh the timestep marker for the last Jacobian evaluation.
            vstate.NSLJ = vstate.n_step;

            // Indicate that the Jacobian is current for this solve.
            vstate.JCUR = 1;

            amrex::Real fac = 0.0_rt;
            for (int i = 1; i <= int_neqs; ++i) {
                fac += (vstate.savf(i) * vstate.ewt(i)) * (vstate.savf(i) * vstate.ewt(i));
            }
            fac = std::sqrt(fac / int_neqs);

            amrex::Real R0 = 1000.0_rt * std::abs(vstate.H) * UROUND * int_neqs * fac;
            if (R0 == 0.0_rt) {
                R0 = 1.0_rt;
            }

            const bool in_jacobian = true;
            for (int j = 1; j <= int_neqs; ++j) {
                const amrex::Real yj = vstate.y(j);

                const amrex::Real R = amrex::max(std::sqrt(UROUND) * std::abs(yj), R0 / vstate.ewt(j));
                vstate.y(j) += R;
                fac = 1.0_rt / R;

                rhs(vstate.tn, state, vstate, vstate.acor, in_jacobian);
                for (int i = 1; i <= int_neqs; ++i) {
                    vstate.jac.set(i, j, (vstate.acor(i) - vstate.savf(i)) * fac);
                }

                vstate.y(j) = yj;
            }

            // Increment the RHS evaluation counter by N.
            vstate.n_rhs += int_neqs;

#ifdef ALLOW_JACOBIAN_CACHING
            // Store the Jacobian if we're caching.
            if (use_jacobian_caching == 1) {
                vstate.jac_save = vstate.jac;
            }
#endif

        }

#ifdef ALLOW_JACOBIAN_CACHING
    }
    else {

        // Load the cached Jacobian.

        // Indicate the Jacobian is not current for this step.
        vstate.JCUR = 0;
        vstate.jac = vstate.jac_save;

    }
#endif

    // Multiply Jacobian by a scalar, add the identity matrix
    // (along the diagonal), and do LU decomposition.

    const amrex::Real hrl1 = vstate.H * vstate.RL1;
    const amrex::Real con = -hrl1;

    vstate.jac.mul(con);
    vstate.jac.add_identity();

    int IER{};

#ifdef NEW_NETWORK_IMPLEMENTATION
    RHS::dgefa(vstate.jac);
    IER = 0;
#else
    if (integrator_rp::linalg_do_pivoting == 1) {
        constexpr bool allow_pivot{true};
        dgefa<int_neqs, allow_pivot>(vstate.jac, vstate.pivot, IER);
    } else {
        constexpr bool allow_pivot{false};
        dgefa<int_neqs, allow_pivot>(vstate.jac, vstate.pivot, IER);
    }
#endif

    if (IER != 0) {
        IERPJ = 1;
    }
}

#endif
