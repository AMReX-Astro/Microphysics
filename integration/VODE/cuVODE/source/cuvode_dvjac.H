#ifndef _cuvode_dvjac_H_
#define _cuvode_dvjac_H_

#include <cuvode_parameters.H>
#include <cuvode_type.H>
#include <cuvode_linpack.H>
#include <vode_rhs.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dvjac (int* pivot, int& IERPJ, dvode_t& vstate)
{

    // dvjac is called by dvnlsd to compute and process the matrix
    // P = I - h*rl1*J , where J is an approximation to the Jacobian
    // that we obtain either through direct evaluation or caching from
    // a previous evaluation. P is then subjected to LU decomposition
    // in preparation for later solution of linear systems with P as
    // coefficient matrix. This is done by DGEFA.

    IERPJ = 0;

    // See whether the Jacobian should be evaluated. Start by basing
    // the decision on whether we're caching the Jacobian.
    int evaluate_jacobian = -vstate.JSV;

    if (vstate.JSV == 1) {
        // Now evaluate the cases where we're caching the Jacobian but aren't
        // going to be using the cached Jacobian.

        // On the first step we don't have a cached Jacobian. Also, after enough
        // steps, we consider the cached Jacobian too old and will want to re-evaluate
        // it, so we look at whether the step of the last Jacobian evaluation (NSLJ)
        // is more than max_steps_between_jacobian_evals steps in the past.
        if (vstate.NST == 0 || vstate.NST > vstate.NSLJ + max_steps_between_jacobian_evals) {
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

        // We want to evaluate the Jacobian -- now the path depends on
        // whether we're using the numerical or analytic Jacobian.

        if (vstate.jacobian == 1) {

            // For the analytic Jacobian, call the user-supplied function.

            // Increment the Jacobian evaluation counter.
            vstate.NJE += 1;

            // Refresh the timestep marker for the last Jacobian evaluation.
            vstate.NSLJ = vstate.NST;

            // Indicate that the Jacobian is current for this solve.
            vstate.JCUR = 1;

            for (int n = 0; n < VODE_NEQS * VODE_NEQS; ++n) {
                vstate.jac[n] = 0.0_rt;
            }

            jac(vstate.tn, vstate, vstate.jac);

            // Store the Jacobian if we're caching.
            if (vstate.JSV == 1) {
                for (int n = 0; n < VODE_NEQS * VODE_NEQS; ++n) {
                    vstate.jac_save[n] = vstate.jac[n];
                }
            }

        }
        else {

            // For the numerical Jacobian, make N calls to the RHS to approximate it.

            // Increment the Jacobian evaluation counter.
            vstate.NJE += 1;

            // Refresh the timestep marker for the last Jacobian evaluation.
            vstate.NSLJ = vstate.NST;

            // Indicate that the Jacobian is current for this solve.
            vstate.JCUR = 1;

            Real fac = 0.0_rt;
            for (int n = 0; n < VODE_NEQS; ++n) {
                fac += (vstate.savf[n] * vstate.savf[n]) * (vstate.ewt[n] * vstate.ewt[n]);
            }
            fac = std::sqrt(fac / VODE_NEQS);

            Real R0 = 1000.0_rt * std::abs(vstate.H) * UROUND * VODE_NEQS * fac;
            if (R0 == 0.0_rt) {
                R0 = 1.0_rt;
            }

            int j1 = 0;
            for (int j = 0; j < VODE_NEQS; ++j) {
                Real yj = vstate.y[j];

                Real R = amrex::max(SRUR * std::abs(yj), R0 / vstate.ewt[j]);
                vstate.y[j] += R;
                fac = 1.0_rt / R;

                f_rhs(vstate.tn, vstate, vstate.acor);
                for (int n = 0; n < VODE_NEQS; ++n) {
                    vstate.jac[n+j1] = (vstate.acor[n] - vstate.savf[n]) * fac;
                }

                vstate.y[j] = yj;

                j1 += VODE_NEQS;
            }

            // Increment the RHS evaluation counter by N.
            vstate.NFE += VODE_NEQS;

            // Store the Jacobian if we're caching.
            if (vstate.JSV == 1) {
                for (int n = 0; n < VODE_NEQS * VODE_NEQS; ++n) {
                    vstate.jac_save[n] = vstate.jac[n];
                }
            }

        }

    }
    else {

        // Load the cached Jacobian.

        // Indicate the Jacobian is not current for this step.
        vstate.JCUR = 0;
        for (int n = 0; n < VODE_NEQS * VODE_NEQS; ++n) {
            vstate.jac[n] = vstate.jac_save[n];
        }

    }

    // Multiply Jacobian by a scalar, add the identity matrix
    // (along the diagonal), and do LU decomposition.

    Real hrl1 = vstate.H * vstate.RL1;
    Real con = -hrl1;
    for (int n = 0; n < VODE_NEQS * VODE_NEQS; ++n) {
        vstate.jac[n] *= con;
    }

    int j = 0;
    for (int n = 0; n < VODE_NEQS; ++n) {
        // Add 1 from the identity matrix; then,
        // advance to the next location on the diagonal,
        // which is in the next row (advance by VODE_NEQS)
        // and over by one column (advance by 1).
        vstate.jac[j] += 1.0_rt;
        j += VODE_NEQS + 1;
    }

    int IER;
    dgefa(vstate.jac, pivot, IER);

    if (IER != 0) {
        IERPJ = 1;
    }

}

#endif
