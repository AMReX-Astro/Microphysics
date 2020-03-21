#ifndef _cuvode_dvnlsd_H_
#define _cuvode_dvnlsd_H_

#include <cuvode_parameters.H>
#include <cuvode_type.H>
#include <cuvode_linpack.H>
#include <cuvode_dvjac.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dvnlsd (Real* pivot, int& NFLAG, dvode_t& vstate)
{

    // dvnlsd is a nonlinear system solver that uses a chord (modified
    // Newton) method with direct linear algebraic system solvers.
    // It then handles the corrector phase of this integration package.

    ! Declare local variables
    real(rt) :: CSCALE, DCON, DEL, DELP
    integer  :: I, IERPJ, M

    // Parameter declarations
    const Real CCMAX = 0.3e0_rt;
    const Real CRDOWN = 0.3e0_rt;
    const Real RDIV  = 2.0e0_rt;
    const int MAXCOR = 3;
    const int MSBP = 20;

    // On the first step, on a change of method order, or after a
    // nonlinear convergence failure with NFLAG = -2, set IPUP = 1
    // to force a Jacobian update.

    if (vstate.JSTART == 0) {
        vstate.NSLP = 0;
    }
    if (NFLAG == 0) {
        vstate.ICF = 0;
    }
    if (NFLAG == -2) {
        vstate.IPUP = 1;
    }
    if (vstate.JSTART == 0) {
        vstate.IPUP = 1;
    }

    // RC is the ratio of new to old values of the coefficient H / EL(2) = h / l1.
    // When RC differs from 1 by more than CCMAX, IPUP is set to 1
    // to force DVJAC to be called, if a Jacobian is involved.
    // In any case, DVJAC is called at least every MSBP steps.

    vstate.DRC = std::abs(vstate.RC - 1.0_rt);
    if (vstate.DRC > CCMAX || vstate.NST >= vstate.NSLP + MSBP) {
        vstate.IPUP = 1;
    }

    // Up to MAXCOR corrector iterations are taken.  A convergence test is
    // made on the r.m.s. norm of each correction, weighted by the error
    // weight array EWT.  The sum of the corrections is accumulated in the
    // array ACOR.  The YH array is not altered in the corrector loop.

    bool converged = false;

    while (true) {

        int M = 0;
        Real DELP = 0.0_rt;

        for (int n = 0; n < VODE_NEQS; ++n) {
            vstate.Y[n] = vstate.YH[0][n];
        }

        f_rhs(vstate.TN, vstate, vstate.savf);
        vstate.NFE += 1;

        if (vstate.IPUP == 1) {

            // If indicated, the matrix P = I - h*rl1*J is reevaluated and
            // preprocessed before starting the corrector iteration.  IPUP is set
            // to 0 as an indicator that this has been done.

            dvjac(pivot, IERPJ, vstate);
            vstate.IPUP = 0;
            vstate.RC = 1.0_rt;
            vstate.DRC = 0.0_rt;
            vstate.CRATE = 1.0_rt;
            vstate.NSLP = vstate.NST;

            // If matrix is singular, take error return to force cut in step size.
            if (IERPJ /= 0) {
                NFLAG = -1;
                vstate.ICF = 2;
                vstate.IPUP = 1;
                return;
            }

        }

        for (int n = 0; n < VODE_NEQS; ++n) {
            vstate.acor[n] = 0.0_rt;
        }

        // Corrector iteration loop.

        while (true) {

            // Compute the corrector error, and solve the linear system with
            // that as right-hand side and  P as coefficient matrix. The
            // correction is scaled by the factor 2/(1+RC) to account for
            // changes in h*rl1 since the last dvjac call.

            for (int n = 0; n < VODE_NEQS; ++n) {
                vstate.Y[n] = (vstate.RL1 * vstate.H) * vstate.SAVF[n] -
                              (vstate.RL1 * vstate.YH[1][n] + vstate.ACOR[n]);
            }

            dgesl(vstate.jac, pivot, vstate.Y);

            if (vstate.RC /= 1.0_rt) {
                Real CSCALE = 2.0_rt / (1.0_rt + vstate.RC);
                for (int n = 0; n < VODE_NEQS; ++n) {
                    vstate.Y[n] *= CSCALE;
                }
            }

#ifdef CLEAN_INTEGRATOR_CORRECTION
            // Clean the correction to Y. Use vstate.Y as scratch space.

            // Find the corrected Y: Yc = Y_previous + Y_delta
            for (int n = 0; n < VODE_NEQS; ++n) {
                vstate.Y[n] = vstate.Y[n] + (vstate.YH[0][n] + vstate.ACOR[n]);
            }

            // Clean the corrected Y: Yc' = clean(Yc)
            clean_state(vstate.Y, vstate.RPAR);

            // Find the cleaned correction: clean(Y_delta) = Yc' - Y_previous
            for (int n = 0; n < VODE_NEQS; ++n) {
                vstate.Y[n] = vstate.Y[n] - (vstate.YH[0][n] + vstate.ACOR[n]);
            }
#endif

            Real DEL = 0.0_rt;
            for (int n = 0; n < VODE_NEQS; ++n) {
                DEL += (vstate.Y[n] * vstate.Y[n]) * (vstate.EWT[n] * vstate.EWT[n]);
            }
            DEL = std::sqrt(DEL / VODE_NEQS);

            for (int n = 0; n < VODE_NEQS; ++n) {
                vstate.acor[n] += vstate.Y[n];
            }

            for (int n = 0; n < VODE_NEQS; ++n) {
                vstate.Y[n] = vstate.YH[0][n] + vstate.ACOR[n];
            }

            // Test for convergence. If M > 0, an estimate of the convergence
            // rate constant is stored in CRATE, and this is used in the test.

            if (M != 0) {
                vstate.CRATE = amrex::max(CRDOWN * vstate.CRATE, DEL / DELP);
            }

            Real DCON = DEL * amrex::min(1.0_rt, vstate.CRATE) / vstate.TQ[3];
            if (DCON <= 1.0_rt) {
                // we converged, exit the outer loop
                converged = true;
                break;
            }

            M += 1;
            if (M == MAXCOR) {
                //! exit the inner correction iteration
                break;
            }

            if (M >= 2 && DEL > RDIV * DELP) {
                // exit the inner correction iteration
                break;
            }

            DELP = DEL;
            f_rhs(vstate.TN, vstate, vstate.SAVF);
            vstate.NFE += 1;

        }

        if (converged) {
            break;
        }

        if (vstate.JCUR == 1) {
            NFLAG = -1;
            vstate.ICF = 2;
            vstate.IPUP = 1;
            return;
        }

        vstate.ICF = 1;
        vstate.IPUP = 1;

    }

    // Return for successful step.
    NFLAG = 0;
    vstate.JCUR = 0;
    vstate.ICF = 0;
    if (M == 0) {
        vstate.ACNRM = DEL;
    }

    if (M > 0) {
        vstate.ACNRM = 0.0_rt;
        for (int n = 0; n < VODE_NEQS; ++n) {
            vstate.ACNRM += (vstate.ACOR[n] * vstate.ACOR[n]) * (vstate.EWT[n] * vstate.EWT[n]);
        }
        vstate.ACNRM = std::sqrt(vstate.ACNRM / VODE_NEQS);
    }

}

#endif
