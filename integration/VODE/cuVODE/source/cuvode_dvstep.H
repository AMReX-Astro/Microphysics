#ifndef _cuvode_dvstep_H_
#define _cuvode_dvstep_H_

#include <cuvode_parameters.H>
#include <cuvode_type.H>
#include <cuvode_dvset.H>
#include <cuvode_dvjust.H>
#include <cuvode_dvnlsd.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void advance_nordsieck (dvode_t& vstate)
{

    // Effectively multiplies the Nordsieck history
    // array by the Pascal triangle matrix.

    for (int k = vstate.NQ; k >= 1; --k) {
        for (int j = k; j <= vstate.NQ; ++j) {
            for (int i = 0; i < VODE_NEQS; ++i) {
                vstate.YH[j-1][i] += vstate.YH[j][i];
            }
        }
    }

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void retract_nordsieck (dvode_t& vstate)
{

    // Undoes the Pascal triangle matrix multiplication
    // implemented in subroutine advance_nordsieck.

    for (int k = vstate.NQ; k >= 1; --k) {
        for (int j = k; j <= vstate.NQ; ++j) {
            for (int i = 0; i < VODE_NEQS; ++i) {
                vstate.YH[j-1][i] -= vstate.YH[j][i];
            }
        }
    }

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dvstep (Real* pivot, dvode_t& vstate)
{

    // dvstep performs one step of the integration of an initial value
    // problem for a system of ordinary differential equations.
    // dvstep calls subroutine dvnlsd for the solution of the nonlinear system
    // arising in the time step.  Thus it is independent of the problem
    // Jacobian structure and the type of nonlinear system solution method.
    // dvstep returns a completion flag KFLAG.

    // A return with KFLAG = -1 or -2 means either abs(H) = HMIN or 10
    // consecutive failures occurred. On a return with KFLAG negative,
    // the values of TN and the YH array are as of the beginning of the last
    // step, and H is the last step size attempted.

    // On a successful return, ETAMAX is reset and ACOR is scaled.

    ! Declare local variables
    real(rt) :: CNQUOT, DDN, DSM, DUP, TOLD
    real(rt) :: ETAQ, ETAQM1, ETAQP1, FLOTL, R
    integer  :: I, I1, I2, IBACK, J, JB, NCF, NFLAG

    // Parameter declarations
    const int KFC = -3;
    const int KFH = -7
    const int MXNCF = 10;
    const Real ADDON = 1.0e-6_rt;
    const Real BIAS1 = 6.0e0_rt;
    const Real BIAS2 = 6.0e0_rt;
    const Real BIAS3 = 10.0e0_rt;
    const Real ETACF = 0.25e0_rt;
    const Real ETAMIN = 0.1e0_rt;
    const Real ETAMXF = 0.2e0_rt;
    const Real ETAMX1 = 1.0e4_rt;
    const Real ETAMX2 = 10.0e0_rt;
    const Real ETAMX3 = 10.0e0_rt;
    const Real ONEPSM = 1.00001e0_rt;
    const Real THRESH = 1.5e0_rt;

    Real ETAQ = 1.0_rt;
    Real ETAQM1 = 1.0_rt;

    vstate.KFLAG = 0;
    Real TOLD = vstate.TN;
    int NCF = 0;
    vstate.JCUR = 0;
    int NFLAG = 0;

    bool do_initialization = false;

    if (vstate.JSTART == 0) {

        // On the first call, the order is set to 1, and other variables are
        // initialized.  ETAMAX is the maximum ratio by which H can be increased
        // in a single step.  It is normally 10, but is larger during the
        // first step to compensate for the small initial H.  If a failure
        // occurs (in corrector convergence or error test), ETAMAX is set to 1
        // for the next increase.

        vstate.NQ = 1;
        vstate.L = 2;
        vstate.NQNYH = vstate.NQ * VODE_NEQS;
        vstate.TAU[0] = vstate.H;
        vstate.PRL1 = 1.0_rt;
        vstate.RC = 0.0_rt;
        vstate.ETAMAX = ETAMX1;
        vstate.NQWAIT = 2;
        vstate.HSCAL = vstate.H;

    }
    else {

        // Take preliminary actions on a normal continuation step (JSTART > 0).
        // If the driver changed H, then ETA must be reset and NEWH set to 1.
        // If a change of order was dictated on the previous step, then
        // it is done here and appropriate adjustments in the history are made.
        // On an order decrease, the history array is adjusted by DVJUST.
        // On an order increase, the history array is augmented by a column.
        //  On a change of step size H, the history array YH is rescaled.

        if (vstate.NEWH == 0) {
            do_initialization = false;
        }
        else {

            if (vstate.NEWQ < vstate.NQ) {
                dvjust(-1, vstate);
                vstate.NQ = vstate.NEWQ;
                vstate.L = vstate.NQ + 1;
                vstate.NQWAIT = vstate.L;
            }
            else if (vstate.NEWQ > vstate.NQ) {
                dvjust(1, vstate);
                vstate.NQ = vstate.NEWQ;
                vstate.L = vstate.NQ + 1;
                vstate.NQWAIT = vstate.L;
            }

            do_initialization = true;
        }

    }

    if (do_initialization) {

        // Rescale the history array for a change in H by a factor of ETA.

        R = 1.0_rt;

        for (int j = 2; j <= vstate.L; ++j) {
            R *= vstate.ETA;
            for (int i = 0; i < VODE_NEQS; ++i) {
                vstate.YH[j-1][i] *= R;
            }
        }

        vstate.H = vstate.HSCAL * vstate.ETA;
        vstate.HSCAL = vstate.H;
        vstate.RC = vstate.RC * vstate.ETA;
        vstate.NQNYH = vstate.NQ * VODE_NEQS;

    }

    // Compute the predicted values by effectively
    // multiplying the YH array by the Pascal triangle matrix.
    // dvset is called to calculate all integration coefficients.
    // RC is the ratio of new to old values of the coefficient H / EL(2) = h / l1.

    while (true) {

        vstate.TN += vstate.H;

        advance_nordsieck(vstate);

        dvset(vstate);

        vstate.RL1 = 1.0_rt / vstate.EL[1];
        vstate.RC *= (vstate.RL1 / vstate.PRL1);
        vstate.PRL1 = vstate.RL1;

        // Call the nonlinear system solver.

        dvnlsd(pivot, NFLAG, vstate);

        if (NFLAG != 0) {

            // The VNLS routine failed to achieve convergence (NFLAG /= 0).
            // The YH array is retracted to its values before prediction.
            // The step size H is reduced and the step is retried, if possible.
            // Otherwise, an error exit is taken.

            NCF += 1;
            vstate.ETAMAX = 1.0_rt;
            vstate.TN = TOLD;

            retract_nordsieck(vstate);

            if (NFLAG < -1) {
                if (NFLAG == -2) {
                    vstate.KFLAG = -3;
                }
                else if (NFLAG == -3) {
                    vstate.KFLAG = -4;
                }
                vstate.JSTART = 1;
                return;
            }

            if (std::abs(vstate.H) <= HMIN * ONEPSM) {
                vstate.KFLAG = -2;
                vstate.JSTART = 1;
                return;
            }

            if (NCF == MXNCF) {
                vstate.KFLAG = -2;
                vstate.JSTART = 1;
                return;
            }

            vstate.ETA = ETACF;
            vstate.ETA = amrex::max(vstate.ETA, HMIN / std::abs(vstate.H));
            NFLAG = -1;

            // Rescale the history array for a change in H by a factor of ETA.
            Real R = 1.0_rt;

            for (int j = 2; j <= vstate.L; ++j) {
                R *= vstate.ETA;
                for (int i = 0; i < VODE_NEQS; ++i) {
                    vstate.YH[j-1][i] *= R;
                }
            }

            vstate.H = vstate.HSCAL * vstate.ETA;
            vstate.HSCAL = vstate.H;
            vstate.RC *= vstate.ETA;
            vstate.NQNYH = vstate.NQ * VODE_NEQS;

            continue;
        }

        // The corrector has converged (NFLAG = 0). The local error test is
        // made and control passes to statement 500 if it fails.

        Real DSM = vstate.ACNRM / vstate.TQ[1];
        if (DSM <= 1.0_rt) {

            // After a successful step, update the YH and TAU arrays and decrement
            // NQWAIT. If NQWAIT is then 1 and NQ .lt. MAXORD, then ACOR is saved
            // for use in a possible order increase on the next step.
            // If ETAMAX = 1 (a failure occurred this step), keep NQWAIT >= 2.

            vstate.KFLAG = 0;
            vstate.NST += 1;
            vstate.HU = vstate.H;
            for (int iback = 1; iback <= vstate.NQ; ++iback) {
                int i = vstate.L - IBACK - 1;
                vstate.TAU[i+1] = vstate.TAU[i];
            }

            vstate.TAU[0] = vstate.H;
            for (int j = 1; j <= vstate.L; ++j) {
                for (int i = 0; i < VODE_NEQS; ++i) {
                    vstate.yh[j-1][i] += vstate.EL[j-1] * vstate.acor[i];
                }
            }

            vstate.NQWAIT -= 1;
            if ((vstate.L != VODE_LMAX) && (vstate.NQWAIT == 1)) {
                for (int n = 0; n < VODE_NEQS; ++n) {
                    vstate.yh[VODE_LMAX-1][n] = vstate.acor[n];
                }
                vstate.CONP = vstate.TQ[4];
            }

            if (vstate.ETAMAX != 1.0_rt) {
                break;
            }

            if (vstate.NQWAIT < 2) {
                vstate.NQWAIT = 2;
            }

            vstate.NEWQ = vstate.NQ;
            vstate.NEWH = 0;
            vstate.ETA = 1.0_rt;
            vstate.HNEW = vstate.H;
            vstate.ETAMAX = ETAMX3;
            if (vstate.NST <= 10) {
                vstate.ETAMAX = ETAMX2;
            }

            Real R = 1.0_rt / vstate.TQ[1];
            for (int n = 0; n < VODE_NEQS; ++n) {
                vstate.acor[n] *= R;
            }

            vstate.JSTART = 1;

            return;

        }

        // The error test failed. KFLAG keeps track of multiple failures.
        // Restore TN and the YH array to their previous values, and prepare
        // to try the step again. Compute the optimum step size for the
        // same order. After repeated failures, H is forced to decrease
        // more rapidly.

        vstate.KFLAG -= 1;
        NFLAG = -2;
        vstate.TN = TOLD;

        retract_nordsieck(vstate);

        if (std::abs(vstate.H) <= HMIN * ONEPSM) {
            vstate.KFLAG = -1;
            vstate.JSTART = 1;
            return;
        }

        vstate.ETAMAX = 1.0_rt;

        if (vstate.KFLAG > KFC) {

            // Compute ratio of new H to current H at the current order.
            Real FLOTL = vstate.L;
            vstate.ETA = 1.0_rt / (std::pow(BIAS2 * DSM, 1.0_rt / FLOTL) + ADDON);
            vstate.ETA = amrex::max(vstate.ETA, HMIN / std::abs(vstate.H), ETAMIN);
            if ((vstate.KFLAG <= -2) && (vstate.ETA > ETAMXF)) {
                vstate.ETA = ETAMXF;
            }

            // Rescale the history array for a change in H by a factor of ETA.
            Real R = 1.0_rt;

            for (int j = 2; j <= vstate.L; ++j) {
                R *= vstate.ETA;
                for (int n = 0; n < VODE_NEQS; ++n) {
                    vstate.YH[j-1][n] *= R;
                }
            }

            vstate.H = vstate.HSCAL * vstate.ETA;
            vstate.HSCAL = vstate.H;
            vstate.RC = vstate.RC * vstate.ETA;
            vstate.NQNYH = vstate.NQ * VODE_NEQS;

            continue;

        }

        // Control reaches this section if 3 or more consecutive failures
        // have occurred. It is assumed that the elements of the YH array
        // have accumulated errors of the wrong order. The order is reduced
        // by one, if possible. Then H is reduced by a factor of 0.1 and
        // the step is retried.  After a total of 7 consecutive failures,
        // an exit is taken with KFLAG = -1.

        if (vstate.KFLAG == KFH) {
            vstate.KFLAG = -1;
            vstate.JSTART = 1;
            return;
        }

        if (vstate.NQ /= 1) {
            vstate.ETA = amrex::max(ETAMIN, HMIN / std::abs(vstate.H));
            dvjust(-1, vstate);
            vstate.L = vstate.NQ;
            vstate.NQ = vstate.NQ - 1;
            vstate.NQWAIT = vstate.L;

            // Rescale the history array for a change in H by a factor of ETA.
            Real R = 1.0_rt;

            for (int j = 2; j <= vstate.L; ++j) {
                R *= vstate.ETA;
                for (int n = 0; n < VODE_NEQS; ++n) {
                    vstate.YH[j-1][n] *= R;
                }
            }

            vstate.H *= vstate.ETA;
            vstate.HSCAL = vstate.H;
            vstate.RC = vstate.RC * vstate.ETA;
            vstate.NQNYH = vstate.NQ * VODE_NEQS;

            continue;

        }

        vstate.ETA = amrex::max(ETAMIN, HMIN / std::abs(vstate.H));
        vstate.H *= vstate.ETA;
        vstate.HSCAL = vstate.H;
        vstate.TAU[0] = vstate.H;
        f_rhs(vstate.TN, vstate, vstate.savf);
        vstate.NFE += 1;
        for (int n = 0; n < VODE_NEQS; ++n) {
            vstate.YH[1][n] = vstate.H * vstate.savf[n];
        }

        vstate.NQWAIT = 10;

    }

    // If NQWAIT = 0, an increase or decrease in order by one is considered.
    // Factors ETAQ, ETAQM1, ETAQP1 are computed by which H could
    // be multiplied at order q, q-1, or q+1, respectively.
    // The largest of these is determined, and the new order and
    // step size set accordingly.
    // A change of H or NQ is made only if H increases by at least a
    // factor of THRESH.  If an order change is considered and rejected,
    // then NQWAIT is set to 2 (reconsider it after 2 steps).

    bool already_set_eta = false;

    // Compute ratio of new H to current H at the current order.
    Real FLOTL = vstate.L;
    ETAQ = 1.0_rt / (std::pow(BIAS2 * DSM, 1.0_rt / FLOTL) + ADDON);
    if (vstate.NQWAIT == 0) {
        vstate.NQWAIT = 2;
        ETAQM1 = 0.0_rt;
        if (vstate.NQ != 1) {
            // Compute ratio of new H to current H at the current order less one.
            Real DDN = 0.0_rt;
            for (int n = 0; n < VODE_NEQS; ++n) {
                DDN += (vstate.YH[vstate.L-1][n] * vstate.YH[vstate.L-1][n]) * (vstate.EWT[n] * vstate.EWT[n]);
            }
            DDN = std::sqrt(DDN / VODE_NEQS) / vstate.TQ[0];
            ETAQM1 = 1.0_rt / (std::pow(BIAS1 * DDN, 1.0_rt / (FLOTL - 1.0_rt)) + ADDON);
        }

        ETAQP1 = 0.0_rt;
        if (vstate.L != VODE_LMAX) {
            // Compute ratio of new H to current H at current order plus one.
            Real CNQUOT = (vstate.TQ[4] / vstate.CONP) * std::pow(vstate.H / vstate.TAU[1], vstate.L);
            for (int n = 0; n < VODE_NEQS; ++n) {
                vstate.savf[n] = vstate.acor[n] - CNQUOT * vstate.YH[VODE_LMAX-1][n];
            }
            Real DUP = 0.0_rt;
            for (int n = 0; n < VODE_NEQS; ++n) {
                DUP += (vstate.SAVF[n] * vstate.SAVF[n]) * (vstate.EWT[n] * vstate.EWT[n]);
            }
            DUP = std::sqrt(DUP / VODE_NEQS) / vstate.TQ[2];
            ETAQP1 = 1.0_rt / (std::pow(BIAS3 * DUP, 1.0_rt / (FLOTL + 1.0_rt)) + ADDON);
        }

        if (ETAQ < ETAQP1) {
            if (ETAQP1 > ETAQM1) {
                vstate.ETA = ETAQP1;
                vstate.NEWQ = vstate.NQ + 1;
                for (int n = 0; n < VODE_NEQS; ++n) {
                    vstate.YH[VODE_LMAX-1][n] = vstate.acor[n];
                }
            }
            else {
                vstate.ETA = ETAQM1;
                vstate.NEWQ = vstate.NQ - 1;
            }
            already_set_eta = true;
        }

        if (ETAQ < ETAQM1 && !already_set_eta) {
            vstate.ETA = ETAQM1;
            vstate.NEWQ = vstate.NQ - 1;
            already_set_eta = true;
        }

    }

    if (!already_set_eta) {
        vstate.ETA = ETAQ;
        vstate.NEWQ = vstate.NQ;
    }

    // Test tentative new H against THRESH and ETAMAX, and HMXI, then exit.
    if (vstate.ETA >= THRESH && vstate.ETAMAX /= 1.0_rt) {
        vstate.ETA = amrex::min(vstate.ETA, vstate.ETAMAX);
        vstate.ETA = vstate.ETA / amrex::max(1.0_rt, std::abs(vstate.H) * HMXI * vstate.ETA);
        vstate.NEWH = 1;
        vstate.HNEW = vstate.H * vstate.ETA;
        vstate.ETAMAX = ETAMX3;
        if (vstate.NST <= 10) {
            vstate.ETAMAX = ETAMX2;
        }
        Real R = 1.0_rt / vstate.TQ[1];
        for (int n = 0; n < VODE_NEQS; ++n) {
            vstate.acor[n] *= R;
        }
        vstate.JSTART = 1;
        return;
    }

    vstate.NEWQ = vstate.NQ;
    vstate.NEWH = 0;
    vstate.ETA = 1.0_rt;
    vstate.HNEW = vstate.H;
    vstate.ETAMAX = ETAMX3;
    if (vstate.NST <= 10) {
        vstate.ETAMAX = ETAMX2;
    }
    Real R = 1.0_rt / vstate.TQ[1];
    for (int n = 0; n < VODE_NEQS; ++n) {
        vstate.acor[n] *= R;
    }
    vstate.JSTART = 1;

}

#endif
