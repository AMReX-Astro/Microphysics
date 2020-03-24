#ifndef _cuvode_dvhin_H_
#define _cuvode_dvhin_H_

#include <vode_rhs.H>

AMREX_GPU_HOST_DEVICE
void dvhin (dvode_t& vstate, Real& H0, int& NITER, int& IER)
{
    // This routine computes the step size, H0, to be attempted on the
    // first step, when the user has not supplied a value for this.

    // First we check that tout - T0 differs significantly from zero. Then
    // an iteration is done to approximate the initial second derivative
    // and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1.
    // A bias factor of 1/2 is applied to the resulting h.
    // The sign of H0 is inferred from the initial values of tout and T0.

    const Real PT1 = 0.1e0_rt;

    NITER = 0;
    Real TDIST = std::abs(vstate.tout - vstate.t);
    Real TROUND = UROUND * amrex::max(std::abs(vstate.t), std::abs(vstate.tout));

    if (TDIST < 2.0_rt * TROUND) {
        // Error return for vstate.tout - vstate.t too small.
        IER = -1;
        return;
    }

    // Set a lower bound on h based on the roundoff level in vstate.t and vstate.tout.
    Real HLB = 100.0_rt * TROUND;

    // Set an upper bound on h based on vstate.tout-vstate.t and the initial Y and YDOT.
    Real HUB = PT1 * TDIST;

    for (int n = 0; n < VODE_NEQS; ++n) {
        Real DELYI = PT1 * std::abs(vstate.yh[0][n]) + vstate.atol[n];
        Real AFI = std::abs(vstate.yh[1][n]);
        if (AFI * HUB > DELYI) {
            HUB = DELYI / AFI;
        }
    }

    // Set initial guess for h as geometric mean of upper and lower bounds.
    int iter = 0;
    Real HG = std::sqrt(HLB * HUB);

    // If the bounds have crossed, exit with the mean value.
    bool do_iterations = true;
    if (HUB < HLB) {
        H0 = HG;
        do_iterations = false;
    }

    Real hnew;

    if (do_iterations) {

        while (true) {

            // Estimate the second derivative as a difference quotient in f.
            Real H = std::copysign(HG, vstate.tout - vstate.t);
            Real T1 = vstate.t + H;
            for (int n = 0; n < VODE_NEQS; ++n) {
                vstate.y[n] = vstate.yh[0][n] + H * vstate.yh[1][n];
            }

            f_rhs(T1, vstate, vstate.acor);

            for (int n = 0; n < VODE_NEQS; ++n) {
                vstate.acor[n] = (vstate.acor[n] - vstate.yh[1][n]) / H;
            }

            Real YDDNRM = 0.0_rt;
            for (int n = 0; n < VODE_NEQS; ++n) {
                YDDNRM += (vstate.acor[n] * vstate.acor[n]) * (vstate.ewt[n] * vstate.ewt[n]);
            }
            YDDNRM = std::sqrt(YDDNRM / VODE_NEQS);

            // Get the corresponding new value of h.
            if (YDDNRM * HUB * HUB > 2.0_rt) {
                hnew = std::sqrt(2.0_rt / YDDNRM);
            }
            else {
                hnew = std::sqrt(HG * HUB);
            }
            iter += 1;

            // Test the stopping conditions.
            // Stop if the new and previous h values differ by a factor of .lt. 2.
            // Stop if four iterations have been done.  Also, stop with previous h
            // if hnew/HG .gt. 2 after first iteration, as this probably means that
            // the second derivative value is bad because of cancellation error.

            if (iter >= 4) {
                break;
            }

            Real HRAT = hnew / HG;
            if ((HRAT > 0.5_rt) && (HRAT < 2.0_rt)) {
                break;
            }

            if ((iter >= 2) && (hnew > 2.0_rt * HG)) {
                hnew = HG;
                break;
            }

            HG = hnew;

        }

        // Iteration done. Apply bounds, bias factor, and sign. Then exit.
        H0 = hnew * 0.5_rt;
        if (H0 < HLB) {
            H0 = HLB;
        }
        if (H0 > HUB) {
            H0 = HUB;
        }

    }

    H0 = std::copysign(H0, vstate.tout - vstate.t);
    NITER = iter;
    IER = 0;

}

#endif
