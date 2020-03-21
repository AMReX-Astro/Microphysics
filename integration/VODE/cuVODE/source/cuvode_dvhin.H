#ifndef _cuvode_dvhin_H_
#define _cuvode_dvhin_H_

AMREX_GPU_HOST_DEVICE
void dvhin (dvode_t& vstate, Real& H0, int& NITER, int& IER)
{
    // This routine computes the step size, H0, to be attempted on the
    // first step, when the user has not supplied a value for this.

    // First we check that TOUT - T0 differs significantly from zero. Then
    // an iteration is done to approximate the initial second derivative
    // and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1.
    // A bias factor of 1/2 is applied to the resulting h.
    // The sign of H0 is inferred from the initial values of TOUT and T0.

    const Real PT1 = 0.1e0_rt;

    NITER = 0;
    Real TDIST = std::abs(vstate.TOUT - vstate.T);
    Real TROUND = UROUND * amrex::max(std::abs(vstate.T), std::abs(vstate.TOUT));

    if (TDIST < 2.0_rt * TROUND) {
        // Error return for vstate.TOUT - vstate.T too small.
        IER = -1;
        return;
    }

    // Set a lower bound on h based on the roundoff level in vstate.T and vstate.TOUT.
    Real HLB = 100.0_rt * TROUND;

    // Set an upper bound on h based on vstate.TOUT-vstate.T and the initial Y and YDOT.
    Real HUB = PT1 * TDIST;

    for (int n = 0; n < VODE_NEQS; ++n) {
        Real DELYI = PT1 * std::abs(vstate.YH[0][n]) + vstate.ATOL[n];
        Real AFI = std::abs(vstate.YH[1][n]);
        if (AFI * HUB > DELYI) {
            HUB = DELYI / AFI;
        }
    }

    // Set initial guess for h as geometric mean of upper and lower bounds.
    int ITER = 0;
    Real HG = std::sqrt(HLB * HUB);

    // If the bounds have crossed, exit with the mean value.
    bool do_iterations = true;
    if (HUB < HLB) {
        H0 = HG;
        do_iterations = false;
    }

    if (do_iterations) {

        while (true) {

            // Estimate the second derivative as a difference quotient in f.
            H = std::copysign(HG, vstate.TOUT - vstate.T);
            Real T1 = vstate.T + H;
            for (int n = 0; n < VODE_NEQS; ++n) {
                vstate.Y[n] = vstate.YH[0][n] + H * vstate.YH[1][n];
            }

            f_rhs(T1, vstate, vstate.ACOR);

            for (int n = 0; n < VODE_NEQS; ++n) {
                vstate.ACOR[n] = (vstate.ACOR[n] - vstate.YH[1][n]) / H;
            }

            Real YDDNRM = 0.0_rt;
            for (int n = 0; n < VODE_NEQS; ++n) {
                YDDNRM += (vstate.ACOR[n] * vstate.ACOR[n]) * (vstate.EWT[n] * vstate.EWT[n]);
            }
            YDDNRM = std::sqrt(YDDNRM / VODE_NEQS);

            // Get the corresponding new value of h.
            Real HNEW;
            if (YDDNRM * HUB * HUB > 2.0_rt) {
                HNEW = std::sqrt(2.0_rt / YDDNRM);
            }
            else {
                HNEW = std::sqrt(HG * HUB);
            }
            ITER += 1;

            // Test the stopping conditions.
            // Stop if the new and previous h values differ by a factor of .lt. 2.
            // Stop if four iterations have been done.  Also, stop with previous h
            // if HNEW/HG .gt. 2 after first iteration, as this probably means that
            // the second derivative value is bad because of cancellation error.

            if (iter >= 4) {
                break;
            }

            Real HRAT = HNEW / HG;
            if ((HRAT > 0.5_rt) && (HRAT < 2.0_rt)) {
                break;
            }

            if ((ITER >= 2) && (HNEW > 2.0_rt * HG)) {
                HNEW = HG;
                break;
            }

            HG = HNEW;

        }

        // Iteration done. Apply bounds, bias factor, and sign. Then exit.
        H0 = HNEW * 0.5_rt;
        if (H0 < HLB) {
            H0 = HLB;
        }
        if (H0 > HUB) {
            H0 = HUB;
        }

    }

    Real H0 = std::copysign(H0, vstate.TOUT - vstate.T);
    NITER = ITER;
    IER = 0;

}

#endif
