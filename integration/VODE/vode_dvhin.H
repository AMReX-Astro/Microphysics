#ifndef _vode_dvhin_H_
#define _vode_dvhin_H_

#include <vode_rhs.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dvhin (burn_t& state, dvode_t& vstate, Real& H0, int& NITER, int& IER)
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

    for (int i = 1; i <= VODE_NEQS; ++i) {
        Real DELYI = PT1 * std::abs(vstate.yh(i,1)) + vstate.atol(i);
        Real AFI = std::abs(vstate.yh(i,2));
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
            for (int i = 1; i <= VODE_NEQS; ++i) {
                vstate.y(i) = vstate.yh(i,1) + H * vstate.yh(i,2);
            }

            const Real t1 = vstate.t + H;

            rhs(t1, state, vstate, vstate.acor);

            for (int i = 1; i <= VODE_NEQS; ++i) {
                vstate.acor(i) = (vstate.acor(i) - vstate.yh(i,2)) / H;
            }

            Real YDDNRM = 0.0_rt;
            for (int i = 1; i <= VODE_NEQS; ++i) {
                YDDNRM += (vstate.acor(i) * vstate.ewt(i)) * (vstate.acor(i) * vstate.ewt(i));
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
