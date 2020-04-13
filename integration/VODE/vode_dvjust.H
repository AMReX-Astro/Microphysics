#ifndef _vode_dvjust_H_
#define _vode_dvjust_H_

#include <vode_type.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dvjust (int IORD, dvode_t& vstate)
{

    // This subroutine adjusts the YH array on reduction of order,
    // and also when the order is increased.

    if ((vstate.NQ == 2) && (IORD != 1)) {
        return;
    }

    int NQM1 = vstate.NQ - 1;
    int NQM2 = vstate.NQ - 2;

    Real HSUM, XI, ALPH0, ALPH1, PROD, XIOLD, T1;

    // Check to see if the order is being increased or decreased.

    if (IORD != 1) {

        // Order decrease.
        for (int j = 1; j <= VODE_LMAX; ++j) {
            vstate.el(j) = 0.0_rt;
        }
        vstate.el(3) = 1.0_rt;
        HSUM = 0.0_rt;
        for (int j = 1; j <= NQM2; ++j) {
            // Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)).
            HSUM += vstate.tau(j);
            XI = HSUM / vstate.HSCAL;
            for (int iback = 1; iback <= j+1; ++iback) {
                int i = (j + 4) - iback;
                vstate.el(i) = vstate.el(i) * XI + vstate.el(i-1);
            }
        }

        // Subtract correction terms from YH array.
        for (int j = 3; j <= vstate.NQ; ++j) {
            for (int i = 1; i <= VODE_NEQS; ++i) {
                vstate.yh(i,j) -= vstate.yh(i,vstate.L) * vstate.el(j);
            }
        }

    }
    else {

        // Order increase.
        for (int j = 1; j <= VODE_LMAX; ++j) {
            vstate.el(j) = 0.0_rt;
        }

        vstate.el(3) = 1.0_rt;
        ALPH0 = -1.0_rt;
        ALPH1 = 1.0_rt;
        PROD = 1.0_rt;
        XIOLD = 1.0_rt;
        HSUM = vstate.HSCAL;

        if (vstate.NQ != 1) {
            for (int j = 1; j <= NQM1; ++j) {
                // Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)).
                HSUM += vstate.tau(j+1);
                XI = HSUM / vstate.HSCAL;
                PROD *= XI;
                ALPH0 -= 1.0_rt / (j + 1);
                ALPH1 += 1.0_rt / XI;
                for (int iback = 1; iback <= j+1; ++iback) {
                    int i = (j + 4) - iback;
                    vstate.el(i) = vstate.el(i) * XIOLD + vstate.el(i-1);
                }
                XIOLD = XI;
            }
        }

        T1 = (-ALPH0 - ALPH1) / PROD;
        // Load column L + 1 in YH array.
        for (int i = 1; i <= VODE_NEQS; ++i) {
            vstate.yh(i,vstate.L+1) = T1 * vstate.yh(i,VODE_LMAX);
        }

        // Add correction terms to YH array.
        for (int j = 3; j <= vstate.NQ + 1; ++j) {
            for (int i = 1; i <= VODE_NEQS; ++i) {
                vstate.yh(i,j) += vstate.el(j) * vstate.yh(i,vstate.L+1);
            }
        }

    }

}

#endif
