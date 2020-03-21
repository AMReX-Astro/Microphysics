#ifndef _cuvode_dvjust_H_
#define _cuvode_dvjust_H_

#include <cuvode_parameters.H>
#include <cuvode_type.H>

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

    // Check to see if the order is being increased or decreased.

    if (IORD != 1) {

        // Order decrease.
        for (int j = 0; j < VODE_LMAX; ++j) {
            vstate.EL[j] = 0.0_rt;
        }
        vstate.EL[2] = 1.0_rt;
        Real HSUM = 0.0_rt;
        for (int j = 1; j <= NQM2; ++j) {
            // Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)).
            HSUM += vstate.TAU[j-1];
            Real XI = HSUM / vstate.HSCAL;
            for (int iback = 1; iback <= j+1; ++iback) {
                int i = (j + 4) - iback - 1;
                vstate.EL[i] = vstate.EL[i] * XI + vstate.EL[i-1];
            }
        }

        // Subtract correction terms from YH array.
        for (int j = 3; j <= vstate.NQ; ++j) {
            for (int i = 0; i < VODE_NEQS; ++i) {
                vstate.YH[j-1][i] = vstate.YH[j-1][i] - vstate.YH[vstate.L-1][i] * vstate.EL[j-1];
            }
        }

    }
    else {

        // Order increase.
        for (int j = 1; j <= VODE_LMAX; ++j) {
            vstate.EL[j-1] = 0.0_rt;
        }

        vstate.EL[2] = 1.0_rt;
        Real ALPH0 = -1.0_rt;
        Real ALPH1 = 1.0_rt;
        Real PROD = 1.0_rt;
        Real XIOLD = 1.0_rt;
        Real HSUM = vstate.HSCAL;

        if (vstate.NQ != 1) {
            for (int j = 1; j <= NQM1; ++j) { 
                // Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)).
                HSUM += vstate.TAU[j];
                Real XI = HSUM / vstate.HSCAL;
                PROD *= XI;
                ALPH0 -= 1.0_rt / (j + 1);
                ALPH1 += 1.0_rt / XI;
                for (int iback = 1; iback <= j+1; ++iback) {
                    int i = (j + 4) - iback - 1;
                    vstate.EL[i] = vstate.EL[i] * XIOLD + vstate.EL[i-1];
                }
                XIOLD = XI;
            }
        }       

        Real T1 = (-ALPH0 - ALPH1) / PROD;
        // Load column L + 1 in YH array.
        for (int n = 0; n < VODE_NEQS; ++n) {
            vstate.YH[vstate.L][n] = T1 * vstate.YH[VODE_LMAX-1][n];
        }

        // Add correction terms to YH array.
        for (int j = 3; j <= vstate.NQ + 1; ++j) {
            for (int n = 0; n < VODE_NEQS; ++n) {
                vstate.YH[j-1][n] = vstate.YH[j-1][n] + vstate.EL[j-1] * vstate.YH[vstate.L][n];
            }
        }

    }

}

#endif
