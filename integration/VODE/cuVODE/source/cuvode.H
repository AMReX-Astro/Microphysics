#ifndef _cuvode_H_
#define _cuvode_H_

#include <cuvode_parameters.H>
#include <cuvode_type.H>
#include <cuvode_dvhin.H>
#include <cuvode_dvstep.H>  
#include <vode_rhs.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dvode (dvode_t& vstate)
{

    // Return if the final time matches the starting time.

    if (vstate.TOUT == vstate.T) {
        return;
    }

    // All remaining initializations, the initial call to F,
    // and the calculation of the initial step size.
    // The error weights in EWT are inverted after being loaded.

    vstate.TN = vstate.T;

    vstate.JSTART = 0;
    vstate.NST = 0;
    vstate.NJE = 0;
    vstate.NSLJ = 0;
    vstate.HU = 0.0_rt;

    // Initial call to F.

    f_rhs(vstate.T, vstate, vstate.yh[1]);

    // We start off with one RHS evaluation.

    vstate.NFE = 1;

    // Load the initial value array in YH.

    for (int n = 0; n < VODE_NEQS; ++n) {
        vstate.YH[0][n] = vstate.Y[n];
    }

    // Load and invert the EWT array. (H is temporarily set to 1.0.)
    vstate.NQ = 1;
    vstate.H = 1.0_rt;

    for (int n = 0; n < VODE_NEQS; ++n) {
        vstate.ewt[n] = vstate.RTOL[n] * std::abs(vstate.YH[0][n]) + vstate.ATOL(n);
        vstate.ewt[n] = 1.0_rt / vstate.ewt[n];
    }

    // Call DVHIN to set initial step size H0 to be attempted.
    int NITER, IER;
    Real H0 = 0.0_rt;
    dvhin(vstate, H0, NITER, IER);
    vstate.NFE += NITER;

    if (IER != 0) {
#ifndef AMREX_USE_GPU
        std::cout << "DVODE: TOUT too close to T to start integration" << std::endl;
#endif
        vstate.ISTATE = -3;
        return;
    }

    // Load H with H0 and scale YH[1] by H0.
    vstate.H = H0;
    for (int n = 0; n < VODE_NEQS; ++n) {
        vstate.YH[1][n] *= H0;
    }

    bool skip_loop_start = true;
    Real pivot[VODE_NEQS];

    // Now do the actual integration as a loop over dvstep.

    while (true) {

       if (!skip_loop_start) {

           // First check for too many steps being taken, update EWT (if not at
           // start of problem), check for too much accuracy being requested, and
           // check for H below the roundoff level in T.

           if (vstate.NST >= vstate.MXSTEP) {
               // The maximum number of steps was taken before reaching TOUT.
#ifndef AMREX_USE_GPU
               std::cout << "DVODE: maximum number of steps taken before reaching TOUT" << std::endl;
#endif
               vstate.ISTATE = -1;

               for (int n = 0; n < VODE_NEQS; ++n) {
                   vstate.Y[n] = vstate.YH[0][n];
               }

               vstate.T = vstate.TN;

               return;

           }

           for (int n = 0; n < VODE_NEQS; ++n) {
               vstate.ewt[n] = vstate.RTOL[n] * std::abs(vstate.YH[0][n]) + vstate.ATOL[n];
               vstate.ewt[n] = 1.0_rt / vstate.ewt[n];
           }

       }
       else {
           skip_loop_start = false;
       }

       Real TOLSF = 0.0_rt;
       for (int n = 0; n < VODE_NEQS; ++n) {
           TOLSF += (vstate.YH[0][n] * vstate.YH[0][n]) * (vstate.EWT[n] * vstate.EWT[n]);
       }
       TOLSF = UROUND * std::sqrt(TOLSF / VODE_NEQS);

       if (TOLSF > 1.0_rt) {

           TOLSF *= 2.0_rt;

           if (vstate.NST == 0) {
#ifndef AMREX_USE_GPU
               std::cout << "DVODE: too much accuracy requested at start of integration" << std::endl;
#endif
               vstate.ISTATE = -3;
               return;
           }

           // Too much accuracy requested for machine precision.
#ifndef AMREX_USE_GPU
           std::cout << "DVODE: too much accuracy requested" << std::endl;
#endif
           vstate.ISTATE = -2;

           for (int n = 0; n < VODE_NEQS; ++n) {
               vstate.Y[n] = vstate.YH[0][n];
           }

           vstate.T = vstate.TN;

           return;

       }

       dvstep(pivot, vstate);

       // Branch on KFLAG. KFLAG can be 0, -1, or -2.

       if (vstate.kflag == -1) {
           // Error test failed repeatedly or with ABS(H) = HMIN.
#ifndef AMREX_USE_GPU
           std::cout << "DVODE: error test failed repeatedly or with abs(H) = HMIN" << std::endl;
#endif
           vstate.ISTATE = -4;

           // Set Y array, T, and optional output.
           for (int n = 0; n < VODE_NEQS; ++n) {
               vstate.Y[n] = vstate.YH[0][n];
           }

           vstate.T = vstate.TN;

           return;

       }
       else if (vstate.kflag == -2) {
           // Convergence failed repeatedly or with ABS(H) = HMIN.
#ifndef AMREX_USE_GPU
           std::cout << "DVODE: corrector convergence failed repeatedly or with abs(H) = HMIN" << std::endl;
#endif
           vstate.ISTATE = -5;

           // Set Y array, T, and optional output.
           for (int n = 0; n < VODE_NEQS; ++n) {
               vstate.Y[n] = vstate.YH[0][n];
           }

           vstate.T = vstate.TN;

           return;

       }

       // Otherwise, we've had a successful return from the integrator (KFLAG = 0).
       // Test for our stopping condition.

       if ((vstate.TN - vstate.TOUT) * vstate.H < 0.0_rt) continue;

       // If TOUT has been reached, interpolate.

       for (int n = 0; n < VODE_NEQS; ++n) {
           vstate.Y[n] = vstate.YH[vstate.L-1][n];
       }

       Real S = (vstate.TOUT - vstate.TN) / vstate.H;

       for (int jb = 1; jb <= vstate.NQ; ++jb) {
           int j = vstate.NQ - jb - 1;
           for (int n = 0; n < VODE_NEQS; ++n) {
               vstate.Y[n] = vstate.YH[j+1][n] + S * vstate.Y[n];
           }
       }

       vstate.T = vstate.TOUT;

       vstate.ISTATE = 2;

       return;

    }

    // We've successfully completed the ingration.

    for (int n = 0; n < VODE_NEQS; ++n) {
        vstate.Y[n] = vstate.YH[0][n];
    }

    vstate.T = vstate.TN;

    vstate.ISTATE = 2;

}

#endif
