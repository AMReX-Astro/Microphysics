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

    // Local variables
    bool IHIT;
    Real EWTI, H0, HMAX, HMX, S;
    Real RH, SIZE, TCRIT, TNEXT, TOLSF, TP;
    int IER, KGO, LENJ, LENP;
    int MBAND, MFA, ML, MU, NITER;
    int NSLAST;

    Array1D<int, 1, VODE_NEQS> pivot;

    // Return if the final time matches the starting time.

    if (vstate.tout == vstate.t) {
        return;
    }

    // All remaining initializations, the initial call to F,
    // and the calculation of the initial step size.
    // The error weights in ewt are inverted after being loaded.

    vstate.tn = vstate.t;

    vstate.JSTART = 0;
    vstate.NST = 0;
    vstate.NJE = 0;
    vstate.NSLJ = 0;
    vstate.HU = 0.0_rt;

    // Initial call to F.

    f_rhs(vstate.t, vstate, vstate.yh[1]);

    // We start off with one RHS evaluation.

    vstate.NFE = 1;

    // Load the initial value array in yh.

    for (int n = 0; n < VODE_NEQS; ++n) {
        vstate.yh[0][n] = vstate.y[n];
    }

    // Load and invert the ewt array. (H is temporarily set to 1.0.)
    vstate.NQ = 1;
    vstate.H = 1.0_rt;

    for (int n = 0; n < VODE_NEQS; ++n) {
        vstate.ewt[n] = vstate.rtol[n] * std::abs(vstate.yh[0][n]) + vstate.atol[n];
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
        vstate.istate = -3;
        return;
    }

    // Load H with H0 and scale yh[1] by H0.
    vstate.H = H0;
    for (int n = 0; n < VODE_NEQS; ++n) {
        vstate.yh[1][n] *= H0;
    }

    bool skip_loop_start = true;
    int pivot[VODE_NEQS];

    // Now do the actual integration as a loop over dvstep.

    while (true) {

       if (!skip_loop_start) {

           // First check for too many steps being taken, update ewt (if not at
           // start of problem), check for too much accuracy being requested, and
           // check for H below the roundoff level in T.

           if (vstate.NST >= vstate.MXSTEP) {
               // The maximum number of steps was taken before reaching TOUT.
#ifndef AMREX_USE_GPU
               std::cout << "DVODE: maximum number of steps taken before reaching TOUT" << std::endl;
#endif
               vstate.istate = -1;

               for (int n = 0; n < VODE_NEQS; ++n) {
                   vstate.y[n] = vstate.yh[0][n];
               }

               vstate.t = vstate.tn;

               return;

           }

           for (int n = 0; n < VODE_NEQS; ++n) {
               vstate.ewt[n] = vstate.rtol[n] * std::abs(vstate.yh[0][n]) + vstate.atol[n];
               vstate.ewt[n] = 1.0_rt / vstate.ewt[n];
           }

       }
       else {
           skip_loop_start = false;
       }

       Real TOLSF = 0.0_rt;
       for (int n = 0; n < VODE_NEQS; ++n) {
           TOLSF += (vstate.yh[0][n] * vstate.yh[0][n]) * (vstate.ewt[n] * vstate.ewt[n]);
       }
       TOLSF = UROUND * std::sqrt(TOLSF / VODE_NEQS);

       if (TOLSF > 1.0_rt) {

           TOLSF *= 2.0_rt;

           if (vstate.NST == 0) {
#ifndef AMREX_USE_GPU
               std::cout << "DVODE: too much accuracy requested at start of integration" << std::endl;
#endif
               vstate.istate = -3;
               return;
           }

           // Too much accuracy requested for machine precision.
#ifndef AMREX_USE_GPU
           std::cout << "DVODE: too much accuracy requested" << std::endl;
#endif
           vstate.istate = -2;

           for (int n = 0; n < VODE_NEQS; ++n) {
               vstate.y[n] = vstate.yh[0][n];
           }

           vstate.t = vstate.tn;

           return;

       }

       dvstep(pivot, vstate);

       // Branch on KFLAG. KFLAG can be 0, -1, or -2.

       if (vstate.kflag == -1) {
           // Error test failed repeatedly or with ABS(H) = HMIN.
#ifndef AMREX_USE_GPU
           std::cout << "DVODE: error test failed repeatedly or with abs(H) = HMIN" << std::endl;
#endif
           vstate.istate = -4;

           // Set Y array, T, and optional output.
           for (int n = 0; n < VODE_NEQS; ++n) {
               vstate.y[n] = vstate.yh[0][n];
           }

           vstate.t = vstate.tn;

           return;

       }
       else if (vstate.kflag == -2) {
           // Convergence failed repeatedly or with ABS(H) = HMIN.
#ifndef AMREX_USE_GPU
           std::cout << "DVODE: corrector convergence failed repeatedly or with abs(H) = HMIN" << std::endl;
#endif
           vstate.istate = -5;

           // Set Y array, T, and optional output.
           for (int n = 0; n < VODE_NEQS; ++n) {
               vstate.y[n] = vstate.yh[0][n];
           }

           vstate.t = vstate.tn;

           return;

       }

       // Otherwise, we've had a successful return from the integrator (KFLAG = 0).
       // Test for our stopping condition.

       if ((vstate.tn - vstate.tout) * vstate.H < 0.0_rt) continue;

       // If TOUT has been reached, interpolate.

       for (int n = 0; n < VODE_NEQS; ++n) {
           vstate.y[n] = vstate.yh[vstate.L-1][n];
       }

       Real S = (vstate.tout - vstate.tn) / vstate.H;

       for (int jb = 1; jb <= vstate.NQ; ++jb) {
           int j = vstate.NQ - jb - 1;
           for (int n = 0; n < VODE_NEQS; ++n) {
               vstate.y[n] = vstate.yh[j+1][n] + S * vstate.y[n];
           }
       }

       vstate.t = vstate.tout;

       vstate.istate = 2;

       return;

    }

    // We've successfully completed the integration.

    for (int n = 0; n < VODE_NEQS; ++n) {
        vstate.y[n] = vstate.yh[0][n];
    }

    vstate.t = vstate.tn;

    vstate.istate = 2;

}

#endif
