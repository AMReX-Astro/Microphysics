#ifndef _vode_dvode_H_
#define _vode_dvode_H_

#include <vode_type.H>
#include <vode_dvhin.H>
#include <vode_dvstep.H>  
#include <vode_rhs.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dvode (burn_t& state, dvode_t& vstate)
{

    // Local variables
    Real H0, S;
    int IER, NITER;

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

    // Initial call to the RHS.

    Array1D<Real, 1, VODE_NEQS> f_init;

    rhs(state, vstate, f_init);

    for (int i = 1; i <= VODE_NEQS; ++i) {
        vstate.yh(i,2) = f_init(i);
    }

    vstate.NFE = 1;

    // Load the initial value array in yh.

    for (int i = 1; i <= VODE_NEQS; ++i) {
        vstate.yh(i,1) = vstate.y(i);
    }

    // Load and invert the ewt array. (H is temporarily set to 1.0.)
    vstate.NQ = 1;
    vstate.H = 1.0_rt;

    for (int i = 1; i <= VODE_NEQS; ++i) {
        vstate.ewt(i) = vstate.rtol(i) * std::abs(vstate.yh(i,1)) + vstate.atol(i);
        vstate.ewt(i) = 1.0_rt / vstate.ewt(i);
    }

    // Call DVHIN to set initial step size H0 to be attempted.
    H0 = 0.0_rt;
    dvhin(state, vstate, H0, NITER, IER);
    vstate.NFE += NITER;

    if (IER != 0) {
#ifndef AMREX_USE_GPU
        std::cout << "DVODE: TOUT too close to T to start integration" << std::endl;
#endif
        vstate.istate = -3;
        return;
    }

    // Load H with H0 and scale yh(:,2) by H0.
    vstate.H = H0;
    for (int i = 1; i <= VODE_NEQS; ++i) {
        vstate.yh(i,2) *= H0;
    }

    bool skip_loop_start = true;

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

               for (int i = 1; i <= VODE_NEQS; ++i) {
                   vstate.y(i) = vstate.yh(i,1);
               }

               vstate.t = vstate.tn;

               return;

           }

           for (int i = 1; i <= VODE_NEQS; ++i) {
               vstate.ewt(i) = vstate.rtol(i) * std::abs(vstate.yh(i,1)) + vstate.atol(i);
               vstate.ewt(i) = 1.0_rt / vstate.ewt(i);
           }

       }
       else {
           skip_loop_start = false;
       }

       Real TOLSF = 0.0_rt;
       for (int i = 1; i <= VODE_NEQS; ++i) {
           TOLSF += (vstate.yh(i,1) * vstate.yh(i,1)) * (vstate.ewt(i) * vstate.ewt(i));
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

           for (int i = 1; i <= VODE_NEQS; ++i) {
               vstate.y(i) = vstate.yh(i,1);
           }

           vstate.t = vstate.tn;

           return;

       }

       dvstep(pivot, state, vstate);

       // Branch on KFLAG. KFLAG can be 0, -1, or -2.

       if (vstate.kflag == -1) {
           // Error test failed repeatedly or with ABS(H) = HMIN.
#ifndef AMREX_USE_GPU
           std::cout << "DVODE: error test failed repeatedly or with abs(H) = HMIN" << std::endl;
#endif
           vstate.istate = -4;

           // Set Y array, T, and optional output.
           for (int i = 1; i <= VODE_NEQS; ++i) {
               vstate.y(i) = vstate.yh(i,1);
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
           for (int i = 1; i <= VODE_NEQS; ++i) {
               vstate.y(i) = vstate.yh(i,1);
           }

           vstate.t = vstate.tn;

           return;

       }

       // Otherwise, we've had a successful return from the integrator (KFLAG = 0).
       // Test for our stopping condition.

       if ((vstate.tn - vstate.tout) * vstate.H < 0.0_rt) continue;

       // If TOUT has been reached, interpolate.

       for (int i = 1; i <= VODE_NEQS; ++i) {
           vstate.y(i) = vstate.yh(i,vstate.L);
       }

       S = (vstate.tout - vstate.tn) / vstate.H;

       for (int jb = 1; jb <= vstate.NQ; ++jb) {
           int j = vstate.NQ - jb;
           for (int i = 1; i <= VODE_NEQS; ++i) {
               vstate.y(i) = vstate.yh(i,j+1) + S * vstate.y(i);
           }
       }

       vstate.t = vstate.tout;

       vstate.istate = 2;

       return;

    }

    // We've successfully completed the integration.

    for (int i = 1; i <= VODE_NEQS; ++i) {
        vstate.y(i) = vstate.yh(i,1);
    }

    vstate.t = vstate.tn;

    vstate.istate = 2;

}

#endif
