#ifndef VODE_DVODE_H
#define VODE_DVODE_H

#include <vode_type.H>
#include <vode_dvhin.H>
#include <vode_dvstep.H>
#ifdef STRANG
#include <integrator_rhs_strang.H>
#endif
#ifdef SDC
#include <integrator_rhs_simplified_sdc.H>
#endif
#ifdef NSE_TABLE
#include <nse_table_check.H>
#endif
#ifdef NSE_NET
#include <nse_check.H>
#endif

template <typename BurnT, typename DvodeT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int dvode (BurnT& state, DvodeT& vstate)
{
    constexpr int int_neqs = integrator_neqs<BurnT>();

    // Local variables
    Real H0{}, S{};
    int IER{}, NITER{};

    // Flag determining if we were successful.

    int istate = IERR_SUCCESS;

    // Return if the final time matches the starting time.

    if (vstate.tout == vstate.t) {
        return istate;
    }

    // All remaining initializations, the initial call to F,
    // and the calculation of the initial step size.
    // The error weights in ewt are inverted after being loaded.

    vstate.tn = vstate.t;

    vstate.NST = 0;
    vstate.NJE = 0;
    vstate.NSLJ = 0;

    // Initial call to the RHS.

    Array1D<Real, 1, int_neqs> f_init;

    rhs(vstate.t, state, vstate, f_init);

    for (int i = 1; i <= int_neqs; ++i) {
        vstate.yh(i,2) = f_init(i);
    }

    vstate.NFE = 1;

    // Load the initial value array in yh.

    for (int i = 1; i <= int_neqs; ++i) {
        vstate.yh(i,1) = vstate.y(i);
    }

    // Load and invert the ewt array. (H is temporarily set to 1.0.)
    vstate.NQ = 1;
    vstate.H = 1.0_rt;

    for (int i = 1; i <= NumSpec; ++i) {
        vstate.ewt(i) = vstate.rtol_spec * std::abs(vstate.yh(i,1)) + vstate.atol_spec;
        vstate.ewt(i) = 1.0_rt / vstate.ewt(i);
    }
    vstate.ewt(NumSpec+1) = vstate.rtol_enuc * std::abs(vstate.yh(NumSpec+1,1)) + vstate.atol_enuc;
    vstate.ewt(NumSpec+1) = 1.0_rt / vstate.ewt(NumSpec+1);

    // Call DVHIN to set initial step size H0 to be attempted.
    H0 = 0.0_rt;
    dvhin(state, vstate, H0, NITER, IER);
    vstate.NFE += NITER;

    if (IER != 0) {
#ifndef AMREX_USE_GPU
        std::cout << "DVODE: TOUT too close to T to start integration" << std::endl;
#endif
        istate = -3;
        return istate;
    }

    // Load H with H0 and scale yh(:,2) by H0.
    vstate.H = H0;
    for (int i = 1; i <= int_neqs; ++i) {
        vstate.yh(i,2) *= H0;
    }

    // Start with the order set to 1, and initialize other variables.
    // ETAMAX is the maximum ratio by which H can be increased
    // in a single step. It is normally 10, but is larger during the
    // first step to compensate for the small initial H. If a failure
    // occurs (in corrector convergence or error test), ETAMAX is set to 1
    // for the next increase.

    vstate.NQ = 1;
    vstate.NEWQ = 1;
    vstate.L = 2;
    vstate.tau(1) = vstate.H;
    vstate.PRL1 = 1.0_rt;
    vstate.RC = 0.0_rt;
    vstate.ETAMAX = 1.0e4_rt;
    vstate.NQWAIT = 2;
    vstate.HSCAL = vstate.H;
    vstate.NEWH = 0;
    vstate.NSLP = 0;
    vstate.IPUP = 1;

    bool skip_loop_start = true;

    // Now do the actual integration as a loop over dvstep.

    while (true) {

       if (!skip_loop_start) {

           // First check for too many steps being taken, update ewt (if not at
           // start of problem), check for too much accuracy being requested, and
           // check for H below the roundoff level in T.

           if (vstate.NST >= ode_max_steps) {
               // The maximum number of steps was taken before reaching TOUT.
#ifndef AMREX_USE_GPU
               std::cout << Font::Bold << FGColor::Red << "DVODE: maximum number of steps taken before reaching TOUT" << ResetDisplay << std::endl;
#endif
               for (int i = 1; i <= int_neqs; ++i) {
                   vstate.y(i) = vstate.yh(i,1);
               }

               vstate.t = vstate.tn;

               return IERR_TOO_MANY_STEPS;

           }

           for (int i = 1; i <= NumSpec; ++i) {
               vstate.ewt(i) = vstate.rtol_spec * std::abs(vstate.yh(i,1)) + vstate.atol_spec;
               vstate.ewt(i) = 1.0_rt / vstate.ewt(i);
           }
           vstate.ewt(NumSpec+1) = vstate.rtol_enuc * std::abs(vstate.yh(NumSpec+1,1)) + vstate.atol_enuc;
           vstate.ewt(NumSpec+1) = 1.0_rt / vstate.ewt(NumSpec+1);

       }
       else {
           skip_loop_start = false;
       }

       Real TOLSF = 0.0_rt;
       for (int i = 1; i <= int_neqs; ++i) {
           TOLSF += (vstate.yh(i,1) * vstate.ewt(i)) * (vstate.yh(i,1) * vstate.ewt(i));
       }
       TOLSF = UROUND * std::sqrt(TOLSF / int_neqs);

       if (TOLSF > 1.0_rt) {

           if (vstate.NST == 0) {
#ifndef AMREX_USE_GPU
               std::cout << Font::Bold << FGColor::Red << "DVODE: too much accuracy requested at start of integration" << ResetDisplay << std::endl;
#endif
               return IERR_TOO_MUCH_ACCURACY_REQUESTED;
           }

           // Too much accuracy requested for machine precision.
#ifndef AMREX_USE_GPU
           std::cout << "DVODE: too much accuracy requested" << std::endl;
#endif
           for (int i = 1; i <= int_neqs; ++i) {
               vstate.y(i) = vstate.yh(i,1);
           }

           vstate.t = vstate.tn;

           return IERR_TOO_MUCH_ACCURACY_REQUESTED;

       }

       int kflag = dvstep(state, vstate);


       // Branch on KFLAG. KFLAG can be 0, -1, or -2.

       if (kflag == -1) {
           // Error test failed repeatedly or with ABS(H) = HMIN.
#ifndef AMREX_USE_GPU
           std::cout << Font::Bold << FGColor::Red << "DVODE: error test failed repeatedly or with abs(H) = HMIN" << ResetDisplay << std::endl;
#endif
           // Set Y array, T, and optional output.
           for (int i = 1; i <= int_neqs; ++i) {
               vstate.y(i) = vstate.yh(i,1);
           }

           vstate.t = vstate.tn;
           return IERR_DT_UNDERFLOW;

       }
       else if (kflag == -2) {
           // Convergence failed repeatedly or with ABS(H) = HMIN.
#ifndef AMREX_USE_GPU
           std::cout << Font::Bold << FGColor::Red << "DVODE: corrector convergence failed repeatedly or with abs(H) = HMIN" << ResetDisplay << std::endl;
#endif
           // Set Y array, T, and optional output.
           for (int i = 1; i <= int_neqs; ++i) {
               vstate.y(i) = vstate.yh(i,1);
           }

           vstate.t = vstate.tn;
           return IERR_CORRECTOR_CONVERGENCE;

       }

#ifdef NSE
       // check if, during the course of integration, we hit NSE, and
       // if so, bail out we rely on the state being consistent after
       // the call to dvstep, even if the step failed.

       // we only do this after MIN_NSE_BAILOUT_STEPS to prevent us
       // from hitting this right at the start when VODE might do so
       // wild exploration.  Also ensure we are not working > tmax,
       // so we don't need to worry about extrapolating back in time.

       if (vstate.NST > MIN_NSE_BAILOUT_STEPS && vstate.tn <= vstate.tout) {
           // first we need to make the burn_t in sync

#ifdef STRANG
           update_thermodynamics(state, vstate);
#endif
#ifdef SDC
           int_to_burn(vstate.tn, vstate, state);
#endif

           if (in_nse(state)) {
               vstate.t = vstate.tn;
               return IERR_ENTERED_NSE;
           }
       }
#endif


       // Otherwise, we've had a successful return from the integrator (kflag = 0).
       // Test for our stopping condition.

       if ((vstate.tn - vstate.tout) * vstate.H < 0.0_rt) continue;

       // If TOUT has been reached, interpolate.

       for (int i = 1; i <= int_neqs; ++i) {
           vstate.y(i) = vstate.yh(i,vstate.L);
       }

       S = (vstate.tout - vstate.tn) / vstate.H;

       for (int jb = 1; jb <= vstate.NQ; ++jb) {
           const int j = vstate.NQ - jb;
           for (int i = 1; i <= int_neqs; ++i) {
               vstate.y(i) = vstate.yh(i,j+1) + S * vstate.y(i);
           }
       }

       vstate.t = vstate.tout;

       return IERR_SUCCESS;

    }
}

#endif
