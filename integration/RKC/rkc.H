#ifndef RKC_H
#define RKC_H

#include <rkc_type.H>
#include <rkc_util.H>
#include <burn_type.H>
#ifdef STRANG
#include <integrator_type_strang.H>
#include <integrator_rhs_strang.H>
#endif
#ifdef SDC
#include <integrator_type_sdc.H>
#include <integrator_rhs_sdc.H>
#endif
#include <circle_theorem.H>
#include <integrator_data.H>

#ifdef NSE_TABLE
#include <nse_table_check.H>
#endif
#ifdef NSE_NET
#include <nse_check.H>
#endif


template <typename BurnT, typename RkcT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void step (BurnT& state, RkcT& rstate, const amrex::Real h, const int m)
{
    constexpr int int_neqs = integrator_neqs<BurnT>();

    // Take a step of size H from T to T+H to get Y(*).

    amrex::Real w0 = 1.0_rt + 2.0_rt / (13.0_rt * m * m);

    amrex::Real temp1 = w0 * w0 - 1.0_rt;
    amrex::Real temp2 = std::sqrt(temp1);
    amrex::Real arg = m * std::log(w0 + temp2);
    amrex::Real w1 = std::sinh(arg) * temp1 /
        (std::cosh(arg) * m * temp2 -  w0 * std::sinh(arg));
    amrex::Real bjm1 = 1.0_rt / amrex::Math::powi<2>(2.0_rt * w0);
    amrex::Real bjm2 = bjm1;

    // Evaluate the first stage.

    for (int i = 1; i <= int_neqs; ++i) {
        rstate.yjm2(i) = rstate.yn(i);
    }
    amrex::Real mus = w1 * bjm1;
    for (int i = 1; i <= int_neqs; ++i) {
        rstate.yjm1(i) = rstate.yn(i) + h * mus * rstate.fn(i);
    }

    amrex::Real thjm2{};
    amrex::Real thjm1 = mus;
    amrex::Real zjm1 = w0;
    amrex::Real zjm2 = 1.0_rt;
    amrex::Real dzjm1 = 1.0_rt;
    amrex::Real dzjm2 = 0.0_rt;
    amrex::Real d2zjm1 = 0.0_rt;
    amrex::Real d2zjm2 = 0.0_rt;

    // Evaluate stages j = 2,...,m.

    for (int j = 2; j <= m; ++j) {
        amrex::Real zj = 2.0_rt * w0 * zjm1 - zjm2;
        amrex::Real dzj = 2.0_rt * w0 * dzjm1 - dzjm2 + 2.0_rt * zjm1;
        amrex::Real d2zj = 2.0_rt * w0 * d2zjm1 - d2zjm2 + 4.0_rt * dzjm1;
        amrex::Real bj = d2zj / amrex::Math::powi<2>(dzj);
        amrex::Real ajm1 = 1.0_rt - zjm1 * bjm1;
        amrex::Real mu = 2.0_rt * w0 * bj / bjm1;
        amrex::Real nu = -bj / bjm2;
        mus = mu * w1 / w0;

        // Use the y array for temporary storage here.

        // for this call
        // we want rstate.yjm1() to be the input y
        // the original Fortran code stored the output temporarily
        // in rstate.y(), but we need that for input

        for (int i = 1; i <= int_neqs; ++i) {
            rstate.y(i) = rstate.yjm1(i);
        }
        RArray1D ydot;
        rhs(rstate.t + h*thjm1, state, rstate, ydot);

        for (int i = 1; i <= int_neqs; ++i) {
            rstate.y(i) = mu * rstate.yjm1(i) +
                          nu * rstate.yjm2(i) +
                          (1.0_rt - mu - nu) * rstate.yn(i) +
                          h * mus * (ydot(i) - ajm1 * rstate.fn(i));
        }
        amrex::Real thj = mu * thjm1 + nu * thjm2 + mus * (1.0_rt - ajm1);

        // Shift the data for the next stage.

        if (j < m) {
            for (int i = 1; i <= int_neqs; ++i) {
                rstate.yjm2(i) = rstate.yjm1(i);
                rstate.yjm1(i) = rstate.y(i);
            }
            thjm2 = thjm1;
            thjm1 = thj;
            bjm2 = bjm1;
            bjm1 = bj;
            zjm2 = zjm1;
            zjm1 = zj;
            dzjm2 = dzjm1;
            dzjm1 = dzj;
            d2zjm2 = d2zjm1;
            d2zjm1 = d2zj;
        }
    }
}


template <typename BurnT, typename RkcT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rkclow (BurnT& state, RkcT& rstate)
{
    constexpr int int_neqs = integrator_neqs<BurnT>();

    // RKC is an interface to RKCLOW where the actual solution takes place.

    const amrex::Real onep1 = 1.1_rt;
    const amrex::Real p4 = 0.4_rt;
    const amrex::Real p8 = 0.8_rt;
    const amrex::Real one3rd = 1.0_rt/3.0_rt;
    const amrex::Real two3rd = 2.0_rt/3.0_rt;

    // Initialize on the first call.

    // note: the original Fortran code only had a single rtol, so
    // I needed to generalize this
    int mmax = static_cast<int>(std::round(std::sqrt(std::max(rstate.rtol_spec,
                                                              rstate.rtol_enuc) / (10.0_rt * UROUND))));
    mmax = std::max(mmax, 2);
    bool newspc = true;
    bool jacatt = false;
    int nstsig = 0;
    for (int i = 1; i <= int_neqs; ++i) {
        rstate.yn(i) = rstate.y(i);
    }

    // we want to call with yn = y as the input and store the output in fn
    rhs(rstate.t, state, rstate, rstate.fn);

    rstate.n_rhs++;
    amrex::Real tdir = std::copysign(1.0_rt, rstate.tout - rstate.t);
    rstate.hmax = std::abs(rstate.tout - rstate.t);

    amrex::Real hmin{};
    amrex::Real sprad{};
    amrex::Real absh{};
    amrex::Real errold{};
    amrex::Real h{};
    amrex::Real hold{};

    // Start of loop for taking one step.
    while (rstate.n_step < integrator_rp::ode_max_steps) {

        // Estimate the spectral radius of the Jacobian
        // when newspc = .true..

        if (newspc) {
            if (integrator_rp::use_circle_theorem) {
                circle_theorem_sprad(rstate.t, state, rstate, sprad);
            } else {
                int ierr = rkcrho(state, rstate, rstate.hmax, sprad);
                if (ierr !=0) {
                    return IERR_SPRAD_CONVERGENCE;
                }
            }
            jacatt = true;
        }


        // Compute an initial step size.

        if (rstate.n_step == 0) {
            absh = rkc_init_dt(state, rstate, rstate.hmax, sprad);
        }

        // Adjust the step size and determine the number of stages m.

        bool last = false;
        if (onep1 * absh >=  std::abs(rstate.tout - rstate.t)) {
            absh = std::abs(rstate.tout - rstate.t);
            last = true;
        }
        int m = 1 + static_cast<int>(std::sqrt(1.54_rt * absh * sprad + 1.0_rt));

        // Limit m to mmax to control the growth of roundoff error.

        if (m > mmax) {
            m = mmax;
            absh = static_cast<amrex::Real>(m * m - 1) / (1.54_rt * sprad);
            last = false;
        }
        rstate.maxm = std::max(m, rstate.maxm);

        // A tentative solution at t+h is returned in
        // y and its slope is evaluated in yjm1(*).

        h = tdir * absh;
        hmin = 10.0_rt * UROUND * std::max(std::abs(rstate.t),
                                           std::abs(rstate.t + h));
        step(state, rstate, h, m);
        rhs(rstate.t+h, state, rstate, rstate.yjm1);
        rstate.n_rhs += m;
        rstate.n_step++;

        // Estimate the local error and compute its weighted RMS norm.

        amrex::Real err{};
        for (int i = 1; i <= int_neqs; ++i) {
            amrex::Real wt{};
            if (i <= NumSpec) {
                wt = rstate.rtol_spec * std::max(std::abs(rstate.yn(i)),
                                                 std::abs(rstate.y(i))) +
                    rstate.atol_spec;
            } else {
                wt = rstate.rtol_enuc * std::max(std::abs(rstate.yn(i)),
                                                 std::abs(rstate.y(i))) +
                    rstate.atol_enuc;
            }

            amrex::Real est = p8 * (rstate.yn(i) - rstate.y(i)) +
                       p4 * h * (rstate.fn(i) + rstate.yjm1(i));
            err += amrex::Math::powi<2>(est / wt);
        }

        err = std::sqrt(err / int_neqs);

        // before we accept or reject the step, let's check if we've entered
        // NSE

#ifdef NSE
       // check if, during the course of integration, we hit NSE, and
       // if so, bail out we rely on the state being consistent after
       // the call to dvstep, even if the step failed.

       // we only do this after MIN_NSE_BAILOUT_STEPS to prevent us
       // from hitting this right at the start when VODE might do so
       // wild exploration.  Also ensure we are not working > tmax,
       // so we don't need to worry about extrapolating back in time.

       if (rstate.nsteps > MIN_NSE_BAILOUT_STEPS && rstate.t <= rstate.tout) {
           // first we need to make the burn_t in sync

#ifdef STRANG
           update_thermodynamics(state, rstate);
#endif
#ifdef SDC
           int_to_burn(rstate.t, rstate, state);
#endif

           if (in_nse(state)) {
               return IERR_ENTERED_NSE;
           }
       }
#endif

        if (err > 1.0_rt) {
            // Step is rejected.
            rstate.nrejct++;
            absh = p8 * absh / std::pow(err, one3rd);
            if (absh < hmin) {
                return IERR_DT_UNDERFLOW;
            } else {
                newspc = ! jacatt;
                continue;
            }
        }

        // Step is accepted.

        rstate.naccpt++;
        rstate.t += h;
        jacatt = false;
        nstsig = nstsig+1 % 25;
        newspc = false;
        if (nstsig == 0) {
            newspc = ! jacatt;
        }

        // Update the data for interpolation stored in work(*).
        for (int i = 1; i <= int_neqs; ++i) {
            amrex::Real ylast = rstate.yn(i);
            amrex::Real yplast = rstate.fn(i);
            rstate.yn(i) = rstate.y(i);
            rstate.fn(i) = rstate.yjm1(i);
            rstate.yjm1(i) = ylast;
            rstate.yjm2(i) = yplast;
        }
        amrex::Real fac = 10.0_rt;
        if (rstate.naccpt == 1) {
            amrex::Real temp2 = std::pow(err, one3rd);
            if (p8 < fac * temp2) {
                fac = p8/temp2;
            }
        } else {
            amrex::Real temp1 = p8 * absh * std::pow(errold, one3rd);
            amrex::Real temp2 = std::abs(hold) * std::pow(err, two3rd);
            if (temp1 < fac * temp2) {
                fac = temp1 / temp2;
            }
        }
        absh = std::max(0.1_rt, fac) * absh;
        absh = std::max(hmin, std::min(rstate.hmax, absh));
        errold = err;
        hold = h;
        h = tdir * absh;

        if (last) {
            return IERR_SUCCESS;
        }
    }

    return IERR_TOO_MANY_STEPS;
}

template <typename BurnT, typename RkcT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rkc (BurnT& state, RkcT& rstate)
{
    const amrex::Real rmax = 0.1_rt;
    const amrex::Real rmin = 10.0_rt * UROUND;

    // Test the input data.
    bool valid = true;

    if ((rtol_spec > rmax) || (rtol_enuc > rmax) ||
        (rtol_spec < rmin) || (rtol_enuc < rmin)) {
        valid = false;
    }

    if (atol_spec < 0.0_rt || atol_enuc < 0.0_rt) {
        valid = false;
    }

    if (! valid) {
        return IERR_BAD_INPUTS;
    }

    // Initialize counters and pointers.

    rstate.n_rhs = 0;
    rstate.n_step = 0;
    rstate.naccpt = 0;
    rstate.nrejct = 0;
    rstate.nfesig = 0;
    rstate.maxm = 0;

    int ierr = rkclow(state, rstate);
    return ierr;
}

#endif
