#ifndef RKC_H
#define RKC_H

#include <rkc_type.H>
#include <rkc_util.H>
#include <burn_type.H>
#ifdef STRANG
#include <integrator_type_strang.H>
#include <integrator_rhs_strang.H>
#endif
#ifdef SIMPLIFIED_SDC
#include <integrator_type_simplified_sdc.H>
#include <integrator_rhs_simplified_sdc.H>
#endif

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void step (BurnT& state, rkc_t& rstate, const Real h, const int m)
{

    // Take a step of size H from T to T+H to get Y(*).

    Real w0 = 1.0_rt + 2.0_rt / (13.0_rt * m * m);

    Real temp1 = w0 * w0 - 1.0_rt;
    Real temp2 = std::sqrt(temp1);
    Real arg = m * std::log(w0 + temp2);
    Real w1 = std::sinh(arg) * temp1 /
        (std::cosh(arg) * m * temp2 -  w0 * std::sinh(arg));
    Real bjm1 = 1.0_rt / std::pow(2.0_rt * w0, 2);
    Real bjm2 = bjm1;

    // Evaluate the first stage.

    for (int i = 1; i <= RKC_NEQS; ++i) {
        rstate.yjm2(i) = rstate.yn(i);
    }
    Real mus = w1 * bjm1;
    for (int i = 1; i <= RKC_NEQS; ++i) {
        rstate.yjm1(i) = rstate.yn(i) + h * mus * rstate.fn(i);
    }

    Real thjm2{};
    Real thjm1 = mus;
    Real zjm1 = w0;
    Real zjm2 = 1.0_rt;
    Real dzjm1 = 1.0_rt;
    Real dzjm2 = 0.0_rt;
    Real d2zjm1 = 0.0_rt;
    Real d2zjm2 = 0.0_rt;

    // Evaluate stages j = 2,...,m.

    for (int j = 2; j <= m; ++j) {
        Real zj = 2.0_rt * w0 * zjm1 - zjm2;
        Real dzj = 2.0_rt * w0 * dzjm1 - dzjm2 + 2.0_rt * zjm1;
        Real d2zj = 2.0_rt * w0 * d2zjm1 - d2zjm2 + 4.0_rt * dzjm1;
        Real bj = d2zj / std::pow(dzj, 2);
        Real ajm1 = 1.0_rt - zjm1 * bjm1;
        Real mu = 2.0_rt * w0 * bj / bjm1;
        Real nu = -bj / bjm2;
        mus = mu * w1 / w0;

        // Use the y array for temporary storage here.

        // for this call
        // we want rstate.yjm1() to be the input y
        // the original Fortran code stored the output temporarily
        // in rstate.y(), but we need that for input

        for (int i = 1; i <= RKC_NEQS; ++i) {
            rstate.y(i) = rstate.yjm1(i);
        }
        RArray1D ydot;
        rhs(rstate.t + h*thjm1, state, rstate, ydot);

        for (int i = 1; i <= RKC_NEQS; ++i) {
            rstate.y(i) = mu * rstate.yjm1(i) +
                          nu * rstate.yjm2(i) +
                          (1.0_rt - mu - nu) * rstate.yn(i) +
                          h * mus * (ydot(i) - ajm1 * rstate.fn(i));
        }
        Real thj = mu * thjm1 + nu * thjm2 + mus * (1.0_rt - ajm1);

        // Shift the data for the next stage.

        if (j < m) {
            for (int i = 1; i <= RKC_NEQS; ++i) {
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


template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rkclow (BurnT& state, rkc_t& rstate)
{


    // RKC is an interface to RKCLOW where the actual solution takes place.

    const Real onep1 = 1.1_rt;
    const Real p4 = 0.4_rt;
    const Real p8 = 0.8_rt;
    const Real one3rd = 1.0_rt/3.0_rt;
    const Real two3rd = 2.0_rt/3.0_rt;

    // Initialize on the first call.

    // note: the original Fortran code only had a single rtol, so
    // I needed to generalize this
    int mmax = std::round(std::sqrt(std::max(rstate.rtol_spec, rstate.rtol_enuc) /
                                    (10.0_rt * UROUND)));
    mmax = std::max(mmax, 2);
    bool newspc = true;
    bool jacatt = false;
    int nstsig = 0;
    for (int i = 1; i <= RKC_NEQS; ++i) {
        rstate.yn(i) = rstate.y(i);
    }

    // we want to call with yn = y as the input and store the output in fn
    rhs(rstate.t, state, rstate, rstate.fn);

    rstate.nfe++;
    Real tdir = std::copysign(1.0_rt, rstate.tend - rstate.t);
    rstate.hmax = std::abs(rstate.tend - rstate.t);

    Real hmin{};
    Real sprad{};
    Real absh{};
    Real errold{};
    Real h{};
    Real hold{};

    // Start of loop for taking one step.
    while (true) {

        // Estimate the spectral radius of the Jacobian
        // when newspc = .true..  A convergence failure
        // in rkcrho is reported by idid = 6.

        if (newspc) {
            sprad = rkcrho(state, rstate, rstate.hmax);
            if (rstate.idid == 6) {
                return rstate.idid;
            }

            jacatt = true;
        }


        // Compute an initial step size.

        if (rstate.nsteps == 0) {
            absh = rkc_init_dt(state, rstate, rstate.hmax, sprad);
        }

        // Adjust the step size and determine the number of stages m.

        bool last = false;
        if (onep1 * absh >=  std::abs(rstate.tend - rstate.t)) {
            absh = std::abs(rstate.tend - rstate.t);
            last = true;
        }
        int m = 1 + static_cast<int>(std::sqrt(1.54_rt * absh * sprad + 1.0_rt));

        // Limit m to mmax to control the growth of roundoff error.

        if (m > mmax) {
            m = mmax;
            absh = static_cast<Real>(m * m - 1) / (1.54_rt * sprad);
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
        rstate.nfe += m;
        rstate.nsteps++;

        // Estimate the local error and compute its weighted RMS norm.

        Real err{};
        for (int i = 1; i <= RKC_NEQS; ++i) {
            Real wt{};
            if (i <= NumSpec) {
                wt = rstate.rtol_spec * std::max(std::abs(rstate.yn(i)),
                                                 std::abs(rstate.y(i))) +
                    rstate.atol_spec;
            } else {
                wt = rstate.rtol_enuc * std::max(std::abs(rstate.yn(i)),
                                                 std::abs(rstate.y(i))) +
                    rstate.atol_enuc;
            }
            if (wt == 0.0_rt) {
                rstate.idid = 3;
                return rstate.idid;
            }
            Real est = p8 * (rstate.yn(i) - rstate.y(i)) +
                       p4 * h * (rstate.fn(i) + rstate.yjm1(i));
            err += std::pow(est / wt, 2);
        }

        err = std::sqrt(err / RKC_NEQS);

        if (err > 1.0_rt) {
            // Step is rejected.
            rstate.nrejct++;
            absh = p8 * absh / std::pow(err, one3rd);
            if (absh < hmin) {
                rstate.idid = 4;
                return rstate.idid;
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
        for (int i = 1; i <= RKC_NEQS; ++i) {
            Real ylast = rstate.yn(i);
            Real yplast = rstate.fn(i);
            rstate.yn(i) = rstate.y(i);
            rstate.fn(i) = rstate.yjm1(i);
            rstate.yjm1(i) = ylast;
            rstate.yjm2(i) = yplast;
        }
        Real fac = 10.0_rt;
        if (rstate.naccpt == 1) {
            Real temp2 = std::pow(err, one3rd);
            if (p8 < fac * temp2) {
                fac = p8/temp2;
            }
        } else {
            Real temp1 = p8 * absh * std::pow(errold, one3rd);
            Real temp2 = std::abs(hold) * std::pow(err, two3rd);
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
            rstate.idid = 1;
            return rstate.idid;
        }
    }
}

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rkc (BurnT& state, rkc_t& rstate)
{

    const Real rmax = 0.1_rt;
    const Real rmin = 10.0_rt * UROUND;

    // the original rkc.f checked on the input value of idid, but we
    // only work for idid = 0, which is a new integration (not a
    // continuation)

    // Test the input data.
    bool valid = true;

    if ((rtol_spec > rmax) || (rtol_enuc > rmax) ||
        (rtol_spec < rmin) || (rtol_spec < rmin)) {
        valid = false;
    }

    if (atol_spec < 0.0_rt || atol_enuc < 0.0_rt) {
        valid = false;
    }

    if (! valid) {
        rstate.idid = 5;  // TODO: remove idid from rstate
        return rstate.idid;
    }

    // Initialize counters and pointers.

    rstate.nfe = 0;
    rstate.nsteps = 0;
    rstate.naccpt = 0;
    rstate.nrejct = 0;
    rstate.nfesig = 0;
    rstate.maxm = 0;

    rkclow(state, rstate);
    return rstate.idid;
}

#endif
