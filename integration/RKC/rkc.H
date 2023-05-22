#ifndef RKC_H
#define RKC_H

#include <rkc_type.H>
#include <burn_type.H>
#ifdef STRANG
#include <integrator_type_strang.H>
#include <integrator_rhs_strang.H>
#endif
#ifdef SIMPLIFIED_SDC
#include <integrator_type_simplified_sdc.H>
#include <integrator_rhs_simplified_sdc.H>
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void step (burn_t& state, rkc_t& rstate, const Real h, const int m)
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
        rstate.vtemp2(i) = rstate.yn(i);
    }
    Real mus = w1 * bjm1;
    for (int i = 1; i <= RKC_NEQS; ++i) {
        rstate.vtemp1(i) = rstate.yn(i) + h * mus * rstate.fn(i);
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
        // we want rstate.vtemp1() to be the input y
        // the original Fortran code stored the output temporarily
        // in rstate.y(), but we need that for input

        for (int i = 1; i <= RKC_NEQS; ++i) {
            rstate.y(i) = rstate.vtemp1(i);
        }
        RArray1D ydot;
        rhs(rstate.t + h*thjm1, state, rstate, ydot);

        for (int i = 1; i <= RKC_NEQS; ++i) {
            rstate.y(i) = mu * rstate.vtemp1(i) +
                          nu * rstate.vtemp2(i) +
                          (1.0_rt - mu - nu) * rstate.yn(i) +
                          h * mus * (ydot(i) - ajm1 * rstate.fn(i));
        }
        Real thj = mu * thjm1 + nu * thjm2 + mus * (1.0_rt - ajm1);

        // Shift the data for the next stage.

        if (j < m) {
            for (int i = 1; i <= RKC_NEQS; ++i) {
                rstate.vtemp2(i) = rstate.vtemp1(i);
                rstate.vtemp1(i) = rstate.y(i);
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real rkcrho (burn_t& state, rkc_t& rstate)
{

    // RKCRHO attempts to compute a close upper bound, SPRAD, on
    // the spectral radius of the Jacobian matrix using a nonlinear
    // power method.  A convergence failure is reported by IDID = 6.

    constexpr int itmax{50};
    constexpr Real onep2{1.2_rt};
    constexpr Real p01{0.01_rt};

    // sprad smaller than small = 1/hmax are not
    // interesting because they do not constrain the step size.

    Real small = 1.0_rt / rstate.hmax;

    // The initial slope is used as guess when nsteps = 0 and
    // thereafter the last computed eigenvector.  Some care
    // is needed to deal with special cases. Approximations to
    // the eigenvector are normalized so that their Euclidean
    // norm has the constant value dynrm.

    if (rstate.nsteps == 0) {
        for (int i = 1; i <= RKC_NEQS; ++i) {
            rstate.vtemp1(i) = rstate.fn(i);
        }
    } else {
        for (int i = 1; i <= RKC_NEQS; ++i) {
            rstate.vtemp1(i) = rstate.sprad(i);
        }
    }

    Real ynrm{};
    Real vnrm{};
    for (int i = 1; i <= RKC_NEQS; ++i) {
        ynrm += std::pow(rstate.yn(i), 2);
        vnrm += std::pow(rstate.vtemp1(i), 2);
    }
    ynrm = std::sqrt(ynrm);
    vnrm = std::sqrt(vnrm);

    Real dynrm{};
    if (ynrm != 0.0_rt && vnrm != 0.0_rt) {
        dynrm = ynrm * std::sqrt(UROUND);
        for (int i = 1; i <= RKC_NEQS; ++i) {
            rstate.vtemp1(i) = rstate.yn(i) + rstate.vtemp1(i) * (dynrm / vnrm);
        }
    } else if (ynrm != 0.0_rt) {
        dynrm = ynrm * std::sqrt(UROUND);
        for (int i = 1; i <= RKC_NEQS; ++i) {
            rstate.vtemp1(i) = rstate.yn(i) + rstate.yn(i) * std::sqrt(UROUND);
        }
    } else if (vnrm != 0.0_rt) {
        dynrm = UROUND;
        for (int i = 1; i <= RKC_NEQS; ++i) {
            rstate.vtemp1(i) = rstate.vtemp1(i) * (dynrm / vnrm);
        }
    } else {
        dynrm = UROUND;
        for (int i = 1; i <= RKC_NEQS; ++i) {
            rstate.vtemp1(i) = dynrm;
        }
    }

    // Now iterate with a nonlinear power method.

    Real sigma{};
    Real sprad{};

    for (int iter = 0; iter < itmax; ++iter) {
        std::cout << "iter = " << iter << " " << sprad << std::endl;

        // we want to call the RHS with vtemp1() as the input, so we need to copy
        // this temporarily into vstate.y()
        RArray1D ysav;
        for (int i = 1; i <= RKC_NEQS; ++i) {
            ysav(i) = rstate.y(i);
            rstate.y(i) = rstate.vtemp1(i);
        }
        rhs(rstate.t, state, rstate, rstate.vtemp2);
        for (int i = 1; i <= RKC_NEQS; ++i) {
            rstate.y(i) = ysav(i);
        }
        rstate.nfesig++;
        Real dfnrm{};
        for (int i = 1; i <= RKC_NEQS; ++i) {
            dfnrm += std::pow(rstate.vtemp2(i) - rstate.fn(i), 2);
        }
        dfnrm = std::sqrt(dfnrm);
        Real sigmal = sigma;
        sigma = dfnrm / dynrm;

        // sprad is a little bigger than the estimate sigma of the
        // spectral radius, so is more likely to be an upper bound.

        sprad = onep2 * sigma;

        if (iter >= 1 &&
            std::abs(sigma - sigmal) <= std::max(sigma, small) * p01) {
            for (int i = 1; i <= RKC_NEQS; ++i) {
                rstate.sprad(i) = rstate.vtemp1(i) - rstate.yn(i);
            }
            return sprad;
        }

        // The next vtemp1(*) is the change in f
        // scaled so that norm(v - yn) = dynrm.

        if (dfnrm != 0.0_rt) {
            for (int i = 1; i <= RKC_NEQS; ++i) {
                rstate.vtemp1(i) = rstate.yn(i) +
                    (rstate.vtemp2(i) - rstate.fn(i)) * (dynrm / dfnrm);
            }
        } else {

            // The new v(*) degenerated to yn(*)--"randomly" perturb
            // current approximation to the eigenvector by changing
            // the sign of one component.

            int index = 1 + (iter + 1) % RKC_NEQS;
            rstate.vtemp1(index) = rstate.yn(index) -
                (rstate.vtemp1(index) - rstate.yn(index));
        }
    }

    // Set flag to report a convergence failure.

    rstate.idid = 6;
    return sprad;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rkclow (burn_t& state, rkc_t& rstate)
{


    // RKC is an interface to RKCLOW where the actual solution takes place.

    const Real onep1 = 1.1_rt;
    const Real p1 = 0.1_rt;
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

    Real hmin = 10.0_rt * UROUND * std::max(std::abs(rstate.t), rstate.hmax);
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
            sprad = rkcrho(state, rstate);
            if (rstate.idid == 6) {
                return rstate.idid;
            }

            jacatt = true;
        }


        // Compute an initial step size.

        if (rstate.nsteps == 0) {
            absh = rstate.hmax;
            if (sprad * absh > 1.0_rt) {
                absh = 1.0_rt / sprad;
            }
            absh = std::max(absh, hmin);

            for (int i = 1; i <= RKC_NEQS; ++i) {
                rstate.vtemp1(i) = rstate.yn(i) + absh * rstate.fn(i);
            }

            // we want to call with vtemp1 as the input, so we need to stash the
            // current y
            RArray1D ysav;
            for (int i = 1; i <= RKC_NEQS; ++i) {
                ysav(i) = rstate.y(i);
                rstate.y(i) = rstate.vtemp1(i);
            }
            rhs(rstate.t+absh, state, rstate, rstate.vtemp2);
            for (int i = 1; i <= RKC_NEQS; ++i) {
                rstate.y(i) = ysav(i);
            }

            rstate.nfe++;
            Real est{};

            // compute the weights using the tolerances
            for (int i = 1; i <= RKC_NEQS; ++i) {
                Real wt{};
                if (i <= NumSpec) {
                    wt = rstate.rtol_spec * std::abs(rstate.yn(i)) + rstate.atol_spec;
                } else {
                    wt = rstate.rtol_enuc * std::abs(rstate.yn(i)) + rstate.atol_enuc;
                }
                if (wt == 0.0_rt) {
                    rstate.idid = 3;
                    return rstate.idid;
                }
                est += std::pow((rstate.vtemp2(i) - rstate.fn(i)) / wt, 2);
            }
            est = absh * std::sqrt(est / RKC_NEQS);

            if (p1 * absh < rstate.hmax * std::sqrt(est)) {
                absh = std::max(p1 * absh / std::sqrt(est), hmin);
            } else {
                absh = rstate.hmax;
            }
        }

        // Adjust the step size and determine the number of stages m.

        bool last = false;
        if (onep1 * absh >=  std::abs(rstate.tend - rstate.t)) {
            absh = std::abs(rstate.tend - rstate.t);
            std::cout << "made it into last!!!!!!" << std::endl;
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
        // y and its slope is evaluated in vtemp1(*).

        h = tdir * absh;
        hmin = 10.0_rt * UROUND * std::max(std::abs(rstate.t),
                                           std::abs(rstate.t + h));
        step(state, rstate, h, m);
        rhs(rstate.t+h, state, rstate, rstate.vtemp1);
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
                       p4 * h * (rstate.fn(i) + rstate.vtemp1(i));
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
            rstate.fn(i) = rstate.vtemp1(i);
            rstate.vtemp1(i) = ylast;
            rstate.vtemp2(i) = yplast;
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
        absh = std::max(p1, fac) * absh;
        absh = std::max(hmin, std::min(rstate.hmax, absh));
        errold = err;
        hold = h;
        h = tdir * absh;

        if (last) {
            std::cout << "setting idid in last" << std::endl;
            rstate.idid = 1;
            return rstate.idid;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rkc (burn_t& state, rkc_t& rstate)
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

    std::cout << "calling rkclow" << std::endl;

    rkclow(state, rstate);

    std::cout << rstate.idid << std::endl;
    return rstate.idid;
}

#endif
