#ifndef RKC_UTIL_H
#define RKC_UTIL_H

#include <rkc_type.H>

template <typename BurnT, typename IntT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real rkc_init_dt (BurnT& state, IntT& rstate, const Real max_dt, const Real sprad)
{

    // estimate the initial timestep
    //
    // this requires that yn hold the initial state and fn the corresponding RHS
    // it also requires temporary space yjm1 and yjm2
    // finally, we require the initial spectral radius, sprad, and the maximum
    // timestep (tend - tstart), max_dt.

    Real hmin = 10.0_rt * UROUND * std::max(std::abs(rstate.t), max_dt);

    Real absh = max_dt;

    if (sprad * absh > 1.0_rt) {
        absh = 1.0_rt / sprad;
    }
    absh = std::max(absh, hmin);

    for (int i = 1; i <= INT_NEQS; ++i) {
        rstate.yjm1(i) = rstate.yn(i) + absh * rstate.fn(i);
    }

    // we want to call with yjm1 as the input, so we need to stash the
    // current y
    RArray1D ysav;
    for (int i = 1; i <= INT_NEQS; ++i) {
        ysav(i) = rstate.y(i);
        rstate.y(i) = rstate.yjm1(i);
    }
    rhs(rstate.t+absh, state, rstate, rstate.yjm2);
    for (int i = 1; i <= INT_NEQS; ++i) {
        rstate.y(i) = ysav(i);
    }

    rstate.nfe++;
    Real est{};

    // compute the weights using the tolerances
    for (int i = 1; i <= INT_NEQS; ++i) {
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
        est += std::pow((rstate.yjm2(i) - rstate.fn(i)) / wt, 2);
    }
    est = absh * std::sqrt(est / INT_NEQS);

    if (0.1_rt * absh < max_dt * std::sqrt(est)) {
        absh = std::max(0.1_rt * absh / std::sqrt(est), hmin);
    } else {
        absh = max_dt;
    }

    return absh;

}

template <typename BurnT, typename IntT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real rkcrho (BurnT& state, IntT& rstate, const Real max_dt)
{

    // RKCRHO attempts to compute a close upper bound, SPRAD, on
    // the spectral radius of the Jacobian matrix using a nonlinear
    // power method.  A convergence failure is reported by IDID = 6.

    // for the start (nsteps == 0), this works off of the initial
    // function evaluation, rstate.fn(i), which has units of dy/dt
    //
    // for other steps, with works off of rstate.sprad, which is the
    // eigenvector from the previous solve.  It also has units of
    // y
    //
    // this means that the caller needs to initialize
    //
    // rstate.yn to the old-timesolution
    // rstate.fn to the RHS called with yn for nsteps == 0
    // rstate.sprad to the previous eigenvector for nsteps > 0
    //
    // max_dt is the maximum timestep, typically tend - tbegin
    //
    // we use rstate.yjm1 and rstate.yjm2 as internal storage

    constexpr int itmax{50};

    // sprad smaller than small = 1/hmax are not
    // interesting because they do not constrain the step size.

    Real small = 1.0_rt / max_dt;

    // The initial slope is used as guess when nsteps = 0 and
    // thereafter the last computed eigenvector.  Some care
    // is needed to deal with special cases. Approximations to
    // the eigenvector are normalized so that their Euclidean
    // norm has the constant value dynrm.

    if (rstate.nsteps == 0) {
        for (int i = 1; i <= RKC_NEQS; ++i) {
            rstate.yjm1(i) = rstate.fn(i);
        }
    } else {
        for (int i = 1; i <= RKC_NEQS; ++i) {
            rstate.yjm1(i) = rstate.sprad(i);
        }
    }

    Real ynrm{};
    Real vnrm{};
    for (int i = 1; i <= RKC_NEQS; ++i) {
        ynrm += std::pow(rstate.yn(i), 2);
        vnrm += std::pow(rstate.yjm1(i), 2);
    }
    ynrm = std::sqrt(ynrm);
    vnrm = std::sqrt(vnrm);

    Real dynrm{};
    if (ynrm != 0.0_rt && vnrm != 0.0_rt) {
        dynrm = ynrm * std::sqrt(UROUND);
        for (int i = 1; i <= RKC_NEQS; ++i) {
            rstate.yjm1(i) = rstate.yn(i) + rstate.yjm1(i) * (dynrm / vnrm);
        }
    } else if (ynrm != 0.0_rt) {
        dynrm = ynrm * std::sqrt(UROUND);
        for (int i = 1; i <= RKC_NEQS; ++i) {
            rstate.yjm1(i) = rstate.yn(i) + rstate.yn(i) * std::sqrt(UROUND);
        }
    } else if (vnrm != 0.0_rt) {
        dynrm = UROUND;
        for (int i = 1; i <= RKC_NEQS; ++i) {
            rstate.yjm1(i) = rstate.yjm1(i) * (dynrm / vnrm);
        }
    } else {
        dynrm = UROUND;
        for (int i = 1; i <= RKC_NEQS; ++i) {
            rstate.yjm1(i) = dynrm;
        }
    }

    // Now iterate with a nonlinear power method.

    Real sigma{};
    Real sprad{};

    for (int iter = 0; iter < itmax; ++iter) {

        // we want to call the RHS with yjm1() as the input, so we need to copy
        // this temporarily into vstate.y()
        RArray1D ysav;
        for (int i = 1; i <= RKC_NEQS; ++i) {
            ysav(i) = rstate.y(i);
            rstate.y(i) = rstate.yjm1(i);
        }
        rhs(rstate.t, state, rstate, rstate.yjm2);
        for (int i = 1; i <= RKC_NEQS; ++i) {
            rstate.y(i) = ysav(i);
        }
        rstate.nfesig++;
        Real dfnrm{};
        for (int i = 1; i <= RKC_NEQS; ++i) {
            dfnrm += std::pow(rstate.yjm2(i) - rstate.fn(i), 2);
        }
        dfnrm = std::sqrt(dfnrm);
        Real sigmal = sigma;
        sigma = dfnrm / dynrm;

        // sprad is a little bigger than the estimate sigma of the
        // spectral radius, so is more likely to be an upper bound.

        sprad = 1.2_rt * sigma;

        if (iter >= 1 &&
            std::abs(sigma - sigmal) <= std::max(sigma, small) * 0.01_rt) {
            for (int i = 1; i <= RKC_NEQS; ++i) {
                rstate.sprad(i) = rstate.yjm1(i) - rstate.yn(i);
            }
            return sprad;
        }

        // The next yjm1(*) is the change in f
        // scaled so that norm(v - yn) = dynrm.

        if (dfnrm != 0.0_rt) {
            for (int i = 1; i <= RKC_NEQS; ++i) {
                rstate.yjm1(i) = rstate.yn(i) +
                    (rstate.yjm2(i) - rstate.fn(i)) * (dynrm / dfnrm);
            }
        } else {

            // The new v(*) degenerated to yn(*)--"randomly" perturb
            // current approximation to the eigenvector by changing
            // the sign of one component.

            int index = 1 + (iter + 1) % RKC_NEQS;
            rstate.yjm1(index) = rstate.yn(index) -
                (rstate.yjm1(index) - rstate.yn(index));
        }
    }

    // Set flag to report a convergence failure.

    rstate.idid = 6;
    return sprad;
}


#endif
