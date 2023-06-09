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
    // it also requires temporary space vtemp1 and vtemp2
    // finally, we require the initial spectral radius, sprad, and the maximum
    // timestep (tend - tstart), max_dt.

    Real hmin = 10.0_rt * UROUND * std::max(std::abs(rstate.t), max_dt);

    Real absh = max_dt;

    if (sprad * absh > 1.0_rt) {
        absh = 1.0_rt / sprad;
    }
    absh = std::max(absh, hmin);

    for (int i = 1; i <= INT_NEQS; ++i) {
        rstate.vtemp1(i) = rstate.yn(i) + absh * rstate.fn(i);
    }

    // we want to call with vtemp1 as the input, so we need to stash the
    // current y
    RArray1D ysav;
    for (int i = 1; i <= INT_NEQS; ++i) {
        ysav(i) = rstate.y(i);
        rstate.y(i) = rstate.vtemp1(i);
    }
    rhs(rstate.t+absh, state, rstate, rstate.vtemp2);
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
        est += std::pow((rstate.vtemp2(i) - rstate.fn(i)) / wt, 2);
    }
    est = absh * std::sqrt(est / INT_NEQS);

    if (0.1_rt * absh < max_dt * std::sqrt(est)) {
        absh = std::max(0.1_rt * absh / std::sqrt(est), hmin);
    } else {
        absh = max_dt;
    }

    return absh;

}

#endif
