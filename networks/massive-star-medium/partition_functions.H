#ifndef PARTITION_FUNCTIONS_H
#define PARTITION_FUNCTIONS_H

#include <numbers>

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Math.H>
#include <AMReX_ConstexprFor.H>

#include <tfactors.H>
#include <fundamental_constants.H>
#include <network_properties.H>
#include <interp_tools.H>

using namespace amrex::literals;
using namespace Species;

namespace part_fun {


    // interpolation routine

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void interpolate_pf(const amrex::Real t9, int idx,
                        const T& temp_array, const T& pf_array,
                        amrex::Real& logpf, amrex::Real& dlogpf_dT9) {

        // if the index is -1, we are either outside of the bounds or
        // we never computed it.  try recomputing here
        if (idx == -1) {
            idx = interp_net::find_index(t9, temp_array);
        }

        // if the index is still -1, then we are out-of-bounds, so
        // just skip the interpolation

        if (idx != -1) {

            // we have idx such that:
            // temp_array[idx] <= t9 < temp_array[idx+1]
            AMREX_ASSERT(t9 >= temp_array(idx) && t9 < temp_array(idx+1));

            // find the derivative: dlog(pf)/dT9 -- this is
            // (log(pf_{i+1}) - log(pf_i)) / (T_{i+1} - T_i)

            dlogpf_dT9 = (pf_array(idx+1) - pf_array(idx)) /
                           (temp_array(idx+1) - temp_array(idx));

            // find log(PF)

            logpf = pf_array(idx) + dlogpf_dT9 * (t9 - temp_array(idx));

        } else {

            // T < the smallest T or >= the largest T in the partition function table
            logpf = 0.0;
            dlogpf_dT9 = 0.0;

        }

    }

    struct pf_cache_t {
        // Store the coefficient and derivative adjacent in memory, as they're
        // always accessed at the same time.
        amrex::Array2D<amrex::Real, 1, NumSpecTotal, 1, 2, amrex::Order::C> data;
        AMREX_GPU_DEVICE AMREX_INLINE
        pf_cache_t() {
            // The entries will be default-initialized to -900,
            // log(partition_function) should never be -900.
            for (int i = 1; i <= NumSpecTotal; ++i) {
                    data(i, 1) = -900.0_rt;
                    data(i, 2) = -900.0_rt;
            }
        }
    };

}

// get the partition function for nucleus inuc.  Here pf_cache
// should already contain the index into the temp / pf arrays

AMREX_GPU_DEVICE AMREX_INLINE
void get_partition_function(const int inuc,
                            [[maybe_unused]] const amrex::Real T9,
                            const part_fun::pf_cache_t& pf_cache,
                            amrex::Real& logpf, amrex::Real& dlogpf_dT9) {

    logpf = 0.0_rt;
    dlogpf_dT9 = 0.0_rt;

    // inuc is the 1-based index for the species

    switch (inuc) {


    default:

        logpf = 0.0_rt;
        dlogpf_dT9 = 0.0_rt;

    }

}


// a version of get_partition_function where we have not cached the
// index into the arrays. Uses tfactors instead of T9 as argument,
// and returns pf and dpf_dT instead of logpf and dlogpf_dT9

AMREX_GPU_DEVICE AMREX_INLINE
void get_partition_function(const int inuc,
                            [[maybe_unused]] const tf_t& tfactors,
                            amrex::Real& pf, amrex::Real& dpf_dT) {

    // create a dummy cache.  This will initialize the indices
    // to -1 and we will force a search for the correct index
    // when we do the interpolation.

    part_fun::pf_cache_t pf_cache;

    amrex::Real logpf, dlogpf_dT9;
    get_partition_function(inuc, tfactors.T9, pf_cache, logpf, dlogpf_dT9);

    // Convert logpf to pf and dlogpf_dT9 to dpf_dT
    pf = std::exp(logpf);
    dpf_dT = pf * dlogpf_dT9 * 1.e-9_rt;
}


AMREX_GPU_DEVICE AMREX_INLINE
void get_partition_function_cached(const int inuc, const amrex::Real T9,
                                   part_fun::pf_cache_t& pf_cache,
                                   amrex::Real& logpf, amrex::Real& dlogpf_dT9) {
    if (pf_cache.data(inuc, 1) != -900.0_rt) {
        // present in cache
        amrex::ignore_unused(T9);
        logpf = pf_cache.data(inuc, 1);
        dlogpf_dT9 = pf_cache.data(inuc, 2);
    } else {
        get_partition_function(inuc, T9, pf_cache, logpf, dlogpf_dT9);
        pf_cache.data(inuc, 1) = logpf;
        pf_cache.data(inuc, 2) = dlogpf_dT9;
    }
}


template <int spec>
AMREX_GPU_DEVICE AMREX_INLINE
constexpr amrex::Real get_spin_state() {

    static_assert(spec >= 1 && spec <= NumSpec);

    // Set the spin states (2J + 1) of the element

    if constexpr (
        spec == He4 || spec == C12 || spec == O16 ||
        spec == Ne20 || spec == Mg24 || spec == Si28 ||
        spec == Si30 || spec == S32 || spec == S34 ||
        spec == Ar36 || spec == Ar38 || spec == Ca40 ||
        spec == Ca42 || spec == Ti44 || spec == Ti46 ||
        spec == Ti48 || spec == Cr48 || spec == Cr50 ||
        spec == Cr52 || spec == Cr54 || spec == Fe52 ||
        spec == Fe54 || spec == Fe56 || spec == Fe58 ||
        spec == Ni56 || spec == Ni58 || spec == Ni60 ||
        spec == Zn60 || spec == Sc42
                  )
    {
        return 1.0_rt;
    }
    else if constexpr (
        spec == N || spec == H1 || spec == Si29 ||
        spec == P29 || spec == P31 || spec == P33 ||
        spec == Fe57
                       )
    {
        return 2.0_rt;
    }
    else if constexpr (
        spec == Mg23 || spec == Si31 || spec == S33 ||
        spec == S35 || spec == Cl33 || spec == Cl35 ||
        spec == Cl37 || spec == Ar37 || spec == K37 ||
        spec == K39 || spec == V47 || spec == Cu59 ||
        spec == Ni57 || spec == Ni59 || spec == Fe55 ||
        spec == Cr53 || spec == Na23
                       )
    {
        return 4.0_rt;
    }
    else if constexpr (
        spec == Co58 || spec == Sc44
                       )
    {
        return 5.0_rt;
    }
    else if constexpr (
        spec == Al27 || spec == Si27 || spec == Mn51 ||
        spec == Mn55 || spec == Cr49 || spec == Ti47
                       )
    {
        return 6.0_rt;
    }
    else if constexpr (
        spec == Mn54
                       )
    {
        return 7.0_rt;
    }
    else if constexpr (
        spec == Sc41 || spec == Sc43 || spec == Sc45 ||
        spec == V49 || spec == V51 || spec == Mn53 ||
        spec == Co55 || spec == Co57 || spec == Co59 ||
        spec == Fe53 || spec == Cr51 || spec == Ti45 ||
        spec == Ca41
                       )
    {
        return 8.0_rt;
    }
    else if constexpr (
        spec == Co56 || spec == V48
                       )
    {
        return 9.0_rt;
    }
    else if constexpr (
        spec == Mn52 || spec == V50
                       )
    {
        return 13.0_rt;
    }

    // Return -1 if we don't recognize the species.
    return -1.0_rt;
}


// Legacy (non-templated) interface

AMREX_GPU_DEVICE AMREX_INLINE
amrex::Real get_spin_state (int spec) {

    amrex::Real s = -1.0_rt;

    // Set the spin states (2J + 1) of the element

    amrex::constexpr_for<1, NumSpec+1>([&] (auto n) {
        if (n == spec) {
            s = get_spin_state<n>();
        }
    });

    return s;
}


#endif
