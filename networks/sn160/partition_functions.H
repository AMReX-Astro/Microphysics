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
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
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
        AMREX_GPU_HOST_DEVICE AMREX_INLINE
        pf_cache_t() {
            // The entries will be default-initialized to -900, log(partition_function)
            // should never be -900.
            for (int i = 1; i <= NumSpecTotal; ++i) {
                    data(i,1) = -900.0_rt;
                    data(i,2) = -900.0_rt;
            }
        }
    };

}

// get the partition function for nucleus inuc.  Here pf_cache
// should already contain the index into the temp / pf arrays

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_partition_function(const int inuc,
                            [[maybe_unused]] const tf_t& tfactors,
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


// a version of get_parition_function where we have not cached the
// index into the arrays. And returns pf and dpf_dT instead of
// logpf and dlogpf_dT9

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_partition_function(const int inuc,
                            [[maybe_unused]] const tf_t& tfactors,
                            amrex::Real& pf, amrex::Real& dpf_dT) {

    // create a dummy cache.  This will initialize the indices
    // to -1 and we will force a search for the correct index
    // when we do the interpolation.

    part_fun::pf_cache_t pf_cache;

    amrex::Real logpf, dlogpf_dT9;
    get_partition_function(inuc, tfactors, pf_cache, logpf, dlogpf_dT9);

    // Convert to logpf to pf and dlogpf_dT9 to dpf_dT
    pf = std::exp(logpf);
    dpf_dT = pf * dlogpf_dT9 * 1.e-9_rt;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_partition_function_cached(const int inuc, const tf_t& tfactors,
                                   part_fun::pf_cache_t& pf_cache,
                                   amrex::Real& logpf, amrex::Real& dlogpf_dT9) {
    if (pf_cache.data(inuc, 1) != -900.0_rt) {
        // present in cache
        amrex::ignore_unused(tfactors);
        logpf = pf_cache.data(inuc, 1);
        dlogpf_dT9 = pf_cache.data(inuc, 2);
    } else {
        get_partition_function(inuc, tfactors, pf_cache, logpf, dlogpf_dT9);
        pf_cache.data(inuc, 1) = logpf;
        pf_cache.data(inuc, 2) = dlogpf_dT9;
    }
}


template <int spec>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr amrex::Real get_spin_state() {

    static_assert(spec >= 1 && spec <= NumSpec);

    // Set the spin states (2J + 1) of the element

    if constexpr (
        spec == He4 || spec == C12 || spec == C14 ||
        spec == O14 || spec == O16 || spec == O18 ||
        spec == Ne18 || spec == Ne20 || spec == Ne22 ||
        spec == Mg24 || spec == Mg26 || spec == Si28 ||
        spec == Si30 || spec == Si32 || spec == S32 ||
        spec == S34 || spec == S36 || spec == Cl34 ||
        spec == Ar36 || spec == Ar38 || spec == Ar40 ||
        spec == Ca40 || spec == Ca42 || spec == Ca44 ||
        spec == Ca46 || spec == Ca48 || spec == Ti44 ||
        spec == Ti46 || spec == Ti48 || spec == Ti50 ||
        spec == V46 || spec == Cr48 || spec == Cr50 ||
        spec == Cr52 || spec == Cr54 || spec == Mn50 ||
        spec == Fe52 || spec == Fe54 || spec == Fe56 ||
        spec == Fe58 || spec == Co54 || spec == Ni56 ||
        spec == Ni58 || spec == Ni60 || spec == Ni62 ||
        spec == Ni64 || spec == Zn60 || spec == Zn62 ||
        spec == Zn64 || spec == Zn66 || spec == Ga62 ||
        spec == Ga64 || spec == Ge64
                  )
    {
        return 1.0_rt;
    }
    else if constexpr (
        spec == N || spec == H1 || spec == He3 ||
        spec == C13 || spec == N13 || spec == N15 ||
        spec == O15 || spec == F19 || spec == Ne19 ||
        spec == Si29 || spec == P29 || spec == P31 ||
        spec == P33 || spec == Fe57 || spec == Ni63
                       )
    {
        return 2.0_rt;
    }
    else if constexpr (
        spec == H2 || spec == Li6 || spec == N14 ||
        spec == F18 || spec == P30 || spec == P32 ||
        spec == Cu58 || spec == Cu62 || spec == Cu64
                       )
    {
        return 3.0_rt;
    }
    else if constexpr (
        spec == Li7 || spec == Be7 || spec == Be9 ||
        spec == B11 || spec == Ne21 || spec == Na21 ||
        spec == Na23 || spec == Mg23 || spec == Si31 ||
        spec == S33 || spec == S35 || spec == Cl33 ||
        spec == Cl35 || spec == Cl37 || spec == Ar37 ||
        spec == K37 || spec == K39 || spec == K41 ||
        spec == Ti51 || spec == V47 || spec == Cr53 ||
        spec == Fe55 || spec == Ni57 || spec == Ni59 ||
        spec == Ni61 || spec == Cu57 || spec == Cu59 ||
        spec == Cu61 || spec == Cu63 || spec == Cu65 ||
        spec == Zn59 || spec == Zn61 || spec == Zn63 ||
        spec == Ga63 || spec == Ge63
                       )
    {
        return 4.0_rt;
    }
    else if constexpr (
        spec == B8 || spec == Cl36 || spec == Sc44 ||
        spec == Co58 || spec == Cu60
                       )
    {
        return 5.0_rt;
    }
    else if constexpr (
        spec == O17 || spec == F17 || spec == Mg25 ||
        spec == Al25 || spec == Al27 || spec == Ti47 ||
        spec == Cr49 || spec == Mn51 || spec == Mn55 ||
        spec == Zn65
                       )
    {
        return 6.0_rt;
    }
    else if constexpr (
        spec == B10 || spec == Na22 || spec == K38 ||
        spec == V52 || spec == Mn54
                       )
    {
        return 7.0_rt;
    }
    else if constexpr (
        spec == Ar39 || spec == Ca41 || spec == Ca43 ||
        spec == Ca45 || spec == Ca47 || spec == Sc43 ||
        spec == Sc45 || spec == Sc47 || spec == Sc49 ||
        spec == Ti45 || spec == Ti49 || spec == V49 ||
        spec == V51 || spec == Cr51 || spec == Mn53 ||
        spec == Fe53 || spec == Co53 || spec == Co55 ||
        spec == Co57 || spec == Co59
                       )
    {
        return 8.0_rt;
    }
    else if constexpr (
        spec == K40 || spec == Sc46 || spec == V48 ||
        spec == Co56
                       )
    {
        return 9.0_rt;
    }
    else if constexpr (
        spec == Al26
                       )
    {
        return 11.0_rt;
    }
    else if constexpr (
        spec == Sc48 || spec == V50 || spec == Mn52
                       )
    {
        return 13.0_rt;
    }

    // Return -1 if we don't recognize the species.
    return -1.0_rt;
}


// Legacy (non-templated) interface

AMREX_GPU_HOST_DEVICE AMREX_INLINE
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
