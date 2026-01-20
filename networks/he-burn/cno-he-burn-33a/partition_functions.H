#ifndef PARTITION_FUNCTIONS_H
#define PARTITION_FUNCTIONS_H

#include <numbers>

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Math.H>

#include <tfactors.H>
#include <fundamental_constants.H>
#include <network_properties.H>
#include <interp_tools.H>

using namespace amrex::literals;
using namespace Species;

namespace part_fun {

    constexpr int npts_1 = 72;

    // this is T9

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> temp_array_1;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> O14_pf_array;
    constexpr amrex::Real O14_pf_threshold_T9 = 10.0;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> O15_pf_array;
    constexpr amrex::Real O15_pf_threshold_T9 = 10.0;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> O16_pf_array;
    constexpr amrex::Real O16_pf_threshold_T9 = 10.0;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> O17_pf_array;
    constexpr amrex::Real O17_pf_threshold_T9 = 2.0;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> O18_pf_array;
    constexpr amrex::Real O18_pf_threshold_T9 = 3.0;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> F17_pf_array;
    constexpr amrex::Real F17_pf_threshold_T9 = 1.0;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> F18_pf_array;
    constexpr amrex::Real F18_pf_threshold_T9 = 1.5;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> F19_pf_array;
    constexpr amrex::Real F19_pf_threshold_T9 = 0.2;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Ne18_pf_array;
    constexpr amrex::Real Ne18_pf_threshold_T9 = 1.0;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Ne19_pf_array;
    constexpr amrex::Real Ne19_pf_threshold_T9 = 0.15;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Ne20_pf_array;
    constexpr amrex::Real Ne20_pf_threshold_T9 = 1.0;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Ne21_pf_array;
    constexpr amrex::Real Ne21_pf_threshold_T9 = 0.2;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Na22_pf_array;
    constexpr amrex::Real Na22_pf_threshold_T9 = 0.4;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Na23_pf_array;
    constexpr amrex::Real Na23_pf_threshold_T9 = 0.3;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Mg22_pf_array;
    constexpr amrex::Real Mg22_pf_threshold_T9 = 0.8;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Mg24_pf_array;
    constexpr amrex::Real Mg24_pf_threshold_T9 = 0.9;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Al27_pf_array;
    constexpr amrex::Real Al27_pf_threshold_T9 = 0.7;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Si28_pf_array;
    constexpr amrex::Real Si28_pf_threshold_T9 = 1.0;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> P31_pf_array;
    constexpr amrex::Real P31_pf_threshold_T9 = 0.9;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> S32_pf_array;
    constexpr amrex::Real S32_pf_threshold_T9 = 1.5;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Cl35_pf_array;
    constexpr amrex::Real Cl35_pf_threshold_T9 = 1.0;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Ar36_pf_array;
    constexpr amrex::Real Ar36_pf_threshold_T9 = 1.0;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> K39_pf_array;
    constexpr amrex::Real K39_pf_threshold_T9 = 2.0;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Ca40_pf_array;
    constexpr amrex::Real Ca40_pf_threshold_T9 = 2.5;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Sc43_pf_array;
    constexpr amrex::Real Sc43_pf_threshold_T9 = 0.01;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Ti44_pf_array;
    constexpr amrex::Real Ti44_pf_threshold_T9 = 0.7;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> V47_pf_array;
    constexpr amrex::Real V47_pf_threshold_T9 = -1.0;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Cr48_pf_array;
    constexpr amrex::Real Cr48_pf_threshold_T9 = 0.5;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Mn51_pf_array;
    constexpr amrex::Real Mn51_pf_threshold_T9 = 0.15;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Fe52_pf_array;
    constexpr amrex::Real Fe52_pf_threshold_T9 = 0.6;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Co55_pf_array;
    constexpr amrex::Real Co55_pf_threshold_T9 = 1.5;

    // this is log(partition function)

    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> Ni56_pf_array;
    constexpr amrex::Real Ni56_pf_threshold_T9 = 1.5;


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
        int index_temp_array_1{-1};

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

    case O14:
        if (tfactors.T9 > part_fun::O14_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::O14_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case O15:
        if (tfactors.T9 > part_fun::O15_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::O15_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case O16:
        if (tfactors.T9 > part_fun::O16_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::O16_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case O17:
        if (tfactors.T9 > part_fun::O17_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::O17_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case O18:
        if (tfactors.T9 > part_fun::O18_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::O18_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case F17:
        if (tfactors.T9 > part_fun::F17_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::F17_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case F18:
        if (tfactors.T9 > part_fun::F18_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::F18_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case F19:
        if (tfactors.T9 > part_fun::F19_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::F19_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Ne18:
        if (tfactors.T9 > part_fun::Ne18_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Ne18_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Ne19:
        if (tfactors.T9 > part_fun::Ne19_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Ne19_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Ne20:
        if (tfactors.T9 > part_fun::Ne20_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Ne20_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Ne21:
        if (tfactors.T9 > part_fun::Ne21_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Ne21_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Na22:
        if (tfactors.T9 > part_fun::Na22_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Na22_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Na23:
        if (tfactors.T9 > part_fun::Na23_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Na23_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Mg22:
        if (tfactors.T9 > part_fun::Mg22_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Mg22_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Mg24:
        if (tfactors.T9 > part_fun::Mg24_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Mg24_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Al27:
        if (tfactors.T9 > part_fun::Al27_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Al27_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Si28:
        if (tfactors.T9 > part_fun::Si28_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Si28_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case P31:
        if (tfactors.T9 > part_fun::P31_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::P31_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case S32:
        if (tfactors.T9 > part_fun::S32_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::S32_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Cl35:
        if (tfactors.T9 > part_fun::Cl35_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Cl35_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Ar36:
        if (tfactors.T9 > part_fun::Ar36_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Ar36_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case K39:
        if (tfactors.T9 > part_fun::K39_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::K39_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Ca40:
        if (tfactors.T9 > part_fun::Ca40_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Ca40_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Sc43:
        if (tfactors.T9 > part_fun::Sc43_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Sc43_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Ti44:
        if (tfactors.T9 > part_fun::Ti44_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Ti44_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case V47:
        if (tfactors.T9 > part_fun::V47_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::V47_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Cr48:
        if (tfactors.T9 > part_fun::Cr48_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Cr48_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Mn51:
        if (tfactors.T9 > part_fun::Mn51_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Mn51_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Fe52:
        if (tfactors.T9 > part_fun::Fe52_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Fe52_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Co55:
        if (tfactors.T9 > part_fun::Co55_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Co55_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case Ni56:
        if (tfactors.T9 > part_fun::Ni56_pf_threshold_T9) {
            part_fun::interpolate_pf(tfactors.T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::Ni56_pf_array, logpf, dlogpf_dT9);
        }
        break;


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
        spec == He4 || spec == C12 || spec == O14 ||
        spec == O16 || spec == O18 || spec == Ne18 ||
        spec == Ne20 || spec == Mg22 || spec == Mg24 ||
        spec == Si28 || spec == S32 || spec == Ar36 ||
        spec == Ca40 || spec == Ti44 || spec == Cr48 ||
        spec == Fe52 || spec == Ni56
                  )
    {
        return 1.0_rt;
    }
    else if constexpr (
        spec == H1 || spec == C13 || spec == N13 ||
        spec == N15 || spec == O15 || spec == F19 ||
        spec == Ne19 || spec == P31
                       )
    {
        return 2.0_rt;
    }
    else if constexpr (
        spec == N14 || spec == F18
                       )
    {
        return 3.0_rt;
    }
    else if constexpr (
        spec == Ne21 || spec == Na23 || spec == Cl35 ||
        spec == K39 || spec == V47
                       )
    {
        return 4.0_rt;
    }
    else if constexpr (
        spec == O17 || spec == F17 || spec == Al27 ||
        spec == Mn51
                       )
    {
        return 6.0_rt;
    }
    else if constexpr (
        spec == Na22
                       )
    {
        return 7.0_rt;
    }
    else if constexpr (
        spec == Sc43 || spec == Co55
                       )
    {
        return 8.0_rt;
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
