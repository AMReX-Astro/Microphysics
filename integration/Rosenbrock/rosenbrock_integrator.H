#ifndef ROSENBROCK_INTEGRATOR_H
#define ROSENBROCK_INTEGRATOR_H

#include <AMReX_Algorithm.H>
#include <AMReX_GpuDevice.H>

#include <cstdint>
#include <cmath>
#include <limits>

#include <network.H>
#include <burn_type.H>
#include <linpack.H>

#ifdef STRANG
#include <integrator_rhs_strang.H>
#endif
#ifdef SDC
#include <integrator_rhs_sdc.H>
#endif

#include <integrator_data.H>
#include <rosenbrock_type.H>
#include <rosenbrock_tableau.H>

namespace rosenbrock {

namespace warp {

using mask_t = unsigned long long;

#ifdef AMREX_USE_GPU
static_assert(amrex::Gpu::Device::warp_size == 32 ||
              amrex::Gpu::Device::warp_size == 64,
              "Rosenbrock warp-coupled integration expects a GPU warp/wavefront size of 32 or 64");
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
constexpr int size ()
{
#ifdef AMREX_USE_GPU
    return amrex::Gpu::Device::warp_size;
#else
    return 1;
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
mask_t active_mask ()
{
    mask_t mask = 1ULL;
#if defined(AMREX_USE_CUDA)
    AMREX_IF_ON_DEVICE((mask = static_cast<mask_t>(__activemask());))
#elif defined(AMREX_USE_HIP)
    AMREX_IF_ON_DEVICE((mask = static_cast<mask_t>(__ballot(1));))
#endif
    return mask;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool is_active (const mask_t mask, const int lane)
{
    return (mask & (mask_t{1} << lane)) != 0;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int active_count (const mask_t mask)
{
    int count = 0;
    for (int lane = 0; lane < size(); ++lane) {
        if (is_active(mask, lane)) {
            count += 1;
        }
    }
    return count;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int first_lane (const mask_t mask)
{
    for (int lane = 0; lane < size(); ++lane) {
        if (is_active(mask, lane)) {
            return lane;
        }
    }
    return 0;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T broadcast (const T value, const int src_lane, const mask_t mask)
{
#if defined(AMREX_USE_CUDA)
    AMREX_IF_ON_DEVICE((return __shfl_sync(static_cast<unsigned>(mask), value, src_lane);))
#elif defined(AMREX_USE_HIP)
    AMREX_IF_ON_DEVICE((return __shfl(value, src_lane);))
#endif
    return value;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T sum (const T value, const mask_t mask)
{
    T result = T(0);
    for (int lane = 0; lane < size(); ++lane) {
        if (is_active(mask, lane)) {
            result += broadcast(value, lane, mask);
        }
    }
    return result;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T max (const T value, const mask_t mask)
{
    T result = broadcast(value, first_lane(mask), mask);
    for (int lane = 0; lane < size(); ++lane) {
        if (is_active(mask, lane)) {
            result = amrex::max(result, broadcast(value, lane, mask));
        }
    }
    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool all (const bool value, const mask_t mask)
{
    int result = 1;
    const int ivalue = value ? 1 : 0;
    for (int lane = 0; lane < size(); ++lane) {
        if (is_active(mask, lane)) {
            result = result && broadcast(ivalue, lane, mask);
        }
    }
    return result != 0;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool any (const bool value, const mask_t mask)
{
    int result = 0;
    const int ivalue = value ? 1 : 0;
    for (int lane = 0; lane < size(); ++lane) {
        if (is_active(mask, lane)) {
            result = result || broadcast(ivalue, lane, mask);
        }
    }
    return result != 0;
}

} // namespace warp

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real& stage (rosenbrock_t<int_neqs>& rstate, const int i, const int n)
{
    if (i == 1) {
        return rstate.ak1(n);
    }
    if (i == 2) {
        return rstate.ak2(n);
    }
    if (i == 3) {
        return rstate.ak3(n);
    }
    if (i == 4) {
        return rstate.ak4(n);
    }
    if (i == 5) {
        return rstate.ak5(n);
    }
    if (i == 6) {
        return rstate.ak6(n);
    }
    if (i == 7) {
        return rstate.ak7(n);
    }
    if (i == 8) {
        return rstate.ak8(n);
    }
    return rstate.work(n);
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real stage (const rosenbrock_t<int_neqs>& rstate, const int i, const int n)
{
    if (i == 1) {
        return rstate.ak1(n);
    }
    if (i == 2) {
        return rstate.ak2(n);
    }
    if (i == 3) {
        return rstate.ak3(n);
    }
    if (i == 4) {
        return rstate.ak4(n);
    }
    if (i == 5) {
        return rstate.ak5(n);
    }
    if (i == 6) {
        return rstate.ak6(n);
    }
    if (i == 7) {
        return rstate.ak7(n);
    }
    if (i == 8) {
        return rstate.ak8(n);
    }
    return rstate.work(n);
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
RArray1D& stage_vector (rosenbrock_t<int_neqs>& rstate, const int i)
{
    if (i == 1) {
        return rstate.ak1;
    }
    if (i == 2) {
        return rstate.ak2;
    }
    if (i == 3) {
        return rstate.ak3;
    }
    if (i == 4) {
        return rstate.ak4;
    }
    if (i == 5) {
        return rstate.ak5;
    }
    if (i == 6) {
        return rstate.ak6;
    }
    if (i == 7) {
        return rstate.ak7;
    }
    if (i == 8) {
        return rstate.ak8;
    }
    return rstate.work;
}

template <typename RosenbrockT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real rtol_for (const RosenbrockT& rstate, const int n)
{
    return n == net_ienuc ? rstate.rtol_enuc : rstate.rtol_spec;
}

template <typename RosenbrockT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real atol_for (const RosenbrockT& rstate, const int n)
{
    return n == net_ienuc ? rstate.atol_enuc : rstate.atol_spec;
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real error_norm (const rosenbrock_t<int_neqs>& rstate)
{
    amrex::Real err = 0.0_rt;
    for (int n = 1; n <= int_neqs; ++n) {
        const amrex::Real sk = atol_for(rstate, n) +
            rtol_for(rstate, n) * amrex::max(std::abs(rstate.y(n)), std::abs(rstate.ynew(n)));
        const amrex::Real term = rstate.work(n) / sk;
        err += term * term;
    }
    return std::sqrt(err / static_cast<amrex::Real>(int_neqs));
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real error_norm_warp (const rosenbrock_t<int_neqs>& rstate,
                             const warp::mask_t mask, const int active_lanes)
{
    amrex::Real err = 0.0_rt;
    for (int n = 1; n <= int_neqs; ++n) {
        const amrex::Real sk = atol_for(rstate, n) +
            rtol_for(rstate, n) * amrex::max(std::abs(rstate.y(n)), std::abs(rstate.ynew(n)));
        const amrex::Real term = rstate.work(n) / sk;
        err += term * term;
    }

    const amrex::Real warp_err = warp::sum(err, mask);
    return std::sqrt(warp_err / static_cast<amrex::Real>(active_lanes * int_neqs));
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool species_change_valid (const rosenbrock_t<int_neqs>& rstate)
{
    bool valid = true;

#ifdef STRANG
    constexpr amrex::Real increase_change_factor = 4.0_rt;
    constexpr amrex::Real decrease_change_factor = 0.25_rt;

    for (int i = 1; i <= NumSpec; ++i) {
        if (std::abs(rstate.y(i)) > integrator_rp::X_reject_buffer * rstate.atol_spec &&
            std::abs(rstate.ynew(i)) > integrator_rp::X_reject_buffer * rstate.atol_spec &&
            (std::abs(rstate.ynew(i)) > increase_change_factor * std::abs(rstate.y(i)) ||
             std::abs(rstate.ynew(i)) < decrease_change_factor * std::abs(rstate.y(i)))) {
            valid = false;
        }
    }
#endif

    return valid;
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real yass_change_norm (const rosenbrock_t<int_neqs>& rstate)
{
    amrex::Real max_change = 0.0_rt;
    amrex::Real max_species = 0.0_rt;

    for (int i = 1; i <= NumSpec && i <= int_neqs; ++i) {
        max_species = amrex::max(max_species, std::abs(rstate.y(i)));
    }

    const amrex::Real active_floor = integrator_rp::yass_floor * max_species;

    for (int i = 1; i <= NumSpec && i <= int_neqs; ++i) {
        const amrex::Real yi = std::abs(rstate.y(i));
        if (yi > active_floor) {
            max_change = amrex::max(max_change, std::abs(rstate.ynew(i) - rstate.y(i)) / yi);
        }
    }

    return max_change / integrator_rp::yass_epsilon;
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real yass_change_norm_warp (const rosenbrock_t<int_neqs>& rstate,
                                   const warp::mask_t mask)
{
    return warp::max(yass_change_norm(rstate), mask);
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool valid_integrator_state (const amrex::Array1D<amrex::Real, 1, int_neqs>& y)
{
    bool valid = true;

    for (int n = 1; n <= int_neqs; ++n) {
        const amrex::Real yn = y(n);
        if (yn != yn || std::abs(yn) == std::numeric_limits<amrex::Real>::infinity()) {
            valid = false;
        }
    }

#ifdef STRANG
    for (int n = 1; n <= NumSpec; ++n) {
        if (y(n) < -integrator_rp::atol_spec) {
            valid = false;
        }

        if (! integrator_rp::use_number_densities && y(n) > 1.0_rt) {
            valid = false;
        }
    }

    if (y(net_ienuc) <= 0.0_rt) {
        valid = false;
    }
#endif

    return valid;
}

template <typename BurnT, int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void evaluate_jacobian (BurnT& state, rosenbrock_t<int_neqs>& rstate, const amrex::Real time)
{
    constexpr bool in_jacobian = true;

    rhs(time, state, rstate, rstate.ak1, in_jacobian);
    rstate.n_rhs += 1;

    if (rstate.jacobian_type == 1) {
        jac(time, state, rstate, rstate.jac);
        rstate.n_jac += 1;
        return;
    }

    constexpr amrex::Real UROUND = std::numeric_limits<amrex::Real>::epsilon();
    amrex::Array1D<amrex::Real, 1, int_neqs> ewt;
    amrex::Array1D<amrex::Real, 1, int_neqs> ybase;

    // rhs() may clean or renormalize all components of rstate.y.
    for (int i = 1; i <= int_neqs; ++i) {
        ybase(i) = rstate.y(i);
    }

    amrex::Real fac = 0.0_rt;
    for (int i = 1; i <= int_neqs; ++i) {
        ewt(i) = 1.0_rt / (rtol_for(rstate, i) * std::abs(rstate.y(i)) + atol_for(rstate, i));
        fac += (rstate.ak1(i) * ewt(i)) * (rstate.ak1(i) * ewt(i));
    }
    fac = std::sqrt(fac / static_cast<amrex::Real>(int_neqs));

    amrex::Real R0 = 1000.0_rt * std::abs(rstate.dt) *
        UROUND * static_cast<amrex::Real>(int_neqs) * fac;
    if (R0 == 0.0_rt) {
        R0 = 1.0_rt;
    }

    for (int j = 1; j <= int_neqs; ++j) {
        const amrex::Real yj = ybase(j);
        const amrex::Real R = amrex::max(std::sqrt(UROUND) * std::abs(yj), R0 / ewt(j));
        for (int n = 1; n <= int_neqs; ++n) {
            rstate.y(n) = ybase(n);
        }
        rstate.y(j) += R;

        rhs(time, state, rstate, rstate.work, in_jacobian);

        fac = 1.0_rt / R;
        for (int i = 1; i <= int_neqs; ++i) {
            rstate.jac.set(i, j, (rstate.work(i) - rstate.ak1(i)) * fac);
        }

        for (int n = 1; n <= int_neqs; ++n) {
            rstate.y(n) = ybase(n);
        }
    }

    rstate.n_rhs += int_neqs;
    rstate.n_jac += 1;
}

template <typename BurnT, int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void evaluate_jacobian_warp (BurnT& state, rosenbrock_t<int_neqs>& rstate,
                             const amrex::Real time, const warp::mask_t mask,
                             const int active_lanes)
{
    constexpr bool in_jacobian = true;

    rhs(time, state, rstate, rstate.ak1, in_jacobian);
    rstate.n_rhs += 1;

    if (rstate.jacobian_type == 1) {
        jac(time, state, rstate, rstate.jac);
        rstate.n_jac += 1;
        return;
    }

    constexpr amrex::Real UROUND = std::numeric_limits<amrex::Real>::epsilon();
    amrex::Array1D<amrex::Real, 1, int_neqs> ewt;
    amrex::Array1D<amrex::Real, 1, int_neqs> ybase;

    // rhs() may clean or renormalize all components of rstate.y.
    for (int i = 1; i <= int_neqs; ++i) {
        ybase(i) = rstate.y(i);
    }

    amrex::Real fac_local = 0.0_rt;
    for (int i = 1; i <= int_neqs; ++i) {
        ewt(i) = 1.0_rt / (rtol_for(rstate, i) * std::abs(rstate.y(i)) + atol_for(rstate, i));
        fac_local += (rstate.ak1(i) * ewt(i)) * (rstate.ak1(i) * ewt(i));
    }

    const int coupled_neqs = active_lanes * int_neqs;
    amrex::Real fac = std::sqrt(warp::sum(fac_local, mask) /
                                static_cast<amrex::Real>(coupled_neqs));

    amrex::Real R0 = 1000.0_rt * std::abs(rstate.dt) *
        UROUND * static_cast<amrex::Real>(coupled_neqs) * fac;
    if (R0 == 0.0_rt) {
        R0 = 1.0_rt;
    }

    for (int j = 1; j <= int_neqs; ++j) {
        const amrex::Real yj = ybase(j);
        const amrex::Real R = amrex::max(std::sqrt(UROUND) * std::abs(yj), R0 / ewt(j));
        for (int n = 1; n <= int_neqs; ++n) {
            rstate.y(n) = ybase(n);
        }
        rstate.y(j) += R;

        rhs(time, state, rstate, rstate.work, in_jacobian);

        fac = 1.0_rt / R;
        for (int i = 1; i <= int_neqs; ++i) {
            rstate.jac.set(i, j, (rstate.work(i) - rstate.ak1(i)) * fac);
        }

        for (int n = 1; n <= int_neqs; ++n) {
            rstate.y(n) = ybase(n);
        }
    }

    rstate.n_rhs += int_neqs;
    rstate.n_jac += 1;
}

template <typename BurnT, int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rhs_at (const amrex::Real time, BurnT& state, rosenbrock_t<int_neqs>& rstate,
             const amrex::Array1D<amrex::Real, 1, int_neqs>& y,
             amrex::Array1D<amrex::Real, 1, int_neqs>& ydot)
{
    amrex::Array1D<amrex::Real, 1, int_neqs> ysave;

    for (int n = 1; n <= int_neqs; ++n) {
        ysave(n) = rstate.y(n);
        rstate.y(n) = y(n);
    }

    rhs(time, state, rstate, ydot);

    for (int n = 1; n <= int_neqs; ++n) {
        rstate.y(n) = ysave(n);
    }
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int decompose (rosenbrock_t<int_neqs>& rstate, const amrex::Real fac)
{
    for (int j = 1; j <= int_neqs; ++j) {
        for (int i = 1; i <= int_neqs; ++i) {
            rstate.jac(i, j) = -rstate.jac(i, j);
        }
        rstate.jac(j, j) += fac;
    }

    int ierr_linpack = 0;
    if (integrator_rp::linalg_do_pivoting == 1) {
        constexpr bool allow_pivot{true};
        dgefa<int_neqs, allow_pivot>(rstate.jac, rstate.pivot, ierr_linpack);
    } else {
        constexpr bool allow_pivot{false};
        dgefa<int_neqs, allow_pivot>(rstate.jac, rstate.pivot, ierr_linpack);
    }

    return ierr_linpack;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real select_real (const bool take, const amrex::Real if_true,
                         const amrex::Real if_false)
{
    static_assert(sizeof(amrex::Real) == sizeof(std::uint64_t) ||
                  sizeof(amrex::Real) == sizeof(std::uint32_t),
                  "select_real expects 32-bit or 64-bit amrex::Real");

    if constexpr (sizeof(amrex::Real) == sizeof(std::uint64_t)) {
        union real_bits {
            amrex::Real r;
            std::uint64_t u;
        };

        real_bits true_bits{};
        real_bits false_bits{};
        real_bits result{};
        true_bits.r = if_true;
        false_bits.r = if_false;

        const std::uint64_t mask = std::uint64_t{0} - static_cast<std::uint64_t>(take);
        result.u = (true_bits.u & mask) | (false_bits.u & ~mask);
        return result.r;
    } else {
        union real_bits {
            amrex::Real r;
            std::uint32_t u;
        };

        real_bits true_bits{};
        real_bits false_bits{};
        real_bits result{};
        true_bits.r = if_true;
        false_bits.r = if_false;

        const std::uint32_t mask = std::uint32_t{0} - static_cast<std::uint32_t>(take);
        result.u = (true_bits.u & mask) | (false_bits.u & ~mask);
        return result.r;
    }
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real select_array_value (const RArray1D& a, const int row)
{
    amrex::Real value = 0.0_rt;
    amrex::constexpr_for<1, int_neqs+1>([&] (auto ii) {
        constexpr int i = ii;
        value = select_real(row == i, a(i), value);
    });
    return value;
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void swap_array_values (RArray1D& a, const int l, const int k,
                        const amrex::Real al, const amrex::Real ak,
                        const bool do_swap)
{
    amrex::constexpr_for<1, int_neqs+1>([&] (auto ii) {
        constexpr int i = ii;
        const amrex::Real old_value = a(i);
        amrex::Real new_value = select_real(l == i, ak, old_value);
        new_value = select_real(k == i, al, new_value);
        a(i) = select_real(do_swap, new_value, old_value);
    });
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real select_jac_value (const ArrayUtil::MathArray2D<jac_t, 1, int_neqs, 1, int_neqs>& a,
                              const int row, const int col)
{
    amrex::Real value = 0.0_rt;
    amrex::constexpr_for<1, int_neqs+1>([&] (auto ii) {
        constexpr int i = ii;
        value = select_real(row == i, a(i, col), value);
    });
    return value;
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void swap_jac_values (ArrayUtil::MathArray2D<jac_t, 1, int_neqs, 1, int_neqs>& a,
                      const int l, const int k, const int col,
                      const amrex::Real al, const amrex::Real ak,
                      const bool do_swap)
{
    amrex::constexpr_for<1, int_neqs+1>([&] (auto ii) {
        constexpr int i = ii;
        const amrex::Real old_value = a(i, col);
        amrex::Real new_value = select_real(l == i, ak, old_value);
        new_value = select_real(k == i, al, new_value);
        a(i, col) = select_real(do_swap, new_value, old_value);
    });
}

template <int int_neqs, bool allow_pivot>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int decompose_warp_uniform_impl (rosenbrock_t<int_neqs>& rstate, const amrex::Real fac)
{
    for (int j = 1; j <= int_neqs; ++j) {
        for (int i = 1; i <= int_neqs; ++i) {
            rstate.jac(i, j) = -rstate.jac(i, j);
        }
        rstate.jac(j, j) += fac;
    }

    int info = 0;
    constexpr int nm1 = int_neqs - 1;

    if constexpr (nm1 >= 1) {
        for (int k = 1; k <= nm1; ++k) {

            int l = k;

            if constexpr (allow_pivot) {
                amrex::Real dmax = std::abs(rstate.jac(k, k));
                amrex::constexpr_for<1, int_neqs+1>([&] (auto ii) {
                    constexpr int i = ii;
                    const amrex::Real ai = std::abs(rstate.jac(i, k));
                    const bool take = (i > k) && ai > dmax;
                    l = take ? i : l;
                    dmax = select_real(take, ai, dmax);
                });

                rstate.pivot(k) = static_cast<short>(l);
            }

            const amrex::Real alk = select_jac_value<int_neqs>(rstate.jac, l, k);
            const bool nonzero_pivot = alk != 0.0e0_rt;
            info = nonzero_pivot ? info : k;

            bool do_swap = false;
            if constexpr (allow_pivot) {
                do_swap = nonzero_pivot && (l != k);

                const amrex::Real akk = rstate.jac(k, k);
                swap_jac_values<int_neqs>(rstate.jac, l, k, k, alk, akk, do_swap);
            }

            const amrex::Real safe_pivot =
                select_real(nonzero_pivot, rstate.jac(k, k), 1.0_rt);
            const amrex::Real t_mult =
                select_real(nonzero_pivot, -1.0e0_rt / safe_pivot, 0.0_rt);
            amrex::constexpr_for<1, int_neqs+1>([&] (auto jj) {
                constexpr int j = jj;
                const amrex::Real ajk = rstate.jac(j, k);
                const amrex::Real new_ajk = ajk * t_mult;
                rstate.jac(j, k) =
                    select_real(nonzero_pivot && (j > k), new_ajk, ajk);
            });

            amrex::constexpr_for<1, int_neqs+1>([&] (auto jj) {
                constexpr int j = jj;
                const bool active_j = j > k;
                amrex::Real t = select_jac_value<int_neqs>(rstate.jac, l, j);

                if constexpr (allow_pivot) {
                    const amrex::Real alj = t;
                    const amrex::Real akj = rstate.jac(k, j);
                    swap_jac_values<int_neqs>(rstate.jac, l, k, j, alj, akj,
                                              do_swap && active_j);
                }

                t = select_real(nonzero_pivot && active_j, t, 0.0_rt);
                amrex::constexpr_for<1, int_neqs+1>([&] (auto ii) {
                    constexpr int i = ii;
                    const amrex::Real old_aij = rstate.jac(i, j);
                    const amrex::Real new_aij = old_aij + t * rstate.jac(i, k);
                    rstate.jac(i, j) =
                        select_real((i > k) && active_j, new_aij, old_aij);
                });
            });
        }
    }

    if constexpr (allow_pivot) {
        rstate.pivot(int_neqs) = static_cast<short>(int_neqs);
    }

    info = (rstate.jac(int_neqs, int_neqs) == 0.0e0_rt) ? int_neqs : info;
    return info;
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int decompose_warp_uniform (rosenbrock_t<int_neqs>& rstate, const amrex::Real fac)
{
    if (integrator_rp::linalg_do_pivoting == 1) {
        constexpr bool allow_pivot{true};
        return decompose_warp_uniform_impl<int_neqs, allow_pivot>(rstate, fac);
    }

    constexpr bool allow_pivot{false};
    return decompose_warp_uniform_impl<int_neqs, allow_pivot>(rstate, fac);
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void solve (rosenbrock_t<int_neqs>& rstate, RArray1D& b)
{
    if (integrator_rp::linalg_do_pivoting == 1) {
        constexpr bool allow_pivot{true};
        dgesl<int_neqs, allow_pivot>(rstate.jac, rstate.pivot, b);
    } else {
        constexpr bool allow_pivot{false};
        dgesl<int_neqs, allow_pivot>(rstate.jac, rstate.pivot, b);
    }
}

template <int int_neqs, bool allow_pivot>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void solve_warp_uniform_impl (const rosenbrock_t<int_neqs>& rstate, RArray1D& b)
{
    constexpr int nm1 = int_neqs - 1;

    if constexpr (nm1 >= 1) {
        for (int k = 1; k <= nm1; ++k) {

            amrex::Real t{};
            if constexpr (allow_pivot) {
                const int l = rstate.pivot(k);
                const amrex::Real bl = select_array_value<int_neqs>(b, l);
                const amrex::Real bk = b(k);
                const bool do_swap = l != k;
                swap_array_values<int_neqs>(b, l, k, bl, bk, do_swap);
                t = bl;
            } else {
                t = b(k);
            }

            amrex::constexpr_for<1, int_neqs+1>([&] (auto jj) {
                constexpr int j = jj;
                const amrex::Real old_bj = b(j);
                const amrex::Real new_bj = old_bj + t * rstate.jac(j, k);
                b(j) = select_real(j > k, new_bj, old_bj);
            });
        }
    }

    for (int kb = 1; kb <= int_neqs; ++kb) {
        const int k = int_neqs + 1 - kb;
        b(k) = b(k) / rstate.jac(k, k);
        const amrex::Real t = -b(k);
        amrex::constexpr_for<1, int_neqs+1>([&] (auto jj) {
            constexpr int j = jj;
            const amrex::Real old_bj = b(j);
            const amrex::Real new_bj = old_bj + t * rstate.jac(j, k);
            b(j) = select_real(j < k, new_bj, old_bj);
        });
    }
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void solve_warp_uniform (const rosenbrock_t<int_neqs>& rstate, RArray1D& b)
{
    if (integrator_rp::linalg_do_pivoting == 1) {
        constexpr bool allow_pivot{true};
        solve_warp_uniform_impl<int_neqs, allow_pivot>(rstate, b);
    } else {
        constexpr bool allow_pivot{false};
        solve_warp_uniform_impl<int_neqs, allow_pivot>(rstate, b);
    }
}

} // namespace rosenbrock

template <typename BurnT, typename Tableau>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rosenbrock_integrator (BurnT& state, rosenbrock_t<integrator_neqs<BurnT>()>& rstate)
{
    using C = Tableau;
    constexpr int int_neqs = integrator_neqs<BurnT>();
    static_assert(C::stages <= 8, "Rosenbrock integrator stores up to 8 stages");

    if (rstate.tout == rstate.t) {
        return IERR_SUCCESS;
    }

    for (int n = 1; n <= int_neqs; ++n) {
        if (rosenbrock::atol_for(rstate, n) <= 0.0_rt ||
            rosenbrock::rtol_for(rstate, n) <= 10.0_rt * std::numeric_limits<amrex::Real>::epsilon()) {
            return IERR_TOO_MUCH_ACCURACY_REQUESTED;
        }
    }

    if (integrator_rp::rosenbrock_timestep_controller == 1 &&
        (integrator_rp::yass_epsilon <= 0.0_rt ||
         integrator_rp::yass_floor < 0.0_rt ||
         integrator_rp::yass_fac_min <= 0.0_rt ||
         integrator_rp::yass_fac_max < integrator_rp::yass_fac_min ||
         integrator_rp::yass_reduction_fac <= 0.0_rt)) {
        return IERR_BAD_INPUTS;
    }

    rstate.n_rhs = 0;
    rstate.n_jac = 0;
    rstate.n_step = 0;

    rstate.dt = rstate.tout - rstate.t;

    const amrex::Real posneg = rstate.tout >= rstate.t ? 1.0_rt : -1.0_rt;
    const amrex::Real hmax = amrex::min(integrator_rp::ode_max_dt,
                                        std::abs(rstate.tout - rstate.t));
    amrex::Real h = amrex::min(std::abs(rstate.dt), hmax) * posneg;
    if (std::abs(h) <= 10.0_rt * std::numeric_limits<amrex::Real>::epsilon()) {
        h = 1.e-6_rt * posneg;
    }

    bool reject = false;
    bool last = false;
    int nsing = 0;
    int n_reject = 0;
    amrex::Real facold = 1.0_rt;
    amrex::Real errold = 1.0_rt;
    amrex::Real hopt = h;
    amrex::Real x = rstate.t;
    amrex::Array1D<amrex::Real, 1, int_neqs> rhs_tmp;

    while (true) {
        if (rstate.n_step > integrator_rp::ode_max_steps) {
            rstate.t = x;
            rstate.dt = h;
            return IERR_TOO_MANY_STEPS;
        }

        if (0.1_rt * std::abs(h) <= std::abs(x) * std::numeric_limits<amrex::Real>::epsilon()) {
            rstate.t = x;
            rstate.dt = h;
            return IERR_DT_UNDERFLOW;
        }

        if (last) {
            rstate.t = x;
            rstate.dt = hopt;
            return IERR_SUCCESS;
        }

        hopt = h;
        if ((x + h * (1.0_rt + timestep_safety_factor) - rstate.tout) * posneg >= 0.0_rt) {
            h = rstate.tout - x;
            last = true;
        }

        while (true) {
            rstate.dt = h;

            if (0.1_rt * std::abs(h) <= std::abs(x) * std::numeric_limits<amrex::Real>::epsilon()) {
                rstate.t = x;
                rstate.dt = h;
                return IERR_DT_UNDERFLOW;
            }

            const amrex::Real fac = 1.0_rt / (h * C::gamma);
            rosenbrock::evaluate_jacobian(state, rstate, x);
            int ierr_linpack = rosenbrock::decompose(rstate, fac);

            if (ierr_linpack != 0) {
                nsing += 1;
                if (nsing >= 5) {
                    rstate.t = x;
                    rstate.dt = h;
                    return IERR_LU_DECOMPOSITION_ERROR;
                }
                h *= 0.5_rt;
                reject = true;
                last = false;
                continue;
            }

            rosenbrock::solve(rstate, rstate.ak1);

            bool valid_stage_state = true;
            for (int istage = 2; istage <= C::stages; ++istage) {
                for (int n = 1; n <= int_neqs; ++n) {
                    rstate.ynew(n) = rstate.y(n);
                    rosenbrock::stage(rstate, istage, n) = 0.0_rt;

                    for (int jstage = 1; jstage < istage; ++jstage) {
                        const amrex::Real akj = rosenbrock::stage(rstate, jstage, n);
                        rstate.ynew(n) += C::a(istage, jstage) * akj;
                        rosenbrock::stage(rstate, istage, n) +=
                            (C::c(istage, jstage) / h) * akj;
                    }
                }

                if (! rosenbrock::valid_integrator_state(rstate.ynew)) {
                    h *= 0.25_rt;
                    reject = true;
                    n_reject += 1;
                    last = false;
                    valid_stage_state = false;
                    break;
                }

                rosenbrock::rhs_at(x + C::ct(istage) * h, state, rstate, rstate.ynew, rhs_tmp);
                rstate.n_rhs += 1;

                for (int n = 1; n <= int_neqs; ++n) {
                    rosenbrock::stage(rstate, istage, n) += rhs_tmp(n);
                }
                rosenbrock::solve(rstate, rosenbrock::stage_vector(rstate, istage));
            }

            if (! valid_stage_state) {
                continue;
            }

            for (int n = 1; n <= int_neqs; ++n) {
                amrex::Real y_update = 0.0_rt;
                amrex::Real err = 0.0_rt;

                for (int istage = 1; istage <= C::stages; ++istage) {
                    const amrex::Real aki = rosenbrock::stage(rstate, istage, n);
                    y_update += C::b(istage) * aki;
                    err += C::e(istage) * aki;
                }

                rstate.ynew(n) = rstate.y(n) + y_update;
                rstate.work(n) = err;
            }

            rstate.n_step += 1;

            const bool valid_state = rosenbrock::valid_integrator_state(rstate.ynew);
            if (! valid_state) {
                reject = true;
                n_reject += 1;
                last = false;
                h *= 0.25_rt;
                continue;
            }

            const bool valid_update = rosenbrock::species_change_valid(rstate);
            constexpr amrex::Real err_min = 1.e-10_rt;
            amrex::Real err = rosenbrock::error_norm(rstate);
            amrex::Real controller_fac = 1.0_rt;

            if (integrator_rp::rosenbrock_timestep_controller == 1) {
                err = rosenbrock::yass_change_norm(rstate);
                controller_fac = amrex::Clamp(1.0_rt / amrex::max(err, err_min),
                                              integrator_rp::yass_fac_min,
                                              integrator_rp::yass_fac_max);
            } else {
                controller_fac = amrex::Clamp(
                    std::pow(1.0_rt / amrex::max(err, err_min),
                             1.0_rt / (integrator_rp::h211b_b * integrator_rp::h211b_k)) *
                    std::pow(1.0_rt / amrex::max(errold, err_min),
                             1.0_rt / (integrator_rp::h211b_b * integrator_rp::h211b_k)) *
                    std::pow(facold, -1.0_rt / integrator_rp::h211b_b),
                    integrator_rp::h211b_fac_min, integrator_rp::h211b_fac_max);
            }
            facold = controller_fac;
            errold = err;

            amrex::Real hnew = h * controller_fac;

            if (err <= 1.0_rt && valid_update) {
                for (int n = 1; n <= int_neqs; ++n) {
                    rstate.y(n) = rstate.ynew(n);
                }
                x += h;

                if (std::abs(hnew) > hmax) {
                    hnew = posneg * hmax;
                }
                if (reject) {
                    hnew = posneg * amrex::min(std::abs(hnew), std::abs(h));
                }
                reject = false;
                n_reject = 0;
                h = hnew;
                break;
            }

            reject = true;
            n_reject += 1;
            last = false;
            if (n_reject >= 2) {
                hnew *= integrator_rp::rosenbrock_timestep_controller == 1 ?
                    integrator_rp::yass_reduction_fac : integrator_rp::h211b_reduction_fac;
            }
            if (valid_update) {
                const amrex::Real reject_fac =
                    integrator_rp::rosenbrock_timestep_controller == 1 ?
                    integrator_rp::yass_reduction_fac : integrator_rp::h211b_reduction_fac;
                h = posneg * amrex::min(std::abs(hnew), reject_fac * std::abs(h));
            } else {
                h *= 0.25_rt;
            }
        }
    }
}

template <typename BurnT, typename Tableau>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rosenbrock_integrator_warp_coupled (BurnT& state, rosenbrock_t<integrator_neqs<BurnT>()>& rstate)
{
    using C = Tableau;
    constexpr int int_neqs = integrator_neqs<BurnT>();
    static_assert(C::stages <= 8, "Rosenbrock integrator stores up to 8 stages");

    const auto mask = rosenbrock::warp::active_mask();
    const int active_lanes = rosenbrock::warp::active_count(mask);
    const int leader = rosenbrock::warp::first_lane(mask);

    const amrex::Real t_ref = rosenbrock::warp::broadcast(rstate.t, leader, mask);
    const amrex::Real tout_ref = rosenbrock::warp::broadcast(rstate.tout, leader, mask);
    const bool same_interval = rosenbrock::warp::all(rstate.t == t_ref && rstate.tout == tout_ref, mask);
    if (! same_interval) {
        return IERR_BAD_INPUTS;
    }

    const int jacobian_type_ref = rosenbrock::warp::broadcast(static_cast<int>(rstate.jacobian_type),
                                                             leader, mask);
    if (! rosenbrock::warp::all(static_cast<int>(rstate.jacobian_type) == jacobian_type_ref, mask)) {
        return IERR_BAD_INPUTS;
    }

    if (tout_ref == t_ref) {
        return IERR_SUCCESS;
    }

    bool bad_tolerance = false;
    for (int n = 1; n <= int_neqs; ++n) {
        if (rosenbrock::atol_for(rstate, n) <= 0.0_rt ||
            rosenbrock::rtol_for(rstate, n) <= 10.0_rt * std::numeric_limits<amrex::Real>::epsilon()) {
            bad_tolerance = true;
        }
    }
    if (rosenbrock::warp::any(bad_tolerance, mask)) {
        return IERR_TOO_MUCH_ACCURACY_REQUESTED;
    }

    if (integrator_rp::rosenbrock_timestep_controller == 1 &&
        (integrator_rp::yass_epsilon <= 0.0_rt ||
         integrator_rp::yass_floor < 0.0_rt ||
         integrator_rp::yass_fac_min <= 0.0_rt ||
         integrator_rp::yass_fac_max < integrator_rp::yass_fac_min ||
         integrator_rp::yass_reduction_fac <= 0.0_rt)) {
        return IERR_BAD_INPUTS;
    }

    rstate.n_rhs = 0;
    rstate.n_jac = 0;
    rstate.n_step = 0;

    rstate.dt = tout_ref - t_ref;

    const amrex::Real posneg = tout_ref >= t_ref ? 1.0_rt : -1.0_rt;
    const amrex::Real hmax = amrex::min(integrator_rp::ode_max_dt,
                                        std::abs(tout_ref - t_ref));
    amrex::Real h = amrex::min(std::abs(rstate.dt), hmax) * posneg;
    if (std::abs(h) <= 10.0_rt * std::numeric_limits<amrex::Real>::epsilon()) {
        h = 1.e-6_rt * posneg;
    }

    bool reject = false;
    bool last = false;
    int nsing = 0;
    int n_reject = 0;
    amrex::Real facold = 1.0_rt;
    amrex::Real errold = 1.0_rt;
    amrex::Real hopt = h;
    amrex::Real x = t_ref;
    amrex::Array1D<amrex::Real, 1, int_neqs> rhs_tmp;

    while (true) {
        if (rosenbrock::warp::any(rstate.n_step > integrator_rp::ode_max_steps, mask)) {
            rstate.t = x;
            rstate.dt = h;
            return IERR_TOO_MANY_STEPS;
        }

        if (0.1_rt * std::abs(h) <= std::abs(x) * std::numeric_limits<amrex::Real>::epsilon()) {
            rstate.t = x;
            rstate.dt = h;
            return IERR_DT_UNDERFLOW;
        }

        if (last) {
            rstate.t = x;
            rstate.dt = hopt;
            return IERR_SUCCESS;
        }

        hopt = h;
        if ((x + h * (1.0_rt + timestep_safety_factor) - tout_ref) * posneg >= 0.0_rt) {
            h = tout_ref - x;
            last = true;
        }

        while (true) {
            rstate.dt = h;

            if (0.1_rt * std::abs(h) <= std::abs(x) * std::numeric_limits<amrex::Real>::epsilon()) {
                rstate.t = x;
                rstate.dt = h;
                return IERR_DT_UNDERFLOW;
            }

            // Across the warp, the coupled Jacobian is block diagonal.  Each
            // lane factors its block, while all adaptive decisions below are
            // made from warp-wide reductions over the concatenated system.
            const amrex::Real fac = 1.0_rt / (h * C::gamma);
            rosenbrock::evaluate_jacobian_warp(state, rstate, x, mask, active_lanes);
            int ierr_linpack = rosenbrock::decompose_warp_uniform(rstate, fac);

            if (rosenbrock::warp::any(ierr_linpack != 0, mask)) {
                nsing += 1;
                if (nsing >= 5) {
                    rstate.t = x;
                    rstate.dt = h;
                    return IERR_LU_DECOMPOSITION_ERROR;
                }
                h *= 0.5_rt;
                reject = true;
                last = false;
                continue;
            }

            rosenbrock::solve_warp_uniform(rstate, rstate.ak1);

            bool valid_stage_state = true;
            for (int istage = 2; istage <= C::stages; ++istage) {
                for (int n = 1; n <= int_neqs; ++n) {
                    rstate.ynew(n) = rstate.y(n);
                    rosenbrock::stage(rstate, istage, n) = 0.0_rt;

                    for (int jstage = 1; jstage < istage; ++jstage) {
                        const amrex::Real akj = rosenbrock::stage(rstate, jstage, n);
                        rstate.ynew(n) += C::a(istage, jstage) * akj;
                        rosenbrock::stage(rstate, istage, n) +=
                            (C::c(istage, jstage) / h) * akj;
                    }
                }

                const bool local_valid_stage = rosenbrock::valid_integrator_state(rstate.ynew);
                if (! rosenbrock::warp::all(local_valid_stage, mask)) {
                    h *= 0.25_rt;
                    reject = true;
                    n_reject += 1;
                    last = false;
                    valid_stage_state = false;
                    break;
                }

                rosenbrock::rhs_at(x + C::ct(istage) * h, state, rstate, rstate.ynew, rhs_tmp);
                rstate.n_rhs += 1;

                for (int n = 1; n <= int_neqs; ++n) {
                    rosenbrock::stage(rstate, istage, n) += rhs_tmp(n);
                }
                rosenbrock::solve_warp_uniform(rstate, rosenbrock::stage_vector(rstate, istage));
            }

            if (! valid_stage_state) {
                continue;
            }

            for (int n = 1; n <= int_neqs; ++n) {
                amrex::Real y_update = 0.0_rt;
                amrex::Real err = 0.0_rt;

                for (int istage = 1; istage <= C::stages; ++istage) {
                    const amrex::Real aki = rosenbrock::stage(rstate, istage, n);
                    y_update += C::b(istage) * aki;
                    err += C::e(istage) * aki;
                }

                rstate.ynew(n) = rstate.y(n) + y_update;
                rstate.work(n) = err;
            }

            rstate.n_step += 1;

            const bool valid_state = rosenbrock::warp::all(rosenbrock::valid_integrator_state(rstate.ynew), mask);
            if (! valid_state) {
                reject = true;
                n_reject += 1;
                last = false;
                h *= 0.25_rt;
                continue;
            }

            const bool valid_update = rosenbrock::warp::all(rosenbrock::species_change_valid(rstate), mask);
            constexpr amrex::Real err_min = 1.e-10_rt;
            amrex::Real err = rosenbrock::error_norm_warp(rstate, mask, active_lanes);
            amrex::Real controller_fac = 1.0_rt;

            if (integrator_rp::rosenbrock_timestep_controller == 1) {
                err = rosenbrock::yass_change_norm_warp(rstate, mask);
                controller_fac = amrex::Clamp(1.0_rt / amrex::max(err, err_min),
                                              integrator_rp::yass_fac_min,
                                              integrator_rp::yass_fac_max);
            } else {
                controller_fac = amrex::Clamp(
                    std::pow(1.0_rt / amrex::max(err, err_min),
                             1.0_rt / (integrator_rp::h211b_b * integrator_rp::h211b_k)) *
                    std::pow(1.0_rt / amrex::max(errold, err_min),
                             1.0_rt / (integrator_rp::h211b_b * integrator_rp::h211b_k)) *
                    std::pow(facold, -1.0_rt / integrator_rp::h211b_b),
                    integrator_rp::h211b_fac_min, integrator_rp::h211b_fac_max);
            }
            facold = controller_fac;
            errold = err;

            amrex::Real hnew = h * controller_fac;

            if (err <= 1.0_rt && valid_update) {
                for (int n = 1; n <= int_neqs; ++n) {
                    rstate.y(n) = rstate.ynew(n);
                }
                x += h;

                if (std::abs(hnew) > hmax) {
                    hnew = posneg * hmax;
                }
                if (reject) {
                    hnew = posneg * amrex::min(std::abs(hnew), std::abs(h));
                }
                reject = false;
                n_reject = 0;
                h = hnew;
                break;
            }

            reject = true;
            n_reject += 1;
            last = false;
            if (n_reject >= 2) {
                hnew *= integrator_rp::rosenbrock_timestep_controller == 1 ?
                    integrator_rp::yass_reduction_fac : integrator_rp::h211b_reduction_fac;
            }
            if (valid_update) {
                const amrex::Real reject_fac =
                    integrator_rp::rosenbrock_timestep_controller == 1 ?
                    integrator_rp::yass_reduction_fac : integrator_rp::h211b_reduction_fac;
                h = posneg * amrex::min(std::abs(hnew), reject_fac * std::abs(h));
            } else {
                h *= 0.25_rt;
            }
        }
    }
}

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rosenbrock_integrator_warp_coupled (BurnT& state, rosenbrock_t<integrator_neqs<BurnT>()>& rstate)
{
    if (integrator_rp::rosenbrock_tableau == 1) {
        return rosenbrock_integrator_warp_coupled<BurnT, rosenbrock::rodas4p_tableau>(state, rstate);
    }
    if (integrator_rp::rosenbrock_tableau == 2) {
        return rosenbrock_integrator_warp_coupled<BurnT, rosenbrock::rodas3p_tableau>(state, rstate);
    }
    if (integrator_rp::rosenbrock_tableau == 3) {
        return rosenbrock_integrator_warp_coupled<BurnT, rosenbrock::ros2s_tableau>(state, rstate);
    }
    if (integrator_rp::rosenbrock_tableau == 4) {
        return rosenbrock_integrator_warp_coupled<BurnT, rosenbrock::ros2_tableau>(state, rstate);
    }
    if (integrator_rp::rosenbrock_tableau == 5) {
        return rosenbrock_integrator_warp_coupled<BurnT, rosenbrock::rosenbrock_euler_tableau>(state, rstate);
    }

    return rosenbrock_integrator_warp_coupled<BurnT, rosenbrock::rodas5p_tableau>(state, rstate);
}

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rosenbrock_integrator (BurnT& state, rosenbrock_t<integrator_neqs<BurnT>()>& rstate)
{
    if (integrator_rp::rosenbrock_tableau == 1) {
        return rosenbrock_integrator<BurnT, rosenbrock::rodas4p_tableau>(state, rstate);
    }
    if (integrator_rp::rosenbrock_tableau == 2) {
        return rosenbrock_integrator<BurnT, rosenbrock::rodas3p_tableau>(state, rstate);
    }
    if (integrator_rp::rosenbrock_tableau == 3) {
        return rosenbrock_integrator<BurnT, rosenbrock::ros2s_tableau>(state, rstate);
    }
    if (integrator_rp::rosenbrock_tableau == 4) {
        return rosenbrock_integrator<BurnT, rosenbrock::ros2_tableau>(state, rstate);
    }
    if (integrator_rp::rosenbrock_tableau == 5) {
        return rosenbrock_integrator<BurnT, rosenbrock::rosenbrock_euler_tableau>(state, rstate);
    }

    return rosenbrock_integrator<BurnT, rosenbrock::rodas5p_tableau>(state, rstate);
}

#endif
