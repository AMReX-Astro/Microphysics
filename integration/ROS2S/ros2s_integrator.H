#ifndef ROS2S_INTEGRATOR_H
#define ROS2S_INTEGRATOR_H

#include <AMReX_Algorithm.H>

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
#include <ros2s_type.H>

namespace ros2s {

struct coefficients {
    static constexpr amrex::Real gamma = 0.292893218813452_rt;
    static constexpr amrex::Real ct2 = 0.585786437626905_rt;
    static constexpr amrex::Real a21 = 2.0000000000000036_rt;
    static constexpr amrex::Real a31 = 6.828427124746214_rt;
    static constexpr amrex::Real a32 = 3.4142135623731007_rt;
    static constexpr amrex::Real c21 = -6.828427124746214_rt;
    static constexpr amrex::Real c31 = -10.949747468305889_rt;
    static constexpr amrex::Real c32 = -7.535533905932761_rt;
    static constexpr amrex::Real b1 = 6.828427124746214_rt;
    static constexpr amrex::Real b2 = 3.414213562373101_rt;
    static constexpr amrex::Real b3 = 1.0_rt;
    static constexpr amrex::Real e1 = -0.23570226039551292_rt;
    static constexpr amrex::Real e2 = -0.23570226039551567_rt;
    static constexpr amrex::Real e3 = -0.13807118745769906_rt;
};

template <typename Ros2sT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real rtol_for (const Ros2sT& ros2s, const int n)
{
    return n == net_ienuc ? ros2s.rtol_enuc : ros2s.rtol_spec;
}

template <typename Ros2sT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real atol_for (const Ros2sT& ros2s, const int n)
{
    return n == net_ienuc ? ros2s.atol_enuc : ros2s.atol_spec;
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real error_norm (const ros2s_t<int_neqs>& ros2s)
{
    amrex::Real err = 0.0_rt;
    for (int n = 1; n <= int_neqs; ++n) {
        const amrex::Real sk = atol_for(ros2s, n) +
            rtol_for(ros2s, n) * amrex::max(std::abs(ros2s.y(n)), std::abs(ros2s.ynew(n)));
        const amrex::Real term = ros2s.work(n) / sk;
        err += term * term;
    }
    return std::sqrt(err / static_cast<amrex::Real>(int_neqs));
}

template <typename BurnT, int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void evaluate_jacobian (BurnT& state, ros2s_t<int_neqs>& ros2s, const amrex::Real time)
{
    constexpr bool in_jacobian = true;

    rhs(time, state, ros2s, ros2s.ak1, in_jacobian);
    ros2s.n_rhs += 1;

    if (ros2s.jacobian_type == 1) {
        jac(time, state, ros2s, ros2s.jac);
        ros2s.n_jac += 1;
        return;
    }

    amrex::Real fac = 0.0_rt;
    for (int n = 1; n <= int_neqs; ++n) {
        const amrex::Real weight = 1.0_rt /
            (rtol_for(ros2s, n) * std::abs(ros2s.y(n)) + atol_for(ros2s, n));
        fac += (ros2s.ak1(n) * weight) * (ros2s.ak1(n) * weight);
    }
    fac = std::sqrt(fac / static_cast<amrex::Real>(int_neqs));

    amrex::Real r0 = 1000.0_rt * std::abs(ros2s.dt) *
        std::numeric_limits<amrex::Real>::epsilon() *
        static_cast<amrex::Real>(int_neqs) * fac;
    if (r0 == 0.0_rt) {
        r0 = 1.0_rt;
    }

    for (int j = 1; j <= int_neqs; ++j) {
        const amrex::Real yj = ros2s.y(j);
        const amrex::Real weight = rtol_for(ros2s, j) * std::abs(yj) + atol_for(ros2s, j);
        const amrex::Real dy = amrex::max(std::sqrt(std::numeric_limits<amrex::Real>::epsilon()) *
                                          std::abs(yj), r0 * weight);
        ros2s.y(j) += dy;

        rhs(time, state, ros2s, ros2s.work, in_jacobian);
        const amrex::Real dyinv = 1.0_rt / dy;
        for (int i = 1; i <= int_neqs; ++i) {
            ros2s.jac(i, j) = (ros2s.work(i) - ros2s.ak1(i)) * dyinv;
        }

        ros2s.y(j) = yj;
    }

    ros2s.n_rhs += int_neqs;
    ros2s.n_jac += 1;
}

template <typename BurnT, int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rhs_at (const amrex::Real time, BurnT& state, ros2s_t<int_neqs>& ros2s,
             const amrex::Array1D<amrex::Real, 1, int_neqs>& y,
             amrex::Array1D<amrex::Real, 1, int_neqs>& ydot)
{
    for (int n = 1; n <= int_neqs; ++n) {
        ros2s.ysave(n) = ros2s.y(n);
        ros2s.y(n) = y(n);
    }

    rhs(time, state, ros2s, ydot);

    for (int n = 1; n <= int_neqs; ++n) {
        ros2s.y(n) = ros2s.ysave(n);
    }
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int decompose (ros2s_t<int_neqs>& ros2s, const amrex::Real fac)
{
    for (int j = 1; j <= int_neqs; ++j) {
        for (int i = 1; i <= int_neqs; ++i) {
            ros2s.jac(i, j) = -ros2s.jac(i, j);
        }
        ros2s.jac(j, j) += fac;
    }

    int ierr_linpack = 0;
    if (integrator_rp::linalg_do_pivoting == 1) {
        constexpr bool allow_pivot{true};
        dgefa<int_neqs, allow_pivot>(ros2s.jac, ros2s.pivot, ierr_linpack);
    } else {
        constexpr bool allow_pivot{false};
        dgefa<int_neqs, allow_pivot>(ros2s.jac, ros2s.pivot, ierr_linpack);
    }

    return ierr_linpack;
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void solve (ros2s_t<int_neqs>& ros2s, RArray1D& b)
{
    if (integrator_rp::linalg_do_pivoting == 1) {
        constexpr bool allow_pivot{true};
        dgesl<int_neqs, allow_pivot>(ros2s.jac, ros2s.pivot, b);
    } else {
        constexpr bool allow_pivot{false};
        dgesl<int_neqs, allow_pivot>(ros2s.jac, ros2s.pivot, b);
    }
}

} // namespace ros2s

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int ros2s_integrator (BurnT& state, ros2s_t<integrator_neqs<BurnT>()>& ros2s)
{
    using C = ros2s::coefficients;
    constexpr int int_neqs = integrator_neqs<BurnT>();

    if (ros2s.tout == ros2s.t) {
        return IERR_SUCCESS;
    }

    for (int n = 1; n <= int_neqs; ++n) {
        if (ros2s::atol_for(ros2s, n) <= 0.0_rt ||
            ros2s::rtol_for(ros2s, n) <= 10.0_rt * std::numeric_limits<amrex::Real>::epsilon()) {
            return IERR_TOO_MUCH_ACCURACY_REQUESTED;
        }
    }

    ros2s.n_rhs = 0;
    ros2s.n_jac = 0;
    ros2s.n_step = 0;
    ros2s.n_accept = 0;
    ros2s.n_reject = 0;

    rhs(ros2s.t, state, ros2s, ros2s.ak1);
    ros2s.n_rhs += 1;
    ros2s.dt = ros2s.tout - ros2s.t;

    const amrex::Real posneg = ros2s.tout >= ros2s.t ? 1.0_rt : -1.0_rt;
    const amrex::Real hmax = amrex::min(integrator_rp::ode_max_dt,
                                        std::abs(ros2s.tout - ros2s.t));
    amrex::Real h = amrex::min(std::abs(ros2s.dt), hmax) * posneg;
    if (std::abs(h) <= 10.0_rt * std::numeric_limits<amrex::Real>::epsilon()) {
        h = 1.e-6_rt * posneg;
    }

    bool reject = false;
    bool last = false;
    int nsing = 0;
    amrex::Real hacc = 0.0_rt;
    amrex::Real erracc = 1.0_rt;
    amrex::Real hopt = h;
    amrex::Real x = ros2s.t;

    while (true) {
        if (ros2s.n_step > integrator_rp::ode_max_steps) {
            ros2s.t = x;
            ros2s.dt = h;
            return IERR_TOO_MANY_STEPS;
        }

        if (0.1_rt * std::abs(h) <= std::abs(x) * std::numeric_limits<amrex::Real>::epsilon()) {
            ros2s.t = x;
            ros2s.dt = h;
            return IERR_DT_UNDERFLOW;
        }

        if (last) {
            ros2s.t = x;
            ros2s.dt = hopt;
            return IERR_SUCCESS;
        }

        hopt = h;
        if ((x + h * (1.0_rt + timestep_safety_factor) - ros2s.tout) * posneg >= 0.0_rt) {
            h = ros2s.tout - x;
            last = true;
        }

        while (true) {
            ros2s.dt = h;
            ros2s::evaluate_jacobian(state, ros2s, x);

            const amrex::Real fac = 1.0_rt / (h * C::gamma);
            if (ros2s::decompose(ros2s, fac) != 0) {
                nsing += 1;
                if (nsing >= 5) {
                    ros2s.t = x;
                    ros2s.dt = h;
                    return IERR_LU_DECOMPOSITION_ERROR;
                }
                h *= 0.5_rt;
                reject = true;
                last = false;
                continue;
            }

            ros2s::solve(ros2s, ros2s.ak1);

            for (int n = 1; n <= int_neqs; ++n) {
                ros2s.ynew(n) = ros2s.y(n) + C::a21 * ros2s.ak1(n);
                ros2s.ak2(n) = (C::c21 / h) * ros2s.ak1(n);
            }
            ros2s::rhs_at(x + C::ct2 * h, state, ros2s, ros2s.ynew, ros2s.rhs_tmp);
            ros2s.n_rhs += 1;
            for (int n = 1; n <= int_neqs; ++n) {
                ros2s.ak2(n) += ros2s.rhs_tmp(n);
            }
            ros2s::solve(ros2s, ros2s.ak2);

            for (int n = 1; n <= int_neqs; ++n) {
                ros2s.ynew(n) = ros2s.y(n) + C::a31 * ros2s.ak1(n) + C::a32 * ros2s.ak2(n);
                ros2s.work(n) = (C::c31 * ros2s.ak1(n) + C::c32 * ros2s.ak2(n)) / h;
            }
            ros2s::rhs_at(x + h, state, ros2s, ros2s.ynew, ros2s.rhs_tmp);
            ros2s.n_rhs += 1;
            for (int n = 1; n <= int_neqs; ++n) {
                ros2s.work(n) += ros2s.rhs_tmp(n);
            }
            ros2s::solve(ros2s, ros2s.work);

            for (int n = 1; n <= int_neqs; ++n) {
                const amrex::Real ak3 = ros2s.work(n);
                ros2s.ynew(n) = ros2s.y(n) + C::b1 * ros2s.ak1(n) +
                    C::b2 * ros2s.ak2(n) + C::b3 * ak3;
                ros2s.work(n) = C::e1 * ros2s.ak1(n) + C::e2 * ros2s.ak2(n) + C::e3 * ak3;
            }

            ros2s.n_step += 1;

            const amrex::Real err = ros2s::error_norm(ros2s);
            const amrex::Real safe = 0.9_rt;
            const amrex::Real fac_min = 0.2_rt;
            const amrex::Real fac_max = 6.0_rt;
            const amrex::Real raw_fac = std::cbrt(err) / safe;
            const amrex::Real lower_fac = 1.0_rt / fac_max;
            const amrex::Real upper_fac = 1.0_rt / fac_min;
            const amrex::Real fac_step = amrex::Clamp(raw_fac, lower_fac, upper_fac);
            amrex::Real hnew = h / fac_step;

            if (err <= 1.0_rt) {
                ros2s.n_accept += 1;
                if (ros2s.n_accept > 1) {
                    const amrex::Real facgus = amrex::Clamp((hacc / h) *
                                                            std::cbrt((err * err) / erracc) /
                                                            safe,
                                                            lower_fac, upper_fac);
                    hnew = h / amrex::max(fac_step, facgus);
                }
                hacc = h;
                erracc = amrex::max(1.e-2_rt, err);

                for (int n = 1; n <= int_neqs; ++n) {
                    ros2s.y(n) = ros2s.ynew(n);
                }
                x += h;

                if (std::abs(hnew) > hmax) {
                    hnew = posneg * hmax;
                }
                if (reject) {
                    hnew = posneg * amrex::min(std::abs(hnew), std::abs(h));
                }
                reject = false;
                h = hnew;
                break;
            }

            reject = true;
            last = false;
            h = hnew;
            if (ros2s.n_accept >= 1) {
                ros2s.n_reject += 1;
            }
        }
    }
}

#endif
