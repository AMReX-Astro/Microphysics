#ifndef ROS2S_INTEGRATOR_H
#define ROS2S_INTEGRATOR_H

#include <AMReX_Algorithm.H>
#ifndef AMREX_USE_GPU
#include <AMReX_ParallelDescriptor.H>
#endif

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
#include <ros2s_tableau.H>

namespace ros2s {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real profile_second ()
{
#ifdef AMREX_USE_GPU
    return 0.0_rt;
#else
    return amrex::ParallelDescriptor::second();
#endif
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real& stage (ros2s_t<int_neqs>& ros2s, const int i, const int n)
{
    if (i == 1) {
        return ros2s.ak1(n);
    }
    if (i == 2) {
        return ros2s.ak2(n);
    }
    if (i == 3) {
        return ros2s.ak3(n);
    }
    if (i == 4) {
        return ros2s.ak4(n);
    }
    return ros2s.work(n);
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real stage (const ros2s_t<int_neqs>& ros2s, const int i, const int n)
{
    if (i == 1) {
        return ros2s.ak1(n);
    }
    if (i == 2) {
        return ros2s.ak2(n);
    }
    if (i == 3) {
        return ros2s.ak3(n);
    }
    if (i == 4) {
        return ros2s.ak4(n);
    }
    return ros2s.work(n);
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
RArray1D& stage_vector (ros2s_t<int_neqs>& ros2s, const int i)
{
    if (i == 1) {
        return ros2s.ak1;
    }
    if (i == 2) {
        return ros2s.ak2;
    }
    if (i == 3) {
        return ros2s.ak3;
    }
    if (i == 4) {
        return ros2s.ak4;
    }
    return ros2s.work;
}

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

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool species_change_valid (const ros2s_t<int_neqs>& ros2s)
{
#ifdef STRANG
    constexpr amrex::Real increase_change_factor = 4.0_rt;
    constexpr amrex::Real decrease_change_factor = 0.25_rt;

    for (int i = 1; i <= NumSpec; ++i) {
        if (std::abs(ros2s.y(i)) > integrator_rp::X_reject_buffer * ros2s.atol_spec &&
            std::abs(ros2s.ynew(i)) > integrator_rp::X_reject_buffer * ros2s.atol_spec &&
            (std::abs(ros2s.ynew(i)) > increase_change_factor * std::abs(ros2s.y(i)) ||
             std::abs(ros2s.ynew(i)) < decrease_change_factor * std::abs(ros2s.y(i)))) {
            return false;
        }
    }
#endif

    return true;
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real yass_change_norm (const ros2s_t<int_neqs>& ros2s)
{
    amrex::Real max_change = 0.0_rt;
    amrex::Real max_species = 0.0_rt;

    for (int i = 1; i <= NumSpec && i <= int_neqs; ++i) {
        max_species = amrex::max(max_species, std::abs(ros2s.y(i)));
    }

    const amrex::Real active_floor = integrator_rp::yass_floor * max_species;

    for (int i = 1; i <= NumSpec && i <= int_neqs; ++i) {
        const amrex::Real yi = std::abs(ros2s.y(i));
        if (yi > active_floor) {
            max_change = amrex::max(max_change, std::abs(ros2s.ynew(i) - ros2s.y(i)) / yi);
        }
    }

    return max_change / integrator_rp::yass_epsilon;
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool valid_integrator_state (const amrex::Array1D<amrex::Real, 1, int_neqs>& y)
{
    for (int n = 1; n <= int_neqs; ++n) {
        const amrex::Real yn = y(n);
        if (yn != yn || std::abs(yn) == std::numeric_limits<amrex::Real>::infinity()) {
            return false;
        }
    }

#ifdef STRANG
    for (int n = 1; n <= NumSpec; ++n) {
        if (y(n) < 0.0_rt) {
            return false;
        }

        if (! integrator_rp::use_number_densities && y(n) > 1.0_rt) {
            return false;
        }
    }

    if (y(net_ienuc) <= 0.0_rt) {
        return false;
    }
#endif

    return true;
}

template <typename BurnT, int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void evaluate_jacobian (BurnT& state, ros2s_t<int_neqs>& ros2s, const amrex::Real time)
{
    constexpr bool in_jacobian = true;

    amrex::Real profile_start = profile_second();
    rhs(time, state, ros2s, ros2s.ak1, in_jacobian);
    ros2s.rhs_walltime += profile_second() - profile_start;
    ros2s.n_rhs += 1;

    if (ros2s.jacobian_type == 1) {
        profile_start = profile_second();
        jac(time, state, ros2s, ros2s.jac);
        ros2s.jacobian_walltime += profile_second() - profile_start;
        ros2s.n_jac += 1;
        return;
    }

    profile_start = profile_second();
    constexpr amrex::Real UROUND = std::numeric_limits<amrex::Real>::epsilon();
    amrex::Array1D<amrex::Real, 1, int_neqs> ewt;

    amrex::Real fac = 0.0_rt;
    for (int i = 1; i <= int_neqs; ++i) {
        ewt(i) = 1.0_rt / (rtol_for(ros2s, i) * std::abs(ros2s.y(i)) + atol_for(ros2s, i));
        fac += (ros2s.ak1(i) * ewt(i)) * (ros2s.ak1(i) * ewt(i));
    }
    fac = std::sqrt(fac / static_cast<amrex::Real>(int_neqs));

    amrex::Real R0 = 1000.0_rt * std::abs(ros2s.dt) *
        UROUND * static_cast<amrex::Real>(int_neqs) * fac;
    if (R0 == 0.0_rt) {
        R0 = 1.0_rt;
    }
    ros2s.jacobian_walltime += profile_second() - profile_start;

    for (int j = 1; j <= int_neqs; ++j) {
        profile_start = profile_second();
        const amrex::Real yj = ros2s.y(j);
        const amrex::Real R = amrex::max(std::sqrt(UROUND) * std::abs(yj), R0 / ewt(j));
        ros2s.y(j) += R;
        ros2s.jacobian_walltime += profile_second() - profile_start;

        const amrex::Real rhs_profile_start = profile_second();
        rhs(time, state, ros2s, ros2s.work, in_jacobian);
        ros2s.rhs_walltime += profile_second() - rhs_profile_start;

        profile_start = profile_second();
        fac = 1.0_rt / R;
        for (int i = 1; i <= int_neqs; ++i) {
            ros2s.jac.set(i, j, (ros2s.work(i) - ros2s.ak1(i)) * fac);
        }

        ros2s.y(j) = yj;
        ros2s.jacobian_walltime += profile_second() - profile_start;
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
    amrex::Array1D<amrex::Real, 1, int_neqs> ysave;

    for (int n = 1; n <= int_neqs; ++n) {
        ysave(n) = ros2s.y(n);
        ros2s.y(n) = y(n);
    }

    const amrex::Real profile_start = profile_second();
    rhs(time, state, ros2s, ydot);
    ros2s.rhs_walltime += profile_second() - profile_start;

    for (int n = 1; n <= int_neqs; ++n) {
        ros2s.y(n) = ysave(n);
    }
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int decompose (ros2s_t<int_neqs>& ros2s, const amrex::Real fac)
{
    const amrex::Real profile_start = profile_second();

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

    ros2s.lu_walltime += profile_second() - profile_start;

    return ierr_linpack;
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void solve (ros2s_t<int_neqs>& ros2s, RArray1D& b)
{
    const amrex::Real profile_start = profile_second();

    if (integrator_rp::linalg_do_pivoting == 1) {
        constexpr bool allow_pivot{true};
        dgesl<int_neqs, allow_pivot>(ros2s.jac, ros2s.pivot, b);
    } else {
        constexpr bool allow_pivot{false};
        dgesl<int_neqs, allow_pivot>(ros2s.jac, ros2s.pivot, b);
    }

    ros2s.linear_solve_walltime += profile_second() - profile_start;
}

} // namespace ros2s

template <typename BurnT, typename Tableau>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rosenbrock_integrator (BurnT& state, ros2s_t<integrator_neqs<BurnT>()>& ros2s)
{
    using C = Tableau;
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

    if (integrator_rp::rosenbrock_timestep_controller == 1 &&
        (integrator_rp::yass_epsilon <= 0.0_rt ||
         integrator_rp::yass_floor < 0.0_rt ||
         integrator_rp::yass_fac_min <= 0.0_rt ||
         integrator_rp::yass_fac_max < integrator_rp::yass_fac_min ||
         integrator_rp::yass_reduction_fac <= 0.0_rt)) {
        return IERR_BAD_INPUTS;
    }

    ros2s.n_rhs = 0;
    ros2s.n_jac = 0;
    ros2s.n_step = 0;
    ros2s.rhs_walltime = 0.0_rt;
    ros2s.jacobian_walltime = 0.0_rt;
    ros2s.lu_walltime = 0.0_rt;
    ros2s.linear_solve_walltime = 0.0_rt;

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
    int n_reject = 0;
    amrex::Real facold = 1.0_rt;
    amrex::Real errold = 1.0_rt;
    amrex::Real hopt = h;
    amrex::Real x = ros2s.t;
    amrex::Array1D<amrex::Real, 1, int_neqs> rhs_tmp;

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

            if (0.1_rt * std::abs(h) <= std::abs(x) * std::numeric_limits<amrex::Real>::epsilon()) {
                ros2s.t = x;
                ros2s.dt = h;
                return IERR_DT_UNDERFLOW;
            }

            const amrex::Real fac = 1.0_rt / (h * C::gamma);
            ros2s::evaluate_jacobian(state, ros2s, x);
            int ierr_linpack = ros2s::decompose(ros2s, fac);

            if (ierr_linpack != 0) {
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

            bool valid_stage_state = true;
            for (int istage = 2; istage <= C::stages; ++istage) {
                for (int n = 1; n <= int_neqs; ++n) {
                    ros2s.ynew(n) = ros2s.y(n);
                    ros2s::stage(ros2s, istage, n) = 0.0_rt;

                    for (int jstage = 1; jstage < istage; ++jstage) {
                        const amrex::Real akj = ros2s::stage(ros2s, jstage, n);
                        ros2s.ynew(n) += C::a(istage, jstage) * akj;
                        ros2s::stage(ros2s, istage, n) +=
                            (C::c(istage, jstage) / h) * akj;
                    }
                }

                if (! ros2s::valid_integrator_state(ros2s.ynew)) {
                    h *= 0.25_rt;
                    reject = true;
                    n_reject += 1;
                    last = false;
                    valid_stage_state = false;
                    break;
                }

                ros2s::rhs_at(x + C::ct(istage) * h, state, ros2s, ros2s.ynew, rhs_tmp);
                ros2s.n_rhs += 1;

                for (int n = 1; n <= int_neqs; ++n) {
                    ros2s::stage(ros2s, istage, n) += rhs_tmp(n);
                }
                ros2s::solve(ros2s, ros2s::stage_vector(ros2s, istage));
            }

            if (! valid_stage_state) {
                continue;
            }

            for (int n = 1; n <= int_neqs; ++n) {
                amrex::Real y_update = 0.0_rt;
                amrex::Real err = 0.0_rt;

                for (int istage = 1; istage <= C::stages; ++istage) {
                    const amrex::Real aki = ros2s::stage(ros2s, istage, n);
                    y_update += C::b(istage) * aki;
                    err += C::e(istage) * aki;
                }

                ros2s.ynew(n) = ros2s.y(n) + y_update;
                ros2s.work(n) = err;
            }

            ros2s.n_step += 1;

            const bool valid_state = ros2s::valid_integrator_state(ros2s.ynew);
            if (! valid_state) {
                reject = true;
                n_reject += 1;
                last = false;
                h *= 0.25_rt;
                continue;
            }

            const bool valid_update = ros2s::species_change_valid(ros2s);
            constexpr amrex::Real err_min = 1.e-10_rt;
            amrex::Real err = ros2s::error_norm(ros2s);
            amrex::Real controller_fac = 1.0_rt;

            if (integrator_rp::rosenbrock_timestep_controller == 1) {
                err = ros2s::yass_change_norm(ros2s);
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
            h = valid_update ? hnew : 0.25_rt * h;
        }
    }
}

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int ros2s_integrator (BurnT& state, ros2s_t<integrator_neqs<BurnT>()>& ros2s)
{
    if (integrator_rp::rosenbrock_tableau == 1) {
        return rosenbrock_integrator<BurnT, ros2s::ros2_tableau>(state, ros2s);
    }
    if (integrator_rp::rosenbrock_tableau == 2) {
        return rosenbrock_integrator<BurnT, ros2s::sandu_a_tableau>(state, ros2s);
    }
    if (integrator_rp::rosenbrock_tableau == 3) {
        return rosenbrock_integrator<BurnT, ros2s::sandu_b_tableau>(state, ros2s);
    }
    if (integrator_rp::rosenbrock_tableau == 4) {
        return rosenbrock_integrator<BurnT, ros2s::sandu_c_tableau>(state, ros2s);
    }
    if (integrator_rp::rosenbrock_tableau == 5) {
        return rosenbrock_integrator<BurnT, ros2s::sandu_d_tableau>(state, ros2s);
    }
    if (integrator_rp::rosenbrock_tableau == 6) {
        return rosenbrock_integrator<BurnT, ros2s::rosenbrock_euler_tableau>(state, ros2s);
    }

    return rosenbrock_integrator<BurnT, ros2s::ros2s_tableau>(state, ros2s);
}

#endif
