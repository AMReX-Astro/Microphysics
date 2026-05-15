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

struct ros2s_tableau {
    static constexpr int stages = 3;
    static constexpr amrex::Real gamma = 0.292893218813452_rt;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real ct (const int i) {
        return i == 2 ? 0.585786437626905_rt : 1.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real a (const int i, const int j) {
        if (i == 2 && j == 1) {
            return 2.0000000000000036_rt;
        }
        if (i == 3 && j == 1) {
            return 6.828427124746214_rt;
        }
        if (i == 3 && j == 2) {
            return 3.4142135623731007_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real c (const int i, const int j) {
        if (i == 2 && j == 1) {
            return -6.828427124746214_rt;
        }
        if (i == 3 && j == 1) {
            return -10.949747468305889_rt;
        }
        if (i == 3 && j == 2) {
            return -7.535533905932761_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real b (const int i) {
        if (i == 1) {
            return 6.828427124746214_rt;
        }
        if (i == 2) {
            return 3.414213562373101_rt;
        }
        return 1.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real e (const int i) {
        if (i == 1) {
            return -0.23570226039551292_rt;
        }
        if (i == 2) {
            return -0.23570226039551567_rt;
        }
        return -0.13807118745769906_rt;
    }
};

struct ros2_tableau {
    static constexpr int stages = 2;
    static constexpr amrex::Real gamma = 1.7071067811865475_rt;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real ct (const int i) {
        (void) i;
        return 1.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real a (const int i, const int j) {
        if (i == 2 && j == 1) {
            return 0.585786437626905_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real c (const int i, const int j) {
        if (i == 2 && j == 1) {
            return -1.1715728752538102_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real b (const int i) {
        if (i == 1) {
            return 0.8786796564403575_rt;
        }
        return 0.2928932188134525_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real e (const int i) {
        (void) i;
        return 0.2928932188134525_rt;
    }
};

struct sandu_a_tableau {
    static constexpr int stages = 3;
    static constexpr amrex::Real gamma = 0.7886751345948129_rt; // (3 + sqrt(3)) / 6

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real ct (const int i) {
        return i == 1 ? 0.0_rt : 1.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real a (const int i, const int j) {
        if ((i == 2 || i == 3) && j == 1) {
            return 1.2679491924311228_rt; // 3 - sqrt(3)
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real c (const int i, const int j) {
        if (i == 2 && j == 1) {
            return 0.9282032302755092_rt; // -6 + 4 * sqrt(3)
        }
        if (i == 3 && (j == 1 || j == 2)) {
            return -0.4641016151377544_rt; // 3 - 2 * sqrt(3)
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real b (const int i) {
        if (i == 1) {
            return 1.2679491924311228_rt; // 3 - sqrt(3)
        }
        if (i == 3) {
            return 1.0_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real e (const int i) {
        return i == 3 ? 1.0_rt : 0.0_rt;
    }
};

struct sandu_b_tableau {
    static constexpr int stages = 3;
    static constexpr amrex::Real gamma = 0.7886751345948129_rt; // (3 + sqrt(3)) / 6

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real ct (const int i) {
        return i == 1 ? 0.0_rt : 1.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real a (const int i, const int j) {
        if ((i == 2 || i == 3) && j == 1) {
            return 1.2679491924311228_rt; // 3 - sqrt(3)
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real c (const int i, const int j) {
        if (i == 3 && j == 1) {
            return -0.5358983848622456_rt; // -4 + 2 * sqrt(3)
        }
        if (i == 3 && j == 2) {
            return -0.7320508075688772_rt; // 1 - sqrt(3)
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real b (const int i) {
        if (i == 1) {
            return 1.2679491924311228_rt; // 3 - sqrt(3)
        }
        if (i == 3) {
            return 1.0_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real e (const int i) {
        return i == 3 ? 1.0_rt : 0.0_rt;
    }
};

struct sandu_c_tableau {
    static constexpr int stages = 3;
    static constexpr amrex::Real gamma = 0.78867513459481275_rt; // (3 + sqrt(3)) / 6

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real ct (const int i) {
        if (i == 3) {
            return 0.6666666666666666_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real a (const int i, const int j) {
        if (i == 3 && j == 1) {
            return -0.4335531486700170_rt; // (1272 - 823 * sqrt(3)) / 354
        }
        if (i == 3 && j == 2) {
            return 1.7222282832648295_rt; // 3 * (-51 + 49 * sqrt(3)) / 59
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real c (const int i, const int j) {
        if (i == 2 && j == 1) {
            return -0.3264246859454366_rt; // -1 + 7 * sqrt(3) / 18
        }
        if (i == 3 && j == 1) {
            return 0.4138899169340995_rt; // (25 - 13 * sqrt(3)) / 6
        }
        if (i == 3 && j == 2) {
            return 1.6076951545867360_rt; // 12 - 6 * sqrt(3)
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real b (const int i) {
        if (i == 1) {
            return -7.1285171234651825_rt; // (-12089 + 5037 * sqrt(3)) / 472
        }
        if (i == 2) {
            return 8.4030361763035106_rt; // 9 * (344 - 135 * sqrt(3)) / 118
        }
        return 0.9509618943233421_rt; // 3 * (3 - sqrt(3)) / 4
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real e (const int i) {
        if (i == 1) {
            return -8.396466315896305_rt;
        }
        if (i == 2) {
            return 8.4030361763035106_rt;
        }
        return 0.9509618943233421_rt;
    }
};

struct sandu_d_tableau {
    static constexpr int stages = 4;
    static constexpr amrex::Real gamma = 0.5_rt;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real ct (const int i) {
        return i == 1 ? 0.0_rt : 1.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real a (const int i, const int j) {
        if (i > 1 && j == 1) {
            return 2.0_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real c (const int i, const int j) {
        if (i == 2 && j == 1) {
            return -1.3333333333333333_rt; // -4 / 3
        }
        if (i == 3 && j == 1) {
            return -3.3333333333333333_rt; // -10 / 3
        }
        if (i == 3 && j == 2) {
            return -2.0_rt;
        }
        if (i == 4 && j == 1) {
            return -0.5_rt;
        }
        if (i == 4 && j == 3) {
            return 1.5_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real b (const int i) {
        if (i == 1) {
            return 2.0_rt;
        }
        if (i == 4) {
            return 1.0_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real e (const int i) {
        if (i == 1) {
            return -0.6666666666666666_rt; // 2 - 8 / 3
        }
        if (i == 2 || i == 3) {
            return -1.0_rt;
        }
        return 1.3333333333333333_rt; // 1 - (-1 / 3)
    }
};

using coefficients = ros2s_tableau;

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
    amrex::Array1D<amrex::Real, 1, int_neqs> ysave;

    for (int n = 1; n <= int_neqs; ++n) {
        ysave(n) = ros2s.y(n);
        ros2s.y(n) = y(n);
    }

    rhs(time, state, ros2s, ydot);

    for (int n = 1; n <= int_neqs; ++n) {
        ros2s.y(n) = ysave(n);
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

    ros2s.n_rhs = 0;
    ros2s.n_jac = 0;
    ros2s.n_step = 0;

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

            const amrex::Real err = ros2s::error_norm(ros2s);
            const bool valid_update = ros2s::species_change_valid(ros2s);
            constexpr amrex::Real err_min = 1.e-10_rt;
            const amrex::Real controller_fac = amrex::Clamp(
                std::pow(1.0_rt / amrex::max(err, err_min),
                         1.0_rt / (integrator_rp::h211b_b * integrator_rp::h211b_k)) *
                std::pow(1.0_rt / amrex::max(errold, err_min),
                         1.0_rt / (integrator_rp::h211b_b * integrator_rp::h211b_k)) *
                std::pow(facold, -1.0_rt / integrator_rp::h211b_b),
                integrator_rp::h211b_fac_min, integrator_rp::h211b_fac_max);
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
                hnew *= integrator_rp::h211b_reduction_fac;
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

    return rosenbrock_integrator<BurnT, ros2s::ros2s_tableau>(state, ros2s);
}

#endif
