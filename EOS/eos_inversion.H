#ifndef EOS_INVERSION_H
#define EOS_INVERSION_H

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <eos_type.H>

namespace eos_inversion {

constexpr int max_newton = 100;

// The input mode determines the targets and independent variables once,
// before entering the Newton loop.
template <eos_input_t Input>
struct inversion_traits {
    static constexpr bool double_iter = Input == eos_input_ps || Input == eos_input_ph;
    static constexpr EOS::eos_iter var =
        Input == eos_input_rh || Input == eos_input_th ? EOS::ienth :
        Input == eos_input_re ? EOS::iener : EOS::ipres;
    static constexpr EOS::eos_iter var2 = Input == eos_input_ps ? EOS::ientr : EOS::ienth;
    static constexpr EOS::eos_iter dvar =
        Input == eos_input_tp || Input == eos_input_th ? EOS::idens : EOS::itemp;
};

template <EOS::eos_iter Var, typename T>
inline constexpr bool has_quantity =
    (Var == EOS::ipres && has_pressure<T>) ||
    (Var == EOS::iener && has_energy<T>) ||
    (Var == EOS::ientr && has_entropy<T>) ||
    (Var == EOS::ienth && has_enthalpy<T>);

struct eos_quantity {
    amrex::Real value;
    amrex::Real dT;
    amrex::Real drho;
};

// Share the property selection between the one- and two-variable updates.
template <EOS::eos_iter Var, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
eos_quantity quantity (const T& state)
{
    static_assert(has_quantity<Var, T>, "EOS state is missing the inversion target");

    if constexpr (Var == EOS::ipres) {
        return {state.p, state.dpdT, state.dpdr};
    } else if constexpr (Var == EOS::iener) {
        return {state.e, state.dedT, state.dedr};
    } else if constexpr (Var == EOS::ientr) {
        return {state.s, state.dsdT, state.dsdr};
    } else if constexpr (Var == EOS::ienth) {
        return {state.h, state.dhdT, state.dhdr};
    }
}


template <EOS::eos_iter Var, EOS::eos_iter DVar, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void single_iter_update (T& state, amrex::Real v_want, bool& converged)
{
    static_assert(DVar == EOS::itemp || DVar == EOS::idens);

    const auto q = quantity<Var>(state);
    amrex::Real x, dvdx, xtol, smallx;

    if constexpr (DVar == EOS::itemp) {
        x = state.T;
        dvdx = q.dT;
        smallx = EOSData::mintemp;
        xtol = eos_rp::eos_ttol;
    } else {
        x = state.rho;
        dvdx = q.drho;
        smallx = EOSData::mindens;
        xtol = eos_rp::eos_dtol;
    }

    // Now do the calculation for the next guess for T/rho
    amrex::Real xnew = x - (q.value - v_want) / dvdx;

    // Don't let the temperature/density change by more than a factor of two
    xnew = amrex::Clamp(xnew, 0.5_rt * x, 2.0_rt * x);

    // Don't let us freeze/evacuate
    xnew = amrex::max(smallx, xnew);

    if constexpr (DVar == EOS::itemp) {
        state.T = xnew;
    } else {
        state.rho = xnew;
    }

    amrex::Real error = std::abs(xnew - x);
    if (error < xtol * x) converged = true;
}


template <EOS::eos_iter Var1, EOS::eos_iter Var2, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void double_iter_update (T& state, amrex::Real v1_want, amrex::Real v2_want,
                         bool& converged)
{
    const auto [v1, dv1dt, dv1dr] = quantity<Var1>(state);
    const auto [v2, dv2dt, dv2dr] = quantity<Var2>(state);

    amrex::Real told = state.T;
    amrex::Real rold = state.rho;

    // Two functions, f and g, to iterate over
    amrex::Real v1i = v1_want - v1;
    amrex::Real v2i = v2_want - v2;

    //
    // 0 = f + dfdr * delr + dfdt * delt
    // 0 = g + dgdr * delr + dgdt * delt
    //

    // note that dfi/dT = - df/dT
    amrex::Real delr = (-v1i * dv2dt + v2i * dv1dt) /
        (dv2dr * dv1dt - dv2dt * dv1dr);

    amrex::Real rnew = rold + delr;

    amrex::Real tnew = told + (v1i - dv1dr * delr) / dv1dt;

    // Don't let the temperature or density change by more
    // than a factor of two
    tnew = amrex::Clamp(tnew, 0.5e0_rt * told, 2.0e0_rt * told);
    rnew = amrex::Clamp(rnew, 0.5e0_rt * rold, 2.0e0_rt * rold);

    // Don't let us freeze or evacuate
    tnew = amrex::max(EOSData::mintemp, tnew);
    rnew = amrex::max(EOSData::mindens, rnew);

    // Store the new temperature and density
    state.rho = rnew;
    state.T = tnew;

    // Compute the errors
    amrex::Real error1 = std::abs(rnew - rold);
    amrex::Real error2 = std::abs(tnew - told);

    if (error1 < eos_rp::eos_dtol * rold &&
        error2 < eos_rp::eos_ttol * told) converged = true;
}


template <eos_input_t Input, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool invert_eos_impl (T& state)
{
    using traits = inversion_traits<Input>;

    // Runtime dispatch instantiates every mode. Discard unsupported bodies so
    // reduced states never instantiate accessors for fields they do not have.
    constexpr bool supported = Input == eos_input_rt ||
        (has_quantity<traits::var, T> &&
         (!traits::double_iter || has_quantity<traits::var2, T>));

    if constexpr (!supported) {
        return false;
    } else {
        amrex::Real v_want{}, v1_want{}, v2_want{};
        if constexpr (Input != eos_input_rt) {
            if constexpr (traits::double_iter) {
                v1_want = quantity<traits::var>(state).value;
                v2_want = quantity<traits::var2>(state).value;
            } else {
                v_want = quantity<traits::var>(state).value;
            }
        }

        // With rho and T supplied, evaluate the EOS only once.
        bool converged = Input == eos_input_rt;

        for (int iter = 1; iter <= max_newton; ++iter) {
            apply_eos(state);

            // Calculate enthalpy the usual way, h = e + p / rho.
            if constexpr (has_enthalpy<T>) {
                state.h = state.e + state.p / state.rho;
                state.dhdr = state.dedr + state.dpdr / state.rho - state.p / (state.rho * state.rho);
                state.dhdT = state.dedT + state.dpdT / state.rho;
            }

            // Evaluate the EOS at the converged guess before finalizing.
            if (converged) {
                break;
            }
            if constexpr (Input != eos_input_rt) {
                if constexpr (traits::double_iter) {
                    double_iter_update<traits::var, traits::var2>(state, v1_want, v2_want, converged);
                } else {
                    single_iter_update<traits::var, traits::dvar>(state, v_want, converged);
                }
            }
        }

        finalize_state(Input, state, v_want, v1_want, v2_want);
        return converged;
    }
}


} // namespace eos_inversion

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool invert_eos (eos_input_t input, T& state)
{
    using namespace eos_inversion;

    switch (input) {
    case eos_input_rt:
        return invert_eos_impl<eos_input_rt>(state);
    case eos_input_rh:
        return invert_eos_impl<eos_input_rh>(state);
    case eos_input_tp:
        return invert_eos_impl<eos_input_tp>(state);
    case eos_input_rp:
        return invert_eos_impl<eos_input_rp>(state);
    case eos_input_re:
        return invert_eos_impl<eos_input_re>(state);
    case eos_input_ps:
        return invert_eos_impl<eos_input_ps>(state);
    case eos_input_ph:
        return invert_eos_impl<eos_input_ph>(state);
    case eos_input_th:
        return invert_eos_impl<eos_input_th>(state);
    }

#ifndef AMREX_USE_GPU
    amrex::Error("Unknown EOS input");
#endif
    return false;
}

#endif
