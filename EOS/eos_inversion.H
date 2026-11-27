#ifndef EOS_INVERSION
#define EOS_INVERSION

#include <AMReX.H>
#include <AMReX_REAL.H>

constexpr int max_newton = 100;

template <typename I, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void prepare_for_iterations (I input, const T& state,
                             bool& single_iter, amrex::Real& v_want, amrex::Real& v1_want, amrex::Real& v2_want,
                             int& var, int& dvar, int& var1, int& var2)
{
    using namespace helmholtz;
    using namespace EOS;

    single_iter = true;

    if (input == eos_input_rt) {

        // Nothing to do here.

    }  else if (input == eos_input_rh) {

        if constexpr (has_enthalpy<T>) {
            v_want = state.h;
            var  = ienth;
            dvar = itemp;
        }

    } else if (input == eos_input_tp) {

        if constexpr (has_pressure<T>) {
            v_want = state.p;
            var  = ipres;
            dvar = idens;
        }

    } else if (input == eos_input_rp) {

        if constexpr (has_pressure<T>) {
            v_want = state.p;
            var  = ipres;
            dvar = itemp;
        }

    } else if (input == eos_input_re) {

        if constexpr (has_energy<T>) {
            v_want = state.e;
            var  = iener;
            dvar = itemp;
        }

    } else if (input == eos_input_ps) {

        if constexpr (has_pressure<T> && has_entropy<T>) {
            single_iter = false;
            v1_want = state.p;
            v2_want = state.s;
            var1 = ipres;
            var2 = ientr;
        }

    } else if (input == eos_input_ph) {

        if constexpr (has_pressure<T> && has_enthalpy<T>) {
            single_iter = false;
            v1_want = state.p;
            v2_want = state.h;
            var1 = ipres;
            var2 = ienth;
        }

    } else if (input == eos_input_th) {

        if constexpr (has_enthalpy<T>) {
            v_want = state.h;
            var  = ienth;
            dvar = idens;
        }

#ifndef AMREX_USE_GPU
    } else {
        amrex::Error("Unknown EOS input");
#endif
    }

}



template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void single_iter_update (T& state, int var, int dvar,
                         amrex::Real v_want, bool& converged)
{

    using namespace EOS;

    amrex::Real x, v = 0.0_rt, dvdx = 0.0_rt, xtol, smallx;

    if (dvar == itemp) {

        x = state.T;
        smallx = EOSData::mintemp;
        xtol = eos_rp::eos_ttol;

        if (var == ipres) {
            if constexpr (has_pressure<T>) {
                v    = state.p;
                dvdx = state.dpdT;
            }
        }
        else if (var == iener) {
            if constexpr (has_energy<T>) {
                v    = state.e;
                dvdx = state.dedT;
            }
        }
        else if (var == ientr) {
            if constexpr (has_entropy<T>) {
                v    = state.s;
                dvdx = state.dsdT;
            }
        }
        else if (var == ienth) {
            if constexpr (has_enthalpy<T>) {
                v    = state.h;
                dvdx = state.dhdT;
            }
        }

    } else {
        // dvar == density

        x = state.rho;
        smallx = EOSData::mindens;
        xtol = eos_rp::eos_dtol;

        if (var == ipres) {
            if constexpr (has_pressure<T>) {
                v    = state.p;
                dvdx = state.dpdr;
            }
        }
        else if (var == iener) {
            if constexpr (has_energy<T>) {
                v    = state.e;
                dvdx = state.dedr;
            }
        }
        else if (var == ientr) {
            if constexpr (has_entropy<T>) {
                v    = state.s;
                dvdx = state.dsdr;
            }
        }
        else if (var == ienth) {
            if constexpr (has_enthalpy<T>) {
                v    = state.h;
                dvdx = state.dhdr;
            }
        }

    }

    // Now do the calculation for the next guess for T/rho
    amrex::Real xnew = x - (v - v_want) / dvdx;

    // Don't let the temperature/density change by more than a factor of two
    xnew = amrex::Clamp(xnew, 0.5_rt * x, 2.0_rt * x);

    // Don't let us freeze/evacuate
    xnew = amrex::max(smallx, xnew);

    // Store the new temperature/density

    if (dvar == itemp) {
        state.T = xnew;
    } else {
        state.rho = xnew;
    }

    // Compute the error from the last iteration

    amrex::Real error = std::abs(xnew - x);

    if (error < xtol * x) converged = true;
}



template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void double_iter_update (T& state, int var1, int var2,
                         amrex::Real v1_want, amrex::Real v2_want, bool& converged)
{
    using namespace EOS;

    amrex::Real v1 = 0.0_rt, dv1dt = 0.0_rt, dv1dr = 0.0_rt, v2 = 0.0_rt, dv2dt = 0.0_rt, dv2dr = 0.0_rt;

    // Figure out which variables we're using

    amrex::Real told = state.T;
    amrex::Real rold = state.rho;

    if (var1 == ipres) {
        if constexpr (has_pressure<T>) {
            v1    = state.p;
            dv1dt = state.dpdT;
            dv1dr = state.dpdr;
        }
    }
    else if (var1 == iener) {
        if constexpr (has_energy<T>) {
            v1    = state.e;
            dv1dt = state.dedT;
            dv1dr = state.dedr;
        }
    }
    else if (var1 == ientr) {
        if constexpr (has_entropy<T>) {
            v1    = state.s;
            dv1dt = state.dsdT;
            dv1dr = state.dsdr;
        }
    }
    else if (var1 == ienth) {
        if constexpr (has_enthalpy<T>) {
            v1    = state.h;
            dv1dt = state.dhdT;
            dv1dr = state.dhdr;
        }
    }

    if (var2 == ipres) {
        if constexpr (has_pressure<T>) {
            v2    = state.p;
            dv2dt = state.dpdT;
            dv2dr = state.dpdr;
        }
    }
    else if (var2 == iener) {
        if constexpr (has_energy<T>) {
            v2    = state.e;
            dv2dt = state.dedT;
            dv2dr = state.dedr;
        }
    }
    else if (var2 == ientr) {
        if constexpr (has_entropy<T>) {
            v2    = state.s;
            dv2dt = state.dsdT;
            dv2dr = state.dsdr;
        }
    }
    else if (var2 == ienth) {
        if constexpr (has_enthalpy<T>) {
            v2    = state.h;
            dv2dt = state.dhdT;
            dv2dr = state.dhdr;
        }
    }

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


template <typename I, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool invert_eos (I input, T& state)
{

    bool single_iter, converged;
    int var{}, dvar{}, var1, var2;
    amrex::Real v_want{}, v1_want{}, v2_want{};

    prepare_for_iterations(input, state, single_iter,
                           v_want, v1_want, v2_want,
                           var, dvar, var1, var2);

    converged = false;

    // Only take a single step if we're coming in with both rho and T;
    // in this call we just want the EOS to fill the other
    // thermodynamic quantities.

    if (input == eos_input_rt) converged = true;

    // Iterate until converged.

    for (int iter = 1; iter <= max_newton; ++iter) {

        // apply the EOS implementation to the state

        apply_eos(state);

        // Calculate enthalpy the usual way, h = e + p / rho.

        if constexpr (has_enthalpy<T>) {
            state.h = state.e + state.p / state.rho;
            state.dhdr = state.dedr + state.dpdr / state.rho - state.p / (state.rho * state.rho);
            state.dhdT = state.dedT + state.dpdT / state.rho;
        }

        if (converged) {
            break;
        } else if (single_iter) {
            single_iter_update(state, var, dvar, v_want, converged);
        } else {
            double_iter_update(state, var1, var2, v1_want, v2_want, converged);
        }

    }

    // do any finalization for the EOS implmentation.

    finalize_state(input, state, v_want, v1_want, v2_want);

    return converged;
}

#endif
