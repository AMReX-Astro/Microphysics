#ifndef HELM_INTERFACE_H
#define HELM_INTERFACE_H

#include <iostream>
#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_ParallelDescriptor.H>
#include <extern_parameters.H>
#include <eos_type.H>
#include <helmholtz_core.H>




template <typename I, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void core_eos (I input, T& state)
{

    using namespace eos_rp;

    static_assert(std::is_same<I, eos_input_t>::value, "input must be an eos_input_t");

    using namespace helmholtz;

    constexpr int max_newton = 100;

    bool single_iter, converged;
    int var{}, dvar{}, var1, var2;
    Real v_want{}, v1_want{}, v2_want{};

    prepare_for_iterations(input, state, single_iter, v_want, v1_want, v2_want, var, dvar, var1, var2);

    converged = false;

    // Only take a single step if we're coming in with both rho and T;
    // in this call we just want the EOS to fill the other thermodynamic quantities.

    if (input == eos_input_rt) converged = true;

    // Iterate until converged.

    for (int iter = 1; iter <= max_newton; ++iter) {

        // Radiation must come first since it initializes the
        // state instead of adding to it.

        apply_radiation(state);

        apply_ions(state);

        apply_electrons(state);

        if (do_coulomb) {
            apply_coulomb_corrections(state);
        }

        // Calculate enthalpy the usual way, h = e + p / rho.

        if constexpr (has_enthalpy<T>::value) {
            state.h = state.e + state.p / state.rho;
            state.dhdr = state.dedr + state.dpdr / state.rho - state.p / (state.rho * state.rho);
            state.dhdT = state.dedT + state.dpdT / state.rho;
        }

        if (converged) {
            break;
        }
        else if (single_iter) {
            single_iter_update(state, var, dvar, v_want, converged);
        }
        else {
            double_iter_update(state, var1, var2, v1_want, v2_want, converged);
        }

    }

    finalize_state(input, state, v_want, v1_want, v2_want);
}

template <typename I, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void core_eos_interface (const I input, T& state, bool use_raw_inputs = false)
{
  static_assert(std::is_same<I, eos_input_t>::value, "input must be an eos_input_t");

  // Input arguments

  bool has_been_reset = false;
  bool use_composition_routine = true;

  // Local variables

#ifndef AMREX_USE_GPU
  if (!EOSData::initialized) {
    amrex::Error("EOS: not initialized");
  }
#endif

  if (use_raw_inputs) {
    use_composition_routine = false;
  }

  if constexpr (has_xn<T>::value) {
      if (use_composition_routine) {
          // Get abar, zbar, etc.
          composition(state);
      }
  }

  // Force the inputs to be valid.
  reset_inputs(input, state, has_been_reset);

  // Allow the user to override any details of the
  // EOS state. This should generally occur right
  // before the actual_eos call.
  eos_override(state);

  // Call the EOS.

  if (!has_been_reset) {
    core_eos(input, state);
  }
}


#endif
