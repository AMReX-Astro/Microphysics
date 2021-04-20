#ifndef _actual_eos_H_
#define _actual_eos_H_

#include <string>
#include <extern_parameters.H>
#include <fundamental_constants.H>
#include <eos_data.H>
#include <cmath>

const std::string eos_name = "tillotson";

inline
void actual_eos_init() {

}

template <typename I>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool is_input_valid (I input)
{
  static_assert(std::is_same<I, eos_input_t>::value, "input must be an eos_input_t");

  bool valid = false;

  if (input == eos_input_rt || input == eos_input_re) {
      // This EOS doesn't have a true dependence on T, so
      // eos_input_re is really the only mode that is consistent
      // with the Tillotson formulation; but since eos_input_rt
      // is commonly needed, we need to support that too.

      valid = true;
  }

  return valid;
}

template <typename I, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_eos (I input, T& state)
{
    static_assert(std::is_same<I, eos_input_t>::value, "input must be an eos_input_t");

    // First make T and e consistent with each other.
    // We use an arbitrary expression for T; although this shouldn't be used for
    // other physics (e.g. chemistry) it is fine for making sure that if we come
    // in with eos_input_rt that we are consistent with eos_input_re.

    if (input == eos_input_re) {
        state.T = amrex::max(state.e / eos_c_v, EOSData::mintemp);
    }
    else if (input == eos_input_rt) {
        state.e = amrex::max(eos_c_v * state.T, EOSData::mine);
    }

    Real eta = state.rho / eos_rho_0;
    Real mu = eta - 1.0_rt;
    Real omega = state.e / (eos_e_0 * eta * eta) + 1.0_rt;
    Real z = (1.0_rt / eta - 1.0_rt);

    Real g_c = eos_la + eos_lb / omega;
    Real g_e = eos_la + eos_lb / omega * std::exp(-eos_beta * z * z);

    Real P_c = g_c * state.e * state.rho + eos_A * mu + eos_B * mu * mu;
    Real P_e = g_e * state.e * state.rho + eos_A * mu * std::exp(-eos_alpha * z - eos_beta * z * z);

    // Floor the pressure since the above expressions can be negative

    P_c = amrex::max(P_c, EOSData::minp);
    P_e = amrex::max(P_e, EOSData::minp);

    Real cs_squared_c = (g_c + 2.0_rt * eos_lb * (omega - 1.0_rt) / (omega * omega)) * state.e +
                        (eos_A + 2.0_rt * eos_B * mu) / eos_rho_0;
    Real cs_squared_e = (g_e + (2.0_rt * eos_lb / omega) * std::exp(-eos_beta * z * z) *
                         ((omega - 1.0_rt) / omega + eos_beta * z / eta)) * state.e +
                        (eos_A / eos_rho_0) * std::exp(-eos_alpha * z - eos_beta * z * z) *
                        (1.0_rt + mu / (eta * eta) * (eos_alpha + 2.0_rt * eos_beta * z));

    // Floor the sound speed as well (and make sure we get the floor dimensionally right)

    cs_squared_c = amrex::max(cs_squared_c, EOSData::minp / state.rho);
    cs_squared_e = amrex::max(cs_squared_e, EOSData::minp / state.rho);

    Real dPde_c = eos_la * state.rho + eos_lb * state.rho / (omega * omega);
    Real dPde_e = eos_la * state.rho + eos_lb * state.rho / (omega * omega) * std::exp(-eos_beta * z * z);

    // At this point we should have valid (rho, e); compute P = P(rho, e).

    if (state.rho >= eos_rho_0 || (state.rho < eos_rho_0 && state.e < eos_e_s)) {
        // Region I, II (compressed state)

        state.p = P_c;
        state.dpdr_e = cs_squared_c;
        state.dpde = dPde_c;
    }
    else if (state.rho < eos_rho_0 && state.e > eos_e_s_prime) {
        // Region IV (expanded state)

        state.p = P_e;
        state.dpdr_e = cs_squared_e;
        state.dpde = dPde_e;
    }
    else {
        // Region III (interpolate between compressed and expanded state)

        Real denom = (eos_e_s_prime - eos_e_s);
        Real xi_c = (state.e - eos_e_s) / denom;
        Real xi_e = (eos_e_s_prime - state.e) / denom;

        state.p = P_c * xi_c + P_e * xi_e;
        state.dpdr_e = cs_squared_c * xi_c + cs_squared_e * xi_e;
        // Note that the interpolation function depends on e, so this term
        // needs to have the derivatives of xi to be accurate.
        state.dpde = dPde_c * xi_c + P_c / denom + dPde_e * xi_e - P_e / denom;
    }

    // Sound speed

    state.cs = std::sqrt(state.dpdr_e);

    // Compute gamma from sound speed

    state.gam1 = state.cs * state.rho / state.p;

    // Fill in heat capacity (c_v from provided data, c_p from c_v)

    state.cv = eos_c_v;
    state.cp = state.gam1 * eos_c_v;

    // For this EOS dpdr == dpdr_e

    state.dpdr = state.dpdr_e;
}

inline
void actual_eos_finalize() {

}

#endif
