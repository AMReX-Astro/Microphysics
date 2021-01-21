#ifndef _actual_eos_H_
#define _actual_eos_H_

// This is the equation of state for zero-temperature white dwarf 
// matter composed of degenerate electrons:
// P = A * (x * (2x**2 - 3)(x**2 + 1)**1/2 + 3 sinh**-1(x))
// 
// where rho = B x**3 and the constants are given by:
//
// A = pi m_e**4 c**5 / (3 h**3) = 6.0 x 10^22 dyne cm**-2
// B = 8 pi m_e**3 c**3 mu_e m_p  / (3 h**3) = 9.8 x 10^5 mu_e g cm**-3
//
// The equation of state comes from Chandrasekhar (1935), and the enthalpy
// is calculated by Hachisu (1986):
//
// h = (8A / B) (1 + x**2)**(1/2)
//
// The internal energy is calculated using the standard relation:
// 
// h = e + P / rho

#include <AMReX.H>
#include <fundamental_constants.H>
#include <eos_type.H>

using namespace amrex;

const std::string eos_name = "ztwd";
  
const Real A = M_PI * std::pow(C::m_e, 4) * std::pow(C::c_light, 5) / (3.0_rt * std::pow(C::hplanck, 3));
const Real B2 = 8.0_rt * M_PI * std::pow(C::m_e, 3) * std::pow(C::c_light, 3) * C::m_p  / (3.0_rt * std::pow(C::hplanck, 3));
const Real iter_tol = 1.e-10_rt;
const int  max_iter = 1000;

inline
void actual_eos_init ()
{
}


template <typename I>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool is_input_valid (I input)
{
  static_assert(std::is_same<I, eos_input_t>::value, "input must be an eos_input_t");

  bool valid = true;

  return valid;
}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real pressure (Real x)
{
    return A * (x * (2.0_rt * x * x - 3.0_rt) * std::sqrt(x * x + 1.0_rt) + 3.0_rt * std::asinh(x));
}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real enthalpy (Real x, Real B)
{
    return (8.0_rt * A / B) * std::sqrt(1.0_rt + x * x);
}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real dpdx (Real x)
{
    return A * ((2.0_rt * x * x - 3.0_rt) * std::sqrt(x * x + 1.0_rt) +
                x * (4.0_rt * x) * std::sqrt(x * x + 1.0_rt) +
                x * x * (2.0_rt * x * x - 3.0_rt) / std::sqrt(x * x + 1.0_rt) +
                3.0_rt / std::sqrt(x * x + 1.0_rt));
}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real dhdx (Real x, Real B)
{
    return enthalpy(x, B) * (x / (x * x + 1.0_rt));
}




AMREX_GPU_HOST_DEVICE AMREX_INLINE
void pres_iter (Real pres, Real& dens, Real B)
{

    // Starting guess for the iteration.

    Real x = 1.0_rt;

    // We are solving the equation:
    // f(x) = p_want - p(x) = 0.
    // Then we can use Newton's method, with dfdx = -dpdx.
    // We iterate until the density change is close enough to zero.

    int iter;

    for (iter = 1; iter <= max_iter; ++iter)
    {
        Real dx = (pres - pressure(x)) / dpdx(x);

        x = x + dx;

        if (std::abs(dx) / x < iter_tol) {
            break;
        }
    }

#ifndef AMREX_USE_CUDA
    if (iter == max_iter) {
        amrex::Error("EOS: pres_iter failed to converge.");
    }
#endif

    dens = B * x * x * x;

}



template <typename I, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_eos (I input, T& state)
{
    static_assert(std::is_same<I, eos_input_t>::value, "input must be an eos_input_t");

    Real dens = state.rho;
    Real temp = state.T;

    Real pres = 1.0_rt;
    if constexpr (has_pressure<T>::value) {
        pres = state.p;
    }

    Real enth = 1.0_rt;
    if constexpr (has_enthalpy<T>::value) {
        enth = state.h;
    }

    Real eint = 1.0_rt;
    if constexpr (has_energy<T>::value) {
        eint = state.e;
    }

    Real entr = 1.0_rt;
    if constexpr (has_entropy<T>::value) {
        entr = state.s;
    }

    Real B = B2 * state.mu_e;

    Real x, dxdr;

    switch (input) {

    //-------------------------------------------------------------------------
    // Now do the calculations. In every case,
    // make sure we have pressure, density, energy, and enthalpy.
    // Relevant equations:
    // rho = B x**3
    // p   = A ( x * (2x**2 - 3)(x**2 + 1)**1/2 + 3 sinh**-1(x) )
    // h   = (8A / B) * (1 + x**2)**1/2
    // e   = h - p / rho
    //-------------------------------------------------------------------------

    case eos_input_rh:

        // dens, enthalpy, and xmass are inputs

        // Solve for the pressure and energy:

        x = std::pow(dens / B, 1.0_rt / 3.0_rt);
        pres = pressure(x);
        eint = enth - pres / dens;

        break;

    case eos_input_rt:

        // dens, temp, and xmass are inputs

        // Solve for the pressure, energy and enthalpy:

        x = std::pow(dens / B, 1.0_rt / 3.0_rt);
        pres = pressure(x);
        enth = enthalpy(x, B);
        eint = enth - pres / dens;

        break;

    case eos_input_tp:

        // temp, pres, and xmass are inputs

        // Solve for the density, energy and enthalpy:

        pres_iter(pres, dens, B);

        x = std::pow(dens / B, 1.0_rt / 3.0_rt);
        enth = enthalpy(x, B);
        eint = enth - pres / dens;

        break;

    case eos_input_rp:

        // dens, pres, and xmass are inputs

        // Solve for the enthalpy and energy:

        x = std::pow(dens / B, 1.0_rt / 3.0_rt);
        enth = enthalpy(x, B);
        eint = enth - pres / dens;

        break;

    case eos_input_re:

        // dens, energy, and xmass are inputs

        // Solve for the pressure and enthalpy:

        x = std::pow(dens / B, 1.0_rt / 3.0_rt);
        pres = pressure(x);
        enth = enthalpy(x, B);

        break;

    case eos_input_ps:

        // pressure, entropy and xmass are inputs

        // Solve for the density, energy and enthalpy:

        pres_iter(pres, dens, B);

        x = std::pow(dens / B, 1.0_rt / 3.0_rt);
        enth = enthalpy(x, B);
        eint = enth - pres / dens;

        break;

    case eos_input_ph:

        // pressure, enthalpy and xmass are inputs

        // Solve for the density and energy:

        x = std::sqrt((B * enth) / (8.0_rt * A) * (B * enth) / (8.0_rt * A) - 1.0_rt);
        dens = B * x * x * x;
        eint = enth - pres / dens;

        break;

    case eos_input_th:

        // temperature, enthalpy and xmass are inputs

        // Solve for the density, energy and pressure:

        x = std::sqrt((B * enth) / (8.0_rt * A) * (B * enth) / (8.0_rt * A) - 1.0_rt);
        dens = B * x * x * x;
        pres = pressure(x);
        eint = enth - pres / dens;

        break;

    default:

#ifndef AMREX_USE_GPU
        amrex::Error("EOS: invalid input.");
#endif

        break;

    }

    //-------------------------------------------------------------------------
    // Now we have all relevant quantities, regardless of the inputs.
    //-------------------------------------------------------------------------

    state.T   = temp;
    state.rho = dens;
    if constexpr (has_enthalpy<T>::value) {
        state.h   = enth;
    }
    if constexpr (has_entropy<T>::value) {
        state.s   = entr;
    }
    if constexpr (has_energy<T>::value) {
        state.e   = eint;
    }
    if constexpr (has_pressure<T>::value) {
        state.p   = pres;
    }

    // All temperature derivatives are zero since the gas is temperature-independent.

    if constexpr (has_pressure<T>::value) {
        state.dpdT = 0.0_rt;
    }
    if constexpr (has_enthalpy<T>::value) {
        state.dhdT = 0.0_rt;
    }
    if constexpr (has_energy<T>::value) {
        state.dedT = 0.0_rt;
    }
    if constexpr (has_entropy<T>::value) {
        state.dsdT = 0.0_rt;
    }

    // Density derivatives are computed using the chain rule, e.g. dpdr = dpdx * dxdr.

    x = std::pow(dens / B, 1.0_rt / 3.0_rt);
    dxdr = (1.0_rt / 3.0_rt) * x / dens;

    Real dpdr = dxdr * dpdx(x);
    Real dhdr = dxdr * dhdx(x, B);

    if constexpr (has_pressure<T>::value) {
        state.dpdr = dpdr;
    }
    if constexpr (has_enthalpy<T>::value) {
        state.dhdr = dhdr;
    }
    if constexpr (has_energy<T>::value) {
        state.dedr = dhdr - dpdr / state.rho + pres / (state.rho * state.rho);
    }
    if constexpr (has_entropy<T>::value) {
        state.dsdr = 0.0_rt;
    }

    // Heat capacities are zero: the gas properties don't change when the temperature changes.

    if constexpr (has_energy<T>::value) {
        state.cv = 0.0_rt;
    }
    if constexpr (has_pressure<T>::value && has_energy<T>::value) {
        state.cp = 0.0_rt;
    }

    // Adiabatic gamma_1 == d(log p) / d(log rho) |_s.

    if constexpr (has_pressure<T>::value && has_energy<T>::value) {
        state.gam1 = state.dpdr * (state.rho / state.p);
    }

#ifdef EXTRA_THERMO
    // Derivatives with respect to A and Z.

    if constexpr (has_pressure<T>::value) {
        state.dpdA = -state.p / state.abar;
        state.dpdZ =  state.p / (1.0_rt + state.zbar);
    }

    if constexpr (has_energy<T>::value) {
        state.dedA = -state.e / state.abar;
        state.dedZ =  state.e / (1.0_rt + state.zbar);
    }
#endif

    // Sound speed.

    if constexpr (has_pressure<T>::value && has_energy<T>::value) {
        state.cs = std::sqrt(state.dpdr);
    }
}



inline
void actual_eos_finalize ()
{
}

#endif
