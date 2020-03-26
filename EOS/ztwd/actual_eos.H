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
  
const Real A = M_PI * std::pow(m_e, 4) * std::pow(c_light, 5) / (3.0_rt * std::pow(hplanck, 3));
const Real B2 = 8.0_rt * M_PI * std::pow(m_e, 3) * std::pow(c_light, 3) * m_p  / (3.0_rt * std::pow(hplanck, 3));
const Real iter_tol = 1.e-10_rt;
const int  max_iter = 1000;

inline
void actual_eos_init ()
{
}


AMREX_GPU_HOST_DEVICE inline
bool is_input_valid(eos_input_t input) {

  bool valid = true;

  return valid;
}



AMREX_GPU_HOST_DEVICE inline
Real pressure (Real x)
{
    return A * (x * (2.0_rt * x * x - 3.0_rt) * std::sqrt(x * x + 1.0_rt) + 3.0_rt * std::asinh(x));
}



AMREX_GPU_HOST_DEVICE inline
Real enthalpy (Real x, Real B)
{
    return (8.0_rt * A / B) * std::sqrt(1.0_rt + x * x);
}



AMREX_GPU_HOST_DEVICE inline
Real dpdx (Real x)
{
    return A * ((2.0_rt * x * x - 3.0_rt) * std::sqrt(x * x + 1.0_rt) +
                x * (4.0_rt * x) * std::sqrt(x * x + 1.0_rt) +
                x * x * (2.0_rt * x * x - 3.0_rt) / std::sqrt(x * x + 1.0_rt) +
                3.0_rt / std::sqrt(x * x + 1.0_rt));
}



AMREX_GPU_HOST_DEVICE inline
Real dhdx (Real x, Real B)
{
    return enthalpy(x, B) * (x / (x * x + 1.0_rt));
}




AMREX_GPU_HOST_DEVICE inline
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



AMREX_GPU_HOST_DEVICE inline
void actual_eos (eos_input_t input, eos_t& state)
{

    Real dens = state.rho;
    Real temp = state.T;
    Real pres = state.p;
    Real enth = state.h;
    Real eint = state.e;
    Real entr = state.s;

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

#ifndef AMREX_USE_CUDA
        amrex::Error("EOS: invalid input.");
#endif

        break;

    }

    //-------------------------------------------------------------------------
    // Now we have all relevant quantities, regardless of the inputs.
    //-------------------------------------------------------------------------

    state.T   = temp;
    state.rho = dens;
    state.h   = enth;
    state.s   = entr;
    state.e   = eint;
    state.p   = pres;

    // All temperature derivatives are zero since the gas is temperature-independent.

    state.dpdT = 0.0_rt;
    state.dhdT = 0.0_rt;
    state.dedT = 0.0_rt;
    state.dsdT = 0.0_rt;

    // Density derivatives are computed using the chain rule, e.g. dpdr = dpdx * dxdr.

    x = std::pow(dens / B, 1.0_rt / 3.0_rt);
    dxdr = (1.0_rt / 3.0_rt) * x / dens;

    state.dpdr = dxdr * dpdx(x);
    state.dhdr = dxdr * dhdx(x, B);
    state.dedr = state.dhdr - state.dpdr / state.rho + state.p / (state.rho * state.rho);
    state.dsdr = 0.0_rt;

    // Heat capacities are zero: the gas properties don't change when the temperature changes.

    state.cv = 0.0_rt;
    state.cp = 0.0_rt;

    // Adiabatic gamma_1 == d(log p) / d(log rho) |_s.

    state.gam1 = state.dpdr * (state.rho / state.p);

#ifdef EXTRA_THERMO
    // Derivatives with respect to A and Z.

    state.dpdA = -state.p / state.abar;
    state.dpdZ =  state.p / (1.0_rt + state.zbar);

    state.dedA = -state.e / state.abar;
    state.dedZ =  state.e / (1.0_rt + state.zbar);
#endif

    // Sound speed.

    state.cs = std::sqrt(state.dpdr);

}



inline
void actual_eos_finalize ()
{
}

#endif
