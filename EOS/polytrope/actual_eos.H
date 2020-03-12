#ifndef _actual_eos_H_
#define _actual_eos_H_

// This is the equation of state for a polytropic fluid:
// P = K rho^gamma
//
// The internal energy is given by a gamma law:
//
// e = (P / rho) * (1 / (gamma - 1))
//
// Unlike the gamma law EOS, e is always a dependent variable
// that is directly determined by the fluid density. This guarantees
// that the fluid always obeys the polytropic relationship.
//
// gamma and K are fixed quantities for the run, and must either be
// supplied by the user or selected from a list of available options.
// Currently, we have fully degenerate ionized gases (both relativistic
// and non-relativistic), where the pressure is supplied by electrons.
//
// Note that here we define the mean number of electrons per ion as:
//
//   1/mu_e = sum_k { X_k Z_k / A_k }
//
// This is assumed to be constant for the degenerate gases.

#include <AMReX.H>
#include <extern_parameters.H>
#include <eos_type.H>
#include <actual_eos_data.H>

const std::string eos_name = "polytrope";

inline
void actual_eos_init ()
{

    // Available pre-defined polytrope options:

    // 1: Non-relativistic, fully degenerate electron gas
    // 2: Relativistic, fully degenerate electron gas 

    if (polytrope_type > 0) {
        mu_e = polytrope_mu_e;
        polytrope = polytrope_type;

        if (polytrope == 1) {
            gamma_const = 5.0_rt / 3.0_rt;
            K_const     = 9.9154e12_rt; // (3 / pi)^(2/3) * h^2 / (20 * m_e * m_p^(5/3))
            K_const     = K_const / std::pow(mu_e, gamma_const);
        }
        else if (polytrope == 2) {
            gamma_const = 4.0_rt / 3.0_rt;
            K_const     = 1.2316e15_rt; // (3 / pi)^(1/3) * h c / (8 * m_p^(4/3))
            K_const     = K_const / std::pow(mu_e, gamma_const);
        }
        else {
            amrex::Error("EOS: Polytrope type currently not defined");
        }
    }
    else if (polytrope_gamma > 0.0_rt && polytrope_K > 0.0_rt) {
        gamma_const = polytrope_gamma;
        K_const     = polytrope_K;
        mu_e        = 2.0_rt; // This will not be used
    }
    else {
        amrex::Error("EOS: Neither polytrope type nor both gamma and K are defined");
    }

    gm1 = gamma_const - 1.0_rt;

    polytrope_index = 1.0_rt / (gamma_const - 1.0_rt);

}



//---------------------------------------------------------------------------
// Public interfaces 
//---------------------------------------------------------------------------

inline
void eos_get_polytrope_parameters (int& polytrope_out, Real& gamma_out, Real& K_out, Real& mu_e_out)
{

    polytrope_out = polytrope;
    gamma_out     = gamma_const;
    K_out         = K_const;
    mu_e_out      = mu_e;

}

inline
void eos_set_polytrope_parameters (int polytrope_in, int gamma_in, int K_in, int mu_e_in)
{

    polytrope   = polytrope_in;
    gamma_const = gamma_in;
    K_const     = K_in;
    mu_e        = mu_e_in;

}



//---------------------------------------------------------------------------
// The main interface
//---------------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE inline
void actual_eos (eos_input_t input, eos_t& state)
{

    Real dens = state.rho;
    Real temp = state.T;
    Real pres = state.p;
    Real enth = state.h;
    Real eint = state.e;
    Real entr = state.s;

    switch (input) {

    //-------------------------------------------------------------------------
    // Now do the calculations. In every case,
    // make sure we have pressure, density, energy, and enthalpy.
    // Relevant equations:
    // h   = e + p / rho = (p / rho) * gamma / (gamma - 1) = e * gamma
    // p   = K * (rho ** gamma) = (gamma - 1) * rho * e
    // rho = (p / K)**(1 / gamma)
    // e   = h - p / rho = (p / rho) / (gamma - 1)         = h / gamma
    //-------------------------------------------------------------------------

    case eos_input_rh:

        // dens, enthalpy, and xmass are inputs

        // Solve for the pressure and energy:

        pres = enth * dens * gm1 / gamma_const;
        eint = enth / gamma_const;

        break;

    case eos_input_rt:

        // dens, temp, and xmass are inputs

        // Solve for the pressure, energy and enthalpy:

        pres = K_const * std::pow(dens, gamma_const);
        enth = pres / dens * gamma_const / gm1;
        eint = enth / gamma_const;

        break;

    case eos_input_tp:

        // temp, pres, and xmass are inputs

        // Solve for the density, energy and enthalpy:

        dens = std::pow(pres / K_const, 1.0_rt / gamma_const);
        enth = pres / dens * gamma_const / gm1;
        eint = enth / gamma_const;

        break;

    case eos_input_rp:

        // dens, pres, and xmass are inputs

        // Solve for the enthalpy and energy:

        enth = (pres / dens) * gamma_const / gm1;
        eint = enth / gamma_const;

        break;

    case eos_input_re:

        // dens, energy, and xmass are inputs

        // Solve for the pressure and enthalpy:

        pres = K_const * std::pow(dens, gamma_const);
        enth = eint * gamma_const;

        break;

    case eos_input_ps:

        // pressure, entropy and xmass are inputs

        // Solve for the density, energy and enthalpy:

        dens = std::pow(pres / K_const, 1.0_rt / gamma_const);
        enth = pres / dens * gamma_const / gm1;
        eint = enth / gamma_const;

        break;

    case eos_input_ph:

        // pressure, enthalpy and xmass are inputs

        // Solve for the density and energy:

        dens = std::pow(pres / K_const, 1.0_rt / gamma_const);
        eint = (pres / dens) * 1.0_rt / gm1;

        break;

    case eos_input_th:

        // temperature, enthalpy and xmass are inputs

        // Solve for the density, energy and pressure:

        eint = enth / gamma_const;
        dens = std::pow(gm1 / gamma_const * enth / K_const, 1.0_rt / gm1);
        pres = gm1 * dens * eint;

        break;

    default:

#if !(defined(AMREX_USE_ACC) || defined(AMREX_USE_CUDA))
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

    // Compute the thermodynamic derivatives and specific heats 
    state.dpdT = 0.0_rt;
    state.dpdr = gamma_const * pres / dens;
    state.dedT = 0.0_rt;
    state.dedr = pres / (dens * dens);
    state.dsdT = 0.0_rt;
    state.dsdr = 0.0_rt;
    state.dhdT = 0.0_rt;
    state.dhdr = state.dedr + gm1 * pres / (dens * dens);

    state.dpde = 0.0_rt;
    state.dpdr_e = gamma_const * pres / dens;
    
    state.cv = state.dedT;
    state.cp = gamma_const * state.cv;

    state.gam1 = gamma_const;

#ifdef EXTRA_THERMO
    // Compute dpdX, dedX, dhdX.

    state.dpdA = - state.p / state.abar;
    state.dpdZ =   state.p / (1.0_rt + state.zbar);

    state.dedA = - state.e / state.abar;
    state.dedZ =   state.e / (1.0_rt + state.zbar);
#endif

    // sound speed
    state.cs = std::sqrt(gamma_const * pres / dens);

}

inline
void actual_eos_finalize ()
{
}

#endif
