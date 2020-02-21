// This is a constant gamma equation of state, using an ideal gas.
//
// The gas may either be completely ionized or completely neutral.
//
// The ratio of specific heats (gamma) is allowed to vary.  NOTE: the
// expression for entropy is only valid for an ideal MONATOMIC gas
// (gamma = 5/3).

#include <actual_eos.H>
#include <extern_parameters.H>
#include <eos.H>
#include <fundamental_constants.H>
#include <cmath>

void actual_eos_cxx_init() {

  // constant ratio of specific heats
  if (eos_gamma <= 0.0) {
    amrex::Error("gamma_const cannot be < 0");
  }

}

void actual_eos_cxx(const eos_input_t input, eos_t& state) {

  // Get the mass of a nucleon from Avogadro's number.
  const Real m_nucleon = 1.0 / n_A;
  const Real fac = 1.0 / std::pow(2.0 * M_PI * hbar * hbar, 1.5);

  if (eos_assume_neutral) {
    state.mu = state.abar;

  } else {
    Real sum = 0.0;
    for (int n = 0; n < NumSpec; n++) {
      sum += (zion[n] + 1.0) * state.xn[n] / aion[n];
    }
    state.mu = 1.0 / sum;
  }

  // For all EOS input modes EXCEPT eos_input_rt, first compute dens
  // and temp as needed from the inputs.

  switch (input)
    {

    case eos_input_rt:

      // dens, temp and xmass are inputs

      // We don't need to do anything here
      break;

    case eos_input_rh:

      // dens, enthalpy, and xmass are inputs

      // Solve for the temperature:
      // h = e + p/rho = (p/rho)*[1 + 1/(gamma-1)] = (p/rho)*gamma/(gamma-1)

      state.T = (state.h * state.mu * m_nucleon / k_B)*(eos_gamma - 1.0)/eos_gamma;

      break;

    case eos_input_tp:

      // temp, pres, and xmass are inputs

      // Solve for the density:
      // p = rho k T / (mu m_nucleon)

      state.rho =  state.p * state.mu * m_nucleon / (k_B * state.T);

      break;

    case eos_input_rp:

      // dens, pres, and xmass are inputs

      // Solve for the temperature:
      // p = rho k T / (mu m_nucleon)

      state.T = state.p * state.mu * m_nucleon / (k_B * state.rho);

      break;

    case eos_input_re:

      // dens, energy, and xmass are inputs

      // Solve for the temperature
      // e = k T / [(mu m_nucleon)*(gamma-1)]

      state.T = state.e * state.mu * m_nucleon * (eos_gamma - 1.0) / k_B;

      break;

    case eos_input_ps:

      // pressure, entropy, and xmass are inputs

      // Solve for the temperature
      // Invert Sackur-Tetrode eqn (below) using
      // rho = p mu m_nucleon / (k T)

      state.T = std::pow(state.p, 2.0/5.0) *
        std::pow(2.0*M_PI*hbar*hbar/(state.mu*m_nucleon), 3.0/5.0) *
        exp(2.0*state.mu*m_nucleon*state.s/(5.0*k_B) - 1.0) / k_B;

      // Solve for the density
      // rho = p mu m_nucleon / (k T)

      state.rho = state.p * state.mu * m_nucleon / (k_B * state.T);

      break;

    case eos_input_ph:

      // pressure, enthalpy and xmass are inputs

      // Solve for temperature and density

      state.rho = state.p / state.h * eos_gamma / (eos_gamma - 1.0);
      state.T = state.p * state.mu * m_nucleon / (k_B * state.rho);

      break;

    case eos_input_th:

      // temperature, enthalpy and xmass are inputs

      // This system is underconstrained.

#ifndef AMREX_USE_GPU
      amrex::Error("EOS: eos_input_th is not a valid input for the gamma law EOS.");
#endif

      break;

    default:

#ifndef AMREX_USE_GPU
      amrex::Error("EOS: invalid input.");
#endif

      break;

    }

  // Now we have the density and temperature (and mass fractions /
  // mu), regardless of the inputs.

  Real Tinv = 1.0 / state.T;
  Real rhoinv = 1.0 / state.rho;

  // Compute the pressure simply from the ideal gas law, and the
  // specific internal energy using the gamma-law EOS relation.
  state.p = state.rho * state.T * (k_B / (state.mu*m_nucleon));
  state.e = state.p/(eos_gamma - 1.0) * rhoinv;

  // enthalpy is h = e + p/rho
  state.h = state.e + state.p * rhoinv;

  // entropy (per gram) of an ideal monoatomic gas (the Sackur-Tetrode equation)
  // NOTE: this expression is only valid for gamma = 5/3.
  state.s = (k_B/(state.mu*m_nucleon))*
    (2.5 + log((std::pow(state.mu*m_nucleon, 2.5) * rhoinv ) *
               std::pow(k_B * state.T, 1.5) * fac ) );

  // Compute the thermodynamic derivatives and specific heats
  state.dpdT = state.p * Tinv;
  state.dpdr = state.p * rhoinv;
  state.dedT = state.e * Tinv;
  state.dedr = 0.0;
  state.dsdT = 1.5 * (k_B / (state.mu * m_nucleon)) * Tinv;
  state.dsdr = - (k_B / (state.mu * m_nucleon)) * rhoinv;
  state.dhdT = state.dedT + state.dpdT * rhoinv;
  state.dhdr = 0.0;

  state.xne = 0.0;
  state.xnp = 0.0;
  state.eta = 0.0;
  state.pele = 0.0;
  state.ppos = 0.0;

  state.cv = state.dedT;
  state.cp = eos_gamma * state.cv;

  state.gam1 = eos_gamma;

  state.dpdr_e = state.dpdr - state.dpdT * state.dedr * (1.0/state.dedT);
  state.dpde = state.dpdT * (1.0/state.dedT);

  // sound speed
  state.cs = std::sqrt(eos_gamma * state.p * rhoinv);

#ifdef EXTRA_THERMO
  state.dpdA = - state.p * (1.0/state.abar);
  state.dedA = - state.e * (1.0/state.abar);

  if (eos_assume_neutral) {
      state.dpdZ = 0.0;
      state.dedZ = 0.0;
  } else {
    state.dpdZ = state.p * (1.0/(1.0 + state.zbar));
    state.dedZ = state.e * (1.0/(1.0 + state.zbar));
  } 
#endif
}

void actual_eos_cxx_finalize() {

}
