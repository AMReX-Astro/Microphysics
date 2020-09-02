#ifndef SDC_TYPE_H
#define SDC_TYPE_H

#include <AMReX_REAL.H>

#include <network.H>
#include <eos_type.H>
#include <extern_parameters.H>

using namespace amrex;

// A generic structure holding data necessary to do a nuclear burn
// in the simplified-SDC formalism.

#if defined(SDC_EVOLVE_ENERGY)

// these indicies represent the order that the conserved state comes
// into the ODE integration from the hydro code.
//
// they also represent the order of the advective sources
//
// integrate rho*X, internal energy, total energy
// carry momentum as an unevolved variable

const int SEDEN = 0;
const int SEINT = 1;
const int SFS   = 2;
const int SRHO  = SFS + nspec;
const int SMX   = SRHO + 1;
const int SMY   = SRHO + 2;
const int SMZ   = SRHO + 3;

const int SVAR  = SMZ+1;
const int SVAR_EVOLVE = SRHO;

#elif defined(SDC_EVOLVE_ENTHALPY)

// integrate rho*X (species masses) and rho*h (enthalpy)
// carry pressure for EOS calls in RHS

const int SFS = 0;
const int SENTH = SFS + nspec;
const int SVAR  = SENTH+1;
const int SVAR_EVOLVE = SVAR;

#endif

struct sdc_t {

  Real y[SVAR];
  Real ydot_a[SVAR];

#if defined(SDC_EVOLVE_ENERGY)
  bool T_from_eden;
#elif defined(SDC_EVOLVE_ENTHALPY)
  // Pressure in case we wish to use it for EOS calls
  Real p0;
  // Density is defined by sum(rho*X) = rho in this method
  Real rho;
#endif

  int i;
  int j;
  int k;

  int n_rhs;
  int n_jac;

  int sdc_iter;
  int num_sdc_iters;

  bool success;
}

#endif
