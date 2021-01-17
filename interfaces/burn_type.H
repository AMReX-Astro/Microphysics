#ifndef _burn_type_H_
#define _burn_type_H_

#include <AMReX_REAL.H>

#include <network.H>
#include <eos_type.H>
#include <extern_parameters.H>

#include <ArrayUtilities.H>

using namespace amrex;

// A generic structure holding data necessary to do a nuclear burn.

// Set the number of independent variables -- this should be
// temperature, enuc + the number of species which participate
// in the evolution equations.

// For Strang evolution, this will be the number of equations we are
// evolving.  For simplified-SDC, we will need neqs when getting the
// reaction sources from the network's RHS.
const int neqs = 2 + NumSpec;

// Indices of the temperature and energy variables in the work arrays.

const int net_itemp = NumSpec + 1;
const int net_ienuc = NumSpec + 2;

// this is the data type that is used to get the ydots from the actual
// RHS of the network, regardless of Strang or SDC
typedef amrex::Array1D<Real, 1, neqs> YdotNetArray1D;

#ifdef SIMPLIFIED_SDC

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
// the following are not evolved
const int SFX   = SFS + NumSpec;
const int SRHO  = SFX + NumAux; // this is SFS + NumSpec if NumAux = 0;
const int SMX   = SRHO + 1;
const int SMY   = SRHO + 2;
const int SMZ   = SRHO + 3;

const int SVAR  = SMZ+1;
const int SVAR_EVOLVE = SFX;

#elif defined(SDC_EVOLVE_ENTHALPY)

// integrate rho*X (species masses) and rho*h (enthalpy)
// carry pressure for EOS calls in RHS

const int SFS = 0;
const int SENTH = SFS + NumSpec;
const int SVAR  = SENTH+1;
const int SVAR_EVOLVE = SVAR;

#endif

#endif

// this is the data type of the dense Jacobian that the network wants.
// it is not the same size as the Jacobian that VODE cares about when
// we are doing simplified-SDC
typedef ArrayUtil::MathArray2D<1, neqs, 1, neqs> JacNetArray2D;

struct burn_t
{

  // this first group are the quantities the network RHS uses

  Real rho;
  Real T;
  Real e;
  Real xn[NumSpec];
#if NAUX_NET > 0
  Real aux[NumAux];
#endif

  // now come the bits that we need for SDC or Strang evolution

#ifdef SIMPLIFIED_SDC

  // y is the input conserved state.  We will keep this state updated
  // in time as we integrate, such that upon output it will be the
  // final conserved state.
  Real y[SVAR];

  // we need to store a copy of the original state as well so we can
  // handle the non-evolved state evolution
  Real rho_orig;

#if defined(SDC_EVOLVE_ENERGY)
  Real umx_orig;
  Real umy_orig;
  Real umz_orig;
#endif

  // ydot_a are the advective terms that will modify the state y due
  // to hydrodynamics over the timestep.
  Real ydot_a[SVAR];

#if defined(SDC_EVOLVE_ENERGY)
  bool T_from_eden;
#elif defined(SDC_EVOLVE_ENTHALPY)
  // Pressure in case we wish to use it for EOS calls
  Real p0;
#endif

  int i;
  int j;
  int k;

  int sdc_iter;
  int num_sdc_iters;

#endif

#ifdef STRANG
  // Strang stuff

  Real cv;
  Real cp;

  // Last temperature we evaluated the EOS at
  Real T_old;
  Real cv_old;
  Real cp_old;

  // Temperature derivatives of specific heat
  Real dcvdT;
  Real dcpdT;

  // Are we doing a self-heating burn?
  bool self_heat;

#endif

#ifdef TRUE_SDC
  Real E_var;
  Real mom[3];
  Real f_source[NumSpec+2];
#endif

  Real y_e;
  Real eta;
  Real abar;
  Real zbar;

  // diagnostics
  int n_rhs, n_jac;

  // Was the burn successful?
  bool success;
};


inline
std::ostream& operator<< (std::ostream& o, burn_t const& burn_state)
{
  o << "rho = " << burn_state.rho << std::endl;
  o << "T =   " << burn_state.T << std::endl;
  o << "xn = ";
  for (int n = 0; n < NumSpec; ++n) {
    o << burn_state.xn[n] << " ";
  }
  o << std::endl;
#if NAUX_NET > 0
  o << "aux = ";
  for (int n = 0; n < NumAux; ++n) {
    o << burn_state.aux[n] << " ";
  }
  o << std::endl;
#endif

#ifdef SIMPLIFIED_SDC
  o << "y = ";
  for (int n = 0; n < SVAR; n++) {
      o << burn_state.y[n] << " ";
  }
  o << std::endl;

  o << "ydot_a = ";
  for (int n = 0; n < SVAR; n++) {
      o << burn_state.ydot_a[n] << " ";
  }
  o << std::endl;
#endif
  return o;
}


// Given an eos type, copy the data relevant to the burn type.

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eos_to_burn (const T& eos_state, burn_t& burn_state)
{
    burn_state.rho  = eos_state.rho;
    burn_state.T    = eos_state.T;
    burn_state.e    = eos_state.e;
    for (int n = 0; n < NumSpec; ++n) {
        burn_state.xn[n] = eos_state.xn[n];
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
        burn_state.aux[n] = eos_state.aux[n];
    }
#endif

#ifdef STRANG
    burn_state.cv   = eos_state.cv;
    burn_state.cp   = eos_state.cp;
#endif

    burn_state.y_e  = eos_state.y_e;
    burn_state.eta  = eos_state.eta;
    burn_state.abar = eos_state.abar;
    burn_state.zbar = eos_state.zbar;
}



// Given a burn type, copy the data relevant to the eos type. 
// Note that when doing simplified SDC integration, we should
// avoid using this interface because the energy includes a
// contribution from the advection term. However this is useful
// for instantaneous RHS evaluations.

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void burn_to_eos (const burn_t& burn_state, T& eos_state)
{
    eos_state.rho  = burn_state.rho;
    eos_state.T    = burn_state.T;
    eos_state.e    = burn_state.e;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = burn_state.xn[n];
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
        eos_state.aux[n] = burn_state.aux[n];
    }
#endif

#ifdef STRANG
    eos_state.cv   = burn_state.cv;
    eos_state.cp   = burn_state.cp;
#endif

    eos_state.y_e  = burn_state.y_e;
    eos_state.eta  = burn_state.eta;
    eos_state.abar = burn_state.abar;
    eos_state.zbar = burn_state.zbar;
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void normalize_abundances_burn (burn_t& state)
{

    Real sum = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        state.xn[n] = amrex::max(small_x, amrex::min(1.0_rt, state.xn[n]));
        sum += state.xn[n];
    }
    for (int n = 0; n < NumSpec; ++n) {
        state.xn[n] /= sum;
    }

}

#endif
