#ifndef INTEGRATOR_DATA_H
#define INTEGRATOR_DATA_H

#include <burn_type.H>

// Define the size of the ODE system that VODE will integrate

#ifdef STRANG
const int INT_NEQS = NumSpec + 1;
#endif

#ifdef SIMPLIFIED_SDC
const int INT_NEQS = SVAR_EVOLVE;
#endif

#ifdef TRUE_SDC
const int INT_NEQS = NumSpec + 2;
#endif

typedef amrex::Array1D<short, 1, INT_NEQS> IArray1D;
typedef amrex::Array1D<Real, 1, INT_NEQS> RArray1D;
typedef ArrayUtil::MathArray2D<1, INT_NEQS, 1, INT_NEQS> RArray2D;

#endif
