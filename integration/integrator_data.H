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

#endif
