PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE
USE_CUDA   = FALSE
USE_ACC    = FALSE

USE_COMPILE_PIC = TRUE
USE_COMPILE_WITH_F2PY = TRUE
SKIP_LINKING = TRUE

USE_REACT = TRUE
USE_EXTRA_THERMO = TRUE

EBASE =

# define the location of the Microphysics top directory
MICROPHYSICS_HOME  := ../

# This sets the EOS directory in Microphysics/EOS -- note: gamma_law will not work,
# you'll need to use gamma_law_general
EOS_DIR     := helmholtz

# This sets the network directory in Microphysics/networks
NETWORK_DIR := aprox13

CONDUCTIVITY_DIR := stellar

# This sets the integrator directory in Microphysics/integration
INTEGRATOR_DIR := VODE90

EXTERN_SEARCH += .

Bpack   := Make.package
Blocs   := .

include $(MICROPHYSICS_HOME)/Make.Microphysics

