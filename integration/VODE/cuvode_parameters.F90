module cuvode_parameters_module

  use burn_type_module, only: neqs
#ifdef SIMPLIFIED_SDC
  use sdc_type_module, only : SVAR_EVOLVE
#endif
  use network, only : nspec

  use amrex_fort_module, only : rt => amrex_real
  implicit none

#ifdef TRUE_SDC
  integer, parameter :: VODE_NEQS = nspec + 2
#elif SIMPLIFIED_SDC
  integer, parameter :: VODE_NEQS = SVAR_EVOLVE
#else
  integer, parameter :: VODE_NEQS = neqs
#endif

  ! Our problem is stiff, so tell ODEPACK that. 21 means stiff, jacobian
  ! function is supplied; 22 means stiff, figure out my jacobian through
  ! differencing.

  ! Negative method flags mean on the GPU we turn off Jacobian caching
  ! to reduce our memory requirements.
  integer, parameter :: MF_ANALYTIC_JAC_NOCACHE = -21, MF_NUMERICAL_JAC_NOCACHE = -22
  integer, parameter :: MF_ANALYTIC_JAC_CACHED = 21, MF_NUMERICAL_JAC_CACHED = 22

  ! For VODE, LMAX = MAXORD + 1, so the following are specific
  ! to our choice of method (see the dvode README for details)

  ! We are using the backward-differentiation formulas (BDF-s)
  ! and the maximum order for BDF mode should be no greater than 5.
  ! If the implicit Adams method is desired, the maximum order
  ! should be no greater than 12. This is the VODE_MAXORD variable.

  ! IMPORTANT:
  ! Input sanitization in VODE has been removed for this for
  ! performance reasons so make sure the following parameter
  ! is set correctly.
  integer, parameter :: VODE_MAXORD = 5

  integer, parameter :: VODE_LMAX = VODE_MAXORD + 1

end module cuvode_parameters_module
