module cuvode_parameters_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  ! User-definable number of equations
  integer, parameter :: VODE_NEQS = 0

  ! Our problem is stiff, so tell ODEPACK that. 21 means stiff, jacobian
  ! function is supplied; 22 means stiff, figure out my jacobian through
  ! differencing.

#ifdef AMREX_USE_CUDA
  ! Negative method flags mean on the GPU we turn off Jacobian caching
  ! to reduce our memory requirements.
  integer, parameter :: MF_ANALYTIC_JAC = -21, MF_NUMERICAL_JAC = -22
#else
  integer, parameter :: MF_ANALYTIC_JAC = 21, MF_NUMERICAL_JAC = 22
#endif

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
