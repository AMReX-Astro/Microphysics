module vode_parameters_module

  use burn_type_module, only: neqs

  implicit none

  integer, parameter :: VODE_NEQS = neqs

  ! Our problem is stiff, so tell ODEPACK that. 21 means stiff, jacobian
  ! function is supplied; 22 means stiff, figure out my jacobian through
  ! differencing.

#ifdef CUDA
  ! Negative method flags mean on the GPU we turn off Jacobian caching
  ! to reduce our memory requirements.
  integer, parameter :: MF_ANALYTIC_JAC = -21, MF_NUMERICAL_JAC = -22
#else
  integer, parameter :: MF_ANALYTIC_JAC = 21, MF_NUMERICAL_JAC = 22
#endif

  ! Tolerance parameters:
  !
  !  itol specifies whether to use an single absolute tolerance for
  !  all variables (1), or to pass an array of absolute tolerances, one
  !  for each variable with a scalar relative tol (2), a scalar absolute
  !  and array of relative tolerances (3), or arrays for both (4).
  !
  !  The error is determined as e(i) = rtol*abs(y(i)) + atol, and must
  !  be > 0.  Since we have some compositions that may be 0 initially,
  !  we will specify both an absolute and a relative tolerance.
  !
  ! We will use arrays for both the absolute and relative tolerances,
  ! since we want to be easier on the temperature than the species.

  integer, parameter :: VODE_ITOL = 4

  ! We want to do a normal computation, and get the output values of y(t)
  ! after stepping though dt.

  integer, PARAMETER :: ITASK = 1

  ! We will override the maximum number of steps, so turn on the
  ! optional arguments flag.

  integer, parameter :: IOPT = 1

  ! Declare a real work array of size 22 + 9*NEQ + 2*NEQ**2 and an
  ! integer work array of size 30 + NEQ. These are VODE constants
  ! that depend on the integration mode we're using -- see dvode.f.

  integer, parameter :: VODE_LRW = 22 + 9*neqs + 2*neqs**2
  integer, parameter :: VODE_LIW = 30 + neqs

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
  integer, parameter :: VODE_LENWM = 2 + 2 * neqs**2

end module vode_parameters_module
