module cuvode_parameters_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  ! User-definable number of equations
  integer, parameter :: VODE_NEQS = 3

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
