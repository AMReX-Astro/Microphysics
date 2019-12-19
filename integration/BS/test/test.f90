program main

  use amrex_fort_module, only : rt => amrex_real
  use stiff_ode

  implicit none

  integer, parameter :: neq = 3

  real(kind=rt) :: y(neq), t, tmax, eps
  integer :: ierr
  
  external f_rhs

  t = 0.0_rt
  
  y(:) = [1.0_rt, 0.0_rt, 0.0_rt]

  tmax = 4.e10_rt
  eps = 1.e-10_rt

  call ode(y, neq, t, tmax, eps, f_rhs, ierr)

  print *, y
  print *, ierr

end program main
