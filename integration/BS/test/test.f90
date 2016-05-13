program main

  use bl_types
  use stiff_ode

  implicit none

  integer, parameter :: neq = 3

  real(kind=dp_t) :: y(neq), t, tmax, eps
  integer :: ierr
  
  external f_rhs

  t = 0.0
  
  y(:) = [1.0_dp_t, 0.0_dp_t, 0.0_dp_t]

  tmax = 4.e10_dp_t
  eps = 1.e-10_dp_t

  call ode(y, neq, t, tmax, eps, f_rhs, ierr)

  print *, y
  print *, ierr

end program main
