program test

  use bl_types
  use microphysics_math_module

  implicit none

  real(dp_t) :: a(2)

  a(1) = 1.0_dp_t
  a(2) = 2.e-16_dp_t

  print *, esum(a, 2)
  print *, a(1) + a(2)

end program test
