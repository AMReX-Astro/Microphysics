program test

  use bl_types
  use microphysics_math_module

  implicit none

  real(dp_t) :: a(2)
  real(dp_t) :: b(3)

  a(1) = 1.0_dp_t
  a(2) = 2.e-16_dp_t

  b(1) = 1.0_dp_t
  b(2) = 3e-16_dp_t
  b(3) = 1.0_dp_t

  print *, esum(a, 2)
  print *, a(1) + a(2)

  print *, esum(b, 3)
  print *, b(1) + b(2) + b(3)

end program test
