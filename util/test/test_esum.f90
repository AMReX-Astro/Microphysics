program test

  use amrex_fort_module, only : rt => amrex_real

  use microphysics_math_module

  implicit none

  real(rt) :: a(2)
  real(rt) :: b(3)

  a(1) = 1.0_rt
  a(2) = 2.e-16_rt

  b(1) = 1.0_rt
  b(2) = 3e-16_rt
  b(3) = 1.0_rt

  print *, esum2(a)
  print *, a(1) + a(2)

  print *, esum3(b)
  print *, b(1) + b(2) + b(3)

end program test
