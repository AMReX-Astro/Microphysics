module microphysics_type_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt), parameter :: ZERO = 0.0_rt
  real(rt), parameter :: HALF = 0.5_rt
  real(rt), parameter :: ONE = 1.0_rt
  real(rt), parameter :: TWO = 2.0_rt
  real(rt), parameter :: THREE = 3.0_rt
  real(rt), parameter :: FOUR = 4.0_rt
  real(rt), parameter :: FIVE = 5.0_rt

  real(rt), parameter :: FIVE3RD = 5.0_rt/3.0_rt

end module microphysics_type_module
