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
  real(rt), parameter :: SIX  = 6.0_rt
  real(rt), parameter :: EIGHT = 8.0_rt
  real(rt), parameter :: NINE = 9.0_rt
  real(rt), parameter :: TEN = 10.0_rt
  real(rt), parameter :: TWELVE = 12.0_rt
  real(rt), parameter :: HUN  = 100.0_rt
  real(rt), parameter :: THOU = 1000.0_rt

  real(rt), parameter :: THIRD = 1.0_rt/3.0_rt
  real(rt), parameter :: FIVE3RD = 5.0_rt/3.0_rt

  real(rt), parameter :: FOURTH = 1.0_rt/4.0_rt

  real(rt), parameter :: FIVE12TH = 5.0_rt/12.0_rt
  real(rt), parameter :: FIVE32ND = 5.0_rt/32.0_rt

  real(rt), parameter :: SIXTH = 1.0_rt / 6.0_rt

  real(rt), parameter :: M_PI = 3.141592653589793238462643383279502884197e0_rt

end module microphysics_type_module
