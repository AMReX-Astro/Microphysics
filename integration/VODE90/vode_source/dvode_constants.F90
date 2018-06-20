module dvode_constants_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), parameter :: ZERO = 0.0_rt
  real(rt), parameter :: ONE  = 1.0_rt
  real(rt), parameter :: HALF = 0.5_rt
  real(rt), parameter :: TWO  = 2.0_rt
  real(rt), parameter :: FOUR = 4.0_rt
  real(rt), parameter :: SIX  = 6.0_rt
  real(rt), parameter :: HUN  = 100.0_rt
  real(rt), parameter :: THOU = 1000.0_rt

end module dvode_constants_module
