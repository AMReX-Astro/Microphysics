module dvode_constants_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt), parameter :: ZERO = 0.0D0
  real(rt), parameter :: ONE  = 1.0D0
  real(rt), parameter :: HALF = 0.5D0
  real(rt), parameter :: TWO  = 2.0D0
  real(rt), parameter :: FOUR = 4.0D0
  real(rt), parameter :: SIX  = 6.0D0
  real(rt), parameter :: HUN  = 100.0D0
  real(rt), parameter :: THOU = 1000.0D0

end module dvode_constants_module
