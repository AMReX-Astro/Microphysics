module dvode_constants_module

  use bl_types, only: dp_t

  implicit none

  real(dp_t), parameter :: ZERO = 0.0D0
  real(dp_t), parameter :: ONE  = 1.0D0
  real(dp_t), parameter :: HALF = 0.5D0
  real(dp_t), parameter :: TWO  = 2.0D0
  real(dp_t), parameter :: FOUR = 4.0D0
  real(dp_t), parameter :: SIX  = 6.0D0
  real(dp_t), parameter :: HUN  = 100.0D0
  real(dp_t), parameter :: THOU = 1000.0D0

end module dvode_constants_module
