subroutine do_nse_F() bind (C, name="do_nse_F")

  ! this routine simply calls the table interpolation with a fixed
  ! input, so we can compare the results with those we get from the
  ! entire Microphysics machinery.

  use nse_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt), parameter :: temp = 1.e9_rt
  real(rt), parameter :: rho = 1.e9_rt
  real(rt), parameter :: ye = 0.46_rt

  real(rt) :: abar, dq, dyedt, X(nspec)

  call nse_interp(temp, rho, ye, abar, dq, dyedt, X)

  print *, "temp:    ", temp
  print *, "rho:     ", rho
  print *, "ye:      ", ye
  print *, "abar:    ", abar
  print *, "be/a:    ", dq
  print *, "dyedt:   ", dyedt
  print *, "X:       ", X(:)

end subroutine do_nse_F
