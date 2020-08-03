subroutine do_nse_F() bind (C, name="do_nse_F")

  ! this routine simply calls the table interpolation with a fixed
  ! input, so we can compare the results with those we get from the
  ! entire Microphysics machinery.

  use nse_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt), parameter :: t9 = 1.0_rt
  real(rt), parameter :: rho = 1.e9_rt
  real(rt), parameter :: ye = 0.46

  real(rt) :: abar, dq, dyedt, X(nspec)

  call nse_interp(t9, rho, ye, abar, dq, dyedt, X)

  print *, "      t9    ", "        rho    ", "       ye     ", &
       "     abar     ", "      be/a    ", "      dyedt  "

  write (*, 41) t9, rho, ye, abar, dq, dyedt, X(:)
41 format (1pe12.3, 5e14.5, 19e14.5)

end subroutine do_nse_F
