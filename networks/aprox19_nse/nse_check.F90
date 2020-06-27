module nse_check_module

  use network

  implicit none

contains

  subroutine in_nse(rho, T, X, nse_check)

    real(rt), intent(in) :: rho, T, X(nse)
    integer, intent(out) :: nse_check

    real(rt) Fe_group, C_group, He_group

    nse_check = 0

    if (rho > rho_nse .and. T > T_nse) then

       ! Ma checks on Fe-group (for our composition, this means Cr48, Fe52, Fe54, Ni56)
       ! and C-group (for us, that is C12, N14)
       ! and He-group (for us, that is H1, He3, He4)

       Fe_group = X(icr48) + X(ife52) + X(ife54) + X(ini56)
       C_group = X(ic12) + X(in14)
       He_group = X(ih1) + X(ihe3) + X(ihe4)

       if (Fe_group + He_group > 0.88 .and. C_group < 0.01) then
          nse_check = 1
       end if
    end if

  end subroutine in_nse

end module nse_check_module
