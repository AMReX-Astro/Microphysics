module nse_check_module

  use network
  use burn_type_module, only : burn_t
  use extern_probin_module, only : rho_nse, T_nse, C_nse, He_Fe_nse

  implicit none

contains

  subroutine in_nse(burn_state, nse_check)

    type(burn_t), intent(in) :: burn_state
    integer, intent(out) :: nse_check

    real(rt) Fe_group, C_group, He_group

    !$gpu

    nse_check = 0

    if (burn_state % rho > rho_nse .and. burn_state % T > T_nse) then

       ! Ma checks on Fe-group (for our composition, this means Cr48, Fe52, Fe54, Ni56)
       ! and C-group (for us, that is C12, N14)
       ! and He-group (for us, that is H1, He3, He4)

       Fe_group = burn_state % xn(icr48) + burn_state % xn(ife52) + &
                  burn_state % xn(ife54) + burn_state % xn(ini56)
       C_group = burn_state % xn(ic12) + burn_state % xn(in14)
       He_group = burn_state % xn(ih1) + burn_state % xn(ihe3) + burn_state % xn(ihe4)

       if (Fe_group + He_group > He_Fe_nse .and. C_group < C_nse) then
          nse_check = 1
       end if
    end if

  end subroutine in_nse

end module nse_check_module
