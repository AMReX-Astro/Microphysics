! rpar are the real quantities that are passed through the VODE call to the
! RHS and Jacobian routines

module rpar_indices

  implicit none

  integer :: n_rpar_comps = 0

  integer :: irp_dens, irp_cv, irp_cp, irp_dedY, irp_dhdY
  integer :: irp_nspec
  integer :: irp_abar, irp_zbar
  integer :: irp_eta, irp_ye
  integer :: irp_self_heat
  integer :: irp_have_rates
  integer :: irp_rates
  integer :: irp_Told, irp_dcvdt, irp_dcpdt

contains

  function get_next_rpar_index(num) result (next)

    implicit none

    ! Return the next starting index for a plotfile quantity,
    ! and increment the counter of plotfile quantities by num.

    integer, intent(in) :: num
    integer             :: next

    if (num < 1) then

       next = -1

    else

       ! Return the value of the first index of this group of data.

       next = n_rpar_comps + 1

       ! Update our local record of how many variables we are using.

       n_rpar_comps = n_rpar_comps + num

    endif

  end function get_next_rpar_index


  subroutine init_rpar_indices()

    use network, only: nspec, nspec_evolve
    use actual_network_data, only: nrates
    use burn_type_module, only: num_rate_groups

    implicit none

    irp_dens      = get_next_rpar_index(1)
    irp_cp        = get_next_rpar_index(1)
    irp_cv        = get_next_rpar_index(1)
    irp_abar      = get_next_rpar_index(1)
    irp_zbar      = get_next_rpar_index(1)
    irp_eta       = get_next_rpar_index(1)
    irp_ye        = get_next_rpar_index(1)
    irp_self_heat = get_next_rpar_index(1)
    irp_nspec     = get_next_rpar_index(nspec - nspec_evolve)
    irp_dhdY      = get_next_rpar_index(nspec)
    irp_dedY      = get_next_rpar_index(nspec)
    irp_have_rates = get_next_rpar_index(1)
    irp_rates     = get_next_rpar_index(num_rate_groups * nrates)
    irp_Told      = get_next_rpar_index(1)
    irp_dcvdt     = get_next_rpar_index(1)
    irp_dcpdt     = get_next_rpar_index(1)

  end subroutine init_rpar_indices

end module rpar_indices
