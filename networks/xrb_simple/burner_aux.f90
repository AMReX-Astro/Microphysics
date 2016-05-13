module rpar_indices

  implicit none

  integer, save :: n_rpar_comps = 0

  integer, save :: irp_dens, irp_cp, irp_dhdX, irp_T9_eos, irp_dTcrit,&
                   irp_rates, irp_drtdt

contains

  function get_next_rpar_index(num) result (next)

    ! return the next starting index for a quantity,
    ! and increment the counter of quantities by num
    integer :: num, next

    next = n_rpar_comps + 1
    n_rpar_comps = n_rpar_comps + num

    return
  end function get_next_rpar_index


  subroutine init_rpar_indices(nrat, nspec)

    integer, intent(in) :: nrat, nspec

    irp_dens  = get_next_rpar_index(1)
    irp_cp    = get_next_rpar_index(1)
    irp_dhdX  = get_next_rpar_index(nspec)
    irp_T9_eos = get_next_rpar_index(1)
    irp_dTcrit = get_next_rpar_index(1)

    !============================================================
    ! only rate-related things below this point
    irp_rates = get_next_rpar_index(nrat)
    irp_drtdt = get_next_rpar_index(nrat)

  end subroutine init_rpar_indices

end module rpar_indices
