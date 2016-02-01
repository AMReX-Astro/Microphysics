
module rpar_indices

  implicit none

  integer, save :: n_rpar_comps = 0

  integer, save :: irp_dens, irp_cp, irp_dhdX, irp_Teos, irp_Tcrit, &
                   irp_rates, irp_drtdt, irp_Y56

contains

  function get_next_rpar_index(num) result (next)

    ! return the next starting index for a plotfile quantity,
    ! and increment the counter of plotfile quantities by num
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
    irp_Teos  = get_next_rpar_index(1)
    irp_Tcrit = get_next_rpar_index(1)

    irp_rates = get_next_rpar_index(nrat)
    irp_drtdt = get_next_rpar_index(nrat)

    irp_Y56 = get_next_rpar_index(1)

  end subroutine init_rpar_indices

end module rpar_indices
