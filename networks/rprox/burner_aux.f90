module rpar_indices

  implicit none

  integer :: n_rpar_comps = 0

  integer :: irp_dens, irp_cp, irp_dhdX, irp_T9_eos, irp_dTcrit,&
             irp_rates, irp_drtdt, &
             irp_dlambCNOdh1, irp_drs1dhe4, irp_drr1dh1, &
             irp_dlambda1dhe4, irp_dlambda2dhe4, &
             irp_delta1, irp_delta2, irp_r56eff, irp_dr56effdt

contains

  function get_next_rpar_index(num) result (next)

    ! return the next starting index for a plotfile quantity,
    ! and increment the counter of plotfile quantities by num
    integer :: num, next

    next = n_rpar_comps + 1
    n_rpar_comps = n_rpar_comps + num

  end function get_next_rpar_index


  subroutine init_rpar_indices(nrat, nspec)

    integer, intent(in) :: nrat, nspec

    irp_dens  = get_next_rpar_index(1)
    irp_cp    = get_next_rpar_index(1)
    irp_dhdX  = get_next_rpar_index(nspec)
    irp_T9_eos = get_next_rpar_index(nspec)
    irp_dTcrit = get_next_rpar_index(nspec)

    !============================================================
    ! only rate-related things below this point
    irp_rates = get_next_rpar_index(nrat)
    irp_drtdt = get_next_rpar_index(nrat)

    irp_dlambCNOdh1 = get_next_rpar_index(1)
    irp_drs1dhe4 = get_next_rpar_index(1)
    irp_drr1dh1 = get_next_rpar_index(1)
    irp_dlambda1dhe4 = get_next_rpar_index(1)
    irp_dlambda2dhe4 = get_next_rpar_index(1)

    irp_delta1 = get_next_rpar_index(1)
    irp_delta2 = get_next_rpar_index(1)

    irp_r56eff = get_next_rpar_index(1)
    irp_dr56effdt = get_next_rpar_index(1)

  end subroutine init_rpar_indices

end module rpar_indices
