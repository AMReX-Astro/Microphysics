module network_indices

  implicit none

  integer, parameter :: ic12_ = 1
  integer, parameter :: io16_ = 2
  integer, parameter :: iash_ = 3

end module network_indices

module rpar_indices

  implicit none

  integer, save :: n_rpar_comps = 0

  integer, save :: irp_dens, irp_cp, irp_dhdx, irp_o16
  integer, save :: irp_rate, irp_dratedt, irp_sc1212, irp_dsc1212dt, irp_xc12tmp

contains

  function get_next_rpar_index(num) result (next)

    ! return the next starting index for a plotfile quantity,
    ! and increment the counter of plotfile quantities by num
    integer :: num, next

    next = n_rpar_comps + 1
    n_rpar_comps = n_rpar_comps + num

  end function get_next_rpar_index


  subroutine init_rpar_indices(nspec)

    integer, intent(in) :: nspec

    irp_dens  = get_next_rpar_index(1)
    irp_cp    = get_next_rpar_index(1)
    irp_dhdX  = get_next_rpar_index(nspec)
    irp_o16   = get_next_rpar_index(1)

    irp_rate      = get_next_rpar_index(1)
    irp_dratedt   = get_next_rpar_index(1)
    irp_sc1212    = get_next_rpar_index(1)
    irp_dsc1212dt = get_next_rpar_index(1)
    irp_xc12tmp   = get_next_rpar_index(1)

  end subroutine init_rpar_indices

end module rpar_indices
