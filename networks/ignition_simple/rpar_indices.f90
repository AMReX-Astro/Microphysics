! rpar are the real quantities that are passed through the VODE call to the
! RHS and Jacobian routines

module rpar_indices

  implicit none

  integer, save :: n_rpar_comps = 0

  integer, save :: irp_dens, irp_cp, irp_dhdX, irp_O16
  integer, save :: irp_rate, irp_dratedt, irp_sc1212, irp_dsc1212dt, irp_xc12tmp
  
contains

  function get_next_rpar_index(num) result (next)

    implicit none
    
    ! Return the next starting index for a plotfile quantity,
    ! and increment the counter of plotfile quantities by num.
    
    integer, intent(in) :: num
    integer             :: next

   ! Return the value of the first index of this group of data.

    next = n_rpar_comps + 1

    ! Update our local record of how many variables we are using.
    
    n_rpar_comps = n_rpar_comps + num

  end function get_next_rpar_index


  subroutine init_rpar_indices(nspec)

    implicit none
    
    integer, intent(in) :: nspec

    irp_dens = get_next_rpar_index(1)
    irp_cp   = get_next_rpar_index(1)
    irp_dhdX = get_next_rpar_index(nspec)
    irp_O16  = get_next_rpar_index(1)

    irp_rate      = get_next_rpar_index(1)
    irp_dratedt   = get_next_rpar_index(1)
    irp_sc1212    = get_next_rpar_index(1)
    irp_dsc1212dt = get_next_rpar_index(1)
    irp_xc12tmp   = get_next_rpar_index(1)    
    
  end subroutine init_rpar_indices

end module rpar_indices
