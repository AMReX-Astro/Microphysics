! rpar are the real quantities that are passed through the VODE call to the
! RHS and Jacobian routines

module rpar_indices

  implicit none

  integer :: n_rpar_comps = 0

  integer :: irp_dens, irp_cv, irp_cp, irp_dedX, irp_dhdX, irp_smallx
  integer :: irp_abar, irp_zbar
  integer :: irp_dydt, irp_dratesdt

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


  subroutine init_rpar_indices(nrates, nspec)

    implicit none
    
    integer, intent(in) :: nrates, nspec

    irp_dens     = get_next_rpar_index(1)
    irp_cp       = get_next_rpar_index(1)
    irp_dhdX     = get_next_rpar_index(nspec)
    irp_smallx   = get_next_rpar_index(1)
    irp_cv       = get_next_rpar_index(1)
    irp_abar     = get_next_rpar_index(1)
    irp_zbar     = get_next_rpar_index(1)
    irp_dedX     = get_next_rpar_index(nspec)
    irp_dydt     = get_next_rpar_index(nspec)
    irp_dratesdt = get_next_rpar_index(nrates)

  end subroutine init_rpar_indices

end module rpar_indices
