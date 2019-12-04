module cuvode_nordsieck_module

  use cuvode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_LIW,   &
                                    VODE_LENWM, VODE_MAXORD, VODE_ITOL
  use cuvode_types_module, only: dvode_t, rwork_t

  use cuvode_constants_module

  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine advance_nordsieck(rwork, vstate)

    ! Effectively multiplies the Nordsieck history
    ! array by the Pascal triangle matrix.

    implicit none

    ! Declare arguments
    type(rwork_t), intent(inout) :: rwork
    type(dvode_t), intent(in   ) :: vstate

    ! Declare local variables
    integer :: k, j, i

    !$gpu

    do k = vstate % NQ, 1, -1
       do j = k, vstate % NQ
          do i = 1, VODE_NEQS
             rwork % YH(i, j) = rwork % YH(i, j) + rwork % YH(i, j+1)
          end do
       end do
    end do

  end subroutine advance_nordsieck


#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine retract_nordsieck(rwork, vstate)

    ! Undoes the Pascal triangle matrix multiplication
    ! implemented in subroutine advance_nordsieck.

    implicit none

    ! Declare arguments
    type(rwork_t), intent(inout) :: rwork
    type(dvode_t), intent(in   ) :: vstate

    ! Declare local variables
    integer :: k, j, i

    !$gpu

    do k = vstate % NQ, 1, -1
       do j = k, vstate % NQ
          do i = 1, VODE_NEQS
             rwork % YH(i, j) = rwork % YH(i, j) - rwork % YH(i, j+1)
          end do
       end do
    end do

  end subroutine retract_nordsieck

end module cuvode_nordsieck_module
