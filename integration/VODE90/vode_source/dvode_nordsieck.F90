module dvode_nordsieck_module

  use vode_type_module, only: rwork_t
  use vode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_LIW,   &
                                    VODE_LENWM, VODE_MAXORD, VODE_ITOL
  use dvode_type_module, only: dvode_t

  use dvode_constants_module

  implicit none

contains

  AMREX_DEVICE subroutine advance_nordsieck(rwork, vstate)

    ! Effectively multiplies the Nordsieck history
    ! array by the Pascal triangle matrix.

    implicit none

    ! Declare arguments
    type(rwork_t), intent(inout) :: rwork
    type(dvode_t), intent(in   ) :: vstate

    ! Declare local variables
    integer :: k, j, i

    do k = vstate % NQ, 1, -1
       do j = k, vstate % NQ
          do i = 1, VODE_NEQS
             rwork % YH(i, j) = rwork % YH(i, j) + rwork % YH(i, j+1)
          end do
       end do
    end do

  end subroutine advance_nordsieck


  AMREX_DEVICE subroutine retract_nordsieck(rwork, vstate)

    ! Undoes the Pascal triangle matrix multiplication
    ! implemented in subroutine advance_nordsieck.

    implicit none

    ! Declare arguments
    type(rwork_t), intent(inout) :: rwork
    type(dvode_t), intent(in   ) :: vstate

    ! Declare local variables
    integer :: k, j, i

    do k = vstate % NQ, 1, -1
       do j = k, vstate % NQ
          do i = 1, VODE_NEQS
             rwork % YH(i, j) = rwork % YH(i, j) - rwork % YH(i, j+1)
          end do
       end do
    end do

  end subroutine retract_nordsieck

end module dvode_nordsieck_module
