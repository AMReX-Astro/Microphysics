module react_zones_module

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module

  implicit none

contains

  subroutine init_state(lo, hi, &
                        state, s_lo, s_hi, ncomp, npts) bind(C, name="init_state")

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), ncomp)
    integer, intent(in) :: npts, ncomp

    integer :: i, j, k, n

    n = 0
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             state(i, j, k, 1) = ONE - real(n, kind=rt)/real(npts**3, kind=rt)
             state(i, j, k, 2) = ZERO
             state(i, j, k, 3) = ZERO

             n = n + 1
          enddo
       enddo
    enddo

  end subroutine init_state


  subroutine do_react(lo, hi, &
                      state, s_lo, s_hi, ncomp, dt) bind(C, name="do_react")

    use cuvode_types_module, only: dvode_t
    use cuvode_module, only: dvode

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), ncomp)
    integer, intent(in), value :: ncomp
    real(rt), intent(in), value :: dt

    integer         :: ii, jj, kk, n

    ! VODE variables
    type (dvode_t) :: dvode_state

    ! istate determines the state of the calculation.  A value of 1 meeans
    ! this is the first call to the problem -- this is what we will want.
    
    integer :: istate
    

    !$gpu

    do ii = lo(1), hi(1)
       do jj = lo(2), hi(2)
          do kk = lo(3), hi(3)

             ! Use an analytic Jacobian
             dvode_state % jacobian = 1

             ! Set the absolute tolerances
             dvode_state % atol(1) = 1.e-8_rt
             dvode_state % atol(2) = 1.e-14_rt
             dvode_state % atol(3) = 1.e-6_rt

             ! Set the relative tolerances
             dvode_state % rtol(1) = 1.e-4_rt
             dvode_state % rtol(2) = 1.e-4_rt
             dvode_state % rtol(3) = 1.e-4_rt

             ! We want VODE to re-initialize each time we call it.
             dvode_state % istate = 1

             ! Take no more than 500 steps.
             dvode_state % MXSTEP = 500

             ! Initialize the integration time and set the final time to dt
             dvode_state % T = ZERO
             dvode_state % TOUT = dt

             ! Initialize the initial conditions
             do n = 1, ncomp
                dvode_state % y(n) = state(ii, jj, kk, n)
             enddo

             ! Call the integration routine.
             call dvode(dvode_state)

             ! Check if the integration failed
             if (dvode_state % istate < 0) then
#ifndef AMREX_USE_CUDA       
                print *, 'ERROR: integration failed'
                print *, 'istate = ', dvode_state % istate
                print *, 'time = ', dvode_state % T
                print *, 'Y start = ', state(ii, jj, kk, :)
                print *, 'Y current = ', dvode_state % y
#endif
                stop
             endif
             
             ! Store the final result
             do n = 1, ncomp
                state(ii, jj, kk, n) = dvode_state % y(n)
             enddo

          enddo
       enddo
    enddo

  end subroutine do_react

end module react_zones_module
