! A variable-sized stepsize driver for the BS integrator,
! based on the NR routine


module stiff_ode

  use bl_types

  implicit none

  real(kind=dp_t), parameter :: dt_min = 1.d-20
  real(kind=dp_t), parameter :: dt_ini = 1.d-10
  real(kind=dp_t), parameter :: SMALL = 1.d-30

  integer, parameter :: MAX_STEPS = 10000

  ! error codes
  integer, parameter :: IERR_DT_TOO_SMALL = -100

contains

  subroutine ode(yinit, neq, t, tmax, eps, f_rhs)
    ! integrate from t to tmax

    integer, intent(in) :: neq
    real(kind=dp_t), intent(inout) :: yinit(neq)
    real(kind=dp_t), intent(inout) :: t
    real(kind=dp_t), intent(in) :: tmax
    real(kind=dp_t), intent(in) :: eps

    external f_rhs

  
    real(kind=dp_t) :: y(neq), yscal(neq)


    ! initialize
    y(:) = yinit(:)
    dt = dtini

    do n = 1, MAX_STEPS

       ! get the scaling
       call f_rhs(t, y, dydt)
       yscal(:) =  abs(y(:)) + abs(dt*dydt(:)) + SMALL
       
       ! make sure we don't overshoot the ending time
       if (t + dt > tmax) dt = tmax - t
       
       ! take a step
       call single_step()
       
       ! finished?
       if (t - tmax >= ZERO) then
          yinit(:) = y(:)
          ierr = 0
          exit
       endif
       
       dt = dt_next

       if (dt < dt_min) then
          ierr = ERR_DT_TOO_SMALL
          exit
       endif

    enddo
    
  end subroutine ode

end module stiff_ode
  
  

  
