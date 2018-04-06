module cj_det_module

  use bl_types, only: dp_t
  use eos_type_module
  use eos_module

  implicit none

contains

  subroutine adiabat(eos_state_fuel, eos_state_ash, q)

    implicit none

    type(eos_t), intent(in) :: eos_state_fuel
    type(eos_t), intent(inout) :: eos_state_ash
    real(dp_t), intent(in) :: q

    real(dp_t), parameter :: tol = 1.e-8_dp_t
    logical :: converged
    real(dp_t) :: f, dfdT, dT
    integer :: iter
    integer, parameter :: max_iter = 25

    ! we want to zero e_1 + q - e_2 + 0.5*(p_1 + p_2)*(v_1 - v_2)
    ! where v = 1/rho

    ! we have a rho_2 (the ash input), so we need to find the T that
    ! makes this all work

    ! initial T guess
    eos_state_ash % T = eos_state_fuel % T

    converged = .false.

    iter = 0
    do while (.not. converged .and. iter < max_iter)
       call eos(eos_input_rt, eos_state_ash)

       f = eos_state_fuel % e + q - eos_state_ash % e + &
            0.5_dp_t * (eos_state_fuel % p + eos_state_ash % p) * &
            (1.0_dp_t/eos_state_fuel % rho - 1.0_dp_t/eos_state_ash % rho)

       dfdT = -eos_state_ash % dedT + 0.5_dp_t * eos_state_ash % dpdT * &
            (1.0_dp_t/eos_state_fuel % rho - 1.0_dp_t/eos_state_ash % rho)

       dT = -f/dfdT
       print *, eos_state_ash % T, abs(dT)

       if (abs(dT) < tol * eos_state_ash % T) then
          converged = .true.
       endif

       eos_state_ash % T = eos_state_ash % T + dT

       iter = iter + 1

    enddo

    if (.not. converged) then
       call bl_error("Error: non-convergence")
    endif

  end subroutine adiabat

end module cj_det_module
