module cj_det_module

  use amrex_fort_module, only : rt => amrex_real

  use eos_type_module
  use eos_module

  implicit none

contains

  subroutine adiabat(eos_state_fuel, eos_state_ash, q, istatus)

    implicit none

    type(eos_t), intent(in) :: eos_state_fuel
    type(eos_t), intent(inout) :: eos_state_ash
    real(rt), intent(in) :: q
    integer, intent(out) :: istatus

    real(rt), parameter :: tol = 1.e-8_rt
    logical :: converged
    real(rt) :: f, dfdT, dT
    integer :: iter
    integer, parameter :: max_iter = 50

    ! we want to zero e_1 + q - e_2 + 0.5*(p_1 + p_2)*(v_1 - v_2)
    ! where v = 1/rho

    ! we have a rho_2 (the ash input), so we need to find the T that
    ! makes this all work

    ! we assume that we come in with a reasonable guess for T

    converged = .false.

    iter = 0
    do while (.not. converged .and. iter < max_iter)
       call eos(eos_input_rt, eos_state_ash)

       f = eos_state_fuel % e + q - eos_state_ash % e + &
            0.5_rt * (eos_state_fuel % p + eos_state_ash % p) * &
            (1.0_rt/eos_state_fuel % rho - 1.0_rt/eos_state_ash % rho)

       dfdT = -eos_state_ash % dedT + 0.5_rt * eos_state_ash % rt * &
            (1.0_rt/eos_state_fuel % rho - 1.0_rt/eos_state_ash % rho)

       dT = -f/dfdT

       if (abs(dT) < tol * eos_state_ash % T) then
          converged = .true.
       endif

       eos_state_ash % T = eos_state_ash % T + dT

       iter = iter + 1

    enddo

    if (.not. converged) then
       istatus = -1
    else
       istatus = 0
    endif

  end subroutine adiabat

  subroutine cj_cond(eos_state_fuel, eos_state_ash, q)

    implicit none

    type(eos_t), intent(in) :: eos_state_fuel
    type(eos_t), intent(inout) :: eos_state_ash
    real(rt), intent(in) :: q

    real(rt), parameter :: tol = 1.e-8_rt
    logical :: converged
    real(rt) :: rho_old, drho
    integer :: iter, istatus
    integer, parameter :: max_iter = 50

    ! iterate, picking the density that corresponds to the CJ point
    call eos(eos_input_rt, eos_state_ash)

    drho = 1.e30_rt

    ! this is the density we find from the tangent point to the
    ! Hugoniot curve
    eos_state_ash % rho = eos_state_fuel % rho * &
         (1.0_rt + (eos_state_ash % p - eos_state_fuel % p) / &
                      (eos_state_ash % gam1 * eos_state_ash % p))

    iter = 0
    converged = .false.
    do while (.not. converged .and. iter < max_iter)

       rho_old = eos_state_ash % rho

       call adiabat(eos_state_fuel, eos_state_ash, q, istatus)

       ! this is the density we find from the tangent point to the
       ! Hugoniot curve
       eos_state_ash % rho = eos_state_fuel % rho * &
            (1.0_rt + (eos_state_ash % p - eos_state_fuel % p) / &
                         (eos_state_ash % gam1 * eos_state_ash % p))

       drho = eos_state_ash % rho - rho_old

       if (abs(drho) < tol * eos_state_ash % rho) then
          converged = .true.
       endif

       iter = iter + 1
    enddo

    if (.not. converged .or. istatus == -1) then
       call bl_error("CJ did not converge")
    endif

  end subroutine cj_cond

end module cj_det_module
