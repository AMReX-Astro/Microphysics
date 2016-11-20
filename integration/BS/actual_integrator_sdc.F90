! Common variables and routines for burners
! that use BS for their integration.

module actual_integrator_module

  implicit none

contains

  subroutine actual_integrator_init()

    use bs_type_module, only: nseq

    implicit none

    nseq = [2, 6, 10, 14, 22, 34, 50, 70]

    !$acc update device(nseq)

  end subroutine actual_integrator_init



  ! Main interface

  subroutine actual_integrator(state_in, state_out, dt, time)

    !$acc routine seq

    use extern_probin_module, only: burner_verbose, &
                                    rtol_spec, rtol_temp, rtol_enuc, &
                                    atol_spec, atol_temp, atol_enuc
    use sdc_type_module, only: sdc_t
    use stiff_ode, only: ode
    use bs_type_module, only: bs_t, sdc_to_bs, bs_to_sdc
    use bl_types, only: dp_t
    use bl_constants_module, only: ZERO

    implicit none

    ! Input arguments

    type(sdc_t), intent(in   ) :: state_in
    type(sdc_t), intent(inout) :: state_out
    real(dp_t),  intent(in   ) :: dt, time

    ! Local variables
    integer :: ierr

    real(kind=dp_t) :: t0, t1                   ! starting and ending time

    type (bs_t) :: bs

    real(dp_t) :: retry_change_factor

    ! BS does not allow for per-equation tolerances, so aggregate them here
    bs % atol(:) = 0.d0
    bs % rtol(:) = max(rtol_spec, rtol_temp, rtol_enuc)

    ! Initialize the integration time.
    t0 = ZERO
    t1 = t0 + dt

    call sdc_to_bs(state_in, bs)

    bs % n_rhs = 0
    bs % n_jac = 0

    bs % self_heat = .false.

    ! Call the integration routine.

    call ode(bs, t0, t1, maxval(bs % rtol), ierr)
    
    ! Store the final data

    call bs_to_sdc(state_out, bs)

  end subroutine actual_integrator

end module actual_integrator_module
