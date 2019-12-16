! Common variables and routines for burners
! that use BS for their integration.

module bs_integrator_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine bs_integrator_init()

    use bs_type_module, only: nseq

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    nseq = [2, 6, 10, 14, 22, 34, 50, 70]

    !$acc update device(nseq)

  end subroutine bs_integrator_init



  ! Main interface

  subroutine bs_integrator(state_in, state_out, dt, time, status)

    !$acc routine seq

    use extern_probin_module, only: burner_verbose
#if defined(SDC_EVOLVE_ENERGY)
    use sdc_type_module, only: sdc_t, SRHO, SEINT, SEDEN, SFS
#elif defined(SDC_EVOLVE_ENTHALPY)
    use sdc_type_module, only: sdc_t, SFS, SENTH
#endif
    use network, only: nspec
    use stiff_ode, only: ode, IERR_NONE
    use bs_type_module, only: bs_t, sdc_to_bs, bs_to_sdc
    use amrex_constants_module, only: ZERO
    use amrex_fort_module, only : rt => amrex_real
    use bs_rpar_indices, only : irp_t0
    use integration_data, only: integration_status_t

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    ! Input arguments

    type(sdc_t), intent(in   ) :: state_in
    type(sdc_t), intent(inout) :: state_out
    real(rt),  intent(in   ) :: dt, time
    type(integration_status_t), intent(inout) :: status

    ! Local variables
    integer :: ierr

    real(rt) :: t0, t1                   ! starting and ending time

    type (bs_t) :: bs

    ! BS does not allow for per-equation tolerances, so aggregate them here
    bs % atol(:) = 0.e0_rt
    bs % rtol(:) = max(status % rtol_spec, status % rtol_temp, status % rtol_enuc)

    ! Start out by assuming a successful burn.

    state_out % success = .true.

    ! Initialize the integration time.
    t0 = ZERO
    t1 = t0 + dt

    call sdc_to_bs(state_in, bs)

    bs % n_rhs = 0
    bs % n_jac = 0

    bs % self_heat = .false.

    ! set the time offset -- we integrate from 0 to dt, so this
    ! is the offset to simulation time
    
    bs % u(irp_t0) = time

    ! Call the integration routine.

    call ode(bs, t0, t1, maxval(bs % rtol), ierr)
    
    ! Store the final data

    call bs_to_sdc(state_out, bs)

    ! If we failed, print out the current state of the integration.

    if (ierr /= IERR_NONE) then

#ifndef CUDA
#if defined(SDC_EVOLVE_ENERGY)
       print *, 'ERROR: integration failed in net'
       print *, 'ierr = ', ierr
       print *, 'time = ', bs % t
       print *, 'dens start = ', state_in % y(SRHO)
       print *, 'eint start = ', state_in % y(SEINT) / state_in % y(SRHO)
       print *, 'xn start = ', state_in % y(SFS:SFS+nspec-1) / state_in % y(SRHO)
       print *, 'dens current = ', state_out % y(SRHO)
       print *, 'eint current = ', state_out % y(SEINT) / state_out % y(SRHO)
       print *, 'xn current = ', state_out % y(SFS:SFS+nspec-1) / state_out % y(SRHO)
       print *, 'energy generated = ', state_out % y(SEDEN) / state_out % y(SRHO) - &
            state_in % y(SEDEN) / state_in % y(SRHO)
#endif
#endif

       state_out % success = .false.
       return

    endif

  end subroutine bs_integrator

end module bs_integrator_module
