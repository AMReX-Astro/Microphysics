! This is the interface to the burner for the simplified SDC case.

module actual_integrator_module

  use amrex_error_module
  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module

  use sdc_type_module
  use vode_type_module

  use cuvode_parameters_module

  implicit none

contains

  subroutine actual_integrator_init()

  end subroutine actual_integrator_init


  subroutine actual_integrator(state_in, state_out, dt, time)

    use vode_rpar_indices
    use vode_rhs_module
    use cuvode_module, only: dvode
    use cuvode_types_module, only: dvode_t, rwork_t
    use extern_probin_module, only: jacobian, burner_verbose, &
                                    rtol_spec, rtol_temp, rtol_enuc, &
                                    atol_spec, atol_temp, atol_enuc, &
                                    burning_mode, retry_burn, &
                                    retry_burn_factor, retry_burn_max_change, &
                                    call_eos_in_rhs, dT_crit, use_jacobian_caching
    use cuvode_parameters_module

    ! Input arguments

    type (sdc_t), intent(in   ) :: state_in
    type (sdc_t), intent(inout) :: state_out
    real(rt),    intent(in   ) :: dt, time

    ! Local variables

    real(rt) :: local_time

    ! Work arrays

    type(rwork_t) :: rwork
    integer    :: iwork(VODE_LIW)

    integer :: MF_JAC

    ! istate determines the state of the calculation.  A value of 1 meeans
    ! this is the first call to the problem -- this is what we will want.

    integer :: istate

    integer :: ipar

    real(rt) :: sum
    real(rt) :: retry_change_factor
    type (dvode_t) :: dvode_state

    !$gpu

    if (jacobian == 1) then ! Analytical
       MF_JAC = MF_ANALYTIC_JAC_CACHED
    else if (jacobian == 2) then ! Numerical
       MF_JAC = MF_NUMERICAL_JAC_CACHED
    else
#ifndef AMREX_USE_CUDA
       call amrex_error("Error: unknown Jacobian mode in actual_integrator.f90.")
#endif
    endif

    if (.not. use_jacobian_caching) then
       MF_JAC = -MF_JAC
    endif

    ! Set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.
    !
    ! **NOTE** if you reduce these tolerances, you probably will need
    ! to (a) decrease dT_crit, (b) increase the maximum number of
    ! steps allowed.

#if defined(SDC_EVOLVE_ENERGY)

    dvode_state % atol(SFS:SFS-1+nspec) = atol_spec ! mass fractions
    dvode_state % atol(SEDEN)           = atol_enuc ! temperature
    dvode_state % atol(SEINT)           = atol_enuc ! energy generated

    dvode_state % rtol(SFS:SFS-1+nspec) = rtol_spec ! mass fractions
    dvode_state % rtol(SEDEN)           = rtol_enuc ! temperature
    dvode_state % rtol(SEINT)           = rtol_enuc ! energy generated

#elif defined(SDC_EVOLVE_ENTHALPY)

    atol(SFS:SFS-1+nspec) = status % atol_spec ! mass fractions
    atol(SENTH)           = status % atol_enuc ! enthalpy

    rtol(SFS:SFS-1+nspec) = status % rtol_spec ! mass fractions
    rtol(SENTH)           = status % rtol_enuc ! enthalpy

#endif

    ! We want VODE to re-initialize each time we call it.

    dvode_state % istate = 1

    ! Initialize work arrays to zero.
    rwork % CONDOPT = ZERO
    rwork % YH   = ZERO
    rwork % WM   = ZERO
    rwork % EWT  = ZERO
    rwork % SAVF = ZERO
    rwork % ACOR = ZERO    
    iwork(:) = 0

    ! Set the maximum number of steps allowed (the VODE default is 500).

    iwork(6) = 150000

    ! Disable printing of messages about T + H == T unless we are in verbose mode.

    if (burner_verbose) then
       iwork(7) = 1
    else
       iwork(7) = 0
    endif

    ! Initialize the integration time.

    dvode_state % T = ZERO
    dvode_state % TOUT = dt

    ! Convert our input sdc state into the form VODE expects

    call sdc_to_vode(state_in, dvode_state % y, dvode_state % rpar)


    ! this is not used but we set it to prevent accessing uninitialzed
    ! data in common routines with the non-SDC integrator
    dvode_state % rpar(irp_self_heat) = -ONE

    ! Set the time offset -- this converts between the local integration
    ! time and the simulation time
    dvode_state % rpar(irp_t0) = time


    ! Call the integration routine.
    call dvode(dvode_state, rwork, iwork, ITASK, IOPT, MF_JAC)


    ! Store the final data
    call vode_to_sdc(time, dvode_state % y, dvode_state % rpar, state_out)

    ! get the number of RHS calls and jac evaluations from the VODE
    ! work arrays
    state_out % n_rhs = iwork(12)
    state_out % n_jac = iwork(13)

  end subroutine actual_integrator

end module actual_integrator_module
