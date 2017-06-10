! Common variables and routines for burners
! that use VODE for their integration.

module actual_integrator_module

  use eos_module, only: eos, eos_input_rt
  use network
  use rpar_indices
  use vode_type_module
  use burn_type_module
  use bl_types

  implicit none
  
contains

  subroutine actual_integrator_init()

    implicit none

  end subroutine actual_integrator_init


  ! Main interface
#ifdef CUDA
  attributes(device) &
#endif
  subroutine actual_integrator(state_in, state_out, dt, time)

    !$acc routine seq
    
    use rpar_indices
    use managed_probin_module, only: cu_jacobian, cu_burner_verbose, &
         cu_rtol_spec, cu_rtol_temp, cu_rtol_enuc, &
         cu_atol_spec, cu_atol_temp, cu_atol_enuc, &
         cu_burning_mode, cu_retry_burn, &
         cu_retry_burn_factor, cu_retry_burn_max_change, &
         cu_call_eos_in_rhs, cu_dt_crit
    use vode_rhs_module, only: f_rhs, jac    
    use actual_rhs_module, only : update_unevolved_species
    use dvode_module, only: dvode
    use eos_type_module, only: eos_t, copy_eos_t
    use dvode_type_module, only: dvode_t
    use bl_constants_module, only: ZERO    

    implicit none

    ! Input arguments

    type (burn_t), intent(in   ) :: state_in
    type (burn_t), intent(inout) :: state_out
    real(dp_t),    intent(in   ) :: dt, time

    ! Local variables

    real(dp_t) :: local_time
    type (eos_t) :: eos_state_in, eos_state_out, eos_state_temp
    type (dvode_t) :: dvode_state

    ! Work arrays

    real(dp_t) :: y(neqs)
    real(dp_t) :: atol(neqs), rtol(neqs)
    integer    :: iwork(LIW)
    real(dp_t) :: rpar(n_rpar_comps)

    integer :: MF_JAC

    ! istate determines the state of the calculation.  A value of 1 meeans
    ! this is the first call to the problem -- this is what we will want.

    integer :: istate

    integer :: ipar(n_ipar_comps)

    real(dp_t) :: sum
    real(dp_t) :: retry_change_factor

    real(dp_t) :: ener_offset

    type(rwork_t) :: rwork

    if (cu_jacobian == 1) then ! Analytical
       MF_JAC = MF_ANALYTIC_JAC
    else if (cu_jacobian == 2) then ! Numerical
       MF_JAC = MF_NUMERICAL_JAC
    else
       stop
       !CUDA
       !call bl_error("Error: unknown Jacobian mode in actual_integrator.f90.")
    endif

    ! Set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.
    !
    ! **NOTE** if you reduce these tolerances, you probably will need
    ! to (a) decrease dT_crit, (b) increase the maximum number of
    ! steps allowed.

    atol(1:nspec_evolve) = cu_atol_spec ! mass fractions
    atol(net_itemp)      = cu_atol_temp ! temperature
    atol(net_ienuc)      = cu_atol_enuc ! energy generated

    rtol(1:nspec_evolve) = cu_rtol_spec ! mass fractions
    rtol(net_itemp)      = cu_rtol_temp ! temperature
    rtol(net_ienuc)      = cu_rtol_enuc ! energy generated

    ! We want VODE to re-initialize each time we call it.

    istate = 1

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

    if (cu_burner_verbose) then
       iwork(7) = 1
    else
       iwork(7) = 0
    endif

    ! Initialize the integration time.

    local_time = ZERO

    ! Convert our input burn state into an EOS type.

    call burn_to_eos(state_in, eos_state_in)

    ! We assume that (rho, T) coming in are valid, do an EOS call
    ! to fill the rest of the thermodynamic variables.

    call eos(eos_input_rt, eos_state_in)

    ! Convert the EOS state data into the form VODE expects.

    call eos_to_vode(eos_state_in, y, rpar)

    ener_offset = eos_state_in % e

    y(net_ienuc) = ener_offset

    ! Pass through whether we are doing self-heating.

    if (cu_burning_mode == 0 .or. cu_burning_mode == 2) then
       rpar(irp_self_heat) = -ONE
    else if (cu_burning_mode == 1 .or. cu_burning_mode == 3) then
       rpar(irp_self_heat) = ONE
    else
       stop
       !CUDA
       !call bl_error("Error: unknown burning_mode in actual_integrator.f90.")
    endif

    ! Copy in the zone size.

    rpar(irp_dx) = state_in % dx

    ! Set the sound crossing time.

    rpar(irp_t_sound) = state_in % dx / eos_state_in % cs

    ! Set the time offset -- this converts between the local integration 
    ! time and the simulation time

    rpar(irp_t0) = time

    ! If we are using the dT_crit functionality and therefore doing a linear
    ! interpolation of the specific heat in between EOS calls, do a second
    ! EOS call here to establish an initial slope.

    rpar(irp_Told) = eos_state_in % T

    if (cu_dT_crit < 1.0d19) then

       call copy_eos_t(eos_state_temp, eos_state_in)
       eos_state_temp % T = eos_state_in % T * (ONE + sqrt(epsilon(ONE)))

       call eos(eos_input_rt, eos_state_temp)

       rpar(irp_dcvdt) = (eos_state_temp % cv - eos_state_in % cv) / (eos_state_temp % T - eos_state_in % T)
       rpar(irp_dcpdt) = (eos_state_temp % cp - eos_state_in % cp) / (eos_state_temp % T - eos_state_in % T)

    endif

    ! Save the initial state.

    rpar(irp_y_init:irp_y_init + neqs - 1) = y

    ! Call the integration routine.

    call dvode(neqs, y, local_time, local_time + dt, &
               ITOL, rtol, atol, ITASK, &
               istate, IOPT, rwork, LRW, iwork, LIW, MF_JAC, &
               rpar, ipar, dvode_state)

    ! If we are using hybrid burning and the energy release was negative (or we failed),
    ! re-run this in self-heating mode.

    if ( cu_burning_mode == 2 .and. &
         (y(net_ienuc) - ener_offset < ZERO .or. &
          istate < 0) ) then

       rpar(irp_self_heat) = ONE

       istate = 1

       rwork % CONDOPT = ZERO
       rwork % YH   = ZERO
       rwork % WM   = ZERO
       rwork % EWT  = ZERO
       rwork % SAVF = ZERO
       rwork % ACOR = ZERO    
       iwork(:) = 0

       iwork(6) = 150000

       local_time = ZERO

       call eos_to_vode(eos_state_in, y, rpar)

       rpar(irp_Told) = eos_state_in % T

       if (cu_dT_crit < 1.0d19) then

          rpar(irp_dcvdt) = (eos_state_temp % cv - eos_state_in % cv) / (eos_state_temp % T - eos_state_in % T)
          rpar(irp_dcpdt) = (eos_state_temp % cp - eos_state_in % cp) / (eos_state_temp % T - eos_state_in % T)

       endif

       y(net_ienuc) = ener_offset

       call dvode(neqs, y, local_time, local_time + dt, &
                  ITOL, rtol, atol, ITASK, &
                  istate, IOPT, rwork, LRW, iwork, LIW, MF_JAC, &
                  rpar, ipar, dvode_state)

    endif

    ! If we still failed, print out the current state of the integration.

    if (istate < 0) then
       
#ifndef CUDA       
       print *, 'ERROR: integration failed in net'
       print *, 'istate = ', istate
       print *, 'time = ', local_time
       print *, 'dens = ', state_in % rho
       print *, 'temp start = ', state_in % T
       print *, 'xn start = ', state_in % xn
       print *, 'temp current = ', y(net_itemp)
       print *, 'xn current = ', y(1:nspec_evolve) * aion(1:nspec_evolve), &
            rpar(irp_nspec:irp_nspec+n_not_evolved-1) * aion(nspec_evolve+1:)
       print *, 'energy generated = ', y(net_ienuc) - ener_offset
#endif
       
       if (.not. cu_retry_burn) then

          stop
          !CUDA
          !call bl_error("ERROR in burner: integration failed")

       else

#ifndef CUDA          
          print *, 'Retrying burn with looser tolerances'
#endif          

          retry_change_factor = ONE

          do while (istate < 0 .and. retry_change_factor <= cu_retry_burn_max_change)

             retry_change_factor = retry_change_factor * cu_retry_burn_factor

             istate = 1

             rwork % CONDOPT = ZERO             
             rwork % YH   = ZERO
             rwork % WM   = ZERO
             rwork % EWT  = ZERO
             rwork % SAVF = ZERO
             rwork % ACOR = ZERO    
             iwork(:) = 0

             atol = atol * cu_retry_burn_factor
             rtol = rtol * cu_retry_burn_factor

             iwork(6) = 150000

             local_time = ZERO

             call eos_to_vode(eos_state_in, y, rpar)

             rpar(irp_Told) = eos_state_in % T

             if (cu_dT_crit < 1.0d19) then

                rpar(irp_dcvdt) = (eos_state_temp % cv - eos_state_in % cv) / (eos_state_temp % T - eos_state_in % T)
                rpar(irp_dcpdt) = (eos_state_temp % cp - eos_state_in % cp) / (eos_state_temp % T - eos_state_in % T)

             endif

             y(net_ienuc) = ener_offset

             call dvode(neqs, y, local_time, local_time + dt, &
                        ITOL, rtol, atol, ITASK, &
                        istate, IOPT, rwork, LRW, iwork, LIW, MF_JAC, &
                        rpar, ipar, dvode_state)

          enddo

          if (retry_change_factor > cu_retry_burn_max_change .and. istate < 0) then

             stop
             !CUDA
             !call bl_error("ERROR in burner: integration failed")

          endif

       endif

    endif

    ! Subtract the energy offset
    y(net_ienuc) = y(net_ienuc) - ener_offset

    ! Store the final data, and then normalize abundances.
    call vode_to_burn(y, rpar, state_out)

    ! get the number of RHS calls and jac evaluations from the VODE
    ! work arrays
    state_out % n_rhs = iwork(12)
    state_out % n_jac = iwork(13)

    if (nspec_evolve < nspec) then
       call update_unevolved_species(state_out)
    endif

    call normalize_abundances_burn(state_out)

#ifndef CUDA    
    if (cu_burner_verbose) then

       ! Print out some integration statistics, if desired.
       print *, 'integration summary: '
       print *, 'dens: ', state_out % rho, ' temp: ', state_out % T, &
                ' energy released: ', state_out % e - state_in % e
       print *, 'number of steps taken: ', iwork(11)
       print *, 'number of f evaluations: ', iwork(12)
    endif
#endif
    
  end subroutine actual_integrator

end module actual_integrator_module
