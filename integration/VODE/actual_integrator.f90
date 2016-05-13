! Common variables and routines for burners
! that use VODE for their integration.

module actual_integrator_module

  use eos_module
  use network
  use rpar_indices
  use vode_convert_module
  use burn_type_module

  implicit none

  ! Our problem is stiff, so tell ODEPACK that. 21 means stiff, jacobian
  ! function is supplied; 22 means stiff, figure out my jacobian through
  ! differencing.

  integer, parameter :: MF_ANALYTIC_JAC = 21, MF_NUMERICAL_JAC = 22

  ! Tolerance parameters:
  !
  !  itol specifies whether to use an single absolute tolerance for
  !  all variables (1), or to pass an array of absolute tolerances, one
  !  for each variable with a scalar relative tol (2), a scalar absolute
  !  and array of relative tolerances (3), or arrays for both (4).
  !
  !  The error is determined as e(i) = rtol*abs(y(i)) + atol, and must
  !  be > 0.  Since we have some compositions that may be 0 initially,
  !  we will specify both an absolute and a relative tolerance.
  !
  ! We will use arrays for both the absolute and relative tolerances,
  ! since we want to be easier on the temperature than the species.

  integer, parameter :: ITOL = 4

  ! We want to do a normal computation, and get the output values of y(t)
  ! after stepping though dt.

  integer, PARAMETER :: ITASK = 1

  ! We will override the maximum number of steps, so turn on the
  ! optional arguments flag.

  integer, parameter :: IOPT = 1

  ! Declare a real work array of size 22 + 9*NEQ + 2*NEQ**2 and an
  ! integer work array of size 30 + NEQ. These are VODE constants
  ! that depend on the integration mode we're using -- see dvode.f.

  integer, parameter :: LRW = 22 + 9*neqs + 2*neqs**2
  integer, parameter :: LIW = 30 + neqs

contains

  subroutine actual_integrator_init()

    implicit none

    call init_rpar_indices()

  end subroutine actual_integrator_init



  ! Main interface

  subroutine actual_integrator(state_in, state_out, dt, time)

    use rpar_indices
    use extern_probin_module, only: jacobian, burner_verbose, &
                                    rtol_spec, rtol_temp, rtol_enuc, &
                                    atol_spec, atol_temp, atol_enuc, &
                                    burning_mode, retry_burn, &
                                    retry_burn_factor, retry_burn_max_change, &
                                    call_eos_in_rhs, dT_crit
    use integration_data, only: ener_scale

    implicit none

    ! Input arguments

    type (burn_t),       intent(in   ) :: state_in
    type (burn_t),       intent(inout) :: state_out
    double precision,    intent(in   ) :: dt, time

    ! Local variables

    double precision :: local_time
    type (eos_t)     :: eos_state_in, eos_state_out, eos_state_temp

    ! Work arrays

    double precision :: y(neqs)
    double precision :: atol(neqs), rtol(neqs)
    double precision :: rwork(LRW)
    integer          :: iwork(LIW)
    double precision :: rpar(n_rpar_comps)

    integer :: MF_JAC

    ! istate determines the state of the calculation.  A value of 1 meeans
    ! this is the first call to the problem -- this is what we will want.

    integer :: istate

    integer :: ipar

    double precision :: sum
    double precision :: retry_change_factor

    EXTERNAL jac, f_rhs

    if (jacobian == 1) then ! Analytical
       MF_JAC = MF_ANALYTIC_JAC
    else if (jacobian == 2) then ! Numerical
       MF_JAC = MF_NUMERICAL_JAC
    else
       call bl_error("Error: unknown Jacobian mode in actual_integrator.f90.")
    endif

    ! Set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.
    !
    ! **NOTE** if you reduce these tolerances, you probably will need
    ! to (a) decrease dT_crit, (b) increase the maximum number of
    ! steps allowed.

    atol(1:nspec)   = atol_spec ! mass fractions
    atol(net_itemp) = atol_temp ! temperature
    atol(net_ienuc) = atol_enuc ! energy generated

    rtol(1:nspec)   = rtol_spec ! mass fractions
    rtol(net_itemp) = rtol_temp ! temperature
    rtol(net_ienuc) = rtol_enuc ! energy generated

    ! We want VODE to re-initialize each time we call it.

    istate = 1

    ! Initialize work arrays to zero.

    rwork(:) = ZERO
    iwork(:) = 0

    ! Set the maximum number of steps allowed (the VODE default is 500).

    iwork(6) = 150000

    ! Initialize the integration time.

    local_time = ZERO

    ! Convert our input burn state into an EOS type.

    call burn_to_eos(state_in, eos_state_in)

    ! We assume that the valid quantities coming in are (rho, e); do an EOS call
    ! to make sure all other variables are consistent.

    call eos(eos_input_burn, eos_state_in)

    ! Send this data back to the burn state in case the energy changed
    ! due to a reset/flooring.

    call eos_to_burn(eos_state_in, state_in)

    ! Convert the EOS state data into the form VODE expects.

    call eos_to_vode(eos_state_in, y, rpar)

    y(net_ienuc) = ZERO

    ! Pass through whether we are doing self-heating.

    if (burning_mode == 0 .or. burning_mode == 2) then
       rpar(irp_self_heat) = -ONE
    else if (burning_mode == 1) then
       rpar(irp_self_heat) = ONE
    else
       call bl_error("Error: unknown burning_mode in actual_integrator.f90.")
    endif

    ! If we are using the dT_crit functionality and therefore doing a linear
    ! interpolation of the specific heat in between EOS calls, do a second
    ! EOS call here to establish an initial slope.

    rpar(irp_Told) = eos_state_in % T

    if (dT_crit < 1.0d19) then

       eos_state_temp = eos_state_in
       eos_state_temp % T = eos_state_in % T * (ONE + sqrt(epsilon(ONE)))

       call eos(eos_input_rt, eos_state_temp)

       rpar(irp_dcvdt) = (eos_state_temp % cv - eos_state_in % cv) / (eos_state_temp % T - eos_state_in % T)
       rpar(irp_dcpdt) = (eos_state_temp % cp - eos_state_in % cp) / (eos_state_temp % T - eos_state_in % T)

    endif

    ! Call the integration routine.

    call dvode(f_rhs, neqs, y, local_time, dt, ITOL, rtol, atol, ITASK, &
               istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_JAC, rpar, ipar)

    ! If we are using hybrid burning and the energy release was negative (or we failed),
    ! re-run this in self-heating mode.

    if ( burning_mode == 2 .and. (y(net_ienuc) < ZERO .or. istate < 0) ) then

       rpar(irp_self_heat) = ONE

       istate = 1

       rwork(:) = ZERO
       iwork(:) = 0

       iwork(6) = 150000

       local_time = ZERO

       call eos_to_vode(eos_state_in, y, rpar)

       rpar(irp_Told) = eos_state_in % T

       if (dT_crit < 1.0d19) then

          rpar(irp_dcvdt) = (eos_state_temp % cv - eos_state_in % cv) / (eos_state_temp % T - eos_state_in % T)
          rpar(irp_dcpdt) = (eos_state_temp % cp - eos_state_in % cp) / (eos_state_temp % T - eos_state_in % T)

       endif

       y(net_ienuc) = ZERO

       call dvode(f_rhs, neqs, y, local_time, dt, ITOL, rtol, atol, ITASK, &
                  istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_JAC, rpar, ipar)

    endif

    ! If we still failed, print out the current state of the integration.

    if (istate < 0) then
       print *, 'ERROR: integration failed in net'
       print *, 'istate = ', istate
       print *, 'time = ', local_time
       print *, 'dens = ', state_in % rho
       print *, 'temp start = ', state_in % T
       print *, 'xn start = ', state_in % xn
       print *, 'temp current = ', y(net_itemp) * temp_scale
       print *, 'xn current = ', y(1:nspec)
       print *, 'energy generated = ', y(net_ienuc)

       if (.not. retry_burn) then

          call bl_error("ERROR in burner: integration failed")

       else

          print *, 'Retrying burn with looser tolerances'

          retry_change_factor = ONE

          do while (istate < 0 .and. retry_change_factor <= retry_burn_max_change)

             retry_change_factor = retry_change_factor * retry_burn_factor

             istate = 1

             rwork(:) = ZERO
             iwork(:) = 0

             atol = atol * retry_burn_factor
             rtol = rtol * retry_burn_factor

             iwork(6) = 150000

             local_time = ZERO

             call eos_to_vode(eos_state_in, y, rpar)

             rpar(irp_Told) = eos_state_in % T

             if (dT_crit < 1.0d19) then

                rpar(irp_dcvdt) = (eos_state_temp % cv - eos_state_in % cv) / (eos_state_temp % T - eos_state_in % T)
                rpar(irp_dcpdt) = (eos_state_temp % cp - eos_state_in % cp) / (eos_state_temp % T - eos_state_in % T)

             endif

             y(net_ienuc) = ZERO

             call dvode(f_rhs, neqs, y, local_time, dt, ITOL, rtol, atol, ITASK, &
                        istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_JAC, rpar, ipar)

          enddo

          if (retry_change_factor > retry_burn_max_change .and. istate < 0) then

             call bl_error("ERROR in burner: integration failed")

          endif

       endif

    endif

    ! Store the final data.

    call vode_to_eos(eos_state_out, y, rpar)

    call normalize_abundances(eos_state_out)

    ! Energy was integrated in the system -- we use this integrated
    ! energy which contains both the reaction energy release and
    ! neutrino losses. The final energy is the initial energy
    ! plus this energy release. Note that we get a new temperature too,
    ! but we will discard it and call the EOS to get a final temperature
    ! consistent with this new energy.

    eos_state_out % e = eos_state_in % e + y(net_ienuc) * ener_scale

    eos_state_out % reset = .true.

    call eos(eos_input_re, eos_state_out)

    call eos_to_burn(eos_state_out, state_out)

    if (burner_verbose) then

       ! Print out some integration statistics, if desired.

       print *, 'integration summary: '
       print *, 'dens: ', state_out % rho, ' temp: ', state_out % T, &
                ' energy released: ', state_out % e - state_in % e
       print *, 'number of steps taken: ', iwork(11)
       print *, 'number of f evaluations: ', iwork(12)

    endif

  end subroutine actual_integrator

end module actual_integrator_module
