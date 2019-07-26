module actual_rhs_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module
  use physical_constants, only: N_AVO
  use network
  use table_rates
  use burn_type_module

  implicit none

  type :: rate_eval_t
     real(rt) :: unscreened_rates(4, nrates)
     real(rt) :: screened_rates(nrates)
     real(rt) :: add_energy(nrat_tabular)
     real(rt) :: add_energy_rate(nrat_tabular)
  end type rate_eval_t
  
contains

  subroutine actual_rhs_init()
    ! STUB FOR MAESTRO'S TEST_REACT. ALL THE INIT IS DONE BY BURNER_INIT
    return
  end subroutine actual_rhs_init


  subroutine update_unevolved_species(state)
    ! STUB FOR INTEGRATOR
    type(burn_t)     :: state

    !$gpu
    
    return
  end subroutine update_unevolved_species


  subroutine zero_rate_eval(rate_eval)

    implicit none

    type(rate_eval_t), intent(inout) :: rate_eval

    !$gpu

    rate_eval % unscreened_rates(i_rate, :) = ZERO
    rate_eval % unscreened_rates(i_drate_dt, :) = ZERO
    rate_eval % unscreened_rates(i_scor, :) = ONE
    rate_eval % unscreened_rates(i_dscor_dt, :) = ZERO
    rate_eval % screened_rates = ZERO
    rate_eval % add_energy = ZERO
    rate_eval % add_energy_rate = ZERO

  end subroutine zero_rate_eval


  subroutine evaluate_rates(state, rate_eval)
    !$acc routine seq

    use reaclib_rates, only: screen_reaclib, reaclib_evaluate
    use screening_module, only: screen5, plasma_state, fill_plasma_state

    implicit none
    
    type(burn_t)     :: state
    type(rate_eval_t), intent(out) :: rate_eval
    type(plasma_state) :: pstate
    real(rt) :: Y(nspec)
    real(rt) :: raw_rates(4, nrates)
    real(rt) :: reactvec(num_rate_groups+2)
    integer :: i, j
    real(rt) :: rhoy, scor, dscor_dt, dscor_dd

    !$gpu

    Y(:) = state % xn(:) * aion_inv(:)
    rhoy = state % rho * state % y_e

    ! Zero out the rates
    call zero_rate_eval(rate_eval)

    ! Calculate Reaclib rates
    call fill_plasma_state(pstate, state % T, state % rho, Y)
    do i = 1, nrat_reaclib
       call reaclib_evaluate(pstate, state % T, i, reactvec)
       rate_eval % unscreened_rates(:,i) = reactvec(1:4)
    end do

    ! Evaluate screening factors
    if (screen_reaclib) then

      call screen5(pstate, 1, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,34) = scor
      rate_eval % unscreened_rates(i_dscor_dt,34) = dscor_dt
      rate_eval % unscreened_rates(i_scor,35) = scor
      rate_eval % unscreened_rates(i_dscor_dt,35) = dscor_dt
      rate_eval % unscreened_rates(i_scor,89) = scor
      rate_eval % unscreened_rates(i_dscor_dt,89) = dscor_dt
      rate_eval % unscreened_rates(i_scor,91) = scor
      rate_eval % unscreened_rates(i_dscor_dt,91) = dscor_dt
      rate_eval % unscreened_rates(i_scor,528) = scor
      rate_eval % unscreened_rates(i_dscor_dt,528) = dscor_dt


      call screen5(pstate, 2, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,36) = scor
      rate_eval % unscreened_rates(i_dscor_dt,36) = dscor_dt


      call screen5(pstate, 3, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,37) = scor
      rate_eval % unscreened_rates(i_dscor_dt,37) = dscor_dt


      call screen5(pstate, 4, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,38) = scor
      rate_eval % unscreened_rates(i_dscor_dt,38) = dscor_dt


      call screen5(pstate, 5, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,39) = scor
      rate_eval % unscreened_rates(i_dscor_dt,39) = dscor_dt


      call screen5(pstate, 6, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,40) = scor
      rate_eval % unscreened_rates(i_dscor_dt,40) = dscor_dt


      call screen5(pstate, 7, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,41) = scor
      rate_eval % unscreened_rates(i_dscor_dt,41) = dscor_dt


      call screen5(pstate, 8, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,42) = scor
      rate_eval % unscreened_rates(i_dscor_dt,42) = dscor_dt
      rate_eval % unscreened_rates(i_scor,63) = scor
      rate_eval % unscreened_rates(i_dscor_dt,63) = dscor_dt


      call screen5(pstate, 9, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,43) = scor
      rate_eval % unscreened_rates(i_dscor_dt,43) = dscor_dt


      call screen5(pstate, 10, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,44) = scor
      rate_eval % unscreened_rates(i_dscor_dt,44) = dscor_dt


      call screen5(pstate, 11, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,45) = scor
      rate_eval % unscreened_rates(i_dscor_dt,45) = dscor_dt
      rate_eval % unscreened_rates(i_scor,66) = scor
      rate_eval % unscreened_rates(i_dscor_dt,66) = dscor_dt


      call screen5(pstate, 12, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,46) = scor
      rate_eval % unscreened_rates(i_dscor_dt,46) = dscor_dt
      rate_eval % unscreened_rates(i_scor,67) = scor
      rate_eval % unscreened_rates(i_dscor_dt,67) = dscor_dt


      call screen5(pstate, 13, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,47) = scor
      rate_eval % unscreened_rates(i_dscor_dt,47) = dscor_dt
      rate_eval % unscreened_rates(i_scor,68) = scor
      rate_eval % unscreened_rates(i_dscor_dt,68) = dscor_dt


      call screen5(pstate, 14, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,48) = scor
      rate_eval % unscreened_rates(i_dscor_dt,48) = dscor_dt
      rate_eval % unscreened_rates(i_scor,69) = scor
      rate_eval % unscreened_rates(i_dscor_dt,69) = dscor_dt


      call screen5(pstate, 15, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,49) = scor
      rate_eval % unscreened_rates(i_dscor_dt,49) = dscor_dt
      rate_eval % unscreened_rates(i_scor,70) = scor
      rate_eval % unscreened_rates(i_dscor_dt,70) = dscor_dt


      call screen5(pstate, 16, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,50) = scor
      rate_eval % unscreened_rates(i_dscor_dt,50) = dscor_dt
      rate_eval % unscreened_rates(i_scor,71) = scor
      rate_eval % unscreened_rates(i_dscor_dt,71) = dscor_dt


      call screen5(pstate, 17, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,51) = scor
      rate_eval % unscreened_rates(i_dscor_dt,51) = dscor_dt
      rate_eval % unscreened_rates(i_scor,72) = scor
      rate_eval % unscreened_rates(i_dscor_dt,72) = dscor_dt


      call screen5(pstate, 18, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,52) = scor
      rate_eval % unscreened_rates(i_dscor_dt,52) = dscor_dt
      rate_eval % unscreened_rates(i_scor,73) = scor
      rate_eval % unscreened_rates(i_dscor_dt,73) = dscor_dt


      call screen5(pstate, 19, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,53) = scor
      rate_eval % unscreened_rates(i_dscor_dt,53) = dscor_dt
      rate_eval % unscreened_rates(i_scor,74) = scor
      rate_eval % unscreened_rates(i_dscor_dt,74) = dscor_dt


      call screen5(pstate, 20, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,54) = scor
      rate_eval % unscreened_rates(i_dscor_dt,54) = dscor_dt
      rate_eval % unscreened_rates(i_scor,75) = scor
      rate_eval % unscreened_rates(i_dscor_dt,75) = dscor_dt


      call screen5(pstate, 21, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,55) = scor
      rate_eval % unscreened_rates(i_dscor_dt,55) = dscor_dt
      rate_eval % unscreened_rates(i_scor,76) = scor
      rate_eval % unscreened_rates(i_dscor_dt,76) = dscor_dt


      call screen5(pstate, 22, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,56) = scor
      rate_eval % unscreened_rates(i_dscor_dt,56) = dscor_dt
      rate_eval % unscreened_rates(i_scor,78) = scor
      rate_eval % unscreened_rates(i_dscor_dt,78) = dscor_dt


      call screen5(pstate, 23, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,57) = scor
      rate_eval % unscreened_rates(i_dscor_dt,57) = dscor_dt
      rate_eval % unscreened_rates(i_scor,80) = scor
      rate_eval % unscreened_rates(i_dscor_dt,80) = dscor_dt


      call screen5(pstate, 24, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,58) = scor
      rate_eval % unscreened_rates(i_dscor_dt,58) = dscor_dt
      rate_eval % unscreened_rates(i_scor,81) = scor
      rate_eval % unscreened_rates(i_dscor_dt,81) = dscor_dt


      call screen5(pstate, 25, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,59) = scor
      rate_eval % unscreened_rates(i_dscor_dt,59) = dscor_dt


      call screen5(pstate, 26, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,60) = scor
      rate_eval % unscreened_rates(i_dscor_dt,60) = dscor_dt
      rate_eval % unscreened_rates(i_scor,90) = scor
      rate_eval % unscreened_rates(i_dscor_dt,90) = dscor_dt


      call screen5(pstate, 27, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,61) = scor
      rate_eval % unscreened_rates(i_dscor_dt,61) = dscor_dt
      rate_eval % unscreened_rates(i_scor,88) = scor
      rate_eval % unscreened_rates(i_dscor_dt,88) = dscor_dt


      call screen5(pstate, 28, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,62) = scor
      rate_eval % unscreened_rates(i_dscor_dt,62) = dscor_dt


      call screen5(pstate, 29, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,64) = scor
      rate_eval % unscreened_rates(i_dscor_dt,64) = dscor_dt


      call screen5(pstate, 30, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,65) = scor
      rate_eval % unscreened_rates(i_dscor_dt,65) = dscor_dt


      call screen5(pstate, 31, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,77) = scor
      rate_eval % unscreened_rates(i_dscor_dt,77) = dscor_dt
      rate_eval % unscreened_rates(i_scor,97) = scor
      rate_eval % unscreened_rates(i_dscor_dt,97) = dscor_dt


      call screen5(pstate, 32, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,79) = scor
      rate_eval % unscreened_rates(i_dscor_dt,79) = dscor_dt
      rate_eval % unscreened_rates(i_scor,98) = scor
      rate_eval % unscreened_rates(i_dscor_dt,98) = dscor_dt


      call screen5(pstate, 33, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,82) = scor
      rate_eval % unscreened_rates(i_dscor_dt,82) = dscor_dt
      rate_eval % unscreened_rates(i_scor,109) = scor
      rate_eval % unscreened_rates(i_dscor_dt,109) = dscor_dt


      call screen5(pstate, 34, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,83) = scor
      rate_eval % unscreened_rates(i_dscor_dt,83) = dscor_dt
      rate_eval % unscreened_rates(i_scor,110) = scor
      rate_eval % unscreened_rates(i_dscor_dt,110) = dscor_dt
      rate_eval % unscreened_rates(i_scor,111) = scor
      rate_eval % unscreened_rates(i_dscor_dt,111) = dscor_dt


      call screen5(pstate, 35, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,84) = scor
      rate_eval % unscreened_rates(i_dscor_dt,84) = dscor_dt
      rate_eval % unscreened_rates(i_scor,112) = scor
      rate_eval % unscreened_rates(i_dscor_dt,112) = dscor_dt


      call screen5(pstate, 36, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,85) = scor
      rate_eval % unscreened_rates(i_dscor_dt,85) = dscor_dt


      call screen5(pstate, 37, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,86) = scor
      rate_eval % unscreened_rates(i_dscor_dt,86) = dscor_dt


      call screen5(pstate, 38, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,87) = scor
      rate_eval % unscreened_rates(i_dscor_dt,87) = dscor_dt


      call screen5(pstate, 39, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,94) = scor
      rate_eval % unscreened_rates(i_dscor_dt,94) = dscor_dt
      rate_eval % unscreened_rates(i_scor,95) = scor
      rate_eval % unscreened_rates(i_dscor_dt,95) = dscor_dt


      call screen5(pstate, 40, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,96) = scor
      rate_eval % unscreened_rates(i_dscor_dt,96) = dscor_dt


      call screen5(pstate, 41, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,99) = scor
      rate_eval % unscreened_rates(i_dscor_dt,99) = dscor_dt
      rate_eval % unscreened_rates(i_scor,100) = scor
      rate_eval % unscreened_rates(i_dscor_dt,100) = dscor_dt


      call screen5(pstate, 42, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,101) = scor
      rate_eval % unscreened_rates(i_dscor_dt,101) = dscor_dt


      call screen5(pstate, 43, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,102) = scor
      rate_eval % unscreened_rates(i_dscor_dt,102) = dscor_dt


      call screen5(pstate, 44, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,103) = scor
      rate_eval % unscreened_rates(i_dscor_dt,103) = dscor_dt
      rate_eval % unscreened_rates(i_scor,104) = scor
      rate_eval % unscreened_rates(i_dscor_dt,104) = dscor_dt


      call screen5(pstate, 45, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,105) = scor
      rate_eval % unscreened_rates(i_dscor_dt,105) = dscor_dt
      rate_eval % unscreened_rates(i_scor,107) = scor
      rate_eval % unscreened_rates(i_dscor_dt,107) = dscor_dt


      call screen5(pstate, 46, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,106) = scor
      rate_eval % unscreened_rates(i_dscor_dt,106) = dscor_dt
      rate_eval % unscreened_rates(i_scor,108) = scor
      rate_eval % unscreened_rates(i_dscor_dt,108) = dscor_dt


      call screen5(pstate, 47, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,113) = scor
      rate_eval % unscreened_rates(i_dscor_dt,113) = dscor_dt
      rate_eval % unscreened_rates(i_scor,114) = scor
      rate_eval % unscreened_rates(i_dscor_dt,114) = dscor_dt


      call screen5(pstate, 48, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,115) = scor
      rate_eval % unscreened_rates(i_dscor_dt,115) = dscor_dt


      call screen5(pstate, 49, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,117) = scor
      rate_eval % unscreened_rates(i_dscor_dt,117) = dscor_dt


      call screen5(pstate, 50, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,119) = scor
      rate_eval % unscreened_rates(i_dscor_dt,119) = dscor_dt
      rate_eval % unscreened_rates(i_scor,120) = scor
      rate_eval % unscreened_rates(i_dscor_dt,120) = dscor_dt


      call screen5(pstate, 51, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,121) = scor
      rate_eval % unscreened_rates(i_dscor_dt,121) = dscor_dt
      rate_eval % unscreened_rates(i_scor,123) = scor
      rate_eval % unscreened_rates(i_dscor_dt,123) = dscor_dt


      call screen5(pstate, 52, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,122) = scor
      rate_eval % unscreened_rates(i_dscor_dt,122) = dscor_dt


      call screen5(pstate, 53, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,127) = scor
      rate_eval % unscreened_rates(i_dscor_dt,127) = dscor_dt
      rate_eval % unscreened_rates(i_scor,129) = scor
      rate_eval % unscreened_rates(i_dscor_dt,129) = dscor_dt


      call screen5(pstate, 54, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,128) = scor
      rate_eval % unscreened_rates(i_dscor_dt,128) = dscor_dt
      rate_eval % unscreened_rates(i_scor,130) = scor
      rate_eval % unscreened_rates(i_dscor_dt,130) = dscor_dt


      call screen5(pstate, 55, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,134) = scor
      rate_eval % unscreened_rates(i_dscor_dt,134) = dscor_dt
      rate_eval % unscreened_rates(i_scor,136) = scor
      rate_eval % unscreened_rates(i_dscor_dt,136) = dscor_dt


      call screen5(pstate, 56, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,135) = scor
      rate_eval % unscreened_rates(i_dscor_dt,135) = dscor_dt
      rate_eval % unscreened_rates(i_scor,137) = scor
      rate_eval % unscreened_rates(i_dscor_dt,137) = dscor_dt


      call screen5(pstate, 57, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,140) = scor
      rate_eval % unscreened_rates(i_dscor_dt,140) = dscor_dt
      rate_eval % unscreened_rates(i_scor,142) = scor
      rate_eval % unscreened_rates(i_dscor_dt,142) = dscor_dt


      call screen5(pstate, 58, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,141) = scor
      rate_eval % unscreened_rates(i_dscor_dt,141) = dscor_dt
      rate_eval % unscreened_rates(i_scor,143) = scor
      rate_eval % unscreened_rates(i_dscor_dt,143) = dscor_dt


      call screen5(pstate, 59, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,144) = scor
      rate_eval % unscreened_rates(i_dscor_dt,144) = dscor_dt
      rate_eval % unscreened_rates(i_scor,146) = scor
      rate_eval % unscreened_rates(i_dscor_dt,146) = dscor_dt


      call screen5(pstate, 60, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,145) = scor
      rate_eval % unscreened_rates(i_dscor_dt,145) = dscor_dt


      call screen5(pstate, 61, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,149) = scor
      rate_eval % unscreened_rates(i_dscor_dt,149) = dscor_dt
      rate_eval % unscreened_rates(i_scor,150) = scor
      rate_eval % unscreened_rates(i_dscor_dt,150) = dscor_dt


      call screen5(pstate, 62, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,154) = scor
      rate_eval % unscreened_rates(i_dscor_dt,154) = dscor_dt
      rate_eval % unscreened_rates(i_scor,156) = scor
      rate_eval % unscreened_rates(i_dscor_dt,156) = dscor_dt


      call screen5(pstate, 63, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,155) = scor
      rate_eval % unscreened_rates(i_dscor_dt,155) = dscor_dt
      rate_eval % unscreened_rates(i_scor,157) = scor
      rate_eval % unscreened_rates(i_dscor_dt,157) = dscor_dt


      call screen5(pstate, 64, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,162) = scor
      rate_eval % unscreened_rates(i_dscor_dt,162) = dscor_dt
      rate_eval % unscreened_rates(i_scor,164) = scor
      rate_eval % unscreened_rates(i_dscor_dt,164) = dscor_dt


      call screen5(pstate, 65, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,163) = scor
      rate_eval % unscreened_rates(i_dscor_dt,163) = dscor_dt
      rate_eval % unscreened_rates(i_scor,165) = scor
      rate_eval % unscreened_rates(i_dscor_dt,165) = dscor_dt


      call screen5(pstate, 66, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,169) = scor
      rate_eval % unscreened_rates(i_dscor_dt,169) = dscor_dt
      rate_eval % unscreened_rates(i_scor,171) = scor
      rate_eval % unscreened_rates(i_dscor_dt,171) = dscor_dt


      call screen5(pstate, 67, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,170) = scor
      rate_eval % unscreened_rates(i_dscor_dt,170) = dscor_dt
      rate_eval % unscreened_rates(i_scor,172) = scor
      rate_eval % unscreened_rates(i_dscor_dt,172) = dscor_dt


      call screen5(pstate, 68, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,175) = scor
      rate_eval % unscreened_rates(i_dscor_dt,175) = dscor_dt
      rate_eval % unscreened_rates(i_scor,177) = scor
      rate_eval % unscreened_rates(i_dscor_dt,177) = dscor_dt


      call screen5(pstate, 69, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,176) = scor
      rate_eval % unscreened_rates(i_dscor_dt,176) = dscor_dt
      rate_eval % unscreened_rates(i_scor,178) = scor
      rate_eval % unscreened_rates(i_dscor_dt,178) = dscor_dt


      call screen5(pstate, 70, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,181) = scor
      rate_eval % unscreened_rates(i_dscor_dt,181) = dscor_dt
      rate_eval % unscreened_rates(i_scor,183) = scor
      rate_eval % unscreened_rates(i_dscor_dt,183) = dscor_dt


      call screen5(pstate, 71, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,182) = scor
      rate_eval % unscreened_rates(i_dscor_dt,182) = dscor_dt
      rate_eval % unscreened_rates(i_scor,184) = scor
      rate_eval % unscreened_rates(i_dscor_dt,184) = dscor_dt


      call screen5(pstate, 72, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,186) = scor
      rate_eval % unscreened_rates(i_dscor_dt,186) = dscor_dt
      rate_eval % unscreened_rates(i_scor,188) = scor
      rate_eval % unscreened_rates(i_dscor_dt,188) = dscor_dt


      call screen5(pstate, 73, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,187) = scor
      rate_eval % unscreened_rates(i_dscor_dt,187) = dscor_dt


      call screen5(pstate, 74, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,192) = scor
      rate_eval % unscreened_rates(i_dscor_dt,192) = dscor_dt
      rate_eval % unscreened_rates(i_scor,194) = scor
      rate_eval % unscreened_rates(i_dscor_dt,194) = dscor_dt


      call screen5(pstate, 75, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,193) = scor
      rate_eval % unscreened_rates(i_dscor_dt,193) = dscor_dt
      rate_eval % unscreened_rates(i_scor,195) = scor
      rate_eval % unscreened_rates(i_dscor_dt,195) = dscor_dt


      call screen5(pstate, 76, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,199) = scor
      rate_eval % unscreened_rates(i_dscor_dt,199) = dscor_dt
      rate_eval % unscreened_rates(i_scor,201) = scor
      rate_eval % unscreened_rates(i_dscor_dt,201) = dscor_dt


      call screen5(pstate, 77, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,200) = scor
      rate_eval % unscreened_rates(i_dscor_dt,200) = dscor_dt
      rate_eval % unscreened_rates(i_scor,202) = scor
      rate_eval % unscreened_rates(i_dscor_dt,202) = dscor_dt


      call screen5(pstate, 78, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,205) = scor
      rate_eval % unscreened_rates(i_dscor_dt,205) = dscor_dt
      rate_eval % unscreened_rates(i_scor,207) = scor
      rate_eval % unscreened_rates(i_dscor_dt,207) = dscor_dt


      call screen5(pstate, 79, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,206) = scor
      rate_eval % unscreened_rates(i_dscor_dt,206) = dscor_dt
      rate_eval % unscreened_rates(i_scor,208) = scor
      rate_eval % unscreened_rates(i_dscor_dt,208) = dscor_dt


      call screen5(pstate, 80, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,210) = scor
      rate_eval % unscreened_rates(i_dscor_dt,210) = dscor_dt
      rate_eval % unscreened_rates(i_scor,212) = scor
      rate_eval % unscreened_rates(i_dscor_dt,212) = dscor_dt


      call screen5(pstate, 81, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,211) = scor
      rate_eval % unscreened_rates(i_dscor_dt,211) = dscor_dt


      call screen5(pstate, 82, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,217) = scor
      rate_eval % unscreened_rates(i_dscor_dt,217) = dscor_dt


      call screen5(pstate, 83, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,218) = scor
      rate_eval % unscreened_rates(i_dscor_dt,218) = dscor_dt
      rate_eval % unscreened_rates(i_scor,219) = scor
      rate_eval % unscreened_rates(i_dscor_dt,219) = dscor_dt


      call screen5(pstate, 84, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,223) = scor
      rate_eval % unscreened_rates(i_dscor_dt,223) = dscor_dt
      rate_eval % unscreened_rates(i_scor,225) = scor
      rate_eval % unscreened_rates(i_dscor_dt,225) = dscor_dt


      call screen5(pstate, 85, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,224) = scor
      rate_eval % unscreened_rates(i_dscor_dt,224) = dscor_dt
      rate_eval % unscreened_rates(i_scor,226) = scor
      rate_eval % unscreened_rates(i_dscor_dt,226) = dscor_dt


      call screen5(pstate, 86, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,229) = scor
      rate_eval % unscreened_rates(i_dscor_dt,229) = dscor_dt
      rate_eval % unscreened_rates(i_scor,230) = scor
      rate_eval % unscreened_rates(i_dscor_dt,230) = dscor_dt


      call screen5(pstate, 87, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,235) = scor
      rate_eval % unscreened_rates(i_dscor_dt,235) = dscor_dt


      call screen5(pstate, 88, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,236) = scor
      rate_eval % unscreened_rates(i_dscor_dt,236) = dscor_dt
      rate_eval % unscreened_rates(i_scor,237) = scor
      rate_eval % unscreened_rates(i_dscor_dt,237) = dscor_dt


      call screen5(pstate, 89, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,243) = scor
      rate_eval % unscreened_rates(i_dscor_dt,243) = dscor_dt
      rate_eval % unscreened_rates(i_scor,245) = scor
      rate_eval % unscreened_rates(i_dscor_dt,245) = dscor_dt


      call screen5(pstate, 90, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,244) = scor
      rate_eval % unscreened_rates(i_dscor_dt,244) = dscor_dt
      rate_eval % unscreened_rates(i_scor,246) = scor
      rate_eval % unscreened_rates(i_dscor_dt,246) = dscor_dt


      call screen5(pstate, 91, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,250) = scor
      rate_eval % unscreened_rates(i_dscor_dt,250) = dscor_dt
      rate_eval % unscreened_rates(i_scor,252) = scor
      rate_eval % unscreened_rates(i_dscor_dt,252) = dscor_dt


      call screen5(pstate, 92, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,251) = scor
      rate_eval % unscreened_rates(i_dscor_dt,251) = dscor_dt
      rate_eval % unscreened_rates(i_scor,253) = scor
      rate_eval % unscreened_rates(i_dscor_dt,253) = dscor_dt


      call screen5(pstate, 93, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,256) = scor
      rate_eval % unscreened_rates(i_dscor_dt,256) = dscor_dt
      rate_eval % unscreened_rates(i_scor,258) = scor
      rate_eval % unscreened_rates(i_dscor_dt,258) = dscor_dt


      call screen5(pstate, 94, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,257) = scor
      rate_eval % unscreened_rates(i_dscor_dt,257) = dscor_dt
      rate_eval % unscreened_rates(i_scor,259) = scor
      rate_eval % unscreened_rates(i_dscor_dt,259) = dscor_dt


      call screen5(pstate, 95, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,262) = scor
      rate_eval % unscreened_rates(i_dscor_dt,262) = dscor_dt
      rate_eval % unscreened_rates(i_scor,264) = scor
      rate_eval % unscreened_rates(i_dscor_dt,264) = dscor_dt


      call screen5(pstate, 96, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,263) = scor
      rate_eval % unscreened_rates(i_dscor_dt,263) = dscor_dt
      rate_eval % unscreened_rates(i_scor,265) = scor
      rate_eval % unscreened_rates(i_dscor_dt,265) = dscor_dt


      call screen5(pstate, 97, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,267) = scor
      rate_eval % unscreened_rates(i_dscor_dt,267) = dscor_dt
      rate_eval % unscreened_rates(i_scor,269) = scor
      rate_eval % unscreened_rates(i_dscor_dt,269) = dscor_dt


      call screen5(pstate, 98, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,268) = scor
      rate_eval % unscreened_rates(i_dscor_dt,268) = dscor_dt
      rate_eval % unscreened_rates(i_scor,270) = scor
      rate_eval % unscreened_rates(i_dscor_dt,270) = dscor_dt


      call screen5(pstate, 99, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,274) = scor
      rate_eval % unscreened_rates(i_dscor_dt,274) = dscor_dt
      rate_eval % unscreened_rates(i_scor,276) = scor
      rate_eval % unscreened_rates(i_dscor_dt,276) = dscor_dt


      call screen5(pstate, 100, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,275) = scor
      rate_eval % unscreened_rates(i_dscor_dt,275) = dscor_dt
      rate_eval % unscreened_rates(i_scor,277) = scor
      rate_eval % unscreened_rates(i_dscor_dt,277) = dscor_dt


      call screen5(pstate, 101, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,281) = scor
      rate_eval % unscreened_rates(i_dscor_dt,281) = dscor_dt
      rate_eval % unscreened_rates(i_scor,283) = scor
      rate_eval % unscreened_rates(i_dscor_dt,283) = dscor_dt


      call screen5(pstate, 102, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,282) = scor
      rate_eval % unscreened_rates(i_dscor_dt,282) = dscor_dt
      rate_eval % unscreened_rates(i_scor,284) = scor
      rate_eval % unscreened_rates(i_dscor_dt,284) = dscor_dt


      call screen5(pstate, 103, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,287) = scor
      rate_eval % unscreened_rates(i_dscor_dt,287) = dscor_dt
      rate_eval % unscreened_rates(i_scor,289) = scor
      rate_eval % unscreened_rates(i_dscor_dt,289) = dscor_dt


      call screen5(pstate, 104, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,288) = scor
      rate_eval % unscreened_rates(i_dscor_dt,288) = dscor_dt
      rate_eval % unscreened_rates(i_scor,290) = scor
      rate_eval % unscreened_rates(i_dscor_dt,290) = dscor_dt


      call screen5(pstate, 105, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,292) = scor
      rate_eval % unscreened_rates(i_dscor_dt,292) = dscor_dt
      rate_eval % unscreened_rates(i_scor,294) = scor
      rate_eval % unscreened_rates(i_dscor_dt,294) = dscor_dt


      call screen5(pstate, 106, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,293) = scor
      rate_eval % unscreened_rates(i_dscor_dt,293) = dscor_dt
      rate_eval % unscreened_rates(i_scor,295) = scor
      rate_eval % unscreened_rates(i_dscor_dt,295) = dscor_dt


      call screen5(pstate, 107, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,299) = scor
      rate_eval % unscreened_rates(i_dscor_dt,299) = dscor_dt
      rate_eval % unscreened_rates(i_scor,300) = scor
      rate_eval % unscreened_rates(i_dscor_dt,300) = dscor_dt


      call screen5(pstate, 108, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,305) = scor
      rate_eval % unscreened_rates(i_dscor_dt,305) = dscor_dt
      rate_eval % unscreened_rates(i_scor,307) = scor
      rate_eval % unscreened_rates(i_dscor_dt,307) = dscor_dt


      call screen5(pstate, 109, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,306) = scor
      rate_eval % unscreened_rates(i_dscor_dt,306) = dscor_dt
      rate_eval % unscreened_rates(i_scor,308) = scor
      rate_eval % unscreened_rates(i_dscor_dt,308) = dscor_dt


      call screen5(pstate, 110, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,312) = scor
      rate_eval % unscreened_rates(i_dscor_dt,312) = dscor_dt
      rate_eval % unscreened_rates(i_scor,314) = scor
      rate_eval % unscreened_rates(i_dscor_dt,314) = dscor_dt


      call screen5(pstate, 111, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,313) = scor
      rate_eval % unscreened_rates(i_dscor_dt,313) = dscor_dt
      rate_eval % unscreened_rates(i_scor,315) = scor
      rate_eval % unscreened_rates(i_dscor_dt,315) = dscor_dt


      call screen5(pstate, 112, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,321) = scor
      rate_eval % unscreened_rates(i_dscor_dt,321) = dscor_dt
      rate_eval % unscreened_rates(i_scor,322) = scor
      rate_eval % unscreened_rates(i_dscor_dt,322) = dscor_dt


      call screen5(pstate, 113, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,327) = scor
      rate_eval % unscreened_rates(i_dscor_dt,327) = dscor_dt
      rate_eval % unscreened_rates(i_scor,329) = scor
      rate_eval % unscreened_rates(i_dscor_dt,329) = dscor_dt


      call screen5(pstate, 114, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,328) = scor
      rate_eval % unscreened_rates(i_dscor_dt,328) = dscor_dt
      rate_eval % unscreened_rates(i_scor,330) = scor
      rate_eval % unscreened_rates(i_dscor_dt,330) = dscor_dt


      call screen5(pstate, 115, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,336) = scor
      rate_eval % unscreened_rates(i_dscor_dt,336) = dscor_dt
      rate_eval % unscreened_rates(i_scor,338) = scor
      rate_eval % unscreened_rates(i_dscor_dt,338) = dscor_dt


      call screen5(pstate, 116, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,337) = scor
      rate_eval % unscreened_rates(i_dscor_dt,337) = dscor_dt
      rate_eval % unscreened_rates(i_scor,339) = scor
      rate_eval % unscreened_rates(i_dscor_dt,339) = dscor_dt


      call screen5(pstate, 117, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,343) = scor
      rate_eval % unscreened_rates(i_dscor_dt,343) = dscor_dt
      rate_eval % unscreened_rates(i_scor,345) = scor
      rate_eval % unscreened_rates(i_dscor_dt,345) = dscor_dt


      call screen5(pstate, 118, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,344) = scor
      rate_eval % unscreened_rates(i_dscor_dt,344) = dscor_dt
      rate_eval % unscreened_rates(i_scor,346) = scor
      rate_eval % unscreened_rates(i_dscor_dt,346) = dscor_dt


      call screen5(pstate, 119, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,349) = scor
      rate_eval % unscreened_rates(i_dscor_dt,349) = dscor_dt
      rate_eval % unscreened_rates(i_scor,351) = scor
      rate_eval % unscreened_rates(i_dscor_dt,351) = dscor_dt


      call screen5(pstate, 120, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,350) = scor
      rate_eval % unscreened_rates(i_dscor_dt,350) = dscor_dt
      rate_eval % unscreened_rates(i_scor,352) = scor
      rate_eval % unscreened_rates(i_dscor_dt,352) = dscor_dt


      call screen5(pstate, 121, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,355) = scor
      rate_eval % unscreened_rates(i_dscor_dt,355) = dscor_dt
      rate_eval % unscreened_rates(i_scor,357) = scor
      rate_eval % unscreened_rates(i_dscor_dt,357) = dscor_dt


      call screen5(pstate, 122, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,356) = scor
      rate_eval % unscreened_rates(i_dscor_dt,356) = dscor_dt
      rate_eval % unscreened_rates(i_scor,358) = scor
      rate_eval % unscreened_rates(i_dscor_dt,358) = dscor_dt


      call screen5(pstate, 123, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,359) = scor
      rate_eval % unscreened_rates(i_dscor_dt,359) = dscor_dt
      rate_eval % unscreened_rates(i_scor,361) = scor
      rate_eval % unscreened_rates(i_dscor_dt,361) = dscor_dt


      call screen5(pstate, 124, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,360) = scor
      rate_eval % unscreened_rates(i_dscor_dt,360) = dscor_dt


      call screen5(pstate, 125, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,364) = scor
      rate_eval % unscreened_rates(i_dscor_dt,364) = dscor_dt
      rate_eval % unscreened_rates(i_scor,366) = scor
      rate_eval % unscreened_rates(i_dscor_dt,366) = dscor_dt


      call screen5(pstate, 126, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,365) = scor
      rate_eval % unscreened_rates(i_dscor_dt,365) = dscor_dt
      rate_eval % unscreened_rates(i_scor,367) = scor
      rate_eval % unscreened_rates(i_dscor_dt,367) = dscor_dt


      call screen5(pstate, 127, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,371) = scor
      rate_eval % unscreened_rates(i_dscor_dt,371) = dscor_dt
      rate_eval % unscreened_rates(i_scor,373) = scor
      rate_eval % unscreened_rates(i_dscor_dt,373) = dscor_dt


      call screen5(pstate, 128, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,372) = scor
      rate_eval % unscreened_rates(i_dscor_dt,372) = dscor_dt
      rate_eval % unscreened_rates(i_scor,374) = scor
      rate_eval % unscreened_rates(i_dscor_dt,374) = dscor_dt


      call screen5(pstate, 129, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,378) = scor
      rate_eval % unscreened_rates(i_dscor_dt,378) = dscor_dt
      rate_eval % unscreened_rates(i_scor,380) = scor
      rate_eval % unscreened_rates(i_dscor_dt,380) = dscor_dt


      call screen5(pstate, 130, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,379) = scor
      rate_eval % unscreened_rates(i_dscor_dt,379) = dscor_dt
      rate_eval % unscreened_rates(i_scor,381) = scor
      rate_eval % unscreened_rates(i_dscor_dt,381) = dscor_dt


      call screen5(pstate, 131, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,384) = scor
      rate_eval % unscreened_rates(i_dscor_dt,384) = dscor_dt
      rate_eval % unscreened_rates(i_scor,386) = scor
      rate_eval % unscreened_rates(i_dscor_dt,386) = dscor_dt


      call screen5(pstate, 132, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,385) = scor
      rate_eval % unscreened_rates(i_dscor_dt,385) = dscor_dt
      rate_eval % unscreened_rates(i_scor,387) = scor
      rate_eval % unscreened_rates(i_dscor_dt,387) = dscor_dt


      call screen5(pstate, 133, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,388) = scor
      rate_eval % unscreened_rates(i_dscor_dt,388) = dscor_dt
      rate_eval % unscreened_rates(i_scor,390) = scor
      rate_eval % unscreened_rates(i_dscor_dt,390) = dscor_dt


      call screen5(pstate, 134, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,389) = scor
      rate_eval % unscreened_rates(i_dscor_dt,389) = dscor_dt


      call screen5(pstate, 135, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,393) = scor
      rate_eval % unscreened_rates(i_dscor_dt,393) = dscor_dt
      rate_eval % unscreened_rates(i_scor,395) = scor
      rate_eval % unscreened_rates(i_dscor_dt,395) = dscor_dt


      call screen5(pstate, 136, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,394) = scor
      rate_eval % unscreened_rates(i_dscor_dt,394) = dscor_dt
      rate_eval % unscreened_rates(i_scor,396) = scor
      rate_eval % unscreened_rates(i_dscor_dt,396) = dscor_dt


      call screen5(pstate, 137, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,401) = scor
      rate_eval % unscreened_rates(i_dscor_dt,401) = dscor_dt


      call screen5(pstate, 138, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,402) = scor
      rate_eval % unscreened_rates(i_dscor_dt,402) = dscor_dt
      rate_eval % unscreened_rates(i_scor,403) = scor
      rate_eval % unscreened_rates(i_dscor_dt,403) = dscor_dt


      call screen5(pstate, 139, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,408) = scor
      rate_eval % unscreened_rates(i_dscor_dt,408) = dscor_dt
      rate_eval % unscreened_rates(i_scor,410) = scor
      rate_eval % unscreened_rates(i_dscor_dt,410) = dscor_dt


      call screen5(pstate, 140, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,409) = scor
      rate_eval % unscreened_rates(i_dscor_dt,409) = dscor_dt
      rate_eval % unscreened_rates(i_scor,411) = scor
      rate_eval % unscreened_rates(i_dscor_dt,411) = dscor_dt


      call screen5(pstate, 141, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,415) = scor
      rate_eval % unscreened_rates(i_dscor_dt,415) = dscor_dt
      rate_eval % unscreened_rates(i_scor,417) = scor
      rate_eval % unscreened_rates(i_dscor_dt,417) = dscor_dt


      call screen5(pstate, 142, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,416) = scor
      rate_eval % unscreened_rates(i_dscor_dt,416) = dscor_dt
      rate_eval % unscreened_rates(i_scor,418) = scor
      rate_eval % unscreened_rates(i_dscor_dt,418) = dscor_dt


      call screen5(pstate, 143, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,422) = scor
      rate_eval % unscreened_rates(i_dscor_dt,422) = dscor_dt
      rate_eval % unscreened_rates(i_scor,423) = scor
      rate_eval % unscreened_rates(i_dscor_dt,423) = dscor_dt


      call screen5(pstate, 144, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,428) = scor
      rate_eval % unscreened_rates(i_dscor_dt,428) = dscor_dt
      rate_eval % unscreened_rates(i_scor,430) = scor
      rate_eval % unscreened_rates(i_dscor_dt,430) = dscor_dt


      call screen5(pstate, 145, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,429) = scor
      rate_eval % unscreened_rates(i_dscor_dt,429) = dscor_dt
      rate_eval % unscreened_rates(i_scor,431) = scor
      rate_eval % unscreened_rates(i_dscor_dt,431) = dscor_dt


      call screen5(pstate, 146, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,437) = scor
      rate_eval % unscreened_rates(i_dscor_dt,437) = dscor_dt
      rate_eval % unscreened_rates(i_scor,439) = scor
      rate_eval % unscreened_rates(i_dscor_dt,439) = dscor_dt


      call screen5(pstate, 147, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,438) = scor
      rate_eval % unscreened_rates(i_dscor_dt,438) = dscor_dt
      rate_eval % unscreened_rates(i_scor,440) = scor
      rate_eval % unscreened_rates(i_dscor_dt,440) = dscor_dt


      call screen5(pstate, 148, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,444) = scor
      rate_eval % unscreened_rates(i_dscor_dt,444) = dscor_dt
      rate_eval % unscreened_rates(i_scor,446) = scor
      rate_eval % unscreened_rates(i_dscor_dt,446) = dscor_dt


      call screen5(pstate, 149, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,445) = scor
      rate_eval % unscreened_rates(i_dscor_dt,445) = dscor_dt
      rate_eval % unscreened_rates(i_scor,447) = scor
      rate_eval % unscreened_rates(i_dscor_dt,447) = dscor_dt


      call screen5(pstate, 150, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,450) = scor
      rate_eval % unscreened_rates(i_dscor_dt,450) = dscor_dt
      rate_eval % unscreened_rates(i_scor,452) = scor
      rate_eval % unscreened_rates(i_dscor_dt,452) = dscor_dt


      call screen5(pstate, 151, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,451) = scor
      rate_eval % unscreened_rates(i_dscor_dt,451) = dscor_dt
      rate_eval % unscreened_rates(i_scor,453) = scor
      rate_eval % unscreened_rates(i_dscor_dt,453) = dscor_dt


      call screen5(pstate, 152, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,457) = scor
      rate_eval % unscreened_rates(i_dscor_dt,457) = dscor_dt
      rate_eval % unscreened_rates(i_scor,459) = scor
      rate_eval % unscreened_rates(i_dscor_dt,459) = dscor_dt


      call screen5(pstate, 153, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,458) = scor
      rate_eval % unscreened_rates(i_dscor_dt,458) = dscor_dt
      rate_eval % unscreened_rates(i_scor,460) = scor
      rate_eval % unscreened_rates(i_dscor_dt,460) = dscor_dt


      call screen5(pstate, 154, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,462) = scor
      rate_eval % unscreened_rates(i_dscor_dt,462) = dscor_dt
      rate_eval % unscreened_rates(i_scor,464) = scor
      rate_eval % unscreened_rates(i_dscor_dt,464) = dscor_dt


      call screen5(pstate, 155, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,463) = scor
      rate_eval % unscreened_rates(i_dscor_dt,463) = dscor_dt


      call screen5(pstate, 156, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,467) = scor
      rate_eval % unscreened_rates(i_dscor_dt,467) = dscor_dt
      rate_eval % unscreened_rates(i_scor,469) = scor
      rate_eval % unscreened_rates(i_dscor_dt,469) = dscor_dt


      call screen5(pstate, 157, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,468) = scor
      rate_eval % unscreened_rates(i_dscor_dt,468) = dscor_dt
      rate_eval % unscreened_rates(i_scor,470) = scor
      rate_eval % unscreened_rates(i_dscor_dt,470) = dscor_dt


      call screen5(pstate, 158, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,474) = scor
      rate_eval % unscreened_rates(i_dscor_dt,474) = dscor_dt
      rate_eval % unscreened_rates(i_scor,476) = scor
      rate_eval % unscreened_rates(i_dscor_dt,476) = dscor_dt


      call screen5(pstate, 159, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,475) = scor
      rate_eval % unscreened_rates(i_dscor_dt,475) = dscor_dt
      rate_eval % unscreened_rates(i_scor,477) = scor
      rate_eval % unscreened_rates(i_dscor_dt,477) = dscor_dt


      call screen5(pstate, 160, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,481) = scor
      rate_eval % unscreened_rates(i_dscor_dt,481) = dscor_dt
      rate_eval % unscreened_rates(i_scor,483) = scor
      rate_eval % unscreened_rates(i_dscor_dt,483) = dscor_dt


      call screen5(pstate, 161, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,482) = scor
      rate_eval % unscreened_rates(i_dscor_dt,482) = dscor_dt
      rate_eval % unscreened_rates(i_scor,484) = scor
      rate_eval % unscreened_rates(i_dscor_dt,484) = dscor_dt


      call screen5(pstate, 162, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,487) = scor
      rate_eval % unscreened_rates(i_dscor_dt,487) = dscor_dt
      rate_eval % unscreened_rates(i_scor,489) = scor
      rate_eval % unscreened_rates(i_dscor_dt,489) = dscor_dt


      call screen5(pstate, 163, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,488) = scor
      rate_eval % unscreened_rates(i_dscor_dt,488) = dscor_dt
      rate_eval % unscreened_rates(i_scor,490) = scor
      rate_eval % unscreened_rates(i_dscor_dt,490) = dscor_dt


      call screen5(pstate, 164, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,492) = scor
      rate_eval % unscreened_rates(i_dscor_dt,492) = dscor_dt
      rate_eval % unscreened_rates(i_scor,494) = scor
      rate_eval % unscreened_rates(i_dscor_dt,494) = dscor_dt


      call screen5(pstate, 165, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,493) = scor
      rate_eval % unscreened_rates(i_dscor_dt,493) = dscor_dt


      call screen5(pstate, 166, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,497) = scor
      rate_eval % unscreened_rates(i_dscor_dt,497) = dscor_dt
      rate_eval % unscreened_rates(i_scor,499) = scor
      rate_eval % unscreened_rates(i_dscor_dt,499) = dscor_dt


      call screen5(pstate, 167, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,498) = scor
      rate_eval % unscreened_rates(i_dscor_dt,498) = dscor_dt
      rate_eval % unscreened_rates(i_scor,500) = scor
      rate_eval % unscreened_rates(i_dscor_dt,500) = dscor_dt


      call screen5(pstate, 168, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,503) = scor
      rate_eval % unscreened_rates(i_dscor_dt,503) = dscor_dt
      rate_eval % unscreened_rates(i_scor,504) = scor
      rate_eval % unscreened_rates(i_dscor_dt,504) = dscor_dt


      call screen5(pstate, 169, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,509) = scor
      rate_eval % unscreened_rates(i_dscor_dt,509) = dscor_dt
      rate_eval % unscreened_rates(i_scor,511) = scor
      rate_eval % unscreened_rates(i_dscor_dt,511) = dscor_dt


      call screen5(pstate, 170, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,510) = scor
      rate_eval % unscreened_rates(i_dscor_dt,510) = dscor_dt
      rate_eval % unscreened_rates(i_scor,512) = scor
      rate_eval % unscreened_rates(i_dscor_dt,512) = dscor_dt


      call screen5(pstate, 171, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,517) = scor
      rate_eval % unscreened_rates(i_dscor_dt,517) = dscor_dt
      rate_eval % unscreened_rates(i_scor,519) = scor
      rate_eval % unscreened_rates(i_dscor_dt,519) = dscor_dt


      call screen5(pstate, 172, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,518) = scor
      rate_eval % unscreened_rates(i_dscor_dt,518) = dscor_dt
      rate_eval % unscreened_rates(i_scor,520) = scor
      rate_eval % unscreened_rates(i_dscor_dt,520) = dscor_dt


      call screen5(pstate, 173, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,524) = scor
      rate_eval % unscreened_rates(i_dscor_dt,524) = dscor_dt
      rate_eval % unscreened_rates(i_scor,526) = scor
      rate_eval % unscreened_rates(i_dscor_dt,526) = dscor_dt


      call screen5(pstate, 174, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,525) = scor
      rate_eval % unscreened_rates(i_dscor_dt,525) = dscor_dt
      rate_eval % unscreened_rates(i_scor,527) = scor
      rate_eval % unscreened_rates(i_dscor_dt,527) = dscor_dt


      call screen5(pstate, 175, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,532) = scor
      rate_eval % unscreened_rates(i_dscor_dt,532) = dscor_dt


      call screen5(pstate, 176, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,533) = scor
      rate_eval % unscreened_rates(i_dscor_dt,533) = dscor_dt
      rate_eval % unscreened_rates(i_scor,534) = scor
      rate_eval % unscreened_rates(i_dscor_dt,534) = dscor_dt


      call screen5(pstate, 177, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,539) = scor
      rate_eval % unscreened_rates(i_dscor_dt,539) = dscor_dt
      rate_eval % unscreened_rates(i_scor,541) = scor
      rate_eval % unscreened_rates(i_dscor_dt,541) = dscor_dt


      call screen5(pstate, 178, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,540) = scor
      rate_eval % unscreened_rates(i_dscor_dt,540) = dscor_dt
      rate_eval % unscreened_rates(i_scor,542) = scor
      rate_eval % unscreened_rates(i_dscor_dt,542) = dscor_dt


      call screen5(pstate, 179, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,548) = scor
      rate_eval % unscreened_rates(i_dscor_dt,548) = dscor_dt
      rate_eval % unscreened_rates(i_scor,550) = scor
      rate_eval % unscreened_rates(i_dscor_dt,550) = dscor_dt


      call screen5(pstate, 180, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,549) = scor
      rate_eval % unscreened_rates(i_dscor_dt,549) = dscor_dt
      rate_eval % unscreened_rates(i_scor,551) = scor
      rate_eval % unscreened_rates(i_dscor_dt,551) = dscor_dt


      call screen5(pstate, 181, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,555) = scor
      rate_eval % unscreened_rates(i_dscor_dt,555) = dscor_dt
      rate_eval % unscreened_rates(i_scor,557) = scor
      rate_eval % unscreened_rates(i_dscor_dt,557) = dscor_dt


      call screen5(pstate, 182, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,556) = scor
      rate_eval % unscreened_rates(i_dscor_dt,556) = dscor_dt
      rate_eval % unscreened_rates(i_scor,558) = scor
      rate_eval % unscreened_rates(i_dscor_dt,558) = dscor_dt


      call screen5(pstate, 183, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,561) = scor
      rate_eval % unscreened_rates(i_dscor_dt,561) = dscor_dt
      rate_eval % unscreened_rates(i_scor,563) = scor
      rate_eval % unscreened_rates(i_dscor_dt,563) = dscor_dt


      call screen5(pstate, 184, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,562) = scor
      rate_eval % unscreened_rates(i_dscor_dt,562) = dscor_dt
      rate_eval % unscreened_rates(i_scor,564) = scor
      rate_eval % unscreened_rates(i_dscor_dt,564) = dscor_dt


      call screen5(pstate, 185, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,568) = scor
      rate_eval % unscreened_rates(i_dscor_dt,568) = dscor_dt
      rate_eval % unscreened_rates(i_scor,570) = scor
      rate_eval % unscreened_rates(i_dscor_dt,570) = dscor_dt


      call screen5(pstate, 186, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,569) = scor
      rate_eval % unscreened_rates(i_dscor_dt,569) = dscor_dt
      rate_eval % unscreened_rates(i_scor,571) = scor
      rate_eval % unscreened_rates(i_dscor_dt,571) = dscor_dt


      call screen5(pstate, 187, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,573) = scor
      rate_eval % unscreened_rates(i_dscor_dt,573) = dscor_dt
      rate_eval % unscreened_rates(i_scor,575) = scor
      rate_eval % unscreened_rates(i_dscor_dt,575) = dscor_dt


      call screen5(pstate, 188, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,574) = scor
      rate_eval % unscreened_rates(i_dscor_dt,574) = dscor_dt
      rate_eval % unscreened_rates(i_scor,576) = scor
      rate_eval % unscreened_rates(i_dscor_dt,576) = dscor_dt


      call screen5(pstate, 189, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,579) = scor
      rate_eval % unscreened_rates(i_dscor_dt,579) = dscor_dt
      rate_eval % unscreened_rates(i_scor,581) = scor
      rate_eval % unscreened_rates(i_dscor_dt,581) = dscor_dt


      call screen5(pstate, 190, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,580) = scor
      rate_eval % unscreened_rates(i_dscor_dt,580) = dscor_dt
      rate_eval % unscreened_rates(i_scor,582) = scor
      rate_eval % unscreened_rates(i_dscor_dt,582) = dscor_dt


      call screen5(pstate, 191, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,586) = scor
      rate_eval % unscreened_rates(i_dscor_dt,586) = dscor_dt
      rate_eval % unscreened_rates(i_scor,588) = scor
      rate_eval % unscreened_rates(i_dscor_dt,588) = dscor_dt


      call screen5(pstate, 192, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,587) = scor
      rate_eval % unscreened_rates(i_dscor_dt,587) = dscor_dt
      rate_eval % unscreened_rates(i_scor,589) = scor
      rate_eval % unscreened_rates(i_dscor_dt,589) = dscor_dt


      call screen5(pstate, 193, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,593) = scor
      rate_eval % unscreened_rates(i_dscor_dt,593) = dscor_dt
      rate_eval % unscreened_rates(i_scor,595) = scor
      rate_eval % unscreened_rates(i_dscor_dt,595) = dscor_dt


      call screen5(pstate, 194, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,594) = scor
      rate_eval % unscreened_rates(i_dscor_dt,594) = dscor_dt
      rate_eval % unscreened_rates(i_scor,596) = scor
      rate_eval % unscreened_rates(i_dscor_dt,596) = dscor_dt


      call screen5(pstate, 195, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,600) = scor
      rate_eval % unscreened_rates(i_dscor_dt,600) = dscor_dt
      rate_eval % unscreened_rates(i_scor,602) = scor
      rate_eval % unscreened_rates(i_dscor_dt,602) = dscor_dt


      call screen5(pstate, 196, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,601) = scor
      rate_eval % unscreened_rates(i_dscor_dt,601) = dscor_dt
      rate_eval % unscreened_rates(i_scor,603) = scor
      rate_eval % unscreened_rates(i_dscor_dt,603) = dscor_dt


      call screen5(pstate, 197, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,605) = scor
      rate_eval % unscreened_rates(i_dscor_dt,605) = dscor_dt
      rate_eval % unscreened_rates(i_scor,607) = scor
      rate_eval % unscreened_rates(i_dscor_dt,607) = dscor_dt


      call screen5(pstate, 198, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,606) = scor
      rate_eval % unscreened_rates(i_dscor_dt,606) = dscor_dt
      rate_eval % unscreened_rates(i_scor,608) = scor
      rate_eval % unscreened_rates(i_dscor_dt,608) = dscor_dt


      call screen5(pstate, 199, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,612) = scor
      rate_eval % unscreened_rates(i_dscor_dt,612) = dscor_dt
      rate_eval % unscreened_rates(i_scor,614) = scor
      rate_eval % unscreened_rates(i_dscor_dt,614) = dscor_dt


      call screen5(pstate, 200, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,613) = scor
      rate_eval % unscreened_rates(i_dscor_dt,613) = dscor_dt
      rate_eval % unscreened_rates(i_scor,615) = scor
      rate_eval % unscreened_rates(i_dscor_dt,615) = dscor_dt


      call screen5(pstate, 201, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,619) = scor
      rate_eval % unscreened_rates(i_dscor_dt,619) = dscor_dt


      call screen5(pstate, 202, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,620) = scor
      rate_eval % unscreened_rates(i_dscor_dt,620) = dscor_dt
      rate_eval % unscreened_rates(i_scor,621) = scor
      rate_eval % unscreened_rates(i_dscor_dt,621) = dscor_dt


      call screen5(pstate, 203, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,626) = scor
      rate_eval % unscreened_rates(i_dscor_dt,626) = dscor_dt
      rate_eval % unscreened_rates(i_scor,628) = scor
      rate_eval % unscreened_rates(i_dscor_dt,628) = dscor_dt


      call screen5(pstate, 204, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,627) = scor
      rate_eval % unscreened_rates(i_dscor_dt,627) = dscor_dt
      rate_eval % unscreened_rates(i_scor,629) = scor
      rate_eval % unscreened_rates(i_dscor_dt,629) = dscor_dt
      rate_eval % unscreened_rates(i_scor,630) = scor
      rate_eval % unscreened_rates(i_dscor_dt,630) = dscor_dt


      call screen5(pstate, 205, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,634) = scor
      rate_eval % unscreened_rates(i_dscor_dt,634) = dscor_dt
      rate_eval % unscreened_rates(i_scor,636) = scor
      rate_eval % unscreened_rates(i_dscor_dt,636) = dscor_dt


      call screen5(pstate, 206, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,635) = scor
      rate_eval % unscreened_rates(i_dscor_dt,635) = dscor_dt
      rate_eval % unscreened_rates(i_scor,637) = scor
      rate_eval % unscreened_rates(i_dscor_dt,637) = dscor_dt


      call screen5(pstate, 207, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,641) = scor
      rate_eval % unscreened_rates(i_dscor_dt,641) = dscor_dt
      rate_eval % unscreened_rates(i_scor,643) = scor
      rate_eval % unscreened_rates(i_dscor_dt,643) = dscor_dt


      call screen5(pstate, 208, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,642) = scor
      rate_eval % unscreened_rates(i_dscor_dt,642) = dscor_dt
      rate_eval % unscreened_rates(i_scor,644) = scor
      rate_eval % unscreened_rates(i_dscor_dt,644) = dscor_dt


      call screen5(pstate, 209, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,649) = scor
      rate_eval % unscreened_rates(i_dscor_dt,649) = dscor_dt
      rate_eval % unscreened_rates(i_scor,650) = scor
      rate_eval % unscreened_rates(i_dscor_dt,650) = dscor_dt


      call screen5(pstate, 210, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,654) = scor
      rate_eval % unscreened_rates(i_dscor_dt,654) = dscor_dt
      rate_eval % unscreened_rates(i_scor,656) = scor
      rate_eval % unscreened_rates(i_dscor_dt,656) = dscor_dt


      call screen5(pstate, 211, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,655) = scor
      rate_eval % unscreened_rates(i_dscor_dt,655) = dscor_dt
      rate_eval % unscreened_rates(i_scor,657) = scor
      rate_eval % unscreened_rates(i_dscor_dt,657) = dscor_dt


      call screen5(pstate, 212, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,662) = scor
      rate_eval % unscreened_rates(i_dscor_dt,662) = dscor_dt
      rate_eval % unscreened_rates(i_scor,664) = scor
      rate_eval % unscreened_rates(i_dscor_dt,664) = dscor_dt


      call screen5(pstate, 213, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,663) = scor
      rate_eval % unscreened_rates(i_dscor_dt,663) = dscor_dt
      rate_eval % unscreened_rates(i_scor,665) = scor
      rate_eval % unscreened_rates(i_dscor_dt,665) = dscor_dt


      call screen5(pstate, 214, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,671) = scor
      rate_eval % unscreened_rates(i_dscor_dt,671) = dscor_dt
      rate_eval % unscreened_rates(i_scor,673) = scor
      rate_eval % unscreened_rates(i_dscor_dt,673) = dscor_dt


      call screen5(pstate, 215, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,672) = scor
      rate_eval % unscreened_rates(i_dscor_dt,672) = dscor_dt
      rate_eval % unscreened_rates(i_scor,674) = scor
      rate_eval % unscreened_rates(i_dscor_dt,674) = dscor_dt


      call screen5(pstate, 216, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,678) = scor
      rate_eval % unscreened_rates(i_dscor_dt,678) = dscor_dt
      rate_eval % unscreened_rates(i_scor,680) = scor
      rate_eval % unscreened_rates(i_dscor_dt,680) = dscor_dt


      call screen5(pstate, 217, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,679) = scor
      rate_eval % unscreened_rates(i_dscor_dt,679) = dscor_dt
      rate_eval % unscreened_rates(i_scor,681) = scor
      rate_eval % unscreened_rates(i_dscor_dt,681) = dscor_dt


      call screen5(pstate, 218, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,685) = scor
      rate_eval % unscreened_rates(i_dscor_dt,685) = dscor_dt
      rate_eval % unscreened_rates(i_scor,687) = scor
      rate_eval % unscreened_rates(i_dscor_dt,687) = dscor_dt


      call screen5(pstate, 219, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,686) = scor
      rate_eval % unscreened_rates(i_dscor_dt,686) = dscor_dt
      rate_eval % unscreened_rates(i_scor,688) = scor
      rate_eval % unscreened_rates(i_dscor_dt,688) = dscor_dt


      call screen5(pstate, 220, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,692) = scor
      rate_eval % unscreened_rates(i_dscor_dt,692) = dscor_dt
      rate_eval % unscreened_rates(i_scor,694) = scor
      rate_eval % unscreened_rates(i_dscor_dt,694) = dscor_dt


      call screen5(pstate, 221, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,693) = scor
      rate_eval % unscreened_rates(i_dscor_dt,693) = dscor_dt
      rate_eval % unscreened_rates(i_scor,695) = scor
      rate_eval % unscreened_rates(i_dscor_dt,695) = dscor_dt


      call screen5(pstate, 222, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,696) = scor
      rate_eval % unscreened_rates(i_dscor_dt,696) = dscor_dt
      rate_eval % unscreened_rates(i_scor,698) = scor
      rate_eval % unscreened_rates(i_dscor_dt,698) = dscor_dt


      call screen5(pstate, 223, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,697) = scor
      rate_eval % unscreened_rates(i_dscor_dt,697) = dscor_dt


      call screen5(pstate, 224, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,701) = scor
      rate_eval % unscreened_rates(i_dscor_dt,701) = dscor_dt
      rate_eval % unscreened_rates(i_scor,703) = scor
      rate_eval % unscreened_rates(i_dscor_dt,703) = dscor_dt


      call screen5(pstate, 225, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,702) = scor
      rate_eval % unscreened_rates(i_dscor_dt,702) = dscor_dt
      rate_eval % unscreened_rates(i_scor,704) = scor
      rate_eval % unscreened_rates(i_dscor_dt,704) = dscor_dt


      call screen5(pstate, 226, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,707) = scor
      rate_eval % unscreened_rates(i_dscor_dt,707) = dscor_dt
      rate_eval % unscreened_rates(i_scor,709) = scor
      rate_eval % unscreened_rates(i_dscor_dt,709) = dscor_dt


      call screen5(pstate, 227, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,708) = scor
      rate_eval % unscreened_rates(i_dscor_dt,708) = dscor_dt
      rate_eval % unscreened_rates(i_scor,710) = scor
      rate_eval % unscreened_rates(i_dscor_dt,710) = dscor_dt


      call screen5(pstate, 228, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,714) = scor
      rate_eval % unscreened_rates(i_dscor_dt,714) = dscor_dt
      rate_eval % unscreened_rates(i_scor,716) = scor
      rate_eval % unscreened_rates(i_dscor_dt,716) = dscor_dt


      call screen5(pstate, 229, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,715) = scor
      rate_eval % unscreened_rates(i_dscor_dt,715) = dscor_dt
      rate_eval % unscreened_rates(i_scor,717) = scor
      rate_eval % unscreened_rates(i_dscor_dt,717) = dscor_dt


      call screen5(pstate, 230, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,721) = scor
      rate_eval % unscreened_rates(i_dscor_dt,721) = dscor_dt
      rate_eval % unscreened_rates(i_scor,723) = scor
      rate_eval % unscreened_rates(i_dscor_dt,723) = dscor_dt


      call screen5(pstate, 231, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,722) = scor
      rate_eval % unscreened_rates(i_dscor_dt,722) = dscor_dt
      rate_eval % unscreened_rates(i_scor,724) = scor
      rate_eval % unscreened_rates(i_dscor_dt,724) = dscor_dt


      call screen5(pstate, 232, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,728) = scor
      rate_eval % unscreened_rates(i_dscor_dt,728) = dscor_dt
      rate_eval % unscreened_rates(i_scor,730) = scor
      rate_eval % unscreened_rates(i_dscor_dt,730) = dscor_dt


      call screen5(pstate, 233, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,729) = scor
      rate_eval % unscreened_rates(i_dscor_dt,729) = dscor_dt
      rate_eval % unscreened_rates(i_scor,731) = scor
      rate_eval % unscreened_rates(i_dscor_dt,731) = dscor_dt


      call screen5(pstate, 234, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,732) = scor
      rate_eval % unscreened_rates(i_dscor_dt,732) = dscor_dt
      rate_eval % unscreened_rates(i_scor,734) = scor
      rate_eval % unscreened_rates(i_dscor_dt,734) = dscor_dt


      call screen5(pstate, 235, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,733) = scor
      rate_eval % unscreened_rates(i_dscor_dt,733) = dscor_dt


      call screen5(pstate, 236, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,737) = scor
      rate_eval % unscreened_rates(i_dscor_dt,737) = dscor_dt
      rate_eval % unscreened_rates(i_scor,739) = scor
      rate_eval % unscreened_rates(i_dscor_dt,739) = dscor_dt


      call screen5(pstate, 237, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,738) = scor
      rate_eval % unscreened_rates(i_dscor_dt,738) = dscor_dt
      rate_eval % unscreened_rates(i_scor,740) = scor
      rate_eval % unscreened_rates(i_dscor_dt,740) = dscor_dt


      call screen5(pstate, 238, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,744) = scor
      rate_eval % unscreened_rates(i_dscor_dt,744) = dscor_dt
      rate_eval % unscreened_rates(i_scor,746) = scor
      rate_eval % unscreened_rates(i_dscor_dt,746) = dscor_dt


      call screen5(pstate, 239, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,745) = scor
      rate_eval % unscreened_rates(i_dscor_dt,745) = dscor_dt
      rate_eval % unscreened_rates(i_scor,747) = scor
      rate_eval % unscreened_rates(i_dscor_dt,747) = dscor_dt


      call screen5(pstate, 240, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,750) = scor
      rate_eval % unscreened_rates(i_dscor_dt,750) = dscor_dt
      rate_eval % unscreened_rates(i_scor,751) = scor
      rate_eval % unscreened_rates(i_dscor_dt,751) = dscor_dt


      call screen5(pstate, 241, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,755) = scor
      rate_eval % unscreened_rates(i_dscor_dt,755) = dscor_dt
      rate_eval % unscreened_rates(i_scor,757) = scor
      rate_eval % unscreened_rates(i_dscor_dt,757) = dscor_dt


      call screen5(pstate, 242, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,756) = scor
      rate_eval % unscreened_rates(i_dscor_dt,756) = dscor_dt
      rate_eval % unscreened_rates(i_scor,758) = scor
      rate_eval % unscreened_rates(i_dscor_dt,758) = dscor_dt


      call screen5(pstate, 243, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,763) = scor
      rate_eval % unscreened_rates(i_dscor_dt,763) = dscor_dt
      rate_eval % unscreened_rates(i_scor,765) = scor
      rate_eval % unscreened_rates(i_dscor_dt,765) = dscor_dt


      call screen5(pstate, 244, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,764) = scor
      rate_eval % unscreened_rates(i_dscor_dt,764) = dscor_dt
      rate_eval % unscreened_rates(i_scor,766) = scor
      rate_eval % unscreened_rates(i_dscor_dt,766) = dscor_dt


      call screen5(pstate, 245, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,770) = scor
      rate_eval % unscreened_rates(i_dscor_dt,770) = dscor_dt
      rate_eval % unscreened_rates(i_scor,772) = scor
      rate_eval % unscreened_rates(i_dscor_dt,772) = dscor_dt


      call screen5(pstate, 246, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,771) = scor
      rate_eval % unscreened_rates(i_dscor_dt,771) = dscor_dt
      rate_eval % unscreened_rates(i_scor,773) = scor
      rate_eval % unscreened_rates(i_dscor_dt,773) = dscor_dt


      call screen5(pstate, 247, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,777) = scor
      rate_eval % unscreened_rates(i_dscor_dt,777) = dscor_dt
      rate_eval % unscreened_rates(i_scor,779) = scor
      rate_eval % unscreened_rates(i_dscor_dt,779) = dscor_dt


      call screen5(pstate, 248, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,778) = scor
      rate_eval % unscreened_rates(i_dscor_dt,778) = dscor_dt
      rate_eval % unscreened_rates(i_scor,780) = scor
      rate_eval % unscreened_rates(i_dscor_dt,780) = dscor_dt


      call screen5(pstate, 249, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,786) = scor
      rate_eval % unscreened_rates(i_dscor_dt,786) = dscor_dt


      call screen5(pstate, 250, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,787) = scor
      rate_eval % unscreened_rates(i_dscor_dt,787) = dscor_dt
      rate_eval % unscreened_rates(i_scor,788) = scor
      rate_eval % unscreened_rates(i_dscor_dt,788) = dscor_dt


      call screen5(pstate, 251, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,792) = scor
      rate_eval % unscreened_rates(i_dscor_dt,792) = dscor_dt
      rate_eval % unscreened_rates(i_scor,794) = scor
      rate_eval % unscreened_rates(i_dscor_dt,794) = dscor_dt


      call screen5(pstate, 252, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,793) = scor
      rate_eval % unscreened_rates(i_dscor_dt,793) = dscor_dt
      rate_eval % unscreened_rates(i_scor,795) = scor
      rate_eval % unscreened_rates(i_dscor_dt,795) = dscor_dt


      call screen5(pstate, 253, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,800) = scor
      rate_eval % unscreened_rates(i_dscor_dt,800) = dscor_dt
      rate_eval % unscreened_rates(i_scor,802) = scor
      rate_eval % unscreened_rates(i_dscor_dt,802) = dscor_dt


      call screen5(pstate, 254, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,801) = scor
      rate_eval % unscreened_rates(i_dscor_dt,801) = dscor_dt
      rate_eval % unscreened_rates(i_scor,803) = scor
      rate_eval % unscreened_rates(i_dscor_dt,803) = dscor_dt


      call screen5(pstate, 255, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,807) = scor
      rate_eval % unscreened_rates(i_dscor_dt,807) = dscor_dt
      rate_eval % unscreened_rates(i_scor,809) = scor
      rate_eval % unscreened_rates(i_dscor_dt,809) = dscor_dt


      call screen5(pstate, 256, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,808) = scor
      rate_eval % unscreened_rates(i_dscor_dt,808) = dscor_dt
      rate_eval % unscreened_rates(i_scor,810) = scor
      rate_eval % unscreened_rates(i_dscor_dt,810) = dscor_dt


      call screen5(pstate, 257, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,814) = scor
      rate_eval % unscreened_rates(i_dscor_dt,814) = dscor_dt
      rate_eval % unscreened_rates(i_scor,816) = scor
      rate_eval % unscreened_rates(i_dscor_dt,816) = dscor_dt


      call screen5(pstate, 258, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,815) = scor
      rate_eval % unscreened_rates(i_dscor_dt,815) = dscor_dt
      rate_eval % unscreened_rates(i_scor,817) = scor
      rate_eval % unscreened_rates(i_dscor_dt,817) = dscor_dt


      call screen5(pstate, 259, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,821) = scor
      rate_eval % unscreened_rates(i_dscor_dt,821) = dscor_dt
      rate_eval % unscreened_rates(i_scor,823) = scor
      rate_eval % unscreened_rates(i_dscor_dt,823) = dscor_dt


      call screen5(pstate, 260, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,822) = scor
      rate_eval % unscreened_rates(i_dscor_dt,822) = dscor_dt
      rate_eval % unscreened_rates(i_scor,824) = scor
      rate_eval % unscreened_rates(i_dscor_dt,824) = dscor_dt


      call screen5(pstate, 261, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,828) = scor
      rate_eval % unscreened_rates(i_dscor_dt,828) = dscor_dt
      rate_eval % unscreened_rates(i_scor,830) = scor
      rate_eval % unscreened_rates(i_dscor_dt,830) = dscor_dt


      call screen5(pstate, 262, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,829) = scor
      rate_eval % unscreened_rates(i_dscor_dt,829) = dscor_dt
      rate_eval % unscreened_rates(i_scor,831) = scor
      rate_eval % unscreened_rates(i_dscor_dt,831) = dscor_dt


      call screen5(pstate, 263, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,833) = scor
      rate_eval % unscreened_rates(i_dscor_dt,833) = dscor_dt
      rate_eval % unscreened_rates(i_scor,835) = scor
      rate_eval % unscreened_rates(i_dscor_dt,835) = dscor_dt


      call screen5(pstate, 264, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,834) = scor
      rate_eval % unscreened_rates(i_dscor_dt,834) = dscor_dt


      call screen5(pstate, 265, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,838) = scor
      rate_eval % unscreened_rates(i_dscor_dt,838) = dscor_dt
      rate_eval % unscreened_rates(i_scor,840) = scor
      rate_eval % unscreened_rates(i_dscor_dt,840) = dscor_dt


      call screen5(pstate, 266, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,839) = scor
      rate_eval % unscreened_rates(i_dscor_dt,839) = dscor_dt
      rate_eval % unscreened_rates(i_scor,841) = scor
      rate_eval % unscreened_rates(i_dscor_dt,841) = dscor_dt


      call screen5(pstate, 267, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,844) = scor
      rate_eval % unscreened_rates(i_dscor_dt,844) = dscor_dt
      rate_eval % unscreened_rates(i_scor,846) = scor
      rate_eval % unscreened_rates(i_dscor_dt,846) = dscor_dt


      call screen5(pstate, 268, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,845) = scor
      rate_eval % unscreened_rates(i_dscor_dt,845) = dscor_dt
      rate_eval % unscreened_rates(i_scor,847) = scor
      rate_eval % unscreened_rates(i_dscor_dt,847) = dscor_dt


      call screen5(pstate, 269, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,851) = scor
      rate_eval % unscreened_rates(i_dscor_dt,851) = dscor_dt
      rate_eval % unscreened_rates(i_scor,853) = scor
      rate_eval % unscreened_rates(i_dscor_dt,853) = dscor_dt


      call screen5(pstate, 270, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,852) = scor
      rate_eval % unscreened_rates(i_dscor_dt,852) = dscor_dt
      rate_eval % unscreened_rates(i_scor,854) = scor
      rate_eval % unscreened_rates(i_dscor_dt,854) = dscor_dt


      call screen5(pstate, 271, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,858) = scor
      rate_eval % unscreened_rates(i_dscor_dt,858) = dscor_dt
      rate_eval % unscreened_rates(i_scor,860) = scor
      rate_eval % unscreened_rates(i_dscor_dt,860) = dscor_dt


      call screen5(pstate, 272, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,859) = scor
      rate_eval % unscreened_rates(i_dscor_dt,859) = dscor_dt
      rate_eval % unscreened_rates(i_scor,861) = scor
      rate_eval % unscreened_rates(i_dscor_dt,861) = dscor_dt


      call screen5(pstate, 273, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,865) = scor
      rate_eval % unscreened_rates(i_dscor_dt,865) = dscor_dt
      rate_eval % unscreened_rates(i_scor,867) = scor
      rate_eval % unscreened_rates(i_dscor_dt,867) = dscor_dt


      call screen5(pstate, 274, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,866) = scor
      rate_eval % unscreened_rates(i_dscor_dt,866) = dscor_dt
      rate_eval % unscreened_rates(i_scor,868) = scor
      rate_eval % unscreened_rates(i_dscor_dt,868) = dscor_dt


      call screen5(pstate, 275, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,870) = scor
      rate_eval % unscreened_rates(i_dscor_dt,870) = dscor_dt
      rate_eval % unscreened_rates(i_scor,872) = scor
      rate_eval % unscreened_rates(i_dscor_dt,872) = dscor_dt


      call screen5(pstate, 276, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,871) = scor
      rate_eval % unscreened_rates(i_dscor_dt,871) = dscor_dt


      call screen5(pstate, 277, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,875) = scor
      rate_eval % unscreened_rates(i_dscor_dt,875) = dscor_dt
      rate_eval % unscreened_rates(i_scor,877) = scor
      rate_eval % unscreened_rates(i_dscor_dt,877) = dscor_dt


      call screen5(pstate, 278, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,876) = scor
      rate_eval % unscreened_rates(i_dscor_dt,876) = dscor_dt
      rate_eval % unscreened_rates(i_scor,878) = scor
      rate_eval % unscreened_rates(i_dscor_dt,878) = dscor_dt


      call screen5(pstate, 279, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,882) = scor
      rate_eval % unscreened_rates(i_dscor_dt,882) = dscor_dt
      rate_eval % unscreened_rates(i_scor,884) = scor
      rate_eval % unscreened_rates(i_dscor_dt,884) = dscor_dt


      call screen5(pstate, 280, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,883) = scor
      rate_eval % unscreened_rates(i_dscor_dt,883) = dscor_dt
      rate_eval % unscreened_rates(i_scor,885) = scor
      rate_eval % unscreened_rates(i_dscor_dt,885) = dscor_dt


      call screen5(pstate, 281, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,889) = scor
      rate_eval % unscreened_rates(i_dscor_dt,889) = dscor_dt
      rate_eval % unscreened_rates(i_scor,890) = scor
      rate_eval % unscreened_rates(i_dscor_dt,890) = dscor_dt


      call screen5(pstate, 282, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,894) = scor
      rate_eval % unscreened_rates(i_dscor_dt,894) = dscor_dt
      rate_eval % unscreened_rates(i_scor,896) = scor
      rate_eval % unscreened_rates(i_dscor_dt,896) = dscor_dt


      call screen5(pstate, 283, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,895) = scor
      rate_eval % unscreened_rates(i_dscor_dt,895) = dscor_dt
      rate_eval % unscreened_rates(i_scor,897) = scor
      rate_eval % unscreened_rates(i_dscor_dt,897) = dscor_dt


      call screen5(pstate, 284, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,902) = scor
      rate_eval % unscreened_rates(i_dscor_dt,902) = dscor_dt
      rate_eval % unscreened_rates(i_scor,904) = scor
      rate_eval % unscreened_rates(i_dscor_dt,904) = dscor_dt


      call screen5(pstate, 285, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,903) = scor
      rate_eval % unscreened_rates(i_dscor_dt,903) = dscor_dt
      rate_eval % unscreened_rates(i_scor,905) = scor
      rate_eval % unscreened_rates(i_dscor_dt,905) = dscor_dt


      call screen5(pstate, 286, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,909) = scor
      rate_eval % unscreened_rates(i_dscor_dt,909) = dscor_dt
      rate_eval % unscreened_rates(i_scor,911) = scor
      rate_eval % unscreened_rates(i_dscor_dt,911) = dscor_dt


      call screen5(pstate, 287, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,910) = scor
      rate_eval % unscreened_rates(i_dscor_dt,910) = dscor_dt
      rate_eval % unscreened_rates(i_scor,912) = scor
      rate_eval % unscreened_rates(i_dscor_dt,912) = dscor_dt


      call screen5(pstate, 288, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,916) = scor
      rate_eval % unscreened_rates(i_dscor_dt,916) = dscor_dt
      rate_eval % unscreened_rates(i_scor,918) = scor
      rate_eval % unscreened_rates(i_dscor_dt,918) = dscor_dt


      call screen5(pstate, 289, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,917) = scor
      rate_eval % unscreened_rates(i_dscor_dt,917) = dscor_dt
      rate_eval % unscreened_rates(i_scor,919) = scor
      rate_eval % unscreened_rates(i_dscor_dt,919) = dscor_dt


      call screen5(pstate, 290, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,922) = scor
      rate_eval % unscreened_rates(i_dscor_dt,922) = dscor_dt
      rate_eval % unscreened_rates(i_scor,923) = scor
      rate_eval % unscreened_rates(i_dscor_dt,923) = dscor_dt


      call screen5(pstate, 291, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,928) = scor
      rate_eval % unscreened_rates(i_dscor_dt,928) = dscor_dt
      rate_eval % unscreened_rates(i_scor,930) = scor
      rate_eval % unscreened_rates(i_dscor_dt,930) = dscor_dt


      call screen5(pstate, 292, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,929) = scor
      rate_eval % unscreened_rates(i_dscor_dt,929) = dscor_dt
      rate_eval % unscreened_rates(i_scor,931) = scor
      rate_eval % unscreened_rates(i_dscor_dt,931) = dscor_dt


      call screen5(pstate, 293, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,935) = scor
      rate_eval % unscreened_rates(i_dscor_dt,935) = dscor_dt
      rate_eval % unscreened_rates(i_scor,937) = scor
      rate_eval % unscreened_rates(i_dscor_dt,937) = dscor_dt


      call screen5(pstate, 294, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,936) = scor
      rate_eval % unscreened_rates(i_dscor_dt,936) = dscor_dt
      rate_eval % unscreened_rates(i_scor,938) = scor
      rate_eval % unscreened_rates(i_dscor_dt,938) = dscor_dt


      call screen5(pstate, 295, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,943) = scor
      rate_eval % unscreened_rates(i_dscor_dt,943) = dscor_dt
      rate_eval % unscreened_rates(i_scor,945) = scor
      rate_eval % unscreened_rates(i_dscor_dt,945) = dscor_dt


      call screen5(pstate, 296, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,944) = scor
      rate_eval % unscreened_rates(i_dscor_dt,944) = dscor_dt
      rate_eval % unscreened_rates(i_scor,946) = scor
      rate_eval % unscreened_rates(i_dscor_dt,946) = dscor_dt


      call screen5(pstate, 297, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,951) = scor
      rate_eval % unscreened_rates(i_dscor_dt,951) = dscor_dt
      rate_eval % unscreened_rates(i_scor,953) = scor
      rate_eval % unscreened_rates(i_dscor_dt,953) = dscor_dt


      call screen5(pstate, 298, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,952) = scor
      rate_eval % unscreened_rates(i_dscor_dt,952) = dscor_dt
      rate_eval % unscreened_rates(i_scor,954) = scor
      rate_eval % unscreened_rates(i_dscor_dt,954) = dscor_dt


      call screen5(pstate, 299, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,958) = scor
      rate_eval % unscreened_rates(i_dscor_dt,958) = dscor_dt
      rate_eval % unscreened_rates(i_scor,960) = scor
      rate_eval % unscreened_rates(i_dscor_dt,960) = dscor_dt


      call screen5(pstate, 300, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,959) = scor
      rate_eval % unscreened_rates(i_dscor_dt,959) = dscor_dt
      rate_eval % unscreened_rates(i_scor,961) = scor
      rate_eval % unscreened_rates(i_dscor_dt,961) = dscor_dt


      call screen5(pstate, 301, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,965) = scor
      rate_eval % unscreened_rates(i_dscor_dt,965) = dscor_dt
      rate_eval % unscreened_rates(i_scor,967) = scor
      rate_eval % unscreened_rates(i_dscor_dt,967) = dscor_dt


      call screen5(pstate, 302, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,966) = scor
      rate_eval % unscreened_rates(i_dscor_dt,966) = dscor_dt
      rate_eval % unscreened_rates(i_scor,968) = scor
      rate_eval % unscreened_rates(i_dscor_dt,968) = dscor_dt


      call screen5(pstate, 303, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,972) = scor
      rate_eval % unscreened_rates(i_dscor_dt,972) = dscor_dt
      rate_eval % unscreened_rates(i_scor,973) = scor
      rate_eval % unscreened_rates(i_dscor_dt,973) = dscor_dt


      call screen5(pstate, 304, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,974) = scor
      rate_eval % unscreened_rates(i_dscor_dt,974) = dscor_dt


      call screen5(pstate, 305, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,976) = scor
      rate_eval % unscreened_rates(i_dscor_dt,976) = dscor_dt
      rate_eval % unscreened_rates(i_scor,978) = scor
      rate_eval % unscreened_rates(i_dscor_dt,978) = dscor_dt


      call screen5(pstate, 306, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,977) = scor
      rate_eval % unscreened_rates(i_dscor_dt,977) = dscor_dt
      rate_eval % unscreened_rates(i_scor,979) = scor
      rate_eval % unscreened_rates(i_dscor_dt,979) = dscor_dt


      call screen5(pstate, 307, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,982) = scor
      rate_eval % unscreened_rates(i_dscor_dt,982) = dscor_dt
      rate_eval % unscreened_rates(i_scor,983) = scor
      rate_eval % unscreened_rates(i_dscor_dt,983) = dscor_dt


      call screen5(pstate, 308, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,984) = scor
      rate_eval % unscreened_rates(i_dscor_dt,984) = dscor_dt


      call screen5(pstate, 309, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,987) = scor
      rate_eval % unscreened_rates(i_dscor_dt,987) = dscor_dt
      rate_eval % unscreened_rates(i_scor,988) = scor
      rate_eval % unscreened_rates(i_dscor_dt,988) = dscor_dt


      call screen5(pstate, 310, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,992) = scor
      rate_eval % unscreened_rates(i_dscor_dt,992) = dscor_dt
      rate_eval % unscreened_rates(i_scor,993) = scor
      rate_eval % unscreened_rates(i_dscor_dt,993) = dscor_dt


      call screen5(pstate, 311, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,994) = scor
      rate_eval % unscreened_rates(i_dscor_dt,994) = dscor_dt


      call screen5(pstate, 312, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,998) = scor
      rate_eval % unscreened_rates(i_dscor_dt,998) = dscor_dt
      rate_eval % unscreened_rates(i_scor,999) = scor
      rate_eval % unscreened_rates(i_dscor_dt,999) = dscor_dt


      call screen5(pstate, 313, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1003) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1003) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1004) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1004) = dscor_dt


      call screen5(pstate, 314, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1006) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1006) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1007) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1007) = dscor_dt


      call screen5(pstate, 315, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1011) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1011) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1012) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1012) = dscor_dt


      call screen5(pstate, 316, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1016) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1016) = dscor_dt


      call screen5(pstate, 317, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1021) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1021) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1022) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1022) = dscor_dt


      call screen5(pstate, 318, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1026) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1026) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1028) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1028) = dscor_dt


      call screen5(pstate, 319, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1027) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1027) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1029) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1029) = dscor_dt


      call screen5(pstate, 320, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1034) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1034) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1036) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1036) = dscor_dt


      call screen5(pstate, 321, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1035) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1035) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1037) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1037) = dscor_dt


      call screen5(pstate, 322, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1041) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1041) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1042) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1042) = dscor_dt


      call screen5(pstate, 323, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1043) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1043) = dscor_dt


      call screen5(pstate, 324, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1047) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1047) = dscor_dt


      call screen5(pstate, 325, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1051) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1051) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1052) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1052) = dscor_dt


      call screen5(pstate, 326, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1058) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1058) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1060) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1060) = dscor_dt


      call screen5(pstate, 327, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1059) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1059) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1061) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1061) = dscor_dt


      call screen5(pstate, 328, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1067) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1067) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1069) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1069) = dscor_dt


      call screen5(pstate, 329, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1068) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1068) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1070) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1070) = dscor_dt


      call screen5(pstate, 330, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1075) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1075) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1077) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1077) = dscor_dt


      call screen5(pstate, 331, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1076) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1076) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1078) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1078) = dscor_dt


      call screen5(pstate, 332, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1082) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1082) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1083) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1083) = dscor_dt


      call screen5(pstate, 333, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1084) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1084) = dscor_dt


      call screen5(pstate, 334, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1088) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1088) = dscor_dt


      call screen5(pstate, 335, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1092) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1092) = dscor_dt


      call screen5(pstate, 336, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1093) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1093) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1094) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1094) = dscor_dt


      call screen5(pstate, 337, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1097) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1097) = dscor_dt


      call screen5(pstate, 338, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1102) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1102) = dscor_dt


      call screen5(pstate, 339, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1103) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1103) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1104) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1104) = dscor_dt


      call screen5(pstate, 340, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1108) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1108) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1109) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1109) = dscor_dt


      call screen5(pstate, 341, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1110) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1110) = dscor_dt


      call screen5(pstate, 342, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1115) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1115) = dscor_dt


      call screen5(pstate, 343, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1119) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1119) = dscor_dt


      call screen5(pstate, 344, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1120) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1120) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1121) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1121) = dscor_dt


      call screen5(pstate, 345, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1126) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1126) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1128) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1128) = dscor_dt


      call screen5(pstate, 346, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1127) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1127) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1129) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1129) = dscor_dt


      call screen5(pstate, 347, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1134) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1134) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1135) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1135) = dscor_dt


      call screen5(pstate, 348, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1136) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1136) = dscor_dt


      call screen5(pstate, 349, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1141) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1141) = dscor_dt


      call screen5(pstate, 350, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1144) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1144) = dscor_dt


      call screen5(pstate, 351, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1145) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1145) = dscor_dt


      call screen5(pstate, 352, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1150) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1150) = dscor_dt


      call screen5(pstate, 353, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1155) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1155) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1156) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1156) = dscor_dt


      call screen5(pstate, 354, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1161) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1161) = dscor_dt
      rate_eval % unscreened_rates(i_scor,1162) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1162) = dscor_dt


      call screen5(pstate, 355, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1163) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1163) = dscor_dt


      call screen5(pstate, 356, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1168) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1168) = dscor_dt

    end if


    ! Compute screened rates
    rate_eval % screened_rates = rate_eval % unscreened_rates(i_rate, :) * &
                                 rate_eval % unscreened_rates(i_scor, :)

  end subroutine evaluate_rates


  subroutine actual_rhs(state)
    
    !$acc routine seq

    use extern_probin_module, only: do_constant_volume_burn
    use burn_type_module, only: net_itemp, net_ienuc
    use sneut_module, only: sneut5
    use temperature_integration_module, only: temperature_rhs

    implicit none

    type(burn_t) :: state
    type(rate_eval_t) :: rate_eval
    real(rt) :: Y(nspec), ydot_nuc(nspec)
    real(rt) :: reactvec(num_rate_groups+2)
    integer :: i, j
    real(rt) :: rhoy, ye, enuc
    real(rt) :: sneut, dsneutdt, dsneutdd, snuda, snudz

    !$gpu

    ! Set molar abundances
    Y(:) = state % xn(:) * aion_inv(:)

    call evaluate_rates(state, rate_eval)

    call rhs_nuc(state, ydot_nuc, Y, rate_eval % screened_rates)
    state % ydot(1:nspec) = ydot_nuc

    ! ion binding energy contributions
    call ener_gener_rate(ydot_nuc, enuc)

    ! additional per-reaction energies
    ! including Q-value modification and electron chemical potential

    ! additional energy generation rates
    ! including gamma heating and reaction neutrino losses (non-thermal)


    ! Get the thermal neutrino losses
    call sneut5(state % T, state % rho, state % abar, state % zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)

    ! Append the energy equation (this is erg/g/s)
    state % ydot(net_ienuc) = enuc - sneut

    ! Append the temperature equation
    call temperature_rhs(state)

  end subroutine actual_rhs


  subroutine rhs_nuc(state, ydot_nuc, Y, screened_rates)

    !$acc routine seq

    implicit none

    type (burn_t), intent(in) :: state
    real(rt), intent(out) :: ydot_nuc(nspec)
    real(rt), intent(in)  :: Y(nspec)
    real(rt), intent(in)  :: screened_rates(nrates)

    !$gpu



    ydot_nuc(jp) = ( &
      screened_rates(k_al22__p_mg21)*Y(jal22) + screened_rates(k_al22__p_na21__weak__wc17)* &
      Y(jal22) + 2.0d0*screened_rates(k_al22__p_p_ne20__weak__wc12)*Y(jal22) &
      + screened_rates(k_al23__p_mg22)*Y(jal23) + &
      screened_rates(k_al23__p_na22__weak__wc12)*Y(jal23) + &
      screened_rates(k_al24__p_mg23)*Y(jal24) + &
      screened_rates(k_al24__p_na23__weak__wc12)*Y(jal24) + &
      screened_rates(k_al25__p_mg24)*Y(jal25) + screened_rates(k_al26__p_mg25)* &
      Y(jal26) + screened_rates(k_al27__p_mg26)*Y(jal27) + &
      screened_rates(k_ar31__p_cl30)*Y(jar31) + 2.0d0* &
      screened_rates(k_ar31__p_p_p29__weak__wc12)*Y(jar31) + &
      screened_rates(k_ar31__p_s30__weak__wc17)*Y(jar31) + &
      screened_rates(k_ar32__p_cl31)*Y(jar32) + &
      screened_rates(k_ar32__p_s31__weak__wc12)*Y(jar32) + &
      screened_rates(k_ar33__p_cl32)*Y(jar33) + &
      screened_rates(k_ar33__p_s32__weak__wc12)*Y(jar33) + &
      screened_rates(k_ar34__p_cl33)*Y(jar34) + screened_rates(k_ar35__p_cl34)* &
      Y(jar35) + screened_rates(k_ar36__p_cl35)*Y(jar36) + &
      screened_rates(k_ar37__p_cl36)*Y(jar37) + screened_rates(k_ar38__p_cl37)* &
      Y(jar38) + screened_rates(k_b8__p_be7)*Y(jb8) + screened_rates(k_ca34__p_k33)* &
      Y(jca34) + screened_rates(k_ca35__p_ar34__weak__wc17)*Y(jca35) + &
      screened_rates(k_ca35__p_k34)*Y(jca35) + 2.0d0* &
      screened_rates(k_ca35__p_p_cl33__weak__wc12)*Y(jca35) + &
      screened_rates(k_ca36__p_ar35__weak__wc12)*Y(jca36) + &
      screened_rates(k_ca36__p_k35)*Y(jca36) + &
      screened_rates(k_ca37__p_ar36__weak__wc12)*Y(jca37) + &
      screened_rates(k_ca37__p_k36)*Y(jca37) + screened_rates(k_ca38__p_k37)*Y(jca38) &
      + screened_rates(k_ca39__p_k38)*Y(jca39) + screened_rates(k_ca40__p_k39)* &
      Y(jca40) + screened_rates(k_ca41__p_k40)*Y(jca41) + &
      screened_rates(k_ca42__p_k41)*Y(jca42) + screened_rates(k_cl29__p_s28)*Y(jcl29) &
      + screened_rates(k_cl30__p_s29)*Y(jcl30) + &
      screened_rates(k_cl31__p_p30__weak__wc12)*Y(jcl31) + &
      screened_rates(k_cl31__p_s30)*Y(jcl31) + &
      screened_rates(k_cl32__p_p31__weak__wc12)*Y(jcl32) + &
      screened_rates(k_cl32__p_s31)*Y(jcl32) + screened_rates(k_cl33__p_s32)*Y(jcl33) &
      + screened_rates(k_cl34__p_s33)*Y(jcl34) + screened_rates(k_cl35__p_s34)* &
      Y(jcl35) + screened_rates(k_cl36__p_s35)*Y(jcl36) + &
      screened_rates(k_co47__p_fe46)*Y(jco47) + screened_rates(k_co48__p_fe47)* &
      Y(jco48) + screened_rates(k_co49__p_fe48)*Y(jco49) + &
      screened_rates(k_co50__p_fe49)*Y(jco50) + &
      screened_rates(k_co50__p_mn49__weak__wc12)*Y(jco50) + &
      screened_rates(k_co51__p_fe50)*Y(jco51) + screened_rates(k_co52__p_fe51)* &
      Y(jco52) + screened_rates(k_co53__p_fe52)*Y(jco53) + &
      screened_rates(k_co54__p_fe53)*Y(jco54) + screened_rates(k_co55__p_fe54)* &
      Y(jco55) + screened_rates(k_co56__p_fe55)*Y(jco56) + &
      screened_rates(k_cr42__p_ti41__weak__wc12)*Y(jcr42) + &
      screened_rates(k_cr42__p_v41)*Y(jcr42) + 3.0d0* &
      screened_rates(k_cr43__p_p_p_ca40__weak__wc12)*Y(jcr43) + 2.0d0* &
      screened_rates(k_cr43__p_p_sc41__weak__wc12)*Y(jcr43) + &
      screened_rates(k_cr43__p_ti42__weak__wc17)*Y(jcr43) + &
      screened_rates(k_cr43__p_v42)*Y(jcr43) + &
      screened_rates(k_cr44__p_ti43__weak__wc12)*Y(jcr44) + &
      screened_rates(k_cr44__p_v43)*Y(jcr44) + &
      screened_rates(k_cr45__p_ti44__weak__wc12)*Y(jcr45) + &
      screened_rates(k_cr45__p_v44)*Y(jcr45) + screened_rates(k_cr46__p_v45)*Y(jcr46) &
      + screened_rates(k_cr47__p_v46)*Y(jcr47) + screened_rates(k_cr48__p_v47)* &
      Y(jcr48) + screened_rates(k_cr49__p_v48)*Y(jcr49) + &
      screened_rates(k_cr50__p_v49)*Y(jcr50) + screened_rates(k_cr51__p_v50)*Y(jcr51) &
      + screened_rates(k_d_be7__p_he4_he4)*Y(jbe7)*Y(jd)*state % rho + &
      screened_rates(k_d_he3__p_he4)*Y(jd)*Y(jhe3)*state % rho + &
      screened_rates(k_f16__p_o15)*Y(jf16) + screened_rates(k_f17__p_o16)*Y(jf17) + &
      screened_rates(k_f18__p_o17)*Y(jf18) + screened_rates(k_f19__p_o18)*Y(jf19) + &
      screened_rates(k_fe45__p_cr44__weak__wc17)*Y(jfe45) + &
      screened_rates(k_fe45__p_mn44)*Y(jfe45) + 3.0d0* &
      screened_rates(k_fe45__p_p_p_ti42__weak__wc12)*Y(jfe45) + 2.0d0* &
      screened_rates(k_fe45__p_p_v43__weak__wc12)*Y(jfe45) + &
      screened_rates(k_fe46__p_cr45__weak__wc12)*Y(jfe46) + &
      screened_rates(k_fe46__p_mn45)*Y(jfe46) + &
      screened_rates(k_fe47__p_cr46__weak__wc12)*Y(jfe47) + &
      screened_rates(k_fe47__p_mn46)*Y(jfe47) + &
      screened_rates(k_fe48__p_cr47__weak__wc12)*Y(jfe48) + &
      screened_rates(k_fe48__p_mn47)*Y(jfe48) + &
      screened_rates(k_fe49__p_cr48__weak__wc12)*Y(jfe49) + &
      screened_rates(k_fe49__p_mn48)*Y(jfe49) + screened_rates(k_fe50__p_mn49)* &
      Y(jfe50) + screened_rates(k_fe51__p_mn50)*Y(jfe51) + &
      screened_rates(k_fe52__p_mn51)*Y(jfe52) + screened_rates(k_fe53__p_mn52)* &
      Y(jfe53) + screened_rates(k_fe54__p_mn53)*Y(jfe54) + &
      screened_rates(k_fe55__p_mn54)*Y(jfe55) + screened_rates(k_fe56__p_mn55)* &
      Y(jfe56) + screened_rates(k_he3__p_d)*Y(jhe3) + 2.0d0* &
      screened_rates(k_he3_be7__p_p_he4_he4)*Y(jbe7)*Y(jhe3)*state % rho + &
      screened_rates(k_he3_he3__p_p_he4)*Y(jhe3)**2*state % rho + &
      screened_rates(k_he4_al22__p_si25)*Y(jal22)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_al23__p_si26)*Y(jal23)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_al24__p_si27)*Y(jal24)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_al25__p_si28)*Y(jal25)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_al26__p_si29)*Y(jal26)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_al27__p_si30)*Y(jal27)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ar31__p_k34)*Y(jar31)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ar32__p_k35)*Y(jar32)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ar33__p_k36)*Y(jar33)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ar34__p_k37)*Y(jar34)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ar35__p_k38)*Y(jar35)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ar36__p_k39)*Y(jar36)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ar37__p_k40)*Y(jar37)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ar38__p_k41)*Y(jar38)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_c12__p_n15)*Y(jc12)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ca34__p_sc37)*Y(jca34)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ca35__p_sc38)*Y(jca35)*Y(jhe4)*state % rho + 2.0d0* &
      screened_rates(k_he4_ca36__p_p_ca38)*Y(jca36)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ca36__p_sc39)*Y(jca36)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ca37__p_sc40)*Y(jca37)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ca38__p_sc41)*Y(jca38)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ca39__p_sc42)*Y(jca39)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ca40__p_sc43)*Y(jca40)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ca41__p_sc44)*Y(jca41)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ca42__p_sc45)*Y(jca42)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ca43__p_sc46)*Y(jca43)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cl29__p_ar32)*Y(jcl29)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cl30__p_ar33)*Y(jcl30)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cl31__p_ar34)*Y(jcl31)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cl32__p_ar35)*Y(jcl32)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cl33__p_ar36)*Y(jcl33)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cl34__p_ar37)*Y(jcl34)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cl35__p_ar38)*Y(jcl35)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cl36__p_ar39)*Y(jcl36)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_co47__p_ni50)*Y(jco47)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_co48__p_ni51)*Y(jco48)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_co49__p_ni52)*Y(jco49)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_co50__p_ni53)*Y(jco50)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_co51__p_ni54)*Y(jco51)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_co52__p_ni55)*Y(jco52)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_co53__p_ni56)*Y(jco53)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cr42__p_mn45)*Y(jcr42)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cr43__p_mn46)*Y(jcr43)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cr44__p_mn47)*Y(jcr44)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cr45__p_mn48)*Y(jcr45)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cr46__p_mn49)*Y(jcr46)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cr47__p_mn50)*Y(jcr47)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cr48__p_mn51)*Y(jcr48)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cr49__p_mn52)*Y(jcr49)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cr50__p_mn53)*Y(jcr50)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cr51__p_mn54)*Y(jcr51)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cr52__p_mn55)*Y(jcr52)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_f16__p_ne19)*Y(jf16)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_f17__p_ne20)*Y(jf17)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_f18__p_ne21)*Y(jf18)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_f19__p_ne22)*Y(jf19)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe45__p_co48)*Y(jfe45)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe46__p_co49)*Y(jfe46)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe47__p_co50)*Y(jfe47)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe48__p_co51)*Y(jfe48)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe49__p_co52)*Y(jfe49)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe50__p_co53)*Y(jfe50)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe51__p_co54)*Y(jfe51)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe52__p_co55)*Y(jfe52)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe53__p_co56)*Y(jfe53)*Y(jhe4)*state % rho + 0.5d0* &
      screened_rates(k_he4_he4__p_li7)*Y(jhe4)**2*state % rho + &
      screened_rates(k_he4_k33__p_ca36)*Y(jhe4)*Y(jk33)*state % rho + &
      screened_rates(k_he4_k34__p_ca37)*Y(jhe4)*Y(jk34)*state % rho + &
      screened_rates(k_he4_k35__p_ca38)*Y(jhe4)*Y(jk35)*state % rho + &
      screened_rates(k_he4_k36__p_ca39)*Y(jhe4)*Y(jk36)*state % rho + &
      screened_rates(k_he4_k37__p_ca40)*Y(jhe4)*Y(jk37)*state % rho + &
      screened_rates(k_he4_k38__p_ca41)*Y(jhe4)*Y(jk38)*state % rho + &
      screened_rates(k_he4_k39__p_ca42)*Y(jhe4)*Y(jk39)*state % rho + &
      screened_rates(k_he4_k40__p_ca43)*Y(jhe4)*Y(jk40)*state % rho + &
      screened_rates(k_he4_k41__p_ca44)*Y(jhe4)*Y(jk41)*state % rho + &
      screened_rates(k_he4_mg21__p_al24)*Y(jhe4)*Y(jmg21)*state % rho + &
      screened_rates(k_he4_mg22__p_al25)*Y(jhe4)*Y(jmg22)*state % rho + &
      screened_rates(k_he4_mg23__p_al26)*Y(jhe4)*Y(jmg23)*state % rho + &
      screened_rates(k_he4_mg24__p_al27)*Y(jhe4)*Y(jmg24)*state % rho + &
      screened_rates(k_he4_mg25__p_al28)*Y(jhe4)*Y(jmg25)*state % rho + &
      screened_rates(k_he4_mn44__p_fe47)*Y(jhe4)*Y(jmn44)*state % rho + &
      screened_rates(k_he4_mn45__p_fe48)*Y(jhe4)*Y(jmn45)*state % rho + &
      screened_rates(k_he4_mn46__p_fe49)*Y(jhe4)*Y(jmn46)*state % rho + &
      screened_rates(k_he4_mn47__p_fe50)*Y(jhe4)*Y(jmn47)*state % rho + &
      screened_rates(k_he4_mn48__p_fe51)*Y(jhe4)*Y(jmn48)*state % rho + &
      screened_rates(k_he4_mn49__p_fe52)*Y(jhe4)*Y(jmn49)*state % rho + &
      screened_rates(k_he4_mn50__p_fe53)*Y(jhe4)*Y(jmn50)*state % rho + &
      screened_rates(k_he4_mn51__p_fe54)*Y(jhe4)*Y(jmn51)*state % rho + &
      screened_rates(k_he4_mn52__p_fe55)*Y(jhe4)*Y(jmn52)*state % rho + &
      screened_rates(k_he4_mn53__p_fe56)*Y(jhe4)*Y(jmn53)*state % rho + &
      screened_rates(k_he4_n12__p_o15)*Y(jhe4)*Y(jn12)*state % rho + &
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*Y(jn13)*state % rho + &
      screened_rates(k_he4_n14__p_o17)*Y(jhe4)*Y(jn14)*state % rho + &
      screened_rates(k_he4_n15__p_o18)*Y(jhe4)*Y(jn15)*state % rho + &
      screened_rates(k_he4_na19__p_mg22)*Y(jhe4)*Y(jna19)*state % rho + &
      screened_rates(k_he4_na20__p_mg23)*Y(jhe4)*Y(jna20)*state % rho + &
      screened_rates(k_he4_na21__p_mg24)*Y(jhe4)*Y(jna21)*state % rho + &
      screened_rates(k_he4_na22__p_mg25)*Y(jhe4)*Y(jna22)*state % rho + &
      screened_rates(k_he4_na23__p_mg26)*Y(jhe4)*Y(jna23)*state % rho + &
      screened_rates(k_he4_ne17__p_na20)*Y(jhe4)*Y(jne17)*state % rho + &
      screened_rates(k_he4_ne18__p_na21)*Y(jhe4)*Y(jne18)*state % rho + &
      screened_rates(k_he4_ne19__p_na22)*Y(jhe4)*Y(jne19)*state % rho + &
      screened_rates(k_he4_ne20__p_na23)*Y(jhe4)*Y(jne20)*state % rho + &
      screened_rates(k_he4_ne21__p_na24)*Y(jhe4)*Y(jne21)*state % rho + &
      screened_rates(k_he4_o14__p_f17)*Y(jhe4)*Y(jo14)*state % rho + &
      screened_rates(k_he4_o15__p_f18)*Y(jhe4)*Y(jo15)*state % rho + &
      screened_rates(k_he4_o16__p_f19)*Y(jhe4)*Y(jo16)*state % rho + &
      screened_rates(k_he4_o17__p_f20)*Y(jhe4)*Y(jo17)*state % rho + &
      screened_rates(k_he4_p26__p_s29)*Y(jhe4)*Y(jp26)*state % rho + &
      screened_rates(k_he4_p27__p_s30)*Y(jhe4)*Y(jp27)*state % rho + &
      screened_rates(k_he4_p28__p_s31)*Y(jhe4)*Y(jp28)*state % rho + &
      screened_rates(k_he4_p29__p_s32)*Y(jhe4)*Y(jp29)*state % rho + &
      screened_rates(k_he4_p30__p_s33)*Y(jhe4)*Y(jp30)*state % rho + &
      screened_rates(k_he4_p31__p_s34)*Y(jhe4)*Y(jp31)*state % rho + &
      screened_rates(k_he4_p32__p_s35)*Y(jhe4)*Y(jp32)*state % rho + &
      screened_rates(k_he4_s28__p_cl31)*Y(jhe4)*Y(js28)*state % rho + &
      screened_rates(k_he4_s29__p_cl32)*Y(jhe4)*Y(js29)*state % rho + &
      screened_rates(k_he4_s30__p_cl33)*Y(jhe4)*Y(js30)*state % rho + &
      screened_rates(k_he4_s31__p_cl34)*Y(jhe4)*Y(js31)*state % rho + &
      screened_rates(k_he4_s32__p_cl35)*Y(jhe4)*Y(js32)*state % rho + &
      screened_rates(k_he4_s33__p_cl36)*Y(jhe4)*Y(js33)*state % rho + &
      screened_rates(k_he4_s34__p_cl37)*Y(jhe4)*Y(js34)*state % rho + &
      screened_rates(k_he4_sc36__p_ti39)*Y(jhe4)*Y(jsc36)*state % rho + &
      screened_rates(k_he4_sc37__p_ti40)*Y(jhe4)*Y(jsc37)*state % rho + &
      screened_rates(k_he4_sc38__p_ti41)*Y(jhe4)*Y(jsc38)*state % rho + &
      screened_rates(k_he4_sc39__p_ti42)*Y(jhe4)*Y(jsc39)*state % rho + &
      screened_rates(k_he4_sc40__p_ti43)*Y(jhe4)*Y(jsc40)*state % rho + &
      screened_rates(k_he4_sc41__p_ti44)*Y(jhe4)*Y(jsc41)*state % rho + &
      screened_rates(k_he4_sc42__p_ti45)*Y(jhe4)*Y(jsc42)*state % rho + &
      screened_rates(k_he4_sc43__p_ti46)*Y(jhe4)*Y(jsc43)*state % rho + &
      screened_rates(k_he4_sc44__p_ti47)*Y(jhe4)*Y(jsc44)*state % rho + &
      screened_rates(k_he4_sc45__p_ti48)*Y(jhe4)*Y(jsc45)*state % rho + &
      screened_rates(k_he4_si24__p_p27)*Y(jhe4)*Y(jsi24)*state % rho + &
      screened_rates(k_he4_si25__p_p28)*Y(jhe4)*Y(jsi25)*state % rho + &
      screened_rates(k_he4_si26__p_p29)*Y(jhe4)*Y(jsi26)*state % rho + &
      screened_rates(k_he4_si27__p_p30)*Y(jhe4)*Y(jsi27)*state % rho + &
      screened_rates(k_he4_si28__p_p31)*Y(jhe4)*Y(jsi28)*state % rho + &
      screened_rates(k_he4_si29__p_p32)*Y(jhe4)*Y(jsi29)*state % rho + &
      screened_rates(k_he4_si30__p_p33)*Y(jhe4)*Y(jsi30)*state % rho + &
      screened_rates(k_he4_ti38__p_v41)*Y(jhe4)*Y(jti38)*state % rho + &
      screened_rates(k_he4_ti39__p_v42)*Y(jhe4)*Y(jti39)*state % rho + &
      screened_rates(k_he4_ti40__p_v43)*Y(jhe4)*Y(jti40)*state % rho + &
      screened_rates(k_he4_ti41__p_v44)*Y(jhe4)*Y(jti41)*state % rho + &
      screened_rates(k_he4_ti42__p_v45)*Y(jhe4)*Y(jti42)*state % rho + &
      screened_rates(k_he4_ti43__p_v46)*Y(jhe4)*Y(jti43)*state % rho + &
      screened_rates(k_he4_ti44__p_v47)*Y(jhe4)*Y(jti44)*state % rho + &
      screened_rates(k_he4_ti45__p_v48)*Y(jhe4)*Y(jti45)*state % rho + &
      screened_rates(k_he4_ti46__p_v49)*Y(jhe4)*Y(jti46)*state % rho + &
      screened_rates(k_he4_ti47__p_v50)*Y(jhe4)*Y(jti47)*state % rho + &
      screened_rates(k_he4_v40__p_cr43)*Y(jhe4)*Y(jv40)*state % rho + &
      screened_rates(k_he4_v41__p_cr44)*Y(jhe4)*Y(jv41)*state % rho + &
      screened_rates(k_he4_v42__p_cr45)*Y(jhe4)*Y(jv42)*state % rho + &
      screened_rates(k_he4_v43__p_cr46)*Y(jhe4)*Y(jv43)*state % rho + &
      screened_rates(k_he4_v44__p_cr47)*Y(jhe4)*Y(jv44)*state % rho + &
      screened_rates(k_he4_v45__p_cr48)*Y(jhe4)*Y(jv45)*state % rho + &
      screened_rates(k_he4_v46__p_cr49)*Y(jhe4)*Y(jv46)*state % rho + &
      screened_rates(k_he4_v47__p_cr50)*Y(jhe4)*Y(jv47)*state % rho + &
      screened_rates(k_he4_v48__p_cr51)*Y(jhe4)*Y(jv48)*state % rho + &
      screened_rates(k_he4_v49__p_cr52)*Y(jhe4)*Y(jv49)*state % rho + &
      screened_rates(k_k33__p_ar32)*Y(jk33) + screened_rates(k_k34__p_ar33)*Y(jk34) + &
      screened_rates(k_k35__p_ar34)*Y(jk35) + screened_rates(k_k35__p_cl34__weak__wc12) &
      *Y(jk35) + screened_rates(k_k36__p_ar35)*Y(jk36) + &
      screened_rates(k_k36__p_cl35__weak__wc12)*Y(jk36) + screened_rates(k_k37__p_ar36) &
      *Y(jk37) + screened_rates(k_k38__p_ar37)*Y(jk38) + &
      screened_rates(k_k39__p_ar38)*Y(jk39) + screened_rates(k_k40__p_ar39)*Y(jk40) + &
      screened_rates(k_mg21__p_na20)*Y(jmg21) + &
      screened_rates(k_mg21__p_ne20__weak__wc12)*Y(jmg21) + &
      screened_rates(k_mg22__p_na21)*Y(jmg22) + screened_rates(k_mg23__p_na22)* &
      Y(jmg23) + screened_rates(k_mg24__p_na23)*Y(jmg24) + &
      screened_rates(k_mg25__p_na24)*Y(jmg25) + screened_rates(k_mn44__p_cr43)* &
      Y(jmn44) + screened_rates(k_mn45__p_cr44)*Y(jmn45) + &
      screened_rates(k_mn46__p_cr45)*Y(jmn46) + &
      screened_rates(k_mn46__p_v45__weak__wc12)*Y(jmn46) + &
      screened_rates(k_mn47__p_cr46)*Y(jmn47) + screened_rates(k_mn48__p_cr47)* &
      Y(jmn48) + screened_rates(k_mn48__p_v47__weak__wc12)*Y(jmn48) + &
      screened_rates(k_mn49__p_cr48)*Y(jmn49) + screened_rates(k_mn50__p_cr49)* &
      Y(jmn50) + screened_rates(k_mn51__p_cr50)*Y(jmn51) + &
      screened_rates(k_mn52__p_cr51)*Y(jmn52) + screened_rates(k_mn53__p_cr52)* &
      Y(jmn53) + screened_rates(k_n13__p_c12)*Y(jn13) + screened_rates(k_n14__p_c13)* &
      Y(jn14) + screened_rates(k_na19__p_ne18)*Y(jna19) + &
      screened_rates(k_na20__p_ne19)*Y(jna20) + screened_rates(k_na21__p_ne20)* &
      Y(jna21) + screened_rates(k_na22__p_ne21)*Y(jna22) + &
      screened_rates(k_na23__p_ne22)*Y(jna23) + &
      screened_rates(k_ne17__p_o16__weak__wc12)*Y(jne17) + &
      screened_rates(k_ne18__p_f17)*Y(jne18) + screened_rates(k_ne19__p_f18)*Y(jne19) &
      + screened_rates(k_ne20__p_f19)*Y(jne20) + screened_rates(k_ne21__p_f20)* &
      Y(jne21) + screened_rates(k_ni48__p_co47)*Y(jni48) + &
      screened_rates(k_ni49__p_co48)*Y(jni49) + &
      screened_rates(k_ni49__p_fe48__weak__wc12)*Y(jni49) + &
      screened_rates(k_ni50__p_co49)*Y(jni50) + &
      screened_rates(k_ni50__p_fe49__weak__wc12)*Y(jni50) + &
      screened_rates(k_ni51__p_co50)*Y(jni51) + &
      screened_rates(k_ni51__p_fe50__weak__wc12)*Y(jni51) + &
      screened_rates(k_ni52__p_co51)*Y(jni52) + &
      screened_rates(k_ni52__p_fe51__weak__wc12)*Y(jni52) + &
      screened_rates(k_ni53__p_co52)*Y(jni53) + &
      screened_rates(k_ni53__p_fe52__weak__wc12)*Y(jni53) + &
      screened_rates(k_ni54__p_co53)*Y(jni54) + screened_rates(k_ni55__p_co54)* &
      Y(jni55) + screened_rates(k_ni56__p_co55)*Y(jni56) + &
      screened_rates(k_o14__p_n13)*Y(jo14) + screened_rates(k_o15__p_n14)*Y(jo15) + &
      screened_rates(k_o16__p_n15)*Y(jo16) + screened_rates(k_p26__p_si25)*Y(jp26) + &
      screened_rates(k_p27__p_al26__weak__wc12)*Y(jp27) + screened_rates(k_p27__p_si26) &
      *Y(jp27) + screened_rates(k_p28__p_al27__weak__wc12)*Y(jp28) + &
      screened_rates(k_p28__p_si27)*Y(jp28) + screened_rates(k_p29__p_si28)*Y(jp29) + &
      screened_rates(k_p30__p_si29)*Y(jp30) + screened_rates(k_p31__p_si30)*Y(jp31) - &
      screened_rates(k_p_al23__si24)*Y(jal23)*Y(jp)*state % rho - &
      screened_rates(k_p_al24__he4_mg21)*Y(jal24)*Y(jp)*state % rho - &
      screened_rates(k_p_al24__si25)*Y(jal24)*Y(jp)*state % rho - &
      screened_rates(k_p_al25__he4_mg22)*Y(jal25)*Y(jp)*state % rho - &
      screened_rates(k_p_al25__si26)*Y(jal25)*Y(jp)*state % rho - &
      screened_rates(k_p_al26__he4_mg23)*Y(jal26)*Y(jp)*state % rho - &
      screened_rates(k_p_al26__si27)*Y(jal26)*Y(jp)*state % rho - &
      screened_rates(k_p_al27__he4_mg24)*Y(jal27)*Y(jp)*state % rho - &
      screened_rates(k_p_al27__si28)*Y(jal27)*Y(jp)*state % rho - &
      screened_rates(k_p_al28__he4_mg25)*Y(jal28)*Y(jp)*state % rho - &
      screened_rates(k_p_al28__si29)*Y(jal28)*Y(jp)*state % rho - &
      screened_rates(k_p_ar32__he4_cl29)*Y(jar32)*Y(jp)*state % rho - &
      screened_rates(k_p_ar32__k33)*Y(jar32)*Y(jp)*state % rho - &
      screened_rates(k_p_ar33__he4_cl30)*Y(jar33)*Y(jp)*state % rho - &
      screened_rates(k_p_ar33__k34)*Y(jar33)*Y(jp)*state % rho - &
      screened_rates(k_p_ar34__he4_cl31)*Y(jar34)*Y(jp)*state % rho - &
      screened_rates(k_p_ar34__k35)*Y(jar34)*Y(jp)*state % rho - &
      screened_rates(k_p_ar35__he4_cl32)*Y(jar35)*Y(jp)*state % rho - &
      screened_rates(k_p_ar35__k36)*Y(jar35)*Y(jp)*state % rho - &
      screened_rates(k_p_ar36__he4_cl33)*Y(jar36)*Y(jp)*state % rho - &
      screened_rates(k_p_ar36__k37)*Y(jar36)*Y(jp)*state % rho - &
      screened_rates(k_p_ar37__he4_cl34)*Y(jar37)*Y(jp)*state % rho - &
      screened_rates(k_p_ar37__k38)*Y(jar37)*Y(jp)*state % rho - &
      screened_rates(k_p_ar38__he4_cl35)*Y(jar38)*Y(jp)*state % rho - &
      screened_rates(k_p_ar38__k39)*Y(jar38)*Y(jp)*state % rho - &
      screened_rates(k_p_ar39__he4_cl36)*Y(jar39)*Y(jp)*state % rho - &
      screened_rates(k_p_ar39__k40)*Y(jar39)*Y(jp)*state % rho - &
      screened_rates(k_p_be7__b8)*Y(jbe7)*Y(jp)*state % rho - &
      screened_rates(k_p_c12__n13)*Y(jc12)*Y(jp)*state % rho - &
      screened_rates(k_p_c13__n14)*Y(jc13)*Y(jp)*state % rho - &
      screened_rates(k_p_ca35__sc36)*Y(jca35)*Y(jp)*state % rho - &
      screened_rates(k_p_ca36__he4_k33)*Y(jca36)*Y(jp)*state % rho - &
      screened_rates(k_p_ca36__sc37)*Y(jca36)*Y(jp)*state % rho - &
      screened_rates(k_p_ca37__he4_k34)*Y(jca37)*Y(jp)*state % rho - &
      screened_rates(k_p_ca37__sc38)*Y(jca37)*Y(jp)*state % rho - &
      screened_rates(k_p_ca38__he4_k35)*Y(jca38)*Y(jp)*state % rho - &
      screened_rates(k_p_ca38__sc39)*Y(jca38)*Y(jp)*state % rho - &
      screened_rates(k_p_ca39__he4_k36)*Y(jca39)*Y(jp)*state % rho - &
      screened_rates(k_p_ca39__sc40)*Y(jca39)*Y(jp)*state % rho - &
      screened_rates(k_p_ca40__he4_k37)*Y(jca40)*Y(jp)*state % rho - &
      screened_rates(k_p_ca40__sc41)*Y(jca40)*Y(jp)*state % rho - &
      screened_rates(k_p_ca41__he4_k38)*Y(jca41)*Y(jp)*state % rho - &
      screened_rates(k_p_ca41__sc42)*Y(jca41)*Y(jp)*state % rho - &
      screened_rates(k_p_ca42__he4_k39)*Y(jca42)*Y(jp)*state % rho - &
      screened_rates(k_p_ca42__sc43)*Y(jca42)*Y(jp)*state % rho - &
      screened_rates(k_p_ca43__he4_k40)*Y(jca43)*Y(jp)*state % rho - &
      screened_rates(k_p_ca43__sc44)*Y(jca43)*Y(jp)*state % rho - &
      screened_rates(k_p_ca44__he4_k41)*Y(jca44)*Y(jp)*state % rho - &
      screened_rates(k_p_ca44__sc45)*Y(jca44)*Y(jp)*state % rho - &
      screened_rates(k_p_cl30__ar31)*Y(jcl30)*Y(jp)*state % rho - &
      screened_rates(k_p_cl31__ar32)*Y(jcl31)*Y(jp)*state % rho - &
      screened_rates(k_p_cl31__he4_s28)*Y(jcl31)*Y(jp)*state % rho - &
      screened_rates(k_p_cl32__ar33)*Y(jcl32)*Y(jp)*state % rho - &
      screened_rates(k_p_cl32__he4_s29)*Y(jcl32)*Y(jp)*state % rho - &
      screened_rates(k_p_cl33__ar34)*Y(jcl33)*Y(jp)*state % rho - &
      screened_rates(k_p_cl33__he4_s30)*Y(jcl33)*Y(jp)*state % rho - &
      screened_rates(k_p_cl34__ar35)*Y(jcl34)*Y(jp)*state % rho - &
      screened_rates(k_p_cl34__he4_s31)*Y(jcl34)*Y(jp)*state % rho - &
      screened_rates(k_p_cl35__ar36)*Y(jcl35)*Y(jp)*state % rho - &
      screened_rates(k_p_cl35__he4_s32)*Y(jcl35)*Y(jp)*state % rho - &
      screened_rates(k_p_cl36__ar37)*Y(jcl36)*Y(jp)*state % rho - &
      screened_rates(k_p_cl36__he4_s33)*Y(jcl36)*Y(jp)*state % rho - &
      screened_rates(k_p_cl37__ar38)*Y(jcl37)*Y(jp)*state % rho - &
      screened_rates(k_p_cl37__he4_s34)*Y(jcl37)*Y(jp)*state % rho - &
      screened_rates(k_p_co47__ni48)*Y(jco47)*Y(jp)*state % rho - &
      screened_rates(k_p_co48__he4_fe45)*Y(jco48)*Y(jp)*state % rho - &
      screened_rates(k_p_co48__ni49)*Y(jco48)*Y(jp)*state % rho - &
      screened_rates(k_p_co49__he4_fe46)*Y(jco49)*Y(jp)*state % rho - &
      screened_rates(k_p_co49__ni50)*Y(jco49)*Y(jp)*state % rho - &
      screened_rates(k_p_co50__he4_fe47)*Y(jco50)*Y(jp)*state % rho - &
      screened_rates(k_p_co50__ni51)*Y(jco50)*Y(jp)*state % rho - &
      screened_rates(k_p_co51__he4_fe48)*Y(jco51)*Y(jp)*state % rho - &
      screened_rates(k_p_co51__ni52)*Y(jco51)*Y(jp)*state % rho - &
      screened_rates(k_p_co52__he4_fe49)*Y(jco52)*Y(jp)*state % rho - &
      screened_rates(k_p_co52__ni53)*Y(jco52)*Y(jp)*state % rho - &
      screened_rates(k_p_co53__he4_fe50)*Y(jco53)*Y(jp)*state % rho - &
      screened_rates(k_p_co53__ni54)*Y(jco53)*Y(jp)*state % rho - &
      screened_rates(k_p_co54__he4_fe51)*Y(jco54)*Y(jp)*state % rho - &
      screened_rates(k_p_co54__ni55)*Y(jco54)*Y(jp)*state % rho - &
      screened_rates(k_p_co55__he4_fe52)*Y(jco55)*Y(jp)*state % rho - &
      screened_rates(k_p_co55__ni56)*Y(jco55)*Y(jp)*state % rho - &
      screened_rates(k_p_co56__he4_fe53)*Y(jco56)*Y(jp)*state % rho - &
      screened_rates(k_p_cr43__he4_v40)*Y(jcr43)*Y(jp)*state % rho - &
      screened_rates(k_p_cr43__mn44)*Y(jcr43)*Y(jp)*state % rho - &
      screened_rates(k_p_cr44__he4_v41)*Y(jcr44)*Y(jp)*state % rho - &
      screened_rates(k_p_cr44__mn45)*Y(jcr44)*Y(jp)*state % rho - &
      screened_rates(k_p_cr45__he4_v42)*Y(jcr45)*Y(jp)*state % rho - &
      screened_rates(k_p_cr45__mn46)*Y(jcr45)*Y(jp)*state % rho - &
      screened_rates(k_p_cr46__he4_v43)*Y(jcr46)*Y(jp)*state % rho - &
      screened_rates(k_p_cr46__mn47)*Y(jcr46)*Y(jp)*state % rho - &
      screened_rates(k_p_cr47__he4_v44)*Y(jcr47)*Y(jp)*state % rho - &
      screened_rates(k_p_cr47__mn48)*Y(jcr47)*Y(jp)*state % rho - &
      screened_rates(k_p_cr48__he4_v45)*Y(jcr48)*Y(jp)*state % rho - &
      screened_rates(k_p_cr48__mn49)*Y(jcr48)*Y(jp)*state % rho - &
      screened_rates(k_p_cr49__he4_v46)*Y(jcr49)*Y(jp)*state % rho - &
      screened_rates(k_p_cr49__mn50)*Y(jcr49)*Y(jp)*state % rho - &
      screened_rates(k_p_cr50__he4_v47)*Y(jcr50)*Y(jp)*state % rho - &
      screened_rates(k_p_cr50__mn51)*Y(jcr50)*Y(jp)*state % rho - &
      screened_rates(k_p_cr51__he4_v48)*Y(jcr51)*Y(jp)*state % rho - &
      screened_rates(k_p_cr51__mn52)*Y(jcr51)*Y(jp)*state % rho - &
      screened_rates(k_p_cr52__he4_v49)*Y(jcr52)*Y(jp)*state % rho - &
      screened_rates(k_p_cr52__mn53)*Y(jcr52)*Y(jp)*state % rho - &
      screened_rates(k_p_d__he3)*Y(jd)*Y(jp)*state % rho - &
      screened_rates(k_p_f17__he4_o14)*Y(jf17)*Y(jp)*state % rho - &
      screened_rates(k_p_f17__ne18)*Y(jf17)*Y(jp)*state % rho - &
      screened_rates(k_p_f18__he4_o15)*Y(jf18)*Y(jp)*state % rho - &
      screened_rates(k_p_f18__ne19)*Y(jf18)*Y(jp)*state % rho - &
      screened_rates(k_p_f19__he4_o16)*Y(jf19)*Y(jp)*state % rho - &
      screened_rates(k_p_f19__ne20)*Y(jf19)*Y(jp)*state % rho - &
      screened_rates(k_p_f20__he4_o17)*Y(jf20)*Y(jp)*state % rho - &
      screened_rates(k_p_f20__ne21)*Y(jf20)*Y(jp)*state % rho - &
      screened_rates(k_p_fe46__co47)*Y(jfe46)*Y(jp)*state % rho - &
      screened_rates(k_p_fe47__co48)*Y(jfe47)*Y(jp)*state % rho - &
      screened_rates(k_p_fe47__he4_mn44)*Y(jfe47)*Y(jp)*state % rho - &
      screened_rates(k_p_fe48__co49)*Y(jfe48)*Y(jp)*state % rho - &
      screened_rates(k_p_fe48__he4_mn45)*Y(jfe48)*Y(jp)*state % rho - &
      screened_rates(k_p_fe49__co50)*Y(jfe49)*Y(jp)*state % rho - &
      screened_rates(k_p_fe49__he4_mn46)*Y(jfe49)*Y(jp)*state % rho - &
      screened_rates(k_p_fe50__co51)*Y(jfe50)*Y(jp)*state % rho - &
      screened_rates(k_p_fe50__he4_mn47)*Y(jfe50)*Y(jp)*state % rho - &
      screened_rates(k_p_fe51__co52)*Y(jfe51)*Y(jp)*state % rho - &
      screened_rates(k_p_fe51__he4_mn48)*Y(jfe51)*Y(jp)*state % rho - &
      screened_rates(k_p_fe52__co53)*Y(jfe52)*Y(jp)*state % rho - &
      screened_rates(k_p_fe52__he4_mn49)*Y(jfe52)*Y(jp)*state % rho - &
      screened_rates(k_p_fe53__co54)*Y(jfe53)*Y(jp)*state % rho - &
      screened_rates(k_p_fe53__he4_mn50)*Y(jfe53)*Y(jp)*state % rho - &
      screened_rates(k_p_fe54__co55)*Y(jfe54)*Y(jp)*state % rho - &
      screened_rates(k_p_fe54__he4_mn51)*Y(jfe54)*Y(jp)*state % rho - &
      screened_rates(k_p_fe55__co56)*Y(jfe55)*Y(jp)*state % rho - &
      screened_rates(k_p_fe55__he4_mn52)*Y(jfe55)*Y(jp)*state % rho - &
      screened_rates(k_p_fe56__he4_mn53)*Y(jfe56)*Y(jp)*state % rho - &
      screened_rates(k_p_he3__he4__weak__bet_pos_)*Y(jhe3)*Y(jp)*state % rho - &
      screened_rates(k_p_he4__d_he3)*Y(jhe4)*Y(jp)*state % rho - 0.5d0* &
      screened_rates(k_p_he4_he4__d_be7)*Y(jhe4)**2*Y(jp)*state % rho**2 - &
      screened_rates(k_p_k33__ca34)*Y(jk33)*Y(jp)*state % rho - &
      screened_rates(k_p_k34__ca35)*Y(jk34)*Y(jp)*state % rho - &
      screened_rates(k_p_k34__he4_ar31)*Y(jk34)*Y(jp)*state % rho - &
      screened_rates(k_p_k35__ca36)*Y(jk35)*Y(jp)*state % rho - &
      screened_rates(k_p_k35__he4_ar32)*Y(jk35)*Y(jp)*state % rho - &
      screened_rates(k_p_k36__ca37)*Y(jk36)*Y(jp)*state % rho - &
      screened_rates(k_p_k36__he4_ar33)*Y(jk36)*Y(jp)*state % rho - &
      screened_rates(k_p_k37__ca38)*Y(jk37)*Y(jp)*state % rho - &
      screened_rates(k_p_k37__he4_ar34)*Y(jk37)*Y(jp)*state % rho - &
      screened_rates(k_p_k38__ca39)*Y(jk38)*Y(jp)*state % rho - &
      screened_rates(k_p_k38__he4_ar35)*Y(jk38)*Y(jp)*state % rho - &
      screened_rates(k_p_k39__ca40)*Y(jk39)*Y(jp)*state % rho - &
      screened_rates(k_p_k39__he4_ar36)*Y(jk39)*Y(jp)*state % rho - &
      screened_rates(k_p_k40__ca41)*Y(jk40)*Y(jp)*state % rho - &
      screened_rates(k_p_k40__he4_ar37)*Y(jk40)*Y(jp)*state % rho - &
      screened_rates(k_p_k41__ca42)*Y(jk41)*Y(jp)*state % rho - &
      screened_rates(k_p_k41__he4_ar38)*Y(jk41)*Y(jp)*state % rho - &
      screened_rates(k_p_li7__he4_he4)*Y(jli7)*Y(jp)*state % rho - &
      screened_rates(k_p_mg21__al22)*Y(jmg21)*Y(jp)*state % rho - &
      screened_rates(k_p_mg22__al23)*Y(jmg22)*Y(jp)*state % rho - &
      screened_rates(k_p_mg22__he4_na19)*Y(jmg22)*Y(jp)*state % rho - &
      screened_rates(k_p_mg23__al24)*Y(jmg23)*Y(jp)*state % rho - &
      screened_rates(k_p_mg23__he4_na20)*Y(jmg23)*Y(jp)*state % rho - &
      screened_rates(k_p_mg24__al25)*Y(jmg24)*Y(jp)*state % rho - &
      screened_rates(k_p_mg24__he4_na21)*Y(jmg24)*Y(jp)*state % rho - &
      screened_rates(k_p_mg25__al26)*Y(jmg25)*Y(jp)*state % rho - &
      screened_rates(k_p_mg25__he4_na22)*Y(jmg25)*Y(jp)*state % rho - &
      screened_rates(k_p_mg26__al27)*Y(jmg26)*Y(jp)*state % rho - &
      screened_rates(k_p_mg26__he4_na23)*Y(jmg26)*Y(jp)*state % rho - &
      screened_rates(k_p_mn44__fe45)*Y(jmn44)*Y(jp)*state % rho - &
      screened_rates(k_p_mn45__fe46)*Y(jmn45)*Y(jp)*state % rho - &
      screened_rates(k_p_mn45__he4_cr42)*Y(jmn45)*Y(jp)*state % rho - &
      screened_rates(k_p_mn46__fe47)*Y(jmn46)*Y(jp)*state % rho - &
      screened_rates(k_p_mn46__he4_cr43)*Y(jmn46)*Y(jp)*state % rho - &
      screened_rates(k_p_mn47__fe48)*Y(jmn47)*Y(jp)*state % rho - &
      screened_rates(k_p_mn47__he4_cr44)*Y(jmn47)*Y(jp)*state % rho - &
      screened_rates(k_p_mn48__fe49)*Y(jmn48)*Y(jp)*state % rho - &
      screened_rates(k_p_mn48__he4_cr45)*Y(jmn48)*Y(jp)*state % rho - &
      screened_rates(k_p_mn49__fe50)*Y(jmn49)*Y(jp)*state % rho - &
      screened_rates(k_p_mn49__he4_cr46)*Y(jmn49)*Y(jp)*state % rho - &
      screened_rates(k_p_mn50__fe51)*Y(jmn50)*Y(jp)*state % rho - &
      screened_rates(k_p_mn50__he4_cr47)*Y(jmn50)*Y(jp)*state % rho - &
      screened_rates(k_p_mn51__fe52)*Y(jmn51)*Y(jp)*state % rho - &
      screened_rates(k_p_mn51__he4_cr48)*Y(jmn51)*Y(jp)*state % rho - &
      screened_rates(k_p_mn52__fe53)*Y(jmn52)*Y(jp)*state % rho - &
      screened_rates(k_p_mn52__he4_cr49)*Y(jmn52)*Y(jp)*state % rho - &
      screened_rates(k_p_mn53__fe54)*Y(jmn53)*Y(jp)*state % rho - &
      screened_rates(k_p_mn53__he4_cr50)*Y(jmn53)*Y(jp)*state % rho - &
      screened_rates(k_p_mn54__fe55)*Y(jmn54)*Y(jp)*state % rho - &
      screened_rates(k_p_mn54__he4_cr51)*Y(jmn54)*Y(jp)*state % rho - &
      screened_rates(k_p_mn55__fe56)*Y(jmn55)*Y(jp)*state % rho - &
      screened_rates(k_p_mn55__he4_cr52)*Y(jmn55)*Y(jp)*state % rho - &
      screened_rates(k_p_n13__o14)*Y(jn13)*Y(jp)*state % rho - &
      screened_rates(k_p_n14__o15)*Y(jn14)*Y(jp)*state % rho - &
      screened_rates(k_p_n15__he4_c12)*Y(jn15)*Y(jp)*state % rho - &
      screened_rates(k_p_n15__o16)*Y(jn15)*Y(jp)*state % rho - &
      screened_rates(k_p_na20__he4_ne17)*Y(jna20)*Y(jp)*state % rho - &
      screened_rates(k_p_na20__mg21)*Y(jna20)*Y(jp)*state % rho - &
      screened_rates(k_p_na21__he4_ne18)*Y(jna21)*Y(jp)*state % rho - &
      screened_rates(k_p_na21__mg22)*Y(jna21)*Y(jp)*state % rho - &
      screened_rates(k_p_na22__he4_ne19)*Y(jna22)*Y(jp)*state % rho - &
      screened_rates(k_p_na22__mg23)*Y(jna22)*Y(jp)*state % rho - &
      screened_rates(k_p_na23__he4_ne20)*Y(jna23)*Y(jp)*state % rho - &
      screened_rates(k_p_na23__mg24)*Y(jna23)*Y(jp)*state % rho - &
      screened_rates(k_p_na24__he4_ne21)*Y(jna24)*Y(jp)*state % rho - &
      screened_rates(k_p_na24__mg25)*Y(jna24)*Y(jp)*state % rho - &
      screened_rates(k_p_ne18__na19)*Y(jne18)*Y(jp)*state % rho - &
      screened_rates(k_p_ne19__he4_f16)*Y(jne19)*Y(jp)*state % rho - &
      screened_rates(k_p_ne19__na20)*Y(jne19)*Y(jp)*state % rho - &
      screened_rates(k_p_ne20__he4_f17)*Y(jne20)*Y(jp)*state % rho - &
      screened_rates(k_p_ne20__na21)*Y(jne20)*Y(jp)*state % rho - &
      screened_rates(k_p_ne21__he4_f18)*Y(jne21)*Y(jp)*state % rho - &
      screened_rates(k_p_ne21__na22)*Y(jne21)*Y(jp)*state % rho - &
      screened_rates(k_p_ne22__he4_f19)*Y(jne22)*Y(jp)*state % rho - &
      screened_rates(k_p_ne22__na23)*Y(jne22)*Y(jp)*state % rho - &
      screened_rates(k_p_ni50__he4_co47)*Y(jni50)*Y(jp)*state % rho - &
      screened_rates(k_p_ni51__he4_co48)*Y(jni51)*Y(jp)*state % rho - &
      screened_rates(k_p_ni52__he4_co49)*Y(jni52)*Y(jp)*state % rho - &
      screened_rates(k_p_ni53__he4_co50)*Y(jni53)*Y(jp)*state % rho - &
      screened_rates(k_p_ni54__he4_co51)*Y(jni54)*Y(jp)*state % rho - &
      screened_rates(k_p_ni55__he4_co52)*Y(jni55)*Y(jp)*state % rho - &
      screened_rates(k_p_ni56__he4_co53)*Y(jni56)*Y(jp)*state % rho - &
      screened_rates(k_p_o15__f16)*Y(jo15)*Y(jp)*state % rho - &
      screened_rates(k_p_o15__he4_n12)*Y(jo15)*Y(jp)*state % rho - &
      screened_rates(k_p_o16__f17)*Y(jo16)*Y(jp)*state % rho - &
      screened_rates(k_p_o16__he4_n13)*Y(jo16)*Y(jp)*state % rho - &
      screened_rates(k_p_o17__f18)*Y(jo17)*Y(jp)*state % rho - &
      screened_rates(k_p_o17__he4_n14)*Y(jo17)*Y(jp)*state % rho - &
      screened_rates(k_p_o18__f19)*Y(jo18)*Y(jp)*state % rho - &
      screened_rates(k_p_o18__he4_n15)*Y(jo18)*Y(jp)*state % rho - &
      screened_rates(k_p_p27__he4_si24)*Y(jp27)*Y(jp)*state % rho - &
      screened_rates(k_p_p27__s28)*Y(jp27)*Y(jp)*state % rho - &
      screened_rates(k_p_p28__he4_si25)*Y(jp28)*Y(jp)*state % rho - &
      screened_rates(k_p_p28__s29)*Y(jp28)*Y(jp)*state % rho - &
      screened_rates(k_p_p29__he4_si26)*Y(jp29)*Y(jp)*state % rho - &
      screened_rates(k_p_p29__s30)*Y(jp29)*Y(jp)*state % rho - &
      screened_rates(k_p_p30__he4_si27)*Y(jp30)*Y(jp)*state % rho - &
      screened_rates(k_p_p30__s31)*Y(jp30)*Y(jp)*state % rho - &
      screened_rates(k_p_p31__he4_si28)*Y(jp31)*Y(jp)*state % rho - &
      screened_rates(k_p_p31__s32)*Y(jp31)*Y(jp)*state % rho - &
      screened_rates(k_p_p32__he4_si29)*Y(jp32)*Y(jp)*state % rho - &
      screened_rates(k_p_p32__s33)*Y(jp32)*Y(jp)*state % rho - &
      screened_rates(k_p_p33__he4_si30)*Y(jp33)*Y(jp)*state % rho - &
      screened_rates(k_p_p33__s34)*Y(jp33)*Y(jp)*state % rho - &
      screened_rates(k_p_p__d__weak__bet_pos_)*Y(jp)**2*state % rho - &
      screened_rates(k_p_p__d__weak__electron_capture)*Y(jp)**2*state % rho**2* &
      state % y_e - screened_rates(k_p_p_ca38__he4_ca36)*Y(jca38)*Y(jp)**2* &
      state % rho**2 - screened_rates(k_p_p_he4__he3_he3)*Y(jhe4)*Y(jp)**2* &
      state % rho**2 - 0.5d0*screened_rates(k_p_p_he4_he4__he3_be7)*Y(jhe4)**2* &
      Y(jp)**2*state % rho**3 - screened_rates(k_p_s28__cl29)*Y(jp)*Y(js28)* &
      state % rho - screened_rates(k_p_s29__cl30)*Y(jp)*Y(js29)*state % rho - &
      screened_rates(k_p_s29__he4_p26)*Y(jp)*Y(js29)*state % rho - &
      screened_rates(k_p_s30__cl31)*Y(jp)*Y(js30)*state % rho - &
      screened_rates(k_p_s30__he4_p27)*Y(jp)*Y(js30)*state % rho - &
      screened_rates(k_p_s31__cl32)*Y(jp)*Y(js31)*state % rho - &
      screened_rates(k_p_s31__he4_p28)*Y(jp)*Y(js31)*state % rho - &
      screened_rates(k_p_s32__cl33)*Y(jp)*Y(js32)*state % rho - &
      screened_rates(k_p_s32__he4_p29)*Y(jp)*Y(js32)*state % rho - &
      screened_rates(k_p_s33__cl34)*Y(jp)*Y(js33)*state % rho - &
      screened_rates(k_p_s33__he4_p30)*Y(jp)*Y(js33)*state % rho - &
      screened_rates(k_p_s34__cl35)*Y(jp)*Y(js34)*state % rho - &
      screened_rates(k_p_s34__he4_p31)*Y(jp)*Y(js34)*state % rho - &
      screened_rates(k_p_s35__cl36)*Y(jp)*Y(js35)*state % rho - &
      screened_rates(k_p_s35__he4_p32)*Y(jp)*Y(js35)*state % rho - &
      screened_rates(k_p_sc37__he4_ca34)*Y(jp)*Y(jsc37)*state % rho - &
      screened_rates(k_p_sc37__ti38)*Y(jp)*Y(jsc37)*state % rho - &
      screened_rates(k_p_sc38__he4_ca35)*Y(jp)*Y(jsc38)*state % rho - &
      screened_rates(k_p_sc38__ti39)*Y(jp)*Y(jsc38)*state % rho - &
      screened_rates(k_p_sc39__he4_ca36)*Y(jp)*Y(jsc39)*state % rho - &
      screened_rates(k_p_sc39__ti40)*Y(jp)*Y(jsc39)*state % rho - &
      screened_rates(k_p_sc40__he4_ca37)*Y(jp)*Y(jsc40)*state % rho - &
      screened_rates(k_p_sc40__ti41)*Y(jp)*Y(jsc40)*state % rho - &
      screened_rates(k_p_sc41__he4_ca38)*Y(jp)*Y(jsc41)*state % rho - &
      screened_rates(k_p_sc41__ti42)*Y(jp)*Y(jsc41)*state % rho - &
      screened_rates(k_p_sc42__he4_ca39)*Y(jp)*Y(jsc42)*state % rho - &
      screened_rates(k_p_sc42__ti43)*Y(jp)*Y(jsc42)*state % rho - &
      screened_rates(k_p_sc43__he4_ca40)*Y(jp)*Y(jsc43)*state % rho - &
      screened_rates(k_p_sc43__ti44)*Y(jp)*Y(jsc43)*state % rho - &
      screened_rates(k_p_sc44__he4_ca41)*Y(jp)*Y(jsc44)*state % rho - &
      screened_rates(k_p_sc44__ti45)*Y(jp)*Y(jsc44)*state % rho - &
      screened_rates(k_p_sc45__he4_ca42)*Y(jp)*Y(jsc45)*state % rho - &
      screened_rates(k_p_sc45__ti46)*Y(jp)*Y(jsc45)*state % rho - &
      screened_rates(k_p_sc46__he4_ca43)*Y(jp)*Y(jsc46)*state % rho - &
      screened_rates(k_p_sc46__ti47)*Y(jp)*Y(jsc46)*state % rho - &
      screened_rates(k_p_si25__he4_al22)*Y(jp)*Y(jsi25)*state % rho - &
      screened_rates(k_p_si25__p26)*Y(jp)*Y(jsi25)*state % rho - &
      screened_rates(k_p_si26__he4_al23)*Y(jp)*Y(jsi26)*state % rho - &
      screened_rates(k_p_si26__p27)*Y(jp)*Y(jsi26)*state % rho - &
      screened_rates(k_p_si27__he4_al24)*Y(jp)*Y(jsi27)*state % rho - &
      screened_rates(k_p_si27__p28)*Y(jp)*Y(jsi27)*state % rho - &
      screened_rates(k_p_si28__he4_al25)*Y(jp)*Y(jsi28)*state % rho - &
      screened_rates(k_p_si28__p29)*Y(jp)*Y(jsi28)*state % rho - &
      screened_rates(k_p_si29__he4_al26)*Y(jp)*Y(jsi29)*state % rho - &
      screened_rates(k_p_si29__p30)*Y(jp)*Y(jsi29)*state % rho - &
      screened_rates(k_p_si30__he4_al27)*Y(jp)*Y(jsi30)*state % rho - &
      screened_rates(k_p_si30__p31)*Y(jp)*Y(jsi30)*state % rho - &
      screened_rates(k_p_ti39__he4_sc36)*Y(jp)*Y(jti39)*state % rho - &
      screened_rates(k_p_ti39__v40)*Y(jp)*Y(jti39)*state % rho - &
      screened_rates(k_p_ti40__he4_sc37)*Y(jp)*Y(jti40)*state % rho - &
      screened_rates(k_p_ti40__v41)*Y(jp)*Y(jti40)*state % rho - &
      screened_rates(k_p_ti41__he4_sc38)*Y(jp)*Y(jti41)*state % rho - &
      screened_rates(k_p_ti41__v42)*Y(jp)*Y(jti41)*state % rho - &
      screened_rates(k_p_ti42__he4_sc39)*Y(jp)*Y(jti42)*state % rho - &
      screened_rates(k_p_ti42__v43)*Y(jp)*Y(jti42)*state % rho - &
      screened_rates(k_p_ti43__he4_sc40)*Y(jp)*Y(jti43)*state % rho - &
      screened_rates(k_p_ti43__v44)*Y(jp)*Y(jti43)*state % rho - &
      screened_rates(k_p_ti44__he4_sc41)*Y(jp)*Y(jti44)*state % rho - &
      screened_rates(k_p_ti44__v45)*Y(jp)*Y(jti44)*state % rho - &
      screened_rates(k_p_ti45__he4_sc42)*Y(jp)*Y(jti45)*state % rho - &
      screened_rates(k_p_ti45__v46)*Y(jp)*Y(jti45)*state % rho - &
      screened_rates(k_p_ti46__he4_sc43)*Y(jp)*Y(jti46)*state % rho - &
      screened_rates(k_p_ti46__v47)*Y(jp)*Y(jti46)*state % rho - &
      screened_rates(k_p_ti47__he4_sc44)*Y(jp)*Y(jti47)*state % rho - &
      screened_rates(k_p_ti47__v48)*Y(jp)*Y(jti47)*state % rho - &
      screened_rates(k_p_ti48__he4_sc45)*Y(jp)*Y(jti48)*state % rho - &
      screened_rates(k_p_ti48__v49)*Y(jp)*Y(jti48)*state % rho - &
      screened_rates(k_p_v41__cr42)*Y(jp)*Y(jv41)*state % rho - &
      screened_rates(k_p_v41__he4_ti38)*Y(jp)*Y(jv41)*state % rho - &
      screened_rates(k_p_v42__cr43)*Y(jp)*Y(jv42)*state % rho - &
      screened_rates(k_p_v42__he4_ti39)*Y(jp)*Y(jv42)*state % rho - &
      screened_rates(k_p_v43__cr44)*Y(jp)*Y(jv43)*state % rho - &
      screened_rates(k_p_v43__he4_ti40)*Y(jp)*Y(jv43)*state % rho - &
      screened_rates(k_p_v44__cr45)*Y(jp)*Y(jv44)*state % rho - &
      screened_rates(k_p_v44__he4_ti41)*Y(jp)*Y(jv44)*state % rho - &
      screened_rates(k_p_v45__cr46)*Y(jp)*Y(jv45)*state % rho - &
      screened_rates(k_p_v45__he4_ti42)*Y(jp)*Y(jv45)*state % rho - &
      screened_rates(k_p_v46__cr47)*Y(jp)*Y(jv46)*state % rho - &
      screened_rates(k_p_v46__he4_ti43)*Y(jp)*Y(jv46)*state % rho - &
      screened_rates(k_p_v47__cr48)*Y(jp)*Y(jv47)*state % rho - &
      screened_rates(k_p_v47__he4_ti44)*Y(jp)*Y(jv47)*state % rho - &
      screened_rates(k_p_v48__cr49)*Y(jp)*Y(jv48)*state % rho - &
      screened_rates(k_p_v48__he4_ti45)*Y(jp)*Y(jv48)*state % rho - &
      screened_rates(k_p_v49__cr50)*Y(jp)*Y(jv49)*state % rho - &
      screened_rates(k_p_v49__he4_ti46)*Y(jp)*Y(jv49)*state % rho - &
      screened_rates(k_p_v50__cr51)*Y(jp)*Y(jv50)*state % rho - &
      screened_rates(k_p_v50__he4_ti47)*Y(jp)*Y(jv50)*state % rho + &
      screened_rates(k_s28__p_p27)*Y(js28) + screened_rates(k_s28__p_si27__weak__wc12)* &
      Y(js28) + screened_rates(k_s29__p_p28)*Y(js29) + &
      screened_rates(k_s29__p_si28__weak__wc12)*Y(js29) + screened_rates(k_s30__p_p29)* &
      Y(js30) + screened_rates(k_s31__p_p30)*Y(js31) + screened_rates(k_s32__p_p31)* &
      Y(js32) + screened_rates(k_s33__p_p32)*Y(js33) + screened_rates(k_s34__p_p33)* &
      Y(js34) + screened_rates(k_sc36__p_ca35)*Y(jsc36) + &
      screened_rates(k_sc37__p_ca36)*Y(jsc37) + screened_rates(k_sc38__p_ca37)* &
      Y(jsc38) + screened_rates(k_sc39__p_ca38)*Y(jsc39) + &
      screened_rates(k_sc40__p_ca39)*Y(jsc40) + &
      screened_rates(k_sc40__p_k39__weak__wc12)*Y(jsc40) + &
      screened_rates(k_sc41__p_ca40)*Y(jsc41) + screened_rates(k_sc42__p_ca41)* &
      Y(jsc42) + screened_rates(k_sc43__p_ca42)*Y(jsc43) + &
      screened_rates(k_sc44__p_ca43)*Y(jsc44) + screened_rates(k_sc45__p_ca44)* &
      Y(jsc45) + screened_rates(k_si24__p_al23)*Y(jsi24) + &
      screened_rates(k_si24__p_mg23__weak__wc12)*Y(jsi24) + &
      screened_rates(k_si25__p_al24)*Y(jsi25) + &
      screened_rates(k_si25__p_mg24__weak__wc12)*Y(jsi25) + &
      screened_rates(k_si26__p_al25)*Y(jsi26) + screened_rates(k_si27__p_al26)* &
      Y(jsi27) + screened_rates(k_si28__p_al27)*Y(jsi28) + &
      screened_rates(k_si29__p_al28)*Y(jsi29) + screened_rates(k_ti38__p_sc37)* &
      Y(jti38) + screened_rates(k_ti39__p_ca38__weak__wc12)*Y(jti39) + &
      screened_rates(k_ti39__p_sc38)*Y(jti39) + &
      screened_rates(k_ti40__p_ca39__weak__wc12)*Y(jti40) + &
      screened_rates(k_ti40__p_sc39)*Y(jti40) + &
      screened_rates(k_ti41__p_ca40__weak__wc12)*Y(jti41) + &
      screened_rates(k_ti41__p_sc40)*Y(jti41) + screened_rates(k_ti42__p_sc41)* &
      Y(jti42) + screened_rates(k_ti43__p_sc42)*Y(jti43) + &
      screened_rates(k_ti44__p_sc43)*Y(jti44) + screened_rates(k_ti45__p_sc44)* &
      Y(jti45) + screened_rates(k_ti46__p_sc45)*Y(jti46) + &
      screened_rates(k_ti47__p_sc46)*Y(jti47) + screened_rates(k_v40__p_ti39)*Y(jv40) &
      + screened_rates(k_v41__p_ti40)*Y(jv41) + screened_rates(k_v42__p_ti41)*Y(jv42) &
      + screened_rates(k_v43__p_ti42)*Y(jv43) + screened_rates(k_v44__p_ti43)*Y(jv44) &
      + screened_rates(k_v45__p_ti44)*Y(jv45) + screened_rates(k_v46__p_ti45)*Y(jv46) &
      + screened_rates(k_v47__p_ti46)*Y(jv47) + screened_rates(k_v48__p_ti47)*Y(jv48) &
      + screened_rates(k_v49__p_ti48)*Y(jv49) &
       )

    ydot_nuc(jd) = ( &
      -screened_rates(k_d_be7__p_he4_he4)*Y(jbe7)*Y(jd)*state % rho - &
      screened_rates(k_d_d__he4)*Y(jd)**2*state % rho - screened_rates(k_d_he3__p_he4)* &
      Y(jd)*Y(jhe3)*state % rho + screened_rates(k_he3__p_d)*Y(jhe3) + 2.0d0* &
      screened_rates(k_he4__d_d)*Y(jhe4) - screened_rates(k_p_d__he3)*Y(jd)*Y(jp)* &
      state % rho + screened_rates(k_p_he4__d_he3)*Y(jhe4)*Y(jp)*state % rho + &
      0.5d0*screened_rates(k_p_he4_he4__d_be7)*Y(jhe4)**2*Y(jp)*state % rho**2 &
      + 0.5d0*screened_rates(k_p_p__d__weak__bet_pos_)*Y(jp)**2*state % rho + &
      0.5d0*screened_rates(k_p_p__d__weak__electron_capture)*Y(jp)**2*state % rho &
      **2*state % y_e &
       )

    ydot_nuc(jhe3) = ( &
      screened_rates(k_be7__he4_he3)*Y(jbe7) - screened_rates(k_d_he3__p_he4)*Y(jd)* &
      Y(jhe3)*state % rho - screened_rates(k_he3__p_d)*Y(jhe3) - &
      screened_rates(k_he3_be7__p_p_he4_he4)*Y(jbe7)*Y(jhe3)*state % rho - &
      screened_rates(k_he3_he3__p_p_he4)*Y(jhe3)**2*state % rho - &
      screened_rates(k_he4_he3__be7)*Y(jhe3)*Y(jhe4)*state % rho + &
      screened_rates(k_p_d__he3)*Y(jd)*Y(jp)*state % rho - &
      screened_rates(k_p_he3__he4__weak__bet_pos_)*Y(jhe3)*Y(jp)*state % rho + &
      screened_rates(k_p_he4__d_he3)*Y(jhe4)*Y(jp)*state % rho + &
      screened_rates(k_p_p_he4__he3_he3)*Y(jhe4)*Y(jp)**2*state % rho**2 + &
      0.25d0*screened_rates(k_p_p_he4_he4__he3_be7)*Y(jhe4)**2*Y(jp)**2* &
      state % rho**3 &
       )

    ydot_nuc(jhe4) = ( &
      screened_rates(k_al22__he4_ne18__weak__wc12)*Y(jal22) + screened_rates(k_al23__he4_na19) &
      *Y(jal23) + screened_rates(k_al24__he4_na20)*Y(jal24) + &
      screened_rates(k_al24__he4_ne20__weak__wc12)*Y(jal24) + &
      screened_rates(k_al25__he4_na21)*Y(jal25) + screened_rates(k_al26__he4_na22)* &
      Y(jal26) + screened_rates(k_al27__he4_na23)*Y(jal27) + &
      screened_rates(k_al28__he4_na24)*Y(jal28) + screened_rates(k_ar32__he4_s28)* &
      Y(jar32) + screened_rates(k_ar33__he4_s29)*Y(jar33) + &
      screened_rates(k_ar34__he4_s30)*Y(jar34) + screened_rates(k_ar35__he4_s31)* &
      Y(jar35) + screened_rates(k_ar36__he4_s32)*Y(jar36) + &
      screened_rates(k_ar37__he4_s33)*Y(jar37) + screened_rates(k_ar38__he4_s34)* &
      Y(jar38) + screened_rates(k_ar39__he4_s35)*Y(jar39) + 2.0d0* &
      screened_rates(k_b8__he4_he4__weak__wc12)*Y(jb8) + screened_rates(k_be7__he4_he3) &
      *Y(jbe7) + 3.0d0*screened_rates(k_c12__he4_he4_he4)*Y(jc12) + 0.5d0* &
      screened_rates(k_c12_c12__he4_ne20)*Y(jc12)**2*state % rho + &
      screened_rates(k_ca35__he4_ar31)*Y(jca35) + screened_rates(k_ca36__he4_ar32)* &
      Y(jca36) + screened_rates(k_ca37__he4_ar33)*Y(jca37) + &
      screened_rates(k_ca38__he4_ar34)*Y(jca38) + screened_rates(k_ca39__he4_ar35)* &
      Y(jca39) + screened_rates(k_ca40__he4_ar36)*Y(jca40) + &
      screened_rates(k_ca41__he4_ar37)*Y(jca41) + screened_rates(k_ca42__he4_ar38)* &
      Y(jca42) + screened_rates(k_ca43__he4_ar39)*Y(jca43) + &
      screened_rates(k_cl30__he4_p26)*Y(jcl30) + screened_rates(k_cl31__he4_p27)* &
      Y(jcl31) + screened_rates(k_cl32__he4_p28)*Y(jcl32) + &
      screened_rates(k_cl32__he4_si28__weak__wc12)*Y(jcl32) + &
      screened_rates(k_cl33__he4_p29)*Y(jcl33) + screened_rates(k_cl34__he4_p30)* &
      Y(jcl34) + screened_rates(k_cl35__he4_p31)*Y(jcl35) + &
      screened_rates(k_cl36__he4_p32)*Y(jcl36) + screened_rates(k_cl37__he4_p33)* &
      Y(jcl37) + screened_rates(k_co48__he4_mn44)*Y(jco48) + &
      screened_rates(k_co49__he4_mn45)*Y(jco49) + screened_rates(k_co50__he4_mn46)* &
      Y(jco50) + screened_rates(k_co51__he4_mn47)*Y(jco51) + &
      screened_rates(k_co52__he4_mn48)*Y(jco52) + screened_rates(k_co53__he4_mn49)* &
      Y(jco53) + screened_rates(k_co54__he4_mn50)*Y(jco54) + &
      screened_rates(k_co55__he4_mn51)*Y(jco55) + screened_rates(k_co56__he4_mn52)* &
      Y(jco56) + screened_rates(k_cr42__he4_ti38)*Y(jcr42) + &
      screened_rates(k_cr43__he4_ti39)*Y(jcr43) + screened_rates(k_cr44__he4_ti40)* &
      Y(jcr44) + screened_rates(k_cr45__he4_ti41)*Y(jcr45) + &
      screened_rates(k_cr46__he4_ti42)*Y(jcr46) + screened_rates(k_cr47__he4_ti43)* &
      Y(jcr47) + screened_rates(k_cr48__he4_ti44)*Y(jcr48) + &
      screened_rates(k_cr49__he4_ti45)*Y(jcr49) + screened_rates(k_cr50__he4_ti46)* &
      Y(jcr50) + screened_rates(k_cr51__he4_ti47)*Y(jcr51) + &
      screened_rates(k_cr52__he4_ti48)*Y(jcr52) + 2.0d0* &
      screened_rates(k_d_be7__p_he4_he4)*Y(jbe7)*Y(jd)*state % rho + 0.5d0* &
      screened_rates(k_d_d__he4)*Y(jd)**2*state % rho + screened_rates(k_d_he3__p_he4)* &
      Y(jd)*Y(jhe3)*state % rho + screened_rates(k_f18__he4_n14)*Y(jf18) + &
      screened_rates(k_f19__he4_n15)*Y(jf19) + screened_rates(k_fe46__he4_cr42)* &
      Y(jfe46) + screened_rates(k_fe47__he4_cr43)*Y(jfe47) + &
      screened_rates(k_fe48__he4_cr44)*Y(jfe48) + screened_rates(k_fe49__he4_cr45)* &
      Y(jfe49) + screened_rates(k_fe50__he4_cr46)*Y(jfe50) + &
      screened_rates(k_fe51__he4_cr47)*Y(jfe51) + screened_rates(k_fe52__he4_cr48)* &
      Y(jfe52) + screened_rates(k_fe53__he4_cr49)*Y(jfe53) + &
      screened_rates(k_fe54__he4_cr50)*Y(jfe54) + screened_rates(k_fe55__he4_cr51)* &
      Y(jfe55) + screened_rates(k_fe56__he4_cr52)*Y(jfe56) + 2.0d0* &
      screened_rates(k_he3_be7__p_p_he4_he4)*Y(jbe7)*Y(jhe3)*state % rho + &
      0.5d0*screened_rates(k_he3_he3__p_p_he4)*Y(jhe3)**2*state % rho - &
      screened_rates(k_he4__d_d)*Y(jhe4) - screened_rates(k_he4_al22__p26)*Y(jal22)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_al22__p_si25)*Y(jal22)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_al23__p27)*Y(jal23)*Y(jhe4)* &
      state % rho - screened_rates(k_he4_al23__p_si26)*Y(jal23)*Y(jhe4)* &
      state % rho - screened_rates(k_he4_al24__p28)*Y(jal24)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al24__p_si27)*Y(jal24)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al25__p29)*Y(jal25)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al25__p_si28)*Y(jal25)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al26__p30)*Y(jal26)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al26__p_si29)*Y(jal26)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al27__p31)*Y(jal27)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al27__p_si30)*Y(jal27)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al28__p32)*Y(jal28)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar31__ca35)*Y(jar31)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar31__p_k34)*Y(jar31)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar32__ca36)*Y(jar32)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar32__p_k35)*Y(jar32)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar33__ca37)*Y(jar33)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar33__p_k36)*Y(jar33)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar34__ca38)*Y(jar34)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar34__p_k37)*Y(jar34)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar35__ca39)*Y(jar35)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar35__p_k38)*Y(jar35)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar36__ca40)*Y(jar36)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar36__p_k39)*Y(jar36)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar37__ca41)*Y(jar37)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar37__p_k40)*Y(jar37)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar38__ca42)*Y(jar38)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar38__p_k41)*Y(jar38)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar39__ca43)*Y(jar39)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_c12__p_n15)*Y(jc12)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca34__p_sc37)*Y(jca34)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca34__ti38)*Y(jca34)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca35__p_sc38)*Y(jca35)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca35__ti39)*Y(jca35)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca36__p_p_ca38)*Y(jca36)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca36__p_sc39)*Y(jca36)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca36__ti40)*Y(jca36)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca37__p_sc40)*Y(jca37)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca37__ti41)*Y(jca37)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca38__p_sc41)*Y(jca38)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca38__ti42)*Y(jca38)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca39__p_sc42)*Y(jca39)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca39__ti43)*Y(jca39)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca40__p_sc43)*Y(jca40)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca40__ti44)*Y(jca40)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca41__p_sc44)*Y(jca41)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca41__ti45)*Y(jca41)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca42__p_sc45)*Y(jca42)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca42__ti46)*Y(jca42)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca43__p_sc46)*Y(jca43)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca43__ti47)*Y(jca43)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca44__ti48)*Y(jca44)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl29__k33)*Y(jcl29)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl29__p_ar32)*Y(jcl29)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl30__k34)*Y(jcl30)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl30__p_ar33)*Y(jcl30)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl31__k35)*Y(jcl31)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl31__p_ar34)*Y(jcl31)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl32__k36)*Y(jcl32)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl32__p_ar35)*Y(jcl32)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl33__k37)*Y(jcl33)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl33__p_ar36)*Y(jcl33)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl34__k38)*Y(jcl34)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl34__p_ar37)*Y(jcl34)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl35__k39)*Y(jcl35)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl35__p_ar38)*Y(jcl35)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl36__k40)*Y(jcl36)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl36__p_ar39)*Y(jcl36)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl37__k41)*Y(jcl37)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_co47__p_ni50)*Y(jco47)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_co48__p_ni51)*Y(jco48)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_co49__p_ni52)*Y(jco49)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_co50__p_ni53)*Y(jco50)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_co51__p_ni54)*Y(jco51)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_co52__p_ni55)*Y(jco52)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_co53__p_ni56)*Y(jco53)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr42__fe46)*Y(jcr42)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr42__p_mn45)*Y(jcr42)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr43__fe47)*Y(jcr43)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr43__p_mn46)*Y(jcr43)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr44__fe48)*Y(jcr44)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr44__p_mn47)*Y(jcr44)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr45__fe49)*Y(jcr45)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr45__p_mn48)*Y(jcr45)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr46__fe50)*Y(jcr46)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr46__p_mn49)*Y(jcr46)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr47__fe51)*Y(jcr47)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr47__p_mn50)*Y(jcr47)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr48__fe52)*Y(jcr48)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr48__p_mn51)*Y(jcr48)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr49__fe53)*Y(jcr49)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr49__p_mn52)*Y(jcr49)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr50__fe54)*Y(jcr50)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr50__p_mn53)*Y(jcr50)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr51__fe55)*Y(jcr51)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr51__p_mn54)*Y(jcr51)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr52__fe56)*Y(jcr52)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr52__p_mn55)*Y(jcr52)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_f16__na20)*Y(jf16)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_f16__p_ne19)*Y(jf16)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_f17__na21)*Y(jf17)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_f17__p_ne20)*Y(jf17)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_f18__na22)*Y(jf18)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_f18__p_ne21)*Y(jf18)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_f19__na23)*Y(jf19)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_f19__p_ne22)*Y(jf19)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_f20__na24)*Y(jf20)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe45__ni49)*Y(jfe45)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe45__p_co48)*Y(jfe45)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe46__ni50)*Y(jfe46)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe46__p_co49)*Y(jfe46)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe47__ni51)*Y(jfe47)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe47__p_co50)*Y(jfe47)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe48__ni52)*Y(jfe48)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe48__p_co51)*Y(jfe48)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe49__ni53)*Y(jfe49)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe49__p_co52)*Y(jfe49)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe50__ni54)*Y(jfe50)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe50__p_co53)*Y(jfe50)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe51__ni55)*Y(jfe51)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe51__p_co54)*Y(jfe51)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe52__ni56)*Y(jfe52)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe52__p_co55)*Y(jfe52)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe53__p_co56)*Y(jfe53)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_he3__be7)*Y(jhe3)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_he4__p_li7)*Y(jhe4)**2*state % rho - 0.5d0* &
      screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**3*state % rho**2 - &
      screened_rates(k_he4_k33__p_ca36)*Y(jhe4)*Y(jk33)*state % rho - &
      screened_rates(k_he4_k33__sc37)*Y(jhe4)*Y(jk33)*state % rho - &
      screened_rates(k_he4_k34__p_ca37)*Y(jhe4)*Y(jk34)*state % rho - &
      screened_rates(k_he4_k34__sc38)*Y(jhe4)*Y(jk34)*state % rho - &
      screened_rates(k_he4_k35__p_ca38)*Y(jhe4)*Y(jk35)*state % rho - &
      screened_rates(k_he4_k35__sc39)*Y(jhe4)*Y(jk35)*state % rho - &
      screened_rates(k_he4_k36__p_ca39)*Y(jhe4)*Y(jk36)*state % rho - &
      screened_rates(k_he4_k36__sc40)*Y(jhe4)*Y(jk36)*state % rho - &
      screened_rates(k_he4_k37__p_ca40)*Y(jhe4)*Y(jk37)*state % rho - &
      screened_rates(k_he4_k37__sc41)*Y(jhe4)*Y(jk37)*state % rho - &
      screened_rates(k_he4_k38__p_ca41)*Y(jhe4)*Y(jk38)*state % rho - &
      screened_rates(k_he4_k38__sc42)*Y(jhe4)*Y(jk38)*state % rho - &
      screened_rates(k_he4_k39__p_ca42)*Y(jhe4)*Y(jk39)*state % rho - &
      screened_rates(k_he4_k39__sc43)*Y(jhe4)*Y(jk39)*state % rho - &
      screened_rates(k_he4_k40__p_ca43)*Y(jhe4)*Y(jk40)*state % rho - &
      screened_rates(k_he4_k40__sc44)*Y(jhe4)*Y(jk40)*state % rho - &
      screened_rates(k_he4_k41__p_ca44)*Y(jhe4)*Y(jk41)*state % rho - &
      screened_rates(k_he4_k41__sc45)*Y(jhe4)*Y(jk41)*state % rho - &
      screened_rates(k_he4_mg21__p_al24)*Y(jhe4)*Y(jmg21)*state % rho - &
      screened_rates(k_he4_mg21__si25)*Y(jhe4)*Y(jmg21)*state % rho - &
      screened_rates(k_he4_mg22__p_al25)*Y(jhe4)*Y(jmg22)*state % rho - &
      screened_rates(k_he4_mg22__si26)*Y(jhe4)*Y(jmg22)*state % rho - &
      screened_rates(k_he4_mg23__p_al26)*Y(jhe4)*Y(jmg23)*state % rho - &
      screened_rates(k_he4_mg23__si27)*Y(jhe4)*Y(jmg23)*state % rho - &
      screened_rates(k_he4_mg24__p_al27)*Y(jhe4)*Y(jmg24)*state % rho - &
      screened_rates(k_he4_mg24__si28)*Y(jhe4)*Y(jmg24)*state % rho - &
      screened_rates(k_he4_mg25__p_al28)*Y(jhe4)*Y(jmg25)*state % rho - &
      screened_rates(k_he4_mg25__si29)*Y(jhe4)*Y(jmg25)*state % rho - &
      screened_rates(k_he4_mg26__si30)*Y(jhe4)*Y(jmg26)*state % rho - &
      screened_rates(k_he4_mn44__co48)*Y(jhe4)*Y(jmn44)*state % rho - &
      screened_rates(k_he4_mn44__p_fe47)*Y(jhe4)*Y(jmn44)*state % rho - &
      screened_rates(k_he4_mn45__co49)*Y(jhe4)*Y(jmn45)*state % rho - &
      screened_rates(k_he4_mn45__p_fe48)*Y(jhe4)*Y(jmn45)*state % rho - &
      screened_rates(k_he4_mn46__co50)*Y(jhe4)*Y(jmn46)*state % rho - &
      screened_rates(k_he4_mn46__p_fe49)*Y(jhe4)*Y(jmn46)*state % rho - &
      screened_rates(k_he4_mn47__co51)*Y(jhe4)*Y(jmn47)*state % rho - &
      screened_rates(k_he4_mn47__p_fe50)*Y(jhe4)*Y(jmn47)*state % rho - &
      screened_rates(k_he4_mn48__co52)*Y(jhe4)*Y(jmn48)*state % rho - &
      screened_rates(k_he4_mn48__p_fe51)*Y(jhe4)*Y(jmn48)*state % rho - &
      screened_rates(k_he4_mn49__co53)*Y(jhe4)*Y(jmn49)*state % rho - &
      screened_rates(k_he4_mn49__p_fe52)*Y(jhe4)*Y(jmn49)*state % rho - &
      screened_rates(k_he4_mn50__co54)*Y(jhe4)*Y(jmn50)*state % rho - &
      screened_rates(k_he4_mn50__p_fe53)*Y(jhe4)*Y(jmn50)*state % rho - &
      screened_rates(k_he4_mn51__co55)*Y(jhe4)*Y(jmn51)*state % rho - &
      screened_rates(k_he4_mn51__p_fe54)*Y(jhe4)*Y(jmn51)*state % rho - &
      screened_rates(k_he4_mn52__co56)*Y(jhe4)*Y(jmn52)*state % rho - &
      screened_rates(k_he4_mn52__p_fe55)*Y(jhe4)*Y(jmn52)*state % rho - &
      screened_rates(k_he4_mn53__p_fe56)*Y(jhe4)*Y(jmn53)*state % rho - &
      screened_rates(k_he4_n12__p_o15)*Y(jhe4)*Y(jn12)*state % rho - &
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*Y(jn13)*state % rho - &
      screened_rates(k_he4_n14__f18)*Y(jhe4)*Y(jn14)*state % rho - &
      screened_rates(k_he4_n14__p_o17)*Y(jhe4)*Y(jn14)*state % rho - &
      screened_rates(k_he4_n15__f19)*Y(jhe4)*Y(jn15)*state % rho - &
      screened_rates(k_he4_n15__p_o18)*Y(jhe4)*Y(jn15)*state % rho - &
      screened_rates(k_he4_na19__al23)*Y(jhe4)*Y(jna19)*state % rho - &
      screened_rates(k_he4_na19__p_mg22)*Y(jhe4)*Y(jna19)*state % rho - &
      screened_rates(k_he4_na20__al24)*Y(jhe4)*Y(jna20)*state % rho - &
      screened_rates(k_he4_na20__p_mg23)*Y(jhe4)*Y(jna20)*state % rho - &
      screened_rates(k_he4_na21__al25)*Y(jhe4)*Y(jna21)*state % rho - &
      screened_rates(k_he4_na21__p_mg24)*Y(jhe4)*Y(jna21)*state % rho - &
      screened_rates(k_he4_na22__al26)*Y(jhe4)*Y(jna22)*state % rho - &
      screened_rates(k_he4_na22__p_mg25)*Y(jhe4)*Y(jna22)*state % rho - &
      screened_rates(k_he4_na23__al27)*Y(jhe4)*Y(jna23)*state % rho - &
      screened_rates(k_he4_na23__p_mg26)*Y(jhe4)*Y(jna23)*state % rho - &
      screened_rates(k_he4_na24__al28)*Y(jhe4)*Y(jna24)*state % rho - &
      screened_rates(k_he4_ne17__mg21)*Y(jhe4)*Y(jne17)*state % rho - &
      screened_rates(k_he4_ne17__p_na20)*Y(jhe4)*Y(jne17)*state % rho - &
      screened_rates(k_he4_ne18__mg22)*Y(jhe4)*Y(jne18)*state % rho - &
      screened_rates(k_he4_ne18__p_na21)*Y(jhe4)*Y(jne18)*state % rho - &
      screened_rates(k_he4_ne19__mg23)*Y(jhe4)*Y(jne19)*state % rho - &
      screened_rates(k_he4_ne19__p_na22)*Y(jhe4)*Y(jne19)*state % rho - &
      screened_rates(k_he4_ne20__c12_c12)*Y(jhe4)*Y(jne20)*state % rho - &
      screened_rates(k_he4_ne20__mg24)*Y(jhe4)*Y(jne20)*state % rho - &
      screened_rates(k_he4_ne20__p_na23)*Y(jhe4)*Y(jne20)*state % rho - &
      screened_rates(k_he4_ne21__mg25)*Y(jhe4)*Y(jne21)*state % rho - &
      screened_rates(k_he4_ne21__p_na24)*Y(jhe4)*Y(jne21)*state % rho - &
      screened_rates(k_he4_ne22__mg26)*Y(jhe4)*Y(jne22)*state % rho - &
      screened_rates(k_he4_o14__ne18)*Y(jhe4)*Y(jo14)*state % rho - &
      screened_rates(k_he4_o14__p_f17)*Y(jhe4)*Y(jo14)*state % rho - &
      screened_rates(k_he4_o15__ne19)*Y(jhe4)*Y(jo15)*state % rho - &
      screened_rates(k_he4_o15__p_f18)*Y(jhe4)*Y(jo15)*state % rho - &
      screened_rates(k_he4_o16__ne20)*Y(jhe4)*Y(jo16)*state % rho - &
      screened_rates(k_he4_o16__p_f19)*Y(jhe4)*Y(jo16)*state % rho - &
      screened_rates(k_he4_o17__ne21)*Y(jhe4)*Y(jo17)*state % rho - &
      screened_rates(k_he4_o17__p_f20)*Y(jhe4)*Y(jo17)*state % rho - &
      screened_rates(k_he4_o18__ne22)*Y(jhe4)*Y(jo18)*state % rho - &
      screened_rates(k_he4_p26__cl30)*Y(jhe4)*Y(jp26)*state % rho - &
      screened_rates(k_he4_p26__p_s29)*Y(jhe4)*Y(jp26)*state % rho - &
      screened_rates(k_he4_p27__cl31)*Y(jhe4)*Y(jp27)*state % rho - &
      screened_rates(k_he4_p27__p_s30)*Y(jhe4)*Y(jp27)*state % rho - &
      screened_rates(k_he4_p28__cl32)*Y(jhe4)*Y(jp28)*state % rho - &
      screened_rates(k_he4_p28__p_s31)*Y(jhe4)*Y(jp28)*state % rho - &
      screened_rates(k_he4_p29__cl33)*Y(jhe4)*Y(jp29)*state % rho - &
      screened_rates(k_he4_p29__p_s32)*Y(jhe4)*Y(jp29)*state % rho - &
      screened_rates(k_he4_p30__cl34)*Y(jhe4)*Y(jp30)*state % rho - &
      screened_rates(k_he4_p30__p_s33)*Y(jhe4)*Y(jp30)*state % rho - &
      screened_rates(k_he4_p31__cl35)*Y(jhe4)*Y(jp31)*state % rho - &
      screened_rates(k_he4_p31__p_s34)*Y(jhe4)*Y(jp31)*state % rho - &
      screened_rates(k_he4_p32__cl36)*Y(jhe4)*Y(jp32)*state % rho - &
      screened_rates(k_he4_p32__p_s35)*Y(jhe4)*Y(jp32)*state % rho - &
      screened_rates(k_he4_p33__cl37)*Y(jhe4)*Y(jp33)*state % rho - &
      screened_rates(k_he4_s28__ar32)*Y(jhe4)*Y(js28)*state % rho - &
      screened_rates(k_he4_s28__p_cl31)*Y(jhe4)*Y(js28)*state % rho - &
      screened_rates(k_he4_s29__ar33)*Y(jhe4)*Y(js29)*state % rho - &
      screened_rates(k_he4_s29__p_cl32)*Y(jhe4)*Y(js29)*state % rho - &
      screened_rates(k_he4_s30__ar34)*Y(jhe4)*Y(js30)*state % rho - &
      screened_rates(k_he4_s30__p_cl33)*Y(jhe4)*Y(js30)*state % rho - &
      screened_rates(k_he4_s31__ar35)*Y(jhe4)*Y(js31)*state % rho - &
      screened_rates(k_he4_s31__p_cl34)*Y(jhe4)*Y(js31)*state % rho - &
      screened_rates(k_he4_s32__ar36)*Y(jhe4)*Y(js32)*state % rho - &
      screened_rates(k_he4_s32__p_cl35)*Y(jhe4)*Y(js32)*state % rho - &
      screened_rates(k_he4_s33__ar37)*Y(jhe4)*Y(js33)*state % rho - &
      screened_rates(k_he4_s33__p_cl36)*Y(jhe4)*Y(js33)*state % rho - &
      screened_rates(k_he4_s34__ar38)*Y(jhe4)*Y(js34)*state % rho - &
      screened_rates(k_he4_s34__p_cl37)*Y(jhe4)*Y(js34)*state % rho - &
      screened_rates(k_he4_s35__ar39)*Y(jhe4)*Y(js35)*state % rho - &
      screened_rates(k_he4_sc36__p_ti39)*Y(jhe4)*Y(jsc36)*state % rho - &
      screened_rates(k_he4_sc36__v40)*Y(jhe4)*Y(jsc36)*state % rho - &
      screened_rates(k_he4_sc37__p_ti40)*Y(jhe4)*Y(jsc37)*state % rho - &
      screened_rates(k_he4_sc37__v41)*Y(jhe4)*Y(jsc37)*state % rho - &
      screened_rates(k_he4_sc38__p_ti41)*Y(jhe4)*Y(jsc38)*state % rho - &
      screened_rates(k_he4_sc38__v42)*Y(jhe4)*Y(jsc38)*state % rho - &
      screened_rates(k_he4_sc39__p_ti42)*Y(jhe4)*Y(jsc39)*state % rho - &
      screened_rates(k_he4_sc39__v43)*Y(jhe4)*Y(jsc39)*state % rho - &
      screened_rates(k_he4_sc40__p_ti43)*Y(jhe4)*Y(jsc40)*state % rho - &
      screened_rates(k_he4_sc40__v44)*Y(jhe4)*Y(jsc40)*state % rho - &
      screened_rates(k_he4_sc41__p_ti44)*Y(jhe4)*Y(jsc41)*state % rho - &
      screened_rates(k_he4_sc41__v45)*Y(jhe4)*Y(jsc41)*state % rho - &
      screened_rates(k_he4_sc42__p_ti45)*Y(jhe4)*Y(jsc42)*state % rho - &
      screened_rates(k_he4_sc42__v46)*Y(jhe4)*Y(jsc42)*state % rho - &
      screened_rates(k_he4_sc43__p_ti46)*Y(jhe4)*Y(jsc43)*state % rho - &
      screened_rates(k_he4_sc43__v47)*Y(jhe4)*Y(jsc43)*state % rho - &
      screened_rates(k_he4_sc44__p_ti47)*Y(jhe4)*Y(jsc44)*state % rho - &
      screened_rates(k_he4_sc44__v48)*Y(jhe4)*Y(jsc44)*state % rho - &
      screened_rates(k_he4_sc45__p_ti48)*Y(jhe4)*Y(jsc45)*state % rho - &
      screened_rates(k_he4_sc45__v49)*Y(jhe4)*Y(jsc45)*state % rho - &
      screened_rates(k_he4_sc46__v50)*Y(jhe4)*Y(jsc46)*state % rho - &
      screened_rates(k_he4_si24__p_p27)*Y(jhe4)*Y(jsi24)*state % rho - &
      screened_rates(k_he4_si24__s28)*Y(jhe4)*Y(jsi24)*state % rho - &
      screened_rates(k_he4_si25__p_p28)*Y(jhe4)*Y(jsi25)*state % rho - &
      screened_rates(k_he4_si25__s29)*Y(jhe4)*Y(jsi25)*state % rho - &
      screened_rates(k_he4_si26__p_p29)*Y(jhe4)*Y(jsi26)*state % rho - &
      screened_rates(k_he4_si26__s30)*Y(jhe4)*Y(jsi26)*state % rho - &
      screened_rates(k_he4_si27__p_p30)*Y(jhe4)*Y(jsi27)*state % rho - &
      screened_rates(k_he4_si27__s31)*Y(jhe4)*Y(jsi27)*state % rho - &
      screened_rates(k_he4_si28__p_p31)*Y(jhe4)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_si28__s32)*Y(jhe4)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_si29__p_p32)*Y(jhe4)*Y(jsi29)*state % rho - &
      screened_rates(k_he4_si29__s33)*Y(jhe4)*Y(jsi29)*state % rho - &
      screened_rates(k_he4_si30__p_p33)*Y(jhe4)*Y(jsi30)*state % rho - &
      screened_rates(k_he4_si30__s34)*Y(jhe4)*Y(jsi30)*state % rho - &
      screened_rates(k_he4_ti38__cr42)*Y(jhe4)*Y(jti38)*state % rho - &
      screened_rates(k_he4_ti38__p_v41)*Y(jhe4)*Y(jti38)*state % rho - &
      screened_rates(k_he4_ti39__cr43)*Y(jhe4)*Y(jti39)*state % rho - &
      screened_rates(k_he4_ti39__p_v42)*Y(jhe4)*Y(jti39)*state % rho - &
      screened_rates(k_he4_ti40__cr44)*Y(jhe4)*Y(jti40)*state % rho - &
      screened_rates(k_he4_ti40__p_v43)*Y(jhe4)*Y(jti40)*state % rho - &
      screened_rates(k_he4_ti41__cr45)*Y(jhe4)*Y(jti41)*state % rho - &
      screened_rates(k_he4_ti41__p_v44)*Y(jhe4)*Y(jti41)*state % rho - &
      screened_rates(k_he4_ti42__cr46)*Y(jhe4)*Y(jti42)*state % rho - &
      screened_rates(k_he4_ti42__p_v45)*Y(jhe4)*Y(jti42)*state % rho - &
      screened_rates(k_he4_ti43__cr47)*Y(jhe4)*Y(jti43)*state % rho - &
      screened_rates(k_he4_ti43__p_v46)*Y(jhe4)*Y(jti43)*state % rho - &
      screened_rates(k_he4_ti44__cr48)*Y(jhe4)*Y(jti44)*state % rho - &
      screened_rates(k_he4_ti44__p_v47)*Y(jhe4)*Y(jti44)*state % rho - &
      screened_rates(k_he4_ti45__cr49)*Y(jhe4)*Y(jti45)*state % rho - &
      screened_rates(k_he4_ti45__p_v48)*Y(jhe4)*Y(jti45)*state % rho - &
      screened_rates(k_he4_ti46__cr50)*Y(jhe4)*Y(jti46)*state % rho - &
      screened_rates(k_he4_ti46__p_v49)*Y(jhe4)*Y(jti46)*state % rho - &
      screened_rates(k_he4_ti47__cr51)*Y(jhe4)*Y(jti47)*state % rho - &
      screened_rates(k_he4_ti47__p_v50)*Y(jhe4)*Y(jti47)*state % rho - &
      screened_rates(k_he4_ti48__cr52)*Y(jhe4)*Y(jti48)*state % rho - &
      screened_rates(k_he4_v40__mn44)*Y(jhe4)*Y(jv40)*state % rho - &
      screened_rates(k_he4_v40__p_cr43)*Y(jhe4)*Y(jv40)*state % rho - &
      screened_rates(k_he4_v41__mn45)*Y(jhe4)*Y(jv41)*state % rho - &
      screened_rates(k_he4_v41__p_cr44)*Y(jhe4)*Y(jv41)*state % rho - &
      screened_rates(k_he4_v42__mn46)*Y(jhe4)*Y(jv42)*state % rho - &
      screened_rates(k_he4_v42__p_cr45)*Y(jhe4)*Y(jv42)*state % rho - &
      screened_rates(k_he4_v43__mn47)*Y(jhe4)*Y(jv43)*state % rho - &
      screened_rates(k_he4_v43__p_cr46)*Y(jhe4)*Y(jv43)*state % rho - &
      screened_rates(k_he4_v44__mn48)*Y(jhe4)*Y(jv44)*state % rho - &
      screened_rates(k_he4_v44__p_cr47)*Y(jhe4)*Y(jv44)*state % rho - &
      screened_rates(k_he4_v45__mn49)*Y(jhe4)*Y(jv45)*state % rho - &
      screened_rates(k_he4_v45__p_cr48)*Y(jhe4)*Y(jv45)*state % rho - &
      screened_rates(k_he4_v46__mn50)*Y(jhe4)*Y(jv46)*state % rho - &
      screened_rates(k_he4_v46__p_cr49)*Y(jhe4)*Y(jv46)*state % rho - &
      screened_rates(k_he4_v47__mn51)*Y(jhe4)*Y(jv47)*state % rho - &
      screened_rates(k_he4_v47__p_cr50)*Y(jhe4)*Y(jv47)*state % rho - &
      screened_rates(k_he4_v48__mn52)*Y(jhe4)*Y(jv48)*state % rho - &
      screened_rates(k_he4_v48__p_cr51)*Y(jhe4)*Y(jv48)*state % rho - &
      screened_rates(k_he4_v49__mn53)*Y(jhe4)*Y(jv49)*state % rho - &
      screened_rates(k_he4_v49__p_cr52)*Y(jhe4)*Y(jv49)*state % rho - &
      screened_rates(k_he4_v50__mn54)*Y(jhe4)*Y(jv50)*state % rho + &
      screened_rates(k_k33__he4_cl29)*Y(jk33) + screened_rates(k_k34__he4_cl30)* &
      Y(jk34) + screened_rates(k_k35__he4_cl31)*Y(jk35) + &
      screened_rates(k_k36__he4_cl32)*Y(jk36) + &
      screened_rates(k_k36__he4_s32__weak__wc12)*Y(jk36) + &
      screened_rates(k_k37__he4_cl33)*Y(jk37) + screened_rates(k_k38__he4_cl34)* &
      Y(jk38) + screened_rates(k_k39__he4_cl35)*Y(jk39) + &
      screened_rates(k_k40__he4_cl36)*Y(jk40) + screened_rates(k_k41__he4_cl37)* &
      Y(jk41) + screened_rates(k_mg21__he4_ne17)*Y(jmg21) + &
      screened_rates(k_mg22__he4_ne18)*Y(jmg22) + screened_rates(k_mg23__he4_ne19)* &
      Y(jmg23) + screened_rates(k_mg24__he4_ne20)*Y(jmg24) + &
      screened_rates(k_mg25__he4_ne21)*Y(jmg25) + screened_rates(k_mg26__he4_ne22)* &
      Y(jmg26) + screened_rates(k_mn44__he4_v40)*Y(jmn44) + &
      screened_rates(k_mn45__he4_v41)*Y(jmn45) + screened_rates(k_mn46__he4_v42)* &
      Y(jmn46) + screened_rates(k_mn47__he4_v43)*Y(jmn47) + &
      screened_rates(k_mn48__he4_v44)*Y(jmn48) + screened_rates(k_mn49__he4_v45)* &
      Y(jmn49) + screened_rates(k_mn50__he4_v46)*Y(jmn50) + &
      screened_rates(k_mn51__he4_v47)*Y(jmn51) + screened_rates(k_mn52__he4_v48)* &
      Y(jmn52) + screened_rates(k_mn53__he4_v49)*Y(jmn53) + &
      screened_rates(k_mn54__he4_v50)*Y(jmn54) + screened_rates(k_na20__he4_f16)* &
      Y(jna20) + screened_rates(k_na20__he4_o16__weak__wc12)*Y(jna20) + &
      screened_rates(k_na21__he4_f17)*Y(jna21) + screened_rates(k_na22__he4_f18)* &
      Y(jna22) + screened_rates(k_na23__he4_f19)*Y(jna23) + &
      screened_rates(k_ne18__he4_o14)*Y(jne18) + screened_rates(k_ne19__he4_o15)* &
      Y(jne19) + screened_rates(k_ne20__he4_o16)*Y(jne20) + &
      screened_rates(k_ne21__he4_o17)*Y(jne21) + screened_rates(k_ni49__he4_fe45)* &
      Y(jni49) + screened_rates(k_ni50__he4_fe46)*Y(jni50) + &
      screened_rates(k_ni51__he4_fe47)*Y(jni51) + screened_rates(k_ni52__he4_fe48)* &
      Y(jni52) + screened_rates(k_ni53__he4_fe49)*Y(jni53) + &
      screened_rates(k_ni54__he4_fe50)*Y(jni54) + screened_rates(k_ni55__he4_fe51)* &
      Y(jni55) + screened_rates(k_ni56__he4_fe52)*Y(jni56) + &
      screened_rates(k_o16__he4_c12)*Y(jo16) + screened_rates(k_p26__he4_al22)* &
      Y(jp26) + screened_rates(k_p27__he4_al23)*Y(jp27) + &
      screened_rates(k_p28__he4_al24)*Y(jp28) + &
      screened_rates(k_p28__he4_mg24__weak__wc12)*Y(jp28) + &
      screened_rates(k_p29__he4_al25)*Y(jp29) + screened_rates(k_p30__he4_al26)* &
      Y(jp30) + screened_rates(k_p31__he4_al27)*Y(jp31) + &
      screened_rates(k_p32__he4_al28)*Y(jp32) + screened_rates(k_p_al24__he4_mg21)* &
      Y(jal24)*Y(jp)*state % rho + screened_rates(k_p_al25__he4_mg22)* &
      Y(jal25)*Y(jp)*state % rho + screened_rates(k_p_al26__he4_mg23)* &
      Y(jal26)*Y(jp)*state % rho + screened_rates(k_p_al27__he4_mg24)* &
      Y(jal27)*Y(jp)*state % rho + screened_rates(k_p_al28__he4_mg25)* &
      Y(jal28)*Y(jp)*state % rho + screened_rates(k_p_ar32__he4_cl29)* &
      Y(jar32)*Y(jp)*state % rho + screened_rates(k_p_ar33__he4_cl30)* &
      Y(jar33)*Y(jp)*state % rho + screened_rates(k_p_ar34__he4_cl31)* &
      Y(jar34)*Y(jp)*state % rho + screened_rates(k_p_ar35__he4_cl32)* &
      Y(jar35)*Y(jp)*state % rho + screened_rates(k_p_ar36__he4_cl33)* &
      Y(jar36)*Y(jp)*state % rho + screened_rates(k_p_ar37__he4_cl34)* &
      Y(jar37)*Y(jp)*state % rho + screened_rates(k_p_ar38__he4_cl35)* &
      Y(jar38)*Y(jp)*state % rho + screened_rates(k_p_ar39__he4_cl36)* &
      Y(jar39)*Y(jp)*state % rho + screened_rates(k_p_ca36__he4_k33)*Y(jca36) &
      *Y(jp)*state % rho + screened_rates(k_p_ca37__he4_k34)*Y(jca37)*Y(jp)* &
      state % rho + screened_rates(k_p_ca38__he4_k35)*Y(jca38)*Y(jp)*state % rho + &
      screened_rates(k_p_ca39__he4_k36)*Y(jca39)*Y(jp)*state % rho + &
      screened_rates(k_p_ca40__he4_k37)*Y(jca40)*Y(jp)*state % rho + &
      screened_rates(k_p_ca41__he4_k38)*Y(jca41)*Y(jp)*state % rho + &
      screened_rates(k_p_ca42__he4_k39)*Y(jca42)*Y(jp)*state % rho + &
      screened_rates(k_p_ca43__he4_k40)*Y(jca43)*Y(jp)*state % rho + &
      screened_rates(k_p_ca44__he4_k41)*Y(jca44)*Y(jp)*state % rho + &
      screened_rates(k_p_cl31__he4_s28)*Y(jcl31)*Y(jp)*state % rho + &
      screened_rates(k_p_cl32__he4_s29)*Y(jcl32)*Y(jp)*state % rho + &
      screened_rates(k_p_cl33__he4_s30)*Y(jcl33)*Y(jp)*state % rho + &
      screened_rates(k_p_cl34__he4_s31)*Y(jcl34)*Y(jp)*state % rho + &
      screened_rates(k_p_cl35__he4_s32)*Y(jcl35)*Y(jp)*state % rho + &
      screened_rates(k_p_cl36__he4_s33)*Y(jcl36)*Y(jp)*state % rho + &
      screened_rates(k_p_cl37__he4_s34)*Y(jcl37)*Y(jp)*state % rho + &
      screened_rates(k_p_co48__he4_fe45)*Y(jco48)*Y(jp)*state % rho + &
      screened_rates(k_p_co49__he4_fe46)*Y(jco49)*Y(jp)*state % rho + &
      screened_rates(k_p_co50__he4_fe47)*Y(jco50)*Y(jp)*state % rho + &
      screened_rates(k_p_co51__he4_fe48)*Y(jco51)*Y(jp)*state % rho + &
      screened_rates(k_p_co52__he4_fe49)*Y(jco52)*Y(jp)*state % rho + &
      screened_rates(k_p_co53__he4_fe50)*Y(jco53)*Y(jp)*state % rho + &
      screened_rates(k_p_co54__he4_fe51)*Y(jco54)*Y(jp)*state % rho + &
      screened_rates(k_p_co55__he4_fe52)*Y(jco55)*Y(jp)*state % rho + &
      screened_rates(k_p_co56__he4_fe53)*Y(jco56)*Y(jp)*state % rho + &
      screened_rates(k_p_cr43__he4_v40)*Y(jcr43)*Y(jp)*state % rho + &
      screened_rates(k_p_cr44__he4_v41)*Y(jcr44)*Y(jp)*state % rho + &
      screened_rates(k_p_cr45__he4_v42)*Y(jcr45)*Y(jp)*state % rho + &
      screened_rates(k_p_cr46__he4_v43)*Y(jcr46)*Y(jp)*state % rho + &
      screened_rates(k_p_cr47__he4_v44)*Y(jcr47)*Y(jp)*state % rho + &
      screened_rates(k_p_cr48__he4_v45)*Y(jcr48)*Y(jp)*state % rho + &
      screened_rates(k_p_cr49__he4_v46)*Y(jcr49)*Y(jp)*state % rho + &
      screened_rates(k_p_cr50__he4_v47)*Y(jcr50)*Y(jp)*state % rho + &
      screened_rates(k_p_cr51__he4_v48)*Y(jcr51)*Y(jp)*state % rho + &
      screened_rates(k_p_cr52__he4_v49)*Y(jcr52)*Y(jp)*state % rho + &
      screened_rates(k_p_f17__he4_o14)*Y(jf17)*Y(jp)*state % rho + &
      screened_rates(k_p_f18__he4_o15)*Y(jf18)*Y(jp)*state % rho + &
      screened_rates(k_p_f19__he4_o16)*Y(jf19)*Y(jp)*state % rho + &
      screened_rates(k_p_f20__he4_o17)*Y(jf20)*Y(jp)*state % rho + &
      screened_rates(k_p_fe47__he4_mn44)*Y(jfe47)*Y(jp)*state % rho + &
      screened_rates(k_p_fe48__he4_mn45)*Y(jfe48)*Y(jp)*state % rho + &
      screened_rates(k_p_fe49__he4_mn46)*Y(jfe49)*Y(jp)*state % rho + &
      screened_rates(k_p_fe50__he4_mn47)*Y(jfe50)*Y(jp)*state % rho + &
      screened_rates(k_p_fe51__he4_mn48)*Y(jfe51)*Y(jp)*state % rho + &
      screened_rates(k_p_fe52__he4_mn49)*Y(jfe52)*Y(jp)*state % rho + &
      screened_rates(k_p_fe53__he4_mn50)*Y(jfe53)*Y(jp)*state % rho + &
      screened_rates(k_p_fe54__he4_mn51)*Y(jfe54)*Y(jp)*state % rho + &
      screened_rates(k_p_fe55__he4_mn52)*Y(jfe55)*Y(jp)*state % rho + &
      screened_rates(k_p_fe56__he4_mn53)*Y(jfe56)*Y(jp)*state % rho + &
      screened_rates(k_p_he3__he4__weak__bet_pos_)*Y(jhe3)*Y(jp)*state % rho - &
      screened_rates(k_p_he4__d_he3)*Y(jhe4)*Y(jp)*state % rho - &
      screened_rates(k_p_he4_he4__d_be7)*Y(jhe4)**2*Y(jp)*state % rho**2 + &
      screened_rates(k_p_k34__he4_ar31)*Y(jk34)*Y(jp)*state % rho + &
      screened_rates(k_p_k35__he4_ar32)*Y(jk35)*Y(jp)*state % rho + &
      screened_rates(k_p_k36__he4_ar33)*Y(jk36)*Y(jp)*state % rho + &
      screened_rates(k_p_k37__he4_ar34)*Y(jk37)*Y(jp)*state % rho + &
      screened_rates(k_p_k38__he4_ar35)*Y(jk38)*Y(jp)*state % rho + &
      screened_rates(k_p_k39__he4_ar36)*Y(jk39)*Y(jp)*state % rho + &
      screened_rates(k_p_k40__he4_ar37)*Y(jk40)*Y(jp)*state % rho + &
      screened_rates(k_p_k41__he4_ar38)*Y(jk41)*Y(jp)*state % rho + 2.0d0* &
      screened_rates(k_p_li7__he4_he4)*Y(jli7)*Y(jp)*state % rho + &
      screened_rates(k_p_mg22__he4_na19)*Y(jmg22)*Y(jp)*state % rho + &
      screened_rates(k_p_mg23__he4_na20)*Y(jmg23)*Y(jp)*state % rho + &
      screened_rates(k_p_mg24__he4_na21)*Y(jmg24)*Y(jp)*state % rho + &
      screened_rates(k_p_mg25__he4_na22)*Y(jmg25)*Y(jp)*state % rho + &
      screened_rates(k_p_mg26__he4_na23)*Y(jmg26)*Y(jp)*state % rho + &
      screened_rates(k_p_mn45__he4_cr42)*Y(jmn45)*Y(jp)*state % rho + &
      screened_rates(k_p_mn46__he4_cr43)*Y(jmn46)*Y(jp)*state % rho + &
      screened_rates(k_p_mn47__he4_cr44)*Y(jmn47)*Y(jp)*state % rho + &
      screened_rates(k_p_mn48__he4_cr45)*Y(jmn48)*Y(jp)*state % rho + &
      screened_rates(k_p_mn49__he4_cr46)*Y(jmn49)*Y(jp)*state % rho + &
      screened_rates(k_p_mn50__he4_cr47)*Y(jmn50)*Y(jp)*state % rho + &
      screened_rates(k_p_mn51__he4_cr48)*Y(jmn51)*Y(jp)*state % rho + &
      screened_rates(k_p_mn52__he4_cr49)*Y(jmn52)*Y(jp)*state % rho + &
      screened_rates(k_p_mn53__he4_cr50)*Y(jmn53)*Y(jp)*state % rho + &
      screened_rates(k_p_mn54__he4_cr51)*Y(jmn54)*Y(jp)*state % rho + &
      screened_rates(k_p_mn55__he4_cr52)*Y(jmn55)*Y(jp)*state % rho + &
      screened_rates(k_p_n15__he4_c12)*Y(jn15)*Y(jp)*state % rho + &
      screened_rates(k_p_na20__he4_ne17)*Y(jna20)*Y(jp)*state % rho + &
      screened_rates(k_p_na21__he4_ne18)*Y(jna21)*Y(jp)*state % rho + &
      screened_rates(k_p_na22__he4_ne19)*Y(jna22)*Y(jp)*state % rho + &
      screened_rates(k_p_na23__he4_ne20)*Y(jna23)*Y(jp)*state % rho + &
      screened_rates(k_p_na24__he4_ne21)*Y(jna24)*Y(jp)*state % rho + &
      screened_rates(k_p_ne19__he4_f16)*Y(jne19)*Y(jp)*state % rho + &
      screened_rates(k_p_ne20__he4_f17)*Y(jne20)*Y(jp)*state % rho + &
      screened_rates(k_p_ne21__he4_f18)*Y(jne21)*Y(jp)*state % rho + &
      screened_rates(k_p_ne22__he4_f19)*Y(jne22)*Y(jp)*state % rho + &
      screened_rates(k_p_ni50__he4_co47)*Y(jni50)*Y(jp)*state % rho + &
      screened_rates(k_p_ni51__he4_co48)*Y(jni51)*Y(jp)*state % rho + &
      screened_rates(k_p_ni52__he4_co49)*Y(jni52)*Y(jp)*state % rho + &
      screened_rates(k_p_ni53__he4_co50)*Y(jni53)*Y(jp)*state % rho + &
      screened_rates(k_p_ni54__he4_co51)*Y(jni54)*Y(jp)*state % rho + &
      screened_rates(k_p_ni55__he4_co52)*Y(jni55)*Y(jp)*state % rho + &
      screened_rates(k_p_ni56__he4_co53)*Y(jni56)*Y(jp)*state % rho + &
      screened_rates(k_p_o15__he4_n12)*Y(jo15)*Y(jp)*state % rho + &
      screened_rates(k_p_o16__he4_n13)*Y(jo16)*Y(jp)*state % rho + &
      screened_rates(k_p_o17__he4_n14)*Y(jo17)*Y(jp)*state % rho + &
      screened_rates(k_p_o18__he4_n15)*Y(jo18)*Y(jp)*state % rho + &
      screened_rates(k_p_p27__he4_si24)*Y(jp27)*Y(jp)*state % rho + &
      screened_rates(k_p_p28__he4_si25)*Y(jp28)*Y(jp)*state % rho + &
      screened_rates(k_p_p29__he4_si26)*Y(jp29)*Y(jp)*state % rho + &
      screened_rates(k_p_p30__he4_si27)*Y(jp30)*Y(jp)*state % rho + &
      screened_rates(k_p_p31__he4_si28)*Y(jp31)*Y(jp)*state % rho + &
      screened_rates(k_p_p32__he4_si29)*Y(jp32)*Y(jp)*state % rho + &
      screened_rates(k_p_p33__he4_si30)*Y(jp33)*Y(jp)*state % rho + 0.5d0* &
      screened_rates(k_p_p_ca38__he4_ca36)*Y(jca38)*Y(jp)**2*state % rho**2 - &
      0.5d0*screened_rates(k_p_p_he4__he3_he3)*Y(jhe4)*Y(jp)**2*state % rho**2 &
      - 0.5d0*screened_rates(k_p_p_he4_he4__he3_be7)*Y(jhe4)**2*Y(jp)**2* &
      state % rho**3 + screened_rates(k_p_s29__he4_p26)*Y(jp)*Y(js29)*state % rho &
      + screened_rates(k_p_s30__he4_p27)*Y(jp)*Y(js30)*state % rho + &
      screened_rates(k_p_s31__he4_p28)*Y(jp)*Y(js31)*state % rho + &
      screened_rates(k_p_s32__he4_p29)*Y(jp)*Y(js32)*state % rho + &
      screened_rates(k_p_s33__he4_p30)*Y(jp)*Y(js33)*state % rho + &
      screened_rates(k_p_s34__he4_p31)*Y(jp)*Y(js34)*state % rho + &
      screened_rates(k_p_s35__he4_p32)*Y(jp)*Y(js35)*state % rho + &
      screened_rates(k_p_sc37__he4_ca34)*Y(jp)*Y(jsc37)*state % rho + &
      screened_rates(k_p_sc38__he4_ca35)*Y(jp)*Y(jsc38)*state % rho + &
      screened_rates(k_p_sc39__he4_ca36)*Y(jp)*Y(jsc39)*state % rho + &
      screened_rates(k_p_sc40__he4_ca37)*Y(jp)*Y(jsc40)*state % rho + &
      screened_rates(k_p_sc41__he4_ca38)*Y(jp)*Y(jsc41)*state % rho + &
      screened_rates(k_p_sc42__he4_ca39)*Y(jp)*Y(jsc42)*state % rho + &
      screened_rates(k_p_sc43__he4_ca40)*Y(jp)*Y(jsc43)*state % rho + &
      screened_rates(k_p_sc44__he4_ca41)*Y(jp)*Y(jsc44)*state % rho + &
      screened_rates(k_p_sc45__he4_ca42)*Y(jp)*Y(jsc45)*state % rho + &
      screened_rates(k_p_sc46__he4_ca43)*Y(jp)*Y(jsc46)*state % rho + &
      screened_rates(k_p_si25__he4_al22)*Y(jp)*Y(jsi25)*state % rho + &
      screened_rates(k_p_si26__he4_al23)*Y(jp)*Y(jsi26)*state % rho + &
      screened_rates(k_p_si27__he4_al24)*Y(jp)*Y(jsi27)*state % rho + &
      screened_rates(k_p_si28__he4_al25)*Y(jp)*Y(jsi28)*state % rho + &
      screened_rates(k_p_si29__he4_al26)*Y(jp)*Y(jsi29)*state % rho + &
      screened_rates(k_p_si30__he4_al27)*Y(jp)*Y(jsi30)*state % rho + &
      screened_rates(k_p_ti39__he4_sc36)*Y(jp)*Y(jti39)*state % rho + &
      screened_rates(k_p_ti40__he4_sc37)*Y(jp)*Y(jti40)*state % rho + &
      screened_rates(k_p_ti41__he4_sc38)*Y(jp)*Y(jti41)*state % rho + &
      screened_rates(k_p_ti42__he4_sc39)*Y(jp)*Y(jti42)*state % rho + &
      screened_rates(k_p_ti43__he4_sc40)*Y(jp)*Y(jti43)*state % rho + &
      screened_rates(k_p_ti44__he4_sc41)*Y(jp)*Y(jti44)*state % rho + &
      screened_rates(k_p_ti45__he4_sc42)*Y(jp)*Y(jti45)*state % rho + &
      screened_rates(k_p_ti46__he4_sc43)*Y(jp)*Y(jti46)*state % rho + &
      screened_rates(k_p_ti47__he4_sc44)*Y(jp)*Y(jti47)*state % rho + &
      screened_rates(k_p_ti48__he4_sc45)*Y(jp)*Y(jti48)*state % rho + &
      screened_rates(k_p_v41__he4_ti38)*Y(jp)*Y(jv41)*state % rho + &
      screened_rates(k_p_v42__he4_ti39)*Y(jp)*Y(jv42)*state % rho + &
      screened_rates(k_p_v43__he4_ti40)*Y(jp)*Y(jv43)*state % rho + &
      screened_rates(k_p_v44__he4_ti41)*Y(jp)*Y(jv44)*state % rho + &
      screened_rates(k_p_v45__he4_ti42)*Y(jp)*Y(jv45)*state % rho + &
      screened_rates(k_p_v46__he4_ti43)*Y(jp)*Y(jv46)*state % rho + &
      screened_rates(k_p_v47__he4_ti44)*Y(jp)*Y(jv47)*state % rho + &
      screened_rates(k_p_v48__he4_ti45)*Y(jp)*Y(jv48)*state % rho + &
      screened_rates(k_p_v49__he4_ti46)*Y(jp)*Y(jv49)*state % rho + &
      screened_rates(k_p_v50__he4_ti47)*Y(jp)*Y(jv50)*state % rho + &
      screened_rates(k_s28__he4_si24)*Y(js28) + screened_rates(k_s29__he4_si25)* &
      Y(js29) + screened_rates(k_s30__he4_si26)*Y(js30) + &
      screened_rates(k_s31__he4_si27)*Y(js31) + screened_rates(k_s32__he4_si28)* &
      Y(js32) + screened_rates(k_s33__he4_si29)*Y(js33) + &
      screened_rates(k_s34__he4_si30)*Y(js34) + screened_rates(k_sc37__he4_k33)* &
      Y(jsc37) + screened_rates(k_sc38__he4_k34)*Y(jsc38) + &
      screened_rates(k_sc39__he4_k35)*Y(jsc39) + &
      screened_rates(k_sc40__he4_ar36__weak__wc12)*Y(jsc40) + &
      screened_rates(k_sc40__he4_k36)*Y(jsc40) + screened_rates(k_sc41__he4_k37)* &
      Y(jsc41) + screened_rates(k_sc42__he4_k38)*Y(jsc42) + &
      screened_rates(k_sc43__he4_k39)*Y(jsc43) + screened_rates(k_sc44__he4_k40)* &
      Y(jsc44) + screened_rates(k_sc45__he4_k41)*Y(jsc45) + &
      screened_rates(k_si25__he4_mg21)*Y(jsi25) + screened_rates(k_si26__he4_mg22)* &
      Y(jsi26) + screened_rates(k_si27__he4_mg23)*Y(jsi27) + &
      screened_rates(k_si28__he4_mg24)*Y(jsi28) + screened_rates(k_si29__he4_mg25)* &
      Y(jsi29) + screened_rates(k_si30__he4_mg26)*Y(jsi30) + &
      screened_rates(k_ti38__he4_ca34)*Y(jti38) + screened_rates(k_ti39__he4_ca35)* &
      Y(jti39) + screened_rates(k_ti40__he4_ca36)*Y(jti40) + &
      screened_rates(k_ti41__he4_ca37)*Y(jti41) + screened_rates(k_ti42__he4_ca38)* &
      Y(jti42) + screened_rates(k_ti43__he4_ca39)*Y(jti43) + &
      screened_rates(k_ti44__he4_ca40)*Y(jti44) + screened_rates(k_ti45__he4_ca41)* &
      Y(jti45) + screened_rates(k_ti46__he4_ca42)*Y(jti46) + &
      screened_rates(k_ti47__he4_ca43)*Y(jti47) + screened_rates(k_ti48__he4_ca44)* &
      Y(jti48) + screened_rates(k_v40__he4_sc36)*Y(jv40) + &
      screened_rates(k_v41__he4_sc37)*Y(jv41) + screened_rates(k_v42__he4_sc38)* &
      Y(jv42) + screened_rates(k_v43__he4_sc39)*Y(jv43) + &
      screened_rates(k_v44__he4_sc40)*Y(jv44) + screened_rates(k_v45__he4_sc41)* &
      Y(jv45) + screened_rates(k_v46__he4_sc42)*Y(jv46) + &
      screened_rates(k_v47__he4_sc43)*Y(jv47) + screened_rates(k_v48__he4_sc44)* &
      Y(jv48) + screened_rates(k_v49__he4_sc45)*Y(jv49) + &
      screened_rates(k_v50__he4_sc46)*Y(jv50) &
       )

    ydot_nuc(jli7) = ( &
      screened_rates(k_be7__li7__weak__electron_capture)*Y(jbe7)*state % rho*state % y_e + &
      0.5d0*screened_rates(k_he4_he4__p_li7)*Y(jhe4)**2*state % rho - &
      screened_rates(k_p_li7__he4_he4)*Y(jli7)*Y(jp)*state % rho &
       )

    ydot_nuc(jbe7) = ( &
      screened_rates(k_b8__p_be7)*Y(jb8) - screened_rates(k_be7__he4_he3)*Y(jbe7) - &
      screened_rates(k_be7__li7__weak__electron_capture)*Y(jbe7)*state % rho* &
      state % y_e - screened_rates(k_d_be7__p_he4_he4)*Y(jbe7)*Y(jd)*state % rho - &
      screened_rates(k_he3_be7__p_p_he4_he4)*Y(jbe7)*Y(jhe3)*state % rho + &
      screened_rates(k_he4_he3__be7)*Y(jhe3)*Y(jhe4)*state % rho - &
      screened_rates(k_p_be7__b8)*Y(jbe7)*Y(jp)*state % rho + 0.5d0* &
      screened_rates(k_p_he4_he4__d_be7)*Y(jhe4)**2*Y(jp)*state % rho**2 + &
      0.25d0*screened_rates(k_p_p_he4_he4__he3_be7)*Y(jhe4)**2*Y(jp)**2* &
      state % rho**3 &
       )

    ydot_nuc(jbe8) = ( &
      screened_rates(k_b8__be8__weak__wc17)*Y(jb8) &
       )

    ydot_nuc(jb8) = ( &
      -screened_rates(k_b8__be8__weak__wc17)*Y(jb8) - &
      screened_rates(k_b8__he4_he4__weak__wc12)*Y(jb8) - screened_rates(k_b8__p_be7)* &
      Y(jb8) + screened_rates(k_p_be7__b8)*Y(jbe7)*Y(jp)*state % rho &
       )

    ydot_nuc(jc12) = ( &
      -screened_rates(k_c12__he4_he4_he4)*Y(jc12) - screened_rates(k_c12_c12__he4_ne20)* &
      Y(jc12)**2*state % rho - screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4) &
      *state % rho - screened_rates(k_he4_c12__p_n15)*Y(jc12)*Y(jhe4)*state % rho &
      + 0.16666666666666667d0*screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**3* &
      state % rho**2 + 2.0d0*screened_rates(k_he4_ne20__c12_c12)*Y(jhe4)* &
      Y(jne20)*state % rho + screened_rates(k_n12__c12__weak__wc12)*Y(jn12) + &
      screened_rates(k_n13__p_c12)*Y(jn13) + screened_rates(k_o16__he4_c12)*Y(jo16) - &
      screened_rates(k_p_c12__n13)*Y(jc12)*Y(jp)*state % rho + &
      screened_rates(k_p_n15__he4_c12)*Y(jn15)*Y(jp)*state % rho &
       )

    ydot_nuc(jc13) = ( &
      screened_rates(k_n13__c13__weak__wc12)*Y(jn13) + screened_rates(k_n14__p_c13)*Y(jn14) &
      - screened_rates(k_p_c13__n14)*Y(jc13)*Y(jp)*state % rho &
       )

    ydot_nuc(jn12) = ( &
      -screened_rates(k_he4_n12__p_o15)*Y(jhe4)*Y(jn12)*state % rho - &
      screened_rates(k_n12__c12__weak__wc12)*Y(jn12) + screened_rates(k_p_o15__he4_n12) &
      *Y(jo15)*Y(jp)*state % rho &
       )

    ydot_nuc(jn13) = ( &
      -screened_rates(k_he4_n13__p_o16)*Y(jhe4)*Y(jn13)*state % rho - &
      screened_rates(k_n13__c13__weak__wc12)*Y(jn13) - screened_rates(k_n13__p_c12)* &
      Y(jn13) + screened_rates(k_o14__p_n13)*Y(jo14) + screened_rates(k_p_c12__n13)* &
      Y(jc12)*Y(jp)*state % rho - screened_rates(k_p_n13__o14)*Y(jn13)* &
      Y(jp)*state % rho + screened_rates(k_p_o16__he4_n13)*Y(jo16)*Y(jp)* &
      state % rho &
       )

    ydot_nuc(jn14) = ( &
      screened_rates(k_f18__he4_n14)*Y(jf18) - screened_rates(k_he4_n14__f18)*Y(jhe4)* &
      Y(jn14)*state % rho - screened_rates(k_he4_n14__p_o17)*Y(jhe4)*Y(jn14)* &
      state % rho - screened_rates(k_n14__p_c13)*Y(jn14) + &
      screened_rates(k_o14__n14__weak__wc12)*Y(jo14) + screened_rates(k_o15__p_n14)* &
      Y(jo15) + screened_rates(k_p_c13__n14)*Y(jc13)*Y(jp)*state % rho - &
      screened_rates(k_p_n14__o15)*Y(jn14)*Y(jp)*state % rho + &
      screened_rates(k_p_o17__he4_n14)*Y(jo17)*Y(jp)*state % rho &
       )

    ydot_nuc(jn15) = ( &
      screened_rates(k_f19__he4_n15)*Y(jf19) + screened_rates(k_he4_c12__p_n15)*Y(jc12)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_n15__f19)*Y(jhe4)*Y(jn15)* &
      state % rho - screened_rates(k_he4_n15__p_o18)*Y(jhe4)*Y(jn15)*state % rho + &
      screened_rates(k_o15__n15__weak__wc12)*Y(jo15) + screened_rates(k_o16__p_n15)* &
      Y(jo16) - screened_rates(k_p_n15__he4_c12)*Y(jn15)*Y(jp)*state % rho - &
      screened_rates(k_p_n15__o16)*Y(jn15)*Y(jp)*state % rho + &
      screened_rates(k_p_o18__he4_n15)*Y(jo18)*Y(jp)*state % rho &
       )

    ydot_nuc(jo14) = ( &
      -screened_rates(k_he4_o14__ne18)*Y(jhe4)*Y(jo14)*state % rho - &
      screened_rates(k_he4_o14__p_f17)*Y(jhe4)*Y(jo14)*state % rho + &
      screened_rates(k_ne18__he4_o14)*Y(jne18) - screened_rates(k_o14__n14__weak__wc12) &
      *Y(jo14) - screened_rates(k_o14__p_n13)*Y(jo14) + &
      screened_rates(k_p_f17__he4_o14)*Y(jf17)*Y(jp)*state % rho + &
      screened_rates(k_p_n13__o14)*Y(jn13)*Y(jp)*state % rho &
       )

    ydot_nuc(jo15) = ( &
      screened_rates(k_f16__p_o15)*Y(jf16) + screened_rates(k_he4_n12__p_o15)*Y(jhe4)* &
      Y(jn12)*state % rho - screened_rates(k_he4_o15__ne19)*Y(jhe4)*Y(jo15)* &
      state % rho - screened_rates(k_he4_o15__p_f18)*Y(jhe4)*Y(jo15)*state % rho + &
      screened_rates(k_ne19__he4_o15)*Y(jne19) - screened_rates(k_o15__n15__weak__wc12) &
      *Y(jo15) - screened_rates(k_o15__p_n14)*Y(jo15) + &
      screened_rates(k_p_f18__he4_o15)*Y(jf18)*Y(jp)*state % rho + &
      screened_rates(k_p_n14__o15)*Y(jn14)*Y(jp)*state % rho - &
      screened_rates(k_p_o15__f16)*Y(jo15)*Y(jp)*state % rho - &
      screened_rates(k_p_o15__he4_n12)*Y(jo15)*Y(jp)*state % rho &
       )

    ydot_nuc(jo16) = ( &
      screened_rates(k_f17__p_o16)*Y(jf17) + screened_rates(k_he4_c12__o16)*Y(jc12)* &
      Y(jhe4)*state % rho + screened_rates(k_he4_n13__p_o16)*Y(jhe4)*Y(jn13)* &
      state % rho - screened_rates(k_he4_o16__ne20)*Y(jhe4)*Y(jo16)*state % rho - &
      screened_rates(k_he4_o16__p_f19)*Y(jhe4)*Y(jo16)*state % rho + &
      screened_rates(k_na20__he4_o16__weak__wc12)*Y(jna20) + &
      screened_rates(k_ne17__p_o16__weak__wc12)*Y(jne17) + &
      screened_rates(k_ne20__he4_o16)*Y(jne20) - screened_rates(k_o16__he4_c12)* &
      Y(jo16) - screened_rates(k_o16__p_n15)*Y(jo16) + &
      screened_rates(k_p_f19__he4_o16)*Y(jf19)*Y(jp)*state % rho + &
      screened_rates(k_p_n15__o16)*Y(jn15)*Y(jp)*state % rho - &
      screened_rates(k_p_o16__f17)*Y(jo16)*Y(jp)*state % rho - &
      screened_rates(k_p_o16__he4_n13)*Y(jo16)*Y(jp)*state % rho &
       )

    ydot_nuc(jo17) = ( &
      screened_rates(k_f17__o17__weak__wc12)*Y(jf17) + screened_rates(k_f18__p_o17)*Y(jf18) &
      + screened_rates(k_he4_n14__p_o17)*Y(jhe4)*Y(jn14)*state % rho - &
      screened_rates(k_he4_o17__ne21)*Y(jhe4)*Y(jo17)*state % rho - &
      screened_rates(k_he4_o17__p_f20)*Y(jhe4)*Y(jo17)*state % rho + &
      screened_rates(k_ne21__he4_o17)*Y(jne21) + screened_rates(k_p_f20__he4_o17)* &
      Y(jf20)*Y(jp)*state % rho - screened_rates(k_p_o17__f18)*Y(jo17)* &
      Y(jp)*state % rho - screened_rates(k_p_o17__he4_n14)*Y(jo17)*Y(jp)* &
      state % rho &
       )

    ydot_nuc(jo18) = ( &
      screened_rates(k_f18__o18__weak__wc12)*Y(jf18) + screened_rates(k_f19__p_o18)*Y(jf19) &
      + screened_rates(k_he4_n15__p_o18)*Y(jhe4)*Y(jn15)*state % rho - &
      screened_rates(k_he4_o18__ne22)*Y(jhe4)*Y(jo18)*state % rho - &
      screened_rates(k_p_o18__f19)*Y(jo18)*Y(jp)*state % rho - &
      screened_rates(k_p_o18__he4_n15)*Y(jo18)*Y(jp)*state % rho &
       )

    ydot_nuc(jf16) = ( &
      -screened_rates(k_f16__p_o15)*Y(jf16) - screened_rates(k_he4_f16__na20)*Y(jf16)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_f16__p_ne19)*Y(jf16)*Y(jhe4) &
      *state % rho + screened_rates(k_na20__he4_f16)*Y(jna20) + &
      screened_rates(k_p_ne19__he4_f16)*Y(jne19)*Y(jp)*state % rho + &
      screened_rates(k_p_o15__f16)*Y(jo15)*Y(jp)*state % rho &
       )

    ydot_nuc(jf17) = ( &
      -screened_rates(k_f17__o17__weak__wc12)*Y(jf17) - screened_rates(k_f17__p_o16)*Y(jf17) &
      - screened_rates(k_he4_f17__na21)*Y(jf17)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_f17__p_ne20)*Y(jf17)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_o14__p_f17)*Y(jhe4)*Y(jo14)*state % rho + &
      screened_rates(k_na21__he4_f17)*Y(jna21) + &
      screened_rates(k_ne17__f17__weak__wc17)*Y(jne17) + screened_rates(k_ne18__p_f17)* &
      Y(jne18) - screened_rates(k_p_f17__he4_o14)*Y(jf17)*Y(jp)*state % rho - &
      screened_rates(k_p_f17__ne18)*Y(jf17)*Y(jp)*state % rho + &
      screened_rates(k_p_ne20__he4_f17)*Y(jne20)*Y(jp)*state % rho + &
      screened_rates(k_p_o16__f17)*Y(jo16)*Y(jp)*state % rho &
       )

    ydot_nuc(jf18) = ( &
      -screened_rates(k_f18__he4_n14)*Y(jf18) - screened_rates(k_f18__o18__weak__wc12)* &
      Y(jf18) - screened_rates(k_f18__p_o17)*Y(jf18) - &
      screened_rates(k_he4_f18__na22)*Y(jf18)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_f18__p_ne21)*Y(jf18)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_n14__f18)*Y(jhe4)*Y(jn14)*state % rho + &
      screened_rates(k_he4_o15__p_f18)*Y(jhe4)*Y(jo15)*state % rho + &
      screened_rates(k_na22__he4_f18)*Y(jna22) + &
      screened_rates(k_ne18__f18__weak__wc12)*Y(jne18) + screened_rates(k_ne19__p_f18)* &
      Y(jne19) - screened_rates(k_p_f18__he4_o15)*Y(jf18)*Y(jp)*state % rho - &
      screened_rates(k_p_f18__ne19)*Y(jf18)*Y(jp)*state % rho + &
      screened_rates(k_p_ne21__he4_f18)*Y(jne21)*Y(jp)*state % rho + &
      screened_rates(k_p_o17__f18)*Y(jo17)*Y(jp)*state % rho &
       )

    ydot_nuc(jf19) = ( &
      -screened_rates(k_f19__he4_n15)*Y(jf19) - screened_rates(k_f19__p_o18)*Y(jf19) - &
      screened_rates(k_he4_f19__na23)*Y(jf19)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_f19__p_ne22)*Y(jf19)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_n15__f19)*Y(jhe4)*Y(jn15)*state % rho + &
      screened_rates(k_he4_o16__p_f19)*Y(jhe4)*Y(jo16)*state % rho + &
      screened_rates(k_na23__he4_f19)*Y(jna23) + &
      screened_rates(k_ne19__f19__weak__wc12)*Y(jne19) + screened_rates(k_ne20__p_f19)* &
      Y(jne20) - screened_rates(k_p_f19__he4_o16)*Y(jf19)*Y(jp)*state % rho - &
      screened_rates(k_p_f19__ne20)*Y(jf19)*Y(jp)*state % rho + &
      screened_rates(k_p_ne22__he4_f19)*Y(jne22)*Y(jp)*state % rho + &
      screened_rates(k_p_o18__f19)*Y(jo18)*Y(jp)*state % rho &
       )

    ydot_nuc(jf20) = ( &
      -screened_rates(k_f20__ne20__weak__wc12)*Y(jf20) - screened_rates(k_he4_f20__na24)* &
      Y(jf20)*Y(jhe4)*state % rho + screened_rates(k_he4_o17__p_f20)*Y(jhe4)* &
      Y(jo17)*state % rho + screened_rates(k_ne21__p_f20)*Y(jne21) - &
      screened_rates(k_p_f20__he4_o17)*Y(jf20)*Y(jp)*state % rho - &
      screened_rates(k_p_f20__ne21)*Y(jf20)*Y(jp)*state % rho &
       )

    ydot_nuc(jne17) = ( &
      -screened_rates(k_he4_ne17__mg21)*Y(jhe4)*Y(jne17)*state % rho - &
      screened_rates(k_he4_ne17__p_na20)*Y(jhe4)*Y(jne17)*state % rho + &
      screened_rates(k_mg21__he4_ne17)*Y(jmg21) - &
      screened_rates(k_ne17__f17__weak__wc17)*Y(jne17) - &
      screened_rates(k_ne17__p_o16__weak__wc12)*Y(jne17) + &
      screened_rates(k_p_na20__he4_ne17)*Y(jna20)*Y(jp)*state % rho &
       )

    ydot_nuc(jne18) = ( &
      screened_rates(k_al22__he4_ne18__weak__wc12)*Y(jal22) - screened_rates(k_he4_ne18__mg22) &
      *Y(jhe4)*Y(jne18)*state % rho - screened_rates(k_he4_ne18__p_na21)* &
      Y(jhe4)*Y(jne18)*state % rho + screened_rates(k_he4_o14__ne18)*Y(jhe4)* &
      Y(jo14)*state % rho + screened_rates(k_mg22__he4_ne18)*Y(jmg22) + &
      screened_rates(k_na19__p_ne18)*Y(jna19) - screened_rates(k_ne18__f18__weak__wc12) &
      *Y(jne18) - screened_rates(k_ne18__he4_o14)*Y(jne18) - &
      screened_rates(k_ne18__p_f17)*Y(jne18) + screened_rates(k_p_f17__ne18)*Y(jf17)* &
      Y(jp)*state % rho + screened_rates(k_p_na21__he4_ne18)*Y(jna21)*Y(jp)* &
      state % rho - screened_rates(k_p_ne18__na19)*Y(jne18)*Y(jp)*state % rho &
       )

    ydot_nuc(jne19) = ( &
      screened_rates(k_he4_f16__p_ne19)*Y(jf16)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ne19__mg23)*Y(jhe4)*Y(jne19)*state % rho - &
      screened_rates(k_he4_ne19__p_na22)*Y(jhe4)*Y(jne19)*state % rho + &
      screened_rates(k_he4_o15__ne19)*Y(jhe4)*Y(jo15)*state % rho + &
      screened_rates(k_mg23__he4_ne19)*Y(jmg23) + &
      screened_rates(k_na19__ne19__weak__bqa_pos_)*Y(jna19) + &
      screened_rates(k_na20__p_ne19)*Y(jna20) - screened_rates(k_ne19__f19__weak__wc12) &
      *Y(jne19) - screened_rates(k_ne19__he4_o15)*Y(jne19) - &
      screened_rates(k_ne19__p_f18)*Y(jne19) + screened_rates(k_p_f18__ne19)*Y(jf18)* &
      Y(jp)*state % rho + screened_rates(k_p_na22__he4_ne19)*Y(jna22)*Y(jp)* &
      state % rho - screened_rates(k_p_ne19__he4_f16)*Y(jne19)*Y(jp)*state % rho - &
      screened_rates(k_p_ne19__na20)*Y(jne19)*Y(jp)*state % rho &
       )

    ydot_nuc(jne20) = ( &
      screened_rates(k_al22__p_p_ne20__weak__wc12)*Y(jal22) + &
      screened_rates(k_al24__he4_ne20__weak__wc12)*Y(jal24) + 0.5d0* &
      screened_rates(k_c12_c12__he4_ne20)*Y(jc12)**2*state % rho + &
      screened_rates(k_f20__ne20__weak__wc12)*Y(jf20) + &
      screened_rates(k_he4_f17__p_ne20)*Y(jf17)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ne20__c12_c12)*Y(jhe4)*Y(jne20)*state % rho - &
      screened_rates(k_he4_ne20__mg24)*Y(jhe4)*Y(jne20)*state % rho - &
      screened_rates(k_he4_ne20__p_na23)*Y(jhe4)*Y(jne20)*state % rho + &
      screened_rates(k_he4_o16__ne20)*Y(jhe4)*Y(jo16)*state % rho + &
      screened_rates(k_mg21__p_ne20__weak__wc12)*Y(jmg21) + &
      screened_rates(k_mg24__he4_ne20)*Y(jmg24) + &
      screened_rates(k_na20__ne20__weak__wc12)*Y(jna20) + &
      screened_rates(k_na21__p_ne20)*Y(jna21) - screened_rates(k_ne20__he4_o16)* &
      Y(jne20) - screened_rates(k_ne20__p_f19)*Y(jne20) + &
      screened_rates(k_p_f19__ne20)*Y(jf19)*Y(jp)*state % rho + &
      screened_rates(k_p_na23__he4_ne20)*Y(jna23)*Y(jp)*state % rho - &
      screened_rates(k_p_ne20__he4_f17)*Y(jne20)*Y(jp)*state % rho - &
      screened_rates(k_p_ne20__na21)*Y(jne20)*Y(jp)*state % rho &
       )

    ydot_nuc(jne21) = ( &
      screened_rates(k_he4_f18__p_ne21)*Y(jf18)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ne21__mg25)*Y(jhe4)*Y(jne21)*state % rho - &
      screened_rates(k_he4_ne21__p_na24)*Y(jhe4)*Y(jne21)*state % rho + &
      screened_rates(k_he4_o17__ne21)*Y(jhe4)*Y(jo17)*state % rho + &
      screened_rates(k_mg25__he4_ne21)*Y(jmg25) + &
      screened_rates(k_na21__ne21__weak__wc12)*Y(jna21) + &
      screened_rates(k_na22__p_ne21)*Y(jna22) - screened_rates(k_ne21__he4_o17)* &
      Y(jne21) - screened_rates(k_ne21__p_f20)*Y(jne21) + &
      screened_rates(k_p_f20__ne21)*Y(jf20)*Y(jp)*state % rho + &
      screened_rates(k_p_na24__he4_ne21)*Y(jna24)*Y(jp)*state % rho - &
      screened_rates(k_p_ne21__he4_f18)*Y(jne21)*Y(jp)*state % rho - &
      screened_rates(k_p_ne21__na22)*Y(jne21)*Y(jp)*state % rho &
       )

    ydot_nuc(jne22) = ( &
      screened_rates(k_he4_f19__p_ne22)*Y(jf19)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ne22__mg26)*Y(jhe4)*Y(jne22)*state % rho + &
      screened_rates(k_he4_o18__ne22)*Y(jhe4)*Y(jo18)*state % rho + &
      screened_rates(k_mg26__he4_ne22)*Y(jmg26) + &
      screened_rates(k_na22__ne22__weak__wc12)*Y(jna22) + &
      screened_rates(k_na23__p_ne22)*Y(jna23) - screened_rates(k_p_ne22__he4_f19)* &
      Y(jne22)*Y(jp)*state % rho - screened_rates(k_p_ne22__na23)*Y(jne22)* &
      Y(jp)*state % rho &
       )

    ydot_nuc(jna19) = ( &
      screened_rates(k_al23__he4_na19)*Y(jal23) - screened_rates(k_he4_na19__al23)*Y(jhe4)* &
      Y(jna19)*state % rho - screened_rates(k_he4_na19__p_mg22)*Y(jhe4)* &
      Y(jna19)*state % rho - screened_rates(k_na19__ne19__weak__bqa_pos_)* &
      Y(jna19) - screened_rates(k_na19__p_ne18)*Y(jna19) + &
      screened_rates(k_p_mg22__he4_na19)*Y(jmg22)*Y(jp)*state % rho + &
      screened_rates(k_p_ne18__na19)*Y(jne18)*Y(jp)*state % rho &
       )

    ydot_nuc(jna20) = ( &
      screened_rates(k_al24__he4_na20)*Y(jal24) + screened_rates(k_he4_f16__na20)*Y(jf16)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_na20__al24)*Y(jhe4)*Y(jna20) &
      *state % rho - screened_rates(k_he4_na20__p_mg23)*Y(jhe4)*Y(jna20)* &
      state % rho + screened_rates(k_he4_ne17__p_na20)*Y(jhe4)*Y(jne17)* &
      state % rho + screened_rates(k_mg21__p_na20)*Y(jmg21) - &
      screened_rates(k_na20__he4_f16)*Y(jna20) - &
      screened_rates(k_na20__he4_o16__weak__wc12)*Y(jna20) - &
      screened_rates(k_na20__ne20__weak__wc12)*Y(jna20) - &
      screened_rates(k_na20__p_ne19)*Y(jna20) + screened_rates(k_p_mg23__he4_na20)* &
      Y(jmg23)*Y(jp)*state % rho - screened_rates(k_p_na20__he4_ne17)* &
      Y(jna20)*Y(jp)*state % rho - screened_rates(k_p_na20__mg21)*Y(jna20)* &
      Y(jp)*state % rho + screened_rates(k_p_ne19__na20)*Y(jne19)*Y(jp)* &
      state % rho &
       )

    ydot_nuc(jna21) = ( &
      screened_rates(k_al22__p_na21__weak__wc17)*Y(jal22) + screened_rates(k_al25__he4_na21)* &
      Y(jal25) + screened_rates(k_he4_f17__na21)*Y(jf17)*Y(jhe4)*state % rho &
      - screened_rates(k_he4_na21__al25)*Y(jhe4)*Y(jna21)*state % rho - &
      screened_rates(k_he4_na21__p_mg24)*Y(jhe4)*Y(jna21)*state % rho + &
      screened_rates(k_he4_ne18__p_na21)*Y(jhe4)*Y(jne18)*state % rho + &
      screened_rates(k_mg21__na21__weak__wc12)*Y(jmg21) + &
      screened_rates(k_mg22__p_na21)*Y(jmg22) - screened_rates(k_na21__he4_f17)* &
      Y(jna21) - screened_rates(k_na21__ne21__weak__wc12)*Y(jna21) - &
      screened_rates(k_na21__p_ne20)*Y(jna21) + screened_rates(k_p_mg24__he4_na21)* &
      Y(jmg24)*Y(jp)*state % rho - screened_rates(k_p_na21__he4_ne18)* &
      Y(jna21)*Y(jp)*state % rho - screened_rates(k_p_na21__mg22)*Y(jna21)* &
      Y(jp)*state % rho + screened_rates(k_p_ne20__na21)*Y(jne20)*Y(jp)* &
      state % rho &
       )

    ydot_nuc(jna22) = ( &
      screened_rates(k_al23__p_na22__weak__wc12)*Y(jal23) + screened_rates(k_al26__he4_na22)* &
      Y(jal26) + screened_rates(k_he4_f18__na22)*Y(jf18)*Y(jhe4)*state % rho &
      - screened_rates(k_he4_na22__al26)*Y(jhe4)*Y(jna22)*state % rho - &
      screened_rates(k_he4_na22__p_mg25)*Y(jhe4)*Y(jna22)*state % rho + &
      screened_rates(k_he4_ne19__p_na22)*Y(jhe4)*Y(jne19)*state % rho + &
      screened_rates(k_mg22__na22__weak__wc12)*Y(jmg22) + &
      screened_rates(k_mg23__p_na22)*Y(jmg23) - screened_rates(k_na22__he4_f18)* &
      Y(jna22) - screened_rates(k_na22__ne22__weak__wc12)*Y(jna22) - &
      screened_rates(k_na22__p_ne21)*Y(jna22) + screened_rates(k_p_mg25__he4_na22)* &
      Y(jmg25)*Y(jp)*state % rho - screened_rates(k_p_na22__he4_ne19)* &
      Y(jna22)*Y(jp)*state % rho - screened_rates(k_p_na22__mg23)*Y(jna22)* &
      Y(jp)*state % rho + screened_rates(k_p_ne21__na22)*Y(jne21)*Y(jp)* &
      state % rho &
       )

    ydot_nuc(jna23) = ( &
      screened_rates(k_al24__p_na23__weak__wc12)*Y(jal24) + screened_rates(k_al27__he4_na23)* &
      Y(jal27) + screened_rates(k_he4_f19__na23)*Y(jf19)*Y(jhe4)*state % rho &
      - screened_rates(k_he4_na23__al27)*Y(jhe4)*Y(jna23)*state % rho - &
      screened_rates(k_he4_na23__p_mg26)*Y(jhe4)*Y(jna23)*state % rho + &
      screened_rates(k_he4_ne20__p_na23)*Y(jhe4)*Y(jne20)*state % rho + &
      screened_rates(k_mg23__na23__weak__wc12)*Y(jmg23) + &
      screened_rates(k_mg24__p_na23)*Y(jmg24) - screened_rates(k_na23__he4_f19)* &
      Y(jna23) - screened_rates(k_na23__p_ne22)*Y(jna23) + &
      screened_rates(k_p_mg26__he4_na23)*Y(jmg26)*Y(jp)*state % rho - &
      screened_rates(k_p_na23__he4_ne20)*Y(jna23)*Y(jp)*state % rho - &
      screened_rates(k_p_na23__mg24)*Y(jna23)*Y(jp)*state % rho + &
      screened_rates(k_p_ne22__na23)*Y(jne22)*Y(jp)*state % rho &
       )

    ydot_nuc(jna24) = ( &
      screened_rates(k_al28__he4_na24)*Y(jal28) + screened_rates(k_he4_f20__na24)*Y(jf20)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_na24__al28)*Y(jhe4)*Y(jna24) &
      *state % rho + screened_rates(k_he4_ne21__p_na24)*Y(jhe4)*Y(jne21)* &
      state % rho + screened_rates(k_mg25__p_na24)*Y(jmg25) - &
      screened_rates(k_p_na24__he4_ne21)*Y(jna24)*Y(jp)*state % rho - &
      screened_rates(k_p_na24__mg25)*Y(jna24)*Y(jp)*state % rho &
       )

    ydot_nuc(jmg21) = ( &
      screened_rates(k_al22__p_mg21)*Y(jal22) - screened_rates(k_he4_mg21__p_al24)*Y(jhe4)* &
      Y(jmg21)*state % rho - screened_rates(k_he4_mg21__si25)*Y(jhe4)* &
      Y(jmg21)*state % rho + screened_rates(k_he4_ne17__mg21)*Y(jhe4)* &
      Y(jne17)*state % rho - screened_rates(k_mg21__he4_ne17)*Y(jmg21) - &
      screened_rates(k_mg21__na21__weak__wc12)*Y(jmg21) - &
      screened_rates(k_mg21__p_na20)*Y(jmg21) - &
      screened_rates(k_mg21__p_ne20__weak__wc12)*Y(jmg21) + &
      screened_rates(k_p_al24__he4_mg21)*Y(jal24)*Y(jp)*state % rho - &
      screened_rates(k_p_mg21__al22)*Y(jmg21)*Y(jp)*state % rho + &
      screened_rates(k_p_na20__mg21)*Y(jna20)*Y(jp)*state % rho + &
      screened_rates(k_si25__he4_mg21)*Y(jsi25) &
       )

    ydot_nuc(jmg22) = ( &
      screened_rates(k_al22__mg22__weak__wc12)*Y(jal22) + screened_rates(k_al23__p_mg22)* &
      Y(jal23) - screened_rates(k_he4_mg22__p_al25)*Y(jhe4)*Y(jmg22)* &
      state % rho - screened_rates(k_he4_mg22__si26)*Y(jhe4)*Y(jmg22)*state % rho &
      + screened_rates(k_he4_na19__p_mg22)*Y(jhe4)*Y(jna19)*state % rho + &
      screened_rates(k_he4_ne18__mg22)*Y(jhe4)*Y(jne18)*state % rho - &
      screened_rates(k_mg22__he4_ne18)*Y(jmg22) - &
      screened_rates(k_mg22__na22__weak__wc12)*Y(jmg22) - &
      screened_rates(k_mg22__p_na21)*Y(jmg22) + screened_rates(k_p_al25__he4_mg22)* &
      Y(jal25)*Y(jp)*state % rho - screened_rates(k_p_mg22__al23)*Y(jmg22)* &
      Y(jp)*state % rho - screened_rates(k_p_mg22__he4_na19)*Y(jmg22)*Y(jp)* &
      state % rho + screened_rates(k_p_na21__mg22)*Y(jna21)*Y(jp)*state % rho + &
      screened_rates(k_si26__he4_mg22)*Y(jsi26) &
       )

    ydot_nuc(jmg23) = ( &
      screened_rates(k_al23__mg23__weak__wc12)*Y(jal23) + screened_rates(k_al24__p_mg23)* &
      Y(jal24) - screened_rates(k_he4_mg23__p_al26)*Y(jhe4)*Y(jmg23)* &
      state % rho - screened_rates(k_he4_mg23__si27)*Y(jhe4)*Y(jmg23)*state % rho &
      + screened_rates(k_he4_na20__p_mg23)*Y(jhe4)*Y(jna20)*state % rho + &
      screened_rates(k_he4_ne19__mg23)*Y(jhe4)*Y(jne19)*state % rho - &
      screened_rates(k_mg23__he4_ne19)*Y(jmg23) - &
      screened_rates(k_mg23__na23__weak__wc12)*Y(jmg23) - &
      screened_rates(k_mg23__p_na22)*Y(jmg23) + screened_rates(k_p_al26__he4_mg23)* &
      Y(jal26)*Y(jp)*state % rho - screened_rates(k_p_mg23__al24)*Y(jmg23)* &
      Y(jp)*state % rho - screened_rates(k_p_mg23__he4_na20)*Y(jmg23)*Y(jp)* &
      state % rho + screened_rates(k_p_na22__mg23)*Y(jna22)*Y(jp)*state % rho + &
      screened_rates(k_si24__p_mg23__weak__wc12)*Y(jsi24) + &
      screened_rates(k_si27__he4_mg23)*Y(jsi27) &
       )

    ydot_nuc(jmg24) = ( &
      screened_rates(k_al24__mg24__weak__wc12)*Y(jal24) + screened_rates(k_al25__p_mg24)* &
      Y(jal25) - screened_rates(k_he4_mg24__p_al27)*Y(jhe4)*Y(jmg24)* &
      state % rho - screened_rates(k_he4_mg24__si28)*Y(jhe4)*Y(jmg24)*state % rho &
      + screened_rates(k_he4_na21__p_mg24)*Y(jhe4)*Y(jna21)*state % rho + &
      screened_rates(k_he4_ne20__mg24)*Y(jhe4)*Y(jne20)*state % rho - &
      screened_rates(k_mg24__he4_ne20)*Y(jmg24) - screened_rates(k_mg24__p_na23)* &
      Y(jmg24) + screened_rates(k_p28__he4_mg24__weak__wc12)*Y(jp28) + &
      screened_rates(k_p_al27__he4_mg24)*Y(jal27)*Y(jp)*state % rho - &
      screened_rates(k_p_mg24__al25)*Y(jmg24)*Y(jp)*state % rho - &
      screened_rates(k_p_mg24__he4_na21)*Y(jmg24)*Y(jp)*state % rho + &
      screened_rates(k_p_na23__mg24)*Y(jna23)*Y(jp)*state % rho + &
      screened_rates(k_si25__p_mg24__weak__wc12)*Y(jsi25) + &
      screened_rates(k_si28__he4_mg24)*Y(jsi28) &
       )

    ydot_nuc(jmg25) = ( &
      screened_rates(k_al25__mg25__weak__wc12)*Y(jal25) + screened_rates(k_al26__p_mg25)* &
      Y(jal26) - screened_rates(k_he4_mg25__p_al28)*Y(jhe4)*Y(jmg25)* &
      state % rho - screened_rates(k_he4_mg25__si29)*Y(jhe4)*Y(jmg25)*state % rho &
      + screened_rates(k_he4_na22__p_mg25)*Y(jhe4)*Y(jna22)*state % rho + &
      screened_rates(k_he4_ne21__mg25)*Y(jhe4)*Y(jne21)*state % rho - &
      screened_rates(k_mg25__he4_ne21)*Y(jmg25) - screened_rates(k_mg25__p_na24)* &
      Y(jmg25) + screened_rates(k_p_al28__he4_mg25)*Y(jal28)*Y(jp)* &
      state % rho - screened_rates(k_p_mg25__al26)*Y(jmg25)*Y(jp)*state % rho - &
      screened_rates(k_p_mg25__he4_na22)*Y(jmg25)*Y(jp)*state % rho + &
      screened_rates(k_p_na24__mg25)*Y(jna24)*Y(jp)*state % rho + &
      screened_rates(k_si29__he4_mg25)*Y(jsi29) &
       )

    ydot_nuc(jmg26) = ( &
      screened_rates(k_al26__mg26__weak__wc12)*Y(jal26) + screened_rates(k_al27__p_mg26)* &
      Y(jal27) - screened_rates(k_he4_mg26__si30)*Y(jhe4)*Y(jmg26)* &
      state % rho + screened_rates(k_he4_na23__p_mg26)*Y(jhe4)*Y(jna23)* &
      state % rho + screened_rates(k_he4_ne22__mg26)*Y(jhe4)*Y(jne22)*state % rho &
      - screened_rates(k_mg26__he4_ne22)*Y(jmg26) - screened_rates(k_p_mg26__al27)* &
      Y(jmg26)*Y(jp)*state % rho - screened_rates(k_p_mg26__he4_na23)* &
      Y(jmg26)*Y(jp)*state % rho + screened_rates(k_si30__he4_mg26)*Y(jsi30) &
       )

    ydot_nuc(jal22) = ( &
      -screened_rates(k_al22__he4_ne18__weak__wc12)*Y(jal22) - &
      screened_rates(k_al22__mg22__weak__wc12)*Y(jal22) - &
      screened_rates(k_al22__p_mg21)*Y(jal22) - &
      screened_rates(k_al22__p_na21__weak__wc17)*Y(jal22) - &
      screened_rates(k_al22__p_p_ne20__weak__wc12)*Y(jal22) - &
      screened_rates(k_he4_al22__p26)*Y(jal22)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al22__p_si25)*Y(jal22)*Y(jhe4)*state % rho + &
      screened_rates(k_p26__he4_al22)*Y(jp26) + screened_rates(k_p_mg21__al22)* &
      Y(jmg21)*Y(jp)*state % rho + screened_rates(k_p_si25__he4_al22)*Y(jp)* &
      Y(jsi25)*state % rho &
       )

    ydot_nuc(jal23) = ( &
      -screened_rates(k_al23__he4_na19)*Y(jal23) - screened_rates(k_al23__mg23__weak__wc12)* &
      Y(jal23) - screened_rates(k_al23__p_mg22)*Y(jal23) - &
      screened_rates(k_al23__p_na22__weak__wc12)*Y(jal23) - &
      screened_rates(k_he4_al23__p27)*Y(jal23)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al23__p_si26)*Y(jal23)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_na19__al23)*Y(jhe4)*Y(jna19)*state % rho + &
      screened_rates(k_p27__he4_al23)*Y(jp27) - screened_rates(k_p_al23__si24)* &
      Y(jal23)*Y(jp)*state % rho + screened_rates(k_p_mg22__al23)*Y(jmg22)* &
      Y(jp)*state % rho + screened_rates(k_p_si26__he4_al23)*Y(jp)*Y(jsi26)* &
      state % rho + screened_rates(k_si24__p_al23)*Y(jsi24) &
       )

    ydot_nuc(jal24) = ( &
      -screened_rates(k_al24__he4_na20)*Y(jal24) - &
      screened_rates(k_al24__he4_ne20__weak__wc12)*Y(jal24) - &
      screened_rates(k_al24__mg24__weak__wc12)*Y(jal24) - &
      screened_rates(k_al24__p_mg23)*Y(jal24) - &
      screened_rates(k_al24__p_na23__weak__wc12)*Y(jal24) - &
      screened_rates(k_he4_al24__p28)*Y(jal24)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al24__p_si27)*Y(jal24)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mg21__p_al24)*Y(jhe4)*Y(jmg21)*state % rho + &
      screened_rates(k_he4_na20__al24)*Y(jhe4)*Y(jna20)*state % rho + &
      screened_rates(k_p28__he4_al24)*Y(jp28) - screened_rates(k_p_al24__he4_mg21)* &
      Y(jal24)*Y(jp)*state % rho - screened_rates(k_p_al24__si25)*Y(jal24)* &
      Y(jp)*state % rho + screened_rates(k_p_mg23__al24)*Y(jmg23)*Y(jp)* &
      state % rho + screened_rates(k_p_si27__he4_al24)*Y(jp)*Y(jsi27)*state % rho &
      + screened_rates(k_si24__al24__weak__wc12)*Y(jsi24) + &
      screened_rates(k_si25__p_al24)*Y(jsi25) &
       )

    ydot_nuc(jal25) = ( &
      -screened_rates(k_al25__he4_na21)*Y(jal25) - screened_rates(k_al25__mg25__weak__wc12)* &
      Y(jal25) - screened_rates(k_al25__p_mg24)*Y(jal25) - &
      screened_rates(k_he4_al25__p29)*Y(jal25)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al25__p_si28)*Y(jal25)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mg22__p_al25)*Y(jhe4)*Y(jmg22)*state % rho + &
      screened_rates(k_he4_na21__al25)*Y(jhe4)*Y(jna21)*state % rho + &
      screened_rates(k_p29__he4_al25)*Y(jp29) - screened_rates(k_p_al25__he4_mg22)* &
      Y(jal25)*Y(jp)*state % rho - screened_rates(k_p_al25__si26)*Y(jal25)* &
      Y(jp)*state % rho + screened_rates(k_p_mg24__al25)*Y(jmg24)*Y(jp)* &
      state % rho + screened_rates(k_p_si28__he4_al25)*Y(jp)*Y(jsi28)*state % rho &
      + screened_rates(k_si25__al25__weak__wc12)*Y(jsi25) + &
      screened_rates(k_si26__p_al25)*Y(jsi26) &
       )

    ydot_nuc(jal26) = ( &
      -screened_rates(k_al26__he4_na22)*Y(jal26) - screened_rates(k_al26__mg26__weak__wc12)* &
      Y(jal26) - screened_rates(k_al26__p_mg25)*Y(jal26) - &
      screened_rates(k_he4_al26__p30)*Y(jal26)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al26__p_si29)*Y(jal26)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mg23__p_al26)*Y(jhe4)*Y(jmg23)*state % rho + &
      screened_rates(k_he4_na22__al26)*Y(jhe4)*Y(jna22)*state % rho + &
      screened_rates(k_p27__p_al26__weak__wc12)*Y(jp27) + &
      screened_rates(k_p30__he4_al26)*Y(jp30) - screened_rates(k_p_al26__he4_mg23)* &
      Y(jal26)*Y(jp)*state % rho - screened_rates(k_p_al26__si27)*Y(jal26)* &
      Y(jp)*state % rho + screened_rates(k_p_mg25__al26)*Y(jmg25)*Y(jp)* &
      state % rho + screened_rates(k_p_si29__he4_al26)*Y(jp)*Y(jsi29)*state % rho &
      + screened_rates(k_si26__al26__weak__wc12)*Y(jsi26) + &
      screened_rates(k_si27__p_al26)*Y(jsi27) &
       )

    ydot_nuc(jal27) = ( &
      -screened_rates(k_al27__he4_na23)*Y(jal27) - screened_rates(k_al27__p_mg26)*Y(jal27) - &
      screened_rates(k_he4_al27__p31)*Y(jal27)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al27__p_si30)*Y(jal27)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mg24__p_al27)*Y(jhe4)*Y(jmg24)*state % rho + &
      screened_rates(k_he4_na23__al27)*Y(jhe4)*Y(jna23)*state % rho + &
      screened_rates(k_p28__p_al27__weak__wc12)*Y(jp28) + &
      screened_rates(k_p31__he4_al27)*Y(jp31) - screened_rates(k_p_al27__he4_mg24)* &
      Y(jal27)*Y(jp)*state % rho - screened_rates(k_p_al27__si28)*Y(jal27)* &
      Y(jp)*state % rho + screened_rates(k_p_mg26__al27)*Y(jmg26)*Y(jp)* &
      state % rho + screened_rates(k_p_si30__he4_al27)*Y(jp)*Y(jsi30)*state % rho &
      + screened_rates(k_si27__al27__weak__wc12)*Y(jsi27) + &
      screened_rates(k_si28__p_al27)*Y(jsi28) &
       )

    ydot_nuc(jal28) = ( &
      -screened_rates(k_al28__he4_na24)*Y(jal28) - screened_rates(k_he4_al28__p32)*Y(jal28)* &
      Y(jhe4)*state % rho + screened_rates(k_he4_mg25__p_al28)*Y(jhe4)* &
      Y(jmg25)*state % rho + screened_rates(k_he4_na24__al28)*Y(jhe4)* &
      Y(jna24)*state % rho + screened_rates(k_p32__he4_al28)*Y(jp32) - &
      screened_rates(k_p_al28__he4_mg25)*Y(jal28)*Y(jp)*state % rho - &
      screened_rates(k_p_al28__si29)*Y(jal28)*Y(jp)*state % rho + &
      screened_rates(k_si29__p_al28)*Y(jsi29) &
       )

    ydot_nuc(jsi24) = ( &
      -screened_rates(k_he4_si24__p_p27)*Y(jhe4)*Y(jsi24)*state % rho - &
      screened_rates(k_he4_si24__s28)*Y(jhe4)*Y(jsi24)*state % rho + &
      screened_rates(k_p_al23__si24)*Y(jal23)*Y(jp)*state % rho + &
      screened_rates(k_p_p27__he4_si24)*Y(jp27)*Y(jp)*state % rho + &
      screened_rates(k_s28__he4_si24)*Y(js28) - &
      screened_rates(k_si24__al24__weak__wc12)*Y(jsi24) - &
      screened_rates(k_si24__p_al23)*Y(jsi24) - &
      screened_rates(k_si24__p_mg23__weak__wc12)*Y(jsi24) &
       )

    ydot_nuc(jsi25) = ( &
      screened_rates(k_he4_al22__p_si25)*Y(jal22)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mg21__si25)*Y(jhe4)*Y(jmg21)*state % rho - &
      screened_rates(k_he4_si25__p_p28)*Y(jhe4)*Y(jsi25)*state % rho - &
      screened_rates(k_he4_si25__s29)*Y(jhe4)*Y(jsi25)*state % rho + &
      screened_rates(k_p26__p_si25)*Y(jp26) + screened_rates(k_p_al24__si25)*Y(jal24) &
      *Y(jp)*state % rho + screened_rates(k_p_p28__he4_si25)*Y(jp28)*Y(jp)* &
      state % rho - screened_rates(k_p_si25__he4_al22)*Y(jp)*Y(jsi25)*state % rho &
      - screened_rates(k_p_si25__p26)*Y(jp)*Y(jsi25)*state % rho + &
      screened_rates(k_s29__he4_si25)*Y(js29) - &
      screened_rates(k_si25__al25__weak__wc12)*Y(jsi25) - &
      screened_rates(k_si25__he4_mg21)*Y(jsi25) - screened_rates(k_si25__p_al24)* &
      Y(jsi25) - screened_rates(k_si25__p_mg24__weak__wc12)*Y(jsi25) &
       )

    ydot_nuc(jsi26) = ( &
      screened_rates(k_he4_al23__p_si26)*Y(jal23)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mg22__si26)*Y(jhe4)*Y(jmg22)*state % rho - &
      screened_rates(k_he4_si26__p_p29)*Y(jhe4)*Y(jsi26)*state % rho - &
      screened_rates(k_he4_si26__s30)*Y(jhe4)*Y(jsi26)*state % rho + &
      screened_rates(k_p26__si26__weak__wc12)*Y(jp26) + screened_rates(k_p27__p_si26)* &
      Y(jp27) + screened_rates(k_p_al25__si26)*Y(jal25)*Y(jp)*state % rho + &
      screened_rates(k_p_p29__he4_si26)*Y(jp29)*Y(jp)*state % rho - &
      screened_rates(k_p_si26__he4_al23)*Y(jp)*Y(jsi26)*state % rho - &
      screened_rates(k_p_si26__p27)*Y(jp)*Y(jsi26)*state % rho + &
      screened_rates(k_s30__he4_si26)*Y(js30) - &
      screened_rates(k_si26__al26__weak__wc12)*Y(jsi26) - &
      screened_rates(k_si26__he4_mg22)*Y(jsi26) - screened_rates(k_si26__p_al25)* &
      Y(jsi26) &
       )

    ydot_nuc(jsi27) = ( &
      screened_rates(k_he4_al24__p_si27)*Y(jal24)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mg23__si27)*Y(jhe4)*Y(jmg23)*state % rho - &
      screened_rates(k_he4_si27__p_p30)*Y(jhe4)*Y(jsi27)*state % rho - &
      screened_rates(k_he4_si27__s31)*Y(jhe4)*Y(jsi27)*state % rho + &
      screened_rates(k_p27__si27__weak__wc12)*Y(jp27) + screened_rates(k_p28__p_si27)* &
      Y(jp28) + screened_rates(k_p_al26__si27)*Y(jal26)*Y(jp)*state % rho + &
      screened_rates(k_p_p30__he4_si27)*Y(jp30)*Y(jp)*state % rho - &
      screened_rates(k_p_si27__he4_al24)*Y(jp)*Y(jsi27)*state % rho - &
      screened_rates(k_p_si27__p28)*Y(jp)*Y(jsi27)*state % rho + &
      screened_rates(k_s28__p_si27__weak__wc12)*Y(js28) + &
      screened_rates(k_s31__he4_si27)*Y(js31) - &
      screened_rates(k_si27__al27__weak__wc12)*Y(jsi27) - &
      screened_rates(k_si27__he4_mg23)*Y(jsi27) - screened_rates(k_si27__p_al26)* &
      Y(jsi27) &
       )

    ydot_nuc(jsi28) = ( &
      screened_rates(k_cl32__he4_si28__weak__wc12)*Y(jcl32) + &
      screened_rates(k_he4_al25__p_si28)*Y(jal25)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mg24__si28)*Y(jhe4)*Y(jmg24)*state % rho - &
      screened_rates(k_he4_si28__p_p31)*Y(jhe4)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_si28__s32)*Y(jhe4)*Y(jsi28)*state % rho + &
      screened_rates(k_p28__si28__weak__wc12)*Y(jp28) + screened_rates(k_p29__p_si28)* &
      Y(jp29) + screened_rates(k_p_al27__si28)*Y(jal27)*Y(jp)*state % rho + &
      screened_rates(k_p_p31__he4_si28)*Y(jp31)*Y(jp)*state % rho - &
      screened_rates(k_p_si28__he4_al25)*Y(jp)*Y(jsi28)*state % rho - &
      screened_rates(k_p_si28__p29)*Y(jp)*Y(jsi28)*state % rho + &
      screened_rates(k_s29__p_si28__weak__wc12)*Y(js29) + &
      screened_rates(k_s32__he4_si28)*Y(js32) - screened_rates(k_si28__he4_mg24)* &
      Y(jsi28) - screened_rates(k_si28__p_al27)*Y(jsi28) &
       )

    ydot_nuc(jsi29) = ( &
      screened_rates(k_he4_al26__p_si29)*Y(jal26)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mg25__si29)*Y(jhe4)*Y(jmg25)*state % rho - &
      screened_rates(k_he4_si29__p_p32)*Y(jhe4)*Y(jsi29)*state % rho - &
      screened_rates(k_he4_si29__s33)*Y(jhe4)*Y(jsi29)*state % rho + &
      screened_rates(k_p29__si29__weak__wc12)*Y(jp29) + screened_rates(k_p30__p_si29)* &
      Y(jp30) + screened_rates(k_p_al28__si29)*Y(jal28)*Y(jp)*state % rho + &
      screened_rates(k_p_p32__he4_si29)*Y(jp32)*Y(jp)*state % rho - &
      screened_rates(k_p_si29__he4_al26)*Y(jp)*Y(jsi29)*state % rho - &
      screened_rates(k_p_si29__p30)*Y(jp)*Y(jsi29)*state % rho + &
      screened_rates(k_s33__he4_si29)*Y(js33) - screened_rates(k_si29__he4_mg25)* &
      Y(jsi29) - screened_rates(k_si29__p_al28)*Y(jsi29) &
       )

    ydot_nuc(jsi30) = ( &
      screened_rates(k_he4_al27__p_si30)*Y(jal27)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mg26__si30)*Y(jhe4)*Y(jmg26)*state % rho - &
      screened_rates(k_he4_si30__p_p33)*Y(jhe4)*Y(jsi30)*state % rho - &
      screened_rates(k_he4_si30__s34)*Y(jhe4)*Y(jsi30)*state % rho + &
      screened_rates(k_p30__si30__weak__wc12)*Y(jp30) + screened_rates(k_p31__p_si30)* &
      Y(jp31) + screened_rates(k_p_p33__he4_si30)*Y(jp33)*Y(jp)*state % rho - &
      screened_rates(k_p_si30__he4_al27)*Y(jp)*Y(jsi30)*state % rho - &
      screened_rates(k_p_si30__p31)*Y(jp)*Y(jsi30)*state % rho + &
      screened_rates(k_s34__he4_si30)*Y(js34) - screened_rates(k_si30__he4_mg26)* &
      Y(jsi30) &
       )

    ydot_nuc(jp26) = ( &
      screened_rates(k_cl30__he4_p26)*Y(jcl30) + screened_rates(k_he4_al22__p26)*Y(jal22)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_p26__cl30)*Y(jhe4)*Y(jp26)* &
      state % rho - screened_rates(k_he4_p26__p_s29)*Y(jhe4)*Y(jp26)*state % rho - &
      screened_rates(k_p26__he4_al22)*Y(jp26) - screened_rates(k_p26__p_si25)*Y(jp26) &
      - screened_rates(k_p26__si26__weak__wc12)*Y(jp26) + &
      screened_rates(k_p_s29__he4_p26)*Y(jp)*Y(js29)*state % rho + &
      screened_rates(k_p_si25__p26)*Y(jp)*Y(jsi25)*state % rho &
       )

    ydot_nuc(jp27) = ( &
      screened_rates(k_cl31__he4_p27)*Y(jcl31) + screened_rates(k_he4_al23__p27)*Y(jal23)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_p27__cl31)*Y(jhe4)*Y(jp27)* &
      state % rho - screened_rates(k_he4_p27__p_s30)*Y(jhe4)*Y(jp27)*state % rho + &
      screened_rates(k_he4_si24__p_p27)*Y(jhe4)*Y(jsi24)*state % rho - &
      screened_rates(k_p27__he4_al23)*Y(jp27) - &
      screened_rates(k_p27__p_al26__weak__wc12)*Y(jp27) - screened_rates(k_p27__p_si26) &
      *Y(jp27) - screened_rates(k_p27__si27__weak__wc12)*Y(jp27) - &
      screened_rates(k_p_p27__he4_si24)*Y(jp27)*Y(jp)*state % rho - &
      screened_rates(k_p_p27__s28)*Y(jp27)*Y(jp)*state % rho + &
      screened_rates(k_p_s30__he4_p27)*Y(jp)*Y(js30)*state % rho + &
      screened_rates(k_p_si26__p27)*Y(jp)*Y(jsi26)*state % rho + &
      screened_rates(k_s28__p_p27)*Y(js28) &
       )

    ydot_nuc(jp28) = ( &
      screened_rates(k_cl32__he4_p28)*Y(jcl32) + screened_rates(k_he4_al24__p28)*Y(jal24)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_p28__cl32)*Y(jhe4)*Y(jp28)* &
      state % rho - screened_rates(k_he4_p28__p_s31)*Y(jhe4)*Y(jp28)*state % rho + &
      screened_rates(k_he4_si25__p_p28)*Y(jhe4)*Y(jsi25)*state % rho - &
      screened_rates(k_p28__he4_al24)*Y(jp28) - &
      screened_rates(k_p28__he4_mg24__weak__wc12)*Y(jp28) - &
      screened_rates(k_p28__p_al27__weak__wc12)*Y(jp28) - screened_rates(k_p28__p_si27) &
      *Y(jp28) - screened_rates(k_p28__si28__weak__wc12)*Y(jp28) - &
      screened_rates(k_p_p28__he4_si25)*Y(jp28)*Y(jp)*state % rho - &
      screened_rates(k_p_p28__s29)*Y(jp28)*Y(jp)*state % rho + &
      screened_rates(k_p_s31__he4_p28)*Y(jp)*Y(js31)*state % rho + &
      screened_rates(k_p_si27__p28)*Y(jp)*Y(jsi27)*state % rho + &
      screened_rates(k_s28__p28__weak__wc12)*Y(js28) + screened_rates(k_s29__p_p28)* &
      Y(js29) &
       )

    ydot_nuc(jp29) = ( &
      screened_rates(k_ar31__p_p_p29__weak__wc12)*Y(jar31) + screened_rates(k_cl33__he4_p29)* &
      Y(jcl33) + screened_rates(k_he4_al25__p29)*Y(jal25)*Y(jhe4)*state % rho &
      - screened_rates(k_he4_p29__cl33)*Y(jhe4)*Y(jp29)*state % rho - &
      screened_rates(k_he4_p29__p_s32)*Y(jhe4)*Y(jp29)*state % rho + &
      screened_rates(k_he4_si26__p_p29)*Y(jhe4)*Y(jsi26)*state % rho - &
      screened_rates(k_p29__he4_al25)*Y(jp29) - screened_rates(k_p29__p_si28)*Y(jp29) &
      - screened_rates(k_p29__si29__weak__wc12)*Y(jp29) - &
      screened_rates(k_p_p29__he4_si26)*Y(jp29)*Y(jp)*state % rho - &
      screened_rates(k_p_p29__s30)*Y(jp29)*Y(jp)*state % rho + &
      screened_rates(k_p_s32__he4_p29)*Y(jp)*Y(js32)*state % rho + &
      screened_rates(k_p_si28__p29)*Y(jp)*Y(jsi28)*state % rho + &
      screened_rates(k_s29__p29__weak__wc12)*Y(js29) + screened_rates(k_s30__p_p29)* &
      Y(js30) &
       )

    ydot_nuc(jp30) = ( &
      screened_rates(k_cl31__p_p30__weak__wc12)*Y(jcl31) + screened_rates(k_cl34__he4_p30)* &
      Y(jcl34) + screened_rates(k_he4_al26__p30)*Y(jal26)*Y(jhe4)*state % rho &
      - screened_rates(k_he4_p30__cl34)*Y(jhe4)*Y(jp30)*state % rho - &
      screened_rates(k_he4_p30__p_s33)*Y(jhe4)*Y(jp30)*state % rho + &
      screened_rates(k_he4_si27__p_p30)*Y(jhe4)*Y(jsi27)*state % rho - &
      screened_rates(k_p30__he4_al26)*Y(jp30) - screened_rates(k_p30__p_si29)*Y(jp30) &
      - screened_rates(k_p30__si30__weak__wc12)*Y(jp30) - &
      screened_rates(k_p_p30__he4_si27)*Y(jp30)*Y(jp)*state % rho - &
      screened_rates(k_p_p30__s31)*Y(jp30)*Y(jp)*state % rho + &
      screened_rates(k_p_s33__he4_p30)*Y(jp)*Y(js33)*state % rho + &
      screened_rates(k_p_si29__p30)*Y(jp)*Y(jsi29)*state % rho + &
      screened_rates(k_s30__p30__weak__wc12)*Y(js30) + screened_rates(k_s31__p_p30)* &
      Y(js31) &
       )

    ydot_nuc(jp31) = ( &
      screened_rates(k_cl32__p_p31__weak__wc12)*Y(jcl32) + screened_rates(k_cl35__he4_p31)* &
      Y(jcl35) + screened_rates(k_he4_al27__p31)*Y(jal27)*Y(jhe4)*state % rho &
      - screened_rates(k_he4_p31__cl35)*Y(jhe4)*Y(jp31)*state % rho - &
      screened_rates(k_he4_p31__p_s34)*Y(jhe4)*Y(jp31)*state % rho + &
      screened_rates(k_he4_si28__p_p31)*Y(jhe4)*Y(jsi28)*state % rho - &
      screened_rates(k_p31__he4_al27)*Y(jp31) - screened_rates(k_p31__p_si30)*Y(jp31) &
      - screened_rates(k_p_p31__he4_si28)*Y(jp31)*Y(jp)*state % rho - &
      screened_rates(k_p_p31__s32)*Y(jp31)*Y(jp)*state % rho + &
      screened_rates(k_p_s34__he4_p31)*Y(jp)*Y(js34)*state % rho + &
      screened_rates(k_p_si30__p31)*Y(jp)*Y(jsi30)*state % rho + &
      screened_rates(k_s31__p31__weak__wc12)*Y(js31) + screened_rates(k_s32__p_p31)* &
      Y(js32) &
       )

    ydot_nuc(jp32) = ( &
      screened_rates(k_cl36__he4_p32)*Y(jcl36) + screened_rates(k_he4_al28__p32)*Y(jal28)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_p32__cl36)*Y(jhe4)*Y(jp32)* &
      state % rho - screened_rates(k_he4_p32__p_s35)*Y(jhe4)*Y(jp32)*state % rho + &
      screened_rates(k_he4_si29__p_p32)*Y(jhe4)*Y(jsi29)*state % rho - &
      screened_rates(k_p32__he4_al28)*Y(jp32) - screened_rates(k_p_p32__he4_si29)* &
      Y(jp32)*Y(jp)*state % rho - screened_rates(k_p_p32__s33)*Y(jp32)* &
      Y(jp)*state % rho + screened_rates(k_p_s35__he4_p32)*Y(jp)*Y(js35)* &
      state % rho + screened_rates(k_s33__p_p32)*Y(js33) &
       )

    ydot_nuc(jp33) = ( &
      screened_rates(k_cl37__he4_p33)*Y(jcl37) - screened_rates(k_he4_p33__cl37)*Y(jhe4)* &
      Y(jp33)*state % rho + screened_rates(k_he4_si30__p_p33)*Y(jhe4)* &
      Y(jsi30)*state % rho - screened_rates(k_p_p33__he4_si30)*Y(jp33)*Y(jp)* &
      state % rho - screened_rates(k_p_p33__s34)*Y(jp33)*Y(jp)*state % rho + &
      screened_rates(k_s34__p_p33)*Y(js34) &
       )

    ydot_nuc(js28) = ( &
      screened_rates(k_ar32__he4_s28)*Y(jar32) + screened_rates(k_cl29__p_s28)*Y(jcl29) - &
      screened_rates(k_he4_s28__ar32)*Y(jhe4)*Y(js28)*state % rho - &
      screened_rates(k_he4_s28__p_cl31)*Y(jhe4)*Y(js28)*state % rho + &
      screened_rates(k_he4_si24__s28)*Y(jhe4)*Y(jsi24)*state % rho + &
      screened_rates(k_p_cl31__he4_s28)*Y(jcl31)*Y(jp)*state % rho + &
      screened_rates(k_p_p27__s28)*Y(jp27)*Y(jp)*state % rho - &
      screened_rates(k_p_s28__cl29)*Y(jp)*Y(js28)*state % rho - &
      screened_rates(k_s28__he4_si24)*Y(js28) - screened_rates(k_s28__p28__weak__wc12)* &
      Y(js28) - screened_rates(k_s28__p_p27)*Y(js28) - &
      screened_rates(k_s28__p_si27__weak__wc12)*Y(js28) &
       )

    ydot_nuc(js29) = ( &
      screened_rates(k_ar33__he4_s29)*Y(jar33) + screened_rates(k_cl29__s29__weak__bqa_pos_)* &
      Y(jcl29) + screened_rates(k_cl30__p_s29)*Y(jcl30) + &
      screened_rates(k_he4_p26__p_s29)*Y(jhe4)*Y(jp26)*state % rho - &
      screened_rates(k_he4_s29__ar33)*Y(jhe4)*Y(js29)*state % rho - &
      screened_rates(k_he4_s29__p_cl32)*Y(jhe4)*Y(js29)*state % rho + &
      screened_rates(k_he4_si25__s29)*Y(jhe4)*Y(jsi25)*state % rho + &
      screened_rates(k_p_cl32__he4_s29)*Y(jcl32)*Y(jp)*state % rho + &
      screened_rates(k_p_p28__s29)*Y(jp28)*Y(jp)*state % rho - &
      screened_rates(k_p_s29__cl30)*Y(jp)*Y(js29)*state % rho - &
      screened_rates(k_p_s29__he4_p26)*Y(jp)*Y(js29)*state % rho - &
      screened_rates(k_s29__he4_si25)*Y(js29) - screened_rates(k_s29__p29__weak__wc12)* &
      Y(js29) - screened_rates(k_s29__p_p28)*Y(js29) - &
      screened_rates(k_s29__p_si28__weak__wc12)*Y(js29) &
       )

    ydot_nuc(js30) = ( &
      screened_rates(k_ar31__p_s30__weak__wc17)*Y(jar31) + screened_rates(k_ar34__he4_s30)* &
      Y(jar34) + screened_rates(k_cl30__s30__weak__bqa_pos_)*Y(jcl30) + &
      screened_rates(k_cl31__p_s30)*Y(jcl31) + screened_rates(k_he4_p27__p_s30)* &
      Y(jhe4)*Y(jp27)*state % rho - screened_rates(k_he4_s30__ar34)*Y(jhe4)* &
      Y(js30)*state % rho - screened_rates(k_he4_s30__p_cl33)*Y(jhe4)*Y(js30) &
      *state % rho + screened_rates(k_he4_si26__s30)*Y(jhe4)*Y(jsi26)*state % rho &
      + screened_rates(k_p_cl33__he4_s30)*Y(jcl33)*Y(jp)*state % rho + &
      screened_rates(k_p_p29__s30)*Y(jp29)*Y(jp)*state % rho - &
      screened_rates(k_p_s30__cl31)*Y(jp)*Y(js30)*state % rho - &
      screened_rates(k_p_s30__he4_p27)*Y(jp)*Y(js30)*state % rho - &
      screened_rates(k_s30__he4_si26)*Y(js30) - screened_rates(k_s30__p30__weak__wc12)* &
      Y(js30) - screened_rates(k_s30__p_p29)*Y(js30) &
       )

    ydot_nuc(js31) = ( &
      screened_rates(k_ar32__p_s31__weak__wc12)*Y(jar32) + screened_rates(k_ar35__he4_s31)* &
      Y(jar35) + screened_rates(k_cl31__s31__weak__wc12)*Y(jcl31) + &
      screened_rates(k_cl32__p_s31)*Y(jcl32) + screened_rates(k_he4_p28__p_s31)* &
      Y(jhe4)*Y(jp28)*state % rho - screened_rates(k_he4_s31__ar35)*Y(jhe4)* &
      Y(js31)*state % rho - screened_rates(k_he4_s31__p_cl34)*Y(jhe4)*Y(js31) &
      *state % rho + screened_rates(k_he4_si27__s31)*Y(jhe4)*Y(jsi27)*state % rho &
      + screened_rates(k_p_cl34__he4_s31)*Y(jcl34)*Y(jp)*state % rho + &
      screened_rates(k_p_p30__s31)*Y(jp30)*Y(jp)*state % rho - &
      screened_rates(k_p_s31__cl32)*Y(jp)*Y(js31)*state % rho - &
      screened_rates(k_p_s31__he4_p28)*Y(jp)*Y(js31)*state % rho - &
      screened_rates(k_s31__he4_si27)*Y(js31) - screened_rates(k_s31__p31__weak__wc12)* &
      Y(js31) - screened_rates(k_s31__p_p30)*Y(js31) &
       )

    ydot_nuc(js32) = ( &
      screened_rates(k_ar33__p_s32__weak__wc12)*Y(jar33) + screened_rates(k_ar36__he4_s32)* &
      Y(jar36) + screened_rates(k_cl32__s32__weak__wc12)*Y(jcl32) + &
      screened_rates(k_cl33__p_s32)*Y(jcl33) + screened_rates(k_he4_p29__p_s32)* &
      Y(jhe4)*Y(jp29)*state % rho - screened_rates(k_he4_s32__ar36)*Y(jhe4)* &
      Y(js32)*state % rho - screened_rates(k_he4_s32__p_cl35)*Y(jhe4)*Y(js32) &
      *state % rho + screened_rates(k_he4_si28__s32)*Y(jhe4)*Y(jsi28)*state % rho &
      + screened_rates(k_k36__he4_s32__weak__wc12)*Y(jk36) + &
      screened_rates(k_p_cl35__he4_s32)*Y(jcl35)*Y(jp)*state % rho + &
      screened_rates(k_p_p31__s32)*Y(jp31)*Y(jp)*state % rho - &
      screened_rates(k_p_s32__cl33)*Y(jp)*Y(js32)*state % rho - &
      screened_rates(k_p_s32__he4_p29)*Y(jp)*Y(js32)*state % rho - &
      screened_rates(k_s32__he4_si28)*Y(js32) - screened_rates(k_s32__p_p31)*Y(js32) &
       )

    ydot_nuc(js33) = ( &
      screened_rates(k_ar37__he4_s33)*Y(jar37) + screened_rates(k_cl33__s33__weak__wc12)* &
      Y(jcl33) + screened_rates(k_cl34__p_s33)*Y(jcl34) + &
      screened_rates(k_he4_p30__p_s33)*Y(jhe4)*Y(jp30)*state % rho - &
      screened_rates(k_he4_s33__ar37)*Y(jhe4)*Y(js33)*state % rho - &
      screened_rates(k_he4_s33__p_cl36)*Y(jhe4)*Y(js33)*state % rho + &
      screened_rates(k_he4_si29__s33)*Y(jhe4)*Y(jsi29)*state % rho + &
      screened_rates(k_p_cl36__he4_s33)*Y(jcl36)*Y(jp)*state % rho + &
      screened_rates(k_p_p32__s33)*Y(jp32)*Y(jp)*state % rho - &
      screened_rates(k_p_s33__cl34)*Y(jp)*Y(js33)*state % rho - &
      screened_rates(k_p_s33__he4_p30)*Y(jp)*Y(js33)*state % rho - &
      screened_rates(k_s33__he4_si29)*Y(js33) - screened_rates(k_s33__p_p32)*Y(js33) &
       )

    ydot_nuc(js34) = ( &
      screened_rates(k_ar38__he4_s34)*Y(jar38) + screened_rates(k_cl34__s34__weak__wc12)* &
      Y(jcl34) + screened_rates(k_cl35__p_s34)*Y(jcl35) + &
      screened_rates(k_he4_p31__p_s34)*Y(jhe4)*Y(jp31)*state % rho - &
      screened_rates(k_he4_s34__ar38)*Y(jhe4)*Y(js34)*state % rho - &
      screened_rates(k_he4_s34__p_cl37)*Y(jhe4)*Y(js34)*state % rho + &
      screened_rates(k_he4_si30__s34)*Y(jhe4)*Y(jsi30)*state % rho + &
      screened_rates(k_p_cl37__he4_s34)*Y(jcl37)*Y(jp)*state % rho + &
      screened_rates(k_p_p33__s34)*Y(jp33)*Y(jp)*state % rho - &
      screened_rates(k_p_s34__cl35)*Y(jp)*Y(js34)*state % rho - &
      screened_rates(k_p_s34__he4_p31)*Y(jp)*Y(js34)*state % rho - &
      screened_rates(k_s34__he4_si30)*Y(js34) - screened_rates(k_s34__p_p33)*Y(js34) &
       )

    ydot_nuc(js35) = ( &
      screened_rates(k_ar39__he4_s35)*Y(jar39) + screened_rates(k_cl36__p_s35)*Y(jcl36) + &
      screened_rates(k_he4_p32__p_s35)*Y(jhe4)*Y(jp32)*state % rho - &
      screened_rates(k_he4_s35__ar39)*Y(jhe4)*Y(js35)*state % rho - &
      screened_rates(k_p_s35__cl36)*Y(jp)*Y(js35)*state % rho - &
      screened_rates(k_p_s35__he4_p32)*Y(jp)*Y(js35)*state % rho &
       )

    ydot_nuc(jcl29) = ( &
      -screened_rates(k_cl29__p_s28)*Y(jcl29) - screened_rates(k_cl29__s29__weak__bqa_pos_)* &
      Y(jcl29) - screened_rates(k_he4_cl29__k33)*Y(jcl29)*Y(jhe4)*state % rho &
      - screened_rates(k_he4_cl29__p_ar32)*Y(jcl29)*Y(jhe4)*state % rho + &
      screened_rates(k_k33__he4_cl29)*Y(jk33) + screened_rates(k_p_ar32__he4_cl29)* &
      Y(jar32)*Y(jp)*state % rho + screened_rates(k_p_s28__cl29)*Y(jp)* &
      Y(js28)*state % rho &
       )

    ydot_nuc(jcl30) = ( &
      screened_rates(k_ar31__p_cl30)*Y(jar31) - screened_rates(k_cl30__he4_p26)*Y(jcl30) - &
      screened_rates(k_cl30__p_s29)*Y(jcl30) - &
      screened_rates(k_cl30__s30__weak__bqa_pos_)*Y(jcl30) - &
      screened_rates(k_he4_cl30__k34)*Y(jcl30)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl30__p_ar33)*Y(jcl30)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_p26__cl30)*Y(jhe4)*Y(jp26)*state % rho + &
      screened_rates(k_k34__he4_cl30)*Y(jk34) + screened_rates(k_p_ar33__he4_cl30)* &
      Y(jar33)*Y(jp)*state % rho - screened_rates(k_p_cl30__ar31)*Y(jcl30)* &
      Y(jp)*state % rho + screened_rates(k_p_s29__cl30)*Y(jp)*Y(js29)* &
      state % rho &
       )

    ydot_nuc(jcl31) = ( &
      screened_rates(k_ar31__cl31__weak__wc12)*Y(jar31) + screened_rates(k_ar32__p_cl31)* &
      Y(jar32) - screened_rates(k_cl31__he4_p27)*Y(jcl31) - &
      screened_rates(k_cl31__p_p30__weak__wc12)*Y(jcl31) - &
      screened_rates(k_cl31__p_s30)*Y(jcl31) - screened_rates(k_cl31__s31__weak__wc12)* &
      Y(jcl31) - screened_rates(k_he4_cl31__k35)*Y(jcl31)*Y(jhe4)*state % rho &
      - screened_rates(k_he4_cl31__p_ar34)*Y(jcl31)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_p27__cl31)*Y(jhe4)*Y(jp27)*state % rho + &
      screened_rates(k_he4_s28__p_cl31)*Y(jhe4)*Y(js28)*state % rho + &
      screened_rates(k_k35__he4_cl31)*Y(jk35) + screened_rates(k_p_ar34__he4_cl31)* &
      Y(jar34)*Y(jp)*state % rho - screened_rates(k_p_cl31__ar32)*Y(jcl31)* &
      Y(jp)*state % rho - screened_rates(k_p_cl31__he4_s28)*Y(jcl31)*Y(jp)* &
      state % rho + screened_rates(k_p_s30__cl31)*Y(jp)*Y(js30)*state % rho &
       )

    ydot_nuc(jcl32) = ( &
      screened_rates(k_ar32__cl32__weak__wc12)*Y(jar32) + screened_rates(k_ar33__p_cl32)* &
      Y(jar33) - screened_rates(k_cl32__he4_p28)*Y(jcl32) - &
      screened_rates(k_cl32__he4_si28__weak__wc12)*Y(jcl32) - &
      screened_rates(k_cl32__p_p31__weak__wc12)*Y(jcl32) - &
      screened_rates(k_cl32__p_s31)*Y(jcl32) - screened_rates(k_cl32__s32__weak__wc12)* &
      Y(jcl32) - screened_rates(k_he4_cl32__k36)*Y(jcl32)*Y(jhe4)*state % rho &
      - screened_rates(k_he4_cl32__p_ar35)*Y(jcl32)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_p28__cl32)*Y(jhe4)*Y(jp28)*state % rho + &
      screened_rates(k_he4_s29__p_cl32)*Y(jhe4)*Y(js29)*state % rho + &
      screened_rates(k_k36__he4_cl32)*Y(jk36) + screened_rates(k_p_ar35__he4_cl32)* &
      Y(jar35)*Y(jp)*state % rho - screened_rates(k_p_cl32__ar33)*Y(jcl32)* &
      Y(jp)*state % rho - screened_rates(k_p_cl32__he4_s29)*Y(jcl32)*Y(jp)* &
      state % rho + screened_rates(k_p_s31__cl32)*Y(jp)*Y(js31)*state % rho &
       )

    ydot_nuc(jcl33) = ( &
      screened_rates(k_ar33__cl33__weak__wc12)*Y(jar33) + screened_rates(k_ar34__p_cl33)* &
      Y(jar34) + screened_rates(k_ca35__p_p_cl33__weak__wc12)*Y(jca35) - &
      screened_rates(k_cl33__he4_p29)*Y(jcl33) - screened_rates(k_cl33__p_s32)* &
      Y(jcl33) - screened_rates(k_cl33__s33__weak__wc12)*Y(jcl33) - &
      screened_rates(k_he4_cl33__k37)*Y(jcl33)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl33__p_ar36)*Y(jcl33)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_p29__cl33)*Y(jhe4)*Y(jp29)*state % rho + &
      screened_rates(k_he4_s30__p_cl33)*Y(jhe4)*Y(js30)*state % rho + &
      screened_rates(k_k37__he4_cl33)*Y(jk37) + screened_rates(k_p_ar36__he4_cl33)* &
      Y(jar36)*Y(jp)*state % rho - screened_rates(k_p_cl33__ar34)*Y(jcl33)* &
      Y(jp)*state % rho - screened_rates(k_p_cl33__he4_s30)*Y(jcl33)*Y(jp)* &
      state % rho + screened_rates(k_p_s32__cl33)*Y(jp)*Y(js32)*state % rho &
       )

    ydot_nuc(jcl34) = ( &
      screened_rates(k_ar34__cl34__weak__wc12)*Y(jar34) + screened_rates(k_ar35__p_cl34)* &
      Y(jar35) - screened_rates(k_cl34__he4_p30)*Y(jcl34) - &
      screened_rates(k_cl34__p_s33)*Y(jcl34) - screened_rates(k_cl34__s34__weak__wc12)* &
      Y(jcl34) - screened_rates(k_he4_cl34__k38)*Y(jcl34)*Y(jhe4)*state % rho &
      - screened_rates(k_he4_cl34__p_ar37)*Y(jcl34)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_p30__cl34)*Y(jhe4)*Y(jp30)*state % rho + &
      screened_rates(k_he4_s31__p_cl34)*Y(jhe4)*Y(js31)*state % rho + &
      screened_rates(k_k35__p_cl34__weak__wc12)*Y(jk35) + &
      screened_rates(k_k38__he4_cl34)*Y(jk38) + screened_rates(k_p_ar37__he4_cl34)* &
      Y(jar37)*Y(jp)*state % rho - screened_rates(k_p_cl34__ar35)*Y(jcl34)* &
      Y(jp)*state % rho - screened_rates(k_p_cl34__he4_s31)*Y(jcl34)*Y(jp)* &
      state % rho + screened_rates(k_p_s33__cl34)*Y(jp)*Y(js33)*state % rho &
       )

    ydot_nuc(jcl35) = ( &
      screened_rates(k_ar35__cl35__weak__wc12)*Y(jar35) + screened_rates(k_ar36__p_cl35)* &
      Y(jar36) - screened_rates(k_cl35__he4_p31)*Y(jcl35) - &
      screened_rates(k_cl35__p_s34)*Y(jcl35) - screened_rates(k_he4_cl35__k39)* &
      Y(jcl35)*Y(jhe4)*state % rho - screened_rates(k_he4_cl35__p_ar38)* &
      Y(jcl35)*Y(jhe4)*state % rho + screened_rates(k_he4_p31__cl35)*Y(jhe4)* &
      Y(jp31)*state % rho + screened_rates(k_he4_s32__p_cl35)*Y(jhe4)*Y(js32) &
      *state % rho + screened_rates(k_k36__p_cl35__weak__wc12)*Y(jk36) + &
      screened_rates(k_k39__he4_cl35)*Y(jk39) + screened_rates(k_p_ar38__he4_cl35)* &
      Y(jar38)*Y(jp)*state % rho - screened_rates(k_p_cl35__ar36)*Y(jcl35)* &
      Y(jp)*state % rho - screened_rates(k_p_cl35__he4_s32)*Y(jcl35)*Y(jp)* &
      state % rho + screened_rates(k_p_s34__cl35)*Y(jp)*Y(js34)*state % rho &
       )

    ydot_nuc(jcl36) = ( &
      screened_rates(k_ar37__p_cl36)*Y(jar37) - screened_rates(k_cl36__he4_p32)*Y(jcl36) - &
      screened_rates(k_cl36__p_s35)*Y(jcl36) - screened_rates(k_he4_cl36__k40)* &
      Y(jcl36)*Y(jhe4)*state % rho - screened_rates(k_he4_cl36__p_ar39)* &
      Y(jcl36)*Y(jhe4)*state % rho + screened_rates(k_he4_p32__cl36)*Y(jhe4)* &
      Y(jp32)*state % rho + screened_rates(k_he4_s33__p_cl36)*Y(jhe4)*Y(js33) &
      *state % rho + screened_rates(k_k40__he4_cl36)*Y(jk40) + &
      screened_rates(k_p_ar39__he4_cl36)*Y(jar39)*Y(jp)*state % rho - &
      screened_rates(k_p_cl36__ar37)*Y(jcl36)*Y(jp)*state % rho - &
      screened_rates(k_p_cl36__he4_s33)*Y(jcl36)*Y(jp)*state % rho + &
      screened_rates(k_p_s35__cl36)*Y(jp)*Y(js35)*state % rho &
       )

    ydot_nuc(jcl37) = ( &
      screened_rates(k_ar37__cl37__weak__wc12)*Y(jar37) + screened_rates(k_ar38__p_cl37)* &
      Y(jar38) - screened_rates(k_cl37__he4_p33)*Y(jcl37) - &
      screened_rates(k_he4_cl37__k41)*Y(jcl37)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_p33__cl37)*Y(jhe4)*Y(jp33)*state % rho + &
      screened_rates(k_he4_s34__p_cl37)*Y(jhe4)*Y(js34)*state % rho + &
      screened_rates(k_k41__he4_cl37)*Y(jk41) - screened_rates(k_p_cl37__ar38)* &
      Y(jcl37)*Y(jp)*state % rho - screened_rates(k_p_cl37__he4_s34)*Y(jcl37) &
      *Y(jp)*state % rho &
       )

    ydot_nuc(jar31) = ( &
      -screened_rates(k_ar31__cl31__weak__wc12)*Y(jar31) - screened_rates(k_ar31__p_cl30)* &
      Y(jar31) - screened_rates(k_ar31__p_p_p29__weak__wc12)*Y(jar31) - &
      screened_rates(k_ar31__p_s30__weak__wc17)*Y(jar31) + &
      screened_rates(k_ca35__he4_ar31)*Y(jca35) - screened_rates(k_he4_ar31__ca35)* &
      Y(jar31)*Y(jhe4)*state % rho - screened_rates(k_he4_ar31__p_k34)* &
      Y(jar31)*Y(jhe4)*state % rho + screened_rates(k_p_cl30__ar31)*Y(jcl30)* &
      Y(jp)*state % rho + screened_rates(k_p_k34__he4_ar31)*Y(jk34)*Y(jp)* &
      state % rho &
       )

    ydot_nuc(jar32) = ( &
      -screened_rates(k_ar32__cl32__weak__wc12)*Y(jar32) - screened_rates(k_ar32__he4_s28)* &
      Y(jar32) - screened_rates(k_ar32__p_cl31)*Y(jar32) - &
      screened_rates(k_ar32__p_s31__weak__wc12)*Y(jar32) + &
      screened_rates(k_ca36__he4_ar32)*Y(jca36) - screened_rates(k_he4_ar32__ca36)* &
      Y(jar32)*Y(jhe4)*state % rho - screened_rates(k_he4_ar32__p_k35)* &
      Y(jar32)*Y(jhe4)*state % rho + screened_rates(k_he4_cl29__p_ar32)* &
      Y(jcl29)*Y(jhe4)*state % rho + screened_rates(k_he4_s28__ar32)*Y(jhe4)* &
      Y(js28)*state % rho + screened_rates(k_k33__p_ar32)*Y(jk33) - &
      screened_rates(k_p_ar32__he4_cl29)*Y(jar32)*Y(jp)*state % rho - &
      screened_rates(k_p_ar32__k33)*Y(jar32)*Y(jp)*state % rho + &
      screened_rates(k_p_cl31__ar32)*Y(jcl31)*Y(jp)*state % rho + &
      screened_rates(k_p_k35__he4_ar32)*Y(jk35)*Y(jp)*state % rho &
       )

    ydot_nuc(jar33) = ( &
      -screened_rates(k_ar33__cl33__weak__wc12)*Y(jar33) - screened_rates(k_ar33__he4_s29)* &
      Y(jar33) - screened_rates(k_ar33__p_cl32)*Y(jar33) - &
      screened_rates(k_ar33__p_s32__weak__wc12)*Y(jar33) + &
      screened_rates(k_ca37__he4_ar33)*Y(jca37) - screened_rates(k_he4_ar33__ca37)* &
      Y(jar33)*Y(jhe4)*state % rho - screened_rates(k_he4_ar33__p_k36)* &
      Y(jar33)*Y(jhe4)*state % rho + screened_rates(k_he4_cl30__p_ar33)* &
      Y(jcl30)*Y(jhe4)*state % rho + screened_rates(k_he4_s29__ar33)*Y(jhe4)* &
      Y(js29)*state % rho + screened_rates(k_k33__ar33__weak__bqa_pos_)*Y(jk33) &
      + screened_rates(k_k34__p_ar33)*Y(jk34) - screened_rates(k_p_ar33__he4_cl30)* &
      Y(jar33)*Y(jp)*state % rho - screened_rates(k_p_ar33__k34)*Y(jar33)* &
      Y(jp)*state % rho + screened_rates(k_p_cl32__ar33)*Y(jcl32)*Y(jp)* &
      state % rho + screened_rates(k_p_k36__he4_ar33)*Y(jk36)*Y(jp)*state % rho &
       )

    ydot_nuc(jar34) = ( &
      -screened_rates(k_ar34__cl34__weak__wc12)*Y(jar34) - screened_rates(k_ar34__he4_s30)* &
      Y(jar34) - screened_rates(k_ar34__p_cl33)*Y(jar34) + &
      screened_rates(k_ca35__p_ar34__weak__wc17)*Y(jca35) + &
      screened_rates(k_ca38__he4_ar34)*Y(jca38) - screened_rates(k_he4_ar34__ca38)* &
      Y(jar34)*Y(jhe4)*state % rho - screened_rates(k_he4_ar34__p_k37)* &
      Y(jar34)*Y(jhe4)*state % rho + screened_rates(k_he4_cl31__p_ar34)* &
      Y(jcl31)*Y(jhe4)*state % rho + screened_rates(k_he4_s30__ar34)*Y(jhe4)* &
      Y(js30)*state % rho + screened_rates(k_k34__ar34__weak__bqa_pos_)*Y(jk34) &
      + screened_rates(k_k35__p_ar34)*Y(jk35) - screened_rates(k_p_ar34__he4_cl31)* &
      Y(jar34)*Y(jp)*state % rho - screened_rates(k_p_ar34__k35)*Y(jar34)* &
      Y(jp)*state % rho + screened_rates(k_p_cl33__ar34)*Y(jcl33)*Y(jp)* &
      state % rho + screened_rates(k_p_k37__he4_ar34)*Y(jk37)*Y(jp)*state % rho &
       )

    ydot_nuc(jar35) = ( &
      -screened_rates(k_ar35__cl35__weak__wc12)*Y(jar35) - screened_rates(k_ar35__he4_s31)* &
      Y(jar35) - screened_rates(k_ar35__p_cl34)*Y(jar35) + &
      screened_rates(k_ca36__p_ar35__weak__wc12)*Y(jca36) + &
      screened_rates(k_ca39__he4_ar35)*Y(jca39) - screened_rates(k_he4_ar35__ca39)* &
      Y(jar35)*Y(jhe4)*state % rho - screened_rates(k_he4_ar35__p_k38)* &
      Y(jar35)*Y(jhe4)*state % rho + screened_rates(k_he4_cl32__p_ar35)* &
      Y(jcl32)*Y(jhe4)*state % rho + screened_rates(k_he4_s31__ar35)*Y(jhe4)* &
      Y(js31)*state % rho + screened_rates(k_k35__ar35__weak__wc12)*Y(jk35) + &
      screened_rates(k_k36__p_ar35)*Y(jk36) - screened_rates(k_p_ar35__he4_cl32)* &
      Y(jar35)*Y(jp)*state % rho - screened_rates(k_p_ar35__k36)*Y(jar35)* &
      Y(jp)*state % rho + screened_rates(k_p_cl34__ar35)*Y(jcl34)*Y(jp)* &
      state % rho + screened_rates(k_p_k38__he4_ar35)*Y(jk38)*Y(jp)*state % rho &
       )

    ydot_nuc(jar36) = ( &
      -screened_rates(k_ar36__he4_s32)*Y(jar36) - screened_rates(k_ar36__p_cl35)*Y(jar36) + &
      screened_rates(k_ca37__p_ar36__weak__wc12)*Y(jca37) + &
      screened_rates(k_ca40__he4_ar36)*Y(jca40) - screened_rates(k_he4_ar36__ca40)* &
      Y(jar36)*Y(jhe4)*state % rho - screened_rates(k_he4_ar36__p_k39)* &
      Y(jar36)*Y(jhe4)*state % rho + screened_rates(k_he4_cl33__p_ar36)* &
      Y(jcl33)*Y(jhe4)*state % rho + screened_rates(k_he4_s32__ar36)*Y(jhe4)* &
      Y(js32)*state % rho + screened_rates(k_k36__ar36__weak__wc12)*Y(jk36) + &
      screened_rates(k_k37__p_ar36)*Y(jk37) - screened_rates(k_p_ar36__he4_cl33)* &
      Y(jar36)*Y(jp)*state % rho - screened_rates(k_p_ar36__k37)*Y(jar36)* &
      Y(jp)*state % rho + screened_rates(k_p_cl35__ar36)*Y(jcl35)*Y(jp)* &
      state % rho + screened_rates(k_p_k39__he4_ar36)*Y(jk39)*Y(jp)*state % rho + &
      screened_rates(k_sc40__he4_ar36__weak__wc12)*Y(jsc40) &
       )

    ydot_nuc(jar37) = ( &
      -screened_rates(k_ar37__cl37__weak__wc12)*Y(jar37) - screened_rates(k_ar37__he4_s33)* &
      Y(jar37) - screened_rates(k_ar37__p_cl36)*Y(jar37) + &
      screened_rates(k_ca41__he4_ar37)*Y(jca41) - screened_rates(k_he4_ar37__ca41)* &
      Y(jar37)*Y(jhe4)*state % rho - screened_rates(k_he4_ar37__p_k40)* &
      Y(jar37)*Y(jhe4)*state % rho + screened_rates(k_he4_cl34__p_ar37)* &
      Y(jcl34)*Y(jhe4)*state % rho + screened_rates(k_he4_s33__ar37)*Y(jhe4)* &
      Y(js33)*state % rho + screened_rates(k_k37__ar37__weak__wc12)*Y(jk37) + &
      screened_rates(k_k38__p_ar37)*Y(jk38) - screened_rates(k_p_ar37__he4_cl34)* &
      Y(jar37)*Y(jp)*state % rho - screened_rates(k_p_ar37__k38)*Y(jar37)* &
      Y(jp)*state % rho + screened_rates(k_p_cl36__ar37)*Y(jcl36)*Y(jp)* &
      state % rho + screened_rates(k_p_k40__he4_ar37)*Y(jk40)*Y(jp)*state % rho &
       )

    ydot_nuc(jar38) = ( &
      -screened_rates(k_ar38__he4_s34)*Y(jar38) - screened_rates(k_ar38__p_cl37)*Y(jar38) + &
      screened_rates(k_ca42__he4_ar38)*Y(jca42) - screened_rates(k_he4_ar38__ca42)* &
      Y(jar38)*Y(jhe4)*state % rho - screened_rates(k_he4_ar38__p_k41)* &
      Y(jar38)*Y(jhe4)*state % rho + screened_rates(k_he4_cl35__p_ar38)* &
      Y(jcl35)*Y(jhe4)*state % rho + screened_rates(k_he4_s34__ar38)*Y(jhe4)* &
      Y(js34)*state % rho + screened_rates(k_k38__ar38__weak__wc12)*Y(jk38) + &
      screened_rates(k_k39__p_ar38)*Y(jk39) - screened_rates(k_p_ar38__he4_cl35)* &
      Y(jar38)*Y(jp)*state % rho - screened_rates(k_p_ar38__k39)*Y(jar38)* &
      Y(jp)*state % rho + screened_rates(k_p_cl37__ar38)*Y(jcl37)*Y(jp)* &
      state % rho + screened_rates(k_p_k41__he4_ar38)*Y(jk41)*Y(jp)*state % rho &
       )

    ydot_nuc(jar39) = ( &
      -screened_rates(k_ar39__he4_s35)*Y(jar39) + screened_rates(k_ca43__he4_ar39)*Y(jca43) &
      - screened_rates(k_he4_ar39__ca43)*Y(jar39)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_cl36__p_ar39)*Y(jcl36)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_s35__ar39)*Y(jhe4)*Y(js35)*state % rho + &
      screened_rates(k_k40__p_ar39)*Y(jk40) - screened_rates(k_p_ar39__he4_cl36)* &
      Y(jar39)*Y(jp)*state % rho - screened_rates(k_p_ar39__k40)*Y(jar39)* &
      Y(jp)*state % rho &
       )

    ydot_nuc(jk33) = ( &
      screened_rates(k_ca34__p_k33)*Y(jca34) + screened_rates(k_he4_cl29__k33)*Y(jcl29)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_k33__p_ca36)*Y(jhe4)*Y(jk33) &
      *state % rho - screened_rates(k_he4_k33__sc37)*Y(jhe4)*Y(jk33)*state % rho - &
      screened_rates(k_k33__ar33__weak__bqa_pos_)*Y(jk33) - &
      screened_rates(k_k33__he4_cl29)*Y(jk33) - screened_rates(k_k33__p_ar32)*Y(jk33) &
      + screened_rates(k_p_ar32__k33)*Y(jar32)*Y(jp)*state % rho + &
      screened_rates(k_p_ca36__he4_k33)*Y(jca36)*Y(jp)*state % rho - &
      screened_rates(k_p_k33__ca34)*Y(jk33)*Y(jp)*state % rho + &
      screened_rates(k_sc37__he4_k33)*Y(jsc37) &
       )

    ydot_nuc(jk34) = ( &
      screened_rates(k_ca34__k34__weak__bqa_pos_)*Y(jca34) + screened_rates(k_ca35__p_k34)* &
      Y(jca35) + screened_rates(k_he4_ar31__p_k34)*Y(jar31)*Y(jhe4)* &
      state % rho + screened_rates(k_he4_cl30__k34)*Y(jcl30)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_k34__p_ca37)*Y(jhe4)*Y(jk34)*state % rho - &
      screened_rates(k_he4_k34__sc38)*Y(jhe4)*Y(jk34)*state % rho - &
      screened_rates(k_k34__ar34__weak__bqa_pos_)*Y(jk34) - &
      screened_rates(k_k34__he4_cl30)*Y(jk34) - screened_rates(k_k34__p_ar33)*Y(jk34) &
      + screened_rates(k_p_ar33__k34)*Y(jar33)*Y(jp)*state % rho + &
      screened_rates(k_p_ca37__he4_k34)*Y(jca37)*Y(jp)*state % rho - &
      screened_rates(k_p_k34__ca35)*Y(jk34)*Y(jp)*state % rho - &
      screened_rates(k_p_k34__he4_ar31)*Y(jk34)*Y(jp)*state % rho + &
      screened_rates(k_sc38__he4_k34)*Y(jsc38) &
       )

    ydot_nuc(jk35) = ( &
      screened_rates(k_ca35__k35__weak__wc17)*Y(jca35) + screened_rates(k_ca36__p_k35)* &
      Y(jca36) + screened_rates(k_he4_ar32__p_k35)*Y(jar32)*Y(jhe4)* &
      state % rho + screened_rates(k_he4_cl31__k35)*Y(jcl31)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_k35__p_ca38)*Y(jhe4)*Y(jk35)*state % rho - &
      screened_rates(k_he4_k35__sc39)*Y(jhe4)*Y(jk35)*state % rho - &
      screened_rates(k_k35__ar35__weak__wc12)*Y(jk35) - screened_rates(k_k35__he4_cl31) &
      *Y(jk35) - screened_rates(k_k35__p_ar34)*Y(jk35) - &
      screened_rates(k_k35__p_cl34__weak__wc12)*Y(jk35) + screened_rates(k_p_ar34__k35) &
      *Y(jar34)*Y(jp)*state % rho + screened_rates(k_p_ca38__he4_k35)* &
      Y(jca38)*Y(jp)*state % rho - screened_rates(k_p_k35__ca36)*Y(jk35)* &
      Y(jp)*state % rho - screened_rates(k_p_k35__he4_ar32)*Y(jk35)*Y(jp)* &
      state % rho + screened_rates(k_sc39__he4_k35)*Y(jsc39) &
       )

    ydot_nuc(jk36) = ( &
      screened_rates(k_ca36__k36__weak__wc12)*Y(jca36) + screened_rates(k_ca37__p_k36)* &
      Y(jca37) + screened_rates(k_he4_ar33__p_k36)*Y(jar33)*Y(jhe4)* &
      state % rho + screened_rates(k_he4_cl32__k36)*Y(jcl32)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_k36__p_ca39)*Y(jhe4)*Y(jk36)*state % rho - &
      screened_rates(k_he4_k36__sc40)*Y(jhe4)*Y(jk36)*state % rho - &
      screened_rates(k_k36__ar36__weak__wc12)*Y(jk36) - screened_rates(k_k36__he4_cl32) &
      *Y(jk36) - screened_rates(k_k36__he4_s32__weak__wc12)*Y(jk36) - &
      screened_rates(k_k36__p_ar35)*Y(jk36) - screened_rates(k_k36__p_cl35__weak__wc12) &
      *Y(jk36) + screened_rates(k_p_ar35__k36)*Y(jar35)*Y(jp)*state % rho + &
      screened_rates(k_p_ca39__he4_k36)*Y(jca39)*Y(jp)*state % rho - &
      screened_rates(k_p_k36__ca37)*Y(jk36)*Y(jp)*state % rho - &
      screened_rates(k_p_k36__he4_ar33)*Y(jk36)*Y(jp)*state % rho + &
      screened_rates(k_sc40__he4_k36)*Y(jsc40) &
       )

    ydot_nuc(jk37) = ( &
      screened_rates(k_ca37__k37__weak__wc12)*Y(jca37) + screened_rates(k_ca38__p_k37)* &
      Y(jca38) + screened_rates(k_he4_ar34__p_k37)*Y(jar34)*Y(jhe4)* &
      state % rho + screened_rates(k_he4_cl33__k37)*Y(jcl33)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_k37__p_ca40)*Y(jhe4)*Y(jk37)*state % rho - &
      screened_rates(k_he4_k37__sc41)*Y(jhe4)*Y(jk37)*state % rho - &
      screened_rates(k_k37__ar37__weak__wc12)*Y(jk37) - screened_rates(k_k37__he4_cl33) &
      *Y(jk37) - screened_rates(k_k37__p_ar36)*Y(jk37) + &
      screened_rates(k_p_ar36__k37)*Y(jar36)*Y(jp)*state % rho + &
      screened_rates(k_p_ca40__he4_k37)*Y(jca40)*Y(jp)*state % rho - &
      screened_rates(k_p_k37__ca38)*Y(jk37)*Y(jp)*state % rho - &
      screened_rates(k_p_k37__he4_ar34)*Y(jk37)*Y(jp)*state % rho + &
      screened_rates(k_sc41__he4_k37)*Y(jsc41) &
       )

    ydot_nuc(jk38) = ( &
      screened_rates(k_ca38__k38__weak__wc12)*Y(jca38) + screened_rates(k_ca39__p_k38)* &
      Y(jca39) + screened_rates(k_he4_ar35__p_k38)*Y(jar35)*Y(jhe4)* &
      state % rho + screened_rates(k_he4_cl34__k38)*Y(jcl34)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_k38__p_ca41)*Y(jhe4)*Y(jk38)*state % rho - &
      screened_rates(k_he4_k38__sc42)*Y(jhe4)*Y(jk38)*state % rho - &
      screened_rates(k_k38__ar38__weak__wc12)*Y(jk38) - screened_rates(k_k38__he4_cl34) &
      *Y(jk38) - screened_rates(k_k38__p_ar37)*Y(jk38) + &
      screened_rates(k_p_ar37__k38)*Y(jar37)*Y(jp)*state % rho + &
      screened_rates(k_p_ca41__he4_k38)*Y(jca41)*Y(jp)*state % rho - &
      screened_rates(k_p_k38__ca39)*Y(jk38)*Y(jp)*state % rho - &
      screened_rates(k_p_k38__he4_ar35)*Y(jk38)*Y(jp)*state % rho + &
      screened_rates(k_sc42__he4_k38)*Y(jsc42) &
       )

    ydot_nuc(jk39) = ( &
      screened_rates(k_ca39__k39__weak__wc12)*Y(jca39) + screened_rates(k_ca40__p_k39)* &
      Y(jca40) + screened_rates(k_he4_ar36__p_k39)*Y(jar36)*Y(jhe4)* &
      state % rho + screened_rates(k_he4_cl35__k39)*Y(jcl35)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_k39__p_ca42)*Y(jhe4)*Y(jk39)*state % rho - &
      screened_rates(k_he4_k39__sc43)*Y(jhe4)*Y(jk39)*state % rho - &
      screened_rates(k_k39__he4_cl35)*Y(jk39) - screened_rates(k_k39__p_ar38)*Y(jk39) &
      + screened_rates(k_p_ar38__k39)*Y(jar38)*Y(jp)*state % rho + &
      screened_rates(k_p_ca42__he4_k39)*Y(jca42)*Y(jp)*state % rho - &
      screened_rates(k_p_k39__ca40)*Y(jk39)*Y(jp)*state % rho - &
      screened_rates(k_p_k39__he4_ar36)*Y(jk39)*Y(jp)*state % rho + &
      screened_rates(k_sc40__p_k39__weak__wc12)*Y(jsc40) + &
      screened_rates(k_sc43__he4_k39)*Y(jsc43) &
       )

    ydot_nuc(jk40) = ( &
      screened_rates(k_ca41__p_k40)*Y(jca41) + screened_rates(k_he4_ar37__p_k40)*Y(jar37)* &
      Y(jhe4)*state % rho + screened_rates(k_he4_cl36__k40)*Y(jcl36)*Y(jhe4)* &
      state % rho - screened_rates(k_he4_k40__p_ca43)*Y(jhe4)*Y(jk40)*state % rho &
      - screened_rates(k_he4_k40__sc44)*Y(jhe4)*Y(jk40)*state % rho - &
      screened_rates(k_k40__he4_cl36)*Y(jk40) - screened_rates(k_k40__p_ar39)*Y(jk40) &
      + screened_rates(k_p_ar39__k40)*Y(jar39)*Y(jp)*state % rho + &
      screened_rates(k_p_ca43__he4_k40)*Y(jca43)*Y(jp)*state % rho - &
      screened_rates(k_p_k40__ca41)*Y(jk40)*Y(jp)*state % rho - &
      screened_rates(k_p_k40__he4_ar37)*Y(jk40)*Y(jp)*state % rho + &
      screened_rates(k_sc44__he4_k40)*Y(jsc44) &
       )

    ydot_nuc(jk41) = ( &
      screened_rates(k_ca41__k41__weak__wc12)*Y(jca41) + screened_rates(k_ca42__p_k41)* &
      Y(jca42) + screened_rates(k_he4_ar38__p_k41)*Y(jar38)*Y(jhe4)* &
      state % rho + screened_rates(k_he4_cl37__k41)*Y(jcl37)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_k41__p_ca44)*Y(jhe4)*Y(jk41)*state % rho - &
      screened_rates(k_he4_k41__sc45)*Y(jhe4)*Y(jk41)*state % rho - &
      screened_rates(k_k41__he4_cl37)*Y(jk41) + screened_rates(k_p_ca44__he4_k41)* &
      Y(jca44)*Y(jp)*state % rho - screened_rates(k_p_k41__ca42)*Y(jk41)* &
      Y(jp)*state % rho - screened_rates(k_p_k41__he4_ar38)*Y(jk41)*Y(jp)* &
      state % rho + screened_rates(k_sc45__he4_k41)*Y(jsc45) &
       )

    ydot_nuc(jca34) = ( &
      -screened_rates(k_ca34__k34__weak__bqa_pos_)*Y(jca34) - screened_rates(k_ca34__p_k33)* &
      Y(jca34) - screened_rates(k_he4_ca34__p_sc37)*Y(jca34)*Y(jhe4)* &
      state % rho - screened_rates(k_he4_ca34__ti38)*Y(jca34)*Y(jhe4)*state % rho &
      + screened_rates(k_p_k33__ca34)*Y(jk33)*Y(jp)*state % rho + &
      screened_rates(k_p_sc37__he4_ca34)*Y(jp)*Y(jsc37)*state % rho + &
      screened_rates(k_ti38__he4_ca34)*Y(jti38) &
       )

    ydot_nuc(jca35) = ( &
      -screened_rates(k_ca35__he4_ar31)*Y(jca35) - screened_rates(k_ca35__k35__weak__wc17)* &
      Y(jca35) - screened_rates(k_ca35__p_ar34__weak__wc17)*Y(jca35) - &
      screened_rates(k_ca35__p_k34)*Y(jca35) - &
      screened_rates(k_ca35__p_p_cl33__weak__wc12)*Y(jca35) + &
      screened_rates(k_he4_ar31__ca35)*Y(jar31)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca35__p_sc38)*Y(jca35)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca35__ti39)*Y(jca35)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ca35__sc36)*Y(jca35)*Y(jp)*state % rho + &
      screened_rates(k_p_k34__ca35)*Y(jk34)*Y(jp)*state % rho + &
      screened_rates(k_p_sc38__he4_ca35)*Y(jp)*Y(jsc38)*state % rho + &
      screened_rates(k_sc36__p_ca35)*Y(jsc36) + screened_rates(k_ti39__he4_ca35)* &
      Y(jti39) &
       )

    ydot_nuc(jca36) = ( &
      -screened_rates(k_ca36__he4_ar32)*Y(jca36) - screened_rates(k_ca36__k36__weak__wc12)* &
      Y(jca36) - screened_rates(k_ca36__p_ar35__weak__wc12)*Y(jca36) - &
      screened_rates(k_ca36__p_k35)*Y(jca36) + screened_rates(k_he4_ar32__ca36)* &
      Y(jar32)*Y(jhe4)*state % rho - screened_rates(k_he4_ca36__p_p_ca38)* &
      Y(jca36)*Y(jhe4)*state % rho - screened_rates(k_he4_ca36__p_sc39)* &
      Y(jca36)*Y(jhe4)*state % rho - screened_rates(k_he4_ca36__ti40)* &
      Y(jca36)*Y(jhe4)*state % rho + screened_rates(k_he4_k33__p_ca36)* &
      Y(jhe4)*Y(jk33)*state % rho - screened_rates(k_p_ca36__he4_k33)* &
      Y(jca36)*Y(jp)*state % rho - screened_rates(k_p_ca36__sc37)*Y(jca36)* &
      Y(jp)*state % rho + screened_rates(k_p_k35__ca36)*Y(jk35)*Y(jp)* &
      state % rho + 0.5d0*screened_rates(k_p_p_ca38__he4_ca36)*Y(jca38)*Y(jp)** &
      2*state % rho**2 + screened_rates(k_p_sc39__he4_ca36)*Y(jp)*Y(jsc39)* &
      state % rho + screened_rates(k_sc36__ca36__weak__bqa_pos_)*Y(jsc36) + &
      screened_rates(k_sc37__p_ca36)*Y(jsc37) + screened_rates(k_ti40__he4_ca36)* &
      Y(jti40) &
       )

    ydot_nuc(jca37) = ( &
      -screened_rates(k_ca37__he4_ar33)*Y(jca37) - screened_rates(k_ca37__k37__weak__wc12)* &
      Y(jca37) - screened_rates(k_ca37__p_ar36__weak__wc12)*Y(jca37) - &
      screened_rates(k_ca37__p_k36)*Y(jca37) + screened_rates(k_he4_ar33__ca37)* &
      Y(jar33)*Y(jhe4)*state % rho - screened_rates(k_he4_ca37__p_sc40)* &
      Y(jca37)*Y(jhe4)*state % rho - screened_rates(k_he4_ca37__ti41)* &
      Y(jca37)*Y(jhe4)*state % rho + screened_rates(k_he4_k34__p_ca37)* &
      Y(jhe4)*Y(jk34)*state % rho - screened_rates(k_p_ca37__he4_k34)* &
      Y(jca37)*Y(jp)*state % rho - screened_rates(k_p_ca37__sc38)*Y(jca37)* &
      Y(jp)*state % rho + screened_rates(k_p_k36__ca37)*Y(jk36)*Y(jp)* &
      state % rho + screened_rates(k_p_sc40__he4_ca37)*Y(jp)*Y(jsc40)*state % rho &
      + screened_rates(k_sc37__ca37__weak__bqa_pos_)*Y(jsc37) + &
      screened_rates(k_sc38__p_ca37)*Y(jsc38) + screened_rates(k_ti41__he4_ca37)* &
      Y(jti41) &
       )

    ydot_nuc(jca38) = ( &
      -screened_rates(k_ca38__he4_ar34)*Y(jca38) - screened_rates(k_ca38__k38__weak__wc12)* &
      Y(jca38) - screened_rates(k_ca38__p_k37)*Y(jca38) + &
      screened_rates(k_he4_ar34__ca38)*Y(jar34)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ca36__p_p_ca38)*Y(jca36)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca38__p_sc41)*Y(jca38)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca38__ti42)*Y(jca38)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_k35__p_ca38)*Y(jhe4)*Y(jk35)*state % rho - &
      screened_rates(k_p_ca38__he4_k35)*Y(jca38)*Y(jp)*state % rho - &
      screened_rates(k_p_ca38__sc39)*Y(jca38)*Y(jp)*state % rho + &
      screened_rates(k_p_k37__ca38)*Y(jk37)*Y(jp)*state % rho - 0.5d0* &
      screened_rates(k_p_p_ca38__he4_ca36)*Y(jca38)*Y(jp)**2*state % rho**2 + &
      screened_rates(k_p_sc41__he4_ca38)*Y(jp)*Y(jsc41)*state % rho + &
      screened_rates(k_sc38__ca38__weak__mo97)*Y(jsc38) + &
      screened_rates(k_sc39__p_ca38)*Y(jsc39) + &
      screened_rates(k_ti39__p_ca38__weak__wc12)*Y(jti39) + &
      screened_rates(k_ti42__he4_ca38)*Y(jti42) &
       )

    ydot_nuc(jca39) = ( &
      -screened_rates(k_ca39__he4_ar35)*Y(jca39) - screened_rates(k_ca39__k39__weak__wc12)* &
      Y(jca39) - screened_rates(k_ca39__p_k38)*Y(jca39) + &
      screened_rates(k_he4_ar35__ca39)*Y(jar35)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca39__p_sc42)*Y(jca39)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca39__ti43)*Y(jca39)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_k36__p_ca39)*Y(jhe4)*Y(jk36)*state % rho - &
      screened_rates(k_p_ca39__he4_k36)*Y(jca39)*Y(jp)*state % rho - &
      screened_rates(k_p_ca39__sc40)*Y(jca39)*Y(jp)*state % rho + &
      screened_rates(k_p_k38__ca39)*Y(jk38)*Y(jp)*state % rho + &
      screened_rates(k_p_sc42__he4_ca39)*Y(jp)*Y(jsc42)*state % rho + &
      screened_rates(k_sc39__ca39__weak__mo97)*Y(jsc39) + &
      screened_rates(k_sc40__p_ca39)*Y(jsc40) + &
      screened_rates(k_ti40__p_ca39__weak__wc12)*Y(jti40) + &
      screened_rates(k_ti43__he4_ca39)*Y(jti43) &
       )

    ydot_nuc(jca40) = ( &
      -screened_rates(k_ca40__he4_ar36)*Y(jca40) - screened_rates(k_ca40__p_k39)*Y(jca40) + &
      screened_rates(k_cr43__p_p_p_ca40__weak__wc12)*Y(jcr43) + &
      screened_rates(k_he4_ar36__ca40)*Y(jar36)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca40__p_sc43)*Y(jca40)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca40__ti44)*Y(jca40)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_k37__p_ca40)*Y(jhe4)*Y(jk37)*state % rho - &
      screened_rates(k_p_ca40__he4_k37)*Y(jca40)*Y(jp)*state % rho - &
      screened_rates(k_p_ca40__sc41)*Y(jca40)*Y(jp)*state % rho + &
      screened_rates(k_p_k39__ca40)*Y(jk39)*Y(jp)*state % rho + &
      screened_rates(k_p_sc43__he4_ca40)*Y(jp)*Y(jsc43)*state % rho + &
      screened_rates(k_sc40__ca40__weak__wc12)*Y(jsc40) + &
      screened_rates(k_sc41__p_ca40)*Y(jsc41) + &
      screened_rates(k_ti41__p_ca40__weak__wc12)*Y(jti41) + &
      screened_rates(k_ti44__he4_ca40)*Y(jti44) &
       )

    ydot_nuc(jca41) = ( &
      -screened_rates(k_ca41__he4_ar37)*Y(jca41) - screened_rates(k_ca41__k41__weak__wc12)* &
      Y(jca41) - screened_rates(k_ca41__p_k40)*Y(jca41) + &
      screened_rates(k_he4_ar37__ca41)*Y(jar37)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca41__p_sc44)*Y(jca41)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca41__ti45)*Y(jca41)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_k38__p_ca41)*Y(jhe4)*Y(jk38)*state % rho - &
      screened_rates(k_p_ca41__he4_k38)*Y(jca41)*Y(jp)*state % rho - &
      screened_rates(k_p_ca41__sc42)*Y(jca41)*Y(jp)*state % rho + &
      screened_rates(k_p_k40__ca41)*Y(jk40)*Y(jp)*state % rho + &
      screened_rates(k_p_sc44__he4_ca41)*Y(jp)*Y(jsc44)*state % rho + &
      screened_rates(k_sc41__ca41__weak__wc12)*Y(jsc41) + &
      screened_rates(k_sc42__p_ca41)*Y(jsc42) + screened_rates(k_ti45__he4_ca41)* &
      Y(jti45) &
       )

    ydot_nuc(jca42) = ( &
      -screened_rates(k_ca42__he4_ar38)*Y(jca42) - screened_rates(k_ca42__p_k41)*Y(jca42) + &
      screened_rates(k_he4_ar38__ca42)*Y(jar38)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca42__p_sc45)*Y(jca42)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca42__ti46)*Y(jca42)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_k39__p_ca42)*Y(jhe4)*Y(jk39)*state % rho - &
      screened_rates(k_p_ca42__he4_k39)*Y(jca42)*Y(jp)*state % rho - &
      screened_rates(k_p_ca42__sc43)*Y(jca42)*Y(jp)*state % rho + &
      screened_rates(k_p_k41__ca42)*Y(jk41)*Y(jp)*state % rho + &
      screened_rates(k_p_sc45__he4_ca42)*Y(jp)*Y(jsc45)*state % rho + &
      screened_rates(k_sc42__ca42__weak__wc12)*Y(jsc42) + &
      screened_rates(k_sc43__p_ca42)*Y(jsc43) + screened_rates(k_ti46__he4_ca42)* &
      Y(jti46) &
       )

    ydot_nuc(jca43) = ( &
      -screened_rates(k_ca43__he4_ar39)*Y(jca43) + screened_rates(k_he4_ar39__ca43)*Y(jar39) &
      *Y(jhe4)*state % rho - screened_rates(k_he4_ca43__p_sc46)*Y(jca43)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_ca43__ti47)*Y(jca43)*Y(jhe4) &
      *state % rho + screened_rates(k_he4_k40__p_ca43)*Y(jhe4)*Y(jk40)*state % rho &
      - screened_rates(k_p_ca43__he4_k40)*Y(jca43)*Y(jp)*state % rho - &
      screened_rates(k_p_ca43__sc44)*Y(jca43)*Y(jp)*state % rho + &
      screened_rates(k_p_sc46__he4_ca43)*Y(jp)*Y(jsc46)*state % rho + &
      screened_rates(k_sc43__ca43__weak__wc12)*Y(jsc43) + &
      screened_rates(k_sc44__p_ca43)*Y(jsc44) + screened_rates(k_ti47__he4_ca43)* &
      Y(jti47) &
       )

    ydot_nuc(jca44) = ( &
      -screened_rates(k_he4_ca44__ti48)*Y(jca44)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_k41__p_ca44)*Y(jhe4)*Y(jk41)*state % rho - &
      screened_rates(k_p_ca44__he4_k41)*Y(jca44)*Y(jp)*state % rho - &
      screened_rates(k_p_ca44__sc45)*Y(jca44)*Y(jp)*state % rho + &
      screened_rates(k_sc44__ca44__weak__wc12)*Y(jsc44) + &
      screened_rates(k_sc45__p_ca44)*Y(jsc45) + screened_rates(k_ti48__he4_ca44)* &
      Y(jti48) &
       )

    ydot_nuc(jsc36) = ( &
      -screened_rates(k_he4_sc36__p_ti39)*Y(jhe4)*Y(jsc36)*state % rho - &
      screened_rates(k_he4_sc36__v40)*Y(jhe4)*Y(jsc36)*state % rho + &
      screened_rates(k_p_ca35__sc36)*Y(jca35)*Y(jp)*state % rho + &
      screened_rates(k_p_ti39__he4_sc36)*Y(jp)*Y(jti39)*state % rho - &
      screened_rates(k_sc36__ca36__weak__bqa_pos_)*Y(jsc36) - &
      screened_rates(k_sc36__p_ca35)*Y(jsc36) + screened_rates(k_v40__he4_sc36)* &
      Y(jv40) &
       )

    ydot_nuc(jsc37) = ( &
      screened_rates(k_he4_ca34__p_sc37)*Y(jca34)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_k33__sc37)*Y(jhe4)*Y(jk33)*state % rho - &
      screened_rates(k_he4_sc37__p_ti40)*Y(jhe4)*Y(jsc37)*state % rho - &
      screened_rates(k_he4_sc37__v41)*Y(jhe4)*Y(jsc37)*state % rho + &
      screened_rates(k_p_ca36__sc37)*Y(jca36)*Y(jp)*state % rho - &
      screened_rates(k_p_sc37__he4_ca34)*Y(jp)*Y(jsc37)*state % rho - &
      screened_rates(k_p_sc37__ti38)*Y(jp)*Y(jsc37)*state % rho + &
      screened_rates(k_p_ti40__he4_sc37)*Y(jp)*Y(jti40)*state % rho - &
      screened_rates(k_sc37__ca37__weak__bqa_pos_)*Y(jsc37) - &
      screened_rates(k_sc37__he4_k33)*Y(jsc37) - screened_rates(k_sc37__p_ca36)* &
      Y(jsc37) + screened_rates(k_ti38__p_sc37)*Y(jti38) + &
      screened_rates(k_v41__he4_sc37)*Y(jv41) &
       )

    ydot_nuc(jsc38) = ( &
      screened_rates(k_he4_ca35__p_sc38)*Y(jca35)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_k34__sc38)*Y(jhe4)*Y(jk34)*state % rho - &
      screened_rates(k_he4_sc38__p_ti41)*Y(jhe4)*Y(jsc38)*state % rho - &
      screened_rates(k_he4_sc38__v42)*Y(jhe4)*Y(jsc38)*state % rho + &
      screened_rates(k_p_ca37__sc38)*Y(jca37)*Y(jp)*state % rho - &
      screened_rates(k_p_sc38__he4_ca35)*Y(jp)*Y(jsc38)*state % rho - &
      screened_rates(k_p_sc38__ti39)*Y(jp)*Y(jsc38)*state % rho + &
      screened_rates(k_p_ti41__he4_sc38)*Y(jp)*Y(jti41)*state % rho - &
      screened_rates(k_sc38__ca38__weak__mo97)*Y(jsc38) - &
      screened_rates(k_sc38__he4_k34)*Y(jsc38) - screened_rates(k_sc38__p_ca37)* &
      Y(jsc38) + screened_rates(k_ti38__sc38__weak__mo97)*Y(jti38) + &
      screened_rates(k_ti39__p_sc38)*Y(jti39) + screened_rates(k_v42__he4_sc38)* &
      Y(jv42) &
       )

    ydot_nuc(jsc39) = ( &
      screened_rates(k_he4_ca36__p_sc39)*Y(jca36)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_k35__sc39)*Y(jhe4)*Y(jk35)*state % rho - &
      screened_rates(k_he4_sc39__p_ti42)*Y(jhe4)*Y(jsc39)*state % rho - &
      screened_rates(k_he4_sc39__v43)*Y(jhe4)*Y(jsc39)*state % rho + &
      screened_rates(k_p_ca38__sc39)*Y(jca38)*Y(jp)*state % rho - &
      screened_rates(k_p_sc39__he4_ca36)*Y(jp)*Y(jsc39)*state % rho - &
      screened_rates(k_p_sc39__ti40)*Y(jp)*Y(jsc39)*state % rho + &
      screened_rates(k_p_ti42__he4_sc39)*Y(jp)*Y(jti42)*state % rho - &
      screened_rates(k_sc39__ca39__weak__mo97)*Y(jsc39) - &
      screened_rates(k_sc39__he4_k35)*Y(jsc39) - screened_rates(k_sc39__p_ca38)* &
      Y(jsc39) + screened_rates(k_ti39__sc39__weak__wc17)*Y(jti39) + &
      screened_rates(k_ti40__p_sc39)*Y(jti40) + screened_rates(k_v43__he4_sc39)* &
      Y(jv43) &
       )

    ydot_nuc(jsc40) = ( &
      screened_rates(k_he4_ca37__p_sc40)*Y(jca37)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_k36__sc40)*Y(jhe4)*Y(jk36)*state % rho - &
      screened_rates(k_he4_sc40__p_ti43)*Y(jhe4)*Y(jsc40)*state % rho - &
      screened_rates(k_he4_sc40__v44)*Y(jhe4)*Y(jsc40)*state % rho + &
      screened_rates(k_p_ca39__sc40)*Y(jca39)*Y(jp)*state % rho - &
      screened_rates(k_p_sc40__he4_ca37)*Y(jp)*Y(jsc40)*state % rho - &
      screened_rates(k_p_sc40__ti41)*Y(jp)*Y(jsc40)*state % rho + &
      screened_rates(k_p_ti43__he4_sc40)*Y(jp)*Y(jti43)*state % rho - &
      screened_rates(k_sc40__ca40__weak__wc12)*Y(jsc40) - &
      screened_rates(k_sc40__he4_ar36__weak__wc12)*Y(jsc40) - &
      screened_rates(k_sc40__he4_k36)*Y(jsc40) - screened_rates(k_sc40__p_ca39)* &
      Y(jsc40) - screened_rates(k_sc40__p_k39__weak__wc12)*Y(jsc40) + &
      screened_rates(k_ti40__sc40__weak__wc17)*Y(jti40) + &
      screened_rates(k_ti41__p_sc40)*Y(jti41) + screened_rates(k_v44__he4_sc40)* &
      Y(jv44) &
       )

    ydot_nuc(jsc41) = ( &
      screened_rates(k_cr43__p_p_sc41__weak__wc12)*Y(jcr43) + &
      screened_rates(k_he4_ca38__p_sc41)*Y(jca38)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_k37__sc41)*Y(jhe4)*Y(jk37)*state % rho - &
      screened_rates(k_he4_sc41__p_ti44)*Y(jhe4)*Y(jsc41)*state % rho - &
      screened_rates(k_he4_sc41__v45)*Y(jhe4)*Y(jsc41)*state % rho + &
      screened_rates(k_p_ca40__sc41)*Y(jca40)*Y(jp)*state % rho - &
      screened_rates(k_p_sc41__he4_ca38)*Y(jp)*Y(jsc41)*state % rho - &
      screened_rates(k_p_sc41__ti42)*Y(jp)*Y(jsc41)*state % rho + &
      screened_rates(k_p_ti44__he4_sc41)*Y(jp)*Y(jti44)*state % rho - &
      screened_rates(k_sc41__ca41__weak__wc12)*Y(jsc41) - &
      screened_rates(k_sc41__he4_k37)*Y(jsc41) - screened_rates(k_sc41__p_ca40)* &
      Y(jsc41) + screened_rates(k_ti41__sc41__weak__wc17)*Y(jti41) + &
      screened_rates(k_ti42__p_sc41)*Y(jti42) + screened_rates(k_v45__he4_sc41)* &
      Y(jv45) &
       )

    ydot_nuc(jsc42) = ( &
      screened_rates(k_he4_ca39__p_sc42)*Y(jca39)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_k38__sc42)*Y(jhe4)*Y(jk38)*state % rho - &
      screened_rates(k_he4_sc42__p_ti45)*Y(jhe4)*Y(jsc42)*state % rho - &
      screened_rates(k_he4_sc42__v46)*Y(jhe4)*Y(jsc42)*state % rho + &
      screened_rates(k_p_ca41__sc42)*Y(jca41)*Y(jp)*state % rho - &
      screened_rates(k_p_sc42__he4_ca39)*Y(jp)*Y(jsc42)*state % rho - &
      screened_rates(k_p_sc42__ti43)*Y(jp)*Y(jsc42)*state % rho + &
      screened_rates(k_p_ti45__he4_sc42)*Y(jp)*Y(jti45)*state % rho - &
      screened_rates(k_sc42__ca42__weak__wc12)*Y(jsc42) - &
      screened_rates(k_sc42__he4_k38)*Y(jsc42) - screened_rates(k_sc42__p_ca41)* &
      Y(jsc42) + screened_rates(k_ti42__sc42__weak__wc12)*Y(jti42) + &
      screened_rates(k_ti43__p_sc42)*Y(jti43) + screened_rates(k_v46__he4_sc42)* &
      Y(jv46) &
       )

    ydot_nuc(jsc43) = ( &
      screened_rates(k_he4_ca40__p_sc43)*Y(jca40)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_k39__sc43)*Y(jhe4)*Y(jk39)*state % rho - &
      screened_rates(k_he4_sc43__p_ti46)*Y(jhe4)*Y(jsc43)*state % rho - &
      screened_rates(k_he4_sc43__v47)*Y(jhe4)*Y(jsc43)*state % rho + &
      screened_rates(k_p_ca42__sc43)*Y(jca42)*Y(jp)*state % rho - &
      screened_rates(k_p_sc43__he4_ca40)*Y(jp)*Y(jsc43)*state % rho - &
      screened_rates(k_p_sc43__ti44)*Y(jp)*Y(jsc43)*state % rho + &
      screened_rates(k_p_ti46__he4_sc43)*Y(jp)*Y(jti46)*state % rho - &
      screened_rates(k_sc43__ca43__weak__wc12)*Y(jsc43) - &
      screened_rates(k_sc43__he4_k39)*Y(jsc43) - screened_rates(k_sc43__p_ca42)* &
      Y(jsc43) + screened_rates(k_ti43__sc43__weak__wc12)*Y(jti43) + &
      screened_rates(k_ti44__p_sc43)*Y(jti44) + screened_rates(k_v47__he4_sc43)* &
      Y(jv47) &
       )

    ydot_nuc(jsc44) = ( &
      screened_rates(k_he4_ca41__p_sc44)*Y(jca41)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_k40__sc44)*Y(jhe4)*Y(jk40)*state % rho - &
      screened_rates(k_he4_sc44__p_ti47)*Y(jhe4)*Y(jsc44)*state % rho - &
      screened_rates(k_he4_sc44__v48)*Y(jhe4)*Y(jsc44)*state % rho + &
      screened_rates(k_p_ca43__sc44)*Y(jca43)*Y(jp)*state % rho - &
      screened_rates(k_p_sc44__he4_ca41)*Y(jp)*Y(jsc44)*state % rho - &
      screened_rates(k_p_sc44__ti45)*Y(jp)*Y(jsc44)*state % rho + &
      screened_rates(k_p_ti47__he4_sc44)*Y(jp)*Y(jti47)*state % rho - &
      screened_rates(k_sc44__ca44__weak__wc12)*Y(jsc44) - &
      screened_rates(k_sc44__he4_k40)*Y(jsc44) - screened_rates(k_sc44__p_ca43)* &
      Y(jsc44) + screened_rates(k_ti44__sc44__weak__wc12)*Y(jti44) + &
      screened_rates(k_ti45__p_sc44)*Y(jti45) + screened_rates(k_v48__he4_sc44)* &
      Y(jv48) &
       )

    ydot_nuc(jsc45) = ( &
      screened_rates(k_he4_ca42__p_sc45)*Y(jca42)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_k41__sc45)*Y(jhe4)*Y(jk41)*state % rho - &
      screened_rates(k_he4_sc45__p_ti48)*Y(jhe4)*Y(jsc45)*state % rho - &
      screened_rates(k_he4_sc45__v49)*Y(jhe4)*Y(jsc45)*state % rho + &
      screened_rates(k_p_ca44__sc45)*Y(jca44)*Y(jp)*state % rho - &
      screened_rates(k_p_sc45__he4_ca42)*Y(jp)*Y(jsc45)*state % rho - &
      screened_rates(k_p_sc45__ti46)*Y(jp)*Y(jsc45)*state % rho + &
      screened_rates(k_p_ti48__he4_sc45)*Y(jp)*Y(jti48)*state % rho - &
      screened_rates(k_sc45__he4_k41)*Y(jsc45) - screened_rates(k_sc45__p_ca44)* &
      Y(jsc45) + screened_rates(k_ti45__sc45__weak__wc12)*Y(jti45) + &
      screened_rates(k_ti46__p_sc45)*Y(jti46) + screened_rates(k_v49__he4_sc45)* &
      Y(jv49) &
       )

    ydot_nuc(jsc46) = ( &
      screened_rates(k_he4_ca43__p_sc46)*Y(jca43)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_sc46__v50)*Y(jhe4)*Y(jsc46)*state % rho - &
      screened_rates(k_p_sc46__he4_ca43)*Y(jp)*Y(jsc46)*state % rho - &
      screened_rates(k_p_sc46__ti47)*Y(jp)*Y(jsc46)*state % rho + &
      screened_rates(k_ti47__p_sc46)*Y(jti47) + screened_rates(k_v50__he4_sc46)* &
      Y(jv50) &
       )

    ydot_nuc(jti38) = ( &
      screened_rates(k_cr42__he4_ti38)*Y(jcr42) + screened_rates(k_he4_ca34__ti38)*Y(jca34)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_ti38__cr42)*Y(jhe4)*Y(jti38) &
      *state % rho - screened_rates(k_he4_ti38__p_v41)*Y(jhe4)*Y(jti38)* &
      state % rho + screened_rates(k_p_sc37__ti38)*Y(jp)*Y(jsc37)*state % rho + &
      screened_rates(k_p_v41__he4_ti38)*Y(jp)*Y(jv41)*state % rho - &
      screened_rates(k_ti38__he4_ca34)*Y(jti38) - screened_rates(k_ti38__p_sc37)* &
      Y(jti38) - screened_rates(k_ti38__sc38__weak__mo97)*Y(jti38) &
       )

    ydot_nuc(jti39) = ( &
      screened_rates(k_cr43__he4_ti39)*Y(jcr43) + screened_rates(k_he4_ca35__ti39)*Y(jca35)* &
      Y(jhe4)*state % rho + screened_rates(k_he4_sc36__p_ti39)*Y(jhe4)* &
      Y(jsc36)*state % rho - screened_rates(k_he4_ti39__cr43)*Y(jhe4)* &
      Y(jti39)*state % rho - screened_rates(k_he4_ti39__p_v42)*Y(jhe4)* &
      Y(jti39)*state % rho + screened_rates(k_p_sc38__ti39)*Y(jp)*Y(jsc38)* &
      state % rho - screened_rates(k_p_ti39__he4_sc36)*Y(jp)*Y(jti39)*state % rho &
      - screened_rates(k_p_ti39__v40)*Y(jp)*Y(jti39)*state % rho + &
      screened_rates(k_p_v42__he4_ti39)*Y(jp)*Y(jv42)*state % rho - &
      screened_rates(k_ti39__he4_ca35)*Y(jti39) - &
      screened_rates(k_ti39__p_ca38__weak__wc12)*Y(jti39) - &
      screened_rates(k_ti39__p_sc38)*Y(jti39) - &
      screened_rates(k_ti39__sc39__weak__wc17)*Y(jti39) + screened_rates(k_v40__p_ti39) &
      *Y(jv40) &
       )

    ydot_nuc(jti40) = ( &
      screened_rates(k_cr44__he4_ti40)*Y(jcr44) + screened_rates(k_he4_ca36__ti40)*Y(jca36)* &
      Y(jhe4)*state % rho + screened_rates(k_he4_sc37__p_ti40)*Y(jhe4)* &
      Y(jsc37)*state % rho - screened_rates(k_he4_ti40__cr44)*Y(jhe4)* &
      Y(jti40)*state % rho - screened_rates(k_he4_ti40__p_v43)*Y(jhe4)* &
      Y(jti40)*state % rho + screened_rates(k_p_sc39__ti40)*Y(jp)*Y(jsc39)* &
      state % rho - screened_rates(k_p_ti40__he4_sc37)*Y(jp)*Y(jti40)*state % rho &
      - screened_rates(k_p_ti40__v41)*Y(jp)*Y(jti40)*state % rho + &
      screened_rates(k_p_v43__he4_ti40)*Y(jp)*Y(jv43)*state % rho - &
      screened_rates(k_ti40__he4_ca36)*Y(jti40) - &
      screened_rates(k_ti40__p_ca39__weak__wc12)*Y(jti40) - &
      screened_rates(k_ti40__p_sc39)*Y(jti40) - &
      screened_rates(k_ti40__sc40__weak__wc17)*Y(jti40) + &
      screened_rates(k_v40__ti40__weak__bqa_pos_)*Y(jv40) + &
      screened_rates(k_v41__p_ti40)*Y(jv41) &
       )

    ydot_nuc(jti41) = ( &
      screened_rates(k_cr42__p_ti41__weak__wc12)*Y(jcr42) + screened_rates(k_cr45__he4_ti41)* &
      Y(jcr45) + screened_rates(k_he4_ca37__ti41)*Y(jca37)*Y(jhe4)* &
      state % rho + screened_rates(k_he4_sc38__p_ti41)*Y(jhe4)*Y(jsc38)* &
      state % rho - screened_rates(k_he4_ti41__cr45)*Y(jhe4)*Y(jti41)*state % rho &
      - screened_rates(k_he4_ti41__p_v44)*Y(jhe4)*Y(jti41)*state % rho + &
      screened_rates(k_p_sc40__ti41)*Y(jp)*Y(jsc40)*state % rho - &
      screened_rates(k_p_ti41__he4_sc38)*Y(jp)*Y(jti41)*state % rho - &
      screened_rates(k_p_ti41__v42)*Y(jp)*Y(jti41)*state % rho + &
      screened_rates(k_p_v44__he4_ti41)*Y(jp)*Y(jv44)*state % rho - &
      screened_rates(k_ti41__he4_ca37)*Y(jti41) - &
      screened_rates(k_ti41__p_ca40__weak__wc12)*Y(jti41) - &
      screened_rates(k_ti41__p_sc40)*Y(jti41) - &
      screened_rates(k_ti41__sc41__weak__wc17)*Y(jti41) + &
      screened_rates(k_v41__ti41__weak__bqa_pos_)*Y(jv41) + &
      screened_rates(k_v42__p_ti41)*Y(jv42) &
       )

    ydot_nuc(jti42) = ( &
      screened_rates(k_cr43__p_ti42__weak__wc17)*Y(jcr43) + screened_rates(k_cr46__he4_ti42)* &
      Y(jcr46) + screened_rates(k_fe45__p_p_p_ti42__weak__wc12)*Y(jfe45) + &
      screened_rates(k_he4_ca38__ti42)*Y(jca38)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_sc39__p_ti42)*Y(jhe4)*Y(jsc39)*state % rho - &
      screened_rates(k_he4_ti42__cr46)*Y(jhe4)*Y(jti42)*state % rho - &
      screened_rates(k_he4_ti42__p_v45)*Y(jhe4)*Y(jti42)*state % rho + &
      screened_rates(k_p_sc41__ti42)*Y(jp)*Y(jsc41)*state % rho - &
      screened_rates(k_p_ti42__he4_sc39)*Y(jp)*Y(jti42)*state % rho - &
      screened_rates(k_p_ti42__v43)*Y(jp)*Y(jti42)*state % rho + &
      screened_rates(k_p_v45__he4_ti42)*Y(jp)*Y(jv45)*state % rho - &
      screened_rates(k_ti42__he4_ca38)*Y(jti42) - screened_rates(k_ti42__p_sc41)* &
      Y(jti42) - screened_rates(k_ti42__sc42__weak__wc12)*Y(jti42) + &
      screened_rates(k_v42__ti42__weak__mo97)*Y(jv42) + screened_rates(k_v43__p_ti42)* &
      Y(jv43) &
       )

    ydot_nuc(jti43) = ( &
      screened_rates(k_cr44__p_ti43__weak__wc12)*Y(jcr44) + screened_rates(k_cr47__he4_ti43)* &
      Y(jcr47) + screened_rates(k_he4_ca39__ti43)*Y(jca39)*Y(jhe4)* &
      state % rho + screened_rates(k_he4_sc40__p_ti43)*Y(jhe4)*Y(jsc40)* &
      state % rho - screened_rates(k_he4_ti43__cr47)*Y(jhe4)*Y(jti43)*state % rho &
      - screened_rates(k_he4_ti43__p_v46)*Y(jhe4)*Y(jti43)*state % rho + &
      screened_rates(k_p_sc42__ti43)*Y(jp)*Y(jsc42)*state % rho - &
      screened_rates(k_p_ti43__he4_sc40)*Y(jp)*Y(jti43)*state % rho - &
      screened_rates(k_p_ti43__v44)*Y(jp)*Y(jti43)*state % rho + &
      screened_rates(k_p_v46__he4_ti43)*Y(jp)*Y(jv46)*state % rho - &
      screened_rates(k_ti43__he4_ca39)*Y(jti43) - screened_rates(k_ti43__p_sc42)* &
      Y(jti43) - screened_rates(k_ti43__sc43__weak__wc12)*Y(jti43) + &
      screened_rates(k_v43__ti43__weak__wc12)*Y(jv43) + screened_rates(k_v44__p_ti43)* &
      Y(jv44) &
       )

    ydot_nuc(jti44) = ( &
      screened_rates(k_cr45__p_ti44__weak__wc12)*Y(jcr45) + screened_rates(k_cr48__he4_ti44)* &
      Y(jcr48) + screened_rates(k_he4_ca40__ti44)*Y(jca40)*Y(jhe4)* &
      state % rho + screened_rates(k_he4_sc41__p_ti44)*Y(jhe4)*Y(jsc41)* &
      state % rho - screened_rates(k_he4_ti44__cr48)*Y(jhe4)*Y(jti44)*state % rho &
      - screened_rates(k_he4_ti44__p_v47)*Y(jhe4)*Y(jti44)*state % rho + &
      screened_rates(k_p_sc43__ti44)*Y(jp)*Y(jsc43)*state % rho - &
      screened_rates(k_p_ti44__he4_sc41)*Y(jp)*Y(jti44)*state % rho - &
      screened_rates(k_p_ti44__v45)*Y(jp)*Y(jti44)*state % rho + &
      screened_rates(k_p_v47__he4_ti44)*Y(jp)*Y(jv47)*state % rho - &
      screened_rates(k_ti44__he4_ca40)*Y(jti44) - screened_rates(k_ti44__p_sc43)* &
      Y(jti44) - screened_rates(k_ti44__sc44__weak__wc12)*Y(jti44) + &
      screened_rates(k_v44__ti44__weak__wc12)*Y(jv44) + screened_rates(k_v45__p_ti44)* &
      Y(jv45) &
       )

    ydot_nuc(jti45) = ( &
      screened_rates(k_cr49__he4_ti45)*Y(jcr49) + screened_rates(k_he4_ca41__ti45)*Y(jca41)* &
      Y(jhe4)*state % rho + screened_rates(k_he4_sc42__p_ti45)*Y(jhe4)* &
      Y(jsc42)*state % rho - screened_rates(k_he4_ti45__cr49)*Y(jhe4)* &
      Y(jti45)*state % rho - screened_rates(k_he4_ti45__p_v48)*Y(jhe4)* &
      Y(jti45)*state % rho + screened_rates(k_p_sc44__ti45)*Y(jp)*Y(jsc44)* &
      state % rho - screened_rates(k_p_ti45__he4_sc42)*Y(jp)*Y(jti45)*state % rho &
      - screened_rates(k_p_ti45__v46)*Y(jp)*Y(jti45)*state % rho + &
      screened_rates(k_p_v48__he4_ti45)*Y(jp)*Y(jv48)*state % rho - &
      screened_rates(k_ti45__he4_ca41)*Y(jti45) - screened_rates(k_ti45__p_sc44)* &
      Y(jti45) - screened_rates(k_ti45__sc45__weak__wc12)*Y(jti45) + &
      screened_rates(k_v45__ti45__weak__wc12)*Y(jv45) + screened_rates(k_v46__p_ti45)* &
      Y(jv46) &
       )

    ydot_nuc(jti46) = ( &
      screened_rates(k_cr50__he4_ti46)*Y(jcr50) + screened_rates(k_he4_ca42__ti46)*Y(jca42)* &
      Y(jhe4)*state % rho + screened_rates(k_he4_sc43__p_ti46)*Y(jhe4)* &
      Y(jsc43)*state % rho - screened_rates(k_he4_ti46__cr50)*Y(jhe4)* &
      Y(jti46)*state % rho - screened_rates(k_he4_ti46__p_v49)*Y(jhe4)* &
      Y(jti46)*state % rho + screened_rates(k_p_sc45__ti46)*Y(jp)*Y(jsc45)* &
      state % rho - screened_rates(k_p_ti46__he4_sc43)*Y(jp)*Y(jti46)*state % rho &
      - screened_rates(k_p_ti46__v47)*Y(jp)*Y(jti46)*state % rho + &
      screened_rates(k_p_v49__he4_ti46)*Y(jp)*Y(jv49)*state % rho - &
      screened_rates(k_ti46__he4_ca42)*Y(jti46) - screened_rates(k_ti46__p_sc45)* &
      Y(jti46) + screened_rates(k_v46__ti46__weak__wc12)*Y(jv46) + &
      screened_rates(k_v47__p_ti46)*Y(jv47) &
       )

    ydot_nuc(jti47) = ( &
      screened_rates(k_cr51__he4_ti47)*Y(jcr51) + screened_rates(k_he4_ca43__ti47)*Y(jca43)* &
      Y(jhe4)*state % rho + screened_rates(k_he4_sc44__p_ti47)*Y(jhe4)* &
      Y(jsc44)*state % rho - screened_rates(k_he4_ti47__cr51)*Y(jhe4)* &
      Y(jti47)*state % rho - screened_rates(k_he4_ti47__p_v50)*Y(jhe4)* &
      Y(jti47)*state % rho + screened_rates(k_p_sc46__ti47)*Y(jp)*Y(jsc46)* &
      state % rho - screened_rates(k_p_ti47__he4_sc44)*Y(jp)*Y(jti47)*state % rho &
      - screened_rates(k_p_ti47__v48)*Y(jp)*Y(jti47)*state % rho + &
      screened_rates(k_p_v50__he4_ti47)*Y(jp)*Y(jv50)*state % rho - &
      screened_rates(k_ti47__he4_ca43)*Y(jti47) - screened_rates(k_ti47__p_sc46)* &
      Y(jti47) + screened_rates(k_v47__ti47__weak__wc12)*Y(jv47) + &
      screened_rates(k_v48__p_ti47)*Y(jv48) &
       )

    ydot_nuc(jti48) = ( &
      screened_rates(k_cr52__he4_ti48)*Y(jcr52) + screened_rates(k_he4_ca44__ti48)*Y(jca44)* &
      Y(jhe4)*state % rho + screened_rates(k_he4_sc45__p_ti48)*Y(jhe4)* &
      Y(jsc45)*state % rho - screened_rates(k_he4_ti48__cr52)*Y(jhe4)* &
      Y(jti48)*state % rho - screened_rates(k_p_ti48__he4_sc45)*Y(jp)* &
      Y(jti48)*state % rho - screened_rates(k_p_ti48__v49)*Y(jp)*Y(jti48)* &
      state % rho - screened_rates(k_ti48__he4_ca44)*Y(jti48) + &
      screened_rates(k_v48__ti48__weak__wc12)*Y(jv48) + screened_rates(k_v49__p_ti48)* &
      Y(jv49) &
       )

    ydot_nuc(jv40) = ( &
      screened_rates(k_he4_sc36__v40)*Y(jhe4)*Y(jsc36)*state % rho - &
      screened_rates(k_he4_v40__mn44)*Y(jhe4)*Y(jv40)*state % rho - &
      screened_rates(k_he4_v40__p_cr43)*Y(jhe4)*Y(jv40)*state % rho + &
      screened_rates(k_mn44__he4_v40)*Y(jmn44) + screened_rates(k_p_cr43__he4_v40)* &
      Y(jcr43)*Y(jp)*state % rho + screened_rates(k_p_ti39__v40)*Y(jp)* &
      Y(jti39)*state % rho - screened_rates(k_v40__he4_sc36)*Y(jv40) - &
      screened_rates(k_v40__p_ti39)*Y(jv40) - &
      screened_rates(k_v40__ti40__weak__bqa_pos_)*Y(jv40) &
       )

    ydot_nuc(jv41) = ( &
      screened_rates(k_cr42__p_v41)*Y(jcr42) + screened_rates(k_he4_sc37__v41)*Y(jhe4)* &
      Y(jsc37)*state % rho + screened_rates(k_he4_ti38__p_v41)*Y(jhe4)* &
      Y(jti38)*state % rho - screened_rates(k_he4_v41__mn45)*Y(jhe4)*Y(jv41)* &
      state % rho - screened_rates(k_he4_v41__p_cr44)*Y(jhe4)*Y(jv41)*state % rho &
      + screened_rates(k_mn45__he4_v41)*Y(jmn45) + screened_rates(k_p_cr44__he4_v41)* &
      Y(jcr44)*Y(jp)*state % rho + screened_rates(k_p_ti40__v41)*Y(jp)* &
      Y(jti40)*state % rho - screened_rates(k_p_v41__cr42)*Y(jp)*Y(jv41)* &
      state % rho - screened_rates(k_p_v41__he4_ti38)*Y(jp)*Y(jv41)*state % rho - &
      screened_rates(k_v41__he4_sc37)*Y(jv41) - screened_rates(k_v41__p_ti40)*Y(jv41) &
      - screened_rates(k_v41__ti41__weak__bqa_pos_)*Y(jv41) &
       )

    ydot_nuc(jv42) = ( &
      screened_rates(k_cr42__v42__weak__wc12)*Y(jcr42) + screened_rates(k_cr43__p_v42)* &
      Y(jcr43) + screened_rates(k_he4_sc38__v42)*Y(jhe4)*Y(jsc38)*state % rho &
      + screened_rates(k_he4_ti39__p_v42)*Y(jhe4)*Y(jti39)*state % rho - &
      screened_rates(k_he4_v42__mn46)*Y(jhe4)*Y(jv42)*state % rho - &
      screened_rates(k_he4_v42__p_cr45)*Y(jhe4)*Y(jv42)*state % rho + &
      screened_rates(k_mn46__he4_v42)*Y(jmn46) + screened_rates(k_p_cr45__he4_v42)* &
      Y(jcr45)*Y(jp)*state % rho + screened_rates(k_p_ti41__v42)*Y(jp)* &
      Y(jti41)*state % rho - screened_rates(k_p_v42__cr43)*Y(jp)*Y(jv42)* &
      state % rho - screened_rates(k_p_v42__he4_ti39)*Y(jp)*Y(jv42)*state % rho - &
      screened_rates(k_v42__he4_sc38)*Y(jv42) - screened_rates(k_v42__p_ti41)*Y(jv42) &
      - screened_rates(k_v42__ti42__weak__mo97)*Y(jv42) &
       )

    ydot_nuc(jv43) = ( &
      screened_rates(k_cr43__v43__weak__wc12)*Y(jcr43) + screened_rates(k_cr44__p_v43)* &
      Y(jcr44) + screened_rates(k_fe45__p_p_v43__weak__wc12)*Y(jfe45) + &
      screened_rates(k_he4_sc39__v43)*Y(jhe4)*Y(jsc39)*state % rho + &
      screened_rates(k_he4_ti40__p_v43)*Y(jhe4)*Y(jti40)*state % rho - &
      screened_rates(k_he4_v43__mn47)*Y(jhe4)*Y(jv43)*state % rho - &
      screened_rates(k_he4_v43__p_cr46)*Y(jhe4)*Y(jv43)*state % rho + &
      screened_rates(k_mn47__he4_v43)*Y(jmn47) + screened_rates(k_p_cr46__he4_v43)* &
      Y(jcr46)*Y(jp)*state % rho + screened_rates(k_p_ti42__v43)*Y(jp)* &
      Y(jti42)*state % rho - screened_rates(k_p_v43__cr44)*Y(jp)*Y(jv43)* &
      state % rho - screened_rates(k_p_v43__he4_ti40)*Y(jp)*Y(jv43)*state % rho - &
      screened_rates(k_v43__he4_sc39)*Y(jv43) - screened_rates(k_v43__p_ti42)*Y(jv43) &
      - screened_rates(k_v43__ti43__weak__wc12)*Y(jv43) &
       )

    ydot_nuc(jv44) = ( &
      screened_rates(k_cr44__v44__weak__wc12)*Y(jcr44) + screened_rates(k_cr45__p_v44)* &
      Y(jcr45) + screened_rates(k_he4_sc40__v44)*Y(jhe4)*Y(jsc40)*state % rho &
      + screened_rates(k_he4_ti41__p_v44)*Y(jhe4)*Y(jti41)*state % rho - &
      screened_rates(k_he4_v44__mn48)*Y(jhe4)*Y(jv44)*state % rho - &
      screened_rates(k_he4_v44__p_cr47)*Y(jhe4)*Y(jv44)*state % rho + &
      screened_rates(k_mn48__he4_v44)*Y(jmn48) + screened_rates(k_p_cr47__he4_v44)* &
      Y(jcr47)*Y(jp)*state % rho + screened_rates(k_p_ti43__v44)*Y(jp)* &
      Y(jti43)*state % rho - screened_rates(k_p_v44__cr45)*Y(jp)*Y(jv44)* &
      state % rho - screened_rates(k_p_v44__he4_ti41)*Y(jp)*Y(jv44)*state % rho - &
      screened_rates(k_v44__he4_sc40)*Y(jv44) - screened_rates(k_v44__p_ti43)*Y(jv44) &
      - screened_rates(k_v44__ti44__weak__wc12)*Y(jv44) &
       )

    ydot_nuc(jv45) = ( &
      screened_rates(k_cr45__v45__weak__wc12)*Y(jcr45) + screened_rates(k_cr46__p_v45)* &
      Y(jcr46) + screened_rates(k_he4_sc41__v45)*Y(jhe4)*Y(jsc41)*state % rho &
      + screened_rates(k_he4_ti42__p_v45)*Y(jhe4)*Y(jti42)*state % rho - &
      screened_rates(k_he4_v45__mn49)*Y(jhe4)*Y(jv45)*state % rho - &
      screened_rates(k_he4_v45__p_cr48)*Y(jhe4)*Y(jv45)*state % rho + &
      screened_rates(k_mn46__p_v45__weak__wc12)*Y(jmn46) + &
      screened_rates(k_mn49__he4_v45)*Y(jmn49) + screened_rates(k_p_cr48__he4_v45)* &
      Y(jcr48)*Y(jp)*state % rho + screened_rates(k_p_ti44__v45)*Y(jp)* &
      Y(jti44)*state % rho - screened_rates(k_p_v45__cr46)*Y(jp)*Y(jv45)* &
      state % rho - screened_rates(k_p_v45__he4_ti42)*Y(jp)*Y(jv45)*state % rho - &
      screened_rates(k_v45__he4_sc41)*Y(jv45) - screened_rates(k_v45__p_ti44)*Y(jv45) &
      - screened_rates(k_v45__ti45__weak__wc12)*Y(jv45) &
       )

    ydot_nuc(jv46) = ( &
      screened_rates(k_cr46__v46__weak__wc12)*Y(jcr46) + screened_rates(k_cr47__p_v46)* &
      Y(jcr47) + screened_rates(k_he4_sc42__v46)*Y(jhe4)*Y(jsc42)*state % rho &
      + screened_rates(k_he4_ti43__p_v46)*Y(jhe4)*Y(jti43)*state % rho - &
      screened_rates(k_he4_v46__mn50)*Y(jhe4)*Y(jv46)*state % rho - &
      screened_rates(k_he4_v46__p_cr49)*Y(jhe4)*Y(jv46)*state % rho + &
      screened_rates(k_mn50__he4_v46)*Y(jmn50) + screened_rates(k_p_cr49__he4_v46)* &
      Y(jcr49)*Y(jp)*state % rho + screened_rates(k_p_ti45__v46)*Y(jp)* &
      Y(jti45)*state % rho - screened_rates(k_p_v46__cr47)*Y(jp)*Y(jv46)* &
      state % rho - screened_rates(k_p_v46__he4_ti43)*Y(jp)*Y(jv46)*state % rho - &
      screened_rates(k_v46__he4_sc42)*Y(jv46) - screened_rates(k_v46__p_ti45)*Y(jv46) &
      - screened_rates(k_v46__ti46__weak__wc12)*Y(jv46) &
       )

    ydot_nuc(jv47) = ( &
      screened_rates(k_cr47__v47__weak__wc12)*Y(jcr47) + screened_rates(k_cr48__p_v47)* &
      Y(jcr48) + screened_rates(k_he4_sc43__v47)*Y(jhe4)*Y(jsc43)*state % rho &
      + screened_rates(k_he4_ti44__p_v47)*Y(jhe4)*Y(jti44)*state % rho - &
      screened_rates(k_he4_v47__mn51)*Y(jhe4)*Y(jv47)*state % rho - &
      screened_rates(k_he4_v47__p_cr50)*Y(jhe4)*Y(jv47)*state % rho + &
      screened_rates(k_mn48__p_v47__weak__wc12)*Y(jmn48) + &
      screened_rates(k_mn51__he4_v47)*Y(jmn51) + screened_rates(k_p_cr50__he4_v47)* &
      Y(jcr50)*Y(jp)*state % rho + screened_rates(k_p_ti46__v47)*Y(jp)* &
      Y(jti46)*state % rho - screened_rates(k_p_v47__cr48)*Y(jp)*Y(jv47)* &
      state % rho - screened_rates(k_p_v47__he4_ti44)*Y(jp)*Y(jv47)*state % rho - &
      screened_rates(k_v47__he4_sc43)*Y(jv47) - screened_rates(k_v47__p_ti46)*Y(jv47) &
      - screened_rates(k_v47__ti47__weak__wc12)*Y(jv47) &
       )

    ydot_nuc(jv48) = ( &
      screened_rates(k_cr48__v48__weak__wc12)*Y(jcr48) + screened_rates(k_cr49__p_v48)* &
      Y(jcr49) + screened_rates(k_he4_sc44__v48)*Y(jhe4)*Y(jsc44)*state % rho &
      + screened_rates(k_he4_ti45__p_v48)*Y(jhe4)*Y(jti45)*state % rho - &
      screened_rates(k_he4_v48__mn52)*Y(jhe4)*Y(jv48)*state % rho - &
      screened_rates(k_he4_v48__p_cr51)*Y(jhe4)*Y(jv48)*state % rho + &
      screened_rates(k_mn52__he4_v48)*Y(jmn52) + screened_rates(k_p_cr51__he4_v48)* &
      Y(jcr51)*Y(jp)*state % rho + screened_rates(k_p_ti47__v48)*Y(jp)* &
      Y(jti47)*state % rho - screened_rates(k_p_v48__cr49)*Y(jp)*Y(jv48)* &
      state % rho - screened_rates(k_p_v48__he4_ti45)*Y(jp)*Y(jv48)*state % rho - &
      screened_rates(k_v48__he4_sc44)*Y(jv48) - screened_rates(k_v48__p_ti47)*Y(jv48) &
      - screened_rates(k_v48__ti48__weak__wc12)*Y(jv48) &
       )

    ydot_nuc(jv49) = ( &
      screened_rates(k_cr49__v49__weak__wc12)*Y(jcr49) + screened_rates(k_cr50__p_v49)* &
      Y(jcr50) + screened_rates(k_he4_sc45__v49)*Y(jhe4)*Y(jsc45)*state % rho &
      + screened_rates(k_he4_ti46__p_v49)*Y(jhe4)*Y(jti46)*state % rho - &
      screened_rates(k_he4_v49__mn53)*Y(jhe4)*Y(jv49)*state % rho - &
      screened_rates(k_he4_v49__p_cr52)*Y(jhe4)*Y(jv49)*state % rho + &
      screened_rates(k_mn53__he4_v49)*Y(jmn53) + screened_rates(k_p_cr52__he4_v49)* &
      Y(jcr52)*Y(jp)*state % rho + screened_rates(k_p_ti48__v49)*Y(jp)* &
      Y(jti48)*state % rho - screened_rates(k_p_v49__cr50)*Y(jp)*Y(jv49)* &
      state % rho - screened_rates(k_p_v49__he4_ti46)*Y(jp)*Y(jv49)*state % rho - &
      screened_rates(k_v49__he4_sc45)*Y(jv49) - screened_rates(k_v49__p_ti48)*Y(jv49) &
       )

    ydot_nuc(jv50) = ( &
      screened_rates(k_cr51__p_v50)*Y(jcr51) + screened_rates(k_he4_sc46__v50)*Y(jhe4)* &
      Y(jsc46)*state % rho + screened_rates(k_he4_ti47__p_v50)*Y(jhe4)* &
      Y(jti47)*state % rho - screened_rates(k_he4_v50__mn54)*Y(jhe4)*Y(jv50)* &
      state % rho + screened_rates(k_mn54__he4_v50)*Y(jmn54) - &
      screened_rates(k_p_v50__cr51)*Y(jp)*Y(jv50)*state % rho - &
      screened_rates(k_p_v50__he4_ti47)*Y(jp)*Y(jv50)*state % rho - &
      screened_rates(k_v50__he4_sc46)*Y(jv50) &
       )

    ydot_nuc(jcr42) = ( &
      -screened_rates(k_cr42__he4_ti38)*Y(jcr42) - screened_rates(k_cr42__p_ti41__weak__wc12)* &
      Y(jcr42) - screened_rates(k_cr42__p_v41)*Y(jcr42) - &
      screened_rates(k_cr42__v42__weak__wc12)*Y(jcr42) + &
      screened_rates(k_fe46__he4_cr42)*Y(jfe46) - screened_rates(k_he4_cr42__fe46)* &
      Y(jcr42)*Y(jhe4)*state % rho - screened_rates(k_he4_cr42__p_mn45)* &
      Y(jcr42)*Y(jhe4)*state % rho + screened_rates(k_he4_ti38__cr42)*Y(jhe4) &
      *Y(jti38)*state % rho + screened_rates(k_p_mn45__he4_cr42)*Y(jmn45)* &
      Y(jp)*state % rho + screened_rates(k_p_v41__cr42)*Y(jp)*Y(jv41)* &
      state % rho &
       )

    ydot_nuc(jcr43) = ( &
      -screened_rates(k_cr43__he4_ti39)*Y(jcr43) - &
      screened_rates(k_cr43__p_p_p_ca40__weak__wc12)*Y(jcr43) - &
      screened_rates(k_cr43__p_p_sc41__weak__wc12)*Y(jcr43) - &
      screened_rates(k_cr43__p_ti42__weak__wc17)*Y(jcr43) - &
      screened_rates(k_cr43__p_v42)*Y(jcr43) - screened_rates(k_cr43__v43__weak__wc12)* &
      Y(jcr43) + screened_rates(k_fe47__he4_cr43)*Y(jfe47) - &
      screened_rates(k_he4_cr43__fe47)*Y(jcr43)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr43__p_mn46)*Y(jcr43)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ti39__cr43)*Y(jhe4)*Y(jti39)*state % rho + &
      screened_rates(k_he4_v40__p_cr43)*Y(jhe4)*Y(jv40)*state % rho + &
      screened_rates(k_mn44__p_cr43)*Y(jmn44) - screened_rates(k_p_cr43__he4_v40)* &
      Y(jcr43)*Y(jp)*state % rho - screened_rates(k_p_cr43__mn44)*Y(jcr43)* &
      Y(jp)*state % rho + screened_rates(k_p_mn46__he4_cr43)*Y(jmn46)*Y(jp)* &
      state % rho + screened_rates(k_p_v42__cr43)*Y(jp)*Y(jv42)*state % rho &
       )

    ydot_nuc(jcr44) = ( &
      -screened_rates(k_cr44__he4_ti40)*Y(jcr44) - screened_rates(k_cr44__p_ti43__weak__wc12)* &
      Y(jcr44) - screened_rates(k_cr44__p_v43)*Y(jcr44) - &
      screened_rates(k_cr44__v44__weak__wc12)*Y(jcr44) + &
      screened_rates(k_fe45__p_cr44__weak__wc17)*Y(jfe45) + &
      screened_rates(k_fe48__he4_cr44)*Y(jfe48) - screened_rates(k_he4_cr44__fe48)* &
      Y(jcr44)*Y(jhe4)*state % rho - screened_rates(k_he4_cr44__p_mn47)* &
      Y(jcr44)*Y(jhe4)*state % rho + screened_rates(k_he4_ti40__cr44)*Y(jhe4) &
      *Y(jti40)*state % rho + screened_rates(k_he4_v41__p_cr44)*Y(jhe4)* &
      Y(jv41)*state % rho + screened_rates(k_mn44__cr44__weak__bqa_pos_)* &
      Y(jmn44) + screened_rates(k_mn45__p_cr44)*Y(jmn45) - &
      screened_rates(k_p_cr44__he4_v41)*Y(jcr44)*Y(jp)*state % rho - &
      screened_rates(k_p_cr44__mn45)*Y(jcr44)*Y(jp)*state % rho + &
      screened_rates(k_p_mn47__he4_cr44)*Y(jmn47)*Y(jp)*state % rho + &
      screened_rates(k_p_v43__cr44)*Y(jp)*Y(jv43)*state % rho &
       )

    ydot_nuc(jcr45) = ( &
      -screened_rates(k_cr45__he4_ti41)*Y(jcr45) - screened_rates(k_cr45__p_ti44__weak__wc12)* &
      Y(jcr45) - screened_rates(k_cr45__p_v44)*Y(jcr45) - &
      screened_rates(k_cr45__v45__weak__wc12)*Y(jcr45) + &
      screened_rates(k_fe46__p_cr45__weak__wc12)*Y(jfe46) + &
      screened_rates(k_fe49__he4_cr45)*Y(jfe49) - screened_rates(k_he4_cr45__fe49)* &
      Y(jcr45)*Y(jhe4)*state % rho - screened_rates(k_he4_cr45__p_mn48)* &
      Y(jcr45)*Y(jhe4)*state % rho + screened_rates(k_he4_ti41__cr45)*Y(jhe4) &
      *Y(jti41)*state % rho + screened_rates(k_he4_v42__p_cr45)*Y(jhe4)* &
      Y(jv42)*state % rho + screened_rates(k_mn45__cr45__weak__bqa_pos_)* &
      Y(jmn45) + screened_rates(k_mn46__p_cr45)*Y(jmn46) - &
      screened_rates(k_p_cr45__he4_v42)*Y(jcr45)*Y(jp)*state % rho - &
      screened_rates(k_p_cr45__mn46)*Y(jcr45)*Y(jp)*state % rho + &
      screened_rates(k_p_mn48__he4_cr45)*Y(jmn48)*Y(jp)*state % rho + &
      screened_rates(k_p_v44__cr45)*Y(jp)*Y(jv44)*state % rho &
       )

    ydot_nuc(jcr46) = ( &
      -screened_rates(k_cr46__he4_ti42)*Y(jcr46) - screened_rates(k_cr46__p_v45)*Y(jcr46) - &
      screened_rates(k_cr46__v46__weak__wc12)*Y(jcr46) + &
      screened_rates(k_fe47__p_cr46__weak__wc12)*Y(jfe47) + &
      screened_rates(k_fe50__he4_cr46)*Y(jfe50) - screened_rates(k_he4_cr46__fe50)* &
      Y(jcr46)*Y(jhe4)*state % rho - screened_rates(k_he4_cr46__p_mn49)* &
      Y(jcr46)*Y(jhe4)*state % rho + screened_rates(k_he4_ti42__cr46)*Y(jhe4) &
      *Y(jti42)*state % rho + screened_rates(k_he4_v43__p_cr46)*Y(jhe4)* &
      Y(jv43)*state % rho + screened_rates(k_mn46__cr46__weak__wc12)*Y(jmn46) + &
      screened_rates(k_mn47__p_cr46)*Y(jmn47) - screened_rates(k_p_cr46__he4_v43)* &
      Y(jcr46)*Y(jp)*state % rho - screened_rates(k_p_cr46__mn47)*Y(jcr46)* &
      Y(jp)*state % rho + screened_rates(k_p_mn49__he4_cr46)*Y(jmn49)*Y(jp)* &
      state % rho + screened_rates(k_p_v45__cr46)*Y(jp)*Y(jv45)*state % rho &
       )

    ydot_nuc(jcr47) = ( &
      -screened_rates(k_cr47__he4_ti43)*Y(jcr47) - screened_rates(k_cr47__p_v46)*Y(jcr47) - &
      screened_rates(k_cr47__v47__weak__wc12)*Y(jcr47) + &
      screened_rates(k_fe48__p_cr47__weak__wc12)*Y(jfe48) + &
      screened_rates(k_fe51__he4_cr47)*Y(jfe51) - screened_rates(k_he4_cr47__fe51)* &
      Y(jcr47)*Y(jhe4)*state % rho - screened_rates(k_he4_cr47__p_mn50)* &
      Y(jcr47)*Y(jhe4)*state % rho + screened_rates(k_he4_ti43__cr47)*Y(jhe4) &
      *Y(jti43)*state % rho + screened_rates(k_he4_v44__p_cr47)*Y(jhe4)* &
      Y(jv44)*state % rho + screened_rates(k_mn47__cr47__weak__wc17)*Y(jmn47) + &
      screened_rates(k_mn48__p_cr47)*Y(jmn48) - screened_rates(k_p_cr47__he4_v44)* &
      Y(jcr47)*Y(jp)*state % rho - screened_rates(k_p_cr47__mn48)*Y(jcr47)* &
      Y(jp)*state % rho + screened_rates(k_p_mn50__he4_cr47)*Y(jmn50)*Y(jp)* &
      state % rho + screened_rates(k_p_v46__cr47)*Y(jp)*Y(jv46)*state % rho &
       )

    ydot_nuc(jcr48) = ( &
      -screened_rates(k_cr48__he4_ti44)*Y(jcr48) - screened_rates(k_cr48__p_v47)*Y(jcr48) - &
      screened_rates(k_cr48__v48__weak__wc12)*Y(jcr48) + &
      screened_rates(k_fe49__p_cr48__weak__wc12)*Y(jfe49) + &
      screened_rates(k_fe52__he4_cr48)*Y(jfe52) - screened_rates(k_he4_cr48__fe52)* &
      Y(jcr48)*Y(jhe4)*state % rho - screened_rates(k_he4_cr48__p_mn51)* &
      Y(jcr48)*Y(jhe4)*state % rho + screened_rates(k_he4_ti44__cr48)*Y(jhe4) &
      *Y(jti44)*state % rho + screened_rates(k_he4_v45__p_cr48)*Y(jhe4)* &
      Y(jv45)*state % rho + screened_rates(k_mn48__cr48__weak__wc12)*Y(jmn48) + &
      screened_rates(k_mn49__p_cr48)*Y(jmn49) - screened_rates(k_p_cr48__he4_v45)* &
      Y(jcr48)*Y(jp)*state % rho - screened_rates(k_p_cr48__mn49)*Y(jcr48)* &
      Y(jp)*state % rho + screened_rates(k_p_mn51__he4_cr48)*Y(jmn51)*Y(jp)* &
      state % rho + screened_rates(k_p_v47__cr48)*Y(jp)*Y(jv47)*state % rho &
       )

    ydot_nuc(jcr49) = ( &
      -screened_rates(k_cr49__he4_ti45)*Y(jcr49) - screened_rates(k_cr49__p_v48)*Y(jcr49) - &
      screened_rates(k_cr49__v49__weak__wc12)*Y(jcr49) + &
      screened_rates(k_fe53__he4_cr49)*Y(jfe53) - screened_rates(k_he4_cr49__fe53)* &
      Y(jcr49)*Y(jhe4)*state % rho - screened_rates(k_he4_cr49__p_mn52)* &
      Y(jcr49)*Y(jhe4)*state % rho + screened_rates(k_he4_ti45__cr49)*Y(jhe4) &
      *Y(jti45)*state % rho + screened_rates(k_he4_v46__p_cr49)*Y(jhe4)* &
      Y(jv46)*state % rho + screened_rates(k_mn49__cr49__weak__wc12)*Y(jmn49) + &
      screened_rates(k_mn50__p_cr49)*Y(jmn50) - screened_rates(k_p_cr49__he4_v46)* &
      Y(jcr49)*Y(jp)*state % rho - screened_rates(k_p_cr49__mn50)*Y(jcr49)* &
      Y(jp)*state % rho + screened_rates(k_p_mn52__he4_cr49)*Y(jmn52)*Y(jp)* &
      state % rho + screened_rates(k_p_v48__cr49)*Y(jp)*Y(jv48)*state % rho &
       )

    ydot_nuc(jcr50) = ( &
      -screened_rates(k_cr50__he4_ti46)*Y(jcr50) - screened_rates(k_cr50__p_v49)*Y(jcr50) + &
      screened_rates(k_fe54__he4_cr50)*Y(jfe54) - screened_rates(k_he4_cr50__fe54)* &
      Y(jcr50)*Y(jhe4)*state % rho - screened_rates(k_he4_cr50__p_mn53)* &
      Y(jcr50)*Y(jhe4)*state % rho + screened_rates(k_he4_ti46__cr50)*Y(jhe4) &
      *Y(jti46)*state % rho + screened_rates(k_he4_v47__p_cr50)*Y(jhe4)* &
      Y(jv47)*state % rho + screened_rates(k_mn50__cr50__weak__wc12)*Y(jmn50) + &
      screened_rates(k_mn51__p_cr50)*Y(jmn51) - screened_rates(k_p_cr50__he4_v47)* &
      Y(jcr50)*Y(jp)*state % rho - screened_rates(k_p_cr50__mn51)*Y(jcr50)* &
      Y(jp)*state % rho + screened_rates(k_p_mn53__he4_cr50)*Y(jmn53)*Y(jp)* &
      state % rho + screened_rates(k_p_v49__cr50)*Y(jp)*Y(jv49)*state % rho &
       )

    ydot_nuc(jcr51) = ( &
      -screened_rates(k_cr51__he4_ti47)*Y(jcr51) - screened_rates(k_cr51__p_v50)*Y(jcr51) + &
      screened_rates(k_fe55__he4_cr51)*Y(jfe55) - screened_rates(k_he4_cr51__fe55)* &
      Y(jcr51)*Y(jhe4)*state % rho - screened_rates(k_he4_cr51__p_mn54)* &
      Y(jcr51)*Y(jhe4)*state % rho + screened_rates(k_he4_ti47__cr51)*Y(jhe4) &
      *Y(jti47)*state % rho + screened_rates(k_he4_v48__p_cr51)*Y(jhe4)* &
      Y(jv48)*state % rho + screened_rates(k_mn51__cr51__weak__wc12)*Y(jmn51) + &
      screened_rates(k_mn52__p_cr51)*Y(jmn52) - screened_rates(k_p_cr51__he4_v48)* &
      Y(jcr51)*Y(jp)*state % rho - screened_rates(k_p_cr51__mn52)*Y(jcr51)* &
      Y(jp)*state % rho + screened_rates(k_p_mn54__he4_cr51)*Y(jmn54)*Y(jp)* &
      state % rho + screened_rates(k_p_v50__cr51)*Y(jp)*Y(jv50)*state % rho &
       )

    ydot_nuc(jcr52) = ( &
      -screened_rates(k_cr52__he4_ti48)*Y(jcr52) + screened_rates(k_fe56__he4_cr52)*Y(jfe56) &
      - screened_rates(k_he4_cr52__fe56)*Y(jcr52)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr52__p_mn55)*Y(jcr52)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ti48__cr52)*Y(jhe4)*Y(jti48)*state % rho + &
      screened_rates(k_he4_v49__p_cr52)*Y(jhe4)*Y(jv49)*state % rho + &
      screened_rates(k_mn52__cr52__weak__wc12)*Y(jmn52) + &
      screened_rates(k_mn53__p_cr52)*Y(jmn53) - screened_rates(k_p_cr52__he4_v49)* &
      Y(jcr52)*Y(jp)*state % rho - screened_rates(k_p_cr52__mn53)*Y(jcr52)* &
      Y(jp)*state % rho + screened_rates(k_p_mn55__he4_cr52)*Y(jmn55)*Y(jp)* &
      state % rho &
       )

    ydot_nuc(jmn44) = ( &
      screened_rates(k_co48__he4_mn44)*Y(jco48) + screened_rates(k_fe45__p_mn44)*Y(jfe45) - &
      screened_rates(k_he4_mn44__co48)*Y(jhe4)*Y(jmn44)*state % rho - &
      screened_rates(k_he4_mn44__p_fe47)*Y(jhe4)*Y(jmn44)*state % rho + &
      screened_rates(k_he4_v40__mn44)*Y(jhe4)*Y(jv40)*state % rho - &
      screened_rates(k_mn44__cr44__weak__bqa_pos_)*Y(jmn44) - &
      screened_rates(k_mn44__he4_v40)*Y(jmn44) - screened_rates(k_mn44__p_cr43)* &
      Y(jmn44) + screened_rates(k_p_cr43__mn44)*Y(jcr43)*Y(jp)*state % rho + &
      screened_rates(k_p_fe47__he4_mn44)*Y(jfe47)*Y(jp)*state % rho - &
      screened_rates(k_p_mn44__fe45)*Y(jmn44)*Y(jp)*state % rho &
       )

    ydot_nuc(jmn45) = ( &
      screened_rates(k_co49__he4_mn45)*Y(jco49) + screened_rates(k_fe45__mn45__weak__wc17)* &
      Y(jfe45) + screened_rates(k_fe46__p_mn45)*Y(jfe46) + &
      screened_rates(k_he4_cr42__p_mn45)*Y(jcr42)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_mn45__co49)*Y(jhe4)*Y(jmn45)*state % rho - &
      screened_rates(k_he4_mn45__p_fe48)*Y(jhe4)*Y(jmn45)*state % rho + &
      screened_rates(k_he4_v41__mn45)*Y(jhe4)*Y(jv41)*state % rho - &
      screened_rates(k_mn45__cr45__weak__bqa_pos_)*Y(jmn45) - &
      screened_rates(k_mn45__he4_v41)*Y(jmn45) - screened_rates(k_mn45__p_cr44)* &
      Y(jmn45) + screened_rates(k_p_cr44__mn45)*Y(jcr44)*Y(jp)*state % rho + &
      screened_rates(k_p_fe48__he4_mn45)*Y(jfe48)*Y(jp)*state % rho - &
      screened_rates(k_p_mn45__fe46)*Y(jmn45)*Y(jp)*state % rho - &
      screened_rates(k_p_mn45__he4_cr42)*Y(jmn45)*Y(jp)*state % rho &
       )

    ydot_nuc(jmn46) = ( &
      screened_rates(k_co50__he4_mn46)*Y(jco50) + screened_rates(k_fe46__mn46__weak__wc12)* &
      Y(jfe46) + screened_rates(k_fe47__p_mn46)*Y(jfe47) + &
      screened_rates(k_he4_cr43__p_mn46)*Y(jcr43)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_mn46__co50)*Y(jhe4)*Y(jmn46)*state % rho - &
      screened_rates(k_he4_mn46__p_fe49)*Y(jhe4)*Y(jmn46)*state % rho + &
      screened_rates(k_he4_v42__mn46)*Y(jhe4)*Y(jv42)*state % rho - &
      screened_rates(k_mn46__cr46__weak__wc12)*Y(jmn46) - &
      screened_rates(k_mn46__he4_v42)*Y(jmn46) - screened_rates(k_mn46__p_cr45)* &
      Y(jmn46) - screened_rates(k_mn46__p_v45__weak__wc12)*Y(jmn46) + &
      screened_rates(k_p_cr45__mn46)*Y(jcr45)*Y(jp)*state % rho + &
      screened_rates(k_p_fe49__he4_mn46)*Y(jfe49)*Y(jp)*state % rho - &
      screened_rates(k_p_mn46__fe47)*Y(jmn46)*Y(jp)*state % rho - &
      screened_rates(k_p_mn46__he4_cr43)*Y(jmn46)*Y(jp)*state % rho &
       )

    ydot_nuc(jmn47) = ( &
      screened_rates(k_co51__he4_mn47)*Y(jco51) + screened_rates(k_fe47__mn47__weak__wc12)* &
      Y(jfe47) + screened_rates(k_fe48__p_mn47)*Y(jfe48) + &
      screened_rates(k_he4_cr44__p_mn47)*Y(jcr44)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_mn47__co51)*Y(jhe4)*Y(jmn47)*state % rho - &
      screened_rates(k_he4_mn47__p_fe50)*Y(jhe4)*Y(jmn47)*state % rho + &
      screened_rates(k_he4_v43__mn47)*Y(jhe4)*Y(jv43)*state % rho - &
      screened_rates(k_mn47__cr47__weak__wc17)*Y(jmn47) - &
      screened_rates(k_mn47__he4_v43)*Y(jmn47) - screened_rates(k_mn47__p_cr46)* &
      Y(jmn47) + screened_rates(k_p_cr46__mn47)*Y(jcr46)*Y(jp)*state % rho + &
      screened_rates(k_p_fe50__he4_mn47)*Y(jfe50)*Y(jp)*state % rho - &
      screened_rates(k_p_mn47__fe48)*Y(jmn47)*Y(jp)*state % rho - &
      screened_rates(k_p_mn47__he4_cr44)*Y(jmn47)*Y(jp)*state % rho &
       )

    ydot_nuc(jmn48) = ( &
      screened_rates(k_co52__he4_mn48)*Y(jco52) + screened_rates(k_fe48__mn48__weak__wc12)* &
      Y(jfe48) + screened_rates(k_fe49__p_mn48)*Y(jfe49) + &
      screened_rates(k_he4_cr45__p_mn48)*Y(jcr45)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_mn48__co52)*Y(jhe4)*Y(jmn48)*state % rho - &
      screened_rates(k_he4_mn48__p_fe51)*Y(jhe4)*Y(jmn48)*state % rho + &
      screened_rates(k_he4_v44__mn48)*Y(jhe4)*Y(jv44)*state % rho - &
      screened_rates(k_mn48__cr48__weak__wc12)*Y(jmn48) - &
      screened_rates(k_mn48__he4_v44)*Y(jmn48) - screened_rates(k_mn48__p_cr47)* &
      Y(jmn48) - screened_rates(k_mn48__p_v47__weak__wc12)*Y(jmn48) + &
      screened_rates(k_p_cr47__mn48)*Y(jcr47)*Y(jp)*state % rho + &
      screened_rates(k_p_fe51__he4_mn48)*Y(jfe51)*Y(jp)*state % rho - &
      screened_rates(k_p_mn48__fe49)*Y(jmn48)*Y(jp)*state % rho - &
      screened_rates(k_p_mn48__he4_cr45)*Y(jmn48)*Y(jp)*state % rho &
       )

    ydot_nuc(jmn49) = ( &
      screened_rates(k_co50__p_mn49__weak__wc12)*Y(jco50) + screened_rates(k_co53__he4_mn49)* &
      Y(jco53) + screened_rates(k_fe49__mn49__weak__wc12)*Y(jfe49) + &
      screened_rates(k_fe50__p_mn49)*Y(jfe50) + screened_rates(k_he4_cr46__p_mn49)* &
      Y(jcr46)*Y(jhe4)*state % rho - screened_rates(k_he4_mn49__co53)*Y(jhe4) &
      *Y(jmn49)*state % rho - screened_rates(k_he4_mn49__p_fe52)*Y(jhe4)* &
      Y(jmn49)*state % rho + screened_rates(k_he4_v45__mn49)*Y(jhe4)*Y(jv45)* &
      state % rho - screened_rates(k_mn49__cr49__weak__wc12)*Y(jmn49) - &
      screened_rates(k_mn49__he4_v45)*Y(jmn49) - screened_rates(k_mn49__p_cr48)* &
      Y(jmn49) + screened_rates(k_p_cr48__mn49)*Y(jcr48)*Y(jp)*state % rho + &
      screened_rates(k_p_fe52__he4_mn49)*Y(jfe52)*Y(jp)*state % rho - &
      screened_rates(k_p_mn49__fe50)*Y(jmn49)*Y(jp)*state % rho - &
      screened_rates(k_p_mn49__he4_cr46)*Y(jmn49)*Y(jp)*state % rho &
       )

    ydot_nuc(jmn50) = ( &
      screened_rates(k_co54__he4_mn50)*Y(jco54) + screened_rates(k_fe50__mn50__weak__wc12)* &
      Y(jfe50) + screened_rates(k_fe51__p_mn50)*Y(jfe51) + &
      screened_rates(k_he4_cr47__p_mn50)*Y(jcr47)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_mn50__co54)*Y(jhe4)*Y(jmn50)*state % rho - &
      screened_rates(k_he4_mn50__p_fe53)*Y(jhe4)*Y(jmn50)*state % rho + &
      screened_rates(k_he4_v46__mn50)*Y(jhe4)*Y(jv46)*state % rho - &
      screened_rates(k_mn50__cr50__weak__wc12)*Y(jmn50) - &
      screened_rates(k_mn50__he4_v46)*Y(jmn50) - screened_rates(k_mn50__p_cr49)* &
      Y(jmn50) + screened_rates(k_p_cr49__mn50)*Y(jcr49)*Y(jp)*state % rho + &
      screened_rates(k_p_fe53__he4_mn50)*Y(jfe53)*Y(jp)*state % rho - &
      screened_rates(k_p_mn50__fe51)*Y(jmn50)*Y(jp)*state % rho - &
      screened_rates(k_p_mn50__he4_cr47)*Y(jmn50)*Y(jp)*state % rho &
       )

    ydot_nuc(jmn51) = ( &
      screened_rates(k_co55__he4_mn51)*Y(jco55) + screened_rates(k_fe51__mn51__weak__wc12)* &
      Y(jfe51) + screened_rates(k_fe52__p_mn51)*Y(jfe52) + &
      screened_rates(k_he4_cr48__p_mn51)*Y(jcr48)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_mn51__co55)*Y(jhe4)*Y(jmn51)*state % rho - &
      screened_rates(k_he4_mn51__p_fe54)*Y(jhe4)*Y(jmn51)*state % rho + &
      screened_rates(k_he4_v47__mn51)*Y(jhe4)*Y(jv47)*state % rho - &
      screened_rates(k_mn51__cr51__weak__wc12)*Y(jmn51) - &
      screened_rates(k_mn51__he4_v47)*Y(jmn51) - screened_rates(k_mn51__p_cr50)* &
      Y(jmn51) + screened_rates(k_p_cr50__mn51)*Y(jcr50)*Y(jp)*state % rho + &
      screened_rates(k_p_fe54__he4_mn51)*Y(jfe54)*Y(jp)*state % rho - &
      screened_rates(k_p_mn51__fe52)*Y(jmn51)*Y(jp)*state % rho - &
      screened_rates(k_p_mn51__he4_cr48)*Y(jmn51)*Y(jp)*state % rho &
       )

    ydot_nuc(jmn52) = ( &
      screened_rates(k_co56__he4_mn52)*Y(jco56) + screened_rates(k_fe52__mn52__weak__wc12)* &
      Y(jfe52) + screened_rates(k_fe53__p_mn52)*Y(jfe53) + &
      screened_rates(k_he4_cr49__p_mn52)*Y(jcr49)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_mn52__co56)*Y(jhe4)*Y(jmn52)*state % rho - &
      screened_rates(k_he4_mn52__p_fe55)*Y(jhe4)*Y(jmn52)*state % rho + &
      screened_rates(k_he4_v48__mn52)*Y(jhe4)*Y(jv48)*state % rho - &
      screened_rates(k_mn52__cr52__weak__wc12)*Y(jmn52) - &
      screened_rates(k_mn52__he4_v48)*Y(jmn52) - screened_rates(k_mn52__p_cr51)* &
      Y(jmn52) + screened_rates(k_p_cr51__mn52)*Y(jcr51)*Y(jp)*state % rho + &
      screened_rates(k_p_fe55__he4_mn52)*Y(jfe55)*Y(jp)*state % rho - &
      screened_rates(k_p_mn52__fe53)*Y(jmn52)*Y(jp)*state % rho - &
      screened_rates(k_p_mn52__he4_cr49)*Y(jmn52)*Y(jp)*state % rho &
       )

    ydot_nuc(jmn53) = ( &
      screened_rates(k_fe53__mn53__weak__wc12)*Y(jfe53) + screened_rates(k_fe54__p_mn53)* &
      Y(jfe54) + screened_rates(k_he4_cr50__p_mn53)*Y(jcr50)*Y(jhe4)* &
      state % rho - screened_rates(k_he4_mn53__p_fe56)*Y(jhe4)*Y(jmn53)* &
      state % rho + screened_rates(k_he4_v49__mn53)*Y(jhe4)*Y(jv49)*state % rho - &
      screened_rates(k_mn53__he4_v49)*Y(jmn53) - screened_rates(k_mn53__p_cr52)* &
      Y(jmn53) + screened_rates(k_p_cr52__mn53)*Y(jcr52)*Y(jp)*state % rho + &
      screened_rates(k_p_fe56__he4_mn53)*Y(jfe56)*Y(jp)*state % rho - &
      screened_rates(k_p_mn53__fe54)*Y(jmn53)*Y(jp)*state % rho - &
      screened_rates(k_p_mn53__he4_cr50)*Y(jmn53)*Y(jp)*state % rho &
       )

    ydot_nuc(jmn54) = ( &
      screened_rates(k_fe55__p_mn54)*Y(jfe55) + screened_rates(k_he4_cr51__p_mn54)*Y(jcr51)* &
      Y(jhe4)*state % rho + screened_rates(k_he4_v50__mn54)*Y(jhe4)*Y(jv50)* &
      state % rho - screened_rates(k_mn54__he4_v50)*Y(jmn54) - &
      screened_rates(k_p_mn54__fe55)*Y(jmn54)*Y(jp)*state % rho - &
      screened_rates(k_p_mn54__he4_cr51)*Y(jmn54)*Y(jp)*state % rho &
       )

    ydot_nuc(jmn55) = ( &
      screened_rates(k_fe55__mn55__weak__wc12)*Y(jfe55) + screened_rates(k_fe56__p_mn55)* &
      Y(jfe56) + screened_rates(k_he4_cr52__p_mn55)*Y(jcr52)*Y(jhe4)* &
      state % rho - screened_rates(k_p_mn55__fe56)*Y(jmn55)*Y(jp)*state % rho - &
      screened_rates(k_p_mn55__he4_cr52)*Y(jmn55)*Y(jp)*state % rho &
       )

    ydot_nuc(jfe45) = ( &
      -screened_rates(k_fe45__mn45__weak__wc17)*Y(jfe45) - &
      screened_rates(k_fe45__p_cr44__weak__wc17)*Y(jfe45) - &
      screened_rates(k_fe45__p_mn44)*Y(jfe45) - &
      screened_rates(k_fe45__p_p_p_ti42__weak__wc12)*Y(jfe45) - &
      screened_rates(k_fe45__p_p_v43__weak__wc12)*Y(jfe45) - &
      screened_rates(k_he4_fe45__ni49)*Y(jfe45)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe45__p_co48)*Y(jfe45)*Y(jhe4)*state % rho + &
      screened_rates(k_ni49__he4_fe45)*Y(jni49) + screened_rates(k_p_co48__he4_fe45)* &
      Y(jco48)*Y(jp)*state % rho + screened_rates(k_p_mn44__fe45)*Y(jmn44)* &
      Y(jp)*state % rho &
       )

    ydot_nuc(jfe46) = ( &
      screened_rates(k_co47__p_fe46)*Y(jco47) - screened_rates(k_fe46__he4_cr42)*Y(jfe46) - &
      screened_rates(k_fe46__mn46__weak__wc12)*Y(jfe46) - &
      screened_rates(k_fe46__p_cr45__weak__wc12)*Y(jfe46) - &
      screened_rates(k_fe46__p_mn45)*Y(jfe46) + screened_rates(k_he4_cr42__fe46)* &
      Y(jcr42)*Y(jhe4)*state % rho - screened_rates(k_he4_fe46__ni50)* &
      Y(jfe46)*Y(jhe4)*state % rho - screened_rates(k_he4_fe46__p_co49)* &
      Y(jfe46)*Y(jhe4)*state % rho + screened_rates(k_ni50__he4_fe46)* &
      Y(jni50) + screened_rates(k_p_co49__he4_fe46)*Y(jco49)*Y(jp)* &
      state % rho - screened_rates(k_p_fe46__co47)*Y(jfe46)*Y(jp)*state % rho + &
      screened_rates(k_p_mn45__fe46)*Y(jmn45)*Y(jp)*state % rho &
       )

    ydot_nuc(jfe47) = ( &
      screened_rates(k_co47__fe47__weak__bqa_pos_)*Y(jco47) + screened_rates(k_co48__p_fe47)* &
      Y(jco48) - screened_rates(k_fe47__he4_cr43)*Y(jfe47) - &
      screened_rates(k_fe47__mn47__weak__wc12)*Y(jfe47) - &
      screened_rates(k_fe47__p_cr46__weak__wc12)*Y(jfe47) - &
      screened_rates(k_fe47__p_mn46)*Y(jfe47) + screened_rates(k_he4_cr43__fe47)* &
      Y(jcr43)*Y(jhe4)*state % rho - screened_rates(k_he4_fe47__ni51)* &
      Y(jfe47)*Y(jhe4)*state % rho - screened_rates(k_he4_fe47__p_co50)* &
      Y(jfe47)*Y(jhe4)*state % rho + screened_rates(k_he4_mn44__p_fe47)* &
      Y(jhe4)*Y(jmn44)*state % rho + screened_rates(k_ni51__he4_fe47)* &
      Y(jni51) + screened_rates(k_p_co50__he4_fe47)*Y(jco50)*Y(jp)* &
      state % rho - screened_rates(k_p_fe47__co48)*Y(jfe47)*Y(jp)*state % rho - &
      screened_rates(k_p_fe47__he4_mn44)*Y(jfe47)*Y(jp)*state % rho + &
      screened_rates(k_p_mn46__fe47)*Y(jmn46)*Y(jp)*state % rho &
       )

    ydot_nuc(jfe48) = ( &
      screened_rates(k_co48__fe48__weak__bqa_pos_)*Y(jco48) + screened_rates(k_co49__p_fe48)* &
      Y(jco49) - screened_rates(k_fe48__he4_cr44)*Y(jfe48) - &
      screened_rates(k_fe48__mn48__weak__wc12)*Y(jfe48) - &
      screened_rates(k_fe48__p_cr47__weak__wc12)*Y(jfe48) - &
      screened_rates(k_fe48__p_mn47)*Y(jfe48) + screened_rates(k_he4_cr44__fe48)* &
      Y(jcr44)*Y(jhe4)*state % rho - screened_rates(k_he4_fe48__ni52)* &
      Y(jfe48)*Y(jhe4)*state % rho - screened_rates(k_he4_fe48__p_co51)* &
      Y(jfe48)*Y(jhe4)*state % rho + screened_rates(k_he4_mn45__p_fe48)* &
      Y(jhe4)*Y(jmn45)*state % rho + screened_rates(k_ni49__p_fe48__weak__wc12) &
      *Y(jni49) + screened_rates(k_ni52__he4_fe48)*Y(jni52) + &
      screened_rates(k_p_co51__he4_fe48)*Y(jco51)*Y(jp)*state % rho - &
      screened_rates(k_p_fe48__co49)*Y(jfe48)*Y(jp)*state % rho - &
      screened_rates(k_p_fe48__he4_mn45)*Y(jfe48)*Y(jp)*state % rho + &
      screened_rates(k_p_mn47__fe48)*Y(jmn47)*Y(jp)*state % rho &
       )

    ydot_nuc(jfe49) = ( &
      screened_rates(k_co49__fe49__weak__bqa_pos_)*Y(jco49) + screened_rates(k_co50__p_fe49)* &
      Y(jco50) - screened_rates(k_fe49__he4_cr45)*Y(jfe49) - &
      screened_rates(k_fe49__mn49__weak__wc12)*Y(jfe49) - &
      screened_rates(k_fe49__p_cr48__weak__wc12)*Y(jfe49) - &
      screened_rates(k_fe49__p_mn48)*Y(jfe49) + screened_rates(k_he4_cr45__fe49)* &
      Y(jcr45)*Y(jhe4)*state % rho - screened_rates(k_he4_fe49__ni53)* &
      Y(jfe49)*Y(jhe4)*state % rho - screened_rates(k_he4_fe49__p_co52)* &
      Y(jfe49)*Y(jhe4)*state % rho + screened_rates(k_he4_mn46__p_fe49)* &
      Y(jhe4)*Y(jmn46)*state % rho + screened_rates(k_ni50__p_fe49__weak__wc12) &
      *Y(jni50) + screened_rates(k_ni53__he4_fe49)*Y(jni53) + &
      screened_rates(k_p_co52__he4_fe49)*Y(jco52)*Y(jp)*state % rho - &
      screened_rates(k_p_fe49__co50)*Y(jfe49)*Y(jp)*state % rho - &
      screened_rates(k_p_fe49__he4_mn46)*Y(jfe49)*Y(jp)*state % rho + &
      screened_rates(k_p_mn48__fe49)*Y(jmn48)*Y(jp)*state % rho &
       )

    ydot_nuc(jfe50) = ( &
      screened_rates(k_co50__fe50__weak__wc12)*Y(jco50) + screened_rates(k_co51__p_fe50)* &
      Y(jco51) - screened_rates(k_fe50__he4_cr46)*Y(jfe50) - &
      screened_rates(k_fe50__mn50__weak__wc12)*Y(jfe50) - &
      screened_rates(k_fe50__p_mn49)*Y(jfe50) + screened_rates(k_he4_cr46__fe50)* &
      Y(jcr46)*Y(jhe4)*state % rho - screened_rates(k_he4_fe50__ni54)* &
      Y(jfe50)*Y(jhe4)*state % rho - screened_rates(k_he4_fe50__p_co53)* &
      Y(jfe50)*Y(jhe4)*state % rho + screened_rates(k_he4_mn47__p_fe50)* &
      Y(jhe4)*Y(jmn47)*state % rho + screened_rates(k_ni51__p_fe50__weak__wc12) &
      *Y(jni51) + screened_rates(k_ni54__he4_fe50)*Y(jni54) + &
      screened_rates(k_p_co53__he4_fe50)*Y(jco53)*Y(jp)*state % rho - &
      screened_rates(k_p_fe50__co51)*Y(jfe50)*Y(jp)*state % rho - &
      screened_rates(k_p_fe50__he4_mn47)*Y(jfe50)*Y(jp)*state % rho + &
      screened_rates(k_p_mn49__fe50)*Y(jmn49)*Y(jp)*state % rho &
       )

    ydot_nuc(jfe51) = ( &
      screened_rates(k_co51__fe51__weak__mo97)*Y(jco51) + screened_rates(k_co52__p_fe51)* &
      Y(jco52) - screened_rates(k_fe51__he4_cr47)*Y(jfe51) - &
      screened_rates(k_fe51__mn51__weak__wc12)*Y(jfe51) - &
      screened_rates(k_fe51__p_mn50)*Y(jfe51) + screened_rates(k_he4_cr47__fe51)* &
      Y(jcr47)*Y(jhe4)*state % rho - screened_rates(k_he4_fe51__ni55)* &
      Y(jfe51)*Y(jhe4)*state % rho - screened_rates(k_he4_fe51__p_co54)* &
      Y(jfe51)*Y(jhe4)*state % rho + screened_rates(k_he4_mn48__p_fe51)* &
      Y(jhe4)*Y(jmn48)*state % rho + screened_rates(k_ni52__p_fe51__weak__wc12) &
      *Y(jni52) + screened_rates(k_ni55__he4_fe51)*Y(jni55) + &
      screened_rates(k_p_co54__he4_fe51)*Y(jco54)*Y(jp)*state % rho - &
      screened_rates(k_p_fe51__co52)*Y(jfe51)*Y(jp)*state % rho - &
      screened_rates(k_p_fe51__he4_mn48)*Y(jfe51)*Y(jp)*state % rho + &
      screened_rates(k_p_mn50__fe51)*Y(jmn50)*Y(jp)*state % rho &
       )

    ydot_nuc(jfe52) = ( &
      screened_rates(k_co52__fe52__weak__wc12)*Y(jco52) + screened_rates(k_co53__p_fe52)* &
      Y(jco53) - screened_rates(k_fe52__he4_cr48)*Y(jfe52) - &
      screened_rates(k_fe52__mn52__weak__wc12)*Y(jfe52) - &
      screened_rates(k_fe52__p_mn51)*Y(jfe52) + screened_rates(k_he4_cr48__fe52)* &
      Y(jcr48)*Y(jhe4)*state % rho - screened_rates(k_he4_fe52__ni56)* &
      Y(jfe52)*Y(jhe4)*state % rho - screened_rates(k_he4_fe52__p_co55)* &
      Y(jfe52)*Y(jhe4)*state % rho + screened_rates(k_he4_mn49__p_fe52)* &
      Y(jhe4)*Y(jmn49)*state % rho + screened_rates(k_ni53__p_fe52__weak__wc12) &
      *Y(jni53) + screened_rates(k_ni56__he4_fe52)*Y(jni56) + &
      screened_rates(k_p_co55__he4_fe52)*Y(jco55)*Y(jp)*state % rho - &
      screened_rates(k_p_fe52__co53)*Y(jfe52)*Y(jp)*state % rho - &
      screened_rates(k_p_fe52__he4_mn49)*Y(jfe52)*Y(jp)*state % rho + &
      screened_rates(k_p_mn51__fe52)*Y(jmn51)*Y(jp)*state % rho &
       )

    ydot_nuc(jfe53) = ( &
      screened_rates(k_co53__fe53__weak__wc12)*Y(jco53) + screened_rates(k_co54__p_fe53)* &
      Y(jco54) - screened_rates(k_fe53__he4_cr49)*Y(jfe53) - &
      screened_rates(k_fe53__mn53__weak__wc12)*Y(jfe53) - &
      screened_rates(k_fe53__p_mn52)*Y(jfe53) + screened_rates(k_he4_cr49__fe53)* &
      Y(jcr49)*Y(jhe4)*state % rho - screened_rates(k_he4_fe53__p_co56)* &
      Y(jfe53)*Y(jhe4)*state % rho + screened_rates(k_he4_mn50__p_fe53)* &
      Y(jhe4)*Y(jmn50)*state % rho + screened_rates(k_p_co56__he4_fe53)* &
      Y(jco56)*Y(jp)*state % rho - screened_rates(k_p_fe53__co54)*Y(jfe53)* &
      Y(jp)*state % rho - screened_rates(k_p_fe53__he4_mn50)*Y(jfe53)*Y(jp)* &
      state % rho + screened_rates(k_p_mn52__fe53)*Y(jmn52)*Y(jp)*state % rho &
       )

    ydot_nuc(jfe54) = ( &
      screened_rates(k_co54__fe54__weak__wc12)*Y(jco54) + screened_rates(k_co55__p_fe54)* &
      Y(jco55) - screened_rates(k_fe54__he4_cr50)*Y(jfe54) - &
      screened_rates(k_fe54__p_mn53)*Y(jfe54) + screened_rates(k_he4_cr50__fe54)* &
      Y(jcr50)*Y(jhe4)*state % rho + screened_rates(k_he4_mn51__p_fe54)* &
      Y(jhe4)*Y(jmn51)*state % rho - screened_rates(k_p_fe54__co55)*Y(jfe54)* &
      Y(jp)*state % rho - screened_rates(k_p_fe54__he4_mn51)*Y(jfe54)*Y(jp)* &
      state % rho + screened_rates(k_p_mn53__fe54)*Y(jmn53)*Y(jp)*state % rho &
       )

    ydot_nuc(jfe55) = ( &
      screened_rates(k_co55__fe55__weak__wc12)*Y(jco55) + screened_rates(k_co56__p_fe55)* &
      Y(jco56) - screened_rates(k_fe55__he4_cr51)*Y(jfe55) - &
      screened_rates(k_fe55__mn55__weak__wc12)*Y(jfe55) - &
      screened_rates(k_fe55__p_mn54)*Y(jfe55) + screened_rates(k_he4_cr51__fe55)* &
      Y(jcr51)*Y(jhe4)*state % rho + screened_rates(k_he4_mn52__p_fe55)* &
      Y(jhe4)*Y(jmn52)*state % rho - screened_rates(k_p_fe55__co56)*Y(jfe55)* &
      Y(jp)*state % rho - screened_rates(k_p_fe55__he4_mn52)*Y(jfe55)*Y(jp)* &
      state % rho + screened_rates(k_p_mn54__fe55)*Y(jmn54)*Y(jp)*state % rho &
       )

    ydot_nuc(jfe56) = ( &
      screened_rates(k_co56__fe56__weak__wc12)*Y(jco56) - screened_rates(k_fe56__he4_cr52)* &
      Y(jfe56) - screened_rates(k_fe56__p_mn55)*Y(jfe56) + &
      screened_rates(k_he4_cr52__fe56)*Y(jcr52)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mn53__p_fe56)*Y(jhe4)*Y(jmn53)*state % rho - &
      screened_rates(k_p_fe56__he4_mn53)*Y(jfe56)*Y(jp)*state % rho + &
      screened_rates(k_p_mn55__fe56)*Y(jmn55)*Y(jp)*state % rho &
       )

    ydot_nuc(jco47) = ( &
      -screened_rates(k_co47__fe47__weak__bqa_pos_)*Y(jco47) - screened_rates(k_co47__p_fe46)* &
      Y(jco47) - screened_rates(k_he4_co47__p_ni50)*Y(jco47)*Y(jhe4)* &
      state % rho + screened_rates(k_ni48__p_co47)*Y(jni48) - &
      screened_rates(k_p_co47__ni48)*Y(jco47)*Y(jp)*state % rho + &
      screened_rates(k_p_fe46__co47)*Y(jfe46)*Y(jp)*state % rho + &
      screened_rates(k_p_ni50__he4_co47)*Y(jni50)*Y(jp)*state % rho &
       )

    ydot_nuc(jco48) = ( &
      -screened_rates(k_co48__fe48__weak__bqa_pos_)*Y(jco48) - &
      screened_rates(k_co48__he4_mn44)*Y(jco48) - screened_rates(k_co48__p_fe47)* &
      Y(jco48) - screened_rates(k_he4_co48__p_ni51)*Y(jco48)*Y(jhe4)* &
      state % rho + screened_rates(k_he4_fe45__p_co48)*Y(jfe45)*Y(jhe4)* &
      state % rho + screened_rates(k_he4_mn44__co48)*Y(jhe4)*Y(jmn44)*state % rho &
      + screened_rates(k_ni48__co48__weak__wc17)*Y(jni48) + &
      screened_rates(k_ni49__p_co48)*Y(jni49) - screened_rates(k_p_co48__he4_fe45)* &
      Y(jco48)*Y(jp)*state % rho - screened_rates(k_p_co48__ni49)*Y(jco48)* &
      Y(jp)*state % rho + screened_rates(k_p_fe47__co48)*Y(jfe47)*Y(jp)* &
      state % rho + screened_rates(k_p_ni51__he4_co48)*Y(jni51)*Y(jp)*state % rho &
       )

    ydot_nuc(jco49) = ( &
      -screened_rates(k_co49__fe49__weak__bqa_pos_)*Y(jco49) - &
      screened_rates(k_co49__he4_mn45)*Y(jco49) - screened_rates(k_co49__p_fe48)* &
      Y(jco49) - screened_rates(k_he4_co49__p_ni52)*Y(jco49)*Y(jhe4)* &
      state % rho + screened_rates(k_he4_fe46__p_co49)*Y(jfe46)*Y(jhe4)* &
      state % rho + screened_rates(k_he4_mn45__co49)*Y(jhe4)*Y(jmn45)*state % rho &
      + screened_rates(k_ni49__co49__weak__wc12)*Y(jni49) + &
      screened_rates(k_ni50__p_co49)*Y(jni50) - screened_rates(k_p_co49__he4_fe46)* &
      Y(jco49)*Y(jp)*state % rho - screened_rates(k_p_co49__ni50)*Y(jco49)* &
      Y(jp)*state % rho + screened_rates(k_p_fe48__co49)*Y(jfe48)*Y(jp)* &
      state % rho + screened_rates(k_p_ni52__he4_co49)*Y(jni52)*Y(jp)*state % rho &
       )

    ydot_nuc(jco50) = ( &
      -screened_rates(k_co50__fe50__weak__wc12)*Y(jco50) - screened_rates(k_co50__he4_mn46)* &
      Y(jco50) - screened_rates(k_co50__p_fe49)*Y(jco50) - &
      screened_rates(k_co50__p_mn49__weak__wc12)*Y(jco50) - &
      screened_rates(k_he4_co50__p_ni53)*Y(jco50)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe47__p_co50)*Y(jfe47)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mn46__co50)*Y(jhe4)*Y(jmn46)*state % rho + &
      screened_rates(k_ni50__co50__weak__wc12)*Y(jni50) + &
      screened_rates(k_ni51__p_co50)*Y(jni51) - screened_rates(k_p_co50__he4_fe47)* &
      Y(jco50)*Y(jp)*state % rho - screened_rates(k_p_co50__ni51)*Y(jco50)* &
      Y(jp)*state % rho + screened_rates(k_p_fe49__co50)*Y(jfe49)*Y(jp)* &
      state % rho + screened_rates(k_p_ni53__he4_co50)*Y(jni53)*Y(jp)*state % rho &
       )

    ydot_nuc(jco51) = ( &
      -screened_rates(k_co51__fe51__weak__mo97)*Y(jco51) - screened_rates(k_co51__he4_mn47)* &
      Y(jco51) - screened_rates(k_co51__p_fe50)*Y(jco51) - &
      screened_rates(k_he4_co51__p_ni54)*Y(jco51)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe48__p_co51)*Y(jfe48)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mn47__co51)*Y(jhe4)*Y(jmn47)*state % rho + &
      screened_rates(k_ni51__co51__weak__wc17)*Y(jni51) + &
      screened_rates(k_ni52__p_co51)*Y(jni52) - screened_rates(k_p_co51__he4_fe48)* &
      Y(jco51)*Y(jp)*state % rho - screened_rates(k_p_co51__ni52)*Y(jco51)* &
      Y(jp)*state % rho + screened_rates(k_p_fe50__co51)*Y(jfe50)*Y(jp)* &
      state % rho + screened_rates(k_p_ni54__he4_co51)*Y(jni54)*Y(jp)*state % rho &
       )

    ydot_nuc(jco52) = ( &
      -screened_rates(k_co52__fe52__weak__wc12)*Y(jco52) - screened_rates(k_co52__he4_mn48)* &
      Y(jco52) - screened_rates(k_co52__p_fe51)*Y(jco52) - &
      screened_rates(k_he4_co52__p_ni55)*Y(jco52)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe49__p_co52)*Y(jfe49)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mn48__co52)*Y(jhe4)*Y(jmn48)*state % rho + &
      screened_rates(k_ni52__co52__weak__wc12)*Y(jni52) + &
      screened_rates(k_ni53__p_co52)*Y(jni53) - screened_rates(k_p_co52__he4_fe49)* &
      Y(jco52)*Y(jp)*state % rho - screened_rates(k_p_co52__ni53)*Y(jco52)* &
      Y(jp)*state % rho + screened_rates(k_p_fe51__co52)*Y(jfe51)*Y(jp)* &
      state % rho + screened_rates(k_p_ni55__he4_co52)*Y(jni55)*Y(jp)*state % rho &
       )

    ydot_nuc(jco53) = ( &
      -screened_rates(k_co53__fe53__weak__wc12)*Y(jco53) - screened_rates(k_co53__he4_mn49)* &
      Y(jco53) - screened_rates(k_co53__p_fe52)*Y(jco53) - &
      screened_rates(k_he4_co53__p_ni56)*Y(jco53)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe50__p_co53)*Y(jfe50)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mn49__co53)*Y(jhe4)*Y(jmn49)*state % rho + &
      screened_rates(k_ni53__co53__weak__wc12)*Y(jni53) + &
      screened_rates(k_ni54__p_co53)*Y(jni54) - screened_rates(k_p_co53__he4_fe50)* &
      Y(jco53)*Y(jp)*state % rho - screened_rates(k_p_co53__ni54)*Y(jco53)* &
      Y(jp)*state % rho + screened_rates(k_p_fe52__co53)*Y(jfe52)*Y(jp)* &
      state % rho + screened_rates(k_p_ni56__he4_co53)*Y(jni56)*Y(jp)*state % rho &
       )

    ydot_nuc(jco54) = ( &
      -screened_rates(k_co54__fe54__weak__wc12)*Y(jco54) - screened_rates(k_co54__he4_mn50)* &
      Y(jco54) - screened_rates(k_co54__p_fe53)*Y(jco54) + &
      screened_rates(k_he4_fe51__p_co54)*Y(jfe51)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mn50__co54)*Y(jhe4)*Y(jmn50)*state % rho + &
      screened_rates(k_ni54__co54__weak__wc12)*Y(jni54) + &
      screened_rates(k_ni55__p_co54)*Y(jni55) - screened_rates(k_p_co54__he4_fe51)* &
      Y(jco54)*Y(jp)*state % rho - screened_rates(k_p_co54__ni55)*Y(jco54)* &
      Y(jp)*state % rho + screened_rates(k_p_fe53__co54)*Y(jfe53)*Y(jp)* &
      state % rho &
       )

    ydot_nuc(jco55) = ( &
      -screened_rates(k_co55__fe55__weak__wc12)*Y(jco55) - screened_rates(k_co55__he4_mn51)* &
      Y(jco55) - screened_rates(k_co55__p_fe54)*Y(jco55) + &
      screened_rates(k_he4_fe52__p_co55)*Y(jfe52)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mn51__co55)*Y(jhe4)*Y(jmn51)*state % rho + &
      screened_rates(k_ni55__co55__weak__wc12)*Y(jni55) + &
      screened_rates(k_ni56__p_co55)*Y(jni56) - screened_rates(k_p_co55__he4_fe52)* &
      Y(jco55)*Y(jp)*state % rho - screened_rates(k_p_co55__ni56)*Y(jco55)* &
      Y(jp)*state % rho + screened_rates(k_p_fe54__co55)*Y(jfe54)*Y(jp)* &
      state % rho &
       )

    ydot_nuc(jco56) = ( &
      -screened_rates(k_co56__fe56__weak__wc12)*Y(jco56) - screened_rates(k_co56__he4_mn52)* &
      Y(jco56) - screened_rates(k_co56__p_fe55)*Y(jco56) + &
      screened_rates(k_he4_fe53__p_co56)*Y(jfe53)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mn52__co56)*Y(jhe4)*Y(jmn52)*state % rho + &
      screened_rates(k_ni56__co56__weak__wc12)*Y(jni56) - &
      screened_rates(k_p_co56__he4_fe53)*Y(jco56)*Y(jp)*state % rho + &
      screened_rates(k_p_fe55__co56)*Y(jfe55)*Y(jp)*state % rho &
       )

    ydot_nuc(jni48) = ( &
      -screened_rates(k_ni48__co48__weak__wc17)*Y(jni48) - screened_rates(k_ni48__p_co47)* &
      Y(jni48) + screened_rates(k_p_co47__ni48)*Y(jco47)*Y(jp)*state % rho &
       )

    ydot_nuc(jni49) = ( &
      screened_rates(k_he4_fe45__ni49)*Y(jfe45)*Y(jhe4)*state % rho - &
      screened_rates(k_ni49__co49__weak__wc12)*Y(jni49) - &
      screened_rates(k_ni49__he4_fe45)*Y(jni49) - screened_rates(k_ni49__p_co48)* &
      Y(jni49) - screened_rates(k_ni49__p_fe48__weak__wc12)*Y(jni49) + &
      screened_rates(k_p_co48__ni49)*Y(jco48)*Y(jp)*state % rho &
       )

    ydot_nuc(jni50) = ( &
      screened_rates(k_he4_co47__p_ni50)*Y(jco47)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe46__ni50)*Y(jfe46)*Y(jhe4)*state % rho - &
      screened_rates(k_ni50__co50__weak__wc12)*Y(jni50) - &
      screened_rates(k_ni50__he4_fe46)*Y(jni50) - screened_rates(k_ni50__p_co49)* &
      Y(jni50) - screened_rates(k_ni50__p_fe49__weak__wc12)*Y(jni50) + &
      screened_rates(k_p_co49__ni50)*Y(jco49)*Y(jp)*state % rho - &
      screened_rates(k_p_ni50__he4_co47)*Y(jni50)*Y(jp)*state % rho &
       )

    ydot_nuc(jni51) = ( &
      screened_rates(k_he4_co48__p_ni51)*Y(jco48)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe47__ni51)*Y(jfe47)*Y(jhe4)*state % rho - &
      screened_rates(k_ni51__co51__weak__wc17)*Y(jni51) - &
      screened_rates(k_ni51__he4_fe47)*Y(jni51) - screened_rates(k_ni51__p_co50)* &
      Y(jni51) - screened_rates(k_ni51__p_fe50__weak__wc12)*Y(jni51) + &
      screened_rates(k_p_co50__ni51)*Y(jco50)*Y(jp)*state % rho - &
      screened_rates(k_p_ni51__he4_co48)*Y(jni51)*Y(jp)*state % rho &
       )

    ydot_nuc(jni52) = ( &
      screened_rates(k_he4_co49__p_ni52)*Y(jco49)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe48__ni52)*Y(jfe48)*Y(jhe4)*state % rho - &
      screened_rates(k_ni52__co52__weak__wc12)*Y(jni52) - &
      screened_rates(k_ni52__he4_fe48)*Y(jni52) - screened_rates(k_ni52__p_co51)* &
      Y(jni52) - screened_rates(k_ni52__p_fe51__weak__wc12)*Y(jni52) + &
      screened_rates(k_p_co51__ni52)*Y(jco51)*Y(jp)*state % rho - &
      screened_rates(k_p_ni52__he4_co49)*Y(jni52)*Y(jp)*state % rho &
       )

    ydot_nuc(jni53) = ( &
      screened_rates(k_he4_co50__p_ni53)*Y(jco50)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe49__ni53)*Y(jfe49)*Y(jhe4)*state % rho - &
      screened_rates(k_ni53__co53__weak__wc12)*Y(jni53) - &
      screened_rates(k_ni53__he4_fe49)*Y(jni53) - screened_rates(k_ni53__p_co52)* &
      Y(jni53) - screened_rates(k_ni53__p_fe52__weak__wc12)*Y(jni53) + &
      screened_rates(k_p_co52__ni53)*Y(jco52)*Y(jp)*state % rho - &
      screened_rates(k_p_ni53__he4_co50)*Y(jni53)*Y(jp)*state % rho &
       )

    ydot_nuc(jni54) = ( &
      screened_rates(k_he4_co51__p_ni54)*Y(jco51)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe50__ni54)*Y(jfe50)*Y(jhe4)*state % rho - &
      screened_rates(k_ni54__co54__weak__wc12)*Y(jni54) - &
      screened_rates(k_ni54__he4_fe50)*Y(jni54) - screened_rates(k_ni54__p_co53)* &
      Y(jni54) + screened_rates(k_p_co53__ni54)*Y(jco53)*Y(jp)*state % rho - &
      screened_rates(k_p_ni54__he4_co51)*Y(jni54)*Y(jp)*state % rho &
       )

    ydot_nuc(jni55) = ( &
      screened_rates(k_he4_co52__p_ni55)*Y(jco52)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe51__ni55)*Y(jfe51)*Y(jhe4)*state % rho - &
      screened_rates(k_ni55__co55__weak__wc12)*Y(jni55) - &
      screened_rates(k_ni55__he4_fe51)*Y(jni55) - screened_rates(k_ni55__p_co54)* &
      Y(jni55) + screened_rates(k_p_co54__ni55)*Y(jco54)*Y(jp)*state % rho - &
      screened_rates(k_p_ni55__he4_co52)*Y(jni55)*Y(jp)*state % rho &
       )

    ydot_nuc(jni56) = ( &
      screened_rates(k_he4_co53__p_ni56)*Y(jco53)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_fe52__ni56)*Y(jfe52)*Y(jhe4)*state % rho - &
      screened_rates(k_ni56__co56__weak__wc12)*Y(jni56) - &
      screened_rates(k_ni56__he4_fe52)*Y(jni56) - screened_rates(k_ni56__p_co55)* &
      Y(jni56) + screened_rates(k_p_co55__ni56)*Y(jco55)*Y(jp)*state % rho - &
      screened_rates(k_p_ni56__he4_co53)*Y(jni56)*Y(jp)*state % rho &
       )


  end subroutine rhs_nuc


  subroutine actual_jac(state)

    !$acc routine seq

    use burn_type_module, only: net_itemp, net_ienuc
    use sneut_module, only: sneut5
    use temperature_integration_module, only: temperature_jac
    use jacobian_sparsity_module, only: get_jac_entry, set_jac_entry, set_jac_zero

    implicit none
    
    type(burn_t) :: state
    type(rate_eval_t) :: rate_eval
    real(rt) :: reactvec(num_rate_groups+2)
    real(rt) :: screened_rates_dt(nrates)
    real(rt) :: Y(nspec), yderivs(nspec)
    real(rt) :: ye, rhoy, b1, scratch
    real(rt) :: sneut, dsneutdt, dsneutdd, snuda, snudz
    integer  :: j, k

    !$gpu

    ! Set molar abundances
    Y(:) = state % xn(:) * aion_inv(:)
    
    call evaluate_rates(state, rate_eval)

    ! Zero out the Jacobian
    call set_jac_zero(state)

    ! Species Jacobian elements with respect to other species
    call jac_nuc(state, Y, rate_eval % screened_rates)

    ! Evaluate the species Jacobian elements with respect to temperature by
    ! calling the RHS using the temperature derivative of the screened rate
    screened_rates_dt = rate_eval % unscreened_rates(i_rate, :) * &
                        rate_eval % unscreened_rates(i_dscor_dt, :) + &
                        rate_eval % unscreened_rates(i_drate_dt, :) * &
                        rate_eval % unscreened_rates(i_scor, :)

    call rhs_nuc(state, yderivs, Y, screened_rates_dt)

    do k = 1, nspec
       call set_jac_entry(state, k, net_itemp, yderivs(k))
    enddo

    ! Energy generation rate Jacobian elements with respect to species
    do j = 1, nspec
       do k = 1, nspec
          call get_jac_entry(state, k, j, yderivs(k))
       enddo
       call ener_gener_rate(yderivs, scratch)
       call set_jac_entry(state, net_ienuc, j, scratch)
    enddo

    ! Account for the thermal neutrino losses
    call sneut5(state % T, state % rho, state % abar, state % zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)

    do j = 1, nspec
       b1 = ((aion(j) - state % abar) * state % abar * snuda + (zion(j) - state % zbar) * state % abar * snudz)
       call get_jac_entry(state, net_ienuc, j, scratch)
       scratch = scratch - b1
       call set_jac_entry(state, net_ienuc, j, scratch)
    enddo

    ! Energy generation rate Jacobian element with respect to temperature
    do k = 1, nspec
       call get_jac_entry(state, k, net_itemp, yderivs(k))
    enddo
    call ener_gener_rate(yderivs, scratch)
    scratch = scratch - dsneutdt    
    call set_jac_entry(state, net_ienuc, net_itemp, scratch)

    ! Temperature Jacobian elements
    call temperature_jac(state)

  end subroutine actual_jac


  subroutine jac_nuc(state, Y, screened_rates)

    !$acc routine seq

    use jacobian_sparsity_module, only: set_jac_entry

    implicit none

    type(burn_t), intent(inout) :: state
    real(rt), intent(in)  :: Y(nspec)
    real(rt), intent(in)  :: screened_rates(nrates)
    real(rt) :: scratch


    !$gpu


    scratch = (&
      -screened_rates(k_p_al23__si24)*Y(jal23)*state % rho - screened_rates(k_p_al24__he4_mg21)* &
      Y(jal24)*state % rho - screened_rates(k_p_al24__si25)*Y(jal24)*state % rho - &
      screened_rates(k_p_al25__he4_mg22)*Y(jal25)*state % rho - &
      screened_rates(k_p_al25__si26)*Y(jal25)*state % rho - &
      screened_rates(k_p_al26__he4_mg23)*Y(jal26)*state % rho - &
      screened_rates(k_p_al26__si27)*Y(jal26)*state % rho - &
      screened_rates(k_p_al27__he4_mg24)*Y(jal27)*state % rho - &
      screened_rates(k_p_al27__si28)*Y(jal27)*state % rho - &
      screened_rates(k_p_al28__he4_mg25)*Y(jal28)*state % rho - &
      screened_rates(k_p_al28__si29)*Y(jal28)*state % rho - &
      screened_rates(k_p_ar32__he4_cl29)*Y(jar32)*state % rho - &
      screened_rates(k_p_ar32__k33)*Y(jar32)*state % rho - &
      screened_rates(k_p_ar33__he4_cl30)*Y(jar33)*state % rho - &
      screened_rates(k_p_ar33__k34)*Y(jar33)*state % rho - &
      screened_rates(k_p_ar34__he4_cl31)*Y(jar34)*state % rho - &
      screened_rates(k_p_ar34__k35)*Y(jar34)*state % rho - &
      screened_rates(k_p_ar35__he4_cl32)*Y(jar35)*state % rho - &
      screened_rates(k_p_ar35__k36)*Y(jar35)*state % rho - &
      screened_rates(k_p_ar36__he4_cl33)*Y(jar36)*state % rho - &
      screened_rates(k_p_ar36__k37)*Y(jar36)*state % rho - &
      screened_rates(k_p_ar37__he4_cl34)*Y(jar37)*state % rho - &
      screened_rates(k_p_ar37__k38)*Y(jar37)*state % rho - &
      screened_rates(k_p_ar38__he4_cl35)*Y(jar38)*state % rho - &
      screened_rates(k_p_ar38__k39)*Y(jar38)*state % rho - &
      screened_rates(k_p_ar39__he4_cl36)*Y(jar39)*state % rho - &
      screened_rates(k_p_ar39__k40)*Y(jar39)*state % rho - screened_rates(k_p_be7__b8)* &
      Y(jbe7)*state % rho - screened_rates(k_p_c12__n13)*Y(jc12)*state % rho - &
      screened_rates(k_p_c13__n14)*Y(jc13)*state % rho - screened_rates(k_p_ca35__sc36)* &
      Y(jca35)*state % rho - screened_rates(k_p_ca36__he4_k33)*Y(jca36)* &
      state % rho - screened_rates(k_p_ca36__sc37)*Y(jca36)*state % rho - &
      screened_rates(k_p_ca37__he4_k34)*Y(jca37)*state % rho - &
      screened_rates(k_p_ca37__sc38)*Y(jca37)*state % rho - &
      screened_rates(k_p_ca38__he4_k35)*Y(jca38)*state % rho - &
      screened_rates(k_p_ca38__sc39)*Y(jca38)*state % rho - &
      screened_rates(k_p_ca39__he4_k36)*Y(jca39)*state % rho - &
      screened_rates(k_p_ca39__sc40)*Y(jca39)*state % rho - &
      screened_rates(k_p_ca40__he4_k37)*Y(jca40)*state % rho - &
      screened_rates(k_p_ca40__sc41)*Y(jca40)*state % rho - &
      screened_rates(k_p_ca41__he4_k38)*Y(jca41)*state % rho - &
      screened_rates(k_p_ca41__sc42)*Y(jca41)*state % rho - &
      screened_rates(k_p_ca42__he4_k39)*Y(jca42)*state % rho - &
      screened_rates(k_p_ca42__sc43)*Y(jca42)*state % rho - &
      screened_rates(k_p_ca43__he4_k40)*Y(jca43)*state % rho - &
      screened_rates(k_p_ca43__sc44)*Y(jca43)*state % rho - &
      screened_rates(k_p_ca44__he4_k41)*Y(jca44)*state % rho - &
      screened_rates(k_p_ca44__sc45)*Y(jca44)*state % rho - screened_rates(k_p_cl30__ar31) &
      *Y(jcl30)*state % rho - screened_rates(k_p_cl31__ar32)*Y(jcl31)*state % rho &
      - screened_rates(k_p_cl31__he4_s28)*Y(jcl31)*state % rho - &
      screened_rates(k_p_cl32__ar33)*Y(jcl32)*state % rho - &
      screened_rates(k_p_cl32__he4_s29)*Y(jcl32)*state % rho - &
      screened_rates(k_p_cl33__ar34)*Y(jcl33)*state % rho - &
      screened_rates(k_p_cl33__he4_s30)*Y(jcl33)*state % rho - &
      screened_rates(k_p_cl34__ar35)*Y(jcl34)*state % rho - &
      screened_rates(k_p_cl34__he4_s31)*Y(jcl34)*state % rho - &
      screened_rates(k_p_cl35__ar36)*Y(jcl35)*state % rho - &
      screened_rates(k_p_cl35__he4_s32)*Y(jcl35)*state % rho - &
      screened_rates(k_p_cl36__ar37)*Y(jcl36)*state % rho - &
      screened_rates(k_p_cl36__he4_s33)*Y(jcl36)*state % rho - &
      screened_rates(k_p_cl37__ar38)*Y(jcl37)*state % rho - &
      screened_rates(k_p_cl37__he4_s34)*Y(jcl37)*state % rho - &
      screened_rates(k_p_co47__ni48)*Y(jco47)*state % rho - &
      screened_rates(k_p_co48__he4_fe45)*Y(jco48)*state % rho - &
      screened_rates(k_p_co48__ni49)*Y(jco48)*state % rho - &
      screened_rates(k_p_co49__he4_fe46)*Y(jco49)*state % rho - &
      screened_rates(k_p_co49__ni50)*Y(jco49)*state % rho - &
      screened_rates(k_p_co50__he4_fe47)*Y(jco50)*state % rho - &
      screened_rates(k_p_co50__ni51)*Y(jco50)*state % rho - &
      screened_rates(k_p_co51__he4_fe48)*Y(jco51)*state % rho - &
      screened_rates(k_p_co51__ni52)*Y(jco51)*state % rho - &
      screened_rates(k_p_co52__he4_fe49)*Y(jco52)*state % rho - &
      screened_rates(k_p_co52__ni53)*Y(jco52)*state % rho - &
      screened_rates(k_p_co53__he4_fe50)*Y(jco53)*state % rho - &
      screened_rates(k_p_co53__ni54)*Y(jco53)*state % rho - &
      screened_rates(k_p_co54__he4_fe51)*Y(jco54)*state % rho - &
      screened_rates(k_p_co54__ni55)*Y(jco54)*state % rho - &
      screened_rates(k_p_co55__he4_fe52)*Y(jco55)*state % rho - &
      screened_rates(k_p_co55__ni56)*Y(jco55)*state % rho - &
      screened_rates(k_p_co56__he4_fe53)*Y(jco56)*state % rho - &
      screened_rates(k_p_cr43__he4_v40)*Y(jcr43)*state % rho - &
      screened_rates(k_p_cr43__mn44)*Y(jcr43)*state % rho - &
      screened_rates(k_p_cr44__he4_v41)*Y(jcr44)*state % rho - &
      screened_rates(k_p_cr44__mn45)*Y(jcr44)*state % rho - &
      screened_rates(k_p_cr45__he4_v42)*Y(jcr45)*state % rho - &
      screened_rates(k_p_cr45__mn46)*Y(jcr45)*state % rho - &
      screened_rates(k_p_cr46__he4_v43)*Y(jcr46)*state % rho - &
      screened_rates(k_p_cr46__mn47)*Y(jcr46)*state % rho - &
      screened_rates(k_p_cr47__he4_v44)*Y(jcr47)*state % rho - &
      screened_rates(k_p_cr47__mn48)*Y(jcr47)*state % rho - &
      screened_rates(k_p_cr48__he4_v45)*Y(jcr48)*state % rho - &
      screened_rates(k_p_cr48__mn49)*Y(jcr48)*state % rho - &
      screened_rates(k_p_cr49__he4_v46)*Y(jcr49)*state % rho - &
      screened_rates(k_p_cr49__mn50)*Y(jcr49)*state % rho - &
      screened_rates(k_p_cr50__he4_v47)*Y(jcr50)*state % rho - &
      screened_rates(k_p_cr50__mn51)*Y(jcr50)*state % rho - &
      screened_rates(k_p_cr51__he4_v48)*Y(jcr51)*state % rho - &
      screened_rates(k_p_cr51__mn52)*Y(jcr51)*state % rho - &
      screened_rates(k_p_cr52__he4_v49)*Y(jcr52)*state % rho - &
      screened_rates(k_p_cr52__mn53)*Y(jcr52)*state % rho - screened_rates(k_p_d__he3)* &
      Y(jd)*state % rho - screened_rates(k_p_f17__he4_o14)*Y(jf17)*state % rho - &
      screened_rates(k_p_f17__ne18)*Y(jf17)*state % rho - screened_rates(k_p_f18__he4_o15) &
      *Y(jf18)*state % rho - screened_rates(k_p_f18__ne19)*Y(jf18)*state % rho - &
      screened_rates(k_p_f19__he4_o16)*Y(jf19)*state % rho - screened_rates(k_p_f19__ne20) &
      *Y(jf19)*state % rho - screened_rates(k_p_f20__he4_o17)*Y(jf20)*state % rho &
      - screened_rates(k_p_f20__ne21)*Y(jf20)*state % rho - screened_rates(k_p_fe46__co47) &
      *Y(jfe46)*state % rho - screened_rates(k_p_fe47__co48)*Y(jfe47)*state % rho &
      - screened_rates(k_p_fe47__he4_mn44)*Y(jfe47)*state % rho - &
      screened_rates(k_p_fe48__co49)*Y(jfe48)*state % rho - &
      screened_rates(k_p_fe48__he4_mn45)*Y(jfe48)*state % rho - &
      screened_rates(k_p_fe49__co50)*Y(jfe49)*state % rho - &
      screened_rates(k_p_fe49__he4_mn46)*Y(jfe49)*state % rho - &
      screened_rates(k_p_fe50__co51)*Y(jfe50)*state % rho - &
      screened_rates(k_p_fe50__he4_mn47)*Y(jfe50)*state % rho - &
      screened_rates(k_p_fe51__co52)*Y(jfe51)*state % rho - &
      screened_rates(k_p_fe51__he4_mn48)*Y(jfe51)*state % rho - &
      screened_rates(k_p_fe52__co53)*Y(jfe52)*state % rho - &
      screened_rates(k_p_fe52__he4_mn49)*Y(jfe52)*state % rho - &
      screened_rates(k_p_fe53__co54)*Y(jfe53)*state % rho - &
      screened_rates(k_p_fe53__he4_mn50)*Y(jfe53)*state % rho - &
      screened_rates(k_p_fe54__co55)*Y(jfe54)*state % rho - &
      screened_rates(k_p_fe54__he4_mn51)*Y(jfe54)*state % rho - &
      screened_rates(k_p_fe55__co56)*Y(jfe55)*state % rho - &
      screened_rates(k_p_fe55__he4_mn52)*Y(jfe55)*state % rho - &
      screened_rates(k_p_fe56__he4_mn53)*Y(jfe56)*state % rho - &
      screened_rates(k_p_he3__he4__weak__bet_pos_)*Y(jhe3)*state % rho - &
      screened_rates(k_p_he4__d_he3)*Y(jhe4)*state % rho - 0.5d0* &
      screened_rates(k_p_he4_he4__d_be7)*Y(jhe4)**2*state % rho**2 - &
      screened_rates(k_p_k33__ca34)*Y(jk33)*state % rho - screened_rates(k_p_k34__ca35)* &
      Y(jk34)*state % rho - screened_rates(k_p_k34__he4_ar31)*Y(jk34)*state % rho &
      - screened_rates(k_p_k35__ca36)*Y(jk35)*state % rho - &
      screened_rates(k_p_k35__he4_ar32)*Y(jk35)*state % rho - &
      screened_rates(k_p_k36__ca37)*Y(jk36)*state % rho - &
      screened_rates(k_p_k36__he4_ar33)*Y(jk36)*state % rho - &
      screened_rates(k_p_k37__ca38)*Y(jk37)*state % rho - &
      screened_rates(k_p_k37__he4_ar34)*Y(jk37)*state % rho - &
      screened_rates(k_p_k38__ca39)*Y(jk38)*state % rho - &
      screened_rates(k_p_k38__he4_ar35)*Y(jk38)*state % rho - &
      screened_rates(k_p_k39__ca40)*Y(jk39)*state % rho - &
      screened_rates(k_p_k39__he4_ar36)*Y(jk39)*state % rho - &
      screened_rates(k_p_k40__ca41)*Y(jk40)*state % rho - &
      screened_rates(k_p_k40__he4_ar37)*Y(jk40)*state % rho - &
      screened_rates(k_p_k41__ca42)*Y(jk41)*state % rho - &
      screened_rates(k_p_k41__he4_ar38)*Y(jk41)*state % rho - &
      screened_rates(k_p_li7__he4_he4)*Y(jli7)*state % rho - &
      screened_rates(k_p_mg21__al22)*Y(jmg21)*state % rho - screened_rates(k_p_mg22__al23) &
      *Y(jmg22)*state % rho - screened_rates(k_p_mg22__he4_na19)*Y(jmg22)* &
      state % rho - screened_rates(k_p_mg23__al24)*Y(jmg23)*state % rho - &
      screened_rates(k_p_mg23__he4_na20)*Y(jmg23)*state % rho - &
      screened_rates(k_p_mg24__al25)*Y(jmg24)*state % rho - &
      screened_rates(k_p_mg24__he4_na21)*Y(jmg24)*state % rho - &
      screened_rates(k_p_mg25__al26)*Y(jmg25)*state % rho - &
      screened_rates(k_p_mg25__he4_na22)*Y(jmg25)*state % rho - &
      screened_rates(k_p_mg26__al27)*Y(jmg26)*state % rho - &
      screened_rates(k_p_mg26__he4_na23)*Y(jmg26)*state % rho - &
      screened_rates(k_p_mn44__fe45)*Y(jmn44)*state % rho - screened_rates(k_p_mn45__fe46) &
      *Y(jmn45)*state % rho - screened_rates(k_p_mn45__he4_cr42)*Y(jmn45)* &
      state % rho - screened_rates(k_p_mn46__fe47)*Y(jmn46)*state % rho - &
      screened_rates(k_p_mn46__he4_cr43)*Y(jmn46)*state % rho - &
      screened_rates(k_p_mn47__fe48)*Y(jmn47)*state % rho - &
      screened_rates(k_p_mn47__he4_cr44)*Y(jmn47)*state % rho - &
      screened_rates(k_p_mn48__fe49)*Y(jmn48)*state % rho - &
      screened_rates(k_p_mn48__he4_cr45)*Y(jmn48)*state % rho - &
      screened_rates(k_p_mn49__fe50)*Y(jmn49)*state % rho - &
      screened_rates(k_p_mn49__he4_cr46)*Y(jmn49)*state % rho - &
      screened_rates(k_p_mn50__fe51)*Y(jmn50)*state % rho - &
      screened_rates(k_p_mn50__he4_cr47)*Y(jmn50)*state % rho - &
      screened_rates(k_p_mn51__fe52)*Y(jmn51)*state % rho - &
      screened_rates(k_p_mn51__he4_cr48)*Y(jmn51)*state % rho - &
      screened_rates(k_p_mn52__fe53)*Y(jmn52)*state % rho - &
      screened_rates(k_p_mn52__he4_cr49)*Y(jmn52)*state % rho - &
      screened_rates(k_p_mn53__fe54)*Y(jmn53)*state % rho - &
      screened_rates(k_p_mn53__he4_cr50)*Y(jmn53)*state % rho - &
      screened_rates(k_p_mn54__fe55)*Y(jmn54)*state % rho - &
      screened_rates(k_p_mn54__he4_cr51)*Y(jmn54)*state % rho - &
      screened_rates(k_p_mn55__fe56)*Y(jmn55)*state % rho - &
      screened_rates(k_p_mn55__he4_cr52)*Y(jmn55)*state % rho - &
      screened_rates(k_p_n13__o14)*Y(jn13)*state % rho - screened_rates(k_p_n14__o15)* &
      Y(jn14)*state % rho - screened_rates(k_p_n15__he4_c12)*Y(jn15)*state % rho - &
      screened_rates(k_p_n15__o16)*Y(jn15)*state % rho - &
      screened_rates(k_p_na20__he4_ne17)*Y(jna20)*state % rho - &
      screened_rates(k_p_na20__mg21)*Y(jna20)*state % rho - &
      screened_rates(k_p_na21__he4_ne18)*Y(jna21)*state % rho - &
      screened_rates(k_p_na21__mg22)*Y(jna21)*state % rho - &
      screened_rates(k_p_na22__he4_ne19)*Y(jna22)*state % rho - &
      screened_rates(k_p_na22__mg23)*Y(jna22)*state % rho - &
      screened_rates(k_p_na23__he4_ne20)*Y(jna23)*state % rho - &
      screened_rates(k_p_na23__mg24)*Y(jna23)*state % rho - &
      screened_rates(k_p_na24__he4_ne21)*Y(jna24)*state % rho - &
      screened_rates(k_p_na24__mg25)*Y(jna24)*state % rho - screened_rates(k_p_ne18__na19) &
      *Y(jne18)*state % rho - screened_rates(k_p_ne19__he4_f16)*Y(jne19)* &
      state % rho - screened_rates(k_p_ne19__na20)*Y(jne19)*state % rho - &
      screened_rates(k_p_ne20__he4_f17)*Y(jne20)*state % rho - &
      screened_rates(k_p_ne20__na21)*Y(jne20)*state % rho - &
      screened_rates(k_p_ne21__he4_f18)*Y(jne21)*state % rho - &
      screened_rates(k_p_ne21__na22)*Y(jne21)*state % rho - &
      screened_rates(k_p_ne22__he4_f19)*Y(jne22)*state % rho - &
      screened_rates(k_p_ne22__na23)*Y(jne22)*state % rho - &
      screened_rates(k_p_ni50__he4_co47)*Y(jni50)*state % rho - &
      screened_rates(k_p_ni51__he4_co48)*Y(jni51)*state % rho - &
      screened_rates(k_p_ni52__he4_co49)*Y(jni52)*state % rho - &
      screened_rates(k_p_ni53__he4_co50)*Y(jni53)*state % rho - &
      screened_rates(k_p_ni54__he4_co51)*Y(jni54)*state % rho - &
      screened_rates(k_p_ni55__he4_co52)*Y(jni55)*state % rho - &
      screened_rates(k_p_ni56__he4_co53)*Y(jni56)*state % rho - &
      screened_rates(k_p_o15__f16)*Y(jo15)*state % rho - screened_rates(k_p_o15__he4_n12)* &
      Y(jo15)*state % rho - screened_rates(k_p_o16__f17)*Y(jo16)*state % rho - &
      screened_rates(k_p_o16__he4_n13)*Y(jo16)*state % rho - screened_rates(k_p_o17__f18)* &
      Y(jo17)*state % rho - screened_rates(k_p_o17__he4_n14)*Y(jo17)*state % rho - &
      screened_rates(k_p_o18__f19)*Y(jo18)*state % rho - screened_rates(k_p_o18__he4_n15)* &
      Y(jo18)*state % rho - screened_rates(k_p_p27__he4_si24)*Y(jp27)*state % rho &
      - screened_rates(k_p_p27__s28)*Y(jp27)*state % rho - &
      screened_rates(k_p_p28__he4_si25)*Y(jp28)*state % rho - screened_rates(k_p_p28__s29) &
      *Y(jp28)*state % rho - screened_rates(k_p_p29__he4_si26)*Y(jp29)*state % rho &
      - screened_rates(k_p_p29__s30)*Y(jp29)*state % rho - &
      screened_rates(k_p_p30__he4_si27)*Y(jp30)*state % rho - screened_rates(k_p_p30__s31) &
      *Y(jp30)*state % rho - screened_rates(k_p_p31__he4_si28)*Y(jp31)*state % rho &
      - screened_rates(k_p_p31__s32)*Y(jp31)*state % rho - &
      screened_rates(k_p_p32__he4_si29)*Y(jp32)*state % rho - screened_rates(k_p_p32__s33) &
      *Y(jp32)*state % rho - screened_rates(k_p_p33__he4_si30)*Y(jp33)*state % rho &
      - screened_rates(k_p_p33__s34)*Y(jp33)*state % rho - 2.0d0* &
      screened_rates(k_p_p__d__weak__bet_pos_)*Y(jp)*state % rho - 2.0d0* &
      screened_rates(k_p_p__d__weak__electron_capture)*Y(jp)*state % rho**2* &
      state % y_e - 2.0d0*screened_rates(k_p_p_ca38__he4_ca36)*Y(jca38)*Y(jp)* &
      state % rho**2 - 2.0d0*screened_rates(k_p_p_he4__he3_he3)*Y(jhe4)*Y(jp)* &
      state % rho**2 - 1.0d0*screened_rates(k_p_p_he4_he4__he3_be7)*Y(jhe4)**2* &
      Y(jp)*state % rho**3 - screened_rates(k_p_s28__cl29)*Y(js28)*state % rho - &
      screened_rates(k_p_s29__cl30)*Y(js29)*state % rho - screened_rates(k_p_s29__he4_p26) &
      *Y(js29)*state % rho - screened_rates(k_p_s30__cl31)*Y(js30)*state % rho - &
      screened_rates(k_p_s30__he4_p27)*Y(js30)*state % rho - screened_rates(k_p_s31__cl32) &
      *Y(js31)*state % rho - screened_rates(k_p_s31__he4_p28)*Y(js31)*state % rho &
      - screened_rates(k_p_s32__cl33)*Y(js32)*state % rho - &
      screened_rates(k_p_s32__he4_p29)*Y(js32)*state % rho - screened_rates(k_p_s33__cl34) &
      *Y(js33)*state % rho - screened_rates(k_p_s33__he4_p30)*Y(js33)*state % rho &
      - screened_rates(k_p_s34__cl35)*Y(js34)*state % rho - &
      screened_rates(k_p_s34__he4_p31)*Y(js34)*state % rho - screened_rates(k_p_s35__cl36) &
      *Y(js35)*state % rho - screened_rates(k_p_s35__he4_p32)*Y(js35)*state % rho &
      - screened_rates(k_p_sc37__he4_ca34)*Y(jsc37)*state % rho - &
      screened_rates(k_p_sc37__ti38)*Y(jsc37)*state % rho - &
      screened_rates(k_p_sc38__he4_ca35)*Y(jsc38)*state % rho - &
      screened_rates(k_p_sc38__ti39)*Y(jsc38)*state % rho - &
      screened_rates(k_p_sc39__he4_ca36)*Y(jsc39)*state % rho - &
      screened_rates(k_p_sc39__ti40)*Y(jsc39)*state % rho - &
      screened_rates(k_p_sc40__he4_ca37)*Y(jsc40)*state % rho - &
      screened_rates(k_p_sc40__ti41)*Y(jsc40)*state % rho - &
      screened_rates(k_p_sc41__he4_ca38)*Y(jsc41)*state % rho - &
      screened_rates(k_p_sc41__ti42)*Y(jsc41)*state % rho - &
      screened_rates(k_p_sc42__he4_ca39)*Y(jsc42)*state % rho - &
      screened_rates(k_p_sc42__ti43)*Y(jsc42)*state % rho - &
      screened_rates(k_p_sc43__he4_ca40)*Y(jsc43)*state % rho - &
      screened_rates(k_p_sc43__ti44)*Y(jsc43)*state % rho - &
      screened_rates(k_p_sc44__he4_ca41)*Y(jsc44)*state % rho - &
      screened_rates(k_p_sc44__ti45)*Y(jsc44)*state % rho - &
      screened_rates(k_p_sc45__he4_ca42)*Y(jsc45)*state % rho - &
      screened_rates(k_p_sc45__ti46)*Y(jsc45)*state % rho - &
      screened_rates(k_p_sc46__he4_ca43)*Y(jsc46)*state % rho - &
      screened_rates(k_p_sc46__ti47)*Y(jsc46)*state % rho - &
      screened_rates(k_p_si25__he4_al22)*Y(jsi25)*state % rho - &
      screened_rates(k_p_si25__p26)*Y(jsi25)*state % rho - &
      screened_rates(k_p_si26__he4_al23)*Y(jsi26)*state % rho - &
      screened_rates(k_p_si26__p27)*Y(jsi26)*state % rho - &
      screened_rates(k_p_si27__he4_al24)*Y(jsi27)*state % rho - &
      screened_rates(k_p_si27__p28)*Y(jsi27)*state % rho - &
      screened_rates(k_p_si28__he4_al25)*Y(jsi28)*state % rho - &
      screened_rates(k_p_si28__p29)*Y(jsi28)*state % rho - &
      screened_rates(k_p_si29__he4_al26)*Y(jsi29)*state % rho - &
      screened_rates(k_p_si29__p30)*Y(jsi29)*state % rho - &
      screened_rates(k_p_si30__he4_al27)*Y(jsi30)*state % rho - &
      screened_rates(k_p_si30__p31)*Y(jsi30)*state % rho - &
      screened_rates(k_p_ti39__he4_sc36)*Y(jti39)*state % rho - &
      screened_rates(k_p_ti39__v40)*Y(jti39)*state % rho - &
      screened_rates(k_p_ti40__he4_sc37)*Y(jti40)*state % rho - &
      screened_rates(k_p_ti40__v41)*Y(jti40)*state % rho - &
      screened_rates(k_p_ti41__he4_sc38)*Y(jti41)*state % rho - &
      screened_rates(k_p_ti41__v42)*Y(jti41)*state % rho - &
      screened_rates(k_p_ti42__he4_sc39)*Y(jti42)*state % rho - &
      screened_rates(k_p_ti42__v43)*Y(jti42)*state % rho - &
      screened_rates(k_p_ti43__he4_sc40)*Y(jti43)*state % rho - &
      screened_rates(k_p_ti43__v44)*Y(jti43)*state % rho - &
      screened_rates(k_p_ti44__he4_sc41)*Y(jti44)*state % rho - &
      screened_rates(k_p_ti44__v45)*Y(jti44)*state % rho - &
      screened_rates(k_p_ti45__he4_sc42)*Y(jti45)*state % rho - &
      screened_rates(k_p_ti45__v46)*Y(jti45)*state % rho - &
      screened_rates(k_p_ti46__he4_sc43)*Y(jti46)*state % rho - &
      screened_rates(k_p_ti46__v47)*Y(jti46)*state % rho - &
      screened_rates(k_p_ti47__he4_sc44)*Y(jti47)*state % rho - &
      screened_rates(k_p_ti47__v48)*Y(jti47)*state % rho - &
      screened_rates(k_p_ti48__he4_sc45)*Y(jti48)*state % rho - &
      screened_rates(k_p_ti48__v49)*Y(jti48)*state % rho - screened_rates(k_p_v41__cr42)* &
      Y(jv41)*state % rho - screened_rates(k_p_v41__he4_ti38)*Y(jv41)*state % rho &
      - screened_rates(k_p_v42__cr43)*Y(jv42)*state % rho - &
      screened_rates(k_p_v42__he4_ti39)*Y(jv42)*state % rho - &
      screened_rates(k_p_v43__cr44)*Y(jv43)*state % rho - &
      screened_rates(k_p_v43__he4_ti40)*Y(jv43)*state % rho - &
      screened_rates(k_p_v44__cr45)*Y(jv44)*state % rho - &
      screened_rates(k_p_v44__he4_ti41)*Y(jv44)*state % rho - &
      screened_rates(k_p_v45__cr46)*Y(jv45)*state % rho - &
      screened_rates(k_p_v45__he4_ti42)*Y(jv45)*state % rho - &
      screened_rates(k_p_v46__cr47)*Y(jv46)*state % rho - &
      screened_rates(k_p_v46__he4_ti43)*Y(jv46)*state % rho - &
      screened_rates(k_p_v47__cr48)*Y(jv47)*state % rho - &
      screened_rates(k_p_v47__he4_ti44)*Y(jv47)*state % rho - &
      screened_rates(k_p_v48__cr49)*Y(jv48)*state % rho - &
      screened_rates(k_p_v48__he4_ti45)*Y(jv48)*state % rho - &
      screened_rates(k_p_v49__cr50)*Y(jv49)*state % rho - &
      screened_rates(k_p_v49__he4_ti46)*Y(jv49)*state % rho - &
      screened_rates(k_p_v50__cr51)*Y(jv50)*state % rho - &
      screened_rates(k_p_v50__he4_ti47)*Y(jv50)*state % rho &
       )
    call set_jac_entry(state, jp, jp, scratch)

    scratch = (&
      screened_rates(k_d_be7__p_he4_he4)*Y(jbe7)*state % rho + screened_rates(k_d_he3__p_he4)* &
      Y(jhe3)*state % rho - screened_rates(k_p_d__he3)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jd, scratch)

    scratch = (&
      screened_rates(k_d_he3__p_he4)*Y(jd)*state % rho + screened_rates(k_he3__p_d) + 2.0d0* &
      screened_rates(k_he3_be7__p_p_he4_he4)*Y(jbe7)*state % rho + 2.0d0* &
      screened_rates(k_he3_he3__p_p_he4)*Y(jhe3)*state % rho - &
      screened_rates(k_p_he3__he4__weak__bet_pos_)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jhe3, scratch)

    scratch = (&
      screened_rates(k_he4_al22__p_si25)*Y(jal22)*state % rho + &
      screened_rates(k_he4_al23__p_si26)*Y(jal23)*state % rho + &
      screened_rates(k_he4_al24__p_si27)*Y(jal24)*state % rho + &
      screened_rates(k_he4_al25__p_si28)*Y(jal25)*state % rho + &
      screened_rates(k_he4_al26__p_si29)*Y(jal26)*state % rho + &
      screened_rates(k_he4_al27__p_si30)*Y(jal27)*state % rho + &
      screened_rates(k_he4_ar31__p_k34)*Y(jar31)*state % rho + &
      screened_rates(k_he4_ar32__p_k35)*Y(jar32)*state % rho + &
      screened_rates(k_he4_ar33__p_k36)*Y(jar33)*state % rho + &
      screened_rates(k_he4_ar34__p_k37)*Y(jar34)*state % rho + &
      screened_rates(k_he4_ar35__p_k38)*Y(jar35)*state % rho + &
      screened_rates(k_he4_ar36__p_k39)*Y(jar36)*state % rho + &
      screened_rates(k_he4_ar37__p_k40)*Y(jar37)*state % rho + &
      screened_rates(k_he4_ar38__p_k41)*Y(jar38)*state % rho + &
      screened_rates(k_he4_c12__p_n15)*Y(jc12)*state % rho + &
      screened_rates(k_he4_ca34__p_sc37)*Y(jca34)*state % rho + &
      screened_rates(k_he4_ca35__p_sc38)*Y(jca35)*state % rho + 2.0d0* &
      screened_rates(k_he4_ca36__p_p_ca38)*Y(jca36)*state % rho + &
      screened_rates(k_he4_ca36__p_sc39)*Y(jca36)*state % rho + &
      screened_rates(k_he4_ca37__p_sc40)*Y(jca37)*state % rho + &
      screened_rates(k_he4_ca38__p_sc41)*Y(jca38)*state % rho + &
      screened_rates(k_he4_ca39__p_sc42)*Y(jca39)*state % rho + &
      screened_rates(k_he4_ca40__p_sc43)*Y(jca40)*state % rho + &
      screened_rates(k_he4_ca41__p_sc44)*Y(jca41)*state % rho + &
      screened_rates(k_he4_ca42__p_sc45)*Y(jca42)*state % rho + &
      screened_rates(k_he4_ca43__p_sc46)*Y(jca43)*state % rho + &
      screened_rates(k_he4_cl29__p_ar32)*Y(jcl29)*state % rho + &
      screened_rates(k_he4_cl30__p_ar33)*Y(jcl30)*state % rho + &
      screened_rates(k_he4_cl31__p_ar34)*Y(jcl31)*state % rho + &
      screened_rates(k_he4_cl32__p_ar35)*Y(jcl32)*state % rho + &
      screened_rates(k_he4_cl33__p_ar36)*Y(jcl33)*state % rho + &
      screened_rates(k_he4_cl34__p_ar37)*Y(jcl34)*state % rho + &
      screened_rates(k_he4_cl35__p_ar38)*Y(jcl35)*state % rho + &
      screened_rates(k_he4_cl36__p_ar39)*Y(jcl36)*state % rho + &
      screened_rates(k_he4_co47__p_ni50)*Y(jco47)*state % rho + &
      screened_rates(k_he4_co48__p_ni51)*Y(jco48)*state % rho + &
      screened_rates(k_he4_co49__p_ni52)*Y(jco49)*state % rho + &
      screened_rates(k_he4_co50__p_ni53)*Y(jco50)*state % rho + &
      screened_rates(k_he4_co51__p_ni54)*Y(jco51)*state % rho + &
      screened_rates(k_he4_co52__p_ni55)*Y(jco52)*state % rho + &
      screened_rates(k_he4_co53__p_ni56)*Y(jco53)*state % rho + &
      screened_rates(k_he4_cr42__p_mn45)*Y(jcr42)*state % rho + &
      screened_rates(k_he4_cr43__p_mn46)*Y(jcr43)*state % rho + &
      screened_rates(k_he4_cr44__p_mn47)*Y(jcr44)*state % rho + &
      screened_rates(k_he4_cr45__p_mn48)*Y(jcr45)*state % rho + &
      screened_rates(k_he4_cr46__p_mn49)*Y(jcr46)*state % rho + &
      screened_rates(k_he4_cr47__p_mn50)*Y(jcr47)*state % rho + &
      screened_rates(k_he4_cr48__p_mn51)*Y(jcr48)*state % rho + &
      screened_rates(k_he4_cr49__p_mn52)*Y(jcr49)*state % rho + &
      screened_rates(k_he4_cr50__p_mn53)*Y(jcr50)*state % rho + &
      screened_rates(k_he4_cr51__p_mn54)*Y(jcr51)*state % rho + &
      screened_rates(k_he4_cr52__p_mn55)*Y(jcr52)*state % rho + &
      screened_rates(k_he4_f16__p_ne19)*Y(jf16)*state % rho + &
      screened_rates(k_he4_f17__p_ne20)*Y(jf17)*state % rho + &
      screened_rates(k_he4_f18__p_ne21)*Y(jf18)*state % rho + &
      screened_rates(k_he4_f19__p_ne22)*Y(jf19)*state % rho + &
      screened_rates(k_he4_fe45__p_co48)*Y(jfe45)*state % rho + &
      screened_rates(k_he4_fe46__p_co49)*Y(jfe46)*state % rho + &
      screened_rates(k_he4_fe47__p_co50)*Y(jfe47)*state % rho + &
      screened_rates(k_he4_fe48__p_co51)*Y(jfe48)*state % rho + &
      screened_rates(k_he4_fe49__p_co52)*Y(jfe49)*state % rho + &
      screened_rates(k_he4_fe50__p_co53)*Y(jfe50)*state % rho + &
      screened_rates(k_he4_fe51__p_co54)*Y(jfe51)*state % rho + &
      screened_rates(k_he4_fe52__p_co55)*Y(jfe52)*state % rho + &
      screened_rates(k_he4_fe53__p_co56)*Y(jfe53)*state % rho + 1.0d0* &
      screened_rates(k_he4_he4__p_li7)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_k33__p_ca36)*Y(jk33)*state % rho + &
      screened_rates(k_he4_k34__p_ca37)*Y(jk34)*state % rho + &
      screened_rates(k_he4_k35__p_ca38)*Y(jk35)*state % rho + &
      screened_rates(k_he4_k36__p_ca39)*Y(jk36)*state % rho + &
      screened_rates(k_he4_k37__p_ca40)*Y(jk37)*state % rho + &
      screened_rates(k_he4_k38__p_ca41)*Y(jk38)*state % rho + &
      screened_rates(k_he4_k39__p_ca42)*Y(jk39)*state % rho + &
      screened_rates(k_he4_k40__p_ca43)*Y(jk40)*state % rho + &
      screened_rates(k_he4_k41__p_ca44)*Y(jk41)*state % rho + &
      screened_rates(k_he4_mg21__p_al24)*Y(jmg21)*state % rho + &
      screened_rates(k_he4_mg22__p_al25)*Y(jmg22)*state % rho + &
      screened_rates(k_he4_mg23__p_al26)*Y(jmg23)*state % rho + &
      screened_rates(k_he4_mg24__p_al27)*Y(jmg24)*state % rho + &
      screened_rates(k_he4_mg25__p_al28)*Y(jmg25)*state % rho + &
      screened_rates(k_he4_mn44__p_fe47)*Y(jmn44)*state % rho + &
      screened_rates(k_he4_mn45__p_fe48)*Y(jmn45)*state % rho + &
      screened_rates(k_he4_mn46__p_fe49)*Y(jmn46)*state % rho + &
      screened_rates(k_he4_mn47__p_fe50)*Y(jmn47)*state % rho + &
      screened_rates(k_he4_mn48__p_fe51)*Y(jmn48)*state % rho + &
      screened_rates(k_he4_mn49__p_fe52)*Y(jmn49)*state % rho + &
      screened_rates(k_he4_mn50__p_fe53)*Y(jmn50)*state % rho + &
      screened_rates(k_he4_mn51__p_fe54)*Y(jmn51)*state % rho + &
      screened_rates(k_he4_mn52__p_fe55)*Y(jmn52)*state % rho + &
      screened_rates(k_he4_mn53__p_fe56)*Y(jmn53)*state % rho + &
      screened_rates(k_he4_n12__p_o15)*Y(jn12)*state % rho + &
      screened_rates(k_he4_n13__p_o16)*Y(jn13)*state % rho + &
      screened_rates(k_he4_n14__p_o17)*Y(jn14)*state % rho + &
      screened_rates(k_he4_n15__p_o18)*Y(jn15)*state % rho + &
      screened_rates(k_he4_na19__p_mg22)*Y(jna19)*state % rho + &
      screened_rates(k_he4_na20__p_mg23)*Y(jna20)*state % rho + &
      screened_rates(k_he4_na21__p_mg24)*Y(jna21)*state % rho + &
      screened_rates(k_he4_na22__p_mg25)*Y(jna22)*state % rho + &
      screened_rates(k_he4_na23__p_mg26)*Y(jna23)*state % rho + &
      screened_rates(k_he4_ne17__p_na20)*Y(jne17)*state % rho + &
      screened_rates(k_he4_ne18__p_na21)*Y(jne18)*state % rho + &
      screened_rates(k_he4_ne19__p_na22)*Y(jne19)*state % rho + &
      screened_rates(k_he4_ne20__p_na23)*Y(jne20)*state % rho + &
      screened_rates(k_he4_ne21__p_na24)*Y(jne21)*state % rho + &
      screened_rates(k_he4_o14__p_f17)*Y(jo14)*state % rho + &
      screened_rates(k_he4_o15__p_f18)*Y(jo15)*state % rho + &
      screened_rates(k_he4_o16__p_f19)*Y(jo16)*state % rho + &
      screened_rates(k_he4_o17__p_f20)*Y(jo17)*state % rho + &
      screened_rates(k_he4_p26__p_s29)*Y(jp26)*state % rho + &
      screened_rates(k_he4_p27__p_s30)*Y(jp27)*state % rho + &
      screened_rates(k_he4_p28__p_s31)*Y(jp28)*state % rho + &
      screened_rates(k_he4_p29__p_s32)*Y(jp29)*state % rho + &
      screened_rates(k_he4_p30__p_s33)*Y(jp30)*state % rho + &
      screened_rates(k_he4_p31__p_s34)*Y(jp31)*state % rho + &
      screened_rates(k_he4_p32__p_s35)*Y(jp32)*state % rho + &
      screened_rates(k_he4_s28__p_cl31)*Y(js28)*state % rho + &
      screened_rates(k_he4_s29__p_cl32)*Y(js29)*state % rho + &
      screened_rates(k_he4_s30__p_cl33)*Y(js30)*state % rho + &
      screened_rates(k_he4_s31__p_cl34)*Y(js31)*state % rho + &
      screened_rates(k_he4_s32__p_cl35)*Y(js32)*state % rho + &
      screened_rates(k_he4_s33__p_cl36)*Y(js33)*state % rho + &
      screened_rates(k_he4_s34__p_cl37)*Y(js34)*state % rho + &
      screened_rates(k_he4_sc36__p_ti39)*Y(jsc36)*state % rho + &
      screened_rates(k_he4_sc37__p_ti40)*Y(jsc37)*state % rho + &
      screened_rates(k_he4_sc38__p_ti41)*Y(jsc38)*state % rho + &
      screened_rates(k_he4_sc39__p_ti42)*Y(jsc39)*state % rho + &
      screened_rates(k_he4_sc40__p_ti43)*Y(jsc40)*state % rho + &
      screened_rates(k_he4_sc41__p_ti44)*Y(jsc41)*state % rho + &
      screened_rates(k_he4_sc42__p_ti45)*Y(jsc42)*state % rho + &
      screened_rates(k_he4_sc43__p_ti46)*Y(jsc43)*state % rho + &
      screened_rates(k_he4_sc44__p_ti47)*Y(jsc44)*state % rho + &
      screened_rates(k_he4_sc45__p_ti48)*Y(jsc45)*state % rho + &
      screened_rates(k_he4_si24__p_p27)*Y(jsi24)*state % rho + &
      screened_rates(k_he4_si25__p_p28)*Y(jsi25)*state % rho + &
      screened_rates(k_he4_si26__p_p29)*Y(jsi26)*state % rho + &
      screened_rates(k_he4_si27__p_p30)*Y(jsi27)*state % rho + &
      screened_rates(k_he4_si28__p_p31)*Y(jsi28)*state % rho + &
      screened_rates(k_he4_si29__p_p32)*Y(jsi29)*state % rho + &
      screened_rates(k_he4_si30__p_p33)*Y(jsi30)*state % rho + &
      screened_rates(k_he4_ti38__p_v41)*Y(jti38)*state % rho + &
      screened_rates(k_he4_ti39__p_v42)*Y(jti39)*state % rho + &
      screened_rates(k_he4_ti40__p_v43)*Y(jti40)*state % rho + &
      screened_rates(k_he4_ti41__p_v44)*Y(jti41)*state % rho + &
      screened_rates(k_he4_ti42__p_v45)*Y(jti42)*state % rho + &
      screened_rates(k_he4_ti43__p_v46)*Y(jti43)*state % rho + &
      screened_rates(k_he4_ti44__p_v47)*Y(jti44)*state % rho + &
      screened_rates(k_he4_ti45__p_v48)*Y(jti45)*state % rho + &
      screened_rates(k_he4_ti46__p_v49)*Y(jti46)*state % rho + &
      screened_rates(k_he4_ti47__p_v50)*Y(jti47)*state % rho + &
      screened_rates(k_he4_v40__p_cr43)*Y(jv40)*state % rho + &
      screened_rates(k_he4_v41__p_cr44)*Y(jv41)*state % rho + &
      screened_rates(k_he4_v42__p_cr45)*Y(jv42)*state % rho + &
      screened_rates(k_he4_v43__p_cr46)*Y(jv43)*state % rho + &
      screened_rates(k_he4_v44__p_cr47)*Y(jv44)*state % rho + &
      screened_rates(k_he4_v45__p_cr48)*Y(jv45)*state % rho + &
      screened_rates(k_he4_v46__p_cr49)*Y(jv46)*state % rho + &
      screened_rates(k_he4_v47__p_cr50)*Y(jv47)*state % rho + &
      screened_rates(k_he4_v48__p_cr51)*Y(jv48)*state % rho + &
      screened_rates(k_he4_v49__p_cr52)*Y(jv49)*state % rho - &
      screened_rates(k_p_he4__d_he3)*Y(jp)*state % rho - 1.0d0* &
      screened_rates(k_p_he4_he4__d_be7)*Y(jhe4)*Y(jp)*state % rho**2 - &
      screened_rates(k_p_p_he4__he3_he3)*Y(jp)**2*state % rho**2 - 1.0d0* &
      screened_rates(k_p_p_he4_he4__he3_be7)*Y(jhe4)*Y(jp)**2*state % rho**3 &
       )
    call set_jac_entry(state, jp, jhe4, scratch)

    scratch = (&
      -screened_rates(k_p_li7__he4_he4)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jli7, scratch)

    scratch = (&
      screened_rates(k_d_be7__p_he4_he4)*Y(jd)*state % rho + 2.0d0* &
      screened_rates(k_he3_be7__p_p_he4_he4)*Y(jhe3)*state % rho - &
      screened_rates(k_p_be7__b8)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jbe7, scratch)

    scratch = (&
      screened_rates(k_b8__p_be7) &
       )
    call set_jac_entry(state, jp, jb8, scratch)

    scratch = (&
      screened_rates(k_he4_c12__p_n15)*Y(jhe4)*state % rho - screened_rates(k_p_c12__n13)*Y(jp) &
      *state % rho &
       )
    call set_jac_entry(state, jp, jc12, scratch)

    scratch = (&
      -screened_rates(k_p_c13__n14)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jc13, scratch)

    scratch = (&
      screened_rates(k_he4_n12__p_o15)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp, jn12, scratch)

    scratch = (&
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*state % rho + screened_rates(k_n13__p_c12) - &
      screened_rates(k_p_n13__o14)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jn13, scratch)

    scratch = (&
      screened_rates(k_he4_n14__p_o17)*Y(jhe4)*state % rho + screened_rates(k_n14__p_c13) - &
      screened_rates(k_p_n14__o15)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jn14, scratch)

    scratch = (&
      screened_rates(k_he4_n15__p_o18)*Y(jhe4)*state % rho - screened_rates(k_p_n15__he4_c12)* &
      Y(jp)*state % rho - screened_rates(k_p_n15__o16)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jn15, scratch)

    scratch = (&
      screened_rates(k_he4_o14__p_f17)*Y(jhe4)*state % rho + screened_rates(k_o14__p_n13) &
       )
    call set_jac_entry(state, jp, jo14, scratch)

    scratch = (&
      screened_rates(k_he4_o15__p_f18)*Y(jhe4)*state % rho + screened_rates(k_o15__p_n14) - &
      screened_rates(k_p_o15__f16)*Y(jp)*state % rho - screened_rates(k_p_o15__he4_n12)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jo15, scratch)

    scratch = (&
      screened_rates(k_he4_o16__p_f19)*Y(jhe4)*state % rho + screened_rates(k_o16__p_n15) - &
      screened_rates(k_p_o16__f17)*Y(jp)*state % rho - screened_rates(k_p_o16__he4_n13)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jo16, scratch)

    scratch = (&
      screened_rates(k_he4_o17__p_f20)*Y(jhe4)*state % rho - screened_rates(k_p_o17__f18)*Y(jp) &
      *state % rho - screened_rates(k_p_o17__he4_n14)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jo17, scratch)

    scratch = (&
      -screened_rates(k_p_o18__f19)*Y(jp)*state % rho - screened_rates(k_p_o18__he4_n15)*Y(jp)* &
      state % rho &
       )
    call set_jac_entry(state, jp, jo18, scratch)

    scratch = (&
      screened_rates(k_f16__p_o15) + screened_rates(k_he4_f16__p_ne19)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp, jf16, scratch)

    scratch = (&
      screened_rates(k_f17__p_o16) + screened_rates(k_he4_f17__p_ne20)*Y(jhe4)*state % rho - &
      screened_rates(k_p_f17__he4_o14)*Y(jp)*state % rho - screened_rates(k_p_f17__ne18)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jf17, scratch)

    scratch = (&
      screened_rates(k_f18__p_o17) + screened_rates(k_he4_f18__p_ne21)*Y(jhe4)*state % rho - &
      screened_rates(k_p_f18__he4_o15)*Y(jp)*state % rho - screened_rates(k_p_f18__ne19)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jf18, scratch)

    scratch = (&
      screened_rates(k_f19__p_o18) + screened_rates(k_he4_f19__p_ne22)*Y(jhe4)*state % rho - &
      screened_rates(k_p_f19__he4_o16)*Y(jp)*state % rho - screened_rates(k_p_f19__ne20)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jf19, scratch)

    scratch = (&
      -screened_rates(k_p_f20__he4_o17)*Y(jp)*state % rho - screened_rates(k_p_f20__ne21)*Y(jp) &
      *state % rho &
       )
    call set_jac_entry(state, jp, jf20, scratch)

    scratch = (&
      screened_rates(k_he4_ne17__p_na20)*Y(jhe4)*state % rho + &
      screened_rates(k_ne17__p_o16__weak__wc12) &
       )
    call set_jac_entry(state, jp, jne17, scratch)

    scratch = (&
      screened_rates(k_he4_ne18__p_na21)*Y(jhe4)*state % rho + screened_rates(k_ne18__p_f17) - &
      screened_rates(k_p_ne18__na19)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jne18, scratch)

    scratch = (&
      screened_rates(k_he4_ne19__p_na22)*Y(jhe4)*state % rho + screened_rates(k_ne19__p_f18) - &
      screened_rates(k_p_ne19__he4_f16)*Y(jp)*state % rho - screened_rates(k_p_ne19__na20) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jne19, scratch)

    scratch = (&
      screened_rates(k_he4_ne20__p_na23)*Y(jhe4)*state % rho + screened_rates(k_ne20__p_f19) - &
      screened_rates(k_p_ne20__he4_f17)*Y(jp)*state % rho - screened_rates(k_p_ne20__na21) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jne20, scratch)

    scratch = (&
      screened_rates(k_he4_ne21__p_na24)*Y(jhe4)*state % rho + screened_rates(k_ne21__p_f20) - &
      screened_rates(k_p_ne21__he4_f18)*Y(jp)*state % rho - screened_rates(k_p_ne21__na22) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jne21, scratch)

    scratch = (&
      -screened_rates(k_p_ne22__he4_f19)*Y(jp)*state % rho - screened_rates(k_p_ne22__na23)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jne22, scratch)

    scratch = (&
      screened_rates(k_he4_na19__p_mg22)*Y(jhe4)*state % rho + screened_rates(k_na19__p_ne18) &
       )
    call set_jac_entry(state, jp, jna19, scratch)

    scratch = (&
      screened_rates(k_he4_na20__p_mg23)*Y(jhe4)*state % rho + screened_rates(k_na20__p_ne19) - &
      screened_rates(k_p_na20__he4_ne17)*Y(jp)*state % rho - &
      screened_rates(k_p_na20__mg21)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jna20, scratch)

    scratch = (&
      screened_rates(k_he4_na21__p_mg24)*Y(jhe4)*state % rho + screened_rates(k_na21__p_ne20) - &
      screened_rates(k_p_na21__he4_ne18)*Y(jp)*state % rho - &
      screened_rates(k_p_na21__mg22)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jna21, scratch)

    scratch = (&
      screened_rates(k_he4_na22__p_mg25)*Y(jhe4)*state % rho + screened_rates(k_na22__p_ne21) - &
      screened_rates(k_p_na22__he4_ne19)*Y(jp)*state % rho - &
      screened_rates(k_p_na22__mg23)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jna22, scratch)

    scratch = (&
      screened_rates(k_he4_na23__p_mg26)*Y(jhe4)*state % rho + screened_rates(k_na23__p_ne22) - &
      screened_rates(k_p_na23__he4_ne20)*Y(jp)*state % rho - &
      screened_rates(k_p_na23__mg24)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jna23, scratch)

    scratch = (&
      -screened_rates(k_p_na24__he4_ne21)*Y(jp)*state % rho - screened_rates(k_p_na24__mg25)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jna24, scratch)

    scratch = (&
      screened_rates(k_he4_mg21__p_al24)*Y(jhe4)*state % rho + screened_rates(k_mg21__p_na20) + &
      screened_rates(k_mg21__p_ne20__weak__wc12) - screened_rates(k_p_mg21__al22)*Y(jp) &
      *state % rho &
       )
    call set_jac_entry(state, jp, jmg21, scratch)

    scratch = (&
      screened_rates(k_he4_mg22__p_al25)*Y(jhe4)*state % rho + screened_rates(k_mg22__p_na21) - &
      screened_rates(k_p_mg22__al23)*Y(jp)*state % rho - &
      screened_rates(k_p_mg22__he4_na19)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jmg22, scratch)

    scratch = (&
      screened_rates(k_he4_mg23__p_al26)*Y(jhe4)*state % rho + screened_rates(k_mg23__p_na22) - &
      screened_rates(k_p_mg23__al24)*Y(jp)*state % rho - &
      screened_rates(k_p_mg23__he4_na20)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jmg23, scratch)

    scratch = (&
      screened_rates(k_he4_mg24__p_al27)*Y(jhe4)*state % rho + screened_rates(k_mg24__p_na23) - &
      screened_rates(k_p_mg24__al25)*Y(jp)*state % rho - &
      screened_rates(k_p_mg24__he4_na21)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jmg24, scratch)

    scratch = (&
      screened_rates(k_he4_mg25__p_al28)*Y(jhe4)*state % rho + screened_rates(k_mg25__p_na24) - &
      screened_rates(k_p_mg25__al26)*Y(jp)*state % rho - &
      screened_rates(k_p_mg25__he4_na22)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jmg25, scratch)

    scratch = (&
      -screened_rates(k_p_mg26__al27)*Y(jp)*state % rho - screened_rates(k_p_mg26__he4_na23)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jmg26, scratch)

    scratch = (&
      screened_rates(k_al22__p_mg21) + screened_rates(k_al22__p_na21__weak__wc17) + 2.0d0* &
      screened_rates(k_al22__p_p_ne20__weak__wc12) + screened_rates(k_he4_al22__p_si25)* &
      Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp, jal22, scratch)

    scratch = (&
      screened_rates(k_al23__p_mg22) + screened_rates(k_al23__p_na22__weak__wc12) + &
      screened_rates(k_he4_al23__p_si26)*Y(jhe4)*state % rho - &
      screened_rates(k_p_al23__si24)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jal23, scratch)

    scratch = (&
      screened_rates(k_al24__p_mg23) + screened_rates(k_al24__p_na23__weak__wc12) + &
      screened_rates(k_he4_al24__p_si27)*Y(jhe4)*state % rho - &
      screened_rates(k_p_al24__he4_mg21)*Y(jp)*state % rho - &
      screened_rates(k_p_al24__si25)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jal24, scratch)

    scratch = (&
      screened_rates(k_al25__p_mg24) + screened_rates(k_he4_al25__p_si28)*Y(jhe4)*state % rho - &
      screened_rates(k_p_al25__he4_mg22)*Y(jp)*state % rho - &
      screened_rates(k_p_al25__si26)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jal25, scratch)

    scratch = (&
      screened_rates(k_al26__p_mg25) + screened_rates(k_he4_al26__p_si29)*Y(jhe4)*state % rho - &
      screened_rates(k_p_al26__he4_mg23)*Y(jp)*state % rho - &
      screened_rates(k_p_al26__si27)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jal26, scratch)

    scratch = (&
      screened_rates(k_al27__p_mg26) + screened_rates(k_he4_al27__p_si30)*Y(jhe4)*state % rho - &
      screened_rates(k_p_al27__he4_mg24)*Y(jp)*state % rho - &
      screened_rates(k_p_al27__si28)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jal27, scratch)

    scratch = (&
      -screened_rates(k_p_al28__he4_mg25)*Y(jp)*state % rho - screened_rates(k_p_al28__si29)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jal28, scratch)

    scratch = (&
      screened_rates(k_he4_si24__p_p27)*Y(jhe4)*state % rho + screened_rates(k_si24__p_al23) + &
      screened_rates(k_si24__p_mg23__weak__wc12) &
       )
    call set_jac_entry(state, jp, jsi24, scratch)

    scratch = (&
      screened_rates(k_he4_si25__p_p28)*Y(jhe4)*state % rho - screened_rates(k_p_si25__he4_al22)* &
      Y(jp)*state % rho - screened_rates(k_p_si25__p26)*Y(jp)*state % rho + &
      screened_rates(k_si25__p_al24) + screened_rates(k_si25__p_mg24__weak__wc12) &
       )
    call set_jac_entry(state, jp, jsi25, scratch)

    scratch = (&
      screened_rates(k_he4_si26__p_p29)*Y(jhe4)*state % rho - screened_rates(k_p_si26__he4_al23)* &
      Y(jp)*state % rho - screened_rates(k_p_si26__p27)*Y(jp)*state % rho + &
      screened_rates(k_si26__p_al25) &
       )
    call set_jac_entry(state, jp, jsi26, scratch)

    scratch = (&
      screened_rates(k_he4_si27__p_p30)*Y(jhe4)*state % rho - screened_rates(k_p_si27__he4_al24)* &
      Y(jp)*state % rho - screened_rates(k_p_si27__p28)*Y(jp)*state % rho + &
      screened_rates(k_si27__p_al26) &
       )
    call set_jac_entry(state, jp, jsi27, scratch)

    scratch = (&
      screened_rates(k_he4_si28__p_p31)*Y(jhe4)*state % rho - screened_rates(k_p_si28__he4_al25)* &
      Y(jp)*state % rho - screened_rates(k_p_si28__p29)*Y(jp)*state % rho + &
      screened_rates(k_si28__p_al27) &
       )
    call set_jac_entry(state, jp, jsi28, scratch)

    scratch = (&
      screened_rates(k_he4_si29__p_p32)*Y(jhe4)*state % rho - screened_rates(k_p_si29__he4_al26)* &
      Y(jp)*state % rho - screened_rates(k_p_si29__p30)*Y(jp)*state % rho + &
      screened_rates(k_si29__p_al28) &
       )
    call set_jac_entry(state, jp, jsi29, scratch)

    scratch = (&
      screened_rates(k_he4_si30__p_p33)*Y(jhe4)*state % rho - screened_rates(k_p_si30__he4_al27)* &
      Y(jp)*state % rho - screened_rates(k_p_si30__p31)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jsi30, scratch)

    scratch = (&
      screened_rates(k_he4_p26__p_s29)*Y(jhe4)*state % rho + screened_rates(k_p26__p_si25) &
       )
    call set_jac_entry(state, jp, jp26, scratch)

    scratch = (&
      screened_rates(k_he4_p27__p_s30)*Y(jhe4)*state % rho + &
      screened_rates(k_p27__p_al26__weak__wc12) + screened_rates(k_p27__p_si26) - &
      screened_rates(k_p_p27__he4_si24)*Y(jp)*state % rho - screened_rates(k_p_p27__s28)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jp27, scratch)

    scratch = (&
      screened_rates(k_he4_p28__p_s31)*Y(jhe4)*state % rho + &
      screened_rates(k_p28__p_al27__weak__wc12) + screened_rates(k_p28__p_si27) - &
      screened_rates(k_p_p28__he4_si25)*Y(jp)*state % rho - screened_rates(k_p_p28__s29)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jp28, scratch)

    scratch = (&
      screened_rates(k_he4_p29__p_s32)*Y(jhe4)*state % rho + screened_rates(k_p29__p_si28) - &
      screened_rates(k_p_p29__he4_si26)*Y(jp)*state % rho - screened_rates(k_p_p29__s30)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jp29, scratch)

    scratch = (&
      screened_rates(k_he4_p30__p_s33)*Y(jhe4)*state % rho + screened_rates(k_p30__p_si29) - &
      screened_rates(k_p_p30__he4_si27)*Y(jp)*state % rho - screened_rates(k_p_p30__s31)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jp30, scratch)

    scratch = (&
      screened_rates(k_he4_p31__p_s34)*Y(jhe4)*state % rho + screened_rates(k_p31__p_si30) - &
      screened_rates(k_p_p31__he4_si28)*Y(jp)*state % rho - screened_rates(k_p_p31__s32)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jp31, scratch)

    scratch = (&
      screened_rates(k_he4_p32__p_s35)*Y(jhe4)*state % rho - screened_rates(k_p_p32__he4_si29)* &
      Y(jp)*state % rho - screened_rates(k_p_p32__s33)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jp32, scratch)

    scratch = (&
      -screened_rates(k_p_p33__he4_si30)*Y(jp)*state % rho - screened_rates(k_p_p33__s34)*Y(jp) &
      *state % rho &
       )
    call set_jac_entry(state, jp, jp33, scratch)

    scratch = (&
      screened_rates(k_he4_s28__p_cl31)*Y(jhe4)*state % rho - screened_rates(k_p_s28__cl29)* &
      Y(jp)*state % rho + screened_rates(k_s28__p_p27) + &
      screened_rates(k_s28__p_si27__weak__wc12) &
       )
    call set_jac_entry(state, jp, js28, scratch)

    scratch = (&
      screened_rates(k_he4_s29__p_cl32)*Y(jhe4)*state % rho - screened_rates(k_p_s29__cl30)* &
      Y(jp)*state % rho - screened_rates(k_p_s29__he4_p26)*Y(jp)*state % rho + &
      screened_rates(k_s29__p_p28) + screened_rates(k_s29__p_si28__weak__wc12) &
       )
    call set_jac_entry(state, jp, js29, scratch)

    scratch = (&
      screened_rates(k_he4_s30__p_cl33)*Y(jhe4)*state % rho - screened_rates(k_p_s30__cl31)* &
      Y(jp)*state % rho - screened_rates(k_p_s30__he4_p27)*Y(jp)*state % rho + &
      screened_rates(k_s30__p_p29) &
       )
    call set_jac_entry(state, jp, js30, scratch)

    scratch = (&
      screened_rates(k_he4_s31__p_cl34)*Y(jhe4)*state % rho - screened_rates(k_p_s31__cl32)* &
      Y(jp)*state % rho - screened_rates(k_p_s31__he4_p28)*Y(jp)*state % rho + &
      screened_rates(k_s31__p_p30) &
       )
    call set_jac_entry(state, jp, js31, scratch)

    scratch = (&
      screened_rates(k_he4_s32__p_cl35)*Y(jhe4)*state % rho - screened_rates(k_p_s32__cl33)* &
      Y(jp)*state % rho - screened_rates(k_p_s32__he4_p29)*Y(jp)*state % rho + &
      screened_rates(k_s32__p_p31) &
       )
    call set_jac_entry(state, jp, js32, scratch)

    scratch = (&
      screened_rates(k_he4_s33__p_cl36)*Y(jhe4)*state % rho - screened_rates(k_p_s33__cl34)* &
      Y(jp)*state % rho - screened_rates(k_p_s33__he4_p30)*Y(jp)*state % rho + &
      screened_rates(k_s33__p_p32) &
       )
    call set_jac_entry(state, jp, js33, scratch)

    scratch = (&
      screened_rates(k_he4_s34__p_cl37)*Y(jhe4)*state % rho - screened_rates(k_p_s34__cl35)* &
      Y(jp)*state % rho - screened_rates(k_p_s34__he4_p31)*Y(jp)*state % rho + &
      screened_rates(k_s34__p_p33) &
       )
    call set_jac_entry(state, jp, js34, scratch)

    scratch = (&
      -screened_rates(k_p_s35__cl36)*Y(jp)*state % rho - screened_rates(k_p_s35__he4_p32)*Y(jp) &
      *state % rho &
       )
    call set_jac_entry(state, jp, js35, scratch)

    scratch = (&
      screened_rates(k_cl29__p_s28) + screened_rates(k_he4_cl29__p_ar32)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp, jcl29, scratch)

    scratch = (&
      screened_rates(k_cl30__p_s29) + screened_rates(k_he4_cl30__p_ar33)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cl30__ar31)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jcl30, scratch)

    scratch = (&
      screened_rates(k_cl31__p_p30__weak__wc12) + screened_rates(k_cl31__p_s30) + &
      screened_rates(k_he4_cl31__p_ar34)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cl31__ar32)*Y(jp)*state % rho - screened_rates(k_p_cl31__he4_s28) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jcl31, scratch)

    scratch = (&
      screened_rates(k_cl32__p_p31__weak__wc12) + screened_rates(k_cl32__p_s31) + &
      screened_rates(k_he4_cl32__p_ar35)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cl32__ar33)*Y(jp)*state % rho - screened_rates(k_p_cl32__he4_s29) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jcl32, scratch)

    scratch = (&
      screened_rates(k_cl33__p_s32) + screened_rates(k_he4_cl33__p_ar36)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cl33__ar34)*Y(jp)*state % rho - screened_rates(k_p_cl33__he4_s30) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jcl33, scratch)

    scratch = (&
      screened_rates(k_cl34__p_s33) + screened_rates(k_he4_cl34__p_ar37)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cl34__ar35)*Y(jp)*state % rho - screened_rates(k_p_cl34__he4_s31) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jcl34, scratch)

    scratch = (&
      screened_rates(k_cl35__p_s34) + screened_rates(k_he4_cl35__p_ar38)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cl35__ar36)*Y(jp)*state % rho - screened_rates(k_p_cl35__he4_s32) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jcl35, scratch)

    scratch = (&
      screened_rates(k_cl36__p_s35) + screened_rates(k_he4_cl36__p_ar39)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cl36__ar37)*Y(jp)*state % rho - screened_rates(k_p_cl36__he4_s33) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jcl36, scratch)

    scratch = (&
      -screened_rates(k_p_cl37__ar38)*Y(jp)*state % rho - screened_rates(k_p_cl37__he4_s34)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jcl37, scratch)

    scratch = (&
      screened_rates(k_ar31__p_cl30) + 2.0d0*screened_rates(k_ar31__p_p_p29__weak__wc12) + &
      screened_rates(k_ar31__p_s30__weak__wc17) + screened_rates(k_he4_ar31__p_k34)* &
      Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp, jar31, scratch)

    scratch = (&
      screened_rates(k_ar32__p_cl31) + screened_rates(k_ar32__p_s31__weak__wc12) + &
      screened_rates(k_he4_ar32__p_k35)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ar32__he4_cl29)*Y(jp)*state % rho - screened_rates(k_p_ar32__k33) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jar32, scratch)

    scratch = (&
      screened_rates(k_ar33__p_cl32) + screened_rates(k_ar33__p_s32__weak__wc12) + &
      screened_rates(k_he4_ar33__p_k36)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ar33__he4_cl30)*Y(jp)*state % rho - screened_rates(k_p_ar33__k34) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jar33, scratch)

    scratch = (&
      screened_rates(k_ar34__p_cl33) + screened_rates(k_he4_ar34__p_k37)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ar34__he4_cl31)*Y(jp)*state % rho - screened_rates(k_p_ar34__k35) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jar34, scratch)

    scratch = (&
      screened_rates(k_ar35__p_cl34) + screened_rates(k_he4_ar35__p_k38)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ar35__he4_cl32)*Y(jp)*state % rho - screened_rates(k_p_ar35__k36) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jar35, scratch)

    scratch = (&
      screened_rates(k_ar36__p_cl35) + screened_rates(k_he4_ar36__p_k39)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ar36__he4_cl33)*Y(jp)*state % rho - screened_rates(k_p_ar36__k37) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jar36, scratch)

    scratch = (&
      screened_rates(k_ar37__p_cl36) + screened_rates(k_he4_ar37__p_k40)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ar37__he4_cl34)*Y(jp)*state % rho - screened_rates(k_p_ar37__k38) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jar37, scratch)

    scratch = (&
      screened_rates(k_ar38__p_cl37) + screened_rates(k_he4_ar38__p_k41)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ar38__he4_cl35)*Y(jp)*state % rho - screened_rates(k_p_ar38__k39) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jar38, scratch)

    scratch = (&
      -screened_rates(k_p_ar39__he4_cl36)*Y(jp)*state % rho - screened_rates(k_p_ar39__k40)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jar39, scratch)

    scratch = (&
      screened_rates(k_he4_k33__p_ca36)*Y(jhe4)*state % rho + screened_rates(k_k33__p_ar32) - &
      screened_rates(k_p_k33__ca34)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jk33, scratch)

    scratch = (&
      screened_rates(k_he4_k34__p_ca37)*Y(jhe4)*state % rho + screened_rates(k_k34__p_ar33) - &
      screened_rates(k_p_k34__ca35)*Y(jp)*state % rho - screened_rates(k_p_k34__he4_ar31)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jk34, scratch)

    scratch = (&
      screened_rates(k_he4_k35__p_ca38)*Y(jhe4)*state % rho + screened_rates(k_k35__p_ar34) + &
      screened_rates(k_k35__p_cl34__weak__wc12) - screened_rates(k_p_k35__ca36)*Y(jp)* &
      state % rho - screened_rates(k_p_k35__he4_ar32)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jk35, scratch)

    scratch = (&
      screened_rates(k_he4_k36__p_ca39)*Y(jhe4)*state % rho + screened_rates(k_k36__p_ar35) + &
      screened_rates(k_k36__p_cl35__weak__wc12) - screened_rates(k_p_k36__ca37)*Y(jp)* &
      state % rho - screened_rates(k_p_k36__he4_ar33)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jk36, scratch)

    scratch = (&
      screened_rates(k_he4_k37__p_ca40)*Y(jhe4)*state % rho + screened_rates(k_k37__p_ar36) - &
      screened_rates(k_p_k37__ca38)*Y(jp)*state % rho - screened_rates(k_p_k37__he4_ar34)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jk37, scratch)

    scratch = (&
      screened_rates(k_he4_k38__p_ca41)*Y(jhe4)*state % rho + screened_rates(k_k38__p_ar37) - &
      screened_rates(k_p_k38__ca39)*Y(jp)*state % rho - screened_rates(k_p_k38__he4_ar35)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jk38, scratch)

    scratch = (&
      screened_rates(k_he4_k39__p_ca42)*Y(jhe4)*state % rho + screened_rates(k_k39__p_ar38) - &
      screened_rates(k_p_k39__ca40)*Y(jp)*state % rho - screened_rates(k_p_k39__he4_ar36)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jk39, scratch)

    scratch = (&
      screened_rates(k_he4_k40__p_ca43)*Y(jhe4)*state % rho + screened_rates(k_k40__p_ar39) - &
      screened_rates(k_p_k40__ca41)*Y(jp)*state % rho - screened_rates(k_p_k40__he4_ar37)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jk40, scratch)

    scratch = (&
      screened_rates(k_he4_k41__p_ca44)*Y(jhe4)*state % rho - screened_rates(k_p_k41__ca42)* &
      Y(jp)*state % rho - screened_rates(k_p_k41__he4_ar38)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jk41, scratch)

    scratch = (&
      screened_rates(k_ca34__p_k33) + screened_rates(k_he4_ca34__p_sc37)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp, jca34, scratch)

    scratch = (&
      screened_rates(k_ca35__p_ar34__weak__wc17) + screened_rates(k_ca35__p_k34) + 2.0d0* &
      screened_rates(k_ca35__p_p_cl33__weak__wc12) + screened_rates(k_he4_ca35__p_sc38)* &
      Y(jhe4)*state % rho - screened_rates(k_p_ca35__sc36)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jca35, scratch)

    scratch = (&
      screened_rates(k_ca36__p_ar35__weak__wc12) + screened_rates(k_ca36__p_k35) + 2.0d0* &
      screened_rates(k_he4_ca36__p_p_ca38)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_ca36__p_sc39)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ca36__he4_k33)*Y(jp)*state % rho - screened_rates(k_p_ca36__sc37) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jca36, scratch)

    scratch = (&
      screened_rates(k_ca37__p_ar36__weak__wc12) + screened_rates(k_ca37__p_k36) + &
      screened_rates(k_he4_ca37__p_sc40)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ca37__he4_k34)*Y(jp)*state % rho - screened_rates(k_p_ca37__sc38) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jca37, scratch)

    scratch = (&
      screened_rates(k_ca38__p_k37) + screened_rates(k_he4_ca38__p_sc41)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ca38__he4_k35)*Y(jp)*state % rho - screened_rates(k_p_ca38__sc39) &
      *Y(jp)*state % rho - screened_rates(k_p_p_ca38__he4_ca36)*Y(jp)**2* &
      state % rho**2 &
       )
    call set_jac_entry(state, jp, jca38, scratch)

    scratch = (&
      screened_rates(k_ca39__p_k38) + screened_rates(k_he4_ca39__p_sc42)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ca39__he4_k36)*Y(jp)*state % rho - screened_rates(k_p_ca39__sc40) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jca39, scratch)

    scratch = (&
      screened_rates(k_ca40__p_k39) + screened_rates(k_he4_ca40__p_sc43)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ca40__he4_k37)*Y(jp)*state % rho - screened_rates(k_p_ca40__sc41) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jca40, scratch)

    scratch = (&
      screened_rates(k_ca41__p_k40) + screened_rates(k_he4_ca41__p_sc44)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ca41__he4_k38)*Y(jp)*state % rho - screened_rates(k_p_ca41__sc42) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jca41, scratch)

    scratch = (&
      screened_rates(k_ca42__p_k41) + screened_rates(k_he4_ca42__p_sc45)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ca42__he4_k39)*Y(jp)*state % rho - screened_rates(k_p_ca42__sc43) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jca42, scratch)

    scratch = (&
      screened_rates(k_he4_ca43__p_sc46)*Y(jhe4)*state % rho - screened_rates(k_p_ca43__he4_k40)* &
      Y(jp)*state % rho - screened_rates(k_p_ca43__sc44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jca43, scratch)

    scratch = (&
      -screened_rates(k_p_ca44__he4_k41)*Y(jp)*state % rho - screened_rates(k_p_ca44__sc45)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jca44, scratch)

    scratch = (&
      screened_rates(k_he4_sc36__p_ti39)*Y(jhe4)*state % rho + screened_rates(k_sc36__p_ca35) &
       )
    call set_jac_entry(state, jp, jsc36, scratch)

    scratch = (&
      screened_rates(k_he4_sc37__p_ti40)*Y(jhe4)*state % rho - screened_rates(k_p_sc37__he4_ca34) &
      *Y(jp)*state % rho - screened_rates(k_p_sc37__ti38)*Y(jp)*state % rho + &
      screened_rates(k_sc37__p_ca36) &
       )
    call set_jac_entry(state, jp, jsc37, scratch)

    scratch = (&
      screened_rates(k_he4_sc38__p_ti41)*Y(jhe4)*state % rho - screened_rates(k_p_sc38__he4_ca35) &
      *Y(jp)*state % rho - screened_rates(k_p_sc38__ti39)*Y(jp)*state % rho + &
      screened_rates(k_sc38__p_ca37) &
       )
    call set_jac_entry(state, jp, jsc38, scratch)

    scratch = (&
      screened_rates(k_he4_sc39__p_ti42)*Y(jhe4)*state % rho - screened_rates(k_p_sc39__he4_ca36) &
      *Y(jp)*state % rho - screened_rates(k_p_sc39__ti40)*Y(jp)*state % rho + &
      screened_rates(k_sc39__p_ca38) &
       )
    call set_jac_entry(state, jp, jsc39, scratch)

    scratch = (&
      screened_rates(k_he4_sc40__p_ti43)*Y(jhe4)*state % rho - screened_rates(k_p_sc40__he4_ca37) &
      *Y(jp)*state % rho - screened_rates(k_p_sc40__ti41)*Y(jp)*state % rho + &
      screened_rates(k_sc40__p_ca39) + screened_rates(k_sc40__p_k39__weak__wc12) &
       )
    call set_jac_entry(state, jp, jsc40, scratch)

    scratch = (&
      screened_rates(k_he4_sc41__p_ti44)*Y(jhe4)*state % rho - screened_rates(k_p_sc41__he4_ca38) &
      *Y(jp)*state % rho - screened_rates(k_p_sc41__ti42)*Y(jp)*state % rho + &
      screened_rates(k_sc41__p_ca40) &
       )
    call set_jac_entry(state, jp, jsc41, scratch)

    scratch = (&
      screened_rates(k_he4_sc42__p_ti45)*Y(jhe4)*state % rho - screened_rates(k_p_sc42__he4_ca39) &
      *Y(jp)*state % rho - screened_rates(k_p_sc42__ti43)*Y(jp)*state % rho + &
      screened_rates(k_sc42__p_ca41) &
       )
    call set_jac_entry(state, jp, jsc42, scratch)

    scratch = (&
      screened_rates(k_he4_sc43__p_ti46)*Y(jhe4)*state % rho - screened_rates(k_p_sc43__he4_ca40) &
      *Y(jp)*state % rho - screened_rates(k_p_sc43__ti44)*Y(jp)*state % rho + &
      screened_rates(k_sc43__p_ca42) &
       )
    call set_jac_entry(state, jp, jsc43, scratch)

    scratch = (&
      screened_rates(k_he4_sc44__p_ti47)*Y(jhe4)*state % rho - screened_rates(k_p_sc44__he4_ca41) &
      *Y(jp)*state % rho - screened_rates(k_p_sc44__ti45)*Y(jp)*state % rho + &
      screened_rates(k_sc44__p_ca43) &
       )
    call set_jac_entry(state, jp, jsc44, scratch)

    scratch = (&
      screened_rates(k_he4_sc45__p_ti48)*Y(jhe4)*state % rho - screened_rates(k_p_sc45__he4_ca42) &
      *Y(jp)*state % rho - screened_rates(k_p_sc45__ti46)*Y(jp)*state % rho + &
      screened_rates(k_sc45__p_ca44) &
       )
    call set_jac_entry(state, jp, jsc45, scratch)

    scratch = (&
      -screened_rates(k_p_sc46__he4_ca43)*Y(jp)*state % rho - screened_rates(k_p_sc46__ti47)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jsc46, scratch)

    scratch = (&
      screened_rates(k_he4_ti38__p_v41)*Y(jhe4)*state % rho + screened_rates(k_ti38__p_sc37) &
       )
    call set_jac_entry(state, jp, jti38, scratch)

    scratch = (&
      screened_rates(k_he4_ti39__p_v42)*Y(jhe4)*state % rho - screened_rates(k_p_ti39__he4_sc36)* &
      Y(jp)*state % rho - screened_rates(k_p_ti39__v40)*Y(jp)*state % rho + &
      screened_rates(k_ti39__p_ca38__weak__wc12) + screened_rates(k_ti39__p_sc38) &
       )
    call set_jac_entry(state, jp, jti39, scratch)

    scratch = (&
      screened_rates(k_he4_ti40__p_v43)*Y(jhe4)*state % rho - screened_rates(k_p_ti40__he4_sc37)* &
      Y(jp)*state % rho - screened_rates(k_p_ti40__v41)*Y(jp)*state % rho + &
      screened_rates(k_ti40__p_ca39__weak__wc12) + screened_rates(k_ti40__p_sc39) &
       )
    call set_jac_entry(state, jp, jti40, scratch)

    scratch = (&
      screened_rates(k_he4_ti41__p_v44)*Y(jhe4)*state % rho - screened_rates(k_p_ti41__he4_sc38)* &
      Y(jp)*state % rho - screened_rates(k_p_ti41__v42)*Y(jp)*state % rho + &
      screened_rates(k_ti41__p_ca40__weak__wc12) + screened_rates(k_ti41__p_sc40) &
       )
    call set_jac_entry(state, jp, jti41, scratch)

    scratch = (&
      screened_rates(k_he4_ti42__p_v45)*Y(jhe4)*state % rho - screened_rates(k_p_ti42__he4_sc39)* &
      Y(jp)*state % rho - screened_rates(k_p_ti42__v43)*Y(jp)*state % rho + &
      screened_rates(k_ti42__p_sc41) &
       )
    call set_jac_entry(state, jp, jti42, scratch)

    scratch = (&
      screened_rates(k_he4_ti43__p_v46)*Y(jhe4)*state % rho - screened_rates(k_p_ti43__he4_sc40)* &
      Y(jp)*state % rho - screened_rates(k_p_ti43__v44)*Y(jp)*state % rho + &
      screened_rates(k_ti43__p_sc42) &
       )
    call set_jac_entry(state, jp, jti43, scratch)

    scratch = (&
      screened_rates(k_he4_ti44__p_v47)*Y(jhe4)*state % rho - screened_rates(k_p_ti44__he4_sc41)* &
      Y(jp)*state % rho - screened_rates(k_p_ti44__v45)*Y(jp)*state % rho + &
      screened_rates(k_ti44__p_sc43) &
       )
    call set_jac_entry(state, jp, jti44, scratch)

    scratch = (&
      screened_rates(k_he4_ti45__p_v48)*Y(jhe4)*state % rho - screened_rates(k_p_ti45__he4_sc42)* &
      Y(jp)*state % rho - screened_rates(k_p_ti45__v46)*Y(jp)*state % rho + &
      screened_rates(k_ti45__p_sc44) &
       )
    call set_jac_entry(state, jp, jti45, scratch)

    scratch = (&
      screened_rates(k_he4_ti46__p_v49)*Y(jhe4)*state % rho - screened_rates(k_p_ti46__he4_sc43)* &
      Y(jp)*state % rho - screened_rates(k_p_ti46__v47)*Y(jp)*state % rho + &
      screened_rates(k_ti46__p_sc45) &
       )
    call set_jac_entry(state, jp, jti46, scratch)

    scratch = (&
      screened_rates(k_he4_ti47__p_v50)*Y(jhe4)*state % rho - screened_rates(k_p_ti47__he4_sc44)* &
      Y(jp)*state % rho - screened_rates(k_p_ti47__v48)*Y(jp)*state % rho + &
      screened_rates(k_ti47__p_sc46) &
       )
    call set_jac_entry(state, jp, jti47, scratch)

    scratch = (&
      -screened_rates(k_p_ti48__he4_sc45)*Y(jp)*state % rho - screened_rates(k_p_ti48__v49)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jti48, scratch)

    scratch = (&
      screened_rates(k_he4_v40__p_cr43)*Y(jhe4)*state % rho + screened_rates(k_v40__p_ti39) &
       )
    call set_jac_entry(state, jp, jv40, scratch)

    scratch = (&
      screened_rates(k_he4_v41__p_cr44)*Y(jhe4)*state % rho - screened_rates(k_p_v41__cr42)* &
      Y(jp)*state % rho - screened_rates(k_p_v41__he4_ti38)*Y(jp)*state % rho + &
      screened_rates(k_v41__p_ti40) &
       )
    call set_jac_entry(state, jp, jv41, scratch)

    scratch = (&
      screened_rates(k_he4_v42__p_cr45)*Y(jhe4)*state % rho - screened_rates(k_p_v42__cr43)* &
      Y(jp)*state % rho - screened_rates(k_p_v42__he4_ti39)*Y(jp)*state % rho + &
      screened_rates(k_v42__p_ti41) &
       )
    call set_jac_entry(state, jp, jv42, scratch)

    scratch = (&
      screened_rates(k_he4_v43__p_cr46)*Y(jhe4)*state % rho - screened_rates(k_p_v43__cr44)* &
      Y(jp)*state % rho - screened_rates(k_p_v43__he4_ti40)*Y(jp)*state % rho + &
      screened_rates(k_v43__p_ti42) &
       )
    call set_jac_entry(state, jp, jv43, scratch)

    scratch = (&
      screened_rates(k_he4_v44__p_cr47)*Y(jhe4)*state % rho - screened_rates(k_p_v44__cr45)* &
      Y(jp)*state % rho - screened_rates(k_p_v44__he4_ti41)*Y(jp)*state % rho + &
      screened_rates(k_v44__p_ti43) &
       )
    call set_jac_entry(state, jp, jv44, scratch)

    scratch = (&
      screened_rates(k_he4_v45__p_cr48)*Y(jhe4)*state % rho - screened_rates(k_p_v45__cr46)* &
      Y(jp)*state % rho - screened_rates(k_p_v45__he4_ti42)*Y(jp)*state % rho + &
      screened_rates(k_v45__p_ti44) &
       )
    call set_jac_entry(state, jp, jv45, scratch)

    scratch = (&
      screened_rates(k_he4_v46__p_cr49)*Y(jhe4)*state % rho - screened_rates(k_p_v46__cr47)* &
      Y(jp)*state % rho - screened_rates(k_p_v46__he4_ti43)*Y(jp)*state % rho + &
      screened_rates(k_v46__p_ti45) &
       )
    call set_jac_entry(state, jp, jv46, scratch)

    scratch = (&
      screened_rates(k_he4_v47__p_cr50)*Y(jhe4)*state % rho - screened_rates(k_p_v47__cr48)* &
      Y(jp)*state % rho - screened_rates(k_p_v47__he4_ti44)*Y(jp)*state % rho + &
      screened_rates(k_v47__p_ti46) &
       )
    call set_jac_entry(state, jp, jv47, scratch)

    scratch = (&
      screened_rates(k_he4_v48__p_cr51)*Y(jhe4)*state % rho - screened_rates(k_p_v48__cr49)* &
      Y(jp)*state % rho - screened_rates(k_p_v48__he4_ti45)*Y(jp)*state % rho + &
      screened_rates(k_v48__p_ti47) &
       )
    call set_jac_entry(state, jp, jv48, scratch)

    scratch = (&
      screened_rates(k_he4_v49__p_cr52)*Y(jhe4)*state % rho - screened_rates(k_p_v49__cr50)* &
      Y(jp)*state % rho - screened_rates(k_p_v49__he4_ti46)*Y(jp)*state % rho + &
      screened_rates(k_v49__p_ti48) &
       )
    call set_jac_entry(state, jp, jv49, scratch)

    scratch = (&
      -screened_rates(k_p_v50__cr51)*Y(jp)*state % rho - screened_rates(k_p_v50__he4_ti47)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jv50, scratch)

    scratch = (&
      screened_rates(k_cr42__p_ti41__weak__wc12) + screened_rates(k_cr42__p_v41) + &
      screened_rates(k_he4_cr42__p_mn45)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp, jcr42, scratch)

    scratch = (&
      3.0d0*screened_rates(k_cr43__p_p_p_ca40__weak__wc12) + 2.0d0* &
      screened_rates(k_cr43__p_p_sc41__weak__wc12) + &
      screened_rates(k_cr43__p_ti42__weak__wc17) + screened_rates(k_cr43__p_v42) + &
      screened_rates(k_he4_cr43__p_mn46)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cr43__he4_v40)*Y(jp)*state % rho - screened_rates(k_p_cr43__mn44) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jcr43, scratch)

    scratch = (&
      screened_rates(k_cr44__p_ti43__weak__wc12) + screened_rates(k_cr44__p_v43) + &
      screened_rates(k_he4_cr44__p_mn47)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cr44__he4_v41)*Y(jp)*state % rho - screened_rates(k_p_cr44__mn45) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jcr44, scratch)

    scratch = (&
      screened_rates(k_cr45__p_ti44__weak__wc12) + screened_rates(k_cr45__p_v44) + &
      screened_rates(k_he4_cr45__p_mn48)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cr45__he4_v42)*Y(jp)*state % rho - screened_rates(k_p_cr45__mn46) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jcr45, scratch)

    scratch = (&
      screened_rates(k_cr46__p_v45) + screened_rates(k_he4_cr46__p_mn49)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cr46__he4_v43)*Y(jp)*state % rho - screened_rates(k_p_cr46__mn47) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jcr46, scratch)

    scratch = (&
      screened_rates(k_cr47__p_v46) + screened_rates(k_he4_cr47__p_mn50)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cr47__he4_v44)*Y(jp)*state % rho - screened_rates(k_p_cr47__mn48) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jcr47, scratch)

    scratch = (&
      screened_rates(k_cr48__p_v47) + screened_rates(k_he4_cr48__p_mn51)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cr48__he4_v45)*Y(jp)*state % rho - screened_rates(k_p_cr48__mn49) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jcr48, scratch)

    scratch = (&
      screened_rates(k_cr49__p_v48) + screened_rates(k_he4_cr49__p_mn52)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cr49__he4_v46)*Y(jp)*state % rho - screened_rates(k_p_cr49__mn50) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jcr49, scratch)

    scratch = (&
      screened_rates(k_cr50__p_v49) + screened_rates(k_he4_cr50__p_mn53)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cr50__he4_v47)*Y(jp)*state % rho - screened_rates(k_p_cr50__mn51) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jcr50, scratch)

    scratch = (&
      screened_rates(k_cr51__p_v50) + screened_rates(k_he4_cr51__p_mn54)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cr51__he4_v48)*Y(jp)*state % rho - screened_rates(k_p_cr51__mn52) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jcr51, scratch)

    scratch = (&
      screened_rates(k_he4_cr52__p_mn55)*Y(jhe4)*state % rho - screened_rates(k_p_cr52__he4_v49)* &
      Y(jp)*state % rho - screened_rates(k_p_cr52__mn53)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jcr52, scratch)

    scratch = (&
      screened_rates(k_he4_mn44__p_fe47)*Y(jhe4)*state % rho + screened_rates(k_mn44__p_cr43) - &
      screened_rates(k_p_mn44__fe45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jmn44, scratch)

    scratch = (&
      screened_rates(k_he4_mn45__p_fe48)*Y(jhe4)*state % rho + screened_rates(k_mn45__p_cr44) - &
      screened_rates(k_p_mn45__fe46)*Y(jp)*state % rho - &
      screened_rates(k_p_mn45__he4_cr42)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jmn45, scratch)

    scratch = (&
      screened_rates(k_he4_mn46__p_fe49)*Y(jhe4)*state % rho + screened_rates(k_mn46__p_cr45) + &
      screened_rates(k_mn46__p_v45__weak__wc12) - screened_rates(k_p_mn46__fe47)*Y(jp)* &
      state % rho - screened_rates(k_p_mn46__he4_cr43)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jmn46, scratch)

    scratch = (&
      screened_rates(k_he4_mn47__p_fe50)*Y(jhe4)*state % rho + screened_rates(k_mn47__p_cr46) - &
      screened_rates(k_p_mn47__fe48)*Y(jp)*state % rho - &
      screened_rates(k_p_mn47__he4_cr44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jmn47, scratch)

    scratch = (&
      screened_rates(k_he4_mn48__p_fe51)*Y(jhe4)*state % rho + screened_rates(k_mn48__p_cr47) + &
      screened_rates(k_mn48__p_v47__weak__wc12) - screened_rates(k_p_mn48__fe49)*Y(jp)* &
      state % rho - screened_rates(k_p_mn48__he4_cr45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jmn48, scratch)

    scratch = (&
      screened_rates(k_he4_mn49__p_fe52)*Y(jhe4)*state % rho + screened_rates(k_mn49__p_cr48) - &
      screened_rates(k_p_mn49__fe50)*Y(jp)*state % rho - &
      screened_rates(k_p_mn49__he4_cr46)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jmn49, scratch)

    scratch = (&
      screened_rates(k_he4_mn50__p_fe53)*Y(jhe4)*state % rho + screened_rates(k_mn50__p_cr49) - &
      screened_rates(k_p_mn50__fe51)*Y(jp)*state % rho - &
      screened_rates(k_p_mn50__he4_cr47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jmn50, scratch)

    scratch = (&
      screened_rates(k_he4_mn51__p_fe54)*Y(jhe4)*state % rho + screened_rates(k_mn51__p_cr50) - &
      screened_rates(k_p_mn51__fe52)*Y(jp)*state % rho - &
      screened_rates(k_p_mn51__he4_cr48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jmn51, scratch)

    scratch = (&
      screened_rates(k_he4_mn52__p_fe55)*Y(jhe4)*state % rho + screened_rates(k_mn52__p_cr51) - &
      screened_rates(k_p_mn52__fe53)*Y(jp)*state % rho - &
      screened_rates(k_p_mn52__he4_cr49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jmn52, scratch)

    scratch = (&
      screened_rates(k_he4_mn53__p_fe56)*Y(jhe4)*state % rho + screened_rates(k_mn53__p_cr52) - &
      screened_rates(k_p_mn53__fe54)*Y(jp)*state % rho - &
      screened_rates(k_p_mn53__he4_cr50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jmn53, scratch)

    scratch = (&
      -screened_rates(k_p_mn54__fe55)*Y(jp)*state % rho - screened_rates(k_p_mn54__he4_cr51)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jmn54, scratch)

    scratch = (&
      -screened_rates(k_p_mn55__fe56)*Y(jp)*state % rho - screened_rates(k_p_mn55__he4_cr52)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jmn55, scratch)

    scratch = (&
      screened_rates(k_fe45__p_cr44__weak__wc17) + screened_rates(k_fe45__p_mn44) + 3.0d0* &
      screened_rates(k_fe45__p_p_p_ti42__weak__wc12) + 2.0d0* &
      screened_rates(k_fe45__p_p_v43__weak__wc12) + screened_rates(k_he4_fe45__p_co48)* &
      Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp, jfe45, scratch)

    scratch = (&
      screened_rates(k_fe46__p_cr45__weak__wc12) + screened_rates(k_fe46__p_mn45) + &
      screened_rates(k_he4_fe46__p_co49)*Y(jhe4)*state % rho - &
      screened_rates(k_p_fe46__co47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jfe46, scratch)

    scratch = (&
      screened_rates(k_fe47__p_cr46__weak__wc12) + screened_rates(k_fe47__p_mn46) + &
      screened_rates(k_he4_fe47__p_co50)*Y(jhe4)*state % rho - &
      screened_rates(k_p_fe47__co48)*Y(jp)*state % rho - &
      screened_rates(k_p_fe47__he4_mn44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jfe47, scratch)

    scratch = (&
      screened_rates(k_fe48__p_cr47__weak__wc12) + screened_rates(k_fe48__p_mn47) + &
      screened_rates(k_he4_fe48__p_co51)*Y(jhe4)*state % rho - &
      screened_rates(k_p_fe48__co49)*Y(jp)*state % rho - &
      screened_rates(k_p_fe48__he4_mn45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jfe48, scratch)

    scratch = (&
      screened_rates(k_fe49__p_cr48__weak__wc12) + screened_rates(k_fe49__p_mn48) + &
      screened_rates(k_he4_fe49__p_co52)*Y(jhe4)*state % rho - &
      screened_rates(k_p_fe49__co50)*Y(jp)*state % rho - &
      screened_rates(k_p_fe49__he4_mn46)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jfe49, scratch)

    scratch = (&
      screened_rates(k_fe50__p_mn49) + screened_rates(k_he4_fe50__p_co53)*Y(jhe4)*state % rho - &
      screened_rates(k_p_fe50__co51)*Y(jp)*state % rho - &
      screened_rates(k_p_fe50__he4_mn47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jfe50, scratch)

    scratch = (&
      screened_rates(k_fe51__p_mn50) + screened_rates(k_he4_fe51__p_co54)*Y(jhe4)*state % rho - &
      screened_rates(k_p_fe51__co52)*Y(jp)*state % rho - &
      screened_rates(k_p_fe51__he4_mn48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jfe51, scratch)

    scratch = (&
      screened_rates(k_fe52__p_mn51) + screened_rates(k_he4_fe52__p_co55)*Y(jhe4)*state % rho - &
      screened_rates(k_p_fe52__co53)*Y(jp)*state % rho - &
      screened_rates(k_p_fe52__he4_mn49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jfe52, scratch)

    scratch = (&
      screened_rates(k_fe53__p_mn52) + screened_rates(k_he4_fe53__p_co56)*Y(jhe4)*state % rho - &
      screened_rates(k_p_fe53__co54)*Y(jp)*state % rho - &
      screened_rates(k_p_fe53__he4_mn50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jfe53, scratch)

    scratch = (&
      screened_rates(k_fe54__p_mn53) - screened_rates(k_p_fe54__co55)*Y(jp)*state % rho - &
      screened_rates(k_p_fe54__he4_mn51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jfe54, scratch)

    scratch = (&
      screened_rates(k_fe55__p_mn54) - screened_rates(k_p_fe55__co56)*Y(jp)*state % rho - &
      screened_rates(k_p_fe55__he4_mn52)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jfe55, scratch)

    scratch = (&
      screened_rates(k_fe56__p_mn55) - screened_rates(k_p_fe56__he4_mn53)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jfe56, scratch)

    scratch = (&
      screened_rates(k_co47__p_fe46) + screened_rates(k_he4_co47__p_ni50)*Y(jhe4)*state % rho - &
      screened_rates(k_p_co47__ni48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jco47, scratch)

    scratch = (&
      screened_rates(k_co48__p_fe47) + screened_rates(k_he4_co48__p_ni51)*Y(jhe4)*state % rho - &
      screened_rates(k_p_co48__he4_fe45)*Y(jp)*state % rho - &
      screened_rates(k_p_co48__ni49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jco48, scratch)

    scratch = (&
      screened_rates(k_co49__p_fe48) + screened_rates(k_he4_co49__p_ni52)*Y(jhe4)*state % rho - &
      screened_rates(k_p_co49__he4_fe46)*Y(jp)*state % rho - &
      screened_rates(k_p_co49__ni50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jco49, scratch)

    scratch = (&
      screened_rates(k_co50__p_fe49) + screened_rates(k_co50__p_mn49__weak__wc12) + &
      screened_rates(k_he4_co50__p_ni53)*Y(jhe4)*state % rho - &
      screened_rates(k_p_co50__he4_fe47)*Y(jp)*state % rho - &
      screened_rates(k_p_co50__ni51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jco50, scratch)

    scratch = (&
      screened_rates(k_co51__p_fe50) + screened_rates(k_he4_co51__p_ni54)*Y(jhe4)*state % rho - &
      screened_rates(k_p_co51__he4_fe48)*Y(jp)*state % rho - &
      screened_rates(k_p_co51__ni52)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jco51, scratch)

    scratch = (&
      screened_rates(k_co52__p_fe51) + screened_rates(k_he4_co52__p_ni55)*Y(jhe4)*state % rho - &
      screened_rates(k_p_co52__he4_fe49)*Y(jp)*state % rho - &
      screened_rates(k_p_co52__ni53)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jco52, scratch)

    scratch = (&
      screened_rates(k_co53__p_fe52) + screened_rates(k_he4_co53__p_ni56)*Y(jhe4)*state % rho - &
      screened_rates(k_p_co53__he4_fe50)*Y(jp)*state % rho - &
      screened_rates(k_p_co53__ni54)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jco53, scratch)

    scratch = (&
      screened_rates(k_co54__p_fe53) - screened_rates(k_p_co54__he4_fe51)*Y(jp)*state % rho - &
      screened_rates(k_p_co54__ni55)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jco54, scratch)

    scratch = (&
      screened_rates(k_co55__p_fe54) - screened_rates(k_p_co55__he4_fe52)*Y(jp)*state % rho - &
      screened_rates(k_p_co55__ni56)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jco55, scratch)

    scratch = (&
      screened_rates(k_co56__p_fe55) - screened_rates(k_p_co56__he4_fe53)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jco56, scratch)

    scratch = (&
      screened_rates(k_ni48__p_co47) &
       )
    call set_jac_entry(state, jp, jni48, scratch)

    scratch = (&
      screened_rates(k_ni49__p_co48) + screened_rates(k_ni49__p_fe48__weak__wc12) &
       )
    call set_jac_entry(state, jp, jni49, scratch)

    scratch = (&
      screened_rates(k_ni50__p_co49) + screened_rates(k_ni50__p_fe49__weak__wc12) - &
      screened_rates(k_p_ni50__he4_co47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jni50, scratch)

    scratch = (&
      screened_rates(k_ni51__p_co50) + screened_rates(k_ni51__p_fe50__weak__wc12) - &
      screened_rates(k_p_ni51__he4_co48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jni51, scratch)

    scratch = (&
      screened_rates(k_ni52__p_co51) + screened_rates(k_ni52__p_fe51__weak__wc12) - &
      screened_rates(k_p_ni52__he4_co49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jni52, scratch)

    scratch = (&
      screened_rates(k_ni53__p_co52) + screened_rates(k_ni53__p_fe52__weak__wc12) - &
      screened_rates(k_p_ni53__he4_co50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jni53, scratch)

    scratch = (&
      screened_rates(k_ni54__p_co53) - screened_rates(k_p_ni54__he4_co51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jni54, scratch)

    scratch = (&
      screened_rates(k_ni55__p_co54) - screened_rates(k_p_ni55__he4_co52)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jni55, scratch)

    scratch = (&
      screened_rates(k_ni56__p_co55) - screened_rates(k_p_ni56__he4_co53)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp, jni56, scratch)

    scratch = (&
      -screened_rates(k_p_d__he3)*Y(jd)*state % rho + screened_rates(k_p_he4__d_he3)*Y(jhe4)* &
      state % rho + 0.5d0*screened_rates(k_p_he4_he4__d_be7)*Y(jhe4)**2*state % rho &
      **2 + 1.0d0*screened_rates(k_p_p__d__weak__bet_pos_)*Y(jp)*state % rho + &
      1.0d0*screened_rates(k_p_p__d__weak__electron_capture)*Y(jp)*state % rho**2 &
      *state % y_e &
       )
    call set_jac_entry(state, jd, jp, scratch)

    scratch = (&
      -screened_rates(k_d_be7__p_he4_he4)*Y(jbe7)*state % rho - 2.0d0*screened_rates(k_d_d__he4)* &
      Y(jd)*state % rho - screened_rates(k_d_he3__p_he4)*Y(jhe3)*state % rho - &
      screened_rates(k_p_d__he3)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jd, jd, scratch)

    scratch = (&
      -screened_rates(k_d_he3__p_he4)*Y(jd)*state % rho + screened_rates(k_he3__p_d) &
       )
    call set_jac_entry(state, jd, jhe3, scratch)

    scratch = (&
      2.0d0*screened_rates(k_he4__d_d) + screened_rates(k_p_he4__d_he3)*Y(jp)*state % rho + 1.0d0 &
      *screened_rates(k_p_he4_he4__d_be7)*Y(jhe4)*Y(jp)*state % rho**2 &
       )
    call set_jac_entry(state, jd, jhe4, scratch)

    scratch = (&
      -screened_rates(k_d_be7__p_he4_he4)*Y(jd)*state % rho &
       )
    call set_jac_entry(state, jd, jbe7, scratch)

    scratch = (&
      screened_rates(k_p_d__he3)*Y(jd)*state % rho - screened_rates(k_p_he3__he4__weak__bet_pos_) &
      *Y(jhe3)*state % rho + screened_rates(k_p_he4__d_he3)*Y(jhe4)*state % rho + &
      2.0d0*screened_rates(k_p_p_he4__he3_he3)*Y(jhe4)*Y(jp)*state % rho**2 + &
      0.5d0*screened_rates(k_p_p_he4_he4__he3_be7)*Y(jhe4)**2*Y(jp)*state % rho &
      **3 &
       )
    call set_jac_entry(state, jhe3, jp, scratch)

    scratch = (&
      -screened_rates(k_d_he3__p_he4)*Y(jhe3)*state % rho + screened_rates(k_p_d__he3)*Y(jp)* &
      state % rho &
       )
    call set_jac_entry(state, jhe3, jd, scratch)

    scratch = (&
      -screened_rates(k_d_he3__p_he4)*Y(jd)*state % rho - screened_rates(k_he3__p_d) - &
      screened_rates(k_he3_be7__p_p_he4_he4)*Y(jbe7)*state % rho - 2.0d0* &
      screened_rates(k_he3_he3__p_p_he4)*Y(jhe3)*state % rho - &
      screened_rates(k_he4_he3__be7)*Y(jhe4)*state % rho - &
      screened_rates(k_p_he3__he4__weak__bet_pos_)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe3, jhe3, scratch)

    scratch = (&
      -screened_rates(k_he4_he3__be7)*Y(jhe3)*state % rho + screened_rates(k_p_he4__d_he3)* &
      Y(jp)*state % rho + screened_rates(k_p_p_he4__he3_he3)*Y(jp)**2*state % rho &
      **2 + 0.5d0*screened_rates(k_p_p_he4_he4__he3_be7)*Y(jhe4)*Y(jp)**2* &
      state % rho**3 &
       )
    call set_jac_entry(state, jhe3, jhe4, scratch)

    scratch = (&
      screened_rates(k_be7__he4_he3) - screened_rates(k_he3_be7__p_p_he4_he4)*Y(jhe3)*state % rho &
       )
    call set_jac_entry(state, jhe3, jbe7, scratch)

    scratch = (&
      screened_rates(k_p_al24__he4_mg21)*Y(jal24)*state % rho + &
      screened_rates(k_p_al25__he4_mg22)*Y(jal25)*state % rho + &
      screened_rates(k_p_al26__he4_mg23)*Y(jal26)*state % rho + &
      screened_rates(k_p_al27__he4_mg24)*Y(jal27)*state % rho + &
      screened_rates(k_p_al28__he4_mg25)*Y(jal28)*state % rho + &
      screened_rates(k_p_ar32__he4_cl29)*Y(jar32)*state % rho + &
      screened_rates(k_p_ar33__he4_cl30)*Y(jar33)*state % rho + &
      screened_rates(k_p_ar34__he4_cl31)*Y(jar34)*state % rho + &
      screened_rates(k_p_ar35__he4_cl32)*Y(jar35)*state % rho + &
      screened_rates(k_p_ar36__he4_cl33)*Y(jar36)*state % rho + &
      screened_rates(k_p_ar37__he4_cl34)*Y(jar37)*state % rho + &
      screened_rates(k_p_ar38__he4_cl35)*Y(jar38)*state % rho + &
      screened_rates(k_p_ar39__he4_cl36)*Y(jar39)*state % rho + &
      screened_rates(k_p_ca36__he4_k33)*Y(jca36)*state % rho + &
      screened_rates(k_p_ca37__he4_k34)*Y(jca37)*state % rho + &
      screened_rates(k_p_ca38__he4_k35)*Y(jca38)*state % rho + &
      screened_rates(k_p_ca39__he4_k36)*Y(jca39)*state % rho + &
      screened_rates(k_p_ca40__he4_k37)*Y(jca40)*state % rho + &
      screened_rates(k_p_ca41__he4_k38)*Y(jca41)*state % rho + &
      screened_rates(k_p_ca42__he4_k39)*Y(jca42)*state % rho + &
      screened_rates(k_p_ca43__he4_k40)*Y(jca43)*state % rho + &
      screened_rates(k_p_ca44__he4_k41)*Y(jca44)*state % rho + &
      screened_rates(k_p_cl31__he4_s28)*Y(jcl31)*state % rho + &
      screened_rates(k_p_cl32__he4_s29)*Y(jcl32)*state % rho + &
      screened_rates(k_p_cl33__he4_s30)*Y(jcl33)*state % rho + &
      screened_rates(k_p_cl34__he4_s31)*Y(jcl34)*state % rho + &
      screened_rates(k_p_cl35__he4_s32)*Y(jcl35)*state % rho + &
      screened_rates(k_p_cl36__he4_s33)*Y(jcl36)*state % rho + &
      screened_rates(k_p_cl37__he4_s34)*Y(jcl37)*state % rho + &
      screened_rates(k_p_co48__he4_fe45)*Y(jco48)*state % rho + &
      screened_rates(k_p_co49__he4_fe46)*Y(jco49)*state % rho + &
      screened_rates(k_p_co50__he4_fe47)*Y(jco50)*state % rho + &
      screened_rates(k_p_co51__he4_fe48)*Y(jco51)*state % rho + &
      screened_rates(k_p_co52__he4_fe49)*Y(jco52)*state % rho + &
      screened_rates(k_p_co53__he4_fe50)*Y(jco53)*state % rho + &
      screened_rates(k_p_co54__he4_fe51)*Y(jco54)*state % rho + &
      screened_rates(k_p_co55__he4_fe52)*Y(jco55)*state % rho + &
      screened_rates(k_p_co56__he4_fe53)*Y(jco56)*state % rho + &
      screened_rates(k_p_cr43__he4_v40)*Y(jcr43)*state % rho + &
      screened_rates(k_p_cr44__he4_v41)*Y(jcr44)*state % rho + &
      screened_rates(k_p_cr45__he4_v42)*Y(jcr45)*state % rho + &
      screened_rates(k_p_cr46__he4_v43)*Y(jcr46)*state % rho + &
      screened_rates(k_p_cr47__he4_v44)*Y(jcr47)*state % rho + &
      screened_rates(k_p_cr48__he4_v45)*Y(jcr48)*state % rho + &
      screened_rates(k_p_cr49__he4_v46)*Y(jcr49)*state % rho + &
      screened_rates(k_p_cr50__he4_v47)*Y(jcr50)*state % rho + &
      screened_rates(k_p_cr51__he4_v48)*Y(jcr51)*state % rho + &
      screened_rates(k_p_cr52__he4_v49)*Y(jcr52)*state % rho + &
      screened_rates(k_p_f17__he4_o14)*Y(jf17)*state % rho + &
      screened_rates(k_p_f18__he4_o15)*Y(jf18)*state % rho + &
      screened_rates(k_p_f19__he4_o16)*Y(jf19)*state % rho + &
      screened_rates(k_p_f20__he4_o17)*Y(jf20)*state % rho + &
      screened_rates(k_p_fe47__he4_mn44)*Y(jfe47)*state % rho + &
      screened_rates(k_p_fe48__he4_mn45)*Y(jfe48)*state % rho + &
      screened_rates(k_p_fe49__he4_mn46)*Y(jfe49)*state % rho + &
      screened_rates(k_p_fe50__he4_mn47)*Y(jfe50)*state % rho + &
      screened_rates(k_p_fe51__he4_mn48)*Y(jfe51)*state % rho + &
      screened_rates(k_p_fe52__he4_mn49)*Y(jfe52)*state % rho + &
      screened_rates(k_p_fe53__he4_mn50)*Y(jfe53)*state % rho + &
      screened_rates(k_p_fe54__he4_mn51)*Y(jfe54)*state % rho + &
      screened_rates(k_p_fe55__he4_mn52)*Y(jfe55)*state % rho + &
      screened_rates(k_p_fe56__he4_mn53)*Y(jfe56)*state % rho + &
      screened_rates(k_p_he3__he4__weak__bet_pos_)*Y(jhe3)*state % rho - &
      screened_rates(k_p_he4__d_he3)*Y(jhe4)*state % rho - &
      screened_rates(k_p_he4_he4__d_be7)*Y(jhe4)**2*state % rho**2 + &
      screened_rates(k_p_k34__he4_ar31)*Y(jk34)*state % rho + &
      screened_rates(k_p_k35__he4_ar32)*Y(jk35)*state % rho + &
      screened_rates(k_p_k36__he4_ar33)*Y(jk36)*state % rho + &
      screened_rates(k_p_k37__he4_ar34)*Y(jk37)*state % rho + &
      screened_rates(k_p_k38__he4_ar35)*Y(jk38)*state % rho + &
      screened_rates(k_p_k39__he4_ar36)*Y(jk39)*state % rho + &
      screened_rates(k_p_k40__he4_ar37)*Y(jk40)*state % rho + &
      screened_rates(k_p_k41__he4_ar38)*Y(jk41)*state % rho + 2.0d0* &
      screened_rates(k_p_li7__he4_he4)*Y(jli7)*state % rho + &
      screened_rates(k_p_mg22__he4_na19)*Y(jmg22)*state % rho + &
      screened_rates(k_p_mg23__he4_na20)*Y(jmg23)*state % rho + &
      screened_rates(k_p_mg24__he4_na21)*Y(jmg24)*state % rho + &
      screened_rates(k_p_mg25__he4_na22)*Y(jmg25)*state % rho + &
      screened_rates(k_p_mg26__he4_na23)*Y(jmg26)*state % rho + &
      screened_rates(k_p_mn45__he4_cr42)*Y(jmn45)*state % rho + &
      screened_rates(k_p_mn46__he4_cr43)*Y(jmn46)*state % rho + &
      screened_rates(k_p_mn47__he4_cr44)*Y(jmn47)*state % rho + &
      screened_rates(k_p_mn48__he4_cr45)*Y(jmn48)*state % rho + &
      screened_rates(k_p_mn49__he4_cr46)*Y(jmn49)*state % rho + &
      screened_rates(k_p_mn50__he4_cr47)*Y(jmn50)*state % rho + &
      screened_rates(k_p_mn51__he4_cr48)*Y(jmn51)*state % rho + &
      screened_rates(k_p_mn52__he4_cr49)*Y(jmn52)*state % rho + &
      screened_rates(k_p_mn53__he4_cr50)*Y(jmn53)*state % rho + &
      screened_rates(k_p_mn54__he4_cr51)*Y(jmn54)*state % rho + &
      screened_rates(k_p_mn55__he4_cr52)*Y(jmn55)*state % rho + &
      screened_rates(k_p_n15__he4_c12)*Y(jn15)*state % rho + &
      screened_rates(k_p_na20__he4_ne17)*Y(jna20)*state % rho + &
      screened_rates(k_p_na21__he4_ne18)*Y(jna21)*state % rho + &
      screened_rates(k_p_na22__he4_ne19)*Y(jna22)*state % rho + &
      screened_rates(k_p_na23__he4_ne20)*Y(jna23)*state % rho + &
      screened_rates(k_p_na24__he4_ne21)*Y(jna24)*state % rho + &
      screened_rates(k_p_ne19__he4_f16)*Y(jne19)*state % rho + &
      screened_rates(k_p_ne20__he4_f17)*Y(jne20)*state % rho + &
      screened_rates(k_p_ne21__he4_f18)*Y(jne21)*state % rho + &
      screened_rates(k_p_ne22__he4_f19)*Y(jne22)*state % rho + &
      screened_rates(k_p_ni50__he4_co47)*Y(jni50)*state % rho + &
      screened_rates(k_p_ni51__he4_co48)*Y(jni51)*state % rho + &
      screened_rates(k_p_ni52__he4_co49)*Y(jni52)*state % rho + &
      screened_rates(k_p_ni53__he4_co50)*Y(jni53)*state % rho + &
      screened_rates(k_p_ni54__he4_co51)*Y(jni54)*state % rho + &
      screened_rates(k_p_ni55__he4_co52)*Y(jni55)*state % rho + &
      screened_rates(k_p_ni56__he4_co53)*Y(jni56)*state % rho + &
      screened_rates(k_p_o15__he4_n12)*Y(jo15)*state % rho + &
      screened_rates(k_p_o16__he4_n13)*Y(jo16)*state % rho + &
      screened_rates(k_p_o17__he4_n14)*Y(jo17)*state % rho + &
      screened_rates(k_p_o18__he4_n15)*Y(jo18)*state % rho + &
      screened_rates(k_p_p27__he4_si24)*Y(jp27)*state % rho + &
      screened_rates(k_p_p28__he4_si25)*Y(jp28)*state % rho + &
      screened_rates(k_p_p29__he4_si26)*Y(jp29)*state % rho + &
      screened_rates(k_p_p30__he4_si27)*Y(jp30)*state % rho + &
      screened_rates(k_p_p31__he4_si28)*Y(jp31)*state % rho + &
      screened_rates(k_p_p32__he4_si29)*Y(jp32)*state % rho + &
      screened_rates(k_p_p33__he4_si30)*Y(jp33)*state % rho + 1.0d0* &
      screened_rates(k_p_p_ca38__he4_ca36)*Y(jca38)*Y(jp)*state % rho**2 - &
      1.0d0*screened_rates(k_p_p_he4__he3_he3)*Y(jhe4)*Y(jp)*state % rho**2 - &
      1.0d0*screened_rates(k_p_p_he4_he4__he3_be7)*Y(jhe4)**2*Y(jp)*state % rho &
      **3 + screened_rates(k_p_s29__he4_p26)*Y(js29)*state % rho + &
      screened_rates(k_p_s30__he4_p27)*Y(js30)*state % rho + &
      screened_rates(k_p_s31__he4_p28)*Y(js31)*state % rho + &
      screened_rates(k_p_s32__he4_p29)*Y(js32)*state % rho + &
      screened_rates(k_p_s33__he4_p30)*Y(js33)*state % rho + &
      screened_rates(k_p_s34__he4_p31)*Y(js34)*state % rho + &
      screened_rates(k_p_s35__he4_p32)*Y(js35)*state % rho + &
      screened_rates(k_p_sc37__he4_ca34)*Y(jsc37)*state % rho + &
      screened_rates(k_p_sc38__he4_ca35)*Y(jsc38)*state % rho + &
      screened_rates(k_p_sc39__he4_ca36)*Y(jsc39)*state % rho + &
      screened_rates(k_p_sc40__he4_ca37)*Y(jsc40)*state % rho + &
      screened_rates(k_p_sc41__he4_ca38)*Y(jsc41)*state % rho + &
      screened_rates(k_p_sc42__he4_ca39)*Y(jsc42)*state % rho + &
      screened_rates(k_p_sc43__he4_ca40)*Y(jsc43)*state % rho + &
      screened_rates(k_p_sc44__he4_ca41)*Y(jsc44)*state % rho + &
      screened_rates(k_p_sc45__he4_ca42)*Y(jsc45)*state % rho + &
      screened_rates(k_p_sc46__he4_ca43)*Y(jsc46)*state % rho + &
      screened_rates(k_p_si25__he4_al22)*Y(jsi25)*state % rho + &
      screened_rates(k_p_si26__he4_al23)*Y(jsi26)*state % rho + &
      screened_rates(k_p_si27__he4_al24)*Y(jsi27)*state % rho + &
      screened_rates(k_p_si28__he4_al25)*Y(jsi28)*state % rho + &
      screened_rates(k_p_si29__he4_al26)*Y(jsi29)*state % rho + &
      screened_rates(k_p_si30__he4_al27)*Y(jsi30)*state % rho + &
      screened_rates(k_p_ti39__he4_sc36)*Y(jti39)*state % rho + &
      screened_rates(k_p_ti40__he4_sc37)*Y(jti40)*state % rho + &
      screened_rates(k_p_ti41__he4_sc38)*Y(jti41)*state % rho + &
      screened_rates(k_p_ti42__he4_sc39)*Y(jti42)*state % rho + &
      screened_rates(k_p_ti43__he4_sc40)*Y(jti43)*state % rho + &
      screened_rates(k_p_ti44__he4_sc41)*Y(jti44)*state % rho + &
      screened_rates(k_p_ti45__he4_sc42)*Y(jti45)*state % rho + &
      screened_rates(k_p_ti46__he4_sc43)*Y(jti46)*state % rho + &
      screened_rates(k_p_ti47__he4_sc44)*Y(jti47)*state % rho + &
      screened_rates(k_p_ti48__he4_sc45)*Y(jti48)*state % rho + &
      screened_rates(k_p_v41__he4_ti38)*Y(jv41)*state % rho + &
      screened_rates(k_p_v42__he4_ti39)*Y(jv42)*state % rho + &
      screened_rates(k_p_v43__he4_ti40)*Y(jv43)*state % rho + &
      screened_rates(k_p_v44__he4_ti41)*Y(jv44)*state % rho + &
      screened_rates(k_p_v45__he4_ti42)*Y(jv45)*state % rho + &
      screened_rates(k_p_v46__he4_ti43)*Y(jv46)*state % rho + &
      screened_rates(k_p_v47__he4_ti44)*Y(jv47)*state % rho + &
      screened_rates(k_p_v48__he4_ti45)*Y(jv48)*state % rho + &
      screened_rates(k_p_v49__he4_ti46)*Y(jv49)*state % rho + &
      screened_rates(k_p_v50__he4_ti47)*Y(jv50)*state % rho &
       )
    call set_jac_entry(state, jhe4, jp, scratch)

    scratch = (&
      2.0d0*screened_rates(k_d_be7__p_he4_he4)*Y(jbe7)*state % rho + 1.0d0* &
      screened_rates(k_d_d__he4)*Y(jd)*state % rho + screened_rates(k_d_he3__p_he4)* &
      Y(jhe3)*state % rho &
       )
    call set_jac_entry(state, jhe4, jd, scratch)

    scratch = (&
      screened_rates(k_d_he3__p_he4)*Y(jd)*state % rho + 2.0d0* &
      screened_rates(k_he3_be7__p_p_he4_he4)*Y(jbe7)*state % rho + 1.0d0* &
      screened_rates(k_he3_he3__p_p_he4)*Y(jhe3)*state % rho - &
      screened_rates(k_he4_he3__be7)*Y(jhe4)*state % rho + &
      screened_rates(k_p_he3__he4__weak__bet_pos_)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jhe3, scratch)

    scratch = (&
      -screened_rates(k_he4__d_d) - screened_rates(k_he4_al22__p26)*Y(jal22)*state % rho - &
      screened_rates(k_he4_al22__p_si25)*Y(jal22)*state % rho - &
      screened_rates(k_he4_al23__p27)*Y(jal23)*state % rho - &
      screened_rates(k_he4_al23__p_si26)*Y(jal23)*state % rho - &
      screened_rates(k_he4_al24__p28)*Y(jal24)*state % rho - &
      screened_rates(k_he4_al24__p_si27)*Y(jal24)*state % rho - &
      screened_rates(k_he4_al25__p29)*Y(jal25)*state % rho - &
      screened_rates(k_he4_al25__p_si28)*Y(jal25)*state % rho - &
      screened_rates(k_he4_al26__p30)*Y(jal26)*state % rho - &
      screened_rates(k_he4_al26__p_si29)*Y(jal26)*state % rho - &
      screened_rates(k_he4_al27__p31)*Y(jal27)*state % rho - &
      screened_rates(k_he4_al27__p_si30)*Y(jal27)*state % rho - &
      screened_rates(k_he4_al28__p32)*Y(jal28)*state % rho - &
      screened_rates(k_he4_ar31__ca35)*Y(jar31)*state % rho - &
      screened_rates(k_he4_ar31__p_k34)*Y(jar31)*state % rho - &
      screened_rates(k_he4_ar32__ca36)*Y(jar32)*state % rho - &
      screened_rates(k_he4_ar32__p_k35)*Y(jar32)*state % rho - &
      screened_rates(k_he4_ar33__ca37)*Y(jar33)*state % rho - &
      screened_rates(k_he4_ar33__p_k36)*Y(jar33)*state % rho - &
      screened_rates(k_he4_ar34__ca38)*Y(jar34)*state % rho - &
      screened_rates(k_he4_ar34__p_k37)*Y(jar34)*state % rho - &
      screened_rates(k_he4_ar35__ca39)*Y(jar35)*state % rho - &
      screened_rates(k_he4_ar35__p_k38)*Y(jar35)*state % rho - &
      screened_rates(k_he4_ar36__ca40)*Y(jar36)*state % rho - &
      screened_rates(k_he4_ar36__p_k39)*Y(jar36)*state % rho - &
      screened_rates(k_he4_ar37__ca41)*Y(jar37)*state % rho - &
      screened_rates(k_he4_ar37__p_k40)*Y(jar37)*state % rho - &
      screened_rates(k_he4_ar38__ca42)*Y(jar38)*state % rho - &
      screened_rates(k_he4_ar38__p_k41)*Y(jar38)*state % rho - &
      screened_rates(k_he4_ar39__ca43)*Y(jar39)*state % rho - &
      screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho - &
      screened_rates(k_he4_c12__p_n15)*Y(jc12)*state % rho - &
      screened_rates(k_he4_ca34__p_sc37)*Y(jca34)*state % rho - &
      screened_rates(k_he4_ca34__ti38)*Y(jca34)*state % rho - &
      screened_rates(k_he4_ca35__p_sc38)*Y(jca35)*state % rho - &
      screened_rates(k_he4_ca35__ti39)*Y(jca35)*state % rho - &
      screened_rates(k_he4_ca36__p_p_ca38)*Y(jca36)*state % rho - &
      screened_rates(k_he4_ca36__p_sc39)*Y(jca36)*state % rho - &
      screened_rates(k_he4_ca36__ti40)*Y(jca36)*state % rho - &
      screened_rates(k_he4_ca37__p_sc40)*Y(jca37)*state % rho - &
      screened_rates(k_he4_ca37__ti41)*Y(jca37)*state % rho - &
      screened_rates(k_he4_ca38__p_sc41)*Y(jca38)*state % rho - &
      screened_rates(k_he4_ca38__ti42)*Y(jca38)*state % rho - &
      screened_rates(k_he4_ca39__p_sc42)*Y(jca39)*state % rho - &
      screened_rates(k_he4_ca39__ti43)*Y(jca39)*state % rho - &
      screened_rates(k_he4_ca40__p_sc43)*Y(jca40)*state % rho - &
      screened_rates(k_he4_ca40__ti44)*Y(jca40)*state % rho - &
      screened_rates(k_he4_ca41__p_sc44)*Y(jca41)*state % rho - &
      screened_rates(k_he4_ca41__ti45)*Y(jca41)*state % rho - &
      screened_rates(k_he4_ca42__p_sc45)*Y(jca42)*state % rho - &
      screened_rates(k_he4_ca42__ti46)*Y(jca42)*state % rho - &
      screened_rates(k_he4_ca43__p_sc46)*Y(jca43)*state % rho - &
      screened_rates(k_he4_ca43__ti47)*Y(jca43)*state % rho - &
      screened_rates(k_he4_ca44__ti48)*Y(jca44)*state % rho - &
      screened_rates(k_he4_cl29__k33)*Y(jcl29)*state % rho - &
      screened_rates(k_he4_cl29__p_ar32)*Y(jcl29)*state % rho - &
      screened_rates(k_he4_cl30__k34)*Y(jcl30)*state % rho - &
      screened_rates(k_he4_cl30__p_ar33)*Y(jcl30)*state % rho - &
      screened_rates(k_he4_cl31__k35)*Y(jcl31)*state % rho - &
      screened_rates(k_he4_cl31__p_ar34)*Y(jcl31)*state % rho - &
      screened_rates(k_he4_cl32__k36)*Y(jcl32)*state % rho - &
      screened_rates(k_he4_cl32__p_ar35)*Y(jcl32)*state % rho - &
      screened_rates(k_he4_cl33__k37)*Y(jcl33)*state % rho - &
      screened_rates(k_he4_cl33__p_ar36)*Y(jcl33)*state % rho - &
      screened_rates(k_he4_cl34__k38)*Y(jcl34)*state % rho - &
      screened_rates(k_he4_cl34__p_ar37)*Y(jcl34)*state % rho - &
      screened_rates(k_he4_cl35__k39)*Y(jcl35)*state % rho - &
      screened_rates(k_he4_cl35__p_ar38)*Y(jcl35)*state % rho - &
      screened_rates(k_he4_cl36__k40)*Y(jcl36)*state % rho - &
      screened_rates(k_he4_cl36__p_ar39)*Y(jcl36)*state % rho - &
      screened_rates(k_he4_cl37__k41)*Y(jcl37)*state % rho - &
      screened_rates(k_he4_co47__p_ni50)*Y(jco47)*state % rho - &
      screened_rates(k_he4_co48__p_ni51)*Y(jco48)*state % rho - &
      screened_rates(k_he4_co49__p_ni52)*Y(jco49)*state % rho - &
      screened_rates(k_he4_co50__p_ni53)*Y(jco50)*state % rho - &
      screened_rates(k_he4_co51__p_ni54)*Y(jco51)*state % rho - &
      screened_rates(k_he4_co52__p_ni55)*Y(jco52)*state % rho - &
      screened_rates(k_he4_co53__p_ni56)*Y(jco53)*state % rho - &
      screened_rates(k_he4_cr42__fe46)*Y(jcr42)*state % rho - &
      screened_rates(k_he4_cr42__p_mn45)*Y(jcr42)*state % rho - &
      screened_rates(k_he4_cr43__fe47)*Y(jcr43)*state % rho - &
      screened_rates(k_he4_cr43__p_mn46)*Y(jcr43)*state % rho - &
      screened_rates(k_he4_cr44__fe48)*Y(jcr44)*state % rho - &
      screened_rates(k_he4_cr44__p_mn47)*Y(jcr44)*state % rho - &
      screened_rates(k_he4_cr45__fe49)*Y(jcr45)*state % rho - &
      screened_rates(k_he4_cr45__p_mn48)*Y(jcr45)*state % rho - &
      screened_rates(k_he4_cr46__fe50)*Y(jcr46)*state % rho - &
      screened_rates(k_he4_cr46__p_mn49)*Y(jcr46)*state % rho - &
      screened_rates(k_he4_cr47__fe51)*Y(jcr47)*state % rho - &
      screened_rates(k_he4_cr47__p_mn50)*Y(jcr47)*state % rho - &
      screened_rates(k_he4_cr48__fe52)*Y(jcr48)*state % rho - &
      screened_rates(k_he4_cr48__p_mn51)*Y(jcr48)*state % rho - &
      screened_rates(k_he4_cr49__fe53)*Y(jcr49)*state % rho - &
      screened_rates(k_he4_cr49__p_mn52)*Y(jcr49)*state % rho - &
      screened_rates(k_he4_cr50__fe54)*Y(jcr50)*state % rho - &
      screened_rates(k_he4_cr50__p_mn53)*Y(jcr50)*state % rho - &
      screened_rates(k_he4_cr51__fe55)*Y(jcr51)*state % rho - &
      screened_rates(k_he4_cr51__p_mn54)*Y(jcr51)*state % rho - &
      screened_rates(k_he4_cr52__fe56)*Y(jcr52)*state % rho - &
      screened_rates(k_he4_cr52__p_mn55)*Y(jcr52)*state % rho - &
      screened_rates(k_he4_f16__na20)*Y(jf16)*state % rho - &
      screened_rates(k_he4_f16__p_ne19)*Y(jf16)*state % rho - &
      screened_rates(k_he4_f17__na21)*Y(jf17)*state % rho - &
      screened_rates(k_he4_f17__p_ne20)*Y(jf17)*state % rho - &
      screened_rates(k_he4_f18__na22)*Y(jf18)*state % rho - &
      screened_rates(k_he4_f18__p_ne21)*Y(jf18)*state % rho - &
      screened_rates(k_he4_f19__na23)*Y(jf19)*state % rho - &
      screened_rates(k_he4_f19__p_ne22)*Y(jf19)*state % rho - &
      screened_rates(k_he4_f20__na24)*Y(jf20)*state % rho - &
      screened_rates(k_he4_fe45__ni49)*Y(jfe45)*state % rho - &
      screened_rates(k_he4_fe45__p_co48)*Y(jfe45)*state % rho - &
      screened_rates(k_he4_fe46__ni50)*Y(jfe46)*state % rho - &
      screened_rates(k_he4_fe46__p_co49)*Y(jfe46)*state % rho - &
      screened_rates(k_he4_fe47__ni51)*Y(jfe47)*state % rho - &
      screened_rates(k_he4_fe47__p_co50)*Y(jfe47)*state % rho - &
      screened_rates(k_he4_fe48__ni52)*Y(jfe48)*state % rho - &
      screened_rates(k_he4_fe48__p_co51)*Y(jfe48)*state % rho - &
      screened_rates(k_he4_fe49__ni53)*Y(jfe49)*state % rho - &
      screened_rates(k_he4_fe49__p_co52)*Y(jfe49)*state % rho - &
      screened_rates(k_he4_fe50__ni54)*Y(jfe50)*state % rho - &
      screened_rates(k_he4_fe50__p_co53)*Y(jfe50)*state % rho - &
      screened_rates(k_he4_fe51__ni55)*Y(jfe51)*state % rho - &
      screened_rates(k_he4_fe51__p_co54)*Y(jfe51)*state % rho - &
      screened_rates(k_he4_fe52__ni56)*Y(jfe52)*state % rho - &
      screened_rates(k_he4_fe52__p_co55)*Y(jfe52)*state % rho - &
      screened_rates(k_he4_fe53__p_co56)*Y(jfe53)*state % rho - &
      screened_rates(k_he4_he3__be7)*Y(jhe3)*state % rho - 2.0d0* &
      screened_rates(k_he4_he4__p_li7)*Y(jhe4)*state % rho - 1.5d0* &
      screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**2*state % rho**2 - &
      screened_rates(k_he4_k33__p_ca36)*Y(jk33)*state % rho - &
      screened_rates(k_he4_k33__sc37)*Y(jk33)*state % rho - &
      screened_rates(k_he4_k34__p_ca37)*Y(jk34)*state % rho - &
      screened_rates(k_he4_k34__sc38)*Y(jk34)*state % rho - &
      screened_rates(k_he4_k35__p_ca38)*Y(jk35)*state % rho - &
      screened_rates(k_he4_k35__sc39)*Y(jk35)*state % rho - &
      screened_rates(k_he4_k36__p_ca39)*Y(jk36)*state % rho - &
      screened_rates(k_he4_k36__sc40)*Y(jk36)*state % rho - &
      screened_rates(k_he4_k37__p_ca40)*Y(jk37)*state % rho - &
      screened_rates(k_he4_k37__sc41)*Y(jk37)*state % rho - &
      screened_rates(k_he4_k38__p_ca41)*Y(jk38)*state % rho - &
      screened_rates(k_he4_k38__sc42)*Y(jk38)*state % rho - &
      screened_rates(k_he4_k39__p_ca42)*Y(jk39)*state % rho - &
      screened_rates(k_he4_k39__sc43)*Y(jk39)*state % rho - &
      screened_rates(k_he4_k40__p_ca43)*Y(jk40)*state % rho - &
      screened_rates(k_he4_k40__sc44)*Y(jk40)*state % rho - &
      screened_rates(k_he4_k41__p_ca44)*Y(jk41)*state % rho - &
      screened_rates(k_he4_k41__sc45)*Y(jk41)*state % rho - &
      screened_rates(k_he4_mg21__p_al24)*Y(jmg21)*state % rho - &
      screened_rates(k_he4_mg21__si25)*Y(jmg21)*state % rho - &
      screened_rates(k_he4_mg22__p_al25)*Y(jmg22)*state % rho - &
      screened_rates(k_he4_mg22__si26)*Y(jmg22)*state % rho - &
      screened_rates(k_he4_mg23__p_al26)*Y(jmg23)*state % rho - &
      screened_rates(k_he4_mg23__si27)*Y(jmg23)*state % rho - &
      screened_rates(k_he4_mg24__p_al27)*Y(jmg24)*state % rho - &
      screened_rates(k_he4_mg24__si28)*Y(jmg24)*state % rho - &
      screened_rates(k_he4_mg25__p_al28)*Y(jmg25)*state % rho - &
      screened_rates(k_he4_mg25__si29)*Y(jmg25)*state % rho - &
      screened_rates(k_he4_mg26__si30)*Y(jmg26)*state % rho - &
      screened_rates(k_he4_mn44__co48)*Y(jmn44)*state % rho - &
      screened_rates(k_he4_mn44__p_fe47)*Y(jmn44)*state % rho - &
      screened_rates(k_he4_mn45__co49)*Y(jmn45)*state % rho - &
      screened_rates(k_he4_mn45__p_fe48)*Y(jmn45)*state % rho - &
      screened_rates(k_he4_mn46__co50)*Y(jmn46)*state % rho - &
      screened_rates(k_he4_mn46__p_fe49)*Y(jmn46)*state % rho - &
      screened_rates(k_he4_mn47__co51)*Y(jmn47)*state % rho - &
      screened_rates(k_he4_mn47__p_fe50)*Y(jmn47)*state % rho - &
      screened_rates(k_he4_mn48__co52)*Y(jmn48)*state % rho - &
      screened_rates(k_he4_mn48__p_fe51)*Y(jmn48)*state % rho - &
      screened_rates(k_he4_mn49__co53)*Y(jmn49)*state % rho - &
      screened_rates(k_he4_mn49__p_fe52)*Y(jmn49)*state % rho - &
      screened_rates(k_he4_mn50__co54)*Y(jmn50)*state % rho - &
      screened_rates(k_he4_mn50__p_fe53)*Y(jmn50)*state % rho - &
      screened_rates(k_he4_mn51__co55)*Y(jmn51)*state % rho - &
      screened_rates(k_he4_mn51__p_fe54)*Y(jmn51)*state % rho - &
      screened_rates(k_he4_mn52__co56)*Y(jmn52)*state % rho - &
      screened_rates(k_he4_mn52__p_fe55)*Y(jmn52)*state % rho - &
      screened_rates(k_he4_mn53__p_fe56)*Y(jmn53)*state % rho - &
      screened_rates(k_he4_n12__p_o15)*Y(jn12)*state % rho - &
      screened_rates(k_he4_n13__p_o16)*Y(jn13)*state % rho - &
      screened_rates(k_he4_n14__f18)*Y(jn14)*state % rho - &
      screened_rates(k_he4_n14__p_o17)*Y(jn14)*state % rho - &
      screened_rates(k_he4_n15__f19)*Y(jn15)*state % rho - &
      screened_rates(k_he4_n15__p_o18)*Y(jn15)*state % rho - &
      screened_rates(k_he4_na19__al23)*Y(jna19)*state % rho - &
      screened_rates(k_he4_na19__p_mg22)*Y(jna19)*state % rho - &
      screened_rates(k_he4_na20__al24)*Y(jna20)*state % rho - &
      screened_rates(k_he4_na20__p_mg23)*Y(jna20)*state % rho - &
      screened_rates(k_he4_na21__al25)*Y(jna21)*state % rho - &
      screened_rates(k_he4_na21__p_mg24)*Y(jna21)*state % rho - &
      screened_rates(k_he4_na22__al26)*Y(jna22)*state % rho - &
      screened_rates(k_he4_na22__p_mg25)*Y(jna22)*state % rho - &
      screened_rates(k_he4_na23__al27)*Y(jna23)*state % rho - &
      screened_rates(k_he4_na23__p_mg26)*Y(jna23)*state % rho - &
      screened_rates(k_he4_na24__al28)*Y(jna24)*state % rho - &
      screened_rates(k_he4_ne17__mg21)*Y(jne17)*state % rho - &
      screened_rates(k_he4_ne17__p_na20)*Y(jne17)*state % rho - &
      screened_rates(k_he4_ne18__mg22)*Y(jne18)*state % rho - &
      screened_rates(k_he4_ne18__p_na21)*Y(jne18)*state % rho - &
      screened_rates(k_he4_ne19__mg23)*Y(jne19)*state % rho - &
      screened_rates(k_he4_ne19__p_na22)*Y(jne19)*state % rho - &
      screened_rates(k_he4_ne20__c12_c12)*Y(jne20)*state % rho - &
      screened_rates(k_he4_ne20__mg24)*Y(jne20)*state % rho - &
      screened_rates(k_he4_ne20__p_na23)*Y(jne20)*state % rho - &
      screened_rates(k_he4_ne21__mg25)*Y(jne21)*state % rho - &
      screened_rates(k_he4_ne21__p_na24)*Y(jne21)*state % rho - &
      screened_rates(k_he4_ne22__mg26)*Y(jne22)*state % rho - &
      screened_rates(k_he4_o14__ne18)*Y(jo14)*state % rho - &
      screened_rates(k_he4_o14__p_f17)*Y(jo14)*state % rho - &
      screened_rates(k_he4_o15__ne19)*Y(jo15)*state % rho - &
      screened_rates(k_he4_o15__p_f18)*Y(jo15)*state % rho - &
      screened_rates(k_he4_o16__ne20)*Y(jo16)*state % rho - &
      screened_rates(k_he4_o16__p_f19)*Y(jo16)*state % rho - &
      screened_rates(k_he4_o17__ne21)*Y(jo17)*state % rho - &
      screened_rates(k_he4_o17__p_f20)*Y(jo17)*state % rho - &
      screened_rates(k_he4_o18__ne22)*Y(jo18)*state % rho - &
      screened_rates(k_he4_p26__cl30)*Y(jp26)*state % rho - &
      screened_rates(k_he4_p26__p_s29)*Y(jp26)*state % rho - &
      screened_rates(k_he4_p27__cl31)*Y(jp27)*state % rho - &
      screened_rates(k_he4_p27__p_s30)*Y(jp27)*state % rho - &
      screened_rates(k_he4_p28__cl32)*Y(jp28)*state % rho - &
      screened_rates(k_he4_p28__p_s31)*Y(jp28)*state % rho - &
      screened_rates(k_he4_p29__cl33)*Y(jp29)*state % rho - &
      screened_rates(k_he4_p29__p_s32)*Y(jp29)*state % rho - &
      screened_rates(k_he4_p30__cl34)*Y(jp30)*state % rho - &
      screened_rates(k_he4_p30__p_s33)*Y(jp30)*state % rho - &
      screened_rates(k_he4_p31__cl35)*Y(jp31)*state % rho - &
      screened_rates(k_he4_p31__p_s34)*Y(jp31)*state % rho - &
      screened_rates(k_he4_p32__cl36)*Y(jp32)*state % rho - &
      screened_rates(k_he4_p32__p_s35)*Y(jp32)*state % rho - &
      screened_rates(k_he4_p33__cl37)*Y(jp33)*state % rho - &
      screened_rates(k_he4_s28__ar32)*Y(js28)*state % rho - &
      screened_rates(k_he4_s28__p_cl31)*Y(js28)*state % rho - &
      screened_rates(k_he4_s29__ar33)*Y(js29)*state % rho - &
      screened_rates(k_he4_s29__p_cl32)*Y(js29)*state % rho - &
      screened_rates(k_he4_s30__ar34)*Y(js30)*state % rho - &
      screened_rates(k_he4_s30__p_cl33)*Y(js30)*state % rho - &
      screened_rates(k_he4_s31__ar35)*Y(js31)*state % rho - &
      screened_rates(k_he4_s31__p_cl34)*Y(js31)*state % rho - &
      screened_rates(k_he4_s32__ar36)*Y(js32)*state % rho - &
      screened_rates(k_he4_s32__p_cl35)*Y(js32)*state % rho - &
      screened_rates(k_he4_s33__ar37)*Y(js33)*state % rho - &
      screened_rates(k_he4_s33__p_cl36)*Y(js33)*state % rho - &
      screened_rates(k_he4_s34__ar38)*Y(js34)*state % rho - &
      screened_rates(k_he4_s34__p_cl37)*Y(js34)*state % rho - &
      screened_rates(k_he4_s35__ar39)*Y(js35)*state % rho - &
      screened_rates(k_he4_sc36__p_ti39)*Y(jsc36)*state % rho - &
      screened_rates(k_he4_sc36__v40)*Y(jsc36)*state % rho - &
      screened_rates(k_he4_sc37__p_ti40)*Y(jsc37)*state % rho - &
      screened_rates(k_he4_sc37__v41)*Y(jsc37)*state % rho - &
      screened_rates(k_he4_sc38__p_ti41)*Y(jsc38)*state % rho - &
      screened_rates(k_he4_sc38__v42)*Y(jsc38)*state % rho - &
      screened_rates(k_he4_sc39__p_ti42)*Y(jsc39)*state % rho - &
      screened_rates(k_he4_sc39__v43)*Y(jsc39)*state % rho - &
      screened_rates(k_he4_sc40__p_ti43)*Y(jsc40)*state % rho - &
      screened_rates(k_he4_sc40__v44)*Y(jsc40)*state % rho - &
      screened_rates(k_he4_sc41__p_ti44)*Y(jsc41)*state % rho - &
      screened_rates(k_he4_sc41__v45)*Y(jsc41)*state % rho - &
      screened_rates(k_he4_sc42__p_ti45)*Y(jsc42)*state % rho - &
      screened_rates(k_he4_sc42__v46)*Y(jsc42)*state % rho - &
      screened_rates(k_he4_sc43__p_ti46)*Y(jsc43)*state % rho - &
      screened_rates(k_he4_sc43__v47)*Y(jsc43)*state % rho - &
      screened_rates(k_he4_sc44__p_ti47)*Y(jsc44)*state % rho - &
      screened_rates(k_he4_sc44__v48)*Y(jsc44)*state % rho - &
      screened_rates(k_he4_sc45__p_ti48)*Y(jsc45)*state % rho - &
      screened_rates(k_he4_sc45__v49)*Y(jsc45)*state % rho - &
      screened_rates(k_he4_sc46__v50)*Y(jsc46)*state % rho - &
      screened_rates(k_he4_si24__p_p27)*Y(jsi24)*state % rho - &
      screened_rates(k_he4_si24__s28)*Y(jsi24)*state % rho - &
      screened_rates(k_he4_si25__p_p28)*Y(jsi25)*state % rho - &
      screened_rates(k_he4_si25__s29)*Y(jsi25)*state % rho - &
      screened_rates(k_he4_si26__p_p29)*Y(jsi26)*state % rho - &
      screened_rates(k_he4_si26__s30)*Y(jsi26)*state % rho - &
      screened_rates(k_he4_si27__p_p30)*Y(jsi27)*state % rho - &
      screened_rates(k_he4_si27__s31)*Y(jsi27)*state % rho - &
      screened_rates(k_he4_si28__p_p31)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_si28__s32)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_si29__p_p32)*Y(jsi29)*state % rho - &
      screened_rates(k_he4_si29__s33)*Y(jsi29)*state % rho - &
      screened_rates(k_he4_si30__p_p33)*Y(jsi30)*state % rho - &
      screened_rates(k_he4_si30__s34)*Y(jsi30)*state % rho - &
      screened_rates(k_he4_ti38__cr42)*Y(jti38)*state % rho - &
      screened_rates(k_he4_ti38__p_v41)*Y(jti38)*state % rho - &
      screened_rates(k_he4_ti39__cr43)*Y(jti39)*state % rho - &
      screened_rates(k_he4_ti39__p_v42)*Y(jti39)*state % rho - &
      screened_rates(k_he4_ti40__cr44)*Y(jti40)*state % rho - &
      screened_rates(k_he4_ti40__p_v43)*Y(jti40)*state % rho - &
      screened_rates(k_he4_ti41__cr45)*Y(jti41)*state % rho - &
      screened_rates(k_he4_ti41__p_v44)*Y(jti41)*state % rho - &
      screened_rates(k_he4_ti42__cr46)*Y(jti42)*state % rho - &
      screened_rates(k_he4_ti42__p_v45)*Y(jti42)*state % rho - &
      screened_rates(k_he4_ti43__cr47)*Y(jti43)*state % rho - &
      screened_rates(k_he4_ti43__p_v46)*Y(jti43)*state % rho - &
      screened_rates(k_he4_ti44__cr48)*Y(jti44)*state % rho - &
      screened_rates(k_he4_ti44__p_v47)*Y(jti44)*state % rho - &
      screened_rates(k_he4_ti45__cr49)*Y(jti45)*state % rho - &
      screened_rates(k_he4_ti45__p_v48)*Y(jti45)*state % rho - &
      screened_rates(k_he4_ti46__cr50)*Y(jti46)*state % rho - &
      screened_rates(k_he4_ti46__p_v49)*Y(jti46)*state % rho - &
      screened_rates(k_he4_ti47__cr51)*Y(jti47)*state % rho - &
      screened_rates(k_he4_ti47__p_v50)*Y(jti47)*state % rho - &
      screened_rates(k_he4_ti48__cr52)*Y(jti48)*state % rho - &
      screened_rates(k_he4_v40__mn44)*Y(jv40)*state % rho - &
      screened_rates(k_he4_v40__p_cr43)*Y(jv40)*state % rho - &
      screened_rates(k_he4_v41__mn45)*Y(jv41)*state % rho - &
      screened_rates(k_he4_v41__p_cr44)*Y(jv41)*state % rho - &
      screened_rates(k_he4_v42__mn46)*Y(jv42)*state % rho - &
      screened_rates(k_he4_v42__p_cr45)*Y(jv42)*state % rho - &
      screened_rates(k_he4_v43__mn47)*Y(jv43)*state % rho - &
      screened_rates(k_he4_v43__p_cr46)*Y(jv43)*state % rho - &
      screened_rates(k_he4_v44__mn48)*Y(jv44)*state % rho - &
      screened_rates(k_he4_v44__p_cr47)*Y(jv44)*state % rho - &
      screened_rates(k_he4_v45__mn49)*Y(jv45)*state % rho - &
      screened_rates(k_he4_v45__p_cr48)*Y(jv45)*state % rho - &
      screened_rates(k_he4_v46__mn50)*Y(jv46)*state % rho - &
      screened_rates(k_he4_v46__p_cr49)*Y(jv46)*state % rho - &
      screened_rates(k_he4_v47__mn51)*Y(jv47)*state % rho - &
      screened_rates(k_he4_v47__p_cr50)*Y(jv47)*state % rho - &
      screened_rates(k_he4_v48__mn52)*Y(jv48)*state % rho - &
      screened_rates(k_he4_v48__p_cr51)*Y(jv48)*state % rho - &
      screened_rates(k_he4_v49__mn53)*Y(jv49)*state % rho - &
      screened_rates(k_he4_v49__p_cr52)*Y(jv49)*state % rho - &
      screened_rates(k_he4_v50__mn54)*Y(jv50)*state % rho - screened_rates(k_p_he4__d_he3) &
      *Y(jp)*state % rho - 2.0d0*screened_rates(k_p_he4_he4__d_be7)*Y(jhe4)* &
      Y(jp)*state % rho**2 - 0.5d0*screened_rates(k_p_p_he4__he3_he3)*Y(jp)**2* &
      state % rho**2 - 1.0d0*screened_rates(k_p_p_he4_he4__he3_be7)*Y(jhe4)* &
      Y(jp)**2*state % rho**3 &
       )
    call set_jac_entry(state, jhe4, jhe4, scratch)

    scratch = (&
      2.0d0*screened_rates(k_p_li7__he4_he4)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jli7, scratch)

    scratch = (&
      screened_rates(k_be7__he4_he3) + 2.0d0*screened_rates(k_d_be7__p_he4_he4)*Y(jd)*state % rho &
      + 2.0d0*screened_rates(k_he3_be7__p_p_he4_he4)*Y(jhe3)*state % rho &
       )
    call set_jac_entry(state, jhe4, jbe7, scratch)

    scratch = (&
      2.0d0*screened_rates(k_b8__he4_he4__weak__wc12) &
       )
    call set_jac_entry(state, jhe4, jb8, scratch)

    scratch = (&
      3.0d0*screened_rates(k_c12__he4_he4_he4) + 1.0d0*screened_rates(k_c12_c12__he4_ne20)* &
      Y(jc12)*state % rho - screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_c12__p_n15)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jc12, scratch)

    scratch = (&
      -screened_rates(k_he4_n12__p_o15)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jn12, scratch)

    scratch = (&
      -screened_rates(k_he4_n13__p_o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jn13, scratch)

    scratch = (&
      -screened_rates(k_he4_n14__f18)*Y(jhe4)*state % rho - screened_rates(k_he4_n14__p_o17)* &
      Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jn14, scratch)

    scratch = (&
      -screened_rates(k_he4_n15__f19)*Y(jhe4)*state % rho - screened_rates(k_he4_n15__p_o18)* &
      Y(jhe4)*state % rho + screened_rates(k_p_n15__he4_c12)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jn15, scratch)

    scratch = (&
      -screened_rates(k_he4_o14__ne18)*Y(jhe4)*state % rho - screened_rates(k_he4_o14__p_f17)* &
      Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jo14, scratch)

    scratch = (&
      -screened_rates(k_he4_o15__ne19)*Y(jhe4)*state % rho - screened_rates(k_he4_o15__p_f18)* &
      Y(jhe4)*state % rho + screened_rates(k_p_o15__he4_n12)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jo15, scratch)

    scratch = (&
      -screened_rates(k_he4_o16__ne20)*Y(jhe4)*state % rho - screened_rates(k_he4_o16__p_f19)* &
      Y(jhe4)*state % rho + screened_rates(k_o16__he4_c12) + &
      screened_rates(k_p_o16__he4_n13)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jo16, scratch)

    scratch = (&
      -screened_rates(k_he4_o17__ne21)*Y(jhe4)*state % rho - screened_rates(k_he4_o17__p_f20)* &
      Y(jhe4)*state % rho + screened_rates(k_p_o17__he4_n14)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jo17, scratch)

    scratch = (&
      -screened_rates(k_he4_o18__ne22)*Y(jhe4)*state % rho + screened_rates(k_p_o18__he4_n15)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jo18, scratch)

    scratch = (&
      -screened_rates(k_he4_f16__na20)*Y(jhe4)*state % rho - screened_rates(k_he4_f16__p_ne19)* &
      Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jf16, scratch)

    scratch = (&
      -screened_rates(k_he4_f17__na21)*Y(jhe4)*state % rho - screened_rates(k_he4_f17__p_ne20)* &
      Y(jhe4)*state % rho + screened_rates(k_p_f17__he4_o14)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jf17, scratch)

    scratch = (&
      screened_rates(k_f18__he4_n14) - screened_rates(k_he4_f18__na22)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_f18__p_ne21)*Y(jhe4)*state % rho + &
      screened_rates(k_p_f18__he4_o15)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jf18, scratch)

    scratch = (&
      screened_rates(k_f19__he4_n15) - screened_rates(k_he4_f19__na23)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_f19__p_ne22)*Y(jhe4)*state % rho + &
      screened_rates(k_p_f19__he4_o16)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jf19, scratch)

    scratch = (&
      -screened_rates(k_he4_f20__na24)*Y(jhe4)*state % rho + screened_rates(k_p_f20__he4_o17)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jf20, scratch)

    scratch = (&
      -screened_rates(k_he4_ne17__mg21)*Y(jhe4)*state % rho - screened_rates(k_he4_ne17__p_na20)* &
      Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jne17, scratch)

    scratch = (&
      -screened_rates(k_he4_ne18__mg22)*Y(jhe4)*state % rho - screened_rates(k_he4_ne18__p_na21)* &
      Y(jhe4)*state % rho + screened_rates(k_ne18__he4_o14) &
       )
    call set_jac_entry(state, jhe4, jne18, scratch)

    scratch = (&
      -screened_rates(k_he4_ne19__mg23)*Y(jhe4)*state % rho - screened_rates(k_he4_ne19__p_na22)* &
      Y(jhe4)*state % rho + screened_rates(k_ne19__he4_o15) + &
      screened_rates(k_p_ne19__he4_f16)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jne19, scratch)

    scratch = (&
      -screened_rates(k_he4_ne20__c12_c12)*Y(jhe4)*state % rho - screened_rates(k_he4_ne20__mg24) &
      *Y(jhe4)*state % rho - screened_rates(k_he4_ne20__p_na23)*Y(jhe4)* &
      state % rho + screened_rates(k_ne20__he4_o16) + screened_rates(k_p_ne20__he4_f17)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jne20, scratch)

    scratch = (&
      -screened_rates(k_he4_ne21__mg25)*Y(jhe4)*state % rho - screened_rates(k_he4_ne21__p_na24)* &
      Y(jhe4)*state % rho + screened_rates(k_ne21__he4_o17) + &
      screened_rates(k_p_ne21__he4_f18)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jne21, scratch)

    scratch = (&
      -screened_rates(k_he4_ne22__mg26)*Y(jhe4)*state % rho + screened_rates(k_p_ne22__he4_f19)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jne22, scratch)

    scratch = (&
      -screened_rates(k_he4_na19__al23)*Y(jhe4)*state % rho - screened_rates(k_he4_na19__p_mg22)* &
      Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jna19, scratch)

    scratch = (&
      -screened_rates(k_he4_na20__al24)*Y(jhe4)*state % rho - screened_rates(k_he4_na20__p_mg23)* &
      Y(jhe4)*state % rho + screened_rates(k_na20__he4_f16) + &
      screened_rates(k_na20__he4_o16__weak__wc12) + screened_rates(k_p_na20__he4_ne17)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jna20, scratch)

    scratch = (&
      -screened_rates(k_he4_na21__al25)*Y(jhe4)*state % rho - screened_rates(k_he4_na21__p_mg24)* &
      Y(jhe4)*state % rho + screened_rates(k_na21__he4_f17) + &
      screened_rates(k_p_na21__he4_ne18)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jna21, scratch)

    scratch = (&
      -screened_rates(k_he4_na22__al26)*Y(jhe4)*state % rho - screened_rates(k_he4_na22__p_mg25)* &
      Y(jhe4)*state % rho + screened_rates(k_na22__he4_f18) + &
      screened_rates(k_p_na22__he4_ne19)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jna22, scratch)

    scratch = (&
      -screened_rates(k_he4_na23__al27)*Y(jhe4)*state % rho - screened_rates(k_he4_na23__p_mg26)* &
      Y(jhe4)*state % rho + screened_rates(k_na23__he4_f19) + &
      screened_rates(k_p_na23__he4_ne20)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jna23, scratch)

    scratch = (&
      -screened_rates(k_he4_na24__al28)*Y(jhe4)*state % rho + screened_rates(k_p_na24__he4_ne21)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jna24, scratch)

    scratch = (&
      -screened_rates(k_he4_mg21__p_al24)*Y(jhe4)*state % rho - screened_rates(k_he4_mg21__si25)* &
      Y(jhe4)*state % rho + screened_rates(k_mg21__he4_ne17) &
       )
    call set_jac_entry(state, jhe4, jmg21, scratch)

    scratch = (&
      -screened_rates(k_he4_mg22__p_al25)*Y(jhe4)*state % rho - screened_rates(k_he4_mg22__si26)* &
      Y(jhe4)*state % rho + screened_rates(k_mg22__he4_ne18) + &
      screened_rates(k_p_mg22__he4_na19)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jmg22, scratch)

    scratch = (&
      -screened_rates(k_he4_mg23__p_al26)*Y(jhe4)*state % rho - screened_rates(k_he4_mg23__si27)* &
      Y(jhe4)*state % rho + screened_rates(k_mg23__he4_ne19) + &
      screened_rates(k_p_mg23__he4_na20)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jmg23, scratch)

    scratch = (&
      -screened_rates(k_he4_mg24__p_al27)*Y(jhe4)*state % rho - screened_rates(k_he4_mg24__si28)* &
      Y(jhe4)*state % rho + screened_rates(k_mg24__he4_ne20) + &
      screened_rates(k_p_mg24__he4_na21)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jmg24, scratch)

    scratch = (&
      -screened_rates(k_he4_mg25__p_al28)*Y(jhe4)*state % rho - screened_rates(k_he4_mg25__si29)* &
      Y(jhe4)*state % rho + screened_rates(k_mg25__he4_ne21) + &
      screened_rates(k_p_mg25__he4_na22)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jmg25, scratch)

    scratch = (&
      -screened_rates(k_he4_mg26__si30)*Y(jhe4)*state % rho + screened_rates(k_mg26__he4_ne22) + &
      screened_rates(k_p_mg26__he4_na23)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jmg26, scratch)

    scratch = (&
      screened_rates(k_al22__he4_ne18__weak__wc12) - screened_rates(k_he4_al22__p26)*Y(jhe4)* &
      state % rho - screened_rates(k_he4_al22__p_si25)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jal22, scratch)

    scratch = (&
      screened_rates(k_al23__he4_na19) - screened_rates(k_he4_al23__p27)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al23__p_si26)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jal23, scratch)

    scratch = (&
      screened_rates(k_al24__he4_na20) + screened_rates(k_al24__he4_ne20__weak__wc12) - &
      screened_rates(k_he4_al24__p28)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al24__p_si27)*Y(jhe4)*state % rho + &
      screened_rates(k_p_al24__he4_mg21)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jal24, scratch)

    scratch = (&
      screened_rates(k_al25__he4_na21) - screened_rates(k_he4_al25__p29)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al25__p_si28)*Y(jhe4)*state % rho + &
      screened_rates(k_p_al25__he4_mg22)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jal25, scratch)

    scratch = (&
      screened_rates(k_al26__he4_na22) - screened_rates(k_he4_al26__p30)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al26__p_si29)*Y(jhe4)*state % rho + &
      screened_rates(k_p_al26__he4_mg23)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jal26, scratch)

    scratch = (&
      screened_rates(k_al27__he4_na23) - screened_rates(k_he4_al27__p31)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al27__p_si30)*Y(jhe4)*state % rho + &
      screened_rates(k_p_al27__he4_mg24)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jal27, scratch)

    scratch = (&
      screened_rates(k_al28__he4_na24) - screened_rates(k_he4_al28__p32)*Y(jhe4)*state % rho + &
      screened_rates(k_p_al28__he4_mg25)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jal28, scratch)

    scratch = (&
      -screened_rates(k_he4_si24__p_p27)*Y(jhe4)*state % rho - screened_rates(k_he4_si24__s28)* &
      Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jsi24, scratch)

    scratch = (&
      -screened_rates(k_he4_si25__p_p28)*Y(jhe4)*state % rho - screened_rates(k_he4_si25__s29)* &
      Y(jhe4)*state % rho + screened_rates(k_p_si25__he4_al22)*Y(jp)*state % rho + &
      screened_rates(k_si25__he4_mg21) &
       )
    call set_jac_entry(state, jhe4, jsi25, scratch)

    scratch = (&
      -screened_rates(k_he4_si26__p_p29)*Y(jhe4)*state % rho - screened_rates(k_he4_si26__s30)* &
      Y(jhe4)*state % rho + screened_rates(k_p_si26__he4_al23)*Y(jp)*state % rho + &
      screened_rates(k_si26__he4_mg22) &
       )
    call set_jac_entry(state, jhe4, jsi26, scratch)

    scratch = (&
      -screened_rates(k_he4_si27__p_p30)*Y(jhe4)*state % rho - screened_rates(k_he4_si27__s31)* &
      Y(jhe4)*state % rho + screened_rates(k_p_si27__he4_al24)*Y(jp)*state % rho + &
      screened_rates(k_si27__he4_mg23) &
       )
    call set_jac_entry(state, jhe4, jsi27, scratch)

    scratch = (&
      -screened_rates(k_he4_si28__p_p31)*Y(jhe4)*state % rho - screened_rates(k_he4_si28__s32)* &
      Y(jhe4)*state % rho + screened_rates(k_p_si28__he4_al25)*Y(jp)*state % rho + &
      screened_rates(k_si28__he4_mg24) &
       )
    call set_jac_entry(state, jhe4, jsi28, scratch)

    scratch = (&
      -screened_rates(k_he4_si29__p_p32)*Y(jhe4)*state % rho - screened_rates(k_he4_si29__s33)* &
      Y(jhe4)*state % rho + screened_rates(k_p_si29__he4_al26)*Y(jp)*state % rho + &
      screened_rates(k_si29__he4_mg25) &
       )
    call set_jac_entry(state, jhe4, jsi29, scratch)

    scratch = (&
      -screened_rates(k_he4_si30__p_p33)*Y(jhe4)*state % rho - screened_rates(k_he4_si30__s34)* &
      Y(jhe4)*state % rho + screened_rates(k_p_si30__he4_al27)*Y(jp)*state % rho + &
      screened_rates(k_si30__he4_mg26) &
       )
    call set_jac_entry(state, jhe4, jsi30, scratch)

    scratch = (&
      -screened_rates(k_he4_p26__cl30)*Y(jhe4)*state % rho - screened_rates(k_he4_p26__p_s29)* &
      Y(jhe4)*state % rho + screened_rates(k_p26__he4_al22) &
       )
    call set_jac_entry(state, jhe4, jp26, scratch)

    scratch = (&
      -screened_rates(k_he4_p27__cl31)*Y(jhe4)*state % rho - screened_rates(k_he4_p27__p_s30)* &
      Y(jhe4)*state % rho + screened_rates(k_p27__he4_al23) + &
      screened_rates(k_p_p27__he4_si24)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jp27, scratch)

    scratch = (&
      -screened_rates(k_he4_p28__cl32)*Y(jhe4)*state % rho - screened_rates(k_he4_p28__p_s31)* &
      Y(jhe4)*state % rho + screened_rates(k_p28__he4_al24) + &
      screened_rates(k_p28__he4_mg24__weak__wc12) + screened_rates(k_p_p28__he4_si25)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jp28, scratch)

    scratch = (&
      -screened_rates(k_he4_p29__cl33)*Y(jhe4)*state % rho - screened_rates(k_he4_p29__p_s32)* &
      Y(jhe4)*state % rho + screened_rates(k_p29__he4_al25) + &
      screened_rates(k_p_p29__he4_si26)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jp29, scratch)

    scratch = (&
      -screened_rates(k_he4_p30__cl34)*Y(jhe4)*state % rho - screened_rates(k_he4_p30__p_s33)* &
      Y(jhe4)*state % rho + screened_rates(k_p30__he4_al26) + &
      screened_rates(k_p_p30__he4_si27)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jp30, scratch)

    scratch = (&
      -screened_rates(k_he4_p31__cl35)*Y(jhe4)*state % rho - screened_rates(k_he4_p31__p_s34)* &
      Y(jhe4)*state % rho + screened_rates(k_p31__he4_al27) + &
      screened_rates(k_p_p31__he4_si28)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jp31, scratch)

    scratch = (&
      -screened_rates(k_he4_p32__cl36)*Y(jhe4)*state % rho - screened_rates(k_he4_p32__p_s35)* &
      Y(jhe4)*state % rho + screened_rates(k_p32__he4_al28) + &
      screened_rates(k_p_p32__he4_si29)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jp32, scratch)

    scratch = (&
      -screened_rates(k_he4_p33__cl37)*Y(jhe4)*state % rho + screened_rates(k_p_p33__he4_si30)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jp33, scratch)

    scratch = (&
      -screened_rates(k_he4_s28__ar32)*Y(jhe4)*state % rho - screened_rates(k_he4_s28__p_cl31)* &
      Y(jhe4)*state % rho + screened_rates(k_s28__he4_si24) &
       )
    call set_jac_entry(state, jhe4, js28, scratch)

    scratch = (&
      -screened_rates(k_he4_s29__ar33)*Y(jhe4)*state % rho - screened_rates(k_he4_s29__p_cl32)* &
      Y(jhe4)*state % rho + screened_rates(k_p_s29__he4_p26)*Y(jp)*state % rho + &
      screened_rates(k_s29__he4_si25) &
       )
    call set_jac_entry(state, jhe4, js29, scratch)

    scratch = (&
      -screened_rates(k_he4_s30__ar34)*Y(jhe4)*state % rho - screened_rates(k_he4_s30__p_cl33)* &
      Y(jhe4)*state % rho + screened_rates(k_p_s30__he4_p27)*Y(jp)*state % rho + &
      screened_rates(k_s30__he4_si26) &
       )
    call set_jac_entry(state, jhe4, js30, scratch)

    scratch = (&
      -screened_rates(k_he4_s31__ar35)*Y(jhe4)*state % rho - screened_rates(k_he4_s31__p_cl34)* &
      Y(jhe4)*state % rho + screened_rates(k_p_s31__he4_p28)*Y(jp)*state % rho + &
      screened_rates(k_s31__he4_si27) &
       )
    call set_jac_entry(state, jhe4, js31, scratch)

    scratch = (&
      -screened_rates(k_he4_s32__ar36)*Y(jhe4)*state % rho - screened_rates(k_he4_s32__p_cl35)* &
      Y(jhe4)*state % rho + screened_rates(k_p_s32__he4_p29)*Y(jp)*state % rho + &
      screened_rates(k_s32__he4_si28) &
       )
    call set_jac_entry(state, jhe4, js32, scratch)

    scratch = (&
      -screened_rates(k_he4_s33__ar37)*Y(jhe4)*state % rho - screened_rates(k_he4_s33__p_cl36)* &
      Y(jhe4)*state % rho + screened_rates(k_p_s33__he4_p30)*Y(jp)*state % rho + &
      screened_rates(k_s33__he4_si29) &
       )
    call set_jac_entry(state, jhe4, js33, scratch)

    scratch = (&
      -screened_rates(k_he4_s34__ar38)*Y(jhe4)*state % rho - screened_rates(k_he4_s34__p_cl37)* &
      Y(jhe4)*state % rho + screened_rates(k_p_s34__he4_p31)*Y(jp)*state % rho + &
      screened_rates(k_s34__he4_si30) &
       )
    call set_jac_entry(state, jhe4, js34, scratch)

    scratch = (&
      -screened_rates(k_he4_s35__ar39)*Y(jhe4)*state % rho + screened_rates(k_p_s35__he4_p32)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, js35, scratch)

    scratch = (&
      -screened_rates(k_he4_cl29__k33)*Y(jhe4)*state % rho - screened_rates(k_he4_cl29__p_ar32)* &
      Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcl29, scratch)

    scratch = (&
      screened_rates(k_cl30__he4_p26) - screened_rates(k_he4_cl30__k34)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl30__p_ar33)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcl30, scratch)

    scratch = (&
      screened_rates(k_cl31__he4_p27) - screened_rates(k_he4_cl31__k35)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl31__p_ar34)*Y(jhe4)*state % rho + &
      screened_rates(k_p_cl31__he4_s28)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcl31, scratch)

    scratch = (&
      screened_rates(k_cl32__he4_p28) + screened_rates(k_cl32__he4_si28__weak__wc12) - &
      screened_rates(k_he4_cl32__k36)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl32__p_ar35)*Y(jhe4)*state % rho + &
      screened_rates(k_p_cl32__he4_s29)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcl32, scratch)

    scratch = (&
      screened_rates(k_cl33__he4_p29) - screened_rates(k_he4_cl33__k37)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl33__p_ar36)*Y(jhe4)*state % rho + &
      screened_rates(k_p_cl33__he4_s30)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcl33, scratch)

    scratch = (&
      screened_rates(k_cl34__he4_p30) - screened_rates(k_he4_cl34__k38)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl34__p_ar37)*Y(jhe4)*state % rho + &
      screened_rates(k_p_cl34__he4_s31)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcl34, scratch)

    scratch = (&
      screened_rates(k_cl35__he4_p31) - screened_rates(k_he4_cl35__k39)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl35__p_ar38)*Y(jhe4)*state % rho + &
      screened_rates(k_p_cl35__he4_s32)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcl35, scratch)

    scratch = (&
      screened_rates(k_cl36__he4_p32) - screened_rates(k_he4_cl36__k40)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl36__p_ar39)*Y(jhe4)*state % rho + &
      screened_rates(k_p_cl36__he4_s33)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcl36, scratch)

    scratch = (&
      screened_rates(k_cl37__he4_p33) - screened_rates(k_he4_cl37__k41)*Y(jhe4)*state % rho + &
      screened_rates(k_p_cl37__he4_s34)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcl37, scratch)

    scratch = (&
      -screened_rates(k_he4_ar31__ca35)*Y(jhe4)*state % rho - screened_rates(k_he4_ar31__p_k34)* &
      Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jar31, scratch)

    scratch = (&
      screened_rates(k_ar32__he4_s28) - screened_rates(k_he4_ar32__ca36)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar32__p_k35)*Y(jhe4)*state % rho + &
      screened_rates(k_p_ar32__he4_cl29)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jar32, scratch)

    scratch = (&
      screened_rates(k_ar33__he4_s29) - screened_rates(k_he4_ar33__ca37)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar33__p_k36)*Y(jhe4)*state % rho + &
      screened_rates(k_p_ar33__he4_cl30)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jar33, scratch)

    scratch = (&
      screened_rates(k_ar34__he4_s30) - screened_rates(k_he4_ar34__ca38)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar34__p_k37)*Y(jhe4)*state % rho + &
      screened_rates(k_p_ar34__he4_cl31)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jar34, scratch)

    scratch = (&
      screened_rates(k_ar35__he4_s31) - screened_rates(k_he4_ar35__ca39)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar35__p_k38)*Y(jhe4)*state % rho + &
      screened_rates(k_p_ar35__he4_cl32)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jar35, scratch)

    scratch = (&
      screened_rates(k_ar36__he4_s32) - screened_rates(k_he4_ar36__ca40)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar36__p_k39)*Y(jhe4)*state % rho + &
      screened_rates(k_p_ar36__he4_cl33)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jar36, scratch)

    scratch = (&
      screened_rates(k_ar37__he4_s33) - screened_rates(k_he4_ar37__ca41)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar37__p_k40)*Y(jhe4)*state % rho + &
      screened_rates(k_p_ar37__he4_cl34)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jar37, scratch)

    scratch = (&
      screened_rates(k_ar38__he4_s34) - screened_rates(k_he4_ar38__ca42)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar38__p_k41)*Y(jhe4)*state % rho + &
      screened_rates(k_p_ar38__he4_cl35)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jar38, scratch)

    scratch = (&
      screened_rates(k_ar39__he4_s35) - screened_rates(k_he4_ar39__ca43)*Y(jhe4)*state % rho + &
      screened_rates(k_p_ar39__he4_cl36)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jar39, scratch)

    scratch = (&
      -screened_rates(k_he4_k33__p_ca36)*Y(jhe4)*state % rho - screened_rates(k_he4_k33__sc37)* &
      Y(jhe4)*state % rho + screened_rates(k_k33__he4_cl29) &
       )
    call set_jac_entry(state, jhe4, jk33, scratch)

    scratch = (&
      -screened_rates(k_he4_k34__p_ca37)*Y(jhe4)*state % rho - screened_rates(k_he4_k34__sc38)* &
      Y(jhe4)*state % rho + screened_rates(k_k34__he4_cl30) + &
      screened_rates(k_p_k34__he4_ar31)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jk34, scratch)

    scratch = (&
      -screened_rates(k_he4_k35__p_ca38)*Y(jhe4)*state % rho - screened_rates(k_he4_k35__sc39)* &
      Y(jhe4)*state % rho + screened_rates(k_k35__he4_cl31) + &
      screened_rates(k_p_k35__he4_ar32)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jk35, scratch)

    scratch = (&
      -screened_rates(k_he4_k36__p_ca39)*Y(jhe4)*state % rho - screened_rates(k_he4_k36__sc40)* &
      Y(jhe4)*state % rho + screened_rates(k_k36__he4_cl32) + &
      screened_rates(k_k36__he4_s32__weak__wc12) + screened_rates(k_p_k36__he4_ar33)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jk36, scratch)

    scratch = (&
      -screened_rates(k_he4_k37__p_ca40)*Y(jhe4)*state % rho - screened_rates(k_he4_k37__sc41)* &
      Y(jhe4)*state % rho + screened_rates(k_k37__he4_cl33) + &
      screened_rates(k_p_k37__he4_ar34)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jk37, scratch)

    scratch = (&
      -screened_rates(k_he4_k38__p_ca41)*Y(jhe4)*state % rho - screened_rates(k_he4_k38__sc42)* &
      Y(jhe4)*state % rho + screened_rates(k_k38__he4_cl34) + &
      screened_rates(k_p_k38__he4_ar35)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jk38, scratch)

    scratch = (&
      -screened_rates(k_he4_k39__p_ca42)*Y(jhe4)*state % rho - screened_rates(k_he4_k39__sc43)* &
      Y(jhe4)*state % rho + screened_rates(k_k39__he4_cl35) + &
      screened_rates(k_p_k39__he4_ar36)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jk39, scratch)

    scratch = (&
      -screened_rates(k_he4_k40__p_ca43)*Y(jhe4)*state % rho - screened_rates(k_he4_k40__sc44)* &
      Y(jhe4)*state % rho + screened_rates(k_k40__he4_cl36) + &
      screened_rates(k_p_k40__he4_ar37)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jk40, scratch)

    scratch = (&
      -screened_rates(k_he4_k41__p_ca44)*Y(jhe4)*state % rho - screened_rates(k_he4_k41__sc45)* &
      Y(jhe4)*state % rho + screened_rates(k_k41__he4_cl37) + &
      screened_rates(k_p_k41__he4_ar38)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jk41, scratch)

    scratch = (&
      -screened_rates(k_he4_ca34__p_sc37)*Y(jhe4)*state % rho - screened_rates(k_he4_ca34__ti38)* &
      Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jca34, scratch)

    scratch = (&
      screened_rates(k_ca35__he4_ar31) - screened_rates(k_he4_ca35__p_sc38)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca35__ti39)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jca35, scratch)

    scratch = (&
      screened_rates(k_ca36__he4_ar32) - screened_rates(k_he4_ca36__p_p_ca38)*Y(jhe4)*state % rho &
      - screened_rates(k_he4_ca36__p_sc39)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca36__ti40)*Y(jhe4)*state % rho + &
      screened_rates(k_p_ca36__he4_k33)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jca36, scratch)

    scratch = (&
      screened_rates(k_ca37__he4_ar33) - screened_rates(k_he4_ca37__p_sc40)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca37__ti41)*Y(jhe4)*state % rho + &
      screened_rates(k_p_ca37__he4_k34)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jca37, scratch)

    scratch = (&
      screened_rates(k_ca38__he4_ar34) - screened_rates(k_he4_ca38__p_sc41)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca38__ti42)*Y(jhe4)*state % rho + &
      screened_rates(k_p_ca38__he4_k35)*Y(jp)*state % rho + 0.5d0* &
      screened_rates(k_p_p_ca38__he4_ca36)*Y(jp)**2*state % rho**2 &
       )
    call set_jac_entry(state, jhe4, jca38, scratch)

    scratch = (&
      screened_rates(k_ca39__he4_ar35) - screened_rates(k_he4_ca39__p_sc42)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca39__ti43)*Y(jhe4)*state % rho + &
      screened_rates(k_p_ca39__he4_k36)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jca39, scratch)

    scratch = (&
      screened_rates(k_ca40__he4_ar36) - screened_rates(k_he4_ca40__p_sc43)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca40__ti44)*Y(jhe4)*state % rho + &
      screened_rates(k_p_ca40__he4_k37)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jca40, scratch)

    scratch = (&
      screened_rates(k_ca41__he4_ar37) - screened_rates(k_he4_ca41__p_sc44)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca41__ti45)*Y(jhe4)*state % rho + &
      screened_rates(k_p_ca41__he4_k38)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jca41, scratch)

    scratch = (&
      screened_rates(k_ca42__he4_ar38) - screened_rates(k_he4_ca42__p_sc45)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca42__ti46)*Y(jhe4)*state % rho + &
      screened_rates(k_p_ca42__he4_k39)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jca42, scratch)

    scratch = (&
      screened_rates(k_ca43__he4_ar39) - screened_rates(k_he4_ca43__p_sc46)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca43__ti47)*Y(jhe4)*state % rho + &
      screened_rates(k_p_ca43__he4_k40)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jca43, scratch)

    scratch = (&
      -screened_rates(k_he4_ca44__ti48)*Y(jhe4)*state % rho + screened_rates(k_p_ca44__he4_k41)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jca44, scratch)

    scratch = (&
      -screened_rates(k_he4_sc36__p_ti39)*Y(jhe4)*state % rho - screened_rates(k_he4_sc36__v40)* &
      Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jsc36, scratch)

    scratch = (&
      -screened_rates(k_he4_sc37__p_ti40)*Y(jhe4)*state % rho - screened_rates(k_he4_sc37__v41)* &
      Y(jhe4)*state % rho + screened_rates(k_p_sc37__he4_ca34)*Y(jp)*state % rho + &
      screened_rates(k_sc37__he4_k33) &
       )
    call set_jac_entry(state, jhe4, jsc37, scratch)

    scratch = (&
      -screened_rates(k_he4_sc38__p_ti41)*Y(jhe4)*state % rho - screened_rates(k_he4_sc38__v42)* &
      Y(jhe4)*state % rho + screened_rates(k_p_sc38__he4_ca35)*Y(jp)*state % rho + &
      screened_rates(k_sc38__he4_k34) &
       )
    call set_jac_entry(state, jhe4, jsc38, scratch)

    scratch = (&
      -screened_rates(k_he4_sc39__p_ti42)*Y(jhe4)*state % rho - screened_rates(k_he4_sc39__v43)* &
      Y(jhe4)*state % rho + screened_rates(k_p_sc39__he4_ca36)*Y(jp)*state % rho + &
      screened_rates(k_sc39__he4_k35) &
       )
    call set_jac_entry(state, jhe4, jsc39, scratch)

    scratch = (&
      -screened_rates(k_he4_sc40__p_ti43)*Y(jhe4)*state % rho - screened_rates(k_he4_sc40__v44)* &
      Y(jhe4)*state % rho + screened_rates(k_p_sc40__he4_ca37)*Y(jp)*state % rho + &
      screened_rates(k_sc40__he4_ar36__weak__wc12) + screened_rates(k_sc40__he4_k36) &
       )
    call set_jac_entry(state, jhe4, jsc40, scratch)

    scratch = (&
      -screened_rates(k_he4_sc41__p_ti44)*Y(jhe4)*state % rho - screened_rates(k_he4_sc41__v45)* &
      Y(jhe4)*state % rho + screened_rates(k_p_sc41__he4_ca38)*Y(jp)*state % rho + &
      screened_rates(k_sc41__he4_k37) &
       )
    call set_jac_entry(state, jhe4, jsc41, scratch)

    scratch = (&
      -screened_rates(k_he4_sc42__p_ti45)*Y(jhe4)*state % rho - screened_rates(k_he4_sc42__v46)* &
      Y(jhe4)*state % rho + screened_rates(k_p_sc42__he4_ca39)*Y(jp)*state % rho + &
      screened_rates(k_sc42__he4_k38) &
       )
    call set_jac_entry(state, jhe4, jsc42, scratch)

    scratch = (&
      -screened_rates(k_he4_sc43__p_ti46)*Y(jhe4)*state % rho - screened_rates(k_he4_sc43__v47)* &
      Y(jhe4)*state % rho + screened_rates(k_p_sc43__he4_ca40)*Y(jp)*state % rho + &
      screened_rates(k_sc43__he4_k39) &
       )
    call set_jac_entry(state, jhe4, jsc43, scratch)

    scratch = (&
      -screened_rates(k_he4_sc44__p_ti47)*Y(jhe4)*state % rho - screened_rates(k_he4_sc44__v48)* &
      Y(jhe4)*state % rho + screened_rates(k_p_sc44__he4_ca41)*Y(jp)*state % rho + &
      screened_rates(k_sc44__he4_k40) &
       )
    call set_jac_entry(state, jhe4, jsc44, scratch)

    scratch = (&
      -screened_rates(k_he4_sc45__p_ti48)*Y(jhe4)*state % rho - screened_rates(k_he4_sc45__v49)* &
      Y(jhe4)*state % rho + screened_rates(k_p_sc45__he4_ca42)*Y(jp)*state % rho + &
      screened_rates(k_sc45__he4_k41) &
       )
    call set_jac_entry(state, jhe4, jsc45, scratch)

    scratch = (&
      -screened_rates(k_he4_sc46__v50)*Y(jhe4)*state % rho + screened_rates(k_p_sc46__he4_ca43)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jsc46, scratch)

    scratch = (&
      -screened_rates(k_he4_ti38__cr42)*Y(jhe4)*state % rho - screened_rates(k_he4_ti38__p_v41)* &
      Y(jhe4)*state % rho + screened_rates(k_ti38__he4_ca34) &
       )
    call set_jac_entry(state, jhe4, jti38, scratch)

    scratch = (&
      -screened_rates(k_he4_ti39__cr43)*Y(jhe4)*state % rho - screened_rates(k_he4_ti39__p_v42)* &
      Y(jhe4)*state % rho + screened_rates(k_p_ti39__he4_sc36)*Y(jp)*state % rho + &
      screened_rates(k_ti39__he4_ca35) &
       )
    call set_jac_entry(state, jhe4, jti39, scratch)

    scratch = (&
      -screened_rates(k_he4_ti40__cr44)*Y(jhe4)*state % rho - screened_rates(k_he4_ti40__p_v43)* &
      Y(jhe4)*state % rho + screened_rates(k_p_ti40__he4_sc37)*Y(jp)*state % rho + &
      screened_rates(k_ti40__he4_ca36) &
       )
    call set_jac_entry(state, jhe4, jti40, scratch)

    scratch = (&
      -screened_rates(k_he4_ti41__cr45)*Y(jhe4)*state % rho - screened_rates(k_he4_ti41__p_v44)* &
      Y(jhe4)*state % rho + screened_rates(k_p_ti41__he4_sc38)*Y(jp)*state % rho + &
      screened_rates(k_ti41__he4_ca37) &
       )
    call set_jac_entry(state, jhe4, jti41, scratch)

    scratch = (&
      -screened_rates(k_he4_ti42__cr46)*Y(jhe4)*state % rho - screened_rates(k_he4_ti42__p_v45)* &
      Y(jhe4)*state % rho + screened_rates(k_p_ti42__he4_sc39)*Y(jp)*state % rho + &
      screened_rates(k_ti42__he4_ca38) &
       )
    call set_jac_entry(state, jhe4, jti42, scratch)

    scratch = (&
      -screened_rates(k_he4_ti43__cr47)*Y(jhe4)*state % rho - screened_rates(k_he4_ti43__p_v46)* &
      Y(jhe4)*state % rho + screened_rates(k_p_ti43__he4_sc40)*Y(jp)*state % rho + &
      screened_rates(k_ti43__he4_ca39) &
       )
    call set_jac_entry(state, jhe4, jti43, scratch)

    scratch = (&
      -screened_rates(k_he4_ti44__cr48)*Y(jhe4)*state % rho - screened_rates(k_he4_ti44__p_v47)* &
      Y(jhe4)*state % rho + screened_rates(k_p_ti44__he4_sc41)*Y(jp)*state % rho + &
      screened_rates(k_ti44__he4_ca40) &
       )
    call set_jac_entry(state, jhe4, jti44, scratch)

    scratch = (&
      -screened_rates(k_he4_ti45__cr49)*Y(jhe4)*state % rho - screened_rates(k_he4_ti45__p_v48)* &
      Y(jhe4)*state % rho + screened_rates(k_p_ti45__he4_sc42)*Y(jp)*state % rho + &
      screened_rates(k_ti45__he4_ca41) &
       )
    call set_jac_entry(state, jhe4, jti45, scratch)

    scratch = (&
      -screened_rates(k_he4_ti46__cr50)*Y(jhe4)*state % rho - screened_rates(k_he4_ti46__p_v49)* &
      Y(jhe4)*state % rho + screened_rates(k_p_ti46__he4_sc43)*Y(jp)*state % rho + &
      screened_rates(k_ti46__he4_ca42) &
       )
    call set_jac_entry(state, jhe4, jti46, scratch)

    scratch = (&
      -screened_rates(k_he4_ti47__cr51)*Y(jhe4)*state % rho - screened_rates(k_he4_ti47__p_v50)* &
      Y(jhe4)*state % rho + screened_rates(k_p_ti47__he4_sc44)*Y(jp)*state % rho + &
      screened_rates(k_ti47__he4_ca43) &
       )
    call set_jac_entry(state, jhe4, jti47, scratch)

    scratch = (&
      -screened_rates(k_he4_ti48__cr52)*Y(jhe4)*state % rho + screened_rates(k_p_ti48__he4_sc45)* &
      Y(jp)*state % rho + screened_rates(k_ti48__he4_ca44) &
       )
    call set_jac_entry(state, jhe4, jti48, scratch)

    scratch = (&
      -screened_rates(k_he4_v40__mn44)*Y(jhe4)*state % rho - screened_rates(k_he4_v40__p_cr43)* &
      Y(jhe4)*state % rho + screened_rates(k_v40__he4_sc36) &
       )
    call set_jac_entry(state, jhe4, jv40, scratch)

    scratch = (&
      -screened_rates(k_he4_v41__mn45)*Y(jhe4)*state % rho - screened_rates(k_he4_v41__p_cr44)* &
      Y(jhe4)*state % rho + screened_rates(k_p_v41__he4_ti38)*Y(jp)*state % rho + &
      screened_rates(k_v41__he4_sc37) &
       )
    call set_jac_entry(state, jhe4, jv41, scratch)

    scratch = (&
      -screened_rates(k_he4_v42__mn46)*Y(jhe4)*state % rho - screened_rates(k_he4_v42__p_cr45)* &
      Y(jhe4)*state % rho + screened_rates(k_p_v42__he4_ti39)*Y(jp)*state % rho + &
      screened_rates(k_v42__he4_sc38) &
       )
    call set_jac_entry(state, jhe4, jv42, scratch)

    scratch = (&
      -screened_rates(k_he4_v43__mn47)*Y(jhe4)*state % rho - screened_rates(k_he4_v43__p_cr46)* &
      Y(jhe4)*state % rho + screened_rates(k_p_v43__he4_ti40)*Y(jp)*state % rho + &
      screened_rates(k_v43__he4_sc39) &
       )
    call set_jac_entry(state, jhe4, jv43, scratch)

    scratch = (&
      -screened_rates(k_he4_v44__mn48)*Y(jhe4)*state % rho - screened_rates(k_he4_v44__p_cr47)* &
      Y(jhe4)*state % rho + screened_rates(k_p_v44__he4_ti41)*Y(jp)*state % rho + &
      screened_rates(k_v44__he4_sc40) &
       )
    call set_jac_entry(state, jhe4, jv44, scratch)

    scratch = (&
      -screened_rates(k_he4_v45__mn49)*Y(jhe4)*state % rho - screened_rates(k_he4_v45__p_cr48)* &
      Y(jhe4)*state % rho + screened_rates(k_p_v45__he4_ti42)*Y(jp)*state % rho + &
      screened_rates(k_v45__he4_sc41) &
       )
    call set_jac_entry(state, jhe4, jv45, scratch)

    scratch = (&
      -screened_rates(k_he4_v46__mn50)*Y(jhe4)*state % rho - screened_rates(k_he4_v46__p_cr49)* &
      Y(jhe4)*state % rho + screened_rates(k_p_v46__he4_ti43)*Y(jp)*state % rho + &
      screened_rates(k_v46__he4_sc42) &
       )
    call set_jac_entry(state, jhe4, jv46, scratch)

    scratch = (&
      -screened_rates(k_he4_v47__mn51)*Y(jhe4)*state % rho - screened_rates(k_he4_v47__p_cr50)* &
      Y(jhe4)*state % rho + screened_rates(k_p_v47__he4_ti44)*Y(jp)*state % rho + &
      screened_rates(k_v47__he4_sc43) &
       )
    call set_jac_entry(state, jhe4, jv47, scratch)

    scratch = (&
      -screened_rates(k_he4_v48__mn52)*Y(jhe4)*state % rho - screened_rates(k_he4_v48__p_cr51)* &
      Y(jhe4)*state % rho + screened_rates(k_p_v48__he4_ti45)*Y(jp)*state % rho + &
      screened_rates(k_v48__he4_sc44) &
       )
    call set_jac_entry(state, jhe4, jv48, scratch)

    scratch = (&
      -screened_rates(k_he4_v49__mn53)*Y(jhe4)*state % rho - screened_rates(k_he4_v49__p_cr52)* &
      Y(jhe4)*state % rho + screened_rates(k_p_v49__he4_ti46)*Y(jp)*state % rho + &
      screened_rates(k_v49__he4_sc45) &
       )
    call set_jac_entry(state, jhe4, jv49, scratch)

    scratch = (&
      -screened_rates(k_he4_v50__mn54)*Y(jhe4)*state % rho + screened_rates(k_p_v50__he4_ti47)* &
      Y(jp)*state % rho + screened_rates(k_v50__he4_sc46) &
       )
    call set_jac_entry(state, jhe4, jv50, scratch)

    scratch = (&
      screened_rates(k_cr42__he4_ti38) - screened_rates(k_he4_cr42__fe46)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr42__p_mn45)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcr42, scratch)

    scratch = (&
      screened_rates(k_cr43__he4_ti39) - screened_rates(k_he4_cr43__fe47)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr43__p_mn46)*Y(jhe4)*state % rho + &
      screened_rates(k_p_cr43__he4_v40)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcr43, scratch)

    scratch = (&
      screened_rates(k_cr44__he4_ti40) - screened_rates(k_he4_cr44__fe48)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr44__p_mn47)*Y(jhe4)*state % rho + &
      screened_rates(k_p_cr44__he4_v41)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcr44, scratch)

    scratch = (&
      screened_rates(k_cr45__he4_ti41) - screened_rates(k_he4_cr45__fe49)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr45__p_mn48)*Y(jhe4)*state % rho + &
      screened_rates(k_p_cr45__he4_v42)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcr45, scratch)

    scratch = (&
      screened_rates(k_cr46__he4_ti42) - screened_rates(k_he4_cr46__fe50)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr46__p_mn49)*Y(jhe4)*state % rho + &
      screened_rates(k_p_cr46__he4_v43)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcr46, scratch)

    scratch = (&
      screened_rates(k_cr47__he4_ti43) - screened_rates(k_he4_cr47__fe51)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr47__p_mn50)*Y(jhe4)*state % rho + &
      screened_rates(k_p_cr47__he4_v44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcr47, scratch)

    scratch = (&
      screened_rates(k_cr48__he4_ti44) - screened_rates(k_he4_cr48__fe52)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr48__p_mn51)*Y(jhe4)*state % rho + &
      screened_rates(k_p_cr48__he4_v45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcr48, scratch)

    scratch = (&
      screened_rates(k_cr49__he4_ti45) - screened_rates(k_he4_cr49__fe53)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr49__p_mn52)*Y(jhe4)*state % rho + &
      screened_rates(k_p_cr49__he4_v46)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcr49, scratch)

    scratch = (&
      screened_rates(k_cr50__he4_ti46) - screened_rates(k_he4_cr50__fe54)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr50__p_mn53)*Y(jhe4)*state % rho + &
      screened_rates(k_p_cr50__he4_v47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcr50, scratch)

    scratch = (&
      screened_rates(k_cr51__he4_ti47) - screened_rates(k_he4_cr51__fe55)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr51__p_mn54)*Y(jhe4)*state % rho + &
      screened_rates(k_p_cr51__he4_v48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcr51, scratch)

    scratch = (&
      screened_rates(k_cr52__he4_ti48) - screened_rates(k_he4_cr52__fe56)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr52__p_mn55)*Y(jhe4)*state % rho + &
      screened_rates(k_p_cr52__he4_v49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jcr52, scratch)

    scratch = (&
      -screened_rates(k_he4_mn44__co48)*Y(jhe4)*state % rho - screened_rates(k_he4_mn44__p_fe47)* &
      Y(jhe4)*state % rho + screened_rates(k_mn44__he4_v40) &
       )
    call set_jac_entry(state, jhe4, jmn44, scratch)

    scratch = (&
      -screened_rates(k_he4_mn45__co49)*Y(jhe4)*state % rho - screened_rates(k_he4_mn45__p_fe48)* &
      Y(jhe4)*state % rho + screened_rates(k_mn45__he4_v41) + &
      screened_rates(k_p_mn45__he4_cr42)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jmn45, scratch)

    scratch = (&
      -screened_rates(k_he4_mn46__co50)*Y(jhe4)*state % rho - screened_rates(k_he4_mn46__p_fe49)* &
      Y(jhe4)*state % rho + screened_rates(k_mn46__he4_v42) + &
      screened_rates(k_p_mn46__he4_cr43)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jmn46, scratch)

    scratch = (&
      -screened_rates(k_he4_mn47__co51)*Y(jhe4)*state % rho - screened_rates(k_he4_mn47__p_fe50)* &
      Y(jhe4)*state % rho + screened_rates(k_mn47__he4_v43) + &
      screened_rates(k_p_mn47__he4_cr44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jmn47, scratch)

    scratch = (&
      -screened_rates(k_he4_mn48__co52)*Y(jhe4)*state % rho - screened_rates(k_he4_mn48__p_fe51)* &
      Y(jhe4)*state % rho + screened_rates(k_mn48__he4_v44) + &
      screened_rates(k_p_mn48__he4_cr45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jmn48, scratch)

    scratch = (&
      -screened_rates(k_he4_mn49__co53)*Y(jhe4)*state % rho - screened_rates(k_he4_mn49__p_fe52)* &
      Y(jhe4)*state % rho + screened_rates(k_mn49__he4_v45) + &
      screened_rates(k_p_mn49__he4_cr46)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jmn49, scratch)

    scratch = (&
      -screened_rates(k_he4_mn50__co54)*Y(jhe4)*state % rho - screened_rates(k_he4_mn50__p_fe53)* &
      Y(jhe4)*state % rho + screened_rates(k_mn50__he4_v46) + &
      screened_rates(k_p_mn50__he4_cr47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jmn50, scratch)

    scratch = (&
      -screened_rates(k_he4_mn51__co55)*Y(jhe4)*state % rho - screened_rates(k_he4_mn51__p_fe54)* &
      Y(jhe4)*state % rho + screened_rates(k_mn51__he4_v47) + &
      screened_rates(k_p_mn51__he4_cr48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jmn51, scratch)

    scratch = (&
      -screened_rates(k_he4_mn52__co56)*Y(jhe4)*state % rho - screened_rates(k_he4_mn52__p_fe55)* &
      Y(jhe4)*state % rho + screened_rates(k_mn52__he4_v48) + &
      screened_rates(k_p_mn52__he4_cr49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jmn52, scratch)

    scratch = (&
      -screened_rates(k_he4_mn53__p_fe56)*Y(jhe4)*state % rho + screened_rates(k_mn53__he4_v49) + &
      screened_rates(k_p_mn53__he4_cr50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jmn53, scratch)

    scratch = (&
      screened_rates(k_mn54__he4_v50) + screened_rates(k_p_mn54__he4_cr51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jmn54, scratch)

    scratch = (&
      screened_rates(k_p_mn55__he4_cr52)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jmn55, scratch)

    scratch = (&
      -screened_rates(k_he4_fe45__ni49)*Y(jhe4)*state % rho - screened_rates(k_he4_fe45__p_co48)* &
      Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jfe45, scratch)

    scratch = (&
      screened_rates(k_fe46__he4_cr42) - screened_rates(k_he4_fe46__ni50)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe46__p_co49)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jfe46, scratch)

    scratch = (&
      screened_rates(k_fe47__he4_cr43) - screened_rates(k_he4_fe47__ni51)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe47__p_co50)*Y(jhe4)*state % rho + &
      screened_rates(k_p_fe47__he4_mn44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jfe47, scratch)

    scratch = (&
      screened_rates(k_fe48__he4_cr44) - screened_rates(k_he4_fe48__ni52)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe48__p_co51)*Y(jhe4)*state % rho + &
      screened_rates(k_p_fe48__he4_mn45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jfe48, scratch)

    scratch = (&
      screened_rates(k_fe49__he4_cr45) - screened_rates(k_he4_fe49__ni53)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe49__p_co52)*Y(jhe4)*state % rho + &
      screened_rates(k_p_fe49__he4_mn46)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jfe49, scratch)

    scratch = (&
      screened_rates(k_fe50__he4_cr46) - screened_rates(k_he4_fe50__ni54)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe50__p_co53)*Y(jhe4)*state % rho + &
      screened_rates(k_p_fe50__he4_mn47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jfe50, scratch)

    scratch = (&
      screened_rates(k_fe51__he4_cr47) - screened_rates(k_he4_fe51__ni55)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe51__p_co54)*Y(jhe4)*state % rho + &
      screened_rates(k_p_fe51__he4_mn48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jfe51, scratch)

    scratch = (&
      screened_rates(k_fe52__he4_cr48) - screened_rates(k_he4_fe52__ni56)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe52__p_co55)*Y(jhe4)*state % rho + &
      screened_rates(k_p_fe52__he4_mn49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jfe52, scratch)

    scratch = (&
      screened_rates(k_fe53__he4_cr49) - screened_rates(k_he4_fe53__p_co56)*Y(jhe4)*state % rho + &
      screened_rates(k_p_fe53__he4_mn50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jfe53, scratch)

    scratch = (&
      screened_rates(k_fe54__he4_cr50) + screened_rates(k_p_fe54__he4_mn51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jfe54, scratch)

    scratch = (&
      screened_rates(k_fe55__he4_cr51) + screened_rates(k_p_fe55__he4_mn52)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jfe55, scratch)

    scratch = (&
      screened_rates(k_fe56__he4_cr52) + screened_rates(k_p_fe56__he4_mn53)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jfe56, scratch)

    scratch = (&
      -screened_rates(k_he4_co47__p_ni50)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jco47, scratch)

    scratch = (&
      screened_rates(k_co48__he4_mn44) - screened_rates(k_he4_co48__p_ni51)*Y(jhe4)*state % rho + &
      screened_rates(k_p_co48__he4_fe45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jco48, scratch)

    scratch = (&
      screened_rates(k_co49__he4_mn45) - screened_rates(k_he4_co49__p_ni52)*Y(jhe4)*state % rho + &
      screened_rates(k_p_co49__he4_fe46)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jco49, scratch)

    scratch = (&
      screened_rates(k_co50__he4_mn46) - screened_rates(k_he4_co50__p_ni53)*Y(jhe4)*state % rho + &
      screened_rates(k_p_co50__he4_fe47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jco50, scratch)

    scratch = (&
      screened_rates(k_co51__he4_mn47) - screened_rates(k_he4_co51__p_ni54)*Y(jhe4)*state % rho + &
      screened_rates(k_p_co51__he4_fe48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jco51, scratch)

    scratch = (&
      screened_rates(k_co52__he4_mn48) - screened_rates(k_he4_co52__p_ni55)*Y(jhe4)*state % rho + &
      screened_rates(k_p_co52__he4_fe49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jco52, scratch)

    scratch = (&
      screened_rates(k_co53__he4_mn49) - screened_rates(k_he4_co53__p_ni56)*Y(jhe4)*state % rho + &
      screened_rates(k_p_co53__he4_fe50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jco53, scratch)

    scratch = (&
      screened_rates(k_co54__he4_mn50) + screened_rates(k_p_co54__he4_fe51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jco54, scratch)

    scratch = (&
      screened_rates(k_co55__he4_mn51) + screened_rates(k_p_co55__he4_fe52)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jco55, scratch)

    scratch = (&
      screened_rates(k_co56__he4_mn52) + screened_rates(k_p_co56__he4_fe53)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jco56, scratch)

    scratch = (&
      screened_rates(k_ni49__he4_fe45) &
       )
    call set_jac_entry(state, jhe4, jni49, scratch)

    scratch = (&
      screened_rates(k_ni50__he4_fe46) + screened_rates(k_p_ni50__he4_co47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jni50, scratch)

    scratch = (&
      screened_rates(k_ni51__he4_fe47) + screened_rates(k_p_ni51__he4_co48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jni51, scratch)

    scratch = (&
      screened_rates(k_ni52__he4_fe48) + screened_rates(k_p_ni52__he4_co49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jni52, scratch)

    scratch = (&
      screened_rates(k_ni53__he4_fe49) + screened_rates(k_p_ni53__he4_co50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jni53, scratch)

    scratch = (&
      screened_rates(k_ni54__he4_fe50) + screened_rates(k_p_ni54__he4_co51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jni54, scratch)

    scratch = (&
      screened_rates(k_ni55__he4_fe51) + screened_rates(k_p_ni55__he4_co52)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jni55, scratch)

    scratch = (&
      screened_rates(k_ni56__he4_fe52) + screened_rates(k_p_ni56__he4_co53)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jhe4, jni56, scratch)

    scratch = (&
      -screened_rates(k_p_li7__he4_he4)*Y(jli7)*state % rho &
       )
    call set_jac_entry(state, jli7, jp, scratch)

    scratch = (&
      1.0d0*screened_rates(k_he4_he4__p_li7)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jli7, jhe4, scratch)

    scratch = (&
      -screened_rates(k_p_li7__he4_he4)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jli7, jli7, scratch)

    scratch = (&
      screened_rates(k_be7__li7__weak__electron_capture)*state % rho*state % y_e &
       )
    call set_jac_entry(state, jli7, jbe7, scratch)

    scratch = (&
      -screened_rates(k_p_be7__b8)*Y(jbe7)*state % rho + 0.5d0*screened_rates(k_p_he4_he4__d_be7) &
      *Y(jhe4)**2*state % rho**2 + 0.5d0*screened_rates(k_p_p_he4_he4__he3_be7)* &
      Y(jhe4)**2*Y(jp)*state % rho**3 &
       )
    call set_jac_entry(state, jbe7, jp, scratch)

    scratch = (&
      -screened_rates(k_d_be7__p_he4_he4)*Y(jbe7)*state % rho &
       )
    call set_jac_entry(state, jbe7, jd, scratch)

    scratch = (&
      -screened_rates(k_he3_be7__p_p_he4_he4)*Y(jbe7)*state % rho + &
      screened_rates(k_he4_he3__be7)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jbe7, jhe3, scratch)

    scratch = (&
      screened_rates(k_he4_he3__be7)*Y(jhe3)*state % rho + 1.0d0* &
      screened_rates(k_p_he4_he4__d_be7)*Y(jhe4)*Y(jp)*state % rho**2 + 0.5d0* &
      screened_rates(k_p_p_he4_he4__he3_be7)*Y(jhe4)*Y(jp)**2*state % rho**3 &
       )
    call set_jac_entry(state, jbe7, jhe4, scratch)

    scratch = (&
      -screened_rates(k_be7__he4_he3) - screened_rates(k_be7__li7__weak__electron_capture)* &
      state % rho*state % y_e - screened_rates(k_d_be7__p_he4_he4)*Y(jd)*state % rho - &
      screened_rates(k_he3_be7__p_p_he4_he4)*Y(jhe3)*state % rho - &
      screened_rates(k_p_be7__b8)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jbe7, jbe7, scratch)

    scratch = (&
      screened_rates(k_b8__p_be7) &
       )
    call set_jac_entry(state, jbe7, jb8, scratch)

    scratch = (&
      screened_rates(k_b8__be8__weak__wc17) &
       )
    call set_jac_entry(state, jbe8, jb8, scratch)

    scratch = (&
      screened_rates(k_p_be7__b8)*Y(jbe7)*state % rho &
       )
    call set_jac_entry(state, jb8, jp, scratch)

    scratch = (&
      screened_rates(k_p_be7__b8)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jb8, jbe7, scratch)

    scratch = (&
      -screened_rates(k_b8__be8__weak__wc17) - screened_rates(k_b8__he4_he4__weak__wc12) - &
      screened_rates(k_b8__p_be7) &
       )
    call set_jac_entry(state, jb8, jb8, scratch)

    scratch = (&
      -screened_rates(k_p_c12__n13)*Y(jc12)*state % rho + screened_rates(k_p_n15__he4_c12)* &
      Y(jn15)*state % rho &
       )
    call set_jac_entry(state, jc12, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho - screened_rates(k_he4_c12__p_n15)* &
      Y(jc12)*state % rho + 0.5d0*screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**2 &
      *state % rho**2 + 2.0d0*screened_rates(k_he4_ne20__c12_c12)*Y(jne20)* &
      state % rho &
       )
    call set_jac_entry(state, jc12, jhe4, scratch)

    scratch = (&
      -screened_rates(k_c12__he4_he4_he4) - 2.0d0*screened_rates(k_c12_c12__he4_ne20)*Y(jc12)* &
      state % rho - screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_c12__p_n15)*Y(jhe4)*state % rho - screened_rates(k_p_c12__n13)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jc12, jc12, scratch)

    scratch = (&
      screened_rates(k_n12__c12__weak__wc12) &
       )
    call set_jac_entry(state, jc12, jn12, scratch)

    scratch = (&
      screened_rates(k_n13__p_c12) &
       )
    call set_jac_entry(state, jc12, jn13, scratch)

    scratch = (&
      screened_rates(k_p_n15__he4_c12)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jc12, jn15, scratch)

    scratch = (&
      screened_rates(k_o16__he4_c12) &
       )
    call set_jac_entry(state, jc12, jo16, scratch)

    scratch = (&
      2.0d0*screened_rates(k_he4_ne20__c12_c12)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jc12, jne20, scratch)

    scratch = (&
      -screened_rates(k_p_c13__n14)*Y(jc13)*state % rho &
       )
    call set_jac_entry(state, jc13, jp, scratch)

    scratch = (&
      -screened_rates(k_p_c13__n14)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jc13, jc13, scratch)

    scratch = (&
      screened_rates(k_n13__c13__weak__wc12) &
       )
    call set_jac_entry(state, jc13, jn13, scratch)

    scratch = (&
      screened_rates(k_n14__p_c13) &
       )
    call set_jac_entry(state, jc13, jn14, scratch)

    scratch = (&
      screened_rates(k_p_o15__he4_n12)*Y(jo15)*state % rho &
       )
    call set_jac_entry(state, jn12, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_n12__p_o15)*Y(jn12)*state % rho &
       )
    call set_jac_entry(state, jn12, jhe4, scratch)

    scratch = (&
      -screened_rates(k_he4_n12__p_o15)*Y(jhe4)*state % rho - &
      screened_rates(k_n12__c12__weak__wc12) &
       )
    call set_jac_entry(state, jn12, jn12, scratch)

    scratch = (&
      screened_rates(k_p_o15__he4_n12)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jn12, jo15, scratch)

    scratch = (&
      screened_rates(k_p_c12__n13)*Y(jc12)*state % rho - screened_rates(k_p_n13__o14)*Y(jn13)* &
      state % rho + screened_rates(k_p_o16__he4_n13)*Y(jo16)*state % rho &
       )
    call set_jac_entry(state, jn13, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_n13__p_o16)*Y(jn13)*state % rho &
       )
    call set_jac_entry(state, jn13, jhe4, scratch)

    scratch = (&
      screened_rates(k_p_c12__n13)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jn13, jc12, scratch)

    scratch = (&
      -screened_rates(k_he4_n13__p_o16)*Y(jhe4)*state % rho - &
      screened_rates(k_n13__c13__weak__wc12) - screened_rates(k_n13__p_c12) - &
      screened_rates(k_p_n13__o14)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jn13, jn13, scratch)

    scratch = (&
      screened_rates(k_o14__p_n13) &
       )
    call set_jac_entry(state, jn13, jo14, scratch)

    scratch = (&
      screened_rates(k_p_o16__he4_n13)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jn13, jo16, scratch)

    scratch = (&
      screened_rates(k_p_c13__n14)*Y(jc13)*state % rho - screened_rates(k_p_n14__o15)*Y(jn14)* &
      state % rho + screened_rates(k_p_o17__he4_n14)*Y(jo17)*state % rho &
       )
    call set_jac_entry(state, jn14, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_n14__f18)*Y(jn14)*state % rho - screened_rates(k_he4_n14__p_o17)* &
      Y(jn14)*state % rho &
       )
    call set_jac_entry(state, jn14, jhe4, scratch)

    scratch = (&
      screened_rates(k_p_c13__n14)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jn14, jc13, scratch)

    scratch = (&
      -screened_rates(k_he4_n14__f18)*Y(jhe4)*state % rho - screened_rates(k_he4_n14__p_o17)* &
      Y(jhe4)*state % rho - screened_rates(k_n14__p_c13) - screened_rates(k_p_n14__o15)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jn14, jn14, scratch)

    scratch = (&
      screened_rates(k_o14__n14__weak__wc12) &
       )
    call set_jac_entry(state, jn14, jo14, scratch)

    scratch = (&
      screened_rates(k_o15__p_n14) &
       )
    call set_jac_entry(state, jn14, jo15, scratch)

    scratch = (&
      screened_rates(k_p_o17__he4_n14)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jn14, jo17, scratch)

    scratch = (&
      screened_rates(k_f18__he4_n14) &
       )
    call set_jac_entry(state, jn14, jf18, scratch)

    scratch = (&
      -screened_rates(k_p_n15__he4_c12)*Y(jn15)*state % rho - screened_rates(k_p_n15__o16)* &
      Y(jn15)*state % rho + screened_rates(k_p_o18__he4_n15)*Y(jo18)*state % rho &
       )
    call set_jac_entry(state, jn15, jp, scratch)

    scratch = (&
      screened_rates(k_he4_c12__p_n15)*Y(jc12)*state % rho - screened_rates(k_he4_n15__f19)* &
      Y(jn15)*state % rho - screened_rates(k_he4_n15__p_o18)*Y(jn15)*state % rho &
       )
    call set_jac_entry(state, jn15, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_c12__p_n15)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jn15, jc12, scratch)

    scratch = (&
      -screened_rates(k_he4_n15__f19)*Y(jhe4)*state % rho - screened_rates(k_he4_n15__p_o18)* &
      Y(jhe4)*state % rho - screened_rates(k_p_n15__he4_c12)*Y(jp)*state % rho - &
      screened_rates(k_p_n15__o16)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jn15, jn15, scratch)

    scratch = (&
      screened_rates(k_o15__n15__weak__wc12) &
       )
    call set_jac_entry(state, jn15, jo15, scratch)

    scratch = (&
      screened_rates(k_o16__p_n15) &
       )
    call set_jac_entry(state, jn15, jo16, scratch)

    scratch = (&
      screened_rates(k_p_o18__he4_n15)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jn15, jo18, scratch)

    scratch = (&
      screened_rates(k_f19__he4_n15) &
       )
    call set_jac_entry(state, jn15, jf19, scratch)

    scratch = (&
      screened_rates(k_p_f17__he4_o14)*Y(jf17)*state % rho + screened_rates(k_p_n13__o14)* &
      Y(jn13)*state % rho &
       )
    call set_jac_entry(state, jo14, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_o14__ne18)*Y(jo14)*state % rho - screened_rates(k_he4_o14__p_f17)* &
      Y(jo14)*state % rho &
       )
    call set_jac_entry(state, jo14, jhe4, scratch)

    scratch = (&
      screened_rates(k_p_n13__o14)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jo14, jn13, scratch)

    scratch = (&
      -screened_rates(k_he4_o14__ne18)*Y(jhe4)*state % rho - screened_rates(k_he4_o14__p_f17)* &
      Y(jhe4)*state % rho - screened_rates(k_o14__n14__weak__wc12) - &
      screened_rates(k_o14__p_n13) &
       )
    call set_jac_entry(state, jo14, jo14, scratch)

    scratch = (&
      screened_rates(k_p_f17__he4_o14)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jo14, jf17, scratch)

    scratch = (&
      screened_rates(k_ne18__he4_o14) &
       )
    call set_jac_entry(state, jo14, jne18, scratch)

    scratch = (&
      screened_rates(k_p_f18__he4_o15)*Y(jf18)*state % rho + screened_rates(k_p_n14__o15)* &
      Y(jn14)*state % rho - screened_rates(k_p_o15__f16)*Y(jo15)*state % rho - &
      screened_rates(k_p_o15__he4_n12)*Y(jo15)*state % rho &
       )
    call set_jac_entry(state, jo15, jp, scratch)

    scratch = (&
      screened_rates(k_he4_n12__p_o15)*Y(jn12)*state % rho - screened_rates(k_he4_o15__ne19)* &
      Y(jo15)*state % rho - screened_rates(k_he4_o15__p_f18)*Y(jo15)*state % rho &
       )
    call set_jac_entry(state, jo15, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_n12__p_o15)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jo15, jn12, scratch)

    scratch = (&
      screened_rates(k_p_n14__o15)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jo15, jn14, scratch)

    scratch = (&
      -screened_rates(k_he4_o15__ne19)*Y(jhe4)*state % rho - screened_rates(k_he4_o15__p_f18)* &
      Y(jhe4)*state % rho - screened_rates(k_o15__n15__weak__wc12) - &
      screened_rates(k_o15__p_n14) - screened_rates(k_p_o15__f16)*Y(jp)*state % rho - &
      screened_rates(k_p_o15__he4_n12)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jo15, jo15, scratch)

    scratch = (&
      screened_rates(k_f16__p_o15) &
       )
    call set_jac_entry(state, jo15, jf16, scratch)

    scratch = (&
      screened_rates(k_p_f18__he4_o15)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jo15, jf18, scratch)

    scratch = (&
      screened_rates(k_ne19__he4_o15) &
       )
    call set_jac_entry(state, jo15, jne19, scratch)

    scratch = (&
      screened_rates(k_p_f19__he4_o16)*Y(jf19)*state % rho + screened_rates(k_p_n15__o16)* &
      Y(jn15)*state % rho - screened_rates(k_p_o16__f17)*Y(jo16)*state % rho - &
      screened_rates(k_p_o16__he4_n13)*Y(jo16)*state % rho &
       )
    call set_jac_entry(state, jo16, jp, scratch)

    scratch = (&
      screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho + screened_rates(k_he4_n13__p_o16)* &
      Y(jn13)*state % rho - screened_rates(k_he4_o16__ne20)*Y(jo16)*state % rho - &
      screened_rates(k_he4_o16__p_f19)*Y(jo16)*state % rho &
       )
    call set_jac_entry(state, jo16, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jo16, jc12, scratch)

    scratch = (&
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jo16, jn13, scratch)

    scratch = (&
      screened_rates(k_p_n15__o16)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jo16, jn15, scratch)

    scratch = (&
      -screened_rates(k_he4_o16__ne20)*Y(jhe4)*state % rho - screened_rates(k_he4_o16__p_f19)* &
      Y(jhe4)*state % rho - screened_rates(k_o16__he4_c12) - screened_rates(k_o16__p_n15) &
      - screened_rates(k_p_o16__f17)*Y(jp)*state % rho - screened_rates(k_p_o16__he4_n13)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jo16, jo16, scratch)

    scratch = (&
      screened_rates(k_f17__p_o16) &
       )
    call set_jac_entry(state, jo16, jf17, scratch)

    scratch = (&
      screened_rates(k_p_f19__he4_o16)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jo16, jf19, scratch)

    scratch = (&
      screened_rates(k_ne17__p_o16__weak__wc12) &
       )
    call set_jac_entry(state, jo16, jne17, scratch)

    scratch = (&
      screened_rates(k_ne20__he4_o16) &
       )
    call set_jac_entry(state, jo16, jne20, scratch)

    scratch = (&
      screened_rates(k_na20__he4_o16__weak__wc12) &
       )
    call set_jac_entry(state, jo16, jna20, scratch)

    scratch = (&
      screened_rates(k_p_f20__he4_o17)*Y(jf20)*state % rho - screened_rates(k_p_o17__f18)* &
      Y(jo17)*state % rho - screened_rates(k_p_o17__he4_n14)*Y(jo17)*state % rho &
       )
    call set_jac_entry(state, jo17, jp, scratch)

    scratch = (&
      screened_rates(k_he4_n14__p_o17)*Y(jn14)*state % rho - screened_rates(k_he4_o17__ne21)* &
      Y(jo17)*state % rho - screened_rates(k_he4_o17__p_f20)*Y(jo17)*state % rho &
       )
    call set_jac_entry(state, jo17, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_n14__p_o17)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jo17, jn14, scratch)

    scratch = (&
      -screened_rates(k_he4_o17__ne21)*Y(jhe4)*state % rho - screened_rates(k_he4_o17__p_f20)* &
      Y(jhe4)*state % rho - screened_rates(k_p_o17__f18)*Y(jp)*state % rho - &
      screened_rates(k_p_o17__he4_n14)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jo17, jo17, scratch)

    scratch = (&
      screened_rates(k_f17__o17__weak__wc12) &
       )
    call set_jac_entry(state, jo17, jf17, scratch)

    scratch = (&
      screened_rates(k_f18__p_o17) &
       )
    call set_jac_entry(state, jo17, jf18, scratch)

    scratch = (&
      screened_rates(k_p_f20__he4_o17)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jo17, jf20, scratch)

    scratch = (&
      screened_rates(k_ne21__he4_o17) &
       )
    call set_jac_entry(state, jo17, jne21, scratch)

    scratch = (&
      -screened_rates(k_p_o18__f19)*Y(jo18)*state % rho - screened_rates(k_p_o18__he4_n15)* &
      Y(jo18)*state % rho &
       )
    call set_jac_entry(state, jo18, jp, scratch)

    scratch = (&
      screened_rates(k_he4_n15__p_o18)*Y(jn15)*state % rho - screened_rates(k_he4_o18__ne22)* &
      Y(jo18)*state % rho &
       )
    call set_jac_entry(state, jo18, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_n15__p_o18)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jo18, jn15, scratch)

    scratch = (&
      -screened_rates(k_he4_o18__ne22)*Y(jhe4)*state % rho - screened_rates(k_p_o18__f19)*Y(jp) &
      *state % rho - screened_rates(k_p_o18__he4_n15)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jo18, jo18, scratch)

    scratch = (&
      screened_rates(k_f18__o18__weak__wc12) &
       )
    call set_jac_entry(state, jo18, jf18, scratch)

    scratch = (&
      screened_rates(k_f19__p_o18) &
       )
    call set_jac_entry(state, jo18, jf19, scratch)

    scratch = (&
      screened_rates(k_p_ne19__he4_f16)*Y(jne19)*state % rho + screened_rates(k_p_o15__f16)* &
      Y(jo15)*state % rho &
       )
    call set_jac_entry(state, jf16, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_f16__na20)*Y(jf16)*state % rho - screened_rates(k_he4_f16__p_ne19)* &
      Y(jf16)*state % rho &
       )
    call set_jac_entry(state, jf16, jhe4, scratch)

    scratch = (&
      screened_rates(k_p_o15__f16)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jf16, jo15, scratch)

    scratch = (&
      -screened_rates(k_f16__p_o15) - screened_rates(k_he4_f16__na20)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_f16__p_ne19)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jf16, jf16, scratch)

    scratch = (&
      screened_rates(k_p_ne19__he4_f16)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jf16, jne19, scratch)

    scratch = (&
      screened_rates(k_na20__he4_f16) &
       )
    call set_jac_entry(state, jf16, jna20, scratch)

    scratch = (&
      -screened_rates(k_p_f17__he4_o14)*Y(jf17)*state % rho - screened_rates(k_p_f17__ne18)* &
      Y(jf17)*state % rho + screened_rates(k_p_ne20__he4_f17)*Y(jne20)*state % rho &
      + screened_rates(k_p_o16__f17)*Y(jo16)*state % rho &
       )
    call set_jac_entry(state, jf17, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_f17__na21)*Y(jf17)*state % rho - screened_rates(k_he4_f17__p_ne20)* &
      Y(jf17)*state % rho + screened_rates(k_he4_o14__p_f17)*Y(jo14)*state % rho &
       )
    call set_jac_entry(state, jf17, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_o14__p_f17)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jf17, jo14, scratch)

    scratch = (&
      screened_rates(k_p_o16__f17)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jf17, jo16, scratch)

    scratch = (&
      -screened_rates(k_f17__o17__weak__wc12) - screened_rates(k_f17__p_o16) - &
      screened_rates(k_he4_f17__na21)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_f17__p_ne20)*Y(jhe4)*state % rho - &
      screened_rates(k_p_f17__he4_o14)*Y(jp)*state % rho - screened_rates(k_p_f17__ne18)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jf17, jf17, scratch)

    scratch = (&
      screened_rates(k_ne17__f17__weak__wc17) &
       )
    call set_jac_entry(state, jf17, jne17, scratch)

    scratch = (&
      screened_rates(k_ne18__p_f17) &
       )
    call set_jac_entry(state, jf17, jne18, scratch)

    scratch = (&
      screened_rates(k_p_ne20__he4_f17)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jf17, jne20, scratch)

    scratch = (&
      screened_rates(k_na21__he4_f17) &
       )
    call set_jac_entry(state, jf17, jna21, scratch)

    scratch = (&
      -screened_rates(k_p_f18__he4_o15)*Y(jf18)*state % rho - screened_rates(k_p_f18__ne19)* &
      Y(jf18)*state % rho + screened_rates(k_p_ne21__he4_f18)*Y(jne21)*state % rho &
      + screened_rates(k_p_o17__f18)*Y(jo17)*state % rho &
       )
    call set_jac_entry(state, jf18, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_f18__na22)*Y(jf18)*state % rho - screened_rates(k_he4_f18__p_ne21)* &
      Y(jf18)*state % rho + screened_rates(k_he4_n14__f18)*Y(jn14)*state % rho + &
      screened_rates(k_he4_o15__p_f18)*Y(jo15)*state % rho &
       )
    call set_jac_entry(state, jf18, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_n14__f18)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jf18, jn14, scratch)

    scratch = (&
      screened_rates(k_he4_o15__p_f18)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jf18, jo15, scratch)

    scratch = (&
      screened_rates(k_p_o17__f18)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jf18, jo17, scratch)

    scratch = (&
      -screened_rates(k_f18__he4_n14) - screened_rates(k_f18__o18__weak__wc12) - &
      screened_rates(k_f18__p_o17) - screened_rates(k_he4_f18__na22)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_f18__p_ne21)*Y(jhe4)*state % rho - &
      screened_rates(k_p_f18__he4_o15)*Y(jp)*state % rho - screened_rates(k_p_f18__ne19)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jf18, jf18, scratch)

    scratch = (&
      screened_rates(k_ne18__f18__weak__wc12) &
       )
    call set_jac_entry(state, jf18, jne18, scratch)

    scratch = (&
      screened_rates(k_ne19__p_f18) &
       )
    call set_jac_entry(state, jf18, jne19, scratch)

    scratch = (&
      screened_rates(k_p_ne21__he4_f18)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jf18, jne21, scratch)

    scratch = (&
      screened_rates(k_na22__he4_f18) &
       )
    call set_jac_entry(state, jf18, jna22, scratch)

    scratch = (&
      -screened_rates(k_p_f19__he4_o16)*Y(jf19)*state % rho - screened_rates(k_p_f19__ne20)* &
      Y(jf19)*state % rho + screened_rates(k_p_ne22__he4_f19)*Y(jne22)*state % rho &
      + screened_rates(k_p_o18__f19)*Y(jo18)*state % rho &
       )
    call set_jac_entry(state, jf19, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_f19__na23)*Y(jf19)*state % rho - screened_rates(k_he4_f19__p_ne22)* &
      Y(jf19)*state % rho + screened_rates(k_he4_n15__f19)*Y(jn15)*state % rho + &
      screened_rates(k_he4_o16__p_f19)*Y(jo16)*state % rho &
       )
    call set_jac_entry(state, jf19, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_n15__f19)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jf19, jn15, scratch)

    scratch = (&
      screened_rates(k_he4_o16__p_f19)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jf19, jo16, scratch)

    scratch = (&
      screened_rates(k_p_o18__f19)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jf19, jo18, scratch)

    scratch = (&
      -screened_rates(k_f19__he4_n15) - screened_rates(k_f19__p_o18) - screened_rates(k_he4_f19__na23)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_f19__p_ne22)*Y(jhe4)*state % rho &
      - screened_rates(k_p_f19__he4_o16)*Y(jp)*state % rho - screened_rates(k_p_f19__ne20) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jf19, jf19, scratch)

    scratch = (&
      screened_rates(k_ne19__f19__weak__wc12) &
       )
    call set_jac_entry(state, jf19, jne19, scratch)

    scratch = (&
      screened_rates(k_ne20__p_f19) &
       )
    call set_jac_entry(state, jf19, jne20, scratch)

    scratch = (&
      screened_rates(k_p_ne22__he4_f19)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jf19, jne22, scratch)

    scratch = (&
      screened_rates(k_na23__he4_f19) &
       )
    call set_jac_entry(state, jf19, jna23, scratch)

    scratch = (&
      -screened_rates(k_p_f20__he4_o17)*Y(jf20)*state % rho - screened_rates(k_p_f20__ne21)* &
      Y(jf20)*state % rho &
       )
    call set_jac_entry(state, jf20, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_f20__na24)*Y(jf20)*state % rho + screened_rates(k_he4_o17__p_f20)* &
      Y(jo17)*state % rho &
       )
    call set_jac_entry(state, jf20, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_o17__p_f20)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jf20, jo17, scratch)

    scratch = (&
      -screened_rates(k_f20__ne20__weak__wc12) - screened_rates(k_he4_f20__na24)*Y(jhe4)* &
      state % rho - screened_rates(k_p_f20__he4_o17)*Y(jp)*state % rho - &
      screened_rates(k_p_f20__ne21)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jf20, jf20, scratch)

    scratch = (&
      screened_rates(k_ne21__p_f20) &
       )
    call set_jac_entry(state, jf20, jne21, scratch)

    scratch = (&
      screened_rates(k_p_na20__he4_ne17)*Y(jna20)*state % rho &
       )
    call set_jac_entry(state, jne17, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_ne17__mg21)*Y(jne17)*state % rho - screened_rates(k_he4_ne17__p_na20) &
      *Y(jne17)*state % rho &
       )
    call set_jac_entry(state, jne17, jhe4, scratch)

    scratch = (&
      -screened_rates(k_he4_ne17__mg21)*Y(jhe4)*state % rho - screened_rates(k_he4_ne17__p_na20)* &
      Y(jhe4)*state % rho - screened_rates(k_ne17__f17__weak__wc17) - &
      screened_rates(k_ne17__p_o16__weak__wc12) &
       )
    call set_jac_entry(state, jne17, jne17, scratch)

    scratch = (&
      screened_rates(k_p_na20__he4_ne17)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jne17, jna20, scratch)

    scratch = (&
      screened_rates(k_mg21__he4_ne17) &
       )
    call set_jac_entry(state, jne17, jmg21, scratch)

    scratch = (&
      screened_rates(k_p_f17__ne18)*Y(jf17)*state % rho + screened_rates(k_p_na21__he4_ne18)* &
      Y(jna21)*state % rho - screened_rates(k_p_ne18__na19)*Y(jne18)*state % rho &
       )
    call set_jac_entry(state, jne18, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_ne18__mg22)*Y(jne18)*state % rho - screened_rates(k_he4_ne18__p_na21) &
      *Y(jne18)*state % rho + screened_rates(k_he4_o14__ne18)*Y(jo14)*state % rho &
       )
    call set_jac_entry(state, jne18, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_o14__ne18)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jne18, jo14, scratch)

    scratch = (&
      screened_rates(k_p_f17__ne18)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jne18, jf17, scratch)

    scratch = (&
      -screened_rates(k_he4_ne18__mg22)*Y(jhe4)*state % rho - screened_rates(k_he4_ne18__p_na21)* &
      Y(jhe4)*state % rho - screened_rates(k_ne18__f18__weak__wc12) - &
      screened_rates(k_ne18__he4_o14) - screened_rates(k_ne18__p_f17) - &
      screened_rates(k_p_ne18__na19)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jne18, jne18, scratch)

    scratch = (&
      screened_rates(k_na19__p_ne18) &
       )
    call set_jac_entry(state, jne18, jna19, scratch)

    scratch = (&
      screened_rates(k_p_na21__he4_ne18)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jne18, jna21, scratch)

    scratch = (&
      screened_rates(k_mg22__he4_ne18) &
       )
    call set_jac_entry(state, jne18, jmg22, scratch)

    scratch = (&
      screened_rates(k_al22__he4_ne18__weak__wc12) &
       )
    call set_jac_entry(state, jne18, jal22, scratch)

    scratch = (&
      screened_rates(k_p_f18__ne19)*Y(jf18)*state % rho + screened_rates(k_p_na22__he4_ne19)* &
      Y(jna22)*state % rho - screened_rates(k_p_ne19__he4_f16)*Y(jne19)* &
      state % rho - screened_rates(k_p_ne19__na20)*Y(jne19)*state % rho &
       )
    call set_jac_entry(state, jne19, jp, scratch)

    scratch = (&
      screened_rates(k_he4_f16__p_ne19)*Y(jf16)*state % rho - screened_rates(k_he4_ne19__mg23)* &
      Y(jne19)*state % rho - screened_rates(k_he4_ne19__p_na22)*Y(jne19)* &
      state % rho + screened_rates(k_he4_o15__ne19)*Y(jo15)*state % rho &
       )
    call set_jac_entry(state, jne19, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_o15__ne19)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jne19, jo15, scratch)

    scratch = (&
      screened_rates(k_he4_f16__p_ne19)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jne19, jf16, scratch)

    scratch = (&
      screened_rates(k_p_f18__ne19)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jne19, jf18, scratch)

    scratch = (&
      -screened_rates(k_he4_ne19__mg23)*Y(jhe4)*state % rho - screened_rates(k_he4_ne19__p_na22)* &
      Y(jhe4)*state % rho - screened_rates(k_ne19__f19__weak__wc12) - &
      screened_rates(k_ne19__he4_o15) - screened_rates(k_ne19__p_f18) - &
      screened_rates(k_p_ne19__he4_f16)*Y(jp)*state % rho - screened_rates(k_p_ne19__na20) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jne19, jne19, scratch)

    scratch = (&
      screened_rates(k_na19__ne19__weak__bqa_pos_) &
       )
    call set_jac_entry(state, jne19, jna19, scratch)

    scratch = (&
      screened_rates(k_na20__p_ne19) &
       )
    call set_jac_entry(state, jne19, jna20, scratch)

    scratch = (&
      screened_rates(k_p_na22__he4_ne19)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jne19, jna22, scratch)

    scratch = (&
      screened_rates(k_mg23__he4_ne19) &
       )
    call set_jac_entry(state, jne19, jmg23, scratch)

    scratch = (&
      screened_rates(k_p_f19__ne20)*Y(jf19)*state % rho + screened_rates(k_p_na23__he4_ne20)* &
      Y(jna23)*state % rho - screened_rates(k_p_ne20__he4_f17)*Y(jne20)* &
      state % rho - screened_rates(k_p_ne20__na21)*Y(jne20)*state % rho &
       )
    call set_jac_entry(state, jne20, jp, scratch)

    scratch = (&
      screened_rates(k_he4_f17__p_ne20)*Y(jf17)*state % rho - screened_rates(k_he4_ne20__c12_c12) &
      *Y(jne20)*state % rho - screened_rates(k_he4_ne20__mg24)*Y(jne20)* &
      state % rho - screened_rates(k_he4_ne20__p_na23)*Y(jne20)*state % rho + &
      screened_rates(k_he4_o16__ne20)*Y(jo16)*state % rho &
       )
    call set_jac_entry(state, jne20, jhe4, scratch)

    scratch = (&
      1.0d0*screened_rates(k_c12_c12__he4_ne20)*Y(jc12)*state % rho &
       )
    call set_jac_entry(state, jne20, jc12, scratch)

    scratch = (&
      screened_rates(k_he4_o16__ne20)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jne20, jo16, scratch)

    scratch = (&
      screened_rates(k_he4_f17__p_ne20)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jne20, jf17, scratch)

    scratch = (&
      screened_rates(k_p_f19__ne20)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jne20, jf19, scratch)

    scratch = (&
      screened_rates(k_f20__ne20__weak__wc12) &
       )
    call set_jac_entry(state, jne20, jf20, scratch)

    scratch = (&
      -screened_rates(k_he4_ne20__c12_c12)*Y(jhe4)*state % rho - screened_rates(k_he4_ne20__mg24) &
      *Y(jhe4)*state % rho - screened_rates(k_he4_ne20__p_na23)*Y(jhe4)* &
      state % rho - screened_rates(k_ne20__he4_o16) - screened_rates(k_ne20__p_f19) - &
      screened_rates(k_p_ne20__he4_f17)*Y(jp)*state % rho - screened_rates(k_p_ne20__na21) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jne20, jne20, scratch)

    scratch = (&
      screened_rates(k_na20__ne20__weak__wc12) &
       )
    call set_jac_entry(state, jne20, jna20, scratch)

    scratch = (&
      screened_rates(k_na21__p_ne20) &
       )
    call set_jac_entry(state, jne20, jna21, scratch)

    scratch = (&
      screened_rates(k_p_na23__he4_ne20)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jne20, jna23, scratch)

    scratch = (&
      screened_rates(k_mg21__p_ne20__weak__wc12) &
       )
    call set_jac_entry(state, jne20, jmg21, scratch)

    scratch = (&
      screened_rates(k_mg24__he4_ne20) &
       )
    call set_jac_entry(state, jne20, jmg24, scratch)

    scratch = (&
      screened_rates(k_al22__p_p_ne20__weak__wc12) &
       )
    call set_jac_entry(state, jne20, jal22, scratch)

    scratch = (&
      screened_rates(k_al24__he4_ne20__weak__wc12) &
       )
    call set_jac_entry(state, jne20, jal24, scratch)

    scratch = (&
      screened_rates(k_p_f20__ne21)*Y(jf20)*state % rho + screened_rates(k_p_na24__he4_ne21)* &
      Y(jna24)*state % rho - screened_rates(k_p_ne21__he4_f18)*Y(jne21)* &
      state % rho - screened_rates(k_p_ne21__na22)*Y(jne21)*state % rho &
       )
    call set_jac_entry(state, jne21, jp, scratch)

    scratch = (&
      screened_rates(k_he4_f18__p_ne21)*Y(jf18)*state % rho - screened_rates(k_he4_ne21__mg25)* &
      Y(jne21)*state % rho - screened_rates(k_he4_ne21__p_na24)*Y(jne21)* &
      state % rho + screened_rates(k_he4_o17__ne21)*Y(jo17)*state % rho &
       )
    call set_jac_entry(state, jne21, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_o17__ne21)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jne21, jo17, scratch)

    scratch = (&
      screened_rates(k_he4_f18__p_ne21)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jne21, jf18, scratch)

    scratch = (&
      screened_rates(k_p_f20__ne21)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jne21, jf20, scratch)

    scratch = (&
      -screened_rates(k_he4_ne21__mg25)*Y(jhe4)*state % rho - screened_rates(k_he4_ne21__p_na24)* &
      Y(jhe4)*state % rho - screened_rates(k_ne21__he4_o17) - &
      screened_rates(k_ne21__p_f20) - screened_rates(k_p_ne21__he4_f18)*Y(jp)*state % rho &
      - screened_rates(k_p_ne21__na22)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jne21, jne21, scratch)

    scratch = (&
      screened_rates(k_na21__ne21__weak__wc12) &
       )
    call set_jac_entry(state, jne21, jna21, scratch)

    scratch = (&
      screened_rates(k_na22__p_ne21) &
       )
    call set_jac_entry(state, jne21, jna22, scratch)

    scratch = (&
      screened_rates(k_p_na24__he4_ne21)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jne21, jna24, scratch)

    scratch = (&
      screened_rates(k_mg25__he4_ne21) &
       )
    call set_jac_entry(state, jne21, jmg25, scratch)

    scratch = (&
      -screened_rates(k_p_ne22__he4_f19)*Y(jne22)*state % rho - screened_rates(k_p_ne22__na23)* &
      Y(jne22)*state % rho &
       )
    call set_jac_entry(state, jne22, jp, scratch)

    scratch = (&
      screened_rates(k_he4_f19__p_ne22)*Y(jf19)*state % rho - screened_rates(k_he4_ne22__mg26)* &
      Y(jne22)*state % rho + screened_rates(k_he4_o18__ne22)*Y(jo18)*state % rho &
       )
    call set_jac_entry(state, jne22, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_o18__ne22)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jne22, jo18, scratch)

    scratch = (&
      screened_rates(k_he4_f19__p_ne22)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jne22, jf19, scratch)

    scratch = (&
      -screened_rates(k_he4_ne22__mg26)*Y(jhe4)*state % rho - screened_rates(k_p_ne22__he4_f19)* &
      Y(jp)*state % rho - screened_rates(k_p_ne22__na23)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jne22, jne22, scratch)

    scratch = (&
      screened_rates(k_na22__ne22__weak__wc12) &
       )
    call set_jac_entry(state, jne22, jna22, scratch)

    scratch = (&
      screened_rates(k_na23__p_ne22) &
       )
    call set_jac_entry(state, jne22, jna23, scratch)

    scratch = (&
      screened_rates(k_mg26__he4_ne22) &
       )
    call set_jac_entry(state, jne22, jmg26, scratch)

    scratch = (&
      screened_rates(k_p_mg22__he4_na19)*Y(jmg22)*state % rho + screened_rates(k_p_ne18__na19)* &
      Y(jne18)*state % rho &
       )
    call set_jac_entry(state, jna19, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_na19__al23)*Y(jna19)*state % rho - screened_rates(k_he4_na19__p_mg22) &
      *Y(jna19)*state % rho &
       )
    call set_jac_entry(state, jna19, jhe4, scratch)

    scratch = (&
      screened_rates(k_p_ne18__na19)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jna19, jne18, scratch)

    scratch = (&
      -screened_rates(k_he4_na19__al23)*Y(jhe4)*state % rho - screened_rates(k_he4_na19__p_mg22)* &
      Y(jhe4)*state % rho - screened_rates(k_na19__ne19__weak__bqa_pos_) - &
      screened_rates(k_na19__p_ne18) &
       )
    call set_jac_entry(state, jna19, jna19, scratch)

    scratch = (&
      screened_rates(k_p_mg22__he4_na19)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jna19, jmg22, scratch)

    scratch = (&
      screened_rates(k_al23__he4_na19) &
       )
    call set_jac_entry(state, jna19, jal23, scratch)

    scratch = (&
      screened_rates(k_p_mg23__he4_na20)*Y(jmg23)*state % rho - &
      screened_rates(k_p_na20__he4_ne17)*Y(jna20)*state % rho - &
      screened_rates(k_p_na20__mg21)*Y(jna20)*state % rho + screened_rates(k_p_ne19__na20) &
      *Y(jne19)*state % rho &
       )
    call set_jac_entry(state, jna20, jp, scratch)

    scratch = (&
      screened_rates(k_he4_f16__na20)*Y(jf16)*state % rho - screened_rates(k_he4_na20__al24)* &
      Y(jna20)*state % rho - screened_rates(k_he4_na20__p_mg23)*Y(jna20)* &
      state % rho + screened_rates(k_he4_ne17__p_na20)*Y(jne17)*state % rho &
       )
    call set_jac_entry(state, jna20, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_f16__na20)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jna20, jf16, scratch)

    scratch = (&
      screened_rates(k_he4_ne17__p_na20)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jna20, jne17, scratch)

    scratch = (&
      screened_rates(k_p_ne19__na20)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jna20, jne19, scratch)

    scratch = (&
      -screened_rates(k_he4_na20__al24)*Y(jhe4)*state % rho - screened_rates(k_he4_na20__p_mg23)* &
      Y(jhe4)*state % rho - screened_rates(k_na20__he4_f16) - &
      screened_rates(k_na20__he4_o16__weak__wc12) - &
      screened_rates(k_na20__ne20__weak__wc12) - screened_rates(k_na20__p_ne19) - &
      screened_rates(k_p_na20__he4_ne17)*Y(jp)*state % rho - &
      screened_rates(k_p_na20__mg21)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jna20, jna20, scratch)

    scratch = (&
      screened_rates(k_mg21__p_na20) &
       )
    call set_jac_entry(state, jna20, jmg21, scratch)

    scratch = (&
      screened_rates(k_p_mg23__he4_na20)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jna20, jmg23, scratch)

    scratch = (&
      screened_rates(k_al24__he4_na20) &
       )
    call set_jac_entry(state, jna20, jal24, scratch)

    scratch = (&
      screened_rates(k_p_mg24__he4_na21)*Y(jmg24)*state % rho - &
      screened_rates(k_p_na21__he4_ne18)*Y(jna21)*state % rho - &
      screened_rates(k_p_na21__mg22)*Y(jna21)*state % rho + screened_rates(k_p_ne20__na21) &
      *Y(jne20)*state % rho &
       )
    call set_jac_entry(state, jna21, jp, scratch)

    scratch = (&
      screened_rates(k_he4_f17__na21)*Y(jf17)*state % rho - screened_rates(k_he4_na21__al25)* &
      Y(jna21)*state % rho - screened_rates(k_he4_na21__p_mg24)*Y(jna21)* &
      state % rho + screened_rates(k_he4_ne18__p_na21)*Y(jne18)*state % rho &
       )
    call set_jac_entry(state, jna21, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_f17__na21)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jna21, jf17, scratch)

    scratch = (&
      screened_rates(k_he4_ne18__p_na21)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jna21, jne18, scratch)

    scratch = (&
      screened_rates(k_p_ne20__na21)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jna21, jne20, scratch)

    scratch = (&
      -screened_rates(k_he4_na21__al25)*Y(jhe4)*state % rho - screened_rates(k_he4_na21__p_mg24)* &
      Y(jhe4)*state % rho - screened_rates(k_na21__he4_f17) - &
      screened_rates(k_na21__ne21__weak__wc12) - screened_rates(k_na21__p_ne20) - &
      screened_rates(k_p_na21__he4_ne18)*Y(jp)*state % rho - &
      screened_rates(k_p_na21__mg22)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jna21, jna21, scratch)

    scratch = (&
      screened_rates(k_mg21__na21__weak__wc12) &
       )
    call set_jac_entry(state, jna21, jmg21, scratch)

    scratch = (&
      screened_rates(k_mg22__p_na21) &
       )
    call set_jac_entry(state, jna21, jmg22, scratch)

    scratch = (&
      screened_rates(k_p_mg24__he4_na21)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jna21, jmg24, scratch)

    scratch = (&
      screened_rates(k_al22__p_na21__weak__wc17) &
       )
    call set_jac_entry(state, jna21, jal22, scratch)

    scratch = (&
      screened_rates(k_al25__he4_na21) &
       )
    call set_jac_entry(state, jna21, jal25, scratch)

    scratch = (&
      screened_rates(k_p_mg25__he4_na22)*Y(jmg25)*state % rho - &
      screened_rates(k_p_na22__he4_ne19)*Y(jna22)*state % rho - &
      screened_rates(k_p_na22__mg23)*Y(jna22)*state % rho + screened_rates(k_p_ne21__na22) &
      *Y(jne21)*state % rho &
       )
    call set_jac_entry(state, jna22, jp, scratch)

    scratch = (&
      screened_rates(k_he4_f18__na22)*Y(jf18)*state % rho - screened_rates(k_he4_na22__al26)* &
      Y(jna22)*state % rho - screened_rates(k_he4_na22__p_mg25)*Y(jna22)* &
      state % rho + screened_rates(k_he4_ne19__p_na22)*Y(jne19)*state % rho &
       )
    call set_jac_entry(state, jna22, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_f18__na22)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jna22, jf18, scratch)

    scratch = (&
      screened_rates(k_he4_ne19__p_na22)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jna22, jne19, scratch)

    scratch = (&
      screened_rates(k_p_ne21__na22)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jna22, jne21, scratch)

    scratch = (&
      -screened_rates(k_he4_na22__al26)*Y(jhe4)*state % rho - screened_rates(k_he4_na22__p_mg25)* &
      Y(jhe4)*state % rho - screened_rates(k_na22__he4_f18) - &
      screened_rates(k_na22__ne22__weak__wc12) - screened_rates(k_na22__p_ne21) - &
      screened_rates(k_p_na22__he4_ne19)*Y(jp)*state % rho - &
      screened_rates(k_p_na22__mg23)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jna22, jna22, scratch)

    scratch = (&
      screened_rates(k_mg22__na22__weak__wc12) &
       )
    call set_jac_entry(state, jna22, jmg22, scratch)

    scratch = (&
      screened_rates(k_mg23__p_na22) &
       )
    call set_jac_entry(state, jna22, jmg23, scratch)

    scratch = (&
      screened_rates(k_p_mg25__he4_na22)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jna22, jmg25, scratch)

    scratch = (&
      screened_rates(k_al23__p_na22__weak__wc12) &
       )
    call set_jac_entry(state, jna22, jal23, scratch)

    scratch = (&
      screened_rates(k_al26__he4_na22) &
       )
    call set_jac_entry(state, jna22, jal26, scratch)

    scratch = (&
      screened_rates(k_p_mg26__he4_na23)*Y(jmg26)*state % rho - &
      screened_rates(k_p_na23__he4_ne20)*Y(jna23)*state % rho - &
      screened_rates(k_p_na23__mg24)*Y(jna23)*state % rho + screened_rates(k_p_ne22__na23) &
      *Y(jne22)*state % rho &
       )
    call set_jac_entry(state, jna23, jp, scratch)

    scratch = (&
      screened_rates(k_he4_f19__na23)*Y(jf19)*state % rho - screened_rates(k_he4_na23__al27)* &
      Y(jna23)*state % rho - screened_rates(k_he4_na23__p_mg26)*Y(jna23)* &
      state % rho + screened_rates(k_he4_ne20__p_na23)*Y(jne20)*state % rho &
       )
    call set_jac_entry(state, jna23, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_f19__na23)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jna23, jf19, scratch)

    scratch = (&
      screened_rates(k_he4_ne20__p_na23)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jna23, jne20, scratch)

    scratch = (&
      screened_rates(k_p_ne22__na23)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jna23, jne22, scratch)

    scratch = (&
      -screened_rates(k_he4_na23__al27)*Y(jhe4)*state % rho - screened_rates(k_he4_na23__p_mg26)* &
      Y(jhe4)*state % rho - screened_rates(k_na23__he4_f19) - &
      screened_rates(k_na23__p_ne22) - screened_rates(k_p_na23__he4_ne20)*Y(jp)* &
      state % rho - screened_rates(k_p_na23__mg24)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jna23, jna23, scratch)

    scratch = (&
      screened_rates(k_mg23__na23__weak__wc12) &
       )
    call set_jac_entry(state, jna23, jmg23, scratch)

    scratch = (&
      screened_rates(k_mg24__p_na23) &
       )
    call set_jac_entry(state, jna23, jmg24, scratch)

    scratch = (&
      screened_rates(k_p_mg26__he4_na23)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jna23, jmg26, scratch)

    scratch = (&
      screened_rates(k_al24__p_na23__weak__wc12) &
       )
    call set_jac_entry(state, jna23, jal24, scratch)

    scratch = (&
      screened_rates(k_al27__he4_na23) &
       )
    call set_jac_entry(state, jna23, jal27, scratch)

    scratch = (&
      -screened_rates(k_p_na24__he4_ne21)*Y(jna24)*state % rho - screened_rates(k_p_na24__mg25)* &
      Y(jna24)*state % rho &
       )
    call set_jac_entry(state, jna24, jp, scratch)

    scratch = (&
      screened_rates(k_he4_f20__na24)*Y(jf20)*state % rho - screened_rates(k_he4_na24__al28)* &
      Y(jna24)*state % rho + screened_rates(k_he4_ne21__p_na24)*Y(jne21)* &
      state % rho &
       )
    call set_jac_entry(state, jna24, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_f20__na24)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jna24, jf20, scratch)

    scratch = (&
      screened_rates(k_he4_ne21__p_na24)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jna24, jne21, scratch)

    scratch = (&
      -screened_rates(k_he4_na24__al28)*Y(jhe4)*state % rho - screened_rates(k_p_na24__he4_ne21)* &
      Y(jp)*state % rho - screened_rates(k_p_na24__mg25)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jna24, jna24, scratch)

    scratch = (&
      screened_rates(k_mg25__p_na24) &
       )
    call set_jac_entry(state, jna24, jmg25, scratch)

    scratch = (&
      screened_rates(k_al28__he4_na24) &
       )
    call set_jac_entry(state, jna24, jal28, scratch)

    scratch = (&
      screened_rates(k_p_al24__he4_mg21)*Y(jal24)*state % rho - screened_rates(k_p_mg21__al22)* &
      Y(jmg21)*state % rho + screened_rates(k_p_na20__mg21)*Y(jna20)*state % rho &
       )
    call set_jac_entry(state, jmg21, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_mg21__p_al24)*Y(jmg21)*state % rho - screened_rates(k_he4_mg21__si25) &
      *Y(jmg21)*state % rho + screened_rates(k_he4_ne17__mg21)*Y(jne17)* &
      state % rho &
       )
    call set_jac_entry(state, jmg21, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ne17__mg21)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmg21, jne17, scratch)

    scratch = (&
      screened_rates(k_p_na20__mg21)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmg21, jna20, scratch)

    scratch = (&
      -screened_rates(k_he4_mg21__p_al24)*Y(jhe4)*state % rho - screened_rates(k_he4_mg21__si25)* &
      Y(jhe4)*state % rho - screened_rates(k_mg21__he4_ne17) - &
      screened_rates(k_mg21__na21__weak__wc12) - screened_rates(k_mg21__p_na20) - &
      screened_rates(k_mg21__p_ne20__weak__wc12) - screened_rates(k_p_mg21__al22)*Y(jp) &
      *state % rho &
       )
    call set_jac_entry(state, jmg21, jmg21, scratch)

    scratch = (&
      screened_rates(k_al22__p_mg21) &
       )
    call set_jac_entry(state, jmg21, jal22, scratch)

    scratch = (&
      screened_rates(k_p_al24__he4_mg21)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmg21, jal24, scratch)

    scratch = (&
      screened_rates(k_si25__he4_mg21) &
       )
    call set_jac_entry(state, jmg21, jsi25, scratch)

    scratch = (&
      screened_rates(k_p_al25__he4_mg22)*Y(jal25)*state % rho - screened_rates(k_p_mg22__al23)* &
      Y(jmg22)*state % rho - screened_rates(k_p_mg22__he4_na19)*Y(jmg22)* &
      state % rho + screened_rates(k_p_na21__mg22)*Y(jna21)*state % rho &
       )
    call set_jac_entry(state, jmg22, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_mg22__p_al25)*Y(jmg22)*state % rho - screened_rates(k_he4_mg22__si26) &
      *Y(jmg22)*state % rho + screened_rates(k_he4_na19__p_mg22)*Y(jna19)* &
      state % rho + screened_rates(k_he4_ne18__mg22)*Y(jne18)*state % rho &
       )
    call set_jac_entry(state, jmg22, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ne18__mg22)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmg22, jne18, scratch)

    scratch = (&
      screened_rates(k_he4_na19__p_mg22)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmg22, jna19, scratch)

    scratch = (&
      screened_rates(k_p_na21__mg22)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmg22, jna21, scratch)

    scratch = (&
      -screened_rates(k_he4_mg22__p_al25)*Y(jhe4)*state % rho - screened_rates(k_he4_mg22__si26)* &
      Y(jhe4)*state % rho - screened_rates(k_mg22__he4_ne18) - &
      screened_rates(k_mg22__na22__weak__wc12) - screened_rates(k_mg22__p_na21) - &
      screened_rates(k_p_mg22__al23)*Y(jp)*state % rho - &
      screened_rates(k_p_mg22__he4_na19)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmg22, jmg22, scratch)

    scratch = (&
      screened_rates(k_al22__mg22__weak__wc12) &
       )
    call set_jac_entry(state, jmg22, jal22, scratch)

    scratch = (&
      screened_rates(k_al23__p_mg22) &
       )
    call set_jac_entry(state, jmg22, jal23, scratch)

    scratch = (&
      screened_rates(k_p_al25__he4_mg22)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmg22, jal25, scratch)

    scratch = (&
      screened_rates(k_si26__he4_mg22) &
       )
    call set_jac_entry(state, jmg22, jsi26, scratch)

    scratch = (&
      screened_rates(k_p_al26__he4_mg23)*Y(jal26)*state % rho - screened_rates(k_p_mg23__al24)* &
      Y(jmg23)*state % rho - screened_rates(k_p_mg23__he4_na20)*Y(jmg23)* &
      state % rho + screened_rates(k_p_na22__mg23)*Y(jna22)*state % rho &
       )
    call set_jac_entry(state, jmg23, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_mg23__p_al26)*Y(jmg23)*state % rho - screened_rates(k_he4_mg23__si27) &
      *Y(jmg23)*state % rho + screened_rates(k_he4_na20__p_mg23)*Y(jna20)* &
      state % rho + screened_rates(k_he4_ne19__mg23)*Y(jne19)*state % rho &
       )
    call set_jac_entry(state, jmg23, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ne19__mg23)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmg23, jne19, scratch)

    scratch = (&
      screened_rates(k_he4_na20__p_mg23)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmg23, jna20, scratch)

    scratch = (&
      screened_rates(k_p_na22__mg23)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmg23, jna22, scratch)

    scratch = (&
      -screened_rates(k_he4_mg23__p_al26)*Y(jhe4)*state % rho - screened_rates(k_he4_mg23__si27)* &
      Y(jhe4)*state % rho - screened_rates(k_mg23__he4_ne19) - &
      screened_rates(k_mg23__na23__weak__wc12) - screened_rates(k_mg23__p_na22) - &
      screened_rates(k_p_mg23__al24)*Y(jp)*state % rho - &
      screened_rates(k_p_mg23__he4_na20)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmg23, jmg23, scratch)

    scratch = (&
      screened_rates(k_al23__mg23__weak__wc12) &
       )
    call set_jac_entry(state, jmg23, jal23, scratch)

    scratch = (&
      screened_rates(k_al24__p_mg23) &
       )
    call set_jac_entry(state, jmg23, jal24, scratch)

    scratch = (&
      screened_rates(k_p_al26__he4_mg23)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmg23, jal26, scratch)

    scratch = (&
      screened_rates(k_si24__p_mg23__weak__wc12) &
       )
    call set_jac_entry(state, jmg23, jsi24, scratch)

    scratch = (&
      screened_rates(k_si27__he4_mg23) &
       )
    call set_jac_entry(state, jmg23, jsi27, scratch)

    scratch = (&
      screened_rates(k_p_al27__he4_mg24)*Y(jal27)*state % rho - screened_rates(k_p_mg24__al25)* &
      Y(jmg24)*state % rho - screened_rates(k_p_mg24__he4_na21)*Y(jmg24)* &
      state % rho + screened_rates(k_p_na23__mg24)*Y(jna23)*state % rho &
       )
    call set_jac_entry(state, jmg24, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_mg24__p_al27)*Y(jmg24)*state % rho - screened_rates(k_he4_mg24__si28) &
      *Y(jmg24)*state % rho + screened_rates(k_he4_na21__p_mg24)*Y(jna21)* &
      state % rho + screened_rates(k_he4_ne20__mg24)*Y(jne20)*state % rho &
       )
    call set_jac_entry(state, jmg24, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ne20__mg24)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmg24, jne20, scratch)

    scratch = (&
      screened_rates(k_he4_na21__p_mg24)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmg24, jna21, scratch)

    scratch = (&
      screened_rates(k_p_na23__mg24)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmg24, jna23, scratch)

    scratch = (&
      -screened_rates(k_he4_mg24__p_al27)*Y(jhe4)*state % rho - screened_rates(k_he4_mg24__si28)* &
      Y(jhe4)*state % rho - screened_rates(k_mg24__he4_ne20) - &
      screened_rates(k_mg24__p_na23) - screened_rates(k_p_mg24__al25)*Y(jp)*state % rho - &
      screened_rates(k_p_mg24__he4_na21)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmg24, jmg24, scratch)

    scratch = (&
      screened_rates(k_al24__mg24__weak__wc12) &
       )
    call set_jac_entry(state, jmg24, jal24, scratch)

    scratch = (&
      screened_rates(k_al25__p_mg24) &
       )
    call set_jac_entry(state, jmg24, jal25, scratch)

    scratch = (&
      screened_rates(k_p_al27__he4_mg24)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmg24, jal27, scratch)

    scratch = (&
      screened_rates(k_si25__p_mg24__weak__wc12) &
       )
    call set_jac_entry(state, jmg24, jsi25, scratch)

    scratch = (&
      screened_rates(k_si28__he4_mg24) &
       )
    call set_jac_entry(state, jmg24, jsi28, scratch)

    scratch = (&
      screened_rates(k_p28__he4_mg24__weak__wc12) &
       )
    call set_jac_entry(state, jmg24, jp28, scratch)

    scratch = (&
      screened_rates(k_p_al28__he4_mg25)*Y(jal28)*state % rho - screened_rates(k_p_mg25__al26)* &
      Y(jmg25)*state % rho - screened_rates(k_p_mg25__he4_na22)*Y(jmg25)* &
      state % rho + screened_rates(k_p_na24__mg25)*Y(jna24)*state % rho &
       )
    call set_jac_entry(state, jmg25, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_mg25__p_al28)*Y(jmg25)*state % rho - screened_rates(k_he4_mg25__si29) &
      *Y(jmg25)*state % rho + screened_rates(k_he4_na22__p_mg25)*Y(jna22)* &
      state % rho + screened_rates(k_he4_ne21__mg25)*Y(jne21)*state % rho &
       )
    call set_jac_entry(state, jmg25, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ne21__mg25)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmg25, jne21, scratch)

    scratch = (&
      screened_rates(k_he4_na22__p_mg25)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmg25, jna22, scratch)

    scratch = (&
      screened_rates(k_p_na24__mg25)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmg25, jna24, scratch)

    scratch = (&
      -screened_rates(k_he4_mg25__p_al28)*Y(jhe4)*state % rho - screened_rates(k_he4_mg25__si29)* &
      Y(jhe4)*state % rho - screened_rates(k_mg25__he4_ne21) - &
      screened_rates(k_mg25__p_na24) - screened_rates(k_p_mg25__al26)*Y(jp)*state % rho - &
      screened_rates(k_p_mg25__he4_na22)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmg25, jmg25, scratch)

    scratch = (&
      screened_rates(k_al25__mg25__weak__wc12) &
       )
    call set_jac_entry(state, jmg25, jal25, scratch)

    scratch = (&
      screened_rates(k_al26__p_mg25) &
       )
    call set_jac_entry(state, jmg25, jal26, scratch)

    scratch = (&
      screened_rates(k_p_al28__he4_mg25)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmg25, jal28, scratch)

    scratch = (&
      screened_rates(k_si29__he4_mg25) &
       )
    call set_jac_entry(state, jmg25, jsi29, scratch)

    scratch = (&
      -screened_rates(k_p_mg26__al27)*Y(jmg26)*state % rho - screened_rates(k_p_mg26__he4_na23)* &
      Y(jmg26)*state % rho &
       )
    call set_jac_entry(state, jmg26, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_mg26__si30)*Y(jmg26)*state % rho + screened_rates(k_he4_na23__p_mg26) &
      *Y(jna23)*state % rho + screened_rates(k_he4_ne22__mg26)*Y(jne22)* &
      state % rho &
       )
    call set_jac_entry(state, jmg26, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ne22__mg26)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmg26, jne22, scratch)

    scratch = (&
      screened_rates(k_he4_na23__p_mg26)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmg26, jna23, scratch)

    scratch = (&
      -screened_rates(k_he4_mg26__si30)*Y(jhe4)*state % rho - screened_rates(k_mg26__he4_ne22) - &
      screened_rates(k_p_mg26__al27)*Y(jp)*state % rho - &
      screened_rates(k_p_mg26__he4_na23)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmg26, jmg26, scratch)

    scratch = (&
      screened_rates(k_al26__mg26__weak__wc12) &
       )
    call set_jac_entry(state, jmg26, jal26, scratch)

    scratch = (&
      screened_rates(k_al27__p_mg26) &
       )
    call set_jac_entry(state, jmg26, jal27, scratch)

    scratch = (&
      screened_rates(k_si30__he4_mg26) &
       )
    call set_jac_entry(state, jmg26, jsi30, scratch)

    scratch = (&
      screened_rates(k_p_mg21__al22)*Y(jmg21)*state % rho + screened_rates(k_p_si25__he4_al22)* &
      Y(jsi25)*state % rho &
       )
    call set_jac_entry(state, jal22, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_al22__p26)*Y(jal22)*state % rho - screened_rates(k_he4_al22__p_si25)* &
      Y(jal22)*state % rho &
       )
    call set_jac_entry(state, jal22, jhe4, scratch)

    scratch = (&
      screened_rates(k_p_mg21__al22)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jal22, jmg21, scratch)

    scratch = (&
      -screened_rates(k_al22__he4_ne18__weak__wc12) - screened_rates(k_al22__mg22__weak__wc12) - &
      screened_rates(k_al22__p_mg21) - screened_rates(k_al22__p_na21__weak__wc17) - &
      screened_rates(k_al22__p_p_ne20__weak__wc12) - screened_rates(k_he4_al22__p26)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_al22__p_si25)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jal22, jal22, scratch)

    scratch = (&
      screened_rates(k_p_si25__he4_al22)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jal22, jsi25, scratch)

    scratch = (&
      screened_rates(k_p26__he4_al22) &
       )
    call set_jac_entry(state, jal22, jp26, scratch)

    scratch = (&
      -screened_rates(k_p_al23__si24)*Y(jal23)*state % rho + screened_rates(k_p_mg22__al23)* &
      Y(jmg22)*state % rho + screened_rates(k_p_si26__he4_al23)*Y(jsi26)* &
      state % rho &
       )
    call set_jac_entry(state, jal23, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_al23__p27)*Y(jal23)*state % rho - screened_rates(k_he4_al23__p_si26)* &
      Y(jal23)*state % rho + screened_rates(k_he4_na19__al23)*Y(jna19)*state % rho &
       )
    call set_jac_entry(state, jal23, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_na19__al23)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jal23, jna19, scratch)

    scratch = (&
      screened_rates(k_p_mg22__al23)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jal23, jmg22, scratch)

    scratch = (&
      -screened_rates(k_al23__he4_na19) - screened_rates(k_al23__mg23__weak__wc12) - &
      screened_rates(k_al23__p_mg22) - screened_rates(k_al23__p_na22__weak__wc12) - &
      screened_rates(k_he4_al23__p27)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al23__p_si26)*Y(jhe4)*state % rho - &
      screened_rates(k_p_al23__si24)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jal23, jal23, scratch)

    scratch = (&
      screened_rates(k_si24__p_al23) &
       )
    call set_jac_entry(state, jal23, jsi24, scratch)

    scratch = (&
      screened_rates(k_p_si26__he4_al23)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jal23, jsi26, scratch)

    scratch = (&
      screened_rates(k_p27__he4_al23) &
       )
    call set_jac_entry(state, jal23, jp27, scratch)

    scratch = (&
      -screened_rates(k_p_al24__he4_mg21)*Y(jal24)*state % rho - screened_rates(k_p_al24__si25)* &
      Y(jal24)*state % rho + screened_rates(k_p_mg23__al24)*Y(jmg23)*state % rho + &
      screened_rates(k_p_si27__he4_al24)*Y(jsi27)*state % rho &
       )
    call set_jac_entry(state, jal24, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_al24__p28)*Y(jal24)*state % rho - screened_rates(k_he4_al24__p_si27)* &
      Y(jal24)*state % rho + screened_rates(k_he4_mg21__p_al24)*Y(jmg21)* &
      state % rho + screened_rates(k_he4_na20__al24)*Y(jna20)*state % rho &
       )
    call set_jac_entry(state, jal24, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_na20__al24)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jal24, jna20, scratch)

    scratch = (&
      screened_rates(k_he4_mg21__p_al24)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jal24, jmg21, scratch)

    scratch = (&
      screened_rates(k_p_mg23__al24)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jal24, jmg23, scratch)

    scratch = (&
      -screened_rates(k_al24__he4_na20) - screened_rates(k_al24__he4_ne20__weak__wc12) - &
      screened_rates(k_al24__mg24__weak__wc12) - screened_rates(k_al24__p_mg23) - &
      screened_rates(k_al24__p_na23__weak__wc12) - screened_rates(k_he4_al24__p28)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_al24__p_si27)*Y(jhe4)*state % rho &
      - screened_rates(k_p_al24__he4_mg21)*Y(jp)*state % rho - &
      screened_rates(k_p_al24__si25)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jal24, jal24, scratch)

    scratch = (&
      screened_rates(k_si24__al24__weak__wc12) &
       )
    call set_jac_entry(state, jal24, jsi24, scratch)

    scratch = (&
      screened_rates(k_si25__p_al24) &
       )
    call set_jac_entry(state, jal24, jsi25, scratch)

    scratch = (&
      screened_rates(k_p_si27__he4_al24)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jal24, jsi27, scratch)

    scratch = (&
      screened_rates(k_p28__he4_al24) &
       )
    call set_jac_entry(state, jal24, jp28, scratch)

    scratch = (&
      -screened_rates(k_p_al25__he4_mg22)*Y(jal25)*state % rho - screened_rates(k_p_al25__si26)* &
      Y(jal25)*state % rho + screened_rates(k_p_mg24__al25)*Y(jmg24)*state % rho + &
      screened_rates(k_p_si28__he4_al25)*Y(jsi28)*state % rho &
       )
    call set_jac_entry(state, jal25, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_al25__p29)*Y(jal25)*state % rho - screened_rates(k_he4_al25__p_si28)* &
      Y(jal25)*state % rho + screened_rates(k_he4_mg22__p_al25)*Y(jmg22)* &
      state % rho + screened_rates(k_he4_na21__al25)*Y(jna21)*state % rho &
       )
    call set_jac_entry(state, jal25, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_na21__al25)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jal25, jna21, scratch)

    scratch = (&
      screened_rates(k_he4_mg22__p_al25)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jal25, jmg22, scratch)

    scratch = (&
      screened_rates(k_p_mg24__al25)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jal25, jmg24, scratch)

    scratch = (&
      -screened_rates(k_al25__he4_na21) - screened_rates(k_al25__mg25__weak__wc12) - &
      screened_rates(k_al25__p_mg24) - screened_rates(k_he4_al25__p29)*Y(jhe4)*state % rho &
      - screened_rates(k_he4_al25__p_si28)*Y(jhe4)*state % rho - &
      screened_rates(k_p_al25__he4_mg22)*Y(jp)*state % rho - &
      screened_rates(k_p_al25__si26)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jal25, jal25, scratch)

    scratch = (&
      screened_rates(k_si25__al25__weak__wc12) &
       )
    call set_jac_entry(state, jal25, jsi25, scratch)

    scratch = (&
      screened_rates(k_si26__p_al25) &
       )
    call set_jac_entry(state, jal25, jsi26, scratch)

    scratch = (&
      screened_rates(k_p_si28__he4_al25)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jal25, jsi28, scratch)

    scratch = (&
      screened_rates(k_p29__he4_al25) &
       )
    call set_jac_entry(state, jal25, jp29, scratch)

    scratch = (&
      -screened_rates(k_p_al26__he4_mg23)*Y(jal26)*state % rho - screened_rates(k_p_al26__si27)* &
      Y(jal26)*state % rho + screened_rates(k_p_mg25__al26)*Y(jmg25)*state % rho + &
      screened_rates(k_p_si29__he4_al26)*Y(jsi29)*state % rho &
       )
    call set_jac_entry(state, jal26, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_al26__p30)*Y(jal26)*state % rho - screened_rates(k_he4_al26__p_si29)* &
      Y(jal26)*state % rho + screened_rates(k_he4_mg23__p_al26)*Y(jmg23)* &
      state % rho + screened_rates(k_he4_na22__al26)*Y(jna22)*state % rho &
       )
    call set_jac_entry(state, jal26, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_na22__al26)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jal26, jna22, scratch)

    scratch = (&
      screened_rates(k_he4_mg23__p_al26)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jal26, jmg23, scratch)

    scratch = (&
      screened_rates(k_p_mg25__al26)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jal26, jmg25, scratch)

    scratch = (&
      -screened_rates(k_al26__he4_na22) - screened_rates(k_al26__mg26__weak__wc12) - &
      screened_rates(k_al26__p_mg25) - screened_rates(k_he4_al26__p30)*Y(jhe4)*state % rho &
      - screened_rates(k_he4_al26__p_si29)*Y(jhe4)*state % rho - &
      screened_rates(k_p_al26__he4_mg23)*Y(jp)*state % rho - &
      screened_rates(k_p_al26__si27)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jal26, jal26, scratch)

    scratch = (&
      screened_rates(k_si26__al26__weak__wc12) &
       )
    call set_jac_entry(state, jal26, jsi26, scratch)

    scratch = (&
      screened_rates(k_si27__p_al26) &
       )
    call set_jac_entry(state, jal26, jsi27, scratch)

    scratch = (&
      screened_rates(k_p_si29__he4_al26)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jal26, jsi29, scratch)

    scratch = (&
      screened_rates(k_p27__p_al26__weak__wc12) &
       )
    call set_jac_entry(state, jal26, jp27, scratch)

    scratch = (&
      screened_rates(k_p30__he4_al26) &
       )
    call set_jac_entry(state, jal26, jp30, scratch)

    scratch = (&
      -screened_rates(k_p_al27__he4_mg24)*Y(jal27)*state % rho - screened_rates(k_p_al27__si28)* &
      Y(jal27)*state % rho + screened_rates(k_p_mg26__al27)*Y(jmg26)*state % rho + &
      screened_rates(k_p_si30__he4_al27)*Y(jsi30)*state % rho &
       )
    call set_jac_entry(state, jal27, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_al27__p31)*Y(jal27)*state % rho - screened_rates(k_he4_al27__p_si30)* &
      Y(jal27)*state % rho + screened_rates(k_he4_mg24__p_al27)*Y(jmg24)* &
      state % rho + screened_rates(k_he4_na23__al27)*Y(jna23)*state % rho &
       )
    call set_jac_entry(state, jal27, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_na23__al27)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jal27, jna23, scratch)

    scratch = (&
      screened_rates(k_he4_mg24__p_al27)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jal27, jmg24, scratch)

    scratch = (&
      screened_rates(k_p_mg26__al27)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jal27, jmg26, scratch)

    scratch = (&
      -screened_rates(k_al27__he4_na23) - screened_rates(k_al27__p_mg26) - &
      screened_rates(k_he4_al27__p31)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_al27__p_si30)*Y(jhe4)*state % rho - &
      screened_rates(k_p_al27__he4_mg24)*Y(jp)*state % rho - &
      screened_rates(k_p_al27__si28)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jal27, jal27, scratch)

    scratch = (&
      screened_rates(k_si27__al27__weak__wc12) &
       )
    call set_jac_entry(state, jal27, jsi27, scratch)

    scratch = (&
      screened_rates(k_si28__p_al27) &
       )
    call set_jac_entry(state, jal27, jsi28, scratch)

    scratch = (&
      screened_rates(k_p_si30__he4_al27)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jal27, jsi30, scratch)

    scratch = (&
      screened_rates(k_p28__p_al27__weak__wc12) &
       )
    call set_jac_entry(state, jal27, jp28, scratch)

    scratch = (&
      screened_rates(k_p31__he4_al27) &
       )
    call set_jac_entry(state, jal27, jp31, scratch)

    scratch = (&
      -screened_rates(k_p_al28__he4_mg25)*Y(jal28)*state % rho - screened_rates(k_p_al28__si29)* &
      Y(jal28)*state % rho &
       )
    call set_jac_entry(state, jal28, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_al28__p32)*Y(jal28)*state % rho + screened_rates(k_he4_mg25__p_al28)* &
      Y(jmg25)*state % rho + screened_rates(k_he4_na24__al28)*Y(jna24)*state % rho &
       )
    call set_jac_entry(state, jal28, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_na24__al28)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jal28, jna24, scratch)

    scratch = (&
      screened_rates(k_he4_mg25__p_al28)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jal28, jmg25, scratch)

    scratch = (&
      -screened_rates(k_al28__he4_na24) - screened_rates(k_he4_al28__p32)*Y(jhe4)*state % rho - &
      screened_rates(k_p_al28__he4_mg25)*Y(jp)*state % rho - &
      screened_rates(k_p_al28__si29)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jal28, jal28, scratch)

    scratch = (&
      screened_rates(k_si29__p_al28) &
       )
    call set_jac_entry(state, jal28, jsi29, scratch)

    scratch = (&
      screened_rates(k_p32__he4_al28) &
       )
    call set_jac_entry(state, jal28, jp32, scratch)

    scratch = (&
      screened_rates(k_p_al23__si24)*Y(jal23)*state % rho + screened_rates(k_p_p27__he4_si24)* &
      Y(jp27)*state % rho &
       )
    call set_jac_entry(state, jsi24, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_si24__p_p27)*Y(jsi24)*state % rho - screened_rates(k_he4_si24__s28)* &
      Y(jsi24)*state % rho &
       )
    call set_jac_entry(state, jsi24, jhe4, scratch)

    scratch = (&
      screened_rates(k_p_al23__si24)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsi24, jal23, scratch)

    scratch = (&
      -screened_rates(k_he4_si24__p_p27)*Y(jhe4)*state % rho - screened_rates(k_he4_si24__s28)* &
      Y(jhe4)*state % rho - screened_rates(k_si24__al24__weak__wc12) - &
      screened_rates(k_si24__p_al23) - screened_rates(k_si24__p_mg23__weak__wc12) &
       )
    call set_jac_entry(state, jsi24, jsi24, scratch)

    scratch = (&
      screened_rates(k_p_p27__he4_si24)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsi24, jp27, scratch)

    scratch = (&
      screened_rates(k_s28__he4_si24) &
       )
    call set_jac_entry(state, jsi24, js28, scratch)

    scratch = (&
      screened_rates(k_p_al24__si25)*Y(jal24)*state % rho + screened_rates(k_p_p28__he4_si25)* &
      Y(jp28)*state % rho - screened_rates(k_p_si25__he4_al22)*Y(jsi25)* &
      state % rho - screened_rates(k_p_si25__p26)*Y(jsi25)*state % rho &
       )
    call set_jac_entry(state, jsi25, jp, scratch)

    scratch = (&
      screened_rates(k_he4_al22__p_si25)*Y(jal22)*state % rho + screened_rates(k_he4_mg21__si25)* &
      Y(jmg21)*state % rho - screened_rates(k_he4_si25__p_p28)*Y(jsi25)* &
      state % rho - screened_rates(k_he4_si25__s29)*Y(jsi25)*state % rho &
       )
    call set_jac_entry(state, jsi25, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_mg21__si25)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsi25, jmg21, scratch)

    scratch = (&
      screened_rates(k_he4_al22__p_si25)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsi25, jal22, scratch)

    scratch = (&
      screened_rates(k_p_al24__si25)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsi25, jal24, scratch)

    scratch = (&
      -screened_rates(k_he4_si25__p_p28)*Y(jhe4)*state % rho - screened_rates(k_he4_si25__s29)* &
      Y(jhe4)*state % rho - screened_rates(k_p_si25__he4_al22)*Y(jp)*state % rho - &
      screened_rates(k_p_si25__p26)*Y(jp)*state % rho - &
      screened_rates(k_si25__al25__weak__wc12) - screened_rates(k_si25__he4_mg21) - &
      screened_rates(k_si25__p_al24) - screened_rates(k_si25__p_mg24__weak__wc12) &
       )
    call set_jac_entry(state, jsi25, jsi25, scratch)

    scratch = (&
      screened_rates(k_p26__p_si25) &
       )
    call set_jac_entry(state, jsi25, jp26, scratch)

    scratch = (&
      screened_rates(k_p_p28__he4_si25)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsi25, jp28, scratch)

    scratch = (&
      screened_rates(k_s29__he4_si25) &
       )
    call set_jac_entry(state, jsi25, js29, scratch)

    scratch = (&
      screened_rates(k_p_al25__si26)*Y(jal25)*state % rho + screened_rates(k_p_p29__he4_si26)* &
      Y(jp29)*state % rho - screened_rates(k_p_si26__he4_al23)*Y(jsi26)* &
      state % rho - screened_rates(k_p_si26__p27)*Y(jsi26)*state % rho &
       )
    call set_jac_entry(state, jsi26, jp, scratch)

    scratch = (&
      screened_rates(k_he4_al23__p_si26)*Y(jal23)*state % rho + screened_rates(k_he4_mg22__si26)* &
      Y(jmg22)*state % rho - screened_rates(k_he4_si26__p_p29)*Y(jsi26)* &
      state % rho - screened_rates(k_he4_si26__s30)*Y(jsi26)*state % rho &
       )
    call set_jac_entry(state, jsi26, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_mg22__si26)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsi26, jmg22, scratch)

    scratch = (&
      screened_rates(k_he4_al23__p_si26)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsi26, jal23, scratch)

    scratch = (&
      screened_rates(k_p_al25__si26)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsi26, jal25, scratch)

    scratch = (&
      -screened_rates(k_he4_si26__p_p29)*Y(jhe4)*state % rho - screened_rates(k_he4_si26__s30)* &
      Y(jhe4)*state % rho - screened_rates(k_p_si26__he4_al23)*Y(jp)*state % rho - &
      screened_rates(k_p_si26__p27)*Y(jp)*state % rho - &
      screened_rates(k_si26__al26__weak__wc12) - screened_rates(k_si26__he4_mg22) - &
      screened_rates(k_si26__p_al25) &
       )
    call set_jac_entry(state, jsi26, jsi26, scratch)

    scratch = (&
      screened_rates(k_p26__si26__weak__wc12) &
       )
    call set_jac_entry(state, jsi26, jp26, scratch)

    scratch = (&
      screened_rates(k_p27__p_si26) &
       )
    call set_jac_entry(state, jsi26, jp27, scratch)

    scratch = (&
      screened_rates(k_p_p29__he4_si26)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsi26, jp29, scratch)

    scratch = (&
      screened_rates(k_s30__he4_si26) &
       )
    call set_jac_entry(state, jsi26, js30, scratch)

    scratch = (&
      screened_rates(k_p_al26__si27)*Y(jal26)*state % rho + screened_rates(k_p_p30__he4_si27)* &
      Y(jp30)*state % rho - screened_rates(k_p_si27__he4_al24)*Y(jsi27)* &
      state % rho - screened_rates(k_p_si27__p28)*Y(jsi27)*state % rho &
       )
    call set_jac_entry(state, jsi27, jp, scratch)

    scratch = (&
      screened_rates(k_he4_al24__p_si27)*Y(jal24)*state % rho + screened_rates(k_he4_mg23__si27)* &
      Y(jmg23)*state % rho - screened_rates(k_he4_si27__p_p30)*Y(jsi27)* &
      state % rho - screened_rates(k_he4_si27__s31)*Y(jsi27)*state % rho &
       )
    call set_jac_entry(state, jsi27, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_mg23__si27)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsi27, jmg23, scratch)

    scratch = (&
      screened_rates(k_he4_al24__p_si27)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsi27, jal24, scratch)

    scratch = (&
      screened_rates(k_p_al26__si27)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsi27, jal26, scratch)

    scratch = (&
      -screened_rates(k_he4_si27__p_p30)*Y(jhe4)*state % rho - screened_rates(k_he4_si27__s31)* &
      Y(jhe4)*state % rho - screened_rates(k_p_si27__he4_al24)*Y(jp)*state % rho - &
      screened_rates(k_p_si27__p28)*Y(jp)*state % rho - &
      screened_rates(k_si27__al27__weak__wc12) - screened_rates(k_si27__he4_mg23) - &
      screened_rates(k_si27__p_al26) &
       )
    call set_jac_entry(state, jsi27, jsi27, scratch)

    scratch = (&
      screened_rates(k_p27__si27__weak__wc12) &
       )
    call set_jac_entry(state, jsi27, jp27, scratch)

    scratch = (&
      screened_rates(k_p28__p_si27) &
       )
    call set_jac_entry(state, jsi27, jp28, scratch)

    scratch = (&
      screened_rates(k_p_p30__he4_si27)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsi27, jp30, scratch)

    scratch = (&
      screened_rates(k_s28__p_si27__weak__wc12) &
       )
    call set_jac_entry(state, jsi27, js28, scratch)

    scratch = (&
      screened_rates(k_s31__he4_si27) &
       )
    call set_jac_entry(state, jsi27, js31, scratch)

    scratch = (&
      screened_rates(k_p_al27__si28)*Y(jal27)*state % rho + screened_rates(k_p_p31__he4_si28)* &
      Y(jp31)*state % rho - screened_rates(k_p_si28__he4_al25)*Y(jsi28)* &
      state % rho - screened_rates(k_p_si28__p29)*Y(jsi28)*state % rho &
       )
    call set_jac_entry(state, jsi28, jp, scratch)

    scratch = (&
      screened_rates(k_he4_al25__p_si28)*Y(jal25)*state % rho + screened_rates(k_he4_mg24__si28)* &
      Y(jmg24)*state % rho - screened_rates(k_he4_si28__p_p31)*Y(jsi28)* &
      state % rho - screened_rates(k_he4_si28__s32)*Y(jsi28)*state % rho &
       )
    call set_jac_entry(state, jsi28, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_mg24__si28)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsi28, jmg24, scratch)

    scratch = (&
      screened_rates(k_he4_al25__p_si28)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsi28, jal25, scratch)

    scratch = (&
      screened_rates(k_p_al27__si28)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsi28, jal27, scratch)

    scratch = (&
      -screened_rates(k_he4_si28__p_p31)*Y(jhe4)*state % rho - screened_rates(k_he4_si28__s32)* &
      Y(jhe4)*state % rho - screened_rates(k_p_si28__he4_al25)*Y(jp)*state % rho - &
      screened_rates(k_p_si28__p29)*Y(jp)*state % rho - screened_rates(k_si28__he4_mg24) - &
      screened_rates(k_si28__p_al27) &
       )
    call set_jac_entry(state, jsi28, jsi28, scratch)

    scratch = (&
      screened_rates(k_p28__si28__weak__wc12) &
       )
    call set_jac_entry(state, jsi28, jp28, scratch)

    scratch = (&
      screened_rates(k_p29__p_si28) &
       )
    call set_jac_entry(state, jsi28, jp29, scratch)

    scratch = (&
      screened_rates(k_p_p31__he4_si28)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsi28, jp31, scratch)

    scratch = (&
      screened_rates(k_s29__p_si28__weak__wc12) &
       )
    call set_jac_entry(state, jsi28, js29, scratch)

    scratch = (&
      screened_rates(k_s32__he4_si28) &
       )
    call set_jac_entry(state, jsi28, js32, scratch)

    scratch = (&
      screened_rates(k_cl32__he4_si28__weak__wc12) &
       )
    call set_jac_entry(state, jsi28, jcl32, scratch)

    scratch = (&
      screened_rates(k_p_al28__si29)*Y(jal28)*state % rho + screened_rates(k_p_p32__he4_si29)* &
      Y(jp32)*state % rho - screened_rates(k_p_si29__he4_al26)*Y(jsi29)* &
      state % rho - screened_rates(k_p_si29__p30)*Y(jsi29)*state % rho &
       )
    call set_jac_entry(state, jsi29, jp, scratch)

    scratch = (&
      screened_rates(k_he4_al26__p_si29)*Y(jal26)*state % rho + screened_rates(k_he4_mg25__si29)* &
      Y(jmg25)*state % rho - screened_rates(k_he4_si29__p_p32)*Y(jsi29)* &
      state % rho - screened_rates(k_he4_si29__s33)*Y(jsi29)*state % rho &
       )
    call set_jac_entry(state, jsi29, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_mg25__si29)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsi29, jmg25, scratch)

    scratch = (&
      screened_rates(k_he4_al26__p_si29)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsi29, jal26, scratch)

    scratch = (&
      screened_rates(k_p_al28__si29)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsi29, jal28, scratch)

    scratch = (&
      -screened_rates(k_he4_si29__p_p32)*Y(jhe4)*state % rho - screened_rates(k_he4_si29__s33)* &
      Y(jhe4)*state % rho - screened_rates(k_p_si29__he4_al26)*Y(jp)*state % rho - &
      screened_rates(k_p_si29__p30)*Y(jp)*state % rho - screened_rates(k_si29__he4_mg25) - &
      screened_rates(k_si29__p_al28) &
       )
    call set_jac_entry(state, jsi29, jsi29, scratch)

    scratch = (&
      screened_rates(k_p29__si29__weak__wc12) &
       )
    call set_jac_entry(state, jsi29, jp29, scratch)

    scratch = (&
      screened_rates(k_p30__p_si29) &
       )
    call set_jac_entry(state, jsi29, jp30, scratch)

    scratch = (&
      screened_rates(k_p_p32__he4_si29)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsi29, jp32, scratch)

    scratch = (&
      screened_rates(k_s33__he4_si29) &
       )
    call set_jac_entry(state, jsi29, js33, scratch)

    scratch = (&
      screened_rates(k_p_p33__he4_si30)*Y(jp33)*state % rho - screened_rates(k_p_si30__he4_al27)* &
      Y(jsi30)*state % rho - screened_rates(k_p_si30__p31)*Y(jsi30)*state % rho &
       )
    call set_jac_entry(state, jsi30, jp, scratch)

    scratch = (&
      screened_rates(k_he4_al27__p_si30)*Y(jal27)*state % rho + screened_rates(k_he4_mg26__si30)* &
      Y(jmg26)*state % rho - screened_rates(k_he4_si30__p_p33)*Y(jsi30)* &
      state % rho - screened_rates(k_he4_si30__s34)*Y(jsi30)*state % rho &
       )
    call set_jac_entry(state, jsi30, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_mg26__si30)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsi30, jmg26, scratch)

    scratch = (&
      screened_rates(k_he4_al27__p_si30)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsi30, jal27, scratch)

    scratch = (&
      -screened_rates(k_he4_si30__p_p33)*Y(jhe4)*state % rho - screened_rates(k_he4_si30__s34)* &
      Y(jhe4)*state % rho - screened_rates(k_p_si30__he4_al27)*Y(jp)*state % rho - &
      screened_rates(k_p_si30__p31)*Y(jp)*state % rho - screened_rates(k_si30__he4_mg26) &
       )
    call set_jac_entry(state, jsi30, jsi30, scratch)

    scratch = (&
      screened_rates(k_p30__si30__weak__wc12) &
       )
    call set_jac_entry(state, jsi30, jp30, scratch)

    scratch = (&
      screened_rates(k_p31__p_si30) &
       )
    call set_jac_entry(state, jsi30, jp31, scratch)

    scratch = (&
      screened_rates(k_p_p33__he4_si30)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsi30, jp33, scratch)

    scratch = (&
      screened_rates(k_s34__he4_si30) &
       )
    call set_jac_entry(state, jsi30, js34, scratch)

    scratch = (&
      screened_rates(k_p_s29__he4_p26)*Y(js29)*state % rho + screened_rates(k_p_si25__p26)* &
      Y(jsi25)*state % rho &
       )
    call set_jac_entry(state, jp26, jp, scratch)

    scratch = (&
      screened_rates(k_he4_al22__p26)*Y(jal22)*state % rho - screened_rates(k_he4_p26__cl30)* &
      Y(jp26)*state % rho - screened_rates(k_he4_p26__p_s29)*Y(jp26)*state % rho &
       )
    call set_jac_entry(state, jp26, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_al22__p26)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp26, jal22, scratch)

    scratch = (&
      screened_rates(k_p_si25__p26)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp26, jsi25, scratch)

    scratch = (&
      -screened_rates(k_he4_p26__cl30)*Y(jhe4)*state % rho - screened_rates(k_he4_p26__p_s29)* &
      Y(jhe4)*state % rho - screened_rates(k_p26__he4_al22) - &
      screened_rates(k_p26__p_si25) - screened_rates(k_p26__si26__weak__wc12) &
       )
    call set_jac_entry(state, jp26, jp26, scratch)

    scratch = (&
      screened_rates(k_p_s29__he4_p26)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp26, js29, scratch)

    scratch = (&
      screened_rates(k_cl30__he4_p26) &
       )
    call set_jac_entry(state, jp26, jcl30, scratch)

    scratch = (&
      -screened_rates(k_p_p27__he4_si24)*Y(jp27)*state % rho - screened_rates(k_p_p27__s28)* &
      Y(jp27)*state % rho + screened_rates(k_p_s30__he4_p27)*Y(js30)*state % rho + &
      screened_rates(k_p_si26__p27)*Y(jsi26)*state % rho &
       )
    call set_jac_entry(state, jp27, jp, scratch)

    scratch = (&
      screened_rates(k_he4_al23__p27)*Y(jal23)*state % rho - screened_rates(k_he4_p27__cl31)* &
      Y(jp27)*state % rho - screened_rates(k_he4_p27__p_s30)*Y(jp27)*state % rho + &
      screened_rates(k_he4_si24__p_p27)*Y(jsi24)*state % rho &
       )
    call set_jac_entry(state, jp27, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_al23__p27)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp27, jal23, scratch)

    scratch = (&
      screened_rates(k_he4_si24__p_p27)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp27, jsi24, scratch)

    scratch = (&
      screened_rates(k_p_si26__p27)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp27, jsi26, scratch)

    scratch = (&
      -screened_rates(k_he4_p27__cl31)*Y(jhe4)*state % rho - screened_rates(k_he4_p27__p_s30)* &
      Y(jhe4)*state % rho - screened_rates(k_p27__he4_al23) - &
      screened_rates(k_p27__p_al26__weak__wc12) - screened_rates(k_p27__p_si26) - &
      screened_rates(k_p27__si27__weak__wc12) - screened_rates(k_p_p27__he4_si24)*Y(jp) &
      *state % rho - screened_rates(k_p_p27__s28)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp27, jp27, scratch)

    scratch = (&
      screened_rates(k_s28__p_p27) &
       )
    call set_jac_entry(state, jp27, js28, scratch)

    scratch = (&
      screened_rates(k_p_s30__he4_p27)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp27, js30, scratch)

    scratch = (&
      screened_rates(k_cl31__he4_p27) &
       )
    call set_jac_entry(state, jp27, jcl31, scratch)

    scratch = (&
      -screened_rates(k_p_p28__he4_si25)*Y(jp28)*state % rho - screened_rates(k_p_p28__s29)* &
      Y(jp28)*state % rho + screened_rates(k_p_s31__he4_p28)*Y(js31)*state % rho + &
      screened_rates(k_p_si27__p28)*Y(jsi27)*state % rho &
       )
    call set_jac_entry(state, jp28, jp, scratch)

    scratch = (&
      screened_rates(k_he4_al24__p28)*Y(jal24)*state % rho - screened_rates(k_he4_p28__cl32)* &
      Y(jp28)*state % rho - screened_rates(k_he4_p28__p_s31)*Y(jp28)*state % rho + &
      screened_rates(k_he4_si25__p_p28)*Y(jsi25)*state % rho &
       )
    call set_jac_entry(state, jp28, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_al24__p28)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp28, jal24, scratch)

    scratch = (&
      screened_rates(k_he4_si25__p_p28)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp28, jsi25, scratch)

    scratch = (&
      screened_rates(k_p_si27__p28)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp28, jsi27, scratch)

    scratch = (&
      -screened_rates(k_he4_p28__cl32)*Y(jhe4)*state % rho - screened_rates(k_he4_p28__p_s31)* &
      Y(jhe4)*state % rho - screened_rates(k_p28__he4_al24) - &
      screened_rates(k_p28__he4_mg24__weak__wc12) - &
      screened_rates(k_p28__p_al27__weak__wc12) - screened_rates(k_p28__p_si27) - &
      screened_rates(k_p28__si28__weak__wc12) - screened_rates(k_p_p28__he4_si25)*Y(jp) &
      *state % rho - screened_rates(k_p_p28__s29)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp28, jp28, scratch)

    scratch = (&
      screened_rates(k_s28__p28__weak__wc12) &
       )
    call set_jac_entry(state, jp28, js28, scratch)

    scratch = (&
      screened_rates(k_s29__p_p28) &
       )
    call set_jac_entry(state, jp28, js29, scratch)

    scratch = (&
      screened_rates(k_p_s31__he4_p28)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp28, js31, scratch)

    scratch = (&
      screened_rates(k_cl32__he4_p28) &
       )
    call set_jac_entry(state, jp28, jcl32, scratch)

    scratch = (&
      -screened_rates(k_p_p29__he4_si26)*Y(jp29)*state % rho - screened_rates(k_p_p29__s30)* &
      Y(jp29)*state % rho + screened_rates(k_p_s32__he4_p29)*Y(js32)*state % rho + &
      screened_rates(k_p_si28__p29)*Y(jsi28)*state % rho &
       )
    call set_jac_entry(state, jp29, jp, scratch)

    scratch = (&
      screened_rates(k_he4_al25__p29)*Y(jal25)*state % rho - screened_rates(k_he4_p29__cl33)* &
      Y(jp29)*state % rho - screened_rates(k_he4_p29__p_s32)*Y(jp29)*state % rho + &
      screened_rates(k_he4_si26__p_p29)*Y(jsi26)*state % rho &
       )
    call set_jac_entry(state, jp29, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_al25__p29)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp29, jal25, scratch)

    scratch = (&
      screened_rates(k_he4_si26__p_p29)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp29, jsi26, scratch)

    scratch = (&
      screened_rates(k_p_si28__p29)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp29, jsi28, scratch)

    scratch = (&
      -screened_rates(k_he4_p29__cl33)*Y(jhe4)*state % rho - screened_rates(k_he4_p29__p_s32)* &
      Y(jhe4)*state % rho - screened_rates(k_p29__he4_al25) - &
      screened_rates(k_p29__p_si28) - screened_rates(k_p29__si29__weak__wc12) - &
      screened_rates(k_p_p29__he4_si26)*Y(jp)*state % rho - screened_rates(k_p_p29__s30)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp29, jp29, scratch)

    scratch = (&
      screened_rates(k_s29__p29__weak__wc12) &
       )
    call set_jac_entry(state, jp29, js29, scratch)

    scratch = (&
      screened_rates(k_s30__p_p29) &
       )
    call set_jac_entry(state, jp29, js30, scratch)

    scratch = (&
      screened_rates(k_p_s32__he4_p29)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp29, js32, scratch)

    scratch = (&
      screened_rates(k_cl33__he4_p29) &
       )
    call set_jac_entry(state, jp29, jcl33, scratch)

    scratch = (&
      screened_rates(k_ar31__p_p_p29__weak__wc12) &
       )
    call set_jac_entry(state, jp29, jar31, scratch)

    scratch = (&
      -screened_rates(k_p_p30__he4_si27)*Y(jp30)*state % rho - screened_rates(k_p_p30__s31)* &
      Y(jp30)*state % rho + screened_rates(k_p_s33__he4_p30)*Y(js33)*state % rho + &
      screened_rates(k_p_si29__p30)*Y(jsi29)*state % rho &
       )
    call set_jac_entry(state, jp30, jp, scratch)

    scratch = (&
      screened_rates(k_he4_al26__p30)*Y(jal26)*state % rho - screened_rates(k_he4_p30__cl34)* &
      Y(jp30)*state % rho - screened_rates(k_he4_p30__p_s33)*Y(jp30)*state % rho + &
      screened_rates(k_he4_si27__p_p30)*Y(jsi27)*state % rho &
       )
    call set_jac_entry(state, jp30, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_al26__p30)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp30, jal26, scratch)

    scratch = (&
      screened_rates(k_he4_si27__p_p30)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp30, jsi27, scratch)

    scratch = (&
      screened_rates(k_p_si29__p30)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp30, jsi29, scratch)

    scratch = (&
      -screened_rates(k_he4_p30__cl34)*Y(jhe4)*state % rho - screened_rates(k_he4_p30__p_s33)* &
      Y(jhe4)*state % rho - screened_rates(k_p30__he4_al26) - &
      screened_rates(k_p30__p_si29) - screened_rates(k_p30__si30__weak__wc12) - &
      screened_rates(k_p_p30__he4_si27)*Y(jp)*state % rho - screened_rates(k_p_p30__s31)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp30, jp30, scratch)

    scratch = (&
      screened_rates(k_s30__p30__weak__wc12) &
       )
    call set_jac_entry(state, jp30, js30, scratch)

    scratch = (&
      screened_rates(k_s31__p_p30) &
       )
    call set_jac_entry(state, jp30, js31, scratch)

    scratch = (&
      screened_rates(k_p_s33__he4_p30)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp30, js33, scratch)

    scratch = (&
      screened_rates(k_cl31__p_p30__weak__wc12) &
       )
    call set_jac_entry(state, jp30, jcl31, scratch)

    scratch = (&
      screened_rates(k_cl34__he4_p30) &
       )
    call set_jac_entry(state, jp30, jcl34, scratch)

    scratch = (&
      -screened_rates(k_p_p31__he4_si28)*Y(jp31)*state % rho - screened_rates(k_p_p31__s32)* &
      Y(jp31)*state % rho + screened_rates(k_p_s34__he4_p31)*Y(js34)*state % rho + &
      screened_rates(k_p_si30__p31)*Y(jsi30)*state % rho &
       )
    call set_jac_entry(state, jp31, jp, scratch)

    scratch = (&
      screened_rates(k_he4_al27__p31)*Y(jal27)*state % rho - screened_rates(k_he4_p31__cl35)* &
      Y(jp31)*state % rho - screened_rates(k_he4_p31__p_s34)*Y(jp31)*state % rho + &
      screened_rates(k_he4_si28__p_p31)*Y(jsi28)*state % rho &
       )
    call set_jac_entry(state, jp31, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_al27__p31)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp31, jal27, scratch)

    scratch = (&
      screened_rates(k_he4_si28__p_p31)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp31, jsi28, scratch)

    scratch = (&
      screened_rates(k_p_si30__p31)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp31, jsi30, scratch)

    scratch = (&
      -screened_rates(k_he4_p31__cl35)*Y(jhe4)*state % rho - screened_rates(k_he4_p31__p_s34)* &
      Y(jhe4)*state % rho - screened_rates(k_p31__he4_al27) - &
      screened_rates(k_p31__p_si30) - screened_rates(k_p_p31__he4_si28)*Y(jp)*state % rho &
      - screened_rates(k_p_p31__s32)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp31, jp31, scratch)

    scratch = (&
      screened_rates(k_s31__p31__weak__wc12) &
       )
    call set_jac_entry(state, jp31, js31, scratch)

    scratch = (&
      screened_rates(k_s32__p_p31) &
       )
    call set_jac_entry(state, jp31, js32, scratch)

    scratch = (&
      screened_rates(k_p_s34__he4_p31)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp31, js34, scratch)

    scratch = (&
      screened_rates(k_cl32__p_p31__weak__wc12) &
       )
    call set_jac_entry(state, jp31, jcl32, scratch)

    scratch = (&
      screened_rates(k_cl35__he4_p31) &
       )
    call set_jac_entry(state, jp31, jcl35, scratch)

    scratch = (&
      -screened_rates(k_p_p32__he4_si29)*Y(jp32)*state % rho - screened_rates(k_p_p32__s33)* &
      Y(jp32)*state % rho + screened_rates(k_p_s35__he4_p32)*Y(js35)*state % rho &
       )
    call set_jac_entry(state, jp32, jp, scratch)

    scratch = (&
      screened_rates(k_he4_al28__p32)*Y(jal28)*state % rho - screened_rates(k_he4_p32__cl36)* &
      Y(jp32)*state % rho - screened_rates(k_he4_p32__p_s35)*Y(jp32)*state % rho + &
      screened_rates(k_he4_si29__p_p32)*Y(jsi29)*state % rho &
       )
    call set_jac_entry(state, jp32, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_al28__p32)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp32, jal28, scratch)

    scratch = (&
      screened_rates(k_he4_si29__p_p32)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp32, jsi29, scratch)

    scratch = (&
      -screened_rates(k_he4_p32__cl36)*Y(jhe4)*state % rho - screened_rates(k_he4_p32__p_s35)* &
      Y(jhe4)*state % rho - screened_rates(k_p32__he4_al28) - &
      screened_rates(k_p_p32__he4_si29)*Y(jp)*state % rho - screened_rates(k_p_p32__s33)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp32, jp32, scratch)

    scratch = (&
      screened_rates(k_s33__p_p32) &
       )
    call set_jac_entry(state, jp32, js33, scratch)

    scratch = (&
      screened_rates(k_p_s35__he4_p32)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp32, js35, scratch)

    scratch = (&
      screened_rates(k_cl36__he4_p32) &
       )
    call set_jac_entry(state, jp32, jcl36, scratch)

    scratch = (&
      -screened_rates(k_p_p33__he4_si30)*Y(jp33)*state % rho - screened_rates(k_p_p33__s34)* &
      Y(jp33)*state % rho &
       )
    call set_jac_entry(state, jp33, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_p33__cl37)*Y(jp33)*state % rho + screened_rates(k_he4_si30__p_p33)* &
      Y(jsi30)*state % rho &
       )
    call set_jac_entry(state, jp33, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_si30__p_p33)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jp33, jsi30, scratch)

    scratch = (&
      -screened_rates(k_he4_p33__cl37)*Y(jhe4)*state % rho - screened_rates(k_p_p33__he4_si30)* &
      Y(jp)*state % rho - screened_rates(k_p_p33__s34)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jp33, jp33, scratch)

    scratch = (&
      screened_rates(k_s34__p_p33) &
       )
    call set_jac_entry(state, jp33, js34, scratch)

    scratch = (&
      screened_rates(k_cl37__he4_p33) &
       )
    call set_jac_entry(state, jp33, jcl37, scratch)

    scratch = (&
      screened_rates(k_p_cl31__he4_s28)*Y(jcl31)*state % rho + screened_rates(k_p_p27__s28)* &
      Y(jp27)*state % rho - screened_rates(k_p_s28__cl29)*Y(js28)*state % rho &
       )
    call set_jac_entry(state, js28, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_s28__ar32)*Y(js28)*state % rho - screened_rates(k_he4_s28__p_cl31)* &
      Y(js28)*state % rho + screened_rates(k_he4_si24__s28)*Y(jsi24)*state % rho &
       )
    call set_jac_entry(state, js28, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_si24__s28)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, js28, jsi24, scratch)

    scratch = (&
      screened_rates(k_p_p27__s28)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, js28, jp27, scratch)

    scratch = (&
      -screened_rates(k_he4_s28__ar32)*Y(jhe4)*state % rho - screened_rates(k_he4_s28__p_cl31)* &
      Y(jhe4)*state % rho - screened_rates(k_p_s28__cl29)*Y(jp)*state % rho - &
      screened_rates(k_s28__he4_si24) - screened_rates(k_s28__p28__weak__wc12) - &
      screened_rates(k_s28__p_p27) - screened_rates(k_s28__p_si27__weak__wc12) &
       )
    call set_jac_entry(state, js28, js28, scratch)

    scratch = (&
      screened_rates(k_cl29__p_s28) &
       )
    call set_jac_entry(state, js28, jcl29, scratch)

    scratch = (&
      screened_rates(k_p_cl31__he4_s28)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, js28, jcl31, scratch)

    scratch = (&
      screened_rates(k_ar32__he4_s28) &
       )
    call set_jac_entry(state, js28, jar32, scratch)

    scratch = (&
      screened_rates(k_p_cl32__he4_s29)*Y(jcl32)*state % rho + screened_rates(k_p_p28__s29)* &
      Y(jp28)*state % rho - screened_rates(k_p_s29__cl30)*Y(js29)*state % rho - &
      screened_rates(k_p_s29__he4_p26)*Y(js29)*state % rho &
       )
    call set_jac_entry(state, js29, jp, scratch)

    scratch = (&
      screened_rates(k_he4_p26__p_s29)*Y(jp26)*state % rho - screened_rates(k_he4_s29__ar33)* &
      Y(js29)*state % rho - screened_rates(k_he4_s29__p_cl32)*Y(js29)*state % rho &
      + screened_rates(k_he4_si25__s29)*Y(jsi25)*state % rho &
       )
    call set_jac_entry(state, js29, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_si25__s29)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, js29, jsi25, scratch)

    scratch = (&
      screened_rates(k_he4_p26__p_s29)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, js29, jp26, scratch)

    scratch = (&
      screened_rates(k_p_p28__s29)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, js29, jp28, scratch)

    scratch = (&
      -screened_rates(k_he4_s29__ar33)*Y(jhe4)*state % rho - screened_rates(k_he4_s29__p_cl32)* &
      Y(jhe4)*state % rho - screened_rates(k_p_s29__cl30)*Y(jp)*state % rho - &
      screened_rates(k_p_s29__he4_p26)*Y(jp)*state % rho - screened_rates(k_s29__he4_si25) &
      - screened_rates(k_s29__p29__weak__wc12) - screened_rates(k_s29__p_p28) - &
      screened_rates(k_s29__p_si28__weak__wc12) &
       )
    call set_jac_entry(state, js29, js29, scratch)

    scratch = (&
      screened_rates(k_cl29__s29__weak__bqa_pos_) &
       )
    call set_jac_entry(state, js29, jcl29, scratch)

    scratch = (&
      screened_rates(k_cl30__p_s29) &
       )
    call set_jac_entry(state, js29, jcl30, scratch)

    scratch = (&
      screened_rates(k_p_cl32__he4_s29)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, js29, jcl32, scratch)

    scratch = (&
      screened_rates(k_ar33__he4_s29) &
       )
    call set_jac_entry(state, js29, jar33, scratch)

    scratch = (&
      screened_rates(k_p_cl33__he4_s30)*Y(jcl33)*state % rho + screened_rates(k_p_p29__s30)* &
      Y(jp29)*state % rho - screened_rates(k_p_s30__cl31)*Y(js30)*state % rho - &
      screened_rates(k_p_s30__he4_p27)*Y(js30)*state % rho &
       )
    call set_jac_entry(state, js30, jp, scratch)

    scratch = (&
      screened_rates(k_he4_p27__p_s30)*Y(jp27)*state % rho - screened_rates(k_he4_s30__ar34)* &
      Y(js30)*state % rho - screened_rates(k_he4_s30__p_cl33)*Y(js30)*state % rho &
      + screened_rates(k_he4_si26__s30)*Y(jsi26)*state % rho &
       )
    call set_jac_entry(state, js30, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_si26__s30)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, js30, jsi26, scratch)

    scratch = (&
      screened_rates(k_he4_p27__p_s30)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, js30, jp27, scratch)

    scratch = (&
      screened_rates(k_p_p29__s30)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, js30, jp29, scratch)

    scratch = (&
      -screened_rates(k_he4_s30__ar34)*Y(jhe4)*state % rho - screened_rates(k_he4_s30__p_cl33)* &
      Y(jhe4)*state % rho - screened_rates(k_p_s30__cl31)*Y(jp)*state % rho - &
      screened_rates(k_p_s30__he4_p27)*Y(jp)*state % rho - screened_rates(k_s30__he4_si26) &
      - screened_rates(k_s30__p30__weak__wc12) - screened_rates(k_s30__p_p29) &
       )
    call set_jac_entry(state, js30, js30, scratch)

    scratch = (&
      screened_rates(k_cl30__s30__weak__bqa_pos_) &
       )
    call set_jac_entry(state, js30, jcl30, scratch)

    scratch = (&
      screened_rates(k_cl31__p_s30) &
       )
    call set_jac_entry(state, js30, jcl31, scratch)

    scratch = (&
      screened_rates(k_p_cl33__he4_s30)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, js30, jcl33, scratch)

    scratch = (&
      screened_rates(k_ar31__p_s30__weak__wc17) &
       )
    call set_jac_entry(state, js30, jar31, scratch)

    scratch = (&
      screened_rates(k_ar34__he4_s30) &
       )
    call set_jac_entry(state, js30, jar34, scratch)

    scratch = (&
      screened_rates(k_p_cl34__he4_s31)*Y(jcl34)*state % rho + screened_rates(k_p_p30__s31)* &
      Y(jp30)*state % rho - screened_rates(k_p_s31__cl32)*Y(js31)*state % rho - &
      screened_rates(k_p_s31__he4_p28)*Y(js31)*state % rho &
       )
    call set_jac_entry(state, js31, jp, scratch)

    scratch = (&
      screened_rates(k_he4_p28__p_s31)*Y(jp28)*state % rho - screened_rates(k_he4_s31__ar35)* &
      Y(js31)*state % rho - screened_rates(k_he4_s31__p_cl34)*Y(js31)*state % rho &
      + screened_rates(k_he4_si27__s31)*Y(jsi27)*state % rho &
       )
    call set_jac_entry(state, js31, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_si27__s31)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, js31, jsi27, scratch)

    scratch = (&
      screened_rates(k_he4_p28__p_s31)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, js31, jp28, scratch)

    scratch = (&
      screened_rates(k_p_p30__s31)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, js31, jp30, scratch)

    scratch = (&
      -screened_rates(k_he4_s31__ar35)*Y(jhe4)*state % rho - screened_rates(k_he4_s31__p_cl34)* &
      Y(jhe4)*state % rho - screened_rates(k_p_s31__cl32)*Y(jp)*state % rho - &
      screened_rates(k_p_s31__he4_p28)*Y(jp)*state % rho - screened_rates(k_s31__he4_si27) &
      - screened_rates(k_s31__p31__weak__wc12) - screened_rates(k_s31__p_p30) &
       )
    call set_jac_entry(state, js31, js31, scratch)

    scratch = (&
      screened_rates(k_cl31__s31__weak__wc12) &
       )
    call set_jac_entry(state, js31, jcl31, scratch)

    scratch = (&
      screened_rates(k_cl32__p_s31) &
       )
    call set_jac_entry(state, js31, jcl32, scratch)

    scratch = (&
      screened_rates(k_p_cl34__he4_s31)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, js31, jcl34, scratch)

    scratch = (&
      screened_rates(k_ar32__p_s31__weak__wc12) &
       )
    call set_jac_entry(state, js31, jar32, scratch)

    scratch = (&
      screened_rates(k_ar35__he4_s31) &
       )
    call set_jac_entry(state, js31, jar35, scratch)

    scratch = (&
      screened_rates(k_p_cl35__he4_s32)*Y(jcl35)*state % rho + screened_rates(k_p_p31__s32)* &
      Y(jp31)*state % rho - screened_rates(k_p_s32__cl33)*Y(js32)*state % rho - &
      screened_rates(k_p_s32__he4_p29)*Y(js32)*state % rho &
       )
    call set_jac_entry(state, js32, jp, scratch)

    scratch = (&
      screened_rates(k_he4_p29__p_s32)*Y(jp29)*state % rho - screened_rates(k_he4_s32__ar36)* &
      Y(js32)*state % rho - screened_rates(k_he4_s32__p_cl35)*Y(js32)*state % rho &
      + screened_rates(k_he4_si28__s32)*Y(jsi28)*state % rho &
       )
    call set_jac_entry(state, js32, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_si28__s32)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, js32, jsi28, scratch)

    scratch = (&
      screened_rates(k_he4_p29__p_s32)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, js32, jp29, scratch)

    scratch = (&
      screened_rates(k_p_p31__s32)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, js32, jp31, scratch)

    scratch = (&
      -screened_rates(k_he4_s32__ar36)*Y(jhe4)*state % rho - screened_rates(k_he4_s32__p_cl35)* &
      Y(jhe4)*state % rho - screened_rates(k_p_s32__cl33)*Y(jp)*state % rho - &
      screened_rates(k_p_s32__he4_p29)*Y(jp)*state % rho - screened_rates(k_s32__he4_si28) &
      - screened_rates(k_s32__p_p31) &
       )
    call set_jac_entry(state, js32, js32, scratch)

    scratch = (&
      screened_rates(k_cl32__s32__weak__wc12) &
       )
    call set_jac_entry(state, js32, jcl32, scratch)

    scratch = (&
      screened_rates(k_cl33__p_s32) &
       )
    call set_jac_entry(state, js32, jcl33, scratch)

    scratch = (&
      screened_rates(k_p_cl35__he4_s32)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, js32, jcl35, scratch)

    scratch = (&
      screened_rates(k_ar33__p_s32__weak__wc12) &
       )
    call set_jac_entry(state, js32, jar33, scratch)

    scratch = (&
      screened_rates(k_ar36__he4_s32) &
       )
    call set_jac_entry(state, js32, jar36, scratch)

    scratch = (&
      screened_rates(k_k36__he4_s32__weak__wc12) &
       )
    call set_jac_entry(state, js32, jk36, scratch)

    scratch = (&
      screened_rates(k_p_cl36__he4_s33)*Y(jcl36)*state % rho + screened_rates(k_p_p32__s33)* &
      Y(jp32)*state % rho - screened_rates(k_p_s33__cl34)*Y(js33)*state % rho - &
      screened_rates(k_p_s33__he4_p30)*Y(js33)*state % rho &
       )
    call set_jac_entry(state, js33, jp, scratch)

    scratch = (&
      screened_rates(k_he4_p30__p_s33)*Y(jp30)*state % rho - screened_rates(k_he4_s33__ar37)* &
      Y(js33)*state % rho - screened_rates(k_he4_s33__p_cl36)*Y(js33)*state % rho &
      + screened_rates(k_he4_si29__s33)*Y(jsi29)*state % rho &
       )
    call set_jac_entry(state, js33, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_si29__s33)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, js33, jsi29, scratch)

    scratch = (&
      screened_rates(k_he4_p30__p_s33)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, js33, jp30, scratch)

    scratch = (&
      screened_rates(k_p_p32__s33)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, js33, jp32, scratch)

    scratch = (&
      -screened_rates(k_he4_s33__ar37)*Y(jhe4)*state % rho - screened_rates(k_he4_s33__p_cl36)* &
      Y(jhe4)*state % rho - screened_rates(k_p_s33__cl34)*Y(jp)*state % rho - &
      screened_rates(k_p_s33__he4_p30)*Y(jp)*state % rho - screened_rates(k_s33__he4_si29) &
      - screened_rates(k_s33__p_p32) &
       )
    call set_jac_entry(state, js33, js33, scratch)

    scratch = (&
      screened_rates(k_cl33__s33__weak__wc12) &
       )
    call set_jac_entry(state, js33, jcl33, scratch)

    scratch = (&
      screened_rates(k_cl34__p_s33) &
       )
    call set_jac_entry(state, js33, jcl34, scratch)

    scratch = (&
      screened_rates(k_p_cl36__he4_s33)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, js33, jcl36, scratch)

    scratch = (&
      screened_rates(k_ar37__he4_s33) &
       )
    call set_jac_entry(state, js33, jar37, scratch)

    scratch = (&
      screened_rates(k_p_cl37__he4_s34)*Y(jcl37)*state % rho + screened_rates(k_p_p33__s34)* &
      Y(jp33)*state % rho - screened_rates(k_p_s34__cl35)*Y(js34)*state % rho - &
      screened_rates(k_p_s34__he4_p31)*Y(js34)*state % rho &
       )
    call set_jac_entry(state, js34, jp, scratch)

    scratch = (&
      screened_rates(k_he4_p31__p_s34)*Y(jp31)*state % rho - screened_rates(k_he4_s34__ar38)* &
      Y(js34)*state % rho - screened_rates(k_he4_s34__p_cl37)*Y(js34)*state % rho &
      + screened_rates(k_he4_si30__s34)*Y(jsi30)*state % rho &
       )
    call set_jac_entry(state, js34, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_si30__s34)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, js34, jsi30, scratch)

    scratch = (&
      screened_rates(k_he4_p31__p_s34)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, js34, jp31, scratch)

    scratch = (&
      screened_rates(k_p_p33__s34)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, js34, jp33, scratch)

    scratch = (&
      -screened_rates(k_he4_s34__ar38)*Y(jhe4)*state % rho - screened_rates(k_he4_s34__p_cl37)* &
      Y(jhe4)*state % rho - screened_rates(k_p_s34__cl35)*Y(jp)*state % rho - &
      screened_rates(k_p_s34__he4_p31)*Y(jp)*state % rho - screened_rates(k_s34__he4_si30) &
      - screened_rates(k_s34__p_p33) &
       )
    call set_jac_entry(state, js34, js34, scratch)

    scratch = (&
      screened_rates(k_cl34__s34__weak__wc12) &
       )
    call set_jac_entry(state, js34, jcl34, scratch)

    scratch = (&
      screened_rates(k_cl35__p_s34) &
       )
    call set_jac_entry(state, js34, jcl35, scratch)

    scratch = (&
      screened_rates(k_p_cl37__he4_s34)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, js34, jcl37, scratch)

    scratch = (&
      screened_rates(k_ar38__he4_s34) &
       )
    call set_jac_entry(state, js34, jar38, scratch)

    scratch = (&
      -screened_rates(k_p_s35__cl36)*Y(js35)*state % rho - screened_rates(k_p_s35__he4_p32)* &
      Y(js35)*state % rho &
       )
    call set_jac_entry(state, js35, jp, scratch)

    scratch = (&
      screened_rates(k_he4_p32__p_s35)*Y(jp32)*state % rho - screened_rates(k_he4_s35__ar39)* &
      Y(js35)*state % rho &
       )
    call set_jac_entry(state, js35, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_p32__p_s35)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, js35, jp32, scratch)

    scratch = (&
      -screened_rates(k_he4_s35__ar39)*Y(jhe4)*state % rho - screened_rates(k_p_s35__cl36)* &
      Y(jp)*state % rho - screened_rates(k_p_s35__he4_p32)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, js35, js35, scratch)

    scratch = (&
      screened_rates(k_cl36__p_s35) &
       )
    call set_jac_entry(state, js35, jcl36, scratch)

    scratch = (&
      screened_rates(k_ar39__he4_s35) &
       )
    call set_jac_entry(state, js35, jar39, scratch)

    scratch = (&
      screened_rates(k_p_ar32__he4_cl29)*Y(jar32)*state % rho + screened_rates(k_p_s28__cl29)* &
      Y(js28)*state % rho &
       )
    call set_jac_entry(state, jcl29, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cl29__k33)*Y(jcl29)*state % rho - screened_rates(k_he4_cl29__p_ar32)* &
      Y(jcl29)*state % rho &
       )
    call set_jac_entry(state, jcl29, jhe4, scratch)

    scratch = (&
      screened_rates(k_p_s28__cl29)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl29, js28, scratch)

    scratch = (&
      -screened_rates(k_cl29__p_s28) - screened_rates(k_cl29__s29__weak__bqa_pos_) - &
      screened_rates(k_he4_cl29__k33)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl29__p_ar32)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcl29, jcl29, scratch)

    scratch = (&
      screened_rates(k_p_ar32__he4_cl29)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl29, jar32, scratch)

    scratch = (&
      screened_rates(k_k33__he4_cl29) &
       )
    call set_jac_entry(state, jcl29, jk33, scratch)

    scratch = (&
      screened_rates(k_p_ar33__he4_cl30)*Y(jar33)*state % rho - screened_rates(k_p_cl30__ar31)* &
      Y(jcl30)*state % rho + screened_rates(k_p_s29__cl30)*Y(js29)*state % rho &
       )
    call set_jac_entry(state, jcl30, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cl30__k34)*Y(jcl30)*state % rho - screened_rates(k_he4_cl30__p_ar33)* &
      Y(jcl30)*state % rho + screened_rates(k_he4_p26__cl30)*Y(jp26)*state % rho &
       )
    call set_jac_entry(state, jcl30, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_p26__cl30)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcl30, jp26, scratch)

    scratch = (&
      screened_rates(k_p_s29__cl30)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl30, js29, scratch)

    scratch = (&
      -screened_rates(k_cl30__he4_p26) - screened_rates(k_cl30__p_s29) - &
      screened_rates(k_cl30__s30__weak__bqa_pos_) - screened_rates(k_he4_cl30__k34)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_cl30__p_ar33)*Y(jhe4)*state % rho &
      - screened_rates(k_p_cl30__ar31)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl30, jcl30, scratch)

    scratch = (&
      screened_rates(k_ar31__p_cl30) &
       )
    call set_jac_entry(state, jcl30, jar31, scratch)

    scratch = (&
      screened_rates(k_p_ar33__he4_cl30)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl30, jar33, scratch)

    scratch = (&
      screened_rates(k_k34__he4_cl30) &
       )
    call set_jac_entry(state, jcl30, jk34, scratch)

    scratch = (&
      screened_rates(k_p_ar34__he4_cl31)*Y(jar34)*state % rho - screened_rates(k_p_cl31__ar32)* &
      Y(jcl31)*state % rho - screened_rates(k_p_cl31__he4_s28)*Y(jcl31)* &
      state % rho + screened_rates(k_p_s30__cl31)*Y(js30)*state % rho &
       )
    call set_jac_entry(state, jcl31, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cl31__k35)*Y(jcl31)*state % rho - screened_rates(k_he4_cl31__p_ar34)* &
      Y(jcl31)*state % rho + screened_rates(k_he4_p27__cl31)*Y(jp27)*state % rho + &
      screened_rates(k_he4_s28__p_cl31)*Y(js28)*state % rho &
       )
    call set_jac_entry(state, jcl31, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_p27__cl31)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcl31, jp27, scratch)

    scratch = (&
      screened_rates(k_he4_s28__p_cl31)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcl31, js28, scratch)

    scratch = (&
      screened_rates(k_p_s30__cl31)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl31, js30, scratch)

    scratch = (&
      -screened_rates(k_cl31__he4_p27) - screened_rates(k_cl31__p_p30__weak__wc12) - &
      screened_rates(k_cl31__p_s30) - screened_rates(k_cl31__s31__weak__wc12) - &
      screened_rates(k_he4_cl31__k35)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl31__p_ar34)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cl31__ar32)*Y(jp)*state % rho - screened_rates(k_p_cl31__he4_s28) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl31, jcl31, scratch)

    scratch = (&
      screened_rates(k_ar31__cl31__weak__wc12) &
       )
    call set_jac_entry(state, jcl31, jar31, scratch)

    scratch = (&
      screened_rates(k_ar32__p_cl31) &
       )
    call set_jac_entry(state, jcl31, jar32, scratch)

    scratch = (&
      screened_rates(k_p_ar34__he4_cl31)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl31, jar34, scratch)

    scratch = (&
      screened_rates(k_k35__he4_cl31) &
       )
    call set_jac_entry(state, jcl31, jk35, scratch)

    scratch = (&
      screened_rates(k_p_ar35__he4_cl32)*Y(jar35)*state % rho - screened_rates(k_p_cl32__ar33)* &
      Y(jcl32)*state % rho - screened_rates(k_p_cl32__he4_s29)*Y(jcl32)* &
      state % rho + screened_rates(k_p_s31__cl32)*Y(js31)*state % rho &
       )
    call set_jac_entry(state, jcl32, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cl32__k36)*Y(jcl32)*state % rho - screened_rates(k_he4_cl32__p_ar35)* &
      Y(jcl32)*state % rho + screened_rates(k_he4_p28__cl32)*Y(jp28)*state % rho + &
      screened_rates(k_he4_s29__p_cl32)*Y(js29)*state % rho &
       )
    call set_jac_entry(state, jcl32, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_p28__cl32)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcl32, jp28, scratch)

    scratch = (&
      screened_rates(k_he4_s29__p_cl32)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcl32, js29, scratch)

    scratch = (&
      screened_rates(k_p_s31__cl32)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl32, js31, scratch)

    scratch = (&
      -screened_rates(k_cl32__he4_p28) - screened_rates(k_cl32__he4_si28__weak__wc12) - &
      screened_rates(k_cl32__p_p31__weak__wc12) - screened_rates(k_cl32__p_s31) - &
      screened_rates(k_cl32__s32__weak__wc12) - screened_rates(k_he4_cl32__k36)*Y(jhe4) &
      *state % rho - screened_rates(k_he4_cl32__p_ar35)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cl32__ar33)*Y(jp)*state % rho - screened_rates(k_p_cl32__he4_s29) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl32, jcl32, scratch)

    scratch = (&
      screened_rates(k_ar32__cl32__weak__wc12) &
       )
    call set_jac_entry(state, jcl32, jar32, scratch)

    scratch = (&
      screened_rates(k_ar33__p_cl32) &
       )
    call set_jac_entry(state, jcl32, jar33, scratch)

    scratch = (&
      screened_rates(k_p_ar35__he4_cl32)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl32, jar35, scratch)

    scratch = (&
      screened_rates(k_k36__he4_cl32) &
       )
    call set_jac_entry(state, jcl32, jk36, scratch)

    scratch = (&
      screened_rates(k_p_ar36__he4_cl33)*Y(jar36)*state % rho - screened_rates(k_p_cl33__ar34)* &
      Y(jcl33)*state % rho - screened_rates(k_p_cl33__he4_s30)*Y(jcl33)* &
      state % rho + screened_rates(k_p_s32__cl33)*Y(js32)*state % rho &
       )
    call set_jac_entry(state, jcl33, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cl33__k37)*Y(jcl33)*state % rho - screened_rates(k_he4_cl33__p_ar36)* &
      Y(jcl33)*state % rho + screened_rates(k_he4_p29__cl33)*Y(jp29)*state % rho + &
      screened_rates(k_he4_s30__p_cl33)*Y(js30)*state % rho &
       )
    call set_jac_entry(state, jcl33, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_p29__cl33)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcl33, jp29, scratch)

    scratch = (&
      screened_rates(k_he4_s30__p_cl33)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcl33, js30, scratch)

    scratch = (&
      screened_rates(k_p_s32__cl33)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl33, js32, scratch)

    scratch = (&
      -screened_rates(k_cl33__he4_p29) - screened_rates(k_cl33__p_s32) - &
      screened_rates(k_cl33__s33__weak__wc12) - screened_rates(k_he4_cl33__k37)*Y(jhe4) &
      *state % rho - screened_rates(k_he4_cl33__p_ar36)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cl33__ar34)*Y(jp)*state % rho - screened_rates(k_p_cl33__he4_s30) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl33, jcl33, scratch)

    scratch = (&
      screened_rates(k_ar33__cl33__weak__wc12) &
       )
    call set_jac_entry(state, jcl33, jar33, scratch)

    scratch = (&
      screened_rates(k_ar34__p_cl33) &
       )
    call set_jac_entry(state, jcl33, jar34, scratch)

    scratch = (&
      screened_rates(k_p_ar36__he4_cl33)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl33, jar36, scratch)

    scratch = (&
      screened_rates(k_k37__he4_cl33) &
       )
    call set_jac_entry(state, jcl33, jk37, scratch)

    scratch = (&
      screened_rates(k_ca35__p_p_cl33__weak__wc12) &
       )
    call set_jac_entry(state, jcl33, jca35, scratch)

    scratch = (&
      screened_rates(k_p_ar37__he4_cl34)*Y(jar37)*state % rho - screened_rates(k_p_cl34__ar35)* &
      Y(jcl34)*state % rho - screened_rates(k_p_cl34__he4_s31)*Y(jcl34)* &
      state % rho + screened_rates(k_p_s33__cl34)*Y(js33)*state % rho &
       )
    call set_jac_entry(state, jcl34, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cl34__k38)*Y(jcl34)*state % rho - screened_rates(k_he4_cl34__p_ar37)* &
      Y(jcl34)*state % rho + screened_rates(k_he4_p30__cl34)*Y(jp30)*state % rho + &
      screened_rates(k_he4_s31__p_cl34)*Y(js31)*state % rho &
       )
    call set_jac_entry(state, jcl34, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_p30__cl34)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcl34, jp30, scratch)

    scratch = (&
      screened_rates(k_he4_s31__p_cl34)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcl34, js31, scratch)

    scratch = (&
      screened_rates(k_p_s33__cl34)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl34, js33, scratch)

    scratch = (&
      -screened_rates(k_cl34__he4_p30) - screened_rates(k_cl34__p_s33) - &
      screened_rates(k_cl34__s34__weak__wc12) - screened_rates(k_he4_cl34__k38)*Y(jhe4) &
      *state % rho - screened_rates(k_he4_cl34__p_ar37)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cl34__ar35)*Y(jp)*state % rho - screened_rates(k_p_cl34__he4_s31) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl34, jcl34, scratch)

    scratch = (&
      screened_rates(k_ar34__cl34__weak__wc12) &
       )
    call set_jac_entry(state, jcl34, jar34, scratch)

    scratch = (&
      screened_rates(k_ar35__p_cl34) &
       )
    call set_jac_entry(state, jcl34, jar35, scratch)

    scratch = (&
      screened_rates(k_p_ar37__he4_cl34)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl34, jar37, scratch)

    scratch = (&
      screened_rates(k_k35__p_cl34__weak__wc12) &
       )
    call set_jac_entry(state, jcl34, jk35, scratch)

    scratch = (&
      screened_rates(k_k38__he4_cl34) &
       )
    call set_jac_entry(state, jcl34, jk38, scratch)

    scratch = (&
      screened_rates(k_p_ar38__he4_cl35)*Y(jar38)*state % rho - screened_rates(k_p_cl35__ar36)* &
      Y(jcl35)*state % rho - screened_rates(k_p_cl35__he4_s32)*Y(jcl35)* &
      state % rho + screened_rates(k_p_s34__cl35)*Y(js34)*state % rho &
       )
    call set_jac_entry(state, jcl35, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cl35__k39)*Y(jcl35)*state % rho - screened_rates(k_he4_cl35__p_ar38)* &
      Y(jcl35)*state % rho + screened_rates(k_he4_p31__cl35)*Y(jp31)*state % rho + &
      screened_rates(k_he4_s32__p_cl35)*Y(js32)*state % rho &
       )
    call set_jac_entry(state, jcl35, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_p31__cl35)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcl35, jp31, scratch)

    scratch = (&
      screened_rates(k_he4_s32__p_cl35)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcl35, js32, scratch)

    scratch = (&
      screened_rates(k_p_s34__cl35)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl35, js34, scratch)

    scratch = (&
      -screened_rates(k_cl35__he4_p31) - screened_rates(k_cl35__p_s34) - screened_rates(k_he4_cl35__k39)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_cl35__p_ar38)*Y(jhe4)*state % rho &
      - screened_rates(k_p_cl35__ar36)*Y(jp)*state % rho - &
      screened_rates(k_p_cl35__he4_s32)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl35, jcl35, scratch)

    scratch = (&
      screened_rates(k_ar35__cl35__weak__wc12) &
       )
    call set_jac_entry(state, jcl35, jar35, scratch)

    scratch = (&
      screened_rates(k_ar36__p_cl35) &
       )
    call set_jac_entry(state, jcl35, jar36, scratch)

    scratch = (&
      screened_rates(k_p_ar38__he4_cl35)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl35, jar38, scratch)

    scratch = (&
      screened_rates(k_k36__p_cl35__weak__wc12) &
       )
    call set_jac_entry(state, jcl35, jk36, scratch)

    scratch = (&
      screened_rates(k_k39__he4_cl35) &
       )
    call set_jac_entry(state, jcl35, jk39, scratch)

    scratch = (&
      screened_rates(k_p_ar39__he4_cl36)*Y(jar39)*state % rho - screened_rates(k_p_cl36__ar37)* &
      Y(jcl36)*state % rho - screened_rates(k_p_cl36__he4_s33)*Y(jcl36)* &
      state % rho + screened_rates(k_p_s35__cl36)*Y(js35)*state % rho &
       )
    call set_jac_entry(state, jcl36, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cl36__k40)*Y(jcl36)*state % rho - screened_rates(k_he4_cl36__p_ar39)* &
      Y(jcl36)*state % rho + screened_rates(k_he4_p32__cl36)*Y(jp32)*state % rho + &
      screened_rates(k_he4_s33__p_cl36)*Y(js33)*state % rho &
       )
    call set_jac_entry(state, jcl36, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_p32__cl36)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcl36, jp32, scratch)

    scratch = (&
      screened_rates(k_he4_s33__p_cl36)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcl36, js33, scratch)

    scratch = (&
      screened_rates(k_p_s35__cl36)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl36, js35, scratch)

    scratch = (&
      -screened_rates(k_cl36__he4_p32) - screened_rates(k_cl36__p_s35) - screened_rates(k_he4_cl36__k40)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_cl36__p_ar39)*Y(jhe4)*state % rho &
      - screened_rates(k_p_cl36__ar37)*Y(jp)*state % rho - &
      screened_rates(k_p_cl36__he4_s33)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl36, jcl36, scratch)

    scratch = (&
      screened_rates(k_ar37__p_cl36) &
       )
    call set_jac_entry(state, jcl36, jar37, scratch)

    scratch = (&
      screened_rates(k_p_ar39__he4_cl36)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl36, jar39, scratch)

    scratch = (&
      screened_rates(k_k40__he4_cl36) &
       )
    call set_jac_entry(state, jcl36, jk40, scratch)

    scratch = (&
      -screened_rates(k_p_cl37__ar38)*Y(jcl37)*state % rho - screened_rates(k_p_cl37__he4_s34)* &
      Y(jcl37)*state % rho &
       )
    call set_jac_entry(state, jcl37, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cl37__k41)*Y(jcl37)*state % rho + screened_rates(k_he4_p33__cl37)* &
      Y(jp33)*state % rho + screened_rates(k_he4_s34__p_cl37)*Y(js34)*state % rho &
       )
    call set_jac_entry(state, jcl37, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_p33__cl37)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcl37, jp33, scratch)

    scratch = (&
      screened_rates(k_he4_s34__p_cl37)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcl37, js34, scratch)

    scratch = (&
      -screened_rates(k_cl37__he4_p33) - screened_rates(k_he4_cl37__k41)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cl37__ar38)*Y(jp)*state % rho - screened_rates(k_p_cl37__he4_s34) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcl37, jcl37, scratch)

    scratch = (&
      screened_rates(k_ar37__cl37__weak__wc12) &
       )
    call set_jac_entry(state, jcl37, jar37, scratch)

    scratch = (&
      screened_rates(k_ar38__p_cl37) &
       )
    call set_jac_entry(state, jcl37, jar38, scratch)

    scratch = (&
      screened_rates(k_k41__he4_cl37) &
       )
    call set_jac_entry(state, jcl37, jk41, scratch)

    scratch = (&
      screened_rates(k_p_cl30__ar31)*Y(jcl30)*state % rho + screened_rates(k_p_k34__he4_ar31)* &
      Y(jk34)*state % rho &
       )
    call set_jac_entry(state, jar31, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_ar31__ca35)*Y(jar31)*state % rho - screened_rates(k_he4_ar31__p_k34)* &
      Y(jar31)*state % rho &
       )
    call set_jac_entry(state, jar31, jhe4, scratch)

    scratch = (&
      screened_rates(k_p_cl30__ar31)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar31, jcl30, scratch)

    scratch = (&
      -screened_rates(k_ar31__cl31__weak__wc12) - screened_rates(k_ar31__p_cl30) - &
      screened_rates(k_ar31__p_p_p29__weak__wc12) - &
      screened_rates(k_ar31__p_s30__weak__wc17) - screened_rates(k_he4_ar31__ca35)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_ar31__p_k34)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jar31, jar31, scratch)

    scratch = (&
      screened_rates(k_p_k34__he4_ar31)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar31, jk34, scratch)

    scratch = (&
      screened_rates(k_ca35__he4_ar31) &
       )
    call set_jac_entry(state, jar31, jca35, scratch)

    scratch = (&
      -screened_rates(k_p_ar32__he4_cl29)*Y(jar32)*state % rho - screened_rates(k_p_ar32__k33)* &
      Y(jar32)*state % rho + screened_rates(k_p_cl31__ar32)*Y(jcl31)*state % rho + &
      screened_rates(k_p_k35__he4_ar32)*Y(jk35)*state % rho &
       )
    call set_jac_entry(state, jar32, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_ar32__ca36)*Y(jar32)*state % rho - screened_rates(k_he4_ar32__p_k35)* &
      Y(jar32)*state % rho + screened_rates(k_he4_cl29__p_ar32)*Y(jcl29)* &
      state % rho + screened_rates(k_he4_s28__ar32)*Y(js28)*state % rho &
       )
    call set_jac_entry(state, jar32, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_s28__ar32)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jar32, js28, scratch)

    scratch = (&
      screened_rates(k_he4_cl29__p_ar32)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jar32, jcl29, scratch)

    scratch = (&
      screened_rates(k_p_cl31__ar32)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar32, jcl31, scratch)

    scratch = (&
      -screened_rates(k_ar32__cl32__weak__wc12) - screened_rates(k_ar32__he4_s28) - &
      screened_rates(k_ar32__p_cl31) - screened_rates(k_ar32__p_s31__weak__wc12) - &
      screened_rates(k_he4_ar32__ca36)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar32__p_k35)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ar32__he4_cl29)*Y(jp)*state % rho - screened_rates(k_p_ar32__k33) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar32, jar32, scratch)

    scratch = (&
      screened_rates(k_k33__p_ar32) &
       )
    call set_jac_entry(state, jar32, jk33, scratch)

    scratch = (&
      screened_rates(k_p_k35__he4_ar32)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar32, jk35, scratch)

    scratch = (&
      screened_rates(k_ca36__he4_ar32) &
       )
    call set_jac_entry(state, jar32, jca36, scratch)

    scratch = (&
      -screened_rates(k_p_ar33__he4_cl30)*Y(jar33)*state % rho - screened_rates(k_p_ar33__k34)* &
      Y(jar33)*state % rho + screened_rates(k_p_cl32__ar33)*Y(jcl32)*state % rho + &
      screened_rates(k_p_k36__he4_ar33)*Y(jk36)*state % rho &
       )
    call set_jac_entry(state, jar33, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_ar33__ca37)*Y(jar33)*state % rho - screened_rates(k_he4_ar33__p_k36)* &
      Y(jar33)*state % rho + screened_rates(k_he4_cl30__p_ar33)*Y(jcl30)* &
      state % rho + screened_rates(k_he4_s29__ar33)*Y(js29)*state % rho &
       )
    call set_jac_entry(state, jar33, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_s29__ar33)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jar33, js29, scratch)

    scratch = (&
      screened_rates(k_he4_cl30__p_ar33)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jar33, jcl30, scratch)

    scratch = (&
      screened_rates(k_p_cl32__ar33)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar33, jcl32, scratch)

    scratch = (&
      -screened_rates(k_ar33__cl33__weak__wc12) - screened_rates(k_ar33__he4_s29) - &
      screened_rates(k_ar33__p_cl32) - screened_rates(k_ar33__p_s32__weak__wc12) - &
      screened_rates(k_he4_ar33__ca37)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar33__p_k36)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ar33__he4_cl30)*Y(jp)*state % rho - screened_rates(k_p_ar33__k34) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar33, jar33, scratch)

    scratch = (&
      screened_rates(k_k33__ar33__weak__bqa_pos_) &
       )
    call set_jac_entry(state, jar33, jk33, scratch)

    scratch = (&
      screened_rates(k_k34__p_ar33) &
       )
    call set_jac_entry(state, jar33, jk34, scratch)

    scratch = (&
      screened_rates(k_p_k36__he4_ar33)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar33, jk36, scratch)

    scratch = (&
      screened_rates(k_ca37__he4_ar33) &
       )
    call set_jac_entry(state, jar33, jca37, scratch)

    scratch = (&
      -screened_rates(k_p_ar34__he4_cl31)*Y(jar34)*state % rho - screened_rates(k_p_ar34__k35)* &
      Y(jar34)*state % rho + screened_rates(k_p_cl33__ar34)*Y(jcl33)*state % rho + &
      screened_rates(k_p_k37__he4_ar34)*Y(jk37)*state % rho &
       )
    call set_jac_entry(state, jar34, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_ar34__ca38)*Y(jar34)*state % rho - screened_rates(k_he4_ar34__p_k37)* &
      Y(jar34)*state % rho + screened_rates(k_he4_cl31__p_ar34)*Y(jcl31)* &
      state % rho + screened_rates(k_he4_s30__ar34)*Y(js30)*state % rho &
       )
    call set_jac_entry(state, jar34, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_s30__ar34)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jar34, js30, scratch)

    scratch = (&
      screened_rates(k_he4_cl31__p_ar34)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jar34, jcl31, scratch)

    scratch = (&
      screened_rates(k_p_cl33__ar34)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar34, jcl33, scratch)

    scratch = (&
      -screened_rates(k_ar34__cl34__weak__wc12) - screened_rates(k_ar34__he4_s30) - &
      screened_rates(k_ar34__p_cl33) - screened_rates(k_he4_ar34__ca38)*Y(jhe4)* &
      state % rho - screened_rates(k_he4_ar34__p_k37)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ar34__he4_cl31)*Y(jp)*state % rho - screened_rates(k_p_ar34__k35) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar34, jar34, scratch)

    scratch = (&
      screened_rates(k_k34__ar34__weak__bqa_pos_) &
       )
    call set_jac_entry(state, jar34, jk34, scratch)

    scratch = (&
      screened_rates(k_k35__p_ar34) &
       )
    call set_jac_entry(state, jar34, jk35, scratch)

    scratch = (&
      screened_rates(k_p_k37__he4_ar34)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar34, jk37, scratch)

    scratch = (&
      screened_rates(k_ca35__p_ar34__weak__wc17) &
       )
    call set_jac_entry(state, jar34, jca35, scratch)

    scratch = (&
      screened_rates(k_ca38__he4_ar34) &
       )
    call set_jac_entry(state, jar34, jca38, scratch)

    scratch = (&
      -screened_rates(k_p_ar35__he4_cl32)*Y(jar35)*state % rho - screened_rates(k_p_ar35__k36)* &
      Y(jar35)*state % rho + screened_rates(k_p_cl34__ar35)*Y(jcl34)*state % rho + &
      screened_rates(k_p_k38__he4_ar35)*Y(jk38)*state % rho &
       )
    call set_jac_entry(state, jar35, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_ar35__ca39)*Y(jar35)*state % rho - screened_rates(k_he4_ar35__p_k38)* &
      Y(jar35)*state % rho + screened_rates(k_he4_cl32__p_ar35)*Y(jcl32)* &
      state % rho + screened_rates(k_he4_s31__ar35)*Y(js31)*state % rho &
       )
    call set_jac_entry(state, jar35, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_s31__ar35)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jar35, js31, scratch)

    scratch = (&
      screened_rates(k_he4_cl32__p_ar35)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jar35, jcl32, scratch)

    scratch = (&
      screened_rates(k_p_cl34__ar35)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar35, jcl34, scratch)

    scratch = (&
      -screened_rates(k_ar35__cl35__weak__wc12) - screened_rates(k_ar35__he4_s31) - &
      screened_rates(k_ar35__p_cl34) - screened_rates(k_he4_ar35__ca39)*Y(jhe4)* &
      state % rho - screened_rates(k_he4_ar35__p_k38)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ar35__he4_cl32)*Y(jp)*state % rho - screened_rates(k_p_ar35__k36) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar35, jar35, scratch)

    scratch = (&
      screened_rates(k_k35__ar35__weak__wc12) &
       )
    call set_jac_entry(state, jar35, jk35, scratch)

    scratch = (&
      screened_rates(k_k36__p_ar35) &
       )
    call set_jac_entry(state, jar35, jk36, scratch)

    scratch = (&
      screened_rates(k_p_k38__he4_ar35)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar35, jk38, scratch)

    scratch = (&
      screened_rates(k_ca36__p_ar35__weak__wc12) &
       )
    call set_jac_entry(state, jar35, jca36, scratch)

    scratch = (&
      screened_rates(k_ca39__he4_ar35) &
       )
    call set_jac_entry(state, jar35, jca39, scratch)

    scratch = (&
      -screened_rates(k_p_ar36__he4_cl33)*Y(jar36)*state % rho - screened_rates(k_p_ar36__k37)* &
      Y(jar36)*state % rho + screened_rates(k_p_cl35__ar36)*Y(jcl35)*state % rho + &
      screened_rates(k_p_k39__he4_ar36)*Y(jk39)*state % rho &
       )
    call set_jac_entry(state, jar36, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_ar36__ca40)*Y(jar36)*state % rho - screened_rates(k_he4_ar36__p_k39)* &
      Y(jar36)*state % rho + screened_rates(k_he4_cl33__p_ar36)*Y(jcl33)* &
      state % rho + screened_rates(k_he4_s32__ar36)*Y(js32)*state % rho &
       )
    call set_jac_entry(state, jar36, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_s32__ar36)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jar36, js32, scratch)

    scratch = (&
      screened_rates(k_he4_cl33__p_ar36)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jar36, jcl33, scratch)

    scratch = (&
      screened_rates(k_p_cl35__ar36)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar36, jcl35, scratch)

    scratch = (&
      -screened_rates(k_ar36__he4_s32) - screened_rates(k_ar36__p_cl35) - &
      screened_rates(k_he4_ar36__ca40)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar36__p_k39)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ar36__he4_cl33)*Y(jp)*state % rho - screened_rates(k_p_ar36__k37) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar36, jar36, scratch)

    scratch = (&
      screened_rates(k_k36__ar36__weak__wc12) &
       )
    call set_jac_entry(state, jar36, jk36, scratch)

    scratch = (&
      screened_rates(k_k37__p_ar36) &
       )
    call set_jac_entry(state, jar36, jk37, scratch)

    scratch = (&
      screened_rates(k_p_k39__he4_ar36)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar36, jk39, scratch)

    scratch = (&
      screened_rates(k_ca37__p_ar36__weak__wc12) &
       )
    call set_jac_entry(state, jar36, jca37, scratch)

    scratch = (&
      screened_rates(k_ca40__he4_ar36) &
       )
    call set_jac_entry(state, jar36, jca40, scratch)

    scratch = (&
      screened_rates(k_sc40__he4_ar36__weak__wc12) &
       )
    call set_jac_entry(state, jar36, jsc40, scratch)

    scratch = (&
      -screened_rates(k_p_ar37__he4_cl34)*Y(jar37)*state % rho - screened_rates(k_p_ar37__k38)* &
      Y(jar37)*state % rho + screened_rates(k_p_cl36__ar37)*Y(jcl36)*state % rho + &
      screened_rates(k_p_k40__he4_ar37)*Y(jk40)*state % rho &
       )
    call set_jac_entry(state, jar37, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_ar37__ca41)*Y(jar37)*state % rho - screened_rates(k_he4_ar37__p_k40)* &
      Y(jar37)*state % rho + screened_rates(k_he4_cl34__p_ar37)*Y(jcl34)* &
      state % rho + screened_rates(k_he4_s33__ar37)*Y(js33)*state % rho &
       )
    call set_jac_entry(state, jar37, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_s33__ar37)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jar37, js33, scratch)

    scratch = (&
      screened_rates(k_he4_cl34__p_ar37)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jar37, jcl34, scratch)

    scratch = (&
      screened_rates(k_p_cl36__ar37)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar37, jcl36, scratch)

    scratch = (&
      -screened_rates(k_ar37__cl37__weak__wc12) - screened_rates(k_ar37__he4_s33) - &
      screened_rates(k_ar37__p_cl36) - screened_rates(k_he4_ar37__ca41)*Y(jhe4)* &
      state % rho - screened_rates(k_he4_ar37__p_k40)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ar37__he4_cl34)*Y(jp)*state % rho - screened_rates(k_p_ar37__k38) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar37, jar37, scratch)

    scratch = (&
      screened_rates(k_k37__ar37__weak__wc12) &
       )
    call set_jac_entry(state, jar37, jk37, scratch)

    scratch = (&
      screened_rates(k_k38__p_ar37) &
       )
    call set_jac_entry(state, jar37, jk38, scratch)

    scratch = (&
      screened_rates(k_p_k40__he4_ar37)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar37, jk40, scratch)

    scratch = (&
      screened_rates(k_ca41__he4_ar37) &
       )
    call set_jac_entry(state, jar37, jca41, scratch)

    scratch = (&
      -screened_rates(k_p_ar38__he4_cl35)*Y(jar38)*state % rho - screened_rates(k_p_ar38__k39)* &
      Y(jar38)*state % rho + screened_rates(k_p_cl37__ar38)*Y(jcl37)*state % rho + &
      screened_rates(k_p_k41__he4_ar38)*Y(jk41)*state % rho &
       )
    call set_jac_entry(state, jar38, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_ar38__ca42)*Y(jar38)*state % rho - screened_rates(k_he4_ar38__p_k41)* &
      Y(jar38)*state % rho + screened_rates(k_he4_cl35__p_ar38)*Y(jcl35)* &
      state % rho + screened_rates(k_he4_s34__ar38)*Y(js34)*state % rho &
       )
    call set_jac_entry(state, jar38, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_s34__ar38)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jar38, js34, scratch)

    scratch = (&
      screened_rates(k_he4_cl35__p_ar38)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jar38, jcl35, scratch)

    scratch = (&
      screened_rates(k_p_cl37__ar38)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar38, jcl37, scratch)

    scratch = (&
      -screened_rates(k_ar38__he4_s34) - screened_rates(k_ar38__p_cl37) - &
      screened_rates(k_he4_ar38__ca42)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar38__p_k41)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ar38__he4_cl35)*Y(jp)*state % rho - screened_rates(k_p_ar38__k39) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar38, jar38, scratch)

    scratch = (&
      screened_rates(k_k38__ar38__weak__wc12) &
       )
    call set_jac_entry(state, jar38, jk38, scratch)

    scratch = (&
      screened_rates(k_k39__p_ar38) &
       )
    call set_jac_entry(state, jar38, jk39, scratch)

    scratch = (&
      screened_rates(k_p_k41__he4_ar38)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar38, jk41, scratch)

    scratch = (&
      screened_rates(k_ca42__he4_ar38) &
       )
    call set_jac_entry(state, jar38, jca42, scratch)

    scratch = (&
      -screened_rates(k_p_ar39__he4_cl36)*Y(jar39)*state % rho - screened_rates(k_p_ar39__k40)* &
      Y(jar39)*state % rho &
       )
    call set_jac_entry(state, jar39, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_ar39__ca43)*Y(jar39)*state % rho + screened_rates(k_he4_cl36__p_ar39) &
      *Y(jcl36)*state % rho + screened_rates(k_he4_s35__ar39)*Y(js35)*state % rho &
       )
    call set_jac_entry(state, jar39, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_s35__ar39)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jar39, js35, scratch)

    scratch = (&
      screened_rates(k_he4_cl36__p_ar39)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jar39, jcl36, scratch)

    scratch = (&
      -screened_rates(k_ar39__he4_s35) - screened_rates(k_he4_ar39__ca43)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ar39__he4_cl36)*Y(jp)*state % rho - screened_rates(k_p_ar39__k40) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jar39, jar39, scratch)

    scratch = (&
      screened_rates(k_k40__p_ar39) &
       )
    call set_jac_entry(state, jar39, jk40, scratch)

    scratch = (&
      screened_rates(k_ca43__he4_ar39) &
       )
    call set_jac_entry(state, jar39, jca43, scratch)

    scratch = (&
      screened_rates(k_p_ar32__k33)*Y(jar32)*state % rho + screened_rates(k_p_ca36__he4_k33)* &
      Y(jca36)*state % rho - screened_rates(k_p_k33__ca34)*Y(jk33)*state % rho &
       )
    call set_jac_entry(state, jk33, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cl29__k33)*Y(jcl29)*state % rho - screened_rates(k_he4_k33__p_ca36)* &
      Y(jk33)*state % rho - screened_rates(k_he4_k33__sc37)*Y(jk33)*state % rho &
       )
    call set_jac_entry(state, jk33, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cl29__k33)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jk33, jcl29, scratch)

    scratch = (&
      screened_rates(k_p_ar32__k33)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk33, jar32, scratch)

    scratch = (&
      -screened_rates(k_he4_k33__p_ca36)*Y(jhe4)*state % rho - screened_rates(k_he4_k33__sc37)* &
      Y(jhe4)*state % rho - screened_rates(k_k33__ar33__weak__bqa_pos_) - &
      screened_rates(k_k33__he4_cl29) - screened_rates(k_k33__p_ar32) - &
      screened_rates(k_p_k33__ca34)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk33, jk33, scratch)

    scratch = (&
      screened_rates(k_ca34__p_k33) &
       )
    call set_jac_entry(state, jk33, jca34, scratch)

    scratch = (&
      screened_rates(k_p_ca36__he4_k33)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk33, jca36, scratch)

    scratch = (&
      screened_rates(k_sc37__he4_k33) &
       )
    call set_jac_entry(state, jk33, jsc37, scratch)

    scratch = (&
      screened_rates(k_p_ar33__k34)*Y(jar33)*state % rho + screened_rates(k_p_ca37__he4_k34)* &
      Y(jca37)*state % rho - screened_rates(k_p_k34__ca35)*Y(jk34)*state % rho - &
      screened_rates(k_p_k34__he4_ar31)*Y(jk34)*state % rho &
       )
    call set_jac_entry(state, jk34, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar31__p_k34)*Y(jar31)*state % rho + screened_rates(k_he4_cl30__k34)* &
      Y(jcl30)*state % rho - screened_rates(k_he4_k34__p_ca37)*Y(jk34)*state % rho &
      - screened_rates(k_he4_k34__sc38)*Y(jk34)*state % rho &
       )
    call set_jac_entry(state, jk34, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cl30__k34)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jk34, jcl30, scratch)

    scratch = (&
      screened_rates(k_he4_ar31__p_k34)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jk34, jar31, scratch)

    scratch = (&
      screened_rates(k_p_ar33__k34)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk34, jar33, scratch)

    scratch = (&
      -screened_rates(k_he4_k34__p_ca37)*Y(jhe4)*state % rho - screened_rates(k_he4_k34__sc38)* &
      Y(jhe4)*state % rho - screened_rates(k_k34__ar34__weak__bqa_pos_) - &
      screened_rates(k_k34__he4_cl30) - screened_rates(k_k34__p_ar33) - &
      screened_rates(k_p_k34__ca35)*Y(jp)*state % rho - screened_rates(k_p_k34__he4_ar31)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk34, jk34, scratch)

    scratch = (&
      screened_rates(k_ca34__k34__weak__bqa_pos_) &
       )
    call set_jac_entry(state, jk34, jca34, scratch)

    scratch = (&
      screened_rates(k_ca35__p_k34) &
       )
    call set_jac_entry(state, jk34, jca35, scratch)

    scratch = (&
      screened_rates(k_p_ca37__he4_k34)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk34, jca37, scratch)

    scratch = (&
      screened_rates(k_sc38__he4_k34) &
       )
    call set_jac_entry(state, jk34, jsc38, scratch)

    scratch = (&
      screened_rates(k_p_ar34__k35)*Y(jar34)*state % rho + screened_rates(k_p_ca38__he4_k35)* &
      Y(jca38)*state % rho - screened_rates(k_p_k35__ca36)*Y(jk35)*state % rho - &
      screened_rates(k_p_k35__he4_ar32)*Y(jk35)*state % rho &
       )
    call set_jac_entry(state, jk35, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar32__p_k35)*Y(jar32)*state % rho + screened_rates(k_he4_cl31__k35)* &
      Y(jcl31)*state % rho - screened_rates(k_he4_k35__p_ca38)*Y(jk35)*state % rho &
      - screened_rates(k_he4_k35__sc39)*Y(jk35)*state % rho &
       )
    call set_jac_entry(state, jk35, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cl31__k35)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jk35, jcl31, scratch)

    scratch = (&
      screened_rates(k_he4_ar32__p_k35)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jk35, jar32, scratch)

    scratch = (&
      screened_rates(k_p_ar34__k35)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk35, jar34, scratch)

    scratch = (&
      -screened_rates(k_he4_k35__p_ca38)*Y(jhe4)*state % rho - screened_rates(k_he4_k35__sc39)* &
      Y(jhe4)*state % rho - screened_rates(k_k35__ar35__weak__wc12) - &
      screened_rates(k_k35__he4_cl31) - screened_rates(k_k35__p_ar34) - &
      screened_rates(k_k35__p_cl34__weak__wc12) - screened_rates(k_p_k35__ca36)*Y(jp)* &
      state % rho - screened_rates(k_p_k35__he4_ar32)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk35, jk35, scratch)

    scratch = (&
      screened_rates(k_ca35__k35__weak__wc17) &
       )
    call set_jac_entry(state, jk35, jca35, scratch)

    scratch = (&
      screened_rates(k_ca36__p_k35) &
       )
    call set_jac_entry(state, jk35, jca36, scratch)

    scratch = (&
      screened_rates(k_p_ca38__he4_k35)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk35, jca38, scratch)

    scratch = (&
      screened_rates(k_sc39__he4_k35) &
       )
    call set_jac_entry(state, jk35, jsc39, scratch)

    scratch = (&
      screened_rates(k_p_ar35__k36)*Y(jar35)*state % rho + screened_rates(k_p_ca39__he4_k36)* &
      Y(jca39)*state % rho - screened_rates(k_p_k36__ca37)*Y(jk36)*state % rho - &
      screened_rates(k_p_k36__he4_ar33)*Y(jk36)*state % rho &
       )
    call set_jac_entry(state, jk36, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar33__p_k36)*Y(jar33)*state % rho + screened_rates(k_he4_cl32__k36)* &
      Y(jcl32)*state % rho - screened_rates(k_he4_k36__p_ca39)*Y(jk36)*state % rho &
      - screened_rates(k_he4_k36__sc40)*Y(jk36)*state % rho &
       )
    call set_jac_entry(state, jk36, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cl32__k36)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jk36, jcl32, scratch)

    scratch = (&
      screened_rates(k_he4_ar33__p_k36)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jk36, jar33, scratch)

    scratch = (&
      screened_rates(k_p_ar35__k36)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk36, jar35, scratch)

    scratch = (&
      -screened_rates(k_he4_k36__p_ca39)*Y(jhe4)*state % rho - screened_rates(k_he4_k36__sc40)* &
      Y(jhe4)*state % rho - screened_rates(k_k36__ar36__weak__wc12) - &
      screened_rates(k_k36__he4_cl32) - screened_rates(k_k36__he4_s32__weak__wc12) - &
      screened_rates(k_k36__p_ar35) - screened_rates(k_k36__p_cl35__weak__wc12) - &
      screened_rates(k_p_k36__ca37)*Y(jp)*state % rho - screened_rates(k_p_k36__he4_ar33)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk36, jk36, scratch)

    scratch = (&
      screened_rates(k_ca36__k36__weak__wc12) &
       )
    call set_jac_entry(state, jk36, jca36, scratch)

    scratch = (&
      screened_rates(k_ca37__p_k36) &
       )
    call set_jac_entry(state, jk36, jca37, scratch)

    scratch = (&
      screened_rates(k_p_ca39__he4_k36)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk36, jca39, scratch)

    scratch = (&
      screened_rates(k_sc40__he4_k36) &
       )
    call set_jac_entry(state, jk36, jsc40, scratch)

    scratch = (&
      screened_rates(k_p_ar36__k37)*Y(jar36)*state % rho + screened_rates(k_p_ca40__he4_k37)* &
      Y(jca40)*state % rho - screened_rates(k_p_k37__ca38)*Y(jk37)*state % rho - &
      screened_rates(k_p_k37__he4_ar34)*Y(jk37)*state % rho &
       )
    call set_jac_entry(state, jk37, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar34__p_k37)*Y(jar34)*state % rho + screened_rates(k_he4_cl33__k37)* &
      Y(jcl33)*state % rho - screened_rates(k_he4_k37__p_ca40)*Y(jk37)*state % rho &
      - screened_rates(k_he4_k37__sc41)*Y(jk37)*state % rho &
       )
    call set_jac_entry(state, jk37, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cl33__k37)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jk37, jcl33, scratch)

    scratch = (&
      screened_rates(k_he4_ar34__p_k37)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jk37, jar34, scratch)

    scratch = (&
      screened_rates(k_p_ar36__k37)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk37, jar36, scratch)

    scratch = (&
      -screened_rates(k_he4_k37__p_ca40)*Y(jhe4)*state % rho - screened_rates(k_he4_k37__sc41)* &
      Y(jhe4)*state % rho - screened_rates(k_k37__ar37__weak__wc12) - &
      screened_rates(k_k37__he4_cl33) - screened_rates(k_k37__p_ar36) - &
      screened_rates(k_p_k37__ca38)*Y(jp)*state % rho - screened_rates(k_p_k37__he4_ar34)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk37, jk37, scratch)

    scratch = (&
      screened_rates(k_ca37__k37__weak__wc12) &
       )
    call set_jac_entry(state, jk37, jca37, scratch)

    scratch = (&
      screened_rates(k_ca38__p_k37) &
       )
    call set_jac_entry(state, jk37, jca38, scratch)

    scratch = (&
      screened_rates(k_p_ca40__he4_k37)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk37, jca40, scratch)

    scratch = (&
      screened_rates(k_sc41__he4_k37) &
       )
    call set_jac_entry(state, jk37, jsc41, scratch)

    scratch = (&
      screened_rates(k_p_ar37__k38)*Y(jar37)*state % rho + screened_rates(k_p_ca41__he4_k38)* &
      Y(jca41)*state % rho - screened_rates(k_p_k38__ca39)*Y(jk38)*state % rho - &
      screened_rates(k_p_k38__he4_ar35)*Y(jk38)*state % rho &
       )
    call set_jac_entry(state, jk38, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar35__p_k38)*Y(jar35)*state % rho + screened_rates(k_he4_cl34__k38)* &
      Y(jcl34)*state % rho - screened_rates(k_he4_k38__p_ca41)*Y(jk38)*state % rho &
      - screened_rates(k_he4_k38__sc42)*Y(jk38)*state % rho &
       )
    call set_jac_entry(state, jk38, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cl34__k38)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jk38, jcl34, scratch)

    scratch = (&
      screened_rates(k_he4_ar35__p_k38)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jk38, jar35, scratch)

    scratch = (&
      screened_rates(k_p_ar37__k38)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk38, jar37, scratch)

    scratch = (&
      -screened_rates(k_he4_k38__p_ca41)*Y(jhe4)*state % rho - screened_rates(k_he4_k38__sc42)* &
      Y(jhe4)*state % rho - screened_rates(k_k38__ar38__weak__wc12) - &
      screened_rates(k_k38__he4_cl34) - screened_rates(k_k38__p_ar37) - &
      screened_rates(k_p_k38__ca39)*Y(jp)*state % rho - screened_rates(k_p_k38__he4_ar35)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk38, jk38, scratch)

    scratch = (&
      screened_rates(k_ca38__k38__weak__wc12) &
       )
    call set_jac_entry(state, jk38, jca38, scratch)

    scratch = (&
      screened_rates(k_ca39__p_k38) &
       )
    call set_jac_entry(state, jk38, jca39, scratch)

    scratch = (&
      screened_rates(k_p_ca41__he4_k38)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk38, jca41, scratch)

    scratch = (&
      screened_rates(k_sc42__he4_k38) &
       )
    call set_jac_entry(state, jk38, jsc42, scratch)

    scratch = (&
      screened_rates(k_p_ar38__k39)*Y(jar38)*state % rho + screened_rates(k_p_ca42__he4_k39)* &
      Y(jca42)*state % rho - screened_rates(k_p_k39__ca40)*Y(jk39)*state % rho - &
      screened_rates(k_p_k39__he4_ar36)*Y(jk39)*state % rho &
       )
    call set_jac_entry(state, jk39, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar36__p_k39)*Y(jar36)*state % rho + screened_rates(k_he4_cl35__k39)* &
      Y(jcl35)*state % rho - screened_rates(k_he4_k39__p_ca42)*Y(jk39)*state % rho &
      - screened_rates(k_he4_k39__sc43)*Y(jk39)*state % rho &
       )
    call set_jac_entry(state, jk39, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cl35__k39)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jk39, jcl35, scratch)

    scratch = (&
      screened_rates(k_he4_ar36__p_k39)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jk39, jar36, scratch)

    scratch = (&
      screened_rates(k_p_ar38__k39)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk39, jar38, scratch)

    scratch = (&
      -screened_rates(k_he4_k39__p_ca42)*Y(jhe4)*state % rho - screened_rates(k_he4_k39__sc43)* &
      Y(jhe4)*state % rho - screened_rates(k_k39__he4_cl35) - &
      screened_rates(k_k39__p_ar38) - screened_rates(k_p_k39__ca40)*Y(jp)*state % rho - &
      screened_rates(k_p_k39__he4_ar36)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk39, jk39, scratch)

    scratch = (&
      screened_rates(k_ca39__k39__weak__wc12) &
       )
    call set_jac_entry(state, jk39, jca39, scratch)

    scratch = (&
      screened_rates(k_ca40__p_k39) &
       )
    call set_jac_entry(state, jk39, jca40, scratch)

    scratch = (&
      screened_rates(k_p_ca42__he4_k39)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk39, jca42, scratch)

    scratch = (&
      screened_rates(k_sc40__p_k39__weak__wc12) &
       )
    call set_jac_entry(state, jk39, jsc40, scratch)

    scratch = (&
      screened_rates(k_sc43__he4_k39) &
       )
    call set_jac_entry(state, jk39, jsc43, scratch)

    scratch = (&
      screened_rates(k_p_ar39__k40)*Y(jar39)*state % rho + screened_rates(k_p_ca43__he4_k40)* &
      Y(jca43)*state % rho - screened_rates(k_p_k40__ca41)*Y(jk40)*state % rho - &
      screened_rates(k_p_k40__he4_ar37)*Y(jk40)*state % rho &
       )
    call set_jac_entry(state, jk40, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar37__p_k40)*Y(jar37)*state % rho + screened_rates(k_he4_cl36__k40)* &
      Y(jcl36)*state % rho - screened_rates(k_he4_k40__p_ca43)*Y(jk40)*state % rho &
      - screened_rates(k_he4_k40__sc44)*Y(jk40)*state % rho &
       )
    call set_jac_entry(state, jk40, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cl36__k40)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jk40, jcl36, scratch)

    scratch = (&
      screened_rates(k_he4_ar37__p_k40)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jk40, jar37, scratch)

    scratch = (&
      screened_rates(k_p_ar39__k40)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk40, jar39, scratch)

    scratch = (&
      -screened_rates(k_he4_k40__p_ca43)*Y(jhe4)*state % rho - screened_rates(k_he4_k40__sc44)* &
      Y(jhe4)*state % rho - screened_rates(k_k40__he4_cl36) - &
      screened_rates(k_k40__p_ar39) - screened_rates(k_p_k40__ca41)*Y(jp)*state % rho - &
      screened_rates(k_p_k40__he4_ar37)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk40, jk40, scratch)

    scratch = (&
      screened_rates(k_ca41__p_k40) &
       )
    call set_jac_entry(state, jk40, jca41, scratch)

    scratch = (&
      screened_rates(k_p_ca43__he4_k40)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk40, jca43, scratch)

    scratch = (&
      screened_rates(k_sc44__he4_k40) &
       )
    call set_jac_entry(state, jk40, jsc44, scratch)

    scratch = (&
      screened_rates(k_p_ca44__he4_k41)*Y(jca44)*state % rho - screened_rates(k_p_k41__ca42)* &
      Y(jk41)*state % rho - screened_rates(k_p_k41__he4_ar38)*Y(jk41)*state % rho &
       )
    call set_jac_entry(state, jk41, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar38__p_k41)*Y(jar38)*state % rho + screened_rates(k_he4_cl37__k41)* &
      Y(jcl37)*state % rho - screened_rates(k_he4_k41__p_ca44)*Y(jk41)*state % rho &
      - screened_rates(k_he4_k41__sc45)*Y(jk41)*state % rho &
       )
    call set_jac_entry(state, jk41, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cl37__k41)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jk41, jcl37, scratch)

    scratch = (&
      screened_rates(k_he4_ar38__p_k41)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jk41, jar38, scratch)

    scratch = (&
      -screened_rates(k_he4_k41__p_ca44)*Y(jhe4)*state % rho - screened_rates(k_he4_k41__sc45)* &
      Y(jhe4)*state % rho - screened_rates(k_k41__he4_cl37) - &
      screened_rates(k_p_k41__ca42)*Y(jp)*state % rho - screened_rates(k_p_k41__he4_ar38)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk41, jk41, scratch)

    scratch = (&
      screened_rates(k_ca41__k41__weak__wc12) &
       )
    call set_jac_entry(state, jk41, jca41, scratch)

    scratch = (&
      screened_rates(k_ca42__p_k41) &
       )
    call set_jac_entry(state, jk41, jca42, scratch)

    scratch = (&
      screened_rates(k_p_ca44__he4_k41)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jk41, jca44, scratch)

    scratch = (&
      screened_rates(k_sc45__he4_k41) &
       )
    call set_jac_entry(state, jk41, jsc45, scratch)

    scratch = (&
      screened_rates(k_p_k33__ca34)*Y(jk33)*state % rho + screened_rates(k_p_sc37__he4_ca34)* &
      Y(jsc37)*state % rho &
       )
    call set_jac_entry(state, jca34, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_ca34__p_sc37)*Y(jca34)*state % rho - screened_rates(k_he4_ca34__ti38) &
      *Y(jca34)*state % rho &
       )
    call set_jac_entry(state, jca34, jhe4, scratch)

    scratch = (&
      screened_rates(k_p_k33__ca34)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca34, jk33, scratch)

    scratch = (&
      -screened_rates(k_ca34__k34__weak__bqa_pos_) - screened_rates(k_ca34__p_k33) - &
      screened_rates(k_he4_ca34__p_sc37)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca34__ti38)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca34, jca34, scratch)

    scratch = (&
      screened_rates(k_p_sc37__he4_ca34)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca34, jsc37, scratch)

    scratch = (&
      screened_rates(k_ti38__he4_ca34) &
       )
    call set_jac_entry(state, jca34, jti38, scratch)

    scratch = (&
      -screened_rates(k_p_ca35__sc36)*Y(jca35)*state % rho + screened_rates(k_p_k34__ca35)* &
      Y(jk34)*state % rho + screened_rates(k_p_sc38__he4_ca35)*Y(jsc38)* &
      state % rho &
       )
    call set_jac_entry(state, jca35, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar31__ca35)*Y(jar31)*state % rho - screened_rates(k_he4_ca35__p_sc38)* &
      Y(jca35)*state % rho - screened_rates(k_he4_ca35__ti39)*Y(jca35)*state % rho &
       )
    call set_jac_entry(state, jca35, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ar31__ca35)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca35, jar31, scratch)

    scratch = (&
      screened_rates(k_p_k34__ca35)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca35, jk34, scratch)

    scratch = (&
      -screened_rates(k_ca35__he4_ar31) - screened_rates(k_ca35__k35__weak__wc17) - &
      screened_rates(k_ca35__p_ar34__weak__wc17) - screened_rates(k_ca35__p_k34) - &
      screened_rates(k_ca35__p_p_cl33__weak__wc12) - screened_rates(k_he4_ca35__p_sc38)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_ca35__ti39)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ca35__sc36)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca35, jca35, scratch)

    scratch = (&
      screened_rates(k_sc36__p_ca35) &
       )
    call set_jac_entry(state, jca35, jsc36, scratch)

    scratch = (&
      screened_rates(k_p_sc38__he4_ca35)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca35, jsc38, scratch)

    scratch = (&
      screened_rates(k_ti39__he4_ca35) &
       )
    call set_jac_entry(state, jca35, jti39, scratch)

    scratch = (&
      -screened_rates(k_p_ca36__he4_k33)*Y(jca36)*state % rho - screened_rates(k_p_ca36__sc37)* &
      Y(jca36)*state % rho + screened_rates(k_p_k35__ca36)*Y(jk35)*state % rho + &
      1.0d0*screened_rates(k_p_p_ca38__he4_ca36)*Y(jca38)*Y(jp)*state % rho**2 &
      + screened_rates(k_p_sc39__he4_ca36)*Y(jsc39)*state % rho &
       )
    call set_jac_entry(state, jca36, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar32__ca36)*Y(jar32)*state % rho - &
      screened_rates(k_he4_ca36__p_p_ca38)*Y(jca36)*state % rho - &
      screened_rates(k_he4_ca36__p_sc39)*Y(jca36)*state % rho - &
      screened_rates(k_he4_ca36__ti40)*Y(jca36)*state % rho + &
      screened_rates(k_he4_k33__p_ca36)*Y(jk33)*state % rho &
       )
    call set_jac_entry(state, jca36, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ar32__ca36)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca36, jar32, scratch)

    scratch = (&
      screened_rates(k_he4_k33__p_ca36)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca36, jk33, scratch)

    scratch = (&
      screened_rates(k_p_k35__ca36)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca36, jk35, scratch)

    scratch = (&
      -screened_rates(k_ca36__he4_ar32) - screened_rates(k_ca36__k36__weak__wc12) - &
      screened_rates(k_ca36__p_ar35__weak__wc12) - screened_rates(k_ca36__p_k35) - &
      screened_rates(k_he4_ca36__p_p_ca38)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca36__p_sc39)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca36__ti40)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ca36__he4_k33)*Y(jp)*state % rho - screened_rates(k_p_ca36__sc37) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca36, jca36, scratch)

    scratch = (&
      0.5d0*screened_rates(k_p_p_ca38__he4_ca36)*Y(jp)**2*state % rho**2 &
       )
    call set_jac_entry(state, jca36, jca38, scratch)

    scratch = (&
      screened_rates(k_sc36__ca36__weak__bqa_pos_) &
       )
    call set_jac_entry(state, jca36, jsc36, scratch)

    scratch = (&
      screened_rates(k_sc37__p_ca36) &
       )
    call set_jac_entry(state, jca36, jsc37, scratch)

    scratch = (&
      screened_rates(k_p_sc39__he4_ca36)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca36, jsc39, scratch)

    scratch = (&
      screened_rates(k_ti40__he4_ca36) &
       )
    call set_jac_entry(state, jca36, jti40, scratch)

    scratch = (&
      -screened_rates(k_p_ca37__he4_k34)*Y(jca37)*state % rho - screened_rates(k_p_ca37__sc38)* &
      Y(jca37)*state % rho + screened_rates(k_p_k36__ca37)*Y(jk36)*state % rho + &
      screened_rates(k_p_sc40__he4_ca37)*Y(jsc40)*state % rho &
       )
    call set_jac_entry(state, jca37, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar33__ca37)*Y(jar33)*state % rho - screened_rates(k_he4_ca37__p_sc40)* &
      Y(jca37)*state % rho - screened_rates(k_he4_ca37__ti41)*Y(jca37)*state % rho &
      + screened_rates(k_he4_k34__p_ca37)*Y(jk34)*state % rho &
       )
    call set_jac_entry(state, jca37, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ar33__ca37)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca37, jar33, scratch)

    scratch = (&
      screened_rates(k_he4_k34__p_ca37)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca37, jk34, scratch)

    scratch = (&
      screened_rates(k_p_k36__ca37)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca37, jk36, scratch)

    scratch = (&
      -screened_rates(k_ca37__he4_ar33) - screened_rates(k_ca37__k37__weak__wc12) - &
      screened_rates(k_ca37__p_ar36__weak__wc12) - screened_rates(k_ca37__p_k36) - &
      screened_rates(k_he4_ca37__p_sc40)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca37__ti41)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ca37__he4_k34)*Y(jp)*state % rho - screened_rates(k_p_ca37__sc38) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca37, jca37, scratch)

    scratch = (&
      screened_rates(k_sc37__ca37__weak__bqa_pos_) &
       )
    call set_jac_entry(state, jca37, jsc37, scratch)

    scratch = (&
      screened_rates(k_sc38__p_ca37) &
       )
    call set_jac_entry(state, jca37, jsc38, scratch)

    scratch = (&
      screened_rates(k_p_sc40__he4_ca37)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca37, jsc40, scratch)

    scratch = (&
      screened_rates(k_ti41__he4_ca37) &
       )
    call set_jac_entry(state, jca37, jti41, scratch)

    scratch = (&
      -screened_rates(k_p_ca38__he4_k35)*Y(jca38)*state % rho - screened_rates(k_p_ca38__sc39)* &
      Y(jca38)*state % rho + screened_rates(k_p_k37__ca38)*Y(jk37)*state % rho - &
      1.0d0*screened_rates(k_p_p_ca38__he4_ca36)*Y(jca38)*Y(jp)*state % rho**2 &
      + screened_rates(k_p_sc41__he4_ca38)*Y(jsc41)*state % rho &
       )
    call set_jac_entry(state, jca38, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar34__ca38)*Y(jar34)*state % rho + &
      screened_rates(k_he4_ca36__p_p_ca38)*Y(jca36)*state % rho - &
      screened_rates(k_he4_ca38__p_sc41)*Y(jca38)*state % rho - &
      screened_rates(k_he4_ca38__ti42)*Y(jca38)*state % rho + &
      screened_rates(k_he4_k35__p_ca38)*Y(jk35)*state % rho &
       )
    call set_jac_entry(state, jca38, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ar34__ca38)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca38, jar34, scratch)

    scratch = (&
      screened_rates(k_he4_k35__p_ca38)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca38, jk35, scratch)

    scratch = (&
      screened_rates(k_p_k37__ca38)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca38, jk37, scratch)

    scratch = (&
      screened_rates(k_he4_ca36__p_p_ca38)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca38, jca36, scratch)

    scratch = (&
      -screened_rates(k_ca38__he4_ar34) - screened_rates(k_ca38__k38__weak__wc12) - &
      screened_rates(k_ca38__p_k37) - screened_rates(k_he4_ca38__p_sc41)*Y(jhe4)* &
      state % rho - screened_rates(k_he4_ca38__ti42)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ca38__he4_k35)*Y(jp)*state % rho - screened_rates(k_p_ca38__sc39) &
      *Y(jp)*state % rho - 0.5d0*screened_rates(k_p_p_ca38__he4_ca36)*Y(jp)**2* &
      state % rho**2 &
       )
    call set_jac_entry(state, jca38, jca38, scratch)

    scratch = (&
      screened_rates(k_sc38__ca38__weak__mo97) &
       )
    call set_jac_entry(state, jca38, jsc38, scratch)

    scratch = (&
      screened_rates(k_sc39__p_ca38) &
       )
    call set_jac_entry(state, jca38, jsc39, scratch)

    scratch = (&
      screened_rates(k_p_sc41__he4_ca38)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca38, jsc41, scratch)

    scratch = (&
      screened_rates(k_ti39__p_ca38__weak__wc12) &
       )
    call set_jac_entry(state, jca38, jti39, scratch)

    scratch = (&
      screened_rates(k_ti42__he4_ca38) &
       )
    call set_jac_entry(state, jca38, jti42, scratch)

    scratch = (&
      -screened_rates(k_p_ca39__he4_k36)*Y(jca39)*state % rho - screened_rates(k_p_ca39__sc40)* &
      Y(jca39)*state % rho + screened_rates(k_p_k38__ca39)*Y(jk38)*state % rho + &
      screened_rates(k_p_sc42__he4_ca39)*Y(jsc42)*state % rho &
       )
    call set_jac_entry(state, jca39, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar35__ca39)*Y(jar35)*state % rho - screened_rates(k_he4_ca39__p_sc42)* &
      Y(jca39)*state % rho - screened_rates(k_he4_ca39__ti43)*Y(jca39)*state % rho &
      + screened_rates(k_he4_k36__p_ca39)*Y(jk36)*state % rho &
       )
    call set_jac_entry(state, jca39, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ar35__ca39)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca39, jar35, scratch)

    scratch = (&
      screened_rates(k_he4_k36__p_ca39)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca39, jk36, scratch)

    scratch = (&
      screened_rates(k_p_k38__ca39)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca39, jk38, scratch)

    scratch = (&
      -screened_rates(k_ca39__he4_ar35) - screened_rates(k_ca39__k39__weak__wc12) - &
      screened_rates(k_ca39__p_k38) - screened_rates(k_he4_ca39__p_sc42)*Y(jhe4)* &
      state % rho - screened_rates(k_he4_ca39__ti43)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ca39__he4_k36)*Y(jp)*state % rho - screened_rates(k_p_ca39__sc40) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca39, jca39, scratch)

    scratch = (&
      screened_rates(k_sc39__ca39__weak__mo97) &
       )
    call set_jac_entry(state, jca39, jsc39, scratch)

    scratch = (&
      screened_rates(k_sc40__p_ca39) &
       )
    call set_jac_entry(state, jca39, jsc40, scratch)

    scratch = (&
      screened_rates(k_p_sc42__he4_ca39)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca39, jsc42, scratch)

    scratch = (&
      screened_rates(k_ti40__p_ca39__weak__wc12) &
       )
    call set_jac_entry(state, jca39, jti40, scratch)

    scratch = (&
      screened_rates(k_ti43__he4_ca39) &
       )
    call set_jac_entry(state, jca39, jti43, scratch)

    scratch = (&
      -screened_rates(k_p_ca40__he4_k37)*Y(jca40)*state % rho - screened_rates(k_p_ca40__sc41)* &
      Y(jca40)*state % rho + screened_rates(k_p_k39__ca40)*Y(jk39)*state % rho + &
      screened_rates(k_p_sc43__he4_ca40)*Y(jsc43)*state % rho &
       )
    call set_jac_entry(state, jca40, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar36__ca40)*Y(jar36)*state % rho - screened_rates(k_he4_ca40__p_sc43)* &
      Y(jca40)*state % rho - screened_rates(k_he4_ca40__ti44)*Y(jca40)*state % rho &
      + screened_rates(k_he4_k37__p_ca40)*Y(jk37)*state % rho &
       )
    call set_jac_entry(state, jca40, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ar36__ca40)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca40, jar36, scratch)

    scratch = (&
      screened_rates(k_he4_k37__p_ca40)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca40, jk37, scratch)

    scratch = (&
      screened_rates(k_p_k39__ca40)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca40, jk39, scratch)

    scratch = (&
      -screened_rates(k_ca40__he4_ar36) - screened_rates(k_ca40__p_k39) - &
      screened_rates(k_he4_ca40__p_sc43)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca40__ti44)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ca40__he4_k37)*Y(jp)*state % rho - screened_rates(k_p_ca40__sc41) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca40, jca40, scratch)

    scratch = (&
      screened_rates(k_sc40__ca40__weak__wc12) &
       )
    call set_jac_entry(state, jca40, jsc40, scratch)

    scratch = (&
      screened_rates(k_sc41__p_ca40) &
       )
    call set_jac_entry(state, jca40, jsc41, scratch)

    scratch = (&
      screened_rates(k_p_sc43__he4_ca40)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca40, jsc43, scratch)

    scratch = (&
      screened_rates(k_ti41__p_ca40__weak__wc12) &
       )
    call set_jac_entry(state, jca40, jti41, scratch)

    scratch = (&
      screened_rates(k_ti44__he4_ca40) &
       )
    call set_jac_entry(state, jca40, jti44, scratch)

    scratch = (&
      screened_rates(k_cr43__p_p_p_ca40__weak__wc12) &
       )
    call set_jac_entry(state, jca40, jcr43, scratch)

    scratch = (&
      -screened_rates(k_p_ca41__he4_k38)*Y(jca41)*state % rho - screened_rates(k_p_ca41__sc42)* &
      Y(jca41)*state % rho + screened_rates(k_p_k40__ca41)*Y(jk40)*state % rho + &
      screened_rates(k_p_sc44__he4_ca41)*Y(jsc44)*state % rho &
       )
    call set_jac_entry(state, jca41, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar37__ca41)*Y(jar37)*state % rho - screened_rates(k_he4_ca41__p_sc44)* &
      Y(jca41)*state % rho - screened_rates(k_he4_ca41__ti45)*Y(jca41)*state % rho &
      + screened_rates(k_he4_k38__p_ca41)*Y(jk38)*state % rho &
       )
    call set_jac_entry(state, jca41, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ar37__ca41)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca41, jar37, scratch)

    scratch = (&
      screened_rates(k_he4_k38__p_ca41)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca41, jk38, scratch)

    scratch = (&
      screened_rates(k_p_k40__ca41)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca41, jk40, scratch)

    scratch = (&
      -screened_rates(k_ca41__he4_ar37) - screened_rates(k_ca41__k41__weak__wc12) - &
      screened_rates(k_ca41__p_k40) - screened_rates(k_he4_ca41__p_sc44)*Y(jhe4)* &
      state % rho - screened_rates(k_he4_ca41__ti45)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ca41__he4_k38)*Y(jp)*state % rho - screened_rates(k_p_ca41__sc42) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca41, jca41, scratch)

    scratch = (&
      screened_rates(k_sc41__ca41__weak__wc12) &
       )
    call set_jac_entry(state, jca41, jsc41, scratch)

    scratch = (&
      screened_rates(k_sc42__p_ca41) &
       )
    call set_jac_entry(state, jca41, jsc42, scratch)

    scratch = (&
      screened_rates(k_p_sc44__he4_ca41)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca41, jsc44, scratch)

    scratch = (&
      screened_rates(k_ti45__he4_ca41) &
       )
    call set_jac_entry(state, jca41, jti45, scratch)

    scratch = (&
      -screened_rates(k_p_ca42__he4_k39)*Y(jca42)*state % rho - screened_rates(k_p_ca42__sc43)* &
      Y(jca42)*state % rho + screened_rates(k_p_k41__ca42)*Y(jk41)*state % rho + &
      screened_rates(k_p_sc45__he4_ca42)*Y(jsc45)*state % rho &
       )
    call set_jac_entry(state, jca42, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar38__ca42)*Y(jar38)*state % rho - screened_rates(k_he4_ca42__p_sc45)* &
      Y(jca42)*state % rho - screened_rates(k_he4_ca42__ti46)*Y(jca42)*state % rho &
      + screened_rates(k_he4_k39__p_ca42)*Y(jk39)*state % rho &
       )
    call set_jac_entry(state, jca42, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ar38__ca42)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca42, jar38, scratch)

    scratch = (&
      screened_rates(k_he4_k39__p_ca42)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca42, jk39, scratch)

    scratch = (&
      screened_rates(k_p_k41__ca42)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca42, jk41, scratch)

    scratch = (&
      -screened_rates(k_ca42__he4_ar38) - screened_rates(k_ca42__p_k41) - &
      screened_rates(k_he4_ca42__p_sc45)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca42__ti46)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ca42__he4_k39)*Y(jp)*state % rho - screened_rates(k_p_ca42__sc43) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca42, jca42, scratch)

    scratch = (&
      screened_rates(k_sc42__ca42__weak__wc12) &
       )
    call set_jac_entry(state, jca42, jsc42, scratch)

    scratch = (&
      screened_rates(k_sc43__p_ca42) &
       )
    call set_jac_entry(state, jca42, jsc43, scratch)

    scratch = (&
      screened_rates(k_p_sc45__he4_ca42)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca42, jsc45, scratch)

    scratch = (&
      screened_rates(k_ti46__he4_ca42) &
       )
    call set_jac_entry(state, jca42, jti46, scratch)

    scratch = (&
      -screened_rates(k_p_ca43__he4_k40)*Y(jca43)*state % rho - screened_rates(k_p_ca43__sc44)* &
      Y(jca43)*state % rho + screened_rates(k_p_sc46__he4_ca43)*Y(jsc46)* &
      state % rho &
       )
    call set_jac_entry(state, jca43, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar39__ca43)*Y(jar39)*state % rho - screened_rates(k_he4_ca43__p_sc46)* &
      Y(jca43)*state % rho - screened_rates(k_he4_ca43__ti47)*Y(jca43)*state % rho &
      + screened_rates(k_he4_k40__p_ca43)*Y(jk40)*state % rho &
       )
    call set_jac_entry(state, jca43, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ar39__ca43)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca43, jar39, scratch)

    scratch = (&
      screened_rates(k_he4_k40__p_ca43)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca43, jk40, scratch)

    scratch = (&
      -screened_rates(k_ca43__he4_ar39) - screened_rates(k_he4_ca43__p_sc46)*Y(jhe4)*state % rho &
      - screened_rates(k_he4_ca43__ti47)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ca43__he4_k40)*Y(jp)*state % rho - screened_rates(k_p_ca43__sc44) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca43, jca43, scratch)

    scratch = (&
      screened_rates(k_sc43__ca43__weak__wc12) &
       )
    call set_jac_entry(state, jca43, jsc43, scratch)

    scratch = (&
      screened_rates(k_sc44__p_ca43) &
       )
    call set_jac_entry(state, jca43, jsc44, scratch)

    scratch = (&
      screened_rates(k_p_sc46__he4_ca43)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca43, jsc46, scratch)

    scratch = (&
      screened_rates(k_ti47__he4_ca43) &
       )
    call set_jac_entry(state, jca43, jti47, scratch)

    scratch = (&
      -screened_rates(k_p_ca44__he4_k41)*Y(jca44)*state % rho - screened_rates(k_p_ca44__sc45)* &
      Y(jca44)*state % rho &
       )
    call set_jac_entry(state, jca44, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_ca44__ti48)*Y(jca44)*state % rho + screened_rates(k_he4_k41__p_ca44)* &
      Y(jk41)*state % rho &
       )
    call set_jac_entry(state, jca44, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_k41__p_ca44)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jca44, jk41, scratch)

    scratch = (&
      -screened_rates(k_he4_ca44__ti48)*Y(jhe4)*state % rho - screened_rates(k_p_ca44__he4_k41)* &
      Y(jp)*state % rho - screened_rates(k_p_ca44__sc45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jca44, jca44, scratch)

    scratch = (&
      screened_rates(k_sc44__ca44__weak__wc12) &
       )
    call set_jac_entry(state, jca44, jsc44, scratch)

    scratch = (&
      screened_rates(k_sc45__p_ca44) &
       )
    call set_jac_entry(state, jca44, jsc45, scratch)

    scratch = (&
      screened_rates(k_ti48__he4_ca44) &
       )
    call set_jac_entry(state, jca44, jti48, scratch)

    scratch = (&
      screened_rates(k_p_ca35__sc36)*Y(jca35)*state % rho + screened_rates(k_p_ti39__he4_sc36)* &
      Y(jti39)*state % rho &
       )
    call set_jac_entry(state, jsc36, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_sc36__p_ti39)*Y(jsc36)*state % rho - screened_rates(k_he4_sc36__v40)* &
      Y(jsc36)*state % rho &
       )
    call set_jac_entry(state, jsc36, jhe4, scratch)

    scratch = (&
      screened_rates(k_p_ca35__sc36)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc36, jca35, scratch)

    scratch = (&
      -screened_rates(k_he4_sc36__p_ti39)*Y(jhe4)*state % rho - screened_rates(k_he4_sc36__v40)* &
      Y(jhe4)*state % rho - screened_rates(k_sc36__ca36__weak__bqa_pos_) - &
      screened_rates(k_sc36__p_ca35) &
       )
    call set_jac_entry(state, jsc36, jsc36, scratch)

    scratch = (&
      screened_rates(k_p_ti39__he4_sc36)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc36, jti39, scratch)

    scratch = (&
      screened_rates(k_v40__he4_sc36) &
       )
    call set_jac_entry(state, jsc36, jv40, scratch)

    scratch = (&
      screened_rates(k_p_ca36__sc37)*Y(jca36)*state % rho - screened_rates(k_p_sc37__he4_ca34)* &
      Y(jsc37)*state % rho - screened_rates(k_p_sc37__ti38)*Y(jsc37)*state % rho + &
      screened_rates(k_p_ti40__he4_sc37)*Y(jti40)*state % rho &
       )
    call set_jac_entry(state, jsc37, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca34__p_sc37)*Y(jca34)*state % rho + screened_rates(k_he4_k33__sc37)* &
      Y(jk33)*state % rho - screened_rates(k_he4_sc37__p_ti40)*Y(jsc37)* &
      state % rho - screened_rates(k_he4_sc37__v41)*Y(jsc37)*state % rho &
       )
    call set_jac_entry(state, jsc37, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_k33__sc37)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc37, jk33, scratch)

    scratch = (&
      screened_rates(k_he4_ca34__p_sc37)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc37, jca34, scratch)

    scratch = (&
      screened_rates(k_p_ca36__sc37)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc37, jca36, scratch)

    scratch = (&
      -screened_rates(k_he4_sc37__p_ti40)*Y(jhe4)*state % rho - screened_rates(k_he4_sc37__v41)* &
      Y(jhe4)*state % rho - screened_rates(k_p_sc37__he4_ca34)*Y(jp)*state % rho - &
      screened_rates(k_p_sc37__ti38)*Y(jp)*state % rho - &
      screened_rates(k_sc37__ca37__weak__bqa_pos_) - screened_rates(k_sc37__he4_k33) - &
      screened_rates(k_sc37__p_ca36) &
       )
    call set_jac_entry(state, jsc37, jsc37, scratch)

    scratch = (&
      screened_rates(k_ti38__p_sc37) &
       )
    call set_jac_entry(state, jsc37, jti38, scratch)

    scratch = (&
      screened_rates(k_p_ti40__he4_sc37)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc37, jti40, scratch)

    scratch = (&
      screened_rates(k_v41__he4_sc37) &
       )
    call set_jac_entry(state, jsc37, jv41, scratch)

    scratch = (&
      screened_rates(k_p_ca37__sc38)*Y(jca37)*state % rho - screened_rates(k_p_sc38__he4_ca35)* &
      Y(jsc38)*state % rho - screened_rates(k_p_sc38__ti39)*Y(jsc38)*state % rho + &
      screened_rates(k_p_ti41__he4_sc38)*Y(jti41)*state % rho &
       )
    call set_jac_entry(state, jsc38, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca35__p_sc38)*Y(jca35)*state % rho + screened_rates(k_he4_k34__sc38)* &
      Y(jk34)*state % rho - screened_rates(k_he4_sc38__p_ti41)*Y(jsc38)* &
      state % rho - screened_rates(k_he4_sc38__v42)*Y(jsc38)*state % rho &
       )
    call set_jac_entry(state, jsc38, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_k34__sc38)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc38, jk34, scratch)

    scratch = (&
      screened_rates(k_he4_ca35__p_sc38)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc38, jca35, scratch)

    scratch = (&
      screened_rates(k_p_ca37__sc38)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc38, jca37, scratch)

    scratch = (&
      -screened_rates(k_he4_sc38__p_ti41)*Y(jhe4)*state % rho - screened_rates(k_he4_sc38__v42)* &
      Y(jhe4)*state % rho - screened_rates(k_p_sc38__he4_ca35)*Y(jp)*state % rho - &
      screened_rates(k_p_sc38__ti39)*Y(jp)*state % rho - &
      screened_rates(k_sc38__ca38__weak__mo97) - screened_rates(k_sc38__he4_k34) - &
      screened_rates(k_sc38__p_ca37) &
       )
    call set_jac_entry(state, jsc38, jsc38, scratch)

    scratch = (&
      screened_rates(k_ti38__sc38__weak__mo97) &
       )
    call set_jac_entry(state, jsc38, jti38, scratch)

    scratch = (&
      screened_rates(k_ti39__p_sc38) &
       )
    call set_jac_entry(state, jsc38, jti39, scratch)

    scratch = (&
      screened_rates(k_p_ti41__he4_sc38)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc38, jti41, scratch)

    scratch = (&
      screened_rates(k_v42__he4_sc38) &
       )
    call set_jac_entry(state, jsc38, jv42, scratch)

    scratch = (&
      screened_rates(k_p_ca38__sc39)*Y(jca38)*state % rho - screened_rates(k_p_sc39__he4_ca36)* &
      Y(jsc39)*state % rho - screened_rates(k_p_sc39__ti40)*Y(jsc39)*state % rho + &
      screened_rates(k_p_ti42__he4_sc39)*Y(jti42)*state % rho &
       )
    call set_jac_entry(state, jsc39, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca36__p_sc39)*Y(jca36)*state % rho + screened_rates(k_he4_k35__sc39)* &
      Y(jk35)*state % rho - screened_rates(k_he4_sc39__p_ti42)*Y(jsc39)* &
      state % rho - screened_rates(k_he4_sc39__v43)*Y(jsc39)*state % rho &
       )
    call set_jac_entry(state, jsc39, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_k35__sc39)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc39, jk35, scratch)

    scratch = (&
      screened_rates(k_he4_ca36__p_sc39)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc39, jca36, scratch)

    scratch = (&
      screened_rates(k_p_ca38__sc39)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc39, jca38, scratch)

    scratch = (&
      -screened_rates(k_he4_sc39__p_ti42)*Y(jhe4)*state % rho - screened_rates(k_he4_sc39__v43)* &
      Y(jhe4)*state % rho - screened_rates(k_p_sc39__he4_ca36)*Y(jp)*state % rho - &
      screened_rates(k_p_sc39__ti40)*Y(jp)*state % rho - &
      screened_rates(k_sc39__ca39__weak__mo97) - screened_rates(k_sc39__he4_k35) - &
      screened_rates(k_sc39__p_ca38) &
       )
    call set_jac_entry(state, jsc39, jsc39, scratch)

    scratch = (&
      screened_rates(k_ti39__sc39__weak__wc17) &
       )
    call set_jac_entry(state, jsc39, jti39, scratch)

    scratch = (&
      screened_rates(k_ti40__p_sc39) &
       )
    call set_jac_entry(state, jsc39, jti40, scratch)

    scratch = (&
      screened_rates(k_p_ti42__he4_sc39)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc39, jti42, scratch)

    scratch = (&
      screened_rates(k_v43__he4_sc39) &
       )
    call set_jac_entry(state, jsc39, jv43, scratch)

    scratch = (&
      screened_rates(k_p_ca39__sc40)*Y(jca39)*state % rho - screened_rates(k_p_sc40__he4_ca37)* &
      Y(jsc40)*state % rho - screened_rates(k_p_sc40__ti41)*Y(jsc40)*state % rho + &
      screened_rates(k_p_ti43__he4_sc40)*Y(jti43)*state % rho &
       )
    call set_jac_entry(state, jsc40, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca37__p_sc40)*Y(jca37)*state % rho + screened_rates(k_he4_k36__sc40)* &
      Y(jk36)*state % rho - screened_rates(k_he4_sc40__p_ti43)*Y(jsc40)* &
      state % rho - screened_rates(k_he4_sc40__v44)*Y(jsc40)*state % rho &
       )
    call set_jac_entry(state, jsc40, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_k36__sc40)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc40, jk36, scratch)

    scratch = (&
      screened_rates(k_he4_ca37__p_sc40)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc40, jca37, scratch)

    scratch = (&
      screened_rates(k_p_ca39__sc40)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc40, jca39, scratch)

    scratch = (&
      -screened_rates(k_he4_sc40__p_ti43)*Y(jhe4)*state % rho - screened_rates(k_he4_sc40__v44)* &
      Y(jhe4)*state % rho - screened_rates(k_p_sc40__he4_ca37)*Y(jp)*state % rho - &
      screened_rates(k_p_sc40__ti41)*Y(jp)*state % rho - &
      screened_rates(k_sc40__ca40__weak__wc12) - &
      screened_rates(k_sc40__he4_ar36__weak__wc12) - screened_rates(k_sc40__he4_k36) - &
      screened_rates(k_sc40__p_ca39) - screened_rates(k_sc40__p_k39__weak__wc12) &
       )
    call set_jac_entry(state, jsc40, jsc40, scratch)

    scratch = (&
      screened_rates(k_ti40__sc40__weak__wc17) &
       )
    call set_jac_entry(state, jsc40, jti40, scratch)

    scratch = (&
      screened_rates(k_ti41__p_sc40) &
       )
    call set_jac_entry(state, jsc40, jti41, scratch)

    scratch = (&
      screened_rates(k_p_ti43__he4_sc40)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc40, jti43, scratch)

    scratch = (&
      screened_rates(k_v44__he4_sc40) &
       )
    call set_jac_entry(state, jsc40, jv44, scratch)

    scratch = (&
      screened_rates(k_p_ca40__sc41)*Y(jca40)*state % rho - screened_rates(k_p_sc41__he4_ca38)* &
      Y(jsc41)*state % rho - screened_rates(k_p_sc41__ti42)*Y(jsc41)*state % rho + &
      screened_rates(k_p_ti44__he4_sc41)*Y(jti44)*state % rho &
       )
    call set_jac_entry(state, jsc41, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca38__p_sc41)*Y(jca38)*state % rho + screened_rates(k_he4_k37__sc41)* &
      Y(jk37)*state % rho - screened_rates(k_he4_sc41__p_ti44)*Y(jsc41)* &
      state % rho - screened_rates(k_he4_sc41__v45)*Y(jsc41)*state % rho &
       )
    call set_jac_entry(state, jsc41, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_k37__sc41)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc41, jk37, scratch)

    scratch = (&
      screened_rates(k_he4_ca38__p_sc41)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc41, jca38, scratch)

    scratch = (&
      screened_rates(k_p_ca40__sc41)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc41, jca40, scratch)

    scratch = (&
      -screened_rates(k_he4_sc41__p_ti44)*Y(jhe4)*state % rho - screened_rates(k_he4_sc41__v45)* &
      Y(jhe4)*state % rho - screened_rates(k_p_sc41__he4_ca38)*Y(jp)*state % rho - &
      screened_rates(k_p_sc41__ti42)*Y(jp)*state % rho - &
      screened_rates(k_sc41__ca41__weak__wc12) - screened_rates(k_sc41__he4_k37) - &
      screened_rates(k_sc41__p_ca40) &
       )
    call set_jac_entry(state, jsc41, jsc41, scratch)

    scratch = (&
      screened_rates(k_ti41__sc41__weak__wc17) &
       )
    call set_jac_entry(state, jsc41, jti41, scratch)

    scratch = (&
      screened_rates(k_ti42__p_sc41) &
       )
    call set_jac_entry(state, jsc41, jti42, scratch)

    scratch = (&
      screened_rates(k_p_ti44__he4_sc41)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc41, jti44, scratch)

    scratch = (&
      screened_rates(k_v45__he4_sc41) &
       )
    call set_jac_entry(state, jsc41, jv45, scratch)

    scratch = (&
      screened_rates(k_cr43__p_p_sc41__weak__wc12) &
       )
    call set_jac_entry(state, jsc41, jcr43, scratch)

    scratch = (&
      screened_rates(k_p_ca41__sc42)*Y(jca41)*state % rho - screened_rates(k_p_sc42__he4_ca39)* &
      Y(jsc42)*state % rho - screened_rates(k_p_sc42__ti43)*Y(jsc42)*state % rho + &
      screened_rates(k_p_ti45__he4_sc42)*Y(jti45)*state % rho &
       )
    call set_jac_entry(state, jsc42, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca39__p_sc42)*Y(jca39)*state % rho + screened_rates(k_he4_k38__sc42)* &
      Y(jk38)*state % rho - screened_rates(k_he4_sc42__p_ti45)*Y(jsc42)* &
      state % rho - screened_rates(k_he4_sc42__v46)*Y(jsc42)*state % rho &
       )
    call set_jac_entry(state, jsc42, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_k38__sc42)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc42, jk38, scratch)

    scratch = (&
      screened_rates(k_he4_ca39__p_sc42)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc42, jca39, scratch)

    scratch = (&
      screened_rates(k_p_ca41__sc42)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc42, jca41, scratch)

    scratch = (&
      -screened_rates(k_he4_sc42__p_ti45)*Y(jhe4)*state % rho - screened_rates(k_he4_sc42__v46)* &
      Y(jhe4)*state % rho - screened_rates(k_p_sc42__he4_ca39)*Y(jp)*state % rho - &
      screened_rates(k_p_sc42__ti43)*Y(jp)*state % rho - &
      screened_rates(k_sc42__ca42__weak__wc12) - screened_rates(k_sc42__he4_k38) - &
      screened_rates(k_sc42__p_ca41) &
       )
    call set_jac_entry(state, jsc42, jsc42, scratch)

    scratch = (&
      screened_rates(k_ti42__sc42__weak__wc12) &
       )
    call set_jac_entry(state, jsc42, jti42, scratch)

    scratch = (&
      screened_rates(k_ti43__p_sc42) &
       )
    call set_jac_entry(state, jsc42, jti43, scratch)

    scratch = (&
      screened_rates(k_p_ti45__he4_sc42)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc42, jti45, scratch)

    scratch = (&
      screened_rates(k_v46__he4_sc42) &
       )
    call set_jac_entry(state, jsc42, jv46, scratch)

    scratch = (&
      screened_rates(k_p_ca42__sc43)*Y(jca42)*state % rho - screened_rates(k_p_sc43__he4_ca40)* &
      Y(jsc43)*state % rho - screened_rates(k_p_sc43__ti44)*Y(jsc43)*state % rho + &
      screened_rates(k_p_ti46__he4_sc43)*Y(jti46)*state % rho &
       )
    call set_jac_entry(state, jsc43, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca40__p_sc43)*Y(jca40)*state % rho + screened_rates(k_he4_k39__sc43)* &
      Y(jk39)*state % rho - screened_rates(k_he4_sc43__p_ti46)*Y(jsc43)* &
      state % rho - screened_rates(k_he4_sc43__v47)*Y(jsc43)*state % rho &
       )
    call set_jac_entry(state, jsc43, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_k39__sc43)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc43, jk39, scratch)

    scratch = (&
      screened_rates(k_he4_ca40__p_sc43)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc43, jca40, scratch)

    scratch = (&
      screened_rates(k_p_ca42__sc43)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc43, jca42, scratch)

    scratch = (&
      -screened_rates(k_he4_sc43__p_ti46)*Y(jhe4)*state % rho - screened_rates(k_he4_sc43__v47)* &
      Y(jhe4)*state % rho - screened_rates(k_p_sc43__he4_ca40)*Y(jp)*state % rho - &
      screened_rates(k_p_sc43__ti44)*Y(jp)*state % rho - &
      screened_rates(k_sc43__ca43__weak__wc12) - screened_rates(k_sc43__he4_k39) - &
      screened_rates(k_sc43__p_ca42) &
       )
    call set_jac_entry(state, jsc43, jsc43, scratch)

    scratch = (&
      screened_rates(k_ti43__sc43__weak__wc12) &
       )
    call set_jac_entry(state, jsc43, jti43, scratch)

    scratch = (&
      screened_rates(k_ti44__p_sc43) &
       )
    call set_jac_entry(state, jsc43, jti44, scratch)

    scratch = (&
      screened_rates(k_p_ti46__he4_sc43)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc43, jti46, scratch)

    scratch = (&
      screened_rates(k_v47__he4_sc43) &
       )
    call set_jac_entry(state, jsc43, jv47, scratch)

    scratch = (&
      screened_rates(k_p_ca43__sc44)*Y(jca43)*state % rho - screened_rates(k_p_sc44__he4_ca41)* &
      Y(jsc44)*state % rho - screened_rates(k_p_sc44__ti45)*Y(jsc44)*state % rho + &
      screened_rates(k_p_ti47__he4_sc44)*Y(jti47)*state % rho &
       )
    call set_jac_entry(state, jsc44, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca41__p_sc44)*Y(jca41)*state % rho + screened_rates(k_he4_k40__sc44)* &
      Y(jk40)*state % rho - screened_rates(k_he4_sc44__p_ti47)*Y(jsc44)* &
      state % rho - screened_rates(k_he4_sc44__v48)*Y(jsc44)*state % rho &
       )
    call set_jac_entry(state, jsc44, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_k40__sc44)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc44, jk40, scratch)

    scratch = (&
      screened_rates(k_he4_ca41__p_sc44)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc44, jca41, scratch)

    scratch = (&
      screened_rates(k_p_ca43__sc44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc44, jca43, scratch)

    scratch = (&
      -screened_rates(k_he4_sc44__p_ti47)*Y(jhe4)*state % rho - screened_rates(k_he4_sc44__v48)* &
      Y(jhe4)*state % rho - screened_rates(k_p_sc44__he4_ca41)*Y(jp)*state % rho - &
      screened_rates(k_p_sc44__ti45)*Y(jp)*state % rho - &
      screened_rates(k_sc44__ca44__weak__wc12) - screened_rates(k_sc44__he4_k40) - &
      screened_rates(k_sc44__p_ca43) &
       )
    call set_jac_entry(state, jsc44, jsc44, scratch)

    scratch = (&
      screened_rates(k_ti44__sc44__weak__wc12) &
       )
    call set_jac_entry(state, jsc44, jti44, scratch)

    scratch = (&
      screened_rates(k_ti45__p_sc44) &
       )
    call set_jac_entry(state, jsc44, jti45, scratch)

    scratch = (&
      screened_rates(k_p_ti47__he4_sc44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc44, jti47, scratch)

    scratch = (&
      screened_rates(k_v48__he4_sc44) &
       )
    call set_jac_entry(state, jsc44, jv48, scratch)

    scratch = (&
      screened_rates(k_p_ca44__sc45)*Y(jca44)*state % rho - screened_rates(k_p_sc45__he4_ca42)* &
      Y(jsc45)*state % rho - screened_rates(k_p_sc45__ti46)*Y(jsc45)*state % rho + &
      screened_rates(k_p_ti48__he4_sc45)*Y(jti48)*state % rho &
       )
    call set_jac_entry(state, jsc45, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca42__p_sc45)*Y(jca42)*state % rho + screened_rates(k_he4_k41__sc45)* &
      Y(jk41)*state % rho - screened_rates(k_he4_sc45__p_ti48)*Y(jsc45)* &
      state % rho - screened_rates(k_he4_sc45__v49)*Y(jsc45)*state % rho &
       )
    call set_jac_entry(state, jsc45, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_k41__sc45)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc45, jk41, scratch)

    scratch = (&
      screened_rates(k_he4_ca42__p_sc45)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc45, jca42, scratch)

    scratch = (&
      screened_rates(k_p_ca44__sc45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc45, jca44, scratch)

    scratch = (&
      -screened_rates(k_he4_sc45__p_ti48)*Y(jhe4)*state % rho - screened_rates(k_he4_sc45__v49)* &
      Y(jhe4)*state % rho - screened_rates(k_p_sc45__he4_ca42)*Y(jp)*state % rho - &
      screened_rates(k_p_sc45__ti46)*Y(jp)*state % rho - screened_rates(k_sc45__he4_k41) - &
      screened_rates(k_sc45__p_ca44) &
       )
    call set_jac_entry(state, jsc45, jsc45, scratch)

    scratch = (&
      screened_rates(k_ti45__sc45__weak__wc12) &
       )
    call set_jac_entry(state, jsc45, jti45, scratch)

    scratch = (&
      screened_rates(k_ti46__p_sc45) &
       )
    call set_jac_entry(state, jsc45, jti46, scratch)

    scratch = (&
      screened_rates(k_p_ti48__he4_sc45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc45, jti48, scratch)

    scratch = (&
      screened_rates(k_v49__he4_sc45) &
       )
    call set_jac_entry(state, jsc45, jv49, scratch)

    scratch = (&
      -screened_rates(k_p_sc46__he4_ca43)*Y(jsc46)*state % rho - screened_rates(k_p_sc46__ti47)* &
      Y(jsc46)*state % rho &
       )
    call set_jac_entry(state, jsc46, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca43__p_sc46)*Y(jca43)*state % rho - screened_rates(k_he4_sc46__v50)* &
      Y(jsc46)*state % rho &
       )
    call set_jac_entry(state, jsc46, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ca43__p_sc46)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jsc46, jca43, scratch)

    scratch = (&
      -screened_rates(k_he4_sc46__v50)*Y(jhe4)*state % rho - screened_rates(k_p_sc46__he4_ca43)* &
      Y(jp)*state % rho - screened_rates(k_p_sc46__ti47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jsc46, jsc46, scratch)

    scratch = (&
      screened_rates(k_ti47__p_sc46) &
       )
    call set_jac_entry(state, jsc46, jti47, scratch)

    scratch = (&
      screened_rates(k_v50__he4_sc46) &
       )
    call set_jac_entry(state, jsc46, jv50, scratch)

    scratch = (&
      screened_rates(k_p_sc37__ti38)*Y(jsc37)*state % rho + screened_rates(k_p_v41__he4_ti38)* &
      Y(jv41)*state % rho &
       )
    call set_jac_entry(state, jti38, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca34__ti38)*Y(jca34)*state % rho - screened_rates(k_he4_ti38__cr42)* &
      Y(jti38)*state % rho - screened_rates(k_he4_ti38__p_v41)*Y(jti38)* &
      state % rho &
       )
    call set_jac_entry(state, jti38, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ca34__ti38)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti38, jca34, scratch)

    scratch = (&
      screened_rates(k_p_sc37__ti38)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti38, jsc37, scratch)

    scratch = (&
      -screened_rates(k_he4_ti38__cr42)*Y(jhe4)*state % rho - screened_rates(k_he4_ti38__p_v41)* &
      Y(jhe4)*state % rho - screened_rates(k_ti38__he4_ca34) - &
      screened_rates(k_ti38__p_sc37) - screened_rates(k_ti38__sc38__weak__mo97) &
       )
    call set_jac_entry(state, jti38, jti38, scratch)

    scratch = (&
      screened_rates(k_p_v41__he4_ti38)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti38, jv41, scratch)

    scratch = (&
      screened_rates(k_cr42__he4_ti38) &
       )
    call set_jac_entry(state, jti38, jcr42, scratch)

    scratch = (&
      screened_rates(k_p_sc38__ti39)*Y(jsc38)*state % rho - screened_rates(k_p_ti39__he4_sc36)* &
      Y(jti39)*state % rho - screened_rates(k_p_ti39__v40)*Y(jti39)*state % rho + &
      screened_rates(k_p_v42__he4_ti39)*Y(jv42)*state % rho &
       )
    call set_jac_entry(state, jti39, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca35__ti39)*Y(jca35)*state % rho + screened_rates(k_he4_sc36__p_ti39)* &
      Y(jsc36)*state % rho - screened_rates(k_he4_ti39__cr43)*Y(jti39)*state % rho &
      - screened_rates(k_he4_ti39__p_v42)*Y(jti39)*state % rho &
       )
    call set_jac_entry(state, jti39, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ca35__ti39)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti39, jca35, scratch)

    scratch = (&
      screened_rates(k_he4_sc36__p_ti39)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti39, jsc36, scratch)

    scratch = (&
      screened_rates(k_p_sc38__ti39)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti39, jsc38, scratch)

    scratch = (&
      -screened_rates(k_he4_ti39__cr43)*Y(jhe4)*state % rho - screened_rates(k_he4_ti39__p_v42)* &
      Y(jhe4)*state % rho - screened_rates(k_p_ti39__he4_sc36)*Y(jp)*state % rho - &
      screened_rates(k_p_ti39__v40)*Y(jp)*state % rho - screened_rates(k_ti39__he4_ca35) - &
      screened_rates(k_ti39__p_ca38__weak__wc12) - screened_rates(k_ti39__p_sc38) - &
      screened_rates(k_ti39__sc39__weak__wc17) &
       )
    call set_jac_entry(state, jti39, jti39, scratch)

    scratch = (&
      screened_rates(k_v40__p_ti39) &
       )
    call set_jac_entry(state, jti39, jv40, scratch)

    scratch = (&
      screened_rates(k_p_v42__he4_ti39)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti39, jv42, scratch)

    scratch = (&
      screened_rates(k_cr43__he4_ti39) &
       )
    call set_jac_entry(state, jti39, jcr43, scratch)

    scratch = (&
      screened_rates(k_p_sc39__ti40)*Y(jsc39)*state % rho - screened_rates(k_p_ti40__he4_sc37)* &
      Y(jti40)*state % rho - screened_rates(k_p_ti40__v41)*Y(jti40)*state % rho + &
      screened_rates(k_p_v43__he4_ti40)*Y(jv43)*state % rho &
       )
    call set_jac_entry(state, jti40, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca36__ti40)*Y(jca36)*state % rho + screened_rates(k_he4_sc37__p_ti40)* &
      Y(jsc37)*state % rho - screened_rates(k_he4_ti40__cr44)*Y(jti40)*state % rho &
      - screened_rates(k_he4_ti40__p_v43)*Y(jti40)*state % rho &
       )
    call set_jac_entry(state, jti40, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ca36__ti40)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti40, jca36, scratch)

    scratch = (&
      screened_rates(k_he4_sc37__p_ti40)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti40, jsc37, scratch)

    scratch = (&
      screened_rates(k_p_sc39__ti40)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti40, jsc39, scratch)

    scratch = (&
      -screened_rates(k_he4_ti40__cr44)*Y(jhe4)*state % rho - screened_rates(k_he4_ti40__p_v43)* &
      Y(jhe4)*state % rho - screened_rates(k_p_ti40__he4_sc37)*Y(jp)*state % rho - &
      screened_rates(k_p_ti40__v41)*Y(jp)*state % rho - screened_rates(k_ti40__he4_ca36) - &
      screened_rates(k_ti40__p_ca39__weak__wc12) - screened_rates(k_ti40__p_sc39) - &
      screened_rates(k_ti40__sc40__weak__wc17) &
       )
    call set_jac_entry(state, jti40, jti40, scratch)

    scratch = (&
      screened_rates(k_v40__ti40__weak__bqa_pos_) &
       )
    call set_jac_entry(state, jti40, jv40, scratch)

    scratch = (&
      screened_rates(k_v41__p_ti40) &
       )
    call set_jac_entry(state, jti40, jv41, scratch)

    scratch = (&
      screened_rates(k_p_v43__he4_ti40)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti40, jv43, scratch)

    scratch = (&
      screened_rates(k_cr44__he4_ti40) &
       )
    call set_jac_entry(state, jti40, jcr44, scratch)

    scratch = (&
      screened_rates(k_p_sc40__ti41)*Y(jsc40)*state % rho - screened_rates(k_p_ti41__he4_sc38)* &
      Y(jti41)*state % rho - screened_rates(k_p_ti41__v42)*Y(jti41)*state % rho + &
      screened_rates(k_p_v44__he4_ti41)*Y(jv44)*state % rho &
       )
    call set_jac_entry(state, jti41, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca37__ti41)*Y(jca37)*state % rho + screened_rates(k_he4_sc38__p_ti41)* &
      Y(jsc38)*state % rho - screened_rates(k_he4_ti41__cr45)*Y(jti41)*state % rho &
      - screened_rates(k_he4_ti41__p_v44)*Y(jti41)*state % rho &
       )
    call set_jac_entry(state, jti41, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ca37__ti41)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti41, jca37, scratch)

    scratch = (&
      screened_rates(k_he4_sc38__p_ti41)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti41, jsc38, scratch)

    scratch = (&
      screened_rates(k_p_sc40__ti41)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti41, jsc40, scratch)

    scratch = (&
      -screened_rates(k_he4_ti41__cr45)*Y(jhe4)*state % rho - screened_rates(k_he4_ti41__p_v44)* &
      Y(jhe4)*state % rho - screened_rates(k_p_ti41__he4_sc38)*Y(jp)*state % rho - &
      screened_rates(k_p_ti41__v42)*Y(jp)*state % rho - screened_rates(k_ti41__he4_ca37) - &
      screened_rates(k_ti41__p_ca40__weak__wc12) - screened_rates(k_ti41__p_sc40) - &
      screened_rates(k_ti41__sc41__weak__wc17) &
       )
    call set_jac_entry(state, jti41, jti41, scratch)

    scratch = (&
      screened_rates(k_v41__ti41__weak__bqa_pos_) &
       )
    call set_jac_entry(state, jti41, jv41, scratch)

    scratch = (&
      screened_rates(k_v42__p_ti41) &
       )
    call set_jac_entry(state, jti41, jv42, scratch)

    scratch = (&
      screened_rates(k_p_v44__he4_ti41)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti41, jv44, scratch)

    scratch = (&
      screened_rates(k_cr42__p_ti41__weak__wc12) &
       )
    call set_jac_entry(state, jti41, jcr42, scratch)

    scratch = (&
      screened_rates(k_cr45__he4_ti41) &
       )
    call set_jac_entry(state, jti41, jcr45, scratch)

    scratch = (&
      screened_rates(k_p_sc41__ti42)*Y(jsc41)*state % rho - screened_rates(k_p_ti42__he4_sc39)* &
      Y(jti42)*state % rho - screened_rates(k_p_ti42__v43)*Y(jti42)*state % rho + &
      screened_rates(k_p_v45__he4_ti42)*Y(jv45)*state % rho &
       )
    call set_jac_entry(state, jti42, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca38__ti42)*Y(jca38)*state % rho + screened_rates(k_he4_sc39__p_ti42)* &
      Y(jsc39)*state % rho - screened_rates(k_he4_ti42__cr46)*Y(jti42)*state % rho &
      - screened_rates(k_he4_ti42__p_v45)*Y(jti42)*state % rho &
       )
    call set_jac_entry(state, jti42, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ca38__ti42)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti42, jca38, scratch)

    scratch = (&
      screened_rates(k_he4_sc39__p_ti42)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti42, jsc39, scratch)

    scratch = (&
      screened_rates(k_p_sc41__ti42)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti42, jsc41, scratch)

    scratch = (&
      -screened_rates(k_he4_ti42__cr46)*Y(jhe4)*state % rho - screened_rates(k_he4_ti42__p_v45)* &
      Y(jhe4)*state % rho - screened_rates(k_p_ti42__he4_sc39)*Y(jp)*state % rho - &
      screened_rates(k_p_ti42__v43)*Y(jp)*state % rho - screened_rates(k_ti42__he4_ca38) - &
      screened_rates(k_ti42__p_sc41) - screened_rates(k_ti42__sc42__weak__wc12) &
       )
    call set_jac_entry(state, jti42, jti42, scratch)

    scratch = (&
      screened_rates(k_v42__ti42__weak__mo97) &
       )
    call set_jac_entry(state, jti42, jv42, scratch)

    scratch = (&
      screened_rates(k_v43__p_ti42) &
       )
    call set_jac_entry(state, jti42, jv43, scratch)

    scratch = (&
      screened_rates(k_p_v45__he4_ti42)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti42, jv45, scratch)

    scratch = (&
      screened_rates(k_cr43__p_ti42__weak__wc17) &
       )
    call set_jac_entry(state, jti42, jcr43, scratch)

    scratch = (&
      screened_rates(k_cr46__he4_ti42) &
       )
    call set_jac_entry(state, jti42, jcr46, scratch)

    scratch = (&
      screened_rates(k_fe45__p_p_p_ti42__weak__wc12) &
       )
    call set_jac_entry(state, jti42, jfe45, scratch)

    scratch = (&
      screened_rates(k_p_sc42__ti43)*Y(jsc42)*state % rho - screened_rates(k_p_ti43__he4_sc40)* &
      Y(jti43)*state % rho - screened_rates(k_p_ti43__v44)*Y(jti43)*state % rho + &
      screened_rates(k_p_v46__he4_ti43)*Y(jv46)*state % rho &
       )
    call set_jac_entry(state, jti43, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca39__ti43)*Y(jca39)*state % rho + screened_rates(k_he4_sc40__p_ti43)* &
      Y(jsc40)*state % rho - screened_rates(k_he4_ti43__cr47)*Y(jti43)*state % rho &
      - screened_rates(k_he4_ti43__p_v46)*Y(jti43)*state % rho &
       )
    call set_jac_entry(state, jti43, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ca39__ti43)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti43, jca39, scratch)

    scratch = (&
      screened_rates(k_he4_sc40__p_ti43)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti43, jsc40, scratch)

    scratch = (&
      screened_rates(k_p_sc42__ti43)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti43, jsc42, scratch)

    scratch = (&
      -screened_rates(k_he4_ti43__cr47)*Y(jhe4)*state % rho - screened_rates(k_he4_ti43__p_v46)* &
      Y(jhe4)*state % rho - screened_rates(k_p_ti43__he4_sc40)*Y(jp)*state % rho - &
      screened_rates(k_p_ti43__v44)*Y(jp)*state % rho - screened_rates(k_ti43__he4_ca39) - &
      screened_rates(k_ti43__p_sc42) - screened_rates(k_ti43__sc43__weak__wc12) &
       )
    call set_jac_entry(state, jti43, jti43, scratch)

    scratch = (&
      screened_rates(k_v43__ti43__weak__wc12) &
       )
    call set_jac_entry(state, jti43, jv43, scratch)

    scratch = (&
      screened_rates(k_v44__p_ti43) &
       )
    call set_jac_entry(state, jti43, jv44, scratch)

    scratch = (&
      screened_rates(k_p_v46__he4_ti43)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti43, jv46, scratch)

    scratch = (&
      screened_rates(k_cr44__p_ti43__weak__wc12) &
       )
    call set_jac_entry(state, jti43, jcr44, scratch)

    scratch = (&
      screened_rates(k_cr47__he4_ti43) &
       )
    call set_jac_entry(state, jti43, jcr47, scratch)

    scratch = (&
      screened_rates(k_p_sc43__ti44)*Y(jsc43)*state % rho - screened_rates(k_p_ti44__he4_sc41)* &
      Y(jti44)*state % rho - screened_rates(k_p_ti44__v45)*Y(jti44)*state % rho + &
      screened_rates(k_p_v47__he4_ti44)*Y(jv47)*state % rho &
       )
    call set_jac_entry(state, jti44, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca40__ti44)*Y(jca40)*state % rho + screened_rates(k_he4_sc41__p_ti44)* &
      Y(jsc41)*state % rho - screened_rates(k_he4_ti44__cr48)*Y(jti44)*state % rho &
      - screened_rates(k_he4_ti44__p_v47)*Y(jti44)*state % rho &
       )
    call set_jac_entry(state, jti44, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ca40__ti44)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti44, jca40, scratch)

    scratch = (&
      screened_rates(k_he4_sc41__p_ti44)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti44, jsc41, scratch)

    scratch = (&
      screened_rates(k_p_sc43__ti44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti44, jsc43, scratch)

    scratch = (&
      -screened_rates(k_he4_ti44__cr48)*Y(jhe4)*state % rho - screened_rates(k_he4_ti44__p_v47)* &
      Y(jhe4)*state % rho - screened_rates(k_p_ti44__he4_sc41)*Y(jp)*state % rho - &
      screened_rates(k_p_ti44__v45)*Y(jp)*state % rho - screened_rates(k_ti44__he4_ca40) - &
      screened_rates(k_ti44__p_sc43) - screened_rates(k_ti44__sc44__weak__wc12) &
       )
    call set_jac_entry(state, jti44, jti44, scratch)

    scratch = (&
      screened_rates(k_v44__ti44__weak__wc12) &
       )
    call set_jac_entry(state, jti44, jv44, scratch)

    scratch = (&
      screened_rates(k_v45__p_ti44) &
       )
    call set_jac_entry(state, jti44, jv45, scratch)

    scratch = (&
      screened_rates(k_p_v47__he4_ti44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti44, jv47, scratch)

    scratch = (&
      screened_rates(k_cr45__p_ti44__weak__wc12) &
       )
    call set_jac_entry(state, jti44, jcr45, scratch)

    scratch = (&
      screened_rates(k_cr48__he4_ti44) &
       )
    call set_jac_entry(state, jti44, jcr48, scratch)

    scratch = (&
      screened_rates(k_p_sc44__ti45)*Y(jsc44)*state % rho - screened_rates(k_p_ti45__he4_sc42)* &
      Y(jti45)*state % rho - screened_rates(k_p_ti45__v46)*Y(jti45)*state % rho + &
      screened_rates(k_p_v48__he4_ti45)*Y(jv48)*state % rho &
       )
    call set_jac_entry(state, jti45, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca41__ti45)*Y(jca41)*state % rho + screened_rates(k_he4_sc42__p_ti45)* &
      Y(jsc42)*state % rho - screened_rates(k_he4_ti45__cr49)*Y(jti45)*state % rho &
      - screened_rates(k_he4_ti45__p_v48)*Y(jti45)*state % rho &
       )
    call set_jac_entry(state, jti45, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ca41__ti45)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti45, jca41, scratch)

    scratch = (&
      screened_rates(k_he4_sc42__p_ti45)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti45, jsc42, scratch)

    scratch = (&
      screened_rates(k_p_sc44__ti45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti45, jsc44, scratch)

    scratch = (&
      -screened_rates(k_he4_ti45__cr49)*Y(jhe4)*state % rho - screened_rates(k_he4_ti45__p_v48)* &
      Y(jhe4)*state % rho - screened_rates(k_p_ti45__he4_sc42)*Y(jp)*state % rho - &
      screened_rates(k_p_ti45__v46)*Y(jp)*state % rho - screened_rates(k_ti45__he4_ca41) - &
      screened_rates(k_ti45__p_sc44) - screened_rates(k_ti45__sc45__weak__wc12) &
       )
    call set_jac_entry(state, jti45, jti45, scratch)

    scratch = (&
      screened_rates(k_v45__ti45__weak__wc12) &
       )
    call set_jac_entry(state, jti45, jv45, scratch)

    scratch = (&
      screened_rates(k_v46__p_ti45) &
       )
    call set_jac_entry(state, jti45, jv46, scratch)

    scratch = (&
      screened_rates(k_p_v48__he4_ti45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti45, jv48, scratch)

    scratch = (&
      screened_rates(k_cr49__he4_ti45) &
       )
    call set_jac_entry(state, jti45, jcr49, scratch)

    scratch = (&
      screened_rates(k_p_sc45__ti46)*Y(jsc45)*state % rho - screened_rates(k_p_ti46__he4_sc43)* &
      Y(jti46)*state % rho - screened_rates(k_p_ti46__v47)*Y(jti46)*state % rho + &
      screened_rates(k_p_v49__he4_ti46)*Y(jv49)*state % rho &
       )
    call set_jac_entry(state, jti46, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca42__ti46)*Y(jca42)*state % rho + screened_rates(k_he4_sc43__p_ti46)* &
      Y(jsc43)*state % rho - screened_rates(k_he4_ti46__cr50)*Y(jti46)*state % rho &
      - screened_rates(k_he4_ti46__p_v49)*Y(jti46)*state % rho &
       )
    call set_jac_entry(state, jti46, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ca42__ti46)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti46, jca42, scratch)

    scratch = (&
      screened_rates(k_he4_sc43__p_ti46)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti46, jsc43, scratch)

    scratch = (&
      screened_rates(k_p_sc45__ti46)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti46, jsc45, scratch)

    scratch = (&
      -screened_rates(k_he4_ti46__cr50)*Y(jhe4)*state % rho - screened_rates(k_he4_ti46__p_v49)* &
      Y(jhe4)*state % rho - screened_rates(k_p_ti46__he4_sc43)*Y(jp)*state % rho - &
      screened_rates(k_p_ti46__v47)*Y(jp)*state % rho - screened_rates(k_ti46__he4_ca42) - &
      screened_rates(k_ti46__p_sc45) &
       )
    call set_jac_entry(state, jti46, jti46, scratch)

    scratch = (&
      screened_rates(k_v46__ti46__weak__wc12) &
       )
    call set_jac_entry(state, jti46, jv46, scratch)

    scratch = (&
      screened_rates(k_v47__p_ti46) &
       )
    call set_jac_entry(state, jti46, jv47, scratch)

    scratch = (&
      screened_rates(k_p_v49__he4_ti46)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti46, jv49, scratch)

    scratch = (&
      screened_rates(k_cr50__he4_ti46) &
       )
    call set_jac_entry(state, jti46, jcr50, scratch)

    scratch = (&
      screened_rates(k_p_sc46__ti47)*Y(jsc46)*state % rho - screened_rates(k_p_ti47__he4_sc44)* &
      Y(jti47)*state % rho - screened_rates(k_p_ti47__v48)*Y(jti47)*state % rho + &
      screened_rates(k_p_v50__he4_ti47)*Y(jv50)*state % rho &
       )
    call set_jac_entry(state, jti47, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca43__ti47)*Y(jca43)*state % rho + screened_rates(k_he4_sc44__p_ti47)* &
      Y(jsc44)*state % rho - screened_rates(k_he4_ti47__cr51)*Y(jti47)*state % rho &
      - screened_rates(k_he4_ti47__p_v50)*Y(jti47)*state % rho &
       )
    call set_jac_entry(state, jti47, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ca43__ti47)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti47, jca43, scratch)

    scratch = (&
      screened_rates(k_he4_sc44__p_ti47)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti47, jsc44, scratch)

    scratch = (&
      screened_rates(k_p_sc46__ti47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti47, jsc46, scratch)

    scratch = (&
      -screened_rates(k_he4_ti47__cr51)*Y(jhe4)*state % rho - screened_rates(k_he4_ti47__p_v50)* &
      Y(jhe4)*state % rho - screened_rates(k_p_ti47__he4_sc44)*Y(jp)*state % rho - &
      screened_rates(k_p_ti47__v48)*Y(jp)*state % rho - screened_rates(k_ti47__he4_ca43) - &
      screened_rates(k_ti47__p_sc46) &
       )
    call set_jac_entry(state, jti47, jti47, scratch)

    scratch = (&
      screened_rates(k_v47__ti47__weak__wc12) &
       )
    call set_jac_entry(state, jti47, jv47, scratch)

    scratch = (&
      screened_rates(k_v48__p_ti47) &
       )
    call set_jac_entry(state, jti47, jv48, scratch)

    scratch = (&
      screened_rates(k_p_v50__he4_ti47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jti47, jv50, scratch)

    scratch = (&
      screened_rates(k_cr51__he4_ti47) &
       )
    call set_jac_entry(state, jti47, jcr51, scratch)

    scratch = (&
      -screened_rates(k_p_ti48__he4_sc45)*Y(jti48)*state % rho - screened_rates(k_p_ti48__v49)* &
      Y(jti48)*state % rho &
       )
    call set_jac_entry(state, jti48, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca44__ti48)*Y(jca44)*state % rho + screened_rates(k_he4_sc45__p_ti48)* &
      Y(jsc45)*state % rho - screened_rates(k_he4_ti48__cr52)*Y(jti48)*state % rho &
       )
    call set_jac_entry(state, jti48, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ca44__ti48)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti48, jca44, scratch)

    scratch = (&
      screened_rates(k_he4_sc45__p_ti48)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jti48, jsc45, scratch)

    scratch = (&
      -screened_rates(k_he4_ti48__cr52)*Y(jhe4)*state % rho - screened_rates(k_p_ti48__he4_sc45)* &
      Y(jp)*state % rho - screened_rates(k_p_ti48__v49)*Y(jp)*state % rho - &
      screened_rates(k_ti48__he4_ca44) &
       )
    call set_jac_entry(state, jti48, jti48, scratch)

    scratch = (&
      screened_rates(k_v48__ti48__weak__wc12) &
       )
    call set_jac_entry(state, jti48, jv48, scratch)

    scratch = (&
      screened_rates(k_v49__p_ti48) &
       )
    call set_jac_entry(state, jti48, jv49, scratch)

    scratch = (&
      screened_rates(k_cr52__he4_ti48) &
       )
    call set_jac_entry(state, jti48, jcr52, scratch)

    scratch = (&
      screened_rates(k_p_cr43__he4_v40)*Y(jcr43)*state % rho + screened_rates(k_p_ti39__v40)* &
      Y(jti39)*state % rho &
       )
    call set_jac_entry(state, jv40, jp, scratch)

    scratch = (&
      screened_rates(k_he4_sc36__v40)*Y(jsc36)*state % rho - screened_rates(k_he4_v40__mn44)* &
      Y(jv40)*state % rho - screened_rates(k_he4_v40__p_cr43)*Y(jv40)*state % rho &
       )
    call set_jac_entry(state, jv40, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_sc36__v40)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv40, jsc36, scratch)

    scratch = (&
      screened_rates(k_p_ti39__v40)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv40, jti39, scratch)

    scratch = (&
      -screened_rates(k_he4_v40__mn44)*Y(jhe4)*state % rho - screened_rates(k_he4_v40__p_cr43)* &
      Y(jhe4)*state % rho - screened_rates(k_v40__he4_sc36) - &
      screened_rates(k_v40__p_ti39) - screened_rates(k_v40__ti40__weak__bqa_pos_) &
       )
    call set_jac_entry(state, jv40, jv40, scratch)

    scratch = (&
      screened_rates(k_p_cr43__he4_v40)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv40, jcr43, scratch)

    scratch = (&
      screened_rates(k_mn44__he4_v40) &
       )
    call set_jac_entry(state, jv40, jmn44, scratch)

    scratch = (&
      screened_rates(k_p_cr44__he4_v41)*Y(jcr44)*state % rho + screened_rates(k_p_ti40__v41)* &
      Y(jti40)*state % rho - screened_rates(k_p_v41__cr42)*Y(jv41)*state % rho - &
      screened_rates(k_p_v41__he4_ti38)*Y(jv41)*state % rho &
       )
    call set_jac_entry(state, jv41, jp, scratch)

    scratch = (&
      screened_rates(k_he4_sc37__v41)*Y(jsc37)*state % rho + screened_rates(k_he4_ti38__p_v41)* &
      Y(jti38)*state % rho - screened_rates(k_he4_v41__mn45)*Y(jv41)*state % rho - &
      screened_rates(k_he4_v41__p_cr44)*Y(jv41)*state % rho &
       )
    call set_jac_entry(state, jv41, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_sc37__v41)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv41, jsc37, scratch)

    scratch = (&
      screened_rates(k_he4_ti38__p_v41)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv41, jti38, scratch)

    scratch = (&
      screened_rates(k_p_ti40__v41)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv41, jti40, scratch)

    scratch = (&
      -screened_rates(k_he4_v41__mn45)*Y(jhe4)*state % rho - screened_rates(k_he4_v41__p_cr44)* &
      Y(jhe4)*state % rho - screened_rates(k_p_v41__cr42)*Y(jp)*state % rho - &
      screened_rates(k_p_v41__he4_ti38)*Y(jp)*state % rho - &
      screened_rates(k_v41__he4_sc37) - screened_rates(k_v41__p_ti40) - &
      screened_rates(k_v41__ti41__weak__bqa_pos_) &
       )
    call set_jac_entry(state, jv41, jv41, scratch)

    scratch = (&
      screened_rates(k_cr42__p_v41) &
       )
    call set_jac_entry(state, jv41, jcr42, scratch)

    scratch = (&
      screened_rates(k_p_cr44__he4_v41)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv41, jcr44, scratch)

    scratch = (&
      screened_rates(k_mn45__he4_v41) &
       )
    call set_jac_entry(state, jv41, jmn45, scratch)

    scratch = (&
      screened_rates(k_p_cr45__he4_v42)*Y(jcr45)*state % rho + screened_rates(k_p_ti41__v42)* &
      Y(jti41)*state % rho - screened_rates(k_p_v42__cr43)*Y(jv42)*state % rho - &
      screened_rates(k_p_v42__he4_ti39)*Y(jv42)*state % rho &
       )
    call set_jac_entry(state, jv42, jp, scratch)

    scratch = (&
      screened_rates(k_he4_sc38__v42)*Y(jsc38)*state % rho + screened_rates(k_he4_ti39__p_v42)* &
      Y(jti39)*state % rho - screened_rates(k_he4_v42__mn46)*Y(jv42)*state % rho - &
      screened_rates(k_he4_v42__p_cr45)*Y(jv42)*state % rho &
       )
    call set_jac_entry(state, jv42, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_sc38__v42)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv42, jsc38, scratch)

    scratch = (&
      screened_rates(k_he4_ti39__p_v42)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv42, jti39, scratch)

    scratch = (&
      screened_rates(k_p_ti41__v42)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv42, jti41, scratch)

    scratch = (&
      -screened_rates(k_he4_v42__mn46)*Y(jhe4)*state % rho - screened_rates(k_he4_v42__p_cr45)* &
      Y(jhe4)*state % rho - screened_rates(k_p_v42__cr43)*Y(jp)*state % rho - &
      screened_rates(k_p_v42__he4_ti39)*Y(jp)*state % rho - &
      screened_rates(k_v42__he4_sc38) - screened_rates(k_v42__p_ti41) - &
      screened_rates(k_v42__ti42__weak__mo97) &
       )
    call set_jac_entry(state, jv42, jv42, scratch)

    scratch = (&
      screened_rates(k_cr42__v42__weak__wc12) &
       )
    call set_jac_entry(state, jv42, jcr42, scratch)

    scratch = (&
      screened_rates(k_cr43__p_v42) &
       )
    call set_jac_entry(state, jv42, jcr43, scratch)

    scratch = (&
      screened_rates(k_p_cr45__he4_v42)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv42, jcr45, scratch)

    scratch = (&
      screened_rates(k_mn46__he4_v42) &
       )
    call set_jac_entry(state, jv42, jmn46, scratch)

    scratch = (&
      screened_rates(k_p_cr46__he4_v43)*Y(jcr46)*state % rho + screened_rates(k_p_ti42__v43)* &
      Y(jti42)*state % rho - screened_rates(k_p_v43__cr44)*Y(jv43)*state % rho - &
      screened_rates(k_p_v43__he4_ti40)*Y(jv43)*state % rho &
       )
    call set_jac_entry(state, jv43, jp, scratch)

    scratch = (&
      screened_rates(k_he4_sc39__v43)*Y(jsc39)*state % rho + screened_rates(k_he4_ti40__p_v43)* &
      Y(jti40)*state % rho - screened_rates(k_he4_v43__mn47)*Y(jv43)*state % rho - &
      screened_rates(k_he4_v43__p_cr46)*Y(jv43)*state % rho &
       )
    call set_jac_entry(state, jv43, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_sc39__v43)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv43, jsc39, scratch)

    scratch = (&
      screened_rates(k_he4_ti40__p_v43)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv43, jti40, scratch)

    scratch = (&
      screened_rates(k_p_ti42__v43)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv43, jti42, scratch)

    scratch = (&
      -screened_rates(k_he4_v43__mn47)*Y(jhe4)*state % rho - screened_rates(k_he4_v43__p_cr46)* &
      Y(jhe4)*state % rho - screened_rates(k_p_v43__cr44)*Y(jp)*state % rho - &
      screened_rates(k_p_v43__he4_ti40)*Y(jp)*state % rho - &
      screened_rates(k_v43__he4_sc39) - screened_rates(k_v43__p_ti42) - &
      screened_rates(k_v43__ti43__weak__wc12) &
       )
    call set_jac_entry(state, jv43, jv43, scratch)

    scratch = (&
      screened_rates(k_cr43__v43__weak__wc12) &
       )
    call set_jac_entry(state, jv43, jcr43, scratch)

    scratch = (&
      screened_rates(k_cr44__p_v43) &
       )
    call set_jac_entry(state, jv43, jcr44, scratch)

    scratch = (&
      screened_rates(k_p_cr46__he4_v43)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv43, jcr46, scratch)

    scratch = (&
      screened_rates(k_mn47__he4_v43) &
       )
    call set_jac_entry(state, jv43, jmn47, scratch)

    scratch = (&
      screened_rates(k_fe45__p_p_v43__weak__wc12) &
       )
    call set_jac_entry(state, jv43, jfe45, scratch)

    scratch = (&
      screened_rates(k_p_cr47__he4_v44)*Y(jcr47)*state % rho + screened_rates(k_p_ti43__v44)* &
      Y(jti43)*state % rho - screened_rates(k_p_v44__cr45)*Y(jv44)*state % rho - &
      screened_rates(k_p_v44__he4_ti41)*Y(jv44)*state % rho &
       )
    call set_jac_entry(state, jv44, jp, scratch)

    scratch = (&
      screened_rates(k_he4_sc40__v44)*Y(jsc40)*state % rho + screened_rates(k_he4_ti41__p_v44)* &
      Y(jti41)*state % rho - screened_rates(k_he4_v44__mn48)*Y(jv44)*state % rho - &
      screened_rates(k_he4_v44__p_cr47)*Y(jv44)*state % rho &
       )
    call set_jac_entry(state, jv44, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_sc40__v44)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv44, jsc40, scratch)

    scratch = (&
      screened_rates(k_he4_ti41__p_v44)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv44, jti41, scratch)

    scratch = (&
      screened_rates(k_p_ti43__v44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv44, jti43, scratch)

    scratch = (&
      -screened_rates(k_he4_v44__mn48)*Y(jhe4)*state % rho - screened_rates(k_he4_v44__p_cr47)* &
      Y(jhe4)*state % rho - screened_rates(k_p_v44__cr45)*Y(jp)*state % rho - &
      screened_rates(k_p_v44__he4_ti41)*Y(jp)*state % rho - &
      screened_rates(k_v44__he4_sc40) - screened_rates(k_v44__p_ti43) - &
      screened_rates(k_v44__ti44__weak__wc12) &
       )
    call set_jac_entry(state, jv44, jv44, scratch)

    scratch = (&
      screened_rates(k_cr44__v44__weak__wc12) &
       )
    call set_jac_entry(state, jv44, jcr44, scratch)

    scratch = (&
      screened_rates(k_cr45__p_v44) &
       )
    call set_jac_entry(state, jv44, jcr45, scratch)

    scratch = (&
      screened_rates(k_p_cr47__he4_v44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv44, jcr47, scratch)

    scratch = (&
      screened_rates(k_mn48__he4_v44) &
       )
    call set_jac_entry(state, jv44, jmn48, scratch)

    scratch = (&
      screened_rates(k_p_cr48__he4_v45)*Y(jcr48)*state % rho + screened_rates(k_p_ti44__v45)* &
      Y(jti44)*state % rho - screened_rates(k_p_v45__cr46)*Y(jv45)*state % rho - &
      screened_rates(k_p_v45__he4_ti42)*Y(jv45)*state % rho &
       )
    call set_jac_entry(state, jv45, jp, scratch)

    scratch = (&
      screened_rates(k_he4_sc41__v45)*Y(jsc41)*state % rho + screened_rates(k_he4_ti42__p_v45)* &
      Y(jti42)*state % rho - screened_rates(k_he4_v45__mn49)*Y(jv45)*state % rho - &
      screened_rates(k_he4_v45__p_cr48)*Y(jv45)*state % rho &
       )
    call set_jac_entry(state, jv45, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_sc41__v45)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv45, jsc41, scratch)

    scratch = (&
      screened_rates(k_he4_ti42__p_v45)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv45, jti42, scratch)

    scratch = (&
      screened_rates(k_p_ti44__v45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv45, jti44, scratch)

    scratch = (&
      -screened_rates(k_he4_v45__mn49)*Y(jhe4)*state % rho - screened_rates(k_he4_v45__p_cr48)* &
      Y(jhe4)*state % rho - screened_rates(k_p_v45__cr46)*Y(jp)*state % rho - &
      screened_rates(k_p_v45__he4_ti42)*Y(jp)*state % rho - &
      screened_rates(k_v45__he4_sc41) - screened_rates(k_v45__p_ti44) - &
      screened_rates(k_v45__ti45__weak__wc12) &
       )
    call set_jac_entry(state, jv45, jv45, scratch)

    scratch = (&
      screened_rates(k_cr45__v45__weak__wc12) &
       )
    call set_jac_entry(state, jv45, jcr45, scratch)

    scratch = (&
      screened_rates(k_cr46__p_v45) &
       )
    call set_jac_entry(state, jv45, jcr46, scratch)

    scratch = (&
      screened_rates(k_p_cr48__he4_v45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv45, jcr48, scratch)

    scratch = (&
      screened_rates(k_mn46__p_v45__weak__wc12) &
       )
    call set_jac_entry(state, jv45, jmn46, scratch)

    scratch = (&
      screened_rates(k_mn49__he4_v45) &
       )
    call set_jac_entry(state, jv45, jmn49, scratch)

    scratch = (&
      screened_rates(k_p_cr49__he4_v46)*Y(jcr49)*state % rho + screened_rates(k_p_ti45__v46)* &
      Y(jti45)*state % rho - screened_rates(k_p_v46__cr47)*Y(jv46)*state % rho - &
      screened_rates(k_p_v46__he4_ti43)*Y(jv46)*state % rho &
       )
    call set_jac_entry(state, jv46, jp, scratch)

    scratch = (&
      screened_rates(k_he4_sc42__v46)*Y(jsc42)*state % rho + screened_rates(k_he4_ti43__p_v46)* &
      Y(jti43)*state % rho - screened_rates(k_he4_v46__mn50)*Y(jv46)*state % rho - &
      screened_rates(k_he4_v46__p_cr49)*Y(jv46)*state % rho &
       )
    call set_jac_entry(state, jv46, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_sc42__v46)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv46, jsc42, scratch)

    scratch = (&
      screened_rates(k_he4_ti43__p_v46)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv46, jti43, scratch)

    scratch = (&
      screened_rates(k_p_ti45__v46)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv46, jti45, scratch)

    scratch = (&
      -screened_rates(k_he4_v46__mn50)*Y(jhe4)*state % rho - screened_rates(k_he4_v46__p_cr49)* &
      Y(jhe4)*state % rho - screened_rates(k_p_v46__cr47)*Y(jp)*state % rho - &
      screened_rates(k_p_v46__he4_ti43)*Y(jp)*state % rho - &
      screened_rates(k_v46__he4_sc42) - screened_rates(k_v46__p_ti45) - &
      screened_rates(k_v46__ti46__weak__wc12) &
       )
    call set_jac_entry(state, jv46, jv46, scratch)

    scratch = (&
      screened_rates(k_cr46__v46__weak__wc12) &
       )
    call set_jac_entry(state, jv46, jcr46, scratch)

    scratch = (&
      screened_rates(k_cr47__p_v46) &
       )
    call set_jac_entry(state, jv46, jcr47, scratch)

    scratch = (&
      screened_rates(k_p_cr49__he4_v46)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv46, jcr49, scratch)

    scratch = (&
      screened_rates(k_mn50__he4_v46) &
       )
    call set_jac_entry(state, jv46, jmn50, scratch)

    scratch = (&
      screened_rates(k_p_cr50__he4_v47)*Y(jcr50)*state % rho + screened_rates(k_p_ti46__v47)* &
      Y(jti46)*state % rho - screened_rates(k_p_v47__cr48)*Y(jv47)*state % rho - &
      screened_rates(k_p_v47__he4_ti44)*Y(jv47)*state % rho &
       )
    call set_jac_entry(state, jv47, jp, scratch)

    scratch = (&
      screened_rates(k_he4_sc43__v47)*Y(jsc43)*state % rho + screened_rates(k_he4_ti44__p_v47)* &
      Y(jti44)*state % rho - screened_rates(k_he4_v47__mn51)*Y(jv47)*state % rho - &
      screened_rates(k_he4_v47__p_cr50)*Y(jv47)*state % rho &
       )
    call set_jac_entry(state, jv47, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_sc43__v47)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv47, jsc43, scratch)

    scratch = (&
      screened_rates(k_he4_ti44__p_v47)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv47, jti44, scratch)

    scratch = (&
      screened_rates(k_p_ti46__v47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv47, jti46, scratch)

    scratch = (&
      -screened_rates(k_he4_v47__mn51)*Y(jhe4)*state % rho - screened_rates(k_he4_v47__p_cr50)* &
      Y(jhe4)*state % rho - screened_rates(k_p_v47__cr48)*Y(jp)*state % rho - &
      screened_rates(k_p_v47__he4_ti44)*Y(jp)*state % rho - &
      screened_rates(k_v47__he4_sc43) - screened_rates(k_v47__p_ti46) - &
      screened_rates(k_v47__ti47__weak__wc12) &
       )
    call set_jac_entry(state, jv47, jv47, scratch)

    scratch = (&
      screened_rates(k_cr47__v47__weak__wc12) &
       )
    call set_jac_entry(state, jv47, jcr47, scratch)

    scratch = (&
      screened_rates(k_cr48__p_v47) &
       )
    call set_jac_entry(state, jv47, jcr48, scratch)

    scratch = (&
      screened_rates(k_p_cr50__he4_v47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv47, jcr50, scratch)

    scratch = (&
      screened_rates(k_mn48__p_v47__weak__wc12) &
       )
    call set_jac_entry(state, jv47, jmn48, scratch)

    scratch = (&
      screened_rates(k_mn51__he4_v47) &
       )
    call set_jac_entry(state, jv47, jmn51, scratch)

    scratch = (&
      screened_rates(k_p_cr51__he4_v48)*Y(jcr51)*state % rho + screened_rates(k_p_ti47__v48)* &
      Y(jti47)*state % rho - screened_rates(k_p_v48__cr49)*Y(jv48)*state % rho - &
      screened_rates(k_p_v48__he4_ti45)*Y(jv48)*state % rho &
       )
    call set_jac_entry(state, jv48, jp, scratch)

    scratch = (&
      screened_rates(k_he4_sc44__v48)*Y(jsc44)*state % rho + screened_rates(k_he4_ti45__p_v48)* &
      Y(jti45)*state % rho - screened_rates(k_he4_v48__mn52)*Y(jv48)*state % rho - &
      screened_rates(k_he4_v48__p_cr51)*Y(jv48)*state % rho &
       )
    call set_jac_entry(state, jv48, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_sc44__v48)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv48, jsc44, scratch)

    scratch = (&
      screened_rates(k_he4_ti45__p_v48)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv48, jti45, scratch)

    scratch = (&
      screened_rates(k_p_ti47__v48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv48, jti47, scratch)

    scratch = (&
      -screened_rates(k_he4_v48__mn52)*Y(jhe4)*state % rho - screened_rates(k_he4_v48__p_cr51)* &
      Y(jhe4)*state % rho - screened_rates(k_p_v48__cr49)*Y(jp)*state % rho - &
      screened_rates(k_p_v48__he4_ti45)*Y(jp)*state % rho - &
      screened_rates(k_v48__he4_sc44) - screened_rates(k_v48__p_ti47) - &
      screened_rates(k_v48__ti48__weak__wc12) &
       )
    call set_jac_entry(state, jv48, jv48, scratch)

    scratch = (&
      screened_rates(k_cr48__v48__weak__wc12) &
       )
    call set_jac_entry(state, jv48, jcr48, scratch)

    scratch = (&
      screened_rates(k_cr49__p_v48) &
       )
    call set_jac_entry(state, jv48, jcr49, scratch)

    scratch = (&
      screened_rates(k_p_cr51__he4_v48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv48, jcr51, scratch)

    scratch = (&
      screened_rates(k_mn52__he4_v48) &
       )
    call set_jac_entry(state, jv48, jmn52, scratch)

    scratch = (&
      screened_rates(k_p_cr52__he4_v49)*Y(jcr52)*state % rho + screened_rates(k_p_ti48__v49)* &
      Y(jti48)*state % rho - screened_rates(k_p_v49__cr50)*Y(jv49)*state % rho - &
      screened_rates(k_p_v49__he4_ti46)*Y(jv49)*state % rho &
       )
    call set_jac_entry(state, jv49, jp, scratch)

    scratch = (&
      screened_rates(k_he4_sc45__v49)*Y(jsc45)*state % rho + screened_rates(k_he4_ti46__p_v49)* &
      Y(jti46)*state % rho - screened_rates(k_he4_v49__mn53)*Y(jv49)*state % rho - &
      screened_rates(k_he4_v49__p_cr52)*Y(jv49)*state % rho &
       )
    call set_jac_entry(state, jv49, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_sc45__v49)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv49, jsc45, scratch)

    scratch = (&
      screened_rates(k_he4_ti46__p_v49)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv49, jti46, scratch)

    scratch = (&
      screened_rates(k_p_ti48__v49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv49, jti48, scratch)

    scratch = (&
      -screened_rates(k_he4_v49__mn53)*Y(jhe4)*state % rho - screened_rates(k_he4_v49__p_cr52)* &
      Y(jhe4)*state % rho - screened_rates(k_p_v49__cr50)*Y(jp)*state % rho - &
      screened_rates(k_p_v49__he4_ti46)*Y(jp)*state % rho - &
      screened_rates(k_v49__he4_sc45) - screened_rates(k_v49__p_ti48) &
       )
    call set_jac_entry(state, jv49, jv49, scratch)

    scratch = (&
      screened_rates(k_cr49__v49__weak__wc12) &
       )
    call set_jac_entry(state, jv49, jcr49, scratch)

    scratch = (&
      screened_rates(k_cr50__p_v49) &
       )
    call set_jac_entry(state, jv49, jcr50, scratch)

    scratch = (&
      screened_rates(k_p_cr52__he4_v49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jv49, jcr52, scratch)

    scratch = (&
      screened_rates(k_mn53__he4_v49) &
       )
    call set_jac_entry(state, jv49, jmn53, scratch)

    scratch = (&
      -screened_rates(k_p_v50__cr51)*Y(jv50)*state % rho - screened_rates(k_p_v50__he4_ti47)* &
      Y(jv50)*state % rho &
       )
    call set_jac_entry(state, jv50, jp, scratch)

    scratch = (&
      screened_rates(k_he4_sc46__v50)*Y(jsc46)*state % rho + screened_rates(k_he4_ti47__p_v50)* &
      Y(jti47)*state % rho - screened_rates(k_he4_v50__mn54)*Y(jv50)*state % rho &
       )
    call set_jac_entry(state, jv50, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_sc46__v50)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv50, jsc46, scratch)

    scratch = (&
      screened_rates(k_he4_ti47__p_v50)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jv50, jti47, scratch)

    scratch = (&
      -screened_rates(k_he4_v50__mn54)*Y(jhe4)*state % rho - screened_rates(k_p_v50__cr51)* &
      Y(jp)*state % rho - screened_rates(k_p_v50__he4_ti47)*Y(jp)*state % rho - &
      screened_rates(k_v50__he4_sc46) &
       )
    call set_jac_entry(state, jv50, jv50, scratch)

    scratch = (&
      screened_rates(k_cr51__p_v50) &
       )
    call set_jac_entry(state, jv50, jcr51, scratch)

    scratch = (&
      screened_rates(k_mn54__he4_v50) &
       )
    call set_jac_entry(state, jv50, jmn54, scratch)

    scratch = (&
      screened_rates(k_p_mn45__he4_cr42)*Y(jmn45)*state % rho + screened_rates(k_p_v41__cr42)* &
      Y(jv41)*state % rho &
       )
    call set_jac_entry(state, jcr42, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cr42__fe46)*Y(jcr42)*state % rho - screened_rates(k_he4_cr42__p_mn45) &
      *Y(jcr42)*state % rho + screened_rates(k_he4_ti38__cr42)*Y(jti38)* &
      state % rho &
       )
    call set_jac_entry(state, jcr42, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ti38__cr42)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr42, jti38, scratch)

    scratch = (&
      screened_rates(k_p_v41__cr42)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr42, jv41, scratch)

    scratch = (&
      -screened_rates(k_cr42__he4_ti38) - screened_rates(k_cr42__p_ti41__weak__wc12) - &
      screened_rates(k_cr42__p_v41) - screened_rates(k_cr42__v42__weak__wc12) - &
      screened_rates(k_he4_cr42__fe46)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr42__p_mn45)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr42, jcr42, scratch)

    scratch = (&
      screened_rates(k_p_mn45__he4_cr42)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr42, jmn45, scratch)

    scratch = (&
      screened_rates(k_fe46__he4_cr42) &
       )
    call set_jac_entry(state, jcr42, jfe46, scratch)

    scratch = (&
      -screened_rates(k_p_cr43__he4_v40)*Y(jcr43)*state % rho - screened_rates(k_p_cr43__mn44)* &
      Y(jcr43)*state % rho + screened_rates(k_p_mn46__he4_cr43)*Y(jmn46)* &
      state % rho + screened_rates(k_p_v42__cr43)*Y(jv42)*state % rho &
       )
    call set_jac_entry(state, jcr43, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cr43__fe47)*Y(jcr43)*state % rho - screened_rates(k_he4_cr43__p_mn46) &
      *Y(jcr43)*state % rho + screened_rates(k_he4_ti39__cr43)*Y(jti39)* &
      state % rho + screened_rates(k_he4_v40__p_cr43)*Y(jv40)*state % rho &
       )
    call set_jac_entry(state, jcr43, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ti39__cr43)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr43, jti39, scratch)

    scratch = (&
      screened_rates(k_he4_v40__p_cr43)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr43, jv40, scratch)

    scratch = (&
      screened_rates(k_p_v42__cr43)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr43, jv42, scratch)

    scratch = (&
      -screened_rates(k_cr43__he4_ti39) - screened_rates(k_cr43__p_p_p_ca40__weak__wc12) - &
      screened_rates(k_cr43__p_p_sc41__weak__wc12) - &
      screened_rates(k_cr43__p_ti42__weak__wc17) - screened_rates(k_cr43__p_v42) - &
      screened_rates(k_cr43__v43__weak__wc12) - screened_rates(k_he4_cr43__fe47)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_cr43__p_mn46)*Y(jhe4)*state % rho &
      - screened_rates(k_p_cr43__he4_v40)*Y(jp)*state % rho - &
      screened_rates(k_p_cr43__mn44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr43, jcr43, scratch)

    scratch = (&
      screened_rates(k_mn44__p_cr43) &
       )
    call set_jac_entry(state, jcr43, jmn44, scratch)

    scratch = (&
      screened_rates(k_p_mn46__he4_cr43)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr43, jmn46, scratch)

    scratch = (&
      screened_rates(k_fe47__he4_cr43) &
       )
    call set_jac_entry(state, jcr43, jfe47, scratch)

    scratch = (&
      -screened_rates(k_p_cr44__he4_v41)*Y(jcr44)*state % rho - screened_rates(k_p_cr44__mn45)* &
      Y(jcr44)*state % rho + screened_rates(k_p_mn47__he4_cr44)*Y(jmn47)* &
      state % rho + screened_rates(k_p_v43__cr44)*Y(jv43)*state % rho &
       )
    call set_jac_entry(state, jcr44, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cr44__fe48)*Y(jcr44)*state % rho - screened_rates(k_he4_cr44__p_mn47) &
      *Y(jcr44)*state % rho + screened_rates(k_he4_ti40__cr44)*Y(jti40)* &
      state % rho + screened_rates(k_he4_v41__p_cr44)*Y(jv41)*state % rho &
       )
    call set_jac_entry(state, jcr44, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ti40__cr44)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr44, jti40, scratch)

    scratch = (&
      screened_rates(k_he4_v41__p_cr44)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr44, jv41, scratch)

    scratch = (&
      screened_rates(k_p_v43__cr44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr44, jv43, scratch)

    scratch = (&
      -screened_rates(k_cr44__he4_ti40) - screened_rates(k_cr44__p_ti43__weak__wc12) - &
      screened_rates(k_cr44__p_v43) - screened_rates(k_cr44__v44__weak__wc12) - &
      screened_rates(k_he4_cr44__fe48)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr44__p_mn47)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cr44__he4_v41)*Y(jp)*state % rho - screened_rates(k_p_cr44__mn45) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr44, jcr44, scratch)

    scratch = (&
      screened_rates(k_mn44__cr44__weak__bqa_pos_) &
       )
    call set_jac_entry(state, jcr44, jmn44, scratch)

    scratch = (&
      screened_rates(k_mn45__p_cr44) &
       )
    call set_jac_entry(state, jcr44, jmn45, scratch)

    scratch = (&
      screened_rates(k_p_mn47__he4_cr44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr44, jmn47, scratch)

    scratch = (&
      screened_rates(k_fe45__p_cr44__weak__wc17) &
       )
    call set_jac_entry(state, jcr44, jfe45, scratch)

    scratch = (&
      screened_rates(k_fe48__he4_cr44) &
       )
    call set_jac_entry(state, jcr44, jfe48, scratch)

    scratch = (&
      -screened_rates(k_p_cr45__he4_v42)*Y(jcr45)*state % rho - screened_rates(k_p_cr45__mn46)* &
      Y(jcr45)*state % rho + screened_rates(k_p_mn48__he4_cr45)*Y(jmn48)* &
      state % rho + screened_rates(k_p_v44__cr45)*Y(jv44)*state % rho &
       )
    call set_jac_entry(state, jcr45, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cr45__fe49)*Y(jcr45)*state % rho - screened_rates(k_he4_cr45__p_mn48) &
      *Y(jcr45)*state % rho + screened_rates(k_he4_ti41__cr45)*Y(jti41)* &
      state % rho + screened_rates(k_he4_v42__p_cr45)*Y(jv42)*state % rho &
       )
    call set_jac_entry(state, jcr45, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ti41__cr45)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr45, jti41, scratch)

    scratch = (&
      screened_rates(k_he4_v42__p_cr45)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr45, jv42, scratch)

    scratch = (&
      screened_rates(k_p_v44__cr45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr45, jv44, scratch)

    scratch = (&
      -screened_rates(k_cr45__he4_ti41) - screened_rates(k_cr45__p_ti44__weak__wc12) - &
      screened_rates(k_cr45__p_v44) - screened_rates(k_cr45__v45__weak__wc12) - &
      screened_rates(k_he4_cr45__fe49)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr45__p_mn48)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cr45__he4_v42)*Y(jp)*state % rho - screened_rates(k_p_cr45__mn46) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr45, jcr45, scratch)

    scratch = (&
      screened_rates(k_mn45__cr45__weak__bqa_pos_) &
       )
    call set_jac_entry(state, jcr45, jmn45, scratch)

    scratch = (&
      screened_rates(k_mn46__p_cr45) &
       )
    call set_jac_entry(state, jcr45, jmn46, scratch)

    scratch = (&
      screened_rates(k_p_mn48__he4_cr45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr45, jmn48, scratch)

    scratch = (&
      screened_rates(k_fe46__p_cr45__weak__wc12) &
       )
    call set_jac_entry(state, jcr45, jfe46, scratch)

    scratch = (&
      screened_rates(k_fe49__he4_cr45) &
       )
    call set_jac_entry(state, jcr45, jfe49, scratch)

    scratch = (&
      -screened_rates(k_p_cr46__he4_v43)*Y(jcr46)*state % rho - screened_rates(k_p_cr46__mn47)* &
      Y(jcr46)*state % rho + screened_rates(k_p_mn49__he4_cr46)*Y(jmn49)* &
      state % rho + screened_rates(k_p_v45__cr46)*Y(jv45)*state % rho &
       )
    call set_jac_entry(state, jcr46, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cr46__fe50)*Y(jcr46)*state % rho - screened_rates(k_he4_cr46__p_mn49) &
      *Y(jcr46)*state % rho + screened_rates(k_he4_ti42__cr46)*Y(jti42)* &
      state % rho + screened_rates(k_he4_v43__p_cr46)*Y(jv43)*state % rho &
       )
    call set_jac_entry(state, jcr46, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ti42__cr46)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr46, jti42, scratch)

    scratch = (&
      screened_rates(k_he4_v43__p_cr46)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr46, jv43, scratch)

    scratch = (&
      screened_rates(k_p_v45__cr46)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr46, jv45, scratch)

    scratch = (&
      -screened_rates(k_cr46__he4_ti42) - screened_rates(k_cr46__p_v45) - &
      screened_rates(k_cr46__v46__weak__wc12) - screened_rates(k_he4_cr46__fe50)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_cr46__p_mn49)*Y(jhe4)*state % rho &
      - screened_rates(k_p_cr46__he4_v43)*Y(jp)*state % rho - &
      screened_rates(k_p_cr46__mn47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr46, jcr46, scratch)

    scratch = (&
      screened_rates(k_mn46__cr46__weak__wc12) &
       )
    call set_jac_entry(state, jcr46, jmn46, scratch)

    scratch = (&
      screened_rates(k_mn47__p_cr46) &
       )
    call set_jac_entry(state, jcr46, jmn47, scratch)

    scratch = (&
      screened_rates(k_p_mn49__he4_cr46)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr46, jmn49, scratch)

    scratch = (&
      screened_rates(k_fe47__p_cr46__weak__wc12) &
       )
    call set_jac_entry(state, jcr46, jfe47, scratch)

    scratch = (&
      screened_rates(k_fe50__he4_cr46) &
       )
    call set_jac_entry(state, jcr46, jfe50, scratch)

    scratch = (&
      -screened_rates(k_p_cr47__he4_v44)*Y(jcr47)*state % rho - screened_rates(k_p_cr47__mn48)* &
      Y(jcr47)*state % rho + screened_rates(k_p_mn50__he4_cr47)*Y(jmn50)* &
      state % rho + screened_rates(k_p_v46__cr47)*Y(jv46)*state % rho &
       )
    call set_jac_entry(state, jcr47, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cr47__fe51)*Y(jcr47)*state % rho - screened_rates(k_he4_cr47__p_mn50) &
      *Y(jcr47)*state % rho + screened_rates(k_he4_ti43__cr47)*Y(jti43)* &
      state % rho + screened_rates(k_he4_v44__p_cr47)*Y(jv44)*state % rho &
       )
    call set_jac_entry(state, jcr47, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ti43__cr47)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr47, jti43, scratch)

    scratch = (&
      screened_rates(k_he4_v44__p_cr47)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr47, jv44, scratch)

    scratch = (&
      screened_rates(k_p_v46__cr47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr47, jv46, scratch)

    scratch = (&
      -screened_rates(k_cr47__he4_ti43) - screened_rates(k_cr47__p_v46) - &
      screened_rates(k_cr47__v47__weak__wc12) - screened_rates(k_he4_cr47__fe51)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_cr47__p_mn50)*Y(jhe4)*state % rho &
      - screened_rates(k_p_cr47__he4_v44)*Y(jp)*state % rho - &
      screened_rates(k_p_cr47__mn48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr47, jcr47, scratch)

    scratch = (&
      screened_rates(k_mn47__cr47__weak__wc17) &
       )
    call set_jac_entry(state, jcr47, jmn47, scratch)

    scratch = (&
      screened_rates(k_mn48__p_cr47) &
       )
    call set_jac_entry(state, jcr47, jmn48, scratch)

    scratch = (&
      screened_rates(k_p_mn50__he4_cr47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr47, jmn50, scratch)

    scratch = (&
      screened_rates(k_fe48__p_cr47__weak__wc12) &
       )
    call set_jac_entry(state, jcr47, jfe48, scratch)

    scratch = (&
      screened_rates(k_fe51__he4_cr47) &
       )
    call set_jac_entry(state, jcr47, jfe51, scratch)

    scratch = (&
      -screened_rates(k_p_cr48__he4_v45)*Y(jcr48)*state % rho - screened_rates(k_p_cr48__mn49)* &
      Y(jcr48)*state % rho + screened_rates(k_p_mn51__he4_cr48)*Y(jmn51)* &
      state % rho + screened_rates(k_p_v47__cr48)*Y(jv47)*state % rho &
       )
    call set_jac_entry(state, jcr48, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cr48__fe52)*Y(jcr48)*state % rho - screened_rates(k_he4_cr48__p_mn51) &
      *Y(jcr48)*state % rho + screened_rates(k_he4_ti44__cr48)*Y(jti44)* &
      state % rho + screened_rates(k_he4_v45__p_cr48)*Y(jv45)*state % rho &
       )
    call set_jac_entry(state, jcr48, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ti44__cr48)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr48, jti44, scratch)

    scratch = (&
      screened_rates(k_he4_v45__p_cr48)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr48, jv45, scratch)

    scratch = (&
      screened_rates(k_p_v47__cr48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr48, jv47, scratch)

    scratch = (&
      -screened_rates(k_cr48__he4_ti44) - screened_rates(k_cr48__p_v47) - &
      screened_rates(k_cr48__v48__weak__wc12) - screened_rates(k_he4_cr48__fe52)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_cr48__p_mn51)*Y(jhe4)*state % rho &
      - screened_rates(k_p_cr48__he4_v45)*Y(jp)*state % rho - &
      screened_rates(k_p_cr48__mn49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr48, jcr48, scratch)

    scratch = (&
      screened_rates(k_mn48__cr48__weak__wc12) &
       )
    call set_jac_entry(state, jcr48, jmn48, scratch)

    scratch = (&
      screened_rates(k_mn49__p_cr48) &
       )
    call set_jac_entry(state, jcr48, jmn49, scratch)

    scratch = (&
      screened_rates(k_p_mn51__he4_cr48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr48, jmn51, scratch)

    scratch = (&
      screened_rates(k_fe49__p_cr48__weak__wc12) &
       )
    call set_jac_entry(state, jcr48, jfe49, scratch)

    scratch = (&
      screened_rates(k_fe52__he4_cr48) &
       )
    call set_jac_entry(state, jcr48, jfe52, scratch)

    scratch = (&
      -screened_rates(k_p_cr49__he4_v46)*Y(jcr49)*state % rho - screened_rates(k_p_cr49__mn50)* &
      Y(jcr49)*state % rho + screened_rates(k_p_mn52__he4_cr49)*Y(jmn52)* &
      state % rho + screened_rates(k_p_v48__cr49)*Y(jv48)*state % rho &
       )
    call set_jac_entry(state, jcr49, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cr49__fe53)*Y(jcr49)*state % rho - screened_rates(k_he4_cr49__p_mn52) &
      *Y(jcr49)*state % rho + screened_rates(k_he4_ti45__cr49)*Y(jti45)* &
      state % rho + screened_rates(k_he4_v46__p_cr49)*Y(jv46)*state % rho &
       )
    call set_jac_entry(state, jcr49, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ti45__cr49)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr49, jti45, scratch)

    scratch = (&
      screened_rates(k_he4_v46__p_cr49)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr49, jv46, scratch)

    scratch = (&
      screened_rates(k_p_v48__cr49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr49, jv48, scratch)

    scratch = (&
      -screened_rates(k_cr49__he4_ti45) - screened_rates(k_cr49__p_v48) - &
      screened_rates(k_cr49__v49__weak__wc12) - screened_rates(k_he4_cr49__fe53)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_cr49__p_mn52)*Y(jhe4)*state % rho &
      - screened_rates(k_p_cr49__he4_v46)*Y(jp)*state % rho - &
      screened_rates(k_p_cr49__mn50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr49, jcr49, scratch)

    scratch = (&
      screened_rates(k_mn49__cr49__weak__wc12) &
       )
    call set_jac_entry(state, jcr49, jmn49, scratch)

    scratch = (&
      screened_rates(k_mn50__p_cr49) &
       )
    call set_jac_entry(state, jcr49, jmn50, scratch)

    scratch = (&
      screened_rates(k_p_mn52__he4_cr49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr49, jmn52, scratch)

    scratch = (&
      screened_rates(k_fe53__he4_cr49) &
       )
    call set_jac_entry(state, jcr49, jfe53, scratch)

    scratch = (&
      -screened_rates(k_p_cr50__he4_v47)*Y(jcr50)*state % rho - screened_rates(k_p_cr50__mn51)* &
      Y(jcr50)*state % rho + screened_rates(k_p_mn53__he4_cr50)*Y(jmn53)* &
      state % rho + screened_rates(k_p_v49__cr50)*Y(jv49)*state % rho &
       )
    call set_jac_entry(state, jcr50, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cr50__fe54)*Y(jcr50)*state % rho - screened_rates(k_he4_cr50__p_mn53) &
      *Y(jcr50)*state % rho + screened_rates(k_he4_ti46__cr50)*Y(jti46)* &
      state % rho + screened_rates(k_he4_v47__p_cr50)*Y(jv47)*state % rho &
       )
    call set_jac_entry(state, jcr50, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ti46__cr50)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr50, jti46, scratch)

    scratch = (&
      screened_rates(k_he4_v47__p_cr50)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr50, jv47, scratch)

    scratch = (&
      screened_rates(k_p_v49__cr50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr50, jv49, scratch)

    scratch = (&
      -screened_rates(k_cr50__he4_ti46) - screened_rates(k_cr50__p_v49) - &
      screened_rates(k_he4_cr50__fe54)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr50__p_mn53)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cr50__he4_v47)*Y(jp)*state % rho - screened_rates(k_p_cr50__mn51) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr50, jcr50, scratch)

    scratch = (&
      screened_rates(k_mn50__cr50__weak__wc12) &
       )
    call set_jac_entry(state, jcr50, jmn50, scratch)

    scratch = (&
      screened_rates(k_mn51__p_cr50) &
       )
    call set_jac_entry(state, jcr50, jmn51, scratch)

    scratch = (&
      screened_rates(k_p_mn53__he4_cr50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr50, jmn53, scratch)

    scratch = (&
      screened_rates(k_fe54__he4_cr50) &
       )
    call set_jac_entry(state, jcr50, jfe54, scratch)

    scratch = (&
      -screened_rates(k_p_cr51__he4_v48)*Y(jcr51)*state % rho - screened_rates(k_p_cr51__mn52)* &
      Y(jcr51)*state % rho + screened_rates(k_p_mn54__he4_cr51)*Y(jmn54)* &
      state % rho + screened_rates(k_p_v50__cr51)*Y(jv50)*state % rho &
       )
    call set_jac_entry(state, jcr51, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cr51__fe55)*Y(jcr51)*state % rho - screened_rates(k_he4_cr51__p_mn54) &
      *Y(jcr51)*state % rho + screened_rates(k_he4_ti47__cr51)*Y(jti47)* &
      state % rho + screened_rates(k_he4_v48__p_cr51)*Y(jv48)*state % rho &
       )
    call set_jac_entry(state, jcr51, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ti47__cr51)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr51, jti47, scratch)

    scratch = (&
      screened_rates(k_he4_v48__p_cr51)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr51, jv48, scratch)

    scratch = (&
      screened_rates(k_p_v50__cr51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr51, jv50, scratch)

    scratch = (&
      -screened_rates(k_cr51__he4_ti47) - screened_rates(k_cr51__p_v50) - &
      screened_rates(k_he4_cr51__fe55)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr51__p_mn54)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cr51__he4_v48)*Y(jp)*state % rho - screened_rates(k_p_cr51__mn52) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr51, jcr51, scratch)

    scratch = (&
      screened_rates(k_mn51__cr51__weak__wc12) &
       )
    call set_jac_entry(state, jcr51, jmn51, scratch)

    scratch = (&
      screened_rates(k_mn52__p_cr51) &
       )
    call set_jac_entry(state, jcr51, jmn52, scratch)

    scratch = (&
      screened_rates(k_p_mn54__he4_cr51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr51, jmn54, scratch)

    scratch = (&
      screened_rates(k_fe55__he4_cr51) &
       )
    call set_jac_entry(state, jcr51, jfe55, scratch)

    scratch = (&
      -screened_rates(k_p_cr52__he4_v49)*Y(jcr52)*state % rho - screened_rates(k_p_cr52__mn53)* &
      Y(jcr52)*state % rho + screened_rates(k_p_mn55__he4_cr52)*Y(jmn55)* &
      state % rho &
       )
    call set_jac_entry(state, jcr52, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cr52__fe56)*Y(jcr52)*state % rho - screened_rates(k_he4_cr52__p_mn55) &
      *Y(jcr52)*state % rho + screened_rates(k_he4_ti48__cr52)*Y(jti48)* &
      state % rho + screened_rates(k_he4_v49__p_cr52)*Y(jv49)*state % rho &
       )
    call set_jac_entry(state, jcr52, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ti48__cr52)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr52, jti48, scratch)

    scratch = (&
      screened_rates(k_he4_v49__p_cr52)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jcr52, jv49, scratch)

    scratch = (&
      -screened_rates(k_cr52__he4_ti48) - screened_rates(k_he4_cr52__fe56)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr52__p_mn55)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cr52__he4_v49)*Y(jp)*state % rho - screened_rates(k_p_cr52__mn53) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr52, jcr52, scratch)

    scratch = (&
      screened_rates(k_mn52__cr52__weak__wc12) &
       )
    call set_jac_entry(state, jcr52, jmn52, scratch)

    scratch = (&
      screened_rates(k_mn53__p_cr52) &
       )
    call set_jac_entry(state, jcr52, jmn53, scratch)

    scratch = (&
      screened_rates(k_p_mn55__he4_cr52)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jcr52, jmn55, scratch)

    scratch = (&
      screened_rates(k_fe56__he4_cr52) &
       )
    call set_jac_entry(state, jcr52, jfe56, scratch)

    scratch = (&
      screened_rates(k_p_cr43__mn44)*Y(jcr43)*state % rho + screened_rates(k_p_fe47__he4_mn44)* &
      Y(jfe47)*state % rho - screened_rates(k_p_mn44__fe45)*Y(jmn44)*state % rho &
       )
    call set_jac_entry(state, jmn44, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_mn44__co48)*Y(jmn44)*state % rho - screened_rates(k_he4_mn44__p_fe47) &
      *Y(jmn44)*state % rho + screened_rates(k_he4_v40__mn44)*Y(jv40)*state % rho &
       )
    call set_jac_entry(state, jmn44, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_v40__mn44)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn44, jv40, scratch)

    scratch = (&
      screened_rates(k_p_cr43__mn44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn44, jcr43, scratch)

    scratch = (&
      -screened_rates(k_he4_mn44__co48)*Y(jhe4)*state % rho - screened_rates(k_he4_mn44__p_fe47)* &
      Y(jhe4)*state % rho - screened_rates(k_mn44__cr44__weak__bqa_pos_) - &
      screened_rates(k_mn44__he4_v40) - screened_rates(k_mn44__p_cr43) - &
      screened_rates(k_p_mn44__fe45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn44, jmn44, scratch)

    scratch = (&
      screened_rates(k_fe45__p_mn44) &
       )
    call set_jac_entry(state, jmn44, jfe45, scratch)

    scratch = (&
      screened_rates(k_p_fe47__he4_mn44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn44, jfe47, scratch)

    scratch = (&
      screened_rates(k_co48__he4_mn44) &
       )
    call set_jac_entry(state, jmn44, jco48, scratch)

    scratch = (&
      screened_rates(k_p_cr44__mn45)*Y(jcr44)*state % rho + screened_rates(k_p_fe48__he4_mn45)* &
      Y(jfe48)*state % rho - screened_rates(k_p_mn45__fe46)*Y(jmn45)*state % rho - &
      screened_rates(k_p_mn45__he4_cr42)*Y(jmn45)*state % rho &
       )
    call set_jac_entry(state, jmn45, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr42__p_mn45)*Y(jcr42)*state % rho - screened_rates(k_he4_mn45__co49)* &
      Y(jmn45)*state % rho - screened_rates(k_he4_mn45__p_fe48)*Y(jmn45)* &
      state % rho + screened_rates(k_he4_v41__mn45)*Y(jv41)*state % rho &
       )
    call set_jac_entry(state, jmn45, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_v41__mn45)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn45, jv41, scratch)

    scratch = (&
      screened_rates(k_he4_cr42__p_mn45)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn45, jcr42, scratch)

    scratch = (&
      screened_rates(k_p_cr44__mn45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn45, jcr44, scratch)

    scratch = (&
      -screened_rates(k_he4_mn45__co49)*Y(jhe4)*state % rho - screened_rates(k_he4_mn45__p_fe48)* &
      Y(jhe4)*state % rho - screened_rates(k_mn45__cr45__weak__bqa_pos_) - &
      screened_rates(k_mn45__he4_v41) - screened_rates(k_mn45__p_cr44) - &
      screened_rates(k_p_mn45__fe46)*Y(jp)*state % rho - &
      screened_rates(k_p_mn45__he4_cr42)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn45, jmn45, scratch)

    scratch = (&
      screened_rates(k_fe45__mn45__weak__wc17) &
       )
    call set_jac_entry(state, jmn45, jfe45, scratch)

    scratch = (&
      screened_rates(k_fe46__p_mn45) &
       )
    call set_jac_entry(state, jmn45, jfe46, scratch)

    scratch = (&
      screened_rates(k_p_fe48__he4_mn45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn45, jfe48, scratch)

    scratch = (&
      screened_rates(k_co49__he4_mn45) &
       )
    call set_jac_entry(state, jmn45, jco49, scratch)

    scratch = (&
      screened_rates(k_p_cr45__mn46)*Y(jcr45)*state % rho + screened_rates(k_p_fe49__he4_mn46)* &
      Y(jfe49)*state % rho - screened_rates(k_p_mn46__fe47)*Y(jmn46)*state % rho - &
      screened_rates(k_p_mn46__he4_cr43)*Y(jmn46)*state % rho &
       )
    call set_jac_entry(state, jmn46, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr43__p_mn46)*Y(jcr43)*state % rho - screened_rates(k_he4_mn46__co50)* &
      Y(jmn46)*state % rho - screened_rates(k_he4_mn46__p_fe49)*Y(jmn46)* &
      state % rho + screened_rates(k_he4_v42__mn46)*Y(jv42)*state % rho &
       )
    call set_jac_entry(state, jmn46, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_v42__mn46)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn46, jv42, scratch)

    scratch = (&
      screened_rates(k_he4_cr43__p_mn46)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn46, jcr43, scratch)

    scratch = (&
      screened_rates(k_p_cr45__mn46)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn46, jcr45, scratch)

    scratch = (&
      -screened_rates(k_he4_mn46__co50)*Y(jhe4)*state % rho - screened_rates(k_he4_mn46__p_fe49)* &
      Y(jhe4)*state % rho - screened_rates(k_mn46__cr46__weak__wc12) - &
      screened_rates(k_mn46__he4_v42) - screened_rates(k_mn46__p_cr45) - &
      screened_rates(k_mn46__p_v45__weak__wc12) - screened_rates(k_p_mn46__fe47)*Y(jp)* &
      state % rho - screened_rates(k_p_mn46__he4_cr43)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn46, jmn46, scratch)

    scratch = (&
      screened_rates(k_fe46__mn46__weak__wc12) &
       )
    call set_jac_entry(state, jmn46, jfe46, scratch)

    scratch = (&
      screened_rates(k_fe47__p_mn46) &
       )
    call set_jac_entry(state, jmn46, jfe47, scratch)

    scratch = (&
      screened_rates(k_p_fe49__he4_mn46)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn46, jfe49, scratch)

    scratch = (&
      screened_rates(k_co50__he4_mn46) &
       )
    call set_jac_entry(state, jmn46, jco50, scratch)

    scratch = (&
      screened_rates(k_p_cr46__mn47)*Y(jcr46)*state % rho + screened_rates(k_p_fe50__he4_mn47)* &
      Y(jfe50)*state % rho - screened_rates(k_p_mn47__fe48)*Y(jmn47)*state % rho - &
      screened_rates(k_p_mn47__he4_cr44)*Y(jmn47)*state % rho &
       )
    call set_jac_entry(state, jmn47, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr44__p_mn47)*Y(jcr44)*state % rho - screened_rates(k_he4_mn47__co51)* &
      Y(jmn47)*state % rho - screened_rates(k_he4_mn47__p_fe50)*Y(jmn47)* &
      state % rho + screened_rates(k_he4_v43__mn47)*Y(jv43)*state % rho &
       )
    call set_jac_entry(state, jmn47, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_v43__mn47)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn47, jv43, scratch)

    scratch = (&
      screened_rates(k_he4_cr44__p_mn47)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn47, jcr44, scratch)

    scratch = (&
      screened_rates(k_p_cr46__mn47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn47, jcr46, scratch)

    scratch = (&
      -screened_rates(k_he4_mn47__co51)*Y(jhe4)*state % rho - screened_rates(k_he4_mn47__p_fe50)* &
      Y(jhe4)*state % rho - screened_rates(k_mn47__cr47__weak__wc17) - &
      screened_rates(k_mn47__he4_v43) - screened_rates(k_mn47__p_cr46) - &
      screened_rates(k_p_mn47__fe48)*Y(jp)*state % rho - &
      screened_rates(k_p_mn47__he4_cr44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn47, jmn47, scratch)

    scratch = (&
      screened_rates(k_fe47__mn47__weak__wc12) &
       )
    call set_jac_entry(state, jmn47, jfe47, scratch)

    scratch = (&
      screened_rates(k_fe48__p_mn47) &
       )
    call set_jac_entry(state, jmn47, jfe48, scratch)

    scratch = (&
      screened_rates(k_p_fe50__he4_mn47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn47, jfe50, scratch)

    scratch = (&
      screened_rates(k_co51__he4_mn47) &
       )
    call set_jac_entry(state, jmn47, jco51, scratch)

    scratch = (&
      screened_rates(k_p_cr47__mn48)*Y(jcr47)*state % rho + screened_rates(k_p_fe51__he4_mn48)* &
      Y(jfe51)*state % rho - screened_rates(k_p_mn48__fe49)*Y(jmn48)*state % rho - &
      screened_rates(k_p_mn48__he4_cr45)*Y(jmn48)*state % rho &
       )
    call set_jac_entry(state, jmn48, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr45__p_mn48)*Y(jcr45)*state % rho - screened_rates(k_he4_mn48__co52)* &
      Y(jmn48)*state % rho - screened_rates(k_he4_mn48__p_fe51)*Y(jmn48)* &
      state % rho + screened_rates(k_he4_v44__mn48)*Y(jv44)*state % rho &
       )
    call set_jac_entry(state, jmn48, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_v44__mn48)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn48, jv44, scratch)

    scratch = (&
      screened_rates(k_he4_cr45__p_mn48)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn48, jcr45, scratch)

    scratch = (&
      screened_rates(k_p_cr47__mn48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn48, jcr47, scratch)

    scratch = (&
      -screened_rates(k_he4_mn48__co52)*Y(jhe4)*state % rho - screened_rates(k_he4_mn48__p_fe51)* &
      Y(jhe4)*state % rho - screened_rates(k_mn48__cr48__weak__wc12) - &
      screened_rates(k_mn48__he4_v44) - screened_rates(k_mn48__p_cr47) - &
      screened_rates(k_mn48__p_v47__weak__wc12) - screened_rates(k_p_mn48__fe49)*Y(jp)* &
      state % rho - screened_rates(k_p_mn48__he4_cr45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn48, jmn48, scratch)

    scratch = (&
      screened_rates(k_fe48__mn48__weak__wc12) &
       )
    call set_jac_entry(state, jmn48, jfe48, scratch)

    scratch = (&
      screened_rates(k_fe49__p_mn48) &
       )
    call set_jac_entry(state, jmn48, jfe49, scratch)

    scratch = (&
      screened_rates(k_p_fe51__he4_mn48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn48, jfe51, scratch)

    scratch = (&
      screened_rates(k_co52__he4_mn48) &
       )
    call set_jac_entry(state, jmn48, jco52, scratch)

    scratch = (&
      screened_rates(k_p_cr48__mn49)*Y(jcr48)*state % rho + screened_rates(k_p_fe52__he4_mn49)* &
      Y(jfe52)*state % rho - screened_rates(k_p_mn49__fe50)*Y(jmn49)*state % rho - &
      screened_rates(k_p_mn49__he4_cr46)*Y(jmn49)*state % rho &
       )
    call set_jac_entry(state, jmn49, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr46__p_mn49)*Y(jcr46)*state % rho - screened_rates(k_he4_mn49__co53)* &
      Y(jmn49)*state % rho - screened_rates(k_he4_mn49__p_fe52)*Y(jmn49)* &
      state % rho + screened_rates(k_he4_v45__mn49)*Y(jv45)*state % rho &
       )
    call set_jac_entry(state, jmn49, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_v45__mn49)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn49, jv45, scratch)

    scratch = (&
      screened_rates(k_he4_cr46__p_mn49)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn49, jcr46, scratch)

    scratch = (&
      screened_rates(k_p_cr48__mn49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn49, jcr48, scratch)

    scratch = (&
      -screened_rates(k_he4_mn49__co53)*Y(jhe4)*state % rho - screened_rates(k_he4_mn49__p_fe52)* &
      Y(jhe4)*state % rho - screened_rates(k_mn49__cr49__weak__wc12) - &
      screened_rates(k_mn49__he4_v45) - screened_rates(k_mn49__p_cr48) - &
      screened_rates(k_p_mn49__fe50)*Y(jp)*state % rho - &
      screened_rates(k_p_mn49__he4_cr46)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn49, jmn49, scratch)

    scratch = (&
      screened_rates(k_fe49__mn49__weak__wc12) &
       )
    call set_jac_entry(state, jmn49, jfe49, scratch)

    scratch = (&
      screened_rates(k_fe50__p_mn49) &
       )
    call set_jac_entry(state, jmn49, jfe50, scratch)

    scratch = (&
      screened_rates(k_p_fe52__he4_mn49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn49, jfe52, scratch)

    scratch = (&
      screened_rates(k_co50__p_mn49__weak__wc12) &
       )
    call set_jac_entry(state, jmn49, jco50, scratch)

    scratch = (&
      screened_rates(k_co53__he4_mn49) &
       )
    call set_jac_entry(state, jmn49, jco53, scratch)

    scratch = (&
      screened_rates(k_p_cr49__mn50)*Y(jcr49)*state % rho + screened_rates(k_p_fe53__he4_mn50)* &
      Y(jfe53)*state % rho - screened_rates(k_p_mn50__fe51)*Y(jmn50)*state % rho - &
      screened_rates(k_p_mn50__he4_cr47)*Y(jmn50)*state % rho &
       )
    call set_jac_entry(state, jmn50, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr47__p_mn50)*Y(jcr47)*state % rho - screened_rates(k_he4_mn50__co54)* &
      Y(jmn50)*state % rho - screened_rates(k_he4_mn50__p_fe53)*Y(jmn50)* &
      state % rho + screened_rates(k_he4_v46__mn50)*Y(jv46)*state % rho &
       )
    call set_jac_entry(state, jmn50, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_v46__mn50)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn50, jv46, scratch)

    scratch = (&
      screened_rates(k_he4_cr47__p_mn50)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn50, jcr47, scratch)

    scratch = (&
      screened_rates(k_p_cr49__mn50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn50, jcr49, scratch)

    scratch = (&
      -screened_rates(k_he4_mn50__co54)*Y(jhe4)*state % rho - screened_rates(k_he4_mn50__p_fe53)* &
      Y(jhe4)*state % rho - screened_rates(k_mn50__cr50__weak__wc12) - &
      screened_rates(k_mn50__he4_v46) - screened_rates(k_mn50__p_cr49) - &
      screened_rates(k_p_mn50__fe51)*Y(jp)*state % rho - &
      screened_rates(k_p_mn50__he4_cr47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn50, jmn50, scratch)

    scratch = (&
      screened_rates(k_fe50__mn50__weak__wc12) &
       )
    call set_jac_entry(state, jmn50, jfe50, scratch)

    scratch = (&
      screened_rates(k_fe51__p_mn50) &
       )
    call set_jac_entry(state, jmn50, jfe51, scratch)

    scratch = (&
      screened_rates(k_p_fe53__he4_mn50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn50, jfe53, scratch)

    scratch = (&
      screened_rates(k_co54__he4_mn50) &
       )
    call set_jac_entry(state, jmn50, jco54, scratch)

    scratch = (&
      screened_rates(k_p_cr50__mn51)*Y(jcr50)*state % rho + screened_rates(k_p_fe54__he4_mn51)* &
      Y(jfe54)*state % rho - screened_rates(k_p_mn51__fe52)*Y(jmn51)*state % rho - &
      screened_rates(k_p_mn51__he4_cr48)*Y(jmn51)*state % rho &
       )
    call set_jac_entry(state, jmn51, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr48__p_mn51)*Y(jcr48)*state % rho - screened_rates(k_he4_mn51__co55)* &
      Y(jmn51)*state % rho - screened_rates(k_he4_mn51__p_fe54)*Y(jmn51)* &
      state % rho + screened_rates(k_he4_v47__mn51)*Y(jv47)*state % rho &
       )
    call set_jac_entry(state, jmn51, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_v47__mn51)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn51, jv47, scratch)

    scratch = (&
      screened_rates(k_he4_cr48__p_mn51)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn51, jcr48, scratch)

    scratch = (&
      screened_rates(k_p_cr50__mn51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn51, jcr50, scratch)

    scratch = (&
      -screened_rates(k_he4_mn51__co55)*Y(jhe4)*state % rho - screened_rates(k_he4_mn51__p_fe54)* &
      Y(jhe4)*state % rho - screened_rates(k_mn51__cr51__weak__wc12) - &
      screened_rates(k_mn51__he4_v47) - screened_rates(k_mn51__p_cr50) - &
      screened_rates(k_p_mn51__fe52)*Y(jp)*state % rho - &
      screened_rates(k_p_mn51__he4_cr48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn51, jmn51, scratch)

    scratch = (&
      screened_rates(k_fe51__mn51__weak__wc12) &
       )
    call set_jac_entry(state, jmn51, jfe51, scratch)

    scratch = (&
      screened_rates(k_fe52__p_mn51) &
       )
    call set_jac_entry(state, jmn51, jfe52, scratch)

    scratch = (&
      screened_rates(k_p_fe54__he4_mn51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn51, jfe54, scratch)

    scratch = (&
      screened_rates(k_co55__he4_mn51) &
       )
    call set_jac_entry(state, jmn51, jco55, scratch)

    scratch = (&
      screened_rates(k_p_cr51__mn52)*Y(jcr51)*state % rho + screened_rates(k_p_fe55__he4_mn52)* &
      Y(jfe55)*state % rho - screened_rates(k_p_mn52__fe53)*Y(jmn52)*state % rho - &
      screened_rates(k_p_mn52__he4_cr49)*Y(jmn52)*state % rho &
       )
    call set_jac_entry(state, jmn52, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr49__p_mn52)*Y(jcr49)*state % rho - screened_rates(k_he4_mn52__co56)* &
      Y(jmn52)*state % rho - screened_rates(k_he4_mn52__p_fe55)*Y(jmn52)* &
      state % rho + screened_rates(k_he4_v48__mn52)*Y(jv48)*state % rho &
       )
    call set_jac_entry(state, jmn52, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_v48__mn52)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn52, jv48, scratch)

    scratch = (&
      screened_rates(k_he4_cr49__p_mn52)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn52, jcr49, scratch)

    scratch = (&
      screened_rates(k_p_cr51__mn52)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn52, jcr51, scratch)

    scratch = (&
      -screened_rates(k_he4_mn52__co56)*Y(jhe4)*state % rho - screened_rates(k_he4_mn52__p_fe55)* &
      Y(jhe4)*state % rho - screened_rates(k_mn52__cr52__weak__wc12) - &
      screened_rates(k_mn52__he4_v48) - screened_rates(k_mn52__p_cr51) - &
      screened_rates(k_p_mn52__fe53)*Y(jp)*state % rho - &
      screened_rates(k_p_mn52__he4_cr49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn52, jmn52, scratch)

    scratch = (&
      screened_rates(k_fe52__mn52__weak__wc12) &
       )
    call set_jac_entry(state, jmn52, jfe52, scratch)

    scratch = (&
      screened_rates(k_fe53__p_mn52) &
       )
    call set_jac_entry(state, jmn52, jfe53, scratch)

    scratch = (&
      screened_rates(k_p_fe55__he4_mn52)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn52, jfe55, scratch)

    scratch = (&
      screened_rates(k_co56__he4_mn52) &
       )
    call set_jac_entry(state, jmn52, jco56, scratch)

    scratch = (&
      screened_rates(k_p_cr52__mn53)*Y(jcr52)*state % rho + screened_rates(k_p_fe56__he4_mn53)* &
      Y(jfe56)*state % rho - screened_rates(k_p_mn53__fe54)*Y(jmn53)*state % rho - &
      screened_rates(k_p_mn53__he4_cr50)*Y(jmn53)*state % rho &
       )
    call set_jac_entry(state, jmn53, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr50__p_mn53)*Y(jcr50)*state % rho - &
      screened_rates(k_he4_mn53__p_fe56)*Y(jmn53)*state % rho + &
      screened_rates(k_he4_v49__mn53)*Y(jv49)*state % rho &
       )
    call set_jac_entry(state, jmn53, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_v49__mn53)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn53, jv49, scratch)

    scratch = (&
      screened_rates(k_he4_cr50__p_mn53)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn53, jcr50, scratch)

    scratch = (&
      screened_rates(k_p_cr52__mn53)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn53, jcr52, scratch)

    scratch = (&
      -screened_rates(k_he4_mn53__p_fe56)*Y(jhe4)*state % rho - screened_rates(k_mn53__he4_v49) - &
      screened_rates(k_mn53__p_cr52) - screened_rates(k_p_mn53__fe54)*Y(jp)*state % rho - &
      screened_rates(k_p_mn53__he4_cr50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn53, jmn53, scratch)

    scratch = (&
      screened_rates(k_fe53__mn53__weak__wc12) &
       )
    call set_jac_entry(state, jmn53, jfe53, scratch)

    scratch = (&
      screened_rates(k_fe54__p_mn53) &
       )
    call set_jac_entry(state, jmn53, jfe54, scratch)

    scratch = (&
      screened_rates(k_p_fe56__he4_mn53)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn53, jfe56, scratch)

    scratch = (&
      -screened_rates(k_p_mn54__fe55)*Y(jmn54)*state % rho - screened_rates(k_p_mn54__he4_cr51)* &
      Y(jmn54)*state % rho &
       )
    call set_jac_entry(state, jmn54, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr51__p_mn54)*Y(jcr51)*state % rho + screened_rates(k_he4_v50__mn54)* &
      Y(jv50)*state % rho &
       )
    call set_jac_entry(state, jmn54, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_v50__mn54)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn54, jv50, scratch)

    scratch = (&
      screened_rates(k_he4_cr51__p_mn54)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn54, jcr51, scratch)

    scratch = (&
      -screened_rates(k_mn54__he4_v50) - screened_rates(k_p_mn54__fe55)*Y(jp)*state % rho - &
      screened_rates(k_p_mn54__he4_cr51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn54, jmn54, scratch)

    scratch = (&
      screened_rates(k_fe55__p_mn54) &
       )
    call set_jac_entry(state, jmn54, jfe55, scratch)

    scratch = (&
      -screened_rates(k_p_mn55__fe56)*Y(jmn55)*state % rho - screened_rates(k_p_mn55__he4_cr52)* &
      Y(jmn55)*state % rho &
       )
    call set_jac_entry(state, jmn55, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr52__p_mn55)*Y(jcr52)*state % rho &
       )
    call set_jac_entry(state, jmn55, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cr52__p_mn55)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jmn55, jcr52, scratch)

    scratch = (&
      -screened_rates(k_p_mn55__fe56)*Y(jp)*state % rho - screened_rates(k_p_mn55__he4_cr52)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(state, jmn55, jmn55, scratch)

    scratch = (&
      screened_rates(k_fe55__mn55__weak__wc12) &
       )
    call set_jac_entry(state, jmn55, jfe55, scratch)

    scratch = (&
      screened_rates(k_fe56__p_mn55) &
       )
    call set_jac_entry(state, jmn55, jfe56, scratch)

    scratch = (&
      screened_rates(k_p_co48__he4_fe45)*Y(jco48)*state % rho + screened_rates(k_p_mn44__fe45)* &
      Y(jmn44)*state % rho &
       )
    call set_jac_entry(state, jfe45, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_fe45__ni49)*Y(jfe45)*state % rho - screened_rates(k_he4_fe45__p_co48) &
      *Y(jfe45)*state % rho &
       )
    call set_jac_entry(state, jfe45, jhe4, scratch)

    scratch = (&
      screened_rates(k_p_mn44__fe45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe45, jmn44, scratch)

    scratch = (&
      -screened_rates(k_fe45__mn45__weak__wc17) - screened_rates(k_fe45__p_cr44__weak__wc17) - &
      screened_rates(k_fe45__p_mn44) - screened_rates(k_fe45__p_p_p_ti42__weak__wc12) - &
      screened_rates(k_fe45__p_p_v43__weak__wc12) - screened_rates(k_he4_fe45__ni49)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_fe45__p_co48)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe45, jfe45, scratch)

    scratch = (&
      screened_rates(k_p_co48__he4_fe45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe45, jco48, scratch)

    scratch = (&
      screened_rates(k_ni49__he4_fe45) &
       )
    call set_jac_entry(state, jfe45, jni49, scratch)

    scratch = (&
      screened_rates(k_p_co49__he4_fe46)*Y(jco49)*state % rho - screened_rates(k_p_fe46__co47)* &
      Y(jfe46)*state % rho + screened_rates(k_p_mn45__fe46)*Y(jmn45)*state % rho &
       )
    call set_jac_entry(state, jfe46, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr42__fe46)*Y(jcr42)*state % rho - screened_rates(k_he4_fe46__ni50)* &
      Y(jfe46)*state % rho - screened_rates(k_he4_fe46__p_co49)*Y(jfe46)* &
      state % rho &
       )
    call set_jac_entry(state, jfe46, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cr42__fe46)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe46, jcr42, scratch)

    scratch = (&
      screened_rates(k_p_mn45__fe46)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe46, jmn45, scratch)

    scratch = (&
      -screened_rates(k_fe46__he4_cr42) - screened_rates(k_fe46__mn46__weak__wc12) - &
      screened_rates(k_fe46__p_cr45__weak__wc12) - screened_rates(k_fe46__p_mn45) - &
      screened_rates(k_he4_fe46__ni50)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe46__p_co49)*Y(jhe4)*state % rho - &
      screened_rates(k_p_fe46__co47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe46, jfe46, scratch)

    scratch = (&
      screened_rates(k_co47__p_fe46) &
       )
    call set_jac_entry(state, jfe46, jco47, scratch)

    scratch = (&
      screened_rates(k_p_co49__he4_fe46)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe46, jco49, scratch)

    scratch = (&
      screened_rates(k_ni50__he4_fe46) &
       )
    call set_jac_entry(state, jfe46, jni50, scratch)

    scratch = (&
      screened_rates(k_p_co50__he4_fe47)*Y(jco50)*state % rho - screened_rates(k_p_fe47__co48)* &
      Y(jfe47)*state % rho - screened_rates(k_p_fe47__he4_mn44)*Y(jfe47)* &
      state % rho + screened_rates(k_p_mn46__fe47)*Y(jmn46)*state % rho &
       )
    call set_jac_entry(state, jfe47, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr43__fe47)*Y(jcr43)*state % rho - screened_rates(k_he4_fe47__ni51)* &
      Y(jfe47)*state % rho - screened_rates(k_he4_fe47__p_co50)*Y(jfe47)* &
      state % rho + screened_rates(k_he4_mn44__p_fe47)*Y(jmn44)*state % rho &
       )
    call set_jac_entry(state, jfe47, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cr43__fe47)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe47, jcr43, scratch)

    scratch = (&
      screened_rates(k_he4_mn44__p_fe47)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe47, jmn44, scratch)

    scratch = (&
      screened_rates(k_p_mn46__fe47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe47, jmn46, scratch)

    scratch = (&
      -screened_rates(k_fe47__he4_cr43) - screened_rates(k_fe47__mn47__weak__wc12) - &
      screened_rates(k_fe47__p_cr46__weak__wc12) - screened_rates(k_fe47__p_mn46) - &
      screened_rates(k_he4_fe47__ni51)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe47__p_co50)*Y(jhe4)*state % rho - &
      screened_rates(k_p_fe47__co48)*Y(jp)*state % rho - &
      screened_rates(k_p_fe47__he4_mn44)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe47, jfe47, scratch)

    scratch = (&
      screened_rates(k_co47__fe47__weak__bqa_pos_) &
       )
    call set_jac_entry(state, jfe47, jco47, scratch)

    scratch = (&
      screened_rates(k_co48__p_fe47) &
       )
    call set_jac_entry(state, jfe47, jco48, scratch)

    scratch = (&
      screened_rates(k_p_co50__he4_fe47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe47, jco50, scratch)

    scratch = (&
      screened_rates(k_ni51__he4_fe47) &
       )
    call set_jac_entry(state, jfe47, jni51, scratch)

    scratch = (&
      screened_rates(k_p_co51__he4_fe48)*Y(jco51)*state % rho - screened_rates(k_p_fe48__co49)* &
      Y(jfe48)*state % rho - screened_rates(k_p_fe48__he4_mn45)*Y(jfe48)* &
      state % rho + screened_rates(k_p_mn47__fe48)*Y(jmn47)*state % rho &
       )
    call set_jac_entry(state, jfe48, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr44__fe48)*Y(jcr44)*state % rho - screened_rates(k_he4_fe48__ni52)* &
      Y(jfe48)*state % rho - screened_rates(k_he4_fe48__p_co51)*Y(jfe48)* &
      state % rho + screened_rates(k_he4_mn45__p_fe48)*Y(jmn45)*state % rho &
       )
    call set_jac_entry(state, jfe48, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cr44__fe48)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe48, jcr44, scratch)

    scratch = (&
      screened_rates(k_he4_mn45__p_fe48)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe48, jmn45, scratch)

    scratch = (&
      screened_rates(k_p_mn47__fe48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe48, jmn47, scratch)

    scratch = (&
      -screened_rates(k_fe48__he4_cr44) - screened_rates(k_fe48__mn48__weak__wc12) - &
      screened_rates(k_fe48__p_cr47__weak__wc12) - screened_rates(k_fe48__p_mn47) - &
      screened_rates(k_he4_fe48__ni52)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe48__p_co51)*Y(jhe4)*state % rho - &
      screened_rates(k_p_fe48__co49)*Y(jp)*state % rho - &
      screened_rates(k_p_fe48__he4_mn45)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe48, jfe48, scratch)

    scratch = (&
      screened_rates(k_co48__fe48__weak__bqa_pos_) &
       )
    call set_jac_entry(state, jfe48, jco48, scratch)

    scratch = (&
      screened_rates(k_co49__p_fe48) &
       )
    call set_jac_entry(state, jfe48, jco49, scratch)

    scratch = (&
      screened_rates(k_p_co51__he4_fe48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe48, jco51, scratch)

    scratch = (&
      screened_rates(k_ni49__p_fe48__weak__wc12) &
       )
    call set_jac_entry(state, jfe48, jni49, scratch)

    scratch = (&
      screened_rates(k_ni52__he4_fe48) &
       )
    call set_jac_entry(state, jfe48, jni52, scratch)

    scratch = (&
      screened_rates(k_p_co52__he4_fe49)*Y(jco52)*state % rho - screened_rates(k_p_fe49__co50)* &
      Y(jfe49)*state % rho - screened_rates(k_p_fe49__he4_mn46)*Y(jfe49)* &
      state % rho + screened_rates(k_p_mn48__fe49)*Y(jmn48)*state % rho &
       )
    call set_jac_entry(state, jfe49, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr45__fe49)*Y(jcr45)*state % rho - screened_rates(k_he4_fe49__ni53)* &
      Y(jfe49)*state % rho - screened_rates(k_he4_fe49__p_co52)*Y(jfe49)* &
      state % rho + screened_rates(k_he4_mn46__p_fe49)*Y(jmn46)*state % rho &
       )
    call set_jac_entry(state, jfe49, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cr45__fe49)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe49, jcr45, scratch)

    scratch = (&
      screened_rates(k_he4_mn46__p_fe49)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe49, jmn46, scratch)

    scratch = (&
      screened_rates(k_p_mn48__fe49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe49, jmn48, scratch)

    scratch = (&
      -screened_rates(k_fe49__he4_cr45) - screened_rates(k_fe49__mn49__weak__wc12) - &
      screened_rates(k_fe49__p_cr48__weak__wc12) - screened_rates(k_fe49__p_mn48) - &
      screened_rates(k_he4_fe49__ni53)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe49__p_co52)*Y(jhe4)*state % rho - &
      screened_rates(k_p_fe49__co50)*Y(jp)*state % rho - &
      screened_rates(k_p_fe49__he4_mn46)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe49, jfe49, scratch)

    scratch = (&
      screened_rates(k_co49__fe49__weak__bqa_pos_) &
       )
    call set_jac_entry(state, jfe49, jco49, scratch)

    scratch = (&
      screened_rates(k_co50__p_fe49) &
       )
    call set_jac_entry(state, jfe49, jco50, scratch)

    scratch = (&
      screened_rates(k_p_co52__he4_fe49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe49, jco52, scratch)

    scratch = (&
      screened_rates(k_ni50__p_fe49__weak__wc12) &
       )
    call set_jac_entry(state, jfe49, jni50, scratch)

    scratch = (&
      screened_rates(k_ni53__he4_fe49) &
       )
    call set_jac_entry(state, jfe49, jni53, scratch)

    scratch = (&
      screened_rates(k_p_co53__he4_fe50)*Y(jco53)*state % rho - screened_rates(k_p_fe50__co51)* &
      Y(jfe50)*state % rho - screened_rates(k_p_fe50__he4_mn47)*Y(jfe50)* &
      state % rho + screened_rates(k_p_mn49__fe50)*Y(jmn49)*state % rho &
       )
    call set_jac_entry(state, jfe50, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr46__fe50)*Y(jcr46)*state % rho - screened_rates(k_he4_fe50__ni54)* &
      Y(jfe50)*state % rho - screened_rates(k_he4_fe50__p_co53)*Y(jfe50)* &
      state % rho + screened_rates(k_he4_mn47__p_fe50)*Y(jmn47)*state % rho &
       )
    call set_jac_entry(state, jfe50, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cr46__fe50)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe50, jcr46, scratch)

    scratch = (&
      screened_rates(k_he4_mn47__p_fe50)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe50, jmn47, scratch)

    scratch = (&
      screened_rates(k_p_mn49__fe50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe50, jmn49, scratch)

    scratch = (&
      -screened_rates(k_fe50__he4_cr46) - screened_rates(k_fe50__mn50__weak__wc12) - &
      screened_rates(k_fe50__p_mn49) - screened_rates(k_he4_fe50__ni54)*Y(jhe4)* &
      state % rho - screened_rates(k_he4_fe50__p_co53)*Y(jhe4)*state % rho - &
      screened_rates(k_p_fe50__co51)*Y(jp)*state % rho - &
      screened_rates(k_p_fe50__he4_mn47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe50, jfe50, scratch)

    scratch = (&
      screened_rates(k_co50__fe50__weak__wc12) &
       )
    call set_jac_entry(state, jfe50, jco50, scratch)

    scratch = (&
      screened_rates(k_co51__p_fe50) &
       )
    call set_jac_entry(state, jfe50, jco51, scratch)

    scratch = (&
      screened_rates(k_p_co53__he4_fe50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe50, jco53, scratch)

    scratch = (&
      screened_rates(k_ni51__p_fe50__weak__wc12) &
       )
    call set_jac_entry(state, jfe50, jni51, scratch)

    scratch = (&
      screened_rates(k_ni54__he4_fe50) &
       )
    call set_jac_entry(state, jfe50, jni54, scratch)

    scratch = (&
      screened_rates(k_p_co54__he4_fe51)*Y(jco54)*state % rho - screened_rates(k_p_fe51__co52)* &
      Y(jfe51)*state % rho - screened_rates(k_p_fe51__he4_mn48)*Y(jfe51)* &
      state % rho + screened_rates(k_p_mn50__fe51)*Y(jmn50)*state % rho &
       )
    call set_jac_entry(state, jfe51, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr47__fe51)*Y(jcr47)*state % rho - screened_rates(k_he4_fe51__ni55)* &
      Y(jfe51)*state % rho - screened_rates(k_he4_fe51__p_co54)*Y(jfe51)* &
      state % rho + screened_rates(k_he4_mn48__p_fe51)*Y(jmn48)*state % rho &
       )
    call set_jac_entry(state, jfe51, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cr47__fe51)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe51, jcr47, scratch)

    scratch = (&
      screened_rates(k_he4_mn48__p_fe51)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe51, jmn48, scratch)

    scratch = (&
      screened_rates(k_p_mn50__fe51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe51, jmn50, scratch)

    scratch = (&
      -screened_rates(k_fe51__he4_cr47) - screened_rates(k_fe51__mn51__weak__wc12) - &
      screened_rates(k_fe51__p_mn50) - screened_rates(k_he4_fe51__ni55)*Y(jhe4)* &
      state % rho - screened_rates(k_he4_fe51__p_co54)*Y(jhe4)*state % rho - &
      screened_rates(k_p_fe51__co52)*Y(jp)*state % rho - &
      screened_rates(k_p_fe51__he4_mn48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe51, jfe51, scratch)

    scratch = (&
      screened_rates(k_co51__fe51__weak__mo97) &
       )
    call set_jac_entry(state, jfe51, jco51, scratch)

    scratch = (&
      screened_rates(k_co52__p_fe51) &
       )
    call set_jac_entry(state, jfe51, jco52, scratch)

    scratch = (&
      screened_rates(k_p_co54__he4_fe51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe51, jco54, scratch)

    scratch = (&
      screened_rates(k_ni52__p_fe51__weak__wc12) &
       )
    call set_jac_entry(state, jfe51, jni52, scratch)

    scratch = (&
      screened_rates(k_ni55__he4_fe51) &
       )
    call set_jac_entry(state, jfe51, jni55, scratch)

    scratch = (&
      screened_rates(k_p_co55__he4_fe52)*Y(jco55)*state % rho - screened_rates(k_p_fe52__co53)* &
      Y(jfe52)*state % rho - screened_rates(k_p_fe52__he4_mn49)*Y(jfe52)* &
      state % rho + screened_rates(k_p_mn51__fe52)*Y(jmn51)*state % rho &
       )
    call set_jac_entry(state, jfe52, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr48__fe52)*Y(jcr48)*state % rho - screened_rates(k_he4_fe52__ni56)* &
      Y(jfe52)*state % rho - screened_rates(k_he4_fe52__p_co55)*Y(jfe52)* &
      state % rho + screened_rates(k_he4_mn49__p_fe52)*Y(jmn49)*state % rho &
       )
    call set_jac_entry(state, jfe52, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cr48__fe52)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe52, jcr48, scratch)

    scratch = (&
      screened_rates(k_he4_mn49__p_fe52)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe52, jmn49, scratch)

    scratch = (&
      screened_rates(k_p_mn51__fe52)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe52, jmn51, scratch)

    scratch = (&
      -screened_rates(k_fe52__he4_cr48) - screened_rates(k_fe52__mn52__weak__wc12) - &
      screened_rates(k_fe52__p_mn51) - screened_rates(k_he4_fe52__ni56)*Y(jhe4)* &
      state % rho - screened_rates(k_he4_fe52__p_co55)*Y(jhe4)*state % rho - &
      screened_rates(k_p_fe52__co53)*Y(jp)*state % rho - &
      screened_rates(k_p_fe52__he4_mn49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe52, jfe52, scratch)

    scratch = (&
      screened_rates(k_co52__fe52__weak__wc12) &
       )
    call set_jac_entry(state, jfe52, jco52, scratch)

    scratch = (&
      screened_rates(k_co53__p_fe52) &
       )
    call set_jac_entry(state, jfe52, jco53, scratch)

    scratch = (&
      screened_rates(k_p_co55__he4_fe52)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe52, jco55, scratch)

    scratch = (&
      screened_rates(k_ni53__p_fe52__weak__wc12) &
       )
    call set_jac_entry(state, jfe52, jni53, scratch)

    scratch = (&
      screened_rates(k_ni56__he4_fe52) &
       )
    call set_jac_entry(state, jfe52, jni56, scratch)

    scratch = (&
      screened_rates(k_p_co56__he4_fe53)*Y(jco56)*state % rho - screened_rates(k_p_fe53__co54)* &
      Y(jfe53)*state % rho - screened_rates(k_p_fe53__he4_mn50)*Y(jfe53)* &
      state % rho + screened_rates(k_p_mn52__fe53)*Y(jmn52)*state % rho &
       )
    call set_jac_entry(state, jfe53, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr49__fe53)*Y(jcr49)*state % rho - screened_rates(k_he4_fe53__p_co56)* &
      Y(jfe53)*state % rho + screened_rates(k_he4_mn50__p_fe53)*Y(jmn50)* &
      state % rho &
       )
    call set_jac_entry(state, jfe53, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cr49__fe53)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe53, jcr49, scratch)

    scratch = (&
      screened_rates(k_he4_mn50__p_fe53)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe53, jmn50, scratch)

    scratch = (&
      screened_rates(k_p_mn52__fe53)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe53, jmn52, scratch)

    scratch = (&
      -screened_rates(k_fe53__he4_cr49) - screened_rates(k_fe53__mn53__weak__wc12) - &
      screened_rates(k_fe53__p_mn52) - screened_rates(k_he4_fe53__p_co56)*Y(jhe4)* &
      state % rho - screened_rates(k_p_fe53__co54)*Y(jp)*state % rho - &
      screened_rates(k_p_fe53__he4_mn50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe53, jfe53, scratch)

    scratch = (&
      screened_rates(k_co53__fe53__weak__wc12) &
       )
    call set_jac_entry(state, jfe53, jco53, scratch)

    scratch = (&
      screened_rates(k_co54__p_fe53) &
       )
    call set_jac_entry(state, jfe53, jco54, scratch)

    scratch = (&
      screened_rates(k_p_co56__he4_fe53)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe53, jco56, scratch)

    scratch = (&
      -screened_rates(k_p_fe54__co55)*Y(jfe54)*state % rho - screened_rates(k_p_fe54__he4_mn51)* &
      Y(jfe54)*state % rho + screened_rates(k_p_mn53__fe54)*Y(jmn53)*state % rho &
       )
    call set_jac_entry(state, jfe54, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr50__fe54)*Y(jcr50)*state % rho + screened_rates(k_he4_mn51__p_fe54)* &
      Y(jmn51)*state % rho &
       )
    call set_jac_entry(state, jfe54, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cr50__fe54)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe54, jcr50, scratch)

    scratch = (&
      screened_rates(k_he4_mn51__p_fe54)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe54, jmn51, scratch)

    scratch = (&
      screened_rates(k_p_mn53__fe54)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe54, jmn53, scratch)

    scratch = (&
      -screened_rates(k_fe54__he4_cr50) - screened_rates(k_fe54__p_mn53) - screened_rates(k_p_fe54__co55) &
      *Y(jp)*state % rho - screened_rates(k_p_fe54__he4_mn51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe54, jfe54, scratch)

    scratch = (&
      screened_rates(k_co54__fe54__weak__wc12) &
       )
    call set_jac_entry(state, jfe54, jco54, scratch)

    scratch = (&
      screened_rates(k_co55__p_fe54) &
       )
    call set_jac_entry(state, jfe54, jco55, scratch)

    scratch = (&
      -screened_rates(k_p_fe55__co56)*Y(jfe55)*state % rho - screened_rates(k_p_fe55__he4_mn52)* &
      Y(jfe55)*state % rho + screened_rates(k_p_mn54__fe55)*Y(jmn54)*state % rho &
       )
    call set_jac_entry(state, jfe55, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr51__fe55)*Y(jcr51)*state % rho + screened_rates(k_he4_mn52__p_fe55)* &
      Y(jmn52)*state % rho &
       )
    call set_jac_entry(state, jfe55, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cr51__fe55)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe55, jcr51, scratch)

    scratch = (&
      screened_rates(k_he4_mn52__p_fe55)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe55, jmn52, scratch)

    scratch = (&
      screened_rates(k_p_mn54__fe55)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe55, jmn54, scratch)

    scratch = (&
      -screened_rates(k_fe55__he4_cr51) - screened_rates(k_fe55__mn55__weak__wc12) - &
      screened_rates(k_fe55__p_mn54) - screened_rates(k_p_fe55__co56)*Y(jp)*state % rho - &
      screened_rates(k_p_fe55__he4_mn52)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe55, jfe55, scratch)

    scratch = (&
      screened_rates(k_co55__fe55__weak__wc12) &
       )
    call set_jac_entry(state, jfe55, jco55, scratch)

    scratch = (&
      screened_rates(k_co56__p_fe55) &
       )
    call set_jac_entry(state, jfe55, jco56, scratch)

    scratch = (&
      -screened_rates(k_p_fe56__he4_mn53)*Y(jfe56)*state % rho + screened_rates(k_p_mn55__fe56)* &
      Y(jmn55)*state % rho &
       )
    call set_jac_entry(state, jfe56, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr52__fe56)*Y(jcr52)*state % rho + screened_rates(k_he4_mn53__p_fe56)* &
      Y(jmn53)*state % rho &
       )
    call set_jac_entry(state, jfe56, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cr52__fe56)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe56, jcr52, scratch)

    scratch = (&
      screened_rates(k_he4_mn53__p_fe56)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jfe56, jmn53, scratch)

    scratch = (&
      screened_rates(k_p_mn55__fe56)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe56, jmn55, scratch)

    scratch = (&
      -screened_rates(k_fe56__he4_cr52) - screened_rates(k_fe56__p_mn55) - &
      screened_rates(k_p_fe56__he4_mn53)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jfe56, jfe56, scratch)

    scratch = (&
      screened_rates(k_co56__fe56__weak__wc12) &
       )
    call set_jac_entry(state, jfe56, jco56, scratch)

    scratch = (&
      -screened_rates(k_p_co47__ni48)*Y(jco47)*state % rho + screened_rates(k_p_fe46__co47)* &
      Y(jfe46)*state % rho + screened_rates(k_p_ni50__he4_co47)*Y(jni50)* &
      state % rho &
       )
    call set_jac_entry(state, jco47, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_co47__p_ni50)*Y(jco47)*state % rho &
       )
    call set_jac_entry(state, jco47, jhe4, scratch)

    scratch = (&
      screened_rates(k_p_fe46__co47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco47, jfe46, scratch)

    scratch = (&
      -screened_rates(k_co47__fe47__weak__bqa_pos_) - screened_rates(k_co47__p_fe46) - &
      screened_rates(k_he4_co47__p_ni50)*Y(jhe4)*state % rho - &
      screened_rates(k_p_co47__ni48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco47, jco47, scratch)

    scratch = (&
      screened_rates(k_ni48__p_co47) &
       )
    call set_jac_entry(state, jco47, jni48, scratch)

    scratch = (&
      screened_rates(k_p_ni50__he4_co47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco47, jni50, scratch)

    scratch = (&
      -screened_rates(k_p_co48__he4_fe45)*Y(jco48)*state % rho - screened_rates(k_p_co48__ni49)* &
      Y(jco48)*state % rho + screened_rates(k_p_fe47__co48)*Y(jfe47)*state % rho + &
      screened_rates(k_p_ni51__he4_co48)*Y(jni51)*state % rho &
       )
    call set_jac_entry(state, jco48, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_co48__p_ni51)*Y(jco48)*state % rho + &
      screened_rates(k_he4_fe45__p_co48)*Y(jfe45)*state % rho + &
      screened_rates(k_he4_mn44__co48)*Y(jmn44)*state % rho &
       )
    call set_jac_entry(state, jco48, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_mn44__co48)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jco48, jmn44, scratch)

    scratch = (&
      screened_rates(k_he4_fe45__p_co48)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jco48, jfe45, scratch)

    scratch = (&
      screened_rates(k_p_fe47__co48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco48, jfe47, scratch)

    scratch = (&
      -screened_rates(k_co48__fe48__weak__bqa_pos_) - screened_rates(k_co48__he4_mn44) - &
      screened_rates(k_co48__p_fe47) - screened_rates(k_he4_co48__p_ni51)*Y(jhe4)* &
      state % rho - screened_rates(k_p_co48__he4_fe45)*Y(jp)*state % rho - &
      screened_rates(k_p_co48__ni49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco48, jco48, scratch)

    scratch = (&
      screened_rates(k_ni48__co48__weak__wc17) &
       )
    call set_jac_entry(state, jco48, jni48, scratch)

    scratch = (&
      screened_rates(k_ni49__p_co48) &
       )
    call set_jac_entry(state, jco48, jni49, scratch)

    scratch = (&
      screened_rates(k_p_ni51__he4_co48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco48, jni51, scratch)

    scratch = (&
      -screened_rates(k_p_co49__he4_fe46)*Y(jco49)*state % rho - screened_rates(k_p_co49__ni50)* &
      Y(jco49)*state % rho + screened_rates(k_p_fe48__co49)*Y(jfe48)*state % rho + &
      screened_rates(k_p_ni52__he4_co49)*Y(jni52)*state % rho &
       )
    call set_jac_entry(state, jco49, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_co49__p_ni52)*Y(jco49)*state % rho + &
      screened_rates(k_he4_fe46__p_co49)*Y(jfe46)*state % rho + &
      screened_rates(k_he4_mn45__co49)*Y(jmn45)*state % rho &
       )
    call set_jac_entry(state, jco49, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_mn45__co49)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jco49, jmn45, scratch)

    scratch = (&
      screened_rates(k_he4_fe46__p_co49)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jco49, jfe46, scratch)

    scratch = (&
      screened_rates(k_p_fe48__co49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco49, jfe48, scratch)

    scratch = (&
      -screened_rates(k_co49__fe49__weak__bqa_pos_) - screened_rates(k_co49__he4_mn45) - &
      screened_rates(k_co49__p_fe48) - screened_rates(k_he4_co49__p_ni52)*Y(jhe4)* &
      state % rho - screened_rates(k_p_co49__he4_fe46)*Y(jp)*state % rho - &
      screened_rates(k_p_co49__ni50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco49, jco49, scratch)

    scratch = (&
      screened_rates(k_ni49__co49__weak__wc12) &
       )
    call set_jac_entry(state, jco49, jni49, scratch)

    scratch = (&
      screened_rates(k_ni50__p_co49) &
       )
    call set_jac_entry(state, jco49, jni50, scratch)

    scratch = (&
      screened_rates(k_p_ni52__he4_co49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco49, jni52, scratch)

    scratch = (&
      -screened_rates(k_p_co50__he4_fe47)*Y(jco50)*state % rho - screened_rates(k_p_co50__ni51)* &
      Y(jco50)*state % rho + screened_rates(k_p_fe49__co50)*Y(jfe49)*state % rho + &
      screened_rates(k_p_ni53__he4_co50)*Y(jni53)*state % rho &
       )
    call set_jac_entry(state, jco50, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_co50__p_ni53)*Y(jco50)*state % rho + &
      screened_rates(k_he4_fe47__p_co50)*Y(jfe47)*state % rho + &
      screened_rates(k_he4_mn46__co50)*Y(jmn46)*state % rho &
       )
    call set_jac_entry(state, jco50, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_mn46__co50)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jco50, jmn46, scratch)

    scratch = (&
      screened_rates(k_he4_fe47__p_co50)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jco50, jfe47, scratch)

    scratch = (&
      screened_rates(k_p_fe49__co50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco50, jfe49, scratch)

    scratch = (&
      -screened_rates(k_co50__fe50__weak__wc12) - screened_rates(k_co50__he4_mn46) - &
      screened_rates(k_co50__p_fe49) - screened_rates(k_co50__p_mn49__weak__wc12) - &
      screened_rates(k_he4_co50__p_ni53)*Y(jhe4)*state % rho - &
      screened_rates(k_p_co50__he4_fe47)*Y(jp)*state % rho - &
      screened_rates(k_p_co50__ni51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco50, jco50, scratch)

    scratch = (&
      screened_rates(k_ni50__co50__weak__wc12) &
       )
    call set_jac_entry(state, jco50, jni50, scratch)

    scratch = (&
      screened_rates(k_ni51__p_co50) &
       )
    call set_jac_entry(state, jco50, jni51, scratch)

    scratch = (&
      screened_rates(k_p_ni53__he4_co50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco50, jni53, scratch)

    scratch = (&
      -screened_rates(k_p_co51__he4_fe48)*Y(jco51)*state % rho - screened_rates(k_p_co51__ni52)* &
      Y(jco51)*state % rho + screened_rates(k_p_fe50__co51)*Y(jfe50)*state % rho + &
      screened_rates(k_p_ni54__he4_co51)*Y(jni54)*state % rho &
       )
    call set_jac_entry(state, jco51, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_co51__p_ni54)*Y(jco51)*state % rho + &
      screened_rates(k_he4_fe48__p_co51)*Y(jfe48)*state % rho + &
      screened_rates(k_he4_mn47__co51)*Y(jmn47)*state % rho &
       )
    call set_jac_entry(state, jco51, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_mn47__co51)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jco51, jmn47, scratch)

    scratch = (&
      screened_rates(k_he4_fe48__p_co51)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jco51, jfe48, scratch)

    scratch = (&
      screened_rates(k_p_fe50__co51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco51, jfe50, scratch)

    scratch = (&
      -screened_rates(k_co51__fe51__weak__mo97) - screened_rates(k_co51__he4_mn47) - &
      screened_rates(k_co51__p_fe50) - screened_rates(k_he4_co51__p_ni54)*Y(jhe4)* &
      state % rho - screened_rates(k_p_co51__he4_fe48)*Y(jp)*state % rho - &
      screened_rates(k_p_co51__ni52)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco51, jco51, scratch)

    scratch = (&
      screened_rates(k_ni51__co51__weak__wc17) &
       )
    call set_jac_entry(state, jco51, jni51, scratch)

    scratch = (&
      screened_rates(k_ni52__p_co51) &
       )
    call set_jac_entry(state, jco51, jni52, scratch)

    scratch = (&
      screened_rates(k_p_ni54__he4_co51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco51, jni54, scratch)

    scratch = (&
      -screened_rates(k_p_co52__he4_fe49)*Y(jco52)*state % rho - screened_rates(k_p_co52__ni53)* &
      Y(jco52)*state % rho + screened_rates(k_p_fe51__co52)*Y(jfe51)*state % rho + &
      screened_rates(k_p_ni55__he4_co52)*Y(jni55)*state % rho &
       )
    call set_jac_entry(state, jco52, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_co52__p_ni55)*Y(jco52)*state % rho + &
      screened_rates(k_he4_fe49__p_co52)*Y(jfe49)*state % rho + &
      screened_rates(k_he4_mn48__co52)*Y(jmn48)*state % rho &
       )
    call set_jac_entry(state, jco52, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_mn48__co52)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jco52, jmn48, scratch)

    scratch = (&
      screened_rates(k_he4_fe49__p_co52)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jco52, jfe49, scratch)

    scratch = (&
      screened_rates(k_p_fe51__co52)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco52, jfe51, scratch)

    scratch = (&
      -screened_rates(k_co52__fe52__weak__wc12) - screened_rates(k_co52__he4_mn48) - &
      screened_rates(k_co52__p_fe51) - screened_rates(k_he4_co52__p_ni55)*Y(jhe4)* &
      state % rho - screened_rates(k_p_co52__he4_fe49)*Y(jp)*state % rho - &
      screened_rates(k_p_co52__ni53)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco52, jco52, scratch)

    scratch = (&
      screened_rates(k_ni52__co52__weak__wc12) &
       )
    call set_jac_entry(state, jco52, jni52, scratch)

    scratch = (&
      screened_rates(k_ni53__p_co52) &
       )
    call set_jac_entry(state, jco52, jni53, scratch)

    scratch = (&
      screened_rates(k_p_ni55__he4_co52)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco52, jni55, scratch)

    scratch = (&
      -screened_rates(k_p_co53__he4_fe50)*Y(jco53)*state % rho - screened_rates(k_p_co53__ni54)* &
      Y(jco53)*state % rho + screened_rates(k_p_fe52__co53)*Y(jfe52)*state % rho + &
      screened_rates(k_p_ni56__he4_co53)*Y(jni56)*state % rho &
       )
    call set_jac_entry(state, jco53, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_co53__p_ni56)*Y(jco53)*state % rho + &
      screened_rates(k_he4_fe50__p_co53)*Y(jfe50)*state % rho + &
      screened_rates(k_he4_mn49__co53)*Y(jmn49)*state % rho &
       )
    call set_jac_entry(state, jco53, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_mn49__co53)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jco53, jmn49, scratch)

    scratch = (&
      screened_rates(k_he4_fe50__p_co53)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jco53, jfe50, scratch)

    scratch = (&
      screened_rates(k_p_fe52__co53)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco53, jfe52, scratch)

    scratch = (&
      -screened_rates(k_co53__fe53__weak__wc12) - screened_rates(k_co53__he4_mn49) - &
      screened_rates(k_co53__p_fe52) - screened_rates(k_he4_co53__p_ni56)*Y(jhe4)* &
      state % rho - screened_rates(k_p_co53__he4_fe50)*Y(jp)*state % rho - &
      screened_rates(k_p_co53__ni54)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco53, jco53, scratch)

    scratch = (&
      screened_rates(k_ni53__co53__weak__wc12) &
       )
    call set_jac_entry(state, jco53, jni53, scratch)

    scratch = (&
      screened_rates(k_ni54__p_co53) &
       )
    call set_jac_entry(state, jco53, jni54, scratch)

    scratch = (&
      screened_rates(k_p_ni56__he4_co53)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco53, jni56, scratch)

    scratch = (&
      -screened_rates(k_p_co54__he4_fe51)*Y(jco54)*state % rho - screened_rates(k_p_co54__ni55)* &
      Y(jco54)*state % rho + screened_rates(k_p_fe53__co54)*Y(jfe53)*state % rho &
       )
    call set_jac_entry(state, jco54, jp, scratch)

    scratch = (&
      screened_rates(k_he4_fe51__p_co54)*Y(jfe51)*state % rho + screened_rates(k_he4_mn50__co54)* &
      Y(jmn50)*state % rho &
       )
    call set_jac_entry(state, jco54, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_mn50__co54)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jco54, jmn50, scratch)

    scratch = (&
      screened_rates(k_he4_fe51__p_co54)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jco54, jfe51, scratch)

    scratch = (&
      screened_rates(k_p_fe53__co54)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco54, jfe53, scratch)

    scratch = (&
      -screened_rates(k_co54__fe54__weak__wc12) - screened_rates(k_co54__he4_mn50) - &
      screened_rates(k_co54__p_fe53) - screened_rates(k_p_co54__he4_fe51)*Y(jp)* &
      state % rho - screened_rates(k_p_co54__ni55)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco54, jco54, scratch)

    scratch = (&
      screened_rates(k_ni54__co54__weak__wc12) &
       )
    call set_jac_entry(state, jco54, jni54, scratch)

    scratch = (&
      screened_rates(k_ni55__p_co54) &
       )
    call set_jac_entry(state, jco54, jni55, scratch)

    scratch = (&
      -screened_rates(k_p_co55__he4_fe52)*Y(jco55)*state % rho - screened_rates(k_p_co55__ni56)* &
      Y(jco55)*state % rho + screened_rates(k_p_fe54__co55)*Y(jfe54)*state % rho &
       )
    call set_jac_entry(state, jco55, jp, scratch)

    scratch = (&
      screened_rates(k_he4_fe52__p_co55)*Y(jfe52)*state % rho + screened_rates(k_he4_mn51__co55)* &
      Y(jmn51)*state % rho &
       )
    call set_jac_entry(state, jco55, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_mn51__co55)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jco55, jmn51, scratch)

    scratch = (&
      screened_rates(k_he4_fe52__p_co55)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jco55, jfe52, scratch)

    scratch = (&
      screened_rates(k_p_fe54__co55)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco55, jfe54, scratch)

    scratch = (&
      -screened_rates(k_co55__fe55__weak__wc12) - screened_rates(k_co55__he4_mn51) - &
      screened_rates(k_co55__p_fe54) - screened_rates(k_p_co55__he4_fe52)*Y(jp)* &
      state % rho - screened_rates(k_p_co55__ni56)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco55, jco55, scratch)

    scratch = (&
      screened_rates(k_ni55__co55__weak__wc12) &
       )
    call set_jac_entry(state, jco55, jni55, scratch)

    scratch = (&
      screened_rates(k_ni56__p_co55) &
       )
    call set_jac_entry(state, jco55, jni56, scratch)

    scratch = (&
      -screened_rates(k_p_co56__he4_fe53)*Y(jco56)*state % rho + screened_rates(k_p_fe55__co56)* &
      Y(jfe55)*state % rho &
       )
    call set_jac_entry(state, jco56, jp, scratch)

    scratch = (&
      screened_rates(k_he4_fe53__p_co56)*Y(jfe53)*state % rho + screened_rates(k_he4_mn52__co56)* &
      Y(jmn52)*state % rho &
       )
    call set_jac_entry(state, jco56, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_mn52__co56)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jco56, jmn52, scratch)

    scratch = (&
      screened_rates(k_he4_fe53__p_co56)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jco56, jfe53, scratch)

    scratch = (&
      screened_rates(k_p_fe55__co56)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jco56, jfe55, scratch)

    scratch = (&
      -screened_rates(k_co56__fe56__weak__wc12) - screened_rates(k_co56__he4_mn52) - &
      screened_rates(k_co56__p_fe55) - screened_rates(k_p_co56__he4_fe53)*Y(jp)* &
      state % rho &
       )
    call set_jac_entry(state, jco56, jco56, scratch)

    scratch = (&
      screened_rates(k_ni56__co56__weak__wc12) &
       )
    call set_jac_entry(state, jco56, jni56, scratch)

    scratch = (&
      screened_rates(k_p_co47__ni48)*Y(jco47)*state % rho &
       )
    call set_jac_entry(state, jni48, jp, scratch)

    scratch = (&
      screened_rates(k_p_co47__ni48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jni48, jco47, scratch)

    scratch = (&
      -screened_rates(k_ni48__co48__weak__wc17) - screened_rates(k_ni48__p_co47) &
       )
    call set_jac_entry(state, jni48, jni48, scratch)

    scratch = (&
      screened_rates(k_p_co48__ni49)*Y(jco48)*state % rho &
       )
    call set_jac_entry(state, jni49, jp, scratch)

    scratch = (&
      screened_rates(k_he4_fe45__ni49)*Y(jfe45)*state % rho &
       )
    call set_jac_entry(state, jni49, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_fe45__ni49)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jni49, jfe45, scratch)

    scratch = (&
      screened_rates(k_p_co48__ni49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jni49, jco48, scratch)

    scratch = (&
      -screened_rates(k_ni49__co49__weak__wc12) - screened_rates(k_ni49__he4_fe45) - &
      screened_rates(k_ni49__p_co48) - screened_rates(k_ni49__p_fe48__weak__wc12) &
       )
    call set_jac_entry(state, jni49, jni49, scratch)

    scratch = (&
      screened_rates(k_p_co49__ni50)*Y(jco49)*state % rho - screened_rates(k_p_ni50__he4_co47)* &
      Y(jni50)*state % rho &
       )
    call set_jac_entry(state, jni50, jp, scratch)

    scratch = (&
      screened_rates(k_he4_co47__p_ni50)*Y(jco47)*state % rho + screened_rates(k_he4_fe46__ni50)* &
      Y(jfe46)*state % rho &
       )
    call set_jac_entry(state, jni50, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_fe46__ni50)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jni50, jfe46, scratch)

    scratch = (&
      screened_rates(k_he4_co47__p_ni50)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jni50, jco47, scratch)

    scratch = (&
      screened_rates(k_p_co49__ni50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jni50, jco49, scratch)

    scratch = (&
      -screened_rates(k_ni50__co50__weak__wc12) - screened_rates(k_ni50__he4_fe46) - &
      screened_rates(k_ni50__p_co49) - screened_rates(k_ni50__p_fe49__weak__wc12) - &
      screened_rates(k_p_ni50__he4_co47)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jni50, jni50, scratch)

    scratch = (&
      screened_rates(k_p_co50__ni51)*Y(jco50)*state % rho - screened_rates(k_p_ni51__he4_co48)* &
      Y(jni51)*state % rho &
       )
    call set_jac_entry(state, jni51, jp, scratch)

    scratch = (&
      screened_rates(k_he4_co48__p_ni51)*Y(jco48)*state % rho + screened_rates(k_he4_fe47__ni51)* &
      Y(jfe47)*state % rho &
       )
    call set_jac_entry(state, jni51, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_fe47__ni51)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jni51, jfe47, scratch)

    scratch = (&
      screened_rates(k_he4_co48__p_ni51)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jni51, jco48, scratch)

    scratch = (&
      screened_rates(k_p_co50__ni51)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jni51, jco50, scratch)

    scratch = (&
      -screened_rates(k_ni51__co51__weak__wc17) - screened_rates(k_ni51__he4_fe47) - &
      screened_rates(k_ni51__p_co50) - screened_rates(k_ni51__p_fe50__weak__wc12) - &
      screened_rates(k_p_ni51__he4_co48)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jni51, jni51, scratch)

    scratch = (&
      screened_rates(k_p_co51__ni52)*Y(jco51)*state % rho - screened_rates(k_p_ni52__he4_co49)* &
      Y(jni52)*state % rho &
       )
    call set_jac_entry(state, jni52, jp, scratch)

    scratch = (&
      screened_rates(k_he4_co49__p_ni52)*Y(jco49)*state % rho + screened_rates(k_he4_fe48__ni52)* &
      Y(jfe48)*state % rho &
       )
    call set_jac_entry(state, jni52, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_fe48__ni52)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jni52, jfe48, scratch)

    scratch = (&
      screened_rates(k_he4_co49__p_ni52)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jni52, jco49, scratch)

    scratch = (&
      screened_rates(k_p_co51__ni52)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jni52, jco51, scratch)

    scratch = (&
      -screened_rates(k_ni52__co52__weak__wc12) - screened_rates(k_ni52__he4_fe48) - &
      screened_rates(k_ni52__p_co51) - screened_rates(k_ni52__p_fe51__weak__wc12) - &
      screened_rates(k_p_ni52__he4_co49)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jni52, jni52, scratch)

    scratch = (&
      screened_rates(k_p_co52__ni53)*Y(jco52)*state % rho - screened_rates(k_p_ni53__he4_co50)* &
      Y(jni53)*state % rho &
       )
    call set_jac_entry(state, jni53, jp, scratch)

    scratch = (&
      screened_rates(k_he4_co50__p_ni53)*Y(jco50)*state % rho + screened_rates(k_he4_fe49__ni53)* &
      Y(jfe49)*state % rho &
       )
    call set_jac_entry(state, jni53, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_fe49__ni53)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jni53, jfe49, scratch)

    scratch = (&
      screened_rates(k_he4_co50__p_ni53)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jni53, jco50, scratch)

    scratch = (&
      screened_rates(k_p_co52__ni53)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jni53, jco52, scratch)

    scratch = (&
      -screened_rates(k_ni53__co53__weak__wc12) - screened_rates(k_ni53__he4_fe49) - &
      screened_rates(k_ni53__p_co52) - screened_rates(k_ni53__p_fe52__weak__wc12) - &
      screened_rates(k_p_ni53__he4_co50)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jni53, jni53, scratch)

    scratch = (&
      screened_rates(k_p_co53__ni54)*Y(jco53)*state % rho - screened_rates(k_p_ni54__he4_co51)* &
      Y(jni54)*state % rho &
       )
    call set_jac_entry(state, jni54, jp, scratch)

    scratch = (&
      screened_rates(k_he4_co51__p_ni54)*Y(jco51)*state % rho + screened_rates(k_he4_fe50__ni54)* &
      Y(jfe50)*state % rho &
       )
    call set_jac_entry(state, jni54, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_fe50__ni54)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jni54, jfe50, scratch)

    scratch = (&
      screened_rates(k_he4_co51__p_ni54)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jni54, jco51, scratch)

    scratch = (&
      screened_rates(k_p_co53__ni54)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jni54, jco53, scratch)

    scratch = (&
      -screened_rates(k_ni54__co54__weak__wc12) - screened_rates(k_ni54__he4_fe50) - &
      screened_rates(k_ni54__p_co53) - screened_rates(k_p_ni54__he4_co51)*Y(jp)* &
      state % rho &
       )
    call set_jac_entry(state, jni54, jni54, scratch)

    scratch = (&
      screened_rates(k_p_co54__ni55)*Y(jco54)*state % rho - screened_rates(k_p_ni55__he4_co52)* &
      Y(jni55)*state % rho &
       )
    call set_jac_entry(state, jni55, jp, scratch)

    scratch = (&
      screened_rates(k_he4_co52__p_ni55)*Y(jco52)*state % rho + screened_rates(k_he4_fe51__ni55)* &
      Y(jfe51)*state % rho &
       )
    call set_jac_entry(state, jni55, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_fe51__ni55)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jni55, jfe51, scratch)

    scratch = (&
      screened_rates(k_he4_co52__p_ni55)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jni55, jco52, scratch)

    scratch = (&
      screened_rates(k_p_co54__ni55)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jni55, jco54, scratch)

    scratch = (&
      -screened_rates(k_ni55__co55__weak__wc12) - screened_rates(k_ni55__he4_fe51) - &
      screened_rates(k_ni55__p_co54) - screened_rates(k_p_ni55__he4_co52)*Y(jp)* &
      state % rho &
       )
    call set_jac_entry(state, jni55, jni55, scratch)

    scratch = (&
      screened_rates(k_p_co55__ni56)*Y(jco55)*state % rho - screened_rates(k_p_ni56__he4_co53)* &
      Y(jni56)*state % rho &
       )
    call set_jac_entry(state, jni56, jp, scratch)

    scratch = (&
      screened_rates(k_he4_co53__p_ni56)*Y(jco53)*state % rho + screened_rates(k_he4_fe52__ni56)* &
      Y(jfe52)*state % rho &
       )
    call set_jac_entry(state, jni56, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_fe52__ni56)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jni56, jfe52, scratch)

    scratch = (&
      screened_rates(k_he4_co53__p_ni56)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jni56, jco53, scratch)

    scratch = (&
      screened_rates(k_p_co55__ni56)*Y(jp)*state % rho &
       )
    call set_jac_entry(state, jni56, jco55, scratch)

    scratch = (&
      -screened_rates(k_ni56__co56__weak__wc12) - screened_rates(k_ni56__he4_fe52) - &
      screened_rates(k_ni56__p_co55) - screened_rates(k_p_ni56__he4_co53)*Y(jp)* &
      state % rho &
       )
    call set_jac_entry(state, jni56, jni56, scratch)


  end subroutine jac_nuc

end module actual_rhs_module
