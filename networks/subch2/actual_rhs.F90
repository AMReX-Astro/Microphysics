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
      rate_eval % unscreened_rates(i_scor,33) = scor
      rate_eval % unscreened_rates(i_dscor_dt,33) = dscor_dt


      call screen5(pstate, 2, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,34) = scor
      rate_eval % unscreened_rates(i_dscor_dt,34) = dscor_dt


      call screen5(pstate, 3, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,35) = scor
      rate_eval % unscreened_rates(i_dscor_dt,35) = dscor_dt


      call screen5(pstate, 4, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,36) = scor
      rate_eval % unscreened_rates(i_dscor_dt,36) = dscor_dt


      call screen5(pstate, 5, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,37) = scor
      rate_eval % unscreened_rates(i_dscor_dt,37) = dscor_dt


      call screen5(pstate, 6, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,38) = scor
      rate_eval % unscreened_rates(i_dscor_dt,38) = dscor_dt
      rate_eval % unscreened_rates(i_scor,70) = scor
      rate_eval % unscreened_rates(i_dscor_dt,70) = dscor_dt


      call screen5(pstate, 7, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,39) = scor
      rate_eval % unscreened_rates(i_dscor_dt,39) = dscor_dt
      rate_eval % unscreened_rates(i_scor,74) = scor
      rate_eval % unscreened_rates(i_dscor_dt,74) = dscor_dt
      rate_eval % unscreened_rates(i_scor,75) = scor
      rate_eval % unscreened_rates(i_dscor_dt,75) = dscor_dt


      call screen5(pstate, 8, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,40) = scor
      rate_eval % unscreened_rates(i_dscor_dt,40) = dscor_dt
      rate_eval % unscreened_rates(i_scor,76) = scor
      rate_eval % unscreened_rates(i_dscor_dt,76) = dscor_dt
      rate_eval % unscreened_rates(i_scor,77) = scor
      rate_eval % unscreened_rates(i_dscor_dt,77) = dscor_dt


      call screen5(pstate, 9, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,41) = scor
      rate_eval % unscreened_rates(i_dscor_dt,41) = dscor_dt


      call screen5(pstate, 10, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,42) = scor
      rate_eval % unscreened_rates(i_dscor_dt,42) = dscor_dt
      rate_eval % unscreened_rates(i_scor,78) = scor
      rate_eval % unscreened_rates(i_dscor_dt,78) = dscor_dt
      rate_eval % unscreened_rates(i_scor,79) = scor
      rate_eval % unscreened_rates(i_dscor_dt,79) = dscor_dt
      rate_eval % unscreened_rates(i_scor,80) = scor
      rate_eval % unscreened_rates(i_dscor_dt,80) = dscor_dt


      call screen5(pstate, 11, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,43) = scor
      rate_eval % unscreened_rates(i_dscor_dt,43) = dscor_dt
      rate_eval % unscreened_rates(i_scor,81) = scor
      rate_eval % unscreened_rates(i_dscor_dt,81) = dscor_dt
      rate_eval % unscreened_rates(i_scor,82) = scor
      rate_eval % unscreened_rates(i_dscor_dt,82) = dscor_dt
      rate_eval % unscreened_rates(i_scor,83) = scor
      rate_eval % unscreened_rates(i_dscor_dt,83) = dscor_dt


      call screen5(pstate, 12, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,44) = scor
      rate_eval % unscreened_rates(i_dscor_dt,44) = dscor_dt


      call screen5(pstate, 13, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,45) = scor
      rate_eval % unscreened_rates(i_dscor_dt,45) = dscor_dt
      rate_eval % unscreened_rates(i_scor,84) = scor
      rate_eval % unscreened_rates(i_dscor_dt,84) = dscor_dt


      call screen5(pstate, 14, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,46) = scor
      rate_eval % unscreened_rates(i_dscor_dt,46) = dscor_dt
      rate_eval % unscreened_rates(i_scor,85) = scor
      rate_eval % unscreened_rates(i_dscor_dt,85) = dscor_dt


      call screen5(pstate, 15, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,47) = scor
      rate_eval % unscreened_rates(i_dscor_dt,47) = dscor_dt


      call screen5(pstate, 16, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,48) = scor
      rate_eval % unscreened_rates(i_dscor_dt,48) = dscor_dt
      rate_eval % unscreened_rates(i_scor,86) = scor
      rate_eval % unscreened_rates(i_dscor_dt,86) = dscor_dt


      call screen5(pstate, 17, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,49) = scor
      rate_eval % unscreened_rates(i_dscor_dt,49) = dscor_dt
      rate_eval % unscreened_rates(i_scor,87) = scor
      rate_eval % unscreened_rates(i_dscor_dt,87) = dscor_dt


      call screen5(pstate, 18, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,50) = scor
      rate_eval % unscreened_rates(i_dscor_dt,50) = dscor_dt


      call screen5(pstate, 19, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,51) = scor
      rate_eval % unscreened_rates(i_dscor_dt,51) = dscor_dt
      rate_eval % unscreened_rates(i_scor,88) = scor
      rate_eval % unscreened_rates(i_dscor_dt,88) = dscor_dt


      call screen5(pstate, 20, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,52) = scor
      rate_eval % unscreened_rates(i_dscor_dt,52) = dscor_dt
      rate_eval % unscreened_rates(i_scor,89) = scor
      rate_eval % unscreened_rates(i_dscor_dt,89) = dscor_dt


      call screen5(pstate, 21, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,53) = scor
      rate_eval % unscreened_rates(i_dscor_dt,53) = dscor_dt


      call screen5(pstate, 22, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,54) = scor
      rate_eval % unscreened_rates(i_dscor_dt,54) = dscor_dt
      rate_eval % unscreened_rates(i_scor,90) = scor
      rate_eval % unscreened_rates(i_dscor_dt,90) = dscor_dt


      call screen5(pstate, 23, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,55) = scor
      rate_eval % unscreened_rates(i_dscor_dt,55) = dscor_dt
      rate_eval % unscreened_rates(i_scor,91) = scor
      rate_eval % unscreened_rates(i_dscor_dt,91) = dscor_dt


      call screen5(pstate, 24, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,56) = scor
      rate_eval % unscreened_rates(i_dscor_dt,56) = dscor_dt


      call screen5(pstate, 25, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,57) = scor
      rate_eval % unscreened_rates(i_dscor_dt,57) = dscor_dt
      rate_eval % unscreened_rates(i_scor,92) = scor
      rate_eval % unscreened_rates(i_dscor_dt,92) = dscor_dt


      call screen5(pstate, 26, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,58) = scor
      rate_eval % unscreened_rates(i_dscor_dt,58) = dscor_dt
      rate_eval % unscreened_rates(i_scor,93) = scor
      rate_eval % unscreened_rates(i_dscor_dt,93) = dscor_dt


      call screen5(pstate, 27, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,59) = scor
      rate_eval % unscreened_rates(i_dscor_dt,59) = dscor_dt


      call screen5(pstate, 28, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,60) = scor
      rate_eval % unscreened_rates(i_dscor_dt,60) = dscor_dt
      rate_eval % unscreened_rates(i_scor,94) = scor
      rate_eval % unscreened_rates(i_dscor_dt,94) = dscor_dt


      call screen5(pstate, 29, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,61) = scor
      rate_eval % unscreened_rates(i_dscor_dt,61) = dscor_dt
      rate_eval % unscreened_rates(i_scor,95) = scor
      rate_eval % unscreened_rates(i_dscor_dt,95) = dscor_dt


      call screen5(pstate, 30, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,62) = scor
      rate_eval % unscreened_rates(i_dscor_dt,62) = dscor_dt


      call screen5(pstate, 31, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,63) = scor
      rate_eval % unscreened_rates(i_dscor_dt,63) = dscor_dt


      call screen5(pstate, 32, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,64) = scor
      rate_eval % unscreened_rates(i_dscor_dt,64) = dscor_dt


      call screen5(pstate, 33, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,65) = scor
      rate_eval % unscreened_rates(i_dscor_dt,65) = dscor_dt
      rate_eval % unscreened_rates(i_scor,66) = scor
      rate_eval % unscreened_rates(i_dscor_dt,66) = dscor_dt


      call screen5(pstate, 34, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,67) = scor
      rate_eval % unscreened_rates(i_dscor_dt,67) = dscor_dt
      rate_eval % unscreened_rates(i_scor,68) = scor
      rate_eval % unscreened_rates(i_dscor_dt,68) = dscor_dt


      call screen5(pstate, 35, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,69) = scor
      rate_eval % unscreened_rates(i_dscor_dt,69) = dscor_dt


      call screen5(pstate, 36, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,71) = scor
      rate_eval % unscreened_rates(i_dscor_dt,71) = dscor_dt
      rate_eval % unscreened_rates(i_scor,72) = scor
      rate_eval % unscreened_rates(i_dscor_dt,72) = dscor_dt


      call screen5(pstate, 37, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,73) = scor
      rate_eval % unscreened_rates(i_dscor_dt,73) = dscor_dt


      call screen5(pstate, 38, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,96) = scor
      rate_eval % unscreened_rates(i_dscor_dt,96) = dscor_dt

    end if


    ! Compute screened rates
    rate_eval % screened_rates = rate_eval % unscreened_rates(i_rate, :) * &
                                 rate_eval % unscreened_rates(i_scor, :)

  end subroutine evaluate_rates


  subroutine actual_rhs(state, ydot)
    
    !$acc routine seq

    use extern_probin_module, only: do_constant_volume_burn
    use burn_type_module, only: net_itemp, net_ienuc, neqs
    use sneut_module, only: sneut5
    use temperature_integration_module, only: temperature_rhs

    implicit none

    type(burn_t), intent(in) :: state
    real(rt), intent(inout) :: ydot(neqs)

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
    ydot(1:nspec) = ydot_nuc

    ! ion binding energy contributions
    call ener_gener_rate(ydot_nuc, enuc)

    ! additional per-reaction energies
    ! including Q-value modification and electron chemical potential

    ! additional energy generation rates
    ! including gamma heating and reaction neutrino losses (non-thermal)


    ! Get the thermal neutrino losses
    call sneut5(state % T, state % rho, state % abar, state % zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)

    ! Append the energy equation (this is erg/g/s)
    ydot(net_ienuc) = enuc - sneut

    ! Append the temperature equation
    call temperature_rhs(state, ydot)

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
      screened_rates(k_ar36__p_cl35)*Y(jar36) + screened_rates(k_c12_ne20__p_p31)*Y(jc12)* &
      Y(jne20)*state % rho + screened_rates(k_c12_o16__p_al27)*Y(jc12)* &
      Y(jo16)*state % rho + screened_rates(k_ca40__p_k39)*Y(jca40) + &
      screened_rates(k_cr48__p_v47)*Y(jcr48) + screened_rates(k_fe52__p_mn51)* &
      Y(jfe52) + screened_rates(k_he4_ar36__p_k39)*Y(jar36)*Y(jhe4)* &
      state % rho + screened_rates(k_he4_ca40__p_sc43)*Y(jca40)*Y(jhe4)* &
      state % rho + screened_rates(k_he4_cr48__p_mn51)*Y(jcr48)*Y(jhe4)* &
      state % rho + screened_rates(k_he4_f18__p_ne21)*Y(jf18)*Y(jhe4)*state % rho &
      + screened_rates(k_he4_fe52__p_co55)*Y(jfe52)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mg24__p_al27)*Y(jhe4)*Y(jmg24)*state % rho + &
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*Y(jn13)*state % rho + &
      screened_rates(k_he4_s32__p_cl35)*Y(jhe4)*Y(js32)*state % rho + &
      screened_rates(k_he4_si28__p_p31)*Y(jhe4)*Y(jsi28)*state % rho + &
      screened_rates(k_he4_ti44__p_v47)*Y(jhe4)*Y(jti44)*state % rho + &
      screened_rates(k_n13__p_c12)*Y(jn13) + screened_rates(k_ni56__p_co55)*Y(jni56) &
      + 0.5d0*screened_rates(k_o16_o16__p_p31)*Y(jo16)**2*state % rho - &
      screened_rates(k_p_al27__c12_o16)*Y(jal27)*Y(jp)*state % rho - &
      screened_rates(k_p_al27__he4_mg24)*Y(jal27)*Y(jp)*state % rho - &
      screened_rates(k_p_al27__si28)*Y(jal27)*Y(jp)*state % rho - &
      screened_rates(k_p_c12__n13)*Y(jc12)*Y(jp)*state % rho - &
      screened_rates(k_p_cl35__ar36)*Y(jcl35)*Y(jp)*state % rho - &
      screened_rates(k_p_cl35__he4_s32)*Y(jcl35)*Y(jp)*state % rho - &
      screened_rates(k_p_co55__he4_fe52)*Y(jco55)*Y(jp)*state % rho - &
      screened_rates(k_p_co55__ni56)*Y(jco55)*Y(jp)*state % rho - &
      screened_rates(k_p_k39__ca40)*Y(jk39)*Y(jp)*state % rho - &
      screened_rates(k_p_k39__he4_ar36)*Y(jk39)*Y(jp)*state % rho - &
      screened_rates(k_p_mn51__fe52)*Y(jmn51)*Y(jp)*state % rho - &
      screened_rates(k_p_mn51__he4_cr48)*Y(jmn51)*Y(jp)*state % rho - &
      screened_rates(k_p_ne21__he4_f18)*Y(jne21)*Y(jp)*state % rho - &
      screened_rates(k_p_o16__he4_n13)*Y(jo16)*Y(jp)*state % rho - &
      screened_rates(k_p_p31__c12_ne20)*Y(jp31)*Y(jp)*state % rho - &
      screened_rates(k_p_p31__he4_si28)*Y(jp31)*Y(jp)*state % rho - &
      screened_rates(k_p_p31__o16_o16)*Y(jp31)*Y(jp)*state % rho - &
      screened_rates(k_p_p31__s32)*Y(jp31)*Y(jp)*state % rho - &
      screened_rates(k_p_sc43__he4_ca40)*Y(jp)*Y(jsc43)*state % rho - &
      screened_rates(k_p_sc43__ti44)*Y(jp)*Y(jsc43)*state % rho - &
      screened_rates(k_p_v47__cr48)*Y(jp)*Y(jv47)*state % rho - &
      screened_rates(k_p_v47__he4_ti44)*Y(jp)*Y(jv47)*state % rho + &
      screened_rates(k_s32__p_p31)*Y(js32) + screened_rates(k_si28__p_al27)*Y(jsi28) &
      + screened_rates(k_ti44__p_sc43)*Y(jti44) &
       )

    ydot_nuc(jhe4) = ( &
      screened_rates(k_ar36__he4_s32)*Y(jar36) + 3.0d0*screened_rates(k_c12__he4_he4_he4)* &
      Y(jc12) + 0.5d0*screened_rates(k_c12_c12__he4_ne20)*Y(jc12)**2* &
      state % rho + screened_rates(k_c12_ne20__he4_si28)*Y(jc12)*Y(jne20)* &
      state % rho + screened_rates(k_c12_o16__he4_mg24)*Y(jc12)*Y(jo16)* &
      state % rho + screened_rates(k_ca40__he4_ar36)*Y(jca40) + &
      screened_rates(k_cl35__he4_p31)*Y(jcl35) + screened_rates(k_co55__he4_mn51)* &
      Y(jco55) + screened_rates(k_cr48__he4_ti44)*Y(jcr48) + &
      screened_rates(k_f18__he4_n14)*Y(jf18) + screened_rates(k_fe52__he4_cr48)* &
      Y(jfe52) - screened_rates(k_he4_al27__p31)*Y(jal27)*Y(jhe4)*state % rho &
      - screened_rates(k_he4_ar36__ca40)*Y(jar36)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar36__p_k39)*Y(jar36)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_c14__o18)*Y(jc14)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca40__p_sc43)*Y(jca40)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca40__ti44)*Y(jca40)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cl35__k39)*Y(jcl35)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr48__fe52)*Y(jcr48)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr48__p_mn51)*Y(jcr48)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_f18__p_ne21)*Y(jf18)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe52__ni56)*Y(jfe52)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe52__p_co55)*Y(jfe52)*Y(jhe4)*state % rho - 0.5d0* &
      screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**3*state % rho**2 - &
      screened_rates(k_he4_k39__sc43)*Y(jhe4)*Y(jk39)*state % rho - &
      screened_rates(k_he4_mg24__c12_o16)*Y(jhe4)*Y(jmg24)*state % rho - &
      screened_rates(k_he4_mg24__p_al27)*Y(jhe4)*Y(jmg24)*state % rho - &
      screened_rates(k_he4_mg24__si28)*Y(jhe4)*Y(jmg24)*state % rho - &
      screened_rates(k_he4_mn51__co55)*Y(jhe4)*Y(jmn51)*state % rho - &
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*Y(jn13)*state % rho - &
      screened_rates(k_he4_n14__f18)*Y(jhe4)*Y(jn14)*state % rho - &
      screened_rates(k_he4_ne20__c12_c12)*Y(jhe4)*Y(jne20)*state % rho - &
      screened_rates(k_he4_ne20__mg24)*Y(jhe4)*Y(jne20)*state % rho - &
      screened_rates(k_he4_o16__ne20)*Y(jhe4)*Y(jo16)*state % rho - &
      screened_rates(k_he4_p31__cl35)*Y(jhe4)*Y(jp31)*state % rho - &
      screened_rates(k_he4_s32__ar36)*Y(jhe4)*Y(js32)*state % rho - &
      screened_rates(k_he4_s32__p_cl35)*Y(jhe4)*Y(js32)*state % rho - &
      screened_rates(k_he4_sc43__v47)*Y(jhe4)*Y(jsc43)*state % rho - &
      screened_rates(k_he4_si28__c12_ne20)*Y(jhe4)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_si28__o16_o16)*Y(jhe4)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_si28__p_p31)*Y(jhe4)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_si28__s32)*Y(jhe4)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_ti44__cr48)*Y(jhe4)*Y(jti44)*state % rho - &
      screened_rates(k_he4_ti44__p_v47)*Y(jhe4)*Y(jti44)*state % rho - &
      screened_rates(k_he4_v47__mn51)*Y(jhe4)*Y(jv47)*state % rho + &
      screened_rates(k_k39__he4_cl35)*Y(jk39) + screened_rates(k_mg24__he4_ne20)* &
      Y(jmg24) + screened_rates(k_mn51__he4_v47)*Y(jmn51) + &
      screened_rates(k_ne20__he4_o16)*Y(jne20) + screened_rates(k_ni56__he4_fe52)* &
      Y(jni56) + screened_rates(k_o16__he4_c12)*Y(jo16) + 0.5d0* &
      screened_rates(k_o16_o16__he4_si28)*Y(jo16)**2*state % rho + &
      screened_rates(k_o18__he4_c14)*Y(jo18) + screened_rates(k_p31__he4_al27)* &
      Y(jp31) + screened_rates(k_p_al27__he4_mg24)*Y(jal27)*Y(jp)*state % rho &
      + screened_rates(k_p_cl35__he4_s32)*Y(jcl35)*Y(jp)*state % rho + &
      screened_rates(k_p_co55__he4_fe52)*Y(jco55)*Y(jp)*state % rho + &
      screened_rates(k_p_k39__he4_ar36)*Y(jk39)*Y(jp)*state % rho + &
      screened_rates(k_p_mn51__he4_cr48)*Y(jmn51)*Y(jp)*state % rho + &
      screened_rates(k_p_ne21__he4_f18)*Y(jne21)*Y(jp)*state % rho + &
      screened_rates(k_p_o16__he4_n13)*Y(jo16)*Y(jp)*state % rho + &
      screened_rates(k_p_p31__he4_si28)*Y(jp31)*Y(jp)*state % rho + &
      screened_rates(k_p_sc43__he4_ca40)*Y(jp)*Y(jsc43)*state % rho + &
      screened_rates(k_p_v47__he4_ti44)*Y(jp)*Y(jv47)*state % rho + &
      screened_rates(k_s32__he4_si28)*Y(js32) + screened_rates(k_sc43__he4_k39)* &
      Y(jsc43) + screened_rates(k_si28__he4_mg24)*Y(jsi28) + &
      screened_rates(k_ti44__he4_ca40)*Y(jti44) + screened_rates(k_v47__he4_sc43)* &
      Y(jv47) &
       )

    ydot_nuc(jc12) = ( &
      -screened_rates(k_c12__he4_he4_he4)*Y(jc12) - screened_rates(k_c12_c12__he4_ne20)* &
      Y(jc12)**2*state % rho - screened_rates(k_c12_ne20__he4_si28)*Y(jc12)* &
      Y(jne20)*state % rho - screened_rates(k_c12_ne20__p_p31)*Y(jc12)* &
      Y(jne20)*state % rho - screened_rates(k_c12_o16__he4_mg24)*Y(jc12)* &
      Y(jo16)*state % rho - screened_rates(k_c12_o16__p_al27)*Y(jc12)*Y(jo16) &
      *state % rho - screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*state % rho + &
      0.16666666666666667d0*screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**3* &
      state % rho**2 + screened_rates(k_he4_mg24__c12_o16)*Y(jhe4)*Y(jmg24)* &
      state % rho + 2.0d0*screened_rates(k_he4_ne20__c12_c12)*Y(jhe4)*Y(jne20)* &
      state % rho + screened_rates(k_he4_si28__c12_ne20)*Y(jhe4)*Y(jsi28)* &
      state % rho + screened_rates(k_n13__p_c12)*Y(jn13) + screened_rates(k_o16__he4_c12)* &
      Y(jo16) + screened_rates(k_p_al27__c12_o16)*Y(jal27)*Y(jp)*state % rho &
      - screened_rates(k_p_c12__n13)*Y(jc12)*Y(jp)*state % rho + &
      screened_rates(k_p_p31__c12_ne20)*Y(jp31)*Y(jp)*state % rho &
       )

    ydot_nuc(jc14) = ( &
      -screened_rates(k_c14__n14__weak__wc12)*Y(jc14) - screened_rates(k_he4_c14__o18)* &
      Y(jc14)*Y(jhe4)*state % rho + screened_rates(k_o18__he4_c14)*Y(jo18) &
       )

    ydot_nuc(jn13) = ( &
      -screened_rates(k_he4_n13__p_o16)*Y(jhe4)*Y(jn13)*state % rho - &
      screened_rates(k_n13__p_c12)*Y(jn13) + screened_rates(k_p_c12__n13)*Y(jc12)* &
      Y(jp)*state % rho + screened_rates(k_p_o16__he4_n13)*Y(jo16)*Y(jp)* &
      state % rho &
       )

    ydot_nuc(jn14) = ( &
      screened_rates(k_c14__n14__weak__wc12)*Y(jc14) + screened_rates(k_f18__he4_n14)* &
      Y(jf18) - screened_rates(k_he4_n14__f18)*Y(jhe4)*Y(jn14)*state % rho &
       )

    ydot_nuc(jo16) = ( &
      -screened_rates(k_c12_o16__he4_mg24)*Y(jc12)*Y(jo16)*state % rho - &
      screened_rates(k_c12_o16__p_al27)*Y(jc12)*Y(jo16)*state % rho + &
      screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mg24__c12_o16)*Y(jhe4)*Y(jmg24)*state % rho + &
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*Y(jn13)*state % rho - &
      screened_rates(k_he4_o16__ne20)*Y(jhe4)*Y(jo16)*state % rho + 2.0d0* &
      screened_rates(k_he4_si28__o16_o16)*Y(jhe4)*Y(jsi28)*state % rho + &
      screened_rates(k_ne20__he4_o16)*Y(jne20) - screened_rates(k_o16__he4_c12)* &
      Y(jo16) - screened_rates(k_o16_o16__he4_si28)*Y(jo16)**2*state % rho - &
      screened_rates(k_o16_o16__p_p31)*Y(jo16)**2*state % rho + &
      screened_rates(k_p_al27__c12_o16)*Y(jal27)*Y(jp)*state % rho - &
      screened_rates(k_p_o16__he4_n13)*Y(jo16)*Y(jp)*state % rho + 2.0d0* &
      screened_rates(k_p_p31__o16_o16)*Y(jp31)*Y(jp)*state % rho &
       )

    ydot_nuc(jo18) = ( &
      screened_rates(k_f18__o18__weak__wc12)*Y(jf18) + screened_rates(k_he4_c14__o18)* &
      Y(jc14)*Y(jhe4)*state % rho - screened_rates(k_o18__he4_c14)*Y(jo18) &
       )

    ydot_nuc(jf18) = ( &
      -screened_rates(k_f18__he4_n14)*Y(jf18) - screened_rates(k_f18__o18__weak__wc12)* &
      Y(jf18) - screened_rates(k_he4_f18__p_ne21)*Y(jf18)*Y(jhe4)*state % rho &
      + screened_rates(k_he4_n14__f18)*Y(jhe4)*Y(jn14)*state % rho + &
      screened_rates(k_p_ne21__he4_f18)*Y(jne21)*Y(jp)*state % rho &
       )

    ydot_nuc(jne20) = ( &
      0.5d0*screened_rates(k_c12_c12__he4_ne20)*Y(jc12)**2*state % rho - &
      screened_rates(k_c12_ne20__he4_si28)*Y(jc12)*Y(jne20)*state % rho - &
      screened_rates(k_c12_ne20__p_p31)*Y(jc12)*Y(jne20)*state % rho - &
      screened_rates(k_he4_ne20__c12_c12)*Y(jhe4)*Y(jne20)*state % rho - &
      screened_rates(k_he4_ne20__mg24)*Y(jhe4)*Y(jne20)*state % rho + &
      screened_rates(k_he4_o16__ne20)*Y(jhe4)*Y(jo16)*state % rho + &
      screened_rates(k_he4_si28__c12_ne20)*Y(jhe4)*Y(jsi28)*state % rho + &
      screened_rates(k_mg24__he4_ne20)*Y(jmg24) - screened_rates(k_ne20__he4_o16)* &
      Y(jne20) + screened_rates(k_p_p31__c12_ne20)*Y(jp31)*Y(jp)*state % rho &
       )

    ydot_nuc(jne21) = ( &
      screened_rates(k_he4_f18__p_ne21)*Y(jf18)*Y(jhe4)*state % rho - &
      screened_rates(k_p_ne21__he4_f18)*Y(jne21)*Y(jp)*state % rho &
       )

    ydot_nuc(jmg24) = ( &
      screened_rates(k_c12_o16__he4_mg24)*Y(jc12)*Y(jo16)*state % rho - &
      screened_rates(k_he4_mg24__c12_o16)*Y(jhe4)*Y(jmg24)*state % rho - &
      screened_rates(k_he4_mg24__p_al27)*Y(jhe4)*Y(jmg24)*state % rho - &
      screened_rates(k_he4_mg24__si28)*Y(jhe4)*Y(jmg24)*state % rho + &
      screened_rates(k_he4_ne20__mg24)*Y(jhe4)*Y(jne20)*state % rho - &
      screened_rates(k_mg24__he4_ne20)*Y(jmg24) + screened_rates(k_p_al27__he4_mg24)* &
      Y(jal27)*Y(jp)*state % rho + screened_rates(k_si28__he4_mg24)*Y(jsi28) &
       )

    ydot_nuc(jal27) = ( &
      screened_rates(k_c12_o16__p_al27)*Y(jc12)*Y(jo16)*state % rho - &
      screened_rates(k_he4_al27__p31)*Y(jal27)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_mg24__p_al27)*Y(jhe4)*Y(jmg24)*state % rho + &
      screened_rates(k_p31__he4_al27)*Y(jp31) - screened_rates(k_p_al27__c12_o16)* &
      Y(jal27)*Y(jp)*state % rho - screened_rates(k_p_al27__he4_mg24)* &
      Y(jal27)*Y(jp)*state % rho - screened_rates(k_p_al27__si28)*Y(jal27)* &
      Y(jp)*state % rho + screened_rates(k_si28__p_al27)*Y(jsi28) &
       )

    ydot_nuc(jsi28) = ( &
      screened_rates(k_c12_ne20__he4_si28)*Y(jc12)*Y(jne20)*state % rho + &
      screened_rates(k_he4_mg24__si28)*Y(jhe4)*Y(jmg24)*state % rho - &
      screened_rates(k_he4_si28__c12_ne20)*Y(jhe4)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_si28__o16_o16)*Y(jhe4)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_si28__p_p31)*Y(jhe4)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_si28__s32)*Y(jhe4)*Y(jsi28)*state % rho + 0.5d0* &
      screened_rates(k_o16_o16__he4_si28)*Y(jo16)**2*state % rho + &
      screened_rates(k_p_al27__si28)*Y(jal27)*Y(jp)*state % rho + &
      screened_rates(k_p_p31__he4_si28)*Y(jp31)*Y(jp)*state % rho + &
      screened_rates(k_s32__he4_si28)*Y(js32) - screened_rates(k_si28__he4_mg24)* &
      Y(jsi28) - screened_rates(k_si28__p_al27)*Y(jsi28) &
       )

    ydot_nuc(jp31) = ( &
      screened_rates(k_c12_ne20__p_p31)*Y(jc12)*Y(jne20)*state % rho + &
      screened_rates(k_cl35__he4_p31)*Y(jcl35) + screened_rates(k_he4_al27__p31)* &
      Y(jal27)*Y(jhe4)*state % rho - screened_rates(k_he4_p31__cl35)*Y(jhe4)* &
      Y(jp31)*state % rho + screened_rates(k_he4_si28__p_p31)*Y(jhe4)* &
      Y(jsi28)*state % rho + 0.5d0*screened_rates(k_o16_o16__p_p31)*Y(jo16)**2* &
      state % rho - screened_rates(k_p31__he4_al27)*Y(jp31) - &
      screened_rates(k_p_p31__c12_ne20)*Y(jp31)*Y(jp)*state % rho - &
      screened_rates(k_p_p31__he4_si28)*Y(jp31)*Y(jp)*state % rho - &
      screened_rates(k_p_p31__o16_o16)*Y(jp31)*Y(jp)*state % rho - &
      screened_rates(k_p_p31__s32)*Y(jp31)*Y(jp)*state % rho + &
      screened_rates(k_s32__p_p31)*Y(js32) &
       )

    ydot_nuc(js32) = ( &
      screened_rates(k_ar36__he4_s32)*Y(jar36) - screened_rates(k_he4_s32__ar36)*Y(jhe4)* &
      Y(js32)*state % rho - screened_rates(k_he4_s32__p_cl35)*Y(jhe4)*Y(js32) &
      *state % rho + screened_rates(k_he4_si28__s32)*Y(jhe4)*Y(jsi28)*state % rho &
      + screened_rates(k_p_cl35__he4_s32)*Y(jcl35)*Y(jp)*state % rho + &
      screened_rates(k_p_p31__s32)*Y(jp31)*Y(jp)*state % rho - &
      screened_rates(k_s32__he4_si28)*Y(js32) - screened_rates(k_s32__p_p31)*Y(js32) &
       )

    ydot_nuc(jcl35) = ( &
      screened_rates(k_ar36__p_cl35)*Y(jar36) - screened_rates(k_cl35__he4_p31)*Y(jcl35) - &
      screened_rates(k_he4_cl35__k39)*Y(jcl35)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_p31__cl35)*Y(jhe4)*Y(jp31)*state % rho + &
      screened_rates(k_he4_s32__p_cl35)*Y(jhe4)*Y(js32)*state % rho + &
      screened_rates(k_k39__he4_cl35)*Y(jk39) - screened_rates(k_p_cl35__ar36)* &
      Y(jcl35)*Y(jp)*state % rho - screened_rates(k_p_cl35__he4_s32)*Y(jcl35) &
      *Y(jp)*state % rho &
       )

    ydot_nuc(jar36) = ( &
      -screened_rates(k_ar36__he4_s32)*Y(jar36) - screened_rates(k_ar36__p_cl35)*Y(jar36) + &
      screened_rates(k_ca40__he4_ar36)*Y(jca40) - screened_rates(k_he4_ar36__ca40)* &
      Y(jar36)*Y(jhe4)*state % rho - screened_rates(k_he4_ar36__p_k39)* &
      Y(jar36)*Y(jhe4)*state % rho + screened_rates(k_he4_s32__ar36)*Y(jhe4)* &
      Y(js32)*state % rho + screened_rates(k_p_cl35__ar36)*Y(jcl35)*Y(jp)* &
      state % rho + screened_rates(k_p_k39__he4_ar36)*Y(jk39)*Y(jp)*state % rho &
       )

    ydot_nuc(jk39) = ( &
      screened_rates(k_ca40__p_k39)*Y(jca40) + screened_rates(k_he4_ar36__p_k39)*Y(jar36)* &
      Y(jhe4)*state % rho + screened_rates(k_he4_cl35__k39)*Y(jcl35)*Y(jhe4)* &
      state % rho - screened_rates(k_he4_k39__sc43)*Y(jhe4)*Y(jk39)*state % rho - &
      screened_rates(k_k39__he4_cl35)*Y(jk39) - screened_rates(k_p_k39__ca40)*Y(jk39) &
      *Y(jp)*state % rho - screened_rates(k_p_k39__he4_ar36)*Y(jk39)*Y(jp)* &
      state % rho + screened_rates(k_sc43__he4_k39)*Y(jsc43) &
       )

    ydot_nuc(jca40) = ( &
      -screened_rates(k_ca40__he4_ar36)*Y(jca40) - screened_rates(k_ca40__p_k39)*Y(jca40) + &
      screened_rates(k_he4_ar36__ca40)*Y(jar36)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca40__p_sc43)*Y(jca40)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca40__ti44)*Y(jca40)*Y(jhe4)*state % rho + &
      screened_rates(k_p_k39__ca40)*Y(jk39)*Y(jp)*state % rho + &
      screened_rates(k_p_sc43__he4_ca40)*Y(jp)*Y(jsc43)*state % rho + &
      screened_rates(k_ti44__he4_ca40)*Y(jti44) &
       )

    ydot_nuc(jsc43) = ( &
      screened_rates(k_he4_ca40__p_sc43)*Y(jca40)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_k39__sc43)*Y(jhe4)*Y(jk39)*state % rho - &
      screened_rates(k_he4_sc43__v47)*Y(jhe4)*Y(jsc43)*state % rho - &
      screened_rates(k_p_sc43__he4_ca40)*Y(jp)*Y(jsc43)*state % rho - &
      screened_rates(k_p_sc43__ti44)*Y(jp)*Y(jsc43)*state % rho - &
      screened_rates(k_sc43__he4_k39)*Y(jsc43) + screened_rates(k_ti44__p_sc43)* &
      Y(jti44) + screened_rates(k_v47__he4_sc43)*Y(jv47) &
       )

    ydot_nuc(jti44) = ( &
      screened_rates(k_cr48__he4_ti44)*Y(jcr48) + screened_rates(k_he4_ca40__ti44)*Y(jca40)* &
      Y(jhe4)*state % rho - screened_rates(k_he4_ti44__cr48)*Y(jhe4)*Y(jti44) &
      *state % rho - screened_rates(k_he4_ti44__p_v47)*Y(jhe4)*Y(jti44)* &
      state % rho + screened_rates(k_p_sc43__ti44)*Y(jp)*Y(jsc43)*state % rho + &
      screened_rates(k_p_v47__he4_ti44)*Y(jp)*Y(jv47)*state % rho - &
      screened_rates(k_ti44__he4_ca40)*Y(jti44) - screened_rates(k_ti44__p_sc43)* &
      Y(jti44) &
       )

    ydot_nuc(jv47) = ( &
      screened_rates(k_cr48__p_v47)*Y(jcr48) + screened_rates(k_he4_sc43__v47)*Y(jhe4)* &
      Y(jsc43)*state % rho + screened_rates(k_he4_ti44__p_v47)*Y(jhe4)* &
      Y(jti44)*state % rho - screened_rates(k_he4_v47__mn51)*Y(jhe4)*Y(jv47)* &
      state % rho + screened_rates(k_mn51__he4_v47)*Y(jmn51) - &
      screened_rates(k_p_v47__cr48)*Y(jp)*Y(jv47)*state % rho - &
      screened_rates(k_p_v47__he4_ti44)*Y(jp)*Y(jv47)*state % rho - &
      screened_rates(k_v47__he4_sc43)*Y(jv47) &
       )

    ydot_nuc(jcr48) = ( &
      -screened_rates(k_cr48__he4_ti44)*Y(jcr48) - screened_rates(k_cr48__p_v47)*Y(jcr48) + &
      screened_rates(k_fe52__he4_cr48)*Y(jfe52) - screened_rates(k_he4_cr48__fe52)* &
      Y(jcr48)*Y(jhe4)*state % rho - screened_rates(k_he4_cr48__p_mn51)* &
      Y(jcr48)*Y(jhe4)*state % rho + screened_rates(k_he4_ti44__cr48)*Y(jhe4) &
      *Y(jti44)*state % rho + screened_rates(k_p_mn51__he4_cr48)*Y(jmn51)* &
      Y(jp)*state % rho + screened_rates(k_p_v47__cr48)*Y(jp)*Y(jv47)* &
      state % rho &
       )

    ydot_nuc(jmn51) = ( &
      screened_rates(k_co55__he4_mn51)*Y(jco55) + screened_rates(k_fe52__p_mn51)*Y(jfe52) + &
      screened_rates(k_he4_cr48__p_mn51)*Y(jcr48)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_mn51__co55)*Y(jhe4)*Y(jmn51)*state % rho + &
      screened_rates(k_he4_v47__mn51)*Y(jhe4)*Y(jv47)*state % rho - &
      screened_rates(k_mn51__he4_v47)*Y(jmn51) - screened_rates(k_p_mn51__fe52)* &
      Y(jmn51)*Y(jp)*state % rho - screened_rates(k_p_mn51__he4_cr48)* &
      Y(jmn51)*Y(jp)*state % rho &
       )

    ydot_nuc(jfe52) = ( &
      -screened_rates(k_fe52__he4_cr48)*Y(jfe52) - screened_rates(k_fe52__p_mn51)*Y(jfe52) + &
      screened_rates(k_he4_cr48__fe52)*Y(jcr48)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe52__ni56)*Y(jfe52)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe52__p_co55)*Y(jfe52)*Y(jhe4)*state % rho + &
      screened_rates(k_ni56__he4_fe52)*Y(jni56) + screened_rates(k_p_co55__he4_fe52)* &
      Y(jco55)*Y(jp)*state % rho + screened_rates(k_p_mn51__fe52)*Y(jmn51)* &
      Y(jp)*state % rho &
       )

    ydot_nuc(jco55) = ( &
      -screened_rates(k_co55__he4_mn51)*Y(jco55) + screened_rates(k_he4_fe52__p_co55)* &
      Y(jfe52)*Y(jhe4)*state % rho + screened_rates(k_he4_mn51__co55)*Y(jhe4) &
      *Y(jmn51)*state % rho + screened_rates(k_ni56__p_co55)*Y(jni56) - &
      screened_rates(k_p_co55__he4_fe52)*Y(jco55)*Y(jp)*state % rho - &
      screened_rates(k_p_co55__ni56)*Y(jco55)*Y(jp)*state % rho &
       )

    ydot_nuc(jni56) = ( &
      screened_rates(k_he4_fe52__ni56)*Y(jfe52)*Y(jhe4)*state % rho - &
      screened_rates(k_ni56__he4_fe52)*Y(jni56) - screened_rates(k_ni56__p_co55)* &
      Y(jni56) + screened_rates(k_p_co55__ni56)*Y(jco55)*Y(jp)*state % rho &
       )


  end subroutine rhs_nuc


  subroutine actual_jac(state, jac)

    !$acc routine seq

    use burn_type_module, only: net_itemp, net_ienuc, neqs
    use sneut_module, only: sneut5
    use temperature_integration_module, only: temperature_jac
    use jacobian_sparsity_module, only: get_jac_entry, set_jac_entry, set_jac_zero

    implicit none
    
    type(burn_t), intent(in) :: state
    real(rt) :: jac(neqs, neqs)

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
    call set_jac_zero(jac)

    ! Species Jacobian elements with respect to other species
    call jac_nuc(state, jac, Y, rate_eval % screened_rates)

    ! Evaluate the species Jacobian elements with respect to temperature by
    ! calling the RHS using the temperature derivative of the screened rate
    screened_rates_dt = rate_eval % unscreened_rates(i_rate, :) * &
                        rate_eval % unscreened_rates(i_dscor_dt, :) + &
                        rate_eval % unscreened_rates(i_drate_dt, :) * &
                        rate_eval % unscreened_rates(i_scor, :)

    call rhs_nuc(state, yderivs, Y, screened_rates_dt)

    do k = 1, nspec
       call set_jac_entry(jac, k, net_itemp, yderivs(k))
    enddo

    ! Energy generation rate Jacobian elements with respect to species
    do j = 1, nspec
       do k = 1, nspec
          call get_jac_entry(jac, k, j, yderivs(k))
       enddo
       call ener_gener_rate(yderivs, scratch)
       call set_jac_entry(jac, net_ienuc, j, scratch)
    enddo

    ! Account for the thermal neutrino losses
    call sneut5(state % T, state % rho, state % abar, state % zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)

    do j = 1, nspec
       b1 = ((aion(j) - state % abar) * state % abar * snuda + (zion(j) - state % zbar) * state % abar * snudz)
       call get_jac_entry(jac, net_ienuc, j, scratch)
       scratch = scratch - b1
       call set_jac_entry(jac, net_ienuc, j, scratch)
    enddo

    ! Energy generation rate Jacobian element with respect to temperature
    do k = 1, nspec
       call get_jac_entry(jac, k, net_itemp, yderivs(k))
    enddo
    call ener_gener_rate(yderivs, scratch)
    scratch = scratch - dsneutdt    
    call set_jac_entry(jac, net_ienuc, net_itemp, scratch)

    ! Temperature Jacobian elements
    call temperature_jac(state, jac)

  end subroutine actual_jac


  subroutine jac_nuc(state, jac, Y, screened_rates)

    !$acc routine seq

    use jacobian_sparsity_module, only: set_jac_entry

    implicit none

    type(burn_t), intent(in) :: state
    real(rt), intent(inout) :: jac(neqs, neqs)

    real(rt), intent(in)  :: Y(nspec)
    real(rt), intent(in)  :: screened_rates(nrates)
    real(rt) :: scratch


    !$gpu


    scratch = (&
      -screened_rates(k_p_al27__c12_o16)*Y(jal27)*state % rho - &
      screened_rates(k_p_al27__he4_mg24)*Y(jal27)*state % rho - &
      screened_rates(k_p_al27__si28)*Y(jal27)*state % rho - screened_rates(k_p_c12__n13)* &
      Y(jc12)*state % rho - screened_rates(k_p_cl35__ar36)*Y(jcl35)*state % rho - &
      screened_rates(k_p_cl35__he4_s32)*Y(jcl35)*state % rho - &
      screened_rates(k_p_co55__he4_fe52)*Y(jco55)*state % rho - &
      screened_rates(k_p_co55__ni56)*Y(jco55)*state % rho - screened_rates(k_p_k39__ca40)* &
      Y(jk39)*state % rho - screened_rates(k_p_k39__he4_ar36)*Y(jk39)*state % rho &
      - screened_rates(k_p_mn51__fe52)*Y(jmn51)*state % rho - &
      screened_rates(k_p_mn51__he4_cr48)*Y(jmn51)*state % rho - &
      screened_rates(k_p_ne21__he4_f18)*Y(jne21)*state % rho - &
      screened_rates(k_p_o16__he4_n13)*Y(jo16)*state % rho - &
      screened_rates(k_p_p31__c12_ne20)*Y(jp31)*state % rho - &
      screened_rates(k_p_p31__he4_si28)*Y(jp31)*state % rho - &
      screened_rates(k_p_p31__o16_o16)*Y(jp31)*state % rho - screened_rates(k_p_p31__s32)* &
      Y(jp31)*state % rho - screened_rates(k_p_sc43__he4_ca40)*Y(jsc43)* &
      state % rho - screened_rates(k_p_sc43__ti44)*Y(jsc43)*state % rho - &
      screened_rates(k_p_v47__cr48)*Y(jv47)*state % rho - &
      screened_rates(k_p_v47__he4_ti44)*Y(jv47)*state % rho &
       )
    call set_jac_entry(jac, jp, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar36__p_k39)*Y(jar36)*state % rho + screened_rates(k_he4_ca40__p_sc43) &
      *Y(jca40)*state % rho + screened_rates(k_he4_cr48__p_mn51)*Y(jcr48)* &
      state % rho + screened_rates(k_he4_f18__p_ne21)*Y(jf18)*state % rho + &
      screened_rates(k_he4_fe52__p_co55)*Y(jfe52)*state % rho + &
      screened_rates(k_he4_mg24__p_al27)*Y(jmg24)*state % rho + &
      screened_rates(k_he4_n13__p_o16)*Y(jn13)*state % rho + &
      screened_rates(k_he4_s32__p_cl35)*Y(js32)*state % rho + &
      screened_rates(k_he4_si28__p_p31)*Y(jsi28)*state % rho + &
      screened_rates(k_he4_ti44__p_v47)*Y(jti44)*state % rho &
       )
    call set_jac_entry(jac, jp, jhe4, scratch)

    scratch = (&
      screened_rates(k_c12_ne20__p_p31)*Y(jne20)*state % rho + screened_rates(k_c12_o16__p_al27)* &
      Y(jo16)*state % rho - screened_rates(k_p_c12__n13)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jp, jc12, scratch)

    scratch = (&
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*state % rho + screened_rates(k_n13__p_c12) &
       )
    call set_jac_entry(jac, jp, jn13, scratch)

    scratch = (&
      screened_rates(k_c12_o16__p_al27)*Y(jc12)*state % rho + 1.0d0* &
      screened_rates(k_o16_o16__p_p31)*Y(jo16)*state % rho - &
      screened_rates(k_p_o16__he4_n13)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jp, jo16, scratch)

    scratch = (&
      screened_rates(k_he4_f18__p_ne21)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jp, jf18, scratch)

    scratch = (&
      screened_rates(k_c12_ne20__p_p31)*Y(jc12)*state % rho &
       )
    call set_jac_entry(jac, jp, jne20, scratch)

    scratch = (&
      -screened_rates(k_p_ne21__he4_f18)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jp, jne21, scratch)

    scratch = (&
      screened_rates(k_he4_mg24__p_al27)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jp, jmg24, scratch)

    scratch = (&
      -screened_rates(k_p_al27__c12_o16)*Y(jp)*state % rho - screened_rates(k_p_al27__he4_mg24)* &
      Y(jp)*state % rho - screened_rates(k_p_al27__si28)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jp, jal27, scratch)

    scratch = (&
      screened_rates(k_he4_si28__p_p31)*Y(jhe4)*state % rho + screened_rates(k_si28__p_al27) &
       )
    call set_jac_entry(jac, jp, jsi28, scratch)

    scratch = (&
      -screened_rates(k_p_p31__c12_ne20)*Y(jp)*state % rho - screened_rates(k_p_p31__he4_si28)* &
      Y(jp)*state % rho - screened_rates(k_p_p31__o16_o16)*Y(jp)*state % rho - &
      screened_rates(k_p_p31__s32)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jp, jp31, scratch)

    scratch = (&
      screened_rates(k_he4_s32__p_cl35)*Y(jhe4)*state % rho + screened_rates(k_s32__p_p31) &
       )
    call set_jac_entry(jac, jp, js32, scratch)

    scratch = (&
      -screened_rates(k_p_cl35__ar36)*Y(jp)*state % rho - screened_rates(k_p_cl35__he4_s32)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jp, jcl35, scratch)

    scratch = (&
      screened_rates(k_ar36__p_cl35) + screened_rates(k_he4_ar36__p_k39)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jp, jar36, scratch)

    scratch = (&
      -screened_rates(k_p_k39__ca40)*Y(jp)*state % rho - screened_rates(k_p_k39__he4_ar36)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jp, jk39, scratch)

    scratch = (&
      screened_rates(k_ca40__p_k39) + screened_rates(k_he4_ca40__p_sc43)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jp, jca40, scratch)

    scratch = (&
      -screened_rates(k_p_sc43__he4_ca40)*Y(jp)*state % rho - screened_rates(k_p_sc43__ti44)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jp, jsc43, scratch)

    scratch = (&
      screened_rates(k_he4_ti44__p_v47)*Y(jhe4)*state % rho + screened_rates(k_ti44__p_sc43) &
       )
    call set_jac_entry(jac, jp, jti44, scratch)

    scratch = (&
      -screened_rates(k_p_v47__cr48)*Y(jp)*state % rho - screened_rates(k_p_v47__he4_ti44)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jp, jv47, scratch)

    scratch = (&
      screened_rates(k_cr48__p_v47) + screened_rates(k_he4_cr48__p_mn51)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jp, jcr48, scratch)

    scratch = (&
      -screened_rates(k_p_mn51__fe52)*Y(jp)*state % rho - screened_rates(k_p_mn51__he4_cr48)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jp, jmn51, scratch)

    scratch = (&
      screened_rates(k_fe52__p_mn51) + screened_rates(k_he4_fe52__p_co55)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jp, jfe52, scratch)

    scratch = (&
      -screened_rates(k_p_co55__he4_fe52)*Y(jp)*state % rho - screened_rates(k_p_co55__ni56)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jp, jco55, scratch)

    scratch = (&
      screened_rates(k_ni56__p_co55) &
       )
    call set_jac_entry(jac, jp, jni56, scratch)

    scratch = (&
      screened_rates(k_p_al27__he4_mg24)*Y(jal27)*state % rho + screened_rates(k_p_cl35__he4_s32) &
      *Y(jcl35)*state % rho + screened_rates(k_p_co55__he4_fe52)*Y(jco55)* &
      state % rho + screened_rates(k_p_k39__he4_ar36)*Y(jk39)*state % rho + &
      screened_rates(k_p_mn51__he4_cr48)*Y(jmn51)*state % rho + &
      screened_rates(k_p_ne21__he4_f18)*Y(jne21)*state % rho + &
      screened_rates(k_p_o16__he4_n13)*Y(jo16)*state % rho + &
      screened_rates(k_p_p31__he4_si28)*Y(jp31)*state % rho + &
      screened_rates(k_p_sc43__he4_ca40)*Y(jsc43)*state % rho + &
      screened_rates(k_p_v47__he4_ti44)*Y(jv47)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_al27__p31)*Y(jal27)*state % rho - screened_rates(k_he4_ar36__ca40)* &
      Y(jar36)*state % rho - screened_rates(k_he4_ar36__p_k39)*Y(jar36)* &
      state % rho - screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho - &
      screened_rates(k_he4_c14__o18)*Y(jc14)*state % rho - &
      screened_rates(k_he4_ca40__p_sc43)*Y(jca40)*state % rho - &
      screened_rates(k_he4_ca40__ti44)*Y(jca40)*state % rho - &
      screened_rates(k_he4_cl35__k39)*Y(jcl35)*state % rho - &
      screened_rates(k_he4_cr48__fe52)*Y(jcr48)*state % rho - &
      screened_rates(k_he4_cr48__p_mn51)*Y(jcr48)*state % rho - &
      screened_rates(k_he4_f18__p_ne21)*Y(jf18)*state % rho - &
      screened_rates(k_he4_fe52__ni56)*Y(jfe52)*state % rho - &
      screened_rates(k_he4_fe52__p_co55)*Y(jfe52)*state % rho - 1.5d0* &
      screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**2*state % rho**2 - &
      screened_rates(k_he4_k39__sc43)*Y(jk39)*state % rho - &
      screened_rates(k_he4_mg24__c12_o16)*Y(jmg24)*state % rho - &
      screened_rates(k_he4_mg24__p_al27)*Y(jmg24)*state % rho - &
      screened_rates(k_he4_mg24__si28)*Y(jmg24)*state % rho - &
      screened_rates(k_he4_mn51__co55)*Y(jmn51)*state % rho - &
      screened_rates(k_he4_n13__p_o16)*Y(jn13)*state % rho - &
      screened_rates(k_he4_n14__f18)*Y(jn14)*state % rho - &
      screened_rates(k_he4_ne20__c12_c12)*Y(jne20)*state % rho - &
      screened_rates(k_he4_ne20__mg24)*Y(jne20)*state % rho - &
      screened_rates(k_he4_o16__ne20)*Y(jo16)*state % rho - &
      screened_rates(k_he4_p31__cl35)*Y(jp31)*state % rho - &
      screened_rates(k_he4_s32__ar36)*Y(js32)*state % rho - &
      screened_rates(k_he4_s32__p_cl35)*Y(js32)*state % rho - &
      screened_rates(k_he4_sc43__v47)*Y(jsc43)*state % rho - &
      screened_rates(k_he4_si28__c12_ne20)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_si28__o16_o16)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_si28__p_p31)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_si28__s32)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_ti44__cr48)*Y(jti44)*state % rho - &
      screened_rates(k_he4_ti44__p_v47)*Y(jti44)*state % rho - &
      screened_rates(k_he4_v47__mn51)*Y(jv47)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jhe4, scratch)

    scratch = (&
      3.0d0*screened_rates(k_c12__he4_he4_he4) + 1.0d0*screened_rates(k_c12_c12__he4_ne20)* &
      Y(jc12)*state % rho + screened_rates(k_c12_ne20__he4_si28)*Y(jne20)* &
      state % rho + screened_rates(k_c12_o16__he4_mg24)*Y(jo16)*state % rho - &
      screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jc12, scratch)

    scratch = (&
      -screened_rates(k_he4_c14__o18)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jc14, scratch)

    scratch = (&
      -screened_rates(k_he4_n13__p_o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jn13, scratch)

    scratch = (&
      -screened_rates(k_he4_n14__f18)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jn14, scratch)

    scratch = (&
      screened_rates(k_c12_o16__he4_mg24)*Y(jc12)*state % rho - screened_rates(k_he4_o16__ne20)* &
      Y(jhe4)*state % rho + screened_rates(k_o16__he4_c12) + 1.0d0* &
      screened_rates(k_o16_o16__he4_si28)*Y(jo16)*state % rho + &
      screened_rates(k_p_o16__he4_n13)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jo16, scratch)

    scratch = (&
      screened_rates(k_o18__he4_c14) &
       )
    call set_jac_entry(jac, jhe4, jo18, scratch)

    scratch = (&
      screened_rates(k_f18__he4_n14) - screened_rates(k_he4_f18__p_ne21)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jf18, scratch)

    scratch = (&
      screened_rates(k_c12_ne20__he4_si28)*Y(jc12)*state % rho - &
      screened_rates(k_he4_ne20__c12_c12)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ne20__mg24)*Y(jhe4)*state % rho + &
      screened_rates(k_ne20__he4_o16) &
       )
    call set_jac_entry(jac, jhe4, jne20, scratch)

    scratch = (&
      screened_rates(k_p_ne21__he4_f18)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jne21, scratch)

    scratch = (&
      -screened_rates(k_he4_mg24__c12_o16)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_mg24__p_al27)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_mg24__si28)*Y(jhe4)*state % rho + &
      screened_rates(k_mg24__he4_ne20) &
       )
    call set_jac_entry(jac, jhe4, jmg24, scratch)

    scratch = (&
      -screened_rates(k_he4_al27__p31)*Y(jhe4)*state % rho + screened_rates(k_p_al27__he4_mg24)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jal27, scratch)

    scratch = (&
      -screened_rates(k_he4_si28__c12_ne20)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_si28__o16_o16)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_si28__p_p31)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_si28__s32)*Y(jhe4)*state % rho + &
      screened_rates(k_si28__he4_mg24) &
       )
    call set_jac_entry(jac, jhe4, jsi28, scratch)

    scratch = (&
      -screened_rates(k_he4_p31__cl35)*Y(jhe4)*state % rho + screened_rates(k_p31__he4_al27) + &
      screened_rates(k_p_p31__he4_si28)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jp31, scratch)

    scratch = (&
      -screened_rates(k_he4_s32__ar36)*Y(jhe4)*state % rho - screened_rates(k_he4_s32__p_cl35)* &
      Y(jhe4)*state % rho + screened_rates(k_s32__he4_si28) &
       )
    call set_jac_entry(jac, jhe4, js32, scratch)

    scratch = (&
      screened_rates(k_cl35__he4_p31) - screened_rates(k_he4_cl35__k39)*Y(jhe4)*state % rho + &
      screened_rates(k_p_cl35__he4_s32)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jcl35, scratch)

    scratch = (&
      screened_rates(k_ar36__he4_s32) - screened_rates(k_he4_ar36__ca40)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar36__p_k39)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jar36, scratch)

    scratch = (&
      -screened_rates(k_he4_k39__sc43)*Y(jhe4)*state % rho + screened_rates(k_k39__he4_cl35) + &
      screened_rates(k_p_k39__he4_ar36)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jk39, scratch)

    scratch = (&
      screened_rates(k_ca40__he4_ar36) - screened_rates(k_he4_ca40__p_sc43)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca40__ti44)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jca40, scratch)

    scratch = (&
      -screened_rates(k_he4_sc43__v47)*Y(jhe4)*state % rho + screened_rates(k_p_sc43__he4_ca40)* &
      Y(jp)*state % rho + screened_rates(k_sc43__he4_k39) &
       )
    call set_jac_entry(jac, jhe4, jsc43, scratch)

    scratch = (&
      -screened_rates(k_he4_ti44__cr48)*Y(jhe4)*state % rho - screened_rates(k_he4_ti44__p_v47)* &
      Y(jhe4)*state % rho + screened_rates(k_ti44__he4_ca40) &
       )
    call set_jac_entry(jac, jhe4, jti44, scratch)

    scratch = (&
      -screened_rates(k_he4_v47__mn51)*Y(jhe4)*state % rho + screened_rates(k_p_v47__he4_ti44)* &
      Y(jp)*state % rho + screened_rates(k_v47__he4_sc43) &
       )
    call set_jac_entry(jac, jhe4, jv47, scratch)

    scratch = (&
      screened_rates(k_cr48__he4_ti44) - screened_rates(k_he4_cr48__fe52)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr48__p_mn51)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jcr48, scratch)

    scratch = (&
      -screened_rates(k_he4_mn51__co55)*Y(jhe4)*state % rho + screened_rates(k_mn51__he4_v47) + &
      screened_rates(k_p_mn51__he4_cr48)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jmn51, scratch)

    scratch = (&
      screened_rates(k_fe52__he4_cr48) - screened_rates(k_he4_fe52__ni56)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe52__p_co55)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jfe52, scratch)

    scratch = (&
      screened_rates(k_co55__he4_mn51) + screened_rates(k_p_co55__he4_fe52)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jco55, scratch)

    scratch = (&
      screened_rates(k_ni56__he4_fe52) &
       )
    call set_jac_entry(jac, jhe4, jni56, scratch)

    scratch = (&
      screened_rates(k_p_al27__c12_o16)*Y(jal27)*state % rho - screened_rates(k_p_c12__n13)* &
      Y(jc12)*state % rho + screened_rates(k_p_p31__c12_ne20)*Y(jp31)*state % rho &
       )
    call set_jac_entry(jac, jc12, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho + 0.5d0* &
      screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**2*state % rho**2 + &
      screened_rates(k_he4_mg24__c12_o16)*Y(jmg24)*state % rho + 2.0d0* &
      screened_rates(k_he4_ne20__c12_c12)*Y(jne20)*state % rho + &
      screened_rates(k_he4_si28__c12_ne20)*Y(jsi28)*state % rho &
       )
    call set_jac_entry(jac, jc12, jhe4, scratch)

    scratch = (&
      -screened_rates(k_c12__he4_he4_he4) - 2.0d0*screened_rates(k_c12_c12__he4_ne20)*Y(jc12)* &
      state % rho - screened_rates(k_c12_ne20__he4_si28)*Y(jne20)*state % rho - &
      screened_rates(k_c12_ne20__p_p31)*Y(jne20)*state % rho - &
      screened_rates(k_c12_o16__he4_mg24)*Y(jo16)*state % rho - &
      screened_rates(k_c12_o16__p_al27)*Y(jo16)*state % rho - &
      screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho - screened_rates(k_p_c12__n13)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jc12, jc12, scratch)

    scratch = (&
      screened_rates(k_n13__p_c12) &
       )
    call set_jac_entry(jac, jc12, jn13, scratch)

    scratch = (&
      -screened_rates(k_c12_o16__he4_mg24)*Y(jc12)*state % rho - &
      screened_rates(k_c12_o16__p_al27)*Y(jc12)*state % rho + &
      screened_rates(k_o16__he4_c12) &
       )
    call set_jac_entry(jac, jc12, jo16, scratch)

    scratch = (&
      -screened_rates(k_c12_ne20__he4_si28)*Y(jc12)*state % rho - &
      screened_rates(k_c12_ne20__p_p31)*Y(jc12)*state % rho + 2.0d0* &
      screened_rates(k_he4_ne20__c12_c12)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jc12, jne20, scratch)

    scratch = (&
      screened_rates(k_he4_mg24__c12_o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jc12, jmg24, scratch)

    scratch = (&
      screened_rates(k_p_al27__c12_o16)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jc12, jal27, scratch)

    scratch = (&
      screened_rates(k_he4_si28__c12_ne20)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jc12, jsi28, scratch)

    scratch = (&
      screened_rates(k_p_p31__c12_ne20)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jc12, jp31, scratch)

    scratch = (&
      -screened_rates(k_he4_c14__o18)*Y(jc14)*state % rho &
       )
    call set_jac_entry(jac, jc14, jhe4, scratch)

    scratch = (&
      -screened_rates(k_c14__n14__weak__wc12) - screened_rates(k_he4_c14__o18)*Y(jhe4)* &
      state % rho &
       )
    call set_jac_entry(jac, jc14, jc14, scratch)

    scratch = (&
      screened_rates(k_o18__he4_c14) &
       )
    call set_jac_entry(jac, jc14, jo18, scratch)

    scratch = (&
      screened_rates(k_p_c12__n13)*Y(jc12)*state % rho + screened_rates(k_p_o16__he4_n13)* &
      Y(jo16)*state % rho &
       )
    call set_jac_entry(jac, jn13, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_n13__p_o16)*Y(jn13)*state % rho &
       )
    call set_jac_entry(jac, jn13, jhe4, scratch)

    scratch = (&
      screened_rates(k_p_c12__n13)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jn13, jc12, scratch)

    scratch = (&
      -screened_rates(k_he4_n13__p_o16)*Y(jhe4)*state % rho - screened_rates(k_n13__p_c12) &
       )
    call set_jac_entry(jac, jn13, jn13, scratch)

    scratch = (&
      screened_rates(k_p_o16__he4_n13)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jn13, jo16, scratch)

    scratch = (&
      -screened_rates(k_he4_n14__f18)*Y(jn14)*state % rho &
       )
    call set_jac_entry(jac, jn14, jhe4, scratch)

    scratch = (&
      screened_rates(k_c14__n14__weak__wc12) &
       )
    call set_jac_entry(jac, jn14, jc14, scratch)

    scratch = (&
      -screened_rates(k_he4_n14__f18)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jn14, jn14, scratch)

    scratch = (&
      screened_rates(k_f18__he4_n14) &
       )
    call set_jac_entry(jac, jn14, jf18, scratch)

    scratch = (&
      screened_rates(k_p_al27__c12_o16)*Y(jal27)*state % rho - screened_rates(k_p_o16__he4_n13)* &
      Y(jo16)*state % rho + 2.0d0*screened_rates(k_p_p31__o16_o16)*Y(jp31)* &
      state % rho &
       )
    call set_jac_entry(jac, jo16, jp, scratch)

    scratch = (&
      screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho + screened_rates(k_he4_mg24__c12_o16)* &
      Y(jmg24)*state % rho + screened_rates(k_he4_n13__p_o16)*Y(jn13)*state % rho &
      - screened_rates(k_he4_o16__ne20)*Y(jo16)*state % rho + 2.0d0* &
      screened_rates(k_he4_si28__o16_o16)*Y(jsi28)*state % rho &
       )
    call set_jac_entry(jac, jo16, jhe4, scratch)

    scratch = (&
      -screened_rates(k_c12_o16__he4_mg24)*Y(jo16)*state % rho - &
      screened_rates(k_c12_o16__p_al27)*Y(jo16)*state % rho + &
      screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jo16, jc12, scratch)

    scratch = (&
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jo16, jn13, scratch)

    scratch = (&
      -screened_rates(k_c12_o16__he4_mg24)*Y(jc12)*state % rho - &
      screened_rates(k_c12_o16__p_al27)*Y(jc12)*state % rho - &
      screened_rates(k_he4_o16__ne20)*Y(jhe4)*state % rho - screened_rates(k_o16__he4_c12) &
      - 2.0d0*screened_rates(k_o16_o16__he4_si28)*Y(jo16)*state % rho - 2.0d0* &
      screened_rates(k_o16_o16__p_p31)*Y(jo16)*state % rho - &
      screened_rates(k_p_o16__he4_n13)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jo16, jo16, scratch)

    scratch = (&
      screened_rates(k_ne20__he4_o16) &
       )
    call set_jac_entry(jac, jo16, jne20, scratch)

    scratch = (&
      screened_rates(k_he4_mg24__c12_o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jo16, jmg24, scratch)

    scratch = (&
      screened_rates(k_p_al27__c12_o16)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jo16, jal27, scratch)

    scratch = (&
      2.0d0*screened_rates(k_he4_si28__o16_o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jo16, jsi28, scratch)

    scratch = (&
      2.0d0*screened_rates(k_p_p31__o16_o16)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jo16, jp31, scratch)

    scratch = (&
      screened_rates(k_he4_c14__o18)*Y(jc14)*state % rho &
       )
    call set_jac_entry(jac, jo18, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_c14__o18)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jo18, jc14, scratch)

    scratch = (&
      -screened_rates(k_o18__he4_c14) &
       )
    call set_jac_entry(jac, jo18, jo18, scratch)

    scratch = (&
      screened_rates(k_f18__o18__weak__wc12) &
       )
    call set_jac_entry(jac, jo18, jf18, scratch)

    scratch = (&
      screened_rates(k_p_ne21__he4_f18)*Y(jne21)*state % rho &
       )
    call set_jac_entry(jac, jf18, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_f18__p_ne21)*Y(jf18)*state % rho + screened_rates(k_he4_n14__f18)* &
      Y(jn14)*state % rho &
       )
    call set_jac_entry(jac, jf18, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_n14__f18)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jf18, jn14, scratch)

    scratch = (&
      -screened_rates(k_f18__he4_n14) - screened_rates(k_f18__o18__weak__wc12) - &
      screened_rates(k_he4_f18__p_ne21)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jf18, jf18, scratch)

    scratch = (&
      screened_rates(k_p_ne21__he4_f18)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jf18, jne21, scratch)

    scratch = (&
      screened_rates(k_p_p31__c12_ne20)*Y(jp31)*state % rho &
       )
    call set_jac_entry(jac, jne20, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_ne20__c12_c12)*Y(jne20)*state % rho - &
      screened_rates(k_he4_ne20__mg24)*Y(jne20)*state % rho + &
      screened_rates(k_he4_o16__ne20)*Y(jo16)*state % rho + &
      screened_rates(k_he4_si28__c12_ne20)*Y(jsi28)*state % rho &
       )
    call set_jac_entry(jac, jne20, jhe4, scratch)

    scratch = (&
      1.0d0*screened_rates(k_c12_c12__he4_ne20)*Y(jc12)*state % rho - &
      screened_rates(k_c12_ne20__he4_si28)*Y(jne20)*state % rho - &
      screened_rates(k_c12_ne20__p_p31)*Y(jne20)*state % rho &
       )
    call set_jac_entry(jac, jne20, jc12, scratch)

    scratch = (&
      screened_rates(k_he4_o16__ne20)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jne20, jo16, scratch)

    scratch = (&
      -screened_rates(k_c12_ne20__he4_si28)*Y(jc12)*state % rho - &
      screened_rates(k_c12_ne20__p_p31)*Y(jc12)*state % rho - &
      screened_rates(k_he4_ne20__c12_c12)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ne20__mg24)*Y(jhe4)*state % rho - &
      screened_rates(k_ne20__he4_o16) &
       )
    call set_jac_entry(jac, jne20, jne20, scratch)

    scratch = (&
      screened_rates(k_mg24__he4_ne20) &
       )
    call set_jac_entry(jac, jne20, jmg24, scratch)

    scratch = (&
      screened_rates(k_he4_si28__c12_ne20)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jne20, jsi28, scratch)

    scratch = (&
      screened_rates(k_p_p31__c12_ne20)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jne20, jp31, scratch)

    scratch = (&
      -screened_rates(k_p_ne21__he4_f18)*Y(jne21)*state % rho &
       )
    call set_jac_entry(jac, jne21, jp, scratch)

    scratch = (&
      screened_rates(k_he4_f18__p_ne21)*Y(jf18)*state % rho &
       )
    call set_jac_entry(jac, jne21, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_f18__p_ne21)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jne21, jf18, scratch)

    scratch = (&
      -screened_rates(k_p_ne21__he4_f18)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jne21, jne21, scratch)

    scratch = (&
      screened_rates(k_p_al27__he4_mg24)*Y(jal27)*state % rho &
       )
    call set_jac_entry(jac, jmg24, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_mg24__c12_o16)*Y(jmg24)*state % rho - &
      screened_rates(k_he4_mg24__p_al27)*Y(jmg24)*state % rho - &
      screened_rates(k_he4_mg24__si28)*Y(jmg24)*state % rho + &
      screened_rates(k_he4_ne20__mg24)*Y(jne20)*state % rho &
       )
    call set_jac_entry(jac, jmg24, jhe4, scratch)

    scratch = (&
      screened_rates(k_c12_o16__he4_mg24)*Y(jo16)*state % rho &
       )
    call set_jac_entry(jac, jmg24, jc12, scratch)

    scratch = (&
      screened_rates(k_c12_o16__he4_mg24)*Y(jc12)*state % rho &
       )
    call set_jac_entry(jac, jmg24, jo16, scratch)

    scratch = (&
      screened_rates(k_he4_ne20__mg24)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jmg24, jne20, scratch)

    scratch = (&
      -screened_rates(k_he4_mg24__c12_o16)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_mg24__p_al27)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_mg24__si28)*Y(jhe4)*state % rho - &
      screened_rates(k_mg24__he4_ne20) &
       )
    call set_jac_entry(jac, jmg24, jmg24, scratch)

    scratch = (&
      screened_rates(k_p_al27__he4_mg24)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jmg24, jal27, scratch)

    scratch = (&
      screened_rates(k_si28__he4_mg24) &
       )
    call set_jac_entry(jac, jmg24, jsi28, scratch)

    scratch = (&
      -screened_rates(k_p_al27__c12_o16)*Y(jal27)*state % rho - &
      screened_rates(k_p_al27__he4_mg24)*Y(jal27)*state % rho - &
      screened_rates(k_p_al27__si28)*Y(jal27)*state % rho &
       )
    call set_jac_entry(jac, jal27, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_al27__p31)*Y(jal27)*state % rho + screened_rates(k_he4_mg24__p_al27)* &
      Y(jmg24)*state % rho &
       )
    call set_jac_entry(jac, jal27, jhe4, scratch)

    scratch = (&
      screened_rates(k_c12_o16__p_al27)*Y(jo16)*state % rho &
       )
    call set_jac_entry(jac, jal27, jc12, scratch)

    scratch = (&
      screened_rates(k_c12_o16__p_al27)*Y(jc12)*state % rho &
       )
    call set_jac_entry(jac, jal27, jo16, scratch)

    scratch = (&
      screened_rates(k_he4_mg24__p_al27)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jal27, jmg24, scratch)

    scratch = (&
      -screened_rates(k_he4_al27__p31)*Y(jhe4)*state % rho - screened_rates(k_p_al27__c12_o16)* &
      Y(jp)*state % rho - screened_rates(k_p_al27__he4_mg24)*Y(jp)*state % rho - &
      screened_rates(k_p_al27__si28)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jal27, jal27, scratch)

    scratch = (&
      screened_rates(k_si28__p_al27) &
       )
    call set_jac_entry(jac, jal27, jsi28, scratch)

    scratch = (&
      screened_rates(k_p31__he4_al27) &
       )
    call set_jac_entry(jac, jal27, jp31, scratch)

    scratch = (&
      screened_rates(k_p_al27__si28)*Y(jal27)*state % rho + screened_rates(k_p_p31__he4_si28)* &
      Y(jp31)*state % rho &
       )
    call set_jac_entry(jac, jsi28, jp, scratch)

    scratch = (&
      screened_rates(k_he4_mg24__si28)*Y(jmg24)*state % rho - &
      screened_rates(k_he4_si28__c12_ne20)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_si28__o16_o16)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_si28__p_p31)*Y(jsi28)*state % rho - &
      screened_rates(k_he4_si28__s32)*Y(jsi28)*state % rho &
       )
    call set_jac_entry(jac, jsi28, jhe4, scratch)

    scratch = (&
      screened_rates(k_c12_ne20__he4_si28)*Y(jne20)*state % rho &
       )
    call set_jac_entry(jac, jsi28, jc12, scratch)

    scratch = (&
      1.0d0*screened_rates(k_o16_o16__he4_si28)*Y(jo16)*state % rho &
       )
    call set_jac_entry(jac, jsi28, jo16, scratch)

    scratch = (&
      screened_rates(k_c12_ne20__he4_si28)*Y(jc12)*state % rho &
       )
    call set_jac_entry(jac, jsi28, jne20, scratch)

    scratch = (&
      screened_rates(k_he4_mg24__si28)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jsi28, jmg24, scratch)

    scratch = (&
      screened_rates(k_p_al27__si28)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jsi28, jal27, scratch)

    scratch = (&
      -screened_rates(k_he4_si28__c12_ne20)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_si28__o16_o16)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_si28__p_p31)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_si28__s32)*Y(jhe4)*state % rho - &
      screened_rates(k_si28__he4_mg24) - screened_rates(k_si28__p_al27) &
       )
    call set_jac_entry(jac, jsi28, jsi28, scratch)

    scratch = (&
      screened_rates(k_p_p31__he4_si28)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jsi28, jp31, scratch)

    scratch = (&
      screened_rates(k_s32__he4_si28) &
       )
    call set_jac_entry(jac, jsi28, js32, scratch)

    scratch = (&
      -screened_rates(k_p_p31__c12_ne20)*Y(jp31)*state % rho - screened_rates(k_p_p31__he4_si28)* &
      Y(jp31)*state % rho - screened_rates(k_p_p31__o16_o16)*Y(jp31)*state % rho - &
      screened_rates(k_p_p31__s32)*Y(jp31)*state % rho &
       )
    call set_jac_entry(jac, jp31, jp, scratch)

    scratch = (&
      screened_rates(k_he4_al27__p31)*Y(jal27)*state % rho - screened_rates(k_he4_p31__cl35)* &
      Y(jp31)*state % rho + screened_rates(k_he4_si28__p_p31)*Y(jsi28)*state % rho &
       )
    call set_jac_entry(jac, jp31, jhe4, scratch)

    scratch = (&
      screened_rates(k_c12_ne20__p_p31)*Y(jne20)*state % rho &
       )
    call set_jac_entry(jac, jp31, jc12, scratch)

    scratch = (&
      1.0d0*screened_rates(k_o16_o16__p_p31)*Y(jo16)*state % rho &
       )
    call set_jac_entry(jac, jp31, jo16, scratch)

    scratch = (&
      screened_rates(k_c12_ne20__p_p31)*Y(jc12)*state % rho &
       )
    call set_jac_entry(jac, jp31, jne20, scratch)

    scratch = (&
      screened_rates(k_he4_al27__p31)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jp31, jal27, scratch)

    scratch = (&
      screened_rates(k_he4_si28__p_p31)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jp31, jsi28, scratch)

    scratch = (&
      -screened_rates(k_he4_p31__cl35)*Y(jhe4)*state % rho - screened_rates(k_p31__he4_al27) - &
      screened_rates(k_p_p31__c12_ne20)*Y(jp)*state % rho - &
      screened_rates(k_p_p31__he4_si28)*Y(jp)*state % rho - &
      screened_rates(k_p_p31__o16_o16)*Y(jp)*state % rho - screened_rates(k_p_p31__s32)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jp31, jp31, scratch)

    scratch = (&
      screened_rates(k_s32__p_p31) &
       )
    call set_jac_entry(jac, jp31, js32, scratch)

    scratch = (&
      screened_rates(k_cl35__he4_p31) &
       )
    call set_jac_entry(jac, jp31, jcl35, scratch)

    scratch = (&
      screened_rates(k_p_cl35__he4_s32)*Y(jcl35)*state % rho + screened_rates(k_p_p31__s32)* &
      Y(jp31)*state % rho &
       )
    call set_jac_entry(jac, js32, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_s32__ar36)*Y(js32)*state % rho - screened_rates(k_he4_s32__p_cl35)* &
      Y(js32)*state % rho + screened_rates(k_he4_si28__s32)*Y(jsi28)*state % rho &
       )
    call set_jac_entry(jac, js32, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_si28__s32)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, js32, jsi28, scratch)

    scratch = (&
      screened_rates(k_p_p31__s32)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, js32, jp31, scratch)

    scratch = (&
      -screened_rates(k_he4_s32__ar36)*Y(jhe4)*state % rho - screened_rates(k_he4_s32__p_cl35)* &
      Y(jhe4)*state % rho - screened_rates(k_s32__he4_si28) - screened_rates(k_s32__p_p31) &
       )
    call set_jac_entry(jac, js32, js32, scratch)

    scratch = (&
      screened_rates(k_p_cl35__he4_s32)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, js32, jcl35, scratch)

    scratch = (&
      screened_rates(k_ar36__he4_s32) &
       )
    call set_jac_entry(jac, js32, jar36, scratch)

    scratch = (&
      -screened_rates(k_p_cl35__ar36)*Y(jcl35)*state % rho - screened_rates(k_p_cl35__he4_s32)* &
      Y(jcl35)*state % rho &
       )
    call set_jac_entry(jac, jcl35, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cl35__k39)*Y(jcl35)*state % rho + screened_rates(k_he4_p31__cl35)* &
      Y(jp31)*state % rho + screened_rates(k_he4_s32__p_cl35)*Y(js32)*state % rho &
       )
    call set_jac_entry(jac, jcl35, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_p31__cl35)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jcl35, jp31, scratch)

    scratch = (&
      screened_rates(k_he4_s32__p_cl35)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jcl35, js32, scratch)

    scratch = (&
      -screened_rates(k_cl35__he4_p31) - screened_rates(k_he4_cl35__k39)*Y(jhe4)*state % rho - &
      screened_rates(k_p_cl35__ar36)*Y(jp)*state % rho - screened_rates(k_p_cl35__he4_s32) &
      *Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jcl35, jcl35, scratch)

    scratch = (&
      screened_rates(k_ar36__p_cl35) &
       )
    call set_jac_entry(jac, jcl35, jar36, scratch)

    scratch = (&
      screened_rates(k_k39__he4_cl35) &
       )
    call set_jac_entry(jac, jcl35, jk39, scratch)

    scratch = (&
      screened_rates(k_p_cl35__ar36)*Y(jcl35)*state % rho + screened_rates(k_p_k39__he4_ar36)* &
      Y(jk39)*state % rho &
       )
    call set_jac_entry(jac, jar36, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_ar36__ca40)*Y(jar36)*state % rho - screened_rates(k_he4_ar36__p_k39)* &
      Y(jar36)*state % rho + screened_rates(k_he4_s32__ar36)*Y(js32)*state % rho &
       )
    call set_jac_entry(jac, jar36, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_s32__ar36)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jar36, js32, scratch)

    scratch = (&
      screened_rates(k_p_cl35__ar36)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jar36, jcl35, scratch)

    scratch = (&
      -screened_rates(k_ar36__he4_s32) - screened_rates(k_ar36__p_cl35) - &
      screened_rates(k_he4_ar36__ca40)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ar36__p_k39)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jar36, jar36, scratch)

    scratch = (&
      screened_rates(k_p_k39__he4_ar36)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jar36, jk39, scratch)

    scratch = (&
      screened_rates(k_ca40__he4_ar36) &
       )
    call set_jac_entry(jac, jar36, jca40, scratch)

    scratch = (&
      -screened_rates(k_p_k39__ca40)*Y(jk39)*state % rho - screened_rates(k_p_k39__he4_ar36)* &
      Y(jk39)*state % rho &
       )
    call set_jac_entry(jac, jk39, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar36__p_k39)*Y(jar36)*state % rho + screened_rates(k_he4_cl35__k39)* &
      Y(jcl35)*state % rho - screened_rates(k_he4_k39__sc43)*Y(jk39)*state % rho &
       )
    call set_jac_entry(jac, jk39, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cl35__k39)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jk39, jcl35, scratch)

    scratch = (&
      screened_rates(k_he4_ar36__p_k39)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jk39, jar36, scratch)

    scratch = (&
      -screened_rates(k_he4_k39__sc43)*Y(jhe4)*state % rho - screened_rates(k_k39__he4_cl35) - &
      screened_rates(k_p_k39__ca40)*Y(jp)*state % rho - screened_rates(k_p_k39__he4_ar36)* &
      Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jk39, jk39, scratch)

    scratch = (&
      screened_rates(k_ca40__p_k39) &
       )
    call set_jac_entry(jac, jk39, jca40, scratch)

    scratch = (&
      screened_rates(k_sc43__he4_k39) &
       )
    call set_jac_entry(jac, jk39, jsc43, scratch)

    scratch = (&
      screened_rates(k_p_k39__ca40)*Y(jk39)*state % rho + screened_rates(k_p_sc43__he4_ca40)* &
      Y(jsc43)*state % rho &
       )
    call set_jac_entry(jac, jca40, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ar36__ca40)*Y(jar36)*state % rho - screened_rates(k_he4_ca40__p_sc43)* &
      Y(jca40)*state % rho - screened_rates(k_he4_ca40__ti44)*Y(jca40)*state % rho &
       )
    call set_jac_entry(jac, jca40, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ar36__ca40)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jca40, jar36, scratch)

    scratch = (&
      screened_rates(k_p_k39__ca40)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jca40, jk39, scratch)

    scratch = (&
      -screened_rates(k_ca40__he4_ar36) - screened_rates(k_ca40__p_k39) - &
      screened_rates(k_he4_ca40__p_sc43)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_ca40__ti44)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jca40, jca40, scratch)

    scratch = (&
      screened_rates(k_p_sc43__he4_ca40)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jca40, jsc43, scratch)

    scratch = (&
      screened_rates(k_ti44__he4_ca40) &
       )
    call set_jac_entry(jac, jca40, jti44, scratch)

    scratch = (&
      -screened_rates(k_p_sc43__he4_ca40)*Y(jsc43)*state % rho - screened_rates(k_p_sc43__ti44)* &
      Y(jsc43)*state % rho &
       )
    call set_jac_entry(jac, jsc43, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca40__p_sc43)*Y(jca40)*state % rho + screened_rates(k_he4_k39__sc43)* &
      Y(jk39)*state % rho - screened_rates(k_he4_sc43__v47)*Y(jsc43)*state % rho &
       )
    call set_jac_entry(jac, jsc43, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_k39__sc43)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jsc43, jk39, scratch)

    scratch = (&
      screened_rates(k_he4_ca40__p_sc43)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jsc43, jca40, scratch)

    scratch = (&
      -screened_rates(k_he4_sc43__v47)*Y(jhe4)*state % rho - screened_rates(k_p_sc43__he4_ca40)* &
      Y(jp)*state % rho - screened_rates(k_p_sc43__ti44)*Y(jp)*state % rho - &
      screened_rates(k_sc43__he4_k39) &
       )
    call set_jac_entry(jac, jsc43, jsc43, scratch)

    scratch = (&
      screened_rates(k_ti44__p_sc43) &
       )
    call set_jac_entry(jac, jsc43, jti44, scratch)

    scratch = (&
      screened_rates(k_v47__he4_sc43) &
       )
    call set_jac_entry(jac, jsc43, jv47, scratch)

    scratch = (&
      screened_rates(k_p_sc43__ti44)*Y(jsc43)*state % rho + screened_rates(k_p_v47__he4_ti44)* &
      Y(jv47)*state % rho &
       )
    call set_jac_entry(jac, jti44, jp, scratch)

    scratch = (&
      screened_rates(k_he4_ca40__ti44)*Y(jca40)*state % rho - screened_rates(k_he4_ti44__cr48)* &
      Y(jti44)*state % rho - screened_rates(k_he4_ti44__p_v47)*Y(jti44)* &
      state % rho &
       )
    call set_jac_entry(jac, jti44, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ca40__ti44)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jti44, jca40, scratch)

    scratch = (&
      screened_rates(k_p_sc43__ti44)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jti44, jsc43, scratch)

    scratch = (&
      -screened_rates(k_he4_ti44__cr48)*Y(jhe4)*state % rho - screened_rates(k_he4_ti44__p_v47)* &
      Y(jhe4)*state % rho - screened_rates(k_ti44__he4_ca40) - &
      screened_rates(k_ti44__p_sc43) &
       )
    call set_jac_entry(jac, jti44, jti44, scratch)

    scratch = (&
      screened_rates(k_p_v47__he4_ti44)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jti44, jv47, scratch)

    scratch = (&
      screened_rates(k_cr48__he4_ti44) &
       )
    call set_jac_entry(jac, jti44, jcr48, scratch)

    scratch = (&
      -screened_rates(k_p_v47__cr48)*Y(jv47)*state % rho - screened_rates(k_p_v47__he4_ti44)* &
      Y(jv47)*state % rho &
       )
    call set_jac_entry(jac, jv47, jp, scratch)

    scratch = (&
      screened_rates(k_he4_sc43__v47)*Y(jsc43)*state % rho + screened_rates(k_he4_ti44__p_v47)* &
      Y(jti44)*state % rho - screened_rates(k_he4_v47__mn51)*Y(jv47)*state % rho &
       )
    call set_jac_entry(jac, jv47, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_sc43__v47)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jv47, jsc43, scratch)

    scratch = (&
      screened_rates(k_he4_ti44__p_v47)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jv47, jti44, scratch)

    scratch = (&
      -screened_rates(k_he4_v47__mn51)*Y(jhe4)*state % rho - screened_rates(k_p_v47__cr48)* &
      Y(jp)*state % rho - screened_rates(k_p_v47__he4_ti44)*Y(jp)*state % rho - &
      screened_rates(k_v47__he4_sc43) &
       )
    call set_jac_entry(jac, jv47, jv47, scratch)

    scratch = (&
      screened_rates(k_cr48__p_v47) &
       )
    call set_jac_entry(jac, jv47, jcr48, scratch)

    scratch = (&
      screened_rates(k_mn51__he4_v47) &
       )
    call set_jac_entry(jac, jv47, jmn51, scratch)

    scratch = (&
      screened_rates(k_p_mn51__he4_cr48)*Y(jmn51)*state % rho + screened_rates(k_p_v47__cr48)* &
      Y(jv47)*state % rho &
       )
    call set_jac_entry(jac, jcr48, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_cr48__fe52)*Y(jcr48)*state % rho - screened_rates(k_he4_cr48__p_mn51) &
      *Y(jcr48)*state % rho + screened_rates(k_he4_ti44__cr48)*Y(jti44)* &
      state % rho &
       )
    call set_jac_entry(jac, jcr48, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_ti44__cr48)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jcr48, jti44, scratch)

    scratch = (&
      screened_rates(k_p_v47__cr48)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jcr48, jv47, scratch)

    scratch = (&
      -screened_rates(k_cr48__he4_ti44) - screened_rates(k_cr48__p_v47) - &
      screened_rates(k_he4_cr48__fe52)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_cr48__p_mn51)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jcr48, jcr48, scratch)

    scratch = (&
      screened_rates(k_p_mn51__he4_cr48)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jcr48, jmn51, scratch)

    scratch = (&
      screened_rates(k_fe52__he4_cr48) &
       )
    call set_jac_entry(jac, jcr48, jfe52, scratch)

    scratch = (&
      -screened_rates(k_p_mn51__fe52)*Y(jmn51)*state % rho - screened_rates(k_p_mn51__he4_cr48)* &
      Y(jmn51)*state % rho &
       )
    call set_jac_entry(jac, jmn51, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr48__p_mn51)*Y(jcr48)*state % rho - screened_rates(k_he4_mn51__co55)* &
      Y(jmn51)*state % rho + screened_rates(k_he4_v47__mn51)*Y(jv47)*state % rho &
       )
    call set_jac_entry(jac, jmn51, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_v47__mn51)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jmn51, jv47, scratch)

    scratch = (&
      screened_rates(k_he4_cr48__p_mn51)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jmn51, jcr48, scratch)

    scratch = (&
      -screened_rates(k_he4_mn51__co55)*Y(jhe4)*state % rho - screened_rates(k_mn51__he4_v47) - &
      screened_rates(k_p_mn51__fe52)*Y(jp)*state % rho - &
      screened_rates(k_p_mn51__he4_cr48)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jmn51, jmn51, scratch)

    scratch = (&
      screened_rates(k_fe52__p_mn51) &
       )
    call set_jac_entry(jac, jmn51, jfe52, scratch)

    scratch = (&
      screened_rates(k_co55__he4_mn51) &
       )
    call set_jac_entry(jac, jmn51, jco55, scratch)

    scratch = (&
      screened_rates(k_p_co55__he4_fe52)*Y(jco55)*state % rho + screened_rates(k_p_mn51__fe52)* &
      Y(jmn51)*state % rho &
       )
    call set_jac_entry(jac, jfe52, jp, scratch)

    scratch = (&
      screened_rates(k_he4_cr48__fe52)*Y(jcr48)*state % rho - screened_rates(k_he4_fe52__ni56)* &
      Y(jfe52)*state % rho - screened_rates(k_he4_fe52__p_co55)*Y(jfe52)* &
      state % rho &
       )
    call set_jac_entry(jac, jfe52, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_cr48__fe52)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jfe52, jcr48, scratch)

    scratch = (&
      screened_rates(k_p_mn51__fe52)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jfe52, jmn51, scratch)

    scratch = (&
      -screened_rates(k_fe52__he4_cr48) - screened_rates(k_fe52__p_mn51) - &
      screened_rates(k_he4_fe52__ni56)*Y(jhe4)*state % rho - &
      screened_rates(k_he4_fe52__p_co55)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jfe52, jfe52, scratch)

    scratch = (&
      screened_rates(k_p_co55__he4_fe52)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jfe52, jco55, scratch)

    scratch = (&
      screened_rates(k_ni56__he4_fe52) &
       )
    call set_jac_entry(jac, jfe52, jni56, scratch)

    scratch = (&
      -screened_rates(k_p_co55__he4_fe52)*Y(jco55)*state % rho - screened_rates(k_p_co55__ni56)* &
      Y(jco55)*state % rho &
       )
    call set_jac_entry(jac, jco55, jp, scratch)

    scratch = (&
      screened_rates(k_he4_fe52__p_co55)*Y(jfe52)*state % rho + screened_rates(k_he4_mn51__co55)* &
      Y(jmn51)*state % rho &
       )
    call set_jac_entry(jac, jco55, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_mn51__co55)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jco55, jmn51, scratch)

    scratch = (&
      screened_rates(k_he4_fe52__p_co55)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jco55, jfe52, scratch)

    scratch = (&
      -screened_rates(k_co55__he4_mn51) - screened_rates(k_p_co55__he4_fe52)*Y(jp)*state % rho - &
      screened_rates(k_p_co55__ni56)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jco55, jco55, scratch)

    scratch = (&
      screened_rates(k_ni56__p_co55) &
       )
    call set_jac_entry(jac, jco55, jni56, scratch)

    scratch = (&
      screened_rates(k_p_co55__ni56)*Y(jco55)*state % rho &
       )
    call set_jac_entry(jac, jni56, jp, scratch)

    scratch = (&
      screened_rates(k_he4_fe52__ni56)*Y(jfe52)*state % rho &
       )
    call set_jac_entry(jac, jni56, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_fe52__ni56)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jni56, jfe52, scratch)

    scratch = (&
      screened_rates(k_p_co55__ni56)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jni56, jco55, scratch)

    scratch = (&
      -screened_rates(k_ni56__he4_fe52) - screened_rates(k_ni56__p_co55) &
       )
    call set_jac_entry(jac, jni56, jni56, scratch)


  end subroutine jac_nuc

end module actual_rhs_module
