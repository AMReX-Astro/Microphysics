module actual_rhs_module

  use microphysics_type_module, only: rt, ZERO, ONE
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
      rate_eval % unscreened_rates(i_scor,2) = scor
      rate_eval % unscreened_rates(i_dscor_dt,2) = dscor_dt


      call screen5(pstate, 2, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,3) = scor
      rate_eval % unscreened_rates(i_dscor_dt,3) = dscor_dt


      call screen5(pstate, 3, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,4) = scor
      rate_eval % unscreened_rates(i_dscor_dt,4) = dscor_dt
      rate_eval % unscreened_rates(i_scor,11) = scor
      rate_eval % unscreened_rates(i_dscor_dt,11) = dscor_dt


      call screen5(pstate, 4, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,5) = scor
      rate_eval % unscreened_rates(i_dscor_dt,5) = dscor_dt
      rate_eval % unscreened_rates(i_scor,12) = scor
      rate_eval % unscreened_rates(i_dscor_dt,12) = dscor_dt


      call screen5(pstate, 5, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,6) = scor
      rate_eval % unscreened_rates(i_dscor_dt,6) = dscor_dt


      call screen5(pstate, 6, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,7) = scor
      rate_eval % unscreened_rates(i_dscor_dt,7) = dscor_dt
      rate_eval % unscreened_rates(i_scor,13) = scor
      rate_eval % unscreened_rates(i_dscor_dt,13) = dscor_dt


      call screen5(pstate, 7, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,8) = scor
      rate_eval % unscreened_rates(i_dscor_dt,8) = dscor_dt
      rate_eval % unscreened_rates(i_scor,14) = scor
      rate_eval % unscreened_rates(i_dscor_dt,14) = dscor_dt


      call screen5(pstate, 8, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,9) = scor
      rate_eval % unscreened_rates(i_dscor_dt,9) = dscor_dt
      rate_eval % unscreened_rates(i_scor,10) = scor
      rate_eval % unscreened_rates(i_dscor_dt,10) = dscor_dt

    end if

    ! Calculate tabular rates
    call tabular_evaluate(rate_table_j_f20_o20, rhoy_table_j_f20_o20, temp_table_j_f20_o20, &
                          num_rhoy_j_f20_o20, num_temp_j_f20_o20, num_vars_j_f20_o20, &
                          rhoy, state % T, reactvec)
    rate_eval % unscreened_rates(:,15) = reactvec(1:4)
    rate_eval % add_energy(1) = reactvec(5)
    rate_eval % add_energy_rate(1)  = reactvec(6)

    call tabular_evaluate(rate_table_j_ne20_f20, rhoy_table_j_ne20_f20, temp_table_j_ne20_f20, &
                          num_rhoy_j_ne20_f20, num_temp_j_ne20_f20, num_vars_j_ne20_f20, &
                          rhoy, state % T, reactvec)
    rate_eval % unscreened_rates(:,16) = reactvec(1:4)
    rate_eval % add_energy(2) = reactvec(5)
    rate_eval % add_energy_rate(2)  = reactvec(6)

    call tabular_evaluate(rate_table_j_o20_f20, rhoy_table_j_o20_f20, temp_table_j_o20_f20, &
                          num_rhoy_j_o20_f20, num_temp_j_o20_f20, num_vars_j_o20_f20, &
                          rhoy, state % T, reactvec)
    rate_eval % unscreened_rates(:,17) = reactvec(1:4)
    rate_eval % add_energy(3) = reactvec(5)
    rate_eval % add_energy_rate(3)  = reactvec(6)

    call tabular_evaluate(rate_table_j_f20_ne20, rhoy_table_j_f20_ne20, temp_table_j_f20_ne20, &
                          num_rhoy_j_f20_ne20, num_temp_j_f20_ne20, num_vars_j_f20_ne20, &
                          rhoy, state % T, reactvec)
    rate_eval % unscreened_rates(:,18) = reactvec(1:4)
    rate_eval % add_energy(4) = reactvec(5)
    rate_eval % add_energy_rate(4)  = reactvec(6)


    ! Compute screened rates
    rate_eval % screened_rates = rate_eval % unscreened_rates(i_rate, :) * &
                                 rate_eval % unscreened_rates(i_scor, :)

  end subroutine evaluate_rates


  subroutine actual_rhs(state)
    
    !$acc routine seq

    use extern_probin_module, only: do_constant_volume_burn, disable_thermal_neutrinos
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
    enuc = enuc + N_AVO * state % ydot(jf20) * rate_eval % add_energy(j_f20_o20)
    enuc = enuc + N_AVO * state % ydot(jne20) * rate_eval % add_energy(j_ne20_f20)
    enuc = enuc + N_AVO * state % ydot(jo20) * rate_eval % add_energy(j_o20_f20)
    enuc = enuc + N_AVO * state % ydot(jf20) * rate_eval % add_energy(j_f20_ne20)

    ! additional energy generation rates
    ! including gamma heating and reaction neutrino losses (non-thermal)
    enuc = enuc + N_AVO * Y(jf20) * rate_eval % add_energy_rate(j_f20_o20)
    enuc = enuc + N_AVO * Y(jne20) * rate_eval % add_energy_rate(j_ne20_f20)
    enuc = enuc + N_AVO * Y(jo20) * rate_eval % add_energy_rate(j_o20_f20)
    enuc = enuc + N_AVO * Y(jf20) * rate_eval % add_energy_rate(j_f20_ne20)


    ! Get the thermal neutrino losses
    if (.not. disable_thermal_neutrinos) then
       call sneut5(state % T, state % rho, state % abar, state % zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)
    else
       sneut = 0.0e0_rt
    end if

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

    real(rt) :: scratch_0
    real(rt) :: scratch_1
    real(rt) :: scratch_2
    real(rt) :: scratch_3
    real(rt) :: scratch_4
    real(rt) :: scratch_5
    real(rt) :: scratch_6
    real(rt) :: scratch_7
    real(rt) :: scratch_8
    real(rt) :: scratch_9
    real(rt) :: scratch_10
    real(rt) :: scratch_11
    real(rt) :: scratch_12
    real(rt) :: scratch_13
    real(rt) :: scratch_14
    real(rt) :: scratch_15
    real(rt) :: scratch_16
    real(rt) :: scratch_17
    real(rt) :: scratch_18
    real(rt) :: scratch_19
    real(rt) :: scratch_20
    real(rt) :: scratch_21
    real(rt) :: scratch_22
    real(rt) :: scratch_23
    real(rt) :: scratch_24
    real(rt) :: scratch_25
    real(rt) :: scratch_26
    real(rt) :: scratch_27
    real(rt) :: scratch_28
    real(rt) :: scratch_29
    real(rt) :: scratch_30
    real(rt) :: scratch_31

    scratch_0 = Y(jhe4)*Y(jmg24)*state % rho
    scratch_1 = screened_rates(k_he4_mg24__p_al27)*scratch_0
    scratch_2 = Y(jal27)*Y(jp)*state % rho
    scratch_3 = screened_rates(k_p_al27__he4_mg24)*scratch_2
    scratch_4 = screened_rates(k_p_al27__si28)*scratch_2
    scratch_5 = scratch_1 - scratch_3 - scratch_4
    scratch_6 = Y(jo16)**2*state % rho
    scratch_7 = screened_rates(k_o16_o16__p_p31)*scratch_6
    scratch_8 = Y(jhe4)*Y(jsi28)*state % rho
    scratch_9 = screened_rates(k_he4_si28__p_p31)*scratch_8
    scratch_10 = Y(jp31)*Y(jp)*state % rho
    scratch_11 = screened_rates(k_p_p31__he4_si28)*scratch_10
    scratch_12 = screened_rates(k_p_p31__s32)*scratch_10
    scratch_13 = -scratch_11 - scratch_12 + 0.5e0_rt*scratch_7 + scratch_9
    scratch_14 = screened_rates(k_ne20__he4_o16)*Y(jne20)
    scratch_15 = Y(jhe4)*state % rho
    scratch_16 = screened_rates(k_he4_al27__p31)*Y(jal27)*scratch_15
    scratch_17 = -scratch_16
    scratch_18 = -scratch_1
    scratch_19 = screened_rates(k_he4_mg24__si28)*scratch_0
    scratch_20 = -scratch_19
    scratch_21 = screened_rates(k_he4_ne20__mg24)*Y(jne20)*scratch_15
    scratch_22 = -scratch_21
    scratch_23 = screened_rates(k_he4_o16__ne20)*Y(jo16)*scratch_15
    scratch_24 = -scratch_23
    scratch_25 = screened_rates(k_o16_o16__he4_si28)*scratch_6
    scratch_26 = screened_rates(k_he4_si28__s32)*scratch_8
    scratch_27 = scratch_11 + 0.5e0_rt*scratch_25 - scratch_26 - scratch_9
    scratch_28 = screened_rates(k_f20__o20)*Y(jf20)
    scratch_29 = screened_rates(k_o20__f20)*Y(jo20)
    scratch_30 = screened_rates(k_ne20__f20)*Y(jne20)
    scratch_31 = screened_rates(k_f20__ne20)*Y(jf20)

    ydot_nuc(jp) = ( &
      scratch_13 + scratch_5 &
       )

    ydot_nuc(jhe4) = ( &
      scratch_14 + scratch_17 + scratch_18 + scratch_20 + scratch_22 + &
      scratch_24 + scratch_27 + scratch_3 &
       )

    ydot_nuc(jo16) = ( &
      scratch_14 + scratch_24 - scratch_25 - scratch_7 &
       )

    ydot_nuc(jo20) = ( &
      scratch_28 - scratch_29 &
       )

    ydot_nuc(jf20) = ( &
      -scratch_28 + scratch_29 + scratch_30 - scratch_31 &
       )

    ydot_nuc(jne20) = ( &
      -scratch_14 + scratch_22 + scratch_23 - scratch_30 + scratch_31 &
       )

    ydot_nuc(jmg24) = ( &
      scratch_18 + scratch_20 + scratch_21 + scratch_3 &
       )

    ydot_nuc(jal27) = ( &
      scratch_17 + scratch_5 &
       )

    ydot_nuc(jsi28) = ( &
      scratch_19 + scratch_27 + scratch_4 &
       )

    ydot_nuc(jp31) = ( &
      scratch_13 + scratch_16 &
       )

    ydot_nuc(js32) = ( &
      scratch_12 + scratch_26 &
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
       b1 = (-state % abar * state % abar * snuda + (zion(j) - state % zbar) * state % abar * snudz)
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

    real(rt) :: scratch_0
    real(rt) :: scratch_1
    real(rt) :: scratch_2
    real(rt) :: scratch_3
    real(rt) :: scratch_4
    real(rt) :: scratch_5
    real(rt) :: scratch_6
    real(rt) :: scratch_7
    real(rt) :: scratch_8
    real(rt) :: scratch_9
    real(rt) :: scratch_10
    real(rt) :: scratch_11
    real(rt) :: scratch_12
    real(rt) :: scratch_13
    real(rt) :: scratch_14
    real(rt) :: scratch_15
    real(rt) :: scratch_16
    real(rt) :: scratch_17
    real(rt) :: scratch_18
    real(rt) :: scratch_19
    real(rt) :: scratch_20
    real(rt) :: scratch_21
    real(rt) :: scratch_22
    real(rt) :: scratch_23
    real(rt) :: scratch_24
    real(rt) :: scratch_25
    real(rt) :: scratch_26
    real(rt) :: scratch_27
    real(rt) :: scratch_28
    real(rt) :: scratch_29
    real(rt) :: scratch_30
    real(rt) :: scratch_31
    real(rt) :: scratch_32
    real(rt) :: scratch_33
    real(rt) :: scratch_34
    real(rt) :: scratch_35
    real(rt) :: scratch_36
    real(rt) :: scratch_37
    real(rt) :: scratch_38
    real(rt) :: scratch_39
    real(rt) :: scratch_40
    real(rt) :: scratch_41
    real(rt) :: scratch_42
    real(rt) :: scratch_43
    real(rt) :: scratch_44
    real(rt) :: scratch_45
    real(rt) :: scratch_46

    !$gpu

    scratch_0 = Y(jp31)*state % rho
    scratch_1 = screened_rates(k_p_p31__he4_si28)*scratch_0
    scratch_2 = screened_rates(k_p_p31__s32)*scratch_0
    scratch_3 = -scratch_1 - scratch_2
    scratch_4 = Y(jal27)*state % rho
    scratch_5 = screened_rates(k_p_al27__he4_mg24)*scratch_4
    scratch_6 = screened_rates(k_p_al27__si28)*scratch_4
    scratch_7 = -scratch_5 - scratch_6
    scratch_8 = Y(jmg24)*state % rho
    scratch_9 = screened_rates(k_he4_mg24__p_al27)*scratch_8
    scratch_10 = Y(jsi28)*state % rho
    scratch_11 = screened_rates(k_he4_si28__p_p31)*scratch_10
    scratch_12 = screened_rates(k_o16_o16__p_p31)*Y(jo16)*state % rho
    scratch_13 = 1.0e0_rt*scratch_12
    scratch_14 = Y(jhe4)*state % rho
    scratch_15 = screened_rates(k_he4_mg24__p_al27)*scratch_14
    scratch_16 = Y(jp)*state % rho
    scratch_17 = screened_rates(k_p_al27__he4_mg24)*scratch_16
    scratch_18 = screened_rates(k_p_al27__si28)*scratch_16
    scratch_19 = -scratch_17 - scratch_18
    scratch_20 = screened_rates(k_he4_si28__p_p31)*scratch_14
    scratch_21 = screened_rates(k_p_p31__he4_si28)*scratch_16
    scratch_22 = screened_rates(k_p_p31__s32)*scratch_16
    scratch_23 = -scratch_21 - scratch_22
    scratch_24 = screened_rates(k_he4_al27__p31)*scratch_4
    scratch_25 = -scratch_24
    scratch_26 = -scratch_9
    scratch_27 = screened_rates(k_he4_mg24__si28)*scratch_8
    scratch_28 = -scratch_27
    scratch_29 = screened_rates(k_he4_ne20__mg24)*Y(jne20)*state % rho
    scratch_30 = -scratch_29
    scratch_31 = screened_rates(k_he4_o16__ne20)*Y(jo16)*state % rho
    scratch_32 = -scratch_31
    scratch_33 = screened_rates(k_he4_si28__s32)*scratch_10
    scratch_34 = -scratch_11 - scratch_33
    scratch_35 = screened_rates(k_he4_o16__ne20)*scratch_14
    scratch_36 = -scratch_35
    scratch_37 = screened_rates(k_o16_o16__he4_si28)*Y(jo16)*state % rho
    scratch_38 = 1.0e0_rt*scratch_37
    scratch_39 = screened_rates(k_he4_ne20__mg24)*scratch_14
    scratch_40 = -scratch_39
    scratch_41 = screened_rates(k_he4_mg24__si28)*scratch_14
    scratch_42 = -scratch_15 - scratch_41
    scratch_43 = screened_rates(k_he4_al27__p31)*scratch_14
    scratch_44 = -scratch_43
    scratch_45 = screened_rates(k_he4_si28__s32)*scratch_14
    scratch_46 = -scratch_20 - scratch_45

    scratch = (&
      scratch_3 + scratch_7 &
       )
    call set_jac_entry(state, jp, jp, scratch)

    scratch = (&
      scratch_11 + scratch_9 &
       )
    call set_jac_entry(state, jp, jhe4, scratch)

    scratch = (&
      scratch_13 &
       )
    call set_jac_entry(state, jp, jo16, scratch)

    scratch = (&
      scratch_15 &
       )
    call set_jac_entry(state, jp, jmg24, scratch)

    scratch = (&
      scratch_19 &
       )
    call set_jac_entry(state, jp, jal27, scratch)

    scratch = (&
      scratch_20 &
       )
    call set_jac_entry(state, jp, jsi28, scratch)

    scratch = (&
      scratch_23 &
       )
    call set_jac_entry(state, jp, jp31, scratch)

    scratch = (&
      scratch_1 + scratch_5 &
       )
    call set_jac_entry(state, jhe4, jp, scratch)

    scratch = (&
      scratch_25 + scratch_26 + scratch_28 + scratch_30 + scratch_32 + &
      scratch_34 &
       )
    call set_jac_entry(state, jhe4, jhe4, scratch)

    scratch = (&
      scratch_36 + scratch_38 &
       )
    call set_jac_entry(state, jhe4, jo16, scratch)

    scratch = (&
      screened_rates(k_ne20__he4_o16) + scratch_40 &
       )
    call set_jac_entry(state, jhe4, jne20, scratch)

    scratch = (&
      scratch_42 &
       )
    call set_jac_entry(state, jhe4, jmg24, scratch)

    scratch = (&
      scratch_17 + scratch_44 &
       )
    call set_jac_entry(state, jhe4, jal27, scratch)

    scratch = (&
      scratch_46 &
       )
    call set_jac_entry(state, jhe4, jsi28, scratch)

    scratch = (&
      scratch_21 &
       )
    call set_jac_entry(state, jhe4, jp31, scratch)

    scratch = (&
      scratch_32 &
       )
    call set_jac_entry(state, jo16, jhe4, scratch)

    scratch = (&
      -2.0e0_rt*scratch_12 + scratch_36 - 2.0e0_rt*scratch_37 &
       )
    call set_jac_entry(state, jo16, jo16, scratch)

    scratch = (&
      screened_rates(k_ne20__he4_o16) &
       )
    call set_jac_entry(state, jo16, jne20, scratch)

    scratch = (&
      -screened_rates(k_o20__f20) &
       )
    call set_jac_entry(state, jo20, jo20, scratch)

    scratch = (&
      screened_rates(k_f20__o20) &
       )
    call set_jac_entry(state, jo20, jf20, scratch)

    scratch = (&
      screened_rates(k_o20__f20) &
       )
    call set_jac_entry(state, jf20, jo20, scratch)

    scratch = (&
      -screened_rates(k_f20__ne20) - screened_rates(k_f20__o20) &
       )
    call set_jac_entry(state, jf20, jf20, scratch)

    scratch = (&
      screened_rates(k_ne20__f20) &
       )
    call set_jac_entry(state, jf20, jne20, scratch)

    scratch = (&
      scratch_30 + scratch_31 &
       )
    call set_jac_entry(state, jne20, jhe4, scratch)

    scratch = (&
      scratch_35 &
       )
    call set_jac_entry(state, jne20, jo16, scratch)

    scratch = (&
      screened_rates(k_f20__ne20) &
       )
    call set_jac_entry(state, jne20, jf20, scratch)

    scratch = (&
      -screened_rates(k_ne20__f20) - screened_rates(k_ne20__he4_o16) + scratch_40 &
       )
    call set_jac_entry(state, jne20, jne20, scratch)

    scratch = (&
      scratch_5 &
       )
    call set_jac_entry(state, jmg24, jp, scratch)

    scratch = (&
      scratch_26 + scratch_28 + scratch_29 &
       )
    call set_jac_entry(state, jmg24, jhe4, scratch)

    scratch = (&
      scratch_39 &
       )
    call set_jac_entry(state, jmg24, jne20, scratch)

    scratch = (&
      scratch_42 &
       )
    call set_jac_entry(state, jmg24, jmg24, scratch)

    scratch = (&
      scratch_17 &
       )
    call set_jac_entry(state, jmg24, jal27, scratch)

    scratch = (&
      scratch_7 &
       )
    call set_jac_entry(state, jal27, jp, scratch)

    scratch = (&
      scratch_25 + scratch_9 &
       )
    call set_jac_entry(state, jal27, jhe4, scratch)

    scratch = (&
      scratch_15 &
       )
    call set_jac_entry(state, jal27, jmg24, scratch)

    scratch = (&
      scratch_19 + scratch_44 &
       )
    call set_jac_entry(state, jal27, jal27, scratch)

    scratch = (&
      scratch_1 + scratch_6 &
       )
    call set_jac_entry(state, jsi28, jp, scratch)

    scratch = (&
      scratch_27 + scratch_34 &
       )
    call set_jac_entry(state, jsi28, jhe4, scratch)

    scratch = (&
      scratch_38 &
       )
    call set_jac_entry(state, jsi28, jo16, scratch)

    scratch = (&
      scratch_41 &
       )
    call set_jac_entry(state, jsi28, jmg24, scratch)

    scratch = (&
      scratch_18 &
       )
    call set_jac_entry(state, jsi28, jal27, scratch)

    scratch = (&
      scratch_46 &
       )
    call set_jac_entry(state, jsi28, jsi28, scratch)

    scratch = (&
      scratch_21 &
       )
    call set_jac_entry(state, jsi28, jp31, scratch)

    scratch = (&
      scratch_3 &
       )
    call set_jac_entry(state, jp31, jp, scratch)

    scratch = (&
      scratch_11 + scratch_24 &
       )
    call set_jac_entry(state, jp31, jhe4, scratch)

    scratch = (&
      scratch_13 &
       )
    call set_jac_entry(state, jp31, jo16, scratch)

    scratch = (&
      scratch_43 &
       )
    call set_jac_entry(state, jp31, jal27, scratch)

    scratch = (&
      scratch_20 &
       )
    call set_jac_entry(state, jp31, jsi28, scratch)

    scratch = (&
      scratch_23 &
       )
    call set_jac_entry(state, jp31, jp31, scratch)

    scratch = (&
      scratch_2 &
       )
    call set_jac_entry(state, js32, jp, scratch)

    scratch = (&
      scratch_33 &
       )
    call set_jac_entry(state, js32, jhe4, scratch)

    scratch = (&
      scratch_45 &
       )
    call set_jac_entry(state, js32, jsi28, scratch)

    scratch = (&
      scratch_22 &
       )
    call set_jac_entry(state, js32, jp31, scratch)


  end subroutine jac_nuc

end module actual_rhs_module
