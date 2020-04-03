module actual_rhs_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module
  use physical_constants, only: N_AVO
  use network
  use table_rates
  use burn_type_module

  implicit none

  ! Indices into rate groups in the rate_eval_t type
  integer, parameter :: i_rate        = 1
  integer, parameter :: i_drate_dt    = 2
  integer, parameter :: i_scor        = 3
  integer, parameter :: i_dscor_dt    = 4

  type :: rate_eval_t
     real(rt) :: unscreened_rates(num_rate_groups, nrates)
     real(rt) :: screened_rates(nrates)
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
    integer :: i, j
    real(rt) :: rhoy
    real(rt) :: rate, drate_dt, edot_nu
    real(rt) :: scor, dscor_dt, dscor_dd

    !$gpu

    Y(:) = state % xn(:) * aion_inv(:)
    rhoy = state % rho * state % y_e

    ! Zero out the rates
    call zero_rate_eval(rate_eval)

    ! Calculate Reaclib rates
    call fill_plasma_state(pstate, state % T, state % rho, Y)
    do i = 1, nrat_reaclib
       call reaclib_evaluate(pstate, state % T, i, rate, drate_dt)
       rate_eval % unscreened_rates(i_rate, i) = rate
       rate_eval % unscreened_rates(i_drate_dt, i) = drate_dt
    end do

    ! Evaluate screening factors
    if (screen_reaclib) then

      call screen5(pstate, 1, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,5) = scor
      rate_eval % unscreened_rates(i_dscor_dt,5) = dscor_dt


      call screen5(pstate, 2, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,6) = scor
      rate_eval % unscreened_rates(i_dscor_dt,6) = dscor_dt


      call screen5(pstate, 3, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,7) = scor
      rate_eval % unscreened_rates(i_dscor_dt,7) = dscor_dt


      call screen5(pstate, 4, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,8) = scor
      rate_eval % unscreened_rates(i_dscor_dt,8) = dscor_dt


      call screen5(pstate, 5, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,9) = scor
      rate_eval % unscreened_rates(i_dscor_dt,9) = dscor_dt


      call screen5(pstate, 6, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,10) = scor
      rate_eval % unscreened_rates(i_dscor_dt,10) = dscor_dt


      call screen5(pstate, 7, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,11) = scor
      rate_eval % unscreened_rates(i_dscor_dt,11) = dscor_dt
      rate_eval % unscreened_rates(i_scor,15) = scor
      rate_eval % unscreened_rates(i_dscor_dt,15) = dscor_dt


      call screen5(pstate, 8, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,12) = scor
      rate_eval % unscreened_rates(i_dscor_dt,12) = dscor_dt


      call screen5(pstate, 9, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,13) = scor
      rate_eval % unscreened_rates(i_dscor_dt,13) = dscor_dt
      rate_eval % unscreened_rates(i_scor,17) = scor
      rate_eval % unscreened_rates(i_dscor_dt,17) = dscor_dt


      call screen5(pstate, 10, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,14) = scor
      rate_eval % unscreened_rates(i_dscor_dt,14) = dscor_dt


      call screen5(pstate, 11, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,16) = scor
      rate_eval % unscreened_rates(i_dscor_dt,16) = dscor_dt


      call screen5(pstate, 12, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,18) = scor
      rate_eval % unscreened_rates(i_dscor_dt,18) = dscor_dt


      call screen5(pstate, 13, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,19) = scor
      rate_eval % unscreened_rates(i_dscor_dt,19) = dscor_dt

    end if


    ! Compute screened rates
    rate_eval % screened_rates = rate_eval % unscreened_rates(i_rate, :) * &
                                 rate_eval % unscreened_rates(i_scor, :)

  end subroutine evaluate_rates


  subroutine actual_rhs(state, ydot)
    
    !$acc routine seq

    use extern_probin_module, only: do_constant_volume_burn, disable_thermal_neutrinos
    use burn_type_module, only: net_itemp, net_ienuc, neqs
    use sneut_module, only: sneut5
    use temperature_integration_module, only: temperature_rhs

    implicit none

    type(burn_t), intent(in) :: state
    real(rt), intent(inout) :: ydot(neqs)

    type(rate_eval_t) :: rate_eval
    real(rt) :: Y(nspec), ydot_nuc(nspec)
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

    ! include reaction neutrino losses (non-thermal)

    ! Get the thermal neutrino losses
    if (.not. disable_thermal_neutrinos) then
       call sneut5(state % T, state % rho, state % abar, state % zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)
    else
       sneut = ZERO
    end if

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
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*Y(jn13)*state % rho + &
      screened_rates(k_he4_o14__p_f17)*Y(jhe4)*Y(jo14)*state % rho - &
      screened_rates(k_p_c12__n13)*Y(jc12)*Y(jp)*state % rho - &
      screened_rates(k_p_c13__n14)*Y(jc13)*Y(jp)*state % rho - &
      screened_rates(k_p_f18__he4_o15)*Y(jf18)*Y(jp)*state % rho - &
      screened_rates(k_p_n13__o14)*Y(jn13)*Y(jp)*state % rho - &
      screened_rates(k_p_n14__o15)*Y(jn14)*Y(jp)*state % rho - &
      screened_rates(k_p_n15__he4_c12)*Y(jn15)*Y(jp)*state % rho - &
      screened_rates(k_p_n15__o16)*Y(jn15)*Y(jp)*state % rho - &
      screened_rates(k_p_o16__f17)*Y(jo16)*Y(jp)*state % rho - &
      screened_rates(k_p_o17__f18)*Y(jo17)*Y(jp)*state % rho - &
      screened_rates(k_p_o17__he4_n14)*Y(jo17)*Y(jp)*state % rho &
       )

    ydot_nuc(jhe4) = ( &
      -screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*state % rho - 0.5e0_rt* &
      screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**3*state % rho**2 - &
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*Y(jn13)*state % rho - &
      screened_rates(k_he4_n14__f18)*Y(jhe4)*Y(jn14)*state % rho - &
      screened_rates(k_he4_o14__p_f17)*Y(jhe4)*Y(jo14)*state % rho + &
      screened_rates(k_p_f18__he4_o15)*Y(jf18)*Y(jp)*state % rho + &
      screened_rates(k_p_n15__he4_c12)*Y(jn15)*Y(jp)*state % rho + &
      screened_rates(k_p_o17__he4_n14)*Y(jo17)*Y(jp)*state % rho &
       )

    ydot_nuc(jc12) = ( &
      -screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*state % rho + &
      0.16666666666666667e0_rt*screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**3* &
      state % rho**2 - screened_rates(k_p_c12__n13)*Y(jc12)*Y(jp)*state % rho + &
      screened_rates(k_p_n15__he4_c12)*Y(jn15)*Y(jp)*state % rho &
       )

    ydot_nuc(jc13) = ( &
      screened_rates(k_n13__c13__weak__wc12)*Y(jn13) - screened_rates(k_p_c13__n14)*Y(jc13)* &
      Y(jp)*state % rho &
       )

    ydot_nuc(jn13) = ( &
      -screened_rates(k_he4_n13__p_o16)*Y(jhe4)*Y(jn13)*state % rho - &
      screened_rates(k_n13__c13__weak__wc12)*Y(jn13) + screened_rates(k_p_c12__n13)* &
      Y(jc12)*Y(jp)*state % rho - screened_rates(k_p_n13__o14)*Y(jn13)* &
      Y(jp)*state % rho &
       )

    ydot_nuc(jn14) = ( &
      -screened_rates(k_he4_n14__f18)*Y(jhe4)*Y(jn14)*state % rho + &
      screened_rates(k_o14__n14__weak__wc12)*Y(jo14) + screened_rates(k_p_c13__n14)* &
      Y(jc13)*Y(jp)*state % rho - screened_rates(k_p_n14__o15)*Y(jn14)* &
      Y(jp)*state % rho + screened_rates(k_p_o17__he4_n14)*Y(jo17)*Y(jp)* &
      state % rho &
       )

    ydot_nuc(jn15) = ( &
      screened_rates(k_o15__n15__weak__wc12)*Y(jo15) - screened_rates(k_p_n15__he4_c12)* &
      Y(jn15)*Y(jp)*state % rho - screened_rates(k_p_n15__o16)*Y(jn15)* &
      Y(jp)*state % rho &
       )

    ydot_nuc(jo14) = ( &
      -screened_rates(k_he4_o14__p_f17)*Y(jhe4)*Y(jo14)*state % rho - &
      screened_rates(k_o14__n14__weak__wc12)*Y(jo14) + screened_rates(k_p_n13__o14)* &
      Y(jn13)*Y(jp)*state % rho &
       )

    ydot_nuc(jo15) = ( &
      -screened_rates(k_o15__n15__weak__wc12)*Y(jo15) + screened_rates(k_p_f18__he4_o15)* &
      Y(jf18)*Y(jp)*state % rho + screened_rates(k_p_n14__o15)*Y(jn14)* &
      Y(jp)*state % rho &
       )

    ydot_nuc(jo16) = ( &
      screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*state % rho + &
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*Y(jn13)*state % rho + &
      screened_rates(k_p_n15__o16)*Y(jn15)*Y(jp)*state % rho - &
      screened_rates(k_p_o16__f17)*Y(jo16)*Y(jp)*state % rho &
       )

    ydot_nuc(jo17) = ( &
      screened_rates(k_f17__o17__weak__wc12)*Y(jf17) - screened_rates(k_p_o17__f18)*Y(jo17)* &
      Y(jp)*state % rho - screened_rates(k_p_o17__he4_n14)*Y(jo17)*Y(jp)* &
      state % rho &
       )

    ydot_nuc(jf17) = ( &
      -screened_rates(k_f17__o17__weak__wc12)*Y(jf17) + screened_rates(k_he4_o14__p_f17)* &
      Y(jhe4)*Y(jo14)*state % rho + screened_rates(k_p_o16__f17)*Y(jo16)* &
      Y(jp)*state % rho &
       )

    ydot_nuc(jf18) = ( &
      screened_rates(k_he4_n14__f18)*Y(jhe4)*Y(jn14)*state % rho - &
      screened_rates(k_p_f18__he4_o15)*Y(jf18)*Y(jp)*state % rho + &
      screened_rates(k_p_o17__f18)*Y(jo17)*Y(jp)*state % rho &
       )


  end subroutine rhs_nuc


  subroutine actual_jac(state, jac)

    !$acc routine seq

    use burn_type_module, only: net_itemp, net_ienuc, neqs, njrows, njcols
    use extern_probin_module, only: disable_thermal_neutrinos
    use sneut_module, only: sneut5
    use temperature_integration_module, only: temperature_jac
    use jacobian_sparsity_module, only: get_jac_entry, set_jac_entry, set_jac_zero

    implicit none
    
    type(burn_t), intent(in) :: state
    real(rt), intent(inout) :: jac(njrows, njcols)

    type(rate_eval_t) :: rate_eval
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
    if (.not. disable_thermal_neutrinos) then
       call sneut5(state % T, state % rho, state % abar, state % zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)

       do j = 1, nspec
          b1 = (-state % abar * state % abar * snuda + (zion(j) - state % zbar) * state % abar * snudz)
          call get_jac_entry(jac, net_ienuc, j, scratch)
          scratch = scratch - b1
          call set_jac_entry(jac, net_ienuc, j, scratch)
       enddo
    endif

    ! Energy generation rate Jacobian element with respect to temperature
    do k = 1, nspec
       call get_jac_entry(jac, k, net_itemp, yderivs(k))
    enddo
    call ener_gener_rate(yderivs, scratch)
    if (.not. disable_thermal_neutrinos) then
       scratch = scratch - dsneutdt
    endif
    call set_jac_entry(jac, net_ienuc, net_itemp, scratch)

    ! Temperature Jacobian elements
    call temperature_jac(state, jac)

  end subroutine actual_jac


  subroutine jac_nuc(state, jac, Y, screened_rates)

    !$acc routine seq

    use jacobian_sparsity_module, only: set_jac_entry

    implicit none

    type(burn_t), intent(in) :: state
    real(rt), intent(inout) :: jac(njrows, njcols)

    real(rt), intent(in)  :: Y(nspec)
    real(rt), intent(in)  :: screened_rates(nrates)
    real(rt) :: scratch


    !$gpu


    scratch = (&
      -screened_rates(k_p_c12__n13)*Y(jc12)*state % rho - screened_rates(k_p_c13__n14)*Y(jc13)* &
      state % rho - screened_rates(k_p_f18__he4_o15)*Y(jf18)*state % rho - &
      screened_rates(k_p_n13__o14)*Y(jn13)*state % rho - screened_rates(k_p_n14__o15)* &
      Y(jn14)*state % rho - screened_rates(k_p_n15__he4_c12)*Y(jn15)*state % rho - &
      screened_rates(k_p_n15__o16)*Y(jn15)*state % rho - screened_rates(k_p_o16__f17)* &
      Y(jo16)*state % rho - screened_rates(k_p_o17__f18)*Y(jo17)*state % rho - &
      screened_rates(k_p_o17__he4_n14)*Y(jo17)*state % rho &
       )
    call set_jac_entry(jac, jp, jp, scratch)

    scratch = (&
      screened_rates(k_he4_n13__p_o16)*Y(jn13)*state % rho + screened_rates(k_he4_o14__p_f17)* &
      Y(jo14)*state % rho &
       )
    call set_jac_entry(jac, jp, jhe4, scratch)

    scratch = (&
      -screened_rates(k_p_c12__n13)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jp, jc12, scratch)

    scratch = (&
      -screened_rates(k_p_c13__n14)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jp, jc13, scratch)

    scratch = (&
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*state % rho - screened_rates(k_p_n13__o14)*Y(jp) &
      *state % rho &
       )
    call set_jac_entry(jac, jp, jn13, scratch)

    scratch = (&
      -screened_rates(k_p_n14__o15)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jp, jn14, scratch)

    scratch = (&
      -screened_rates(k_p_n15__he4_c12)*Y(jp)*state % rho - screened_rates(k_p_n15__o16)*Y(jp)* &
      state % rho &
       )
    call set_jac_entry(jac, jp, jn15, scratch)

    scratch = (&
      screened_rates(k_he4_o14__p_f17)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jp, jo14, scratch)

    scratch = (&
      -screened_rates(k_p_o16__f17)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jp, jo16, scratch)

    scratch = (&
      -screened_rates(k_p_o17__f18)*Y(jp)*state % rho - screened_rates(k_p_o17__he4_n14)*Y(jp)* &
      state % rho &
       )
    call set_jac_entry(jac, jp, jo17, scratch)

    scratch = (&
      -screened_rates(k_p_f18__he4_o15)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jp, jf18, scratch)

    scratch = (&
      screened_rates(k_p_f18__he4_o15)*Y(jf18)*state % rho + screened_rates(k_p_n15__he4_c12)* &
      Y(jn15)*state % rho + screened_rates(k_p_o17__he4_n14)*Y(jo17)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho - 1.5e0_rt* &
      screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**2*state % rho**2 - &
      screened_rates(k_he4_n13__p_o16)*Y(jn13)*state % rho - &
      screened_rates(k_he4_n14__f18)*Y(jn14)*state % rho - &
      screened_rates(k_he4_o14__p_f17)*Y(jo14)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jhe4, scratch)

    scratch = (&
      -screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jc12, scratch)

    scratch = (&
      -screened_rates(k_he4_n13__p_o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jn13, scratch)

    scratch = (&
      -screened_rates(k_he4_n14__f18)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jn14, scratch)

    scratch = (&
      screened_rates(k_p_n15__he4_c12)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jn15, scratch)

    scratch = (&
      -screened_rates(k_he4_o14__p_f17)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jo14, scratch)

    scratch = (&
      screened_rates(k_p_o17__he4_n14)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jo17, scratch)

    scratch = (&
      screened_rates(k_p_f18__he4_o15)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jf18, scratch)

    scratch = (&
      -screened_rates(k_p_c12__n13)*Y(jc12)*state % rho + screened_rates(k_p_n15__he4_c12)* &
      Y(jn15)*state % rho &
       )
    call set_jac_entry(jac, jc12, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho + 0.5e0_rt* &
      screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**2*state % rho**2 &
       )
    call set_jac_entry(jac, jc12, jhe4, scratch)

    scratch = (&
      -screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho - screened_rates(k_p_c12__n13)*Y(jp)* &
      state % rho &
       )
    call set_jac_entry(jac, jc12, jc12, scratch)

    scratch = (&
      screened_rates(k_p_n15__he4_c12)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jc12, jn15, scratch)

    scratch = (&
      -screened_rates(k_p_c13__n14)*Y(jc13)*state % rho &
       )
    call set_jac_entry(jac, jc13, jp, scratch)

    scratch = (&
      -screened_rates(k_p_c13__n14)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jc13, jc13, scratch)

    scratch = (&
      screened_rates(k_n13__c13__weak__wc12) &
       )
    call set_jac_entry(jac, jc13, jn13, scratch)

    scratch = (&
      screened_rates(k_p_c12__n13)*Y(jc12)*state % rho - screened_rates(k_p_n13__o14)*Y(jn13)* &
      state % rho &
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
      -screened_rates(k_he4_n13__p_o16)*Y(jhe4)*state % rho - &
      screened_rates(k_n13__c13__weak__wc12) - screened_rates(k_p_n13__o14)*Y(jp)* &
      state % rho &
       )
    call set_jac_entry(jac, jn13, jn13, scratch)

    scratch = (&
      screened_rates(k_p_c13__n14)*Y(jc13)*state % rho - screened_rates(k_p_n14__o15)*Y(jn14)* &
      state % rho + screened_rates(k_p_o17__he4_n14)*Y(jo17)*state % rho &
       )
    call set_jac_entry(jac, jn14, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_n14__f18)*Y(jn14)*state % rho &
       )
    call set_jac_entry(jac, jn14, jhe4, scratch)

    scratch = (&
      screened_rates(k_p_c13__n14)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jn14, jc13, scratch)

    scratch = (&
      -screened_rates(k_he4_n14__f18)*Y(jhe4)*state % rho - screened_rates(k_p_n14__o15)*Y(jp)* &
      state % rho &
       )
    call set_jac_entry(jac, jn14, jn14, scratch)

    scratch = (&
      screened_rates(k_o14__n14__weak__wc12) &
       )
    call set_jac_entry(jac, jn14, jo14, scratch)

    scratch = (&
      screened_rates(k_p_o17__he4_n14)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jn14, jo17, scratch)

    scratch = (&
      -screened_rates(k_p_n15__he4_c12)*Y(jn15)*state % rho - screened_rates(k_p_n15__o16)* &
      Y(jn15)*state % rho &
       )
    call set_jac_entry(jac, jn15, jp, scratch)

    scratch = (&
      -screened_rates(k_p_n15__he4_c12)*Y(jp)*state % rho - screened_rates(k_p_n15__o16)*Y(jp)* &
      state % rho &
       )
    call set_jac_entry(jac, jn15, jn15, scratch)

    scratch = (&
      screened_rates(k_o15__n15__weak__wc12) &
       )
    call set_jac_entry(jac, jn15, jo15, scratch)

    scratch = (&
      screened_rates(k_p_n13__o14)*Y(jn13)*state % rho &
       )
    call set_jac_entry(jac, jo14, jp, scratch)

    scratch = (&
      -screened_rates(k_he4_o14__p_f17)*Y(jo14)*state % rho &
       )
    call set_jac_entry(jac, jo14, jhe4, scratch)

    scratch = (&
      screened_rates(k_p_n13__o14)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jo14, jn13, scratch)

    scratch = (&
      -screened_rates(k_he4_o14__p_f17)*Y(jhe4)*state % rho - &
      screened_rates(k_o14__n14__weak__wc12) &
       )
    call set_jac_entry(jac, jo14, jo14, scratch)

    scratch = (&
      screened_rates(k_p_f18__he4_o15)*Y(jf18)*state % rho + screened_rates(k_p_n14__o15)* &
      Y(jn14)*state % rho &
       )
    call set_jac_entry(jac, jo15, jp, scratch)

    scratch = (&
      screened_rates(k_p_n14__o15)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jo15, jn14, scratch)

    scratch = (&
      -screened_rates(k_o15__n15__weak__wc12) &
       )
    call set_jac_entry(jac, jo15, jo15, scratch)

    scratch = (&
      screened_rates(k_p_f18__he4_o15)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jo15, jf18, scratch)

    scratch = (&
      screened_rates(k_p_n15__o16)*Y(jn15)*state % rho - screened_rates(k_p_o16__f17)*Y(jo16)* &
      state % rho &
       )
    call set_jac_entry(jac, jo16, jp, scratch)

    scratch = (&
      screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho + screened_rates(k_he4_n13__p_o16)* &
      Y(jn13)*state % rho &
       )
    call set_jac_entry(jac, jo16, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jo16, jc12, scratch)

    scratch = (&
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jo16, jn13, scratch)

    scratch = (&
      screened_rates(k_p_n15__o16)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jo16, jn15, scratch)

    scratch = (&
      -screened_rates(k_p_o16__f17)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jo16, jo16, scratch)

    scratch = (&
      -screened_rates(k_p_o17__f18)*Y(jo17)*state % rho - screened_rates(k_p_o17__he4_n14)* &
      Y(jo17)*state % rho &
       )
    call set_jac_entry(jac, jo17, jp, scratch)

    scratch = (&
      -screened_rates(k_p_o17__f18)*Y(jp)*state % rho - screened_rates(k_p_o17__he4_n14)*Y(jp)* &
      state % rho &
       )
    call set_jac_entry(jac, jo17, jo17, scratch)

    scratch = (&
      screened_rates(k_f17__o17__weak__wc12) &
       )
    call set_jac_entry(jac, jo17, jf17, scratch)

    scratch = (&
      screened_rates(k_p_o16__f17)*Y(jo16)*state % rho &
       )
    call set_jac_entry(jac, jf17, jp, scratch)

    scratch = (&
      screened_rates(k_he4_o14__p_f17)*Y(jo14)*state % rho &
       )
    call set_jac_entry(jac, jf17, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_o14__p_f17)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jf17, jo14, scratch)

    scratch = (&
      screened_rates(k_p_o16__f17)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jf17, jo16, scratch)

    scratch = (&
      -screened_rates(k_f17__o17__weak__wc12) &
       )
    call set_jac_entry(jac, jf17, jf17, scratch)

    scratch = (&
      -screened_rates(k_p_f18__he4_o15)*Y(jf18)*state % rho + screened_rates(k_p_o17__f18)* &
      Y(jo17)*state % rho &
       )
    call set_jac_entry(jac, jf18, jp, scratch)

    scratch = (&
      screened_rates(k_he4_n14__f18)*Y(jn14)*state % rho &
       )
    call set_jac_entry(jac, jf18, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_n14__f18)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jf18, jn14, scratch)

    scratch = (&
      screened_rates(k_p_o17__f18)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jf18, jo17, scratch)

    scratch = (&
      -screened_rates(k_p_f18__he4_o15)*Y(jp)*state % rho &
       )
    call set_jac_entry(jac, jf18, jf18, scratch)


  end subroutine jac_nuc

end module actual_rhs_module
