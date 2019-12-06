module actual_rhs_module

  use microphysics_type_module
  use amrex_constants_module
  use physical_constants, only: N_AVO
  use network
  use reaclib_rates, only: screen_reaclib, reaclib_evaluate
  use table_rates
  use screening_module, only: screen5, plasma_state, fill_plasma_state
  use sneut_module, only: sneut5
  use temperature_integration_module, only: temperature_rhs, temperature_jac
  use burn_type_module

  implicit none

  type :: rate_eval_t
     real(rt) :: unscreened_rates(4, nrates)
     real(rt) :: screened_rates(nrates)
     real(rt) :: dqweak(nrat_tabular)
     real(rt) :: epart(nrat_tabular)
  end type rate_eval_t
  
contains

  subroutine actual_rhs_init()
    ! STUB FOR MAESTRO'S TEST_REACT. ALL THE INIT IS DONE BY BURNER_INIT
    return
  end subroutine actual_rhs_init
  
  subroutine update_unevolved_species(state)
    ! STUB FOR INTEGRATOR
    type(burn_t)     :: state
    return
  end subroutine update_unevolved_species

  subroutine ener_gener_rate(dydt, enuc)
    ! Computes the instantaneous energy generation rate
    !$acc routine seq

    implicit none

    real(rt) :: dydt(nspec), enuc

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * mion(:)) * enuc_conv2

  end subroutine ener_gener_rate

  subroutine evaluate_rates(state, rate_eval)
    !$acc routine seq

    implicit none
    
    type(burn_t)     :: state
    type(rate_eval_t), intent(out) :: rate_eval
    type(plasma_state) :: pstate
    real(rt) :: Y(nspec)
    real(rt) :: raw_rates(4, nrates)
    real(rt) :: reactvec(num_rate_groups+2)
    integer :: i, j
    real(rt) :: dens, temp, rhoy, scor, dscor_dt, dscor_dd

    Y(:) = state%xn(:) * aion_inv(:)
    dens = state%rho
    temp = state%T
    rhoy = dens * state%y_e

    ! Calculate Reaclib rates
    call fill_plasma_state(pstate, temp, dens, Y)
    do i = 1, nrat_reaclib
       call reaclib_evaluate(pstate, temp, i, reactvec)
       rate_eval % unscreened_rates(:,i) = reactvec(1:4)
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

  subroutine actual_rhs(state)
    
    !$acc routine seq

    use extern_probin_module, only: do_constant_volume_burn
    use burn_type_module, only: net_itemp, net_ienuc

    implicit none

    type(burn_t) :: state
    type(rate_eval_t) :: rate_eval
    type(plasma_state) :: pstate
    real(rt) :: Y(nspec)
    real(rt) :: reactvec(num_rate_groups+2)
    integer :: i, j
    real(rt) :: dens, temp, rhoy, ye, enuc
    real(rt) :: sneut, dsneutdt, dsneutdd, snuda, snudz

    ! Set molar abundances
    Y(:) = state%xn(:) * aion_inv(:)

    dens = state%rho
    temp = state%T

    call evaluate_rates(state, rate_eval)

    call rhs_nuc(state, state % ydot(1:nspec), Y, rate_eval % screened_rates)

    ! ion binding energy contributions
    call ener_gener_rate(state%ydot(1:nspec), enuc)
    
    ! weak Q-value modification dqweak (density and temperature dependent)
    
    ! weak particle energy generation rates from gamma heating and neutrino loss
    ! (does not include plasma neutrino losses)


    ! Get the neutrino losses
    call sneut5(temp, dens, state%abar, state%zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)

    ! Append the energy equation (this is erg/g/s)
    state%ydot(net_ienuc) = enuc - sneut

    ! Append the temperature equation
    call temperature_rhs(state)
    
    ! write(*,*) '______________________________'
    ! do i = 1, nspec+2
    !    write(*,*) 'state%ydot(',i,'): ',state%ydot(i)
    ! end do
  end subroutine actual_rhs

  subroutine rhs_nuc(state, ydot_nuc, Y, screened_rates)

    !$acc routine seq

    type(burn_t),   intent(inout) :: state
    real(rt), intent(out) :: ydot_nuc(nspec)
    real(rt), intent(in)  :: Y(nspec)
    real(rt), intent(in)  :: screened_rates(nrates)



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

  
  subroutine actual_jac(state)

    !$acc routine seq

    use burn_type_module, only: net_itemp, net_ienuc
    
    implicit none
    
    type(burn_t) :: state
    type(rate_eval_t) :: rate_eval
    type(plasma_state) :: pstate
    real(rt) :: reactvec(num_rate_groups+2)
    real(rt) :: screened_rates_dt(nrates)
    real(rt) :: Y(nspec)
    real(rt) :: dens, temp, ye, rhoy, b1
    real(rt) :: sneut, dsneutdt, dsneutdd, snuda, snudz
    integer :: i, j

    dens = state%rho
    temp = state%T

    ! Set molar abundances
    Y(:) = state%xn(:) * aion_inv(:)
    
    call evaluate_rates(state, rate_eval)
    
    ! Species Jacobian elements with respect to other species
    call jac_nuc(state, Y, rate_eval % screened_rates)

    ! Species Jacobian elements with respect to energy generation rate
    state%jac(1:nspec, net_ienuc) = 0.0e0_rt

    ! Evaluate the species Jacobian elements with respect to temperature by
    ! calling the RHS using the temperature derivative of the screened rate
    screened_rates_dt = rate_eval % unscreened_rates(i_rate, :) * &
         rate_eval % unscreened_rates(i_dscor_dt, :) + &
         rate_eval % unscreened_rates(i_drate_dt, :) * &
         rate_eval % unscreened_rates(i_scor, :)

    call rhs_nuc(state, state%jac(1:nspec, net_itemp), Y, screened_rates_dt)
    
    ! Energy generation rate Jacobian elements with respect to species
    do j = 1, nspec
       call ener_gener_rate(state % jac(1:nspec,j), state % jac(net_ienuc,j))
    enddo

    ! Account for the thermal neutrino losses
    call sneut5(temp, dens, state%abar, state%zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)
    do j = 1, nspec
       b1 = (-state%abar * state%abar * snuda + (zion(j) - state%zbar) * state%abar * snudz)
       state % jac(net_ienuc,j) = state % jac(net_ienuc,j) - b1
    enddo

    ! Energy generation rate Jacobian element with respect to energy generation rate
    state%jac(net_ienuc, net_ienuc) = 0.0e0_rt

    ! Energy generation rate Jacobian element with respect to temperature
    call ener_gener_rate(state%jac(1:nspec, net_itemp), state%jac(net_ienuc, net_itemp))
    state%jac(net_ienuc, net_itemp) = state%jac(net_ienuc, net_itemp) - dsneutdt

    ! Add dqweak and epart contributions!!!

    ! Temperature Jacobian elements
    call temperature_jac(state)

  end subroutine actual_jac

  subroutine jac_nuc(state, Y, screened_rates)

    !$acc routine seq
    
    type(burn_t),   intent(inout) :: state
    real(rt), intent(in)  :: Y(nspec)
    real(rt), intent(in)  :: screened_rates(nrates)



    state % jac(jp,jp) = ( &
      -screened_rates(k_p_c12__n13)*Y(jc12)*state % rho - screened_rates(k_p_c13__n14)*Y(jc13)* &
      state % rho - screened_rates(k_p_f18__he4_o15)*Y(jf18)*state % rho - &
      screened_rates(k_p_n13__o14)*Y(jn13)*state % rho - screened_rates(k_p_n14__o15)* &
      Y(jn14)*state % rho - screened_rates(k_p_n15__he4_c12)*Y(jn15)*state % rho - &
      screened_rates(k_p_n15__o16)*Y(jn15)*state % rho - screened_rates(k_p_o16__f17)* &
      Y(jo16)*state % rho - screened_rates(k_p_o17__f18)*Y(jo17)*state % rho - &
      screened_rates(k_p_o17__he4_n14)*Y(jo17)*state % rho &
       )

    state % jac(jp,jhe4) = ( &
      screened_rates(k_he4_n13__p_o16)*Y(jn13)*state % rho + screened_rates(k_he4_o14__p_f17)* &
      Y(jo14)*state % rho &
       )

    state % jac(jp,jc12) = ( &
      -screened_rates(k_p_c12__n13)*Y(jp)*state % rho &
       )

    state % jac(jp,jc13) = ( &
      -screened_rates(k_p_c13__n14)*Y(jp)*state % rho &
       )

    state % jac(jp,jn13) = ( &
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*state % rho - screened_rates(k_p_n13__o14)*Y(jp) &
      *state % rho &
       )

    state % jac(jp,jn14) = ( &
      -screened_rates(k_p_n14__o15)*Y(jp)*state % rho &
       )

    state % jac(jp,jn15) = ( &
      -screened_rates(k_p_n15__he4_c12)*Y(jp)*state % rho - screened_rates(k_p_n15__o16)*Y(jp)* &
      state % rho &
       )

    state % jac(jp,jo14) = ( &
      screened_rates(k_he4_o14__p_f17)*Y(jhe4)*state % rho &
       )

    state % jac(jp,jo15) = ( &
      0.0e0_rt &
       )

    state % jac(jp,jo16) = ( &
      -screened_rates(k_p_o16__f17)*Y(jp)*state % rho &
       )

    state % jac(jp,jo17) = ( &
      -screened_rates(k_p_o17__f18)*Y(jp)*state % rho - screened_rates(k_p_o17__he4_n14)*Y(jp)* &
      state % rho &
       )

    state % jac(jp,jf17) = ( &
      0.0e0_rt &
       )

    state % jac(jp,jf18) = ( &
      -screened_rates(k_p_f18__he4_o15)*Y(jp)*state % rho &
       )

    state % jac(jhe4,jp) = ( &
      screened_rates(k_p_f18__he4_o15)*Y(jf18)*state % rho + screened_rates(k_p_n15__he4_c12)* &
      Y(jn15)*state % rho + screened_rates(k_p_o17__he4_n14)*Y(jo17)*state % rho &
       )

    state % jac(jhe4,jhe4) = ( &
      -screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho - 1.5e0_rt* &
      screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**2*state % rho**2 - &
      screened_rates(k_he4_n13__p_o16)*Y(jn13)*state % rho - &
      screened_rates(k_he4_n14__f18)*Y(jn14)*state % rho - &
      screened_rates(k_he4_o14__p_f17)*Y(jo14)*state % rho &
       )

    state % jac(jhe4,jc12) = ( &
      -screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho &
       )

    state % jac(jhe4,jc13) = ( &
      0.0e0_rt &
       )

    state % jac(jhe4,jn13) = ( &
      -screened_rates(k_he4_n13__p_o16)*Y(jhe4)*state % rho &
       )

    state % jac(jhe4,jn14) = ( &
      -screened_rates(k_he4_n14__f18)*Y(jhe4)*state % rho &
       )

    state % jac(jhe4,jn15) = ( &
      screened_rates(k_p_n15__he4_c12)*Y(jp)*state % rho &
       )

    state % jac(jhe4,jo14) = ( &
      -screened_rates(k_he4_o14__p_f17)*Y(jhe4)*state % rho &
       )

    state % jac(jhe4,jo15) = ( &
      0.0e0_rt &
       )

    state % jac(jhe4,jo16) = ( &
      0.0e0_rt &
       )

    state % jac(jhe4,jo17) = ( &
      screened_rates(k_p_o17__he4_n14)*Y(jp)*state % rho &
       )

    state % jac(jhe4,jf17) = ( &
      0.0e0_rt &
       )

    state % jac(jhe4,jf18) = ( &
      screened_rates(k_p_f18__he4_o15)*Y(jp)*state % rho &
       )

    state % jac(jc12,jp) = ( &
      -screened_rates(k_p_c12__n13)*Y(jc12)*state % rho + screened_rates(k_p_n15__he4_c12)* &
      Y(jn15)*state % rho &
       )

    state % jac(jc12,jhe4) = ( &
      -screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho + 0.5e0_rt* &
      screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**2*state % rho**2 &
       )

    state % jac(jc12,jc12) = ( &
      -screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho - screened_rates(k_p_c12__n13)*Y(jp)* &
      state % rho &
       )

    state % jac(jc12,jc13) = ( &
      0.0e0_rt &
       )

    state % jac(jc12,jn13) = ( &
      0.0e0_rt &
       )

    state % jac(jc12,jn14) = ( &
      0.0e0_rt &
       )

    state % jac(jc12,jn15) = ( &
      screened_rates(k_p_n15__he4_c12)*Y(jp)*state % rho &
       )

    state % jac(jc12,jo14) = ( &
      0.0e0_rt &
       )

    state % jac(jc12,jo15) = ( &
      0.0e0_rt &
       )

    state % jac(jc12,jo16) = ( &
      0.0e0_rt &
       )

    state % jac(jc12,jo17) = ( &
      0.0e0_rt &
       )

    state % jac(jc12,jf17) = ( &
      0.0e0_rt &
       )

    state % jac(jc12,jf18) = ( &
      0.0e0_rt &
       )

    state % jac(jc13,jp) = ( &
      -screened_rates(k_p_c13__n14)*Y(jc13)*state % rho &
       )

    state % jac(jc13,jhe4) = ( &
      0.0e0_rt &
       )

    state % jac(jc13,jc12) = ( &
      0.0e0_rt &
       )

    state % jac(jc13,jc13) = ( &
      -screened_rates(k_p_c13__n14)*Y(jp)*state % rho &
       )

    state % jac(jc13,jn13) = ( &
      screened_rates(k_n13__c13__weak__wc12) &
       )

    state % jac(jc13,jn14) = ( &
      0.0e0_rt &
       )

    state % jac(jc13,jn15) = ( &
      0.0e0_rt &
       )

    state % jac(jc13,jo14) = ( &
      0.0e0_rt &
       )

    state % jac(jc13,jo15) = ( &
      0.0e0_rt &
       )

    state % jac(jc13,jo16) = ( &
      0.0e0_rt &
       )

    state % jac(jc13,jo17) = ( &
      0.0e0_rt &
       )

    state % jac(jc13,jf17) = ( &
      0.0e0_rt &
       )

    state % jac(jc13,jf18) = ( &
      0.0e0_rt &
       )

    state % jac(jn13,jp) = ( &
      screened_rates(k_p_c12__n13)*Y(jc12)*state % rho - screened_rates(k_p_n13__o14)*Y(jn13)* &
      state % rho &
       )

    state % jac(jn13,jhe4) = ( &
      -screened_rates(k_he4_n13__p_o16)*Y(jn13)*state % rho &
       )

    state % jac(jn13,jc12) = ( &
      screened_rates(k_p_c12__n13)*Y(jp)*state % rho &
       )

    state % jac(jn13,jc13) = ( &
      0.0e0_rt &
       )

    state % jac(jn13,jn13) = ( &
      -screened_rates(k_he4_n13__p_o16)*Y(jhe4)*state % rho - &
      screened_rates(k_n13__c13__weak__wc12) - screened_rates(k_p_n13__o14)*Y(jp)* &
      state % rho &
       )

    state % jac(jn13,jn14) = ( &
      0.0e0_rt &
       )

    state % jac(jn13,jn15) = ( &
      0.0e0_rt &
       )

    state % jac(jn13,jo14) = ( &
      0.0e0_rt &
       )

    state % jac(jn13,jo15) = ( &
      0.0e0_rt &
       )

    state % jac(jn13,jo16) = ( &
      0.0e0_rt &
       )

    state % jac(jn13,jo17) = ( &
      0.0e0_rt &
       )

    state % jac(jn13,jf17) = ( &
      0.0e0_rt &
       )

    state % jac(jn13,jf18) = ( &
      0.0e0_rt &
       )

    state % jac(jn14,jp) = ( &
      screened_rates(k_p_c13__n14)*Y(jc13)*state % rho - screened_rates(k_p_n14__o15)*Y(jn14)* &
      state % rho + screened_rates(k_p_o17__he4_n14)*Y(jo17)*state % rho &
       )

    state % jac(jn14,jhe4) = ( &
      -screened_rates(k_he4_n14__f18)*Y(jn14)*state % rho &
       )

    state % jac(jn14,jc12) = ( &
      0.0e0_rt &
       )

    state % jac(jn14,jc13) = ( &
      screened_rates(k_p_c13__n14)*Y(jp)*state % rho &
       )

    state % jac(jn14,jn13) = ( &
      0.0e0_rt &
       )

    state % jac(jn14,jn14) = ( &
      -screened_rates(k_he4_n14__f18)*Y(jhe4)*state % rho - screened_rates(k_p_n14__o15)*Y(jp)* &
      state % rho &
       )

    state % jac(jn14,jn15) = ( &
      0.0e0_rt &
       )

    state % jac(jn14,jo14) = ( &
      screened_rates(k_o14__n14__weak__wc12) &
       )

    state % jac(jn14,jo15) = ( &
      0.0e0_rt &
       )

    state % jac(jn14,jo16) = ( &
      0.0e0_rt &
       )

    state % jac(jn14,jo17) = ( &
      screened_rates(k_p_o17__he4_n14)*Y(jp)*state % rho &
       )

    state % jac(jn14,jf17) = ( &
      0.0e0_rt &
       )

    state % jac(jn14,jf18) = ( &
      0.0e0_rt &
       )

    state % jac(jn15,jp) = ( &
      -screened_rates(k_p_n15__he4_c12)*Y(jn15)*state % rho - screened_rates(k_p_n15__o16)* &
      Y(jn15)*state % rho &
       )

    state % jac(jn15,jhe4) = ( &
      0.0e0_rt &
       )

    state % jac(jn15,jc12) = ( &
      0.0e0_rt &
       )

    state % jac(jn15,jc13) = ( &
      0.0e0_rt &
       )

    state % jac(jn15,jn13) = ( &
      0.0e0_rt &
       )

    state % jac(jn15,jn14) = ( &
      0.0e0_rt &
       )

    state % jac(jn15,jn15) = ( &
      -screened_rates(k_p_n15__he4_c12)*Y(jp)*state % rho - screened_rates(k_p_n15__o16)*Y(jp)* &
      state % rho &
       )

    state % jac(jn15,jo14) = ( &
      0.0e0_rt &
       )

    state % jac(jn15,jo15) = ( &
      screened_rates(k_o15__n15__weak__wc12) &
       )

    state % jac(jn15,jo16) = ( &
      0.0e0_rt &
       )

    state % jac(jn15,jo17) = ( &
      0.0e0_rt &
       )

    state % jac(jn15,jf17) = ( &
      0.0e0_rt &
       )

    state % jac(jn15,jf18) = ( &
      0.0e0_rt &
       )

    state % jac(jo14,jp) = ( &
      screened_rates(k_p_n13__o14)*Y(jn13)*state % rho &
       )

    state % jac(jo14,jhe4) = ( &
      -screened_rates(k_he4_o14__p_f17)*Y(jo14)*state % rho &
       )

    state % jac(jo14,jc12) = ( &
      0.0e0_rt &
       )

    state % jac(jo14,jc13) = ( &
      0.0e0_rt &
       )

    state % jac(jo14,jn13) = ( &
      screened_rates(k_p_n13__o14)*Y(jp)*state % rho &
       )

    state % jac(jo14,jn14) = ( &
      0.0e0_rt &
       )

    state % jac(jo14,jn15) = ( &
      0.0e0_rt &
       )

    state % jac(jo14,jo14) = ( &
      -screened_rates(k_he4_o14__p_f17)*Y(jhe4)*state % rho - &
      screened_rates(k_o14__n14__weak__wc12) &
       )

    state % jac(jo14,jo15) = ( &
      0.0e0_rt &
       )

    state % jac(jo14,jo16) = ( &
      0.0e0_rt &
       )

    state % jac(jo14,jo17) = ( &
      0.0e0_rt &
       )

    state % jac(jo14,jf17) = ( &
      0.0e0_rt &
       )

    state % jac(jo14,jf18) = ( &
      0.0e0_rt &
       )

    state % jac(jo15,jp) = ( &
      screened_rates(k_p_f18__he4_o15)*Y(jf18)*state % rho + screened_rates(k_p_n14__o15)* &
      Y(jn14)*state % rho &
       )

    state % jac(jo15,jhe4) = ( &
      0.0e0_rt &
       )

    state % jac(jo15,jc12) = ( &
      0.0e0_rt &
       )

    state % jac(jo15,jc13) = ( &
      0.0e0_rt &
       )

    state % jac(jo15,jn13) = ( &
      0.0e0_rt &
       )

    state % jac(jo15,jn14) = ( &
      screened_rates(k_p_n14__o15)*Y(jp)*state % rho &
       )

    state % jac(jo15,jn15) = ( &
      0.0e0_rt &
       )

    state % jac(jo15,jo14) = ( &
      0.0e0_rt &
       )

    state % jac(jo15,jo15) = ( &
      -screened_rates(k_o15__n15__weak__wc12) &
       )

    state % jac(jo15,jo16) = ( &
      0.0e0_rt &
       )

    state % jac(jo15,jo17) = ( &
      0.0e0_rt &
       )

    state % jac(jo15,jf17) = ( &
      0.0e0_rt &
       )

    state % jac(jo15,jf18) = ( &
      screened_rates(k_p_f18__he4_o15)*Y(jp)*state % rho &
       )

    state % jac(jo16,jp) = ( &
      screened_rates(k_p_n15__o16)*Y(jn15)*state % rho - screened_rates(k_p_o16__f17)*Y(jo16)* &
      state % rho &
       )

    state % jac(jo16,jhe4) = ( &
      screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho + screened_rates(k_he4_n13__p_o16)* &
      Y(jn13)*state % rho &
       )

    state % jac(jo16,jc12) = ( &
      screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho &
       )

    state % jac(jo16,jc13) = ( &
      0.0e0_rt &
       )

    state % jac(jo16,jn13) = ( &
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*state % rho &
       )

    state % jac(jo16,jn14) = ( &
      0.0e0_rt &
       )

    state % jac(jo16,jn15) = ( &
      screened_rates(k_p_n15__o16)*Y(jp)*state % rho &
       )

    state % jac(jo16,jo14) = ( &
      0.0e0_rt &
       )

    state % jac(jo16,jo15) = ( &
      0.0e0_rt &
       )

    state % jac(jo16,jo16) = ( &
      -screened_rates(k_p_o16__f17)*Y(jp)*state % rho &
       )

    state % jac(jo16,jo17) = ( &
      0.0e0_rt &
       )

    state % jac(jo16,jf17) = ( &
      0.0e0_rt &
       )

    state % jac(jo16,jf18) = ( &
      0.0e0_rt &
       )

    state % jac(jo17,jp) = ( &
      -screened_rates(k_p_o17__f18)*Y(jo17)*state % rho - screened_rates(k_p_o17__he4_n14)* &
      Y(jo17)*state % rho &
       )

    state % jac(jo17,jhe4) = ( &
      0.0e0_rt &
       )

    state % jac(jo17,jc12) = ( &
      0.0e0_rt &
       )

    state % jac(jo17,jc13) = ( &
      0.0e0_rt &
       )

    state % jac(jo17,jn13) = ( &
      0.0e0_rt &
       )

    state % jac(jo17,jn14) = ( &
      0.0e0_rt &
       )

    state % jac(jo17,jn15) = ( &
      0.0e0_rt &
       )

    state % jac(jo17,jo14) = ( &
      0.0e0_rt &
       )

    state % jac(jo17,jo15) = ( &
      0.0e0_rt &
       )

    state % jac(jo17,jo16) = ( &
      0.0e0_rt &
       )

    state % jac(jo17,jo17) = ( &
      -screened_rates(k_p_o17__f18)*Y(jp)*state % rho - screened_rates(k_p_o17__he4_n14)*Y(jp)* &
      state % rho &
       )

    state % jac(jo17,jf17) = ( &
      screened_rates(k_f17__o17__weak__wc12) &
       )

    state % jac(jo17,jf18) = ( &
      0.0e0_rt &
       )

    state % jac(jf17,jp) = ( &
      screened_rates(k_p_o16__f17)*Y(jo16)*state % rho &
       )

    state % jac(jf17,jhe4) = ( &
      screened_rates(k_he4_o14__p_f17)*Y(jo14)*state % rho &
       )

    state % jac(jf17,jc12) = ( &
      0.0e0_rt &
       )

    state % jac(jf17,jc13) = ( &
      0.0e0_rt &
       )

    state % jac(jf17,jn13) = ( &
      0.0e0_rt &
       )

    state % jac(jf17,jn14) = ( &
      0.0e0_rt &
       )

    state % jac(jf17,jn15) = ( &
      0.0e0_rt &
       )

    state % jac(jf17,jo14) = ( &
      screened_rates(k_he4_o14__p_f17)*Y(jhe4)*state % rho &
       )

    state % jac(jf17,jo15) = ( &
      0.0e0_rt &
       )

    state % jac(jf17,jo16) = ( &
      screened_rates(k_p_o16__f17)*Y(jp)*state % rho &
       )

    state % jac(jf17,jo17) = ( &
      0.0e0_rt &
       )

    state % jac(jf17,jf17) = ( &
      -screened_rates(k_f17__o17__weak__wc12) &
       )

    state % jac(jf17,jf18) = ( &
      0.0e0_rt &
       )

    state % jac(jf18,jp) = ( &
      -screened_rates(k_p_f18__he4_o15)*Y(jf18)*state % rho + screened_rates(k_p_o17__f18)* &
      Y(jo17)*state % rho &
       )

    state % jac(jf18,jhe4) = ( &
      screened_rates(k_he4_n14__f18)*Y(jn14)*state % rho &
       )

    state % jac(jf18,jc12) = ( &
      0.0e0_rt &
       )

    state % jac(jf18,jc13) = ( &
      0.0e0_rt &
       )

    state % jac(jf18,jn13) = ( &
      0.0e0_rt &
       )

    state % jac(jf18,jn14) = ( &
      screened_rates(k_he4_n14__f18)*Y(jhe4)*state % rho &
       )

    state % jac(jf18,jn15) = ( &
      0.0e0_rt &
       )

    state % jac(jf18,jo14) = ( &
      0.0e0_rt &
       )

    state % jac(jf18,jo15) = ( &
      0.0e0_rt &
       )

    state % jac(jf18,jo16) = ( &
      0.0e0_rt &
       )

    state % jac(jf18,jo17) = ( &
      screened_rates(k_p_o17__f18)*Y(jp)*state % rho &
       )

    state % jac(jf18,jf17) = ( &
      0.0e0_rt &
       )

    state % jac(jf18,jf18) = ( &
      -screened_rates(k_p_f18__he4_o15)*Y(jp)*state % rho &
       )

    
  end subroutine jac_nuc

end module actual_rhs_module
