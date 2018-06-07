module actual_rhs_module

  use bl_types
  use bl_constants_module
  use physical_constants, only: N_AVO
  use network
  use reaclib_rates, only: screen_reaclib, reaclib_evaluate
  use table_rates
  use screening_module, only: plasma_state, fill_plasma_state
  use sneut_module, only: sneut5
  use temperature_integration_module, only: temperature_rhs, temperature_jac
  use burn_type_module

  implicit none

  type :: rate_eval_t
     double precision :: unscreened_rates(4, nrates)
     double precision :: screened_rates(nrates)
     double precision :: dqweak(nrat_tabular)
     double precision :: epart(nrat_tabular)
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

  subroutine evaluate_rates(state, rate_eval)
    !$acc routine seq
    type(burn_t)     :: state
    type(rate_eval_t), intent(out) :: rate_eval
    type(plasma_state) :: pstate
    double precision :: Y(nspec)
    double precision :: raw_rates(4, nrates)
    double precision :: reactvec(num_rate_groups+2)
    integer :: i, j
    double precision :: dens, temp, rhoy

    Y(:) = state%xn(:) * aion_inv(:)
    dens = state%rho
    temp = state%T
    rhoy = dens*state%y_e
    
    ! Calculate Reaclib rates
    call fill_plasma_state(pstate, temp, dens, Y)
    do i = 1, nrat_reaclib
       call reaclib_evaluate(pstate, temp, i, reactvec)
       rate_eval % unscreened_rates(:,i) = reactvec(1:4)
    end do


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
    double precision :: Y(nspec)
    double precision :: ydot_nuc(nspec)
    double precision :: reactvec(num_rate_groups+2)
    integer :: i, j
    double precision :: dens, temp, rhoy, ye, enuc
    double precision :: sneut, dsneutdt, dsneutdd, snuda, snudz

    ! Set molar abundances
    Y(:) = state%xn(:) * aion_inv(:)

    dens = state%rho
    temp = state%T

    call evaluate_rates(state, rate_eval)

    call rhs_nuc(ydot_nuc, Y, rate_eval % screened_rates, dens)
    state%ydot(1:nspec) = ydot_nuc

    ! ion binding energy contributions
    call ener_gener_rate(ydot_nuc, enuc)
    
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

  subroutine rhs_nuc(ydot_nuc, Y, screened_rates, dens)

    !$acc routine seq

    
    double precision, intent(out) :: ydot_nuc(nspec)
    double precision, intent(in)  :: Y(nspec)
    double precision, intent(in)  :: screened_rates(nrates)
    double precision, intent(in)  :: dens



    ydot_nuc(jhe4) = ( &
      screened_rates(k_ar36__he4_s32)*Y(jar36) + 3.0d0*screened_rates(k_c12__he4_he4_he4)* &
      Y(jc12) + 0.5d0*screened_rates(k_c12_c12__he4_ne20)*Y(jc12)**2*dens + &
      screened_rates(k_c12_ne20__he4_si28)*Y(jc12)*Y(jne20)*dens + &
      screened_rates(k_c12_o16__he4_mg24)*Y(jc12)*Y(jo16)*dens + &
      screened_rates(k_ca40__he4_ar36)*Y(jca40) + screened_rates(k_cr48__he4_ti44)* &
      Y(jcr48) + screened_rates(k_fe52__he4_cr48)*Y(jfe52) - &
      screened_rates(k_he4_ar36__ca40)*Y(jar36)*Y(jhe4)*dens - &
      screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*dens - &
      screened_rates(k_he4_ca40__ti44)*Y(jca40)*Y(jhe4)*dens - &
      screened_rates(k_he4_cr48__fe52)*Y(jcr48)*Y(jhe4)*dens - &
      screened_rates(k_he4_fe52__ni56)*Y(jfe52)*Y(jhe4)*dens - 0.5d0* &
      screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**3*dens**2 - &
      screened_rates(k_he4_mg24__c12_o16)*Y(jhe4)*Y(jmg24)*dens - &
      screened_rates(k_he4_mg24__si28)*Y(jhe4)*Y(jmg24)*dens - &
      screened_rates(k_he4_ne20__c12_c12)*Y(jhe4)*Y(jne20)*dens - &
      screened_rates(k_he4_ne20__mg24)*Y(jhe4)*Y(jne20)*dens - &
      screened_rates(k_he4_ni56__zn60)*Y(jhe4)*Y(jni56)*dens - &
      screened_rates(k_he4_o16__ne20)*Y(jhe4)*Y(jo16)*dens - &
      screened_rates(k_he4_s32__ar36)*Y(jhe4)*Y(js32)*dens - &
      screened_rates(k_he4_si28__c12_ne20)*Y(jhe4)*Y(jsi28)*dens - &
      screened_rates(k_he4_si28__o16_o16)*Y(jhe4)*Y(jsi28)*dens - &
      screened_rates(k_he4_si28__s32)*Y(jhe4)*Y(jsi28)*dens - &
      screened_rates(k_he4_ti44__cr48)*Y(jhe4)*Y(jti44)*dens + &
      screened_rates(k_mg24__he4_ne20)*Y(jmg24) + screened_rates(k_ne20__he4_o16)* &
      Y(jne20) + screened_rates(k_ni56__he4_fe52)*Y(jni56) + &
      screened_rates(k_o16__he4_c12)*Y(jo16) + 0.5d0* &
      screened_rates(k_o16_o16__he4_si28)*Y(jo16)**2*dens + &
      screened_rates(k_s32__he4_si28)*Y(js32) + screened_rates(k_si28__he4_mg24)* &
      Y(jsi28) + screened_rates(k_ti44__he4_ca40)*Y(jti44) + &
      screened_rates(k_zn60__he4_ni56)*Y(jzn60) &
       )

    ydot_nuc(jc12) = ( &
      -screened_rates(k_c12__he4_he4_he4)*Y(jc12) - screened_rates(k_c12_c12__he4_ne20)* &
      Y(jc12)**2*dens - screened_rates(k_c12_ne20__he4_si28)*Y(jc12)* &
      Y(jne20)*dens - screened_rates(k_c12_o16__he4_mg24)*Y(jc12)*Y(jo16)* &
      dens - screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*dens + &
      0.16666666666666667d0*screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**3*dens &
      **2 + screened_rates(k_he4_mg24__c12_o16)*Y(jhe4)*Y(jmg24)*dens + &
      2.0d0*screened_rates(k_he4_ne20__c12_c12)*Y(jhe4)*Y(jne20)*dens + &
      screened_rates(k_he4_si28__c12_ne20)*Y(jhe4)*Y(jsi28)*dens + &
      screened_rates(k_o16__he4_c12)*Y(jo16) &
       )

    ydot_nuc(jo16) = ( &
      -screened_rates(k_c12_o16__he4_mg24)*Y(jc12)*Y(jo16)*dens + &
      screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*dens + &
      screened_rates(k_he4_mg24__c12_o16)*Y(jhe4)*Y(jmg24)*dens - &
      screened_rates(k_he4_o16__ne20)*Y(jhe4)*Y(jo16)*dens + 2.0d0* &
      screened_rates(k_he4_si28__o16_o16)*Y(jhe4)*Y(jsi28)*dens + &
      screened_rates(k_ne20__he4_o16)*Y(jne20) - screened_rates(k_o16__he4_c12)* &
      Y(jo16) - screened_rates(k_o16_o16__he4_si28)*Y(jo16)**2*dens &
       )

    ydot_nuc(jne20) = ( &
      0.5d0*screened_rates(k_c12_c12__he4_ne20)*Y(jc12)**2*dens - &
      screened_rates(k_c12_ne20__he4_si28)*Y(jc12)*Y(jne20)*dens - &
      screened_rates(k_he4_ne20__c12_c12)*Y(jhe4)*Y(jne20)*dens - &
      screened_rates(k_he4_ne20__mg24)*Y(jhe4)*Y(jne20)*dens + &
      screened_rates(k_he4_o16__ne20)*Y(jhe4)*Y(jo16)*dens + &
      screened_rates(k_he4_si28__c12_ne20)*Y(jhe4)*Y(jsi28)*dens + &
      screened_rates(k_mg24__he4_ne20)*Y(jmg24) - screened_rates(k_ne20__he4_o16)* &
      Y(jne20) &
       )

    ydot_nuc(jmg24) = ( &
      screened_rates(k_c12_o16__he4_mg24)*Y(jc12)*Y(jo16)*dens - &
      screened_rates(k_he4_mg24__c12_o16)*Y(jhe4)*Y(jmg24)*dens - &
      screened_rates(k_he4_mg24__si28)*Y(jhe4)*Y(jmg24)*dens + &
      screened_rates(k_he4_ne20__mg24)*Y(jhe4)*Y(jne20)*dens - &
      screened_rates(k_mg24__he4_ne20)*Y(jmg24) + screened_rates(k_si28__he4_mg24)* &
      Y(jsi28) &
       )

    ydot_nuc(jsi28) = ( &
      screened_rates(k_c12_ne20__he4_si28)*Y(jc12)*Y(jne20)*dens + &
      screened_rates(k_he4_mg24__si28)*Y(jhe4)*Y(jmg24)*dens - &
      screened_rates(k_he4_si28__c12_ne20)*Y(jhe4)*Y(jsi28)*dens - &
      screened_rates(k_he4_si28__o16_o16)*Y(jhe4)*Y(jsi28)*dens - &
      screened_rates(k_he4_si28__s32)*Y(jhe4)*Y(jsi28)*dens + 0.5d0* &
      screened_rates(k_o16_o16__he4_si28)*Y(jo16)**2*dens + &
      screened_rates(k_s32__he4_si28)*Y(js32) - screened_rates(k_si28__he4_mg24)* &
      Y(jsi28) &
       )

    ydot_nuc(js32) = ( &
      screened_rates(k_ar36__he4_s32)*Y(jar36) - screened_rates(k_he4_s32__ar36)*Y(jhe4)* &
      Y(js32)*dens + screened_rates(k_he4_si28__s32)*Y(jhe4)*Y(jsi28)*dens &
      - screened_rates(k_s32__he4_si28)*Y(js32) &
       )

    ydot_nuc(jar36) = ( &
      -screened_rates(k_ar36__he4_s32)*Y(jar36) + screened_rates(k_ca40__he4_ar36)*Y(jca40) &
      - screened_rates(k_he4_ar36__ca40)*Y(jar36)*Y(jhe4)*dens + &
      screened_rates(k_he4_s32__ar36)*Y(jhe4)*Y(js32)*dens &
       )

    ydot_nuc(jca40) = ( &
      -screened_rates(k_ca40__he4_ar36)*Y(jca40) + screened_rates(k_he4_ar36__ca40)*Y(jar36) &
      *Y(jhe4)*dens - screened_rates(k_he4_ca40__ti44)*Y(jca40)*Y(jhe4)* &
      dens + screened_rates(k_ti44__he4_ca40)*Y(jti44) &
       )

    ydot_nuc(jti44) = ( &
      screened_rates(k_cr48__he4_ti44)*Y(jcr48) + screened_rates(k_he4_ca40__ti44)*Y(jca40)* &
      Y(jhe4)*dens - screened_rates(k_he4_ti44__cr48)*Y(jhe4)*Y(jti44)* &
      dens - screened_rates(k_ti44__he4_ca40)*Y(jti44) &
       )

    ydot_nuc(jcr48) = ( &
      -screened_rates(k_cr48__he4_ti44)*Y(jcr48) + screened_rates(k_fe52__he4_cr48)*Y(jfe52) &
      - screened_rates(k_he4_cr48__fe52)*Y(jcr48)*Y(jhe4)*dens + &
      screened_rates(k_he4_ti44__cr48)*Y(jhe4)*Y(jti44)*dens &
       )

    ydot_nuc(jfe52) = ( &
      -screened_rates(k_fe52__he4_cr48)*Y(jfe52) + screened_rates(k_he4_cr48__fe52)*Y(jcr48) &
      *Y(jhe4)*dens - screened_rates(k_he4_fe52__ni56)*Y(jfe52)*Y(jhe4)* &
      dens + screened_rates(k_ni56__he4_fe52)*Y(jni56) &
       )

    ydot_nuc(jni56) = ( &
      screened_rates(k_he4_fe52__ni56)*Y(jfe52)*Y(jhe4)*dens - &
      screened_rates(k_he4_ni56__zn60)*Y(jhe4)*Y(jni56)*dens - &
      screened_rates(k_ni56__he4_fe52)*Y(jni56) + screened_rates(k_zn60__he4_ni56)* &
      Y(jzn60) &
       )

    ydot_nuc(jzn60) = ( &
      screened_rates(k_he4_ni56__zn60)*Y(jhe4)*Y(jni56)*dens - &
      screened_rates(k_zn60__he4_ni56)*Y(jzn60) &
       )


  end subroutine rhs_nuc

  
  subroutine actual_jac(state)

    !$acc routine seq

    use burn_type_module, only: net_itemp, net_ienuc
    
    implicit none
    
    type(burn_t) :: state
    type(rate_eval_t) :: rate_eval
    type(plasma_state) :: pstate
    double precision :: reactvec(num_rate_groups+2)
    double precision :: screened_rates_dt(nrates)
    double precision :: dfdy_nuc(nspec, nspec)
    double precision :: Y(nspec)
    double precision :: dens, temp, ye, rhoy, b1
    double precision :: sneut, dsneutdt, dsneutdd, snuda, snudz
    integer :: i, j

    dens = state%rho
    temp = state%T

    ! Set molar abundances
    Y(:) = state%xn(:) * aion_inv(:)
    
    call evaluate_rates(state, rate_eval)
    
    ! Species Jacobian elements with respect to other species
    call jac_nuc(dfdy_nuc, Y, rate_eval % screened_rates, dens)
    state%jac(1:nspec, 1:nspec) = dfdy_nuc

    ! Species Jacobian elements with respect to energy generation rate
    state%jac(1:nspec, net_ienuc) = 0.0d0

    ! Evaluate the species Jacobian elements with respect to temperature by
    ! calling the RHS using the temperature derivative of the screened rate
    screened_rates_dt = rate_eval % unscreened_rates(i_rate, :) * &
         rate_eval % unscreened_rates(i_dscor_dt, :) + &
         rate_eval % unscreened_rates(i_drate_dt, :) * &
         rate_eval % unscreened_rates(i_scor, :)

    call rhs_nuc(state%jac(1:nspec, net_itemp), Y, screened_rates_dt, dens)
    
    ! Energy generation rate Jacobian elements with respect to species
    do j = 1, nspec
       call ener_gener_rate(state % jac(1:nspec,j), state % jac(net_ienuc,j))
    enddo

    ! Account for the thermal neutrino losses
    call sneut5(temp, dens, state%abar, state%zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)
    do j = 1, nspec
       b1 = ((aion(j) - state%abar) * state%abar * snuda + (zion(j) - state%zbar) * state%abar * snudz)
       state % jac(net_ienuc,j) = state % jac(net_ienuc,j) - b1
    enddo

    ! Energy generation rate Jacobian element with respect to energy generation rate
    state%jac(net_ienuc, net_ienuc) = 0.0d0

    ! Energy generation rate Jacobian element with respect to temperature
    call ener_gener_rate(state%jac(1:nspec, net_itemp), state%jac(net_ienuc, net_itemp))
    state%jac(net_ienuc, net_itemp) = state%jac(net_ienuc, net_itemp) - dsneutdt

    ! Add dqweak and epart contributions!!!

    ! Temperature Jacobian elements
    call temperature_jac(state)

  end subroutine actual_jac

  subroutine jac_nuc(dfdy_nuc, Y, screened_rates, dens)

    !$acc routine seq
    

    double precision, intent(out) :: dfdy_nuc(nspec, nspec)
    double precision, intent(in)  :: Y(nspec)
    double precision, intent(in)  :: screened_rates(nrates)
    double precision, intent(in)  :: dens



    dfdy_nuc(jhe4,jhe4) = ( &
      -screened_rates(k_he4_ar36__ca40)*Y(jar36)*dens - screened_rates(k_he4_c12__o16)* &
      Y(jc12)*dens - screened_rates(k_he4_ca40__ti44)*Y(jca40)*dens - &
      screened_rates(k_he4_cr48__fe52)*Y(jcr48)*dens - screened_rates(k_he4_fe52__ni56) &
      *Y(jfe52)*dens - 1.5d0*screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**2* &
      dens**2 - screened_rates(k_he4_mg24__c12_o16)*Y(jmg24)*dens - &
      screened_rates(k_he4_mg24__si28)*Y(jmg24)*dens - &
      screened_rates(k_he4_ne20__c12_c12)*Y(jne20)*dens - &
      screened_rates(k_he4_ne20__mg24)*Y(jne20)*dens - screened_rates(k_he4_ni56__zn60) &
      *Y(jni56)*dens - screened_rates(k_he4_o16__ne20)*Y(jo16)*dens - &
      screened_rates(k_he4_s32__ar36)*Y(js32)*dens - &
      screened_rates(k_he4_si28__c12_ne20)*Y(jsi28)*dens - &
      screened_rates(k_he4_si28__o16_o16)*Y(jsi28)*dens - &
      screened_rates(k_he4_si28__s32)*Y(jsi28)*dens - screened_rates(k_he4_ti44__cr48)* &
      Y(jti44)*dens &
       )

    dfdy_nuc(jhe4,jc12) = ( &
      3.0d0*screened_rates(k_c12__he4_he4_he4) + 1.0d0*screened_rates(k_c12_c12__he4_ne20)* &
      Y(jc12)*dens + screened_rates(k_c12_ne20__he4_si28)*Y(jne20)*dens + &
      screened_rates(k_c12_o16__he4_mg24)*Y(jo16)*dens - screened_rates(k_he4_c12__o16) &
      *Y(jhe4)*dens &
       )

    dfdy_nuc(jhe4,jo16) = ( &
      screened_rates(k_c12_o16__he4_mg24)*Y(jc12)*dens - screened_rates(k_he4_o16__ne20)* &
      Y(jhe4)*dens + screened_rates(k_o16__he4_c12) + 1.0d0* &
      screened_rates(k_o16_o16__he4_si28)*Y(jo16)*dens &
       )

    dfdy_nuc(jhe4,jne20) = ( &
      screened_rates(k_c12_ne20__he4_si28)*Y(jc12)*dens - screened_rates(k_he4_ne20__c12_c12)* &
      Y(jhe4)*dens - screened_rates(k_he4_ne20__mg24)*Y(jhe4)*dens + &
      screened_rates(k_ne20__he4_o16) &
       )

    dfdy_nuc(jhe4,jmg24) = ( &
      -screened_rates(k_he4_mg24__c12_o16)*Y(jhe4)*dens - screened_rates(k_he4_mg24__si28)* &
      Y(jhe4)*dens + screened_rates(k_mg24__he4_ne20) &
       )

    dfdy_nuc(jhe4,jsi28) = ( &
      -screened_rates(k_he4_si28__c12_ne20)*Y(jhe4)*dens - screened_rates(k_he4_si28__o16_o16) &
      *Y(jhe4)*dens - screened_rates(k_he4_si28__s32)*Y(jhe4)*dens + &
      screened_rates(k_si28__he4_mg24) &
       )

    dfdy_nuc(jhe4,js32) = ( &
      -screened_rates(k_he4_s32__ar36)*Y(jhe4)*dens + screened_rates(k_s32__he4_si28) &
       )

    dfdy_nuc(jhe4,jar36) = ( &
      screened_rates(k_ar36__he4_s32) - screened_rates(k_he4_ar36__ca40)*Y(jhe4)*dens &
       )

    dfdy_nuc(jhe4,jca40) = ( &
      screened_rates(k_ca40__he4_ar36) - screened_rates(k_he4_ca40__ti44)*Y(jhe4)*dens &
       )

    dfdy_nuc(jhe4,jti44) = ( &
      -screened_rates(k_he4_ti44__cr48)*Y(jhe4)*dens + screened_rates(k_ti44__he4_ca40) &
       )

    dfdy_nuc(jhe4,jcr48) = ( &
      screened_rates(k_cr48__he4_ti44) - screened_rates(k_he4_cr48__fe52)*Y(jhe4)*dens &
       )

    dfdy_nuc(jhe4,jfe52) = ( &
      screened_rates(k_fe52__he4_cr48) - screened_rates(k_he4_fe52__ni56)*Y(jhe4)*dens &
       )

    dfdy_nuc(jhe4,jni56) = ( &
      -screened_rates(k_he4_ni56__zn60)*Y(jhe4)*dens + screened_rates(k_ni56__he4_fe52) &
       )

    dfdy_nuc(jhe4,jzn60) = ( &
      screened_rates(k_zn60__he4_ni56) &
       )

    dfdy_nuc(jc12,jhe4) = ( &
      -screened_rates(k_he4_c12__o16)*Y(jc12)*dens + 0.5d0*screened_rates(k_he4_he4_he4__c12)* &
      Y(jhe4)**2*dens**2 + screened_rates(k_he4_mg24__c12_o16)*Y(jmg24)*dens &
      + 2.0d0*screened_rates(k_he4_ne20__c12_c12)*Y(jne20)*dens + &
      screened_rates(k_he4_si28__c12_ne20)*Y(jsi28)*dens &
       )

    dfdy_nuc(jc12,jc12) = ( &
      -screened_rates(k_c12__he4_he4_he4) - 2.0d0*screened_rates(k_c12_c12__he4_ne20)*Y(jc12)* &
      dens - screened_rates(k_c12_ne20__he4_si28)*Y(jne20)*dens - &
      screened_rates(k_c12_o16__he4_mg24)*Y(jo16)*dens - screened_rates(k_he4_c12__o16) &
      *Y(jhe4)*dens &
       )

    dfdy_nuc(jc12,jo16) = ( &
      -screened_rates(k_c12_o16__he4_mg24)*Y(jc12)*dens + screened_rates(k_o16__he4_c12) &
       )

    dfdy_nuc(jc12,jne20) = ( &
      -screened_rates(k_c12_ne20__he4_si28)*Y(jc12)*dens + 2.0d0* &
      screened_rates(k_he4_ne20__c12_c12)*Y(jhe4)*dens &
       )

    dfdy_nuc(jc12,jmg24) = ( &
      screened_rates(k_he4_mg24__c12_o16)*Y(jhe4)*dens &
       )

    dfdy_nuc(jc12,jsi28) = ( &
      screened_rates(k_he4_si28__c12_ne20)*Y(jhe4)*dens &
       )

    dfdy_nuc(jc12,js32) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jar36) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jca40) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jti44) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jcr48) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jfe52) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jni56) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jzn60) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jhe4) = ( &
      screened_rates(k_he4_c12__o16)*Y(jc12)*dens + screened_rates(k_he4_mg24__c12_o16)* &
      Y(jmg24)*dens - screened_rates(k_he4_o16__ne20)*Y(jo16)*dens + 2.0d0* &
      screened_rates(k_he4_si28__o16_o16)*Y(jsi28)*dens &
       )

    dfdy_nuc(jo16,jc12) = ( &
      -screened_rates(k_c12_o16__he4_mg24)*Y(jo16)*dens + screened_rates(k_he4_c12__o16)* &
      Y(jhe4)*dens &
       )

    dfdy_nuc(jo16,jo16) = ( &
      -screened_rates(k_c12_o16__he4_mg24)*Y(jc12)*dens - screened_rates(k_he4_o16__ne20)* &
      Y(jhe4)*dens - screened_rates(k_o16__he4_c12) - 2.0d0* &
      screened_rates(k_o16_o16__he4_si28)*Y(jo16)*dens &
       )

    dfdy_nuc(jo16,jne20) = ( &
      screened_rates(k_ne20__he4_o16) &
       )

    dfdy_nuc(jo16,jmg24) = ( &
      screened_rates(k_he4_mg24__c12_o16)*Y(jhe4)*dens &
       )

    dfdy_nuc(jo16,jsi28) = ( &
      2.0d0*screened_rates(k_he4_si28__o16_o16)*Y(jhe4)*dens &
       )

    dfdy_nuc(jo16,js32) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jar36) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jca40) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jti44) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jcr48) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jfe52) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jni56) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jzn60) = ( &
      0.0d0 &
       )

    dfdy_nuc(jne20,jhe4) = ( &
      -screened_rates(k_he4_ne20__c12_c12)*Y(jne20)*dens - screened_rates(k_he4_ne20__mg24)* &
      Y(jne20)*dens + screened_rates(k_he4_o16__ne20)*Y(jo16)*dens + &
      screened_rates(k_he4_si28__c12_ne20)*Y(jsi28)*dens &
       )

    dfdy_nuc(jne20,jc12) = ( &
      1.0d0*screened_rates(k_c12_c12__he4_ne20)*Y(jc12)*dens - &
      screened_rates(k_c12_ne20__he4_si28)*Y(jne20)*dens &
       )

    dfdy_nuc(jne20,jo16) = ( &
      screened_rates(k_he4_o16__ne20)*Y(jhe4)*dens &
       )

    dfdy_nuc(jne20,jne20) = ( &
      -screened_rates(k_c12_ne20__he4_si28)*Y(jc12)*dens - screened_rates(k_he4_ne20__c12_c12) &
      *Y(jhe4)*dens - screened_rates(k_he4_ne20__mg24)*Y(jhe4)*dens - &
      screened_rates(k_ne20__he4_o16) &
       )

    dfdy_nuc(jne20,jmg24) = ( &
      screened_rates(k_mg24__he4_ne20) &
       )

    dfdy_nuc(jne20,jsi28) = ( &
      screened_rates(k_he4_si28__c12_ne20)*Y(jhe4)*dens &
       )

    dfdy_nuc(jne20,js32) = ( &
      0.0d0 &
       )

    dfdy_nuc(jne20,jar36) = ( &
      0.0d0 &
       )

    dfdy_nuc(jne20,jca40) = ( &
      0.0d0 &
       )

    dfdy_nuc(jne20,jti44) = ( &
      0.0d0 &
       )

    dfdy_nuc(jne20,jcr48) = ( &
      0.0d0 &
       )

    dfdy_nuc(jne20,jfe52) = ( &
      0.0d0 &
       )

    dfdy_nuc(jne20,jni56) = ( &
      0.0d0 &
       )

    dfdy_nuc(jne20,jzn60) = ( &
      0.0d0 &
       )

    dfdy_nuc(jmg24,jhe4) = ( &
      -screened_rates(k_he4_mg24__c12_o16)*Y(jmg24)*dens - screened_rates(k_he4_mg24__si28)* &
      Y(jmg24)*dens + screened_rates(k_he4_ne20__mg24)*Y(jne20)*dens &
       )

    dfdy_nuc(jmg24,jc12) = ( &
      screened_rates(k_c12_o16__he4_mg24)*Y(jo16)*dens &
       )

    dfdy_nuc(jmg24,jo16) = ( &
      screened_rates(k_c12_o16__he4_mg24)*Y(jc12)*dens &
       )

    dfdy_nuc(jmg24,jne20) = ( &
      screened_rates(k_he4_ne20__mg24)*Y(jhe4)*dens &
       )

    dfdy_nuc(jmg24,jmg24) = ( &
      -screened_rates(k_he4_mg24__c12_o16)*Y(jhe4)*dens - screened_rates(k_he4_mg24__si28)* &
      Y(jhe4)*dens - screened_rates(k_mg24__he4_ne20) &
       )

    dfdy_nuc(jmg24,jsi28) = ( &
      screened_rates(k_si28__he4_mg24) &
       )

    dfdy_nuc(jmg24,js32) = ( &
      0.0d0 &
       )

    dfdy_nuc(jmg24,jar36) = ( &
      0.0d0 &
       )

    dfdy_nuc(jmg24,jca40) = ( &
      0.0d0 &
       )

    dfdy_nuc(jmg24,jti44) = ( &
      0.0d0 &
       )

    dfdy_nuc(jmg24,jcr48) = ( &
      0.0d0 &
       )

    dfdy_nuc(jmg24,jfe52) = ( &
      0.0d0 &
       )

    dfdy_nuc(jmg24,jni56) = ( &
      0.0d0 &
       )

    dfdy_nuc(jmg24,jzn60) = ( &
      0.0d0 &
       )

    dfdy_nuc(jsi28,jhe4) = ( &
      screened_rates(k_he4_mg24__si28)*Y(jmg24)*dens - screened_rates(k_he4_si28__c12_ne20)* &
      Y(jsi28)*dens - screened_rates(k_he4_si28__o16_o16)*Y(jsi28)*dens - &
      screened_rates(k_he4_si28__s32)*Y(jsi28)*dens &
       )

    dfdy_nuc(jsi28,jc12) = ( &
      screened_rates(k_c12_ne20__he4_si28)*Y(jne20)*dens &
       )

    dfdy_nuc(jsi28,jo16) = ( &
      1.0d0*screened_rates(k_o16_o16__he4_si28)*Y(jo16)*dens &
       )

    dfdy_nuc(jsi28,jne20) = ( &
      screened_rates(k_c12_ne20__he4_si28)*Y(jc12)*dens &
       )

    dfdy_nuc(jsi28,jmg24) = ( &
      screened_rates(k_he4_mg24__si28)*Y(jhe4)*dens &
       )

    dfdy_nuc(jsi28,jsi28) = ( &
      -screened_rates(k_he4_si28__c12_ne20)*Y(jhe4)*dens - screened_rates(k_he4_si28__o16_o16) &
      *Y(jhe4)*dens - screened_rates(k_he4_si28__s32)*Y(jhe4)*dens - &
      screened_rates(k_si28__he4_mg24) &
       )

    dfdy_nuc(jsi28,js32) = ( &
      screened_rates(k_s32__he4_si28) &
       )

    dfdy_nuc(jsi28,jar36) = ( &
      0.0d0 &
       )

    dfdy_nuc(jsi28,jca40) = ( &
      0.0d0 &
       )

    dfdy_nuc(jsi28,jti44) = ( &
      0.0d0 &
       )

    dfdy_nuc(jsi28,jcr48) = ( &
      0.0d0 &
       )

    dfdy_nuc(jsi28,jfe52) = ( &
      0.0d0 &
       )

    dfdy_nuc(jsi28,jni56) = ( &
      0.0d0 &
       )

    dfdy_nuc(jsi28,jzn60) = ( &
      0.0d0 &
       )

    dfdy_nuc(js32,jhe4) = ( &
      -screened_rates(k_he4_s32__ar36)*Y(js32)*dens + screened_rates(k_he4_si28__s32)* &
      Y(jsi28)*dens &
       )

    dfdy_nuc(js32,jc12) = ( &
      0.0d0 &
       )

    dfdy_nuc(js32,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(js32,jne20) = ( &
      0.0d0 &
       )

    dfdy_nuc(js32,jmg24) = ( &
      0.0d0 &
       )

    dfdy_nuc(js32,jsi28) = ( &
      screened_rates(k_he4_si28__s32)*Y(jhe4)*dens &
       )

    dfdy_nuc(js32,js32) = ( &
      -screened_rates(k_he4_s32__ar36)*Y(jhe4)*dens - screened_rates(k_s32__he4_si28) &
       )

    dfdy_nuc(js32,jar36) = ( &
      screened_rates(k_ar36__he4_s32) &
       )

    dfdy_nuc(js32,jca40) = ( &
      0.0d0 &
       )

    dfdy_nuc(js32,jti44) = ( &
      0.0d0 &
       )

    dfdy_nuc(js32,jcr48) = ( &
      0.0d0 &
       )

    dfdy_nuc(js32,jfe52) = ( &
      0.0d0 &
       )

    dfdy_nuc(js32,jni56) = ( &
      0.0d0 &
       )

    dfdy_nuc(js32,jzn60) = ( &
      0.0d0 &
       )

    dfdy_nuc(jar36,jhe4) = ( &
      -screened_rates(k_he4_ar36__ca40)*Y(jar36)*dens + screened_rates(k_he4_s32__ar36)* &
      Y(js32)*dens &
       )

    dfdy_nuc(jar36,jc12) = ( &
      0.0d0 &
       )

    dfdy_nuc(jar36,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jar36,jne20) = ( &
      0.0d0 &
       )

    dfdy_nuc(jar36,jmg24) = ( &
      0.0d0 &
       )

    dfdy_nuc(jar36,jsi28) = ( &
      0.0d0 &
       )

    dfdy_nuc(jar36,js32) = ( &
      screened_rates(k_he4_s32__ar36)*Y(jhe4)*dens &
       )

    dfdy_nuc(jar36,jar36) = ( &
      -screened_rates(k_ar36__he4_s32) - screened_rates(k_he4_ar36__ca40)*Y(jhe4)*dens &
       )

    dfdy_nuc(jar36,jca40) = ( &
      screened_rates(k_ca40__he4_ar36) &
       )

    dfdy_nuc(jar36,jti44) = ( &
      0.0d0 &
       )

    dfdy_nuc(jar36,jcr48) = ( &
      0.0d0 &
       )

    dfdy_nuc(jar36,jfe52) = ( &
      0.0d0 &
       )

    dfdy_nuc(jar36,jni56) = ( &
      0.0d0 &
       )

    dfdy_nuc(jar36,jzn60) = ( &
      0.0d0 &
       )

    dfdy_nuc(jca40,jhe4) = ( &
      screened_rates(k_he4_ar36__ca40)*Y(jar36)*dens - screened_rates(k_he4_ca40__ti44)* &
      Y(jca40)*dens &
       )

    dfdy_nuc(jca40,jc12) = ( &
      0.0d0 &
       )

    dfdy_nuc(jca40,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jca40,jne20) = ( &
      0.0d0 &
       )

    dfdy_nuc(jca40,jmg24) = ( &
      0.0d0 &
       )

    dfdy_nuc(jca40,jsi28) = ( &
      0.0d0 &
       )

    dfdy_nuc(jca40,js32) = ( &
      0.0d0 &
       )

    dfdy_nuc(jca40,jar36) = ( &
      screened_rates(k_he4_ar36__ca40)*Y(jhe4)*dens &
       )

    dfdy_nuc(jca40,jca40) = ( &
      -screened_rates(k_ca40__he4_ar36) - screened_rates(k_he4_ca40__ti44)*Y(jhe4)*dens &
       )

    dfdy_nuc(jca40,jti44) = ( &
      screened_rates(k_ti44__he4_ca40) &
       )

    dfdy_nuc(jca40,jcr48) = ( &
      0.0d0 &
       )

    dfdy_nuc(jca40,jfe52) = ( &
      0.0d0 &
       )

    dfdy_nuc(jca40,jni56) = ( &
      0.0d0 &
       )

    dfdy_nuc(jca40,jzn60) = ( &
      0.0d0 &
       )

    dfdy_nuc(jti44,jhe4) = ( &
      screened_rates(k_he4_ca40__ti44)*Y(jca40)*dens - screened_rates(k_he4_ti44__cr48)* &
      Y(jti44)*dens &
       )

    dfdy_nuc(jti44,jc12) = ( &
      0.0d0 &
       )

    dfdy_nuc(jti44,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jti44,jne20) = ( &
      0.0d0 &
       )

    dfdy_nuc(jti44,jmg24) = ( &
      0.0d0 &
       )

    dfdy_nuc(jti44,jsi28) = ( &
      0.0d0 &
       )

    dfdy_nuc(jti44,js32) = ( &
      0.0d0 &
       )

    dfdy_nuc(jti44,jar36) = ( &
      0.0d0 &
       )

    dfdy_nuc(jti44,jca40) = ( &
      screened_rates(k_he4_ca40__ti44)*Y(jhe4)*dens &
       )

    dfdy_nuc(jti44,jti44) = ( &
      -screened_rates(k_he4_ti44__cr48)*Y(jhe4)*dens - screened_rates(k_ti44__he4_ca40) &
       )

    dfdy_nuc(jti44,jcr48) = ( &
      screened_rates(k_cr48__he4_ti44) &
       )

    dfdy_nuc(jti44,jfe52) = ( &
      0.0d0 &
       )

    dfdy_nuc(jti44,jni56) = ( &
      0.0d0 &
       )

    dfdy_nuc(jti44,jzn60) = ( &
      0.0d0 &
       )

    dfdy_nuc(jcr48,jhe4) = ( &
      -screened_rates(k_he4_cr48__fe52)*Y(jcr48)*dens + screened_rates(k_he4_ti44__cr48)* &
      Y(jti44)*dens &
       )

    dfdy_nuc(jcr48,jc12) = ( &
      0.0d0 &
       )

    dfdy_nuc(jcr48,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jcr48,jne20) = ( &
      0.0d0 &
       )

    dfdy_nuc(jcr48,jmg24) = ( &
      0.0d0 &
       )

    dfdy_nuc(jcr48,jsi28) = ( &
      0.0d0 &
       )

    dfdy_nuc(jcr48,js32) = ( &
      0.0d0 &
       )

    dfdy_nuc(jcr48,jar36) = ( &
      0.0d0 &
       )

    dfdy_nuc(jcr48,jca40) = ( &
      0.0d0 &
       )

    dfdy_nuc(jcr48,jti44) = ( &
      screened_rates(k_he4_ti44__cr48)*Y(jhe4)*dens &
       )

    dfdy_nuc(jcr48,jcr48) = ( &
      -screened_rates(k_cr48__he4_ti44) - screened_rates(k_he4_cr48__fe52)*Y(jhe4)*dens &
       )

    dfdy_nuc(jcr48,jfe52) = ( &
      screened_rates(k_fe52__he4_cr48) &
       )

    dfdy_nuc(jcr48,jni56) = ( &
      0.0d0 &
       )

    dfdy_nuc(jcr48,jzn60) = ( &
      0.0d0 &
       )

    dfdy_nuc(jfe52,jhe4) = ( &
      screened_rates(k_he4_cr48__fe52)*Y(jcr48)*dens - screened_rates(k_he4_fe52__ni56)* &
      Y(jfe52)*dens &
       )

    dfdy_nuc(jfe52,jc12) = ( &
      0.0d0 &
       )

    dfdy_nuc(jfe52,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jfe52,jne20) = ( &
      0.0d0 &
       )

    dfdy_nuc(jfe52,jmg24) = ( &
      0.0d0 &
       )

    dfdy_nuc(jfe52,jsi28) = ( &
      0.0d0 &
       )

    dfdy_nuc(jfe52,js32) = ( &
      0.0d0 &
       )

    dfdy_nuc(jfe52,jar36) = ( &
      0.0d0 &
       )

    dfdy_nuc(jfe52,jca40) = ( &
      0.0d0 &
       )

    dfdy_nuc(jfe52,jti44) = ( &
      0.0d0 &
       )

    dfdy_nuc(jfe52,jcr48) = ( &
      screened_rates(k_he4_cr48__fe52)*Y(jhe4)*dens &
       )

    dfdy_nuc(jfe52,jfe52) = ( &
      -screened_rates(k_fe52__he4_cr48) - screened_rates(k_he4_fe52__ni56)*Y(jhe4)*dens &
       )

    dfdy_nuc(jfe52,jni56) = ( &
      screened_rates(k_ni56__he4_fe52) &
       )

    dfdy_nuc(jfe52,jzn60) = ( &
      0.0d0 &
       )

    dfdy_nuc(jni56,jhe4) = ( &
      screened_rates(k_he4_fe52__ni56)*Y(jfe52)*dens - screened_rates(k_he4_ni56__zn60)* &
      Y(jni56)*dens &
       )

    dfdy_nuc(jni56,jc12) = ( &
      0.0d0 &
       )

    dfdy_nuc(jni56,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jni56,jne20) = ( &
      0.0d0 &
       )

    dfdy_nuc(jni56,jmg24) = ( &
      0.0d0 &
       )

    dfdy_nuc(jni56,jsi28) = ( &
      0.0d0 &
       )

    dfdy_nuc(jni56,js32) = ( &
      0.0d0 &
       )

    dfdy_nuc(jni56,jar36) = ( &
      0.0d0 &
       )

    dfdy_nuc(jni56,jca40) = ( &
      0.0d0 &
       )

    dfdy_nuc(jni56,jti44) = ( &
      0.0d0 &
       )

    dfdy_nuc(jni56,jcr48) = ( &
      0.0d0 &
       )

    dfdy_nuc(jni56,jfe52) = ( &
      screened_rates(k_he4_fe52__ni56)*Y(jhe4)*dens &
       )

    dfdy_nuc(jni56,jni56) = ( &
      -screened_rates(k_he4_ni56__zn60)*Y(jhe4)*dens - screened_rates(k_ni56__he4_fe52) &
       )

    dfdy_nuc(jni56,jzn60) = ( &
      screened_rates(k_zn60__he4_ni56) &
       )

    dfdy_nuc(jzn60,jhe4) = ( &
      screened_rates(k_he4_ni56__zn60)*Y(jni56)*dens &
       )

    dfdy_nuc(jzn60,jc12) = ( &
      0.0d0 &
       )

    dfdy_nuc(jzn60,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jzn60,jne20) = ( &
      0.0d0 &
       )

    dfdy_nuc(jzn60,jmg24) = ( &
      0.0d0 &
       )

    dfdy_nuc(jzn60,jsi28) = ( &
      0.0d0 &
       )

    dfdy_nuc(jzn60,js32) = ( &
      0.0d0 &
       )

    dfdy_nuc(jzn60,jar36) = ( &
      0.0d0 &
       )

    dfdy_nuc(jzn60,jca40) = ( &
      0.0d0 &
       )

    dfdy_nuc(jzn60,jti44) = ( &
      0.0d0 &
       )

    dfdy_nuc(jzn60,jcr48) = ( &
      0.0d0 &
       )

    dfdy_nuc(jzn60,jfe52) = ( &
      0.0d0 &
       )

    dfdy_nuc(jzn60,jni56) = ( &
      screened_rates(k_he4_ni56__zn60)*Y(jhe4)*dens &
       )

    dfdy_nuc(jzn60,jzn60) = ( &
      -screened_rates(k_zn60__he4_ni56) &
       )

    
  end subroutine jac_nuc

end module actual_rhs_module
