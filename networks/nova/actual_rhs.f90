module actual_rhs_module

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module
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



    ydot_nuc(jp) = ( &
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*Y(jn13)*dens + &
      screened_rates(k_he4_o14__p_f17)*Y(jhe4)*Y(jo14)*dens - &
      screened_rates(k_p_c12__n13)*Y(jc12)*Y(jp)*dens - screened_rates(k_p_c13__n14)* &
      Y(jc13)*Y(jp)*dens - screened_rates(k_p_f18__he4_o15)*Y(jf18)* &
      Y(jp)*dens - screened_rates(k_p_n13__o14)*Y(jn13)*Y(jp)*dens - &
      screened_rates(k_p_n14__o15)*Y(jn14)*Y(jp)*dens - &
      screened_rates(k_p_n15__he4_c12)*Y(jn15)*Y(jp)*dens - &
      screened_rates(k_p_n15__o16)*Y(jn15)*Y(jp)*dens - screened_rates(k_p_o16__f17)* &
      Y(jo16)*Y(jp)*dens - screened_rates(k_p_o17__f18)*Y(jo17)*Y(jp)* &
      dens - screened_rates(k_p_o17__he4_n14)*Y(jo17)*Y(jp)*dens &
       )

    ydot_nuc(jhe4) = ( &
      -screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*dens - 0.5d0* &
      screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**3*dens**2 - &
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*Y(jn13)*dens - &
      screened_rates(k_he4_n14__f18)*Y(jhe4)*Y(jn14)*dens - &
      screened_rates(k_he4_o14__p_f17)*Y(jhe4)*Y(jo14)*dens + &
      screened_rates(k_p_f18__he4_o15)*Y(jf18)*Y(jp)*dens + &
      screened_rates(k_p_n15__he4_c12)*Y(jn15)*Y(jp)*dens + &
      screened_rates(k_p_o17__he4_n14)*Y(jo17)*Y(jp)*dens &
       )

    ydot_nuc(jc12) = ( &
      -screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*dens + 0.16666666666666667d0* &
      screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**3*dens**2 - &
      screened_rates(k_p_c12__n13)*Y(jc12)*Y(jp)*dens + &
      screened_rates(k_p_n15__he4_c12)*Y(jn15)*Y(jp)*dens &
       )

    ydot_nuc(jc13) = ( &
      screened_rates(k_n13__c13)*Y(jn13) - screened_rates(k_p_c13__n14)*Y(jc13)*Y(jp)*dens &
       )

    ydot_nuc(jn13) = ( &
      -screened_rates(k_he4_n13__p_o16)*Y(jhe4)*Y(jn13)*dens - screened_rates(k_n13__c13)* &
      Y(jn13) + screened_rates(k_p_c12__n13)*Y(jc12)*Y(jp)*dens - &
      screened_rates(k_p_n13__o14)*Y(jn13)*Y(jp)*dens &
       )

    ydot_nuc(jn14) = ( &
      -screened_rates(k_he4_n14__f18)*Y(jhe4)*Y(jn14)*dens + screened_rates(k_o14__n14)* &
      Y(jo14) + screened_rates(k_p_c13__n14)*Y(jc13)*Y(jp)*dens - &
      screened_rates(k_p_n14__o15)*Y(jn14)*Y(jp)*dens + &
      screened_rates(k_p_o17__he4_n14)*Y(jo17)*Y(jp)*dens &
       )

    ydot_nuc(jn15) = ( &
      screened_rates(k_o15__n15)*Y(jo15) - screened_rates(k_p_n15__he4_c12)*Y(jn15)*Y(jp)* &
      dens - screened_rates(k_p_n15__o16)*Y(jn15)*Y(jp)*dens &
       )

    ydot_nuc(jo14) = ( &
      -screened_rates(k_he4_o14__p_f17)*Y(jhe4)*Y(jo14)*dens - screened_rates(k_o14__n14)* &
      Y(jo14) + screened_rates(k_p_n13__o14)*Y(jn13)*Y(jp)*dens &
       )

    ydot_nuc(jo15) = ( &
      -screened_rates(k_o15__n15)*Y(jo15) + screened_rates(k_p_f18__he4_o15)*Y(jf18)*Y(jp) &
      *dens + screened_rates(k_p_n14__o15)*Y(jn14)*Y(jp)*dens &
       )

    ydot_nuc(jo16) = ( &
      screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*dens + screened_rates(k_he4_n13__p_o16) &
      *Y(jhe4)*Y(jn13)*dens + screened_rates(k_p_n15__o16)*Y(jn15)*Y(jp) &
      *dens - screened_rates(k_p_o16__f17)*Y(jo16)*Y(jp)*dens &
       )

    ydot_nuc(jo17) = ( &
      screened_rates(k_f17__o17)*Y(jf17) - screened_rates(k_p_o17__f18)*Y(jo17)*Y(jp)*dens &
      - screened_rates(k_p_o17__he4_n14)*Y(jo17)*Y(jp)*dens &
       )

    ydot_nuc(jf17) = ( &
      -screened_rates(k_f17__o17)*Y(jf17) + screened_rates(k_he4_o14__p_f17)*Y(jhe4)* &
      Y(jo14)*dens + screened_rates(k_p_o16__f17)*Y(jo16)*Y(jp)*dens &
       )

    ydot_nuc(jf18) = ( &
      screened_rates(k_he4_n14__f18)*Y(jhe4)*Y(jn14)*dens - screened_rates(k_p_f18__he4_o15) &
      *Y(jf18)*Y(jp)*dens + screened_rates(k_p_o17__f18)*Y(jo17)*Y(jp)* &
      dens &
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



    dfdy_nuc(jp,jp) = ( &
      -screened_rates(k_p_c12__n13)*Y(jc12)*dens - screened_rates(k_p_c13__n14)*Y(jc13)*dens &
      - screened_rates(k_p_f18__he4_o15)*Y(jf18)*dens - screened_rates(k_p_n13__o14)* &
      Y(jn13)*dens - screened_rates(k_p_n14__o15)*Y(jn14)*dens - &
      screened_rates(k_p_n15__he4_c12)*Y(jn15)*dens - screened_rates(k_p_n15__o16)* &
      Y(jn15)*dens - screened_rates(k_p_o16__f17)*Y(jo16)*dens - &
      screened_rates(k_p_o17__f18)*Y(jo17)*dens - screened_rates(k_p_o17__he4_n14)* &
      Y(jo17)*dens &
       )

    dfdy_nuc(jp,jhe4) = ( &
      screened_rates(k_he4_n13__p_o16)*Y(jn13)*dens + screened_rates(k_he4_o14__p_f17)* &
      Y(jo14)*dens &
       )

    dfdy_nuc(jp,jc12) = ( &
      -screened_rates(k_p_c12__n13)*Y(jp)*dens &
       )

    dfdy_nuc(jp,jc13) = ( &
      -screened_rates(k_p_c13__n14)*Y(jp)*dens &
       )

    dfdy_nuc(jp,jn13) = ( &
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*dens - screened_rates(k_p_n13__o14)*Y(jp)* &
      dens &
       )

    dfdy_nuc(jp,jn14) = ( &
      -screened_rates(k_p_n14__o15)*Y(jp)*dens &
       )

    dfdy_nuc(jp,jn15) = ( &
      -screened_rates(k_p_n15__he4_c12)*Y(jp)*dens - screened_rates(k_p_n15__o16)*Y(jp)*dens &
       )

    dfdy_nuc(jp,jo14) = ( &
      screened_rates(k_he4_o14__p_f17)*Y(jhe4)*dens &
       )

    dfdy_nuc(jp,jo15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jp,jo16) = ( &
      -screened_rates(k_p_o16__f17)*Y(jp)*dens &
       )

    dfdy_nuc(jp,jo17) = ( &
      -screened_rates(k_p_o17__f18)*Y(jp)*dens - screened_rates(k_p_o17__he4_n14)*Y(jp)*dens &
       )

    dfdy_nuc(jp,jf17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jp,jf18) = ( &
      -screened_rates(k_p_f18__he4_o15)*Y(jp)*dens &
       )

    dfdy_nuc(jhe4,jp) = ( &
      screened_rates(k_p_f18__he4_o15)*Y(jf18)*dens + screened_rates(k_p_n15__he4_c12)* &
      Y(jn15)*dens + screened_rates(k_p_o17__he4_n14)*Y(jo17)*dens &
       )

    dfdy_nuc(jhe4,jhe4) = ( &
      -screened_rates(k_he4_c12__o16)*Y(jc12)*dens - 1.5d0*screened_rates(k_he4_he4_he4__c12)* &
      Y(jhe4)**2*dens**2 - screened_rates(k_he4_n13__p_o16)*Y(jn13)*dens - &
      screened_rates(k_he4_n14__f18)*Y(jn14)*dens - screened_rates(k_he4_o14__p_f17)* &
      Y(jo14)*dens &
       )

    dfdy_nuc(jhe4,jc12) = ( &
      -screened_rates(k_he4_c12__o16)*Y(jhe4)*dens &
       )

    dfdy_nuc(jhe4,jc13) = ( &
      0.0d0 &
       )

    dfdy_nuc(jhe4,jn13) = ( &
      -screened_rates(k_he4_n13__p_o16)*Y(jhe4)*dens &
       )

    dfdy_nuc(jhe4,jn14) = ( &
      -screened_rates(k_he4_n14__f18)*Y(jhe4)*dens &
       )

    dfdy_nuc(jhe4,jn15) = ( &
      screened_rates(k_p_n15__he4_c12)*Y(jp)*dens &
       )

    dfdy_nuc(jhe4,jo14) = ( &
      -screened_rates(k_he4_o14__p_f17)*Y(jhe4)*dens &
       )

    dfdy_nuc(jhe4,jo15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jhe4,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jhe4,jo17) = ( &
      screened_rates(k_p_o17__he4_n14)*Y(jp)*dens &
       )

    dfdy_nuc(jhe4,jf17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jhe4,jf18) = ( &
      screened_rates(k_p_f18__he4_o15)*Y(jp)*dens &
       )

    dfdy_nuc(jc12,jp) = ( &
      -screened_rates(k_p_c12__n13)*Y(jc12)*dens + screened_rates(k_p_n15__he4_c12)*Y(jn15)* &
      dens &
       )

    dfdy_nuc(jc12,jhe4) = ( &
      -screened_rates(k_he4_c12__o16)*Y(jc12)*dens + 0.5d0*screened_rates(k_he4_he4_he4__c12)* &
      Y(jhe4)**2*dens**2 &
       )

    dfdy_nuc(jc12,jc12) = ( &
      -screened_rates(k_he4_c12__o16)*Y(jhe4)*dens - screened_rates(k_p_c12__n13)*Y(jp)*dens &
       )

    dfdy_nuc(jc12,jc13) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jn13) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jn14) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jn15) = ( &
      screened_rates(k_p_n15__he4_c12)*Y(jp)*dens &
       )

    dfdy_nuc(jc12,jo14) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jo15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jo17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jf17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jf18) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc13,jp) = ( &
      -screened_rates(k_p_c13__n14)*Y(jc13)*dens &
       )

    dfdy_nuc(jc13,jhe4) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc13,jc12) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc13,jc13) = ( &
      -screened_rates(k_p_c13__n14)*Y(jp)*dens &
       )

    dfdy_nuc(jc13,jn13) = ( &
      screened_rates(k_n13__c13) &
       )

    dfdy_nuc(jc13,jn14) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc13,jn15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc13,jo14) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc13,jo15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc13,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc13,jo17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc13,jf17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc13,jf18) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn13,jp) = ( &
      screened_rates(k_p_c12__n13)*Y(jc12)*dens - screened_rates(k_p_n13__o14)*Y(jn13)*dens &
       )

    dfdy_nuc(jn13,jhe4) = ( &
      -screened_rates(k_he4_n13__p_o16)*Y(jn13)*dens &
       )

    dfdy_nuc(jn13,jc12) = ( &
      screened_rates(k_p_c12__n13)*Y(jp)*dens &
       )

    dfdy_nuc(jn13,jc13) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn13,jn13) = ( &
      -screened_rates(k_he4_n13__p_o16)*Y(jhe4)*dens - screened_rates(k_n13__c13) - &
      screened_rates(k_p_n13__o14)*Y(jp)*dens &
       )

    dfdy_nuc(jn13,jn14) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn13,jn15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn13,jo14) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn13,jo15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn13,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn13,jo17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn13,jf17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn13,jf18) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn14,jp) = ( &
      screened_rates(k_p_c13__n14)*Y(jc13)*dens - screened_rates(k_p_n14__o15)*Y(jn14)*dens &
      + screened_rates(k_p_o17__he4_n14)*Y(jo17)*dens &
       )

    dfdy_nuc(jn14,jhe4) = ( &
      -screened_rates(k_he4_n14__f18)*Y(jn14)*dens &
       )

    dfdy_nuc(jn14,jc12) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn14,jc13) = ( &
      screened_rates(k_p_c13__n14)*Y(jp)*dens &
       )

    dfdy_nuc(jn14,jn13) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn14,jn14) = ( &
      -screened_rates(k_he4_n14__f18)*Y(jhe4)*dens - screened_rates(k_p_n14__o15)*Y(jp)*dens &
       )

    dfdy_nuc(jn14,jn15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn14,jo14) = ( &
      screened_rates(k_o14__n14) &
       )

    dfdy_nuc(jn14,jo15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn14,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn14,jo17) = ( &
      screened_rates(k_p_o17__he4_n14)*Y(jp)*dens &
       )

    dfdy_nuc(jn14,jf17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn14,jf18) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn15,jp) = ( &
      -screened_rates(k_p_n15__he4_c12)*Y(jn15)*dens - screened_rates(k_p_n15__o16)*Y(jn15)* &
      dens &
       )

    dfdy_nuc(jn15,jhe4) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn15,jc12) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn15,jc13) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn15,jn13) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn15,jn14) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn15,jn15) = ( &
      -screened_rates(k_p_n15__he4_c12)*Y(jp)*dens - screened_rates(k_p_n15__o16)*Y(jp)*dens &
       )

    dfdy_nuc(jn15,jo14) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn15,jo15) = ( &
      screened_rates(k_o15__n15) &
       )

    dfdy_nuc(jn15,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn15,jo17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn15,jf17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn15,jf18) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo14,jp) = ( &
      screened_rates(k_p_n13__o14)*Y(jn13)*dens &
       )

    dfdy_nuc(jo14,jhe4) = ( &
      -screened_rates(k_he4_o14__p_f17)*Y(jo14)*dens &
       )

    dfdy_nuc(jo14,jc12) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo14,jc13) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo14,jn13) = ( &
      screened_rates(k_p_n13__o14)*Y(jp)*dens &
       )

    dfdy_nuc(jo14,jn14) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo14,jn15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo14,jo14) = ( &
      -screened_rates(k_he4_o14__p_f17)*Y(jhe4)*dens - screened_rates(k_o14__n14) &
       )

    dfdy_nuc(jo14,jo15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo14,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo14,jo17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo14,jf17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo14,jf18) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo15,jp) = ( &
      screened_rates(k_p_f18__he4_o15)*Y(jf18)*dens + screened_rates(k_p_n14__o15)*Y(jn14)* &
      dens &
       )

    dfdy_nuc(jo15,jhe4) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo15,jc12) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo15,jc13) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo15,jn13) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo15,jn14) = ( &
      screened_rates(k_p_n14__o15)*Y(jp)*dens &
       )

    dfdy_nuc(jo15,jn15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo15,jo14) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo15,jo15) = ( &
      -screened_rates(k_o15__n15) &
       )

    dfdy_nuc(jo15,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo15,jo17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo15,jf17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo15,jf18) = ( &
      screened_rates(k_p_f18__he4_o15)*Y(jp)*dens &
       )

    dfdy_nuc(jo16,jp) = ( &
      screened_rates(k_p_n15__o16)*Y(jn15)*dens - screened_rates(k_p_o16__f17)*Y(jo16)*dens &
       )

    dfdy_nuc(jo16,jhe4) = ( &
      screened_rates(k_he4_c12__o16)*Y(jc12)*dens + screened_rates(k_he4_n13__p_o16)*Y(jn13) &
      *dens &
       )

    dfdy_nuc(jo16,jc12) = ( &
      screened_rates(k_he4_c12__o16)*Y(jhe4)*dens &
       )

    dfdy_nuc(jo16,jc13) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jn13) = ( &
      screened_rates(k_he4_n13__p_o16)*Y(jhe4)*dens &
       )

    dfdy_nuc(jo16,jn14) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jn15) = ( &
      screened_rates(k_p_n15__o16)*Y(jp)*dens &
       )

    dfdy_nuc(jo16,jo14) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jo15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jo16) = ( &
      -screened_rates(k_p_o16__f17)*Y(jp)*dens &
       )

    dfdy_nuc(jo16,jo17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jf17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jf18) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo17,jp) = ( &
      -screened_rates(k_p_o17__f18)*Y(jo17)*dens - screened_rates(k_p_o17__he4_n14)*Y(jo17)* &
      dens &
       )

    dfdy_nuc(jo17,jhe4) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo17,jc12) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo17,jc13) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo17,jn13) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo17,jn14) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo17,jn15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo17,jo14) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo17,jo15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo17,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo17,jo17) = ( &
      -screened_rates(k_p_o17__f18)*Y(jp)*dens - screened_rates(k_p_o17__he4_n14)*Y(jp)*dens &
       )

    dfdy_nuc(jo17,jf17) = ( &
      screened_rates(k_f17__o17) &
       )

    dfdy_nuc(jo17,jf18) = ( &
      0.0d0 &
       )

    dfdy_nuc(jf17,jp) = ( &
      screened_rates(k_p_o16__f17)*Y(jo16)*dens &
       )

    dfdy_nuc(jf17,jhe4) = ( &
      screened_rates(k_he4_o14__p_f17)*Y(jo14)*dens &
       )

    dfdy_nuc(jf17,jc12) = ( &
      0.0d0 &
       )

    dfdy_nuc(jf17,jc13) = ( &
      0.0d0 &
       )

    dfdy_nuc(jf17,jn13) = ( &
      0.0d0 &
       )

    dfdy_nuc(jf17,jn14) = ( &
      0.0d0 &
       )

    dfdy_nuc(jf17,jn15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jf17,jo14) = ( &
      screened_rates(k_he4_o14__p_f17)*Y(jhe4)*dens &
       )

    dfdy_nuc(jf17,jo15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jf17,jo16) = ( &
      screened_rates(k_p_o16__f17)*Y(jp)*dens &
       )

    dfdy_nuc(jf17,jo17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jf17,jf17) = ( &
      -screened_rates(k_f17__o17) &
       )

    dfdy_nuc(jf17,jf18) = ( &
      0.0d0 &
       )

    dfdy_nuc(jf18,jp) = ( &
      -screened_rates(k_p_f18__he4_o15)*Y(jf18)*dens + screened_rates(k_p_o17__f18)*Y(jo17)* &
      dens &
       )

    dfdy_nuc(jf18,jhe4) = ( &
      screened_rates(k_he4_n14__f18)*Y(jn14)*dens &
       )

    dfdy_nuc(jf18,jc12) = ( &
      0.0d0 &
       )

    dfdy_nuc(jf18,jc13) = ( &
      0.0d0 &
       )

    dfdy_nuc(jf18,jn13) = ( &
      0.0d0 &
       )

    dfdy_nuc(jf18,jn14) = ( &
      screened_rates(k_he4_n14__f18)*Y(jhe4)*dens &
       )

    dfdy_nuc(jf18,jn15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jf18,jo14) = ( &
      0.0d0 &
       )

    dfdy_nuc(jf18,jo15) = ( &
      0.0d0 &
       )

    dfdy_nuc(jf18,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jf18,jo17) = ( &
      screened_rates(k_p_o17__f18)*Y(jp)*dens &
       )

    dfdy_nuc(jf18,jf17) = ( &
      0.0d0 &
       )

    dfdy_nuc(jf18,jf18) = ( &
      -screened_rates(k_p_f18__he4_o15)*Y(jp)*dens &
       )

    
  end subroutine jac_nuc

end module actual_rhs_module
