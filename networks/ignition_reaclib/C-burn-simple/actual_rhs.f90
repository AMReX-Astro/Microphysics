module actual_rhs_module

  use bl_types
  use bl_constants_module
  use physical_constants, only: N_AVO
  use actual_burner_data, only: c_light
  use network
  use actual_network_data, only: nrat_reaclib, nrat_tabular
  use net_rates, only: screen_reaclib, reaclib_evaluate
  use table_rates
  use screening_module, only: plasma_state, fill_plasma_state
  use burn_type_module

  implicit none

contains

  subroutine actual_rhs(state)
    
    !$acc routine seq

    use extern_probin_module, only: do_constant_volume_burn
    use burn_type_module, only: num_rate_groups, net_itemp, net_ienuc

    implicit none

    type(burn_t) :: state
    type(plasma_state) :: pstate
    double precision, dimension(nspec) :: Y
    double precision :: reactvec(num_rate_groups+2)
    double precision :: dqweak(nrat_tabular)
    double precision :: epart(nrat_tabular)
    integer :: i, j
    double precision :: dens, temp, rhoy, ye

    double precision :: scratch_0
    double precision :: scratch_1
    double precision :: scratch_2
    double precision :: scratch_3
    double precision :: scratch_4
    double precision :: scratch_5
    double precision :: scratch_6
    double precision :: scratch_7
    
    ! ! Notify
    ! write(*,*) '(Executing subroutine actual_rhs)'

    ! ! Print out the composition
    ! do i = 1, nspec
    !    write(*,*) 'state%xn(', i, '): ', state%xn(i)
    ! end do

    ! Set molar abundances and enforce them to be positive
    do i = 1, nspec
       Y(i) = max(0.0d0, state%xn(i)/aion(i))
    end do

    dens = state%rho
    temp = state%T
    ye   = state%y_e
    rhoy = dens*ye

    ! For now this doesn't have dqweak or epart, need to add those
    ! Calculate Reaclib rates
    call fill_plasma_state(pstate, temp, dens, Y(1:nspec))
    do i = 1, nrat_reaclib
       call reaclib_evaluate(pstate, temp, i, reactvec)
       state%rates(:,i) = reactvec(1:4)
    end do

    ! Include only if there are tabular rates
    do i = 1, nrat_tabular
       call tabular_evaluate(table_meta(i), rhoy, temp, reactvec)
       j = i + nrat_reaclib
       state%rates(:,j) = reactvec(1:4)
       dqweak(i) = reactvec(5)
       epart(i)  = reactvec(6)
    end do

    scratch_0 = state%rates(i_rate, k_n_p)*state%rates(i_scor, k_n_p)*Y(jn)
    scratch_1 = Y(jc12)**2*dens
    scratch_2 = state%rates(i_rate, k_c12_c12n_mg23)*state%rates(i_scor, k_c12_c12n_mg23)*scratch_1
    scratch_3 = 0.5d0*scratch_2
    scratch_4 = state%rates(i_rate, k_c12_c12p_na23)*state%rates(i_scor, k_c12_c12p_na23)*scratch_1
    scratch_5 = 0.5d0*scratch_4
    scratch_6 = state%rates(i_rate, k_c12_c12a_ne20)*state%rates(i_scor, k_c12_c12a_ne20)*scratch_1
    scratch_7 = 0.5d0*scratch_6

    state%ydot(jn) = ( &
      -scratch_0 + scratch_3 &
       )

    state%ydot(jp) = ( &
      scratch_0 + scratch_5 &
       )

    state%ydot(jhe4) = ( &
      scratch_7 &
       )

    state%ydot(jc12) = ( &
      -scratch_2 - scratch_4 - scratch_6 &
       )

    state%ydot(jne20) = ( &
      scratch_7 &
       )

    state%ydot(jna23) = ( &
      scratch_5 &
       )

    state%ydot(jmg23) = ( &
      scratch_3 &
       )

    
    state%ydot(net_ienuc) = 0.0d0
    ! ion binding energy contributions
    do i = 1, nspec
       state%ydot(net_ienuc) = state%ydot(net_ienuc) - N_AVO * mion(i) * c_light * c_light * state%ydot(i)
    end do
    
    ! weak Q-value modification dqweak (density and temperature dependent)
    
    ! weak particle energy generation rates from gamma heating and neutrino loss
    ! (does not include plasma neutrino losses)

    ! Let temperature be a constant
    state%ydot(net_itemp) = 0.0d0
    
    ! write(*,*) '______________________________'
    ! do i = 1, nspec+2
    !    write(*,*) 'state%ydot(',i,'): ',state%ydot(i)
    ! end do
  end subroutine actual_rhs

  subroutine actual_jac(state)

    !$acc routine seq

    use burn_type_module, only: num_rate_groups, net_itemp, net_ienuc
    
    implicit none
    
    type(burn_t) :: state
    type(plasma_state) :: pstate
    double precision :: reactvec(num_rate_groups+2)
    double precision :: Y(nspec)
    double precision :: dens, temp, ye, rhoy
    integer :: i, j

    double precision :: scratch_0
    double precision :: scratch_1
    double precision :: scratch_2
    double precision :: scratch_3
    double precision :: scratch_4
    double precision :: scratch_5
    double precision :: scratch_6
    double precision :: scratch_7
    double precision :: scratch_8

    ! ! Notify
    ! write(*,*) '(Executing subroutine actual_jac)'

    ! ! Print out the composition
    ! do i = 1, nspec
    !    write(*,*) 'state%xn(', i, '): ', state%xn(i)
    ! end do
    
    dens = state%rho
    temp = state%T
    ye   = state%y_e
    rhoy = dens*ye

    ! Set molar abundances and enforce them to be positive
    do i = 1, nspec
       Y(i) = max(0.0d0, state%xn(i)/aion(i))
    end do

    state%jac(:,:) = ZERO
    
    if (.not. state%have_rates) then
       ! Calculate Reaclib rates
       call fill_plasma_state(pstate, temp, dens, Y(1:nspec))
       do i = 1, nrat_reaclib
          call reaclib_evaluate(pstate, temp, i, reactvec)
          state%rates(:,i) = reactvec(1:4)
       end do

       ! Include only if there are tabular rates
       do i = 1, nrat_tabular
          call tabular_evaluate(table_meta(i), rhoy, temp, reactvec)
          j = i + nrat_reaclib
          state%rates(:,j) = reactvec(1:4)
       end do
    end if

    scratch_0 = state%rates(i_rate, k_n_p)*state%rates(i_scor, k_n_p)
    scratch_1 = state%rates(i_rate, k_c12_c12n_mg23)*state%rates(i_scor, k_c12_c12n_mg23)
    scratch_2 = 1.0d0*Y(jc12)*dens
    scratch_3 = scratch_1*scratch_2
    scratch_4 = state%rates(i_rate, k_c12_c12p_na23)*state%rates(i_scor, k_c12_c12p_na23)
    scratch_5 = scratch_2*scratch_4
    scratch_6 = state%rates(i_rate, k_c12_c12a_ne20)*state%rates(i_scor, k_c12_c12a_ne20)*Y(jc12)* &
      dens
    scratch_7 = 1.0d0*scratch_6
    scratch_8 = 2.0d0*Y(jc12)*dens
    
    ! state%jac(j, i) = d(YDOT(j))/dY(i)
    
      state%jac(jn,jn) = ( &
        -scratch_0 &
         )

      state%jac(jn,jp) = ( &
        0.0d0 &
         )

      state%jac(jn,jhe4) = ( &
        0.0d0 &
         )

      state%jac(jn,jc12) = ( &
        scratch_3 &
         )

      state%jac(jn,jne20) = ( &
        0.0d0 &
         )

      state%jac(jn,jna23) = ( &
        0.0d0 &
         )

      state%jac(jn,jmg23) = ( &
        0.0d0 &
         )

      state%jac(jp,jn) = ( &
        scratch_0 &
         )

      state%jac(jp,jp) = ( &
        0.0d0 &
         )

      state%jac(jp,jhe4) = ( &
        0.0d0 &
         )

      state%jac(jp,jc12) = ( &
        scratch_5 &
         )

      state%jac(jp,jne20) = ( &
        0.0d0 &
         )

      state%jac(jp,jna23) = ( &
        0.0d0 &
         )

      state%jac(jp,jmg23) = ( &
        0.0d0 &
         )

      state%jac(jhe4,jn) = ( &
        0.0d0 &
         )

      state%jac(jhe4,jp) = ( &
        0.0d0 &
         )

      state%jac(jhe4,jhe4) = ( &
        0.0d0 &
         )

      state%jac(jhe4,jc12) = ( &
        scratch_7 &
         )

      state%jac(jhe4,jne20) = ( &
        0.0d0 &
         )

      state%jac(jhe4,jna23) = ( &
        0.0d0 &
         )

      state%jac(jhe4,jmg23) = ( &
        0.0d0 &
         )

      state%jac(jc12,jn) = ( &
        0.0d0 &
         )

      state%jac(jc12,jp) = ( &
        0.0d0 &
         )

      state%jac(jc12,jhe4) = ( &
        0.0d0 &
         )

      state%jac(jc12,jc12) = ( &
        -scratch_1*scratch_8 - scratch_4*scratch_8 - 2.0d0*scratch_6 &
         )

      state%jac(jc12,jne20) = ( &
        0.0d0 &
         )

      state%jac(jc12,jna23) = ( &
        0.0d0 &
         )

      state%jac(jc12,jmg23) = ( &
        0.0d0 &
         )

      state%jac(jne20,jn) = ( &
        0.0d0 &
         )

      state%jac(jne20,jp) = ( &
        0.0d0 &
         )

      state%jac(jne20,jhe4) = ( &
        0.0d0 &
         )

      state%jac(jne20,jc12) = ( &
        scratch_7 &
         )

      state%jac(jne20,jne20) = ( &
        0.0d0 &
         )

      state%jac(jne20,jna23) = ( &
        0.0d0 &
         )

      state%jac(jne20,jmg23) = ( &
        0.0d0 &
         )

      state%jac(jna23,jn) = ( &
        0.0d0 &
         )

      state%jac(jna23,jp) = ( &
        0.0d0 &
         )

      state%jac(jna23,jhe4) = ( &
        0.0d0 &
         )

      state%jac(jna23,jc12) = ( &
        scratch_5 &
         )

      state%jac(jna23,jne20) = ( &
        0.0d0 &
         )

      state%jac(jna23,jna23) = ( &
        0.0d0 &
         )

      state%jac(jna23,jmg23) = ( &
        0.0d0 &
         )

      state%jac(jmg23,jn) = ( &
        0.0d0 &
         )

      state%jac(jmg23,jp) = ( &
        0.0d0 &
         )

      state%jac(jmg23,jhe4) = ( &
        0.0d0 &
         )

      state%jac(jmg23,jc12) = ( &
        scratch_3 &
         )

      state%jac(jmg23,jne20) = ( &
        0.0d0 &
         )

      state%jac(jmg23,jna23) = ( &
        0.0d0 &
         )

      state%jac(jmg23,jmg23) = ( &
        0.0d0 &
         )


  end subroutine actual_jac

end module actual_rhs_module
