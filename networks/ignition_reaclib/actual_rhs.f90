module actual_rhs_module

  use bl_types
  use bl_constants_module
  use physical_constants, only: N_AVO
  use network
  use net_rates, only: nreact, screen_reaclib, rate_evaluate
  use screening_module, only: plasma_state, fill_plasma_state
  use burn_type_module

  implicit none

contains

  subroutine actual_rhs(state)
    
    !$acc routine seq

    use extern_probin_module, only: do_constant_volume_burn
    use burn_type_module, only: num_rate_groups, net_itemp, net_ienuc
    use table_rates, only: num_tables

    implicit none

    type(burn_t) :: state
    type(plasma_state) :: pstate
    double precision, dimension(nspec) :: Y
    double precision :: reactvec(num_rate_groups+2)
    double precision :: dqweak(num_tables)
    double precision :: epart(num_tables)
    integer :: i
    double precision :: dens, temp, rhoy, ye

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
    ! Calculate rates
    call fill_plasma_state(pstate, temp, dens, Y(1:nspec))
    do i = 1, nreact
       call rate_evaluate(pstate, rhoy, temp, i, reactvec)
       state%rates(:,i) = reactvec(1:4)
    end do

    state%ydot(jn) = ( &
       - Y(jn) * state%rates(i_scor, k_n_p) * state%rates(i_rate, k_n_p) &
       + 5.00000000000000d-01 * dens * Y(jc12)**2 * state%rates(i_scor, k_c12_c12n_mg23) * state%rates(i_rate, k_c12_c12n_mg23) &
       )

    state%ydot(jp) = ( &
       + 5.00000000000000d-01 * dens * Y(jc12)**2 * state%rates(i_scor, k_c12_c12p_na23) * state%rates(i_rate, k_c12_c12p_na23) &
       + Y(jn) * state%rates(i_scor, k_n_p) * state%rates(i_rate, k_n_p) &
       )

    state%ydot(jhe4) = ( &
       + 5.00000000000000d-01 * dens * Y(jc12)**2 * state%rates(i_scor, k_c12_c12a_ne20) * state%rates(i_rate, k_c12_c12a_ne20) &
       )

    state%ydot(jc12) = ( &
       - 2 * 5.00000000000000d-01 * dens * Y(jc12)**2 * state%rates(i_scor, k_c12_c12a_ne20) * state%rates(i_rate, k_c12_c12a_ne20) &
       - 2 * 5.00000000000000d-01 * dens * Y(jc12)**2 * state%rates(i_scor, k_c12_c12n_mg23) * state%rates(i_rate, k_c12_c12n_mg23) &
       - 2 * 5.00000000000000d-01 * dens * Y(jc12)**2 * state%rates(i_scor, k_c12_c12p_na23) * state%rates(i_rate, k_c12_c12p_na23) &
       )

    state%ydot(jne20) = ( &
       + 5.00000000000000d-01 * dens * Y(jc12)**2 * state%rates(i_scor, k_c12_c12a_ne20) * state%rates(i_rate, k_c12_c12a_ne20) &
       )

    state%ydot(jna23) = ( &
       + 5.00000000000000d-01 * dens * Y(jc12)**2 * state%rates(i_scor, k_c12_c12p_na23) * state%rates(i_rate, k_c12_c12p_na23) &
       )

    state%ydot(jmg23) = ( &
       + 5.00000000000000d-01 * dens * Y(jc12)**2 * state%rates(i_scor, k_c12_c12n_mg23) * state%rates(i_rate, k_c12_c12n_mg23) &
       )

    ! Let temperature be a constant
    state%ydot(net_itemp) = 0.0d0
    
    state%ydot(net_ienuc) = 0.0d0
    ! ion binding energy contributions
    do i = 1, nspec
       state%ydot(net_ienuc) = state%ydot(net_ienuc) + N_AVO * bion(i) * state%ydot(i)
    end do
    
    ! weak Q-value modification dqweak (density and temperature dependent)
    
    ! weak particle energy generation rates from gamma heating and neutrino loss
    ! (does not include plasma neutrino losses)

    ! write(*,*) '______________________________'
    ! do i = 1, nspec+1
    !    write(*,*) 'state%ydot(',i,'): ',state%ydot(i)
    ! end do
  end subroutine actual_rhs

  subroutine actual_jac(state)

    !$acc routine seq

    use burn_type_module, only: num_rate_groups, net_itemp, net_ienuc
    use table_rates, only: num_tables
    
    implicit none
    
    type(burn_t) :: state
    type(plasma_state) :: pstate
    double precision :: reactvec(num_rate_groups+2)
    double precision :: Y(nspec)
    double precision :: dens, temp, ye, rhoy
    integer :: i, j

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

    if (.not. state%have_rates) then
       ! Calculate rates
       call fill_plasma_state(pstate, temp, dens, Y(1:nspec))
       do i = 1, nreact
          call rate_evaluate(pstate, rhoy, temp, i, reactvec)
          state%rates(:,i) = reactvec(1:4)
       end do
    end if
    
    ! state%jac(j, i) = d(YDOT(j))/dY(i)
    state%jac(jn,jn) = ( &
         -    state%rates(i_scor, k_n_p) * state%rates(i_rate, k_n_p) &
         )

    state%jac(jn,jp) = ( &
         + 0.0d0 &
         )

    state%jac(jn,jhe4) = ( &
         + 0.0d0 &
         )

    state%jac(jn,jc12) = ( &
         + 5.00000000000000d-01 * dens * 2*Y(jc12) * state%rates(i_scor, k_c12_c12n_mg23) * state%rates(i_rate, k_c12_c12n_mg23) &
         )

    state%jac(jn,jne20) = ( &
         + 0.0d0 &
         )

    state%jac(jn,jna23) = ( &
         + 0.0d0 &
         )

    state%jac(jn,jmg23) = ( &
         + 0.0d0 &
         )

    state%jac(jp,jn) = ( &
         +    state%rates(i_scor, k_n_p) * state%rates(i_rate, k_n_p) &
         )

    state%jac(jp,jp) = ( &
         + 0.0d0 &
         )

    state%jac(jp,jhe4) = ( &
         + 0.0d0 &
         )

    state%jac(jp,jc12) = ( &
         + 5.00000000000000d-01 * dens * 2*Y(jc12) * state%rates(i_scor, k_c12_c12p_na23) * state%rates(i_rate, k_c12_c12p_na23) &
         )

    state%jac(jp,jne20) = ( &
         + 0.0d0 &
         )

    state%jac(jp,jna23) = ( &
         + 0.0d0 &
         )

    state%jac(jp,jmg23) = ( &
         + 0.0d0 &
         )

    state%jac(jhe4,jn) = ( &
         + 0.0d0 &
         )

    state%jac(jhe4,jp) = ( &
         + 0.0d0 &
         )

    state%jac(jhe4,jhe4) = ( &
         + 0.0d0 &
         )

    state%jac(jhe4,jc12) = ( &
         + 5.00000000000000d-01 * dens * 2*Y(jc12) * state%rates(i_scor, k_c12_c12a_ne20) * state%rates(i_rate, k_c12_c12a_ne20) &
         )

    state%jac(jhe4,jne20) = ( &
         + 0.0d0 &
         )

    state%jac(jhe4,jna23) = ( &
         + 0.0d0 &
         )

    state%jac(jhe4,jmg23) = ( &
         + 0.0d0 &
         )

    state%jac(jc12,jn) = ( &
         + 0.0d0 &
         )

    state%jac(jc12,jp) = ( &
         + 0.0d0 &
         )

    state%jac(jc12,jhe4) = ( &
         + 0.0d0 &
         )

    state%jac(jc12,jc12) = ( &
         - 2 * 5.00000000000000d-01 * dens * 2*Y(jc12) * state%rates(i_scor, k_c12_c12a_ne20) * state%rates(i_rate, k_c12_c12a_ne20) &
         - 2 * 5.00000000000000d-01 * dens * 2*Y(jc12) * state%rates(i_scor, k_c12_c12n_mg23) * state%rates(i_rate, k_c12_c12n_mg23) &
         - 2 * 5.00000000000000d-01 * dens * 2*Y(jc12) * state%rates(i_scor, k_c12_c12p_na23) * state%rates(i_rate, k_c12_c12p_na23) &
         )

    state%jac(jc12,jne20) = ( &
         + 0.0d0 &
         )

    state%jac(jc12,jna23) = ( &
         + 0.0d0 &
         )

    state%jac(jc12,jmg23) = ( &
         + 0.0d0 &
         )

    state%jac(jne20,jn) = ( &
         + 0.0d0 &
         )

    state%jac(jne20,jp) = ( &
         + 0.0d0 &
         )

    state%jac(jne20,jhe4) = ( &
         + 0.0d0 &
         )

    state%jac(jne20,jc12) = ( &
         + 5.00000000000000d-01 * dens * 2*Y(jc12) * state%rates(i_scor, k_c12_c12a_ne20) * state%rates(i_rate, k_c12_c12a_ne20) &
         )

    state%jac(jne20,jne20) = ( &
         + 0.0d0 &
         )

    state%jac(jne20,jna23) = ( &
         + 0.0d0 &
         )

    state%jac(jne20,jmg23) = ( &
         + 0.0d0 &
         )

    state%jac(jna23,jn) = ( &
         + 0.0d0 &
         )

    state%jac(jna23,jp) = ( &
         + 0.0d0 &
         )

    state%jac(jna23,jhe4) = ( &
         + 0.0d0 &
         )

    state%jac(jna23,jc12) = ( &
         + 5.00000000000000d-01 * dens * 2*Y(jc12) * state%rates(i_scor, k_c12_c12p_na23) * state%rates(i_rate, k_c12_c12p_na23) &
         )

    state%jac(jna23,jne20) = ( &
         + 0.0d0 &
         )

    state%jac(jna23,jna23) = ( &
         + 0.0d0 &
         )

    state%jac(jna23,jmg23) = ( &
         + 0.0d0 &
         )

    state%jac(jmg23,jn) = ( &
         + 0.0d0 &
         )

    state%jac(jmg23,jp) = ( &
         + 0.0d0 &
         )

    state%jac(jmg23,jhe4) = ( &
         + 0.0d0 &
         )

    state%jac(jmg23,jc12) = ( &
         + 5.00000000000000d-01 * dens * 2*Y(jc12) * state%rates(i_scor, k_c12_c12n_mg23) * state%rates(i_rate, k_c12_c12n_mg23) &
         )

    state%jac(jmg23,jne20) = ( &
         + 0.0d0 &
         )

    state%jac(jmg23,jna23) = ( &
         + 0.0d0 &
         )

    state%jac(jmg23,jmg23) = ( &
         + 0.0d0 &
         )

  end subroutine actual_jac

end module actual_rhs_module
