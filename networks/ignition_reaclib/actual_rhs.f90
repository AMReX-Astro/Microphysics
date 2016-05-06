module actual_rhs_module

  use bl_types
  use bl_constants_module
  use physical_constants, only: N_AVO
  use network
  use net_rates, only: nreact, screen_reaclib, rate_evaluate
  use screening_module, only: plasma_state, fill_plasma_state

  implicit none

contains

  subroutine actual_rhs(state)
    
    !$acc routine seq

    use extern_probin_module, only: do_constant_volume_burn

    implicit none

    type(burn_t) :: state
    type(plasma_pstate) :: pstate
    double precision, dimension(nspec) :: Y
    integer :: i
    double precision :: dens, temp, rhoy

    ! Set molar abundances and enforce them to be positive
    do i = 1, nspec
       Y(i) = max(0.0d0, state%xn(i)/aion(i))
    end do

    dens = state%rho
    temp = state%T
    ye   = state%y_e
    rhoy = dens*ye
    
    ! Calculate rates
    call fill_plasma_pstate(pstate, temp, dens, Y(1:nspec))
    do i = 1, nreact
       call rate_evaluate(pstate, rhoy, temp, i, state%rates(:,i))
    end do
    state%rates(:,:) = state%rates(1:4,:)

    write(*,*) "RHS Time: ", T
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

    state%ydot(jenuc) = 0.0d0
    ! ion binding energy contributions
    do i = 1, nspec
       state%ydot(jenuc) = state%ydot(jenuc) + N_AVO * ebind(i) * state%ydot(i)
    end do

    ! weak Q-value modification dqweak (density and temperature dependent)
    
    ! weak particle energy generation rates from gamma heating and neutrino loss
    ! (does not include plasma neutrino losses)
    
    write(*,*) '______________________________'
    do i = 1, nspec+1
       write(*,*) 'state%ydot(',i,'): ',state%ydot(i)
    end do
  end subroutine actual_rhs

  subroutine actual_jac(state)

    !$acc routine seq

    implicit none
    
    type(burn_t) :: state
    type(plasma_pstate) :: pstate
    double precision :: Y(nspec)
    double precision :: dens, temp, ye, rhoy
    integer :: i, j

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
       call fill_plasma_pstate(pstate, temp, dens, Y(1:nspec))
       do i = 1, nreact
          call rate_evaluate(pstate, rhoy, temp, i, state%rates(:,i))
       end do
    end if
    
    ! DJAC(j, i) = d(YDOT(j))/dY(i)
    DJAC(jn,jn) = ( &
         -    state%rates(i_scor, k_n_p) * state%rates(i_rate, k_n_p) &
         )

    DJAC(jn,jp) = ( &
         + 0.0d0 &
         )

    DJAC(jn,jhe4) = ( &
         + 0.0d0 &
         )

    DJAC(jn,jc12) = ( &
         + 5.00000000000000d-01 * dens * 2*Y(jc12) * state%rates(i_scor, k_c12_c12n_mg23) * state%rates(i_rate, k_c12_c12n_mg23) &
         )

    DJAC(jn,jne20) = ( &
         + 0.0d0 &
         )

    DJAC(jn,jna23) = ( &
         + 0.0d0 &
         )

    DJAC(jn,jmg23) = ( &
         + 0.0d0 &
         )

    DJAC(jp,jn) = ( &
         +    state%rates(i_scor, k_n_p) * state%rates(i_rate, k_n_p) &
         )

    DJAC(jp,jp) = ( &
         + 0.0d0 &
         )

    DJAC(jp,jhe4) = ( &
         + 0.0d0 &
         )

    DJAC(jp,jc12) = ( &
         + 5.00000000000000d-01 * dens * 2*Y(jc12) * state%rates(i_scor, k_c12_c12p_na23) * state%rates(i_rate, k_c12_c12p_na23) &
         )

    DJAC(jp,jne20) = ( &
         + 0.0d0 &
         )

    DJAC(jp,jna23) = ( &
         + 0.0d0 &
         )

    DJAC(jp,jmg23) = ( &
         + 0.0d0 &
         )

    DJAC(jhe4,jn) = ( &
         + 0.0d0 &
         )

    DJAC(jhe4,jp) = ( &
         + 0.0d0 &
         )

    DJAC(jhe4,jhe4) = ( &
         + 0.0d0 &
         )

    DJAC(jhe4,jc12) = ( &
         + 5.00000000000000d-01 * dens * 2*Y(jc12) * state%rates(i_scor, k_c12_c12a_ne20) * state%rates(i_rate, k_c12_c12a_ne20) &
         )

    DJAC(jhe4,jne20) = ( &
         + 0.0d0 &
         )

    DJAC(jhe4,jna23) = ( &
         + 0.0d0 &
         )

    DJAC(jhe4,jmg23) = ( &
         + 0.0d0 &
         )

    DJAC(jc12,jn) = ( &
         + 0.0d0 &
         )

    DJAC(jc12,jp) = ( &
         + 0.0d0 &
         )

    DJAC(jc12,jhe4) = ( &
         + 0.0d0 &
         )

    DJAC(jc12,jc12) = ( &
         - 2 * 5.00000000000000d-01 * dens * 2*Y(jc12) * state%rates(i_scor, k_c12_c12a_ne20) * state%rates(i_rate, k_c12_c12a_ne20) &
         - 2 * 5.00000000000000d-01 * dens * 2*Y(jc12) * state%rates(i_scor, k_c12_c12n_mg23) * state%rates(i_rate, k_c12_c12n_mg23) &
         - 2 * 5.00000000000000d-01 * dens * 2*Y(jc12) * state%rates(i_scor, k_c12_c12p_na23) * state%rates(i_rate, k_c12_c12p_na23) &
         )

    DJAC(jc12,jne20) = ( &
         + 0.0d0 &
         )

    DJAC(jc12,jna23) = ( &
         + 0.0d0 &
         )

    DJAC(jc12,jmg23) = ( &
         + 0.0d0 &
         )

    DJAC(jne20,jn) = ( &
         + 0.0d0 &
         )

    DJAC(jne20,jp) = ( &
         + 0.0d0 &
         )

    DJAC(jne20,jhe4) = ( &
         + 0.0d0 &
         )

    DJAC(jne20,jc12) = ( &
         + 5.00000000000000d-01 * dens * 2*Y(jc12) * state%rates(i_scor, k_c12_c12a_ne20) * state%rates(i_rate, k_c12_c12a_ne20) &
         )

    DJAC(jne20,jne20) = ( &
         + 0.0d0 &
         )

    DJAC(jne20,jna23) = ( &
         + 0.0d0 &
         )

    DJAC(jne20,jmg23) = ( &
         + 0.0d0 &
         )

    DJAC(jna23,jn) = ( &
         + 0.0d0 &
         )

    DJAC(jna23,jp) = ( &
         + 0.0d0 &
         )

    DJAC(jna23,jhe4) = ( &
         + 0.0d0 &
         )

    DJAC(jna23,jc12) = ( &
         + 5.00000000000000d-01 * dens * 2*Y(jc12) * state%rates(i_scor, k_c12_c12p_na23) * state%rates(i_rate, k_c12_c12p_na23) &
         )

    DJAC(jna23,jne20) = ( &
         + 0.0d0 &
         )

    DJAC(jna23,jna23) = ( &
         + 0.0d0 &
         )

    DJAC(jna23,jmg23) = ( &
         + 0.0d0 &
         )

    DJAC(jmg23,jn) = ( &
         + 0.0d0 &
         )

    DJAC(jmg23,jp) = ( &
         + 0.0d0 &
         )

    DJAC(jmg23,jhe4) = ( &
         + 0.0d0 &
         )

    DJAC(jmg23,jc12) = ( &
         + 5.00000000000000d-01 * dens * 2*Y(jc12) * state%rates(i_scor, k_c12_c12n_mg23) * state%rates(i_rate, k_c12_c12n_mg23) &
         )

    DJAC(jmg23,jne20) = ( &
         + 0.0d0 &
         )

    DJAC(jmg23,jna23) = ( &
         + 0.0d0 &
         )

    DJAC(jmg23,jmg23) = ( &
         + 0.0d0 &
         )
  end subroutine actual_jac

end module actual_rhs_module
