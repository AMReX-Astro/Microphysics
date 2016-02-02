module actual_rhs_module

  implicit none

contains

  subroutine actual_rhs(neq, time, y, ydot, rpar)

    use bl_types
    use bl_constants_module
    use network
    use vode_data
    use rpar_indices
    use screening_module, only: screenz
    use actual_burner_data
    use actual_burner_module, only: ener_gener_rate

    implicit none

    integer          :: neq
    double precision :: time
    double precision :: y(neq), ydot(neq)
    double precision :: rpar(n_rpar_comps)

    double precision :: temp, T9, T9a, dT9dt, dT9adt

    double precision :: scratch, dscratchdt
    double precision :: rate, dratedt, sc1212, dsc1212dt, xc12tmp

    double precision :: dens

    double precision :: a, b, dadt, dbdt

    temp    = y(net_itemp)
    dens    = rpar(irp_dens)

    ! call the screening routine
    call screenz(temp,dens,6.0d0,6.0d0,12.0d0,12.0d0,y(1:nspec),aion,zion,nspec,sc1212,dsc1212dt)

    ! compute some often used temperature constants
    T9     = temp/1.d9
    dT9dt  = ONE/1.d9
    T9a    = T9/(1.0d0 + 0.0396d0*T9)
    dT9adt = (T9a / T9 - (T9a / (1.0d0 + 0.0396d0*T9)) * 0.0396d0) * dT9dt

    ! compute the CF88 rate
    scratch    = T9a**THIRD
    dscratchdt = THIRD * T9a**(-TWO3RD) * dT9adt

    a       = 4.27d26*T9a**(FIVE*SIXTH)*T9**(-1.5d0)
    dadt    = (FIVE * SIXTH) * (a/T9a) * dT9adt - 1.5d0 * (a/T9) * dT9dt

    b       = dexp(-84.165d0/scratch - 2.12d-3*T9*T9*T9)
    dbdt    = (84.165d0 * dscratchdt/ scratch**TWO - THREE * 2.12d-3 * T9 * T9 * dT9dt) * b

    rate    = a *  b
    dratedt = dadt * b + a * dbdt

    ! The change in number density of C12 is
    ! d(n12)/dt = - 2 * 1/2 (n12)**2 <sigma v>
    !
    ! where <sigma v> is the average of the relative velocity times the cross
    ! section for the reaction, and the factor accounting for the total number
    ! of particle pairs has a 1/2 because we are considering a reaction involving 
    ! identical particles (see Clayton p. 293).  Finally, the -2 means that for
    ! each reaction, we lose 2 carbon nuclei.
    !
    ! The corresponding Mg24 change is
    ! d(n24)/dt = + 1/2 (n12)**2 <sigma v>
    !
    ! note that no factor of 2 appears here, because we create only 1 Mg nuclei.
    !
    ! Switching over to mass fractions, using n = rho X N_A/A, where N_A is
    ! Avagadro's number, and A is the mass number of the nucleon, we get
    !
    ! d(X12)/dt = -2 *1/2 (X12)**2 rho N_A <sigma v> / A12
    !
    ! d(X24)/dt = + 1/2 (X12)**2 rho N_A <sigma v> (A24/A12**2)
    !
    ! these are equal and opposite.
    !
    ! The quantity [N_A <sigma v>] is what is tabulated in Caughlin and Fowler.

    ! we will always refer to the species by integer indices that come from
    ! the network module -- this makes things robust to a shuffling of the 
    ! species ordering

    xc12tmp = max(y(ic12) * aion(ic12), ZERO)
    ydot(ic12)  = -TWELFTH * dens * sc1212 * rate * xc12tmp**2
    ydot(io16)  = ZERO
    ydot(img24) = -ydot(ic12)

    ! These get sent to the Jacobian

    rpar(irp_rate)      = rate
    rpar(irp_dratedt)   = dratedt
    rpar(irp_sc1212)    = sc1212
    rpar(irp_dsc1212dt) = dsc1212dt
    rpar(irp_xc12tmp)   = xc12tmp

    ! Convert back to molar form

    ydot(1:nspec) = ydot(1:nspec) / aion

    call ener_gener_rate(ydot(1:nspec), ydot(net_ienuc))

    call temperature_rhs(neq, y, ydot, rpar)

  end subroutine actual_rhs



  subroutine actual_jac(neq, time, y, pd, rpar)

    use network
    use rpar_indices
    use bl_constants_module, only: ZERO, SIXTH, TWELFTH
    use vode_data
    use actual_burner_data
    use actual_burner_module, only: ener_gener_rate

    implicit none

    integer          :: neq
    double precision :: time
    double precision :: y(neq), pd(neq, neq)
    double precision :: rpar(n_rpar_comps)

    double precision :: dens
    double precision :: rate, dratedt, scorr, dscorrdt, xc12tmp

    integer :: j

    ! Get data from the rpar array

    dens     = rpar(irp_dens)

    rate     = rpar(irp_rate)     
    dratedt  = rpar(irp_dratedt)  
    scorr    = rpar(irp_sc1212)   
    dscorrdt = rpar(irp_dsc1212dt)
    xc12tmp  = rpar(irp_xc12tmp)    

    ! initialize
    pd(:,:)  = ZERO

    ! carbon jacobian elements
    pd(ic12, ic12)  = -SIXTH * dens * scorr * rate * xc12tmp

    ! add the temperature derivatives: df(y_i) / dT
    pd(ic12, net_itemp)  = -TWELFTH * ( dens * rate * xc12tmp**2 * dscorrdt + &
                            dens * scorr * xc12tmp**2 * dratedt )

    ! Energy generation rate Jacobian elements with respect to species

    do j = 1, nspec
       call ener_gener_rate(pd(1:nspec,j), pd(net_ienuc,j))
    enddo

    ! Jacobian elements with respect to temperature

    call ener_gener_rate(pd(1:nspec,net_itemp), pd(net_ienuc,net_itemp))

    call temperature_jac(neq, y, pd, rpar)

  end subroutine actual_jac

end module actual_rhs_module
