module actual_rhs_module

  use amrex_constants_module
  use network
  use burn_type_module
  use temperature_integration_module, only: temperature_rhs, temperature_jac
  use rate_type_module

  implicit none

contains

  subroutine actual_rhs_init()

    implicit none

  end subroutine actual_rhs_init

#ifdef CUDA
  attributes(device) &
#endif
  subroutine actual_rhs(state)

    !$acc routine seq

    use managed_probin_module, only: cu_do_constant_volume_burn

    implicit none

    type (burn_t)    :: state
    type (rate_t)    :: rr

    double precision :: temp, T9, T9a, dT9dt, dT9adt

    double precision :: scratch, dscratchdt
    double precision :: rate, dratedt, sc1212, dsc1212dt, xc12tmp

    double precision :: dens

    double precision :: a, b, dadt, dbdt

    double precision :: y(nspec)


    ! We enforce that X(O16) remains constant, and that X(Mg24) always mirrors changes in X(C12).
    call update_unevolved_species(state)

    call evaluate_rates(state, rr)


    ! Now get the data from the state.

    temp = state % T
    dens = state % rho
    y    = state % xn * aion_inv

    rate      = rr % rates(1,1)
    dratedt   = rr % rates(2,1)
    sc1212    = rr % rates(3,1)
    dsc1212dt = rr % rates(4,1)

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

    xc12tmp = max(state % xn(ic12), ZERO)
    state % ydot(ic12)  = -TWELFTH * dens * sc1212 * rate * xc12tmp**2

    ! Convert back to molar form

    state % ydot(ic12) = state % ydot(ic12) * aion_inv(ic12)

    call ener_gener_rate(state % ydot(ic12), state % ydot(net_ienuc))

    ! Do the temperature equation explicitly here since
    ! the generic form doesn't work when nspec_evolve < nspec.

    if (state % self_heat) then

       if (cu_do_constant_volume_burn) then
          state % ydot(net_itemp) = state % ydot(net_ienuc) / state % cv

       else
          state % ydot(net_itemp) = state % ydot(net_ienuc) / state % cp

       endif
    endif

  end subroutine actual_rhs

#ifdef CUDA
  attributes(device) &
#endif
  subroutine actual_jac(state)

    !$acc routine seq

    use managed_probin_module, only: cu_do_constant_volume_burn

    implicit none

    type (burn_t)    :: state
    type (rate_t)    :: rr

    double precision :: dens
    double precision :: rate, dratedt, scorr, dscorrdt, xc12tmp

    double precision :: cvInv, cpInv

    call evaluate_rates(state, rr)

    ! Get data from the state

    dens     = state % rho

    rate     = rr % rates(1,1)
    dratedt  = rr % rates(2,1)
    scorr    = rr % rates(3,1)
    dscorrdt = rr % rates(4,1)
    xc12tmp  = max(state % xn(ic12), ZERO)

    ! initialize
    state % jac(:,:) = ZERO

    ! carbon jacobian elements
    state % jac(ic12, ic12)  = -SIXTH * dens * scorr * rate * xc12tmp

    ! add the temperature derivatives: df(y_i) / dT
    state % jac(ic12, net_itemp)  = -TWELFTH * ( dens * rate * xc12tmp**2 * dscorrdt + &
                                     dens * scorr * xc12tmp**2 * dratedt )

    ! Convert back to molar form
    ! Note that the factor of 1/A cancels in the (C12,C12) Jacobian element,
    ! so this conversion is necessarily only for the temperature derivative.

    state % jac(ic12,net_itemp) = state % jac(ic12,net_itemp) * aion_inv(ic12)

    ! Energy generation rate Jacobian elements with respect to species

    call ener_gener_rate(state % jac(ic12,ic12), state % jac(net_ienuc,ic12))

    ! Jacobian elements with respect to temperature

    call ener_gener_rate(state % jac(ic12,net_itemp), state % jac(net_ienuc,net_itemp))

    if (state % self_heat) then

       if (cu_do_constant_volume_burn) then

          cvInv = ONE / state % cv

          ! d(itemp)/d(yi)

          state % jac(net_itemp,ic12) = state % jac(net_ienuc,ic12) * cvInv

          ! d(itemp)/d(temp)

          state % jac(net_itemp, net_itemp) = state % jac(net_ienuc,net_itemp) * cvInv

       else

          cpInv = ONE / state % cp

          ! d(itemp)/d(yi)

          state % jac(net_itemp,ic12) = state % jac(net_ienuc,ic12) * cpInv

          ! d(itemp)/d(temp)

          state % jac(net_itemp,net_itemp) = state % jac(net_ienuc,net_itemp) * cpInv

       endif

    endif

  end subroutine actual_jac


#ifdef CUDA
  attributes(device) &
#endif
  subroutine evaluate_rates(state, rr)

    !$acc routine seq

    use screening_module, only: screenz

    implicit none

    type (burn_t), intent(in) :: state
    type (rate_t), intent(out) :: rr

    double precision :: temp, T9, T9a, dT9dt, dT9adt

    double precision :: scratch, dscratchdt
    double precision :: rate, dratedt, sc1212, dsc1212dt, xc12tmp

    double precision :: dens

    double precision :: a, b, dadt, dbdt

    double precision :: y(nspec)

    temp = state % T
    dens = state % rho
    y    = state % xn * aion_inv

    ! call the screening routine
    call screenz(temp,dens,6.0d0,6.0d0,12.0d0,12.0d0,y,sc1212,dsc1212dt)

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

    ! These get sent to the Jacobian

    rr % rates(1,:)  = rate
    rr % rates(2,:)  = dratedt
    rr % rates(3,:)  = sc1212
    rr % rates(4,:)  = dsc1212dt

    rr % T_eval = temp

  end subroutine evaluate_rates



  ! Computes the instantaneous energy generation rate
#ifdef CUDA
  attributes(device) &
#endif
  subroutine ener_gener_rate(dydt, enuc)

    use network

    !$acc routine seq

    implicit none

    double precision :: dydt(nspec_evolve), enuc

    ! This is basically e = m c**2

    ! Note that since we don't explicitly evolve Mg24
    ! in this network, we need to explicitly add its
    ! contribution in this routine. We can factor out
    ! the common factor of dydt(ic12), we just need to
    ! account for a factor of aion(ic12) / aion(img24)
    ! for the second term to make the expression work.

    enuc = dydt(ic12) * (mion(ic12) - mion(img24) / 2) * enuc_conv2

  end subroutine ener_gener_rate

#ifdef CUDA
  attributes(device) &
#endif  
  subroutine update_unevolved_species(state)

    !$acc routine seq

    implicit none

    type (burn_t)    :: state

    state % xn(iMg24) = ONE - state % xn(iC12) - state % xn(iO16)

  end subroutine update_unevolved_species

end module actual_rhs_module
