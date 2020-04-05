module vode_type_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module
  use cuvode_parameters_module, only : VODE_NEQS
  use cuvode_types_module, only : dvode_t
  use network, only : nspec, aion, aion_inv

  use vode_rpar_indices
  use sdc_type_module

  use extern_probin_module, only: renormalize_abundances

  implicit none

  ! this should be larger than any reasonable temperature we will encounter   
  real (kind=rt), parameter :: MAX_TEMP = 1.0e11_rt          

  public

contains

  subroutine clean_state(time, vode_state)

    use eos_type_module, only : eos_t, eos_input_re, eos_input_rt
    use eos_module, only : eos

    real(rt), intent(in) :: time
    type(dvode_t), intent(inout) :: vode_state

    real(rt) :: max_e, ke

    type (eos_t) :: eos_state

    !$gpu

    ! update rho, rho*u, etc.
    call fill_unevolved_variables(time, vode_state)

    ! Ensure that mass fractions always stay positive.
    vode_state % y(SFS:SFS+nspec-1) = &
         max(min(vode_state % y(SFS:SFS+nspec-1), vode_state % rpar(irp_SRHO)), &
             vode_state % rpar(irp_SRHO) * 1.e-200_rt)

    ! renormalize abundances as necessary
    if (renormalize_abundances) then
       call renormalize_species(time, vode_state)
    endif

#ifdef SDC_EVOLVE_ENERGY

    ! Ensure that internal energy never goes above the maximum limit
    ! provided by the EOS. Same for the internal energy implied by the
    ! total energy (which we get by subtracting kinetic energy).
    eos_state % rho = vode_state % rpar(irp_SRHO)
    eos_state % T = MAX_TEMP
    eos_state % xn = vode_state % y(SFS:SFS+nspec-1) / vode_state % rpar(irp_SRHO)

    call eos(eos_input_rt, eos_state)

    max_e = eos_state % e

    vode_state % y(SEINT) = min(vode_state % rpar(irp_SRHO) * max_e, vode_state % y(SEINT))

    ke = vode_state % y(SEDEN) - &
         HALF * sum(vode_state % rpar(irp_SMX:irp_SMZ)**2) / vode_state % rpar(irp_SRHO)

    vode_state % y(SEDEN) = min(vode_state % rpar(irp_SRHO) * max_e + ke, vode_state % y(SEDEN))

#endif

  end subroutine clean_state


  subroutine fill_unevolved_variables(time, vode_state)

    real(rt), intent(in) :: time
    type(dvode_t), intent(inout) :: vode_state

    !$gpu

#if defined(SDC_EVOLVE_ENERGY)

    ! we are always integrating from t = 0, so there is no offset
    ! time needed here.  The indexing of irp_ydot_a is based on
    ! the indices in sdc_type_module
    vode_state % rpar(irp_SRHO) = vode_state % rpar(irp_u_init-1+irp_SRHO) + &
         vode_state % rpar(irp_ydot_a-1+SRHO) * time

    vode_state % rpar(irp_SMX) = vode_state % rpar(irp_u_init-1+irp_SMX) + &
         vode_state % rpar(irp_ydot_a-1+SMX) * time
    vode_state % rpar(irp_SMY) = vode_state % rpar(irp_u_init-1+irp_SMY) + &
         vode_state % rpar(irp_ydot_a-1+SMY) * time
    vode_state % rpar(irp_SMZ) = vode_state % rpar(irp_u_init-1+irp_SMZ) + &
         vode_state % rpar(irp_ydot_a-1+SMZ) * time

#elif defined(SDC_EVOLVE_ENTHALPY)

    ! Keep density consistent with the partial densities.
    vode_state % rpar(irp_SRHO) = sum(vode_state % y(SFS:SFS - 1 + nspec))

#endif

  end subroutine fill_unevolved_variables

  subroutine renormalize_species(time, vode_state)

    real(rt), intent(in) :: time
    type(dvode_t), intent(inout) :: vode_state

    real(rt) :: nspec_sum

    !$gpu

    ! We only renormalize species when evolving energy because
    ! when we evolve enthalpy, we define the density as
    ! the sum of the partial densities rho*X for each species.
#ifdef SDC_EVOLVE_ENERGY

    ! update rho, rho*u, etc.
    call fill_unevolved_variables(time, vode_state)

    nspec_sum = sum(vode_state % y(SFS:SFS-1+nspec)) / vode_state % rpar(irp_SRHO)

    vode_state % y(SFS:SFS-1+nspec) = vode_state % y(SFS:SFS-1+nspec) / nspec_sum

#endif

  end subroutine renormalize_species


  ! Given a burn state, fill the rpar and integration state data.

  subroutine sdc_to_vode(sdc, vode_state)

    type (sdc_t) :: sdc
    type(dvode_t), intent(inout) :: vode_state

    !$gpu

    vode_state % y(:) = sdc % y(1:SVAR_EVOLVE)

    ! advective sources
    vode_state % rpar(irp_ydot_a:irp_ydot_a-1+SVAR) = sdc % ydot_a(:)

#if defined(SDC_EVOLVE_ENERGY)

    ! unevolved state variables
    vode_state % rpar(irp_SRHO) = sdc % y(SRHO)
    vode_state % rpar(irp_SMX:irp_SMZ) = sdc % y(SMX:SMZ)

    ! initial state for unevolved variables
    vode_state % rpar(irp_u_init-1+irp_SRHO) = sdc % y(SRHO)
    vode_state % rpar(irp_u_init-1+irp_SMX:irp_u_init-1+irp_SMZ) = sdc % y(SMX:SMZ)

    ! other parameters
    if (sdc % T_from_eden) then
       vode_state % rpar(irp_T_from_eden) = ONE
    else
       vode_state % rpar(irp_T_from_eden) = -ONE
    endif

#elif defined(SDC_EVOLVE_ENTHALPY)

    vode_state % rpar(irp_p0)   = sdc % p0
    vode_state % rpar(irp_SRHO) = sdc % rho

#endif

#ifdef NONAKA_PLOT

    ! bookkeeping information
    vode_state % rpar(irp_i) = sdc % i
    vode_state % rpar(irp_j) = sdc % j
    vode_state % rpar(irp_k) = sdc % k
    vode_state % rpar(irp_iter) = sdc % sdc_iter

#endif

  end subroutine sdc_to_vode

  subroutine vode_to_sdc(time, vode_state, sdc)

    real(rt), intent(in) :: time
    type (sdc_t) :: sdc
    type(dvode_t), intent(inout) :: vode_state

    !$gpu

    sdc % y(1:SVAR_EVOLVE) = vode_state % y(:)

    ! unevolved state variables
    call fill_unevolved_variables(time, vode_state)

#if defined(SDC_EVOLVE_ENERGY)

    sdc % y(SRHO) = vode_state % rpar(irp_SRHO)
    sdc % y(SMX:SMZ) = vode_state % rpar(irp_SMX:irp_SMZ)

#elif defined(SDC_EVOLVE_ENTHALPY)

    sdc % p0  = vode_state % rpar(irp_p0)
    sdc % rho = vode_state % rpar(irp_SRHO)

#endif

#ifdef NONAKA_PLOT

    sdc % i = vode_state % rpar(irp_i)
    sdc % j = vode_state % rpar(irp_j)
    sdc % k = vode_state % rpar(irp_k)
    sdc % sdc_iter = vode_state % rpar(irp_iter)

#endif

  end subroutine vode_to_sdc


  subroutine rhs_to_vode(time, burn_state, ydot_react, vode_state, ydot)

    use burn_type_module, only : burn_t, net_ienuc, neqs

    real(rt), intent(in) :: time
    type(dvode_t), intent(inout) :: vode_state
    type(burn_t), intent(in) :: burn_state
    real(rt), intent(in) :: ydot_react(neqs)
    real(rt), intent(out) :: ydot(SVAR_EVOLVE)

    !$gpu

    call fill_unevolved_variables(time, vode_state)

    ! ydot_react has just the contribution to the RHS from the
    ! reaction network.  Note that these are in terms of dY/dt

    ! start with the contribution from the non-reacting sources
    ydot(:) = vode_state % rpar(irp_ydot_a:irp_ydot_a-1+SVAR_EVOLVE)

    ! add in the reacting terms -- here we convert from dY/dt to dX/dt
    ydot(SFS:SFS-1+nspec) = ydot(SFS:SFS-1+nspec) + &
         vode_state % rpar(irp_SRHO) * aion(1:nspec) * ydot_react(1:nspec)

#if defined(SDC_EVOLVE_ENERGY)

    ydot(SEINT) = ydot(SEINT) + vode_state % rpar(irp_SRHO) * ydot_react(net_ienuc)
    ydot(SEDEN) = ydot(SEDEN) + vode_state % rpar(irp_SRHO) * ydot_react(net_ienuc)

#elif defined(SDC_EVOLVE_ENTHALPY)

    ydot(SENTH) = ydot(SENTH) + vode_state % rpar(irp_SRHO) * ydot_react(net_ienuc)

#endif

  end subroutine rhs_to_vode


  subroutine jac_to_vode(time, jac_react, vode_state, jac)

    ! this is only used with an analytic Jacobian


    ! we come in with burn_state % jac being the Jacobian of the reacting system
    ! but we need to convert it to the SDC system

    use burn_type_module, only : burn_t, net_ienuc, net_itemp, copy_burn_t, neqs
    use eos_type_module, only : eos_input_re, eos_input_rh, eos_t
    use eos_module, only : eos
    use eos_composition_module, only : eos_xderivs_t, composition_derivatives
    use actual_rhs_module

    real(rt), intent(in) :: time
    type(dvode_t), intent(inout) :: vode_state
    real(rt), intent(inout) :: jac_react(neqs, neqs)
    real(rt)    :: jac(SVAR_EVOLVE,SVAR_EVOLVE)

    integer :: m, n
    integer, parameter :: iwrho = 1, iwfs=2, iwT = iwfs+nspec, iwvar = 2+nspec

    real(rt) :: dRdw(SVAR_EVOLVE, iwvar)
    real(rt) :: dwdU(iwvar, SVAR_EVOLVE)
    type(burn_t) :: burn_state
    type(burn_t) :: burn_state_pert
    type(eos_t) :: eos_state
    type(eos_xderivs_t) :: eos_xderivs
    real(rt) :: K
    real(rt), parameter :: eps = 1.e-8_rt
    real(rt) :: ydot(neqs), ydot_pert(neqs)

#if defined(SDC_EVOLVE_ENTHALPY)
    real(rt) :: jac_temp(SVAR_EVOLVE+1, SVAR_EVOLVE+1)
    integer :: SRHO_EXTRA
#endif

    real(rt), parameter :: smallK = 1.e-15_rt

    !$gpu

    ! burn_state % jac has the derivatives with respect to the native
    ! network variables, X, T. e.  It does not have derivatives with
    ! respect to density, so we'll have to compute those ourselves.

    ! The Jacobian from the nets is in terms of dYdot/dY, but we want
    ! it was dXdot/dX, so convert here.
    do n = 1, nspec
       jac_react(n,:) = jac_react(n,:) * aion(n)
       jac_react(:,n) = jac_react(:,n) * aion_inv(n)
    enddo

    ! also fill the ydot -- we can't assume that it is valid on input
    call vode_to_burn(time, vode_state, burn_state)
    call actual_rhs(burn_state, ydot)

    ! at this point, our Jacobian should be entirely in terms of X,
    ! not Y.  Let's now fix the rhs terms themselves to be in terms of
    ! dX/dt and not dY/dt.
    ydot(1:nspec) = ydot(1:nspec) * aion(1:nspec)

#if defined(SDC_EVOLVE_ENERGY)

    ! Our jacobian, dR/dw has the form:
    !
    !  SFS         / d(rho X1dot)/drho  d(rho X1dot)/dX1   d(rho X1dit)/dX2   ...  d(rho X1dot)/dT \
    !              | d(rho X2dot)/drho  d(rho X2dot)/dX1   d(rho X2dot)/dX2   ...  d(rho X2dot)/dT |
    !  SFS-1+nspec |   ...                                                                         |
    !  SEINT       | d(rho Edot)/drho   d(rho Edot)/dX1    d(rho Edot)/dX2    ...  d(rho Edot)/dT  |
    !  SEDEN       \ d(rho Edot)/drho   d(rho Edot)/dX1    d(rho Edot)/dX2    ...  d(rho Edot)/dT  /

    dRdw(:,:) = ZERO

    ! now perturb density and call the RHS to compute the derivative wrt rho
    ! species rates come back in terms of molar fractions
    call copy_burn_t(burn_state_pert, burn_state)
    burn_state_pert % rho = burn_state % rho * (ONE + eps)

    burn_state_pert % i = burn_state % i
    burn_state_pert % j = burn_state % j
    burn_state_pert % k = burn_state % k

    call actual_rhs(burn_state_pert, ydot_pert)

    ! make the rates dX/dt and not dY/dt
    ydot_pert(1:nspec) = ydot_pert(1:nspec) * aion(1:nspec)

    ! fill the column of dRdw corresponding to the derivative
    ! with respect to rho
    do m = 1, nspec
       ! d( d(rho X_m)/dt)/drho
       dRdw(SFS-1+m, iwrho) = ydot(m) + &
            vode_state % rpar(irp_SRHO) * (ydot_pert(m) - ydot(m))/(eps * burn_state % rho)
    enddo

    ! d( d(rho e)/dt)/drho
    dRdw(SEINT, iwrho) = ydot(net_ienuc) + &
         vode_state % rpar(irp_SRHO) * (ydot_pert(net_ienuc) - ydot(net_ienuc))/(eps * burn_state % rho)

    ! d( d(rho E)/dt)/drho
    dRdw(SEDEN, iwrho) = ydot(net_ienuc) + &
         vode_state % rpar(irp_SRHO) * (ydot_pert(net_ienuc) - ydot(net_ienuc))/(eps * burn_state % rho)

    ! fill the columns of dRdw corresponding to each derivative
    ! with respect to species mass fraction
    do n = 1, nspec
       do m = 1, nspec
          ! d( d(rho X_m)/dt)/dX_n
          dRdw(SFS-1+m, iwfs-1+n) = vode_state % rpar(irp_SRHO) * jac_react(m, n)
       enddo

       ! d( d(rho e)/dt)/dX_n
       dRdw(SEINT, iwfs-1+n) = vode_state % rpar(irp_SRHO) * jac_react(net_ienuc, n)

       ! d( d(rho E)/dt)/dX_n
       dRdw(SEDEN, iwfs-1+n) = vode_state % rpar(irp_SRHO) * jac_react(net_ienuc, n)

    enddo

    ! now fill the column corresponding to derivatives with respect to
    ! temperature -- this column is iwT

    ! d( d(rho X_m)/dt)/dT
    do m = 1, nspec
       dRdw(SFS-1+m, iwT) = vode_state % rpar(irp_SRHO) * jac_react(m, net_itemp)
    enddo

    ! d( d(rho e)/dt)/dT
    dRdw(SEINT, iwT) = vode_state % rpar(irp_SRHO) * jac_react(net_ienuc, net_itemp)

    ! d( d(rho E)/dt)/dT
    dRdw(SEDEN, iwT) = vode_state % rpar(irp_SRHO) * jac_react(net_ienuc, net_itemp)

    ! that completes dRdw

    ! construct dwdU

    dwdU(:, :) = ZERO

    ! kinetic energy, K = 1/2 |U|^2
    K = 0.5_rt * sum(vode_state % rpar(irp_SMX:irp_SMZ)**2)/vode_state % rpar(irp_SRHO)**2

    if (K == ZERO) then
       K = smallK
    end if

    ! density row (iwrho)
    dwdU(iwrho, SEINT) = -1.0_rt/K
    dwdU(iwrho, SEDEN) = 1.0_rt/K

    ! species rows
    do m = 1, nspec
       dwdU(iwfs-1+m, SFS-1+m) = 1.0_rt/vode_state % rpar(irp_SRHO)
       dwdU(iwfs-1+m, SEINT) = burn_state % xn(m) / (burn_state % rho * K)
       dwdU(iwfs-1+m, SEDEN) = -burn_state % xn(m) / (burn_state % rho * K)
    end do


    eos_state % rho = vode_state % rpar(irp_SRHO)
    eos_state % T = 1.e4   ! initial guess
    eos_state % xn(:) = vode_state % y(SFS:SFS-1+nspec)/vode_state % rpar(irp_SRHO)
    eos_state % e = vode_state % y(SEINT) / vode_state % rpar(irp_SRHO)

    call eos(eos_input_re, eos_state)

    call composition_derivatives(eos_state, eos_xderivs)

    ! temperature row
    dwdU(iwT, SFS:SFS-1+nspec) = -eos_xderivs % dedX(1:nspec)/ (eos_state % rho * eos_state % dedT)
    dwdU(iwT, SEINT) = - (sum(eos_state % xn * eos_xderivs % dedX) - &
         eos_state % rho * eos_state % dedr - eos_state % e - K) / (eos_state % rho * eos_state % dedT * K)
    dwdU(iwT, SEDEN) = (sum(eos_state % xn * eos_xderivs % dedX) - &
         eos_state % rho * eos_state % dedr - eos_state % e) / (eos_state % rho * eos_state % dedT * K)

    ! compute J = dR/dw dw/dU

    ! CUDA Fortran doesn't support the matmul intrinsic :(
    jac(:,:) = matmul(dRdw, dwdU)

    ! call dgemm(1, 1, SVAR_EVOLVE, SVAR_EVOLVE, iwvar, ONE, &
    !            dRdw, SVAR_EVOLVE, &
    !            dwdU, iwvar, &
    !            ZERO, &
    !            jac, SVAR_EVOLVE)


#elif defined(SDC_EVOLVE_ENTHALPY)

    ! Our R source has components for species and enthalpy only.  But
    ! we will extend it here to include the mass density too to ensure
    ! that we have a square matrix in dU/dw that we can take the
    ! inverse of to use below.  When we compute the final Jacobian, we will
    ! discard the density row.

    ! Our jacobian, dR/dw has the form:
    !
    !  SFS         / d(rho X1dot)/drho  d(rho X1dot)/dX1   d(rho X1dit)/dX2   ...  d(rho X1dot)/dT \
    !              | d(rho X2dot)/drho  d(rho X2dot)/dX1   d(rho X2dot)/dX2   ...  d(rho X2dot)/dT |
    !  SFS-1+nspec |   ...                                                                         |
    !  SENTH       | d(rho h)/drho      d(rho h)/dX1       d(rho h)/dX2       ...  d(rho h)/dT     |
    !  SRHO_EXTRA  \ 0                  0                  0                       0               /

    SRHO_EXTRA = SENTH + 1

    dRdw(:,:) = ZERO

    ! now perturb density and call the RHS to compute the derivative wrt rho
    ! species rates come back in terms of molar fractions
    call copy_burn_t(burn_state_pert, burn_state)
    burn_state_pert % rho = burn_state % rho * (ONE + eps)

    burn_state_pert % i = burn_state % i
    burn_state_pert % j = burn_state % j
    burn_state_pert % k = burn_state % k

    call actual_rhs(burn_state_pert, ydot_pert)

    ! make the rates dX/dt and not dY/dt
    ydot_pert(1:nspec) = ydot_pert(1:nspec) * aion(1:nspec)

    ! fill the column of dRdw corresponding to the derivative
    ! with respect to rho
    do m = 1, nspec
       ! d( d(rho X_m)/dt)/drho
       dRdw(SFS-1+m, iwrho) = ydot(m) + &
            vode_state % rpar(irp_SRHO) * (ydot_pert(m) - ydot(m))/(eps * burn_state % rho)
    enddo

    ! d( d(rho h)/dt)/drho
    dRdw(SENTH, iwrho) = ydot(net_ienuc) + &
         vode_state % rpar(irp_SRHO) * (ydot_pert(net_ienuc) - ydot(net_ienuc))/(eps * burn_state % rho)

    ! d( d(rho)/dt)/drho
    dRdw(SRHO_EXTRA, iwrho) = ZERO

    ! fill the columns of dRdw corresponding to each derivative
    ! with respect to species mass fraction
    do n = 1, nspec
       do m = 1, nspec
          ! d( d(rho X_m)/dt)/dX_n
          dRdw(SFS-1+m, iwfs-1+n) = vode_state % rpar(irp_SRHO) * jac_react(m, n)
       enddo

       ! d( d(rho h)/dt)/dX_n
       dRdw(SENTH, iwfs-1+n) = vode_state % rpar(irp_SRHO) * jac_react(net_ienuc, n)

       ! d( d(rho)/dt)/dX_n
       dRdw(SRHO_EXTRA, iwfs-1+n) = ZERO

    enddo

    ! now fill the column corresponding to derivatives with respect to
    ! temperature -- this column is iwT

    ! d( d(rho X_m)/dt)/dT
    do m = 1, nspec
       dRdw(SFS-1+m, iwT) = vode_state % rpar(irp_SRHO) * jac_react(m, net_itemp)
    enddo

    ! d( d(rho h)/dt)/dT
    dRdw(SENTH, iwT) = vode_state % rpar(irp_SRHO) * jac_react(net_ienuc, net_itemp)

    ! d( d(rho)/dt)/dT
    dRdw(SRHO_EXTRA, iwT) = ZERO

    ! that completes dRdw

    ! construct dwdU.  Here we take U = (rho X, rho h, rho)^T

    dwdU(:, :) = ZERO

    ! density row (iwrho)
    dwdU(iwrho, SRHO_EXTRA) = 1.0_rt

    ! species rows
    do m = 1, nspec
       dwdU(iwfs-1+m, SFS-1+m) = 1.0_rt/vode_state % rpar(irp_SRHO)
       dwdU(iwfs-1+m, SRHO_EXTRA) = -burn_state % xn(m) / burn_state % rho
    end do


    eos_state % rho = vode_state % rpar(irp_SRHO)
    eos_state % T = 1.e4   ! initial guess
    eos_state % xn(:) = y(SFS:SFS-1+nspec)/vode_state % rpar(irp_SRHO)
    eos_state % h = y(SENTH)/vode_state % rpar(irp_SRHO)

    call eos(eos_input_rh, eos_state)

    call composition_derivatives(eos_state, eos_xderivs)

    ! temperature row
    dwdU(iwT, SFS:SFS-1+nspec) = -eos_xderivs % dhdX(1:nspec)/ (eos_state % rho * eos_state % dedT)
    dwdU(iwT, SENTH) = ONE/(eos_state % rho * eos_state % dhdT)
    dwdU(iwT, SRHO_EXTRA) = (sum(eos_state % xn * eos_xderivs % dhdX) - &
         eos_state % rho * eos_state % dhdr - eos_state % h) / (eos_state % rho * eos_state % dhdT)

    ! compute J_temp = dR/dw dw/dU

    ! CUDA Fortran doesn't support the matmul intrinsic :(
    jac_temp(:,:) = matmul(dRdw, dwdU)
    ! call dgemm(1, 1, SVAR_EVOLVE+1, SVAR_EVOLVE+1, iwvar, ONE, &
    !            dRdw, SVAR_EVOLVE+1, &
    !            dwdU, iwvar, &
    !            ZERO, &
    !            jac_temp, SVAR_EVOLVE+1)

    ! now discard the SRHO_EXTRA row / column
    do n = 1, SVAR_EVOLVE
       do m = 1, SVAR_EVOLVE
          jac(m, n) = jac_temp(m, n)
       end do
    end do

#endif

  end subroutine jac_to_vode


  subroutine vode_to_burn(time, vode_state, burn_state)

    use eos_type_module, only : eos_t, eos_input_re, eos_input_rt, eos_input_rp, eos_input_rh
    use eos_type_module, only : eos_get_small_temp, eos_get_max_temp
    use eos_module, only : eos
    use burn_type_module, only : eos_to_burn, burn_t

#ifdef SDC_EVOLVE_ENTHALPY
    use meth_params_module, only: use_tfromp
#endif

    type (burn_t) :: burn_state
    real(rt), intent(in) :: time
    type(dvode_t), intent(inout) :: vode_state

    type(eos_t) :: eos_state

    real(rt) :: rhoInv, min_temp, max_temp

    !$gpu

    ! update rho, rho*u, etc.
    call fill_unevolved_variables(time, vode_state)

    rhoInv = ONE / vode_state % rpar(irp_SRHO)

    eos_state % rho = vode_state % rpar(irp_SRHO)
    eos_state % xn  = vode_state % y(SFS:SFS+nspec-1) * rhoInv

#if defined(SDC_EVOLVE_ENERGY)

    if (vode_state % rpar(irp_T_from_eden) > ZERO) then
       eos_state % e = (vode_state % y(SEDEN) - &
            HALF*rhoInv*sum(vode_state % rpar(irp_SMX:irp_SMZ)**2)) * rhoInv
    else
       eos_state % e = vode_state % y(SEINT) * rhoInv
    endif

#elif defined(SDC_EVOLVE_ENTHALPY)

    if (use_tfromp) then
       ! NOT SURE IF THIS IS VALID
       eos_state % p = vode_state % rpar(irp_p0)
    else
       eos_state % h = vode_state % y(SENTH) * rhoInv
    endif

#endif

    ! Give the temperature an initial guess -- use the geometric mean
    ! of the minimum and maximum temperatures.

    call eos_get_small_temp(min_temp)
    call eos_get_max_temp(max_temp)

    eos_state % T = sqrt(min_temp * max_temp)

#if defined(SDC_EVOLVE_ENERGY)

    call eos(eos_input_re, eos_state)

#elif defined(SDC_EVOLVE_ENTHALPY)

    if (use_tfromp) then
       ! NOT SURE IF THIS IS VALID
       ! used to be an Abort statement
       call eos(eos_input_rp, eos_state)
    else
       call eos(eos_input_rh, eos_state)
    endif

#endif

    call eos_to_burn(eos_state, burn_state)

    burn_state % time = time

    if (vode_state % rpar(irp_self_heat) > ZERO) then
       burn_state % self_heat = .true.
    else
       burn_state % self_heat = .false.
    endif

#ifdef SDC_EVOLVE_ENTHALPY

    burn_state % p0 = vode_state % rpar(irp_p0)

#endif

#ifdef NONAKA_PLOT

    burn_state % i = vode_state % rpar(irp_i)
    burn_state % j = vode_state % rpar(irp_j)
    burn_state % k = vode_state % rpar(irp_k)

#endif

  end subroutine vode_to_burn

end module vode_type_module
