#ifndef VODE_TYPE_SIMPLIFIED_SDC_H
#define VODE_TYPE_SIMPLIFIED_SDC

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void clean_state(const Real time, burn_t& state, dvode_t& vode_state)
{

    // update rho, rho*u, etc.

    fill_unevolved_variables(time, state, vode_state);

    // Ensure that mass fractions always stay positive.

   for (int n = 1; n <= NumSpec; ++n) {
     // we use 1-based indexing, so we need to offset SFS
     vode_state.y[SFS+n] = amrex::max(amrex::min(vode_state.y[SFS+n], state.rho),
                                      state.rho * 1.e-200_rt);
    }


    // renormalize abundances as necessary

    if (renormalize_abundances) {then
        renormalize_species(time, vode_state);
    }

#ifdef SDC_EVOLVE_ENERGY

    // Ensure that internal energy never goes above the maximum limit
    // provided by the EOS. Same for the internal energy implied by the
    // total energy (which we get by subtracting kinetic energy).

    eos_t eos_state;
    eos_state.rho = state.y[SRHO];
    eos_state.T = MAX_TEMP;
    for (int n = 0; n < NumSpec; n++) {
      // vode_state.y uses 1-based indexing
      eos_state.xn[n] = vode_state.y[SFS+n] / state.y[SRHO];
    }
#if NAUX > 0
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = state.y[SFX+n] / state.y[SRHO];
    }
#endif

    eos(eos_input_rt, eos_state);

    Real max_e = eos_state.e;

    vode_state.y[SEINT+1] = amrex::min(state.y[SRHO] * max_e, vode_state.y[SEINT+1]);

    Real ke = vode_state.y[SEDEN+1] -
      0.5_rt * (state.y[SMX] * state.y[SMX] +
                state.y[SMY] * state.y[SMY] +
                state.y[SMZ] * state.y[SMZ]) / state.y[SRHO];

    vode_state.y[SEDEN+1] = amrex::min(state.y[SRHO] * max_e + ke, vode_state.y[SEDEN+1]);

#endif

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void fill_unevolved_variables(const Real time, burn_t& state, dvode_t& vode_state)
{

  // some quantities are only advected.  Here we update those state
  // variables in burn_t's y[] array due only to advection.

#if defined(SDC_EVOLVE_ENERGY)

  // we are always integrating from t = 0, so there is no offset
  // time needed here.  The indexing of ydot_a is based on
  // the indices in burn_t and is 0-based
  state.y[SRHO] = state.rho_orig + state.ydot_a[SRHO] * time;

  // for consistency
  state.rho = state.y[SRHO];

  state.y[SMX] = state.umx_orig + state.ydot_a[SMX] * time;
  state.y[SMY] = state.umy_orig + state.ydot_a[SMY] * time;
  state.y[SMZ] = state.umz_orig + state.ydot_a[SMZ] * time;

#elif defined(SDC_EVOLVE_ENTHALPY)

  // Keep density consistent with the partial densities.
  state.rho = 0.0_rt;
  for (int n = 1; n <= NumSpec; n++) {
    // use 1-based indexing into vode_state.y
    state.rho += vode_state.y[SFS+n];
  }
#endif


#if NAUX > 0
      amrex::Error("Error: fill_unevolved_variables does not handle aux yet");
#endif

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void renormalize_species(const Real time, burn_t& state, dvode_t& vode_state)
{

    // We only renormalize species when evolving energy because
    // when we evolve enthalpy, we define the density as
    // the sum of the partial densities rho*X for each species.

#ifdef SDC_EVOLVE_ENERGY

    // update rho, rho*u, etc.

    fill_unevolved_variables(time, state, vode_state);

    Real nspec_sum = 0.0_rt;
    for (int n = 1; n <= NumSpec; n++) {
      // use 1-based indexing
      nspec_sum += vode_state.y[SFS+n];
    }
    nspec_sum /= state.y[SRHO];

    for (int n = 1; n <= NumSpec; n++) {
      vode_state.y[SFS+n] /= nspec_sum;
    }

#endif

}


// Given a burn state, fill the rpar and integration state data.
// this is indended only to be called once -- at the very start of the integration

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void burn_to_vode(burn_t& state, dvode_t& vode_state)
{

  for (int n = 0; n < SVAR_EVOLVE; n++) {
    // vode_state uses 1-based indexing
    vode_state.y[n+1] = state.y[n];
  }

    // store the original rho and rho U
#if defined(SDC_EVOLVE_ENERGY)

    state.rho_orig = state.y[SRHO];

    state.umx_orig = state.y[SMX];
    state.umy_orig = state.y[SMY];
    state.umz_orig = state.y[SMZ];

#else

    state.rho_orig = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
      state.rho_orig += state.y[SFS+n];
    }

#endif

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void vode_to_burn(const Real time, dvode_t& vode_state, burn_t& state)
{

  // this makes burn_t represent the current integration state.  The
  // main thing this really does is compute the current temperature

  for (int n = 0; n < SVAR_EVOLVE; n++) {
    // note vode_state uses 1-based indexing
    state.y[n] = vode_state.y[n+1];
  }

  // update rho, rho*u, etc.

  fill_unevolved_variables(time, state, vode_state);

  Real rhoInv = 1.0_rt / state.rho;

  eos_t eos_state;
  eos_state.rho = state.rho;
  for (int n = 0; n < NumState; n++) {
    // vode_state uses 1-based indexing
    eos_state.xn[n] = vode_state.y[SFS+1+n] * rhoInv;
  }

#if defined(SDC_EVOLVE_ENERGY)

  if (state.T_from_eden) {
    eos_state.e = (vode_state.y[SEDEN+1] -
                   0.5_rt * rhoInv * (state.y[SMX] * state.y[SMX] +
                                      state.y[SMY] * state.y[SMY] +
                                      state.y[SMZ] * state.y[SMZ])) * rhoInv;
  } else {
    eos_state.e = vode_state.y[SEINT-1] * rhoInv;
  }

#elif defined(SDC_EVOLVE_ENTHALPY)

  if (use_tfromp) {
    // under this assumption, we are assuming that there is no base state evolution
    eos_state.p = state.p0;
  } else {
    eos_state.h = vode_state.y[SENTH+1] * rhoInv;
  }

#endif

  // Give the temperature an initial guess -- use the geometric mean
  // of the minimum and maximum temperatures.

  // TODO: we should just change this to use state.T

  Real min_temp;
  eos_get_small_temp(&min_temp);

  Real max_temp;
  eos_get_max_temp(&max_temp);

  eos_state.T = std::sqrt(min_temp * max_temp);

#if defined(SDC_EVOLVE_ENERGY)

  eos(eos_input_re, eos_state);

#elif defined(SDC_EVOLVE_ENTHALPY)

  if (use_tfromp) {
    // NOT SURE IF THIS IS VALID
    // used to be an Abort statement
    eos(eos_input_rp, eos_state);
  } else {
    eos(eos_input_rh, eos_state);
  }

#endif

  eos_to_burn(eos_state, state);

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rhs_to_vode(const Real time,
                 burn_t& state, YdotNetArray1D& ydot_react, dvode_t& vode_state,
                 RArray1D& ydot)
{


  fill_unevolved_variables(time, state, vode_state);

  // ydot_react has just the contribution to the RHS from the
  // reaction network.  Note that these are in terms of dY/dt

  // start with the contribution from the non-reacting sources

  // Note: both ydot and ydot_react are 1-based
  for (int n = 0; n < SVAR_EVOLVE; n++) {
    ydot(n+1) = state.ydot_a[n];
  }

  // add in the reacting terms -- here we convert from dY/dt to dX/dt
  // The species derivatives are the first NumSpec components of
  // ydot_react

  for (int n = 1; n <= NumSpec; n++) {
    ydot(SFS+n) += state.rho * aion[n-1] * ydot_react(n);
  }

#if defined(SDC_EVOLVE_ENERGY)

  ydot(SEINT+1) += state.rho * ydot_react(net_ienuc);
  ydot(SEDEN+1) += state.rho * ydot_react(net_ienuc);

#elif defined(SDC_EVOLVE_ENTHALPY)

  ydot(SENTH+1) += state.rho * ydot_react(net_ienuc);

#endif

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void jac_to_vode(time, jac_react, vode_state, jac)
{

  // this is only used with an analytic Jacobian


    // we come in with burn_state % jac being the Jacobian of the reacting system
    // but we need to convert it to the SDC system

    use burn_type_module, only : burn_t, net_ienuc, net_itemp, copy_burn_t, neqs
    use eos_type_module, only : eos_input_re, eos_input_rh, eos_t
    use eos_module, only : eos
    use eos_composition_module, only : eos_xderivs_t, composition_derivatives
    use actual_rhs_module

    real(rt), intent(in) :: time
    type(dvode_t), intent(inout) :: vode_state
    real(rt), intent(inout) :: jac_react(neqs, neqs)
    real(rt)    :: jac(SVAR_EVOLVE,SVAR_EVOLVE)

    integer :: m, n, k
#if defined(SDC_EVOLVE_ENERGY)
    integer, parameter :: iwrho = 1, iwfs=2, iwK = iwfs+nspec, iwT = iwK+1, iwvar = 3+nspec
#else
    integer, parameter :: iwrho = 1, iwfs=2, iwT = iwfs+nspec, iwvar = 2+nspec
#endif

    // SVAR_EVOLVE doesn't include rho, but we will include it here in
    // the intermediate this affects both the Castro
    // (SDC_EVOLVE_ENERGY) and MAESTROeX (SDC_EVOLVE_ENTHALPY) systems.
    real(rt) :: dRdw(SVAR_EVOLVE+1, iwvar)
    real(rt) :: dwdU(iwvar, SVAR_EVOLVE+1)

    type(burn_t) :: burn_state
    type(burn_t) :: burn_state_pert
    type(eos_t) :: eos_state
    type(eos_xderivs_t) :: eos_xderivs
    real(rt) :: kineng
    real(rt), parameter :: eps = 1.e-8_rt
    real(rt) :: ydot(neqs), ydot_pert(neqs)

    integer :: SRHO_EXTRA

    real(rt), parameter :: smallK = 1.e-15_rt

    //$gpu

    // burn_state % jac has the derivatives with respect to the native
    // network variables, X, T. e.  It does not have derivatives with
    // respect to density, so we'll have to compute those ourselves.

    // The Jacobian from the nets is in terms of dYdot/dY, but we want
    // it was dXdot/dX, so convert here.
    do n = 1, nspec
       jac_react(n,:) = jac_react(n,:) * aion(n)
       jac_react(:,n) = jac_react(:,n) * aion_inv(n)
    enddo

    // also fill the ydot -- we can't assume that it is valid on input
    call vode_to_burn(time, vode_state, burn_state)
    call actual_rhs(burn_state, ydot)

    // at this point, our Jacobian should be entirely in terms of X,
    // not Y.  Let's now fix the rhs terms themselves to be in terms of
    // dX/dt and not dY/dt.
    ydot(1:nspec) = ydot(1:nspec) * aion(1:nspec)

    SRHO_EXTRA = SVAR_EVOLVE + 1

    dRdw(:,:) = ZERO
    dwdU(:, :) = ZERO

#if defined(SDC_EVOLVE_ENERGY)

    // The system we integrate has the form (rho X_k, rho E, rho e), but we will temporarily augment
    // This with rho, giving U = (rho, rho X_k, rho E, rho e).
    //
    // The intermediate state, w, has the form w = (rho, X_k, K, T), where K is 1/2 |U|^2

    // First compute dR/dw using the Jacobian that comes from the
    // network.  Note: this doesn't include density derivatives, so we
    // compute those via differencing.

    // dR/dw has the form:
    //
    //  SFS         / d(rho X1dot)/drho  d(rho X1dot)/dX1   d(rho X1dit)/dX2   ...  0   d(rho X1dot)/dT \
    //              | d(rho X2dot)/drho  d(rho X2dot)/dX1   d(rho X2dot)/dX2   ...  0   d(rho X2dot)/dT |
    //  SFS-1+nspec |   ...                                                         0                   |
    //  SEINT       | d(rho Edot)/drho   d(rho Edot)/dX1    d(rho Edot)/dX2    ...  0   d(rho Edot)/dT  |
    //  SEDEN       | d(rho Edot)/drho   d(rho Edot)/dX1    d(rho Edot)/dX2    ...  0   d(rho Edot)/dT  |
    //  SRHO_EXTRA  \        0                  0                  0                0          0       /
    //
    //                                                                              ^
    //                                                                              K derivatives

    // now perturb density and call the RHS to compute the derivative wrt rho
    // species rates come back in terms of molar fractions
    call copy_burn_t(burn_state_pert, burn_state)
    burn_state_pert % rho = burn_state % rho * (ONE + eps)

    burn_state_pert % i = burn_state % i
    burn_state_pert % j = burn_state % j
    burn_state_pert % k = burn_state % k

    call actual_rhs(burn_state_pert, ydot_pert)

    // make the rates dX/dt and not dY/dt
    ydot_pert(1:nspec) = ydot_pert(1:nspec) * aion(1:nspec)

    // fill the column of dRdw corresponding to the derivative
    // with respect to rho
    do m = 1, nspec
       // d( d(rho X_m)/dt)/drho
       dRdw(SFS-1+m, iwrho) = ydot(m) + &
            vode_state % rpar(irp_SRHO) * (ydot_pert(m) - ydot(m))/(eps * burn_state % rho)
    enddo

    // d( d(rho e)/dt)/drho
    dRdw(SEINT, iwrho) = ydot(net_ienuc) + &
         vode_state % rpar(irp_SRHO) * (ydot_pert(net_ienuc) - ydot(net_ienuc))/(eps * burn_state % rho)

    // d( d(rho E)/dt)/drho
    dRdw(SEDEN, iwrho) = ydot(net_ienuc) + &
         vode_state % rpar(irp_SRHO) * (ydot_pert(net_ienuc) - ydot(net_ienuc))/(eps * burn_state % rho)

    // fill the columns of dRdw corresponding to each derivative
    // with respect to species mass fraction
    do n = 1, nspec
       do m = 1, nspec
          // d( d(rho X_m)/dt)/dX_n
          dRdw(SFS-1+m, iwfs-1+n) = vode_state % rpar(irp_SRHO) * jac_react(m, n)
       enddo

       // d( d(rho e)/dt)/dX_n
       dRdw(SEINT, iwfs-1+n) = vode_state % rpar(irp_SRHO) * jac_react(net_ienuc, n)

       // d( d(rho E)/dt)/dX_n
       dRdw(SEDEN, iwfs-1+n) = vode_state % rpar(irp_SRHO) * jac_react(net_ienuc, n)

    enddo

    // now fill the column corresponding to derivatives with respect to
    // temperature -- this column is iwT

    // d( d(rho X_m)/dt)/dT
    do m = 1, nspec
       dRdw(SFS-1+m, iwT) = vode_state % rpar(irp_SRHO) * jac_react(m, net_itemp)
    enddo

    // d( d(rho e)/dt)/dT
    dRdw(SEINT, iwT) = vode_state % rpar(irp_SRHO) * jac_react(net_ienuc, net_itemp)

    // d( d(rho E)/dt)/dT
    dRdw(SEDEN, iwT) = vode_state % rpar(irp_SRHO) * jac_react(net_ienuc, net_itemp)

    // for the K derivatives, dRdw(:, iwK), and the rho sources,
    // dRdw(SRHO_EXTRA, :), we don't need to do anything, because these
    // are already zeroed out

    // that completes dRdw


    // construct dwdU

    // kinetic energy, K = 1/2 |U|^2
    kineng = 0.5_rt * sum(vode_state % rpar(irp_SMX:irp_SMZ)**2)/vode_state % rpar(irp_SRHO)**2

    // density row (iwrho)
    dwdU(iwrho, SRHO_EXTRA) = 1.0_rt

    // species rows
    do m = 1, nspec
       dwdU(iwfs-1+m, SFS-1+m) = 1.0_rt/vode_state % rpar(irp_SRHO)
       dwdU(iwfs-1+m, SRHO_EXTRA) = -burn_state % xn(m) / burn_state % rho
    end do

    // K row
    dwdU(iwK, SRHO_EXTRA) = -kineng / burn_state % rho
    dwdU(iwK, SEINT) = -1.0_rt / burn_state % rho
    dwdU(iwK, SEDEN) = 1.0_rt / burn_state % rho

    // T row
    eos_state % rho = vode_state % rpar(irp_SRHO)
    eos_state % T = 1.e4   // initial guess
    eos_state % xn(:) = vode_state % y(SFS:SFS-1+nspec)/vode_state % rpar(irp_SRHO)
    eos_state % e = vode_state % y(SEINT) / vode_state % rpar(irp_SRHO)

    call eos(eos_input_re, eos_state)

    call composition_derivatives(eos_state, eos_xderivs)

    // temperature row
    dwdU(iwT, SFS:SFS-1+nspec) = -eos_xderivs % dedX(1:nspec)/ (eos_state % rho * eos_state % dedT)
    dwdU(iwT, SEINT) = 1.0_rt / (eos_state % rho * eos_state % dedT)
    dwdU(iwT, SEDEN) = 0.0_rt
    dwdU(iwT, SRHO_EXTRA) = &
         (sum(eos_state % xn * eos_xderivs % dedX) - eos_state % rho * eos_state % dedr - eos_state % e) / &
         (eos_state % rho * eos_state % dedT)

#elif defined(SDC_EVOLVE_ENTHALPY)

    // Our R source has components for species and enthalpy only.  But
    // we will extend it here to include the mass density too to ensure
    // that we have a square matrix in dU/dw that we can take the
    // inverse of to use below.  When we compute the final Jacobian, we will
    // discard the density row.

    // Our jacobian, dR/dw has the form:
    //
    //  SFS         / d(rho X1dot)/drho  d(rho X1dot)/dX1   d(rho X1dit)/dX2   ...  d(rho X1dot)/dT \
    //              | d(rho X2dot)/drho  d(rho X2dot)/dX1   d(rho X2dot)/dX2   ...  d(rho X2dot)/dT |
    //  SFS-1+nspec |   ...                                                                         |
    //  SENTH       | d(rho h)/drho      d(rho h)/dX1       d(rho h)/dX2       ...  d(rho h)/dT     |
    //  SRHO_EXTRA  \ 0                  0                  0                       0               /


    // now perturb density and call the RHS to compute the derivative wrt rho
    // species rates come back in terms of molar fractions
    call copy_burn_t(burn_state_pert, burn_state)
    burn_state_pert % rho = burn_state % rho * (ONE + eps)

    burn_state_pert % i = burn_state % i
    burn_state_pert % j = burn_state % j
    burn_state_pert % k = burn_state % k

    call actual_rhs(burn_state_pert, ydot_pert)

    // make the rates dX/dt and not dY/dt
    ydot_pert(1:nspec) = ydot_pert(1:nspec) * aion(1:nspec)

    // fill the column of dRdw corresponding to the derivative
    // with respect to rho
    do m = 1, nspec
       // d( d(rho X_m)/dt)/drho
       dRdw(SFS-1+m, iwrho) = ydot(m) + &
            vode_state % rpar(irp_SRHO) * (ydot_pert(m) - ydot(m))/(eps * burn_state % rho)
    enddo

    // d( d(rho h)/dt)/drho
    dRdw(SENTH, iwrho) = ydot(net_ienuc) + &
         vode_state % rpar(irp_SRHO) * (ydot_pert(net_ienuc) - ydot(net_ienuc))/(eps * burn_state % rho)

    // d( d(rho)/dt)/drho
    dRdw(SRHO_EXTRA, iwrho) = ZERO

    // fill the columns of dRdw corresponding to each derivative
    // with respect to species mass fraction
    do n = 1, nspec
       do m = 1, nspec
          // d( d(rho X_m)/dt)/dX_n
          dRdw(SFS-1+m, iwfs-1+n) = vode_state % rpar(irp_SRHO) * jac_react(m, n)
       enddo

       // d( d(rho h)/dt)/dX_n
       dRdw(SENTH, iwfs-1+n) = vode_state % rpar(irp_SRHO) * jac_react(net_ienuc, n)

       // d( d(rho)/dt)/dX_n
       dRdw(SRHO_EXTRA, iwfs-1+n) = ZERO

    enddo

    // now fill the column corresponding to derivatives with respect to
    // temperature -- this column is iwT

    // d( d(rho X_m)/dt)/dT
    do m = 1, nspec
       dRdw(SFS-1+m, iwT) = vode_state % rpar(irp_SRHO) * jac_react(m, net_itemp)
    enddo

    // d( d(rho h)/dt)/dT
    dRdw(SENTH, iwT) = vode_state % rpar(irp_SRHO) * jac_react(net_ienuc, net_itemp)

    // d( d(rho)/dt)/dT
    dRdw(SRHO_EXTRA, iwT) = ZERO

    // that completes dRdw

    // construct dwdU.  Here we take U = (rho X, rho h, rho)^T

    // density row (iwrho)
    dwdU(iwrho, SRHO_EXTRA) = 1.0_rt

    // species rows
    do m = 1, nspec
       dwdU(iwfs-1+m, SFS-1+m) = 1.0_rt/vode_state % rpar(irp_SRHO)
       dwdU(iwfs-1+m, SRHO_EXTRA) = -burn_state % xn(m) / burn_state % rho
    end do


    eos_state % rho = vode_state % rpar(irp_SRHO)
    eos_state % T = 1.e4   // initial guess
    eos_state % xn(:) = vode_state % y(SFS:SFS-1+nspec)/vode_state % rpar(irp_SRHO)
    eos_state % h = vode_state % y(SENTH)/vode_state % rpar(irp_SRHO)

    call eos(eos_input_rh, eos_state)

    call composition_derivatives(eos_state, eos_xderivs)

    // temperature row
    dwdU(iwT, SFS:SFS-1+nspec) = -eos_xderivs % dhdX(1:nspec)/ (eos_state % rho * eos_state % dedT)
    dwdU(iwT, SENTH) = ONE/(eos_state % rho * eos_state % dhdT)
    dwdU(iwT, SRHO_EXTRA) = (sum(eos_state % xn * eos_xderivs % dhdX) - &
         eos_state % rho * eos_state % dhdr - eos_state % h) / (eos_state % rho * eos_state % dhdT)

#endif


    // compute J = dR/dw dw/dU

    // J is SVAR_EVOLVE x SVAR_EVOLVE, which will call m x n
    //
    // J = dR/dw dw/dU
    //
    //   dR/dw is SVAR_EVOLVE+1 x iwvar, which we call m x k
    //   dw/dU is iwvar x SVAR_EVOLVE+1, which we call k x n
    //

    // we need to cut out the density (SRHO_EXTRA) row and column of
    // the Jacobian, since that is not in our full SVAR_EVOLVE state
    do n = 1, SVAR_EVOLVE
       if (n == SRHO_EXTRA) cycle
       do m = 1, SVAR_EVOLVE
          if (m == SRHO_EXTRA) cycle

          jac(m, n) = 0.0_rt
          do k = 1, iwvar
             jac(m, n) = jac(m,n) + dRdw(m, k) * dwdU(k, n)
          end do
       end do
    end do

}



#endif

