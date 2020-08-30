module eos_composition_module

  use eos_type_module, only : eos_t
#ifdef NSE_THERMO
  use network, only : nspec, iye, iabar
#else
  use network, only: nspec
#endif
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  type :: eos_xderivs_t
    real(rt) :: dedX(nspec)
    real(rt) :: dpdX(nspec)
    real(rt) :: dhdX(nspec)
 end type eos_xderivs_t

contains

  ! Given a set of mass fractions, calculate quantities that depend
  ! on the composition like abar and zbar.

  subroutine composition(state)

    use amrex_constants_module, only: ONE
    use network, only: aion_inv, zion

    implicit none

    type (eos_t), intent(inout) :: state

    !$gpu

    ! Calculate abar, the mean nucleon number,
    ! zbar, the mean proton number,
    ! mu, the mean molecular weight,
    ! mu_e, the mean number of nucleons per electron, and
    ! y_e, the electron fraction.

#ifdef NSE_THERMO
    state % mu_e = ONE / state % aux(iye)
    state % y_e = state % aux(iye)
    state % abar = state% aux(iabar)
    state % zbar = state % abar * state % y_e

#else
    state % y_e = sum(state % xn(:) * zion(:) * aion_inv(:))
    state % mu_e = ONE / state % y_e

    state % abar = ONE / (sum(state % xn(:) * aion_inv(:)))
    state % zbar = state % abar * state % y_e
#endif

  end subroutine composition


  ! Compute thermodynamic derivatives with respect to xn(:)

  subroutine composition_derivatives(state, state_xderivs)

    use amrex_constants_module, only: ZERO
#ifdef EXTRA_THERMO
    use network, only: aion, aion_inv, zion
#endif

    implicit none

    type (eos_t), intent(in) :: state
    type (eos_xderivs_t), intent(out) :: state_xderivs

    !$gpu

#ifdef EXTRA_THERMO
    state_xderivs % dpdX(:) = state % dpdA * (state % abar * aion_inv(:))   &
                                        * (aion(:) - state % abar) &
                            + state % dpdZ * (state % abar * aion_inv(:))   &
                                        * (zion(:) - state % zbar)

    state_xderivs % dEdX(:) = state % dedA * (state % abar * aion_inv(:))   &
                                        * (aion(:) - state % abar) &
                            + state % dedZ * (state % abar * aion_inv(:))   &
                                        * (zion(:) - state % zbar)

    if (state % dPdr .ne. ZERO) then

       state_xderivs % dhdX(:) = state_xderivs % dedX(:) &
            + (state % p / state % rho**2 - state % dedr) &
            *  state_xderivs % dPdX(:) / state % dPdr

    endif
#endif

  end subroutine composition_derivatives

end module eos_composition_module
