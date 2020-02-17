module temperature_integration_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  logical, save, allocatable :: self_heat

#ifdef AMREX_USE_CUDA
  attributes(managed) :: self_heat
#endif

  !$acc declare create(self_heat)

  public

contains

  ! Sets up the temperature equation. This should be called from
  ! within the actual_rhs routine but is provided here as a convenience
  ! since most networks will use the same temperature ODE.

  subroutine temperature_rhs(state, ydot)

    !$acc routine seq

    use amrex_constants_module, only: ZERO, ONE
    use network, only: nspec
    use burn_type_module
    use jacobian_sparsity_module, only: get_jac_entry, set_jac_entry
    use extern_probin_module, only: dT_crit, call_eos_in_rhs

    implicit none

    type (burn_t) :: state
    real(rt), intent(inout) :: ydot(neqs)

    real(rt) :: cx, cx_inv

    !$gpu

    if (state % self_heat) then

       ! Set up the temperature ODE.  For constant pressure, Dp/Dt = 0, we
       ! evolve :
       !    dT/dt = (1/c_p) [ Hnuc ]
       !
       ! For constant volume, div{U} = 0, and we evolve:
       !    dT/dt = (1/c_v) [ Hnuc ]
       !
       ! See low Mach paper III, including Eq. A3 for details.
       ! Note that we no longer include the chemical potential (dE/dX or dH/dX)
       ! terms because we believe they analytically should vanish.
       !
       ! The burn_t "cx" field is the correct specific heat based on
       ! do_constant_volume_burn.

       if (.not. call_eos_in_rhs .and. dT_crit < 1.0e19_rt) then
          cx = state % cx + (state % T - state % T_old) * state % dcxdt

       else
          cx = state % cx

       endif

       cx_inv = ONE / cx
       ydot(net_itemp) = ydot(net_ienuc) * cx_inv

    endif

  end subroutine temperature_rhs



  ! Sets up the temperature entries in the Jacobian. This should be called from
  ! within the actual_jac routine but is provided here as a convenience
  ! since most networks will use the same temperature ODE.

  subroutine temperature_jac(state, jac)

    !$acc routine seq

    use amrex_constants_module, only: ZERO, ONE
    use network, only: nspec
    use burn_type_module
    use jacobian_sparsity_module, only: get_jac_entry, set_jac_entry
    use extern_probin_module, only: dT_crit, call_eos_in_rhs

    implicit none

    type (burn_t) :: state
    real(rt) :: jac(neqs, neqs)
    real(rt) :: scratch, cspec, cspec_inv

    integer :: k

    !$gpu

    ! Temperature Jacobian elements

    if (state % self_heat) then

       if (.not. call_eos_in_rhs .and. dT_crit < 1.0e19_rt) then

          cspec = state % cx + (state % T - state % T_old) * state % dxvdt

       else

          cspec = state % cx

       endif

       cspec_inv = ONE / cspec

       ! d(itemp)/d(yi)

       do k = 1, nspec_evolve
          call get_jac_entry(jac, net_ienuc, k, scratch)
          scratch = scratch * cspecInv
          call set_jac_entry(jac, net_itemp, k, scratch)
       enddo

       ! d(itemp)/d(temp) -- we get this from the equation for d (denuc / dt) / dT
       ! since dT/dt = 1/c_x denuc/dt in our formalism

       call get_jac_entry(jac, net_ienuc, net_itemp, scratch)
       scratch = scratch * cspec_inv
       call set_jac_entry(jac, net_itemp, net_itemp, scratch)

       ! d(itemp)/d(enuc)

       scratch = ZERO
       call set_jac_entry(jac, net_itemp, net_ienuc, scratch)

    endif

  end subroutine temperature_jac



  subroutine temperature_rhs_init()

    use extern_probin_module, only: burning_mode
    use amrex_error_module, only: amrex_error

    implicit none

    ! Provide a default value, then consult the burning_mode.

    allocate(self_heat)
    self_heat = .true.

    if (burning_mode == 0 .or. burning_mode == 2) then
       self_heat = .false.
    else if (burning_mode == 1 .or. burning_mode == 3) then
       self_heat = .true.
    else
       call amrex_error("Error: unknown burning_mode in temperature_rhs_init()")
    end if

    !$acc update device(self_heat)

  end subroutine temperature_rhs_init

end module temperature_integration_module
