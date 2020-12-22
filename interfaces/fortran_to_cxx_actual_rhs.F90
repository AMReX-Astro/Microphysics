module fortran_to_cxx_actual_rhs_module

  use amrex_fort_module, only: rt => amrex_real
  use network, only: nspec, naux
  use burn_type_module, only: burn_t, neqs

  implicit none

contains

  subroutine fortran_to_cxx_actual_rhs(rho, T, e, xn, abar, zbar, y_e, eta, &
#if NAUX_NET > 0
                                       aux, &
#endif
#ifndef SIMPLIFIED_SDC
                                       cv, cp, self_heat, T_old, dcvdT, dcpdT, &
#endif
                                       ydot) bind(C, name='fortran_to_cxx_actual_rhs')

    use actual_rhs_module, only: actual_rhs

    implicit none

    real(rt), intent(in), value :: rho, T, e, abar, zbar, y_e, eta
    real(rt), intent(in) :: xn(nspec)
#if NAUX_NET > 0
    real(rt), intent(in) :: aux(naux)
#endif
#ifndef SIMPLIFIED_SDC
    real(rt), intent(in), value :: cv, cp, T_old, dcvdT, dcpdT
    integer, intent(in), value :: self_heat
#endif
    real(rt), intent(inout) :: ydot(neqs)

    type (burn_t) :: state

    state % rho = rho
    state % T = T
    state % e = e
    state % abar = abar
    state % zbar = zbar
    state % y_e = y_e
    state % eta = eta
    state % xn(:) = xn(:)
#if NAUX_NET > 0
    state % aux(:) = aux(:)
#endif

#ifndef SIMPLIFIED_SDC
    state % cv = cv
    state % cp = cp

    if (self_heat == 1) then
       state % self_heat = .true.
    else
       state % self_heat = .false.
    end if

    state % T_old = T_old
    state % dcvdT = dcvdT
    state % dcpdT = dcpdT
#endif

    call actual_rhs(state, ydot)

  end subroutine fortran_to_cxx_actual_rhs



  subroutine fortran_to_cxx_actual_jac(rho, T, e, xn, abar, zbar, y_e, eta, &
#if NAUX_NET > 0
                                       aux, &
#endif
#ifndef SIMPLIFIED_SDC
                                       cv, cp, self_heat, T_old, dcvdT, dcpdT, &
#endif
                                       jac) bind(C, name='fortran_to_cxx_actual_jac')

    use actual_rhs_module, only: actual_jac

    implicit none

    real(rt), intent(in), value :: rho, T, e, abar, zbar, y_e, eta
    real(rt), intent(in) :: xn(nspec)
#if NAUX_NET > 0
    real(rt), intent(in) :: aux(naux)
#endif
#ifndef SIMPLIFIED_SDC
    real(rt), intent(in), value :: cv, cp, T_old, dcvdT, dcpdT
    integer, intent(in), value :: self_heat
#endif
    real(rt), intent(inout) :: jac(neqs, neqs)

    type (burn_t) :: state

    state % rho = rho
    state % T = T
    state % e = e
    state % abar = abar
    state % zbar = zbar
    state % y_e = y_e
    state % eta = eta
    state % xn(:) = xn(:)
#if NAUX_NET > 0
    state % aux(:) = aux(:)
#endif

#ifndef SIMPLIFIED_SDC
    state % cv = cv
    state % cp = cp

    if (self_heat == 1) then
       state % self_heat = .true.
    else
       state % self_heat = .false.
    end if

    state % T_old = T_old
    state % dcvdT = dcvdT
    state % dcpdT = dcpdT
#endif

    call actual_jac(state, jac)

  end subroutine fortran_to_cxx_actual_jac

end module fortran_to_cxx_actual_rhs_module
