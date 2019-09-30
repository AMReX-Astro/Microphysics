module rhs_module

contains
  
  ! The f_rhs routine provides the right-hand-side for the DVODE solver.
  subroutine f_rhs(time, y, ydot, rpar)

    use sdc_parameters_module, only: SDC_NEQS
    use sdc_rpar_indices, only: n_rpar_comps
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(IN) :: time, y(SDC_NEQS)
    real(rt), intent(INOUT) :: rpar(n_rpar_comps)
    real(rt), intent(INOUT) :: ydot(SDC_NEQS)

    !$gpu

    ydot(1) = -0.04d0 * y(1) + 1.d4 * y(2) * y(3)
    ydot(2) = 0.04d0 * y(1) - 1.d4 * y(2) * y(3) - 3.e7 * y(2)**2
    ydot(3) = 3.e7 * y(2)**2

  end subroutine f_rhs


  ! Analytical Jacobian
  subroutine jac(time, y, ml, mu, pd, nrpd, rpar)

    use sdc_parameters_module, only: SDC_NEQS
    use sdc_rpar_indices, only: n_rpar_comps
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer, intent(IN   ) :: ml, mu, nrpd
    real(rt), intent(in) :: time
    real(rt), intent(INOUT) :: y(SDC_NEQS), rpar(n_rpar_comps)
    real(rt), intent(  OUT) :: pd(SDC_NEQS,SDC_NEQS)

    !$gpu

    pd(:,:) = 0.0d0

    pd(1,1) = -0.04d0
    pd(1,2) = 1.d4 * y(3)
    pd(1,3) = 1.d4 * y(2)
    pd(2,1) = 0.0d0
    pd(2,3) = -pd(1,3)
    pd(3,2) = 6.d7 * y(2)
    pd(2,2) = -pd(1,2) - pd(3,2)

  end subroutine jac

end module rhs_module
