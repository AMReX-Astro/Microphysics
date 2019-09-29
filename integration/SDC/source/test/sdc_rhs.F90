module rhs_module

contains
  
  ! The f_rhs routine provides the right-hand-side for the DVODE solver.
  subroutine f_rhs(time, y, ydot, rpar)

    use cuvode_parameters_module, only: VODE_NEQS
    use vode_rpar_indices, only: n_rpar_comps
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(INOUT) :: time, y(VODE_NEQS)
    real(rt), intent(INOUT) :: rpar(n_rpar_comps)
    real(rt), intent(INOUT) :: ydot(VODE_NEQS)

    !$gpu

    YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
    YDOT(3) = 3.D7*Y(2)*Y(2)
    YDOT(2) = -YDOT(1) - YDOT(3)

  end subroutine f_rhs


  ! Analytical Jacobian
  subroutine jac(time, y, ml, mu, pd, nrpd, rpar)

    use cuvode_parameters_module, only: VODE_NEQS
    use vode_rpar_indices, only: n_rpar_comps
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer   , intent(IN   ) :: ml, mu, nrpd
    real(rt), intent(INOUT) :: y(VODE_NEQS), rpar(n_rpar_comps), time
    real(rt), intent(  OUT) :: pd(VODE_NEQS,VODE_NEQS)

    !$gpu

    PD(1,1) = -.04D0
    PD(1,2) = 1.D4*Y(3)
    PD(1,3) = 1.D4*Y(2)
    PD(2,1) = .04D0
    PD(2,3) = -PD(1,3)
    PD(3,2) = 6.D7*Y(2)
    PD(2,2) = -PD(1,2) - PD(3,2)

  end subroutine jac

end module rhs_module
