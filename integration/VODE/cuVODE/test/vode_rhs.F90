module vode_rhs_module

contains
  
  ! The f_rhs routine provides the right-hand-side for the DVODE solver.
  subroutine f_rhs(time, dvode_state, ydot)

    use cuvode_parameters_module, only: VODE_NEQS
    use amrex_fort_module, only: rt => amrex_real
    use cuvode_types_module, only: dvode_t

    implicit none

    real(rt),      intent(in   ) :: time
    type(dvode_t), intent(in   ) :: dvode_state
    real(rt),      intent(inout) :: ydot(VODE_NEQS)

    !$gpu

    ydot(1) = -.04e0_rt * dvode_state % Y(1) + 1.e4_rt * dvode_state % Y(2) * dvode_state % Y(3)
    ydot(3) = 3.e7_rt * dvode_state % Y(2) * dvode_state % Y(2)
    ydot(2) = -ydot(1) - ydot(3)

  end subroutine f_rhs


  ! Analytical Jacobian
  subroutine jac(time, dvode_state, ml, mu, pd, nrpd)

    use cuvode_parameters_module, only: VODE_NEQS
    use cuvode_types_module, only: dvode_t
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer ,      intent(in   ) :: ml, mu, nrpd
    type(dvode_t), intent(in   ) :: dvode_state
    real(rt),      intent(in   ) :: time
    real(rt),      intent(inout) :: pd(VODE_NEQS,VODE_NEQS)

    !$gpu

    pd(1,1) = -.04e0_rt
    pd(1,2) = 1.e4_rt * dvode_state % Y(3)
    pd(1,3) = 1.e4_rt * dvode_state % Y(2)
    pd(2,1) = .04e0_rt
    pd(2,3) = -pd(1,3)
    pd(3,2) = 6.e7_rt * dvode_state % Y(2)
    pd(2,2) = -pd(1,2) - pd(3,2)

  end subroutine jac

end module vode_rhs_module
