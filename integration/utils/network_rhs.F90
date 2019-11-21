module network_rhs_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public

contains

  subroutine network_rhs(state, reference_time)

    use actual_rhs_module, only: actual_rhs
    use burn_type_module, only: burn_t

#ifdef NONAKA_PLOT
    use nonaka_plot_module, only: nonaka_rhs
#endif

    implicit none

    type(burn_t), intent(inout) :: state
    real(rt),     intent(in)    :: reference_time

    !$gpu

    call actual_rhs(state)

#ifdef NONAKA_PLOT
    call nonaka_rhs(state, reference_time)
#endif

  end subroutine network_rhs


  subroutine network_jac(state, jac, reference_time)

    use actual_rhs_module, only: actual_jac
    use burn_type_module, only: burn_t, neqs

    implicit none

    type(burn_t), intent(inout) :: state
    real(rt) :: jac(neqs, neqs)
    real(rt),     intent(in)    :: reference_time

    !$gpu

    call actual_jac(state, jac)

  end subroutine network_jac

end module network_rhs_module
