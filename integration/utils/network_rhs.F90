module network_rhs_module

  implicit none

  public

contains

  subroutine network_rhs(state)

    use actual_rhs_module, only: actual_rhs
    use burn_type_module, only: burn_t

#ifdef NONAKA_PLOT
    use nonaka_plot_module, only: nonaka_rhs
#endif

    implicit none

    type(burn_t), intent(inout) :: state

    !$gpu

    call actual_rhs(state)

#ifdef NONAKA_PLOT
    call nonaka_rhs(state)
#endif

  end subroutine network_rhs


  subroutine network_jac(state)

    use actual_rhs_module, only: actual_jac
    use burn_type_module, only: burn_t

    implicit none

    type(burn_t), intent(inout) :: state

    !$gpu

    call actual_jac(state)

  end subroutine network_jac

end module network_rhs_module
