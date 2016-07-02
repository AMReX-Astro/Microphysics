module vode_type_module

  implicit none

contains

  subroutine clean_state(y, rpar)

    use bl_types, only: dp_t
    use actual_network, only: aion, nspec, nspec_evolve
    use integration_data, only: aionInv
    use burn_type_module, only: neqs
    use rpar_indices, only: n_rpar_comps

    implicit none

    real(dp_t) :: y(neqs), rpar(n_rpar_comps)

    ! Ensure that mass fractions always stay positive.

    y(1:nspec_evolve) = max(y(1:nspec_evolve) * aion(1:nspec_evolve), 1.d-200) * aionInv(1:nspec_evolve)

  end subroutine clean_state



  subroutine renormalize_species(y, rpar)

    use bl_types, only: dp_t
    use actual_network, only: aion, nspec, nspec_evolve
    use integration_data, only: aionInv
    use burn_type_module, only: neqs
    use rpar_indices, only: n_rpar_comps, irp_nspec

    implicit none

    real(dp_t) :: y(neqs), rpar(n_rpar_comps)

    real(dp_t) :: nspec_sum

    nspec_sum = sum(y(1:nspec_evolve) * aion(1:nspec_evolve)) + &
                sum(rpar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) * aion(nspec_evolve+1:nspec))

    y(1:nspec_evolve) = y(1:nspec_evolve) / nspec_sum
    rpar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) = rpar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) / nspec_sum

  end subroutine renormalize_species

end module vode_type_module
