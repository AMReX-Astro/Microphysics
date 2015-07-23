! the burner_aux module contains thermodynamic state variables that are
! needed in the RHS and Jacobian routines

module burner_aux_module

  use bl_types
  use network, only : nspec

  implicit none

  double precision :: dens_pass
  double precision :: c_p_pass
  double precision :: dhdx_pass(nspec)
  double precision :: X_O16_pass

  !$OMP THREADPRIVATE(dens_pass, c_p_pass, dhdx_pass, X_O16_pass)

end module burner_aux_module
