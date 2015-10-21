! the burner_aux module contains thermodynamic state variables that are
! needed in the RHS and Jacobian routines

module burner_aux_module

  use bl_types
  use network, only : nspec

  implicit none

  real(kind=dp_t), save :: sdc_rhoX_pass(nspec)
  real(kind=dp_t), save :: sdc_rhoh_pass
  real(kind=dp_t), save :: p0_pass  

end module burner_aux_module
