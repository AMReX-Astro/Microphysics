! the burner_aux module contains thermodynamic state variables that are
! needed in the RHS and Jacobian routines

module burner_aux_module

  use amrex_fort_module, only : rt => amrex_real
  use network, only : nspec

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt), save :: sdc_rhoX_pass(nspec)
  real(rt), save :: sdc_rhoh_pass
  real(rt), save :: p0_pass  

end module burner_aux_module
