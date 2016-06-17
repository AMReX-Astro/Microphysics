module integration_data

  use network, only: nspec
  use bl_constants_module, only: ONE
  use bl_types, only: dp_t

  implicit none

  ! We may want to scale certain quantities during the integration.

  real(dp_t), save :: dens_scale = ONE
  real(dp_t), save :: temp_scale = ONE
  real(dp_t), save :: ener_scale = ONE

  real(dp_t), save :: inv_dens_scale = ONE
  real(dp_t), save :: inv_temp_scale = ONE
  real(dp_t), save :: inv_ener_scale = ONE

  !$acc declare create(dens_scale, temp_scale, ener_scale)
  !$acc declare create(inv_dens_scale, inv_temp_scale, inv_ener_scale)

  ! Multiplicative inverse of atomic numbers.

  real(dp_t), save :: aionInv(nspec)

  !$acc declare create(aionInv)

end module integration_data
