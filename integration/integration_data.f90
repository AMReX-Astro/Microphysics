module integration_data

  use network, only: nspec
  use bl_constants_module, only: ONE

  implicit none

  ! We may want to scale certain quantities during the integration.

  double precision, save :: dens_scale = ONE
  double precision, save :: temp_scale = ONE
  double precision, save :: ener_scale = ONE

  !$acc declare create(dens_scale, temp_scale, ener_scale)

end module integration_data
