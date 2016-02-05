module vode_data

  use network, only: nspec
  use bl_constants_module, only: ONE

  implicit none

  ! We may want to scale the temperature or density during the integration.

  double precision, save :: dens_scale = ONE
  double precision, save :: temp_scale = ONE  
  
end module vode_data
