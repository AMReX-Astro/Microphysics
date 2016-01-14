module vode_data

  use network, only: nspec
  use bl_constants_module, only: ONE

  implicit none

  ! Indices of the temperature and energy variables in the work array

  integer, parameter :: net_itemp = nspec + 1
  integer, parameter :: net_ienuc = nspec + 2

  ! We may want to scale the temperature or density during the integration.

  double precision, save :: dens_scale = ONE
  double precision, save :: temp_scale = ONE  
  
end module vode_data
