module vode_indices

  use network, only: nspec

  implicit none

  ! Indices of the temperature and energy variables in the work array

  integer, parameter :: net_itemp = nspec + 1
  integer, parameter :: net_ienuc = nspec + 2

end module vode_indices
