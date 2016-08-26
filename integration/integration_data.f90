module integration_data

  use network, only: nspec
  use bl_types, only: dp_t

  implicit none

  ! Multiplicative inverse of atomic numbers.

  real(dp_t), save :: aionInv(nspec)

  !$acc declare create(aionInv)

end module integration_data
