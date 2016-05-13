module actual_network_data

  implicit none

  integer, parameter :: nspec = 3
  integer, parameter :: nspec_evolve = 2
  integer, parameter :: naux  = 0

  integer, parameter :: ifuel_  = 1
  integer, parameter :: iash_   = 2
  integer, parameter :: iinert_ = 3

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  double precision, save :: aion(nspec), zion(nspec), ebin(nspec)

end module actual_network_data
