module actual_network_data

  implicit none

  integer, parameter :: nspec = 3
  integer, parameter :: naux  = 0

  integer, parameter :: ic12_  = 1
  integer, parameter :: io16_  = 2
  integer, parameter :: iash_  = 3

  character (len=16), save :: spec_names(nspec) 
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  double precision, save :: aion(nspec), zion(nspec), bion(nspec)
  double precision, save :: nion(nspec), mion(nspec), wion(nspec)

end module actual_network_data
