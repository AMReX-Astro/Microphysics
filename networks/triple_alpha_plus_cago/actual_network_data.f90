module actual_network_data

  implicit none

  integer, parameter :: nspec = 4
  integer, parameter :: nspec_evolve = 4
  integer, parameter :: naux  = 0

  integer, parameter :: ihe4_  = 1
  integer, parameter :: ic12_  = 2
  integer, parameter :: io16_  = 3
  integer, parameter :: ife56_ = 4

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  double precision, save :: aion(nspec), zion(nspec), ebin(nspec)

  character (len=22), parameter :: network_name = "triple_alpha_plus_cago"

end module actual_network_data
