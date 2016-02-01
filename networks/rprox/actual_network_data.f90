module actual_network_data

  implicit none

  integer, parameter :: nspec = 10
  integer, parameter :: naux  = 0

  integer, parameter :: ic12  = 1
  integer, parameter :: io14  = 2
  integer, parameter :: io15  = 3
  integer, parameter :: io16  = 4
  integer, parameter :: if17  = 5
  integer, parameter :: img22 = 6
  integer, parameter :: is30  = 7
  integer, parameter :: ini56 = 8
  integer, parameter :: ihe4  = 9
  integer, parameter :: ih1   = 10

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  double precision, save :: aion(nspec), zion(nspec), ebin(nspec)

  character (len=32), parameter :: network_name = "rprox"

end module actual_network_data
