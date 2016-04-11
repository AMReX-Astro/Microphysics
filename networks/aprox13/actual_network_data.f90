module actual_network_data

  implicit none

  integer, parameter :: nspec = 13
  integer, parameter :: nspec_evolve = 13
  integer, parameter :: naux  = 0

  integer, parameter :: ihe4  = 1
  integer, parameter :: ic12  = 2
  integer, parameter :: io16  = 3
  integer, parameter :: ine20 = 4
  integer, parameter :: img24 = 5
  integer, parameter :: isi28 = 6
  integer, parameter :: is32  = 7
  integer, parameter :: iar36 = 8
  integer, parameter :: ica40 = 9
  integer, parameter :: iti44 = 10
  integer, parameter :: icr48 = 11
  integer, parameter :: ife52 = 12
  integer, parameter :: ini56 = 13

  double precision, save :: aion(nspec), zion(nspec), nion(nspec)
  double precision, save :: bion(nspec), mion(nspec), wion(nspec)

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  character (len=32), save :: network_name = "aprox13"

end module actual_network_data
