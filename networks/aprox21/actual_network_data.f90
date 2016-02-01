module actual_network_data

  implicit none

  integer, parameter :: nspec = 21
  integer, parameter :: naux  = 0
  
  integer, parameter :: ih1   = 1
  integer, parameter :: ihe3  = 2
  integer, parameter :: ihe4  = 3
  integer, parameter :: ic12  = 4
  integer, parameter :: in14  = 5
  integer, parameter :: io16  = 6
  integer, parameter :: ine20 = 7
  integer, parameter :: img24 = 8
  integer, parameter :: isi28 = 9
  integer, parameter :: is32  = 10
  integer, parameter :: iar36 = 11
  integer, parameter :: ica40 = 12
  integer, parameter :: iti44 = 13
  integer, parameter :: icr48 = 14
  integer, parameter :: icr56 = 15
  integer, parameter :: ife52 = 16
  integer, parameter :: ife54 = 17
  integer, parameter :: ife56 = 18
  integer, parameter :: ini56 = 19
  integer, parameter :: ineut = 20
  integer, parameter :: iprot = 21

  double precision, save :: aion(nspec), zion(nspec), nion(nspec)
  double precision, save :: bion(nspec), mion(nspec), wion(nspec)

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  character (len=32), save :: network_name = "aprox21"

end module actual_network_data
