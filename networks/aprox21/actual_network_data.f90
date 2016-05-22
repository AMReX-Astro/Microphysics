module actual_network_data

  implicit none

  integer, parameter :: nspec = 21
  integer, parameter :: nspec_evolve = 21
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

  ! Some fundamental physical constants

  double precision, parameter :: avo = 6.0221417930d23
  double precision, parameter :: c_light = 2.99792458d10

  double precision, parameter :: ev2erg  = 1.60217648740d-12
  double precision, parameter :: mev2erg = ev2erg*1.0d6
  double precision, parameter :: mev2gr  = mev2erg/c_light**2

  double precision, parameter :: mn = 1.67492721184d-24
  double precision, parameter :: mp = 1.67262163783d-24
  double precision, parameter :: me = 9.1093821545d-28

  ! Conversion factor for the nuclear energy generation rate.

  double precision, parameter :: enuc_conv2 = -avo*c_light*c_light

end module actual_network_data
