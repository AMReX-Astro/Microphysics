module actual_network_data

  implicit none

  integer, parameter :: nspec = 19
  integer, parameter :: nspec_evolve = 19
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
  integer, parameter :: ife52 = 15
  integer, parameter :: ife54 = 16
  integer, parameter :: ini56 = 17
  integer, parameter :: ineut = 18
  integer, parameter :: iprot = 19

  double precision, save :: aion(nspec), zion(nspec), nion(nspec)
  double precision, save :: bion(nspec), mion(nspec), wion(nspec)

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  character (len=32), save :: network_name = "aprox19"

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

  ! Rates data

  integer, parameter :: nrates = 100

  integer, parameter :: ir3a   = 1
  integer, parameter :: irg3a  = 2
  integer, parameter :: ircag  = 3
  integer, parameter :: ir1212 = 4
  integer, parameter :: ir1216 = 5
  integer, parameter :: ir1616 = 6
  integer, parameter :: iroga  = 7
  integer, parameter :: iroag  = 8
  integer, parameter :: irnega = 9
  integer, parameter :: irneag = 10
  integer, parameter :: irmgga = 11
  integer, parameter :: irmgag = 12
  integer, parameter :: irsiga = 13
  integer, parameter :: irmgap = 14
  integer, parameter :: iralpa = 15
  integer, parameter :: iralpg = 16
  integer, parameter :: irsigp = 17
  integer, parameter :: irsiag = 18
  integer, parameter :: irsga  = 19
  integer, parameter :: irsiap = 20
  integer, parameter :: irppa  = 21
  integer, parameter :: irppg  = 22
  integer, parameter :: irsgp  = 23
  integer, parameter :: irsag  = 24
  integer, parameter :: irarga = 25
  integer, parameter :: irsap  = 26
  integer, parameter :: irclpa = 27
  integer, parameter :: irclpg = 28
  integer, parameter :: irargp = 29
  integer, parameter :: irarag = 30
  integer, parameter :: ircaga = 31
  integer, parameter :: irarap = 32
  integer, parameter :: irkpa  = 33
  integer, parameter :: irkpg  = 34
  integer, parameter :: ircagp = 35
  integer, parameter :: ircaag = 36
  integer, parameter :: irtiga = 37
  integer, parameter :: ircaap = 38
  integer, parameter :: irscpa = 39
  integer, parameter :: irscpg = 40
  integer, parameter :: irtigp = 41
  integer, parameter :: irtiag = 42
  integer, parameter :: ircrga = 43
  integer, parameter :: irtiap = 44
  integer, parameter :: irvpa  = 45
  integer, parameter :: irvpg  = 46
  integer, parameter :: ircrgp = 47
  integer, parameter :: ircrag = 48
  integer, parameter :: irfega = 49
  integer, parameter :: ircrap = 50
  integer, parameter :: irmnpa = 51
  integer, parameter :: irmnpg = 52
  integer, parameter :: irfegp = 53
  integer, parameter :: irfeag = 54
  integer, parameter :: irniga = 55
  integer, parameter :: irfeap = 56
  integer, parameter :: ircopa = 57
  integer, parameter :: ircopg = 58
  integer, parameter :: irnigp = 59
  integer, parameter :: ir52ng = 60
  integer, parameter :: ir53gn = 61
  integer, parameter :: ir53ng = 62
  integer, parameter :: ir54gn = 63
  integer, parameter :: irfepg = 64
  integer, parameter :: ircogp = 65
  integer, parameter :: irheng = 66
  integer, parameter :: irhegn = 67
  integer, parameter :: irhng  = 68
  integer, parameter :: irdgn  = 69
  integer, parameter :: irdpg  = 70
  integer, parameter :: irhegp = 71
  integer, parameter :: irpen   = 72
  integer, parameter :: irnep   = 73
  integer, parameter :: irn56ec = 74
  integer, parameter :: irpp    = 75
  integer, parameter :: ir33    = 76
  integer, parameter :: irhe3ag = 77
  integer, parameter :: ircpg  = 78
  integer, parameter :: irnpg  = 79
  integer, parameter :: ifa    = 80
  integer, parameter :: ifg    = 81
  integer, parameter :: iropg  = 82
  integer, parameter :: irnag  = 83
  integer, parameter :: irr1   = 84
  integer, parameter :: irs1   = 85
  integer, parameter :: irt1   = 86
  integer, parameter :: iru1   = 87
  integer, parameter :: irv1   = 88
  integer, parameter :: irw1   = 89
  integer, parameter :: irx1   = 90
  integer, parameter :: ir1f54 = 91
  integer, parameter :: ir2f54 = 92
  integer, parameter :: ir3f54 = 93
  integer, parameter :: ir4f54 = 94
  integer, parameter :: ir5f54 = 95
  integer, parameter :: ir6f54 = 96
  integer, parameter :: ir7f54 = 97
  integer, parameter :: ir8f54 = 98

  integer, parameter :: iralf1 = 99
  integer, parameter :: iralf2 = 100

  character (len=16), save :: ratenames(nrates)

end module actual_network_data
