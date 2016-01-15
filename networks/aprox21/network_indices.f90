module network_indices

  implicit none

  ! ids of the species
  
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

  ! set the id numbers of the reaction rates

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

  ! for fe54 photodisintegration
  integer, parameter :: ir52ng = 60
  integer, parameter :: ir53gn = 61
  integer, parameter :: ir53ng = 62
  integer, parameter :: ir54gn = 63
  integer, parameter :: irfepg = 64
  integer, parameter :: ircogp = 65

  ! for he4 photodisintegration
  integer, parameter :: irheng = 66
  integer, parameter :: irhegn = 67
  integer, parameter :: irhng  = 68
  integer, parameter :: irdgn  = 69
  integer, parameter :: irdpg  = 70
  integer, parameter :: irhegp = 71

  ! weak reactions
  integer, parameter :: irpen   = 72
  integer, parameter :: irnep   = 73
  integer, parameter :: irn56ec = 74

  ! ppchain
  integer, parameter :: irpp    = 75
  integer, parameter :: ir33    = 76
  integer, parameter :: irhe3ag = 77

  ! cno cycles
  integer, parameter :: ircpg  = 78
  integer, parameter :: irnpg  = 79
  integer, parameter :: ifa    = 80
  integer, parameter :: ifg    = 81
  integer, parameter :: iropg  = 82
  integer, parameter :: irnag  = 83

  ! for reactions to fe56 
  integer, parameter :: ir54ng   = 84
  integer, parameter :: ir55gn   = 85
  integer, parameter :: ir55ng   = 86
  integer, parameter :: ir56gn   = 87
  integer, parameter :: irfe54ap = 88
  integer, parameter :: irco57pa = 89
  integer, parameter :: irfe56pg = 90
  integer, parameter :: irco57gp = 91

  ! the equilibrium links
  integer, parameter :: irr1   = 92
  integer, parameter :: irs1   = 93
  integer, parameter :: irt1   = 94
  integer, parameter :: iru1   = 95
  integer, parameter :: irv1   = 96
  integer, parameter :: irw1   = 97
  integer, parameter :: irx1   = 98

  integer, parameter :: ir1f54 = 99
  integer, parameter :: ir2f54 = 100
  integer, parameter :: ir3f54 = 101
  integer, parameter :: ir4f54 = 102
  integer, parameter :: ir5f54 = 103
  integer, parameter :: ir6f54 = 104
  integer, parameter :: ir7f54 = 105
  integer, parameter :: ir8f54 = 106

  integer, parameter :: iralf1 = 107
  integer, parameter :: iralf2 = 108

  integer, parameter :: irfe56_aux1 = 109
  integer, parameter :: irfe56_aux2 = 110
  integer, parameter :: irfe56_aux3 = 111
  integer, parameter :: irfe56_aux4 = 112

  ! rpar data

  integer :: irp_dydt
  integer :: irp_rates
  integer :: irp_drdy1
  integer :: irp_drdy2

end module network_indices
