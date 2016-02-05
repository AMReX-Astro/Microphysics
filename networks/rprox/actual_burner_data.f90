module actual_burner_data

  implicit none

  integer, parameter :: nrates    = 18

  integer, parameter :: irlambCNO = 1
  integer, parameter :: irag15o   = 2
  integer, parameter :: irr1      = 3
  integer, parameter :: irag16o   = 4
  integer, parameter :: irpg16o   = 5
  integer, parameter :: irpg17f   = 6
  integer, parameter :: irgp17f   = 7
  integer, parameter :: irlambda2 = 8
  integer, parameter :: irap14o   = 9
  integer, parameter :: irs1      = 10
  integer, parameter :: irlambda1 = 11
  integer, parameter :: ir3a      = 12
  integer, parameter :: irpg12c   = 13
  integer, parameter :: irwk14o   = 14
  integer, parameter :: irwk17f   = 15
  integer, parameter :: irwk15o   = 16
  integer, parameter :: irLweak   = 17
  integer, parameter :: irla2     = 18

  character (len=10), save :: reac_names(nrates)

  integer, parameter :: dlambCNOdh1   = 1
  integer, parameter :: drs1dhe4      = 2
  integer, parameter :: drr1dh1       = 3
  integer, parameter :: dlambda1dhe4  = 4
  integer, parameter :: dlambda2dhe4  = 5
  integer, parameter :: delta1        = 6
  integer, parameter :: delta2        = 7
  integer, parameter :: r56eff        = 8
  integer, parameter :: dr56effdt     = 9

end module actual_burner_data
