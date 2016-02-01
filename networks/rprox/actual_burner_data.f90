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

  integer, save :: irp_rates
  integer, save :: irp_drtdt
  integer, save :: irp_dlambCNOdh1
  integer, save :: irp_drs1dhe4
  integer, save :: irp_drr1dh1
  integer, save :: irp_dlambda1dhe4
  integer, save :: irp_dlambda2dhe4
  integer, save :: irp_delta1
  integer, save :: irp_delta2
  integer, save :: irp_r56eff
  integer, save :: irp_dr56effdt

end module actual_burner_data
