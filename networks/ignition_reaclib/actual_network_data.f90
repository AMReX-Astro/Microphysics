module actual_network_data

  implicit none

  ! Evolution and auxiliary
  integer, parameter :: nspec_evolve = 1
  integer, parameter :: naux  = 0

  ! Number of nuclear species in the network
  integer, parameter :: nspec = 7

  ! Number of entries in reactvec returned by rate_evaluate
  integer, parameter :: nreactvec = 6

  ! Binding Energies Per Nucleon (MeV)
  double precision :: ebind_per_nucleon(nspec)

  ! Nucleon mass number A
  double precision :: aion(nspec)

  ! Nucleon atomic number Z
  double precision :: zion(nspec)

  ! Nucleon neutron number N
  double precision :: nion(nspec)
  
  ! Binding Energies (ergs)
  double precision :: bion(nspec)

  ! Nuclides
  integer, parameter :: jn   = 1
  integer, parameter :: jp   = 2
  integer, parameter :: jhe4   = 3
  integer, parameter :: jc12   = 4
  integer, parameter :: jne20   = 5
  integer, parameter :: jna23   = 6
  integer, parameter :: jmg23   = 7
  ! Energy Generation Rate
  integer, parameter :: jenuc   = 8

  ! Reactions
  integer, parameter :: k_c12_c12a_ne20   = 1
  integer, parameter :: k_c12_c12n_mg23   = 2
  integer, parameter :: k_c12_c12p_na23   = 3
  integer, parameter :: k_n_p   = 4

  ! reactvec indices
  integer, parameter :: i_rate        = 1
  integer, parameter :: i_drate_dt    = 2
  integer, parameter :: i_scor        = 3
  integer, parameter :: i_dscor_dt    = 4
  integer, parameter :: i_dqweak      = 5
  integer, parameter :: i_epart       = 6

  character (len=16), save :: spec_names(nspec) 
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  double precision :: aion(nspec), zion(nspec), bion(nspec)
  double precision :: nion(nspec), mion(nspec), wion(nspec)

  !$acc declare create(aion, zion, bion, nion, mion, wion)

end module actual_network_data
