module actual_burner_data

  implicit none

  integer, parameter :: nrates = 2

  integer, parameter :: ir3a_   = 1
  integer, parameter :: ircago_ = 2

  character (len=10), save :: reac_names(nrates)

end module actual_burner_data
