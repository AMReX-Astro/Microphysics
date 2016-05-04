module actual_burner_data

  implicit none

  integer, parameter :: nrates = 2

  integer, parameter :: ir3a_   = 1
  integer, parameter :: ircago_ = 2

  character (len=10), save :: reac_names(nrates)

contains

  subroutine ener_gener_rate(dydt, enuc)

    use network

    implicit none

    double precision :: dydt(nspec), enuc

    enuc = sum(dydt(:) * ebin(:))

  end subroutine ener_gener_rate

end module actual_burner_data
