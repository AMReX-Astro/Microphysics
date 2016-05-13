module actual_burner_data

  implicit none

  integer, parameter :: nrates = 1

contains

  subroutine ener_gener_rate(dydt, ebin, enuc)

    use network

    implicit none

    double precision :: dydt(nspec), ebin(nspec), enuc

    enuc = sum(dydt(:) * aion(:) * ebin(:))

  end subroutine ener_gener_rate

end module actual_burner_data
