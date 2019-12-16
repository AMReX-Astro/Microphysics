! This module contains the dydt routine which change in species with respect
! to time based on the reaction rates.

module dydt_module

  use amrex_constants_module
  use network

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine dydt(ymol, rates, ydot)

    !$acc routine seq

    use amrex_fort_module, only : rt => amrex_real
    real(rt)        , intent(IN   ) :: ymol(nspec), rates(nrates)
    real(rt)        , intent(  OUT) :: ydot(nspec_evolve)

    !$gpu

    ydot(:) = ZERO

    ! He4 reactions
    ydot(ihe4) = - THREE * ymol(ihe4)**THREE * rates(ir3a)              &
                 - ONE * ymol(ic12) * ymol(ihe4) * rates(ircago)

    ! C12 reactions
    ydot(ic12) =   ONE * ymol(ihe4)**THREE * rates(ir3a)                &
                 - ONE * ymol(ic12) * ymol(ihe4) * rates(ircago)

    ! O16 reactions
    ydot(io16) =  ONE * ymol(ic12) * ymol(ihe4) * rates(ircago)

  end subroutine dydt

end module dydt_module
