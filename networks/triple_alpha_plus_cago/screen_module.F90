! This module contains the screen routine which applies screening corrections
! to the reaction rates.  The input composition must be in terms of molar
! fractions.  The triple-alpha reaction goes through a compound nucleus
! channel and therefore screening must be applied to both reactions.
!
! A call is made to screenz which lives in screen.f to apply the screening
! to a single reaction.
module screen_module

  use amrex_constants_module
  use amrex_fort_module, only : rt => amrex_real
  use network

  implicit none

contains

  subroutine screen(temp, dens, ymol, rates, dratesdt)

    !$acc routine seq

    use screening_module, only: screen5, plasma_state, fill_plasma_state

    real(rt), intent(IN   ) :: temp, dens, ymol(nspec)
    real(rt), intent(INOUT) :: rates(nrates)
    real(rt), intent(INOUT) :: dratesdt(nrates)

    real(rt) :: scorr1, dscorr1dt, dscorr1dd
    real(rt) :: scorr2, dscorr2dt, dscorr2dd
    real(rt) ::  scorr,  dscorrdt, dscorrdd

    real(rt) :: rates_in(nrates), dratesdt_in(nrates)
    integer :: jscr

    type (plasma_state) :: state

    !$gpu
    
    rates_in    = rates
    dratesdt_in = dratesdt

    call fill_plasma_state(state, temp, dens, ymol)

    jscr = 1
    call screen5(state,jscr,scorr1,dscorr1dt,dscorr1dd)

    jscr = jscr + 1
    call screen5(state,jscr,scorr2,dscorr2dt,dscorr2dd)

    scorr    = scorr1 * scorr2
    dscorrdt = dscorr1dt * scorr2 + scorr1 * dscorr2dt

    rates(ir3a)    = rates_in(ir3a) * scorr
    dratesdt(ir3a) = dratesdt_in(ir3a) * scorr +                        &
                     rates_in(ir3a) * dscorrdt

    ! C12 + alpha --> O16
    jscr = jscr + 1
    call screen5(state,jscr,scorr,dscorrdt,dscorrdd)

    rates(ircago)    = rates_in(ircago) * scorr
    dratesdt(ircago) = dratesdt_in(ircago) * scorr + &
                        rates_in(ircago) * dscorrdt

    return

  end subroutine screen

end module screen_module
