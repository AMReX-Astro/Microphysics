! Aaron Jackson 2009
!
! This subroutine checks for local ignition conditions.
! If found, we calculate the desired detonation point from these
! ignition conditions.  We compare the new detonation point with previous
! detonation points and make sure we are not detonating too close to a
! previous one.
!
! If ignition_test is .true. then, detX, detY, detZ contain the detonation
! coords.
! If ignition_test is .false. then they do not matter.

subroutine bn_paraSpark(x, y, z, dens, pres, phfa, flame, c12, ne22, &
                        detX, detY, detZ, ignition_test)

  use Burn_data, ONLY : pbIgnRho, pbIgnRhoFact, pbIgnPhfa, pbIgnNum, &
                        pbIgnX, pbIgnY, pbIgnZ, pbIgnSep, pbIgnDist

  implicit none

#include "constants.h"  
#include "Flash.h"
  
  real, intent(in)     :: x, y, z, dens, pres, phfa, flame, c12, ne22
  real, intent(out)    :: detX, detY, detZ
  logical, intent(out) :: ignition_test

  real, parameter :: ye12 = 0.5
  real, parameter :: ye16 = 0.5
  real, parameter :: ye22 = 10.0/22.0

  integer :: i
  real :: r, costheta, sintheta, detD, test_dens, yei

  ignition_test = .false.

  detX = 0.0e0
  detY = 0.0e0
  detZ = 0.0e0

!  if ( flame >= pbIgnPhfa .and. flame < 0.99 ) then

  ! only if we are at the flame front edge 
  if ( flame > 0.001 .and. flame < 0.01 ) then

! **********************************************************************
!    This density is an estimate of the unburned density from
! Hansen & Kawaler.  This estimate uses the material composition and the
! pressure. This correspondence is good for dens >~ 1.5e5 (degenerate)
! and 0.92 is adjusted for fit.
!    Since we are in the flame, and our ddt-density is in the
! fuel, this estimate produces more consistant results with placing
! the detonation points by hand.
! **********************************************************************
!     yei = c12*ye12 + ne22*ye22 + (1.0-c12-ne22)*ye16
!     yei = 1.0e0/yei

!     test_dens = 0.92*yei*sqrt( (pres/1.243e15)**(6.0/4.0) + &
!                                (pres/1.004e13)**(6.0/5.0) )

! **********************************************************************
! This density we test for is the local density in the flame.
! However, as material is burned, it becomes less dense and we could meet
! these conditions sooner than we want.
! **********************************************************************
     test_dens = dens

     ! check that rhofact > dens >= rho
     ! rhofact needs to be less than 1
     if ( test_dens > pbIgnRhoFact .and.  &
          test_dens <= pbIgnRho ) then

        ignition_test = .true.

!**** This will need to be changed to handle other geometries!
        ! calculate the detonation coordinates
        r = sqrt(x**2 + y**2 + z**2)
        costheta = y/r
        sintheta = x/r
        r = r + pbIgnDist
        detX = r*sintheta
        detY = r*costheta
        detZ = 0.0e0

        ! now check that we are not igniting near something we've
        ! already detonated
        if (pbIgnNum > 0) then

           do i = 1, pbIgnNum

              ! distance to previous det point i
              detD = sqrt( ( detX - pbIgnX(i) )**2 + &
                           ( detY - pbIgnY(i) )**2 + &
                           ( detZ - pbIgnZ(i) )**2 )

              ! if we're too close, then we fail spark test
              if ( detD <= pbIgnSep ) then

                 ignition_test = .false.
                 detX = 0.e0
                 detY = 0.e0
                 detZ = 0.e0

                 ! this detonation is too close to a previous one
                 ! fail ignition_test and exit loop.
                 exit

              endif ! end separation check

           enddo

        endif ! end previous detonation number check

     endif ! end density check

  endif ! end in-flame check

  return
  
end subroutine bn_paraSpark
