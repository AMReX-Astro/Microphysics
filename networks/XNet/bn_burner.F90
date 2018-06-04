!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/XNet/bn_burner
!!
!! NAME
!!  
!!  bn_burner
!!
!! SYNOPSIS
!! 
!!  bn_burner(
!!       real(in)  :: tstep,
!!       real(in)  :: temp,
!!       real(in)  :: density,
!!       real(in)  :: xIn(:),
!!       real(out) :: xOut(:),
!!       real(out) :: sdotRate)
!!
!!  
!! DESCRIPTION
!!
!!  Routine bn_burner drives the nuclear burning network  
!!     given time step tstep, temperature temp, density density, and 
!!     composition xIn, this routine returns the burned composition xOut
!!     and the energy generation rate sdotRate.
!!
!! ARGUMENTS
!!
!!  tstep:    time step 
!!  temp:     temperature
!!  density:  density
!!  xIn:      composition in
!!  xOut:     composition out
!!  sdotRate: energy generation rate
!!  burnedZone: mask for which zones to burn
!!  kstep:    maximum number of burning timesteps
!!
!! NOTES
!!
!!  Within the network setup process
!!
!!  In this nuclearBurn directory, there are additional routines
!!  routine bn_azbar computes composition variables from the mass fractions; they are
!!              stored in Burn_dataEOS
!!  routine bn_sneutx computes neutrino losses
!!
!!
!!***

subroutine bn_burner(tstep,temp,density,xIn,xOut,sdotRate,burnedZone,kstep)

  use abundances, ONLY : y, ystart
  use conditions, ONLY : tdel, t, t9, rho, ye
  use controls, ONLY : szbatch, nzbatch, lzactive
  use thermo_data, ONLY : tstart, tstop, tdelstart, t9start, rhostart, &
    yestart, th, t9h, rhoh, yeh, nh
  use timers, ONLY : timer_burner, xnet_wtime
  use xnet_constants, ONLY : avn, epmev
  use xnet_data, ONLY : zz, aa, be
  use xnet_interface, ONLY : full_net

  implicit none

  ! arguments
  real, intent(in)                                      :: tstep
  real, intent(in), dimension(:)                        :: temp, density
  real, intent(in), dimension(:,:)                      :: xIn
  real, intent(out), dimension(size(xIn,1),size(xIn,2)) :: xOut
  real, intent(out), dimension(size(temp))              :: sdotRate
  logical, intent(in), dimension(:)                     :: burnedZone
  integer, intent(out)                                  :: kstep

  !..local varaibles      
  real, parameter ::  conv = avn*epmev
  integer :: i, numzones

  timer_burner = timer_burner - xnet_wtime()

  ! Active zone mask
  numzones = size(temp)
  szbatch = 1
  nzbatch = numzones
  lzactive = burnedZone

  ! Set the thermo history data for a constant conditions burn
  ! TODO: Call Microphysics abar/zbar routines
  do i = 1, numzones
  !   xmass = xIn(:,i)
  !   call bn_azbar()
  !   yestart(i) = bye
     ystart(:,i) = xIn(:,i) / aa(:)
  end do
  tstart    = 0.0
  tstop     = tstep
  tdelstart = 0.0
  t9start   = temp*1.0e-9
  rhostart  = density
  th(1,:)   = tstart  ; th(2,:)   = tstop
  t9h(1,:)  = t9start ; t9h(2,:)  = t9start
  rhoh(1,:) = density ; rhoh(2,:) = density
  yeh(1,:)  = yestart ; yeh(2,:)  = yestart
  nh        = 2

  ! Load initial abundances, time and timestep
  tdel = 0.0
  t    = tstart
  y    = ystart
  t9   = t9start
  rho  = rhostart
  ye   = yestart

  ! Evolve abundance from tstart to tstop
  if (any(burnedZone)) then
     call full_net(kstep)
  else
     kstep = 0
  end if

  do i = 1, numzones
     if (burnedZone(i)) then

        !..the average energy generated over the time step
        sdotRate(i) = sum((y(:,i) - ystart(:,i)) * be(:)) * conv / tstep

        !TODO: Call Microphysics neutrino routines
        !..take into account neutrino losses
        !btemp = t9(i)
        !bden  = rho(i)
        !xmass = xIn(:,i)
        !call bn_azbar()
        !call bn_sneutx()
        !sdotRate(i) = sdotRate(i) - sneut

        !..update the composition
        xOut(:,i) = y(:,i) * aa(:)
     else
        sdotRate(i) = 0.0e0
        xOut(:,i) = xIn(:,i)
     end if
  end do

  timer_burner = timer_burner + xnet_wtime()

  return
end subroutine bn_burner
