! Common variables and routines for burners
! that use XNet for their integration.

module actual_integrator_module

   use eos_type_module
   use eos_module
   use network
   use burn_type_module
   use bl_types
   use bl_constants_module

   implicit none

contains

   subroutine actual_integrator(state_in, state_out, dt, time)

      !$acc routine seq

      use abundances, ONLY : y, ystart
      use conditions, ONLY : tdel, t, t9, rho, ye
      use controls, ONLY : szbatch, nzbatch, lzactive
      use thermo_data, ONLY : tstart, tstop, tdelstart, t9start, rhostart, &
         yestart, th, t9h, rhoh, yeh, nh
      use timers, ONLY : timer_burner, xnet_wtime
      use xnet_constants, ONLY : avn, epmev
      use xnet_interface, ONLY : full_net

      implicit none

      ! Input arguments

      type (burn_t), intent(in   ) :: state_in
      type (burn_t), intent(inout) :: state_out
      real(dp_t),    intent(in   ) :: dt, time

      ! Local variables

      real(dp_t), parameter ::  conv = avn*epmev
      integer :: i, numzones, kstep

      timer_burner = timer_burner - xnet_wtime()

      ! Active zone mask
      numzones = 1 ! Microphysics only passes 1 zone to burner
      szbatch = 1
      nzbatch = numzones
      lzactive = .false.
      lzactive(1) = .true.

      ! Set the thermo history data for a constant conditions burn
      do i = 1, numzones
         yestart(i) = state_in % y_e
         ystart(:,i) = state_in % xn(:) * aion_inv(:)
      end do
      tstart    = 0.0
      tstop     = dt
      tdelstart = 0.0
      t9start   = state_in % T * 1.0e-9
      rhostart  = state_in % rho
      th(1,:)   = tstart  ; th(2,:)   = tstop
      t9h(1,:)  = t9start ; t9h(2,:)  = t9start
      rhoh(1,:) = rhostart; rhoh(2,:) = rhostart
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
      if (any(lzactive)) then
         call full_net(kstep)
      else
         kstep = 0
      end if

      do i = 1, numzones
         if (lzactive(i)) then

            ! Update the state
            state_out % rho = rho
            state_out % T = t9
            state_out % y_e = ye
            state_out % xn(:) = y(:,i) * aion(:)

            state_out % time = time + dt

            state_out % e = state_in % e + sum((y(:,i) - ystart(:,i)) * bion(:)) * conv

            !TODO: Call Microphysics neutrino routines
            !..take into account neutrino losses
            !call bn_sneutx()
            !sdotRate(i) = sdotRate(i) - sneut
         else
            state_out % e = state_in % e
            state_out % xn(:) = state_in % xn(:)
         end if
      end do

      timer_burner = timer_burner + xnet_wtime()

      return
   end subroutine actual_integrator

end module actual_integrator_module
