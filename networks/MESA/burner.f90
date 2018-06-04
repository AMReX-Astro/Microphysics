module burner_module

   use amrex_constants_module
   use amrex_error_module
   use amrex_fort_module, only : rt => amrex_real
   use eos_module
   use network

contains

   subroutine burner(dens, temp, Xin, dt, Xout, rho_omegadot, rho_Hnuc)

      implicit none
   
      ! INPUT:
      real(rt), intent(in) :: dens, temp, dt
      real(rt), intent(in) :: Xin(nspec)
   
      ! OUTPUT:
      real(rt), intent(out) :: Xout(nspec)
      real(rt), intent(out) :: rho_omegadot(nspec), rho_Hnuc

      ! LOCAL:
      integer :: i
      real(rt) :: dX, burn_ergs
      logical, save :: firstCall = .true.

      if (firstCall) then

         if (.not. network_initialized) then
            call amrex_error("ERROR in burner: must initialize network first")
         endif

         firstCall = .false.
      endif

      ! call the MESA network
      call Do_One_Zone_Burn(dens, temp, dt, Xin, burn_ergs, Xout)

      ! calculate rho_omegadot
      rho_omegadot(:) = 0.e0_rt
      do i=1,nspec
         dX = Xout(i) - Xin(i)
         rho_omegadot(i) = dens * dX / dt
      enddo

      ! calculate rho_Hnuc (burn_ergs has units of ergs/g)
      rho_Hnuc = dens * burn_ergs / dt

      return

   end subroutine burner

end module burner_module
