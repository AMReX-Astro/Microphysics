module burner_module

   use bl_types
   use bl_constants_module
   use bl_error_module
   use eos_module
   use network

contains

   subroutine burner(dens, temp, Xin, dt, Xout, rho_omegadot, rho_Hnuc)

      implicit none
   
      ! INPUT:
      real(kind=dp_t), intent(in) :: dens, temp, dt
      real(kind=dp_t), intent(in) :: Xin(nspec)
   
      ! OUTPUT:
      real(kind=dp_t), intent(out) :: Xout(nspec)
      real(kind=dp_t), intent(out) :: rho_omegadot(nspec), rho_Hnuc

      ! LOCAL:
      integer :: i
      real(kind=dp_t) :: dX, burn_ergs
      logical, save :: firstCall = .true.

      if (firstCall) then

         if (.not. network_initialized) then
            call bl_error("ERROR in burner: must initialize network first")
         endif

         firstCall = .false.
      endif

      ! call the MESA network
      call Do_One_Zone_Burn(dens, temp, dt, Xin, burn_ergs, Xout)

      ! calculate rho_omegadot
      rho_omegadot(:) = 0.e0_dp_t
      do i=1,nspec
         dX = Xout(i) - Xin(i)
         rho_omegadot(i) = dens * dX / dt
      enddo

      ! calculate rho_Hnuc (burn_ergs has units of ergs/g)
      rho_Hnuc = dens * burn_ergs / dt

      return

   end subroutine burner

end module burner_module
