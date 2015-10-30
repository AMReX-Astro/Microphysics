!
! routine to shutdown the MESA net 
!

subroutine shutdown_mesa_net()

   use network, only: chem_id, net_iso
   use net_utils, only: which_rates

   implicit none

   deallocate(chem_id, net_iso, which_rates)


end subroutine shutdown_mesa_net

