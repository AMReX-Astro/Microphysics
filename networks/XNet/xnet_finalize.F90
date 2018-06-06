subroutine xnet_finalize()

   use gpu_controls, ONLY : gpu_finalize
   use xnet_interface, ONLY : jacobian_finalize

   implicit none

   call jacobian_finalize
   call gpu_finalize

   return
end subroutine xnet_finalize
