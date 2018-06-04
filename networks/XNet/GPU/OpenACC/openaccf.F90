module openaccf
   use, intrinsic :: iso_c_binding
   implicit none

   interface
      subroutine acc_map_data(hostptr,devptr,bytes) &
            bind(c,name="acc_map_data")
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: hostptr
         type(c_ptr), value :: devptr
         integer(c_size_t), value :: bytes
      end subroutine acc_map_data

      subroutine acc_unmap_data(hostptr) &
            bind(c,name="acc_unmap_data")
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: hostptr
      end subroutine acc_unmap_data

      type(c_ptr) function acc_deviceptr(hostptr) &
            bind(c,name="acc_deviceptr")
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: hostptr
      end function acc_deviceptr

      type(c_ptr) function acc_hostptr(devptr) &
            bind(c,name="acc_hostptr")
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: devptr
      end function acc_hostptr
   end interface

end module