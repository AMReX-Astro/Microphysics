module integration_data

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  type :: integration_status_t

     real(rt) :: atol_spec, atol_enuc, atol_temp
     real(rt) :: rtol_spec, rtol_enuc, rtol_temp

  end type integration_status_t

end module integration_data
