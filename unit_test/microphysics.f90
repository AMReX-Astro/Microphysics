module microphysics_module

  use network
  use eos_module, only : eos_init
  use actual_rhs_module, only : actual_rhs_init
  use actual_burner_module, only : actual_burner_init

  
  implicit none

contains

  subroutine microphysics_init(small_temp, small_dens)

    double precision, optional :: small_temp
    double precision, optional :: small_dens


    if (present(small_temp) .and. present(small_dens)) then
       call eos_init(small_temp=small_temp, small_dens=small_dens)
    else if (present(small_temp)) then
       call eos_init(small_temp=small_temp)
    else if (present(small_dens)) then
       call eos_init(small_dens=small_dens)
    else
       call eos_init()
    endif

    call network_init()
    call actual_rhs_init()
    call actual_burner_init()
    
  end subroutine microphysics_init

  subroutine microphysics_finalize()

    !call network_finalize()

  end subroutine microphysics_finalize

end module microphysics_module
    
