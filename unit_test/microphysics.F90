module microphysics_module

  use BoxLib
  use backtrace_module, only : set_fpe_trap
  use network
  use eos_module, only : eos_init, eos_finalize
  use actual_rhs_module, only : actual_rhs_init
  use managed_probin_module, only: managed_probin_init, managed_probin_finalize
#ifndef SDC
  use actual_burner_module, only : actual_burner_init
#endif

  implicit none

contains

  subroutine microphysics_init(small_temp, small_dens)

    implicit none
    
    double precision, optional :: small_temp
    double precision, optional :: small_dens


    !call boxlib_initialize()

    !call set_fpe_trap(.true., .true., .true.)

    call managed_probin_init()
    
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
#ifndef SDC
    call actual_burner_init()
#endif

  end subroutine microphysics_init

  subroutine microphysics_finalize()

    implicit none

    call eos_finalize()
    call network_finalize()
    call managed_probin_finalize()

  end subroutine microphysics_finalize

end module microphysics_module
