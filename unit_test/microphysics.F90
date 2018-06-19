module microphysics_module

  use BoxLib
  use backtrace_module, only : set_fpe_trap
  use network
  use eos_module, only : eos_init
  use actual_rhs_module, only : actual_rhs_init
#ifndef SDC
  use actual_burner_module, only : actual_burner_init
#endif

  implicit none

contains

  subroutine microphysics_init(small_temp, small_dens)

    double precision, optional :: small_temp
    double precision, optional :: small_dens


    !call boxlib_initialize()

    !call set_fpe_trap(.true., .true., .true.)

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

    use eos_module, only: eos_finalize
#ifdef USE_SCREENING
    use screening_module, only: screening_finalize
    call screening_finalize()
#endif
    call eos_finalize()
    call network_finalize()

  end subroutine microphysics_finalize

end module microphysics_module
