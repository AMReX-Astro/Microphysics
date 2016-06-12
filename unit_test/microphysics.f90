module microphysics_module

  use network
  use eos_module, only : eos_init
  use actual_rhs_module, only : actual_rhs_init
  use actual_burner_module, only : actual_burner_init

  
  implicit none

contains

  subroutine microphysics_init()

    call eos_init()

    call network_init()
    call actual_rhs_init()
    call actual_burner_init()
    
  end subroutine microphysics_init

  subroutine microphysics_finalize()

    !call network_finalize()

  end subroutine microphysics_finalize

end module microphysics_module
    
