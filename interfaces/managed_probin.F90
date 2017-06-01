module managed_probin_module

  use bl_types
  
  implicit none

  logical, managed, allocatable, save :: cu_do_constant_volume_burn
  logical, managed, allocatable, save :: cu_call_eos_in_rhs
  logical, managed, allocatable, save :: cu_centered_diff_jac  
  real(dp_t), managed, allocatable, save :: cu_dt_crit
  real(dp_t), managed, allocatable, save :: cu_tmax

  !$acc declare &
  !$acc create(cu_do_constant_volume_burn, cu_call_eos_in_rhs, cu_dt_crit) &
  !$acc create(cu_centered_diff_jac, cu_tmax)
  
contains

  subroutine managed_probin_init()
    use extern_probin_module, only: do_constant_volume_burn, dT_crit, call_eos_in_rhs, &
         centered_diff_jac
    use probin_module, only: tmax

    implicit none

    ! Allocate and set managed memory probin parameters
    allocate(cu_do_constant_volume_burn)
    allocate(cu_call_eos_in_rhs)
    allocate(cu_centered_diff_jac)
    allocate(cu_dt_crit)
    allocate(cu_tmax)

    cu_do_constant_volume_burn = do_constant_volume_burn
    cu_call_eos_in_rhs = call_eos_in_rhs
    cu_centered_diff_jac = centered_diff_jac
    cu_dt_crit = dT_crit
    cu_tmax = tmax
  end subroutine managed_probin_init

  subroutine managed_probin_finalize()
    
    implicit none

    ! Deallocate managed memory probin parameters
    deallocate(cu_do_constant_volume_burn)
    deallocate(cu_call_eos_in_rhs)
    deallocate(cu_centered_diff_jac)
    deallocate(cu_dt_crit)
    deallocate(cu_tmax)
  end subroutine managed_probin_finalize

end module managed_probin_module
