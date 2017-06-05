module managed_probin_module

  use bl_types
  
  implicit none
  
  logical, &
#ifdef CUDA       
       managed, &
#endif       
       allocatable, save :: cu_do_constant_volume_burn
  logical, &
#ifdef CUDA       
       managed, &
#endif       
       allocatable, save :: cu_call_eos_in_rhs
  logical, &
#ifdef CUDA       
       managed, &
#endif       
       allocatable, save :: cu_centered_diff_jac
  logical, &
#ifdef CUDA       
       managed, &
#endif              
       allocatable, save :: cu_renormalize_abundances
  logical, &
#ifdef CUDA       
       managed, &
#endif
       allocatable, save :: cu_integrate_temperature
  logical, &
#ifdef CUDA       
       managed, &
#endif       
       allocatable, save :: cu_integrate_energy
  
  integer, &
#ifdef CUDA       
       managed, &
#endif       
       allocatable, save :: cu_burning_mode
  
  real(dp_t), &
#ifdef CUDA       
       managed, &
#endif       
       allocatable, save :: cu_burning_mode_factor  
  real(dp_t), &
#ifdef CUDA       
       managed, &
#endif       
       allocatable, save :: cu_dt_crit
  real(dp_t), &
#ifdef CUDA       
       managed, &
#endif       
       allocatable, save :: cu_tmax

  !$acc declare &
  !$acc create(cu_do_constant_volume_burn, cu_call_eos_in_rhs, cu_dt_crit) &
  !$acc create(cu_centered_diff_jac, cu_tmax, cu_renormalize_abundances) &
  !$acc create(cu_burning_mode, cu_burning_mode_factor) &
  !$acc create(cu_integrate_temperature, cu_integrate_energy)
  
contains

  subroutine managed_probin_init()
    use extern_probin_module, only: do_constant_volume_burn, dT_crit, call_eos_in_rhs, &
         centered_diff_jac, renormalize_abundances, &
         burning_mode, burning_mode_factor, &
         integrate_temperature, integrate_energy

    use probin_module, only: tmax

    implicit none

    ! Allocate and set managed memory probin parameters
    allocate(cu_do_constant_volume_burn)
    allocate(cu_call_eos_in_rhs)
    allocate(cu_centered_diff_jac)
    allocate(cu_dt_crit)
    allocate(cu_tmax)
    allocate(cu_renormalize_abundances)
    allocate(cu_burning_mode)
    allocate(cu_burning_mode_factor)
    allocate(cu_integrate_temperature)
    allocate(cu_integrate_energy)
    

    cu_do_constant_volume_burn = do_constant_volume_burn
    cu_call_eos_in_rhs = call_eos_in_rhs
    cu_centered_diff_jac = centered_diff_jac
    cu_dt_crit = dT_crit
    cu_tmax = tmax
    cu_renormalize_abundances = renormalize_abundances
    cu_burning_mode = burning_mode
    cu_burning_mode_factor = burning_mode_factor
    cu_integrate_temperature = integrate_temperature
    cu_integrate_energy = integrate_energy
  end subroutine managed_probin_init

  subroutine managed_probin_finalize()
    
    implicit none

    ! Deallocate managed memory probin parameters
    deallocate(cu_do_constant_volume_burn)
    deallocate(cu_call_eos_in_rhs)
    deallocate(cu_centered_diff_jac)
    deallocate(cu_dt_crit)
    deallocate(cu_tmax)
    deallocate(cu_renormalize_abundances)
    deallocate(cu_burning_mode)
    deallocate(cu_burning_mode_factor)
    deallocate(cu_integrate_temperature)
    deallocate(cu_integrate_energy)
  end subroutine managed_probin_finalize

end module managed_probin_module
