module managed_probin_module

  use bl_types
  
  implicit none

!   logical, &
! #ifdef CUDA       
!        managed, &
! #endif       
!        allocatable, save :: cu_use_tables
  
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
  
  integer, &
#ifdef CUDA       
       managed, &
#endif              
       allocatable, save :: cu_jacobian
  
  logical, &
#ifdef CUDA       
       managed, &
#endif              
       allocatable, save :: cu_burner_verbose, cu_retry_burn
  
  real(dp_t), &
#ifdef CUDA       
       managed, &
#endif              
       allocatable, save :: cu_rtol_spec, cu_rtol_temp, cu_rtol_enuc, &
       cu_atol_spec, cu_atol_temp, cu_atol_enuc, cu_retry_burn_factor, cu_retry_burn_max_change

  !$acc declare &
  !$acc create(cu_do_constant_volume_burn, cu_call_eos_in_rhs, cu_dt_crit) &
  !$acc create(cu_centered_diff_jac, cu_tmax, cu_renormalize_abundances) &
  !$acc create(cu_burning_mode, cu_burning_mode_factor) &
  !$acc create(cu_integrate_temperature, cu_integrate_energy) &
  !$acc create(cu_jacobian, cu_burner_verbose, cu_rtol_spec, cu_rtol_temp, cu_rtol_enuc) &
  !$acc create(cu_atol_spec, cu_atol_temp, cu_atol_enuc) &
  !$acc create(cu_retry_burn, cu_retry_burn_factor, cu_retry_burn_max_change)
  
contains

  subroutine managed_probin_init()
    use extern_probin_module, only: do_constant_volume_burn, dT_crit, call_eos_in_rhs, &
         centered_diff_jac, renormalize_abundances, &
         burning_mode, burning_mode_factor, &
         integrate_temperature, integrate_energy, &
         jacobian, burner_verbose, &
         rtol_spec, rtol_temp, rtol_enuc, &
         atol_spec, atol_temp, atol_enuc, &
         retry_burn, &
         retry_burn_factor, retry_burn_max_change!, &
!         use_tables
    use probin_module, only: tmax

    implicit none

    ! Allocate and set managed memory probin parameters
!    allocate(cu_use_tables)
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
    allocate(cu_jacobian)
    allocate(cu_burner_verbose)
    allocate(cu_rtol_spec)
    allocate(cu_rtol_temp)
    allocate(cu_rtol_enuc)
    allocate(cu_atol_spec)
    allocate(cu_atol_temp)
    allocate(cu_atol_enuc)
    allocate(cu_retry_burn)
    allocate(cu_retry_burn_factor)
    allocate(cu_retry_burn_max_change)

!    cu_use_tables = use_tables
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
    cu_jacobian = jacobian
    cu_burner_verbose = burner_verbose
    cu_rtol_spec = rtol_spec
    cu_rtol_temp = rtol_temp
    cu_rtol_enuc = rtol_enuc
    cu_atol_spec = atol_spec
    cu_atol_temp = atol_temp
    cu_atol_enuc = atol_enuc
    cu_retry_burn = retry_burn
    cu_retry_burn_factor = retry_burn_factor
    cu_retry_burn_max_change = retry_burn_max_change
  end subroutine managed_probin_init

  subroutine managed_probin_finalize()
    
    implicit none

    ! Deallocate managed memory probin parameters
!    deallocate(cu_use_tables)
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
    deallocate(cu_jacobian)
    deallocate(cu_burner_verbose)
    deallocate(cu_rtol_spec)
    deallocate(cu_rtol_temp)
    deallocate(cu_rtol_enuc)
    deallocate(cu_atol_spec)
    deallocate(cu_atol_temp)
    deallocate(cu_atol_enuc)
    deallocate(cu_retry_burn)
    deallocate(cu_retry_burn_factor)
    deallocate(cu_retry_burn_max_change)
  end subroutine managed_probin_finalize

end module managed_probin_module
