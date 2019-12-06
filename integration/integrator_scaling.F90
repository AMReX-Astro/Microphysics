module integrator_scaling_module

  use extern_probin_module, only: reactions_density_scale, reactions_energy_scale, reactions_temperature_scale
  use microphysics_type_module, only: rt, ONE

  implicit none

  real(rt), allocatable :: dens_scale, ener_scale, temp_scale
  real(rt), allocatable :: inv_dens_scale, inv_ener_scale, inv_temp_scale

#ifdef AMREX_USE_CUDA
  attributes(managed) :: dens_scale, ener_scale, temp_scale
  attributes(managed) :: inv_dens_scale, inv_ener_scale, inv_temp_scale
#endif

contains

  subroutine integrator_scaling_init

    implicit none

    allocate(dens_scale, ener_scale, temp_scale)
    allocate(inv_dens_scale, inv_ener_scale, inv_temp_scale)

    dens_scale = reactions_density_scale
    ener_scale = reactions_energy_scale
    temp_scale = reactions_temperature_scale

    inv_dens_scale = ONE / dens_scale
    inv_ener_scale = ONE / ener_scale
    inv_temp_scale = ONE / temp_scale

  end subroutine integrator_scaling_init
  
end module integrator_scaling_module
