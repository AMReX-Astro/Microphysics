program cj_det

  use bl_types, only: dp_t

  use actual_rhs_module
  use eos_type_module
  use eos_module
  use network
  use probin_module, only: smallx
  use cj_det_module

  implicit none

  type(eos_t) :: eos_state_fuel, eos_state_ash

  real(dp_t) :: q_burn
  real(dp_t), parameter :: rho_min_fac = 0.9_dp_t, rho_max_fac = 10.0_dp_t
  integer, parameter :: npts_ad = 150

  ! set the unburned (fuel) state
  eos_state_fuel % rho = 1.e7
  eos_state_fuel % T = 1.e8
  eos_state_fuel % xn(:) = smallx
  eos_state_fuel % xn(1) = 1.0 - (nspec - 1)*smallx

  call eos(eos_input_rt, eos_state_fuel)


  ! set the ash composition
  eos_state_ash = eos_state_fuel
  eos_state_ash % xn(:) = smallx
  eos_state_ash % xn(nspec) = 1.0 - (nspec - 1)*smallx

  ! get the q value
  call ener_gener_rate(eos_state_ash % xn(:) - eos_state_fuel % xn(:), q_burn)

  ! store the shock adiabat and the detonation adiabat
  rho_min = rho_min_fac * eos_state_fuel % rho
  rho_max = rho_max_fac * eos_state_fuel % rho
  dlogrho = (log10(rho_max) - log10(rho_min))/(npts-1)

  do n = 0, npts_ad-1
     eos_state_ash % rho = 10.0_dp_t**(dlog10(rho_min) + n*dlogrho)

     call adiabat(eos_state_fuel, eos_state_ash, 0.0_rt)
     call adiabat(eos_state_fuel, eos_state_ash, q_burn)

     print *, rho2, p2_shock, p2_det

  enddo

  ! solve for the detonation conditions

end program cj_det

