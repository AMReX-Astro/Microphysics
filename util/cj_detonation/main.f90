program cj_det

  use bl_types, only: dp_t

  use actual_rhs_module
  use eos_type_module
  use eos_module
  use network
  use probin_module, only: smallx

  implicit none

  type(eos_t) :: eos_state_fuel, eos_state_ash

  real(dp_t) :: q_burn

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


  ! solve for the detonation conditions

end program cj_det
