program cj_det

  use amrex_fort_module, only : rt => amrex_real

  use actual_rhs_module
  use eos_type_module
  use eos_module
  use microphysics_module
  use network
  use runtime_init_module
  use probin_module, only: smallx
  use cj_det_module

  implicit none

  type(eos_t) :: eos_state_fuel, eos_state_ash

  real(rt) :: q_burn
  real(rt), parameter :: rho_min_fac = 0.9_rt, rho_max_fac = 10.0_rt
  integer, parameter :: npts_ad = 150

  real(rt) :: rho_min, rho_max, dlogrho, p2_shock, p2_det, D_cj, rho_cj, p_cj, cs_det
  integer :: n
  integer, parameter :: npts = 100

  integer :: lun, istatus

  ! runtime
  call runtime_init(.true.)

  call microphysics_init()

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

  ! get the q value -- we need the change in molar fractions
  call ener_gener_rate(eos_state_ash % xn(:)/aion(:) - eos_state_fuel % xn(:)/aion(:), q_burn)

  ! store the shock adiabat and the detonation adiabat
  rho_min = rho_min_fac * eos_state_fuel % rho
  rho_max = rho_max_fac * eos_state_fuel % rho
  dlogrho = (log10(rho_max) - log10(rho_min))/(npts-1)

  ! initial guess
  eos_state_ash % T = eos_state_fuel % T


  ! now let's get the CJ velocity
  call cj_cond(eos_state_fuel, eos_state_ash, q_burn)

  ! we get this from the mass flux: rho_1 v_1 = j
  D_cj = (1.0_rt / eos_state_fuel % rho) * sqrt( &
       (eos_state_ash % p - eos_state_fuel % p) / &
       (1.0_rt/eos_state_fuel % rho - 1.0_rt/eos_state_ash % rho))

  rho_cj = eos_state_ash % rho
  p_cj = eos_state_ash % p

  ! output info along with points on the Hugoniot

  open(newunit=lun, file="hugoniot.txt", status="unknown")

  write(lun, *) "# initial rho = ", eos_state_fuel % rho
  write(lun, *) "# initial p = ", eos_state_fuel % p
  write(lun, *) "# ash rho = ", eos_state_ash % rho
  write(lun, *) "# ash p = ", eos_state_ash % p
  write(lun, *) "# CJ speed = ", D_cj

  ! test
  cs_det = sqrt(eos_state_ash % gam1 * eos_state_ash % p / eos_state_ash % rho)

  ! this checks that the slope of the Rayleigh line is rho_2 * cs_2
  !print *, eos_state_ash % rho * cs_det, eos_state_fuel % rho * D_cj

  do n = 0, npts_ad-1

     eos_state_ash % rho = 10.0_rt**(dlog10(rho_min) + n*dlogrho)

     call adiabat(eos_state_fuel, eos_state_ash, 0.0_rt, istatus)
     p2_shock = eos_state_ash % p

     if (istatus == -1) then
        exit
     endif

     call adiabat(eos_state_fuel, eos_state_ash, q_burn, istatus)
     p2_det = eos_state_ash % p

     if (istatus == -1) then
        exit
     endif

     write(lun, *) eos_state_ash % rho, p2_shock, p2_det

  enddo

  close(unit=lun)



end program cj_det

