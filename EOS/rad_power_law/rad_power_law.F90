! This is an artificial equation of state used primarily for radiation tests.
!
! It is defined by the relationship:
! c_v = K * rho**m * T**(-n)
! where K, m, and n are user-defined constant parameters.
!
! Ignoring the integration constant, we thus have the relationships:
! e = K * rho**m * T**(1-n) / (1 - n)
! T = ((1 - n) * e * rho**(-m) / K)**(1/(1-n))
!
! Consequently the only input modes supported are eos_input_rt and eos_input_re.
! Pressure and Gamma_1 are not defined, so this EOS cannot be used for hydro.

module actual_eos_module

  use amrex_fort_module, only: rt => amrex_real
  use eos_type_module, only: eos_t, eos_input_rt, eos_input_re

  implicit none

  character (len=64), parameter :: eos_name = "rad_power_law"

  real(rt), allocatable, save :: const_c_v, c_v_exp_m, c_v_exp_n

#ifdef AMREX_USE_CUDA
  attributes(managed) :: const_c_v, c_v_exp_m, c_v_exp_n
#endif

contains

  subroutine actual_eos_init

    use extern_probin_module, only: eos_const_c_v, eos_c_v_exp_m, eos_c_v_exp_n
    use amrex_error_module, only: amrex_error

    implicit none

    allocate(const_c_v)
    allocate(c_v_exp_m)
    allocate(c_v_exp_n)

    const_c_v = eos_const_c_v
    c_v_exp_m = eos_c_v_exp_m
    c_v_exp_n = eos_c_v_exp_n

    if (const_c_v .le. 0.e0_rt) then
       call amrex_error("eos_const_c_v must be > 0")
    end if

    if (c_v_exp_n .eq. 1.0e0_rt) then
       call amrex_error("eos_c_v_exp_n == 1 is unsupported")
    end if

  end subroutine actual_eos_init



  subroutine actual_eos_finalize()

    implicit none

    deallocate(const_c_v)
    deallocate(c_v_exp_m)
    deallocate(c_v_exp_n)

  end subroutine actual_eos_finalize



  subroutine is_input_valid(input, valid)
    implicit none
    integer, intent(in) :: input
    logical, intent(out) :: valid

    !$gpu

    valid = .true.

    if (input == eos_input_rh .or. &
        input == eos_input_tp .or. &
        input == eos_input_rp .or. &
        input == eos_input_ps .or. &
        input == eos_input_ph .or. &
        input == eos_input_th) then
       valid = .false.
    end if
  end subroutine is_input_valid


  subroutine actual_eos(input, state)

#ifndef AMREX_USE_GPU
    use amrex_error_module, only: amrex_error
#endif

    implicit none

    integer,      intent(in   ) :: input
    type (eos_t), intent(inout) :: state

    !$gpu

    select case (input)

    case (eos_input_rt)

       state % cv = const_c_v * state % rho**c_v_exp_m * state % T**(-c_v_exp_n)
       state % e = const_c_v * state % rho**c_v_exp_m * state % T**(1 - c_v_exp_n) / (1 - c_v_exp_n)

    case (eos_input_re)

       state % T = ((1 - c_v_exp_n) * state % e * state % rho**(-c_v_exp_m) / const_c_v)**(1.0_rt / (1.0_rt - c_v_exp_n))
       state % cv = const_c_v * state % rho**c_v_exp_m * state % T**(-c_v_exp_n)

    case default

#ifndef AMREX_USE_GPU
       call amrex_error('EOS: invalid input.')
#endif

    end select

    ! Set some data to nonsense values so that things intentionally go wrong
    ! if this EOS is somehow used for hydro.

    state % p    = -1.e0_rt
    state % gam1 = -1.e0_rt
    state % cs   = -1.e0_rt
    state % s    = -1.e0_rt
    state % h    = -1.e0_rt

    ! Give dpdr a value for the purposes of the composition_derivatives routine.

    state % dPdr = 0.e0_rt

  end subroutine actual_eos

end module actual_eos_module
