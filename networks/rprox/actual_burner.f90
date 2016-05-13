! This module contains a version of Wallace and Woosley's (ApJS 45,389
! (1981)) rprox reaction network burner.

module actual_burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use eos_type_module
  use network
  use actual_burner_data
  use burn_type_module

contains

  subroutine actual_burner_init()

    use integrator_module, only: integrator_init
    use integration_data, only: temp_scale

    implicit none

    temp_scale = 1.0d9

    reac_names(irlambCNO) = "rlambdaCNO"
    reac_names(irag15o)   = "rag15o"
    reac_names(irr1)      = "rr1"
    reac_names(irag16o)   = "rag16o"
    reac_names(irpg16o)   = "rpg16o"
    reac_names(irpg17f)   = "rpg17f"
    reac_names(irgp17f)   = "rgp17f"
    reac_names(irlambda2) = "rlambda2"
    reac_names(irap14o)   = "rap14o"
    reac_names(irs1)      = "rs1"
    reac_names(irlambda1) = "rlambda1"
    reac_names(ir3a)      = "r3a"
    reac_names(irpg12c)   = "rpg12c"
    reac_names(irwk14o)   = "wk14o"
    reac_names(irwk17f)   = "wk17f"
    reac_names(irwk15o)   = "wk15o"
    reac_names(irLweak)   = "Lweak"
    reac_names(irla2)     = "la2"

    call integrator_init()

  end subroutine actual_burner_init



  subroutine actual_burner(state_in, state_out, dt, time)  

    use integrator_module, only: integrator

    implicit none

    type (burn_t),    intent(in   ) :: state_in
    type (burn_t),    intent(inout) :: state_out
    double precision, intent(in   ) :: dt, time    

    double precision :: T9

    T9 = state_in % T / 10.0**9

    ! Only burn if 0.2 < T9 < 2.5 or X(H1) > 0.05.
    ! The last restriction is a kludge based on the last paragraph of WW81.

    if ((T9 .gt. 0.2d0 .and. T9 .lt. 2.5d0) .or. state_in % xn(ih1) > 0.05d0) then

       ! Call the integration routine.

       call integrator(state_in, state_out, dt, time)

    endif
    
  end subroutine actual_burner

end module actual_burner_module
