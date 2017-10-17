module sdc_type_module

  use actual_network, only: nspec
  use bl_types, only: dp_t

  implicit none

  ! A generic structure holding data necessary to do a nuclear burn
  ! in the SDC formalism.

#if (SDC_METHOD == 1)
  ! these indicies represent the order that the conserved state comes
  ! into the ODE integration from the hydro code.
  !
  ! they also represent the order of the advective sources
  !
  ! integrate rho*X, internal energy, total energy
  ! carry momentum as an unevolved variable

  integer, parameter :: SEDEN = 1
  integer, parameter :: SEINT = 2
  integer, parameter :: SFS   = 3
  integer, parameter :: SRHO  = SFS + nspec
  integer, parameter :: SMX   = SRHO + 1
  integer, parameter :: SMY   = SRHO + 2
  integer, parameter :: SMZ   = SRHO + 3

  integer, parameter :: SVAR  = SMZ
  integer, parameter :: SVAR_EVOLVE = SRHO - 1
#elif (SDC_METHOD == 2)
  ! integrate rho*X (species masses) and rho*h (enthalpy)
  ! carry pressure for EOS calls in RHS

  integer, parameter :: SSPEC = 1
  integer, parameter :: SENTH = SSPEC + nspec
  integer, parameter :: SVAR  = SENTH
#endif

  type :: sdc_t

     real(dp_t) :: y(SVAR)
     real(dp_t) :: ydot_a(SVAR)

#if (SDC_METHOD == 1)
     logical :: T_from_eden
#elif (SDC_METHOD == 2)
     ! Pressure in case we wish to use it for EOS calls
     real(dp_t) :: p0
#endif

     integer :: i
     integer :: j
     integer :: k

     integer :: n_rhs
     integer :: n_jac

     integer :: sdc_iter

  end type sdc_t

end module sdc_type_module
