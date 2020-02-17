! rpar is a set of real quantities that we use to convert auxiliary data back
! and forth between the Microphysics data format and the specific data format
! used by the integrators.

! SDC TODO: we need to store the advective sources here.  Also some bits
! are not needed in the SDC implementation (like the cv and cp stuff)

module vode_rpar_indices

#ifdef SIMPLIFIED_SDC
  use sdc_type_module, only: SVAR, SVAR_EVOLVE
#elif TRUE_SDC
  use actual_network, only: nspec, nspec_evolve
#else
  use actual_network, only: nspec, nspec_evolve
  use burn_type_module, only: neqs
#endif

  use amrex_fort_module, only : rt => amrex_real
  implicit none


#ifdef SIMPLIFIED_SDC

#if defined(SDC_EVOLVE_ENERGY)

  ! Note: we require these components to be first, to allow for offset
  ! indexing with irp_ydot_a and irp_u_init
  integer, parameter :: irp_SRHO = 1
  integer, parameter :: irp_SMX  = 2
  integer, parameter :: irp_SMY  = 3
  integer, parameter :: irp_SMZ  = 4

  ! these are the SVAR advective terms (including those that we do
  ! not explicitly evolve).  Note: these need to be indexed using
  ! the indicies defined in sdc_type_module
  integer, parameter :: irp_ydot_a = 5

  ! these are the SVAR - SVAR_EVOLVE initial values of the unevolved
  ! components of the conserved state.  Note: these need to be indexed
  ! using the components defined above (irp_SRHO, irp_SMX, ...)
  integer, parameter :: irp_u_init = irp_ydot_a + SVAR

  integer, parameter :: irp_T_from_eden = irp_u_init + SVAR - SVAR_EVOLVE
  integer, parameter :: irp_self_heat = irp_T_from_eden + 1

#elif defined(SDC_EVOLVE_ENTHALPY)

  ! Note: we require these components to be first, to allow for offset
  ! indexing with irp_ydot_a
  integer, parameter :: irp_SRHO = 1

  ! These are the SVAR advective terms.  Note: these need to be
  ! indexed using the indices defined in sdc_type_module
  integer, parameter :: irp_ydot_a = 2

  ! This is the pressure component, carried in case we
  ! wish to call the EOS using pressure as an input
  integer, parameter :: irp_p0 = irp_ydot_a + SVAR

  ! These are various bookkeeping parameters
  integer, parameter :: irp_self_heat = irp_p0 + 1

#endif

  integer, parameter :: irp_t0 = irp_self_heat + 1

#ifdef NONAKA_PLOT

  integer, parameter :: irp_i = irp_t0 + 1
  integer, parameter :: irp_j = irp_t0 + 2
  integer, parameter :: irp_k = irp_t0 + 3
  integer, parameter :: irp_iter = irp_k + 1

  integer, parameter :: n_rpar_comps = irp_iter

#else

  integer, parameter :: n_rpar_comps = irp_t0

#endif

#elif TRUE_SDC
  ! f_source is function we are zeroing.  There are nspec_evolve + 2
  ! components (density and energy), since those are the unknowns for
  ! the nonlinear system
  integer, parameter :: irp_f_source = 1

  ! dt is the timestep (1 component)
  integer, parameter :: irp_dt = irp_f_source + nspec_evolve + 2

  ! mom is the momentum (3 components)
  integer, parameter :: irp_mom = irp_dt + 1

  ! evar is the other energy variable (rho e if we are solving for rho E, and vice versa)
  integer, parameter :: irp_evar = irp_mom + 3

  ! the temperature -- used as a guess in the EOS
  integer, parameter :: irp_temp = irp_evar + 1

  ! the unevolved species -- note: unevolved here means not reacting
  integer, parameter :: irp_spec = irp_temp + 1  ! nspec - nspec_evolve components

  integer, parameter :: n_rpar_comps = nspec_evolve + 8 + (nspec - nspec_evolve)

#else
  integer, parameter :: n_not_evolved = nspec - nspec_evolve

  integer, parameter :: irp_dens = 1
  integer, parameter :: irp_cx = irp_dens + 1
  integer, parameter :: irp_nspec = irp_cx + 1
  integer, parameter :: irp_abar = irp_nspec + n_not_evolved
  integer, parameter :: irp_zbar = irp_abar + 1
  integer, parameter :: irp_eta = irp_zbar + 1
  integer, parameter :: irp_ye = irp_eta + 1
  integer, parameter :: irp_cs = irp_ye + 1
  integer, parameter :: irp_dx = irp_cs + 1
  integer, parameter :: irp_t_sound = irp_dx + 1
  integer, parameter :: irp_y_init = irp_t_sound + 1
  integer, parameter :: irp_self_heat = irp_y_init + neqs
  integer, parameter :: irp_Told = irp_self_heat + 1
  integer, parameter :: irp_dcxdt = irp_Told + 1
  integer, parameter :: irp_t0 = irp_dcxdt + 1

  integer, parameter :: n_rpar_comps = irp_t0 + 1
#endif


end module vode_rpar_indices
