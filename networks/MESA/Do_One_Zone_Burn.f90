!###################################################################
! One Zone Burn Routine:
!
! based on MESA version 4088  mesa/net/test/one_zone_burn.f
!
! Ryan Orvedahl Aug 2, 2012
!
! revised:
!       Sept 24, 2012
!               -remove initialize_module_ptrs, shutdown_module_ptrs
!       Mar 31, 2013
!               -remove extraneous comments/lines of code
!               -fix the energy calculation
!       ~April 2013
!               -update to version 4849
!       Jul 2, 2013
!               -update to version 4942 (this may have taken place earlier)
!               -remove set_rattab_rates and net_setup_tables calls
!               -get energy to agree with mesa/net/test/ test_burn.f runs
!                  by returning burn_ergs, not burn_ergs_total
!
!###################################################################
subroutine Do_One_Zone_Burn(density, temperature, tstop, x_mesa_in, &
                    burn_erg_out, ending_x)
   
   ! Local
   use net_utils,   only: decsol_choice, solver_choice, handle_net, &
                      num_reactions, species, &
                      x_previous, x_initial, xin, &
                      dxdt, d_dxdt_dT, d_dxdt_dRho, d_dxdt_dx, &
                      peak_time, peak_abundance, eps_nuc_total, &
                      burn_rho, burn_temp, burn_neu_total, eps_neu_total, &
                      eta, data_output_min_t, theta_e_for_graboske_et_al, &
                      screening_mode, category_factors, rate_factors, &
                      burn_ergs_total, burn_ergs
   
   ! MESA
   use net_lib,     only: num_categories, net_1_zone_burn
   use rates_def,   only: std_reaction_neuQs, std_reaction_Qs
   use screen_def,  only: extended_screening

   ! MAESTRO
   use network,     only: nspec
   use eos_module,  only: eos
   use eos_type_module, only: eos_t

   ! AMReX
   use amrex_error_module, only: amrex_error

   implicit none

   ! INPUT:
   real(rt), intent(in) :: density, temperature, tstop
   real(rt), intent(in) :: x_mesa_in(nspec)

   ! OUTPUT:
   real(rt), intent(out) :: burn_erg_out
   real(rt), intent(out) :: ending_x(nspec)

   ! LOCAL:
   type (eos_t) :: eos_state
   integer, save :: eos_input = 1
   logical, save :: do_eos_diag = .false.

   real(rt) :: burn_tend, burn_rtol, burn_atol

   real(rt), dimension(:), pointer :: d_eps_nuc_dx

   integer :: handle, max_steps
   
   real(rt) :: logRho, logT, Rho, T
   
   ! args to control the solver -- see num/public/num_isolve.dek
   real(rt) :: h 
   real(rt) :: max_step_size ! maximal step size.

   ! absolute and relative error tolerances
   real(rt) :: rtol(1) ! relative error tolerance(s)
   real(rt) :: atol(1) ! absolute error tolerance(s)
   integer :: itol ! switch for rtol and atol
   
   integer :: nfcn    ! number of function evaluations
   integer :: njac    ! number of jacobian evaluations
   integer :: nstep   ! number of computed steps
   integer :: naccpt  ! number of accepted steps
   integer :: nrejct  ! number of rejected steps
   
   integer :: ierr, iout, info, caller_id

   logical, save :: clip = .true.

   integer, parameter :: num_times = 1
   
   real(rt) :: dt, time_doing_net
   real(rt), dimension(:), pointer :: times, dxdt_source_term
   real(rt), dimension(:,:), pointer :: log10Ts_f, log10Rhos_f, etas_f
   real(rt), dimension(:), pointer :: log10Ts_f1, log10Rhos_f1, etas_f1

   ! interface for the burn_solout routine so it can be passed
   !  to the MESA burner net_1_zone_burn
   interface
      subroutine burn_solout(nr, xold, x, n, y, rwork_y, iwork_y, &
              interp_y, lrpar,rpar, lipar, ipar, irtrn)
         ! nr is the step number.
         ! x is the current x value; xold is the previous x value.
         ! y is the current y value.
         ! irtrn negative means terminate integration.
         ! rwork_y and iwork_y hold info for interp_y
         ! note that these are not the same as the rwork and iwork arrays for
         !   the solver.
         use amrex_fort_module, only : rt => amrex_real
         use num_lib,   only: safe_log10
         use net_lib,   only: chem_isos, del_Mn, del_Mp, num_categories, &
                          net_work_size, net_get
         use rates_def, only: num_rvs, std_reaction_neuQs, std_reaction_Qs
         use const_def, only: Qconv
         use chem_lib,  only: composition_info
         use net_utils, only: burn_ergs, burn_ergs_total, burn_neu_total, &
                          burn_rho, burn_temp, eps_nuc_total, eps_neu_total, &
                          eta, data_output_min_t, theta_e_for_graboske_et_al, &
                          screening_mode, x_previous, x_initial, xin, &
                          peak_abundance, peak_time, dxdt, d_dxdt_dRho, &
                          d_dxdt_dT, d_dxdt_dx, species, num_reactions, &
                          handle_net
         use network,   only: chem_id
         use amrex_error_module, only: amrex_error
         integer, intent(in) :: nr, n, lrpar, lipar
         real(rt), intent(in) :: xold, x
         real(rt), intent(inout) :: y(n)
         ! y can be modified if necessary to keep it in valid range of 
         !  possible solutions.
         real(rt), intent(inout), target :: rwork_y(*)
         integer, intent(inout), target :: iwork_y(*)
         real(rt), intent(inout), pointer :: rpar(:)
         integer, intent(inout), pointer :: ipar(:)
         interface
            include 'num_interp_y.dek'
         end interface
         integer, intent(out) :: irtrn ! < 0 causes solver to return to 
                                       ! calling program.
      end subroutine burn_solout
   end interface
      
   ! set defaults
   burn_rtol = 5d-10
   burn_atol = 5d-10
   max_steps = 15000

   burn_neu_total = 0
   eps_neu_total = 0      
   eps_nuc_total = 0

   data_output_min_t = 1d-6

   screening_mode = extended_screening
   theta_e_for_graboske_et_al = 1

   !--- initialize burn parameters with incoming values
   burn_tend = tstop
   burn_rho = density
   burn_temp = temperature
   handle = handle_net
   ending_x = x_mesa_in ! ending_x is a copy in order to preserve x_mesa_in

   ierr = 0

   Rho = burn_rho
   T = burn_temp     
   logT = log10(T)
   logRho = log10(burn_rho)

   ! setup eos_state for eos call:
   eos_state%rho = Rho
   eos_state%T = T
   eos_state%xn = ending_x

   !------------------------------------------------------------ 
   ! eos call to get eta
   !------------------------------------------------------------ 
   call eos(eos_input, eos_state)
   eta = eos_state%eta
   
   allocate(rate_factors(num_reactions), category_factors(num_categories), &
         times(num_times), log10Ts_f1(4*num_times), log10Rhos_f1(4*num_times), &
         etas_f1(4*num_times), stat=info)
   if (info /= 0) then
      call amrex_error("error allocating in initialize_module_ptrs")
   end if

   rate_factors(:) = 1
   category_factors(:) = 1

   ! size of times and second index of log10*s_f arrays is num_times = 1 
   times(1) = dt

   log10Ts_f(1:4,1:num_times) => log10Ts_f1(1:4*num_times)
   log10Rhos_f(1:4,1:num_times) => log10Rhos_f1(1:4*num_times)
   etas_f(1:4,1:num_times) => etas_f1(1:4*num_times)

   log10Ts_f(1,1) = logT
   log10Ts_f(2:4,1) = 0

   log10Rhos_f(1,1) = logRho
   log10Rhos_f(2:4,1) = 0

   etas_f(1,1) = eta
   etas_f(2:4,1) = 0

   dxdt_source_term => null()
   
   allocate( xin(species), &
      d_eps_nuc_dx(species), dxdt(species), d_dxdt_dRho(species), &
      d_dxdt_dT(species), d_dxdt_dx(species, species), &
      stat=ierr)

   allocate( &
       x_initial(species), x_previous(species), &
       peak_abundance(species), peak_time(species), stat=ierr)
   if (ierr /= 0) then
      call amrex_error("allocate failed for Do_One_Zone_Burn (MESA)")
   end if

   peak_abundance(:) = 0
   peak_time(:) = 0

   max_step_size = 0
   iout = 1
   itol = 0
   
   h = 1d-14
   
   rtol(:) = burn_rtol
   atol(:) = burn_atol
   
   xin = 0
   xin(:) = x_mesa_in(:)

   xin(:) = xin(:)/sum(xin(:))
      
   x_initial(1:species) = xin(1:species)
   x_previous(1:species) = xin(1:species)
   caller_id = 0

   time_doing_net = -1

   ! set initial value
   burn_ergs_total = 0.0d0
   !--------------------------------------------------------------
   !
   ! MESA one-zone-burn routine:
   !
   !--------------------------------------------------------------
   call net_1_zone_burn( &
         handle, solver_choice, species, num_reactions, 0.0_rt, burn_tend, &
         xin, clip, num_times, times, log10Ts_f1, log10Rhos_f1, etas_f1, &
         dxdt_source_term, rate_factors, category_factors, std_reaction_Qs, &
         std_reaction_neuQs, screening_mode, theta_e_for_graboske_et_al, & 
         h, max_step_size, max_steps, rtol, atol, itol, & 
         decsol_choice, caller_id, burn_solout, iout, ending_x, & 
         nfcn, njac, nstep, naccpt, nrejct, time_doing_net, ierr)

   if (ierr /= 0) then
      call amrex_error("net_1_zone_burn ierr")
   end if

   !--------------------------------------------------------------
   ! Calculate burn_ergs_total:
   !--------------------------------------------------------------
   burn_erg_out = burn_ergs

   ! deallocate the pointers
   deallocate(rate_factors, category_factors, times, log10Ts_f1, &
              log10Rhos_f1, etas_f1)

   deallocate( &
       x_initial, x_previous, &
       peak_abundance, peak_time)

   deallocate( &
      xin, d_eps_nuc_dx, dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx)

   return

end subroutine Do_One_Zone_Burn

