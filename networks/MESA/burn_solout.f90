!
! burn_solout routine to calculate energy generation rate
!
!   this routine is called after every step of the burn
!   and as a result the actual call to this subroutine 
!   is buried deep inside the MESA source code
!
!   this routine has little to no modification compared to
!   the MESA source version
!

subroutine burn_solout(step, told, time, n, x, rwork_y, iwork_y, &
                       interp_y, lrpar, rpar, lipar, ipar, irtrn)

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
   use amrex_fort_module, only : rt => amrex_real

   implicit none

   integer, intent(in) :: step, n, lrpar, lipar
   real(rt), intent(in) :: told, time
   real(rt), intent(inout) :: x(n)
   ! x can be modified if necessary to keep it in valid range of 
   ! possible solutions.
   !real(rt), intent(inout), target :: rpar(lrpar), rwork_y(*)
   !integer, intent(inout), target :: ipar(lipar), iwork_y(*)
   real(rt), intent(inout), target :: rwork_y(*)
   integer, intent(inout), target :: iwork_y(*)
   real(rt), intent(inout), pointer :: rpar(:)
   integer, intent(inout), pointer :: ipar(:)
   interface
      include 'num_interp_y.dek'
   end interface
   integer, intent(out) :: irtrn ! < 0 causes solver to return to 
                                 !calling program.

   real(rt) :: logRho, logT
   real(rt) :: d_eps_nuc_dRho, &
          d_eps_nuc_dT, eps_nuc_categories(num_rvs, num_categories)
   real(rt), dimension(:, :), pointer :: reaction_eps_nuc
   real(rt), dimension(:, :), pointer :: rate_screened, rate_raw

   real(rt) :: Rho, T, xsum, d_eps_nuc_dx(species), dx, &
          dt, xh, xhe, mass_correction

   integer :: info, i, j, lwork, cid
   real(rt), dimension(species) :: dabar_dx, dzbar_dx, dmc_dx
   real(rt), pointer :: work(:), rate_factors(:), category_factors(:)

   real(rt) :: abar, zbar, eps_nuc, ye, z2bar

   logical, parameter :: reuse_given_rates = .false.

   ! this is new to version 5118: (make sure all interfaces referencing 
   !  this routine have been updated too)
   !real(rt) :: d_eta_dlnT, d_eta_dlnRho

   irtrn = 0
   if (time == 0) return

   x(1:n) = max(0d0, min(1d0, x(1:n)))
   xsum = sum(x(1:n))
   x(1:n) = x(1:n)/xsum

   info = 0
   lwork = net_work_size(handle_net, info)
   if (info /= 0) call amrex_error("error in burn_solout (MESA) routine")

   allocate(work(lwork),  &
      rate_factors(num_reactions), category_factors(num_categories), &
      rate_screened(num_rvs, num_reactions),  &
      reaction_eps_nuc(num_rvs, num_reactions),  &
      rate_raw(num_rvs, num_reactions), stat=info)

   if (info /= 0) call amrex_error("error allocating in burn_solout (MESA)")

   T = burn_temp
   Rho = burn_rho
   logT = log10(T)
   logRho = log10(Rho)
   xin = x

   call composition_info( &
         species, chem_id, xin, xh, xhe, abar, zbar, z2bar, ye, &
         mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)

   rate_factors(:) = 1
   category_factors(:) = 1
   ! this is version 4942:
   call net_get(handle_net, species, num_reactions,  &
         xin, T, logT, Rho, logRho,  &
         abar, zbar, z2bar, ye, eta, rate_factors, category_factors,  &
         std_reaction_Qs, std_reaction_neuQs, &
         eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
         dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
         screening_mode, theta_e_for_graboske_et_al,     &
         rate_screened, rate_raw, reuse_given_rates,  &
         reaction_eps_nuc, eps_nuc_categories, eps_neu_total, &
         lwork, work, info)

   ! this is version 5118:
   !call net_get(handle_net, species, num_reactions,  &
         !xin, T, logT, Rho, logRho,  &
         !abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
         !rate_factors, category_factors,  &
         !std_reaction_Qs, std_reaction_neuQs, &
         !eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
         !dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
         !screening_mode, theta_e_for_graboske_et_al,     &
         !rate_screened, rate_raw, reuse_given_rates,  &
         !reaction_eps_nuc, eps_nuc_categories, eps_neu_total, &
         !lwork, work, info)

   if (info /= 0) then
      call amrex_error("error bad return from net_get (MESA)")
   end if

   deallocate(work, rate_factors, category_factors, rate_screened, &
      reaction_eps_nuc, rate_raw)

   dt = time - told

   ! set burn_ergs according to change from initial abundances
   burn_neu_total = burn_neu_total + eps_neu_total*dt
   burn_ergs = 0
   eps_nuc = 0
   xsum = 0
   do i=1,species
      cid = chem_id(i)
      dx = x(i) - x_initial(i)
      xsum = xsum + x(i)

      ! units of burn_ergs are MeV/g
      burn_ergs = burn_ergs +  &
         (chem_isos% binding_energy(cid) - chem_isos% Z(cid)*del_Mp -  &
            chem_isos% N(cid)*del_Mn)*dx/chem_isos% W(cid)

      dx = x(i) - x_previous(i)
      eps_nuc = eps_nuc +  &
         (chem_isos% binding_energy(cid) - chem_isos% Z(cid)*del_Mp -  &
            chem_isos% N(cid)*del_Mn)*dx/chem_isos% W(cid)
   end do

   ! Qconv converts from MeV/g to ergs/g
   eps_nuc = eps_nuc*Qconv/dt - eps_neu_total
   burn_ergs = burn_ergs*Qconv - burn_neu_total

   burn_ergs_total = burn_ergs_total + burn_ergs
   eps_nuc_total = eps_nuc_total + eps_nuc

   !--- DEBUG
   if (.false.) then
       write(*,'(a40,i12)') 'step',step
       write(*,'(a40,99(1pd26.16))') 'eps_nuc',eps_nuc
       !write(*,'(a40,99(1pd26.16))') 'burn_ergs',burn_ergs
       !write(*,'(a40,99(1pd26.16))') 'burn_ergs_total',burn_ergs_total
       write(*,'(a40,99(1pd26.16))') 'eps_nuc_total',eps_nuc_total
       write(*,*)
   endif
   !--- DEBUG

   x_previous(1:species) = x(1:species)

   if (time < data_output_min_t) return

   do j=1, species
      if (x(j) > peak_abundance(j)) then
         peak_abundance(j) = x(j)
         peak_time(j) = time
      end if
   end do

end subroutine burn_solout
