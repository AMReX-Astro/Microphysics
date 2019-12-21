! a simple code to check the analytic Jacobian via numerical 
! differencing

program testjacobian

  use amrex_constants_module
  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real
  use network
  use eos_module
  use burner_module
  use burner_aux_module

  implicit none

  real(rt) :: dens, temp
  real(rt), dimension(nspec) :: Xin
  real(rt), dimension(nspec_advance+1) :: y, ydot
  real(rt), dimension(nspec_advance+1) :: yp, ym
  real(rt), dimension(nspec_advance+1) :: ydotp, ydotm
  real(rt), dimension(nspec_advance+1,nspec_advance+1) :: pd
  real(rt) :: enucdot

  real(rt) :: rpar
  integer :: ipar

  integer :: ic12, io16, img24
  integer :: i, j, n

  real(rt), parameter :: delta = 0.001e0_rt
  real(rt) :: num_jac

  call network_init()
  call eos_init()

  ic12 = network_species_index("carbon-12")
  io16 = network_species_index("oxygen-16")
  img24 = network_species_index("magnesium-24")

  if (ic12 < 0 .or. io16 < 0 .or. img24 < 0) then
     call amrex_error("ERROR: species index not defined")
  endif
  
  dens = 2.6e9_rt
  temp = 7.e8_rt

  Xin(ic12) = 0.5_rt
  Xin(io16) = 0.5_rt
  Xin(img24) = 0.0_rt


  den_eos(1) = dens
  temp_eos(1) = temp
  xn_eos(1,:) = Xin(:)
  
  call eos(eos_input_rt, den_eos, temp_eos, &
           npts, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           .false.)


  print *, 'evaluating the RHS...'

  ! load the state
  y(1) = Xin(ic12)
  y(nspec_advance+1) = temp

  ! set the burner_aux variables
  dens_pass = dens
  c_p_pass = cp_eos(1)
  dhdx_pass(:) = dhdX_eos(1,:)
  X_O16_pass = Xin(io16)
  

  call f_rhs(nspec_advance+1, ZERO, y, ydot, rpar, ipar)

  call jac(nspec_advance+1, ZERO, y, 0, 0, pd, nspec_advance+1, rpar, ipar)

999 format(1x, "df(",i1,")/dy(",i1,")", g20.10, g20.10)

  do j = 1, nspec_advance+1

     yp(:) = y(:)
     ym(:) = y(:)

     yp(j) = (1.e0_rt + delta)*y(j)
     call f_rhs(nspec_advance+1, ZERO, yp, ydotp, rpar, ipar)
     
     ym(j) = (1.e0_rt - delta)*y(j)
     call f_rhs(nspec_advance+1, ZERO, ym, ydotm, rpar, ipar)        

     do i = 1, nspec_advance+1
        
        num_jac = (ydotp(i) - ydotm(i))/(yp(j) - ym(j))

        write (*,999) i, j, num_jac, pd(i,j)

     enddo
  enddo

end program testjacobian
