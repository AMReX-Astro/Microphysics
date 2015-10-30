! a simple code to check the analytic Jacobian via numerical 
! differencing

program testjacobian

  use bl_types
  use bl_constants_module
  use bl_error_module
  use network
  use eos_module
  use eos_type_module
  use burner_module
  use network_indices
  use rpar_indices

  implicit none

  real(kind=dp_t) :: dens, temp
  integer, parameter :: neq = nevolve+1
  real(kind=dp_t), dimension(nspec) :: Xin
  real(kind=dp_t), dimension(neq) :: y, ydot
  real(kind=dp_t), dimension(neq) :: yp, ym
  real(kind=dp_t), dimension(neq) :: ydotp, ydotm
  real(kind=dp_t), dimension(neq,neq) :: pd
  real(kind=dp_t) :: enucdot

  real(kind=dp_t), allocatable :: rpar(:)
  integer :: ipar

  type(eos_t) :: eos_state

  integer :: i, j, n

  real(kind=dp_t), parameter :: delta = 0.001d0
  real(kind=dp_t), parameter :: SMALL = 1.d-12
  real(kind=dp_t) :: num_jac

  character(len=16) :: namei,namej

  call network_init()
  call eos_init()

  allocate(rpar(n_rpar_comps))

  dens = 1.0e6_dp_t
  temp = 2.e8_dp_t

  Xin = 0.0d0
  Xin(ihe4_) = 0.5d0
  Xin(ic12_) = 0.5d0

  rpar(irp_dens) = dens
  rpar(irp_Teos) = temp
  rpar(irp_Y56) = ZERO

  eos_state%rho = dens
  eos_state%T   = temp
  eos_state%xn  = Xin

  call eos(eos_input_rt,eos_state)

  rpar(irp_cp) = eos_state%cp
  rpar(irp_dhdX:irp_dhdX+nspec-1) = eos_state%dhdX

 print *, 'evaluating the RHS...'

  ! load the state
  y(1:nevolve) = Xin(1:nevolve)
  y(neq) = temp

  call f_rhs(neq, ZERO, y, ydot, rpar, ipar)

  call jac(neq, ZERO, y, 0, 0, pd, neq, rpar, ipar)

888 format(a,"-derivatives that don't match:")
999 format(5x, "df(",a,")/dy(",a,")", g18.10, g18.10, g18.10)

  do j = 1, neq

     yp(:) = y(:)
     ym(:) = y(:)

     yp(j) = (1.d0 + delta)*y(j) 
     call f_rhs(neq, ZERO, yp, ydotp, rpar, ipar)
     
     ym(j) = (1.d0 - delta)*y(j) 
     call f_rhs(neq, ZERO, ym, ydotm, rpar, ipar)        
     
     if (j==neq) then
        namej = "T"     
     else
        namej = short_spec_names(j)
     endif

     write(*,888) trim(namej)

     do i = 1, neq
        
        num_jac = (ydotp(i) - ydotm(i))/(yp(j) - ym(j) + SMALL)

        if (i==neq) then
           namei = "T"
        else
           namei = short_spec_names(i)
        endif

        ! only dump the ones that don't match
        if (num_jac /= ZERO) then
           write (*,999) trim(namei), &
                trim(namej), num_jac, pd(i,j), (num_jac-pd(i,j))/num_jac
        else
           write (*,999) trim(namei), &
                trim(namej), num_jac, pd(i,j)
        endif
        
     enddo
  enddo


end program testjacobian
