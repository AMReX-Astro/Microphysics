! a simple code to check the analytic Jacobian via numerical 
! differencing

subroutine test_jacobian() bind(C)

  use network
  use network_indices
  use eos_module
  use burner_module
  use actual_burner_module
  use rpar_indices
  use vode_module
  use vode_data
  use extern_probin_module, only: jacobian
  
  implicit none

  type (eos_t) :: state

  character (len=32) :: probin_file
  integer :: probin_pass(32)
  integer :: i, j

  integer          :: ipar
  double precision :: y(NEQ), ydot(NEQ), yp(NEQ), ym(NEQ), ydotp(NEQ), ydotm(NEQ), pd(NEQ,NEQ)
  double precision, allocatable :: rpar(:)

  double precision, parameter :: delta = 1.d-6
  double precision, parameter :: SMALL = 1.d-12
  double precision :: num_jac

  integer :: old_jacobian
  
  character(len=16) :: namei,namej  
  
  probin_file = "probin"
  do i = 1, len(trim(probin_file))
     probin_pass(i) = ichar(probin_file(i:i))
  enddo

  call runtime_init(probin_pass(1:len(trim(probin_file))), len(trim(probin_file)))

  call network_init()
  call burner_init()
  call eos_init()

  ! Allocate rpar. We need to do this after initializing the network
  ! since that is what sets up the rpar indices.

  allocate(rpar(n_rpar_comps))

  ! Set up EOS state
  
  state % rho       = 1.4311401611205835d7
  state % T         = 4.6993994016410122d9
  
  state % xn(ihe4)  = 4.2717633762309063d-3
  state % xn(ic12)  = 2.4502021307478711d-5
  state % xn(io16)  = 1.2059146851610723d-4
  state % xn(ine20) = 5.4419551339421394d-7
  state % xn(img24) = 2.5178594678377961d-4
  state % xn(isi28) = 3.5998829467937532d-1
  state % xn(is32)  = 2.7075529188304326d-1
  state % xn(iar36) = 9.1747472911892503d-2
  state % xn(ica40) = 8.0560189657331735d-2
  state % xn(iti44) = 6.1369127564250370d-4
  state % xn(icr48) = 2.5528582259065832d-3
  state % xn(ife52) = 1.9491916518179594d-2
  state % xn(ini56) = 1.6962109761781674d-1  
  
  call normalize_abundances(state)

  call eos(eos_input_rt, state)

  ! Set up the integration state

  call eos_to_vode(state, y, rpar)
  
  rpar(irp_self_heat) = ONE
  
  ! Evaluate the analytical Jacobian. Note that we
  ! need to call f_rhs first because that will fill rpar
  ! with the rates that the Jacobian needs.

  old_jacobian = jacobian
  jacobian = 1
  
  call f_rhs(NEQ, ZERO, y, ydot, rpar, ipar)  
  call jac(NEQ, ZERO, y, 0, 0, pd, NEQ, rpar, ipar)  
  
888 format(a,"-derivatives that don't match:")
999 format(5x, "df(",a,")/dy(",a,")", g18.10, g18.10, g18.10)

  ! Now evaluate a numerical estimate of the Jacobian
  ! using the RHS.

  jacobian = old_jacobian
  
  do j = 1, NEQ

     yp(:) = y(:)
     ym(:) = y(:)

     yp(j) = (ONE + delta) * y(j) 
     call f_rhs(NEQ, ZERO, yp, ydotp, rpar, ipar)
     
     ym(j) = (ONE - delta) * y(j) 
     call f_rhs(NEQ, ZERO, ym, ydotm, rpar, ipar)        

     if (j <= nspec) then
        namej = short_spec_names(j)
     else if (j==net_ienuc) then
        namej = "e"
     else if (j==net_itemp) then
        namej = "T"
     endif

     write(*,888) trim(namej)

     do i = 2, NEQ
        
        num_jac = (ydotp(i) - ydotm(i))/(yp(j) - ym(j) + SMALL)

        if (i <= nspec) then
           namei = short_spec_names(i)
        else if (i==net_ienuc) then
           namei = "e"
        else if (i==net_itemp) then
           namei = "T"
        endif

        ! only dump the ones that don't match
        if (num_jac /= pd(i,j)) then
           if (num_jac /= ZERO) then
              write (*,999) trim(namei), &
                   trim(namej), num_jac, pd(i,j), (num_jac-pd(i,j))/num_jac
           else
              write (*,999) trim(namei), &
                   trim(namej), num_jac, pd(i,j)
           endif
        endif

     enddo
  enddo  
  
end subroutine test_jacobian
