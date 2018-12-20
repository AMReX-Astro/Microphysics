program evaluate_rhs

  use amrex_fort_module, only : rt => amrex_real

  use actual_network
  use actual_rhs_module
  use reaclib_rates, only: init_reaclib, net_screening_init
  use table_rates, only: init_tabular
  use eos_type_module
  use eos_module
  use burn_type_module

  type(burn_t) :: state
  type(eos_t)  :: eos_state
  integer      :: i, j
  real(rt)   :: ynum, yden
  character(len=20) :: FMT

  ! Initialize Network
  call actual_network_init()
  call init_reaclib()
  call init_tabular()
  call net_screening_init()

  ! Initialize EOS
  call eos_init()
  
  ! Fill burn state
  write(*,*) 'State Density (g/cm^3): '
  read(*,*) state%rho
  write(*,*) 'State Temperature (K): '
  read(*,*) state%T
  do i = 1, nspec
     write(*,*) 'Mass Fraction (', spec_names(i), '): '
     read(*,*) state%xn(i)
  end do

  ! Call EOS to get thermo variables, abar, zbar, ye, etc.
  call burn_to_eos(state, eos_state)
  call eos(eos_input_rt, eos_state)
  call eos_to_burn(eos_state, state)
  
  ! Evaluate RHS
  call actual_rhs(state)
  
  ! Evaluate Jacobian
  call actual_jac(state)

  write(*,'(A25I25)') 'nspec_evolve: ', nspec_evolve
  write(*,'(A25E30.16E5)') 'Density (g/cm^3): ', state%rho
  write(*,'(A25E30.16E5)') 'Temperature (K): ', state%T
  write(*,'(A25E30.16E5)') 'Ye: ', state%y_e
  write(*,*) 'RHS Evaluation'
  ! Print RHS
  do i = 1, nspec_evolve
     write(*,'(A5A5A3E30.16E5)') 'ydot(', short_spec_names(i), '): ', state%ydot(i)
     write(*,'(A5A5A3E30.16E5)') 'xdot(', short_spec_names(i), '): ', state%ydot(i)*aion(i)
  end do
  write(*,'(A10E30.16E5)') 'dot temp: ', state%ydot(net_itemp)
  write(*,'(A10E30.16E5)') 'dot enuc: ', state%ydot(net_ienuc)
  write(*,*) '--------------------'
  write(*,*) 'Jacobian Evaluation: d(dYi/dt)/dYj'
  ! Print Jacobian
  write(FMT, '("(", I0, "E30.16E5)")') net_ienuc
  do i = 1, net_ienuc
     write(*,FMT) ( state%jac(i, j), j = 1, net_ienuc )
  end do
  write(*,*) '--------------------'
  write(*,*) 'd(dY(1:nspec_evolve)/dt)/dY(1:nspec_evolve)'
  write(FMT, '("(", I0, "E30.16E5)")') nspec_evolve
  do i = 1, nspec_evolve
     write(*,FMT) ( state%jac(i, j), j = 1, nspec_evolve )
  end do
  
  write(*,*) 'd(dY(net_ienuc)/dt)/dY(1:nspec_evolve)'
  write(*,FMT) ( state%jac(net_ienuc, j), j = 1, nspec_evolve )
  
  write(*,*) 'd(dY(1:nspec_evolve)/dt)/dY(net_ienuc)'
  write(*,FMT) ( state%jac(j, net_ienuc), j = 1, nspec_evolve )
  
  write(*,*) 'd(dY(net_itemp)/dt)/dY(1:nspec_evolve)'
  write(*,FMT) ( state%jac(net_itemp, j), j = 1, nspec_evolve )
  
  write(*,*) 'd(dY(1:nspec_evolve)/dt)/dY(net_itemp)'
  write(*,FMT) ( state%jac(j, net_itemp), j = 1, nspec_evolve )
  
  write(FMT, '("(", I0, "E30.16E5)")') 1
  
  write(*,*) 'd(dY(net_ienuc)/dt)/dY(net_ienuc)'
  write(*,FMT) state%jac(net_ienuc, net_ienuc)
  
  write(*,*) 'd(dY(net_ienuc)/dt)/dY(net_itemp)'
  write(*,FMT) state%jac(net_ienuc, net_itemp)

  write(*,*) 'd(dY(net_itemp)/dt)/dY(net_itemp)'
  write(*,FMT) state%jac(net_itemp, net_itemp)
  
  write(*,*) 'd(dY(net_itemp)/dt)/dY(net_ienuc)'
  write(*,FMT) state%jac(net_itemp, net_ienuc)
  
end program evaluate_rhs
