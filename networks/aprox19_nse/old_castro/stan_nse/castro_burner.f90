module castro_burner_module

  use bl_types
  use bl_error_module
  use network
  use nse_module
  use burning_module
  use interpolate_module

contains

  subroutine burner (dens, temp, Xin, ein, tau, dt, time, Xout, eout)

    use probdata_module, only : taudiv, bea_fuel, abar_fuel
    implicit none

    real(kind=dp_t), intent(in   ) :: dens,temp,Xin(nspec),ein,tau,dt,time
    real(kind=dp_t), intent(  out) :: eout,Xout(nspec)


    real(kind=dp_t) :: enuc, dyedt,ye,abar,dq,npalpha,siplusca,iron_group 
    real(kind=dp_t) :: t9 
    real(kind=dp_t) :: dXdt, xk
    real(kind=dp_t) :: dummy,rho_fuel
    real(kind=dp_t) :: deltaq,dbeadx
    real(kind=dp_t) :: bea_network,abar_network
    real(kind=dp_t) :: tauo5

    integer, save :: ic12,io16,ihe4,ine20,img24,isi28,ife56,iye,ibea, iabar
    
    logical, save :: firstCall = .true.

    integer :: i

    !$OMP CRITICAL (stan_burner_init)
    if (firstCall) then
       if (.NOT. network_initialized) then
          call bl_error("ERROR in burner: must initialize network first")
       endif
       ic12  = network_species_index("carbon-12")
       io16  = network_species_index("oxygen-16")
       ihe4  = network_species_index("helium-4")
       ine20 = network_species_index("neon-20")
       img24 = network_species_index("magnesium-24")
       isi28 = network_species_index("silicon-28")
       ife56 = network_species_index("iron-56")
       iye   = network_species_index("Ye")
       ibea  = network_species_index("BeA")
       iabar = network_species_index("Abar")

       if (ic12 < 0) then
          call bl_error("ERROR in burner: species undefined")
       endif
       firstCall = .false.
    endif
    !$OMP END CRITICAL (stan_burner_init)

    call networkburn(log10(dens),Xin(ic12),rho_fuel,dummy,dummy,&
                     dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy)
    t9 = temp/1.d9

    tauo5 = tau

    !flame (including carbon, oxygen burning, silicon burning et al.)
    if (Xin(ic12) .ge. 0.01 .and. t9 .gt. 2.0d0 .and. rho_fuel .gt. 6.5d0) then

       if(t9 .gt. 3.0d0) tauo5 = tau/taudiv

       if (Xin(ic12) .ge. 0.02d0) then
                dXdt = -0.5d0/tauo5
       else
!                dXdt = -min(0.5d0/tauo5,Xin(ic12)/dt)
                dXdt = -Xin(ic12)/dt
       endif

       call networkburn(log10(dens),Xin(ic12),dummy,Xout(ihe4),dummy,&
                  Xout(io16),Xout(ine20),Xout(img24),Xout(isi28),&
                  Xout(ife56),bea_network,abar_network,dummy)

       dbeadx = -(bea_network-bea_fuel)/0.5d0
       enuc = 1.602d-6*6.023d23*dbeadx*dXdt
!       enuc = max(0.0d0,enuc)
       eout = ein+ enuc*dt

       Xout(ic12) = Xin(ic12)+dXdt*dt
       Xout(ic12) = max(Xout(ic12),1.d-5)
       Xout(iabar)= abar_fuel+((0.5d0 - Xout(ic12))/0.5d0)*(abar_network - abar_fuel)   
       Xout(ibea) = bea_fuel +((0.5d0 - Xout(ic12))/0.5d0)*(bea_network - bea_fuel)    
       Xout(iye)  = Xin(iye)

    ! nse  
    else if((Xin(ife56)+Xin(ihe4) .gt. 0.88d0 .or. log10(dens) .gt. 8.5d0) .and. &
                     Xin(ic12) .lt. 0.01d0 .and. temp .gt. 3.d9) then
       ye = Xin(iye)
       call interp(t9,dens,ye,abar,dq,dyedt,npalpha,siplusca,iron_group)
       deltaq = dq-Xin(ibea)
       deltaq = 0.3d0*deltaq           !underrelaxation
       Xout(ibea) = Xin(ibea) + deltaq
       enuc = deltaq*1.602d-6*6.023d23/dt
       eout = ein+ enuc*dt

       Xout(iye)   = Xin(iye)-dt*dyedt
       Xout(iabar) = abar
       if(Xin(ic12) .lt. 1.d-3) then
          Xout(ic12) = Xin(ic12)
       else 
          Xout(ic12)  = 1.d-3
       endif
       Xout(io16)  = 1.d-20
       Xout(ihe4)  = npalpha
       Xout(ine20) = 1.d-20
       Xout(img24) = 1.d-20
       Xout(isi28) = siplusca
       Xout(ife56) = iron_group

    else
       enuc= 0.d0
       eout = ein
       do i =1, nspec 
         Xout(i) = Xin(i)
       enddo
    endif
 
  end subroutine burner

end module castro_burner_module
