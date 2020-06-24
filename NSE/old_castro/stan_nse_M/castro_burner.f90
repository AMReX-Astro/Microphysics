module castro_burner_module

  use bl_types
  use bl_error_module
  use network
  use nse_module
  use burning_module

contains

  subroutine burner (dens, temp, Xin, ein, tau, dt, time, Xout, eout)

    use probdata_module, only :  bea_fuel, abar_fuel
    implicit none

    real(kind=dp_t), intent(in   ) :: dens,temp,Xin(nspec+naux),ein,tau,dt,time
    real(kind=dp_t), intent(  out) :: eout,Xout(nspec+naux)


    real(kind=dp_t) :: enuc, dyedt,ye,abar,dq,npalpha,siplusca,iron_group 
    real(kind=dp_t) :: t9 
    real(kind=dp_t) :: dXdt, xk, totmass
    real(kind=dp_t) :: dummy,denslog10,dx12
    real(kind=dp_t) :: deltaq,dbeadx
    real(kind=dp_t) :: bea_network,abar_network

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

    t9 = temp/1.d9
    denslog10 = log10(dens)

    !flame (including carbon, oxygen burning, silicon burning et al.)
    if (Xin(ic12) .ge. 0.01 .and. t9 .gt. 2.0d0 .and. denslog10 .gt. 6.5d0) then

!    if have substantial carbon at high t, assume this is a flame

       dXdt = -0.5d0/tau

!    burn out small residual carbon in zones that should be in nse

       if (Xin(ic12).lt.0.02d0.and.denslog10.gt.8.0d0) dXdt = -Xin(ic12)/dt

       call networkburn(denslog10,Xin(ic12),dummy,Xout(ihe4),dummy,&
                  Xout(io16),Xout(ine20),Xout(img24),Xout(isi28),&
                  Xout(ife56),bea_network,abar_network,dummy)

!    dbeadx is the total change in binding energy per nucleon when 
!    0.5 mass fraction of c12 burns (and a corresponding amount of o16)

       dbeadx = -(bea_network-bea_fuel)/0.5d0

!    apply a fraction of that change based upon how much c12 burns this step,
!    to get the energy generation rate. evolve c12 abundance. 
!    dXdt is negative and so too, usually is dbeadx

       enuc = 1.602d-6*6.023d23*dbeadx*dXdt
       eout = ein+ enuc*dt

       dx12 = dXdt*dt
       Xout(ic12) = Xin(ic12)+dx12
       Xout(ic12) = max(Xout(ic12),1.d-5)
       
!    evolve other species - sum of incoming mass fractions should be unity
!    so should the sum of mass fractions from the burn interpolation 

       Xout(ihe4)  = Xin(ihe4)  - Xout(ihe4) *(2.0d0*dx12)
       Xout(io16)  = Xin(io16)  - (Xout(io16) - 0.50d0)*(2.0d0*dx12)
       Xout(ine20) = Xin(ine20) - Xout(ine20)*(2.0d0*dx12)
       Xout(img24) = Xin(img24) - Xout(img24)*(2.0d0*dx12)
       Xout(isi28) = Xin(isi28) - Xout(isi28)*(2.0d0*dx12)
       Xout(ife56) = Xin(ife56) - Xout(ife56)*(2.0d0*dx12)


! protect against underflow to negative abundance

       Xout(ihe4)  = max(0.0d0,Xout(ihe4))
       Xout(io16)  = max(0.0d0,Xout(io16)) 
       Xout(ine20) = max(0.0d0,Xout(ine20))
       Xout(img24) = max(0.0d0,Xout(img24))
       Xout(isi28) = max(0.0d0,Xout(isi28))
       Xout(ife56) = max(0.0d0,Xout(ife56))

! make sure sume of mass fractions is 1. e.g., at low dens some of the 
! ashes are 12c even when the burning is over

       totmass = Xout(ihe4)+Xout(ic12)+Xout(io16)+Xout(ine20) &
               + Xout(img24)+Xout(isi28)+Xout(ife56)
       Xout(ihe4) = Xout(ihe4)/totmass
       Xout(ic12) = Xout(ic12)/totmass
       Xout(io16) = Xout(io16)/totmass
       Xout(ine20) = Xout(ine20)/totmass
       Xout(img24) = Xout(img24)/totmass
       Xout(isi28) = Xout(isi28)/totmass
       Xout(ife56) = Xout(ife56)/totmass       

! generate abar from actual composition. reciprocal sof sum of y's

       totmass = Xout(ihe4)/4.0d0+Xout(ic12)/12.0d0+Xout(io16)/16.0d0 &
         +Xout(ine20)/20.0d0+Xout(img24)/24.0d0+Xout(isi28)/28.0d0    &
         +Xout(ife56)/56.0d0
       Xout(iabar)= 1.0d0/totmass

! this is approximate. base current be/a upon current mass fraction
! and burning at the present density. if carbon is mixed into
! ashes of burning that happened at a higher density, this won't
! be exactly right. it is not clear that this quantity is used anywhere
! execpt in the transition to nse

       Xout(ibea) = bea_fuel+((0.5d0-Xout(ic12))/0.5d0)*(bea_network-bea_fuel)

! another option would be this, but bea could inprinciple grow too large
! from mixing carbon into a zone that had already burned

!       Xout(ibea) = Xin(ibea)+dbeadx*dXdt*dt

! no weak interactions

       Xout(iye)  = Xin(iye)

!   check for nse. the density and temperature must both be high as 
!   well as the iron abundace. otherwise this could turn e.g. si
!   zones artificially into fe zones
  
    else if(((Xin(ife56)+Xin(ihe4)).gt.0.88d0.or.denslog10.gt. 8.5d0) .and. &
                     Xin(ic12).lt.0.01d0.and.temp.gt.3.d9.and.denslog10.gt.8.0d0) then
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
       Xout(:) = Xin(:)
    endif
 
  end subroutine burner

end module castro_burner_module
