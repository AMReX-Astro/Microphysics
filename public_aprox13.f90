      program drive_aprox13
      include 'implno.dek'
      include 'const.dek'
      include 'timers.dek'
      include 'burn_common.dek'
      include 'network.dek'

! this program exercises the aprox13 network

! declare
      character*40     summary
      integer          i,j,nok,nbad

      double precision tstart,tstep,conserv,tin,din,ein,vin,zin,xin(18), &
                       tout,dout,eout,xout(18),edum



! initialize the network and eos
      call init_aprox13
      call read_helm_table


! keep coming back to here
20    continue

      call net_input(tstart,tstep,tin,din,vin,zin,ein,xin)


! start the clock
      call zsecond(timzer)


! a message
        write(6,*)
        write(6,*) 'starting integration'
        write(6,*)

! burn it

        call burner(tstart,tstep, &
                    tin,din,vin,zin,ein,xin, &
                    tout,dout,eout,xout, &
                    conserv,nok,nbad)


!      edum = ev2erg*1.0d6*avo * sum((xout(1:ionmax) - xin(1:ionmax))/aion(1:ionmax)*bion(1:ionmax))
!      write(6,123) edum
!      write(6,123) edum - sneut*tstep
!123   format(1x,1pe14.6)



! output a summary of the integration

       call net_summary(tstep,tin,din,ein, &
                        tout,dout,eout,conserv, &
                        nbad,nok,xout)


! back for another input point
      goto 20
      end








      subroutine burner(beg,tstep, &
                        tin,din,vin,zin,ein,xin, &
                        tout,dout,eout,xout, &
                        conserv,nok,nbad)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'



! input:
! beg       starting time
! tstep     time over which to integrate
! tin       initial temperature
! din       initial density
! vin       initial velocity                                                                                           
! zin       initial position
! ein       initial internal energy
! xin(1:19) initial composition vector

! output:
! tout       final temperature
! dout       final density
! eout       final internal energy
! xout(1:19) final composition vector
! conserv    1 - sum of mass fractions
! nok        number of good time steps taken
! nbad       number of bad timesteps attempted


! declare the pass
      integer          nok,nbad
      double precision beg,tstep,tin,din,vin,zin,ein,xin(*), &
                       tout,dout,eout,xout(*),conserv


! local variables
      integer          i
      double precision abar,zbar,wbar,ye,xcess

! for the integration driver
      integer          kount
      double precision stptry,stpmin,tend,ys2(abignet*nzmax), &
                       odescal,tol

! usually adequate
      parameter        (tol     = 1.0d-6, &
                        odescal = 1.0d-8)

! for very accurate integrations
!      parameter        (tol     = 1.0d-10, &
!                        odescal = 1.0d-12)


      external         aprox13,saprox13,baprox13,daprox13


!      external         forder_ma28
!      external         forder_umf
!      external         forder_y12m
!      external         forder_ludcmp
!      external         forder_leqs
!      external         forder_lapack
!      external         forder_gift
!      external         forder_biconj

!      external         rosen_ma28
!      external         rosen_umf
!      external         rosen_y12m
!      external         rosen_ludcmp
!      external         rosen_leqs
!      external         rosen_lapack
!      external         rosen_gift
!      external         rosen_biconj

      external         stifbs_ma28
!      external         stifbs_umf
!      external         stifbs_y12m
!      external         stifbs_ludcmp
!      external         stifbs_leqs
!      external         stifbs_lapack
!      external         stifbs_gift
!      external         stifbs_biconj





! set the initial condition


! load the mass fractions
       xmass(ionbeg:ionend) = xin(ionbeg:ionend)


! get abar, zbar and a few other composition variables
       call azbar(xmass(ionbeg),aion(ionbeg),zion(ionbeg),wion(ionbeg),ionmax, &
                  ymass(ionbeg),abar,zbar,wbar,ye,xcess)


! stuff the initial conditions into ys2
       ys2(ionbeg:ionend) = ymass(ionbeg:ionend)
       ys2(iener) = ein
       ys2(itemp) = tin
       ys2(iden)  = din
       ys2(ivelx) = vin
       ys2(iposx) = zin


! single step (tend=tstep), hydrostatic, or expansion ending times.
! the variable tstep has two meanings here. tstep in single step mode
! is the size of the time step to try. tstep in hydrostatic or expansion
! mode is the ending integration time. the integration driver really
! gets some exercise if tstep is large in single step mode.

      tend = tstep
      if (one_step) then
       stptry = tstep
       stpmin = tstep * 1.0d-20
      else
       stptry = max(beg * 1.0d-10,1.0d-16)
       stpmin = stptry * 1.0d-12
      end if



! integrate the aprox13 network
       call netint(beg,stptry,stpmin,tend,ys2, &
                   tol,neqs,nok,nbad,kount,odescal, &
!                  aprox13,saprox13,baprox13,forder_ma28)
!                  aprox13,saprox13,baprox13,forder_umf)
!                  aprox13,saprox13,baprox13,forder_y12m)
!                  aprox13,daprox13,baprox13,forder_ludcmp)
!                  aprox13,daprox13,baprox13,forder_leqs)
!                  aprox13,daprox13,baprox13,forder_lapack)
!                  aprox13,daprox13,baprox13,forder_gift)
!                  aprox13,saprox13,baprox13,forder_biconj)
!                  aprox13,saprox13,baprox13,rosen_ma28)
!                  aprox13,saprox13,baprox13,rosen_umf)
!                  aprox13,saprox13,baprox13,rosen_y12m)
!                  aprox13,daprox13,baprox13,rosen_ludcmp)
!                  aprox13,daprox13,baprox13,rosen_leqs)
!                  aprox13,daprox13,baprox13,rosen_lapack)
!                  aprox13,daprox13,baprox13,rosen_gift)
!                  aprox13,saprox13,baprox13,rosen_biconj)
                  aprox13,saprox13,baprox13,stifbs_ma28)
!                  aprox13,saprox13,baprox13,stifbs_umf)
!                  aprox13,saprox13,baprox13,stifbs_y12m)
!                  aprox13,daprox13,baprox13,stifbs_ludcmp)
!                  aprox13,daprox13,baprox13,stifbs_leqs)
!                  aprox13,daprox13,baprox13,stifbs_lapack)
!                  aprox13,daprox13,baprox13,stifbs_gift)
!                  aprox13,saprox13,baprox13,stifbs_biconj)




! set the output
      do i=ionbeg,ionend
       xout(i) = ys2(i) * aion(i)
      enddo
      tout = ys2(itemp)
      dout = ys2(iden)
      eout = ys2(iener)

      conserv = 0.0d0
      do i=ionbeg,ionend
       conserv = conserv + xout(i)
      enddo
      conserv = 1.0d0 - conserv

      return
      end



!---------------------------------------------------------------------







!---------------------------------------------------------------------
! this file contains aprox13 network

! routine aprox13 sets up the odes
! routine rhs evaluates the right hand sides
! routine daprox13 sets up the dense aprox13 jacobian
! routine baprox13 builds the nonzero locations for saprox13
! routine saprox13 sets up the sparse aprox13 jacobian
! routine aprox13rat generates the reaction rates for routine aprox13
! routine aprox13tab generates the raw rates using table interpolation
! routine screen_aprox13 applies screening corrections to the raw rates
! routine init_aprox13 initializes the aprox13 network





!---------------------------------------------------------------------




!---------------------------------------------------------------------
      subroutine net_input(tstart,tstep,tin,din,vin,zin,ein,xin)
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
      include 'burn_common.dek'
      include 'network.dek'
      include 'cjdet.dek'


! declare the pass
      double precision tstart,tstep,tin,din,vin,zin,ein,xin(*)


! local variables
      character*80     string,word
      integer          i,j,k,ibtype,ictype,igues,kkase,ians,getnam
      double precision xneut,xh1,xhe4,xc12,xc13,xn14,xo16,xne20,xne22, &
                       xsi28,xfe52,xfe54,xfe56,xni56,zye,sum,abar,zbar, &
                       wbar,xcess,ye,ye_orig,xmup,xmun,qdum,a,z,xelem, &
                       andgrev,value


! bigbang specifics
      double precision fac,f1,zeta3
      parameter        (zeta3 = 1.20205690315732d0)


! popular format statements
01    format(1x,a,a,a)
02    format(1x,a,'=',1pe10.3,' ',a,'=',1pe10.3,' ', &
                a,'=',1pe10.3,' ',a,'=',1pe10.3,' ', &
                a,'=',1pe10.3)
 03   format(a)
 04   format(1x,a,'=',i2,' ',a,'=',i2,' ', &
                a,'=',i2,' ',a,'=',i2,' ', &
                a,'=',i2)



! initialize the common block variables
      call net_initialize


! inititailize local variables
      ibtype    = 0
      ictype    = 0
      tstart    = 0.0d0
      tstep     = 0.0d0
      bpres     = 0.0d0
      tin       = 0.0d0
      din       = 0.0d0
      vin       = 0.0d0
      zin       = 0.0d0
      zye       = 0.0d0
      xin(1:ionmax) = 1.0d-30


!---------------------------------------------------------------------------



! get the burn type
 10   write(6,01) 'give burning mode:'
      write(6,01) '     ibtype = 0 = stop'
      write(6,01) '              1 = onestep'
      write(6,01) '              2 = hydrostatic'
      write(6,01) '              3 = expansion'
      write(6,01) '              4 = self-heat at constant density'
      write(6,01) '              5 = self heat at constant pressure'
      write(6,01) '              6 = self-heat pressure-temp trajectory'
      write(6,01) '              7 = big bang '
      write(6,01) '              8 = detonation'
      write(6,01) '              9 = temp-den trajectory'

      read(5,*)  ibtype
      if (ibtype .lt. 0 .or. ibtype .gt. 9) goto 10

! set the burn type logical
      if (ibtype .eq. 0) then
       stop 'normal termination'
      else if (ibtype .eq. 1) then
       one_step = .true.
      else if (ibtype .eq. 2) then
       hydrostatic = .true.
      else if (ibtype .eq. 3) then
       expansion = .true.
      else if (ibtype .eq. 4) then
       self_heat_const_den = .true.
      else if (ibtype .eq. 5) then
       self_heat_const_pres = .true.
      else if (ibtype .eq. 6) then
       pt_hist = .true.
      else if (ibtype .eq. 7) then
       bbang = .true.
      else if (ibtype .eq. 8) then
       detonation = .true.
      else if (ibtype .eq. 9) then
       trho_hist = .true.
      else
       goto 10
      end if





! general options
 11   write(6,*)
      write(6,04) 'set general options:'
      write(6,04) 'screen_on',screen_on
      write(6,04) 'use_tables',use_tables
      write(6,04) 'weak_on',weak_on
      write(6,04) 'ffn_on',ffn_on
      write(6,04) 'pure_network',pure_network
      write(6,04) 'nse_analysis',nse_analysis
      write(6,04) 'allow_nse_evol',allow_nse_evol
      write(6,04) 'iprint_files',iprint_files
      write(6,04) 'iprint_screen',iprint_screen
      write(6,02) 'sthreshold',sthreshold,' set > 1 to disable'
      write(6,*)
      write(6,01) 'if these are ok, enter 1, otherwise enter 0 =>'

      read(5,*) ians
      if (ians .lt. 0 .or. ians .gt. 1) goto 11

      if (ians .eq. 0) then
 12    write(6,01) 'give the 9 integer and one real vector =>'

       read(5,*) screen_on, use_tables, weak_on, ffn_on, &
                 pure_network, nse_analysis, allow_nse_evol, &
                 iprint_files, iprint_screen, &
                 sthreshold

       if (screen_on .lt. 0 .or. screen_on .gt. 1) goto 12
       if (use_tables .lt. 0 .or. use_tables .gt. 1) goto 12
       if (weak_on .lt. 0 .or. weak_on .gt. 1) goto 12
       if (ffn_on .lt. 0 .or. ffn_on .gt. 1) goto 12
       if (pure_network .lt. 0 .or. pure_network .gt. 1) goto 12
       if (nse_analysis .lt. 0 .or. nse_analysis .gt. 1) goto 12
       if (iprint_files .lt. 0 .or. iprint_files .gt. 1) goto 12
       if (iprint_screen .lt. 0 .or. iprint_screen .gt. 1) goto 12
       goto 11
      end if



! get the bigbang parameters; set default to wmap 2008 (5 year) values
      if (bbang) then

       eta1    = 6.23e-10
       xnnu    = 3.0d0
       hubble  = 70.5d0
       cmbtemp = 2.725d0

 13    write(6,*)
       write(6,02) 'bigbang parameters:'
       write(6,02) 'eta',eta1
       write(6,02) 'number of neutrino families',xnnu
       write(6,02) 'hubble constant',hubble
       write(6,02) 'present cmb temperature',cmbtemp

       write(6,01) 'if these are ok, enter 1, otherwise enter 0 =>'
       read(5,*) ians
       if (ians .lt. 0 .or. ians .gt. 1) goto 13

       if (ians .eq. 0) then
        write(6,01) 'give eta, xnu, hubble, and cmbtemp  =>'
        read(5,*) eta1, xnnu, hubble, cmbtemp
        goto 13
       end if
      end if



! get an alternative the stopping condition; when the
! mass fraction of a given isotope falls below a given level


 14   write(6,*)
      write(6,*) 'stop when an isotope falls below a given abundance?', &
                  ' 1=yes 0=no'
      read(5,*)  ians
      if (ians .lt. 0 .or. ians .gt. 1) goto 14

      if (ians .eq. 0) then
       name_stop = 'he4 '
       xmass_stop = -1.0d30
      end if

 15   if (ians .eq. 1) then
       write(6,*) 'give the name of the isotope and the mass fraction'
       write(6,*) 'for example: c12 0.50'

       read(5,03) string
       j = 1
       i = getnam(string,word,j)
       name_stop = word(1:5)

       i = getnam(string,word,j)
       xmass_stop = value(word)

       write(6,*) name_stop,xmass_stop
      end if


! check that the name_stop isotope is in the network
      do i=1,ionmax
       if (ionam(i) .eq. name_stop) then
        id_stop = i
        goto 16
       end if
      enddo
      write(6,*)
      write(6,*) 'name_stop>',name_stop,'< not in network'
      write(6,*)
      if (ians .eq. 1) goto 15
      stop ' bad name for stopping isotope'
 16   continue




! get the initial thermodynamics
      write(6,*)
      if (self_heat_const_pres) then
       write(6,01) 'give the ending time, temperature, pressure =>'
       read(5,*)  tstep,tin,bpres

      else if (bbang) then
       write(6,01) 'give the ending time, initial temperature =>'
       read(5,*)  tstep,tin

      else if (.not. (trho_hist .or. pt_hist)) then
       write(6,01) 'give the ending time, temperature, density =>'
       read(5,*)  tstep,tin,din
      end if

! limit the temperature since the rates are invalid much above t9=100
       tin = min(1.0d11,tin)



! get the composition
      if (.not. bbang) then
 20    write(6,01) 'give initial composition:'
       write(6,01) '     ictype = 0 = leave alone; read from file'
       write(6,01) '              1 = solar abundances'
       write(6,01) '              2 = nse'
       write(6,01) '              3 = specify initial composition'

       read(5,*) ictype
       if (ictype .lt. 0 .or. ictype .gt. 3) goto 20

       if (ictype .eq. 3) then
        write(6,01) &
        'n h1 he4 c12 c13 n14 o16 ne20 ne22 si28 fe52 fe54 fe56 ni56 =>'
        read(5,*) xneut,xh1,xhe4,xc12,xc13,xn14,xo16,xne20,xne22,xsi28, &
                  xfe52,xfe54,xfe56,xni56
       end if
      end if


! get the output root file name
      write(6,*)  ' '
      write(6,01) 'give output root name, <cr> for default "foo_"=>'
      read(5,03) hfile
      if (hfile(1:2) .eq. '  ')  hfile = 'foo_'


!---------------------------------------------------------------------------



! set some more variables based on the burn type

! adiabatic expansion
! psi =  1 is an adiabatic expansion, -1 in an adiabatic implosion

      if (expansion) then
       psi       = 1.0d0
!       psi       = -1.0d0
       den0      = din
       temp0     = tin
       temp_stop = 1.0d7
!       temp_stop = 1.0d10
       if ( (psi .ge. 1.0  .and. temp_stop .ge. tin)  .or. &
            (psi .le. -1.0 .and. temp_stop .le. tin)) &
          stop 'bad adiabatic temp_stop in routine burner'



! big bang
      else if (bbang) then

! set the initial n and p abundances; equation 3 of wagoner et al 1967
       fac = exp((mn - mp)*clight**2/(kerg*tin))
       xneut = 1.0d0/(1.0d0 + fac)
       xh1   = 1.0d0 - xneut

! set the density from the temperature and eta1
       f1  = 30.0d0 * zeta3/pi**4 * asol/(kerg*avo)
       din = f1 * eta1 * tin**3


! thermodynamic profile being given
      else if (trho_hist .or. pt_hist) then
       write(6,*) 'give the trajectory file =>'
       read(5,03) trho_file
      end if


!---------------------------------------------------------------------------



! read the thermodynamic trajectory and initial abundances
! transfer the info stored in xsum and zsum from the update2 call

      if (trho_hist) then
       call update2(tstart,tin,din)
       xin(1:ionmax) = xsum(1:ionmax)
       tstart     = zwork1(1)
       tstep      = zwork1(2)
       zye        = zwork1(3)
      end if


      if (pt_hist) then
       call update3(tstart,tin,bpres)
       xin(1:ionmax) = xsum(1:ionmax)
       tstart     = zwork1(1)
       tstep      = zwork1(2)
       zye        = zwork1(3)
      end if



! massage the input composition, includes possible changes to the
! the abundances read in from the trho_hist file

! solar abundances
      if (ictype .eq. 1) then
       do i=1,ionmax
        xin(i) = andgrev(ionam(i),z,a,xelem)
       enddo
       if (iprot .ne. 0) xin(iprot) = andgrev('h1   ',z,a,xelem)


! put it in nse
      else if (ictype .eq. 2) then
       if (zye .eq. 0.0) zye   = 0.5d0
       igues = 1
       call nse(tin,din,zye,igues,1,1,xin,xmun,xmup,0)


! set the composition variables
      else if (ictype .eq. 3 .or. bbang) then
       if (ineut .ne. 0) xin(ineut) = xneut
       if (ih1   .ne. 0) xin(ih1)   = xh1
       if (iprot .ne. 0) xin(iprot) = xh1
       if (ih1 .ne. 0 .and. iprot .ne. 0) xin(iprot) = 0.0d0
       if (ihe4  .ne. 0) xin(ihe4)  = xhe4
       if (ic12  .ne. 0) xin(ic12)  = xc12
       if (ic13  .ne. 0) xin(ic13)  = xc13
       if (in14  .ne. 0) xin(in14)  = xn14
       if (io16  .ne. 0) xin(io16)  = xo16
       if (ine20 .ne. 0) xin(ine20) = xne20
       if (ine22 .ne. 0) xin(ine22) = xne22
       if (isi28 .ne. 0) xin(isi28) = xsi28
       if (ife52 .ne. 0) xin(ife52) = xfe52
       if (ife54 .ne. 0) xin(ife54) = xfe54
       if (ife56 .ne. 0) xin(ife56) = xfe56
       if (ini56 .ne. 0) xin(ini56) = xni56


! hardcode something here

!if (ih1 .ne. 0)   xin(ih1)=       7.0572558936810803E-01
!if (ih2 .ne. 0)   xin(ih2)=       4.8010000000000003E-05
!if (ihe3 .ne. 0)  xin(ihe3)=      2.9291000000000001E-05
!if (ihe4 .ne. 0)  xin(ihe4)=      2.7521000000000001E-01
!if (ili7 .ne. 0)  xin(ili7)=      9.3489999999999999E-09
!if (ic12 .ne. 0)  xin(ic12)=      3.0324000000000002E-03
!if (ic13 .ne. 0)  xin(ic13)=      3.6501000000000002E-05
!if (in14 .ne. 0)  xin(in14)=      1.1049000000000000E-03
!if (in15 .ne. 0)  xin(in15)=      4.3633999999999996E-06
!if (io16 .ne. 0)  xin(io16)=      9.5917999999999993E-03
!if (io17 .ne. 0)  xin(io17)=      3.8873000000000000E-06
!if (io18 .ne. 0)  xin(io18)=      2.1673000000000001E-05
!if (if19 .ne. 0)  xin(if19)=      4.0515000000000000E-07
!if (ine20 .ne. 0) xin(ine20)=     1.6188999999999999E-03
!if (img24 .ne. 0) xin(img24)=     3.5722704328918775E-03


!if (ihe4 .ne. 0)  xin(ihe4)=     5.45516e-07 
!if (ic12 .ne. 0)  xin(ic12)=     0.491254d0
!if (io16 .ne. 0)  xin(io16)=      0.494312d0
!if (ine20 .ne. 0) xin(ine20)=     0.0142452d0
!if (img24 .ne. 0) xin(img24)=     0.000187270d0
!if (isi28 .ne. 0) xin(isi28)=     9.08493e-07
!if (is32 .ne. 0) xin(is32)=      4.30614e-11
!if (iar36 .ne. 0) xin(iar36)=     5.21367e-16
!if (ica40 .ne. 0) xin(ica40)=     1.06910e-20
!if (iti44 .ne. 0) xin(iti44)=     1.00000e-20
!if (icr48 .ne. 0) xin(icr48)=     1.00000e-20
!if (ife52 .ne. 0) xin(ife52)=     1.00000e-20
!if (ini56 .ne. 0) xin(ini56)=     1.00000e-20


      end if


! write out the input composition so far
!      write(6,02) (ionam(i),xin(i), i=1,ionmax)
!      read(5,*)


! normalize the composition
      do i=1,ionmax
       xin(i) = min(1.0d0,max(xin(i),1.0d-30))
      end do
      sum = 0.0d0
       do i=1,ionmax
        sum = sum + xin(i)
       enddo
      sum = 1.0d0/sum
      do i=1,ionmax
       xin(i) = min(1.0d0,max(xin(i) * sum,1.0d-30))
      enddo

!      write(6,*) 'post norm', ionmax,xin(ih1)
!      write(6,02) (ionam(i),xin(i), i=1,ionmax)
!      read(5,*)

!---------------------------------------------------------------------------


! get the ye of the initial compositon
        call azbar(xin,aion,zion,wion,ionmax, &
                   zwork1,abar,zbar,wbar,ye_orig,xcess)

!       write(6,123) abar,zbar
!123    format(1x,1p2e12.3)
!       read(5,*)



! modify the composition if ye_orig is less than 0.55
!        if (ye_orig .le. 0.55) then
!
! set the mass fraction of fe58 to set the desired ye
!         ye_want = 0.495d0
!         ye_want = 0.50d0
!         if (ye_want .eq. 0.5) then
!          xin(ife58) = 0.0d0
!         else
!          xin(ife58) = (ye_orig - ye_want) /
!     1                  (ye_orig - zion(ife58)/aion(ife58))
!         end if
!
! reset the mass fractions of everything else
!         sum = 1.0d0 - xin(ife58)
!         do i=1,ionmax
!          if (i .ne. ife58) xin(i) = xin(i) * sum
!         enddo
!        end if


!---------------------------------------------------------------------------


! modify for a detonation
! get the chapman-jouget solution
       if (detonation) then
        kkase = 1
         mach  = 0.0d0
        do i=1,ionmax
         xmass_up(i) = xin(i)
        enddo
        temp_up = tin
        den_up  = din
        call cjsolve(kkase,xmass_up,temp_up,den_up,mach, &
                    qburn_cj,xmass_cj,ener_up,pres_up,cs_up, &
                    vel_det,vel_cj,temp_cj,den_cj,ener_cj,pres_cj,cs_cj)


        write(6,*)  ' '
        write(6,63) 'cj state (should be sonic with vel_mat = cs_cj):'
        write(6,61) 'temp_cj',temp_cj,'den_cj ',den_cj, &
                    'pres_cj',pres_cj
        write(6,61) 'cs_cj  ',cs_cj, &
                    'vel_mat',vel_cj,'vel_det',vel_det
        write(6,61) 'mach_cj',vel_cj/cs_cj,'qburn_cj',qburn_cj

 63     format(1x,a)
 61     format(1x,a7,'=',1pe10.3,' ',a7,'=',1pe10.3,' ', &
                 a7,'=',1pe10.3,' ',a4,'=',1pe10.3)


        write(6,*) ' '
        write(6,*) 'top 10 cj nse mass fractions:'
        call indexx(ionmax,xmass_cj,izwork1)
        write(6,02) (ionam(izwork1(i)), &
                   xmass_cj(izwork1(i)), i=ionmax,ionmax-9,-1)


! get shock solution
        kkase = 4
        mach_sh = vel_det/cs_up
        call cjsolve(kkase,xmass_up,temp_up,den_up,mach_sh, &
                    qdum,xmass_up,ener_up,pres_up,cs_up, &
                    vel_det,vel_sh,temp_sh,den_sh,ener_sh,pres_sh,cs_sh)


! reset the initial conditions for znd detonations
        tin      = temp_sh
        din      = den_sh
        vin      = vel_sh
        zin      = 1.0e-16*vel_sh
        den_stop = 1.00d0 * den_cj

        write(6,*)
        write(6,*) 'resetting initial conditions for a detonation to:'
        write(6,64) 'tin=',tin,' din=',din,' vin=',vin,' zin=',zin
 64     format(1x,4(a,1pe12.4) )
       end if


!---------------------------------------------------------------------------


! get the abundance variables for the final mixture
        call azbar(xin,aion,zion,wion,ionmax, &
                   zwork1,abar,zbar,wbar,ye_orig,xcess)


! get the thermodynamic state
      temp_row(1) = tin
      den_row(1)  = din
      ptot_row(1) = bpres
      abar_row(1) = abar
      zbar_row(1) = zbar
      jlo_eos = 1
      jhi_eos = 1

!      write(6,*) tin,abar,zbar

      if (self_heat_const_pres .or. pt_hist) then
       den_row(1)  = bpres * abar/(avo * kerg * tin)
       call invert_helm_pt
       din = den_row(1)

!       write(6,778) bpres,din
!       read(5,*)

      else
       call helmeos
       bpres = ptot_row(1)
      endif

      ein   = etot_row(1)


!---------------------------------------------------------------------------


! write out the final input
        write(6,*)
        write(6,02) 'tstart',tstart,'tstep',tstep
        write(6,02) 'tin',tin,'din',din,'bpres',bpres,'ein',ein

! largest mass fractions
        call indexx(ionmax,xin,izwork1)
        j = min(20,ionmax)
        k = max(ionmax-19,1)
        write(6,*) j,' largest mass fractions'
        do i=ionmax,k,-1
         if (xin(izwork1(i)) .gt. 1.0e-12) &
            write(6,02) ionam(izwork1(i)),xin(izwork1(i))
        end do

! nonconservation, abar, zbar of the mixture
        sum = 0.0d0
         do i=1,ionmax
          sum = sum + xin(i)
         enddo
        write(6,02) '1-sum',1.0d0 - sum
        write(6,02) 'abar',abar,'zbar',zbar,'ye',zbar/abar
        write(6,*)

!        read(5,*)



! there is probably a better place for this
! if requested, adjust the number of equations being solved
      if (pure_network .eq. 1) then
       neqs  = ionmax
       btemp = tin
       bden  = din
      end if

      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------
!
! this routine contains auxillary network routine

! routines for a tree construction to mark nonzero matrix locations
! routine screen6 computes screening factors
! routine screen5 computes screening factors
! routine snupp computes neutrino loss rates for the pp chain
! routine snucno computes neutrino loss rates for the cno cycles
! routine sneut5 computes neutrino loss rates
! routine ifermi12 does an inverse fermi integral of order 1/2
! routine zfermim12 does an inverse fermi integral of order -1/2

! routine ecapnuc02 computes electron capture rates
! routine ecapnuc computes electron capture rates
! routine mazurek computes ni56 electron capture rates
! routine time_scales computes various timescales
! routine ener_gener_rate computes the instantaneous energy generation rate








!---------------------------------------------------------------------
      subroutine net_output(kount,x,y,derivs)
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
      include 'burn_common.dek'
      include 'network.dek'
      include 'cjdet.dek'

! writes the output

! declare the pass
      external         derivs
      integer          kount
      double precision x,y(*)


! local variables
      character*8      atim
      character*9      adat
      character*80     string
      integer          k,kk,j,lop,ilop,jrem,kb,ke,nn,lenstr
      double precision sum,xcons,ycons,yex,ydum(abignet), &
                       dydt_dum(nzmax*abignet),xdum(abignet), &
                       abar,zbar,wbar,ye,xcess,zero,tdum,ddum,pdum, &
                       ener,denerdt,zc12,xc12,ff, &
                       chem_pot(nzmax*abignet),chem_sum, &
                       ydum_sav(nzmax*abignet),posx,velx
      parameter        (zero = 0.0d0)


! for nse
      integer          igues
      double precision xmun,xmup,t9,tau_nse,tau_qse,taud


! popular format statements
01    format(1x,'*',t13,a,t33,a,t47,a,t61,a,t75,a,t89,a, &
                    t103,a,t117,a,t131,a,t145,a,t159,a)
03    format(a30,i4.4,a2,i8,a)
04    format(1x,i6,1pe20.12,1p15e14.6)
05    format(1x,i6,1pe20.12,1p12e14.6)
07    format(1x,'* ',a,5(a,1pe11.3))



!      write(6,*) kount,neqs,nzone

! initialize the files with their headers
      if (kount .eq. 1) then



! for every spatial zone
       do k=1,max(1,nzone)
        kk = neqs*(k-1)

! logical unit 22 records the energetics
        write(string,03) hfile,0,'_z',k,'.dat'
        call sqeeze(string)
        call today(adat,atim)
        open (unit=22, file=string, status='unknown')


! logical unit 23 records the thermodynamics
        write(string,03) hfile,1,'_z',k,'.dat'
        call sqeeze(string)
        open (unit=23, file=string, status='unknown')


         write(22,01) adat,atim
         write(23,01) adat,atim

        if (one_step) then
         write(22,07) 'one_step:','  btemp=',btemp,' bden=',bden
         write(23,07) 'one_step:','  btemp=',btemp,' bden=',bden

        else if (hydrostatic) then
         write(22,07) 'hydrostatic:','  btemp=',btemp,' bden=',bden
         write(23,07) 'hydrostatic:','  btemp=',btemp,' bden=',bden

        else if (expansion) then
         write(22,07) 'expansion:','  temp0=',temp0,' den0=',den0, &
                     ' temp_stop=',temp_stop
         write(23,07) 'expansion:','  temp0=',temp0,' den0=',den0, &
                     ' temp_stop=',temp_stop

        else if (self_heat_const_den) then
         write(22,07) 'self_heat:','   temp0=',y(itemp+kk), &
                                  '   den0=',y(iden+kk)
         write(23,07) 'self_heat:','   temp0=',y(itemp+kk), &
                                  '   den0=',y(iden+kk)

        else if (self_heat_const_pres) then
         write(22,07) 'self_heat:','   temp0=',y(itemp+kk), &
                                  '   den0=',y(iden+kk)
         write(23,07) 'self_heat:','   temp0=',y(itemp+kk), &
                                  '   den0=',y(iden+kk)

        else if (bbang) then
         write(22,07) 'big bang:','   eta=',eta1, &
                                  '   Nnu=',xnnu, &
                                  '   H0=',hubble, &
                                  '   Tcmb=',cmbtemp
         write(23,07) 'big_bang:','   eta=',eta1, &
                                  '   Nnu=',xnnu, &
                                  '   H0=',hubble, &
                                  '   Tcmb=',cmbtemp

        else if (detonation) then
         write(22,07) 'detonation:',' temp0=',temp_up, &
                                  '   den0=',den_up, &
                                  '   pres0=',pres_up, &
                                  '   mach=',mach_sh
         write(22,07) '           ',' temp_sh=',temp_sh, &
                                  '   den_sh=',den_sh, &
                                  '   pres_sh=',pres_sh, &
                                  '   vel_sh=',vel_sh, &
                                  '   cs_sh=',cs_sh
         write(22,07) '           ',' temp_cj=',temp_cj, &
                                  '   den_cj=',den_cj, &
                                  '   pres_cj=',pres_cj, &
                                  '   vel_cj=',vel_cj, &
                                  '   cs_cj=',cs_cj

         write(23,07) 'detonation:','  temp0=',temp_up, &
                                  '   den0=',den_up, &
                                  '   pres0=',pres_up, &
                                  '   mach=',mach_sh
         write(23,07) '           ',' temp_sh=',temp_sh, &
                                  '   den_sh=',den_sh, &
                                  '   pres_sh=',pres_sh, &
                                  '   vel_sh=',vel_sh, &
                                  '   cs_sh=',cs_sh
         write(23,07) '           ',' temp_cj=',temp_cj, &
                                  '   den_cj=',den_cj, &
                                  '   pres_cj=',pres_cj, &
                                  '   vel_cj=',vel_cj, &
                                  '   cs_cj=',cs_cj

        else if (trho_hist) then
         call update2(zero,tdum,ddum)
         write(22,07) 'trho_hist:','  mass interior =',mint, &
                                   '  shell mass =',mshell
         write(23,07) 'trho_hist:',' mass interior =',mint, &
                                   '  shell mass =',mshell

        else if (pt_hist) then
         call update3(zero,tdum,pdum)
         write(22,07) 'pt_hist:  ','  mass interior =',mint, &
                                   '  shell mass =',mshell
         write(23,07) 'pt_hist:  ',' mass interior =',mint, &
                                   '  shell mass =',mshell
        end if


        write(22,01) 'time','temp','den','ener','sdot','sneut', &
                     's-snu','ye','1-sum'

        write(23,01) 'time','pos','vel','temp','den','pres','ener', &
                     'entr','cs'


! close up the files
        close(unit=22)
        close(unit=23)


! end of spatial loop
       enddo



! if we are doing an nse analysis, we'll write out another file

       if (nse_analysis .eq. 1) then

! for every spatial zone
       do k=1,max(1,nzone)
        kk = neqs*(k-1)

! logical unit 25 records the nse analysis
        write(string,03) hfile,0,'_z',k,'_nse.dat'
        call sqeeze(string)
        call today(adat,atim)
        open (unit=25, file=string, status='unknown')

        write(25,01) adat,atim

        if (one_step) then
         write(25,07) 'one_step:','  btemp=',btemp,' bden=',bden

        else if (hydrostatic) then
         write(25,07) 'hydrostatic:','  btemp=',btemp,' bden=',bden

        else if (expansion) then
         write(25,07) 'expansion:','  temp0=',temp0,' den0=',den0, &
                     ' temp_stop=',temp_stop

        else if (self_heat_const_den) then
         write(25,07) 'self_heat:','   temp0=',y(itemp+kk), &
                                  '   den0=',y(iden+kk)

        else if (self_heat_const_pres) then
         write(25,07) 'self_heat:','   temp0=',y(itemp+kk), &
                                  '   den0=',y(iden+kk)

        else if (bbang) then
         write(25,07) 'big bang:','   eta=',eta1, &
                                  '   Nnu=',xnnu, &
                                  '   H0=',hubble, &
                                  '   Tcmb=',cmbtemp

        else if (detonation) then
         write(25,07) 'detonation:','  temp0=',temp_up, &
                                  '   den0=',den_up, &
                                  '   pres0=',pres_up, &
                                  '   mach=',mach_sh
         write(25,07) '           ',' temp_sh=',temp_sh, &
                                  '   den_sh=',den_sh, &
                                  '   pres_sh=',pres_sh, &
                                  '   vel_sh=',vel_sh, &
                                  '   cs_sh=',cs_sh
         write(25,07) '           ',' temp_cj=',temp_cj, &
                                  '   den_cj=',den_cj, &
                                  '   pres_cj=',pres_cj, &
                                  '   vel_cj=',vel_cj, &
                                  '   cs_cj=',cs_cj

        else if (trho_hist) then
         call update2(zero,tdum,ddum)
         write(25,07) 'trho_hist:','  mass interior =',mint, &
                                   '  shell mass =',mshell

        else if (pt_hist) then
         call update3(zero,tdum,pdum)
         write(25,07) 'pt_hist:  ','  mass interior =',mint, &
                                   '  shell mass =',mshell
        end if


        write(25,01) 'time','temp','den','ye','tqse','tnse','delta', &
                     '1-sum'


! close up the files
        close(unit=25)

! end of spatial loop
       enddo

       end if



! done writing thermodynamic headers



! for every spatial zone
       do k=1,max(1,nzone)
        kk = neqs*(k-1)


! write out the isotopic mass fractions in blocks of 8
! lop is how many groups of 8 exist; jrem is the remainder
        lop  = ionmax/8
        jrem  = ionmax - 8*lop
        do ilop = 1,lop+1
         kb = 1 + 8*(ilop-1)
         ke = 8 + 8*(ilop-1)
         if (ilop .eq. lop+1  .and. jrem .eq. 0) goto 50
         if (ilop .eq. lop+1) ke = ionmax


! logical unit 34 records the abundance evolution
! open the output file
         write(string,03) hfile,ilop+1,'_z',k,'.dat'
         call sqeeze(string)
         open (unit=34, file=string, status='unknown')

         write(34,01) adat,atim


        if (one_step) then
         write(34,07) 'one_step:','  btemp=',btemp,' bden=',bden

        else if (hydrostatic) then
         write(34,07) 'hydrostatic:','  btemp=',btemp,' bden=',bden

        else if (expansion) then
         write(34,07) 'expansion:','  temp0=',temp0,' den0=',den0, &
                                  ' temp_stop=',temp_stop

        else if (self_heat_const_den) then
         write(34,07) 'self_heat:','   temp0=',y(itemp+kk), &
                                  '   den0=',y(iden+kk)

        else if (bbang) then
         write(34,07) 'big bang:','   eta=',eta1, &
                                  '   Nnu=',xnnu, &
                                  '   H0=',hubble, &
                                  '   Tcmb=',cmbtemp

        else if (detonation) then
         write(34,07) 'detonation:','  temp0=',temp_up, &
                                  '   den0=',den_up, &
                                  '   pres0=',pres_up, &
                                  '   mach=',mach_sh
         write(34,07) '           ',' temp_sh=',temp_sh, &
                                  '   den_sh=',den_sh, &
                                  '   pres_sh=',pres_sh, &
                                  '   vel_sh=',vel_sh, &
                                  '   cs_sh=',cs_sh
         write(34,07) '           ',' temp_cj=',temp_cj, &
                                  '   den_cj=',den_cj, &
                                  '   pres_cj=',pres_cj, &
                                  '   vel_cj=',vel_cj, &
                                  '   cs_cj=',cs_cj


        else if (trho_hist) then
         call update2(zero,tdum,ddum)
         write(34,07) 'trho_hist:','  mass interior =',mint, &
                                   '  shell mass =',mshell

        else if (pt_hist) then
         call update3(zero,tdum,ddum)
         write(34,07) 'pt_hist:  ','  mass interior =',mint, &
                                   '  shell mass =',mshell

        end if

        write(34,01) 'time',(ionam(nn), nn=kb,ke)

        close(unit=34)
 50     continue
       enddo

! end of the spatial loop
      enddo



! if we are doing an nse analysis, we'll write out another
! set of abundance file

       if (nse_analysis .eq. 1) then


! for every spatial zone
       do k=1,max(1,nzone)
        kk = neqs*(k-1)

! write out the isotopic mass fractions in blocks of 8
! lop is how many groups of 8 exist; jrem is the remainder
        lop  = ionmax/8
        jrem  = ionmax - 8*lop
        do ilop = 1,lop+1
         kb = 1 + 8*(ilop-1)
         ke = 8 + 8*(ilop-1)
         if (ilop .eq. lop+1  .and. jrem .eq. 0) goto 60
         if (ilop .eq. lop+1) ke = ionmax


! logical unit 35 records the abundance evolution
! open the output file
         write(string,03) hfile,ilop+1,'_z',k,'_nse.dat'
         call sqeeze(string)
         open (unit=35, file=string, status='unknown')

         write(35,01) adat,atim


        if (one_step) then
         write(35,07) 'one_step:','  btemp=',btemp,' bden=',bden

        else if (hydrostatic) then
         write(35,07) 'hydrostatic:','  btemp=',btemp,' bden=',bden

        else if (expansion) then
         write(35,07) 'expansion:','  temp0=',temp0,' den0=',den0, &
                                  ' temp_stop=',temp_stop

        else if (self_heat_const_den) then
         write(35,07) 'self_heat:','   temp0=',y(itemp+kk), &
                                  '   den0=',y(iden+kk)

        else if (bbang) then
         write(35,07) 'big bang:','   eta=',eta1, &
                                  '   Nnu=',xnnu, &
                                  '   H0=',hubble, &
                                  '   Tcmb=',cmbtemp

        else if (detonation) then
         write(35,07) 'detonation:','  temp0=',temp_up, &
                                  '   den0=',den_up, &
                                  '   pres0=',pres_up, &
                                  '   mach=',mach_sh
         write(35,07) '           ',' temp_sh=',temp_sh, &
                                  '   den_sh=',den_sh, &
                                  '   pres_sh=',pres_sh, &
                                  '   vel_sh=',vel_sh, &
                                  '   cs_sh=',cs_sh
         write(35,07) '           ',' temp_cj=',temp_cj, &
                                  '   den_cj=',den_cj, &
                                  '   pres_cj=',pres_cj, &
                                  '   vel_cj=',vel_cj, &
                                  '   cs_cj=',cs_cj


        else if (trho_hist) then
         call update2(zero,tdum,ddum)
         write(35,07) 'trho_hist:','  mass interior =',mint, &
                                   '  shell mass =',mshell

        else if (pt_hist) then
         call update3(zero,tdum,ddum)
         write(35,07) 'pt_hist:  ','  mass interior =',mint, &
                                   '  shell mass =',mshell

        end if

        write(35,01) 'time',(ionam(nn), nn=kb,ke)

        close(unit=35)
 60     continue
       enddo

! end of the spatial loop and nse analyis test if
      enddo
      end if

!       write(6,*) 'wrote mass fraction headers'

! end of the file initialization
      end if

!      write(6,*) 'done with initialization'







! normal execution starts here

! for any time point

! for every spatial zone
      do k=1,max(1,nzone)
       kk = neqs*(k-1)

! open the files in append mode (f77) or position mode (f90)

! energetics file
       write(string,03) hfile,0,'_z',k,'.dat'
       call sqeeze(string)
!       open (unit=22, file=string, status='old', access='append')
       open (unit=22, file=string, status='old', position='append')


! thermodynamics file
       write(string,03) hfile,1,'_z',k,'.dat'
       call sqeeze(string)
!       open (unit=23, file=string, status='old', access='append')
       open (unit=23, file=string, status='old', position='append')


! form the mass fractions
       do j=1,ionmax
        xdum(j) = min(1.0d0,max(y(j+kk)*aion(j),1.0d-30))
       enddo


! mass conservation
       sum = 0.0d0
       do j=1,ionmax
        sum = sum + xdum(j)
       enddo
       sum = 1.0d0 - sum
       xcons = sum


! y sum
!       sum = 0.0d0
!       do j=1,ionmax
!        if (zion(j) .gt. 2.0) then
!         sum = sum + max(y(j+kk),1.0d-30)
!        endif
!       enddo
!       ycons = sum


! get ye using normalized mass fractions
       sum = 0.0d0
       do j=1,ionmax
        sum = sum + xdum(j)
       enddo
       sum = 1.0d0/sum
       do j=1,ionmax
        xdum(j) = min(1.0d0,max(sum*xdum(j),1.0d-30))
       enddo


! get abar, zbar and a few other composition variables
       call azbar(xdum(ionbeg),aion(ionbeg),zion(ionbeg),wion(ionbeg),ionmax, &
                  ydum(ionbeg),abar,zbar,wbar,yex,xcess)




! get the right hand sides, exact energy generation rate and so on
       if (nse_on .eq. 0) then
        call derivs(x,y,dydt_dum)
        if (pure_network .eq. 0) then
         ener = y(iener + kk)
         denerdt = dydt_dum(iener + kk)
        else
         ener = 0.0d0
         denerdt = 0.0d0
        end if
       else
        sdot    = 0.0d0
        sneut   = 0.0d0
        ener    = 0.0d0
        denerdt = 0.0d0
       end if


! call an eos
       if (pure_network .eq. 0) then
        temp_row(1) = y(itemp+kk)
        den_row(1)  = y(iden+kk)
       else
        temp_row(1) = btemp
        den_row(1)  = bden
       end if
       if (trho_hist) call update2(x,temp_row(1),den_row(1))
       abar_row(1) = abar
       zbar_row(1) = zbar
       jlo_eos = 1
       jhi_eos = 1

      if (pt_hist) then
       call update3(x,temp_row(1),bpres)
       den_row(1)  = bpres * abar/(avo * kerg * temp_row(1))
       call invert_helm_pt
      else
       call helmeos
!       call eosfxt
      end if


! figure some time scales
       call time_scales(temp_row(1),den_row(1),taud,tau_nse,tau_qse)


! compute the chemical potentials
!       do j=1,ionmax
!        chem_pot(j) = abar*((zion(j) - zbar)*deionz_row(1) &
!                          + (aion(j) - abar)*deiona_row(1))
!       end do
!       sum = 0.0d0
!       do j=1,ionmax
!        sum = sum + chem_pot(j) * dydt_dum(j)
!       end do
!       chem_sum = sum



! and write what we found


! total c12+c12 rate, mass fraction of c12, function
!       zc12 = ratdum(ir1212n) + ratdum(ir1212p) + ratdum(ir1212a)
!       xc12 = y(ic12)*aion(ic12)
!       ff   = sdot/(y(ic12)**2 * zc12) * 2.0d0/3.0d0

       write(22,05) kount,x,temp_row(1),den_row(1), &
                    ener,sdot,sneut,denerdt,yex,xcons


!                    chem_sum,chem_sum/denerdt
!     2              xc12,zc12/den_row(1),ff


       if (iposx .eq. 0) then 
        posx = 0.0d0
       else
        posx = y(iposx+kk)
       end if
       if (ivelx .eq. 0) then 
        velx = 0.0d0
       else
        velx = y(ivelx+kk)
       end if

       write(23,05) kount,x,posx,velx, &
                   temp_row(1),den_row(1),ptot_row(1), &
                   ener,stot_row(1),cs_row(1)



! close up the files
       close(unit=22)
       close(unit=23)

! end of spatial loop
      end do


!      write(6,*) 'done with thermo file'




! for every spatial zone
      do k=1,max(1,nzone)
       kk = neqs*(k-1)

! write out the isotopic mass fractions in blocks of 8
! lop is how many groups of 8 exist; jrem is the remainder
       lop  = ionmax/8
       jrem  = ionmax - 8*lop
       do ilop = 1,lop+1
        kb = 1 + 8*(ilop-1)
        ke = 8 + 8*(ilop-1)
        if (ilop .eq. lop+1  .and. jrem .eq. 0) goto 70
        if (ilop .eq. lop+1) ke = ionmax

! open the output file in append mode (f77) or position mode (f90)
! abundance evolution file
        write(string,03) hfile,ilop+1,'_z',k,'.dat'
        call sqeeze(string)
!        open (unit=34, file=string, status='old', access='append')
        open (unit=34, file=string, status='old', position='append')

        write(34,04) kount,x,(y(nn+kk)*aion(nn), nn=kb,ke)
!        write(34,04) kount,x,(y(nn+kk), nn=kb,ke)

        close(unit=34)
70      continue
       enddo


! end of spatial zone loop
      enddo

!      write(6,*) 'done with mass fractions file'




! start of nse analysis

      if (nse_analysis .eq. 1) then

! for every spatial zone
       do k=1,max(1,nzone)
        kk = neqs*(k-1)


! open the files in append mode (f77) or position mode (f90)

! nse analysis file
       write(string,03) hfile,0,'_z',k,'_nse.dat'
       call sqeeze(string)
!       open (unit=25, file=string, status='old', access='append')
       open (unit=25, file=string, status='old', position='append')


! form the mass fractions
       do j=1,ionmax
        xdum(j) = min(1.0d0,max(y(j+kk)*aion(j),1.0d-30))
       enddo



! normalized mass fractions
       sum = 0.0d0
       do j=1,ionmax
        sum = sum + xdum(j)
       enddo
       xcons = 1.0d0 - sum
       sum = 1.0d0/sum
       do j=1,ionmax
        xdum(j) = min(1.0d0,max(sum*xdum(j),1.0d-30))
       enddo


! get abar, zbar and a few other composition variables
       call azbar(xdum(ionbeg),aion(ionbeg),zion(ionbeg),wion(ionbeg),ionmax, &
                  ydum(ionbeg),abar,zbar,wbar,yex,xcess)



! set the temperature and density
       if (pure_network .eq. 0) then
        temp_row(1) = y(itemp+kk)
        den_row(1)  = y(iden+kk)
       else
        temp_row(1) = btemp
        den_row(1)  = bden
       end if
       if (trho_hist) call update2(x,temp_row(1),den_row(1))
       if (pt_hist) then
        call update3(x,temp_row(1),bpres)
        den_row(1)  = bpres * abar/(avo * kerg * temp_row(1))
        call invert_helm_pt
       end if


! with the temperature, density, and ye
! compute the nse state if the temperature is high enough

       if (temp_row(1) .gt. 2.0e9) then
        igues = 1
        call nse(temp_row(1),den_row(1),yex,igues,1,1,xsum,xmun,xmup,0)
       else
        do j=1,ionmax
         xsum(j) = 1.0e20
        enddo
       end if

! figure delta on the top 20 nse mass fractions
       call indexx(ionmax,xsum(ionbeg),izwork1(ionbeg))
       sum = 0.0d0
       kb  = 0
       do j = ionmax, max(1,ionmax-19), -1
        if (xsum(izwork1(j)) .ge. 1.0e-6) then
         kb = kb + 1
         tdum = (xsum(izwork1(j)) - xdum(izwork1(j)))/xsum(izwork1(j))
!         tdum = (xsum(izwork1(j)) - xdum(izwork1(j)))**2
         sum  = sum + tdum
        end if
       enddo
       sum = sum/float(kb)
!       sum = sqrt(sum/kb)


! figure the time scales
       call time_scales(temp_row(1),den_row(1),taud,tau_nse,tau_qse)



! write out what we got
       write(25,05) kount,x,temp_row(1),den_row(1),yex, &
                   tau_qse,tau_nse,sum,xcons

! close up the files
       close(unit=25)


! write out the isotopic mass fractions in blocks of 8
! lop is how many groups of 8 exist; jrem is the remainder
       lop  = ionmax/8
       jrem  = ionmax - 8*lop
       do ilop = 1,lop+1
        kb = 1 + 8*(ilop-1)
        ke = 8 + 8*(ilop-1)
        if (ilop .eq. lop+1  .and. jrem .eq. 0) goto 80
        if (ilop .eq. lop+1) ke = ionmax


! open the output file in append mode (f77) or position mode (f90)
! abundance evolution file

        write(string,03) hfile,ilop+1,'_z',k,'_nse.dat'
        call sqeeze(string)
!        open (unit=35, file=string, status='old', access='append')
        open (unit=35, file=string, status='old', position='append')

        write(35,04) kount,x,(xsum(nn+kk), nn=kb,ke)

        close(unit=35)
80      continue
       enddo


! end of spatial zone loop
      enddo

! end of the nse analysis if
      end if



      return
      end
!---------------------------------------------------------------------











!---------------------------------------------------------------------
! reaction rate library

! torch rates
! li7(t,n)   a(an,g)    be9(p,d)    be9(p,n)    b10(a,n)   b11(a,n)
! n14(p,a)   c11(p,g)   c12(a,n)    c13(a,n)    c13(p,n)   c14(a,g)
! c14(p,n)   c14(p,g)   o16(p,a)    n14(p,n)    n14(a,n)   n15(p,n)
! n15(a,n)   n15(a,g)   o14(a,g)    o17(a,g)    o17(a,n)   o18(a,g)
! o18(a,n)   ne20(p,a)  f18(p,g)    f19(p,g)    f19(p,n)   f19(a,p)
! na22(n,a)  ne20(p,g)  na23(p,a)   ne20(n,g)   ne21(p,g)  ne21(a,g)
! ne22(p,g)  ne22(a,g)  na22(n,p)   ne22(a,n)   na21(p,g)  mg24(p,a)
! ne21(a,n)  na22(p,g)  na23(p,g)   na23(p,n)   mg24(p,g)  al27(p,a)
! mg25(p,g)  mg25(a,p)  mg25(a,g)   mg25(a,n)   mg26(p,g)  mg26(a,g)
! mg26(a,n)  al25(p,g)  al26(p,g)   al27(a,n)   si27(p,g)  si28(p,g)
! si29(p,g)  si30(p,g)

! bigbang rates:
! n(e-nu)p   p(e-,nu)n  d(p,n)      d(n,g)      d(d,p)     d(d,n)
! t(p,n)     d(d,g)     t(p,g)      t(d,n)      t(t,2n)    he3(d,p)
! he3(t,d)   he3(t,np)  he4(np,g)   he4(d,g)    he4(t,n)   li6(p,he3)
! li6(n,g)   li7(d,n)   lit(t,2n)   li7(he3,np) li6(p,g)   li7(p,n)
! be7(d,p)   be7(t,np)  be7(3he,2p) li6(a,g)    li7(a,n)   be9(p,g)
! b10(p,a)   li7(a,g)   b11(p,a)    be7(a,g)    b11(p,n)   b8(a,p)
! b10(p,g)   c11(n,a)   be9(a,n)    b11(p,g)    b11(a,p)

! pp123 rates:
! p(p,e+nu)  p(n,g)     d(p,g)      he3(n,g)    he3+he3    he3(a,g)
! be7(e-,nu) be7(p,g)   li7(p,g)    li7(p,a)    b8(e+,nu)

! cno rates:
! c12(p,g)   n13(e-nu)  c13(p,g)    n14(p,g)    o15(e-nu)  n14(a,g)
! n15(p,g)   n15(p,a)   o16(p,g)    o17(p,a)    o17(p,g)   o18(p,a)
! o18(p,g)   f17(e-nu)  f18(e-nu)   f19(p,a)

! hot cno rates
! n13(p,g)   o14(e-nu)  o14(a,p)    o15(a,g)    f17(p,g)   ne18(e-nu)
! f18(p,a)   ne18(a,p)  ne19(p,g)   ne19(e-nu)  si26(a,p)

! alfa chain rates:
! a(aa,g)    c12(a,g)   c12+c12     c12+o16     o16+o16    o16(a,g)
! ne20(a,g)  ne20(a,g)  mg24(a,g)   mg24(a,p)   al27(p,g)  si28(a,g)
! si28(a,p)  p31(p,g)   s32(a,g)    s32(a,p)    cl35(p,g)  ar36(a,g)
! ar36(a,p)  k39(p,g)   ca40(a,g)   ca40(a,p)   sc43(p,g)  ti44(a,g)
! ti44(a,p)  v47(p,g)   cr48(a,g)   cr(a,p)     mn51(p,g)  fe52(a,g)
! fe52(a,p)  co55(p,g)

! photodisintegration rates:
! fe52(n,g) fe53(n,g)  fe54(p,g)












!---------------------------------------------------------------------







!---------------------------------------------------------------------
      subroutine net_initialize
      include 'implno.dek'
      include 'network.dek'

! initializes quantities

! local variables
      integer   i


! general options
      screen_on      = 1
      use_tables     = 1
      weak_on        = 1
      ffn_on         = 0
      pure_network   = 0
      nse_analysis   = 0
      allow_nse_evol = 0


! printing information
      iprint_files  = 1
      iprint_screen = 1


! inititailize the burn type logicals
      one_step             = .false.
      hydrostatic          = .false.
      expansion            = .false.
      self_heat_const_den  = .false.
      self_heat_const_pres = .false.
      pt_hist              = .false.
      bbang                = .false.
      detonation           = .false.
      trho_hist            = .false.


! adiabatic expanion off
      psi       = 0.0d0
      temp_stop = 1.0d30


! mass fractions above sthreshold are written to the summary file
      sthreshold = 1.0d30

      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------
      subroutine net_summary(tstep,tin,din,ein,tout,dout,eout,conserv, &
                             nbad,nok,xout)
      include 'implno.dek'
      include 'timers.dek'
      include 'vector_eos.dek'
      include 'burn_common.dek'
      include 'network.dek'


! writes out a summary of the network run

! declare the pass
      integer          nbad,nok
      double precision tstep,tin,din,ein,tout,dout,eout,conserv, &
                       xout(*)

! local variables
      character*80     summary
      integer          i,j,k,lenstr,ioff
      double precision abar,zbar,wbar,ye,xcess


! popular format statements
 01   format(a,'summary.dat')
 02   format(1x,a,'=',1pe10.3,' ',a,'=',1pe10.3,' ', &
                a,'=',1pe10.3,' ',a,'=',1pe10.3,' ', &
                a,'=',1pe10.3)
 03   format(1x,a,1pe20.12)
 04   format(1x,a,':',/, &
             1x,3(a,1pe20.12),/, &
             1x,3(a,1pe20.12),/, &
             1x,2(a,1pe11.3),2(a,i5))
 08   format(1x,a,1pe10.3,a)
 09   format(1x,a,i2,a)



! construct the file name and open it
       write(summary,01) hfile(1:lenstr(hfile,80))
       call sqeeze(summary)
       open(unit=41,file=summary,status='unknown')


       write(6,*) ' '
       write(6,04) netname, &
                   ' tin =',tin,' din =',din,' ein =',ein, &
                   ' tout=',tout,' dout=',dout,' eout=',eout, &
                   ' dener=',(eout - ein),' sum =',conserv, &
                   ' nbad=',nbad,' nok=',nok
       write(6,*) ' '

       write(41,*) ' '
       write(41,04) netname, &
                   ' tin =',tin,' din =',din,' ein =',ein, &
                   ' tout=',tout,' dout=',dout,' eout=',eout, &
                   ' enuc=',(eout - ein)/tstep,' sum =',conserv, &
                   ' nbad=',nbad,' nok=',nok
       write(41,*) ' '



! write out the biggest mass fractions
       call indexx(ionmax,xout(ionbeg),izwork1(ionbeg))
       ioff = ionbeg - 1

       if (sthreshold .le. 1  .and. sthreshold .gt. 0.0) then
        do i=ionmax,1,-1
         if (xout(izwork1(i)+ioff) .lt. sthreshold) then
          k = i + 1
          write(6,08)  'mass fractions larger than ',sthreshold
          write(41,08) 'mass fractions larger than ',sthreshold
          goto 20
         end if
        end do
       else
        j = min(20,ionmax)
        k = max(ionmax-19,ionbeg)
        write(6,09)  'top ',j,' mass fractions:'
        write(41,09) 'top ',j,' mass fractions:'
       end if

 20   continue


       write(6,02) (ionam(izwork1(i)+ioff),xout(izwork1(i)+ioff), i=ionmax,k,-1)
       if (iprot .ne. 0 .and. ineut .ne. 0) then
        write(6,02) ionam(iprot),xout(iprot), &
                    ionam(ineut),xout(ineut), &
                    ionam(ihe4),xout(ihe4)
       end if
       write(6,*) ' '

       write(41,02) (ionam(izwork1(i)+ioff),xout(izwork1(i)+ioff), i=ionmax,k,-1)
       if (iprot .ne. 0 .and. ineut .ne. 0) then
        write(41,02) ionam(iprot),xout(iprot), &
                     ionam(ineut),xout(ineut), &
                     ionam(ihe4),xout(ihe4)
       end if
       write(41,*) ' '



! end the clock
      call zsecond(timtot)
      timtot = timtot - timzer
      call timlap(timtot,hours,minuts,secs,msecs)
      write(6,100) hours,minuts,secs,msecs
      write(41,100) hours,minuts,secs,msecs
 100  format(1x,'cpu time : ',i2.2,' hrs  ',i2.2,' min  ', &
                              i2.2,' sec  ',i6,' usec',/,/)


! close up shop
      close(unit=41)
      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------
! this file contains routines that sort, search and select parts of arrays:
!
! index and rank makers:
! routine indexx constructs a sort index for a real array



      subroutine indexx(n,arr,indx)
      include 'implno.dek'
!
! indexes an array arr(1:n). that is it outputs the array indx(1:n) such
! that arr(indx(j)) is in ascending order for j=1...n. the input quantities
! are not changed.
!
! declare
      integer          n,indx(n),m,nstack
      parameter        (m=7, nstack = 50)
      integer          i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
      double precision arr(n),a
!
! initialize
      do 11 j=1,n
       indx(j) = j
11    continue
      jstack = 0
      l      = 1
      ir     = n
!
! insertion sort when subbarray small enough
1     if (ir - l .lt. m) then
       do 13 j=l+1,ir
        indxt = indx(j)
        a     = arr(indxt)
        do 12 i=j-1,l,-1
         if (arr(indx(i)) .le. a) go to 2
         indx(i+1) = indx(i)
12      continue
        i = l - 1
2       indx(i+1) = indxt
13     continue
!
! pop stack and begin a new round of partitioning
       if (jstack .eq. 0) return
       ir     = istack(jstack)
       l      = istack(jstack-1)
       jstack = jstack - 2
!
! choose median of left, center and right elements as partitioning element
! also rearrange so that a(l+1) < a(l) < a(ir)
      else
       k         = (l + ir)/2
       itemp     = indx(k)
       indx(k)   = indx(l+1)
       indx(l+1) = itemp

       if (arr(indx(l)) .gt. arr(indx(ir))) then
        itemp    = indx(l)
        indx(l)  = indx(ir)
        indx(ir) = itemp
       end if


       if(arr(indx(l+1)).gt.arr(indx(ir)))then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
       endif
       if(arr(indx(l)).gt.arr(indx(l+1)))then
        itemp=indx(l)
        indx(l)=indx(l+1)
        indx(l+1)=itemp
       endif

!
! initialize pointers for partitioning
       i     = l + 1
       j     = ir
       indxt = indx(l+1)
       a     = arr(indxt)
3      continue
       i = i + 1
       if (arr(indx(i)) .lt. a) go to 3
4      continue
       j = j - 1
       if (arr(indx(j)) .gt. a) go to 4
       if (j .lt. i) go to 5
       itemp   = indx(i)
       indx(i) = indx(j)
       indx(j) = itemp
       go to 3
!
5      indx(l+1) = indx(j)
       indx(j)   = indxt
       jstack    = jstack + 2
!
! push pointers to larger subarray on stack
       if (jstack .gt. nstack) stop 'jstack > nstack in routine indexx'
       if (ir - i + 1  .ge.  j - l) then
        istack(jstack)   = ir
        istack(jstack-1) = i
        ir               = j - 1
       else
        istack(jstack)   = j-1
        istack(jstack-1) = l
        l                = i
       end if
      end if
      go to 1
      end
!---------------------------------------------------------------------








!---------------------------------------------------------------------





