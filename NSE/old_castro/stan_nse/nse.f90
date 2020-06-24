module nse_module

  use bl_types
  use bl_constants_module
  use bl_error_module


contains

   subroutine nsetable_init()
   use bl_error_module
   use parallel
   implicit none

   double precision  ttlog(46221),ddlog(46221),yetab(46221), &
                     abartab(46221),ebtab(46221),wratetab(46221),helium(46221), &
                     sica(46221), irongroup(46221)
   common /table_nse/  ttlog,ddlog,yetab,abartab,ebtab,wratetab,helium,sica,irongroup
   save /table_nse/

   integer          :: j,ir,it,iye
   integer          :: ntemp,nden,nye

   logical          first1
   data             first1/.false./

   if ( first1 ) return

   first1 = .true.

!.. open the table
   open(unit=3,file='nse_table.dat',status='old',err=1010)
   GO TO 1011
 1010 continue
   call bl_error('EOSINIT: Failed to open nse_table.dat')
 1011 continue

!  set table parameters
   ntemp=71
   nden=31
   nye=21

!  read in table
   do ir = 1,nden
   do it = 1,ntemp
      do iye = 1,nye
         j =(ir-1)*ntemp*nye + (it-1)*nye + iye
         read (3,5)  ttlog(j),ddlog(j),yetab(j),helium(j),sica(j),irongroup(j),abartab(j),&
                    ebtab(j),wratetab(j)
 5 format(2f12.5,1pe12.5,6e12.5)
      end do
   end do
   end do

      
   close(3)

   if ( parallel_ioprocessor() ) then
        write(6,*)
        write(6,*) 'finished reading nse table'
        write(6,*)
   end if

   end subroutine nsetable_init


  subroutine interp (t9,rho,ye,abar,dq,dyedt,he,si_ca,fe)
  use bl_error_module
  use parallel
  implicit none

  double precision t9,rho,ye,abar,dq,dyedt,he,si_ca,fe
  
   double precision  ttlog(46221),ddlog(46221),yetab(46221), &
                     abartab(46221),ebtab(46221),wratetab(46221),helium(46221), &
                     sica(46221), irongroup(46221)
   common /table_nse/  ttlog,ddlog,yetab,abartab,ebtab,wratetab,helium,sica,irongroup
   save /table_nse/

  double precision tlog,rholog, yet
  double precision t0,r0,x0,td,rd,xd,omtd,omrd,omxd,t932,t9a
  double precision t9a13,t9a56,r1212,pene1,pene2,x12,x12sq
  integer it1,it2,ir1,ir2,ic1,ic2,it1r1c1,it1r1c2,it1r2c1
  integer it1r2c2,it2r1c1,it2r1c2,it2r2c1,it2r2c2


  tlog = log10(1.0d9*t9)  
  rholog = log10(rho)
  yet = ye
 

  if (tlog.lt.9.0d0) tlog = 9.0d0
  if (tlog.gt.10.4d0) tlog = 10.4d0

  it1 = int((tlog -9.0d0)*50.0d0 - 1.d-6)
  it1 = it1 + 1 
  it2 = it1 + 1

  if (rholog.lt.7.0d0) rholog = 7.0d0
  if (rholog.gt.10.0d0) rholog = 10.0d0

  ir1 = int((rholog -7.0d0)*10.0d0 - 1.d-6)
  ir1 = ir1 + 1
  ir2 = ir1+1

  if (yet.lt.0.40d0) yet = 0.40d0
  if (yet.gt.0.50d0) yet = 0.50d0
      
  ic1 = int((0.50d0-yet)/0.005d0 - 1.0d-6)
  ic1 = ic1 + 1
  ic2 = ic1 + 1


! write (6,*) it1,it2,ir1,ir2,ic1,ic2

! find the eight interpolation points in the 1D arrays

  it1r1c1 = (ir1-1)*71*21 + (it1-1)*21 + ic1 
  it1r1c2 = (ir1-1)*71*21 + (it1-1)*21 + ic2
  it1r2c1 = (ir2-1)*71*21 + (it1-1)*21 + ic1
  it1r2c2 = (ir2-1)*71*21 + (it1-1)*21 + ic2
  it2r1c1 = (ir1-1)*71*21 + (it2-1)*21 + ic1
  it2r1c2 = (ir1-1)*71*21 + (it2-1)*21 + ic2
  it2r2c1 = (ir2-1)*71*21 + (it2-1)*21 + ic1
  it2r2c2 = (ir2-1)*71*21 + (it2-1)*21 + ic2

  t0 = 9.0d0 + real(it1-1)*0.02d0
  r0 = 7.0d0 + real(ir1-1)*0.10d0
  x0 = 0.50d0 - real(ic1-1)*0.005d0

  td = (tlog - t0)/0.02d0
  rd = (rholog - r0)/0.10d0
  xd = (x0-yet)/0.005d0
  xd = max(0.0d0,xd)
  omtd = 1.0d0 - td
  omrd = 1.0d0 - rd
  omxd = 1.0d0 - xd


      abar = abartab(it1r1c1)*omtd*omrd*omxd    &
            + abartab(it1r1c2)*omtd*omrd*xd     &
            + abartab(it1r2c1)*omtd*rd*omxd     &
            + abartab(it1r2c2)*omtd*rd*xd       &
            + abartab(it2r1c1)*td*omrd*omxd     &
            + abartab(it2r1c2)*td*omrd*xd       &
            + abartab(it2r2c1)*td*rd*omxd       &
            + abartab(it2r2c2)*td*rd*xd        

      dq   = ebtab(it1r1c1)*omtd*omrd*omxd      &
            + ebtab(it1r1c2)*omtd*omrd*xd       &
            + ebtab(it1r2c1)*omtd*rd*omxd       &
            + ebtab(it1r2c2)*omtd*rd*xd         &
            + ebtab(it2r1c1)*td*omrd*omxd       &
            + ebtab(it2r1c2)*td*omrd*xd         &
            + ebtab(it2r2c1)*td*rd*omxd         &
            + ebtab(it2r2c2)*td*rd*xd

      dyedt = wratetab(it1r1c1)*omtd*omrd*omxd  &
            + wratetab(it1r1c2)*omtd*omrd*xd    &
            + wratetab(it1r2c1)*omtd*rd*omxd    &
            + wratetab(it1r2c2)*omtd*rd*xd      &
            + wratetab(it2r1c1)*td*omrd*omxd    &
            + wratetab(it2r1c2)*td*omrd*xd      &
            + wratetab(it2r2c1)*td*rd*omxd      &
            + wratetab(it2r2c2)*td*rd*xd

      he    = helium(it1r1c1)*omtd*omrd*omxd  &
            + helium(it1r1c2)*omtd*omrd*xd    &
            + helium(it1r2c1)*omtd*rd*omxd    &
            + helium(it1r2c2)*omtd*rd*xd      &
            + helium(it2r1c1)*td*omrd*omxd    &
            + helium(it2r1c2)*td*omrd*xd      &
            + helium(it2r2c1)*td*rd*omxd      &
            + helium(it2r2c2)*td*rd*xd

      si_ca = sica(it1r1c1)*omtd*omrd*omxd  &
            + sica(it1r1c2)*omtd*omrd*xd    &
            + sica(it1r2c1)*omtd*rd*omxd    &
            + sica(it1r2c2)*omtd*rd*xd      &
            + sica(it2r1c1)*td*omrd*omxd    &
            + sica(it2r1c2)*td*omrd*xd      &
            + sica(it2r2c1)*td*rd*omxd      &
            + sica(it2r2c2)*td*rd*xd


      fe    = irongroup(it1r1c1)*omtd*omrd*omxd  &
            + irongroup(it1r1c2)*omtd*omrd*xd    &
            + irongroup(it1r2c1)*omtd*rd*omxd    &
            + irongroup(it1r2c2)*omtd*rd*xd      &
            + irongroup(it2r1c1)*td*omrd*omxd    &
            + irongroup(it2r1c2)*td*omrd*xd      &
            + irongroup(it2r2c1)*td*rd*omxd      &
            + irongroup(it2r2c2)*td*rd*xd

   end subroutine interp

end module nse_module
