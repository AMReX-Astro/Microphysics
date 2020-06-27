      program usetable

      implicit real*8 (a-h,o-z)

c     begin initialization

      common /table/ ttlog(46221),ddlog(46221),yetab(46221),
     1     abartab(46221),ebtab(46221),wratetab(46221),helium(46221)

      
c     set table parameters

      ntemp=71
      nden=31
      nye=21

c     read in table

      do 30, irho = 1,nden
      do 20, it9 = 1,ntemp
      do 10, iye = 1,nye

      j = (irho-1)*ntemp*nye + (it9-1)*nye + iye  
      read (2,5) ttlog(j),ddlog(j),yetab(j),helium(j),abartab(j),
     1     ebtab(j),wratetab(j)
 5     format(2f12.5,1pe12.5,4e12.5)

 10   continue
 20   continue
 30   continue
 
c     end initialization 

c     enter test points

 9    write (6,*) 'enter 1 for log10 entry; 2 for actual values'
      read (5,*) ichoice
      if (ichoice.eq.1) go to 7
      if (ichoice.eq.2) go to 39
      call exit(1)

 39   write (6,390)
 390  format('Enter t9, rho, ye')
      read (5,*) t9,rho,ye
      go to 8

 7    write (6,391)
 391  format('Enter log 10 t9, log 10 rho, ye')
      read (5,*) xlogt9,xlogrho,ye
      t9 = 10.0d0**(xlogt9-9.0d0)
      rho = 10.0d0**xlogrho

 8    call interp (t9,rho,ye,abar,dq,dyedt)

      write (6,40)
 40   format ('      t9    ','      rho    ','     ye     ',
     1      '   abar     ','    be/a    ','    dyedt  ')
      write (6,41) t9,rho,ye,abar,dq,dyedt
 41   format (1pe12.3,5e12.3)
      go to 9

      call exit(1)
      end


      subroutine interp (t9,rho,ye,abar,dq,dyedt)

      implicit real*8 (a-h,o-z)


      common /table/ ttlog(46221),ddlog(46221),yetab(46221),
     1     abartab(46221),ebtab(46221),wratetab(46221),helium(46221)

 
      tlog = log10(1.0d9*t9)
      rholog = log10(rho)
      yet = ye
 
c      ntemp=71
c      nden=31
c      nye=21

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


c      write (6,*) it1,it2,ir1,ir2,ic1,ic2

c     find the eight interpolation points in the 1D arrays

      it1r1c1 = (ir1-1)*71*21 + (it1-1)*21 + ic1 
      it1r1c2 = (ir1-1)*71*21 + (it1-1)*21 + ic2
      it1r2c1 = (ir2-1)*71*21 + (it1-1)*21 + ic1
      it1r2c2 = (ir2-1)*71*21 + (it1-1)*21 + ic2
      it2r1c1 = (ir1-1)*71*21 + (it2-1)*21 + ic1
      it2r1c2 = (ir1-1)*71*21 + (it2-1)*21 + ic2
      it2r2c1 = (ir2-1)*71*21 + (it2-1)*21 + ic1
      it2r2c2 = (ir2-1)*71*21 + (it2-1)*21 + ic2

c     debugging information
c
c      write (6,50) wratetab(it1r1c1), wratetab(it1r1c2), 
c     1     wratetab(it1r2c1), wratetab(it1r2c2), wratetab(it2r1c1), 
c     2     wratetab(it2r1c2),wratetab(it2r2c1), wratetab(it2r2c2)
c 50    format(1pe12.3,7e12.3)

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
c     write (6,*) td,rd,xd,omtd,omrd,omxd

      abar = abartab(it1r1c1)*omtd*omrd*omxd
     1       + abartab(it1r1c2)*omtd*omrd*xd
     2       + abartab(it1r2c1)*omtd*rd*omxd
     3       + abartab(it1r2c2)*omtd*rd*xd
     4       + abartab(it2r1c1)*td*omrd*omxd 
     5       + abartab(it2r1c2)*td*omrd*xd
     6       + abartab(it2r2c1)*td*rd*omxd
     7       + abartab(it2r2c2)*td*rd*xd

      dq   = ebtab(it1r1c1)*omtd*omrd*omxd
     1       + ebtab(it1r1c2)*omtd*omrd*xd
     2       + ebtab(it1r2c1)*omtd*rd*omxd
     3       + ebtab(it1r2c2)*omtd*rd*xd
     4       + ebtab(it2r1c1)*td*omrd*omxd 
     5       + ebtab(it2r1c2)*td*omrd*xd
     6       + ebtab(it2r2c1)*td*rd*omxd
     7       + ebtab(it2r2c2)*td*rd*xd

      dyedt = wratetab(it1r1c1)*omtd*omrd*omxd
     1       + wratetab(it1r1c2)*omtd*omrd*xd
     2       + wratetab(it1r2c1)*omtd*rd*omxd
     3       + wratetab(it1r2c2)*omtd*rd*xd
     4       + wratetab(it2r1c1)*td*omrd*omxd 
     5       + wratetab(it2r1c2)*td*omrd*xd
     6       + wratetab(it2r2c1)*td*rd*omxd
     7       + wratetab(it2r2c2)*td*rd*xd

      return
      end
