#ifndef ROCK4_H
#define ROCK4_H

#include <rock4_type.H>
#include <burn_type.H>
#ifdef STRANG
#include <integrator_type_strang.H>
#include <integrator_rhs_strang.H>
#endif
#ifdef SIMPLIFIED_SDC
#include <integrator_type_simplified_sdc.H>
#include <integrator_rhs_simplified_sdc.H>
#endif

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rock4 (BurnT& state, rkc_t& rstate)
{
    // initializations

    // TODO: get the initial set size

    // Test the initial step size and tolerances.

    if (h.gt.abs(tend-t)) then
        write(6,*)'initial step is longer than the integration interval'
        arret=.true.
        idid=-1
        return
      end if  
      if (h.lt.10.d0*uround) then
        write(6,*)'initial step-size is too small'
        idid=-1
        arret=.true.
        return
      end if 

c -------- Call to the core integrator.-----------
      call rockfcore(neqn,t,tend,h,y,f,work,work(n1),work(n2),
     & work(n3),work(n4),work(n5),work(n6),work(n7),ms,atol,rtol,
     & ntol,recf,fpa,fpb,fpbe,iwork,arret,uround,idid)
      return
      end
c ----------------------------------------------
c     End of subroutine ROCK4.
c ----------------------------------------------
c
       subroutine rockfcore(neqn,t,tend,h,y,f,work,yn,fn,fnt,yjm1,
     & yjm2,yjm3,yjm4,ms,atol,rtol,ntol,recf,fpa,fpb,fpbe,iwork,
     & arret,uround,idid)
c ----------------------------------------------
c    Core integrator for ROCK4.
c ---------------------------------------------- 
c             Declarations
c ----------------------------------------------
       double precision y(*),yn(*),fn(*),fnt(neqn),work(*),
     & yjm1(neqn),yjm2(neqn),yjm3(neqn),yjm4(neqn),recf(4382),
     & atol(*),rtol(*),fpa(50,6),fpb(50,4),fpbe(50,5),err,errp,
     & eigmax,t,tend,h,hnew,hp,facmax,fac,facp,rho,uround,told
       integer ms(50),mp(2),iwork(12),neqn,mdeg,i,ntol,
     & nrho,mdego,nrej,idid
       logical last,reject,arret
       external f
c *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** 
c             initializations
c *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***  

      iwork(5)=0
      iwork(6)=100000000
      iwork(7)=0
      iwork(8)=0
      iwork(9)=0
      iwork(10)=0      
      iwork(11)=0
      iwork(12)=0

      facmax=5.d0
      told=0.d0
      nrej=0
      hp=h
      mdego=0
      last=.false.
      reject=.false. 
      err=0.d0
      nrho=0
      idid=1
c -------- initialization of the integration step.-------- 
10    if (rstate.nsteps != 0) {
        do i=1,neqn
          yn(i)=y(i)
          fn(i)=fnt(i)
        end do 
      } else {
        do i=1,neqn
          yn(i)=y(i)
        end do
      call f(neqn,t,yn,fn) 
      }
      rstate.nsteps++;
      errp=err

          // Step-size is adjusted.
20    if(1.01d0*h.ge.abs(tend-t)) then
        h=dabs(tend-t)
        last=.true.
      end if
      if (h.lt.10.d0*uround) then
        write(6,*)' tolerances are too small'
        idid=-2
        arret=.true.
        return
      end if
          // Spectral radius.
      if (nrho == 0) {
          // Computed internally by rockfrho
            call rockfrho(neqn,t,y,f,yn,fn,work,yjm1,
     &                   yjm2,eigmax,uround,idid,iwork)
            rstate.max_rho = std::max(rstate.max_rho, std::round(eigmax));
            rstate.min_rho = std::min(rstate.min_rho, std::round(eigmax));
      }

c -------- The number of stages.-------- 
       mdeg=sqrt((3.d0+h*eigmax)/0.353d0)+1
      if (mdeg.gt.152) then
        h=0.8d0*(152.d0**2*0.353d0-3.d0)/eigmax
        mdeg=152
        last=.false.
      end if
      mdeg=max(mdeg,5)-4
      if (mdeg.ne.mdego) then
        call mdegre(mdeg,mp,ms)
      end if
      if (mdeg+4 > rstate.max_stages) {
          rstate.max_stages = mdeg + 4;
      }

c -------- Computation of an integration step.--------
      call rfstep(neqn,t,h,y,f,yn,fn,fnt,yjm1,yjm2,yjm3,yjm4,
     & mdeg,mp,err,atol,rtol,ntol,recf,fpa,fpb,fpbe)
      mdego=mdeg
      rstate.nsteps++;
      rstate.nfe += mdeg + 4;

c -------- Error control procedure.--------
      fac=(1.d0/err)**(0.25d0)
      if (errp.ne.0.d0.and..not.reject) then
        facp=errp**(0.25d0)*fac**2*(h/hp)
        fac=dmin1(fac,facp)
      end if
      if (reject) then
        facmax=1.d0
      end if  
      fac=dmin1(facmax,dmax1(0.1d0,0.8d0*fac))
      hnew=h*fac
c -------- Accepted step.--------
      if (err.lt.1.d0) then
        rstate.naccpt++;
        facmax=2.d0
        t=t+h
        if (reject) then
          hnew=dmin1(hnew,h)
          if (tend.lt.t)  hnew=dmax1(hnew,h)
          reject= .false.
          nrej=0
        end if
        hp=h
        h=hnew
        nrho=nrho+1
        nrho=mod(nrho+1,10)
        if (last) then 
          return
        else
          goto 10
        end if
                   else            
c -------- Rejected step.--------
        rstate.nrejct++;
        reject= .true.
        last=.false.
        h= 0.8d0*hnew
        if (rstate.nsteps == 0) {
            h = 0.1e0 * h;
        }
        if (told.eq.t) then 
          nrej=nrej+1
          if (nrej.eq.10) h=1.0d-5
        end if
        told=t      
c -------The spectral radius is recomputed.-------- 
c        after a step failure
        if (nrho.ne.0) then 
          nrho=0
                       else
          nrho=1
        end if
        goto 20
      end if        
      return      
      end
c ----------------------------------------------
c     End of subroutine rockcore.
c ---------------------------------------------- 
c        
      subroutine rfstep(neqn,t,h,y,f,yn,fn,fnt,yjm1,yjm2,yjm3,
     & yjm4,mdeg,mp,err,atol,rtol,ntol,recf,fpa,fpb,fpbe)
c ----------------------------------------------
c  Solut. at t+h by an explicit (mdeg+4)-stages formula.
c-----------------------------------------------
c ---------------------------------------------- 
c             Declarations
c ----------------------------------------------
       double precision y(neqn),yn(neqn),fn(neqn),fnt(neqn),
     & yjm2(neqn),yjm3(neqn),yjm4(neqn),recf(4382),fpa(50,6),
     & yjm1(neqn),fpb(50,4),fpbe(50,5),atol(*),rtol(*),err,t,
     & h,ci1,ci2,ci3,temp1,temp2,temp3,temp4,temp5,ato,rto
       integer mp(2),neqn,mdeg,mr,mz,i,j,ntol
       external f
c *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** 
c             initializations
c *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***
      err=0.d0
      mz=mp(1)
      mr=mp(2)
c -------- First stage.--------
      temp1=h*recf(mr)
      ci1=t+temp1
      ci2=t+temp1
      ci3=t
      do i=1,neqn
        yjm2(i)=yn(i)
        yjm1(i)=yn(i)+temp1*fn(i)
        if (mdeg.lt.2) then
          y(i)=yjm1(i)
        end if
      end do
c --------- Stage for j=2..mdeg.--------
      do i=2,mdeg
        temp1=h*recf(mr+2*(i-2)+1)
        temp3=-recf(mr+2*(i-2)+2)
        temp2=1.d0-temp3
        call f(neqn,ci1,yjm1,y)
        ci1=temp1+temp2*ci2+temp3*ci3
        do j=1,neqn
          y(j)=temp1*y(j)+temp2*yjm1(j)+temp3*yjm2(j)
c -------- Shift the value "y" for the next stage.--------
          if (i.lt.mdeg) then
            yjm2(j)=yjm1(j)
            yjm1(j)=y(j)
          end if 
        end do
        ci3=ci2
        ci2=ci1
      end do          
c -------- The finishing procedure (4-stage method).
c -------- Stage 1.--------  
      temp1=h*fpa(mz,1)
      call f(neqn,ci1,y,yjm1)
      do j=1,neqn
        yjm3(j)=y(j)+temp1*yjm1(j)
      end do
c -------- Stage 2.--------    
      ci2=ci1+temp1
      temp1=h*fpa(mz,2)
      temp2=h*fpa(mz,3)
      call f(neqn,ci2,yjm3,yjm2)
      do j=1,neqn
        yjm4(j)=y(j)+temp1*yjm1(j)+temp2*yjm2(j)
      end do
c -------- Stage 3.--------
      ci2=ci1+temp1+temp2
      temp1=h*fpa(mz,4)
      temp2=h*fpa(mz,5)
      temp3=h*fpa(mz,6)
      call f(neqn,ci2,yjm4,yjm3)
      do j=1,neqn
        fnt(j)=y(j)+temp1*yjm1(j)+temp2*yjm2(j)+temp3*yjm3(j)
      end do
c -------- Stage 4.--------
      ci2=ci1+temp1+temp2+temp3
      temp1=h*fpb(mz,1)
      temp2=h*fpb(mz,2)
      temp3=h*fpb(mz,3)
      temp4=h*fpb(mz,4)
      call f(neqn,ci2,fnt,yjm4)
      do j=1,neqn
        y(j)=y(j)+temp1*yjm1(j)+temp2*yjm2(j)+temp3*yjm3(j)
     &  +temp4*yjm4(j)
      end do
c -------- Error evaluation (embedded method of order 3).--------
      temp1=h*fpbe(mz,1)-temp1
      temp2=h*fpbe(mz,2)-temp2
      temp3=h*fpbe(mz,3)-temp3
      temp4=h*fpbe(mz,4)-temp4
      temp5=h*fpbe(mz,5)
      call f(neqn,t+h,y,fnt)
c -------- Atol and rtol are scalar.--------
      if (ntol.eq.0) then
        ato=atol(1)
        rto=rtol(1)
        do j=1,neqn
          ci1=dabs(y(j))*rto
          err=err+((temp1*yjm1(j)+temp2*yjm2(j)+temp3*yjm3(j)
     &   +temp4*yjm4(j)+temp5*fnt(j))/(ato+ci1))**2
        end do
c -------- Atol and rtol are array.--------
                     else
        do j=1,neqn
          ci1=dabs(y(j))*rtol(j)
          err=err+((temp1*yjm1(j)+temp2*yjm2(j)+temp3*yjm3(j)
     &   +temp4*yjm4(j)+temp5*fnt(j))/(atol(j)+ci1))**2
        end do
      end if
      err=sqrt(err/neqn)
      return 
      end
c ----------------------------------------------
c     End of subroutine rfstep.
c ----------------------------------------------   
c    
      subroutine mdegre(mdeg,mp,ms)
c-------------------------------------------------------------          
c       Find the optimal degree.
c       MP(1): pointer which select the degree in ms(i)\1,2,..
c             such that mdeg<=ms(i).
c       MP(2): pointer which gives the corresponding position
c       of a_1 in the data recf for the selected degree.        
c------------------------------------------------------------- 
c ---------------------------------------------- 
c             Declarations
c ----------------------------------------------
      integer ms(50),mp(2),mdeg,i
c *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** 
c             initializations 
c *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***     
      mp(2)=1
c -------- Find the degree.--------      
      do i=1,50
        if ((ms(i)/mdeg).ge.1) then
          mdeg=ms(i)
          mp(1)=i
          return
        end if
      mp(2)=mp(2)+ms(i)*2-1
      end do
      return
      end   
c ----------------------------------------------
c     End of subroutine mdegr.
c ----------------------------------------------
c   
      subroutine rockfrho(neqn,t,y,f,yn,fn,work,z,fz,eigmax,
     &                   uround,idid,iwork)
c------------------------------------------------------------ 
c     Rockfrho compute eigmax, a close upper bound of the
c     spectral radius of the Jacobian matrix using a 
c     power method (J.N. Franklin (matrix theory)). 
c     The algorithm used is a small change (initial vector
c     and stopping criteria) of that of
c     Sommeijer-Shampine-Verwer, implemented in RKC.
c-------------------------------------------------------------
c             Declarations
c-------------------------------------------------------------
       double precision y(neqn),yn(neqn),fn(neqn),z(neqn),
     & work(*),fz(neqn),t,eigmax,eigmaxo,sqrtu,uround,znor,
     & ynor,quot,dzyn,dfzfn,safe
       integer iwork(12),neqn,n8,i,iter,maxiter,nind,ntest,
     & ind,idid  
       parameter (maxiter=50)
       parameter (safe=1.2d0)
       external f
c *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** 
c             initializations
c *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***     
      sqrtu=sqrt(uround)
      ynor=0.d0
      znor=0.d0
      n8=7*neqn
c ------ The initial vectors for the power method are yn --------
c       and yn+c*f(v_n), where vn=f(yn) a perturbation of yn 
c       (if first step) or a perturbation of the last computed
c       eigenvector (otherwise)
c
      if (rstate.nsteps == 0) {
        do i=1,neqn
          fz(i)=fn(i)
        end do
        call f(neqn,t,fz,z)
        rstate.nfesig++;
            } else {
        do i=1,neqn
          z(i)=work(n8+i)
        end do
      }
c ------ Perturbation.--------
      do i=1,neqn
        ynor=ynor+yn(i)**2
        znor=znor+z(i)**2
      end do
      ynor=sqrt(ynor)
      znor=sqrt(znor)
c ------ Normalization of the vector z so that --------
c        the difference z-yn lie in a circle 
c        around yn (i.e has a constant modulus).
c
      if (ynor.ne.0.d0.and.znor.ne.0.d0) then
        dzyn=ynor*sqrtu
        quot=dzyn/znor
        do i=1,neqn
          z(i)=yn(i)+z(i)*quot
        end do
      elseif(ynor.ne.0.d0) then
        dzyn=ynor*sqrtu
        do i=1,neqn
          z(i)=yn(i)+yn(i)*sqrtu
        end do
      elseif(znor.ne.0.d0) then
          dzyn=uround
          quot=dzyn/znor
          do i=1,neqn
            z(i)=z(i)*quot
          end do
      else
        dzyn=uround
        do i=1,neqn
          z(i)=dzyn
        end do
      end if
c ------ Start the power method.--------
      eigmax=0.d0
      do iter=1,maxiter
        call f(neqn,t,z,fz)
        rstate.nfesig++;
        dfzfn=0.d0
        do i=1,neqn
          dfzfn=dfzfn+(fz(i)-fn(i))**2
        end do
        dfzfn=dsqrt(dfzfn)
        eigmaxo=eigmax
        eigmax=dfzfn/dzyn
        eigmax=safe*eigmax
c ------ The stopping criteria is based on a 
c        relative error between two succesive
c        estimation ``eigmax'' of the spectral
c        radius. 
c
        if (iter.ge.2.and.dabs(eigmax-eigmaxo)
     &    .le.(eigmax*0.05d0)) then
c ----- The last eigenvector is stored.-------- 
          do i=1,neqn
            work(n8+i)=z(i)-yn(i)
          end do
          return
        end if
c ----- The next z is defined by -------- 
c       z_new=yn+coef*(fz-fn) where
c       coef is choosen so that
c       norm(z_new-yn)=norm(z_old-yn)
c
        if (dfzfn.ne.0.d0) then
          quot=dzyn/dfzfn
          do i=1,neqn
            z(i)=yn(i)+(fz(i)-fn(i))*quot
          end do
        else
c ----- The new z is defined by an arbitrary --------
c       perturbation of the current approximation
c       of the eigenvector.
c
          nind=neqn
          ntest=0
          ind=1+mod(iter,nind)
          if (z(ind).ne.yn(ind).or.ntest.eq.10) then
            z(ind)=yn(ind)-(z(ind)-yn(ind))
          else 
            nind=neqn+ind
            ntest=ntest+1
          end if
        end if
      end do
      write(6,*) 'convergence failure in the 
     & spectral radius computation'
      idid=-3
      return
      end
         
      
        
        
              
      
     


