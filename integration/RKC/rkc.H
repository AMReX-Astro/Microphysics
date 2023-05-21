#ifndef RKC_H
#define RKC_H

#include <rkc_type.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rkc (burn_t& state, rkc_t rstate)
{

    const Real rmax = 0.1_rt;
    const Real rmin = 10.0_rt * UROUND;

    // the original rkc.f checked on the input value of idid, but we
    // only work for idid = 0, which is a new integration (not a
    // continuation)

    // Test the input data.
    bool valid = true;

    if ((rtol_spec > rmax) || (rtol_enuc > rmax) ||
        (rtol_spec < rmin) || (rtol_spec < rmin)) {
        valid = false;
    }

    if (atol_spec < 0.0_rt || atol_enuc < 0.0_rt) {
        valid = false;
    }

    if (! valid) {
        rstate.idid = 5;  // TODO: remove idid from rstate
        return rstate.idid;
    }

    // Initialize counters and pointers.

    rstate.nfe = 0;
    rstate.nsteps = 0;
    rstate.naccpt = 0;
    rstate.nrejct = 0;
    rstate.nfesig = 0;
    rstate.maxm = 0;

    rkclow(state, rstate);

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rkclow (burn_t& state, rkc_t rstate)
{


    // RKC is an interface to RKCLOW where the actual solution takes place.

    const Real onep1 = 1.1_rt;
    const Real onep54 = 1.54_rt;
    const Real p1 = 0.1_rt;
    const Real p4 = 0.4_rt;
    const Real p8 = 0.8_rt;
    const Real one3rd = 1.0_rt/3.0_rt;
    const Real two3rd = 2.0_rt/3.0_rt;

    // Initialize on the first call.

    rstate.mmax = std::round(std::sqrt(rtol/(10.0_rt * uround)));
    rstate.mmax = std::max(rstate.mmax, 2);
    bool newspc = true;
    bool jacatt = false;
    int nstsig = 0;
    for (int i = 1; i <= RKC_NEQS; ++i) {
        rstate.work.yn(i) = rstate.y(i);
    }

    call f(neqn,t,yn,fn);

    rstate.nfe++;
    int tdir = std::copysign(1.0_rt, tend - t);
    rstate.hmax = std::abs(tend - t);

    int hmin = 10.0_rt * uround * std::max(std::abs(t), hmax);

c------------------------------------      
c  Start of loop for taking one step.  
c------------------------------------
20    continue      
c----------------------------------------------
c  Estimate the spectral radius of the Jacobian 
c  when newspc = .true..  A convergence failure 
c  in rkcrho is reported by idid = 6.
c----------------------------------------------
      if(newspc) then
        if(info(2) .eq. 1) then
          sprad = spcrad(neqn,t,yn)
        else
          call rkcrho(neqn,t,f,yn,fn,vtemp1,vtemp2,work,sprad,idid)
          if(idid .eq. 6) return
        endif
        jacatt = .true.
      endif
c-------------------------------
c  Compute an initial step size.
c-------------------------------        
      if(nsteps .eq. 0) then
        absh = hmax
        if(sprad*absh .gt. one) absh = one/sprad
        absh = max(absh,hmin)
        do 30 i = 1,neqn
          vtemp1(i) = yn(i) + absh*fn(i)
30      continue   
        call f(neqn,t+absh,vtemp1,vtemp2)
        nfe = nfe + 1
        est = zero
        at = atol(1)
        do 40 i = 1,neqn
          if(array) at = atol(i)
          wt = at + rtol*abs(yn(i))
          if(wt .eq. zero) then
            idid = 3
            return
          endif
          est = est + ((vtemp2(i) - fn(i))/wt)**2
40      continue          
        est = absh*sqrt(est/neqn)
        if(p1*absh .lt. hmax*sqrt(est)) then
          absh = max(p1*absh/sqrt(est), hmin)
        else 
          absh = hmax
        endif
      endif        
c------------------------------------------------------------
c  Adjust the step size and determine the number of stages m.
c------------------------------------------------------------      
      last = .false.
      if(onep1*absh .ge. abs(tend - t)) then
        absh = abs(tend - t)
        last = .true.
      endif
      m = 1 + int(sqrt(onep54*absh*sprad + one))
c----------------------------------------------------------
c  Limit m to mmax to control the growth of roundoff error.
c----------------------------------------------------------
      if(m .gt. mmax) then
        m = mmax
        absh = (m**2 - 1)/(onep54*sprad)
        last = .false.
      endif
      maxm = max(m,maxm)
c--------------------------------------------
c  A tentative solution at t+h is returned in
c  y and its slope is evaluated in vtemp1(*).
c--------------------------------------------
      h = tdir*absh
      hmin = ten*uround*max(abs(t),abs(t + h))
      call step(neqn,f,t,yn,fn,h,m,y,vtemp1,vtemp2)
      call f(neqn,t+h,y,vtemp1)
      nfe = nfe + m
      nsteps = nsteps + 1
c-------------------------------------------------------------
c  Estimate the local error and compute its weighted RMS norm.
c-------------------------------------------------------------
      err = zero
      at = atol(1)
      do 50 i = 1, neqn
        if(array) at = atol(i)
        wt = at + rtol*max(abs(y(i)),abs(yn(i)))
        if(wt .eq. zero) then
          idid = 3
          return
        endif
        est = p8*(yn(i) - y(i)) + p4*h*(fn(i) + vtemp1(i))
        err = err + (est/wt)**2
50    continue
      err = sqrt(err/neqn)
c
      if(err .gt. one) then
c-------------------
c  Step is rejected.
c-------------------
        nrejct = nrejct + 1
        absh = p8*absh/(err**one3rd)
        if(absh .lt. hmin) then
          idid = 4
          return
        else
          newspc = .not. jacatt
          goto 20
        endif
      endif
c-------------------
c  Step is accepted.
c-------------------
      naccpt = naccpt + 1
      t = t + h
      jacatt = info(3) .eq. 1
      nstsig = mod(nstsig+1,25)
      newspc = .false.
      if((info(2) .eq. 1) .or. (nstsig .eq. 0)) newspc = .not. jacatt
c------------------------------------------------------
c  Update the data for interpolation stored in work(*).
c------------------------------------------------------      
      work(1) = h
      work(2) = t
      do 60 i = 1, neqn
         ylast = yn(i)
         yplast = fn(i)
         yn(i) = y(i)
         fn(i) = vtemp1(i)
         vtemp1(i) = ylast
         vtemp2(i) = yplast
60    continue
      fac = ten
      if(naccpt .eq. 1) then
        temp2 = err**one3rd
        if(p8 .lt. fac*temp2) fac = p8/temp2
      else
        temp1 = p8*absh*errold**one3rd
        temp2 = abs(hold)*err**two3rd
        if(temp1 .lt. fac*temp2) fac = temp1/temp2
      endif
      absh = max(p1,fac)*absh
      absh = max(hmin,min(hmax,absh))
      errold = err
      hold = h
      h = tdir*absh
      if(last) then
        idid = 1
        return
      elseif(info(1) .eq. 0) then
        idid = 2
        return
      else
        goto 20
      endif
      end

      subroutine step(neqn,f,t,yn,fn,h,m,y,yjm1,yjm2)
c--------------------------------------------------
c  Take a step of size H from T to T+H to get Y(*).
c--------------------------------------------------
      integer          neqn,m
      double precision t,yn(neqn),fn(neqn),h,
     &                 y(neqn),yjm1(neqn),yjm2(neqn)
      external         f
c      
      double precision one,two,four,c13,zero
      parameter       (one=1d0,two=2d0,four=4d0,c13=13d0,zero=0d0)
      integer          i,j
      double precision ajm1,arg,bj,bjm1,bjm2,dzj,dzjm1,dzjm2,
     &                 d2zj,d2zjm1,d2zjm2,mu,mus,nu,
     &                 temp1,temp2,thj,thjm1,thjm2,w0,w1,
     &                 zj,zjm1,zjm2
c
      w0 = one + two/(c13*m**2)
      temp1 = w0**2 - one
      temp2 = sqrt(temp1)
      arg = m*log(w0 + temp2)
      w1 = sinh(arg)*temp1 / (cosh(arg)*m*temp2 - w0*sinh(arg))
      bjm1 = one/(two*w0)**2
      bjm2 = bjm1
c---------------------------
c  Evaluate the first stage.
c---------------------------
      do 10 i = 1, neqn
        yjm2(i) = yn(i)
10    continue
      mus = w1*bjm1
      do 20 i = 1, neqn
        yjm1(i) = yn(i) + h*mus*fn(i)
20    continue
      thjm2  = zero
      thjm1  = mus
      zjm1   = w0
      zjm2   = one
      dzjm1  = one
      dzjm2  = zero
      d2zjm1 = zero
      d2zjm2 = zero
c------------------------------
c  Evaluate stages j = 2,...,m.
c------------------------------
      do 50 j = 2, m
        zj   =   two*w0*zjm1 - zjm2
        dzj  =   two*w0*dzjm1 - dzjm2 + two*zjm1
        d2zj =   two*w0*d2zjm1 - d2zjm2 + four*dzjm1
        bj   =   d2zj/dzj**2
        ajm1 =   one - zjm1*bjm1
        mu   =   two*w0*bj/bjm1
        nu   = - bj/bjm2
        mus  =   mu*w1/w0
c---------------------------------------------
c  Use the y array for temporary storage here.         
c---------------------------------------------
        call f(neqn,t + h*thjm1,yjm1,y)
        do 30 i = 1, neqn
          y(i) = mu*yjm1(i) + nu*yjm2(i) + (one - mu - nu)*yn(i) +
     &           h*mus*(y(i) - ajm1*fn(i))
30      continue
        thj = mu*thjm1 + nu*thjm2 + mus*(one - ajm1)
c------------------------------------
c  Shift the data for the next stage.
c------------------------------------
        if(j .lt. m) then
          do 40 i = 1, neqn
            yjm2(i) = yjm1(i)
            yjm1(i) = y(i)
40        continue
          thjm2  = thjm1
          thjm1  = thj
          bjm2   = bjm1
          bjm1   = bj
          zjm2   = zjm1
          zjm1   = zj
          dzjm2  = dzjm1
          dzjm1  = dzj
          d2zjm2 = d2zjm1
          d2zjm1 = d2zj
        endif
50    continue
      return
      end


      subroutine rkcrho(neqn,t,f,yn,fn,v,fv,work,sprad,idid)
c---------------------------------------------------------------
c  RKCRHO attempts to compute a close upper bound, SPRAD, on
c  the spectral radius of the Jacobian matrix using a nonlinear
c  power method.  A convergence failure is reported by IDID = 6.
c---------------------------------------------------------------
      integer          neqn,idid
      double precision t,yn(neqn),fn(neqn),v(neqn),fv(neqn),work(*),
     &                 sprad
      external         f
c      
      integer          itmax
      parameter       (itmax=50)
      double precision zero,one,onep2,p01
      parameter       (zero=0d0,one=1d0,onep2=1.2d0,p01=0.01d0)
      integer          i,iter,index,ptr5
      double precision uround,sqrtu,ynrm,sigma,sigmal,
     &                 dynrm,dfnrm,vnrm,small
      integer          nfe,nsteps,naccpt,nrejct,nfesig,maxm
      common /rkcdid/  nfe,nsteps,naccpt,nrejct,nfesig,maxm
c
c------------------------------------------------------------
c  hmax = work(6).  sprad smaller than small = 1/hmax are not
c  interesting because they do not constrain the step size.  
c------------------------------------------------------------
      small = one/work(6)
c---------------------------------------------------------
c  The initial slope is used as guess when nsteps = 0 and
c  thereafter the last computed eigenvector.  Some care
c  is needed to deal with special cases. Approximations to
c  the eigenvector are normalized so that their Euclidean
c  norm has the constant value dynrm.
c---------------------------------------------------------
      ptr5 = nint(work(7)) + 4*neqn
      if(nsteps .eq. 0) then
        do 10 i = 1,neqn
          v(i) = fn(i)
10      continue
      else
        do 20 i = 1,neqn          
          v(i) = work(ptr5+i-1)
20      continue
      endif
      ynrm = zero
      vnrm = zero
      do 30 i = 1,neqn
        ynrm = ynrm + yn(i)**2     
        vnrm = vnrm + v(i)**2
30    continue
      ynrm = sqrt(ynrm)
      vnrm = sqrt(vnrm)
      if(ynrm .ne. zero .and. vnrm .ne. zero) then
        dynrm = ynrm*sqrtu
        do 40 i = 1,neqn
          v(i) = yn(i) + v(i)*(dynrm/vnrm)
40      continue          
      elseif(ynrm .ne. zero) then
        dynrm = ynrm*sqrtu
        do 50 i = 1, neqn
          v(i) = yn(i) + yn(i)*sqrtu
50      continue
      elseif(vnrm .ne. zero) then
        dynrm = uround
        do 60 i = 1,neqn
          v(i) = v(i)*(dynrm/vnrm)
60      continue
      else          
        dynrm = uround
        do 70 i = 1,neqn
          v(i) = dynrm
70      continue          
      endif
c--------------------------------------------
c  Now iterate with a nonlinear power method.
c--------------------------------------------
      sigma = zero
      do 110 iter = 1, itmax
        call f(neqn,t,v,fv)
        nfesig = nfesig + 1
        dfnrm = zero
        do 80 i = 1, neqn
          dfnrm = dfnrm + (fv(i) - fn(i))**2
80      continue   
        dfnrm = sqrt(dfnrm)
        sigmal = sigma
        sigma = dfnrm/dynrm
c----------------------------------------------------------
c  sprad is a little bigger than the estimate sigma of the 
c  spectral radius, so is more likely to be an upper bound.
c----------------------------------------------------------
        sprad = onep2*sigma
        if(iter .ge. 2 .and.
     &     abs(sigma - sigmal) .le. max(sigma,small)*p01) then
          do 90 i = 1,neqn          
            work(ptr5+i-1) = v(i) - yn(i)
90        continue
          return
        endif
c--------------------------------------
c  The next v(*) is the change in f
c  scaled so that norm(v - yn) = dynrm.
c--------------------------------------
        if(dfnrm .ne. zero) then
          do 100 i = 1,neqn
            v(i) = yn(i) + (fv(i) - fn(i))*(dynrm/dfnrm)
100       continue           
        else
c-------------------------------------------------------
c  The new v(*) degenerated to yn(*)--"randomly" perturb 
c  current approximation to the eigenvector by changing 
c  the sign of one component.
c-------------------------------------------------------
          index = 1 + mod(iter,neqn)
          v(index) = yn(index) - (v(index) - yn(index))
        endif
110   continue
c-------------------------------------------
c  Set flag to report a convergence failure.
c-------------------------------------------
      idid = 6
      return 
      end

#endif
