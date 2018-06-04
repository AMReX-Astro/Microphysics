!***************************************************************************************************
! jacobian_dense.f90 10/18/17
! The routines in this file assume a dense Jacobian and use a dense linear algebra package.
!
! The bulk of the computational cost of the network (60-95%) is the solving of the matrix equation.
! Careful selection of the matrix solver is therefore very important to fast computation. For
! networks from a few dozen up to a couple hundred species, hand tuned dense solvers such as those
! supplied by the hardware manufacturer (often LAPACK) or third-parties like NAG, IMSL, etc. are
! fastest. However for larger matrices, sparse solvers are faster.
!*******************************************************************************

Module jacobian_data
  !-------------------------------------------------------------------------------------------------
  ! The Jacobian matrix for the solver.
  !-------------------------------------------------------------------------------------------------
  Use, Intrinsic :: iso_c_binding, Only: C_DOUBLE, C_INT, C_INTPTR_T, C_PTR, C_SIZE_T
  Use xnet_types, Only: dp
  Implicit None

  ! CPU data (pointers for pinned memory)
  Real(dp), Pointer :: jac(:,:,:), rhs(:,:)
  Integer, Pointer :: indx(:,:), info(:)
  !$omp threadprivate(jac,rhs,indx,info)

  ! C pointers for pinned memory
  Type(C_PTR) :: hjac, hrhs, hinfo, hindx
  !$omp threadprivate(hjac,hrhs,hinfo,hindx)

  ! Device pointers for arrays
  Type(C_PTR) :: djac, drhs, dindx, dinfo
  Integer(C_INTPTR_T), Pointer :: djacf(:,:,:), drhsf(:,:)
  !$omp threadprivate(djac,drhs,dindx,dinfo,djacf,drhsf)

  ! Arrays of pointers to each device array batch element address
  Type(C_PTR), Allocatable, Target :: djaci(:), drhsi(:)
  !$omp threadprivate(djaci,drhsi)

  ! Host and device addresses for the arrays of device pointers
  Type(C_PTR) :: hdjac_array, hdrhs_array
  Type(C_PTR) :: djac_array, drhs_array
  !$omp threadprivate(hdjac_array,hdrhs_array,djac_array,drhs_array)

  ! Array size parameters
  Real(C_DOUBLE), Parameter :: ddummy = 0.0
  Integer(C_INT), Parameter :: idummy = 0
  Integer(C_INTPTR_T), Parameter :: cptr_dummy = 0
  Integer(C_SIZE_T), Parameter :: sizeof_double = sizeof(ddummy)
  Integer(C_SIZE_T), Parameter :: sizeof_int = sizeof(idummy)
  Integer(C_SIZE_T), Parameter :: sizeof_cptr = sizeof(cptr_dummy)
  Integer(C_SIZE_T) :: sizeof_jac, sizeof_rhs, sizeof_indx, sizeof_info

  ! Parameters for GPU array dimensions
  Integer :: msize ! Size of linear system to be solved

End Module jacobian_data

Subroutine read_jacobian_data(data_dir)
  !-------------------------------------------------------------------------------------------------
  ! Initializes the Jacobian data.
  !-------------------------------------------------------------------------------------------------
  Use, Intrinsic :: iso_c_binding, Only: c_f_pointer, c_loc
  Use controls, Only: iheat, nzbatchmx
  Use cublasf, Only: cublasSetVector
  Use cudaf, Only: cudaHostAlloc, cudaHostAllocDefault, cudaMalloc
  Use jacobian_data, Only: msize, sizeof_double, sizeof_int, sizeof_cptr, sizeof_jac, sizeof_rhs, &
    sizeof_indx, sizeof_info, jac, rhs, indx, info, hjac, hrhs, hindx, hinfo, djac, drhs, dindx, &
    dinfo, djacf, drhsf, djaci, drhsi, hdjac_array, hdrhs_array, djac_array, drhs_array
  Use nuc_number, Only: ny
  Use openaccf, Only: acc_map_data
  Implicit None

  ! Input variables
  Character(*), Intent(in) :: data_dir

  ! Local variables
  Integer :: istat, i

  ! Calculate array sizes
  If ( iheat > 0 ) Then
    msize = ny + 1
  Else
    msize = ny
  EndIf
  sizeof_jac = msize*msize*nzbatchmx*sizeof_double
  sizeof_rhs = msize*nzbatchmx*sizeof_double
  sizeof_indx = msize*nzbatchmx*sizeof_int
  sizeof_info = nzbatchmx*sizeof_int

  !$omp parallel default(shared) private(istat,i)

  ! Allocate CPU memory (pinned for asynchronous host<->device copy)
  istat = cudaHostAlloc(hjac, sizeof_jac, cudaHostAllocDefault)
  istat = cudaHostAlloc(hrhs, sizeof_rhs, cudaHostAllocDefault)
  istat = cudaHostAlloc(hindx, sizeof_indx, cudaHostAllocDefault)
  istat = cudaHostAlloc(hinfo, sizeof_info, cudaHostAllocDefault)
  Call c_f_pointer(hjac, jac, (/msize,msize,nzbatchmx/))
  Call c_f_pointer(hrhs, rhs, (/msize,nzbatchmx/))
  Call c_f_pointer(hindx, indx, (/msize,nzbatchmx/))
  Call c_f_pointer(hinfo, info, (/nzbatchmx/))

  ! Allocate GPU memory for arrays
  istat = cudaMalloc(djac, msize*msize*nzbatchmx*sizeof_double)
  istat = cudaMalloc(drhs, msize*nzbatchmx*sizeof_double)
  istat = cudaMalloc(dindx, msize*nzbatchmx*sizeof_int)
  istat = cudaMalloc(dinfo, nzbatchmx*sizeof_int)

  ! Setup fortran pointers to device-addresses of device arrays
  Call c_f_pointer(djac, djacf, (/msize,msize,nzbatchmx/))
  Call c_f_pointer(drhs, drhsf, (/msize,nzbatchmx/))

  ! Setup arrays of pointers to each device array batch element address
  Allocate (djaci(nzbatchmx))
  Allocate (drhsi(nzbatchmx))
  do i = 1, nzbatchmx
    ! Get the device-addresses for this batch element
    djaci(i) = c_loc(djacf(1,1,i))
    drhsi(i) = c_loc(drhsf(1,i))
  end do

  ! Get the host-addresses for the arrays of device-pointers
  hdjac_array = c_loc(djaci(1))
  hdrhs_array = c_loc(drhsi(1))

  ! Allocate GPU memory for batched GPU operations and copy array of pointers to device
  istat = cudaMalloc(djac_array, nzbatchmx*sizeof_cptr)
  istat = cudaMalloc(drhs_array, nzbatchmx*sizeof_cptr)
  istat = cublasSetVector(nzbatchmx, sizeof_cptr, hdjac_array, 1, djac_array, 1)
  istat = cublasSetVector(nzbatchmx, sizeof_cptr, hdrhs_array, 1, drhs_array, 1)

  ! Map CUDA allocated arrays to OpenACC
  Call acc_map_data(hjac,djac,sizeof_jac)
  Call acc_map_data(hrhs,drhs,sizeof_rhs)

  !$omp end parallel

  Return
End Subroutine read_jacobian_data

Subroutine jacobian_finalize
  !-------------------------------------------------------------------------------------------------
  ! Free the page-locked and device memory used in the dense solver.
  !-------------------------------------------------------------------------------------------------
  Use cudaf, Only: cudaFree, cudaFreeHost
  Use jacobian_data
  Use openaccf, Only: acc_unmap_data

  ! Local variables
  Integer :: istat

  !$omp parallel default(shared) private(istat)

  call acc_unmap_data(hjac)
  call acc_unmap_data(hrhs)
  istat = cudaFree(djac)
  istat = cudaFree(drhs)
  istat = cudaFree(dindx)
  istat = cudaFree(dinfo)
  Deallocate (djaci)
  Deallocate (drhsi)

  istat = cudaFree(djac_array)
  istat = cudaFree(drhs_array)
  istat = cudaFreeHost(hjac)
  istat = cudaFreeHost(hrhs)
  istat = cudaFreeHost(hindx)
  istat = cudaFreeHost(hinfo)

  !$omp end parallel

End Subroutine jacobian_finalize

Subroutine jacobian_build(diag,mult,mask_in)
  !-------------------------------------------------------------------------------------------------
  ! This routine calculates the reaction Jacobian matrix, dYdot/dY, and augments by multiplying all
  ! elements by mult and adding diag to the diagonal elements.
  !-------------------------------------------------------------------------------------------------
  Use abundances, Only: yt
  Use conditions, Only: cv
  Use controls, Only: iheat, idiag, lun_diag, nzbatchmx, szbatch, lzactive
  Use cross_sect_data, Only: dcsect1dt9, dcsect2dt9, dcsect3dt9
  Use jacobian_data, Only: jac, msize
  Use nuc_number, Only: ny
  Use nuclear_data, Only: mex
  Use reac_rate_data, Only: a1, a2, a3, b1, b2, b3, la, le, mu1, mu2, mu3, &
    & n11, n21, n22, n31, n32, n33
  Use timers, Only: xnet_wtime, start_timer, stop_timer, timer_jacob
  Use xnet_types, Only: dp
  Implicit None

  ! Input variables
  Real(dp), Intent(in) :: diag(:), mult

  ! Optional variables
  Logical, Optional, Target, Intent(in) :: mask_in(:)

  ! Local variables
  Integer :: i, j, i0, i1, i11, i21, i22, i31, i32, i33, j1, la1, le1, la2, le2, la3, le3, izb, izone
  Real(dp) :: dydotdt9(ny), dt9dotdy(msize), s1, s2, s3, r1, r2, r3, y11, y21, y22, y31, y32, y33
  Real(dp) :: dr1dt9, dr2dt9, dr3dt9
  Logical, Pointer :: mask(:)

  If ( present(mask_in) ) Then
    mask => mask_in
  Else
    mask => lzactive
  EndIf

  start_timer = xnet_wtime()
  timer_jacob = timer_jacob - start_timer

  !$acc kernels &
  !$acc copyin(diag,mask,b1,b2,b3, &
  !$acc dcsect1dt9,dcsect2dt9,dcsect3dt9,yt,cv) &
  !$acc present(jac, &
  !$acc mex,la,le,mu1,mu2,mu3,n11,n21,n22,n31,n32,n33,a1,a2,a3)

  ! Build the Jacobian
  !$acc loop independent
  Do izb = 1, nzbatchmx
    If ( mask(izb) ) Then
      jac(:,:,izb) = 0.0
      !$acc loop independent
      Do i0 = 1, ny
        la1 = la(1,i0) ; la2 = la(2,i0) ; la3 = la(3,i0)
        le1 = le(1,i0) ; le2 = le(2,i0) ; le3 = le(3,i0)
        Do j1 = la1, le1
            jac(i0,n11(j1),izb) = jac(i0,n11(j1),izb) + b1(j1,izb)
        EndDo
        Do j1 = la2, le2
            jac(i0,n21(j1),izb) = jac(i0,n21(j1),izb) + b2(j1,izb)*yt(n22(j1),izb)
            jac(i0,n22(j1),izb) = jac(i0,n22(j1),izb) + b2(j1,izb)*yt(n21(j1),izb)
        EndDo
        Do j1 = la3, le3
            jac(i0,n31(j1),izb) = jac(i0,n31(j1),izb) + b3(j1,izb)*yt(n32(j1),izb)*yt(n33(j1),izb)
            jac(i0,n32(j1),izb) = jac(i0,n32(j1),izb) + b3(j1,izb)*yt(n31(j1),izb)*yt(n33(j1),izb)
            jac(i0,n33(j1),izb) = jac(i0,n33(j1),izb) + b3(j1,izb)*yt(n31(j1),izb)*yt(n32(j1),izb)
        EndDo
      EndDo
    EndIf
  EndDo

  If ( iheat > 0 ) Then
    !$acc loop independent
    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        !$acc loop independent
        Do i0 = 1, ny
          la1 = la(1,i0) ; la2 = la(2,i0) ; la3 = la(3,i0)
          le1 = le(1,i0) ; le2 = le(2,i0) ; le3 = le(3,i0)
          s1 = 0.0
          Do j1 = la1, le1
            dr1dt9 = a1(j1)*dcsect1dt9(mu1(j1),izb)
            i11 = n11(j1)
            y11 = yt(i11,izb)
            s1 = s1 + dr1dt9 * y11
          EndDo
          s2 = 0.0
          Do j1 = la2, le2
            dr2dt9 = a2(j1)*dcsect2dt9(mu2(j1),izb)
            i21 = n21(j1)
            i22 = n22(j1)
            y21 = yt(i21,izb)
            y22 = yt(i22,izb)
            s2 = s2 + dr2dt9 * y21 * y22
          EndDo
          s3 = 0.0
          Do j1 = la3, le3
            dr3dt9 = a3(j1)*dcsect3dt9(mu3(j1),izb)
            i31 = n31(j1)
            i32 = n32(j1)
            i33 = n33(j1)
            y31 = yt(i31,izb)
            y32 = yt(i32,izb)
            y33 = yt(i33,izb)
            s3 = s3 + dr3dt9 * y31 * y32 * y33
          EndDo
          jac(i0,ny+1,izb)= s1 + s2 + s3
        EndDo
        !$acc loop independent
        Do i1 = 1, msize
          s1 = 0.0
          Do i0 = 1, ny
            s1 = s1 + mex(i0)*jac(i0,i1,izb)
          EndDo
          jac(ny+1,i1,izb) = -s1 / cv(izb)
        EndDo
        ! The BLAS version of jac(ny+1,i,izb) = -sum(mex*jac(1:ny,i,izb))/cv(izb)
!       Call dgemv('T',ny,msize,1.0/cv(izb),jac(1,1,izb),msize,mex,1,0.0,dt9dotdy,1)
!       jac(ny+1,:,izb) = -dt9dotdy
      EndIf
    EndDo

    ! This also works, but could be inefficient if there are few active zones in batch
!   Call dgemv('T',ny,msize*nzbatchmx,1.0,jac(1,1,1),msize,mex,1,0.0,dt9dotdy,1)
!   ForAll ( izb = 1:nzbatchmx, j1 = 1:msize, mask(izb) )
!     jac(ny+1,j1,izb) = -dt9dotdy(j1,izb) / cv(izb)
!   EndForAll
  EndIf

  ! Apply the externally provided factors
  !$acc loop independent
  Do izb = 1, nzbatchmx
    If ( mask(izb) ) Then
      jac(:,:,izb) = mult * jac(:,:,izb)
      !$acc loop independent
      Do i0 = 1, msize
        jac(i0,i0,izb) = jac(i0,i0,izb) + diag(izb)
      EndDo
    EndIf
  EndDo
  !$acc end kernels

  If ( idiag >= 5 ) Then
    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        izone = izb + szbatch - 1
        Write(lun_diag,"(a9,i5,2es14.7)") 'JAC_BUILD',izone,diag(izb),mult
        Write(lun_diag,"(14es9.1)") ((jac(i,j,izb),j=1,msize),i=1,msize)
      EndIf
    EndDo
  EndIf

  stop_timer = xnet_wtime()
  timer_jacob = timer_jacob + stop_timer

  Return
End Subroutine jacobian_build

Subroutine jacobian_solve(kstep,yrhs,dy,t9rhs,dt9,mask_in)
  !-------------------------------------------------------------------------------------------------
  ! This routine solves the system of equations composed of the Jacobian and RHS vector.
  !-------------------------------------------------------------------------------------------------
  Use controls, Only: idiag, iheat, lun_diag, nzbatchmx, szbatch, lzactive
  Use cublasf, Only: cublasDgetrfBatched, cublasDgetrsBatched, cublasGetVectorAsync, &
    & cublasSetMatrixAsync, cublasSetVectorAsync
  Use cudaf, Only: cudaStreamSynchronize
  Use gpu_controls, Only: handle, stream
  Use jacobian_data
  Use nuc_number, Only: ny
  Use timers, Only: xnet_wtime, start_timer, stop_timer, timer_solve
  Use xnet_types, Only: dp
  Implicit None

  ! Input variables
  Integer, Intent(in) :: kstep
  Real(dp), Intent(in) :: yrhs(:,:)
  Real(dp), Intent(in) :: t9rhs(:)

  ! Output variables
  Real(dp), Intent(out) :: dy(size(yrhs,1),size(yrhs,2))
  Real(dp), Intent(out) :: dt9(size(t9rhs))

  ! Optional variables
  Logical, Optional, Target, Intent(in) :: mask_in(:)

  ! Local variables
  Integer :: i, izb, izone, istat
  Logical, Pointer :: mask(:)

  If ( present(mask_in) ) Then
    mask => mask_in
  Else
    mask => lzactive
  EndIf

  start_timer = xnet_wtime()
  timer_solve = timer_solve - start_timer

  Do izb = 1, nzbatchmx
    If ( mask(izb) ) Then
      rhs(1:ny,izb) = yrhs(:,izb)
      If ( iheat > 0 ) rhs(ny+1,izb) = t9rhs(izb)
    EndIf
  EndDo

  ! Copy the system to the GPU
  istat = cublasSetVectorAsync(msize*nzbatchmx, sizeof_double, hrhs, 1, drhs, 1, stream)

  ! Solve the linear system
  istat = cublasDgetrfBatched(handle, msize, djac_array, msize, dindx, dinfo, nzbatchmx)
  istat = cublasDgetrsBatched(handle, 0, msize, 1, djac_array, msize, dindx, drhs_array, msize, hinfo, nzbatchmx)

  ! Copy the solution back to the CPU
  istat = cublasGetVectorAsync(msize*nzbatchmx, sizeof_double, drhs, 1, hrhs, 1, stream)
  istat = cudaStreamSynchronize(stream)
  Do izb = 1, nzbatchmx
    If ( mask(izb) ) Then
      dy(:,izb) = rhs(1:ny,izb)
      If ( iheat > 0 ) dt9(izb) = rhs(ny+1,izb)
    EndIf
  EndDo

  If ( idiag >= 4 ) Then
    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        izone = izb + szbatch - 1
        Write(lun_diag,"(a,i5)") 'JAC_SOLVE',izone
        Write(lun_diag,"(14es10.3)") (dy(i,izb),i=1,ny)
        If ( iheat > 0 ) Write(lun_diag,"(es10.3)") dt9(izb)
      EndIf
    EndDo
  EndIf

  stop_timer = xnet_wtime()
  timer_solve = timer_solve + stop_timer

  Return
End Subroutine jacobian_solve

Subroutine jacobian_decomp(kstep,mask_in)
  !-------------------------------------------------------------------------------------------------
  ! This routine performs the LU matrix decomposition for the Jacobian.
  !-------------------------------------------------------------------------------------------------
  Use controls, Only: idiag, lun_diag, nzbatchmx, szbatch, lzactive
  Use cublasf, Only: cublasDgetrfBatched, cublasSetMatrixAsync
  Use gpu_controls, Only: handle
  Use jacobian_data
  Use timers, Only: xnet_wtime, start_timer, stop_timer, timer_solve
  Use xnet_types, Only: dp
  Implicit None

  ! Input variables
  Integer, Intent(in) :: kstep

  ! Optional variables
  Logical, Optional, Target, Intent(in) :: mask_in(:)

  ! Local variables
  Integer :: i, j, izb, izone, istat
  Logical, Pointer :: mask(:)

  If ( present(mask_in) ) Then
    mask => mask_in
  Else
    mask => lzactive
  EndIf

  start_timer = xnet_wtime()
  timer_solve = timer_solve - start_timer

  ! Calculate the LU decomposition
  istat = cublasDgetrfBatched(handle, msize, djac_array, msize, dindx, dinfo, nzbatchmx)

  If ( idiag >= 4 ) Then
    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        izone = izb + szbatch - 1
        Write(lun_diag,"(a3,i5)") 'LUD',izone
        Write(lun_diag,"(14es9.1)") ((jac(i,j,izb),j=1,msize),i=1,msize)
      EndIf
    EndDo
  EndIf

  stop_timer = xnet_wtime()
  timer_solve = timer_solve + stop_timer

  Return
End Subroutine jacobian_decomp

Subroutine jacobian_bksub(kstep,yrhs,dy,t9rhs,dt9,mask_in)
  !-------------------------------------------------------------------------------------------------
  ! This routine performs back-substitution for a LU matrix and the RHS vector.
  !-------------------------------------------------------------------------------------------------
  Use controls, Only: idiag, iheat, lun_diag, nzbatchmx, szbatch, lzactive
  Use cublasf, Only: cublasDgetrsBatched, cublasGetVectorAsync, cublasSetVectorAsync
  Use cudaf, Only: cudaStreamSynchronize
  Use gpu_controls, Only: handle, stream
  Use jacobian_data
  Use nuc_number, Only: ny
  Use timers, Only: xnet_wtime, start_timer, stop_timer, timer_solve
  Use xnet_types, Only: dp
  Implicit None

  ! Input variables
  Integer, Intent(in) :: kstep
  Real(dp), Intent(in) :: yrhs(:,:)
  Real(dp), Intent(in) :: t9rhs(:)

  ! Output variables
  Real(dp), Intent(out) :: dy(size(yrhs,1),size(yrhs,2))
  Real(dp), Intent(out) :: dt9(size(t9rhs))

  ! Optional variables
  Logical, Optional, Target, Intent(in) :: mask_in(:)

  ! Local variables
  Integer :: i, izb, izone, istat
  Logical, Pointer :: mask(:)

  If ( present(mask_in) ) Then
    mask => mask_in
  Else
    mask => lzactive
  EndIf

  start_timer = xnet_wtime()
  timer_solve = timer_solve - start_timer

  Do izb = 1, nzbatchmx
    If ( mask(izb) ) Then
      rhs(1:ny,izb) = yrhs(:,izb)
      If ( iheat > 0 ) rhs(ny+1,izb) = t9rhs(izb)
    EndIf
  EndDo

  ! Copy the RHS to the GPU
  istat = cublasSetVectorAsync(msize*nzbatchmx, sizeof_double, hrhs, 1, drhs, 1, stream)

  ! Solve the LU-decomposed triangular system via back-substitution
  istat = cublasDgetrsBatched(handle, 0, msize, 1, djac_array, msize, dindx, drhs_array, msize, hinfo, nzbatchmx)

  ! Copy the solution back to the CPU
  istat = cublasGetVectorAsync(msize*nzbatchmx, sizeof_double, drhs, 1, hrhs, 1, stream)

  istat = cudaStreamSynchronize(stream)
  Do izb = 1, nzbatchmx
    If ( mask(izb) ) Then
      dy(:,izb) = rhs(1:ny,izb)
      If ( iheat > 0 ) dt9(izb) = rhs(ny+1,izb)
    EndIf
  EndDo

  If ( idiag >= 4 ) Then
    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        izone = izb + szbatch - 1
        Write(lun_diag,"(a,i5)") 'BKSUB', izone
        Write(lun_diag,"(14es10.3)") (dy(i,izb),i=1,ny)
        If ( iheat > 0 ) Write(lun_diag,"(es10.3)") dt9(izb)
      EndIf
    EndDo
  EndIf

  stop_timer = xnet_wtime()
  timer_solve = timer_solve + stop_timer

  Return
End Subroutine jacobian_bksub

Subroutine yderiv_and_jacobian_build(diag,mult,mask_in)
  Use abundances, Only: yt, ydot
  Use conditions, Only: cv, t9t, t9dot
  Use controls, Only: idiag, iheat, lun_diag, nzbatchmx, szbatch, lzactive
  Use cross_sect_data, Only: csect1, csect2, csect3, n1i, n2i, n3i, &
    & dcsect1dt9, dcsect2dt9, dcsect3dt9
  Use jacobian_data, Only: jac, msize
  Use nuc_number, Only: ny
  Use nuclear_data, Only: nname, mex
  Use reac_rate_data, Only: a1, a2, a3, b1, b2, b3, la, le, mu1, mu2, mu3, &
    & n11, n21, n22, n31, n32, n33
  Use timers, Only: xnet_wtime, start_timer, stop_timer, timer_deriv, timer_jacob
  Use xnet_types, Only: dp
  Implicit None

  ! Input variables
  Real(dp), Intent(in) :: diag(:), mult

  ! Optional variables
  Logical, Optional, Target, Intent(in) :: mask_in(:)

  ! Local variables
  Integer :: i0, i1, i11, i21, i22, i31, i32, i33, j1, la1, le1, la2, le2, la3, le3, izb, izone
  Real(dp) :: dydotdt9(ny), dt9dotdy(msize), s1, s2, s3, r1, r2, r3, y11, y21, y22, y31, y32, y33
  Real(dp) :: dr1dt9, dr2dt9, dr3dt9
  Logical, Pointer :: mask(:)

  If ( present(mask_in) ) Then
    mask => mask_in
  Else
    mask => lzactive
  EndIf

  start_timer = xnet_wtime()
  timer_deriv = timer_deriv - start_timer
  timer_jacob = timer_jacob - start_timer

  !$acc kernels &
  !$acc copyin(diag,mask, &
  !$acc csect1,csect2,csect3,dcsect1dt9,dcsect2dt9,dcsect3dt9, &
  !$acc yt,ydot,t9dot,cv) &
  !$acc present(jac, &
  !$acc mex,la,le,mu1,mu2,mu3,n11,n21,n22,n31,n32,n33,a1,a2,a3)
  
  !$acc loop independent
  Do izb = 1, nzbatchmx
    If ( mask(izb) ) Then

      ! Calculate Ydot for each nucleus, summing over the reactions which affect it.
      !$acc loop independent
      Do i0 = 1, ny
        jac(i0,:,izb) = 0.0
        la1 = la(1,i0) ; la2 = la(2,i0) ; la3 = la(3,i0)
        le1 = le(1,i0) ; le2 = le(2,i0) ; le3 = le(3,i0)

        ! Sum over the reactions with 1 reactant
        s1 = 0.0
        Do j1 = la1, le1
          r1 = a1(j1)*csect1(mu1(j1),izb)
          i11 = n11(j1)
          y11 = yt(i11,izb)
          s1 = s1 + r1 * y11
          jac(i0,i11,izb) = jac(i0,i11,izb) + r1
        EndDo

        ! Sum over the reactions with 2 reactants
        s2 = 0.0
        Do j1 = la2, le2
          r2 = a2(j1)*csect2(mu2(j1),izb)
          i21 = n21(j1)
          i22 = n22(j1)
          y21 = yt(i21,izb)
          y22 = yt(i22,izb)
          s2 = s2 + r2 * y21 * y22
          jac(i0,i21,izb) = jac(i0,i21,izb) + r2 * y22
          jac(i0,i22,izb) = jac(i0,i22,izb) + r2 * y21
        EndDo

        ! Sum over the reactions with 3 reactants
        s3 = 0.0
        Do j1 = la3, le3
          r3 = a3(j1)*csect3(mu3(j1),izb)
          i31 = n31(j1)
          i32 = n32(j1)
          i33 = n33(j1)
          y31 = yt(i31,izb)
          y32 = yt(i32,izb)
          y33 = yt(i33,izb)
          s3 = s3 + r3 * y31 * y32 * y33
          jac(i0,i31,izb) = jac(i0,i31,izb) + r3 * y32 * y33
          jac(i0,i32,izb) = jac(i0,i32,izb) + r3 * y33 * y31
          jac(i0,i33,izb) = jac(i0,i33,izb) + r3 * y31 * y32
        EndDo

        ! Sum the 3 components of Ydot
        ydot(i0,izb) = s1 + s2 + s3
      EndDo

      If ( iheat > 0 ) Then

        ! Surprisingly, this seems to perform better than the DGEMV below
        s1 = 0.0
        Do i0 = 1, ny
          s1 = s1 + mex(i0)*ydot(i0,izb)
        EndDo
        t9dot(izb) = -s1 / cv(izb)

        !$acc loop independent
        Do i0 = 1, ny
          la1 = la(1,i0) ; la2 = la(2,i0) ; la3 = la(3,i0)
          le1 = le(1,i0) ; le2 = le(2,i0) ; le3 = le(3,i0)
          s1 = 0.0
          Do j1 = la1, le1
            dr1dt9 = a1(j1)*dcsect1dt9(mu1(j1),izb)
            i11 = n11(j1)
            y11 = yt(i11,izb)
            s1 = s1 + dr1dt9 * y11
          EndDo
          s2 = 0.0
          Do j1 = la2, le2
            dr2dt9 = a2(j1)*dcsect2dt9(mu2(j1),izb)
            i21 = n21(j1)
            i22 = n22(j1)
            y21 = yt(i21,izb)
            y22 = yt(i22,izb)
            s2 = s2 + dr2dt9 * y21 * y22
          EndDo
          s3 = 0.0
          Do j1 = la3, le3
            dr3dt9 = a3(j1)*dcsect3dt9(mu3(j1),izb)
            i31 = n31(j1)
            i32 = n32(j1)
            i33 = n33(j1)
            y31 = yt(i31,izb)
            y32 = yt(i32,izb)
            y33 = yt(i33,izb)
            s3 = s3 + dr3dt9 * y31 * y32 * y33
          EndDo
          jac(i0,ny+1,izb)= s1 + s2 + s3
        EndDo

        !$acc loop independent
        Do i1 = 1, msize
          s1 = 0.0
          Do i0 = 1, ny
            s1 = s1 + mex(i0)*jac(i0,i1,izb)
          EndDo
          jac(ny+1,i1,izb) = -s1 / cv(izb)
        EndDo
        ! The BLAS version of jac(ny+1,i,izb) = -sum(mex*jac(1:ny,i,izb))/cv(izb)
!       Call dgemv('T',ny,msize,1.0/cv(izb),jac(1,1,izb),msize,mex,1,0.0,dt9dotdy,1)
!       jac(ny+1,:,izb) = -dt9dotdy
      EndIf
    EndIf
  EndDo

  ! This also works, but could be inefficient if there are few active zones in batch
! If ( iheat > 0 ) Then
!   Call dgemv('T',ny,nzbatchmx,1.0,ydot,ny,mex,1,0.0,t9dot,1)
!   Call dgemv('T',ny,msize*nzbatchmx,1.0,jac(1,1,1),msize,mex,1,0.0,dt9dotdy,1)
!   Do izb = 1, nzbatchmx
!     If ( mask(izb) ) Then
!       t9dot(izb) = -t9dot(izb)/ cv(izb)
!       Do i1 = 1, msize
!         jac(ny+1,i1,izb) = -dt9dotdy(i1,izb) / cv(izb)
!       EndDo
!     EndIf
!   EndDo
! EndIf

  ! Apply the externally provided factors
  !$acc loop independent
  Do izb = 1, nzbatchmx
    If ( mask(izb) ) Then
      jac(:,:,izb) = mult * jac(:,:,izb)
      !$acc loop independent
      Do i0 = 1, msize
        jac(i0,i0,izb) = jac(i0,i0,izb) + diag(izb)
      EndDo
    EndIf
  EndDo
  !$acc end kernels

  stop_timer = xnet_wtime()
  timer_deriv = timer_deriv + stop_timer
  timer_jacob = timer_jacob + stop_timer

  Return
End Subroutine yderiv_and_jacobian_build
