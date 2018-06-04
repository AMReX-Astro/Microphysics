!***************************************************************************************************
! gpu_control.f90 10/18/17
! This file contains modules and subroutines to control the GPU execution of XNet.
!***************************************************************************************************

Module gpu_controls
  Use, Intrinsic :: iso_c_binding, Only: C_INT, C_PTR
  Implicit None

  ! CUDA/CUBLAS management
  Type(C_PTR) :: handle, stream, event
  !$omp threadprivate(handle,stream,event)

  Integer(C_INT) :: deviceCount
  Integer :: mydevice

Contains

  Subroutine gpu_init
    Use controls, Only: lun_stdout, myid
    Use cublasf, Only: cublasCreate_v2, cublasSetStream_v2, CUBLAS_STATUS_SUCCESS
    Use cudaf, Only: cudaGetDeviceCount, cudaDeviceProp, cudaGetDeviceProperties, cudaSetDevice, &
      & cudaStreamCreateWithFlags, cudaStreamNonBlocking, cudaStreamDefault, cudaSuccess, &
      & cudaEventCreateWithFlags, cudaEventDefault, cudaDeviceSynchronize
    Use openacc, Only: acc_init, acc_set_device_num, acc_device_nvidia
    Implicit None

    ! Local variables
    Integer :: istat
    Type(cudaDeviceProp) :: deviceProp

    ! Initialize GPU
    istat = cudaGetDeviceCount(deviceCount)
    If ( istat /= cudaSuccess ) Write(lun_stdout,*) "cudaGetDeviceCount, istat", istat

    If ( deviceCount > 0 ) Then
      mydevice = mod(myid,deviceCount)
    Else
      Write(lun_stdout,*) 'No CUDA capable device found'
    EndIf

    istat = cudaGetDeviceProperties(deviceProp,mydevice)
    If ( istat /= cudaSuccess ) Write(lun_stdout,*) "cudaGetDeviceProperties, istat", istat

!    Write(lun_stdout,'(3(a,i2),3(a,i1))') "Rank: ",myid,", Device: ",mydevice+1," (of ",deviceCount, &
!      "), CC: ",deviceProp%major,".",deviceProp%minor,", ComputeMode: ",deviceProp%computeMode

    !$omp parallel default(shared) private(istat)

    istat = cudaSetDevice(mydevice)
    If ( istat /= cudaSuccess ) Write(lun_stdout,*) "cudaSetDevice, istat", istat
    
    ! Create cublas handles
    istat = cublasCreate_v2(handle)
    If ( istat /= CUBLAS_STATUS_SUCCESS ) Write(lun_stdout,*) 'cublasCreate_v2, istat', istat

    ! Create CUDA streams
!   istat = cudaStreamCreateWithFlags(stream, cudaStreamDefault)
    istat = cudaStreamCreateWithFlags(stream, cudaStreamNonBlocking)
    if (istat /= cudaSuccess) Write(lun_stdout,*) "cudaStreamCreateWithFlags, istat", istat

    ! Associate each cublas handle with a CUDA stream
    istat = cublasSetStream_v2(handle, stream)
    if ( istat /= CUBLAS_STATUS_SUCCESS ) Write(lun_stdout,*) 'cublasSetStream_v2, istat', istat

    istat = cudaEventCreateWithFlags(event, cudaEventDefault)

    call acc_set_device_num( mydevice, acc_device_nvidia )
    call acc_init( acc_device_nvidia )

    !$acc enter data

    !$omp end parallel

    istat = cudaDeviceSynchronize()
    If ( istat /= cudaSuccess ) Write(lun_stdout,*) "cudaDeviceSynchronize, istat", istat

    Return
  End Subroutine gpu_init

  Subroutine gpu_finalize
    Use cublasf, Only: cublasDestroy_v2
    Use cudaf, Only: cudaEventDestroy, cudaStreamDestroy
    Implicit None

    ! Local variables
    Integer :: istat

    !$omp parallel default(shared) private(istat)

    !$acc exit data

    istat = cudaEventDestroy(event)
    istat = cudaStreamDestroy(stream)
    istat = cublasDestroy_v2(handle)

    !$omp end parallel

    Return
  End Subroutine gpu_finalize

End Module gpu_controls
