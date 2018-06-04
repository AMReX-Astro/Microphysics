!***************************************************************************************************
! types.f90 10/18/17
! This file contains the module defining different fortran types.
!***************************************************************************************************

Module xnet_types
  !-------------------------------------------------------------------------------------------------
  ! These are explicitly defined type parameters for floating-point (real) number precision, with
  ! double precision declared as Real(dp) instead of Real(8).
  !-------------------------------------------------------------------------------------------------
  Use, Intrinsic :: iso_fortran_env, Only: int64, real32, real64, real128
  Implicit None
  Integer, Parameter :: i8 = int64   ! 64-bit integer
  Integer, Parameter :: sp = real32  ! Single precision
  Integer, Parameter :: dp = real64  ! Double precision
  Integer, Parameter :: qp = real128 ! Quad percision
End Module xnet_types