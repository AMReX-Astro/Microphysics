!***************************************************************************************************
! terminate_mpi.f90 10/18/17
! This file contains routines for gracefully stopping XNet exeuction in the event of an error.
!***************************************************************************************************

Subroutine xnet_terminate(c_diagnostic,i_diagnostic)
  !-------------------------------------------------------------------------------------------------
  ! This routine gracefully exits XNet with a diagnostic statement
  !-------------------------------------------------------------------------------------------------
  Use Driver_interface, Only: Driver_abortFlash
  Implicit None

  ! Input variables
  Character(*), Intent(in) :: c_diagnostic
  Integer, Intent(in), Optional :: i_diagnostic

  Call Driver_abortFlash(c_diagnostic)

  Return
End Subroutine xnet_terminate