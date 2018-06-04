!***************************************************************************************************
! jacobian_bcast_dense.f90 10/18/17
! Needed for MPI execution with dense matrix solver.
! This routine broadcasts the jacobian data between MPI tasks.
!***************************************************************************************************

Subroutine jacobian_bcast(data_dir)
  !-------------------------------------------------------------------------------------------------
  ! This routine distributes Jacobian data for dense solver.
  !-------------------------------------------------------------------------------------------------
  Use xnet_interface, Only: read_jacobian_data
  Implicit None

  ! Input variables
  Character(*), Intent(in) :: data_dir

  ! For dense solvers, no additional data must be read or broadcast, only allocations performed
  Call read_jacobian_data(data_dir)

  Return
End Subroutine jacobian_bcast
