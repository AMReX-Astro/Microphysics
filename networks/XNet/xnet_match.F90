!***************************************************************************************************
! match.f90 10/18/17
! This file containes the data structures and input routines for matching reactions which involve
! the same nuclei (forward and reverse reactions as well as reactions with multiple components).
!***************************************************************************************************

Module match_data
  !-------------------------------------------------------------------------------------------------
  ! This module contains the data necessary to match up reactions.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None
  Integer                   :: mflx                      ! Number of unique reaction pairs
  Integer, Allocatable      :: nflx(:,:)                 ! The nuclei in each unique reaction
  Integer, Allocatable      :: iwflx(:)                  ! Reaction pair weak flag
  Integer, Allocatable      :: ifl1(:), ifl2(:), ifl3(:) ! Maps reaction rate arrays to reaction pair arrays
  Real(dp), Allocatable     :: qflx(:)                   ! Reaction pair Q values
  Character(4), Allocatable :: descx(:)                  ! Descriptor for the reaction pair
End Module match_data

Subroutine read_match_data(data_dir)
  !-------------------------------------------------------------------------------------------------
  ! This routine reads in the reaction matching data and allocates the necessary arrays.
  !-------------------------------------------------------------------------------------------------
  Use controls, Only: idiag, lun_diag, lun_stdout, getNewUnit
  Use cross_sect_data, Only: nreac
  Use match_data, Only: mflx, nflx, ifl1, ifl2, ifl3, iwflx, qflx, descx
  Implicit None

  ! Input variables
  Character(*), Intent(in) :: data_dir

  ! Local variables
  Integer :: i, lun_match, nr(3)

  ! Open and read the matching data arrays
  Open(getNewUnit(lun_match), file=trim(data_dir)//"/match_data", form="unformatted", status="old")
  Read(lun_match) mflx, nr
  Allocate (ifl1(nr(1)),ifl2(nr(2)),ifl3(nr(3)))
  Allocate (nflx(7,mflx),qflx(mflx),iwflx(mflx),descx(mflx))
  Read(lun_match) ifl1, ifl2, ifl3
  Read(lun_match) nflx,qflx,iwflx,descx
  Close(lun_match)

  ! Make sure match_data agrees with cross_sect_data
  Do i = 1, 3
    If ( nr(i) /= nreac(i) ) Write(lun_stdout,*) 'NR mismatch',i,nr(i),nreac(i)
  EndDo

  !$omp parallel default(shared)
  If ( idiag > 0 ) Write(lun_diag,*) 'Match',mflx
  !$omp end parallel

  Return
End Subroutine read_match_data
