!***************************************************************************************************
! data_distribute_mpi.f90 10/18/17
! Needed for MPI execution. These routines broadcast the nuclear and network data between MPI tasks.
!***************************************************************************************************

Subroutine control_bcast(data_dir)
  !-------------------------------------------------------------------------------------------------
  ! This routine reads and distrbutes the control file which contains the parameters which control
  ! the action of the network. control also indicates the relative directory from which the nuclear
  ! data should be loaded, as well as the names of the files containing the initial abundances and
  ! thermodynamic trajectories.
  !-------------------------------------------------------------------------------------------------
  Use controls
  Use xnet_interface, Only: net_preprocess
  Use xnet_types, Only: dp
  Use mpi
  Implicit None

  ! Output variables
  Character(80), Intent(out) :: data_dir

  ! Local variables
  Character(80) :: control_char(6)
  Real(dp) :: control_real(9)
  Integer :: control_int(13)
  Integer :: ierr

  ! PE0 ...
  If ( myid == 0 ) Then

    Call read_controls(data_dir)
    IF ( iprocess > 0 ) CALL net_preprocess( lun_stdout, data_dir, data_dir )

    ! Load Control passing arrays
    control_int(1)  = szone     ; control_int(2)  = nzone  ; control_int(3)  = iweak0
    control_int(4)  = iscrn     ; control_int(5)  = isolv  ; control_int(6)  = kstmx
    control_int(7)  = kitmx     ; control_int(8)  = ijac   ; control_int(9)  = iconvc
    control_int(10) = idiag     ; control_int(11) = itsout ; control_int(12) = iheat
    control_int(13) = ineutrino

    control_real(1) = changemx ; control_real(2) = yacc      ; control_real(3) = tolm
    control_real(4) = tolc     ; control_real(5) = ymin      ; control_real(6) = tdel_maxmult
    control_real(7) = t9nse    ; control_real(8) = changemxt ; control_real(9) = tolt9

    control_char(1:3) = descript       ; control_char(4) = data_dir
    control_char(5)   = bin_file_base  ; control_char(6) = ev_file_base
  EndIf

  ! All PE
  ! Broadcast network control parameters
  Call mpi_bcast(control_int,13,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(control_real,9,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(control_char,6*80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(output_nuc,14*5,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

  ! Unpack network control parameters
  If ( myid /= 0 ) Then
    szone     = control_int(1)  ; nzone  = control_int(2)  ; iweak0 = control_int(3)
    iscrn     = control_int(4)  ; isolv  = control_int(5)  ; kstmx  = control_int(6)
    kitmx     = control_int(7)  ; ijac   = control_int(8)  ; iconvc = control_int(9)
    idiag     = control_int(10) ; itsout = control_int(11) ; iheat  = control_int(12)
    ineutrino = control_int(13)

    changemx = control_real(1) ; yacc      = control_real(2) ; tolm         = control_real(3)
    tolc     = control_real(4) ; ymin      = control_real(5) ; tdel_maxmult = control_real(6)
    t9nse    = control_real(7) ; changemxt = control_real(8) ; tolt9        = control_real(9)

    descript       = control_char(1:3) ; data_dir         = control_char(4)
    bin_file_base  = control_char(5)   ; ev_file_base     = control_char(6)
    Allocate (inab_file(nzone),thermo_file(nzone))
  EndIf

  ! Broadcast input datafiles
  Call mpi_bcast(inab_file,nzone*80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(thermo_file,nzone*80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

  !$omp parallel default(shared)
  If ( idiag > 0 ) Write(lun_diag,*) 'CBcast',idiag,kitmx,kstmx
  !$omp end parallel

  Return
End Subroutine control_bcast

Subroutine netdata_bcast(data_dir,data_desc)
  !-------------------------------------------------------------------------------------------------
  ! This routine handles nuclear data I/O for MPI versions, broadcasting the necessary nuclear and
  ! reaction data from PE0 to the production PEs.
  !-------------------------------------------------------------------------------------------------
  Use controls, Only: idiag, iheat, iscrn, lun_diag, myid, nzbatchmx
  Use nuc_number, Only: ny
  Use nuclear_data
  Use part_funct_data
  Use cross_sect_data
  Use ffn_data
! Use nnu_data                                                                                  !NNU
! Use neutrino_data                                                                             !NNU
  Use reac_rate_data
  Use screening_data
  Use xnet_interface, Only: read_nuclear_data, read_reaction_data
  Use mpi
  Implicit None

  ! Input variables
  Character(*), Intent(in) :: data_dir

  ! Output variables
  Character(80), Intent(out) :: data_desc

  ! Local variables
  Integer :: matrix_shape(5)
  Integer :: nbc, ierr, nr1, nr2, nr3

  ! On PE0 ...
  If ( myid == 0 ) Then

    ! ... read the nuclear and reaction data
    Call read_nuclear_data(data_dir,data_desc)
    Call read_reaction_data(data_dir)

  EndIf

  ! Share data for nuc_number module
  Call mpi_bcast(ny,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  ! Share data description
  Call mpi_bcast(data_desc,80,MPI_CHARACTER, 0,MPI_COMM_WORLD,ierr)

  ! Share data for the nuclear_data module
  If ( myid /= 0 ) Allocate(nname(0:ny),aa(ny),zz(ny),nn(ny),be(ny),mex(ny),mm(ny),zz2(ny),zz53(ny),zzi(ny))
  Call mpi_bcast(aa,ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(zz,ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(nn,ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(be,ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(mm,ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(mex,ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(mex_n,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(mex_p,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(ineut,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(iprot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(inmin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(inmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(izmin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(izmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(zz2,ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(zz53,ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(zzi,ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  nbc = 5*(ny+1)
  Call mpi_bcast(nname,nbc,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

  ! Share data for the part_funct_data module
  If ( myid /= 0 ) Allocate(g(ng,ny),angm(0:ny))
  Call mpi_bcast(t9i,ng,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(g,ng*ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(angm,ny+1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

  ! Share data for the cross_sect_data module
  Call mpi_bcast(nreac,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  nr1 = nreac(1) ; nr2 = nreac(2) ; nr3 = nreac(3)
  If ( myid /= 0 ) Then
    Allocate (iffn(nr1),innu(nr1))
    Allocate (n1i(4,nr1),iwk1(nr1),ires1(nr1),irev1(nr1),rc1(7,nr1),q1(nr1))
    Allocate (n2i(5,nr2),iwk2(nr2),ires2(nr2),irev2(nr2),rc2(7,nr2),q2(nr2))
    Allocate (n3i(6,nr3),iwk3(nr3),ires3(nr3),irev3(nr3),rc3(7,nr3),q3(nr3))
  EndIf
  nbc = nreac(1)
  Call mpi_bcast(rc1,7*nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(n1i,4*nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(iwk1,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(ires1,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(irev1,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(q1,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(iffn,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(innu,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  nbc = nreac(2)
  Call mpi_bcast(rc2,7*nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(n2i,5*nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(iwk2,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(ires2,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(irev2,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(q2,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  nbc = nreac(3)
  Call mpi_bcast(rc3,7*nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(n3i,6*nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(iwk3,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(ires3,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(irev3,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(q3,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

  ! Share the data for the reac_rate_data module
  Call mpi_bcast(nan,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  If ( myid /= 0 ) Then
    Allocate (la(3,ny),le(3,ny))
    Allocate (mu1(nan(1)),a1(nan(1)),n11(nan(1)))
    Allocate (mu2(nan(2)),a2(nan(2)),n21(nan(2)),n22(nan(2)))
    Allocate (mu3(nan(3)),a3(nan(3)),n31(nan(3)),n32(nan(3)),n33(nan(3)))
  Endif
  nbc = 3*ny
  Call mpi_bcast(la,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(le,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  nbc = nan(1)
  Call mpi_bcast(mu1,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(a1,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(n11,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  nbc = nan(2)
  Call mpi_bcast(mu2,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(a2,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(n21,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(n22,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  nbc = nan(3)
  Call mpi_bcast(mu3,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(a3,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(n31,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(n32,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(n33,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  ! Share the data for the screening_data module
  If ( iscrn > 0 ) Then
    If ( myid /= 0 ) Then
      Allocate (zseq(0:izmax+2),zseq53(0:izmax+2),zseqi(0:izmax+2))
      Allocate (z21(nr2),z22(nr2),iz21(nr2),iz22(nr2),iz2c(nr2),z2c(nr2),zeta2w(nr2),zeta2i(nr2))
      Allocate (z31(nr3),z32(nr3),z33(nr3),iz31(nr3),iz32(nr3),iz33(nr3),iz3c(nr3),z3c(nr3),zeta3w(nr3),zeta3i(nr3))
    EndIf
    nbc = izmax+3
    Call mpi_bcast(zseq,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(zseq53,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(zseqi,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    nbc = nreac(2)
    Call mpi_bcast(z21,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(z22,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(iz21,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(iz22,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(iz2c,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(z2c,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(zeta2w,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(zeta2i,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    nbc = nreac(3)
    Call mpi_bcast(z31,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(z32,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(z33,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(iz31,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(iz32,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(iz33,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(iz3c,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(z3c,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(zeta3w,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(zeta3i,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  EndIf

  ! Share the data for the nnu_data
! Call mpi_bcast(nnnu,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)                                      !NNU
! If ( nnnu > 0 ) Then                                                                          !NNU
!   If ( myid /= 0 ) Then                                                                       !NNU
!     Allocate (sigmanu(nnnu,7))                                                                !NNU
!     !$omp parallel default(shared)                                                            !NNU
!     Allocate (rcsnu(nnnu,4))                                                                  !NNU
!     !$omp end parallel                                                                        !NNU
!   EndIf                                                                                       !NNU
!   Call mpi_bcast(sigmanu,7*nnnu,MPI_REAL8,0,MPI_COMM_WORLD,ierr)                              !NNU
! ElseIf ( myid /= 0 ) Then                                                                     !NNU
!   Allocate (sigmanu(1,7))                                                                     !NNU
!   Allocate (rcsnu(1,7))                                                                       !NNU
! EndIf                                                                                         !NNU

  ! Share the data for the ffn_data
  Call mpi_bcast(nffn,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  If ( nffn > 0) Then
    If ( myid /= 0 ) Then
      Allocate (ffnsum(nffn,ngrid),ffnenu(nffn,ngrid))
    Endif
    nbc = nffn*ngrid
    Call mpi_bcast(ffnsum,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    Call mpi_bcast(ffnenu,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  ElseIf ( myid /= 0 ) Then
    Allocate(ffnsum(1,ngrid),ffnenu(1,ngrid))
  Endif

  ! Allocate threadprivate arrays
  !$omp parallel default(shared)
  If ( myid /= 0 ) Then
    Allocate (gg(0:ny,nzbatchmx))
    Allocate (csect1(nr1,nzbatchmx),h1(nr1,nzbatchmx),rpf1(nr1))
    Allocate (csect2(nr2,nzbatchmx),h2(nr2,nzbatchmx),rpf2(nr2))
    Allocate (csect3(nr3,nzbatchmx),h3(nr3,nzbatchmx))
    Allocate (b1(nan(1),nzbatchmx))
    Allocate (b2(nan(2),nzbatchmx))
    Allocate (b3(nan(3),nzbatchmx))
    If ( iheat > 0 ) Then
      Allocate (dlngdt9(0:ny,nzbatchmx))
      Allocate (dcsect1dt9(nr1,nzbatchmx),dh1dt9(nr1,nzbatchmx),dlnrpf1dt9(nr1))
      Allocate (dcsect2dt9(nr2,nzbatchmx),dh2dt9(nr2,nzbatchmx),dlnrpf2dt9(nr2))
      Allocate (dcsect3dt9(nr3,nzbatchmx),dh3dt9(nr3,nzbatchmx))
    EndIf
  EndIf
  If ( idiag > 0 ) Write(lun_diag,*) 'NBcast',ny,nffn,nreac
  !$omp end parallel

  Return
End Subroutine netdata_bcast

Subroutine match_bcast(data_dir)
  !-------------------------------------------------------------------------------------------------
  ! This routine reads in the match_data modules on PE0 and broadcasts the data to other processors.
  !-------------------------------------------------------------------------------------------------
  Use controls, Only: idiag, lun_diag, myid
  Use cross_sect_data, Only: nreac
  Use match_data, Only: mflx, nflx, ifl1, ifl2, ifl3, iwflx, qflx, descx
  Use xnet_interface, Only: read_match_data
  Use mpi
  Implicit None

  ! Input variables
  Character(*), Intent(in) :: data_dir

  ! Local variables
  Integer :: ierr

  ! The control PE...
  If ( myid == 0 ) Then

    ! ... reads in the flux matching data
    Call read_match_data(data_dir)

  EndIf

  ! Share the match data
  Call mpi_bcast(mflx,1,       MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  If ( myid /= 0 ) Then
      Allocate (ifl1(nreac(1)),ifl2(nreac(2)),ifl3(nreac(3)))
      Allocate (nflx(7,mflx),qflx(mflx),iwflx(mflx),descx(mflx))
  EndIf
  Call mpi_bcast(ifl1,nreac(1),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(ifl2,nreac(2),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(ifl3,nreac(3),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(nflx,7*mflx,  MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(qflx,mflx,    MPI_REAL8,  0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(iwflx,mflx,   MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  Call mpi_bcast(descx,4*mflx, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

  !$omp parallel default(shared)
  If ( idiag > 0 ) Write(lun_diag,*) 'MBcast',mflx
  !$omp end parallel

  Return
End Subroutine match_bcast
