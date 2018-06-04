!***************************************************************************************************
! util.f90 10/18/17
! This file contains various utility routines that are common throughout XNet.
!***************************************************************************************************

Subroutine benuc(y,enb,enm,ytot,ztot,atot)
  !-------------------------------------------------------------------------------------------------
  ! This routine finds moments of the abundance distribution useful for hydrodynamics, including the
  ! total abundance, electron fraction, binding energy, and mass excess energy and outputs in mol/g
  ! and ergs/g.
  !-------------------------------------------------------------------------------------------------
  Use xnet_constants, Only: epmev, avn
  Use nuclear_data, Only: zz, nn, aa, be, mex_p, mex_n
  Use xnet_types, Only: dp
  Implicit None

  ! Input variables
  Real(dp), Intent(in) :: y(:)

  ! Output variables
  Real(dp), Intent(out) :: enb ! Binding energy [ergs g^{-1}]
  Real(dp), Intent(out) :: enm ! Mass excess [ergs g^{-1}]
  Real(dp), Intent(out) :: ytot, ztot, atot ! Abundance moments

  ytot = sum(y)
  ztot = sum(y * zz)
  atot = sum(y * aa)
  enb  = sum(y * be)
  enm  = mex_p*ztot + mex_n*(1.0-ztot) - enb

  ! Change units from MeV/nucleon to erg/g
  enb = epmev * avn * enb
  enm = epmev * avn * enm

  Return
End Subroutine benuc

Subroutine norm(yy)
  !-------------------------------------------------------------------------------------------------
  ! This routine renormalizes the abundances to guarantee mass conservation.
  !-------------------------------------------------------------------------------------------------
  Use nuclear_data, Only: aa
  Use xnet_types, Only: dp
  Implicit None

  ! Input/Output variables
  Real(dp), Intent(inout) :: yy(:)

  ! Local variables
  Real(dp) :: xtot, rxt

  ! Renormalize total mass fraction to 1
  xtot = sum(yy * aa)
  rxt  = 1.0 / xtot
  yy   = yy * rxt

  Return
End Subroutine norm

Subroutine ye_norm(yy,ye)
  !-------------------------------------------------------------------------------------------------
  ! This routine renormalizes the abundances to guarantee mass and charge conservation.
  !-------------------------------------------------------------------------------------------------
  Use nuclear_data, Only: zz, nn, aa
  Use xnet_types, Only: dp
  Implicit None

  ! Input variables
  Real(dp), Intent(in) :: ye

  ! Input/Output variables
  Real(dp), Intent(inout) :: yy(:)

  ! Local variables
  Real(dp) :: zy, nny, zny, zzy, beta, alph

  ! Calculate moments needed to calculate normalization coefficients
  nny  = sum(nn * yy)
  zy   = sum(zz * yy)
  zny  = sum(nn * zz * yy / aa)
  zzy  = sum(zz * zz * yy / aa)
  beta = (ye*nny - zny) / (nny*zzy - zy*zny)
  alph = (1.0 - beta*zy) / nny
  yy   = yy * (alph*nn + beta*zz) / aa

  Return
End Subroutine ye_norm

Subroutine index_from_name(nuc_name,inuc)
  !-------------------------------------------------------------------------------------------------
  ! This routine takes a nuclear name and finds the corresponding index for the current data set.
  ! inuc = 0      indicates a blank name.
  ! inuc = ny + 1 indicates that the nuclear name is not found in the current set.
  !-------------------------------------------------------------------------------------------------
  Use nuc_number, Only: ny
  Use nuclear_data, Only: nname
  Implicit None

  ! Input variables
  Character(*), Intent(in) :: nuc_name

  ! Output variables
  Integer, Intent(out) :: inuc

  ! Local variables
  Integer :: i, name_len

  name_len = len_trim(nuc_name)
  If ( name_len > 0 ) then
    Do i = 1, ny
      If ( trim(adjustl(nuc_name)) == trim(adjustl(nname(i))) ) Exit
    EndDo
    inuc = i
  Else
    inuc = 0
  EndIf

  Return
End Subroutine index_from_name

Subroutine name_ordered(base_string,n,nmax)
  !-------------------------------------------------------------------------------------------------
  ! This routine appends the integer n (padded with "0"s up to the size of the integer nmax) to the
  ! character variable base_string.
  !-------------------------------------------------------------------------------------------------
  Implicit None

  ! Input variables
  Integer, Intent(in) :: n, nmax

  ! Input/Output variables
  Character(*), Intent(inout) :: base_string

  ! Local variables
  Integer :: nmax_digits
  Character(1) :: nmax_digits_string
  Character(6) :: n_format
  Character(9) :: n_string

  ! Find character length of imax
  Write(nmax_digits_string,"(i1)") int(log10(real(nmax))) + 1

  ! Construct format spec and write n as zero padded string
  n_format = "(i"//nmax_digits_string//"."//nmax_digits_string//")"
  Write(n_string,n_format) n

  ! Append n_string to base_string
  base_string = trim(base_string)//trim(n_string)

  Return
End Subroutine name_ordered

Subroutine string_lc(string)
  !-------------------------------------------------------------------------------------------------
  ! This routine converts an ASCII string to all lower case.
  !-------------------------------------------------------------------------------------------------
  Implicit None

  ! Input/Output variables
  Character(*), Intent(inout) :: string

  ! Local variables
  Integer, Parameter :: lc_a_ascii=iachar('a')
  Integer, Parameter :: uc_a_ascii=iachar('A')
  Integer, Parameter :: lc_z_ascii=iachar('z')
  Integer, Parameter :: uc_z_ascii=iachar('Z')
  Integer :: i, x

  Do i = 1, len_trim(string)
    x = iachar(string(i:i))
    If ( x >= uc_a_ascii .and. x <= uc_z_ascii ) Then
      x = x + (lc_a_ascii - uc_a_ascii)
      string(i:i) = achar(x)
    EndIf
  EndDo

  Return
End Subroutine string_lc
