module cleaning_module

contains

  function sanitize(val) result (new_val)

    implicit none

    double precision, intent(in) :: val
    double precision :: new_val

    if (abs(val) < 1.d-99) then
       new_val = 0.0d0
    else
       new_val = val
    endif

  end function sanitize

end module cleaning_module

program rewrite

  ! the NSE table has some very small numbers of the form 5.045-100 --
  ! this is how Fortran writes numbers with exponents < -99 -- it
  ! drops the "E".  Here we reset those small values to 0 so the table
  ! can be read in C++

  use cleaning_module

  implicit none

  integer, parameter :: ntemp = 71
  integer, parameter :: nden = 31
  integer, parameter :: nye = 21

  integer, parameter :: nspec = 19

  integer :: uold, unew
  double precision :: ttlog, ddlog, yetab
  double precision :: the, tsi, tfe
  double precision :: abartab, ebtab, wratetab
  double precision :: massfractab(nspec)

  integer :: irho, it9, iye
  integer :: k

  ! read in table
  open(newunit=uold, file="nse19.tbl.orig")

  open(newunit=unew, file="nse19.tbl")

5 format(2f12.5, 1pe12.5, 6e12.5, 19e12.5)

  do irho = 1, nden
     do it9 = 1, ntemp
        do iye = 1, nye

           ! read in the existing table
           read (uold, 5) ttlog, ddlog, yetab, &
                the, tsi, tfe, &
                abartab, ebtab, wratetab, &
                (massfractab(k), k=1, nspec)

           ! sanitize and output
           write (unew, 5) sanitize(ttlog), sanitize(ddlog), sanitize(yetab), &
                sanitize(the), sanitize(tsi), sanitize(tfe), &
                sanitize(abartab), sanitize(ebtab), sanitize(wratetab), &
                (sanitize(massfractab(k)), k=1, nspec)

        end do
     end do
  end do

end program rewrite
