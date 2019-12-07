module util_module

  use amrex_constants_module
  use amrex_fort_module, only : rt => amrex_real

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine get_xn(xn_zone)

    use network,       only: nspec, spec_names
    use extern_probin_module, only: xin_file

    use amrex_fort_module, only : rt => amrex_real
    real(rt), intent(  out) :: xn_zone(:,:)

    integer         :: un, i
    real(rt) :: summ, usr_in
    character (len=4096) :: line

    ! read in an inputs file containing the mass fractions.
    ! each species is on its own line.
    ! Allow for comment lines with '#' in the first column
    open(newunit=un, file=xin_file, status='old')

    summ = ZERO

    i = 1
    do while (i <= nspec)
       ! read the line into a character buffer
       read (un,'(a)') line
       if (index(line, '#') == 1) cycle

       read (line,*) xn_zone(i,:)

       i = i + 1
    enddo

    ! need to add check on the species name, and check that we sum to 1

    close(un)
  end subroutine get_xn

end module util_module
