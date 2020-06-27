program usetable

  use nse_module

  implicit none

  integer :: irho, it9, iye
  integer :: j, k

  integer :: ichoice
  double precision :: t9, rho, ye, X(nspec)
  double precision ::xlogt9, xlogrho
  double precision :: abar, dq, dyedt


  ! enter test points

  call init_nse()

  do while (.true.)

     print *, "enter 1 for log10 entry; 2 for actual values; anything else to exit"

     read (*,*) ichoice

     if (ichoice == 1) then

        print *, "Enter log10(T), log10(rho), ye"
        read (*,*) xlogt9, xlogrho, ye

        t9 = 10.0d0**(xlogt9-9.0d0)
        rho = 10.0d0**xlogrho

     else if (ichoice == 2) then

        print *, "Enter T9, rho, ye"
        read (*,*) t9, rho, ye

     else

        call exit(1)

     end if

     call interp (t9, rho, ye, abar, dq, dyedt, X)

     print *, "      t9    ", "        rho    ", "       ye     ", &
              "     abar     ", "      be/a    ", "      dyedt  "

     write (*,41) t9, rho, ye, abar, dq, dyedt, X(:)
41   format (1pe12.3, 5e14.5, 19e14.5)

  end do

end program usetable
