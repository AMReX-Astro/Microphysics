program test

  use burning_module

  implicit none

  double precision :: rholog,rhologf,dummy,c12,drho
  double precision, parameter :: rhologstart = 7.0d0, &
                                rhologend = 10.0d0
  integer, parameter :: nsteps = 1000

  drho = (rhologend - rhologstart)/dble(nsteps+1)

  rholog = rhologstart
  c12 = 0.5d0
  do while (rholog <= rhologend)
     call networkburn(rholog,c12,rhologf,dummy,dummy,dummy,dummy,dummy, &
          dummy,dummy,dummy,dummy,dummy)
     print *, rholog, rhologf
     rholog = rholog + drho
  end do
  
end program test


