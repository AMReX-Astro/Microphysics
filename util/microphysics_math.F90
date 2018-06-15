
module microphysics_math_module
  
  use amrex_error_module

  implicit none

  public

contains

#ifdef CUDA
  attributes(device) &
#endif
  function esum10(array) result(esum)

    !$acc routine seq

    use bl_error_module, only: bl_error
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: array(:)
    real(rt) :: esum

    integer :: i, j, k, km

    ! Note that for performance reasons we are not
    ! initializing the unused values in this array.

    real(rt) :: partials(0:9)
    real(rt) :: x, y, z, hi, lo

    ! j keeps track of how many entries in partials are actually used.
    ! The algorithm we model this off of, written in Python, simply
    ! deletes array entries at the end of every outer loop iteration.
    ! The Fortran equivalent to this might be to just zero them out,
    ! but this results in a huge performance hit given how often
    ! this routine is called during in a burn. So we opt instead to
    ! just track how many of the values are meaningful, which j does
    ! automatically, and ignore any data in the remaining slots.

    j = 0

    ! The first partial is just the first term.
    partials(j) = array(1)

    do i = 2, 10

       km = j
       j = 0

       x = array(i)

       do k = 0, km
          y = partials(k)

          if (abs(x) < abs(y)) then
             ! Swap x, y
             z = y
             y = x
             x = z
          endif

          hi = x + y
          lo = y - (hi - x)

          if (lo .ne. 0.0_rt) then
             partials(j) = lo
             j = j + 1
          endif

          x = hi

       enddo

       partials(j) = x

    enddo

    esum = sum(partials(0:j))

  end function esum10


#ifdef CUDA
  attributes(device) &
#endif
  function esum12(array) result(esum)

    !$acc routine seq

    use amrex_error_module, only: amrex_error
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    real(rt), intent(in) :: array(:)
    real(rt) :: esum

    integer :: i, j, k, km

    ! Note that for performance reasons we are not
    ! initializing the unused values in this array.

    real(rt) :: partials(0:11)
    real(rt) :: x, y, z, hi, lo

    ! j keeps track of how many entries in partials are actually used.
    ! The algorithm we model this off of, written in Python, simply
    ! deletes array entries at the end of every outer loop iteration.
    ! The Fortran equivalent to this might be to just zero them out,
    ! but this results in a huge performance hit given how often
    ! this routine is called during in a burn. So we opt instead to
    ! just track how many of the values are meaningful, which j does
    ! automatically, and ignore any data in the remaining slots.

    j = 0

    ! The first partial is just the first term.
    partials(j) = array(1)

    do i = 2, 12

       km = j
       j = 0

       x = array(i)

       do k = 0, km
          y = partials(k)

          if (abs(x) < abs(y)) then
             ! Swap x, y
             z = y
             y = x
             x = z
          endif

          hi = x + y
          lo = y - (hi - x)

          if (lo .ne. 0.0_rt) then
             partials(j) = lo
             j = j + 1
          endif

          x = hi

       enddo

       partials(j) = x

    enddo

    esum = sum(partials(0:j))

  end function esum12


#ifdef CUDA
  attributes(device) &
#endif
  function esum13(array) result(esum)

    !$acc routine seq

    use bl_error_module, only: bl_error
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: array(:)
    real(rt) :: esum

    integer :: i, j, k, km

    ! Note that for performance reasons we are not
    ! initializing the unused values in this array.

    real(rt) :: partials(0:12)
    real(rt) :: x, y, z, hi, lo

    ! j keeps track of how many entries in partials are actually used.
    ! The algorithm we model this off of, written in Python, simply
    ! deletes array entries at the end of every outer loop iteration.
    ! The Fortran equivalent to this might be to just zero them out,
    ! but this results in a huge performance hit given how often
    ! this routine is called during in a burn. So we opt instead to
    ! just track how many of the values are meaningful, which j does
    ! automatically, and ignore any data in the remaining slots.

    j = 0

    ! The first partial is just the first term.
    partials(j) = array(1)

    do i = 2, 13

       km = j
       j = 0

       x = array(i)

       do k = 0, km
          y = partials(k)

          if (abs(x) < abs(y)) then
             ! Swap x, y
             z = y
             y = x
             x = z
          endif

          hi = x + y
          lo = y - (hi - x)

          if (lo .ne. 0.0_rt) then
             partials(j) = lo
             j = j + 1
          endif

          x = hi

       enddo

       partials(j) = x

    enddo

#ifndef ACC
#ifndef CUDA    
    if (j > n - 1) then
       call amrex_error("Error: too many partials created in esum.")
    endif
#endif
#endif

    esum = sum(partials(0:j))

  end function esum13


#ifdef CUDA
  attributes(device) &
#endif
  function esum15(array) result(esum)

    !$acc routine seq

    use bl_error_module, only: bl_error
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: array(:)
    real(rt) :: esum

    integer :: i, j, k, km

    ! Note that for performance reasons we are not
    ! initializing the unused values in this array.

    real(rt) :: partials(0:14)
    real(rt) :: x, y, z, hi, lo

    ! j keeps track of how many entries in partials are actually used.
    ! The algorithm we model this off of, written in Python, simply
    ! deletes array entries at the end of every outer loop iteration.
    ! The Fortran equivalent to this might be to just zero them out,
    ! but this results in a huge performance hit given how often
    ! this routine is called during in a burn. So we opt instead to
    ! just track how many of the values are meaningful, which j does
    ! automatically, and ignore any data in the remaining slots.

    j = 0

    ! The first partial is just the first term.
    partials(j) = array(1)

    do i = 2, 15

       km = j
       j = 0

       x = array(i)

       do k = 0, km
          y = partials(k)

          if (abs(x) < abs(y)) then
             ! Swap x, y
             z = y
             y = x
             x = z
          endif

          hi = x + y
          lo = y - (hi - x)

          if (lo .ne. 0.0_rt) then
             partials(j) = lo
             j = j + 1
          endif

          x = hi

       enddo

       partials(j) = x

    enddo

    esum = sum(partials(0:j))

  end function esum15


#ifdef CUDA
  attributes(device) &
#endif
  function esum17(array) result(esum)

    !$acc routine seq

    use bl_error_module, only: bl_error
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: array(:)
    real(rt) :: esum

    integer :: i, j, k, km

    ! Note that for performance reasons we are not
    ! initializing the unused values in this array.

    real(rt) :: partials(0:16)
    real(rt) :: x, y, z, hi, lo

    ! j keeps track of how many entries in partials are actually used.
    ! The algorithm we model this off of, written in Python, simply
    ! deletes array entries at the end of every outer loop iteration.
    ! The Fortran equivalent to this might be to just zero them out,
    ! but this results in a huge performance hit given how often
    ! this routine is called during in a burn. So we opt instead to
    ! just track how many of the values are meaningful, which j does
    ! automatically, and ignore any data in the remaining slots.

    j = 0

    ! The first partial is just the first term.
    partials(j) = array(1)

    do i = 2, 17

       km = j
       j = 0

       x = array(i)

       do k = 0, km
          y = partials(k)

          if (abs(x) < abs(y)) then
             ! Swap x, y
             z = y
             y = x
             x = z
          endif

          hi = x + y
          lo = y - (hi - x)

          if (lo .ne. 0.0_rt) then
             partials(j) = lo
             j = j + 1
          endif

          x = hi

       enddo

       partials(j) = x

    enddo

    esum = sum(partials(0:j))

  end function esum17


#ifdef CUDA
  attributes(device) &
#endif
  function esum20(array) result(esum)

    !$acc routine seq

    use bl_error_module, only: bl_error
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: array(:)
    real(rt) :: esum

    integer :: i, j, k, km

    ! Note that for performance reasons we are not
    ! initializing the unused values in this array.

    real(rt) :: partials(0:19)
    real(rt) :: x, y, z, hi, lo

    ! j keeps track of how many entries in partials are actually used.
    ! The algorithm we model this off of, written in Python, simply
    ! deletes array entries at the end of every outer loop iteration.
    ! The Fortran equivalent to this might be to just zero them out,
    ! but this results in a huge performance hit given how often
    ! this routine is called during in a burn. So we opt instead to
    ! just track how many of the values are meaningful, which j does
    ! automatically, and ignore any data in the remaining slots.

    j = 0

    ! The first partial is just the first term.
    partials(j) = array(1)

    do i = 2, 20

       km = j
       j = 0

       x = array(i)

       do k = 0, km
          y = partials(k)

          if (abs(x) < abs(y)) then
             ! Swap x, y
             z = y
             y = x
             x = z
          endif

          hi = x + y
          lo = y - (hi - x)

          if (lo .ne. 0.0_rt) then
             partials(j) = lo
             j = j + 1
          endif

          x = hi

       enddo

       partials(j) = x

    enddo

    esum = sum(partials(0:j))

  end function esum20


#ifdef CUDA
  attributes(device) &
#endif
  function esum25(array) result(esum)

    !$acc routine seq

    use bl_error_module, only: bl_error
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: array(:)
    real(rt) :: esum

    integer :: i, j, k, km

    ! Note that for performance reasons we are not
    ! initializing the unused values in this array.

    real(rt) :: partials(0:24)
    real(rt) :: x, y, z, hi, lo

    ! j keeps track of how many entries in partials are actually used.
    ! The algorithm we model this off of, written in Python, simply
    ! deletes array entries at the end of every outer loop iteration.
    ! The Fortran equivalent to this might be to just zero them out,
    ! but this results in a huge performance hit given how often
    ! this routine is called during in a burn. So we opt instead to
    ! just track how many of the values are meaningful, which j does
    ! automatically, and ignore any data in the remaining slots.

    j = 0

    ! The first partial is just the first term.
    partials(j) = array(1)

    do i = 2, 25

       km = j
       j = 0

       x = array(i)

       do k = 0, km
          y = partials(k)

          if (abs(x) < abs(y)) then
             ! Swap x, y
             z = y
             y = x
             x = z
          endif

          hi = x + y
          lo = y - (hi - x)

          if (lo .ne. 0.0_rt) then
             partials(j) = lo
             j = j + 1
          endif

          x = hi

       enddo

       partials(j) = x

    enddo

    esum = sum(partials(0:j))

  end function esum25


#ifdef CUDA
  attributes(device) &
#endif
  function esum26(array) result(esum)

    !$acc routine seq

    use bl_error_module, only: bl_error
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: array(:)
    real(rt) :: esum

    integer :: i, j, k, km

    ! Note that for performance reasons we are not
    ! initializing the unused values in this array.

    real(rt) :: partials(0:25)
    real(rt) :: x, y, z, hi, lo

    ! j keeps track of how many entries in partials are actually used.
    ! The algorithm we model this off of, written in Python, simply
    ! deletes array entries at the end of every outer loop iteration.
    ! The Fortran equivalent to this might be to just zero them out,
    ! but this results in a huge performance hit given how often
    ! this routine is called during in a burn. So we opt instead to
    ! just track how many of the values are meaningful, which j does
    ! automatically, and ignore any data in the remaining slots.

    j = 0

    ! The first partial is just the first term.
    partials(j) = array(1)

    do i = 2, 26

       km = j
       j = 0

       x = array(i)

       do k = 0, km
          y = partials(k)

          if (abs(x) < abs(y)) then
             ! Swap x, y
             z = y
             y = x
             x = z
          endif

          hi = x + y
          lo = y - (hi - x)

          if (lo .ne. 0.0_rt) then
             partials(j) = lo
             j = j + 1
          endif

          x = hi

       enddo

       partials(j) = x

    enddo

    esum = sum(partials(0:j))

  end function esum26


#ifdef CUDA
  attributes(device) &
#endif
  function esum3(array) result(esum)

    !$acc routine seq

    use bl_error_module, only: bl_error
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: array(:)
    real(rt) :: esum

    integer :: i, j, k, km

    ! Note that for performance reasons we are not
    ! initializing the unused values in this array.

    real(rt) :: partials(0:2)
    real(rt) :: x, y, z, hi, lo

    ! j keeps track of how many entries in partials are actually used.
    ! The algorithm we model this off of, written in Python, simply
    ! deletes array entries at the end of every outer loop iteration.
    ! The Fortran equivalent to this might be to just zero them out,
    ! but this results in a huge performance hit given how often
    ! this routine is called during in a burn. So we opt instead to
    ! just track how many of the values are meaningful, which j does
    ! automatically, and ignore any data in the remaining slots.

    j = 0

    ! The first partial is just the first term.
    partials(j) = array(1)

    do i = 2, 3

       km = j
       j = 0

       x = array(i)

       do k = 0, km
          y = partials(k)

          if (abs(x) < abs(y)) then
             ! Swap x, y
             z = y
             y = x
             x = z
          endif

          hi = x + y
          lo = y - (hi - x)

          if (lo .ne. 0.0_rt) then
             partials(j) = lo
             j = j + 1
          endif

          x = hi

       enddo

       partials(j) = x

    enddo

    esum = sum(partials(0:j))

  end function esum3


#ifdef CUDA
  attributes(device) &
#endif
  function esum4(array) result(esum)

    !$acc routine seq

    use bl_error_module, only: bl_error
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: array(:)
    real(rt) :: esum

    integer :: i, j, k, km

    ! Note that for performance reasons we are not
    ! initializing the unused values in this array.

    real(rt) :: partials(0:3)
    real(rt) :: x, y, z, hi, lo

    ! j keeps track of how many entries in partials are actually used.
    ! The algorithm we model this off of, written in Python, simply
    ! deletes array entries at the end of every outer loop iteration.
    ! The Fortran equivalent to this might be to just zero them out,
    ! but this results in a huge performance hit given how often
    ! this routine is called during in a burn. So we opt instead to
    ! just track how many of the values are meaningful, which j does
    ! automatically, and ignore any data in the remaining slots.

    j = 0

    ! The first partial is just the first term.
    partials(j) = array(1)

    do i = 2, 4

       km = j
       j = 0

       x = array(i)

       do k = 0, km
          y = partials(k)

          if (abs(x) < abs(y)) then
             ! Swap x, y
             z = y
             y = x
             x = z
          endif

          hi = x + y
          lo = y - (hi - x)

          if (lo .ne. 0.0_rt) then
             partials(j) = lo
             j = j + 1
          endif

          x = hi

       enddo

       partials(j) = x

    enddo

    esum = sum(partials(0:j))

  end function esum4


#ifdef CUDA
  attributes(device) &
#endif
  function esum5(array) result(esum)

    !$acc routine seq

    use bl_error_module, only: bl_error
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: array(:)
    real(rt) :: esum

    integer :: i, j, k, km

    ! Note that for performance reasons we are not
    ! initializing the unused values in this array.

    real(rt) :: partials(0:4)
    real(rt) :: x, y, z, hi, lo

    ! j keeps track of how many entries in partials are actually used.
    ! The algorithm we model this off of, written in Python, simply
    ! deletes array entries at the end of every outer loop iteration.
    ! The Fortran equivalent to this might be to just zero them out,
    ! but this results in a huge performance hit given how often
    ! this routine is called during in a burn. So we opt instead to
    ! just track how many of the values are meaningful, which j does
    ! automatically, and ignore any data in the remaining slots.

    j = 0

    ! The first partial is just the first term.
    partials(j) = array(1)

    do i = 2, 5

       km = j
       j = 0

       x = array(i)

       do k = 0, km
          y = partials(k)

          if (abs(x) < abs(y)) then
             ! Swap x, y
             z = y
             y = x
             x = z
          endif

          hi = x + y
          lo = y - (hi - x)

          if (lo .ne. 0.0_rt) then
             partials(j) = lo
             j = j + 1
          endif

          x = hi

       enddo

       partials(j) = x

    enddo

    esum = sum(partials(0:j))

  end function esum5


#ifdef CUDA
  attributes(device) &
#endif
  function esum6(array) result(esum)

    !$acc routine seq

    use bl_error_module, only: bl_error
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: array(:)
    real(rt) :: esum

    integer :: i, j, k, km

    ! Note that for performance reasons we are not
    ! initializing the unused values in this array.

    real(rt) :: partials(0:5)
    real(rt) :: x, y, z, hi, lo

    ! j keeps track of how many entries in partials are actually used.
    ! The algorithm we model this off of, written in Python, simply
    ! deletes array entries at the end of every outer loop iteration.
    ! The Fortran equivalent to this might be to just zero them out,
    ! but this results in a huge performance hit given how often
    ! this routine is called during in a burn. So we opt instead to
    ! just track how many of the values are meaningful, which j does
    ! automatically, and ignore any data in the remaining slots.

    j = 0

    ! The first partial is just the first term.
    partials(j) = array(1)

    do i = 2, 6

       km = j
       j = 0

       x = array(i)

       do k = 0, km
          y = partials(k)

          if (abs(x) < abs(y)) then
             ! Swap x, y
             z = y
             y = x
             x = z
          endif

          hi = x + y
          lo = y - (hi - x)

          if (lo .ne. 0.0_rt) then
             partials(j) = lo
             j = j + 1
          endif

          x = hi

       enddo

       partials(j) = x

    enddo

    esum = sum(partials(0:j))

  end function esum6


#ifdef CUDA
  attributes(device) &
#endif
  function esum7(array) result(esum)

    !$acc routine seq

    use bl_error_module, only: bl_error
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: array(:)
    real(rt) :: esum

    integer :: i, j, k, km

    ! Note that for performance reasons we are not
    ! initializing the unused values in this array.

    real(rt) :: partials(0:6)
    real(rt) :: x, y, z, hi, lo

    ! j keeps track of how many entries in partials are actually used.
    ! The algorithm we model this off of, written in Python, simply
    ! deletes array entries at the end of every outer loop iteration.
    ! The Fortran equivalent to this might be to just zero them out,
    ! but this results in a huge performance hit given how often
    ! this routine is called during in a burn. So we opt instead to
    ! just track how many of the values are meaningful, which j does
    ! automatically, and ignore any data in the remaining slots.

    j = 0

    ! The first partial is just the first term.
    partials(j) = array(1)

    do i = 2, 7

       km = j
       j = 0

       x = array(i)

       do k = 0, km
          y = partials(k)

          if (abs(x) < abs(y)) then
             ! Swap x, y
             z = y
             y = x
             x = z
          endif

          hi = x + y
          lo = y - (hi - x)

          if (lo .ne. 0.0_rt) then
             partials(j) = lo
             j = j + 1
          endif

          x = hi

       enddo

       partials(j) = x

    enddo

    esum = sum(partials(0:j))

  end function esum7


#ifdef CUDA
  attributes(device) &
#endif
  function esum8(array) result(esum)

    !$acc routine seq

    use bl_error_module, only: bl_error
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: array(:)
    real(rt) :: esum

    integer :: i, j, k, km

    ! Note that for performance reasons we are not
    ! initializing the unused values in this array.

    real(rt) :: partials(0:7)
    real(rt) :: x, y, z, hi, lo

    ! j keeps track of how many entries in partials are actually used.
    ! The algorithm we model this off of, written in Python, simply
    ! deletes array entries at the end of every outer loop iteration.
    ! The Fortran equivalent to this might be to just zero them out,
    ! but this results in a huge performance hit given how often
    ! this routine is called during in a burn. So we opt instead to
    ! just track how many of the values are meaningful, which j does
    ! automatically, and ignore any data in the remaining slots.

    j = 0

    ! The first partial is just the first term.
    partials(j) = array(1)

    do i = 2, 8

       km = j
       j = 0

       x = array(i)

       do k = 0, km
          y = partials(k)

          if (abs(x) < abs(y)) then
             ! Swap x, y
             z = y
             y = x
             x = z
          endif

          hi = x + y
          lo = y - (hi - x)

          if (lo .ne. 0.0_rt) then
             partials(j) = lo
             j = j + 1
          endif

          x = hi

       enddo

       partials(j) = x

    enddo

    esum = sum(partials(0:j))

  end function esum8


#ifdef CUDA
  attributes(device) &
#endif
  function esum9(array) result(esum)

    !$acc routine seq

    use bl_error_module, only: bl_error
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: array(:)
    real(rt) :: esum

    integer :: i, j, k, km

    ! Note that for performance reasons we are not
    ! initializing the unused values in this array.

    real(rt) :: partials(0:8)
    real(rt) :: x, y, z, hi, lo

    ! j keeps track of how many entries in partials are actually used.
    ! The algorithm we model this off of, written in Python, simply
    ! deletes array entries at the end of every outer loop iteration.
    ! The Fortran equivalent to this might be to just zero them out,
    ! but this results in a huge performance hit given how often
    ! this routine is called during in a burn. So we opt instead to
    ! just track how many of the values are meaningful, which j does
    ! automatically, and ignore any data in the remaining slots.

    j = 0

    ! The first partial is just the first term.
    partials(j) = array(1)

    do i = 2, 9

       km = j
       j = 0

       x = array(i)

       do k = 0, km
          y = partials(k)

          if (abs(x) < abs(y)) then
             ! Swap x, y
             z = y
             y = x
             x = z
          endif

          hi = x + y
          lo = y - (hi - x)

          if (lo .ne. 0.0_rt) then
             partials(j) = lo
             j = j + 1
          endif

          x = hi

       enddo

       partials(j) = x

    enddo

    esum = sum(partials(0:j))

  end function esum9


end module microphysics_math_module
