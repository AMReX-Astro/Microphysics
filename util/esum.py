import os
import re

esum_template = """
#ifdef CUDA
  attributes(device) &
#endif
  function esum@NUM@(array) result(esum)

    !$acc routine seq

    use bl_error_module, only: bl_error
    use bl_types, only: dp_t

    implicit none

    real(dp_t), intent(in) :: array(:)
    real(dp_t) :: esum

    integer :: i, j, k, km

    ! Note that for performance reasons we are not
    ! initializing the unused values in this array.

    real(dp_t) :: partials(0:@NUM_MINUS_ONE@)
    real(dp_t) :: x, y, z, hi, lo

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

    do i = 2, n

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

          if (lo .ne. 0.0_dp_t) then
             partials(j) = lo
             j = j + 1
          endif

          x = hi

       enddo

       partials(j) = x

    enddo

    esum = sum(partials(0:j))

  end function esum@NUM@

"""

module_start = """
module microphysics_math_module

  implicit none

  public

contains
"""

module_end = """
end module microphysics_math_module
"""

file_extensions = [".f90", ".F90"]

esum_needed = []

top_dir = os.getcwd()

# find the files with proper extension
potential_files = []
for dir_name, _, files in os.walk(top_dir):
    for ext in file_extensions:
        potential_files += [os.path.join(dir_name, f) for f in files
                            if f.endswith(ext) and os.path.isfile(os.path.join(dir_name, f))]

# now search those files for the ones that contains esum
esum_files = []
for pf in potential_files:
    with open(pf, "r") as fh:
        for line in fh:
            if "esum" in line:
                esum_files.append(pf)
                break


# now do the real work (probably could have done that above :)
esum_re = re.compile("(esum).*?(\\()((?:[a-z][a-z0-9_]*)).*?(\\d+)(\\))")

for ef in esum_files:
    print("working on {}".format(ef))

    # back it up
    os.rename(ef, "{}_orig".format(ef))

    with open(ef, "w") as fnew, open("{}_orig".format(ef), "r") as fold:
        for line in fold:
            eout = esum_re.search(line.strip())
            if eout:
                full_esum = eout.group(0)
                var = eout.group(3)
                num = eout.group(4)

                line = line.replace(full_esum, "esum{}({})".format(num, var))

                if not num in esum_needed:
                    esum_needed.append(num)

            fnew.write(line)


with open("esum_module.f90", "w") as ef:

    ef.write(module_start)

    for num in sorted(set(esum_needed)):
        ef.write(esum_template.replace("@NUM@", num).replace("@NUM_MINUS_ONE@", str(int(num)-1)))

    ef.write(module_end)

