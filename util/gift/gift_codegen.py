"""
Python port of public_gift.f from https://cococubed.com/code_pages/gift.shtml

Modified to also output Python and AMReX C++ code.

Original header:
-----------------------------------------------------------------------
    g i f t    is a program generator, which outputs fortran code for
               g_auss i_nversion of similar linear systems of equations
               by a f_ormal t_reatment considering only the non-zero
               elements of the corresponding matrices.

    authors:
      original version:    h. falk,  r.g. wolff           1983
      modified by          e. mueller              april  1996
                           f. x. timmes            august 2000

    notes: mueller added the code for a vectorized version.
           fxt removed this vectorization for better performance on
           parallel machines, and minimized the number of divides.
           for n equations, there are now only n divides.
           if you want the vectorized version, contact fxt.

           to use this program you will need to modify the code in
           routine setmat to reflect the sparsity pattern of your matrix,
           and set the logical size of your matrix a few lines below.
           you shouldn't have to do much else.

           the output of running this program is set to gift_test.f.
           to use a gift generated routine, it is important to note
           that the size of the input matrix is a(n,n+1). the extra
           column is where the right hand side is stored on input,
           and contains the solution to a*x=b on output.

----------------------------------------------------------------------
"""

import argparse
from typing import Any, Callable

import numpy as np
import numpy.typing as npt

FORMATS: dict[str, dict[int, str | Callable[..., str]]] = {}
FORMATS["fortran"] = {
    1001: """\
      subroutine {:s}(ab,n1,n2)
      implicit none
      integer n1,n2
      double precision ab(n1,n2)""",
    3001: """\
      subroutine {:s}{:3d}(ab,n1,n2)
      implicit none
      integer n1,n2
      double precision ab(n1,n2),c,tmp({:3d})""",
    3002: "      call {:s}{:3d}(ab,n1,n2)",
    # gfortran and ifort appear to strip lines using only X, so no leading space here
    6003: "\n        tmp({0:3d})   = 1.0d0/ab({0:3d},{0:3d})",
    1003: "        c           = ab({0:3d},{1:3d}) * tmp({1:4d})",
    1004: "        ab({0:3d},{1:3d}) = ab({0:3d},{1:3d}) - ab({2:3d},{1:3d}) * c",
    1005: "        ab({0:3d},{1:3d}) = ab({0:3d},{1:3d}) - ab({2:3d},{1:3d}) * c",
    # no spaces on empty lines here either
    1007: """\n\nc backsubstitution
        tmp({0:3d})     = 1.0d0/ab({0:3d},{0:3d})
        ab({0:3d},{1:3d}) =   ab({0:3d},{1:3d}) * tmp({0:4d})""",
    1008: "        ab({0:3d},{1:3d}) = ( ab({0:3d},{1:3d})",
    1009: "     &                  - ab({0:3d},{1:3d}) * ab({1:3d},{2:3d})",
    2001: "     &                                              ) *tmp({:4d})",
    4001: "     &                                              )",
    2002: "      return\n      end",
}

FORMATS["python"] = {
    7001: "import numpy as np",
    1001: """

def {:s}(ab, n1, n2):
    del n1, n2  # unused""",
    3001: """

def {:s}{:3d}(ab):
    tmp = np.zeros({:d}, dtype=np.double)""",
    3002: "    {:s}{:3d}(ab)",
    6003: lambda i: f"\n    tmp[{i-1:3d}]     = 1.0e0 / ab[{i-1:3d}, {i-1:3d}]",
    1003: lambda i,j: f"    c            = ab[{i-1:3d}, {j-1:3d}] * tmp[{j-1:3d}]",
    1004: lambda i,j,k: f"    ab[{i-1:3d}, {j-1:3d}] = ab[{i-1:3d}, {j-1:3d}] - ab[{k-1:3d}, {j-1:3d}] * c",
    1005: lambda i,j,k: f"    ab[{i-1:3d}, {j-1:3d}] = ab[{i-1:3d}, {j-1:3d}] - ab[{k-1:3d}, {j-1:3d}] * c",
    1007: lambda i,j: f"""\n
    # backsubstitution
    tmp[{i-1:3d}]     = 1.0e0 / ab[{i-1:3d}, {i-1:3d}]
    ab[{i-1:3d}, {j-1:3d}] =   ab[{i-1:3d}, {j-1:3d}] * tmp[{i-1:3d}]""",
    1008: lambda i,j: f"    ab[{i-1:3d}, {j-1:3d}] = ( ab[{i-1:3d}, {j-1:3d}]",
    1009: lambda i,j,k: f"                     - ab[{i-1:3d}, {j-1:3d}] * ab[{j-1:3d}, {k-1:3d}]",
    2001: lambda i: f"                                                   ) * tmp[{i-1:3d}]",
    4001: "                                                   )",
}

FORMATS["cpp"] = {
    7001: """\
#ifndef GIFT_TEST_H
#define GIFT_TEST_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
""",
    7002: "\n#endif // GIFT_TEST_H",
    1001: """\

template <int neqs>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void {:s}(amrex::Array2D<amrex::Real, 1, neqs, 1, neqs+1> &ab)
{{""",
    3001: """\

template <int neqs>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void {:s}{:3d}(amrex::Array2D<amrex::Real, 1, neqs, 1, neqs+1> &ab)
{{
    amrex::Array1D<amrex::Real, 1, {:d}> tmp;
    amrex::Real c;""",
    3002: "    {:s}{:3d}(ab);",
    6003: "\n    tmp({0:3d})    = 1.0e0/ab({0:3d},{0:3d});",
    1003: "    c           = ab({0:3d},{1:3d}) * tmp({1:4d});",
    1004: "    ab({0:3d},{1:3d}) = ab({0:3d},{1:3d}) - ab({2:3d},{1:3d}) * c;",
    1005: "    ab({0:3d},{1:3d}) = ab({0:3d},{1:3d}) - ab({2:3d},{1:3d}) * c;",
    1007: """\n
    // backsubstitution
    tmp({0:3d})    = 1.0e0/ab({0:3d},{0:3d});
    ab({0:3d},{1:3d}) =   ab({0:3d},{1:3d}) * tmp({0:4d});""",
    1008: "    ab({0:3d},{1:3d}) = ( ab({0:3d},{1:3d})",
    1009: "                    - ab({0:3d},{1:3d}) * ab({1:3d},{2:3d})",
    2001: "                                                ) * tmp({:4d});",
    4001: "                                                );",
    2002: "}}",
}

EXTENSIONS = {"fortran": ".f", "python": ".py", "cpp": ".H"}


def gift(a: npt.NDArray[np.bool_], b: npt.NDArray[np.bool_], lang: str) -> None:
    # pylint: disable=too-many-locals, too-many-statements, too-many-branches
    # pylint: disable=too-many-nested-blocks

    assert a.ndim == 2
    n = a.shape[0]
    assert a.shape[1] == n

    assert b.ndim == 1
    assert b.shape[0] == n

    ctop = "gift_test"
    csub = "gsub_test_"

    nplus = n + 1

    # write out a map of the matrix
    npri = min(n, 132)
    print("  ", f"map of matrix n={n:4d}\n", sep="")

    for i in range(npri):
        print(" ", *('x' if a[i, j] else '.' for j in range(npri)), sep="")

    # now start writing the fortran code

    ilabel = 1  # do-loop label counter

    itlsub = 1_000_000  # maximum number of do-loops per subroutine
    nosub1 = 100  # number of first subroutine

    numst = 0  # statement number counter within a do-loop
    lsize = 1_000_000  # maximum number of statements per do-loop
    icount = 0  # block number counter within a do-loop
    nblk = 1_000_000  # maximum number of blocks per do-loop
    totalst = 0  # total number of statements in the program

    nosub = nosub1  # initialize subroutine number counter

    string = ctop + EXTENSIONS[lang]
    f = open(string, mode="w")  # pylint: disable=consider-using-with

    def fwrite(fmt_num: int, *args: Any) -> None:
        if fmt_num in FORMATS[lang]:
            fmt = FORMATS[lang][fmt_num]
            if isinstance(fmt, str):
                fmt = fmt.format
            print(fmt(*args), file=f)

    fwrite(7001)  # write file header

    fwrite(3001, csub, nosub, n)
    for jj in range(2, n + 1):
        j = n + 2 - jj
        fwrite(6003, j)

        for ii in range(1, j):
            i = j - 1 + 1 - ii

            if a[i - 1, j - 1]:
                if b[j - 1]:
                    if not b[i - 1]:
                        b[i - 1] = b[j - 1]
                    if icount >= nblk and numst > lsize:
                        ilabel += 1  # advance do-loop label

                        if ilabel % itlsub == 0:
                            fwrite(2002)  # write return / end
                            nosub += 1  # advance subroutine label
                            fwrite(3001, csub, nosub, n)  # write subroutine header
                        totalst += numst
                        numst = 0
                        icount = 0

                    icount += 1

                    fwrite(1003, i, j)

                    numst += 1
                    if b[i - 1]:
                        fwrite(1004, i, nplus, j)
                        numst += 1

                for k in range(1, j):
                    if not a[j - 1, k - 1]:
                        continue
                    if not a[i - 1, k - 1]:
                        a[i - 1, k - 1] = a[j - 1, k - 1]
                    if not a[i - 1, k - 1]:
                        continue
                    fwrite(1005, i, k, j)
                    numst += 1

            for k in range(j, n + 1):
                a[i - 1, k - 1] = False

    ilabel += 1  # advance do-loop counter

    if ilabel % itlsub == 0:
        fwrite(2002)  # write return / end
        nosub += 1  # advance subroutine label
        fwrite(3001, csub, nosub, n)

    icount = 0
    totalst += numst
    numst = 0

    # backsubstitution
    ione = 1
    fwrite(1007, ione, nplus)

    for l in range(2, n + 1):  # noqa: E741
        if icount >= nblk and numst > lsize:
            ilabel += 1  # advance do-loop counter

            if ilabel % itlsub == 0:
                fwrite(2001)  # write return / end
                nosub += 1  # advance subroutine label
                fwrite(3001, csub, nosub, n)  # write subroutine header

            icount = 0
            totalst += numst
            numst = 0

        icount += 1
        fwrite(1008, l, nplus)
        numst += 1

        istmt = 0
        for k in range(1, l):
            if not a[l - 1, k - 1]:
                continue
            istmt += 1

            if istmt % 18 == 0:  # limits continuation lines
                fwrite(4001)  # writes closing ")"
                numst += 1
                icount += 1
                fwrite(1008, l, nplus)  # new statement
                numst += 1

            fwrite(1009, l, k, nplus)  # continuation card
            numst += 1

        fwrite(2001, l)
        numst += 1

    totalst += numst

    fwrite(2002)  # write return / end

    fwrite(1001, ctop)  # write master subroutine header
    for k in range(nosub1, nosub + 1):
        fwrite(3002, csub, k)  # write call to sub-routine
    fwrite(2002)  # write return / end

    fwrite(7002)  # write file footer

    f.close()

    print("  ")
    print(f"  total number of do-loops:     {ilabel:12d}")
    print(f"  total number of subroutines:  {nosub - nosub1 + 1:12d}")
    print(f"  subroutines with prefix: {csub}")
    print(f"  and subroutine gift have been written to {string:20s}")
    print(f"  last subroutine is {csub}{nosub:12d}")
    print(f"  total statements: {totalst}")


def gen_matrices(n: int) -> tuple[npt.NDArray[np.bool_], npt.NDArray[np.bool_]]:
    """Sets the nonzero patterns of the dense matrix a and the dense right hand side b.

    a has shape(n, n), b has shape (n).
    Returns (a, b).
    """

    # initialize the matrix and right hand side
    a = np.zeros((n, n), dtype=bool)
    b = np.ones(n, dtype=bool)

    # set the pattern for a simple tridiagonal
    band_size = 3
    for i in range(band_size//2 + 1):
        flat_diag = np.ones(n-i, dtype=bool)
        a += np.diag(flat_diag, k=i)
        a += np.diag(flat_diag, k=-i)

    return a, b


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "language",
        choices=list(FORMATS) + [[]],
        nargs="*",
        help="language(s) to generate code for",
    )
    args = parser.parse_args()

    # aprox13 (output by test_linear_algebra)
    a = np.genfromtxt("""
   *  He4  C12  O16 Ne20 Mg24 Si28  S32 Ar36 Ca40 Ti44 Cr48 Fe52 Ni56 enuc
 He4    1    1    1    1    1    1    1    1    1    1    1    1    1    1
 C12    1    1    1    1    1    1    0    0    0    0    0    0    0    1
 O16    1    1    1    1    1    1    1    0    0    0    0    0    0    1
Ne20    1    1    1    1    1    0    0    0    0    0    0    0    0    1
Mg24    1    1    1    1    1    1    0    0    0    0    0    0    0    1
Si28    1    1    1    0    1    1    1    0    0    0    0    0    0    1
 S32    1    0    1    0    0    1    1    1    0    0    0    0    0    1
Ar36    1    0    0    0    0    0    1    1    1    0    0    0    0    1
Ca40    1    0    0    0    0    0    0    1    1    1    0    0    0    1
Ti44    1    0    0    0    0    0    0    0    1    1    1    0    0    1
Cr48    1    0    0    0    0    0    0    0    0    1    1    1    0    1
Fe52    1    0    0    0    0    0    0    0    0    0    1    1    1    1
Ni56    1    0    0    0    0    0    0    0    0    0    0    1    1    1
enuc    1    1    1    1    1    1    1    1    1    1    1    1    1    1
    """.splitlines(), dtype=int)[1:, 1:].astype(bool)
    b = np.ones(a.shape[0], dtype=bool)

    # try moving enuc to the front
    # a = np.roll(np.roll(a, 1, axis=0), 1, axis=1)

    # non-sparse matrix for comparison
    # a = np.ones_like(a, dtype=bool)

    # a, b = gen_matrices(10)
    if not args.language:
        args.language = list(FORMATS)
    for language in args.language:
        gift(a.copy(), b.copy(), language)


if __name__ == "__main__":
    main()
