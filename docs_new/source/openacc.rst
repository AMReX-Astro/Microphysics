****
GPUs
****

General Ideas
=============

Common Compiler Errors
======================

PGI
---

To get more detailed about what kernels are launching, do:

::

    export PGI_ACC_NOTIFY=1

Compilation errors
^^^^^^^^^^^^^^^^^^

-  *Multiple exit statements*

   ::

       PGF90-S-0155-Compiler failed to translate accelerator region (see -Minfo messages): 
       Unexpected flow graph (../../integration/BS/stiff_ode.F90: 1)
       single_step_rosen:

   This results when you have multiple exit statements from a
   do-loop. You need to consolidate any error/convergence checking in
   a loop and have at most a single exit statement.

   Note: sometimes a do while can help here, but there is a sense
   that while-loops do not perform optimally on GPUs.

-  *ACC routines not in Fortran modules*

   ::

       PGF90-S-0155-Procedures called in a compute region must have acc routine information: 
         dgefa (../../integration/BS/stiff_ode.F90: 711)
       PGF90-S-0155-Accelerator region ignored; see -Minfo messages  
         (../../integration/BS/stiff_ode.F90)

   This occurs when a subroutine relies on another routine that is not part
   of a Fortran 90 module. In this case, even if that routine already has

   ::

       !$acc routine seq

   we need to mark the *calling* routine as well, with:

   ::

       !$acc routine(dgesl) seq

   (e.g., for the Fortran routine dgesl).

Runtime errors
^^^^^^^^^^^^^^

-  *Multi-d array copies*

   ::

       Unhandled builtin: 601 (pgf90_mzero8)
       PGF90-F-0000-Internal compiler error. Unhandled builtin function.       
         0 (../../networks/triple_alpha_plus_cago/actual_rhs.f90: 146)
       PGF90/x86-64 Linux 16.5-0: compilation aborted

   This error results from doing a multi-d array copy (with Fortran
   notation) in GPU code. The fix is to explicitly write out a loop over
   rows.

-  *Illegal memory access*

   ::

       call to cuMemcpyDtoHAsync returned error 700: Illegal address during kernel execution
       call to cuMemFreeHost returned error 700: Illegal address during kernel execution

   This indicates that you went out of bounds in memory access or,
   sometimes it seems, generated some NaNs.

Debugging
=========

cuda-gdb
--------

Basic debugging can be done using cuda-gdb. This will work just
like gdb and can give you the name of a routine where a crash
occurred, but generally doesn’t produce line numbers.
