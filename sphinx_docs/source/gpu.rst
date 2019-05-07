****
GPUs
****

General Ideas
=============

The current strategy for StarKiller Microphysics is to use CUDA
Fortran for all microphysics routines necessary to evaluate the EOS
and integrate a reaction network so that the high level wrappers `eos`
and `burner` can be called from a global CUDA kernel in which each
cell in a grid is assigned to a GPU thread.

This global kernel loops through the grid, creating thread-local
variables from the derived types in `Microphysics/interfaces` with
initial conditions from the grid. These interface types are passed to
the microphysics routines by each thread and the desired results are
saved back to the grid.

Reaction Network Integration
============================

There are currently two ways to integrate reaction networks using
GPUs: CUDA Fortran VODE and CVODE.

We recommend using CUDA Fortran VODE for the networks in Microphysics
for production simulations at this time.

VODE90 - CUDA Fortran VODE
--------------------------

This integrator is located in `Microphysics/integration/VODE90` and
follows the paradigm described above where it is designed to be called
from a single GPU thread and uses only local thread memory. It is a
port of the legacy Fortran 77 `dvode` code to CUDA Fortran.

By necessity, this port reorganized the local data layout in `VODE` to
eliminate `common` blocks in favor of derived types and eliminate
assumed-size or assumed-shape arrays in favor of arrays with sizes
explicitly known as compile time parameters. This allows the CUDA
Fortran compiler to avoid costly global memory allocations at runtime.

To reduce memory, this port also eliminated support for VODE's
explicit Adams-Moulton integration, as we use only the implicit BDF
method provided by VODE for reaction networks. This code has been
tested against the original `dvode` integrator in
`Microphysics/integration/VODE` and yields the same results.

Standalone Test
^^^^^^^^^^^^^^^

There is a standalone test problem for the CUDA Fortran VODE
integrator that integrates a grid of cells, with the original VODE
test problem in each cell.

To run this problem, see the test setup in
`Microphysics/integration/VODE90/cuVODE/test` and the Readme file
located there.

Reaction Network Test
^^^^^^^^^^^^^^^^^^^^^

There is also a test setup that exercises this integrator (as well as
other integrators in Microphysics) on our reaction networks while
optionally calling the EOS in the right hand side evaluation. This
test thus closely emulates the environment of a production simulation.

To run this test, see the setup in `Microphysics/unit_test/test_react`
and the included Readme file.

CVODE
-----

Microphysics also includes a test setup and interfaces for using CVODE
to integrate an entire grid of cells together either serially, on the
GPU using CVODE's Krylov linear solver, or on the GPU using cuSOLVER's
batched sparse QR linear solver.

This approach does not follow the GPU paradigm of one thread per cell
as in the case of CUDA Fortran VODE. The CVODE interface packs the
network ODEs from an entire grid into a single vector of
ODEs. Logically, the linear system that arises when solving the BDF
equation for such a system is block-diagonal. All cells are integrated
together in lock-step as if there were only one system of ODEs.

CVODE implements GPU parallelism by expressing vector operations as
operations on `NVector` data structures that exist in GPU
memory. Arithmetic operations on `NVector` objects are implemented as
CUDA kernels operating on the GPU-resident `NVector` data.

This interface is under development and testing and is not ready for
production use.

Reaction Network Test
^^^^^^^^^^^^^^^^^^^^^

There is a test setup that exercises the CVODE interface in
`Microphysics/unit_test/test_cvode_react` along with a Readme file
with details for compiling and running.

It is recommended to experiment with the different linear solvers to
determine which is best for a given network. For `aprox13` the sparse
direct solver in the cuSOLVER toolkit is preferred.

At present, this setup requires a modified version of CVODE
incorporating the interface to the batched sparse QR linear solver in
the cuSOLVER toolkit. The modified version of CVODE is located here:

`<https://github.com/dwillcox/cvode-4.0.0-development>`_

Note that because this integrator interface effectively couples many
independent cells together, the error norm check in the nonlinear
solver may not be sensitive to relatively large convergence errors in
only a small subset of the grid cells. As a result, it can occur that
this integration scheme will not yield results identical to
integrating each cell independently of the others.

Common Compiler Errors
======================

PGI - OpenACC
-------------

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
