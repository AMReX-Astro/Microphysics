# Runge-Kutta-Chebyshev Integrator

This is an implementation of the Runge-Kutta-Chebyshev (RKC)
integrator.  This is a port of the `rkc.f` routine from Sommeijer et
al. 1997, as downloaded from Netlib.

Here's the original comments:

```
ABSTRACT: RKC integrates initial value problems for systems of first
order ordinary differential equations.  It is based on a family of
explicit Runge-Kutta-Chebyshev formulas of order two.  The stability
of members of the family increases quadratically in the number of
stages m. An estimate of the spectral radius is used at each step to
select the smallest m resulting in a stable integration. RKC is
appropriate for the solution to modest accuracy of mildly stiff
problems with eigenvalues of Jacobians that are close to the negative
real axis.  For such problems it has the advantages of explicit
one-step methods and very low storage. If it should turn out that RKC
is using m far beyond 100, the problem is not mildly stiff and
alternative methods should be considered.  Answers can be obtained
cheaply anywhere in the interval of integration by means of a
continuous extension evaluated in the subroutine RKCINT.

The initial value problems arising from semi-discretization of
diffusion-dominated parabolic partial differential equations and of
reaction-diffusion equations, especially in two and three spatial
variables, exemplify the problems for which RKC was designed.  Two
example programs, ExA and ExB, are provided that show how to use RKC.

USAGE: RKC integrates a system of NEQN first order ordinary
differential equations specified by a subroutine F from T to TEND.
The initial values at T are input in Y(*).  On all returns from RKC,
Y(*) is an approximate solution at T.  In the computation of Y(*), the
local error has been controlled at each step to satisfy a relative
error tolerance RTOL and absolute error tolerances ATOL(*).  The array
INFO(*) specifies the way the problem is to be solved.  WORK(*) is a
work array. IDID reports success or the reason the computation has
been terminated.

FIRST CALL TO RK

You must provide storage in your calling program for the arrays in the
call list -- Y(NEQN), INFO(4), WORK(8+5*NEQN).  If INFO(2) = 0, you
can reduce the storage for the work array to WORK(8+4*NEQN).  ATOL may
be a scalar or an array.  If it is an array, you must provide storage
for ATOL(NEQN).  You must declare F in an external statement, supply
the subroutine F and the function SPCRAD, and initialize the following
quantities:

  NEQN:  The number of differential equations. Integer.

  T:     The initial point of the integration. Double precision.
         Must be a variable.

  TEND:  The end of the interval of integration.  Double precision.
         TEND may be less than T.

  Y(*):  The initial value of the solution.  Double precision array
         of length NEQN.

  F:     The name of a subroutine for evaluating the differential
         equation.  It must have the form

           subroutine f(neqn,t,y,dy)
           integer          neqn
           double precision t,y(neqn),dy(neqn)
           dy(1)    = ...
           ...
           dy(neqn) = ...
           return
           end

RTOL,
ATOL(*):  At each step of the integration the local error is controlled
          so that its RMS norm is no larger than tolerances RTOL, ATOL(*).
          RTOL is a double precision scalar. ATOL(*) is either a double
          precision scalar or a double precision array of length NEQN.
          RKC is designed for the solution of problems to modest accuracy.
          Because it is based on a method of order 2, it is relatively
          expensive to achieve high accuracy.

          RTOL is a relative error tolerance.  You must ask for some
          relative accuracy, but you cannot ask for too much for the
          precision available.  Accordingly, it is required that
          0.1 >= RTOL >= 10*uround. (See below for the machine and
          precision dependent quantity uround.)

          ATOL is an absolute error tolerance that can be either a
          scalar or an array.  When it is an array, the tolerances are
          applied to corresponding components of the solution and when
          it is a scalar, it is applied to all components.  A scalar
          tolerance is reasonable only when all solution components are
          scaled to be of comparable size.  A scalar tolerance saves a
          useful amount of storage and is convenient.  Use INFO(*) to
          tell RKC whether ATOL is a scalar or an array.

          The absolute error tolerances ATOL(*) must satisfy ATOL(i) >= 0
          for i = 1,...,NEQN.  ATOL(j)= 0 specifies a pure relative error
          test on component j of the solution, so it is an error if this
          component vanishes in the course of the integration.

          If all is going well, reducing the tolerances by a factor of
          0.1 will reduce the error in the computed solution by a factor
          of roughly 0.2.

INFO(*)   Integer array of length 4 that specifies how the problem
          is to be solved.

INFO(1):  RKC integrates the initial value problem from T to TEND.
          This is done by computing approximate solutions at points
          chosen automatically throughout [T, TEND].  Ordinarily RKC
          returns at each step with an approximate solution. These
          approximations show how y behaves throughout the interval.
          The subroutine RKCINT can be used to obtain answers anywhere
          in the span of a step very inexpensively. This makes it
          possible to obtain answers at specific points in [T, TEND]
          and to obtain many answers very cheaply when attempting to
          locating where some function of the solution has a zero
          (event location).  Sometimes you will be interested only in
          a solution at TEND, so you can suppress the returns at each
          step along the way if you wish.

INFO(1)  = 0 Return after each step on the way to TEND with a
             solution Y(*) at the output value of T.

         = 1 Compute a solution Y(*) at TEND only.

INFO(2):  RKC needs an estimate of the spectral radius of the Jacobian.
          You must provide a function that must be called SPCRAD and
          have the form

            double precision function spcrad(neqn,t,y)
            integer          neqn
            double precision t,y(neqn)

            spcrad = < expression depending on info(2) >

            return
            end

          You can provide a dummy function and let RKC compute the
          estimate. Sometimes it is convenient for you to compute in
          SPCRAD a reasonably close upper bound on the spectral radius,
          using, e.g., Gershgorin's theorem.  This may be faster and/or
          more reliable than having RKC compute one.

INFO(2)  = 0  RKC is to compute the estimate internally.
              Assign any value to SPCRAD.

         = 1  SPCRAD returns an upper bound on the spectral
              radius of the Jacobian of f at (t,y).

INFO(3):  If you know that the Jacobian is constant, you should say so.

INFO(3)  = 0  The Jacobian may not be constant.

         = 1  The Jacobian is constant.

INFO(4):  You must tell RKC whether ATOL is a scalar or an array.

INFO(4)  = 0  ATOL is a double precision scalar.

         = 1  ATOL is a double precision array of length NEQN.

WORK(*):  Work array.  Double precision array of length at least
          8 + 5*NEQN if INFO(2) = 0 and otherwise, 8 + 4*NEQN.

IDID:     Set IDID = 0 to initialize the integration.



RETURNS FROM RKC

T:        The integration has advanced to T.

Y(*):     The solution at T.

IDID:     The value of IDID reports what happened.

                        SUCCESS

  IDID     = 1 T = TEND, so the integration is complete.

           = 2 Took a step to the output value of T.  To continue on
               towards TEND, just call RKC again.   WARNING:  Do not
               alter any argument between calls.

               The last step, HLAST, is returned as WORK(1). RKCINT
               can be used to approximate the solution anywhere in
               [T-HLAST, T] very inexpensively using data in WORK(*).

               The work can be monitored by inspecting data in RKCDID.

                        FAILURE

           = 3 Improper error control: For some j, ATOL(j) = 0
               and Y(j) = 0.

           = 4 Unable to achieve the desired accuracy with the
               precision available.  A severe lack of smoothness in
               the solution y(t) or the function f(t,y) is likely.

           = 5 Invalid input parameters:  NEQN <= 0, RTOL > 0.1,
               RTOL < 10*UROUND, or ATOL(i) < 0 for some i.

           = 6 The method used by RKC to estimate the spectral
               radius of the Jacobian failed to converge.

RKCDID is a labelled common block that communicates statistics
       about the integration process:
       common /rkcdid/   nfe,nsteps,naccpt,nrejct,nfesig,maxm

       The integer counters are:

      NFE      number of evaluations of F used
                 to integrate the initial value problem
      NSTEPS   number of integration steps
      NACCPT   number of accepted steps
      NREJCT   number of rejected steps
      NFESIG   number of evaluations of F used
                 to estimate the spectral radius
      MAXM     maximum number of stages used

      This data can be used to monitor the work and terminate a run
      that proves to be unacceptably expensive.  Also, if MAXM should
      be far beyond 100, the problem is too expensive for RKC and
      alternative methods should be considered.

 CAUTION: MACHINE/PRECISION ISSUES

   UROUND (the machine precision) is the smallest number such that
   1 + UROUND > 1, where 1 is a floating point number in the working
   precision. UROUND is set in a parameter statement in RKC. Its
   value depends on both the precision and the machine used, so it
   must be set appropriately.  UROUND is the only constant in RKC
   that depends on the precision.

   This version of RKC is written in double precision. It can be changed
   to single precision by replacing DOUBLE PRECISION in the declarations
   by REAL and changing the type of the floating point constants set in
   PARAMETER statements from double precision to real.

Authors: B.P. Sommeijer and J.G. Verwer
         Centre for Mathematics and Computer Science (CWI)
         Kruislaan 413
         1098 SJ  Amsterdam
         The Netherlands
         e-mail: bsom@cwi.nl

         L.F. Shampine
         Mathematics Department
         Southern Methodist University
         Dallas, Texas 75275-0156
         USA
         e-mail: lshampin@mail.smu.edu

Details of the methods used and the performance of RKC can be
found in

       B.P. Sommeijer, L.F. Shampine and J.G. Verwer
       RKC: an Explicit Solver for Parabolic PDEs.
            Technical Report MAS-R9715, CWI, Amsterdam, 1997

This source code for RKC and some examples, as well as the
reference solution to the second example can also be obtained
by anonymous ftp from the address ftp://ftp.cwi.nl/pub/bsom/rkc
```
