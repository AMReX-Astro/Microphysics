# ROCK4 Integrator

Here are the original comments from the `rock4.f` source code:

```

   Numerical solution of a (mildly) stiff system
   of first order differential equations. ROCK4 is
   based on a family of 4th order explicit Runge-
   Kutta methods with nearly optimal stability domain
   on the negative real axis. The numerical method is
   based on a three-term recursion relation, and a
   composition of two sub-methods.
   The size (along the negative axis) of the stability
   domains increases quadratically with the stage number.

   Intended for problems of large dimensions with
   eigenvalues of the Jacobian close to the negative
   real axis. Typically for problems originating from
   paraboliPDEs.


   Author: A. Abdulle
           Universite de Geneve, Dept. de mathematiques
           Ch-1211 Geneve 24, Switzerland
           e-mail: assyr.abdulle@math.unige.ch

   Version of April 2002  minor bug corrected

   The analysis of the ROCK4 method is described in:

   [ ] A. Abdulle
       Fourth order Chebyshev methods with
       recurrence relation
       To appear in SISC
       http://www.unige.ch/math/biblio/preprint/liste.html

    Input parameters
    ----------------
    NEQN:       Number of differential equations of the system
                (integer).

    T:          Initial point of integration (double precision).

    TEND:       End of the interval of integration,
                may be less than t (double precision)

    H:          Initial step size guess
                (usually between 1d-4 and 1d-6).

    Y(NEQN):    Initial value of the solution
                (double precision array of length neqn).

    F:          Name (external) of subroutine computing the value
                of f(x,y). Must have the form

                  subroutine f(neqn,t,y,dy)
                  double precision y(neqn),dy(neqn)
                  integer neqn
                  dy(1)=...
                  ...
                  dy(neqn)=...
                  return
                  end

                Implementation:
                for stability issues when the problem
                is originating from paraboliPDEs, transforming
                inhomogeneous boundary conditions in homogeneous ones
                (by adding the appropriate function to the right-hand side)
                may increase the performance of the code.


    ATOL(*) :   Absolute and relative error tolerances
    RTOL(*)     can be both scalar (double precision)
                or vectors of length neqn (double precision).

    RHO:        Name (external) of a function (double precision)
                giving the spectral radius of the Jacobian
                matrix of f at (t,y). Must have the form

                  double precision function rho(neqn,t,y)
                  double precision y(neqn),t
                  integer neqn
                  ...
                  rho=...
                  return
                  end

                N.b. Gerschgorin's theorem can be helpful. If the
                Jacobian is known to be constant it should be
                specified by setting iwork(2)=1 (see below).

                ROCK4 can also compute this estimate. In that
                case, provide a dummy function rho(neqn,t,y)
                and set iwork(1)=0 (see below).

                If it is possible to give an estimate of
                the spectral radius, it should be preferred to
                the estimate computed internally by ROCK4.

    IWORK(*):   Integer array of length 12 that gives information
                on how the problem is to be solved and communicates
                statistics about the integration process.

    IWORK(1):   =0 ROCK4 attempts to compute the spectral radius
                   internally. Define a dummy function

                   double precision function rho(neqn,t,y)
                   double precision y(neqn),t
                   integer neqn
                   rho=0.d0
                   return
                   end

                =1 RHO returns an upper bound of the spectral
                   radius  of the Jacobian matrix of f at (t,y).

    IWORK(2):   =0 The Jacobian is not constant.
                =1 The Jacobian is constant
                   the function rho is called only once.

    IWORK(3):   =0 Return and solution at tend.
                =1 The code returns after each step t_i chosen
                   automatically between [t,tend] (solution
                   at t_i is in y(*) ).
                   To continue call ROCK4 again without changing
                   any arguments.

    IWORK(4):   =0 Atol and rtol are scalar.
                =1 Atol and rtol are array of length neqn.

    WORK(*):      Workspace of length 8*neqn if iwork(1)=0,
                  otherwise of length 7*neqn.
                  Work(1),..,work(7*neqn) serve as
                  working space for the solution of
                  the ode.
                  Work(7*neqn+1),..,work(8*neqn)
                  serve as working space for the
                  internal computation of the
                  spectral radius of the Jacobian.

    IDID:         Report on successfulness upon return
                  (integer).

    Output parameters
    -----------------
    T:          T-value for which the solution has been computed
                (after successful return t=tend).

    Y(NEQN):    Numerical solution at tend.

    IDID:       Reports what happened upon return

    IDID        =1 Successful computation t=tend.
                =2 Successful computation of one step.
                   To continue call ROCK4 again without
                   altering any arguments.
                =-1 invalid input parameters.
                =-2 Stepsize becomes to small.
                =-3 The method used in ROCK4 to estimate
                    the spectral radius did not converge.

    IWORK(5)    =Number of function evaluations.
    IWORK(6)    =Number of steps.
    IWORK(7)    =Number of accepted steps.
    IWORK(8)    =Number of rejected steps.
    IWORK(9)    =Number of evaluations of f used
                to estimate the spectral radius
                (equal to zero if iwork(1)=1).
    IWORK(10)   =Maximum number of stages used.
    IWORK(11)   =Maximum value of the estimated
                 bound for the spectral radius
                 (rounded to the nearest integer).
    IWORK(12)   =Minimum value of the estimated
                 bound for the spectral radius
                 (rounded to the nearest integer).


   Caution:     The variable UROUND (the rounding unit) is set to
   -------      1.0d-16 and may depends on the machines.

*** *** *** *** *** *** *** *** *** *** *** *** *** *** ***
        Numerical method
*** *** *** *** *** *** *** *** *** *** *** *** *** *** ***
     The nearly optimal stability polynomial is computed
     as a product:  R_s(z)=P_{s-4}(z)*w(z).
     We realize this polynomial as a composition of two
     Runge-Kutta methods W(P), where the stability
     polynomial of the method P is P_{s-4}(z) and the stability
     polynomial of the method W is w(z). The first s-4 stages
     possess a three-term recurrence formula.

*** *** *** *** *** *** *** *** *** *** *** *** *** *** ***
        Stability functions and three-term recurrence relation
*** *** *** *** *** *** *** *** *** *** *** *** *** *** ***

    The stability functions: R_j(z)=P_{j-4}(z)  (internal j<=ms-4)
                            R_{ms}=P_{s-4}(z)*w(z) ( absolute)
                            w(z)=1+a_1*z+a_2*z^2+a_3*z^3+a_4*z^4

    P_j(z) orthogonal with respect to w(z)^2/sqrt{1-x^2}

    Recurrence formula:

    P_j(z)=(a_j*z-b_j)*P_{j-1}(z)-c_j*P_{j-2}(z)
                j=1..ms-4 b_1=-1,c_1=0

    Normalization: P_(0)=1  =>b_{j}=-(1+c_{j})

    Runge-Kutta formula:

    g_j(z)=a_{j}*f(g_{j-1})-b_{j}*g_{j-1}-c_{j}*g_{j-2} j=1..ms-4

    Data (rec. param.): rec(i)= a_1,a_2,c_2,a_3,c_3,.,a_(ms-4),
                         c_(ms-4)  for ms=1,3,5,..

    The finishing procedure:

    g_{s-3}=g_{s-4} + h*a_{2,1}*f(g_{s-4}
    g_{s-2}=g_{s-4} + h*(a_{3,1}*f(g_{s-4}+a_{3,2}*f(g_{s-3})
    g_{s-1}=g_{s-4} + h*(a_{4,1}*f(g_{s-4}
                    + a_{4,2}*f(g_{s-3}a_{4,3}*f(g_{s-2})
    y1:=g_{s}=g_{s-4} h*(b_1*f(g_{s-4}+b_2*f(g_{s-3}+
    b_3*f(g_{s-2})+b_4*f(g_{s-1})

    Embedded method:

    y1e:=g_{s}=g_{s-4} h*(be_1*f(g_{s-4}+be_2*f(g_{s-3}+
    be_3*f(g_{s-2})+be_4*f(g_{s-1}+be_5*f(y_1))

    Data (finish. proced.):fpa(i,j)=a_ij,fpb(i)=b_i,fpbe(i)=be_i
                            i=1,.,4 j=1,.,i-1  for ms=1,3,5,..

    Chosen degrees: s=5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,
    26,28,30,32,34,36,38,40,42,45,48,51,54,57,60,63,67,71,75,
    80,85,90,96,102,109,116,124,133,142,152  ms=s-4
```
