.. _sec:burn_cell_sdc:

*****************
``burn_cell_sdc``
*****************

.. index:: burn_cell_sdc

``burn_cell_sdc`` is a simple one-zone burn analogous to :ref:`sec:burn_cell` that
exercises the simplified-SDC code paths in the integrator.
The system that is evolved
has the form:

.. math::

   \begin{align*}
      \frac{d(\rho X_k)}{dt} &= \Adv{\rho X_k}^{n+1/2} + \rho \dot{\omega}_k(\rho, X_k, T) \\
      \frac{d(\rho e)}{dt} &= \Adv{\rho e}^{n+1/2} + \rho \epsilon(\rho, X_k, T)
   \end{align*}

with density constructed as needed via:

$$\rho(t) = \rho^n + \Adv{\rho}^{n+1/2} (t - t^n)$$

In this system, we now need to specify the advective terms for the unit test.

.. note::

   This test can also be used with NSE to test how the integration
   bails out into NSE if the ``in_nse()`` test is true.




Getting Started
===============

The ``burn_cell_sdc`` code is located in
``Microphysics/unit_test/burn_cell_sdc``.   By default, ``USE_SIMPLIFIED_SDC=TRUE``
is set in the ``GNUmakefile``.


Setting the thermodynamics
--------------------------

The parameters that affect the thermodynamics are:

* ``unit_test.density`` : the initial density

* ``unit_test.temperature`` : the initial temperature

* ``unit_test.rhoe`` : the initial $(rho e)$.  If this is not set (or
  set to be $< 0$), then it will be computed from the temperature
  using the EOS.

* ``unit_test.small_temp`` : the low temperature cutoff used in the equation of state

* ``unit_test.small_dens`` : the low density cutoff used in the equation of state

As with ``burn_cell``, the composition can be set either by setting each mass fraction explicitly via the
parameters, ``unit_test.X1``, ``unit_test.X2``, ..., or initialized to be
uniform via ``unit_test.uniform_xn=1``.

Aux composition
---------------

Advective terms
---------------

But default, the advective terms are set to zero.  In this mode, ``burn_cell_sdc``
largely mimics ``burn_cell`` (although with a slightly different ODE system
integrated).

The advective terms can be set




Controlling time
----------------

The test will run unit a time ``unit_test.tmax``, outputting the state
at regular intervals.  The parameters controlling the output are:

* ``unit_test.tmax`` : the end point of integration.

* ``unit_test.tfirst`` : the first time interval to output.

* ``unit_test.nsteps`` : the number of steps to divide the integration into,
  logarithmically-spaced.

If there is only a single step, ``unit_test.nsteps = 1``, then we integrate
from $[0, \mathrm{tmax}]$.

If there are multiple steps, then the first output will be at a time
$\mathrm{tmax} / \mathrm{nsteps}$, and the steps will be
logarithmically-spaced afterwards.


Integration parameters
----------------------

The tolerances, choice of Jacobian, and other integration parameters
can be set via the usual Microphysics runtime parameters, e.g.
``integrator.atol_spec``.


Rerunning a burn fail
---------------------


Building and Running the Code
=============================

The code can be built simply as:

.. prompt:: bash

   make

and the network and integrator can be changed using the normal
Microphysics build system parameters.

.. important::

   You need to do a ``make clean`` before rebuilding with a different
   network or integrator.


To run the code, in the ``burn_cell_sdc`` directory run::

   ./main3d.gnu.ex inputs

where ``inputs`` is the name of your inputs file.
