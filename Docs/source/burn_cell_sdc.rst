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

The composition can be set either by specifying individual mass fractions
or setting ``unit_test.uniform_xn`` as described in :ref:`sec:defining_unit_test_composition`.

Aux composition
---------------

When built with ``USE_NSE_TABLE=TRUE`` (see :ref:`tabulated_nse`) or with
``USE_AUX_THERMO=TRUE`` (see :ref:`aux_eos_comp`) then the auxiliary
composition, $Y_e$, $\bar{A}$, $\langle
B/A\rangle$, is defined.

The auxiliary composition can either be initialized directly via
``unit_test.Aux1``, ``unit_test.Aux2``, and ``unit_test.Aux3``, or
their values can be computed at initialization from the mass fractions
if ``unit_test.recompute_aux=1`` is set.



Advective terms
---------------

But default, the advective terms are set to zero.  In this mode,
``burn_cell_sdc`` largely mimics ``burn_cell`` (although with a
slightly different ODE system integrated).

The advective terms can be set manually via

* ``unit_test.Adv_rho`` : $\Adv{\rho}^{n+1/2}$
* ``unit_test.Adv_rhoe`` : $\Adv{\rho e}^{n+1/2}$
* ``unit_test.Adv_X1``, ``unit_test.Adv_X2``, ... : $\Adv{\rho X_1}^{n+1/2}$, $\Adv{\rho X_2}^{n+1/2}$, ... (up to $X_{35}$)
* ``unit_test.Adv_Aux1``, ``unit_test.Adv_Aux2``, ``unit_test.Adv_Aux3`` : $\Adv{\rho \alpha_1}^{n+1/2}$, $\Adv{\rho \alpha_2}^{n+1/2}$, $\Adv{\rho \alpha_3}^{n+1/2}$




Controlling time
----------------

The integration time and output frequency are controlled by the
same set of parameters as in ``burn_cell``. see :ref:`sec:burn_cell_time`.


Integration parameters
----------------------

The tolerances, choice of Jacobian, and other integration parameters
can be set via the usual Microphysics runtime parameters, e.g.
``integrator.atol_spec``.


Rerunning a burn fail
---------------------

.. index:: parse_integration_failure.py, USE_GPU_PRINTF

When a network integration encounters a failure, it will output the
entire burn state to ``stdout`` (for GPU builds, this needs to be
enabled explicitly by building with ``USE_GPU_PRINTF``).

The script ``unit_test/burn_cell_sdc/parse_integration_failure.py``
can be used to parse the output (copy and paste the full error into a
file) and produce the runtime parameter settings needed to reproduce
the burn.  This is especially important with SDC, since it will contain
all of the advective terms.


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
