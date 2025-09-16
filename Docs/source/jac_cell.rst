.. _sec:jac_cell:

************
``jac_cell``
************

.. index:: jac_cell, numerical_jacobian.H

``jac_cell`` is used to test the accuracy of an analytic Jacobian
provided by a network, but comparing to a finite-difference
approximation.  Given a thermodynamic state, ``jac_cell`` evaluates
the analytic Jacobian as well as a numerical approximation to it, and
then outputs (to ``stdout``) the Jacobian elements side-by-side for
each method of computing the Jacobian.

.. note::

   Some integrators, like VODE, have their own numerical Jacobian
   routines.  This test uses the implementation in
   ``integration/util/numerical_jacobian.H``, which closely follows the ideas
   of the VODE numerical approximation, as described in :cite:`lsode`.


Getting Started
===============

The ``jac_cell`` code is located in
``Microphysics/unit_test/jac_cell``.  An inputs file which sets the
default parameters for your thermodynamic state is needed to run the
test.

Setting the thermodynamics
--------------------------

The parameters that affect the thermodynamics are:

* ``unit_test.density`` : the initial density

* ``unit_test.temperature`` : the initial temperature

The composition can be set either by specifying individual mass fractions
or setting ``unit_test.uniform_xn`` as described in :ref:`sec:defining_unit_test_composition`.

If the values don't sum to ``1`` initially, then the test will do a
normalization.  This normalization can be disabled by setting:

::

    unit_test.skip_initial_normalization = 1


Building and Running the Code
=============================

The code can be built simply as:

.. prompt:: bash

   make

.. note::

   By default, this will build with the ``aprox13`` network, which
   uses the C++ templating method (:ref:`sec:templated_rhs`) of
   building the analytic Jacobian at compile-time.  The network can be
   changed via the build parameter ``NETWORK_DIR``.

.. important::

   The screening routines provided by Microphysics do not return the
   derivative of the screening factor with respect to composition.  This
   means that the analytic Jacobian provided by the network will ignore
   this contribution, but the numerical Jacobian will implicitly include
   it.

   To compare just the network Jacobian elements, it is suggested that
   you build with ``SCREEN_METHOD=null``.


The build process will automatically create links in the build
directory to any required EOS table.

To run the code, in the ``jac_cell`` directory run::

   ./main3d.gnu.ex inputs

where ``inputs`` is the provided inputs file.  You may edit the
thermodynamic state in that file prior to running.


Output
======

All output is directed to ``stdout``.  Example output is shown below
(this output includes screening, highlighting the difference that the
screening composition derivatives make):


.. literalinclude:: ../../unit_test/jac_cell/ci-benchmarks/jac_cell_aprox13.out
