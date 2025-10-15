.. _sec:eos_cell:

************
``eos_cell``
************

.. index:: eos_cell

``eos_cell`` simply calls the equation of state on an input density, temperature,
and composition, and then outputs the full thermodynamic state.  This is mainly
used to understand the thermodynamics one might encounter in a simulation
when using a particular EOS.

Getting Started
===============

The ``eos_cell`` code is located in
``Microphysics/unit_test/eos_cell``.  An inputs file which sets the
default parameters for your thermodynamic state is needed to run the
test.

Setting the thermodynamics
--------------------------

The parameters that affect the thermodynamics are:

* ``unit_test.density`` : the initial density

* ``unit_test.temperature`` : the initial temperature

The composition can be set either by specifying individual mass fractions
or setting ``unit_test.uniform_xn`` as described in :ref:`sec:defining_unit_test_composition`.


Building and Running the Code
=============================

The code can be built simply as:

.. prompt:: bash

   make

.. note::

   Even though there are no reactions, a network is still required,
   and can be set via the ``NETWORK_DIR`` build variable.  By default,
   the ``aprox13`` network is used.

   The network choice serves only to set the composition, and a
   ``general_null`` network may also be used.

The build process will automatically create links in the build
directory to any required EOS table.

To run the code, in the ``eos_cell`` directory run::

   ./main3d.gnu.ex inputs_eos

where ``inputs_eos`` is the provided inputs file.  You may edit the
thermodynamic state in that file prior to running.


Output
======

All output is directed to ``stdout`` and simply lists the entries in the
full ``eos_t`` datatype.  Example output is shown below:

.. literalinclude:: ../../unit_test/eos_cell/ci-benchmarks/eos_helmholtz.out
