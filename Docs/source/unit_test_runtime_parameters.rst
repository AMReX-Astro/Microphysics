***********************************
Unit Test Common Runtime Parameters
***********************************

There are a number of runtime parameters that are common to all (or most) of the unit tests.
These are defined in the top-level ``unit_test/_parameters`` file.

Thermodynamics
==============

The equation of state enforces minimum density and temperatures, which must be set
upon initialization.  These are controlled by the following runtime parameters:

* ``unit_test.small_temp`` : the low temperature cutoff used in the equation of state

* ``unit_test.small_dens`` : the low density cutoff used in the equation of state


.. _sec:defining_unit_test_composition:

Defining composition
====================

Most of the unit tests require a composition to be defined (for the
initial mass-fractions, $X_k$).  There are a few ways this can be done
(depending on the test).


* One-zone (``*_cell``) tests (see :ref:`sec:one_zone_tests`) usually do one of:

  * *Explicitly setting the individual mass fractions.*  This is
    controlled by the parameters ``unit_test.X1``, ``unit_test.X2``, ..., ``unit_test.X35``,
    e.g.:

    ::

        unit_test.X1 = 0.5
        unit_test.X2 = 0.2
        unit_test.X3 = 0.2
        unit_test.X4 = 0.1

    While many of the tests will renormalize the abundances, the user
    should take care to ensure that the mass fractions sum to unity.

  * *Setting the composition to be uniform.*  This is controlled by
    ``unit_test.uniform_xn``.  If this is set to ``1``, then each mass fraction
    is initialized to ``1 / NumSpec``.

* Comprehensive tests (see :ref:`sec:comprehensive_tests`) need many different compositions, since they are creating a cube

  of varying thermodynamic properties, and thus require a prescription
  to create the composition.  This is done by setting ``unit_test.primary_species_1``,
  ``unit_test.primary_species_2``, and ``unit_test.primary_species_3`` to one of the
  *names* of the species in the network.

  The function ``setup_composition()`` is then used to set limits on
  the species abundances (it takes a parameter which is the index into
  the cube of data that is being initialized) which is then used by
  ``get_xn()`` to create the individual mass fractions.  Both of these
  routines are contained in ``react_util.H``.
