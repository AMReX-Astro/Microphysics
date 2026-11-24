*********************
Fundamental Constants
*********************

Fundamental constants are defined in ``constants/fundamental_constants.H``.
This file contains physical constants and unit conversion factors, such as
Boltzmann’s constant, $k_B$, and is **automatically generated** by the
Python script ``write_fundamental_constants.py``.

All constant values are obtained from ``scipy.constants`` to ensure consistency
between ``Microphysics`` and ``pynucastro``, the latter is used for
generating the reaction networks used by ``Microphysics``.

There are two exceptions:

1. Constants under the ``Legacy`` namespace represent outdated values
   historically used in networks such as ``aprox13``, ``aprox19``,
   and ``aprox21``. These constants are retained for backward compatibility
   and should **NOT** be used in other places.
2. Newton’s gravitational constant continues to use an older value to preserve
   **hydrostatic equilibrium** in existing initial models used by ``Castro``.

.. note::
   The consistency of physical constants between ``pynucastro`` and
   ``Microphysics`` is particularly important when using inverse rates derived
   from detailed balance. Discrepancies in constant values can affect the
   agreement between mass fractions obtained from direct integration and those
   computed from nuclear statistical equilibrium (NSE).

   This a subtle point, but important when using NSE evolution mode
   (i.e., when compiled with ``USE_NSE_NET=TRUE`` and ``USE_NSE_TABLE=TRUE``).
   See :ref:`ch:nse` for more details.
