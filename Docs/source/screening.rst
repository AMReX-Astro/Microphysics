.. _sec:screening:

***************************
Screening of Reaction Rates
***************************

.. index:: SCREEN_METHOD

Introduction
------------

Plasma screening is the enhancement of nuclear reaction rates, :math:`R`
due to the Coulomb coupling of the surrounding plasma and electrons
and ions. The enhancement of the reaction rates is done via the form:

.. math::
   R_{\mathrm{scr}} = R \exp{(h)}

where :math:`R_{\mathrm{scr}}` is the screened reaction rate and :math:`h`
characterizes the magnitude of the screening.
Plasma screening can be broken up to different regimes depending
on the Coulomb coupling parameter, :math:`\Gamma`,
of the reaction rate reactants. Generally, :math:`\Gamma \ll 1` and
:math:`\Gamma \gtrsim 1` correspond to weak and strong screening regimes,
respectively. :math:`\Gamma` is defined as:

.. math::
   \Gamma = \alpha(Z_1, Z_2) \Gamma_e

where :math:`\alpha(Z_1, Z_2)` characterizes the Coulomb strength of
the reactants, its definition can vary slightly depending on the
screening routine. :math:`\Gamma_e` is the Coulomb coupling parameter
depending only on the thermodynamic conditions.

.. math::
   \Gamma_e= e^2 \frac{\sqrt[3]{4 \pi n_e / 3}}{k_B T}

where :math:`e` is the electron charge and :math:`n_e` is the
electron number density.


Screening Options
-----------------

The screening enhancement factor can be can be computed using
several different methods, controlled by the make parameter ``SCREEN_METHOD``.
For example,

.. prompt:: bash

   make SCREEN_METHOD=screen5

Any of the available screening methods can be used with any reaction network.

The options are:

* ``screen5`` :

  This is the screening routine from the Kepler stellar evolution code
  and is the default used with the distributed versions of the "aprox"
  family of reaction networks. In the weak screenng regime,
  :math:`\Gamma < 0.3`, it uses screening described in
  :cite:`graboske:1973`. In the strong screening regime,
  :math:`\Gamma > 0.8`, it uses screening described in
  :cite:`jancovici:1977, alastuey:1978, itoh:1979`.
  For the intermediate screening regime, :math:`0.3 < \Gamma < 0.8`,
  a weighted blending between the weak and strong screening are used.
  The overall procedure is described in :cite:`Wallace:1982`.

  This is the default screening method.

* ``chugunov2007`` :

  This implements the screening of :cite:`chugunov:2007`, following
  :cite:`yakovlev:2006` to extend to binary mixtures. It is suitable
  for :math:`\Gamma \lesssim 600`.

* ``chugunov2009`` :

  This implements the screening of :cite:`chugunov:2009`.  The main
  difference is that the 2007 one calculates an effective coupling
  parameter based on the two ions and then treats it as a
  one-component plasma, while this version (2009) treats it fully as a
  multi-component plasma (which is significantly more expensive)

  This includes the portion in the appendix that blends in the weak
  screening limit.

* ``chabrier1998`` :

  This implements the screening of :cite:`Chabrier_1998` as well as
  the quantum corrections for strong screening according to screen5,
  which is suggested in the appendix of :cite:`Calder_2007`.
  This screening is compatible with NSE calculations unlike ``screen5``,
  ``chugunov2007``, and ``chugunov2009``. This screening is valid in the
  weak screening regime, :math:`\Gamma < 0.1`, and strong screening regime,
  :math:`1 \lesssim \Gamma \lesssim 160`.

.. index:: screening.enable_debye_huckel_skip, screening.debye_huckel_skip_threshold

* ``debye_huckel`` :

  This is just the Debye-Hückel weak-screening limit from
  :cite:`chugunov:2009`.

  While it can be used on its own (by building with
  ``SCREEN_METHOD=debye_huckel``, it is really meant to be used as a
  test to determine whether a more extensive screening approximation
  should be used.  By setting ``screening.enable_debye_huckel_skip``,
  we first compute this weak-screening approximation and then, if it
  is larger than ``screening.debye_huckel_skip_threshold``, the full
  screening factor is computed (using the method specified via
  ``SCREEN_METHOD``).

* ``null`` :

  This disables screening by always returning 1 for the screening
  enhancement factor.

Runtime Options
---------------

.. index:: screening.enable_chabrier1998_quantum_corr

* ``screening.enable_chabrier1998_quantum_corr = 1`` in the input file
  enables an additional quantum correction term added to the screening
  factor when ``SCREEN_METHOD=chabrier1998``. This is disabled by
  default since ``chabrier1998`` is often used along with
  ``USE_NSE_NET=TRUE``, and the NSE solver doesn't include quantum
  corrections.
