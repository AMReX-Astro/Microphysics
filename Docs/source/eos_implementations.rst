****************************
Available Equations of State
****************************

.. index:: eos_t, EOS_DIR

The following equations of state are available in Microphysics.

.. note::

   The EOS is chosen at compile-time via the ``EOS_DIR`` make
   variable.

.. note::

   Except where noted, each of these EOSs will provide the full
   thermodynamic data (including all derivatives) in the ``eos_t``
   type.


``breakout``
============

The ``breakout`` EOS is essentially the same as ``gamma_law``, but it gets
its composition information from the auxiliary data.  In particular,
it expects an auxiliary quantity named ``invmu`` which is the inverse
of the mean molecular weight:

.. math::

   \frac{1}{\mu} = \sum_k \frac{Z_k X_k}{A_k}

The ``general_null`` network provides this when used with the ``breakout.net``
network inputs.

``gamma_law``
=============

.. index:: eos.eos_gamma, eos.assume_neutral

``gamma_law`` represents a gamma law gas, with
equation of state:

.. math:: p = (\gamma - 1) \rho e.

:math:`\gamma` is specified by the runtime parameter ``eos.eos_gamma``. For
an ideal gas, this represents the ratio of specific heats. The gas is
assumed to be ideal, with the pressure given by

.. math:: p = \frac{\rho k T}{\mu m_u}

where :math:`k` is Boltzmann’s constant and :math:`\mu` is the mean molecular
weight, calculated from the composition, :math:`X_k`. This EOS assumes
the gas is either completely neutral (``eos.assume_neutral = 1``),
giving:

.. math:: \mu^{-1} = \sum_k \frac{X_k}{A_k}

or completely ionized (``eos.assume_neutral = 0``), giving:

.. math:: \mu^{-1} = \sum_k \left ( 1 + Z_k \right ) \frac{X_k}{A_k}

The entropy comes from the Sackur-Tetrode equation. Because of the
complex way that composition enters into the entropy, the entropy
formulation here is only correct for a :math:`\gamma = 5/3` gas.


``helmholtz``
=============

``helmholtz`` contains a general, publicly available stellar
equation of state based on the Helmholtz free energy, with
contributions from ions, radiation, and electron degeneracy, as
described in :cite:`timmes:1999, timmes:2000, flash`.

.. note::

   Our implementation of the ``helmholtz`` EOS has been modified
   extensively from the original Fortran source.  It has been
   made threadsafe and makes heavy use of C++ templating to optimize
   the evaluation of thermodynamic quantities.

The ``helmholtz`` EOS has the ability to perform a Newton-Raphson
iteration so that if we call the EOS with, e.g., density and energy,
and iterate over temperature until we find the temperature
that matches this density–energy combination. If we cannot find an
appropriate temperature, we will reset it to ``small_temp``, which
needs to be set in the equation of state wrapper module in the code
calling this.

.. index:: eos.use_eos_coulomb, eos.eos_input_is_constant, eos.eos_ttol, eos.eos_dtol, eos.prad_limiter_rho_c, eos.prad_limiter_delta_rho

The following runtime parameters affect the EOS:

* ``eos.use_eos_coulomb`` : do we include Coulomb corrections?  This
  is enabled by default.  Coulomb corrections can cause problems in
  some regimes, because the implementation in ``helmholtz`` doesn't
  have the correct asymptotic behavior and can lead to negative
  pressures or energies.

* ``eos.eos_input_is_constant`` : when inverting the EOS for find the
  density and/or temperature that match the inputs, there is a choice
  of whether to update the inputs to match the final density /
  temperature, respecting thermodynamic consistency.  If
  ``eos_input_is_constant=1`` is set (the default), then we leave the
  input thermodynamic quantities unchanged, respecting energy
  conservation.

* ``eos.eos_ttol``, ``eos.eos_dtol`` : these are the tolerances
  for temperature and density used by the Newton solver when
  inverting the EOS.

* ``eos.prad_limiter_rho_c``, ``eos.prad_limiter_delta_rho`` : by
  default, radiation pressure is included in the optically-thick, LTE
  limit (with $p_\gamma = (1/3)a T^4$).  At low densities, this can
  cause issues, leading to an artificially high soundspeed dominated
  by radiation when, in fact, we should be optically thin.  These
  parameters allow us turn off the radiation component smoothly,
  starting at a density ``eos.prad_limiter_rho_c`` and transitioning
  via a $\tanh$ profile to zero over a scale
  ``eos.prad_limiter_delta_rho``.

We thank Frank Timmes for permitting us to modify his code and
publicly release it in this repository.

``metal_chem``
==============

This is a multi-gamma equation of state for metal ISM chemistry.

``multigamma``
==============

``multigamma`` is an ideal gas equation of state where each
species can have a different value of :math:`\gamma`. This mainly affects
how the internal energy is constructed as each species, represented
with a mass fraction :math:`X_k` will have its contribution to the total
specific internal energy take the form of :math:`e = p/\rho/(\gamma_k -  1)`.
The main thermodynamic quantities take the form:

.. math::

   \begin{aligned}
   p &= \frac{\rho k T}{m_u} \sum_k \frac{X_k}{A_k} \\
   e &= \frac{k T}{m_u} \sum_k \frac{1}{\gamma_k - 1} \frac{X_k}{A_k} \\
   h &= \frac{k T}{m_u} \sum_k \frac{\gamma_k}{\gamma_k - 1} \frac{X_k}{A_k}\end{aligned}

We recognize that the usual astrophysical :math:`\bar{A}^{-1} = \sum_k
X_k/A_k`, but now we have two other sums that involve different
:math:`\gamma_k` weightings.

The specific heats are constructed as usual,

.. math::

   \begin{aligned}
   c_v &= \left . \frac{\partial e}{\partial T} \right |_\rho =
       \frac{k}{m_u} \sum_k \frac{1}{\gamma_k - 1} \frac{X_k}{A_k} \\
   c_p &= \left . \frac{\partial h}{\partial T} \right |_p =
       \frac{k}{m_u} \sum_k \frac{\gamma_k}{\gamma_k - 1} \frac{X_k}{A_k}\end{aligned}

and it can be seen that the specific gas constant, :math:`R \equiv c_p -
c_v` is independent of the :math:`\gamma_i`, and is simply :math:`R =
k/m_u\bar{A}` giving the usual relation that :math:`p = R\rho T`.
Furthermore, we can show

.. math::

   \Gamma_1 \equiv \left . \frac{\partial \log p}{\partial \log \rho} \right |_s =
      \left ( \sum_k \frac{\gamma_k}{\gamma_k - 1} \frac{X_k}{A_k} \right ) \bigg /
      \left ( \sum_k \frac{1}{\gamma_k - 1} \frac{X_k}{A_k} \right ) =
   \frac{c_p}{c_v} \equiv \gamma_\mathrm{effective}

and :math:`p = \rho e (\gamma_\mathrm{effective} - 1)`.

This equation of state takes several runtime parameters that can set
the :math:`\gamma_i` for a specific species. The parameters are:

.. index:: eos.eos_gamma_default

-  ``eos.eos_gamma_default``: the default :math:`\gamma` to apply for all
   species

-  ``eos.species_X_name`` and ``eos.species_X_gamma``: set the
   :math:`\gamma_i` for the species whose name is given as
   ``eos.species_X_name`` to the value provided by ``eos.species_X_gamma``.
   Here, ``X`` can be one of the letters: ``a``, ``b``, or
   ``c``, allowing us to specify custom :math:`\gamma_i` for up to three
   different species.



``polytrope``
=============

.. index:: eos.polytrope_K, eos.polytrope_gamma, eos.polytrope_type, eos.polytrope_mu_e

``polytrope`` represents a polytropic fluid, with equation of
state:

.. math:: p = K \rho^\gamma.

The gas is also assumed to obey the above gamma law relation
connecting the pressure and internal energy. Therefore :math:`\rho` is the
only independent variable; there is no temperature dependence. The
user either selects from a set of predefined options reflecting
physical polytropes (e.g. a non-relativistic, fully degenerate
electron gas) or inputs their own values for :math:`K` and :math:`\gamma`
via ``eos.polytrope_K`` and ``eos.polytrope_gamma``.

The runtime parameter ``eos.polytrope_type`` selects the pre-defined
polytropic relations. The options are:

-  ``eos.polytrope_type = 1``: sets :math:`\gamma = 5/3` and

   .. math:: K = \left ( \frac{3}{\pi} \right)^{2/3} \frac{h^2}{20 m_e m_p^{5/3}} \frac{1}{\mu_e^{5/3}}

   where :math:`mu_e` is the mean molecular weight per electron, specified via ``eos.polytrope_mu_e``

   This is the form appropriate for a non-relativistic
   fully-degenerate electron gas.

-  ``eos.polytrope_type = 2``: sets :math:`\gamma = 4/3` and

   .. math:: K = \left ( \frac{3}{\pi} \right)^{1/3} \frac{hc}{8 m_p^{4/3}} \frac{1}{\mu_e^{4/3}}

   This is the form appropriate for a relativistic fully-degenerate
   electron gas.


``primordial_chem``
===================

This is a version of the multi-gamma equation of state that models primordial chemistry.

``rad_power_law``
=================

This is an artificial equation of state for radiation transport test problems.  It uses
a parameterization of the specific heat at constant volume:

.. math::

   c_v = A \rho^m T^{-n}

and energy:

.. math::

   e = \frac{A}{1 - n} \rho^m T^{1-n}

where the runtime parameters provide the constants:

* ``eos.eos_const_c_v`` $= A$

* ``eos.eos_c_v_exp_m`` $= m$

* ``eos.eos_c_v_exp_n`` $= n$


``tillotson``
=============

This is an equation of state for hypervelocity impacts based on :cite:`tillotson:1962`.


``ztwd``
========

``ztwd`` is the zero-temperature degenerate electron equation
of state of Chandrasekhar (1935), which is designed to describe
white dward material. The pressure satisfies the equation:

.. math:: p(x) = A \left( x(2x^2-3)(x^2 + 1)^{1/2} + 3\, \text{sinh}^{-1}(x) \right),

with :math:`A = \pi m_e^4 c^5 / (3 h^3)`. Here :math:`x` is a dimensionless
measure of the Fermi momentum, with :math:`\rho = B x^3` and :math:`B = 8\pi \mu_e
m_p m_e^3 c^3 / (3h^3)`, where :math:`\mu_e` is the mean molecular weight
per electron and :math:`h` is the Planck constant.

The enthalpy was worked out by Hachisu (1986):

.. math:: h(x) = \frac{8A}{B}\left(x^2 + 1\right)^{1/2}.

(note the unfortunate notation here, but this :math:`h` is specific
enthalpy). The specific internal energy satisfies the standard
relationship to the specific enthalpy:

.. math:: e = h - p / \rho.

Since the pressure-density relationship does not admit a closed-form
solution for the density in terms of the pressure, if we call the EOS
with pressure as a primary input then we do Newton-Raphson iteration
to find the density that matches this pressure.
