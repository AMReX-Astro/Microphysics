******************
Equations of State
******************

In this chapter on equations of state, we list the EOS routines
available for your use, and then we describe the correct structure of
an EOS module in case you want to build your own.

Available Equations of State
============================

The following equations of state are available in Microphysics.
Except where noted, each of these EOSs will provide the full
thermodynamic data (including all derivatives) in the ``eos_t``
type.

gamma_law
---------

``gamma_law`` represents a gamma law gas, with
equation of state:

.. math:: p = (\gamma - 1) \rho e.

:math:`\gamma` is specified by the runtime parameter ``eos_gamma``. For
an ideal gas, this represents the ratio of specific heats. The gas is
assumed to be ideal, with the pressure given by

.. math:: p = \frac{\rho k T}{\mu m_u}

where :math:`k` is Boltzmann’s constant and :math:`\mu` is the mean molecular
weight, calculated from the composition, :math:`X_k`. This EOS assumes
the gas is either completely neutral (``assume_neutral = T``),
giving:

.. math:: \mu^{-1} = \sum_k \frac{X_k}{A_k}

or completely ionized (``assume_neutral = F``), giving:

.. math:: \mu^{-1} = \sum_k \left ( 1 + Z_k \right ) \frac{X_k}{A_k}

The entropy comes from the Sackur-Tetrode equation. Because of the
complex way that composition enters into the entropy, the entropy
formulation here is only correct for a :math:`\gamma = 5/3` gas.

Note that the implementation provided in Microphysics is the same as
the version shipped with MAESTRO, but more general than the
``gamma_law`` EOS provided with CASTRO. CASTRO’s default EOS only
fills the thermodynamic information in ``eos_t`` that is required
by the hydrodynamics module in CASTRO.

polytrope
---------

``polytrope`` represents a polytropic fluid, with equation of
state:

.. math:: p = K \rho^\gamma.

The gas is also assumed to obey the above gamma law relation
connecting the pressure and internal energy. Therefore :math:`\rho` is the
only independent variable; there is no temperature dependence. The
user either selects from a set of predefined options reflecting
physical polytropes (e.g. a non-relativistic, fully degenerate
electron gas) or inputs their own values for :math:`K` and :math:`\gamma`
via ``polytrope_K`` and ``polytrope_gamma``.

The runtime parameter ``polytrope_type`` selects the pre-defined
polytropic relations. The options are:

-  ``polytrope_type = 1``: sets :math:`\gamma = 5/3` and

   .. math:: K = \left ( \frac{3}{\pi} \right)^{2/3} \frac{h^2}{20 m_e m_p^{5/3}} \frac{1}{\mu_e^{5/3}}

   where :math:`mu_e` is the mean molecular weight per electron, specified via ``polytrope_mu_e``

   This is the form appropriate for a non-relativistic
   fully-degenerate electron gas.

-  ``polytrope_type = 2``: sets :math:`\gamma = 4/3` and

   .. math:: K = \left ( \frac{3}{\pi} \right)^{1/3} \frac{hc}{8 m_p^{4/3}} \frac{1}{\mu_e^{4/3}}

   This is the form appropriate for a relativistic fully-degenerate
   electron gas.

ztwd
----

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

multigamma
----------

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

-  ``eos_gamma_default``: the default :math:`\gamma` to apply for all
   species

-  ``species_X_name`` and ``species_X_gamma``: set the
   :math:`\gamma_i` for the species whose name is given as
   ``species_X_name`` to the value provided by ``species_X_gamma``.
   Here, ``X`` can be one of the letters: ``a``, ``b``, or
   ``c``, allowing us to specify custom :math:`\gamma_i` for up to three
   different species.

helmholtz
---------

``helmholtz`` contains a general, publicly available stellar
equation of state based on the Helmholtz free energy, with
contributions from ions, radiation, and electron degeneracy, as
described in :cite:`timmes:1999`, :cite:`timmes:2000`, :cite:`flash`.

We have modified this EOS a bit to fit within the context of our
codes. The vectorization is explicitly thread-safe for use with OpenMP
and OpenACC. In addition, we have added the ability to perform a
Newton-Raphson iteration so that if we call the EOS with density and
energy (say), then we will iterate over temperature until we find the
temperature that matches this density–energy combination. If we
cannot find an appropriate temperature, we will reset it to
``small_temp``, which needs to be set in the equation of state wrapper
module in the code calling this. However, there is a choice of whether
to update the energy to match this temperature, respecting
thermodynamic consistency, or to leave the energy alone, respecting
energy conservation. This is controlled through the
``eos.eos_input_is_constant`` parameter in your inputs file.

We thank Frank Timmes for permitting us to modify his code and
publicly release it in this repository.

stellarcollapse
---------------

``stellarcollapse`` is the equation of state module provided
on http://stellarcollapse.org. It is designed
to be used for core-collapse supernovae and is compatible with a
large number of equations of state that are designed to describe
matter near nuclear density. You will need to download an
appropriate interpolation table from that site to use this.

Interface and Modes
===================

The EOS is called as:

.. code:: c++

   eos(mode, eos_type)

where *mode* determines which thermodynamic quantities are inputs,
and is one of:

* ``eos_input_rt`` : density and temperature are inputs

* ``eos_input_rh`` : density and specific enthalpy are inputs

* ``eos_input_tp`` : temperature and pressure are inputs

* ``eos_input_rp`` : density and pressure are inputs

* ``eos_input_re`` : density and specific internal energy are inputs

* ``eos_input_ps`` : pressure and entropy are inputs

* ``eos_input_ph`` : pressure and specific enthalpy are inputs

* ``eos_input_th`` : temperature and specific enthalpy are inputs

The *eos_type* passed in is one of

* ``eos_t`` : provides access to all available thermodynamic information,
  including derivatives.

* ``eos_re_t`` : only provides the energy-based thermodynamic information, including
  energy derivatives.

* ``eos_rep_t`` : expands on ``eos_re_t`` to include pressure information

In general, you should use the type that has the smallest set of
information needed, since we optimize out needless quantities at
compile type (via C++ templating) for ``eos_re_t`` and ``eos_rep_t``.

.. note::

   All of these modes require composition as an input.  Usually this is
   via the set of mass fractions, ``eos_t.xn[]``, but if ``USE_AUX_THERMO``
   is set to ``TRUE``, then we instead use the auxiliary quantities
   stored in ``eos_t.aux[]``.

.. _aux_eos_comp:

Auxiliary Composition
---------------------


With ``USE_AUX_THERMO=TRUE``, we interpret the composition from the auxiliary variables.
The auxiliary variables are

* ``eos_state.aux[iye]`` : electron fraction, defined as

  .. math::

     Y_e = \sum_k \frac{X_k Z_k}{A_k}

* ``eos_state.aux[iabar]`` : the average mass of the nuclei, :math:`\bar{A}`, defined as:

  .. math::

     \frac{1}{\bar{A}} = \sum_k \frac{X_k}{A_k}

  In many stellar evolutions texts, this would be written as :math:`\mu_I`.

* ``eos_state.aux[ibea]`` : the binding energy per nucleon (units of
  MeV), defined as

  .. math::

     \left \langle \frac{B}{A} \right \rangle  = \sum_k \frac{X_k B_k}{A_k}

  where :math:`B_k` is the binding energy of nucleus :math:`k`

Given a composition of mass fractions, the function
``set_aux_comp_from_X(state_t& state)`` will initialize these
auxiliary quantities.

The equation of state also needs :math:`\bar{Z}` which is easily computed as

.. math::

   \bar{Z} = \bar{A} Y_e


Initialization and Cutoff Values
================================

Input Validation
================

The EOS will make sure that the inputs are within an acceptable range,
(e.g., ``small_temp`` :math:`< T <` ``maxT``). If they are not, then it
resets them silently—no error is thrown.

If you are calling the EOS with ``eos_input_re``, and if :math:`e <
10^{-200}`, then it calls the EOS with ``eos_input_rt`` with T =
max ( ``small_temp``, T ).

User’s are encourage to do their own validation of inputs before calling
the EOS.

EOS Structure
=============

Each EOS should have two main routines by which it interfaces to the
rest of CASTRO. At the beginning of the simulation,
``actual_eos_init`` will perform any initialization steps and save
EOS variables (mainly ``smallt``, the temperature floor, and
``smalld``, the density floor). These variables are stored in the
main EOS module of the code calling these routines. This would be the
appropriate time for, say, loading an interpolation table into memory.

The main evaluation routine is called ``actual_eos``. It should
accept an eos_input and an eos_t state; see Section :ref:`data_structures`.



