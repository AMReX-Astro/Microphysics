**********************
Transport Coefficients
**********************

Thermal Conductivity
====================

The thermal conductivities, $k_\mathrm{th}$, are provided to allow for
modeling of thermal diffusion, for instance in an energy equation:

.. math::

   \frac{\partial (\rho e)}{\partial t} = \nabla \cdot k_\mathrm{th} \nabla T

Thermal conductivities are provided by the ``conductivity/``
directory.  The main interface has the form:

.. code:: c++

   template <typename T>
   AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
   void conductivity (T& state)

where currently, the input type needs to be the full EOS type, ``eos_t``.  The conductivity
is then available as ``eos_t eos_state.conductivity``.

.. important::

   It is assumed that the state is thermodynamically consistent
   before calling the conductivity routine.
   It may be necessary to do an EOS
   call first, to enforce the consistency.

There are several implementations of the conductivity currently:

* ``constant`` :

  This simply sets the conductivity to a constant value set via
  ``conductivity.const_conductivity``.

* ``constant_opacity`` :

  This is intended for creating a conductivity from an opacity, for
  instance, electron scattering.  The opacity is taken as constant
  and set via ``conductivity.const_opacity``, and then the conductivity
  is set as:

  .. math::

     k_\mathrm{th} = \frac{16 \sigma_\mathrm{SB} T^3}{3 \kappa \rho}

  where $\sigma_\mathrm{SB}$ is the Stefan-Boltzmann constant.

* ``powerlaw`` :

  This models a temperature-dependent conductivity of the form:

  .. math::

     k_\mathrm{th} = k_0 T^\nu

  where $k_0$ is set via ``conductivity.cond_coeff`` and $\nu$ is set via
  ``conductivity.cond_exponent``.

* ``stellar`` :

  This is a general stellar conductivity based on :cite:`timmes:2000b`.
  It combines radiative opacities and thermal conductivities into a
  single expression.  This conductivity is suitable for modeling
  laminar flames.

.. index:: CONDUCTIVITY_DIR

These can be set via the ``CONDUCTIVITY_DIR`` make variable.


Radiation Opacities
===================

For radiation transport, we provide simple opacity routines that return
the Planck and Rosseland mean opacities.

The interface for these opacities is:

.. code:: c++

   AMREX_GPU_HOST_DEVICE AMREX_INLINE
   void actual_opacity (amrex::Real& kp, amrex::Real& kr,
                        amrex::Real rho, amrex::Real temp, amrex::Real rhoYe, amrex::Real nu,
                        bool get_Planck_mean, bool get_Rosseland_mean)

where the boolean ``get_Planck_mean`` and ``get_Rosseland_mean`` specify where those
opacity terms are filled upon return.

There are 2 interfaces, which are selected via the make variable ``OPACITY_DIR``.

* ``breakout`` :

  This returns a simple electron scattering opacity with a crude approximation
  connecting the Planck and Rosseland means.

* ``rad_power_law`` :

  This constructs a power-law opacity.  For the Planck mean, it is:

  .. math::

     \kappa_p = \kappa_{p,0} \rho^{m_p} T^{-{n_p}} \nu^{p_p}

  where $\kappa_{p,0}$ is set via ``opacity.const_kappa_p``, $m_p$ is set via ``opacity.kappa_p_exp_m``,
  $n_p$ is set via ``opacity.kappa_p_exp_n``, and $p_p$ is set via ``opacity.kappa_p_exp_p``.

  The Rosseland mean has a scattering component, $\kappa_s$, which is computed as:

  .. math::

     \kappa_s = \kappa_{s,0} \rho^{m_s} T^{-{n_s}} \nu^{p_s}

  where $\kappa_{s,0}$ is set via ``opacity.const_scatter``, $m_s$ is set via ``opacity.scatter_exp_m``,
  $n_s$ is set via ``opacity.scatter_n``, and $p_s$ is set via ``opacity.scatter_p``.  The Rosseland
  mean is then computed as:

  .. math::

     \kappa'_r = \kappa_{r,0} \rho^{m_r} T^{-{n_r}} \nu^{p_r}

  where $\kappa_{r,0}$ is set via ``opacity.const_kappa_r``, $m_r$ is set via ``opacity.kappa_r_exp_m``,
  $n_r$ is set via ``opacity.kappa_r_exp_n``, and $p_r$ is set via ``opacity.kappa_r_exp_p``, and then
  combined with scattering as:

  .. math::

     \kappa_r = \max \{ \kappa'_r + \kappa_s, \kappa_\mathrm{floor} \}

  where $\kappa_\mathrm{floor}$ is set via ``opacity.kappa_floor``,
