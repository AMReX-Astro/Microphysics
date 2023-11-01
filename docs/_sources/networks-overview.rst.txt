*****************************
Overview of Reaction Networks
*****************************

Types of Networks
=================

There are several different types of reaction networks in Microphysics.  Largely these can
be split into:

* *hardcoded networks*: the structure of the network rates linking the
  nuclei was written by hand. Many of these networks are *approximate
  networks* where the links between two nuclei includes steps
  involving nuclei that are not directly represented in the network.
  Since the rates are hardcoded, it can be difficult to update the
  individual rates.

  A common example is to combine :math:`(\alpha,\gamma)` and
  :math:`(\alpha,p)(p,\gamma)` channels without explicitly including
  the protons or nuclei produced by the :math:`(\alpha,p)` reaction in
  the network.  This reduces the size of the network at the expense of
  making the individual rates more complex.

  Many of these networks came from the collection of networks that
  `Frank Timmes provides
  <https://cococubed.com/code_pages/burn.shtml>`_, and subsequently
  converted into C++ and adapted to our framework.

* *automatically-generated networks*: these networks are automatically
  generated simply by providing a list of nuclei and having
  `pynucastro <https://github.com/pynucastro/pynucastro>`_ find all of
  the linking rates.  These can easily be regenerated to reflect new rates.

Network Requirements and Structure
==================================

A network both defines the composition advected by the hydro code as
well as describes the burning processes between those isotopes.
Evolving the species in a network requires an integrator. The design
of Microphysics decouples the integrator from the network, allowing
for the ability to swap integrators as desired. We discuss the
integrators in a later section.

A network is defined by a ``.net`` file which provides a list of species
and some data about each species (its name and some isotopic data).

An example is ``Microphysics/networks/iso7/iso7.net``:

.. literalinclude:: ../../networks/iso7/iso7.net

Lines beginning with ``#`` are comments.  In this example, there are 7
nuclei that are evolved in the network.  We also have 3 additional
nuclei prefixed with ``__extra_`` -- these nuclei will not be know to
the ODE system but we will be able to access their properties (like
:math:`A` and :math:`Z`) while constructing the righthand side of the
network.

At build
time, a file ``network_properties.H`` is automatically generated which contains
a number of variables, including:

* ``NumSpec`` : the number of species evolved in the network

* ``NumSpecExtra`` : the number of evolved and "extra" species, as
  prefixed with ``__extra_`` in the ``.net`` file.


* ``NumAux`` : the number of auxiliary quantities needed by the network (these are not evolved).

* ``aion[NumSpec]`` : the atomic weight (in atomic mass units) of the species.  This includes the extra species.

* ``zion[NumSpec]`` : the atomic number of the species.  This includes the extra species.

* ``spec_names[NumSpec]`` : a descriptive name of the species (e.g. "hydrogen-1").  This includes the extra species.

* ``short_spec_names[NumSpec]`` : a shortened version of the species name (e.g. "H1").  This includes the extra species.

* ``aux_names[NumAux]``: the names of the auxiliary quantities

* ``short_aux_names[NumAux]`` : a shortened version of the auxiliary name

.. note::

   A convention adopted in Microphysics is that each network is
   responsible for determining the energy release from a change in
   composition. Most networks will provide an array of the species
   binding energies and a routine to compute the energy yield from the
   reaction rates.

There are two primary files within each network directory.

* ``actual_network.H`` (as well as ``actual_network_data.H`` and ``actual_network_data.cpp``):

   This header defines data for the network such as enumerators for the reaction rates
   and the binding energies. The header must define the following function, which can be
   used to initialize any of the network's internal data at runtime.

   * ``actual_network_init()``

* ``actual_rhs.H`` (as well as ``actual_rhs_data.H`` and ``actual_rhs_data.cpp``):

   This header defines the functions which are used in the integrator for a burn:

   * ``actual_rhs_init()``

   * ``actual_rhs(state, rhs)``

   * ``actual_jac(state, jac)``

   This supplies an interface for computing the right-hand-side of the
   network, the time-derivative of each species (and the temperature
   and nuclear energy release), as well as the analytic Jacobian.
   Both ``actual_rhs`` and ``actual_jac`` take as arguments a burn_t
   state and (respectively) the time-derivatives and Jacobian
   elements to fill in.

   Note: some networks do not provide an analytic Jacobian and instead
   rely on the numerical difference-approximation to the Jacobian. In
   this case, the interface ``actual_jac`` is still needed to compile.

Notice that these modules have initialization routines:

* ``actual_network_init()``

* ``actual_rhs_init()``

These must be called upon initialization. These should be not called
within OpenMP parallel regions, because in general they will modify
global data.

Note, depending on the network, some of these may do nothing, but
these interfaces are all required for maximum flexibility.
