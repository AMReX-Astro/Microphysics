.. _neutrino_loss:

***************
Neutrino Losses
***************

.. index:: neutrino_cooling.include_recomb

We model thermal neutrino losses (plasma, photo-, pair-,
recombination, and Bremsstrahlung) using the method described in
:cite:`itoh:1996`.  This neutrino loss term is included in all of the
reaction networks by default, and modifies the energy equation to have
the form (for Strang splitting):

.. math::

   \frac{de}{dt} = \epsilon - \epsilon_\nu

where $\epsilon_\nu$ are the thermal neutrino losses.

.. note::

   Thermal neutrino losses can be disabled at compile time by setting the make
   variable ``USE_NEUTRINOS = FALSE``.

The interface for the neutrino loss function is:

.. code:: c++

   template <int do_derivatives>
   AMREX_GPU_HOST_DEVICE AMREX_INLINE
   void sneut5(const number_t& temp, const amrex::Real den,
                   const number_t& abar, const number_t& zbar,
                   amrex::Real& pair, amrex::Real& phot,
                   amrex::Real& plas, amrex::Real& brem)

Here, the template parameter, ``do_derivatives``, can be used to disable the code
that computes the derivatives of the neutrino loss, for example, if a numerical Jacobian
is used.  The output is

* ``snu`` : $\epsilon_\nu$, the neutrino loss in erg/g/s

* ``dsnudt`` : $d\epsilon_\nu/dT$

* ``dsnudd`` : $d\epsilon_\nu/d\rho$

* ``dsnuda`` : $d\epsilon_\nu/d\bar{A}$

* ``dsnudz`` : $d\epsilon_\nu/d\bar{Z}$

* ``pair`` : contribution from pair neutrino loss in erg/g/s

* ``phot`` : contribution from photo-neutrino loss in erg/g/s

* ``plas`` : contribution from plasma neutrino loss in erg/g/s

* ``brem`` : contribution from Bremsstrahlung neutrino loss in erg/g/s

Additionally, a second method is provided to compute the neutrino losses based on the simplified approximations in :cite:`kippenhahn:1990`.

The interface for the Kippenhahn neutrino loss function is:

.. code:: c++

   template <int do_derivatives>
   AMREX_GPU_HOST_DEVICE AMREX_INLINE
   void kipp(const amrex::Real& temp, const amrex::Real& rho, const amrex::Real& abar,
             const amrex::Real& zbar, amrex::Real& snu, amrex::Real& dsnudt,
             amrex::Real& dsnudrho, amrex::Real& dsnudz, amrex::Real& dsnuda,
             amrex::Real& pair, amrex::Real& phot, amrex::Real& plas, amrex::Real& brem)

To use this method, set the make variable ``NEUTRINO_METHOD = kipp`` during compilation.
Output is the same as for the previous method.

.. note::

   By default, we do not include the recombination terms when calculating the total losses since its effect is negligible.
   This is controlled by the ``include_recomb`` parameter defined in ``Microphysics/neutrinos/_parameters``.
   To include the recombination terms, set ``neutrino_cooling.include_recomb = 1`` in the inputs file.
