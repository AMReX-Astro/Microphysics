.. _neutrino_loss:

***************
Neutrino Losses
***************

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
   void sneut5(const amrex::Real temp, const amrex::Real den,
               const amrex::Real abar, const amrex::Real zbar,
               amrex::Real& snu, amrex::Real& dsnudt, amrex::Real& dsnudd,
               amrex::Real& dsnuda, amrex::Real& dsnudz)

Here, the template parameter, ``do_derivatives``, can be used to disable the code
the computes the derivatives of the neutrino loss, for example, if a numerical Jacobian
is used.  The output is

* ``snu`` : $\epsilon_\nu$, the neutrino loss in erg/g/s

* ``dsnudt`` : $d\epsilon_\nu/dT$

* ``dsnudd`` : $d\epsilon_\nu/d\rho$

* ``dsnuda`` : $d\epsilon_\nu/d\bar{A}$

* ``dsnudz`` : $d\epsilon_\nu/d\bar{Z}$
