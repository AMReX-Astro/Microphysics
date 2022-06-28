***************************
Screening of Reaction Rates
***************************

Screening of reaction rates can be computed using several different methods,
controlled by the make parameter ``SCREEN_METHOD``.  For example,

.. prompt:: bash

   make SCREEN_METHOD=screen5

Any of the available screening methods can be used with any reaction network.

The options are:

* ``screen5`` :

  This is the screening routine from the Kepler stellar evolution code
  and is the default used with the distributed versions of the "aprox"
  family of reaction networks.  It uses the screening described in
  :cite:`graboske:1973` for the weak limit and :cite:`alastuey:1978`,
  :cite:`itoh:1979` for the strong limit.

  This is the default screening method.

* ``chugunov2007`` :

  This implements the screening of :cite:`chugunov:2007`.

* ``chugunov2009``

  This implements the screening of :cite:`chugunov:2009`.  This extends the
  work of :cite:`chugunov:2007` to binary mixtures.

