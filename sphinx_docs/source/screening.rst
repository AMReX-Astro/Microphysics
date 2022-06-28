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

  This implements the screening of :cite:`chugunov:2007`, following
  :cite:`yakovlev:2006` to extend to binary mixtures.

* ``chugunov2009``

  This implements the screening of :cite:`chugunov:2009`.  The main
  difference is that the 2007 one calculates an effective coupling
  parameter based on the two ions and then treats it as a
  one-component plasma, while this version (2009) treats it fully as a
  multi-component plasma (which is significantly more expensive)

  This includes the portion in the appendix that blends in the weak
  screening limit.
