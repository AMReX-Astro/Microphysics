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
  :cite:`graboske:1973` for the weak limit and :cite:`jancovici:1977`,
  :cite:`alastuey:1978`, :cite:`itoh:1979` for the strong limit. The
  overall procedure is described in :cite:`Wallace:1982`.

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

* ``chabrier1998``:

  This implements the screening of :cite:`Chabrier_1998` as well as the quantum corrections for strong screening according to screen5. This is suggested in the appendix of :cite:`Calder_2007`. This screening is compatible with NSE calculations unlike ``screen5``, ``chugunov2007``, and ``chugunov2009``. This screening should be valid in the weak screening regime, :math:`\Gamma < 0.1`, and strong screening regime, :math:`1 \lesssim \Gamma \lesssim 160`.

* ``null`` :

  This disables screening by always returning 1 for the screening
  enhancement factor.

Runtime Options
----------------
* ``screening.enable_chabrier1998_quantum_corr = 1`` in the input file enables an additional quantum correction term added to the screening factor when ``SCREEN_METHOD=chabrier1998``. This is disabled by default since ``chabrier1998`` is often used along with ``USE_NSE_NET=TRUE``, and the NSE solver doesn't include quantum corrections.
