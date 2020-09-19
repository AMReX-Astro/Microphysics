# aprox19

This is a rewrite of Frank Timmes' aprox19 network, based on aprox19.tbz
as downloaded from his website on 2015-12-31. It follows the same reasoning
as the aprox13 network in this repository, so see the README in that directory
for information on our integration strategy. This directory is a copy of the
files in the aprox13/ directory, with appropriate changes for the number of
species and reaction rates.

We thank Frank for allowing us to redistribute these routines.

## NSE option

``aprox19`` can model nuclear statistical equilibrium at high
densities.  This is accomplished by switching from integrating the ODE
system of 19 nuclei to doing a table look-up on the results from a
network with 125 nuclei.

The NSE table was provided by Stan Woosley:

from Stan:

> Probably all the code needs for hydrodynamics and energy generation
> is abar and Ye (given rho and T) so to follow explosions you can use
> the smaller table. It also gives crude measures of nucleosyntheis in
> n+p+helium, si-ca, fe group. Probaly good enough for NSE
>
> If you are following bufk nucleosynthesis you can use the larger table
> which gives 125 isotopes packaged into 19 (CASTRO used to use the 19
> iso network of KEPLER, maybe it no longer does?). These abundances can
> be interpolated by a simple extension of the usetable.f routine to
> more variables

The mapping of the nuclei from the large network to the 19 we carry is:

* 1: is H1  is 0
* 2: is mostly 3He
* 3: 4He
* 4: 12C
* 5: 14N
* 6: 16O
* 7: 20Ne
* 8: 24Mg
* 9: 28Si
* 10: 32S
* 11: 36Ar
* 12: 40Ca
* 13: 44Ti
* 14: 48Cr
* 15: 52Fe
* 16: all oher iron group except 56Ni
* 17: 56Ni
* 18: neutrons
* 19: protons

More from Stan:

> The Ye of those 19 isotopes may not agree exactly with the Ye of the
> table and the latter should be used in the hydro and EOS. In fact
> since the A of abund 16 is not given, Ye is indeterminate.
>
> Carrying all these isotopes merely adds to the memory so long as you
> are in nse, but if you are following a continuum of burning from He
> (or H) to NSE you may need them.  In fact, I think all burning
> stages can be handled with tables e.g. for carbon burning makes
> tables of compsition for carbon burned at constant T and rho and
> interpolate in X12, T and rho, though some thought must e given to
> the initial carbon abundance.
>
> Anyway the above packing is for a calculation with 125 isotopes. I can
> send you a larger table with all 125 species and you can pack it as
> you wish. The file is 158 MB Actually I just sent it by dropbox


Here is a sample dYe/dt is calculated using Fuller rates for this
composition NSE for T9, rho, Ye = 3.16, 1.e8, .5 Its mostly 56Ni 28.0
56.0 9.920E-01

```
  3.16E+00  1.00E+08  5.00E-01
   0.0   1.0   6.149E-18   1.0   1.0   3.141E-05
   1.0   2.0   5.188E-21   1.0   3.0   4.965E-31
   2.0   3.0   1.538E-19   2.0   4.0   5.488E-07
   2.0   5.0   8.226E-28   3.0   5.0   7.983E-17
   6.0  12.0   3.216E-15   8.0  16.0   5.311E-14
  10.0  20.0   1.008E-16  11.0  23.0   1.294E-22
  12.0  24.0   3.484E-12  12.0  25.0   6.399E-20
  12.0  26.0   2.960E-23  13.0  26.0   4.577E-17
  13.0  27.0   8.624E-17  13.0  28.0   2.992E-24
  14.0  28.0   1.272E-06  14.0  29.0   5.002E-13
  14.0  30.0   1.200E-16  15.0  30.0   2.119E-11
  15.0  31.0   3.205E-12  16.0  31.0   2.333E-09
  15.0  32.0   5.157E-19  16.0  32.0   6.497E-06
  15.0  33.0   6.483E-24  16.0  33.0   9.448E-12
  16.0  34.0   2.171E-14  17.0  35.0   3.794E-11
  17.0  36.0   1.442E-17  18.0  36.0   1.039E-05
  17.0  37.0   1.709E-21  18.0  37.0   2.516E-11
  18.0  38.0   2.728E-13  18.0  39.0   4.245E-22
  19.0  39.0   4.852E-10  18.0  40.0   1.721E-27
  19.0  40.0   4.656E-17  20.0  40.0   6.916E-05
  19.0  41.0   5.917E-22  20.0  41.0   6.967E-11
  19.0  42.0   1.339E-29  20.0  42.0   1.026E-13
  19.0  43.0   6.402E-35  20.0  43.0   2.629E-20
  21.0  43.0   2.583E-12  20.0  44.0   9.057E-24
  21.0  44.0   7.455E-17  22.0  44.0   4.369E-07
  21.0  45.0   4.493E-19  22.0  45.0   1.028E-10
  20.0  46.0   7.435E-36  21.0  46.0   5.855E-25
  22.0  46.0   2.656E-11  21.0  47.0   6.112E-29
  22.0  47.0   2.003E-16  23.0  47.0   1.693E-09
  20.0  48.0   6.858E-49  21.0  48.0   1.448E-35
  22.0  48.0   4.188E-19  23.0  48.0   5.716E-13
  24.0  48.0   3.612E-05  21.0  49.0   3.609E-40
  22.0  49.0   1.653E-25  23.0  49.0   9.805E-15
  24.0  49.0   1.105E-07  22.0  50.0   3.375E-29
  23.0  50.0   7.061E-20  24.0  50.0   4.772E-08
  22.0  51.0   1.131E-38  23.0  51.0   6.795E-23
  24.0  51.0   1.069E-12  25.0  51.0   2.601E-06
  22.0  52.0   5.116E-47  23.0  52.0   7.728E-31
  24.0  52.0   1.153E-14  25.0  52.0   1.621E-09
  26.0  52.0   6.167E-03  23.0  53.0   5.212E-38
  24.0  53.0   1.332E-21  25.0  53.0   9.123E-11
  26.0  53.0   2.544E-05  23.0  54.0   1.793E-48
  24.0  54.0   6.539E-27  25.0  54.0   1.499E-16
  26.0  54.0   3.758E-05  27.0  54.0   6.603E-06
  24.0  55.0   1.819E-36  25.0  55.0   1.318E-20
  26.0  55.0   6.657E-10  27.0  55.0   1.012E-03
  24.0  56.0   2.376E-44  25.0  56.0   6.359E-29
  26.0  56.0   6.959E-13  27.0  56.0   1.230E-07
  28.0  56.0   9.920E-01  25.0  57.0   6.007E-36
  26.0  57.0   5.610E-20  27.0  57.0   5.570E-10
  28.0  57.0   5.882E-04  25.0  58.0   3.779E-45
  26.0  58.0   3.792E-25  27.0  58.0   8.540E-16
  28.0  58.0   2.178E-05  26.0  59.0   3.766E-34
  27.0  59.0   3.999E-20  28.0  59.0   1.732E-10
  29.0  59.0   8.068E-07  26.0  60.0   5.413E-41
  27.0  60.0   7.227E-28  28.0  60.0   2.320E-13
  30.0  60.0   9.784E-07  26.0  61.0   1.015E-51
  27.0  61.0   8.329E-34  28.0  61.0   6.271E-20
  26.0  62.0   1.008E-59  27.0  62.0   4.070E-43
  28.0  62.0   1.802E-24  27.0  63.0   3.144E-50
  28.0  63.0   8.278E-33  27.0  64.0   1.316E-60
  28.0  64.0   1.132E-38  28.0  65.0   3.709E-48
  28.0  66.0   4.697E-55
``
which is summarized in the table as

``
     9.50000     8.00000 5.00000E-01 3.19627E-05 8.77542E-05 9.99880E-01 5.58734E+01 8.64229E+00 4.96722E-04 0.00000E+00 1.58959E-19 5.48846E-07 3.21588E-15 0.00000E+00 5.31077E-14 1.00820E-16 3.48367E-12 1.27462E-06 6.49684E-06 1.03892E-05 6.91567E-05 4.38678E-07 3.88804E-05 6.19209E-03 1.66798E-03 9.91981E-01 6.14888E-18 3.14138E-05
``

regarding the term dYedt:

What is tabulated is wrate which is the sum of all electron capture
and positron decay rates times the appropriate abundances minus a
similar rate for beta decay and positron capture

 wrate=rectot+rpdtot-redtot-rpctot

So if electron capture dominates wrate is positive

The inference is that it is a positive term that is subtracted from Ye