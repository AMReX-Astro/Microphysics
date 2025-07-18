# `he-burn-19a`

This is based on the original `subch_simple` network.  It tries to
create an aprox13-like alpha-chain + include a bypass rate for
C12(a, g)O16 discussed in Shen & Bildsten (2008).

Shen & Bildsten discuss the sequences:

* C14(a, g)O18(a, g)Ne22 at high temperatures (T > 1 GK).  We don't
  consider this.

* N14(a, g)F18(a, p)Ne21 is the one they consider important, since it
  produces protons that are then available for C12(p, g)N13(a, p)O16.
  This sequence would leave Ne21 as an endpoint.

  We don't include this sequence, but instead allow for n14 to be
  present as an end-state of H burning, and connect it directly to
  Ne20 via N14(1.5a,g)Ne20, which is an approximation used in
  `aprox19`.  We also include a reverse sequence from there:
  O16(pp,a)N14.

For the C12+C12, C12+O16, and O16+O16 rates, we also need to include
C12(C12,n)Mg23(n, g)Mg24, O16(O16, n)S31(n, g)S32, O16(C12, n)Si27(n,
g)Si28.  Since the neutron captures on those intermediate nuclei are
so fast, we leave those out and take the forward rate to just be the
first rate.  We do not include reverse rates for these processes.

The `he-burn-19a` network simplifies subch_full by including the
following approximations:

* Approximate rates by creating an effective rate for
  (a, g) and (a, p)(p, g) assuming proton equilibrium. This allows
  us to remove 6 intermediate nuclei. This approximation is effective
  when temperatures are below 2.5 GK.

* Omitting the reverse of C12+C12, C12+O16, and O16+O16 -- these
  reverse rates are not in `aprox13`

* Omitting the C12+Ne20 rates and their reverses

* (a,g) links between Na23 and Al27 and between Al27 and P31 -- these
  are not in `aprox13`
