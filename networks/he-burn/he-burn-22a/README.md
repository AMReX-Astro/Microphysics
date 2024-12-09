# subch_simple

The original subch_full network tries to recreate an aprox13 alpha-chain
+ including a bypass rate for c12(a, g)o16 discussed in Shen & Bildsten (2008).

Shen & Bildsten discuss the sequences:

* c14(a, g)o18(a, g)ne22 at high temperatures (T > 1 GK).  We don't consider
  this.

* n14(a, g)f18(a, p)ne21 is the one they consider important, since it
  produces protons that are then available for c12(p, g)n13(a, p)o16.

  Note that, this leaves ne21 as an endpoint, which we need to connect by
  including na22. However, adding na22 only makes it as another endpoint,
  failing to resolve the issue, resulting in a network with a second
  endpoint other than ni56.

For the c12+c12, c12+o16, and o16+o16 rates, we also need to include
c12(c12,n)mg23(n, g)mg24, o16(o16, n)s31(n, g)s32,
o16(c12, n)si27(n, g)si28.  Since the neutron captures on those
intermediate nuclei are so fast, we leave those out and take the
forward rate to just be the first rate.  We do not include reverse
rates for these processes.

The subch_simple network simplifies subch_full by including the
following approximations:

* Approximate rates by creating an effective rate for
  (a, g) and (a, p)(p, g) assuming proton equilibrium. This allows
  us to remove 6 intermediate nuclei. This approximation is effective
  when temperatures are below 2.5 GK.

* The reverse of C12+C12, C12+O16, and O16+O16 -- these reverse rates
  are not in aprox13

* The C12+Ne20 rates and their reverses

* (a,g) links between Na23 and Al27 and between Al27 and P31 -- these
  are not in aprox13
