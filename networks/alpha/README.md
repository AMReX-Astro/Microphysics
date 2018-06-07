# Sub-Ch reaction network

This is a He burning network with the links from Shen & Bildsten (2009) to
bypass the 12C(a,g)16O rate in favor of 12C(p,g)13N(a,p)16O (where the
p comes from 14N(a,g)18F(a,p)21Ne.  So to get this rate going, we need to
have some N14 in the envelope.

This network was generated using pynucastro and the ReacLib rate library and
can be reproduced via:

```
$ python subch.py
```

# Updates

The subch_all_rates network incorporates all Reaclib rates linking
the nuclei in the subch network for comparison.
