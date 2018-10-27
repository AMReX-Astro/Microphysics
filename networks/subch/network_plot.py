# C-burning with A=23 URCA rate module generator

import pynucastro as pyna
from pynucastro.networks import StarKillerNetwork

library_file = "20180319default2"
mylibrary = pyna.rates.Library(library_file)

subCh = pyna.rates.Library()

all_reactants = [(("he4", "he4", "he4"), ("c12")),
        (("c12", "he4"), ("o16")),
        (("n14", "he4"), ("f18")),
        (("f18", "he4"), ("ne21", "p")),
        (("c12", "p"), ("n13")),
        (("n13", "he4"), ("o16", "p")),
        (("o16", "he4"), ("ne20")),
        (("c14", "he4"), ("o18"))]

for r, p in all_reactants:
    if not isinstance(p, tuple):
        p = p,
    rfilter = pyna.rates.RateFilter(reactants=r, products=p)
    _library = mylibrary.filter(rfilter)
    subCh += _library

rc = pyna.RateCollection(libraries=[subCh])

comp = pyna.Composition(rc.get_nuclei())
comp.set_all(0.1)
comp.set_nuc("he4", 0.95)
comp.normalize()

rc.plot(outfile="subch.png", rho=1.e6, T=3.e8, comp=comp)
