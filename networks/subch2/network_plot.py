# C-burning with A=23 URCA rate module generator

import pynucastro as pyna

library_file = "20180319default2"
mylibrary = pyna.rates.Library(library_file)

all_reactants = ["p",
        "he4", "c12", "o16", "ne20", "mg24", "si28", "s32",
        "ar36", "ca40", "ti44", "cr48", "fe52", "ni56",
        "al27", "p31", "cl35", "k39", "sc43", "v47", "mn51", "co55",
        "c14", "n13", "n14", "o18", "f18", "ne21" ]

subCh = mylibrary.linking_nuclei(all_reactants)

rc = pyna.RateCollection(libraries=[subCh])

comp = pyna.Composition(rc.get_nuclei())
comp.set_all(0.1)
comp.set_nuc("he4", 0.95)
comp.normalize()

rc.plot(outfile="subch2.pdf", rho=1.e6, T=1.e9, comp=comp, hide_xalpha=True,
        size=(1500, 450), node_size=500, node_font_size=11, node_color="#337dff", node_shape="s",
        Z_range=(1,29), rotated=True)
