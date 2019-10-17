# C-burning with A=23 URCA rate module generator

import pynucastro as pyna
from pynucastro.networks import StarKillerNetwork

library_file = "20180319default2"
mylibrary = pyna.rates.Library(library_file)

subCh = pyna.rates.Library()

all_reactants = ["p",
        "he4", "c12", "o16", "ne20", "mg24", "si28", "s32",
        "ar36", "ca40", "ti44", "cr48", "fe52", "ni56",
        "al27", "p31", "cl35", "k39", "sc43", "v47", "mn51", "co55",
        "c14", "n13", "n14", "o18", "f18", "ne21" ]

subCh = mylibrary.linking_nuclei(all_reactants)

net = StarKillerNetwork(libraries=[subCh])
net.write_network()
