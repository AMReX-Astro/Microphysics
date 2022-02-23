# C-burning with A=23 URCA rate module generator

import pynucastro as pyna
from pynucastro.networks import StarKillerCxxNetwork

library_file = "20180319default2"

all_reactants = ["n", "p",
                 "he4", "c12", "o16", "ne20", "mg24", "si28", "s32",
                 "ar36", "ca40", "ti44", "cr48", "fe52", "ni56",
                 "al27", "p31", "cl35", "k39", "sc43", "v47", "mn51", "co55",
                 "c14", "n13", "n14", "o18", "f18", "ne21",
                 "mg23", "na23", "si27", "s31"]

def doit():
    reaclib_library = pyna.rates.Library(library_file)


    subch_library = reaclib_library.linking_nuclei(all_reactants)

    # generate a report about any missing rates that we might want to include

    subch_library.validate(reaclib_library)

    net = StarKillerCxxNetwork(libraries=[subch_library])
    net.write_network()

if __name__ == "__main__":
    doit()
