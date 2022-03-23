# C-burning with A=23 URCA rate module generator

import pynucastro as pyna
from pynucastro.networks import StarKillerCxxNetwork

all_reactants = ["n", "p",
                 "he4", "c12", "o16", "ne20", "mg24", "si28", "s32",
                 "ar36", "ca40", "ti44", "cr48", "fe52", "ni56",
                 "al27", "p31", "cl35", "k39", "sc43", "v47", "mn51", "co55",
                 "c14", "n13", "n14", "o18", "f18", "ne21",
                 "mg23", "na23", "si27", "s31"]

def get_subch2_library(validate=False):

    reaclib_library = pyna.ReacLibLibrary()

    subch_library = reaclib_library.linking_nuclei(all_reactants)

    if validate:
        subch_library.validate(reaclib_library)

    all_reactants.remove("n")

    no_neutron_library = reaclib_library.linking_nuclei(all_reactants)

    return subch_library, no_neutron_library


def doit():

    subch_library, no_neutron_library = get_subch2_library()

    n_rates = [r for _, r in (subch_library - no_neutron_library)._rates.items()]

    net = StarKillerCxxNetwork(libraries=[subch_library], symmetric_screening=True)
    net.write_network()

if __name__ == "__main__":
    doit()
