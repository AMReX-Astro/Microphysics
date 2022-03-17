# C-burning with A=23 URCA rate module generator

import pynucastro as pyna
from pynucastro.networks import StarKillerCxxNetwork

all_reactants = ["he4", "c12", "o16"]

def get_3alpha_library(validate=False):

    reaclib_library = pyna.ReacLibLibrary()

    three_alpha_library = reaclib_library.linking_nuclei(all_reactants)

    return three_alpha_library

def doit():

    three_alpha_library = get_3alpha_library()

    net = StarKillerCxxNetwork(libraries=[three_alpha_library], symmetric_screening=True)
    net.write_network()

if __name__ == "__main__":
    doit()
