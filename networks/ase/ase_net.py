# Test a fast reaction cycle for dynamic equilibrium
# S32(g,p)P31(g,p)Si30(g,n)Si29(g,n)Si28(a,g)32S

import pynucastro as pyna
from pynucastro.networks import StarKillerCxxNetwork


reaclib_lib = pyna.ReacLibLibrary()

nucs = ["p", "n", "he4",
        # "c12", "n13", "n14",
        # "o16", "f18", "ne20",
        # "ne21", "na22", "na23",
        # "mg24", "al27" , "si28",
        # "p31", "s32", "ar36",
        # "ca40", "ti44", "cr48",
        # "fe52", "ni56"]
        
        "si28", "si29",
        "si30", "p31", "s32"]

ase = reaclib_lib.linking_nuclei(nucs)

net = StarKillerCxxNetwork(libraries=[ase])
net.write_network()
