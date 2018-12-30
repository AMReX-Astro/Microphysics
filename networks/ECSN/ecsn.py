# C-burning with A=23 URCA rate module generator

import pynucastro as pyna
from pynucastro.networks import StarKillerNetwork

library_file = "20180319default2"
mylibrary = pyna.rates.Library(library_file)

all_nuclei = ["p", "he4",
              "ne20", "o20", "f20",
              "mg24", "al27", "o16",
              "si28", "s32", "p31"]

ecsn_library = mylibrary.linking_nuclei(all_nuclei, with_reverse=True)

print(ecsn_library)

net = StarKillerNetwork(libraries=[ecsn_library])
net.write_network()
