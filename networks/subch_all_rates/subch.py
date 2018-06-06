# This builds a subch network with all the
# Reaclib rates linking the specified nuclei.

import pynucastro as pyna
from pynucastro.networks import StarKillerNetwork

library_file = "20180319default2"
mylibrary = pyna.rates.Library(library_file)

all_nuclei = ['p', 'he4', 'c12', 'c14', 'n13', 'n14', 'o16',
              'o18', 'f18', 'ne20', 'ne21']

subCh = mylibrary.linking_nuclei(all_nuclei, with_reverse=True)

net = StarKillerNetwork(libraries=[subCh])
net.write_network()
