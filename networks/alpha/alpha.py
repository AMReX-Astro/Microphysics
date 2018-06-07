# This builds a subch network with all the
# Reaclib rates linking the specified nuclei.

import pynucastro as pyna
from pynucastro.networks import StarKillerNetwork

library_file = "20180319default2"
mylibrary = pyna.rates.Library(library_file)

all_nuclei = ['he4', 'c12', 'o16', 'ne20', 'mg24', 'si28',
              's32', 'ar36', 'ca40', 'ti44', 'cr48', 'fe52', 'ni56', 'zn60']

alpha = mylibrary.linking_nuclei(all_nuclei, with_reverse=True)

net = StarKillerNetwork(libraries=[alpha])
net.write_network()
