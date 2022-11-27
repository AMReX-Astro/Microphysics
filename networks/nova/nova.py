import pynucastro as pyna
from pynucastro.networks import AmrexAstroCxxNetwork, RateCollection, Composition

mylibrary = pyna.rates.ReacLibLibrary()

all_nuclei = ["p", "he4", "c12", "c13", "n13", "n14", "n15", "o14", "o15", "o16", "o17", "f17", "f18"]

nova_library = mylibrary.linking_nuclei(all_nuclei, with_reverse=False)

print(nova_library)

net = AmrexAstroCxxNetwork(libraries=[nova_library])
net.write_network()

# make a plot

rc = RateCollection(libraries=[nova_library])

comp = Composition(rc.get_nuclei())
comp.set_solar_like()

rc.plot(outfile="nova.png", rho=1.e4, T=9.e7, comp=comp)

