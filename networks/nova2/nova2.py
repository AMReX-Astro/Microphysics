import pynucastro as pyna
from pynucastro.networks import AmrexAstroCxxNetwork

mylibrary = pyna.rates.ReacLibLibrary()

all_nuclei = ["p", "h2", "he3", "he4", "be7", "b8",
              "c12", "c13", "n13", "n14", "n15",
              "o14", "o15", "o16", "o17", "f17", "f18"]

nova_library = mylibrary.linking_nuclei(all_nuclei, with_reverse=False)

print(nova_library)

net = AmrexAstroCxxNetwork(libraries=[nova_library])
net.write_network()

rc = pyna.RateCollection(libraries=[nova_library])

comp = pyna.Composition(rc.get_nuclei())
comp.set_solar_like()

rc.plot(outfile="nova.png", rho=1.e4, T=9.e7, comp=comp)

