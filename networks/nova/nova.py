import pynucastro as pyna
from pynucastro.networks import AmrexAstroCxxNetwork

mylibrary = pyna.rates.ReacLibLibrary()

all_nuclei = ["p", "h2", "he3", "he4", "be7", "b8",
              "c12", "c13", "n13", "n14", "n15",
              "o14", "o15", "o16", "o17", "f17", "f18"]

nova_library = mylibrary.linking_nuclei(all_nuclei, with_reverse=False)

# pp_chain = ["p", "h2", "he3", "he4", "be7", "b8"]
# pp_library = mylibrary.linking_nuclei(pp_chain, with_reverse=False)

net = AmrexAstroCxxNetwork(libraries=[nova_library])
net.write_network()

rc = pyna.RateCollection(libraries=[nova_library])
# filter = pyna.RateCollection(libraries=[pp_library])

comp = pyna.Composition(rc.get_nuclei())
comp.set_solar_like()

rho = 1.e3
T = 1.e8

edge_labels = {(pyna.Nucleus("he4"), pyna.Nucleus("c12")):
               r"$\alpha(\alpha\alpha,\gamma){}^{12}\mathrm{C}$"}

rc.plot(rho, T, comp, rotated=False, outfile="nova.png", N_range=(-1,10) , Z_range=(0,10),  hide_xalpha=True, hide_xp=True,
        curved_edges=True, edge_labels=edge_labels, node_size=300, node_font_size=10)
