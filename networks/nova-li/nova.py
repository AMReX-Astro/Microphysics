import pynucastro as pyna
from pynucastro.networks import AmrexAstroCxxNetwork

rl = pyna.ReacLibLibrary()

all_nuclei = ["p", "h2", "he3", "he4",
              "li6", "li7", "be7", "be8", "b8",
              "c12", "c13", "n13", "n14", "n15",
              "o14", "o15", "o16", "o17", "o18",
              "f17", "f18", "f19",
              "ne18", "ne19", "ne20", "ne21"]

nova_library = rl.linking_nuclei(all_nuclei, with_reverse=False)

tl = pyna.TabularLibrary(ordering=["ffn", "langanke"])
#weak_library = tl.linking_nuclei(all_nuclei)


net = AmrexAstroCxxNetwork(libraries=[nova_library])  # + weak_library])
net.write_network()

rc = pyna.RateCollection(libraries=[nova_library])

rc.plot(rotated=False, outfile="nova.png",
        hide_xalpha=True, hide_xp=True,
        curved_edges=False,
        node_size=400, node_font_size=10)

net.summary()
