import pynucastro as pyna
from pynucastro.networks import AmrexAstroCxxNetwork

rl = pyna.ReacLibLibrary()

# note: adding Li6 makes this network stiff

all_nuclei = ["p", "h2", "he3", "he4",
              "li7", "be7", "be8", "b8", "b11",
              "c12", "c13", "n13", "n14", "n15",
              "o14", "o15", "o16", "o17", "o18",
              "f17", "f18", "f19",
              "ne18", "ne19", "ne20", "ne21"]

nova_library = rl.linking_nuclei(all_nuclei, with_reverse=True)

tl = pyna.TabularLibrary(ordering=["ffn", "langanke"])
#weak_library = tl.linking_nuclei(all_nuclei)


net = AmrexAstroCxxNetwork(libraries=[nova_library])  # + weak_library])

# start of burning state
rho = 1.7e3
T = 7.e7
comp = pyna.Composition(net.unique_nuclei)
comp.set_solar_like()
state1 = (rho, T, comp)

# hotter state
rho = 1.e3
T = 1.e8
comp = pyna.Composition(net.unique_nuclei)
comp.set_equal()
state2 = (rho, T, comp)

cutoff_ratio = 1.e-50

unimportant_rates = net.find_unimportant_rates([state1, state2], cutoff_ratio)
print("unimportant rates: ")
for k, v in unimportant_rates.items():
    print(k, v)

net.remove_rates(unimportant_rates)

net.write_network()

net.plot(rotated=False, outfile="nova.png",
         hide_xalpha=True, hide_xp=True,
         curved_edges=False,
         node_size=400, node_font_size=10)

net.summary()
