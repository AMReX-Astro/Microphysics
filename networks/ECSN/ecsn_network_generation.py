import pynucastro as pyna
from pynucastro.networks import AmrexAstroCxxNetwork

mylibrary = pyna.rates.ReacLibLibrary()

data_list = mylibrary.get_rates()

all_nuclei = ["p", "he4", "ne20", "o20", "f20",
              "mg24", "al27", "o16", "si28",
              "s32", "p31"]

escn_library = mylibrary.linking_nuclei(all_nuclei, with_reverse=True)

sl = pyna.SuzukiLibrary()
tabular_rates = sl.linking_nuclei(["f20", "o20", "ne20"])

rc = pyna.RateCollection(libraries=[escn_library])

comp = pyna.Composition(rc.get_nuclei())
comp.set_nuc("o16", 0.5)
comp.set_nuc("ne20", 0.3)
comp.set_nuc("mg24", 0.1)
comp.set_nuc("o20", 1.e-5)
comp.set_nuc("f20", 1.e-5)
comp.set_nuc("p", 1.e-5)
comp.set_nuc("he4", 1.e-2)
comp.set_nuc("al27", 1.e-2)
comp.set_nuc("si28", 1.e-2)
comp.set_nuc("s32", 1.e-2)
comp.set_nuc("p31", 1.e-2)
comp.normalize()

# Only select the rates which are >= 1.e-20 at 1.e9 K and 7.e9 g/cm3
# Do not select weak rates from "20180319default2"
new_rate_list = []
rho = 7.e9
T = 1.e9
ydots = rc.evaluate_rates(rho=rho, T=T, composition=comp)
for rate in rc.rates:
    if ydots[rate] >= 1.e-20 and rate.weak == False:
        new_rate_list.append(rate)


wd_net = AmrexAstroCxxNetwork(rates=new_rate_list, libraries=[tabular_rates])
wd_net.write_network()

wd_net.plot(rho, T, comp, outfile="ECSN.png",
            hide_xalpha=True, curved_edges=True,
            node_size=600, node_font_size=11)
