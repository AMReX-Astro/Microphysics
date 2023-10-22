# C-burning with A=23 URCA rate module generator

import pynucastro as pyna

rl = pyna.ReacLibLibrary()
rl_lib = rl.linking_nuclei(["n", "p", "he4", "c12", "n13", "o16",
                              "ne20", "na23", "mg23", "mg24"])

# let's remove some unimportant reverse rates (we are not hot
# enough to drive these
rates_to_remove = ["na23(p,n)mg23",
                   "mg24(a,c12)o16",
                   "mg24(g,a)ne20",
                   "mg24(,p)na23",
                   "mg24(,n)mg23",
                   "c12(g,aa)he4",
                   "n13(g,p)c12",
                   "o16(p,a)n13",
                   "o16(g,a)c12",
                   "ne20(g,a)o16",
                   "ne20(a,n)mg23",
                   "na23(p,c12)c12",
                   "ne20(a,c12)c12"]


for r in rates_to_remove:
    print(f"removing {r} : {rl.get_rate_by_.Q}")
    rl_lib.remove_rate(r)

tl = pyna.TabularLibrary()
tl_rates = tl.get_rate_by_name(["na23(,)ne23",
                                "ne23(,)na23",
                                "mg23(,)na23",
                                "n(,)p",
                                "p(,)n"])
print("tl_rates: ")
print(tl_rates)

tl_lib = pyna.Library(rates=tl_rates)

print(tl_lib)

all_lib = rl_lib + tl_lib
print(all_lib)
# we will have duplicate rates -- we want to remove any ReacLib rates
# that we have tabular rates for

dupes = all_lib.find_duplicate_links()

rates_to_remove = []
for d in dupes:
    for r in d:
        if isinstance(r, pyna.rates.ReacLibRate):
            rates_to_remove.append(r)

for r in rates_to_remove:
    all_lib.remove_rate(r)


urca_net = pyna.AmrexAstroCxxNetwork(libraries=[all_lib])

comp = pyna.Composition(urca_net.get_nuclei())
comp.set_all(0.1)
comp.set_nuc("c12", 0.5)
comp.set_nuc("o16", 0.5)
comp.normalize()

urca_net.plot(outfile="urca_medium.png", rho=1.e9, T=6.e8, comp=comp,
              rotated=True, hide_xalpha=True, curved_edges=True,
              size=(1500, 450),
              node_size=500, node_font_size=11, node_color="#337dff", node_shape="s")

urca_net.write_network()
