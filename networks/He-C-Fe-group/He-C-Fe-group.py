import pynucastro as pyna
from pynucastro.rates import ReacLibRate, TabularRate

DO_DERIVED_RATES = True

reaclib_lib = pyna.ReacLibLibrary()
weak_lib = pyna.TabularLibrary()

# these are the nuclei we have in subch_simple

all_reactants = ["p",
                 "he4", "c12", "o16", "ne20", "mg24", "si28", "s32",
                 "ar36", "ca40", "ti44", "cr48", "fe52", "ni56",
                 "cu59", "zn60",
                 "al27", "p31", "cl35", "k39", "sc43", "v47", "mn51", "co55",
                 "n13", "n14", "f18", "ne21", "na22", "na23"]

# create a library of ReacLib rates

core_lib = reaclib_lib.linking_nuclei(all_reactants)

# in this list, we have the reactants, the actual reactants,
# and modified products that we will use instead

other_rates = [("c12(c12,n)mg23", "mg24"),
               ("o16(o16,n)s31", "s32"),
               ("o16(c12,n)si27", "si28")]

for r, mp in other_rates:
    _r = reaclib_lib.get_rate_by_name(r)
    _r.modify_products(mp)
    core_lib.add_rate(_r)

# finally, the aprox nets don't include the reverse rates for
# C12+C12, C12+O16, and O16+O16, so remove those

for r in core_lib.get_rates():
    if sorted(r.products) in [[pyna.Nucleus("c12"), pyna.Nucleus("c12")],
                              [pyna.Nucleus("c12"), pyna.Nucleus("o16")],
                              [pyna.Nucleus("o16"), pyna.Nucleus("o16")]]:
        core_lib.remove_rate(r)

# C12+Ne20 and reverse
# (a,g) links between Na23 and Al27
# (a,g) links between Al27 and P31

rates_to_remove = ["p31(p,c12)ne20",
                   "si28(a,c12)ne20",
                   "ne20(c12,p)p31",
                   "ne20(c12,a)si28",
                   "na23(a,g)al27",
                   "al27(g,a)na23",
                   "al27(a,g)p31",
                   "p31(g,a)al27"]

for r in rates_to_remove:
    print("removing: ", r)
    _r = core_lib.get_rate_by_name(r)
    core_lib.remove_rate(_r)

# now create a list of iron group nuclei and find both the
# ReacLib and weak / tabular rates linking these.

iron_peak = ["n", "p", "he4",
             "mn51", "mn55",
             "fe52", "fe53", "fe54", "fe55", "fe56",
             "co55", "co56", "co57",
             "ni56", "ni57", "ni58",
             "cu59", "zn60"]

iron_reaclib = reaclib_lib.linking_nuclei(iron_peak)

iron_weak_lib = weak_lib.linking_nuclei(iron_peak)

# add the libraries

all_lib = core_lib + iron_reaclib + iron_weak_lib

if DO_DERIVED_RATES:
    rates_to_derive = []
    for r in all_lib.get_rates():
        if r.reverse:
            # this rate was computed using detailed balance, regardless
            # of whether Q < 0 or not.  We want to remove it and then
            # recompute it
            rates_to_derive.append(r)

    # now for each of those derived rates, look to see if the pair exists

    for r in rates_to_derive:
        fr = all_lib.get_rate_by_nuclei(r.products, r.reactants)
        if fr:
            print(f"modifying {r} from {fr}")
            all_lib.remove_rate(r)
            d = pyna.DerivedRate(rate=fr, compute_Q=False, use_pf=True)
            all_lib.add_rate(d)

# we will have duplicate rates -- we want to remove any ReacLib rates
# that we have tabular rates for

dupes = all_lib.find_duplicate_links()

rates_to_remove = []
for d in dupes:
    for r in d:
        if isinstance(r, ReacLibRate):
            rates_to_remove.append(r)

for r in rates_to_remove:
    all_lib.remove_rate(r)

# combine all three libraries into a single network

net = pyna.AmrexAstroCxxNetwork(libraries=[all_lib],
                                symmetric_screening=False)


# now we approximate some (alpha, p)(p, gamma) links

net.make_ap_pg_approx(intermediate_nuclei=["cl35", "k39", "sc43", "v47"])
net.remove_nuclei(["cl35", "k39", "sc43", "v47"])

# let's make a figure

T = 6.e9
rho = 9.e7
comp = pyna.Composition(net.unique_nuclei)
comp.set_all(1.0)
comp.normalize()

fig = net.plot(rho=rho, T=T, comp=comp,
               rotated=True, curved_edges=True, hide_xalpha=True,
               size=(1800, 900),
               node_size=500, node_shape="s", node_font_size=10)

fig.savefig("He-C-Fe-group.png")

print(f"number of nuclei = {len(net.unique_nuclei)}")
print(f"number of ReacLib rates = {len(net.reaclib_rates)}")
print(f"number of tabular rates = {len(net.tabular_rates)}")

net.write_network()
