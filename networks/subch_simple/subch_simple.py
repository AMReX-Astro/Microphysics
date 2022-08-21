# an approximate network for He/C burning with key rates
# to bypass the C12(a,g)O16 rate.  This version uses some
# (a,p)(p,g) approximations.

import pynucastro as pyna
from pynucastro.networks import StarKillerCxxNetwork

def get_library():

    reaclib_lib = pyna.ReacLibLibrary()

    all_reactants = ["p",
                     "he4", "c12", "o16", "ne20", "mg24", "si28", "s32",
                     "ar36", "ca40", "ti44", "cr48", "fe52", "ni56",
                     "al27", "p31", "cl35", "k39", "sc43", "v47", "mn51", "co55",
                     "n13", "n14", "f18", "ne21", "na22", "na23"]

    subch = reaclib_lib.linking_nuclei(all_reactants)

    # in this list, we have the reactants, the actual reactants,
    # and modified products that we will use instead
    other_rates = [(("c12", "c12"), ("mg23", "n"), ("mg24")),
                   (("o16", "o16"), ("s31", "n"), ("s32")),
                   (("c12", "o16"), ("si27", "n"), ("si28"))]

    for r, p, mp in other_rates:
        rfilter = pyna.rates.RateFilter(reactants=r, products=p)
        _library = reaclib_lib.filter(rfilter)
        r = _library.get_rates()[0]
        r.modify_products(mp)
        subch += _library

    return subch

def doit():

    subch = get_library()

    # these are the rates that we are going to allow to be optionally
    # zeroed
    r1 = subch.get_rate("p_c12__n13")
    r2 = subch.get_rate("he4_n13__p_o16")

    net = StarKillerCxxNetwork(libraries=[subch], symmetric_screening=True, disable_rate_params=[r1, r2])
    net.make_ap_pg_approx(intermediate_nuclei=["cl35", "k39", "sc43", "v47", "mn51", "co55"])
    net.remove_nuclei(["cl35", "k39", "sc43", "v47", "mn51", "co55"])

    # finally, the aprox nets don't include the reverse rates for
    # C12+C12, C12+O16, and O16+O16, so remove those

    rates_to_remove = []
    for r in net.rates:
        if sorted(r.products) in [[pyna.Nucleus("c12"), pyna.Nucleus("c12")],
                                  [pyna.Nucleus("c12"), pyna.Nucleus("o16")],
                                  [pyna.Nucleus("o16"), pyna.Nucleus("o16")]]:
            rates_to_remove.append(r)

        # Q = 1.9
        # if (sorted(r.reactants) == sorted([pyna.Nucleus("p"), pyna.Nucleus("p31")]) and
        #     sorted(r.products) == sorted([pyna.Nucleus("he4"), pyna.Nucleus("si28")])):
        #     rates_to_remove.append(r)

        # Q = -10.1
        if (sorted(r.reactants) == sorted([pyna.Nucleus("p"), pyna.Nucleus("p31")]) and
            sorted(r.products) == sorted([pyna.Nucleus("c12"), pyna.Nucleus("ne20")])):
            rates_to_remove.append(r)

        # Q = -12.0
        if (sorted(r.reactants) == sorted([pyna.Nucleus("he4"), pyna.Nucleus("si28")]) and
            sorted(r.products) == sorted([pyna.Nucleus("c12"), pyna.Nucleus("ne20")])):
            rates_to_remove.append(r)

        # Q = -10.1
        if (sorted(r.products) == sorted([pyna.Nucleus("p"), pyna.Nucleus("p31")]) and
            sorted(r.reactants) == sorted([pyna.Nucleus("c12"), pyna.Nucleus("ne20")])):
            rates_to_remove.append(r)

        # Q = -12.0
        if (sorted(r.products) == sorted([pyna.Nucleus("he4"), pyna.Nucleus("si28")]) and
            sorted(r.reactants) == sorted([pyna.Nucleus("c12"), pyna.Nucleus("ne20")])):
            rates_to_remove.append(r)

        # if (sorted(r.reactants) == sorted([pyna.Nucleus("p"), pyna.Nucleus("al27")]) and
        #     sorted(r.products) == sorted([pyna.Nucleus("he4"), pyna.Nucleus("mg24")])):
        #     rates_to_remove.append(r)

        # If (sorted(r.reactants) == sorted([pyna.Nucleus("p"), pyna.Nucleus("ne21")]) and
        #     sorted(r.products) == sorted([pyna.Nucleus("he4"), pyna.Nucleus("f18")])):
        #     rates_to_remove.append(r)

        # if (sorted(r.reactants) == sorted([pyna.Nucleus("p"), pyna.Nucleus("o16")]) and
        #     sorted(r.products) == sorted([pyna.Nucleus("he4"), pyna.Nucleus("n13")])):
        #     rates_to_remove.append(r)

    for r in rates_to_remove:
        print("removing: ", r)

    net.remove_rates(rates_to_remove)

    print(f"number of nuclei: {len(net.unique_nuclei)}")
    print(f"number of rates: {len(net.rates)}")

    comp = pyna.Composition(net.get_nuclei())
    comp.set_all(0.1)
    comp.set_nuc("he4", 0.95)
    comp.normalize()

    net.plot(outfile="subch_simple.png", rho=1.e6, T=1.e9, comp=comp,
             rotated=True, hide_xalpha=True, curved_edges=True,
             size=(1500, 450),
             node_size=500, node_font_size=11, node_color="#337dff", node_shape="s",
             Z_range=(1,29))

    net.write_network()


if __name__ == "__main__":
    doit()
