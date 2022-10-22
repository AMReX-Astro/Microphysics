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
    other_rates = [("c12(c12,n)mg23", "mg24"),
                   ("o16(o16,n)s31", "s32"),
                   ("o16(c12,n)si27", "si28")]

    for r, mp in other_rates:
        _r = reaclib_lib.get_rate_by_name(r)
        _r.modify_products(mp)
        subch += pyna.Library(rates=[_r])

    return subch

def doit():

    subch = get_library()

    # these are the rates that we are going to allow to be optionally
    # zeroed
    r1 = subch.get_rate("p_c12__n13")
    r2 = subch.get_rate("he4_n13__p_o16")

    # finally, the aprox nets don't include the reverse rates for
    # C12+C12, C12+O16, and O16+O16, so remove those

    for r in subch.get_rates():
        if sorted(r.products) in [[pyna.Nucleus("c12"), pyna.Nucleus("c12")],
                                  [pyna.Nucleus("c12"), pyna.Nucleus("o16")],
                                  [pyna.Nucleus("o16"), pyna.Nucleus("o16")]]:
            subch.remove_rate(r)

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
        _r = subch.get_rate_by_name(r)
        subch.remove_rate(_r)

    # at this point we have a library with all the rates that we want.
    # We can create the network now.

    net = StarKillerCxxNetwork(libraries=[subch], symmetric_screening=False,
                               disable_rate_params=[r1, r2])
    net.make_ap_pg_approx(intermediate_nuclei=["cl35", "k39", "sc43", "v47", "mn51", "co55"])
    net.remove_nuclei(["cl35", "k39", "sc43", "v47", "mn51", "co55"])

    print(f"number of nuclei: {len(net.unique_nuclei)}")
    print(f"number of rates: {len(net.rates)}")

    # We are going to want to replace the reverse rates with rates
    # derived via detailed balance.  We'll do this only for the rates
    # that have reverse, and we'll use the RatePair functionality for
    # this


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
