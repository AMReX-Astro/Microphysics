# create the core library used by the he-burn group of networks
# they can then adjust these via various approximations

import pynucastro as pyna
from pynucastro.rates import ReacLibRate


def get_core_library(*,
                     include_n14_sequence=False,
                     include_zn=False,
                     include_iron_peak=False,
                     include_low_ye=False,
                     do_detailed_balance=False,
                     extra_nuclei=None):

    reaclib_lib = pyna.ReacLibLibrary()

    if extra_nuclei is None:
        extra_nuclei = []

    nuclei = ["p",
              "he4", "c12", "n13", "o16",
              "ne20", "na23", "mg24", "si28", "s32",
              "ar36", "ca40", "ti44", "cr48",
              "fe52", "ni56",
              "al27", "p31", "cl35", "k39", "sc43", "v47",
              "mn51", "co55"]

    nuclei += extra_nuclei

    if include_n14_sequence:
        nuclei += ["n14", "f18", "ne21", "na22"]

    if include_zn:
        nuclei += ["cu59", "zn60"]

    core_lib = reaclib_lib.linking_nuclei(nuclei)

    # in this list, we have the reactants, the actual reactants,
    # and modified products that we will use instead
    other_rates = [("c12(c12,n)mg23", "mg24"),
                   ("o16(o16,n)s31", "s32"),
                   ("o16(c12,n)si27", "si28")]

    for r, mp in other_rates:
        _r = reaclib_lib.get_rate_by_name(r)
        _r.modify_products(mp)
        core_lib += pyna.Library(rates=[_r])

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

    iron_peak = []
    all_lib = core_lib

    if include_iron_peak:
        # now create a list of iron group nuclei and find both the
        # ReacLib and weak / tabular rates linking these.

        iron_peak = ["n", "p", "he4",
                     "mn51",
                     "fe52", "fe53", "fe54", "fe55", "fe56",
                     "co55", "co56", "co57",
                     "ni56", "ni57", "ni58"]

        if include_zn:
            iron_peak += ["cu59", "zn60"]

        if include_low_ye:
            iron_peak += ["mn55"]


        all_lib += reaclib_lib.linking_nuclei(iron_peak)

    weak_lib = pyna.TabularLibrary()
    iron_weak_lib = weak_lib.linking_nuclei(iron_peak + nuclei)
    all_lib += iron_weak_lib

    if do_detailed_balance:
        rates_to_derive = core_lib.backward().get_rates()

        # now for each of those derived rates, look to see if the pair exists

        for r in rates_to_derive:
            fr = core_lib.get_rate_by_nuclei(r.products, r.reactants)
            if fr:
                print(f"modifying {r} from {fr}")
                core_lib.remove_rate(r)
                d = pyna.DerivedRate(rate=fr, compute_Q=False, use_pf=True)
                core_lib.add_rate(d)

    # we may have duplicate rates -- we want to remove any ReacLib rates
    # that we have tabular rates for

    dupes = all_lib.find_duplicate_links()

    rates_to_remove = []
    for d in dupes:
        for r in d:
            if isinstance(r, ReacLibRate):
                rates_to_remove.append(r)

    for r in rates_to_remove:
        all_lib.remove_rate(r)

    return all_lib
