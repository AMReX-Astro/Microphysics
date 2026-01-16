import pynucastro as pyna

reaclib_library = pyna.ReacLibLibrary()
fwd_reactions = reaclib_library.forward_for_detailed_balance()

nuclei = ["p", "he4", "fe52", "ni56", "co55"]

fwd_rates_lib = fwd_reactions.linking_nuclei(nuclist=nuclei,
                                             with_reverse=False)

derived = []
for r in fwd_rates_lib.get_rates():
    d = pyna.DerivedRate(r, use_pf=True)
    derived.append(d)

der_rates_lib = pyna.Library(rates=derived)
full_lib = fwd_rates_lib + der_rates_lib
net = pyna.AmrexAstroCxxNetwork(libraries=[full_lib])

net.write_network()
