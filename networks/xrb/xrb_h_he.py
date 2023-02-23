import pynucastro as pyna

rl = pyna.ReacLibLibrary()

h_burn = rl.linking_nuclei(["h1", "he4",
                            "c12", "c13",
                            "n13", "n14", "n15",
                            "o14", "o15", "o16",
                            "f17", "f18",
                            "ne18", "ne19", "ne20", "ne21", "ne22",
                            "na20", "na21", "na22",
                            "mg21", "mg22", "mg23", "mg24"],
                           with_reverse=False)


rc = pyna.StarKillerCxxNetwork(libraries=[h_burn], inert_nuclei=["fe56"])

rc.write_network()

comp = pyna.Composition(rc.get_nuclei())
comp.set_solar_like()

rc.plot(outfile="xrb_h_he.png", rho=1.e6, T=1.e8, comp=comp)
