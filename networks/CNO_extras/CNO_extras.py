import pynucastro as pyna

rl = pyna.ReacLibLibrary()

h_burn = rl.linking_nuclei(["h1", "he4",
                            "c12", "c13",
                            "n13", "n14", "n15",
                            "o14", "o15", "o16","o17","o18",
                            "f17", "f18","f19",
                            "ne18", "ne19", "ne20",
                            "mg22", "mg24"],
                           with_reverse=False)


rc = pyna.StarKillerCxxNetwork(libraries=[h_burn], inert_nuclei=["fe56"])

rc.write_network()

comp = pyna.Composition(rc.get_nuclei())
comp.set_solar_like()

rc.plot(outfile="cno_extras.png", rho=1.e6, T=1.e8, comp=comp, Z_range=[1,13], N_range=[1,13])
rc.plot(outfile="cno_extras_hide_alpha.png", rho=1.e6, T=1.e8, comp=comp, Z_range=[1,13], N_range=[1,13],
                    rotated=True, highlight_filter_function=lambda r: r.Q > 0,
                    curved_edges=True, hide_xalpha=True)
