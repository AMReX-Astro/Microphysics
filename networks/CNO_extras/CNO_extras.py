import pynucastro as pyna

rl = pyna.ReacLibLibrary()

h_burn = rl.linking_nuclei(["h1", "he4",
                            "c12", "c13",
                            "n13", "n14", "n15",
                            "o14", "o15", "o16","o17","o18",
                            "f17", "f18","f19",
                            "ne18", "ne19", "ne20",
                            "mg22", "mg24"],
                           with_reverse=True)


rc = pyna.AmrexAstroCxxNetwork(libraries=[h_burn], inert_nuclei=["fe56"])

rc.write_network()

comp = pyna.Composition(rc.get_nuclei())
comp.set_solar_like()

rho = 1.e6
T = 1.e8

rc.plot(rho, T, comp, outfile="cno_extras.png",
        Z_range=[1, 13], N_range=[1, 13])
rc.plot(outfile="cno_extras_hide_alpha.png",
        Z_range=[1, 13], N_range=[1, 13],
        rotated=True,
        hide_xalpha=True)
