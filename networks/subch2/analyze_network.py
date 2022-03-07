import pynucastro as pyna

import subch2

subch_library = subch2.get_subch2_library()

# make a network chart based on the conditions inbetween the He/C layer

rho = 1.6e6
T = 1.7e9

rc = pyna.RateCollection(libraries=[subch_library])

comp = pyna.Composition(rc.get_nuclei())
comp.set_all(0.1)
comp.normalize()

print(f"number of nuclei = {len(comp.X)}")
print(f"number of rates = {len(rc.rates)}")

rc.plot_network_chart(rho=rho, T=T, comp=comp,
                      size=(800,2400), outfile="subch_network_chart.pdf")
