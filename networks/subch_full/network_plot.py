# C-burning with A=23 URCA rate module generator

import pynucastro as pyna
import subch2

subch_lib = subch2.get_library()

rc = pyna.RateCollection(libraries=[subch_lib])

comp = pyna.Composition(rc.get_nuclei())
comp.set_all(0.1)
comp.set_nuc("he4", 0.95)
comp.normalize()

rc.plot(outfile="subch2.pdf", rho=1.e6, T=1.e9, comp=comp, hide_xalpha=True,
        size=(1500, 650), node_size=500, node_font_size=11, node_color="#337dff", node_shape="s",
        Z_range=(1,29), rotated=True, curved_edges=True)
