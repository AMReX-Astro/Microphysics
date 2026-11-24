Compute the Chapman-Jouguet detonation speed.

The thermodynamic state for the fuel and ash are hardcoded
into main.cpp right now.

You can set the network and EOS in the `GNUmakefile`.

At the moment, this seems to work well with the `aprox13` network.
It does not seem as robust with larger nets, so some work needs to
be done there.


```
main3d.gnu.ex inputs
```

This will solve for the CJ speed and then try to compute points on
the Hugoniot curve -- it will eventually give an error saying it
doesn't converge.

The output (`hugoniot.txt`) can then be plotted using the `cj_plot.py`
script.

