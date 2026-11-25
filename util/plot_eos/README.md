# `plot_eos`

A simple driver that outputs EOS data as rows of format:

```
 T  rho  var
```

where `var` is a quantity from the EOS.

This can then be used for plotting.


## Thermodynamic conditions

The number of density and temperature points and their limits can be
set as parameters.

If an EOS quantity is not provided, you can add it by editing
the `string_to_member` function.

The composition is uniformly set across all species mass fractions.


## Plotting

### GNUplot

For GNUplot, run as:

```
./main3d.gnu.ex plot.skip_after_row=1 > eosdata.out
```

and then you can plot via:

```
plot "eosdata.out" using 1:2:3 with image
```

### python / matplotlib
