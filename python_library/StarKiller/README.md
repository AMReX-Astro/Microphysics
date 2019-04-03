# StarKiller

There are two components to the python interfaces - this StarKiller
python package and the compiled python library that wraps the fortran
subroutines.

This is developed and tested with Python 3.

## Step 1: Install the StarKiller python package

```
python3 setup.py install --user
```

This only needs to be done once.

## Step 2: Compile the Fortran and Python wrappers

To build the python library for this, go to the
`Microphysics/python_library` directory, and `make`. This should
create a *.so library file and print "SUCCESS."

The compilation defaults are to use aprox13 with the Helmholtz
EOS. The network and EOS can be specified on the `make` line as usual
for Microphysics compilation using `make NETWORK_DIR=ABC EOS_DIR=XYZ`.

For now, this needs to be redone any time a different network or EOS
is desired.

Add the `Microphysics/python_library` path to your `PYTHONPATH`
environment variable so python can find the library file.

## Step 3: Copy files to the working directory

- Copy `Microphysics/python_library/helm_table.dat` to your current
  working directory after building the library.

- Also copy (or symlink) any *.dat files from your network directory
  to your working directory. These files contain tabulated or Reaclib
  rates (if using a pynucastro network).

- Copy the desired "probin" initialization file to your working
  directory. The default is "StarKiller/examples/probin_aprox13".

## Step 4: Try a test problem

As an example of the EOS interface, try running `python3 call_eos.py -h`
to see input options.

The `call_eos.py` script is in the
`StarKiller/examples` directory.
