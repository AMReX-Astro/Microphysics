# Using the python interface to Microphysics

To build the python library for this, go to the
`Microphysics/python_library` directory, and do e.g. `make
NETWORK_DIR=subch`. This should create a *.so library file and print
"SUCCESS."

Add `Microphysics/python_library` to your `PYTHONPATH` environment
variable and copy `Microphysics/python_library/helm_table.dat` to this
directory after building the library.

Also copy (or symlink) any *.dat files from your network directory
containing tabulated or Reaclib rates (if using a pynucastro network).