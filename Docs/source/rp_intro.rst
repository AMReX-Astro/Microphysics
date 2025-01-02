******************
Runtime parameters
******************

The behavior of the network and EOS are controlled by many runtime
parameters. These parameters are defined in plain-text files
``_parameters`` located in the different directories that hold the
microphysics code. At compile time, a script in the AMReX build
system, findparams.py, locates all of the ``_parameters`` files that
are needed for the given choice of network, integrator, and EOS, and
assembles all of the runtime parameters into a set of header files
(using the ``write_probin.py`` script).


Parameter definitions take the form::

    # comment describing the parameter
    name              data-type       default-value      priority

Here, the priority is simply an integer. When two directories
define the same parameter, but with different defaults, the version of
the parameter with the highest priority takes precedence. This allows
specific implementations to override the general parameter defaults.

The documentation below is automatically generated, using the comments
in the ``_parameters`` files.  The parameters are grouped by the
namespace under which they live, and parameters that only apply to
specific build configurations or have their defaults overridden are
noted in separate tables.


.. toctree::

   runtime_parameters

