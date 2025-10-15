******************
Runtime Parameters
******************

Defining a Parameter
====================

The behavior of the network and EOS are controlled by many runtime
parameters. These parameters are defined in plain-text files
``_parameters`` located in the different directories that hold the
microphysics code. At compile time, a script in the AMReX build
system, ``findparams.py``, locates all of the ``_parameters`` files
that are needed for the given choice of network, integrator, and EOS,
and other modules, and assembles all of the runtime parameters into a
set of header files (using the ``write_probin.py`` script).

A ``_parameter`` file starts with a namespace declaration like:

.. code::

   @namespace: integrator

Followed by a list of parameters.  Each parameter defined after the namespace
declaration would be prefixed by the namespace in an inputs file or commandline,
like ``integrator.parameter``.

The individual parameter definitions take the form::

    # comment describing the parameter
    name              data-type       default-value      priority

where "data-type" can be: ``real``, ``bool``, ``int``, ``string``.

The "priority" is simply an integer that is used to determine
what happens when two different ``_paramerter`` files define the same
parameter but with different defaults.  In this case, the version of
the parameter with the highest priority takes precedence. This allows
specific implementations to override the general parameter defaults.

Known Parameters
================

The documentation below is automatically generated, using the comments
in the ``_parameters`` files.  The parameters are grouped by the
namespace under which they live, and parameters that only apply to
specific build configurations or have their defaults overridden are
noted in separate tables.


.. toctree::

   runtime_parameters
