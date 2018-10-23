Runtime parameter values
========================

.. raw:: latex

   \small

.. table:: BS parameters.

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Maximum number of     | 10000                 |
   |                       | steps to use in the   |                       |
   |    \endfoot           | ODE integration       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``ode_max_steps``     |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``ode_method``        | use an implementation | 1                     |
   |                       | of the Bulirsch-Stoer |                       |
   |                       | semi-implicit         |                       |
   |                       | extrapolation method  |                       |
   |                       | (1) or a Rosenbrock   |                       |
   |                       | method (2)            |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Floor to use for the  | 1.d-6                 |
   |                       | ODE scaling vector    |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``ode_scale_floor``   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``safety_factor``     | when constructing the | 1.d9                  |
   |                       | intermediate steps in |                       |
   |                       | the stiff ODE         |                       |
   |                       | integration by how    |                       |
   |                       | much do we allow the  |                       |
   |                       | state variables to    |                       |
   |                       | change over a dt      |                       |
   |                       | before giving up on   |                       |
   |                       | the step and retrying |                       |
   |                       | with a smaller step?  |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Which choice to use   | 2                     |
   |                       | for the ODE scaling   |                       |
   |    \rowcolor{tableSha | 1:                    |                       |
   | de}                   | :math:`|y| + |dy/dt|` |                       |
   |                       | ;                     |                       |
   | ``scaling_method``    | 2:                    |                       |
   |                       | :math:`\max(|y|, K)`  |                       |
   |                       | with :math:`K =`      |                       |
   |                       | constant              |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``use_timestep_estima | use the VODE          | .false.               |
   | tor``                 | algorithm’s initial   |                       |
   |                       | timestep estimator?   |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: VBDF parameters.

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | minimum allowable     | 1.d-24                |
   |                       | timestep              |                       |
   |    \endfoot           |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``dt_min``            |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``jac_age``           | number of times we    | 50                    |
   |                       | can use the Jacobian  |                       |
   |                       | before rebuilding     |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | number of times we    | 20                    |
   |                       | use the same Newton   |                       |
   |    \rowcolor{tableSha | iteration matrix      |                       |
   | de}                   | before rebuilding     |                       |
   |                       |                       |                       |
   | ``p_age``             |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``reuse_jac``         | reuse the Jacobian?   | .false.               |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: breakout parameters.

   +--------------------------+--+------+
   |                          |  |      |
   +--------------------------+--+------+
   | Table —continued         |  |      |
   +--------------------------+--+------+
   |                          |  |      |
   +--------------------------+--+------+
   |                          |  |      |
   +--------------------------+--+------+
   | .. raw:: latex           |  | 0.d0 |
   |                          |  |      |
   |    \endfoot              |  |      |
   |                          |  |      |
   | .. raw:: latex           |  |      |
   |                          |  |      |
   |    \hline                |  |      |
   |                          |  |      |
   | .. raw:: latex           |  |      |
   |                          |  |      |
   |    \endlastfoot          |  |      |
   |                          |  |      |
   | .. raw:: latex           |  |      |
   |                          |  |      |
   |    \rowcolor{tableShade} |  |      |
   |                          |  |      |
   | ``eos_gamma``            |  |      |
   +--------------------------+--+------+

.. raw:: latex

   \small

.. table:: burn_cell parameters.

   +----------------------------+--+------+
   | [table: burn_cell runtime] |  |      |
   +============================+==+======+
   |                            |  |      |
   +----------------------------+--+------+
   | Table —continued           |  |      |
   +----------------------------+--+------+
   |                            |  |      |
   +----------------------------+--+------+
   |                            |  |      |
   +----------------------------+--+------+
   | .. raw:: latex             |  | ""   |
   |                            |  |      |
   |    \endfoot                |  |      |
   |                            |  |      |
   | .. raw:: latex             |  |      |
   |                            |  |      |
   |    \hline                  |  |      |
   |                            |  |      |
   | .. raw:: latex             |  |      |
   |                            |  |      |
   |    \endlastfoot            |  |      |
   |                            |  |      |
   | .. raw:: latex             |  |      |
   |                            |  |      |
   |    \rowcolor{tableShade}   |  |      |
   |                            |  |      |
   | ``run_prefix``             |  |      |
   +----------------------------+--+------+
   | ``small_dens``             |  | 1.e5 |
   +----------------------------+--+------+
   | .. raw:: latex             |  | 1.e5 |
   |                            |  |      |
   |    \rowcolor{tableShade}   |  |      |
   |                            |  |      |
   | ``small_temp``             |  |      |
   +----------------------------+--+------+

.. raw:: latex

   \small

.. table:: cj_detonation parameters.

   +--------------------------------+--+--------+
   | [table: cj_detonation runtime] |  |        |
   +================================+==+========+
   |                                |  |        |
   +--------------------------------+--+--------+
   | Table —continued               |  |        |
   +--------------------------------+--+--------+
   |                                |  |        |
   +--------------------------------+--+--------+
   |                                |  |        |
   +--------------------------------+--+--------+
   | .. raw:: latex                 |  | 1.e-10 |
   |                                |  |        |
   |    \endfoot                    |  |        |
   |                                |  |        |
   | .. raw:: latex                 |  |        |
   |                                |  |        |
   |    \hline                      |  |        |
   |                                |  |        |
   | .. raw:: latex                 |  |        |
   |                                |  |        |
   |    \endlastfoot                |  |        |
   |                                |  |        |
   | .. raw:: latex                 |  |        |
   |                                |  |        |
   |    \rowcolor{tableShade}       |  |        |
   |                                |  |        |
   | ``smallx``                     |  |        |
   +--------------------------------+--+--------+

.. raw:: latex

   \small

.. table:: gamma_law_general parameters.

   +------------------------------------+--+-----------+
   | [table: gamma_law_general runtime] |  |           |
   +====================================+==+===========+
   |                                    |  |           |
   +------------------------------------+--+-----------+
   | Table —continued                   |  |           |
   +------------------------------------+--+-----------+
   |                                    |  |           |
   +------------------------------------+--+-----------+
   |                                    |  |           |
   +------------------------------------+--+-----------+
   | .. raw:: latex                     |  | .true.    |
   |                                    |  |           |
   |    \endfoot                        |  |           |
   |                                    |  |           |
   | .. raw:: latex                     |  |           |
   |                                    |  |           |
   |    \hline                          |  |           |
   |                                    |  |           |
   | .. raw:: latex                     |  |           |
   |                                    |  |           |
   |    \endlastfoot                    |  |           |
   |                                    |  |           |
   | .. raw:: latex                     |  |           |
   |                                    |  |           |
   |    \rowcolor{tableShade}           |  |           |
   |                                    |  |           |
   | ``eos_assume_neutral``             |  |           |
   +------------------------------------+--+-----------+
   | ``eos_gamma``                      |  | 5.d0/3.d0 |
   +------------------------------------+--+-----------+

.. raw:: latex

   \small

.. table:: helmholtz parameters.

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Force the EOS output  | .false.               |
   |                       | quantities to match   |                       |
   |    \endfoot           | input                 |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``eos_input_is_consta |                       |                       |
   | nt``                  |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``use_eos_coulomb``   | use the Coulomb       | .true.                |
   |                       | corrections           |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: integration parameters.

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | The maximum           | 1.0d11                |
   |                       | temperature for       |                       |
   |    \endfoot           | reactions in the      |                       |
   |                       | integration.          |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``MAX_TEMP``          |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``SMALL_X_SAFE``      | The absolute cutoff   | 1.0d-30               |
   |                       | for species – note    |                       |
   |                       | that this might be    |                       |
   |                       | larger than small_x,  |                       |
   |                       | but the issue is that |                       |
   |                       | we need to prevent    |                       |
   |                       | underflow issues and  |                       |
   |                       | keep mass fractions   |                       |
   |                       | positive in the       |                       |
   |                       | integrator. You may   |                       |
   |                       | have to increase the  |                       |
   |                       | floor to, e.g. 1.d-20 |                       |
   |                       | if your rates are     |                       |
   |                       | large.                |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1.d-6                 |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``atol_enuc``         |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``atol_spec``         |                       | 1.d-12                |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1.d-6                 |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``atol_temp``         |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``burner_verbose``    | Should we print out   | .false.               |
   |                       | diagnostic output     |                       |
   |                       | after the solve?      |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Integration mode: if  | 1                     |
   |                       | 0, a hydrostatic burn |                       |
   |    \rowcolor{tableSha | (temperature and      |                       |
   | de}                   | density remain        |                       |
   |                       | constant), and if 1,  |                       |
   | ``burning_mode``      | a self-heating burn   |                       |
   |                       | (temperature/energy   |                       |
   |                       | evolve with the       |                       |
   |                       | burning). If 2, a     |                       |
   |                       | hybrid approach       |                       |
   |                       | presented by Raskin   |                       |
   |                       | et al. (2010): do     |                       |
   |                       | hydrostatic           |                       |
   |                       | everywhere, but if    |                       |
   |                       | the hydrostatic burn  |                       |
   |                       | gives us a negative   |                       |
   |                       | energy change, redo   |                       |
   |                       | the burn in           |                       |
   |                       | self-heating mode. If |                       |
   |                       | 3, do normal          |                       |
   |                       | self-heating, but     |                       |
   |                       | limit all values of   |                       |
   |                       | the RHS by the same   |                       |
   |                       | factor :math:`L` such |                       |
   |                       | that                  |                       |
   |                       | :math:`\dot{e} = f_s  |                       |
   |                       | e / t_s`,             |                       |
   |                       | where :math:`\dot{e}` |                       |
   |                       | is the energy         |                       |
   |                       | injection rate,       |                       |
   |                       | :math:`e` is the      |                       |
   |                       | internal energy of    |                       |
   |                       | the zone, :math:`t_s` |                       |
   |                       | is the sound crossing |                       |
   |                       | time, and :math:`f_s` |                       |
   |                       | is a safety factor.   |                       |
   |                       | :math:`L` is computed |                       |
   |                       | as min(1,             |                       |
   |                       | :math:`f_s (e / \dot{ |                       |
   |                       | e}) / t_s`).          |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``burning_mode_factor | If we’re using        | 1.d-1                 |
   | ``                    | burning_mode == 3,    |                       |
   |                       | this is the safety    |                       |
   |                       | factor :math:`f_s` to |                       |
   |                       | use.                  |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Do we call the EOS    | .false.               |
   |                       | each time we enter    |                       |
   |    \rowcolor{tableSha | the EOS? This is      |                       |
   | de}                   | expensive, but more   |                       |
   |                       | accurate. Otherwise,  |                       |
   | ``call_eos_in_rhs``   | we instead call the   |                       |
   |                       | EOS at the start of   |                       |
   |                       | the integration and   |                       |
   |                       | freeze the            |                       |
   |                       | thermodynamics        |                       |
   |                       | throughout the RHS    |                       |
   |                       | evalulation. This     |                       |
   |                       | only affects the      |                       |
   |                       | temperature           |                       |
   |                       | integration (which is |                       |
   |                       | the input to the rate |                       |
   |                       | evaluation). In       |                       |
   |                       | particular, since we  |                       |
   |                       | calculate the         |                       |
   |                       | composition factors   |                       |
   |                       | either way, this      |                       |
   |                       | determines whether    |                       |
   |                       | we’re updating the    |                       |
   |                       | thermodynamic         |                       |
   |                       | derivatives and other |                       |
   |                       | quantities (cp and    |                       |
   |                       | cv) as we go.         |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``centered_diff_jac`` | one-sided numerical   | .false.               |
   |                       | jacobian (.False.) or |                       |
   |                       | centered-difference   |                       |
   |                       | Jacobian (.true.).    |                       |
   |                       | Note: the             |                       |
   |                       | centered-difference   |                       |
   |                       | requires twice as     |                       |
   |                       | many RHS calls        |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | If we want to call    | 1.0d20                |
   |                       | the EOS in general,   |                       |
   |    \rowcolor{tableSha | but don’t want to     |                       |
   | de}                   | overdo it, we can set |                       |
   |                       | a fraction dT_crit    |                       |
   | ``dT_crit``           | such that we only do  |                       |
   |                       | the EOS call if the   |                       |
   |                       | temperature has       |                       |
   |                       | changed by a relative |                       |
   |                       | fraction :math:`>`    |                       |
   |                       | dT_crit. If we use    |                       |
   |                       | this option, we will  |                       |
   |                       | do a linear fit to    |                       |
   |                       | c_v and c_p in        |                       |
   |                       | between EOS calls.    |                       |
   |                       | This will work        |                       |
   |                       | regardless of         |                       |
   |                       | call_eos_in_rhs.      |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``do_constant_volume_ | When evolving the     | .false.               |
   | burn``                | temperature, should   |                       |
   |                       | we assume a constant  |                       |
   |                       | pressure (default) or |                       |
   |                       | a constant volume     |                       |
   |                       | (do_constant_volume_b |                       |
   |                       | urn                   |                       |
   |                       | = T)?                 |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Allow the energy      | .true.                |
   |                       | integration to be     |                       |
   |    \rowcolor{tableSha | disabled by setting   |                       |
   | de}                   | the RHS to zero.      |                       |
   |                       |                       |                       |
   | ``integrate_energy``  |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``integrate_temperatu | Allow the temperature | .true.                |
   | re``                  | integration to be     |                       |
   |                       | disabled by setting   |                       |
   |                       | the RHS to zero.      |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Whether to use an     | 1                     |
   |                       | analytical or         |                       |
   |    \rowcolor{tableSha | numerical Jacobian. 1 |                       |
   | de}                   | == Analytical 2 ==    |                       |
   |                       | Numerical             |                       |
   | ``jacobian``          |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``react_boost``       | boost the reaction    | -1.d0                 |
   |                       | rates by a factor > 1 |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Whether to            | .false.               |
   |                       | renormalize the mass  |                       |
   |    \rowcolor{tableSha | fractions at each     |                       |
   | de}                   | step in the evolution |                       |
   |                       | so that they sum to   |                       |
   | ``renormalize_abundan | unity.                |                       |
   | ces``                 |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``retry_burn``        | If we fail to find a  | .false.               |
   |                       | solution consistent   |                       |
   |                       | with the tolerances,  |                       |
   |                       | do we want to try     |                       |
   |                       | again with a looser   |                       |
   |                       | tolerance?            |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | If we do retry a      | 1.25d0                |
   |                       | burn, by what factor  |                       |
   |    \rowcolor{tableSha | should we loosen the  |                       |
   | de}                   | tolerance?            |                       |
   |                       |                       |                       |
   | ``retry_burn_factor`` |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``retry_burn_max_chan | What is the maximum   | 1.0d2                 |
   | ge``                  | factor we can         |                       |
   |                       | increase the original |                       |
   |                       | tolerances by?        |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1.d-6                 |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``rtol_enuc``         |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``rtol_spec``         | Tolerances for the    | 1.d-12                |
   |                       | solver (relative and  |                       |
   |                       | absolute), for the    |                       |
   |                       | species, temperature, |                       |
   |                       | and energy equations. |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1.d-6                 |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``rtol_temp``         |                       |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: kpp parameters.

   +--------------------------+--+-------+
   |                          |  |       |
   +--------------------------+--+-------+
   | Table —continued         |  |       |
   +--------------------------+--+-------+
   |                          |  |       |
   +--------------------------+--+-------+
   |                          |  |       |
   +--------------------------+--+-------+
   | .. raw:: latex           |  | 10.d0 |
   |                          |  |       |
   |    \endfoot              |  |       |
   |                          |  |       |
   | .. raw:: latex           |  |       |
   |                          |  |       |
   |    \hline                |  |       |
   |                          |  |       |
   | .. raw:: latex           |  |       |
   |                          |  |       |
   |    \endlastfoot          |  |       |
   |                          |  |       |
   | .. raw:: latex           |  |       |
   |                          |  |       |
   |    \rowcolor{tableShade} |  |       |
   |                          |  |       |
   | ``A_burn``               |  |       |
   +--------------------------+--+-------+

.. raw:: latex

   \small

.. table:: multigamma parameters.

   +--------------------------+--+-----+
   |                          |  |     |
   +--------------------------+--+-----+
   | Table —continued         |  |     |
   +--------------------------+--+-----+
   |                          |  |     |
   +--------------------------+--+-----+
   |                          |  |     |
   +--------------------------+--+-----+
   | .. raw:: latex           |  | 1.4 |
   |                          |  |     |
   |    \endfoot              |  |     |
   |                          |  |     |
   | .. raw:: latex           |  |     |
   |                          |  |     |
   |    \hline                |  |     |
   |                          |  |     |
   | .. raw:: latex           |  |     |
   |                          |  |     |
   |    \endlastfoot          |  |     |
   |                          |  |     |
   | .. raw:: latex           |  |     |
   |                          |  |     |
   |    \rowcolor{tableShade} |  |     |
   |                          |  |     |
   | ``eos_gamma_default``    |  |     |
   +--------------------------+--+-----+
   | ``species_a_gamma``      |  | 1.4 |
   +--------------------------+--+-----+
   | .. raw:: latex           |  | ""  |
   |                          |  |     |
   |    \rowcolor{tableShade} |  |     |
   |                          |  |     |
   | ``species_a_name``       |  |     |
   +--------------------------+--+-----+
   | ``species_b_gamma``      |  | 1.4 |
   +--------------------------+--+-----+
   | .. raw:: latex           |  | ""  |
   |                          |  |     |
   |    \rowcolor{tableShade} |  |     |
   |                          |  |     |
   | ``species_b_name``       |  |     |
   +--------------------------+--+-----+
   | ``species_c_gamma``      |  | 1.4 |
   +--------------------------+--+-----+
   | .. raw:: latex           |  | ""  |
   |                          |  |     |
   |    \rowcolor{tableShade} |  |     |
   |                          |  |     |
   | ``species_c_name``       |  |     |
   +--------------------------+--+-----+

.. raw:: latex

   \small

.. table:: networks parameters.

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | cutoff for species    | 1.d-30                |
   |                       | mass fractions        |                       |
   |    \endfoot           |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``small_x``           |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``use_c12ag_deboer17` | Should we use Deboer  | .false.               |
   | `                     | + 2017 rate for       |                       |
   |                       | c12(a,g)o16?          |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Should we use rate    | .false.               |
   |                       | tables if they are    |                       |
   |    \rowcolor{tableSha | present in the        |                       |
   | de}                   | network?              |                       |
   |                       |                       |                       |
   | ``use_tables``        |                       |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: polytrope parameters.

   +--------------------------+--+-------+
   |                          |  |       |
   +--------------------------+--+-------+
   | Table —continued         |  |       |
   +--------------------------+--+-------+
   |                          |  |       |
   +--------------------------+--+-------+
   |                          |  |       |
   +--------------------------+--+-------+
   | .. raw:: latex           |  | 0.0d0 |
   |                          |  |       |
   |    \endfoot              |  |       |
   |                          |  |       |
   | .. raw:: latex           |  |       |
   |                          |  |       |
   |    \hline                |  |       |
   |                          |  |       |
   | .. raw:: latex           |  |       |
   |                          |  |       |
   |    \endlastfoot          |  |       |
   |                          |  |       |
   | .. raw:: latex           |  |       |
   |                          |  |       |
   |    \rowcolor{tableShade} |  |       |
   |                          |  |       |
   | ``polytrope_K``          |  |       |
   +--------------------------+--+-------+
   | ``polytrope_gamma``      |  | 0.0d0 |
   +--------------------------+--+-------+
   | .. raw:: latex           |  | 2.0d0 |
   |                          |  |       |
   |    \rowcolor{tableShade} |  |       |
   |                          |  |       |
   | ``polytrope_mu_e``       |  |       |
   +--------------------------+--+-------+
   | ``polytrope_type``       |  | 0     |
   +--------------------------+--+-------+

.. raw:: latex

   \small

.. table:: powerlaw parameters.

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | reaction thresholds   | 1.0d0                 |
   |                       | (for the power law)   |                       |
   |    \endfoot           |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``T_burn_ref``        |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``burning_mode``      | override the default  | 0                     |
   |                       | burning mode with a   |                       |
   |                       | higher priority       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1.0d0                 |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``f_act``             |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``jacobian``          | override the default  | 2                     |
   |                       | Jacobian mode with a  |                       |
   |                       | higher priority       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | exponent for the      | 4.d0                  |
   |                       | temperature           |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``nu``                |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``rho_burn_ref``      |                       | 1.0d0                 |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the coefficient for   | 1.d0                  |
   |                       | the reaction rate     |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``rtilde``            |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``specific_q_burn``   | reaction specific     | 10.d0                 |
   |                       | q-value (in erg/g)    |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: rprox parameters.

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1.0e-8                |
   |                       |                       |                       |
   |    \endfoot           |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``atol_enuc``         |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``atol_spec``         | override the default  | 1.0e-11               |
   |                       | tolerances for        |                       |
   |                       | backwards             |                       |
   |                       | compatibility         |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1.0e-8                |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``atol_temp``         |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``burning_mode``      | override the default  | 1                     |
   |                       | burning mode with a   |                       |
   |                       | higher priority       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | we typically run this | .false.               |
   |                       | network for           |                       |
   |    \rowcolor{tableSha | constant-pressure     |                       |
   | de}                   | burns                 |                       |
   |                       |                       |                       |
   | ``do_constant_volume_ |                       |                       |
   | burn``                |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``jacobian``          | override so that the  | 1                     |
   |                       | default is an         |                       |
   |                       | analytical Jacobian   |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1.0e-8                |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``rtol_enuc``         |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``rtol_spec``         |                       | 1.0e-12               |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1.0e-8                |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``rtol_temp``         |                       |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: stellarcollapse parameters.

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | name of the HDF5 file | ""                    |
   |                       | containing tabulated  |                       |
   |    \endfoot           | data                  |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``eos_file``          |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``use_energy_shift``  |                       | .false.               |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: test_eos parameters.

   +---------------------------+--+-----------+
   | [table: test_eos runtime] |  |           |
   +===========================+==+===========+
   |                           |  |           |
   +---------------------------+--+-----------+
   | Table —continued          |  |           |
   +---------------------------+--+-----------+
   |                           |  |           |
   +---------------------------+--+-----------+
   |                           |  |           |
   +---------------------------+--+-----------+
   | .. raw:: latex            |  | 1.d9      |
   |                           |  |           |
   |    \endfoot               |  |           |
   |                           |  |           |
   | .. raw:: latex            |  |           |
   |                           |  |           |
   |    \hline                 |  |           |
   |                           |  |           |
   | .. raw:: latex            |  |           |
   |                           |  |           |
   |    \endlastfoot           |  |           |
   |                           |  |           |
   | .. raw:: latex            |  |           |
   |                           |  |           |
   |    \rowcolor{tableShade}  |  |           |
   |                           |  |           |
   | ``dens_max``              |  |           |
   +---------------------------+--+-----------+
   | ``dens_min``              |  | 1.d6      |
   +---------------------------+--+-----------+
   | .. raw:: latex            |  | 0.1d0     |
   |                           |  |           |
   |    \rowcolor{tableShade}  |  |           |
   |                           |  |           |
   | ``metalicity_max``        |  |           |
   +---------------------------+--+-----------+
   | ``small_dens``            |  | 1.e-4     |
   +---------------------------+--+-----------+
   | .. raw:: latex            |  | 1.e4      |
   |                           |  |           |
   |    \rowcolor{tableShade}  |  |           |
   |                           |  |           |
   | ``small_temp``            |  |           |
   +---------------------------+--+-----------+
   | ``temp_max``              |  | 1.d12     |
   +---------------------------+--+-----------+
   | .. raw:: latex            |  | 1.d6      |
   |                           |  |           |
   |    \rowcolor{tableShade}  |  |           |
   |                           |  |           |
   | ``temp_min``              |  |           |
   +---------------------------+--+-----------+
   | ``test_set``              |  | "gr0_3d"  |
   +---------------------------+--+-----------+
   | .. raw:: latex            |  | "uniform" |
   |                           |  |           |
   |    \rowcolor{tableShade}  |  |           |
   |                           |  |           |
   | ``xin_file``              |  |           |
   +---------------------------+--+-----------+

.. raw:: latex

   \small

.. table:: test_react parameters.

   +-----------------------------+--+-----------+
   | [table: test_react runtime] |  |           |
   +=============================+==+===========+
   |                             |  |           |
   +-----------------------------+--+-----------+
   | Table —continued            |  |           |
   +-----------------------------+--+-----------+
   |                             |  |           |
   +-----------------------------+--+-----------+
   |                             |  |           |
   +-----------------------------+--+-----------+
   | .. raw:: latex              |  | 1.d9      |
   |                             |  |           |
   |    \endfoot                 |  |           |
   |                             |  |           |
   | .. raw:: latex              |  |           |
   |                             |  |           |
   |    \hline                   |  |           |
   |                             |  |           |
   | .. raw:: latex              |  |           |
   |                             |  |           |
   |    \endlastfoot             |  |           |
   |                             |  |           |
   | .. raw:: latex              |  |           |
   |                             |  |           |
   |    \rowcolor{tableShade}    |  |           |
   |                             |  |           |
   | ``dens_max``                |  |           |
   +-----------------------------+--+-----------+
   | ``dens_min``                |  | 1.d6      |
   +-----------------------------+--+-----------+
   | .. raw:: latex              |  | 1         |
   |                             |  |           |
   |    \rowcolor{tableShade}    |  |           |
   |                             |  |           |
   | ``do_acc``                  |  |           |
   +-----------------------------+--+-----------+
   | ``run_prefix``              |  | ""        |
   +-----------------------------+--+-----------+
   | .. raw:: latex              |  | 1.e5      |
   |                             |  |           |
   |    \rowcolor{tableShade}    |  |           |
   |                             |  |           |
   | ``small_dens``              |  |           |
   +-----------------------------+--+-----------+
   | ``small_temp``              |  | 1.e5      |
   +-----------------------------+--+-----------+
   | .. raw:: latex              |  | 1.d15     |
   |                             |  |           |
   |    \rowcolor{tableShade}    |  |           |
   |                             |  |           |
   | ``temp_max``                |  |           |
   +-----------------------------+--+-----------+
   | ``temp_min``                |  | 1.d6      |
   +-----------------------------+--+-----------+
   | .. raw:: latex              |  | "gr0_3d"  |
   |                             |  |           |
   |    \rowcolor{tableShade}    |  |           |
   |                             |  |           |
   | ``test_set``                |  |           |
   +-----------------------------+--+-----------+
   | ``tmax``                    |  | 0.1d0     |
   +-----------------------------+--+-----------+
   | .. raw:: latex              |  | "uniform" |
   |                             |  |           |
   |    \rowcolor{tableShade}    |  |           |
   |                             |  |           |
   | ``xin_file``                |  |           |
   +-----------------------------+--+-----------+

.. raw:: latex

   \small

.. table:: test_sdc parameters.

   +---------------------------+--+-----------+
   | [table: test_sdc runtime] |  |           |
   +===========================+==+===========+
   |                           |  |           |
   +---------------------------+--+-----------+
   | Table —continued          |  |           |
   +---------------------------+--+-----------+
   |                           |  |           |
   +---------------------------+--+-----------+
   |                           |  |           |
   +---------------------------+--+-----------+
   | .. raw:: latex            |  | 1.d9      |
   |                           |  |           |
   |    \endfoot               |  |           |
   |                           |  |           |
   | .. raw:: latex            |  |           |
   |                           |  |           |
   |    \hline                 |  |           |
   |                           |  |           |
   | .. raw:: latex            |  |           |
   |                           |  |           |
   |    \endlastfoot           |  |           |
   |                           |  |           |
   | .. raw:: latex            |  |           |
   |                           |  |           |
   |    \rowcolor{tableShade}  |  |           |
   |                           |  |           |
   | ``dens_max``              |  |           |
   +---------------------------+--+-----------+
   | ``dens_min``              |  | 1.d6      |
   +---------------------------+--+-----------+
   | .. raw:: latex            |  | 1         |
   |                           |  |           |
   |    \rowcolor{tableShade}  |  |           |
   |                           |  |           |
   | ``do_acc``                |  |           |
   +---------------------------+--+-----------+
   | ``run_prefix``            |  | ""        |
   +---------------------------+--+-----------+
   | .. raw:: latex            |  | 1.e5      |
   |                           |  |           |
   |    \rowcolor{tableShade}  |  |           |
   |                           |  |           |
   | ``small_dens``            |  |           |
   +---------------------------+--+-----------+
   | ``small_temp``            |  | 1.e5      |
   +---------------------------+--+-----------+
   | .. raw:: latex            |  | 1.d15     |
   |                           |  |           |
   |    \rowcolor{tableShade}  |  |           |
   |                           |  |           |
   | ``temp_max``              |  |           |
   +---------------------------+--+-----------+
   | ``temp_min``              |  | 1.d6      |
   +---------------------------+--+-----------+
   | .. raw:: latex            |  | "gr0_3d"  |
   |                           |  |           |
   |    \rowcolor{tableShade}  |  |           |
   |                           |  |           |
   | ``test_set``              |  |           |
   +---------------------------+--+-----------+
   | ``tmax``                  |  | 0.1d0     |
   +---------------------------+--+-----------+
   | .. raw:: latex            |  | "uniform" |
   |                           |  |           |
   |    \rowcolor{tableShade}  |  |           |
   |                           |  |           |
   | ``xin_file``              |  |           |
   +---------------------------+--+-----------+

.. raw:: latex

   \small

.. table:: triple_alpha_plus_cago parameters.

   +-----------------------+-----------------------+-----------------------+
   | [table: triple_alpha_ |                       |                       |
   | plus_cago runtime]    |                       |                       |
   +=======================+=======================+=======================+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1.0e-8                |
   |                       |                       |                       |
   |    \endfoot           |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``atol_enuc``         |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``atol_spec``         | override the default  | 1.0e-12               |
   |                       | tolerances for        |                       |
   |                       | backwards             |                       |
   |                       | compatibility         |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1.0e-8                |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``atol_temp``         |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``burning_mode``      | override the default  | 1                     |
   |                       | burning mode with a   |                       |
   |                       | higher priority       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | we typically run this | .false.               |
   |                       | network for           |                       |
   |    \rowcolor{tableSha | constant-pressure     |                       |
   | de}                   | burns                 |                       |
   |                       |                       |                       |
   | ``do_constant_volume_ |                       |                       |
   | burn``                |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``jacobian``          | override so that the  | 1                     |
   |                       | default is an         |                       |
   |                       | analytical Jacobian   |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1.0e-6                |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``rtol_enuc``         |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``rtol_spec``         |                       | 1.0e-12               |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1.0e-6                |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``rtol_temp``         |                       |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: xrb_simple parameters.

   +-----------------------+-----------------------+-----------------------+
   | [table: xrb_simple ru |                       |                       |
   | ntime]                |                       |                       |
   +=======================+=======================+=======================+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1.0e-8                |
   |                       |                       |                       |
   |    \endfoot           |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``atol_enuc``         |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``atol_spec``         | override the default  | 1.0e-11               |
   |                       | tolerances for        |                       |
   |                       | backwards             |                       |
   |                       | compatibility         |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1.0e-8                |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``atol_temp``         |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``burning_mode``      | override the default  | 1                     |
   |                       | burning mode with a   |                       |
   |                       | higher priority       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | we typically run this | .false.               |
   |                       | network for           |                       |
   |    \rowcolor{tableSha | constant-pressure     |                       |
   | de}                   | burns                 |                       |
   |                       |                       |                       |
   | ``do_constant_volume_ |                       |                       |
   | burn``                |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``jacobian``          | override so that the  | 2                     |
   |                       | default is a          |                       |
   |                       | numerical Jacobian;   |                       |
   |                       | we don’t yet have an  |                       |
   |                       | analytical Jacobian   |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1.0e-8                |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``rtol_enuc``         |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | ``rtol_spec``         |                       | 1.0e-12               |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1.0e-8                |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   |                       |                       |                       |
   | ``rtol_temp``         |                       |                       |
   +-----------------------+-----------------------+-----------------------+
