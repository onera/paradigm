.. _gnum:


================
Global numbering
================


.. _gen_gnum:

Global numbering generation
---------------------------

**ParaDiGM** relies heavily on global ids.
If your code does not use global ids, these can be generated from geometric data.
This global numbering is achieved by encoding Cartesian coordinates along the Morton space-filling curve.
.. Alternatively, from parents & nuplets...

.. Either way, the first step consists in creating an instance of ``PDM_gen_gnum_t`` (or :class:`Pypdm.Pypdm.GlobalNumbering` in Python).

.. .. code:: c

..   PDM_gen_gnum_t *gen_gnum = PDM_gnum_create(dim,
..                                              n_part,
..                                              merge,
..                                              tolerance,
..                                              comm,
..                                              owner);

.. The ``dim``, ``merge`` and ``tolerance`` arguments are only relevant if you want the global numbering to be based on geometric data.

*TODO: show guided example...*






.. doxygenfile:: pdm_gnum.h
   :project: paradigm



.. _gnum_location:

Location from global ids
------------------------

.. doxygenfile:: pdm_gnum_location.h
   :project: paradigm



.. _global_reduction:

Global reduction operations
---------------------------

.. _global_mean:

Global mean
^^^^^^^^^^^

.. doxygenfile:: pdm_global_mean.h
   :project: paradigm


.. _global_reduce:

Global reduction
^^^^^^^^^^^^^^^^

.. doxygenfile:: pdm_global_reduce.h
   :project: paradigm
