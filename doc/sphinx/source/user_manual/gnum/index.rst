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

Python
^^^^^^

.. autoclass:: Pypdm.Pypdm.GlobalNumbering
  :members:

The following example shows how to build a global numbering from a set of geometric coordinates.

.. code:: python

   import mpi4py.MPI as MPI
   import Pypdm.Pypdm as PDM

   # First, create a GlobalNumbering instance and set some parameters
   gen_gnum = PDM.GlobalNumbering(dim,
                                  n_part,
                                  merge,
                                  tolerance,
                                  MPI.COMM_WORLD)

   # Then, provide the coordinates array for each partition
   # (coords is a list of numpy arrays of type double)
   # (here we omit the optional char_length argument)
   for ipart in range(n_part):
     gen_gnum.set_from_coords(ipart,
                              coords[ipart],
                              None)

   # Once all partitions have been set, build the global numbering
   gen_gnum.compute()

   # Finally, retrieve the computed global id arrays
   for ipart in range(n_part):
     gnum[ipart] = gen_gnum.get(ipart)


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
