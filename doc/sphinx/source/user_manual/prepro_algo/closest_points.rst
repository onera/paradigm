.. _closest_points:

Closest points
==============

C API
-----

Initialization
""""""""""""""

.. doxygenfunction:: PDM_closest_points_create

Point clouds definition
"""""""""""""""""""""""

.. doxygenfunction:: PDM_closest_points_n_part_cloud_set

.. doxygenfunction:: PDM_closest_points_src_cloud_set

.. doxygenfunction:: PDM_closest_points_tgt_cloud_set


Computation
"""""""""""

.. doxygenfunction:: PDM_closest_points_compute

.. doxygenfunction:: PDM_closest_points_dump_times


Results
"""""""

.. doxygenfunction:: PDM_closest_points_part_to_part_get

.. doxygenfunction:: PDM_closest_points_get


Finalization
""""""""""""

.. doxygenfunction:: PDM_closest_points_free



Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  .. todo::
    ...

.. ifconfig:: enable_fortran_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)


Python API
----------

.. ifconfig:: enable_python_doc == 'ON'

  Initialization
  """"""""""""""

  .. autoclass:: Pypdm.Pypdm.ClosestPoints

  Point clouds definition
  """""""""""""""""""""""

  .. autofunction:: Pypdm.Pypdm.ClosestPoints.n_part_cloud_set

  .. autofunction:: Pypdm.Pypdm.ClosestPoints.src_cloud_set

  .. autofunction:: Pypdm.Pypdm.ClosestPoints.tgt_cloud_set


  Computation
  """""""""""

  .. autofunction:: Pypdm.Pypdm.ClosestPoints.compute

  .. autofunction:: Pypdm.Pypdm.ClosestPoints.dump_times


  Results
  """""""

  .. autofunction:: Pypdm.Pypdm.ClosestPoints.part_to_part_get

  .. autofunction:: Pypdm.Pypdm.ClosestPoints.points_get

.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)
