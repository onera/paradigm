.. _prepro_algo:

###################################
Pre-/co-/post-processing algorithms
###################################


.. .. include:: mesh_location.rst

Point cloud localization inside a mesh
======================================

C API
-----

.. .. doxygenfile:: pdm_mesh_location.h

Initialization
""""""""""""""

.. doxygenfunction:: PDM_mesh_location_create

Source mesh definition
""""""""""""""""""""""

.. doxygenfunction:: PDM_mesh_location_mesh_n_part_set

.. doxygenfunction:: PDM_mesh_location_part_set

.. doxygenfunction:: PDM_mesh_location_nodal_part_set

.. doxygenfunction:: PDM_mesh_location_part_set_2d

.. doxygenfunction:: PDM_mesh_location_nodal_part_set_2d

Target point clouds definition
""""""""""""""""""""""""""""""

.. doxygenfunction:: PDM_mesh_location_n_part_cloud_set

.. doxygenfunction:: PDM_mesh_location_cloud_set

Location computation
""""""""""""""""""""

.. doxygenfunction:: PDM_mesh_location_method_set

.. doxygenfunction:: PDM_mesh_location_tolerance_set

.. doxygenfunction:: PDM_mesh_location_compute

Results
"""""""
.. doxygenfunction:: PDM_mesh_location_n_located_get

.. doxygenfunction:: PDM_mesh_location_located_get

.. doxygenfunction:: PDM_mesh_location_n_unlocated_get

.. doxygenfunction:: PDM_mesh_location_unlocated_get

.. doxygenfunction:: PDM_mesh_location_points_in_elt_get

.. doxygenfunction:: PDM_mesh_location_point_location_get

.. doxygenfunction:: PDM_mesh_location_cell_vertex_get

.. doxygenfunction:: PDM_mesh_location_part_to_part_get

Finalization
""""""""""""

.. doxygenfunction:: PDM_mesh_location_free



Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  Initialization
  """"""""""""""

  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_create_

  Source mesh definition
  """"""""""""""""""""""

  .. f:subroutine:: pdm_mesh_location_mesh_n_part_set(mloc, n_part)

    Set the number of partitions of the source mesh

    :param mesh_loc: Pointer to PDM_mesh_location object
    :param n_part:   Number of partitions


  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_part_set_

  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_nodal_part_set

  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_part_set_2d_

  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_nodal_part_set_2d

  Target point clouds definition
  """"""""""""""""""""""""""""""

  .. .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_n_part_cloud_set

  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_cloud_set_

  Location computation
  """"""""""""""""""""

  .. .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_method_set

  .. .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_compute

  Results
  """""""
  .. .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_n_located_get

  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_located_get_

  .. .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_n_unlocated_get

  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_unlocated_get_

  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_points_in_elt_get_

  .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_point_location_get_

  .. .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_cell_vertex_get_

  .. .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_part_to_part_get_

  Finalization
  """"""""""""

  .. .. f:autosubroutine:: pdm_mesh_location/pdm_mesh_location_free

.. ifconfig:: enable_fortran_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)


Python API
----------

.. ifconfig:: enable_python_doc == 'ON'

  Initialization
  """"""""""""""

  .. autoclass:: Pypdm.Pypdm.MeshLocation

  Source mesh definition
  """"""""""""""""""""""

  .. autofunction:: Pypdm.Pypdm.MeshLocation.mesh_n_part_set

  .. autofunction:: Pypdm.Pypdm.MeshLocation.part_set

  .. autofunction:: Pypdm.Pypdm.MeshLocation.nodal_part_set

  .. autofunction:: Pypdm.Pypdm.MeshLocation.part_set_2d

  .. autofunction:: Pypdm.Pypdm.MeshLocation.nodal_part_set_2d

  Target point clouds definition
  """"""""""""""""""""""""""""""

  .. autofunction:: Pypdm.Pypdm.MeshLocation.n_part_cloud_set

  .. autofunction:: Pypdm.Pypdm.MeshLocation.cloud_set

  Location computation
  """"""""""""""""""""

  .. autofunction:: Pypdm.Pypdm.MeshLocation.tolerance_set

  .. autofunction:: Pypdm.Pypdm.MeshLocation.method_set

  .. autofunction:: Pypdm.Pypdm.MeshLocation.compute

  Results
  """""""

  .. autofunction:: Pypdm.Pypdm.MeshLocation.located_get

  .. autofunction:: Pypdm.Pypdm.MeshLocation.unlocated_get

  .. autofunction:: Pypdm.Pypdm.MeshLocation.points_in_elt_get

  .. autofunction:: Pypdm.Pypdm.MeshLocation.point_location_get

  .. autofunction:: Pypdm.Pypdm.MeshLocation.cell_vertex_get

  .. autofunction:: Pypdm.Pypdm.MeshLocation.part_to_part_get

.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)


Other
=====
* nearest neighbors search
* distance from surface
* inside cloud surf
* mesh intersection
* extraction of mesh partitions
* iso-surfaces & slices
* overlay
* mesh adaptation/remeshing
