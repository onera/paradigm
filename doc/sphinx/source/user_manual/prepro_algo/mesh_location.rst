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

.. doxygenfunction:: PDM_mesh_mesh_global_data_set

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

.. doxygenfunction:: pdm_mesh_location_free




Fortran API
-----------




Python API
----------

.. autoclass:: Pypdm.Pypdm.MeshLocation
  :members:
