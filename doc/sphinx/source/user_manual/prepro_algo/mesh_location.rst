.. _mesh_location:

Mesh location
=============


C API
-----

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

    :param c_ptr   mesh_loc [in]: C pointer to PDM_mesh_location_t object
    :param integer n_part   [in]:   Number of partitions


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

  .. py:class:: MeshLocation

    Python structure to perform mesh location operations. Once initialized, all the following
    methods apply to a :class:`MeshLocation` instance.

  .. rubric:: Initialization

  .. autofunction:: Pypdm.Pypdm.MeshLocation.__init__

  .. rubric:: Instance attributes

  .. autoattribute:: Pypdm.Pypdm.MeshLocation.tolerance
  .. autoattribute:: Pypdm.Pypdm.MeshLocation.method

  .. rubric:: Methods summary
  
  .. autosummary::
    ~Pypdm.Pypdm.MeshLocation.mesh_n_part_set
    ~Pypdm.Pypdm.MeshLocation.part_set
    ~Pypdm.Pypdm.MeshLocation.nodal_part_set
    ~Pypdm.Pypdm.MeshLocation.part_set_2d
    ~Pypdm.Pypdm.MeshLocation.nodal_part_set_2d
    ~Pypdm.Pypdm.MeshLocation.n_part_cloud_set
    ~Pypdm.Pypdm.MeshLocation.cloud_set
    ~Pypdm.Pypdm.MeshLocation.compute
    ~Pypdm.Pypdm.MeshLocation.located_get
    ~Pypdm.Pypdm.MeshLocation.unlocated_get
    ~Pypdm.Pypdm.MeshLocation.location_get
    ~Pypdm.Pypdm.MeshLocation.points_in_elt_get
    ~Pypdm.Pypdm.MeshLocation.point_location_get
    ~Pypdm.Pypdm.MeshLocation.cell_vertex_get
    ~Pypdm.Pypdm.MeshLocation.part_to_part_get


  .. rubric:: Source mesh definition

  .. autofunction:: Pypdm.Pypdm.MeshLocation.mesh_n_part_set

  .. autofunction:: Pypdm.Pypdm.MeshLocation.part_set

  .. autofunction:: Pypdm.Pypdm.MeshLocation.nodal_part_set

  .. autofunction:: Pypdm.Pypdm.MeshLocation.part_set_2d

  .. autofunction:: Pypdm.Pypdm.MeshLocation.nodal_part_set_2d

  .. rubric:: Target point clouds definition

  .. autofunction:: Pypdm.Pypdm.MeshLocation.n_part_cloud_set

  .. autofunction:: Pypdm.Pypdm.MeshLocation.cloud_set

  .. rubric:: Location computation

  .. autofunction:: Pypdm.Pypdm.MeshLocation.compute

  .. rubric:: Results

  .. autofunction:: Pypdm.Pypdm.MeshLocation.located_get

  .. autofunction:: Pypdm.Pypdm.MeshLocation.unlocated_get

  .. autofunction:: Pypdm.Pypdm.MeshLocation.location_get

  .. autofunction:: Pypdm.Pypdm.MeshLocation.points_in_elt_get

  .. autofunction:: Pypdm.Pypdm.MeshLocation.point_location_get

  .. autofunction:: Pypdm.Pypdm.MeshLocation.cell_vertex_get

  .. autofunction:: Pypdm.Pypdm.MeshLocation.part_to_part_get

.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)
