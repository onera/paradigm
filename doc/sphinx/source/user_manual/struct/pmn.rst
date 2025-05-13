.. _pmn:

Partitioned Nodal Mesh Structure
================================

C API
-----

Initialization
""""""""""""""

.. doxygenfunction:: PDM_part_mesh_nodal_create

Setters
"""""""

.. doxygenfunction:: PDM_part_mesh_nodal_coord_set

.. doxygenfunction:: PDM_part_mesh_nodal_vtx_gnum_set

.. doxygenfunction:: PDM_part_mesh_nodal_coord_from_parent_set

.. doxygenfunction:: PDM_part_mesh_nodal_section_add

.. doxygenfunction:: PDM_part_mesh_nodal_section_std_set

.. doxygenfunction:: PDM_part_mesh_nodal_section_std_ho_set

.. doxygenfunction:: PDM_part_mesh_nodal_add_part_mesh_nodal_elmts

.. doxygenfunction:: PDM_part_mesh_nodal_section_poly2d_set

.. doxygenfunction:: PDM_part_mesh_nodal_section_poly3d_set

.. doxygenfunction:: PDM_part_mesh_nodal_cell3d_cellface_add

.. doxygenfunction:: PDM_part_mesh_nodal_face2d_faceedge_add

.. doxygenfunction:: PDM_part_mesh_nodal_cells_cellvtx_add

.. doxygenfunction:: PDM_part_mesh_nodal_faces_facevtx_add

.. doxygenfunction:: PDM_part_mesh_nodal_part_comm_graph_set

Compute
"""""""

.. doxygenfunction:: PDM_part_mesh_nodal_section_elt_extents_compute

.. doxygenfunction:: PDM_part_mesh_nodal_section_elt_center_compute

.. doxygenfunction:: PDM_part_mesh_nodal_g_num_in_section_compute

.. doxygenfunction:: PDM_part_mesh_nodal_compute_sections_idx

.. doxygenfunction:: PDM_part_mesh_nodal_compute_straddling_entities

Getters
"""""""

.. doxygenfunction:: PDM_part_mesh_nodal_n_part_get

.. doxygenfunction:: PDM_part_mesh_nodal_mesh_dimension_get

.. doxygenfunction:: PDM_part_mesh_nodal_n_vtx_get

.. doxygenfunction:: PDM_part_mesh_nodal_vtx_coord_get

.. doxygenfunction:: PDM_part_mesh_nodal_vtx_g_num_get

.. doxygenfunction:: PDM_part_mesh_nodal_n_section_in_geom_kind_get

.. doxygenfunction:: PDM_part_mesh_nodal_sections_id_in_geom_kind_get

.. doxygenfunction:: PDM_part_mesh_nodal_section_elt_type_get

.. doxygenfunction:: PDM_part_mesh_nodal_section_in_geom_kind_elt_type_get

.. doxygenfunction:: PDM_part_mesh_nodal_section_n_elt_get

.. doxygenfunction:: PDM_part_mesh_nodal_section_std_get

.. doxygenfunction:: PDM_part_mesh_nodal_section_std_ho_get

.. doxygenfunction:: PDM_part_mesh_nodal_section_parent_num_get

.. doxygenfunction:: PDM_part_mesh_nodal_g_num_get

.. doxygenfunction:: PDM_part_mesh_nodal_section_elt_center_get

.. doxygenfunction:: PDM_part_mesh_nodal_section_poly2d_get

.. doxygenfunction:: PDM_part_mesh_nodal_section_poly3d_get

.. doxygenfunction:: PDM_part_mesh_nodal_section_poly3d_cell_vtx_connect_get

.. doxygenfunction:: PDM_part_mesh_nodal_n_elmts_get

.. doxygenfunction:: PDM_part_mesh_nodal_g_num_get_from_part

.. doxygenfunction:: PDM_part_mesh_nodal_is_set_coord_from_parent

.. doxygenfunction:: PDM_part_mesh_nodal_section_g_num_get

.. doxygenfunction:: PDM_part_mesh_nodal_num_elmt_parent_to_local_get

.. doxygenfunction:: PDM_part_mesh_nodal_vertices_parent_get

.. doxygenfunction:: PDM_part_mesh_nodal_vertices_g_num_parent_get

.. doxygenfunction:: PDM_part_mesh_nodal_section_id_and_geom_kind_get

.. doxygenfunction:: PDM_part_mesh_nodal_section_id_from_geom_kind_get

.. doxygenfunction:: PDM_part_mesh_nodal_n_section_get

.. doxygenfunction:: PDM_part_mesh_nodal_sections_id_get

.. doxygenfunction:: PDM_part_mesh_nodal_group_get

.. doxygenfunction:: PDM_part_mesh_nodal_n_group_get

.. doxygenfunction:: PDM_part_mesh_nodal_part_mesh_nodal_elmts_get

.. doxygenfunction:: PDM_part_mesh_nodal_principal_geom_kind_get

.. doxygenfunction:: PDM_part_mesh_nodal_cell_vtx_connect_get

.. doxygenfunction:: PDM_part_mesh_nodal_part_comm_graph_get

Free
""""

.. doxygenfunction:: PDM_part_mesh_nodal_free

.. doxygenfunction:: PDM_part_mesh_nodal_partial_free

.. doxygenfunction:: PDM_part_mesh_nodal_reset

.. doxygenfunction:: PDM_part_mesh_nodal_section_elt_center_reset

Debug
"""""

.. doxygenfunction:: PDM_part_mesh_nodal_dump_vtk

Example
"""""""

It is necessary to take care of a subtlety relative to the internal structure
of this mesh container. Some functions return local section indices, while others
expect global indices as input. Here's how to proceed :

.. code-block:: c

  int *i_surface_sections = PDM_part_mesh_nodal_sections_id_in_geom_kind_get(pmn, PDM_GEOMETRY_KIND_SURFACIC);

  int i_global_surface_0_section = PDM_part_mesh_nodal_section_id_from_geom_kind_get(pmn,
                                                                                     PDM_GEOMETRY_KIND_SURFACIC,
                                                                                     i_surface_sections[0]);

Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  .. f:autosubroutine:: PDM_part_mesh_nodal_create
  .. f:autosubroutine:: PDM_part_mesh_nodal_coord_set
  .. f:autosubroutine:: PDM_part_mesh_nodal_vtx_gnum_set
  .. f:autofunction::   PDM_part_mesh_nodal_section_add
  .. f:autosubroutine:: PDM_part_mesh_nodal_section_std_set
  .. f:autosubroutine:: PDM_part_mesh_nodal_cells_cellvtx_add
  .. f:autosubroutine:: PDM_part_mesh_nodal_faces_facevtx_add
  .. f:autosubroutine:: PDM_part_mesh_nodal_cell3d_cellface_add
  .. f:autosubroutine:: PDM_part_mesh_nodal_face2d_faceedge_add
  .. f:autosubroutine:: PDM_part_mesh_nodal_principal_geom_kind_get
  .. f:autofunction::   PDM_part_mesh_nodal_n_section_in_geom_kind_get
  .. f:autosubroutine:: PDM_part_mesh_nodal_sections_id_in_geom_kind_get
  .. f:autosubroutine:: PDM_part_mesh_nodal_section_in_geom_kind_elt_type_get
  .. f:autosubroutine:: PDM_part_mesh_nodal_section_n_elt_get
  .. f:autosubroutine:: PDM_part_mesh_nodal_section_elt_type_get
  .. f:autosubroutine:: PDM_part_mesh_nodal_section_std_get
  .. f:autosubroutine:: PDM_part_mesh_nodal_cell_vtx_connect_get
  .. f:autosubroutine:: PDM_part_mesh_nodal_n_vtx_get
  .. f:autosubroutine:: PDM_part_mesh_nodal_vtx_coord_get
  .. f:autosubroutine:: PDM_part_mesh_nodal_vtx_g_num_get
  .. f:autosubroutine:: PDM_part_mesh_nodal_free

.. ifconfig:: enable_fortran_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)

Python API
----------

.. ifconfig:: enable_python_doc == 'ON'


  .. py:class:: PartMeshNodal

    Python class to store meshes defined by nodal connectivity. Once initialized, all the following
    methods apply to a :class:`PartMeshNodal` instance.

    .. rubric:: Initialization

    .. autofunction:: Pypdm.Pypdm.PartMeshNodal.__cinit__

    .. rubric:: Methods summary

    .. autosummary::
      :nosignatures:

      ~Pypdm.Pypdm.PartMeshNodal.set_coordinates
      ~Pypdm.Pypdm.PartMeshNodal.add_section
      ~Pypdm.Pypdm.PartMeshNodal.set_section
      ~Pypdm.Pypdm.PartMeshNodal.n_group_set
      ~Pypdm.Pypdm.PartMeshNodal.group_set
      ~Pypdm.Pypdm.PartMeshNodal.get_sections

    .. automethod:: Pypdm.Pypdm.PartMeshNodal.set_coordinates
    .. automethod:: Pypdm.Pypdm.PartMeshNodal.add_section
    .. automethod:: Pypdm.Pypdm.PartMeshNodal.set_section
    .. automethod:: Pypdm.Pypdm.PartMeshNodal.n_group_set
    .. automethod:: Pypdm.Pypdm.PartMeshNodal.group_set
    .. automethod:: Pypdm.Pypdm.PartMeshNodal.get_sections


.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)


