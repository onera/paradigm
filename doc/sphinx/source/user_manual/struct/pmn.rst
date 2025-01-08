.. _pmn:

Partitionned Nodal Mesh Structure
=================================

C API
-----

Initialization
""""""""""""""

.. doxygenfunction:: PDM_part_mesh_nodal_create

Setters
"""""""

.. doxygenfunction:: PDM_part_mesh_nodal_coord_set

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

Compute
"""""""

.. doxygenfunction:: PDM_part_mesh_nodal_section_elt_extents_compute

.. doxygenfunction:: PDM_part_mesh_nodal_section_elt_center_compute

.. doxygenfunction:: PDM_part_mesh_nodal_g_num_in_section_compute

.. doxygenfunction:: PDM_part_mesh_nodal_compute_sections_idx

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

  .. f:subroutine:: PDM_part_mesh_nodal_section_n_elt_get

    Get number of section elements

    :param c_ptr   pmn       [in]:  C pointer to \ref PDM_part_mesh_nodal_t instance
    :param integer i_section [in]:  Section identifier
    :param integer i_part    [in]:  Partition identifier
    :param integer n_elt     [out]: Number of elements

  .. f:subroutine:: PDM_part_mesh_nodal_section_std_get

    Return standard section description

    :param  c_ptr          pmn                 [in]:  C pointer to \ref PDM_part_mesh_nodal_t instance
    :param  integer        i_section           [in]:  Section identifier
    :param  integer        i_part              [in]:  Partition identifier
    :param  pdm_l_num_s(:) connec              [out]: Connectivity
    :param  pdm_g_num_s(:) numabs              [out]: Global ids
    :param  pdm_l_num_s(:) parent_num          [out]: Parent local ids or *null()*
    :param  pdm_g_num_s(:) parent_entity_g_num [out]: Parent global ids or *null()*
    :param  integer        ownership           [in]:  Data ownership

  .. f:subroutine:: PDM_part_mesh_nodal_vtx_g_num_get

    Return global ids of vertices

    :param  c_ptr          pmn           [in]:  C pointer to \ref PDM_part_mesh_nodal_t instance
    :param  integer        i_part        [in]:  Partition identifier
    :param  pdm_g_num_s(:) vtx_ln_to_gn  [out]: Global ids of vertices (size = ``n_vtx``)

  .. f:subroutine:: PDM_part_mesh_nodal_n_vtx_get

    Return number of vertices

    :param  c_ptr   pmn       [in]:  C pointer to \ref PDM_part_mesh_nodal_t instance
    :param  integer i_part    [in]:  Partition identifier
    :param  integer n_vtx     [out]: Number of vertices

  .. f:subroutine:: PDM_part_mesh_nodal_section_elt_type_get

    Return type of section

    :param  c_ptr   pmn       [in]:  C pointer to \ref PDM_part_mesh_nodal_t instance
    :param  integer i_section [in]:  Section identifier
    :param  integer elt_t     [out]: Type of section

.. ifconfig:: enable_fortran_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)

Python API
----------

.. ifconfig:: enable_python_doc == 'ON'


  .. py:class:: PartMeshNodal

    Interface to build face from Element->Vtx connectivity.
    Once initialized, all the following
    methods apply to a :class:`PartMeshNodal` instance.

    .. rubric:: Initialization

    .. autoclass:: Pypdm.Pypdm.PartMeshNodal.__init__

    .. rubric:: Methods summary

    .. autosummary::
      :nosignatures:

      ~Pypdm.Pypdm.PartMeshNodal.dim_get
      ~Pypdm.Pypdm.PartMeshNodal.get_sections
      ~Pypdm.Pypdm.PartMeshNodal.vtx_g_num_get

    .. automethod:: Pypdm.Pypdm.PartMeshNodal.dim_get
    .. automethod:: Pypdm.Pypdm.PartMeshNodal.get_sections
    .. automethod:: Pypdm.Pypdm.PartMeshNodal.vtx_g_num_get


.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)


