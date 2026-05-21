.. _part_mesh_nodal:

Partitioned Nodal Mesh
======================

Description
"""""""""""

**Part Mesh Nodal** is a partitioned mesh data structure in which elements of different types are stored in separate containers called *sections*.
These sections can be addressed either globally or by spatial dimension (*geometry kind*).

.. If a mesh is composed of elements of multiple dimensions (e.g. a volume mesh with its surface boundary), each dimension is stored in its own container as well (one **Part Mesh Nodal Elmts** instance per dimension).

If a **Part Mesh Nodal** is constructed from a soup of mixed-type elements, the indirection from the sections to the original soup is stored in the data structure (*parent_num*).

The structure can also store group information and the inter-partition communication graph (:ref:`Part Comm Graph <part_comm_graph>`) for each geometry kind as well as for the vertices.


API
"""

.. dropdown:: Initialization


  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_mesh_nodal_create



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_mesh_nodal_create

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. py:class:: PartMeshNodal

          .. automethod:: Pypdm.Pypdm.PartMeshNodal.__init__

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Set vertices


  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_mesh_nodal_coord_set
      .. doxygenfunction:: PDM_part_mesh_nodal_vtx_gnum_set
      .. .. doxygenfunction:: PDM_part_mesh_nodal_coord_from_parent_set


    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_mesh_nodal_coord_set
        .. f:autosubroutine:: PDM_part_mesh_nodal_vtx_gnum_set

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. automethod:: Pypdm.Pypdm.PartMeshNodal.set_coordinates

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Set elements

  The mesh elements can be defined in multiple ways.

  If only a soup of mixed elements is available (of same geometry kind), **Part Mesh Nodal** can divide them into appropriate sections automatically:

  .. dropdown:: Mixed setters

    .. tab-set::
      :sync-group: language

      .. tab-item:: C
        :sync: C

        .. doxygenfunction:: PDM_part_mesh_nodal_cell3d_cellface_add
        .. doxygenfunction:: PDM_part_mesh_nodal_face2d_faceedge_add
        .. doxygenfunction:: PDM_part_mesh_nodal_cells_cellvtx_add
        .. doxygenfunction:: PDM_part_mesh_nodal_faces_facevtx_add



      .. tab-item:: Fortran
        :sync: Fortran

        .. ifconfig:: enable_fortran_doc == 'ON'

          .. f:autosubroutine:: PDM_part_mesh_nodal_cell3d_cellface_add
          .. f:autosubroutine:: PDM_part_mesh_nodal_face2d_faceedge_add
          .. f:autosubroutine:: PDM_part_mesh_nodal_cells_cellvtx_add
          .. f:autosubroutine:: PDM_part_mesh_nodal_faces_facevtx_add

        .. ifconfig:: enable_fortran_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)

      .. tab-item:: Python
        :sync: Python

        .. ifconfig:: enable_python_doc == 'ON'

          .. automethod:: Pypdm.Pypdm.PartMeshNodal.set_sections_from_cell_face
          .. automethod:: Pypdm.Pypdm.PartMeshNodal.set_sections_from_face_edge

        .. ifconfig:: enable_python_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)



  |

  If the mesh elements are already separated by types, the sections can be specified one by one:

  .. dropdown:: Section setters

    .. tab-set::
      :sync-group: language

      .. tab-item:: C
        :sync: C

        .. doxygenfunction:: PDM_part_mesh_nodal_section_add

        .. doxygenfunction:: PDM_part_mesh_nodal_section_std_set
        .. doxygenfunction:: PDM_part_mesh_nodal_section_std_ho_set
        .. doxygenfunction:: PDM_part_mesh_nodal_section_poly2d_set
        .. doxygenfunction:: PDM_part_mesh_nodal_section_poly3d_set



      .. tab-item:: Fortran
        :sync: Fortran

        .. ifconfig:: enable_fortran_doc == 'ON'

          .. f:autofunction::   PDM_part_mesh_nodal_section_add

          .. f:autosubroutine:: PDM_part_mesh_nodal_section_std_set

        .. ifconfig:: enable_fortran_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



      .. tab-item:: Python
        :sync: Python

        .. ifconfig:: enable_python_doc == 'ON'

          .. automethod:: Pypdm.Pypdm.PartMeshNodal.add_section
          .. automethod:: Pypdm.Pypdm.PartMeshNodal.set_std_section
          .. automethod:: Pypdm.Pypdm.PartMeshNodal.set_poly2d_section
          .. automethod:: Pypdm.Pypdm.PartMeshNodal.set_poly3d_section

        .. ifconfig:: enable_python_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)


  |

  Finally, if a **Part Mesh Nodal Elmts** is already available, it can simply be added to the structure:

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_mesh_nodal_add_part_mesh_nodal_elmts




.. dropdown:: Set groups


  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_mesh_nodal_n_group_set
      .. doxygenfunction:: PDM_part_mesh_nodal_group_set



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_mesh_nodal_n_group_set
        .. f:autosubroutine:: PDM_part_mesh_nodal_group_set

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. automethod:: Pypdm.Pypdm.PartMeshNodal.n_group_set
        .. automethod:: Pypdm.Pypdm.PartMeshNodal.group_set

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Set inter-partition communication graph


  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_mesh_nodal_part_comm_graph_set
      .. doxygenfunction:: PDM_part_mesh_nodal_part_comm_graph_vtx_set


    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_mesh_nodal_part_comm_graph_set
        .. f:autosubroutine:: PDM_part_mesh_nodal_part_comm_graph_vtx_set


      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)




.. dropdown:: Computations


  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_mesh_nodal_section_elt_extents_compute
      .. doxygenfunction:: PDM_part_mesh_nodal_section_elt_center_compute
      .. doxygenfunction:: PDM_part_mesh_nodal_g_num_in_section_compute
      .. doxygenfunction:: PDM_part_mesh_nodal_compute_sections_idx
      .. doxygenfunction:: PDM_part_mesh_nodal_compute_straddling_entities




.. dropdown:: Get vertices


  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_mesh_nodal_n_vtx_get
      .. doxygenfunction:: PDM_part_mesh_nodal_vtx_coord_get
      .. doxygenfunction:: PDM_part_mesh_nodal_vtx_g_num_get


    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_mesh_nodal_n_vtx_get
        .. f:autosubroutine:: PDM_part_mesh_nodal_vtx_coord_get
        .. f:autosubroutine:: PDM_part_mesh_nodal_vtx_g_num_get

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. automethod:: Pypdm.Pypdm.PartMeshNodal.coord_get
        .. automethod:: Pypdm.Pypdm.PartMeshNodal.vtx_g_num_get

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Get elements

  The mesh elements can be accessed either all together or by section.
  The *parent_num* indirection is used to link the two representations.

  .. dropdown:: Mixed getters

    .. tab-set::
      :sync-group: language

      .. tab-item:: C
        :sync: C

        .. doxygenfunction:: PDM_part_mesh_nodal_principal_geom_kind_get
        .. doxygenfunction:: PDM_part_mesh_nodal_n_elmts_get
        .. doxygenfunction:: PDM_part_mesh_nodal_cell_vtx_connect_get



      .. tab-item:: Fortran
        :sync: Fortran

        .. ifconfig:: enable_fortran_doc == 'ON'

          .. f:autofunction::   PDM_part_mesh_nodal_principal_geom_kind_get
          .. f:autosubroutine:: PDM_part_mesh_nodal_n_elmts_get
          .. f:autosubroutine:: PDM_part_mesh_nodal_cell_vtx_connect_get

        .. ifconfig:: enable_fortran_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



  .. dropdown:: Section getters

    .. tab-set::
      :sync-group: language

      .. tab-item:: C
        :sync: C

        .. doxygenfunction:: PDM_part_mesh_nodal_n_section_get
        .. doxygenfunction:: PDM_part_mesh_nodal_sections_id_get
        .. doxygenfunction:: PDM_part_mesh_nodal_section_elt_type_get
        .. doxygenfunction:: PDM_part_mesh_nodal_g_num_get
        .. doxygenfunction:: PDM_part_mesh_nodal_section_parent_num_get
        .. doxygenfunction:: PDM_part_mesh_nodal_section_std_get
        .. doxygenfunction:: PDM_part_mesh_nodal_section_std_ho_get
        .. doxygenfunction:: PDM_part_mesh_nodal_section_poly2d_get
        .. doxygenfunction:: PDM_part_mesh_nodal_section_poly3d_get
        .. doxygenfunction:: PDM_part_mesh_nodal_section_poly3d_cell_vtx_connect_get
        .. doxygenfunction:: PDM_part_mesh_nodal_section_elt_center_get



      .. tab-item:: Fortran
        :sync: Fortran

        .. ifconfig:: enable_fortran_doc == 'ON'

          .. .. f:autosubroutine:: PDM_part_mesh_nodal_n_section_get
          .. .. f:autosubroutine:: PDM_part_mesh_nodal_sections_id_get
          .. .. f:autosubroutine:: PDM_part_mesh_nodal_g_num_get
          .. .. f:autosubroutine:: PDM_part_mesh_nodal_section_parent_num_get
          .. f:autosubroutine:: PDM_part_mesh_nodal_section_std_get
          .. .. f:autosubroutine:: PDM_part_mesh_nodal_section_std_ho_get
          .. .. f:autosubroutine:: PDM_part_mesh_nodal_section_poly2d_get
          .. .. f:autosubroutine:: PDM_part_mesh_nodal_section_poly3d_get
          .. .. f:autosubroutine:: PDM_part_mesh_nodal_section_poly3d_cell_vtx_connect_get
          .. .. f:autosubroutine:: PDM_part_mesh_nodal_section_elt_center_get


        .. ifconfig:: enable_fortran_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



      .. tab-item:: Python
        :sync: Python

        .. ifconfig:: enable_python_doc == 'ON'

          .. automethod:: Pypdm.Pypdm.PartMeshNodal.get_sections

        .. ifconfig:: enable_python_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Get groups


  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_mesh_nodal_n_group_get
      .. doxygenfunction:: PDM_part_mesh_nodal_group_get



    .. tab-item:: Fortran
      :sync: Fortran

      .. f:autosubroutine:: PDM_part_mesh_nodal_n_group_get
      .. f:autosubroutine:: PDM_part_mesh_nodal_group_get



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. automethod:: Pypdm.Pypdm.PartMeshNodal.get_n_group
        .. automethod:: Pypdm.Pypdm.PartMeshNodal.get_group

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Get inter-partition communication graph


  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_mesh_nodal_part_comm_graph_get
      .. doxygenfunction:: PDM_part_mesh_nodal_part_comm_graph_vtx_get


    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_mesh_nodal_part_comm_graph_get
        .. f:autosubroutine:: PDM_part_mesh_nodal_part_comm_graph_vtx_get


      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)




.. dropdown:: Finalization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_mesh_nodal_partial_free
      .. doxygenfunction:: PDM_part_mesh_nodal_free



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_mesh_nodal_partial_free
        .. f:autosubroutine:: PDM_part_mesh_nodal_free

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      |python_gc|
