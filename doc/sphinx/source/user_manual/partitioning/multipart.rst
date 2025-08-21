.. _multipart:

Multipart
=========

Description
"""""""""""

**Multipart** is a service for partitioning (possibly multi-domain) block-distributed meshes.
Meshes of dimension 0 (point clouds), 1, 2 and 3 are supported.

The partitioning can be performed using different graph-splitting or geometric methods:

.. doxygenenum:: PDM_split_dual_t



API
"""

.. dropdown:: Initialization


  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_multipart_create



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_multipart_create

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.MultiPart.__init__
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)



.. dropdown:: Set inputs

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_multipart_dmesh_nodal_set

      .. doxygenfunction:: PDM_multipart_dmesh_set

      .. doxygenfunction:: PDM_multipart_domain_interface_shared_set

    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_multipart_dmesh_nodal_set

        .. f:autosubroutine:: PDM_multipart_dmesh_set

        .. f:autosubroutine:: PDM_multipart_domain_interface_shared_set

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.MultiPart.dmesh_nodal_set
          :noindex:
        .. autofunction:: Pypdm.Pypdm.MultiPart.dmesh_set
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Renumbering options

  .. todo::
    List available renumbering methods

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_multipart_set_reordering_options
      .. doxygenfunction:: PDM_multipart_set_reordering_options_vtx



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_multipart_set_reordering_options
        .. f:autosubroutine:: PDM_multipart_set_reordering_options_vtx

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.MultiPart.reordering_set
          :noindex:
        .. autofunction:: Pypdm.Pypdm.MultiPart.reordering_vtx_set
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Perform partitioning


  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_multipart_compute
      .. doxygenfunction:: PDM_multipart_stat_get



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_multipart_compute

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.MultiPart.compute
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Get outputs


  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_multipart_get_part_mesh_nodal

      .. doxygenfunction:: PDM_multipart_part_n_entity_get
      .. doxygenfunction:: PDM_multipart_part_connectivity_get
      .. doxygenfunction:: PDM_multipart_part_ln_to_gn_get
      .. doxygenfunction:: PDM_multipart_part_vtx_coord_get
      .. doxygenfunction:: PDM_multipart_group_get
      .. doxygenfunction:: PDM_multipart_part_graph_comm_get
      .. doxygenfunction:: PDM_multipart_part_ghost_infomation_get
      .. doxygenfunction:: PDM_multipart_partition_color_get
      .. doxygenfunction:: PDM_multipart_part_hyperplane_color_get
      .. doxygenfunction:: PDM_multipart_part_thread_color_get



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_multipart_get_part_mesh_nodal

        .. f:autosubroutine:: PDM_multipart_part_connectivity_get
        .. f:autosubroutine:: PDM_multipart_part_ln_to_gn_get
        .. f:autosubroutine:: PDM_multipart_part_vtx_coord_get
        .. f:autosubroutine:: PDM_multipart_group_get
        .. f:autosubroutine:: PDM_multipart_part_graph_comm_get
        .. f:autosubroutine:: PDM_multipart_part_ghost_infomation_get
        .. f:autosubroutine:: PDM_multipart_partition_color_get

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.MultiPart.part_mesh_nodal_get
          :noindex:

        .. autofunction:: Pypdm.Pypdm.MultiPart.n_entity_get
          :noindex:
        .. autofunction:: Pypdm.Pypdm.MultiPart.connectivity_get
          :noindex:
        .. autofunction:: Pypdm.Pypdm.MultiPart.ln_to_gn_get
          :noindex:
        .. autofunction:: Pypdm.Pypdm.MultiPart.vtx_coord_get
          :noindex:
        .. autofunction:: Pypdm.Pypdm.MultiPart.graph_comm_get
          :noindex:
        .. autofunction:: Pypdm.Pypdm.MultiPart.ghost_information_get
          :noindex:
        .. autofunction:: Pypdm.Pypdm.MultiPart.color_get
          :noindex:
        .. autofunction:: Pypdm.Pypdm.MultiPart.hyper_plane_color_get
          :noindex:
        .. autofunction:: Pypdm.Pypdm.MultiPart.thread_color_get
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Finalization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_multipart_free



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_multipart_free

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      |python_gc|
