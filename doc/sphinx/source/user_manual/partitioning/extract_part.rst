.. _extract_part:

Extract part
============

Description
"""""""""""

**Extract part** is a service for extraction of mesh regions based on arbitrary criteria.
Meshes of dimension 0 (point clouds), 1, 2 and 3 are supported.
Group information associated with extracted parts is preserved.
Facilities for data transfer between input and extracted parts are provided as well.

Three extraction modes are available :

.. doxygenenum:: PDM_extract_part_kind_t

API
"""

.. dropdown:: Initialization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_extract_part_create



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_create

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. py:class:: ExtractPart

          .. automethod:: Pypdm.Pypdm.ExtractPart.__init__

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)



.. dropdown:: Define input mesh

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      One can define the input mesh either as a :ref:`Part Mesh Nodal <pmn>` instance:

      .. doxygenfunction:: PDM_extract_part_part_nodal_set

      Or by providing primitive arrays:

      .. doxygenfunction:: PDM_extract_part_part_set
      .. doxygenfunction:: PDM_extract_part_n_group_set
      .. doxygenfunction:: PDM_extract_part_part_group_set



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        One can define the input mesh either as a :ref:`Part Mesh Nodal <pmn>` instance:

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_part_nodal_set

        Or by providing primitive arrays:

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_part_set
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_n_group_set
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_part_group_set

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. automethod:: Pypdm.Pypdm.ExtractPart.part_set
        .. automethod:: Pypdm.Pypdm.ExtractPart.n_group_set
        .. automethod:: Pypdm.Pypdm.ExtractPart.group_set


      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Define extraction

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_extract_part_selected_lnum_set
      .. doxygenfunction:: PDM_extract_part_target_set
      .. doxygenfunction:: PDM_extract_part_renum_method_set



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_selected_lnum_set
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_target_set

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. automethod:: Pypdm.Pypdm.ExtractPart.selected_lnum_set
        .. automethod:: Pypdm.Pypdm.ExtractPart.target_set


      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Compute extraction

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_extract_part_compute



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_compute

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. automethod:: Pypdm.Pypdm.ExtractPart.compute


      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Get extraction

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      If the input mesh was defined as a :ref:`Part Mesh Nodal <pmn>` structure,
      the extracted mesh is also retrieved in the form of a :ref:`Part Mesh Nodal <pmn>`:

      .. doxygenfunction:: PDM_extract_part_part_mesh_nodal_get

      Otherwise, the extracted mesh can be retrieved using the following accessors:

      .. doxygenfunction:: PDM_extract_part_n_entity_get
      .. doxygenfunction:: PDM_extract_part_connectivity_get
      .. doxygenfunction:: PDM_extract_part_vtx_coord_get
      .. doxygenfunction:: PDM_extract_part_ln_to_gn_get
      .. doxygenfunction:: PDM_extract_part_parent_ln_to_gn_get
      .. doxygenfunction:: PDM_extract_part_group_get
      .. doxygenfunction:: PDM_extract_part_color_get
      .. doxygenfunction:: PDM_extract_part_init_location_get



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        If the input mesh was defined as a :ref:`Part Mesh Nodal <pmn>` instance,
        the extracted mesh is also retrieved in the form of a :ref:`Part Mesh Nodal <pmn>`:

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_part_mesh_nodal_get

        Otherwise, the extracted mesh can be retrieved using the following accessors:

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_n_entity_get
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_connectivity_get
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_vtx_coord_get
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_ln_to_gn_get
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_parent_ln_to_gn_get
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_group_get

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. automethod:: Pypdm.Pypdm.ExtractPart.n_entity_get
        .. automethod:: Pypdm.Pypdm.ExtractPart.connectivity_get
        .. automethod:: Pypdm.Pypdm.ExtractPart.vtx_coord_get
        .. automethod:: Pypdm.Pypdm.ExtractPart.ln_to_gn_get
        .. automethod:: Pypdm.Pypdm.ExtractPart.parent_ln_to_gn_get
        .. automethod:: Pypdm.Pypdm.ExtractPart.group_get


      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Data transfer

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      Application-specific data can be can be transferred between the extracted and input meshes.

      If the extraction was performed in PDM_EXTRACT_PART_KIND_LOCAL mode, the transfer is performed locally using the *parent* indirection.

      .. doxygenfunction:: PDM_extract_part_parent_lnum_get

      Otherwise, the transfer is performed by means of a :ref:`Part-to-part <ptp>` instance for each type of entity.

      .. doxygenfunction:: PDM_extract_part_part_to_part_get
      .. doxygenfunction:: PDM_extract_part_part_to_part_group_get



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        Application-specific data can be can be transferred between the extracted and input meshes.

        If the extraction was performed in PDM_EXTRACT_PART_KIND_LOCAL mode, the transfer is performed locally using the *parent* indirection.

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_parent_lnum_get

        Otherwise, the transfer is performed by means of a :ref:`Part-to-part <ptp>` instance for each type of entity.

        .. note::

          *Direct* exchanges go from extraction to input and *reverse* exchanges go from input to extraction.

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_part_to_part_get
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_part_to_part_group_get

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. automethod:: Pypdm.Pypdm.ExtractPart.part_to_part_get
        .. automethod:: Pypdm.Pypdm.ExtractPart.part_to_part_group_get


      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Finalization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_extract_part_partial_free
      .. doxygenfunction:: PDM_extract_part_free



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_partial_free
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_free

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      |python_gc|


Examples
""""""""

.. dropdown:: Nodal mesh

  .. tab-set::
    :sync-group: language

    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        The following example shows how to perform an extraction (local or with redistribution) with a nodal mesh (extract from ``test/pdm_t_extract_part_nodal_f.F90``).

        Mesh generation
        ~~~~~~~~~~~~~~~

        We start by generating a partitioned 3d mesh.

        .. figure:: ../../../../../doc/images/extract_part/01_initial_mesh.png
          :width: 200
          :alt: Initial partitioned mesh

          Initial partitioned mesh (colored by MPI rank)

        .. literalinclude:: ../../../../../test/pdm_t_extract_part_nodal_f.F90
          :language: fortran
          :dedent: 2
          :lines: 24-165


        Extraction
        ~~~~~~~~~~

        We then extract the cells that cross the :math:`(x = 0)` plane.

        .. figure:: ../../../../../doc/images/extract_part/02_extraction_local.png
          :width: 200
          :alt: Extracted mesh

          Extracted mesh (colored by MPI rank)

        .. literalinclude:: ../../../../../test/pdm_t_extract_part_nodal_f.F90
          :language: fortran
          :dedent: 2
          :lines: 168-250


        Data transfer
        ~~~~~~~~~~~~~

        Finally, we transfer fields from the initial mesh to the extraction.

        Let's create a dummy field ...

        .. figure:: ../../../../../doc/images/extract_part/03_initial_field.png
          :width: 200
          :alt: Initial field

          Field on the initial mesh

        .. literalinclude:: ../../../../../test/pdm_t_extract_part_nodal_f.F90
          :language: fortran
          :dedent: 2
          :lines: 255-257

        ... and transfer it to the extraction  (see ``test/pdm_t_extract_part_nodal_f.F90`` for the detailed implementation of the ``data_transfer`` procedure).

        .. figure:: ../../../../../doc/images/extract_part/04_extracted_field.png
          :width: 200
          :alt: Extracted field

          Field transferred to the extracted mesh

        .. literalinclude:: ../../../../../test/pdm_t_extract_part_nodal_f.F90
          :language: fortran
          :dedent: 2
          :lines: 266-285




      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)

