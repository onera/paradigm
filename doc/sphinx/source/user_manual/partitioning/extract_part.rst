.. _extract_part:

Extract part
============


.. tabs::

  .. tab:: C API

    .. ifconfig:: fake_bool == 'ON'

      Initialization
      """"""""""""""

        .. doxygenfunction:: PDM_extract_part_create


      Define input mesh
      """""""""""""""""

        One can define the input mesh either as a :ref:`Part Mesh Nodal <pmn>` structure:

        .. doxygenfunction:: PDM_extract_part_part_nodal_set

        Or by providing primitive arrays:

        .. doxygenfunction:: PDM_extract_part_part_set
        .. doxygenfunction:: PDM_extract_part_n_group_set
        .. doxygenfunction:: PDM_extract_part_part_group_set


      Define extraction
      """""""""""""""""

        .. doxygenfunction:: PDM_extract_part_selected_lnum_set
        .. doxygenfunction:: PDM_extract_part_target_set
        .. doxygenfunction:: PDM_extract_part_renum_method_set


      Compute extraction
      """"""""""""""""""

        .. doxygenfunction:: PDM_extract_part_compute


      Get extraction
      """"""""""""""

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


      Data transfer
      """""""""""""

        Application-specific data can be can be transferred between the extracted and input meshes.

        If the extraction was performed in PDM_EXTRACT_PART_KIND_LOCAL mode, the transfer is performed locally using the *parent* indirection.

        .. doxygenfunction:: PDM_extract_part_parent_lnum_get

        Otherwise, the transfer is performed by means of a :ref:`Part-to-part <ptp>` instance for each type of entity.

        .. doxygenfunction:: PDM_extract_part_part_to_part_get
        .. doxygenfunction:: PDM_extract_part_part_to_part_group_get


      Finalization
      """"""""""""

        .. doxygenfunction:: PDM_extract_part_partial_free
        .. doxygenfunction:: PDM_extract_part_free





  .. tab:: Fortran API

    .. ifconfig:: enable_fortran_doc == 'ON'

      Initialization
      """"""""""""""

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_create

      Define input mesh
      """""""""""""""""

        One can define the input mesh either as a :ref:`Part Mesh Nodal <pmn>` structure:

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_part_nodal_set

        Or by providing primitive arrays:

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_part_set
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_n_group_set
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_part_group_set

      Define extraction
      """""""""""""""""

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_selected_lnum_set
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_target_set

      Compute extraction
      """"""""""""""""""

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_compute

      Get extraction
      """"""""""""""

        If the input mesh was defined as a :ref:`Part Mesh Nodal <pmn>` structure,
        the extracted mesh is also retrieved in the form of a :ref:`Part Mesh Nodal <pmn>`:

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_part_mesh_nodal_get

        Otherwise, the extracted mesh can be retrieved using the following accessors:

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_n_entity_get
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_connectivity_get
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_vtx_coord_get
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_ln_to_gn_get
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_parent_ln_to_gn_get
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_group_get

      Data transfer
      """""""""""""

        Application-specific data can be can be transferred between the extracted and input meshes.

        If the extraction was performed in PDM_EXTRACT_PART_KIND_LOCAL mode, the transfer is performed locally using the *parent* indirection.

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_parent_lnum_get

        Otherwise, the transfer is performed by means of a :ref:`Part-to-part <ptp>` instance for each type of entity.

        .. note::

          Note that *direct* exchanges go from extraction to input and *reverse* exchanges go from input to extraction.

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_part_to_part_get

      Finalization
      """"""""""""

        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_partial_free
        .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_free

    .. ifconfig:: enable_fortran_doc == 'OFF'

      .. warning::
        Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



  .. tab:: Python API

    .. ifconfig:: enable_python_doc == 'ON'

      .. py:class:: ExtractPart

        Python structure to perform extraction operations. Once initialized, all the following
        methods apply to a :class:`ExtractPart` instance.


        .. rubric:: Initialization

        .. autofunction:: Pypdm.Pypdm.ExtractPart.__init__


        .. rubric:: Define input mesh

        .. autofunction:: Pypdm.Pypdm.ExtractPart.part_set
        .. autofunction:: Pypdm.Pypdm.ExtractPart.part_n_group_set
        .. autofunction:: Pypdm.Pypdm.ExtractPart.part_group_set


        .. rubric:: Define extraction

        .. autofunction:: Pypdm.Pypdm.ExtractPart.selected_lnum_set
        .. autofunction:: Pypdm.Pypdm.ExtractPart.target_set


        .. rubric:: Compute extraction

        .. autofunction:: Pypdm.Pypdm.ExtractPart.compute


        .. rubric:: Get extraction

        .. autofunction:: Pypdm.Pypdm.ExtractPart.n_entity_get
        .. autofunction:: Pypdm.Pypdm.ExtractPart.connectivity_get
        .. autofunction:: Pypdm.Pypdm.ExtractPart.vtx_coord_get
        .. autofunction:: Pypdm.Pypdm.ExtractPart.ln_to_gn_get
        .. autofunction:: Pypdm.Pypdm.ExtractPart.parent_ln_to_gn_get
        .. autofunction:: Pypdm.Pypdm.ExtractPart.extract_part_group_get


        .. rubric:: Data transfer

        .. autofunction:: Pypdm.Pypdm.ExtractPart.part_to_part_get
        .. autofunction:: Pypdm.Pypdm.ExtractPart.part_to_part_group_get



    .. ifconfig:: enable_python_doc == 'OFF'

      .. warning::
        Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)
