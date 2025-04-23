.. _extract_part:

Extract part
============


C API
-----

.. todo::
  ...





Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  Initialization
  """"""""""""""

    .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_create

  Input mesh definition
  """"""""""""""""""""""

    One can define the input mesh either as a :ref:`Part Mesh Nodal <pmn>` structure:

    .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_part_nodal_set

    Or by by providing primitive arrays:

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

    Otherwise, the extracted mesh can be retrived using the following accessors:

    .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_n_entity_get
    .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_connectivity_get
    .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_vtx_coord_get
    .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_ln_to_gn_get
    .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_parent_ln_to_gn_get
    .. f:autosubroutine:: pdm_extract_part/PDM_extract_part_group_get

  Data transfer
  """""""""""""

    Application-specific data can also be transferred from the extracted mesh to the input mesh, and vice versa.

    If the extraction was performed in ``PDM_EXTRACT_PART_KIND_LOCAL`` mode, the transfer is performed locally using the *parent* indirection.

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



Python API
----------

.. ifconfig:: enable_python_doc == 'ON'

  .. todo::
    ...


.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)
