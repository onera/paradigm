.. _part_extension:

Part extension
==============

Description
"""""""""""

**Part extension** is a service for computing extended partitions by fetching ghost cells topology and geometry.
If provided, group information associated with extended entities is fetched as well.

Extension to an arbitrary depth from vertices, edges or faces is supported :

.. _PDM_extend_type_t:

.. doxygenenum:: PDM_extend_type_t

API
"""

.. dropdown:: Initialization


  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_extension_create



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_extension_create

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. py:class:: PartExtension

          .. autofunction:: Pypdm.Pypdm.PartExtension.__init__

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Set input mesh


  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_extension_connectivity_set
      .. doxygenfunction:: PDM_part_extension_vtx_coord_set
      .. doxygenfunction:: PDM_part_extension_ln_to_gn_set
      .. doxygenfunction:: PDM_part_extension_part_bound_graph_set
      .. doxygenfunction:: PDM_part_extension_group_set



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_extension_connectivity_set
        .. f:autosubroutine:: PDM_part_extension_vtx_coord_set
        .. f:autosubroutine:: PDM_part_extension_ln_to_gn_set
        .. f:autosubroutine:: PDM_part_extension_part_bound_graph_set
        .. f:autosubroutine:: PDM_part_extension_group_set

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.PartExtension.connectivity_set
        .. autofunction:: Pypdm.Pypdm.PartExtension.vtx_coord_set
        .. autofunction:: Pypdm.Pypdm.PartExtension.ln_to_gn_set
        .. autofunction:: Pypdm.Pypdm.PartExtension.part_bound_graph_set
        .. autofunction:: Pypdm.Pypdm.PartExtension.group_set

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Compute extended partitions

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_extension_compute



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_extension_compute

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.PartExtension.compute

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Get extended entities


  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_extension_connectivity_get
      .. doxygenfunction:: PDM_part_extension_vtx_coord_get
      .. doxygenfunction:: PDM_part_extension_ln_to_gn_get
      .. doxygenfunction:: PDM_part_extension_group_get



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_extension_connectivity_get
        .. f:autosubroutine:: PDM_part_extension_vtx_coord_get
        .. f:autosubroutine:: PDM_part_extension_ln_to_gn_get
        .. f:autosubroutine:: PDM_part_extension_group_get

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.PartExtension.connectivity_get
        .. autofunction:: Pypdm.Pypdm.PartExtension.vtx_coord_get
        .. autofunction:: Pypdm.Pypdm.PartExtension.ln_to_gn_get
        .. autofunction:: Pypdm.Pypdm.PartExtension.group_get

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Finalization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_extension_free



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_extension_free

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      |python_gc|
