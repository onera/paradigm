.. _ptp:

Part to Part
============

Description
"""""""""""

**Part to Part** is a service for managing MPI data transfers between two arbitrarily partitioned sets of entities.


API
"""

.. dropdown:: Initialization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_to_part_create
      .. doxygenfunction:: PDM_part_to_part_create_from_num2_triplet



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_to_part_create

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.PartToPart.__init__

        .. autofunction:: Pypdm.Pypdm.PartToPart.from_triplet

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Information on Part 2 side

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_to_part_ref_lnum2_get
      .. doxygenfunction:: PDM_part_to_part_unref_lnum2_get
      .. doxygenfunction:: PDM_part_to_part_gnum1_come_from_get



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_to_part_ref_lnum2_get
        .. f:autosubroutine:: PDM_part_to_part_unref_lnum2_get
        .. f:autosubroutine:: PDM_part_to_part_gnum1_come_from_get

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.PartToPart.get_referenced_lnum2

        .. autofunction:: Pypdm.Pypdm.PartToPart.get_unreferenced_lnum2

        .. autofunction:: Pypdm.Pypdm.PartToPart.get_gnum1_come_from


      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Exchange

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_to_part_iexch
      .. doxygenfunction:: PDM_part_to_part_iexch_wait
      .. doxygenfunction:: PDM_part_to_part_reverse_iexch
      .. doxygenfunction:: PDM_part_to_part_reverse_iexch_wait


    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_to_part_iexch
        .. f:autosubroutine:: PDM_part_to_part_iexch_wait

        .. f:autosubroutine:: PDM_part_to_part_reverse_iexch
        .. f:autosubroutine:: PDM_part_to_part_reverse_iexch_wait

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.PartToPart.iexch
        .. autofunction:: Pypdm.Pypdm.PartToPart.wait
        .. autofunction:: Pypdm.Pypdm.PartToPart.reverse_iexch
        .. autofunction:: Pypdm.Pypdm.PartToPart.reverse_wait

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)



.. dropdown:: Finalization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_to_part_free



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_to_part_free

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      |python_gc|
