.. _ptp:

Part to Part
============


API
"""

.. dropdown:: Initialization

  .. tab-set::

    .. tab-item:: C

      .. doxygenfunction:: PDM_part_to_part_create
      .. doxygenfunction:: PDM_part_to_part_create_from_num2_triplet



    .. tab-item:: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_to_part_create

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.PartToPart.__init__
          :noindex:

        .. autofunction:: Pypdm.Pypdm.PartToPart.from_triplet
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Information on Part 2 side

  .. tab-set::

    .. tab-item:: C

      .. doxygenfunction:: PDM_part_to_part_ref_lnum2_get
      .. doxygenfunction:: PDM_part_to_part_unref_lnum2_get
      .. doxygenfunction:: PDM_part_to_part_gnum1_come_from_get



    .. tab-item:: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_to_part_ref_lnum2_get
        .. f:autosubroutine:: PDM_part_to_part_unref_lnum2_get
        .. f:autosubroutine:: PDM_part_to_part_gnum1_come_from_get

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.PartToPart.get_referenced_lnum2
          :noindex:

        .. autofunction:: Pypdm.Pypdm.PartToPart.get_unreferenced_lnum2
          :noindex:

        .. autofunction:: Pypdm.Pypdm.PartToPart.get_gnum1_come_from
          :noindex:


      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Exchange

  .. tab-set::

    .. tab-item:: C

      .. doxygenfunction:: PDM_part_to_part_iexch
      .. doxygenfunction:: PDM_part_to_part_iexch_wait
      .. doxygenfunction:: PDM_part_to_part_reverse_iexch
      .. doxygenfunction:: PDM_part_to_part_reverse_iexch_wait


    .. tab-item:: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_part_to_part_iexch

        .. f:subroutine:: pdm_part_to_part_iexch_wait(ptp, request)

          Finalize a non-blocking exchange (Part1→Part2)

          :p c_ptr ptp [in]:       Part-to-Part instance
          :p integer request [in]: Request

        .. f:autosubroutine:: PDM_part_to_part_reverse_iexch

        .. f:subroutine:: pdm_part_to_part_reverse_iexch_wait(ptp, request)

          Finalize a non-blocking exchange (Part2→Part1)

          :p c_ptr ptp [in]:       Part-to-Part instance
          :p integer request [in]: Request


      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.PartToPart.iexch
          :noindex:
        .. autofunction:: Pypdm.Pypdm.PartToPart.wait
          :noindex:
        .. autofunction:: Pypdm.Pypdm.PartToPart.reverse_iexch
          :noindex:
        .. autofunction:: Pypdm.Pypdm.PartToPart.reverse_wait
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)



.. dropdown:: Finalization

  .. tab-set::

    .. tab-item:: C

      .. doxygenfunction:: PDM_part_to_part_free



    .. tab-item:: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:subroutine:: pdm_part_to_part_free(ptp)

          Free a Part-to-Part structure

          :p c_ptr ptp [inout]: Part-to-part instance

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python

      The instance is automatically freed by Python's garbage collector when it is no longer referenced
