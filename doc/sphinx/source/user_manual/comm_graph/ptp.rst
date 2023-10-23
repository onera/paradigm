.. _ptp:

Part to Part
============

C API
-----

Initialization
""""""""""""""

.. doxygenfunction:: PDM_part_to_part_create

.. doxygenfunction:: PDM_part_to_part_create_from_num2_triplet


Information on Part 2 side
""""""""""""""""""""""""""

.. doxygenfunction:: PDM_part_to_part_ref_lnum2_get

.. doxygenfunction:: PDM_part_to_part_unref_lnum2_get

.. doxygenfunction:: PDM_part_to_part_gnum1_come_from_get


Exchange
""""""""

.. doxygenfunction:: PDM_part_to_part_iexch

.. doxygenfunction:: PDM_part_to_part_iexch_wait

.. doxygenfunction:: PDM_part_to_part_reverse_iexch

.. doxygenfunction:: PDM_part_to_part_reverse_iexch_wait

.. .. doxygenfunction:: PDM_part_to_part_issend

.. .. doxygenfunction:: PDM_part_to_part_issend_wait

.. .. doxygenfunction:: PDM_part_to_part_reverse_issend

.. .. doxygenfunction:: PDM_part_to_part_reverse_issend_wait

.. .. doxygenfunction:: PDM_part_to_part_irecv

.. .. doxygenfunction:: PDM_part_to_part_irecv_wait

.. .. doxygenfunction:: PDM_part_to_part_reverse_irecv

.. .. doxygenfunction:: PDM_part_to_part_reverse_irecv_wait


Finalization
""""""""""""

.. doxygenfunction:: PDM_part_to_part_free



Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  Initialization
  """"""""""""""

  .. f:autosubroutine PDM_part_to_part_create

  Information on Part 2 side
  """"""""""""""""""""""""""

  .. f:autosubroutine PDM_part_to_part_ref_lnum2_get

  .. f:autosubroutine PDM_part_to_part_unref_lnum2_get

  .. f:autosubroutine PDM_part_to_part_gnum1_come_from_get

  Exchange
  """"""""

  .. f:autosubroutine PDM_part_to_part_iexch

  .. f:autosubroutine PDM_part_to_part_iexch_wait

  .. f:autosubroutine PDM_part_to_part_reverse_iexch

  .. f:autosubroutine PDM_part_to_part_reverse_iexch_wait

  .. f:autosubroutine PDM_part_to_part_issend

  .. f:autosubroutine PDM_part_to_part_issend_wait

  .. f:autosubroutine PDM_part_to_part_irecv_raw

  .. f:autosubroutine PDM_part_to_part_irecv_wait_raw

  Finalization
  """"""""""""

  .. f:autosubroutine PDM_part_to_part_free

.. ifconfig:: enable_fortran_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



Python API
----------

.. ifconfig:: enable_python_doc == 'ON'

  Initialization
  """"""""""""""

  .. autoclass:: Pypdm.Pypdm.PartToPart


  Information on Part 2 side
  """"""""""""""""""""""""""

  .. autofunction:: Pypdm.Pypdm.PartToPart.get_referenced_lnum2

  .. autofunction:: Pypdm.Pypdm.PartToPart.get_unreferenced_lnum2

  .. autofunction:: Pypdm.Pypdm.PartToPart.get_gnum1_come_from


  Exchange
  """"""""

  .. autofunction:: Pypdm.Pypdm.PartToPart.iexch

  .. autofunction:: Pypdm.Pypdm.PartToPart.wait

  .. autofunction:: Pypdm.Pypdm.PartToPart.reverse_iexch

  .. autofunction:: Pypdm.Pypdm.PartToPart.reverse_wait


.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)
