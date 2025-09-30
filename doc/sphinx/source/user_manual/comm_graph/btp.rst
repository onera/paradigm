.. _block_to_part:

Block to Part
=============

Description
"""""""""""

**Block to Part** is a service for managing MPI data transfers from block-distributed to partitioned sets of entities.

API
"""

.. dropdown:: Initialization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_block_to_part_create
      .. doxygenfunction:: PDM_block_to_part_create_from_sparse_block


    .. tab-item:: Fortran
      :sync: Fortran

      .. f:autosubroutine:: PDM_block_to_part_create


    .. tab-item:: Python
      :sync: Python

        .. py:class:: BlockToPart

          .. automethod:: Pypdm.Pypdm.BlockToPart.__init__





.. dropdown:: Exchange

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_block_to_part_exch
      .. doxygenfunction:: PDM_block_to_part_exch_in_place


    .. tab-item:: Fortran
      :sync: Fortran

      .. f:subroutine:: pdm_block_to_part_exch(btp, t_stride, block_stride, block_data, part_stride, part_data)

        Exchange data from blocks to partitions. Output arrays are allocated by **ParaDiGM**.

        :p c_ptr btp [in]: Block-to-Part instance
        :p integer t_stride [in]: Stride type
        :p pointer(integer(pdm_l_num_s)) block_stride [in]: Stride for ``block_data``
        :p pointer block_data [in]: Block-distributed data
        :p pdm_pointer_array_t part_stride [out]: Stride for ``part_data``
        :p pdm_pointer_array_t part_data [out]: Partitioned data


      .. f:subroutine:: pdm_block_to_part_exch_in_place(btp, t_stride, block_stride, block_data, part_stride, part_data)

        Exchange data from blocks to partitions. Output arrays are allocated by the user *before* the call to this subroutine.

        :p c_ptr btp [in]: Block-to-Part instance
        :p integer t_stride [in]: Stride type
        :p pointer(integer(pdm_l_num_s)) block_stride [in]: Stride for ``block_data``
        :p pointer block_data [in]: Block-distributed data
        :p pdm_pointer_array_t part_stride [inout]: Stride for ``part_data``
        :p pdm_pointer_array_t part_data [inout]: Partitioned data



    .. tab-item:: Python
      :sync: Python

        .. automethod:: Pypdm.Pypdm.BlockToPart.exchange_field
        .. automethod:: Pypdm.Pypdm.BlockToPart.exchange_field_inplace




.. dropdown:: Finalization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

        .. doxygenfunction:: PDM_block_to_part_free


    .. tab-item:: Fortran
      :sync: Fortran

        .. f:autosubroutine:: PDM_block_to_part_free



    .. tab-item:: Python
      :sync: Python

      |python_gc|

