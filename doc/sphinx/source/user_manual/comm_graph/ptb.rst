.. _part_to_block:

Part to Block
=============

Description
"""""""""""

**Part to Block** is a service for managing MPI data transfers from partitioned to block-distributed sets of entities.

API
"""

.. dropdown:: Initialization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_to_block_create
      .. doxygenfunction:: PDM_part_to_block_create_from_distrib
      .. doxygenfunction:: PDM_part_to_block_geom_create


    .. tab-item:: Fortran
      :sync: Fortran

      .. f:autosubroutine:: PDM_part_to_block_create
      .. f:autosubroutine:: PDM_part_to_block_create_from_distrib



    .. tab-item:: Python
      :sync: Python

        .. py:class:: PartToBlock

          .. autofunction:: Pypdm.Pypdm.PartToBlock.__init__



.. dropdown:: Information on block-distributed frame

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_to_block_n_elt_block_get
      .. doxygenfunction:: PDM_part_to_block_block_gnum_get
      .. doxygenfunction:: PDM_part_to_block_block_gnum_count_get
      .. doxygenfunction:: PDM_part_to_block_distrib_index_get


    .. tab-item:: Fortran
      :sync: Fortran

      .. f:autofunction::   PDM_part_to_block_n_elt_block_get
      .. f:autosubroutine:: PDM_part_to_block_block_gnum_get
      .. f:autosubroutine:: PDM_part_to_block_distrib_index_get

    .. tab-item:: Python
      :sync: Python

      .. autofunction:: Pypdm.Pypdm.PartToBlock.getBlockGnumCopy



.. dropdown:: Exchange

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_to_block_exch
      .. doxygenfunction:: PDM_part_to_block_reverse_exch
      .. doxygenfunction:: PDM_part_to_block_iexch
      .. doxygenfunction:: PDM_part_to_block_iexch_wait
      .. doxygenfunction:: PDM_part_to_block_reverse_iexch
      .. doxygenfunction:: PDM_part_to_block_reverse_iexch_wait


    .. tab-item:: Fortran
      :sync: Fortran

      .. f:subroutine:: pdm_part_to_block_exch(ptb, t_stride, cst_stride, part_stride, part_data, block_stride, block_data)

        Exchange data from partitions to blocks. Output arrays are allocated by **ParaDiGM**.

        :p c_ptr ptb [in]: Part-to-Block instance
        :p integer t_stride [in]: Stride type
        :p integer cst_stride [in]: Constant stride value
        :p pdm_pointer_array_t part_stride [in]: Stride for ``part_data``
        :p pdm_pointer_array_t part_data [in]: Partitioned data
        :p pointer(integer(pdm_l_num_s)) block_stride [out]: Stride for ``block_data``
        :p pointer block_data [out]: Block-distributed data


    .. tab-item:: Python
      :sync: Python

        .. autofunction:: Pypdm.Pypdm.PartToBlock.exchange_field



.. dropdown:: Finalization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_part_to_block_free


    .. tab-item:: Fortran
      :sync: Fortran

      .. f:autosubroutine:: PDM_part_to_block_free


    .. tab-item:: Python
      :sync: Python

      |python_gc|


Examples
""""""""