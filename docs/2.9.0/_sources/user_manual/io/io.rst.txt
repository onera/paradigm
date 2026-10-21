.. _io:

IO
==

Description
"""""""""""

**IO** is a service for managing parallel file reading and writing.

API
"""

.. dropdown:: Create a directory

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_io_mkdir


    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_io_mkdir

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



.. dropdown:: Open a file

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_io_open

      .. doxygenenum:: PDM_io_type_t
      .. doxygenenum:: PDM_io_suff_t
      .. doxygenenum:: PDM_io_kind_t
      .. doxygenenum:: PDM_io_mod_t
      .. doxygenenum:: PDM_io_fmt_t
      .. doxygenenum:: PDM_io_backup_t


    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_io_open

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



.. dropdown:: General file information

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_io_file_name_get
      .. doxygenfunction:: PDM_io_comm_get
      .. doxygenfunction:: PDM_io_dump


    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_io_comm_get
        .. f:autosubroutine:: PDM_io_dump

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



.. dropdown:: Manage cursor in a file

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_io_seek
      .. doxygenfunction:: PDM_io_tell


    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_io_seek
        .. f:autosubroutine:: PDM_io_tell

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



.. dropdown:: Read data from a file

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_io_global_read
      .. doxygenfunction:: PDM_io_par_interlaced_read
      .. doxygenfunction:: PDM_io_par_block_read


    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_io_global_read
        .. f:autosubroutine:: PDM_io_par_interlaced_read
        .. f:autosubroutine:: PDM_io_par_block_read

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



.. dropdown:: Write data in a file

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_io_global_write
      .. doxygenfunction:: PDM_io_par_interlaced_write
      .. doxygenfunction:: PDM_io_par_block_write


    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_io_global_write
        .. f:autosubroutine:: PDM_io_par_interlaced_write
        .. f:autosubroutine:: PDM_io_par_block_write

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



.. dropdown:: Little/big endian swapping

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_io_swap_endian_on
      .. doxygenfunction:: PDM_io_swap_endian_off
      .. doxygenfunction:: PDM_io_swap_endian

      .. doxygenenum:: PDM_io_endian_t


    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_io_swap_endian_on
        .. f:autosubroutine:: PDM_io_swap_endian_off
        .. f:autosubroutine:: PDM_io_swap_endian

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



.. dropdown:: Manage data and format

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_io_n_data_get
      .. doxygenfunction:: PDM_io_fmt_data_set

      .. doxygenenum:: PDM_io_endian_t


    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_io_n_data_get
        .. f:autosubroutine:: PDM_io_fmt_data_set

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



.. dropdown:: Close a file

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_io_close


    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_io_close

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



.. dropdown:: Get elapsed and CPU times

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_io_get_timer_fichier
      .. doxygenfunction:: PDM_io_timer_swap_endian_get
      .. doxygenfunction:: PDM_io_timer_distrib_get
      .. doxygenfunction:: PDM_io_timer_total_get


    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_io_get_timer_fichier
        .. f:autosubroutine:: PDM_io_timer_swap_endian_get
        .. f:autosubroutine:: PDM_io_timer_distrib_get
        .. f:autosubroutine:: PDM_io_timer_total_get

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



.. dropdown:: Finalization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_io_free


    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_io_free

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



Examples
""""""""

.. dropdown:: Write

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. code-block:: c

        // Initialize write
        PDM_io_file_t *pdm_file = NULL;
        PDM_l_num_t    ierr;
        PDM_io_open(filename,
                    PDM_IO_FMT_BIN,         // Binary output file
                    PDM_IO_SUFF_MAN,        // Manual suffix
                    "",                     // Empty suffix
                    PDM_IO_BACKUP_OFF,      // No backup
                    PDM_IO_KIND_MPI_SIMPLE, // Simple MPI-IO access
                    PDM_IO_MOD_WRITE,       // Write mode
                    PDM_IO_NATIVE,          // Native endian
                    comm,
                    -1.,
                    &pdm_file,
                    &ierr);

        // To write a keyword in a mesh file
        // buffer is a char* containing the keyword
        PDM_io_global_write(pdm_file,
              (PDM_l_num_t) sizeof(char),
              (PDM_l_num_t) s_buffer,
                            buffer);

        // To make each MPI rank write a buffer
        PDM_io_par_interlaced_write(pdm_file,
                                    PDM_STRIDE_VAR_INTERLACED,
                    (PDM_l_num_t *) &size, // Buffer size in number of characters
                    (PDM_l_num_t  ) sizeof(char),
                    (PDM_l_num_t  ) 1,
                                    &i_rank_gnum, // MPI rank write order
                     (const void *) buffer);

        // End write
        PDM_io_close(pdm_file);
        PDM_io_free (pdm_file);


    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. code-block:: fortran

          call PDM_io_open(filename,               &
                           PDM_IO_FMT_BIN,         & ! Binary output file
                           PDM_IO_SUFF_MAN,        & ! Manual suffix
                           "",                     & ! Empty suffix
                           PDM_IO_BACKUP_OFF,      & ! No backup
                           PDM_IO_KIND_MPI_SIMPLE, & ! Simple MPI-IO access
                           PDM_IO_MOD_WRITE,       & ! Write mode
                           PDM_IO_NATIVE,          & ! Native endian
                           comm,                   &
                           -1.,                    &
                           pdm_file,               &
                           ierr)

          ! To write a keyword in a mesh file
          ! buffer is a character(*) containing the keyword
          call PDM_io_global_write(pdm_file, &
                                   4,        &
                                   s_buffer, &
                                   buffer)

          ! To make each MPI rank write a buffer
          call PDM_io_par_interlaced_write(pdm_file,                  &
                                           PDM_STRIDE_VAR_INTERLACED, &
                                           size,                      & ! Buffer size in number of characters
                                           4,                         & ! Number of bytes per character
                                           1,                         &
                                           indirection,               & ! MPI rank write order
                                           buffer)

          ! End write
          call PDM_io_close(pdm_file)
          call PDM_io_free (pdm_file)

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)

