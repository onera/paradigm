.. _mpi_io:

MPI-IO
======

C API
-----

Enumerators
~~~~~~~~~~~

.. doxygenenum:: PDM_io_type_t

.. doxygenenum:: PDM_io_suff_t

.. doxygenenum:: PDM_io_kind_t

.. doxygenenum:: PDM_io_mod_t

.. doxygenenum:: PDM_io_fmt_t

.. doxygenenum:: PDM_io_backup_t

.. doxygenenum:: PDM_io_endian_t

.. doxygenenum:: PDM_io_seek_t

Initialization
~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_io_open

.. doxygenfunction:: PDM_io_swap_endian

.. doxygenfunction:: PDM_io_fmt_data_set

.. doxygenfunction:: PDM_io_mkdir

IO information
~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_io_file_name_get

.. doxygenfunction:: PDM_io_seek

.. doxygenfunction:: PDM_io_tell

.. doxygenfunction:: PDM_io_dump

.. doxygenfunction:: PDM_io_comm_get

.. doxygenfunction:: PDM_io_swap_endian_on

.. doxygenfunction:: PDM_io_swap_endian_off

.. doxygenfunction:: PDM_io_n_data_get

Read
~~~~

.. doxygenfunction:: PDM_io_global_read

.. doxygenfunction:: PDM_io_par_interlaced_read

.. doxygenfunction:: PDM_io_par_block_read

Write
~~~~~

.. doxygenfunction:: PDM_io_global_write

.. doxygenfunction:: PDM_io_par_interlaced_write

.. doxygenfunction:: PDM_io_par_block_write

Finalize
~~~~~~~~

.. doxygenfunction:: PDM_io_close

.. doxygenfunction:: PDM_io_free

Timers
~~~~~~

.. doxygenfunction:: PDM_io_get_timer_fichier

.. doxygenfunction:: PDM_io_timer_swap_endian_get

.. doxygenfunction:: PDM_io_timer_distrib_get

.. doxygenfunction:: PDM_io_timer_total_get

Example of usage
~~~~~~~~~~~~~~~~

.. code-block:: c

  // Initialize write
  PDM_io_file_t *writer = NULL;
  PDM_l_num_t    ierr;
  PDM_io_open(filename,
              PDM_IO_FMT_BIN, // Binary output file
              PDM_IO_SUFF_MAN,
              "",
              PDM_IO_BACKUP_OFF,
              PDM_IO_KIND_MPI_SIMPLE,
              PDM_IO_MOD_WRITE,
              PDM_IO_NATIVE,
              comm,
              -1.,
              &writer,
              &ierr);

  // To write a keyword in a mesh file
  // buffer is a char* containing the keyword
  PDM_io_global_write(writer,
                      (PDM_l_num_t) sizeof(char),
                      (PDM_l_num_t) s_buffer,
                      buffer);

  // To make each MPI rank write a buffer
  PDM_io_par_interlaced_write(writer,
                              PDM_STRIDE_VAR_INTERLACED,
              (PDM_l_num_t *) &size, // Buffer size in number of characters
                (PDM_l_num_t) sizeof(char),
                (PDM_l_num_t) one,
                              &i_rank_gnum, // MPI rank write order
               (const void *) buffer);

  // End write
  PDM_io_close(writer);
  PDM_io_free(writer);

Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  Initialization
  ~~~~~~~~~~~~~~

  .. f:subroutine:: PDM_io_open(def)

    Open a file for parallel access

    :param character        nom                [in]:  File name
    :param integer          fmt                [in]:  Text of Binary format
    :param integer          suff_t             [in]:  Type of suffix (manual or automatic)
    :param character        suff_u             [in]:  Suffix (if manual)
    :param integer          s_backup           [in]:  Activates the backup of a pre-existing file in write mode
    :param integer          acces              [in]:  Type (parallel with MPI-IO, parallel without MPI-IO, sequential)
    :param integer          mode               [in]:  Access mode (read, write, read/write)
    :param integer          endian             [in]:  Endian type (big, little or native)
    :param integer          comm               [in]:  Communicator associated to the file
    :param real             prop_noeuds_actifs [in]:  Proportion of active nodes
    :param c_ptr            unite              [out]: Unit of the file
    :param integer          ierr               [out]: Indicates whether the file is of type PDM_io or not (for read-only opening only)

  .. f:subroutine:: PDM_io_swap_endian(def)

    Swap endian pour conversion little endian <-> big endian

    :param integer taille_donnee [in]:  Size of a unit piece of data
    :param integer n_donnees     [in]:  Amount of data
    :param c_ptr   donnees       [in]:  Data
    :param c_ptr   resultats     [out]: Result

  .. f:subroutine:: PDM_io_fmt_data_set(def)

    Defines the format of the individual data for text output

    :param c_ptr     fichier    [in]: Pointer to \ref PDM_io_file_t object
    :param integer   n_char_fmt [in]: Number of characters in the format
    :param integer   data_type  [in]: Type of data
    :param character fmt        [in]: Format

  .. f:subroutine:: PDM_io_mkdir(def)

    Create a directory

    :param character path [in]:  Path to new directory
    :param integer   code [out]: 0 if successful, -1 else

  IO information
  ~~~~~~~~~~~~~~

  .. f:subroutine:: PDM_io_seek(def)

    Set the file position indicator

    :param c_ptr   fichier [in]: Pointer to \ref PDM_io_file_t object
    :param integer offset  [in]: Address
    :param integer seek    [in]: Origin type

  .. f:subroutine:: PDM_io_tell(def)

    Return the current file position

    :param c_ptr   fichier [in]:  Pointer to \ref PDM_io_file_t object
    :param integer offset  [out]: Current position in file

  .. f:subroutine:: PDM_io_dump(def)

    Shows file information

    :param c_ptr fichier [in]: Pointer to \ref PDM_io_file_t object

  .. f:subroutine:: PDM_io_comm_get(def)

    Returns the file communicator

    :param c_ptr   fichier [in]:  Pointer to \ref PDM_io_file_t object
    :param integer f_comm  [out]: MPI communicator

  .. f:subroutine:: PDM_io_swap_endian_on(def)

    Activate endian swap

    :param c_ptr fichier [in]: Pointer to \ref PDM_io_file_t object

  .. f:subroutine:: PDM_io_swap_endian_off(def)

    Deactivate endian swap

    :param c_ptr fichier [in]: Pointer to \ref PDM_io_file_t object

  .. f:subroutine:: PDM_io_n_data_get(def)

    Calculating the total size of a data field

    :param c_ptr      fichier         [in]:  Pointer to \ref PDM_io_file_t object
    :param integer    t_n_composantes [in]:  Type of component sizes (PDM_STRIDE_CST_INTERLACED or PDM_STRIDE_VAR_INTERLACED)
    :param integer(:) n_composantes   [in]:  Number of components for each data
    :param integer    n_donnees       [in]:  Number of data
    :param integer(:) indirection     [in]:  Data redistribution direction
    :param integer    taille          [out]: Total size of a data field

  Read
  ~~~~

  .. f:subroutine:: PDM_io_global_read(def)

    Global read: the master process alone accesses the
    file and redistributes the information to all the communicator's processes

    :param c_ptr   fichier       [in]:  Pointer to \ref PDM_io_file_t object
    :param integer taille_donnee [in]:  Size of a unit piece of data
    :param integer n_donnees     [in]:  Amount of data to be read
    :param c_ptr   donnees       [out]: Read data

  .. f:subroutine:: PDM_io_par_interlaced_read(def)

    Parallel reading of data blocks followed by
    redistribution of the data according to indirection

    :param c_ptr      fichier         [in]:  Pointer to \ref PDM_io_file_t object
    :param integer    t_n_composantes [in]:  Type of component sizes (PDM_STRIDE_CST_INTERLACED or PDM_STRIDE_VAR_INTERLACED)
    :param integer(:) n_composantes   [in]:  Number of components for each piece of data
    :param integer    taille_donnee   [in]:  Unit size of a piece of data
    :param integer    n_donnees       [in]:  Number of data items to be read
    :param integer(:) indirection     [in]:  Indirection of data redistribution
    :param c_ptr      donnees         [out]: Read data

  .. f:subroutine:: PDM_io_par_block_read(def)

    Parallel reading of data blocks
    The blocks must be arranged in ascending order
    according to the numbering of the processes

    :param c_ptr      fichier         [in]:  Pointer to \ref PDM_io_file_t object
    :param integer    t_n_composantes [in]:  Component size type (PDM_STRIDE_CST_INTERLACED or PDM_STRIDE_VAR_INTERLACED)
    :param integer(:) n_composantes   [in]:  Number of components for each data item
    :param integer    taille_donnee   [in]:  Unit size of a piece of data
    :param integer    n_donnees       [in]:  Number of data items to be read
    :param integer    debut_bloc      [in]:  Relative address of start of block
    :param c_ptr      donnees         [out]: Read data

  Write
  ~~~~~

  .. f:subroutine:: PDM_io_global_write(def)

    Global write: The master process has sole access to the file

    :param c_ptr   fichier       [in]: Pointer to \ref PDM_io_file_t object
    :param integer taille_donnee [in]: Size of a unit piece of data
    :param integer n_donnees     [in]: Amount of data to write
    :param c_ptr   donnees       [in]: Data to write

  .. f:subroutine:: PDM_io_par_interlaced_write(def)

    Data sorted according to indirection, then parallel write of data blocks

    :param c_ptr      fichier         [in]: Pointer to \ref PDM_io_file_t object
    :param integer    t_n_composantes [in]: Type of component sizes (PDM_STRIDE_CST_INTERLACED or PDM_STRIDE_VAR_INTERLACED)
    :param integer(:) n_composantes   [in]: Number of components for each data item
    :param integer    taille_donnee   [in]: Unit size of the data
    :param integer    n_donnees       [in]: Number of data items to be written
    :param integer(:) indirection     [in]: Data redistribution direction
    :param c_ptr      donnees         [in]: Data to be written

  .. f:subroutine:: PDM_io_par_block_write(def)

    Parallel writing of data blocks
    Blocks must be arranged in ascending order according
    to numbering of the processes

    :param c_ptr      fichier         [in]: Pointer to \ref PDM_io_file_t object
    :param integer    t_n_composantes [in]: Type of component sizes (PDM_STRIDE_CST_INTERLACED or PDM_STRIDE_VAR_INTERLACED)
    :param integer(:) n_composantes   [in]: Number of components for each data item
    :param integer    taille_donnee   [in]: Unit size of the data
    :param integer    n_donnees       [in]: Number of data to read
    :param integer    debut_bloc      [in]: Relative address of start of block
    :param c_ptr      donnees         [in]: Data to be written

  Finalize
  ~~~~~~~~

  .. f:subroutine:: PDM_io_close(def)

    Closing the file without destroying the PDM_io structure
    associated with unit

    :param c_ptr fichier [in]: Pointer to \ref PDM_io_file_t object

  .. f:subroutine:: PDM_io_free(def)

    Free of the PDM_io structure associated with the unit

    :param c_ptr fichier [in]: Pointer to \ref PDM_io_file_t object

  Timers
  ~~~~~~

  .. f:subroutine:: PDM_io_get_timer_fichier(def)

    Returns the cumulative files access time

    :param c_ptr fichier   [in]:  Pointer to \ref PDM_io_file_t object
    :param real  t_cpu     [out]: CPU time
    :param real  t_elapsed [out]: Elapsed time

  .. f:subroutine:: PDM_io_timer_swap_endian_get(def)

    Returns the cumulative time for data swap

    :param c_ptr fichier   [in]:  Pointer to \ref PDM_io_file_t object
    :param real  t_cpu     [out]: CPU time
    :param real  t_elapsed [out]: Elapsed time

  .. f:subroutine:: PDM_io_timer_distrib_get(def)

    Returns the cumulative time for data distribution

    :param c_ptr fichier   [in]:  Pointer to \ref PDM_io_file_t object
    :param real  t_cpu     [out]: CPU time
    :param real  t_elapsed [out]: Elapsed time

  .. f:subroutine:: PDM_io_timer_total_get(def)

    Returns the total cumulative time

    :param c_ptr fichier   [in]:  Pointer to \ref PDM_io_file_t object
    :param real  t_cpu     [out]: CPU time
    :param real  t_elapsed [out]: Elapsed time

.. ifconfig:: enable_fortran_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)

