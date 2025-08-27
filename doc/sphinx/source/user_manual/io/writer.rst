.. _writer:

Writer
======

Description
"""""""""""

**Writer** is a service for writing meshes and associated fields in parallel using MPI-IO.


API
"""

.. dropdown:: Initialization


  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_writer_create



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_writer_create

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.Writer.__init__

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)

  |

  **Writer** supports binary or ASCII writing:

  .. doxygenenum:: PDM_writer_fmt_fic_t

  |

  Meshes with time-dependent geometry and topology are also supported:

  .. doxygenenum:: PDM_writer_topology_t


.. dropdown:: Define geometry

  If a :ref:`Part Mesh Nodal <part_mesh_nodal>` instance is available, it can be used directly to define the geometry.

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_writer_geom_create_from_mesh_nodal



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_writer_geom_create_from_mesh_nodal

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)




  |

  Otherwise, one geometry instance must be created for each dimension present in the mesh.

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_writer_geom_create



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_writer_geom_create

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.Writer.geom_create

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)


  |


  .. dropdown:: Set vertices

    .. tab-set::
      :sync-group: language

      .. tab-item:: C
        :sync: C

        .. doxygenfunction:: PDM_writer_geom_coord_set
        .. doxygenfunction:: PDM_writer_geom_coord_from_parent_set



      .. tab-item:: Fortran
        :sync: Fortran

        .. ifconfig:: enable_fortran_doc == 'ON'

          .. f:autosubroutine:: PDM_writer_geom_coord_set
          .. f:autosubroutine:: PDM_writer_geom_coord_from_parent_set

        .. ifconfig:: enable_fortran_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



      .. tab-item:: Python
        :sync: Python

        .. ifconfig:: enable_python_doc == 'ON'

          .. autofunction:: Pypdm.Pypdm.Writer.geom_coord_set


        .. ifconfig:: enable_python_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)





  .. dropdown:: Set elements

    The mesh elements can be defined in multiple ways.

    If only a soup of mixed elements is available, it can be set using the appropriate function.

    .. dropdown:: Mixed setters

      .. tab-set::
        :sync-group: language

        .. tab-item:: C
          :sync: C

          .. doxygenfunction:: PDM_writer_geom_cell3d_cellface_add
          .. doxygenfunction:: PDM_writer_geom_cell2d_cellface_add
          .. doxygenfunction:: PDM_writer_geom_faces_facesom_add



        .. tab-item:: Fortran
          :sync: Fortran

          .. ifconfig:: enable_fortran_doc == 'ON'

            .. f:autosubroutine:: PDM_writer_geom_cell3d_cellface_add
            .. f:autosubroutine:: PDM_writer_geom_cell2d_cellface_add
            .. f:autosubroutine:: PDM_writer_geom_faces_facesom_add

          .. ifconfig:: enable_fortran_doc == 'OFF'

            .. warning::
              Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



        .. tab-item:: Python
          :sync: Python

          .. ifconfig:: enable_python_doc == 'ON'

            .. autofunction:: Pypdm.Pypdm.Writer.geom_cell3d_cellface_add
            .. autofunction:: Pypdm.Pypdm.Writer.geom_cell2d_cellface_add
            .. autofunction:: Pypdm.Pypdm.Writer.geom_faces_facevtx_add


    |

    If the mesh elements are already separated by types, sections can be added and set one by one:

    .. dropdown:: Section setters

      .. tab-set::
        :sync-group: language

        .. tab-item:: C
          :sync: C

          .. doxygenfunction:: PDM_writer_geom_bloc_add
          .. doxygenfunction:: PDM_writer_geom_bloc_std_set
          .. doxygenfunction:: PDM_writer_geom_bloc_poly2d_set
          .. doxygenfunction:: PDM_writer_geom_bloc_poly3d_set



        .. tab-item:: Fortran
          :sync: Fortran

          .. ifconfig:: enable_fortran_doc == 'ON'

            .. f:autosubroutine:: PDM_writer_geom_bloc_add
            .. f:autosubroutine:: PDM_writer_geom_bloc_std_set
            .. f:autosubroutine:: PDM_writer_geom_bloc_poly2d_set
            .. f:autosubroutine:: PDM_writer_geom_bloc_poly3d_set

          .. ifconfig:: enable_fortran_doc == 'OFF'

            .. warning::
              Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



        .. tab-item:: Python
          :sync: Python

          .. ifconfig:: enable_python_doc == 'ON'

            .. autofunction:: Pypdm.Pypdm.Writer.geom_block_add
            .. autofunction:: Pypdm.Pypdm.Writer.geom_block_std_set


          .. ifconfig:: enable_python_doc == 'OFF'

            .. warning::
              Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)



.. dropdown:: Define variables

  **Writer** supports both *local* variables (i.e. vertex- of element-based fields), and *global* variables (i.e. that depend solely on time).

  .. dropdown:: Local variables

    .. tab-set::
      :sync-group: language

      .. tab-item:: C
        :sync: C

        .. doxygenfunction:: PDM_writer_var_create
        .. doxygenfunction:: PDM_writer_name_map_add
        .. doxygenfunction:: PDM_writer_var_set



      .. tab-item:: Fortran
        :sync: Fortran

        .. ifconfig:: enable_fortran_doc == 'ON'

          .. f:autosubroutine:: PDM_writer_var_create
          .. f:autosubroutine:: PDM_writer_name_map_add
          .. f:autosubroutine:: PDM_writer_var_set

        .. ifconfig:: enable_fortran_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



      .. tab-item:: Python
        :sync: Python

        .. ifconfig:: enable_python_doc == 'ON'

          .. autofunction:: Pypdm.Pypdm.Writer.var_create
          .. autofunction:: Pypdm.Pypdm.Writer.name_map_add
          .. autofunction:: Pypdm.Pypdm.Writer.var_set

        .. ifconfig:: enable_python_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)

    .. doxygenenum:: PDM_writer_var_loc_t

  .. dropdown:: Global variables

    .. tab-set::
      :sync-group: language

      .. tab-item:: C
        :sync: C

        .. doxygenfunction:: PDM_writer_cst_global_var_create
        .. doxygenfunction:: PDM_writer_cst_global_var_set



      .. tab-item:: Fortran
        :sync: Fortran

        .. ifconfig:: enable_fortran_doc == 'ON'

          .. f:autosubroutine:: PDM_writer_cst_global_var_create
          .. f:autosubroutine:: PDM_writer_cst_global_var_set

        .. ifconfig:: enable_fortran_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)




.. dropdown:: Manage time steps

  All writes must be framed by the beginning and end of a time step.

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_writer_step_beg
      .. doxygenfunction:: PDM_writer_step_end



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_writer_step_beg
        .. f:autosubroutine:: PDM_writer_step_end

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.Writer.step_beg
        .. autofunction:: Pypdm.Pypdm.Writer.step_end

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Write

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_writer_geom_write
      .. doxygenfunction:: PDM_writer_var_write



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_writer_geom_write
        .. f:autosubroutine:: PDM_writer_var_write

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.Writer.geom_write
        .. autofunction:: Pypdm.Pypdm.Writer.var_write

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)



.. dropdown:: Finalization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_writer_var_free
      .. doxygenfunction:: PDM_writer_geom_free
      .. doxygenfunction:: PDM_writer_free



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_writer_var_free
        .. f:autosubroutine:: PDM_writer_geom_free
        .. f:autosubroutine:: PDM_writer_free

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      |python_gc|


Examples
""""""""