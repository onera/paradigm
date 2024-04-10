.. _writer:

Writer
======

C API
-----

Enumerators
~~~~~~~~~~~

.. doxygenenum:: PDM_writer_status_t

.. doxygenenum:: PDM_writer_topology_t

.. doxygenenum:: PDM_writer_fmt_fic_t

.. doxygenenum:: PDM_writer_var_dim_t

.. doxygenenum:: PDM_writer_var_loc_t

The Writer feature allows several types of elements for the input mesh, described by :

.. doxygenenum:: PDM_writer_elt_geom_t

Refer to the table below for the numbering convention of the standard elements in the list.

.. list-table:: Numbering convention used in ParaDiGM for standard elements
  :widths: 50 50

  * - .. figure:: ../../../../images/pdm_mesh_nodal_point.svg
        :alt: PDM_MESH_NODAL_POINT

        ``PDM_MESH_NODAL_POINT``

    - .. figure:: ../../../../images/pdm_mesh_nodal_bar2.svg
        :alt: PDM_MESH_NODAL_BAR2

        ``PDM_MESH_NODAL_BAR2``


  * - .. figure:: ../../../../images/pdm_mesh_nodal_tria3.svg
        :alt: PDM_MESH_NODAL_TRIA3

        ``PDM_MESH_NODAL_TRIA3``


    - .. figure:: ../../../../images/pdm_mesh_nodal_quad4.svg
        :alt: PDM_MESH_NODAL_QUAD4

        ``PDM_MESH_NODAL_QUAD4``


  * - .. figure:: ../../../../images/pdm_mesh_nodal_tetra4.svg
        :alt: PDM_MESH_NODAL_TETRA4

        ``PDM_MESH_NODAL_TETRA4``


    - .. figure:: ../../../../images/pdm_mesh_nodal_pyram5.svg
        :alt: PDM_MESH_NODAL_PYRAM5

        ``PDM_MESH_NODAL_PYRAM5``


  * - .. figure:: ../../../../images/pdm_mesh_nodal_prism6.svg
        :alt: PDM_MESH_NODAL_PRISM6

        ``PDM_MESH_NODAL_PRISM6``


    - .. figure:: ../../../../images/pdm_mesh_nodal_hexa8.svg
        :alt: PDM_MESH_NODAL_HEXA8

        ``PDM_MESH_NODAL_HEXA8``

Initialization
~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_writer_create

Step
~~~~

.. doxygenfunction:: PDM_writer_step_beg

.. doxygenfunction:: PDM_writer_is_open_step

.. doxygenfunction:: PDM_writer_step_end

Geometry
~~~~~~~~

.. doxygenfunction:: PDM_writer_geom_create

.. doxygenfunction:: PDM_writer_geom_create_from_mesh_nodal

.. doxygenfunction:: PDM_writer_geom_set_from_mesh_nodal

.. doxygenfunction:: PDM_writer_geom_coord_set

.. doxygenfunction:: PDM_writer_geom_coord_from_parent_set

.. doxygenfunction:: PDM_writer_geom_bloc_add

.. doxygenfunction:: PDM_writer_geom_bloc_std_set

.. doxygenfunction:: PDM_writer_geom_bloc_poly2d_set

.. doxygenfunction:: PDM_writer_geom_bloc_poly3d_set

.. doxygenfunction:: PDM_writer_geom_cell3d_cellface_add

.. doxygenfunction:: PDM_writer_geom_cell2d_cellface_add

.. doxygenfunction:: PDM_writer_geom_faces_facesom_add

.. doxygenfunction:: PDM_writer_geom_write

.. doxygenfunction:: PDM_writer_geom_data_reset

Variable & Co
~~~~~~~~~~~~~

.. doxygenfunction:: PDM_writer_cst_global_var_create

.. doxygenfunction:: PDM_writer_cst_global_var_set

.. doxygenfunction:: PDM_writer_var_create

.. doxygenfunction:: PDM_writer_var_write

.. doxygenfunction:: PDM_writer_var_set

.. doxygenfunction:: PDM_writer_fmt_add

.. doxygenfunction:: PDM_writer_name_map_add

Finalize
~~~~~~~~

.. doxygenfunction:: PDM_writer_fmt_free

.. doxygenfunction:: PDM_writer_var_data_free

.. doxygenfunction:: PDM_writer_var_free

.. doxygenfunction:: PDM_writer_geom_data_free

.. doxygenfunction:: PDM_writer_geom_free

.. doxygenfunction:: PDM_writer_free

Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  Initialization
  ~~~~~~~~~~~~~~

  .. f:subroutine:: PDM_writer_create(def)

    Create a structure for parallel writing of geometry and associated variables

    :param c_ptr            cs                 [out]: Pointer to Writer instance
    :param character        fmt                [in]:  Output format
    :param integer          fmt_fic            [in]:  Binary or ASCII
    :param integer          topologie          [in]:  Indicates whether the mesh is mobile
    :param integer          st_reprise         [in]:  Finalizes previous outputs before restart
    :param character        rep_sortie         [in]:  Output repository
    :param character        nom_sortie         [in]:  Output filename
    :param integer          f_comm             [in]:  MPI communicator
    :param integer          acces              [in]:  Access type
    :param real             prop_noeuds_actifs [in]:  Amount of active nodes:
                                                        *  -1 : all active
                                                        *   1 : one process per node
                                                        * 0 < val < 1 : one process per active node
    :param character        options            [in]:  Complementary options for the format structured as
                                                      ("name_1 = val_1 : ... : name_n = val_n")

  Step
  ~~~~

  .. f:subroutine:: PDM_writer_step_beg(def)

    Begin a time step

    :param c_ptr   cs            [in]: Pointer to Writer instance
    :param real    physical_time [in]: Time

  .. f:subroutine:: PDM_writer_step_end(def)

    Increment end

    :param c_ptr cs [in]: Pointer to Writer instance

  Geometry
  ~~~~~~~~

  .. f:subroutine:: PDM_writer_geom_create(def)

    Create a new geometry in the writer structure

    :param c_ptr     cs       [in]:  Pointer to Writer instance
    :param integer   id_geom  [out]: Identifier of the geometry within the writer instance
    :param character nom_geom [in]:  Name of the geometry
    :param integer   n_part   [in]:  Number of partitions

  .. f:subroutine:: PDM_writer_geom_coord_set(def)

    Define the coordinates of the current partition

    :param c_ptr                   cs       [in]: Pointer to Writer instance
    :param integer                 id_geom  [in]: Geometry identifier
    :param integer                 id_part  [in]: Partition identifier
    :param integer                 n_som    [in]: Number of vertices
    :param real(8)(:,:)            coords   [in]: Coordinates (shape = [3, ``n_som``])
    :param integer(pdm_g_num_s)(:) numabs   [in]: Vertex global numbering (size = ``n_som``)
    :param integer                 owner    [in]: Ownership

  .. f:subroutine:: PDM_writer_geom_coord_from_parent_set(def)

    Definition of the coordinates of the vertices
    in the current partition from a parent set

    :param c_ptr                   cs            [in]: Pointer to Writer instance
    :param integer                 id_geom       [in]: Geometry identifier
    :param integer                 id_part       [in]: Partition identifier
    :param integer                 n_som         [in]: Number of vertices
    :param integer                 n_som_parent  [in]: Number of parent vertices
    :param integer(pdm_g_num_s)(:) numabs        [in]: Vertex global numbering (size = ``n_som``)
    :param integer(pdm_l_num_s)(:) num_parent    [in]: Vertex parent local numbering (size = ``n_som``)
    :param real(8)(:,:)            coords_parent [in]: Coordinates of parent vertices (shape = [3, ``n_som_parent``])
    :param integer(pdm_g_num_s)(:) numabs_parent [in]: Vertex parent global numbering (size = ``n_som_parent``)
    :param integer                 owner         [in]: Ownership

  .. f:subroutine:: PDM_writer_geom_bloc_add(def)

    Add a block of elements of a given type

    :param c_ptr   cs       [in]:  Pointer to Writer instance
    :param integer id_geom  [in]:  Geometry identifier
    :param integer t_elt    [in]:  Element type
    :param integer owner    [in]:  Ownership
    :param integer id_block [out]: Block identifier

  .. f:subroutine:: PDM_writer_geom_bloc_std_set(def)

    Set in the given geometry a block of elements of a given type

    :param c_ptr                   cs       [in]: Pointer to Writer instance
    :param integer                 id_geom  [in]: Geometry identifier
    :param integer                 id_bloc  [in]: Block identifier
    :param integer                 id_part  [in]: Partition identifier
    :param integer                 n_elt    [in]: Number of elements
    :param integer(pdm_l_num_s)(:) connec   [in]: Element->Vertex connectivity
    :param integer(pdm_g_num_s)(:) numabs   [in]: Element global numbering

  .. f:subroutine:: PDM_writer_geom_bloc_poly2d_set(def)

    Add a block of polygons to the current partition

    :param c_ptr                   cs         [in]: Pointer to Writer instance
    :param integer                 id_geom    [in]: Geometry identifier
    :param integer                 id_bloc    [in]: Block identifier
    :param integer                 id_part    [in]: Partition identifier
    :param integer                 n_elt      [in]: Number of elements
    :param integer(pdm_l_num_s)(:) connec_idx [in]: Index of the Element->Vertex connectivity (size = ``n_elt`` + 1)
    :param integer(pdm_l_num_s)(:) connec     [in]: Element->Vertex connectivity (size = ``connec_idx(n_elt)``)
    :param integer(pdm_g_num_s)(:) numabs     [in]: Element global numbering (size = ``n_elt``)

  .. f:subroutine:: PDM_writer_geom_bloc_poly3d_set(def)

    Add a block of polyhedra to the current partition

    :param c_ptr                   cs          [in]: Pointer to Writer instance
    :param integer                 id_geom     [in]: Geometry identifier
    :param integer                 id_bloc     [in]: Block identifier
    :param integer                 id_part     [in]: Partition identifier
    :param integer                 n_elt       [in]: Number of elements
    :param integer                 n_face      [in]: Number of faces
    :param integer(pdm_l_num_s)(:) facsom_idx  [in]: Index of the Face->Vertex connectivity (size = ``n_face`` + 1)
    :param integer(pdm_l_num_s)(:) facsom      [in]: Face->Vertex connectivity (size = ``facsom_idx(n_face)``)
    :param integer(pdm_l_num_s)(:) cellfac_idx [in]: Index of the Cell->Face connectivity (size = ``n_elt`` + 1)
    :param integer(pdm_l_num_s)(:) cellfac     [in]: Cell->Face connectivity (size = ``cellfac_idx(n_elt)``)
    :param integer(pdm_g_num_s)(:) numabs      [in]: Cell global numbering (size = ``n_elt``)

  .. f:subroutine:: PDM_writer_geom_cell3d_cellface_add(def)

    Add 3D cells described in terms of faces
    This function determines element types and creates
    blocks grouping elements of the same type.
    It returns the indirection to the new arrangement order of the cells.

    :param c_ptr                   cs            [in]: Pointer to Writer instance
    :param integer                 id_geom       [in]: Geometry identifier
    :param integer                 id_part       [in]: Partition identifier
    :param integer                 n_cell        [in]: Number of 3D cells
    :param integer                 n_face        [in]: Number of faces
    :param integer(pdm_l_num_s)(:) face_som_idx  [in]: Index of the Face->Vertex connectivity (size = ``n_face`` + 1)
    :param integer(pdm_l_num_s)(:) face_som_nb   [in]: Number of vertices per face (optional)
    :param integer(pdm_l_num_s)(:) face_som      [in]: Face->Vertex connectivity (size = ``face_som_idx(n_face)``)
    :param integer(pdm_l_num_s)(:) cell_face_idx [in]: Index of the Cell->Face connectivity (size = ``n_cell`` + 1)
    :param integer(pdm_l_num_s)(:) cell_face_nb  [in]: Number of faces per cell (optional)
    :param integer(pdm_l_num_s)(:) cell_face     [in]: Cell->Face connectivity (size = ``cell_face_idx(n_cell)``)
    :param integer(pdm_g_num_s)(:) numabs        [in]: Cell global numbering (size = ``n_cell``)

  .. f:subroutine:: PDM_writer_geom_cell2d_cellface_add(def)

    Add of 2D cells described in terms of faces
    This function determines element types and creates
    blocks grouping elements of the same type.
    It returns the indirection to the new arrangement order of the cells.

    :param c_ptr                   cs            [in]: Pointer to Writer instance
    :param integer                 id_geom       [in]: Geometry identifier
    :param integer                 id_part       [in]: Partition identifier
    :param integer                 n_cell        [in]: Number of 2D cells
    :param integer                 n_face        [in]: Number of faces
    :param integer(pdm_l_num_s)(:) face_som_idx  [in]: Index of the Face->Vertex connectivity (unused)
    :param integer(pdm_l_num_s)(:) face_som_nb   [in]: Number of Vertices per Face (unused)
    :param integer(pdm_l_num_s)(:) face_som      [in]: Face->Vertex connectivity (size = 2 * ``n_face``)
    :param integer(pdm_l_num_s)(:) cell_face_idx [in]: Index of the Cell->Face connectivity (size = ``n_cell`` + 1)
    :param integer(pdm_l_num_s)(:) cell_face_nb  [in]: Number of Faces per Cell (optional)
    :param integer(pdm_l_num_s)(:) cell_face     [in]: Cell->Face connectivity (size = ``cell_face_idx(n_cell)``)
    :param integer(pdm_g_num_s)(:) numabs        [in]: Cell global numbering (size = ``n_cell``)

  .. f:subroutine:: PDM_writer_geom_faces_facesom_add(def)

    Add of faces described in nodal fashion
    This function determines element types and creates
    blocks grouping elements of the same type.
    It returns the indirection to the new arrangement order of the cells.

    :param c_ptr                   cs           [in]: Pointer to Writer instance
    :param integer                 id_geom      [in]: Geometry identifier
    :param integer                 id_part      [in]: Partition identifier
    :param integer                 n_face       [in]: Number of faces
    :param integer(pdm_l_num_s)(:) face_som_idx [in]: Index of the Face->Vertex connectivity (size = ``n_face`` + 1)
    :param integer(pdm_l_num_s)(:) face_som_nb  [in]: Number of Vertices per Face (optional)
    :param integer(pdm_l_num_s)(:) face_som     [in]: Face->Vertex connectivity (size = ``face_som_idx(n_face)``)
    :param integer(pdm_g_num_s)(:) numabs       [in]: Face global numbering (size = ``n_face``)

  .. f:subroutine:: PDM_writer_geom_write(def)

    Write of current mesh

    :param c_ptr   cs      [in]: Pointer to Writer instance
    :param integer id_geom [in]: Geometry identifier

  .. f:subroutine:: PDM_writer_geom_data_reset(def)

    Reset of data describing the current mesh

    :param c_ptr   cs      [in]: Pointer to Writer instance
    :param integer id_geom [in]: Geometry identifier

  Variable & Co
  ~~~~~~~~~~~~~

  .. f:subroutine:: PDM_writer_cst_global_var_create(def)

    Create a global constant variable

    :param c_ptr            cs      [in]:  Pointer to Writer instance
    :param integer          id_var  [out]: Variable identifier
    :param character        nom_var [in]:  Variable name
    :param real             val_var [in]:  Variable value

  .. f:subroutine:: PDM_writer_cst_global_var_set(def)

    Create a global constant variable

    :param c_ptr            cs      [in]: Pointer to Writer instance
    :param integer          id_var  [in]: Variable identifier
    :param real             val_var [in]: Variable value

  .. f:subroutine:: PDM_writer_var_create(def)

    Create a variable

    :param c_ptr     cs           [in]:  Pointer to Writer instance
    :param integer   id_var       [out]: Variable identifier
    :param integer   st_dep_temps [in]:  Indicates whether the variable is time dependent
    :param integer   dim          [in]:  Variable's dimension
    :param integer   loc          [in]:  Variable's location
    :param character nom_var      [in]:  Name of the variable

  .. f:subroutine:: PDM_writer_var_write(def)

    Write variable values

    :param c_ptr     cs     [in]: Pointer to Writer instance
    :param integer   id_var [in]: Variable identifier

  .. f:subroutine:: PDM_writer_var_set(def)

    Update variable values

    .. warning:: the values defined for the elements must be defined in the order in which the blocks are defined!

    :param c_ptr            cs      [in]: Pointer to Writer instance
    :param integer          id_var  [in]: Variable identifier
    :param integer          id_geom [in]: Geometry identifier
    :param integer          id_part [in]: Partition identifier
    :param real(8)(:)       val     [in]: Variable values

  .. f:subroutine:: PDM_writer_fmt_add(def)

    Define a new format writer

    .. warning:: has not been tested, not sure about procedure pointer interoperability

    :param character   name            [in]: Name
    :param procedure() create_fct      [in]: Customize \ref PDM_writer_create function for the new format  (or NULL)
    :param procedure() free_fct        [in]: Customize \ref PDM_writer_free function for the new format (or NULL)
    :param procedure() beg_step_fct    [in]: Customize \ref PDM_writer_step_beg function for the new format (or NULL)
    :param procedure() end_step_fct    [in]: Customize \ref PDM_writer_step_end function for the new format (or NULL)
    :param procedure() geom_create_fct [in]: Customize \ref PDM_writer_geom_create function for the new format (or NULL)
    :param procedure() geom_write_fct  [in]: Customize \ref PDM_writer_geom_write function for the new format
    :param procedure() geom_free_fct   [in]: Customize \ref PDM_writer_geom_free function for the new format (or NULL)
    :param procedure() var_create_fct  [in]: Customize \ref PDM_writer_var_create function for the new format (or NULL)
    :param procedure() var_write_fct   [in]: Customize \ref PDM_writer_var_write function for the new format
    :param procedure() var_free_fct    [in]: Customize \ref PDM_writer_var_free function for the new format (or NULL)

  .. f:subroutine:: PDM_writer_name_map_add(def)

    Variable name mapping

    :param c_ptr     cs          [in]: Pointer to Writer instance
    :param character public_name [in]: Public variable name
    :param character pivate_name [in]: Private variable name

  Finalize
  ~~~~~~~~

  .. f:subroutine:: PDM_writer_fmt_free(def)

    Free format

  .. f:subroutine:: PDM_writer_var_data_free(def)

    Free variable data arrays

    :param c_ptr   cs     [in]: Pointer to Writer instance
    :param integer id_var [in]: Variable identifier

  .. f:subroutine:: PDM_writer_var_free(def)

    Free variable

    :param c_ptr   cs     [in]: Pointer to Writer instance
    :param integer id_var [in]: Variable identifier

  .. f:subroutine:: PDM_writer_geom_data_free(def)

    Free data describing the current mesh
    Indirections on absolute numbering are retained

    :param c_ptr   cs      [in]: Pointer to Writer instance
    :param integer id_geom [in]: Geometry identifier

  .. f:subroutine:: PDM_writer_geom_free(def)

    Free data describing the current mesh

    :param c_ptr   cs      [in]: Pointer to Writer instance
    :param integer id_geom [in]: Geometry identifier

  .. f:subroutine:: PDM_writer_free(def)

    Free a writer structure

    :param c_ptr   cs      [in]: Pointer to Writer instance

Python API
----------

.. ifconfig:: enable_python_doc == 'ON'

  .. py:class:: Writer

    Python object to perform mesh and associated data write.
    Once initialized, all the following
    methods apply to a :class:`Writer` instance.

    .. rubric:: Initialization

    .. autofunction:: Pypdm.Pypdm.Writer.__init__

    .. rubric:: Methods summary

    .. autosummary::
      :nosignatures:

      ~Pypdm.Pypdm.Writer.geom_create
      ~Pypdm.Pypdm.Writer.geom_cell2d_cellface_add
      ~Pypdm.Pypdm.Writer.geom_cell3d_cellface_add
      ~Pypdm.Pypdm.Writer.geom_coord_set
      ~Pypdm.Pypdm.Writer.geom_faces_facevtx_add
      ~Pypdm.Pypdm.Writer.geom_block_add
      ~Pypdm.Pypdm.Writer.geom_block_std_set
      ~Pypdm.Pypdm.Writer.geom_write
      ~Pypdm.Pypdm.Writer.geom_data_free
      ~Pypdm.Pypdm.Writer.geom_free
      ~Pypdm.Pypdm.Writer.var_create
      ~Pypdm.Pypdm.Writer.name_map_add
      ~Pypdm.Pypdm.Writer.var_write
      ~Pypdm.Pypdm.Writer.var_set
      ~Pypdm.Pypdm.Writer.var_data_free
      ~Pypdm.Pypdm.Writer.var_free
      ~Pypdm.Pypdm.Writer.step_beg
      ~Pypdm.Pypdm.Writer.step_end

    .. rubric:: Step

    .. automethod:: Pypdm.Pypdm.Writer.step_beg

    .. automethod:: Pypdm.Pypdm.Writer.step_end

    .. rubric:: Geometry

    .. automethod:: Pypdm.Pypdm.Writer.geom_create

    .. automethod:: Pypdm.Pypdm.Writer.geom_cell2d_cellface_add

    .. automethod:: Pypdm.Pypdm.Writer.geom_cell3d_cellface_add

    .. automethod:: Pypdm.Pypdm.Writer.geom_coord_set

    .. automethod:: Pypdm.Pypdm.Writer.geom_faces_facevtx_add

    .. automethod:: Pypdm.Pypdm.Writer.geom_block_add

    .. automethod:: Pypdm.Pypdm.Writer.geom_block_std_set

    .. automethod:: Pypdm.Pypdm.Writer.geom_write

    .. rubric:: Variable & Co

    .. automethod:: Pypdm.Pypdm.Writer.var_write

    .. automethod:: Pypdm.Pypdm.Writer.var_set

    .. automethod:: Pypdm.Pypdm.Writer.var_free

    .. rubric:: Finalize

    .. automethod:: Pypdm.Pypdm.Writer.geom_data_free

    .. automethod:: Pypdm.Pypdm.Writer.geom_free

    .. automethod:: Pypdm.Pypdm.Writer.var_data_free
