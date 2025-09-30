module pdm_writer

  use iso_c_binding
  use pdm

  implicit none

  !
  ! Statut

  integer, parameter :: PDM_WRITER_OFF = 0
  integer, parameter :: PDM_WRITER_ON  = 1

  !
  ! Type de topologie

  integer, parameter :: PDM_WRITER_TOPO_CST        = 0
  integer, parameter :: PDM_WRITER_TOPO_DEFORMABLE = 1
  integer, parameter :: PDM_WRITER_TOPO_VARIABLE   = 2

  !
  ! Type d'elements géometriques

  integer, parameter :: PDM_WRITER_POINT    = 0
  integer, parameter :: PDM_WRITER_BAR2     = 1
  integer, parameter :: PDM_WRITER_TRIA3    = 2
  integer, parameter :: PDM_WRITER_QUAD4    = 3
  integer, parameter :: PDM_WRITER_POLY_2D  = 4
  integer, parameter :: PDM_WRITER_TETRA4   = 5
  integer, parameter :: PDM_WRITER_PYRAMID5 = 6
  integer, parameter :: PDM_WRITER_PRISM6   = 7
  integer, parameter :: PDM_WRITER_HEXA8    = 8
  integer, parameter :: PDM_WRITER_POLY_3D  = 9

  !
  ! Format de sortie

  integer, parameter ::  PDM_WRITER_FMT_ENSIGHT = 0

  !
  ! Format du fichier

  integer, parameter ::  PDM_WRITER_FMT_BIN   = 0
  integer, parameter ::  PDM_WRITER_FMT_ASCII = 1

  !
  ! Dimension géométrique de la sortie

  integer, parameter ::  PDM_WRITER_DIM_2 = 0
  integer, parameter ::  PDM_WRITER_DIM_3 = 1

  !
  ! Dim des variables

  integer, parameter ::  PDM_WRITER_VAR_CSTE         = 0
  integer, parameter ::  PDM_WRITER_VAR_SCALAIRE     = 1
  integer, parameter ::  PDM_WRITER_VAR_VECTOR       = 3
  integer, parameter ::  PDM_WRITER_VAR_TENSEUR_SYM  = 6
  integer, parameter ::  PDM_WRITER_VAR_TENSEUR_ASYM = 9

  !
  ! Localisation des variables

  integer, parameter ::  PDM_WRITER_VAR_VERTICES   = 0
  integer, parameter ::  PDM_WRITER_VAR_ELEMENTS   = 1
  integer, parameter ::  PDM_WRITER_VAR_PARTICULES = 2



  interface





  !>
  !! \brief Free format
  !!
  !!

  subroutine PDM_writer_fmt_free () &
  bind (c, name='PDM_writer_fmt_free')
  end subroutine PDM_writer_fmt_free


  end interface


  contains


  subroutine PDM_writer_create(cs,                 &
                               fmt,                &
                               fmt_fic,            &
                               topologie,          &
                               st_reprise,         &
                               rep_sortie,         &
                               nom_sortie,         &
                               f_comm,             &
                               acces,              &
                               prop_noeuds_actifs, &
                               options)
    ! Create a new PDM_writer_t instance
    implicit none

    type(c_ptr),      intent(out) :: cs                 ! Pointer to PDM_writer_t instance
    character(len=*), intent(in)  :: fmt                ! Format
    integer,          intent(in)  :: fmt_fic            ! Binary or ASCII
    integer,          intent(in)  :: topologie          ! Indicates whether the mesh is mobile
    integer,          intent(in)  :: st_reprise         ! Finalizes previous outputs before restart
    character(len=*), intent(in)  :: rep_sortie         ! Output repository
    character(len=*), intent(in)  :: nom_sortie         ! Output filename
    integer,          intent(in)  :: f_comm             ! MPI communicator
    integer,          intent(in)  :: acces              ! Access type
    double precision, intent(in)  :: prop_noeuds_actifs ! Amount of active nodes (-1: all active, 1: one process per node, 0 < val < 1: some processes per active node)
    character(len=*), intent(in)  :: options            ! Complementary options for the format structured as ("name_1 = val_1 : ... : name_n = val_n")

    integer(c_int)                :: c_fmt_fic
    integer(c_int)                :: c_topologie
    integer(c_int)                :: c_st_reprise
    integer(c_int)                :: c_comm
    integer(c_int)                :: c_acces
    real(c_double)                :: c_prop_noeuds_actifs

    interface
      function PDM_writer_create_c(fmt,                &
                                   fmt_fic,            &
                                   topologie,          &
                                   st_reprise,         &
                                   rep_sortie,         &
                                   nom_sortie,         &
                                   pdm_mpi_comm,       &
                                   acces,              &
                                   prop_noeuds_actifs, &
                                   options)            &
      result (cs)                                      &
      bind (c, name='PDM_writer_create')
        use iso_c_binding
        implicit none

        type(c_ptr)           :: cs
        character(c_char)     :: fmt
        integer(c_int), value :: fmt_fic
        integer(c_int), value :: topologie
        integer(c_int), value :: st_reprise
        character(c_char)     :: rep_sortie
        character(c_char)     :: nom_sortie
        integer(c_int), value :: pdm_mpi_comm
        integer(c_int), value :: acces
        real(c_double), value :: prop_noeuds_actifs
        character(c_char)     :: options

      end function PDM_writer_create_c
    end interface

    c_comm               = PDM_MPI_Comm_f2c(f_comm)
    c_fmt_fic            = fmt_fic
    c_topologie          = topologie
    c_st_reprise         = st_reprise
    c_acces              = acces
    c_prop_noeuds_actifs = prop_noeuds_actifs

    cs = PDM_writer_create_c(trim(fmt)//C_NULL_CHAR,        &
                             c_fmt_fic,                     &
                             c_topologie,                   &
                             c_st_reprise,                  &
                             trim(rep_sortie)//C_NULL_CHAR, &
                             trim(nom_sortie)//C_NULL_CHAR, &
                             c_comm,                        &
                             c_acces,                       &
                             c_prop_noeuds_actifs,          &
                             trim(options)//C_NULL_CHAR)

  end subroutine PDM_writer_create



  subroutine PDM_writer_step_beg(cs,            &
                                 physical_time)
    ! Begin a time step
    implicit none

    type(c_ptr),      intent(in) :: cs            ! Pointer to PDM_writer_t instance
    double precision, intent(in) :: physical_time ! Time


    interface
      subroutine PDM_writer_step_beg_c(cs,            &
                                       physical_time) &
      bind (c, name="PDM_writer_step_beg")

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        real(c_double), value :: physical_time

      end subroutine PDM_writer_step_beg_c
    end interface

    call PDM_writer_step_beg_c(cs, &
                               physical_time)

  end subroutine PDM_writer_step_beg



  subroutine PDM_writer_step_end(cs)

    ! End a time step
    implicit none

    type(c_ptr), intent(in) :: cs ! Pointer to PDM_writer_t instance

    interface
      subroutine PDM_writer_step_end_c(cs) &
        bind (c, name="PDM_writer_step_end")

        use iso_c_binding
        implicit none

        type(c_ptr), value :: cs

      end subroutine PDM_writer_step_end_c
    end interface

    call PDM_writer_step_end_c(cs)

  end subroutine PDM_writer_step_end



  subroutine PDM_writer_geom_create(cs,       &
                                    id_geom,  &
                                    nom_geom, &
                                    n_part)
    ! Create a new geometry in the writer structure
    implicit none

    type(c_ptr),      intent(in)  :: cs       ! Pointer to PDM_writer_t instance
    integer,          intent(out) :: id_geom  ! Geometry identifier
    character(len=*), intent(in)  :: nom_geom ! Name of the geometry
    integer,          intent(in)  :: n_part   ! Number of partitions

    interface
      function PDM_writer_geom_create_c(cs,       &
                                        nom_geom, &
                                        n_part)   &
      result (id_geom)                            &
      bind (c, name='PDM_writer_geom_create')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int)        :: id_geom
        character(c_char)     :: nom_geom(*)
        integer(c_int), value :: n_part

      end function PDM_writer_geom_create_c
    end interface

    id_geom = PDM_writer_geom_create_c(cs,                          &
                                       trim(nom_geom)//C_NULL_CHAR, &
                                       n_part)

  end subroutine PDM_writer_geom_create



  subroutine PDM_writer_geom_create_from_mesh_nodal(cs,       &
                                                    id_geom,  &
                                                    nom_geom, &
                                                    mesh)
    ! Create a geometry from a nodal mesh structure
    use iso_c_binding
    implicit none

    type(c_ptr),      intent(in)  :: cs        ! Pointer to PDM_writer_t instance
    integer,          intent(out) :: id_geom   ! Geometry identifier
    character(len=*), intent(in)  :: nom_geom  ! Geometry name
    type(c_ptr),      intent(in)  :: mesh      ! Pointer to PDM_part_mesh_nodal_t instance

    interface
      function PDM_writer_geom_create_from_mesh_nodal_c(cs,       &
                                                        nom_geom, &
                                                        mesh)     &
        result (id_geom)                                          &
        bind (c, name='PDM_writer_geom_create_from_mesh_nodal')
          use iso_c_binding
          implicit none

          type(c_ptr),    value :: cs
          integer(c_int)        :: id_geom
          character(c_char)     :: nom_geom(*)
          type(c_ptr),    value :: mesh

      end function PDM_writer_geom_create_from_mesh_nodal_c
    end interface

    id_geom = PDM_writer_geom_create_from_mesh_nodal_c(cs, &
                                                       trim(nom_geom)//C_NULL_CHAR, &
                                                       mesh)

  end subroutine PDM_writer_geom_create_from_mesh_nodal



  subroutine PDM_writer_geom_coord_set(cs,      &
                                       id_geom, &
                                       id_part, &
                                       n_som,   &
                                       coords,  &
                                       numabs,  &
                                       owner)
    ! Define the coordinates of the current partition
    implicit none

    type(c_ptr),                   intent(in) :: cs          ! Pointer to PDM_writer_t instance
    integer,                       intent(in) :: id_geom     ! Geometry identifier
    integer,                       intent(in) :: id_part     ! Partition identifier
    integer,                       intent(in) :: n_som       ! Number of vertices
    real(8),              pointer, intent(in) :: coords(:,:) ! Coordinates (shape = [3, n_som])
    integer(pdm_g_num_s), pointer, intent(in) :: numabs(:)   ! Vertex global IDs
    integer,                       intent(in) :: owner       ! Ownership

    type(c_ptr)                               :: c_coords
    type(c_ptr)                               :: c_numabs

    interface
      subroutine PDM_writer_geom_coord_set_c (cs,      &
                                              id_geom, &
                                              id_part, &
                                              n_som,   &
                                              coords,  &
                                              numabs,  &
                                              owner)   &
      bind (c, name='PDM_writer_geom_coord_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom
        integer(c_int), value :: id_part
        integer(c_int), value :: n_som
        integer(c_int), value :: owner
        type(c_ptr),    value :: coords
        type(c_ptr),    value :: numabs

      end subroutine PDM_writer_geom_coord_set_c
    end interface

    c_coords = C_NULL_PTR
    if (associated(coords)) then
      c_coords = c_loc(coords)
    endif

    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs = c_loc(numabs)
    endif


    call PDM_writer_geom_coord_set_c(cs,       &
                                     id_geom,  &
                                     id_part,  &
                                     n_som,    &
                                     c_coords, &
                                     c_numabs, &
                                     owner)

  end subroutine PDM_writer_geom_coord_set



  subroutine PDM_writer_geom_coord_from_parent_set(cs,            &
                                                   id_geom,       &
                                                   id_part,       &
                                                   n_som,         &
                                                   n_som_parent,  &
                                                   numabs,        &
                                                   num_parent,    &
                                                   coords_parent, &
                                                   numabs_parent, &
                                                   owner)

    ! Definition of the coordinates of the vertices in the current partition from a parent set
    implicit none

    type(c_ptr),                   intent(in) :: cs                 ! Pointer to PDM_writer_t instance
    integer,                       intent(in) :: id_geom            ! Geometry identifier
    integer,                       intent(in) :: id_part            ! Partition identifier
    integer,                       intent(in) :: n_som              ! Number of vertices
    integer,                       intent(in) :: n_som_parent       ! Number of parent vertices
    integer(pdm_g_num_s), pointer, intent(in) :: numabs(:)          ! Vertex global IDs (size = n_som)
    integer(pdm_l_num_s), pointer, intent(in) :: num_parent(:)      ! Vertex parent local IDs (size = n_som)
    real(8),              pointer, intent(in) :: coords_parent(:,:) ! Coordinates of parent vertices (shape = [3, n_som_parent])
    integer(pdm_g_num_s), pointer, intent(in) :: numabs_parent(:)   ! Vertex parent global IDs (size = n_som_parent)
    integer,                       intent(in) :: owner              ! Ownership

    type(c_ptr)                               :: c_numabs
    type(c_ptr)                               :: c_num_parent
    type(c_ptr)                               :: c_coords_parent
    type(c_ptr)                               :: c_numabs_parent

    interface
      subroutine PDM_writer_geom_coord_from_parent_set_c(cs,            &
                                                         id_geom,       &
                                                         id_part,       &
                                                         n_som,         &
                                                         n_som_parent,  &
                                                         numabs,        &
                                                         num_parent,    &
                                                         coords_parent, &
                                                         numabs_parent, &
                                                         owner) &
      bind (c, name='PDM_writer_geom_coord_from_parent_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom
        integer(c_int), value :: id_part
        integer(c_int), value :: n_som
        integer(c_int), value :: n_som_parent
        integer(c_int), value :: owner
        type(c_ptr),    value :: numabs
        type(c_ptr),    value :: num_parent
        type(c_ptr),    value :: coords_parent
        type(c_ptr),    value :: numabs_parent

      end subroutine PDM_writer_geom_coord_from_parent_set_c
    end interface


    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs        = c_loc(numabs)
    endif

    c_num_parent = C_NULL_PTR
    if (associated(num_parent)) then
      c_num_parent    = c_loc(num_parent)
    endif

    c_coords_parent = C_NULL_PTR
    if (associated(coords_parent)) then
      c_coords_parent = c_loc(coords_parent)
    endif

    c_numabs_parent = C_NULL_PTR
    if (associated(numabs_parent)) then
      c_numabs_parent = c_loc(numabs_parent)
    endif


    call PDM_writer_geom_coord_from_parent_set_c(cs,              &
                                                 id_geom,         &
                                                 id_part,         &
                                                 n_som,           &
                                                 n_som_parent,    &
                                                 c_numabs,        &
                                                 c_num_parent,    &
                                                 c_coords_parent, &
                                                 c_numabs_parent, &
                                                 owner)

  end subroutine PDM_writer_geom_coord_from_parent_set



  subroutine PDM_writer_geom_bloc_add(cs,      &
                                      id_geom, &
                                      t_elt,   &
                                      owner,   &
                                      id_bloc)
    ! Add a section of elements of a given type
    implicit none

    type(c_ptr), intent(in)  :: cs      ! Pointer to PDM_writer_t instance
    integer,     intent(in)  :: id_geom ! Geometry identifier
    integer,     intent(in)  :: t_elt   ! Element type
    integer,     intent(in)  :: owner   ! Section ownership
    integer,     intent(out) :: id_bloc ! Section identifier

    interface
      function PDM_writer_geom_bloc_add_c(cs,      &
                                          id_geom, &
                                          t_elt,   &
                                          owner)   &
      result (id_bloc)                             &
      bind (c, name='PDM_writer_geom_bloc_add')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom
        integer(c_int), value :: t_elt
        integer(c_int), value :: owner
        integer(c_int)        :: id_bloc

      end function PDM_writer_geom_bloc_add_c
    end interface

    id_bloc = PDM_writer_geom_bloc_add_c(cs,      &
                                         id_geom, &
                                         t_elt,   &
                                         owner)

  end subroutine PDM_writer_geom_bloc_add



  subroutine PDM_writer_geom_bloc_std_set(cs,      &
                                          id_geom, &
                                          id_bloc, &
                                          id_part, &
                                          n_elt,   &
                                          connec,  &
                                          numabs)
    ! Set in the given geometry a section of elements of a given type
    implicit none

    type(c_ptr),                   intent(in) :: cs        ! Pointer to PDM_writer_t instance
    integer,                       intent(in) :: id_geom   ! Geometry identifier
    integer,                       intent(in) :: id_bloc   ! Section identifier
    integer,                       intent(in) :: id_part   ! Partition identifier
    integer,                       intent(in) :: n_elt     ! Number of elements
    integer(pdm_l_num_s), pointer, intent(in) :: connec(:) ! Element->Vertex connectivity
    integer(pdm_g_num_s), pointer, intent(in) :: numabs(:) ! Element global IDs

    type(c_ptr)                               :: c_connec
    type(c_ptr)                               :: c_numabs

    interface
      subroutine PDM_writer_geom_bloc_std_set_c(cs,      &
                                                id_geom, &
                                                id_bloc, &
                                                id_part, &
                                                n_elt,   &
                                                connec,  &
                                                numabs)  &
      bind (c, name='PDM_writer_geom_bloc_std_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom
        integer(c_int), value :: id_bloc
        integer(c_int), value :: id_part
        integer(c_int), value :: n_elt
        type(c_ptr),    value :: connec
        type(c_ptr),    value :: numabs

      end subroutine PDM_writer_geom_bloc_std_set_c
    end interface

    c_connec = C_NULL_PTR
    if (associated(connec)) then
      c_connec = c_loc(connec)
    endif

    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs = c_loc(numabs)
    endif


    call PDM_writer_geom_bloc_std_set_c(cs,       &
                                        id_geom,  &
                                        id_bloc,  &
                                        id_part,  &
                                        n_elt,    &
                                        c_connec, &
                                        c_numabs)

  end subroutine PDM_writer_geom_bloc_std_set



  subroutine PDM_writer_geom_bloc_poly2d_set(cs,         &
                                             id_geom,    &
                                             id_bloc,    &
                                             id_part,    &
                                             n_elt,      &
                                             connec_idx, &
                                             connec,     &
                                             numabs)
    ! Add a section of polygons to the current partition
    implicit none

    type(c_ptr),                   intent(in) :: cs            ! Pointer to PDM_writer_t instance
    integer,                       intent(in) :: id_geom       ! Geometry identifier
    integer,                       intent(in) :: id_bloc       ! Section identifier
    integer,                       intent(in) :: id_part       ! Partition identifier
    integer,                       intent(in) :: n_elt         ! Number of elements
    integer(pdm_l_num_s), pointer, intent(in) :: connec_idx(:) ! Index of the Element->Vertex connectivity (size = n_elt+1)
    integer(pdm_l_num_s), pointer, intent(in) :: connec(:)     ! Element->Vertex connectivity (size = connec_idx(n_elt+1))
    integer(pdm_g_num_s), pointer, intent(in) :: numabs(:)     ! Element global IDs (size = n_elt)

    type(c_ptr)                               :: c_connec_idx
    type(c_ptr)                               :: c_connec
    type(c_ptr)                               :: c_numabs

    interface
      subroutine PDM_writer_geom_bloc_poly2d_set_c(cs,         &
                                                   id_geom,    &
                                                   id_bloc,    &
                                                   id_part,    &
                                                   n_elt,      &
                                                   connec_idx, &
                                                   connec,     &
                                                   numabs)     &
      bind (c, name='PDM_writer_geom_bloc_poly2d_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom
        integer(c_int), value :: id_bloc
        integer(c_int), value :: id_part
        integer(c_int), value :: n_elt
        type(c_ptr),    value :: connec_idx
        type(c_ptr),    value :: connec
        type(c_ptr),    value :: numabs

      end subroutine PDM_writer_geom_bloc_poly2d_set_c
    end interface

    c_connec_idx = C_NULL_PTR
    if (associated(connec_idx)) then
      c_connec_idx = c_loc(connec_idx)
    endif

    c_connec = C_NULL_PTR
    if (associated(connec)) then
      c_connec = c_loc(connec)
    endif

    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs = c_loc(numabs)
    endif


    call PDM_writer_geom_bloc_poly2d_set_c(cs,           &
                                           id_geom,      &
                                           id_bloc,      &
                                           id_part,      &
                                           n_elt,        &
                                           c_connec_idx, &
                                           c_connec,     &
                                           c_numabs)

  end subroutine PDM_writer_geom_bloc_poly2d_set



  subroutine PDM_writer_geom_bloc_poly3d_set(cs,          &
                                             id_geom,     &
                                             id_bloc,     &
                                             id_part,     &
                                             n_elt,       &
                                             n_face,      &
                                             facsom_idx,  &
                                             facsom,      &
                                             cellfac_idx, &
                                             cellfac,     &
                                             numabs)
    ! Add a section of polyhedra to the current partition
    implicit none

    type(c_ptr),                   intent(in) :: cs             ! Pointer to PDM_writer_t instance
    integer,                       intent(in) :: id_geom        ! Geometry identifier
    integer,                       intent(in) :: id_bloc        ! Section identifier
    integer,                       intent(in) :: id_part        ! Partition identifier
    integer,                       intent(in) :: n_elt          ! Number of elements
    integer,                       intent(in) :: n_face         ! Number of faces
    integer(pdm_l_num_s), pointer, intent(in) :: facsom_idx(:)  ! Index of the Face->Vertex connectivity (size = n_face + 1)
    integer(pdm_l_num_s), pointer, intent(in) :: facsom(:)      ! Face->Vertex connectivity (size = facsom_idx(n_face+1))
    integer(pdm_l_num_s), pointer, intent(in) :: cellfac_idx(:) ! Index of the Cell->Face connectivity (size = n_elt+1)
    integer(pdm_l_num_s), pointer, intent(in) :: cellfac(:)     ! Cell->Face connectivity (size = cellfac_idx(n_elt+1))
    integer(pdm_g_num_s), pointer, intent(in) :: numabs(:)      ! Cell global IDs (size = n_elt)

    type(c_ptr)                               :: c_facsom_idx
    type(c_ptr)                               :: c_facsom
    type(c_ptr)                               :: c_cellfac_idx
    type(c_ptr)                               :: c_cellfac
    type(c_ptr)                               :: c_numabs

    interface
      subroutine PDM_writer_geom_bloc_poly3d_set_c(cs,          &
                                                   id_geom,     &
                                                   id_bloc,     &
                                                   id_part,     &
                                                   n_elt,       &
                                                   n_face,      &
                                                   facsom_idx,  &
                                                   facsom,      &
                                                   cellfac_idx, &
                                                   cellfac,     &
                                                   numabs)      &
      bind (c, name='PDM_writer_geom_bloc_poly3d_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom
        integer(c_int), value :: id_bloc
        integer(c_int), value :: id_part
        integer(c_int), value :: n_elt
        integer(c_int), value :: n_face
        type(c_ptr),    value :: facsom_idx
        type(c_ptr),    value :: facsom
        type(c_ptr),    value :: cellfac_idx
        type(c_ptr),    value :: cellfac
        type(c_ptr),    value :: numabs

      end subroutine PDM_writer_geom_bloc_poly3d_set_c
    end interface

    c_facsom_idx = C_NULL_PTR
    if (associated(facsom_idx)) then
      c_facsom_idx  = c_loc(facsom_idx)
    endif

    c_facsom = C_NULL_PTR
    if (associated(facsom)) then
      c_facsom      = c_loc(facsom)
    endif

    c_cellfac_idx = C_NULL_PTR
    if (associated(cellfac_idx)) then
      c_cellfac_idx = c_loc(cellfac_idx)
    endif

    c_cellfac = C_NULL_PTR
    if (associated(cellfac)) then
      c_cellfac     = c_loc(cellfac)
    endif

    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs      = c_loc(numabs)
    endif


    call PDM_writer_geom_bloc_poly3d_set_c(cs,            &
                                           id_geom,       &
                                           id_bloc,       &
                                           id_part,       &
                                           n_elt,         &
                                           n_face,        &
                                           c_facsom_idx,  &
                                           c_facsom,      &
                                           c_cellfac_idx, &
                                           c_cellfac,     &
                                           c_numabs)

  end subroutine PDM_writer_geom_bloc_poly3d_set



  subroutine PDM_writer_geom_cell3d_cellface_add(cs,            &
                                                 id_geom,       &
                                                 id_part,       &
                                                 n_cell,        &
                                                 n_face,        &
                                                 face_som_idx,  &
                                                 face_som_nb,   &
                                                 face_som,      &
                                                 cell_face_idx, &
                                                 cell_face_nb,  &
                                                 cell_face,     &
                                                 numabs)
    ! Add 3D cells described in terms of faces
    implicit none

    type(c_ptr),                   intent(in) :: cs               ! Pointer to PDM_writer_t instance
    integer,                       intent(in) :: id_geom          ! Geometry identifier
    integer,                       intent(in) :: id_part          ! Partition identifier
    integer,                       intent(in) :: n_cell           ! Number of 3D cells
    integer,                       intent(in) :: n_face           ! Number of faces
    integer(pdm_l_num_s), pointer, intent(in) :: face_som_idx(:)  ! Index of the Face->Vertex connectivity (size = n_face + 1)
    integer(pdm_l_num_s), pointer, intent(in) :: face_som_nb(:)   ! Number of vertices per face (optional, can be set to *null()*)
    integer(pdm_l_num_s), pointer, intent(in) :: face_som(:)      ! Face->Vertex connectivity (size = face_som_idx(n_face+1))
    integer(pdm_l_num_s), pointer, intent(in) :: cell_face_idx(:) ! Index of the Cell->Face connectivity (size = n_cell + 1)
    integer(pdm_l_num_s), pointer, intent(in) :: cell_face_nb(:)  ! Number of faces per cell (optional, can be set to *null()*)
    integer(pdm_l_num_s), pointer, intent(in) :: cell_face(:)     ! Cell->Face connectivity (size = cell_face_idx(n_cell+1))
    integer(pdm_g_num_s), pointer, intent(in) :: numabs(:)        ! Cell global IDs (size = n_cell)

    type(c_ptr)                               :: c_face_som_idx
    type(c_ptr)                               :: c_face_som_nb
    type(c_ptr)                               :: c_face_som
    type(c_ptr)                               :: c_cell_face_idx
    type(c_ptr)                               :: c_cell_face_nb
    type(c_ptr)                               :: c_cell_face
    type(c_ptr)                               :: c_numabs

    interface
      subroutine PDM_writer_geom_cell3d_cellface_add_c(cs,            &
                                                       id_geom,       &
                                                       id_part,       &
                                                       n_cell,        &
                                                       n_face,        &
                                                       face_som_idx,  &
                                                       face_som_nb,   &
                                                       face_som,      &
                                                       cell_face_idx, &
                                                       cell_face_nb,  &
                                                       cell_face,     &
                                                       numabs)        &
      bind (c, name='PDM_writer_geom_cell3d_cellface_add')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom
        integer(c_int), value :: id_part
        integer(c_int), value :: n_cell
        integer(c_int), value :: n_face
        type(c_ptr),    value :: face_som_idx
        type(c_ptr),    value :: face_som_nb
        type(c_ptr),    value :: face_som
        type(c_ptr),    value :: cell_face_idx
        type(c_ptr),    value :: cell_face_nb
        type(c_ptr),    value :: cell_face
        type(c_ptr),    value :: numabs

      end subroutine PDM_writer_geom_cell3d_cellface_add_c
    end interface

    c_face_som_idx = C_NULL_PTR
    if (associated(face_som_idx)) then
        c_face_som_idx  = c_loc(face_som_idx)
    endif

    c_face_som = C_NULL_PTR
    if (associated(face_som)) then
        c_face_som = c_loc(face_som)
    endif

    c_cell_face_idx = C_NULL_PTR
    if (associated(cell_face_idx)) then
        c_cell_face_idx = c_loc(cell_face_idx)
    endif

    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
        c_numabs = c_loc(numabs)
    endif

    c_cell_face = C_NULL_PTR
    if (associated(cell_face)) then
        c_cell_face = c_loc(cell_face)
    endif

    c_face_som_nb = C_NULL_PTR
    if (associated(face_som_nb)) then
      c_face_som_nb = c_loc(face_som_nb)
    endif

    c_cell_face_nb = C_NULL_PTR
    if (associated(cell_face_nb)) then
      c_cell_face_nb = c_loc(cell_face_nb)
    endif

    call PDM_writer_geom_cell3d_cellface_add_c(cs,              &
                                               id_geom,         &
                                               id_part,         &
                                               n_cell,          &
                                               n_face,          &
                                               c_face_som_idx,  &
                                               c_face_som_nb,   &
                                               c_face_som,      &
                                               c_cell_face_idx, &
                                               c_cell_face_nb,  &
                                               c_cell_face,     &
                                               c_numabs)

  end subroutine PDM_writer_geom_cell3d_cellface_add



  subroutine PDM_writer_geom_cell2d_cellface_add(cs,            &
                                                 id_geom,       &
                                                 id_part,       &
                                                 n_cell,        &
                                                 n_face,        &
                                                 face_som_idx,  &
                                                 face_som_nb,   &
                                                 face_som,      &
                                                 cell_face_idx, &
                                                 cell_face_nb,  &
                                                 cell_face,     &
                                                 numabs)
    ! Add 2D cells described in terms of faces
    implicit none

    type(c_ptr),                   intent(in) :: cs               ! Pointer to PDM_writer_t instance
    integer,                       intent(in) :: id_geom          ! Geometry identifier
    integer,                       intent(in) :: id_part          ! Partition identifier
    integer,                       intent(in) :: n_cell           ! Number of 2D cells
    integer,                       intent(in) :: n_face           ! Number of faces
    integer(pdm_l_num_s), pointer, intent(in) :: face_som_idx(:)  ! Index of the Face->Vertex connectivity (unused)
    integer(pdm_l_num_s), pointer, intent(in) :: face_som_nb(:)   ! Number of vertices per face (unused)
    integer(pdm_l_num_s), pointer, intent(in) :: face_som(:)      ! Face->Vertex connectivity (size = 2 * n_face)
    integer(pdm_l_num_s), pointer, intent(in) :: cell_face_idx(:) ! Index of the Cell->Face connectivity (size = n_cell + 1)
    integer(pdm_l_num_s), pointer, intent(in) :: cell_face_nb(:)  ! Number of faces per cell (optional, can be set to *null()*)
    integer(pdm_l_num_s), pointer, intent(in) :: cell_face(:)     ! Cell->Face connectivity (size = cell_face_idx(n_cell+1))
    integer(pdm_g_num_s), pointer, intent(in) :: numabs(:)        ! Cell global IDs (size = n_cell)

    type(c_ptr)                               :: c_face_som_idx
    type(c_ptr)                               :: c_face_som_nb
    type(c_ptr)                               :: c_face_som
    type(c_ptr)                               :: c_cell_face_idx
    type(c_ptr)                               :: c_cell_face_nb
    type(c_ptr)                               :: c_cell_face
    type(c_ptr)                               :: c_numabs

    interface
      subroutine PDM_writer_geom_cell2d_cellface_add_c(cs,            &
                                                       id_geom,       &
                                                       id_part,       &
                                                       n_cell,        &
                                                       n_face,        &
                                                       face_som_idx,  &
                                                       face_som_nb,   &
                                                       face_som,      &
                                                       cell_face_idx, &
                                                       cell_face_nb,  &
                                                       cell_face,     &
                                                       numabs)        &
      bind (c, name='PDM_writer_geom_cell2d_cellface_add')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom
        integer(c_int), value :: id_part
        integer(c_int), value :: n_cell
        integer(c_int), value :: n_face
        type(c_ptr),    value :: face_som_idx
        type(c_ptr),    value :: face_som_nb
        type(c_ptr),    value :: face_som
        type(c_ptr),    value :: cell_face_idx
        type(c_ptr),    value :: cell_face_nb
        type(c_ptr),    value :: cell_face
        type(c_ptr),    value :: numabs

      end subroutine PDM_writer_geom_cell2d_cellface_add_c
    end interface

    c_face_som_idx = C_NULL_PTR
    if (associated(face_som_idx)) then
      c_face_som_idx = c_loc(face_som_idx)
    endif

    c_face_som = C_NULL_PTR
    if (associated(face_som)) then
      c_face_som = c_loc(face_som)
    endif

    c_cell_face_idx = C_NULL_PTR
    if (associated(cell_face_idx)) then
      c_cell_face_idx = c_loc(cell_face_idx)
    endif

    c_face_som_nb = C_NULL_PTR
    if (associated(face_som_nb)) then
      c_face_som_nb = C_NULL_PTR
    endif

    c_cell_face = C_NULL_PTR
    if (associated(cell_face)) then
      c_cell_face = c_loc(cell_face)
    endif

    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs = c_loc(numabs)
    endif

    c_cell_face_nb = C_NULL_PTR
    if (associated(cell_face_nb)) then
      c_cell_face_nb = c_loc(cell_face_nb)
    endif

    call PDM_writer_geom_cell2d_cellface_add_c(cs,              &
                                               id_geom,         &
                                               id_part,         &
                                               n_cell,          &
                                               n_face,          &
                                               c_face_som_idx,  &
                                               c_face_som_nb,   &
                                               c_face_som,      &
                                               c_cell_face_idx, &
                                               c_cell_face_nb,  &
                                               c_cell_face,     &
                                               c_numabs)

  end subroutine PDM_writer_geom_cell2d_cellface_add



  subroutine PDM_writer_geom_faces_facesom_add(cs,           &
                                               id_geom,      &
                                               id_part,      &
                                               n_face,       &
                                               face_som_idx, &
                                               face_som_nb,  &
                                               face_som,     &
                                               numabs)
    ! Add faces described in nodal fashion
    implicit none

    type(c_ptr),                   intent(in) :: cs              ! Pointer to PDM_writer_t instance
    integer,                       intent(in) :: id_geom         ! Geometry identifier
    integer,                       intent(in) :: id_part         ! Partition identifier
    integer,                       intent(in) :: n_face          ! Number of faces
    integer(pdm_l_num_s), pointer, intent(in) :: face_som_idx(:) ! Index of the Face->Vertex connectivity (size = n_face + 1)
    integer(pdm_l_num_s), pointer, intent(in) :: face_som_nb(:)  ! Number of vertices per face (not used)
    integer(pdm_l_num_s), pointer, intent(in) :: face_som(:)     ! Face->Vertex connectivity (size = face_som_idx(n_face+1))
    integer(pdm_g_num_s), pointer, intent(in) :: numabs(:)       ! Face global IDs (size = n_face)

    type(c_ptr)                               :: c_face_som_idx
    type(c_ptr)                               :: c_face_som_nb
    type(c_ptr)                               :: c_face_som
    type(c_ptr)                               :: c_numabs

    interface
      subroutine PDM_writer_geom_faces_facesom_add_c(cs,           &
                                                     id_geom,      &
                                                     id_part,      &
                                                     n_face,       &
                                                     face_som_idx, &
                                                     face_som_nb,  &
                                                     face_som,     &
                                                     numabs)       &
      bind (c, name='PDM_writer_geom_faces_facesom_add')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom
        integer(c_int), value :: id_part
        integer(c_int), value :: n_face
        type(c_ptr),    value :: face_som_idx
        type(c_ptr),    value :: face_som_nb
        type(c_ptr),    value :: face_som
        type(c_ptr),    value :: numabs

      end subroutine PDM_writer_geom_faces_facesom_add_c
    end interface

    c_face_som_idx = C_NULL_PTR
    if (associated(face_som_idx)) then
      c_face_som_idx = c_loc(face_som_idx)
    endif

    c_face_som_nb = C_NULL_PTR
    ! Temporary fix before API change
    ! if (associated(face_som_nb)) then
    !   c_face_som_nb   = c_loc(face_som_nb)
    ! endif

    c_face_som = C_NULL_PTR
    if (associated(face_som)) then
      c_face_som = c_loc(face_som)
    endif

    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs = c_loc(numabs)
    endif


    call PDM_writer_geom_faces_facesom_add_c(cs,             &
                                             id_geom,        &
                                             id_part,        &
                                             n_face,         &
                                             c_face_som_idx, &
                                             c_face_som_nb,  &
                                             c_face_som,     &
                                             c_numabs)

  end subroutine PDM_writer_geom_faces_facesom_add



  subroutine PDM_writer_geom_write(cs,      &
                                   id_geom)
    ! Write current geometry
    implicit none

    type(c_ptr), intent(in) :: cs      ! Pointer to PDM_writer_t instance
    integer,     intent(in) :: id_geom ! Geometry identifier

    interface
      subroutine PDM_writer_geom_write_c(cs,      &
                                         id_geom) &
        bind (c, name='PDM_writer_geom_write')

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom

      end subroutine PDM_writer_geom_write_c
    end interface

    call PDM_writer_geom_write_c(cs,      &
                                 id_geom)

  end subroutine PDM_writer_geom_write



  subroutine PDM_writer_var_create(cs,         &
                                   id_var,     &
                                   st_dep_tps, &
                                   dim,        &
                                   loc,        &
                                   nom_var)
    ! Create a variable
    implicit none

    type(c_ptr),      intent(in)  :: cs         ! Pointer to PDM_writer_t instance
    integer,          intent(out) :: id_var     ! Variable identifier
    integer,          intent(in)  :: st_dep_tps ! Indicates whether the variable is time dependent
    integer,          intent(in)  :: dim        ! Variable's dimension
    integer,          intent(in)  :: loc        ! Variable's location
    character(len=*), intent(in)  :: nom_var    ! Name of the variable

    interface
      function PDM_writer_var_create_c(cs,         &
                                       st_dep_tps, &
                                       dim,        &
                                       loc,        &
                                       nom_var)    &
      result (id_var)                              &
      bind (c, name='PDM_writer_var_create')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int)        :: id_var
        integer(c_int), value :: st_dep_tps
        integer(c_int), value :: dim
        integer(c_int), value :: loc
        character(c_char)     :: nom_var(*)

      end function PDM_writer_var_create_c
    end interface

    id_var = PDM_writer_var_create_c(cs,                         &
                                     st_dep_tps,                 &
                                     dim,                        &
                                     loc,                        &
                                     trim(nom_var)//C_NULL_CHAR)

  end subroutine PDM_writer_var_create



  subroutine PDM_writer_cst_global_var_create(cs,      &
                                              id_var,  &
                                              nom_var, &
                                              val_var)
    ! Create a global constant variable
    implicit none

    type(c_ptr),      intent(in)  :: cs      ! Pointer to PDM_writer_t instance
    integer,          intent(out) :: id_var  ! Variable identifier
    character(len=*), intent(in)  :: nom_var ! Variable name
    real(c_double),   intent(in)  :: val_var ! Variable value

    interface
      function PDM_writer_cst_global_var_create_c(cs,         &
                                                  nom_var,    &
                                                  var_val)    &
      result (id_var)                                         &
      bind (c, name='PDM_writer_cst_global_var_create')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int)        :: id_var
        character(c_char)     :: nom_var(*)
        real(c_double)        :: var_val

      end function PDM_writer_cst_global_var_create_c
    end interface

    id_var = PDM_writer_cst_global_var_create_c(cs,                   &
                                                nom_var//C_NULL_CHAR, &
                                                val_var)

  end subroutine PDM_writer_cst_global_var_create



  subroutine PDM_writer_cst_global_var_set(cs,      &
                                           id_var,  &
                                           val_var)
    ! Set a global constant variable
    implicit none
    type(c_ptr),      intent(in)  :: cs      ! Pointer to PDM_writer_t instance
    integer,          intent(out) :: id_var  ! Variable identifier
    real(c_double),   intent(in)  :: val_var ! Variable value

    interface
      subroutine PDM_writer_cst_global_var_set_c(cs,      &
                                                 id_var,  &
                                                 var_val) &
      bind (c, name='PDM_writer_cst_global_var_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int)        :: id_var
        real(c_double)        :: var_val

      end subroutine PDM_writer_cst_global_var_set_c
    end interface

    call PDM_writer_cst_global_var_set_c(cs,      &
                                         id_var,  &
                                         val_var)

  end subroutine PDM_writer_cst_global_var_set



  subroutine PDM_writer_name_map_add (cs,           &
                                      public_name,  &
                                      private_name)
    ! Variable name mapping
    implicit none

    type(c_ptr),      intent(in) :: cs           ! Pointer to PDM_writer_t instance
    character(len=*), intent(in) :: private_name ! Public variable name
    character(len=*), intent(in) :: public_name  ! Private variable name

    interface
      subroutine PDM_writer_name_map_add_c(cs,           &
                                           public_name,  &
                                           private_name) &
      bind (c, name='PDM_writer_name_map_add')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        character(c_char)     :: public_name(*)
        character(c_char)     :: private_name(*)

      end subroutine PDM_writer_name_map_add_c
    end interface

    call PDM_writer_name_map_add_c(cs,                              &
                                   trim(public_name)//C_NULL_CHAR,  &
                                   trim(private_name)//C_NULL_CHAR)

  end subroutine PDM_writer_name_map_add



  subroutine PDM_writer_var_set(cs,      &
                                id_var,  &
                                id_geom, &
                                id_part, &
                                val)
    ! Update variable values
    !
    ! ..note:: The values are hard-copied (bufferization)
    implicit none

    type(c_ptr),             intent(in) :: cs      ! Pointer to PDM_writer_t instance
    integer,                 intent(in) :: id_var  ! Variable identifier
    integer,                 intent(in) :: id_geom ! Geometry identifier
    integer,                 intent(in) :: id_part ! Partition identifier
    real(c_double), pointer, intent(in) :: val(:)  ! Variable values

    type(c_ptr)                         :: c_val

    interface
      subroutine PDM_writer_var_set_c(cs,      &
                                      id_var,  &
                                      id_geom, &
                                      id_part, &
                                      val)     &
      bind (c, name='PDM_writer_var_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_var
        integer(c_int), value :: id_geom
        integer(c_int), value :: id_part
        type(c_ptr),    value :: val

      end subroutine PDM_writer_var_set_c
    end interface

    c_val = C_NULL_PTR
    if (associated(val)) then
      c_val = c_loc(val)
    endif

    call PDM_writer_var_set_c(cs,      &
                              id_var,  &
                              id_geom, &
                              id_part, &
                              c_val)

  end subroutine PDM_writer_var_set


  subroutine PDM_writer_var_write(cs,     &
                                  id_var)
    ! Write variable values
    implicit none

    type(c_ptr), intent(in) :: cs     ! Pointer to PDM_writer_t instance
    integer,     intent(in) :: id_var ! Variable identifier

    interface
      subroutine PDM_writer_var_write_c(cs,     &
                                        id_var) &
        bind (c, name='PDM_writer_var_write')

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_var

      end subroutine PDM_writer_var_write_c
    end interface

    call PDM_writer_var_write_c(cs,     &
                                id_var)

  end subroutine PDM_writer_var_write



  !>
  !! \brief Add a writer format
  !!
  !! Define a new format writer
  !! WARNING: has not been tested, not sure about procedure pointer interoperability
  !!
  !! \param [in] name            Name
  !! \param [in] create_fct      Customize \ref PDM_writer_create function for the new format  (or NULL)
  !! \param [in] free_fct        Customize \ref PDM_writer_free function for the new format (or NULL)
  !! \param [in] beg_step_fct    Customize \ref PDM_writer_step_beg function for the new format (or NULL)
  !! \param [in] end_step_fct    Customize \ref PDM_writer_step_end function for the new format (or NULL)
  !! \param [in] geom_create_fct Customize \ref PDM_writer_geom_create function for the new format (or NULL)
  !! \param [in] geom_write_fct  Customize \ref PDM_writer_geom_write function for the new format
  !! \param [in] geom_free_fct   Customize \ref PDM_writer_geom_free function for the new format (or NULL)
  !! \param [in] var_create_fct  Customize \ref PDM_writer_var_create function for the new format (or NULL)
  !! \param [in] var_write_fct   Customize \ref PDM_writer_var_write function for the new format
  !! \param [in] var_free_fct    Customize \ref PDM_writer_var_free function for the new format (or NULL)
  !!
  !!

  subroutine PDM_writer_fmt_add (name,            &
                                 create_fct,      &
                                 free_fct,        &
                                 beg_step_fct,    &
                                 end_step_fct,    &
                                 geom_create_fct, &
                                 geom_write_fct,  &
                                 geom_free_fct,   &
                                 var_create_fct,  &
                                 var_write_fct,   &
                                 var_free_fct)
    use iso_c_binding
    implicit none

    character (len=*)    :: name
    procedure(), pointer :: create_fct
    procedure(), pointer :: free_fct
    procedure(), pointer :: beg_step_fct
    procedure(), pointer :: end_step_fct
    procedure(), pointer :: geom_create_fct
    procedure(), pointer :: geom_write_fct
    procedure(), pointer :: geom_free_fct
    procedure(), pointer :: var_create_fct
    procedure(), pointer :: var_write_fct
    procedure(), pointer :: var_free_fct

    type(c_funptr)       :: c_create_fct
    type(c_funptr)       :: c_free_fct
    type(c_funptr)       :: c_beg_step_fct
    type(c_funptr)       :: c_end_step_fct
    type(c_funptr)       :: c_geom_create_fct
    type(c_funptr)       :: c_geom_write_fct
    type(c_funptr)       :: c_geom_free_fct
    type(c_funptr)       :: c_var_create_fct
    type(c_funptr)       :: c_var_write_fct
    type(c_funptr)       :: c_var_free_fct

    interface
      subroutine PDM_writer_fmt_add_c (name,            &
                                       create_fct,      &
                                       free_fct,        &
                                       beg_step_fct,    &
                                       end_step_fct,    &
                                       geom_create_fct, &
                                       geom_write_fct,  &
                                       geom_free_fct,   &
                                       var_create_fct,  &
                                       var_write_fct,   &
                                       var_free_fct)    &
      bind (c, name='PDM_writer_fmt_add')
        use iso_c_binding
        implicit none

        character(c_char) :: name(*)
        type(c_funptr)    :: create_fct
        type(c_funptr)    :: free_fct
        type(c_funptr)    :: beg_step_fct
        type(c_funptr)    :: end_step_fct
        type(c_funptr)    :: geom_create_fct
        type(c_funptr)    :: geom_write_fct
        type(c_funptr)    :: geom_free_fct
        type(c_funptr)    :: var_create_fct
        type(c_funptr)    :: var_write_fct
        type(c_funptr)    :: var_free_fct

      end subroutine PDM_writer_fmt_add_c
    end interface

    c_create_fct      = c_funloc(create_fct)
    c_free_fct        = c_funloc(free_fct)
    c_beg_step_fct    = c_funloc(beg_step_fct)
    c_end_step_fct    = c_funloc(end_step_fct)
    c_geom_create_fct = c_funloc(geom_create_fct)
    c_geom_write_fct  = c_funloc(geom_write_fct)
    c_geom_free_fct   = c_funloc(geom_free_fct)
    c_var_create_fct  = c_funloc(var_create_fct)
    c_var_write_fct   = c_funloc(var_write_fct)
    c_var_free_fct    = c_funloc(var_free_fct)

    call PDM_writer_fmt_add_c (trim(name)//C_NULL_CHAR, &
                               c_create_fct,      &
                               c_free_fct,        &
                               c_beg_step_fct,    &
                               c_end_step_fct,    &
                               c_geom_create_fct, &
                               c_geom_write_fct,  &
                               c_geom_free_fct,   &
                               c_var_create_fct,  &
                               c_var_write_fct,   &
                               c_var_free_fct)

  end subroutine PDM_writer_fmt_add



  subroutine PDM_writer_var_data_free(cs,     &
                                      id_var)
    ! Free variable data arrays
    implicit none

    type(c_ptr), intent(in) :: cs     ! Pointer to PDM_writer_t instance
    integer,     intent(in) :: id_var ! Variable identifier

    interface
      subroutine PDM_writer_var_data_free_c(cs,     &
                                            id_var) &
        bind (c, name='PDM_writer_var_data_free')

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_var

      end subroutine PDM_writer_var_data_free_c
    end interface

    call PDM_writer_var_data_free_c(cs,     &
                                    id_var)

  end subroutine PDM_writer_var_data_free



  subroutine PDM_writer_var_free(cs,     &
                                 id_var)
    ! Free the current variable
    implicit none

    type(c_ptr), intent(in) :: cs     ! Pointer to PDM_writer_t instance
    integer,     intent(in) :: id_var ! Variable identifier

    interface
      subroutine PDM_writer_var_free_c(cs,     &
                                       id_var) &
        bind (c, name='PDM_writer_var_free')

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_var

      end subroutine PDM_writer_var_free_c
    end interface

    call PDM_writer_var_free_c(cs,     &
                               id_var)

  end subroutine PDM_writer_var_free



  subroutine PDM_writer_geom_data_reset(cs,      &
                                        id_geom)
    ! Reset data describing the current geometry
    implicit none

    type(c_ptr), intent(in) :: cs      ! Pointer to PDM_writer_t instance
    integer,     intent(in) :: id_geom ! Geometry identifier

    interface
      subroutine PDM_writer_geom_data_reset_c(cs,      &
                                              id_geom) &
        bind (c, name='PDM_writer_geom_data_reset')

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom

      end subroutine PDM_writer_geom_data_reset_c
    end interface

    call PDM_writer_geom_data_reset_c(cs,      &
                                      id_geom)

  end subroutine PDM_writer_geom_data_reset



  subroutine PDM_writer_geom_data_free(cs,      &
                                       id_geom)
    ! Free data describing the current geometry
    !
    ! Indirections on global IDs are retained
    implicit none

    type(c_ptr), intent(in) :: cs      ! Pointer to PDM_writer_t instance
    integer,     intent(in) :: id_geom ! Geometry identifier

    interface
      subroutine PDM_writer_geom_data_free_c(cs,      &
                                             id_geom) &
        bind (c, name='PDM_writer_geom_data_free')

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom

      end subroutine PDM_writer_geom_data_free_c
    end interface

    call PDM_writer_geom_data_free_c(cs,      &
                                     id_geom)

  end subroutine PDM_writer_geom_data_free



  subroutine PDM_writer_geom_free(cs,      &
                                  id_geom)
    ! Free the current geometry
    implicit none

    type(c_ptr), intent(in) :: cs      ! Pointer to PDM_writer_t instance
    integer,     intent(in) :: id_geom ! Geometry identifier

    interface
      subroutine PDM_writer_geom_free_c(cs,      &
                                        id_geom) &
        bind (c, name='PDM_writer_geom_free')

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom

      end subroutine PDM_writer_geom_free_c
    end interface

    call PDM_writer_geom_free_c(cs,      &
                                id_geom)

  end subroutine PDM_writer_geom_free



  subroutine PDM_writer_free(cs)
    ! Free a PDM_writer_t instance
    implicit none

    type(c_ptr), intent(inout) :: cs ! Pointer to PDM_writer_t instance

    interface
      subroutine PDM_writer_free_c(cs) &
        bind (c, name='PDM_writer_free')
        use iso_c_binding
        implicit none
        type(c_ptr), value :: cs
      end subroutine PDM_writer_free_c
    end interface

    call PDM_writer_free_c(cs)

  end subroutine PDM_writer_free

end module pdm_writer
