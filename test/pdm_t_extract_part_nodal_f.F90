#include "pdm_configf.h"

program extract_part_nodal_f

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use pdm_generate_mesh
  use pdm_gnum
  use pdm_part_mesh_nodal
  use pdm_extract_part
  use pdm_part_to_part
  use pdm_io
  use pdm_writer
  use iso_c_binding

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !---------------------------------------------------------------
  integer, parameter            :: comm = MPI_COMM_WORLD             ! MPI Communicator
  integer                       :: i_rank, ierr

  integer                       :: i_arg
  character(len=99)             :: arg

  integer                       :: extract_kind                      ! Extraction kind (LOCAL or REEQUILIBRATE)
  integer                       :: part_method                       ! Re-partitioning method kind (if REEQUILIBRATE)
  logical                       :: visu                              ! Enable output for visualization

  integer                       :: dim                               ! Mesh dimension
  integer                       :: n_part                            ! Number of partitions on current MPI rank

  integer                       :: n_vtx                             ! Number of vertices
  double precision,     pointer :: vtx_coord(:,:)  => null()         ! Vertex coordinates
  integer(pdm_g_num_s), pointer :: vtx_ln_to_gn(:) => null()         ! Vertex global IDs
  integer                       :: n_tet                             ! Number of tetra
  integer(pdm_l_num_s), pointer :: tet_vtx(:)      => null()         ! Tetra->vtx connectivity
  integer(pdm_g_num_s), pointer :: tet_ln_to_gn(:) => null()         ! Tetra global IDs
  integer(pdm_l_num_s), pointer :: parent_num(:)   => null()         ! Not used
  integer(pdm_g_num_s), pointer :: parent_gnum(:)  => null()         ! Not used
  integer(pdm_l_num_s), pointer :: tet_vtx_idx(:)  => null()         ! Not used

  type(c_ptr)                   :: mesh = C_NULL_PTR                 ! Initial PartMeshNodal (ParaDiGM object)
  integer                       :: id_tet                            ! ID of tetra section

  type(c_ptr)                   :: extrp = C_NULL_PTR                ! ExtractPart (ParaDiGM object)
  integer                       :: n_selected                        ! Number of extracted elements
  integer(pdm_l_num_s), pointer :: selected(:) => null()             ! List of extracted elements

  type(c_ptr)                   :: extract_mesh = C_NULL_PTR         ! Extracted PartMeshNodal (ParaDiGM object)

  integer                       :: extract_n_vtx                     ! Number of extracted vertices
  double precision,     pointer :: extract_vtx_coord(:,:)  => null() ! Extracted vertex coordinates
  integer(pdm_g_num_s), pointer :: extract_vtx_ln_to_gn(:) => null() ! Extracted vertex global IDs
  integer                       :: extract_n_tet                     ! Number of extracted tetra
  integer(pdm_l_num_s), pointer :: extract_tet_vtx(:)      => null() ! Extracted tetra->vtx connectivity
  integer(pdm_g_num_s), pointer :: extract_tet_ln_to_gn(:) => null() ! Extracted tetra global IDs

  integer(pdm_l_num_s), pointer :: extract_vtx_parent(:)   => null() ! Extracted vertices IDs
  double precision,     pointer :: vtx_field(:)            => null() ! Field at vertices
  double precision,     pointer :: tet_field(:)            => null() ! Field at tetrahedra
  double precision,     pointer :: extract_vtx_field(:)    => null() ! Field at extracted vertices
  double precision,     pointer :: extract_tet_field(:)    => null() ! Field at extracted tetrahedra

  integer                       :: i_tet, i_vtx, i
  double precision              :: min_x, max_x
  !---------------------------------------------------------------

  extract_kind = PDM_EXTRACT_PART_KIND_LOCAL
  part_method  = PDM_SPLIT_DUAL_WITH_HILBERT
  visu         = .false.

  !----------------------------------------
  ! Parse command line arguments
  i_arg = 1
  do while (i_arg <= command_argument_count())
    call get_command_argument(i_arg, arg)

    select case(arg)

      case ("-reequilibrate")
        extract_kind = PDM_EXTRACT_PART_KIND_REEQUILIBRATE

      case ("-part_method")
        i_arg = i_arg + 1
        call get_command_argument(i_arg, arg)
        read(arg, *) part_method

      case ("-visu")
        visu = .true.

    endselect

    i_arg = i_arg + 1
  enddo
  !----------------------------------------

  !----------------------------------------
  ! Initialize MPI
  call mpi_init(ierr)
  call mpi_comm_rank(comm, i_rank, ierr)
  !----------------------------------------


  !----------------------------------------
  ! Generate mesh
  call PDM_generate_mesh_ball_simplified(comm,        & ! <- MPI communicator
                                         n_vtx,       & ! -> Local number of vertices
                                         n_tet,       & ! -> Local number of tetra
                                         vtx_coord,   & ! -> Coordinated of local vertices
                                         tet_vtx_idx, & ! -> Index of tetra->vtx connectivity (not used here)
                                         tet_vtx)       ! -> Local tetra->vtx connectivity
  call PDM_fortran_free_c(c_loc(tet_vtx_idx))

  ! Generate global IDs of vtx and tetra (from coordinates)
  call generate_gnums(comm,         & ! <- MPI communicator
                      n_vtx,        & ! <- Local number of vertices
                      n_tet,        & ! <- Local number of tetra
                      vtx_coord,    & ! <- Coordinated of local vertices
                      tet_vtx,      & ! <- Local tetra->vtx connectivity
                      vtx_ln_to_gn, & ! -> Global IDs of local vertices
                      tet_ln_to_gn)   ! -> Global IDs of local tetra
  !----------------------------------------


  !----------------------------------------
  ! Create a PartMeshNodal instance to group all mesh data into a single structure.
  n_part = 1 ! Number of partitions (subdomains) on current MPI rank
  dim    = 3 ! Mesh dimension (volume)
  call PDM_part_mesh_nodal_create(mesh,   & ! -> PartMeshNodal instance
                                  dim,    & ! <- Mesh dimension
                                  n_part, & ! <- Number of subdomains on current MPI rank
                                  comm)     ! <- MPI communicator
  ! Set vertices
  call PDM_part_mesh_nodal_coord_set(mesh,               & ! <- PartMeshNodal instance
                                     0,                  & ! <- ID of current subdomain (i_part)
                                     n_vtx,              & ! <- Number of vertices in current subdomain
                                     vtx_coord,          & ! <- Coordinates of vertices in current subdomain
                                     PDM_OWNERSHIP_USER)   ! <- Ownership

  call PDM_part_mesh_nodal_vtx_gnum_set(mesh,               & ! <- PartMeshNodal instance
                                        0,                  & ! <- ID of current subdomain (i_part)
                                        vtx_ln_to_gn,       & ! <- Global IDs of vertices in current subdomain
                                        PDM_OWNERSHIP_USER)   ! <- Ownership

  ! Set tetrahedra
  id_tet = PDM_part_mesh_nodal_section_add(mesh,                  & ! <- PartMeshNodal instance
                                           PDM_MESH_NODAL_TETRA4)   ! <- Tetrahedra


  call PDM_part_mesh_nodal_section_std_set(mesh,               & ! <- PartMeshNodal instance
                                           id_tet,             & ! <- Tetra section ID
                                           0,                  & ! <- ID of current subdomain (i_part)
                                           n_tet,              & ! <- Number of tetrahedra in current subdomain
                                           tet_vtx,            & ! <- Connectivity tetra->vtx in current subdomain
                                           tet_ln_to_gn,       & ! <- Global IDs of tetra in current subdomain
                                           parent_num,         & ! <- null()
                                           parent_gnum,        & ! <- null()
                                           PDM_OWNERSHIP_USER)   ! <- Ownership
  !----------------------------------------


  !----------------------------------------
  ! Create ExtractPart instance
  call PDM_extract_part_create(extrp,              & ! -> ExtractPart instance
                               dim,                & ! <- Dimension
                               n_part,             & ! <- Number of partitions per rank (input)
                               n_part,             & ! <- Number of partitions per rank (output)
                               extract_kind,       & ! <- Extraction kind
                               part_method,        & ! <- Re-partitioning method (not used in LOCAL mode)
                               .true.,             & ! <- Generate global IDs restricted to extraction
                               PDM_OWNERSHIP_KEEP, & ! <- Ownership of extraction
                               comm)                 ! <- MPI communicator

  ! Set input mesh
  call PDM_extract_part_part_nodal_set(extrp, mesh)

  ! Select elements to extract
  allocate(selected(n_tet))
  n_selected = 0
  do i_tet = 1, n_tet
    min_x =  1.d30
    max_x = -1.d30
    do i = 1, 4
      i_vtx = tet_vtx(4*(i_tet-1) + i)
      min_x = min(min_x, vtx_coord(1,i_vtx))
      max_x = max(max_x, vtx_coord(1,i_vtx))
    enddo

    ! Select current tetra if it crosses the (x = 0) plane
    if (min_x <= 0.d0 .and. max_x >= 0.d0) then
      n_selected = n_selected + 1
      selected(n_selected) = i_tet
    endif
  enddo

  call PDM_extract_part_selected_lnum_set(extrp,              & ! <- ExtractPart instance
                                          0,                  & ! <- ID of current subdomain (i_part)
                                          n_selected,         & ! <- Local number of extracted elements
                                          selected,           & ! <- Local IDs of extracted elements
                                          PDM_OWNERSHIP_USER)   ! <- Ownership

  ! Realize extraction
  call PDM_extract_part_compute(extrp)


  ! Retrieve extracted mesh
  call PDM_extract_part_part_mesh_nodal_get(extrp,              & ! <- ExtractPart instance
                                            extract_mesh,       & ! -> Extracted mesh (PartMeshNodal)
                                            PDM_OWNERSHIP_USER)   ! <- Ownership (user takes ownerhip of 'extract_mesh')

  !  Vertices
  call PDM_part_mesh_nodal_n_vtx_get(extract_mesh, & ! <- Extracted mesh (PartMeshNodal)
                                     0,            & ! <- ID of current subdomain (i_part)
                                     extract_n_vtx)  ! -> Number of extracted vertices in current subdomain

  call PDM_part_mesh_nodal_vtx_coord_get(extract_mesh,       & ! <- Extracted mesh (PartMeshNodal)
                                         0,                  & ! <- ID of current subdomain (i_part)
                                         extract_vtx_coord,  & ! -> Coordinates of extracted vertices
                                         PDM_OWNERSHIP_KEEP)   ! <- Ownership ('extract_mesh' keeps ownership)

  call PDM_part_mesh_nodal_vtx_g_num_get(extract_mesh,         & ! <- Extracted mesh (PartMeshNodal)
                                         0,                    & ! <- ID of current subdomain (i_part)
                                         extract_vtx_ln_to_gn, & ! -> Global IDs of extracted vertices
                                         PDM_OWNERSHIP_KEEP)     ! <- Ownership ('extract_mesh' keeps ownership)

  !  Tetrahedra
  call PDM_part_mesh_nodal_section_n_elt_get(extract_mesh,  & ! <- Extracted mesh (PartMeshNodal)
                                             id_tet,        & ! <- Tetra section ID
                                             0,             & ! <- ID of current subdomain (i_part)
                                             extract_n_tet)   ! -> Number of extracted tetra in current subdomain

  if (extract_kind == PDM_EXTRACT_PART_KIND_LOCAL .and. extract_n_tet /= n_selected) then
    print *, "extract_n_tet =", extract_n_tet, " but expected", n_selected
    STOP
  endif

  call PDM_part_mesh_nodal_section_std_get(extract_mesh,         & ! <- Extracted mesh (PartMeshNodal)
                                           id_tet,               & ! <- Tetra section ID
                                           0,                    & ! <- ID of current subdomain (i_part)
                                           extract_tet_vtx,      & ! -> Extracted connectivity tetra->vtx in current subdomain
                                           extract_tet_ln_to_gn, & ! -> Global IDs of extracted tetra in current subdomain
                                           parent_num,           & ! -> null()
                                           parent_gnum,          & ! -> null()
                                           PDM_OWNERSHIP_KEEP)     ! <- Ownership ('extract_mesh' keeps ownership)
  !----------------------------------------


  !----------------------------------------
  ! Transfer data from initial mesh to extraction

  !  Field at vertices
  allocate(vtx_field(n_vtx))
  allocate(extract_vtx_field(extract_n_vtx))

  vtx_field(:) = cos(4*(vtx_coord(1,:) + vtx_coord(2,:) + vtx_coord(3,:)))

  if (extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) then
    ! Local transfer
    call PDM_extract_part_parent_lnum_get(extrp,               & ! <- ExtractPart instance
                                          0,                   & ! <- ID of current subdomain (i_part)
                                          PDM_MESH_ENTITY_VTX, & ! <- Vertices
                                          extract_n_vtx,       & ! -> Number of extracted vertices in current subdomain
                                          extract_vtx_parent,  & ! -> Local IDs of extracted vertices in initial mesh
                                          PDM_OWNERSHIP_KEEP)    ! <- Ownership ('extrp' keeps ownership)

    extract_vtx_field(1:extract_n_vtx) = vtx_field(extract_vtx_parent(1:extract_n_vtx))
  else
    ! Reequilibrate mode => parallel transfer
    call data_transfer(extrp,               & ! <- ExtractPart instance
                       PDM_MESH_ENTITY_VTX, & ! <- Vertices
                       vtx_field,           & ! <- Field at vertices (local to current subdomain)
                       extract_n_vtx,       & ! <- Number of extracted vertices in current subdomain
                       extract_vtx_field)     ! -> Field at extracted vertices
  endif


  !  Field at tetrahedra : take average of field values at vertices
  allocate(tet_field(n_tet))
  allocate(extract_tet_field(extract_n_tet))

  tet_field(:) = 0.25d0*(vtx_field(tet_vtx(1::4)) + &
                         vtx_field(tet_vtx(2::4)) + &
                         vtx_field(tet_vtx(3::4)) + &
                         vtx_field(tet_vtx(4::4)))

  if (extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) then
    ! Local transfer
    extract_tet_field(1:n_selected) = tet_field(selected(1:n_selected))
  else
    ! Reequilibrate mode => parallel transfer
    call data_transfer(extrp,                & ! <- ExtractPart instance
                       PDM_MESH_ENTITY_CELL, & ! <- Cells (tetra)
                       tet_field,            & ! <- Field at tetra (local to current subdomain)
                       extract_n_tet,        & ! <- Number of extracted tetra in current subdomain
                       extract_tet_field)      ! -> Field at extracted tetra
  endif
  !----------------------------------------


  !----------------------------------------
  ! Visu
  if (visu) then
    call visu_pmn(comm,      & ! <- MPI communicator
                  mesh,      & ! <- Mesh (PartMeshNodal)
                  vtx_field, & ! <- Field at vertices (local to current subdomain)
                  tet_field, & ! <- Field at tetra (local to current subdomain)
                  "init")      ! <- Output name

    call visu_pmn(comm,              & ! <- MPI communicator
                  extract_mesh,      & ! <- Extracted mesh (PartMeshNodal)
                  extract_vtx_field, & ! <- Field at extracted vertices
                  extract_tet_field, & ! <- Field at extracted tetra
                  "extract")           ! <- Output name
  endif
  !----------------------------------------


  !----------------------------------------
  ! Free memory
  deallocate(vtx_field,         &
             tet_field,         &
             extract_vtx_field, &
             extract_tet_field)
  deallocate(selected)
  call PDM_part_mesh_nodal_free(extract_mesh)
  call PDM_extract_part_free(extrp)
  call PDM_part_mesh_nodal_free(mesh)

  call PDM_fortran_free_c(c_loc(vtx_coord))
  call PDM_fortran_free_c(c_loc(tet_vtx))
  call PDM_fortran_free_c(c_loc(vtx_ln_to_gn))
  call PDM_fortran_free_c(c_loc(tet_ln_to_gn))
  !----------------------------------------

  if (i_rank == 0) then
    print *, "The End :)"
  endif

  call mpi_finalize(ierr)


contains

  !--------------------------------------------------------------------------------
  ! Auxiliary routines to make code lighter and easier to read
  subroutine gnum_from_coord(comm,         &
                             n_elt,        &
                             elt_coord,    &
                             elt_ln_to_gn)
    ! Generate global IDs from coordinates (assume 1 part per rank)
    implicit none

    integer,                       intent(in)  :: comm            ! <- MPI communicator
    integer,                       intent(in)  :: n_elt           ! <- Local number of elements
    double precision,     pointer, intent(in)  :: elt_coord(:,:)  ! <- Element coordinates
    integer(pdm_g_num_s), pointer, intent(out) :: elt_ln_to_gn(:) ! -> Element global IDs

    type(c_ptr)                                :: gen_gnum
    double precision,     pointer              :: char_length(:)

    allocate(char_length(n_elt))
    char_length(:) = 1.d0

    call PDM_gnum_create(gen_gnum,           & ! -> C pointer to PDM_gen_gnum_t object
                         3,                  & ! <- Spatial dimension (not mesh element dimension!)
                         1,                  & ! <- Number of subdomains on current rank
                         PDM_TRUE,           & ! <- Coincident points get assigned the same global ID
                         1.d-6,              & ! <- Relative geometric tolerance
                         comm,               & ! <- MPI communicator
                         PDM_OWNERSHIP_USER)   ! <- Ownership (user owns result)

    call PDM_gnum_set_from_coords(gen_gnum,    & ! <- C pointer to PDM_gen_gnum_t object
                                  0,           & ! <- ID of current subdomain (i_part)
                                  n_elt,       & ! <- Local number of elements in current subdomain
                                  elt_coord,   & ! <- Coordinates
                                  char_length)   ! <- Element-based characteristic length

    call PDM_gnum_compute(gen_gnum)
    deallocate(char_length)

    call PDM_gnum_get(gen_gnum,     & ! <- C pointer to PDM_gen_gnum_t object
                      0,            & ! <- ID of current subdomain (i_part)
                      elt_ln_to_gn)   ! -> Global IDs

    call PDM_gnum_free(gen_gnum)

  end subroutine gnum_from_coord


  subroutine generate_gnums(comm,         &
                            n_vtx,        &
                            n_tet,        &
                            vtx_coord,    &
                            tet_vtx,      &
                            vtx_ln_to_gn, &
                            tet_ln_to_gn)
    ! Generate global IDs for vertices and tetrahedra (assume 1 part per rank)
    implicit none

    integer,                       intent(in)  :: comm            ! <- MPI communicator
    integer,                       intent(in)  :: n_vtx           ! <- Local number of vertices
    integer,                       intent(in)  :: n_tet           ! <- Local number of tetra
    double precision,     pointer, intent(in)  :: vtx_coord(:,:)  ! <- Coordinates of vertices
    integer(pdm_l_num_s), pointer, intent(in)  :: tet_vtx(:)      ! <- Tetra->vtx connectivity
    integer(pdm_g_num_s), pointer, intent(out) :: vtx_ln_to_gn(:) ! -> Vertex global IDs
    integer(pdm_g_num_s), pointer, intent(out) :: tet_ln_to_gn(:) ! -> Tetra global IDs

    double precision,     pointer              :: tet_coord(:,:)
    integer                                    :: i_tet, i

    ! Vertices
    call gnum_from_coord(comm, n_vtx, vtx_coord, vtx_ln_to_gn)

    ! Tetra
    allocate(tet_coord(3, n_tet))
    do i_tet = 1, n_tet
      tet_coord(:,i_tet) = 0.d0
      do i = 1, 4
        tet_coord(:,i_tet) = tet_coord(:,i_tet) + 0.25d0*vtx_coord(:,tet_vtx(4*(i_tet-1)+i))
      enddo
    enddo
    call gnum_from_coord(comm, n_tet, tet_coord, tet_ln_to_gn)
    deallocate(tet_coord)

  end subroutine generate_gnums


  subroutine visu_pmn(comm,      &
                      mesh,      &
                      vtx_field, &
                      elt_field, &
                      name)
    ! Export PartMeshNodal in Ensight format
    ! WARNING: does not work if 'mesh' contains multiple geometry kinds (eg. surface *and* volume)
    ! See below how to handle this scenario
    implicit none

    integer,                   intent(in) :: comm         ! <- MPI communicator
    type(c_ptr),               intent(in) :: mesh         ! <- Mesh structure (PartMeshNodal)
    double precision, pointer, intent(in) :: vtx_field(:) ! <- Field at local vertices
    double precision, pointer, intent(in) :: elt_field(:) ! <- Field at local elements
    character(len = *),        intent(in) :: name         ! <- Output file name

    type(c_ptr)                           :: wrt
    integer                               :: id_geom
    integer                               :: id_var_elt_part
    integer                               :: id_var_vtx_field
    integer                               :: id_var_elt_field
    integer                               :: n_elt
    double precision, pointer             :: val_elt_part(:) => null()

    integer                               :: i_rank, err

    call mpi_comm_rank(comm, i_rank, err)

    call PDM_writer_create(wrt,                    & ! -> Writer instance
                           "Ensight",              & ! <- Format
                           PDM_WRITER_FMT_BIN,     & ! <- Binary files
                           PDM_WRITER_TOPO_CST,    & ! <- Topology is constant over time
                           PDM_WRITER_OFF,         & ! <- Write from scratch
                           "extract_part_nodal_f", & ! <- Output directory
                           name,                   & ! <- Output file
                           comm,                   & ! <- MPI communicator
                           PDM_IO_KIND_MPI_SIMPLE, & ! <- MPI IO strategy (type of file access)
                           1.d0,                   & ! <- MPI IO strategy (proportion of active nodes)
                           "")                       ! <- No additional options

    ! WARNING: If multiple geometry kinds (eg. surface *and* volume),
    ! only the one of interest should be passed :
    !   PDM_writer_geom_create
    !   PDM_writer_geom_coord_set
    !   PDM_writer_geom_bloc_add
    !   PDM_writer_geom_bloc_std_set
    !   ...
    call PDM_writer_geom_create_from_mesh_nodal(wrt,     & ! <- Writer instance
                                                id_geom, & ! -> Geometry identifier
                                                name,    & ! <- Geometry name
                                                mesh)      ! <- Mesh structure

    ! Define a variable to visualize the partitioning
    call PDM_writer_var_create(wrt,                     & ! <- Writer instance
                               id_var_elt_part,         & ! -> Variable identifier
                               PDM_WRITER_OFF,          & ! <- Not time-dependent
                               PDM_WRITER_VAR_SCALAIRE, & ! <- Scalar variable
                               PDM_WRITER_VAR_ELEMENTS, & ! <- Element-based variable
                               "i_part")                  ! <- Variable name

    if (associated(vtx_field)) then
      ! Define a variable to visualize the field at vertices
      call PDM_writer_var_create(wrt,                     & ! <- Writer instance
                                 id_var_vtx_field,        & ! -> Variable identifier
                                 PDM_WRITER_OFF,          & ! <- Not time-dependent
                                 PDM_WRITER_VAR_SCALAIRE, & ! <- Scalar variable
                                 PDM_WRITER_VAR_VERTICES, & ! <- Vertex-based variable
                                 "vtx_field")               ! <- Variable name
    endif

    if (associated(elt_field)) then
      ! Define a variable to visualize the field at elements
      call PDM_writer_var_create(wrt,                     & ! <- Writer instance
                                 id_var_elt_field,        & ! -> Variable identifier
                                 PDM_WRITER_OFF,          & ! <- Not time-dependent
                                 PDM_WRITER_VAR_SCALAIRE, & ! <- Scalar variable
                                 PDM_WRITER_VAR_ELEMENTS, & ! <- Element-based variable
                                 "elt_field")               ! <- Variable name
    endif

    ! Begin a time-step (only one here)
    call PDM_writer_step_beg(wrt,  & ! <- Writer instance
                             0.d0)   ! <- Time

    ! Write the geometry
    call PDM_writer_geom_write(wrt,     & ! <- Writer instance
                               id_geom)   ! <- Geometry identifier

    ! Write 'i_part' variable
    call PDM_part_mesh_nodal_section_n_elt_get(mesh,  & ! <- Mesh structure
                                               0,     & ! <- Section identifier
                                               0,     & ! <- ID of current subdomain (i_part)
                                               n_elt)   ! -> Number of elements
    allocate(val_elt_part(n_elt))
    val_elt_part(:) = i_rank

    call PDM_writer_var_set(wrt,             & ! <- Writer instance
                            id_var_elt_part, & ! <- Variable identifier
                            id_geom,         & ! <- Geometry identifier
                            0,               & ! <- ID of current subdomain (i_part)
                            val_elt_part)      ! <- Variable values

    call PDM_writer_var_write(wrt,             & ! <- Writer instance
                              id_var_elt_part)   ! <- Variable identifier
    deallocate(val_elt_part)


    if (associated(vtx_field)) then
      ! Write 'vtx_field'
      call PDM_writer_var_set(wrt,              & ! <- Writer instance
                              id_var_vtx_field, & ! <- Variable identifier
                              id_geom,          & ! <- Geometry identifier
                              0,                & ! <- ID of current subdomain (i_part)
                              vtx_field)          ! <- Variable values

      call PDM_writer_var_write(wrt,              & ! <- Writer instance
                                id_var_vtx_field)   ! <- Variable identifier
    endif


    if (associated(elt_field)) then
      ! Write 'elt_field'
      call PDM_writer_var_set(wrt,              & ! <- Writer instance
                              id_var_elt_field, & ! <- Variable identifier
                              id_geom,          & ! <- Geometry identifier
                              0,                & ! <- ID of current subdomain (i_part)
                              elt_field)          ! <- Variable values

      call PDM_writer_var_write(wrt,              & ! <- Writer instance
                                id_var_elt_field)   ! <- Variable identifier
    endif

    ! End time-step
    call PDM_writer_step_end(wrt)

    ! Free Writer object
    call PDM_writer_free(wrt)

  end subroutine visu_pmn
  !--------------------------------------------------------------------------------


  !--------------------------------------------------------------------------------
  ! Transfer fields from initial mesh to extraction (for REEQUILIBRATE mode)
  subroutine data_transfer(extrp,         &
                           entity_type,   &
                           field,         &
                           n_extract,     &
                           extract_field)
    implicit none

    type(c_ptr),                   intent(in)  :: extrp             ! <- ExtractPart instance
    integer,                       intent(in)  :: entity_type       ! <- Type of mesh entity
    double precision,     pointer, intent(in)  :: field(:)          ! <- Field on local, initial mesh entities
    integer,                       intent(in)  :: n_extract         ! <- Local number of extracted entities
    double precision,     pointer, intent(out) :: extract_field(:)  ! -> Field on local, extracted mesh entities

    type(c_ptr)                                :: ptp
    type(PDM_pointer_array_t), pointer         :: pa_field
    type(PDM_pointer_array_t), pointer         :: pa_extract_field
    double precision,          pointer         :: tmp_extract_field(:)
    integer                                    :: request

    ! Initialize to null pointers
    ptp               =  C_NULL_PTR
    pa_field          => null()
    pa_extract_field  => null()
    tmp_extract_field => null()

    ! Get PartToPart instance
    call PDM_extract_part_part_to_part_get(extrp,              & ! <- ExtractPart instance
                                           entity_type,        & ! <- Type of mesh entity
                                           ptp,                & ! -> PartToPart instance
                                           PDM_OWNERSHIP_KEEP)   ! <- 'extrp' keeps ownership of 'ptp'

    ! Pointer array for input field
    call PDM_pointer_array_create(pa_field,        & ! -> Pointer array
                                  1,               & ! <- Number of subdomains on current rank
                                  PDM_TYPE_DOUBLE)   ! <- Data type

    call PDM_pointer_array_part_set(pa_field, & ! <- Pointer array
                                    0,        & ! <- ID of current subdomain (i_part)
                                    field)      ! <- Field

    ! Exchange
    call PDM_part_to_part_reverse_iexch(ptp,                                   & ! <- PartToPart instance
                                        PDM_MPI_COMM_KIND_P2P,                 & ! <- MPI communication kind
                                        PDM_STRIDE_CST_INTERLACED,             & ! <- Stride type
                                        PDM_PART_TO_PART_DATA_DEF_ORDER_PART2, & ! <- Data ordering
                                        1,                                     & ! <- Value of constant stride
                                        null(),                                & ! <- null()
                                        pa_field,                              & ! <- Field
                                        null(),                                & ! -> null()
                                        pa_extract_field,                      & ! -> Extracted field
                                        request)                                 ! -> MPI request

    ! Finalize exchange
    call PDM_part_to_part_reverse_iexch_wait(ptp,     & ! <- PartToPart instance
                                             request)   ! <- MPI request

    call PDM_pointer_array_part_get(pa_extract_field,   & ! <- Pointer array
                                    0,                  & ! <- ID of current subdomain (i_part)
                                    tmp_extract_field)    ! <- Extracted field

    ! Copy data to avoid ownership conflicts
    extract_field(1:n_extract) = tmp_extract_field(1:n_extract)

    ! Free pointer arrays
    call PDM_pointer_array_free(pa_field)
    call PDM_pointer_array_free(pa_extract_field)

  end subroutine data_transfer
  !--------------------------------------------------------------------------------

end program extract_part_nodal_f