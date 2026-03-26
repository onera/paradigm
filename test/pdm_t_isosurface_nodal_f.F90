#include "pdm_configf.h"

program isosurface_nodal_f

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use pdm_generate_mesh
  use pdm_gnum
  use pdm_part_mesh_nodal
  use pdm_isosurface
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

  type(c_ptr)                   :: mesh = C_NULL_PTR                 ! Initial PartMeshNodal (ParaDiGM object)
  integer                       :: dim                               ! Mesh dimension
  integer                       :: n_part                            ! Number of partitions on current MPI rank
  integer(kind=pdm_g_num_s)     :: n_seg=10, n_layer=0

  integer                       :: n_vtx                             ! Number of vertices
  double precision,     pointer :: vtx_coord(:,:)  => null()         ! Vertex coordinates
  integer                       :: n_tet                             ! Number of tetra
  integer(pdm_l_num_s), pointer :: tet_vtx_idx(:)  => null()         ! Tetra->vtx connectivity
  integer(pdm_l_num_s), pointer :: tet_vtx(:)      => null()         ! Tetra->vtx connectivity
  integer(pdm_g_num_s), pointer :: vtx_ln_to_gn(:) => null()         ! Vertex global IDs
  integer(pdm_g_num_s), pointer :: tet_ln_to_gn(:) => null()         ! Tetra global IDs


  type(c_ptr)                   :: isos              = C_NULL_PTR ! ExtractPart (ParaDiGM object)
  integer                       :: id_iso                         ! Number of extracted elements
  double precision, pointer     :: plane_equation(:) => null()    ! Coefficient (a, b, c, d) of plane equation a*x + b*y + c*z - d = 0
  double precision, pointer     :: isovalues(:)      => null()    ! List of isovalues to compute

  integer                       :: iso_n_vtx                      ! Number of vertices in the isosurface
  double precision,     pointer :: iso_vtx_coord(:,:)   => null() ! Isosurface vertex coordinates
  integer(pdm_g_num_s), pointer :: iso_vtx_ln_to_gn(:)  => null() ! Isosurface vertex global IDs
  integer                       :: iso_n_face                     ! Number of face in the isosurface
  integer(pdm_l_num_s), pointer :: iso_face_vtx_idx(:)  => null() ! Isosurface face->vtx connectivity index
  integer(pdm_l_num_s), pointer :: iso_face_vtx(:)      => null() ! Isosurface face->vtx connectivity
  integer(pdm_g_num_s), pointer :: iso_face_ln_to_gn(:) => null() ! Isosurface face global IDs

  integer(pdm_l_num_s), pointer :: piso_vtx_parent_vtx_idx(:) => null() ! vtx->vtx connectivity index between isosurface mesh and parent mesh
  integer(pdm_l_num_s), pointer :: piso_vtx_parent_vtx(:)     => null() ! vtx->vtx connectivity between isosurface mesh and parent mesh
  double precision,     pointer :: pvtx_parent_weight(:)      => null() ! interpolation weights of parent mesh

  integer(pdm_l_num_s), pointer :: piso_face_parent_cell_idx(:)    => null() ! face->cell connectivity index between isosurface mesh and parent mesh
  integer(pdm_l_num_s), pointer :: piso_face_parent_cell(:)        => null() ! face->cell connectivity between isosurface mesh and parent mesh

  double precision,     pointer :: vtx_field(:)      => null() ! Field at vertices
  double precision,     pointer :: tet_field(:)      => null() ! Field at tetrahedra
  double precision,     pointer :: iso_field(:)      => null() ! Field at extracted vertices
  double precision,     pointer :: iso_face_field(:) => null() ! Field at extracted faces

  integer                       :: i_vtx, i_parent, i_vtx_parent, i_face, i_face_parent


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

  n_part = 1 ! Number of partitions (subdomains) on current MPI rank
  dim    = 3 ! Mesh dimension (volume)

  mesh = PDM_generate_mesh_ball (comm,                   &
                                  PDM_MESH_NODAL_TETRA4, &
                                  1,                     &
                                  C_NULL_PTR,            &
                                  1.0d0,                 &
                                  0.0d0,                 &
                                  0.0d0,                 &
                                  0.0d0,                 &
                                  0.0d0,                 &
                                  n_seg,                 &
                                  n_seg,                 &
                                  n_seg,                 &
                                  n_layer,               &
                                  0.0d0,                 &
                                  n_part,                &
                                  part_method) 

  call PDM_part_mesh_nodal_n_vtx_get(mesh, &
                                     0,    &
                                     n_vtx)

  call PDM_part_mesh_nodal_vtx_coord_get(mesh,      &
                                         0,         &
                                         vtx_coord, &
                                         PDM_OWNERSHIP_USER)

  call PDM_part_mesh_nodal_cell_vtx_connect_get(mesh,                      &
                                                PDM_GEOMETRY_KIND_VOLUMIC, &
                                                0,                         &
                                                tet_vtx_idx,               &
                                                tet_vtx)

  call PDM_part_mesh_nodal_n_elmts_get(mesh,                      &
                                       PDM_GEOMETRY_KIND_VOLUMIC, &
                                       0,                         &
                                       n_tet)

  call generate_gnums(comm,         & ! <- MPI communicator
                      n_vtx,        & ! <- Local number of vertices
                      n_tet,        & ! <- Local number of tetra
                      vtx_coord,    & ! <- Coordinated of local vertices
                      tet_vtx,      & ! <- Local tetra->vtx connectivity
                      vtx_ln_to_gn, & ! -> Global IDs of local vertices
                      tet_ln_to_gn)   ! -> Global IDs of local tetra

  !----------------------------------------
  ! Create ExtractPart instance
  call PDM_isosurface_create(comm, & ! <- MPI communicator
                             3,    & ! <- Dimension
                             isos)   ! -> IsoSurface instance

  call PDM_isosurface_part_mesh_nodal_set(isos, mesh)

  call PDM_isosurface_redistribution_set (isos,         & ! <- IsoSurface instance
                                          extract_kind, & ! <- Reditribution kind
                                          part_method)    ! <- Re-partitioning method (not used in LOCAL mode)

  if (extract_kind == PDM_EXTRACT_PART_KIND_REEQUILIBRATE) then
    call PDM_isosurface_n_part_out_set(isos, n_part)
  end if

  ! Set up IsoSurface
  !  Planar slice in the mesh
  allocate(plane_equation(4))
  plane_equation = [-1.0d0, -1.0d0, 1.0d0, 0.0d0]
  allocate(isovalues(1))
  isovalues = [0.85d0]

  call PDM_isosurface_add (isos,                       & ! <- IsoSurface instance
                           PDM_ISO_SURFACE_KIND_PLANE, & ! <- Kind of isosurface
                           1,                          & ! <- Number of iso-value
                           isovalues,                  & ! <- List of iso-values
                           id_iso)                       ! -> ID of the isosurface

  call PDM_isosurface_equation_set(isos, &          ! <- IsoSurface instance
                                   id_iso, &        ! <- ID of the isosurface
                                   plane_equation); ! <- Coefficient of the plane equation (a*x + b*y + c*z - d = isovalues)


  ! Generate isosurface
  if (extract_kind == PDM_EXTRACT_PART_KIND_REEQUILIBRATE) then

    call PDM_isosurface_part_to_part_enable(isos,                & ! <- IsoSurface instance
                                            id_iso,              & ! <- ID of the isosurface
                                            PDM_MESH_ENTITY_VTX, & ! <- Mesh entity of the part_to_part
                                            0);

    call PDM_isosurface_part_to_part_enable(isos,                 & ! <- IsoSurface instance
                                            id_iso,               & ! <- ID of the isosurface
                                            PDM_MESH_ENTITY_FACE, & ! <- Mesh entity of the part_to_part
                                            0);    
  end if
  call PDM_isosurface_compute(isos, & ! <- IsoSurface instance
                              id_iso) ! <- ID of the isosurface (-1 to compute all isosurface at once)

  ! Retrieve isosurface connectivities
  call PDM_isosurface_pconnectivity_get(isos,                           & ! <- IsoSurface instance
                                        id_iso,                         & ! <- ID of the isosurface
                                        0,                              & ! <- ID of current part
                                        PDM_CONNECTIVITY_TYPE_FACE_VTX, & ! <- Face->vtx
                                        iso_n_face,                     & ! -> Number of isosurface faces
                                        iso_face_vtx_idx,               & ! -> Index for isosurface face->vtx connectivity
                                        iso_face_vtx,                   & ! -> Isosurface face->vtx connectivity
                                        PDM_OWNERSHIP_USER)               ! <- Ownership
  ! Coordinates
  call PDM_isosurface_pvtx_coord_get(isos,               & ! <- IsoSurface instance
                                     id_iso,             & ! <- ID of the isosurface
                                     0,                  & ! <- ID of current parts
                                     iso_n_vtx,          & ! -> Number of isosurface vertices
                                     iso_vtx_coord,      & ! -> Coordinates of isosurface vertices
                                     PDM_OWNERSHIP_USER)   ! <- Ownership

  ! Global IDs (at isosurface)
  call PDM_isosurface_ln_to_gn_get(isos,                 & ! <- IsoSurface instance
                                   id_iso,               & ! <- ID of the isosurface
                                   0,                    & ! <- ID of current part
                                   PDM_MESH_ENTITY_FACE, & ! <- Faces
                                   iso_n_face,           & ! -> Number of isosurface faces
                                   iso_face_ln_to_gn,    & ! -> Global IDs of isosurface faces
                                   PDM_OWNERSHIP_USER)     ! <- Ownership

  call PDM_isosurface_ln_to_gn_get(isos,                & ! <- IsoSurface instance
                                   id_iso,              & ! <- ID of the isosurface
                                   0,                   & ! <- ID of current part
                                   PDM_MESH_ENTITY_VTX, & ! <- Vertices
                                   iso_n_vtx,           & ! -> Number of isosurface vertices
                                   iso_vtx_ln_to_gn,    & ! -> Global IDs of isosurface vertices
                                   PDM_OWNERSHIP_USER)    ! <- Ownership


  !----------------------------------------
  ! Transfer data from initial mesh to isosurface mesh

  !  Field at vertices
  allocate(vtx_field(n_vtx))
  allocate(iso_field(iso_n_vtx))

  vtx_field(:) = cos(4*(vtx_coord(1,:) + vtx_coord(2,:) + vtx_coord(3,:)))

  if (extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) then
    ! Local transfer

    call  PDM_isosurface_pparent_weight_get(isos,                       & ! <- IsoSurface instance
                                            id_iso,                     & ! <- ID of the isosurface
                                            0,                          & ! <- ID of current part
                                            PDM_MESH_ENTITY_VTX,        & ! <- Entity mesh of isosurface : vertices
                                            iso_n_vtx,                  & ! -> Number of isosurface vertices
                                            piso_vtx_parent_vtx_idx,    & ! -> Index of connectivity vtx->vtx between isosurface and global mesh
                                            pvtx_parent_weight, &         ! -> Interpolation weights 
                                            PDM_OWNERSHIP_USER)           ! <- Ownership

    call PDM_isosurface_plocal_parent_get(isos,                    & ! <- IsoSurface instance
                                          id_iso,                  & ! <- ID of the isosurface
                                          0,                       & ! <- ID of current part
                                          PDM_MESH_ENTITY_VTX,     & ! <- Entity mesh of isosurface : vertices
                                          iso_n_vtx,               & ! -> Number of isosurface vertices
                                          piso_vtx_parent_vtx_idx, & ! -> Index of connectivity vtx->vtx between isosurface and global mesh
                                          piso_vtx_parent_vtx,     & ! -> Connectivity vtx->vtx between isosurface and global mesh
                                          PDM_OWNERSHIP_USER)        ! <- Ownership

    do i_vtx=1, iso_n_vtx
      iso_field(i_vtx) = 0.0d0
      do i_vtx_parent=piso_vtx_parent_vtx_idx(i_vtx), piso_vtx_parent_vtx_idx(i_vtx+1)-1
        i_parent = piso_vtx_parent_vtx(i_vtx_parent+1)
        iso_field(i_vtx) = iso_field(i_vtx) + pvtx_parent_weight(i_vtx_parent+1) * vtx_field(i_parent)
      end do
    enddo


  else
    ! Reequilibrate mode => parallel transfer
    call data_transfer(isos,                & ! <- IsoSurface instance
                       PDM_MESH_ENTITY_VTX, & ! <- Vertices
                       vtx_field,           & ! <- Field at vertices (local to current subdomain)
                       iso_n_vtx,           & ! <- Number of isosurface vertices in current subdomain
                       iso_field)             ! -> Field at isosurface vertices
  endif


  !  Field at tetrahedra : take average of field values at vertices
  allocate(tet_field(n_tet))
  allocate(iso_face_field(iso_n_face))

  tet_field(:) = 0.25d0*(vtx_field(tet_vtx(1::4)) + &
                         vtx_field(tet_vtx(2::4)) + &
                         vtx_field(tet_vtx(3::4)) + &
                         vtx_field(tet_vtx(4::4)))

  if (extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) then
    ! Local transfer
    call PDM_isosurface_plocal_parent_get(isos,                      & ! <- IsoSurface instance
                                          id_iso,                    & ! <- ID of the isosurface
                                          0,                         & ! <- ID of current part
                                          PDM_MESH_ENTITY_FACE,      & ! <- Entity mesh of isosurface : faces
                                          iso_n_face,                & ! -> Number of isosurface faces
                                          piso_face_parent_cell_idx, & ! -> Index of connectivity face->cell between isosurface and global mesh
                                          piso_face_parent_cell,     & ! -> Connectivity face->cell between isosurface and global mesh
                                          PDM_OWNERSHIP_USER)          ! <- Ownership

    do i_face=1, iso_n_face
      iso_face_field(i_face) = 0.0d0
      do i_face_parent=piso_face_parent_cell_idx(i_face), piso_face_parent_cell_idx(i_face+1)-1
        i_parent = piso_face_parent_cell(i_face_parent+1)
        iso_face_field(i_face) = iso_face_field(i_face) + tet_field(i_parent)
      end do 
    enddo
  else
    ! Reequilibrate mode => parallel transfer
    call data_transfer(isos,                 & ! <- IsoSurface instance
                       PDM_MESH_ENTITY_FACE, & ! <- Faces
                       tet_field,            & ! <- Field at tetra (local to current subdomain)
                       iso_n_face,           & ! <- Number of isosurface faces in current subdomain
                       iso_face_field)         ! -> Field at isosurface faces

  endif
  !----------------------------------------

  !----------------------------------------
  ! Visu
  if (visu) then
    call visu_mesh(comm,        &
                  "init",       &
                  3,            &
                  n_tet,        &
                  n_vtx,        &
                  tet_vtx_idx,  &
                  tet_vtx,      &
                  vtx_coord,    &
                  tet_ln_to_gn, &
                  vtx_ln_to_gn, &
                  vtx_field,    &
                  tet_field)

    call visu_mesh(comm,             &
                  "isosurface",      &
                  2,                 &
                  iso_n_face,        &
                  iso_n_vtx,         &
                  iso_face_vtx_idx,  &
                  iso_face_vtx,      &
                  iso_vtx_coord,     &
                  iso_face_ln_to_gn, &
                  iso_vtx_ln_to_gn,  &
                  iso_field,         &
                  iso_face_field)
  endif
  !----------------------------------------


  !----------------------------------------
  ! Free memory
  deallocate(vtx_field, &
             tet_field, &
             iso_field, &
             iso_face_field)
  deallocate(plane_equation, isovalues)
  call PDM_isosurface_free(isos)
  call PDM_part_mesh_nodal_free(mesh)

  call PDM_fortran_free_c(c_loc(vtx_coord))
  call PDM_fortran_free_c(c_loc(tet_vtx_idx))
  call PDM_fortran_free_c(c_loc(tet_vtx))
  call PDM_fortran_free_c(c_loc(vtx_ln_to_gn))
  call PDM_fortran_free_c(c_loc(tet_ln_to_gn))

  call PDM_fortran_free_c(c_loc(iso_vtx_coord))
  call PDM_fortran_free_c(c_loc(iso_face_vtx_idx))
  call PDM_fortran_free_c(c_loc(iso_face_vtx))
  call PDM_fortran_free_c(c_loc(iso_vtx_ln_to_gn))
  call PDM_fortran_free_c(c_loc(iso_face_ln_to_gn))  
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


  subroutine visu_mesh(comm,         &
                       name,         &
                       dim,          &
                       n_elt,        &
                       n_vtx,        &
                       elt_vtx_idx,  &
                       elt_vtx,      &
                       vtx_coord,    &
                       elt_ln_to_gn, &
                       vtx_ln_to_gn, &
                       vtx_field,    &
                       elt_field)

    implicit none

    integer,                            intent(in) :: comm            ! <- MPI communicator
    integer,                            intent(in) :: dim             ! <- Mesh dimension
    character(len = *),                 intent(in) :: name            ! <- Output file name
    integer,                            intent(in) :: n_elt           ! <- Number of element
    integer,                            intent(in) :: n_vtx           ! <- Number of vertices
    integer,                   pointer, intent(in) :: elt_vtx_idx(:)  ! <- Element->Vertex connectivity index
    integer,                   pointer, intent(in) :: elt_vtx(:)      ! <- Element->Vertex connectivity
    double precision,          pointer, intent(in) :: vtx_coord(:,:)  ! <- Vertex coordinates
    integer(kind=pdm_g_num_s), pointer, intent(in) :: elt_ln_to_gn(:) ! <- Element global numerotation
    integer(kind=pdm_g_num_s), pointer, intent(in) :: vtx_ln_to_gn(:) ! <- Vertex global numerotation
    double precision,          pointer, intent(in) :: vtx_field(:)    ! <- Field at local vertices
    double precision,          pointer, intent(in) :: elt_field(:)    ! <- Field at local elements

    type(c_ptr)               :: wrt
    integer                   :: id_geom, id_block
    integer                   :: id_var_elt_part
    integer                   :: id_var_vtx_field
    integer                   :: id_var_elt_field
    double precision, pointer :: val_elt_part(:) => null()
    integer                   :: i_rank, err

    call mpi_comm_rank(comm, i_rank, err)

    call PDM_writer_create(wrt,                    & ! -> Writer instance
                           "Ensight",              & ! <- Format
                           PDM_WRITER_FMT_BIN,     & ! <- Binary files
                           PDM_WRITER_TOPO_CST,    & ! <- Topology is constant over time
                           PDM_WRITER_OFF,         & ! <- Write from scratch
                           "isosurface_nodal_f",   & ! <- Output directory
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

    call pdm_writer_geom_create(wrt,     & ! <- Writer instance
                                id_geom, & ! -> ID of the geometry
                                name,    & ! <- Name of the geometry
                                1)         ! <- Number of partition


    call pdm_writer_geom_coord_set(wrt,              & ! <- Writer instance
                                   id_geom,          & ! <- ID of the geometry
                                   0,                & ! <- ID of current part
                                   n_vtx,            & ! <- Number of vertices
                                   vtx_coord,        & ! <- Vertices coordiantes
                                   vtx_ln_to_gn,     & ! <- Vertex global numerotation
                                   PDM_OWNERSHIP_USER) ! <- Ownership

    if (dim == 2) then
     call pdm_writer_geom_faces_facesom_add(wrt,         & ! <- Writer instance
                                            id_geom,     & ! <- ID of the geometry
                                            0,           & ! <- ID of current part
                                            n_elt,       & ! <- Number of elements
                                            elt_vtx_idx, & ! <- Element->Vertex connectivity index
                                            null(),      &
                                            elt_vtx,     & ! <- Element->Vertex connectivity
                                            elt_ln_to_gn)  ! <- Element global numerotation
    else
      call PDM_writer_geom_bloc_add(wrt,                   & ! <- Writer instance
                                    id_geom,               & ! <- ID of the geometry
                                    PDM_MESH_NODAL_TETRA4, & ! <- Type of element
                                    PDM_OWNERSHIP_USER,    & ! <- Ownership
                                    id_block)                ! -> ID of the block

      call PDM_writer_geom_bloc_std_set(wrt,        & ! <- Writer instance
                                        id_geom,    & ! <- ID of the geometry
                                        id_block,   & ! <- ID of the block
                                        0,          & ! <- ID of current part
                                        n_elt,      & ! <- Number of elements
                                        elt_vtx,    & ! <- Element->Vertex connectivity
                                        elt_ln_to_gn) ! <- Element global numerotation
    end if


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
                               id_geom)   ! <- ID of the geometry

    ! Write 'i_part' variable
    allocate(val_elt_part(n_elt))
    val_elt_part(:) = i_rank

    call PDM_writer_var_set(wrt,             & ! <- Writer instance
                            id_var_elt_part, & ! <- Variable identifier
                            id_geom,         & ! <- ID of the geometry
                            0,               & ! <- ID of current part
                            val_elt_part)      ! <- Variable values

    call PDM_writer_var_write(wrt,             & ! <- Writer instance
                              id_var_elt_part)   ! <- Variable identifier
    deallocate(val_elt_part)


    if (associated(vtx_field)) then
      ! Write 'vtx_field'
      call PDM_writer_var_set(wrt,              & ! <- Writer instance
                              id_var_vtx_field, & ! <- Variable identifier
                              id_geom,          & ! <- ID of the geometry
                              0,                & ! <- ID of current part
                              vtx_field)          ! <- Variable values

      call PDM_writer_var_write(wrt,              & ! <- Writer instance
                                id_var_vtx_field)   ! <- Variable identifier
    endif


    if (associated(elt_field)) then
      ! Write 'elt_field'

      call PDM_writer_var_set(wrt,              & ! <- Writer instance
                              id_var_elt_field, & ! <- Variable identifier
                              id_geom,          & ! <- ID of the geometry
                              0,                & ! <- ID of current part
                              elt_field)          ! <- Variable values

      call PDM_writer_var_write(wrt,              & ! <- Writer instance
                                id_var_elt_field)   ! <- Variable identifier
    endif

    ! End time-step
    call PDM_writer_step_end(wrt)

    ! Free Writer object
    call PDM_writer_free(wrt)




  end subroutine visu_mesh


  !--------------------------------------------------------------------------------
  ! Transfer fields from initial mesh to isosurface mesh (for REEQUILIBRATE mode)
  subroutine data_transfer(isos,          &
                           entity_type,   &
                           field,         &
                           iso_n_elt,     &
                           iso_field)
    implicit none

    type(c_ptr),               intent(in)    :: isos             ! <- Isosurface instance
    integer,                   intent(in)    :: entity_type      ! <- Type of mesh entity
    double precision, pointer, intent(in)    :: field(:)         ! <- Field on local, initial mesh entities
    integer,                   intent(in)    :: iso_n_elt        ! <- Local number of isosurface entities
    double precision, pointer, intent(out)   :: iso_field(:)     ! -> Field on local, isosurface mesh entities

    type(c_ptr)                              :: ptp
    type(PDM_pointer_array_t), pointer       :: pa_field
    type(PDM_pointer_array_t), pointer       :: pa_iso_field
    double precision,          pointer       :: tmp_iso_field(:)
    integer                                  :: request
    integer(pdm_l_num_s),      pointer       :: pelt_parent_idx(:)    => null()
    double precision,          pointer       :: pelt_parent_weight(:) => null()
    double precision,          pointer       :: parent_field(:)
    integer                                  :: i_elt, i_elt_parent


    ! Initialize to null pointers
    ptp           =  C_NULL_PTR
    pa_field      => null()
    pa_iso_field  => null()
    tmp_iso_field => null()

    ! Get PartToPart instance
    call PDM_isosurface_part_to_part_get(isos,               & ! <- ExtractPart instance
                                         0,                  & ! <- ID of current part
                                         entity_type,        & ! <- Type of mesh entity
                                         ptp,                & ! -> PartToPart instance
                                         PDM_OWNERSHIP_KEEP)   ! <- 'isos' keeps ownership of 'ptp'

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
                                        pa_iso_field,                          & ! -> Extracted field
                                        request)                                 ! -> MPI request

    ! Finalize exchange
    call PDM_part_to_part_reverse_iexch_wait(ptp,     & ! <- PartToPart instance
                                             request)   ! <- MPI request

    call PDM_pointer_array_part_get (pa_iso_field, &
                                     0,        &
                                     parent_field)

    call  PDM_isosurface_pparent_weight_get(isos,               & ! <- IsoSurface instance
                                            0,                  & ! <- ID of the isosurface
                                            0,                  & ! <- ID of current part
                                            entity_type,        & ! <- Entity mesh of isosurface : vertices
                                            iso_n_elt,          & ! -> Number of extracted element
                                            pelt_parent_idx,    & ! -> Index of connectivity element->parent between isosurface and global mesh
                                            pelt_parent_weight, & ! -> Interpolation weights 
                                            PDM_OWNERSHIP_KEEP)   ! <- Ownership

    do i_elt=1, iso_n_elt
      iso_field(i_elt) = 0.0d0
      do i_elt_parent=pelt_parent_idx(i_elt), pelt_parent_idx(i_elt+1)-1
        iso_field(i_elt) = iso_field(i_elt) + pelt_parent_weight(i_elt_parent+1) * parent_field(i_elt_parent+1)
      end do
    enddo



    ! Free pointer arrays
    call PDM_pointer_array_free(pa_field)
    call PDM_pointer_array_free(pa_iso_field)

  end subroutine data_transfer
  !--------------------------------------------------------------------------------

end program isosurface_nodal_f