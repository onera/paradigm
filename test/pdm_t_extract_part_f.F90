#include "pdm_configf.h"

program extract_part_f

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use pdm_generate_mesh
  use pdm_extract_part
  use pdm_part_to_part
  use pdm_io
  use pdm_writer
  use iso_c_binding

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !----------------------------------------
  ! Define derived type to store part data
  type part_t

    ! Number of entities
    integer                       :: n_cell = 0
    integer                       :: n_face = 0
    integer                       :: n_edge = 0
    integer                       :: n_vtx  = 0

    ! Connectivities
    integer(pdm_l_num_s), pointer :: cell_face_idx(:) => null()
    integer(pdm_l_num_s), pointer :: cell_face(:)     => null()
    integer(pdm_l_num_s), pointer :: face_vtx_idx(:)  => null()
    integer(pdm_l_num_s), pointer :: face_vtx(:)      => null()
    integer(pdm_l_num_s), pointer :: face_edge_idx(:) => null()
    integer(pdm_l_num_s), pointer :: face_edge(:)     => null()
    integer(pdm_l_num_s), pointer :: edge_vtx(:)      => null()

    ! Coordinates
    double precision,     pointer :: vtx_coord(:,:)   => null()
    double precision,     pointer :: cell_center(:,:) => null()

    ! Global IDs
    integer(pdm_g_num_s), pointer :: cell_ln_to_gn(:) => null()
    integer(pdm_g_num_s), pointer :: face_ln_to_gn(:) => null()
    integer(pdm_g_num_s), pointer :: edge_ln_to_gn(:) => null()
    integer(pdm_g_num_s), pointer :: vtx_ln_to_gn(:)  => null()

    ! Parent global IDs
    integer(pdm_g_num_s), pointer :: cell_parent_ln_to_gn(:) => null()
    integer(pdm_g_num_s), pointer :: face_parent_ln_to_gn(:) => null()
    integer(pdm_g_num_s), pointer :: vtx_parent_ln_to_gn(:)  => null()

    ! Fields
    double precision,     pointer :: vtx_field(:)  => null()
    double precision,     pointer :: cell_field(:) => null()

    ! Selected cells
    integer                       :: n_selected  = 0
    integer(pdm_l_num_s), pointer :: selected(:) => null()

  end type part_t
  !----------------------------------------


  !---------------------------------------------------------------
  integer, parameter            :: comm = MPI_COMM_WORLD           ! MPI Communicator
  integer                       :: i_rank, ierr

  integer                       :: i_arg
  character(len=99)             :: arg

  integer                       :: extract_kind                    ! Extraction kind (LOCAL or REEQUILIBRATE)
  integer                       :: n_part_in                       ! Number of initial parts per MPI rank
  integer                       :: n_part_out                      ! Number of output  parts per MPI rank
  integer                       :: part_method                     ! Re-partitioning method kind (if REEQUILIBRATE)
  logical                       :: visu                            ! Enable output for visualization
  integer(pdm_g_num_s)          :: n_subdiv                        ! Number of subdivisions (mesh density)

  type(part_t), allocatable     :: ini_parts(:)                    ! Initial parts
  type(part_t), allocatable     :: ext_parts(:)                    ! Extracted parts

  type(c_ptr)                   :: extrp = C_NULL_PTR              ! ExtractPart (ParaDiGM object)

  integer                       :: i_part, i_cell, i_face, i_vtx
  integer                       :: idx_face, idx_vtx, i
  double precision              :: min_x, max_x

  integer                       :: extract_n_vtx
  integer(pdm_l_num_s), pointer :: extract_vtx_parent(:) => null() ! Extracted vertices IDs
  !---------------------------------------------------------------

  extract_kind = PDM_EXTRACT_PART_KIND_LOCAL
  n_part_in    = 1
  n_part_out   = 1
  part_method  = PDM_SPLIT_DUAL_WITH_HILBERT
  visu         = .false.
  n_subdiv     = 10

  !----------------------------------------
  ! Parse command line arguments
  i_arg = 1
  do while (i_arg <= command_argument_count())
    call get_command_argument(i_arg, arg)

    select case(arg)

      case ("-reequilibrate")
        extract_kind = PDM_EXTRACT_PART_KIND_REEQUILIBRATE

      case ("-n_part_in")
        i_arg = i_arg + 1
        call get_command_argument(i_arg, arg)
        read(arg, *) n_part_in

      case ("-n_part_out")
        i_arg = i_arg + 1
        call get_command_argument(i_arg, arg)
        read(arg, *) n_part_out

      case ("-part_method")
        i_arg = i_arg + 1
        call get_command_argument(i_arg, arg)
        read(arg, *) part_method

      case ("-n")
        i_arg = i_arg + 1
        call get_command_argument(i_arg, arg)
        read(arg, *) n_subdiv

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
  ! Generate a simple 3D mesh
  call generate_mesh(comm,        & ! <- MPI communicator
                     n_subdiv,    & ! <- Number of subdivisions (mesh density)
                     n_part_in,   & ! <- Number of parts
                     part_method, & ! <- Partitioning method
                     ini_parts)     ! -> Mesh parts
  !----------------------------------------


  !----------------------------------------
  ! Create ExtractPart instance
  call PDM_extract_part_create(extrp,              & ! -> ExtractPart instance
                               3,                  & ! <- Dimension
                               n_part_in,          & ! <- Number of parts per rank (input)
                               n_part_out,         & ! <- Number of parts per rank (output)
                               extract_kind,       & ! <- Extraction kind
                               part_method,        & ! <- Re-partitioning method (not used in LOCAL mode)
                               .true.,             & ! <- Generate global IDs restricted to extraction
                               PDM_OWNERSHIP_KEEP, & ! <- Ownership of extraction
                               comm)                 ! <- MPI communicator


  ! Set input mesh
  do i_part = 1, n_part_in
    call pdm_extract_part_part_set(extrp,                           & ! <- ExtractPart instance
                                   i_part-1,                        & ! <- ID of current part
                                   ini_parts(i_part)%n_cell,        & ! <- Number of cells
                                   ini_parts(i_part)%n_face,        & ! <- Number of faces
                                   ini_parts(i_part)%n_edge,        & ! <- Number of edges
                                   ini_parts(i_part)%n_vtx,         & ! <- Number of vertices
                                   ini_parts(i_part)%cell_face_idx, & ! <- Index for cell->face connectivity
                                   ini_parts(i_part)%cell_face,     & ! <- Cell->face connectivity
                                   ini_parts(i_part)%face_edge_idx, & ! <- Index for face->edge connectivity
                                   ini_parts(i_part)%face_edge,     & ! <- Face->edge connectivity
                                   ini_parts(i_part)%edge_vtx,      & ! <- Edge->vtx connectivity
                                   ini_parts(i_part)%face_vtx_idx,  & ! <- Index for face->vtx connectivity
                                   ini_parts(i_part)%face_vtx,      & ! <- Face->vtx connectivity
                                   ini_parts(i_part)%cell_ln_to_gn, & ! <- Cell global IDs
                                   ini_parts(i_part)%face_ln_to_gn, & ! <- Face global IDs
                                   ini_parts(i_part)%edge_ln_to_gn, & ! <- Edge global IDs
                                   ini_parts(i_part)%vtx_ln_to_gn,  & ! <- Vertex global IDs
                                   ini_parts(i_part)%vtx_coord)       ! <- Vertex coordinates
  enddo

  ! Select cells to extract : the ones that cross the (x = 0) plane
  do i_part = 1, n_part_in
    ini_parts(i_part)%n_selected = 0
    allocate(ini_parts(i_part)%selected(ini_parts(i_part)%n_cell))

    do i_cell = 1, ini_parts(i_part)%n_cell
      min_x =  1.d30
      max_x = -1.d30
      do idx_face = ini_parts(i_part)%cell_face_idx(i_cell)+1, ini_parts(i_part)%cell_face_idx(i_cell+1)
        i_face = abs(ini_parts(i_part)%cell_face(idx_face))

        do idx_vtx = ini_parts(i_part)%face_vtx_idx(i_face)+1, ini_parts(i_part)%face_vtx_idx(i_face+1)
          i_vtx = ini_parts(i_part)%face_vtx(idx_vtx)

          min_x = min(min_x, ini_parts(i_part)%vtx_coord(1,i_vtx))
          max_x = max(max_x, ini_parts(i_part)%vtx_coord(1,i_vtx))
        enddo ! End loop on current face's vertices
      enddo ! End loop on current cell's faces

      if (min_x <= 0.d0 .and. max_x >= 0.d0) then
        ini_parts(i_part)%n_selected = ini_parts(i_part)%n_selected + 1
        ini_parts(i_part)%selected(ini_parts(i_part)%n_selected) = i_cell
      endif
    enddo ! End loop on cells

    call PDM_extract_part_selected_lnum_set(extrp,                        & ! <- ExtractPart instance
                                            i_part-1,                     & ! <- ID of current part
                                            ini_parts(i_part)%n_selected, & ! <- Local number of extracted elements
                                            ini_parts(i_part)%selected)     ! <- Local IDs of extracted elements
  enddo ! End loop on parts

  ! Realize extraction
  call PDM_extract_part_compute(extrp)

  ! Retrieve extracted mesh
  allocate(ext_parts(n_part_out))
  do i_part = 1, n_part_out
    ! Connectivities
    call PDM_extract_part_connectivity_get(extrp,                           & ! <- ExtractPart instance
                                           i_part-1,                        & ! <- ID of current part
                                           PDM_CONNECTIVITY_TYPE_CELL_FACE, & ! <- Cell->face
                                           ext_parts(i_part)%n_cell,        & ! -> Number of extracted cells
                                           ext_parts(i_part)%cell_face,     & ! -> Extracted cell->face connectivity
                                           ext_parts(i_part)%cell_face_idx, & ! -> Index for extracted cell->face connectivity
                                           PDM_OWNERSHIP_USER)                ! <- Ownership

    call PDM_extract_part_connectivity_get(extrp,                           & ! <- ExtractPart instance
                                           i_part-1,                        & ! <- ID of current part
                                           PDM_CONNECTIVITY_TYPE_FACE_VTX,  & ! <- Face->vtx
                                           ext_parts(i_part)%n_face,        & ! -> Number of extracted faces
                                           ext_parts(i_part)%face_vtx,      & ! -> Extracted face->vtx connectivity
                                           ext_parts(i_part)%face_vtx_idx,  & ! -> Index for extracted face->vtx connectivity
                                           PDM_OWNERSHIP_USER)                ! <- Ownership

    ! Coordinates
    call PDM_extract_part_vtx_coord_get(extrp,                       & ! <- ExtractPart instance
                                        i_part-1,                    & ! <- ID of current parts
                                        ext_parts(i_part)%n_vtx,     & ! -> Number of extracted vertices
                                        ext_parts(i_part)%vtx_coord, & ! -> Coordinates of extracted vertices
                                        PDM_OWNERSHIP_USER)            ! <- Ownership

    ! Global IDs (in extracted mesh)
    call PDM_extract_part_ln_to_gn_get(extrp,                           & ! <- ExtractPart instance
                                       i_part-1,                        & ! <- ID of current part
                                       PDM_MESH_ENTITY_CELL,            & ! <- Cells
                                       ext_parts(i_part)%n_cell,        & ! -> Number of extracted cells
                                       ext_parts(i_part)%cell_ln_to_gn, & ! -> Global IDs of extracted cells
                                       PDM_OWNERSHIP_USER)                ! <- Ownership

    call PDM_extract_part_ln_to_gn_get(extrp,                           & ! <- ExtractPart instance
                                       i_part-1,                        & ! <- ID of current part
                                       PDM_MESH_ENTITY_FACE,            & ! <- Faces
                                       ext_parts(i_part)%n_face,        & ! -> Number of extracted faces
                                       ext_parts(i_part)%face_ln_to_gn, & ! -> Global IDs of extracted faces
                                       PDM_OWNERSHIP_USER)                ! <- Ownership

    call PDM_extract_part_ln_to_gn_get(extrp,                           & ! <- ExtractPart instance
                                       i_part-1,                        & ! <- ID of current part
                                       PDM_MESH_ENTITY_VTX,             & ! <- Vertices
                                       ext_parts(i_part)%n_vtx,         & ! -> Number of extracted vertices
                                       ext_parts(i_part)%vtx_ln_to_gn,  & ! -> Global IDs of extracted vertices
                                       PDM_OWNERSHIP_USER)                ! <- Ownership

    ! Global IDs (in initial mesh)
    call PDM_extract_part_parent_ln_to_gn_get(extrp,                                  & ! <- ExtractPart instance
                                              i_part-1,                               & ! <- ID of current part
                                              PDM_MESH_ENTITY_CELL,                   & ! <- Cells
                                              ext_parts(i_part)%n_cell,               & ! -> Number of extracted cells
                                              ext_parts(i_part)%cell_parent_ln_to_gn, & ! -> Parent global IDs of extracted cells
                                              PDM_OWNERSHIP_USER)                       ! <- Ownership

    call PDM_extract_part_parent_ln_to_gn_get(extrp,                                  & ! <- ExtractPart instance
                                              i_part-1,                               & ! <- ID of current part
                                              PDM_MESH_ENTITY_FACE,                   & ! <- Faces
                                              ext_parts(i_part)%n_face,               & ! -> Number of extracted faces
                                              ext_parts(i_part)%face_parent_ln_to_gn, & ! -> Parent global IDs of extracted faces
                                              PDM_OWNERSHIP_USER)                       ! <- Ownership

    call PDM_extract_part_parent_ln_to_gn_get(extrp,                                  & ! <- ExtractPart instance
                                              i_part-1,                               & ! <- ID of current part
                                              PDM_MESH_ENTITY_VTX,                    & ! <- Vertices
                                              ext_parts(i_part)%n_vtx,                & ! -> Number of extracted vertices
                                              ext_parts(i_part)%vtx_parent_ln_to_gn,  & ! -> Parent global IDs of extracted vertices
                                              PDM_OWNERSHIP_USER)                       ! <- Ownership

    if (extract_kind == PDM_EXTRACT_PART_KIND_LOCAL .and. ext_parts(i_part)%n_cell /= ini_parts(i_part)%n_selected) then
      print *, "extracted n_cell =", ext_parts(i_part)%n_cell, " but expected", ini_parts(i_part)%n_selected
      STOP
    endif
  enddo ! End loop on parts
  !----------------------------------------


  !----------------------------------------
  ! Transfer data from initial mesh to extraction
  !  Create dummy field on initial mesh
  do i_part = 1, n_part_in
    allocate(ini_parts(i_part)%vtx_field(ini_parts(i_part)%n_vtx))
    call eval_field(ini_parts(i_part)%n_vtx,     &
                    ini_parts(i_part)%vtx_coord, &
                    ini_parts(i_part)%vtx_field)
  enddo

  do i_part = 1, n_part_in
    allocate(ini_parts(i_part)%cell_field(ini_parts(i_part)%n_cell))
    call eval_field(ini_parts(i_part)%n_cell,      &
                    ini_parts(i_part)%cell_center, &
                    ini_parts(i_part)%cell_field)
  enddo

  !  Transfer field at vertices
  do i_part = 1, n_part_out
    allocate(ext_parts(i_part)%vtx_field(ext_parts(i_part)%n_vtx))
  enddo

  if (extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) then
    ! Local transfer
    do i_part = 1, n_part_out
      call PDM_extract_part_parent_lnum_get(extrp,               & ! <- ExtractPart instance
                                            i_part-1,            & ! <- ID of current part
                                            PDM_MESH_ENTITY_VTX, & ! <- Vertices
                                            extract_n_vtx,       & ! -> Number of extracted vertices in current subdomain
                                            extract_vtx_parent,  & ! -> Local IDs of extracted vertices in initial mesh
                                            PDM_OWNERSHIP_KEEP)    ! <- Ownership ('extrp' keeps ownership)
      do i = 1, ext_parts(i_part)%n_vtx
        ext_parts(i_part)%vtx_field(i) = ini_parts(i_part)%vtx_field(extract_vtx_parent(i))
      enddo
    enddo
  else
    ! Reequilibrate mode => transfer in parallel
    call data_transfer(extrp,               & ! <-  ExtractPart instance
                       PDM_MESH_ENTITY_VTX, & ! <-  Vertices
                       n_part_in,           & ! <-  Number of initial parts
                       ini_parts,           & ! <-  Initial parts
                       n_part_out,          & ! <-  Number of extracted parts
                       ext_parts)             ! <-> Extracted parts
  endif


  !  Transfer field at cells
  do i_part = 1, n_part_out
    allocate(ext_parts(i_part)%cell_field(ext_parts(i_part)%n_cell))
  enddo

  if (extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) then
    ! Local transfer
    do i_part = 1, n_part_out
      do i = 1, ext_parts(i_part)%n_cell
        ext_parts(i_part)%cell_field(i) = ini_parts(i_part)%cell_field(ini_parts(i_part)%selected(i))
      enddo
    enddo
  else
    ! Reequilibrate mode => transfer in parallel
    call data_transfer(extrp,                & ! <-  ExtractPart instance
                       PDM_MESH_ENTITY_CELL, & ! <-  Cells
                       n_part_in,            & ! <-  Number of initial parts
                       ini_parts,            & ! <-  Initial parts
                       n_part_out,           & ! <-  Number of extracted parts
                       ext_parts)              ! <-> Extracted parts
  endif
  !----------------------------------------


  !----------------------------------------
  ! Visu
  if (visu) then
    call visu_parts(comm,      & ! <- MPI communicator
                    n_part_in, & ! <- Number of parts
                    ini_parts, & ! <- Mesh parts
                    "init")      ! <- Output name

    call visu_parts(comm,       & ! <- MPI communicator
                    n_part_out, & ! <- Number of parts
                    ext_parts,  & ! <- Mesh parts
                    "extract")    ! <- Output name
  endif
  !----------------------------------------


  !----------------------------------------
  ! Free memory
  call PDM_extract_part_free(extrp)
  do i_part = 1, n_part_in
    call free_part(ini_parts(i_part))
  enddo
  deallocate(ini_parts)

  do i_part = 1, n_part_out
    call free_part(ext_parts(i_part))
  enddo
  deallocate(ext_parts)
  !----------------------------------------

  if (i_rank == 0) then
    print *, "The End :)"
  endif

  call mpi_finalize(ierr)


contains

  !--------------------------------------------------------------------------------
  ! Auxiliary routines to make code lighter and easier to read
  subroutine data_transfer(extrp,       &
                           entity_type, &
                           n_part_in,   &
                           ini_parts,   &
                           n_part_out,  &
                           ext_parts)
    ! Transfer fields from initial mesh to extraction (for REEQUILIBRATE mode)
    implicit none

    type(c_ptr),         intent(in)    :: extrp        ! <-  ExtractPart instance
    integer,             intent(in)    :: entity_type  ! <-  Type of mesh entity
    integer,             intent(in)    :: n_part_in    ! <-  Number of initial parts
    type(part_t),        intent(in)    :: ini_parts(:) ! <-  Initial parts
    integer,             intent(in)    :: n_part_out   ! <-  Number of extracted parts
    type(part_t),        intent(inout) :: ext_parts(:) ! <-> Extracted parts

    type(c_ptr)                        :: ptp
    type(PDM_pointer_array_t), pointer :: pa_field
    type(PDM_pointer_array_t), pointer :: pa_extract_field
    double precision,          pointer :: tmp_extract_field(:)
    integer                            :: request

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
                                  n_part_in,       & ! <- Number of parts on current rank
                                  PDM_TYPE_DOUBLE)   ! <- Data type

    do i_part = 1, n_part_in
      if (entity_type == PDM_MESH_ENTITY_VTX) then
        call PDM_pointer_array_part_set(pa_field,                     & ! <- Pointer array
                                        i_part-1,                     & ! <- ID of current part
                                        ini_parts(i_part)%vtx_field)    ! <- Field
      else ! PDM_MESH_ENTITY_CELL
        call PDM_pointer_array_part_set(pa_field,                     & ! <- Pointer array
                                        i_part-1,                     & ! <- ID of current part
                                        ini_parts(i_part)%cell_field)   ! <- Field
      endif
    enddo

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

    do i_part = 1, n_part_out
      call PDM_pointer_array_part_get(pa_extract_field,   & ! <- Pointer array
                                      i_part-1,           & ! <- ID of current part
                                      tmp_extract_field)    ! <- Extracted field

      ! Copy data to avoid ownership conflicts
      if (entity_type == PDM_MESH_ENTITY_VTX) then
        ext_parts(i_part)%vtx_field(:)  = tmp_extract_field(:)
      else ! PDM_MESH_ENTITY_CELL
        ext_parts(i_part)%cell_field(:) = tmp_extract_field(:)
      endif
    enddo

    ! Free pointer arrays
    call PDM_pointer_array_free(pa_field)
    call PDM_pointer_array_free(pa_extract_field)

  end subroutine data_transfer



  subroutine free_part(part)
    ! Deallocate a `part_t` instance
    implicit none

    type(part_t), intent(inout) :: part ! <-> Mesh part

    call PDM_fortran_free_c(c_loc(part%cell_face_idx))
    call PDM_fortran_free_c(c_loc(part%cell_face))
    call PDM_fortran_free_c(c_loc(part%face_vtx_idx))
    call PDM_fortran_free_c(c_loc(part%face_vtx))
    call PDM_fortran_free_c(c_loc(part%vtx_coord))
    call PDM_fortran_free_c(c_loc(part%cell_ln_to_gn))
    call PDM_fortran_free_c(c_loc(part%face_ln_to_gn))
    call PDM_fortran_free_c(c_loc(part%vtx_ln_to_gn))
    if (associated(part%cell_parent_ln_to_gn)) call PDM_fortran_free_c(c_loc(part%cell_parent_ln_to_gn))
    if (associated(part%face_parent_ln_to_gn)) call PDM_fortran_free_c(c_loc(part%face_parent_ln_to_gn))
    if (associated(part%vtx_parent_ln_to_gn))  call PDM_fortran_free_c(c_loc(part%vtx_parent_ln_to_gn))
    if (associated(part%selected))             deallocate(part%selected)
    if (associated(part%cell_center))          deallocate(part%cell_center)
    deallocate(part%vtx_field)
    deallocate(part%cell_field)

  end subroutine free_part



  subroutine visu_parts(comm,   &
                        n_part, &
                        parts,  &
                        name)
    ! Export mesh parts in Ensight format
    implicit none

    integer,                   intent(in) :: comm     ! <- MPI communicator
    integer,                   intent(in) :: n_part   ! <- Number of parts per MPI rank
    type(part_t),              intent(in) :: parts(:) ! <- Mesh parts
    character(len = *),        intent(in) :: name     ! <- Output file name

    type(c_ptr)                           :: wrt
    integer                               :: id_geom
    integer                               :: id_var_elt_part
    integer                               :: id_var_vtx_field
    integer                               :: id_var_elt_field
    double precision, pointer             :: val_elt_part(:)

    integer                               :: i_rank, err


    call mpi_comm_rank(comm, i_rank, err)

    call PDM_writer_create(wrt,                    & ! -> Writer instance
                           "Ensight",              & ! <- Format
                           PDM_WRITER_FMT_BIN,     & ! <- Binary files
                           PDM_WRITER_TOPO_CST,    & ! <- Topology is constant over time
                           PDM_WRITER_OFF,         & ! <- Write from scratch
                           "extract_part_f",       & ! <- Output directory
                           name,                   & ! <- Output file
                           comm,                   & ! <- MPI communicator
                           PDM_IO_KIND_MPI_SIMPLE, & ! <- MPI IO strategy (type of file access)
                           1.d0,                   & ! <- MPI IO strategy (proportion of active nodes)
                           "")                       ! <- No additional options

    ! Create a geometry
    call PDM_writer_geom_create(wrt,     & ! <- Writer instance
                                id_geom, & ! -> Geometry identifier
                                name,    & ! <- Geometry name
                                n_part)    ! <- Number of parts

    ! Set geometry
    do i_part = 1, n_part
      call PDM_writer_geom_coord_set(wrt,                        & ! <- Writer instance
                                     id_geom,                    & ! <- Geometry identifier
                                     i_part-1,                   & ! <- ID of current part
                                     parts(i_part)%n_vtx,        & ! <- Number of vertices
                                     parts(i_part)%vtx_coord,    & ! <- Vertex coordinates
                                     parts(i_part)%vtx_ln_to_gn, & ! <- Vertex global IDs
                                     PDM_OWNERSHIP_USER)           ! <- Ownership

      call PDM_writer_geom_cell3d_cellface_add(wrt,                         & ! <- Writer instance
                                               id_geom,                     & ! <- Geometry identifier
                                               i_part-1,                    & ! <- ID of current part
                                               parts(i_part)%n_cell,        & ! <- Number of cells
                                               parts(i_part)%n_face,        & ! <- Number of faces
                                               parts(i_part)%face_vtx_idx,  & ! <- Index for face->vtx connectivity
                                               null(),                      & ! <- not used
                                               parts(i_part)%face_vtx,      & ! <- Face->vtx connectivity
                                               parts(i_part)%cell_face_idx, & ! <- Index for cell->face connectivity
                                               null(),                      & ! <- not used
                                               parts(i_part)%cell_face,     & ! <- Cell->face connectivity
                                               parts(i_part)%cell_ln_to_gn)   ! <- Cell global IDs
    enddo


    ! Define a variable to visualize the partitioning
    call PDM_writer_var_create(wrt,                     & ! <- Writer instance
                               id_var_elt_part,         & ! -> Variable identifier
                               PDM_WRITER_OFF,          & ! <- Not time-dependent
                               PDM_WRITER_VAR_SCALAIRE, & ! <- Scalar variable
                               PDM_WRITER_VAR_ELEMENTS, & ! <- Element-based variable
                               "i_part")                  ! <- Variable name

    ! Define a variable to visualize the field at vertices
    call PDM_writer_var_create(wrt,                     & ! <- Writer instance
                               id_var_vtx_field,        & ! -> Variable identifier
                               PDM_WRITER_OFF,          & ! <- Not time-dependent
                               PDM_WRITER_VAR_SCALAIRE, & ! <- Scalar variable
                               PDM_WRITER_VAR_VERTICES, & ! <- Vertex-based variable
                               "vtx_field")               ! <- Variable name

    ! Define a variable to visualize the field at cells
    call PDM_writer_var_create(wrt,                      & ! <- Writer instance
                               id_var_elt_field,         & ! -> Variable identifier
                               PDM_WRITER_OFF,           & ! <- Not time-dependent
                               PDM_WRITER_VAR_SCALAIRE,  & ! <- Scalar variable
                               PDM_WRITER_VAR_ELEMENTS,  & ! <- Vertex-based variable
                               "elt_field")                ! <- Variable name


    ! Begin a time-step (only one here)
    call PDM_writer_step_beg(wrt,  & ! <- Writer instance
                             0.d0)   ! <- Time

    ! Write the geometry
    call PDM_writer_geom_write(wrt,     & ! <- Writer instance
                               id_geom)   ! <- Geometry identifier

    ! Set variables
    do i_part = 1, n_part
      allocate(val_elt_part(parts(i_part)%n_cell))
      val_elt_part(:) = n_part*i_rank + i_part

      call PDM_writer_var_set(wrt,             & ! <- Writer instance
                              id_var_elt_part, & ! <- Variable identifier
                              id_geom,         & ! <- Geometry identifier
                              i_part-1,        & ! <- ID of current part
                              val_elt_part)      ! <- Variable values
      deallocate(val_elt_part)

      call PDM_writer_var_set(wrt,                     & ! <- Writer instance
                              id_var_vtx_field,        & ! <- Variable identifier
                              id_geom,                 & ! <- Geometry identifier
                              i_part-1,                & ! <- ID of current part
                              parts(i_part)%vtx_field)   ! <- Variable values

      call PDM_writer_var_set(wrt,                      & ! <- Writer instance
                              id_var_elt_field,         & ! <- Variable identifier
                              id_geom,                  & ! <- Geometry identifier
                              i_part-1,                 & ! <- ID of current part
                              parts(i_part)%cell_field)   ! <- Variable values
    enddo

    ! Write variables
    call PDM_writer_var_write(wrt,              & ! <- Writer instance
                              id_var_elt_part)    ! <- Variable identifier

    call PDM_writer_var_write(wrt,              & ! <- Writer instance
                              id_var_vtx_field)   ! <- Variable identifier

    call PDM_writer_var_write(wrt,              & ! <- Writer instance
                              id_var_elt_field)   ! <- Variable identifier

    ! End time-step
    call PDM_writer_step_end(wrt)

    ! Free Writer object
    call PDM_writer_free(wrt)

  end subroutine visu_parts



  subroutine eval_field(n_pts, &
                        coord, &
                        field)
    ! Evaluate a dummy field at a set of points
    implicit none

    integer,                   intent(in)    :: n_pts
    double precision, pointer, intent(in)    :: coord(:,:)
    double precision, pointer, intent(inout) :: field(:)

    field(1:n_pts) = cos(4*(coord(1,1:n_pts) - coord(2,1:n_pts) - coord(3,1:n_pts)))

  end subroutine eval_field



  subroutine generate_mesh(comm,        &
                           n_subdiv,    &
                           n_part,      &
                           part_method, &
                           parts)
    ! Generate a simple 3D mesh (ball)
    implicit none

    integer,                   intent(in ) :: comm        ! <- MPI Communicator
    integer(pdm_g_num_s),      intent(in ) :: n_subdiv    ! <- Number of subdivisions
    integer,                   intent(in ) :: n_part      ! <- Number of parts per MPI rank
    integer,                   intent(in ) :: part_method ! <- Partitioning method
    type(part_t), allocatable, intent(out) :: parts(:)    ! -> Mesh parts

    integer(pdm_g_num_s)                   :: n_layer = 0
    integer(pdm_l_num_s),      pointer     :: pn_vtx(:)
    integer(pdm_l_num_s),      pointer     :: pn_edge(:)
    integer(pdm_l_num_s),      pointer     :: pn_face(:)
    integer(pdm_l_num_s),      pointer     :: pn_cell(:)
    integer(pdm_l_num_s),      pointer     :: pn_surface(:)
    type(PDM_pointer_array_t), pointer     :: pvtx_coord             => null()
    type(PDM_pointer_array_t), pointer     :: pedge_vtx              => null()
    type(PDM_pointer_array_t), pointer     :: pface_edge_idx         => null()
    type(PDM_pointer_array_t), pointer     :: pface_edge             => null()
    type(PDM_pointer_array_t), pointer     :: pface_vtx              => null()
    type(PDM_pointer_array_t), pointer     :: pcell_face_idx         => null()
    type(PDM_pointer_array_t), pointer     :: pcell_face             => null()
    type(PDM_pointer_array_t), pointer     :: pvtx_ln_to_gn          => null()
    type(PDM_pointer_array_t), pointer     :: pedge_ln_to_gn         => null()
    type(PDM_pointer_array_t), pointer     :: pface_ln_to_gn         => null()
    type(PDM_pointer_array_t), pointer     :: pcell_ln_to_gn         => null()
    type(PDM_pointer_array_t), pointer     :: psurface_face_idx      => null()
    type(PDM_pointer_array_t), pointer     :: psurface_face          => null()
    type(PDM_pointer_array_t), pointer     :: psurface_face_ln_to_gn => null()

    double precision,          pointer     :: face_center(:,:)
    integer                                :: i_cell, i_face, i_vtx
    integer                                :: idx_face, idx_vtx
    integer                                :: n_connect

    call PDM_generate_mesh_ball_ngon(comm,                   &
                                     PDM_MESH_NODAL_TETRA4,  &
                                     1,                      &
                                     C_NULL_PTR,             &
                                     1.d0,                   &
                                     0.d0,                   &
                                     0.d0,                   &
                                     0.d0,                   &
                                     0.d0,                   &
                                     n_subdiv,               &
                                     n_subdiv,               &
                                     n_subdiv,               &
                                     n_layer,                &
                                     1.d0,                   &
                                     n_part,                 &
                                     part_method,            &
                                     pn_vtx,                 &
                                     pn_edge,                &
                                     pn_face,                &
                                     pn_cell,                &
                                     pvtx_coord,             &
                                     pedge_vtx,              &
                                     pface_edge_idx,         &
                                     pface_edge,             &
                                     pface_vtx,              &
                                     pcell_face_idx,         &
                                     pcell_face,             &
                                     pvtx_ln_to_gn,          &
                                     pedge_ln_to_gn,         &
                                     pface_ln_to_gn,         &
                                     pcell_ln_to_gn,         &
                                     pn_surface,             &
                                     psurface_face_idx,      &
                                     psurface_face,          &
                                     psurface_face_ln_to_gn)

    allocate(parts(n_part))
    do i_part = 1, n_part

      parts(i_part)%n_cell = pn_cell(i_part)
      parts(i_part)%n_face = pn_face(i_part)
      parts(i_part)%n_vtx  = pn_vtx (i_part)

      call PDM_pointer_array_part_get(pcell_face_idx,              &
                                      i_part-1,                    &
                                      parts(i_part)%cell_face_idx)
      call PDM_pointer_array_part_get(pcell_face,                  &
                                      i_part-1,                    &
                                      parts(i_part)%cell_face)

      call PDM_pointer_array_part_get(pface_edge_idx,              &
                                      i_part-1,                    &
                                      parts(i_part)%face_vtx_idx)
      call PDM_pointer_array_part_get(pface_vtx,                   &
                                      i_part-1,                    &
                                      parts(i_part)%face_vtx)

      call PDM_pointer_array_part_get(pvtx_coord,                  &
                                      i_part-1,                    &
                                      PDM_STRIDE_CST_INTERLACED,   &
                                      3,                           &
                                      parts(i_part)%vtx_coord)

      call PDM_pointer_array_part_get(pcell_ln_to_gn,              &
                                      i_part-1,                    &
                                      parts(i_part)%cell_ln_to_gn)
      call PDM_pointer_array_part_get(pface_ln_to_gn,              &
                                      i_part-1,                    &
                                      parts(i_part)%face_ln_to_gn)
      call PDM_pointer_array_part_get(pvtx_ln_to_gn,               &
                                      i_part-1,                    &
                                      parts(i_part)%vtx_ln_to_gn)


      allocate(face_center(3, parts(i_part)%n_face))
      do i_face = 1, parts(i_part)%n_face
        face_center(1:3,i_face) = 0.d0
        do idx_vtx = parts(i_part)%face_vtx_idx(i_face)+1, parts(i_part)%face_vtx_idx(i_face+1)
          i_vtx = parts(i_part)%face_vtx(idx_vtx)
          face_center(1:3,i_face) = face_center(1:3,i_face) + parts(i_part)%vtx_coord(1:3,i_vtx)
        enddo
        n_connect = (parts(i_part)%face_vtx_idx(i_face+1) - ini_parts(i_part)%face_vtx_idx(i_face))
        face_center(1:3,i_face) = face_center(1:3,i_face) / dble(n_connect)
      enddo


      allocate(parts(i_part)%cell_center(3, parts(i_part)%n_cell))
      do i_cell = 1, parts(i_part)%n_cell
        parts(i_part)%cell_center(1:3,i_cell) = 0.d0
        do idx_face = parts(i_part)%cell_face_idx(i_cell)+1, parts(i_part)%cell_face_idx(i_cell+1)
          i_face = abs(parts(i_part)%cell_face(idx_face))
          parts(i_part)%cell_center(1:3,i_cell) = parts(i_part)%cell_center(1:3,i_cell) + face_center(1:3,i_face)
        enddo
        n_connect = (parts(i_part)%cell_face_idx(i_cell+1) - parts(i_part)%cell_face_idx(i_cell))
        parts(i_part)%cell_center(1:3,i_cell) = parts(i_part)%cell_center(1:3,i_cell) / dble(n_connect)

      enddo
      deallocate(face_center)
    enddo

    call PDM_fortran_free_c(c_loc(pn_cell))
    call PDM_fortran_free_c(c_loc(pn_face))
    call PDM_fortran_free_c(c_loc(pn_edge))
    call PDM_fortran_free_c(c_loc(pn_vtx))
    call PDM_fortran_free_c(c_loc(pn_surface))

    call PDM_pointer_array_free(pface_edge)
    call PDM_pointer_array_free(pedge_vtx)
    call PDM_pointer_array_free(pedge_ln_to_gn)

    call PDM_pointer_array_free(psurface_face_idx)
    call PDM_pointer_array_free(psurface_face)
    call PDM_pointer_array_free(psurface_face_ln_to_gn)

  end subroutine generate_mesh

end program extract_part_f
