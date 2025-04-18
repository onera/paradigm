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
  use pdm_io
  use pdm_writer
  use iso_c_binding

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !---------------------------------------------------------------
  integer, parameter            :: comm = MPI_COMM_WORLD         ! MPI Communicator
  integer                       :: i_rank, ierr

  integer                       :: dim                           ! Mesh dimension
  integer                       :: n_part                        ! Number of partitions on current MPI rank

  integer                       :: n_vtx                         ! Number of vertices
  double precision,     pointer :: vtx_coord(:,:)  => null()     ! Vertex coordinates
  integer(pdm_g_num_s), pointer :: vtx_ln_to_gn(:) => null()     ! Vertex global IDs
  ! integer                       :: n_tri                         ! Number of triangles
  ! integer(pdm_g_num_s), pointer :: tri_vtx(:)      => null()     ! Triangle->vtx connectivity
  ! integer(pdm_g_num_s), pointer :: tri_ln_to_gn(:) => null()     ! Triangle global IDs
  integer                       :: n_tet                         ! Number of tetra
  integer(pdm_l_num_s), pointer :: tet_vtx(:)      => null()     ! Tetra->vtx connectivity
  integer(pdm_g_num_s), pointer :: tet_ln_to_gn(:) => null()     ! Tetra global IDs
  integer(pdm_l_num_s), pointer :: parent_num(:)   => null()     ! Not used
  integer(pdm_g_num_s), pointer :: parent_gnum(:)  => null()     ! Not used
  integer(pdm_l_num_s), pointer :: tet_vtx_idx(:)  => null()     ! Not used

  type(c_ptr)                   :: mesh = C_NULL_PTR             ! Initial PartMeshNodal (ParaDiGM object)
  ! integer                       :: id_tri                        ! ID of triangle section
  integer                       :: id_tet                        ! ID of tetra section

  type(c_ptr)                   :: extrp = C_NULL_PTR            ! ExtractPart (ParaDiGM object)
  integer                       :: n_selected                    ! Number of extracted elements
  integer(pdm_l_num_s), pointer :: selected(:) => null()         ! List of extracted elements

  type(c_ptr)                   :: extracted_mesh = C_NULL_PTR   ! Extracted PartMeshNodal (ParaDiGM object)

  integer                       :: i_tet, i_vtx, i
  double precision              :: min_x, max_x
  !---------------------------------------------------------------


  !----------------------------------------
  ! Initialize MPI
  call mpi_init(ierr)
  call mpi_comm_rank(comm, i_rank, ierr)
  !----------------------------------------


  !----------------------------------------
  ! Generate mesh
  call PDM_generate_mesh_ball_simplified(comm,        &
                                         n_vtx,       &
                                         n_tet,       &
                                         vtx_coord,   &
                                         tet_vtx_idx, &
                                         tet_vtx)
  call PDM_fortran_free_c(c_loc(tet_vtx_idx))

  call generate_gnums(comm,         &
                      n_vtx,        &
                      n_tet,        &
                      vtx_coord,    &
                      tet_vtx,      &
                      vtx_ln_to_gn, &
                      tet_ln_to_gn)
  !----------------------------------------


  !----------------------------------------
  ! Create a PartMeshNodal instance to group all mesh data into a single structure.
  n_part = 1 ! Number of partitions (subdomains) on current MPI rank
  dim    = 3 ! Mesh dimension
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

  ! Add sections
  ! !   Triangles
  ! id_tri = PDM_part_mesh_nodal_section_add(mesh,                 & ! <- PartMeshNodal instance
  !                                          PDM_MESH_NODAL_TRIA3)   ! <- Triangles

  ! call PDM_part_mesh_nodal_section_std_set(mesh,               & ! <- PartMeshNodal instance
  !                                          id_tri,             & ! <- Triangle section ID
  !                                          0,                  & ! <- ID of current subdomain (i_part)
  !                                          n_tri,              & ! <- Number of triangles in current subdomain
  !                                          tri_vtx,            & ! <- Connectivity triangle->vtx in current subdomain
  !                                          tri_ln_to_gn,       & ! <- Global IDs of triangle in current subdomain
  !                                          parent_num,         & ! <- null()
  !                                          parent_gnum,        & ! <- null()
  !                                          PDM_OWNERSHIP_USER)   ! <- Ownership


  !   Tetrahedra
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
  call PDM_extract_part_create(extrp,                       & ! -> ExtractPart instance
                               dim,                         & ! <- Dimension
                               n_part,                      & ! <- Number of partitions per rank (input)
                               n_part,                      & ! <- Number of partitions per rank (output)
                               PDM_EXTRACT_PART_KIND_LOCAL, & ! <- Local extraction
                               PDM_SPLIT_DUAL_WITH_HILBERT, & ! <- Partitioning method (not used in LOCAL mode)
                               .true.,                      & ! <- Generate global IDs restricted to extraction
                               PDM_OWNERSHIP_KEEP,          & ! <- Ownership of extraction
                               comm)                          ! <- MPI communicator

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
                                          PDM_OWNERSHIP_USER)   ! <- Ownership (since paradigm-2.6.0)

  ! Realize extraction
  call PDM_extract_part_compute(extrp)


  ! Retrieve extracted mesh
  call PDM_extract_part_part_mesh_nodal_get(extrp,              & ! <- ExtractPart instance
                                            extracted_mesh,     & ! -> Extracted mesh (PartMeshNodal)
                                            PDM_OWNERSHIP_USER)   ! <- Ownership
  !----------------------------------------


  !----------------------------------------
  ! Visu
  call visu_pmn(comm, mesh,           "init")
  call visu_pmn(comm, extracted_mesh, "extract")
  !----------------------------------------


  !----------------------------------------
  ! Free memory
  deallocate(selected)
  call PDM_part_mesh_nodal_free(extracted_mesh)
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

    integer,                       intent(in)  :: comm
    integer,                       intent(in)  :: n_elt
    double precision,     pointer, intent(in)  :: elt_coord(:,:)
    integer(pdm_g_num_s), pointer, intent(out) :: elt_ln_to_gn(:)

    type(c_ptr)                                :: gen_gnum
    double precision,     pointer              :: char_length(:)

    char_length => null()

    call pdm_gnum_create(gen_gnum,           &
                         3,                  &
                         1,                  &
                         0,                  &
                         1.d-6,              &
                         comm,               &
                         PDM_OWNERSHIP_USER)

    call pdm_gnum_set_from_coords(gen_gnum,    &
                                  0,           &
                                  n_elt,       &
                                  elt_coord,   &
                                  char_length)

    call pdm_gnum_compute(gen_gnum)

    call pdm_gnum_get(gen_gnum,     &
                      0,            &
                      elt_ln_to_gn)

    call pdm_gnum_free(gen_gnum)

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

    integer,                       intent(in)  :: comm
    integer,                       intent(in)  :: n_vtx
    integer,                       intent(in)  :: n_tet
    double precision,     pointer, intent(in)  :: vtx_coord(:,:)
    integer(pdm_l_num_s), pointer, intent(in)  :: tet_vtx(:)
    integer(pdm_g_num_s), pointer, intent(out) :: vtx_ln_to_gn(:)
    integer(pdm_g_num_s), pointer, intent(out) :: tet_ln_to_gn(:)

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


  subroutine visu_pmn(comm, &
                      mesh, &
                      name)
    ! Export PartMeshNodal in Ensight format
    implicit none

    integer,            intent(in) :: comm
    type(c_ptr),        intent(in) :: mesh
    character(len = *), intent(in) :: name

    type(c_ptr)                    :: wrt
    integer                        :: id_geom
    integer                        :: id_var_part
    integer                        :: n_elt
    double precision, pointer      :: val_part(:) => null()

    integer                        :: i_rank, err

    call mpi_comm_rank(comm, i_rank, err)

    call PDM_writer_create(wrt,                      &
                           "Ensight",                &
                           PDM_WRITER_FMT_BIN,       &
                           PDM_WRITER_TOPO_VARIABLE, &
                           PDM_WRITER_OFF,           &
                           "extract_part_nodal_f",   &
                           name,                     &
                           comm,                     &
                           PDM_IO_KIND_MPI_SIMPLE,   &
                           1.d0,                     &
                           "")

    call PDM_writer_geom_create_from_mesh_nodal(wrt,     &
                                                id_geom, &
                                                name,    &
                                                mesh)

    call PDM_writer_var_create(wrt,                     &
                               id_var_part,             &
                               PDM_WRITER_OFF,          &
                               PDM_WRITER_VAR_SCALAIRE, &
                               PDM_WRITER_VAR_ELEMENTS, &
                               "i_part")

    call PDM_writer_step_beg(wrt, 0.d0)

    call PDM_writer_geom_write(wrt, id_geom)

    call PDM_part_mesh_nodal_section_n_elt_get(mesh, 0, 0, n_elt)
    allocate(val_part(n_elt))
    val_part(:) = i_rank

    call PDM_writer_var_set(wrt,         &
                            id_var_part, &
                            id_geom,     &
                            0,           &
                            val_part)

    call pdm_writer_var_write(wrt, id_var_part)
    deallocate(val_part)

    call PDM_writer_step_end(wrt)

    call PDM_writer_free(wrt)

    ! call PDM_part_mesh_nodal_dump_vtk(mesh, PDM_GEOMETRY_KIND_VOLUMIC, name)

  end subroutine visu_pmn
  !--------------------------------------------------------------------------------


end program extract_part_nodal_f