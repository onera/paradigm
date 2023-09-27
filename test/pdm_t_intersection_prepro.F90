#include "pdm_configf.h"

program testf

  use iso_c_binding
  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use pdm_io
  use pdm_writer
  use pdm_sphere_surf_gen
  use pdm_fortran
  use pdm_generate_mesh
  use pdm_mesh_intersection

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  integer, parameter                 :: comm                           = MPI_COMM_WORLD
  character(len=99)                  :: arg
  integer                            :: i

  integer                            :: order                          = 1
  type(c_ptr)                        :: ho_ordering                    = C_NULL_PTR
  integer                            :: part_method                    = PDM_SPLIT_DUAL_WITH_HILBERT

  integer                            :: cube_elt_type                  = PDM_MESH_NODAL_HEXA8
  double precision                   :: cube_xmin                      = -0.5 
  double precision                   :: cube_ymin                      = -0.5
  double precision                   :: cube_zmin                      = -0.5
  double precision                   :: cube_lengthx                   = 1.
  double precision                   :: cube_lengthy                   = 1.
  double precision                   :: cube_lengthz                   = 1.
  integer(kind=pdm_g_num_s)          :: cube_n                         = 10
  integer(kind=pdm_g_num_s)          :: cube_n_x                       = -1
  integer(kind=pdm_g_num_s)          :: cube_n_y                       = -1
  integer(kind=pdm_g_num_s)          :: cube_n_z                       = -1
  integer                            :: cube_n_part                    = 1
  integer                            :: cube_dim                       = 3

  integer                            :: sphere_elt_type                = PDM_MESH_NODAL_TRIA3
  double precision                   :: sphere_radius                  = 0.35
  double precision                   :: sphere_center_x                = 0.
  double precision                   :: sphere_center_y                = 0.
  double precision                   :: sphere_center_z                = 0.
  integer(kind=pdm_g_num_s)          :: sphere_n                       = 10
  integer(kind=pdm_g_num_s)          :: sphere_n_u                     = -1
  integer(kind=pdm_g_num_s)          :: sphere_n_v                     = -1
  integer                            :: sphere_n_part                  = 1
  integer                            :: sphere_dim                     = 2

  integer(kind=pdm_l_num_s), pointer :: cube_pn_vtx(:)                 => null()
  integer(kind=pdm_l_num_s), pointer :: cube_pn_edge(:)                => null()
  integer(kind=pdm_l_num_s), pointer :: cube_pn_face(:)                => null()
  integer(kind=pdm_l_num_s), pointer :: cube_pn_cell(:)                => null()
  integer(kind=pdm_l_num_s), pointer :: cube_pn_surface(:)             => null()
  integer(kind=pdm_l_num_s), pointer :: cube_pn_ridge(:)               => null()
  type(PDM_pointer_array_t), pointer :: cube_pvtx_coord                => null()
  type(PDM_pointer_array_t), pointer :: cube_pedge_vtx                 => null()
  type(PDM_pointer_array_t), pointer :: cube_pface_edge_idx            => null()
  type(PDM_pointer_array_t), pointer :: cube_pface_edge                => null()
  type(PDM_pointer_array_t), pointer :: cube_pface_vtx                 => null()
  type(PDM_pointer_array_t), pointer :: cube_pcell_face_idx            => null()
  type(PDM_pointer_array_t), pointer :: cube_pcell_face                => null()
  type(PDM_pointer_array_t), pointer :: cube_pvtx_ln_to_gn             => null()
  type(PDM_pointer_array_t), pointer :: cube_pedge_ln_to_gn            => null()
  type(PDM_pointer_array_t), pointer :: cube_pface_ln_to_gn            => null()
  type(PDM_pointer_array_t), pointer :: cube_pcell_ln_to_gn            => null()
  type(PDM_pointer_array_t), pointer :: cube_psurface_face_idx         => null()
  type(PDM_pointer_array_t), pointer :: cube_psurface_face             => null()
  type(PDM_pointer_array_t), pointer :: cube_psurface_face_ln_to_gn    => null()
  type(PDM_pointer_array_t), pointer :: cube_pridge_edge_idx           => null()
  type(PDM_pointer_array_t), pointer :: cube_pridge_edge               => null()
  type(PDM_pointer_array_t), pointer :: cube_pridge_edge_ln_to_gn      => null()


  integer(kind=pdm_l_num_s), pointer :: sphere_pn_vtx(:)               => null()
  integer(kind=pdm_l_num_s), pointer :: sphere_pn_edge(:)              => null()
  integer(kind=pdm_l_num_s), pointer :: sphere_pn_face(:)              => null()
  type(PDM_pointer_array_t), pointer :: sphere_pvtx_coord              => null()
  type(PDM_pointer_array_t), pointer :: sphere_pedge_vtx               => null()
  type(PDM_pointer_array_t), pointer :: sphere_pface_edge_idx          => null()
  type(PDM_pointer_array_t), pointer :: sphere_pface_edge              => null()
  type(PDM_pointer_array_t), pointer :: sphere_pface_vtx               => null()
  type(PDM_pointer_array_t), pointer :: sphere_pvtx_ln_to_gn           => null()
  type(PDM_pointer_array_t), pointer :: sphere_pedge_ln_to_gn          => null()
  type(PDM_pointer_array_t), pointer :: sphere_pface_ln_to_gn          => null()
 
  integer(kind=pdm_l_num_s)          :: cube_extr_mesh_n_cell            = -1
  integer(kind=pdm_l_num_s), pointer :: cube_extr_mesh_cell_face(:)      => null()
  integer(kind=pdm_l_num_s), pointer :: cube_extr_mesh_cell_face_idx(:)  => null()
  integer(kind=pdm_l_num_s)          :: cube_extr_mesh_n_face            = -1
  integer(kind=pdm_l_num_s), pointer :: cube_extr_mesh_face_vtx(:)       => null()
  integer(kind=pdm_l_num_s), pointer :: cube_extr_mesh_face_vtx_idx(:)   => null()
  integer(kind=pdm_l_num_s)          :: sphere_extr_mesh_n_face          = -1
  integer(kind=pdm_l_num_s), pointer :: sphere_extr_mesh_face_vtx(:)     => null()
  integer(kind=pdm_l_num_s), pointer :: sphere_extr_mesh_face_vtx_idx(:) => null()

  ! Writer
  type(c_ptr)                       :: wrt
  integer                           :: id_geom_sphere
  integer                           :: id_geom_cube
  integer                           :: id_block

  ! Mesh intersecion
  type(c_ptr)                       :: mi
  double precision, parameter        :: project_coeff = 0.0 ! Unused

  ! MPI
  integer                           :: code
  integer                           :: i_rank
  integer                           :: n_rank

  integer                           :: ipart

  integer(pdm_l_num_s), pointer :: unassociated_ptr(:)            => null()

  double precision,     pointer :: ipart_cube_vtx_coord(:,:)      => null()
  integer(pdm_g_num_s), pointer :: ipart_cube_vtx_ln_to_gn(:)     => null()
  double precision,     pointer :: ipart_sphere_vtx_coord(:,:)    => null()
  integer(pdm_g_num_s), pointer :: ipart_sphere_vtx_ln_to_gn(:)   => null()

  integer(pdm_l_num_s), pointer :: ipart_cube_face_edge_idx(:)    => null()
  integer(pdm_l_num_s), pointer :: ipart_cube_face_vtx(:)         => null()
  integer(pdm_l_num_s), pointer :: ipart_sphere_face_edge_idx(:)  => null()
  integer(pdm_l_num_s), pointer :: ipart_sphere_face_vtx(:)       => null()

  integer(pdm_l_num_s), pointer :: ipart_cube_cell_face_idx(:)    => null()
  integer(pdm_l_num_s), pointer :: ipart_cube_cell_face(:)        => null()
  integer(pdm_g_num_s), pointer :: ipart_cube_cell_ln_to_gn(:)    => null()
  integer(pdm_g_num_s), pointer :: ipart_cube_face_ln_to_gn(:)    => null()
  integer(pdm_g_num_s), pointer :: ipart_sphere_face_ln_to_gn(:)  => null()

  integer(pdm_l_num_s), pointer :: box_cube_box_sphere_idx(:)     => null()
  integer(pdm_l_num_s), pointer :: box_cube_box_sphere(:)         => null()

  integer(pdm_l_num_s), pointer :: cptr_int_null(:)               => null()
  integer(pdm_g_num_s), pointer :: cptr_gnum_null(:)              => null()

  type(c_ptr)                   :: cube_extr_mesh                 = C_NULL_PTR
  type(c_ptr)                   :: sphere_extr_mesh               = C_NULL_PTR

  !-----------------------------------------------------------
  !                Read command line arguments
  !-----------------------------------------------------------

  i = 1
  do while (i <= command_argument_count())
    call get_command_argument(i, arg)
    select case(arg)
      case ("-n_cube")
        i = i+1
        call get_command_argument(i, arg)
        read(arg, *) cube_n
      case ("-n_sphere")
        i = i+1
        call get_command_argument(i, arg)
        read(arg, *) sphere_n
      case ("-n_part_cube")
        i = i+1
        call get_command_argument(i, arg)
        read(arg, *) cube_n_part
      case ("-n_part_sphere")
        i = i+1
        call get_command_argument(i, arg)
        read(arg, *) sphere_n_part
      case ("-parmetis")
        part_method = PDM_SPLIT_DUAL_WITH_PARMETIS
      case ("-pt-scotch")
        part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH
      case ("-hilbert")
        part_method = PDM_SPLIT_DUAL_WITH_HILBERT
      case default
        print *, "Invalid command argument ", arg
        stop
    end select

    i = i + 1
  enddo

  cube_n_x = cube_n
  cube_n_y = cube_n
  cube_n_z = cube_n

  sphere_n_u = 2 * sphere_n
  sphere_n_v = sphere_n

  !-----------------------------------------------------------
  !                  MPI init
  !-----------------------------------------------------------

  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)

  !-----------------------------------------------------------
  !                Generate partitionned cube
  !-----------------------------------------------------------

  print *, "*** 1 ***"

  call PDM_generate_mesh_parallelepiped_ngon(comm,                        &
                                             cube_elt_type,               & 
                                             order,                       &
                                             ho_ordering,                 &
                                             cube_xmin,                   &
                                             cube_ymin,                   &
                                             cube_zmin,                   &
                                             cube_lengthx,                &
                                             cube_lengthy,                &
                                             cube_lengthz,                &
                                             cube_n_x,                    &
                                             cube_n_y,                    &
                                             cube_n_z,                    &
                                             cube_n_part,                 &
                                             part_method,                 &
                                             cube_pn_vtx,                 &
                                             cube_pn_edge,                &
                                             cube_pn_face,                &
                                             cube_pn_cell,                &
                                             cube_pvtx_coord,             &
                                             cube_pedge_vtx,              &
                                             cube_pface_edge_idx,         &
                                             cube_pface_edge,             &
                                             cube_pface_vtx,              &
                                             cube_pcell_face_idx,         &
                                             cube_pcell_face,             &
                                             cube_pvtx_ln_to_gn,          &
                                             cube_pedge_ln_to_gn,         &
                                             cube_pface_ln_to_gn,         &
                                             cube_pcell_ln_to_gn,         &
                                             cube_pn_surface,             &
                                             cube_psurface_face_idx,      &
                                             cube_psurface_face,          &
                                             cube_psurface_face_ln_to_gn, &
                                             cube_pn_ridge,               &
                                             cube_pridge_edge_idx,        &
                                             cube_pridge_edge,            &
                                             cube_pridge_edge_ln_to_gn)

  !-----------------------------------------------------------
  !                Generate partitionned sphere
  !-----------------------------------------------------------
  print *, "*** 2 ***"

  call PDM_generate_mesh_sphere_ngon(comm,                  &
                                     sphere_elt_type,       &
                                     order,                 &
                                     ho_ordering,           &
                                     sphere_radius,         &
                                     sphere_center_x,       &
                                     sphere_center_y,       &
                                     sphere_center_z,       &
                                     sphere_n_u,            &
                                     sphere_n_v,            &
                                     sphere_n_part,         &
                                     part_method,           &
                                     sphere_pn_vtx,         &
                                     sphere_pn_edge,        &
                                     sphere_pn_face,        &
                                     sphere_pvtx_coord,     & 
                                     sphere_pedge_vtx,      &
                                     sphere_pface_edge_idx, &
                                     sphere_pface_edge,     &
                                     sphere_pface_vtx,      &
                                     sphere_pvtx_ln_to_gn,  &
                                     sphere_pedge_ln_to_gn, &
                                     sphere_pface_ln_to_gn)

  !-----------------------------------------------------------
  !                   Writer ensight
  !-----------------------------------------------------------

  print *, "*** 3 ***"

  call PDM_writer_create(wrt,                    &
                         "Ensight",              &
                         PDM_WRITER_FMT_ASCII,   &
                         PDM_WRITER_TOPO_CST,    &
                         PDM_WRITER_OFF,         &
                         "intersec_cube_sphere", &
                         "intersec_cube_sphere", &
                         comm,                   &
                         PDM_IO_KIND_MPI_SIMPLE, &
                         1.d0,                   &
                         "")

  print *, "*** 4 ***"

  call PDM_writer_geom_create(wrt,          &
                              id_geom_cube, &
                              "cube",       &
                              cube_n_part)

  call PDM_writer_geom_create(wrt,            &
                              id_geom_sphere, &
                              "sphere",       &
                              sphere_n_part)

  call PDM_writer_step_beg(wrt, 0.d0)


  print *, "*** 5 ***"

  !  Write geometry
  do ipart = 1, cube_n_part
    call PDM_pointer_array_part_get(cube_pvtx_ln_to_gn, &
                                    ipart-1,       &
                                    ipart_cube_vtx_ln_to_gn)

    call PDM_pointer_array_part_get(cube_pvtx_coord,           &
                                    ipart-1,                   &
                                    PDM_STRIDE_CST_INTERLACED, &
                                    3,                         &
                                    ipart_cube_vtx_coord)
  
    call PDM_writer_geom_coord_set(wrt,                        &
                                   id_geom_cube,               &
                                   ipart-1,                    &
                                   cube_pn_vtx(ipart),         &
                                   ipart_cube_vtx_coord,       &
                                   ipart_cube_vtx_ln_to_gn,    &
                                   PDM_OWNERSHIP_USER)

    call PDM_pointer_array_part_get(cube_pcell_ln_to_gn,       &
                                    ipart-1,                   &
                                    ipart_cube_cell_ln_to_gn)

    call PDM_pointer_array_part_get(cube_pface_edge_idx,       &
                                    ipart-1,                   &
                                    ipart_cube_face_edge_idx)

    call PDM_pointer_array_part_get(cube_pface_vtx,            &
                                    ipart-1,                   &
                                    ipart_cube_face_vtx)

    call PDM_pointer_array_part_get(cube_pcell_face_idx,       &
                                    ipart-1,                   &
                                    ipart_cube_cell_face_idx)

    call PDM_pointer_array_part_get(cube_pcell_face,           &
                                    ipart-1,                   &
                                    ipart_cube_cell_face)

    call PDM_writer_geom_cell3d_cellface_add (wrt,                      &
                                              id_geom_cube,             &
                                              ipart-1,                  &
                                              cube_pn_cell(ipart),      &
                                              cube_pn_face(ipart),      &
                                              ipart_cube_face_edge_idx, &
                                              unassociated_ptr,         &
                                              ipart_cube_face_vtx,      &
                                              ipart_cube_cell_face_idx, &
                                              unassociated_ptr,         &
                                              ipart_cube_cell_face,     &
                                              ipart_cube_cell_ln_to_gn)

  end do

  print *, "*** 6 ***"

  call PDM_writer_geom_bloc_add(wrt,                  &
                                id_geom_sphere,       &
                                PDM_MESH_NODAL_TRIA3, & ! PDM_MESH_NODAL_TRIA3
                                PDM_OWNERSHIP_USER,   &
                                id_block)

  do ipart = 1, sphere_n_part
    call PDM_pointer_array_part_get(sphere_pvtx_ln_to_gn,      &
                                    ipart-1,                   &
                                    ipart_sphere_vtx_ln_to_gn)

    call PDM_pointer_array_part_get(sphere_pvtx_coord,         &
                                    ipart-1,                   &
                                    PDM_STRIDE_CST_INTERLACED, &
                                    3,                         &
                                    ipart_sphere_vtx_coord)
  
    call PDM_writer_geom_coord_set(wrt,                        &
                                   id_geom_sphere,             &
                                   ipart-1,                    &
                                   sphere_pn_vtx(ipart),       &
                                   ipart_sphere_vtx_coord,     &
                                   ipart_sphere_vtx_ln_to_gn,  &
                                   PDM_OWNERSHIP_USER)

    call PDM_pointer_array_part_get(sphere_pface_ln_to_gn,  &
                                    ipart-1,                &
                                    ipart_sphere_face_ln_to_gn)

    call PDM_pointer_array_part_get(sphere_pface_edge_idx, &
                                    ipart-1,               &
                                    ipart_sphere_face_edge_idx)

    call PDM_pointer_array_part_get(sphere_pface_vtx,      &
                                    ipart-1,               &    
                                    ipart_sphere_face_vtx)

    call PDM_writer_geom_bloc_std_set(wrt,            &
                                      id_geom_sphere, &
                                      id_block,       &
                                      ipart-1,        &
                                      sphere_pn_face(ipart), &
                                      ipart_sphere_face_vtx,       &
                                      ipart_sphere_face_ln_to_gn)


  end do

  print *, "*** 7 ***"

  call PDM_writer_geom_write(wrt,             &
                             id_geom_cube)

  call PDM_writer_geom_write(wrt,             &
                             id_geom_sphere)

  call PDM_writer_step_end(wrt)

  print *, "*** 8 ***"

  !-----------------------------------------------------------
  !        Search for intersection candidates
  !           result : List of face for each cell 
  !-----------------------------------------------------------

  call PDM_mesh_intersection_create (mi,                                    &
                                     PDM_MESH_INTERSECTION_KIND_PREPROCESS, &
                                     cube_dim,                              &
                                     sphere_dim,                            &
                                     project_coeff,                         &
                                     comm,                                  &
                                     PDM_OWNERSHIP_USER)

  call PDM_mesh_intersection_n_part_set (mi,          &
                                         0,           &  ! Attention commence a 0
                                         cube_n_part)

  call PDM_mesh_intersection_n_part_set (mi,            &
                                         1,             &  ! Attention commence a 0
                                         sphere_n_part)


  do ipart = 1, cube_n_part
    call PDM_pointer_array_part_get(cube_pvtx_ln_to_gn, &
                                    ipart-1,       &
                                    ipart_cube_vtx_ln_to_gn)

    call PDM_pointer_array_part_get(cube_pvtx_coord,           &
                                    ipart-1,                   &
                                    PDM_STRIDE_CST_INTERLACED, &
                                    3,                         &
                                    ipart_cube_vtx_coord)
  
    call PDM_pointer_array_part_get(cube_pcell_ln_to_gn,       &
                                    ipart-1,                   &
                                    ipart_cube_cell_ln_to_gn)

    call PDM_pointer_array_part_get(cube_pface_ln_to_gn,       &
                                    ipart-1,                   &
                                    ipart_cube_face_ln_to_gn)

    call PDM_pointer_array_part_get(cube_pface_edge_idx,       &
                                    ipart-1,                   &
                                    ipart_cube_face_edge_idx)

    call PDM_pointer_array_part_get(cube_pface_vtx,            &
                                    ipart-1,                   &
                                    ipart_cube_face_vtx)

    call PDM_pointer_array_part_get(cube_pcell_face_idx,       &
                                    ipart-1,                   &
                                    ipart_cube_cell_face_idx)

    call PDM_pointer_array_part_get(cube_pcell_face,           &
                                    ipart-1,                   &
                                    ipart_cube_cell_face)

    call PDM_mesh_intersection_part_set (mi,                            &
                                         0,                             &
                                         ipart-1,                       &
                                         cube_pn_cell(ipart),           &
                                         cube_pn_face(ipart),           &
                                         0,                             &
                                         cube_pn_vtx(ipart),            &
                                         ipart_cube_cell_face_idx,      &
                                         ipart_cube_cell_face,          &
                                         cptr_int_null,                 &
                                         cptr_int_null,                 &
                                         cptr_int_null,                 &
                                         ipart_cube_face_edge_idx,      & ! Same as ipart_cube_face_vtx_idx
                                         ipart_cube_face_vtx,           &
                                         ipart_cube_cell_ln_to_gn,      &
                                         ipart_cube_face_ln_to_gn,      &
                                         cptr_gnum_null,                &
                                         ipart_cube_vtx_ln_to_gn,       &
                                         ipart_cube_vtx_coord)      
  enddo 

  do ipart = 1, sphere_n_part
    call PDM_pointer_array_part_get(sphere_pvtx_ln_to_gn,      &
                                    ipart-1,                   &
                                    ipart_sphere_vtx_ln_to_gn)

    call PDM_pointer_array_part_get(sphere_pvtx_coord,         &
                                    ipart-1,                   &
                                    PDM_STRIDE_CST_INTERLACED, &
                                    3,                         &
                                    ipart_sphere_vtx_coord)
  
    call PDM_pointer_array_part_get(sphere_pface_ln_to_gn,  &
                                    ipart-1,                &
                                    ipart_sphere_face_ln_to_gn)

    call PDM_pointer_array_part_get(sphere_pface_edge_idx, &
                                    ipart-1,               &
                                    ipart_sphere_face_edge_idx)

    call PDM_pointer_array_part_get(sphere_pface_vtx,      &
                                    ipart-1,               &
                                    ipart_sphere_face_vtx)

    call PDM_mesh_intersection_part_set (mi,                            &
                                         0,                             &
                                         ipart-1,                       &
                                         0,                             &
                                         sphere_pn_face(ipart),         &
                                         0,                             &
                                         sphere_pn_vtx(ipart),          &
                                         cptr_int_null,                 &
                                         cptr_int_null,                 &
                                         cptr_int_null,                 &
                                         cptr_int_null,                 &
                                         cptr_int_null,                 &
                                         ipart_sphere_face_edge_idx,    & ! Same as ipart_cube_face_vtx_idx
                                         ipart_sphere_face_vtx,         &
                                         cptr_gnum_null,                &
                                         ipart_sphere_face_ln_to_gn,    &
                                         cptr_gnum_null,                &
                                         ipart_sphere_vtx_ln_to_gn,     &
                                         ipart_sphere_vtx_coord)      

  enddo 

  call  PDM_mesh_intersection_compute (mi)

  !-----------------------------------------------------------
  !  Get the intersections preprocessing results
  !-----------------------------------------------------------

  call PDM_mesh_intersection_preprocessing_get (mi,                      &
                                                box_cube_box_sphere_idx, & 
                                                box_cube_box_sphere,     &
                                                cube_extr_mesh,          &
                                                sphere_extr_mesh)

  !-----------------------------------------------------------
  !  Free mesh intersection structure
  !-----------------------------------------------------------

  call PDM_mesh_intersection_free(mi)

  !-----------------------------------------------------------
  !  Get connectivities of extract meshes             
  !-----------------------------------------------------------

  call PDM_extract_part_connectivity_get (cube_extr_mesh,                  &
                                          0,                               & ! A unique partition
                                          PDM_CONNECTIVITY_TYPE_CELL_FACE, &
                                          cube_extr_mesh_n_cell,           &
                                          cube_extr_mesh_cell_face,        &
                                          cube_extr_mesh_cell_face_idx,    &
                                          PDM_OWNERSHIP_KEEP)

  call PDM_extract_part_connectivity_get (cube_extr_mesh,                  &
                                          0,                               & ! A unique partition
                                          PDM_CONNECTIVITY_TYPE_FACE_VTX,  &
                                          cube_extr_mesh_n_face,           &
                                          cube_extr_mesh_face_vtx,         &
                                          cube_extr_mesh_face_vtx_idx,     &
                                          PDM_OWNERSHIP_KEEP)

  call PDM_extract_part_connectivity_get (sphere_extr_mesh,                &
                                          0,                               & ! A unique partition
                                          PDM_CONNECTIVITY_TYPE_FACE_VTX,  &
                                          sphere_extr_mesh_n_face,         &
                                          sphere_extr_mesh_face_vtx,       &
                                          sphere_extr_mesh_face_vtx_idx,   &
                                          PDM_OWNERSHIP_KEEP)

  !-----------------------------------------------------------
  !    Push cube results to the initial partitionning 
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  !    Push sphere results to the initial partitionning 
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  !    Free extract meshes 
  !-----------------------------------------------------------

  call PDM_fortran_free_c(c_loc(box_cube_box_sphere_idx))
  call PDM_fortran_free_c(c_loc(box_cube_box_sphere))

  call PDM_extract_part_free(cube_extr_mesh)
  call PDM_extract_part_free(sphere_extr_mesh)

  !-----------------------------------------------------------
  !    Finalize writer 
  !-----------------------------------------------------------

  call PDM_writer_free(wrt)

  print *, "*** 9 ***"

  !-----------------------------------------------------------
  !    Free structures 
  !-----------------------------------------------------------

  call PDM_fortran_free_c(c_loc(cube_pn_vtx))
  call PDM_fortran_free_c(c_loc(cube_pn_edge))
  call PDM_fortran_free_c(c_loc(cube_pn_face))
  call PDM_fortran_free_c(c_loc(cube_pn_cell))
  call PDM_fortran_free_c(c_loc(cube_pn_surface))
  call PDM_fortran_free_c(c_loc(cube_pn_ridge))

  call PDM_pointer_array_free(cube_pvtx_coord)
  call PDM_pointer_array_free(cube_pedge_vtx)
  call PDM_pointer_array_free(cube_pface_edge_idx)
  call PDM_pointer_array_free(cube_pface_edge)
  call PDM_pointer_array_free(cube_pface_vtx)
  call PDM_pointer_array_free(cube_pcell_face_idx)
  call PDM_pointer_array_free(cube_pcell_face)
  call PDM_pointer_array_free(cube_pvtx_ln_to_gn)
  call PDM_pointer_array_free(cube_pedge_ln_to_gn)
  call PDM_pointer_array_free(cube_pface_ln_to_gn)
  call PDM_pointer_array_free(cube_pcell_ln_to_gn)
  call PDM_pointer_array_free(cube_psurface_face_idx)
  call PDM_pointer_array_free(cube_psurface_face)
  call PDM_pointer_array_free(cube_psurface_face_ln_to_gn)
  call PDM_pointer_array_free(cube_pridge_edge_idx)
  call PDM_pointer_array_free(cube_pridge_edge)
  call PDM_pointer_array_free(cube_pridge_edge_ln_to_gn)

  call PDM_fortran_free_c(c_loc(sphere_pn_vtx))
  call PDM_fortran_free_c(c_loc(sphere_pn_edge))
  call PDM_fortran_free_c(c_loc(sphere_pn_face))

  call PDM_pointer_array_free(sphere_pvtx_coord)
  call PDM_pointer_array_free(sphere_pedge_vtx)
  call PDM_pointer_array_free(sphere_pface_edge_idx)
  call PDM_pointer_array_free(sphere_pface_edge)
  call PDM_pointer_array_free(sphere_pface_vtx)
  call PDM_pointer_array_free(sphere_pvtx_ln_to_gn)
  call PDM_pointer_array_free(sphere_pedge_ln_to_gn)
  call PDM_pointer_array_free(sphere_pface_ln_to_gn)

  print *, "*** 10 ***"

 
  if (i_rank .eq. 0) then
    write(*, *) "-- End"
  end if

  call mpi_finalize(code)

end program testf
