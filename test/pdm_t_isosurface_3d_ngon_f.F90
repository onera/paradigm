#include "pdm_configf.h"

subroutine PDM_my_function(x, &
                           y, &
                           z, &
                           value) &
bind(c)
  use iso_c_binding
  implicit none

  real(c_double), value :: x
  real(c_double), value :: y
  real(c_double), value :: z
  type(c_ptr), value :: value

  double precision, pointer :: f_value(:)

  call c_f_pointer(value, f_value, [1])

  f_value(1) = 0.5d0*(x*x*x*x + y*y*y*y + z*z*z*z) - 8.0d0*(x*x + y*y + z*z) + 60.0d0

end subroutine



program tp_mesh_location

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use pdm_pointer_array
  use pdm_generate_mesh
  use pdm_dcube_gen
  use pdm_isosurface
  use pdm_vtk
  use pdm_part_to_part
  use pdm_writer_wrapper
  use iso_c_binding

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif


  !---------------------------------------------------------------
  integer,              parameter    :: comm = MPI_COMM_WORLD

  ! mesh generation
  integer(pdm_g_num_s), parameter    :: n_vtx_seg = 50
  integer(c_int),       parameter    :: n_part = 3
  integer(c_int),       parameter    :: n_part_out = 1


  integer                            :: order                          = 1
  type(c_ptr)                        :: ho_ordering                    = C_NULL_PTR  


  integer(pdm_l_num_s),      pointer :: n_vtx(:)      => null()
  integer(pdm_l_num_s),      pointer :: n_edge(:)     => null()
  integer(pdm_l_num_s),      pointer :: n_face(:)     => null()
  integer(pdm_l_num_s),      pointer :: n_cell(:)     => null()
  integer(pdm_l_num_s),      pointer :: n_surface(:)  => null()
  integer(pdm_l_num_s),      pointer :: n_ridge(:)    => null()
  type(PDM_pointer_array_t), pointer :: vtx_coord     => null()
  type(PDM_pointer_array_t), pointer :: edge_vtx      => null()
  type(PDM_pointer_array_t), pointer :: face_edge_idx => null()
  type(PDM_pointer_array_t), pointer :: face_edge     => null()
  type(PDM_pointer_array_t), pointer :: face_vtx      => null()
  type(PDM_pointer_array_t), pointer :: cell_face_idx => null()
  type(PDM_pointer_array_t), pointer :: cell_face     => null()
  type(PDM_pointer_array_t), pointer :: surface_face_idx => null()
  type(PDM_pointer_array_t), pointer :: surface_face     => null()
  type(PDM_pointer_array_t), pointer :: ridge_edge_idx   => null()
  type(PDM_pointer_array_t), pointer :: ridge_edge       => null()  
  type(PDM_pointer_array_t), pointer :: vtx_ln_to_gn     => null()
  type(PDM_pointer_array_t), pointer :: edge_ln_to_gn    => null()
  type(PDM_pointer_array_t), pointer :: face_ln_to_gn    => null()
  type(PDM_pointer_array_t), pointer :: cell_ln_to_gn    => null()
  type(PDM_pointer_array_t), pointer :: surface_ln_to_gn => null()
  type(PDM_pointer_array_t), pointer :: ridge_ln_to_gn   => null()


  ! mesh location structure
  type(c_ptr)                        :: isos = C_NULL_PTR
  integer                            :: id_isosurface
  double precision, pointer          :: plane_equation(:)
  double precision, pointer          :: plane_isovalues(:)

  ! Writer
  integer :: visu = 0
  character(len=256) :: filename
  character(len=17)  :: fmt = "(A18 I1 A8 I1 A4)"

  ! src mesh definition
  logical,              parameter    :: nodal = .false.
  integer(pdm_g_num_s),      pointer :: ipart_cell_ln_to_gn(:) => null()
  integer(pdm_l_num_s),      pointer :: ipart_cell_face_idx(:) => null()
  integer(pdm_l_num_s),      pointer :: ipart_cell_face(:)     => null()  
  integer(pdm_g_num_s),      pointer :: ipart_face_ln_to_gn(:) => null()
  integer(pdm_l_num_s),      pointer :: ipart_face_edge_idx(:) => null()
  integer(pdm_l_num_s),      pointer :: ipart_face_edge(:)     => null()
  ! integer(pdm_l_num_s),      pointer :: ipart_face_vtx(:)      => null()
  integer(pdm_g_num_s),      pointer :: ipart_edge_ln_to_gn(:) => null()
  integer(pdm_l_num_s),      pointer :: ipart_edge_vtx(:)      => null()
  integer(pdm_g_num_s),      pointer :: ipart_vtx_ln_to_gn(:)  => null()
  double precision,          pointer :: ipart_vtx_coord(:,:)   => null()

  integer(pdm_g_num_s),      pointer :: ipart_surface_ln_to_gn(:) => null()
  integer(pdm_l_num_s),      pointer :: ipart_surface_face_idx(:) => null()
  integer(pdm_l_num_s),      pointer :: ipart_surface_face(:)     => null()


  integer :: isos_n_face
  integer :: isos_n_vtx

  double precision,          pointer :: isos_vtx_coord(:,:)   => null()  
  integer(pdm_g_num_s),      pointer :: isos_vtx_ln_to_gn(:)  => null()  
  integer,                   pointer :: isos_face_vtx_idx(:)  => null()
  integer,                   pointer :: isos_face_vtx(:)      => null()
  integer(pdm_g_num_s),      pointer :: isos_face_ln_to_gn(:) => null()  
  integer,                   pointer :: isos_face_color(:) => null()  



  type(c_ptr)           :: ptp
  integer :: n_part1, n_part2


  integer :: i_rank, ierr
  integer :: i_part


  interface

    subroutine PDM_my_function(x, &
                               y, &
                               z, &
                               value) &
    bind(c)
      use iso_c_binding
      implicit none

      real(c_double), value :: x
      real(c_double), value :: y
      real(c_double), value :: z
      type(c_ptr), value :: value

    end subroutine

  end interface


  !---------------------------------------------------------------

  ! Initialize MPI
  call mpi_init (ierr)
  call mpi_comm_rank (comm, i_rank, ierr)


  ! Generate partitioned source mesh
  call PDM_generate_mesh_parallelepiped_ngon (comm,                         &
                                               PDM_MESH_NODAL_HEXA8,       &
                                               order,                       &
                                               ho_ordering,                 &
                                               -5.0d0,                      &
                                               -5.0d0,                      &
                                               -5.0d0,                      &
                                               10.d0,                        &
                                               10.d0,                        &
                                               10.d0,                        &
                                               n_vtx_seg,                   &
                                               n_vtx_seg,                   &
                                               n_vtx_seg,                   &
                                               n_part,                      &
                                               PDM_SPLIT_DUAL_WITH_HILBERT, &
                                               n_vtx,                       &
                                               n_edge,                      &
                                               n_face,                      &
                                               n_cell,                      &
                                               vtx_coord,                   &
                                               edge_vtx,                    &
                                               face_edge_idx,               &
                                               face_edge,                   &
                                               face_vtx,                    &
                                               cell_face_idx,               &
                                               cell_face,                   &
                                               vtx_ln_to_gn,                &
                                               edge_ln_to_gn,               &
                                               face_ln_to_gn,               &
                                               cell_ln_to_gn,               &
                                               n_surface,                   &
                                               surface_face_idx,            &
                                               surface_face,                &
                                               surface_ln_to_gn,            &
                                               n_ridge,                     &
                                               ridge_edge_idx,              &
                                               ridge_edge,                  &
                                               ridge_ln_to_gn)


  call pdm_isosurface_create (comm, &
                              3,    &
                              isos)

  call PDM_isosurface_n_part_set (isos, &
                                  n_part)

  do i_part=1, n_part

    ! CELL FACE CONNECTIVITY
    call PDM_pointer_array_part_get (cell_face_idx, &
                                     i_part-1,      &
                                     ipart_cell_face_idx)

    call PDM_pointer_array_part_get (cell_face, &
                                     i_part-1,  &
                                     ipart_cell_face)

    call PDM_isosurface_pconnectivity_set (isos,                            &
                                           i_part-1,                        &
                                           PDM_CONNECTIVITY_TYPE_CELL_FACE, &
                                           n_cell(i_part),                  &
                                           ipart_cell_face_idx,             &
                                           ipart_cell_face)

    ! FACE EDGE CONNECTIVITY
    call PDM_pointer_array_part_get (face_edge_idx, &
                                     i_part-1,      &
                                     ipart_face_edge_idx)

    call PDM_pointer_array_part_get (face_edge, &
                                     i_part-1,  &
                                     ipart_face_edge)

    call PDM_isosurface_pconnectivity_set (isos,                            &
                                           i_part-1,                        &
                                           PDM_CONNECTIVITY_TYPE_FACE_EDGE, &
                                           n_face(i_part),                  &
                                           ipart_face_edge_idx,             &
                                           ipart_face_edge)

    ! EDGE VTX CONNECTIVITY
    call PDM_pointer_array_part_get (edge_vtx, &
                                     i_part-1, &
                                     ipart_edge_vtx)

    call PDM_isosurface_pconnectivity_set (isos,                           &
                                           i_part-1,                       &
                                           PDM_CONNECTIVITY_TYPE_EDGE_VTX, &
                                           n_edge(i_part),                 &
                                           null(),                         &
                                           ipart_edge_vtx)

    ! VTX COORDS
    call PDM_pointer_array_part_get (vtx_coord,                 &
                                     i_part-1,                  &
                                     PDM_STRIDE_CST_INTERLACED, &
                                     3,                         &
                                     ipart_vtx_coord)

    call PDM_isosurface_pvtx_coord_set (isos,          &
                                        i_part-1,      &
                                        n_vtx(i_part), &
                                        ipart_vtx_coord)


    ! CELL LN_TO_GN
    call PDM_pointer_array_part_get (cell_ln_to_gn, &
                                     i_part-1,      &
                                     ipart_cell_ln_to_gn)

    call PDM_isosurface_ln_to_gn_set (isos,                 &
                                      i_part-1,             &
                                      PDM_MESH_ENTITY_CELL, &
                                      ipart_cell_ln_to_gn)


    ! FACE LN_TO_GN
    call PDM_pointer_array_part_get (face_ln_to_gn, &
                                     i_part-1,      &
                                     ipart_face_ln_to_gn)

    call PDM_isosurface_ln_to_gn_set (isos,                 &
                                      i_part-1,             &
                                      PDM_MESH_ENTITY_FACE, &
                                      ipart_face_ln_to_gn)


    ! EDGE LN_TO_GN
    call PDM_pointer_array_part_get(edge_ln_to_gn, &
                                    i_part-1,      &
                                    ipart_edge_ln_to_gn)

    call PDM_isosurface_ln_to_gn_set (isos,                 &
                                      i_part-1,             &
                                      PDM_MESH_ENTITY_EDGE, &
                                      ipart_edge_ln_to_gn)


    ! VTX LN_TO_GN
    call PDM_pointer_array_part_get(vtx_ln_to_gn, &
                                    i_part-1,     &
                                    ipart_vtx_ln_to_gn)


    call PDM_isosurface_ln_to_gn_set (isos,                &
                                      i_part-1,            &
                                      PDM_MESH_ENTITY_VTX, &
                                      ipart_vtx_ln_to_gn)

    ! GROUPS
    call PDM_pointer_array_part_get(surface_face_idx, &
                                    i_part-1,         &
                                    ipart_surface_face_idx)

    call PDM_pointer_array_part_get(surface_face, &
                                    i_part-1,     &
                                    ipart_surface_face)

    call PDM_pointer_array_part_get(surface_ln_to_gn, &
                                    i_part-1,         &
                                    ipart_surface_ln_to_gn)


    call PDM_isosurface_n_group_set (isos,                  &
                                      PDM_MESH_ENTITY_FACE, &
                                      n_surface(i_part))

    call PDM_isosurface_pgroup_set (isos,                   &
                                    i_part-1,               &
                                    PDM_MESH_ENTITY_FACE,   &
                                    ipart_surface_face_idx, &
                                    ipart_surface_face,     &
                                    ipart_surface_ln_to_gn)
  end do




  allocate(plane_equation(3))
  plane_equation = [1.0d0, -1.0d0, 0.0d0]
  allocate(plane_isovalues(1))
  plane_isovalues = [0.0d0]


  ! call PDM_isosurface_add (isos,                        &
  !                          PDM_ISO_SURFACE_KIND_PLANE, &
  !                          1,                           &
  !                          plane_isovalues,             &
  !                          id_isosurface)

  ! call PDM_isosurface_equation_set(isos,          &
  !                                  id_isosurface, &
  !                                  plane_equation);


  call PDM_isosurface_add(isos, &
                          PDM_ISO_SURFACE_KIND_FUNCTION, &
                          1, &
                          plane_isovalues, &
                          id_isosurface)

  call PDM_isosurface_field_function_set(isos, &
                                         id_isosurface, &
                                         PDM_my_function);




  call PDM_isosurface_redistribution_set (isos, &
                                          PDM_EXTRACT_PART_KIND_REEQUILIBRATE, &
                                          1)


  call PDM_isosurface_n_part_out_set (isos,       &
                                      n_part_out)



  call PDM_isosurface_part_to_part_enable(isos,                &
                                          id_isosurface,       &
                                          PDM_MESH_ENTITY_VTX, &
                                          0);

  call PDM_isosurface_part_to_part_enable(isos,                 &
                                          id_isosurface,        &
                                          PDM_MESH_ENTITY_EDGE, &
                                          0);

  call PDM_isosurface_part_to_part_enable(isos,                 &
                                          id_isosurface,        &
                                          PDM_MESH_ENTITY_FACE, &
                                          0);





  call PDM_isosurface_compute (isos, &
                               id_isosurface)


  call PDM_isosurface_part_to_part_get (isos,          &
                                        id_isosurface, &
                                        PDM_MESH_ENTITY_VTX,   &
                                        ptp,           &
                                        PDM_OWNERSHIP_USER)

  call PDM_part_to_part_n_part_get (ptp,     &
                                    n_part1, &
                                    n_part2)


  !  Write geometry
  if (visu == 1) then
    do i_part = 1, n_part_out

    
      isos_n_face = PDM_isosurface_pconnectivity_get (isos,                           &
                                                      id_isosurface,                  &
                                                      i_part-1,                       &
                                                      PDM_CONNECTIVITY_TYPE_FACE_VTX, &
                                                      isos_face_vtx_idx,              &
                                                      isos_face_vtx,                  &
                                                      PDM_OWNERSHIP_KEEP)

      isos_n_face = PDM_isosurface_ln_to_gn_get (isos,                 &
                                                 id_isosurface,        &
                                                 i_part-1,             &
                                                 PDM_MESH_ENTITY_FACE, &
                                                 isos_face_ln_to_gn,   &
                                                 PDM_OWNERSHIP_KEEP)


      isos_n_vtx = PDM_isosurface_pvtx_coord_get (isos,           &
                                                  id_isosurface,  &
                                                  i_part-1,       &
                                                  isos_vtx_coord, &
                                                  PDM_OWNERSHIP_KEEP)

      isos_n_vtx = PDM_isosurface_ln_to_gn_get (isos,                &
                                                id_isosurface,       &
                                                i_part-1,            &
                                                PDM_MESH_ENTITY_VTX, &
                                                isos_vtx_ln_to_gn,   &
                                                PDM_OWNERSHIP_KEEP)

      ! filename = "isosurface_function.vtk" 
      write(filename, fmt) "isosurface_i_part_", i_part, "_i_rank_", i_rank, ".vtk"
      ! write(*,*) ">", trim(filename), "<"

      call PDM_vtk_write_polydata (trim(filename), &
                                   isos_n_vtx, &
                                   isos_vtx_coord, &
                                   isos_vtx_ln_to_gn, &
                                   isos_n_face, &
                                   isos_face_vtx_idx, &
                                   isos_face_vtx, &
                                   isos_face_ln_to_gn, &
                                   isos_face_color)

    end do
  end if




end program tp_mesh_location
