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
  double precision :: a, b, c, d, e

  call c_f_pointer(value, f_value, [1])

  a = 0.0d0
  b = -3.5d0
  c = 0.0d0
  d = 11.0d0
  e = 43.0d0

  ! a = 0.5d0
  ! b = 0.0d0
  ! c = -8.0d0
  ! d = 0.0d0
  ! e = 60.0d0

  f_value(1) = a*(x*x*x*x + y*y*y*y + z*z*z*z) &
             + b*(x*x*x + y*y*y + z*z*z)       &
             + c*(x*x + y*y + z*z)             &
             + d*(x + y + z)                   &
             + e

end subroutine PDM_my_function


program isosurface_3d_ngon

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

  ! Parsing option
  integer :: i_arg
  character(len=99) :: arg

  ! mesh generation
  integer(pdm_g_num_s), parameter :: n_vtx_seg   = 49
  integer                         :: order       = 1
  type(c_ptr)                     :: ho_ordering = C_NULL_PTR  
  integer(c_int)                  :: n_part      = 4
  ! part mesh definition
  integer(pdm_l_num_s), pointer :: n_vtx(:)                  => null()
  integer(pdm_l_num_s), pointer :: n_edge(:)                 => null()
  integer(pdm_l_num_s), pointer :: n_face(:)                 => null()
  integer(pdm_l_num_s), pointer :: n_cell(:)                 => null()
  integer(pdm_l_num_s), pointer :: n_surface(:)              => null()
  integer(pdm_l_num_s), pointer :: n_ridge(:)                => null()
  integer(pdm_l_num_s), pointer :: ipart_cell_face_idx(:)    => null()
  integer(pdm_l_num_s), pointer :: ipart_cell_face(:)        => null()
  integer(pdm_g_num_s), pointer :: ipart_cell_ln_to_gn(:)    => null()
  integer(pdm_l_num_s), pointer :: ipart_face_edge_idx(:)    => null()
  integer(pdm_l_num_s), pointer :: ipart_face_edge(:)        => null()
  integer(pdm_g_num_s), pointer :: ipart_face_ln_to_gn(:)    => null()
  integer(pdm_l_num_s), pointer :: ipart_edge_vtx(:)         => null()
  integer(pdm_g_num_s), pointer :: ipart_edge_ln_to_gn(:)    => null()
  double precision,     pointer :: ipart_vtx_coord(:,:)      => null()
  integer(pdm_g_num_s), pointer :: ipart_vtx_ln_to_gn(:)     => null()
  integer(pdm_l_num_s), pointer :: ipart_surface_face_idx(:) => null()
  integer(pdm_l_num_s), pointer :: ipart_surface_face(:)     => null()
  integer(pdm_g_num_s), pointer :: ipart_surface_ln_to_gn(:) => null()
  ! pointer array
  type(PDM_pointer_array_t), pointer :: cell_face_idx    => null()
  type(PDM_pointer_array_t), pointer :: cell_face        => null()
  type(PDM_pointer_array_t), pointer :: cell_ln_to_gn    => null()
  type(PDM_pointer_array_t), pointer :: face_edge_idx    => null()
  type(PDM_pointer_array_t), pointer :: face_edge        => null()
  type(PDM_pointer_array_t), pointer :: face_ln_to_gn    => null()
  type(PDM_pointer_array_t), pointer :: edge_vtx         => null()
  type(PDM_pointer_array_t), pointer :: edge_ln_to_gn    => null()
  type(PDM_pointer_array_t), pointer :: face_vtx         => null()
  type(PDM_pointer_array_t), pointer :: vtx_coord        => null()
  type(PDM_pointer_array_t), pointer :: vtx_ln_to_gn     => null()
  type(PDM_pointer_array_t), pointer :: surface_face_idx => null()
  type(PDM_pointer_array_t), pointer :: surface_face     => null()
  type(PDM_pointer_array_t), pointer :: surface_ln_to_gn => null()
  type(PDM_pointer_array_t), pointer :: ridge_edge_idx   => null()
  type(PDM_pointer_array_t), pointer :: ridge_edge       => null()
  type(PDM_pointer_array_t), pointer :: ridge_ln_to_gn   => null()

  ! Mesh field definition
  type(PDM_pointer_array_t), pointer :: array_field    => null()
  double precision,          pointer :: ipart_field(:) => null()

  ! Isosurface structure
  type(c_ptr)               :: isos = C_NULL_PTR
  integer                   :: iso1, iso2, iso3, n_iso
  integer                   :: i_iso
  double precision, pointer :: plane_equation(:)
  double precision, pointer :: isovalues1(:), isovalues2(:), isovalues3(:)
  integer(c_int)            :: n_part_out = 1
  integer                   :: extract_kind = PDM_EXTRACT_PART_KIND_LOCAL

  ! Isosurface data
  type(PDM_pointer_array_t), pointer :: array_isos_vtx_coord     => null()
  type(PDM_pointer_array_t), pointer :: array_isos_vtx_ln_to_gn  => null()
  type(PDM_pointer_array_t), pointer :: array_isos_face_vtx_idx  => null()
  type(PDM_pointer_array_t), pointer :: array_isos_face_vtx      => null()
  type(PDM_pointer_array_t), pointer :: array_isos_face_ln_to_gn => null()
  double precision,          pointer :: isos_vtx_coord(:,:)      => null()
  integer(pdm_g_num_s),      pointer :: isos_vtx_ln_to_gn(:)     => null()
  integer,                   pointer :: isos_face_vtx_idx(:)     => null()
  integer,                   pointer :: isos_face_vtx(:)         => null()
  integer(pdm_g_num_s),      pointer :: isos_face_ln_to_gn(:)    => null()
  integer,                   pointer :: isos_n_vtx(:)            => null()
  integer,                   pointer :: isos_n_face(:)           => null()

  ! Isosurface field exchange
  type(c_ptr)                        :: ptp
  type(PDM_pointer_array_t), pointer :: array_iso_field           => null()
  double precision,          pointer :: ipart_iso_field(:)        => null()
  double precision,          pointer :: interp_iso_field(:)       => null()
  type(my_field_t),          pointer :: interp_array_iso_field(:) => null()
  integer                            :: request_vtx = -1
  integer,                   pointer :: pvtx_parent_idx(:)
  integer,                   pointer :: pvtx_parent(:)
  double precision,          pointer :: pvtx_parent_weight(:)
  integer                            :: i_parent

  ! Writer
  logical            :: visu = .true.
  character(len=256) :: filename
  character(len=8)   :: fmt = "(A11 I1)"

  integer :: i_rank, ierr
  integer :: i_part
  integer :: i_vtx, i_vtx_parent

  interface
    subroutine PDM_my_function(x,     &
                               y,     &
                               z,     &
                               value) &
    bind(c)
      use iso_c_binding
      implicit none

      real(c_double), value :: x
      real(c_double), value :: y
      real(c_double), value :: z
      type(c_ptr), value :: value

    end subroutine PDM_my_function
  end interface


  !---------------------------------------------------------------

  ! Initialize MPI
  call mpi_init (ierr)
  call mpi_comm_rank (comm, i_rank, ierr)


  !----------------------------------------
  ! Parse command line arguments
  i_arg = 1
  do while (i_arg <= command_argument_count())
    call get_command_argument(i_arg, arg)

    select case(arg)

      case ("-reequilibrate")
        extract_kind = PDM_EXTRACT_PART_KIND_REEQUILIBRATE

      case ("-n_part")
        i_arg = i_arg + 1
        call get_command_argument(i_arg, arg)
        read(arg, *) n_part

      case ("-n_part_out")
        i_arg = i_arg + 1
        call get_command_argument(i_arg, arg)
        read(arg, *) n_part_out

      case ("-visu")
        visu = .true.

    endselect

    i_arg = i_arg + 1
  enddo





  ! Generate partitioned source mesh
  call PDM_generate_mesh_parallelepiped_ngon (comm,                         &
                                              PDM_MESH_NODAL_TETRA4,        &
                                              order,                        &
                                              ho_ordering,                  &
                                              -5.0d0,                       &
                                              -5.0d0,                       &
                                              -5.0d0,                       &
                                              10.d0,                        &
                                              10.d0,                        &
                                              10.d0,                        &
                                              n_vtx_seg,                    &
                                              n_vtx_seg,                    &
                                              n_vtx_seg,                    &
                                              n_part,                       &
                                              PDM_SPLIT_DUAL_WITH_PARMETIS, &
                                              n_vtx,                        &
                                              n_edge,                       &
                                              n_face,                       &
                                              n_cell,                       &
                                              vtx_coord,                    &
                                              edge_vtx,                     &
                                              face_edge_idx,                &
                                              face_edge,                    &
                                              face_vtx,                     &
                                              cell_face_idx,                &
                                              cell_face,                    &
                                              vtx_ln_to_gn,                 &
                                              edge_ln_to_gn,                &
                                              face_ln_to_gn,                &
                                              cell_ln_to_gn,                &
                                              n_surface,                    &
                                              surface_face_idx,             &
                                              surface_face,                 &
                                              surface_ln_to_gn,             &
                                              n_ridge,                      &
                                              ridge_edge_idx,               &
                                              ridge_edge,                   &
                                              ridge_ln_to_gn)




  call pdm_isosurface_create (comm, &
                              3,    &
                              isos)

  call PDM_isosurface_n_part_set (isos, &
                                  n_part)

  call PDM_pointer_array_create (array_field, &
                                 n_part,      &
                                 PDM_TYPE_DOUBLE)

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
    ! FIELD
    allocate(ipart_field(n_vtx(i_part)))
    call PDM_field(n_vtx(i_part), &
                   ipart_vtx_coord, &
                   ipart_field)

        
    call PDM_pointer_array_part_set(array_field, & ! <- Pointer array
                                    i_part-1,    & ! <- ID of current part
                                    ipart_field )  ! <- Field

  end do

  allocate(plane_equation(3))
  plane_equation = [1.0d0, 1.0d0, 1.0d0]
  allocate(isovalues1(2))
  allocate(isovalues2(1))
  allocate(isovalues3(1))
  isovalues1 = [-4.5d0]
  isovalues2 = [0.0d0]
  isovalues3 = [0.0d0]

  n_iso = 3

  call PDM_isosurface_add (isos,                       &
                           PDM_ISO_SURFACE_KIND_PLANE, &
                           1,                          &
                           isovalues1,                 &
                           iso1)

  call PDM_isosurface_equation_set(isos, &
                                   iso1, &
                                   plane_equation);


  call PDM_isosurface_add(isos,                          &
                          PDM_ISO_SURFACE_KIND_FUNCTION, &
                          1,                             &
                          isovalues2,                    &
                          iso2)

  call PDM_isosurface_field_function_set(isos, &
                                         iso2, &
                                         PDM_my_function);



  call PDM_isosurface_add(isos,                       &
                          PDM_ISO_SURFACE_KIND_FIELD, &
                          1,                          &
                          isovalues3,                 &
                          iso3)

  do i_part=1, n_part

    call PDM_pointer_array_part_get(array_field, &
                                    i_part-1,    &
                                    ipart_field)

    call PDM_isosurface_pfield_set(isos,     &
                                   iso3,     &
                                   i_part-1, &
                                   ipart_field);

  end do


  call PDM_isosurface_redistribution_set (isos,         &
                                          extract_kind, &
                                          PDM_SPLIT_DUAL_WITH_PARMETIS)

  if (extract_kind == PDM_EXTRACT_PART_KIND_REEQUILIBRATE) then
    call PDM_isosurface_n_part_out_set (isos, &
                                        n_part_out)
  else
    n_part_out = n_part
  end if



  do i_iso=1, n_iso


    call PDM_isosurface_part_to_part_enable(isos,                &
                                            i_iso-1,             &
                                            PDM_MESH_ENTITY_VTX, &
                                            0);


    call PDM_isosurface_compute (isos, &
                                 i_iso-1)
  end do

  allocate(isos_n_face(n_part_out))
  allocate(isos_n_vtx(n_part_out))
  allocate(interp_array_iso_field(1))
  interp_array_iso_field(1)%name = "field"


  !  Write geometry
  do i_iso = 1, n_iso

    call PDM_isosurface_part_to_part_get (isos,                &
                                          i_iso-1,             &
                                          PDM_MESH_ENTITY_VTX, &
                                          ptp,                 &
                                          PDM_OWNERSHIP_USER)

    call PDM_pointer_array_create (array_isos_face_vtx_idx, &
                                    n_part_out,             &
                                    PDM_TYPE_INT)

    call PDM_pointer_array_create (array_isos_face_vtx, &
                                    n_part_out,         &
                                    PDM_TYPE_INT)

    call PDM_pointer_array_create (array_isos_face_ln_to_gn, &
                                    n_part_out,              &
                                    PDM_TYPE_G_NUM)

    call PDM_pointer_array_create (array_isos_vtx_coord, &
                                    n_part_out,          &
                                    PDM_TYPE_DOUBLE)

    call PDM_pointer_array_create (array_isos_vtx_ln_to_gn, &
                                    n_part_out,             &
                                    PDM_TYPE_G_NUM)

    call PDM_pointer_array_create (interp_array_iso_field(1)%pa, &
                                    n_part_out,                  &
                                    PDM_TYPE_DOUBLE)

    do i_part = 1, n_part_out

       call PDM_isosurface_pconnectivity_get (isos,                           &
                                              i_iso-1,                        &
                                              i_part-1,                       &
                                              PDM_CONNECTIVITY_TYPE_FACE_VTX, &
                                              isos_n_face(i_part),            &
                                              isos_face_vtx_idx,              &
                                              isos_face_vtx,                  &
                                              PDM_OWNERSHIP_KEEP)

       call PDM_isosurface_ln_to_gn_get (isos,                 &
                                         i_iso-1,              &
                                         i_part-1,             &
                                         PDM_MESH_ENTITY_FACE, &
                                         isos_n_face(i_part),  &
                                         isos_face_ln_to_gn,   &
                                         PDM_OWNERSHIP_KEEP)


       call PDM_isosurface_pvtx_coord_get (isos,               &
                                           i_iso-1,            &
                                           i_part-1,           &
                                           isos_n_vtx(i_part), &
                                           isos_vtx_coord,     &
                                           PDM_OWNERSHIP_KEEP)

       call PDM_isosurface_ln_to_gn_get (isos,                &
                                         i_iso-1,             &
                                         i_part-1,            &
                                         PDM_MESH_ENTITY_VTX, &
                                         isos_n_vtx(i_part),  &
                                         isos_vtx_ln_to_gn,   &
                                         PDM_OWNERSHIP_KEEP)

      call PDM_pointer_array_part_set(array_isos_face_vtx_idx, & ! <- Pointer array
                                      i_part-1,                & ! <- ID of current part
                                      isos_face_vtx_idx )        ! <- Field

      call PDM_pointer_array_part_set(array_isos_face_vtx, & ! <- Pointer array
                                      i_part-1,            & ! <- ID of current part
                                      isos_face_vtx)         ! <- Field
      
      call PDM_pointer_array_part_set(array_isos_face_ln_to_gn, & ! <- Pointer array
                                      i_part-1,                 & ! <- ID of current part
                                      isos_face_ln_to_gn )        ! <- Field

      call PDM_pointer_array_part_set(array_isos_vtx_coord, & ! <- Pointer array
                                      i_part-1,             & ! <- ID of current part
                                      isos_vtx_coord)         ! <- Field

      call PDM_pointer_array_part_set(array_isos_vtx_ln_to_gn, & ! <- Pointer array
                                      i_part-1,                & ! <- ID of current part
                                      isos_vtx_ln_to_gn )        ! <- Field

      if (extract_kind == PDM_EXTRACT_PART_KIND_REEQUILIBRATE) then

        call PDM_part_to_part_reverse_iexch (ptp,                           &
                                     PDM_MPI_COMM_KIND_COLLECTIVE,          &
                                     PDM_STRIDE_CST_INTERLACED,             &
                                     PDM_PART_TO_PART_DATA_DEF_ORDER_PART2, &
                                     1,                                     &
                                     null(),                                &
                                     array_field,                           &
                                     null(),                                &
                                     array_iso_field,                       &
                                     request_vtx)

        call PDM_part_to_part_reverse_iexch_wait(ptp, request_vtx)

        call PDM_pointer_array_part_get (array_iso_field, &
                                         i_part-1,        &
                                         ipart_iso_field)



         call PDM_isosurface_pparent_weight_get(isos, &
                                                i_iso-1,               &
                                                i_part-1,              &
                                                PDM_MESH_ENTITY_VTX,   &
                                                isos_n_vtx(i_part),    &
                                                pvtx_parent_idx,            &
                                                pvtx_parent_weight,         &
                                                PDM_OWNERSHIP_KEEP)

        allocate(interp_iso_field(isos_n_vtx(i_part)))
        do i_vtx=1, isos_n_vtx(i_part)
          interp_iso_field(i_vtx) = 0.0d0
          do i_vtx_parent=pvtx_parent_idx(i_vtx), pvtx_parent_idx(i_vtx+1)-1
            interp_iso_field(i_vtx) = interp_iso_field(i_vtx) + pvtx_parent_weight(i_vtx_parent+1) * ipart_iso_field(i_vtx_parent+1)
          end do 
        enddo

        call PDM_pointer_array_free (array_iso_field)

      else
        call PDM_pointer_array_part_get (array_field, &
                                         i_part-1,    &
                                         ipart_field)

         call  PDM_isosurface_pparent_weight_get(isos,                &
                                                 i_iso-1,             &
                                                 i_part-1,            &
                                                 PDM_MESH_ENTITY_VTX, &
                                                 isos_n_vtx(i_part),  &
                                                 pvtx_parent_idx,     &
                                                 pvtx_parent_weight,  &
                                                 PDM_OWNERSHIP_KEEP)

         call PDM_isosurface_plocal_parent_get(isos,                &
                                               i_iso-1,             &
                                               i_part-1,            &
                                               PDM_MESH_ENTITY_VTX, &
                                               isos_n_vtx(i_part),  &
                                               pvtx_parent_idx,     &
                                               pvtx_parent,         &
                                               PDM_OWNERSHIP_KEEP)


        allocate(interp_iso_field(isos_n_vtx(i_part)))
        do i_vtx=1, isos_n_vtx(i_part)
          interp_iso_field(i_vtx) = 0.0d0
          do i_vtx_parent=pvtx_parent_idx(i_vtx), pvtx_parent_idx(i_vtx+1)-1
            i_parent = pvtx_parent(i_vtx_parent+1)
            interp_iso_field(i_vtx) = interp_iso_field(i_vtx) + pvtx_parent_weight(i_vtx_parent+1) * ipart_field(i_parent)
          end do 
        enddo

      end if


      call PDM_pointer_array_part_set(interp_array_iso_field(1)%pa, & ! <- Pointer array
                                      i_part-1,                     & ! <- ID of current part
                                      interp_iso_field )              ! <- Field


    end do
    if (visu) then
      write(filename, fmt) "isosurface_", i_iso-1
      call writer_wrapper(comm,                     &
                          i_rank,                   &
                          ".",                      &
                          filename,                 &
                          n_part_out,               &
                          isos_n_vtx,               &
                          array_isos_vtx_coord,     &
                          array_isos_vtx_ln_to_gn,  &
                          isos_n_face,              &
                          array_isos_face_vtx_idx,  &
                          array_isos_face_vtx,      &
                          array_isos_face_ln_to_gn, &
                          vtx_field=interp_array_iso_field)
    end if

    call PDM_pointer_array_free (array_isos_face_vtx_idx)
    call PDM_pointer_array_free (array_isos_face_vtx)
    call PDM_pointer_array_free (array_isos_face_ln_to_gn)
    call PDM_pointer_array_free (array_isos_vtx_coord)
    call PDM_pointer_array_free (array_isos_vtx_ln_to_gn)
    call PDM_pointer_array_free (interp_array_iso_field(1)%pa)
  end do

  call PDM_pointer_array_free (array_field)


  contains

  subroutine PDM_field(n_vtx,      &
                       vtx_coords, &
                       field)
    use iso_c_binding
    implicit none

    integer,          intent(in)           :: n_vtx
    double precision, intent(in),  pointer :: vtx_coords(:,:)
    double precision, intent(out), pointer :: field(:)

    double precision, dimension(n_vtx) :: x, y, z, center

    x = vtx_coords(1,:)
    y = vtx_coords(2,:)
    z = vtx_coords(3,:)

    center = 3.0

    field(:) = (x(:)-center)*(x(:)-center) + (y(:)-center)*(y(:)-center) + (z(:)-center)*(z(:)-center) - 2.0d0

  end subroutine PDM_field


end program isosurface_3d_ngon
