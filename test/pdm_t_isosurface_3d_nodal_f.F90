#include "pdm_configf.h"

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


program isosurface_3d_nodal

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use pdm_pointer_array
  use pdm_generate_mesh
  use pdm_dcube_gen
  use pdm_multipart
  use pdm_dcube_nodal_gen
  use pdm_part_mesh_nodal
  use pdm_isosurface
  use pdm_vtk
  use pdm_part_to_part
  use pdm_writer_wrapper
  use pdm_dmesh_nodal
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
  integer, parameter :: n_vtx_seg = 20
  integer(c_int)     :: n_part = 1

  ! part mesh definition
  type(c_ptr)               :: pmesh_nodal
  integer(pdm_l_num_s)      :: n_vtx
  double precision, pointer :: ipart_vtx_coord(:,:) => null()

  ! dist mesh definition
  type(c_ptr)               :: dmesh_nodal
  integer                   :: dn_vtx
  double precision, pointer :: dvtx_coords(:,:)

  ! Mesh field definition
  type(PDM_pointer_array_t), pointer :: array_field    => null() 
  double precision,          pointer :: ipart_field(:) => null()
  double precision,          pointer :: dfield(:)      => null()

  ! Isosurface structure
  type(c_ptr)               :: isos = C_NULL_PTR
  integer                   :: iso1, iso2, iso3, n_iso
  integer                   :: i_iso
  double precision, pointer :: plane_equation(:)
  double precision, pointer :: isovalues1(:), isovalues2(:), isovalues3(:)
  integer(c_int)            :: n_part_out = 1
  integer                   :: extract_kind = PDM_EXTRACT_PART_KIND_LOCAL


  ! Isosurface data
  !  part
  type(PDM_pointer_array_t), pointer :: array_isos_vtx_coord     => null()
  type(PDM_pointer_array_t), pointer :: array_isos_vtx_ln_to_gn  => null()
  type(PDM_pointer_array_t), pointer :: array_isos_face_vtx_idx  => null()
  type(PDM_pointer_array_t), pointer :: array_isos_face_vtx      => null()
  type(PDM_pointer_array_t), pointer :: array_isos_face_ln_to_gn => null()
  integer,                   pointer :: isos_face_vtx_idx(:)     => null()
  integer,                   pointer :: isos_face_vtx(:)         => null()
  integer(pdm_g_num_s),      pointer :: isos_face_ln_to_gn(:)    => null()  
  double precision,          pointer :: isos_vtx_coord(:,:)      => null()  
  integer(pdm_g_num_s),      pointer :: isos_vtx_ln_to_gn(:)     => null()  
  integer,                   pointer :: isos_n_vtx(:)            => null()
  integer,                   pointer :: isos_n_face(:)           => null()
  !  dist
  integer, pointer                   :: isos_dface_vtx_idx(:)    => null()
  integer(kind=pdm_g_num_s), pointer :: isos_dface_vtx(:)        => null()
  double precision,          pointer :: isos_dvtx_coord(:,:)     => null()
  integer                            :: isos_dn_vtx, isos_dn_face

  ! Isosurface field exchange
  type(c_ptr)                        :: ptp
  type(PDM_pointer_array_t), pointer :: array_iso_field           => null()  
  double precision,          pointer :: ipart_iso_field(:)        => null()
  double precision,          pointer :: interp_iso_field(:)       => null()
  type(my_field_t),          pointer :: interp_array_iso_field(:) => null()  
  integer :: request_vtx = -1
  !  part
  integer,          pointer :: pvtx_parent_idx(:)
  integer,          pointer :: pvtx_parent(:)
  double precision, pointer :: pvtx_parent_weight(:)
  integer                   :: i_parent
  !  dist
  integer,          pointer :: dparent_idx(:)
  double precision, pointer :: dparent_weight(:)

  ! Writer
  logical            :: visu = .false.
  character(len=256) :: filename
  character(len=8)   :: fmt = "(A11 I1)"

  integer :: i_rank, ierr
  integer :: i_part
  integer :: i_vtx, i_vtx_parent


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
      type(c_ptr),    value :: value

    end subroutine PDM_my_function
  end interface


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


  !---------------------------------------------------------------

  ! Initialize MPI
  call mpi_init (ierr)
  call mpi_comm_rank (comm, i_rank, ierr)


  ! Generate partitioned source mesh
  call mesh_gen(comm,                 &
                n_vtx_seg,            &
                n_vtx_seg,            &
                n_vtx_seg,            &
                n_part,               &
                PDM_MESH_NODAL_HEXA8, &
                1,                    &
                dmesh_nodal,          &
                pmesh_nodal)




  call pdm_isosurface_create (comm, &
                              3,    &
                              isos)



  if (n_part > 0) then
    call PDM_isosurface_part_mesh_nodal_set(isos, &
                                            pmesh_nodal)

    call PDM_pointer_array_create (array_field, &
                                   n_part,      &
                                   PDM_TYPE_DOUBLE)
    do i_part=1, n_part

      call PDM_part_mesh_nodal_n_vtx_get(pmesh_nodal, &
                                         i_part-1,    &
                                         n_vtx)

      call PDM_part_mesh_nodal_vtx_coord_get(pmesh_nodal,    &
                                            i_part-1,        &
                                            ipart_vtx_coord, &
                                            PDM_OWNERSHIP_KEEP)

      allocate(ipart_field(n_vtx))
      call PDM_field(n_vtx,           &
                     ipart_vtx_coord, &
                     ipart_field)

          
      call PDM_pointer_array_part_set(array_field, & ! <- Pointer array
                                      i_part-1,    & ! <- ID of current part
                                      ipart_field )  ! <- Field

    end do
  else 

    call PDM_pointer_array_create (array_field, &
                                   1,           &
                                   PDM_TYPE_DOUBLE)

    call PDM_isosurface_dmesh_nodal_set(isos, &
                                        dmesh_nodal)

    dn_vtx = pdm_dmesh_nodal_n_vtx_get (dmesh_nodal)

    allocate(dfield(dn_vtx))

    call PDM_DMesh_nodal_vtx_get (dmesh_nodal, &
                                  dvtx_coords, &
                                  PDM_OWNERSHIP_KEEP)

    call PDM_field(dn_vtx, &
                   dvtx_coords, &
                   dfield)

    call PDM_pointer_array_part_set(array_field, & ! <- Pointer array
                                    0,           & ! <- ID of current part
                                    dfield )       ! <- Field    

  end if



  allocate(plane_equation(3))
  plane_equation = [1.0d0, 1.0d0, 1.0d0]
  allocate(isovalues1(1))
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

  if (n_part > 0) then
    do i_part=1, n_part

      call PDM_pointer_array_part_get(array_field, &
                                      i_part-1,    &
                                      ipart_field)

      call PDM_isosurface_pfield_set(isos,     &
                                     iso3,     &
                                     i_part-1, &
                                     ipart_field);

    end do
  else 
    call PDM_isosurface_dfield_set(isos, &
                                   iso3, &
                                   dfield);
  end if


  if (n_part > 0) then
    call PDM_isosurface_redistribution_set (isos,         &
                                            extract_kind, &
                                            PDM_SPLIT_DUAL_WITH_PARMETIS)
  
    if (extract_kind == PDM_EXTRACT_PART_KIND_REEQUILIBRATE)then
      call PDM_isosurface_n_part_out_set (isos,       &
                                          n_part_out)
    else
      n_part_out = n_part
    end if
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

    if (n_part_out > 0) then
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
          call PDM_part_to_part_reverse_iexch (ptp,                                    &
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
                                                  i_iso-1,              &
                                                  i_part-1,             &
                                                  PDM_MESH_ENTITY_VTX,  &
                                                  isos_n_vtx(i_part),   &
                                                  pvtx_parent_idx,      &
                                                  pvtx_parent_weight,   &
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

      write(filename, fmt) "isosurface_", i_iso-1
      if (visu) then
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
      endif
      call PDM_pointer_array_free (array_isos_face_vtx_idx)
      call PDM_pointer_array_free (array_isos_face_vtx)
      call PDM_pointer_array_free (array_isos_face_ln_to_gn)
      call PDM_pointer_array_free (array_isos_vtx_coord)
      call PDM_pointer_array_free (array_isos_vtx_ln_to_gn)
      call PDM_pointer_array_free (interp_array_iso_field(1)%pa)

    else 
      call PDM_isosurface_part_to_part_get (isos,                &
                                            i_iso-1,             &
                                            PDM_MESH_ENTITY_VTX, &
                                            ptp,                 &
                                            PDM_OWNERSHIP_USER)

       call PDM_isosurface_dconnectivity_get (isos,                           &
                                              i_iso-1,                        &
                                              PDM_CONNECTIVITY_TYPE_FACE_VTX, &
                                              isos_dn_face,                   &
                                              isos_dface_vtx_idx,             &
                                              isos_dface_vtx,                 &
                                              PDM_OWNERSHIP_KEEP)

       call PDM_isosurface_dvtx_coord_get (isos,            &
                                           i_iso-1,         &
                                           isos_dn_vtx,     &
                                           isos_dvtx_coord, &
                                           PDM_OWNERSHIP_KEEP)

       call  PDM_isosurface_dparent_weight_get(isos,                &
                                               i_iso-1,             &
                                               PDM_MESH_ENTITY_VTX, &
                                               isos_dn_vtx,         &
                                               dparent_idx,         &
                                               dparent_weight,      &
                                               PDM_OWNERSHIP_KEEP)

      call PDM_part_to_part_reverse_iexch (ptp,                                   &
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
                                       0,               &
                                       ipart_iso_field)

      allocate(interp_iso_field(isos_dn_vtx))
      do i_vtx=1, isos_dn_vtx
        interp_iso_field(i_vtx) = 0.0d0
        do i_vtx_parent=dparent_idx(i_vtx), dparent_idx(i_vtx+1)-1
          interp_iso_field(i_vtx) = interp_iso_field(i_vtx) + dparent_weight(i_vtx_parent+1) * ipart_iso_field(i_vtx_parent+1)
        end do 
      enddo


      call PDM_pointer_array_free (array_iso_field)

    end if

  end do

  call PDM_pointer_array_free (array_field)





  contains


  subroutine mesh_gen(comm,        &
                      n_x,         &
                      n_y,         &
                      n_z,         &
                      n_part,      &
                      elt_type,    &
                      order,       &
                      dmesh_nodal, &
                      pmesh_nodal)
    implicit none

    integer :: comm
    integer :: n_x
    integer :: n_y
    integer :: n_z
    integer :: n_part
    integer :: elt_type
    integer :: order
    type(c_ptr) :: dmesh_nodal
    type(c_ptr) :: pmesh_nodal

    type(c_ptr) :: dcube
    integer :: n_domain = 1
    integer(kind = PDM_l_num_s), pointer :: n_part_domain(:)
    double precision, pointer :: part_fraction(:) => null()
    type(c_ptr) :: multipart

    allocate(n_part_domain(n_domain))
    n_part_domain(1) = n_part

    call PDM_dcube_nodal_gen_create(dcube,    &
                                    comm,     &
                                    n_x,      &
                                    n_y,      &
                                    n_z,      &
                                    10.0d0,   &
                                    -5.0d0,   &
                                    -5.0d0,   &
                                    -5.0d0,   &
                                    elt_type, &
                                    order,    &
                                    PDM_OWNERSHIP_KEEP)

    call PDM_dcube_nodal_gen_build(dcube, &
                                   dmesh_nodal)


    ! call PDM_dmesh_nodal_generate_distribution(dmesh_nodal)

    if (n_part > 0) then
      call PDM_multipart_create(multipart,                   &
                                n_domain,                    &
                                n_part_domain,               &
                                PDM_FALSE,                   &
                                PDM_SPLIT_DUAL_WITH_PARMETIS, &
                                PDM_PART_SIZE_HOMOGENEOUS,   &
                                part_fraction,               &
                                comm,                        &
                                PDM_OWNERSHIP_KEEP)

      call PDM_multipart_dmesh_nodal_set(multipart, &
                                          0,        &
                                          dmesh_nodal)

      call PDM_multipart_compute(multipart)

      call PDM_multipart_get_part_mesh_nodal(multipart,   &
                                             0,           &
                                             pmesh_nodal, &
                                             PDM_OWNERSHIP_KEEP);
    end if

  end subroutine mesh_gen


  subroutine PDM_field(n_vtx,      &
                       vtx_coords, &
                       field)
    implicit none

    integer, intent(in) :: n_vtx
    double precision, intent(in), pointer :: vtx_coords(:,:)
    double precision, intent(out), pointer :: field(:)

    double precision, dimension(n_vtx) :: x, y, z, center

    x = vtx_coords(1,:)
    y = vtx_coords(2,:)
    z = vtx_coords(3,:)

    center = 3.0


    ! do i_vtx=1, n_vtx
      ! field(i_vtx) = vtx_coords(1,i_vtx)*vtx_coords(1,i_vtx) + vtx_coords(2,i_vtx)*vtx_coords(2,i_vtx) + vtx_coords(3,i_vtx)*vtx_coords(3,i_vtx)
    ! end do
    field(:) = (x(:)-center)*(x(:)-center) + (y(:)-center)*(y(:)-center) + (z(:)-center)*(z(:)-center) - 2.0d0

  end subroutine PDM_field


end program isosurface_3d_nodal
