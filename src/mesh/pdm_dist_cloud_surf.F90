#include "pdm_configf.h"

module pdm_dist_cloud_surf

  use pdm
  use pdm_pointer_array

  implicit none

  interface

    function PDM_dist_cloud_surf_cloud_n_part_get_cf(dcs, &
                                                     i_point_cloud) &
      result (n_part) &
      bind(c, name='PDM_dist_cloud_surf_cloud_n_part_get')

      use iso_c_binding
      use pdm

      implicit none

      type(c_ptr),    value :: dcs
      integer(c_int), value :: i_point_cloud
      integer(c_int)        :: n_part

    end function PDM_dist_cloud_surf_cloud_n_part_get_cf


    subroutine PDM_dist_cloud_surf_cloud_dim_get_cf(dcs,           &
                                                    i_point_cloud, &
                                                    i_part,        &
                                                    n_points)      &

      bind (c, name = 'PDM_dist_cloud_surf_cloud_dim_get')

      use iso_c_binding
      implicit none

      type(c_ptr),    value :: dcs
      integer(c_int), value :: i_point_cloud
      integer(c_int), value :: i_part
      integer(c_int)        :: n_points

    end subroutine PDM_dist_cloud_surf_cloud_dim_get_cf


  end interface


  contains


  subroutine pdm_dist_cloud_surf_create(dcs,           &
                                        mesh_nature,   &
                                        n_point_cloud, &
                                        f_comm,        &
                                        owner)

    ! Create a structure to compute distance between point clouds and a surface mesh
    implicit none

    integer, intent(in) :: mesh_nature    ! Nature of the mesh
    integer, intent(in) :: n_point_cloud  ! Number of point clouds
    integer, intent(in) :: f_comm         ! MPI communicator
    integer, intent(in) :: owner          ! Ownership of results

    type (c_ptr)        :: dcs            ! PDM_dist_cloud_surf_t instance

    integer(c_int) :: c_mesh_nature
    integer(c_int) :: c_n_point_cloud
    type (c_ptr)   :: c_comm
    integer(c_int) :: c_owner

    interface
      function pdm_dist_cloud_surf_create_cf(mesh_nature,   &
                                             n_point_cloud, &
                                             comm,          &
                                             owner)         &
        result(dcs)                                         &
        bind (c, name = 'PDM_dist_cloud_surf_create')
        use iso_c_binding
        implicit none
        integer(c_int), value :: mesh_nature
        integer(c_int), value :: n_point_cloud
        type(c_ptr), value    :: comm
        integer(c_int), value :: owner
        type (c_ptr)          :: dcs
      end function pdm_dist_cloud_surf_create_cf
    end interface

    c_comm = PDM_MPI_Comm_f2c(f_comm)

    c_mesh_nature   = mesh_nature
    c_n_point_cloud = n_point_cloud
    c_owner         = owner

    dcs = pdm_dist_cloud_surf_create_cf(c_mesh_nature,   &
                                        c_n_point_cloud, &
                                        c_comm,          &
                                        c_owner)

  end subroutine pdm_dist_cloud_surf_create



  subroutine pdm_dist_cloud_surf_n_part_cloud_set(dcs,           &
                                                  i_point_cloud, &
                                                  n_part)
    ! Set the number of partitions of a point cloud
    implicit none

    type(c_ptr), intent(in) :: dcs           ! PDM_dist_cloud_surf_t instance
    integer,     intent(in) :: i_point_cloud ! Point cloud identifier
    integer,     intent(in) :: n_part        ! Number of partitions

    interface
      subroutine pdm_dist_cloud_surf_n_part_cloud_set_cf(dcs,           &
                                                         i_point_cloud, &
                                                         n_part)        &
      bind (c, name = 'PDM_dist_cloud_surf_n_part_cloud_set')
        use iso_c_binding
        implicit none
        type(c_ptr),    value :: dcs
        integer(c_int), value :: i_point_cloud
        integer(c_int), value :: n_part
      end subroutine pdm_dist_cloud_surf_n_part_cloud_set_cf
    end interface

    call pdm_dist_cloud_surf_n_part_cloud_set_cf(dcs,           &
                                                 i_point_cloud, &
                                                 n_part)

  end subroutine pdm_dist_cloud_surf_n_part_cloud_set



  subroutine pdm_dist_cloud_surf_get(dcs,                   &
                                     i_point_cloud,         &
                                     i_part,                &
                                     closest_elt_distance,  &
                                     closest_elt_projected, &
                                     closest_elt_gnum)

    ! Get mesh distance
    implicit none

    type(c_ptr),          intent(in) :: dcs                        ! PDM_dist_cloud_surf object
    integer,              intent(in) :: i_point_cloud              ! Point cloud identifier
    integer,              intent(in) :: i_part                     ! Partition identifier
    real(8),              pointer    :: closest_elt_distance(:)    ! Distance
    real(8),              pointer    :: closest_elt_projected(:,:) ! Projected point coordinates
    integer(pdm_g_num_s), pointer    :: closest_elt_gnum(:)        ! Global ID of the closest element

    integer(c_int)                   :: c_i_point_cloud
    integer(c_int)                   :: c_i_part
    type(c_ptr)                      :: c_closest_elt_distance
    type(c_ptr)                      :: c_closest_elt_projected
    type(c_ptr)                      :: c_closest_elt_gnum
    integer                          :: n_points

    interface
      subroutine pdm_dist_cloud_surf_get_cf(dcs,                   &
                                            i_point_cloud,         &
                                            i_part,                &
                                            closest_elt_distance,  &
                                            closest_elt_projected, &
                                            closest_elt_gnum)      &
        bind (c, name = 'PDM_dist_cloud_surf_get')

        use iso_c_binding

        implicit none

        type (c_ptr), value       :: dcs
        integer(c_int), value     :: i_point_cloud
        integer(c_int), value     :: i_part
        type(c_ptr)               :: closest_elt_distance
        type(c_ptr)               :: closest_elt_projected
        type(c_ptr)               :: closest_elt_gnum

      end subroutine pdm_dist_cloud_surf_get_cf
    end interface

    c_i_point_cloud = i_point_cloud
    c_i_part        = i_part

    call PDM_dist_cloud_surf_cloud_dim_get_cf (dcs,             &
                                            c_i_point_cloud, &
                                            c_i_part,        &
                                            n_points)

    call pdm_dist_cloud_surf_get_cf (dcs,                     &
                                     c_i_point_cloud,         &
                                     c_i_part,                &
                                     c_closest_elt_distance,  &
                                     c_closest_elt_projected, &
                                     c_closest_elt_gnum)

    call c_f_pointer(c_closest_elt_distance, &
                     closest_elt_distance,   &
                     [n_points])

    call c_f_pointer(c_closest_elt_projected, &
                     closest_elt_projected,   &
                     [3,n_points])

    call c_f_pointer(c_closest_elt_gnum, &
                     closest_elt_gnum,   &
                     [n_points])

  end subroutine pdm_dist_cloud_surf_get



  subroutine pdm_dist_cloud_surf_distri_data(dcs,           &
                                             i_point_cloud, &
                                             stride,        &
                                             surf_data,     &
                                             cloud_data)
    ! Distribute data from the surface mesh to a point cloud
    implicit none

    type(c_ptr),               intent(in) :: dcs           ! PDM_dist_cloud_surf_t instance
    integer,                   intent(in) :: i_point_cloud ! Point cloud identifier
    integer,                   intent(in) :: stride        ! Stride
    type(PDM_pointer_array_t), pointer    :: surf_data     ! Surface mesh data (send)
    type(PDM_pointer_array_t), pointer    :: cloud_data    ! Point cloud data  (recv)

    integer(c_int)                        :: c_i_point_cloud
    integer(c_int)                        :: c_i_part
    integer(c_int)                        :: c_stride
    type(c_ptr)                           :: c_cloud_data
    integer                               :: n_points
    integer                               :: n_part
    integer                               :: i_part
    integer, allocatable                  :: length_data(:)

    interface
      subroutine pdm_dist_cloud_surf_distri_data_cf(dcs,           &
                                                    i_point_cloud, &
                                                    stride,        &
                                                    surf_data,     &
                                                    cloud_data)    &
        bind (c, name = 'PDM_dist_cloud_surf_distri_data')
        use iso_c_binding
        implicit none
        type (c_ptr),   value :: dcs
        integer(c_int), value :: i_point_cloud
        integer(c_int), value :: stride
        type(c_ptr),    value :: surf_data
        type(c_ptr)           :: cloud_data
      end subroutine pdm_dist_cloud_surf_distri_data_cf
    end interface

    c_i_point_cloud = i_point_cloud
    c_stride        = stride

    n_part = PDM_dist_cloud_surf_cloud_n_part_get_cf (dcs, i_point_cloud)

    call pdm_dist_cloud_surf_distri_data_cf (dcs,                    &
                                             c_i_point_cloud,        &
                                             c_stride,               &
                                             c_loc(surf_data%cptr),  &
                                             c_cloud_data)

    allocate( length_data(n_part) )

    do i_part = 1, n_part

      c_i_part = i_part-1

      call PDM_dist_cloud_surf_cloud_dim_get_cf (dcs,             &
                                              c_i_point_cloud, &
                                              c_i_part,        &
                                              n_points)
      length_data(i_part) = n_points*stride

    end do

    call PDM_pointer_array_create (cloud_data,         &
                                   n_part,             &
                                   surf_data%type,     &
                                   c_cloud_data,       &
                                   length_data,        &
                                   PDM_OWNERSHIP_KEEP, &
                                   surf_data%s_data)

    deallocate( length_data )

  end subroutine pdm_dist_cloud_surf_distri_data



  subroutine pdm_dist_cloud_surf_cloud_set(dcs,           &
                                           i_point_cloud, &
                                           i_part,        &
                                           n_points,      &
                                           coords,        &
                                           gnum)

    ! Set a point cloud
    implicit none

    type(c_ptr),               intent(in) :: dcs           ! PDM_dist_cloud_surf_t instance
    integer,                   intent(in) :: i_point_cloud ! Point cloud identifier
    integer,                   intent(in) :: i_part        ! Partition identifier
    integer,                   intent(in) :: n_points      ! Number of points
    real(8),                   pointer    :: coords(:,:)   ! Point coordinates (shape : [3, n_points])
    integer(kind=pdm_g_num_s), pointer    :: gnum(:)       ! Point global ids (shape : [n_points])

    integer(c_int) :: c_i_point_cloud
    integer(c_int) :: c_i_part
    integer(c_int) :: c_n_points

    type(c_ptr)    :: c_coords
    type(c_ptr)    :: c_gnum

    interface
      subroutine pdm_dist_cloud_surf_cloud_set_cf(dcs,           &
                                                  i_point_cloud, &
                                                  i_part,        &
                                                  n_points,      &
                                                  coords,        &
                                                  gnum)          &
        bind (c, name = 'PDM_dist_cloud_surf_cloud_set')
        use iso_c_binding
        implicit none
        type (c_ptr),   value :: dcs
        integer(c_int), value :: i_point_cloud
        integer(c_int), value :: i_part
        integer(c_int), value :: n_points
        type(c_ptr),    value :: coords
        type(c_ptr),    value :: gnum
      end subroutine pdm_dist_cloud_surf_cloud_set_cf
    end interface


    c_i_point_cloud = i_point_cloud
    c_i_part        = i_part
    c_n_points      = n_points

    c_coords = c_loc(coords)
    c_gnum   = c_loc(gnum)

    call pdm_dist_cloud_surf_cloud_set_cf(dcs,             &
                                          c_i_point_cloud, &
                                          c_i_part,        &
                                          c_n_points,      &
                                          c_coords,        &
                                          c_gnum)

  end subroutine pdm_dist_cloud_surf_cloud_set



  subroutine pdm_dist_cloud_surf_nodal_mesh_set(dcs, &
                                                mesh_nodal)

    ! Set the nodal mesh
    implicit none

    type(c_ptr), intent(in) :: dcs        ! PDM_dist_cloud_surf_t instance
    type(c_ptr), intent(in) :: mesh_nodal ! PDM_part_mesh_nodal_t instance

    interface
      subroutine pdm_dist_cloud_surf_nodal_mesh_set_cf(dcs,        &
                                                       mesh_nodal) &
        bind (c, name = 'PDM_dist_cloud_surf_nodal_mesh_set')
        use iso_c_binding
        implicit none
        type(c_ptr), value :: dcs
        type(c_ptr), value :: mesh_nodal
      end subroutine pdm_dist_cloud_surf_nodal_mesh_set_cf
    end interface

    call pdm_dist_cloud_surf_nodal_mesh_set_cf(dcs, mesh_nodal)

  end subroutine pdm_dist_cloud_surf_nodal_mesh_set



  subroutine pdm_dist_cloud_surf_surf_mesh_global_data_set(dcs, &
                                                           n_part)
    ! Set the number of partitions of the mesh
    implicit none

    type(c_ptr), intent(in) :: dcs    ! PDM_dist_cloud_surf_t instance
    integer,     intent(in) :: n_part ! Number of partitions

    interface
      subroutine pdm_dist_cloud_surf_surf_mesh_global_data_set_cf(dcs,    &
                                                                  n_part) &
        bind (c, name = 'PDM_dist_cloud_surf_surf_mesh_global_data_set')
        use iso_c_binding
        implicit none
        type(c_ptr),    value :: dcs
        integer(c_int), value :: n_part
      end subroutine pdm_dist_cloud_surf_surf_mesh_global_data_set_cf
    end interface

    call pdm_dist_cloud_surf_surf_mesh_global_data_set_cf(dcs, &
                                                          n_part)

  end subroutine pdm_dist_cloud_surf_surf_mesh_global_data_set



  subroutine pdm_dist_cloud_surf_surf_mesh_part_set(dcs,           &
                                                    i_part,        &
                                                    n_face,        &
                                                    face_vtx_idx,  &
                                                    face_vtx,      &
                                                    face_ln_to_gn, &
                                                    n_vtx,         &
                                                    coords,        &
                                                    vtx_ln_to_gn)
    ! Set a part of a surface mesh
    implicit none

    type(c_ptr),               intent(in) :: dcs              ! Pointer to \ref PDM_dist_cloud_surf object
    integer,                   intent(in) :: i_part           ! Partition to define
    integer,                   intent(in) :: n_face           ! Number of faces
    integer,                   pointer    :: face_vtx_idx(:)  ! Index in the face -> vertex connectivity
    integer,                   pointer    :: face_vtx(:)      ! face -> vertex connectivity
    integer(kind=pdm_g_num_s), pointer    :: face_ln_to_gn(:) ! Local face numbering to global face numbering
    integer,                   intent(in) :: n_vtx            ! Number of vertices
    real(8),                   pointer    :: coords(:,:)      ! Coordinates
    integer(kind=pdm_g_num_s), pointer    :: vtx_ln_to_gn(:)  ! Local vertex numbering to global vertex numbering

    integer(c_int)     :: c_i_part
    integer(c_int)     :: c_n_face
    type(c_ptr)        :: c_face_vtx_idx
    type(c_ptr)        :: c_face_vtx
    type(c_ptr)        :: c_face_ln_to_gn
    integer(c_int)     :: c_n_vtx
    type(c_ptr)        :: c_coords
    type(c_ptr)        :: c_vtx_ln_to_gn

    interface
      subroutine pdm_dist_cloud_surf_surf_mesh_part_set_cf(dcs,           &
                                                           i_part,        &
                                                           n_face,        &
                                                           face_vtx_idx,  &
                                                           face_vtx,      &
                                                           face_ln_to_gn, &
                                                           n_vtx,         &
                                                           coords,        &
                                                           vtx_ln_to_gn)  &
        bind (c, name = 'PDM_dist_cloud_surf_surf_mesh_part_set')
        use iso_c_binding

        implicit none

        type (c_ptr),   value  :: dcs
        integer(c_int), value  :: i_part
        integer(c_int), value  :: n_face
        type(c_ptr),    value  :: face_vtx_idx
        type(c_ptr),    value  :: face_vtx
        type(c_ptr),    value  :: face_ln_to_gn
        integer(c_int), value  :: n_vtx
        type(c_ptr),    value  :: coords
        type(c_ptr),    value  :: vtx_ln_to_gn

      end subroutine pdm_dist_cloud_surf_surf_mesh_part_set_cf
    end interface


    c_i_part = i_part
    c_n_face = n_face
    c_n_vtx  = n_vtx

    c_face_vtx_idx  = c_loc(face_vtx_idx)
    c_face_vtx      = c_loc(face_vtx)
    c_face_ln_to_gn = c_loc(face_ln_to_gn)
    c_coords        = c_loc(coords)
    c_vtx_ln_to_gn  = c_loc(vtx_ln_to_gn)

    call pdm_dist_cloud_surf_surf_mesh_part_set_cf(dcs,             &
                                                   c_i_part,        &
                                                   c_n_face,        &
                                                   c_face_vtx_idx,  &
                                                   c_face_vtx,      &
                                                   c_face_ln_to_gn, &
                                                   c_n_vtx,         &
                                                   c_coords,        &
                                                   c_vtx_ln_to_gn)

  end subroutine pdm_dist_cloud_surf_surf_mesh_part_set



  subroutine pdm_dist_cloud_surf_compute(dcs)
    ! Compute distance
    implicit none

    type(c_ptr), intent(in) :: dcs ! PDM_dist_cloud_surf_t instance

    interface
      subroutine pdm_dist_cloud_surf_compute_cf(dcs) &
        bind (c, name = 'PDM_dist_cloud_surf_compute')
        use iso_c_binding
        implicit none
        type (c_ptr), value :: dcs
      end subroutine pdm_dist_cloud_surf_compute_cf
    end interface

    call pdm_dist_cloud_surf_compute_cf(dcs)

  end subroutine pdm_dist_cloud_surf_compute



  subroutine pdm_dist_cloud_surf_dump_times(dcs)
    ! Dump elapsed and CPU times
    implicit none

    type(c_ptr), intent(in) :: dcs

    interface
      subroutine pdm_dist_cloud_surf_dump_times_cf(dcs) &
        bind (c, name = 'PDM_dist_cloud_surf_dump_times')
        use iso_c_binding
        implicit none
        type(c_ptr), value :: dcs
      end subroutine pdm_dist_cloud_surf_dump_times_cf
    end interface

    call pdm_dist_cloud_surf_dump_times_cf(dcs)

  end subroutine pdm_dist_cloud_surf_dump_times



  subroutine pdm_dist_cloud_surf_free(dcs)
    ! Free a distance mesh structure
    implicit none

    type(c_ptr), intent(inout) :: dcs ! PDM_dist_cloud_surf_t instance

    interface
      subroutine pdm_dist_cloud_surf_free_cf(dcs) &
        bind (c, name = 'PDM_dist_cloud_surf_free')
        use iso_c_binding
        implicit none
        type (c_ptr), value :: dcs
      end subroutine pdm_dist_cloud_surf_free_cf
    end interface

    call pdm_dist_cloud_surf_free_cf(dcs)

  end subroutine pdm_dist_cloud_surf_free

end module pdm_dist_cloud_surf
