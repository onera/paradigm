#include "pdm_configf.h"

module pdm_mesh_location

  use pdm
  use iso_c_binding

  implicit none

  !!
  !! Enum type PDM_mesh_location_method_t
  !!

  integer(c_int), parameter :: PDM_MESH_LOCATION_OCTREE         = 0 ! Use point octree
  integer(c_int), parameter :: PDM_MESH_LOCATION_DBBTREE        = 1 ! Use bounding-box tree
  integer(c_int), parameter :: PDM_MESH_LOCATION_LOCATE_ALL_TGT = 2 ! Locate all target points
  ! integer(c_int), parameter :: PDM_MESH_LOCATION_DOCTREE        = 3 !


  interface PDM_mesh_location_cell_vertex_get;
    module procedure PDM_mesh_location_cell_vertex_get_cptr
    module procedure PDM_mesh_location_cell_vertex_get_f
  end interface

  interface

    !>
    !!
    !! \brief Get a point cloud
    !!
    !! \param [in]   mloc            Pointer to \ref PDM_mesh_location object
    !! \param [in]   i_point_cloud   Index of point cloud
    !! \param [in]   i_part          Index of partition
    !! \param [out]  n_points        Number of points
    !! \param [out]  coords          Point coordinates
    !! \param [out]  gnum            Point global number
    !!
    !!

    subroutine PDM_mesh_location_cloud_get (mloc, &
                                            i_point_cloud, &
                                            i_part, &
                                            n_points, &
                                            coords, &
                                            gnum) &
     bind (c, name = 'PDM_mesh_location_cloud_get')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc
      integer(c_int), value :: i_point_cloud
      integer(c_int), value :: i_part
      integer(c_int), value :: n_points
      type(c_ptr)           :: coords
      type(c_ptr)           :: gnum

    end subroutine PDM_mesh_location_cloud_get


    !>
    !!
    !! \brief get cell vertex connectivity
    !!
    !! \param [in]   mloc                  Pointer to \ref PDM_mesh_location object
    !! \param [in]   i_part                Index of partition of the cloud
    !! \param [out]  cell_vtx_idx          Index in (size = n_elt + 1)
    !! \param [out]  cell_vtx              Cell vertex connectivity
    !!
    !!

    subroutine PDM_mesh_location_cell_vertex_get_cf(mloc, &
                                                    i_part, &
                                                    cell_vtx_idx, &
                                                    cell_vtx) &
     bind (c, name = 'PDM_mesh_location_cell_vertex_get')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc
      integer(c_int), value :: i_part
      type(c_ptr)           :: cell_vtx_idx
      type(c_ptr)           :: cell_vtx

    end subroutine PDM_mesh_location_cell_vertex_get_cf

    
    function PDM_mesh_location_n_located_get_cf(mloc, &
                                                i_point_cloud, &
                                                i_part) &
                                                result (n_located) &
      bind (c, name = 'PDM_mesh_location_n_located_get')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: mloc
      integer(c_int), value :: i_point_cloud
      integer(c_int), value :: i_part
      integer(c_int)        :: n_located

    end function PDM_mesh_location_n_located_get_cf


    !>
    !!
    !! \brief Get the number of cells
    !!
    !! \param [in]  mloc     Pointer to \ref PDM_mesh_location object
    !! \param [in]  i_part   Index of partition of the mesh
    !!
    !! \return Number of cells
    !!

    function PDM_mesh_location_n_cell_get (mloc, &
                                           i_part) &
                                           result(n_cell) &
      bind (c, name = 'PDM_mesh_location_n_cell_get')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc
      integer(c_int), value :: i_part
      integer(c_int)        :: n_cell

    end function PDM_mesh_location_n_cell_get


    function PDM_mesh_location_n_unlocated_get_cf(mloc, &
                                                  i_point_cloud, &
                                                  i_part) &
                                                  result(n_unlocated) &
      bind (c, name = 'PDM_mesh_location_n_unlocated_get')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: mloc
      integer(c_int), value :: i_point_cloud
      integer(c_int), value :: i_part
      integer(c_int)        :: n_unlocated

    end function PDM_mesh_location_n_unlocated_get_cf


    function PDM_mesh_location_mesh_nodal_get (mloc) &
                                                 result(mesh_nodal) &
      bind (c, name = 'PDM_mesh_location_mesh_nodal_get')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc
      type (c_ptr) :: mesh_nodal

    end function PDM_mesh_location_mesh_nodal_get


  end interface


  contains


  subroutine PDM_mesh_location_create(mloc,          &
                                      n_point_cloud, &
                                      f_comm,        &
                                      owner)
  ! Create a structure to compute the location of point clouds inside a mesh
  implicit none

  type(c_ptr), intent(out) :: mloc          ! C pointer to PDM_mesh_location_t object
  integer,     intent(in)  :: n_point_cloud ! Number of point clouds
  integer,     intent(in)  :: f_comm        ! Fortran MPI communicator
  integer,     intent(in)  :: owner         ! Ownership
  type(c_ptr)              :: c_comm

  interface
    function PDM_mesh_location_create_cf (n_point_cloud, &
                                          comm,          &
                                          owner ) &
                                          result(mloc) &
      bind (c, name = 'PDM_mesh_location_create')
      use iso_c_binding
      implicit none
      integer(c_int), value :: n_point_cloud
      type(c_ptr), value    :: comm
      integer(c_int), value :: owner
      type(c_ptr)           :: mloc
    end function PDM_mesh_location_create_cf
  end interface

  c_comm = PDM_MPI_Comm_f2c(f_comm)

  mloc = PDM_mesh_location_create_cf(n_point_cloud, &
                                     c_comm,        &
                                     owner)

  end subroutine PDM_mesh_location_create



  subroutine PDM_mesh_location_n_part_cloud_set(mloc, &
                                                i_point_cloud, &
                                                n_part)
    ! Set the number of partitions of a point cloud
    implicit none

    type (c_ptr), intent(in) :: mloc          ! C pointer to PDM_mesh_location_t object
    integer,      intent(in) :: i_point_cloud ! Point cloud identifier
    integer,      intent(in) :: n_part        ! Number of partitions

    interface
      subroutine PDM_mesh_location_n_part_cloud_set_cf(mloc, &
                                                       i_point_cloud, &
                                                       n_part) &
      bind (c, name = 'PDM_mesh_location_n_part_cloud_set')
        use iso_c_binding
        implicit none
        type (c_ptr),   value :: mloc
        integer(c_int), value :: i_point_cloud
        integer(c_int), value :: n_part
      end subroutine PDM_mesh_location_n_part_cloud_set_cf
    end interface

    call PDM_mesh_location_n_part_cloud_set_cf(mloc, &
                                               i_point_cloud, &
                                               n_part)

  end subroutine PDM_mesh_location_n_part_cloud_set



  subroutine PDM_mesh_location_cloud_set(mloc, &
                                         i_point_cloud, &
                                         i_part, &
                                         n_points, &
                                         coords, &
                                         gnum)
    ! Set a point cloud
    implicit none

    type (c_ptr), intent(in)           :: mloc          ! C pointer to PDM_mesh_location_t object
    integer, intent(in)                :: i_point_cloud ! Point cloud identifier
    integer, intent(in)                :: i_part        ! Partition identifier
    integer, intent(in)                :: n_points      ! Number of points
    real(8),                   pointer :: coords(:,:)   ! Point coordinates (shape = [3, ``n_points``])
    integer(kind=pdm_g_num_s), pointer :: gnum(:)       ! Point global ids (size = ``n_points``)

    type(c_ptr)                        :: c_coords
    type(c_ptr)                        :: c_gnum

    interface
      subroutine PDM_mesh_location_cloud_set_cf(mloc, &
                                                i_point_cloud, &
                                                i_part, &
                                                n_points, &
                                                coords, &
                                                gnum) &
      bind (c, name = 'PDM_mesh_location_cloud_set')
        use iso_c_binding
        implicit none
        type (c_ptr),   value :: mloc
        integer(c_int), value :: i_point_cloud
        integer(c_int), value :: i_part
        integer(c_int), value :: n_points
        type(c_ptr),    value :: coords
        type(c_ptr),    value :: gnum
      end subroutine PDM_mesh_location_cloud_set_cf
    end interface

    c_coords = C_NULL_PTR
    if (associated(coords)) then
      c_coords = c_loc(coords)
    endif
      
    c_gnum = C_NULL_PTR  
    if (associated(gnum)) then
      c_gnum = c_loc(gnum)
    endif  

    call PDM_mesh_location_cloud_set_cf(mloc,          &
                                        i_point_cloud, &
                                        i_part,        &
                                        n_points,      &
                                        c_coords,      &
                                        c_gnum)

  end subroutine PDM_mesh_location_cloud_set



  subroutine PDM_mesh_location_shared_nodal_mesh_set(mloc, &
                                                     mesh_nodal_id)
    ! Set the Part Mesh Nodal instance
    implicit none

    type (c_ptr), intent(in) :: mloc          ! Mesh location instance
    type (c_ptr), intent(in) :: mesh_nodal_id ! Part Mesh Nodal instance

    interface
      subroutine PDM_mesh_location_shared_nodal_mesh_set_cf(mloc, &
                                                            mesh_nodal_id) &
        bind (c, name = 'PDM_mesh_location_shared_nodal_mesh_set')
        use iso_c_binding
        implicit none
        type (c_ptr), value :: mloc
        type (c_ptr), value :: mesh_nodal_id
      end subroutine PDM_mesh_location_shared_nodal_mesh_set_cf
    end interface

    call PDM_mesh_location_shared_nodal_mesh_set_cf(mloc, &
                                                    mesh_nodal_id)

  end subroutine PDM_mesh_location_shared_nodal_mesh_set



  subroutine PDM_mesh_location_mesh_n_part_set(mloc, &
                                               n_part)
    ! Set the number of partitions of the source mesh
    implicit none

    type (c_ptr), intent(in) :: mloc   ! C pointer to PDM_mesh_location_t object
    integer,      intent(in) :: n_part ! Number of partitions

    interface
      subroutine PDM_mesh_location_mesh_n_part_set_cf(mloc, &
                                                      n_part) &
      bind (c, name = 'PDM_mesh_location_mesh_n_part_set')
        use iso_c_binding
        implicit none
        type (c_ptr), value :: mloc
        integer(c_int), value :: n_part
      end subroutine PDM_mesh_location_mesh_n_part_set_cf
    end interface

    call PDM_mesh_location_mesh_n_part_set_cf(mloc, &
                                              n_part)

  end subroutine PDM_mesh_location_mesh_n_part_set



  subroutine PDM_mesh_location_part_set(mloc,          &
                                        i_part,        &
                                        n_cell,        &
                                        cell_face_idx, &
                                        cell_face,     &
                                        cell_ln_to_gn, &
                                        n_face,        &
                                        face_vtx_idx,  &
                                        face_vtx,      &
                                        face_ln_to_gn, &
                                        n_vtx,         &
                                        coords,        &
                                        vtx_ln_to_gn)
    ! Set a *volume* mesh partition
    implicit none

    type (c_ptr), value                :: mloc             ! C pointer to PDM_mesh_location_t object
    integer, intent(in)                :: i_part           ! Partition identifier
    integer, intent(in)                :: n_cell           ! Number of cells
    integer(kind=pdm_l_num_s), pointer :: cell_face_idx(:) ! Index for cell -> face connectivity (size = ``n_cell`` + 1)
    integer(kind=pdm_l_num_s), pointer :: cell_face(:)     ! Cell -> face connectivity (size = ``cell_face_idx(n_cell+1)``)
    integer(kind=pdm_g_num_s), pointer :: cell_ln_to_gn(:) ! Cell global ids (size = ``n_cell``)
    integer, intent(in)                :: n_face           ! Number of faces
    integer(kind=pdm_l_num_s), pointer :: face_vtx_idx(:)  ! Index for face -> vertex connectivity (size = ``n_face`` + 1)
    integer(kind=pdm_l_num_s), pointer :: face_vtx(:)      ! Face -> vertex connectivity (size = ``face_vtx_idx(n_face+1)``)
    integer(kind=pdm_g_num_s), pointer :: face_ln_to_gn(:) ! Face global ids (size = ``n_face``)
    integer, intent(in)                :: n_vtx            ! Number of vertices
    real(8),                   pointer :: coords(:,:)      ! Vertex coordinates (shape = [3, ``n_vtx``])
    integer(kind=pdm_g_num_s), pointer :: vtx_ln_to_gn(:)  ! Vertex global ids (size = ``n_vtx``)

    type(c_ptr)                        :: c_cell_face_idx
    type(c_ptr)                        :: c_cell_face
    type(c_ptr)                        :: c_cell_ln_to_gn
    type(c_ptr)                        :: c_face_vtx_idx
    type(c_ptr)                        :: c_face_vtx
    type(c_ptr)                        :: c_face_ln_to_gn
    type(c_ptr)                        :: c_coords
    type(c_ptr)                        :: c_vtx_ln_to_gn

    interface
      subroutine PDM_mesh_location_part_set_cf(mloc, &
                                               i_part, &
                                               n_cell, &
                                               cell_face_idx, &
                                               cell_face, &
                                               cell_ln_to_gn, &
                                               n_face, &
                                               face_vtx_idx, &
                                               face_vtx, &
                                               face_ln_to_gn, &
                                               n_vtx, &
                                               coords, &
                                               vtx_ln_to_gn) &
      bind (c, name = 'PDM_mesh_location_part_set')
        use iso_c_binding
        implicit none
        type (c_ptr),   value :: mloc
        integer(c_int), value :: i_part
        integer(c_int), value :: n_cell
        type(c_ptr),    value :: cell_face_idx
        type(c_ptr),    value :: cell_face
        type(c_ptr),    value :: cell_ln_to_gn
        integer(c_int), value :: n_face
        type(c_ptr),    value :: face_vtx_idx
        type(c_ptr),    value :: face_vtx
        type(c_ptr),    value :: face_ln_to_gn
        integer(c_int), value :: n_vtx
        type(c_ptr),    value :: coords
        type(c_ptr),    value :: vtx_ln_to_gn
      end subroutine PDM_mesh_location_part_set_cf
    end interface

    c_cell_face_idx = C_NULL_PTR
    if (associated(cell_face_idx)) then
      c_cell_face_idx = c_loc(cell_face_idx)
    endif
      
    c_cell_face = C_NULL_PTR
    if (associated(cell_face)) then
      c_cell_face = c_loc(cell_face)
    endif
      
    c_cell_ln_to_gn = C_NULL_PTR
    if (associated(cell_ln_to_gn)) then
      c_cell_ln_to_gn = c_loc(cell_ln_to_gn)
    endif
      
    c_face_vtx_idx = C_NULL_PTR
    if (associated(face_vtx_idx)) then
      c_face_vtx_idx = c_loc(face_vtx_idx)
    endif
      
    c_face_vtx = C_NULL_PTR
    if (associated(face_vtx)) then
      c_face_vtx = c_loc(face_vtx)
    endif
      
    c_face_ln_to_gn = C_NULL_PTR
    if (associated(face_ln_to_gn)) then
      c_face_ln_to_gn = c_loc(face_ln_to_gn)
    endif
      
    c_coords = C_NULL_PTR
    if (associated(coords)) then
      c_coords = c_loc(coords)
    endif
      
    c_vtx_ln_to_gn = C_NULL_PTR
    if (associated(vtx_ln_to_gn)) then
      c_vtx_ln_to_gn = c_loc(vtx_ln_to_gn)
    endif
      

    call PDM_mesh_location_part_set_cf(mloc,            &
                                       i_part,          &
                                       n_cell,          &
                                       c_cell_face_idx, &
                                       c_cell_face,     &
                                       c_cell_ln_to_gn, &
                                       n_face,          &
                                       c_face_vtx_idx,  &
                                       c_face_vtx,      &
                                       c_face_ln_to_gn, &
                                       n_vtx,           &
                                       c_coords,        &
                                       c_vtx_ln_to_gn)

  end subroutine PDM_mesh_location_part_set



  subroutine PDM_mesh_location_nodal_part_set(mloc,          &
                                              i_part,        &
                                              n_cell,        &
                                              cell_vtx_idx,  &
                                              cell_vtx,      &
                                              cell_ln_to_gn, &
                                              n_vtx,         &
                                              coords,        &
                                              vtx_ln_to_gn)
    ! Set a *volume* mesh partition defined by nodal connectivity
    !
    ! The mesh is assumed to contain only standard elements
    ! (tetrahedra, pyramids, prisms, hexahedra).
    implicit none

    type (c_ptr), value                :: mloc             ! C pointer to PDM_mesh_location_t object
    integer, intent(in)                :: i_part           ! Partition identifier
    integer, intent(in)                :: n_cell           ! Number of cells
    integer(kind=pdm_l_num_s), pointer :: cell_vtx_idx(:)  ! Index for cell -> face connectivity (size = ``n_cell`` + 1)
    integer(kind=pdm_l_num_s), pointer :: cell_vtx(:)      ! Cell -> face connectivity (size = ``cell_face_idx(n_cell+1)``)
    integer(kind=pdm_g_num_s), pointer :: cell_ln_to_gn(:) ! Cell global ids (size = ``n_cell``)
    integer, intent(in)                :: n_vtx            ! Number of vertices
    real(8),                   pointer :: coords(:,:)      ! Vertex coordinates (shape = [3, ``n_vtx``])
    integer(kind=pdm_g_num_s), pointer :: vtx_ln_to_gn(:)  ! Vertex global ids (size = ``n_vtx``)

    type(c_ptr)                        :: c_cell_vtx_idx
    type(c_ptr)                        :: c_cell_vtx
    type(c_ptr)                        :: c_cell_ln_to_gn
    type(c_ptr)                        :: c_coords
    type(c_ptr)                        :: c_vtx_ln_to_gn

    interface
      subroutine PDM_mesh_location_nodal_part_set_cf(mloc, &
                                                     i_part, &
                                                     n_cell, &
                                                     cell_vtx_idx, &
                                                     cell_vtx, &
                                                     cell_ln_to_gn, &
                                                     n_vtx, &
                                                     coords, &
                                                     vtx_ln_to_gn) &
      bind (c, name = 'PDM_mesh_location_nodal_part_set')
        use iso_c_binding
        implicit none
        type (c_ptr),   value :: mloc
        integer(c_int), value :: i_part
        integer(c_int), value :: n_cell
        type(c_ptr),    value :: cell_vtx_idx
        type(c_ptr),    value :: cell_vtx
        type(c_ptr),    value :: cell_ln_to_gn
        integer(c_int), value :: n_vtx
        type(c_ptr),    value :: coords
        type(c_ptr),    value :: vtx_ln_to_gn
      end subroutine PDM_mesh_location_nodal_part_set_cf
    end interface

    c_cell_vtx_idx = C_NULL_PTR
    if (associated(cell_vtx_idx)) then
      c_cell_vtx_idx = c_loc(cell_vtx_idx)
    endif
      
    c_cell_vtx = C_NULL_PTR
    if (associated(cell_vtx)) then
      c_cell_vtx = c_loc(cell_vtx)
    endif
      
    c_cell_ln_to_gn = C_NULL_PTR
    if (associated(cell_ln_to_gn)) then
      c_cell_ln_to_gn = c_loc(cell_ln_to_gn)
    endif
      
    c_coords = C_NULL_PTR
    if (associated(coords)) then
      c_coords = c_loc(coords)
    endif
      
    c_vtx_ln_to_gn = C_NULL_PTR
    if (associated(vtx_ln_to_gn)) then
      c_vtx_ln_to_gn = c_loc(vtx_ln_to_gn)
    endif    

    call PDM_mesh_location_nodal_part_set_cf(mloc,            &
                                             i_part,          &
                                             n_cell,          &
                                             c_cell_vtx_idx,  &
                                             c_cell_vtx,      &
                                             c_cell_ln_to_gn, &
                                             n_vtx,           &
                                             c_coords,        &
                                             c_vtx_ln_to_gn)

  end subroutine PDM_mesh_location_nodal_part_set



  subroutine PDM_mesh_location_part_set_2d(mloc,          &
                                           i_part,        &
                                           n_face,        &
                                           face_edge_idx, &
                                           face_edge,     &
                                           face_ln_to_gn, &
                                           n_edge,        &
                                           edge_vtx,      &
                                           n_vtx,         &
                                           coords,        &
                                           vtx_ln_to_gn)
    ! Set a *surface* mesh partition
    implicit none

    type (c_ptr), value                :: mloc             ! C pointer to PDM_mesh_location_t object
    integer, intent(in)                :: i_part           ! Partition identifier
    integer, intent(in)                :: n_face           ! Number of faces
    integer(kind=pdm_l_num_s), pointer :: face_edge_idx(:) ! Index for face -> edge connectivity (size = ``n_face`` + 1)
    integer(kind=pdm_l_num_s), pointer :: face_edge(:)     ! Face -> edge connectivity (size = ``face_edge_idx(n_face+1)``)
    integer(kind=pdm_g_num_s), pointer :: face_ln_to_gn(:) ! Face global ids (size = ``n_face``)
    integer, intent(in)                :: n_edge           ! Number of edges
    integer(kind=pdm_l_num_s), pointer :: edge_vtx(:)      ! Edge -> vertex connectivity (size = 2 * ``n_edge``)
    integer, intent(in)                :: n_vtx            ! Number of vertices
    real(8),                   pointer :: coords(:,:)      ! Vertex coordinates (shape = [3, ``n_vtx``])
    integer(kind=pdm_g_num_s), pointer :: vtx_ln_to_gn(:)  ! Vertex global ids (size = ``n_vtx``)

    type(c_ptr)                        :: c_face_edge_idx
    type(c_ptr)                        :: c_face_edge
    type(c_ptr)                        :: c_face_ln_to_gn
    type(c_ptr)                        :: c_edge_vtx
    type(c_ptr)                        :: c_coords
    type(c_ptr)                        :: c_vtx_ln_to_gn

    interface
      subroutine PDM_mesh_location_part_set_2d_cf(mloc, &
                                                  i_part, &
                                                  n_face, &
                                                  face_edge_idx, &
                                                  face_edge, &
                                                  face_ln_to_gn, &
                                                  n_edge, &
                                                  edge_vtx, &
                                                  n_vtx, &
                                                  coords, &
                                                  vtx_ln_to_gn) &
      bind (c, name = 'PDM_mesh_location_part_set_2d')
        use iso_c_binding
        implicit none
        type (c_ptr),   value :: mloc
        integer(c_int), value :: i_part
        integer(c_int), value :: n_face
        type(c_ptr),    value :: face_edge_idx
        type(c_ptr),    value :: face_edge
        type(c_ptr),    value :: face_ln_to_gn
        integer(c_int), value :: n_edge
        type(c_ptr),    value :: edge_vtx
        integer(c_int), value :: n_vtx
        type(c_ptr),    value :: coords
        type(c_ptr),    value :: vtx_ln_to_gn
      end subroutine PDM_mesh_location_part_set_2d_cf
    end interface

    c_face_edge_idx = C_NULL_PTR
    if (associated (face_edge_idx)) then
      c_face_edge_idx = c_loc(face_edge_idx)
    endif
      
    c_face_edge = C_NULL_PTR
    if (associated (face_edge)) then
      c_face_edge = c_loc(face_edge)
    endif
      
    c_face_ln_to_gn = C_NULL_PTR
    if (associated (face_ln_to_gn)) then
      c_face_ln_to_gn = c_loc(face_ln_to_gn)
    endif
      
    c_edge_vtx = C_NULL_PTR
    if (associated (edge_vtx)) then
      c_edge_vtx = c_loc(edge_vtx)
    endif
      
    c_coords = C_NULL_PTR
    if (associated (coords)) then
      c_coords = c_loc(coords)
    endif
      
    c_vtx_ln_to_gn = C_NULL_PTR
    if (associated (vtx_ln_to_gn)) then
      c_vtx_ln_to_gn = c_loc(vtx_ln_to_gn)
    endif   

    call PDM_mesh_location_part_set_2d_cf(mloc, &
                                          i_part, &
                                          n_face, &
                                          c_face_edge_idx, &
                                          c_face_edge, &
                                          c_face_ln_to_gn, &
                                          n_edge, &
                                          c_edge_vtx, &
                                          n_vtx, &
                                          c_coords, &
                                          c_vtx_ln_to_gn)

  end subroutine PDM_mesh_location_part_set_2d



  subroutine PDM_mesh_location_nodal_part_set_2d(mloc,          &
                                                 i_part,        &
                                                 n_face,        &
                                                 face_vtx_idx,  &
                                                 face_vtx,      &
                                                 face_ln_to_gn, &
                                                 n_vtx,         &
                                                 coords,        &
                                                 vtx_ln_to_gn)
    ! Set a *surface* mesh partition defined by nodal connectivity
    implicit none

    type (c_ptr), value                :: mloc             ! C pointer to PDM_mesh_location_t object
    integer, intent(in)                :: i_part           ! Partition identifier
    integer, intent(in)                :: n_face           ! Number of faces
    integer(kind=pdm_l_num_s), pointer :: face_vtx_idx(:)  ! Index for face -> vertex connectivity (size = ``n_face`` + 1)
    integer(kind=pdm_l_num_s), pointer :: face_vtx(:)      ! Face -> vertex connectivity (size = ``face_vtx_idx(n_face+1)``)
    integer(kind=pdm_g_num_s), pointer :: face_ln_to_gn(:) ! Face global ids (size = ``n_face``)
    integer, intent(in)                :: n_vtx            ! Number of vertices
    double precision,          pointer :: coords(:,:)      ! Vertex coordinates (shape = [3, ``n_vtx``])
    integer(kind=pdm_g_num_s), pointer :: vtx_ln_to_gn(:)  ! Vertex global ids (size = ``n_vtx``)

    type(c_ptr)                        :: c_face_vtx_idx
    type(c_ptr)                        :: c_face_vtx
    type(c_ptr)                        :: c_face_ln_to_gn
    type(c_ptr)                        :: c_coords
    type(c_ptr)                        :: c_vtx_ln_to_gn

    interface
      subroutine PDM_mesh_location_nodal_part_set_2d_cf(mloc, &
                                                        i_part, &
                                                        n_face, &
                                                        face_vtx_idx, &
                                                        face_vtx, &
                                                        face_ln_to_gn, &
                                                        n_vtx, &
                                                        coords, &
                                                        vtx_ln_to_gn) &
      bind (c, name = 'PDM_mesh_location_nodal_part_set_2d')
        use iso_c_binding
        implicit none
        type (c_ptr),   value :: mloc
        integer(c_int), value :: i_part
        integer(c_int), value :: n_face
        type(c_ptr),    value :: face_vtx_idx
        type(c_ptr),    value :: face_vtx
        type(c_ptr),    value :: face_ln_to_gn
        integer(c_int), value :: n_vtx
        type(c_ptr),    value :: coords
        type(c_ptr),    value :: vtx_ln_to_gn
      end subroutine PDM_mesh_location_nodal_part_set_2d_cf
    end interface

    c_face_vtx_idx = C_NULL_PTR
    if (associated(face_vtx_idx)) then
      c_face_vtx_idx = c_loc(face_vtx_idx)
    endif
      
    c_face_vtx = C_NULL_PTR
    if (associated(face_vtx)) then
      c_face_vtx = c_loc(face_vtx)
    endif
      
    c_face_ln_to_gn = C_NULL_PTR
    if (associated(face_ln_to_gn)) then
      c_face_ln_to_gn = c_loc(face_ln_to_gn)
    endif
      
    c_coords = C_NULL_PTR
    if (associated(coords)) then
      c_coords = c_loc(coords)
    endif
      
    c_vtx_ln_to_gn = C_NULL_PTR
    if (associated(vtx_ln_to_gn)) then
      c_vtx_ln_to_gn = c_loc(vtx_ln_to_gn)
    endif
      

    call PDM_mesh_location_nodal_part_set_2d_cf(mloc, &
                                                i_part, &
                                                n_face, &
                                                c_face_vtx_idx, &
                                                c_face_vtx, &
                                                c_face_ln_to_gn, &
                                                n_vtx, &
                                                c_coords, &
                                                c_vtx_ln_to_gn)

  end subroutine PDM_mesh_location_nodal_part_set_2d



  subroutine PDM_mesh_location_tolerance_set(mloc, &
                                             tol)
    ! Set the tolerance for bounding boxes
    implicit none

    type (c_ptr), intent(in) :: mloc ! C pointer to PDM_mesh_location_t object
    real(8),      intent(in) :: tol  ! Relative tolerance

    interface
      subroutine PDM_mesh_location_tolerance_set_cf(mloc, &
                                                    tol) &
      bind (c, name = 'PDM_mesh_location_tolerance_set')
        use iso_c_binding
        implicit none
        type (c_ptr),   value :: mloc
        real(c_double), value :: tol
      end subroutine PDM_mesh_location_tolerance_set_cf
    end interface

    call PDM_mesh_location_tolerance_set_cf(mloc, &
                                            tol)

  end subroutine PDM_mesh_location_tolerance_set



  subroutine PDM_mesh_location_method_set(mloc, &
                                          method)
    ! Set the method for computing location (preconditioning stage)
    !
    ! Admissible values are :
    !   - ``PDM_MESH_LOCATION_OCTREE`` : Use point octree (default method)
    !   - ``PDM_MESH_LOCATION_DBBTREE`` : Use bounding-box tree
    !   - ``PDM_MESH_LOCATION_LOCATE_ALL_TGT`` : All target points are guaranteed to be located

    implicit none

    type (c_ptr), intent(in) :: mloc   ! C pointer to PDM_mesh_location_t object
    integer,      intent(in) :: method ! Preconditioning method

    interface
      subroutine PDM_mesh_location_method_set_cf(mloc, &
                                                 method) &
      bind (c, name = 'PDM_mesh_location_method_set')
        use iso_c_binding
        implicit none
        type (c_ptr),   value :: mloc
        integer(c_int), value :: method
      end subroutine PDM_mesh_location_method_set_cf
    end interface

    call PDM_mesh_location_method_set_cf(mloc, &
                                         method)

  end subroutine PDM_mesh_location_method_set



  subroutine PDM_mesh_location_compute(mloc)
    ! Compute point location
    implicit none

    type (c_ptr), intent(in) :: mloc ! C pointer to PDM_mesh_location_t object

    interface
      subroutine PDM_mesh_location_compute_cf(mloc) &
        bind (c, name = 'PDM_mesh_location_compute')
        use iso_c_binding
        implicit none
        type (c_ptr), value :: mloc
      end subroutine PDM_mesh_location_compute_cf
    end interface

    call PDM_mesh_location_compute_cf(mloc)

  end subroutine PDM_mesh_location_compute



  subroutine PDM_mesh_location_dump_times(mloc)
    ! Dump elapsed and CPU times
    implicit none

    type (c_ptr), intent(in) :: mloc ! C pointer to PDM_mesh_location_t object

    interface
      subroutine PDM_mesh_location_dump_times_cf(mloc) &
        bind (c, name = 'PDM_mesh_location_dump_times')
        use iso_c_binding
        implicit none
        type (c_ptr), value :: mloc
      end subroutine PDM_mesh_location_dump_times_cf
    end interface

    call PDM_mesh_location_dump_times_cf(mloc)

  end subroutine PDM_mesh_location_dump_times



  function PDM_mesh_location_n_located_get(mloc, &
                                           i_point_cloud, &
                                           i_part) &
                                           result(n_located)
    ! Get the number of located points
    implicit none

    type (c_ptr), intent(in) :: mloc          ! C pointer to PDM_mesh_location_t object
    integer,      intent(in) :: i_point_cloud ! Point cloud identifier
    integer,      intent(in) :: i_part        ! Partition identifier
    integer                  :: n_located     ! Number of located points

    n_located = PDM_mesh_location_n_located_get_cf(mloc, &
                                                   i_point_cloud, &
                                                   i_part)

  end function PDM_mesh_location_n_located_get



  subroutine PDM_mesh_location_located_get(mloc,          &
                                           i_point_cloud, &
                                           i_part,        &
                                           located)
    ! Get the list of located points
    implicit none

    type (c_ptr), value :: mloc          ! C pointer to PDM_mesh_location_t object
    integer, intent(in) :: i_point_cloud ! Point cloud identifier
    integer, intent(in) :: i_part        ! Partition identifier
    integer, pointer    :: located(:)    ! List of located points

    type(c_ptr)         :: c_located
    integer(c_int)      :: n_located

    interface
      function PDM_mesh_location_located_get_cf(mloc, &
                                                i_point_cloud, &
                                                i_part) &
      result(located) &
        bind (c, name = 'PDM_mesh_location_located_get')
        use iso_c_binding
        implicit none
        type (c_ptr),   value :: mloc
        integer(c_int), value :: i_point_cloud
        integer(c_int), value :: i_part
        type(c_ptr)           :: located
      end function PDM_mesh_location_located_get_cf
    end interface

    n_located = PDM_mesh_location_n_located_get_cf(mloc,          &
                                                   i_point_cloud, &
                                                   i_part)

    c_located = PDM_mesh_location_located_get_cf(mloc,          &
                                                 i_point_cloud, &
                                                 i_part)

    call c_f_pointer(c_located,   &
                     located,     &
                     [n_located])

  end subroutine PDM_mesh_location_located_get



  function PDM_mesh_location_n_unlocated_get(mloc, &
                                             i_point_cloud, &
                                             i_part) &
                                             result(n_unlocated)
    ! Get the number of unlocated points
    implicit none

    type (c_ptr), intent(in) :: mloc          ! C pointer to PDM_mesh_location_t object
    integer,      intent(in) :: i_point_cloud ! Point cloud identifier
    integer,      intent(in) :: i_part        ! Partition identifier
    integer                  :: n_unlocated   ! Number of unlocated points

    n_unlocated = PDM_mesh_location_n_unlocated_get_cf(mloc, &
                                                       i_point_cloud, &
                                                       i_part)

  end function PDM_mesh_location_n_unlocated_get



  subroutine PDM_mesh_location_unlocated_get(mloc,          &
                                             i_point_cloud, &
                                             i_part,        &
                                             unlocated)
    ! Get the list of unlocated points
    implicit none

    type (c_ptr), value :: mloc          ! C pointer to PDM_mesh_location_t object
    integer, intent(in) :: i_point_cloud ! Point cloud identifier
    integer, intent(in) :: i_part        ! Partition identifier
    integer, pointer    :: unlocated(:)  ! List of unlocated points

    type(c_ptr)         :: c_unlocated
    integer(c_int)      :: n_unlocated

    interface
      function PDM_mesh_location_unlocated_get_cf(mloc, &
                                                  i_point_cloud, &
                                                  i_part) &
      result(unlocated) &
        bind (c, name = 'PDM_mesh_location_unlocated_get')
        use iso_c_binding
        implicit none
        type (c_ptr),   value :: mloc
        integer(c_int), value :: i_point_cloud
        integer(c_int), value :: i_part
        type(c_ptr)           :: unlocated
      end function PDM_mesh_location_unlocated_get_cf
    end interface

    n_unlocated = PDM_mesh_location_n_unlocated_get_cf(mloc,          &
                                                       i_point_cloud, &
                                                       i_part)

    c_unlocated = PDM_mesh_location_unlocated_get_cf(mloc,          &
                                                     i_point_cloud, &
                                                     i_part)

    call c_f_pointer(c_unlocated,   &
                     unlocated,     &
                     [n_unlocated])

  end subroutine PDM_mesh_location_unlocated_get



  subroutine PDM_mesh_location_point_location_get(mloc, &
                                                  i_point_cloud, &
                                                  i_part, &
                                                  location, &
                                                  dist2, &
                                                  projected_coords)
    ! Get point location
    !
    ! .. note::
    !   The results are related to located points only
    implicit none

    type (c_ptr), value                :: mloc                  ! C pointer to PDM_mesh_location_t object
    integer, intent(in)                :: i_point_cloud         ! Point cloud identifier
    integer, intent(in)                :: i_part                ! Partition identifier
    integer(kind=pdm_g_num_s), pointer :: location(:)           ! Global id of nearest mesh element for located points (size = *n_located*)
    real(8),                   pointer :: dist2(:)              ! Signed squared distance from nearest element (negative if the point is located inside that element) (size = *n_located*)
    real(8),                   pointer :: projected_coords(:,:) ! Cartesian coordinates of projection onto the nearest element (identity if the point is located inside that element)  (shape = [3, *n_located*])

    type(c_ptr)                        :: c_location
    type(c_ptr)                        :: c_dist2
    type(c_ptr)                        :: c_projected_coords
    integer(c_int)                     :: n_located

    interface
      subroutine PDM_mesh_location_point_location_get_cf(mloc, &
                                                         i_point_cloud, &
                                                         i_part, &
                                                         location, &
                                                         dist2, &
                                                         projected_coords) &
      bind (c, name = 'PDM_mesh_location_point_location_get')
        use iso_c_binding
        implicit none
        type (c_ptr), value :: mloc
        integer(c_int), value :: i_point_cloud
        integer(c_int), value :: i_part
        type(c_ptr)           :: location
        type(c_ptr)           :: dist2
        type(c_ptr)           :: projected_coords
      end subroutine PDM_mesh_location_point_location_get_cf
    end interface

    n_located = PDM_mesh_location_n_located_get_cf(mloc,          &
                                                   i_point_cloud, &
                                                   i_part)

    c_location         = C_NULL_PTR
    c_dist2            = C_NULL_PTR
    c_projected_coords = C_NULL_PTR

    call PDM_mesh_location_point_location_get_cf(mloc, &
                                                 i_point_cloud, &
                                                 i_part, &
                                                 c_location, &
                                                 c_dist2, &
                                                 c_projected_coords)

    call c_f_pointer(c_location, &
                     location,   &
                     [n_located])

    call c_f_pointer(c_dist2, &
                     dist2,   &
                     [n_located])

    call c_f_pointer(c_projected_coords, &
                     projected_coords,   &
                     [3,n_located])

  end subroutine PDM_mesh_location_point_location_get



  subroutine PDM_mesh_location_points_in_elt_get(mloc, &
                                                 i_point_cloud, &
                                                 i_part, &
                                                 elt_pts_inside_idx, &
                                                 points_gnum, &
                                                 points_coords, &
                                                 points_uvw, &
                                                 points_weights_idx, &
                                                 points_weights, &
                                                 points_dist2, &
                                                 points_projected_coords)
    ! Get location data for points located in elements
    implicit none

    type (c_ptr), value                :: mloc                         ! C pointer to PDM_mesh_location_t object
    integer, intent(in)                :: i_point_cloud                ! Point cloud identifier
    integer, intent(in)                :: i_part                       ! Partition identifier
    integer(kind=pdm_l_num_s), pointer :: elt_pts_inside_idx(:)        ! Index for element -> points mapping (size = *n_elt* + 1)
    integer(kind=pdm_g_num_s), pointer :: points_gnum(:)               ! Located points global ids (size = ``elt_pts_inside_idx(n_elt+1)``)
    real(8),                   pointer :: points_coords(:,:)           ! Located points cartesian coordinates (shape = [3, ``elt_pts_inside_idx(n_elt+1)``])
    real(8),                   pointer :: points_uvw(:,:)              ! Located points parametric coordinates (shape = [3, ``elt_pts_inside_idx(n_elt+1)``])
    integer(kind=pdm_l_num_s), pointer :: points_weights_idx(:)        ! Index for interpolation weights (size = ``elt_pts_inside_idx(n_elt+1)`` + 1)
    real(8),                   pointer :: points_weights(:)            ! Interpolation weights (size = ``points_weights_idx(elt_pts_inside_idx(n_elt+1)+1)``)
    real(8),                   pointer :: points_dist2(:)              ! Signed squared distance element-points (< 0 if the point is inside) (size = ``elt_pts_inside_idx(n_elt+1)``)
    real(8),                   pointer :: points_projected_coords(:,:) ! Cartesian coordinates of projection on element (identity if the point is inside) (shape = [3, ``elt_pts_inside_idx(n_elt+1)``])

    type(c_ptr)                        :: c_elt_pts_inside_idx
    type(c_ptr)                        :: c_points_gnum
    type(c_ptr)                        :: c_points_coords
    type(c_ptr)                        :: c_points_uvw
    type(c_ptr)                        :: c_points_weights_idx
    type(c_ptr)                        :: c_points_weights
    type(c_ptr)                        :: c_points_dist2
    type(c_ptr)                        :: c_points_projected_coords
    integer(c_int)                     :: n_elt
    integer                            :: n_pts_t

    interface
      subroutine PDM_mesh_location_points_in_elt_get_cf(mloc, &
                                                        i_point_cloud, &
                                                        i_part, &
                                                        elt_pts_inside_idx, &
                                                        points_gnum, &
                                                        points_coords, &
                                                        points_uvw, &
                                                        points_weights_idx, &
                                                        points_weights, &
                                                        points_dist2, &
                                                        points_projected_coords) &

        bind (c, name = 'PDM_mesh_location_points_in_elt_get')
        use iso_c_binding
        implicit none
        type (c_ptr),   value :: mloc
        integer(c_int), value :: i_point_cloud
        integer(c_int), value :: i_part
        type(c_ptr)           :: elt_pts_inside_idx
        type(c_ptr)           :: points_gnum
        type(c_ptr)           :: points_coords
        type(c_ptr)           :: points_uvw
        type(c_ptr)           :: points_weights_idx
        type(c_ptr)           :: points_weights
        type(c_ptr)           :: points_dist2
        type(c_ptr)           :: points_projected_coords
      end subroutine PDM_mesh_location_points_in_elt_get_cf
    end interface

    n_elt = pdm_mesh_location_n_cell_get(mloc,     &
                                         i_part)

    c_elt_pts_inside_idx      = C_NULL_PTR
    c_points_gnum             = C_NULL_PTR
    c_points_coords           = C_NULL_PTR
    c_points_uvw              = C_NULL_PTR
    c_points_weights_idx      = C_NULL_PTR
    c_points_weights          = C_NULL_PTR
    c_points_dist2            = C_NULL_PTR
    c_points_projected_coords = C_NULL_PTR
    
    call PDM_mesh_location_points_in_elt_get_cf(mloc, &
                                                i_point_cloud, &
                                                i_part, &
                                                c_elt_pts_inside_idx, &
                                                c_points_gnum, &
                                                c_points_coords, &
                                                c_points_uvw, &
                                                c_points_weights_idx, &
                                                c_points_weights, &
                                                c_points_dist2, &
                                                c_points_projected_coords)

    call c_f_pointer(c_elt_pts_inside_idx, &
                     elt_pts_inside_idx,   &
                     [n_elt + 1])

    n_pts_t = elt_pts_inside_idx(n_elt+1)

    call c_f_pointer(c_points_gnum, &
                     points_gnum,   &
                     [n_pts_t])

    call c_f_pointer(c_points_coords, &
                     points_coords,   &
                     [3,n_pts_t])

    call c_f_pointer(c_points_uvw, &
                     points_uvw,   &
                     [3,n_pts_t])

    call c_f_pointer(c_points_weights_idx, &
                     points_weights_idx,   &
                     [n_pts_t+1])

    call c_f_pointer(c_points_weights, &
                     points_weights,   &
                     [points_weights_idx(n_pts_t+1)])

    call c_f_pointer(c_points_gnum, &
                     points_gnum,   &
                     [elt_pts_inside_idx(n_elt+1)])

    call c_f_pointer(c_points_dist2, &
                     points_dist2,   &
                     [n_pts_t])

    call c_f_pointer(c_points_projected_coords, &
                     points_projected_coords,   &
                     [3,n_pts_t])

  end subroutine PDM_mesh_location_points_in_elt_get


  subroutine PDM_mesh_location_cell_vertex_get_cptr(mloc,         &
                                                    i_part,       &
                                                    cell_vtx_idx, &
                                                    cell_vtx)

    use iso_c_binding

    implicit none


    type (c_ptr),        value :: mloc
    integer(c_int), intent(in) :: i_part
    type(c_ptr)                :: cell_vtx_idx
    type(c_ptr)                :: cell_vtx

    call PDM_mesh_location_cell_vertex_get_cf(mloc,         &
                                              i_part,       &
                                              cell_vtx_idx, &
                                              cell_vtx)

  end subroutine PDM_mesh_location_cell_vertex_get_cptr


  subroutine PDM_mesh_location_cell_vertex_get_f(mloc,         &
                                                 i_part,       &
                                                 cell_vtx_idx, &
                                                 cell_vtx)
    ! Get the cell→vertex connectivity used for internal computations
    !
    ! .. note::
    !   For non-standard elements, this connectivity is built by ParaDiGM and is necessary to associate
    !   the `points_weights` array (returned by \ref PDM_mesh_location_points_in_elt_get)
    !   to the appropriate mesh vertices.
    use iso_c_binding

    implicit none


    type (c_ptr),           value :: mloc            ! C pointer to PDM_mesh_location_t object
    integer(c_int),    intent(in) :: i_part          ! Partition identifier
    integer(pdm_l_num_s), pointer :: cell_vtx_idx(:) ! Index for cell -> vertex connectivity
    integer(pdm_l_num_s), pointer :: cell_vtx(:)     ! Cell -> vertex connectivity
    type(c_ptr)                   :: c_cell_vtx_idx
    type(c_ptr)                   :: c_cell_vtx
    integer                       :: n_cell

    c_cell_vtx_idx = C_NULL_PTR
    c_cell_vtx     = C_NULL_PTR

    call PDM_mesh_location_cell_vertex_get_cf(mloc,           &
                                              i_part,         &
                                              c_cell_vtx_idx, &
                                              c_cell_vtx)

    n_cell = pdm_mesh_location_n_cell_get(mloc, i_part)

    call c_f_pointer(c_cell_vtx_idx, cell_vtx_idx, [n_cell+1])
    call c_f_pointer(c_cell_vtx,     cell_vtx,     [cell_vtx_idx(n_cell+1)])

  end subroutine PDM_mesh_location_cell_vertex_get_f



  subroutine PDM_mesh_location_part_to_part_get(mloc,   &
                                                icloud, &
                                                ptp,    &
                                                owner)
    ! Get part_to_part object to exchange data between the source mesh and a target point cloud
    implicit none

    type (c_ptr),   value :: mloc   ! C pointer to PDM_mesh_location_t object
    integer(c_int), value :: icloud ! Point cloud identifier
    type (c_ptr)          :: ptp    ! Pointer to PDM_part_to_part object
    integer(c_int), value :: owner  ! Ownership for ``ptp``

    interface
      subroutine PDM_mesh_location_part_to_part_get_cf(mloc,   &
                                                       icloud, &
                                                       ptp,    &
                                                       owner)  &
      bind (c, name = 'PDM_mesh_location_part_to_part_get')
        use iso_c_binding
        implicit none
        type (c_ptr),   value :: mloc
        integer(c_int), value :: icloud
        type (c_ptr)          :: ptp
        integer(c_int), value :: owner
      end subroutine PDM_mesh_location_part_to_part_get_cf
    end interface

    call PDM_mesh_location_part_to_part_get_cf(mloc,   &
                                               icloud, &
                                               ptp,    &
                                               owner)

  end subroutine PDM_mesh_location_part_to_part_get



  subroutine PDM_mesh_location_free(mloc)
    ! Free a Mesh Location structure
    implicit none

    type (c_ptr), intent(inout) :: mloc

    interface
      subroutine PDM_mesh_location_free_cf(mloc) &
        bind (c, name = 'PDM_mesh_location_free')
        use iso_c_binding
        implicit none
        type (c_ptr), value :: mloc
      end subroutine PDM_mesh_location_free_cf
    end interface

    call PDM_mesh_location_free_cf(mloc)

  end subroutine PDM_mesh_location_free

end module pdm_mesh_location
