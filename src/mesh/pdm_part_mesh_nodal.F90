!
! File:   pdm_part_mesh_nodal.F90
! Author: ndelling
!
! Created on March 8, 2023, 9:34 PM
!

module PDM_part_mesh_nodal

  use pdm
  use iso_c_binding

  implicit none

  contains

!> Create a PDM_part_mesh_nodal structure
!!
!! @param[in]   mesh_dimension   Mesh dimension
!! @param[in]   n_part           Number of partition on the current process
!! @param[in]   f_comm           MPI communicator
!! @param[out]  mesh             Pointer to \ref PDM_part_mesh_nodal object
!!

  subroutine PDM_part_mesh_nodal_create (mesh, mesh_dimension, n_part, f_comm)
    use iso_c_binding

    implicit none

    type(c_ptr)          :: mesh
    integer, intent(in)  :: mesh_dimension
    integer, intent(in)  :: n_part
    integer, intent(in)  :: f_comm

    integer(c_int)       :: c_comm

    interface
      function PDM_part_mesh_nodal_create_c (mesh_dimension, n_part, c_comm) result(mesh) bind(c, name='PDM_part_mesh_nodal_create')
        use iso_c_binding
        implicit none
        integer(c_int), value :: mesh_dimension
        integer(c_int), value :: n_part
        integer(c_int), value :: c_comm
        type(c_ptr)           :: mesh
      end function PDM_part_mesh_nodal_create_c
    end interface

    c_comm = PDM_MPI_Comm_f2c(f_comm)

    mesh = PDM_part_mesh_nodal_create_c (mesh_dimension, n_part, c_comm)

  end subroutine PDM_part_mesh_nodal_create

!> \brief Free partially a PDM_part_mesh_nodal structure
!!
!! @param[in]  mesh       Pointer to \ref PDM_part_mesh_nodal object
!!

  subroutine PDM_part_mesh_nodal_partial_free (mesh)
    use iso_c_binding

    implicit none

    type(c_ptr), value     :: mesh

    interface
      subroutine PDM_part_mesh_nodal_partial_free_c (mesh) &
        bind(c, name='PDM_part_mesh_nodal_partial_free')

        use iso_c_binding

        implicit none

        type(c_ptr), value     :: mesh
      end subroutine PDM_part_mesh_nodal_partial_free_c
    end interface

    call PDM_part_mesh_nodal_partial_free_c (mesh)

  end subroutine PDM_part_mesh_nodal_partial_free

!> Free a PDM_part_mesh_nodal structure
!!
!! @param[in]  mesh       Pointer to \ref PDM_part_mesh_nodal object
!!

  subroutine PDM_part_mesh_nodal_free (mesh)
    use iso_c_binding

    implicit none

    type(c_ptr), value :: mesh

    interface
      subroutine PDM_part_mesh_nodal_free_c (mesh) &
        bind(c, name='PDM_part_mesh_nodal_free')

        use iso_c_binding

        implicit none

        type(c_ptr), value :: mesh
      end subroutine PDM_part_mesh_nodal_free_c
    end interface

    call PDM_part_mesh_nodal_free_c (mesh)

  end subroutine PDM_part_mesh_nodal_free

!> Define partition vertices
!!
!! @param[in]  mesh     Pointer to \ref PDM_part_mesh_nodal object
!! @param[in]  id_part  Partition identifier
!! @param[in]  n_vtx    Number of vertices
!! @param[in]  coords   Interlaced coordinates (size = 3 * \ref n_vtx)
!! @param[in]  numabs   Global numbering
!!

  subroutine PDM_part_mesh_nodal_coord_set (mesh, id_part, n_vtx, coords, numabs, owner)
    use iso_c_binding

    implicit none

    type(c_ptr), value                  :: mesh
    integer, intent(in)                 :: id_part
    integer, intent(in)                 :: n_vtx
    integer, intent(in)                 :: owner
    double precision, pointer           :: coords(:,:)
    integer (pdm_g_num_s), pointer      :: numabs(:)
    type(c_ptr) :: c_coords = C_NULL_PTR
    type(c_ptr) :: c_numabs = C_NULL_PTR

    interface
      subroutine PDM_part_mesh_nodal_coord_set_c (mesh, id_part, n_vtx, coords, numabs, owner) &
        bind(c, name='PDM_part_mesh_nodal_coord_set')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer (c_int), intent(in), value :: n_vtx
        integer (c_int), intent(in), value :: owner
        type (c_ptr),                value :: coords
        type (c_ptr),                value :: numabs
      end subroutine PDM_part_mesh_nodal_coord_set_c
    end interface

    c_coords = c_loc (coords)
    c_numabs = c_loc (numabs)

    call  PDM_part_mesh_nodal_coord_set_c (mesh, id_part, n_vtx, c_coords, c_numabs, owner)

  end subroutine PDM_part_mesh_nodal_coord_set

!> Define standard 3D cells by cell-vertex connectivity
!!
!! @param[in]  mesh          Pointer to \ref PDM_part_mesh_nodal object
!! @param[in]  id_part       Partition identifier
!! @param[in]  n_cell        Number of cells
!! @param[in]  cell_vtx_idx  Index of cell vertex connectivity
!! @param[in]  cell_vtx      Cell vertex connectivity
!! @param[in]  numabs        Global numbering
!! @param[in]  owner         Ownership

  subroutine PDM_part_mesh_nodal_cells_cellvtx_add (mesh, id_part, n_cell, cell_vtx_idx, cell_vtx, numabs, owner)
    use iso_c_binding

    implicit none

    type(c_ptr), value                  :: mesh
    integer, intent(in)                 :: id_part
    integer, intent(in)                 :: n_cell
    integer, intent(in)                 :: owner
    integer (pdm_l_num_s), pointer      :: cell_vtx_idx(:)
    integer (pdm_l_num_s), pointer      :: cell_vtx(:)
    integer (pdm_g_num_s), pointer      :: numabs(:)
    type(c_ptr) :: c_cell_vtx_idx = C_NULL_PTR
    type(c_ptr) :: c_cell_vtx = C_NULL_PTR
    type(c_ptr) :: c_numabs = C_NULL_PTR

    interface
      subroutine PDM_part_mesh_nodal_cells_cellvtx_add_c (mesh, id_part, n_cell, cell_vtx_idx, cell_vtx, &
        numabs, owner) bind(c, name='PDM_part_mesh_nodal_cells_cellvtx_add')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer (c_int), intent(in), value :: n_cell
        integer (c_int), intent(in), value :: owner
        type (c_ptr),                value :: cell_vtx_idx
        type (c_ptr),                value :: cell_vtx
        type (c_ptr),                value :: numabs
      end subroutine PDM_part_mesh_nodal_cells_cellvtx_add_c
    end interface

    c_cell_vtx_idx = c_loc (cell_vtx_idx)
    c_cell_vtx = c_loc (cell_vtx)
    c_numabs = c_loc (numabs)

    call  PDM_part_mesh_nodal_cells_cellvtx_add_c (mesh, id_part, n_cell, c_cell_vtx_idx, c_cell_vtx, &
                                                   c_numabs, owner)

  end subroutine PDM_part_mesh_nodal_cells_cellvtx_add

!> Define faces by face-vertex connectivity
!!
!! @param[in]  mesh          Pointer to \ref PDM_part_mesh_nodal object
!! @param[in]  id_part       Partition identifier
!! @param[in]  n_face        Number of faces
!! @param[in]  face_vtx_idx  Index of face vertex connectivity
!! @param[in]  face_vtx      Face vertex connectivity
!! @param[in]  numabs        Global numbering
!! @param[in]  owner         Ownership

  subroutine PDM_part_mesh_nodal_faces_facevtx_add (mesh, id_part, n_face, face_vtx_idx, face_vtx, numabs, owner)
    use iso_c_binding

    implicit none

    type(c_ptr), value                  :: mesh
    integer, intent(in)                 :: id_part
    integer, intent(in)                 :: n_face
    integer, intent(in)                 :: owner
    integer (pdm_l_num_s), pointer      :: face_vtx_idx(:)
    integer (pdm_l_num_s), pointer      :: face_vtx(:)
    integer (pdm_g_num_s), pointer      :: numabs(:)
    type(c_ptr) :: c_face_vtx_idx = C_NULL_PTR
    type(c_ptr) :: c_face_vtx = C_NULL_PTR
    type(c_ptr) :: c_numabs = C_NULL_PTR

    interface
      subroutine PDM_part_mesh_nodal_faces_facevtx_add_c (mesh, id_part, n_face, face_vtx_idx, face_vtx, &
        numabs, owner) bind(c, name='PDM_part_mesh_nodal_faces_facevtx_add')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer (c_int), intent(in), value :: n_face
        integer (c_int), intent(in), value :: owner
        type (c_ptr),                value :: face_vtx_idx
        type (c_ptr),                value :: face_vtx
        type (c_ptr),                value :: numabs
      end subroutine PDM_part_mesh_nodal_faces_facevtx_add_c
    end interface

    c_face_vtx_idx = c_loc (face_vtx_idx)
    c_face_vtx = c_loc (face_vtx)
    c_numabs = c_loc (numabs)

    call  PDM_part_mesh_nodal_faces_facevtx_add_c (mesh, id_part, n_face, c_face_vtx_idx, c_face_vtx, &
                                                   c_numabs, owner)

  end subroutine PDM_part_mesh_nodal_faces_facevtx_add

!!> Define 2D cells by cell-edge connectivity
!!!
!!! @param[in]  mesh          Pointer to \ref PDM_part_mesh_nodal object
!!! @param[in]  id_part       Partition identifier
!!! @param[in]  n_elt         Number of polyhedra
!!! @param[in]  n_edge        Number of edges used to describe polyhedra
!!! @param[in]  edge_vtx_idx  Index of edge vertex connectivity
!!! @param[in]  edge_vtx      Edge vertex connectivity
!!! @param[in]  cell_edge_idx Index of cell edge connectivity
!!! @param[in]  cell_edge     Cell edge connectivity
!!! @param[in]  numabs        Global numbering
!!! @param[in]  owner         Ownership
!
!  subroutine PDM_part_mesh_nodal_cell2d_celledge_add (mesh, id_part, n_elt, n_edge, edge_vtx_idx, edge_vtx, &
!                                                      cell_edge_idx, cell_edge, numabs, owner)
!    use iso_c_binding
!
!    implicit none
!
!    type(c_ptr), value                  :: mesh
!    integer, intent(in)                 :: id_part
!    integer, intent(in)                 :: n_elt
!    integer, intent(in)                 :: n_edge
!    integer, intent(in)                 :: owner
!    integer (pdm_l_num_s), pointer      :: edge_vtx_idx(:)
!    integer (pdm_l_num_s), pointer      :: edge_vtx(:)
!    integer (pdm_l_num_s), pointer      :: cell_edge_idx(:)
!    integer (pdm_l_num_s), pointer      :: cell_edge(:)
!    integer (pdm_g_num_s), pointer      :: numabs(:)
!    type(c_ptr) :: c_edge_vtx_idx = C_NULL_PTR
!    type(c_ptr) :: c_edge_vtx = C_NULL_PTR
!    type(c_ptr) :: c_cell_edge_idx = C_NULL_PTR
!    type(c_ptr) :: c_cell_edge = C_NULL_PTR
!    type(c_ptr) :: c_numabs = C_NULL_PTR
!
!    interface
!      subroutine PDM_part_mesh_nodal_cell2d_celledge_add_c (mesh, id_part, n_elt, n_edge, edge_vtx_idx, edge_vtx, &
!                                                            cell_edge_idx, cell_edge, numabs, owner) &
!        bind(c, name='PDM_part_mesh_nodal_cell2d_celledge_add')
!
!        use iso_c_binding
!        use pdm
!
!        implicit none
!
!        type(c_ptr),                 value :: mesh
!        integer (c_int), intent(in), value :: id_part
!        integer (c_int), intent(in), value :: n_elt
!        integer (c_int), intent(in), value :: n_edge
!        integer (c_int), intent(in), value :: owner
!        type (c_ptr),                value :: edge_vtx_idx
!        type (c_ptr),                value :: edge_vtx
!        type (c_ptr),                value :: cell_edge_idx
!        type (c_ptr),                value :: cell_edge
!        type (c_ptr),                value :: numabs
!      end subroutine PDM_part_mesh_nodal_cell2d_celledge_add_c
!    end interface
!
!    c_edge_vtx_idx = c_loc (edge_vtx_idx)
!    c_edge_vtx = c_loc (edge_vtx)
!    c_cell_edge_idx = c_loc (cell_edge_idx)
!    c_cell_edge = c_loc (cell_edge)
!    c_numabs = c_loc (numabs)
!
!    call  PDM_part_mesh_nodal_cell2d_celledge_add_c (mesh, id_part, n_elt, n_edge, c_edge_vtx_idx, c_edge_vtx, &
!                                                     c_cell_edge_idx, c_cell_edge, c_numabs, owner)
!
!  end subroutine PDM_part_mesh_nodal_cell2d_celledge_add

!> Define 3D cells by cell-face connectivity
!!
!! @param[in]  mesh          Pointer to \ref PDM_part_mesh_nodal object
!! @param[in]  id_part       Partition identifier
!! @param[in]  n_elt         Number of polyhedra
!! @param[in]  n_face        Number of faces used to describe polyhedra
!! @param[in]  face_vtx_idx  Index of face vertex connectivity
!! @param[in]  face_vtx      Face vertex connectivity
!! @param[in]  face_ln_to_gn Global face numbering
!! @param[in]  cell_face_idx Index of cell face connectivity
!! @param[in]  cell_face     Cell face connectivity
!! @param[in]  numabs        Global numbering
!! @param[in]  owner         Ownership

  subroutine PDM_part_mesh_nodal_cell3d_cellface_add (mesh, id_part, n_elt, n_face, face_vtx_idx, face_vtx, &
                                                      face_ln_to_gn, cell_face_idx, cell_face, numabs, owner)
    use iso_c_binding

    implicit none

    type(c_ptr), value                  :: mesh
    integer, intent(in)                 :: id_part
    integer, intent(in)                 :: n_elt
    integer, intent(in)                 :: n_face
    integer, intent(in)                 :: owner
    integer (pdm_l_num_s), pointer      :: face_vtx_idx(:)
    integer (pdm_l_num_s), pointer      :: face_vtx(:)
    integer (pdm_g_num_s), pointer      :: face_ln_to_gn(:)
    integer (pdm_l_num_s), pointer      :: cell_face_idx(:)
    integer (pdm_l_num_s), pointer      :: cell_face(:)
    integer (pdm_g_num_s), pointer      :: numabs(:)
    type(c_ptr) :: c_face_vtx_idx = C_NULL_PTR
    type(c_ptr) :: c_face_vtx = C_NULL_PTR
    type(c_ptr) :: c_face_ln_to_gn = C_NULL_PTR
    type(c_ptr) :: c_cell_face_idx = C_NULL_PTR
    type(c_ptr) :: c_cell_face = C_NULL_PTR
    type(c_ptr) :: c_numabs = C_NULL_PTR

    interface
      subroutine PDM_part_mesh_nodal_cell3d_cellface_add_c (mesh, id_part, n_elt, n_face, face_vtx_idx, face_vtx, &
                                                            face_ln_to_gn, cell_face_idx, cell_face, numabs, owner) &
        bind(c, name='PDM_part_mesh_nodal_cell3d_cellface_add')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer (c_int), intent(in), value :: n_elt
        integer (c_int), intent(in), value :: n_face
        integer (c_int), intent(in), value :: owner
        type (c_ptr),                value :: face_vtx_idx
        type (c_ptr),                value :: face_vtx
        type (c_ptr),                value :: face_ln_to_gn
        type (c_ptr),                value :: cell_face_idx
        type (c_ptr),                value :: cell_face
        type (c_ptr),                value :: numabs
      end subroutine PDM_part_mesh_nodal_cell3d_cellface_add_c
    end interface

    c_face_vtx_idx = c_loc (face_vtx_idx)
    c_face_vtx = c_loc (face_vtx)
    c_face_ln_to_gn = c_loc (face_ln_to_gn)
    c_cell_face_idx = c_loc (cell_face_idx)
    c_cell_face = c_loc (cell_face)
    c_numabs = c_loc (numabs)

    call  PDM_part_mesh_nodal_cell3d_cellface_add_c (mesh, id_part, n_elt, n_face, c_face_vtx_idx, c_face_vtx, &
                                                     c_face_ln_to_gn, c_cell_face_idx, c_cell_face, c_numabs, owner)

  end subroutine PDM_part_mesh_nodal_cell3d_cellface_add

!!> Get cell-vertex connectivity
!!!
!!! @param[in]  mesh          Pointer to \ref PDM_part_mesh_nodal object
!!! @param[in]  id_part       Partition identifier
!!! @param[out] cell_vtx_idx  Index of cell vertex connectivity
!!! @param[out] cell_face     Cell vertex connectivity
!
!  subroutine PDM_part_mesh_nodal_cell_vtx_connectivity_get (mesh, id_part, cell_vtx_idx, cell_vtx)
!    use iso_c_binding
!
!    implicit none
!
!    type(c_ptr), value             :: mesh
!    integer, intent(in)            :: id_part
!    integer (pdm_l_num_s), pointer :: cell_vtx_idx(:)
!    integer (pdm_l_num_s), pointer :: cell_vtx(:)
!    type(c_ptr)    :: c_cell_vtx_idx = C_NULL_PTR
!    type(c_ptr)    :: c_cell_vtx = C_NULL_PTR
!    integer(c_int) :: n_elt
!
!    interface
!
!      subroutine PDM_part_mesh_nodal_cell_vtx_connectivity_get_c (mesh, id_part, cell_vtx_idx, cell_vtx) &
!        bind(c, name='PDM_part_mesh_nodal_cell_vtx_connectivity_get')
!
!        use iso_c_binding
!        use pdm
!
!        implicit none
!
!        type(c_ptr),                 value :: mesh
!        integer (c_int), intent(in), value :: id_part
!        type (c_ptr)                       :: cell_vtx_idx
!        type (c_ptr)                       :: cell_vtx
!      end subroutine
!
!      function PDM_part_mesh_nodal_n_cell_get_c (mesh, id_part) result (n_elt) &
!        bind(c, name='PDM_part_mesh_nodal_n_cell_get')
!
!        use iso_c_binding
!        use pdm
!
!        implicit none
!
!        type(c_ptr),                 value :: mesh
!        integer (c_int), intent(in), value :: id_part
!        integer(c_int)                     :: n_elt
!      end function
!
!    end interface
!
!    n_elt = PDM_part_mesh_nodal_n_cell_get_c(mesh, id_part)
!
!    call  PDM_part_mesh_nodal_cell_vtx_connectivity_get_c (mesh, id_part, c_cell_vtx_idx, c_cell_vtx)
!
!    call c_f_pointer(c_cell_vtx_idx,   &
!                     cell_vtx_idx,     &
!                     [n_elt+1])
!
!    call c_f_pointer(c_cell_vtx,   &
!                     cell_vtx,     &
!                     [cell_vtx_idx(n_elt+1)])
!
!  end subroutine PDM_part_mesh_nodal_cell_vtx_connectivity_get

end module PDM_part_mesh_nodal
