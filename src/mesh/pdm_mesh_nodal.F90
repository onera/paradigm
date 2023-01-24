!
! File:   pdm_mesh_nodal.F90
! Author: equemera
!
! Created on July 10, 2017, 1:34 PM
!

module pdm_mesh_nodal

  use pdm
  use iso_c_binding

  implicit none

  contains


!> Create a Mesh nodal structure
!!
!! @param[in]   n_part   Number of partition on the current process
!! @param[in]   f_comm   MPI communicator
!! @param[out]  mesh     Pointer to \ref PDM_Mesh_nodal object
!!

  subroutine PDM_mesh_nodal_create (mesh, n_part, f_comm)
    use iso_c_binding

    implicit none

    type(c_ptr)          :: mesh
    integer, intent(in)  :: n_part
    integer, intent(in)  :: f_comm

    integer(c_int)       :: c_comm

    interface
      function PDM_mesh_nodal_create_c (n_part, c_comm) result(mesh) bind(c, name='PDM_Mesh_nodal_create')
        use iso_c_binding
        implicit none
        integer(c_int), value :: n_part
        integer(c_int), value :: c_comm
        type(c_ptr)           :: mesh
      end function PDM_mesh_nodal_create_c
    end interface

    c_comm = PDM_MPI_Comm_f2c(f_comm)

    mesh = PDM_mesh_nodal_create_c (n_part, c_comm)

  end subroutine PDM_mesh_nodal_create


!> \brief Free partially a nodal mesh structure
!!
!! @param[in]  mesh       Pointer to \ref PDM_Mesh_nodal object
!!

  subroutine PDM_mesh_nodal_partial_free (mesh)
    use iso_c_binding

    implicit none

    type(c_ptr), value     :: mesh

    interface
      subroutine PDM_mesh_nodal_partial_free_c (mesh) &
        bind(c, name='PDM_Mesh_nodal_partial_free')

        use iso_c_binding

        implicit none

        type(c_ptr), value     :: mesh
      end subroutine PDM_mesh_nodal_partial_free_c
    end interface

    call PDM_mesh_nodal_partial_free_c (mesh)

  end subroutine PDM_mesh_nodal_partial_free


!> Free a nodal mesh structure
!!
!! @param[in]  mesh       Pointer to \ref PDM_Mesh_nodal object
!!

  subroutine PDM_mesh_nodal_free (mesh)
    use iso_c_binding

    implicit none

    type(c_ptr), value :: mesh

    interface
      subroutine PDM_mesh_nodal_free_c (mesh) &
        bind(c, name='PDM_Mesh_nodal_free')

        use iso_c_binding

        implicit none

        type(c_ptr), value :: mesh
      end subroutine PDM_mesh_nodal_free_c
    end interface

    call PDM_mesh_nodal_free_c (mesh)

  end subroutine PDM_mesh_nodal_free


!> Define partition vertices
!!
!! @param[in]  mesh     Pointer to \ref PDM_Mesh_nodal object
!! @param[in]  id_part  Partition identifier
!! @param[in]  n_vtx    Number of vertices
!! @param[in]  coords   Interlaced coordinates (size = 3 * \ref n_vtx)
!! @param[in]  numabs   Global numbering
!!

  subroutine PDM_mesh_nodal_coord_set (mesh, id_part, n_vtx, coords, numabs, owner)
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
      subroutine PDM_mesh_nodal_coord_set_c (mesh, id_part, n_vtx, coords, numabs, owner) &
        bind(c, name='PDM_Mesh_nodal_coord_set')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer (c_int), intent(in), value :: n_vtx
        integer (c_int), intent(in), value :: owner
        type (c_ptr),                value :: coords
        type (c_ptr),                value :: numabs
      end subroutine PDM_mesh_nodal_coord_set_c
    end interface

    c_coords = c_loc (coords)
    c_numabs = c_loc (numabs)

    call  PDM_mesh_nodal_coord_set_c (mesh, id_part, n_vtx, c_coords, c_numabs, owner)

  end subroutine

!> Define standard 3D cells by cell-vertex connectivity
!!
!! @param[in]  mesh          Pointer to \ref PDM_Mesh_nodal object
!! @param[in]  id_part       Partition identifier
!! @param[in]  n_cell        Number of cells
!! @param[in]  cell_vtx_idx  Index of cell vertex connectivity
!! @param[in]  cell_vtx_nb   Number of vertices for each cell
!! @param[in]  cell_vtx      Cell vertex connectivity
!! @param[in]  numabs        Global numbering
!! @param[in]  owner         Ownership

  subroutine PDM_Mesh_nodal_cells_cellvtx_add (mesh, id_part, n_cell, cell_vtx_idx, cell_vtx_nb, cell_vtx, numabs, owner)
    use iso_c_binding

    implicit none

    type(c_ptr), value                  :: mesh
    integer, intent(in)                 :: id_part
    integer, intent(in)                 :: n_cell
    integer, intent(in)                 :: owner
    integer (pdm_l_num_s), pointer      :: cell_vtx_idx(:)
    integer (pdm_l_num_s), pointer      :: cell_vtx_nb(:)
    integer (pdm_l_num_s), pointer      :: cell_vtx(:)
    integer (pdm_g_num_s), pointer      :: numabs(:)
    type(c_ptr) :: c_cell_vtx_idx = C_NULL_PTR
    type(c_ptr) :: c_cell_vtx_nb = C_NULL_PTR
    type(c_ptr) :: c_cell_vtx = C_NULL_PTR
    type(c_ptr) :: c_numabs = C_NULL_PTR

    interface
      subroutine PDM_Mesh_nodal_cells_cellvtx_add_c (mesh, id_part, n_cell, cell_vtx_idx, cell_vtx_nb, cell_vtx, numabs, owner) &
        bind(c, name='PDM_Mesh_nodal_cells_cellvtx_add')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer (c_int), intent(in), value :: n_cell
        integer (c_int), intent(in), value :: owner
        type (c_ptr),                value :: cell_vtx_idx
        type (c_ptr),                value :: cell_vtx_nb
        type (c_ptr),                value :: cell_vtx
        type (c_ptr),                value :: numabs
      end subroutine PDM_Mesh_nodal_cells_cellvtx_add_c
    end interface

    c_cell_vtx_idx = c_loc (cell_vtx_idx)
    c_cell_vtx_nb = c_loc (cell_vtx_nb)
    c_cell_vtx = c_loc (cell_vtx)
    c_numabs = c_loc (numabs)

    call  PDM_Mesh_nodal_cells_cellvtx_add_c (mesh, id_part, n_cell, c_cell_vtx_idx, c_cell_vtx_nb, c_cell_vtx, c_numabs, owner)

  end subroutine

!> Define faces by face-vertex connectivity
!!
!! @param[in]  mesh          Pointer to \ref PDM_Mesh_nodal object
!! @param[in]  id_part       Partition identifier
!! @param[in]  n_face        Number of faces
!! @param[in]  face_vtx_idx  Index of face vertex connectivity
!! @param[in]  face_vtx_nb   Number of vertices for each face
!! @param[in]  face_vtx      Face vertex connectivity
!! @param[in]  numabs        Global numbering
!! @param[in]  owner         Ownership

  subroutine PDM_Mesh_nodal_faces_facevtx_add (mesh, id_part, n_face, face_vtx_idx, face_vtx_nb, face_vtx, numabs, owner)
    use iso_c_binding

    implicit none

    type(c_ptr), value                  :: mesh
    integer, intent(in)                 :: id_part
    integer, intent(in)                 :: n_face
    integer, intent(in)                 :: owner
    integer (pdm_l_num_s), pointer      :: face_vtx_idx(:)
    integer (pdm_l_num_s), pointer      :: face_vtx_nb(:)
    integer (pdm_l_num_s), pointer      :: face_vtx(:)
    integer (pdm_g_num_s), pointer      :: numabs(:)
    type(c_ptr) :: c_face_vtx_idx = C_NULL_PTR
    type(c_ptr) :: c_face_vtx_nb = C_NULL_PTR
    type(c_ptr) :: c_face_vtx = C_NULL_PTR
    type(c_ptr) :: c_numabs = C_NULL_PTR

    interface
      subroutine PDM_Mesh_nodal_faces_facevtx_add_c (mesh, id_part, n_face, face_vtx_idx, face_vtx_nb, face_vtx, numabs, owner) &
        bind(c, name='PDM_Mesh_nodal_faces_facevtx_add')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer (c_int), intent(in), value :: n_face
        integer (c_int), intent(in), value :: owner
        type (c_ptr),                value :: face_vtx_idx
        type (c_ptr),                value :: face_vtx_nb
        type (c_ptr),                value :: face_vtx
        type (c_ptr),                value :: numabs
      end subroutine PDM_Mesh_nodal_faces_facevtx_add_c
    end interface

    c_face_vtx_idx = c_loc (face_vtx_idx)
    c_face_vtx_nb = c_loc (face_vtx_nb)
    c_face_vtx = c_loc (face_vtx)
    c_numabs = c_loc (numabs)

    call  PDM_Mesh_nodal_faces_facevtx_add_c (mesh, id_part, n_face, c_face_vtx_idx, c_face_vtx_nb, c_face_vtx, c_numabs, owner)

  end subroutine

end module PDM_mesh_nodal
