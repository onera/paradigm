#include "pdm_configf.h"

module pdm_part_connectivity_transform

  use pdm
  use iso_c_binding
  implicit none


  contains

    subroutine PDM_combine_connectivity(n_entity1,           &
                                        entity1_entity2_idx, &
                                        entity1_entity2,     &
                                        entity2_entity3_idx, &
                                        entity2_entity3,     &
                                        entity1_entity3_idx, &
                                        entity1_entity3)

      ! Combine connectivity between entity1_entity2 and entity2_entity3 to have entity1_entity3
      implicit none

      integer,              intent(in) :: n_entity1              ! Number of entity1
      integer(pdm_l_num_s), pointer    :: entity1_entity2_idx(:) ! Connectivity index between entity1 and entity2 (size = n_entity1+1)
      integer(pdm_l_num_s), pointer    :: entity1_entity2(:)     ! Connectivity between entity1 and entity2 (size = entity1_entity2_idx(n_entity1+1) )
      integer(pdm_l_num_s), pointer    :: entity2_entity3_idx(:) ! Connectivity index between entity2 and entity3 (size = n_entity2+1)
      integer(pdm_l_num_s), pointer    :: entity2_entity3(:)     ! Connectivity between entity2 and entity3 (size = entity1_entity2_idx(n_entity2+1) )
      integer(pdm_l_num_s), pointer    :: entity1_entity3_idx(:) ! Connectivity index between entity1 and entity3 (size = n_entity1+1)
      integer(pdm_l_num_s), pointer    :: entity1_entity3(:)     ! Connectivity between entity1 and entity3 (size = entity1_entity2_idx(n_entity1+1) )

      integer(c_int) :: c_n_entity1

      type(c_ptr)    :: c_entity1_entity2_idx
      type(c_ptr)    :: c_entity1_entity2
      type(c_ptr)    :: c_entity2_entity3_idx
      type(c_ptr)    :: c_entity2_entity3
      type(c_ptr)    :: c_entity1_entity3_idx
      type(c_ptr)    :: c_entity1_entity3

      interface
        subroutine PDM_combine_connectivity_cf(n_entity1,           &
                                               entity1_entity2_idx, &
                                               entity1_entity2,     &
                                               entity2_entity3_idx, &
                                               entity2_entity3,     &
                                               entity1_entity3_idx, &
                                               entity1_entity3)     &

        bind (c, name = 'PDM_combine_connectivity')
          use iso_c_binding
          implicit none
          integer(c_int), value :: n_entity1
          type(c_ptr),    value :: entity1_entity2_idx
          type(c_ptr),    value :: entity1_entity2
          type(c_ptr),    value :: entity2_entity3_idx
          type(c_ptr),    value :: entity2_entity3
          type(c_ptr)           :: entity1_entity3_idx
          type(c_ptr)           :: entity1_entity3
        end subroutine PDM_combine_connectivity_cf
      end interface

      c_n_entity1 = n_entity1

      c_entity1_entity2_idx = C_NULL_PTR
      if (associated(entity1_entity2_idx)) then
        c_entity1_entity2_idx = c_loc(entity1_entity2_idx)
      endif

      c_entity1_entity2 = C_NULL_PTR
      if (associated(entity1_entity2)) then
        c_entity1_entity2 = c_loc(entity1_entity2)
      endif

      c_entity2_entity3_idx = C_NULL_PTR
      if (associated(entity2_entity3_idx)) then
        c_entity2_entity3_idx = c_loc(entity2_entity3_idx)
      endif

      c_entity2_entity3 = C_NULL_PTR
      if (associated(entity2_entity3)) then
        c_entity2_entity3 = c_loc(entity2_entity3)
      endif

      c_entity1_entity3_idx = C_NULL_PTR
      c_entity1_entity3     = C_NULL_PTR

      call PDM_combine_connectivity_cf(c_n_entity1,           &
                                       c_entity1_entity2_idx, &
                                       c_entity1_entity2,     &
                                       c_entity2_entity3_idx, &
                                       c_entity2_entity3,     &
                                       c_entity1_entity3_idx, &
                                       c_entity1_entity3)

      call c_f_pointer(c_entity1_entity3_idx, &
                       entity1_entity3_idx,   &
                       [n_entity1+1])

      call c_f_pointer(c_entity1_entity3, &
                       entity1_entity3,   &
                       [entity1_entity3_idx(n_entity1+1)])

    end subroutine PDM_combine_connectivity



    subroutine PDM_connectivity_transpose(n_entity1,           &
                                          n_entity2,           &
                                          entity1_entity2_idx, &
                                          entity1_entity2,     &
                                          entity2_entity1_idx, &
                                          entity2_entity1)

      ! Transpose connectivity entity1_entity2 to have entity2_entity1
      implicit none

      integer,              intent(in) :: n_entity1              ! Number of entity1
      integer,              intent(in) :: n_entity2              ! Number of entity2
      integer(pdm_l_num_s), pointer    :: entity1_entity2_idx(:) ! Connectivity index between entity1 and entity2 (size = n_entity1+1)
      integer(pdm_l_num_s), pointer    :: entity1_entity2(:)     ! Connectivity between entity1 and entity2 (size = entity1_entity2_idx(n_entity1+1) )
      integer(pdm_l_num_s), pointer    :: entity2_entity1_idx(:) ! Connectivity index between entity2 and entity1 (size = n_entity2+1)
      integer(pdm_l_num_s), pointer    :: entity2_entity1(:)     ! Connectivity between entity2 and entity1 (size = entity1_entity2_idx(n_entity2+1) )

      integer(c_int) :: c_n_entity1
      integer(c_int) :: c_n_entity2

      type(c_ptr)    :: c_entity1_entity2_idx
      type(c_ptr)    :: c_entity1_entity2
      type(c_ptr)    :: c_entity2_entity1_idx
      type(c_ptr)    :: c_entity2_entity1

      interface
        subroutine PDM_connectivity_transpose_cf(n_entity1,           &
                                                 n_entity2,           &
                                                 entity1_entity2_idx, &
                                                 entity1_entity2,     &
                                                 entity2_entity1_idx, &
                                                 entity2_entity1)     &
        bind (c, name = 'PDM_connectivity_transpose')
          use iso_c_binding
          implicit none
          integer(c_int), value :: n_entity1
          integer(c_int), value :: n_entity2
          type(c_ptr),    value :: entity1_entity2_idx
          type(c_ptr),    value :: entity1_entity2
          type(c_ptr)           :: entity2_entity1_idx
          type(c_ptr)           :: entity2_entity1
        end subroutine PDM_connectivity_transpose_cf
      end interface

      c_n_entity1 = n_entity1
      c_n_entity2 = n_entity2

      c_entity1_entity2_idx = c_loc(entity1_entity2_idx)
      c_entity1_entity2     = c_loc(entity1_entity2)
      c_entity2_entity1_idx = C_NULL_PTR
      c_entity2_entity1     = C_NULL_PTR

      call PDM_connectivity_transpose_cf(c_n_entity1,           &
                                         c_n_entity2,           &
                                         c_entity1_entity2_idx, &
                                         c_entity1_entity2,     &
                                         c_entity2_entity1_idx, &
                                         c_entity2_entity1)

      call c_f_pointer(c_entity2_entity1_idx, &
                       entity2_entity1_idx,   &
                       [n_entity2+1])

      call c_f_pointer(c_entity2_entity1, &
                       entity2_entity1,   &
                       [entity2_entity1_idx(n_entity2+1)])

    end subroutine PDM_connectivity_transpose




    subroutine PDM_compute_face_vtx_from_face_and_edge(n_face,        &
                                                       face_edge_idx, &
                                                       face_edge,     &
                                                       edge_vtx,      &
                                                       face_vtx)
      ! Combine face->edge with edge-vtx to create face-vtx preserving orientation
      implicit none

      integer,              intent(in) :: n_face           ! Number of faces
      integer(pdm_l_num_s), pointer    :: face_edge_idx(:) ! Index of face->edge connectivity
      integer(pdm_l_num_s), pointer    :: face_edge(:)     ! Face->edge connectivity
      integer(pdm_l_num_s), pointer    :: edge_vtx(:)      ! Edge->vtx connectivity
      integer(pdm_l_num_s), pointer    :: face_vtx(:)      ! Face->vtx connectivity

      integer(c_int) :: c_n_face

      type(c_ptr)    :: c_face_edge_idx     = C_NULL_PTR
      type(c_ptr)    :: c_face_edge         = C_NULL_PTR
      type(c_ptr)    :: c_edge_vtx          = C_NULL_PTR
      type(c_ptr)    :: c_face_vtx          = C_NULL_PTR

      interface
        subroutine PDM_compute_face_vtx_from_face_and_edge_cf(n_face,        &
                                                              face_edge_idx, &
                                                              face_edge,     &
                                                              edge_vtx,      &
                                                              face_vtx)      &
        bind (c, name = 'PDM_compute_face_vtx_from_face_and_edge')
          use iso_c_binding
          implicit none
          integer(c_int), value :: n_face
          type(c_ptr),    value :: face_edge_idx
          type(c_ptr),    value :: face_edge
          type(c_ptr),    value :: edge_vtx
          type(c_ptr)           :: face_vtx
        end subroutine PDM_compute_face_vtx_from_face_and_edge_cf
      end interface

      c_n_face = n_face

      c_face_edge_idx = c_loc(face_edge_idx)
      c_face_edge     = c_loc(face_edge)
      c_edge_vtx      = c_loc(edge_vtx)

      call PDM_compute_face_vtx_from_face_and_edge_cf(c_n_face,        &
                                                      c_face_edge_idx, &
                                                      c_face_edge,     &
                                                      c_edge_vtx,      &
                                                      c_face_vtx)

      call c_f_pointer(c_face_vtx, &
                       face_vtx,   &
                       [face_edge_idx(n_face+1)])

    end subroutine PDM_compute_face_vtx_from_face_and_edge

end module pdm_part_connectivity_transform
