!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2023  ONERA
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 3 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------

#include "pdm_configf.h"

module pdm_extract_part

  use pdm
  use iso_c_binding

  implicit none

  contains

  ! Not wrapped :
  ! PDM_extract_part_entity_center_set
  ! PDM_extract_part_color_get
  ! PDM_extract_part_init_location_get
  ! PDM_extract_part_part_to_part_group_get
  ! PDM_extract_part_renum_method_set
  ! PDM_extract_part_part_mesh_get

   subroutine PDM_extract_part_create(extrp,              &
                                      dim,                &
                                      n_part_in,          &
                                      n_part_out,         &
                                      extract_kind,       &
                                      split_dual_method,  &
                                      compute_child_gnum, &
                                      ownership,          &
                                      comm)
    ! Build an Extract Part instance
    implicit none

    type(c_ptr)                       :: extrp              ! C pointer to PDM_extract_part_t instance
    integer, intent(in)               :: dim                ! Mesh dimension
    integer, intent(in)               :: n_part_in          ! Number of input  partitions
    integer, intent(in)               :: n_part_out         ! Number of output partitions
    integer, intent(in)               :: extract_kind       ! Extraction kind (local/reequilibrate/from target)
    integer, intent(in)               :: split_dual_method  ! Repartitioning method (used only in PDM_EXTRACT_PART_KIND_REEQUILIBRATE mode)
    logical, intent(in)               :: compute_child_gnum ! Enable generation of new global IDs for extraction
    integer, intent(in)               :: ownership          ! Ownership of extraction
    integer, intent(in)               :: comm               ! MPI communicator

    integer(c_int)                    :: c_compute_child_gnum
    integer(c_int)                    :: c_comm

    interface
      function pdm_extract_part_create_c (dim,                &
                                          n_part_in,          &
                                          n_part_out,         &
                                          extract_kind,       &
                                          split_dual_method,  &
                                          compute_child_gnum, &
                                          ownership,          &
                                          comm)               &
      result (extrp)                                          &
      bind(c, name='PDM_extract_part_create')
        use iso_c_binding
        implicit none

        type(c_ptr)           :: extrp
        integer(c_int), value :: dim
        integer(c_int), value :: n_part_in
        integer(c_int), value :: n_part_out
        integer(c_int), value :: extract_kind
        integer(c_int), value :: split_dual_method
        integer(c_int), value :: compute_child_gnum
        integer(c_int), value :: ownership
        integer(c_int), value :: comm

      end function pdm_extract_part_create_c
    end interface

    if (compute_child_gnum) then
      c_compute_child_gnum = 1
    else
      c_compute_child_gnum = 0
    endif

    c_comm = PDM_MPI_Comm_f2c(comm)

    extrp = pdm_extract_part_create_c (dim,                  &
                                       n_part_in,            &
                                       n_part_out,           &
                                       extract_kind,         &
                                       split_dual_method,    &
                                       c_compute_child_gnum, &
                                       ownership,            &
                                       c_comm)

  end subroutine PDM_extract_part_create


  subroutine PDM_extract_part_part_set (extrp,                 &
                                        i_part,                &
                                        n_cell,                &
                                        n_face,                &
                                        n_edge,                &
                                        n_vtx,                 &
                                        cell_face_idx,         &
                                        cell_face,             &
                                        face_edge_idx,         &
                                        face_edge,             &
                                        edge_vtx,              &
                                        face_vtx_idx,          &
                                        face_vtx,              &
                                        cell_ln_to_gn,         &
                                        face_ln_to_gn,         &
                                        edge_ln_to_gn,         &
                                        vtx_ln_to_gn,          &
                                        vtx_coord)
    ! Set partition
    implicit none

    type(c_ptr), value            :: extrp            ! C pointer to PDM_extract_part_t instance
    integer, intent(in)           :: i_part           ! Partition identifer
    integer, intent(in)           :: n_cell           ! Number of cells
    integer, intent(in)           :: n_face           ! Number of faces
    integer, intent(in)           :: n_edge           ! Number of edges
    integer, intent(in)           :: n_vtx            ! Number of vertices
    integer(pdm_l_num_s), pointer :: cell_face_idx(:) ! Index for cell→face connectivity (shape = [n_cell+1])
    integer(pdm_l_num_s), pointer :: cell_face(:)     ! Cell→face connectivity (shape = [cell_face_idx(n_cell+1)])
    integer(pdm_l_num_s), pointer :: face_edge_idx(:) ! Index for face→edge connectivity (shape = [n_face+1])
    integer(pdm_l_num_s), pointer :: face_edge(:)     ! Face→edge connectivity (shape = [face_edge_idx(n_face+1)])
    integer(pdm_l_num_s), pointer :: edge_vtx(:)      ! Edge→vtx connectivity (shape = [2*n_edge])
    integer(pdm_l_num_s), pointer :: face_vtx_idx(:)  ! Index for face→vtx connectivity (shape = [n_face+1])
    integer(pdm_l_num_s), pointer :: face_vtx(:)      ! Face→vtx connectivity (shape = [face_vtx_idx(n_face+1)])
    integer(pdm_g_num_s), pointer :: cell_ln_to_gn(:) ! Cell global IDs (shape = [n_cell])
    integer(pdm_g_num_s), pointer :: face_ln_to_gn(:) ! Face global IDs (shape = [n_face])
    integer(pdm_g_num_s), pointer :: edge_ln_to_gn(:) ! Edge global IDs (shape = [n_edge])
    integer(pdm_g_num_s), pointer :: vtx_ln_to_gn(:)  ! Vertex global IDs (shape = [n_vtx])
    real(8),              pointer :: vtx_coord(:,:)   ! Vertex coordinates (shape = [3, n_vtx])
    type(c_ptr) :: c_cell_face_idx
    type(c_ptr) :: c_cell_face
    type(c_ptr) :: c_face_edge_idx
    type(c_ptr) :: c_face_edge
    type(c_ptr) :: c_edge_vtx
    type(c_ptr) :: c_face_vtx_idx
    type(c_ptr) :: c_face_vtx
    type(c_ptr) :: c_cell_ln_to_gn
    type(c_ptr) :: c_face_ln_to_gn
    type(c_ptr) :: c_edge_ln_to_gn
    type(c_ptr) :: c_vtx_ln_to_gn
    type(c_ptr) :: c_vtx_coord

    interface
      subroutine pdm_extract_part_part_set_c (extrp,                    &
                                              i_part,                   &
                                              n_cell,                   &
                                              n_face,                   &
                                              n_edge,                   &
                                              n_vtx,                    &
                                              cell_face_idx,            &
                                              cell_face,                &
                                              face_edge_idx,            &
                                              face_edge,                &
                                              edge_vtx,                 &
                                              face_vtx_idx,             &
                                              face_vtx,                 &
                                              cell_ln_to_gn,            &
                                              face_ln_to_gn,            &
                                              edge_ln_to_gn,            &
                                              vtx_ln_to_gn,             &
                                              vtx_coord)                &
      bind (c, name='PDM_extract_part_part_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part
        integer(c_int), value :: n_cell
        integer(c_int), value :: n_face
        integer(c_int), value :: n_edge
        integer(c_int), value :: n_vtx
        type(c_ptr),    value :: cell_face_idx
        type(c_ptr),    value :: cell_face
        type(c_ptr),    value :: face_edge_idx
        type(c_ptr),    value :: face_edge
        type(c_ptr),    value :: edge_vtx
        type(c_ptr),    value :: face_vtx_idx
        type(c_ptr),    value :: face_vtx
        type(c_ptr),    value :: cell_ln_to_gn
        type(c_ptr),    value :: face_ln_to_gn
        type(c_ptr),    value :: edge_ln_to_gn
        type(c_ptr),    value :: vtx_ln_to_gn
        type(c_ptr),    value :: vtx_coord

      end subroutine pdm_extract_part_part_set_c
    end interface

    c_cell_face_idx = C_NULL_PTR
    if (associated(cell_face_idx)) then
      c_cell_face_idx = c_loc(cell_face_idx)
    end if

    c_cell_face = C_NULL_PTR
    if (associated(cell_face)) then
      c_cell_face = c_loc(cell_face)
    end if

    c_face_edge_idx = C_NULL_PTR
    if (associated(face_edge_idx)) then
      c_face_edge_idx = c_loc(face_edge_idx)
    end if

    c_face_edge = C_NULL_PTR
    if (associated(face_edge)) then
      c_face_edge = c_loc(face_edge)
    end if

    c_edge_vtx = C_NULL_PTR
    if (associated(edge_vtx)) then
      c_edge_vtx = c_loc(edge_vtx)
    end if

    c_face_vtx_idx = C_NULL_PTR
    if (associated(face_vtx_idx)) then
      c_face_vtx_idx = c_loc(face_vtx_idx)
    end if

    c_face_vtx = C_NULL_PTR
    if (associated(face_vtx)) then
      c_face_vtx = c_loc(face_vtx)
    end if

    c_cell_ln_to_gn = C_NULL_PTR
    if (associated(cell_ln_to_gn)) then
      c_cell_ln_to_gn = c_loc(cell_ln_to_gn)
    end if

    c_face_ln_to_gn = C_NULL_PTR
    if (associated(face_ln_to_gn)) then
      c_face_ln_to_gn = c_loc(face_ln_to_gn)
    end if

    c_edge_ln_to_gn = C_NULL_PTR
    if (associated(edge_ln_to_gn)) then
      c_edge_ln_to_gn = c_loc(edge_ln_to_gn)
    end if

    c_vtx_ln_to_gn = C_NULL_PTR
    if (associated(vtx_ln_to_gn)) then
      c_vtx_ln_to_gn = c_loc(vtx_ln_to_gn)
    end if

    c_vtx_coord = C_NULL_PTR
    if (associated(vtx_coord)) then
      c_vtx_coord = c_loc(vtx_coord)
    end if

    call pdm_extract_part_part_set_c (extrp,                           &
                                      i_part,                          &
                                      n_cell,                          &
                                      n_face,                          &
                                      n_edge,                          &
                                      n_vtx,                           &
                                      c_cell_face_idx,                 &
                                      c_cell_face,                     &
                                      c_face_edge_idx,                 &
                                      c_face_edge,                     &
                                      c_edge_vtx,                      &
                                      c_face_vtx_idx,                  &
                                      c_face_vtx,                      &
                                      c_cell_ln_to_gn,                 &
                                      c_face_ln_to_gn,                 &
                                      c_edge_ln_to_gn,                 &
                                      c_vtx_ln_to_gn,                  &
                                      c_vtx_coord)

  end subroutine PDM_extract_part_part_set


  subroutine PDM_extract_part_n_group_set (extrp,      &
                                           bound_type, &
                                           n_group)
    ! Set number of groups
    implicit none

    type(c_ptr), value            :: extrp      ! C pointer to PDM_extract_part_t instance
    integer, intent(in)           :: bound_type ! Kind of group
    integer, intent(in)           :: n_group    ! Number of groups

    interface
      subroutine pdm_extract_part_n_group_set_c (extrp,      &
                                                 bound_type, &
                                                 n_group)    &
      bind (c, name='PDM_extract_part_n_group_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: bound_type
        integer(c_int), value :: n_group

      end subroutine pdm_extract_part_n_group_set_c
    end interface

    call pdm_extract_part_n_group_set_c (extrp,      &
                                         bound_type, &
                                         n_group)

  end subroutine PDM_extract_part_n_group_set


  subroutine PDM_extract_part_part_group_set (extrp,                 &
                                              i_part,                &
                                              i_group,               &
                                              bound_type,            &
                                              n_group_entity,        &
                                              group_entity,          &
                                              group_entity_ln_to_gn)
    ! Set partition group
    implicit none

    type(c_ptr), value            :: extrp                    ! C pointer to PDM_extract_part_t instance
    integer, intent(in)           :: i_part                   ! Partition identifier
    integer, intent(in)           :: i_group                  ! Group identifier
    integer, intent(in)           :: bound_type               ! Kind of group
    integer, intent(in)           :: n_group_entity           ! Number of entities in current group
    integer(pdm_l_num_s), pointer :: group_entity(:)          ! Local IDs of entities in group (shape = [n_group_entity])
    integer(pdm_g_num_s), pointer :: group_entity_ln_to_gn(:) ! Group-specific global IDs of entities in group (shape = [n_group_entity])

    interface
      subroutine pdm_extract_part_part_group_set_c (extrp,                 &
                                                    i_part,                &
                                                    i_group,               &
                                                    bound_type,            &
                                                    n_group_entity,        &
                                                    group_entity,          &
                                                    group_entity_ln_to_gn) &
      bind (c, name='PDM_extract_part_part_group_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part
        integer(c_int), value :: i_group
        integer(c_int), value :: bound_type
        integer(c_int), value :: n_group_entity
        type(c_ptr),    value :: group_entity
        type(c_ptr),    value :: group_entity_ln_to_gn

      end subroutine pdm_extract_part_part_group_set_c
    end interface

    call pdm_extract_part_part_group_set_c (extrp,                        &
                                            i_part,                       &
                                            i_group,                      &
                                            bound_type,                   &
                                            n_group_entity,               &
                                            c_loc(group_entity),          &
                                            c_loc(group_entity_ln_to_gn))

  end subroutine PDM_extract_part_part_group_set


  subroutine PDM_extract_part_part_nodal_set(extrp, &
                                             pmn)
    ! Set PDM_part_mesh_nodal_t
    implicit none

    type(c_ptr), intent(in) :: extrp ! C pointer to PDM_extract_part_t instance
    type(c_ptr), intent(in) :: pmn   ! C pointer to PDM_part_mesh_nodal_t instance

    interface
      subroutine PDM_extract_part_part_nodal_set_c(extrp, &
                                                   pmn)   &
      bind (c, name='PDM_extract_part_part_nodal_set')
        use iso_c_binding
        implicit none

        type(c_ptr), value :: extrp
        type(c_ptr), value :: pmn

      end subroutine PDM_extract_part_part_nodal_set_c
    end interface

    call PDM_extract_part_part_nodal_set_c(extrp, &
                                           pmn)

  end subroutine PDM_extract_part_part_nodal_set


  subroutine PDM_extract_part_compute(extrp)
    ! Compute extraction
    implicit none

    type(c_ptr) :: extrp ! C pointer to PDM_extract_part_t instance

    interface
      subroutine PDM_extract_part_compute_c(extrp) &
      bind (c, name='PDM_extract_part_compute')
        use iso_c_binding
        implicit none
        type(c_ptr), value :: extrp
      end subroutine PDM_extract_part_compute_c
    end interface

    call PDM_extract_part_compute_c(extrp)

  end subroutine PDM_extract_part_compute


  subroutine PDM_extract_part_group_get (extrp,                                &
                                         bound_type,                           &
                                         i_part,                               &
                                         i_group,                              &
                                         n_extract_group_entity,               &
                                         extract_group_entity,                 &
                                         extract_group_entity_ln_to_gn,        &
                                         extract_group_entity_parent_ln_to_gn, &
                                         ownership)
    ! Get partition group
    implicit none

    type(c_ptr), value            :: extrp                                   ! C pointer to PDM_extract_part_t instance
    integer, intent(in)           :: bound_type                              ! Kind of group
    integer, intent(in)           :: i_part                                  ! Partition identifier
    integer, intent(in)           :: i_group                                 ! Group identifier
    integer                       :: n_extract_group_entity                  ! Number of entities in current group
    integer(pdm_l_num_s), pointer :: extract_group_entity(:)                 ! Local IDs of entities in group (shape = [n_extract_group_entity])
    integer(pdm_g_num_s), pointer :: extract_group_entity_ln_to_gn(:)        ! Group-specific global IDs (in extraction) of entities in group (shape = [n_extract_group_entity])
    integer(pdm_g_num_s), pointer :: extract_group_entity_parent_ln_to_gn(:) ! Group-specific global IDs of entities in group (shape = [n_extract_group_entity])
    integer, intent(in)           :: ownership                               ! Ownership

    type(c_ptr)                   :: c_extract_group_entity
    type(c_ptr)                   :: c_extract_group_entity_ln_to_gn
    type(c_ptr)                   :: c_extract_group_entity_parent_ln_to_gn

    interface
      subroutine pdm_extract_part_group_get_c (extrp,                                &
                                               bound_type,                           &
                                               i_part,                               &
                                               i_group,                              &
                                               n_extract_group_entity,               &
                                               extract_group_entity,                 &
                                               extract_group_entity_ln_to_gn,        &
                                               extract_group_entity_parent_ln_to_gn, &
                                               ownership)                            &
      bind (c, name='PDM_extract_part_group_get')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: bound_type
        integer(c_int), value :: i_part
        integer(c_int), value :: i_group
        integer(c_int)        :: n_extract_group_entity
        type(c_ptr)           :: extract_group_entity
        type(c_ptr)           :: extract_group_entity_ln_to_gn
        type(c_ptr)           :: extract_group_entity_parent_ln_to_gn
        integer(c_int), value :: ownership

      end subroutine pdm_extract_part_group_get_c
    end interface

    c_extract_group_entity                 = C_NULL_PTR
    c_extract_group_entity_ln_to_gn        = C_NULL_PTR
    c_extract_group_entity_parent_ln_to_gn = C_NULL_PTR

    call pdm_extract_part_group_get_c (extrp,                                  &
                                       bound_type,                             &
                                       i_part,                                 &
                                       i_group,                                &
                                       n_extract_group_entity,                 &
                                       c_extract_group_entity,                 &
                                       c_extract_group_entity_ln_to_gn,        &
                                       c_extract_group_entity_parent_ln_to_gn, &
                                       ownership)

    call c_f_pointer(c_extract_group_entity, &
                     extract_group_entity,   &
                     [n_extract_group_entity])

    call c_f_pointer(c_extract_group_entity_ln_to_gn, &
                     extract_group_entity_ln_to_gn,   &
                     [n_extract_group_entity])

    call c_f_pointer(c_extract_group_entity_parent_ln_to_gn, &
                     extract_group_entity_parent_ln_to_gn,   &
                     [n_extract_group_entity])

  end subroutine PDM_extract_part_group_get


  subroutine PDM_extract_part_selected_lnum_set (extrp,        &
                                                 i_part_in,    &
                                                 n_entity,     &
                                                 extract_lnum, &
                                                 ownership)
    ! Select local entities to extract.
    !
    ! (Use only in PDM_EXTRACT_PART_KIND_LOCAL or PDM_EXTRACT_PART_KIND_REEQUILIBRATE mode)
    implicit none

    type(c_ptr), value                   :: extrp           ! C pointer to PDM_extract_part_t instance
    integer, intent(in)                  :: i_part_in       ! Partition identifier
    integer, intent(in)                  :: n_entity        ! Number of entities to extract
    integer(kind = PDM_l_num_s), pointer :: extract_lnum(:) ! Local IDs of entities to extract (shape = [n_entity])
    integer, intent(in)                  :: ownership       ! Ownership

    type(c_ptr)                          :: c_extract_lnum

    interface
      subroutine pdm_extract_part_selected_lnum_set_c (extrp,        &
                                                       i_part_in,    &
                                                       n_entity,     &
                                                       extract_lnum, &
                                                       ownership)    &
      bind (c, name='PDM_extract_part_selected_lnum_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part_in
        integer(c_int), value :: n_entity
        type(c_ptr),    value :: extract_lnum
        integer(c_int), value :: ownership

      end subroutine pdm_extract_part_selected_lnum_set_c
    end interface

    c_extract_lnum = C_NULL_PTR
    if (associated(extract_lnum)) then
      c_extract_lnum = c_loc(extract_lnum)
    end if

    call pdm_extract_part_selected_lnum_set_c (extrp,             &
                                               i_part_in,         &
                                               n_entity,          &
                                               c_extract_lnum,    &
                                               ownership)

  end subroutine PDM_extract_part_selected_lnum_set


  subroutine PDM_extract_part_target_set (extrp,           &
                                          i_part,          &
                                          n_target,        &
                                          target_gnum,     &
                                          target_location, &
                                          ownership)
    ! Set the target entities.
    !
    ! (Use only in PDM_EXTRACT_PART_KIND_FROM_TARGET mode)
    implicit none

    type(c_ptr), value                   :: extrp              ! C pointer to PDM_extract_part_t instance
    integer, intent(in)                  :: i_part             ! Partition identifier
    integer, intent(in)                  :: n_target           ! Number of target entities
    integer(kind = PDM_g_num_s), pointer :: target_gnum(:)     ! Global IDs of target entities (shape = [n_target])
    integer(kind = PDM_l_num_s), pointer :: target_location(:) ! Initial location of target entities (shape = [3*n_target] or *null()*)
    integer, intent(in)                  :: ownership          ! Ownership

    type(c_ptr)                          :: c_target_gnum
    type(c_ptr)                          :: c_target_location

    interface
      subroutine pdm_extract_part_target_set_c (extrp,           &
                                                i_part,          &
                                                n_target,        &
                                                target_gnum,     &
                                                target_location, &
                                                ownership)       &
      bind (c, name='PDM_extract_part_target_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part
        integer(c_int), value :: n_target
        type(c_ptr),    value :: target_gnum
        type(c_ptr),    value :: target_location
        integer(c_int), value :: ownership
      end subroutine pdm_extract_part_target_set_c
    end interface

    c_target_gnum = C_NULL_PTR
    if (associated(target_gnum)) then
      c_target_gnum = c_loc(target_gnum)
    end if

    c_target_location = C_NULL_PTR
    if (associated(target_location)) then
      c_target_location = c_loc(target_location)
    end if

    call pdm_extract_part_target_set_c (extrp,             &
                                        i_part,            &
                                        n_target,          &
                                        c_target_gnum,     &
                                        c_target_location, &
                                        ownership)

  end subroutine PDM_extract_part_target_set


  subroutine PDM_extract_part_n_entity_get (extrp,       &
                                            i_part_out,  &
                                            entity_type, &
                                            n_entity)
    ! Get the number of entities of a given type in extraction
    implicit none

    type(c_ptr), value                :: extrp       ! C pointer to PDM_extract_part_t instance
    integer, intent(in)               :: i_part_out  ! Partition identifier
    integer, intent(in)               :: entity_type ! Type of entity
    integer                           :: n_entity    ! Number of entities

    interface
      function pdm_extract_part_n_entity_get_c (extrp,       &
                                                i_part_out,  &
                                                entity_type) &
      result (n_entity)                                      &
      bind (c, name='PDM_extract_part_n_entity_get')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part_out
        integer(c_int), value :: entity_type
        integer(c_int)        :: n_entity

      end function pdm_extract_part_n_entity_get_c
    end interface

    n_entity = pdm_extract_part_n_entity_get_c (extrp,       &
                                                i_part_out,  &
                                                entity_type)

  end subroutine PDM_extract_part_n_entity_get


  subroutine PDM_extract_part_connectivity_get (extrp,             &
                                                i_part_out,        &
                                                connectivity_type, &
                                                n_entity,          &
                                                connect,           &
                                                connect_idx,       &
                                                ownership)
    ! Get connectivity in extraction
    implicit none

    type(c_ptr), value                   :: extrp             ! C pointer to PDM_extract_part_t instance
    integer, intent(in)                  :: i_part_out        ! Partition identifier
    integer, intent(in)                  :: connectivity_type ! Type of connectivity
    integer, intent(in)                  :: ownership         ! Ownership
    integer                              :: n_entity          ! Number of leading entities
    integer(kind = PDM_l_num_s), pointer :: connect(:)        ! Connectivity (shape = [connect_idx(n_entity+1)])
    integer(kind = PDM_l_num_s), pointer :: connect_idx(:)    ! Connectivity index (shape = [n_entity+1])

    type(c_ptr)                          :: c_connect
    type(c_ptr)                          :: c_connect_idx

    interface
      function pdm_extract_part_connectivity_get_c (extrp,             &
                                                    i_part_out,        &
                                                    connectivity_type, &
                                                    connect,           &
                                                    connect_idx,       &
                                                    ownership)         &
      result (n_entity)                                                &
      bind (c, name='PDM_extract_part_connectivity_get')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part_out
        integer(c_int), value :: connectivity_type
        integer(c_int), value :: ownership
        integer(c_int)        :: n_entity
        type(c_ptr)           :: connect
        type(c_ptr)           :: connect_idx

      end function pdm_extract_part_connectivity_get_c
    end interface

    c_connect     = C_NULL_PTR
    c_connect_idx = C_NULL_PTR

    n_entity = pdm_extract_part_connectivity_get_c (extrp,             &
                                                    i_part_out,        &
                                                    connectivity_type, &
                                                    c_connect,         &
                                                    c_connect_idx,     &
                                                    ownership)

    call c_f_pointer(c_connect_idx, &
                     connect_idx,   &
                     [n_entity+1])

    call c_f_pointer(c_connect, &
                     connect,   &
                     [connect_idx(n_entity+1)])

  end subroutine PDM_extract_part_connectivity_get


  subroutine PDM_extract_part_ln_to_gn_get (extrp,            &
                                            i_part_out,       &
                                            entity_type,      &
                                            n_entity,         &
                                            pentity_ln_to_gn, &
                                            ownership)
    ! Get global IDs of entities in extraction
    implicit none

    type(c_ptr), value                   :: extrp               ! C pointer to PDM_extract_part_t instance
    integer, intent(in)                  :: i_part_out          ! Partition identifier
    integer, intent(in)                  :: entity_type         ! Type of entity
    integer, intent(in)                  :: ownership           ! Ownership
    integer                              :: n_entity            ! Number of entities
    integer(kind = PDM_g_num_s), pointer :: pentity_ln_to_gn(:) ! Global IDs

    type(c_ptr)                          :: c_pentity_ln_to_gn

    interface
      function pdm_extract_part_ln_to_gn_get_c (extrp,            &
                                                i_part_out,       &
                                                entity_type,      &
                                                pentity_ln_to_gn, &
                                                ownership)        &
      result (n_entity)                                           &
      bind (c, name='PDM_extract_part_ln_to_gn_get')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part_out
        integer(c_int), value :: entity_type
        integer(c_int), value :: ownership
        integer(c_int)        :: n_entity
        type(c_ptr)           :: pentity_ln_to_gn

      end function pdm_extract_part_ln_to_gn_get_c
    end interface

    c_pentity_ln_to_gn = C_NULL_PTR

    n_entity = pdm_extract_part_ln_to_gn_get_c (extrp,              &
                                                i_part_out,         &
                                                entity_type,        &
                                                c_pentity_ln_to_gn, &
                                                ownership)

    call c_f_pointer(c_pentity_ln_to_gn, &
                     pentity_ln_to_gn,   &
                     [n_entity])

  end subroutine PDM_extract_part_ln_to_gn_get


  subroutine PDM_extract_part_parent_ln_to_gn_get (extrp,            &
                                                   i_part_out,       &
                                                   entity_type,      &
                                                   n_entity,         &
                                                   parent_ln_to_gn, &
                                                   ownership)
    ! Get parent global IDs of entities in extraction
    implicit none

    type(c_ptr), value                   :: extrp              ! C pointer to PDM_extract_part_t instance
    integer, intent(in)                  :: i_part_out         ! Partition identifier
    integer, intent(in)                  :: entity_type        ! Type of entity
    integer                              :: n_entity           ! Number of entities
    integer(kind = PDM_g_num_s), pointer :: parent_ln_to_gn(:) ! Parent global IDs
    integer, intent(in)                  :: ownership          ! Ownership

    type(c_ptr)                          :: c_parent_ln_to_gn = C_NULL_PTR

    interface
      function pdm_extract_part_parent_ln_to_gn_get_c (extrp,           &
                                                       i_part_out,      &
                                                       entity_type,     &
                                                       parent_ln_to_gn, &
                                                       ownership)       &
      result (n_entity)                                                 &
      bind (c, name='PDM_extract_part_parent_ln_to_gn_get')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part_out
        integer(c_int), value :: entity_type
        integer(c_int), value :: ownership
        integer(c_int)        :: n_entity
        type(c_ptr)           :: parent_ln_to_gn

      end function pdm_extract_part_parent_ln_to_gn_get_c
    end interface

    n_entity = pdm_extract_part_parent_ln_to_gn_get_c (extrp,              &
                                                       i_part_out,         &
                                                       entity_type,        &
                                                       c_parent_ln_to_gn,  &
                                                       ownership)

    call c_f_pointer(c_parent_ln_to_gn, &
                     parent_ln_to_gn,   &
                     [n_entity])

  end subroutine PDM_extract_part_parent_ln_to_gn_get


  subroutine PDM_extract_part_parent_lnum_get (extrp,              &
                                               i_part_out,         &
                                               entity_type,        &
                                               n_entity,           &
                                               parent_entity_lnum, &
                                               ownership)
    ! Get local IDs of parent entities.
    !
    ! (Use only in PDM_EXTRACT_PART_KIND_LOCAL mode)
    implicit none

    type(c_ptr), value                   :: extrp                 ! C pointer to PDM_extract_part_t instance
    integer, intent(in)                  :: i_part_out            ! Partition identifier
    integer, intent(in)                  :: entity_type           ! Type of entity
    integer, intent(in)                  :: ownership             ! Ownership
    integer                              :: n_entity              ! Number of entities
    integer(kind = PDM_l_num_s), pointer :: parent_entity_lnum(:) ! Parent local IDs

    type(c_ptr)                          :: c_parent_entity_lnum

    interface
      function pdm_extract_part_parent_lnum_get_c (extrp,              &
                                                   i_part_out,         &
                                                   entity_type,        &
                                                   parent_entity_lnum, &
                                                   ownership)          &
      result (n_entity)                                                &
      bind (c, name='PDM_extract_part_parent_lnum_get')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part_out
        integer(c_int), value :: entity_type
        integer(c_int), value :: ownership
        integer(c_int)        :: n_entity
        type(c_ptr)           :: parent_entity_lnum

      end function pdm_extract_part_parent_lnum_get_c
    end interface

    c_parent_entity_lnum = C_NULL_PTR

    n_entity = pdm_extract_part_parent_lnum_get_c (extrp,                &
                                                   i_part_out,           &
                                                   entity_type,          &
                                                   c_parent_entity_lnum, &
                                                   ownership)

    call c_f_pointer(c_parent_entity_lnum, &
                     parent_entity_lnum,   &
                     [n_entity])

  end subroutine PDM_extract_part_parent_lnum_get


  subroutine PDM_extract_part_part_to_part_get (extrp,       &
                                                entity_type, &
                                                ptp,         &
                                                ownership)
    ! Get the Part-to-Part instance for a given entity type
    implicit none

    type(c_ptr)                          :: extrp       ! C pointer to PDM_extract_part_t instance
    integer, intent(in)                  :: entity_type ! Type of entity
    type(c_ptr)                          :: ptp         ! C pointer to PDM_part_to_part_t instance
    integer, intent(in)                  :: ownership   ! Ownership

    interface

      subroutine PDM_extract_part_part_to_part_get_c (extrp,       &
                                                       entity_type, &
                                                       ptp,         &
                                                       ownership)   &
        bind (c, name='PDM_extract_part_part_to_part_get')
        use iso_c_binding
        implicit none

        type(c_ptr), value           :: extrp
        integer(c_int), value        :: entity_type
        type(c_ptr)                  :: ptp
        integer(c_int), value        :: ownership
      end subroutine PDM_extract_part_part_to_part_get_c

    end interface

    integer(c_int) :: c_entity_type
    integer(c_int) :: c_ownership

    c_entity_type = entity_type
    c_ownership   = ownership

    call PDM_extract_part_part_to_part_get_c (extrp, c_entity_type, ptp, c_ownership)

  end subroutine PDM_extract_part_part_to_part_get


  subroutine PDM_extract_part_part_to_part_group_get(extrp,      &
                                                     bound_type, &
                                                     i_group,    &
                                                     ptp,        &
                                                     ownership)
    ! Get the Part-to-Part instance for a given group
    implicit none

    type(c_ptr)                          :: extrp      ! C pointer to PDM_extract_part_t instance
    integer, intent(in)                  :: bound_type ! Type of group
    integer, intent(in)                  :: i_group    ! Group identifier
    type(c_ptr)                          :: ptp        ! C pointer to PDM_part_to_part_t instance
    integer, intent(in)                  :: ownership  ! Ownership

    interface
      subroutine PDM_extract_part_part_to_part_group_get_c(extrp,      &
                                                           bound_type, &
                                                           i_group,    &
                                                           ptp,        &
                                                           ownership)  &
      bind (c, name='PDM_extract_part_part_to_part_group_get')
        use iso_c_binding
        implicit none
        type(c_ptr), value           :: extrp
        integer(c_int), value        :: bound_type
        integer(c_int), value        :: i_group
        type(c_ptr)                  :: ptp
        integer(c_int), value        :: ownership
      end subroutine PDM_extract_part_part_to_part_group_get_c
    end interface

    call PDM_extract_part_part_to_part_group_get_c(extrp,      &
                                                   bound_type, &
                                                   i_group,    &
                                                   ptp,        &
                                                   ownership)

  end subroutine PDM_extract_part_part_to_part_group_get


  subroutine PDM_extract_part_vtx_coord_get (extrp,            &
                                             i_part_out,       &
                                             n_entity,         &
                                             pvtx_coord,       &
                                             ownership)
    ! Get vertex coordinates in extraction
    implicit none

    type(c_ptr), value                   :: extrp           ! C pointer to PDM_extract_part_t instance
    integer, intent(in)                  :: i_part_out      ! Partition identifier
    integer, intent(in)                  :: ownership       ! Ownership
    integer                              :: n_entity        ! Number of vertices
    real(8), pointer                     :: pvtx_coord(:,:) ! Coordinates of vertices (shape = [3, n_entity])

    type(c_ptr)                          :: c_pvtx_coord = C_NULL_PTR

    interface
      function pdm_extract_part_vtx_coord_get_c (extrp,            &
                                                 i_part_out,       &
                                                 pvtx_coord,       &
                                                 ownership)        &
      result (n_entity)                                            &
      bind (c, name='PDM_extract_part_vtx_coord_get')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part_out
        integer(c_int), value :: ownership
        integer(c_int)        :: n_entity
        type(c_ptr)           :: pvtx_coord

      end function pdm_extract_part_vtx_coord_get_c
    end interface

    n_entity = pdm_extract_part_vtx_coord_get_c (extrp,              &
                                                 i_part_out,         &
                                                 c_pvtx_coord,       &
                                                 ownership)

    call c_f_pointer(c_pvtx_coord, &
                     pvtx_coord,   &
                     [3,n_entity])

  end subroutine PDM_extract_part_vtx_coord_get


  subroutine PDM_extract_part_part_mesh_nodal_get(extrp,     &
                                                  pmn,       &
                                                  ownership)
    ! Retrieve the extracted PDM_part_mesh_nodal_t
    implicit none

    type(c_ptr), intent(in)  :: extrp     ! C pointer to PDM_extract_part_t instance
    type(c_ptr), intent(out) :: pmn       ! C pointer to PDM_part_mesh_nodal_t instance
    integer,     intent(in)  :: ownership ! Ownership

    interface
      subroutine PDM_extract_part_part_mesh_nodal_get_c(extrp,     &
                                                        pmn,       &
                                                        ownership) &
      bind (c, name='PDM_extract_part_part_mesh_nodal_get')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        type(c_ptr)           :: pmn
        integer(c_int), value :: ownership

      end subroutine PDM_extract_part_part_mesh_nodal_get_c
    end interface

    call PDM_extract_part_part_mesh_nodal_get_c(extrp, &
                                                pmn,   &
                                                ownership)

  end subroutine PDM_extract_part_part_mesh_nodal_get


  subroutine PDM_extract_part_partial_free(extrp)
    ! Free all resulting data if not owner
    !
    ! (It is not necessary to call this subroutine before calling PDM_extract_part_free)
    implicit none

    type(c_ptr) :: extrp ! C pointer to PDM_extract_part_t instance

    interface
      subroutine PDM_extract_part_partial_free_c(extrp) &
      bind (c, name='PDM_extract_part_partial_free')
        use iso_c_binding
        implicit none
        type(c_ptr), value :: extrp
      end subroutine PDM_extract_part_partial_free_c
    end interface

    call PDM_extract_part_partial_free_c(extrp)

  end subroutine PDM_extract_part_partial_free


  subroutine PDM_extract_part_free(extrp)
    ! Free the structure
    implicit none

    type(c_ptr) :: extrp ! C pointer to PDM_extract_part_t instance

    interface
      subroutine PDM_extract_part_free_c(extrp) &
      bind (c, name='PDM_extract_part_free')
        use iso_c_binding
        implicit none
        type(c_ptr), value :: extrp
      end subroutine PDM_extract_part_free_c
    end interface

    call PDM_extract_part_free_c(extrp)

  end subroutine PDM_extract_part_free


end module pdm_extract_part
