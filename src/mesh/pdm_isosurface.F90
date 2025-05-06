!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2020  ONERA
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

module pdm_isosurface

  use pdm
  use iso_c_binding

  implicit none

  integer(c_int), parameter :: PDM_ISO_SURFACE_KIND_FIELD         = 0 ! Use point octree
  integer(c_int), parameter :: PDM_ISO_SURFACE_KIND_PLANE        = 1 ! Use bounding-box tree
  integer(c_int), parameter :: PDM_ISO_SURFACE_KIND_SPHERE = 2 ! Locate all target points
  integer(c_int), parameter :: PDM_ISO_SURFACE_KIND_ELLIPSE = 3 ! Locate all target points
  integer(c_int), parameter :: PDM_ISO_SURFACE_KIND_QUADRIC = 4 ! Locate all target points
  integer(c_int), parameter :: PDM_ISO_SURFACE_KIND_HEART = 5 ! Locate all target points
  integer(c_int), parameter :: PDM_ISO_SURFACE_KIND_FUNCTION = 6 ! Locate all target points
  integer(c_int), parameter :: PDM_ISO_SURFACE_KIND_MAX = 7 ! Locate all target points


  interface

    subroutine PDM_isosurface_tolerance_set(isos,&
                                            tol) &
    bind(c, name = 'PDM_isosurface_set_tolerance')

      use iso_c_binding
      implicit none

      type(c_ptr),    value :: isos
      real(c_double), value :: tol

    end subroutine PDM_isosurface_tolerance_set


    subroutine PDM_isosurface_n_part_set (isos,   &
                                          n_part) &
    bind (c, name = 'PDM_isosurface_n_part_set')
      ! Set the number of partitions of a point cloud
      use iso_c_binding
      implicit none

      type (c_ptr),   value :: isos          ! C pointer to PDM_mesh_location_t object
      integer(c_int), value :: n_part        ! Number of partitions

    end subroutine PDM_isosurface_n_part_set


    subroutine PDM_isosurface_n_group_set (isos,        &
                                           entity_type, &
                                           n_group)     &
    bind(c, name = 'PDM_isosurface_n_group_set')
      ! Set the number of group
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: isos        ! PDM_isosurface_t instance
      integer(c_int), value :: entity_type ! Type of mesh entity
      integer(c_int), value :: n_group     ! Number of group

    end subroutine PDM_isosurface_n_group_set


    subroutine PDM_isosurface_part_mesh_set (isos,  &
                                             pmesh) &
    bind(c, name = 'PDM_isosurface_part_mesh_set')
      ! Set partioned mesh
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: isos  ! PDM_isosurface_t instance
      type(c_ptr),    value :: pmesh ! PDM_part_mesh_t instance

    end subroutine PDM_isosurface_part_mesh_set


    subroutine PDM_isosurface_part_mesh_nodal_set (isos, &
                                                   pmn)  &
    bind(c, name = 'PDM_isosurface_part_mesh_nodal_set')
      ! Set nodal mesh
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: isos ! PDM_isosurface_t instance
      type(c_ptr),    value :: pmn  ! PDM_part_mesh_nodal_t instance

    end subroutine PDM_isosurface_part_mesh_nodal_set


    subroutine PDM_isosurface_dmesh_set (isos,  &
                                         dmesh) &
    bind(c, name = 'PDM_isosurface_dmesh_set')

      use iso_c_binding
      implicit none

      type(c_ptr),    value :: isos
      type(c_ptr),    value :: dmesh

    end subroutine PDM_isosurface_dmesh_set


    subroutine PDM_isosurface_dmesh_nodal_set (isos, &
                                               dmn)  &
    bind(c, name = 'PDM_isosurface_dmesh_nodal_set')

      use iso_c_binding
      implicit none

      type(c_ptr),    value :: isos
      type(c_ptr),    value :: dmn

    end subroutine PDM_isosurface_dmesh_nodal_set


    ! subroutine PDM_isosurface_field_function_set (isos, &
    !                                               id_isosurface, &
    !                                               func)
    !   bind(c, name = 'PDM_isosurface_field_function_set')

    !   use iso_c_binding

    !   implicit none

    !   type(c_ptr),    value :: isos
    !   integer(c_int), value :: id_isosurface
    !   type(c_ptr),    value :: func

    ! end subroutine PDM_isosurface_field_function_set


    subroutine PDM_isosurface_redistribution_set (isos,         &
                                                  extract_kind, &
                                                  part_method)  &
    bind(c, name = 'PDM_isosurface_redistribution_set')

      use iso_c_binding
      implicit none

      type(c_ptr),    value :: isos
      integer(c_int), value :: extract_kind
      integer(c_int), value :: part_method

    end subroutine PDM_isosurface_redistribution_set


    subroutine PDM_isosurface_reset (isos,          &
                                     id_isosurface) &
    bind(c, name = 'PDM_isosurface_reset')

      use iso_c_binding
      implicit none

      type(c_ptr),    value :: isos
      integer(c_int), value :: id_isosurface

    end subroutine PDM_isosurface_reset


    subroutine PDM_isosurface_n_part_out_set (isos,       &
                                              n_part_out) &
    bind(c, name = 'PDM_isosurface_n_part_out_set')

      use iso_c_binding
      implicit none

      type(c_ptr),    value :: isos
      integer(c_int), value :: n_part_out

    end subroutine PDM_isosurface_n_part_out_set


    subroutine PDM_isosurface_compute (isos,          &
                                       id_isosurface) &
    bind(c, name = 'PDM_isosurface_compute')

      use iso_c_binding
      implicit none

      type(c_ptr),    value :: isos
      integer(c_int), value :: id_isosurface

    end subroutine PDM_isosurface_compute


    subroutine PDM_isosurface_dump_times (isos) &
    bind(c, name = 'PDM_isosurface_dump_times')

      use iso_c_binding
      implicit none

      type(c_ptr),    value :: isos

    end subroutine PDM_isosurface_dump_times


    subroutine PDM_isosurface_part_to_part_enable (isos, &
                                                   id_isosurface, &
                                                   entity_type, &
                                                   unify_parent_info) &
    bind(c, name='PDM_isosurface_part_to_part_enable')
      ! Enable construction of a communication graph between source mesh entities and iso-surface entities.
      use iso_c_binding
      implicit none

      type(c_ptr)    :: isos
      integer(c_int) :: id_isosurface
      integer(c_int) :: entity_type
      integer(c_int) :: unify_parent_info

    end subroutine PDM_isosurface_part_to_part_enable


    subroutine PDM_isosurface_part_to_part_get (isos, &
                                                id_isosurface, &
                                                entity_type, &
                                                ptp, &
                                                ownership) &
    bind(c, name='PDM_isosurface_part_to_part_get')
      ! Get \ref PDM_part_to_part_t instance to exchange data between source mesh entities and iso-surface entities.
      use iso_c_binding
      implicit none

      type(c_ptr)    :: isos
      integer(c_int) :: id_isosurface
      integer(c_int) :: entity_type
      type(c_ptr)    :: ptp
      integer(c_int) :: ownership

    end subroutine PDM_isosurface_part_to_part_get


    subroutine PDM_isosurface_free (isos)&
     bind (c, name = 'PDM_isosurface_free')

      use iso_c_binding
      implicit none

      type (c_ptr), value :: isos

    end subroutine PDM_isosurface_free

  end interface


  contains

  subroutine PDM_isosurface_create (comm,           &
                                    mesh_dimension, &
                                    isos)

    ! Create a PDM_isosurface_t instance
    use iso_c_binding
    implicit none

    integer, intent(in) :: comm           ! MPI communicator
    integer, intent(in) :: mesh_dimension ! Dimension of source mesh (2 or 3)
    type(c_ptr)         :: isos           ! isosurface instance

    integer(c_int)      :: c_comm

    interface
      function PDM_isosurface_create_cf (comm,           &
                                         mesh_dimension) &
                                         result(isos)    &
        bind (c, name = 'PDM_isosurface_create')

        use iso_c_binding
        implicit none

        integer(c_int), value :: comm
        integer(c_int), value :: mesh_dimension
        type(c_ptr)           :: isos

      end function PDM_isosurface_create_cf
    end interface

    c_comm = PDM_MPI_Comm_f2c(comm)

    isos = PDM_isosurface_create_cf(c_comm, &
                                    mesh_dimension)

  end subroutine PDM_isosurface_create



  subroutine PDM_isosurface_pconnectivity_set (isos,              &
                                               i_part,            &
                                               connectivity_type, &
                                               n_entity,          &
                                               connect_idx,       &
                                               connect)

    ! Set connectivity
    use iso_c_binding
    implicit none

    type(c_ptr), intent(in) :: isos              ! PDM_isosurface_t instance
    integer,     intent(in) :: i_part            ! Partition identifier
    integer,     intent(in) :: connectivity_type ! Type of connectivity
    integer,     intent(in) :: n_entity          ! Local number of leading entities
    integer, pointer        :: connect_idx(:)    ! Index for connectivity
    integer, pointer        :: connect(:)        ! Connectivity

    type(c_ptr) :: c_connect_idx
    type(c_ptr) :: c_connect

    interface
      subroutine PDM_isosurface_pconnectivity_set_cf (isos,              &
                                                      i_part,            &
                                                      connectivity_type, &
                                                      n_entity,          &
                                                      connect_idx,       &
                                                      connect)           &
      bind(c, name = 'PDM_isosurface_pconnectivity_set')

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: isos
        integer(c_int), value :: i_part
        integer(c_int), value :: connectivity_type
        integer(c_int), value :: n_entity
        type(c_ptr),    value :: connect_idx
        type(c_ptr),    value :: connect

      end subroutine PDM_isosurface_pconnectivity_set_cf
    end interface

    c_connect_idx = C_NULL_PTR
    if (associated(connect_idx)) then
      c_connect_idx = c_loc(connect_idx)
    endif

    c_connect = C_NULL_PTR
    if (associated(connect)) then
      c_connect = c_loc(connect)
    endif

    call PDM_isosurface_pconnectivity_set_cf(isos,              &
                                             i_part,            &
                                             connectivity_type, &
                                             n_entity,          &
                                             c_connect_idx,     &
                                             c_connect)

  end subroutine PDM_isosurface_pconnectivity_set


  subroutine PDM_isosurface_pvtx_coord_set(isos,   &
                                           i_part, &
                                           n_vtx,  &
                                           vtx_coord)

    ! Set vertex coordinates
    use iso_c_binding
    implicit none

    type(c_ptr), intent(in) :: isos           ! PDM_isosurface_t instance
    integer,     intent(in) :: i_part         ! Partition identifier
    integer,     intent(in) :: n_vtx          ! Local number of vertices
    real(8), pointer        :: vtx_coord(:,:) ! Vertex coordinates (shape = [3, n_vtx])

    type(c_ptr) :: c_vtx_coord

    interface
      subroutine PDM_isosurface_pvtx_coord_set_cf (isos,     &
                                                  i_part,    &
                                                  n_vtx,     &
                                                  vtx_coord) &
      bind(c, name = 'PDM_isosurface_pvtx_coord_set')

        use iso_c_binding

        implicit none

        type(c_ptr),    value :: isos
        integer(c_int), value :: i_part
        integer(c_int), value :: n_vtx
        type(c_ptr),    value :: vtx_coord

      end subroutine PDM_isosurface_pvtx_coord_set_cf
    end interface

    c_vtx_coord = C_NULL_PTR
    if (associated(vtx_coord)) then
      c_vtx_coord = c_loc(vtx_coord)
    endif

    call PDM_isosurface_pvtx_coord_set_cf(isos,   &
                                          i_part, &
                                          n_vtx,  &
                                          c_vtx_coord)

  end subroutine PDM_isosurface_pvtx_coord_set

  subroutine PDM_isosurface_ln_to_gn_set(isos,        &
                                         i_part,      &
                                         entity_type, &
                                         ln_to_gn)

    ! Set global ids
    use iso_c_binding
    implicit none

    type(c_ptr), intent(in)            :: isos           ! PDM_isosurface_t instance
    integer,     intent(in)            :: i_part         ! Partition identifier
    integer,     intent(in)            :: entity_type    ! Type of mesh entity
    integer(kind=pdm_g_num_s), pointer :: ln_to_gn(:)    ! Global Ids

    type(c_ptr) :: c_ln_to_gn

    interface
      subroutine PDM_isosurface_ln_to_gn_set_cf (isos,        &
                                                 i_part,      &
                                                 entity_type, &
                                                 ln_to_gn)    &
      bind(c, name = 'PDM_isosurface_ln_to_gn_set')

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: isos
        integer(c_int), value :: i_part
        integer(c_int), value :: entity_type
        type(c_ptr),    value :: ln_to_gn

      end subroutine PDM_isosurface_ln_to_gn_set_cf
    end interface

    c_ln_to_gn = C_NULL_PTR
    if (associated(ln_to_gn)) then
      c_ln_to_gn = c_loc(ln_to_gn)
    endif    

    call PDM_isosurface_ln_to_gn_set_cf(isos,        &
                                        i_part,      &
                                        entity_type, &
                                        c_ln_to_gn)

  end subroutine PDM_isosurface_ln_to_gn_set


  subroutine PDM_isosurface_pgroup_set(isos,             &
                                       i_part,           &
                                       entity_type,      &
                                       group_entity_idx, &
                                       group_entity,     &
                                       group_entity_ln_to_gn)

    ! Set group description
    use iso_c_binding
    implicit none

    type(c_ptr), intent(in)            :: isos                     ! PDM_isosurface_t instance
    integer,     intent(in)            :: i_part                   ! Partition identifier
    integer,     intent(in)            :: entity_type              ! Type of mesh entity
    integer,                   pointer :: group_entity_idx(:)      ! Index for group->entity connectivity (size = n_group + 1)
    integer,                   pointer :: group_entity(:)          ! Group->entity connectivity (size = group_entity_idx[n_group])
    integer(kind=pdm_g_num_s), pointer :: group_entity_ln_to_gn(:) ! Group->entity connectivity (group-specific global ids, size = group_entity_idx[n_group])

    type(c_ptr) :: c_group_entity_idx
    type(c_ptr) :: c_group_entity
    type(c_ptr) :: c_group_entity_ln_to_gn

    interface

      subroutine PDM_isosurface_pgroup_set_cf (isos,                  &
                                               i_part,                &
                                               entity_type,           &
                                               group_entity_idx,      &
                                               group_entity,          &
                                               group_entity_ln_to_gn) &
      bind(c, name = 'PDM_isosurface_pgroup_set')

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: isos
        integer(c_int), value :: i_part
        integer(c_int), value :: entity_type
        type(c_ptr),    value :: group_entity_idx
        type(c_ptr),    value :: group_entity
        type(c_ptr),    value :: group_entity_ln_to_gn

      end subroutine PDM_isosurface_pgroup_set_cf
    end interface

    c_group_entity_idx = C_NULL_PTR
    if (associated(group_entity_idx)) then
      c_group_entity_idx = c_loc(group_entity_idx)
    endif    

    c_group_entity = C_NULL_PTR
    if (associated(group_entity)) then
      c_group_entity = c_loc(group_entity)
    endif    

    c_group_entity_ln_to_gn = C_NULL_PTR
    if (associated(group_entity_ln_to_gn)) then
      c_group_entity_ln_to_gn = c_loc(group_entity_ln_to_gn)
    endif            

    call PDM_isosurface_pgroup_set_cf(isos,               &
                                      i_part,             &
                                      entity_type,        &
                                      c_group_entity_idx, &
                                      c_group_entity,     &
                                      c_group_entity_ln_to_gn)

  end subroutine PDM_isosurface_pgroup_set


  subroutine PDM_isosurface_dconnectivity_set(isos,              &
                                              connectivity_type, &
                                              dconnect_idx,      &
                                              dconnect)


    ! Set block-distributed connectivity
    use iso_c_binding
    implicit none

    type(c_ptr), intent(in)            :: isos              ! PDM_isosurface_t instance
    integer,     intent(in)            :: connectivity_type ! Type of connectivity
    integer,                   pointer :: dconnect_idx(:)   ! Index for connectivity
    integer(kind=pdm_g_num_s), pointer :: dconnect(:)       ! Connectivity

    type(c_ptr) :: c_dconnect_idx
    type(c_ptr) :: c_dconnect

    interface
      subroutine PDM_isosurface_dconnectivity_set_cf (isos,              &
                                                      connectivity_type, &
                                                      dconnect_idx,      &
                                                      dconnect)          &
      bind(c, name = 'PDM_isosurface_dconnectivity_set')

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: isos
        integer(c_int), value :: connectivity_type
        type(c_ptr),    value :: dconnect_idx
        type(c_ptr),    value :: dconnect

      end subroutine PDM_isosurface_dconnectivity_set_cf
    end interface

    c_dconnect_idx = C_NULL_PTR
    if (associated(dconnect_idx)) then
      c_dconnect_idx = c_loc(dconnect_idx)
    endif    

    c_dconnect = C_NULL_PTR
    if (associated(dconnect)) then
      c_dconnect = c_loc(dconnect)
    endif    

    call PDM_isosurface_dconnectivity_set_cf(isos,              &
                                             connectivity_type, &
                                             c_dconnect_idx,    &
                                             c_dconnect)

  end subroutine PDM_isosurface_dconnectivity_set

  subroutine PDM_isosurface_dvtx_coord_set(isos, &
                                           dvtx_coord)

    ! Set block-distributed vertex coordinates
    use iso_c_binding
    implicit none

    type(c_ptr), intent(in) :: isos
    real(8), pointer :: dvtx_coord(:,:)

    type(c_ptr) :: c_dvtx_coord

    interface
      subroutine PDM_isosurface_dvtx_coord_set_cf (isos,       &
                                                   dvtx_coord) &
      bind(c, name = 'PDM_isosurface_dvtx_coord_set')

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: isos
        type(c_ptr),    value :: dvtx_coord

      end subroutine PDM_isosurface_dvtx_coord_set_cf
    end interface

    c_dvtx_coord = C_NULL_PTR
    if (associated(dvtx_coord)) then
      c_dvtx_coord = c_loc(dvtx_coord)
    endif    

    call PDM_isosurface_dvtx_coord_set_cf(isos, &
                                          c_dvtx_coord)

  end subroutine PDM_isosurface_dvtx_coord_set

  subroutine PDM_isosurface_distrib_set(isos,        &
                                        entity_type, &
                                        distrib)
    ! Set entity block distribution index 
    use iso_c_binding
    implicit none

    type(c_ptr), intent(in) :: isos ! PDM_isosurface_t instance
    integer, intent(in) :: entity_type ! Type of mesh entity
    integer(kind=pdm_g_num_s), pointer :: distrib(:) ! Block-distribution (size = nrank + 1)

    type(c_ptr) :: c_distrib

    interface
      subroutine PDM_isosurface_distrib_set_cf (isos,        &
                                                entity_type, &
                                                distrib)     &
      bind(c, name = 'PDM_isosurface_distrib_set')

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: isos
        integer(c_int), value :: entity_type
        type(c_ptr),    value :: distrib

      end subroutine PDM_isosurface_distrib_set_cf
    end interface

    c_distrib = C_NULL_PTR
    if (associated(distrib)) then
      c_distrib = c_loc(distrib)
    endif   

    call PDM_isosurface_distrib_set_cf(isos,        &
                                       entity_type, &
                                       c_distrib)     

  end subroutine PDM_isosurface_distrib_set

  subroutine PDM_isosurface_dgroup_set(isos,              &
                                       entity_type,       &
                                       dgroup_entity_idx, &
                                       dgroup_entity)
    ! Set block-distributed group description
    use iso_c_binding
    implicit none

    type(c_ptr), intent(in)            :: isos                 ! PDM_isosurface_t instance
    integer,     intent(in)            :: entity_type          ! Type of mesh entity
    integer,                   pointer :: dgroup_entity_idx(:) ! Index for group->entity connectivity (size = ``n_group`` + 1)
    integer(kind=pdm_g_num_s), pointer :: dgroup_entity(:)     ! Group->entity connectivity (size = ``dgroup_entity_idx(n_group+1)``)

    type(c_ptr) :: c_dgroup_entity_idx
    type(c_ptr) :: c_dgroup_entity

    interface
      subroutine PDM_isosurface_dgroup_set_cf (isos,              &
                                               entity_type,       &
                                               dgroup_entity_idx, &
                                               dgroup_entity)     &
      bind(c, name = 'PDM_isosurface_dgroup_set')

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: isos
        integer(c_int), value :: entity_type
        type(c_ptr),    value :: dgroup_entity_idx
        type(c_ptr),    value :: dgroup_entity

      end subroutine PDM_isosurface_dgroup_set_cf
    end interface

    c_dgroup_entity_idx = C_NULL_PTR
    if (associated(dgroup_entity_idx)) then
      c_dgroup_entity_idx = c_loc(dgroup_entity_idx)
    endif   

    c_dgroup_entity = C_NULL_PTR
    if (associated(dgroup_entity)) then
      c_dgroup_entity = c_loc(dgroup_entity)
    endif  

    call PDM_isosurface_dgroup_set_cf (isos,                &
                                       entity_type,         &
                                       c_dgroup_entity_idx, &
                                       c_dgroup_entity)     

  end subroutine PDM_isosurface_dgroup_set


  subroutine PDM_isosurface_add (isos,        &
                                 kind,          &
                                 n_isovalues,   &
                                 isovalues,     &
                                 id_isosurface)
    ! Add a requested set of iso-surfaces
    use iso_c_binding
    implicit none

    type(c_ptr), intent(in) :: isos           ! PDM_isosurface_t instance
    integer,     intent(in) :: kind           ! Iso-surface kind (discrete field, slice equation or function pointer)
    integer,     intent(in) :: n_isovalues    ! Number os iso-values to capture
    real(8), pointer        :: isovalues(:)   ! Iso-values to capture (size = ``n_isovalues``)
    integer, intent(out)    :: id_isosurface  ! Iso-surface identifier

    type(c_ptr) :: c_isovalues

    interface
      function PDM_isosurface_add_cf (isos,                 &
                                      kind,                 &
                                      n_isovalues,          &
                                      isovalues)            &
                                      result(id_isosurface) &
      bind(c, name = 'PDM_isosurface_add')

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: isos
        integer(c_int), value :: kind
        integer(c_int), value :: n_isovalues
        type(c_ptr),    value :: isovalues
        integer(c_int)        :: id_isosurface

      end function PDM_isosurface_add_cf
    end interface

    c_isovalues = C_NULL_PTR
    if (associated(isovalues)) then
      c_isovalues = c_loc(isovalues)
    endif  

    id_isosurface =  PDM_isosurface_add_cf (isos,        &
                                            kind,        &
                                            n_isovalues, &
                                            c_isovalues)

  end subroutine PDM_isosurface_add

  subroutine PDM_isosurface_isovalue_set (isos,          &
                                          id_isosurface, &
                                          n_isovalues,   &
                                          isovalues)
    ! Reset isovalues of given isosurface
    use iso_c_binding
    implicit  none

    type(c_ptr), intent(in) :: isos          ! PDM_isosurface_t instance
    integer,     intent(in) :: id_isosurface ! Iso-surface identifier
    integer,     intent(in) :: n_isovalues   ! Number of iso-values to capture
    real(8), pointer        :: isovalues(:)  ! Iso-values to capture (size = ``n_isovalues``)

    type(c_ptr) :: c_isovalues

    interface
      subroutine PDM_isosurface_isovalues_set_cf (isos,          &
                                                  id_isosurface, &
                                                  n_isovalues,   &
                                                  isovalues)     &
      bind(c, name = 'PDM_isosurface_isovalues_set')

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: isos
        integer(c_int), value :: id_isosurface
        integer(c_int), value :: n_isovalues
        type(c_ptr),    value :: isovalues

      end subroutine PDM_isosurface_isovalues_set_cf
    end interface

    c_isovalues = C_NULL_PTR
    if (associated(isovalues)) then
      c_isovalues = c_loc(isovalues)
    endif  

    call PDM_isosurface_isovalues_set_cf (isos,          &
                                          id_isosurface, &
                                          n_isovalues,   &
                                          c_isovalues)    

  end subroutine PDM_isosurface_isovalue_set

  subroutine PDM_isosurface_equation_set (isos,          &
                                          id_isosurface, &
                                          coeff)
    ! Set source field equation
    use iso_c_binding
    implicit none

    type(c_ptr), intent(in) :: isos          ! PDM_isosurface_t instance
    integer,     intent(in) :: id_isosurface ! Iso-surface identifier
    real(8), pointer        :: coeff(:)      ! Equation coefficients

    type(c_ptr) :: c_coeff

    interface
      subroutine PDM_isosurface_equation_set_cf (isos,          &
                                                 id_isosurface, &
                                                 coeff)         &
      bind(c, name = 'PDM_isosurface_equation_set')

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: isos
        integer(c_int), value :: id_isosurface
        type(c_ptr),    value :: coeff

      end subroutine PDM_isosurface_equation_set_cf
    end interface

    c_coeff = C_NULL_PTR
    if (associated(coeff)) then
      c_coeff = c_loc(coeff)
    endif      

    call PDM_isosurface_equation_set_cf (isos,          &
                                         id_isosurface, &
                                         c_coeff)

  end subroutine PDM_isosurface_equation_set


  subroutine PDM_isosurface_pfield_set (isos,          &
                                        id_isosurface, &
                                        i_part,        &
                                        field)
    ! Set field values.
    use iso_c_binding
    implicit none

    type(c_ptr), intent(in) :: isos          ! PDM_isosurface_t instance
    integer,     intent(in) :: id_isosurface ! Iso-surface identifier
    integer,     intent(in) :: i_part        ! Partition identifier
    real(8), pointer        :: field(:)      ! Field values (size = ``n_vtx``)

    type(c_ptr) :: c_field

    interface
      subroutine PDM_isosurface_pfield_set_cf (isos,          &
                                               id_isosurface, &
                                               i_part,        &
                                               field)         &
      bind(c, name = 'PDM_isosurface_pfield_set')

        use iso_c_binding

        implicit none

        type(c_ptr),    value :: isos
        integer(c_int), value :: id_isosurface
        integer(c_int), value :: i_part
        type(c_ptr),    value :: field

      end subroutine PDM_isosurface_pfield_set_cf
    end interface    

    c_field = C_NULL_PTR
    if (associated(field)) then
      c_field = c_loc(field)
    endif      

    call PDM_isosurface_pfield_set_cf (isos,          &
                                       id_isosurface, &
                                       i_part,        &
                                       c_field)

  end subroutine PDM_isosurface_pfield_set


  subroutine PDM_isosurface_dfield_set (isos,          &
                                        id_isosurface, &
                                        dfield)
    ! Set block-distributed field values.
    use iso_c_binding
    implicit none

    type(c_ptr), intent(in) :: isos      ! PDM_isosurface_t instance
    integer, intent(in) :: id_isosurface ! Iso-surface identifier
    real(8), pointer :: dfield(:)        ! Field values (size = ``dn_vtx``)

    type(c_ptr) :: c_dfield

    interface
      subroutine PDM_isosurface_dfield_set_cf (isos,          &
                                               id_isosurface, &
                                               dfield)        &
      bind(c, name = 'PDM_isosurface_dfield_set')

        use iso_c_binding

        implicit none

        type(c_ptr),    value :: isos
        integer(c_int), value :: id_isosurface
        type(c_ptr),    value :: dfield

      end subroutine PDM_isosurface_dfield_set_cf
    end interface

    c_dfield = C_NULL_PTR
    if (associated(dfield)) then
      c_dfield = c_loc(dfield)
    endif      

    call PDM_isosurface_dfield_set_cf (isos,          &
                                       id_isosurface, &
                                       c_dfield)

  end subroutine PDM_isosurface_dfield_set


  function PDM_isosurface_pconnectivity_get (isos,              &
                                             id_isosurface,     &
                                             i_part,            &
                                             connectivity_type, &
                                             connect_idx,       &
                                             connect,           &
                                             ownership)         &
                                             result(n_entity)
    ! Get the iso-surfaces mesh connectivity
    use iso_c_binding
    implicit none

    type(c_ptr), value   :: isos              ! PDM_isosurface_t instance
    integer, intent(in)  :: id_isosurface     ! Iso-surface identifier
    integer, intent(in)  :: i_part            ! Partition identifier
    integer, intent(in)  :: connectivity_type ! Type of connectivity
    integer, pointer     :: connect_idx(:)    ! Index for connectivity
    integer, pointer     :: connect(:)        ! Connectivity
    integer, intent(in)  :: ownership         ! Ownership
    integer              :: n_entity          ! Number of entity

    integer(c_int) :: c_id_isosurface
    integer(c_int) :: c_i_part
    integer(c_int) :: c_connectivity_type
    type(c_ptr)    :: c_connect_idx
    type(c_ptr)    :: c_connect
    integer(c_int) :: c_ownership
    integer(c_int) :: c_n_entity
    
    interface
      function PDM_isosurface_pconnectivity_get_cf(isos,              &
                                                   id_isosurface,     &
                                                   i_part,            &
                                                   connectivity_type, &
                                                   connect_idx,       &
                                                   connect,           &
                                                   ownership)         &
                                                   result(n_entity)   &
      bind(c, name='PDM_isosurface_pconnectivity_get')

        use iso_c_binding
        implicit none

        type(c_ptr)    :: isos
        integer(c_int) :: id_isosurface
        integer(c_int) :: i_part
        integer(c_int) :: connectivity_type
        type(c_ptr)    :: connect_idx
        type(c_ptr)    :: connect
        integer(c_int) :: ownership
        integer(c_int) :: n_entity

      end function PDM_isosurface_pconnectivity_get_cf
    end interface

    c_id_isosurface     = id_isosurface
    c_i_part            = i_part
    c_connectivity_type = connectivity_type
    c_ownership         = ownership

    c_n_entity = PDM_isosurface_pconnectivity_get_cf(isos,                &
                                                     c_id_isosurface,     &
                                                     c_i_part,            &
                                                     c_connectivity_type, &
                                                     c_connect_idx,       &
                                                     c_connect,           &
                                                     c_ownership)

    n_entity = c_n_entity

    call c_f_pointer(c_connect_idx, &
                     connect_idx,   &
                     [n_entity+1])

    call c_f_pointer(c_connect, &
                     connect,   &
                     [connect_idx(n_entity+1)])    

  end function PDM_isosurface_pconnectivity_get


  function PDM_isosurface_pvtx_coord_get (isos,          &
                                          id_isosurface, &
                                          i_part,        &
                                          vtx_coord,     &
                                          ownership)     &
                                          result(n_vtx)
    ! Get coordinates of iso-surface vertices.
    use iso_c_binding
    implicit none

    type(c_ptr), value   :: isos            ! PDM_isosurface_t instance
    integer, intent(in)  :: id_isosurface   ! Iso-surface identifier
    integer, intent(in)  :: i_part          ! Partition identifier
    real(8), pointer     :: vtx_coord(:,:)  ! Vertex coordinates (shape = [3, n_vtx])
    integer, intent(in)  :: ownership       ! Ownership
    integer              :: n_vtx           ! Number of vertices

    integer(c_int) :: c_id_isosurface
    integer(c_int) :: c_i_part
    type(c_ptr)    :: c_vtx_coord
    integer(c_int) :: c_ownership
    integer(c_int) :: c_n_vtx
    
    interface
      function PDM_isosurface_pvtx_coord_get_cf(isos,          &
                                                id_isosurface, &
                                                i_part,        &
                                                vtx_coord,     &
                                                ownership)     &
                                                result(n_vtx)  &
      bind(c, name='PDM_isosurface_pvtx_coord_get')

        use iso_c_binding
        implicit none

        type(c_ptr)    :: isos
        integer(c_int) :: id_isosurface
        integer(c_int) :: i_part
        type(c_ptr)    :: vtx_coord
        integer(c_int) :: ownership
        integer(c_int) :: n_vtx

      end function PDM_isosurface_pvtx_coord_get_cf
    end interface

    c_id_isosurface = id_isosurface
    c_i_part        = i_part
    c_ownership     = ownership

    c_n_vtx = PDM_isosurface_pvtx_coord_get_cf(isos,            &
                                               c_id_isosurface, &
                                               c_i_part,        &
                                               c_vtx_coord,     &
                                               c_ownership)

    n_vtx = c_n_vtx

    call c_f_pointer(c_vtx_coord, &
                     vtx_coord,   &
                     [3, n_vtx])    

  end function PDM_isosurface_pvtx_coord_get


  function PDM_isosurface_ln_to_gn_get (isos,          &
                                        id_isosurface, &
                                        i_part,        &
                                        entity_type,   &
                                        ln_to_gn,      &
                                        ownership)     &
                                        result(n_entity)
    ! Get global ids of iso-surface entities.
    use iso_c_binding
    implicit none

    type(c_ptr), value                 :: isos            ! PDM_isosurface_t instance
    integer, intent(in)                :: id_isosurface   ! Iso-surface identifier
    integer, intent(in)                :: i_part          ! Partition identifier
    integer, intent(in)                :: entity_type     ! Type of mesh entity
    integer(kind=pdm_g_num_s), pointer :: ln_to_gn(:)     ! Global ids
    integer, intent(in)                :: ownership       ! Ownership
    integer                            :: n_entity        ! Number of vertices

    integer(c_int) :: c_id_isosurface
    integer(c_int) :: c_i_part
    integer(c_int) :: c_entity_type
    type(c_ptr)    :: c_ln_to_gn
    integer(c_int) :: c_ownership
    integer(c_int) :: c_n_entity
    
    interface
      function PDM_isosurface_ln_to_gn_get_cf(isos,            &
                                              id_isosurface,   &
                                              i_part,          &
                                              entity_type,     &
                                              ln_to_gn,        &
                                              ownership)       &
                                              result(n_entity) &
      bind(c, name='PDM_isosurface_ln_to_gn_get')

        use iso_c_binding
        implicit none

        type(c_ptr)    :: isos
        integer(c_int) :: id_isosurface
        integer(c_int) :: i_part
        integer(c_int) :: entity_type
        type(c_ptr)    :: ln_to_gn
        integer(c_int) :: ownership
        integer(c_int) :: n_entity

      end function PDM_isosurface_ln_to_gn_get_cf
    end interface

    c_id_isosurface = id_isosurface
    c_i_part        = i_part
    c_entity_type   = entity_type
    c_ownership     = ownership

    c_n_entity = PDM_isosurface_ln_to_gn_get_cf(isos,            &
                                                c_id_isosurface, &
                                                c_i_part,        &
                                                c_entity_type,   &
                                                c_ln_to_gn,      &
                                                c_ownership)

    n_entity = c_n_entity

    call c_f_pointer(c_ln_to_gn, &
                     ln_to_gn,   &
                     [n_entity])    

  end function PDM_isosurface_ln_to_gn_get


  function PDM_isosurface_pgroup_get (isos,                  &
                                      id_isosurface,         &
                                      i_part,                &
                                      entity_type,           &
                                      group_entity_idx,      &
                                      group_entity,          &
                                      group_entity_ln_to_gn, &
                                      ownership)             &
                                      result(n_group)
    ! Get group description
    use iso_c_binding
    implicit none

    type(c_ptr)                        :: isos                     ! PDM_isosurface_t instance
    integer, intent(in)                :: id_isosurface            ! Iso-surface identifier
    integer, intent(in)                :: i_part                   ! Partition identifier
    integer, intent(in)                :: entity_type              ! Type of mesh entity
    integer,                   pointer :: group_entity_idx(:)      ! Index for group->entity connectivity
    integer,                   pointer :: group_entity(:)          ! Group->entity connectivity
    integer(kind=pdm_g_num_s), pointer :: group_entity_ln_to_gn(:) ! Group->entity connectivity (group-specific global ids, size = group_entity_idx[n_group])
    integer, intent(in)                :: ownership                ! Ownership
    integer                            :: n_group                  ! Number of groups

    integer(c_int) :: c_id_isosurface
    integer(c_int) :: c_i_part
    integer(c_int) :: c_entity_type
    type(c_ptr)    :: c_group_entity_idx
    type(c_ptr)    :: c_group_entity
    type(c_ptr)    :: c_group_entity_ln_to_gn
    integer(c_int) :: c_ownership
    integer(c_int) :: c_n_group

    interface
      function PDM_isosurface_pgroup_get_cf (isos,                  &
                                             id_isosurface,         &
                                             i_part,                &
                                             entity_type,           &
                                             group_entity_idx,      &
                                             group_entity,          &
                                             group_entity_ln_to_gn, &
                                             ownership)             &
                                             result(n_group)        &
      bind(c, name='PDM_isosurface_pgroup_get')

        use iso_c_binding
        implicit none

        type(c_ptr)    :: isos
        integer(c_int) :: id_isosurface
        integer(c_int) :: i_part
        integer(c_int) :: entity_type
        type(c_ptr)    :: group_entity_idx
        type(c_ptr)    :: group_entity
        type(c_ptr)    :: group_entity_ln_to_gn
        integer(c_int) :: ownership
        integer(c_int) :: n_group

      end function PDM_isosurface_pgroup_get_cf
    end interface

    c_id_isosurface = id_isosurface
    c_i_part = i_part
    c_entity_type = entity_type
    c_ownership = ownership

    c_n_group = PDM_isosurface_pgroup_get_cf (isos,                    &
                                              c_id_isosurface,         &
                                              c_i_part,                &
                                              c_entity_type,           &
                                              c_group_entity_idx,      &
                                              c_group_entity,          &
                                              c_group_entity_ln_to_gn, &
                                              c_ownership)


    n_group = c_n_group

    call c_f_pointer(c_group_entity_idx, &
                     group_entity_idx,   &
                     [n_group+1])    

    call c_f_pointer(c_group_entity, &
                     group_entity,   &
                     [group_entity_idx(n_group+1)])    

    call c_f_pointer(c_group_entity_ln_to_gn, &
                     group_entity_ln_to_gn,   &
                     [group_entity_idx(n_group+1)])    


  end function PDM_isosurface_pgroup_get


  function PDM_isosurface_dconnectivity_get (isos,                &
                                             id_isosurface,       &
                                             connectivity_type,   &
                                             dconnect_idx,        &
                                             dconnect,            &
                                             ownership)           &
                                             result(n_entity)
    ! Get iso-surface block-distributed mesh connectivity.
    use iso_c_binding
    implicit none

    type(c_ptr)                        :: isos              ! PDM_isosurface_t instance
    integer, intent(in)                :: id_isosurface     ! Iso-surface identifier
    integer, intent(in)                :: connectivity_type ! Type of connectivity
    integer,                   pointer :: dconnect_idx(:)   ! Connectivity index
    integer(kind=pdm_g_num_s), pointer :: dconnect(:)       ! Connectivity
    integer, intent(in)                :: ownership         ! Ownership
    integer                            :: n_entity          ! Number of leading entity

    integer(c_int) :: c_id_isosurface
    integer(c_int) :: c_connectivity_type
    type(c_ptr)    :: c_dconnect_idx
    type(c_ptr)    :: c_dconnect
    integer(c_int) :: c_ownership
    integer(c_int) :: c_n_entity

    interface
      function PDM_isosurface_dconnectivity_get_cf (isos,                &
                                                    id_isosurface,       &
                                                    connectivity_type,   &
                                                    dconnect_idx,        &
                                                    dconnect,            &
                                                    ownership)           &
                                                    result(n_entity)     &
      bind(c, name='PDM_isosurface_dconnectivity_get')

        use iso_c_binding
        implicit none

        type(c_ptr)    :: isos
        integer(c_int) :: id_isosurface
        integer(c_int) :: connectivity_type
        type(c_ptr)    :: dconnect_idx
        type(c_ptr)    :: dconnect
        integer(c_int) :: ownership
        integer(c_int) :: n_entity

      end function PDM_isosurface_dconnectivity_get_cf
    end interface

    c_id_isosurface     = id_isosurface
    c_connectivity_type = connectivity_type
    c_ownership         = ownership

    c_n_entity = PDM_isosurface_dconnectivity_get_cf (isos,                &
                                                      c_id_isosurface,     &
                                                      c_connectivity_type, &
                                                      c_dconnect_idx,      &
                                                      c_dconnect,          &
                                                      c_ownership)           

    n_entity = c_n_entity

    call c_f_pointer(c_dconnect_idx, &
                     dconnect_idx,   &
                     [n_entity + 1])    

    call c_f_pointer(c_dconnect, &
                     dconnect,   &
                     [dconnect_idx(n_entity + 1)])

  end function PDM_isosurface_dconnectivity_get


  function PDM_isosurface_dparent_weight_get (isos,          &
                                             id_isosurface,  &
                                             entity_type,    &
                                             dparent_idx,    &
                                             dparent_weight, &
                                             ownership)      &
                                             result(n_iso_entity)
    ! Get iso-surface parent interpolation weight for iso entities.
    use iso_c_binding
    implicit none

    type(c_ptr)         :: isos              ! PDM_isosurface_t instance
    integer, intent(in) :: id_isosurface     ! Iso-surface identifier
    integer, intent(in) :: entity_type       ! Type of mesh entity
    integer, pointer    :: dparent_idx(:)    ! Index for parent weights
    real(8), pointer    :: dparent_weight(:) ! Parent weight
    integer, intent(in) :: ownership         ! Ownership
    integer             :: n_iso_entity      ! Number of iso-surface entities

    integer(c_int) :: c_id_isosurface
    integer(c_int) :: c_entity_type
    type(c_ptr)    :: c_dparent_idx
    type(c_ptr)    :: c_dparent_weight
    integer(c_int) :: c_ownership
    integer(c_int) :: c_n_iso_entity

    interface
      function PDM_isosurface_dparent_weight_get_cf (isos,                &
                                                     id_isosurface,       &
                                                     entity_type,         &
                                                     dparent_idx,         &
                                                     dparent_weight,      &
                                                     ownership)           &
                                                     result(n_iso_entity) &
      bind(c, name='PDM_isosurface_dparent_weight_get')

        use iso_c_binding
        implicit none

        type(c_ptr)    :: isos
        integer(c_int) :: id_isosurface
        integer(c_int) :: entity_type
        type(c_ptr)    :: dparent_idx
        type(c_ptr)    :: dparent_weight
        integer(c_int) :: ownership
        integer(c_int) :: n_iso_entity

      end function PDM_isosurface_dparent_weight_get_cf
    end interface

    c_id_isosurface = id_isosurface
    c_entity_type   = entity_type
    c_ownership     = ownership

    c_n_iso_entity = PDM_isosurface_dparent_weight_get_cf (isos,             &
                                                           c_id_isosurface,  &
                                                           c_entity_type,    &
                                                           c_dparent_idx,    &
                                                           c_dparent_weight, &
                                                           c_ownership)           

    n_iso_entity = c_n_iso_entity

    call c_f_pointer(c_dparent_idx, &
                     dparent_idx,   &
                     [n_iso_entity + 1])    

    call c_f_pointer(c_dparent_weight, &
                     dparent_weight,   &
                     [dparent_idx(n_iso_entity + 1)])

  end function PDM_isosurface_dparent_weight_get

  function PDM_isosurface_dvtx_coord_get (isos,          &
                                          id_isosurface, &
                                          dvtx_coord,    &
                                          ownership)     &
                                          result(dn_vtx)
    ! Get coordinates of block-distributed iso-surface vertices.
    use iso_c_binding
    implicit none

    type(c_ptr)         :: isos            ! PDM_isosurface_t instance
    integer, intent(in) :: id_isosurface   ! Iso-surface identifier
    real(8), pointer    :: dvtx_coord(:,:) ! Vertex coordinates (shape = [3, dn_vtx])
    integer, intent(in) :: ownership       ! Ownership
    integer             :: dn_vtx          ! Number of vertices

    integer(c_int) :: c_id_isosurface
    type(c_ptr)    :: c_dvtx_coord
    integer(c_int) :: c_ownership
    integer(c_int) :: c_dn_vtx

    interface
      function PDM_isosurface_dvtx_coord_get_cf (isos,          &
                                                 id_isosurface, &
                                                 dvtx_coord,    &
                                                 ownership)     &
                                                 result(dn_vtx) &
      bind(c, name='PDM_isosurface_dvtx_coord_get')

        use iso_c_binding
        implicit none

        type(c_ptr)    :: isos
        integer(c_int) :: id_isosurface
        type(c_ptr)    :: dvtx_coord
        integer(c_int) :: ownership
        integer(c_int) :: dn_vtx

      end function PDM_isosurface_dvtx_coord_get_cf
    end interface

    c_id_isosurface = id_isosurface
    c_ownership     = ownership

    c_dn_vtx = PDM_isosurface_dvtx_coord_get_cf (isos,            &
                                                 c_id_isosurface, &
                                                 c_dvtx_coord,    &
                                                 c_ownership)

    dn_vtx = c_dn_vtx

    call c_f_pointer(c_dvtx_coord, &
                     dvtx_coord,   &
                     [3, dn_vtx])

  end function PDM_isosurface_dvtx_coord_get


  function PDM_isosurface_distrib_get (isos,          &
                                       id_isosurface, &
                                       entity_type,   &
                                       distribution)  &
                                       result(n_entity)
    ! Get block distribution.
    use iso_c_binding
    implicit none

    type(c_ptr)                        :: isos            ! PDM_isosurface_t instance
    integer, intent(in)                :: id_isosurface   ! Iso-surface identifier
    integer, intent(in)                :: entity_type     ! Type of mesh entity
    integer(kind=pdm_g_num_s), pointer :: distribution(:) ! Entity distribution
    integer                            :: n_entity        ! Number of entity

    integer(c_int) :: c_id_isosurface
    integer(c_int) :: c_entity_type
    type(c_ptr)    :: c_distribution
    integer(c_int) :: c_n_entity

    interface
      function PDM_isosurface_distrib_get_cf (isos,            &
                                              id_isosurface,   &
                                              entity_type,     &
                                              distribution)    &
                                              result(n_entity) &
      bind(c, name='PDM_isosurface_distrib_get')

        use iso_c_binding
        implicit none

        type(c_ptr)    :: isos
        integer(c_int) :: id_isosurface
        integer(c_int) :: entity_type
        type(c_ptr)    :: distribution
        integer(c_int) :: n_entity

      end function
    end interface

    c_id_isosurface = id_isosurface
    c_entity_type   = entity_type

    c_n_entity = PDM_isosurface_distrib_get_cf (isos,            &
                                                c_id_isosurface, &
                                                c_entity_type,   &
                                                c_distribution)

    n_entity = c_n_entity

    call c_f_pointer(c_distribution, &
                     distribution,   &
                     [n_entity])

  end function PDM_isosurface_distrib_get


  function PDM_isosurface_dgroup_get (isos,              &
                                      id_isosurface,     &
                                      entity_type,       &
                                      dgroup_entity_idx, &
                                      dgroup_entity,     &
                                      ownership)         &
                                      result(n_group)
    ! Get block-distributed group description.
    use iso_c_binding
    implicit none

    type(c_ptr)                        :: isos                 ! PDM_isosurface_t instance
    integer, intent(in)                :: id_isosurface        ! Iso-surface identifier
    integer, intent(in)                :: entity_type          ! Type of mesh entity
    integer,                   pointer :: dgroup_entity_idx(:) ! Index for group→entity connectivity (size = n_group+1)
    integer(kind=pdm_g_num_s), pointer :: dgroup_entity(:)     ! Group->entity connectivity (group-specific global ids, size = group_entity_idx[n_group+1])
    integer, intent(in)                :: ownership            ! Ownership
    integer                            :: n_group              ! Number of groups

    integer(c_int) :: c_id_isosurface
    integer(c_int) :: c_entity_type
    type(c_ptr)    :: c_dgroup_entity_idx
    type(c_ptr)    :: c_dgroup_entity
    integer(c_int) :: c_ownership
    integer(c_int) :: c_n_group

    interface
      function PDM_isosurface_dgroup_get_cf (isos,              &
                                             id_isosurface,     &
                                             entity_type,       &
                                             dgroup_entity_idx, &
                                             dgroup_entity,     &
                                             ownership)         &
                                             result(n_group)    &
      bind(c, name='PDM_isosurface_dgroup_get')

        use iso_c_binding
        implicit none

        type(c_ptr)    :: isos
        integer(c_int) :: id_isosurface
        integer(c_int) :: entity_type
        type(c_ptr)    :: dgroup_entity_idx
        type(c_ptr)    :: dgroup_entity
        integer(c_int) :: ownership
        integer(c_int) :: n_group

      end function PDM_isosurface_dgroup_get_cf
    end interface

    c_id_isosurface = id_isosurface
    c_entity_type   = entity_type
    c_ownership     = ownership

    c_n_group = PDM_isosurface_dgroup_get_cf (isos,                &
                                              c_id_isosurface,     &
                                              c_entity_type,       &
                                              c_dgroup_entity_idx, &
                                              c_dgroup_entity,     &
                                              c_ownership)

    n_group = c_n_group

    call c_f_pointer(c_dgroup_entity_idx, &
                     dgroup_entity_idx,   &
                     [n_group+1])

    call c_f_pointer(c_dgroup_entity, &
                     dgroup_entity,   &
                     [dgroup_entity_idx(n_group+1)])

  end function PDM_isosurface_dgroup_get


  function PDM_isosurface_pisovalue_entity_idx_get (isos,                &
                                                    id_isosurface,       &
                                                    i_part,              &
                                                    entity_type,         &
                                                    isovalue_entity_idx, &
                                                    ownership)           &
                                                    result(n_isovalues)
    ! Get isovalue.
    use iso_c_binding
    implicit none

    type(c_ptr)         :: isos                    ! PDM_isosurface_t instance
    integer, intent(in) :: id_isosurface           ! Iso-surface identifier
    integer, intent(in) :: i_part                  ! Partition identifier
    integer, intent(in) :: entity_type             ! Type of mesh entity
    integer, pointer    :: isovalue_entity_idx (:) ! Index for isovalue->entity connectivity (size = n_isovalue+1)
    integer, intent(in) :: ownership               ! Ownership
    integer             :: n_isovalues             ! Number of isovalues

    integer(c_int) :: c_id_isosurface
    integer(c_int) :: c_i_part
    integer(c_int) :: c_entity_type
    type(c_ptr)    :: c_isovalue_entity_idx
    integer(c_int) :: c_ownership
    integer(c_int) :: c_n_isovalues

    interface
      function PDM_isosurface_pisovalue_entity_idx_get_cf (isos,                &
                                                           id_isosurface,       &
                                                           i_part,              &
                                                           entity_type,         &
                                                           isovalue_entity_idx, &
                                                           ownership)           &
                                                           result(n_isovalues)  &
      bind(c, name='PDM_isosurface_pisovalue_entity_idx_get')

        use iso_c_binding
        implicit none

        type(c_ptr) :: isos
        integer(c_int) :: id_isosurface
        integer(c_int) :: i_part
        integer(c_int) :: entity_type
        type(c_ptr) :: isovalue_entity_idx
        integer(c_int) :: ownership
        integer(c_int) :: n_isovalues

      end function PDM_isosurface_pisovalue_entity_idx_get_cf
    end interface

    c_id_isosurface = id_isosurface
    c_i_part = i_part
    c_entity_type = entity_type
    c_ownership = ownership

    c_n_isovalues = PDM_isosurface_pisovalue_entity_idx_get_cf (isos,                  &
                                                                c_id_isosurface,       &
                                                                c_i_part,              &
                                                                c_entity_type,         &
                                                                c_isovalue_entity_idx, &
                                                                c_ownership)

    n_isovalues = c_n_isovalues

    call c_f_pointer(c_isovalue_entity_idx, &
                     isovalue_entity_idx,   &
                     [n_isovalues+1])

  end function PDM_isosurface_pisovalue_entity_idx_get


  function PDM_isosurface_disovalue_entity_get (isos,                 &
                                                id_isosurface,        &
                                                entity_type,          &
                                                disovalue_entity_idx, &
                                                disovalue_entity,     &
                                                ownership)            &
                                                result(n_isovalues)
    ! Get distributed isovalue→entity.
    use iso_c_binding
    implicit none

    type(c_ptr)                        :: isos                    ! PDM_isosurface_t instance
    integer, intent(in)                :: id_isosurface           ! Iso-surface identifier
    integer, intent(in)                :: entity_type             ! Type of mesh entity
    integer,                   pointer :: disovalue_entity_idx(:) ! Index for isovalue->entity connectivity (size = [n_isovalue+1])
    integer(kind=pdm_g_num_s), pointer :: disovalue_entity(:)     ! Isovalue→entity connectivity (size = disovalue_entity_idx[n_isovalue+1])
    integer, intent(in)                :: ownership               ! Ownership
    integer                            :: n_isovalues             ! Number of isovalues

    integer(c_int) :: c_id_isosurface
    integer(c_int) :: c_entity_type
    type(c_ptr)    :: c_disovalue_entity_idx
    type(c_ptr)    :: c_disovalue_entity
    integer(c_int) :: c_ownership
    integer(c_int) :: c_n_isovalues

    interface
      function PDM_isosurface_disovalue_entity_get_cf (isos,                 &
                                                       id_isosurface,        &
                                                       entity_type,          &
                                                       disovalue_entity_idx, &
                                                       disovalue_entity,     &
                                                       ownership)            &
                                                       result(n_isovalues)   &
      bind(c, name='PDM_isosurface_disovalue_entity_get')

        use iso_c_binding
        implicit none

        type(c_ptr)    :: isos
        integer(c_int) :: id_isosurface
        integer(c_int) :: entity_type
        type(c_ptr)    :: disovalue_entity_idx
        type(c_ptr)    :: disovalue_entity
        integer(c_int) :: ownership
        integer(c_int) :: n_isovalues

      end function PDM_isosurface_disovalue_entity_get_cf
    end interface

    c_id_isosurface = id_isosurface
    c_entity_type   = entity_type
    c_ownership     = ownership

    c_n_isovalues = PDM_isosurface_disovalue_entity_get_cf (isos,                   &
                                                            c_id_isosurface,        &
                                                            c_entity_type,          &
                                                            c_disovalue_entity_idx, &
                                                            c_disovalue_entity,     &
                                                            c_ownership)

    n_isovalues = c_n_isovalues

    call c_f_pointer(c_disovalue_entity_idx, &
                     disovalue_entity_idx,   &
                     [n_isovalues+1])

    call c_f_pointer(c_disovalue_entity, &
                     disovalue_entity,   &
                     [disovalue_entity_idx(n_isovalues+1)])

  end function PDM_isosurface_disovalue_entity_get


  function PDM_isosurface_plocal_parent_get (isos,              &
                                             id_isosurface,     &
                                             i_part,            &
                                             entity_type,       &
                                             entity_parent_idx, &
                                             entity_parent,     &
                                             ownership)         &
                                             result(n_entity)
    ! Get local parents of iso-surface entities.
    use iso_c_binding
    implicit none

    type(c_ptr)         :: isos                 ! PDM_isosurface_t instance
    integer, intent(in) :: id_isosurface        ! Iso-surface identifier
    integer, intent(in) :: i_part               ! Partition identifier
    integer, intent(in) :: entity_type          ! Type of mesh entity
    integer, pointer    :: entity_parent_idx(:) ! Index for isosurface entity->parent connectivity (size = \p n_entity + 1)
    integer, pointer    :: entity_parent(:)     ! Isosurface entity->parent connectivity (size = \p entity_parent_idx[\p n_entity])
    integer, intent(in) :: ownership            ! Ownership
    integer             :: n_entity             ! Number of entities

    integer(c_int) :: c_id_isosurface
    integer(c_int) :: c_i_part
    integer(c_int) :: c_entity_type
    type(c_ptr)    :: c_entity_parent_idx
    type(c_ptr)    :: c_entity_parent
    integer(c_int) :: c_ownership
    integer(c_int) :: c_n_entity

    interface
      function PDM_isosurface_plocal_parent_get_cf (isos,              &
                                                    id_isosurface,     &
                                                    i_part,            &
                                                    entity_type,       &
                                                    entity_parent_idx, &
                                                    entity_parent,     &
                                                    ownership)         &
                                                    result(n_entity)   &
      bind(c, name='PDM_isosurface_plocal_parent_get')

        use iso_c_binding
        implicit none

        type(c_ptr)    :: isos
        integer(c_int) :: id_isosurface
        integer(c_int) :: i_part
        integer(c_int) :: entity_type
        type(c_ptr)    :: entity_parent_idx
        type(c_ptr)    :: entity_parent
        integer(c_int) :: ownership
        integer(c_int) :: n_entity

      end function PDM_isosurface_plocal_parent_get_cf
    end interface

    c_id_isosurface = id_isosurface
    c_i_part        = i_part
    c_entity_type   = entity_type
    c_ownership     = ownership

    c_n_entity = PDM_isosurface_plocal_parent_get_cf (isos,                &
                                                      c_id_isosurface,     &
                                                      c_i_part,            &
                                                      c_entity_type,       &
                                                      c_entity_parent_idx, &
                                                      c_entity_parent,     &
                                                      c_ownership)

    n_entity = c_n_entity

    call c_f_pointer(c_entity_parent_idx, &
                     entity_parent_idx,   &
                     [n_entity+1])

    call c_f_pointer(c_entity_parent, &
                     entity_parent,   &
                     [entity_parent_idx(n_entity+1)])

  end function PDM_isosurface_plocal_parent_get


  function PDM_isosurface_pparent_weight_get (isos, &
                                              id_isosurface, &
                                              i_part, &
                                              entity_type, &
                                              parent_idx, &
                                              parent_weight, &
                                              ownership) &
                                              result(n_iso_entity)
    ! Get interpolation weights of iso-surface entities.
    use iso_c_binding
    implicit none

    type(c_ptr)         :: isos             ! PDM_isosurface_t instance
    integer, intent(in) :: id_isosurface    ! Iso-surface identifier
    integer, intent(in) :: i_part           ! Partition identifier
    integer, intent(in) :: entity_type      ! Type of mesh entity
    integer, pointer    :: parent_idx(:)    ! Index for interpolation weights
    real(8), pointer    :: parent_weight(:) ! Interpolation weights
    integer, intent(in) :: ownership        ! Ownership
    integer             :: n_iso_entity     ! Number of iso-surface entity

    integer(c_int) :: c_id_isosurface
    integer(c_int) :: c_i_part
    integer(c_int) :: c_entity_type
    type(c_ptr)    :: c_parent_idx
    type(c_ptr)    :: c_parent_weight
    integer(c_int) :: c_ownership
    integer(c_int) :: c_n_iso_entity


    interface
      function PDM_isosurface_pparent_weight_get_cf (isos, &
                                                     id_isosurface, &
                                                     i_part, &
                                                     entity_type, &
                                                     parent_idx, &
                                                     parent_weight, &
                                                     ownership) &
                                                     result(n_iso_entity) &
      bind(c, name='PDM_isosurface_pparent_weight_get')

        use iso_c_binding
        implicit none

        type(c_ptr)    :: isos
        integer(c_int) :: id_isosurface
        integer(c_int) :: i_part
        integer(c_int) :: entity_type
        type(c_ptr)    :: parent_idx
        type(c_ptr)    :: parent_weight
        integer(c_int) :: ownership
        integer(c_int) :: n_iso_entity

      end function PDM_isosurface_pparent_weight_get_cf
    end interface

    c_id_isosurface = id_isosurface
    c_i_part        = i_part
    c_entity_type   = entity_type
    c_ownership     = ownership

    c_n_iso_entity = PDM_isosurface_pparent_weight_get_cf (isos, &
                                                           c_id_isosurface, &
                                                           c_i_part, &
                                                           c_entity_type, &
                                                           c_parent_idx, &
                                                           c_parent_weight, &
                                                           c_ownership)

    n_iso_entity = c_n_iso_entity

    call c_f_pointer(c_parent_idx, &
                     parent_idx,   &
                     [n_iso_entity+1])

    call c_f_pointer(c_parent_weight, &
                     parent_weight,   &
                     [parent_idx(n_iso_entity+1)])

  end function PDM_isosurface_pparent_weight_get




end module pdm_isosurface
