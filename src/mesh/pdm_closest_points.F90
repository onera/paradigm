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

module pdm_closest_points

  use iso_c_binding
  use pdm

  implicit none


  interface

  !>
  !!
  !! \brief  Get the number of target points in a partition
  !!
  !! \param [in]  cls     Pointer to \ref PDM_closest_points_t object
  !! \param [in]  i_part  Index of partition of the target cloud
  !!
  !! \return   Number of target point in the partition \ref i_part
  !!
  !!

  function PDM_closest_points_n_tgt_get (cls,     &
                                        i_part)  &
    result (n_tgt)                                &

    bind (c, name='PDM_closest_points_n_tgt_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value :: cls
    integer(c_int), value :: i_part
    integer(c_int)        :: n_tgt

  end function PDM_closest_points_n_tgt_get


  !>
  !!
  !! \brief  Get the number of source points in a partition
  !!
  !! \param [in]  cls     Pointer to \ref PDM_closest_points_t object
  !! \param [in]  i_part  Index of partition of the target cloud
  !!
  !! \return   Number of source point in the partition \ref i_part
  !!
  !!

  function PDM_closest_points_n_src_get (cls,     &
                                        i_part)  &
    result (n_src)                                &

    bind (c, name='PDM_closest_points_n_src_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value :: cls
    integer(c_int), value :: i_part
    integer(c_int)        :: n_src

  end function PDM_closest_points_n_src_get


  !>
  !!
  !! \brief  Get the number of closest points
  !!
  !! \param [in]  cls     Pointer to \ref PDM_closest_points_t object
  !!
  !! \return   Number of closest points
  !!
  !!

  function PDM_closest_points_n_closest_get (cls)  &
    result (n_closest)                             &

    bind (c, name='PDM_closest_points_n_closest_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value :: cls
    integer(c_int)        :: n_closest

  end function PDM_closest_points_n_closest_get


  !>
  !!
  !! \brief Disable reverse results computation
  !!
  !! \param [inout] cls        Pointer to \ref PDM_closest_point_t object
  !!
  !!

  subroutine PDM_closest_points_reverse_results_disable (cls) &
    bind (c, name = 'PDM_closest_points_reverse_results_disable')

    use iso_c_binding
    implicit none

    type(c_ptr), value :: cls

  end subroutine PDM_closest_points_reverse_results_disable


  end interface


  contains


  subroutine PDM_closest_points_create(cls,       &
                                       f_comm,    &
                                       n_closest, &
                                       owner)
    ! Create a structure to look for the closest points of a point cloud (target cloud) in an other point cloud (source cloud)
    implicit none

    type(c_ptr), intent(out) :: cls       ! PDM_closest_points instance
    integer,     intent(in)  :: f_comm    ! MPI communicator
    integer,     intent(in)  :: n_closest ! Number of closest source points to find for each target point
    integer,     intent(in)  :: owner     ! Ownership

    type(c_ptr)           :: c_comm

    interface
      function PDM_closest_points_create_cf(comm,      &
                                            n_closest, &
                                            owner)     &
      result(cls)                                      &
      bind (c, name = 'PDM_closest_points_create')
        use iso_c_binding
        implicit none
        type(c_ptr)           :: cls
        type(c_ptr), value    :: comm
        integer(c_int), value :: n_closest
        integer(c_int), value :: owner
      end function PDM_closest_points_create_cf
    end interface

    c_comm = PDM_MPI_Comm_f2c(f_comm)

    cls = PDM_closest_points_create_cf(c_comm,    &
                                        n_closest, &
                                        owner)

  end subroutine PDM_closest_points_create



  subroutine PDM_closest_points_n_part_cloud_set(cls, &
                                                  n_part_cloud_src, &
                                                  n_part_cloud_tgt)
    ! Set the number of partitions of both point clouds
    implicit none

    type(c_ptr), intent(in) :: cls              ! PDM_closest_points instance
    integer,     intent(in) :: n_part_cloud_src ! Number of partitions in the source cloud
    integer,     intent(in) :: n_part_cloud_tgt ! Number of partitions in the target cloud

    interface
      subroutine PDM_closest_points_n_part_cloud_set_cf(cls, &
                                                        n_part_cloud_src, &
                                                        n_part_cloud_tgt) &
      bind (c, name = 'PDM_closest_points_n_part_cloud_set')
        use iso_c_binding
        implicit none
        type(c_ptr), value :: cls
        integer(c_int), value :: n_part_cloud_src
        integer(c_int), value :: n_part_cloud_tgt
      end subroutine PDM_closest_points_n_part_cloud_set_cf
    end interface

    call PDM_closest_points_n_part_cloud_set_cf(cls, &
                                                n_part_cloud_src, &
                                                n_part_cloud_tgt)

  end subroutine PDM_closest_points_n_part_cloud_set



  subroutine PDM_closest_points_tgt_cloud_set(cls,      &
                                              i_part,   &
                                              n_points, &
                                              coords,   &
                                              gnum)
    ! Set the target point cloud
    implicit none

    type(c_ptr),            intent(in) :: cls         ! PDM_closest_points instance
    integer,                intent(in) :: i_part      ! Partition identifier
    integer,                intent(in) :: n_points    ! Number of target points
    real(8),                   pointer :: coords(:,:) ! Coordinates of target points (shape = [3, ``n_points``])
    integer(kind=pdm_g_num_s), pointer :: gnum(:)     ! Global IDs of target points

    type(c_ptr)                        :: c_coords
    type(c_ptr)                        :: c_gnum

    interface
      subroutine PDM_closest_points_tgt_cloud_set_cf(cls, &
                                                      i_part, &
                                                      n_points, &
                                                      coords, &
                                                      gnum) &
        bind (c, name = 'PDM_closest_points_tgt_cloud_set')
        use iso_c_binding
        implicit none
        type(c_ptr),    value :: cls
        integer(c_int), value :: i_part
        integer(c_int), value :: n_points
        type(c_ptr),    value :: coords
        type(c_ptr),    value :: gnum
      end subroutine PDM_closest_points_tgt_cloud_set_cf
    end interface

    c_coords = C_NULL_PTR
    if (associated(coords)) then
      c_coords = c_loc(coords)
    end if

    c_gnum = C_NULL_PTR
    if (associated(coords)) then
      c_gnum = c_loc(gnum)
    end if

    call PDM_closest_points_tgt_cloud_set_cf(cls,      &
                                              i_part,   &
                                              n_points, &
                                              c_coords, &
                                              c_gnum)

  end subroutine PDM_closest_points_tgt_cloud_set



  subroutine PDM_closest_points_src_cloud_set(cls,      &
                                              i_part,   &
                                              n_points, &
                                              coords,   &
                                              gnum)
    ! Set the source point cloud
    implicit none

    type(c_ptr),            intent(in) :: cls         ! PDM_closest_points instance
    integer,                intent(in) :: i_part      ! Partition identifier
    integer,                intent(in) :: n_points    ! Number of source points
    real(8),                   pointer :: coords(:,:) ! Coordinates of source points (shape = [3, ``n_points``])
    integer(kind=pdm_g_num_s), pointer :: gnum(:)     ! Global IDs of source points

    type(c_ptr)                        :: c_coords
    type(c_ptr)                        :: c_gnum

    interface
      subroutine PDM_closest_points_src_cloud_set_cf(cls, &
                                                      i_part, &
                                                      n_points, &
                                                      coords, &
                                                      gnum) &
        bind (c, name = 'PDM_closest_points_src_cloud_set')
        use iso_c_binding
        implicit none
        type(c_ptr),    value :: cls
        integer(c_int), value :: i_part
        integer(c_int), value :: n_points
        type(c_ptr),    value :: coords
        type(c_ptr),    value :: gnum
      end subroutine PDM_closest_points_src_cloud_set_cf
    end interface

    c_coords = C_NULL_PTR
    if (associated(coords)) then
      c_coords = c_loc(coords)
    end if

    c_gnum = C_NULL_PTR
    if (associated(coords)) then
      c_gnum = c_loc(gnum)
    end if

    call PDM_closest_points_src_cloud_set_cf(cls,      &
                                              i_part,   &
                                              n_points, &
                                              c_coords, &
                                              c_gnum)

  end subroutine PDM_closest_points_src_cloud_set

  

  subroutine PDM_closest_points_compute(cls)
    ! Look for closest points
    implicit none

    type(c_ptr), intent(in) :: cls ! PDM_closest_points instance

    interface
      subroutine PDM_closest_points_compute_cf(cls) &
        bind (c, name = 'PDM_closest_points_compute')
        use iso_c_binding
        implicit none
        type(c_ptr), value :: cls
      end subroutine PDM_closest_points_compute_cf
    end interface

    call PDM_closest_points_compute_cf(cls)

  end subroutine PDM_closest_points_compute



  subroutine PDM_closest_points_dump_times(cls)
    ! Dump elapsed and CPU times
    implicit none

    type(c_ptr), intent(in) :: cls ! PDM_closest_points instance

    interface
      subroutine PDM_closest_points_dump_times_cf(cls) &
        bind (c, name = 'PDM_closest_points_dump_times')
        use iso_c_binding
        implicit none
        type(c_ptr), value :: cls
      end subroutine PDM_closest_points_dump_times_cf
    end interface

    call PDM_closest_points_dump_times_cf(cls)

  end subroutine PDM_closest_points_dump_times



  subroutine PDM_closest_points_get(cls, &
                                    i_part_tgt, &
                                    closest_src_gnum, &
                                    closest_src_distance)

    ! Get closest source points global IDs and distance
    implicit none

    type(c_ptr),               intent(in) :: cls                     ! PDM_closest_points instance
    integer,                   intent(in) :: i_part_tgt              ! Partition identifier
    integer(kind=pdm_g_num_s), pointer    :: closest_src_gnum(:)     ! Global IDs of closest source points 
    real(8),                   pointer    :: closest_src_distance(:) ! Squared distance of closest source points 

    type(c_ptr)                           :: c_closest_src_gnum
    type(c_ptr)                           :: c_closest_src_distance
    integer(c_int)                        :: n_tgt, n_closest

    interface
      subroutine PDM_closest_points_get_cf(cls, &
                                           i_part_tgt, &
                                           closest_src_gnum, &
                                           closest_src_distance) &
        bind (c, name = 'PDM_closest_points_get')
        use iso_c_binding
        implicit none
        type(c_ptr),    value :: cls
        integer(c_int), value :: i_part_tgt
        type(c_ptr)           :: closest_src_gnum
        type(c_ptr)           :: closest_src_distance
      end subroutine PDM_closest_points_get_cf
    end interface

    c_closest_src_gnum     = C_NULL_PTR
    c_closest_src_distance = C_NULL_PTR

    call PDM_closest_points_get_cf(cls,                    &
                                   i_part_tgt,             &
                                   c_closest_src_gnum,     &
                                   c_closest_src_distance)

    n_tgt = pdm_closest_points_n_tgt_get(cls,        &
                                         i_part_tgt)

    n_closest = pdm_closest_points_n_closest_get(cls)

    call c_f_pointer(c_closest_src_gnum, &
                     closest_src_gnum,   &
                     [n_tgt*n_closest])

    call c_f_pointer(c_closest_src_distance, &
                     closest_src_distance,   &
                     [n_tgt*n_closest])

  end subroutine PDM_closest_points_get



  subroutine PDM_closest_points_tgt_in_src_get(cls, &
                                               i_part_src, &
                                               tgt_in_src_idx, &
                                               tgt_in_src)

    ! Get source->target correspondence 
    implicit none

    type(c_ptr),            intent(in) :: cls                ! PDM_closest_points instance
    integer,                intent(in) :: i_part_src         ! Partition identifier
    integer(kind=pdm_l_num_s), pointer :: tgt_in_src_idx(:)  ! Index for source->target correspondence
    integer(kind=pdm_g_num_s), pointer :: tgt_in_src(:)      ! Source->target correspondence (global IDs)

    type(c_ptr)                        :: c_tgt_in_src_idx
    type(c_ptr)                        :: c_tgt_in_src
    integer(c_int)                     :: n_src

    interface
      subroutine PDM_closest_points_tgt_in_src_get_cf(cls, &
                                                      i_part_src, &
                                                      tgt_in_src_idx, &
                                                      tgt_in_src) &
        bind (c, name = 'PDM_closest_points_tgt_in_src_get')
        use iso_c_binding
        implicit none
        type(c_ptr),    value :: cls
        integer(c_int), value :: i_part_src
        type(c_ptr)           :: tgt_in_src_idx
        type(c_ptr)           :: tgt_in_src
      end subroutine PDM_closest_points_tgt_in_src_get_cf
    end interface

    c_tgt_in_src_idx = C_NULL_PTR
    c_tgt_in_src     = C_NULL_PTR
    call PDM_closest_points_tgt_in_src_get_cf(cls,                  &
                                              i_part_src,           &
                                              c_tgt_in_src_idx,     &
                                              c_tgt_in_src)

    n_src = pdm_closest_points_n_src_get(cls,        &
                                         i_part_src)

    call c_f_pointer(c_tgt_in_src_idx, &
                     tgt_in_src_idx,   &
                     [n_src+1])

    call c_f_pointer(c_tgt_in_src, &
                     tgt_in_src,   &
                     [tgt_in_src_idx(n_src+1)])

  end subroutine PDM_closest_points_tgt_in_src_get



  subroutine PDM_closest_points_tgt_in_src_dist_get(cls, &
                                                    i_part_src, &
                                                    tgt_in_src_idx, &
                                                    tgt_in_src_dist)

    ! Get source->target distance 
    implicit none

    type(c_ptr),            intent(in) :: cls                ! PDM_closest_points instance
    integer,                intent(in) :: i_part_src         ! Partition identifier
    integer(kind=pdm_l_num_s), pointer :: tgt_in_src_idx(:)  ! Index for source->target correspondence
    real(8),                   pointer :: tgt_in_src_dist(:) ! Source->target squared distance

    type(c_ptr)                        :: c_tgt_in_src_idx
    type(c_ptr)                        :: c_tgt_in_src_dist
    integer(c_int)                     :: n_src

    interface
      subroutine PDM_closest_points_tgt_in_src_dist_get_cf(cls, &
                                                           i_part_src, &
                                                           tgt_in_src_idx, &
                                                           tgt_in_src_dist) &
        bind (c, name = 'PDM_closest_points_tgt_in_src_dist_get')
        use iso_c_binding
        implicit none
        type(c_ptr),    value :: cls
        integer(c_int), value :: i_part_src
        type(c_ptr)           :: tgt_in_src_idx
        type(c_ptr)           :: tgt_in_src_dist
      end subroutine PDM_closest_points_tgt_in_src_dist_get_cf
    end interface

    c_tgt_in_src_idx  = C_NULL_PTR
    c_tgt_in_src_dist = C_NULL_PTR
    call PDM_closest_points_tgt_in_src_dist_get_cf(cls,                  &
                                                   i_part_src,           &
                                                   c_tgt_in_src_idx,     &
                                                   c_tgt_in_src_dist)

    n_src = pdm_closest_points_n_src_get(cls,        &
                                         i_part_src)

    call c_f_pointer(c_tgt_in_src_idx, &
                     tgt_in_src_idx,   &
                     [n_src+1])

    call c_f_pointer(c_tgt_in_src_dist, &
                     tgt_in_src_dist,   &
                     [tgt_in_src_idx(n_src+1)])

  end subroutine PDM_closest_points_tgt_in_src_dist_get



  subroutine PDM_closest_points_part_to_part_get(cls,   &
                                                 ptp,   &
                                                 owner)
    ! Get part_to_part object to exchange data between the source and target point clouds
    implicit none

    type (c_ptr),   value :: cls    ! PDM_closest_points instance
    type (c_ptr)          :: ptp    ! Pointer to PDM_part_to_part object
    integer(c_int), value :: owner  ! Ownership for ``ptp``

    interface
      subroutine PDM_closest_points_part_to_part_get_cf(cls,   &
                                                        ptp,   &
                                                        owner) &
      bind (c, name = 'PDM_closest_points_part_to_part_get')
        use iso_c_binding
        implicit none
        type (c_ptr),   value :: cls
        type (c_ptr)          :: ptp
        integer(c_int), value :: owner
      end subroutine PDM_closest_points_part_to_part_get_cf
    end interface

    call PDM_closest_points_part_to_part_get_cf(cls,   &
                                                ptp,   &
                                                owner)

  end subroutine PDM_closest_points_part_to_part_get
  
  
  
  subroutine PDM_closest_points_free(cls)
    ! Free a PDM_closest_points instance
    implicit none

    type(c_ptr), intent(inout) :: cls

    interface
      subroutine PDM_closest_points_free_cf(cls) &
        bind (c, name = 'PDM_closest_points_free')
        use iso_c_binding
        implicit none
        type(c_ptr), value :: cls
      end subroutine PDM_closest_points_free_cf
    end interface

    call PDM_closest_points_free_cf(cls)

  end subroutine PDM_closest_points_free

end module pdm_closest_points
