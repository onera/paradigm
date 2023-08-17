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

module pdm_part_connectivity_transform

  use pdm
  implicit none

  interface PDM_combine_connectivity ; module procedure &
  PDM_combine_connectivity_
  end interface

  interface PDM_connectivity_transpose ; module procedure &
  PDM_connectivity_transpose_
  end interface

  private :: PDM_combine_connectivity_
  private :: PDM_connectivity_transpose_

  interface

    !!>
    !!
    !! \brief PDM_combine_connectivity

    subroutine PDM_combine_connectivity_cf (n_entity1,           &
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

      type(c_ptr), value    :: entity1_entity2_idx
      type(c_ptr), value    :: entity1_entity2
      type(c_ptr), value    :: entity2_entity3_idx
      type(c_ptr), value    :: entity2_entity3

      type(c_ptr)           :: entity1_entity3_idx
      type(c_ptr)           :: entity1_entity3

    end subroutine PDM_combine_connectivity_cf

    !!>
    !!
    !! \brief PDM_connectivity_transpose

    subroutine PDM_connectivity_transpose_cf (n_entity1,           &
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

      type(c_ptr), value    :: entity1_entity2_idx
      type(c_ptr), value    :: entity1_entity2

      type(c_ptr)           :: entity2_entity1_idx
      type(c_ptr)           :: entity2_entity1

    end subroutine PDM_connectivity_transpose_cf

  end interface

  contains

    !!>
    !!
    !! \brief Combine entity1->entity2 with entity2->entity3 to create entity1->entity3
    !!
    !! \param [in]   n_entity1
    !! \param [in]   entity1_entity2_idx
    !! \param [in]   entity1_entity2
    !! \param [in]   entity2_entity3_idx
    !! \param [in]   entity2_entity3
    !! \param [out]  entity1_entity3_idx
    !! \param [out]  entity1_entity3
    !!

    subroutine PDM_combine_connectivity_ (n_entity1,           &
                                          entity1_entity2_idx, &
                                          entity1_entity2,     &
                                          entity2_entity3_idx, &
                                          entity2_entity3,     &
                                          entity1_entity3_idx, &
                                          entity1_entity3)

      use iso_c_binding

      implicit none

      integer, intent(in) :: n_entity1

      integer, pointer    :: entity1_entity2_idx(:)
      integer, pointer    :: entity1_entity2(:)
      integer, pointer    :: entity2_entity3_idx(:)
      integer, pointer    :: entity2_entity3(:)
      integer, pointer    :: entity1_entity3_idx(:)
      integer, pointer    :: entity1_entity3(:)

      integer(c_int) :: c_n_entity1

      type(c_ptr)    :: c_entity1_entity2_idx     = C_NULL_PTR
      type(c_ptr)    :: c_entity1_entity2         = C_NULL_PTR
      type(c_ptr)    :: c_entity2_entity3_idx     = C_NULL_PTR
      type(c_ptr)    :: c_entity2_entity3         = C_NULL_PTR
      type(c_ptr)    :: c_entity1_entity3_idx     = C_NULL_PTR
      type(c_ptr)    :: c_entity1_entity3         = C_NULL_PTR

      c_n_entity1 = n_entity1

      c_entity1_entity2_idx = c_loc(entity1_entity2_idx)
      c_entity1_entity2     = c_loc(entity1_entity2)
      c_entity2_entity3_idx = c_loc(entity2_entity3_idx)
      c_entity2_entity3     = c_loc(entity2_entity3)

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

    end subroutine PDM_combine_connectivity_

    !!>
    !!
    !! \brief Combine entity1->entity2 with entity2->entity3 to create entity1->entity3
    !!
    !! \param [in]   n_entity1
    !! \param [in]   entity1_entity2_idx
    !! \param [in]   entity1_entity2
    !! \param [in]   entity2_entity3_idx
    !! \param [in]   entity2_entity3
    !! \param [out]  entity1_entity3_idx
    !! \param [out]  entity1_entity3
    !!

    subroutine PDM_connectivity_transpose_ (n_entity1,           &
                                            n_entity2,           &
                                            entity1_entity2_idx, &
                                            entity1_entity2,     &
                                            entity2_entity1_idx, &
                                            entity2_entity1)

      use iso_c_binding

      implicit none

      integer, intent(in) :: n_entity1
      integer, intent(in) :: n_entity2

      integer, pointer    :: entity1_entity2_idx(:)
      integer, pointer    :: entity1_entity2(:)
      integer, pointer    :: entity2_entity1_idx(:)
      integer, pointer    :: entity2_entity1(:)

      integer(c_int) :: c_n_entity1
      integer(c_int) :: c_n_entity2

      type(c_ptr)    :: c_entity1_entity2_idx     = C_NULL_PTR
      type(c_ptr)    :: c_entity1_entity2         = C_NULL_PTR
      type(c_ptr)    :: c_entity2_entity1_idx     = C_NULL_PTR
      type(c_ptr)    :: c_entity2_entity1         = C_NULL_PTR

      c_n_entity1 = n_entity1
      c_n_entity2 = n_entity2

      c_entity1_entity2_idx = c_loc(entity1_entity2_idx)
      c_entity1_entity2     = c_loc(entity1_entity2)

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

    end subroutine PDM_connectivity_transpose_

end module pdm_part_connectivity_transform
