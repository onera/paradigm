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

module pdm_gnum

  use pdm

  implicit none

  interface PDM_gnum_create ; module procedure &
  pdm_gnum_create_
  end interface

  private :: pdm_gnum_create_

  interface

  !>
  !!
  !! \brief Build a global numbering structure
  !!
  !! \param [in]   dim          Spatial dimension
  !! \param [in]   n_part       Number of local partitions
  !! \param [in]   merge        Merge double points or not
  !! \param [in]   tolerance    Geometric tolerance (if merge double points is activated)
  !! \param [in]   comm         PDM_MPI communicator
  !!
  !! \return     Pointer to \ref PDM_gen_gnum object
  !!

  function PDM_gnum_create_cf (dim,        &
                               n_part,     &
                               merge,      &
                               tolerance,  &
                               comm,       &
                               owner)      &
  result (gen_gnum) &

  bind (c, name = 'PDM_gnum_create')

  use iso_c_binding
  implicit none

  integer(c_int), value :: dim
  integer(c_int), value :: n_part
  integer(c_int), value :: merge
  real(c_double), value :: tolerance
  integer(c_int), value :: comm
  integer(c_int), value :: owner

  type (c_ptr) :: gen_gnum

  end function PDM_gnum_create_cf


  !>
  !!
  !! \brief Set from coordinates
  !!
  !! \param [in]   gen_gnum     Pointer to \ref PDM_gen_gnum object
  !! \param [in]   i_part       Current partition
  !! \param [in]   n_elts       Number of elements
  !! \param [in]   coords       Coordinates (size = 3 * \ref n_elts)
  !! \param [in]   char_length  Characteristic length (or NULL)
  !!                            (used if merge double points is activated)
  !!
  !!

  subroutine PDM_gnum_set_from_coords (gen_gnum,   &
                                       i_part,     &
                                       n_elts,     &
                                       coords,     &
                                       char_length)&
  bind (c, name = 'PDM_gnum_set_from_coords')

  use iso_c_binding
  implicit none

  type(c_ptr), value :: gen_gnum

  integer(c_int), value :: i_part
  integer(c_int), value :: n_elts
  type(c_ptr), value :: coords
  type(c_ptr), value :: char_length

  end subroutine PDM_gnum_set_from_coords


  !>
  !!
  !! \brief Set Parent global numbering
  !!
  !! \param [in]   gen_gnum     Pointer to \ref PDM_gen_gnum object
  !! \param [in]   i_part       Current partition
  !! \param [in]   n_elts       Number of elements
  !! \param [in]   parent_gnum  Parent global numbering (size = \ref n_elts)
  !!
  !!

  subroutine PDM_gnum_set_from_parents (gen_gnum,   &
                                        i_part,     &
                                        n_elts,     &
                                        parent_gnum)&
  bind (c, name = 'PDM_gnum_set_from_parents')

  use iso_c_binding
  implicit none

  type(c_ptr), value :: gen_gnum

  integer(c_int), value :: i_part
  integer(c_int), value :: n_elts
  type(c_ptr), value :: parent_gnum

  end subroutine PDM_gnum_set_from_parents


  !>
  !!
  !! \brief Compute
  !!
  !! \param [in]   gen_gnum         Pointer to \ref PDM_gen_gnum object
  !!
  !!

  subroutine PDM_gnum_compute (gen_gnum) &
  bind (c, name = 'PDM_gnum_compute')

  use iso_c_binding
  implicit none

  type(c_ptr), value :: gen_gnum

  end subroutine PDM_gnum_compute


  !>
  !!
  !! \brief Get global ids for a given partition
  !!
  !! \param [in]   gen_gnum     Pointer to \ref PDM_gen_gnum object
  !! \param [in]   i_part       Current partition
  !!
  !! \return     Array of global ids
  !!
  !!

  function PDM_gnum_get (gen_gnum, i_part) &
  result (g_nums) &
  bind (c, name = 'PDM_gnum_get')

  use iso_c_binding
  implicit none

  type(c_ptr) :: g_nums
  type(c_ptr), value :: gen_gnum
  integer(c_int), value :: i_part

  end function PDM_gnum_get


  !>
  !!
  !! \brief Free
  !!
  !! \param [in]   gen_gnum         Pointer to \ref PDM_gen_gnum object
  !!
  !!

  subroutine PDM_gnum_free (gen_gnum)     &
  bind (c, name = 'PDM_gnum_free')
  use iso_c_binding
  implicit none

  type (c_ptr)  , value :: gen_gnum

  end subroutine PDM_gnum_free

  end interface




  contains



  !>
  !!
  !! \brief Build a global numbering structure
  !!
  !! \param [out]  gen_gnum     Pointer to \ref PDM_gen_gnum object
  !! \param [in]   dim          Spatial dimension
  !! \param [in]   n_part       Number of local partitions
  !! \param [in]   merge        Merge double points or not
  !! \param [in]   tolerance    Geometric tolerance (if merge double points is activated)
  !! \param [in]   comm         PDM_MPI communicator
  !!

  subroutine PDM_gnum_create_ (gen_gnum,  &
                               dim,       &
                               n_part,    &
                               merge,     &
                               tolerance, &
                               f_comm,    &
                               owner)

  use iso_c_binding
  implicit none

  integer          :: dim
  integer          :: n_part
  integer          :: merge
  double precision :: tolerance
  integer          :: f_comm
  integer          :: owner

  type (c_ptr) :: gen_gnum

  integer(c_int) :: c_dim
  integer(c_int) :: c_n_part
  integer(c_int) :: c_merge
  real(c_double) :: c_tolerance
  integer(c_int) :: c_comm
  integer(c_int) :: c_owner


  c_comm = PDM_MPI_Comm_f2c(f_comm)

  c_dim       = dim
  c_n_part    = n_part
  c_merge     = merge
  c_tolerance = tolerance
  c_owner     = owner

  gen_gnum = PDM_gnum_create_cf (c_dim,       &
                                 c_n_part,    &
                                 c_merge,     &
                                 c_tolerance, &
                                 c_comm,      &
                                 c_owner)

  end subroutine PDM_gnum_create_

  end module pdm_gnum
