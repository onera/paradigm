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

module pdm_dcube_nodal_gen

  use pdm

  implicit none

  interface PDM_dcube_nodal_gen_create ; module procedure &
  PDM_dcube_nodal_gen_create_
  end interface

  interface PDM_dcube_nodal_gen_build ; module procedure &
  PDM_dcube_nodal_gen_build_
  end interface

  interface PDM_dcube_nodal_gen_dmesh_nodal_get ; module procedure &
  PDM_dcube_nodal_gen_dmesh_nodal_get_
  end interface

  private :: PDM_dcube_nodal_gen_create_
  private :: PDM_dcube_nodal_gen_build_
  private :: PDM_dcube_nodal_gen_dmesh_nodal_get_

  interface

    !>
    !!
    !! \brief Create a distributed nodal cube mesh
    !!
    !! \param [in]   comm           Communicator
    !! \param [in]   n_x            Number of elements in x-direction
    !! \param [in]   n_y            Number of elements in y-direction
    !! \param [in]   n_z            Number of elements in z-direction
    !! \param [in]   length         Length of cube side
    !! \param [in]   xmin           Minimal x-coordinate
    !! \param [in]   ymin           Minimal y-coordinate
    !! \param [in]   zmin           Minimal z-coordinate
    !! \param [in]   elt_type       Element type
    !! \param [in]   order          Element order
    !! \param [in]   ownership      Object ownership
    !!

    function PDM_dcube_nodal_gen_create_cf (comm,      &
                                            n_x,       &
                                            n_y,       &
                                            n_z,       &
                                            length,    &
                                            xmin,      &
                                            ymin,      &
                                            zmin,      &
                                            elt_type,  &
                                            order,     &
                                            ownership) &
      result(dcube) &
      bind (c, name = 'PDM_dcube_nodal_gen_create')

      use iso_c_binding
      implicit none

      integer(c_int), value :: comm
      integer(c_int), value :: n_x, n_y, n_z
      integer(c_int), value :: elt_type, order, ownership

      real(c_double), value :: length
      real(c_double), value :: xmin, ymin, zmin

      type (c_ptr) :: dcube

    end function PDM_dcube_nodal_gen_create_cf

    !>
    !!
    !! \brief Build a \ref PDM_dcube_nodal_t structure
    !!
    !! \param [in] dcube           Pointer to \ref PDM_dcube_nodal_t object
    !!

    function PDM_dcube_nodal_gen_build_cf(dcube) &

      result(dmn) &
      bind (c, name = 'PDM_dcube_nodal_gen_build')

      use iso_c_binding
      implicit none

      type (c_ptr), value  :: dcube
      type (c_ptr)         :: dmn

    end function PDM_dcube_nodal_gen_build_cf

    !>
    !!
    !! \brief Get the \ref PDM_dmesh_nodal_t associated to a \ref PDM_dcube_nodal_t
    !!
    !! \param [in] dcube           Pointer to \ref PDM_dcube_nodal_t object
    !!

    function PDM_dcube_nodal_gen_dmesh_nodal_get_cf(dcube) &

      result(dmn) &
      bind (c, name = 'PDM_dcube_nodal_gen_dmesh_nodal_get')

      use iso_c_binding
      implicit none

      type (c_ptr), value  :: dcube
      type (c_ptr)         :: dmn

    end function PDM_dcube_nodal_gen_dmesh_nodal_get_cf

    !>
    !!
    !! \brief Free the structure
    !!
    !! \param [in]   dcube   Pointer to \ref PDM_dcube_nodal_t object
    !!

    subroutine PDM_dcube_nodal_gen_free (dcube) &
    bind(c, name='PDM_dcube_nodal_gen_free')

      use iso_c_binding
      implicit none

      type (c_ptr), value :: dcube
    end subroutine PDM_dcube_nodal_gen_free

  end interface

  contains

    !>
    !!
    !! \brief Create a distributed nodal cube mesh
    !!
    !! \param [in]   comm           Communicator
    !! \param [in]   n_x            Number of elements in x-direction
    !! \param [in]   n_y            Number of elements in y-direction
    !! \param [in]   n_z            Number of elements in z-direction
    !! \param [in]   length         Length of cube side
    !! \param [in]   xmin           Minimal x-coordinate
    !! \param [in]   ymin           Minimal y-coordinate
    !! \param [in]   zmin           Minimal z-coordinate
    !! \param [in]   elt_type       Element type
    !! \param [in]   order          Element order
    !! \param [in]   ownership      Object ownership
    !!
    !! \param [out] dcube           Pointer to \ref PDM_dcube_nodal_t object
    !!

    subroutine PDM_dcube_nodal_gen_create_(dcube,     &
                                           f_comm,    &
                                           n_x,       &
                                           n_y,       &
                                           n_z,       &
                                           length,    &
                                           xmin,      &
                                           ymin,      &
                                           zmin,      &
                                           elt_type,  &
                                           order,     &
                                           ownership)

      use iso_c_binding
      implicit none

      integer,              intent(in) :: f_comm
      integer,              intent(in) :: n_x, n_y, n_z
      integer,              intent(in) :: elt_type, order, ownership

      double precision,     intent(in) :: length
      double precision,     intent(in) :: xmin, ymin, zmin

      integer(c_int) :: c_comm
      integer(c_int) :: c_n_x, c_n_y, c_n_z
      integer(c_int) :: c_elt_type, c_order, c_ownership

      real(c_double) :: c_length
      real(c_double) :: c_xmin, c_ymin, c_zmin

      type (c_ptr)   :: dcube

      c_comm = PDM_MPI_Comm_f2c(f_comm)

      c_n_x       = n_x
      c_n_y       = n_y
      c_n_z       = n_z
      c_length    = length
      c_xmin      = xmin
      c_ymin      = ymin
      c_zmin      = zmin
      c_elt_type  = elt_type
      c_order     = order
      c_ownership = ownership

      dcube = PDM_dcube_nodal_gen_create_cf(c_comm,     &
                                            c_n_x,      &
                                            c_n_y,      &
                                            c_n_z,      &
                                            c_length,   &
                                            c_xmin,     &
                                            c_ymin,     &
                                            c_zmin,     &
                                            c_elt_type, &
                                            c_order,    &
                                            c_ownership)

    end subroutine PDM_dcube_nodal_gen_create_

    !>
    !!
    !! \brief Build a \ref PDM_dcube_nodal_t structure
    !!
    !! \param [in] dcube           Pointer to \ref PDM_dcube_nodal_t object
    !!
    !! \param [out] dmn            Pointer to \ref PDM_mesh_nodal_t object
    !!

    subroutine PDM_dcube_nodal_gen_build_(dcube, &
                                          dmn)

      use iso_c_binding
      implicit none

      type(c_ptr), value   :: dcube

      type (c_ptr)         :: dmn

      dmn = PDM_dcube_nodal_gen_build_cf(dcube)

    end subroutine PDM_dcube_nodal_gen_build_

    !>
    !!
    !! \brief Get the \ref PDM_dmesh_nodal_t associated to a \ref PDM_dcube_nodal_t
    !!
    !! \param [in] dcube           Pointer to \ref PDM_dcube_nodal_t object
    !!
    !! \param [out] dmn            Pointer to \ref PDM_mesh_nodal_t object
    !!

    subroutine PDM_dcube_nodal_gen_dmesh_nodal_get_(dcube, &
                                                    dmn)

      use iso_c_binding
      implicit none

      type(c_ptr), value   :: dcube

      type (c_ptr)         :: dmn

      dmn = PDM_dcube_nodal_gen_dmesh_nodal_get_cf(dcube)

    end subroutine PDM_dcube_nodal_gen_dmesh_nodal_get_

end module pdm_dcube_nodal_gen
