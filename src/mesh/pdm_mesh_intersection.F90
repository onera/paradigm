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

module pdm_mesh_intersection

  use pdm
  use pdm_fortran
  use iso_c_binding
  use pdm_pointer_array

  implicit none

  interface 

  end interface

contains

! PDM_mesh_intersection_t*
! PDM_mesh_intersection_create
! (
!  const PDM_mesh_intersection_kind_t intersection_kind,
!  const int                          dim_mesh_a,
!  const int                          dim_mesh_b,
!  const double                       project_coeff,
!        PDM_MPI_Comm                 comm,
!  const PDM_ownership_t              owner
! );

! /**
!  *
!  * \brief Set global data of a mesh
!  *
!  * \param [in]   mi             Pointer to \ref PDM_mesh_intersection object
!  * \param [in]   i_mesh         Mesh identifier
!  * \param [in]   n_part         Number of partitions
!  *
!  */

! void
! PDM_mesh_intersection_n_part_set
! (
!   PDM_mesh_intersection_t *mi,
!   const int                i_mesh,
!   const int                n_part
! );

! void
! PDM_mesh_intersection_compute
! (
!   PDM_mesh_intersection_t  *mi
! );

! /**
!  *
!  * \brief Set the tolerance for bounding boxes
!  *
!  * \param [in]   mi              Pointer to \ref PDM_mesh_intersection object
!  * \param [in]   tol             Tolerance
!  *
!  */
! void
! PDM_mesh_intersection_tolerance_set
! (
!        PDM_mesh_intersection_t *mi,
!  const double                   tol
! );

! void
! PDM_mesh_intersection_part_set
! (
!   PDM_mesh_intersection_t  *mi,
!   int                       i_mesh,
!   int                       i_part,
!   int                       n_cell,
!   int                       n_face,
!   int                       n_edge,
!   int                       n_vtx,
!   int                      *cell_face_idx,
!   int                      *cell_face,
!   int                      *face_edge_idx,
!   int                      *face_edge,
!   int                      *edge_vtx,
!   int                      *face_vtx_idx,
!   int                      *face_vtx,
!   PDM_g_num_t              *cell_ln_to_gn,
!   PDM_g_num_t              *face_ln_to_gn,
!   PDM_g_num_t              *edge_ln_to_gn,
!   PDM_g_num_t              *vtx_ln_to_gn,
!   double                   *vtx_coord
! );


! void
! PDM_mesh_intersection_free
! (
!  PDM_mesh_intersection_t* mi
! );


! /**
!  * \brief Get part_to_part object to exchange data between the intersected meshes
!  *
!  * \param [in ] mi         Pointer to \ref PDM_mesh_intersection_t object
!  * \param [out] ptp        Pointer to \ref PDM_part_to_part_t object
!  * \param [in ] ownership  Ownership for ptp
!  *
!  */

! void
! PDM_mesh_intersection_part_to_part_get
! (
!  PDM_mesh_intersection_t  *mi,
!  PDM_part_to_part_t      **ptp,
!  PDM_ownership_t           ownership
!  );


! void
! PDM_mesh_intersection_result_from_a_get
! (
!        PDM_mesh_intersection_t  *mi,
!  const int                       ipart,
!        int                     **elt_a_elt_b_idx,
!        PDM_g_num_t             **elt_a_elt_b,
!        double                  **elt_a_elt_b_volume
! );

! void
! PDM_mesh_intersection_result_from_b_get
! (
!        PDM_mesh_intersection_t  *mi,
!  const int                       ipart,
!        double                  **elt_b_elt_a_volume
!  );

! void
! PDM_mesh_intersection_elt_volume_get
! (
!        PDM_mesh_intersection_t  *mi,
!  const int                       imesh,
!  const int                       ipart,
!        double                  **elt_volume
!  );

end module pdm_mesh_intersection


