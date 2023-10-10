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

program pdm_t_mesh_partitioning_sol_f

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use pdm_multipart
  use pdm_vtk
  use pdm_dcube_nodal_gen
  use pdm_dmesh_nodal
  use pdm_part_mesh_nodal
  use pdm_mesh_nodal
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  integer (c_int)                       :: i
  integer (c_int)                       :: fe = 0

  ! MPI
  integer                               :: code
  integer                               :: i_rank
  integer                               :: n_rank
  integer, parameter                    :: comm = MPI_COMM_WORLD

  ! PDM_dcube_nodal_gen_create
  integer                               :: n_x, n_y, n_z
  integer                               :: elt_type, order
  double precision                      :: length
  double precision                      :: xmin, ymin, zmin
  type (c_ptr)                          :: dcube

  ! PDM_dcube_nodal_gen_dmesh_nodal_get
  type (c_ptr)                          :: dmn

  ! PDM_multipart_create
  type (c_ptr)                          :: mpart
  integer (c_int)                       :: n_zone = 1
  integer(kind=PDM_l_num_s), pointer    :: n_part(:) => null()
  integer (c_int)                       :: part_method
  double precision,          pointer    :: part_fraction(:) => null()

  integer (c_int)                       :: i_zone    = 0
  integer (c_int)                       :: i_section = 0
  integer (c_int)                       :: i_part    = 0

  ! PDM_multipart_set_reordering_options
  integer(kind=PDM_l_num_s), pointer    :: renum_cell_properties(:) => null()

  ! PDM_multipart_get_part_mesh_nodal
  type (c_ptr)                          :: pmn

  ! PDM_part_mesh_nodal_section_n_elt_get
  integer (c_int)                       :: n_elt = -1

  ! PDM_part_mesh_nodal_section_std_get
  integer(kind=PDM_l_num_s), pointer  :: elt_vtx(:)
  integer (pdm_g_num_s), pointer      :: elt_ln_to_gn(:)
  integer(kind=PDM_l_num_s), pointer  :: parent_num(:)
  integer (pdm_g_num_s), pointer      :: parent_entity_g_num(:)

  ! PDM_multipart_part_vtx_coord_get
  double precision, pointer           :: coords(:,:)
  integer(c_int)                      :: n_vtx = -1

  ! PDM_part_mesh_nodal_vtx_g_num_get
  integer (pdm_g_num_s), pointer      :: vtx_ln_to_gn(:)

  ! FV
  integer(kind=PDM_g_num_s), pointer :: edge_ln_to_gn(:)
  integer(c_int)                     :: n_edge
  integer(kind=PDM_l_num_s), pointer :: edge_vtx(:)
  integer(kind=PDM_l_num_s), pointer :: edge_vtx_idx(:)

  integer(kind=PDM_g_num_s), pointer :: face_ln_to_gn(:)
  integer(c_int)                     :: n_face
  integer(kind=PDM_l_num_s), pointer :: face_edge(:)
  integer(kind=PDM_l_num_s), pointer :: face_edge_idx(:)

  integer(kind=PDM_g_num_s), pointer :: cell_ln_to_gn(:)
  integer(c_int)                     :: n_cell
  integer(kind=PDM_l_num_s), pointer :: cell_face(:)
  integer(kind=PDM_l_num_s), pointer :: cell_face_idx(:)
  !-----------------------------------------------------------

  ! Initialize MPI environment
  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)

  ! Generate block-distributed parallelepided mesh
  n_x      = 10
  n_y      = 10
  n_z      = 10
  length   = 1.
  xmin     = 0.
  ymin     = 0.
  zmin     = 0.
  elt_type = PDM_MESH_NODAL_TETRA4
  order    = 1
  call PDM_dcube_nodal_gen_create(dcube,     &
                                  comm,      &
                                  n_x,       &
                                  n_y,       &
                                  n_z,       &
                                  length,    &
                                  xmin,      &
                                  ymin,      &
                                  zmin,      &
                                  elt_type,  &
                                  order,     &
                                  PDM_OWNERSHIP_USER)

  call PDM_dcube_nodal_gen_build(dcube, dmn)

  call PDM_dcube_nodal_gen_dmesh_nodal_get(dcube, dmn)

  call PDM_dmesh_nodal_generate_distribution(dmn)

  call PDM_dcube_nodal_gen_free(dcube)

  ! Create partitioning object
  allocate(n_part(n_zone))

  do i = 1, n_zone
    n_part(i) = 1
  end do

  allocate(part_fraction(n_part(i_zone+1)))

  do i = 1, n_part(i_zone+1)
    part_fraction(i) = 1
  end do

  part_method = PDM_SPLIT_DUAL_WITH_HILBERT
  call PDM_multipart_create(mpart,                     &
                            n_zone,                    &
                            n_part,                    &
                            PDM_FALSE,                 &
                            part_method,               &
                            PDM_PART_SIZE_HOMOGENEOUS, &
                            part_fraction,             &
                            comm,                      &
                            PDM_OWNERSHIP_KEEP)

  call PDM_multipart_set_reordering_options(mpart,                      &
                                            i_zone,                     &
                                            "PDM_PART_RENUM_CELL_NONE", &
                                            renum_cell_properties,      &
                                            "PDM_PART_RENUM_FACE_NONE")

  call PDM_multipart_register_dmesh_nodal(mpart,  &
                                          i_zone, &
                                          dmn)

  call PDM_multipart_run_ppart(mpart)

  ! Get mesh arrrays in FE structure
  if (fe .eq. 1) then
    call PDM_multipart_get_part_mesh_nodal(mpart,  &
                                           i_zone, &
                                           pmn,    &
                                           PDM_OWNERSHIP_USER)

    call PDM_part_mesh_nodal_section_n_elt_get(pmn,       &
                                               i_section, &
                                               i_part,    &
                                               n_elt)

    call PDM_part_mesh_nodal_section_std_get(pmn,                 &
                                             i_section,           &
                                             i_part,              &
                                             elt_vtx,             &
                                             elt_ln_to_gn,        &
                                             parent_num,          &
                                             parent_entity_g_num, &
                                             PDM_OWNERSHIP_KEEP)

    call PDM_multipart_part_vtx_coord_get(mpart,              &
                                          i_zone,             &
                                          i_part,             &
                                          coords,             &
                                          PDM_OWNERSHIP_USER, &
                                          n_vtx)

    call PDM_part_mesh_nodal_vtx_g_num_get(pmn,    &
                                           i_part, &
                                           vtx_ln_to_gn)

    call PDM_part_mesh_nodal_free(pmn)
  end if

  ! Get mesh arrrays in FV structure
  if (fe .eq. 0) then
    call PDM_multipart_part_ln_to_gn_get(mpart,                  &
                                         i_zone,                 &
                                         i_part,                 &
                                         PDM_MESH_ENTITY_VERTEX, &
                                         vtx_ln_to_gn,           &
                                         PDM_OWNERSHIP_KEEP,     &
                                         n_vtx)

    call PDM_multipart_part_vtx_coord_get(mpart,             &
                                         i_zone,             &
                                         i_part,             &
                                         coords,             &
                                         PDM_OWNERSHIP_KEEP, &
                                         n_vtx)

    call PDM_multipart_part_ln_to_gn_get(mpart,                &
                                         i_zone,               &
                                         i_part,               &
                                         PDM_MESH_ENTITY_EDGE, &
                                         edge_ln_to_gn,        &
                                         PDM_OWNERSHIP_KEEP,   &
                                         n_edge)

    call PDM_multipart_part_connectivity_get(mpart,                          &
                                             i_zone,                         &
                                             i_part,                         &
                                             PDM_CONNECTIVITY_TYPE_EDGE_VTX, &
                                             edge_vtx,                       &
                                             edge_vtx_idx,                   &
                                             PDM_OWNERSHIP_KEEP,             &
                                             n_edge)

    call PDM_multipart_part_ln_to_gn_get(mpart,                &
                                         i_zone,               &
                                         i_part,               &
                                         PDM_MESH_ENTITY_FACE, &
                                         face_ln_to_gn,        &
                                         PDM_OWNERSHIP_KEEP,   &
                                         n_face)

    call PDM_multipart_part_connectivity_get(mpart,                           &
                                             i_zone,                          &
                                             i_part,                          &
                                             PDM_CONNECTIVITY_TYPE_FACE_EDGE, &
                                             face_edge,                       &
                                             face_edge_idx,                   &
                                             PDM_OWNERSHIP_KEEP,              &
                                             n_face)

    call PDM_multipart_part_ln_to_gn_get(mpart,                &
                                         i_zone,               &
                                         i_part,               &
                                         PDM_MESH_ENTITY_CELL, &
                                         cell_ln_to_gn,        &
                                         PDM_OWNERSHIP_KEEP,   &
                                         n_cell)

    call PDM_multipart_part_connectivity_get(mpart,                           &
                                             i_zone,                          &
                                             i_part,                          &
                                             PDM_CONNECTIVITY_TYPE_CELL_FACE, &
                                             cell_face,                       &
                                             cell_face_idx,                   &
                                             PDM_OWNERSHIP_KEEP,              &
                                             n_cell)
  end if

  ! free
  deallocate(n_part, &
             part_fraction)
  call PDM_DMesh_nodal_free(dmn)
  call PDM_multipart_free(mpart)

  ! Finalize MPI environment
  call mpi_finalize(code)

end program pdm_t_mesh_partitioning_sol_f
