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
  use pdm_mesh_nodal
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  integer (c_int)                       :: i

  ! MPI
  integer                               :: code
  integer                               :: i_rank
  integer                               :: n_rank
  integer, parameter                    :: comm = MPI_COMM_WORLD

  ! PDM_dcube_nodal_gen_create
  integer                               :: n_x, n_y, n_z
  integer                               :: elt_type, order, ownership
  double precision                      :: length
  double precision                      :: xmin, ymin, zmin
  type (c_ptr)                          :: dcube

  ! PDM_dcube_nodal_gen_dmesh_nodal_get
  type (c_ptr)                          :: dmn

  ! PDM_multipart_create
  type (c_ptr)                          :: mpart
  integer (c_int)                       :: n_zone = 1
  integer(kind=PDM_l_num_s), pointer    :: n_part(:) => null()
  integer (c_int)                       :: merge_blocks
  integer (c_int)                       :: split_method
  integer (c_int)                       :: part_size_method
  double precision,          pointer    :: part_fraction(:) => null()

  integer (c_int)                       :: i_zone    = 0
  integer (c_int)                       :: i_section = 0
  integer (c_int)                       :: i_part    = 0

  ! PDM_multipart_set_reordering_options
  integer(kind=PDM_l_num_s), pointer    :: renum_cell_properties(:) => null()

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
  elt_type = 4 ! TO DO : how to use PDM_MESH_NODAL_TETRA4 in Fortran ?
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
                                  ownership) ! TO DO PDM_OWNERSHIP_USER

  call PDM_dcube_nodal_gen_build(dcube, dmn)

  call PDM_dcube_nodal_gen_dmesh_nodal_get(dcube, dmn)

  ! call PDM_dmesh_nodal_generate_distribution(dmn)
  ! TO DO : interface Fortran de dmesh_nodal

  call PDM_dcube_nodal_gen_free(dcube)

  ! Create partitioning object
  allocate(n_part(n_zone))

  do i = 1, n_zone
    n_part(i) = 1
  end do

  call PDM_multipart_create(mpart,            &
                            n_zone,           &
                            n_part,           &
                            merge_blocks,     &
                            split_method,     &
                            part_size_method, &
                            part_fraction,    &
                            comm,             &
                            ownership) ! TO DO PDM_OWNERSHIP_KEEP

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

  ! Get mesh arrrays in FV structure

  ! free
  ! call PDM_DMesh_nodal_free(dmn)
  call PDM_multipart_free(mpart)

  ! Finalize MPI environment
  call mpi_finalize(code)

end program pdm_t_mesh_partitioning_sol_f
