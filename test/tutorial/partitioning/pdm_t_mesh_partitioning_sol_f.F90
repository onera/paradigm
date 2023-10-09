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
  ! use pdm_dcube_nodal_gen
  use pdm_mesh_nodal
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  integer                               :: code
  integer                               :: i_rank
  integer                               :: n_rank
  integer, parameter                    :: comm = MPI_COMM_WORLD
  !-----------------------------------------------------------

  ! Initialize MPI environment
  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)

  ! TO DO

  ! Finalize MPI environment
  call mpi_finalize(code)

end program pdm_t_mesh_partitioning_sol_f
