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

program testf

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use pdm_multipart
  use pdm_dcube_gen
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  ! MPI
  integer,                  parameter :: comm = MPI_COMM_WORLD
  integer                             :: code
  integer                             :: i_rank
  integer                             :: n_rank
  ! UTIL
  integer                             :: i
  ! MULTIPART
  type(c_ptr)                         :: multipart = C_NULL_PTR
  integer(c_int)                      :: split_method
  integer(c_int)                      :: n_part = 1
  integer(c_int)                      :: n_zone = 1
  integer(kind=PDM_l_num_s), pointer  :: n_part_zones(:)  => null()
  double precision,          pointer  :: part_fraction(:) => null()
  !-----------------------------------------------------------

  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)

  ! Initialize multipart

#ifdef PDM_HAVE_PARMETIS
  split_method  = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  split_method  = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#endif
#endif

  allocate(n_part_zones(n_zone))

  do i = 1, n_zone
    n_part_zones(i) = n_part
  end do

  if (i_rank .eq. 0) then
    write(*, *) "PDM_multipart_create"
  end if

  call PDM_multipart_create(multipart, &
                            n_zone, &
                            n_part_zones, &
                            PDM_FALSE, &
                            split_method, &
                            PDM_PART_SIZE_HOMOGENEOUS, &
                            part_fraction, &
                            comm, &
                            PDM_OWNERSHIP_KEEP)

  ! Generate Mesh

  ! Run

  ! Free
  if (i_rank .eq. 0) then
    write(*, *) "PDM_multipart_free"
  end if

  call PDM_multipart_free(multipart)

  if (i_rank .eq. 0) then
    write(*, *) "-- End"
  end if

  call mpi_finalize(code)

end program testf
