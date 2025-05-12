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


program pdm_t_gnum_f

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use PDM_gnum
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  integer,                   parameter :: comm  = MPI_COMM_WORLD
  integer,                   parameter :: n_pts = 10

  type(c_ptr)                          :: gen_gnum = C_NULL_PTR
  double precision,          pointer   :: pts_coord(:,:)  => null()
  integer(kind=pdm_g_num_s), pointer   :: pts_ln_to_gn(:) => null()
  double precision,          pointer   :: char_length(:)  => null()

  integer                              :: code
  integer                              :: i_rank
  integer                              :: n_rank
  !-----------------------------------------------------------


  ! Initialize MPI
  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)


  ! Generate random point cloud
  allocate(pts_coord(3,n_pts))
  allocate(char_length(n_pts))

  call random_seed()
  call random_number(pts_coord)

  ! Define point-wise characteristic length
  char_length(:) = 1.d-6


  ! Create PDM_gen_gnum object
  call pdm_gnum_create(gen_gnum,           & ! -> GenGnum object
                       3,                  & ! <- Dimension
                       1,                  & ! <- Number of parts on current MPI rank
                       0,                  & ! <- Merge duplicate points
                       1.d-3,              & ! <- Global tolerance
                       comm,               & ! <- MPI communicator
                       PDM_OWNERSHIP_USER)   ! <- Ownership of result

  ! Set coordinates
  call pdm_gnum_set_from_coords(gen_gnum,    & ! <- GenGnum object
                                0,           & ! <- Part identifier
                                n_pts,       & ! <- Number of points
                                pts_coord,   & ! <- Coordinates
                                char_length)   ! <- Characteristic length


  ! Compute global numbering
  call pdm_gnum_compute(gen_gnum)


  ! Get global IDs
  call pdm_gnum_get(gen_gnum,     & ! <- GenGnum object
                    0,            & ! <- Part identifier
                    pts_ln_to_gn)   ! -> Global IDs

  ! Free memory
  call pdm_gnum_free(gen_gnum)

  deallocate(pts_coord)
  deallocate(char_length)



  ! Free C-allocated memory
  call pdm_fortran_free_c(c_loc(pts_ln_to_gn))


  if (i_rank .eq. 0) then
    write(*, *) "-- End"
  end if

  call mpi_finalize(code)


end program pdm_t_gnum_f
