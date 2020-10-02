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

program testf

  use pdm
  use mpi
  use pdm_block_to_part
  use iso_c_binding

  implicit none

  integer :: code
  integer :: i_rank
  integer :: n_rank

  integer (kind = pdm_g_num_s), pointer :: block_distrib_index(:)
  type(c_ptr)                           :: cptr_block_distrib_index
  integer (kind = pdm_g_num_s), pointer :: gnum_elt(:)
  type(c_ptr), pointer                  :: cptr_gnum_elt(:)
  type(c_ptr)                           :: cptr_cptr_gnum_elt

  integer(c_int), parameter :: fComm = MPI_COMM_WORLD

  type(c_ptr)              :: btp

  integer(c_int), pointer   :: n_elt(:)
  type(c_ptr)              :: cptr_n_elt

  integer(c_int), parameter :: n_part = 1

  integer(c_int), parameter :: n_elt_case = 5

  integer(c_int), parameter :: t_stride = PDM_STRIDE_VAR

  type(c_ptr)          :: cptr_block_stride !  c_ptr of C array containing block stride
                                            ! (NULL for this case)
  type(c_ptr)          :: cptr_block_data  ! c_ptr of C array containing block data
  integer(c_int), pointer :: block_data(:)

  integer(c_int), pointer :: part_stride(:) ! Fortran array containing data for part 0
                                            ! (Unused for this case with a constant stride)
  type(c_ptr), pointer :: cptr_part_stride(:) ! Fortran array of c_ptr containing c_loc (part_stride)
                                              ! for each partition (for this case n_part = 1)
                                              ! (Unused for this case with a constant stride)
  type(c_ptr)          :: cptr_cptr_part_stride ! c_loc about cptr_part_stride

  integer(c_size_t), parameter :: s_data = 8 ! size of type of the exchanged data
                                             ! (for this case exchange data are 'real*8')

  integer(c_int), pointer :: part_data(:) ! Fortran array containing data for part 0
  type(c_ptr), pointer :: cptr_part_data(:) ! Fortran array of c_ptr containing c_loc (part_data)
                                            ! for each partition (for this case n_part = 1)
  type(c_ptr)          :: cptr_cptr_part_data ! c_loc about cptr_part_data

  !
  ! Init
  !

  call mpi_init(code)
  call mpi_comm_rank(mpi_comm_world, i_rank, code)
  call mpi_comm_size(mpi_comm_world, n_rank, code)

  if (n_rank .ne. 1) then
    print *,'Error : 2 MPI processes are mandatory'
    call mpi_finalize(code)
    stop
  end if

  allocate(block_distrib_index(2))
  block_distrib_index(1) = 1
  block_distrib_index(2) = 6

  cptr_block_distrib_index = c_loc(block_distrib_index)

  allocate(cptr_gnum_elt(n_part))
  allocate(gnum_elt(n_elt_case))
  gnum_elt(1) = 3
  gnum_elt(2) = 5
  gnum_elt(3) = 2
  gnum_elt(4) = 4
  gnum_elt(5) = 1

  cptr_gnum_elt(1) = c_loc(gnum_elt)
  cptr_cptr_gnum_elt = c_loc(cptr_gnum_elt)

  allocate(n_elt(n_part))
  n_elt(1) = n_elt_case
  cptr_n_elt = c_loc(n_elt)

  btp = PDM_block_to_part_create (blockDistribIdx, &
                                  gnum_elt, &
                                  n_elt, &
                                  n_part, &
                                  fcomm)

  call PDM_block_to_part_exch2 (btp, &
                                s_data, &
                                t_stride, &
                                block_stride, &
                                block_data, &
                                part_stride, &
                                part_data)

  btp = PDM_block_to_part_free (btp)

  part_stride => NULL()
  cptr_part_stride => NULL()
  cptr_cptr_part_stride = C_NULL_PTR
  cptr_block_stride = C_NULL_PTR
  cptr_block_data = C_NULL_PTR

  allocate(part_data(n_elt_case))
  part_data(1) = 13
  part_data(2) = 15
  part_data(3) = 12
  part_data(4) = 14
  part_data(5) = 11

  allocate(cptr_part_data(n_part))
  cptr_part_data(1) = c_loc(part_data)

  cptr_cptr_part_data = c_loc(cptr_part_data)

  size_highest_block = PDM_part_to_block_exch (ptb,&
                                               s_data, &
                                               t_stride, &
                                               cst_stride, &
                                               cptr_cptr_part_stride, &
                                               cptr_cptr_part_data, &
                                               cptr_block_stride, &
                                               cptr_block_data)

  call c_f_pointer(cptr_block_data, block_data, [n_elt_case])
  print *, block_data(1), block_data(2), block_data(3), block_data(4), block_data(5)

  ptb = PDM_part_to_block_free (ptb)

  deallocate(cptr_part_data)
  deallocate(part_data)
  deallocate(n_elt)
  deallocate(cptr_gnum_elt)
  deallocate(gnum_elt)

  call mpi_finalize(code)

end program testf
