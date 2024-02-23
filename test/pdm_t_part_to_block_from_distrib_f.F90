!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2024  ONERA
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
  use pdm_pointer_array
  use pdm_part_to_block
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  integer, parameter                    :: comm = MPI_COMM_WORLD
  integer, parameter                    :: n_low  = 10
  integer, parameter                    :: n_high = 20

  character(len=99)                     :: arg

  integer                               :: code
  integer                               :: i_rank
  integer                               :: n_rank

  integer                               :: n_part
  integer(pdm_l_num_s),      pointer    :: n_elt(:) => null()
  type(PDM_pointer_array_t), pointer    :: gnum_elt => null()
  integer(pdm_l_num_s)                  :: dn_elt
  integer(pdm_g_num_s),      pointer    :: ln_to_gn(:) => null()

  type(c_ptr)                           :: ptb = C_NULL_PTR
  integer(pdm_l_num_s),      pointer    :: all_dn_elt(:) => null()
  integer(pdm_g_num_s),      pointer    :: data_distrib_idx(:) => null()

  type(PDM_pointer_array_t), pointer    :: part_data => null()
  integer(pdm_l_num_s),      pointer    :: data(:)   => null()

  integer(pdm_l_num_s),      pointer    :: block_stride(:) => null()
  integer(pdm_l_num_s),      pointer    :: block_data(:)   => null()

  integer(pdm_l_num_s)                  :: expected

  integer                               :: i, j, k
  !-----------------------------------------------------------


  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)

  !  Default values
  n_part = 1

  !  Read command line arguments
  i = 1
  do while (i <= command_argument_count())
    call get_command_argument(i, arg)
    select case(arg)
      case ("-n_part")
        i = i+1
        call get_command_argument(i, arg)
        read(arg, *) n_part
      case default
        print *, "Invalid command argument ", arg
        stop
    end select

    i = i + 1
  enddo


  !  Define partitions

  allocate(n_elt(n_part))
  dn_elt = 0
  do i = 1, n_part
    n_elt(i) = randint(n_low, n_high)
    dn_elt = dn_elt + n_elt(i)
  enddo

  !  Create distribution
  allocate(all_dn_elt(n_rank))
  call MPI_Allgather(dn_elt,     1, MPI_INT, &
                     all_dn_elt, 1, MPI_INT, &
                     comm, code)

  allocate(data_distrib_idx(n_rank+1))
  data_distrib_idx(1) = 0
  do i = 1, n_rank
    data_distrib_idx(i+1) = data_distrib_idx(i) + all_dn_elt(i)
  enddo
  deallocate(all_dn_elt)


  if (i_rank == 0) then
    write (*,*) "data_distrib_idx =", data_distrib_idx
  endif


  !  Create global ids (== block)
  call PDM_pointer_array_create(gnum_elt,       &
                                n_part,         &
                                PDM_TYPE_G_NUM)

  k = 0
  do i = 1, n_part
    allocate(ln_to_gn(n_elt(i)))
    do j = 1, n_elt(i)
      ln_to_gn(j) = data_distrib_idx(i_rank+1) + k + 1
      k = k + 1
    enddo

    call PDM_pointer_array_part_set(gnum_elt, &
                                    i-1,      &
                                    ln_to_gn)

    ! write (*,*) "ln_to_gn =", ln_to_gn
  enddo



  !  Create Part-to-block object
  call PDM_part_to_block_create_from_distrib(ptb,                                &
                                             PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC, & ! t_distrib
                                             PDM_PART_TO_BLOCK_POST_CLEANUP,     & ! t_post
                                             1.d0,                               & ! partActiveNode
                                             gnum_elt,                           &
                                             data_distrib_idx,                   &
                                             n_elt,                              &
                                             n_part,                             &
                                             comm)


  !  Exchange
  call PDM_pointer_array_create(part_data,    &
                                n_part,       &
                                PDM_TYPE_INT)

  do i = 1, n_part
    call PDM_pointer_array_part_get(gnum_elt, &
                                    i-1,      &
                                    ln_to_gn)
    allocate(data(n_elt(i)))
    do j = 1, n_elt(i)
      data(j) = 2*int(ln_to_gn(j), kind=pdm_l_num_s)
    enddo

    call PDM_pointer_array_part_set(part_data, &
                                    i-1,       &
                                    data)

    ! write (*,*) "data =", data
  enddo



  call PDM_part_to_block_exch(ptb,                       &
                              PDM_STRIDE_CST_INTERLACED, & ! t_stride
                              1,                         & ! cst_stride
                              null(),                    &
                              part_data,                 &
                              block_stride,              &
                              block_data)



  do i = 1, dn_elt
    expected = 2*(int(data_distrib_idx(i_rank+1), kind=pdm_l_num_s + i)
    if (block_data(i) /= expected) then
      write (*,*) data_distrib_idx(i_rank+1) + i, "expected", expected, " but got", block_data(i)
    endif
  enddo


  !  Free memory
  call PDM_part_to_block_free(ptb)
  ! Leaks if n_part > 1 but invalid deallocate in intel (╯°□°）╯︵ ┻━┻
  ! do i = 1, n_part
  !   call PDM_pointer_array_part_get(gnum_elt, &
  !                                   i-1,      &
  !                                   ln_to_gn)
  !   deallocate(ln_to_gn)
  !   call PDM_pointer_array_part_get(part_data, &
  !                                   i-1,       &
  !                                   data)
  !   deallocate(data)
  ! enddo
  call PDM_pointer_array_free(gnum_elt)
  call PDM_pointer_array_free(part_data)
  deallocate(n_elt, data_distrib_idx)
  call PDM_fortran_free_c(c_loc(block_data))

  if (i_rank .eq. 0) then
    write(*, *) "-- End"
  end if

  call MPI_finalize(code)

contains

function randint(low, high) result (n)
  implicit none
  integer, intent(in) :: low, high
  integer             :: n
  real                :: r

  call random_number(r)
  n = low + floor((high - low) * r)

end function randint

end program testf

