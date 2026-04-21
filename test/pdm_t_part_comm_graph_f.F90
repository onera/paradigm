
#include "pdm_configf.h"

program test_part_comm_graph

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif

  use pdm_pointer_array
  use pdm_part_comm_graph
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif


  ! ==============================================================
  !  This test emulates the situation illustrated below :
  !  A mesh is partitioned across 2 MPI ranks, each with 2 parts.
  !  We are interested in the inter-partition connections between
  !  the vertices of this mesh, represented by the dotted lines.
  !  The numbers correspond to the local IDs of the vertices
  !  within their respective parts.
  !
  !               -------- Rank 0 -------     ----- Rank 1 ----
  !
  !        {            1 ————— 2 ————— 3 · · 1 ————— 2 ————— 3
  !        {            |       |       |     |       |       |
  !        {            |       |       |     |       |       |
  ! part 1 {            |       |       |     |       |       |
  !        {            4 ————— 5 ————— 6 · · 4 ————— 5 ————— 6
  !        {          ·       · |       |     |       |       |
  !        {        ·       ·   |       |     |       |       |
  !               1 ————— 2     |       |     |       |       |
  !               |       |     7 ————— 8 · · 7 ————— 8 ————— 9
  !        {      |       |   · |       | ·   ·       ·       ·
  !        {      |       | ·   |       |   · ·       ·       ·
  !        {      3 ————— 4     |       |     1 ————— 2 ————— 3
  !        {      |       |     9 ————— 10    |       |       |
  !        {      |       |   ·       ·   ·   |       |       |
  ! part 2 {      |       | ·       ·       · |       |       |
  !        {      5 ————— 6 ————— 7 · · · · · 4 ————— 5 ————— 6
  !        {      |       |       |           |       |       |
  !        {      |       |       |           |       |       |
  !        {      |       |       |           |       |       |
  !        {      8 ————— 9 ————— 10· · · · · 7 ————— 8 ————— 9
  !
  ! ==============================================================


  !-----------------------------------------------
  ! Dummy derived type to hold partition data
  type my_part_t

  integer(pdm_l_num_s), pointer :: entity_graph(:)
  real(8),              pointer :: data(:)
  real(8),              pointer :: send_data(:)

  end type my_part_t
  !-----------------------------------------------


  !--------------------------------------------------------------
  ! Variables
  integer, parameter                     :: comm = MPI_COMM_WORLD
  integer                                :: i_rank, n_rank, err
  integer                                :: i_arg
  character(len=99)                      :: arg

  logical                                :: verbose
  integer                                :: funit
  character                              :: strnum


  integer, parameter                     :: n_part = 2
  integer(pdm_l_num_s),      pointer     :: pn_entity(:)
  integer(pdm_l_num_s),      pointer     :: pn_entity_graph(:)
  type(pdm_pointer_array_t), pointer     :: pentity_graph
  type(my_part_t),           allocatable :: parts(:)
  integer                                :: i_part

  type(c_ptr)                            :: pcg

  integer                                :: stride
  type(pdm_pointer_array_t), pointer     :: send_stride
  type(pdm_pointer_array_t), pointer     :: send_data
  type(pdm_pointer_array_t), pointer     :: recv_stride
  type(pdm_pointer_array_t), pointer     :: recv_data
  real(8),                   pointer     :: data(:)
  real(8)                                :: expected, diff
  integer                                :: request
  type(pdm_pointer_array_t), pointer     :: part_data
  integer                                :: i, j, i_entity, i_try
  !--------------------------------------------------------------

  verbose = .false.

  nullify(pn_entity,       &
          pn_entity_graph, &
          pentity_graph)

  pcg = C_NULL_PTR


  !----------------------------------------
  ! Parse command line arguments
  i_arg = 1
  do while (i_arg <= command_argument_count())
    call get_command_argument(i_arg, arg)
    select case(arg)
      case ("-v")
        verbose = .true.
    endselect
    i_arg = i_arg + 1
  enddo
  !----------------------------------------


  !----------------------------------------
  ! Initialize MPI
  call mpi_init(err)
  call mpi_comm_rank(comm, i_rank, err)
  call mpi_comm_size(comm, n_rank, err)

  if (n_rank /= 2) then
    print *, "This test is supposed to run with exactly 2 MPI ranks"
    call mpi_finalize(err)
    stop
  endif

  if (verbose) then
    write (strnum, '(i1)') i_rank
    open(unit=funit, file="part_comm_graph_f_"//strnum//".log", action='write')
  endif
  !----------------------------------------


  !----------------------------------------
  ! Define the inter-partition communication graph
  allocate(pn_entity(n_part),       &
           pn_entity_graph(n_part), &
           parts(n_part))

  if (i_rank == 0) then
    ! Rank 0
    !   part 1
    pn_entity(1)       = 10
    pn_entity_graph(1) = 10

    allocate(parts(1)%entity_graph(4 * pn_entity_graph(1)))
    parts(1)%entity_graph = [3, 1, 1, 1, &
                             4, 0, 2, 1, &
                             5, 0, 2, 2, &
                             6, 1, 1, 4, &
                             7, 0, 2, 4, &
                             8, 1, 1, 7, &
                             8, 1, 2, 1, &
                             9, 0, 2, 6, &
                            10, 0, 2, 7, &
                            10, 1, 2, 4]

    !   part 2
    pn_entity(2)       = 10
    pn_entity_graph(2) = 7

    allocate(parts(2)%entity_graph(4 * pn_entity_graph(2)))
    parts(2)%entity_graph = [1, 0, 1, 4,  &
                             2, 0, 1, 5,  &
                             4, 0, 1, 7,  &
                             6, 0, 1, 9,  &
                             7, 0, 1, 10, &
                             7, 1, 2, 4,  &
                            10, 1, 2, 7]

  else
    ! Rank 1
    !   part 1
    pn_entity(1)       = 9
    pn_entity_graph(1) = 6

    allocate(parts(1)%entity_graph(4 * pn_entity_graph(1)))
    parts(1)%entity_graph = [1, 0, 1, 3, &
                             4, 0, 1, 6, &
                             7, 0, 1, 8, &
                             7, 1, 2, 1, &
                             8, 1, 2, 2, &
                             9, 1, 2, 3]

    !   part 2
    pn_entity(2)       = 9
    pn_entity_graph(2) = 7
    allocate(parts(2)%entity_graph(4 * pn_entity_graph(2)))
    parts(2)%entity_graph = [1, 0, 1, 8,  &
                             1, 1, 1, 7,  &
                             2, 1, 1, 8,  &
                             3, 1, 1, 9,  &
                             4, 0, 2, 7,  &
                             4, 0, 1, 10, &
                             7, 0, 2, 10]
  endif

  ! Set the pointer_array
  call pdm_pointer_array_create(pentity_graph, &
                                n_part,        &
                                PDM_TYPE_INT)

  do i_part = 1, n_part
    call pdm_pointer_array_part_set(pentity_graph,              &
                                    i_part-1,                   &
                                    parts(i_part)%entity_graph)
  enddo
  !----------------------------------------

  !----------------------------------------
  ! Create a Part Comm Graph instance
  call pdm_part_comm_graph_create(pcg,                &
                                  n_part,             &
                                  pn_entity_graph,    &
                                  pentity_graph,      &
                                  PDM_OWNERSHIP_USER, &
                                  comm)
  !----------------------------------------

  !----------------------------------------
  ! Exchange data
  stride = 2
  nullify(send_stride, &
          send_data,   &
          recv_stride, &
          recv_data)

  call pdm_pointer_array_create(send_data,      &
                                n_part,         &
                                PDM_TYPE_DOUBLE)

  if (verbose) then
    write (funit, *) "SEND"
  endif

  do i_part = 1, n_part

    ! Data for all entities in partition
    allocate(parts(i_part)%data(pn_entity(i_part) * stride))
    do i_entity = 1, pn_entity(i_part)
      do j = 1, stride
        parts(i_part)%data((i_entity-1)*stride+j) = 100.0d0*i_entity + &
                                                     10.0d0*i_rank   + &
                                                      1.0d0*i_part   + &
                                                      0.1d0*(j-1)
      enddo
    enddo


    ! Data for entities on partition boundary
    if (verbose) then
      write (funit, *) "part ", i_part
    endif

    allocate(parts(i_part)%send_data(pn_entity_graph(i_part) * stride))
    do i = 1, pn_entity_graph(i_part)
      i_entity = parts(i_part)%entity_graph(4*(i-1)+1)
      do j = 1, stride
        parts(i_part)%send_data((i-1)*stride+j) = parts(i_part)%data((i_entity-1)*stride+j)
      enddo

      if (verbose) then
        write (funit, *) i_entity, i_rank, i_part
        write (funit, *) "     ", parts(i_part)%send_data((i-1)*stride+1:i*stride)
      endif

    enddo

    call pdm_pointer_array_part_set(send_data,               &
                                    i_part-1,                &
                                    parts(i_part)%send_data)
  enddo

  ! Two rounds of echanges (blocking, then non-blocking)
  do i_try = 1, 2

    if (i_try == 1) then
      ! Blocking exchange
      call pdm_part_comm_graph_exch(pcg,                       &
                                    PDM_STRIDE_CST_INTERLACED, &
                                    stride,                    &
                                    send_stride,               &
                                    send_data,                 &
                                    recv_stride,               &
                                    recv_data)
    else
      ! Initiate non-blocking exchange
      call pdm_part_comm_graph_iexch(pcg,                       &
                                     PDM_MPI_COMM_KIND_P2P,     &
                                     PDM_STRIDE_CST_INTERLACED, &
                                     stride,                    &
                                     send_stride,               &
                                     send_data,                 &
                                     recv_stride,               &
                                     recv_data,                 &
                                     request)

      ! Wait for exchange to finish
      call pdm_part_comm_graph_exch_wait(pcg, request)
    endif


    ! Check the received data
    if (verbose) then
      write (funit, *) "RECV"
    endif

    do i_part = 1, n_part

      if (verbose) then
        write (funit, *) "part ", i_part
      endif

      call pdm_pointer_array_part_get(recv_data, &
                                      i_part-1,  &
                                      data)

      do i = 1, pn_entity_graph(i_part)

        if (verbose) then
          write (funit, *) parts(i_part)%entity_graph(4*(i-1)+4), parts(i_part)%entity_graph(4*(i-1)+2), parts(i_part)%entity_graph(4*(i-1)+3)
          write (funit, *) "     ", data((i-1)*stride+1:i*stride)
        endif

        do j = 1, stride
          expected = 100.0d0*parts(i_part)%entity_graph(4*(i-1)+4) + &
                      10.0d0*parts(i_part)%entity_graph(4*(i-1)+2) + &
                       1.0d0*parts(i_part)%entity_graph(4*(i-1)+3) + &
                       0.1d0*(j-1)

          diff = abs(data((i-1)*stride+j) - expected)
          if (diff > 1.e-9) then
            print *, "Error : (", i_rank, i_part, i, j, ") expected ", expected, " but received ", data((i-1)*stride+j), " from", parts(i_part)%entity_graph(4*(i-1)+2:4*i)
            stop
          endif

        enddo
      enddo
    enddo

    call pdm_pointer_array_free(recv_data)
  enddo
  !----------------------------------------


  !----------------------------------------
  ! Reduction
  nullify(part_data)
  call pdm_pointer_array_create(part_data,      &
                                n_part,         &
                                PDM_TYPE_DOUBLE)
  do i_part = 1, n_part
    call pdm_pointer_array_part_set(part_data,          &
                                    i_part-1,           &
                                    parts(i_part)%data)
  enddo

  call pdm_part_comm_graph_all_reduce(pcg,       &
                                      stride,    &
                                      MPI_MAX,   &
                                      part_data)

  if (verbose) then
    write (funit, *) "REDUCE"
    do i_part = 1, n_part
      write (funit, *) "part ", i_part

      call pdm_pointer_array_part_get(part_data, &
                                      i_part-1,  &
                                      data)

      do i = 1, pn_entity_graph(i_part)
        i_entity = parts(i_part)%entity_graph(4*(i-1)+1)
        write (funit, *) i_entity
        write (funit, *) "     ", data((i_entity-1)*stride+1:i_entity*stride)
      enddo
    enddo
  endif
  !----------------------------------------


  !----------------------------------------
  ! Free memory
  call pdm_pointer_array_free(part_data)
  call pdm_pointer_array_free(send_data)
  call pdm_part_comm_graph_free(pcg)

  do i_part = 1, n_part
    deallocate(parts(i_part)%entity_graph)
    deallocate(parts(i_part)%data)
    deallocate(parts(i_part)%send_data)
  enddo
  call pdm_pointer_array_free(pentity_graph)
  deallocate(pn_entity,       &
             pn_entity_graph, &
             parts)
  !----------------------------------------

  if (i_rank == 0) then
    print *, "The End :)"
  endif

  call mpi_finalize(err)


end program test_part_comm_graph
