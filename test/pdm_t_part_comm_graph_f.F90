
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
  !  The numbers correspond to the local IDs of the vertices and edges.
  !  within their respective parts.
  !
  !               -------- Rank 0 -------     ----- Rank 1 ----
  !
  !        {            1 ————— 2 ————— 3 · · 1 ————— 2 ————— 3
  !        {            |   1   |   2   |     |   1   |   2   |
  !        {            |3     4|      5|     |3     4|      5|
  ! part 1 {            |   6   |   7   |     |   6   |   7   |
  !        {            4 ————— 5 ————— 6 · · 4 ————— 5 ————— 6
  !        {          ·       · |       |     |       |       |
  !        {        ·       ·   |8     9|     |8      |9    10|
  !               1 ————— 2     |  10   |     |  11   |  12   |
  !               |   1   |     7 ————— 8 · · 7 ————— 8 ————— 9
  !        {      |2     3|   · |       | ·   ·       ·       ·
  !        {      |   4   | ·   |11   12|   · ·       ·       ·
  !        {      3 ————— 4     |  13   |     1 ————— 2 ————— 3
  !        {      |       |     9 ————— 10    |   1   |   2   |
  !        {      |5     6|   ·       ·   ·   |3     4|      5|
  ! part 2 {      |   7   | ·       ·       · |   6   |   7   |
  !        {      5 ————— 6 ————— 7 · · · · · 4 ————— 5 ————— 6
  !        {      |       |   8   |           |       |       |
  !        {      |9      |10   11|           |8      |9    10|
  !        {      |  12   |  13   |           |  11   |  12   |
  !        {      8 ————— 9 ————— 10· · · · · 7 ————— 8 ————— 9
  !
  ! ==============================================================


  !-----------------------------------------------
  ! Dummy derived type to hold partition data
  type my_part_t

  integer(pdm_l_num_s), pointer :: vtx_graph(:)
  real(8),              pointer :: data(:)
  real(8),              pointer :: send_data(:)

  integer(pdm_l_num_s), pointer :: edge_vtx_idx(:)
  integer(pdm_l_num_s), pointer :: edge_vtx(:)

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
  integer(pdm_l_num_s),      pointer     :: pn_vtx(:)
  integer(pdm_l_num_s),      pointer     :: pn_vtx_graph(:)
  type(pdm_pointer_array_t), pointer     :: pvtx_graph
  type(my_part_t),           allocatable :: parts(:)
  integer                                :: i_part

  type(c_ptr)                            :: pgc_vtx
  integer(pdm_l_num_s),      pointer     :: is_owner(:)

  integer                                :: stride
  type(pdm_pointer_array_t), pointer     :: send_stride
  type(pdm_pointer_array_t), pointer     :: send_data
  type(pdm_pointer_array_t), pointer     :: recv_stride
  type(pdm_pointer_array_t), pointer     :: recv_data
  real(8),                   pointer     :: data(:)
  real(8)                                :: expected, diff
  integer                                :: request
  type(pdm_pointer_array_t), pointer     :: part_data
  integer                                :: i, j, i_vtx, i_exch

  integer(pdm_l_num_s),      pointer     :: pn_edge(:)
  type(pdm_pointer_array_t), pointer     :: pedge_vtx_idx
  type(pdm_pointer_array_t), pointer     :: pedge_vtx
  type(c_ptr)                            :: pcg_edge
  integer                                :: n_edge_graph
  integer(pdm_l_num_s),      pointer     :: edge_graph(:)
  !--------------------------------------------------------------

  verbose = .false.

  nullify(pn_vtx,       &
          pn_vtx_graph, &
          pvtx_graph)

  pgc_vtx = C_NULL_PTR


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
  allocate(pn_vtx(n_part),       &
           pn_vtx_graph(n_part), &
           parts(n_part))

  if (i_rank == 0) then
    ! Rank 0
    !   part 1
    pn_vtx(1)       = 10
    pn_vtx_graph(1) = 10

    allocate(parts(1)%vtx_graph(4 * pn_vtx_graph(1)))
    parts(1)%vtx_graph = [3, 1, 1, 1, &
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
    pn_vtx(2)       = 10
    pn_vtx_graph(2) = 7

    allocate(parts(2)%vtx_graph(4 * pn_vtx_graph(2)))
    parts(2)%vtx_graph = [1, 0, 1, 4,  &
                          2, 0, 1, 5,  &
                          4, 0, 1, 7,  &
                          6, 0, 1, 9,  &
                          7, 0, 1, 10, &
                          7, 1, 2, 4,  &
                         10, 1, 2, 7]

  else
    ! Rank 1
    !   part 1
    pn_vtx(1)       = 9
    pn_vtx_graph(1) = 6

    allocate(parts(1)%vtx_graph(4 * pn_vtx_graph(1)))
    parts(1)%vtx_graph = [1, 0, 1, 3, &
                          4, 0, 1, 6, &
                          7, 0, 1, 8, &
                          7, 1, 2, 1, &
                          8, 1, 2, 2, &
                          9, 1, 2, 3]

    !   part 2
    pn_vtx(2)       = 9
    pn_vtx_graph(2) = 7
    allocate(parts(2)%vtx_graph(4 * pn_vtx_graph(2)))
    parts(2)%vtx_graph = [1, 0, 1, 8,  &
                          1, 1, 1, 7,  &
                          2, 1, 1, 8,  &
                          3, 1, 1, 9,  &
                          4, 0, 2, 7,  &
                          4, 0, 1, 10, &
                          7, 0, 2, 10]
  endif

  ! Set the pointer_array
  call pdm_pointer_array_create(pvtx_graph,   &
                                n_part,       &
                                PDM_TYPE_INT)

  do i_part = 1, n_part
    call pdm_pointer_array_part_set(pvtx_graph,              &
                                    i_part-1,                &
                                    parts(i_part)%vtx_graph)
  enddo
  !----------------------------------------

  !----------------------------------------
  ! Create a Part Comm Graph instance
  call pdm_part_comm_graph_create(pgc_vtx,            &
                                  n_part,             &
                                  pn_vtx_graph,       &
                                  pvtx_graph,         &
                                  PDM_OWNERSHIP_USER, &
                                  comm)
  !----------------------------------------


  !----------------------------------------
  ! Inspect owners/ghosts
  if (verbose) then
    write (funit, *) "OWNER"
  endif

  do i_part = 1, n_part
    call pdm_part_comm_graph_owner_get(pgc_vtx,  &
                                       i_part-1, &
                                       is_owner)

    if (verbose) then
      write (funit, *) "part ", i_part
      do i = 1, pn_vtx_graph(i_part)
        write (funit, *) "  ", parts(i_part)%vtx_graph(4*(i-1)+1), " : ", is_owner(i)
      enddo
    endif

    if (.not.check_owners(i_rank, i_part, is_owner)) then
      print *, "Error, incorrect `is_owner` for rank", i_rank, " part", i_part
      stop
    endif
  enddo
  !----------------------------------------


  !----------------------------------------
  ! Exchange constant-stride data
  stride = 2
  nullify(send_stride, &
          send_data,   &
          recv_stride, &
          recv_data)

  call pdm_pointer_array_create(send_data,      &
                                n_part,         &
                                PDM_TYPE_DOUBLE)

  do i_part = 1, n_part

    ! Data for all entities in partition
    allocate(parts(i_part)%data(pn_vtx(i_part) * stride))
    do i_vtx = 1, pn_vtx(i_part)
      do j = 1, stride
        parts(i_part)%data((i_vtx-1)*stride+j) = 100.0d0*i_vtx  + &
                                                  10.0d0*i_rank + &
                                                   1.0d0*i_part + &
                                                   0.1d0*(j-1)
      enddo
    enddo


    ! Allocate data for entities on partition boundary
    allocate(parts(i_part)%send_data(pn_vtx_graph(i_part) * stride))
    call pdm_pointer_array_part_set(send_data,               &
                                    i_part-1,                &
                                    parts(i_part)%send_data)
  enddo

  ! Four rounds of echanges (blocking, then non-blocking, then persistent (twice))
  do i_exch = 1, 4

    if (verbose) then
      write (funit, *) "SEND"
    endif

    ! Prepare send_data (multiply part_data) by i_exch
    do i_part = 1, n_part
      if (verbose) then
        write (funit, *) "part ", i_part
      endif

      call pdm_pointer_array_part_get(send_data, &
                                      i_part-1,  &
                                      data)
       do i = 1, pn_vtx_graph(i_part)
        i_vtx = parts(i_part)%vtx_graph(4*(i-1)+1)
        do j = 1, stride
          data((i-1)*stride+j) = i_exch * parts(i_part)%data((i_vtx-1)*stride+j)
        enddo

        if (verbose) then
          write (funit, *) parts(i_part)%vtx_graph(4*(i-1)+1), i_rank, i_part
          write (funit, *) "     ", data((i-1)*stride+1:i*stride)
        endif
      enddo
    enddo

    if (i_exch == 1) then
      if (verbose) then
        write (funit, *) "Blocking exchange"
      endif
      ! Blocking exchange
      call pdm_part_comm_graph_exch(pgc_vtx,                   &
                                    PDM_STRIDE_CST_INTERLACED, &
                                    stride,                    &
                                    send_stride,               &
                                    send_data,                 &
                                    recv_stride,               &
                                    recv_data)

    else if (i_exch == 2) then
      if (verbose) then
        write (funit, *) "Non-blocking exchange"
      endif
      ! Initiate non-blocking exchange
      call pdm_part_comm_graph_iexch(pgc_vtx,                   &
                                     PDM_MPI_COMM_KIND_P2P,     &
                                     PDM_STRIDE_CST_INTERLACED, &
                                     stride,                    &
                                     send_stride,               &
                                     send_data,                 &
                                     recv_stride,               &
                                     recv_data,                 &
                                     request)

      ! Do stuff here to cover MPI communications...

      ! Wait for exchange to finish
      call pdm_part_comm_graph_exch_wait(pgc_vtx, request)

    else
      if (verbose) then
        write (funit, *) "Persistent exchange"
      endif

      if (i_exch == 3) then
        ! Prepare persistent exchange
        call pdm_part_comm_graph_exch_init(pgc_vtx,                   &
                                           PDM_MPI_COMM_KIND_P2P,     &
                                           PDM_STRIDE_CST_INTERLACED, &
                                           stride,                    &
                                           send_stride,               &
                                           send_data,                 &
                                           recv_stride,               &
                                           recv_data,                 &
                                           request)
      endif

      ! We can use the same persistent channel multiple times (with the same send/recv buffers)

      ! Start exchange
      call pdm_part_comm_graph_exch_start(pgc_vtx, request)

      ! Do stuff here to cover MPI communications...

      ! Wait for exchange to finish
      call pdm_part_comm_graph_exch_wait(pgc_vtx, request)

      if (i_exch == 4) then
        ! Free the persistent exchange
        call pdm_part_comm_graph_exch_free(pgc_vtx, request)
      endif
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

      do i = 1, pn_vtx_graph(i_part)

        if (verbose) then
          write (funit, *) parts(i_part)%vtx_graph(4*(i-1)+4), parts(i_part)%vtx_graph(4*(i-1)+2), parts(i_part)%vtx_graph(4*(i-1)+3)
          write (funit, *) "     ", data((i-1)*stride+1:i*stride)
        endif

        do j = 1, stride
          expected = 100.0d0*parts(i_part)%vtx_graph(4*(i-1)+4) + &
                      10.0d0*parts(i_part)%vtx_graph(4*(i-1)+2) + &
                       1.0d0*parts(i_part)%vtx_graph(4*(i-1)+3) + &
                       0.1d0*(j-1)
          expected = expected * i_exch

          diff = abs(data((i-1)*stride+j) - expected)
          if (diff > 1.e-9) then
            print *, "Error i_exch", i_exch, ": (", i_rank, i_part, i, j, ") expected ", expected, " but received ", data((i-1)*stride+j), " from", parts(i_part)%vtx_graph(4*(i-1)+2:4*i)
            stop
          endif

        enddo
      enddo
    enddo

    if (i_exch /= 3) then
      call pdm_pointer_array_free(recv_data)
    endif
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

  call pdm_part_comm_graph_all_reduce(pgc_vtx,   &
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

      do i = 1, pn_vtx_graph(i_part)
        i_vtx = parts(i_part)%vtx_graph(4*(i-1)+1)
        write (funit, *) i_vtx
        write (funit, *) "     ", data((i_vtx-1)*stride+1:i_vtx*stride)
      enddo
    enddo
  endif
  !----------------------------------------



  !----------------------------------------
  ! Create edge->vtx connectivity
  allocate(pn_edge(n_part))
  if (i_rank == 0) then
    ! Rank 0
    !   part 1
    pn_edge(1) = 13
    allocate(parts(1)%edge_vtx(2 * pn_edge(1)))
    parts(1)%edge_vtx = [1, 2,  &
                         2, 3,  &
                         1, 4,  &
                         2, 5,  &
                         3, 6,  &
                         4, 5,  &
                         5, 6,  &
                         5, 7,  &
                         6, 8,  &
                         7, 8,  &
                         7, 9,  &
                         8, 10, &
                         9, 10]

    !   part 2
    pn_edge(2) = 13
    allocate(parts(2)%edge_vtx(2 * pn_edge(2)))
    parts(2)%edge_vtx = [1, 2,  &
                         1, 3,  &
                         2, 4,  &
                         3, 4,  &
                         3, 5,  &
                         4, 6,  &
                         5, 6,  &
                         6, 7,  &
                         5, 8,  &
                         6, 9,  &
                         7, 10, &
                         8, 9,  &
                         9, 10]

  else
    ! Rank 1
    !   parts 1 and 2 are identical
    do i_part = 1, n_part
      pn_edge(i_part) = 12
      allocate(parts(i_part)%edge_vtx(2 * pn_edge(i_part)))
      parts(i_part)%edge_vtx = [1, 2, &
                                2, 3, &
                                1, 4, &
                                2, 5, &
                                3, 6, &
                                4, 5, &
                                5, 6, &
                                4, 7, &
                                5, 8, &
                                6, 9, &
                                7, 8, &
                                8, 9]
    enddo
  endif

  do i_part = 1, n_part
    allocate(parts(i_part)%edge_vtx_idx(pn_edge(i_part)+1))
    do i = 0, pn_edge(i_part)
      parts(i_part)%edge_vtx_idx(i+1) = 2*i
    enddo
  enddo


  nullify(pedge_vtx_idx, &
          pedge_vtx)
  call pdm_pointer_array_create(pedge_vtx_idx, &
                                n_part,        &
                                PDM_TYPE_INT)

  call pdm_pointer_array_create(pedge_vtx,    &
                                n_part,       &
                                PDM_TYPE_INT)

  do i_part = 1, n_part
    call pdm_pointer_array_part_set(pedge_vtx_idx,              &
                                    i_part-1,                   &
                                    parts(i_part)%edge_vtx_idx)

    call pdm_pointer_array_part_set(pedge_vtx,              &
                                    i_part-1,               &
                                    parts(i_part)%edge_vtx)
  enddo
  !----------------------------------------



  !----------------------------------------
  ! Create edge Part Comm Graph for the vtx Part Comm Graph
  call pdm_part_comm_graph_entity1_to_part_comm_graph_entity2(pgc_vtx,       &
                                                              pn_vtx,        &
                                                              pn_edge,       &
                                                              pedge_vtx_idx, &
                                                              pedge_vtx,     &
                                                              pcg_edge)

  if (verbose) then
    write (funit, *) "EDGE PCG"
  endif

  do i_part = 1, n_part
    n_edge_graph = pdm_part_comm_graph_n_entity_get(pcg_edge, i_part-1)

    call pdm_part_comm_graph_entity_graph_get(pcg_edge,           &
                                              i_part-1,           &
                                              edge_graph,         &
                                              PDM_OWNERSHIP_KEEP)

    if (verbose) then
      write (funit, *) "part", i_part
      do i = 1, n_edge_graph
        write (funit, *) "    ", edge_graph(4*(i-1)+1:4*i)
      enddo
    endif

    if (.not.check_edge_graph(i_rank, i_part, edge_graph)) then
      print *, "Error, incorrect `edge_graph` for rank", i_rank, " part", i_part
      stop
    endif
  enddo
  !----------------------------------------



  !----------------------------------------
  ! Free memory
  call pdm_pointer_array_free(part_data)
  call pdm_pointer_array_free(send_data)
  call pdm_part_comm_graph_free(pgc_vtx)
  call pdm_part_comm_graph_free(pcg_edge)

  do i_part = 1, n_part
    deallocate(parts(i_part)%vtx_graph)
    deallocate(parts(i_part)%data)
    deallocate(parts(i_part)%send_data)
    deallocate(parts(i_part)%edge_vtx_idx)
    deallocate(parts(i_part)%edge_vtx)
  enddo
  call pdm_pointer_array_free(pvtx_graph)
  deallocate(pn_vtx,       &
             pn_vtx_graph, &
             parts,        &
             pn_edge)
  call pdm_pointer_array_free(pedge_vtx_idx)
  call pdm_pointer_array_free(pedge_vtx)
  !----------------------------------------

  if (i_rank == 0) then
    print *, "The End :)"
  endif

  call mpi_finalize(err)



  contains


  function check_array_eq(a, b, n) &
  result (equal)

    implicit none
    integer(pdm_l_num_s) :: a(:)
    integer(pdm_l_num_s) :: b(:)
    integer              :: n
    logical              :: equal
    integer              :: i

    equal = .true.
    do i = 1, n
      if (a(i) /= b(i)) then
        equal = .false.
        print *, ""
        return
      endif
    enddo

  end function check_array_eq



  function check_owners(i_rank, i_part, is_owner) &
  result (ok)

    implicit none

    integer              :: i_rank
    integer              :: i_part
    integer(pdm_l_num_s) :: is_owner(:)
    logical              :: ok
    integer(pdm_l_num_s) :: exp_is_owner_01(10)
    integer(pdm_l_num_s) :: exp_is_owner_02(7)
    integer(pdm_l_num_s) :: exp_is_owner_11(6)
    integer(pdm_l_num_s) :: exp_is_owner_12(7)

    exp_is_owner_01 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    exp_is_owner_02 = [0, 0, 0, 0, 0, 0, 1]
    exp_is_owner_11 = [0, 0, 0, 0, 1, 1]
    exp_is_owner_12 = [0, 0, 0, 0, 0, 0, 0]

    if (i_rank == 0) then
      if (i_part == 1) then
        ok = check_array_eq(is_owner, exp_is_owner_01, size(exp_is_owner_01))
      else
        ok = check_array_eq(is_owner, exp_is_owner_02, size(exp_is_owner_02))
      endif
    else
      if (i_part == 1) then
        ok = check_array_eq(is_owner, exp_is_owner_11, size(exp_is_owner_11))
      else
        ok = check_array_eq(is_owner, exp_is_owner_12, size(exp_is_owner_12))
      endif
    endif

  end function check_owners



  function check_edge_graph(i_rank, i_part, edge_graph) &
  result (ok)

    implicit none

    integer              :: i_rank
    integer              :: i_part
    integer(pdm_l_num_s) :: edge_graph(:)
    logical              :: ok
    integer(pdm_l_num_s) :: exp_edge_graph_01(28)
    integer(pdm_l_num_s) :: exp_edge_graph_02(20)
    integer(pdm_l_num_s) :: exp_edge_graph_11(16)
    integer(pdm_l_num_s) :: exp_edge_graph_12(16)

    exp_edge_graph_01 = [ 6, 0, 2, 1, &
                          8, 0, 2, 3, &
                         11, 0, 2, 6, &
                         13, 0, 2, 8, &
                          5, 1, 1, 3, &
                          9, 1, 1, 8, &
                         12, 1, 2, 3]

    exp_edge_graph_02 = [ 1, 0, 1,  6, &
                          3, 0, 1,  8, &
                          6, 0, 1, 11, &
                          8, 0, 1, 13, &
                         11, 1, 2,  8]

    exp_edge_graph_11 = [ 3, 0, 1, 5, &
                          8, 0, 1, 9, &
                         11, 1, 2, 1, &
                         12, 1, 2, 2]

    exp_edge_graph_12 = [3, 0, 1, 12, &
                         8, 0, 2, 11, &
                         1, 1, 1, 11, &
                         2, 1, 1, 12]

    if (i_rank == 0) then
      if (i_part == 1) then
        ok = check_array_eq(edge_graph, exp_edge_graph_01, size(exp_edge_graph_01))
      else
        ok = check_array_eq(edge_graph, exp_edge_graph_02, size(exp_edge_graph_02))
      endif
    else
      if (i_part == 1) then
        ok = check_array_eq(edge_graph, exp_edge_graph_11, size(exp_edge_graph_11))
      else
        ok = check_array_eq(edge_graph, exp_edge_graph_12, size(exp_edge_graph_12))
      endif
    endif

  end function check_edge_graph


end program test_part_comm_graph
