
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
    integer(pdm_l_num_s), pointer :: send_stride(:)

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

  type(c_ptr)                            :: pcg_vtx
  integer(pdm_l_num_s),      pointer     :: is_owner(:)

  integer                                :: i

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

  pcg_vtx = C_NULL_PTR


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
  call pdm_part_comm_graph_create(pcg_vtx,            &
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
    call pdm_part_comm_graph_owner_get(pcg_vtx,  &
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
  call exchanges_constant_stride(pcg_vtx, parts, pn_vtx_graph)
  !----------------------------------------

  !----------------------------------------
  ! Exchange variable-stride data
  call exchanges_variable_stride(pcg_vtx, parts, pn_vtx_graph)
  !----------------------------------------


  !----------------------------------------
  ! Reduction
  call allreduce(pcg_vtx, parts)
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
  call pdm_part_comm_graph_entity1_to_part_comm_graph_entity2(pcg_vtx,       &
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


  if (verbose) then
    close(funit)
  endif



  !----------------------------------------
  ! Free memory
  call pdm_part_comm_graph_free(pcg_vtx)
  call pdm_part_comm_graph_free(pcg_edge)

  do i_part = 1, n_part
    deallocate(parts(i_part)%vtx_graph)
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

    !--------------------------------------------------------------
    integer(pdm_l_num_s) :: a(:)
    integer(pdm_l_num_s) :: b(:)
    integer              :: n
    logical              :: equal
    integer              :: i
    !--------------------------------------------------------------

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

    !--------------------------------------------------------------
    integer              :: i_rank
    integer              :: i_part
    integer(pdm_l_num_s) :: is_owner(:)
    logical              :: ok
    integer(pdm_l_num_s) :: exp_is_owner_01(10)
    integer(pdm_l_num_s) :: exp_is_owner_02(7)
    integer(pdm_l_num_s) :: exp_is_owner_11(6)
    integer(pdm_l_num_s) :: exp_is_owner_12(7)
    !--------------------------------------------------------------

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

    !--------------------------------------------------------------
    integer              :: i_rank
    integer              :: i_part
    integer(pdm_l_num_s) :: edge_graph(:)
    logical              :: ok
    integer(pdm_l_num_s) :: exp_edge_graph_01(28)
    integer(pdm_l_num_s) :: exp_edge_graph_02(20)
    integer(pdm_l_num_s) :: exp_edge_graph_11(16)
    integer(pdm_l_num_s) :: exp_edge_graph_12(16)
    !--------------------------------------------------------------

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



  subroutine allreduce(pcg, parts)
    ! Test in-place reduction
    implicit none

    !--------------------------------------------------------------
    integer, parameter                  :: stride = 2

    type(c_ptr),          intent(in)    :: pcg
    type(my_part_t),      intent(inout) :: parts(n_part)

    type(pdm_pointer_array_t), pointer  :: part_data
    integer                             :: i_part, i_vtx, j, i_op
    integer                             :: op
    !--------------------------------------------------------------

    nullify(part_data)

    if (verbose) then
      write (funit, *) "ALL REDUCE"
    endif

    ! Allocate part data
    call pdm_pointer_array_create(part_data,      &
                                  n_part,         &
                                  PDM_TYPE_DOUBLE)

    do i_part = 1, n_part

      allocate(parts(i_part)%data(pn_vtx(i_part) * stride))

      call pdm_pointer_array_part_set(part_data,          &
                                      i_part-1,           &
                                      parts(i_part)%data)
    enddo


    ! Test all possible reduction operations (MIN, MAX, SUM)
    do i_op = 1, 3

      ! Reset part data
      if (verbose .and. (i_op == 1)) then
        write (funit, *) "BEFORE REDUCTION"
      endif

      do i_part = 1, n_part

        if (verbose .and. (i_op == 1)) then
          write (funit, *) "part ", i_part
        endif

        do i_vtx = 1, pn_vtx(i_part)
          do j = 1, stride
            parts(i_part)%data((i_vtx-1)*stride+j) = 100.0d0*i_vtx  + &
                                                      10.0d0*i_rank + &
                                                       1.0d0*i_part + &
                                                       0.1d0*(j-1)
          enddo

          if (verbose .and. (i_op == 1)) then
            write (funit, *) i_vtx, " :", parts(i_part)%data((i_vtx-1)*stride+1:i_vtx*stride)
          endif
        enddo
      enddo

      ! Perform reduction operation
      if (i_op == 1) then
        op = MPI_MIN
      else if (i_op == 2) then
        op = MPI_MAX
      else
        op = MPI_SUM
      endif

      if (verbose) then
        write (funit, *) "i_op", i_op
      endif

      call pdm_part_comm_graph_allreduce(pcg_vtx,   &
                                         stride,    &
                                         op,        &
                                         .false.,   &
                                         part_data)

      ! Check result of reduction operation
      if (verbose) then
        write (funit, *) "AFTER REDUCTION"
      endif

      do i_part = 1, n_part

        if (verbose) then
          write (funit, *) "part ", i_part
        endif

        do i_vtx = 1, pn_vtx(i_part)
          if (verbose) then
            write (funit, *) i_vtx, " :", parts(i_part)%data((i_vtx-1)*stride+1:i_vtx*stride)
          endif

        enddo
      enddo

    enddo

    do i_part = 1, n_part
      deallocate(parts(i_part)%data)
    enddo

    call pdm_pointer_array_free(part_data)

  end subroutine allreduce



  subroutine exchanges_constant_stride(pcg, parts, pn_vtx_graph)
    ! Test different modes of exchanges with constant stride
    implicit none

    !--------------------------------------------------------------
    integer, parameter                  :: cst_stride = 3

    type(c_ptr),          intent(in)    :: pcg
    type(my_part_t),      intent(inout) :: parts(n_part)
    integer(pdm_l_num_s), intent(in)    :: pn_vtx_graph(n_part)
    integer                             :: i_part, i_vtx, i, j, k

    integer                             :: send_size
    type(pdm_pointer_array_t), pointer  :: send_stride ! not used here
    type(pdm_pointer_array_t), pointer  :: send_data
    type(pdm_pointer_array_t), pointer  :: recv_stride ! not used here
    type(pdm_pointer_array_t), pointer  :: recv_data
    integer                             :: i_exch
    integer                             :: request

    real(8),                   pointer  :: data(:)
    real(8)                             :: expected, diff
    !--------------------------------------------------------------

    if (verbose) then
      write (funit, *) "EXCHANGE CONSTANT STRIDE"
    endif

    nullify(send_stride, &
            send_data,   &
            recv_stride, &
            recv_data)

    call pdm_pointer_array_create(send_data,       &
                                  n_part,          &
                                  PDM_TYPE_DOUBLE)

    ! Allocate send data
    do i_part = 1, n_part

      send_size = cst_stride * pn_vtx_graph(i_part)

      allocate(parts(i_part)%send_data(send_size))

      call pdm_pointer_array_part_set(send_data,               &
                                      i_part-1,                &
                                      parts(i_part)%send_data)

    enddo


    ! Four rounds of echanges (blocking, then non-blocking, then persistent (twice))
    do i_exch = 1, 4

      ! Prepare send_data (multiply part_data) by i_exch
      if (verbose) then
        write (funit, *) "SEND", i_exch
      endif

      do i_part = 1, n_part
        if (verbose) then
          write (funit, *) "part ", i_part
        endif

        call pdm_pointer_array_part_get(send_data, &
                                        i_part-1,  &
                                        data)
        k = 1
        do i = 1, pn_vtx_graph(i_part)
          i_vtx = parts(i_part)%vtx_graph(4*(i-1)+1)
          do j = 1, cst_stride
            data(k) = 100.0d0*i_vtx  + &
                       10.0d0*i_rank + &
                        1.0d0*i_part + &
                        0.1d0*(j-1)
            data(k) = i_exch * data(k)
            k = k + 1
          enddo

          if (verbose) then
            write (funit, *) parts(i_part)%vtx_graph(4*(i-1)+1), i_rank, i_part
            write (funit, *) "      ", data((i-1)*cst_stride+1:i*cst_stride)
          endif
        enddo
      enddo


      if (i_exch == 1) then
        if (verbose) then
          write (funit, *) "Blocking exchange"
        endif
        ! Blocking exchange
        call pdm_part_comm_graph_exch(pcg_vtx,                   &
                                      PDM_STRIDE_CST_INTERLACED, &
                                      cst_stride,                &
                                      send_stride,               &
                                      send_data,                 &
                                      recv_stride,               &
                                      recv_data)

      else if (i_exch == 2) then
        if (verbose) then
          write (funit, *) "Non-blocking exchange"
        endif
        ! Initiate non-blocking exchange
        call pdm_part_comm_graph_iexch(pcg_vtx,                   &
                                       PDM_MPI_COMM_KIND_P2P,     &
                                       PDM_STRIDE_CST_INTERLACED, &
                                       cst_stride,                &
                                       send_stride,               &
                                       send_data,                 &
                                       recv_stride,               &
                                       recv_data,                 &
                                       request)

        ! Do stuff here to cover MPI communications...

        ! Wait for exchange to finish
        call pdm_part_comm_graph_exch_wait(pcg_vtx, request)

      else
        if (verbose) then
          write (funit, *) "Persistent exchange"
        endif

        if (i_exch == 3) then
          ! Prepare persistent exchange
          call pdm_part_comm_graph_exch_init(pcg_vtx,                   &
                                             PDM_MPI_COMM_KIND_P2P,     &
                                             PDM_STRIDE_CST_INTERLACED, &
                                             cst_stride,                &
                                             send_stride,               &
                                             send_data,                 &
                                             recv_stride,               &
                                             recv_data,                 &
                                             request)
        endif

        ! We can use the same persistent channel multiple times (with the same send/recv buffers)

        ! Start exchange
        call pdm_part_comm_graph_exch_start(pcg_vtx, request)

        ! Do stuff here to cover MPI communications...

        ! Wait for exchange to finish
        call pdm_part_comm_graph_exch_wait(pcg_vtx, request)

        if (i_exch == 4) then
          ! Free the persistent exchange
          call pdm_part_comm_graph_exch_free(pcg_vtx, request)
        endif
      endif


      ! Check the received data
      if (verbose) then
        write (funit, *) "RECV", i_exch
      endif

      do i_part = 1, n_part

        if (verbose) then
          write (funit, *) "part ", i_part
        endif

        call pdm_pointer_array_part_get(recv_data, &
                                        i_part-1,  &
                                        data)

        k = 1
        do i = 1, pn_vtx_graph(i_part)

          if (verbose) then
            write (funit, *) "  from", parts(i_part)%vtx_graph(4*(i-1)+4), parts(i_part)%vtx_graph(4*(i-1)+2), parts(i_part)%vtx_graph(4*(i-1)+3)
            write (funit, *) "      ", data((i-1)*cst_stride+1:i*cst_stride)
          endif

          do j = 1, cst_stride
            expected = 100.0d0*parts(i_part)%vtx_graph(4*(i-1)+4) + &
                        10.0d0*parts(i_part)%vtx_graph(4*(i-1)+2) + &
                         1.0d0*parts(i_part)%vtx_graph(4*(i-1)+3) + &
                         0.1d0*(j-1)
            expected = expected * i_exch

            diff = abs(data(k) - expected)
            if (diff > 1.e-9) then
              print *, "Error i_exch", i_exch, ": (", i_rank, i_part, i, j, ") expected ", expected, " but received ", data(k), " from", parts(i_part)%vtx_graph(4*(i-1)+2:4*i)
              stop
            endif

            k = k + 1

          enddo
        enddo

      enddo

      if (i_exch /= 3) then
        call pdm_pointer_array_free(recv_data)
      endif

    enddo

    do i_part = 1, n_part
      deallocate(parts(i_part)%send_data)
    enddo

    call pdm_pointer_array_free(send_data)

  end subroutine exchanges_constant_stride




  subroutine exchanges_variable_stride(pcg, parts, pn_vtx_graph)
    ! Test different modes of exchanges with variable stride
    implicit none

    !--------------------------------------------------------------
    type(c_ptr),          intent(in)    :: pcg
    type(my_part_t),      intent(inout) :: parts(n_part)
    integer(pdm_l_num_s), intent(in)    :: pn_vtx_graph(n_part)
    integer                             :: i_part, i_vtx, i, j, k, idx

    integer                             :: send_size
    type(pdm_pointer_array_t), pointer  :: send_stride
    type(pdm_pointer_array_t), pointer  :: send_data
    type(pdm_pointer_array_t), pointer  :: recv_stride
    type(pdm_pointer_array_t), pointer  :: recv_data
    integer                             :: i_exch
    integer                             :: cst_stride = 0 ! not used here
    integer                             :: request

    integer(pdm_l_num_s),      pointer  :: stride(:)
    real(8),                   pointer  :: data(:)
    real(8)                             :: expected, diff
    !--------------------------------------------------------------

    if (verbose) then
      write (funit, *) "EXCHANGE VARIABLE STRIDE"
    endif

    nullify(send_stride, &
            send_data,   &
            recv_stride, &
            recv_data)

    call pdm_pointer_array_create(send_stride,  &
                                  n_part,       &
                                  PDM_TYPE_INT)

    call pdm_pointer_array_create(send_data,       &
                                  n_part,          &
                                  PDM_TYPE_DOUBLE)

    ! Create send_stride and allocate send data
    do i_part = 1, n_part

      allocate(parts(i_part)%send_stride(pn_vtx_graph(i_part)))

      send_size = 0
      do i = 1, pn_vtx_graph(i_part)
        i_vtx = parts(i_part)%vtx_graph(4*(i-1)+1)
        parts(i_part)%send_stride(i) = 1 + modulo(i_vtx, 4)
        send_size = send_size + parts(i_part)%send_stride(i)
      enddo

      allocate(parts(i_part)%send_data(send_size))

      call pdm_pointer_array_part_set(send_stride,               &
                                      i_part-1,                  &
                                      parts(i_part)%send_stride)

      call pdm_pointer_array_part_set(send_data,               &
                                      i_part-1,                &
                                      parts(i_part)%send_data)

    enddo


    ! Two rounds of echanges (blocking, then non-blocking)
    do i_exch = 1, 2

      ! Prepare send_data (multiply part_data) by i_exch
      if (verbose) then
        write (funit, *) "SEND", i_exch
      endif

      do i_part = 1, n_part
        if (verbose) then
          write (funit, *) "part ", i_part
        endif

        call pdm_pointer_array_part_get(send_data, &
                                        i_part-1,  &
                                        data)
        k = 1
        do i = 1, pn_vtx_graph(i_part)
          idx = k
          i_vtx = parts(i_part)%vtx_graph(4*(i-1)+1)
          do j = 1, parts(i_part)%send_stride(i)
            data(k) = 100.0d0*i_vtx  + &
                       10.0d0*i_rank + &
                        1.0d0*i_part + &
                        0.1d0*(j-1)
            data(k) = i_exch * data(k)
            k = k + 1
          enddo

          if (verbose) then
            write (funit, *) parts(i_part)%vtx_graph(4*(i-1)+1), i_rank, i_part
            write (funit, *) "      stride =", parts(i_part)%send_stride(i)
            write (funit, *) "      data   =", data(idx:k-1)
          endif
        enddo
      enddo


      if (i_exch == 1) then
        if (verbose) then
          write (funit, *) "Blocking exchange"
        endif
        ! Blocking exchange
        call pdm_part_comm_graph_exch(pcg_vtx,                   &
                                      PDM_STRIDE_VAR_INTERLACED, &
                                      cst_stride,                &
                                      send_stride,               &
                                      send_data,                 &
                                      recv_stride,               &
                                      recv_data)

      else if (i_exch == 2) then
        if (verbose) then
          write (funit, *) "Non-blocking exchange"
        endif
        ! Initiate non-blocking exchange
        call pdm_part_comm_graph_iexch(pcg_vtx,                   &
                                       PDM_MPI_COMM_KIND_P2P,     &
                                       PDM_STRIDE_VAR_INTERLACED, &
                                       cst_stride,                &
                                       send_stride,               &
                                       send_data,                 &
                                       recv_stride,               &
                                       recv_data,                 &
                                       request)

        ! Do stuff here to cover MPI communications...

        ! Wait for exchange to finish
        call pdm_part_comm_graph_exch_wait(pcg_vtx, request)

      endif


      ! Check the received data
      if (verbose) then
        write (funit, *) "RECV", i_exch
      endif

      do i_part = 1, n_part

        if (verbose) then
          write (funit, *) "part ", i_part
        endif

        call pdm_pointer_array_part_get(recv_stride, &
                                        i_part-1,    &
                                        stride)


        call pdm_pointer_array_part_get(recv_data, &
                                        i_part-1,  &
                                        data)

        k = 1
        do i = 1, pn_vtx_graph(i_part)

          if (verbose) then
            write (funit, *) "  from", parts(i_part)%vtx_graph(4*(i-1)+4), parts(i_part)%vtx_graph(4*(i-1)+2), parts(i_part)%vtx_graph(4*(i-1)+3)
            write (funit, *) "      stride =", stride(i)
            write (funit, *) "      data   =", data(k:k+stride(i)-1)
          endif

          do j = 1, stride(i)
            expected = 100.0d0*parts(i_part)%vtx_graph(4*(i-1)+4) + &
                        10.0d0*parts(i_part)%vtx_graph(4*(i-1)+2) + &
                         1.0d0*parts(i_part)%vtx_graph(4*(i-1)+3) + &
                         0.1d0*(j-1)
            expected = expected * i_exch

            diff = abs(data(k) - expected)
            if (diff > 1.e-9) then
              print *, "Error i_exch", i_exch, ": (", i_rank, i_part, i, j, ") expected ", expected, " but received ", data(k), " from", parts(i_part)%vtx_graph(4*(i-1)+2:4*i)
              stop
            endif

            k = k + 1

          enddo
        enddo

      enddo


      call pdm_pointer_array_free(recv_stride)
      call pdm_pointer_array_free(recv_data)

    enddo

    do i_part = 1, n_part
      deallocate(parts(i_part)%send_stride, &
                 parts(i_part)%send_data)
    enddo

    call pdm_pointer_array_free(send_stride)
    call pdm_pointer_array_free(send_data)

  end subroutine exchanges_variable_stride


end program test_part_comm_graph
