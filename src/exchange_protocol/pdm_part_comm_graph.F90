#include "pdm_configf.h"

module pdm_part_comm_graph

  use iso_c_binding
  use pdm
  use pdm_pointer_array


  interface

    ! --- Accessors ---
    subroutine pdm_part_comm_graph_entity_graph_get_cf(pcg,          &
                                                       i_part,       &
                                                       entity_graph, &
                                                       ownership)    &
    bind(c, name="PDM_part_comm_graph_entity_graph_get")
      use iso_c_binding
      implicit none
      type(c_ptr),    value :: pcg
      integer(c_int), value :: i_part
      type(c_ptr)           :: entity_graph
      integer(c_int), value :: ownership
    end subroutine pdm_part_comm_graph_entity_graph_get_cf


    function pdm_part_comm_graph_nuplet_size_cf(pcg) &
    result(size) &
    bind(c, name="PDM_part_comm_graph_nuplet_size")
      use iso_c_binding
      implicit none
      type(c_ptr),   value :: pcg
      integer(c_int)       :: size
    end function pdm_part_comm_graph_nuplet_size_cf


    function pdm_part_comm_graph_n_entity_get_cf(pcg,    &
                                                 i_part) &
    result(n_entity) &
    bind(c, name="PDM_part_comm_graph_n_entity_get")
      use iso_c_binding
      implicit none
      type(c_ptr),    value :: pcg
      integer(c_int), value :: i_part
      integer(c_int)        :: n_entity
    end function pdm_part_comm_graph_n_entity_get_cf


    function pdm_part_comm_graph_n_part_get_cf(pcg) &
    result(n_part) &
    bind(c, name="PDM_part_comm_graph_n_part_get")
      use iso_c_binding
      type(c_ptr), value :: pcg
      integer(c_int)     :: n_part
    end function pdm_part_comm_graph_n_part_get_cf


  end interface


  private :: setup_recv_pa

  contains


  subroutine PDM_part_comm_graph_create(pcg,             &
                                        n_part,          &
                                        pn_entity_graph, &
                                        pentity_graph,   &
                                        ownership,       &
                                        comm)
    ! Build a Part Comm Graph instance
    implicit none

    type(c_ptr),                        intent(out) :: pcg                ! PDM_part_comm_graph_t instance
    integer,                            intent(in)  :: n_part             ! Number of parts on current process
    integer(pdm_l_num_s),      pointer, intent(in)  :: pn_entity_graph(:) ! Number of part boundary entities (size = ``n_part``)
    type(pdm_pointer_array_t), pointer, intent(in)  :: pentity_graph      ! Inter-part communication graph (size = 4*``pn_entity_graph`` : [local entity ID (1-based), connected rank (0-based), connected part (1-based), ID of connected entity (1-based)])
    integer,                            intent(in)  :: ownership          ! Ownership
    integer,                            intent(in)  :: comm               ! MPI communicator

    type(c_ptr)                                     :: c_comm

    interface
      function pdm_part_comm_graph_create_cf(n_part,          &
                                             pn_entity_graph, &
                                             pentity_graph,   &
                                             ownership,       &
                                             comm)            &
      result (pcg) &
      bind(c, name="PDM_part_comm_graph_create")
        use iso_c_binding
        implicit none
        integer(c_int), value :: n_part
        type(c_ptr),    value :: pn_entity_graph
        type(c_ptr),    value :: pentity_graph
        integer(c_int), value :: ownership
        type(c_ptr),    value :: comm
        type(c_ptr)           :: pcg
      end function pdm_part_comm_graph_create_cf
    end interface

    c_comm = PDM_MPI_Comm_f2c(comm)

    pcg = pdm_part_comm_graph_create_cf(n_part,                    &
                                        c_loc(pn_entity_graph),    &
                                        c_loc(pentity_graph%cptr), &
                                        ownership,                 &
                                        c_comm)

  end subroutine pdm_part_comm_graph_create



  subroutine PDM_part_comm_graph_with_nuplet_create(pcg,              &
                                                    n_part,           &
                                                    pn_entity_graph,  &
                                                    pentity_graph,    &
                                                    ownership_graph,  &
                                                    nuplet_size,      &
                                                    pentity_nuplet,   &
                                                    ownership_nuplet, &
                                                    is_signed,        &
                                                    comm)
    ! Build a Part Comm Graph instance using additional information represented as a n-uplet
    implicit none

    type(c_ptr),                        intent(out) :: pcg                ! PDM_part_comm_graph_t instance
    integer,                            intent(in)  :: n_part             ! Number of parts on current process
    integer(pdm_l_num_s),      pointer, intent(in)  :: pn_entity_graph(:) ! Number of part boundary entities (size = ``n_part``)
    type(pdm_pointer_array_t), pointer, intent(in)  :: pentity_graph      ! Inter-part communication graph (size = 4*``pn_entity_graph`` : [local entity ID (1-based), connected rank (0-based), connected part (1-based), ID of connected entity (1-based)])
    integer,                            intent(in)  :: ownership_graph    ! Ownership for ``pentity_graph``
    integer,                            intent(in)  :: nuplet_size        ! N-uplet size
    type(pdm_pointer_array_t), pointer, intent(in)  :: pentity_nuplet     ! Additional nuplets (size = ``nuplet_size * pn_entity_graph``)
    integer,                            intent(in)  :: ownership_nuplet   ! Ownership for ``pentity_nuplet``
    logical,                            intent(in)  :: is_signed          ! Use signed nuplets
    integer,                            intent(in)  :: comm               ! MPI communicator

    integer                                         :: c_is_signed
    type(c_ptr)                                     :: c_comm

    interface
      function pdm_part_comm_graph_with_nuplet_create_cf(n_part,           &
                                                         pn_entity_graph,  &
                                                         pentity_graph,    &
                                                         ownership_graph,  &
                                                         nuplet_size,      &
                                                         pentity_nuplet,   &
                                                         ownership_nuplet, &
                                                         is_signed,        &
                                                         comm)             &
      result (pcg) &
      bind(c, name="PDM_part_comm_graph_with_nuplet_create")
        use iso_c_binding
        implicit none
        integer(c_int), value :: n_part
        type(c_ptr),    value :: pn_entity_graph
        type(c_ptr),    value :: pentity_graph
        integer(c_int), value :: ownership_graph
        integer(c_int), value :: nuplet_size
        type(c_ptr),    value :: pentity_nuplet
        integer(c_int), value :: ownership_nuplet
        integer(c_int), value :: is_signed
        type(c_ptr),    value :: comm
        type(c_ptr)           :: pcg
      end function PDM_part_comm_graph_with_nuplet_create_cf
    end interface

    if (is_signed) then
      c_is_signed = 1
    else
      c_is_signed = 0
    endif

    c_comm = PDM_MPI_Comm_f2c(comm)

    pcg = pdm_part_comm_graph_with_nuplet_create_cf(n_part,                     &
                                                    c_loc(pn_entity_graph),     &
                                                    c_loc(pentity_graph%cptr),  &
                                                    ownership_graph,            &
                                                    nuplet_size,                &
                                                    c_loc(pentity_nuplet%cptr), &
                                                    ownership_nuplet,           &
                                                    c_is_signed,                &
                                                    c_comm)

  end subroutine PDM_part_comm_graph_with_nuplet_create



  function PDM_part_comm_graph_is_signed(pcg) result(is_signed)
    ! Return .true. if nuplet description is signed, else .false.
    implicit none

    type(c_ptr), intent(in) :: pcg       ! PDM_part_comm_graph_t instance
    logical                 :: is_signed ! Is the nuplet signed?

    integer                 :: c_is_signed

    interface
      function pdm_part_comm_graph_is_signed_cf(pcg) result(is_signed) &
      bind(c, name="PDM_part_comm_graph_is_signed")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: pcg
        integer(c_int)     :: is_signed
      end function pdm_part_comm_graph_is_signed_cf
    end interface

    c_is_signed = pdm_part_comm_graph_is_signed_cf(pcg)

    if (c_is_signed == 1) then
      is_signed = .true.
    else
      is_signed = .false.
    endif

  end function PDM_part_comm_graph_is_signed



  function PDM_part_comm_graph_nuplet_size(pcg) result(size)
    ! Get nuplet size
    implicit none

    type(c_ptr), intent(in) :: pcg  ! PDM_part_comm_graph_t instance
    integer                 :: size ! Size of nuplet

    size = pdm_part_comm_graph_nuplet_size_cf(pcg)

  end function PDM_part_comm_graph_nuplet_size



  function PDM_part_comm_graph_n_entity_get(pcg,    &
                                            i_part) &
  result(n_entity)
    ! Get number of graph entities
    implicit none

    type(c_ptr), intent(in) :: pcg      ! PDM_part_comm_graph_t instance
    integer,     intent(in) :: i_part   ! Partition identifier
    integer                 :: n_entity ! Number of graph entities

    n_entity = pdm_part_comm_graph_n_entity_get_cf(pcg, i_part)

  end function PDM_part_comm_graph_n_entity_get



  subroutine PDM_part_comm_graph_entity_nuplet_get(pcg,           &
                                                   i_part,        &
                                                   entity_nuplet, &
                                                   ownership)
    ! Get entity nuplets
    implicit none

    type(c_ptr),                   intent(in)  :: pcg              ! PDM_part_comm_graph_t instance
    integer,                       intent(in)  :: i_part           ! Partition identifier
    integer(pdm_l_num_s), pointer, intent(out) :: entity_nuplet(:) ! Entity nuplets (size = nuplet_size * n_entity_graph)
    integer,                       intent(in)  :: ownership        ! Ownership

    integer(c_int)                             :: n_entity
    integer(c_int)                             :: nuplet_size
    type(c_ptr)                                :: c_entity_nuplet

    interface
      subroutine pdm_part_comm_graph_entity_nuplet_get_cf(pcg,           &
                                                          i_part,        &
                                                          entity_nuplet, &
                                                          ownership)     &
      bind(c, name="PDM_part_comm_graph_entity_nuplet_get")
        use iso_c_binding
        implicit none
        type(c_ptr),    value :: pcg
        integer(c_int), value :: i_part
        type(c_ptr)           :: entity_nuplet
        integer(c_int), value :: ownership
      end subroutine pdm_part_comm_graph_entity_nuplet_get_cf
    end interface

    nuplet_size = pdm_part_comm_graph_nuplet_size_cf(pcg)
    n_entity    = pdm_part_comm_graph_n_entity_get_cf(pcg, i_part)

    call pdm_part_comm_graph_entity_nuplet_get_cf(pcg,             &
                                                  i_part,          &
                                                  c_entity_nuplet, &
                                                  ownership)

    call c_f_pointer(c_entity_nuplet,          &
                     entity_nuplet,            &
                     [nuplet_size * n_entity])

  end subroutine PDM_part_comm_graph_entity_nuplet_get



  subroutine PDM_part_comm_graph_entity_graph_get(pcg,          &
                                                  i_part,       &
                                                  entity_graph, &
                                                  ownership)
    ! Get entity graph
    implicit none

    type(c_ptr),                   intent(in)  :: pcg             ! PDM_part_comm_graph_t instance
    integer,                       intent(in)  :: i_part          ! Partition identifier
    integer(pdm_l_num_s), pointer, intent(out) :: entity_graph(:) ! Entity graph (size = 4 * n_entity_graph)
    integer,                       intent(in)  :: ownership       ! Ownership

    integer(c_int)                             :: n_entity
    type(c_ptr)                                :: c_entity_graph

    interface
      subroutine pdm_part_comm_graph_entity_graph_get_cf(pcg,          &
                                                         i_part,       &
                                                         entity_graph, &
                                                         ownership)    &
      bind(c, name="PDM_part_comm_graph_entity_graph_get")
        use iso_c_binding
        implicit none
        type(c_ptr),    value :: pcg
        integer(c_int), value :: i_part
        type(c_ptr)           :: entity_graph
        integer(c_int), value :: ownership
      end subroutine pdm_part_comm_graph_entity_graph_get_cf
    end interface

    n_entity = pdm_part_comm_graph_n_entity_get_cf(pcg, i_part)

    call pdm_part_comm_graph_entity_graph_get_cf(pcg,            &
                                                 i_part,         &
                                                 c_entity_graph, &
                                                 ownership)

    call c_f_pointer(c_entity_graph,  &
                     entity_graph,    &
                     [4 *  n_entity])

  end subroutine PDM_part_comm_graph_entity_graph_get



  subroutine PDM_part_comm_graph_owner_get(pcg,    &
                                           i_part, &
                                           owner)
    ! Get the owner array computed inside the structure, useful to manage reduction of array for example
    implicit none
    type(c_ptr),                   intent(in)  :: pcg      ! PDM_part_comm_graph_t instance
    integer,                       intent(in)  :: i_part   ! Partition identifier
    integer(pdm_l_num_s), pointer, intent(out) :: owner(:) ! Owner status (size = n_entity_graph)

    integer(c_int)                             :: n_entity
    type(c_ptr)                                :: c_owner

    interface
      function pdm_part_comm_graph_owner_get_cf(pcg, i_part) result(res) &
      bind(c, name="PDM_part_comm_graph_owner_get")
        use iso_c_binding
        type(c_ptr),    value :: pcg
        integer(c_int), value :: i_part
        type(c_ptr)           :: res
      end function
    end interface

    n_entity = pdm_part_comm_graph_n_entity_get_cf(pcg, i_part)

    c_owner = pdm_part_comm_graph_owner_get_cf(pcg, i_part)

    call c_f_pointer(c_owner, owner, [n_entity])

  end subroutine PDM_part_comm_graph_owner_get



  subroutine PDM_part_comm_graph_all_reduce(pcg,    &
                                            stride, &
                                            op,     &
                                            pdata)
    ! Inplace reduction value on current graph. Allow synchronization.
    ! Only MPI_DOUBLE and MPI_INT data types are supported
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
    use mpi
#endif
    implicit none
#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
    include "mpif.h"
#endif
    type(c_ptr),               intent(in) :: pcg    ! PDM_part_comm_graph_t instance
    integer,                   intent(in) :: stride ! Constant stride
    integer,                   intent(in) :: op     ! Reduction operation kind (MPI_SUM/MPI_MIN/MPI_MAX)
    type(pdm_pointer_array_t), pointer    :: pdata  ! Buffer of data to synchronise (size = n_entity)

    type(c_ptr)                           :: c_datatype
    type(c_ptr)                           :: c_op

    interface
      subroutine pdm_part_comm_graph_all_reduce_cf(pcg,      &
                                                   datatype, &
                                                   stride,   &
                                                   op,       &
                                                   pdata)    &
      bind(c, name="PDM_part_comm_graph_all_reduce")
        use iso_c_binding
        implicit none
        type(c_ptr),    value :: pcg
        type(c_ptr),    value :: datatype
        integer(c_int), value :: stride
        type(c_ptr),    value :: op
        type(c_ptr),    value :: pdata
      end subroutine
    end interface

    if (pdata%type == PDM_TYPE_INT) then
      c_datatype = PDM_MPI_Type_f2c(MPI_INT)
    else if (pdata%type == PDM_TYPE_DOUBLE) then
      c_datatype = PDM_MPI_Type_f2c(MPI_DOUBLE)
    else
      print *, "PDM_part_comm_graph_all_reduce: data type ", pdata%type, " is not supported"
      stop
    end if

    c_op = PDM_MPI_Op_f2c(op)

    call pdm_part_comm_graph_all_reduce_cf(pcg,               &
                                           c_datatype,        &
                                           stride,            &
                                           c_op,              &
                                           c_loc(pdata%cptr))

  end subroutine PDM_part_comm_graph_all_reduce



  subroutine PDM_part_comm_graph_exch(pcg,         &
                                      t_stride,    &
                                      cst_stride,  &
                                      send_stride, &
                                      send_data,   &
                                      recv_stride, &
                                      recv_data)
    ! Exchange data
    implicit none

    type(c_ptr),               intent(in)  :: pcg         ! PDM_part_comm_graph_t instance
    integer,                   intent(in)  :: t_stride    ! Type of stride
    integer,                   intent(in)  :: cst_stride  ! Constant stride
    type(pdm_pointer_array_t), pointer     :: send_stride ! Stride of send data
    type(pdm_pointer_array_t), pointer     :: send_data   ! Send data
    type(pdm_pointer_array_t), pointer     :: recv_stride ! Stride of recv data
    type(pdm_pointer_array_t), pointer     :: recv_data   ! Recv data

    integer(c_size_t)                      :: c_s_data
    type(c_ptr)                            :: c_send_stride
    type(c_ptr)                            :: c_recv_stride
    type(c_ptr)                            :: c_recv_data

    interface
      subroutine pdm_part_comm_graph_exch_cf(pcg,         &
                                             s_data,      &
                                             t_stride,    &
                                             cst_stride,  &
                                             send_stride, &
                                             send_data,   &
                                             recv_stride, &
                                             recv_data)   &
      bind(c, name="PDM_part_comm_graph_exch")
        use iso_c_binding
        implicit none
        type(c_ptr),       value :: pcg
        integer(c_size_t), value :: s_data
        integer(c_int),    value :: t_stride
        integer(c_int),    value :: cst_stride
        type(c_ptr),       value :: send_stride
        type(c_ptr),       value :: send_data
        type(c_ptr)              :: recv_stride
        type(c_ptr)              :: recv_data
      end subroutine pdm_part_comm_graph_exch_cf
    end interface

    c_s_data = send_data%s_data

    c_send_stride = C_NULL_PTR
    if (associated(send_stride)) then
      c_send_stride = c_loc(send_stride%cptr)
    endif

    c_recv_stride = C_NULL_PTR
    c_recv_data   = C_NULL_PTR

    call pdm_part_comm_graph_exch_cf(pcg,                   &
                                     c_s_data,              &
                                     t_stride,              &
                                     cst_stride,            &
                                     c_send_stride,         &
                                     c_loc(send_data%cptr), &
                                     c_recv_stride,         &
                                     c_recv_data)

    call setup_recv_pa(pcg,              &
                       t_stride,         &
                       cst_stride,       &
                       send_data%type,   &
                       send_data%s_data, &
                       c_recv_stride,    &
                       c_recv_data,      &
                       recv_stride,      &
                       recv_data)

  end subroutine PDM_part_comm_graph_exch



  subroutine PDM_part_comm_graph_iexch(pcg,         &
                                       k_comm,      &
                                       t_stride,    &
                                       cst_stride,  &
                                       send_stride, &
                                       send_data,   &
                                       recv_stride, &
                                       recv_data,   &
                                       request)
    ! Initiate a non-blocking exchange
    implicit none

    type(c_ptr),               intent(in)  :: pcg         ! PDM_part_comm_graph_t instance
    integer,                   intent(in)  :: k_comm      ! Kind of MPI communication
    integer,                   intent(in)  :: t_stride    ! Type of stride
    integer,                   intent(in)  :: cst_stride  ! Constant stride
    type(pdm_pointer_array_t), pointer     :: send_stride ! Stride of send data
    type(pdm_pointer_array_t), pointer     :: send_data   ! Send data
    type(pdm_pointer_array_t), pointer     :: recv_stride ! Stride of recv data
    type(pdm_pointer_array_t), pointer     :: recv_data   ! Recv data
    integer,                   intent(out) :: request     ! Request ID

    integer(c_size_t)                      :: c_s_data
    type(c_ptr)                            :: c_send_stride
    type(c_ptr)                            :: c_recv_stride
    type(c_ptr)                            :: c_recv_data

    interface
      function pdm_part_comm_graph_iexch_cf(pcg,         &
                                            k_comm,      &
                                            s_data,      &
                                            t_stride,    &
                                            cst_stride,  &
                                            send_stride, &
                                            send_data,   &
                                            recv_stride, &
                                            recv_data)   &
      result (request)                                   &
      bind(c, name="PDM_part_comm_graph_iexch")
        use iso_c_binding
        implicit none
        type(c_ptr),       value :: pcg
        integer(c_int),    value :: k_comm
        integer(c_size_t), value :: s_data
        integer(c_int),    value :: t_stride
        integer(c_int),    value :: cst_stride
        type(c_ptr),       value :: send_stride
        type(c_ptr),       value :: send_data
        type(c_ptr)              :: recv_stride
        type(c_ptr)              :: recv_data
        integer(c_int)           :: request
      end function pdm_part_comm_graph_iexch_cf
    end interface

    c_s_data = send_data%s_data

    c_send_stride = C_NULL_PTR
    if (associated(send_stride)) then
      c_send_stride = c_loc(send_stride%cptr)
    endif

    c_recv_stride = C_NULL_PTR
    c_recv_data   = C_NULL_PTR

    request = pdm_part_comm_graph_iexch_cf(pcg,                   &
                                           k_comm,                &
                                           c_s_data,              &
                                           t_stride,              &
                                           cst_stride,            &
                                           c_send_stride,         &
                                           c_loc(send_data%cptr), &
                                           c_recv_stride,         &
                                           c_recv_data)

    call setup_recv_pa(pcg,              &
                       t_stride,         &
                       cst_stride,       &
                       send_data%type,   &
                       send_data%s_data, &
                       c_recv_stride,    &
                       c_recv_data,      &
                       recv_stride,      &
                       recv_data)

  end subroutine PDM_part_comm_graph_iexch



  subroutine PDM_part_comm_graph_exch_wait(pcg,     &
                                           request)
    ! Wait for a non-blocking exchange to finish
    implicit none

    type(c_ptr), intent(in) :: pcg     ! PDM_part_comm_graph_t instance
    integer,     intent(in) :: request ! Request ID

    interface
      subroutine pdm_part_comm_graph_exch_wait_cf(pcg,     &
                                                  request) &
      bind(c, name="PDM_part_comm_graph_exch_wait")
        use iso_c_binding
        implicit none
        type(c_ptr),    value :: pcg
        integer(c_int), value :: request
      end subroutine pdm_part_comm_graph_exch_wait_cf
    end interface

    call pdm_part_comm_graph_exch_wait_cf(pcg,     &
                                          request)

  end subroutine PDM_part_comm_graph_exch_wait



  subroutine PDM_part_comm_graph_reorder(pcg,        &
                                         old_to_new)
    ! Reorder the graph entities
    implicit none

    type(c_ptr),               intent(in) :: pcg        ! PDM_part_comm_graph_t instance
    type(pdm_pointer_array_t), pointer    :: old_to_new ! Permutation table (0-based)

    interface
      subroutine pdm_part_comm_graph_reorder_cf(pcg,        &
                                                old_to_new) &
      bind(c, name="PDM_part_comm_graph_reorder")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: pcg
        type(c_ptr), value :: old_to_new
      end subroutine pdm_part_comm_graph_reorder_cf
    end interface

    call pdm_part_comm_graph_reorder_cf(pcg,                    &
                                        c_loc(old_to_new%cptr))

  end subroutine PDM_part_comm_graph_reorder



  subroutine PDM_part_comm_graph_free(pcg)
    ! Free a PDM_part_comm_graph instance
    implicit none

    type(c_ptr), intent(inout) :: pcg ! PDM_part_comm_graph_t instance

    interface
      subroutine pdm_part_comm_graph_free_cf(pcg) &
      bind(c, name="PDM_part_comm_graph_free")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: pcg
      end subroutine pdm_part_comm_graph_free_cf
    end interface

    call pdm_part_comm_graph_free_cf(pcg)

  end subroutine PDM_part_comm_graph_free



  ! --- Auxiliary procedures ---

  subroutine setup_recv_pa(pcg,           &
                           t_stride,      &
                           cst_stride,    &
                           data_type,     &
                           s_data,        &
                           c_recv_stride, &
                           c_recv_data,   &
                           recv_stride,   &
                           recv_data)
    ! Setup the recv* pointer_arrays for (i)exch routines
    implicit none

    type(c_ptr),               intent(in) :: pcg
    integer,                   intent(in) :: t_stride
    integer,                   intent(in) :: cst_stride
    integer,                   intent(in) :: data_type
    integer,                   intent(in) :: s_data
    type(c_ptr),               intent(in) :: c_recv_stride
    type(c_ptr),               intent(in) :: c_recv_data
    type(pdm_pointer_array_t), pointer    :: recv_stride
    type(pdm_pointer_array_t), pointer    :: recv_data

    integer(c_int)                        :: n_part, i_part, i_entity
    integer(pdm_l_num_s),      pointer    :: length_stride(:)
    integer(pdm_l_num_s),      pointer    :: length_data(:)
    integer(pdm_l_num_s),      pointer    :: stride(:)

    n_part = pdm_part_comm_graph_n_part_get_cf(pcg)

    allocate(length_data(n_part))

    if (t_stride == PDM_STRIDE_VAR_INTERLACED) then
      ! Variable stride
      allocate(length_stride(n_part))
      do i_part = 1, n_part
        length_stride(i_part) = pdm_part_comm_graph_n_entity_get_cf(pcg, i_part-1)
      enddo

      call pdm_pointer_array_create(recv_stride,        &
                                    n_part,             &
                                    PDM_TYPE_INT,       &
                                    c_recv_stride,      &
                                    length_stride,      &
                                    PDM_OWNERSHIP_KEEP)

      do i_part = 1, n_part
        call PDM_pointer_array_part_get(recv_stride, &
                                        i_part-1,    &
                                        stride)
        length_data(i_part) = 0
        do i_entity = 1, length_stride(i_part)
          length_data(i_part) = length_data(i_part) + stride(i_entity)
        enddo
      enddo

    else
      ! Constant stride
      do i_part = 1, n_part
        length_data(i_part) = cst_stride * pdm_part_comm_graph_n_entity_get_cf(pcg, i_part-1)
      enddo

    endif


    if (data_type == PDM_TYPE_CPTR) then
      call pdm_pointer_array_create(recv_data,          &
                                    n_part,             &
                                    data_type,          &
                                    c_recv_data,        &
                                    length_data,        &
                                    PDM_OWNERSHIP_KEEP)
    else
      call pdm_pointer_array_create(recv_data,          &
                                    n_part,             &
                                    data_type,          &
                                    c_recv_data,        &
                                    length_data,        &
                                    PDM_OWNERSHIP_KEEP, &
                                    s_data)
    endif

    if (t_stride == PDM_STRIDE_VAR_INTERLACED) then
      deallocate(length_stride)
    endif
    deallocate(length_data)

  end subroutine setup_recv_pa


end module pdm_part_comm_graph
