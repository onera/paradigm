#include "pdm_configf.h"

module pdm_partitioning_algorithm

use pdm
use pdm_pointer_array

implicit none

interface

  subroutine PDM_compute_graph_comm_entity_ownerhip_single_part_cf (n_entity,          &
                                                                    entity_ln_to_gn,   &
                                                                    n_owned_entity,    &
                                                                    lnum_owned_entity, &
                                                                    comm)              &
    bind(c, name='PDM_compute_graph_comm_entity_ownerhip_single_part')

    use iso_c_binding
    implicit none

    integer(c_int), value :: n_entity
    type(c_ptr),    value :: entity_ln_to_gn
    integer(c_int)        :: n_owned_entity
    type(c_ptr)           :: lnum_owned_entity
    integer(c_int), value :: comm

  end subroutine PDM_compute_graph_comm_entity_ownerhip_single_part_cf

end interface

contains

  !>
  !!
  !! \brief Get the list of owned entities on the current process
  !!
  !! \param [in]  n_part            Number of partitions
  !! \param [in]  n_entity          Number of entities
  !! \param [in]  entity_ln_to_gn   Entity local numbering to global numbering (size = n_entity)
  !! \param [out] n_owned_entity    Number of owned entities
  !! \param [out] lnum_owned_entity Owned entity local numbering (size = n_owned_entity)
  !! \param [in]  comm              MPI communicator
  !!
  !!

  subroutine PDM_compute_graph_comm_entity_ownerhip (n_part,            &
                                                     n_entity,          &
                                                     entity_ln_to_gn,   &
                                                     n_owned_entity,    &
                                                     lnum_owned_entity, &
                                                     f_comm)

    use iso_c_binding
    implicit none

    integer, intent(in)                :: n_part
    integer, pointer                   :: n_entity(:)
    type(PDM_pointer_array_t), pointer :: entity_ln_to_gn
    integer, pointer                   :: n_owned_entity(:)
    type(PDM_pointer_array_t), pointer :: lnum_owned_entity
    integer, intent(in)                :: f_comm

    integer(c_int)                     :: c_comm
    type(c_ptr)                        :: c_n_owned_entity
    type(c_ptr)                        :: c_lnum_owned_entity

    interface

      subroutine PDM_compute_graph_comm_entity_ownerhip_c (n_part,            &
                                                           n_entity,          &
                                                           entity_ln_to_gn,   &
                                                           n_owned_entity,    &
                                                           lnum_owned_entity, &
                                                           comm)              &
        bind(c, name='PDM_compute_graph_comm_entity_ownerhip')

        use iso_c_binding
        implicit none

        integer(c_int), value :: n_part
        type(c_ptr),    value :: n_entity
        type(c_ptr),    value :: entity_ln_to_gn
        type(c_ptr)           :: n_owned_entity
        type(c_ptr)           :: lnum_owned_entity
        integer(c_int), value :: comm

      end subroutine PDM_compute_graph_comm_entity_ownerhip_c

    end interface

    c_comm = PDM_MPI_Comm_f2c(f_comm)

    call PDM_compute_graph_comm_entity_ownerhip_c (n_part,                      &
                                                   c_loc(n_entity),             &
                                                   c_loc(entity_ln_to_gn%cptr), &
                                                   c_n_owned_entity,            &
                                                   c_lnum_owned_entity,         &
                                                   c_comm)

    call c_f_pointer(c_n_owned_entity, &
                     n_owned_entity,   &
                     [n_part])

    call PDM_pointer_array_create (lnum_owned_entity,   &
                                   n_part,              &
                                   PDM_TYPE_INT,        &
                                   c_lnum_owned_entity, &
                                   n_owned_entity,      &
                                   PDM_OWNERSHIP_USER)

  end subroutine PDM_compute_graph_comm_entity_ownerhip

  !>
  !!
  !! \brief Get the list of owned entities on the current process
  !!
  !! \param [in]  n_entity          Number of entities
  !! \param [in]  entity_ln_to_gn   Entity local numbering to global numbering (size = n_entity)
  !! \param [out] n_owned_entity    Number of owned entities
  !! \param [out] lnum_owned_entity Owned entity local numbering (size = n_owned_entity)
  !! \param [in]  comm              MPI communicator
  !!
  !!

  subroutine PDM_compute_graph_comm_entity_ownerhip_single_part (n_entity,          &
                                                                 entity_ln_to_gn,   &
                                                                 n_owned_entity,    &
                                                                 lnum_owned_entity, &
                                                                 f_comm)

    use iso_c_binding
    implicit none

    integer, intent(in)           :: n_entity
    integer(pdm_g_num_s), pointer :: entity_ln_to_gn(:)
    integer                       :: n_owned_entity
    integer, pointer              :: lnum_owned_entity(:)
    integer, intent(in)           :: f_comm

    integer(c_int)                :: c_comm
    type(c_ptr)                   :: c_entity_ln_to_gn
    type(c_ptr)                   :: c_lnum_owned_entity

    c_entity_ln_to_gn = C_NULL_PTR
    if (associated(entity_ln_to_gn)) then
      c_entity_ln_to_gn = c_loc(entity_ln_to_gn)
    end if

    c_comm = PDM_MPI_Comm_f2c(f_comm)

    call PDM_compute_graph_comm_entity_ownerhip_single_part_cf (n_entity,            &
                                                                c_entity_ln_to_gn,   &
                                                                n_owned_entity,      &
                                                                c_lnum_owned_entity, &
                                                                c_comm)

    call c_f_pointer(c_lnum_owned_entity, &
                     lnum_owned_entity,   &
                     [n_owned_entity])

  end subroutine PDM_compute_graph_comm_entity_ownerhip_single_part

end module pdm_partitioning_algorithm
