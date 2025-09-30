!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2022  ONERA
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

module pdm_domain_interface

  use pdm
  use pdm_pointer_array

  implicit none

  interface PDM_domain_interface_create ; module procedure  &
    PDM_domain_interface_create_
  end interface

  interface PDM_domain_interface_set ; module procedure  &
    PDM_domain_interface_set_
  end interface

  interface PDM_domain_interface_get ; module procedure  &
    PDM_domain_interface_get_
  end interface

  interface PDM_domain_interface_n_interface_get ; module procedure  &
    PDM_domain_interface_n_interface_get_
  end interface

interface

  !>
  !! \brief Create a set of domain interfaces
  !!
  !! \param [in] n_interface     Number of interfaces to create
  !! \param [in] n_domain        Number of domains
  !! \param [in] mult_intrf      type of interfaces
  !! \param [in] owner           Data ownership
  !! \param [in] comm            MPI communicator
  !!
  !! \return Pointer to \ref PDM_domain_interface object
  !!
  function PDM_domain_interface_create_c (n_interface, &
                                          n_domain,    &
                                          mult_intrf,  &
                                          owner,       &
                                          comm)        &
  result (dom_intrf)                                   &
  bind (c, name='PDM_domain_interface_create')
    use iso_c_binding
    implicit none

    integer(c_int), value :: n_interface
    integer(c_int), value :: n_domain
    integer(c_int), value :: mult_intrf
    integer(c_int), value :: owner
    integer(c_int), value :: comm
    type(c_ptr)           :: dom_intrf

  end function PDM_domain_interface_create_c

  !>
  !! \brief  Set interfaces connectivity
  !!
  !! \param [in] dom_intrf          Pointer to \ref PDm
  !! \param [in] interface_kind     Type of connectivity (vtx/vtx, edge/edge, face/face)
  !! \param [in] interface_dn       Size of the interfaces
  !! \param [in] interface_ids      Index of the entities
  !! \param [in] interface_dom      Index of the domains
  !!
  subroutine PDM_domain_interface_set_c (dom_intrf,      &
                                         interface_kind, &
                                         interface_dn,   &
                                         interface_ids,  &
                                         interface_dom)  &
  bind (c, name='PDM_domain_interface_set')
    use iso_c_binding
    implicit none

    type(c_ptr),    value :: dom_intrf      !< Set of domain interfaces
    integer(c_int), value :: interface_kind !< Type of connected entities
    type(c_ptr),    value :: interface_dn   !< Size of the interfaces
    type(c_ptr),    value :: interface_ids  !< Indexes of the connected entities
    type(c_ptr),    value :: interface_dom  !< Indexes of the connected domains

  end subroutine PDM_domain_interface_set_c

  !>
  !! \brief Get interface connectivity
  !!
  !! \param [in]  dom_intrf
  !! \param [out] interface_kind    Type of connected entities
  !! \param [out] interface_dn      Size of the interfaces
  !! \param [out] interface_ids     Indexes of the connected entities
  !! \param [out] interface_dom     Indexes of the connected entities
  !!
  subroutine PDM_domain_interface_get_c (dom_intrf,      &
                                         interface_kind, &
                                         interface_dn,   &
                                         interface_ids,  &
                                         interface_dom)  &
  bind (c, name='PDM_domain_interface_get')
    use iso_c_binding
    implicit none

    type(c_ptr), value :: dom_intrf      !< Set of domain interfaces
    integer(c_int)     :: interface_kind !< Type of connected entities
    type(c_ptr)        :: interface_dn   !< Size of the interfaces
    type(c_ptr)        :: interface_ids  !< Indexes of the connected entities
    type(c_ptr)        :: interface_dom  !< Indexes of the connected domains

  end subroutine PDM_domain_interface_get_c

  !>
  !! \brief Get number of interfaces
  !!
  !! \param [in] dom_intrf
  !! \param [out] n_interface   Number of interfaces
  !!
  function PDM_domain_interface_n_interface_get_c (dom_intrf)   &
  result (n_interface) &
  bind (c, name='PDM_domain_interface_n_interface_get')
    use iso_c_binding
    implicit none

    type(c_ptr), value :: dom_intrf    !< Set od domain interfaces
    integer(c_int)     :: n_interface  !< Number of defined interfaces

  end function PDM_domain_interface_n_interface_get_c

end interface

contains

  !>
  !! \brief Create a set of domain interfaces
  !!
  !! \param [in] n_interface     Number of interfaces to create
  !! \param [in] n_domain        Number of domains
  !! \param [in] mult_intrf      type of interfaces
  !! \param [in] owner           Data ownership
  !! \param [in] comm            MPI communicator
  !!
  subroutine PDM_domain_interface_create_ (dom_intrf,   &
                                           n_interface, &
                                           n_domain,    &
                                           mult_intrf,  &
                                           owner,       &
                                           comm)

    use iso_c_binding
    implicit none

    type(c_ptr)         :: dom_intrf   !< Set of domain interfaces
    integer, intent(in) :: n_interface !< Number of interfaces
    integer, intent(in) :: n_domain    !< Number of domains
    integer, intent(in) :: mult_intrf  !< Type of the interfaces
    integer, intent(in) :: owner       !< Data ownership
    integer, intent(in) :: comm        !< MPI communicator

    integer(c_int)      :: c_n_intrf   !< C-style number of interfaces
    integer(c_int)      :: c_n_dom     !< C-style number of domains
    integer(c_int)      :: c_mult_int  !< C-style type of interface
    integer(c_int)      :: c_owner     !< C-style ownership
    integer(c_int)      :: c_comm      !< C-style MPI communicator

    c_n_intrf  = n_interface
    c_n_dom    = n_domain
    c_mult_int = mult_intrf
    c_owner    = owner
    c_comm     = PDM_MPI_Comm_f2c(comm)

    dom_intrf = PDM_domain_interface_create_c (c_n_intrf,  &
                                               c_n_dom,    &
                                               c_mult_int, &
                                               c_owner,    &
                                               c_comm)

  end subroutine PDM_domain_interface_create_

  !>
  !! \brief  Set interfaces connectivity
  !!
  !! \param [in] interface_kind     Type of connectivity (vtx/vtx, edge/edge, face/face)
  !! \param [in] interface_dn       Size of the interfaces
  !! \param [in] interface_ids      Index of the entities
  !! \param [in] interface_dom      Index of the domains
  !!

  subroutine PDM_domain_interface_set_ (dom_intrf,      &
                                        interface_kind, &
                                        interface_dn,   &
                                        interface_ids,  &
                                        interface_dom)

    use iso_c_binding
    implicit none

    type(c_ptr),               value      :: dom_intrf      !< Set of domain interfaces
    integer,                   intent(in) :: interface_kind !< Type of the connected entities
    integer(c_int), pointer,   intent(in) :: interface_dn   !< Size of the interfaces
    type(PDM_pointer_array_t), intent(in) :: interface_ids  !< Indexes of the entities
    type(PDM_pointer_array_t), intent(in) :: interface_dom  !< Indexes of the domains

    call PDM_domain_interface_set_c (dom_intrf,                 &
                                     interface_kind,            &
                                     c_loc(interface_dn),       &
                                     c_loc(interface_ids%cptr), &
                                     c_loc(interface_dom%cptr))

  end subroutine

  !>
  !! \brief Get interface connectivity
  !!
  !! \param [in]  dom_intrf
  !! \param [out] interface_kind !< Type of connected entities
  !! \param [out] interface_dn   !< Size of the interfaces
  !! \param [out] interface_ids  !< Indexes of the connected entities
  !! \param [out] interface_dom  !< Indexes of the connected entities
  !!
  subroutine PDM_domain_interface_get_ (dom_intrf,      &
                                        interface_kind, &
                                        interface_dn,   &
                                        interface_ids,  &
                                        interface_dom)

    use iso_c_binding
    implicit none

    type(c_ptr), value                                  :: dom_intrf       !< Set of domain interfaces
    integer,                              intent(out)   :: interface_kind  !< Type of connected entities
    integer(kind = PDM_l_num_s), pointer, intent(inout) :: interface_dn(:) !< Size of the interfaces
    type(PDM_pointer_array_t),   pointer, intent(inout) :: interface_ids   !< Indexes of the connected entities
    type(PDM_pointer_array_t),   pointer, intent(inout) :: interface_dom   !< Indexes of the connected domains

    integer        :: n_interface     !< F-Style interface number
    integer(c_int) :: c_intrf_numb    !< C-style interface number
    integer(c_int) :: c_intrf_kind    !< C-style type of connected entities
    type(c_ptr)    :: c_interface_dn  !< C-style size of the connected interfaces
    type(c_ptr)    :: c_interface_ids !< C-style indexes of the connected entities
    type(c_ptr)    :: c_interface_dom !< C-style indexes of the connected domains


    call PDM_domain_interface_get_c (dom_intrf, &
                                     c_intrf_kind, &
                                     c_interface_dn, &
                                     c_interface_ids, &
                                     c_interface_dom)

    c_intrf_numb = PDM_domain_interface_n_interface_get_c (dom_intrf)

    n_interface = c_intrf_numb
    interface_kind = c_intrf_kind

    call c_f_pointer(c_interface_dn,   &
                     interface_dn, &
                     [c_intrf_numb])

    call PDM_pointer_array_create (interface_ids,      &
                                   n_interface,        &
                                   PDM_TYPE_G_NUM,     &
                                   c_interface_ids,    &
                                   interface_dn,       &
                                   PDM_OWNERSHIP_USER)

    call PDM_pointer_array_create (interface_dom,      &
                                   n_interface,        &
                                   PDM_TYPE_INT,       &
                                   c_interface_dom,    &
                                   interface_dn,       &
                                   PDM_OWNERSHIP_USER)

  end subroutine PDM_domain_interface_get_

  !>
  !! \brief Get interface number
  !!
  !! \param [in]  dom_intrf
  !! \param [out] n_interface    !< number of interfaces
  !!
  subroutine PDM_domain_interface_n_interface_get_ (dom_intrf,      &
                                                    n_interface)

    use iso_c_binding
    implicit none

    type(c_ptr), value                                  :: dom_intrf      !< Set of domain interfaces
    integer,                              intent(out)   :: n_interface    !< Type of connected entities

    integer(c_int) :: c_intrf_numb    !< C-style interface number


    c_intrf_numb = PDM_domain_interface_n_interface_get_c (dom_intrf)

    n_interface = c_intrf_numb

  end subroutine PDM_domain_interface_n_interface_get_


end module pdm_domain_interface
