!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2023  ONERA
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

module pdm_dbbtree

  use pdm
  use pdm_fortran
  use iso_c_binding
  use pdm_pointer_array

  implicit none

  interface 

	  !>
	  !! \brief Free a \ref PDM_dbbtree_t structure
	  !!
	  !! \param [in]  dbbt   Pointer to a distributed bounding box tree
	  !!
	  !! \return      NULL
	  !!
	  !!

    function PDM_dbbtree_free (dbbt) &  
    result (dbbt_null)               &
    bind (c, name='PDM_dbbtree_free')

      use iso_c_binding
      implicit none

      type(c_ptr), value :: dbbt

      type(c_ptr)        :: dbbt_null

    end function PDM_dbbtree_free 

  end interface

  interface PDM_dbbtree_create
    module procedure PDM_dbbtree_create_
  end interface 

  interface PDM_dbbtree_boxes_set
    module procedure PDM_dbbtree_boxes_set_
    module procedure PDM_dbbtree_boxes_set_with_init_location
  end interface PDM_dbbtree_boxes_set

  ! interface PDM_dbbtree_intersect_boxes_set
  !   module procedure PDM_dbbtree_intersect_boxes_set_
  !   module procedure PDM_dbbtree_intersect_boxes_with_init_location_set
  ! end interface PDM_dbbtree_intersect_boxes_set

  type PDM_dbbtree_t

    type(c_ptr) :: c_dbbtree              = C_NULL_PTR
    type(c_ptr) :: bs_boxes_in_tree       = C_NULL_PTR

    type(c_ptr) :: c_boxes_in_tree_intersected_boxes_idx = C_NULL_PTR   ! A liberer dans le free
    type(c_ptr) :: c_boxes_in_tree_intersected_boxes_lnum = C_NULL_PTR  ! A liberer dans le free

    type(c_ptr) :: bs_intersected_boxes   = C_NULL_PTR

  end type PDM_dbbtree_t

contains

    !>
    !! \brief Return an intialized \ref PDM_dbbtree_t structure
    !!
    !! This function returns an initialized \ref PDM_dbbtree_t structure
    !!
    !! \param [in]  comm             Associated communicator
    !! \param [in]  dim              boxes dimension
    !! \param [in]  global_extents   Globals of elements to storage into the tree
    !!                               (automatic computation if NULL)
    !!
    !! \return      A new initialized \ref PDM_dbbtree_t structure
    !!
    
    function  PDM_dbbtree_create_(comm,          &
                                  dim,           & 
                                  global_extents)&
    result (dbbt)   

    use iso_c_binding
    implicit none

    integer,          intent(in)  :: dim
    integer,          intent(in)  :: comm
    double precision, pointer     :: global_extents(:)

    type(PDM_dbbtree_t), pointer  :: dbbt        

    integer(c_int)                :: c_comm

    interface 

      function PDM_dbbtree_create_c (comm,          &
                                     dim,           & 
                                     global_extents)&
      result (dbbt)                                 &
      bind (c, name='PDM_dbbtree_create')
  
        use iso_c_binding
        implicit none

        integer(c_int), value :: comm
        integer(c_int), value :: dim
        type(c_ptr),    value :: global_extents

        type(c_ptr)           :: dbbt          

      end function PDM_dbbtree_create_c    

    end interface

    allocate (dbbt)

    c_comm = PDM_MPI_Comm_f2c(comm)

    dbbt%c_dbbtree = PDM_dbbtree_create_c (c_comm, dim, c_loc(global_extents))

  end function PDM_dbbtree_create_

  !>
  !! \brief Assign a set of boxes to an empty \ref PDM_dbbtree_t structure.
  !!
  !! This function assigns a set of boxes to an empty \ref PDM_dbbtree_t structure.
  !!
  !! \param [in]  dbbt     Pointer to a distributed bounding box tree
  !! \param [in]  n_part   Number of partitions
  !! \param [in]  nElts    Number of elements of each partition
  !! \param [in]  extents  Extents of each element of each partition
  !! \param [in]  gNum     Global number of each element of each partition
  !!
  !! \return associated \ref PDM_box_set_t structure distributed according to
  !! the tree location
  !!

  subroutine PDM_dbbtree_boxes_set_ (dbbt,    &
                                     n_part,  &
                                     nElts,   &
                                     extents, &
                                     gNum     &
                                     )        
    use iso_c_binding
    implicit none

    type(PDM_dbbtree_t), pointer        :: dbbt          
    integer, intent(in)                 :: n_part
    integer(pdm_l_num_s), pointer       :: nElts(:)
    type (PDM_pointer_array_t), pointer :: extents
    type (PDM_pointer_array_t), pointer :: gNum

    interface 
      function PDM_dbbtree_boxes_set_c (dbbt,    &
                                        n_part,  &
                                        nElts,   &
                                        extents, &
                                        gNum     &
                                        )        &
      result (bs)                                &
      bind (c, name='PDM_dbbtree_boxes_set')
  
        use iso_c_binding
        implicit none
        type (c_ptr), value :: dbbt
        integer(c_int)     :: n_part
        type (c_ptr), value :: nElts
        type (c_ptr), value :: extents
        type (c_ptr), value :: gNum
      
        type (c_ptr)        :: bs
      end function PDM_dbbtree_boxes_set_c
    end interface

    if (.not. associated(dbbt)) then
      print*, "Error PDM_dbbtree_create : dbb is already associated ! "
      call exit
    endif

    dbbt%bs_boxes_in_tree = PDM_dbbtree_boxes_set_c (dbbt%c_dbbtree,      &
                                                     n_part,              &
                                                     c_loc(nElts),        &  
                                                     c_loc(extents%cptr), &
                                                     c_loc(gNum%cptr))
  end subroutine PDM_dbbtree_boxes_set_

  !>
  !! \brief Assign a set of boxes to an empty \ref PDM_dbbtree_t structure.
  !!
  !! This function assigns a set of boxes to an empty \ref PDM_dbbtree_t structure.
  !!
  !! \param [in]  dbbt           Pointer to a distributed bounding box tree
  !! \param [in]  n_part         Number of partitions
  !! \param [in]  nElts          Number of elements of each partition
  !! \param [in]  init_location  Init location of each element of each partition (triplet rank/ipart/ielt)
  !! \param [in]  extents        Extents of each element of each partition
  !! \param [in]  gNum           Global number of each element of each partition
  !!
  !! \return associated \ref PDM_box_set_t structure distributed according to
  !! the tree location
  !!
  !!

  subroutine PDM_dbbtree_boxes_set_with_init_location (dbbt,    &
                                                       n_part,  &
                                                       nElts,   &
                                                       init_location, &                                                         
                                                       extents, &
                                                       gNum     &
                                                       )        

    use iso_c_binding
    implicit none

    type(PDM_dbbtree_t), pointer        :: dbbt          
    integer, intent(in)                 :: n_part
    integer(pdm_l_num_s), pointer       :: nElts(:)
    type (PDM_pointer_array_t), pointer :: init_location
    type (PDM_pointer_array_t), pointer :: extents
    type (PDM_pointer_array_t), pointer :: gNum
      

    interface 
      function PDM_dbbtree_boxes_set_with_init_location_c (dbbt,          &
                                                           n_part,        &
                                                           nElts,         &
                                                           init_location, &                                                         
                                                           extents,       &
                                                           gNum           &
                                                           )              &
      result (bs)                                                         &
      bind (c, name='PDM_dbbtree_boxes_set_with_init_location')
  
        use iso_c_binding
        implicit none

        type (c_ptr), value :: dbbt
        integer(c_int)      :: n_part
        type (c_ptr), value :: nElts
        type (c_ptr), value :: init_location
        type (c_ptr), value :: extents
        type (c_ptr), value :: gNum
      
        type (c_ptr)        :: bs

      end function PDM_dbbtree_boxes_set_with_init_location_c

    end interface

    if (.not. associated(dbbt)) then
      print*, "Error PDM_dbbtree_create : dbb is already associated ! "
      call exit 
    endif

    dbbt%bs_boxes_in_tree = PDM_dbbtree_boxes_set_with_init_location_c (dbbt%c_dbbtree,            &
                                                                        n_part,                    &
                                                                        c_loc(nElts),              &
                                                                        c_loc(init_location%cptr), &
                                                                        c_loc(extents%cptr),       &
                                                                        c_loc(gNum%cptr))

  end subroutine PDM_dbbtree_boxes_set_with_init_location


  !>
  !! \brief Assign boxes to intersect to the tree.
  !!
  !! This function  assigns boxes to intersect to the tree.
  !!
  !! \param [in]  dbbt               Pointer to a distributed bounding box tree
  !! \param [in]  n_part             Number of partitions
  !! \param [in]  nElts              Number of elements of each partition
  !! \param [in]  extents            Extents of each element of each partition
  !! \param [in]  gNum               Global number of each element of each partition
  !! \param [out] C_ARRAY_box_index  Index array on associated tree bounding boxes (Call pdm_fortran_free_c to free it)
  !! \param [out] C_ARRAY_box_g_num  Array of intersecting bounding boxes (Call pdm_fortran_free_c to free it)
  !!
  !! \return associated \ref PDM_box_set_t structure distributed according
  !! to the tree intersection
  !!
  !!/

  function PDM_dbbtree_intersect_boxes_set_ (dbbt,             &
                                             n_part,           &
                                             nElts,            &
                                             extents,          &
                                             gNum,             &
                                             C_ARRAY_box_index,&
                                             C_ARRAY_box_l_num &
                                             )                 &
  result (bs) 

    use iso_c_binding
    implicit none

    type(PDM_dbbtree_t), pointer        :: dbbt          
    integer, intent(in)                 :: n_part
    integer(pdm_l_num_s), pointer       :: nElts(:)
    type (PDM_pointer_array_t), pointer :: extents
    type (PDM_pointer_array_t), pointer :: gNum
    integer(pdm_l_num_s), pointer       :: C_ARRAY_box_index(:)
    integer(pdm_l_num_s), pointer       :: C_ARRAY_box_l_num(:)
      
    type (c_ptr)                        :: bs

    type(c_ptr)                         :: c_box_index
    type(c_ptr)                         :: c_box_l_num

    interface 
      function PDM_dbbtree_intersect_boxes_set_c (dbbt,     &
                                                  n_part,   &
                                                  nElts,    &
                                                  extents,  &
                                                  gNum,     &
                                                  box_index,&
                                                  box_l_num &
                                                  )         &
      result (bs)                                           &
      bind (c, name='PDM_dbbtree_intersect_boxes_set')
  
        use iso_c_binding
        implicit none
        type (c_ptr), value :: dbbt
        integer(c_int)     :: n_part
        type (c_ptr), value :: nElts
        type (c_ptr), value :: extents
        type (c_ptr), value :: gNum
        type (c_ptr)        :: box_index    
        type (c_ptr)        :: box_l_num    

        type (c_ptr)        :: bs
      end function PDM_dbbtree_intersect_boxes_set_c
    end interface

    if (.not. associated(dbbt)) then
      print*, "Error PDM_dbbtree_create : dbb is already associated ! "
      call exit 
    endif

    dbbt%bs_intersected_boxes = PDM_dbbtree_intersect_boxes_set (dbbt%c_dbbtree,      &
                                                                 n_part,              &
                                                                 c_loc(nElts),        &  
                                                                 c_loc(extents%cptr), &
                                                                 c_loc(gNum%cptr),    &
                                                                 c_box_index,         &
                                                                 c_box_l_num)

!    c_f_pointer 

  end function PDM_dbbtree_intersect_boxes_set_

! PDM_box_set_t *
! PDM_dbbtree_intersect_boxes_set
! (
!  PDM_dbbtree_t      *dbbt,
!  const int           n_part,
!  const int          *nElts,
!  const double      **extents,
!  const PDM_g_num_t **gNum,
!  int                *box_index[],
!  int                *box_l_num[]
!  );

  !>
  !! \brief Assign boxes to intersect to the tree.
  !!
  !! This function  assigns boxes to intersect to the tree.
  !!
  !! \param [in]  dbbt               Pointer to a distributed bounding box tree
  !! \param [in]  n_part             Number of partitions
  !! \param [in]  nElts              Number of elements of each partition
  !! \param [in]  extents            Extents of each element of each partition
  !! \param [in]  init_location      Init location of each element of each partition (triplet rank/ipart/ielt)
  !! \param [in]  gNum               Global number of each element of each partition
  !! \param [out] C_ARRAY_box_index  Index array on associated tree bounding boxes (Call pdm_fortran_free_c to free it)
  !! \param [out] C_ARRAY_box_g_num  Array of intersecting bounding boxes (Call pdm_fortran_free_c to free it)
  !!
  !! \return associated \ref PDM_box_set_t structure distributed according
  !! to the tree intersection
  !!
  !!/



! PDM_box_set_t *
! PDM_dbbtree_intersect_boxes_with_init_location_set
! (
!  PDM_dbbtree_t      *dbbt,
!  const int           n_part,
!  const int          *nElts,
!  const int         **init_location,
!  const double      **extents,
!  const PDM_g_num_t **gNum,
!  int                *box_index[],
!  int                *box_l_num[]
! );

end module pdm_dbbtree