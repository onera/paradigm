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

  interface PDM_dbbtree_intersect_boxes_set
     module procedure PDM_dbbtree_intersect_boxes_set_
     module procedure PDM_dbbtree_intersect_boxes_with_init_location_set
  end interface PDM_dbbtree_intersect_boxes_set

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
	!!
    
  function PDM_dbbtree_create_(comm,          &
                               dim,           & 
                               global_extents)&
    result (dbbt)   

    integer,          intent(in)  :: dim
    integer,          intent(in)  :: comm
    double precision, intent(in)  :: global_extents

    type(c_ptr)                   :: dbbt          

    integer(c_int)                :: c_comm

    use iso_c_binding
	  implicit none

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

    c_comm = PDM_MPI_Comm_f2c(comm)

    dbbt = PDM_dbbtree_create_c (c_comm, dim, c_loc(global_extents))

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

  function PDM_dbbtree_boxes_set_ (dbbt,    &
                                   n_part,  &
                                   nElts,   &
                                   extents, &
                                   gNum     &
                                   )        &
  result (bs) 

    use iso_c_binding
	  implicit none

    type (cptr), intent(in)             :: dbbt
    integer, intent(in)                 :: n_part
    integer(pdm_l_num_s), pointer       :: nElts(:)
    type (PDM_pointer_array_t), pointer :: extents
    type (PDM_pointer_array_t), pointer :: gNum
      
    type (cptr)                         :: bs

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

        type (cptr), value :: dbbt
        integer(c_int)     :: n_part
        type (cptr), value :: nElts
        type (cptr), value :: extents
        type (cptr), value :: gNum
      
        type (cptr)        :: bs

      end function PDM_dbbtree_boxes_set_c

    end interface

    bs = PDM_dbbtree_boxes_set_c (dbbt, n_part, c_loc(nElts), )

  end function PDM_dbbtree_boxes_set_

PDM_box_set_t *
PDM_dbbtree_boxes_set
(
 PDM_dbbtree_t      *dbbt,
 const int           n_part,
 const int          *nElts,
 const double      **extents,
 const PDM_g_num_t **gNum
);


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

PDM_box_set_t *
PDM_dbbtree_boxes_set_with_init_location
(
 PDM_dbbtree_t      *dbbt,
 const int           n_part,
 const int          *nElts,
 const int         **init_location,
 const double      **extents,
 const PDM_g_num_t **gNum
);


  !>
  !! \brief Assign boxes to intersect to the tree.
  !!
  !! This function  assigns boxes to intersect to the tree.
  !!
  !! \param [in]  dbbt       Pointer to a distributed bounding box tree
  !! \param [in]  n_part     Number of partitions
  !! \param [in]  nElts      Number of elements of each partition
  !! \param [in]  extents    Extents of each element of each partition
  !! \param [in]  gNum       Global number of each element of each partition
  !! \param [out] box_index  Pointer to the index array on associated tree bounding boxeq
  !! \param [out] box_g_num  Pointer to the list of intersecting bounding boxes
  !!
  !! \return associated \ref PDM_box_set_t structure distributed according
  !! to the tree intersection
  !!
  !!/


PDM_box_set_t *
PDM_dbbtree_intersect_boxes_with_init_location_set
(
 PDM_dbbtree_t      *dbbt,
 const int           n_part,
 const int          *nElts,
 const int         **init_location,
 const double      **extents,
 const PDM_g_num_t **gNum,
 int                *box_index[],
 int                *box_l_num[]
);


  !>
  !! \brief Assign boxes to intersect to the tree.
  !!
  !! This function  assigns boxes to intersect to the tree.
  !!
  !! \param [in]  dbbt           Pointer to a distributed bounding box tree
  !! \param [in]  n_part         Number of partitions
  !! \param [in]  nElts          Number of elements of each partition
  !! \param [in]  extents        Extents of each element of each partition
  !! \param [in]  init_location  Init location of each element of each partition (triplet rank/ipart/ielt)
  !! \param [in]  gNum           Global number of each element of each partition
  !! \param [out] box_index      Pointer to the index array on associated tree bounding boxeq
  !! \param [out] box_g_num      Pointer to the list of intersecting bounding boxes
  !!
  !! \return associated \ref PDM_box_set_t structure distributed according
  !! to the tree intersection
  !!
  !!/

PDM_box_set_t *
PDM_dbbtree_intersect_boxes_set
(
 PDM_dbbtree_t      *dbbt,
 const int           n_part,
 const int          *nElts,
 const double      **extents,
 const PDM_g_num_t **gNum,
 int                *box_index[],
 int                *box_l_num[]
 );

end module pdm_dbbtree