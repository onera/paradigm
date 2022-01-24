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

#include "pdm_configf.h"

module pdm_part

  use pdm

  implicit none

  integer, parameter :: PDM_part_SPLIT_PARMETIS = 1
  integer, parameter :: PDM_part_SPLIT_PTSCOTCH = 2
  integer, parameter :: PDM_part_SPLIT_HILBERT  = 3

  integer, parameter :: PDM_PART_RENUM_FACE_RANDOM        = 1
  integer, parameter :: PDM_PART_RENUM_FACE_NONE          = 2
  integer, parameter :: PDM_PART_RENUM_FACE_LEXICOGRAPHIC = 3

  integer, parameter :: PDM_part_RENUM_CELL_HILBERT = 1
  integer, parameter :: PDM_part_RENUM_CELL_RANDOM  = 2
  integer, parameter :: PDM_part_RENUM_CELL_NONE    = 3
  integer, parameter :: PDM_part_RENUM_CELL_CUTHILL = 4

  interface pdm_part_create ; module procedure  &
    pdm_part_create_
  end interface

  interface pdm_part_part_dim_get ; module procedure  &
    pdm_part_part_dim_get_
  end interface

  interface pdm_part_part_val_get ; module procedure  &
  pdm_part_part_val_get_
  end interface

  interface pdm_part_time_get ; module procedure  &
  pdm_part_time_get_
  end interface

  interface pdm_part_stat_get ; module procedure  &
  pdm_part_stat_get_
  end interface

  interface pdm_part_coarse_mesh_create ; module procedure  &
    pdm_part_coarse_mesh_create_
  end interface

  interface pdm_part_renum_method_cell_idx_get ; module procedure  &
    pdm_part_renum_method_cell_idx_get_
  end interface

  interface pdm_part_renum_method_face_idx_get ; module procedure  &
    pdm_part_renum_method_face_idx_get_
  end interface

  interface pdm_part_renum_method_cell_name_get ; module procedure  &
    pdm_part_renum_method_cell_name_get_
  end interface

  interface pdm_part_renum_method_face_name_get ; module procedure  &
    pdm_part_renum_method_face_name_get_
  end interface

interface

!>
!!
!! \brief Build a initial partitioning
!!
!!  Build a initial partitioning from :
!!      - Cell block distribution with implicit global numbering
!!         (the first cell is the first cell of the first process and
!!          the latest cell is the latest cell of the latest process)
!!      - Face block distribution with implicit global numbering
!!      - Vertex block distribution with implicit global numbering
!!  To repart an existing partition use \ref PDM_part_repart function
!!
!! \param [in]   comm                   MPI Comminicator
!! \param [in]   split_method           Split method
!! \param [in]   renum_cell_method      Cell renumbering method
!! \param [in]   renum_face_method      Cell renumbering method
!! \param [in]   renum_properties_cell  For cache blocking [ n_cell_per_cache_wanted, isAsynchrone, isVectorisation ] \ref PDM_renum_cacheblocking
!! \param [in]   renum_face_method      Face renumbering method
!! \param [in]   renum_properties_face  NOT USED
!! \param [in]   n_part                 Number of partition to build on this process
!! \param [in]   dn_cell                Number of distributed cells
!! \param [in]   dn_face                Number of distributed faces
!! \param [in]   dn_vtx                 Number of distributed vertices
!! \param [in]   n_face_group           Number of face groups
!! \param [in]   dcell_faceIdx          Distributed cell face connectivity index or NULL
!!                                      (size : dn_cell + 1, numbering : 0 to n-1)
!! \param [in]   dcell_face             Distributed cell face connectivity or NULL
!!                                      (size : dface_vtx_idx[dn_cell], numbering : 1 to n)
!! \param [in]   dcell_tag              Cell tag (size : n_cell) or NULL
!! \param [in]   dcell_weight           Cell weight (size : n_cell) or NULL
!! \param [in]   dcell_part             Distributed cell partitioning
!!                                      (size = dn_cell) or NULL (No partitioning if != NULL)
!! \param [in]   dface_cell             Distributed face cell connectivity or NULL
!!                                      (size : 2 * dn_face, numbering : 1 to n)
!! \param [in]   dface_vtx_idx          Distributed face to vertex connectivity index
!!                                      (size : dn_face + 1, numbering : 0 to n-1)
!! \param [in]   dface_vtx              Distributed face to vertex connectivity
!!                                      (size : dface_vtx_idx[dn_face], numbering : 1 to n)
!! \param [in]   dface_tag              Distributed face tag (size : dn_face)
!!                                      or NULL
!! \param [in]   dvtx_coord             Distributed vertex coordinates
!!                                      (size : 3*dn_vtx)
!! \param [in]   dvtx_tag               Distributed vertex tag (size : dn_vtx) or NULL
!! \param [in]   dface_group_idx        Index of distributed faces list of each group
!!                                      (size = n_face_group + 1) or NULL
!! \param [in]   dface_group            Distributed faces list of each group
!!                                      (size = dface_group[dface_group_idx[n_face_group]], numbering : 1 to n)
!!                                      or NULL
!!
!! \return    Pointer to \ref PDM_part object
!!

  function pdm_part_create_cf (comm, &
                               split_method, &
                               renum_cell_method, &
                               renum_face_method, &
                               n_property_cell, &
                               renum_properties_cell, &
                               n_property_face, &
                               renum_properties_face, &
                               n_part, &
                               dn_cell, &
                               dn_face, &
                               dn_vtx, &
                               n_face_group, &
                               dcell_faceIdx, &
                               dcell_face, &
                               dcell_tag, &
                               dcell_weight, &
                               have_dcell_part, &
                               dcell_part, &
                               dface_cell, &
                               dface_vtx_idx, &
                               dface_vtx, &
                               dface_tag, &
                               dvtx_coord, &
                               dvtx_tag, &
                               dface_group_idx, &
                               dface_group) &
  result (ppart) &
  bind (c, name='PDM_part_create')

  use iso_c_binding
  implicit none

  type(c_ptr)            :: ppart
  integer(c_int), value  :: comm
  integer(c_int), value  :: split_method
  character(kind=c_char) :: renum_cell_method
  character(kind=c_char) :: renum_face_method
  integer(c_int), value  :: n_property_cell
  type(c_ptr),    value  :: renum_properties_cell
  integer(c_int), value  :: n_property_face
  type(c_ptr),    value  :: renum_properties_face
  integer(c_int), value  :: n_part
  integer(c_int), value  :: dn_cell
  integer(c_int), value  :: dn_face
  integer(c_int), value  :: dn_vtx
  integer(c_int), value  :: n_face_group
  type(c_ptr),    value  :: dcell_faceIdx
  type(c_ptr),    value  :: dcell_face
  type(c_ptr),    value  :: dcell_tag
  type(c_ptr),    value  :: dcell_weight
  integer(c_int), value  :: have_dcell_part
  type(c_ptr),    value  :: dcell_part
  type(c_ptr),    value  :: dface_cell
  type(c_ptr),    value  :: dface_vtx_idx
  type(c_ptr),    value  :: dface_vtx
  type(c_ptr),    value  :: dface_tag
  type(c_ptr),    value  :: dvtx_coord
  type(c_ptr),    value  :: dvtx_tag
  type(c_ptr),    value  :: dface_group_idx
  type(c_ptr),    value  :: dface_group

  end function pdm_part_create_cf


 !>
 !!
 !! \brief Return a mesh partition dimensions
 !!
 !! \param [in]   ppart               Pointer to \ref PDM_part object
 !! \param [in]   i_part              Current partition
 !! \param [out]  n_cell              Number of cells
 !! \param [out]  n_face              Number of faces
 !! \param [out]  n_face_part_bound   Number of partitioning boundary faces
 !! \param [out]  n_vtx               Number of vertices
 !! \param [out]  n_proc              Number of processus
 !! \param [out]  n_total_part        Number of partitions
 !! \param [out]  scell_face          Size of cell-face connectivity
 !! \param [out]  sface_vtx           Size of face-vertex connectivity
 !! \param [out]  sFacePartBound      Size of face_part_bound array
 !! \param [out]  sface_group         Size of face_group array
 !! \param [out]  n_face_group        Number of face groups
 !!

 subroutine pdm_part_part_dim_get_cf (ppart, &
                                      i_part, &
                                      n_cell, &
                                      n_face, &
                                      n_face_part_bound, &
                                      n_vtx, &
                                      n_proc, &
                                      n_total_part, &
                                      scell_face, &
                                      sface_vtx, &
                                      sface_group, &
                                      n_face_group) &
  bind (c, name="PDM_part_part_dim_get")
     use pdm
     use iso_c_binding

     implicit none

     type(c_ptr),    value :: ppart
     integer(c_int), value :: i_part
     integer(c_int)         :: n_cell
     integer(c_int)         :: n_face
     integer(c_int)         :: n_face_part_bound
     integer(c_int)         :: n_vtx
     integer(c_int)         :: n_proc
     integer(c_int)         :: n_total_part
     integer(c_int)         :: scell_face
     integer(c_int)         :: sface_vtx
     integer(c_int)         :: sface_group
     integer(c_int)         :: n_face_group

   end subroutine pdm_part_part_dim_get_cf


 !>
 !!
 !! \brief Return a mesh partition
 !!
 !! \param [in]   ppart                     Pointer to \ref PDM_part object
 !! \param [in]   i_part                    Current partition
 !! \param [out]  cell_tag                  Cell tag (size = n_cell)
 !! \param [out]  cell_face_idx             Cell to face connectivity index (size = n_cell + 1, numbering : 0 to n-1)
 !! \param [out]  cell_face                 Cell to face connectivity (size = cell_face_idx[n_cell] = lcell_face
 !!                                                                   numbering : 1 to n)
 !! \param [out]  cell_ln_to_gn             Cell local numbering to global numbering (size = n_cell, numbering : 1 to n)
 !! \param [out]  face_tag                  Face tag (size = n_face)
 !! \param [out]  face_cell                 Face to cell connectivity  (size = 2 * n_face, numbering : 1 to n)
 !! \param [out]  face_vtx_idx              Face to Vertex connectivity index (size = n_face + 1, numbering : 0 to n-1)
 !! \param [out]  face_vtx                  Face to Vertex connectivity (size = faceVertexIdx[n_face], numbering : 1 to n)
 !! \param [out]  face_ln_to_gn             Face local numbering to global numbering (size = n_face, numbering : 1 to n)
 !! \param [out]  face_part_bound_proc_idx  Partitioning boundary faces block distribution from processus (size = n_proc + 1)
 !! \param [out]  face_part_bound_part_idx  Partitioning boundary faces block distribution from partition (size = n_total_part + 1)
 !! \param [out]  face_part_bound           Partitioning boundary faces (size = 4 * n_face_part_bound)
 !!                                          sorted by processus, sorted by partition in each processus, and
 !!                                          sorted by absolute face number in each partition
 !!                                         For each face :
 !!                                           - Face local number (numbering : 1 to n)
 !!                                           - Connected process (numbering : 0 to n-1)
 !!                                           - Connected Partition
 !!                                             on the connected process (numbering :1 to n)
 !!                                           - Connected face local number
 !!                                             in the connected partition (numbering :1 to n)
 !! \param [out]  vtx_tag                   Vertex tag (size = nVertex)
 !! \param [out]  vtx                       Vertex coordinates (size = 3 * nVertex)
 !! \param [out]  vtx_ln_to_gn              Vertex local numbering to global numbering (size = n_vtx, numbering : 1 to n)
 !! \param [out]  face_group_idx            Face group index (size = n_face_group + 1, numbering : 1 to n-1)
 !! \param [out]  face_group                Faces for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
 !! \param [out]  face_group_ln_to_gn       Faces global numbering for each group
 !!                                         (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
 !!

 subroutine pdm_part_part_val_get_cf (ppart, &
                                      i_part, &
                                      cell_tag, &
                                      cell_face_idx, &
                                      cell_face, &
                                      cell_ln_to_gn, &
                                      face_tag, &
                                      face_cell, &
                                      face_vtx_idx, &
                                      face_vtx, &
                                      face_ln_to_gn, &
                                      face_part_bound_proc_idx, &
                                      face_part_bound_part_idx, &
                                      face_part_bound, &
                                      vtx_tag, &
                                      vtx, &
                                      vtx_ln_to_gn, &
                                      face_group_idx, &
                                      face_group, &
                                      face_group_ln_to_gn) &
 bind (c, name='PDM_part_part_val_get')

   use pdm
   use iso_c_binding

   implicit none

   type(c_ptr),    value :: ppart
   integer(c_int), value :: i_part
   type(c_ptr)           :: cell_tag
   type(c_ptr)           :: cell_face_idx
   type(c_ptr)           :: cell_face
   type(c_ptr)           :: cell_ln_to_gn
   type(c_ptr)           :: face_tag
   type(c_ptr)           :: face_cell
   type(c_ptr)           :: face_vtx_idx
   type(c_ptr)           :: face_vtx
   type(c_ptr)           :: face_ln_to_gn
   type(c_ptr)           :: face_part_bound_proc_idx
   type(c_ptr)           :: face_part_bound_part_idx
   type(c_ptr)           :: face_part_bound
   type(c_ptr)           :: vtx_tag
   type(c_ptr)           :: vtx
   type(c_ptr)           :: vtx_ln_to_gn
   type(c_ptr)           :: face_group_idx
   type(c_ptr)           :: face_group
   type(c_ptr)           :: face_group_ln_to_gn

 end subroutine pdm_part_part_val_get_cf

 !>
 !!
 !! \brief Free ppart
 !!
 !! \param [in]   ppart               Pointer to \ref PDM_part object
 !!
 !!

 subroutine pdm_part_free (ppart) &
 bind (c, name='PDM_part_free')

   use iso_c_binding
   implicit none

   type(c_ptr), value :: ppart

 end subroutine pdm_part_free

 !>
 !!
 !! \brief Return times
 !!
 !! \param [in]   ppart       Pointer to \ref PDM_part object
 !! \param [out]  elapsed     elapsed times (size = 4)
 !! \param [out]  cpu         cpu times (size = 4)
 !! \param [out]  cpu_user    user cpu times (size = 4)
 !! \param [out]  cpu_sys     system cpu times (size = 4)
 !!
 !!

 subroutine pdm_part_time_get_cf (ppart,    &
                                  elapsed,  &
                                  cpu,      &
                                  cpu_user, &
                                  cpu_sys)  &
 bind (c, name='PDM_part_time_get')

   use iso_c_binding
   implicit none

   type(c_ptr), value :: ppart
   type(c_ptr)        :: elapsed
   type(c_ptr)        :: cpu
   type(c_ptr)        :: cpu_user
   type(c_ptr)        :: cpu_sys

 end subroutine pdm_part_time_get_cf

 !>
 !!
 !! \brief Return statistic
 !!
 !! \param [in]   ppart                          Pointer to \ref PDM_part object
 !! \param [out]  cells_average                  average of cells number
 !! \param [out]  cells_median                   median of cells number
 !! \param [out]  cells_std_deviation            standard deviation of cells number
 !! \param [out]  cells_min                      minimum of cells nummber
 !! \param [out]  cells_max                      maximum of cells nummber
 !! \param [out]  bound_part_faces_average       average of partitioning boundary faces
 !! \param [out]  bound_part_faces_median        median of partitioning boundary faces
 !! \param [out]  bound_part_faces_std_deviation standard deviation of partitioning boundary faces
 !! \param [out]  bound_part_faces_min           minimum of partitioning boundary faces
 !! \param [out]  bound_part_faces_max           maximum of partitioning boundary faces
 !!
 !!

 subroutine pdm_part_stat_get_cf (ppart,                          &
                                  cells_average,                  &
                                  cells_median,                   &
                                  cells_std_deviation,            &
                                  cells_min,                      &
                                  cells_max,                      &
                                  bound_part_faces_average,       &
                                  bound_part_faces_median,        &
                                  bound_part_faces_std_deviation, &
                                  bound_part_faces_min,           &
                                  bound_part_faces_max,           &
                                  bound_part_faces_sum)           &
 bind (c, name='PDM_part_stat_get')

   use pdm
   use iso_c_binding
   implicit none

   type(c_ptr), value :: ppart
   integer(c_int)     :: cells_average
   integer(c_int)     :: cells_median
   real(c_double)     :: cells_std_deviation
   integer(c_int)     :: cells_min
   integer(c_int)     :: cells_max
   integer(c_int)     :: bound_part_faces_average
   integer(c_int)     :: bound_part_faces_median
   real(c_double)     :: bound_part_faces_std_deviation
   integer(c_int)     :: bound_part_faces_min
   integer(c_int)     :: bound_part_faces_max
   integer(c_int)     :: bound_part_faces_sum

 end subroutine pdm_part_stat_get_cf


 !================================================================================
 !
 ! \brief Get the number of renumbering cell methods
 !
 ! \param [out]    Number of methods
 !
 !================================================================================

 subroutine pdm_part_n_renum_method_cell_get (n_method)
   use pdm
   implicit none
   integer      :: n_method

 end subroutine pdm_part_n_renum_method_cell_get



 !================================================================================
 !
 ! \brief Get the number of renumbering cell methods
 !
 ! \param [out]    Number of methods
 !
 !================================================================================

 subroutine pdm_part_n_renum_method_face_get (n_method)
   use pdm
   implicit none
   integer      :: n_method

 end subroutine pdm_part_n_renum_method_face_get

end interface

private :: pdm_part_create_ ,&
           pdm_part_part_dim_get_ ,&
           pdm_part_part_val_get_ ,&
           pdm_part_time_get_ ,&
           pdm_part_stat_get_ ,&
           pdm_part_coarse_mesh_create_ ,&
           pdm_part_renum_method_cell_idx_get_ ,&
           pdm_part_renum_method_face_idx_get_ ,&
           pdm_part_renum_method_cell_name_get_ ,&
           pdm_part_renum_method_face_name_get_

contains

!>
!!
!! \brief Build a initial partitioning
!!
!!  Build a initial partitioning from :
!!      - Cell block distribution with implicit global numbering
!!         (the first cell is the first cell of the first process and
!!          the latest cell is the latest cell of the latest process)
!!      - Face block distribution with implicit global numbering
!!      - Vertex block distribution with implicit global numbering
!!  To repart an existing partition use \ref PDM_part_repart function
!!
!! \param [out]  ppart                  Pointer to \ref PDM_part object
!! \param [in]   f_comm                 MPI Comminicator
!! \param [in]   split_method           Split method
!! \param [in]   renum_cell_method      Cell renumbering method
!! \param [in]   renum_face_method      Cell renumbering method
!! \param [in]   renum_properties_cell  For cache blocking [ n_cell_per_cache_wanted, isAsynchrone, isVectorisation ] \ref PDM_renum_cacheblocking
!! \param [in]   renum_face_method      Face renumbering method
!! \param [in]   renum_properties_face  NOT USED
!! \param [in]   n_part                 Number of partition to build on this process
!! \param [in]   dn_cell                Number of distributed cells
!! \param [in]   dn_face                Number of distributed faces
!! \param [in]   dn_vtx                 Number of distributed vertices
!! \param [in]   n_face_group           Number of face groups
!! \param [in]   dcell_faceIdx          Distributed cell face connectivity index or NULL
!!                                      (size : dn_cell + 1, numbering : 0 to n-1)
!! \param [in]   dcell_face             Distributed cell face connectivity or NULL
!!                                      (size : dface_vtx_idx[dn_cell], numbering : 1 to n)
!! \param [in]   dcell_tag              Cell tag (size : n_cell) or NULL
!! \param [in]   dcell_weight           Cell weight (size : n_cell) or NULL
!! \param [in]   dcell_part             Distributed cell partitioning
!!                                      (size = dn_cell) or NULL (No partitioning if != NULL)
!! \param [in]   dface_cell             Distributed face cell connectivity or NULL
!!                                      (size : 2 * dn_face, numbering : 1 to n)
!! \param [in]   dface_vtx_idx          Distributed face to vertex connectivity index
!!                                      (size : dn_face + 1, numbering : 0 to n-1)
!! \param [in]   dface_vtx              Distributed face to vertex connectivity
!!                                      (size : dface_vtx_idx[dn_face], numbering : 1 to n)
!! \param [in]   dface_tag              Distributed face tag (size : dn_face)
!!                                      or NULL
!! \param [in]   dvtx_coord             Distributed vertex coordinates
!!                                      (size : 3*dn_vtx)
!! \param [in]   dvtx_tag               Distributed vertex tag (size : dn_vtx) or NULL
!! \param [in]   dface_group_idx        Index of distributed faces list of each group
!!                                      (size = n_face_group + 1) or NULL
!! \param [in]   dface_group            Distributed faces list of each group
!!                                      (size = dface_group[dface_group_idx[n_face_group]], numbering : 1 to n)
!!                                      or NULL
!!
!!
 subroutine pdm_part_create_(ppart, &
                             f_comm, &
                             split_method,  &
                             renum_cell_method, &
                             renum_face_method, &
                             nPropertyCell, &
                             renum_properties_cell, &
                             nPropertyFace, &
                             renum_properties_face, &
                             nPart, &
                             dNCell, &
                             dNFace, &
                             dNVtx,&
                             nFaceGroup, &
                             have_dCellFace,&
                             dCellFaceIdx,&
                             dCellFace,&
                             have_dCellTag,&
                             dCellTag, &
                             have_dCellWeight,&
                             dCellWeight,&
                             have_dCellPart,&
                             dCellPart,&
                             have_dFaceCell,&
                             dFaceCell,&
                             dFaceVtxIdx,&
                             dFaceVtx,&
                             have_dFaceTag,&
                             dFaceTag,&
                             dVtxCoord,&
                             have_dVtxTag,&
                             dVtxTag,&
                             dFaceGroupIdx,&
                             dFaceGroup)

    use pdm
    use iso_c_binding

    implicit none

    type(c_ptr)                        :: ppart
    integer                            :: f_comm
    integer                            :: split_method
    character (len=*)                  :: renum_cell_method
    character (len=*)                  :: renum_face_method
    integer                            :: nPropertyCell
    integer(kind=PDM_l_num_s), pointer :: renum_properties_cell(:)
    integer                            :: nPropertyFace
    integer(kind=PDM_l_num_s), pointer :: renum_properties_face(:)
    integer                            :: nPart
    integer                            :: dNCell
    integer                            :: dNFace
    integer                            :: dNVtx
    integer                            :: nFaceGroup
    logical                            :: have_dCellFace
    integer(kind=PDM_l_num_s), pointer :: dCellFaceIdx(:)
    integer(kind=PDM_g_num_s), pointer :: dCellFace(:)
    logical                            :: have_dCellTag
    integer(kind=PDM_l_num_s), pointer :: dCellTag(:)
    logical                            :: have_dCellWeight
    integer(kind=PDM_l_num_s), pointer :: dCellWeight(:)
    logical                            :: have_dCellPart
    integer(kind=PDM_l_num_s), pointer :: dCellPart(:)
    logical                            :: have_dFaceCell
    integer(kind=PDM_g_num_s), pointer :: dFaceCell(:)
    integer(kind=PDM_l_num_s), pointer :: dFaceVtxIdx(:)
    integer(kind=PDM_g_num_s), pointer :: dFaceVtx(:)
    logical                            :: have_dFaceTag
    integer(kind=PDM_l_num_s), pointer :: dFaceTag(:)
    double precision,          pointer :: dVtxCoord(:)
    logical                            :: have_dVtxTag
    integer(kind=PDM_l_num_s), pointer :: dVtxTag(:)
    integer(kind=PDM_l_num_s), pointer :: dFaceGroupIdx(:)
    integer(kind=PDM_g_num_s), pointer :: dFaceGroup(:)

    integer(kind=c_int)                :: c_split_method
    integer(kind=c_int)                :: c_n_property_cell
    type(c_ptr)                        :: c_renum_properties_cell = C_NULL_PTR
    integer(kind=c_int)                :: c_n_property_face
    type(c_ptr)                        :: c_renum_properties_face = C_NULL_PTR
    integer(kind=c_int)                :: c_n_part
    integer(kind=c_int)                :: c_dn_cell
    integer(kind=c_int)                :: c_dn_face
    integer(kind=c_int)                :: c_dn_vtx
    integer(kind=c_int)                :: c_n_face_group
    type(c_ptr)                        :: c_dcell_faceIdx   = C_NULL_PTR
    type(c_ptr)                        :: c_dcell_face      = C_NULL_PTR
    type(c_ptr)                        :: c_dcell_tag       = C_NULL_PTR
    type(c_ptr)                        :: c_dcell_weight    = C_NULL_PTR
    integer(c_int)                     :: c_have_dcell_part
    type(c_ptr)                        :: c_dcell_part      = C_NULL_PTR
    type(c_ptr)                        :: c_dface_cell      = C_NULL_PTR
    type(c_ptr)                        :: c_dface_vtx_idx   = C_NULL_PTR
    type(c_ptr)                        :: c_dface_vtx       = C_NULL_PTR
    type(c_ptr)                        :: c_dface_tag       = C_NULL_PTR
    type(c_ptr)                        :: c_dvtx_coord      = C_NULL_PTR
    type(c_ptr)                        :: c_dvtx_tag        = C_NULL_PTR
    type(c_ptr)                        :: c_dface_group_idx = C_NULL_PTR
    type(c_ptr)                        :: c_dface_group     = C_NULL_PTR
    integer(c_int)                     :: c_comm

    if (have_dCellPart) then
      c_have_dcell_part = 1
    else
      c_have_dcell_part = 0
    endif

    if (have_dCellTag) then
      c_dcell_tag = c_loc(dCellTag)
    endif

    if (.true.) then!have_dCellPart) then
      c_dcell_part = c_loc(dCellPart)
    endif

    if (have_dFaceCell) then
      c_dface_cell = c_loc(dFaceCell)
    endif

    if (have_dFaceTag) then
      c_dface_tag = c_loc(dFaceTag)
    endif

    if (have_dVtxTag) then
      c_dvtx_tag = c_loc(dVtxTag)
    endif

    c_split_method    = split_method
    c_n_property_cell = nPropertyCell
    c_n_property_face = nPropertyFace
    c_n_part          = nPart
    c_dn_cell         = dnCell
    c_dn_face         = dnFace
    c_dn_vtx          = dnVtx
    c_n_face_group    = nFaceGroup

    c_renum_properties_cell = c_loc(renum_properties_cell)
    c_renum_properties_face = c_loc(renum_properties_face)

    c_dcell_faceIdx   = c_loc(dCellFaceIdx )
    c_dcell_face      = c_loc(dCellFace    )
    c_dcell_weight    = c_loc(dCellWeight  )
    c_dface_vtx_idx   = c_loc(dFaceVtxIdx  )
    c_dface_vtx       = c_loc(dFaceVtx     )
    c_dvtx_coord      = c_loc(dVtxCoord    )
    c_dface_group_idx = c_loc(dFaceGroupIdx)
    c_dface_group     = c_loc(dFaceGroup   )

    c_comm = PDM_MPI_Comm_f2c(f_comm)

    ppart = pdm_part_create_cf (c_comm, &
                                c_split_method, &
                                renum_cell_method//C_NULL_CHAR, &
                                renum_face_method//C_NULL_CHAR, &
                                c_n_property_cell, &
                                c_renum_properties_cell, &
                                c_n_property_face, &
                                c_renum_properties_face, &
                                c_n_part, &
                                c_dn_cell, &
                                c_dn_face, &
                                c_dn_vtx, &
                                c_n_face_group, &
                                c_dcell_faceIdx, &
                                c_dcell_face, &
                                c_dcell_tag, &
                                c_dcell_weight, &
                                c_have_dcell_part, &
                                c_dcell_part, &
                                c_dface_cell, &
                                c_dface_vtx_idx, &
                                c_dface_vtx, &
                                c_dface_tag, &
                                c_dvtx_coord, &
                                c_dvtx_tag, &
                                c_dface_group_idx, &
                                c_dface_group)

  end subroutine pdm_part_create_



  !>
 !!
 !! \brief Return a mesh partition dimensions
 !!
 !! \param [in]   ppart               Pointer to \ref PDM_part object
 !! \param [in]   i_part              Current partition
 !! \param [out]  n_cell              Number of cells
 !! \param [out]  n_face              Number of faces
 !! \param [out]  n_face_part_bound   Number of partitioning boundary faces
 !! \param [out]  n_vtx               Number of vertices
 !! \param [out]  n_proc              Number of processus
 !! \param [out]  n_total_part        Number of partitions
 !! \param [out]  scell_face          Size of cell-face connectivity
 !! \param [out]  sface_vtx           Size of face-vertex connectivity
 !! \param [out]  sFacePartBound      Size of face_part_bound array
 !! \param [out]  sface_group         Size of face_group array
 !! \param [out]  n_face_group        Number of face groups
 !!

 subroutine pdm_part_part_dim_get_ (ppart, &
                                    i_part, &
                                    n_cell, &
                                    n_face, &
                                    n_face_part_bound, &
                                    n_vtx, &
                                    n_proc, &
                                    n_total_part, &
                                    scell_face, &
                                    sface_vtx, &
                                    sface_group, &
                                    n_face_group)
   use pdm
   use iso_c_binding

   implicit none

   type(c_ptr), value :: ppart
   integer            :: i_part
   integer            :: n_cell
   integer            :: n_face
   integer            :: n_face_part_bound
   integer            :: n_vtx
   integer            :: n_proc
   integer            :: n_total_part
   integer            :: scell_face
   integer            :: sface_vtx
   integer            :: sface_group
   integer            :: n_face_group

   integer(c_int)     :: c_i_part
   integer(c_int)     :: c_n_cell
   integer(c_int)     :: c_n_face
   integer(c_int)     :: c_n_face_part_bound
   integer(c_int)     :: c_n_vtx
   integer(c_int)     :: c_n_proc
   integer(c_int)     :: c_n_total_part
   integer(c_int)     :: c_scell_face
   integer(c_int)     :: c_sface_vtx
   integer(c_int)     :: c_sface_group
   integer(c_int)     :: c_n_face_group

   c_i_part = i_part

   call pdm_part_part_dim_get_cf(ppart, &
                                 c_i_part, &
                                 c_n_cell, &
                                 c_n_face, &
                                 c_n_face_part_bound, &
                                 c_n_vtx, &
                                 c_n_proc, &
                                 c_n_total_part, &
                                 c_scell_face, &
                                 c_sface_vtx, &
                                 c_sface_group, &
                                 c_n_face_group)

   n_cell            = c_n_cell
   n_face            = c_n_face
   n_face_part_bound = c_n_face_part_bound
   n_vtx             = c_n_vtx
   n_proc            = c_n_proc
   n_total_part      = c_n_total_part
   scell_face        = c_scell_face
   sface_vtx         = c_sface_vtx
   sface_group       = c_sface_group
   n_face_group      = c_n_face_group

 end subroutine pdm_part_part_dim_get_


 !>
 !!
 !! \brief Return a mesh partition
 !!
 !! \param [in]   ppart                     Pointer to \ref PDM_part object
 !! \param [in]   i_part                    Current partition
 !! \param [out]  cell_tag                  Cell tag (size = n_cell)
 !! \param [out]  cell_face_idx             Cell to face connectivity index (size = n_cell + 1, numbering : 0 to n-1)
 !! \param [out]  cell_face                 Cell to face connectivity (size = cell_face_idx[n_cell] = lcell_face
 !!                                                                   numbering : 1 to n)
 !! \param [out]  cell_ln_to_gn             Cell local numbering to global numbering (size = n_cell, numbering : 1 to n)
 !! \param [out]  face_tag                  Face tag (size = n_face)
 !! \param [out]  face_cell                 Face to cell connectivity  (size = 2 * n_face, numbering : 1 to n)
 !! \param [out]  face_vtx_idx              Face to Vertex connectivity index (size = n_face + 1, numbering : 0 to n-1)
 !! \param [out]  face_vtx                  Face to Vertex connectivity (size = faceVertexIdx[n_face], numbering : 1 to n)
 !! \param [out]  face_ln_to_gn             Face local numbering to global numbering (size = n_face, numbering : 1 to n)
 !! \param [out]  face_part_bound_proc_idx  Partitioning boundary faces block distribution from processus (size = n_proc + 1)
 !! \param [out]  face_part_bound_part_idx  Partitioning boundary faces block distribution from partition (size = n_total_part + 1)
 !! \param [out]  face_part_bound           Partitioning boundary faces (size = 4 * n_face_part_bound)
 !!                                          sorted by processus, sorted by partition in each processus, and
 !!                                          sorted by absolute face number in each partition
 !!                                         For each face :
 !!                                           - Face local number (numbering : 1 to n)
 !!                                           - Connected process (numbering : 0 to n-1)
 !!                                           - Connected Partition
 !!                                             on the connected process (numbering :1 to n)
 !!                                           - Connected face local number
 !!                                             in the connected partition (numbering :1 to n)
 !! \param [out]  vtx_tag                   Vertex tag (size = nVertex)
 !! \param [out]  vtx                       Vertex coordinates (size = 3 * nVertex)
 !! \param [out]  vtx_ln_to_gn              Vertex local numbering to global numbering (size = n_vtx, numbering : 1 to n)
 !! \param [out]  face_group_idx            Face group index (size = n_face_group + 1, numbering : 1 to n-1)
 !! \param [out]  face_group                Faces for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
 !! \param [out]  face_group_ln_to_gn       Faces global numbering for each group
 !!                                         (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
 !!

 subroutine pdm_part_part_val_get_ (ppart,                    &
                                    i_part,                   &
                                    cell_tag,                 &
                                    cell_face_idx,            &
                                    cell_face,                &
                                    cell_ln_to_gn,            &
                                    face_tag,                 &
                                    face_cell,                &
                                    face_vtx_idx,             &
                                    face_vtx,                 &
                                    face_ln_to_gn,            &
                                    face_part_bound_proc_idx, &
                                    face_part_bound_part_idx, &
                                    face_part_bound,          &
                                    vtx_tag,                  &
                                    vtx,                      &
                                    vtx_ln_to_gn,             &
                                    face_group_idx,           &
                                    face_group,               &
                                    face_group_ln_to_gn)

   use pdm
   use iso_c_binding

   implicit none

   type (c_ptr), value                   :: ppart
   integer                               :: i_part
   integer (kind = PDM_l_num_s), pointer :: cell_tag(:)
   integer (kind = PDM_l_num_s), pointer :: cell_face_idx(:)
   integer (kind = PDM_l_num_s), pointer :: cell_face(:)
   integer (kind = PDM_g_num_s), pointer :: cell_ln_to_gn(:)
   integer (kind = PDM_l_num_s), pointer :: face_tag(:)
   integer (kind = PDM_l_num_s), pointer :: face_cell(:)
   integer (kind = PDM_l_num_s), pointer :: face_vtx_idx(:)
   integer (kind = PDM_l_num_s), pointer :: face_vtx(:)
   integer (kind = PDM_g_num_s), pointer :: face_ln_to_gn(:)
   integer (kind = PDM_l_num_s), pointer :: face_part_bound_proc_idx(:)
   integer (kind = PDM_l_num_s), pointer :: face_part_bound_part_idx(:)
   integer (kind = PDM_l_num_s), pointer :: face_part_bound(:)
   integer (kind = PDM_l_num_s), pointer :: vtx_tag(:)
   double precision,             pointer :: vtx(:,:)
   integer (kind = PDM_g_num_s), pointer :: vtx_ln_to_gn(:)
   integer (kind = PDM_l_num_s), pointer :: face_group_idx(:)
   integer (kind = PDM_l_num_s), pointer :: face_group(:)
   integer (kind = PDM_g_num_s), pointer :: face_group_ln_to_gn(:)

   integer(c_int)                        :: c_i_part
   type(c_ptr)                           :: c_cell_tag                 = C_NULL_PTR
   type(c_ptr)                           :: c_cell_face_idx            = C_NULL_PTR
   type(c_ptr)                           :: c_cell_face                = C_NULL_PTR
   type(c_ptr)                           :: c_cell_ln_to_gn            = C_NULL_PTR
   type(c_ptr)                           :: c_face_tag                 = C_NULL_PTR
   type(c_ptr)                           :: c_face_cell                = C_NULL_PTR
   type(c_ptr)                           :: c_face_vtx_idx             = C_NULL_PTR
   type(c_ptr)                           :: c_face_vtx                 = C_NULL_PTR
   type(c_ptr)                           :: c_face_ln_to_gn            = C_NULL_PTR
   type(c_ptr)                           :: c_face_part_bound_proc_idx = C_NULL_PTR
   type(c_ptr)                           :: c_face_part_bound_part_idx = C_NULL_PTR
   type(c_ptr)                           :: c_face_part_bound          = C_NULL_PTR
   type(c_ptr)                           :: c_vtx_tag                  = C_NULL_PTR
   type(c_ptr)                           :: c_vtx                      = C_NULL_PTR
   type(c_ptr)                           :: c_vtx_ln_to_gn             = C_NULL_PTR
   type(c_ptr)                           :: c_face_group_idx           = C_NULL_PTR
   type(c_ptr)                           :: c_face_group               = C_NULL_PTR
   type(c_ptr)                           :: c_face_group_ln_to_gn      = C_NULL_PTR
   integer(c_int)                        :: n_cell
   integer(c_int)                        :: n_face
   integer(c_int)                        :: n_face_part_bound
   integer(c_int)                        :: n_vtx
   integer(c_int)                        :: n_proc
   integer(c_int)                        :: n_total_part
   integer(c_int)                        :: scell_face
   integer(c_int)                        :: sface_vtx
   integer(c_int)                        :: sface_group
   integer(c_int)                        :: n_face_group

   c_i_part = i_part

   call pdm_part_part_dim_get_cf(ppart,             &
                                 c_i_part,          &
                                 n_cell,            &
                                 n_face,            &
                                 n_face_part_bound, &
                                 n_vtx,             &
                                 n_proc,            &
                                 n_total_part,      &
                                 scell_face,        &
                                 sface_vtx,         &
                                 sface_group,       &
                                 n_face_group)

   call pdm_part_part_val_get_cf(ppart,                     &
                                 c_i_part,                   &
                                 c_cell_tag,                 &
                                 c_cell_face_idx,            &
                                 c_cell_face,                &
                                 c_cell_ln_to_gn,            &
                                 c_face_tag,                 &
                                 c_face_cell,                &
                                 c_face_vtx_idx,             &
                                 c_face_vtx,                 &
                                 c_face_ln_to_gn,            &
                                 c_face_part_bound_proc_idx, &
                                 c_face_part_bound_part_idx, &
                                 c_face_part_bound,          &
                                 c_vtx_tag,                  &
                                 c_vtx,                      &
                                 c_vtx_ln_to_gn,             &
                                 c_face_group_idx,           &
                                 c_face_group,               &
                                 c_face_group_ln_to_gn)

   call c_f_pointer(c_cell_tag, &
                    cell_tag,   &
                    [n_cell])

   call c_f_pointer(c_cell_face_idx, &
                    cell_face_idx,   &
                    [n_cell+1])

   call c_f_pointer(c_cell_face, &
                    cell_face,   &
                    [scell_face])

   call c_f_pointer(c_cell_ln_to_gn, &
                    cell_ln_to_gn,   &
                    [n_cell])

   call c_f_pointer(c_face_tag, &
                    face_tag,   &
                    [n_face])

   call c_f_pointer(c_face_cell, &
                    face_cell,   &
                    [2*n_face])

   call c_f_pointer(c_face_vtx_idx, &
                    face_vtx_idx,   &
                    [n_face+1])

   call c_f_pointer(c_face_vtx, &
                    face_vtx,   &
                    [sface_vtx])

   call c_f_pointer(c_face_ln_to_gn, &
                    face_ln_to_gn,   &
                    [sface_vtx])

   call c_f_pointer(c_face_part_bound_proc_idx, &
                    face_part_bound_proc_idx,   &
                    [n_proc+1])

   call c_f_pointer(c_face_part_bound_part_idx, &
                    face_part_bound_part_idx,   &
                    [n_proc+1])

   call c_f_pointer(c_face_part_bound, &
                    face_part_bound,   &
                    [4*n_face_part_bound])

   call c_f_pointer(c_vtx_tag, &
                    vtx_tag,   &
                    [n_vtx])

   call c_f_pointer(c_vtx, &
                    vtx,   &
                    [3, n_vtx])

   call c_f_pointer(c_vtx_ln_to_gn, &
                    vtx_ln_to_gn,   &
                    [n_vtx])

   call c_f_pointer(c_face_group_idx, &
                    face_group_idx,   &
                    [n_face_group+1])

   call c_f_pointer(c_face_group, &
                    face_group,   &
                    [sface_group])

   call c_f_pointer(c_face_group_ln_to_gn, &
                    face_group_ln_to_gn,   &
                    [sface_group])

 end subroutine pdm_part_part_val_get_



 !>
 !!
 !! \brief Return times
 !!
 !! \param [in]   ppart       Pointer to \ref PDM_part object
 !! \param [out]  elapsed     elapsed times (size = 4)
 !! \param [out]  cpu         cpu times (size = 4)
 !! \param [out]  cpu_user    user cpu times (size = 4)
 !! \param [out]  cpu_sys     system cpu times (size = 4)
 !!
 !!

 subroutine pdm_part_time_get_ (ppart,    &
                                elapsed,  &
                                cpu,      &
                                cpu_user, &
                                cpu_sys)

   use iso_c_binding
   implicit none

   type(c_ptr), value :: ppart
   double precision   :: elapsed(4)
   double precision   :: cpu(4)
   double precision   :: cpu_user(4)
   double precision   :: cpu_sys(4)

   type(c_ptr)               :: c_elapsed  = C_NULL_PTR
   type(c_ptr)               :: c_cpu      = C_NULL_PTR
   type(c_ptr)               :: c_cpu_user = C_NULL_PTR
   type(c_ptr)               :: c_cpu_sys  = C_NULL_PTR
   double precision, pointer :: ptr(:)     => null()


   call pdm_part_time_get_cf (ppart,    &
                              c_elapsed,  &
                              c_cpu,      &
                              c_cpu_user, &
                              c_cpu_sys)

   call c_f_pointer(c_elapsed, &
                    ptr,       &
                    [4])
   elapsed(1:4) = ptr(1:4)

   call c_f_pointer(c_cpu, &
                    ptr,   &
                    [4])
   cpu(1:4) = ptr(1:4)

   call c_f_pointer(c_cpu_user, &
                    ptr,        &
                    [4])
   cpu_user(1:4) = ptr(1:4)

   call c_f_pointer(c_cpu_sys, &
                    ptr,       &
                    [4])
   cpu_sys(1:4) = ptr(1:4)

 end subroutine pdm_part_time_get_



 !>
 !!
 !! \brief Return statistic
 !!
 !! \param [in]   ppart                          Pointer to \ref PDM_part object
 !! \param [out]  cells_average                  average of cells number
 !! \param [out]  cells_median                   median of cells number
 !! \param [out]  cells_std_deviation            standard deviation of cells number
 !! \param [out]  cells_min                      minimum of cells nummber
 !! \param [out]  cells_max                      maximum of cells nummber
 !! \param [out]  bound_part_faces_average       average of partitioning boundary faces
 !! \param [out]  bound_part_faces_median        median of partitioning boundary faces
 !! \param [out]  bound_part_faces_std_deviation standard deviation of partitioning boundary faces
 !! \param [out]  bound_part_faces_min           minimum of partitioning boundary faces
 !! \param [out]  bound_part_faces_max           maximum of partitioning boundary faces
 !!
 !!

 subroutine pdm_part_stat_get_ (ppart,                          &
                                cells_average,                  &
                                cells_median,                   &
                                cells_std_deviation,            &
                                cells_min,                      &
                                cells_max,                      &
                                bound_part_faces_average,       &
                                bound_part_faces_median,        &
                                bound_part_faces_std_deviation, &
                                bound_part_faces_min,           &
                                bound_part_faces_max,           &
                                bound_part_faces_sum)

   use pdm
   use iso_c_binding
   implicit none

   type(c_ptr), value        :: ppart
   integer(kind=PDM_l_num_s) :: cells_average
   integer(kind=PDM_l_num_s) :: cells_median
   double precision          :: cells_std_deviation
   integer(kind=PDM_l_num_s) :: cells_min
   integer(kind=PDM_l_num_s) :: cells_max
   integer(kind=PDM_l_num_s) :: bound_part_faces_average
   integer(kind=PDM_l_num_s) :: bound_part_faces_median
   double precision          :: bound_part_faces_std_deviation
   integer(kind=PDM_l_num_s) :: bound_part_faces_min
   integer(kind=PDM_l_num_s) :: bound_part_faces_max
   integer(kind=PDM_l_num_s) :: bound_part_faces_sum

   integer(c_int)     :: c_cells_average
   integer(c_int)     :: c_cells_median
   real(c_double)     :: c_cells_std_deviation
   integer(c_int)     :: c_cells_min
   integer(c_int)     :: c_cells_max
   integer(c_int)     :: c_bound_part_faces_average
   integer(c_int)     :: c_bound_part_faces_median
   real(c_double)     :: c_bound_part_faces_std_deviation
   integer(c_int)     :: c_bound_part_faces_min
   integer(c_int)     :: c_bound_part_faces_max
   integer(c_int)     :: c_bound_part_faces_sum

   call pdm_part_stat_get_cf(ppart,                          &
                             c_cells_average,                  &
                             c_cells_median,                   &
                             c_cells_std_deviation,            &
                             c_cells_min,                      &
                             c_cells_max,                      &
                             c_bound_part_faces_average,       &
                             c_bound_part_faces_median,        &
                             c_bound_part_faces_std_deviation, &
                             c_bound_part_faces_min,           &
                             c_bound_part_faces_max,           &
                             c_bound_part_faces_sum)

   cells_average                  = c_cells_average
   cells_median                   = c_cells_median
   cells_std_deviation            = c_cells_std_deviation
   cells_min                      = c_cells_min
   cells_max                      = c_cells_max
   bound_part_faces_average       = c_bound_part_faces_average
   bound_part_faces_median        = c_bound_part_faces_median
   bound_part_faces_std_deviation = c_bound_part_faces_std_deviation
   bound_part_faces_min           = c_bound_part_faces_min
   bound_part_faces_max           = c_bound_part_faces_max
   bound_part_faces_sum           = c_bound_part_faces_sum

 end subroutine pdm_part_stat_get_



 !================================================================================
 !
 ! \brief Return an initialized coarse mesh object
 !
 ! \param [out]  cmId              Coarse mesh identifier
 !
 ! \param [in]   pt_comm           Communicator
 ! \param [in]   method            Choice between (1 for ParMETIS or 2 for PT-Scotch)
 ! \param [in]   nPart             Number of partitions
 ! \param [in]   nTPart            Total number of partitions
 ! \param [in]   have_cellTag      Presence d'un tableau de tags pour les cellules
 ! \param [in]   have_faceTag      Presence d'un tableau de tags pour les faces
 ! \param [in]   have_vtxTag       Presence d'un tableau de tags pour les sommets
 ! \param [in]   have_cellWeight   Presence d'un tableau de poids pour les cellules
 ! \param [in]   have_faceWeight   Presence d'un tableau de poids pour les faces
 ! \param [in]   have_faceGroup    Presence des tableaux de groupes de faces
 !================================================================================


subroutine pdm_part_coarse_mesh_create_ (cmId, &
                                        comm, &
                                        method, &
                                        renum_cell_method, &
                                        renum_face_method, &
                                        nPropertyCell, &
                                        renum_properties_cell, &
                                        nPropertyFace, &
                                        renum_properties_face, &
                                        nPart, &
                                        nTPart, &
                                        nFaceGroup,&
                                        have_cellTag,&
                                        have_faceTag,&
                                        have_vtxTag,&
                                        have_cellWeight, &
                                        have_faceWeight, &
                                        have_faceGroup)

    use pdm

    implicit none

    integer                     ::  cmId
    integer                     ::  comm
    character (len=*)           ::  method
    integer                     ::  nPart
    integer                     ::  nTPart
    integer                     ::  nFaceGroup
    integer                     ::  have_cellTag
    integer                     ::  have_faceTag
    integer                     ::  have_vtxTag
    integer                     ::  have_cellWeight
    integer                     ::  have_faceWeight
    integer                     ::  have_faceGroup
    character (len=*)           ::  renum_cell_method
    character (len=*)           ::  renum_face_method
    integer                     ::  nPropertyCell
    integer                     ::  renum_properties_cell
    integer                     ::  nPropertyFace
    integer                     ::  renum_properties_face

    integer                     ::  l_method
    integer                     ::  l_renum_cell_method
    integer                     ::  l_renum_face_method

    l_method = len(method)
    l_renum_cell_method = len(renum_cell_method)
    l_renum_face_method = len(renum_face_method)

    ! call pdm_part_coarse_mesh_create_cf (cmId, &
    !                                      comm, &
    !                                      method, &
    !                                      l_method, &
    !                                      renum_cell_method, &
    !                                      l_renum_cell_method, &
    !                                      renum_face_method, &
    !                                      l_renum_face_method, &
    !                                      nPropertyCell, &
    !                                      renum_properties_cell, &
    !                                      nPropertyFace, &
    !                                      renum_properties_face, &
    !                                      nPart, &
    !                                      nTPart, &
    !                                      nFaceGroup,&
    !                                      have_cellTag,&
    !                                      have_faceTag,&
    !                                      have_vtxTag,&
    !                                      have_cellWeight, &
    !                                      have_faceWeight, &
    !                                      have_faceGroup)

 end subroutine pdm_part_coarse_mesh_create_


 !================================================================================
 !
 ! \brief Get index of a renumbering face method
 !
 ! \param [in]       name   Name of the method
 ! \param [in, out]  idx    Index of method -1 otherwise
 !
 !================================================================================

  subroutine pdm_part_renum_method_face_idx_get_ (name, &
                                                  idx)

     use pdm

    implicit none

    character (len=*) :: name
    integer           :: idx

    integer           :: l_name

    l_name = len(name)

    ! call pdm_part_renum_method_face_idx_get_cf (name, l_name, idx)

  end subroutine pdm_part_renum_method_face_idx_get_


 !================================================================================
 !
 ! \brief Get index of a renumbering cell method
 !
 ! \param [in]       name   Name of the method
 ! \param [in, out]  idx    Index of method -1 otherwise
 !
 !================================================================================

 subroutine pdm_part_renum_method_cell_idx_get_ (name, &
                                                 idx)

   use pdm

   implicit none

   character (len=*) :: name
   integer           :: idx

   integer           :: l_name

   l_name = len(name)

   ! call pdm_part_renum_method_cell_idx_get_cf (name, l_name, idx)

 end subroutine pdm_part_renum_method_cell_idx_get_

 !================================================================================
 !
 ! \brief Get name of the face renumbering method
 !
 ! \param [in]  idx     Index of the method
 ! \param [in, out]     Name  of the method, '' otherwize
 !
 !================================================================================


subroutine pdm_part_renum_method_face_name_get_ (idx, &
                                                name)

   use pdm
   implicit none

   character (len = *) :: name
   integer  :: idx

   integer  :: l_name
   l_name = len(name)

   ! call pdm_part_renum_method_face_name_get_cf (name, l_name, idx)


end subroutine pdm_part_renum_method_face_name_get_


 !================================================================================
 !
 ! \brief Get name of the face renumbering method
 !
 ! \param [in]  idx     Index of the method
 ! \param [in, out]     Name  of the method, '' otherwize
 !
 !================================================================================


subroutine pdm_part_renum_method_cell_name_get_ (idx, &
                                                name)

   use pdm
   implicit none

   character (len = *) :: name
   integer  :: idx

   integer  :: l_name
   l_name = len(name)

   ! call pdm_part_renum_method_cell_name_get_cf (name, l_name, idx)


end subroutine pdm_part_renum_method_cell_name_get_

end module pdm_part
