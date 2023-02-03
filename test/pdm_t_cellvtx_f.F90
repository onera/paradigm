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


program testf

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use PDM_mesh_nodal
  use iso_c_binding

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  integer,                   parameter :: f_comm  = MPI_COMM_WORLD
  integer,                   parameter :: elements(5) = (/ 3,4,3,0,4 /)

  type(c_ptr)                          :: mesh = C_NULL_PTR
  double precision,          pointer   :: vtx_coord(:,:)   => null()
  integer(kind=pdm_g_num_s), pointer   :: vtx_ln_to_gn(:)  => null()
  integer(kind=pdm_g_num_s), pointer   :: cell_ln_to_gn(:) => null()
  integer,                   pointer   :: cell_vtx_idx(:)  => null()
  integer,                   pointer   :: cell_vtx_nb(:)   => null()
  integer,                   pointer   :: cell_vtx(:)      => null()
  integer,                   pointer   :: cell_face_idx(:) => null()
  integer,                   pointer   :: cell_face_nb(:)  => null()
  integer,                   pointer   :: cell_face(:)     => null()
  integer(kind=pdm_g_num_s), pointer   :: face_ln_to_gn(:) => null()
  integer,                   pointer   :: face_vtx_idx(:)  => null()
  integer,                   pointer   :: face_vtx_nb(:)   => null()
  integer,                   pointer   :: face_vtx(:)      => null()
  integer,                   pointer   :: nums(:,:)        => null()

  integer                              :: code
  integer                              :: i_rank
  integer                              :: n_rank
  integer                              :: n_vtx
  integer                              :: n_face
  integer                              :: n_cell

  integer                              :: i, ifac, jfac, isom, jsom, tmp(6)
  integer                              :: numf, adrs, cfac(2), csom(2)
  integer                              :: fid = 13
  character                            :: strnum
  !-----------------------------------------------------------

  ! Init
  call mpi_init(code)
  call mpi_comm_rank(f_comm, i_rank, code)
  call mpi_comm_size(f_comm, n_rank, code)

  if (n_rank .ne. 2) then
    print *,'Error : 2 MPI processes are mandatory'
    call mpi_finalize(code)
    stop
  end if

  ! Read mesh
  if (i_rank .eq. 0) then
    write(*, *) "-- Read mesh"
  end if

  write (strnum, '(i1)') i_rank+1
  fid = fid + i_rank
  open(unit=fid, file="meshes/mixed_elements_cellvtx."//strnum, action='read')

  do i = 1,26
    read(fid,*)
  end do

  read(fid,*) n_vtx
  read(fid,*) n_face
  read(fid,*) n_cell
  read(fid,*)

  read(fid,*)
  allocate(vtx_ln_to_gn(n_vtx))
  allocate(vtx_coord(3,n_vtx))
  do i = 1,n_vtx
    read(fid,*) vtx_ln_to_gn(i), vtx_coord(1:3,i)
  end do

  read(fid,*)
  allocate(face_ln_to_gn(n_face))
  allocate(face_vtx_idx(n_face+1))
  allocate(face_vtx_nb(n_face))
  allocate(face_vtx(4*n_face))
  face_vtx_idx(1) = 1
  do i = 1,n_face
    read(fid,*) face_ln_to_gn(i), face_vtx_nb(i), tmp(1:face_vtx_nb(i))
    face_vtx_idx(i+1) = face_vtx_idx(i) + face_vtx_nb(i)
    face_vtx(face_vtx_idx(i):face_vtx_idx(i+1)-1) = tmp(1:face_vtx_nb(i))
    do isom = 1,face_vtx_nb(i)
      face_vtx(face_vtx_idx(i)+isom-1) = minloc(abs(vtx_ln_to_gn-face_vtx(face_vtx_idx(i)+isom-1)),1)
    end do
  end do

  read(fid,*)
  allocate(cell_face_idx(n_cell+1))
  allocate(cell_face_nb(n_cell))
  allocate(cell_face(6*n_cell))
  allocate(cell_ln_to_gn(n_cell))
  cell_face_idx(1) = 1
  do i = 1,n_cell
    read(fid,*) cell_ln_to_gn(i), cell_face_nb(i), tmp(1:cell_face_nb(i))
    cell_face_idx(i+1) = cell_face_idx(i) + cell_face_nb(i)
    cell_face(cell_face_idx(i):cell_face_idx(i+1)-1) = tmp(1:cell_face_nb(i))
    do ifac = 1,cell_face_nb(i)
      cell_face(cell_face_idx(i)+ifac-1) = minloc(abs(face_ln_to_gn-cell_face(cell_face_idx(i)+ifac-1)),1)
    end do
  end do
  close(fid)

  ! Convert connectivity
  if (i_rank .eq. 0) then
    write(*, *) "-- Convert connectivity to cell -> vtx"
  end if

  allocate(cell_vtx_idx(n_cell+1))
  allocate(cell_vtx_nb(n_cell))
  allocate(nums(8,n_cell))
  cell_vtx_nb = 0
  nums = 0

  do i = 1,n_cell
    do ifac = 1,cell_face_nb(i)
      numf = cell_face(cell_face_idx(i)+ifac-1)
      adrs = face_vtx_idx(numf)
      do isom = 1,face_vtx_nb(numf)
          if (any(nums(1:cell_vtx_nb(i),i) == face_vtx(adrs+isom-1))) cycle
          cell_vtx_nb(i) = cell_vtx_nb(i) + 1
          nums(cell_vtx_nb(i),i) = face_vtx(adrs+isom-1)
      end do
    end do
  end do

  deallocate(nums)
  allocate(cell_vtx(sum(cell_vtx_nb)))
  cell_vtx_idx(1) = 1
  cell_vtx = -1

  do i = 1,n_cell
    do ifac = 1,cell_face_nb(i)
      numf = cell_face(cell_face_idx(i)+ifac-1)
      if (face_vtx_nb(numf) == elements(cell_vtx_nb(i)-3)) exit
    end do

    adrs = face_vtx_idx(numf)
    cell_vtx(cell_vtx_idx(i):cell_vtx_idx(i)+face_vtx_nb(numf)-1) = face_vtx(adrs:adrs+face_vtx_nb(numf)-1)

    do isom = 1,cell_vtx_nb(i)-face_vtx_nb(numf)
      do ifac = 1,cell_face_nb(i)
        cfac(1) = cell_face(cell_face_idx(i)+ifac-1)
        if (cfac(1) == numf) cycle
        csom(1) = face_vtx_idx(cfac(1))
        if (all(face_vtx(csom(1):csom(1)+face_vtx_nb(cfac(1))-1) /= cell_vtx(cell_vtx_idx(i)+isom-1))) cycle
        do jfac = ifac+1,cell_face_nb(i)
          cfac(2) = cell_face(cell_face_idx(i)+jfac-1)
          if (cfac(2) == numf) cycle
          csom(2) = face_vtx_idx(cfac(2))
          if (all(face_vtx(csom(2):csom(2)+face_vtx_nb(cfac(2))-1) /= cell_vtx(cell_vtx_idx(i)+isom-1))) cycle
          do jsom = 1,face_vtx_nb(cfac(2))
            if (face_vtx(csom(2)+jsom-1) == cell_vtx(cell_vtx_idx(i)+isom-1)) cycle
            if (all(face_vtx(csom(1):csom(1)+face_vtx_nb(cfac(1))-1) /= face_vtx(csom(2)+jsom-1))) cycle
            cell_vtx(cell_vtx_idx(i)+face_vtx_nb(numf)+isom-1) = face_vtx(csom(2)+jsom-1)
            exit
          end do
          if (cell_vtx(cell_vtx_idx(i)+face_vtx_nb(numf)+isom-1) /= -1) exit
        end do
        if (cell_vtx(cell_vtx_idx(i)+face_vtx_nb(numf)+isom-1) /= -1) exit
      end do
    end do

    cell_vtx_idx(i+1) = cell_vtx_idx(i) + cell_vtx_nb(i)

  end do

  ! Convert Fortran -> C
  cell_vtx_idx = cell_vtx_idx - 1
  face_vtx_idx = face_vtx_idx - 1
  cell_face_idx = cell_face_idx - 1

  if (i_rank .eq. 0) then
    write(*, *) "-- Create nodal mesh"
  end if

  ! Création de l'objet "MAILLAGE NODAL"
  call pdm_mesh_nodal_create                    & !
  (mesh,                                        & !- IDENTIFICATEUR OBJET MAILLAGE NODAL
   1,                                           & !- NOMBRE DE PARTITIONS DU MAILLAGE NODAL SUR LE PROCESSUS COURANT
   f_comm)                                        !- COMMUNICATEUR MPI

  if (i_rank .eq. 0) then
    write(*, *) "-- Set vertices"
  end if

  ! Définition des sommets du maillage nodal
  call pdm_mesh_nodal_coord_set                 &
  (mesh,                                        & !- IDENTIFICATEUR OBJET MAILLAGE NODAL
   0,                                           & !- INDICE DE PARTITION DU MAILLAGE NODAL
   n_vtx,                                       & !- NOMBRE DE SOMMETS
   vtx_coord,                                   & !- COORDONNEES DES SOMMETS
   vtx_ln_to_gn,                                & !- NUMEROTATION ABSOLUE DES SOMMETS
   PDM_OWNERSHIP_USER)                            !- OWNERSHIP

  if (i_rank .eq. 0) then
    write(*, *) "-- Set connectivity"
  end if

  ! Définition de la connectivité du maillage nodal
  call pdm_mesh_nodal_cells_cellvtx_add         &
  (mesh,                                        & !- IDENTIFICATEUR OBJET MAILLAGE NODAL
   0,                                           & !- INDICE DE PARTITION DU MAILLAGE NODAL
   n_cell,                                      & !- NOMBRE DE CELLULES
   cell_vtx_idx,                                & !- ADRESSES DES NUMEROS DE SOMMETS PAR CELLULE
   cell_vtx_nb,                                 & !- NOMBRES DE SOMMETS PAR CELLULE
   cell_vtx,                                    & !- NUMEROS DE SOMMETS PAR CELLULE
   cell_ln_to_gn,                               & !- NUMEROTATION ABSOLUE DES CELLULES
   PDM_OWNERSHIP_KEEP)                            !- OWNERSHIP

  ! Free memory
  deallocate(vtx_coord)
  deallocate(vtx_ln_to_gn)
  deallocate(face_ln_to_gn)
  deallocate(face_vtx_idx)
  deallocate(face_vtx_nb)
  deallocate(face_vtx)
  deallocate(cell_face_idx)
  deallocate(cell_face_nb)
  deallocate(cell_face)
  deallocate(cell_ln_to_gn)
  deallocate(cell_vtx_idx)
  deallocate(cell_vtx_nb)
  deallocate(cell_vtx)
  call pdm_mesh_nodal_free(mesh)

  if (i_rank .eq. 0) then
    write(*, *) "-- End"
  end if

  call mpi_finalize(code)

end program testf
