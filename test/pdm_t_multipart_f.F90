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
  use pdm_multipart
  use pdm_dcube_gen
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------
  integer :: code
  !-----------------


  ! Initialize MPI
  call mpi_init(code)

  ! Run tests
  call run_test(use_dpart_id=.false., &
                n_part=1)

  call run_test(use_dpart_id=.true., &
                n_part=1)

  ! Finalize MPI
  call mpi_finalize(code)


  contains


  subroutine run_test(use_dpart_id, &
                      n_part)

    implicit none


    !-----------------------------------------------------------
    logical, intent(in) :: use_dpart_id
    integer, intent(in) :: n_part
    !-----------------------------------------------------------
    ! MPI
    integer,                  parameter   :: comm = MPI_COMM_WORLD
    integer                               :: code
    integer                               :: i_rank
    integer                               :: n_rank
    ! MULTIPART
    type(c_ptr)                           :: multipart
    integer(c_int)                        :: split_method
    integer(c_int),           parameter   :: n_domain = 1
    ! integer(c_int),           parameter   :: i_domain = -1
    integer(kind=PDM_l_num_s), pointer    :: n_part_domains(:)
    double precision,          pointer    :: part_fraction(:)
    integer(kind=PDM_l_num_s), pointer    :: renum_cell_properties(:)
    ! MESH
    type(c_ptr)                           :: dcube
    integer(pdm_g_num_s), parameter       :: n_vtx_seg = 4
    double precision,     parameter       :: length = 5.
    double precision,     parameter       :: zero_x = 1.
    double precision,     parameter       :: zero_y = 1.
    double precision,     parameter       :: zero_z = 1.
    integer                               :: n_face_group
    integer                               :: dn_cell
    integer                               :: dn_face
    integer                               :: dn_vtx
    integer                               :: sface_vtx
    integer                               :: sface_group
    integer (kind = pdm_g_num_s), pointer :: dface_cell(:)
    integer (kind = pdm_l_num_s), pointer :: dcell_face_idx(:)
    integer (kind = pdm_g_num_s), pointer :: dcell_face(:)
    integer (kind = pdm_l_num_s), pointer :: dface_vtx_idx(:)
    integer (kind = pdm_g_num_s), pointer :: dface_vtx(:)
    double precision,             pointer :: dvtx_coord(:,:)
    integer (kind = pdm_l_num_s), pointer :: dface_group_idx(:)
    integer (kind = pdm_g_num_s), pointer :: dface_group(:)

    integer(pdm_l_num_s),         pointer :: dpart_id(:)

    integer(c_int)                        :: n_cell
    integer(pdm_g_num_s),         pointer :: cell_ln_to_gn(:)
    !-----------------------------------------------------------

    ! Initializations
    multipart = C_NULL_PTR
    n_part_domains        => null()
    part_fraction         => null()
    renum_cell_properties => null()

    dcube = C_NULL_PTR

    n_face_group = -1
    dn_cell      = -1
    dn_face      = -1
    dn_vtx       = -1
    sface_vtx    = -1
    sface_group  = -1
    dface_cell      => null()
    dcell_face_idx  => null()
    dcell_face      => null()
    dface_vtx_idx   => null()
    dface_vtx       => null()
    dvtx_coord      => null()
    dface_group_idx => null()
    dface_group     => null()

    dpart_id        => null()
    n_cell = 0
    cell_ln_to_gn   => null()


    call mpi_comm_rank(comm, i_rank, code)
    call mpi_comm_size(comm, n_rank, code)

    if (i_rank .eq. 0) then
      write(*, *) "-- Run test --"
      write(*, *) "  use_dpart_id :", use_dpart_id
      write(*, *) "  n_part       :", n_part
    end if


    ! Initialize multipart
    split_method = PDM_SPLIT_DUAL_WITH_HILBERT

    allocate(n_part_domains(n_domain))

    n_part_domains(1:n_domain) = n_part

    if (i_rank .eq. 0) then
      write(*, *) "PDM_multipart_create"
    end if

    call PDM_multipart_create(multipart, &
                              n_domain, &
                              n_part_domains, &
                              PDM_FALSE, &
                              split_method, &
                              PDM_PART_SIZE_HOMOGENEOUS, &
                              part_fraction, &
                              comm, &
                              PDM_OWNERSHIP_KEEP)

    ! Reordering options
    ! if (i_rank .eq. 0) then
    !   write(*, *) "PDM_multipart_set_reordering_options"
    ! end if

    ! call PDM_multipart_set_reordering_options(multipart, &
    !                                           i_domain, &
    !                                           "PDM_PART_RENUM_CELL_CUTHILL", &
    !                                           renum_cell_properties, &
    !                                           "PDM_PART_RENUM_FACE_LEXICOGRAPHIC")

    ! Generate Mesh (case : n_domain = 1)
    ! > dcube
    if (i_rank .eq. 0) then
      write(*, *) "> dcube"
    end if

    call pdm_dcube_gen_init(dcube,              &
                            comm,               &
                            n_vtx_seg,          &
                            length,             &
                            zero_x,             &
                            zero_y,             &
                            zero_z,             &
                            PDM_OWNERSHIP_KEEP)

    call pdm_dcube_gen_dim_get(dcube,           &
                              n_face_group,    &
                              dn_cell,         &
                              dn_face,         &
                              dn_vtx,          &
                              sface_vtx,       &
                              sface_group)

    call pdm_dcube_gen_data_get(dcube,           &
                                dface_cell,      &
                                dface_vtx_idx,   &
                                dface_vtx,       &
                                dvtx_coord,      &
                                dface_group_idx, &
                                dface_group)


    call PDM_multipart_block_set(multipart, &
                                 0, &
                                 dn_cell, &
                                 dn_face, &
                                 dn_vtx, &
                                 n_face_group, &
                                 dcell_face_idx, &
                                 dcell_face, &
                                 dface_cell, &
                                 dface_vtx_idx, &
                                 dface_vtx, &
                                 dvtx_coord, &
                                 dface_group_idx, &
                                 dface_group)


    ! Set dpart_id
    if (use_dpart_id) then
      allocate(dpart_id(dn_cell))
      dpart_id(1:dn_cell) = i_rank ! similar to PDM_SPLIT_DUAL_WITH_IMPLICIT

      call PDM_multipart_dpart_id_set(multipart, &
                                      0,         &
                                      dpart_id)
    endif


    ! Run
    call PDM_multipart_compute(multipart)

    ! Get
    call PDM_multipart_part_ln_to_gn_get(multipart,            &
                                         0,                    &
                                         0,                    &
                                         PDM_MESH_ENTITY_CELL, &
                                         cell_ln_to_gn,        &
                                         PDM_OWNERSHIP_KEEP,   &
                                         n_cell)

    ! Free
    deallocate(n_part_domains)
    call pdm_dcube_gen_free(dcube)

    if (i_rank .eq. 0) then
      write(*, *) "PDM_multipart_free"
    end if
    if (use_dpart_id) then
      deallocate(dpart_id)
    endif
    call PDM_multipart_free(multipart)

    if (i_rank .eq. 0) then
      write(*, *) "-- End"
      write(*, *) ""
    end if

  end subroutine run_test

end program testf
