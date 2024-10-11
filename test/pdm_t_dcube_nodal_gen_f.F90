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
  use pdm_dcube_nodal_gen
  use pdm_dmesh_nodal
  use pdm_multipart
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  integer,              parameter       :: comm = MPI_COMM_WORLD
  integer,              parameter       :: n_vtx_seg = 10
  double precision,     parameter       :: length = 5.
  double precision,     parameter       :: zero_x = 1.
  double precision,     parameter       :: zero_y = 1.
  double precision,     parameter       :: zero_z = 1.

  type(c_ptr)                           :: dcube
  type(c_ptr)                           :: dmn
  type(c_ptr)                           :: mpart

  integer                               :: n_face_group = 0
  integer                               :: dn_cell
  integer                               :: dn_face
  integer                               :: dn_vtx
  integer                               :: n_domains
  integer (kind = PDM_l_num_s), pointer :: n_part_domains(:)

  integer (kind = pdm_g_num_s), pointer :: dcell_vtx(:)       => null()
  integer (kind = pdm_g_num_s), pointer :: dface_vtx(:)       => null()
  double precision,             pointer :: dvtx_coord(:,:)    => null()
  integer (kind = pdm_l_num_s), pointer :: dface_group_idx(:) => null()
  integer (kind = pdm_g_num_s), pointer :: dface_group(:)     => null()

  integer (kind = pdm_g_num_s)          :: dims(4)
  integer (kind = pdm_l_num_s), pointer :: id_sections(:)

  integer (kind = pdm_l_num_s), pointer :: renum_cell_properties(:) => null()
  double precision,             pointer :: part_fraction(:)         => null()

  integer (c_int)                       :: n_vtx
  integer (c_int)                       :: n_face
  integer (c_int)                       :: n_cell

  double precision, pointer             :: coords(:,:)      => null()
  integer (kind = pdm_g_num_s), pointer :: vtx_ln_to_gn(:)  => null()
  integer (kind = pdm_g_num_s), pointer :: face_ln_to_gn(:) => null()
  integer (kind = pdm_l_num_s), pointer :: face_vtx(:)      => null()
  integer (kind = pdm_l_num_s), pointer :: face_vtx_idx(:)  => null()
  integer (kind = pdm_g_num_s), pointer :: cell_ln_to_gn(:) => null()
  integer (kind = pdm_l_num_s), pointer :: cell_face(:)     => null()
  integer (kind = pdm_l_num_s), pointer :: cell_face_idx(:) => null()

  integer                               :: code
  integer                               :: i_rank
  integer                               :: n_rank
  integer                               :: order
  integer                               :: i_domain
  integer                               :: i_part
  integer                               :: n_section
  integer                               :: id_section
  !-----------------------------------------------------------

  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)

  order = 1
  call pdm_dcube_nodal_gen_create (dcube,                &
                                   comm,                 &
                                   n_vtx_seg,            &
                                   n_vtx_seg,            &
                                   n_vtx_seg,            &
                                   length,               &
                                   zero_x,               &
                                   zero_y,               &
                                   zero_z,               &
                                   PDM_MESH_NODAL_HEXA8, &
                                   order,                &
                                   PDM_OWNERSHIP_KEEP)

  call pdm_dcube_nodal_gen_build (dcube, &
                                  dmn)

  call pdm_dcube_nodal_gen_dmesh_nodal_get (dcube, &
                                            dmn)

  call pdm_dmesh_nodal_section_g_dims_get (dmn,     &
                                           dims(1), &
                                           dims(2), &
                                           dims(3), &
                                           dims(4))

  dn_vtx = pdm_dmesh_nodal_n_vtx_get (dmn)

  call pdm_dmesh_nodal_vtx_get (dmn,        &
                                dvtx_coord, &
                                PDM_OWNERSHIP_KEEP)

  n_section = pdm_dmesh_nodal_n_section_get (dmn, &
                                             PDM_GEOMETRY_KIND_VOLUMIC)

  call pdm_dmesh_nodal_sections_id_get (dmn,                       &
                                        PDM_GEOMETRY_KIND_VOLUMIC, &
                                        id_sections)

  dn_cell = pdm_dmesh_nodal_section_n_elt_get (dmn,                       &
                                               PDM_GEOMETRY_KIND_VOLUMIC, &
                                               id_sections(1))

  call pdm_dmesh_nodal_section_std_get (dmn,                       &
                                        PDM_GEOMETRY_KIND_VOLUMIC, &
                                        id_sections(1),            &
                                        dcell_vtx,                 &
                                        PDM_OWNERSHIP_KEEP)

  n_section = pdm_dmesh_nodal_n_section_get (dmn, &
                                             PDM_GEOMETRY_KIND_SURFACIC)

  call pdm_dmesh_nodal_sections_id_get (dmn,                        &
                                        PDM_GEOMETRY_KIND_SURFACIC, &
                                        id_sections)

  dn_face = pdm_dmesh_nodal_section_n_elt_get (dmn,                        &
                                               PDM_GEOMETRY_KIND_SURFACIC, &
                                               id_sections(1))

  call pdm_dmesh_nodal_section_std_get (dmn,                        &
                                        PDM_GEOMETRY_KIND_SURFACIC, &
                                        id_sections(1),             &
                                        dface_vtx,                  &
                                        PDM_OWNERSHIP_KEEP)

  call pdm_dmesh_nodal_section_group_elmt_get (dmn,                        &
                                               PDM_GEOMETRY_KIND_SURFACIC, &
                                               n_face_group,               &
                                               dface_group_idx,            &
                                               dface_group,                &
                                               PDM_OWNERSHIP_KEEP)

  call pdm_dmesh_nodal_create (dmn,     &
                               comm,    &
                               3,       &
                               dims(4), &
                               dims(1), &
                               dims(2), &
                               dims(3))

  call pdm_dmesh_nodal_coord_set (dmn,        &
                                  dn_vtx,     &
                                  dvtx_coord, &
                                  PDM_OWNERSHIP_USER)

  call pdm_dmesh_nodal_section_add (dmn,                       &
                                    PDM_GEOMETRY_KIND_VOLUMIC, &
                                    PDM_MESH_NODAL_HEXA8,      &
                                    id_section)

  call pdm_dmesh_nodal_section_std_set (dmn,                       &
                                        PDM_GEOMETRY_KIND_VOLUMIC, &
                                        id_section,                &
                                        dn_cell,                   &
                                        dcell_vtx,                 &
                                        PDM_OWNERSHIP_USER)

  call pdm_dmesh_nodal_section_add (dmn,                        &
                                    PDM_GEOMETRY_KIND_SURFACIC, &
                                    PDM_MESH_NODAL_QUAD4,       &
                                    id_section)

  call pdm_dmesh_nodal_section_std_set (dmn,                        &
                                        PDM_GEOMETRY_KIND_SURFACIC, &
                                        id_section,                 &
                                        dn_face,                    &
                                        dface_vtx,                  &
                                        PDM_OWNERSHIP_USER)

  call pdm_dmesh_nodal_section_group_elmt_set (dmn,                        &
                                               PDM_GEOMETRY_KIND_SURFACIC, &
                                               n_face_group,               &
                                               dface_group_idx,            &
                                               dface_group,                &
                                               PDM_OWNERSHIP_USER)

  call pdm_dmesh_nodal_generate_distribution (dmn)

  allocate( n_part_domains(1) ) ; n_part_domains(1) = 1
  n_domains = 1
  i_domain  = 0
  i_part    = 0

  call pdm_multipart_create(mpart,                       &
                            n_domains,                   &
                            n_part_domains,              &
                            PDM_FALSE,                   &
                            PDM_SPLIT_DUAL_WITH_HILBERT, &
                            PDM_PART_SIZE_HOMOGENEOUS,   &
                            part_fraction,               &
                            comm,                        &
                            PDM_OWNERSHIP_KEEP)

  call pdm_multipart_dmesh_nodal_set (mpart,    &
                                      i_domain, &
                                      dmn)

  call pdm_multipart_set_reordering_options (mpart,                      &
                                             i_domain,                   &
                                             "PDM_PART_RENUM_CELL_NONE", &
                                             renum_cell_properties,      &
                                             "PDM_PART_RENUM_FACE_NONE")

  call pdm_multipart_compute (mpart)

  call pdm_multipart_part_ln_to_gn_get (mpart,               &
                                        i_domain,            &
                                        i_part,              &
                                        PDM_MESH_ENTITY_VTX, &
                                        vtx_ln_to_gn,        &
                                        PDM_OWNERSHIP_KEEP,  &
                                        n_vtx)

  call pdm_multipart_part_vtx_coord_get (mpart,              &
                                         i_domain,           &
                                         i_part,             &
                                         coords,             &
                                         PDM_OWNERSHIP_KEEP, &
                                         n_vtx)

  call pdm_multipart_part_ln_to_gn_get (mpart,                &
                                        i_domain,             &
                                        i_part,               &
                                        PDM_MESH_ENTITY_FACE, &
                                        face_ln_to_gn,        &
                                        PDM_OWNERSHIP_KEEP,   &
                                        n_face)

  call pdm_multipart_part_connectivity_get (mpart,                          &
                                            i_domain,                       &
                                            i_part,                         &
                                            PDM_CONNECTIVITY_TYPE_FACE_VTX, &
                                            face_vtx_idx,                   &
                                            face_vtx,                       &
                                            PDM_OWNERSHIP_KEEP,             &
                                            n_face)

  call pdm_multipart_part_ln_to_gn_get (mpart,                &
                                        i_domain,             &
                                        i_part,               &
                                        PDM_MESH_ENTITY_CELL, &
                                        cell_ln_to_gn,        &
                                        PDM_OWNERSHIP_KEEP,   &
                                        n_cell)

  call pdm_multipart_part_connectivity_get (mpart,                           &
                                            i_domain,                        &
                                            i_part,                          &
                                            PDM_CONNECTIVITY_TYPE_CELL_FACE, &
                                            cell_face_idx,                   &
                                            cell_face,                       &
                                            PDM_OWNERSHIP_KEEP,              &
                                            n_cell)

  call pdm_multipart_free (mpart)
  call pdm_dmesh_nodal_free (dmn)
  call pdm_dcube_nodal_gen_free (dcube)

  deallocate( n_part_domains )

  call mpi_finalize(code)

end program testf
