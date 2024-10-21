#ifndef PDM_ISOSURFACE_TEST_UTILS_H
#define PDM_ISOSURFACE_TEST_UTILS_H

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2023       ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_isosurface.h"
#include "pdm_isosurface_priv.h"

#include "pdm_multipart.h"
#include "pdm_multipart_priv.h"
/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/**
 *
 * \brief Read arguments from the command line
 *
 */

void
PDM_isosurface_test_utils_read_args
(
  int                    argc,
  char                 **argv,
  int                   *n_part,
  char                 **mesh_name,
  char                 **sol_name,
  int                   *visu,
  int                   *n_isovalues,
  double               **isovalues,
  PDM_Mesh_nodal_elt_t  *elt_type,
  int                   *randomize,
  PDM_g_num_t           *n_vtx_seg,
  int                   *use_part_mesh,
  int                   *generate_edges,
  int                   *local
);


/**
 *
 * \brief Mesh generation for ngon cases
 *
 */

void
PDM_isosurface_test_utils_gen_mesh
(
  PDM_MPI_Comm          comm,
  const char           *filename,
  int                   n_part,
  PDM_g_num_t           n_vtx_seg,
  int                   randomize,
  PDM_Mesh_nodal_elt_t  elt_type,
  int                   generate_edges,
  PDM_multipart_t     **mpart,
  PDM_part_mesh_t      *pmesh,
  PDM_dmesh_t         **out_dmesh
);


/**
 *
 * \brief Mesh generation for nodal cases
 *
 */

void
PDM_isosurface_test_utils_gen_mesh_nodal
(
  PDM_MPI_Comm            comm,
  const char             *filename,
  int                     n_part,
  PDM_g_num_t             n_vtx_seg,
  int                     randomize,
  PDM_Mesh_nodal_elt_t    elt_type,
  PDM_part_mesh_nodal_t **out_pmn,
  PDM_dmesh_nodal_t     **out_dmn
);


/**
 *
 * \brief Compute iso field from vertex coordinates
 *
 */

void
PDM_isosurface_test_utils_compute_iso_field
(
  int     n_vtx,
  double *vtx_coord,
  double *vtx_field
);


/**
 *
 * \brief Compute interpolate field from vertex coordinates
 *
 */

void
PDM_isosurface_test_utils_compute_itp_field
(
  int     n_vtx,
  double *vtx_coord,
  double *vtx_field
);


/**
 *
 * \brief Interpolation source field onto dist iso mesh
 *
 */

void
PDM_isosurface_test_utils_dist_interpolation
(
  PDM_isosurface_t *isos,
  int               i_iso,
  double           *itp_dfield_vtx,
  double           *itp_dfield_face,
  double           *itp_dfield_cell,
  double          **iso_itp_dfield_vtx,
  double          **iso_itp_dfield_edge,
  double          **iso_itp_dfield_face
);


/**
 *
 * \brief Interpolation source field onto part iso mesh
 *
 */

void
PDM_isosurface_test_utils_part_interpolation
(
  PDM_isosurface_t *isos,
  int               i_iso,
  int               n_part,
  int               local,
  double          **itp_field_vtx,
  double          **itp_field_face,
  double          **itp_field_cell,
  double         ***iso_itp_field_vtx,
  double         ***iso_itp_field_edge,
  double         ***iso_itp_field_face
);



/**
 *
 * \brief Write vtk output from isosurface dist result
 *
 */

void
PDM_isosurface_test_utils_dist_vtk
(
  PDM_isosurface_t *isos,
  int               i_iso,
  double           *iso_vtx_fld,
  double           *iso_edge_fld,
  double           *iso_face_fld,
  PDM_MPI_Comm      comm
);


/**
 *
 * \brief Write vtk output from isosurface part result
 *
 */

void
PDM_isosurface_test_utils_part_vtk
(
  PDM_isosurface_t  *isos,
  int                id_iso,
  int                n_part,
  double           **iso_vtx_fld,
  double           **iso_edge_fld,
  double           **iso_face_fld,
  PDM_MPI_Comm       comm
);


#ifdef  __cplusplus
}
#endif

#endif // PDM_ISOSURFACE_TEST_UTILS_H
