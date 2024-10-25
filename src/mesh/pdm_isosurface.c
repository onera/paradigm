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
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"

#include "pdm_error.h"
#include "pdm_logging.h"

#include "pdm_array.h"
#include "pdm_binary_search.h"

#include "pdm_mesh_nodal.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_extract_part.h"
#include "pdm_vtk.h"
#include "pdm_part_mesh_nodal_to_part_mesh.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_triangulate.h"

#include "pdm_distrib.h"

#include "pdm_part_to_block.h"
#include "pdm_multipart.h"

#include "pdm_isosurface.h"
#include "pdm_isosurface_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

static const char *_isosurface_timer_step_name[ISO_TIMER_N_STEPS] = {
  "dist to part ", // ISO_DIST_TO_PART
  "initial field", // ISO_COMPUTE_INIT_FIELD
  "extract mesh ", // ISO_EXTRACT
  "ngonize      ", // ISO_NGONIZE
  "extract field", // ISO_COMPUTE_EXTRACT_FIELD
  "contouring   ", // ISO_CONTOURING
  "part to dist ", // ISO_PART_TO_DIST
  "exch protocol", // ISO_BUILD_EXCH_PROTOCOL
  "total        "  // ISO_TIMER_TOTAL
};


void
isosurface_timer_start
(
  PDM_isosurface_t         *isos,
  _isosurface_timer_step_t  step
)
{
  PDM_MPI_Barrier(isos->comm); // keep?
  isos->t_start[step] = PDM_MPI_Wtime();
}


void
isosurface_timer_end
(
  PDM_isosurface_t         *isos,
  _isosurface_timer_step_t  step
)
{
  isos->t_end[step] = PDM_MPI_Wtime();

  double delta_t = isos->t_end[step] - isos->t_start[step];

  isos->times_current[step]  = delta_t;
  isos->times_cumul  [step] += delta_t;
}


static void
_dump_delta_time
(
       PDM_MPI_Comm  comm,
       double        delta,
 const char         *prefix
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  double delta_min;
  PDM_MPI_Allreduce(&delta, &delta_min, 1, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);

  double delta_max;
  PDM_MPI_Allreduce(&delta, &delta_max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  if (i_rank == 0) {
    printf("PDM_isosurface timer : %s (min/max in seconds) : %8.2e / %8.2e\n",
           prefix,
           delta_min, delta_max);
  }
}




static inline int
_iso_surface_kind_n_coeff
(
  PDM_iso_surface_kind_t kind
)
{
  switch (kind) {
    case PDM_ISO_SURFACE_KIND_PLANE:
      return 4;
    case PDM_ISO_SURFACE_KIND_SPHERE:
      return 4;
    case PDM_ISO_SURFACE_KIND_ELLIPSE:
      return 6;
    case PDM_ISO_SURFACE_KIND_QUADRIC:
      return 10;
    default:
      return 0;
  }
}



static
inline
double
_plane_field
(
  const double x,
  const double y,
  const double z,
  double *plane_equation
)
{
  return  plane_equation[0] * x
        + plane_equation[1] * y
        + plane_equation[2] * z
        - plane_equation[3];
}



static
inline
double
_sphere_field
(
  const double x,
  const double y,
  const double z,
  double *sphere_equation
)
{
  return   pow(x-sphere_equation[0], 2.)
         + pow(y-sphere_equation[1], 2.)
         + pow(z-sphere_equation[2], 2.)
         - pow(  sphere_equation[3], 2.) ;
}



static
inline
double
_ellipse_field
(
  const double x,
  const double y,
  const double z,
  double *ellipse_equation
)
{
  return   pow( (x-ellipse_equation[0]) / ellipse_equation[3], 2.)
         + pow( (y-ellipse_equation[1]) / ellipse_equation[4], 2.)
         + pow( (z-ellipse_equation[2]) / ellipse_equation[5], 2.)
         -         ellipse_equation[6];
}



static
inline
double
_quadric_field
(
  const double x,
  const double y,
  const double z,
  double *quadric_equation
)
{
  return   quadric_equation[6] * pow( (x-quadric_equation[0]) / quadric_equation[3], 2.)
         + quadric_equation[7] * pow( (y-quadric_equation[1]) / quadric_equation[4], 2.)
         + quadric_equation[8] * pow( (z-quadric_equation[2]) / quadric_equation[5], 2.)
         - quadric_equation[9];
}



static
inline
double
_heart_field
(
  const double x,
  const double y,
  const double z,
  double *heart_equation
)
{
  PDM_UNUSED(heart_equation);
  double a=1.;
  double b=2.;
  // return pow(x*x + y*y -1,3.) - x*x*y*y*y; // 2D
  return pow( pow(       x,2.)
            + pow((1.+b)*y,2.)
            + pow(       z,2.)
            -            1.   ,3.)
        -   x*x*z*z*z
        - a*y*y*z*z*z;
        // +           1.;
}

int
isosurface_is_nodal
(
 PDM_isosurface_t *isos
)
{
  return PDM_ABS(isos->entry_mesh_type) == 3;
}


/**
 * \brief Check if edges were provided (once for all isosurfaces)
 *
 */
static void
_do_we_have_edges
(
 PDM_isosurface_t *isos
)
{
  if (isos->we_have_edges >= 0) {
    // We have already determined if we had edges
    return;
  }

  if (isosurface_is_nodal(isos)) {
    // Irrelevant for nodal mesh
    return;
  }

  int i_have_edges    = 0;
  int i_have_face_vtx = 0;

  if (isos->entry_is_part == 0) {
    // Block-distributed
    i_have_edges    = (isos->distrib_edge  != NULL) && (isos->dface_edge_idx != NULL) && (isos->dface_edge != NULL) && (isos->dedge_vtx != NULL);
    i_have_face_vtx = (isos->dface_vtx_idx != NULL) && (isos->dface_vtx != NULL);
  }
  else {
    // Partitioned
    for (int i_part = 0; i_part < isos->n_part; i_part++) {

      if (isos->n_face[i_part] == 0) {
        continue;
      }

      if (isos->face_edge_idx[i_part] != NULL &&
          isos->face_edge    [i_part] != NULL &&
          isos->edge_vtx     [i_part] != NULL) {
        i_have_edges = 1;
      }

      if (isos->face_vtx_idx[i_part] != NULL &&
          isos->face_vtx    [i_part] != NULL) {
        i_have_face_vtx = 1;
      }
    }
  }

  PDM_MPI_Allreduce(&i_have_edges, &isos->we_have_edges, 1, PDM_MPI_INT, PDM_MPI_MAX, isos->comm);

  if (!isos->we_have_edges) {
    int we_have_face_vtx;
    PDM_MPI_Allreduce(&i_have_face_vtx, &we_have_face_vtx, 1, PDM_MPI_INT, PDM_MPI_MAX, isos->comm);
    if (!we_have_face_vtx) {
      PDM_error(__FILE__, __LINE__, 0, "Either face->vtx or {face->edge, edge->vtx} connectivities must be provided\n");
    }
  }

  // log_trace("isos->we_have_edges = %d\n", isos->we_have_edges);
}


/**
 * \brief Perform implicit partitioning (once for all isosurfaces)
 *
 */
static void
_dist_to_part
(
  PDM_isosurface_t *isos
 )
{
  if (isos->dist_to_part_computed) {
    // We have already switched to paritioned view
    return;
  }

  // We only need to build the block-to-part instance if there is at least
  // one isosurface from discrete field
  PDM_bool_t build_btp_vtx = PDM_FALSE;
  for (int i_iso = 0; i_iso < isos->n_isosurface; i_iso++) {
    if (isos->isosurfaces[i_iso].kind == PDM_ISO_SURFACE_KIND_FIELD) {
      build_btp_vtx = PDM_TRUE;
      break;
    }
  }


  isos->dist_to_part_computed = 1;

  if (isos->entry_is_part != 0) {
    PDM_error(__FILE__, __LINE__, 0, "Expected block-distributed but got partitioned\n");
  }

  int i_rank;
  PDM_MPI_Comm_rank(isos->comm, &i_rank);

  isos->n_part     = 1;
  isos->iso_n_part = 1;

  if (isosurface_is_nodal(isos)) {
    // Nodal
    assert(isos->pmesh_nodal == NULL);

    // TODO: optimize?
    PDM_multipart_t *mpart = PDM_multipart_create(1,
                                                  &isos->n_part,
                                                  PDM_FALSE,
                                                  PDM_SPLIT_DUAL_WITH_IMPLICIT,
                                                  PDM_PART_SIZE_HOMOGENEOUS,
                                                  NULL,
                                                  isos->comm,
                                                  PDM_OWNERSHIP_KEEP);

    PDM_multipart_dmesh_nodal_set(mpart, 0, isos->dmesh_nodal);

    PDM_multipart_compute(mpart);

    PDM_multipart_get_part_mesh_nodal(mpart,
                                      0,
                                      &isos->pmesh_nodal,
                                      PDM_OWNERSHIP_USER);

    PDM_multipart_free(mpart);

    // (Re)create block_to_part for transferring vtx data from block to part (TODO: extract from multipart?)
    if (build_btp_vtx) {
      int          n_vtx        = PDM_part_mesh_nodal_n_vtx_get    (isos->pmesh_nodal, 0);
      PDM_g_num_t *vtx_ln_to_gn = PDM_part_mesh_nodal_vtx_g_num_get(isos->pmesh_nodal, 0);

      const PDM_g_num_t *distrib_vtx      = PDM_DMesh_nodal_distrib_vtx_get(isos->dmesh_nodal);
      const PDM_g_num_t *pvtx_ln_to_gn[1] = {vtx_ln_to_gn};

      PDM_block_to_part_t *btp_vtx = PDM_block_to_part_create(distrib_vtx,
                                                              pvtx_ln_to_gn,
                                                              &n_vtx,
                                                              1,
                                                              isos->comm);

      isos->btp_vtx = btp_vtx;
    }
  }

  else {
    // Ngon
    int n_cell = 0;
    int n_face = 0;
    int n_edge = 0;
    int n_vtx  = 0;

    PDM_g_num_t *cell_ln_to_gn = NULL;
    PDM_g_num_t *face_ln_to_gn = NULL;
    PDM_g_num_t *edge_ln_to_gn = NULL;
    PDM_g_num_t *vtx_ln_to_gn  = NULL;

    int *cell_face_idx = NULL;
    int *cell_face     = NULL;
    int *face_edge_idx = NULL;
    int *face_edge     = NULL;
    int *face_vtx_idx  = NULL;
    int *face_vtx      = NULL;
    int *edge_vtx      = NULL;

    if (isos->entry_mesh_dim == 3) {
      n_cell = isos->distrib_cell[i_rank+1] - isos->distrib_cell[i_rank];

      PDM_malloc(cell_ln_to_gn, n_cell, PDM_g_num_t);
      for (int i = 0; i < n_cell; i++) {
        cell_ln_to_gn[i] = isos->distrib_cell[i_rank] + i + 1;
      }

      PDM_part_dconnectivity_to_pconnectivity_sort_single_part(isos->comm,
                                                               isos->distrib_cell,
                                                               isos->dcell_face_idx,
                                                               isos->dcell_face,
                                                               n_cell,
                                                               cell_ln_to_gn,
                                                               &n_face,
                                                               &face_ln_to_gn,
                                                               &cell_face_idx,
                                                               &cell_face);
    }
    else if (isos->entry_mesh_dim == 2) {
      n_face = isos->distrib_face[i_rank+1] - isos->distrib_face[i_rank];

      PDM_malloc(face_ln_to_gn, n_face, PDM_g_num_t);
      for (int i = 0; i < n_face; i++) {
        face_ln_to_gn[i] = isos->distrib_face[i_rank] + i + 1;
      }
    }

    if (isos->we_have_edges) {
      PDM_part_dconnectivity_to_pconnectivity_sort_single_part(isos->comm,
                                                               isos->distrib_face,
                                                               isos->dface_edge_idx,
                                                               isos->dface_edge,
                                                               n_face,
                                                               face_ln_to_gn,
                                                               &n_edge,
                                                               &edge_ln_to_gn,
                                                               &face_edge_idx,
                                                               &face_edge);

      int dn_edge = isos->distrib_edge[i_rank+1] - isos->distrib_edge[i_rank];
      int *dedge_vtx_idx = PDM_array_new_idx_from_const_stride_int(2, dn_edge);
      int *edge_vtx_idx  = NULL;
      PDM_part_dconnectivity_to_pconnectivity_sort_single_part(isos->comm,
                                                               isos->distrib_edge,
                                                               dedge_vtx_idx,
                                                               isos->dedge_vtx,
                                                               n_edge,
                                                               edge_ln_to_gn,
                                                               &n_vtx,
                                                               &vtx_ln_to_gn,
                                                               &edge_vtx_idx,
                                                               &edge_vtx);
      PDM_free(dedge_vtx_idx);
      PDM_free(edge_vtx_idx );
    }
    else {
      PDM_part_dconnectivity_to_pconnectivity_sort_single_part(isos->comm,
                                                               isos->distrib_face,
                                                               isos->dface_vtx_idx,
                                                               isos->dface_vtx,
                                                               n_face,
                                                               face_ln_to_gn,
                                                               &n_vtx,
                                                               &vtx_ln_to_gn,
                                                               &face_vtx_idx,
                                                               &face_vtx);
    }

    const PDM_g_num_t *pvtx_ln_to_gn[1] = {vtx_ln_to_gn};
    PDM_block_to_part_t *btp_vtx = PDM_block_to_part_create(isos->distrib_vtx,
                                                            pvtx_ln_to_gn,
                                                            &n_vtx,
                                                            1,
                                                            isos->comm);

    int one = 1;
    PDM_block_to_part_exch(btp_vtx,
                           3*sizeof(double),
                           PDM_STRIDE_CST_INTERLACED,
                           &one,
                (void   *) isos->dvtx_coord,
                           NULL,
                (void ***) &isos->vtx_coord);

    // Surfaces
    if (isos->entry_mesh_dim==3) {
      const PDM_g_num_t *pface_ln_to_gn[1] = {face_ln_to_gn};
      PDM_part_distgroup_to_partgroup(isos->comm,
                                      isos->distrib_face,
                                      isos->n_group_face,
                                      isos->dgroup_face_idx,
                                      isos->dgroup_face,
                                      isos->n_part,
                                      &n_face,
                                      pface_ln_to_gn,
                                      &isos->group_face_idx,
                                      &isos->group_face,
                                      &isos->group_face_gnum);
    }

    if (build_btp_vtx) {
      /* Store in struct */
      isos->btp_vtx = btp_vtx; // useful to keep for transferring discrete fields
    }
    else {
      PDM_block_to_part_free(btp_vtx);
    }

    PDM_malloc(isos->n_cell       , 1, int          );
    PDM_malloc(isos->n_face       , 1, int          );
    PDM_malloc(isos->n_edge       , 1, int          );
    PDM_malloc(isos->n_vtx        , 1, int          );
    PDM_malloc(isos->cell_gnum    , 1, PDM_g_num_t *);
    PDM_malloc(isos->face_gnum    , 1, PDM_g_num_t *);
    PDM_malloc(isos->edge_gnum    , 1, PDM_g_num_t *);
    PDM_malloc(isos->vtx_gnum     , 1, PDM_g_num_t *);
    PDM_malloc(isos->cell_face_idx, 1, int         *);
    PDM_malloc(isos->cell_face    , 1, int         *);
    PDM_malloc(isos->face_edge_idx, 1, int         *);
    PDM_malloc(isos->face_edge    , 1, int         *);
    PDM_malloc(isos->face_vtx_idx , 1, int         *);
    PDM_malloc(isos->face_vtx     , 1, int         *);
    PDM_malloc(isos->edge_vtx     , 1, int         *);

    isos->n_cell       [0] = n_cell;
    isos->n_face       [0] = n_face;
    isos->n_edge       [0] = n_edge;
    isos->n_vtx        [0] = n_vtx;
    isos->cell_gnum    [0] = cell_ln_to_gn;
    isos->face_gnum    [0] = face_ln_to_gn;
    isos->edge_gnum    [0] = edge_ln_to_gn;
    isos->vtx_gnum     [0] = vtx_ln_to_gn;
    isos->cell_face_idx[0] = cell_face_idx;
    isos->cell_face    [0] = cell_face;
    isos->face_edge_idx[0] = face_edge_idx;
    isos->face_edge    [0] = face_edge;
    isos->face_vtx_idx [0] = face_vtx_idx;
    isos->face_vtx     [0] = face_vtx;
    isos->edge_vtx     [0] = edge_vtx;
  }
}


static void
_compute_iso_field
(
  PDM_isosurface_t *isos,
  int               id_isosurface,
  int               use_extract
)
{
  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];

  int n_part = isos->n_part;
  if (use_extract) {
    n_part = isos->iso_n_part;
  }

  if (_iso->kind == PDM_ISO_SURFACE_KIND_FIELD) {
    if (isos->entry_is_part == 1) {
      if (_iso->field==NULL) {
        PDM_error(__FILE__, __LINE__, 0, "Field seems not to be defined for iso with id %d\n", id_isosurface);
      }
    }
    else {
      assert(isos->entry_is_part == 0);
      if (_iso->dfield==NULL) {
        PDM_error(__FILE__, __LINE__, 0, "Field seems not to be defined for iso with id %d\n", id_isosurface);
      }
    }

    if (use_extract) {

      if (isos->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) {
        PDM_malloc(isos->extract_field, n_part, double *);
        for (int i_part = 0; i_part < n_part; i_part++) {
          int n_extract_vtx = PDM_extract_part_n_entity_get(isos->extrp,
                                                            i_part,
                                                            PDM_MESH_ENTITY_VTX);

          PDM_malloc(isos->extract_field[i_part], n_extract_vtx, double);

          for (int i = 0; i < n_extract_vtx; i++) {
            int i_parent = isos->extract_vtx_lnum[i_part][i] - 1;
            isos->extract_field[i_part][i] = _iso->field[i_part][i_parent];
          }
        }
      } // End if LOCAL

      else if (isos->extract_kind == PDM_EXTRACT_PART_KIND_REEQUILIBRATE) {
        PDM_part_to_part_t *ptp = NULL;
        PDM_extract_part_part_to_part_get(isos->extrp,
                                          PDM_MESH_ENTITY_VTX,
                                         &ptp,
                                          PDM_OWNERSHIP_KEEP);
        assert(ptp != NULL);

        int request = -1;
        PDM_part_to_part_reverse_iexch(ptp,
                                       PDM_MPI_COMM_KIND_P2P ,
                                       PDM_STRIDE_CST_INTERLACED,
                                       PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                       1,
                                       sizeof(double),
                                       NULL,
                      (const void  **) _iso->field,
                                       NULL,
                            (void ***) &isos->extract_field,
                                       &request);

        PDM_part_to_part_reverse_iexch_wait(ptp, request);
      } // End if REEQUILIBRATE

      else {
        PDM_error(__FILE__, __LINE__, 0, "Invalid extract_kind %d\n", isos->extract_kind);
      }

    } // End if use_extract

    else {
      // Initial mesh (before extraction)
      if (isos->entry_is_part == 0) {
        // Transfer discrete field from block to part
        assert(isos->btp_vtx != NULL);
        assert(_iso->dfield  != NULL);

        int one = 1;
        PDM_block_to_part_exch(isos->btp_vtx,
                               sizeof(double),
                               PDM_STRIDE_CST_INTERLACED,
                               &one,
                    (void   *) _iso->dfield,
                               NULL,
                    (void ***) &_iso->field);
      }
    }

    return;
  }

  if (!use_extract) {
    if (_iso->field==NULL) {
      PDM_malloc(_iso->field, n_part, double *);
      for (int i = 0; i < n_part; i++) {
        _iso->field[i] = NULL;
      }
    }
  }

  if (isos->entry_is_part == 0) {
    // Block-distributed
    assert(isos->dist_to_part_computed);
  }

  double **field     = NULL;
  int     *n_vtx     = NULL;
  double **vtx_coord = NULL;
  PDM_malloc(n_vtx    , n_part, int     );
  PDM_malloc(vtx_coord, n_part, double *);

  if (use_extract) {
    // We work on extracted mesh
    if (isos->extract_field == NULL) {
      PDM_malloc(isos->extract_field, n_part, double *);
    }
    field = isos->extract_field;
    for (int i_part = 0; i_part < n_part; i_part++) {
      if (isosurface_is_nodal(isos)) {
        n_vtx    [i_part] = PDM_part_mesh_nodal_n_vtx_get    (isos->extract_pmesh_nodal, i_part);
        vtx_coord[i_part] = PDM_part_mesh_nodal_vtx_coord_get(isos->extract_pmesh_nodal, i_part);
      }
      else {
        n_vtx[i_part] = PDM_extract_part_vtx_coord_get(isos->extrp,
                                                       i_part,
                                                       &vtx_coord[i_part],
                                                       PDM_OWNERSHIP_KEEP);
      }
    }
  } // End if use_extract

  else {
    // We work on initial mesh (before extraction)
    field = _iso->field;
    if (isosurface_is_nodal(isos)) {
      for (int i_part = 0; i_part < n_part; i_part++) {
        n_vtx    [i_part] = PDM_part_mesh_nodal_n_vtx_get    (isos->pmesh_nodal, i_part);
        vtx_coord[i_part] = PDM_part_mesh_nodal_vtx_coord_get(isos->pmesh_nodal, i_part);
      }
    }
    else {
      for (int i_part = 0; i_part < n_part; i_part++) {
        n_vtx    [i_part] = isos->n_vtx    [i_part];
        vtx_coord[i_part] = isos->vtx_coord[i_part];
      }
    }
  }

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_malloc(field[i_part], n_vtx[i_part], double);
  }

  /* Fill */
  if (_iso->kind == PDM_ISO_SURFACE_KIND_FUNCTION) {
    // User-provided field function
    if (_iso->field_function_python != NULL) {
      // Python
      for (int i_part = 0; i_part < n_part; i_part++) {
        for (int i_vtx = 0; i_vtx < n_vtx[i_part]; i_vtx++) {
          field[i_part][i_vtx] = _iso->field_function_python(isos->python_object,
                                                             id_isosurface,
                                                             vtx_coord[i_part][3*i_vtx  ],
                                                             vtx_coord[i_part][3*i_vtx+1],
                                                             vtx_coord[i_part][3*i_vtx+2]);
        }
      }
    }
    else {
      // C
      assert(_iso->field_function != NULL);
      for (int i_part = 0; i_part < n_part; i_part++) {
        for (int i_vtx = 0; i_vtx < n_vtx[i_part]; i_vtx++) {
          _iso->field_function(vtx_coord[i_part][3*i_vtx  ],
                               vtx_coord[i_part][3*i_vtx+1],
                               vtx_coord[i_part][3*i_vtx+2],
                               &field   [i_part][  i_vtx  ]);
        }
      }
    }
  }

  else {
    // Built-in field function
    double (*field_function) (const double, const double, const double, double *) = NULL;

    switch (_iso->kind) {
      case PDM_ISO_SURFACE_KIND_PLANE: {
        field_function = &_plane_field;
        break;
      }
      case PDM_ISO_SURFACE_KIND_SPHERE: {
        field_function = &_sphere_field;
        break;
      }
      case PDM_ISO_SURFACE_KIND_ELLIPSE: {
        field_function = &_ellipse_field;
        break;
      }
      case PDM_ISO_SURFACE_KIND_QUADRIC: {
        field_function = &_quadric_field;
        break;
      }
      case PDM_ISO_SURFACE_KIND_HEART: {
        field_function = &_heart_field;
        break;
      }
      default: {
        PDM_error(__FILE__, __LINE__, 0, "Invalid isosurface type %d for id_isosurface %d.\n", _iso->kind, id_isosurface);
      }
    }

    for (int i_part = 0; i_part < n_part; i_part++) {
      for (int i_vtx = 0; i_vtx < n_vtx[i_part]; i_vtx++) {
        field[i_part][i_vtx] = field_function(vtx_coord[i_part][3*i_vtx  ],
                                              vtx_coord[i_part][3*i_vtx+1],
                                              vtx_coord[i_part][3*i_vtx+2],
                                              _iso->eq_coeffs);
      }
    }
  }

  PDM_free(n_vtx    );
  PDM_free(vtx_coord);
}


extern inline int
_isosurface_is_at_0_level(
  const double v,
  const double tol
)
{
  return (PDM_ABS(v) <= tol);
}


extern inline int
_isosurface_cross_0_level
(
  const double v0,
  const double v1,
  const double tol
)
{
  return (PDM_ABS(v0) > tol) && (PDM_ABS(v1) > tol) && (v0*v1 < 0);
}


extern inline int
_isosurface_cross_any_level
(
  const double v0,
  const double v1,
  const int    n_isovalues,
  const double isovalues[],
  const double tol
)
{
  int n_crossings = 0;
  for (int i = 0; i < n_isovalues; i++) {
    n_crossings += _isosurface_cross_0_level(v0 - isovalues[i], v1 - isovalues[i], tol);
  }

  return n_crossings;
}


static inline int
_sign
(
  const double v,
  const double tol
)
{
  return (v > tol);
}


extern inline int
_isosurface_cross_0_level_ngon
(
  const double v0,
  const double v1,
  const double tol
)
{
  return _sign(v0, tol) != _sign(v1, tol);
}


extern inline int
_isosurface_cross_any_level_ngon
(
  const double v0,
  const double v1,
  const int    n_isovalues,
  const double isovalues[],
  const double tol
)
{
  int n_crossings = 0;
  for (int i = 0; i < n_isovalues; i++) {
    n_crossings += _isosurface_cross_0_level_ngon(v0 - isovalues[i], v1 - isovalues[i], tol);
  }

  return n_crossings;
}




/**
 * \brief Convert group info into tag
 */
static int
_convert_group_info_to_tag
(
  PDM_part_mesh_nodal_t  *pmn,
  int                     i_part,
  int                     n_elt,
  PDM_geometry_kind_t     geom_kind,
  int                   **elt_bnd_tag_out
)
{
  int dbg = 0;

  *elt_bnd_tag_out = PDM_array_zeros_int(n_elt);
  int *elt_bnd_tag = *elt_bnd_tag_out;
  int n_group = PDM_part_mesh_nodal_n_group_get(pmn, geom_kind);
  for (int i_group=0; i_group<n_group; ++i_group) {
    int          n_group_elmt    = 0;
    int         *group_elmt_lnum = 0;
    PDM_g_num_t *group_elmt_gnum = 0;
    PDM_part_mesh_nodal_group_get(pmn,
                                  geom_kind,
                                  i_part,
                                  i_group,
                                 &n_group_elmt,
                                 &group_elmt_lnum,
                                 &group_elmt_gnum,
                                  PDM_OWNERSHIP_KEEP);
    for (int i_elmt=0; i_elmt<n_group_elmt; ++i_elmt) {
      assert(elt_bnd_tag[group_elmt_lnum[i_elmt]-1]==0);
      elt_bnd_tag[group_elmt_lnum[i_elmt]-1] = i_group+1;
    }
  }
  if (dbg==1) {
    PDM_log_trace_array_int(elt_bnd_tag, n_elt, "elt_bnd_tag ::");
  }
  return n_group;
}


/**
 * \brief Extract nodal elements of interest
 */
static void
_extract_nodal
(
 PDM_isosurface_t *isos,
 int               id_isosurface
)
{
  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];

  isos->extrp = PDM_extract_part_create(isos->entry_mesh_dim,
                                        isos->n_part,
                                        isos->iso_n_part,
                                        isos->extract_kind,
                                        isos->part_method,
                                        PDM_FALSE,
                                        PDM_OWNERSHIP_KEEP,
                                        isos->comm);

  PDM_extract_part_part_nodal_set(isos->extrp, isos->pmesh_nodal);

  // Get principal pmne
  PDM_geometry_kind_t geom_kind = PDM_part_mesh_nodal_principal_geom_kind_get(isos->pmesh_nodal);

  PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(isos->pmesh_nodal, geom_kind);

  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  // Select elements of interest, ie. crossing at least one isovalue
  int  *n_extract    = PDM_array_zeros_int(isos->n_part);
  int **extract_lnum = NULL;
  PDM_malloc(extract_lnum, isos->n_part, int *);

  for (int i_part = 0; i_part < isos->n_part; i_part++) {

    int n_elt_tot = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne, i_part);
    PDM_malloc(extract_lnum[i_part], n_elt_tot, int);

    int i_parent = -1;
    for (int i_section = 0; i_section < n_section; i_section++) {

      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                              sections_id[i_section],
                                                              i_part);

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, sections_id[i_section]);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                 sections_id[i_section],
                                                                 i_part,
                                                                 PDM_OWNERSHIP_BAD_VALUE);

      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        /* Polygons */
        int *connec_idx = NULL;
        int *connec     = NULL;
        PDM_part_mesh_nodal_elmts_section_poly2d_get(pmne,
                                                     sections_id[i_section],
                                                     i_part,
                                                     &connec_idx,
                                                     &connec,
                                                     PDM_OWNERSHIP_BAD_VALUE);

        for (int i_elt = 0; i_elt < n_elt; i_elt++) {

          int is_selected = 0;
          int n_vtx_elt = connec_idx[i_elt+1] - connec_idx[i_elt];
          int *_connec = connec + connec_idx[i_elt];
          for (int i = 0; i < n_vtx_elt; i++) {
            int i_vtx0 = _connec[ i             ] - 1;
            int i_vtx1 = _connec[(i+1)%n_vtx_elt] - 1;

            double val0 = _iso->field[i_part][i_vtx0];
            double val1 = _iso->field[i_part][i_vtx1];
            if (_isosurface_cross_any_level(val0, val1, _iso->n_isovalues, _iso->isovalues, isos->ISOSURFACE_EPS)) {
              is_selected = 1;
              break;
            }
            else {
              for (int i_isovalue = 0; i_isovalue < _iso->n_isovalues; i_isovalue++) {
                if (_isosurface_is_at_0_level(val0-_iso->isovalues[i_isovalue], isos->ISOSURFACE_EPS) ||
                    _isosurface_is_at_0_level(val1-_iso->isovalues[i_isovalue], isos->ISOSURFACE_EPS)) {
                  is_selected = 1;
                  break;
                }
              }
              if (is_selected) {
                break;
              }
            }
          }

          if (is_selected) {
            if (parent_num != NULL) {
              i_parent = parent_num[i_elt];
            }
            else {
              i_parent++;
            }

            extract_lnum[i_part][n_extract[i_part]++] = i_parent + 1;
          }

        } // End loop on elements

      }

      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        /* Polyhedra */
        int          _n_face              = 0;
        PDM_g_num_t *_face_ln_to_gn       = NULL;
        int         *face_vtx_idx         = NULL;
        int         *face_vtx             = NULL;
        PDM_g_num_t *_cell_ln_to_gn       = NULL;
        int         *cell_face_idx        = NULL;
        int         *cell_face            = NULL;
        int         *_parent_num          = NULL;
        PDM_g_num_t *_parent_entity_g_num = NULL;
        PDM_part_mesh_nodal_elmts_section_poly3d_get(pmne,
                                                     sections_id[i_section],
                                                     i_part,
                                                     &_n_face,
                                                     &_face_ln_to_gn,
                                                     &face_vtx_idx,
                                                     &face_vtx,
                                                     &_cell_ln_to_gn,
                                                     &cell_face_idx,
                                                     &cell_face,
                                                     &_parent_num,
                                                     &_parent_entity_g_num,
                                                     PDM_OWNERSHIP_BAD_VALUE);

        for (int i_cell = 0; i_cell < n_elt; i_cell++) {
          int is_selected = 0;
          for (int idx_face = cell_face_idx[i_cell]; idx_face < cell_face_idx[i_cell+1]; idx_face++) {
            int i_face = PDM_ABS(cell_face[idx_face]) - 1;
            int face_vtx_n = face_vtx_idx[i_face+1] - face_vtx_idx[i_face];
            int *fv = face_vtx + face_vtx_idx[i_face];
            for (int i = 0; i < face_vtx_n; i++) {
              int i_vtx0 = fv[ i              ] - 1;
              int i_vtx1 = fv[(i+1)%face_vtx_n] - 1;
              double val0 = _iso->field[i_part][i_vtx0];
              double val1 = _iso->field[i_part][i_vtx1];
              if (_isosurface_cross_any_level(val0, val1, _iso->n_isovalues, _iso->isovalues, isos->ISOSURFACE_EPS)) {
                is_selected = 1;
                break;
              }
              else {
                for (int i_isovalue = 0; i_isovalue < _iso->n_isovalues; i_isovalue++) {
                  if (_isosurface_is_at_0_level(val0-_iso->isovalues[i_isovalue], isos->ISOSURFACE_EPS) ||
                      _isosurface_is_at_0_level(val1-_iso->isovalues[i_isovalue], isos->ISOSURFACE_EPS)) {
                    is_selected = 1;
                    break;
                  }
                }
                if (is_selected) {
                  break;
                }
              }
            }

            if (is_selected) continue;

          } // End of loop on faces of current cell

          if (is_selected) {
            extract_lnum[i_part][n_extract[i_part]++] = i_cell + 1;
          }
        } // End of loop on cells

      }

      else {
        /* Standard elements */
        int         *connec              = NULL;
        PDM_g_num_t *g_num               = NULL;
        int         *_parent_num         = NULL;
        PDM_g_num_t *parent_entity_g_num = NULL;
        int          order               = -1;
        const char  *ho_ordering         = NULL;
        PDM_part_mesh_nodal_elmts_section_std_ho_get(pmne,
                                                     sections_id[i_section],
                                                     i_part,
                                                     &connec,
                                                     &g_num,
                                                     &_parent_num,
                                                     &parent_entity_g_num,
                                                     &order,
                                                     &ho_ordering,
                                                     PDM_OWNERSHIP_BAD_VALUE);

        assert(order == 1);

        const int *elt_edge_vtx = NULL;
        int n_elt_edge = PDM_edge_vtx_per_elmt(t_elt,
                                               &elt_edge_vtx);
        if (elt_edge_vtx == NULL) {
          // Skip irrelevant sections
          continue;
        }

        int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);

        for (int i_elt = 0; i_elt < n_elt; i_elt++) {

          int is_selected = 0;
          int *_connec = connec + n_vtx_elt*i_elt;

          for (int i = 0; i < n_elt_edge; i++) {
            int i_vtx0 = _connec[elt_edge_vtx[2*i  ]] - 1;
            int i_vtx1 = _connec[elt_edge_vtx[2*i+1]] - 1;

            double val0 = _iso->field[i_part][i_vtx0];
            double val1 = _iso->field[i_part][i_vtx1];
            if (_isosurface_cross_any_level(val0, val1, _iso->n_isovalues, _iso->isovalues, isos->ISOSURFACE_EPS)) {
              is_selected = 1;
              break;
            }
            else {
              for (int i_isovalue = 0; i_isovalue < _iso->n_isovalues; i_isovalue++) {
                if (_isosurface_is_at_0_level(val0-_iso->isovalues[i_isovalue], isos->ISOSURFACE_EPS) ||
                    _isosurface_is_at_0_level(val1-_iso->isovalues[i_isovalue], isos->ISOSURFACE_EPS)) {
                  is_selected = 1;
                  break;
                }
              }
              if (is_selected) {
                break;
              }
            }
          }

          if (is_selected) {
            if (parent_num != NULL) {
              i_parent = parent_num[i_elt];
            }
            else {
              i_parent++;
            }

            extract_lnum[i_part][n_extract[i_part]++] = i_parent + 1;
          }

        } // End loop on elements
      }

    } // End loop on sections

    PDM_realloc(extract_lnum[i_part], extract_lnum[i_part], n_extract[i_part], int);

    PDM_extract_part_selected_lnum_set(isos->extrp,
                                       i_part,
                                       n_extract   [i_part],
                                       extract_lnum[i_part]);

  } // End loop on parts


  PDM_extract_part_compute(isos->extrp);

  PDM_extract_part_part_mesh_nodal_get(isos->extrp, &isos->extract_pmesh_nodal, PDM_OWNERSHIP_KEEP);

  if (isos->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) {
    if (isos->entry_mesh_dim==3) {
      isos->extract_cell_lnum = extract_lnum;
      PDM_malloc(isos->extract_face_lnum, isos->iso_n_part, int *);
      for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {
        PDM_extract_part_parent_lnum_get(isos->extrp,
                                         i_part,
                                         PDM_MESH_ENTITY_FACE,
                                         &isos->extract_face_lnum[i_part],
                                         PDM_OWNERSHIP_USER);
      }
    }
    else {
      isos->extract_face_lnum = extract_lnum;
    }
  }
  else {
    for (int i_part = 0; i_part < isos->n_part; i_part++) {
      PDM_free(extract_lnum[i_part]);
    }
    PDM_free(extract_lnum);
  }
  PDM_free(n_extract);

  // for (int i_part = 0; i_part < isos->n_part; i_part++) {
  //   PDM_free(extract_lnum[i_part]);
  // }
  // PDM_free(n_extract   );
  // PDM_free(extract_lnum);
}


/**
 * \brief Extract mesh entities of interest
 */
static void
_extract_ngon
(
 PDM_isosurface_t *isos,
 int               id_isosurface
)
{
  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];

  isos->extrp = PDM_extract_part_create(isos->entry_mesh_dim,
                                        isos->n_part,
                                        isos->iso_n_part,
                                        isos->extract_kind,
                                        isos->part_method,
                                        PDM_FALSE,
                                        PDM_OWNERSHIP_KEEP,
                                        isos->comm);

  // Select entities of interest, ie. crossing at least one isovalue
  int  *n_extract    = PDM_array_zeros_int(isos->n_part);
  int **extract_lnum = NULL;
  PDM_malloc(extract_lnum, isos->n_part, int *);

  if (isos->entry_mesh_dim == 2) {
    // 2D
    for (int i_part = 0; i_part < isos->n_part; i_part++) {
      PDM_malloc(extract_lnum[i_part], isos->n_face[i_part], int);
    }

    if (isos->we_have_edges) {
      // face->vtx
      for (int i_part = 0; i_part < isos->n_part; i_part++) {
        for (int i_face = 0; i_face < isos->n_face[i_part]; i_face++) {
          int is_selected = 0;
          for (int idx_edge = isos->face_edge_idx[i_part][i_face]; idx_edge < isos->face_edge_idx[i_part][i_face+1]; idx_edge++) {
            int i_edge = PDM_ABS(isos->face_edge[i_part][idx_edge]) - 1;
            int i_vtx0 = isos->edge_vtx[i_part][2*i_edge  ] - 1;
            int i_vtx1 = isos->edge_vtx[i_part][2*i_edge+1] - 1;
            double val0 = _iso->field[i_part][i_vtx0];
            double val1 = _iso->field[i_part][i_vtx1];
            if (_isosurface_cross_any_level_ngon(val0, val1, _iso->n_isovalues, _iso->isovalues, isos->ISOSURFACE_EPS)) {
              is_selected = 1;
              break;
            }
          } // End of loop on edges of current face

          if (is_selected) {
            extract_lnum[i_part][n_extract[i_part]++] = i_face + 1;
          }
        } // End of loop on faces
      } // End of loop on parts
    }

    else {
      // face->edge->vtx
      for (int i_part = 0; i_part < isos->n_part; i_part++) {
        for (int i_face = 0; i_face < isos->n_face[i_part]; i_face++) {
          int is_selected = 0;
          int face_vtx_n = isos->face_vtx_idx[i_part][i_face+1] - isos->face_vtx_idx[i_part][i_face];
          int *fv = isos->face_vtx[i_part] + isos->face_vtx_idx[i_part][i_face];
          for (int i = 0; i < face_vtx_n; i++) {
            int i_vtx0 = fv[ i              ] - 1;
            int i_vtx1 = fv[(i+1)%face_vtx_n] - 1;
            double val0 = _iso->field[i_part][i_vtx0];
            double val1 = _iso->field[i_part][i_vtx1];
            if (_isosurface_cross_any_level_ngon(val0, val1, _iso->n_isovalues, _iso->isovalues, isos->ISOSURFACE_EPS)) {
              is_selected = 1;
              break;
            }
          }

          if (is_selected) {
            extract_lnum[i_part][n_extract[i_part]++] = i_face + 1;
          }
        } // End of loop on faces
      } // End of loop on parts
    }
  } // End 2D

  else {
    // 3D
    assert(isos->entry_mesh_dim == 3);
    for (int i_part = 0; i_part < isos->n_part; i_part++) {
      PDM_malloc(extract_lnum[i_part], isos->n_cell[i_part], int);
    }
    if (isos->we_have_edges == 0) {
      // cell->face->vtx
      for (int i_part = 0; i_part < isos->n_part; i_part++) {
        for (int i_cell = 0; i_cell < isos->n_cell[i_part]; i_cell++) {
          int is_selected = 0;
          for (int idx_face = isos->cell_face_idx[i_part][i_cell]; idx_face < isos->cell_face_idx[i_part][i_cell+1]; idx_face++) {
            int i_face = PDM_ABS(isos->cell_face[i_part][idx_face]) - 1;
            int face_vtx_n = isos->face_vtx_idx[i_part][i_face+1] - isos->face_vtx_idx[i_part][i_face];
            int *fv = isos->face_vtx[i_part] + isos->face_vtx_idx[i_part][i_face];
            for (int i = 0; i < face_vtx_n; i++) {
              int i_vtx0 = fv[ i              ] - 1;
              int i_vtx1 = fv[(i+1)%face_vtx_n] - 1;
              double val0 = _iso->field[i_part][i_vtx0];
              double val1 = _iso->field[i_part][i_vtx1];
              if (_isosurface_cross_any_level_ngon(val0, val1, _iso->n_isovalues, _iso->isovalues, isos->ISOSURFACE_EPS)) {
                is_selected = 1;
                break;
              }
            }

            if (is_selected) continue;

          } // End of loop on faces of current cell

          if (is_selected) {
            extract_lnum[i_part][n_extract[i_part]++] = i_cell + 1;
          }
        } // End of loop on cells
      } // End of loop on parts
    }

    else {
      // cell->face->edge->vtx
      for (int i_part = 0; i_part < isos->n_part; i_part++) {
        for (int i_cell = 0; i_cell < isos->n_cell[i_part]; i_cell++) {
          int is_selected = 0;
          for (int idx_face = isos->cell_face_idx[i_part][i_cell]; idx_face < isos->cell_face_idx[i_part][i_cell+1]; idx_face++) {
            int i_face = PDM_ABS(isos->cell_face[i_part][idx_face]) - 1;
            for (int idx_edge = isos->face_edge_idx[i_part][i_face]; idx_edge < isos->face_edge_idx[i_part][i_face+1]; idx_edge++) {
              int i_edge = PDM_ABS(isos->face_edge[i_part][idx_edge]) - 1;
              int i_vtx0 = isos->edge_vtx[i_part][2*i_edge  ] - 1;
              int i_vtx1 = isos->edge_vtx[i_part][2*i_edge+1] - 1;
              double val0 = _iso->field[i_part][i_vtx0];
              double val1 = _iso->field[i_part][i_vtx1];
              if (_isosurface_cross_any_level_ngon(val0, val1, _iso->n_isovalues, _iso->isovalues, isos->ISOSURFACE_EPS)) {
                is_selected = 1;
                break;
              }
            } // End of loop on edges of current face

            if (is_selected) continue;
          } // End of loop on faces of current cell

          if (is_selected) {
            extract_lnum[i_part][n_extract[i_part]++] = i_cell + 1;
          }
        } // End of loop on cells
      } // End of loop on parts
    }

  } // End 3D

  for (int i_part = 0; i_part < isos->n_part; i_part++) {
    PDM_realloc(extract_lnum[i_part], extract_lnum[i_part], n_extract[i_part], int);
  }

  if (isos->entry_mesh_dim==3) {
    PDM_extract_part_n_group_set(isos->extrp, PDM_BOUND_TYPE_FACE, isos->n_group_face);
  }

  for (int i_part = 0; i_part < isos->n_part; i_part++) {
    PDM_extract_part_part_set(isos->extrp,
                              i_part,
                              isos->n_cell       [i_part],
                              isos->n_face       [i_part],
                              isos->n_edge       [i_part],
                              isos->n_vtx        [i_part],
                              isos->cell_face_idx[i_part],
                              isos->cell_face    [i_part],
                              isos->face_edge_idx[i_part],
                              isos->face_edge    [i_part],
                              isos->edge_vtx     [i_part],
                              isos->face_vtx_idx [i_part],
                              isos->face_vtx     [i_part],
                              isos->cell_gnum    [i_part],
                              isos->face_gnum    [i_part],
                              isos->edge_gnum    [i_part],
                              isos->vtx_gnum     [i_part],
                              isos->vtx_coord    [i_part]);

    PDM_extract_part_selected_lnum_set(isos->extrp,
                                       i_part,
                                       n_extract   [i_part],
                                       extract_lnum[i_part]);

    if (isos->entry_mesh_dim==3) {
      for (int i_group = 0; i_group < isos->n_group_face; i_group++) {
        PDM_extract_part_part_group_set(isos->extrp,
                                        i_part,
                                        i_group,
                                        PDM_BOUND_TYPE_FACE,
                                        isos->group_face_idx[i_part][i_group+1] - isos->group_face_idx[i_part][i_group],
                                        isos->group_face     [i_part] + isos->group_face_idx[i_part][i_group],
                                        isos->group_face_gnum[i_part] + isos->group_face_idx[i_part][i_group]);
      }
    }
  }

  PDM_extract_part_compute(isos->extrp);

  // Retrieve extracted mesh and store it as a PDM_part_mesh_t instance
  PDM_extract_part_part_mesh_get(isos->extrp,
                                 &isos->extract_pmesh,
                                 PDM_FALSE);

  if (isos->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) {
    if (isos->entry_mesh_dim==3) {
      isos->extract_cell_lnum = extract_lnum;
      for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {
        PDM_extract_part_parent_lnum_get(isos->extrp,
                                         i_part,
                                         PDM_MESH_ENTITY_FACE,
                                         &isos->extract_face_lnum[i_part],
                                         PDM_OWNERSHIP_USER);
      }
    }
    else {
      isos->extract_face_lnum = extract_lnum;
    }
  }
  else {
    for (int i_part = 0; i_part < isos->n_part; i_part++) {
      PDM_free(extract_lnum[i_part]);
    }
    PDM_free(extract_lnum);
  }
  PDM_free(n_extract);
}


/**
 * \brief Extract mesh regions that contain the isosurface
 */
static void
_extract
(
 PDM_isosurface_t *isos,
 int               id_isosurface
)
{
  if (isos->entry_is_part == 0) {
    // Block-distributed
    assert(isos->dist_to_part_computed);
    isos->extract_kind = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    isos->part_method  = PDM_SPLIT_DUAL_WITH_HILBERT;
  }

  if (isosurface_is_nodal(isos)) {
    _extract_nodal(isos, id_isosurface);
  }
  else {
    _extract_ngon(isos, id_isosurface);
  }

  if (isos->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) {
    PDM_malloc(isos->extract_vtx_lnum, isos->iso_n_part, int *);
    for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {
      PDM_extract_part_parent_lnum_get(isos->extrp,
                                       i_part,
                                       PDM_MESH_ENTITY_VTX,
                                       &isos->extract_vtx_lnum[i_part],
                                       PDM_OWNERSHIP_USER);
    }
  }
}


/**
 * \brief Convert nodal multi-sections to ngon
 */
static void
_ngonize
(
  PDM_isosurface_t *isos
)
{
  assert(isos->extrp != NULL);

  PDM_part_mesh_nodal_t *extract_pmn = NULL;
  PDM_extract_part_part_mesh_nodal_get(isos->extrp,
                                       &extract_pmn,
                                       PDM_OWNERSHIP_KEEP);

  isos->extract_pmesh_nodal = extract_pmn;

  /* Inspect nodal sections and check wether we only have simplices */
  int n_section = PDM_part_mesh_nodal_n_section_get(extract_pmn);

  isos->all_simplices = 1;
  for (int i_section = 0; i_section < n_section; i_section++) {
    PDM_Mesh_nodal_elt_t elt_type = PDM_part_mesh_nodal_section_elt_type_get(extract_pmn,
                                                                             i_section);

    if (PDM_Mesh_nodal_elt_dim_get(elt_type) < 2) {
      continue;
    }
    if (elt_type != PDM_MESH_NODAL_TRIA3 &&
        elt_type != PDM_MESH_NODAL_TETRA4) {
      isos->all_simplices = 0;
      break;
    }
  }

  // We assume all ranks have the same sections, so no need for Allreduce

  if (isos->all_simplices) {
    if (isos->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) {// && isos->entry_mesh_dim==2) {
      if (isos->entry_mesh_dim==2) {
        for (int i_part = 0; i_part < isos->n_part; i_part++) {
          PDM_free(isos->extract_face_lnum[i_part]);
        }
      }
      PDM_free(isos->extract_face_lnum);
    }
    if (isos->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL && isos->entry_mesh_dim==3) {
      for (int i_part = 0; i_part < isos->n_part; i_part++) {
        PDM_free(isos->extract_cell_lnum[i_part]);
      }
      PDM_free(isos->extract_cell_lnum);
    }
    PDM_malloc(isos->extract_n_vtx    , isos->iso_n_part, int          );
    PDM_malloc(isos->extract_vtx_coord, isos->iso_n_part, double      *);
    PDM_malloc(isos->extract_vtx_gnum , isos->iso_n_part, PDM_g_num_t *);

    PDM_malloc(isos->extract_n_tri      , isos->iso_n_part, int          );
    PDM_malloc(isos->extract_tri_vtx    , isos->iso_n_part, int         *);
    PDM_malloc(isos->extract_tri_gnum   , isos->iso_n_part, PDM_g_num_t *);
    PDM_malloc(isos->extract_face_lnum  , isos->iso_n_part, int         *);
    PDM_malloc(isos->extract_tri_n_group, isos->iso_n_part, int          );
    PDM_malloc(isos->extract_tri_tag    , isos->iso_n_part, int         *);

    PDM_malloc(isos->extract_n_tetra   , isos->iso_n_part, int          );
    PDM_malloc(isos->extract_tetra_vtx , isos->iso_n_part, int         *);
    PDM_malloc(isos->extract_tetra_gnum, isos->iso_n_part, PDM_g_num_t *);
    PDM_malloc(isos->extract_cell_lnum,  isos->iso_n_part, int         *);

    for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {

      isos->extract_n_vtx    [i_part] = PDM_part_mesh_nodal_n_vtx_get    (extract_pmn, i_part);
      isos->extract_vtx_coord[i_part] = PDM_part_mesh_nodal_vtx_coord_get(extract_pmn, i_part);
      isos->extract_vtx_gnum [i_part] = PDM_part_mesh_nodal_vtx_g_num_get(extract_pmn, i_part);

      // We assume there is at most one section of TETRA and TRIA
      isos->extract_n_tri      [i_part] = 0;
      isos->extract_n_tetra    [i_part] = 0;
      isos->extract_tri_vtx    [i_part] = NULL;
      isos->extract_tri_gnum   [i_part] = NULL;
      isos->extract_face_lnum  [i_part] = NULL;
      isos->extract_tri_n_group[i_part] = 0;
      isos->extract_tri_tag    [i_part] = NULL;
      isos->extract_tetra_vtx  [i_part] = NULL;
      isos->extract_tetra_gnum [i_part] = NULL;
      isos->extract_cell_lnum  [i_part] = NULL;
      for (int i_section = 0; i_section < n_section; i_section++) {
        PDM_Mesh_nodal_elt_t elt_type = PDM_part_mesh_nodal_section_elt_type_get(extract_pmn,
                                                                                 i_section);
        int n_elt = PDM_part_mesh_nodal_section_n_elt_get(extract_pmn,
                                                          i_section,
                                                          i_part);

        int         *connec;
        PDM_g_num_t *gnum;
        int         *parent_num;
        PDM_g_num_t *parent_entity_g_num;
        PDM_part_mesh_nodal_section_std_get(extract_pmn,
                                            i_section,
                                            i_part,
                                            &connec,
                                            &gnum,
                                            &parent_num,
                                            &parent_entity_g_num,
                                            PDM_OWNERSHIP_BAD_VALUE);

        if (elt_type == PDM_MESH_NODAL_TRIA3) {
          isos->extract_n_tri    [i_part] = n_elt;
          isos->extract_tri_vtx  [i_part] = connec;
          isos->extract_tri_gnum [i_part] = parent_entity_g_num;
          isos->extract_face_lnum[i_part] = parent_num;
          if (isos->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) {
            PDM_extract_part_parent_lnum_get(isos->extrp,
                                             i_part,
                                             PDM_MESH_ENTITY_FACE,
                                             &isos->extract_face_lnum[i_part],
                                             PDM_OWNERSHIP_USER);
          }
          isos->extract_tri_n_group[i_part] = _convert_group_info_to_tag(extract_pmn, i_part, n_elt,
                                                                         PDM_GEOMETRY_KIND_SURFACIC,
                                                                         &isos->extract_tri_tag[i_part]);

        }
        if (elt_type == PDM_MESH_NODAL_TETRA4) {
          isos->extract_n_tetra   [i_part] = n_elt;
          isos->extract_tetra_vtx [i_part] = connec;
          isos->extract_tetra_gnum[i_part] = parent_entity_g_num;
          isos->extract_cell_lnum [i_part] = parent_num;
          if (isos->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) {
            PDM_extract_part_parent_lnum_get(isos->extrp,
                                             i_part,
                                             PDM_MESH_ENTITY_CELL,
                                             &isos->extract_cell_lnum[i_part],
                                             PDM_OWNERSHIP_USER);
          }
        }
      } // End loop on sections

    } // End loop on parts

  }
  else {
    // We have elements other than simplices, we need to ngonize
    // PDM_error(__FILE__, __LINE__, 0, "Nodal not implemented yet for elements other than TRIA3 and TETRA4\n");
    isos->entry_mesh_type = 1 * PDM_SIGN(isos->entry_mesh_type); // we are in fact ngon from now on

    PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm = PDM_part_mesh_nodal_to_part_mesh_create(extract_pmn,
                                                                                            PDM_FALSE,
                                                                                            PDM_OWNERSHIP_USER);

    if (isos->entry_mesh_dim == 3) {
      isos->we_have_edges = 1;
      PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm, PDM_CONNECTIVITY_TYPE_CELL_FACE);

      PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm, PDM_CONNECTIVITY_TYPE_FACE_EDGE);

      PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm, PDM_CONNECTIVITY_TYPE_EDGE_VTX);

      PDM_part_mesh_nodal_to_part_mesh_g_nums_enable(pmn_to_pm, PDM_MESH_ENTITY_CELL);
    }
    else {
      isos->we_have_edges = 0;
      PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm, PDM_CONNECTIVITY_TYPE_FACE_VTX);
    }

    PDM_part_mesh_nodal_to_part_mesh_groups_enable(pmn_to_pm, PDM_BOUND_TYPE_FACE);

    PDM_part_mesh_nodal_to_part_mesh_g_nums_enable(pmn_to_pm, PDM_MESH_ENTITY_FACE); // can we avoid that in 3D
    PDM_part_mesh_nodal_to_part_mesh_g_nums_enable(pmn_to_pm, PDM_MESH_ENTITY_VTX);

    PDM_part_mesh_nodal_to_part_mesh_compute(pmn_to_pm);

    PDM_part_mesh_nodal_to_part_mesh_part_mesh_get(pmn_to_pm,
                                                  &isos->extract_pmesh,
                                                   PDM_OWNERSHIP_USER);



    PDM_part_mesh_nodal_to_part_mesh_free(pmn_to_pm);
  }

}


static void
_triangulate
(
  PDM_isosurface_t *isos
)
{
  assert(isos->extrp != NULL);

  PDM_malloc(isos->extract_n_vtx    , isos->iso_n_part, int          );
  PDM_malloc(isos->extract_vtx_coord, isos->iso_n_part, double      *);
  PDM_malloc(isos->extract_vtx_gnum , isos->iso_n_part, PDM_g_num_t *);

  PDM_malloc(isos->extract_n_tri      , isos->iso_n_part, int          );
  PDM_malloc(isos->extract_tri_vtx    , isos->iso_n_part, int         *);
  PDM_malloc(isos->extract_tri_gnum   , isos->iso_n_part, PDM_g_num_t *);
  PDM_malloc(isos->extract_tri_n_group, isos->iso_n_part, int          );
  PDM_malloc(isos->extract_tri_tag    , isos->iso_n_part, int         *);
  if (isos->extract_kind != PDM_EXTRACT_PART_KIND_LOCAL) {
    PDM_malloc(isos->extract_face_lnum, isos->iso_n_part, int         *);
  }
  PDM_malloc(isos->extract_tri_face,    isos->iso_n_part, int         *);

  PDM_malloc(isos->extract_n_tetra   , isos->iso_n_part, int          );
  PDM_malloc(isos->extract_tetra_vtx , isos->iso_n_part, int         *);
  PDM_malloc(isos->extract_tetra_gnum, isos->iso_n_part, PDM_g_num_t *);

  for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {

    isos->extract_n_tetra   [i_part] = 0;
    isos->extract_tetra_vtx [i_part] = NULL;
    isos->extract_tetra_gnum[i_part] = NULL;

    isos->extract_n_vtx[i_part] = PDM_extract_part_vtx_coord_get(isos->extrp,
                                                                 i_part,
                                                                &isos->extract_vtx_coord[i_part],
                                                                 PDM_OWNERSHIP_BAD_VALUE);
    PDM_part_mesh_entity_ln_to_gn_get(isos->extract_pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                     &isos->extract_vtx_gnum[i_part],
                                      PDM_OWNERSHIP_USER);

    PDM_g_num_t *face_parent_gnum = NULL;
    int n_face = PDM_part_mesh_n_entity_get(isos->extract_pmesh,
                                            i_part,
                                            PDM_MESH_ENTITY_FACE);
    PDM_part_mesh_entity_ln_to_gn_get(isos->extract_pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                     &face_parent_gnum,
                                      PDM_OWNERSHIP_KEEP);



    PDM_bool_t owner_face_vtx = PDM_FALSE;
    int *face_vtx_idx = NULL;
    int *face_vtx     = NULL;
    PDM_part_mesh_connectivity_get(isos->extract_pmesh,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                  &face_vtx,
                                  &face_vtx_idx,
                                   PDM_OWNERSHIP_BAD_VALUE);
    if (face_vtx_idx == NULL) {
      owner_face_vtx = PDM_TRUE;

      int *face_edge = NULL;
      PDM_part_mesh_connectivity_get(isos->extract_pmesh,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                    &face_edge,
                                    &face_vtx_idx,
                                     PDM_OWNERSHIP_BAD_VALUE);

      int *edge_vtx_idx = NULL;
      int *edge_vtx     = NULL;
      PDM_part_mesh_connectivity_get(isos->extract_pmesh,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                    &edge_vtx,
                                    &edge_vtx_idx,
                                     PDM_OWNERSHIP_BAD_VALUE);

      PDM_compute_face_vtx_from_face_and_edge(n_face,
                                              face_vtx_idx,
                                              face_edge,
                                              edge_vtx,
                                             &face_vtx);
    }


    int *face_tria_idx = NULL;

    isos->extract_n_tri[i_part] = PDM_triangulate_faces(n_face,
                                                        face_vtx_idx,
                                                        face_vtx,
                                                        isos->extract_vtx_coord[i_part],
                                                       &face_tria_idx,
                                                       &isos->extract_tri_vtx[i_part],
                                                        NULL);

    if (owner_face_vtx == PDM_TRUE) {
      PDM_free(face_vtx);
    }


    int *face_tag = PDM_array_zeros_int(n_face);

    int          n_group_face;
    int         *group_face;
    PDM_g_num_t *group_face_ln_to_gn;
    for (int i_group=0; i_group<isos->n_group_face; ++i_group) {
      PDM_part_mesh_bound_get(isos->extract_pmesh,
                              i_part,
                              i_group,
                              PDM_BOUND_TYPE_FACE,
                             &n_group_face,
                             &group_face,
                             &group_face_ln_to_gn,
                              PDM_OWNERSHIP_KEEP);

      for (int i_elmt=0; i_elmt<n_group_face; ++i_elmt) {
        assert(face_tag[group_face[i_elmt]-1]==0);
        face_tag[group_face[i_elmt]-1] = i_group+1;
      }
    }

    isos->extract_tri_n_group[i_part] = isos->n_group_face;

    int *extract_face_lnum = NULL;
    if (isos->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) {
      extract_face_lnum = isos->extract_face_lnum[i_part];
    }

    PDM_malloc(isos->extract_tri_gnum [i_part], isos->extract_n_tri[i_part], PDM_g_num_t);
    PDM_malloc(isos->extract_face_lnum[i_part], isos->extract_n_tri[i_part], int        );
    PDM_malloc(isos->extract_tri_tag  [i_part], isos->extract_n_tri[i_part], int        );
    PDM_malloc(isos->extract_tri_face [i_part], isos->extract_n_tri[i_part], int        );

    for (int i_face = 0; i_face < n_face; i_face++) {
      for (int i_tri = face_tria_idx[i_face]; i_tri < face_tria_idx[i_face+1]; i_tri++) {
        isos->extract_tri_gnum [i_part][i_tri] = face_parent_gnum[i_face];
        isos->extract_face_lnum[i_part][i_tri] = i_face + 1;
        isos->extract_tri_tag  [i_part][i_tri] = face_tag[i_face];
        isos->extract_tri_face [i_part][i_tri] = i_face;
      }
    }

    if (isos->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) {
      for (int i_tri = 0; i_tri < isos->extract_n_tri[i_part]; i_tri++) {
        int i_face = isos->extract_face_lnum[i_part][i_tri] - 1;
        isos->extract_face_lnum[i_part][i_tri] = extract_face_lnum[i_face];
      }
      PDM_free(extract_face_lnum);
    }

    PDM_free(face_tria_idx);
    PDM_free(face_tag);
  }
}



/**
 * \brief Block-distribute isosurface mesh vertices
 */
static void
_part_to_dist_vtx
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  _isosurface_t *_iso = &isos->isosurfaces[id_iso];

  // > Set vtx in block frame
  PDM_part_to_block_t *ptb_vtx = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                          1.,
                                                          _iso->iso_entity_gnum[PDM_MESH_ENTITY_VTX],
                                                          NULL,
                                                          _iso->iso_n_entity[PDM_MESH_ENTITY_VTX],
                                                          isos->iso_n_part,
                                                          isos->comm);
  // > Get vtx coords
  double *dvtx_coord = NULL;
  PDM_part_to_block_exch(ptb_vtx,
                         sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         3,
                         NULL,
              (void **)  _iso->iso_vtx_coord,
                         NULL,
              (void **) &dvtx_coord);


  // > Get vtx parent
  int **pvtx_parent_strd = NULL;
  PDM_malloc(pvtx_parent_strd, isos->iso_n_part, int *);
  for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
    pvtx_parent_strd[i_part] = PDM_array_new_size_from_idx_int(_iso->iso_entity_parent_idx[PDM_MESH_ENTITY_VTX][i_part], _iso->iso_n_entity[PDM_MESH_ENTITY_VTX][i_part]);
  }

  int          *dvtx_parent_strd = NULL;
  PDM_g_num_t  *dvtx_parent_gnum = NULL;
  PDM_part_to_block_exch(ptb_vtx,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         pvtx_parent_strd,
              (void **)  _iso->iso_entity_parent_gnum[PDM_MESH_ENTITY_VTX],
                        &dvtx_parent_strd,
              (void **) &dvtx_parent_gnum);
  PDM_free(dvtx_parent_strd);

  double  *dvtx_parent_weight = NULL;
  PDM_part_to_block_exch(ptb_vtx,
                         sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         pvtx_parent_strd,
              (void **)  _iso->iso_entity_parent_wght[PDM_MESH_ENTITY_VTX],
                        &dvtx_parent_strd,
              (void **) &dvtx_parent_weight);


  int dn_vtx = PDM_part_to_block_n_elt_block_get(ptb_vtx);
  int *dvtx_parent_idx = PDM_array_new_idx_from_sizes_int(dvtx_parent_strd, dn_vtx);
  for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
    PDM_free(pvtx_parent_strd[i_part]);
  }
  PDM_free(pvtx_parent_strd);
  PDM_free(dvtx_parent_strd);

  _iso->iso_dn_entity          [PDM_MESH_ENTITY_VTX] = dn_vtx;
  _iso->iso_dentity_parent_idx [PDM_MESH_ENTITY_VTX] = dvtx_parent_idx;
  _iso->iso_dentity_parent_gnum[PDM_MESH_ENTITY_VTX] = dvtx_parent_gnum;
  _iso->iso_dentity_parent_wght[PDM_MESH_ENTITY_VTX] = dvtx_parent_weight;
  _iso->iso_dvtx_coord         = dvtx_coord;
  PDM_part_to_block_free(ptb_vtx);

}


/**
 * \brief Block-distribute isosurface mesh vertices
 */
static void
_part_to_dist_edge_group
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  int dbg = 0;

  _isosurface_t *_iso = &isos->isosurfaces[id_iso];

  // > Allocate tmp result
  int         *_iso_dedge_group_idx  = NULL;
  PDM_g_num_t *_iso_dedge_group_gnum = NULL;
  PDM_malloc(_iso_dedge_group_idx , _iso->iso_n_edge_group+1, int        );
  PDM_malloc(_iso_dedge_group_gnum, _iso->iso_n_edge_group  , PDM_g_num_t);
  _iso_dedge_group_idx[0] = 0;


  // > Exchange info with a part_to_block for each group
  int          *n_elt_group     = NULL;
  PDM_g_num_t **    _group_gnum = NULL;
  PDM_g_num_t **_elt_group_gnum = NULL;
  PDM_malloc(n_elt_group    , isos->iso_n_part, int          );
  PDM_malloc(_elt_group_gnum, isos->iso_n_part, PDM_g_num_t *);
  PDM_malloc(    _group_gnum, isos->iso_n_part, PDM_g_num_t *);

  for (int i_group=0; i_group<_iso->iso_n_edge_group; ++i_group) {

    // > Prepare ptb
    for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
      int i_beg_group = _iso->iso_edge_group_idx[i_part][i_group  ];
      int i_end_group = _iso->iso_edge_group_idx[i_part][i_group+1];
      n_elt_group[i_part] = i_end_group-i_beg_group;

      // > Prepare element group gnum
      PDM_malloc(_elt_group_gnum[i_part], n_elt_group[i_part], PDM_g_num_t);

      int i_write = 0;
      for (int i_elmt=i_beg_group; i_elmt<i_end_group; ++i_elmt) {
        int elmt_lnum = _iso->iso_edge_group_lnum[i_part][i_elmt];
        _elt_group_gnum[i_part][i_write++] = _iso->iso_entity_gnum[PDM_MESH_ENTITY_EDGE][i_part][elmt_lnum-1];
      }

      _group_gnum[i_part] = &_iso->iso_edge_group_gnum[i_part][i_beg_group];
    } // End loop on partitions


    // > Exchange with part_to_block
    PDM_part_to_block_t *ptb_group = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                              PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                              1.,
                                                              _group_gnum,
                                                              NULL,
                                                              n_elt_group,
                                                              isos->iso_n_part,
                                                              isos->comm);

    int dn_elt_group = PDM_part_to_block_n_elt_block_get(ptb_group);
    _iso_dedge_group_idx[i_group+1] = _iso_dedge_group_idx[i_group] + dn_elt_group;

    PDM_realloc(_iso_dedge_group_gnum,
                _iso_dedge_group_gnum,
                _iso_dedge_group_idx[i_group+1],
                PDM_g_num_t);

    PDM_g_num_t *rvcd_gnum = NULL;
    PDM_part_to_block_exch(ptb_group,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST_INTERLACED,
                           1,
                           NULL,
                (void **)  _elt_group_gnum,
                           NULL,
                (void **) &rvcd_gnum);
    memcpy(&_iso_dedge_group_gnum[_iso_dedge_group_idx[i_group]], rvcd_gnum, dn_elt_group * sizeof(PDM_g_num_t));

    PDM_part_to_block_free(ptb_group);
    PDM_free(rvcd_gnum);
    for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
      PDM_free(_elt_group_gnum[i_part]);
    }
  } // End loop on groups
  PDM_free(n_elt_group);
  PDM_free(_elt_group_gnum);
  PDM_free(    _group_gnum);


  /**
   * Set result
   */
  PDM_free(_group_gnum);


  if (dbg==1) {
    int size_group_gnum = _iso_dedge_group_idx[_iso->iso_n_edge_group];
    log_trace("iso_n_edge_group = %d\n", _iso->iso_n_edge_group);
    PDM_log_trace_array_int (_iso_dedge_group_idx , _iso->iso_n_edge_group+1, "iso_dedge_group_idx  ::");
    PDM_log_trace_array_long(_iso_dedge_group_gnum, size_group_gnum         , "iso_dedge_group_gnum ::");
  }


  /**
   * Set result
   */
  _iso->iso_dedge_group_idx  = _iso_dedge_group_idx;
  _iso->iso_dedge_group_gnum = _iso_dedge_group_gnum;
}


/**
 * \brief Block-distribute isosurface mesh vertices
 */
static void
_part_to_dist_elt
(
  PDM_isosurface_t *isos,
  int               id_iso,
  int              *n_elt,
  PDM_g_num_t     **elt_gnum,
  int             **elt_vtx_idx,
  int             **elt_vtx,
  int             **elt_parent_idx,
  PDM_g_num_t     **elt_parent_gnum,
  double          **elt_parent_wght,
  int              *out_dn_elt,
  int             **out_delt_vtx_idx,
  PDM_g_num_t     **out_delt_vtx,
  int             **out_delt_parent_idx,
  PDM_g_num_t     **out_delt_parent_gnum,
  double          **out_delt_parent_wght
)
{
  /**
   * TODO: if unify_parent
   * PDM_PART_TO_BLOCK_POST_CLEANUP -> PDM_PART_TO_BLOCK_POST_MERGE
   * Go through received data :
   *   - select one connectivity example (should be the same on all procs)
   *   - unique parents (must have various from procs)
   * For ngon algo one data should be received for each entity
   */

  int dbg = 0;

  _isosurface_t *_iso = &isos->isosurfaces[id_iso];

  // > Set elt in block frame
  PDM_part_to_block_t *ptb_elt = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                          1.,
                                                          elt_gnum,
                                                          NULL,
                                                          n_elt,
                                                          isos->iso_n_part,
                                                          isos->comm);

  /**
   * Count how many received for block elements
   */
  int  dn_elt           = PDM_part_to_block_n_elt_block_get     (ptb_elt);
  int *block_gnum_count = PDM_part_to_block_block_gnum_count_get(ptb_elt);
  if (dbg==1) {
    for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
      PDM_log_trace_array_long(elt_gnum[i_part], n_elt[i_part], "elt_gnum ::");
    }

    log_trace("dn_elt = %d\n", dn_elt);
    PDM_log_trace_array_int(block_gnum_count, dn_elt, "block_gnum_count ::");
  }
  *out_dn_elt = dn_elt;


  /**
   * Get element connectivity
   */

  // > Translate connectivity into gid
  int         **_elt_vtx_strd = NULL;
  PDM_g_num_t **_elt_vtx      = NULL;
  PDM_malloc(_elt_vtx_strd, isos->iso_n_part, int         *);
  PDM_malloc(_elt_vtx     , isos->iso_n_part, PDM_g_num_t *);

  for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
    int size_elt_vtx = 0;
    if (elt_vtx_idx!=NULL) { // Face
      size_elt_vtx = elt_vtx_idx[i_part][n_elt[i_part]];
      _elt_vtx_strd[i_part] = PDM_array_new_size_from_idx_int(elt_vtx_idx[i_part], n_elt[i_part]);
    }
    else { // Edge
      size_elt_vtx = 2*n_elt[i_part];
      _elt_vtx_strd[i_part] = PDM_array_const_int(n_elt[i_part], 2);
    }
    PDM_malloc(_elt_vtx[i_part], size_elt_vtx, PDM_g_num_t);
    for (int i_elt=0; i_elt<size_elt_vtx; ++i_elt) {
      int elt_lnum = elt_vtx[i_part][i_elt];
      _elt_vtx[i_part][i_elt] = _iso->iso_entity_gnum[PDM_MESH_ENTITY_VTX][i_part][elt_lnum-1];
    }
    int idx = 0;
    for (int i_elt=0; i_elt<n_elt[i_part]; ++i_elt) {
      // log_trace("elt "PDM_FMT_G_NUM" : ", elt_gnum[i_part][i_elt]);
      // PDM_log_trace_array_long(_elt_vtx[i_part] + idx, _elt_vtx_strd[i_part][i_elt], "");
      idx += _elt_vtx_strd[i_part][i_elt];
    }
  }

  // > Exchange
  int         *delt_vtx_strd = NULL;
  PDM_g_num_t *delt_vtx      = NULL;
  int s_block_data = PDM_part_to_block_exch(ptb_elt,
                                            sizeof(PDM_g_num_t),
                                            PDM_STRIDE_VAR_INTERLACED,
                                            1,
                                            _elt_vtx_strd,
                                 (void **)  _elt_vtx,
                                           &delt_vtx_strd,
                                 (void **) &delt_vtx);

  if (dbg) {
    log_trace("s_block_data = %d\n", s_block_data);
    int rcvd_size = 0;
    for (int i_rcvd=0; i_rcvd<dn_elt; ++i_rcvd) {
      rcvd_size+=delt_vtx_strd[i_rcvd];
    }
    PDM_log_trace_array_int (delt_vtx_strd, dn_elt   , "delt_vtx_strd :: ");
    PDM_log_trace_array_long(delt_vtx     , rcvd_size, "delt_vtx      :: ");
  }
  if (elt_vtx_idx!=NULL) { // Face
    *out_delt_vtx_idx = PDM_array_new_idx_from_sizes_int(delt_vtx_strd, dn_elt);
  }
  *out_delt_vtx = delt_vtx;

  for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
    PDM_free(_elt_vtx_strd[i_part]);
    PDM_free(_elt_vtx     [i_part]);
  }
  PDM_free(_elt_vtx_strd);
  PDM_free(_elt_vtx);
  PDM_free(delt_vtx_strd);


  /**
   * Get parent
   */

  // > Prepare stride
  int **_elt_parent_strd = NULL;
  PDM_malloc(_elt_parent_strd, isos->iso_n_part, int *);
  for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
    _elt_parent_strd[i_part] = PDM_array_new_size_from_idx_int(elt_parent_idx[i_part], n_elt[i_part]);
  }

  // > Exchange
  int         *delt_parent_strd = NULL;
  PDM_g_num_t *delt_parent_gnum= NULL;
  s_block_data = PDM_part_to_block_exch(ptb_elt,
                                        sizeof(PDM_g_num_t),
                                        PDM_STRIDE_VAR_INTERLACED,
                                        1,
                                       _elt_parent_strd,
                             (void **)  elt_parent_gnum,
                                       &delt_parent_strd,
                             (void **) &delt_parent_gnum);
  PDM_free(delt_parent_strd);


  double  *delt_parent_weight = NULL;
  PDM_part_to_block_exch(ptb_elt,
                         sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         _elt_parent_strd,
              (void **)  elt_parent_wght,
                        &delt_parent_strd,
              (void **) &delt_parent_weight);

  if (dbg) {
    log_trace("s_block_data = %d\n", s_block_data);
    int rcvd_size = 0;
    for (int i_rcvd=0; i_rcvd<dn_elt; ++i_rcvd) {
      rcvd_size+=delt_parent_strd[i_rcvd];
    }
    PDM_log_trace_array_int (delt_parent_strd, dn_elt   , "delt_parent_strd :: ");
    PDM_log_trace_array_long(delt_parent_gnum, rcvd_size, "delt_parent_gnum :: ");
  }

  for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
    PDM_free(_elt_parent_strd[i_part]);
  }
  PDM_free(_elt_parent_strd);



  *out_delt_parent_idx  = PDM_array_new_idx_from_sizes_int(delt_parent_strd, dn_elt);
  *out_delt_parent_gnum = delt_parent_gnum;
  *out_delt_parent_wght = delt_parent_weight;

  PDM_free(delt_parent_strd);

  PDM_part_to_block_free(ptb_elt);
}


// TODO: factorize as partgroup_to_distgroup?
static void
_part_to_dist_isovalue_entity
(
  PDM_MPI_Comm   comm,
  int            n_part,
  int            n_isovalues,
  int          **isovalue_entity_idx,
  PDM_g_num_t  **entity_gnum,
  int          **disovalue_entity_idx,
  PDM_g_num_t  **disovalue_entity
)
{
  PDM_malloc(*disovalue_entity_idx, n_isovalues+1, int);
  int *_disovalue_entity_idx = *disovalue_entity_idx;
  _disovalue_entity_idx[0] = 0;

  PDM_part_to_block_t **ptb = NULL;
  PDM_malloc(ptb, n_isovalues, PDM_part_to_block_t *);

  for (int i_isovalue = 0; i_isovalue < n_isovalues; i_isovalue++) {

    int          *pn_entity = NULL;
    PDM_g_num_t **part_gnum = NULL;
    double      **weight    = NULL;
    PDM_malloc(pn_entity, n_part, int          );
    PDM_malloc(part_gnum, n_part, PDM_g_num_t *);
    PDM_malloc(weight,    n_part, double      *);

    for (int i_part = 0; i_part < n_part; i_part++) {
      pn_entity[i_part] = isovalue_entity_idx[i_part][i_isovalue+1] - isovalue_entity_idx[i_part][i_isovalue];

      PDM_malloc(part_gnum[i_part], pn_entity[i_part], PDM_g_num_t);
      PDM_malloc(weight   [i_part], pn_entity[i_part], double     );
      for (int i = 0; i < pn_entity[i_part]; i++) {
        int i_entity = isovalue_entity_idx[i_part][i_isovalue] + i;
        part_gnum[i_part][i] = entity_gnum[i_part][i_entity];
        weight   [i_part][i] = 1.;
      }
    }

    ptb[i_isovalue] = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                               PDM_PART_TO_BLOCK_POST_CLEANUP,
                                               1.,
                                               part_gnum,
                                               weight,
                                               pn_entity,
                                               n_part,
                                               comm);

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_free(part_gnum[i_part]);
      PDM_free(weight   [i_part]);
    }
    PDM_free(pn_entity);
    PDM_free(part_gnum);
    PDM_free(weight   );

    int n_elt_block = PDM_part_to_block_n_elt_block_get(ptb[i_isovalue]);

    _disovalue_entity_idx[i_isovalue+1] = _disovalue_entity_idx[i_isovalue] + n_elt_block;

  } // End loop on isovalues

  PDM_malloc(*disovalue_entity, _disovalue_entity_idx[n_isovalues], PDM_g_num_t);
  PDM_g_num_t *_disovalue_entity = *disovalue_entity;

  for (int i_isovalue = 0; i_isovalue < n_isovalues; i_isovalue++) {

    int          n_elt_block = PDM_part_to_block_n_elt_block_get(ptb[i_isovalue]);
    PDM_g_num_t *block_gnum  = PDM_part_to_block_block_gnum_get (ptb[i_isovalue]);

    for (int i = 0; i < n_elt_block; i++) {
      _disovalue_entity[_disovalue_entity_idx[i_isovalue] + i] = block_gnum[i];
    }

    PDM_part_to_block_free(ptb[i_isovalue]);
  }
  PDM_free(ptb);
}


/**
 * \brief Block-distribute isosurface mesh
 */
static void
_part_to_dist
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  _isosurface_t *_iso = &isos->isosurfaces[id_iso];

  /** Vertices */
  _part_to_dist_vtx(isos, id_iso);
  _part_to_dist_isovalue_entity(isos->comm,
                                isos->iso_n_part,
                                _iso->n_isovalues,
                                _iso->isovalue_entity_idx  [PDM_MESH_ENTITY_VTX],
                                _iso->iso_entity_gnum      [PDM_MESH_ENTITY_VTX],
                                &_iso->disovalue_entity_idx[PDM_MESH_ENTITY_VTX],
                                &_iso->disovalue_entity    [PDM_MESH_ENTITY_VTX]);
  _iso->iso_owner_dvtx_coord                            = PDM_OWNERSHIP_KEEP;
  _iso->iso_owner_dparent_idx     [PDM_MESH_ENTITY_VTX] = PDM_OWNERSHIP_KEEP;
  _iso->iso_owner_dparent         [PDM_MESH_ENTITY_VTX] = PDM_OWNERSHIP_KEEP;
  _iso->iso_owner_dparent_wght    [PDM_MESH_ENTITY_VTX] = PDM_OWNERSHIP_KEEP;
  _iso->iso_owner_disovalue_entity[PDM_MESH_ENTITY_VTX] = PDM_OWNERSHIP_KEEP;

  /** Edges */
  _part_to_dist_elt(isos,
                    id_iso,
                    _iso->iso_n_entity           [PDM_MESH_ENTITY_EDGE          ],
                    _iso->iso_entity_gnum        [PDM_MESH_ENTITY_EDGE          ],
                    NULL,
                    _iso->iso_connec             [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                    _iso->iso_entity_parent_idx  [PDM_MESH_ENTITY_EDGE          ],
                    _iso->iso_entity_parent_gnum [PDM_MESH_ENTITY_EDGE          ],
                    _iso->iso_entity_parent_wght [PDM_MESH_ENTITY_EDGE          ],
                   &_iso->iso_dn_entity          [PDM_MESH_ENTITY_EDGE          ],
                    NULL,
                   &_iso->iso_dconnec            [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                   &_iso->iso_dentity_parent_idx [PDM_MESH_ENTITY_EDGE          ],
                   &_iso->iso_dentity_parent_gnum[PDM_MESH_ENTITY_EDGE          ],
                   &_iso->iso_dentity_parent_wght[PDM_MESH_ENTITY_EDGE          ]);
  _part_to_dist_isovalue_entity(isos->comm,
                                isos->iso_n_part,
                                _iso->n_isovalues,
                                _iso->isovalue_entity_idx  [PDM_MESH_ENTITY_EDGE],
                                _iso->iso_entity_gnum      [PDM_MESH_ENTITY_EDGE],
                                &_iso->disovalue_entity_idx[PDM_MESH_ENTITY_EDGE],
                                &_iso->disovalue_entity    [PDM_MESH_ENTITY_EDGE]);

  _iso->iso_owner_dconnec         [PDM_CONNECTIVITY_TYPE_EDGE_VTX] = PDM_OWNERSHIP_KEEP;
  _iso->iso_owner_dparent_idx     [PDM_MESH_ENTITY_EDGE          ] = PDM_OWNERSHIP_KEEP;
  _iso->iso_owner_dparent         [PDM_MESH_ENTITY_EDGE          ] = PDM_OWNERSHIP_KEEP;
  _iso->iso_owner_dparent_wght    [PDM_MESH_ENTITY_EDGE          ] = PDM_OWNERSHIP_KEEP;
  _iso->iso_owner_disovalue_entity[PDM_MESH_ENTITY_EDGE          ] = PDM_OWNERSHIP_KEEP;

  _iso->iso_owner_dedge_bnd = PDM_OWNERSHIP_KEEP;


  if (isos->entry_mesh_dim==3) {
    /** Edge groups */
    _part_to_dist_edge_group(isos, id_iso);


    /** Faces */
    _part_to_dist_elt(isos,
                      id_iso,
                      _iso->iso_n_entity           [PDM_MESH_ENTITY_FACE          ],
                      _iso->iso_entity_gnum        [PDM_MESH_ENTITY_FACE          ],
                      _iso->iso_connec_idx         [PDM_CONNECTIVITY_TYPE_FACE_VTX],
                      _iso->iso_connec             [PDM_CONNECTIVITY_TYPE_FACE_VTX],
                      _iso->iso_entity_parent_idx  [PDM_MESH_ENTITY_FACE          ],
                      _iso->iso_entity_parent_gnum [PDM_MESH_ENTITY_FACE          ],
                      _iso->iso_entity_parent_wght [PDM_MESH_ENTITY_FACE          ],
                     &_iso->iso_dn_entity          [PDM_MESH_ENTITY_FACE          ],
                     &_iso->iso_dconnec_idx        [PDM_CONNECTIVITY_TYPE_FACE_VTX],
                     &_iso->iso_dconnec            [PDM_CONNECTIVITY_TYPE_FACE_VTX],
                     &_iso->iso_dentity_parent_idx [PDM_MESH_ENTITY_FACE          ],
                     &_iso->iso_dentity_parent_gnum[PDM_MESH_ENTITY_FACE          ],
                     &_iso->iso_dentity_parent_wght[PDM_MESH_ENTITY_FACE          ]);
      _part_to_dist_isovalue_entity(isos->comm,
                                    isos->iso_n_part,
                                    _iso->n_isovalues,
                                    _iso->isovalue_entity_idx  [PDM_MESH_ENTITY_FACE],
                                    _iso->iso_entity_gnum      [PDM_MESH_ENTITY_FACE],
                                    &_iso->disovalue_entity_idx[PDM_MESH_ENTITY_FACE],
                                    &_iso->disovalue_entity    [PDM_MESH_ENTITY_FACE]);

    _iso->iso_owner_dconnec         [PDM_CONNECTIVITY_TYPE_FACE_VTX] = PDM_OWNERSHIP_KEEP;
    _iso->iso_owner_dparent_idx     [PDM_MESH_ENTITY_FACE          ] = PDM_OWNERSHIP_KEEP;
    _iso->iso_owner_dparent         [PDM_MESH_ENTITY_FACE          ] = PDM_OWNERSHIP_KEEP;
    _iso->iso_owner_dparent_wght    [PDM_MESH_ENTITY_FACE          ] = PDM_OWNERSHIP_KEEP;
    _iso->iso_owner_disovalue_entity[PDM_MESH_ENTITY_FACE          ] = PDM_OWNERSHIP_KEEP;
  }

}


/* Build part_to_part to link isosurface entities with their 'parent' source entities in user frame */
static void
_build_ptp_part
(
  PDM_isosurface_t    *isos,
  int                  id_iso,
  PDM_mesh_entities_t  entity_type
)
{
  _isosurface_t *_iso = &isos->isosurfaces[id_iso];

  if (_iso->compute_ptp[entity_type] == PDM_FALSE) {
    return;
  }

  if (entity_type == PDM_MESH_ENTITY_FACE && isos->entry_mesh_dim < 3) {
    return;
  }

  assert(isos->entry_is_part == 1);
  assert(isos->extract_kind != PDM_EXTRACT_PART_KIND_LOCAL);

  PDM_extract_part_t *extrp = isos->extrp;
  assert(extrp != NULL);

  int          *n_entity           = _iso->iso_n_entity          [entity_type];
  PDM_g_num_t **entity_gnum        = _iso->iso_entity_gnum       [entity_type];
  int         **entity_parent_idx  = _iso->iso_entity_parent_idx [entity_type];
  int         **entity_parent_lnum = _iso->iso_entity_parent_lnum[entity_type];
  int          *n_parent           = NULL;

  if (isosurface_is_nodal(isos)) {
    PDM_malloc(n_parent, isos->n_part, int);
  }

  PDM_mesh_entities_t parent_entity_type = PDM_MESH_ENTITY_MAX;
  switch (entity_type) {
    case PDM_MESH_ENTITY_VTX: {
      parent_entity_type = PDM_MESH_ENTITY_VTX;
      if (isosurface_is_nodal(isos)) {
        for (int i_part = 0; i_part < isos->n_part; i_part++) {
          n_parent[i_part] = PDM_part_mesh_nodal_n_vtx_get(isos->pmesh_nodal, i_part);
        }
      }
      else {
        n_parent = isos->n_vtx;
      }
      break;
    }

    case PDM_MESH_ENTITY_EDGE: {
      parent_entity_type = PDM_MESH_ENTITY_FACE;
      if (isosurface_is_nodal(isos)) {
        for (int i_part = 0; i_part < isos->n_part; i_part++) {
          n_parent[i_part] = PDM_part_mesh_nodal_n_elmts_get(isos->pmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, i_part);
        }
      }
      else {
        n_parent = isos->n_face;
      }
      break;
    }

    case PDM_MESH_ENTITY_FACE: {
      parent_entity_type = PDM_MESH_ENTITY_CELL;
      if (isosurface_is_nodal(isos)) {
        for (int i_part = 0; i_part < isos->n_part; i_part++) {
          n_parent[i_part] = PDM_part_mesh_nodal_n_elmts_get(isos->pmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, i_part);
        }
      }
      else {
        n_parent = isos->n_cell;
      }
      break;
    }

    default : {
      PDM_error(__FILE__, __LINE__, 0, "Invalid entity_type %d\n", entity_type);
    }
  }


  /* Get init location of extracted parent entities */
  int **entity_parent_triplet_idx   = NULL;
  int **entity_parent_init_location = NULL;
  PDM_malloc(entity_parent_triplet_idx,   isos->iso_n_part, int *);
  PDM_malloc(entity_parent_init_location, isos->iso_n_part, int *);
  for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {

    int *parent_init_location = NULL;
    PDM_extract_part_init_location_get(isos->extrp,
                                       i_part,
                                       parent_entity_type,
                                      &parent_init_location,
                                       PDM_OWNERSHIP_BAD_VALUE);

    entity_parent_triplet_idx[i_part] = PDM_array_new_idx_from_const_stride_int(3, entity_parent_idx[i_part][n_entity[i_part]]);
    PDM_malloc(entity_parent_init_location[i_part], entity_parent_idx[i_part][n_entity[i_part]] * 3, int);

    for (int i_entity = 0; i_entity < n_entity[i_part]; i_entity++) {

      for (int i = entity_parent_idx[i_part][i_entity]; i < entity_parent_idx[i_part][i_entity+1]; i++) {
        int i_parent = entity_parent_lnum[i_part][i] - 1;
        for (int j = 0; j < 3; j++) {
          entity_parent_init_location[i_part][3*i+j] = parent_init_location[3*i_parent+j];
        }
      }

    }
  }


  /* Create ptp from isosurface entities to source entities */
  _iso->iso_ptp[entity_type] = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) entity_gnum,
                                                                         (const int         * ) n_entity,
                                                                                                isos->iso_n_part,
                                                                         (const int         * ) n_parent,
                                                                                                isos->n_part,
                                                                         (const int         **) entity_parent_idx,
                                                                         (const int         **) entity_parent_triplet_idx,
                                                                         (const int         **) entity_parent_init_location,
                                                                                                isos->comm);
  for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {
    PDM_free(entity_parent_triplet_idx  [i_part]);
    PDM_free(entity_parent_init_location[i_part]);
  }
  PDM_free(entity_parent_triplet_idx  );
  PDM_free(entity_parent_init_location);

  _iso->iso_owner_ptp[entity_type] = PDM_OWNERSHIP_KEEP;

  if (isosurface_is_nodal(isos)) {
    PDM_free(n_parent); // no worries, it is deep-copied in part_to_part creation
  }
}


/* Build part_to_part (distributed output) to link isosurface entities
   with their 'parent' source entities in user frame */
static void
_build_ptp_dist
(
  PDM_isosurface_t    *isos,
  int                  id_iso,
  PDM_mesh_entities_t  entity_type
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(isos->comm, &i_rank);
  PDM_MPI_Comm_size(isos->comm, &n_rank);

  _isosurface_t *_iso = &isos->isosurfaces[id_iso];

  if (_iso->compute_ptp[entity_type] == PDM_FALSE) {
    return;
  }

  if (entity_type == PDM_MESH_ENTITY_FACE && isos->entry_mesh_dim < 3) {
    return;
  }

  int          iso_n_entity    = _iso->iso_dn_entity          [entity_type];
  int         *iso_parent_idx  = _iso->iso_dentity_parent_idx [entity_type];
  PDM_g_num_t *iso_parent_gnum = _iso->iso_dentity_parent_gnum[entity_type];

  // > Isosurface utils for ptp
  PDM_g_num_t *iso_distrib = NULL;
  iso_distrib = PDM_compute_entity_distribution(isos->comm, iso_n_entity);
  PDM_g_num_t *iso_gnum = PDM_array_new_arange_gnum(iso_distrib[i_rank]+1, iso_distrib[i_rank+1]+1, 1);
  PDM_free(iso_distrib);


  // > Entry mesh utils for ptp
  PDM_g_num_t *parent_distrib = NULL;
  switch (entity_type) {
    case PDM_MESH_ENTITY_VTX: {
      if (isosurface_is_nodal(isos)) {
        parent_distrib = PDM_dmesh_nodal_vtx_distrib_get(isos->dmesh_nodal);
      }
      else {
        parent_distrib = isos->distrib_vtx;
      }
      break;
    }

    case PDM_MESH_ENTITY_EDGE: {
      assert(!isosurface_is_nodal(isos));
      parent_distrib = isos->distrib_face;
      break;
    }

    case PDM_MESH_ENTITY_FACE: {
      assert(!isosurface_is_nodal(isos));
      parent_distrib = isos->distrib_cell;
      break;
    }

    default : {
      PDM_error(__FILE__, __LINE__, 0, "Invalid entity_type %d\n", entity_type);
    }
  }

  /**
   * Transform parent_gnum element onto triplet
   */

  // > Multiply by 3 iso_parent_idx because of ptp_from_triplet
  int *iso_parent_triplet_idx = NULL;
  PDM_malloc(iso_parent_triplet_idx, iso_n_entity+1, int);
  for (int i_entity=0; i_entity<iso_n_entity+1; ++i_entity) {
    iso_parent_triplet_idx[i_entity] = 3*iso_parent_idx[i_entity];
  }

  // > Build triplet from gnum
  int  n_elt2           = parent_distrib[i_rank+1]-parent_distrib[i_rank];
  int  n_parent         = iso_parent_idx[iso_n_entity];
  int *iso_parent_trplt = NULL;
  PDM_malloc(iso_parent_trplt, 3*n_parent, int);
  for (int i_parent=0; i_parent<n_parent; ++i_parent) {
    PDM_g_num_t gnum = iso_parent_gnum[i_parent]-1;
    int iproc = PDM_binary_search_gap_long(gnum,
                                           parent_distrib,
                                           n_rank + 1);
    int lnum  = gnum - parent_distrib[iproc];
    iso_parent_trplt[3*i_parent  ] = iproc;
    iso_parent_trplt[3*i_parent+1] = 0;
    iso_parent_trplt[3*i_parent+2] = lnum;
  }

  _iso->iso_ptp[entity_type] = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **)&iso_gnum,
                                                                         (const int         * )&iso_n_entity,
                                                                                                1,
                                                                         (const int         * )&n_elt2,
                                                                                                1,
                                                                         (const int         **)&iso_parent_triplet_idx,
                                                                         (const int         **) NULL,
                                                                         (const int         **)&iso_parent_trplt,
                                                                                                isos->comm);
  _iso->iso_owner_ptp[entity_type] = PDM_OWNERSHIP_KEEP;

  PDM_free(iso_gnum);
  PDM_free(iso_parent_triplet_idx);
  PDM_free(iso_parent_trplt);
}


/* Build part_to_part (distributed output) to link isosurface entities
   with their 'parent' source entities in user frame (dmesh_nodal) */
static void
_build_ptp_dist_nodal
(
  PDM_isosurface_t    *isos,
  int                  id_iso,
  PDM_mesh_entities_t  entity_type
)
{
  _isosurface_t *_iso = &isos->isosurfaces[id_iso];
  if (_iso->compute_ptp[entity_type] == PDM_FALSE) {
    return;
  }

  if (entity_type == PDM_MESH_ENTITY_VTX) {
    _build_ptp_dist(isos, id_iso, entity_type);
    return;
  }

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(isos->comm, &i_rank);
  PDM_MPI_Comm_size(isos->comm, &n_rank);


  assert(isos->dmesh_nodal != NULL);

  if (entity_type == PDM_MESH_ENTITY_FACE && isos->entry_mesh_dim < 3) {
    return;
  }

  int          iso_n_entity    = _iso->iso_dn_entity          [entity_type];
  int         *iso_parent_idx  = _iso->iso_dentity_parent_idx [entity_type];
  PDM_g_num_t *iso_parent_gnum = _iso->iso_dentity_parent_gnum[entity_type];

  // > Isosurface utils for ptp
  PDM_g_num_t *iso_distrib = NULL;
  iso_distrib = PDM_compute_entity_distribution(isos->comm, iso_n_entity);
  PDM_g_num_t *iso_gnum = PDM_array_new_arange_gnum(iso_distrib[i_rank]+1, iso_distrib[i_rank+1]+1, 1);
  PDM_free(iso_distrib);

  PDM_geometry_kind_t geom_kind_parent = PDM_GEOMETRY_KIND_MAX;

  if (entity_type == PDM_MESH_ENTITY_EDGE) {
    geom_kind_parent = PDM_GEOMETRY_KIND_SURFACIC;
  }
  else if (entity_type == PDM_MESH_ENTITY_FACE) {
    geom_kind_parent = PDM_GEOMETRY_KIND_VOLUMIC;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Invalid entity_type %d\n", entity_type);
  }
  int n_section = PDM_DMesh_nodal_n_section_get(isos->dmesh_nodal, geom_kind_parent);

  PDM_g_num_t *shifted_section_distrib = NULL;
  PDM_malloc(shifted_section_distrib, n_section+1, PDM_g_num_t);
  shifted_section_distrib[0] = 0;

  // int *sections_id = PDM_DMesh_nodal_sections_id_get(isos->dmesh_nodal, geom_kind);
  int n_elt2 = 0;
  for (int i_section = 0; i_section < n_section; i_section++) {
    const PDM_g_num_t *section_distrib = PDM_DMesh_nodal_distrib_section_get(isos->dmesh_nodal,
                                                                             geom_kind_parent,
                                                                             i_section); // section_id[i_section] ?

    shifted_section_distrib[i_section+1] = shifted_section_distrib[i_section] + section_distrib[n_rank];

    n_elt2 += section_distrib[i_rank+1] - section_distrib[i_rank];
  }


  /**
   * Transform parent_gnum element into triplet
   */

  // > Multiply by 3 iso_parent_idx because of ptp_from_triplet
  int *iso_parent_triplet_idx = NULL;
  PDM_malloc(iso_parent_triplet_idx, iso_n_entity+1, int);
  for (int i_entity=0; i_entity<iso_n_entity+1; ++i_entity) {
    iso_parent_triplet_idx[i_entity] = 3*iso_parent_idx[i_entity];
  }

  // > Build triplet from gnum
  int  n_parent         = iso_parent_idx[iso_n_entity];
  int *iso_parent_trplt = NULL;
  PDM_malloc(iso_parent_trplt, 3*n_parent, int);
  for (int i_parent=0; i_parent<n_parent; ++i_parent) {
    // First, get section
    PDM_g_num_t gnum = iso_parent_gnum[i_parent]-1;
    int i_section = PDM_binary_search_gap_long(gnum,
                                               shifted_section_distrib,
                                               n_section + 1);

    // Second, get proc
    const PDM_g_num_t *section_distrib = PDM_DMesh_nodal_distrib_section_get(isos->dmesh_nodal,
                                                                             geom_kind_parent,
                                                                             i_section); // section_id[i_section] ?
    gnum -= shifted_section_distrib[i_section];
    int iproc = PDM_binary_search_gap_long(gnum,
                                           section_distrib,
                                           n_rank + 1);

    // Third, get lnum in that section of that proc
    int lnum = gnum - section_distrib[iproc];
    for (int j_section = 0; j_section < i_section; j_section++) {
      const PDM_g_num_t *_section_distrib = PDM_DMesh_nodal_distrib_section_get(isos->dmesh_nodal,
                                                                                geom_kind_parent,
                                                                                i_section); // section_id[i_section] ?
      lnum += _section_distrib[iproc+1] - _section_distrib[iproc];
    }

    iso_parent_trplt[3*i_parent  ] = iproc;
    iso_parent_trplt[3*i_parent+1] = 0;
    iso_parent_trplt[3*i_parent+2] = lnum;
  }
  PDM_free(shifted_section_distrib);

  _iso->iso_ptp[entity_type] = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **)&iso_gnum,
                                                                         (const int         * )&iso_n_entity,
                                                                                                1,
                                                                         (const int         * )&n_elt2,
                                                                                                1,
                                                                         (const int         **)&iso_parent_triplet_idx,
                                                                         (const int         **) NULL,
                                                                         (const int         **)&iso_parent_trplt,
                                                                                                isos->comm);
  _iso->iso_owner_ptp[entity_type] = PDM_OWNERSHIP_KEEP;

  PDM_free(iso_gnum);
  PDM_free(iso_parent_triplet_idx);
  PDM_free(iso_parent_trplt);
}



static void
_free_iso_entity
(
  PDM_isosurface_t    *isos,
  PDM_mesh_entities_t  entity_type,
  int                  id_iso
)
{
  if (entity_type == PDM_MESH_ENTITY_FACE && isos->entry_mesh_dim < 3) {
    return;
  }

  _isosurface_t *_iso = &isos->isosurfaces[id_iso];

  /* Vertices */
  if (entity_type == PDM_MESH_ENTITY_VTX) {
    /* Partitioned */
    for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {
      if (_iso->iso_owner_vtx_coord[i_part] == PDM_OWNERSHIP_KEEP) {
        PDM_free(_iso->iso_vtx_coord[i_part]);
      }
    }
    PDM_free(_iso->iso_vtx_coord      );

    /* Block-distributed */
    if (_iso->iso_owner_dvtx_coord == PDM_OWNERSHIP_KEEP) {
      PDM_free(_iso->iso_dvtx_coord);
    }
    if (_iso->iso_owner_dparent_wght[PDM_MESH_ENTITY_VTX] == PDM_OWNERSHIP_KEEP) {
      PDM_free(_iso->iso_dentity_parent_wght);
    }
  }

  /* Edges */
  else if (entity_type == PDM_MESH_ENTITY_EDGE) {
    /* Partitioned */
    for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {
      if (_iso->iso_owner_edge_bnd[i_part] == PDM_OWNERSHIP_KEEP) {
        PDM_free(_iso->iso_edge_group_idx [i_part]);
        PDM_free(_iso->iso_edge_group_lnum[i_part]);
        PDM_free(_iso->iso_edge_group_gnum[i_part]);
      }

      if (_iso->iso_owner_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX][i_part] == PDM_OWNERSHIP_KEEP) {
        PDM_free(_iso->iso_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX][i_part]);
      }
    }
    PDM_free(_iso->iso_edge_group_idx );
    PDM_free(_iso->iso_edge_group_lnum);
    PDM_free(_iso->iso_edge_group_gnum);
    PDM_free(_iso->iso_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX]);

    /* Block-distributed */
    if (isos->entry_mesh_dim == 3 &&
        _iso->iso_owner_dedge_bnd == PDM_OWNERSHIP_KEEP) {
      PDM_free(_iso->iso_dedge_group_idx );
      PDM_free(_iso->iso_dedge_group_gnum);
    }

    if (_iso->iso_owner_dconnec[PDM_CONNECTIVITY_TYPE_EDGE_VTX] == PDM_OWNERSHIP_KEEP) {
      PDM_free(_iso->iso_dconnec[PDM_CONNECTIVITY_TYPE_EDGE_VTX]);
    }
  }

  /* Faces */
  else if (entity_type == PDM_MESH_ENTITY_FACE) {
    /* Partitioned */
    for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {
      if (_iso->iso_owner_connec[PDM_CONNECTIVITY_TYPE_FACE_VTX][i_part] == PDM_OWNERSHIP_KEEP) {
        PDM_free(_iso->iso_connec_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX][i_part]);
        PDM_free(_iso->iso_connec    [PDM_CONNECTIVITY_TYPE_FACE_VTX][i_part]);
      }
    }
    PDM_free(_iso->iso_connec_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX]);
    PDM_free(_iso->iso_connec    [PDM_CONNECTIVITY_TYPE_FACE_VTX]);

    /* Block-distributed */
    if (_iso->iso_owner_dconnec[PDM_CONNECTIVITY_TYPE_FACE_VTX] == PDM_OWNERSHIP_KEEP) {
      PDM_free(_iso->iso_dconnec_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX]);
      PDM_free(_iso->iso_dconnec    [PDM_CONNECTIVITY_TYPE_FACE_VTX]);
    }
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Invalid entity_type %d\n", entity_type);
  }


  /* Generic for all entities */

  /*   Partitioned */
  for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {
    if (_iso->iso_owner_gnum[entity_type][i_part] == PDM_OWNERSHIP_KEEP) {
      PDM_free(_iso->iso_entity_gnum[entity_type][i_part]);
    }
    if (isos->entry_is_part == 0) {
      PDM_free(_iso->iso_entity_parent_gnum[entity_type][i_part]);
    }
    if (_iso->iso_owner_parent_idx[entity_type][i_part] == PDM_OWNERSHIP_KEEP) {
      PDM_free(_iso->iso_entity_parent_idx [entity_type][i_part]);
    }
    if (_iso->iso_owner_parent_lnum[entity_type][i_part] == PDM_OWNERSHIP_KEEP) {
      PDM_free(_iso->iso_entity_parent_lnum[entity_type][i_part]);
    }
    if (_iso->iso_owner_parent_wght[entity_type][i_part] == PDM_OWNERSHIP_KEEP) {
      PDM_free(_iso->iso_entity_parent_wght[entity_type][i_part]);
    }
    if (_iso->iso_owner_isovalue_entity_idx[entity_type][i_part] == PDM_OWNERSHIP_KEEP) {
      PDM_free(_iso->isovalue_entity_idx[entity_type][i_part]);
    }
  }
  PDM_free(_iso->iso_n_entity          [entity_type]);
  PDM_free(_iso->iso_entity_gnum       [entity_type]);
  PDM_free(_iso->iso_entity_parent_idx [entity_type]);
  PDM_free(_iso->iso_entity_parent_lnum[entity_type]);
  PDM_free(_iso->iso_entity_parent_gnum[entity_type]);
  PDM_free(_iso->iso_entity_parent_wght[entity_type]);
  PDM_free(_iso->isovalue_entity_idx   [entity_type]);

  if (_iso->  compute_ptp[entity_type] == PDM_TRUE &&
      _iso->iso_owner_ptp[entity_type] == PDM_OWNERSHIP_KEEP) {
    PDM_part_to_part_free(_iso->iso_ptp[entity_type]);
  }


  /*   Block-distributed */
  if (_iso->iso_owner_dparent_idx[entity_type] == PDM_OWNERSHIP_KEEP) {
    PDM_free(_iso->iso_dentity_parent_idx[entity_type]);
  }
  if (_iso->iso_owner_dparent[entity_type] == PDM_OWNERSHIP_KEEP) {
    PDM_free(_iso->iso_dentity_parent_gnum[entity_type]);
  }
  if (_iso->iso_owner_dparent_wght[entity_type] == PDM_OWNERSHIP_KEEP) {
    PDM_free(_iso->iso_dentity_parent_wght[entity_type]);
  }

  if (_iso->iso_owner_disovalue_entity[entity_type] == PDM_OWNERSHIP_KEEP) {
    PDM_free(_iso->disovalue_entity_idx[entity_type]);
    PDM_free(_iso->disovalue_entity    [entity_type]);
  }
}


static void
_free_field
(
  PDM_isosurface_t *isos,
  int               id_iso,
  int               partial
)
{
  _isosurface_t *_iso = &isos->isosurfaces[id_iso];

  if (_iso->field != NULL) {
    for (int i_part=0; i_part<isos->n_part; ++i_part) {
      if (_iso->kind!=PDM_ISO_SURFACE_KIND_FIELD ||
          isos->entry_is_part==0) {
        PDM_free(_iso->field[i_part]);
      }
    }
    if (!partial) {
      PDM_free(_iso->field);
    }
  }
}


static void
_reset_downer
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  _isosurface_t *_iso = &isos->isosurfaces[id_iso];

  _iso->iso_owner_dvtx_coord         = PDM_OWNERSHIP_BAD_VALUE;
  _iso->iso_owner_dedge_bnd          = PDM_OWNERSHIP_BAD_VALUE;

  for (int i_entity=0; i_entity<PDM_MESH_ENTITY_MAX; ++i_entity) {
    _iso->iso_owner_dparent_idx [i_entity] = PDM_OWNERSHIP_BAD_VALUE;
    _iso->iso_owner_dparent     [i_entity] = PDM_OWNERSHIP_BAD_VALUE;
    _iso->iso_owner_dparent_wght[i_entity] = PDM_OWNERSHIP_BAD_VALUE;
  }
  for (int i_connec=0; i_connec<PDM_CONNECTIVITY_TYPE_MAX; ++i_connec) {
    _iso->iso_owner_dconnec[i_connec] = PDM_OWNERSHIP_BAD_VALUE;
  }
}

static void
_free_owner
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  _isosurface_t *_iso = &isos->isosurfaces[id_iso];

  if (_iso->is_computed == PDM_FALSE) {
    return;
  }

  // > Partitionned
  PDM_free(_iso->iso_owner_vtx_coord        );
  PDM_free(_iso->iso_owner_edge_bnd         );
  for (int i_entity=0; i_entity<PDM_MESH_ENTITY_MAX; ++i_entity) {
    PDM_free(_iso->iso_owner_gnum               [i_entity]);
    PDM_free(_iso->iso_owner_parent_lnum        [i_entity]);
    PDM_free(_iso->iso_owner_parent_idx         [i_entity]);
    PDM_free(_iso->iso_owner_parent_wght        [i_entity]);
    PDM_free(_iso->iso_owner_isovalue_entity_idx[i_entity]);
  }
  PDM_free(_iso->iso_owner_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX]);
  PDM_free(_iso->iso_owner_connec[PDM_CONNECTIVITY_TYPE_FACE_VTX]);
}


static void
_free_extracted_parts
(
  PDM_isosurface_t *isos
)
{
  int is_local   = (isos->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL);
  int is_2d_ngon = (isos->entry_mesh_dim == 2) && (!isosurface_is_nodal(isos) || isos->all_simplices==0);

  for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
    if (is_local) {
      PDM_free(isos->extract_vtx_lnum [i_part]);
      if (isos->entry_mesh_dim == 3) {
        PDM_free(isos->extract_cell_lnum[i_part]);
      }
    }

    if (is_local && (is_2d_ngon || (isosurface_is_nodal(isos) && isos->entry_mesh_dim>1))) {
      PDM_free(isos->extract_face_lnum[i_part]);
    }
    if (!is_local && isos->entry_mesh_dim==2 && !isos->all_simplices) {
      PDM_free(isos->extract_face_lnum[i_part]);
    }


    if (is_2d_ngon) {
      PDM_free(isos->extract_tri_vtx [i_part]);
      PDM_free(isos->extract_tri_gnum[i_part]);
      PDM_free(isos->extract_tri_face[i_part]);
    }

    if (isos->extract_tri_tag != NULL) {
      PDM_free(isos->extract_tri_tag[i_part]);
    }
    if (isos->extract_field != NULL) {
      PDM_free(isos->extract_field[i_part]);
    }
  }
  PDM_free(isos->extract_n_vtx    );
  PDM_free(isos->extract_vtx_coord);
  PDM_free(isos->extract_vtx_gnum );
  PDM_free(isos->extract_vtx_lnum );

  PDM_free(isos->extract_n_tri      );
  PDM_free(isos->extract_tri_vtx    );
  PDM_free(isos->extract_tri_gnum   );
  PDM_free(isos->extract_face_lnum  );
  PDM_free(isos->extract_tri_n_group);
  PDM_free(isos->extract_tri_tag    );
  PDM_free(isos->extract_tri_face   );

  PDM_free(isos->extract_n_tetra   );
  PDM_free(isos->extract_tetra_vtx );
  PDM_free(isos->extract_tetra_gnum);
  PDM_free(isos->extract_cell_lnum );

  PDM_free(isos->extract_field);
}



static void
_isosurface_reset
(
 PDM_isosurface_t *isos,
 int               id_isosurface
)
{
  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];

  if (_iso->is_computed == PDM_FALSE) {
    return;
  }

  _free_iso_entity(isos, PDM_MESH_ENTITY_VTX,  id_isosurface);
  _free_iso_entity(isos, PDM_MESH_ENTITY_EDGE, id_isosurface);
  _free_iso_entity(isos, PDM_MESH_ENTITY_FACE, id_isosurface);
  _free_field     (isos, id_isosurface, isos->entry_is_part);
  // > Distributed
  if (isos->entry_is_part==0) {
    _reset_downer(isos, id_isosurface);
  }
  _free_owner(isos, id_isosurface);

  _iso->is_computed = PDM_FALSE;
}


static void
_compute_entity_weights
(
 PDM_isosurface_t   *isos,
 int                 id_iso,
 PDM_mesh_entities_t entity_type
)
{
  _isosurface_t *_iso = &isos->isosurfaces[id_iso];

  // > Compute weight
  double **parent_wght = NULL;
  PDM_malloc(parent_wght, isos->iso_n_part, double *);
  for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {

    int  n_entity   = _iso->iso_n_entity         [entity_type][i_part];
    int *parent_idx = _iso->iso_entity_parent_idx[entity_type][i_part];
    int  n_parent   = parent_idx[n_entity];
    
    PDM_malloc(parent_wght[i_part], n_parent, double);
    
    for (int i_entity=0; i_entity<n_entity; ++i_entity) {
      
      int entity_n_parent = parent_idx[i_entity+1]-parent_idx[i_entity];
      for (int i_parent=parent_idx[i_entity  ]; 
               i_parent<parent_idx[i_entity+1]; ++i_parent) {

        parent_wght[i_part][i_parent] = 1./((double) entity_n_parent);

      }
    }

  } // end loop on partitions

  _iso->iso_owner_dparent_wght[entity_type] = PDM_OWNERSHIP_KEEP;
  _iso->iso_entity_parent_wght[entity_type] = parent_wght;

}


static void
_compute_missing_entity_weights
(
 PDM_isosurface_t *isos,
 int               id_iso
)
{
  _compute_entity_weights(isos, id_iso, PDM_MESH_ENTITY_EDGE);
  if (isos->entry_mesh_dim==3) {
    _compute_entity_weights(isos, id_iso, PDM_MESH_ENTITY_FACE);
  }
}


static void
_isosurface_compute
(
 PDM_isosurface_t *isos,
 int               id_isosurface
)
{
  int i_rank;
  PDM_MPI_Comm_rank(isos->comm, &i_rank);

  int init_entry_mesh_type = isos->entry_mesh_type;

  int dbg = 0;

  if (dbg==1) {
    log_trace("PDM_isosurface:: compute isosurface n°%d\n", id_isosurface);
  }

  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];

  if (_iso->is_computed == PDM_TRUE) {
    // TODO: Warning or Error?
    if (i_rank == 0) {
      printf("WARNING - PDM_isosurface_compute : id_isosurface %d will be reset before being re-computed\n", id_isosurface);
    }
    _isosurface_reset(isos, id_isosurface);
  }

  isosurface_timer_start(isos, ISO_TIMER_TOTAL);

  assert(_iso->is_computed == PDM_FALSE);

  /* Check if edges were provided by the user */
  _do_we_have_edges(isos);

  if (isos->entry_is_part == 0) {
    /* Implicit partitioning */
    isosurface_timer_start(isos, ISO_TIMER_DIST_TO_PART);
    _dist_to_part(isos);
    isosurface_timer_end(isos, ISO_TIMER_DIST_TO_PART);
  }

  /* Evaluate field */
  isosurface_timer_start(isos, ISO_TIMER_COMPUTE_INIT_FIELD);
  _compute_iso_field(isos, id_isosurface, 0);
  isosurface_timer_end(isos, ISO_TIMER_COMPUTE_INIT_FIELD);

  /* Extract elements of interest */
  isosurface_timer_start(isos, ISO_TIMER_EXTRACT);
  _extract(isos, id_isosurface);
  isosurface_timer_end(isos, ISO_TIMER_EXTRACT);

  isosurface_timer_start(isos, ISO_TIMER_NGONIZE);
  if (isosurface_is_nodal(isos)) {
    /* If nodal with sections other than TRIA3 and TETRA4, fall back to ngon */
    _ngonize(isos);
  }
  if (!isosurface_is_nodal(isos) && isos->entry_mesh_dim == 2) {
    /* If ngon 2d, triangulate polygons and use Marching Triangles contouring algorithm */
    _triangulate(isos);
  }
  isosurface_timer_end(isos, ISO_TIMER_NGONIZE);

  /* Evaluate field on extracted mesh */
  isosurface_timer_start(isos, ISO_TIMER_COMPUTE_EXTRACT_FIELD);
  _compute_iso_field(isos, id_isosurface, 1);
  isosurface_timer_end(isos, ISO_TIMER_COMPUTE_EXTRACT_FIELD);

  /* Build isosurface mesh */
  isosurface_timer_start(isos, ISO_TIMER_CONTOURING);
  if (isosurface_is_nodal(isos) || isos->entry_mesh_dim == 2) {
    PDM_isosurface_marching_algo(isos,
                                 id_isosurface);
  }
  else {
    PDM_isosurface_ngon_algo(isos,
                             id_isosurface);
  }
  PDM_part_mesh_free(isos->extract_pmesh);
  isosurface_timer_end(isos, ISO_TIMER_CONTOURING);

  isos->entry_mesh_type = init_entry_mesh_type;

  // > Compute isosurface entity weights
  _compute_missing_entity_weights(isos, id_isosurface);

  if (isos->entry_is_part == 0) {
    /* Block-distribute the isosurface */
    isosurface_timer_start(isos, ISO_TIMER_PART_TO_DIST);
    _part_to_dist(isos, id_isosurface);
    for (int i_entity = PDM_MESH_ENTITY_FACE; i_entity < PDM_MESH_ENTITY_MAX; i_entity++) {
      if (isosurface_is_nodal(isos)) {
        _build_ptp_dist_nodal(isos,
                              id_isosurface,
                              i_entity);
      }
      else {
        _build_ptp_dist(isos,
                        id_isosurface,
                        i_entity);
      }
    }
    isosurface_timer_end(isos, ISO_TIMER_PART_TO_DIST);
  }
  else {
    // Partitioned
    isosurface_timer_start(isos, ISO_TIMER_BUILD_EXCH_PROTOCOL);
    if (isos->extract_kind != PDM_EXTRACT_PART_KIND_LOCAL) {
      for (int i_entity = PDM_MESH_ENTITY_FACE; i_entity < PDM_MESH_ENTITY_MAX; i_entity++) {
        _build_ptp_part(isos,
                        id_isosurface,
                        i_entity);
      }
    }
    isosurface_timer_end(isos, ISO_TIMER_BUILD_EXCH_PROTOCOL);
  }

  // > Free mesh extraction arrays
  _free_extracted_parts(isos);
  PDM_extract_part_free(isos->extrp);

  _iso->is_computed = PDM_TRUE;

  isosurface_timer_end(isos, ISO_TIMER_TOTAL);
}


/*=============================================================================
 * Public function prototypes
 *============================================================================*/


PDM_isosurface_t *
PDM_isosurface_create
(
 PDM_MPI_Comm             comm,
 int                      mesh_dimension
)
{
  PDM_isosurface_t *isos = NULL;
  PDM_malloc(isos, 1, PDM_isosurface_t);
  memset(isos, 0, sizeof(PDM_isosurface_t));

  // > Save communicator
  isos->comm = comm;

  // > Init tolerance
  isos->ISOSURFACE_EPS = 0.;

  // > Entry mesh information
  isos->entry_is_part   = -1;
  isos->entry_mesh_type = 0;
  isos->entry_mesh_dim  = mesh_dimension;

  // > Isosurface mesh information
  // isos->iso_elt_type = elt_type;
  isos->extract_kind = PDM_EXTRACT_PART_KIND_LOCAL;
  isos->n_part = 1;

  isos->we_have_edges = -1;

  isos->all_simplices = 0;


  // > Reset timers
  memset(isos->times_current, 0, ISO_TIMER_N_STEPS*sizeof(double));
  memset(isos->times_cumul,   0, ISO_TIMER_N_STEPS*sizeof(double));

  return isos;
}


void
PDM_isosurface_n_group_set
(
  PDM_isosurface_t    *isos,
  PDM_mesh_entities_t  entity_type,
  int                  n_group
)
{
  switch (entity_type) {
    case PDM_MESH_ENTITY_FACE:
      isos->n_group_face = n_group;
      break;
    // case PDM_MESH_ENTITY_EDGE:
    //   isos->n_group_edge = n_group;
    //   break;
    // case PDM_MESH_ENTITY_VTX:
    //   isos->n_group_vtx  = n_group;
    //   break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "invalid entity_type (%d) for isosurface boundary.\n", entity_type);
      break;
  }
}


void
PDM_isosurface_set_tolerance
(
  PDM_isosurface_t *isos,
  double            tolerance
)
{
  if (tolerance < 0.) {
    PDM_error(__FILE__, __LINE__, 0, "Invalid tolerance : %f (must be >= 0)\n", tolerance);
  }

  isos->ISOSURFACE_EPS = tolerance;

  for (int id_isosurface = 0; id_isosurface < isos->n_isosurface; id_isosurface++) {
    PDM_ISOSURFACE_CHECK_ISOVALUES_NOT_TOO_CLOSE(isos, id_isosurface);
  }
}


void
PDM_isosurface_set_isovalues
(
  PDM_isosurface_t *isos,
  int               id_isosurface,
  int               n_isovalues,
  double           *isovalues
)
{
  PDM_ISOSURFACE_CHECK_ID(isos, id_isosurface);

  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];

  _iso->n_isovalues = n_isovalues;
  PDM_realloc(_iso->isovalues, _iso->isovalues, n_isovalues, double);
  for (int i=0; i<n_isovalues; ++i) {
    _iso->isovalues[i] = isovalues[i];
  }

  PDM_ISOSURFACE_CHECK_ISOVALUES_NOT_TOO_CLOSE(isos, id_isosurface);
}


int
PDM_isosurface_add
(
  PDM_isosurface_t       *isos,
  PDM_iso_surface_kind_t  kind,
  int                     n_isovalues,
  double                 *isovalues
)
{
  int id_isosurface = isos->n_isosurface;
  isos->n_isosurface++;

  PDM_realloc(isos->isosurfaces, isos->isosurfaces, isos->n_isosurface, _isosurface_t);

  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];

  memset(_iso, 0, sizeof(_isosurface_t));

  _iso->kind        = kind;
  _iso->is_computed = PDM_FALSE;

  PDM_isosurface_set_isovalues(isos,
                               id_isosurface,
                               n_isovalues,
                               isovalues);

  for (int i = 0; i < PDM_MESH_ENTITY_MAX; i++) {
    _iso->compute_ptp  [i] = PDM_FALSE;
    _iso->iso_owner_ptp[i] = PDM_OWNERSHIP_BAD_VALUE;

    _iso->iso_owner_dparent_idx[i] = PDM_OWNERSHIP_BAD_VALUE;
    _iso->iso_owner_dparent    [i] = PDM_OWNERSHIP_BAD_VALUE;
  }

  for (int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; i++) {
    _iso->iso_owner_dconnec[i] = PDM_OWNERSHIP_BAD_VALUE;
  }

  return id_isosurface;
}


void
PDM_isosurface_equation_set
(
 PDM_isosurface_t *isos,
 int               id_isosurface,
 double           *coeff
)
{
  PDM_ISOSURFACE_CHECK_ID(isos, id_isosurface);

  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];
  int n_coeff = _iso_surface_kind_n_coeff(_iso->kind);

  if (_iso->eq_coeffs == NULL) {
    PDM_malloc(_iso->eq_coeffs, n_coeff, double);
  }

  for (int i_coeff=0; i_coeff<n_coeff; ++i_coeff) {
    _iso->eq_coeffs[i_coeff] = coeff[i_coeff];
  }
}


void
PDM_isosurface_field_function_set
(
 PDM_isosurface_t                *isos,
 int                              id_isosurface,
 PDM_isosurface_field_function_t  func
)
{
  PDM_ISOSURFACE_CHECK_ID(isos, id_isosurface);
  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];

  if (_iso->kind == PDM_ISO_SURFACE_KIND_FUNCTION) {
    _iso->field_function = func;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Isosurface n°%d doesn't support PDM_isosurface_field_function_set method cause its kind is %d.\n", id_isosurface, _iso->kind);
  }

}


void
isosurface_field_function_set_python
(
  PDM_isosurface_t                       *isos,
  int                                     id_isosurface,
  PDM_isosurface_field_function_python_t  func
)
{
  PDM_ISOSURFACE_CHECK_ID(isos, id_isosurface);
  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];

  if (_iso->kind == PDM_ISO_SURFACE_KIND_FUNCTION) {
    _iso->field_function_python = func;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Isosurface n°%d doesn't support PDM_isosurface_field_function_set method cause its kind is %d.\n", id_isosurface, _iso->kind);
  }
}

void
isosurface_python_object_set
(
  PDM_isosurface_t *isos,
  void             *python_object
)
{
  isos->python_object = python_object;
}


void
PDM_isosurface_reset
(
 PDM_isosurface_t *isos,
 int               id_isosurface
)
{
  PDM_ISOSURFACE_CHECK_ID(isos, id_isosurface);

  if (id_isosurface < 0) { // Reset all isosurfaces
    for (int i = 0; i < isos->n_isosurface; i++) {
      _isosurface_reset(isos, i);
    }
  }
  else {
    _isosurface_reset(isos, id_isosurface);
  }
}


void
PDM_isosurface_compute
(
 PDM_isosurface_t *isos,
 int               id_isosurface
)
{
  if (isos->iso_n_part == 0) {
    if (isos->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) {
      isos->iso_n_part = isos->n_part;
    }
    else { // REEQUILIBRATE
      isos->iso_n_part = 1;
    }
  }
  else {
    if (isos->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL &&
        isos->iso_n_part   != isos->n_part) {
      PDM_error(__FILE__, __LINE__, 0,
                "In PDM_EXTRACT_PART_KIND_LOCAL mode n_part_out (%d) must be equal to n_part_in (%d)\n",
                isos->iso_n_part, isos->n_part);
    }
  }

  if (id_isosurface < 0) { // Compute all isosurfaces
    for (int i = 0; i < isos->n_isosurface; i++) {
      _isosurface_compute(isos, i);
    }
  }
  else {
    PDM_ISOSURFACE_CHECK_ID(isos, id_isosurface);
    _isosurface_compute(isos, id_isosurface);
  }
}


void
PDM_isosurface_dump_times
(
 PDM_isosurface_t *isos
)
{
  /* TODO: choose between current and cumulative elapsed times? */
  int i_rank;
  PDM_MPI_Comm_rank(isos->comm, &i_rank);

  for (_isosurface_timer_step_t step = 0; step < ISO_TIMER_N_STEPS; step++) {

    /* Skip irrelevant steps */
    if (isos->entry_is_part != 0) {
      // Partitioned
      if (step == ISO_TIMER_DIST_TO_PART ||
          step == ISO_TIMER_PART_TO_DIST) {
        continue;
      }
    }

    if (!isosurface_is_nodal(isos)) {
      // Ngon
      if (step == ISO_TIMER_NGONIZE) {
        continue;
      }
    }

    /* Compute min/max and print */
    _dump_delta_time(isos->comm,
                     isos->times_cumul[step],
                     _isosurface_timer_step_name[step]);
  }

}


void
PDM_isosurface_free
(
  PDM_isosurface_t  *isos
)
{
  if (isos->entry_is_part == 0) {
    // Block-distributed
    PDM_block_to_part_free(isos->btp_vtx);

    if (isosurface_is_nodal(isos)) {
      PDM_part_mesh_nodal_free(isos->pmesh_nodal);
    }
    else {
      PDM_free(isos->cell_face    [0]);
      PDM_free(isos->cell_face_idx[0]);
      PDM_free(isos->face_edge    [0]);
      PDM_free(isos->face_edge_idx[0]);
      PDM_free(isos->face_vtx     [0]);
      PDM_free(isos->face_vtx_idx [0]);
      PDM_free(isos->edge_vtx     [0]);
      PDM_free(isos->vtx_coord    [0]);
      PDM_free(isos->cell_gnum    [0]);
      PDM_free(isos->face_gnum    [0]);
      PDM_free(isos->edge_gnum    [0]);
      PDM_free(isos->vtx_gnum     [0]);
      if (isos->entry_mesh_dim==3) {
        PDM_free(isos->group_face_idx [0]);
        PDM_free(isos->group_face     [0]);
        PDM_free(isos->group_face_gnum[0]);
      }
    }
  }
  PDM_free(isos->n_cell);
  PDM_free(isos->n_face);
  PDM_free(isos->n_edge);
  PDM_free(isos->n_vtx);
  PDM_free(isos->cell_face);
  PDM_free(isos->cell_face_idx);
  PDM_free(isos->face_edge);
  PDM_free(isos->face_edge_idx);
  PDM_free(isos->face_vtx);
  PDM_free(isos->face_vtx_idx);
  PDM_free(isos->edge_vtx);
  PDM_free(isos->vtx_coord);
  PDM_free(isos->cell_gnum);
  PDM_free(isos->face_gnum);
  PDM_free(isos->edge_gnum);
  PDM_free(isos->vtx_gnum);
  PDM_free(isos->group_face_idx);
  PDM_free(isos->group_face);
  PDM_free(isos->group_face_gnum);

  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    _isosurface_t *_iso = &isos->isosurfaces[id_iso];

    PDM_free(_iso->isovalues);
    PDM_free(_iso->eq_coeffs);

    _free_field(isos, id_iso, 0);

    if (_iso->is_computed == PDM_TRUE) {
      _free_iso_entity(isos, PDM_MESH_ENTITY_VTX , id_iso);
      _free_iso_entity(isos, PDM_MESH_ENTITY_EDGE, id_iso);
      _free_iso_entity(isos, PDM_MESH_ENTITY_FACE, id_iso);
    }

    _free_owner(isos, id_iso);
  }
  PDM_free(isos->isosurfaces);

  PDM_free(isos);
}


void
PDM_isosurface_enable_part_to_part
(
  PDM_isosurface_t     *isos,
  int                   id_isosurface,
  PDM_mesh_entities_t   entity_type,
  int                   unify_parent_info
)
{
  /**
   * TODO:
   *  - unify_parent_info : allow demanding user to ask to get all parent over all procs (not urgent)
   *      - build additional ptp to get info
   */
  PDM_ISOSURFACE_CHECK_ID(isos, id_isosurface);

  if (unify_parent_info!=0) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: unify_parent_info option not implemented yet.\n");
  }

  if (id_isosurface >= isos->n_isosurface) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_enable_part_to_part : Invalid id_isosurface %d (n_isosurface = %d)\n", id_isosurface, isos->n_isosurface);
  }

  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];
  _iso->compute_ptp[entity_type] = PDM_TRUE;
}


void
PDM_isosurface_part_to_part_get
(
 PDM_isosurface_t     *isos,
 int                   id_isosurface,
 PDM_mesh_entities_t   entity_type,
 PDM_part_to_part_t  **ptp,
 PDM_ownership_t       ownership
)
{
  PDM_ISOSURFACE_CHECK_ID      (isos, id_isosurface);
  PDM_ISOSURFACE_CHECK_COMPUTED(isos, id_isosurface);
  
  PDM_ISOSURFACE_CHECK_ENTITY_TYPE(entity_type);

  _isosurface_t *_iso = &isos->isosurfaces[id_isosurface];

  if (_iso->compute_ptp[entity_type] == PDM_FALSE) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: part_to_part for entity %d of isosurface %d is not computed.\n", entity_type, id_isosurface);
  }

  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    _iso->iso_owner_ptp[entity_type] = ownership;
  }

  *ptp = _iso->iso_ptp[entity_type];
}


#ifdef  __cplusplus
}
#endif
