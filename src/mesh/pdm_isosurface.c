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
#include "pdm_mesh_nodal.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_extract_part.h"
#include "pdm_vtk.h"
#include "pdm_part_mesh_nodal_to_pmesh.h"

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

static inline int
_is_nodal
(
 PDM_isosurface_t *isos
)
{
  return PDM_ABS(isos->entry_mesh_type) == 3;
}

static void
_do_we_have_edges
(
 PDM_isosurface_t *isos
)
{
  if (isos->we_have_edges >= 0) {
    return;
  }

  if (_is_nodal(isos)) {
    return;
  }

  int i_have_edges    = 0;
  int i_have_face_vtx = 0;

  if (isos->is_dist_or_part == 0) {
    // Block-distributed
    i_have_edges    = (isos->distrib_edge != NULL) && (isos->dface_edge_idx != NULL) && (isos->dface_edge != NULL) && (isos->dedge_vtx != NULL);
    i_have_face_vtx = (isos->dface_vtx_idx != NULL) && (isos->dface_vtx != NULL);
  }
  else {
    // Partitioned
    for (int i_part = 0; i_part < isos->n_part; i_part++) {

      if (isos->n_face[i_part] == 0) {
        continue;
      }

      if (isos->n_edge[i_part] > 0) {
        if (isos->face_edge_idx[i_part] != NULL &&
            isos->face_edge    [i_part] != NULL &&
            isos->edge_vtx     [i_part] != NULL) {
          i_have_edges = 1;
        }
      }

      if (isos->face_vtx_idx[i_part] != NULL &&
          isos->face_vtx    [i_part] != NULL) {
        i_have_edges = 1;
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
    return;
  }

  isos->dist_to_part_computed = 1;

  if (isos->is_dist_or_part != 0) {
    PDM_error(__FILE__, __LINE__, 0, "Expected block-distributed but got partitioned\n");
  }

  int i_rank;
  PDM_MPI_Comm_rank(isos->comm, &i_rank);

  isos->n_part = 1;

  if (_is_nodal(isos)) {
    // Nodal
    //TODO dmesh_nodal to part_mesh_nodal
    PDM_error(__FILE__, __LINE__, 0, "_dist_to_part Nodal not yet implemented\n");
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
      free(dedge_vtx_idx);
      free(edge_vtx_idx );
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
    const PDM_g_num_t *pface_ln_to_gn[1] = {face_ln_to_gn};
    PDM_g_num_t **group_face_ln_to_gn = NULL;
    PDM_part_distgroup_to_partgroup(isos->comm,
                                    isos->distrib_face,
                                    isos->n_dgroup_face,
                                    isos->dgroup_face_idx,
                                    isos->dgroup_face,
                                    isos->n_part,
                                    &n_face,
                                    pface_ln_to_gn,
                                    &isos->group_face_idx,
                                    &isos->group_face,
                                    &group_face_ln_to_gn);
    free(group_face_ln_to_gn[0]);
    free(group_face_ln_to_gn);

    /* Store in struct */
    isos->btp_vtx = btp_vtx; // useful to keep for transferring discrete fields

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
    PDM_malloc(isos->n_group_face , 1, int          );

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
    isos->n_group_face [0] = isos->n_dgroup_face;
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
  if (isos->is_dist_or_part == 0) {
    // Block-distributed
    if (isos->field[id_isosurface] == NULL) {
      PDM_malloc(isos->field[id_isosurface], isos->n_part, double *);
    }
    assert(isos->dist_to_part_computed);
  }

  if (isos->kind[id_isosurface] == PDM_ISO_SURFACE_KIND_FIELD) {
    if (use_extract) {
      if (isos->extract_kind == PDM_EXTRACT_PART_KIND_REEQUILIBRATE) {
        PDM_part_to_part_t *ptp = NULL;
        PDM_extract_part_part_to_part_get(isos->extrp[id_isosurface],
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
                      (const void  **) isos->field[id_isosurface],
                                       NULL,
                            (void ***) &isos->extract_field[id_isosurface],
                                       &request);

        PDM_part_to_part_reverse_iexch_wait(ptp, request);
      }
      else if (isos->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) {
        PDM_malloc(isos->extract_field[id_isosurface], isos->n_part, double *);

        for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {
          int *parent = NULL;
          int n_vtx = PDM_extract_part_parent_lnum_get(isos->extrp[id_isosurface],
                                                       i_part,
                                                       PDM_MESH_ENTITY_VTX,
                                                       &parent,
                                                       PDM_OWNERSHIP_KEEP);

          PDM_malloc(isos->extract_field[id_isosurface][i_part], n_vtx, double);
          for (int i = 0; i < n_vtx; i++) {
            isos->extract_field[id_isosurface][i_part][i] = isos->field[id_isosurface][i_part][parent[i] - 1];
          }
        }
      }
      else {
        PDM_error(__FILE__, __LINE__, 0, "Invalid extract_kind %d\n", (int) isos->extract_kind);
      }

    }
    else {
      if (isos->is_dist_or_part == 0) {
        // Transfer discrete field from block to part
        assert(isos->btp_vtx != NULL);
        assert(isos->dfield[id_isosurface] != NULL);

        int one = 1;
        PDM_block_to_part_exch(isos->btp_vtx,
                               sizeof(double),
                               PDM_STRIDE_CST_INTERLACED,
                               &one,
                    (void   *) isos->dfield[id_isosurface],
                               NULL,
                    (void ***) &isos->field[id_isosurface]);
      }
    }

    return;
  }

  int     *n_vtx     = NULL;
  double **vtx_coord = NULL;
  PDM_malloc(n_vtx    , isos->n_part, int     );
  PDM_malloc(vtx_coord, isos->n_part, double *);
  double **field     = NULL;
  if (use_extract) {
    if (isos->extract_field[id_isosurface] == NULL) {
      PDM_malloc(isos->extract_field[id_isosurface], isos->n_part, double *);
    }
    field = isos->extract_field[id_isosurface];
    for (int i_part = 0; i_part < isos->n_part; i_part++) {
      n_vtx[i_part] = PDM_extract_part_vtx_coord_get(isos->extrp[id_isosurface],
                                                     i_part,
                                                     &vtx_coord[i_part],
                                                     PDM_OWNERSHIP_KEEP);
    }
  }
  else {
    field = isos->field[id_isosurface];
    if (_is_nodal(isos)) {
      for (int i_part = 0; i_part < isos->n_part; i_part++) {
        n_vtx    [i_part] = PDM_part_mesh_nodal_n_vtx_get    (isos->pmesh_nodal, i_part);
        vtx_coord[i_part] = PDM_part_mesh_nodal_vtx_coord_get(isos->pmesh_nodal, i_part);
      }
    }
    else {
      for (int i_part = 0; i_part < isos->n_part; i_part++) {
        n_vtx    [i_part] = isos->n_vtx    [i_part];
        vtx_coord[i_part] = isos->vtx_coord[i_part];
      }
    }
  }

  for (int i_part = 0; i_part < isos->n_part; i_part++) {
    PDM_malloc(field[i_part], n_vtx[i_part], double);
  }

  /* Fill */
  if (isos->kind[id_isosurface] == PDM_ISO_SURFACE_KIND_FUNCTION) {
    for (int i_part = 0; i_part < isos->n_part; i_part++) {
      for (int i_vtx = 0; i_vtx < n_vtx[i_part]; i_vtx++) {
        isos->field_function[id_isosurface](vtx_coord[i_part][3*i_vtx  ],
                                            vtx_coord[i_part][3*i_vtx+1],
                                            vtx_coord[i_part][3*i_vtx+2],
                                            &field[i_part][i_vtx]);
      }
    }
  }
  else {
    double (*field_function) (const double, const double, const double, double *);

    switch (isos->kind[id_isosurface]) {
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
        PDM_error(__FILE__, __LINE__, 0, "Invalid isosurface type %d for id_isosurface %d.\n", isos->kind[id_isosurface], id_isosurface);
      }
    }

    for (int i_part = 0; i_part < isos->n_part; i_part++) {
      for (int i_vtx = 0; i_vtx < n_vtx[i_part]; i_vtx++) {
        field[i_part][i_vtx] = field_function(vtx_coord[i_part][3*i_vtx  ],
                                              vtx_coord[i_part][3*i_vtx+1],
                                              vtx_coord[i_part][3*i_vtx+2],
                                              isos->eq_coeffs[id_isosurface]);
      }
    }
  }

  free(n_vtx    );
  free(vtx_coord);
}



// --->> migrer dans priv?
static const double ISOSURFACE_EPS = 1e-6;

static inline int
_sign
(
  const double v
)
{
  // if (v < -ISOSURFACE_EPS) {
  //   return -1;
  // }
  // else if (v > ISOSURFACE_EPS) {
  //   return 1;
  // }
  // else {
  //   return 0;
  // }
  return (v > ISOSURFACE_EPS);
}


static inline int
_cross_0_level_ngon
(
  const double v0,
  const double v1
)
{
  return _sign(v0) != _sign(v1);
}


static inline int
_cross_any_level_ngon
(
  const double v0,
  const double v1,
  const int    n_isovalues,
  const double isovalues[]
)
{
  int n_crossings = 0;
  for (int i = 0; i < n_isovalues; i++) {
    n_crossings += _cross_0_level_ngon(v0 - isovalues[i], v1 - isovalues[i]);
  }

  return n_crossings;
}
// <<---

static void
_extract
(
 PDM_isosurface_t *isos,
 int               id_isosurface
)
{
  if (isos->is_dist_or_part == 0) {
    // Block-distributed
    assert(isos->dist_to_part_computed);
    isos->extract_kind = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    isos->part_method  = PDM_SPLIT_DUAL_WITH_HILBERT;
  }

  PDM_extract_part_t *extrp = PDM_extract_part_create(isos->entry_mesh_dim,
                                                      isos->n_part,
                                                      isos->n_part,
                                                      isos->extract_kind,
                                                      isos->part_method,
                                                      PDM_FALSE,
                                                      PDM_OWNERSHIP_KEEP,
                                                      isos->comm);
  isos->extrp[id_isosurface] = extrp;


  int  *n_extract    = PDM_array_zeros_int(isos->n_part);
  int **extract_lnum = NULL;
  PDM_malloc(extract_lnum, isos->n_part, int *);

  int          *pn_cell        = NULL;
  int          *pn_face        = NULL;
  int          *pn_edge        = NULL;
  int          *pn_vtx         = NULL;
  int         **pcell_face_idx = NULL;
  int         **pcell_face     = NULL;
  int         **pface_edge_idx = NULL;
  int         **pface_edge     = NULL;
  int         **pface_vtx_idx  = NULL;
  int         **pface_vtx      = NULL;
  int         **pedge_vtx      = NULL;
  double      **pvtx_coord     = NULL;
  PDM_g_num_t **pcell_ln_to_gn = NULL;
  PDM_g_num_t **pface_ln_to_gn = NULL;
  PDM_g_num_t **pedge_ln_to_gn = NULL;
  PDM_g_num_t **pvtx_ln_to_gn  = NULL;

  if (_is_nodal(isos)) {
    // Nodal
    PDM_malloc(pn_cell       , isos->n_part, int          );
    PDM_malloc(pn_face       , isos->n_part, int          );
    PDM_malloc(pn_edge       , isos->n_part, int          );
    PDM_malloc(pn_vtx        , isos->n_part, int          );
    PDM_malloc(pcell_ln_to_gn, isos->n_part, PDM_g_num_t *);
    PDM_malloc(pface_ln_to_gn, isos->n_part, PDM_g_num_t *);
    PDM_malloc(pvtx_ln_to_gn , isos->n_part, PDM_g_num_t *);

    if (isos->entry_mesh_dim == 2) {
      for (int i_part = 0; i_part < isos->n_part; i_part++) {
        pn_cell[i_part] = 0;
      }
    }

    // TODO...
    // PDM_extract_part_part_nodal_set...
    PDM_error(__FILE__, __LINE__, 0, "_extract Nodal not implemented yet\n");
  }

  else {
    // Ngon

    pn_cell        = isos->n_cell;
    pn_face        = isos->n_face;
    pn_edge        = isos->n_edge;
    pn_vtx         = isos->n_vtx;
    pcell_face_idx = isos->cell_face_idx;
    pcell_face     = isos->cell_face;
    pface_edge_idx = isos->face_edge_idx;
    pface_edge     = isos->face_edge;
    pface_vtx_idx  = isos->face_vtx_idx;
    pface_vtx      = isos->face_vtx;
    pedge_vtx      = isos->edge_vtx;
    pvtx_coord     = isos->vtx_coord;
    pcell_ln_to_gn = isos->cell_gnum;
    pface_ln_to_gn = isos->face_gnum;
    pedge_ln_to_gn = isos->edge_gnum;
    pvtx_ln_to_gn  = isos->vtx_gnum;


    if (isos->entry_mesh_dim == 2) {
      // 2D
      for (int i_part = 0; i_part < isos->n_part; i_part++) {
        PDM_malloc(extract_lnum[i_part], pn_face[i_part], int);
      }

      if (isos->we_have_edges) {
        for (int i_part = 0; i_part < isos->n_part; i_part++) {
          for (int i_face = 0; i_face < pn_face[i_part]; i_face++) {
            int is_selected = 0;
            for (int idx_edge = pface_edge_idx[i_part][i_face]; idx_edge < pface_edge_idx[i_part][i_face+1]; idx_edge++) {
              int i_edge = PDM_ABS(pface_edge[i_part][idx_edge]) - 1;
              int i_vtx0 = pedge_vtx[i_part][2*i_edge  ] - 1;
              int i_vtx1 = pedge_vtx[i_part][2*i_edge+1] - 1;
              double val0 = isos->field[id_isosurface][i_part][i_vtx0];
              double val1 = isos->field[id_isosurface][i_part][i_vtx1];
              if (_cross_any_level_ngon(val0, val1, isos->n_isovalues[id_isosurface], isos->isovalues[id_isosurface])) {
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
        for (int i_part = 0; i_part < isos->n_part; i_part++) {
          for (int i_face = 0; i_face < pn_face[i_part]; i_face++) {
            int is_selected = 0;
            int face_vtx_n = pface_vtx_idx[i_part][i_face+1] - pface_vtx_idx[i_part][i_face];
            int *fv = pface_vtx[i_part] + pface_vtx_idx[i_part][i_face];
            for (int i = 0; i < face_vtx_n; i++) {
              int i_vtx0 = fv[ i              ] - 1;
              int i_vtx1 = fv[(i+1)%face_vtx_n] - 1;
              double val0 = isos->field[id_isosurface][i_part][i_vtx0];
              double val1 = isos->field[id_isosurface][i_part][i_vtx1];
              if (_cross_any_level_ngon(val0, val1, isos->n_isovalues[id_isosurface], isos->isovalues[id_isosurface])) {
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
        PDM_malloc(extract_lnum[i_part], pn_cell[i_part], int);
      }
      if (isos->we_have_edges == 0) {
        for (int i_part = 0; i_part < isos->n_part; i_part++) {
          for (int i_cell = 0; i_cell < pn_cell[i_part]; i_cell++) {
            int is_selected = 0;
            for (int idx_face = pcell_face_idx[i_part][i_cell]; idx_face < pcell_face_idx[i_part][i_cell+1]; idx_face++) {
              int i_face = PDM_ABS(pcell_face[i_part][idx_face]) - 1;
              int face_vtx_n = pface_vtx_idx[i_part][i_face+1] - pface_vtx_idx[i_part][i_face];
              int *fv = pface_vtx[i_part] + pface_vtx_idx[i_part][i_face];
              for (int i = 0; i < face_vtx_n; i++) {
                int i_vtx0 = fv[ i              ] - 1;
                int i_vtx1 = fv[(i+1)%face_vtx_n] - 1;
                double val0 = isos->field[id_isosurface][i_part][i_vtx0];
                double val1 = isos->field[id_isosurface][i_part][i_vtx1];
                if (_cross_any_level_ngon(val0, val1, isos->n_isovalues[id_isosurface], isos->isovalues[id_isosurface])) {
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
        for (int i_part = 0; i_part < isos->n_part; i_part++) {
          for (int i_cell = 0; i_cell < pn_cell[i_part]; i_cell++) {
            int is_selected = 0;
            for (int idx_face = pcell_face_idx[i_part][i_cell]; idx_face < pcell_face_idx[i_part][i_cell+1]; idx_face++) {
              int i_face = PDM_ABS(pcell_face[i_part][idx_face]) - 1;
              for (int idx_edge = pface_edge_idx[i_part][i_face]; idx_edge < pface_edge_idx[i_part][i_face+1]; idx_edge++) {
                int i_edge = PDM_ABS(pface_edge[i_part][idx_edge]) - 1;
                int i_vtx0 = pedge_vtx[i_part][2*i_edge  ] - 1;
                int i_vtx1 = pedge_vtx[i_part][2*i_edge+1] - 1;
                double val0 = isos->field[id_isosurface][i_part][i_vtx0];
                double val1 = isos->field[id_isosurface][i_part][i_vtx1];
                if (_cross_any_level_ngon(val0, val1, isos->n_isovalues[id_isosurface], isos->isovalues[id_isosurface])) {
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
      extract_lnum[i_part] = realloc(extract_lnum[i_part], sizeof(int) * n_extract[i_part]);
    }

  } // End Ngon

  PDM_extract_part_n_group_set(extrp, PDM_BOUND_TYPE_FACE, isos->n_group_face[0]);

  for (int i_part = 0; i_part < isos->n_part; i_part++) {
    PDM_extract_part_part_set(extrp,
                              i_part,
                              pn_cell       [i_part],
                              pn_face       [i_part],
                              pn_edge       [i_part],
                              pn_vtx        [i_part],
                              pcell_face_idx[i_part],
                              pcell_face    [i_part],
                              pface_edge_idx[i_part],
                              pface_edge    [i_part],
                              pedge_vtx     [i_part],
                              pface_vtx_idx [i_part],
                              pface_vtx     [i_part],
                              pcell_ln_to_gn[i_part],
                              pface_ln_to_gn[i_part],
                              pedge_ln_to_gn[i_part],
                              pvtx_ln_to_gn [i_part],
                              pvtx_coord    [i_part]);

    PDM_extract_part_selected_lnum_set(extrp,
                                       i_part,
                                       n_extract   [i_part],
                                       extract_lnum[i_part]);

    for (int i_group = 0; i_group < isos->n_group_face[0]; i_group++) {
      PDM_extract_part_part_group_set(extrp,
                                      i_part,
                                      i_group,
                                      PDM_BOUND_TYPE_FACE,
                                      isos->group_face_idx[i_part][i_group+1] - isos->group_face_idx[i_part][i_group],
                                      isos->group_face[i_part]      + isos->group_face_idx[i_part][i_group],
                                      isos->group_face_gnum[i_part] + isos->group_face_idx[i_part][i_group]);
    }
  }


  PDM_extract_part_compute(extrp);


  if (_is_nodal(isos)) {
    // for (int i_part = 0; i_part < isos->n_part; i_part++) {
    //   free(pn_cell       [i_part]);
    //   free(pn_face       [i_part]);
    //   free(pn_edge       [i_part]);
    //   free(pn_vtx        [i_part]);
    //   free(pcell_ln_to_gn[i_part]);
    //   free(pface_ln_to_gn[i_part]);
    //   free(pedge_ln_to_gn[i_part]);
    //   free(pvtx_ln_to_gn [i_part]);
    // }
    // free(pn_cell       );
    // free(pn_face       );
    // free(pn_edge       );
    // free(pn_vtx        );
    // free(pcell_ln_to_gn);
    // free(pface_ln_to_gn);
    // free(pedge_ln_to_gn);
    // free(pvtx_ln_to_gn );
  }

  for (int i_part = 0; i_part < isos->n_part; i_part++) {
    free(extract_lnum[i_part]);
  }
  free(n_extract   );
  free(extract_lnum);
}


/**
 * \brief Convert nodal multi-sections to ngon
 */
static void
_ngonize
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  if (!_is_nodal(isos)) {
    // Already ngon
    return;
  }

  assert(isos->extrp[id_iso] != NULL);

  PDM_part_mesh_nodal_elmts_t *extract_pmne = NULL;
  PDM_extract_part_part_mesh_nodal_get(isos->extrp[id_iso],
                                       &extract_pmne,
                                       PDM_OWNERSHIP_KEEP);

  PDM_part_mesh_nodal_t *extract_pmn = PDM_part_mesh_nodal_create(isos->entry_mesh_dim,
                                                                  isos->n_part,
                                                                  isos->comm);


  for (int i_part = 0; i_part < isos->n_part; i_part++) {
    int extract_n_vtx = PDM_extract_part_n_entity_get(isos->extrp[id_iso],
                                                      i_part,
                                                      PDM_MESH_ENTITY_VTX);
    double *extract_vtx_coord = NULL;
    PDM_extract_part_vtx_coord_get(isos->extrp[id_iso],
                                   i_part,
                                   &extract_vtx_coord,
                                   PDM_OWNERSHIP_KEEP);

    PDM_g_num_t *extract_vtx_ln_to_gn = NULL;
    PDM_extract_part_parent_ln_to_gn_get(isos->extrp[id_iso],
                                         i_part,
                                         PDM_MESH_ENTITY_VTX,
                                         &extract_vtx_ln_to_gn,
                                         PDM_OWNERSHIP_KEEP);

    PDM_part_mesh_nodal_coord_set(extract_pmn,
                                  i_part,
                                  extract_n_vtx,
                                  extract_vtx_coord,
                                  extract_vtx_ln_to_gn,
                                  PDM_OWNERSHIP_USER);
  }

  PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(extract_pmn,
                                                extract_pmne);

  isos->extract_pmesh_nodal = extract_pmn;

  /* Inspect nodal sections and check wether we only have simplices */
  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (extract_pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(extract_pmne);

  int all_simplices = 1;
  for (int i_section = 0; i_section < n_section; i_section++) {
    PDM_Mesh_nodal_elt_t elt_type = PDM_part_mesh_nodal_elmts_section_type_get(extract_pmne,
                                                                               sections_id[i_section]);

    if (elt_type != PDM_MESH_NODAL_TRIA3 &&
        elt_type != PDM_MESH_NODAL_TETRA4) {
      all_simplices = 0;
      break;
    }
  }

  // We assume all ranks have the same sections, so no need for Allreduce

  if (all_simplices) {
    // TODO: decompose into tria_vtx, tetra_vtx, ...
  }
  else {
    // We have elements other than simplices, we need to ngonize
    isos->we_have_edges   = 1; // sure about this hack?
    isos->entry_mesh_type = 1 * PDM_SIGN(isos->entry_mesh_type); // we are in fact ngon from now on
    PDM_part_mesh_t *pmesh = PDM_part_mesh_nodal_to_part_mesh(extract_pmn,
                                                              PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
                                                              PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);
    // TODO: unpack pmesh either into isos->extrp[id_iso] or into isos directly
    PDM_error(__FILE__, __LINE__, 0, "Work left to do\n");
  }

}


/**
 * \brief Block-distrubute isosurface mesh
 */
static void
_part_to_dist
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  PDM_UNUSED(isos);
  PDM_UNUSED(id_iso);

  // Block to parts!
  // TODO...
  PDM_error(__FILE__, __LINE__, 0, "_part_to_dist not yet implemented, but everything OK so far :D\n");
}


static void
_free_iso_vtx
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  if (isos->iso_owner_vtx_coord[id_iso]!=NULL) {
    for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
      if (isos->iso_owner_vtx_coord[id_iso][i_part]==PDM_OWNERSHIP_KEEP) {
        free(isos->iso_vtx_coord[id_iso][i_part]);
      }
      if (isos->iso_owner_gnum[id_iso][i_part][PDM_MESH_ENTITY_VTX]==PDM_OWNERSHIP_KEEP) {
        free(isos->iso_vtx_gnum[id_iso][i_part]);
      }
      free(isos->iso_vtx_lparent_idx   [id_iso][i_part]);
      free(isos->iso_vtx_lparent       [id_iso][i_part]);
      if (isos->iso_owner_vtx_parent_weight[id_iso][i_part]==PDM_OWNERSHIP_KEEP) {
        free(isos->iso_vtx_parent_weight[id_iso][i_part]);
      }
      free(isos->isovalue_vtx_idx[id_iso][i_part]);
    }
  }
  if (isos->iso_n_vtx[id_iso]!=NULL) {
    free(isos->iso_n_vtx            [id_iso]);
  }
  if (isos->iso_vtx_coord[id_iso]!=NULL) {
    free(isos->iso_vtx_coord[id_iso]);
  }
  if (isos->iso_vtx_gnum[id_iso]!=NULL) {
    free(isos->iso_vtx_gnum[id_iso]);
  }
  if (isos->iso_vtx_lparent_idx[id_iso]!=NULL) {
    free(isos->iso_vtx_lparent_idx[id_iso]);
  }
  if (isos->iso_vtx_lparent[id_iso]!=NULL) {
    free(isos->iso_vtx_lparent[id_iso]);
  }
  if (isos->iso_vtx_parent_weight[id_iso]!=NULL) {
    free(isos->iso_vtx_parent_weight[id_iso]);
  }
  if (isos->isovalue_vtx_idx[id_iso]!=NULL) {
    free(isos->isovalue_vtx_idx[id_iso]);
  }

  isos->iso_n_vtx            [id_iso] = NULL;
  isos->iso_vtx_coord        [id_iso] = NULL;
  isos->iso_vtx_gnum         [id_iso] = NULL;
  isos->iso_vtx_lparent_idx  [id_iso] = NULL;
  isos->iso_vtx_lparent      [id_iso] = NULL;
  isos->iso_vtx_parent_weight[id_iso] = NULL;
  isos->isovalue_vtx_idx     [id_iso] = NULL;
}


static void
_free_iso_edge
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  if (isos->iso_owner_connec[id_iso]!=NULL) {
    for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
      if (isos->iso_owner_connec[id_iso][i_part][PDM_CONNECTIVITY_TYPE_EDGE_VTX]==PDM_OWNERSHIP_KEEP) {
        free(isos->iso_edge_vtx[id_iso][i_part]);
      }
      if (isos->iso_owner_gnum[id_iso][i_part][PDM_MESH_ENTITY_EDGE]==PDM_OWNERSHIP_KEEP) {
        free(isos->iso_edge_gnum[id_iso][i_part]);
      }
      free(isos->iso_edge_lparent_idx[id_iso][i_part]);
      free(isos->iso_edge_lparent    [id_iso][i_part]);
      free(isos->iso_edge_group_idx  [id_iso][i_part]);
      free(isos->iso_edge_group_lnum [id_iso][i_part]);
      free(isos->iso_edge_group_gnum [id_iso][i_part]);
      free(isos->isovalue_edge_idx   [id_iso][i_part]);
    }
  }
  if (isos->iso_n_edge[id_iso]!=NULL) {
    free(isos->iso_n_edge[id_iso]);
  }
  if (isos->iso_edge_vtx[id_iso]!=NULL) {
    free(isos->iso_edge_vtx[id_iso]);
  }
  if (isos->iso_edge_gnum[id_iso]!=NULL) {
    free(isos->iso_edge_gnum[id_iso]);
  }
  if (isos->iso_edge_lparent_idx[id_iso]!=NULL) {
    free(isos->iso_edge_lparent_idx[id_iso]);
  }
  if (isos->iso_edge_lparent[id_iso]!=NULL) {
    free(isos->iso_edge_lparent[id_iso]);
  }
  if (isos->iso_edge_group_idx[id_iso]!=NULL) {
    free(isos->iso_edge_group_idx[id_iso]);
  }
  if (isos->iso_edge_group_lnum[id_iso]!=NULL) {
    free(isos->iso_edge_group_lnum[id_iso]);
  }
  if (isos->iso_edge_group_gnum[id_iso]!=NULL) {
    free(isos->iso_edge_group_gnum[id_iso]);
  }
  if (isos->isovalue_edge_idx[id_iso]!=NULL) {
    free(isos->isovalue_edge_idx[id_iso]);
  }

  isos->iso_n_edge          [id_iso] = NULL;
  isos->iso_edge_vtx        [id_iso] = NULL;
  isos->iso_edge_gnum       [id_iso] = NULL;
  isos->iso_edge_lparent_idx[id_iso] = NULL;
  isos->iso_edge_lparent    [id_iso] = NULL;
  isos->iso_n_edge_group    [id_iso] = 0;
  isos->iso_edge_group_idx  [id_iso] = NULL;
  isos->iso_edge_group_lnum [id_iso] = NULL;
  isos->iso_edge_group_gnum [id_iso] = NULL;
  isos->isovalue_edge_idx   [id_iso] = NULL;
}


static void
_free_iso_face
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  if (isos->iso_owner_connec[id_iso]!=NULL) {
    for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
      if (isos->iso_owner_connec[id_iso][i_part][PDM_CONNECTIVITY_TYPE_FACE_VTX]==PDM_OWNERSHIP_KEEP) {
        free(isos->iso_face_vtx_idx   [id_iso][i_part]);
        free(isos->iso_face_vtx       [id_iso][i_part]);
      }
      if (isos->iso_owner_gnum[id_iso][i_part][PDM_MESH_ENTITY_FACE]==PDM_OWNERSHIP_KEEP) {
        free(isos->iso_face_gnum      [id_iso][i_part]);
      }
      free(isos->iso_face_lparent_idx[id_iso][i_part]);
      free(isos->iso_face_lparent    [id_iso][i_part]);
      free(isos->isovalue_face_idx   [id_iso][i_part]);
    }
  }

  if (isos->iso_n_face[id_iso]!=NULL) {
    free(isos->iso_n_face[id_iso]);
  }
  if (isos->iso_face_gnum[id_iso]!=NULL) {
    free(isos->iso_face_gnum[id_iso]);
  }
  if (isos->iso_face_vtx_idx[id_iso]!=NULL) {
    free(isos->iso_face_vtx_idx[id_iso]);
  }
  if (isos->iso_face_vtx[id_iso]!=NULL) {
    free(isos->iso_face_vtx[id_iso]);
  }
  if (isos->iso_face_lparent_idx[id_iso]!=NULL) {
    free(isos->iso_face_lparent_idx[id_iso]);
  }
  if (isos->iso_face_lparent[id_iso]!=NULL) {
    free(isos->iso_face_lparent[id_iso]);
  }
  if (isos->isovalue_face_idx[id_iso]!=NULL) {
    free(isos->isovalue_face_idx[id_iso]);
  }

  isos->iso_n_face          [id_iso] = NULL;
  isos->iso_face_vtx_idx    [id_iso] = NULL;
  isos->iso_face_vtx        [id_iso] = NULL;
  isos->iso_face_gnum       [id_iso] = NULL;
  isos->iso_face_lparent_idx[id_iso] = NULL;
  isos->iso_face_lparent    [id_iso] = NULL;
  isos->isovalue_face_idx   [id_iso] = NULL;
}


static void
_free_field
(
  PDM_isosurface_t *isos,
  int               id_iso,
  int               partial
)
{
  if (isos->field[id_iso]!=NULL) {
    for (int i_part=0; i_part<isos->n_part; ++i_part) {
      if (isos->kind[id_iso]!=PDM_ISO_SURFACE_KIND_FIELD && isos->field[id_iso][i_part]!=NULL) {
        free(isos->field[id_iso][i_part]);
      }
    }
    if (!partial) free(isos->field[id_iso]);
  }
  if (!partial) isos->field[id_iso] = NULL;

  if (isos->extract_field[id_iso] != NULL) {
    for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {
      if (isos->extract_field[id_iso][i_part] != NULL) {
        free(isos->extract_field[id_iso][i_part]);
      }
    }
    free(isos->extract_field[id_iso]);
    isos->extract_field[id_iso] = NULL;
  }
}


static void
_free_owner
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  free(isos->iso_owner_vtx_coord        [id_iso]);
  free(isos->iso_owner_vtx_parent_weight[id_iso]);
  for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
    if (isos->iso_owner_gnum[id_iso]!=NULL) {
      free(isos->iso_owner_gnum   [id_iso][i_part]);
    }
    if (isos->iso_owner_connec[id_iso]!=NULL) {
      free(isos->iso_owner_connec [id_iso][i_part]);
    }
    if (isos->iso_owner_lparent[id_iso]!=NULL) {
      free(isos->iso_owner_lparent[id_iso][i_part]);
    }
  }
  if (isos->iso_owner_gnum[id_iso]!=NULL) {
    free(isos->iso_owner_gnum   [id_iso]);
  }
  if (isos->iso_owner_connec[id_iso]!=NULL) {
    free(isos->iso_owner_connec [id_iso]);
  }
  if (isos->iso_owner_lparent[id_iso]!=NULL) {
    free(isos->iso_owner_lparent[id_iso]);
  }
  if (isos->iso_owner_edge_bnd[id_iso]!=NULL) {
    free(isos->iso_owner_edge_bnd[id_iso]);
  }
  isos->iso_owner_vtx_coord        [id_iso] = NULL;
  isos->iso_owner_vtx_parent_weight[id_iso] = NULL;
  isos->iso_owner_gnum             [id_iso] = NULL;
  isos->iso_owner_connec           [id_iso] = NULL;
  isos->iso_owner_lparent          [id_iso] = NULL;
  isos->iso_owner_edge_bnd         [id_iso] = NULL;
  isos->iso_owner_ptp              [id_iso] = NULL;
}



static void
_isosurface_reset
(
 PDM_isosurface_t *isos,
 int               id_isosurface
)
{
  // > Distributed
  if (isos->is_dist_or_part==0) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_reset not implemented for dist_entry\n");
  }
  // > Partitioned
  else if (isos->is_dist_or_part==1) {
    for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
      _free_iso_vtx (isos, id_isosurface);
      _free_iso_edge(isos, id_isosurface);
      _free_iso_face(isos, id_isosurface);
      _free_owner(isos, id_isosurface);
      _free_field(isos, id_isosurface, 1);
    }
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Isosurface is_dist_or_part = %d is invalid.\n", isos->is_dist_or_part);
  }

  PDM_extract_part_free(isos->extrp[id_isosurface]);
}


static void
_isosurface_compute
(
 PDM_isosurface_t *isos,
 int               id_isosurface
)
{
  int debug = 1;

  if (debug==1) {
    log_trace("PDM_isosurface:: compute isosurface n°%d\n", id_isosurface);
  }

  /* Check if edges were provided by the user */
  _do_we_have_edges(isos);

  if (isos->is_dist_or_part == 0) {
    /* Implicit partitioning */
    _dist_to_part(isos);
  }

  /* Evaluate field */
  _compute_iso_field(isos, id_isosurface, 0);

  /* Extract elements of interest */
  _extract(isos, id_isosurface);

  /* Evaluate field on extracted mesh */
  _compute_iso_field(isos, id_isosurface, 1); // hide in '_extract'?

  /* If nodal with sections other than TRIA3 and TETRA4, fall back to ngon */
  _ngonize(isos, id_isosurface);

  /* Build isosurface mesh */
  if (_is_nodal(isos)) {
    PDM_isosurface_marching_algo(isos,
                                 id_isosurface);
  }
  else {
    PDM_isosurface_ngon_algo(isos,
                             id_isosurface);
  }

  // if (isos->is_dist_or_part==0) { // Distributed entry
  //   if (isos->entry_mesh_type<0) {
  //     PDM_error(__FILE__, __LINE__, 0, "Isosurface is_dist_or_part = %d incoherent with isos->entry_mesh_type = %d < 0.\n", isos->is_dist_or_part, isos->entry_mesh_type);
  //   } else if (isos->entry_mesh_type==1) { // Dist mesh alamano

  //   } else if (isos->entry_mesh_type==2) { // Dist mesh

  //   } else if (isos->entry_mesh_type==3) { // Dist mesh nodal

  //   } else {
  //     PDM_error(__FILE__, __LINE__, 0, "Isosurface isos->entry_mesh_type = %d is invalid for distributed entry.\n", isos->entry_mesh_type);
  //   }
  // } else if (isos->is_dist_or_part==1) { // Partitioned entry
  //   if (isos->entry_mesh_type>0) {
  //     PDM_error(__FILE__, __LINE__, 0, "Isosurface is_dist_or_part = %d incoherent with isos->entry_mesh_type = %d > 0.\n", isos->is_dist_or_part, isos->entry_mesh_type);
  //   }
  //   else if (isos->entry_mesh_type == -1 ||
  //            isos->entry_mesh_type == -2) {
  //     PDM_isosurface_ngon_algo(isos,
  //                              id_isosurface);
  //   }

  //   else if (isos->entry_mesh_type==-3) { // Part mesh nodal

  //     if (debug == 1) {
  //       int mesh_dim = PDM_part_mesh_nodal_mesh_dimension_get(isos->pmesh_nodal);
  //       if (mesh_dim>1) {
  //         PDM_part_mesh_nodal_dump_vtk(isos->pmesh_nodal,
  //                                      PDM_GEOMETRY_KIND_SURFACIC,
  //                                      "pmn_surfacic_entry_mesh");
  //       }
  //       if (mesh_dim>2) {
  //         PDM_part_mesh_nodal_dump_vtk(isos->pmesh_nodal,
  //                                      PDM_GEOMETRY_KIND_VOLUMIC,
  //                                      "pmn_volumic_entry_mesh");
  //       }
  //     }

  //     /*
  //      * Compute field once for all to avoid multiple if in code
  //      * and if user function is costful, would be costful once
  //      */
  //     PDM_isosurface_marching_algo(isos,
  //                                  id_isosurface);

  //   } else {
  //     PDM_error(__FILE__, __LINE__, 0, "Isosurface isos->entry_mesh_type = %d is invalid for partitioned entry.\n", isos->entry_mesh_type);
  //   }
  // } else {
  //   PDM_error(__FILE__, __LINE__, 0, "Isosurface is_dist_or_part = %d is invalid.\n", isos->is_dist_or_part);
  // }

  if (isos->is_dist_or_part == 0) {
    /* Block-distribute the isosurface */
    _part_to_dist(isos, id_isosurface);
  }

}

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


void
_check_entry_mesh_coherence
(
 PDM_isosurface_t *isos,
 int               entry_mesh_type
)
{
  if (isos->entry_mesh_type==0) {
    isos->entry_mesh_type=entry_mesh_type;
  }
  else if (isos->entry_mesh_type!=entry_mesh_type) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t:entry_mesh_type already set to %d.\n", isos->entry_mesh_type);
  }
}


PDM_isosurface_t *
PDM_isosurface_create
(
 PDM_MPI_Comm             comm,
 int                      mesh_dimension
 // PDM_Mesh_nodal_elt_t     elt_type
)
{
  PDM_isosurface_t *isos = NULL;
  PDM_malloc(isos, 1, PDM_isosurface_t);

  // > Save communicator
  isos->comm = comm;

  // > Entry mesh information
  isos->is_dist_or_part = -1; 
  isos->entry_mesh_type =  0; 
  isos->entry_mesh_dim  =  mesh_dimension;

  // // > Isosurface mesh information
  // isos->iso_elt_type = elt_type; 
  isos->extract_kind = PDM_EXTRACT_PART_KIND_LOCAL; 
  isos->extrp        = NULL;

  // > Link with entry mesh
  isos->compute_ptp = NULL;

  // > Isosurfaces informations
  isos->n_isosurface = 0;
  isos->n_isovalues     = NULL;
  isos->  isovalues     = NULL;
  isos->kind            = NULL;
  isos->field_function  = NULL;
  isos->eq_coeffs       = NULL;
  isos->use_gradient    = NULL;
  isos->iso_func        = NULL;

  // > Distributed entry data
  isos->dist_to_part_computed = 0;
  isos->btp_vtx      = NULL;
  isos->dmesh        = NULL;
  isos->dmesh_nodal  = NULL;

  isos->distrib_cell = NULL;
  isos->distrib_face = NULL;
  isos->distrib_edge = NULL;
  isos->distrib_vtx  = NULL;

  isos->dvtx_coord     = NULL;
  isos->dcell_face_idx = NULL;
  isos->dcell_face     = NULL;
  isos->dface_edge_idx = NULL;
  isos->dface_edge     = NULL;
  isos->dface_vtx_idx  = NULL;
  isos->dface_vtx      = NULL;
  isos->dedge_vtx      = NULL;

  isos->dgroup_face_idx = NULL;
  isos->dgroup_face     = NULL;
  // isos->dgroup_edge_idx = NULL;
  // isos->dgroup_edge     = NULL;
  // isos->dgroup_vtx_idx  = NULL;
  // isos->dgroup_vtx      = NULL;

  isos->dfield    = NULL;
  isos->dgradient = NULL;

  // > Partitioned entry data
  isos->pmesh       = NULL;
  isos->pmesh_nodal = NULL;

  isos->n_part = -1;
  isos->n_cell = NULL;
  isos->n_face = NULL;
  isos->n_edge = NULL;
  isos->n_vtx  = NULL;

  isos->cell_gnum = NULL;
  isos->face_gnum = NULL;
  isos->edge_gnum = NULL;
  isos-> vtx_gnum = NULL;

  isos->vtx_coord     = NULL;
  isos->cell_face_idx = NULL;
  isos->cell_face     = NULL;
  isos->face_edge_idx = NULL;
  isos->face_edge     = NULL;
  isos->face_vtx_idx  = NULL;
  isos->face_vtx      = NULL;
  isos->edge_vtx      = NULL;

  isos->n_group_face    = NULL;
  isos->group_face_idx  = NULL;
  isos->group_face      = NULL;
  isos->group_face_gnum = NULL;
  // isos->n_group_edge   = NULL;
  // isos->group_edge_idx = NULL;
  // isos->group_edge     = NULL;
  // isos->n_group_vtx    = NULL;
  // isos->group_vtx_idx  = NULL;
  // isos->group_vtx      = NULL;

  isos->we_have_edges = -1;

  isos->field    = NULL;
  isos->gradient = NULL;

  isos->extract_field = NULL;


  // > Partitioned output data
  isos->iso_n_vtx             = NULL;
  isos->iso_vtx_coord         = NULL;
  isos->iso_vtx_gnum          = NULL;
  isos->iso_vtx_lparent_idx   = NULL;
  isos->iso_vtx_lparent       = NULL;
  isos->iso_vtx_parent_weight = NULL;
  isos->isovalue_vtx_idx      = NULL;

  isos->iso_n_edge            = NULL;
  isos->iso_edge_vtx          = NULL;
  isos->iso_edge_gnum         = NULL;
  isos->iso_edge_lparent_idx  = NULL;
  isos->iso_edge_lparent      = NULL;
  isos->iso_n_edge_group      = NULL;
  isos->iso_edge_group_idx    = NULL;
  isos->iso_edge_group_lnum   = NULL;
  isos->iso_edge_group_gnum   = NULL;
  isos->isovalue_edge_idx     = NULL;

  isos->iso_n_face            = NULL;
  isos->iso_face_vtx_idx      = NULL;
  isos->iso_face_vtx          = NULL;
  isos->iso_face_gnum         = NULL;
  isos->iso_face_lparent_idx  = NULL;
  isos->iso_face_lparent      = NULL;
  isos->isovalue_face_idx     = NULL;

  // > Part_to_part between iso entities and entry mesh entites
  isos->iso_ptp_vtx  = NULL;
  isos->iso_ptp_edge = NULL;
  isos->iso_ptp_face = NULL;

  // > Owners
  isos->iso_owner_vtx_coord         = NULL;
  isos->iso_owner_vtx_parent_weight = NULL;
  isos->iso_owner_gnum              = NULL;
  isos->iso_owner_connec            = NULL;
  isos->iso_owner_lparent           = NULL;
  isos->iso_owner_edge_bnd          = NULL;
  isos->iso_owner_ptp               = NULL;


  return isos;
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
  // TODO: check that difference between isovalues>ISOSURFACE_EPS
  
  int id_isosurface = isos->n_isosurface;
  isos->n_isosurface++;
  
  isos->kind           = realloc(isos->kind          , sizeof(PDM_iso_surface_kind_t          ) * isos->n_isosurface);
  isos->eq_coeffs      = realloc(isos->eq_coeffs     , sizeof(double                         *) * isos->n_isosurface);
  isos->field_function = realloc(isos->field_function, sizeof(PDM_isosurface_field_function_t ) * isos->n_isosurface);
  
  isos->field_function = realloc(isos->field_function, sizeof(PDM_isosurface_field_function_t ) * isos->n_isosurface);
  
  isos->compute_ptp = realloc(isos->compute_ptp, sizeof(int *) * isos->n_isosurface);
  isos->compute_ptp[id_isosurface] = PDM_array_zeros_int(PDM_MESH_ENTITY_MAX);

  // TODO: status if its ok to allocate d and p variable (because set functions order not imposed ?)
  isos->dfield         = realloc(isos->dfield        , sizeof(double                         *) * isos->n_isosurface);
  isos->field          = realloc(isos->field         , sizeof(double                        **) * isos->n_isosurface);
  isos->field[id_isosurface] = NULL;

  isos->extract_field = realloc(isos->extract_field, sizeof(double **) * isos->n_isosurface);
  isos->extract_field[id_isosurface] = NULL;

  if (isos->is_dist_or_part == 1) {
    // Partitioned
    int n_part = 0;
    if (isos->entry_mesh_type==-1) {
      n_part = isos->n_part;
      if (n_part==-1) {
        PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: entry_mesh_type = %d but n_part isn't defined.\n", isos->entry_mesh_type);
      }
    }
    else if (isos->entry_mesh_type==-2) {
      n_part = PDM_part_mesh_n_part_get(isos->pmesh);
    }
    else if (isos->entry_mesh_type==-3) {
      n_part = PDM_part_mesh_nodal_n_part_get(isos->pmesh_nodal);
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: Impossible to defined manually field without setting mesh first.\n", isos->entry_mesh_type);
    }
    PDM_malloc(isos->field[id_isosurface], n_part, double *);
  }

  isos->n_isovalues    = realloc(isos->n_isovalues   , sizeof(int     ) * isos->n_isosurface);
  isos->isovalues      = realloc(isos->isovalues     , sizeof(double *) * isos->n_isosurface);
  isos->use_gradient   = realloc(isos->use_gradient  , sizeof(int    *) * isos->n_isosurface);

  isos->extrp = realloc(isos->extrp, sizeof(PDM_extract_part_t *) * isos->n_isosurface);
  isos->extrp[id_isosurface] = NULL;


  // > Partitioned iso vertices
  isos->iso_n_vtx             = realloc(isos->iso_n_vtx            , sizeof(int           *) * isos->n_isosurface);
  isos->iso_vtx_coord         = realloc(isos->iso_vtx_coord        , sizeof(double       **) * isos->n_isosurface);
  isos->iso_vtx_gnum          = realloc(isos->iso_vtx_gnum         , sizeof(PDM_g_num_t  **) * isos->n_isosurface);
  isos->iso_vtx_lparent_idx   = realloc(isos->iso_vtx_lparent_idx  , sizeof(int          **) * isos->n_isosurface);
  isos->iso_vtx_lparent       = realloc(isos->iso_vtx_lparent      , sizeof(int          **) * isos->n_isosurface);
  isos->iso_vtx_parent_weight = realloc(isos->iso_vtx_parent_weight, sizeof(double       **) * isos->n_isosurface);
  isos->isovalue_vtx_idx      = realloc(isos->isovalue_vtx_idx     , sizeof(int          **) * isos->n_isosurface);
  isos->iso_n_vtx            [id_isosurface] = NULL;
  isos->iso_vtx_coord        [id_isosurface] = NULL;
  isos->iso_vtx_gnum         [id_isosurface] = NULL;
  isos->iso_vtx_lparent_idx  [id_isosurface] = NULL;
  isos->iso_vtx_lparent      [id_isosurface] = NULL;
  isos->iso_vtx_parent_weight[id_isosurface] = NULL;
  isos->isovalue_vtx_idx     [id_isosurface] = NULL;

  // > Partitioned iso edges
  isos->iso_n_edge           = realloc(isos->iso_n_edge          , sizeof(int           *) * isos->n_isosurface);
  isos->iso_edge_vtx         = realloc(isos->iso_edge_vtx        , sizeof(int          **) * isos->n_isosurface);
  isos->iso_edge_gnum        = realloc(isos->iso_edge_gnum       , sizeof(PDM_g_num_t  **) * isos->n_isosurface);
  isos->iso_edge_lparent_idx = realloc(isos->iso_edge_lparent_idx, sizeof(int          **) * isos->n_isosurface);
  isos->iso_edge_lparent     = realloc(isos->iso_edge_lparent    , sizeof(int          **) * isos->n_isosurface);
  isos->iso_n_edge_group     = realloc(isos->iso_n_edge_group    , sizeof(int            ) * isos->n_isosurface);
  isos->iso_edge_group_idx   = realloc(isos->iso_edge_group_idx  , sizeof(int          **) * isos->n_isosurface);
  isos->iso_edge_group_lnum  = realloc(isos->iso_edge_group_lnum , sizeof(int          **) * isos->n_isosurface);
  isos->iso_edge_group_gnum  = realloc(isos->iso_edge_group_gnum , sizeof(PDM_g_num_t  **) * isos->n_isosurface);
  isos->isovalue_edge_idx    = realloc(isos->isovalue_edge_idx   , sizeof(int          **) * isos->n_isosurface);
  isos->iso_n_edge          [id_isosurface] = NULL;
  isos->iso_edge_vtx        [id_isosurface] = NULL;
  isos->iso_edge_gnum       [id_isosurface] = NULL;
  isos->iso_edge_lparent_idx[id_isosurface] = NULL;
  isos->iso_edge_lparent    [id_isosurface] = NULL;
  isos->iso_n_edge_group    [id_isosurface] = 0;
  isos->iso_edge_group_idx  [id_isosurface] = NULL;
  isos->iso_edge_group_lnum [id_isosurface] = NULL;
  isos->iso_edge_group_gnum [id_isosurface] = NULL;
  isos->isovalue_edge_idx   [id_isosurface] = NULL;

  // > Partitioned iso faces
  isos->iso_n_face           = realloc(isos->iso_n_face          , sizeof(int           *) * isos->n_isosurface);
  isos->iso_face_vtx_idx     = realloc(isos->iso_face_vtx_idx    , sizeof(int          **) * isos->n_isosurface);
  isos->iso_face_vtx         = realloc(isos->iso_face_vtx        , sizeof(int          **) * isos->n_isosurface);
  isos->iso_face_gnum        = realloc(isos->iso_face_gnum       , sizeof(PDM_g_num_t  **) * isos->n_isosurface);
  isos->iso_face_lparent_idx = realloc(isos->iso_face_lparent_idx, sizeof(int          **) * isos->n_isosurface);
  isos->iso_face_lparent     = realloc(isos->iso_face_lparent    , sizeof(int          **) * isos->n_isosurface);
  isos->isovalue_face_idx    = realloc(isos->isovalue_face_idx   , sizeof(int          **) * isos->n_isosurface);
  isos->iso_n_face          [id_isosurface] = NULL;
  isos->iso_face_vtx_idx    [id_isosurface] = NULL;
  isos->iso_face_vtx        [id_isosurface] = NULL;
  isos->iso_face_gnum       [id_isosurface] = NULL;
  isos->iso_face_lparent_idx[id_isosurface] = NULL;
  isos->iso_face_lparent    [id_isosurface] = NULL;
  isos->isovalue_face_idx   [id_isosurface] = NULL;

  // > Part_to_part between iso entities and entry mesh entites
  isos->iso_ptp_vtx  = realloc(isos->iso_ptp_vtx , sizeof(PDM_part_to_part_t *) * isos->n_isosurface);
  isos->iso_ptp_edge = realloc(isos->iso_ptp_edge, sizeof(PDM_part_to_part_t *) * isos->n_isosurface);
  isos->iso_ptp_face = realloc(isos->iso_ptp_face, sizeof(PDM_part_to_part_t *) * isos->n_isosurface);

  // > Partitioned owners
  isos->iso_owner_vtx_coord         = realloc(isos->iso_owner_vtx_coord         , sizeof(PDM_ownership_t  *) * isos->n_isosurface);
  isos->iso_owner_vtx_parent_weight = realloc(isos->iso_owner_vtx_parent_weight , sizeof(PDM_ownership_t  *) * isos->n_isosurface);
  isos->iso_owner_gnum              = realloc(isos->iso_owner_gnum              , sizeof(PDM_ownership_t **) * isos->n_isosurface);
  isos->iso_owner_connec            = realloc(isos->iso_owner_connec            , sizeof(PDM_ownership_t **) * isos->n_isosurface);
  isos->iso_owner_lparent           = realloc(isos->iso_owner_lparent           , sizeof(PDM_ownership_t **) * isos->n_isosurface);
  isos->iso_owner_edge_bnd          = realloc(isos->iso_owner_edge_bnd          , sizeof(PDM_ownership_t **) * isos->n_isosurface);
  isos->iso_owner_ptp               = realloc(isos->iso_owner_ptp               , sizeof(PDM_ownership_t  *) * isos->n_isosurface);
  isos->iso_owner_vtx_coord        [id_isosurface] = NULL;
  isos->iso_owner_vtx_parent_weight[id_isosurface] = NULL;
  isos->iso_owner_gnum             [id_isosurface] = NULL;
  isos->iso_owner_connec           [id_isosurface] = NULL;
  isos->iso_owner_lparent          [id_isosurface] = NULL;
  isos->iso_owner_edge_bnd         [id_isosurface] = NULL;
  isos->iso_owner_ptp              [id_isosurface] = NULL;

  isos->kind       [id_isosurface] = kind;
  isos->n_isovalues[id_isosurface] = n_isovalues;
  PDM_malloc(isos->isovalues[id_isosurface], n_isovalues, double);
  for (int i=0; i<n_isovalues; ++i) {
    isos->isovalues[id_isosurface][i] = isovalues[i];
  }

  if (kind==PDM_ISO_SURFACE_KIND_FIELD) {
    isos->eq_coeffs[id_isosurface] = NULL;
  }
  else if (kind==PDM_ISO_SURFACE_KIND_PLANE) {
    int n_coeff = 4;
    PDM_malloc(isos->eq_coeffs[id_isosurface], n_coeff, double);
  }
  else if (kind==PDM_ISO_SURFACE_KIND_SPHERE) {
    int n_coeff = 4;
    PDM_malloc(isos->eq_coeffs[id_isosurface], n_coeff, double);
  }
  else if (kind==PDM_ISO_SURFACE_KIND_ELLIPSE) {
    int n_coeff = 6;
    PDM_malloc(isos->eq_coeffs[id_isosurface], n_coeff, double);
  }
  else if (kind==PDM_ISO_SURFACE_KIND_QUADRIC) {
    int n_coeff = 10;
    PDM_malloc(isos->eq_coeffs[id_isosurface], n_coeff, double);
  }
  else if (kind==PDM_ISO_SURFACE_KIND_HEART) {
    isos->eq_coeffs[id_isosurface] = NULL;
  }
  else if (kind==PDM_ISO_SURFACE_KIND_FUNCTION) {
    isos->eq_coeffs[id_isosurface] = NULL;
  }

  return id_isosurface;
}


void
PDM_isosurface_equation_set
(
 PDM_isosurface_t *isos,
 int               id_isosurface,
 double           *coeff,
 int               use_gradient
)
{
  if        (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_PLANE) {
    int n_coeff = 4;
    for (int i_coeff=0; i_coeff<n_coeff; ++i_coeff) {
      isos->eq_coeffs[id_isosurface][i_coeff] = coeff[i_coeff];
    }
  }
  else if (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_SPHERE) {
    int n_coeff = 4;
    for (int i_coeff=0; i_coeff<n_coeff; ++i_coeff) {
      isos->eq_coeffs[id_isosurface][i_coeff] = coeff[i_coeff];
    }
  }
  else if (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_ELLIPSE) {
    int n_coeff = 6;
    for (int i_coeff=0; i_coeff<n_coeff; ++i_coeff) {
      isos->eq_coeffs[id_isosurface][i_coeff] = coeff[i_coeff];
    }
  }
  else if (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_QUADRIC) {
    int n_coeff = 10;
    for (int i_coeff=0; i_coeff<n_coeff; ++i_coeff) {
      isos->eq_coeffs[id_isosurface][i_coeff] = coeff[i_coeff];
    }
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Isosurface n°%d doesn't support PDM_isosurface_equation_set method cause its kind is %d.\n", id_isosurface, isos->kind[id_isosurface]);
  }

  isos->use_gradient[id_isosurface] = use_gradient;

}


void
PDM_isosurface_field_function_set
(
 PDM_isosurface_t                *isos,
 int                              id_isosurface,
 PDM_isosurface_field_function_t  func
)
{
  if (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_FUNCTION) {
    isos->field_function[id_isosurface] = func;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Isosurface n°%d doesn't support PDM_isosurface_field_function_set method cause its kind is %d.\n", id_isosurface, isos->kind[id_isosurface]);
  }

}


void
PDM_isosurface_reset
(
 PDM_isosurface_t *isos,
 int               id_isosurface
)
{
  if (id_isosurface < 0) {
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
  if (id_isosurface < 0) {
    for (int i = 0; i < isos->n_isosurface; i++) {
      _isosurface_compute(isos, i);
    }
  }
  else {
    _isosurface_compute(isos, id_isosurface);
  }
}


void
PDM_isosurface_dump_times
(
 PDM_isosurface_t *isos
)
{
  PDM_UNUSED(isos);
}


void
PDM_isosurface_free
(
  PDM_isosurface_t  *isos
)
{
  if (isos->is_dist_or_part == 0) {
    // Block-distributed
    PDM_block_to_part_free(isos->btp_vtx);

    if (isos->cell_face      [0] != NULL) free(isos->cell_face      [0]);
    if (isos->cell_face_idx  [0] != NULL) free(isos->cell_face_idx  [0]);
    if (isos->face_edge      [0] != NULL) free(isos->face_edge      [0]);
    if (isos->face_edge_idx  [0] != NULL) free(isos->face_edge_idx  [0]);
    if (isos->face_vtx       [0] != NULL) free(isos->face_vtx       [0]);
    if (isos->face_vtx_idx   [0] != NULL) free(isos->face_vtx_idx   [0]);
    if (isos->edge_vtx       [0] != NULL) free(isos->edge_vtx       [0]);
    if (isos->vtx_coord      [0] != NULL) free(isos->vtx_coord      [0]);
    if (isos->cell_gnum      [0] != NULL) free(isos->cell_gnum      [0]);
    if (isos->face_gnum      [0] != NULL) free(isos->face_gnum      [0]);
    if (isos->edge_gnum      [0] != NULL) free(isos->edge_gnum      [0]);
    if (isos->vtx_gnum       [0] != NULL) free(isos->vtx_gnum       [0]);
    if (isos->group_face_idx [0] != NULL) free(isos->group_face_idx [0]);
    if (isos->group_face     [0] != NULL) free(isos->group_face     [0]);
    if (isos->group_face_gnum[0] != NULL) free(isos->group_face_gnum[0]);
  }
  if (isos->n_cell          != NULL) free(isos->n_cell         );
  if (isos->n_face          != NULL) free(isos->n_face         );
  if (isos->n_edge          != NULL) free(isos->n_edge         );
  if (isos->n_vtx           != NULL) free(isos->n_vtx          );
  if (isos->cell_face       != NULL) free(isos->cell_face      );
  if (isos->cell_face_idx   != NULL) free(isos->cell_face_idx  );
  if (isos->face_edge       != NULL) free(isos->face_edge      );
  if (isos->face_edge_idx   != NULL) free(isos->face_edge_idx  );
  if (isos->face_vtx        != NULL) free(isos->face_vtx       );
  if (isos->face_vtx_idx    != NULL) free(isos->face_vtx_idx   );
  if (isos->edge_vtx        != NULL) free(isos->edge_vtx       );
  if (isos->vtx_coord       != NULL) free(isos->vtx_coord      );
  if (isos->cell_gnum       != NULL) free(isos->cell_gnum      );
  if (isos->face_gnum       != NULL) free(isos->face_gnum      );
  if (isos->edge_gnum       != NULL) free(isos->edge_gnum      );
  if (isos->vtx_gnum        != NULL) free(isos->vtx_gnum       );
  if (isos->n_group_face    != NULL) free(isos->n_group_face   );
  if (isos->group_face_idx  != NULL) free(isos->group_face_idx );
  if (isos->group_face      != NULL) free(isos->group_face     );
  if (isos->group_face_gnum != NULL) free(isos->group_face_gnum);


  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    if (isos->n_isovalues[id_iso]>0) {
      free(isos->isovalues[id_iso]);
      if (isos->kind[id_iso]==PDM_ISO_SURFACE_KIND_FIELD) {
        if (isos->is_dist_or_part==0) { // is distributed
          free(isos->dfield[id_iso]);
        }
      } else  {
        free(isos->eq_coeffs[id_iso]);
      }
    }
  }

  free(isos->n_isovalues);
  free(isos->isovalues);
  free(isos->eq_coeffs);
  free(isos->use_gradient);
  free(isos->field_function);
  free(isos->dfield);
  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    _free_field(isos, id_iso, 0);
  }
  free(isos->field);
  free(isos->kind);

  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    _free_iso_vtx (isos, id_iso);
    _free_iso_edge(isos, id_iso);
    _free_iso_face(isos, id_iso);
    PDM_extract_part_free(isos->extrp[id_iso]);
  }
  free(isos->extract_field);
  free(isos->extrp);

  free(isos->iso_n_vtx);
  free(isos->iso_vtx_coord);
  free(isos->iso_vtx_gnum);
  free(isos->iso_vtx_lparent_idx);
  free(isos->iso_vtx_lparent);
  free(isos->iso_vtx_parent_weight);
  free(isos->isovalue_vtx_idx);
  
  free(isos->iso_n_edge);
  free(isos->iso_edge_vtx);
  free(isos->iso_edge_gnum);
  free(isos->iso_edge_lparent_idx);
  free(isos->iso_edge_lparent);
  free(isos->iso_n_edge_group);
  free(isos->iso_edge_group_idx);
  free(isos->iso_edge_group_lnum);
  free(isos->iso_edge_group_gnum);
  free(isos->isovalue_edge_idx);

  free(isos->iso_n_face);
  free(isos->iso_face_vtx_idx);
  free(isos->iso_face_vtx);
  free(isos->iso_face_gnum);
  free(isos->iso_face_lparent_idx);
  free(isos->iso_face_lparent);
  free(isos->isovalue_face_idx);


  // > Part_to_part between iso entities and entry mesh entites
  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    if (isos->  compute_ptp[id_iso][PDM_MESH_ENTITY_VTX]==1 &&
        isos->iso_owner_ptp[id_iso][PDM_MESH_ENTITY_VTX]==PDM_OWNERSHIP_KEEP) {
      PDM_part_to_part_free(isos->iso_ptp_vtx[id_iso]);
    }
    else if (isos->  compute_ptp[id_iso][PDM_MESH_ENTITY_EDGE]==1 &&
             isos->iso_owner_ptp[id_iso][PDM_MESH_ENTITY_EDGE]==PDM_OWNERSHIP_KEEP) {
      PDM_part_to_part_free(isos->iso_ptp_edge[id_iso]);
    }
    else if (isos->  compute_ptp[id_iso][PDM_MESH_ENTITY_FACE]==1 &&
             isos->iso_owner_ptp[id_iso][PDM_MESH_ENTITY_FACE]==PDM_OWNERSHIP_KEEP) {
      PDM_part_to_part_free(isos->iso_ptp_face[id_iso]);
    }
  }

  free(isos->iso_ptp_vtx);
  free(isos->iso_ptp_edge);
  free(isos->iso_ptp_face);



  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    _free_owner(isos, id_iso);
  }
  free(isos->iso_owner_vtx_coord);
  free(isos->iso_owner_vtx_parent_weight);
  free(isos->iso_owner_gnum);
  free(isos->iso_owner_connec);
  free(isos->iso_owner_lparent);
  free(isos->iso_owner_edge_bnd);
  free(isos->iso_owner_ptp);



  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    free(isos->compute_ptp[id_iso]);
  }
  free(isos->compute_ptp);


  free(isos);
}


#ifdef  __cplusplus
}
#endif
