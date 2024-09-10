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
                                      isos->n_dgroup_face,
                                      isos->dgroup_face_idx,
                                      isos->dgroup_face,
                                      isos->n_part,
                                      &n_face,
                                      pface_ln_to_gn,
                                      &isos->group_face_idx,
                                      &isos->group_face,
                                      &isos->group_face_gnum);
    }

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

  if (isos->kind[id_isosurface] == PDM_ISO_SURFACE_KIND_FIELD) {
    if (use_extract) {
      if (!_is_nodal(isos)) {
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
        // > For now nodal has not extract part so we fake it
        PDM_malloc(isos->extract_field[id_isosurface], isos->n_part, double *);
        for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {
          PDM_malloc(isos->extract_field[id_isosurface][i_part], isos->extract_n_vtx[i_part], double);
          for (int i = 0; i < isos->extract_n_vtx[i_part]; i++) {
            isos->extract_field[id_isosurface][i_part][i] = isos->field[id_isosurface][i_part][i];
          }
        }
        return;
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

  if (isos->is_dist_or_part == 0) {
    // Block-distributed
    if (isos->field[id_isosurface] == NULL) {
      PDM_malloc(isos->field[id_isosurface], isos->n_part, double *);
    }
    assert(isos->dist_to_part_computed);
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
      if (!_is_nodal(isos)) {
        n_vtx[i_part] = PDM_extract_part_vtx_coord_get(isos->extrp[id_isosurface],
                                                       i_part,
                                                       &vtx_coord[i_part],
                                                       PDM_OWNERSHIP_KEEP);
      }
      else {
        n_vtx    [i_part] = isos->extract_n_vtx    [i_part];
        vtx_coord[i_part] = isos->extract_vtx_coord[i_part];
      }
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

  PDM_free(n_vtx    );
  PDM_free(vtx_coord);
}



// --->> migrer dans priv?
static inline int
_sign
(
  const double v, 
  const double tol
)
{
  // if (v < -tol) {
  //   return -1;
  // }
  // else if (v > tol) {
  //   return 1;
  // }
  // else {
  //   return 0;
  // }
  return (v > tol);
}


static inline int
_cross_0_level_ngon
(
  const double v0,
  const double v1,
  const double tol
)
{
  return _sign(v0, tol) != _sign(v1, tol);
}


static inline int
_cross_any_level_ngon
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
    n_crossings += _cross_0_level_ngon(v0 - isovalues[i], v1 - isovalues[i], tol);
  }

  return n_crossings;
}
// <<---



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
  int debug = 0;

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
  if (debug==1) {
    PDM_log_trace_array_int(elt_bnd_tag, n_elt, "elt_bnd_tag ::");
  }
  return n_group;
}


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
    isos->iso_n_part   = 1;
  }


  if (_is_nodal(isos)) {
    // // Nodal
    // PDM_malloc(pn_cell       , isos->n_part, int          );
    // PDM_malloc(pn_face       , isos->n_part, int          );
    // PDM_malloc(pn_edge       , isos->n_part, int          );
    // PDM_malloc(pn_vtx        , isos->n_part, int          );
    // PDM_malloc(pcell_ln_to_gn, isos->n_part, PDM_g_num_t *);
    // PDM_malloc(pface_ln_to_gn, isos->n_part, PDM_g_num_t *);
    // PDM_malloc(pvtx_ln_to_gn , isos->n_part, PDM_g_num_t *);

    // if (isos->entry_mesh_dim == 2) {
    //   for (int i_part = 0; i_part < isos->n_part; i_part++) {
    //     pn_cell[i_part] = 0;
    //   }
    // }

    printf("WARNING: active sub mesh extraction not implemented for nodal. Algo may be slow...\n");

    if (isos->extract_n_vtx==NULL) {
      PDM_malloc(isos->extract_n_vtx    , isos->n_part, int          );
      PDM_malloc(isos->extract_vtx_coord, isos->n_part, double      *);
      PDM_malloc(isos->extract_vtx_gnum , isos->n_part, PDM_g_num_t *);
      PDM_malloc(isos->extract_vtx_lnum , isos->n_part, int         *);

      PDM_malloc(isos->extract_n_tri      , isos->n_part, int          );
      PDM_malloc(isos->extract_tri_vtx    , isos->n_part, int         *);
      PDM_malloc(isos->extract_tri_gnum   , isos->n_part, PDM_g_num_t *);
      PDM_malloc(isos->extract_tri_lnum   , isos->n_part, int         *);
      PDM_malloc(isos->extract_tri_n_group, isos->n_part, int          );
      PDM_malloc(isos->extract_tri_tag    , isos->n_part, int         *);

      PDM_malloc(isos->extract_n_tetra   , isos->n_part, int          );
      PDM_malloc(isos->extract_tetra_vtx , isos->n_part, int         *);
      PDM_malloc(isos->extract_tetra_gnum, isos->n_part, PDM_g_num_t *);
      PDM_malloc(isos->extract_tetra_lnum, isos->n_part, int         *);
      memset(isos->extract_n_tetra   , 0, isos->n_part*sizeof(int          ));
      memset(isos->extract_tetra_vtx , 0, isos->n_part*sizeof(int         *));
      memset(isos->extract_tetra_gnum, 0, isos->n_part*sizeof(PDM_g_num_t *));
      memset(isos->extract_tetra_lnum, 0, isos->n_part*sizeof(int         *));
    }

    int n_section = PDM_part_mesh_nodal_n_section_get(isos->pmesh_nodal);
    for (int i_part = 0; i_part < isos->n_part; i_part++) {

      isos->extract_n_vtx    [i_part] = PDM_part_mesh_nodal_n_vtx_get    (isos->pmesh_nodal, i_part);
      isos->extract_vtx_coord[i_part] = PDM_part_mesh_nodal_vtx_coord_get(isos->pmesh_nodal, i_part);
      isos->extract_vtx_gnum [i_part] = PDM_part_mesh_nodal_vtx_g_num_get(isos->pmesh_nodal, i_part);
      isos->extract_vtx_lnum [i_part] = PDM_array_new_arange_int(1, isos->extract_n_vtx[i_part], 1);

      for (int i_section = 0; i_section < n_section; i_section++) {
        PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_section_elt_type_get(isos->pmesh_nodal, i_section);

        int n_elt = PDM_part_mesh_nodal_section_n_elt_get(isos->pmesh_nodal, i_section, i_part);

        int         *connec;
        PDM_g_num_t *gnum;
        int         *parent_num;
        PDM_g_num_t *parent_entity_g_num;
        PDM_part_mesh_nodal_section_std_get(isos->pmesh_nodal,
                                            i_section,
                                            i_part,
                                            &connec,
                                            &gnum,
                                            &parent_num,
                                            &parent_entity_g_num,
                                            PDM_OWNERSHIP_BAD_VALUE);
        switch (t_elt) {

          case PDM_MESH_NODAL_POINT: {
            break;
          }
          case PDM_MESH_NODAL_BAR2: {
            break;
          }
          case PDM_MESH_NODAL_TRIA3: {
            isos->extract_n_tri      [i_part] = n_elt;
            isos->extract_tri_vtx    [i_part] = connec;
            isos->extract_tri_gnum   [i_part] = gnum;
            isos->extract_tri_lnum   [i_part] = PDM_array_new_arange_int(1, n_elt, 1);
            isos->extract_tri_n_group[i_part] = _convert_group_info_to_tag(isos->pmesh_nodal, i_part, n_elt,
                                                                           PDM_GEOMETRY_KIND_SURFACIC,
                                                                          &isos->extract_tri_tag[i_part]);
            break;
          }
          case PDM_MESH_NODAL_TETRA4: {
            isos->extract_n_tetra    [i_part] = n_elt;
            isos->extract_tetra_vtx  [i_part] = connec;
            isos->extract_tetra_gnum [i_part] = gnum;
            isos->extract_tetra_lnum [i_part] = PDM_array_new_arange_int(1, n_elt, 1);
            break;
          }
          default:{
            PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: only full tri and tetra nodal meshes are managed for now (section with element type = %d)\n", t_elt);
          }
        }
      }  // End loop on sections
    } // End loop on partitions
    // TODO...
    // PDM_extract_part_part_nodal_set...
    // PDM_error(__FILE__, __LINE__, 0, "_extract Nodal not implemented yet\n");
  }
  else { // Ngon

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
              if (_cross_any_level_ngon(val0, val1, isos->n_isovalues[id_isosurface], isos->isovalues[id_isosurface], isos->ISOSURFACE_EPS)) {
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
              if (_cross_any_level_ngon(val0, val1, isos->n_isovalues[id_isosurface], isos->isovalues[id_isosurface], isos->ISOSURFACE_EPS)) {
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
                if (_cross_any_level_ngon(val0, val1, isos->n_isovalues[id_isosurface], isos->isovalues[id_isosurface], isos->ISOSURFACE_EPS)) {
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
                if (_cross_any_level_ngon(val0, val1, isos->n_isovalues[id_isosurface], isos->isovalues[id_isosurface], isos->ISOSURFACE_EPS)) {
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
      PDM_extract_part_n_group_set(extrp, PDM_BOUND_TYPE_FACE, isos->n_group_face[0]);
    }

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

      if (isos->entry_mesh_dim==3) {
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
    }

    PDM_extract_part_compute(extrp);

    for (int i_part = 0; i_part < isos->n_part; i_part++) {
      PDM_free(extract_lnum[i_part]);
    }
    PDM_free(n_extract   );
    PDM_free(extract_lnum);
  } // End Ngon
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
 * \brief Block-distribute isosurface mesh vertices
 */
static void
_part_to_dist_vtx
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  // > Set vtx in block frame
  PDM_part_to_block_t *ptb_vtx = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                          1.,
                                                          isos->iso_entity_gnum[PDM_MESH_ENTITY_VTX][id_iso],
                                                          NULL,
                                                          isos->iso_n_entity[PDM_MESH_ENTITY_VTX][id_iso],
                                                          isos->iso_n_part, // will fail if n_part>1 because of extract_part ?
                                                          isos->comm);
  // > Get vtx coords
  double *dvtx_coord = NULL;
  PDM_part_to_block_exch(ptb_vtx,
                         sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         3,
                         NULL,
              (void **)  isos->iso_vtx_coord[id_iso],
                         NULL,
              (void **) &dvtx_coord);
  

  // > Get vtx parent
  int **pvtx_parent_strd = NULL;
  PDM_malloc(pvtx_parent_strd, isos->iso_n_part, int *);
  for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
    pvtx_parent_strd[i_part] = PDM_array_new_size_from_idx_int(isos->iso_entity_parent_idx[PDM_MESH_ENTITY_VTX][id_iso][i_part], isos->iso_n_entity[PDM_MESH_ENTITY_VTX][id_iso][i_part]);
  }

  int          *dvtx_parent_strd = NULL;
  PDM_g_num_t  *dvtx_parent_gnum = NULL;
  PDM_part_to_block_exch(ptb_vtx,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         pvtx_parent_strd,
              (void **)  isos->iso_entity_parent_gnum[PDM_MESH_ENTITY_VTX][id_iso],
                        &dvtx_parent_strd,
              (void **) &dvtx_parent_gnum);
  PDM_free(dvtx_parent_strd);

  double  *dvtx_parent_weight = NULL;
  PDM_part_to_block_exch(ptb_vtx,
                         sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         pvtx_parent_strd,
              (void **)  isos->iso_vtx_parent_weight[id_iso],
                        &dvtx_parent_strd,
              (void **) &dvtx_parent_weight);


  int dn_vtx = PDM_part_to_block_n_elt_block_get(ptb_vtx);
  // PDM_log_trace_array_int(dvtx_parent_strd, dn_vtx, "dvtx_parent_strd");
  int *dvtx_parent_idx = PDM_array_new_idx_from_sizes_int(dvtx_parent_strd, dn_vtx);
  for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
    PDM_free(pvtx_parent_strd[i_part]);
  }
  PDM_free(pvtx_parent_strd);
  PDM_free(dvtx_parent_strd);

  isos->iso_dn_entity          [PDM_MESH_ENTITY_VTX][id_iso] = dn_vtx;
  isos->iso_dentity_parent_idx [PDM_MESH_ENTITY_VTX][id_iso] = dvtx_parent_idx;
  isos->iso_dentity_parent_gnum[PDM_MESH_ENTITY_VTX][id_iso] = dvtx_parent_gnum;
  isos->iso_dvtx_coord        [id_iso] = dvtx_coord;
  isos->iso_dvtx_parent_weight[id_iso] = dvtx_parent_weight;
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
  int debug = 0;

  if (_is_nodal(isos)==1) {
    // TODO: need to manage incoherent group across procs
    PDM_error(__FILE__, __LINE__, 0, "Not implemented yet.\n");
  }


  // > Allocate tmp result
  int         *_iso_dedge_group_idx  = NULL;
  PDM_g_num_t *_iso_dedge_group_gnum = NULL;
  PDM_malloc(_iso_dedge_group_idx , isos->iso_n_edge_group[id_iso]+1, int        );
  PDM_malloc(_iso_dedge_group_gnum, isos->iso_n_edge_group[id_iso]  , PDM_g_num_t);
  _iso_dedge_group_idx[0] = 0;


  // > Exchange info with a part_to_block for each group
  int          *n_elt_group     = NULL;
  PDM_g_num_t **    _group_gnum = NULL;
  PDM_g_num_t **_elt_group_gnum = NULL;
  PDM_malloc(n_elt_group    , isos->n_part, int          );
  PDM_malloc(_elt_group_gnum, isos->n_part, PDM_g_num_t *);
  PDM_malloc(    _group_gnum, isos->n_part, PDM_g_num_t *);

  for (int i_group=0; i_group<isos->iso_n_edge_group[id_iso]; ++i_group) {

    // > Prepare ptb
    for (int i_part=0; i_part<isos->n_part; ++i_part) {
      int i_beg_group = isos->iso_edge_group_idx[id_iso][i_part][i_group  ];
      int i_end_group = isos->iso_edge_group_idx[id_iso][i_part][i_group+1];
      n_elt_group[i_part] = i_end_group-i_beg_group;

      // > Prepare element group gnum
      PDM_malloc(_elt_group_gnum[i_part], n_elt_group[i_part], PDM_g_num_t);

      int i_write = 0;
      for (int i_elmt=i_beg_group; i_elmt<i_end_group; ++i_elmt) {
        int elmt_lnum = isos->iso_edge_group_lnum[id_iso][i_part][i_elmt];
        _elt_group_gnum[i_part][i_write++] = isos->iso_entity_gnum[PDM_MESH_ENTITY_EDGE][id_iso][i_part][elmt_lnum-1];
      }

      _group_gnum[i_part] = &isos->iso_edge_group_gnum[id_iso][i_part][i_beg_group];
    } // End loop on partitions


    // > Exchange with part_to_block
    PDM_part_to_block_t *ptb_group = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                              PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                              1.,
                                                              _group_gnum,
                                                              NULL,
                                                              n_elt_group,
                                                              isos->iso_n_part, // will fail if n_part>1 because of extract_part ?
                                                              isos->comm);

    int dn_elt_group = PDM_part_to_block_n_elt_block_get(ptb_group);
    _iso_dedge_group_idx[i_group+1] = _iso_dedge_group_idx[i_group] + dn_elt_group;

    PDM_realloc(_iso_dedge_group_gnum,
                _iso_dedge_group_gnum,
                _iso_dedge_group_idx [i_group+1],
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
    for (int i_part=0; i_part<isos->n_part; ++i_part) {
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
    

  if (debug==1) {
    int size_group_gnum = _iso_dedge_group_idx[isos->iso_n_edge_group[id_iso]];
    log_trace("iso_n_edge_group = %d\n", isos->iso_n_edge_group[id_iso]);
    PDM_log_trace_array_int (_iso_dedge_group_idx , isos->iso_n_edge_group[id_iso]+1, "iso_dedge_group_idx  ::");
    PDM_log_trace_array_long(_iso_dedge_group_gnum, size_group_gnum                 , "iso_dedge_group_gnum ::");
  }


  /**
   * Set result
   */
  isos->iso_dedge_group_idx [id_iso] = _iso_dedge_group_idx;
  isos->iso_dedge_group_gnum[id_iso] = _iso_dedge_group_gnum;
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
  int              *out_dn_elt,
  int             **out_delt_vtx_idx,
  PDM_g_num_t     **out_delt_vtx,
  int             **out_delt_parent_idx,
  PDM_g_num_t     **out_delt_parent_gnum
)
{
  int debug = 0;


  // > Set elt in block frame
  PDM_part_to_block_t *ptb_elt = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_MERGE,
                                                          1.,
                                                          elt_gnum,
                                                          NULL,
                                                          n_elt,
                                                          isos->iso_n_part, // will fail if n_part>1 because of extract_part ?
                                                          isos->comm);

  /**
   * Count how many received for block elements
   */
  int  dn_elt           = PDM_part_to_block_n_elt_block_get(ptb_elt);
  int *block_gnum_count = PDM_part_to_block_block_gnum_count_get(ptb_elt);
  if (debug==1) {
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
      _elt_vtx[i_part][i_elt] = isos->iso_entity_gnum[PDM_MESH_ENTITY_VTX][id_iso][i_part][elt_lnum-1];
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

  if (debug) {
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
  PDM_malloc(_elt_parent_strd, isos->iso_n_part, int         *);
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

  if (debug) {
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


  /**
   * Go through received data :
   *   - select one connectivity example (should be the same on all procs)
   *   - unique parents (must have various from procs)
   * For ngon algo one data should be received for each entity
   */
  // int i_read = 0;
  // int i_read_conn = 0;
  // int i_write = 0;
  // int i_write_conn = 0;
  for (int i_elt=0; i_elt<dn_elt; ++i_elt) {
    int n_src = block_gnum_count[i_elt];
    if (_is_nodal(isos)) {
      assert(n_src==1);
    }

    if (n_src==1) {
      continue;
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "Not implemented yet.\n");
      // // > Go through each received data
      // for (int i_src=0; i_src<n_src; ++i_src) {
        
      //   // > Copy first connectivity received
      //   if (i_src==0) {
      //     for (int i_strd=0; i_strd<_delt_vtx_strd[i_read]; ++i_strd)
      //       _delt_vtx[i_write_conn] = _delt_vtx[i_read_conn]
      //       i_write_conn++;
      //       i_read_conn++;
      //     }
      //   }
      //   else {
      //     i_read_conn+=_delt_vtx_strd[i_read];
      //   }

      //   // > Count parent

      //   // > Unique parent

      //   // > Fill with unique parent

      // }
    }

  }
  // _delt_vtx_idx = PDM_array_new_idx_from_sizes_int(_delt_vtx_strd, dn_elt);

  *out_delt_parent_idx  = PDM_array_new_idx_from_sizes_int(delt_parent_strd, dn_elt);
  *out_delt_parent_gnum = delt_parent_gnum;
  PDM_free(delt_parent_strd);

  PDM_part_to_block_free(ptb_elt);
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

  /** Vertices */
  _part_to_dist_vtx(isos, id_iso);
  isos->iso_owner_dvtx_coord                             [id_iso] = PDM_OWNERSHIP_KEEP;
  isos->iso_owner_dvtx_parent_weight                     [id_iso] = PDM_OWNERSHIP_KEEP;
  isos->iso_owner_dparent           [PDM_MESH_ENTITY_VTX][id_iso] = PDM_OWNERSHIP_KEEP;


  /** Edges */
  _part_to_dist_elt(isos,
                    id_iso,
                    isos->iso_n_entity[PDM_MESH_ENTITY_EDGE][id_iso],
                    isos->iso_entity_gnum[PDM_MESH_ENTITY_EDGE][id_iso],
                    NULL,
                    isos->iso_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX][id_iso],
                    isos->iso_entity_parent_idx [PDM_MESH_ENTITY_EDGE][id_iso],
                    isos->iso_entity_parent_gnum[PDM_MESH_ENTITY_EDGE][id_iso],
                   &isos->iso_dn_entity[PDM_MESH_ENTITY_EDGE][id_iso],
                    NULL,
                   &isos->iso_dconnec[PDM_CONNECTIVITY_TYPE_EDGE_VTX][id_iso],
                   &isos->iso_dentity_parent_idx [PDM_MESH_ENTITY_EDGE][id_iso],
                   &isos->iso_dentity_parent_gnum[PDM_MESH_ENTITY_EDGE][id_iso]);
  isos->iso_owner_dconnec[PDM_CONNECTIVITY_TYPE_EDGE_VTX][id_iso] = PDM_OWNERSHIP_KEEP;
  isos->iso_owner_dparent[PDM_MESH_ENTITY_EDGE          ][id_iso] = PDM_OWNERSHIP_KEEP;

  isos->iso_owner_dedge_bnd[id_iso] = PDM_OWNERSHIP_KEEP;


  if (isos->entry_mesh_dim==3) {
    /** Edge groups */
    _part_to_dist_edge_group(isos, id_iso);


    /** Faces */
    _part_to_dist_elt(isos,
                      id_iso,
                      isos->iso_n_entity   [PDM_MESH_ENTITY_FACE          ][id_iso],
                      isos->iso_entity_gnum[PDM_MESH_ENTITY_FACE          ][id_iso],
                      isos->iso_connec_idx [PDM_CONNECTIVITY_TYPE_FACE_VTX][id_iso],
                      isos->iso_connec     [PDM_CONNECTIVITY_TYPE_FACE_VTX][id_iso],
                      isos->iso_entity_parent_idx [PDM_MESH_ENTITY_FACE][id_iso],
                      isos->iso_entity_parent_gnum[PDM_MESH_ENTITY_FACE][id_iso],
                     &isos->iso_dn_entity[PDM_MESH_ENTITY_FACE][id_iso],
                     &isos->iso_dconnec_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX][id_iso],
                     &isos->iso_dconnec    [PDM_CONNECTIVITY_TYPE_FACE_VTX][id_iso],
                     &isos->iso_dentity_parent_idx [PDM_MESH_ENTITY_FACE][id_iso],
                     &isos->iso_dentity_parent_gnum[PDM_MESH_ENTITY_FACE][id_iso]);
    isos->iso_owner_dconnec[PDM_CONNECTIVITY_TYPE_FACE_VTX][id_iso] = PDM_OWNERSHIP_KEEP;
    isos->iso_owner_dparent[PDM_MESH_ENTITY_FACE          ][id_iso] = PDM_OWNERSHIP_KEEP;
  }

}


static void
_get_gnums_from_pmne
(
  PDM_part_mesh_nodal_elmts_t *pmne,
  int                          i_part,
  PDM_g_num_t                 *gnum
)
{
  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  int idx = 0;
  for (int i_section = 0; i_section < n_section; i_section++) {

    int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                            sections_id[i_section],
                                                            i_part);

    int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                               sections_id[i_section],
                                                               i_part,
                                                               PDM_OWNERSHIP_BAD_VALUE);

    PDM_g_num_t *elt_gnum = PDM_part_mesh_nodal_elmts_g_num_get(pmne,
                                                                sections_id[i_section],
                                                                i_part,
                                                                PDM_OWNERSHIP_BAD_VALUE);

    for (int i = 0; i < n_elt; i++) {
      if (parent_num == NULL) {
        idx++;
      }
      else {
        idx = parent_num[i];
      }

      gnum[idx] = elt_gnum[i];
    }
  }
}


/* Build part_to_part to link isosurface entities with their 'parent' source entities in user frame */
static void
_build_ptp
(
  PDM_isosurface_t    *isos,
  int                  id_iso,
  PDM_mesh_entities_t  entity_type
)
{
  /**
   * TODO:
   * - retrieve the init location of extracted parents to build ptp from triplets
   *   and compose child -> extracted_parent -> init_parent
   *   => we already used {extracted_parent -> init_parent} in extract_part, is there a way
   *      to get it again without having to rely on a ptp_reverse_issend?
   */

  if (isos->compute_ptp[entity_type][id_iso] == PDM_FALSE) {
    return;
  }

  if (entity_type == PDM_MESH_ENTITY_FACE && isos->entry_mesh_dim < 3) {
    return;
  }

  assert(isos->is_dist_or_part == 1);
  assert(isos->extract_kind != PDM_EXTRACT_PART_KIND_LOCAL);

  PDM_extract_part_t *extrp = isos->extrp[id_iso];
  assert(extrp != NULL);

  int          *n_entity           = isos->iso_n_entity          [entity_type][id_iso];
  PDM_g_num_t **entity_gnum        = isos->iso_entity_gnum       [entity_type][id_iso];
  int         **entity_parent_idx  = isos->iso_entity_parent_idx [entity_type][id_iso];
  int         **entity_parent_lnum = isos->iso_entity_parent_lnum[entity_type][id_iso];
  PDM_g_num_t **entity_parent_gnum = isos->iso_entity_parent_gnum[entity_type][id_iso];
  int          *n_parent           = NULL;
  PDM_g_num_t **parent_gnum        = NULL;

  if (_is_nodal(isos)) {
    PDM_malloc(n_parent,    isos->n_part, int          );
    PDM_malloc(parent_gnum, isos->n_part, PDM_g_num_t *);
  }

  switch (entity_type) {
    case PDM_MESH_ENTITY_VTX: {
      if (_is_nodal(isos)) {
        for (int i_part = 0; i_part < isos->n_part; i_part++) {
          n_parent   [i_part] = PDM_part_mesh_nodal_n_vtx_get    (isos->pmesh_nodal, i_part);
          parent_gnum[i_part] = PDM_part_mesh_nodal_vtx_g_num_get(isos->pmesh_nodal, i_part);
        }
      }
      else {
        n_parent    = isos->n_vtx;
        parent_gnum = isos->vtx_gnum;
      }
      break;
    }

    case PDM_MESH_ENTITY_EDGE: {
      if (_is_nodal(isos)) {
        PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(isos->pmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);
        for (int i_part = 0; i_part < isos->n_part; i_part++) {
          n_parent[i_part] = PDM_part_mesh_nodal_n_elmts_get(isos->pmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, i_part);
          PDM_malloc(parent_gnum[i_part], n_parent[i_part], PDM_g_num_t);
          _get_gnums_from_pmne(pmne, i_part, parent_gnum[i_part]);
        }
      }
      else {
        n_parent    = isos->n_face;
        parent_gnum = isos->face_gnum;
      }
      break;
    }

    case PDM_MESH_ENTITY_FACE: {
      if (_is_nodal(isos)) {
        PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(isos->pmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC);
        for (int i_part = 0; i_part < isos->n_part; i_part++) {
          n_parent[i_part] = PDM_part_mesh_nodal_n_elmts_get(isos->pmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, i_part);
          PDM_malloc(parent_gnum[i_part], n_parent[i_part], PDM_g_num_t);
          _get_gnums_from_pmne(pmne, i_part, parent_gnum[i_part]);
        }
      }
      else {
        n_parent    = isos->n_cell;
        parent_gnum = isos->cell_gnum;
      }
      break;
    }

    default : {
      PDM_error(__FILE__, __LINE__, 0, "Invalid entity_type %d\n", entity_type);
    }
  }

  PDM_part_to_part_t *ptp = NULL;

  // ptp = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) extrp->target_gnum,
  //                                                 (const int         * ) extrp->n_target,
  //                                                                        isos->iso_n_part,
  //                                                 (const int         * ) n_parent,
  //                                                                        isos->n_part,
  //                                                 (const int         **) entity_parent_idx,
  //                                                                        NULL,
  //                                                 (const int         **) entity_parent_init_location,
  //                                                                        extrp->comm);

  // log_trace("id_iso = %d, entity_type %d\n", id_iso, entity_type);
  // for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {
  //   for (int i = 0; i < n_entity[i_part]; i++) {
  //     log_trace(PDM_FMT_G_NUM" : ", entity_gnum[i_part][i]);
  //     PDM_log_trace_array_long(entity_parent_gnum[i_part] + entity_parent_idx[i_part][i],
  //                              entity_parent_idx[i_part][i+1] - entity_parent_idx[i_part][i],
  //                              "");
  //   }
  // }


  ptp = PDM_part_to_part_create((const PDM_g_num_t **) entity_gnum,
                                (const int         * ) n_entity,
                                                       isos->iso_n_part,
                                (const PDM_g_num_t **) parent_gnum,
                                (const int         * ) n_parent,
                                                       isos->n_part,
                                (const int         **) entity_parent_idx,
                                (const PDM_g_num_t **) entity_parent_gnum,
                                                       isos->comm);

  isos->iso_owner_ptp[entity_type][id_iso] = PDM_OWNERSHIP_KEEP;
  isos->iso_ptp[entity_type][id_iso] = ptp;


  if (_is_nodal(isos)) {
    PDM_free(n_parent); // no worries, it is deep-copied in part_to_part creation
  }
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

  /* Vertices */
  if (entity_type == PDM_MESH_ENTITY_VTX) {
    /* Partitioned */
    for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {
      if (isos->iso_owner_vtx_coord[id_iso][i_part] == PDM_OWNERSHIP_KEEP) {
        PDM_free(isos->iso_vtx_coord[id_iso][i_part]);
      }
      if (isos->iso_owner_vtx_parent_weight[id_iso][i_part] == PDM_OWNERSHIP_KEEP) {
        PDM_free(isos->iso_vtx_parent_weight[id_iso][i_part]);
      }
    }
    PDM_free(isos->iso_vtx_coord        [id_iso]);
    PDM_free(isos->iso_vtx_parent_weight[id_iso]);

    /* Block-distributed */
    if (isos->iso_owner_dvtx_coord[id_iso] == PDM_OWNERSHIP_KEEP) {
      PDM_free(isos->iso_dvtx_coord[id_iso]);
    }
    if (isos->iso_owner_dvtx_parent_weight[id_iso] == PDM_OWNERSHIP_KEEP) {
      PDM_free(isos->iso_dvtx_parent_weight[id_iso]);
    }
  }

  /* Edges */
  else if (entity_type == PDM_MESH_ENTITY_EDGE) {
    /* Partitioned */
    for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {
      if (isos->entry_mesh_dim == 3 &&
          isos->iso_owner_edge_bnd[id_iso][i_part] == PDM_OWNERSHIP_KEEP) {
        PDM_free(isos->iso_edge_group_idx [id_iso][i_part]);
        PDM_free(isos->iso_edge_group_lnum[id_iso][i_part]);
        PDM_free(isos->iso_edge_group_gnum[id_iso][i_part]);
      }

      if (isos->iso_owner_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX][id_iso][i_part] == PDM_OWNERSHIP_KEEP) {
        PDM_free(isos->iso_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX][id_iso][i_part]);
      }
    }
    PDM_free(isos->iso_edge_group_idx [id_iso]);
    PDM_free(isos->iso_edge_group_lnum[id_iso]);
    PDM_free(isos->iso_edge_group_gnum[id_iso]);
    PDM_free(isos->iso_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX][id_iso]);

    /* Block-distributed */
    if (isos->entry_mesh_dim == 3 &&
        isos->iso_owner_dedge_bnd[id_iso] == PDM_OWNERSHIP_KEEP) {
      PDM_free(isos->iso_dedge_group_idx [id_iso]);
      PDM_free(isos->iso_dedge_group_gnum[id_iso]);
    }

    if (isos->iso_owner_dconnec[PDM_CONNECTIVITY_TYPE_EDGE_VTX][id_iso] == PDM_OWNERSHIP_KEEP) {
      PDM_free(isos->iso_dconnec[PDM_CONNECTIVITY_TYPE_EDGE_VTX][id_iso]);
    }
  }

  /* Faces */
  else if (entity_type == PDM_MESH_ENTITY_FACE) {
    /* Partitioned */
    for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {
      if (isos->iso_owner_connec[PDM_CONNECTIVITY_TYPE_FACE_VTX][id_iso][i_part] == PDM_OWNERSHIP_KEEP) {
        PDM_free(isos->iso_connec_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX][id_iso][i_part]);
        PDM_free(isos->iso_connec    [PDM_CONNECTIVITY_TYPE_FACE_VTX][id_iso][i_part]);
      }
    }
    PDM_free(isos->iso_connec_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX][id_iso]);
    PDM_free(isos->iso_connec    [PDM_CONNECTIVITY_TYPE_FACE_VTX][id_iso]);

    /* Block-distributed */
    if (isos->iso_owner_dconnec[PDM_CONNECTIVITY_TYPE_FACE_VTX][id_iso] == PDM_OWNERSHIP_KEEP) {
      PDM_free(isos->iso_dconnec_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX][id_iso]);
      PDM_free(isos->iso_dconnec    [PDM_CONNECTIVITY_TYPE_FACE_VTX][id_iso]);
    }
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Invalid entity_type %d\n", entity_type);
  }


  /* Generic for all entities */

  /*   Partitioned */
  for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {
    if (isos->iso_owner_gnum[entity_type][id_iso][i_part] == PDM_OWNERSHIP_KEEP) {
      PDM_free(isos->iso_entity_gnum[entity_type][id_iso][i_part]);
    }
    if (isos->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL &&
        isos->iso_owner_parent_lnum[entity_type][id_iso][i_part] == PDM_OWNERSHIP_KEEP) {
      PDM_free(isos->iso_entity_parent_idx [entity_type][id_iso][i_part]);
      PDM_free(isos->iso_entity_parent_lnum[entity_type][id_iso][i_part]);
    }
    if (isos->extract_kind==PDM_EXTRACT_PART_KIND_REEQUILIBRATE) {
      PDM_free(isos->iso_entity_parent_idx [entity_type][id_iso][i_part]);
      PDM_free(isos->iso_entity_parent_gnum[entity_type][id_iso][i_part]);
    }

    if (isos->iso_owner_isovalue_entity_idx[entity_type][id_iso][i_part] == PDM_OWNERSHIP_KEEP) {
      PDM_free(isos->isovalue_entity_idx[entity_type][id_iso][i_part]);
    }
  }
  PDM_free(isos->iso_n_entity          [entity_type][id_iso]);
  PDM_free(isos->iso_entity_gnum       [entity_type][id_iso]);
  PDM_free(isos->iso_entity_parent_idx [entity_type][id_iso]);
  PDM_free(isos->iso_entity_parent_lnum[entity_type][id_iso]);
  PDM_free(isos->iso_entity_parent_gnum[entity_type][id_iso]);
  PDM_free(isos->isovalue_entity_idx   [entity_type][id_iso]);

  if (isos->  compute_ptp[entity_type][id_iso] == PDM_TRUE &&
      isos->iso_owner_ptp[entity_type][id_iso] == PDM_OWNERSHIP_KEEP) {
    PDM_part_to_part_free(isos->iso_ptp[entity_type][id_iso]);
  }


  /*   Block-distributed */
  if (isos->iso_owner_dparent[entity_type][id_iso] == PDM_OWNERSHIP_KEEP) {
    PDM_free(isos->iso_dentity_parent_idx [entity_type][id_iso]);
    PDM_free(isos->iso_dentity_parent_gnum[entity_type][id_iso]);
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
  if (isos->field[id_iso]!=NULL) {
    for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
      if (isos->kind[id_iso]!=PDM_ISO_SURFACE_KIND_FIELD ||
          isos->is_dist_or_part==0) {
        PDM_free(isos->field[id_iso][i_part]);
      }
    }
    if (!partial) {
      PDM_free(isos->field[id_iso]);
    }
  }

  if (isos->extract_field[id_iso] != NULL) {
    for (int i_part = 0; i_part < isos->iso_n_part; i_part++) {
      PDM_free(isos->extract_field[id_iso][i_part]);
    }
    PDM_free(isos->extract_field[id_iso]);
  }


}


static void
_reset_downer
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  isos->iso_owner_dvtx_coord        [id_iso] = PDM_OWNERSHIP_BAD_VALUE;
  isos->iso_owner_dvtx_parent_weight[id_iso] = PDM_OWNERSHIP_BAD_VALUE;
  isos->iso_owner_dedge_bnd         [id_iso] = PDM_OWNERSHIP_BAD_VALUE;

  for (int i_entity=0; i_entity<PDM_MESH_ENTITY_MAX; ++i_entity) {
    isos->iso_owner_dparent[i_entity][id_iso] = PDM_OWNERSHIP_BAD_VALUE;
  }
  for (int i_connec=0; i_connec<PDM_CONNECTIVITY_TYPE_MAX; ++i_connec) {
    isos->iso_owner_dconnec[i_connec][id_iso] = PDM_OWNERSHIP_BAD_VALUE;
  }

}


static void
_free_owner
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  if (isos->is_computed[id_iso] == PDM_FALSE) {
    return;
  }

  // > Partitionned
  PDM_free(isos->iso_owner_vtx_coord        [id_iso]);
  PDM_free(isos->iso_owner_vtx_parent_weight[id_iso]);
  PDM_free(isos->iso_owner_edge_bnd         [id_iso]);
  for (int i_entity=0; i_entity<PDM_MESH_ENTITY_MAX; ++i_entity) {
    PDM_free(isos->iso_owner_gnum               [i_entity][id_iso]);
    PDM_free(isos->iso_owner_parent_lnum        [i_entity][id_iso]);
    PDM_free(isos->iso_owner_isovalue_entity_idx[i_entity][id_iso]);
  }
  PDM_free(isos->iso_owner_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX][id_iso]);
  PDM_free(isos->iso_owner_connec[PDM_CONNECTIVITY_TYPE_FACE_VTX][id_iso]);
}


static void
_free_nodal_extract_parts
(
  PDM_isosurface_t *isos
)
{
  for (int i_part=0; i_part<isos->n_part; ++i_part) {
    PDM_free(isos->extract_vtx_lnum[i_part]);
    PDM_free(isos->extract_tri_lnum[i_part]);
    PDM_free(isos->extract_tri_tag [i_part]);
    if (isos->entry_mesh_dim==3) {
      PDM_free(isos->extract_tetra_lnum[i_part]);
    }
  }
}



static void
_isosurface_reset
(
 PDM_isosurface_t *isos,
 int               id_isosurface
)
{
  if (isos->is_computed[id_isosurface] == PDM_FALSE) {
    return;
  }

  _free_iso_entity(isos, PDM_MESH_ENTITY_VTX,  id_isosurface);
  _free_iso_entity(isos, PDM_MESH_ENTITY_EDGE, id_isosurface);
  _free_iso_entity(isos, PDM_MESH_ENTITY_FACE, id_isosurface);
  _free_field     (isos, id_isosurface, isos->is_dist_or_part);
  // > Distributed
  if (isos->is_dist_or_part==0) {
    _reset_downer(isos, id_isosurface);
  }
  _free_owner(isos, id_isosurface);

  PDM_extract_part_free(isos->extrp[id_isosurface]);

  isos->is_computed[id_isosurface] = PDM_FALSE;
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

  int debug = 0;

  if (debug==1) {
    log_trace("PDM_isosurface:: compute isosurface n°%d\n", id_isosurface);
  }

  if (isos->is_computed[id_isosurface] == PDM_TRUE) {
    // TODO: Warning or Error?
    if (i_rank == 0) {
      printf("WARNING - PDM_isosurface_compute : id_isosurface %d will be reset before being re-computed\n", id_isosurface);
    }
    _isosurface_reset(isos, id_isosurface);
  }

  assert(isos->is_computed[id_isosurface] == PDM_FALSE);

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

  // Comment because not implemented yet and i need to test full tetra/tri
  // /* If nodal with sections other than TRIA3 and TETRA4, fall back to ngon */
  // if (_is_nodal(isos)) {
  //   _ngonize(isos, id_isosurface);
  // }

  /* Build isosurface mesh */
  if (_is_nodal(isos)) {
    PDM_isosurface_marching_algo(isos,
                                 id_isosurface);
  }
  else {
    PDM_isosurface_ngon_algo(isos,
                             id_isosurface);
  }

  if (isos->is_dist_or_part == 0) {
    /* Block-distribute the isosurface */
    _part_to_dist(isos, id_isosurface);
  }
  else {
    // Partitioned
    if (isos->extract_kind != PDM_EXTRACT_PART_KIND_LOCAL) {
      for (int i_entity = 0; i_entity < PDM_MESH_ENTITY_MAX; i_entity++) {
        _build_ptp(isos,
                   id_isosurface,
                   i_entity);
      }
    }
  }

  // > Free mesh extraction arrays for nodal
  if (_is_nodal(isos)) {
    _free_nodal_extract_parts(isos);
  }

  isos->is_computed[id_isosurface] = PDM_TRUE;
}


static
void
_assert_isovalues_not_to_close
(
  PDM_isosurface_t *isos,
  int               id_isosurface,
  int               n_isovalues,
  double           *isovalues
)
{
  // > Check that difference between isovalues>ISOSURFACE_EPS
  for (int i_iso=0; i_iso<n_isovalues; ++i_iso) {
    for (int j_iso=i_iso+1; j_iso<n_isovalues; ++j_iso) {
      if (PDM_ABS(isovalues[i_iso]-isovalues[j_iso])<=isos->ISOSURFACE_EPS) {
        PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: isovalue %d = %f too close from isovalue %d = %f (%.2e<=%.2e) for isosurface with id %d.\n",
                                                        i_iso, isovalues[i_iso],
                                                        j_iso, isovalues[j_iso],
                                                        PDM_ABS(isovalues[i_iso]-isovalues[j_iso]), isos->ISOSURFACE_EPS,
                                                        id_isosurface);
      }
    }
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
  memset(isos, 0, sizeof(PDM_isosurface_t));

  // > Save communicator
  isos->comm = comm;

  // > Init tolerance
  isos->ISOSURFACE_EPS = 0.;

  // > Entry mesh information
  isos->is_dist_or_part = -1; 
  isos->entry_mesh_type =  0; 
  isos->entry_mesh_dim  =  mesh_dimension;

  // // > Isosurface mesh information
  // isos->iso_elt_type = elt_type; 
  isos->extract_kind = PDM_EXTRACT_PART_KIND_LOCAL; 
  isos->n_part = 1;

  isos->we_have_edges = -1;

  return isos;
}


void
PDM_isosurface_set_tolerance
(
  PDM_isosurface_t *isos,
  double            tolerance
)
{
  if (tolerance < 0.) {
    PDM_error(__FILE__, __LINE__, 0, "Invalid tolerance : %f (must be >= 0\n", tolerance);
  }

  isos->ISOSURFACE_EPS = tolerance;
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
  
  _assert_isovalues_not_to_close(isos, id_isosurface, n_isovalues, isovalues);

  isos->n_isovalues[id_isosurface] = n_isovalues;
  PDM_realloc(isos->isovalues[id_isosurface], isos->isovalues[id_isosurface], n_isovalues, double);
  for (int i=0; i<n_isovalues; ++i) {
    isos->isovalues[id_isosurface][i] = isovalues[i];
  }
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

  PDM_realloc(isos->is_computed, isos->is_computed, isos->n_isosurface, PDM_bool_t);
  isos->is_computed[id_isosurface] = PDM_FALSE;
  
  _assert_isovalues_not_to_close(isos, id_isosurface, n_isovalues, isovalues);

  PDM_realloc(isos->kind          , isos->kind          , isos->n_isosurface, PDM_iso_surface_kind_t          );
  PDM_realloc(isos->eq_coeffs     , isos->eq_coeffs     , isos->n_isosurface, double                         *);
  PDM_realloc(isos->field_function, isos->field_function, isos->n_isosurface, PDM_isosurface_field_function_t );
  
  PDM_realloc(isos->field_function, isos->field_function, isos->n_isosurface, PDM_isosurface_field_function_t);
  
  PDM_realloc(isos->dfield, isos->dfield, isos->n_isosurface, double  *);
  PDM_realloc(isos->field , isos->field , isos->n_isosurface, double **);
  isos->field[id_isosurface] = NULL;

  PDM_realloc(isos->extract_field, isos->extract_field, isos->n_isosurface, double **);
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
    memset(isos->field[id_isosurface], 0, n_part*sizeof(double *));
  }

  PDM_realloc(isos->n_isovalues , isos->n_isovalues , isos->n_isosurface, int     );
  PDM_realloc(isos->isovalues   , isos->isovalues   , isos->n_isosurface, double *);
  PDM_realloc(isos->use_gradient, isos->use_gradient, isos->n_isosurface, int     );

  PDM_realloc(isos->extrp, isos->extrp, isos->n_isosurface, PDM_extract_part_t *);
  isos->extrp[id_isosurface] = NULL;



  for (int i = 0; i < PDM_MESH_ENTITY_MAX; i++) {
    PDM_realloc(isos->iso_n_entity          [i], isos->iso_n_entity          [i], isos->n_isosurface, int          *);
    PDM_realloc(isos->iso_entity_gnum       [i], isos->iso_entity_gnum       [i], isos->n_isosurface, PDM_g_num_t **);
    PDM_realloc(isos->iso_entity_parent_idx [i], isos->iso_entity_parent_idx [i], isos->n_isosurface, int         **);
    PDM_realloc(isos->iso_entity_parent_lnum[i], isos->iso_entity_parent_lnum[i], isos->n_isosurface, int         **);
    PDM_realloc(isos->iso_entity_parent_gnum[i], isos->iso_entity_parent_gnum[i], isos->n_isosurface, PDM_g_num_t **);
    PDM_realloc(isos->isovalue_entity_idx   [i], isos->isovalue_entity_idx   [i], isos->n_isosurface, int         **);
    PDM_realloc(isos->compute_ptp           [i], isos->compute_ptp           [i], isos->n_isosurface, PDM_bool_t    );
    isos->iso_n_entity          [i][id_isosurface] = NULL;
    isos->iso_entity_gnum       [i][id_isosurface] = NULL;
    isos->iso_entity_parent_idx [i][id_isosurface] = NULL;
    isos->iso_entity_parent_lnum[i][id_isosurface] = NULL;
    isos->iso_entity_parent_gnum[i][id_isosurface] = NULL;
    isos->isovalue_entity_idx   [i][id_isosurface] = NULL;
    isos->compute_ptp           [i][id_isosurface] = PDM_FALSE;
  }

  // > Vertices
  PDM_realloc(isos->iso_vtx_coord        , isos->iso_vtx_coord        , isos->n_isosurface, double **);
  PDM_realloc(isos->iso_vtx_parent_weight, isos->iso_vtx_parent_weight, isos->n_isosurface, double **);
  isos->iso_vtx_coord        [id_isosurface] = NULL;
  isos->iso_vtx_parent_weight[id_isosurface] = NULL;

  // > Edges
  PDM_realloc(isos->iso_n_edge_group   , isos->iso_n_edge_group   , isos->n_isosurface, int            );
  PDM_realloc(isos->iso_edge_group_idx , isos->iso_edge_group_idx , isos->n_isosurface, int          **);
  PDM_realloc(isos->iso_edge_group_lnum, isos->iso_edge_group_lnum, isos->n_isosurface, int          **);
  PDM_realloc(isos->iso_edge_group_gnum, isos->iso_edge_group_gnum, isos->n_isosurface, PDM_g_num_t  **);
  isos->iso_n_edge_group   [id_isosurface] = 0;
  isos->iso_edge_group_idx [id_isosurface] = NULL;
  isos->iso_edge_group_lnum[id_isosurface] = NULL;
  isos->iso_edge_group_gnum[id_isosurface] = NULL;

  PDM_realloc(isos->iso_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX],
              isos->iso_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX], isos->n_isosurface, int **);
  isos->iso_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX][id_isosurface] = NULL;

  // > Faces
  if (isos->entry_mesh_dim == 3) {
    PDM_realloc(isos->iso_connec_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX],
                isos->iso_connec_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX], isos->n_isosurface, int **);
    PDM_realloc(isos->iso_connec    [PDM_CONNECTIVITY_TYPE_FACE_VTX],
                isos->iso_connec    [PDM_CONNECTIVITY_TYPE_FACE_VTX], isos->n_isosurface, int **);
    isos->iso_connec_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX][id_isosurface] = NULL;
    isos->iso_connec    [PDM_CONNECTIVITY_TYPE_FACE_VTX][id_isosurface] = NULL;
  }

  // > Part_to_part between iso entities and entry mesh entites
  PDM_realloc(isos->iso_ptp[PDM_MESH_ENTITY_VTX ], isos->iso_ptp[PDM_MESH_ENTITY_VTX ], isos->n_isosurface, PDM_part_to_part_t *);
  PDM_realloc(isos->iso_ptp[PDM_MESH_ENTITY_EDGE], isos->iso_ptp[PDM_MESH_ENTITY_EDGE], isos->n_isosurface, PDM_part_to_part_t *);
  if (isos->entry_mesh_dim == 3) {
    PDM_realloc(isos->iso_ptp[PDM_MESH_ENTITY_FACE], isos->iso_ptp[PDM_MESH_ENTITY_FACE], isos->n_isosurface, PDM_part_to_part_t *);
  }

  // > Partitioned owners
  PDM_realloc(isos->iso_owner_vtx_coord        , isos->iso_owner_vtx_coord        , isos->n_isosurface, PDM_ownership_t  *);
  PDM_realloc(isos->iso_owner_vtx_parent_weight, isos->iso_owner_vtx_parent_weight, isos->n_isosurface, PDM_ownership_t  *);
  PDM_realloc(isos->iso_owner_edge_bnd         , isos->iso_owner_edge_bnd         , isos->n_isosurface, PDM_ownership_t  *);
  isos->iso_owner_vtx_coord        [id_isosurface] = NULL;
  isos->iso_owner_vtx_parent_weight[id_isosurface] = NULL;
  isos->iso_owner_edge_bnd         [id_isosurface] = NULL;

  for (int i_entity = 0; i_entity < PDM_MESH_ENTITY_MAX; i_entity++) {
    PDM_realloc(isos->iso_owner_gnum               [i_entity], isos->iso_owner_gnum               [i_entity], isos->n_isosurface, PDM_ownership_t *);
    PDM_realloc(isos->iso_owner_parent_lnum        [i_entity], isos->iso_owner_parent_lnum        [i_entity], isos->n_isosurface, PDM_ownership_t *);
    PDM_realloc(isos->iso_owner_isovalue_entity_idx[i_entity], isos->iso_owner_isovalue_entity_idx[i_entity], isos->n_isosurface, PDM_ownership_t *);
    PDM_realloc(isos->iso_owner_ptp                [i_entity], isos->iso_owner_ptp                [i_entity], isos->n_isosurface, PDM_ownership_t);

    isos->iso_owner_gnum               [i_entity][id_isosurface] = NULL;
    isos->iso_owner_parent_lnum        [i_entity][id_isosurface] = NULL;
    isos->iso_owner_isovalue_entity_idx[i_entity][id_isosurface] = NULL;
    isos->iso_owner_ptp                [i_entity][id_isosurface] = PDM_OWNERSHIP_BAD_VALUE;
  }

  PDM_realloc(isos->iso_owner_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX], isos->iso_owner_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX], isos->n_isosurface, PDM_ownership_t *);
  PDM_realloc(isos->iso_owner_connec[PDM_CONNECTIVITY_TYPE_FACE_VTX], isos->iso_owner_connec[PDM_CONNECTIVITY_TYPE_FACE_VTX], isos->n_isosurface, PDM_ownership_t *);
  isos->iso_owner_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX][id_isosurface] = NULL;
  isos->iso_owner_connec[PDM_CONNECTIVITY_TYPE_FACE_VTX][id_isosurface] = NULL;




  for (int i = 0; i < PDM_MESH_ENTITY_MAX; i++) {
    PDM_realloc(isos->iso_dn_entity          [i], isos->iso_dn_entity          [i], isos->n_isosurface, int           );
    PDM_realloc(isos->iso_dentity_parent_idx [i], isos->iso_dentity_parent_idx [i], isos->n_isosurface, int          *);
    PDM_realloc(isos->iso_dentity_parent_gnum[i], isos->iso_dentity_parent_gnum[i], isos->n_isosurface, PDM_g_num_t  *);
    isos->iso_dn_entity          [i][id_isosurface] = 0;
    isos->iso_dentity_parent_idx [i][id_isosurface] = NULL;
    isos->iso_dentity_parent_gnum[i][id_isosurface] = NULL;
  }

  // > Distributed iso vertices
  PDM_realloc(isos->iso_dvtx_coord        , isos->iso_dvtx_coord        , isos->n_isosurface, double       *);
  PDM_realloc(isos->iso_dvtx_parent_weight, isos->iso_dvtx_parent_weight, isos->n_isosurface, double       *);
  isos->iso_dvtx_coord        [id_isosurface] = NULL;
  isos->iso_dvtx_parent_weight[id_isosurface] = NULL;

  // > Distributed iso edges
  PDM_realloc(isos->iso_dedge_group_idx , isos->iso_dedge_group_idx , isos->n_isosurface, int          *);
  PDM_realloc(isos->iso_dedge_group_gnum, isos->iso_dedge_group_gnum, isos->n_isosurface, PDM_g_num_t  *);
  isos->iso_dedge_group_idx [id_isosurface] = NULL;
  isos->iso_dedge_group_gnum[id_isosurface] = NULL;

  PDM_realloc(isos->iso_dconnec[PDM_CONNECTIVITY_TYPE_EDGE_VTX], isos->iso_dconnec[PDM_CONNECTIVITY_TYPE_EDGE_VTX], isos->n_isosurface, PDM_g_num_t  *);
  isos->iso_dconnec[PDM_CONNECTIVITY_TYPE_EDGE_VTX][id_isosurface] = NULL;

  // > Distributed iso faces
  PDM_realloc(isos->iso_dconnec_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX], isos->iso_dconnec_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX], isos->n_isosurface, int          *);
  PDM_realloc(isos->iso_dconnec    [PDM_CONNECTIVITY_TYPE_FACE_VTX], isos->iso_dconnec    [PDM_CONNECTIVITY_TYPE_FACE_VTX], isos->n_isosurface, PDM_g_num_t  *);
  isos->iso_dconnec_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX][id_isosurface] = NULL;
  isos->iso_dconnec    [PDM_CONNECTIVITY_TYPE_FACE_VTX][id_isosurface] = NULL;

  // > Distributed owners
  PDM_realloc(isos->iso_owner_dvtx_coord        , isos->iso_owner_dvtx_coord        , isos->n_isosurface, PDM_ownership_t);
  PDM_realloc(isos->iso_owner_dvtx_parent_weight, isos->iso_owner_dvtx_parent_weight, isos->n_isosurface, PDM_ownership_t);
  PDM_realloc(isos->iso_owner_dedge_bnd         , isos->iso_owner_dedge_bnd         , isos->n_isosurface, PDM_ownership_t);
  isos->iso_owner_dvtx_coord        [id_isosurface] = PDM_OWNERSHIP_BAD_VALUE;
  isos->iso_owner_dvtx_parent_weight[id_isosurface] = PDM_OWNERSHIP_BAD_VALUE;
  isos->iso_owner_dedge_bnd         [id_isosurface] = PDM_OWNERSHIP_BAD_VALUE;

  for (int i_entity = 0; i_entity < PDM_MESH_ENTITY_MAX; i_entity++) {
    PDM_realloc(isos->iso_owner_dparent[i_entity], isos->iso_owner_dparent[i_entity], isos->n_isosurface, PDM_ownership_t);
    isos->iso_owner_dparent[i_entity][id_isosurface] = PDM_OWNERSHIP_BAD_VALUE;
  }

  for (int i_connectivity = 0; i_connectivity < PDM_CONNECTIVITY_TYPE_MAX; i_connectivity++) {
    PDM_realloc(isos->iso_owner_dconnec[i_connectivity], isos->iso_owner_dconnec[i_connectivity], isos->n_isosurface, PDM_ownership_t);
    isos->iso_owner_dconnec[i_connectivity][id_isosurface] = PDM_OWNERSHIP_BAD_VALUE;
  }


  isos->kind       [id_isosurface] = kind;
  isos->n_isovalues[id_isosurface] = n_isovalues;
  PDM_malloc(isos->isovalues[id_isosurface], n_isovalues, double);
  for (int i=0; i<n_isovalues; ++i) {
    isos->isovalues[id_isosurface][i] = isovalues[i];
  }

  int n_coeff = _iso_surface_kind_n_coeff(kind);
  PDM_malloc(isos->eq_coeffs[id_isosurface], n_coeff, double);

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
  PDM_ISOSURFACE_CHECK_ID(isos, id_isosurface);

  int n_coeff = _iso_surface_kind_n_coeff(isos->kind[id_isosurface]);

  for (int i_coeff=0; i_coeff<n_coeff; ++i_coeff) {
    isos->eq_coeffs[id_isosurface][i_coeff] = coeff[i_coeff];
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
  PDM_ISOSURFACE_CHECK_ID(isos, id_isosurface);

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
  PDM_ISOSURFACE_CHECK_ID(isos, id_isosurface);

  if (id_isosurface==-1) { // Reset all isosurface
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
  if (id_isosurface==-1) { // Compute all isosurface
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

    if (_is_nodal(isos)) {
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
  PDM_free(isos->n_group_face);
  PDM_free(isos->group_face_idx);
  PDM_free(isos->group_face);
  PDM_free(isos->group_face_gnum);

  // > Free mesh extraction arrays for nodal
  PDM_free(isos->extract_n_vtx    );
  PDM_free(isos->extract_vtx_coord);
  PDM_free(isos->extract_vtx_gnum );
  PDM_free(isos->extract_vtx_lnum );

  PDM_free(isos->extract_n_tri      );
  PDM_free(isos->extract_tri_vtx    );
  PDM_free(isos->extract_tri_gnum   );
  PDM_free(isos->extract_tri_lnum   );
  PDM_free(isos->extract_tri_n_group);
  PDM_free(isos->extract_tri_tag    );

  PDM_free(isos->extract_n_tetra   );
  PDM_free(isos->extract_tetra_vtx );
  PDM_free(isos->extract_tetra_gnum);
  PDM_free(isos->extract_tetra_lnum);


  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    if (isos->n_isovalues[id_iso]>0) {
      PDM_free(isos->isovalues[id_iso]);
      PDM_free(isos->eq_coeffs[id_iso]);
    }
  }

  PDM_free(isos->n_isovalues);
  PDM_free(isos->isovalues);
  PDM_free(isos->eq_coeffs);
  PDM_free(isos->use_gradient);
  PDM_free(isos->field_function);
  PDM_free(isos->dfield);
  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    _free_field(isos, id_iso, 0);
  }
  PDM_free(isos->field);
  PDM_free(isos->kind);


  /**
   * Partitionned iso
   */
  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    if (isos->is_computed[id_iso] == PDM_TRUE) {
      _free_iso_entity(isos, PDM_MESH_ENTITY_VTX , id_iso);
      _free_iso_entity(isos, PDM_MESH_ENTITY_EDGE, id_iso);
      _free_iso_entity(isos, PDM_MESH_ENTITY_FACE, id_iso);
      PDM_extract_part_free(isos->extrp[id_iso]);
    }
  }
  PDM_free(isos->extract_field);
  PDM_free(isos->extrp);


  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    _free_owner(isos, id_iso);
  }
  for (int i = 0; i < PDM_MESH_ENTITY_MAX; i++) {
    PDM_free(isos->iso_n_entity          [i]);
    PDM_free(isos->iso_entity_gnum       [i]);
    PDM_free(isos->iso_entity_parent_idx [i]);
    PDM_free(isos->iso_entity_parent_lnum[i]);
    PDM_free(isos->iso_entity_parent_gnum[i]);
    PDM_free(isos->isovalue_entity_idx   [i]);
  }
  for (int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; i++) {
    PDM_free(isos->iso_connec_idx  [i]);
    PDM_free(isos->iso_connec      [i]);
    PDM_free(isos->iso_owner_connec[i]);
  }
  PDM_free(isos->iso_vtx_coord);
  PDM_free(isos->iso_vtx_parent_weight);

  PDM_free(isos->iso_n_edge_group);
  PDM_free(isos->iso_edge_group_idx);
  PDM_free(isos->iso_edge_group_lnum);
  PDM_free(isos->iso_edge_group_gnum);


  // > Part_to_part between iso entities and entry mesh entites
  for (PDM_mesh_entities_t i_entity = 0; i_entity < PDM_MESH_ENTITY_MAX; i_entity++) {

    PDM_free(isos->iso_owner_gnum               [i_entity]);
    PDM_free(isos->iso_owner_parent_lnum        [i_entity]);
    PDM_free(isos->iso_owner_isovalue_entity_idx[i_entity]);

    PDM_free(isos->compute_ptp  [i_entity]);
    PDM_free(isos->iso_owner_ptp[i_entity]);
    PDM_free(isos->iso_ptp      [i_entity]);
  }


  /**
   * Distributed iso
   */
  for (int i = 0; i < PDM_MESH_ENTITY_MAX; i++) {
    PDM_free(isos->iso_dn_entity          [i]);
    PDM_free(isos->iso_dentity_parent_idx [i]);
    PDM_free(isos->iso_dentity_parent_gnum[i]);
    PDM_free(isos->iso_owner_dparent      [i]);
  }
  for (int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; i++) {
    PDM_free(isos->iso_dconnec_idx  [i]);
    PDM_free(isos->iso_dconnec      [i]);
    PDM_free(isos->iso_owner_dconnec[i]);
  }

  PDM_free(isos->iso_dvtx_coord);
  PDM_free(isos->iso_dvtx_parent_weight);

  PDM_free(isos->iso_dedge_group_idx);
  PDM_free(isos->iso_dedge_group_gnum);


  PDM_free(isos->iso_owner_vtx_coord);
  PDM_free(isos->iso_owner_vtx_parent_weight);
  PDM_free(isos->iso_owner_edge_bnd);

  PDM_free(isos->iso_owner_dvtx_coord);
  PDM_free(isos->iso_owner_dvtx_parent_weight);
  PDM_free(isos->iso_owner_dedge_bnd);

  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    PDM_free(isos->compute_ptp[id_iso]);
  }
  PDM_free(isos->compute_ptp);

  PDM_free(isos->is_computed);


  PDM_free(isos);
}


#ifdef  __cplusplus
}
#endif
