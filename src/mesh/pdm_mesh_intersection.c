/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_mesh_intersection_priv.h"
#include "pdm_mesh_intersection.h"
#include "pdm_timer.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_sort.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_binary_search.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_priv.h"
#include "pdm_vtk.h"


#ifdef __cplusplus
extern "C"
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

static
void
_compute_extents_3d
(
 int      n_cell,
 int     *cell_face_idx,
 int     *cell_face,
 int     *face_vtx_idx,
 int     *face_vtx,
 double  *vtx_coord,
 double  *box_extents,
 double  *global_extents
)
{
  const double tolerance   = 1.e-12;
  const double eps_extents = 1.e-7;
  const int dim = 3;

  /* Loop over cell */
  for(int i_cell = 0; i_cell < n_cell; ++i_cell ) {

    double *_extents = box_extents + 6 * i_cell;

    /* Init */
    for (int i_dim = 0; i_dim < 3; i_dim++) {
      _extents[  i_dim] =  HUGE_VAL;
      _extents[3+i_dim] = -HUGE_VAL;
    }

    /* Loop over face and vtx to compute bbox */
    for(int idx_face = cell_face_idx[i_cell]; idx_face < cell_face_idx[i_cell+1]; ++idx_face) {

      int i_face = PDM_ABS(cell_face[idx_face])-1;

      for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
        int i_vtx = face_vtx[idx_vtx]-1;

        for (int i_dim = 0; i_dim < 3; i_dim++) {
          double x = vtx_coord[3*i_vtx + i_dim];

          if (x < _extents[i_dim]) {
            _extents[i_dim] = x;
          }
          if (x > _extents[3+i_dim]) {
            _extents[3+i_dim] = x;
          }
        }
      }
    }

    double delta = 0.;
    for (int i_dim = 0; i_dim < 3; i_dim++) {
      double x = _extents[3+i_dim] - _extents[i_dim];

      if (delta < x) {
        delta = x;
      }
    }

    if (delta > eps_extents) {
      delta *= tolerance;
    } else {
      delta = eps_extents;
    }

    for (int i_dim = 0; i_dim < 3; i_dim++) {
      _extents[i_dim]   -= delta;
      _extents[3+i_dim] += delta;
    }

    for (int k = 0; k < dim; k++) {
      global_extents[k]       = PDM_MIN(_extents[k    ], global_extents[k    ]);
      global_extents[dim + k] = PDM_MAX(_extents[dim+k], global_extents[dim+k]);
    }

  } /* End loop cell */
}


static
void
_compute_extents_2d_from_face_vtx
(
 int      n_face,
 int     *face_vtx_idx,
 int     *face_vtx,
 double  *vtx_coord,
 double  *box_extents,
 double  *global_extents
)
{
  const double tolerance   = 1.e-12;
  const double eps_extents = 1.e-7;
  const int dim = 3;

  /* Loop over face */
  for(int i_face = 0; i_face < n_face; ++i_face ) {

    double *_extents = box_extents + 6 * i_face;

    /* Init */
    for (int i_dim = 0; i_dim < 3; i_dim++) {
      _extents[  i_dim] =  HUGE_VAL;
      _extents[3+i_dim] = -HUGE_VAL;
    }

    /* Loop over face and vtx to compute bbox */
    for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
      int i_vtx = face_vtx[idx_vtx]-1;

      for (int i_dim = 0; i_dim < 3; i_dim++) {
        double x = vtx_coord[3*i_vtx + i_dim];

        if (x < _extents[i_dim]) {
          _extents[i_dim] = x;
        }
        if (x > _extents[3+i_dim]) {
          _extents[3+i_dim] = x;
        }
      }
    }

    double delta = 0.;
    for (int i_dim = 0; i_dim < 3; i_dim++) {
      double x = _extents[3+i_dim] - _extents[i_dim];

      if (delta < x) {
        delta = x;
      }
    }

    if (delta > eps_extents) {
      delta *= tolerance;
    } else {
      delta = eps_extents;
    }

    for (int i_dim = 0; i_dim < 3; i_dim++) {
      _extents[i_dim]   -= delta;
      _extents[3+i_dim] += delta;
    }

    for (int k = 0; k < dim; k++) {
      global_extents[k]       = PDM_MIN(_extents[k    ], global_extents[k    ]);
      global_extents[dim + k] = PDM_MAX(_extents[dim+k], global_extents[dim+k]);
    }

  } /* End loop cell */
}

static
void
_compute_part_mesh_extents
(
  PDM_part_mesh_t   *mesh,
  int                dim_mesh,
  double            *global_extents,
  double          ***extents_out
)
{
  int n_part = mesh->n_part;
  double **extents = malloc(n_part * sizeof(double *));
  if(dim_mesh == 3) {
    for(int i_part = 0; i_part < n_part; ++i_part) {

      int    *cell_face     = NULL;
      int    *cell_face_idx = NULL;
      int    *face_vtx_idx  = NULL;
      int    *face_vtx      = NULL;
      double *vtx_coord     = NULL;

      int n_cell = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_CELL);

      PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_CELL_FACE, &cell_face, &cell_face_idx, PDM_OWNERSHIP_USER);
      PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_VTX , &face_vtx , &face_vtx_idx, PDM_OWNERSHIP_USER);

      PDM_part_mesh_vtx_coord_get(mesh, i_part, &vtx_coord, PDM_OWNERSHIP_USER); // Il faudrait un unchanged

      extents[i_part] = malloc(6 * n_cell * sizeof(double));

      _compute_extents_3d(n_cell, cell_face_idx, cell_face, face_vtx_idx, face_vtx, vtx_coord, extents[i_part], global_extents);
    }
  } else if(dim_mesh == 2) {
    for(int i_part = 0; i_part < n_part; ++i_part) {

      int    *face_vtx      = NULL;
      int    *face_vtx_idx  = NULL;
      int    *face_edge_idx = NULL;
      int    *face_edge     = NULL;
      int    *edge_vtx_idx  = NULL;
      int    *edge_vtx      = NULL;
      double *vtx_coord     = NULL;

      // A gerer le cas mixte face_vtx ou face_edge + edge_vtx

      int n_face = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_FACE);
      PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_VTX , &face_vtx , &face_vtx_idx, PDM_OWNERSHIP_USER);

      PDM_part_mesh_vtx_coord_get(mesh, i_part, &vtx_coord, PDM_OWNERSHIP_USER); // Il faudrait un unchanged

      if(face_vtx != NULL) {
        _compute_extents_2d_from_face_vtx(n_face, face_vtx_idx, face_vtx, vtx_coord, extents[i_part], global_extents);
      } else {
        PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_EDGE , &face_edge, &face_edge_idx, PDM_OWNERSHIP_USER);
        assert(face_edge != NULL);
        PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX  , &edge_vtx, &edge_vtx_idx, PDM_OWNERSHIP_USER);
        assert(edge_vtx_idx == NULL);
        int n_edge = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_EDGE);

        edge_vtx_idx = malloc((n_edge+1) * sizeof(int));
        for(int i_edge = 0; i_edge < n_edge+1; ++i_edge){
          edge_vtx_idx[i_edge] = 2 * i_edge;
        }
        _compute_extents_3d(n_face, face_edge_idx, face_edge, edge_vtx_idx, edge_vtx, vtx_coord, extents[i_part], global_extents);
        // _compute_extents_2d_from_face_edge(n_face, face_edge_idx, face_edge, edge_vtx, vtx_coord, extents[i_part], global_extents);
        free(edge_vtx_idx);
      }
    }

  } else {
    int    *edge_vtx_idx  = NULL;
    int    *edge_vtx      = NULL;
    double *vtx_coord     = NULL;
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_edge = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_EDGE);

      PDM_part_mesh_vtx_coord_get(mesh, i_part, &vtx_coord, PDM_OWNERSHIP_USER); // Il faudrait un unchanged

      PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX  , &edge_vtx, &edge_vtx_idx, PDM_OWNERSHIP_USER);
      assert(edge_vtx_idx == NULL);
      edge_vtx_idx = malloc((n_edge+1) * sizeof(int));
      for(int i_edge = 0; i_edge < n_edge+1; ++i_edge){
        edge_vtx_idx[i_edge] = 2 * i_edge;
      }

      _compute_extents_2d_from_face_vtx(n_edge, edge_vtx_idx, edge_vtx, vtx_coord, extents[i_part], global_extents);

      free(edge_vtx_idx);
    }
  }
  *extents_out = extents;
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_mesh_intersection_t*
PDM_mesh_intersection_create
(
 const PDM_mesh_intersection_kind_t intersection_kind,
 const int                          dim_mesh_a,
 const int                          dim_mesh_b,
 const int                          n_part_mesh_a,
 const int                          n_part_mesh_b,
 const double                       project_coeff,
       PDM_MPI_Comm                 comm
)
{
  PDM_mesh_intersection_t *mi = (PDM_mesh_intersection_t *) malloc(sizeof(PDM_mesh_intersection_t));

  mi->comm = comm;
  mi->intersect_kind = intersection_kind;
  mi->n_part_mesh_a  = n_part_mesh_a;
  mi->n_part_mesh_b  = n_part_mesh_b;
  mi->dim_mesh_a     = dim_mesh_a;
  mi->dim_mesh_b     = dim_mesh_b;
  mi->project_coef   = project_coeff;

  mi->mesh_a = PDM_part_mesh_create(n_part_mesh_a, comm);
  mi->mesh_b = PDM_part_mesh_create(n_part_mesh_b, comm);

  return mi;
}

void
PDM_mesh_intersection_compute
(
  PDM_mesh_intersection_t  *mi
)
{
  /*
   * Compute extents of mesh_a and mesh_b
   */
  double **extents_mesh_a = NULL;
  double **extents_mesh_b = NULL;
  double global_extents[6] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL, HUGE_VAL,  HUGE_VAL,  HUGE_VAL};
  double mesh_global_extents[2][6] = {{ HUGE_VAL,  HUGE_VAL,  HUGE_VAL,
                                       -HUGE_VAL, -HUGE_VAL, -HUGE_VAL},
                                      {HUGE_VAL,  HUGE_VAL,  HUGE_VAL,
                                       -HUGE_VAL, -HUGE_VAL, -HUGE_VAL}};
  double g_mesh_global_extents[2][6];
  double g_global_extents        [6];
  _compute_part_mesh_extents(mi->mesh_a, mi->dim_mesh_a, mesh_global_extents[0], &extents_mesh_a);
  _compute_part_mesh_extents(mi->mesh_b, mi->dim_mesh_b, mesh_global_extents[1], &extents_mesh_b);

  /*
   * Global extents exchange
   */
  const intdim = 3;
  for(int i_mesh = 0; i_mesh < 2; ++i_mesh) {
    PDM_MPI_Allreduce(mesh_global_extents[i_mesh], g_mesh_global_extents[i_mesh], dim,
                      PDM_MPI_DOUBLE, PDM_MPI_MIN, mi->comm);
    PDM_MPI_Allreduce(mesh_global_extents[i_mesh]+dim, g_mesh_global_extents[i_mesh]+dim, dim,
                      PDM_MPI_DOUBLE, PDM_MPI_MAX, mi->comm);
  }

  /* Union or intersection of global extents */
  for(int i_mesh = 0; i_mesh < 2; ++i_mesh) {
    for (int k = 0; k < 3; k++) {
      // Union
      // global_extents[k]     = PDM_MIN(mesh_global_extents[i_mesh][k  ], global_extents[k  ]);
      // global_extents[3 + k] = PDM_MAX(mesh_global_extents[i_mesh][3+k], global_extents[3+k]);
      // Intersection
      global_extents[k]     = PDM_MAX(g_mesh_global_extents[i_mesh][k  ], global_extents[k  ]);
      global_extents[3 + k] = PDM_MIN(g_mesh_global_extents[i_mesh][3+k], global_extents[3+k]);
    }
  }
  for(int i = 0; i < 6; ++i) {
    g_global_extents[i] = global_extents[i];
  }
  double max_range = -HUGE_VAL;
  double min_range =  HUGE_VAL;

  for (int k = 0; k < dim; k++) {
    max_range = PDM_MAX(max_range, (g_global_extents[dim+k] - g_global_extents[k]));
    min_range = PDM_MIN(min_range, (g_global_extents[dim+k] - g_global_extents[k]));
  }

  for (int k = 0; k < dim; k++) {
    g_global_extents[k]     += -max_range * 1.1e-3; // On casse la symetrie !
    g_global_extents[dim+k] +=  max_range * 1e-3;
  }


  for(int i_part = 0; i_part < mi->n_part_mesh_a; ++i_part) {
    free(extents_mesh_a[i_part]);
  }
  for(int i_part = 0; i_part < mi->n_part_mesh_b; ++i_part) {
    free(extents_mesh_b[i_part]);
  }
  free(extents_mesh_a);
  free(extents_mesh_b);

}

void
PDM_mesh_intersection_part_set
(
  PDM_mesh_intersection_t  *mi,
  PDM_ol_mesh_t             i_mesh,
  int                       i_part,
  int                       n_cell,
  int                       n_face,
  int                       n_edge,
  int                       n_vtx,
  int                      *cell_face_idx,
  int                      *cell_face,
  int                      *face_edge_idx,
  int                      *face_edge,
  int                      *edge_vtx,
  int                      *face_vtx_idx,
  int                      *face_vtx,
  PDM_g_num_t              *cell_ln_to_gn,
  PDM_g_num_t              *face_ln_to_gn,
  PDM_g_num_t              *edge_ln_to_gn,
  PDM_g_num_t              *vtx_ln_to_gn,
  double                   *vtx_coord
)
{
  PDM_part_mesh_t* mesh = mi->mesh_a;
  if(i_mesh == PDM_OL_MESH_B) {
    mesh = mi->mesh_b;
  }

  PDM_part_mesh_n_entity_set(mesh, PDM_MESH_ENTITY_CELL  , i_part, n_cell);
  PDM_part_mesh_n_entity_set(mesh, PDM_MESH_ENTITY_FACE  , i_part, n_face);
  PDM_part_mesh_n_entity_set(mesh, PDM_MESH_ENTITY_EDGE  , i_part, n_edge);
  PDM_part_mesh_n_entity_set(mesh, PDM_MESH_ENTITY_VERTEX, i_part, n_vtx );

  PDM_part_mesh_connectivity_set(mesh, PDM_CONNECTIVITY_TYPE_CELL_FACE, i_part, cell_face, cell_face_idx, PDM_OWNERSHIP_USER);
  PDM_part_mesh_connectivity_set(mesh, PDM_CONNECTIVITY_TYPE_FACE_EDGE, i_part, face_edge, face_edge_idx, PDM_OWNERSHIP_USER);
  PDM_part_mesh_connectivity_set(mesh, PDM_CONNECTIVITY_TYPE_FACE_VTX , i_part, face_vtx , face_vtx_idx , PDM_OWNERSHIP_USER);
  PDM_part_mesh_connectivity_set(mesh, PDM_CONNECTIVITY_TYPE_EDGE_VTX , i_part, edge_vtx , NULL         , PDM_OWNERSHIP_USER);

  PDM_part_mesh_entity_ln_to_gn_set(mesh, i_part, PDM_MESH_ENTITY_CELL  , cell_ln_to_gn, PDM_OWNERSHIP_USER);
  PDM_part_mesh_entity_ln_to_gn_set(mesh, i_part, PDM_MESH_ENTITY_FACE  , face_ln_to_gn, PDM_OWNERSHIP_USER);
  PDM_part_mesh_entity_ln_to_gn_set(mesh, i_part, PDM_MESH_ENTITY_EDGE  , edge_ln_to_gn, PDM_OWNERSHIP_USER);
  PDM_part_mesh_entity_ln_to_gn_set(mesh, i_part, PDM_MESH_ENTITY_VERTEX, vtx_ln_to_gn , PDM_OWNERSHIP_USER);

  PDM_part_mesh_vtx_coord_set(mesh, i_part, vtx_coord, PDM_OWNERSHIP_USER);
}



void
PDM_mesh_intersection_free
(
 PDM_mesh_intersection_t* mi
)
{

  PDM_part_mesh_free(mi->mesh_a);
  PDM_part_mesh_free(mi->mesh_b);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
