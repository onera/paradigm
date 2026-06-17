/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_part_mesh_geom.h"
#include "pdm_part_mesh_priv.h"
#include "pdm_error.h"
#include "pdm_predicate.h"
#include "pdm_mem_tool.h"
#include "pdm_logging.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


/*=============================================================================
 * Public function definitions
 *============================================================================*/


void
PDM_compute_entity_center
(
  int     n_entity,
  int    *connect_idx,
  int    *connect,
  double *coord,
  double *entity_center
)
{
  for (int i_entity = 0; i_entity < n_entity; i_entity++) {
    double *c = &entity_center[3*i_entity];
    c[0] = c[1] = c[2] = 0.;

    double normalization = 1./(connect_idx[i_entity+1] - connect_idx[i_entity]);

    for (int i = connect_idx[i_entity]; i < connect_idx[i_entity+1]; i++) {
      int idx = PDM_ABS(connect[i]) - 1;
      for (int j = 0; j < 3; j++) {
        c[j] += coord[3*idx+j];
      }
    }

    for (int j = 0; j < 3; j++) {
      c[j] *= normalization;
    }
  }
}

void
PDM_compute_dual_volume_ngon_2d
(
  int       n_part,
  int      *n_face,
  int      *n_edge,
  int      *n_vtx,
  int     **face_edge_idx,
  int     **face_edge,
  int     **edge_vtx,
  double  **vtx_coord,
  double ***out_vtx_volume
)
{
  double **vtx_volume = NULL;
  PDM_malloc(vtx_volume, n_part, double *);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_calloc(vtx_volume[i_part], n_vtx[i_part], double);

    // Compute edge centers
    double *edge_center = NULL;
    PDM_malloc(edge_center, 3*n_edge[i_part], double);
    for (int i_edge = 0; i_edge < n_edge[i_part]; i_edge++) {
      int i_vtx0 = edge_vtx[i_part][2*i_edge  ] - 1;
      int i_vtx1 = edge_vtx[i_part][2*i_edge+1] - 1;

      for (int i = 0; i < 3; i++) {
        edge_center[3*i_edge+i] = 0.5*(vtx_coord[i_part][3*i_vtx0+i] + vtx_coord[i_part][3*i_vtx1+i]);
      }
    }

    // Sum up contributions from each triplet (face, edge, vtx)
    for (int i_face = 0; i_face < n_face[i_part]; i_face++) {

      double face_center[3];
      PDM_compute_entity_center(1,
                                &face_edge_idx[i_part][i_face],
                                face_edge     [i_part],
                                edge_center,
                                face_center);

      for (int idx_edge = face_edge_idx[i_part][i_face]; idx_edge < face_edge_idx[i_part][i_face+1]; idx_edge++) {

        int    i_edge = PDM_ABS (face_edge[i_part][idx_edge]) - 1;
        double sign   = PDM_SIGN(face_edge[i_part][idx_edge]);

        for (int idx_vtx = 2*i_edge; idx_vtx < 2*(i_edge+1); idx_vtx++) {
          int i_vtx = edge_vtx[i_part][idx_vtx] - 1;

          // Assume 2D in xy-plane
          double area = PDM_predicate_orient2d(&vtx_coord[i_part][3*i_vtx],
                                               &edge_center[3*i_edge],
                                                face_center);

          vtx_volume[i_part][i_vtx] += 0.5 * sign * area;
          sign = -sign;

        } // End loop on vertices
      } // End loop on edges
    } // End loop on faces

    PDM_free(edge_center);

  } // End loop on parts

  *out_vtx_volume = vtx_volume;
}

void
PDM_compute_dual_volume_ngon_3d
(
  int       n_part,
  int      *n_cell,
  int      *n_face,
  int      *n_edge,
  int      *n_vtx,
  int     **cell_face_idx,
  int     **cell_face,
  int     **face_edge_idx,
  int     **face_edge,
  int     **edge_vtx,
  double  **vtx_coord,
  double ***out_vtx_volume
)
{
  double **vtx_volume = NULL;
  PDM_malloc(vtx_volume, n_part, double *);

  const double one_sixth = 1./6.;

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_calloc(vtx_volume[i_part], n_vtx[i_part], double);

    double *edge_center = NULL;
    double *face_center = NULL;
    PDM_malloc(edge_center, 3*n_edge[i_part], double);
    PDM_malloc(face_center, 3*n_face[i_part], double);

    // Compute edge centers
    for (int i_edge = 0; i_edge < n_edge[i_part]; i_edge++) {
      int i_vtx0 = edge_vtx[i_part][2*i_edge  ] - 1;
      int i_vtx1 = edge_vtx[i_part][2*i_edge+1] - 1;

      for (int i = 0; i < 3; i++) {
        edge_center[3*i_edge+i] = 0.5*(vtx_coord[i_part][3*i_vtx0+i] + vtx_coord[i_part][3*i_vtx1+i]);
      }
    }

    // Compute face centers
    PDM_compute_entity_center(n_face       [i_part],
                              face_edge_idx[i_part],
                              face_edge    [i_part],
                              edge_center,
                              face_center);

    // Sum up contributions from each quadruplet (cell, face, edge, vtx)
    for (int i_cell = 0; i_cell < n_cell[i_part]; i_cell++) {

      double cell_center[3];
      PDM_compute_entity_center(1,
                                &cell_face_idx[i_part][i_cell],
                                cell_face    [i_part],
                                face_center,
                                cell_center);

      for (int idx_face = cell_face_idx[i_part][i_cell]; idx_face < cell_face_idx[i_part][i_cell+1]; idx_face++) {

        int    i_face    = PDM_ABS (cell_face[i_part][idx_face]) - 1;
        double sign_face = PDM_SIGN(cell_face[i_part][idx_face]);

        for (int idx_edge = face_edge_idx[i_part][i_face]; idx_edge < face_edge_idx[i_part][i_face+1]; idx_edge++) {

          int    i_edge    = PDM_ABS (face_edge[i_part][idx_edge]) - 1;
          double sign_edge = PDM_SIGN(face_edge[i_part][idx_edge]);

          for (int idx_vtx = 2*i_edge; idx_vtx < 2*(i_edge+1); idx_vtx++) {

            int i_vtx = edge_vtx[i_part][idx_vtx] - 1;

            double vol = PDM_predicate_orient3d(&vtx_coord[i_part][3*i_vtx],
                                                &face_center[3*i_face],
                                                &edge_center[3*i_edge],
                                                 cell_center);

            vtx_volume[i_part][i_vtx] += one_sixth * sign_face * sign_edge * vol;
            sign_edge = -sign_edge;

          } // End loop on vertices
        } // End loop on edges
      } // End loop on faces
    } // End loop on cells

    PDM_free(edge_center);
    PDM_free(face_center);

  } // End loop on parts

  *out_vtx_volume = vtx_volume;
}



void
PDM_part_mesh_dual_volume_compute
(
  PDM_part_mesh_t   *pm,
  PDM_bool_t         synchronize,
  double          ***out_dual_vol
)
{

  /* Compute mesh highest dimension */
  int mesh_dimension = 2;
  int tn_cell = 0;
  if (pm->pn_entity[PDM_MESH_ENTITY_CELL] != NULL) {
    for (int i_part = 0; i_part < pm->n_part; i_part++) {
      tn_cell += pm->pn_entity[PDM_MESH_ENTITY_CELL][i_part];
    }
  }

  int max_tn_cell = 0;
  PDM_MPI_Allreduce(&tn_cell, &max_tn_cell, 1,
                    PDM_MPI_INT, PDM_MPI_MAX, pm->comm);

  if (max_tn_cell > 0) {
    mesh_dimension = 3;
  }

  int     *n_cell        = NULL;
  int     *n_face        = NULL;
  int     *n_edge        = NULL;
  int     *n_vtx         = NULL;
  int    **cell_face_idx = NULL;
  int    **cell_face     = NULL;
  int    **face_edge_idx = NULL;
  int    **face_edge     = NULL;
  int    **edge_vtx      = NULL;
  double **vtx_coord     = NULL;
  if (mesh_dimension == 3) {
    PDM_malloc(n_cell,        pm->n_part, int  );
    PDM_malloc(cell_face_idx, pm->n_part, int *);
    PDM_malloc(cell_face,     pm->n_part, int *);
  }
  PDM_malloc(n_face       , pm->n_part, int     );
  PDM_malloc(n_edge       , pm->n_part, int     );
  PDM_malloc(n_vtx        , pm->n_part, int     );
  PDM_malloc(face_edge_idx, pm->n_part, int    *);
  PDM_malloc(face_edge    , pm->n_part, int    *);
  PDM_malloc(edge_vtx     , pm->n_part, int    *);
  PDM_malloc(vtx_coord    , pm->n_part, double *);
  for (int i_part = 0; i_part < pm->n_part; i_part++) {
    if (mesh_dimension == 3) {
      n_cell[i_part] = PDM_part_mesh_n_entity_get(pm,
                                                  i_part,
                                                  PDM_MESH_ENTITY_CELL);
      PDM_part_mesh_connectivity_get(pm,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                     &cell_face    [i_part],
                                     &cell_face_idx[i_part],
                                     PDM_OWNERSHIP_BAD_VALUE);
    }

    n_vtx[i_part] = PDM_part_mesh_n_entity_get(pm,
                                               i_part,
                                               PDM_MESH_ENTITY_VTX);
    PDM_part_mesh_vtx_coord_get(pm,
                                i_part,
                                &vtx_coord[i_part],
                                PDM_OWNERSHIP_BAD_VALUE);

    n_face[i_part] = PDM_part_mesh_n_entity_get(pm,
                                                i_part,
                                                PDM_MESH_ENTITY_FACE);
    PDM_part_mesh_connectivity_get(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                   &face_edge    [i_part],
                                   &face_edge_idx[i_part],
                                   PDM_OWNERSHIP_BAD_VALUE);

    n_edge[i_part] = PDM_part_mesh_n_entity_get(pm,
                                                i_part,
                                                PDM_MESH_ENTITY_EDGE);
    int *edge_vtx_idx = NULL;
    PDM_part_mesh_connectivity_get(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                   &edge_vtx    [i_part],
                                   &edge_vtx_idx,
                                   PDM_OWNERSHIP_BAD_VALUE);
    // PDM_free(edge_vtx_idx);
  }

  if (mesh_dimension == 2) {
    PDM_compute_dual_volume_ngon_2d(pm->n_part,
                                    n_face,
                                    n_edge,
                                    n_vtx,
                                    face_edge_idx,
                                    face_edge,
                                    edge_vtx,
                                    vtx_coord,
                                    out_dual_vol);
  } else {
    PDM_compute_dual_volume_ngon_3d(pm->n_part,
                                    n_cell,
                                    n_face,
                                    n_edge,
                                    n_vtx,
                                    cell_face_idx,
                                    cell_face,
                                    face_edge_idx,
                                    face_edge,
                                    edge_vtx,
                                    vtx_coord,
                                    out_dual_vol);
  }

  if (mesh_dimension == 3) {
    PDM_free(n_cell);
    PDM_free(cell_face_idx);
    PDM_free(cell_face);
  }
  PDM_free(n_face);
  PDM_free(n_edge);
  PDM_free(n_vtx);
  PDM_free(face_edge_idx);
  PDM_free(face_edge);
  PDM_free(edge_vtx);
  PDM_free(vtx_coord);

  if (synchronize) {
    if(pm->pcg[PDM_MESH_ENTITY_VTX] == NULL) {
      PDM_part_mesh_part_comm_graph_compute_from_gnum(pm, PDM_MESH_ENTITY_VTX);
    }

    // Synchro volume :
    PDM_part_comm_graph_all_reduce(pm->pcg[PDM_MESH_ENTITY_VTX],
                                   PDM_MPI_DOUBLE,
                                   1,
                                   PDM_MPI_SUM,
            (unsigned char **)    *out_dual_vol);
  }

}



#ifdef __cplusplus
}
#endif /* __cplusplus */
