
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_array.h"
#include "pdm_part_comm_graph.h"
#include "pdm_part_comm_graph_algorithm.h"
#include "pdm_mem_tool.h"
#include "pdm_error.h"
#include "pdm_laplacian_smoothing.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
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

void PDM_laplacian_smoothing_idw_weights_compute
(
  int       n_part,
  double  **p_vtx_coord,
  int      *p_n_edge,
  int     **p_edge_vtx,
  int       exponent,
  double ***out_p_edge_weight
)
{
  double **p_edge_weight = NULL;
  PDM_malloc(p_edge_weight, n_part, double *);
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_malloc(p_edge_weight[i_part], p_n_edge[i_part], double);
    for (int i_edge = 0; i_edge < p_n_edge[i_part]; i_edge++) {
      int i_vtx = p_edge_vtx[i_part][2*i_edge  ] - 1;
      int j_vtx = p_edge_vtx[i_part][2*i_edge+1] - 1;
      double v[3];
      for (int i = 0; i < 3; i++) {
        v[i] = p_vtx_coord[i_part][3*j_vtx+i] - p_vtx_coord[i_part][3*i_vtx+i];
      }
      p_edge_weight[i_part][i_edge] = PDM_MODULE(v);
      p_edge_weight[i_part][i_edge] = 1./pow(p_edge_weight[i_part][i_edge], exponent);
    }
  }
  *out_p_edge_weight = p_edge_weight;
}

//void PDM_laplacian_smoothing_beltrami_weights_compute
//(
//  int       n_part,
//  int      *p_n_vtx,
//  double  **p_vtx_coord,
//  int      *p_n_elt,
//  int     **p_elt_vtx_idx,
//  int     **p_elt_vtx,
//  int      *p_n_edge,
//  int     **p_edge_vtx,
//  double ***out_p_edge_weight
//)
//{
//  PDM_UNUSED(n_part);
//  PDM_UNUSED(p_n_vtx);
//  PDM_UNUSED(p_vtx_coord);
//  PDM_UNUSED(p_n_elt);
//  PDM_UNUSED(p_elt_vtx_idx);
//  PDM_UNUSED(p_elt_vtx);
//  PDM_UNUSED(p_n_edge);
//  PDM_UNUSED(p_edge_vtx);
//  PDM_UNUSED(out_p_edge_weight);
//  PDM_error(__FILE__, __LINE__, 0, "Beltrami weights not implemented yet\n");
//}

void
PDM_laplacian_smoothing_fields
(
  const PDM_MPI_Comm            comm,
        int                     n_part,
        int                    *p_n_vtx,
        PDM_part_comm_graph_t  *pcg_vtx,
        int                    *p_n_vtx_frozen,
        int                   **p_vtx_frozen,
        int                    *p_n_edge,
        int                   **p_edge_vtx,
        double                **p_edge_weight,
        PDM_part_comm_graph_t  *pcg_edge,
        double                  damping,
        int                     n_iter,
        double                  tol,
        int                     stride,
        double                **p_vtx_field
)
{
  // Generate edge part comm graph if not provided
  int has_pcg_edge = 1;
  if (pcg_edge == NULL) {
    has_pcg_edge = 0;
    int **p_edge_vtx_idx = NULL;
    PDM_malloc(p_edge_vtx_idx, n_part, int *);
    for (int i_part = 0; i_part < n_part; i_part++) {
      p_edge_vtx_idx[i_part] = PDM_array_new_idx_from_const_stride_int(2, p_n_edge[i_part]);
    }
    PDM_part_comm_graph_entity1_to_part_comm_graph_entity2(pcg_vtx,
                                                           p_n_vtx,
                                                           p_n_edge,
                                                           p_edge_vtx_idx,
                                                           p_edge_vtx,
                                                           &pcg_edge);
    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_free(p_edge_vtx_idx[i_part]);
    }
    PDM_free(p_edge_vtx_idx);
  }

  // Compute default edge weights if not provided
  int has_edge_weight = 1;
  if (p_edge_weight == NULL) {
    has_edge_weight = 0;
    PDM_malloc(p_edge_weight, n_part, double *);
    for (int i_part = 0; i_part < n_part; i_part++) {
      p_edge_weight[i_part] = PDM_array_const_double(p_n_edge[i_part], 1.);
    }
  }

  // Prepare to disable ghost edges without in-place modification of user edge weights
  // Could be a function of part_comm_graph ? with a simple getter
  int  *p_n_edge_lap = NULL;
  int **p_edge_lap   = NULL;
  PDM_malloc(p_n_edge_lap, n_part, int  );
  PDM_malloc(p_edge_lap,   n_part, int *);
  for (int i_part = 0; i_part < n_part; i_part++) {
    const int *edge_bound_owner = PDM_part_comm_graph_owner_get(pcg_edge, i_part);

    int *edge_bound = NULL;
    int n_edge_bound = PDM_part_comm_graph_entity_graph_get(pcg_edge,
                                                            i_part,
                                                            &edge_bound,
                                                            PDM_OWNERSHIP_BAD_VALUE);

    PDM_malloc(p_edge_lap[i_part], p_n_edge[i_part], int);
    p_n_edge_lap[i_part] = 0;
    for (int i_bnd = 0; i_bnd < n_edge_bound; i_bnd++) {
      if (edge_bound_owner[i_bnd] == 0) {
        int i_edge = edge_bound[4*i_bnd] - 1;
        p_edge_lap[i_part][p_n_edge_lap[i_part]++] = i_edge+1;
      }
    }
    PDM_realloc(p_edge_lap[i_part], p_edge_lap[i_part], p_n_edge_lap[i_part], int);
  }

  // Pre-compute vertex weights
  double **p_vtx_weight = NULL;
  PDM_malloc(p_vtx_weight, n_part, double *);
  for (int i_part = 0; i_part < n_part; i_part++) {
    p_vtx_weight[i_part] = PDM_array_const_double(p_n_vtx[i_part], 0.);
    for (int i_edge_lap = 0; i_edge_lap < p_n_edge_lap[i_part]; i_edge_lap++) {
      int i_edge = p_edge_lap[i_part][i_edge_lap] - 1;
      for (int i = 0; i < 2; i++) {
        int i_vtx = p_edge_vtx[i_part][2*i_edge+i] - 1;
        p_vtx_weight[i_part][i_vtx] += p_edge_weight[i_part][i_edge];
      }
    }
  }

  // Synchronize vertex weights
  PDM_part_comm_graph_all_reduce(pcg_vtx,
                                 PDM_MPI_DOUBLE,
                                 PDM_MPI_SUM,
              (unsigned char **) p_vtx_weight);

  // Internal field arrays
  double **field_tmp = NULL;
  PDM_malloc(field_tmp, n_part, double *);
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_malloc(field_tmp[i_part], stride * p_n_vtx[i_part], double);
    memcpy(field_tmp[i_part], p_vtx_field[i_part], sizeof(double) * stride * p_n_vtx[i_part]);
  }

  // Iterations
  int    iter = 0;
  double eps  = HUGE_VAL;

  double **field_current = p_vtx_field;
  double **field_prev    = field_tmp;

  while (iter < n_iter && eps > tol) {

    for (int i_part = 0; i_part < n_part; i_part++) {

      // Init fields
      for (int i_val = 0; i_val < stride * p_n_vtx[i_part]; i_val++) {
        field_current[i_part][i_val] = 0.;
      }

      // Laplacian with damping
      for (int i_edge_lap = 0; i_edge_lap < p_n_edge_lap[i_part]; i_edge_lap++) {
        int i_edge = p_edge_lap[i_part][i_edge_lap] - 1;
        for (int i = 0; i < 2; i++) {
          int i_vtx = p_edge_vtx[i_part][2*i_edge+ i     ] - 1;
          int j_vtx = p_edge_vtx[i_part][2*i_edge+(i+1)%2] - 1;
          for (int i_stride = 0; i_stride < stride; i_stride++) {
            field_current[i_part][stride*i_vtx+i_stride] += damping * p_edge_weight[i_part][i_edge] * field_prev[i_part][stride*j_vtx+i_stride];
          }
        }
      }

    }

    // Synchronize fields
    PDM_part_comm_graph_all_reduce_strided(pcg_vtx,
                                           PDM_MPI_DOUBLE,
                                           PDM_MPI_SUM,
                                           stride,
                        (unsigned char **) field_current);

    // Normalize fields
    for (int i_part = 0; i_part < n_part; i_part++) {
      for (int i_vtx = 0; i_vtx < p_n_vtx[i_part]; i_vtx++) {
        for (int i_stride = 0; i_stride < stride; i_stride++) {
          field_current[i_part][stride*i_vtx+i_stride] /= p_vtx_weight[i_part][i_vtx];
        }
      }
    }

    // Reset previous values on vtx group
    if (p_n_vtx_frozen  != NULL) {
      assert(p_vtx_frozen != NULL);
      for (int i_part = 0; i_part < n_part; i_part++) {
        for (int i_vtx_frozen = 0; i_vtx_frozen < p_n_vtx_frozen[i_part]; i_vtx_frozen++) {
          int i_vtx = p_vtx_frozen[i_part][i_vtx_frozen]-1;
          for (int i_stride = 0; i_stride < stride; i_stride++) {
            field_current[i_part][stride*i_vtx+i_stride] = field_prev[i_part][stride*i_vtx+i_stride];
          }
        }
      }
    }

    // Swap current and previous fields
    double **tmp_swap = field_prev;
    field_prev        = field_current;
    field_current     = tmp_swap;

    // Check for convergence
    if (tol > 0.) {
      double _eps = 0.;
      for (int i_part = 0; i_part < n_part; i_part++) {
        for (int i_vtx = 0; i_vtx < p_n_vtx[i_part]; i_vtx++) {
          for (int i_stride = 0; i_stride < stride; i_stride++) {
            _eps = PDM_MAX(_eps, (field_current[i_part][stride*i_vtx+i_stride]-field_prev[i_part][stride*i_vtx+i_stride])/PDM_MAX(field_prev[i_part][stride*i_vtx+i_stride], 1e-16));
          }
        }
      }
      PDM_MPI_Allreduce(&_eps, &eps, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);
    }

    iter++;

  } // end while iterations

  // Update field
  if (field_current != p_vtx_field) {
    for (int i_part = 0; i_part < n_part; i_part++) {
    memcpy(p_vtx_field[i_part], field_current[i_part], sizeof(double) * stride * p_n_vtx[i_part]);
    }
  }

  // Free
  if (has_pcg_edge == 0) {
    PDM_part_comm_graph_free(pcg_edge);
  }
  if (has_edge_weight == 0) {
    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_free(p_edge_weight[i_part]);
    }
    PDM_free(p_edge_weight);
  }
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_free(field_tmp   [i_part]);
    PDM_free(p_vtx_weight[i_part]);
    PDM_free(p_edge_lap  [i_part]);
  }
  PDM_free(field_tmp   );
  PDM_free(p_vtx_weight);
  PDM_free(p_edge_lap  );
  PDM_free(p_n_edge_lap);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
