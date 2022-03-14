/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_part.h"
#include "pdm_timer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_array.h"
#include "pdm_gnum.h"
#include "pdm_extract_part.h"
#include "pdm_extract_part_priv.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_plane.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

static
void
_face_center_2d
(
  int       n_part_in,
  int      *n_extract,
  int     **extract_lnum,
  int     **pface_edge_idx,
  int     **pface_edge,
  int     **pedge_vtx,
  double  **pvtx_coord,
  double ***face_center
)
{
  for(int i_part = 0; i_part < n_part_in; ++i_part) {
    assert(pface_edge    [i_part] != NULL);
    assert(pface_edge_idx[i_part] != NULL);
    assert(pedge_vtx     [i_part] != NULL);
    assert(pvtx_coord    [i_part] != NULL);
  }

  double** entity_center = malloc(n_part_in * sizeof(double * ));
  for(int i_part = 0; i_part < n_part_in; ++i_part) {
    entity_center[i_part] = (double *) malloc(3 * n_extract[i_part] * sizeof(double));

    int    *_pface_edge     = pface_edge    [i_part];
    int    *_pface_edge_idx = pface_edge_idx[i_part];
    int    *_pedge_vtx      = pedge_vtx     [i_part];
    double *_pvtx_coord     = pvtx_coord    [i_part];

    for(int idx_face = 0; idx_face < n_extract[i_part]; ++idx_face) {

      int i_face = extract_lnum[i_part][idx_face];
      entity_center[i_part][3*i_face  ] = 0.;
      entity_center[i_part][3*i_face+1] = 0.;
      entity_center[i_part][3*i_face+2] = 0.;

      double inv = 1./((double)  _pface_edge_idx[i_face+1] - _pface_edge_idx[i_face]);

      for(int idx_edge = _pface_edge_idx[i_face]; idx_edge < _pface_edge_idx[i_face+1]; ++idx_edge) {
        int i_edge = _pface_edge[idx_edge];
        int i_vtx1 = _pedge_vtx[2*i_edge  ] - 1;
        int i_vtx2 = _pedge_vtx[2*i_edge+1] - 1;

        entity_center[i_part][3*i_face  ] += 0.5 * (_pvtx_coord[i_vtx1] + _pvtx_coord[i_vtx2]);
        entity_center[i_part][3*i_face+1] += 0.5 * (_pvtx_coord[i_vtx1] + _pvtx_coord[i_vtx2]);
        entity_center[i_part][3*i_face+2] += 0.5 * (_pvtx_coord[i_vtx1] + _pvtx_coord[i_vtx2]);

      }
      entity_center[i_part][3*i_face  ] = entity_center[i_part][3*i_face  ] * inv;
      entity_center[i_part][3*i_face+1] = entity_center[i_part][3*i_face+1] * inv;
      entity_center[i_part][3*i_face+2] = entity_center[i_part][3*i_face+2] * inv;
    }
  }

  *face_center = entity_center;
}


static
void
_cell_center_3d
(
  int       n_part_in,
  int      *n_extract,
  int     **extract_lnum,
  int     **pcell_face_idx,
  int     **pcell_face,
  int     **pface_edge_idx,
  int     **pface_edge,
  int     **pface_vtx_idx,
  int     **pface_vtx,
  int     **pedge_vtx,
  double  **pvtx_coord,
  double ***cell_center
)
{
  int from_edge = 0;
  int from_face = 0;
  for(int i_part = 0; i_part < n_part_in; ++i_part) {
    if(pface_edge    [i_part] != NULL) {
      from_edge = 1;
    }
    if(pface_vtx    [i_part] != NULL) {
      from_face = 1;
    }
    assert(pvtx_coord    [i_part] != NULL);
  }

  double** entity_center = malloc(n_part_in * sizeof(double * ));

  if(from_face == 1) {
    for(int i_part = 0; i_part < n_part_in; ++i_part) {
      entity_center[i_part] = (double *) malloc(3 * n_extract[i_part] * sizeof(double));

      int    *_pface_edge     = pface_edge    [i_part];
      int    *_pface_edge_idx = pface_edge_idx[i_part];
      int    *_pedge_vtx      = pedge_vtx     [i_part];
      double *_pvtx_coord     = pvtx_coord    [i_part];

      for(int idx_cell = 0; idx_cell < n_extract[i_part]; ++idx_cell) {
        int i_cell = extract_lnum[i_part][idx_cell];
      }
    }
  } else if( from_edge == 1) {

  }


  *cell_center = entity_center;
}



/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_extract_part_t*
PDM_extract_part_create
(
 const int                    dim,
 const int                    n_part_in,
 const int                    n_part_out,
       PDM_split_dual_t       split_dual_method,
       PDM_ownership_t        ownership,
       PDM_MPI_Comm           comm
)
{
  PDM_extract_part_t *extrp = (PDM_extract_part_t *) malloc(sizeof(PDM_extract_part_t));

  extrp->dim               = dim;
  extrp->n_part_in         = n_part_in;
  extrp->n_part_out        = n_part_out;
  extrp->split_dual_method = split_dual_method;
  extrp->ownership         = ownership;
  extrp->comm              = comm;

  extrp->n_cell         = (int          *) malloc(n_part_in * sizeof(int          ));
  extrp->n_face         = (int          *) malloc(n_part_in * sizeof(int          ));
  extrp->n_edge         = (int          *) malloc(n_part_in * sizeof(int          ));
  extrp->n_vtx          = (int          *) malloc(n_part_in * sizeof(int          ));
  extrp->n_extract      = (int          *) malloc(n_part_in * sizeof(int          ));

  extrp->pcell_face     = (int         **) malloc(n_part_in * sizeof(int         *));
  extrp->pcell_face_idx = (int         **) malloc(n_part_in * sizeof(int         *));
  extrp->pface_edge     = (int         **) malloc(n_part_in * sizeof(int         *));
  extrp->pface_edge_idx = (int         **) malloc(n_part_in * sizeof(int         *));
  extrp->pedge_vtx      = (int         **) malloc(n_part_in * sizeof(int         *));

  extrp->extract_lnum   = (int         **) malloc(n_part_in * sizeof(int         *));

  extrp->cell_ln_to_gn  = (PDM_g_num_t **) malloc(n_part_in * sizeof(PDM_g_num_t *));
  extrp->face_ln_to_gn  = (PDM_g_num_t **) malloc(n_part_in * sizeof(PDM_g_num_t *));
  extrp->edge_ln_to_gn  = (PDM_g_num_t **) malloc(n_part_in * sizeof(PDM_g_num_t *));
  extrp->vtx_ln_to_gn   = (PDM_g_num_t **) malloc(n_part_in * sizeof(PDM_g_num_t *));

  extrp->pface_vtx_idx  = (int         **) malloc(n_part_in * sizeof(int         *));
  extrp->pface_vtx      = (int         **) malloc(n_part_in * sizeof(int         *));

  extrp->pvtx_coord     = (double      **) malloc(n_part_in * sizeof(double      *));

  return extrp;
}



void
PDM_extract_part_compute
(
  PDM_extract_part_t        *extrp
)
{
  assert(extrp->dim >= 2);

  int          *pn_entity    = 0;
  PDM_g_num_t **entity_g_num = NULL;
  if(extrp->dim == 3) {
    pn_entity    = extrp->n_cell;
    entity_g_num = extrp->cell_ln_to_gn;
  } else {
    pn_entity    = extrp->n_face;
    entity_g_num = extrp->face_ln_to_gn;
  }

  PDM_UNUSED(pn_entity);

  /*
   *  Create array selected in gnum
   */
  PDM_g_num_t** entity_extract_g_num = (PDM_g_num_t **) malloc( extrp->n_part_in * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    entity_extract_g_num[i_part] = (PDM_g_num_t *) malloc( extrp->n_extract[i_part] * sizeof(PDM_g_num_t));
    for(int i_entity = 0; i_entity < extrp->n_extract[i_part]; ++i_entity) {
      entity_extract_g_num[i_part][i_entity] = entity_g_num[i_part][extrp->extract_lnum[i_part][i_entity]];
    }
    if(0 == 1) {
      PDM_log_trace_array_long(entity_extract_g_num[i_part], extrp->n_extract[i_part], "entity_extract_g_num ::" );
    }
  }


  /*
   * Calcul des coordonnées to setup hilbert ordering (independant of parallelism )
   */
  double **entity_center = NULL;
  if(extrp->split_dual_method == PDM_SPLIT_DUAL_WITH_HILBERT) {
    if(extrp->dim == 2) {
      // Compute entity_center with face_edge + edge_vtx
      _face_center_2d(extrp->n_part_in,
                      extrp->n_extract,
                      extrp->extract_lnum,
                      extrp->pface_edge_idx,
                      extrp->pface_edge,
                      extrp->pedge_vtx,
                      extrp->pvtx_coord,
                      &entity_center);
    } else {  // dim == 3

      _cell_center_3d(extrp->n_part_in,
                      extrp->n_extract,
                      extrp->extract_lnum,
                      extrp->pcell_face_idx,
                      extrp->pcell_face,
                      extrp->pface_edge_idx,
                      extrp->pface_edge,
                      extrp->pface_vtx_idx,
                      extrp->pface_vtx,
                      extrp->pedge_vtx,
                      extrp->pvtx_coord,
                      &entity_center);

    }
  }


  /*
   * Si scotch ou metis, on doit calculer le dcell_face + pdm_gnum (pour avoir les faces child)
   *   Puis dconnectiviy_transpose puis combine
   */



  /*
   * Global numering computation
   */
  // PDM_g_num_t **child_selected_g_num = (PDM_g_num_t **) malloc( n_part_zones * sizeof(PDM_g_num_t *));
  // for (int i_part = 0; i_part < n_part_zones; i_part++){
  //   PDM_gnum_set_from_coords(gnum_extract, i_part, pn_select_cell[i_part], tmp_extract_cell_center[i_part], NULL);
  // }

  // PDM_gnum_compute(gnum_extract);

  // for (int i_part = 0; i_part < n_part_zones; i_part++){
  //   child_selected_g_num[i_part] = PDM_gnum_get(gnum_extract, i_part);
  //   // PDM_log_trace_array_long(child_selected_g_num[i_part], pn_select_cell[i_part], "child_selected_g_num : ");
  // }
  // PDM_gnum_free(gnum_extract);




  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    free(entity_extract_g_num[i_part]);
  }
  free(entity_extract_g_num);

}


void
PDM_extract_part_part_set
(
  PDM_extract_part_t       *extrp,
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
  extrp->n_cell        [i_part] = n_cell;
  extrp->n_face        [i_part] = n_face;
  extrp->n_edge        [i_part] = n_edge;
  extrp->n_vtx         [i_part] = n_vtx;
  extrp->pcell_face    [i_part] = cell_face;
  extrp->pcell_face_idx[i_part] = cell_face_idx;
  extrp->pface_edge    [i_part] = face_edge;
  extrp->pface_edge_idx[i_part] = face_edge_idx;
  extrp->pedge_vtx     [i_part] = edge_vtx;
  extrp->cell_ln_to_gn [i_part] = cell_ln_to_gn;
  extrp->face_ln_to_gn [i_part] = face_ln_to_gn;
  extrp->edge_ln_to_gn [i_part] = edge_ln_to_gn;
  extrp->vtx_ln_to_gn  [i_part] = vtx_ln_to_gn;
  extrp->pface_vtx_idx [i_part] = face_vtx_idx;
  extrp->pface_vtx     [i_part] = face_vtx;
  extrp->pvtx_coord    [i_part] = vtx_coord;
}

void
PDM_extract_part_selected_lnum_set
(
  PDM_extract_part_t       *extrp,
  int                       i_part,
  int                       n_extract,
  int                      *extract_lnum
)
{
  extrp->n_extract   [i_part] = n_extract;
  extrp->extract_lnum[i_part] = extract_lnum;
}

void
PDM_extract_part_free
(
  PDM_extract_part_t  *extrp
)
{
  free(extrp->n_cell        );
  free(extrp->n_face        );
  free(extrp->n_edge        );
  free(extrp->n_vtx         );
  free(extrp->n_extract     );

  free(extrp->pcell_face    );
  free(extrp->pcell_face_idx);

  free(extrp->pface_edge    );
  free(extrp->pface_edge_idx);

  free(extrp->pface_vtx     );
  free(extrp->pface_vtx_idx );
  free(extrp->pedge_vtx     );
  free(extrp->extract_lnum  );

  free(extrp->cell_ln_to_gn );
  free(extrp->face_ln_to_gn );
  free(extrp->edge_ln_to_gn );
  free(extrp->vtx_ln_to_gn  );

  free(extrp->pvtx_coord     );

  free(extrp);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
