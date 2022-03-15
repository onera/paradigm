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
      entity_center[i_part][3*idx_face  ] = 0.;
      entity_center[i_part][3*idx_face+1] = 0.;
      entity_center[i_part][3*idx_face+2] = 0.;

      double inv = 1./((double) _pface_edge_idx[i_face+1] - _pface_edge_idx[i_face]);

      for(int idx_edge = _pface_edge_idx[i_face]; idx_edge < _pface_edge_idx[i_face+1]; ++idx_edge) {
        int i_edge = PDM_ABS(_pface_edge[idx_edge])-1;
        int i_vtx1 = _pedge_vtx[2*i_edge  ] - 1;
        int i_vtx2 = _pedge_vtx[2*i_edge+1] - 1;

        entity_center[i_part][3*idx_face  ] += 0.5 * (_pvtx_coord[i_vtx1] + _pvtx_coord[i_vtx2]);
        entity_center[i_part][3*idx_face+1] += 0.5 * (_pvtx_coord[i_vtx1] + _pvtx_coord[i_vtx2]);
        entity_center[i_part][3*idx_face+2] += 0.5 * (_pvtx_coord[i_vtx1] + _pvtx_coord[i_vtx2]);

      }
      entity_center[i_part][3*idx_face  ] = entity_center[i_part][3*idx_face  ] * inv;
      entity_center[i_part][3*idx_face+1] = entity_center[i_part][3*idx_face+1] * inv;
      entity_center[i_part][3*idx_face+2] = entity_center[i_part][3*idx_face+2] * inv;
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

      int    *_pcell_face     = pcell_face    [i_part];
      int    *_pcell_face_idx = pcell_face_idx[i_part];
      int    *_pface_vtx      = pface_vtx     [i_part];
      int    *_pface_vtx_idx  = pface_vtx_idx [i_part];
      double *_pvtx_coord     = pvtx_coord    [i_part];

      // PDM_log_trace_array_int(extract_lnum[i_part], n_extract[i_part], "extract_lnum ::");
      for(int idx_cell = 0; idx_cell < n_extract[i_part]; ++idx_cell) {
        int i_cell = extract_lnum[i_part][idx_cell];

        entity_center[i_part][3*idx_cell  ] = 0.;
        entity_center[i_part][3*idx_cell+1] = 0.;
        entity_center[i_part][3*idx_cell+2] = 0.;

        double inv = 1./((double)  _pcell_face_idx[idx_cell+1] - _pcell_face_idx[idx_cell]);

        double fcx = 0;
        double fcy = 0;
        double fcz = 0;
        for(int idx_face = _pcell_face_idx[i_cell]; idx_face < _pcell_face_idx[i_cell+1]; ++idx_face) {
          int i_face = PDM_ABS(_pcell_face[idx_face])-1;

          double inv2 = 1./((double)  _pface_vtx_idx[i_face+1] - _pface_vtx_idx[i_face]);

          for(int idx_vtx = _pface_vtx_idx[i_face]; idx_vtx < _pface_vtx_idx[i_face+1]; ++idx_vtx) {
            int i_vtx = _pface_vtx[idx_vtx]-1;
            fcx += _pvtx_coord[i_vtx];
            fcy += _pvtx_coord[i_vtx];
            fcz += _pvtx_coord[i_vtx];
          }
          fcx = fcx * inv2;
          fcy = fcy * inv2;
          fcz = fcz * inv2;

          entity_center[i_part][3*idx_cell  ] += fcx;
          entity_center[i_part][3*idx_cell+1] += fcy;
          entity_center[i_part][3*idx_cell+2] += fcz;
        }

        entity_center[i_part][3*idx_cell  ] = entity_center[i_part][3*idx_cell  ] * inv;
        entity_center[i_part][3*idx_cell+1] = entity_center[i_part][3*idx_cell+1] * inv;
        entity_center[i_part][3*idx_cell+2] = entity_center[i_part][3*idx_cell+2] * inv;
      } /* End cell */
    }
  } else if( from_edge == 1) {
    for(int i_part = 0; i_part < n_part_in; ++i_part) {
      entity_center[i_part] = (double *) malloc(3 * n_extract[i_part] * sizeof(double));

      int    *_pcell_face     = pcell_face    [i_part];
      int    *_pcell_face_idx = pcell_face_idx[i_part];
      int    *_pface_edge     = pface_edge    [i_part];
      int    *_pface_edge_idx = pface_edge_idx[i_part];
      int    *_pedge_vtx      = pedge_vtx     [i_part];
      double *_pvtx_coord     = pvtx_coord    [i_part];

      for(int idx_cell = 0; idx_cell < n_extract[i_part]; ++idx_cell) {
        int i_cell = extract_lnum[i_part][idx_cell];

        entity_center[i_part][3*idx_cell  ] = 0.;
        entity_center[i_part][3*idx_cell+1] = 0.;
        entity_center[i_part][3*idx_cell+2] = 0.;

        double inv = 1./((double)  _pcell_face_idx[idx_cell+1] - _pcell_face_idx[idx_cell]);

        double fcx = 0;
        double fcy = 0;
        double fcz = 0;
        for(int idx_face = _pcell_face_idx[i_cell]; idx_face < _pcell_face_idx[i_cell+1]; ++idx_face) {
          int i_face = PDM_ABS(_pcell_face[idx_face])-1;

          double inv2 = 1./((double)  _pface_edge_idx[i_face+1] - _pface_edge_idx[i_face]);

          for(int idx_edge = _pface_edge_idx[i_face]; idx_edge < _pface_edge_idx[i_face+1]; ++idx_edge) {
            int i_edge = PDM_ABS(_pface_edge[idx_edge])-1;
            int i_vtx1 = _pedge_vtx[2*i_edge  ] - 1;
            int i_vtx2 = _pedge_vtx[2*i_edge+1] - 1;
            fcx += 0.5 * (_pvtx_coord[i_vtx1] + _pvtx_coord[i_vtx2]);
            fcy += 0.5 * (_pvtx_coord[i_vtx1] + _pvtx_coord[i_vtx2]);
            fcz += 0.5 * (_pvtx_coord[i_vtx1] + _pvtx_coord[i_vtx2]);
          }
          fcx = fcx * inv2;
          fcy = fcy * inv2;
          fcz = fcz * inv2;

          entity_center[i_part][3*idx_cell  ] += fcx;
          entity_center[i_part][3*idx_cell+1] += fcy;
          entity_center[i_part][3*idx_cell+2] += fcz;
        }

        entity_center[i_part][3*idx_cell  ] = entity_center[i_part][3*idx_cell  ] * inv;
        entity_center[i_part][3*idx_cell+1] = entity_center[i_part][3*idx_cell+1] * inv;
        entity_center[i_part][3*idx_cell+2] = entity_center[i_part][3*idx_cell+2] * inv;
      } /* End cell */
    }
  }

  *cell_center = entity_center;
}



static
void
extract_entity1_entity2
(
int            n_part,
int           *n_entity1,
int           *n_entity2,
int           *n_extract_entity1,
int          **extract_entity1_lnum,
int          **entity1_entity2_idx,
int          **entity1_entity2,
PDM_g_num_t  **entity2_ln_to_gn,
int         ***extract_entity1_entity2_n,
PDM_g_num_t ***extract_entity1_entity2,
int          **n_extract_entity2,
int         ***extract_entity2_lnum,
int         ***idx_visited
)
{
  int         **_extract_entity1_entity2_n = malloc(n_part * sizeof(int         *));
  int         **_extract_entity2_lnum      = malloc(n_part * sizeof(int         *));
  PDM_g_num_t **_extract_entity1_entity2   = malloc(n_part * sizeof(PDM_g_num_t *));
  int          *_n_extract_entity2         = malloc(n_part * sizeof(int          ));
  int         **_idx_visited               = malloc(n_part * sizeof(int         *));

  for(int i_part = 0; i_part < n_part; ++i_part) {
    int n_extract = n_extract_entity1[i_part];
    int         *_pentity1_entity2     = entity1_entity2    [i_part];
    int         *_pentity1_entity2_idx = entity1_entity2_idx[i_part];
    PDM_g_num_t *_pentity2_ln_to_gn    = entity2_ln_to_gn   [i_part];
    int          _pn_entity1           = n_entity1          [i_part];
    int          _pn_entity2           = n_entity2          [i_part];

    _extract_entity1_entity2  [i_part] = (PDM_g_num_t  *) malloc( _pentity1_entity2_idx[_pn_entity1] * sizeof(PDM_g_num_t));
    _extract_entity1_entity2_n[i_part] = (int          *) malloc( (n_extract)     * sizeof(int        ));
    _extract_entity2_lnum     [i_part] = (int          *) malloc( (_pn_entity2+1) * sizeof(int        ));
    _idx_visited              [i_part] = (int          *) malloc( _pn_entity2     * sizeof(int        ));

    int         *is_visited = (int *) malloc( _pn_entity2 * sizeof(int));
    for(int i = 0; i < _pn_entity2; ++i) {
      is_visited[i] = 0;
    }

    int idx_write = 0;
    _n_extract_entity2[i_part] = 0;
    _extract_entity1_entity2_n[i_part][0] = 0;
    for(int idx_entity = 0; idx_entity < n_extract; ++idx_entity) {
      int i_entity = extract_entity1_lnum[i_part][idx_entity];

      int n_tmp = _pentity1_entity2_idx[i_entity+1] - _pentity1_entity2_idx[i_entity];
      _extract_entity1_entity2_n[i_part][idx_entity] = n_tmp;

      for(int idx_entity2 = _pentity1_entity2_idx[i_entity]; idx_entity2 < _pentity1_entity2_idx[i_entity+1]; ++idx_entity2) {
        int i_entity2 = PDM_ABS(_pentity1_entity2[idx_entity2])-1;
        _extract_entity1_entity2[i_part][idx_write++] = _pentity2_ln_to_gn[i_entity2];

        if(is_visited[i_entity2] == 0) {
          int idx = _n_extract_entity2[i_part]++;
          _extract_entity2_lnum[i_part][idx] = i_entity2;
          is_visited[i_entity2] = 1;
          _idx_visited[i_part][idx] = idx_write-1;
        }
      }
    }
    _extract_entity1_entity2  [i_part] = realloc(_extract_entity1_entity2[i_part], idx_write                  * sizeof(PDM_g_num_t));
    _extract_entity2_lnum     [i_part] = realloc(_extract_entity2_lnum   [i_part], _n_extract_entity2[i_part] * sizeof(int        ));
    _idx_visited              [i_part] = realloc(_idx_visited            [i_part], _n_extract_entity2[i_part] * sizeof(int        ));

    free(is_visited);
  }

  *n_extract_entity2           = _n_extract_entity2;
  *extract_entity2_lnum        = _extract_entity2_lnum;
  *extract_entity1_entity2_n   = _extract_entity1_entity2_n;
  *extract_entity1_entity2     = _extract_entity1_entity2;
  *idx_visited                 = _idx_visited;
}


static
void
extract_and_renum_entity1_entity2
(
PDM_MPI_Comm           comm,
PDM_part_to_block_t   *ptb_entity1,
int                    n_part,
int                   *n_entity1,
int                   *n_entity2,
int                   *n_extract_entity1,
int                  **extract_entity1_lnum,
int                  **entity1_entity2_idx,
int                  **entity1_entity2,
PDM_g_num_t          **entity2_ln_to_gn,
int                 ***selected_entity1_entity2_n,
PDM_g_num_t         ***selected_entity1_entity2,
int                  **n_extract_entity2,
int                 ***extract_entity2_lnum
)
{
  int **extract_entity1_entity2_idx = NULL;
  extract_entity1_entity2(n_part,
                          n_entity1,
                          n_entity2,
                          n_extract_entity1,
                          extract_entity1_lnum,
                          entity1_entity2_idx,
                          entity1_entity2,
                          entity2_ln_to_gn,
                          selected_entity1_entity2_n,
                          selected_entity1_entity2,
                          n_extract_entity2,
                          extract_entity2_lnum,
                          &extract_entity1_entity2_idx);

  int         **_selected_entity1_entity2_n = *selected_entity1_entity2_n;
  PDM_g_num_t **_selected_entity1_entity2   = *selected_entity1_entity2;

  // PDM_gen_gnum_t* gnum_extract_entity2 = PDM_gnum_create(3, 1, PDM_FALSE,
  //                                                        1.e-6,
  //                                                        comm,
  //                                                        PDM_OWNERSHIP_KEEP);

  // PDM_gnum_set_from_parents(gnum_extract_entity2, 0, dn_entity1_entity2, dequi_entity1_entity2);
  // PDM_gnum_compute(gnum_extract_entity2);

  // PDM_g_num_t *child_entity1_entity2_gnum = PDM_gnum_get(gnum_extract_entity2, 0);

  int *_n_extract_entity2 = *n_extract_entity2;

  PDM_gen_gnum_t* gnum_extract_entity2 = PDM_gnum_create(3,
                                                         n_part,
                                                         PDM_FALSE,
                                                         1.e-6,
                                                         comm,
                                                         PDM_OWNERSHIP_USER);
  //
  int *n_connect_tot = (int *) malloc( n_part * sizeof(int));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    n_connect_tot[i_part] = 0;
    for(int i = 0; i < n_extract_entity1[i_part]; ++i) {
      n_connect_tot[i_part] += _selected_entity1_entity2_n[i_part][i];
    }
    // TODO --> compute idx

    PDM_gnum_set_from_parents(gnum_extract_entity2, i_part, n_connect_tot[i_part], _selected_entity1_entity2[i_part]);
  }

  PDM_gnum_compute(gnum_extract_entity2);

  PDM_g_num_t **_child_entity1_entity2           = malloc(n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t **_child_entity2_ln_to_gn          = malloc(n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t **_extract_parent_entity2_ln_to_gn = malloc(n_part * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    _child_entity1_entity2[i_part] = PDM_gnum_get(gnum_extract_entity2, i_part);

    if(0 == 1) {
      PDM_log_trace_array_long(_child_entity1_entity2[i_part], n_connect_tot[i_part], "_child_entity1_entity2 :: ");
      PDM_log_trace_array_long(extract_entity1_entity2_idx[i_part], _n_extract_entity2[i_part], "extract_entity1_entity2_idx :: ");
    }

    // Reforme a temporary ln_to_gn to setup recurence
    _child_entity2_ln_to_gn         [i_part] = (PDM_g_num_t *) malloc( _n_extract_entity2[i_part] * sizeof(PDM_g_num_t));
    _extract_parent_entity2_ln_to_gn[i_part] = (PDM_g_num_t *) malloc( _n_extract_entity2[i_part] * sizeof(PDM_g_num_t));
    for(int i = 0; i < _n_extract_entity2[i_part]; ++i) {
      int idx = extract_entity1_entity2_idx[i_part][i];
      _child_entity2_ln_to_gn         [i_part][i] = _child_entity1_entity2   [i_part][idx];
      _extract_parent_entity2_ln_to_gn[i_part][i] = _selected_entity1_entity2[i_part][idx]; // To exchange to keep link with original part
    }

    free(extract_entity1_entity2_idx[i_part]);
  }
  free(extract_entity1_entity2_idx);
  free(n_connect_tot);
  PDM_gnum_free(gnum_extract_entity2);

  int         *dequi_entity1_entity2_n = NULL;
  PDM_g_num_t *dequi_entity1_entity2   = NULL;
  int dn_entity1_entity2 = PDM_part_to_block_exch(ptb_entity1,
                                                  sizeof(PDM_g_num_t),
                                                  PDM_STRIDE_VAR_INTERLACED,
                                                  -1,
                                                  _selected_entity1_entity2_n,
                                   (void **)      _child_entity1_entity2,
                                                  &dequi_entity1_entity2_n,
                                   (void **)      &dequi_entity1_entity2);


  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(_child_entity1_entity2 [i_part]);
  }
  free(_child_entity1_entity2 );


  PDM_part_to_block_t *ptb_entity2_equi = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                   PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                                   1.,
                                                                   _child_entity2_ln_to_gn,
                                                                   NULL,
                                                                   _n_extract_entity2,
                                                                   n_part,
                                                                   comm);


  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(_child_entity2_ln_to_gn[i_part]);
  }
  free(_child_entity2_ln_to_gn);

  int          dn_entity2_equi        = PDM_part_to_block_n_elt_block_get(ptb_entity2_equi);
  PDM_g_num_t *dextract_entity2_gnum  = PDM_part_to_block_block_gnum_get (ptb_entity2_equi);

  if(0 == 1) {
    PDM_log_trace_array_long(dextract_entity2_gnum, dn_entity2_equi, "dextract_entity2_gnum :: ");
  }



  /*
   *  Exchange parent to keep link
   */

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(_extract_parent_entity2_ln_to_gn[i_part]);
  }
  free(_extract_parent_entity2_ln_to_gn);



  PDM_part_to_block_free(ptb_entity2_equi);


  free(dequi_entity1_entity2_n);
  free(dequi_entity1_entity2);



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

  PDM_gen_gnum_t* gnum_extract = PDM_gnum_create(3, extrp->n_part_in, PDM_FALSE,
                                                 1.e-6,
                                                 extrp->comm,
                                                 PDM_OWNERSHIP_KEEP);
;

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

    for (int i_part = 0; i_part < extrp->n_part_in; i_part++){
      PDM_gnum_set_from_coords(gnum_extract, i_part, extrp->n_extract[i_part], entity_center[i_part], NULL);
    }
  } else {

    for (int i_part = 0; i_part < extrp->n_part_in; i_part++){
      PDM_gnum_set_from_parents(gnum_extract, i_part, extrp->n_extract[i_part], entity_extract_g_num[i_part]);
    }

  }

  /*
   * Global numering computation
   */
  PDM_g_num_t **child_selected_g_num = (PDM_g_num_t **) malloc( extrp->n_part_in * sizeof(PDM_g_num_t *));
  PDM_gnum_compute(gnum_extract);

  for (int i_part = 0; i_part < extrp->n_part_in; i_part++){
    child_selected_g_num[i_part] = PDM_gnum_get(gnum_extract, i_part);
    // PDM_log_trace_array_long(child_selected_g_num[i_part], extrp->n_extract[i_part], "child_selected_g_num : ");
  }

  /*
   *  Remake equilibrate block -> Block is not partial because we use child_gnum
   */
  double **weight = (double **) malloc( extrp->n_part_in * sizeof(double *));
  for (int i_part = 0; i_part < extrp->n_part_in; i_part++){
    weight        [i_part] = malloc(extrp->n_extract[i_part] * sizeof(double));
    for(int i = 0; i < extrp->n_extract[i_part]; ++i) {
      weight[i_part][i] = 1.;
    }
  }

  PDM_part_to_block_t *ptb_equi = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                           1.,
                                                           child_selected_g_num,
                                                           weight,
                                                           extrp->n_extract,
                                                           extrp->n_part_in,
                                                           extrp->comm);


  for (int i_part = 0; i_part < extrp->n_part_in; i_part++){
    free(weight[i_part]);
  }
  free(weight);
  free(child_selected_g_num);
  PDM_gnum_free(gnum_extract);

  /*
   *  exch parent_num + cell_face
   */
  // PDM_part_to_block_exch(ptb_equi,
  //                        sizeof(PDM_g_num_t),
  //                        PDM_STRIDE_VAR_INTERLACED,
  //                        -1,
  //                        pentity1_entity2_n,
  //                        pentity1_entity2,
  //                        &dequi_entity1_entity2_n,
  //                        &dequi_entity1_entity2);

  /*
   * On ramene avec un part_to_block les faces aussi
   * Mais on doit faire un dconnectivity_to_p_connectivty après
   */


  int          dn_entity_equi = PDM_part_to_block_n_elt_block_get(ptb_equi);
  PDM_g_num_t *dextract_gnum  = PDM_part_to_block_block_gnum_get (ptb_equi);

  if(extrp->split_dual_method == PDM_SPLIT_DUAL_WITH_HILBERT) {

    /*
     *
     */


    for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
      free(entity_center[i_part]);
    }
    free(entity_center);
  } else {
    /*
     * Si scotch ou metis, on doit calculer le dcell_face + pdm_gnum (pour avoir les faces child)
     *   Puis dconnectiviy_transpose puis combine
     */

    abort();

  }

  /*
   * Extraction des connectivités
   */
  if(extrp->dim == 3) {


    PDM_g_num_t **selected_cell_face    = NULL;
    int         **selected_cell_face_n  = NULL;
    int          *n_extract_face        = NULL;
    int         **extract_face_lnum     = NULL;
    int         **idx_face_in_cell_face = NULL;
    extract_and_renum_entity1_entity2(extrp->comm,
                                      ptb_equi,
                                      extrp->n_part_in,
                                      extrp->n_cell,
                                      extrp->n_face,
                                      extrp->n_extract,
                                      extrp->extract_lnum,
                                      extrp->pcell_face_idx,
                                      extrp->pcell_face,
                                      extrp->face_ln_to_gn,
                                      &selected_cell_face_n,
                                      &selected_cell_face,
                                      &n_extract_face,
                                      &extract_face_lnum);

    for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
      free(selected_cell_face_n[i_part]);
      free(selected_cell_face  [i_part]);
      free(extract_face_lnum   [i_part]);
    }
    free(selected_cell_face_n);
    free(selected_cell_face  );
    free(extract_face_lnum   );
    free(n_extract_face      );

  }


  PDM_part_to_block_free(ptb_equi);

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
