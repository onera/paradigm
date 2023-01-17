
/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_morton.h"
#include "pdm_array.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_mpi.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_logging.h"
#include "pdm_distant_neighbor.h"
#include "pdm_part_connectivity_transform.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mesh_interpolate_priv.h"
#include "pdm_mesh_interpolate.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
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


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a structure that compute a global mean
 *
 * \param [in]   n_part       Number of local partitions
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Pointer to \ref PDM_mesh_interpolate object
 */
PDM_mesh_interpolate_t*
PDM_mesh_interpolate_create
(
 const int            n_domain,
 const int           *n_part,
 const int           *n_group,
 const int            interp_kind,
 const PDM_MPI_Comm   comm
)
{
  PDM_mesh_interpolate_t* mi = malloc(sizeof(PDM_mesh_interpolate_t));

  mi->n_domain    = n_domain;
  mi->n_part      = malloc( n_domain * sizeof(int)); // Make a copy to avoid pb in cython
  for(int i = 0; i < mi->n_domain; ++i) {
    mi->n_part[i] = n_part[i];
  }
  mi->comm        = comm;

  mi->n_part_idx    = (int * ) malloc( (n_domain + 1) * sizeof(int));
  mi->n_part_g_idx  = (int * ) malloc( (n_domain + 1) * sizeof(int));
  mi->parts = malloc(n_domain * sizeof(_part_t *));


  mi->n_part_idx  [0] = 0;
  mi->n_part_g_idx[0] = 0;
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    mi->parts[i_domain] = malloc( n_part[i_domain] * sizeof(_part_t));
    mi->n_part_idx[i_domain+1] = mi->n_part_idx[i_domain] + mi->n_part[i_domain];

    int n_part_l = n_part[i_domain];
    int n_part_g = -100;
    PDM_MPI_Allreduce(&n_part_l, &n_part_g, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
    mi->n_part_g_idx[i_domain+1] = mi->n_part_idx[i_domain] + n_part_g;

  }

  /* Graph comm  */
  mi->entity_part_bound_proc_idx = malloc( PDM_MESH_ENTITY_MAX * sizeof(int ***));
  mi->entity_part_bound_part_idx = malloc( PDM_MESH_ENTITY_MAX * sizeof(int ***));
  mi->entity_part_bound          = malloc( PDM_MESH_ENTITY_MAX * sizeof(int ***));
  mi->graph_comm_is_defined      = malloc( PDM_MESH_ENTITY_MAX * sizeof(int *  ));

  for(int i_kind = 0; i_kind < PDM_MESH_ENTITY_MAX; ++i_kind) {
    mi->entity_part_bound_proc_idx[i_kind] = malloc(n_domain * sizeof(int **));
    mi->entity_part_bound_part_idx[i_kind] = malloc(n_domain * sizeof(int **));
    mi->entity_part_bound         [i_kind] = malloc(n_domain * sizeof(int **));
    mi->graph_comm_is_defined     [i_kind] = 0;
    for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
      mi->entity_part_bound_proc_idx[i_kind][i_domain] = malloc(n_part[i_domain] * sizeof(int *));
      mi->entity_part_bound_part_idx[i_kind][i_domain] = malloc(n_part[i_domain] * sizeof(int *));
      mi->entity_part_bound         [i_kind][i_domain] = malloc(n_part[i_domain] * sizeof(int *));
      for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
        mi->entity_part_bound_proc_idx[i_kind][i_domain][i_part] = NULL;
        mi->entity_part_bound_part_idx[i_kind][i_domain][i_part] = NULL;
        mi->entity_part_bound         [i_kind][i_domain][i_part] = NULL;
      }
    }
  }

  mi->pdi = NULL;



  return mi;
}

void
PDM_mesh_interpolate_compute
(
  PDM_mesh_interpolate_t *mi
)
{

  /*
   * Compute graph comm from gnum with PDM_part_generate_entity_graph_comm is not provided
   */
  if(mi->graph_comm_is_defined[PDM_MESH_ENTITY_VERTEX] == 0) {
    abort();
    // TODO : compute graph comm from other graph comm OR recumpute from gnum : PDM_part_generate_entity_graph_comm
  }


  /* Deduce graph with all graphe inside same domain and between domain */
  int ***vtx_part_bound_proc_idx = mi->entity_part_bound_proc_idx[PDM_MESH_ENTITY_VERTEX];
  int ***vtx_part_bound_part_idx = mi->entity_part_bound_part_idx[PDM_MESH_ENTITY_VERTEX];
  int ***vtx_part_bound          = mi->entity_part_bound         [PDM_MESH_ENTITY_VERTEX];

  /* Si multidomain on fait un shift et tt roule */
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
    n_part_loc_all_domain += mi->n_part[i_domain];
  }


  int         **pdi_neighbor_idx         = NULL;
  int         **pdi_neighbor             = NULL;
  int           n_composed_interface     = 0;
  int          *composed_interface_idx   = NULL;
  int          *composed_interface       = NULL;
  PDM_g_num_t  *composed_ln_to_gn_sorted = NULL;

  if(mi->pdi != NULL) {
    int **pn_vtx = malloc(mi->n_domain * sizeof(int *));
    for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
      pn_vtx[i_domain] = malloc(mi->n_part[i_domain] * sizeof(int));
      for(int i_part = 0; i_part < mi->n_part[i_domain]; ++i_part) {
        pn_vtx[i_domain][i_part] = mi->parts[i_domain][i_part].n_vtx;
      }
    }

    PDM_part_domain_interface_as_graph(mi->pdi,
                                       PDM_BOUND_TYPE_VTX,
                                       pn_vtx,
                                       NULL,
                                       &pdi_neighbor_idx,
                                       &pdi_neighbor,
                                       &n_composed_interface,
                                       &composed_interface_idx,
                                       &composed_interface,
                                       &composed_ln_to_gn_sorted);
    for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
      free(pn_vtx[i_domain]);
    }
    free(pn_vtx);
  }


  int shift_part   = 0;
  int shift_part_g = 0;

  int **neighbor_idx       = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int *));
  int **neighbor_desc      = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int *));
  int **neighbor_interface = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int *));
  int *pn_vtx              = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int  ));
  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {

    int n_part_total = mi->n_part_g_idx[i_domain+1] - mi->n_part_g_idx[i_domain];

    /* First loop to count */
    for(int i_part = 0; i_part < mi->n_part[i_domain]; ++i_part) {

      int n_vtx = mi->parts[i_domain][i_part].n_vtx;
      pn_vtx[i_part+shift_part] = n_vtx;

      neighbor_idx[i_part+shift_part] = malloc((n_vtx+1) * sizeof(int));
      int* _neighbor_idx   = neighbor_idx[i_part+shift_part];
      int* _vtx_part_bound = vtx_part_bound[i_domain][i_part];

      int* _neighbor_n = PDM_array_zeros_int(n_vtx);

      int n_part_entity_bound_tot = vtx_part_bound_part_idx[i_domain][i_part][n_part_total];
      for(int idx_entity = 0; idx_entity < n_part_entity_bound_tot; ++idx_entity) {
        int i_entity = _vtx_part_bound[4*idx_entity]-1;
        _neighbor_n[i_entity] += 1;
      }

      /*
       * Add comming from interface
       */
      if(pdi_neighbor_idx != NULL) {
        for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
          _neighbor_n[i_entity] += pdi_neighbor_idx[i_part+shift_part][i_entity+1] - pdi_neighbor_idx[i_part+shift_part][i_entity];
        }
      }

      /* Compute index */
      _neighbor_idx[0] = 0;
      for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
        _neighbor_idx[i_entity+1] = _neighbor_idx[i_entity] + _neighbor_n[i_entity];
        _neighbor_n[i_entity] = 0;
      }

      PDM_log_trace_array_int(_neighbor_idx, n_vtx, "_neighbor_idx ::");
      PDM_log_trace_array_int(_neighbor_n, n_vtx, "_neighbor_n ::");


      neighbor_desc     [i_part+shift_part] = (int *) malloc( 3 * _neighbor_idx[n_vtx] * sizeof(int) );
      neighbor_interface[i_part+shift_part] = (int *) malloc(     _neighbor_idx[n_vtx] * sizeof(int) );
      int* _neighbor_desc      = neighbor_desc     [i_part+shift_part];
      int* _neighbor_interface = neighbor_interface[i_part+shift_part];

      /* Fill */
      for(int idx_entity = 0; idx_entity < vtx_part_bound_part_idx[i_domain][i_part][n_part_total]; ++idx_entity) {
        int i_entity = _vtx_part_bound[4*idx_entity]-1;
        int idx_write = _neighbor_idx[i_entity] + _neighbor_n[i_entity]++;
        _neighbor_desc[3*idx_write  ] = _vtx_part_bound[4*idx_entity+1];                // i_proc_opp;
        _neighbor_desc[3*idx_write+1] = _vtx_part_bound[4*idx_entity+2]+shift_part_g-1; // i_part_opp
        _neighbor_desc[3*idx_write+2] = _vtx_part_bound[4*idx_entity+3]-1;              // i_entity_opp
        _neighbor_interface[idx_write] = -40000;
      }

      if(pdi_neighbor_idx != NULL) {
        for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
          for(int idx = pdi_neighbor_idx[i_part+shift_part][i_entity]; idx < pdi_neighbor_idx[i_part+shift_part][i_entity+1]; ++idx) {
            int idx_write = _neighbor_idx[i_entity] + _neighbor_n[i_entity]++;
            _neighbor_desc[3*idx_write  ]  = pdi_neighbor[i_part+shift_part][4*idx  ];
            _neighbor_desc[3*idx_write+1]  = pdi_neighbor[i_part+shift_part][4*idx+1];
            _neighbor_desc[3*idx_write+2]  = pdi_neighbor[i_part+shift_part][4*idx+2];
            _neighbor_interface[idx_write] = pdi_neighbor[i_part+shift_part][4*idx+3];
          }
        }
      }

      free(_neighbor_n);
    }

    shift_part   += mi->n_part              [i_domain];
    shift_part_g += n_part_total;
  }

  /*
   * Create vtx_cell
   */
  shift_part = 0;
  int **pvtx_cell_n   = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int *));
  int **pvtx_cell_idx = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int *));
  int **pvtx_cell     = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int *));
  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < mi->n_part[i_domain]; ++i_part) {

      /* Compute cell_edge */
      int* cell_vtx_idx = NULL;
      int* cell_vtx     = NULL;
      PDM_combine_connectivity(mi->parts[i_domain][i_part].n_cell,
                               mi->parts[i_domain][i_part].cell_face_idx,
                               mi->parts[i_domain][i_part].cell_face,
                               mi->parts[i_domain][i_part].face_vtx_idx,
                               mi->parts[i_domain][i_part].face_vtx,
                               &cell_vtx_idx,
                               &cell_vtx);

      PDM_connectivity_transpose(mi->parts[i_domain][i_part].n_cell,
                                 mi->parts[i_domain][i_part].n_vtx,
                                 cell_vtx_idx,
                                 cell_vtx,
                                 &pvtx_cell_idx[i_part+shift_part],
                                 &pvtx_cell    [i_part+shift_part]);

      free(cell_vtx_idx);
      free(cell_vtx);

      pvtx_cell_n[i_part+shift_part] = malloc(mi->parts[i_domain][i_part].n_vtx * sizeof(int));
      for(int i_vtx = 0; i_vtx < mi->parts[i_domain][i_part].n_vtx; ++i_vtx) {
        pvtx_cell_n[i_part+shift_part][i_vtx] = pvtx_cell_idx[i_part+shift_part][i_vtx+1] - pvtx_cell_idx[i_part+shift_part][i_vtx];
      }

    }
    shift_part   += mi->n_part              [i_domain];
  }


  /*
   * Create distant_neighbor
   */
  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(mi->comm,
                                                           n_part_loc_all_domain,
                                                           pn_vtx,
                                                           neighbor_idx,
                                                           neighbor_desc);

  int **pvtx_cell_opp_n = NULL;
  int **pvtx_cell_opp   = NULL;
  PDM_distant_neighbor_exch(dn,
                            sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            pvtx_cell_n,
                  (void **) pvtx_cell,
                           &pvtx_cell_opp_n,
                 (void ***)&pvtx_cell_opp);

  PDM_distant_neighbor_free(dn);


  shift_part = 0;
  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < mi->n_part[i_domain]; ++i_part) {
      free(neighbor_idx      [i_part+shift_part]);
      free(neighbor_desc     [i_part+shift_part]);
      free(neighbor_interface[i_part+shift_part]);

      if(mi->pdi != NULL) {
        free(pdi_neighbor_idx[i_part+shift_part]);
        free(pdi_neighbor    [i_part+shift_part]);
      }

      free(pvtx_cell_n  [i_part+shift_part]);
      free(pvtx_cell_idx[i_part+shift_part]);
      free(pvtx_cell    [i_part+shift_part]);

      free(pvtx_cell_opp_n[i_part+shift_part]);
      free(pvtx_cell_opp  [i_part+shift_part]);
    }
    shift_part   += mi->n_part              [i_domain];
  }
  free(neighbor_idx);
  free(neighbor_desc);
  free(neighbor_interface);
  free(pn_vtx);
  free(pvtx_cell_n  );
  free(pvtx_cell_idx);
  free(pvtx_cell    );
  free(pvtx_cell_opp_n);
  free(pvtx_cell_opp  );

  if(mi->pdi != NULL) {
    free(pdi_neighbor_idx);
    free(pdi_neighbor    );

    free(composed_interface_idx  );
    free(composed_interface      );
    free(composed_ln_to_gn_sorted);
  }

  /* Compute weight */
  /* Create protocol only between join - Distant neighbor */
}



void
PDM_mesh_interpolate_part_set
(
  PDM_mesh_interpolate_t   *mi,
  int                       i_domain,
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
  mi->parts[i_domain][i_part].n_cell            = n_cell;
  mi->parts[i_domain][i_part].n_face            = n_face;
  mi->parts[i_domain][i_part].n_edge            = n_edge;
  mi->parts[i_domain][i_part].n_vtx             = n_vtx;

  mi->parts[i_domain][i_part].cell_face_idx = cell_face_idx;
  mi->parts[i_domain][i_part].cell_face     = cell_face;

  mi->parts[i_domain][i_part].face_edge_idx = face_edge_idx;
  mi->parts[i_domain][i_part].face_edge     = face_edge;

  mi->parts[i_domain][i_part].face_vtx_idx  = face_vtx_idx;
  mi->parts[i_domain][i_part].face_vtx      = face_vtx;

  mi->parts[i_domain][i_part].edge_vtx      = edge_vtx;

  mi->parts[i_domain][i_part].cell_ln_to_gn = cell_ln_to_gn;
  mi->parts[i_domain][i_part].face_ln_to_gn = face_ln_to_gn;
  mi->parts[i_domain][i_part].edge_ln_to_gn = edge_ln_to_gn;
  mi->parts[i_domain][i_part].vtx_ln_to_gn  = vtx_ln_to_gn;

  mi->parts[i_domain][i_part].vtx = vtx_coord;
}


void
PDM_mesh_interpolate_graph_comm_set
(
  PDM_mesh_interpolate_t   *mi,
  int                       i_domain,
  int                       i_part,
  PDM_mesh_entities_t       mesh_entity,
  int                      *entity_part_bound_proc_idx,
  int                      *entity_part_bound_part_idx,
  int                      *entity_part_bound
)
{
  mi->entity_part_bound_proc_idx[mesh_entity][i_domain][i_part] = entity_part_bound_proc_idx;
  mi->entity_part_bound_part_idx[mesh_entity][i_domain][i_part] = entity_part_bound_part_idx;
  mi->entity_part_bound         [mesh_entity][i_domain][i_part] = entity_part_bound;
  mi->graph_comm_is_defined[mesh_entity] = 1;
}


void
PDM_mesh_interpolate_part_domain_interface_shared_set
(
  PDM_mesh_interpolate_t      *mi,
  PDM_part_domain_interface_t *pdi
)
{
  mi->pdi = pdi;
}


void
PDM_mesh_interpolate_free
(
 PDM_mesh_interpolate_t *mi
)
{

  for(int i_kind = 0; i_kind < PDM_MESH_ENTITY_MAX; ++i_kind) {
    mi->graph_comm_is_defined     [i_kind] = 0;
    for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
      free(mi->entity_part_bound_proc_idx[i_kind][i_domain]);
      free(mi->entity_part_bound_part_idx[i_kind][i_domain]);
      free(mi->entity_part_bound         [i_kind][i_domain]);
    }
    free(mi->entity_part_bound_proc_idx[i_kind]);
    free(mi->entity_part_bound_part_idx[i_kind]);
    free(mi->entity_part_bound         [i_kind]);
  }
  free(mi->entity_part_bound_proc_idx);
  free(mi->entity_part_bound_part_idx);
  free(mi->entity_part_bound         );
  free(mi->graph_comm_is_defined     );


  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
    free(mi->parts[i_domain]);
  }
  free(mi->parts);
  free(mi->n_part_idx);
  free(mi->n_part_g_idx);
  free(mi->n_part);


  free(mi);
}


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
