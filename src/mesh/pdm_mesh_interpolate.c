
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

static
void
_cell_center_3d
(
  int       n_part_in,
  int      *pn_cell,
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
      entity_center[i_part] = (double *) malloc(3 * pn_cell[i_part] * sizeof(double));

      int    *_pcell_face     = pcell_face    [i_part];
      int    *_pcell_face_idx = pcell_face_idx[i_part];
      int    *_pface_vtx      = pface_vtx     [i_part];
      int    *_pface_vtx_idx  = pface_vtx_idx [i_part];
      double *_pvtx_coord     = pvtx_coord    [i_part];

      // PDM_log_trace_array_int(extract_lnum[i_part], pn_cell[i_part], "extract_lnum ::");
      for(int i_cell = 0; i_cell < pn_cell[i_part]; ++i_cell) {

        entity_center[i_part][3*i_cell  ] = 0.;
        entity_center[i_part][3*i_cell+1] = 0.;
        entity_center[i_part][3*i_cell+2] = 0.;

        double inv = 1./((double)  _pcell_face_idx[i_cell+1] - _pcell_face_idx[i_cell]);

        double fcx = 0;
        double fcy = 0;
        double fcz = 0;
        for(int idx_face = _pcell_face_idx[i_cell]; idx_face < _pcell_face_idx[i_cell+1]; ++idx_face) {
          int i_face = PDM_ABS(_pcell_face[idx_face])-1;

          double inv2 = 1./((double)  _pface_vtx_idx[i_face+1] - _pface_vtx_idx[i_face]);

          for(int idx_vtx = _pface_vtx_idx[i_face]; idx_vtx < _pface_vtx_idx[i_face+1]; ++idx_vtx) {
            int i_vtx = _pface_vtx[idx_vtx]-1;
            fcx += _pvtx_coord[3*i_vtx  ];
            fcy += _pvtx_coord[3*i_vtx+1];
            fcz += _pvtx_coord[3*i_vtx+2];
          }
          fcx = fcx * inv2;
          fcy = fcy * inv2;
          fcz = fcz * inv2;

          entity_center[i_part][3*i_cell  ] += fcx;
          entity_center[i_part][3*i_cell+1] += fcy;
          entity_center[i_part][3*i_cell+2] += fcz;
        }

        entity_center[i_part][3*i_cell  ] = entity_center[i_part][3*i_cell  ] * inv;
        entity_center[i_part][3*i_cell+1] = entity_center[i_part][3*i_cell+1] * inv;
        entity_center[i_part][3*i_cell+2] = entity_center[i_part][3*i_cell+2] * inv;
      } /* End cell */
    }
  } else if( from_edge == 1) {
    for(int i_part = 0; i_part < n_part_in; ++i_part) {
      entity_center[i_part] = (double *) malloc(3 * pn_cell[i_part] * sizeof(double));

      int    *_pcell_face     = pcell_face    [i_part];
      int    *_pcell_face_idx = pcell_face_idx[i_part];
      int    *_pface_edge     = pface_edge    [i_part];
      int    *_pface_edge_idx = pface_edge_idx[i_part];
      int    *_pedge_vtx      = pedge_vtx     [i_part];
      double *_pvtx_coord     = pvtx_coord    [i_part];

      for(int i_cell = 0; i_cell < pn_cell[i_part]; ++i_cell) {

        entity_center[i_part][3*i_cell  ] = 0.;
        entity_center[i_part][3*i_cell+1] = 0.;
        entity_center[i_part][3*i_cell+2] = 0.;

        double inv = 1./((double)  _pcell_face_idx[i_cell+1] - _pcell_face_idx[i_cell]);

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
            fcx += 0.5 * (_pvtx_coord[3*i_vtx1  ] + _pvtx_coord[3*i_vtx2  ]);
            fcy += 0.5 * (_pvtx_coord[3*i_vtx1+1] + _pvtx_coord[3*i_vtx2+1]);
            fcz += 0.5 * (_pvtx_coord[3*i_vtx1+2] + _pvtx_coord[3*i_vtx2+2]);
          }
          fcx = fcx * inv2;
          fcy = fcy * inv2;
          fcz = fcz * inv2;

          entity_center[i_part][3*i_cell  ] += fcx;
          entity_center[i_part][3*i_cell+1] += fcy;
          entity_center[i_part][3*i_cell+2] += fcz;
        }

        entity_center[i_part][3*i_cell  ] = entity_center[i_part][3*i_cell  ] * inv;
        entity_center[i_part][3*i_cell+1] = entity_center[i_part][3*i_cell+1] * inv;
        entity_center[i_part][3*i_cell+2] = entity_center[i_part][3*i_cell+2] * inv;
      } /* End cell */
    }
  }

  *cell_center = entity_center;
}

static
void
_prepare_cell_center
(
  PDM_mesh_interpolate_t* mi
)
{
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
    n_part_loc_all_domain += mi->n_part[i_domain];
  }
  int     *pn_cell        = malloc(n_part_loc_all_domain * sizeof(int    *));
  int    **pcell_face_idx = malloc(n_part_loc_all_domain * sizeof(int    *));
  int    **pcell_face     = malloc(n_part_loc_all_domain * sizeof(int    *));
  int    **pface_edge_idx = malloc(n_part_loc_all_domain * sizeof(int    *));
  int    **pface_edge     = malloc(n_part_loc_all_domain * sizeof(int    *));
  int    **pface_vtx_idx  = malloc(n_part_loc_all_domain * sizeof(int    *));
  int    **pface_vtx      = malloc(n_part_loc_all_domain * sizeof(int    *));
  int    **pedge_vtx      = malloc(n_part_loc_all_domain * sizeof(int    *));
  double **pvtx_coord     = malloc(n_part_loc_all_domain * sizeof(double *));

  int shift_part = 0;
  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < mi->n_part[i_domain]; ++i_part) {
      pn_cell       [i_part+shift_part] = mi->parts[i_domain][i_part].n_cell;
      pcell_face_idx[i_part+shift_part] = mi->parts[i_domain][i_part].cell_face_idx;
      pcell_face    [i_part+shift_part] = mi->parts[i_domain][i_part].cell_face;
      pface_edge_idx[i_part+shift_part] = mi->parts[i_domain][i_part].face_edge_idx;
      pface_edge    [i_part+shift_part] = mi->parts[i_domain][i_part].face_edge;
      pface_vtx_idx [i_part+shift_part] = mi->parts[i_domain][i_part].face_vtx_idx;
      pface_vtx     [i_part+shift_part] = mi->parts[i_domain][i_part].face_vtx;
      pedge_vtx     [i_part+shift_part] = mi->parts[i_domain][i_part].edge_vtx;
      pvtx_coord    [i_part+shift_part] = mi->parts[i_domain][i_part].vtx;
    }
    shift_part += mi->n_part[i_domain];
  }

  assert(mi->cell_center == NULL);
  _cell_center_3d(n_part_loc_all_domain,
                  pn_cell,
                  pcell_face_idx,
                  pcell_face,
                  pface_edge_idx,
                  pface_edge,
                  pface_vtx_idx,
                  pface_vtx,
                  pedge_vtx,
                  pvtx_coord,
                  &mi->cell_center);

  free(pn_cell);
  free(pcell_face_idx);
  free(pcell_face    );
  free(pface_edge_idx);
  free(pface_edge    );
  free(pface_vtx_idx );
  free(pface_vtx     );
  free(pedge_vtx     );
  free(pvtx_coord    );


}


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
  PDM_UNUSED(interp_kind);
  PDM_mesh_interpolate_t* mi = malloc(sizeof(PDM_mesh_interpolate_t));

  mi->n_domain    = n_domain;
  mi->n_part      = malloc( n_domain * sizeof(int)); // Make a copy to avoid pb in cython
  mi->n_group     = malloc( n_domain * sizeof(int)); // Make a copy to avoid pb in cython
  for(int i = 0; i < mi->n_domain; ++i) {
    mi->n_part [i] = n_part [i];
    mi->n_group[i] = n_group[i];
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


  mi->group_entity_idx = malloc( PDM_GEOMETRY_KIND_MAX * sizeof(int ***));
  mi->group_entity     = malloc( PDM_GEOMETRY_KIND_MAX * sizeof(int ***));
  mi->group_is_defined = malloc( PDM_GEOMETRY_KIND_MAX * sizeof(int *  ));

  for(int i_kind = 0; i_kind < PDM_GEOMETRY_KIND_MAX; ++i_kind) {
    mi->group_entity_idx[i_kind] = malloc(n_domain * sizeof(int **));
    mi->group_entity    [i_kind] = malloc(n_domain * sizeof(int **));
    mi->group_is_defined[i_kind] = 0;
    for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
      mi->group_entity_idx[i_kind][i_domain] = malloc(n_part[i_domain] * sizeof(int *));
      mi->group_entity    [i_kind][i_domain] = malloc(n_part[i_domain] * sizeof(int *));
      for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
        mi->group_entity_idx[i_kind][i_domain][i_part] = NULL;
        mi->group_entity    [i_kind][i_domain][i_part] = NULL;
      }
    }
  }

  mi->cell_center =  NULL;

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

  _prepare_cell_center(mi);


  /* Deduce graph with all graphe inside same domain and between domain */
  // int ***vtx_part_bound_proc_idx = mi->entity_part_bound_proc_idx[PDM_MESH_ENTITY_VERTEX];
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

    int is_describe_vtx  = 0; //PDM_part_domain_interface_exist_get(mi->pdi, PDM_BOUND_TYPE_VTX );
    int is_describe_edge = PDM_part_domain_interface_exist_get(mi->pdi, PDM_BOUND_TYPE_EDGE);
    int is_describe_face = PDM_part_domain_interface_exist_get(mi->pdi, PDM_BOUND_TYPE_FACE);

    int is_describe_vtx_l  = is_describe_vtx;
    int is_describe_edge_l = is_describe_edge;
    int is_describe_face_l = is_describe_face;
    PDM_MPI_Allreduce(&is_describe_vtx_l , &is_describe_vtx , 1, PDM_MPI_INT, PDM_MPI_MAX, mi->comm);
    PDM_MPI_Allreduce(&is_describe_edge_l, &is_describe_edge, 1, PDM_MPI_INT, PDM_MPI_MAX, mi->comm);
    PDM_MPI_Allreduce(&is_describe_face_l, &is_describe_face, 1, PDM_MPI_INT, PDM_MPI_MAX, mi->comm);


    int **pn_vtx = malloc(mi->n_domain * sizeof(int *));
    for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
      pn_vtx[i_domain] = malloc(mi->n_part[i_domain] * sizeof(int));
      for(int i_part = 0; i_part < mi->n_part[i_domain]; ++i_part) {
        pn_vtx[i_domain][i_part] = mi->parts[i_domain][i_part].n_vtx;
      }
    }

    int          **pn_face        = malloc(mi->n_domain * sizeof(int          *));
    PDM_g_num_t ***pface_ln_to_gn = malloc(mi->n_domain * sizeof(PDM_g_num_t **));
    PDM_g_num_t ***pvtx_ln_to_gn  = malloc(mi->n_domain * sizeof(PDM_g_num_t **));
    int         ***pface_vtx_idx  = malloc(mi->n_domain * sizeof(int         **));
    int         ***pface_vtx      = malloc(mi->n_domain * sizeof(int         **));
    for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
      pn_face       [i_domain] = malloc(mi->n_part[i_domain] * sizeof(int          ));
      pvtx_ln_to_gn [i_domain] = malloc(mi->n_part[i_domain] * sizeof(PDM_g_num_t *));
      pface_ln_to_gn[i_domain] = malloc(mi->n_part[i_domain] * sizeof(PDM_g_num_t *));
      pface_vtx_idx [i_domain] = malloc(mi->n_part[i_domain] * sizeof(int         *));
      pface_vtx     [i_domain] = malloc(mi->n_part[i_domain] * sizeof(int         *));
      for(int i_part = 0; i_part < mi->n_part[i_domain]; ++i_part) {
        pn_face       [i_domain][i_part] = mi->parts[i_domain][i_part].n_face;
        pvtx_ln_to_gn [i_domain][i_part] = mi->parts[i_domain][i_part].vtx_ln_to_gn;
        pface_ln_to_gn[i_domain][i_part] = mi->parts[i_domain][i_part].face_ln_to_gn;
        pface_vtx_idx [i_domain][i_part] = mi->parts[i_domain][i_part].face_vtx_idx;
        pface_vtx     [i_domain][i_part] = mi->parts[i_domain][i_part].face_vtx;
      }
    }

    PDM_part_domain_interface_face2vtx(mi->pdi,
                                       mi->n_part,
                                       pn_face,
                                       pface_ln_to_gn,
                                       pn_vtx,
                                       pvtx_ln_to_gn,
                                       pface_vtx_idx,
                                       pface_vtx);

    printf("is_describe_vtx_l  = %i \n", is_describe_vtx_l );
    printf("is_describe_edge_l = %i \n", is_describe_edge_l);
    printf("is_describe_face_l = %i \n", is_describe_face_l);

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
      free(pn_face       [i_domain]);
      free(pvtx_ln_to_gn [i_domain]);
      free(pface_ln_to_gn[i_domain]);
      free(pface_vtx_idx [i_domain]);
      free(pface_vtx     [i_domain]);
      free(pn_vtx        [i_domain]);
    }
    free(pn_vtx);
    free(pn_face);
    free(pvtx_ln_to_gn );
    free(pface_ln_to_gn);
    free(pface_vtx_idx );
    free(pface_vtx     );
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

      if(0 == 1) {
        PDM_log_trace_array_int(_neighbor_n, n_vtx, "_neighbor_n ::");
      }

      /* Compute index */
      _neighbor_idx[0] = 0;
      for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
        _neighbor_idx[i_entity+1] = _neighbor_idx[i_entity] + _neighbor_n[i_entity];
        _neighbor_n[i_entity] = 0;
      }

      if(0 == 1) {
        PDM_log_trace_array_int(_neighbor_idx, n_vtx, "_neighbor_idx ::");
      }

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
   * Create also vtx_face_bound :
   *   - Vertex can be on 2 boundaries
   */
  shift_part = 0;
  int **vtx_face_bound_idx   = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int * ));
  int **vtx_face_bound_n     = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int * ));
  int **vtx_face_bound       = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int * ));
  int **vtx_face_bound_group = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int * ));
  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {

    int n_group = mi->n_group[i_domain];

    // int *
    for(int i_part = 0; i_part < mi->n_part[i_domain]; ++i_part) {

      int n_vtx = mi->parts[i_domain][i_part].n_vtx;
      vtx_face_bound_idx[i_part+shift_part] = malloc( (n_vtx + 1) * sizeof(int));
      vtx_face_bound_n  [i_part+shift_part] = PDM_array_zeros_int(n_vtx);

      int *_vtx_face_bound_idx = (int *) vtx_face_bound_idx[i_part+shift_part];
      int *_vtx_face_bound_n   = (int *) vtx_face_bound_n  [i_part+shift_part];


      int *face_group_idx = mi->group_entity_idx[PDM_BOUND_TYPE_FACE][i_domain][i_part];
      int *face_group     = mi->group_entity    [PDM_BOUND_TYPE_FACE][i_domain][i_part];

      int* face_vtx_idx = mi->parts[i_domain][i_part].face_vtx_idx;
      int* face_vtx     = mi->parts[i_domain][i_part].face_vtx;

      for(int i_group = 0; i_group < n_group; ++i_group) {
        for(int idx_face = face_group_idx[i_group]; idx_face < face_group_idx[i_group+1]; ++idx_face) {
          int i_face = face_group[idx_face]-1;
          for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
            int i_vtx = face_vtx[idx_vtx];
            _vtx_face_bound_n[i_vtx] += 1;
          }
        }
      }

      _vtx_face_bound_idx[0] = 0;
      for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx ) {
        _vtx_face_bound_idx[i_vtx+1] = _vtx_face_bound_idx[i_vtx] + _vtx_face_bound_n[i_vtx];
        _vtx_face_bound_n[i_vtx] = 0;
      }
      vtx_face_bound      [i_part+shift_part] = malloc(_vtx_face_bound_idx[n_vtx] * sizeof(int));
      vtx_face_bound_group[i_part+shift_part] = malloc(_vtx_face_bound_idx[n_vtx] * sizeof(int));
      int *_vtx_face_bound       =vtx_face_bound      [i_part+shift_part];
      int *_vtx_face_bound_group =vtx_face_bound_group[i_part+shift_part];

      for(int i_group = 0; i_group < n_group; ++i_group) {
        for(int idx_face = face_group_idx[i_group]; idx_face < face_group_idx[i_group+1]; ++idx_face) {
          int i_face = face_group[idx_face]-1;
          for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
            int i_vtx = face_vtx[idx_vtx];
            int idx_write = _vtx_face_bound_idx[i_vtx] + _vtx_face_bound_n[i_vtx]++;
            _vtx_face_bound      [idx_write] = idx_face; // - face_group_idx[i_group]
            _vtx_face_bound_group[idx_write] = i_group;
          }
        }
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

  /* Prepare coordinates to send */
  int    **pvtx_cell_coords_n =  malloc(n_part_loc_all_domain * sizeof(int    *));
  double **pvtx_cell_coords   =  malloc(n_part_loc_all_domain * sizeof(double *));
  shift_part = 0;
  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
    /* First loop to count */
    for(int i_part = 0; i_part < mi->n_part[i_domain]; ++i_part) {
      int n_vtx = pn_vtx[i_part+shift_part];
      int* _neighbor_idx   = neighbor_idx[i_part+shift_part];

      pvtx_cell_coords_n[i_part+shift_part] = PDM_array_zeros_int(n_vtx);
      int *_pvtx_cell_coords_n = pvtx_cell_coords_n[i_part+shift_part];
      int *_pvtx_cell_n        = pvtx_cell_n       [i_part+shift_part];
      int *_pvtx_cell_idx      = pvtx_cell_idx     [i_part+shift_part];
      int *_pvtx_cell          = pvtx_cell         [i_part+shift_part];

      double*  _pcell_center = mi->cell_center[i_part+shift_part];

      /* Count */
      int n_vtx_cell_to_send = 0;
      for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
        int n_neight = _neighbor_idx[i_entity+1] - _neighbor_idx[i_entity];
        if(n_neight > 0 ) {
          n_vtx_cell_to_send += _pvtx_cell_n[i_entity];
          _pvtx_cell_coords_n[i_entity] = _pvtx_cell_n[i_entity];
        }
      }

      pvtx_cell_coords[i_part+shift_part] = malloc( 3 * n_vtx_cell_to_send * sizeof(double));
      double *_pvtx_cell_coords = pvtx_cell_coords[i_part+shift_part];

      n_vtx_cell_to_send = 0;
      for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
        int n_neight = _neighbor_idx[i_entity+1] - _neighbor_idx[i_entity];
        if(n_neight > 0 ) {
          for(int idx_cell = _pvtx_cell_idx[i_entity]; idx_cell < _pvtx_cell_idx[i_entity+1]; ++idx_cell) {
            int  i_cell = PDM_ABS(_pvtx_cell[idx_cell])-1;

            _pvtx_cell_coords[3*n_vtx_cell_to_send  ] = _pcell_center[3*i_cell  ];
            _pvtx_cell_coords[3*n_vtx_cell_to_send+1] = _pcell_center[3*i_cell+1];
            _pvtx_cell_coords[3*n_vtx_cell_to_send+2] = _pcell_center[3*i_cell+2];

          }
          n_vtx_cell_to_send += _pvtx_cell_n[i_entity];
        }
      }
    }
    shift_part += mi->n_part[i_domain];
  }


  PDM_distant_neighbor_free(dn);

  shift_part = 0;
  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < mi->n_part[i_domain]; ++i_part) {
      free(pvtx_cell_coords  [i_part+shift_part]);
      free(pvtx_cell_coords_n[i_part+shift_part]);
    }
    shift_part   += mi->n_part[i_domain];
  }
  free(pvtx_cell_coords  );
  free(pvtx_cell_coords_n);


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

      free(vtx_face_bound_idx  [i_part+shift_part]);
      free(vtx_face_bound_n    [i_part+shift_part]);
      free(vtx_face_bound      [i_part+shift_part]);
      free(vtx_face_bound_group[i_part+shift_part]);

    }
    shift_part   += mi->n_part[i_domain];
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
  free(vtx_face_bound_idx  );
  free(vtx_face_bound_n    );
  free(vtx_face_bound      );
  free(vtx_face_bound_group);

  if(mi->pdi != NULL) {
    free(pdi_neighbor_idx);
    free(pdi_neighbor    );

    free(composed_interface_idx  );
    free(composed_interface      );
    free(composed_ln_to_gn_sorted);
  }

  /* Compute weight */



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
  mi->graph_comm_is_defined     [mesh_entity] = 1;
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

/*
 * A voir avec Julien on passe group par group ?
 */
// void
// PDM_mesh_interpolate_part_group_set
// (
//   PDM_mesh_interpolate_t   *mi,
//   int                       i_domain,
//   int                       i_part,
//   int                       i_group,
//   PDM_mesh_entities_t       entity_kind,
//   int                       n_entity_group_entity,
//   int                      *group_entity
// );

void
PDM_mesh_interpolate_part_group_set
(
  PDM_mesh_interpolate_t   *mi,
  int                       i_domain,
  int                       i_part,
  PDM_bound_type_t          bound_type,
  int                      *entity_bound_idx,
  int                      *entity_bound
)
{
  // assert(i_group < mi->n_group[i_domain]);

  mi->group_entity_idx[bound_type][i_domain][i_part] = entity_bound_idx;
  mi->group_entity    [bound_type][i_domain][i_part] = entity_bound;
  mi->group_is_defined[bound_type] = 1;
}

// Ou calculer in interne ?
// void
// PDM_mesh_interpolate_dual_coord_set
// (
//   PDM_mesh_interpolate_t   *mi,
//   int                       i_domain,
//   int                       i_part,
//   double                   *dual_center_coords
// )
// {

// }

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

  for(int i_kind = 0; i_kind < PDM_GEOMETRY_KIND_MAX; ++i_kind) {
    mi->group_is_defined     [i_kind] = 0;
    for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
      free(mi->group_entity_idx[i_kind][i_domain]);
      free(mi->group_entity    [i_kind][i_domain]);
    }
    free(mi->group_entity_idx[i_kind]);
    free(mi->group_entity    [i_kind]);
  }
  free(mi->group_entity_idx);
  free(mi->group_entity);
  free(mi->group_is_defined     );

  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
    free(mi->parts      [i_domain]);
    free(mi->cell_center[i_domain]);
  }
  free(mi->parts);
  free(mi->n_part_idx);
  free(mi->n_part_g_idx);
  free(mi->n_part);
  free(mi->n_group);
  free(mi->cell_center);


  free(mi);
}


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
