
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
#include "pdm_vtk.h"
#include "pdm_order.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_part_domain_interface_priv.h"
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

// static
// void
// triplet_to_quadruplet
// (
//   int   size,
//   int  *triplet,
//   int  *array,
//   int **quadruplet
// )
// {
//   int *_quadruplet = malloc(4 * size * sizeof(int));
//   for(int i = 0; i < size; ++i) {
//     _quadruplet[4*i  ] = triplet[3*i  ];
//     _quadruplet[4*i+1] = triplet[3*i+1];
//     _quadruplet[4*i+2] = triplet[3*i+2];
//     _quadruplet[4*i+3] = array  [i];
//   }


//   *quadruplet = _quadruplet;
// }

// static inline
// int
// _is_same_quadruplet
// (
// int iproc1, int ipart1, int ielt1, int iinterf1,
// int iproc2, int ipart2, int ielt2, int iinterf2
// )
// {
//   if(iproc1 == iproc2){
//     if(ipart1 == ipart2){
//       if(ielt1 == ielt2){
//         if(iinterf1 == iinterf2){
//           return 1;
//         }
//       }
//     }
//   }
//   return 0;
// }

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

        for(int idx_face = _pcell_face_idx[i_cell]; idx_face < _pcell_face_idx[i_cell+1]; ++idx_face) {
          int i_face = PDM_ABS(_pcell_face[idx_face])-1;

          double inv2 = 1./((double)  _pface_vtx_idx[i_face+1] - _pface_vtx_idx[i_face]);

          double fcx = 0;
          double fcy = 0;
          double fcz = 0;
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
  int     *pn_cell        = malloc(mi->n_part_loc_all_domain * sizeof(int    *));
  int    **pcell_face_idx = malloc(mi->n_part_loc_all_domain * sizeof(int    *));
  int    **pcell_face     = malloc(mi->n_part_loc_all_domain * sizeof(int    *));
  int    **pface_edge_idx = malloc(mi->n_part_loc_all_domain * sizeof(int    *));
  int    **pface_edge     = malloc(mi->n_part_loc_all_domain * sizeof(int    *));
  int    **pface_vtx_idx  = malloc(mi->n_part_loc_all_domain * sizeof(int    *));
  int    **pface_vtx      = malloc(mi->n_part_loc_all_domain * sizeof(int    *));
  int    **pedge_vtx      = malloc(mi->n_part_loc_all_domain * sizeof(int    *));
  double **pvtx_coord     = malloc(mi->n_part_loc_all_domain * sizeof(double *));

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
  _cell_center_3d(mi->n_part_loc_all_domain,
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

  if(1 == 1) {

    int i_rank;
    PDM_MPI_Comm_rank(mi->comm, &i_rank);

    shift_part = 0;
    for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < mi->n_part[i_domain]; ++i_part) {
        char filename[999];
        sprintf(filename, "vtx_coords_%i_%i.vtk", i_rank, i_part+shift_part);
        PDM_vtk_write_point_cloud(filename,
                                  mi->parts[i_domain][i_part].n_vtx,
                                  mi->parts[i_domain][i_part].vtx,
                                  NULL,
                                  NULL);

        sprintf(filename, "cell_center_%i_%i.vtk", i_rank, i_part+shift_part);
        PDM_vtk_write_point_cloud(filename,
                                  mi->parts[i_domain][i_part].n_cell,
                                  mi->cell_center[i_part+shift_part],
                                  NULL,
                                  NULL);
      }
      shift_part += mi->n_part[i_domain];
    }
  }

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

static
void
_prepare_vtx_cell
(
 PDM_mesh_interpolate_t*  mi
)
{

  /*
   * Create vtx_cell
   */
  int shift_part = 0;
  mi->pvtx_cell_n   = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int *));
  mi->pvtx_cell_idx = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int *));
  mi->pvtx_cell     = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int *));
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
                                 &mi->pvtx_cell_idx[i_part+shift_part],
                                 &mi->pvtx_cell    [i_part+shift_part]);

      free(cell_vtx_idx);
      free(cell_vtx);

      mi->pvtx_cell_n[i_part+shift_part] = malloc(mi->parts[i_domain][i_part].n_vtx * sizeof(int));
      for(int i_vtx = 0; i_vtx < mi->parts[i_domain][i_part].n_vtx; ++i_vtx) {
        mi->pvtx_cell_n[i_part+shift_part][i_vtx] = mi->pvtx_cell_idx[i_part+shift_part][i_vtx+1] - mi->pvtx_cell_idx[i_part+shift_part][i_vtx];
      }
    }
    shift_part += mi->n_part[i_domain];
  }
}





static
void
_warm_up_distant_neighbor
(
 PDM_mesh_interpolate_t*  mi
)
{
  /* Deduce graph with all graphe inside same domain and between domain */
  // int ***vtx_part_bound_proc_idx = mi->entity_part_bound_proc_idx[PDM_MESH_ENTITY_VERTEX];
  int ***vtx_part_bound_part_idx = mi->entity_part_bound_part_idx[PDM_MESH_ENTITY_VERTEX];
  int ***vtx_part_bound          = mi->entity_part_bound         [PDM_MESH_ENTITY_VERTEX];


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

  mi->neighbor_idx       = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int *));
  mi->neighbor_desc      = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int *));
  mi->neighbor_interface = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int *));
  mi->pn_vtx             = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int  ));
  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {

    int n_part_total = mi->n_part_g_idx[i_domain+1] - mi->n_part_g_idx[i_domain];

    /* First loop to count */
    for(int i_part = 0; i_part < mi->n_part[i_domain]; ++i_part) {

      int n_vtx = mi->parts[i_domain][i_part].n_vtx;
      mi->pn_vtx[i_part+shift_part] = n_vtx;

      mi->neighbor_idx[i_part+shift_part] = malloc((n_vtx+1) * sizeof(int));
      int* _neighbor_idx   = mi->neighbor_idx[i_part+shift_part];
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

      mi->neighbor_desc     [i_part+shift_part] = (int *) malloc( 3 * _neighbor_idx[n_vtx] * sizeof(int) );
      mi->neighbor_interface[i_part+shift_part] = (int *) malloc(     _neighbor_idx[n_vtx] * sizeof(int) );
      int* _neighbor_desc      = mi->neighbor_desc     [i_part+shift_part];
      int* _neighbor_interface = mi->neighbor_interface[i_part+shift_part];

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

  shift_part = 0;
  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < mi->n_part[i_domain]; ++i_part) {

      if(mi->pdi != NULL) {
        free(pdi_neighbor_idx[i_part+shift_part]);
        free(pdi_neighbor    [i_part+shift_part]);
      }


    }
    shift_part   += mi->n_part[i_domain];
  }


  if(mi->pdi != NULL) {
    free(pdi_neighbor_idx);
    free(pdi_neighbor    );

    free(composed_interface_idx  );
    free(composed_interface      );
    free(composed_ln_to_gn_sorted);
  }

  /*
   * Pruned the connectivity  - Useless ause  the  interessed pruned was to sort cell
   */
  // int i_rank;
  // PDM_MPI_Comm_rank(mi->comm, &i_rank);
  // for(int i_part = 0; i_part < mi->n_part_loc_all_domain; ++i_part) {


  //   int *neighbor_idx       = mi->neighbor_idx      [i_part];
  //   int *neighbor_desc      = mi->neighbor_desc     [i_part];
  //   int *neighbor_interface = mi->neighbor_interface[i_part];

  //   int n_vtx = mi->pn_vtx[i_part];
  //   int s_tot = neighbor_idx[n_vtx];

  //   int *_quad_cell_cell_extended = NULL;
  //   triplet_to_quadruplet(s_tot,
  //                         mi->neighbor_desc     [i_part],
  //                         mi->neighbor_interface[i_part],
  //                         &_quad_cell_cell_extended);


  //   int* order     = (int * ) malloc( s_tot * sizeof(int));
  //   PDM_order_lnum_s(_quad_cell_cell_extended, 4, order, s_tot);

  //   int* pruned_neighbor_idx  = malloc((n_vtx+1)   * sizeof(int));
  //   int* pruned_neighbor      = malloc(3 * s_tot   * sizeof(int));
  //   int* pruned_neighbor_itrf = malloc(    s_tot   * sizeof(int));
  //   int* index_border         = malloc(    s_tot   * sizeof(int));

  //   int idx_unique = 0;
  //   int last_proc  = -1;
  //   int last_part  = -1;
  //   int last_elmt  = -1;
  //   int last_inte  = -40000;
  //   pruned_neighbor_idx[0] = 0;
  //   pruned_neighbor_idx[1] = 0;
  //   for(int i = 0; i < s_tot; ++i) {
  //     int old_order = order[i];

  //     int curr_proc = _quad_cell_cell_extended[4*old_order  ];
  //     int curr_part = _quad_cell_cell_extended[4*old_order+1];
  //     int curr_cell = _quad_cell_cell_extended[4*old_order+2];
  //     int curr_inte = _quad_cell_cell_extended[4*old_order+3];
  //     int is_same = _is_same_quadruplet(last_proc, last_part, last_elmt, last_inte,
  //                                       curr_proc, curr_part, curr_cell, curr_inte);

  //     int is_local = (curr_proc == i_rank) && (curr_part == i_part+shift_part) && (curr_inte == -40000);
  //     if(is_same == 0 && !is_local){ // N'est pas le meme
  //       // idx_unique++;
  //       last_proc = curr_proc;
  //       last_part = curr_part;
  //       last_elmt = curr_cell;
  //       last_inte = curr_inte;

  //       pruned_neighbor     [3*idx_unique  ] = curr_proc;
  //       pruned_neighbor     [3*idx_unique+1] = curr_part;
  //       pruned_neighbor     [3*idx_unique+2] = curr_cell;
  //       pruned_neighbor_itrf[idx_unique    ] = curr_inte;
  //       index_border        [old_order     ] = idx_unique;
  //       idx_unique++;

  //       /* Increment the new counter */
  //       pruned_neighbor_idx[1]++;
  //     } else {

  //       last_proc = curr_proc;
  //       last_part = curr_part;
  //       last_elmt = curr_cell;
  //       last_inte = curr_inte;

  //       index_border[old_order] = idx_unique-1;

  //     }
  //   }

  //   PDM_log_trace_array_int(index_border, s_tot, "index_border ::");

  //   free(pruned_neighbor_idx );
  //   free(pruned_neighbor     );
  //   free(pruned_neighbor_itrf);
  //   free(index_border        );

  //   // free(mi->neighbor_idx      [i_part]);
  //   // free(mi->neighbor_desc     [i_part]);
  //   // free(mi->neighbor_interface[i_part]);

  //   free(order);
  // }


}

static
void
_create_bnd_graph
(
 PDM_mesh_interpolate_t* mi
)
{
  /*
   * Create also vtx_face_bound :
   *   - Vertex can be on 2 boundaries
   */
  int shift_part = 0;
  mi->vtx_face_bound_idx    = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int * ));
  mi->vtx_face_bound_n      = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int * ));
  mi->vtx_face_bound        = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int * ));
  mi->vtx_face_bound_group  = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int * ));
  mi->vtx_face_bound_coords = malloc(mi->n_part_g_idx[mi->n_domain] * sizeof(int * ));

  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {

    int n_group = mi->n_group[i_domain];

    for(int i_part = 0; i_part < mi->n_part[i_domain]; ++i_part) {

      int n_vtx = mi->parts[i_domain][i_part].n_vtx;
      mi->vtx_face_bound_idx[i_part+shift_part] = malloc( (n_vtx + 1) * sizeof(int));
      mi->vtx_face_bound_n  [i_part+shift_part] = PDM_array_zeros_int(n_vtx);

      int *_vtx_face_bound_idx = (int *) mi->vtx_face_bound_idx[i_part+shift_part];
      int *_vtx_face_bound_n   = (int *) mi->vtx_face_bound_n  [i_part+shift_part];

      int  *n_face_group = mi->n_group_entity[PDM_BOUND_TYPE_FACE][i_domain][i_part];
      int **face_group   = mi->group_entity  [PDM_BOUND_TYPE_FACE][i_domain][i_part];

      int* face_vtx_idx = mi->parts[i_domain][i_part].face_vtx_idx;
      int* face_vtx     = mi->parts[i_domain][i_part].face_vtx;

      for(int i_group = 0; i_group < n_group; ++i_group) {
        for(int idx_face = 0; idx_face < n_face_group[i_group]; ++idx_face) {
          int i_face = face_group[i_group][idx_face]-1;
          for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
            int i_vtx = face_vtx[idx_vtx]-1;
            _vtx_face_bound_n[i_vtx] += 1;
          }
        }
      }

      _vtx_face_bound_idx[0] = 0;
      for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx ) {
        _vtx_face_bound_idx[i_vtx+1] = _vtx_face_bound_idx[i_vtx] + _vtx_face_bound_n[i_vtx];
        _vtx_face_bound_n[i_vtx] = 0;
      }
      mi->vtx_face_bound       [i_part+shift_part] = malloc(    _vtx_face_bound_idx[n_vtx] * sizeof(int   ));
      mi->vtx_face_bound_group [i_part+shift_part] = malloc(    _vtx_face_bound_idx[n_vtx] * sizeof(int   ));
      mi->vtx_face_bound_coords[i_part+shift_part] = malloc(3 * _vtx_face_bound_idx[n_vtx] * sizeof(double));
      int    *_vtx_face_bound        = mi->vtx_face_bound       [i_part+shift_part];
      int    *_vtx_face_bound_group  = mi->vtx_face_bound_group [i_part+shift_part];
      double *_vtx_face_bound_coords = mi->vtx_face_bound_coords[i_part+shift_part];
      double *_pvtx_coord            = mi->parts[i_domain][i_part].vtx;

      for(int i_group = 0; i_group < n_group; ++i_group) {
        for(int idx_face = 0; idx_face < n_face_group[i_group]; ++idx_face) {
          int i_face = face_group[i_group][idx_face]-1;

          double center_face[3] = {0., 0., 0.};
          double pond = 1./((double)  face_vtx_idx[i_face+1] - face_vtx_idx[i_face]);
          for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
            int i_vtx = face_vtx[idx_vtx]-1;
            for(int k = 0; k < 3; ++k) {
              center_face[k] += _pvtx_coord[3*i_vtx+k];
            }
          }
          for(int k = 0; k < 3; ++k) {
            center_face[k] = center_face[k] * pond;
          }

          for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
            int i_vtx = face_vtx[idx_vtx]-1;
            int idx_write = _vtx_face_bound_idx[i_vtx] + _vtx_face_bound_n[i_vtx]++;
            _vtx_face_bound      [idx_write] = idx_face;
            _vtx_face_bound_group[idx_write] = i_group;
            for(int k = 0; k < 3; ++k) {
              _vtx_face_bound_coords[3*idx_write+k] = center_face[k];
            }
          }
        }
      }
    }
    shift_part   += mi->n_part              [i_domain];
  }
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

  /* Si multidomain on fait un shift et tt roule */
  mi->n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
    mi->n_part_loc_all_domain += mi->n_part[i_domain];
  }

  mi->pn_vtx              = NULL;
  mi->pvtx_cell_n         = NULL;
  mi->pvtx_cell_idx       = NULL;
  mi->pvtx_cell           = NULL;

  mi->neighbor_idx       = NULL;
  mi->neighbor_desc      = NULL;
  mi->neighbor_interface = NULL;

  mi->vtx_face_bound_idx    = NULL;
  mi->vtx_face_bound_n      = NULL;
  mi->vtx_face_bound        = NULL;
  mi->vtx_face_bound_group  = NULL;
  mi->vtx_face_bound_coords = NULL;

  mi->dn = NULL;

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


  mi->n_group_entity   = malloc( PDM_GEOMETRY_KIND_MAX * sizeof(int  ***));
  mi->group_entity     = malloc( PDM_GEOMETRY_KIND_MAX * sizeof(int ****));
  mi->group_is_defined = malloc( PDM_GEOMETRY_KIND_MAX * sizeof(int *   ));

  for(int i_kind = 0; i_kind < PDM_GEOMETRY_KIND_MAX; ++i_kind) {
    mi->n_group_entity  [i_kind] = malloc(n_domain * sizeof(int  **));
    mi->group_entity    [i_kind] = malloc(n_domain * sizeof(int ***));
    mi->group_is_defined[i_kind] = 0;
    for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
      mi->n_group_entity  [i_kind][i_domain] = malloc(n_part[i_domain] * sizeof(int  *));
      mi->group_entity    [i_kind][i_domain] = malloc(n_part[i_domain] * sizeof(int **));
      for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
        mi->n_group_entity  [i_kind][i_domain][i_part] = malloc(n_group[i_domain] * sizeof(int  ));
        mi->group_entity    [i_kind][i_domain][i_part] = malloc(n_group[i_domain] * sizeof(int *));
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

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(mi->comm, &i_rank);
  PDM_MPI_Comm_size(mi->comm, &n_rank);

  /*
   * Compute graph comm from gnum with PDM_part_generate_entity_graph_comm is not provided
   */
  if(mi->graph_comm_is_defined[PDM_MESH_ENTITY_VERTEX] == 0) {
    abort();
    // TODO : compute graph comm from other graph comm OR recumpute from gnum : PDM_part_generate_entity_graph_comm
  }

  _prepare_cell_center     (mi);
  _prepare_vtx_cell        (mi);
  _create_bnd_graph        (mi);
  _warm_up_distant_neighbor(mi);

  int **pvtx_cell_n   = mi->pvtx_cell_n;
  int **pvtx_cell_idx = mi->pvtx_cell_idx;
  int **pvtx_cell     = mi->pvtx_cell;

  /*
   * Create distant_neighbor
   */
  assert(mi->dn == NULL);
  mi->dn = PDM_distant_neighbor_create(mi->comm,
                                       mi->n_part_loc_all_domain,
                                       mi->pn_vtx,
                                       mi->neighbor_idx,
                                       mi->neighbor_desc);

  /* Prepare coordinates to send */
  int    **pvtx_cell_coords_n =  malloc(mi->n_part_loc_all_domain * sizeof(int    *));
  double **pvtx_cell_coords   =  malloc(mi->n_part_loc_all_domain * sizeof(double *));
  int shift_part = 0;
  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
    /* First loop to count */
    for(int i_part = 0; i_part < mi->n_part[i_domain]; ++i_part) {
      int n_vtx          = mi->pn_vtx      [i_part+shift_part];
      int* _neighbor_idx = mi->neighbor_idx[i_part+shift_part];

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
            n_vtx_cell_to_send++;
          }
        }
      }
    }
    shift_part += mi->n_part[i_domain];
  }

  int    **pvtx_cell_coords_opp_n = NULL;
  double **pvtx_cell_coords_opp   = NULL;
  PDM_distant_neighbor_exch(mi->dn,
                            3 * sizeof(double),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            pvtx_cell_coords_n,
                  (void **) pvtx_cell_coords,
                           &pvtx_cell_coords_opp_n,
                 (void ***)&pvtx_cell_coords_opp);

  /* Same but for face_group */
  int    **pvtx_face_bound_coords_opp_n = NULL;
  double **pvtx_face_bound_coords_opp   = NULL;
  PDM_distant_neighbor_exch(mi->dn,
                            3 * sizeof(double),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            mi->vtx_face_bound_n,
                  (void **) mi->vtx_face_bound_coords,
                           &pvtx_face_bound_coords_opp_n,
                 (void ***)&pvtx_face_bound_coords_opp);

  /*
   * Count receive
   */
  for(int i_part = 0; i_part < mi->n_part_loc_all_domain; ++i_part){

    int nrecv = 0;
    int n_vtx = mi->pn_vtx[i_part];
    int* _neighbor_idx       = mi->neighbor_idx      [i_part];
    int* _neighbor_interface = mi->neighbor_interface[i_part];

    for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
      for(int idx_entity = _neighbor_idx[i_entity]; idx_entity < _neighbor_idx[i_entity+1]; ++idx_entity) {

        if(_neighbor_interface[idx_entity] != -40000) {
          int  i_interface = PDM_ABS(_neighbor_interface[idx_entity])-1;
          for(int idx_recv = 0; idx_recv < pvtx_cell_coords_opp_n[i_part][idx_entity]; ++idx_recv){
            for(int k = 0; k < 3; ++k) {
              pvtx_cell_coords_opp[i_part][3*(nrecv+idx_recv)+k] += PDM_SIGN(_neighbor_interface[idx_entity]) * mi->pdi->translation_vect[i_interface][k];
            }
          }
        }
        nrecv += pvtx_cell_coords_opp_n[i_part][idx_entity];
      }
    }

    if(1 == 1) {
      printf("nrecv = %i \n",  nrecv);
      char filename[999];
      sprintf(filename, "opp_coords_%i_%i.vtk", i_rank, i_part);
      PDM_vtk_write_point_cloud(filename,
                                nrecv,
                                pvtx_cell_coords_opp[i_part],
                                NULL,
                                NULL);
    }

    /* Bnd */
    int nrecv_bnd = 0;

    for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
      for(int idx_entity = _neighbor_idx[i_entity]; idx_entity < _neighbor_idx[i_entity+1]; ++idx_entity) {

        if(_neighbor_interface[idx_entity] != -40000) {
          int  i_interface = PDM_ABS(_neighbor_interface[idx_entity])-1;
          for(int idx_recv = 0; idx_recv < pvtx_face_bound_coords_opp_n[i_part][idx_entity]; ++idx_recv){
            for(int k = 0; k < 3; ++k) {
              pvtx_face_bound_coords_opp[i_part][3*(nrecv_bnd+idx_recv)+k] += PDM_SIGN(_neighbor_interface[idx_entity]) * mi->pdi->translation_vect[i_interface][k];
            }
          }
        }
        nrecv_bnd += pvtx_face_bound_coords_opp_n[i_part][idx_entity];
      }
    }

    if(1 == 1) {
      printf("nrecv_bnd = %i \n",  nrecv_bnd);
      char filename[999];
      sprintf(filename, "opp_coords_bnd_%i_%i.vtk", i_rank, i_part);
      PDM_vtk_write_point_cloud(filename,
                                nrecv_bnd,
                                pvtx_face_bound_coords_opp[i_part],
                                NULL,
                                NULL);
    }


  }





  /* Free */
  for(int i_part = 0; i_part < mi->n_part_loc_all_domain; ++i_part){
    free(pvtx_cell_coords_opp_n[i_part]);
    free(pvtx_cell_coords_opp  [i_part]);

    free(pvtx_face_bound_coords_opp_n[i_part]);
    free(pvtx_face_bound_coords_opp  [i_part]);

    free(pvtx_cell_coords  [i_part]);
    free(pvtx_cell_coords_n[i_part]);
  }
  free(pvtx_cell_coords_opp_n);
  free(pvtx_cell_coords_opp  );
  free(pvtx_cell_coords  );
  free(pvtx_cell_coords_n);
  free(pvtx_face_bound_coords_opp_n  );
  free(pvtx_face_bound_coords_opp);

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


void
PDM_mesh_interpolate_part_group_set
(
  PDM_mesh_interpolate_t   *mi,
  int                       i_domain,
  int                       i_part,
  int                       i_group,
  PDM_bound_type_t          bound_type,
  int                       n_group_entity,
  int                      *group_entity
)
{
  assert(i_group < mi->n_group[i_domain]);

  mi->n_group_entity  [bound_type][i_domain][i_part][i_group] = n_group_entity;
  mi->group_entity    [bound_type][i_domain][i_part][i_group] = group_entity;
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

  PDM_distant_neighbor_free(mi->dn);

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
      for(int i_part = 0; i_part < mi->n_part[i_domain]; ++i_part) {
        free(mi->n_group_entity[i_kind][i_domain][i_part]);
        free(mi->group_entity  [i_kind][i_domain][i_part]);
      }
      free(mi->n_group_entity[i_kind][i_domain]);
      free(mi->group_entity  [i_kind][i_domain]);
    }
    free(mi->n_group_entity[i_kind]);
    free(mi->group_entity  [i_kind]);
  }
  free(mi->n_group_entity);
  free(mi->group_entity);
  free(mi->group_is_defined);


  int shift_part = 0;
  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < mi->n_part[i_domain]; ++i_part) {
      free(mi->neighbor_idx      [i_part+shift_part]);
      free(mi->neighbor_desc     [i_part+shift_part]);
      free(mi->neighbor_interface[i_part+shift_part]);

      free(mi->vtx_face_bound_idx   [i_part+shift_part]);
      free(mi->vtx_face_bound_n     [i_part+shift_part]);
      free(mi->vtx_face_bound       [i_part+shift_part]);
      free(mi->vtx_face_bound_group [i_part+shift_part]);
      free(mi->vtx_face_bound_coords[i_part+shift_part]);

      free(mi->pvtx_cell_n  [i_part+shift_part]);
      free(mi->pvtx_cell_idx[i_part+shift_part]);
      free(mi->pvtx_cell    [i_part+shift_part]);
    }
    shift_part += mi->n_part[i_domain];
  }
  free(mi->neighbor_idx      );
  free(mi->neighbor_desc     );
  free(mi->neighbor_interface);
  free(mi->pvtx_cell_n);
  free(mi->pvtx_cell_idx);
  free(mi->pvtx_cell);

  free(mi->vtx_face_bound_idx   );
  free(mi->vtx_face_bound_n     );
  free(mi->vtx_face_bound       );
  free(mi->vtx_face_bound_group );
  free(mi->vtx_face_bound_coords);
  free(mi->pn_vtx);



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
