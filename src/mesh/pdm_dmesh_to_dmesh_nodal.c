
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
#include "pdm_error.h"
#include "pdm_binary_search.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_elmts_utils.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_logging.h"
#include "pdm_dmesh.h"
#include "pdm_gnum.h"
#include "pdm_distrib.h"
#include "pdm_unique.h"
#include "pdm_quick_sort.h"
#include "pdm_para_graph_dual.h"
#include "pdm_array.h"
#include "pdm_order.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_part_mesh_nodal_to_pmesh.h"
#include "pdm_dmesh_to_dmesh_nodal.h"
#include "pdm_dmesh_to_dmesh_nodal_priv.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_sort.h"
#include "pdm_vtk.h"

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

static
PDM_dmesh_nodal_t*
_dmesh_to_dmesh_nodal
(
 PDM_MPI_Comm    comm,
 PDM_g_num_t    *distrib_cell,
 PDM_g_num_t    *distrib_face,
 PDM_g_num_t    *distrib_edge,
 PDM_g_num_t    *distrib_vtx,
 PDM_g_num_t    *dcell_face,
 int            *dcell_face_idx,
 PDM_g_num_t    *dface_edge,
 int            *dface_edge_idx,
 PDM_g_num_t    *dface_vtx,
 int            *dface_vtx_idx,
 PDM_g_num_t    *dedge_vtx,
 int            *n_bound,
 int           **dbound_idx,
 PDM_g_num_t   **dbound,
 int            *n_blk_gnum,
 PDM_g_num_t   **blk_entity_gnum,
 PDM_g_num_t   **blk_elmt_gnum
)
{

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int is_3d = 0;
  int is_2d = 0;
  int mesh_dimension = 3;
  if(distrib_cell != NULL) {
    is_3d = 1;
    mesh_dimension = 3;
  } else {
    assert(distrib_face != NULL);
    is_2d = 1;
    mesh_dimension = 2;
  }
  printf("_dmesh_to_dmesh_nodal --> is_2d = %i | is_3d = %i \n", is_2d, is_3d);

  PDM_g_num_t  n_vtx  = 0;
  PDM_g_num_t  n_cell = 0;
  PDM_g_num_t  n_face = 0;
  PDM_g_num_t  n_edge = 0;

  if(distrib_cell != NULL) {
    n_cell = distrib_cell[n_rank];
  }

  if(distrib_face != NULL) {
    n_face = distrib_face[n_rank];
  }

  if(distrib_edge != NULL) {
    n_edge = distrib_edge[n_rank];
  }

  if(distrib_vtx != NULL) {
    n_vtx = distrib_vtx[n_rank];
  }

  PDM_dmesh_nodal_t* dmn = PDM_DMesh_nodal_create(comm,
                                                  mesh_dimension,
                                                  n_vtx,
                                                  n_cell,
                                                  n_face,
                                                  n_edge);

  if(is_3d) {

  } else {

    /* Create implicit partitionning */
    int dn_face = distrib_face[i_rank+1] - distrib_face[i_rank];
    PDM_g_num_t* pface_ln_to_gn = malloc(dn_face * sizeof(PDM_g_num_t));
    for(int i_face = 0; i_face < dn_face; ++i_face) {
      pface_ln_to_gn[i_face] = distrib_face[i_rank] + i_face + 1;
    }

    int dn_edge = distrib_edge[i_rank+1] - distrib_edge[i_rank];
    int* dedge_vtx_idx = malloc( (dn_edge+1) * sizeof(int));
    for(int i_edge = 0; i_edge < dn_edge+1; ++i_edge) {
      dedge_vtx_idx[i_edge] = 2 * i_edge;
    }

    int pn_vtx = 0;
    int         *pface_vtx_idx = NULL;
    int         *pface_vtx     = NULL;
    PDM_g_num_t *pvtx_ln_to_gn = NULL;

    /* Translate */
    if(dface_vtx == NULL){
      int pn_edge = 0;
      int         *pface_edge_idx = NULL;
      int         *pface_edge     = NULL;
      PDM_g_num_t *pedge_ln_to_gn = NULL;

      PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                               distrib_face,
                                                               dface_edge_idx,
                                                               dface_edge,
                                                               dn_face,
                                                               pface_ln_to_gn,
                                                               &pn_edge,
                                                               &pedge_ln_to_gn,
                                                               &pface_edge_idx,
                                                               &pface_edge);

      pn_vtx = 0;
      int         *pedge_vtx_idx = NULL;
      int         *pedge_vtx     = NULL;

      PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                               distrib_edge,
                                                               dedge_vtx_idx,
                                                               dedge_vtx,
                                                               pn_edge,
                                                               pedge_ln_to_gn,
                                                               &pn_vtx,
                                                               &pvtx_ln_to_gn,
                                                               &pedge_vtx_idx,
                                                               &pedge_vtx);

      PDM_compute_face_vtx_from_face_and_edge(dn_face,
                                              pface_edge_idx,
                                              pface_edge,
                                              pedge_vtx,
                                              &pface_vtx);

      // PDM_log_trace_connectivity_int(pedge_vtx_idx, pedge_vtx, pn_edge, "pedge_vtx ::");

      pface_vtx_idx = pface_edge_idx;

      free(pface_edge    );
      free(pedge_ln_to_gn);
      free(pedge_vtx_idx);
      free(pedge_vtx    );
    } else {
      PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                               distrib_face,
                                                               dface_vtx_idx,
                                                               dface_vtx,
                                                               dn_face,
                                                               pface_ln_to_gn,
                                                               &pn_vtx,
                                                               &pvtx_ln_to_gn,
                                                               &pface_vtx_idx,
                                                               &pface_vtx);
    }

    /*
     * Reconstruction of section for surfacic
     */
    PDM_g_num_t          *section_n    = malloc((dn_face+1) * sizeof(PDM_g_num_t         )); // Suralloc
    PDM_Mesh_nodal_elt_t *section_kind = malloc( dn_face    * sizeof(PDM_Mesh_nodal_elt_t)); // Suralloc

    int n_section_tot    = 0;
    int n_section_tri    = 0;
    int n_section_quad   = 0;
    int n_section_poly2d = 0;
    int ln_vtx_old = -1;
    for(int i_face = 0; i_face < dn_face; ++i_face) {
      int ln_vtx = pface_vtx_idx[i_face+1] - pface_vtx_idx[i_face];
      if(ln_vtx_old == ln_vtx){
        section_n[n_section_tot-1]++;
        continue;
      }
      if(ln_vtx == 3) {
        section_kind[n_section_tot] = PDM_MESH_NODAL_TRIA3;
        n_section_tri++;
      } else if(ln_vtx == 4){
        section_kind[n_section_tot] = PDM_MESH_NODAL_QUAD4;
        n_section_quad++;
      } else {
        section_kind[n_section_tot] = PDM_MESH_NODAL_POLY_2D;
        n_section_poly2d++;
      }
      section_n[n_section_tot] = 1;
      n_section_tot++;
      ln_vtx_old = ln_vtx;
    }
    section_n    = realloc(section_n   , (n_section_tot+1) * sizeof(PDM_g_num_t         ));
    section_kind = realloc(section_kind,  n_section_tot    * sizeof(PDM_Mesh_nodal_elt_t));

    /*
     * Rebuild sections
     */
    int *section_kind_n = malloc(n_rank * sizeof(int));
    PDM_MPI_Allgather(&n_section_tot, 1, PDM_MPI_INT, section_kind_n, 1, PDM_MPI_INT, comm);

    int *section_kind_idx = malloc((n_rank+1) * sizeof(int));
    section_kind_idx[0] = 0;
    for(int i = 0; i < n_rank; ++i) {
      section_kind_idx[i+1] = section_kind_idx[i] + section_kind_n[i];
    }


    int         *g_section_kind = (int         *) malloc(section_kind_idx[n_rank] * sizeof(int        ));
    PDM_g_num_t *g_section_n      = (PDM_g_num_t *) malloc(section_kind_idx[n_rank] * sizeof(PDM_g_num_t));

    PDM_MPI_Allgatherv(section_kind, n_section_tot, PDM_MPI_INT,
                       g_section_kind, section_kind_n, section_kind_idx, PDM_MPI_INT, comm);
    PDM_MPI_Allgatherv(section_n, n_section_tot, PDM__PDM_MPI_G_NUM,
                       g_section_n, section_kind_n, section_kind_idx, PDM__PDM_MPI_G_NUM, comm);

    if(0 == 0) {
      PDM_log_trace_array_long(section_n, n_section_tot, "section_n ::");
      PDM_log_trace_array_long(distrib_face, n_rank+1, "distrib_face ::");
      PDM_log_trace_connectivity_int(section_kind_idx, g_section_kind, n_rank, "g_section_kind : ");
      PDM_log_trace_connectivity_int(section_kind_idx, g_section_n   , n_rank, "g_section_n      : ");
    }

    /*
     * Post-treament
     */
    int n_section_all_rank = section_kind_idx[n_rank];
    int         *post_section_kind = (int         *) malloc(section_kind_idx[n_rank] * sizeof(int        ));
    PDM_g_num_t *post_section_n    = (PDM_g_num_t *) malloc(section_kind_idx[n_rank] * sizeof(PDM_g_num_t));

    int n_section_post = 0;
    PDM_Mesh_nodal_elt_t first = PDM_MESH_NODAL_N_ELEMENT_TYPES;
    for(int i = 0; i < n_section_all_rank; ++i) {
      if(first != (PDM_Mesh_nodal_elt_t) g_section_kind[i]) {
        first = g_section_kind[i];
        post_section_kind[n_section_post] = first;
        post_section_n   [n_section_post] = g_section_n[i];
        n_section_post++;
      } else {
        post_section_n[n_section_post-1] += g_section_n[i];
      }
    }
    post_section_kind = realloc(post_section_kind, n_section_post * sizeof(int        ));
    post_section_n    = realloc(post_section_n   , n_section_post * sizeof(PDM_g_num_t));

    free(g_section_kind);
    free(g_section_n);
    free(section_kind_idx);
    free(section_kind_n);
    free(section_n);
    free(section_kind);

    if(0 == 0) {
      PDM_log_trace_array_int (post_section_kind, n_section_post, "post_section_kind ::");
      PDM_log_trace_array_long(post_section_n   , n_section_post, "post_section_n    ::");
    }

    /*
     * Compute global shift
     */
    PDM_g_num_t* post_section_idx = malloc((n_section_post+1) * sizeof(PDM_g_num_t));
    post_section_idx[0] = 0;
    for(int i = 0; i < n_section_post; ++i) {
      post_section_idx[i+1] = post_section_idx[i] + post_section_n[i];
    }

    /*
     * Reanrange in global section the local one
     */
    int* local_post_section_n   = malloc((n_section_post  ) * sizeof(int));
    int* local_post_section_idx = malloc((n_section_post+1) * sizeof(int));
    for(int i = 0; i < n_section_post; ++i) {
      local_post_section_n[i] = 0;
    }

    for(int i_face = 0; i_face < dn_face; ++i_face) {
      PDM_g_num_t gnum = distrib_face[i_rank] + i_face + 1;
      int t_section = PDM_binary_search_gap_long(gnum-1, post_section_idx, n_section_post+1);
      local_post_section_n[t_section]++;
    }

    local_post_section_idx[0] = 0;
    for(int i = 0; i < n_section_post; ++i) {
      local_post_section_idx[i+1] = local_post_section_idx[i] + local_post_section_n[i];
    }

    if(1 == 1) {
      PDM_log_trace_array_int (local_post_section_idx, n_section_post+1, "local_post_section_idx ::");
    }

    /*
     * Requilibrate all block
     */
    for(int i_section = 0; i_section < n_section_post; ++i_section) {

      int beg = local_post_section_idx[i_section];
      int end = local_post_section_idx[i_section+1];
      int nl_elmt = end - beg;

      PDM_g_num_t* distrib_elmt = PDM_compute_uniform_entity_distribution(comm, post_section_n[i_section]);
      PDM_g_num_t* ln_to_gn = malloc(local_post_section_n[i_section] * sizeof(PDM_g_num_t));

      for(int i = 0; i < local_post_section_n[i_section]; ++i) {
        ln_to_gn[i] = distrib_face[i_rank] + local_post_section_idx[i_section] + i + 1;
      }

      PDM_part_to_block_t* ptb = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                       PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                                       1.,
                                                                       &ln_to_gn,
                                                                       distrib_elmt,
                                                                       &local_post_section_n[i_section],
                                                                       1,
                                                                       comm);
      int n_face_vtx_tot = pface_vtx_idx[end] - pface_vtx_idx[beg];
      int         *send_face_vtx_n = malloc(nl_elmt        * sizeof(int        ));
      PDM_g_num_t *send_face_vtx   = malloc(n_face_vtx_tot * sizeof(PDM_g_num_t));
      int         *blk_face_vtx_n  = NULL;
      PDM_g_num_t *blk_face_vtx    = NULL;

      int idx_write = 0;
      for(int i = 0; i < nl_elmt; ++i) {
        send_face_vtx_n[i] = pface_vtx_idx[beg+i+1] - pface_vtx_idx[beg+i];
        for(int j = pface_vtx_idx[beg+i]; j < pface_vtx_idx[beg+i+1]; ++j) {
          int i_vtx = pface_vtx[j];
          send_face_vtx[idx_write++] = pvtx_ln_to_gn[i_vtx-1];
        }
      }

      PDM_part_to_block_exch(ptb,
                             sizeof(PDM_g_num_t),
                             PDM_STRIDE_VAR_INTERLACED,
                             -1,
               (int  **)     &send_face_vtx_n,
               (void **)     &send_face_vtx,
               (int  **)     &blk_face_vtx_n,
               (void **)     &blk_face_vtx);

      free(send_face_vtx_n);
      free(send_face_vtx);


      PDM_Mesh_nodal_elt_t t_elt = (PDM_Mesh_nodal_elt_t) post_section_kind[i_section];
      int id_section = PDM_DMesh_nodal_section_add(dmn,
                                                   PDM_GEOMETRY_KIND_SURFACIC,
                                                   t_elt);

      int dn_elmt = distrib_elmt[i_rank+1] - distrib_elmt[i_rank];
      if(t_elt == PDM_MESH_NODAL_POLY_2D) {

        int *blk_face_vtx_idx = PDM_array_new_idx_from_sizes_int(blk_face_vtx_n, dn_elmt);
        PDM_DMesh_nodal_section_poly2d_set(dmn,
                                           PDM_GEOMETRY_KIND_SURFACIC,
                                           id_section,
                                           dn_elmt,
                                           blk_face_vtx_idx,
                                           blk_face_vtx,
                                           PDM_OWNERSHIP_KEEP);
      } else {
        PDM_DMesh_nodal_section_std_set(dmn,
                                        PDM_GEOMETRY_KIND_SURFACIC,
                                        id_section,
                                        dn_elmt,
                                        blk_face_vtx,
                                        PDM_OWNERSHIP_KEEP);
      }

      free(blk_face_vtx_n);
      // free(blk_face_vtx);

      PDM_part_to_block_free(ptb);
      free(distrib_elmt);
      free(ln_to_gn);

    }

    free(post_section_kind);
    free(post_section_n   );
    free(post_section_idx );
    free(local_post_section_n  );
    free(local_post_section_idx);

    if(1 == 0) {
      printf("n_section_tri    = %i\n", n_section_tri   );
      printf("n_section_quad   = %i\n", n_section_quad  );
      printf("n_section_poly2d = %i\n", n_section_poly2d);
    }

    /*
     * Recuperation des bords
     */
    int          n_edge_bound   = n_bound   [PDM_BOUND_TYPE_EDGE];
    int         *edge_bound_idx = dbound_idx[PDM_BOUND_TYPE_EDGE];
    PDM_g_num_t *edge_bound     = dbound    [PDM_BOUND_TYPE_EDGE];
    // int n_section_bar = n_edge_bound;

    int n_edge_bnd_tot = edge_bound_idx[n_edge_bound];
    int          pn_vtx_bnd = 0;
    int         *pedge_vtx_bnd_idx = NULL;
    int         *pedge_bnd_vtx     = NULL;
    PDM_g_num_t *pvtx_bnd_ln_to_gn = NULL;

    PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                             distrib_edge,
                                                             dedge_vtx_idx,
                                                             dedge_vtx,
                                                             n_edge_bnd_tot,
                                                             edge_bound,
                                                             &pn_vtx_bnd,
                                                             &pvtx_bnd_ln_to_gn,
                                                             &pedge_vtx_bnd_idx,
                                                             &pedge_bnd_vtx);

    /*
     * On ne veut pas changer la numerotation absolu des bar
     */

    PDM_Mesh_nodal_elt_t t_elt = PDM_MESH_NODAL_BAR2;
    int id_section = PDM_DMesh_nodal_section_add(dmn,
                                                 PDM_GEOMETRY_KIND_RIDGE,
                                                 t_elt);

    PDM_g_num_t *pedge_bnd_vtx_gnum = malloc(2 * n_edge_bnd_tot * sizeof(PDM_g_num_t));
    for(int i = 0; i < 2 * n_edge_bnd_tot; ++i) {
      pedge_bnd_vtx_gnum[i] = pvtx_bnd_ln_to_gn[pedge_bnd_vtx[i]-1];
    }
    PDM_DMesh_nodal_section_std_set(dmn,
                                    PDM_GEOMETRY_KIND_RIDGE,
                                    id_section,
                                    n_edge_bnd_tot,
                                    pedge_bnd_vtx_gnum,
                                    PDM_OWNERSHIP_KEEP);


    PDM_g_num_t* distrib_bar = PDM_compute_entity_distribution(comm, n_edge_bnd_tot);
    int         *out_edge_bound_idx = malloc((n_edge_bound+1) * sizeof(int        ));
    PDM_g_num_t *out_edge_bound     = malloc( n_edge_bnd_tot  * sizeof(PDM_g_num_t));

    out_edge_bound_idx[0] = edge_bound_idx[0];
    for(int i_group = 0; i_group < n_edge_bound; ++i_group) {
      out_edge_bound_idx[i_group+1] = edge_bound_idx[i_group+1];
      for(int idx_edge = edge_bound_idx[i_group]; idx_edge < edge_bound_idx[i_group+1]; ++idx_edge) {
        out_edge_bound[idx_edge] = distrib_bar[i_rank] + idx_edge + 1;
      }
    }

    PDM_DMesh_nodal_section_group_elmt_set(dmn,
                                           PDM_GEOMETRY_KIND_RIDGE,
                                           n_edge_bound,
                                           out_edge_bound_idx,
                                           out_edge_bound,
                                           PDM_OWNERSHIP_KEEP);

    PDM_log_trace_connectivity_long(out_edge_bound_idx, out_edge_bound, n_edge_bound, "out_edge_bound ::");

    /*
     *  Il faut créer un table qui permet d'updater les numero de faces/edges en numero d'element
     *     - dedge_to_elmts_n
     *     - dedge_to_elmts
     *  Donc c'est un ptb partial
     *    -> Si on a une liste de edge_ln_to_gn
     *  PDM_block_to_part_create_from_sparse_block
     *
     */
    PDM_part_to_block_t* ptb_edge = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                             PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                             1.,
                                                             &edge_bound,
                                                             NULL,
                                                             &n_edge_bnd_tot,
                                                             1,
                                                             comm);
    PDM_g_num_t *gnum_edge = PDM_part_to_block_block_gnum_get(ptb_edge);
    int n_blk_edge = PDM_part_to_block_n_elt_block_get(ptb_edge);

    PDM_g_num_t *blk_out_edge_bound = NULL;
    PDM_part_to_block_exch(ptb_edge,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST_INTERLACED,
                           1,
                           NULL,
                (void **)  &out_edge_bound,
                           NULL,
                (void **)  &blk_out_edge_bound);


    PDM_log_trace_array_long(gnum_edge, n_blk_edge, "gnum_edge ::");
    PDM_log_trace_array_long(blk_out_edge_bound, n_blk_edge, "blk_out_edge_bound ::");

    /* Store it */
    assert(n_blk_gnum     [PDM_BOUND_TYPE_EDGE] == 0);
    assert(blk_entity_gnum[PDM_BOUND_TYPE_EDGE] == NULL);
    assert(blk_elmt_gnum  [PDM_BOUND_TYPE_EDGE] == NULL);

    n_blk_gnum     [PDM_BOUND_TYPE_EDGE] = n_blk_edge;
    blk_entity_gnum[PDM_BOUND_TYPE_EDGE] = malloc(n_blk_edge * sizeof(PDM_g_num_t));
    blk_elmt_gnum  [PDM_BOUND_TYPE_EDGE] = blk_out_edge_bound;

    for(int i = 0; i < n_blk_edge; ++i) {
      blk_elmt_gnum  [PDM_BOUND_TYPE_EDGE][i] = gnum_edge[i];
    }

    /*
     *  Translation face -> elmts
     */
    PDM_block_to_part_t *btp = PDM_block_to_part_create_from_sparse_block(gnum_edge,
                                                                          n_blk_edge,
                                              (const PDM_g_num_t **)      &edge_bound,
                                                                          &n_edge_bnd_tot,
                                                                          1,
                                                                          comm);

    PDM_part_to_block_free(ptb_edge);

    int stride_one = 1;
    PDM_g_num_t **tmp_edge_to_elmt = NULL;
    PDM_block_to_part_exch(btp,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST_INTERLACED,
                           &stride_one,
                           blk_out_edge_bound,
                           NULL,
                (void ***) &tmp_edge_to_elmt);
    PDM_g_num_t *edge_to_elmt = tmp_edge_to_elmt[0];
    free(tmp_edge_to_elmt);

    PDM_log_trace_array_long(edge_to_elmt, n_edge_bnd_tot, "edge_to_elmt ::");


    free(edge_to_elmt);
    // free(blk_out_edge_bound);

    PDM_block_to_part_free(btp);


    free(distrib_bar);

    free(pedge_vtx_bnd_idx);
    free(pedge_bnd_vtx    );
    free(pvtx_bnd_ln_to_gn);
    free(dedge_vtx_idx);

    free(pface_ln_to_gn);
    free(pvtx_ln_to_gn);
    free(pface_vtx_idx);
    free(pface_vtx    );

  }


  return dmn;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_dmesh_to_dmesh_nodal_t*
PDM_dmesh_to_dmesh_nodal_create
(
 const int             n_mesh,
 const PDM_MPI_Comm    comm
)
{

  PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn = (PDM_dmesh_to_dmesh_nodal_t *) malloc(sizeof(PDM_dmesh_to_dmesh_nodal_t));

  dm_to_dmn->comm              = comm;
  dm_to_dmn->results_is_getted = PDM_FALSE;
  dm_to_dmn->n_mesh            = n_mesh;

  dm_to_dmn->dcell_face     = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->dcell_face_idx = malloc(n_mesh * sizeof(int         *));
  dm_to_dmn->dface_edge     = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->dface_edge_idx = malloc(n_mesh * sizeof(int         *));
  dm_to_dmn->dface_vtx      = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->dface_vtx_idx  = malloc(n_mesh * sizeof(int         *));
  dm_to_dmn->dedge_vtx      = malloc(n_mesh * sizeof(PDM_g_num_t *));

  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    dm_to_dmn->dcell_face    [i_mesh] = NULL;
    dm_to_dmn->dcell_face_idx[i_mesh] = NULL;
    dm_to_dmn->dface_edge    [i_mesh] = NULL;
    dm_to_dmn->dface_edge_idx[i_mesh] = NULL;
    dm_to_dmn->dface_vtx     [i_mesh] = NULL;
    dm_to_dmn->dface_vtx_idx [i_mesh] = NULL;
    dm_to_dmn->dedge_vtx     [i_mesh] = NULL;
  }

  dm_to_dmn->distrib_cell = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->distrib_face = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->distrib_edge = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->distrib_vtx  = malloc(n_mesh * sizeof(PDM_g_num_t *));

  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    dm_to_dmn->distrib_cell[i_mesh] = NULL;
    dm_to_dmn->distrib_face[i_mesh] = NULL;
    dm_to_dmn->distrib_edge[i_mesh] = NULL;
    dm_to_dmn->distrib_vtx [i_mesh] = NULL;
  }

  dm_to_dmn->n_bound         = malloc( n_mesh * sizeof(int          *) );
  dm_to_dmn->dbound          = malloc( n_mesh * sizeof(PDM_g_num_t **) );
  dm_to_dmn->dbound_idx      = malloc( n_mesh * sizeof(int         **) );


  dm_to_dmn->n_blk_gnum      = malloc( n_mesh * sizeof(int          *) );
  dm_to_dmn->blk_entity_gnum = malloc( n_mesh * sizeof(PDM_g_num_t **) );
  dm_to_dmn->blk_elmt_gnum   = malloc( n_mesh * sizeof(PDM_g_num_t **) );

  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {

    dm_to_dmn->n_bound   [i_mesh] = malloc( PDM_BOUND_TYPE_MAX * sizeof(int          ) );
    dm_to_dmn->dbound    [i_mesh] = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_g_num_t *) );
    dm_to_dmn->dbound_idx[i_mesh] = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         *) );

    dm_to_dmn->n_blk_gnum     [i_mesh] = malloc( PDM_BOUND_TYPE_MAX * sizeof(int          ) );
    dm_to_dmn->blk_entity_gnum[i_mesh] = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_g_num_t *) );
    dm_to_dmn->blk_elmt_gnum  [i_mesh] = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_g_num_t *) );

    for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i ) {
      dm_to_dmn->n_bound   [i_mesh][i] = 0;
      dm_to_dmn->dbound    [i_mesh][i] = NULL;
      dm_to_dmn->dbound_idx[i_mesh][i] = NULL;

      dm_to_dmn->n_blk_gnum     [i_mesh][i] = 0;
      dm_to_dmn->blk_entity_gnum[i_mesh][i] = NULL;
      dm_to_dmn->blk_elmt_gnum  [i_mesh][i] = NULL;
    }
  }

  dm_to_dmn->dmn = malloc( n_mesh * sizeof(PDM_dmesh_nodal_t *) );
  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    dm_to_dmn->dmn[i_mesh] = NULL;
  }

  return dm_to_dmn;
}

void
PDM_dmesh_to_dmesh_nodal_compute
(
 PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn
)
{
  for(int i_mesh = 0; i_mesh < dm_to_dmn->n_mesh; ++i_mesh) {
    dm_to_dmn->dmn[i_mesh] = _dmesh_to_dmesh_nodal(dm_to_dmn->comm,
                                                   dm_to_dmn->distrib_cell   [i_mesh],
                                                   dm_to_dmn->distrib_face   [i_mesh],
                                                   dm_to_dmn->distrib_edge   [i_mesh],
                                                   dm_to_dmn->distrib_vtx    [i_mesh],
                                                   dm_to_dmn->dcell_face     [i_mesh],
                                                   dm_to_dmn->dcell_face_idx [i_mesh],
                                                   dm_to_dmn->dface_edge     [i_mesh],
                                                   dm_to_dmn->dface_edge_idx [i_mesh],
                                                   dm_to_dmn->dface_vtx      [i_mesh],
                                                   dm_to_dmn->dface_vtx_idx  [i_mesh],
                                                   dm_to_dmn->dedge_vtx      [i_mesh],
                                                   dm_to_dmn->n_bound        [i_mesh],
                                                   dm_to_dmn->dbound_idx     [i_mesh],
                                                   dm_to_dmn->dbound         [i_mesh],
                                                   dm_to_dmn->n_blk_gnum     [i_mesh],
                                                   dm_to_dmn->blk_entity_gnum[i_mesh],
                                                   dm_to_dmn->blk_elmt_gnum  [i_mesh]);
  }


}

void
PDM_dmesh_to_dmesh_nodal_set_dmesh
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        PDM_dmesh_t                *dm
)
{
  PDM_UNUSED(dm_to_dmn);
  PDM_UNUSED(i_mesh);
  PDM_UNUSED(dm);
  abort();
}


void
PDM_dmesh_to_dmesh_nodal_dmesh_nodal_get
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        PDM_dmesh_nodal_t         **dmn,
        PDM_ownership_t             ownership
)
{
  if(ownership == PDM_OWNERSHIP_USER) {
    dm_to_dmn->results_is_getted = PDM_TRUE;
  }
  *dmn = dm_to_dmn->dmn[i_mesh];
}

void
PDM_dmesh_to_dmesh_nodal_distribution_set
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        PDM_g_num_t                *distrib_cell,
        PDM_g_num_t                *distrib_face,
        PDM_g_num_t                *distrib_edge,
        PDM_g_num_t                *distrib_vtx
)
{
  dm_to_dmn->distrib_cell[i_mesh] = distrib_cell;
  dm_to_dmn->distrib_face[i_mesh] = distrib_face;
  dm_to_dmn->distrib_edge[i_mesh] = distrib_edge;
  dm_to_dmn->distrib_vtx [i_mesh] = distrib_vtx;
}

void
PDM_dmesh_to_dmesh_nodal_connectivity_set
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        int                        *dcell_face_idx,
        PDM_g_num_t                *dcell_face,
        int                        *dface_edge_idx,
        PDM_g_num_t                *dface_edge,
        PDM_g_num_t                *dedge_vtx,
        int                        *dface_vtx_idx,
        PDM_g_num_t                *dface_vtx
)
{
  dm_to_dmn->dcell_face_idx[i_mesh] = dcell_face_idx;
  dm_to_dmn->dcell_face    [i_mesh] = dcell_face;
  dm_to_dmn->dface_edge    [i_mesh] = dface_edge;
  dm_to_dmn->dface_edge_idx[i_mesh] = dface_edge_idx;
  dm_to_dmn->dedge_vtx     [i_mesh] = dedge_vtx;
  dm_to_dmn->dface_vtx_idx [i_mesh] = dface_vtx_idx;
  dm_to_dmn->dface_vtx     [i_mesh] = dface_vtx;
}

void
PDM_dmesh_to_dmesh_nodal_group_set
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        PDM_bound_type_t            bound_type,
        int                         n_group,
        int                        *dbound_idx,
        PDM_g_num_t                *dbound
)
{
  dm_to_dmn->n_bound   [i_mesh][bound_type] = n_group;
  dm_to_dmn->dbound_idx[i_mesh][bound_type] = dbound_idx;
  dm_to_dmn->dbound    [i_mesh][bound_type] = dbound;
}


void
PDM_dmesh_to_dmesh_nodal_free
(
 PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn
)
{
  if(dm_to_dmn == NULL) {
    return;
  }

  free(dm_to_dmn->dcell_face    );
  free(dm_to_dmn->dcell_face_idx);
  free(dm_to_dmn->dface_edge    );
  free(dm_to_dmn->dface_edge_idx);
  free(dm_to_dmn->dface_vtx     );
  free(dm_to_dmn->dface_vtx_idx );
  free(dm_to_dmn->dedge_vtx     );

  free(dm_to_dmn->distrib_cell);
  free(dm_to_dmn->distrib_face);
  free(dm_to_dmn->distrib_edge);
  free(dm_to_dmn->distrib_vtx );

  for(int i_mesh = 0; i_mesh < dm_to_dmn->n_mesh; ++i_mesh ) {
    free(dm_to_dmn->n_bound   [i_mesh]);
    free(dm_to_dmn->dbound    [i_mesh]);
    free(dm_to_dmn->dbound_idx[i_mesh]);

    for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i){
      if(dm_to_dmn->blk_entity_gnum[i_mesh][i] != NULL){
        free(dm_to_dmn->blk_entity_gnum[i_mesh][i]);
      }
      if(dm_to_dmn->blk_elmt_gnum[i_mesh][i] != NULL){
        free(dm_to_dmn->blk_elmt_gnum[i_mesh][i]);
      }
    }

    free(dm_to_dmn->n_blk_gnum     [i_mesh]);
    free(dm_to_dmn->blk_entity_gnum[i_mesh]);
    free(dm_to_dmn->blk_elmt_gnum  [i_mesh]);

  }

  free(dm_to_dmn->n_blk_gnum     );
  free(dm_to_dmn->blk_entity_gnum);
  free(dm_to_dmn->blk_elmt_gnum  );

  free(dm_to_dmn->n_bound      );
  free(dm_to_dmn->dbound       );
  free(dm_to_dmn->dbound_idx   );

  if(dm_to_dmn->results_is_getted == PDM_FALSE) {
    for(int i_mesh = 0; i_mesh < dm_to_dmn->n_mesh; ++i_mesh) {
      PDM_DMesh_nodal_free(dm_to_dmn->dmn[i_mesh]);
    }
  }

  free(dm_to_dmn->dmn);

  free(dm_to_dmn);
  dm_to_dmn = NULL;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
