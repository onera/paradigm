#ifndef __PDM_PART_EXTENSION_PRIV_H__
#define __PDM_PART_EXTENSION_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

// #include "pdm_multipart.h"
#include "pdm_part_priv.h"
#include "pdm_part_domain_interface.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _pdm_part_extension_t
 * \brief  Distributed cube
 *
 * _dcube_t define a distributed mesh of a cube
 *
 */

struct _pdm_part_extension_t {
  PDM_MPI_Comm      comm;            /*!< MPI communicator                          */
  PDM_ownership_t   owner;           /*!< Which have the responsabilities of results*/

  int               dim;

  PDM_extend_type_t extend_type;
  int               depth;

  int             n_domain;
  int            *n_part;
  int            *n_part_idx;
  int            *n_part_g_idx;
  int            *lpart_to_dom;

  _part_t  **parts;
  PDM_part_domain_interface_t  *pdi;
  int                           n_interface;

  int user_defined_bnd_graph;

  PDM_bool_t has_connectivity[PDM_CONNECTIVITY_TYPE_MAX];

  /* Store for each depth / each domain / each part */
  int **neighbor_idx;
  int **neighbor_desc;
  int **neighbor_interface;
  int  *n_entity_bound;
  int  *n_cell;

  /* Graph of cell */
  int **dist_neighbor_cell_n;
  int **dist_neighbor_cell_idx;
  int **dist_neighbor_cell_desc;
  int **dist_neighbor_cell_interface;
  int **unique_order_dist_neighbor_cell;
  int  *n_unique_order_dist_neighbor_cell;

  int **cell_cell_idx;
  int **cell_cell;

  /* Management of interface */
  PDM_g_num_t **opp_interface_and_gnum_face; /* For each entity in border contains the opposite gnum and interface */
  PDM_g_num_t **opp_interface_and_gnum_edge;
  PDM_g_num_t **opp_interface_and_gnum_vtx;
  int         **cur_interface_face;          /* Works with opp_interface_and_gnum_*, sort in the smae way and give the local number     */
  int         **cur_interface_edge;
  int         **cur_interface_vtx;
  int         **cur_sens_face;          /* Works with opp_interface_and_gnum_*, sort in the smae way and give the sens of entity to receive     */
  int         **cur_sens_edge;
  int         **cur_sens_vtx;
  int          *n_cur_interface_face;
  int          *n_cur_interface_edge;
  int          *n_cur_interface_vtx;

  /* This one is only on the border and contains only border cells */
  int ***cell_cell_extended_idx;
  int ***cell_cell_extended_n;
  int ***cell_cell_extended;
  int ***cell_cell_interface;
  int  **cell_cell_extended_idx2;
  int  **cell_cell_extended2;
  int  **cell_cell_path_itrf_idx;
  int  **cell_cell_path_itrf;
  int ***unique_order_cell_cell_extended;
  int  **n_unique_order_cell_cell_extended;
  int  **border_cell_list;

  int  **cell_cell_extended_pruned_idx;
  int  **cell_cell_extended_pruned;
  int  **cell_cell_interface_pruned;

  int  *n_tot_part_by_domain;

  int **entity_cell_idx;
  int **entity_cell_n;
  int **entity_cell;

  int **entity_cell_opp_idx;
  int **entity_cell_opp_n;
  int **entity_cell_opp;

  /* Internal results */
  int **face_face_extended_idx;
  int **face_face_extended;
  int **face_face_interface;
  int **face_face_path_itrf_idx;
  int **face_face_path_itrf;

  int **edge_edge_extended_idx;
  int **edge_edge_extended;
  int **edge_edge_interface;
  int **edge_edge_path_itrf_idx;
  int **edge_edge_path_itrf;

  int **vtx_vtx_extended_idx;
  int **vtx_vtx_extended;
  int **vtx_vtx_interface;
  int **vtx_vtx_path_itrf_idx;
  int **vtx_vtx_path_itrf;

  /* Results */
  int          *n_cell_border;
  PDM_g_num_t **border_cell_ln_to_gn;
  PDM_g_num_t **border_cell_ln_to_gn_ancstr;
  int         **border_cell_face_idx;
  int         **border_cell_face;

  int          *n_face_border;
  PDM_g_num_t **border_face_ln_to_gn;
  PDM_g_num_t **border_face_ln_to_gn_ancstr;
  int         **border_face_edge_idx;
  int         **border_face_edge;
  int         **border_face_vtx_idx;
  int         **border_face_vtx;
  int         **border_face_group_idx;
  int         **border_face_group;
  PDM_g_num_t **border_face_group_ln_to_gn;

  int          *n_edge_border;
  PDM_g_num_t **border_edge_ln_to_gn;
  PDM_g_num_t **border_edge_ln_to_gn_ancstr;
  int         **border_edge_vtx_idx;
  int         **border_edge_vtx;
  int         **border_edge_group_idx;
  int         **border_edge_group;
  PDM_g_num_t **border_edge_group_ln_to_gn;

  int          *n_vtx_border;
  PDM_g_num_t **border_vtx_ln_to_gn;
  PDM_g_num_t **border_vtx_ln_to_gn_ancstr;
  double      **border_vtx;
  
  /* Ownership */
  PDM_ownership_t ***ownership_border_ln_to_gn;
  PDM_ownership_t ***ownership_border_ln_to_gn_ancstr;
  PDM_ownership_t ***ownership_border_connectivity;
  PDM_ownership_t  **ownership_border_vtx_coord;
  PDM_ownership_t ***ownership_border_group;
  PDM_ownership_t ***ownership_border_graph;
  PDM_ownership_t ***ownership_border_path_itrf;

  /* Shift by domain for all entities */
  PDM_g_num_t *shift_by_domain_cell;
  PDM_g_num_t *shift_by_domain_face;
  PDM_g_num_t *shift_by_domain_edge;
  PDM_g_num_t *shift_by_domain_vtx;
  PDM_g_num_t *shift_by_domain_edge_group;
  PDM_g_num_t *shift_by_domain_face_group;

  /* Composed interface */
  int  n_composed_interface;
  int *composed_interface_idx;
  int *composed_interface;
  PDM_g_num_t *composed_ln_to_gn_sorted;

  int **pdi_neighbor_idx;
  int **pdi_neighbor;

  /* New manner to do part_extension */
  int compute_kind; // 0 : old // 1 : new
  int ln_part_tot;

  PDM_domain_interface_t  *dom_itrf;

  /* TODO - delete */
  // PDM_part_to_block_t    **ptb_itrf[PDM_BOUND_TYPE_MAX];
  // PDM_g_num_t            **opp_gnum[PDM_BOUND_TYPE_MAX];
  // int                    **opp_sens[PDM_BOUND_TYPE_MAX];

  /* New manner */
  int          dentity_itrf_n_blk              [PDM_BOUND_TYPE_MAX];
  PDM_g_num_t *dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_MAX];
  // PDM_g_num_t *dentity_itrf_blk_ancstr         [PDM_BOUND_TYPE_MAX];
  // int         *dentity_itrf_blk_path_itrf_strd [PDM_BOUND_TYPE_MAX];
  int         *dentity_itrf_blk_path_itrf      [PDM_BOUND_TYPE_MAX];
  int         *dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_MAX];
  PDM_g_num_t *dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_MAX];
  int         *dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_MAX];



  int have_edge;
  int have_face;

  int **pinit_entity_bound_to_pentity_bound_idx;
  int **pinit_entity_bound_to_pentity_bound_triplet;
  int **pinit_entity_bound_to_pentity_bound_interface;

  int owner_vtx_part_bound;

};

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_EXTENSION_PRIV_H__ */
