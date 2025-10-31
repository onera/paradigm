
/*----------------------------------------------------------------------------
 *  Standar headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_array.h"
#include "pdm_block_to_part.h"
#include "pdm_dgeom_elem.h"
#include "pdm_distrib.h"
#include "pdm_error.h"
#include "pdm_hilbert.h"
#include "pdm_mem_tool.h"
#include "pdm_part_comm_graph.h"
#include "pdm_part_comm_graph_priv.h"
#include "pdm_part_geom.h"
#include "pdm_part_to_block.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_priv.h"
#include "pdm_sort.h"
#include "pdm.h"

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

#define CHECK_ELT_SIZE(dimension, i_elt, elt_vtx_n) \
  if (dimension==1 && elt_vtx_n != 2) { \
    PDM_error(__FILE__, __LINE__, 0, "1D elements with more than 2 vertices are not supported, element %d has %d\n", i_elt, elt_vtx_n); \
  } \
  else if (dimension==2 && elt_vtx_n != 3 && elt_vtx_n != 4) { \
    PDM_error(__FILE__, __LINE__, 0, "Only 2D elements with 3 or 4 vertices are supported for new, element %d has %d\n", i_elt, elt_vtx_n); \
  }

/*============================================================================
 * Private function definitions
 *============================================================================*/

static inline void
_compute_edge_normal
(
  int     elt_vtx_n,
  int    *elt_vtx,
  double *vtx_coord,
  double *normal
)
{
  PDM_UNUSED(elt_vtx_n);

  int i_vtx0 = elt_vtx[0] - 1;
  int i_vtx1 = elt_vtx[1] - 1;

  normal[0] = vtx_coord[3*i_vtx1+1] - vtx_coord[3*i_vtx0+1];
  normal[1] = vtx_coord[3*i_vtx0  ] - vtx_coord[3*i_vtx1  ];
  normal[2] = 0.;
}


static inline void
_compute_face_normal
(
  int     elt_vtx_n,
  int    *elt_vtx,
  double *vtx_coord,
  double *normal
)
{
  // TODO: generalise to other types of elements besides TRI and QUAD
  for (int i = 0; i < 3; i++) {
    normal[i] = 0.;
  }

  for (int i_tri = 0; i_tri < elt_vtx_n-2; i_tri++) {
    int i_vtx0 = elt_vtx[i_tri      ] - 1;
    int i_vtx1 = elt_vtx[i_tri+1    ] - 1;
    int i_vtx2 = elt_vtx[elt_vtx_n-1] - 1;

    double vec_1[3], vec_2[3];
    for (int i = 0; i < 3; i++) {
      vec_1[i] = vtx_coord[3*i_vtx1+i] - vtx_coord[3*i_vtx0+i];
      vec_2[i] = vtx_coord[3*i_vtx2+i] - vtx_coord[3*i_vtx0+i];
    }

    double tri_normal[3];
    PDM_CROSS_PRODUCT(tri_normal, vec_1, vec_2);

    for (int i = 0; i < 3; i++) {
      normal[i] += 0.5*tri_normal[i];
    }
  }
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

void
PDM_dcompute_cell_center
(
  const PDM_MPI_Comm  comm,
  const int           dn_cell,
  const int          *dcell_face_idx,
  const PDM_g_num_t  *dcell_face,
  const int          *dface_vtx_idx,
  const PDM_g_num_t  *dface_vtx,
  const PDM_g_num_t  *distrib_face,
  const double       *dvtx_coord,
  const PDM_g_num_t  *distrib_vtx,
  double             *cell_center
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  int dn_face = (int) (distrib_face[i_rank+1] - distrib_face[i_rank]);

  PDM_g_num_t *dface_ln_to_gn = NULL;
  PDM_malloc(dface_ln_to_gn, dn_face, PDM_g_num_t);
  for (int i = 0; i < dn_face; i++) {
    dface_ln_to_gn[i] = distrib_face[i_rank] + i + 1;
  }

  PDM_g_num_t *pvtx_ln_to_gn;
  int         *pface_vtx_idx;
  int         *pface_vtx;
  int          pn_vtx;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_face,
                                                           dface_vtx_idx,
                                                           dface_vtx,
                                                           dn_face,
                                                           dface_ln_to_gn,
                                                           &pn_vtx,
                                                           &pvtx_ln_to_gn,
                                                           &pface_vtx_idx,
                                                           &pface_vtx);
  PDM_free(dface_ln_to_gn);

  double** tmp_pvtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        distrib_vtx,
                                        dvtx_coord,
                                        &pn_vtx,
                                        (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                        &tmp_pvtx_coord);
  double *pvtx_coord = tmp_pvtx_coord[0];
  PDM_free(tmp_pvtx_coord);
  PDM_free(pvtx_ln_to_gn);


  /* Compute face centers */
  double *dface_center = NULL;
  PDM_malloc(dface_center, dn_face * 3, double);
  for (int i = 0; i < dn_face; i++) {
    for (int k = 0; k < 3; k++) {
      dface_center[3*i + k] = 0.;
    }

    double normalization = 1. / (double) (pface_vtx_idx[i+1] - pface_vtx_idx[i]);
    for (int j = pface_vtx_idx[i]; j < pface_vtx_idx[i+1]; j++) {
      int ivtx = pface_vtx[j] - 1;

      for (int k = 0; k < 3; k++) {
        dface_center[3*i + k] += pvtx_coord[3*ivtx + k];
      }
    }

    for (int k = 0; k < 3; k++) {
      dface_center[3*i + k] *= normalization;
    }
  }
  PDM_free(pvtx_coord);
  PDM_free(pface_vtx_idx);
  PDM_free(pface_vtx);

  /* Compute cell centers */
  PDM_compute_center_from_descending_connectivity (dcell_face_idx,
                                                   dcell_face,
                                                   dn_cell,
                                                   distrib_face,
                                                   cell_center,
                                                   dface_center,
                                                   comm);
  PDM_free(dface_center);
}


void
PDM_part_entity_geom
(
 PDM_part_geom_t     method,
 const int           n_part,
 const PDM_MPI_Comm  comm,
 const PDM_g_num_t   dn_entity,
 const double       *dentity_coord,
 const double       *dentity_weight,
       int          *dentity_part
)
{
  PDM_UNUSED(method);
  assert (method == PDM_PART_GEOM_HILBERT);

  const int dim = 3;

  /** TRAITEMENT HILBERT FVM **/
  PDM_hilbert_code_t *hilbert_codes = NULL;
  PDM_malloc(hilbert_codes, dn_entity, PDM_hilbert_code_t);

  /** Initialisation **/

  double extents[2*dim]; /** DIM x 2**/

  /** Get EXTENTS **/
  PDM_hilbert_get_coord_extents_par(dim, dn_entity, dentity_coord, extents, comm);

  /** Hilbert Coordinates Computation **/
  PDM_hilbert_encode_coords(dim, PDM_HILBERT_CS, extents, dn_entity, dentity_coord, hilbert_codes);

  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  int n_total_part;
  PDM_MPI_Allreduce ((void *) &n_part, &n_total_part, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  PDM_hilbert_code_t *hilbert_codes_idx;
  PDM_malloc(hilbert_codes_idx, n_total_part+1, PDM_hilbert_code_t);

  double *weight = NULL;
  PDM_malloc(weight, dn_entity, double);
  if (dentity_weight != NULL) {
    for(int i = 0; i < dn_entity; ++i) {
      weight [i] = dentity_weight [i];
    }
  }
  else {
    for(int i = 0; i < dn_entity; ++i) {
      weight [i] = 1.;
    }
  }

  PDM_hilbert_build_rank_index (dim,
                                n_total_part,
                                dn_entity,
                                hilbert_codes,
                                weight,
                                NULL,
                                hilbert_codes_idx,
                                comm);

  PDM_free(weight);

  /** Remplissage de cell_parts -> en fct des codes Hilbert **/

  for(int i = 0; i < dn_entity; ++i) {
    size_t quantile = PDM_hilbert_quantile_search(n_total_part,
                                                hilbert_codes[i],
                                                hilbert_codes_idx);
    dentity_part[i] = (int) quantile;

  }

  PDM_free(hilbert_codes_idx);
  PDM_free(hilbert_codes);
}


void
PDM_part_geom
(
 PDM_part_geom_t     method,
 const int           n_part,
 const PDM_MPI_Comm  comm,
 const int           dn_cell,
 const int          *dcell_face_idx,
 const PDM_g_num_t  *dcell_face,
 const int          *dcell_weight,
 const int          *dface_vtx_idx,
 const PDM_g_num_t  *dface_vtx,
 const PDM_g_num_t  *distrib_face,
 const double       *dvtx_coord,
 const PDM_g_num_t  *distrib_vtx,
       int          *dcell_part
)
{
  assert (method == PDM_PART_GEOM_HILBERT);
  /*
   * cell center computation
   */
  double *barycenter_coords;
  PDM_malloc(barycenter_coords, dn_cell * 3, double );
  PDM_dcompute_cell_center (comm,
                            dn_cell,
                            dcell_face_idx,
                            dcell_face,
                            dface_vtx_idx,
                            dface_vtx,
                            distrib_face,
                            dvtx_coord,
                            distrib_vtx,
                            barycenter_coords);

  double *dcell_weight_d = NULL;
  if(dcell_weight != NULL) {
    PDM_malloc(dcell_weight_d, dn_cell, double);
    for(int i = 0; i < dn_cell; ++i) {
      dcell_weight_d[i] = dcell_weight[i];
    }
  }

  PDM_part_entity_geom(method,
                       n_part,
                       comm,
                       dn_cell,
                       barycenter_coords,
                       dcell_weight_d,
                       dcell_part);

  if(dcell_weight != NULL) {
    PDM_free(dcell_weight_d);
  }

  PDM_free(barycenter_coords);
}


void
PDM_part_geom_0d
(
 PDM_part_geom_t     method,
 const int           n_part,
 const PDM_MPI_Comm  comm,
 const int           dn_vtx,
 const double       *dvtx_coord,
 const double       *dvtx_weight,
       int          *dvtx_part
)
{
  PDM_part_entity_geom(method,
                       n_part,
                       comm,
                       dn_vtx,
                       dvtx_coord,
                       dvtx_weight,
                       dvtx_part);
}

void
PDM_part_geom_1d
(
 PDM_part_geom_t     method,
 const int           n_part,
 const PDM_MPI_Comm  comm,
 const int           dn_edge,
 const int           dn_vtx,
 const PDM_g_num_t  *dedge_vtx,
 const double       *dvtx_coord,
 const double       *dedge_weight,
       int          *dedge_part
)
{
  PDM_g_num_t *distrib_vtx = PDM_compute_entity_distribution(comm, dn_vtx);

  int *dedge_vtx_idx = PDM_array_new_idx_from_const_stride_int(2, dn_edge);

  double *dedge_center;
  PDM_malloc(dedge_center, dn_edge * 3, double);
  PDM_compute_center_from_descending_connectivity(dedge_vtx_idx,
                                                  dedge_vtx,
                                                  dn_edge,
                                                  distrib_vtx,
                                                  dedge_center,
                                  (double *)      dvtx_coord,
                                                  comm);

  PDM_part_entity_geom(method,
                       n_part,
                       comm,
                       dn_edge,
                       dedge_center,
                       dedge_weight,
                       dedge_part);

  PDM_free(distrib_vtx);
  PDM_free(dedge_center);
  PDM_free(dedge_vtx_idx);
}

void
PDM_part_geom_2d
(
 PDM_part_geom_t     method,
 const int           n_part,
 const PDM_MPI_Comm  comm,
 const int           dn_face,
 const int           dn_edge,
 const int           dn_vtx,
 const int          *dface_vtx_idx,
 const PDM_g_num_t  *dface_vtx,
 const int          *dface_edge_idx,
 const PDM_g_num_t  *dface_edge,
 const PDM_g_num_t  *dedge_vtx,
 const double       *dvtx_coord,
 const double       *dface_weight,
       int          *dface_part
)
{

  PDM_g_num_t *distrib_vtx = PDM_compute_entity_distribution(comm, dn_vtx);
  double *dface_center;
  PDM_malloc(dface_center, dn_face * 3, double);

  if(dface_vtx_idx != NULL) {
    PDM_compute_center_from_descending_connectivity(dface_vtx_idx,
                                                    dface_vtx,
                                                    dn_face,
                                                    distrib_vtx,
                                                    dface_center,
                                    (double *)      dvtx_coord,
                                                    comm);

  } else {
    assert(dface_edge_idx != NULL);
    int *dedge_vtx_idx = PDM_array_new_idx_from_const_stride_int(2, dn_edge);

    double *dedge_center;
    PDM_malloc(dedge_center, dn_edge * 3, double);

    PDM_compute_center_from_descending_connectivity(dedge_vtx_idx,
                                                    dedge_vtx,
                                                    dn_edge,
                                                    distrib_vtx,
                                                    dedge_center,
                                    (double *)      dvtx_coord,
                                                    comm);

    PDM_g_num_t *distrib_edge = PDM_compute_entity_distribution(comm, dn_edge);
    PDM_compute_center_from_descending_connectivity(dface_edge_idx,
                                                    dface_edge,
                                                    dn_face,
                                                    distrib_edge,
                                                    dface_center,
                                    (double *)      dedge_center,
                                                    comm);


    PDM_free(dedge_vtx_idx);
    PDM_free(dedge_center);
    PDM_free(distrib_edge);
  }

  PDM_part_entity_geom(method,
                       n_part,
                       comm,
                       dn_face,
                       dface_center,
                       dface_weight,
                       dface_part);

  PDM_free(distrib_vtx);
  PDM_free(dface_center);
}


void
PDM_dreorder_from_coords
(
 PDM_part_geom_t  method,
 int              dim,
 PDM_g_num_t     *distrib_vtx,
 double          *dcoords,
 PDM_g_num_t     *ln_to_gn,
 PDM_MPI_Comm     comm
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  assert (method == PDM_PART_GEOM_HILBERT);

  int dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];

  PDM_hilbert_code_t *hilbert_codes;
  PDM_malloc(hilbert_codes, dn_vtx, PDM_hilbert_code_t);

  /** Initialisation **/
  double extents[2*dim]; /** DIM x 2**/

  /** Get EXTENTS **/
  PDM_hilbert_get_coord_extents_par(dim, dn_vtx, dcoords, extents, comm);

  /** Hilbert Coordinates Computation **/
  PDM_hilbert_encode_coords(dim, PDM_HILBERT_CS, extents, dn_vtx, dcoords, hilbert_codes);

  PDM_hilbert_code_t *hilbert_codes_idx;
  PDM_malloc(hilbert_codes_idx, n_rank+1, PDM_hilbert_code_t);

  double *weight = PDM_array_const_double(dn_vtx, 1.);
  PDM_hilbert_build_rank_index (dim,
                                n_rank,
                                dn_vtx,
                                hilbert_codes,
                                weight,
                                NULL, // No need order
                                hilbert_codes_idx,
                                comm);
  PDM_free(weight);

  /** Remplissage de cell_parts -> en fct des codes Hilbert **/
  for(int i = 0; i < dn_vtx; ++i) {
    size_t quantile = PDM_hilbert_quantile_search(n_rank,
                                                  hilbert_codes[i],
                                                  hilbert_codes_idx);
    ln_to_gn [i] = (int) (quantile + 1); // Because ln_to_gn of part_to_block begin at 1
  }

  // part_to_block avec ln_to_gn 1 2 3 4 .... pdm_assembly_partition
  // Puis on échange les hilbert_codes, qu'on retrie localement

  PDM_g_num_t *distrib_rank = NULL;
  PDM_malloc(distrib_rank, n_rank+1, PDM_g_num_t);
  for(int i = 0; i < n_rank+1; ++i) {
    distrib_rank[i] = i;
  }

  /*
   * Each proc get all the entities affected to its partitions
   */
  PDM_part_to_block_t *ptb =
   PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                             PDM_PART_TO_BLOCK_POST_NOTHING,
                             1.,
                             &ln_to_gn,
                             distrib_rank,
                             &dn_vtx,
                             1,
                             comm);
  PDM_free(distrib_rank);

  const int n_vtx_block = PDM_part_to_block_n_elt_block_get (ptb);

  // log_trace("n_vtx_block = %i | dn_vtx = %i \n", n_vtx_block, dn_vtx);

  double *blk_hilbert_codes = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(double),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void **) &hilbert_codes,
                          NULL,
                (void **) &blk_hilbert_codes);

  /* Resend */
  for(int i = 0; i < dn_vtx; ++i) {
    ln_to_gn [i] = distrib_vtx[i_rank] + i + 1; // Donc correspond a la numeration absolu initiale
  }

  PDM_g_num_t* blk_ln_to_gn;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void **) &ln_to_gn,
                          NULL,
                (void **) &blk_ln_to_gn);

  PDM_free(hilbert_codes_idx);
  PDM_free(hilbert_codes);
  PDM_part_to_block_free(ptb);

  /* Reorder locally */
  assert (sizeof(double) == sizeof(PDM_hilbert_code_t));
  int *hilbert_order = NULL;
  PDM_malloc(hilbert_order, n_vtx_block, int);
  for (int i = 0; i < n_vtx_block; ++i) {
    hilbert_order [i] = i;
  }
  PDM_sort_double (blk_hilbert_codes, hilbert_order, n_vtx_block);
  PDM_free(blk_hilbert_codes);


  /* Apply order to blk_ln_to_gn */
  PDM_g_num_t *sorted_blk_ln_to_gn = NULL;
  PDM_malloc(sorted_blk_ln_to_gn, n_vtx_block, PDM_g_num_t);
  for(int i = 0; i < n_vtx_block; ++i) {
    sorted_blk_ln_to_gn[i] = blk_ln_to_gn[hilbert_order[i]];
  }
  PDM_free(blk_ln_to_gn);
  PDM_free(hilbert_order);

  PDM_g_num_t* distrib_blk_vtx = PDM_compute_entity_distribution(comm, n_vtx_block);
  PDM_block_to_part_t* btp = PDM_block_to_part_create(distrib_blk_vtx,
                              (const PDM_g_num_t **)  &ln_to_gn,
                                                      &dn_vtx,
                                                      1,
                                                      comm);

  int stride_one = 1;
  PDM_block_to_part_exch_in_place(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
                         &stride_one,
              (void *)   sorted_blk_ln_to_gn,
                         NULL,
              (void **) &ln_to_gn);
  PDM_block_to_part_free(btp);
  PDM_free(sorted_blk_ln_to_gn);
  PDM_free(distrib_blk_vtx);
}


void
PDM_dreorder_from_length
(
 int              dim,
 PDM_g_num_t     *distrib_in,
 double          *length,
 PDM_g_num_t     *ln_to_gn,
 PDM_MPI_Comm     comm
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_length = distrib_in[i_rank+1] - distrib_in[i_rank];

  PDM_hilbert_code_t *tmp_hilbert_codes = NULL;
  PDM_malloc(tmp_hilbert_codes, dn_length, PDM_hilbert_code_t);

  for (int i = 0; i < dn_length; ++i) {
    tmp_hilbert_codes [i] = length [i];
  }

  ///** Calcul des index des codes Hilbert **/
  int *hilbert_order = NULL;
  PDM_malloc(hilbert_order, dn_length, int);

  for (int i = 0; i < dn_length; ++i) {
    hilbert_order [i] = i;
  }

  assert (sizeof(double) == sizeof(PDM_hilbert_code_t));
  PDM_sort_double (tmp_hilbert_codes, hilbert_order, dn_length);

  PDM_free(tmp_hilbert_codes);

  PDM_hilbert_code_t *hilbert_codes_idx = NULL;
  PDM_malloc(hilbert_codes_idx, n_rank+1, PDM_hilbert_code_t);

  double *weight = PDM_array_const_double(dn_length, 1.);
  PDM_hilbert_build_rank_index (dim,
                                n_rank,
                                dn_length,
                                length,
                                weight,
                                hilbert_order,
                                hilbert_codes_idx,
                                comm);


  PDM_free(weight);

  /** Remplissage de cell_parts -> en fct des codes Hilbert **/
  for(int i = 0; i < dn_length; ++i) {
    size_t quantile = PDM_hilbert_quantile_search(n_rank,
                                                  length[i],
                                                  hilbert_codes_idx);
    ln_to_gn [i] = (int) (quantile + 1); // Because ln_to_gn of part_to_block begin at 1
  }

  // part_to_block avec ln_to_gn 1 2 3 4 .... pdm_assembly_partition
  // Puis on échange les hilbert_codes, qu'on retrie localement

  PDM_g_num_t *distrib_rank = NULL;
  PDM_malloc(distrib_rank, n_rank+1, PDM_g_num_t);
  for(int i = 0; i < n_rank+1; ++i) {
    distrib_rank[i] = i;
  }

  /*
   * Each proc get all the entities affected to its partitions
   */
  PDM_part_to_block_t *ptb =
   PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                             PDM_PART_TO_BLOCK_POST_NOTHING,
                             1.,
                             &ln_to_gn,
                             distrib_rank,
                             &dn_length,
                             1,
                             comm);
  PDM_free(distrib_rank);

  const int n_vtx_block = PDM_part_to_block_n_elt_block_get (ptb);

  // log_trace("n_vtx_block = %i | dn_length = %i \n", n_vtx_block, dn_length);

  double *blk_hilbert_codes = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(double),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void **) &length,
                          NULL,
                (void **) &blk_hilbert_codes);

  /* Resend */
  for(int i = 0; i < dn_length; ++i) {
    ln_to_gn [i] = distrib_in[i_rank] + i + 1; // Donc correspond a la numeration absolu initiale
  }

  PDM_g_num_t* blk_ln_to_gn;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void **) &ln_to_gn,
                          NULL,
                (void **) &blk_ln_to_gn);

  PDM_free(hilbert_codes_idx);
  PDM_part_to_block_free(ptb);

  /* Reorder locally */
  assert (sizeof(double) == sizeof(PDM_hilbert_code_t));
  PDM_realloc(hilbert_order ,hilbert_order ,  n_vtx_block ,int);
  for (int i = 0; i < n_vtx_block; ++i) {
    hilbert_order [i] = i;
  }
  PDM_sort_double (blk_hilbert_codes, hilbert_order, n_vtx_block);
  //PDM_log_trace_array_double(blk_hilbert_codes, n_vtx_block, "tmp_edge_length : ");
  PDM_free(blk_hilbert_codes);


  /* Apply order to blk_ln_to_gn */
  PDM_g_num_t *sorted_blk_ln_to_gn = NULL;
  PDM_malloc(sorted_blk_ln_to_gn, n_vtx_block, PDM_g_num_t);
  for(int i = 0; i < n_vtx_block; ++i) {
    sorted_blk_ln_to_gn[i] = blk_ln_to_gn[hilbert_order[i]];
  }
  PDM_free(blk_ln_to_gn);
  PDM_free(hilbert_order);

  PDM_g_num_t* distrib_blk_vtx = PDM_compute_entity_distribution(comm, n_vtx_block);
  PDM_block_to_part_t* btp = PDM_block_to_part_create(distrib_blk_vtx,
                              (const PDM_g_num_t **)  &ln_to_gn,
                                                      &dn_length,
                                                      1,
                                                      comm);

  int stride_one = 1;
  PDM_block_to_part_exch_in_place(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
                         &stride_one,
              (void *)   sorted_blk_ln_to_gn,
                         NULL,
              (void **) &ln_to_gn);
  PDM_block_to_part_free(btp);
  PDM_free(sorted_blk_ln_to_gn);
  PDM_free(distrib_blk_vtx);
}


void
PDM_part_geom_edge_center
(
  int       n_part,
  int      *n_selected,
  int     **selected_lnum,
  int     **pedge_vtx,
  double  **pvtx_coord,
  double ***edge_center
)
{
  for(int i_part = 0; i_part < n_part; ++i_part) {
    assert(pedge_vtx [i_part] != NULL);
    assert(pvtx_coord[i_part] != NULL);
  }

  double **entity_center;
  PDM_malloc(entity_center, n_part, double *);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_malloc(entity_center[i_part], 3 * n_selected[i_part], double);

    double *_pvtx_coord = pvtx_coord[i_part];
    int    *_pedge_vtx  = pedge_vtx [i_part];

    for(int idx_edge = 0; idx_edge < n_selected[i_part]; ++idx_edge) {
      int i_edge = idx_edge;
      if (selected_lnum != NULL) {
        i_edge = selected_lnum[i_part][idx_edge]-1;
      }
      int i_vtx1 = _pedge_vtx[2*i_edge  ]-1;
      int i_vtx2 = _pedge_vtx[2*i_edge+1]-1;
      entity_center[i_part][3*idx_edge  ] = 0.5 * (_pvtx_coord[3*i_vtx1  ] + _pvtx_coord[3*i_vtx2  ]);
      entity_center[i_part][3*idx_edge+1] = 0.5 * (_pvtx_coord[3*i_vtx1+1] + _pvtx_coord[3*i_vtx2+1]);
      entity_center[i_part][3*idx_edge+2] = 0.5 * (_pvtx_coord[3*i_vtx1+2] + _pvtx_coord[3*i_vtx2+2]);
    }
  }
  *edge_center = entity_center;
}


void
PDM_part_geom_face_center_from_edge
(
  int       n_part,
  int      *n_selected,
  int     **selected_lnum,
  int     **pface_edge_idx,
  int     **pface_edge,
  int     **pedge_vtx,
  double  **pvtx_coord,
  double ***face_center
)
{
  for(int i_part = 0; i_part < n_part; ++i_part) {
    assert(pface_edge    [i_part] != NULL);
    assert(pface_edge_idx[i_part] != NULL);
    assert(pedge_vtx     [i_part] != NULL);
    assert(pvtx_coord    [i_part] != NULL);
  }

  double **entity_center;
  PDM_malloc(entity_center, n_part, double * );
  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_malloc(entity_center[i_part], 3 * n_selected[i_part], double);

    int    *_pface_edge     = pface_edge    [i_part];
    int    *_pface_edge_idx = pface_edge_idx[i_part];
    int    *_pedge_vtx      = pedge_vtx     [i_part];
    double *_pvtx_coord     = pvtx_coord    [i_part];

    for(int idx_face = 0; idx_face < n_selected[i_part]; ++idx_face) {

      int i_face = idx_face;
      if (selected_lnum != NULL) {
        i_face = selected_lnum[i_part][idx_face]-1;
      }
      entity_center[i_part][3*idx_face  ] = 0.;
      entity_center[i_part][3*idx_face+1] = 0.;
      entity_center[i_part][3*idx_face+2] = 0.;

      double inv = 1./((double) _pface_edge_idx[i_face+1] - _pface_edge_idx[i_face]);

      for(int idx_edge = _pface_edge_idx[i_face]; idx_edge < _pface_edge_idx[i_face+1]; ++idx_edge) {
        int i_edge = PDM_ABS(_pface_edge[idx_edge])-1;
        int i_vtx1 = _pedge_vtx[2*i_edge  ] - 1;
        int i_vtx2 = _pedge_vtx[2*i_edge+1] - 1;

        entity_center[i_part][3*idx_face  ] += 0.5 * (_pvtx_coord[3*i_vtx1  ] + _pvtx_coord[3*i_vtx2  ]);
        entity_center[i_part][3*idx_face+1] += 0.5 * (_pvtx_coord[3*i_vtx1+1] + _pvtx_coord[3*i_vtx2+1]);
        entity_center[i_part][3*idx_face+2] += 0.5 * (_pvtx_coord[3*i_vtx1+2] + _pvtx_coord[3*i_vtx2+2]);

      }
      entity_center[i_part][3*idx_face  ] = entity_center[i_part][3*idx_face  ] * inv;
      entity_center[i_part][3*idx_face+1] = entity_center[i_part][3*idx_face+1] * inv;
      entity_center[i_part][3*idx_face+2] = entity_center[i_part][3*idx_face+2] * inv;
    }
  }

  *face_center = entity_center;
}


void
PDM_part_geom_face_center_from_vtx
(
  int       n_part,
  int      *n_selected,
  int     **selected_lnum,
  int     **pface_vtx_idx,
  int     **pface_vtx,
  double  **pvtx_coord,
  double ***face_center
)
{
  for(int i_part = 0; i_part < n_part; ++i_part) {
    assert(pface_vtx    [i_part] != NULL);
    assert(pface_vtx_idx[i_part] != NULL);
    assert(pvtx_coord   [i_part] != NULL);
  }

  double **entity_center;
  PDM_malloc(entity_center, n_part, double * );
  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_malloc(entity_center[i_part], 3 * n_selected[i_part], double);

    int    *_pface_vtx     = pface_vtx    [i_part];
    int    *_pface_vtx_idx = pface_vtx_idx[i_part];
    double *_pvtx_coord    = pvtx_coord   [i_part];

    for(int idx_face = 0; idx_face < n_selected[i_part]; ++idx_face) {

      int i_face = idx_face;
      if (selected_lnum != NULL) {
        i_face = selected_lnum[i_part][idx_face]-1;
      }
      entity_center[i_part][3*idx_face  ] = 0.;
      entity_center[i_part][3*idx_face+1] = 0.;
      entity_center[i_part][3*idx_face+2] = 0.;

      double inv = 1./((double) _pface_vtx_idx[i_face+1] - _pface_vtx_idx[i_face]);

      for(int idx_vtx = _pface_vtx_idx[i_face]; idx_vtx < _pface_vtx_idx[i_face+1]; ++idx_vtx) {
        int i_vtx = _pface_vtx[idx_vtx] - 1;

        entity_center[i_part][3*idx_face  ] += _pvtx_coord[3*i_vtx  ];
        entity_center[i_part][3*idx_face+1] += _pvtx_coord[3*i_vtx+1];
        entity_center[i_part][3*idx_face+2] += _pvtx_coord[3*i_vtx+2];

      }
      entity_center[i_part][3*idx_face  ] = entity_center[i_part][3*idx_face  ] * inv;
      entity_center[i_part][3*idx_face+1] = entity_center[i_part][3*idx_face+1] * inv;
      entity_center[i_part][3*idx_face+2] = entity_center[i_part][3*idx_face+2] * inv;
    }
  }

  *face_center = entity_center;
}


void
PDM_part_geom_cell_center
(
  int       n_part,
  int      *n_selected,
  int     **selected_lnum,
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
  int from_edge = (pface_edge != NULL && pface_edge_idx != NULL);
  int from_face = (pface_vtx  != NULL && pface_vtx_idx  != NULL);

  for (int i_part = 0; i_part < n_part; ++i_part) {
    if (n_selected[i_part] > 0) {
      if (from_edge) {
        if (pface_edge[i_part] == NULL || pface_edge_idx[i_part] == NULL) {
          from_edge = 0;
        }
      }
      if (from_face) {
        if (pface_vtx [i_part] == NULL || pface_vtx_idx [i_part] == NULL) {
          from_face = 0;
        }
      }
    }
  }

  if (!from_edge && !from_face) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_geom_cell_center: either face->vtx or face->edge connectivity must be provided\n");
  }

  double **entity_center;
  PDM_malloc(entity_center, n_part, double *);

  if(from_face == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_malloc(entity_center[i_part], 3 * n_selected[i_part], double);

      int    *_pcell_face     = pcell_face    [i_part];
      int    *_pcell_face_idx = pcell_face_idx[i_part];
      int    *_pface_vtx      = pface_vtx     [i_part];
      int    *_pface_vtx_idx  = pface_vtx_idx [i_part];
      double *_pvtx_coord     = pvtx_coord    [i_part];

      // PDM_log_trace_array_int(selected_lnum[i_part], n_selected[i_part], "selected_lnum ::");
      for(int idx_cell = 0; idx_cell < n_selected[i_part]; ++idx_cell) {
        int i_cell = idx_cell;
        if (selected_lnum != NULL) {
          i_cell = selected_lnum[i_part][idx_cell]-1;
        }
        entity_center[i_part][3*idx_cell  ] = 0.;
        entity_center[i_part][3*idx_cell+1] = 0.;
        entity_center[i_part][3*idx_cell+2] = 0.;

        double inv = 1./((double) _pcell_face_idx[idx_cell+1] - _pcell_face_idx[idx_cell]);

        for(int idx_face = _pcell_face_idx[i_cell]; idx_face < _pcell_face_idx[i_cell+1]; ++idx_face) {
          int i_face = PDM_ABS(_pcell_face[idx_face])-1;

          double fcx = 0;
          double fcy = 0;
          double fcz = 0;
          double inv2 = 1./((double) _pface_vtx_idx[i_face+1] - _pface_vtx_idx[i_face]);

          for(int idx_vtx = _pface_vtx_idx[i_face]; idx_vtx < _pface_vtx_idx[i_face+1]; ++idx_vtx) {
            int i_vtx = _pface_vtx[idx_vtx]-1;
            fcx += _pvtx_coord[3*i_vtx  ];
            fcy += _pvtx_coord[3*i_vtx+1];
            fcz += _pvtx_coord[3*i_vtx+2];
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

  else if( from_edge == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_malloc(entity_center[i_part], 3 * n_selected[i_part], double);

      int    *_pcell_face     = pcell_face    [i_part];
      int    *_pcell_face_idx = pcell_face_idx[i_part];
      int    *_pface_edge     = pface_edge    [i_part];
      int    *_pface_edge_idx = pface_edge_idx[i_part];
      int    *_pedge_vtx      = pedge_vtx     [i_part];
      double *_pvtx_coord     = pvtx_coord    [i_part];

      for(int idx_cell = 0; idx_cell < n_selected[i_part]; ++idx_cell) {
        int i_cell = idx_cell;
        if (selected_lnum != NULL) {
          i_cell = selected_lnum[i_part][idx_cell]-1;
        }

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
            fcx += 0.5 * (_pvtx_coord[3*i_vtx1  ] + _pvtx_coord[3*i_vtx2  ]);
            fcy += 0.5 * (_pvtx_coord[3*i_vtx1+1] + _pvtx_coord[3*i_vtx2+1]);
            fcz += 0.5 * (_pvtx_coord[3*i_vtx1+2] + _pvtx_coord[3*i_vtx2+2]);
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


void
PDM_part_geom_vtx_normal_compute
(
  PDM_MPI_Comm             comm,
  int                      n_part,
  int                      dimension,
  int                     *n_selected_elt,
  int                    **selected_elt,
  int                    **elt_vtx_idx,
  int                    **elt_vtx,
  PDM_part_comm_graph_t   *pcg_elt,
  int                     *n_selected_vtx,
  int                    **selected_vtx,
  double                 **vtx_coord,
  PDM_part_comm_graph_t   *pcg_vtx,
  double                ***out_selected_vtx_normal
)
{
  /**
   * Make sure both pcg have the same communicator and n_part
   */
  if (pcg_elt->n_part != pcg_vtx->n_part ||
               n_part != pcg_vtx->n_part ||
      pcg_elt->n_part !=          n_part) {
    PDM_error(__FILE__, __LINE__, 0, "n_part, pcg_elt->n_part and pcg_vtx->n_part are different (%d, %d and %d)\n",
      n_part,
      pcg_elt->n_part,
      pcg_vtx->n_part
    );
  }
  int is_same_comm = 0;
  MPI_Comm_compare(pcg_elt->comm, pcg_vtx->comm, &is_same_comm);
  if (is_same_comm != MPI_IDENT) {
    PDM_error(__FILE__, __LINE__, 0, "pcg_elt and pcg_vtx has different comm\n");
  }
  MPI_Comm_compare(comm, pcg_vtx->comm, &is_same_comm);
  if (is_same_comm != MPI_IDENT) {
    PDM_error(__FILE__, __LINE__, 0, "comm and pcg_vtx->comm are different\n");
  }

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  /* Compute normals */
  void (*_compute_elt_normal) (int, int *, double *, double *) = NULL;
  if (dimension == 1) {
    _compute_elt_normal = &_compute_edge_normal;
  }
  else if (dimension == 2) {
    _compute_elt_normal = &_compute_face_normal;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Invalid dimension (expected 1 or 2, got %d)\n", dimension);
  }


  PDM_malloc(*out_selected_vtx_normal, n_part, double *);

  int    **send_stride = NULL;
  double **send_normal = NULL;
  PDM_malloc(send_stride, n_part, int    *);
  PDM_malloc(send_normal, n_part, double *);

  int **all_vtx_to_selected_vtx = NULL;
  PDM_malloc(all_vtx_to_selected_vtx, n_part, int *);

  for (int i_part = 0; i_part < n_part; i_part++) {

    (*out_selected_vtx_normal)[i_part] = PDM_array_zeros_double(n_selected_vtx[i_part] * 3);

    double *vtx_normal = (*out_selected_vtx_normal)[i_part];

    // Get indirection from all vertices to selected vertices
    int n_elt = 0;
    int n_vtx = 0;
    // This is a bit complicated for just getting the total number of vtx => just pass it as an input?
    if (selected_elt == NULL) {
      // Account for all entities1
      n_elt = n_selected_elt[i_part];
      for (int i_elt = 0; i_elt < n_elt; i_elt++) {
        int elt_vtx_n = elt_vtx_idx[i_part][i_elt+1]-elt_vtx_idx[i_part][i_elt];
        CHECK_ELT_SIZE(dimension, i_elt, elt_vtx_n)
      }
      for (int idx_vtx = 0; idx_vtx < elt_vtx_idx[i_part][n_selected_elt[i_part]]; idx_vtx++) {
        n_vtx = PDM_MAX(n_vtx, elt_vtx[i_part][idx_vtx]);
      }
    }
    else {
      // Only account for subset of entities1
      for (int idx_elt = 0; idx_elt < n_selected_elt[i_part]; idx_elt++) {
        int i_elt = selected_elt[i_part][idx_elt] - 1;
        n_elt = PDM_MAX(n_elt, selected_elt[i_part][idx_elt]);
        int elt_vtx_n = elt_vtx_idx[i_part][i_elt+1]-elt_vtx_idx[i_part][i_elt];
        CHECK_ELT_SIZE(dimension, i_elt, elt_vtx_n)
        for (int idx_vtx = elt_vtx_idx[i_part][i_elt]; idx_vtx < elt_vtx_idx[i_part][i_elt+1]; idx_vtx++) {
          n_vtx = PDM_MAX(n_vtx, elt_vtx[i_part][idx_vtx]);
        }
      }
    }

    if (selected_vtx == NULL) {
      if (n_selected_vtx[i_part] < n_vtx) {
        PDM_error(
          __FILE__, __LINE__, 0,
          "Element connectivity has vtx id > n_selected_vtx (%d > %d) for partition n°%d\n",
          n_vtx,
          n_selected_vtx[i_part],
          i_part
        );
      }
      n_vtx = PDM_MAX(n_vtx, n_selected_vtx[i_part]);
    }
    else {
      for (int idx_vtx = 0; idx_vtx < n_selected_vtx[i_part]; idx_vtx++) {
        n_vtx = PDM_MAX(n_vtx, selected_vtx[i_part][idx_vtx]);
      }
    }

    int *graph_vtx = NULL;
    int n_graph_vtx = PDM_part_comm_graph_entity_graph_get(pcg_vtx,
                                                           i_part,
                                                          &graph_vtx,
                                                           PDM_OWNERSHIP_BAD_VALUE);

    for (int idx_vtx = 0; idx_vtx < n_graph_vtx; idx_vtx++) {
      n_vtx = PDM_MAX(n_vtx, graph_vtx[4*idx_vtx]);
    }

    all_vtx_to_selected_vtx[i_part] = PDM_array_const_int(n_vtx, -1);
    for (int idx_vtx = 0; idx_vtx < n_selected_vtx[i_part]; idx_vtx++) {
      int i_vtx = (selected_vtx == NULL) ? idx_vtx : selected_vtx[i_part][idx_vtx] - 1;
      all_vtx_to_selected_vtx[i_part][i_vtx] = idx_vtx;
    }


    int *elt_is_ghost = PDM_array_zeros_int(n_elt);
    const int *graph_elt_owner = PDM_part_comm_graph_owner_get(pcg_elt, i_part);

    int *graph_elt = NULL;
    int n_graph_elt = PDM_part_comm_graph_entity_graph_get(pcg_elt,
                                                           i_part,
                                                          &graph_elt,
                                                           PDM_OWNERSHIP_BAD_VALUE);

    for (int idx_elt = 0; idx_elt < n_graph_elt; idx_elt++) {
      if (!graph_elt_owner[idx_elt]) {
        int i_elt = graph_elt[4*idx_elt] - 1;
        if (i_elt < n_elt) {
          elt_is_ghost[i_elt] = 1;
        }
      }
    }


    // Compute element normals and add contribution to incident vertices
    // At this point, the normal vectors are scaled by the elements measure
    for (int idx_elt = 0; idx_elt < n_selected_elt[i_part]; idx_elt++) {

      int i_elt = (selected_elt == NULL) ? idx_elt : selected_elt[i_part][idx_elt] - 1;

      // Ignore if ghost element
      if (elt_is_ghost[i_elt]) {
        continue;
      }

      double elt_normal[3];
      _compute_elt_normal(elt_vtx_idx[i_part][i_elt+1] - elt_vtx_idx[i_part][i_elt],
                          &elt_vtx[i_part][elt_vtx_idx[i_part][i_elt]],
                          vtx_coord[i_part],
                          elt_normal);

      for (int idx_vtx = elt_vtx_idx[i_part][i_elt]; idx_vtx < elt_vtx_idx[i_part][i_elt+1]; idx_vtx++) {
        int i_vtx = elt_vtx[i_part][idx_vtx] - 1;
        int i_selected_vtx = all_vtx_to_selected_vtx[i_part][i_vtx];
        if (i_selected_vtx < 0) {
          // PDM_error(__FILE__, __LINE__, 0, "Vertex %d of part %d is not selected (rank %d).\n", i_vtx, i_part, i_rank); //juste ignore it ?
        }
        else {

          for (int i = 0; i < 3; i++) {
            vtx_normal[3*i_selected_vtx+i] += elt_normal[i]; // TODO: scale by fraction of dual measure?
          }
        }
      }
    } // End loop on elements
    PDM_free(elt_is_ghost);

    // Prepare send buffer for inter-partition synchronization
    PDM_malloc(send_stride[i_part], n_graph_vtx,     int   );
    PDM_malloc(send_normal[i_part], n_graph_vtx * 3, double);
    int idx_write = 0;
    for (int idx_vtx = 0; idx_vtx < n_graph_vtx; idx_vtx++) {
      int i_vtx = graph_vtx[4*idx_vtx] - 1;
      int i_selected_vtx = all_vtx_to_selected_vtx[i_part][i_vtx];
      if (i_selected_vtx < 0) {
        send_stride[i_part][idx_vtx] = 0;
      }
      else {
        send_stride[i_part][idx_vtx] = 1;
        for (int i = 0; i < 3; i++) {
          send_normal[i_part][idx_write++] = vtx_normal[3*i_selected_vtx+i];
        }
      }
    }

  } // End loop on parts


  /* Synchronize part boundaries */
  int    **recv_stride = NULL;
  double **recv_normal = NULL;
  PDM_part_comm_graph_exch(pcg_vtx,
                           sizeof(double) * 3,
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           send_stride,
                (void  **) send_normal,
                           &recv_stride,
                (void ***) &recv_normal);


  for (int i_part = 0; i_part < n_part; i_part++) {

    double *vtx_normal = (*out_selected_vtx_normal)[i_part];

    int *graph_vtx = NULL;
    int graph_vtx_n = PDM_part_comm_graph_entity_graph_get(pcg_vtx,
                                                           i_part,
                                                          &graph_vtx,
                                                           PDM_OWNERSHIP_BAD_VALUE);

    // Add remote contributions
    int idx_read = 0;
    for (int idx_vtx = 0; idx_vtx < graph_vtx_n; idx_vtx++) {
      int i_vtx = graph_vtx[4*idx_vtx] - 1;
      int i_selected_vtx = all_vtx_to_selected_vtx[i_part][i_vtx];
      if (i_selected_vtx >= 0) {
        if (recv_stride[i_part][idx_vtx] != 1) {
          PDM_error(__FILE__, __LINE__, 0, "Inconsistent group_vtx between ranks %d and %d, part %d and %d\n",
                     i_rank, graph_vtx[4*idx_vtx+1],
                     i_part, graph_vtx[4*idx_vtx+2]);
        }
        for (int i = 0; i < 3; i++) {
          vtx_normal[3*i_selected_vtx+i] += recv_normal[i_part][idx_read++];
        }
      }
    }

    // Normalize (or don't if you want to be scaled by dual measure)
    for (int i_vtx = 0; i_vtx < n_selected_vtx[i_part]; i_vtx++) {
      double magnitude = PDM_MODULE(&vtx_normal[3*i_vtx]);
      if (magnitude <= 0) { // can happen if vertice is not referenced by edge
        // PDM_error(__FILE__, __LINE__, 0, "Singular vertex %d in part %d of rank %d\n", selected_vtx[i_part][i_vtx], i_part, i_rank);
      }
      else {
        double inv_magnitude = 1./magnitude;
        for (int i = 0; i < 3; i++) {
          vtx_normal[3*i_vtx+i] *= inv_magnitude;
        }
      }
    }

    PDM_free(all_vtx_to_selected_vtx[i_part]);

    PDM_free(send_stride[i_part]);
    PDM_free(send_normal[i_part]);
    PDM_free(recv_stride[i_part]);
    PDM_free(recv_normal[i_part]);
  } // End loop on parts


  // Free memory
  PDM_free(all_vtx_to_selected_vtx);

  PDM_free(send_stride);
  PDM_free(send_normal);
  PDM_free(recv_stride);
  PDM_free(recv_normal);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
