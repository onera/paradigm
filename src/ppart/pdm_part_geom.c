
/*----------------------------------------------------------------------------
 *  Standar headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_error.h"
#include "pdm_part_geom.h"
#include "pdm_hilbert.h"
#include "pdm_sort.h"
#include "pdm_geom_elem.h"
#include "pdm_binary_search.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_distrib.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dgeom_elem.h"
#include "pdm_array.h"

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
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Compute cell center of elements
 *
 * \param [in]  comm         MPI Communicator
 * \param [in]  dn_cell       Number of cells in the current process
 * \param [in]  dn_face       Number of faces in the current process
 * \param [in]  dn_vtx        Number of vertices in the current process
 * \param [in]  dcell_face_idx Index of cell_face
 * \param [in]  dcell_face    cell face connectivity in the current process
 * \param [in]  dface_vtx_idx  Index of face_vtx
 * \param [in]  dface_vtx     face vertex connectivity in the current process
 * \param [in]  distrib_face    face distribution
 * \param [in]  dvtx_coord    coordinates of vertices
 * \param [in]  distrib_vtx     Vertex distribution
 *
 * \param [out] cell_center   Cell centers
 *
 */

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

  /*PDM_log_trace_array_long (distrib_face, n_rank+1, "distrib_face : ");
    PDM_log_trace_array_long (distrib_vtx,  n_rank+1, "distrib_vtx  : ");*/

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

  /*PDM_log_trace_connectivity_long(dface_vtx_idx, dface_vtx, dn_face, "dface_vtx : ");
  PDM_log_trace_connectivity_int (pface_vtx_idx, pface_vtx, dn_face, "pface_vtx : ");
  for (int i = 0; i < pn_vtx; i++) {
    log_trace("vtx %d -> "PDM_FMT_G_NUM"\n", i+1, pvtx_ln_to_gn[i]);
  }*/

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
/*=============================================================================
 * Public function definitions
 *============================================================================*/

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

/**
 *
 * \brief Perform geomtric patitionning
 *
 * \param [in]   method         Geometric method
 * \param [in]   n_part          Number of partition to build on this process
 * \param [in]   comm           Communicator
 * \param [in]   dn_cell         Number of distributed cells
 * \param [in]   dn_face         Number of distributed faces
 * \param [in]   dn_vtx          Number of distributed vertices
 * \param [in]   dcell_face_idx   Distributed cell face connectivity index or NULL
 *                              (size : dn_cell + 1, numbering : 0 to n-1)
 * \param [in]   dcell_face      Distributed cell face connectivity or NULL
 *                              (size : dface_vtx_idx[dn_cell], numbering : 1 to n)
 * \param [in]   dcell_weight    Cell weight (size : n_cell) or NULL
 * \param [in]   dface_vtx_idx    Distributed face to vertex connectivity index
 *                              (size : dn_face + 1, numbering : 0 to n-1)
 * \param [in]   dface_vtx       Distributed face to vertex connectivity
 *                              (size : dface_vtx_idx[dn_face], numbering : 1 to n)
 * \param [in]   dvtx_coord      Distributed vertex coordinates
 *                              (size : 3*dn_vtx)
 * \param [inout]   dcell_part      Distributed cell partitioning
 *                              (size = dn_cell)
 *
 */

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


#ifdef __cplusplus
}
#endif /* __cplusplus */
