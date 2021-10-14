
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


static void
_prepare_connectivity
(
const PDM_MPI_Comm    comm,
const PDM_g_num_t     dn_cell,
const int            *dcell_face_idx,
const PDM_g_num_t    *dcell_face,
const int            *dface_vtx_idx,
const PDM_g_num_t    *dface_vtx,
const PDM_g_num_t    *dface_proc,
int                **face_vtx,
int                **face_vtx_idx,
int                  *sizeface_vtx_idx,
const PDM_g_num_t    *dvtx_proc,
const double         *dvtx_coord,
      double       **lvtx_coord
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /* Offset Distribution to work with block_to_part */
  PDM_g_num_t* dface_proc_loc = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * (n_rank + 1));
  PDM_g_num_t* dvtx_proc_loc  = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * (n_rank + 1));
  for (int i = 0; i < n_rank+1; i++) {
    dface_proc_loc[i] = dface_proc[i]-1;
    dvtx_proc_loc[i]  = dvtx_proc[i]-1;
  }


  /* Verbose */
  if(0 == 1){
    PDM_printf("dface_proc : \n");
    for (int i = 0; i < n_rank+1; i++) {
      PDM_printf("[%i] -> %i \n ",i,  dface_proc[i]);
      PDM_printf("[%i] -> %i \n ",i,  dface_proc_loc[i]);
    }
    for (int i = 0; i < dn_cell; i++) {
      int nFac = dcell_face_idx[i+1] - dcell_face_idx[i];
      printf("nFac[%i] :  %i \n ", i, nFac);
      for (int i_vtx = dcell_face_idx[i]; i_vtx < dcell_face_idx[i+1]; i_vtx++) {
        printf("dcell_face[%i] :  "PDM_FMT_G_NUM" \n ", i_vtx, dcell_face[i_vtx]);
      }
    }
  }


  /* -------------------------------------------------------------- */
  PDM_block_to_part_t *ptb = PDM_block_to_part_create (dface_proc_loc,
                                                       &dcell_face,
                                                       &dcell_face_idx[dn_cell],
                                                        1,
                                                        comm);

  /* Echange sur les nombre de sommet de chaque face */
  int dn_face = dface_proc[i_rank+1] - dface_proc[i_rank];

  int* n_vtx_face_loc = (int *) malloc (sizeof(int) * dn_face              );
  int* n_vtx_face     = (int *) malloc (sizeof(int) * dcell_face_idx[dn_cell]);

  for (int i = 0; i < dn_face; i++) {
    n_vtx_face_loc[i] = dface_vtx_idx[i+1] - dface_vtx_idx[i];
  }

  int stride_one = 1;

  PDM_block_to_part_exch (ptb,
                          sizeof(PDM_l_num_t),
                          PDM_STRIDE_CST,
                          &stride_one,
                 (void *) n_vtx_face_loc,
                          NULL,
                (void **) &n_vtx_face);

  /* Verbose */
  if(0 == 1){
    for (int i = 0; i < dcell_face_idx[dn_cell]; i++) {
      printf("n_vtx_face :  %i \n ", n_vtx_face[i]);
    }
  }

  int ndcell_face_tot = 0;
  for (int i = 0; i < dcell_face_idx[dn_cell]; i++) {
    ndcell_face_tot += n_vtx_face[i];
  }
  // printf("ndcell_face_tot : %i \n", ndcell_face_tot);

  /* Alloc field */
  PDM_g_num_t* lface_vtx    = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * ndcell_face_tot);

  PDM_block_to_part_exch (           ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     n_vtx_face_loc,
                          (void *)   dface_vtx,
                                    &n_vtx_face,
                          (void **) &lface_vtx);

  /* Verbose */
  if(1 == 1){
    for (int i = 0; i < dcell_face_idx[dn_cell]; i++) {
      printf("n_vtx_face :  %i \n ", n_vtx_face[i]);
    }
    for (int i = 0; i < ndcell_face_tot; i++) {
      printf("lface_vtx[%i] :  "PDM_FMT_G_NUM" \n ",i, lface_vtx[i]);
    }
  }

  PDM_block_to_part_free(ptb);
  free (dface_proc_loc);
  /* -------------------------------------------------------------- */


  /* -------------------------------------------------------------- */
  ptb = PDM_block_to_part_create (dvtx_proc_loc,
                                  (const PDM_g_num_t **) &lface_vtx,
                                  &ndcell_face_tot,
                                  1,
                                  comm);

  int stride = 3;


  // int dn_vtx = dvtx_proc[i_rank+1] - dvtx_proc[i_rank];
  // for (int i = 0; i < dn_vtx; i++) {
  //   printf("dvtx_coord[%i] :  %12.5e \n ",i, dvtx_coord[i]);
  // }

  *lvtx_coord = (double *) malloc (sizeof(double) * 3 * ndcell_face_tot);
  PDM_block_to_part_exch (ptb,
                          sizeof(double),
                          PDM_STRIDE_CST,
                          &stride,
                          (void *) dvtx_coord,
                          NULL,
                          (void **) lvtx_coord);

  /* Verbose */
  if(0 == 1){
    for (int i = 0; i < ndcell_face_tot; i++) {
      printf("lvtx_coord[%i] :  %12.5e \n ",i, (*lvtx_coord)[i]);
    }
  }

  PDM_block_to_part_free(ptb);
  free (dvtx_proc_loc);
  /* -------------------------------------------------------------- */

  *sizeface_vtx_idx = dcell_face_idx[dn_cell] + 1 ;

  *face_vtx    = (int *) malloc( sizeof(int *) *   ndcell_face_tot              );
  *face_vtx_idx = (int *) malloc( sizeof(int *) * ( dcell_face_idx[dn_cell] + 1 ) );

  for (int i = 0; i < ndcell_face_tot; i++) {
    (*face_vtx)[i] = i+1;
  }

  (*face_vtx_idx)[0] = 0;
  for (int i = 0; i < dcell_face_idx[dn_cell]; i++) {
    (*face_vtx_idx)[i+1] = (*face_vtx_idx)[i] + n_vtx_face[i];
  }

  /* Free memory */
  free(n_vtx_face);
  free(n_vtx_face_loc);
  free(lface_vtx);

}





static
void
_compute_center_from_descending_connectivity
(
  const int            *dentity1_entity2_idx,
  const PDM_g_num_t    *dentity1_entity2,
  const int             dn_entity1,
  const PDM_g_num_t    *dentity2_distrib,
  double               *dentity1_coord,
  double               *dentity2_coord,
  PDM_MPI_Comm          comm
)
{
  int *dentity1_entity2_sgn = malloc(dentity1_entity2_idx[dn_entity1] * sizeof(int));
  PDM_g_num_t *dentity1_entity2_abs = malloc(dentity1_entity2_idx[dn_entity1] * sizeof(PDM_g_num_t));
  for(int i = 0; i < dentity1_entity2_idx[dn_entity1]; ++i) {
    dentity1_entity2_sgn[i] = PDM_SIGN(dentity1_entity2[i]);
    dentity1_entity2_abs[i] = PDM_ABS (dentity1_entity2[i]);
  }
  PDM_block_to_part_t *btp_entity1_coord = PDM_block_to_part_create (dentity2_distrib,
                                              (const PDM_g_num_t **) &dentity1_entity2_abs,
                                                                     &dentity1_entity2_idx[dn_entity1],
                                                                     1,
                                                                     comm);
  /*for(int i = 0; i < dentity1_entity2_idx[dn_entity1]; ++i) {
    dentity1_entity2[i]     = dentity1_entity2[i]*dentity1_entity2_sgn[i];
    }*/
  free(dentity1_entity2_sgn);
  free(dentity1_entity2_abs);

  int strid_one = 1;
  double **tmp_entity1_entity2_coord;
  PDM_block_to_part_exch2 (btp_entity1_coord,
                           3 * sizeof(double),
                           PDM_STRIDE_CST,
                           &strid_one,
                  (void *) dentity2_coord,
                           NULL,
                (void ***) &tmp_entity1_entity2_coord);
  double *dentity1_entity2_coord = tmp_entity1_entity2_coord[0];
  free(tmp_entity1_entity2_coord);
  PDM_block_to_part_free(btp_entity1_coord);

  for(int i_entity1 = 0; i_entity1 < dn_entity1; ++i_entity1) {
    dentity1_coord[3*i_entity1  ] = 0.;
    dentity1_coord[3*i_entity1+1] = 0.;
    dentity1_coord[3*i_entity1+2] = 0.;
    int n_entity2_per_entity1 = dentity1_entity2_idx[i_entity1+1] - dentity1_entity2_idx[i_entity1];
    double inv = 1./n_entity2_per_entity1;
    for(int idx_entity2 = dentity1_entity2_idx[i_entity1]; idx_entity2 < dentity1_entity2_idx[i_entity1+1]; ++idx_entity2) {
      dentity1_coord[3*i_entity1  ] += dentity1_entity2_coord[3*idx_entity2  ];
      dentity1_coord[3*i_entity1+1] += dentity1_entity2_coord[3*idx_entity2+1];
      dentity1_coord[3*i_entity1+2] += dentity1_entity2_coord[3*idx_entity2+2];
    }
    dentity1_coord[3*i_entity1  ] = dentity1_coord[3*i_entity1  ]*inv;
    dentity1_coord[3*i_entity1+1] = dentity1_coord[3*i_entity1+1]*inv;
    dentity1_coord[3*i_entity1+2] = dentity1_coord[3*i_entity1+2]*inv;

  }
  free(dentity1_entity2_coord);


}






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
 * \param [in]  dface_proc    face distribution
 * \param [in]  dvtx_coord    coordinates of vertices
 * \param [in]  dvtx_proc     Vertex distribution
 *
 * \param [out] cell_center   Cell centers
 *
 */

void
PDM_dcompute_cell_center_old
(
  const PDM_MPI_Comm  comm,
  const int           dn_cell,
  const int          *dcell_face_idx,
  const PDM_g_num_t  *dcell_face,
  const int          *dface_vtx_idx,
  const PDM_g_num_t  *dface_vtx,
  const PDM_g_num_t  *dface_proc,
  const double       *dvtx_coord,
  const PDM_g_num_t  *dvtx_proc,
  double             *cell_center
)
{

  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /* Local connectivity (Allocate in _prepareconnectivity */
  int sizeface_vtx_idx = 0;

  int *face_vtx;
  int *face_vtx_idx;

  /* Local coordinates */
  double *lvtx_coord;

  /*
   * Get face connectivities from other process
   */
  _prepare_connectivity( comm,
                         dn_cell,
                         dcell_face_idx,
                         dcell_face,
                         dface_vtx_idx,
                         dface_vtx,
                         dface_proc,
                        &face_vtx,
                        &face_vtx_idx,
                        &sizeface_vtx_idx,
                         dvtx_proc,
                         dvtx_coord,
                        &lvtx_coord);


  /* Verbose */
  if(0 == 1 ){
    PDM_printf("face_vtx : \n");
    for (int i1 = 0; i1 < sizeface_vtx_idx; i1++) {
      for (int j = face_vtx_idx[i1]; j < face_vtx_idx[i1+1]; j++)
        PDM_printf(" %i", face_vtx[j]);
      PDM_printf("\n");
    }
  }
  free (face_vtx);

  /*
   * Compute cell centers
   */
  for (int i_cell = 0; i_cell < dn_cell; i_cell++) {

    cell_center[3*i_cell    ] = 0.;
    cell_center[3*i_cell + 1] = 0.;
    cell_center[3*i_cell + 2] = 0.;

    int n_vtx_on_cell = 0;

    /* Loop on all faces */
    for (int i_face = dcell_face_idx[i_cell  ]; i_face < dcell_face_idx[i_cell+1]; i_face++) {
      n_vtx_on_cell += face_vtx_idx[i_face+1] - face_vtx_idx[i_face  ];
      /* Loop on all Vtx of current faces */
      for (int i_vtx = face_vtx_idx[i_face  ]; i_vtx < face_vtx_idx[i_face+1]; i_vtx++){

        cell_center[3*i_cell    ] += lvtx_coord[3*i_vtx    ];
        cell_center[3*i_cell + 1] += lvtx_coord[3*i_vtx + 1];
        cell_center[3*i_cell + 2] += lvtx_coord[3*i_vtx + 2];

      }
    }

    /* Ponderate */
    cell_center[3*i_cell    ] = cell_center[3*i_cell    ]/n_vtx_on_cell;
    cell_center[3*i_cell + 1] = cell_center[3*i_cell + 1]/n_vtx_on_cell;
    cell_center[3*i_cell + 2] = cell_center[3*i_cell + 2]/n_vtx_on_cell;
  }
  free (lvtx_coord);
  free (face_vtx_idx);

  /* Verbose */
  if(0 == 1){
    for (int i_cell = 0; i_cell < dn_cell; i_cell++) {
      PDM_printf(" --------- %i / %i \n ", i_cell);
      PDM_printf(" %12.5e %12.5e %12.5e \n", cell_center[3*i_cell    ], cell_center[3*i_cell + 1], cell_center[3*i_cell + 2]);
    }
  }

}


void
PDM_dcompute_cell_center
(
  const PDM_MPI_Comm  comm,
  const int           dn_cell,
  const int          *dcell_face_idx,
  const PDM_g_num_t  *dcell_face,
  const int          *dface_vtx_idx,
  const PDM_g_num_t  *dface_vtx,
  const PDM_g_num_t  *dface_proc,
  const double       *dvtx_coord,
  const PDM_g_num_t  *dvtx_proc,
  double             *cell_center
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  PDM_g_num_t *distrib_face = malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  PDM_g_num_t *distrib_vtx  = malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  for (int i = 0; i <= n_rank; i++) {
    distrib_face[i] = dface_proc[i] - 1;
    distrib_vtx[i]  = dvtx_proc[i]  - 1;
  }

  /*PDM_log_trace_array_long (distrib_face, n_rank+1, "distrib_face : ");
    PDM_log_trace_array_long (distrib_vtx,  n_rank+1, "distrib_vtx  : ");*/

  int dn_face = (int) (dface_proc[i_rank+1] - dface_proc[i_rank]);

  PDM_g_num_t *dface_ln_to_gn = malloc (sizeof(PDM_g_num_t) * dn_face);
  for (int i = 0; i < dn_face; i++) {
    dface_ln_to_gn[i] = dface_proc[i_rank] + i;
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
  free (dface_ln_to_gn);

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
  free(tmp_pvtx_coord);
  free (pvtx_ln_to_gn);
  free (distrib_vtx);


  /* Compute face centers */
  double *dface_center = malloc (sizeof(double) * dn_face * 3);
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
  free (pvtx_coord);
  free (pface_vtx_idx);
  free (pface_vtx);

  /* Compute cell centers */
  _compute_center_from_descending_connectivity (dcell_face_idx,
                                                dcell_face,
                                                dn_cell,
                                                distrib_face,
                                                cell_center,
                                                dface_center,
                                                comm);
  free (dface_center);
  free (distrib_face);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/


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
 const PDM_g_num_t  *dface_proc,
 const double       *dvtx_coord,
 const PDM_g_num_t  *dvtx_proc,
 int                *dcell_part
)
{

  assert (method == PDM_PART_GEOM_HILBERT);

  const int dim = 3;

  double *barycenter_coords = (double *) malloc (dn_cell * 3 * sizeof(double ));

  PDM_hilbert_code_t *hilbert_codes     = (PDM_hilbert_code_t *) malloc (dn_cell * sizeof(PDM_hilbert_code_t));
  PDM_hilbert_code_t *tmp_hilbert_codes = (PDM_hilbert_code_t *) malloc (dn_cell * sizeof(PDM_hilbert_code_t));

  /*
   * cell center computation
   */
  PDM_dcompute_cell_center (comm,
                            dn_cell,
                            dcell_face_idx,
                            dcell_face,
                            dface_vtx_idx,
                            dface_vtx,
                            dface_proc,
                            dvtx_coord,
                            dvtx_proc,
                            barycenter_coords);


  /** TRAITEMENT HILBERT FVM **/

	/** Initialisation **/

  double extents[2*dim]; /** DIM x 2**/

	/** Get EXTENTS **/
  PDM_hilbert_get_coord_extents_par(dim, dn_cell, barycenter_coords, extents, comm);

  // int i_rank;
  // PDM_MPI_Comm_rank  (comm, &i_rank);
  // char filename[999];
  // sprintf(filename, "barycenter_coords_%2.2d.vtk", i_rank);
  // PDM_vtk_write_point_cloud(filename,
  //                           dn_cell,
  //                           barycenter_coords,
  //                           NULL,
  //                           NULL);
  // PDM_log_trace_array_double(extents, 6, "extents :: ");
  // PDM_log_trace_array_double(barycenter_coords, 3*dn_cell, "barycenter_coords :: ");

	/** Hilbert Coordinates Computation **/
  PDM_hilbert_encode_coords(dim, PDM_HILBERT_CS, extents, dn_cell, barycenter_coords, hilbert_codes);

  for (int i = 0; i < dn_cell; ++i) {
    tmp_hilbert_codes [i] = hilbert_codes [i];
	}

  ///** Calcul des index des codes Hilbert **/

  int * hilbert_order = (int * ) malloc (dn_cell * sizeof(int));

  for (int i = 0; i < dn_cell; ++i) {
    hilbert_order [i] = i;
  }

  assert (sizeof(double) == sizeof(PDM_hilbert_code_t));
  PDM_sort_double (tmp_hilbert_codes, hilbert_order, dn_cell);

  free(tmp_hilbert_codes);

  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  PDM_hilbert_code_t *hilbert_codes_idx = (PDM_hilbert_code_t *) malloc ((n_rank+1)*n_part * sizeof(PDM_hilbert_code_t));

  int * weight = (int *) malloc (dn_cell * sizeof(int));
  if (dcell_weight != NULL) {
    for(int i = 0; i < dn_cell; ++i) {
		  weight [i] = dcell_weight [i];
    }
  }
  else {
    for(int i = 0; i < dn_cell; ++i) {
		  weight [i] = 1;
    }
  }

  int n_total_part;
  PDM_MPI_Allreduce ((void *) &n_part, &n_total_part, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);

  PDM_hilbert_build_rank_index (dim,
                                n_total_part,
                                dn_cell,
                                hilbert_codes,
                                weight,
                                hilbert_order,
                                hilbert_codes_idx,
                                comm);

  free(weight);

  /** Remplissage de cell_parts -> en fct des codes Hilbert **/

  for(int i = 0; i < dn_cell; ++i) {
    size_t quantile = PDM_hilbert_quantile_search(n_total_part,
                                                hilbert_codes[i],
                                                hilbert_codes_idx);
    dcell_part [i] = (int) quantile;

  }

  free(barycenter_coords);
  free(hilbert_codes_idx);
  free(hilbert_order);
  free(hilbert_codes);

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

  PDM_hilbert_code_t *hilbert_codes     = (PDM_hilbert_code_t *) malloc (dn_vtx * sizeof(PDM_hilbert_code_t));
  PDM_hilbert_code_t *tmp_hilbert_codes = (PDM_hilbert_code_t *) malloc (dn_vtx * sizeof(PDM_hilbert_code_t));

  /** Initialisation **/
  double extents[2*dim]; /** DIM x 2**/

  /** Get EXTENTS **/
  PDM_hilbert_get_coord_extents_par(dim, dn_vtx, dcoords, extents, comm);

  /** Hilbert Coordinates Computation **/
  PDM_hilbert_encode_coords(dim, PDM_HILBERT_CS, extents, dn_vtx, dcoords, hilbert_codes);

  for (int i = 0; i < dn_vtx; ++i) {
    tmp_hilbert_codes [i] = hilbert_codes [i];
  }

  ///** Calcul des index des codes Hilbert **/
  int * hilbert_order = (int * ) malloc (dn_vtx * sizeof(int));

  for (int i = 0; i < dn_vtx; ++i) {
    hilbert_order [i] = i;
  }

  assert (sizeof(double) == sizeof(PDM_hilbert_code_t));
  PDM_sort_double (tmp_hilbert_codes, hilbert_order, dn_vtx);

  free(tmp_hilbert_codes);

  PDM_hilbert_code_t *hilbert_codes_idx = (PDM_hilbert_code_t *) malloc ((n_rank+1) * sizeof(PDM_hilbert_code_t));

  int *weight = (int *) malloc (dn_vtx * sizeof(int));
  for(int i = 0; i < dn_vtx; ++i) {
    weight [i] = 1;
  }

  PDM_hilbert_build_rank_index (dim,
                                n_rank,
                                dn_vtx,
                                hilbert_codes,
                                weight,
                                hilbert_order,
                                hilbert_codes_idx,
                                comm);


  free(weight);

  /** Remplissage de cell_parts -> en fct des codes Hilbert **/
  for(int i = 0; i < dn_vtx; ++i) {
    size_t quantile = PDM_hilbert_quantile_search(n_rank,
                                                  hilbert_codes[i],
                                                  hilbert_codes_idx);
    ln_to_gn [i] = (int) (quantile + 1); // Because ln_to_gn of part_to_block begin at 1
  }

  // part_to_block avec ln_to_gn 1 2 3 4 .... pdm_assembly_partition
  // Puis on échange les hilbert_codes, qu'on retrie localement

  PDM_g_num_t* distrib_rank = (PDM_g_num_t *) malloc( (n_rank+1) * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_rank+1; ++i) {
    distrib_rank[i] = i;
  }

  /*
   * Each proc get all the entities affected to its partitions
   */
  PDM_part_to_block_t *ptb =
   PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                             PDM_PART_TO_BLOCK_POST_NOTHING,
                             1.,
                             &ln_to_gn,
                             distrib_rank,
                             &dn_vtx,
                             1,
                             comm);
  free(distrib_rank);

  const int n_vtx_block = PDM_part_to_block_n_elt_block_get (ptb);

  // log_trace("n_vtx_block = %i | dn_vtx = %i \n", n_vtx_block, dn_vtx);

  double *blk_hilbert_codes = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(double),
                          PDM_STRIDE_CST,
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
                          PDM_STRIDE_CST,
                          1,
                          NULL,
                (void **) &ln_to_gn,
                          NULL,
                (void **) &blk_ln_to_gn);

  free(hilbert_codes_idx);
  free(hilbert_codes);
  PDM_part_to_block_free(ptb);

  /* Reorder locally */
  assert (sizeof(double) == sizeof(PDM_hilbert_code_t));
  hilbert_order = (int * ) realloc(hilbert_order,  n_vtx_block * sizeof(int));
  for (int i = 0; i < n_vtx_block; ++i) {
    hilbert_order [i] = i;
  }
  PDM_sort_double (blk_hilbert_codes, hilbert_order, n_vtx_block);
  free(blk_hilbert_codes);


  /* Apply order to blk_ln_to_gn */
  PDM_g_num_t* sorted_blk_ln_to_gn = malloc( n_vtx_block * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_vtx_block; ++i) {
    sorted_blk_ln_to_gn[i] = blk_ln_to_gn[hilbert_order[i]];
  }
  free(blk_ln_to_gn);
  free(hilbert_order);

  PDM_g_num_t* distrib_blk_vtx = PDM_compute_entity_distribution(comm, n_vtx_block);
  PDM_block_to_part_t* btp = PDM_block_to_part_create(distrib_blk_vtx,
                              (const PDM_g_num_t **)  &ln_to_gn,
                                                      &dn_vtx,
                                                      1,
                                                      comm);

  int stride_one = 1;
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST,
                         &stride_one,
              (void *)   sorted_blk_ln_to_gn,
                         NULL,
              (void **) &ln_to_gn);
  PDM_block_to_part_free(btp);
  free(sorted_blk_ln_to_gn);
  free(distrib_blk_vtx);
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

  // PDM_hilbert_code_t *hilbert_codes     = (PDM_hilbert_code_t *) malloc (dn_length * sizeof(PDM_hilbert_code_t));
  PDM_hilbert_code_t *tmp_hilbert_codes = (PDM_hilbert_code_t *) malloc (dn_length * sizeof(PDM_hilbert_code_t));

  for (int i = 0; i < dn_length; ++i) {
    tmp_hilbert_codes [i] = length [i];
  }

  ///** Calcul des index des codes Hilbert **/
  int * hilbert_order = (int * ) malloc (dn_length * sizeof(int));

  for (int i = 0; i < dn_length; ++i) {
    hilbert_order [i] = i;
  }

  assert (sizeof(double) == sizeof(PDM_hilbert_code_t));
  PDM_sort_double (tmp_hilbert_codes, hilbert_order, dn_length);

  free(tmp_hilbert_codes);

  PDM_hilbert_code_t *hilbert_codes_idx = (PDM_hilbert_code_t *) malloc ((n_rank+1) * sizeof(PDM_hilbert_code_t));

  int *weight = (int *) malloc (dn_length * sizeof(int));
  for(int i = 0; i < dn_length; ++i) {
    weight [i] = 1;
  }

  PDM_hilbert_build_rank_index (dim,
                                n_rank,
                                dn_length,
                                length,
                                weight,
                                hilbert_order,
                                hilbert_codes_idx,
                                comm);


  free(weight);

  /** Remplissage de cell_parts -> en fct des codes Hilbert **/
  for(int i = 0; i < dn_length; ++i) {
    size_t quantile = PDM_hilbert_quantile_search(n_rank,
                                                  length[i],
                                                  hilbert_codes_idx);
    ln_to_gn [i] = (int) (quantile + 1); // Because ln_to_gn of part_to_block begin at 1
  }

  // part_to_block avec ln_to_gn 1 2 3 4 .... pdm_assembly_partition
  // Puis on échange les hilbert_codes, qu'on retrie localement

  PDM_g_num_t* distrib_rank = (PDM_g_num_t *) malloc( (n_rank+1) * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_rank+1; ++i) {
    distrib_rank[i] = i;
  }

  /*
   * Each proc get all the entities affected to its partitions
   */
  PDM_part_to_block_t *ptb =
   PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                             PDM_PART_TO_BLOCK_POST_NOTHING,
                             1.,
                             &ln_to_gn,
                             distrib_rank,
                             &dn_length,
                             1,
                             comm);
  free(distrib_rank);

  const int n_vtx_block = PDM_part_to_block_n_elt_block_get (ptb);

  // log_trace("n_vtx_block = %i | dn_length = %i \n", n_vtx_block, dn_length);

  double *blk_hilbert_codes = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(double),
                          PDM_STRIDE_CST,
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
                          PDM_STRIDE_CST,
                          1,
                          NULL,
                (void **) &ln_to_gn,
                          NULL,
                (void **) &blk_ln_to_gn);

  free(hilbert_codes_idx);
  PDM_part_to_block_free(ptb);

  /* Reorder locally */
  assert (sizeof(double) == sizeof(PDM_hilbert_code_t));
  hilbert_order = (int * ) realloc(hilbert_order,  n_vtx_block * sizeof(int));
  for (int i = 0; i < n_vtx_block; ++i) {
    hilbert_order [i] = i;
  }
  PDM_sort_double (blk_hilbert_codes, hilbert_order, n_vtx_block);
  free(blk_hilbert_codes);


  /* Apply order to blk_ln_to_gn */
  PDM_g_num_t* sorted_blk_ln_to_gn = malloc( n_vtx_block * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_vtx_block; ++i) {
    sorted_blk_ln_to_gn[i] = blk_ln_to_gn[hilbert_order[i]];
  }
  free(blk_ln_to_gn);
  free(hilbert_order);

  PDM_g_num_t* distrib_blk_vtx = PDM_compute_entity_distribution(comm, n_vtx_block);
  PDM_block_to_part_t* btp = PDM_block_to_part_create(distrib_blk_vtx,
                              (const PDM_g_num_t **)  &ln_to_gn,
                                                      &dn_length,
                                                      1,
                                                      comm);

  int stride_one = 1;
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST,
                         &stride_one,
              (void *)   sorted_blk_ln_to_gn,
                         NULL,
              (void **) &ln_to_gn);
  PDM_block_to_part_free(btp);
  free(sorted_blk_ln_to_gn);
  free(distrib_blk_vtx);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
