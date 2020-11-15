
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_error.h"
#include "pdm_dmesh_nodal_to_dmesh_priv.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_part_to_block.h"
#include "pdm_logging.h"
#include "pdm_dmesh.h"
#include "pdm_quick_sort.h"
#include "pdm_para_graph_dual.h"

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

/**
 * \def _compute_keys
 */
static
void
_compute_keys
(
const int          n_face_elt_tot,
const int         *dcell_face_vtx_idx,
const PDM_g_num_t *dcell_face_vtx,
      PDM_g_num_t *ln_to_gn,
      PDM_g_num_t  key_mod
)
{
  for(int i_face = 0; i_face < n_face_elt_tot; ++i_face ) {
    PDM_g_num_t key = 0;
    for(int idx = dcell_face_vtx_idx[i_face]; idx < dcell_face_vtx_idx[i_face+1]; ++idx) {
      key += dcell_face_vtx[idx];
    }
    // min_vtx =
    ln_to_gn[i_face] = key % key_mod;
  }
}

static void
_make_absolute_face_numbering(_pdm_dmesh_nodal_t* mesh)
{

  PDM_g_num_t n_face_proc = mesh->dn_face;
  PDM_g_num_t beg_num_abs;

  PDM_MPI_Scan(&n_face_proc, &beg_num_abs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, mesh->pdm_mpi_comm);
  beg_num_abs -= n_face_proc;

  /** Compute the distribution of elements amont proc **/
  mesh->face_distrib = (PDM_g_num_t *) malloc((mesh->n_rank+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t _dn_face = (PDM_g_num_t) mesh->dn_face;
  PDM_MPI_Allgather((void *) &_dn_face,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&mesh->face_distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    mesh->pdm_mpi_comm);

  // mesh->face_distrib[0] = 1;
  mesh->face_distrib[0] = 0;

  for (int i = 1; i < mesh->n_rank+1; i++) {
    mesh->face_distrib[i] +=  mesh->face_distrib[i-1];
  }

  if (0 == 1) {
    printf("beg_num_abs::Face : "PDM_FMT_G_NUM" \n", beg_num_abs);
    printf("mesh->face_distrib : "PDM_FMT_G_NUM,  mesh->face_distrib[0]);
    for (int i = 1; i < mesh->n_rank+1; i++) {
      printf(" "PDM_FMT_G_NUM, mesh->face_distrib[i]);
    }
    printf("\n");
  }
}

static
void
_generate_query_entitiy
(
_pdm_dmesh_nodal_t *mesh,
int                 n_face_elt_tot,
PDM_g_num_t        *delmt_face_cell,
int                *dcell_face_vtx_idx,
PDM_g_num_t        *dcell_face_vtx
)
{
  /*
   * We are now all information flatten - we only need to compute hash_keys for each faces
   */
  PDM_g_num_t* ln_to_gn = (PDM_g_num_t*) malloc( n_face_elt_tot * sizeof(PDM_g_num_t));
  PDM_g_num_t key_mod = 4 * mesh->n_vtx_abs;
  // printf("key_mod::%i \n", key_mod);

  // Option des clés
  _compute_keys(n_face_elt_tot,
                dcell_face_vtx_idx,
                dcell_face_vtx,
                ln_to_gn,
                key_mod);

  if(0 == 1) {
    PDM_log_trace_array_long(ln_to_gn, n_face_elt_tot , "ln_to_gn:: ");
    PDM_log_trace_array_int(dcell_face_vtx_idx  , n_face_elt_tot , "dcell_face_vtx_idx:: ");
    PDM_log_trace_array_long(dcell_face_vtx, dcell_face_vtx_idx[n_face_elt_tot] , "dcell_face_vtx:: ");
  }

  /*
   * Prepare exchange by computing stride
   */
  int* dcell_face_vtx_n = (int        *) malloc( n_face_elt_tot * sizeof(int        ));
  int* stride_one       = (int        *) malloc( n_face_elt_tot * sizeof(int        ));
  for(int i_face = 0; i_face < n_face_elt_tot; ++i_face) {
    dcell_face_vtx_n[i_face] = dcell_face_vtx_idx[i_face+1] - dcell_face_vtx_idx[i_face];
    stride_one[i_face]       = 1;
  }
  free(dcell_face_vtx_idx);

  /*
   * Setup part_to_block to filter all keys
   */
  double* weight = (double *) malloc( n_face_elt_tot * sizeof(double));
  for(int i = 0; i < n_face_elt_tot; ++i) {
    weight[i] = 1.;
  }

  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &ln_to_gn,
                                                      &weight,
                                                      &n_face_elt_tot,
                                                      1,
                                                      mesh->pdm_mpi_comm);
  free(weight);
  /*
   * Exchange data
   */
  int*         blk_tot_face_vtx_n = NULL;
  PDM_g_num_t* blk_tot_face_vtx   = NULL;

  int blk_tot_face_vtx_size = PDM_part_to_block_exch(         ptb,
                                                              sizeof(PDM_g_num_t),
                                                              PDM_STRIDE_VAR,
                                                              -1,
                                                              &dcell_face_vtx_n,
                                                    (void **) &dcell_face_vtx,
                                                              &blk_tot_face_vtx_n,
                                                    (void **) &blk_tot_face_vtx);

  int* blk_n_face_per_key = NULL;
  int* blk_face_vtx_n     = NULL;
  int blk_face_vtx_n_size = PDM_part_to_block_exch(         ptb,
                                                            sizeof(int),
                                                            PDM_STRIDE_VAR,
                                                            -1,
                                                            &stride_one,
                                                  (void **) &dcell_face_vtx_n,
                                                            &blk_n_face_per_key,
                                                  (void **) &blk_face_vtx_n);
  free(dcell_face_vtx_n);

  int*         blk_elmt_face_cell_stri = NULL;
  PDM_g_num_t* blk_elmt_face_cell      = NULL;
  int blk_face_cell_size = PDM_part_to_block_exch(         ptb,
                                                           sizeof(PDM_g_num_t),
                                                           PDM_STRIDE_VAR,
                                                           -1,
                                                           &stride_one,
                                                 (void **) &delmt_face_cell,
                                                           &blk_elmt_face_cell_stri,
                                                 (void **) &blk_elmt_face_cell);

  /*
   *  Get the size of the current process bloc
   */
  int blk_size = PDM_part_to_block_n_elt_block_get(ptb);

  if( 0 == 1 ) {
    PDM_log_trace_array_int(blk_tot_face_vtx_n, blk_size             , "blk_tot_face_vtx_n:: ");
    PDM_log_trace_array_long(blk_tot_face_vtx , blk_tot_face_vtx_size, "blk_tot_face_vtx:: "  );

    PDM_log_trace_array_int(blk_n_face_per_key, blk_size         , "blk_n_face_per_key:: ");
    PDM_log_trace_array_int(blk_face_vtx_n    , blk_face_vtx_n_size, "blk_face_vtx_n:: ");

    PDM_log_trace_array_long(blk_elmt_face_cell, blk_face_cell_size, "blk_elmt_face_cell:: ");
  }

  PDM_part_to_block_free(ptb);
  free(stride_one);
  free(blk_elmt_face_cell_stri); // Same as blk_n_face_per_key
  free(blk_tot_face_vtx_n);

  /*
   * Get the max number of vertex of faces
   */
  int* blk_face_vtx_idx  = (int        *) malloc( (blk_face_vtx_n_size+1) * sizeof(int        ));
  int n_max_face_per_key = -1;
  for(int i_face = 0; i_face < blk_size; ++i_face) {
    n_max_face_per_key = PDM_MAX(n_max_face_per_key, blk_n_face_per_key[i_face]);
  }

  int n_max_vtx       = -1;
  blk_face_vtx_idx[0] = 0;
  for(int i_face = 0; i_face < blk_face_vtx_n_size; ++i_face) {
    n_max_vtx          = PDM_MAX(n_max_vtx         , blk_face_vtx_n    [i_face]);
    blk_face_vtx_idx[i_face+1] = blk_face_vtx_idx[i_face] + blk_face_vtx_n[i_face];
  }

  // PDM_log_trace_array_int(blk_face_vtx_idx, blk_face_vtx_n_size, "blk_face_vtx_idx:: ");
  /*
   * We need to identify each uniques faces :
   *      - We have multiple packet to treat
   *      - The connectivity can be sorted in place
   *      - Multiple case can occur :
   *           - Alone face normaly boundary
   *           - Multiple faces
   *           - Same face, we remove and replace by the first
   */
  PDM_g_num_t* loc_face_vtx_1 = (PDM_g_num_t *) malloc( n_max_vtx          * sizeof(PDM_g_num_t) );
  PDM_g_num_t* loc_face_vtx_2 = (PDM_g_num_t *) malloc( n_max_vtx          * sizeof(PDM_g_num_t) );
  int*         already_treat  = (int         *) malloc( n_max_face_per_key * sizeof(int        ) );

  /*
   * Allocate Memory - face_vtx - face_cell
   */
  mesh->_dface_vtx     = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  blk_tot_face_vtx_size  );
  mesh->_dface_vtx_idx = (int         *) malloc( sizeof(int        ) * (blk_face_vtx_n_size+1 ));
  mesh->_dface_cell    = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  blk_face_cell_size * 2 );

  printf("blk_face_cell_size::%i\n", blk_face_cell_size);

  /*
   * Init global numbering
   */
  int i_abs_face = 0;
  mesh->_dface_vtx_idx[0] = 0;

  int idx = 0;
  int idx_face_vtx = 0;
  for(int i_key = 0; i_key < blk_size; ++i_key) {
    // printf("Number of conflicting keys :: %i \n", blk_n_face_per_key[i_key]);

    int n_conflict_faces = blk_n_face_per_key[i_key];

    /* Reset */
    for(int j = 0; j < n_conflict_faces; ++j) {
      already_treat[j] = -1;
    }

    /* Loop over all faces in conflict */
    for(int i_face = 0; i_face < n_conflict_faces; ++i_face) {
      // printf("Number of vtx on faces %i :: %i with index [%i] \n", i_face, blk_face_vtx_n[idx+i_face], idx+i_face);

      int n_vtx_face_1 = blk_face_vtx_n[idx+i_face];
      int beg_1        = blk_face_vtx_idx[idx+i_face];
      if(already_treat[i_face] != 1) {

        PDM_g_num_t key_1 = 0;
        for(int j = 0; j < n_vtx_face_1; ++j) {
          loc_face_vtx_1[j] = blk_tot_face_vtx[beg_1+j];
          key_1 += loc_face_vtx_1[j];
        }
        PDM_quick_sort_long(loc_face_vtx_1, 0, n_vtx_face_1-1);

        for(int i_face_next = i_face+1; i_face_next < n_conflict_faces; ++i_face_next) {

          int n_vtx_face_2 = blk_face_vtx_n[idx+i_face_next];

          if(n_vtx_face_1 == n_vtx_face_2 ) {

            int beg_2 = blk_face_vtx_idx[idx+i_face_next];
            PDM_g_num_t key_2 = 0;
            for(int j = 0; j < n_vtx_face_1; ++j) {
              loc_face_vtx_2[j] = blk_tot_face_vtx[beg_2+j];
              key_2 += loc_face_vtx_2[j];
            }
            PDM_quick_sort_long(loc_face_vtx_2, 0, n_vtx_face_2-1);

            assert(key_1 == key_2);

            int is_same_face = 1;

            for(int i_vtx = 0; i_vtx < n_vtx_face_1; ++i_vtx) {
              if(loc_face_vtx_1[i_vtx] != loc_face_vtx_2[i_vtx]) {
                is_same_face = -1;
                break;
              }
            }

            if(is_same_face == 1 ){
              // printf(" It's a match ! \n");

              mesh->_dface_cell[2*i_abs_face  ] = blk_elmt_face_cell[idx+i_face     ];
              mesh->_dface_cell[2*i_abs_face+1] = blk_elmt_face_cell[idx+i_face_next];

              for(int i_vtx = 0; i_vtx < n_vtx_face_1; ++i_vtx) {
                mesh->_dface_vtx[idx_face_vtx++] = blk_tot_face_vtx[beg_1+i_vtx];
              }
              mesh->_dface_vtx_idx[i_abs_face+1] = mesh->_dface_vtx_idx[i_abs_face] + n_vtx_face_1;
              i_abs_face++;

              already_treat[i_face]      = 1;
              already_treat[i_face_next] = 1;
            }


          } /* End if same number of vertex */
        } /* End loop next face */
      } /* End loop already treated */

      /* If a face is not found a match inside the pool, it's a boundary condition */
      if(already_treat[i_face] != 1) {

        mesh->_dface_cell[2*i_abs_face  ] = blk_elmt_face_cell[idx+i_face];
        mesh->_dface_cell[2*i_abs_face+1] = 0;

        for(int i_vtx = 0; i_vtx < n_vtx_face_1; ++i_vtx) {
          mesh->_dface_vtx[idx_face_vtx++] = blk_tot_face_vtx[beg_1+i_vtx];
        }
        mesh->_dface_vtx_idx[i_abs_face+1] = mesh->_dface_vtx_idx[i_abs_face] + n_vtx_face_1;
        i_abs_face++;

        already_treat[i_face]      = 1;

      }

    } /* End loop face in conflict */

    idx += n_conflict_faces;

  }

  /*
   * Free all unused structure
   */
  free(loc_face_vtx_1);
  free(loc_face_vtx_2);
  free(already_treat);
  free(delmt_face_cell);
  free(dcell_face_vtx);
  free(ln_to_gn);
  free(blk_tot_face_vtx);
  free(blk_n_face_per_key);
  free(blk_face_vtx_n);
  free(blk_elmt_face_cell);

  if( 0 == 1 ){
    printf("i_abs_face::%i \n", i_abs_face);
    PDM_log_trace_array_int(mesh->_dface_vtx_idx, i_abs_face+1                    , "mesh->_dface_vtx_idx:: ");
    PDM_log_trace_array_long(mesh->_dface_vtx   , mesh->_dface_vtx_idx[i_abs_face], "_dface_vtx:: ");
    PDM_log_trace_array_long(mesh->_dface_cell  , 2*i_abs_face                    , "_dface_cell:: ");
  }

  /*
   * Fill up structure
   */
  mesh->dn_face = i_abs_face;

  /*
   * Realloc
   */
  mesh->_dface_vtx_idx = (int *        ) realloc((mesh->_dface_vtx_idx), (mesh->dn_face + 1) * sizeof(int * ) );
  mesh->_dface_vtx     = (PDM_g_num_t *) realloc((mesh->_dface_vtx    ), mesh->_dface_vtx_idx[mesh->dn_face] * sizeof(PDM_g_num_t * ));
  mesh->_dface_cell    = (PDM_g_num_t *) realloc((mesh->_dface_cell   ), mesh->dn_face * 2 * sizeof(PDM_g_num_t * ));

  /*
   * Generate absolute numerotation of faces
   */
  _make_absolute_face_numbering(mesh);

  /*
   * Rebuild cell face
   */
  int n_face_cell = 2*mesh->dn_face;

  int*         part_stri_face_cell = (int         *) malloc( sizeof(int        ) * n_face_cell );
  PDM_g_num_t* dface_cell_tmp      = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * n_face_cell );
  PDM_g_num_t* ln_to_gn_elem       = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * n_face_cell );

  int idx_g = 0;
  assert( mesh->dn_face == mesh->face_distrib[mesh->i_rank+1] - mesh->face_distrib[mesh->i_rank]);
  // for (int i = mesh->face_distrib[mesh->i_rank]; i < mesh->face_distrib[mesh->i_rank+1]; i++) {
  for (int i_face = 0; i_face < mesh->dn_face; ++i_face) {

    PDM_g_num_t g_num_face = (PDM_g_num_t) i_face + mesh->face_distrib[mesh->i_rank] + 1;

    dface_cell_tmp[idx_g] = mesh->_dface_cell[2*i_face];
    ln_to_gn_elem [idx_g] = g_num_face;
    part_stri_face_cell[idx_g] = 1;
    idx_g++;
    if(mesh->_dface_cell[2*i_face+1] != 0){
      dface_cell_tmp[idx_g] = mesh->_dface_cell[2*i_face+1];
      ln_to_gn_elem [idx_g] = g_num_face;
      part_stri_face_cell[idx_g] = 1;
      idx_g++;
    } else {
      dface_cell_tmp[idx_g] = mesh->_dface_cell[2*i_face];
      ln_to_gn_elem [idx_g] = g_num_face;
      part_stri_face_cell[idx_g] = 0;
      idx_g++;
    }
  }
  n_face_cell = idx_g; // Adapt size

  if(0 == 1 ){
    printf("n_face_cell::%i\n", n_face_cell);
    PDM_log_trace_array_int(part_stri_face_cell, n_face_cell, "part_stri_face_cell:: ");
    PDM_log_trace_array_long(ln_to_gn_elem     , n_face_cell, "ln_to_gn_elem:: ");
    PDM_log_trace_array_long(dface_cell_tmp    , n_face_cell, "dface_cell_tmp:: ");
  }

  /*
   *  Use part_to_block with the cell numbering
   */
  PDM_part_to_block_t *ptb2 = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       &dface_cell_tmp,
                                                       NULL,
                                                       &n_face_cell,
                                                       1,
                                                       mesh->pdm_mpi_comm);

  int         *blk_cell_face_n = NULL;
  PDM_g_num_t *blk_cell_face   = NULL;

  int blk_cell_face_size = PDM_part_to_block_exch(          ptb2,
                                                            sizeof(PDM_g_num_t),
                                                            PDM_STRIDE_VAR,
                                                            -1,
                                                            &part_stri_face_cell,
                                                  (void **) &ln_to_gn_elem,
                                                            &blk_cell_face_n,
                                                  (void **) &blk_cell_face);
  PDM_UNUSED(blk_cell_face_size);

  /*
   *  Get the size of the current process bloc
   */
  int delmt_tot = PDM_part_to_block_n_elt_block_get(ptb2);
  mesh->dn_cell = delmt_tot;

  /*
   * Free
   */
  PDM_part_to_block_free(ptb2);
  free(part_stri_face_cell);
  free(dface_cell_tmp     );
  free(ln_to_gn_elem      );

  /*
   * Allcoate
   */
  assert(mesh->dcell_face_idx == NULL);
  mesh->dcell_face_idx = (int * ) malloc( (delmt_tot + 1) * sizeof(int) );

  mesh->dcell_face_idx[0] = 0;
  for(int i = 0; i < delmt_tot; i++){
    mesh->dcell_face_idx[i+1] = mesh->dcell_face_idx[i] + blk_cell_face_n[i];
  }

  // PDM_log_trace_array_int (blk_cell_face_n, delmt_tot         , "blk_cell_face_n:: ");
  // PDM_log_trace_array_long(blk_cell_face  , blk_cell_face_size, "blk_cell_face:: ");

  assert(mesh->dcell_face == NULL);
  mesh->dcell_face = blk_cell_face;


  /* Compress connectivity in place */
  PDM_para_graph_compress_connectivity(mesh->dn_cell,
                                       mesh->dcell_face_idx,
                                       blk_cell_face_n,
                                       mesh->dcell_face);


  if( 1 == 1 ){
    printf("mesh->dn_cell ::%i\n", mesh->dn_cell );
    PDM_log_trace_array_int(mesh->dcell_face_idx, mesh->dn_cell+1                    , "mesh->dcell_face_idx:: ");
    PDM_log_trace_array_long(mesh->dcell_face   , mesh->dcell_face_idx[mesh->dn_cell], "dcell_face:: ");
  }

  /*
   *  Realloc
   */
  mesh->dcell_face = (PDM_g_num_t *) realloc( mesh->dcell_face, mesh->dcell_face_idx[mesh->dn_cell] * sizeof(PDM_g_num_t));


  free(blk_cell_face_n);
}


static
_pdm_dmesh_t*
_generate_faces_from_dmesh_nodal
(
  _pdm_dmesh_nodal_t *dmesh_nodal
)
{

  PDM_DMesh_nodal_t* dmn = (PDM_DMesh_nodal_t *) dmesh_nodal;
  PDM_UNUSED(dmesh_nodal);

  int n_face_elt_tot     = 0;
  int n_sum_vtx_face_tot = 0;

  PDM_dmesh_nodal_decompose_faces_get_size(dmn, &n_face_elt_tot, &n_sum_vtx_face_tot);

  PDM_g_num_t* delmt_face_cell    = (PDM_g_num_t*) malloc(  n_face_elt_tot     * sizeof(PDM_g_num_t));
  int*         dcell_face_vtx_idx = (int        *) malloc( (n_face_elt_tot +1) * sizeof(int        ));
  PDM_g_num_t* dcell_face_vtx     = (PDM_g_num_t*) malloc(  n_sum_vtx_face_tot * sizeof(PDM_g_num_t));

  dcell_face_vtx_idx[0] = 0;
  PDM_sections_decompose_faces(dmn,
                               dcell_face_vtx_idx,
                               dcell_face_vtx,
                               delmt_face_cell,
                               NULL, NULL);

  /*
   *  Create empty dmesh
   */
  PDM_dmesh_t* dm = PDM_dmesh_create(-1, -1, -1, -1, -1);

  printf("_generate_query_entitiy \n");
  _generate_query_entitiy(dmesh_nodal,
                          n_face_elt_tot,
                          delmt_face_cell,
                          dcell_face_vtx_idx,
                          dcell_face_vtx);
  printf("_generate_query_entitiy end \n");

  return (_pdm_dmesh_t *) dm;
}



static
_pdm_dmesh_t*
_generate_edges_from_dmesh_nodal
(
  _pdm_dmesh_nodal_t *dmesh_nodal
)
{

  PDM_DMesh_nodal_t* dmn = (PDM_DMesh_nodal_t *) dmesh_nodal;
  PDM_UNUSED(dmesh_nodal);

  int n_edge_elt_tot     = 0;
  int n_sum_vtx_edge_tot = 0;

  PDM_dmesh_nodal_decompose_edges_get_size(dmn, &n_edge_elt_tot, &n_sum_vtx_edge_tot);

  PDM_g_num_t* delmt_edge_cell    = (PDM_g_num_t*) malloc(  n_edge_elt_tot     * sizeof(PDM_g_num_t));
  int*         dcell_edge_vtx_idx = (int        *) malloc( (n_edge_elt_tot +1) * sizeof(int        ));
  PDM_g_num_t* dcell_edge_vtx     = (PDM_g_num_t*) malloc(  n_sum_vtx_edge_tot * sizeof(PDM_g_num_t));

  dcell_edge_vtx_idx[0] = 0;
  PDM_sections_decompose_edges(dmn,
                               dcell_edge_vtx_idx,
                               dcell_edge_vtx,
                               delmt_edge_cell,
                               NULL, NULL);

  /*
   *  Create empty dmesh
   */
  PDM_dmesh_t* dm = PDM_dmesh_create(-1, -1, -1, -1, -1);

  _generate_query_entitiy(dmesh_nodal,
                          n_edge_elt_tot,
                          delmt_edge_cell,
                          dcell_edge_vtx_idx,
                          dcell_edge_vtx);

  return (_pdm_dmesh_t *) dm;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_dmesh_nodal_to_dmesh_t*
PDM_dmesh_nodal_to_dmesh_create
(
const int             n_mesh,
const PDM_MPI_Comm    comm,
const PDM_ownership_t owner
)
{
  _pdm_dmesh_nodal_to_dmesh_t *dmn_to_dm = (_pdm_dmesh_nodal_to_dmesh_t *) malloc(sizeof(_pdm_dmesh_nodal_to_dmesh_t));

  dmn_to_dm->comm              = comm;
  dmn_to_dm->owner             = owner;
  dmn_to_dm->results_is_getted = PDM_FALSE;
  dmn_to_dm->n_mesh            = n_mesh;

  dmn_to_dm->dmesh_nodal     = (_pdm_dmesh_nodal_t **) malloc(sizeof(_pdm_dmesh_nodal_t *));
  dmn_to_dm->dmesh           = (_pdm_dmesh_t       **) malloc(sizeof(_pdm_dmesh_t       *));

  return (PDM_dmesh_nodal_to_dmesh_t *) dmn_to_dm;
}

/**
 * \brief  Add dmesh_nodal
 */
void
PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal
(
        PDM_dmesh_nodal_to_dmesh_t *dmn_to_dm,
  const int                         i_mesh,
        PDM_DMesh_nodal_t          *dmn
)
{
  _pdm_dmesh_nodal_to_dmesh_t* _dmn_to_dm = (_pdm_dmesh_nodal_to_dmesh_t *) dmn_to_dm;
  _dmn_to_dm->dmesh_nodal[i_mesh] = (_pdm_dmesh_nodal_t*) dmn;
}


void
PDM_dmesh_nodal_to_dmesh_get_dmesh
(
        PDM_dmesh_nodal_to_dmesh_t  *dmn_to_dm,
  const int                          i_mesh,
        PDM_dmesh_t                **dm
)
{
  _pdm_dmesh_nodal_to_dmesh_t* _dmn_to_dm = (_pdm_dmesh_nodal_to_dmesh_t *) dmn_to_dm;

  *dm = (PDM_dmesh_t *) _dmn_to_dm->dmesh[i_mesh];

  _dmn_to_dm->results_is_getted = PDM_TRUE;
}


/**
 * \brief  Free
 */
void
PDM_dmesh_nodal_to_dmesh_free
(
  PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm
)
{
  printf("PDM_dmesh_nodal_to_dmesh_free\n");
  _pdm_dmesh_nodal_to_dmesh_t* _dmn_to_dm = (_pdm_dmesh_nodal_to_dmesh_t *) dmn_to_dm;

  if(( _dmn_to_dm->owner == PDM_OWNERSHIP_KEEP ) ||
     ( _dmn_to_dm->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !_dmn_to_dm->results_is_getted)){

    for(int i_mesh = 0; i_mesh < _dmn_to_dm->n_mesh; ++i_mesh) {
      PDM_dmesh_free( (PDM_dmesh_t*)_dmn_to_dm->dmesh[i_mesh]);
    }
  }

  free(_dmn_to_dm->dmesh_nodal);
  free(_dmn_to_dm->dmesh);
}



void
PDM_dmesh_nodal_to_dmesh_compute
(
  PDM_dmesh_nodal_to_dmesh_t*                dmn_to_dm,
  const PDM_dmesh_nodal_to_dmesh_transform_t transform_kind
)
{
  _pdm_dmesh_nodal_to_dmesh_t* _dmn_to_dm = (_pdm_dmesh_nodal_to_dmesh_t *) dmn_to_dm;

  for(int i_mesh = 0; i_mesh < _dmn_to_dm->n_mesh; ++i_mesh) {

    switch (transform_kind) {

      case PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE:
        {
          _dmn_to_dm->dmesh[i_mesh] = _generate_faces_from_dmesh_nodal(_dmn_to_dm->dmesh_nodal[i_mesh]);
        }
        break;

      case PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE:
        {
          _dmn_to_dm->dmesh[i_mesh] = _generate_edges_from_dmesh_nodal(_dmn_to_dm->dmesh_nodal[i_mesh]);
        }
        break;
    }
  }

  // Boundary management

  // Join management


}
