
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
#include "pdm_binary_search.h"
#include "pdm_dmesh_nodal_to_dmesh_priv.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
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
const int         *delmt_face_vtx_idx,
const PDM_g_num_t *delmt_face_vtx,
      PDM_g_num_t *ln_to_gn,
      PDM_g_num_t  key_mod
)
{
  for(int i_face = 0; i_face < n_face_elt_tot; ++i_face ) {
    PDM_g_num_t key = 0;
    for(int idx = delmt_face_vtx_idx[i_face]; idx < delmt_face_vtx_idx[i_face+1]; ++idx) {
      key += delmt_face_vtx[idx];
    }
    // min_vtx =
    ln_to_gn[i_face] = key % key_mod;
  }
}

static
PDM_g_num_t*
_make_absolute_entity_numbering
(
       int          dn_entity,
 const PDM_MPI_Comm comm
)
{
  int n_rank;

  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_g_num_t n_entity_proc = dn_entity;
  PDM_g_num_t beg_num_abs;

  PDM_MPI_Scan(&n_entity_proc, &beg_num_abs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  beg_num_abs -= n_entity_proc;

  /** Compute the distribution of elements amont proc **/
  PDM_g_num_t *entity_distrib = (PDM_g_num_t *) malloc((n_rank+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t _dn_face = (PDM_g_num_t) dn_entity;
  PDM_MPI_Allgather((void *) &_dn_face,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&entity_distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    comm);

  // entity_distrib[0] = 1;
  entity_distrib[0] = 0;

  for (int i = 1; i < n_rank+1; i++) {
    entity_distrib[i] +=  entity_distrib[i-1];
  }

  if (0 == 1) {
    printf("beg_num_abs::Face : "PDM_FMT_G_NUM" \n", beg_num_abs);
    printf("entity_distrib : "PDM_FMT_G_NUM,  entity_distrib[0]);
    for (int i = 1; i < n_rank+1; i++) {
      printf(" "PDM_FMT_G_NUM, entity_distrib[i]);
    }
    printf("\n");
  }
  return entity_distrib;
}

static
void
_generate_entitiy_connectivity
(
PDM_dmesh_nodal_t   *mesh,
int                  n_entity_elt_tot,
PDM_g_num_t         *delmt_entity,
int                 *delmt_entity_vtx_idx,
PDM_g_num_t         *delmt_entity_vtx,
int                 *dn_entity,
PDM_g_num_t        **entity_distrib,
int                **dentity_vtx_idx,
PDM_g_num_t        **dentity_vtx,
int                **delmt_entity_out_idx,
PDM_g_num_t        **delmt_entity_out,
PDM_g_num_t        **dentity_elmt,
int                **dentity_elmt_idx
)
{
  // PDM_g_num_t       **dentity_elmt_idx --> face : dface_elemt :
  // PDM_g_num_t       **dentity_elmt_idx --> edge : dedge_elemt : --> Il faut un idx
  /*
   * We are now all information flatten - we only need to compute hash_keys for each entitys
   */
  PDM_g_num_t* ln_to_gn = (PDM_g_num_t*) malloc( n_entity_elt_tot * sizeof(PDM_g_num_t));
  PDM_g_num_t key_mod = 4 * mesh->n_vtx_abs;

  // Option des clés
  _compute_keys(n_entity_elt_tot,
                delmt_entity_vtx_idx,
                delmt_entity_vtx,
                ln_to_gn,
                key_mod);

  if(1 == 1) {
    PDM_log_trace_array_long(ln_to_gn, n_entity_elt_tot , "ln_to_gn:: ");
    PDM_log_trace_array_int(delmt_entity_vtx_idx  , n_entity_elt_tot , "delmt_entity_vtx_idx:: ");
    PDM_log_trace_array_long(delmt_entity_vtx, delmt_entity_vtx_idx[n_entity_elt_tot] , "delmt_entity_vtx:: ");
  }

  /*
   * Prepare exchange by computing stride
   */
  int* delmt_entity_vtx_n = (int        *) malloc( n_entity_elt_tot * sizeof(int        ));
  int* stride_one         = (int        *) malloc( n_entity_elt_tot * sizeof(int        ));
  for(int i_entity = 0; i_entity < n_entity_elt_tot; ++i_entity) {
    delmt_entity_vtx_n[i_entity] = delmt_entity_vtx_idx[i_entity+1] - delmt_entity_vtx_idx[i_entity];
    stride_one[i_entity]       = 1;
  }
  free(delmt_entity_vtx_idx);

  /*
   * Setup part_to_block to filter all keys
   */
  double* weight = (double *) malloc( n_entity_elt_tot * sizeof(double));
  for(int i = 0; i < n_entity_elt_tot; ++i) {
    weight[i] = 1.;
  }

  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &ln_to_gn,
                                                      &weight,
                                                      &n_entity_elt_tot,
                                                      1,
                                                      mesh->pdm_mpi_comm);
  free(weight);
  /*
   * Exchange data
   */
  int*         blk_tot_entity_vtx_n = NULL;
  PDM_g_num_t* blk_tot_entity_vtx   = NULL;

  int blk_tot_entity_vtx_size = PDM_part_to_block_exch(         ptb,
                                                              sizeof(PDM_g_num_t),
                                                              PDM_STRIDE_VAR,
                                                              -1,
                                                              &delmt_entity_vtx_n,
                                                    (void **) &delmt_entity_vtx,
                                                              &blk_tot_entity_vtx_n,
                                                    (void **) &blk_tot_entity_vtx);

  int* blk_n_entity_per_key = NULL;
  int* blk_entity_vtx_n     = NULL;
  int blk_entity_vtx_n_size = PDM_part_to_block_exch(         ptb,
                                                            sizeof(int),
                                                            PDM_STRIDE_VAR,
                                                            -1,
                                                            &stride_one,
                                                  (void **) &delmt_entity_vtx_n,
                                                            &blk_n_entity_per_key,
                                                  (void **) &blk_entity_vtx_n);
  free(delmt_entity_vtx_n);

  int*         blk_elmt_entity_elmt_stri = NULL;
  PDM_g_num_t* blk_elmt_entity_elmt      = NULL;
  int blk_entity_elmt_size = PDM_part_to_block_exch(         ptb,
                                                           sizeof(PDM_g_num_t),
                                                           PDM_STRIDE_VAR,
                                                           -1,
                                                           &stride_one,
                                                 (void **) &delmt_entity,
                                                           &blk_elmt_entity_elmt_stri,
                                                 (void **) &blk_elmt_entity_elmt);

  /*
   *  Get the size of the current process bloc
   */
  int blk_size = PDM_part_to_block_n_elt_block_get(ptb);

  if( 0 == 1 ) {
    PDM_log_trace_array_int(blk_tot_entity_vtx_n, blk_size             , "blk_tot_entity_vtx_n:: ");
    PDM_log_trace_array_long(blk_tot_entity_vtx , blk_tot_entity_vtx_size, "blk_tot_entity_vtx:: "  );

    PDM_log_trace_array_int(blk_n_entity_per_key, blk_size         , "blk_n_entity_per_key:: ");
    PDM_log_trace_array_int(blk_entity_vtx_n    , blk_entity_vtx_n_size, "blk_entity_vtx_n:: ");

    PDM_log_trace_array_long(blk_elmt_entity_elmt, blk_entity_elmt_size, "blk_elmt_entity_elmt:: ");
  }

  PDM_part_to_block_free(ptb);
  free(stride_one);
  free(blk_elmt_entity_elmt_stri); // Same as blk_n_entity_per_key
  free(blk_tot_entity_vtx_n);

  /*
   * Get the max number of vertex of entitys
   */
  int* blk_entity_vtx_idx  = (int        *) malloc( (blk_entity_vtx_n_size+1) * sizeof(int        ));
  int n_max_entity_per_key = -1;
  int n_tot_entity_per_key = 0;
  for(int i_entity = 0; i_entity < blk_size; ++i_entity) {
    n_max_entity_per_key = PDM_MAX(n_max_entity_per_key, blk_n_entity_per_key[i_entity]);
    n_tot_entity_per_key += blk_n_entity_per_key[i_entity];
  }

  int n_max_vtx       = -1;
  blk_entity_vtx_idx[0] = 0;
  for(int i_entity = 0; i_entity < blk_entity_vtx_n_size; ++i_entity) {
    n_max_vtx          = PDM_MAX(n_max_vtx         , blk_entity_vtx_n    [i_entity]);
    blk_entity_vtx_idx[i_entity+1] = blk_entity_vtx_idx[i_entity] + blk_entity_vtx_n[i_entity];
  }

  // PDM_log_trace_array_int(blk_entity_vtx_idx, blk_entity_vtx_n_size, "blk_entity_vtx_idx:: ");
  /*
   * We need to identify each uniques entitys :
   *      - We have multiple packet to treat
   *      - The connectivity can be sorted in place
   *      - Multiple case can occur :
   *           - Alone entity normaly boundary
   *           - Multiple entitys
   *           - Same entity, we remove and replace by the first
   */
  PDM_g_num_t* loc_entity_vtx_1 = (PDM_g_num_t *) malloc(  n_max_vtx               * sizeof(PDM_g_num_t) );
  PDM_g_num_t* loc_entity_vtx_2 = (PDM_g_num_t *) malloc(  n_max_vtx               * sizeof(PDM_g_num_t) );
  int*         already_treat    = (int         *) malloc(  n_max_entity_per_key    * sizeof(int        ) );
  int*         same_entity_idx  = (int         *) malloc( (n_max_entity_per_key+1) * sizeof(int        ) );
  int*         sens_entity      = (int         *) malloc(  n_max_entity_per_key    * sizeof(int        ) );

  /*
   * Allocate Memory - entity_vtx - entity_elmt
   */
  *dentity_vtx      = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  blk_tot_entity_vtx_size  );
  *dentity_vtx_idx  = (int         *) malloc( sizeof(int        ) * (blk_entity_vtx_n_size+1 ));
  // *dentity_elmt     = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  blk_entity_elmt_size * 2 );
  *dentity_elmt     = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  n_tot_entity_per_key );
  *dentity_elmt_idx = (int         *) malloc( sizeof(int)         * (blk_entity_elmt_size+1)  );

  PDM_g_num_t *_dentity_vtx      = *dentity_vtx;
  int         *_dentity_vtx_idx  = *dentity_vtx_idx;
  PDM_g_num_t *_dentity_elmt     = *dentity_elmt;
  int         *_dentity_elmt_idx = *dentity_elmt_idx;

  printf("blk_entity_elmt_size::%i\n", blk_entity_elmt_size);
  printf("n_tot_entity_per_key::%i\n", n_tot_entity_per_key);

  /*
   * Init global numbering
   */
  int i_abs_entity = 0;
  _dentity_vtx_idx[0] = 0;

  int idx = 0;
  int idx_entity_vtx = 0;
  _dentity_elmt_idx[0] = 0;
  for(int i_key = 0; i_key < blk_size; ++i_key) {
    // printf(" --- Number of conflicting keys :: %i \n", blk_n_entity_per_key[i_key]);

    int n_conflict_entitys = blk_n_entity_per_key[i_key];

    /* Reset */
    for(int j = 0; j < n_conflict_entitys; ++j) {
      already_treat[j] = -1;
    }

    /* Loop over all entitys in conflict */
    for(int i_entity = 0; i_entity < n_conflict_entitys; ++i_entity) {
      // printf("Number of vtx on entitys %i :: %i with index [%i] \n", i_entity, blk_entity_vtx_n[idx+i_entity], idx+i_entity);

      int n_vtx_entity_1 = blk_entity_vtx_n  [idx+i_entity];
      int beg_1          = blk_entity_vtx_idx[idx+i_entity];
      int idx_next_same_entity = 0;
      sens_entity[idx_next_same_entity] = 1;
      same_entity_idx[idx_next_same_entity++] = i_entity;

      if(already_treat[i_entity] != 1) {

        PDM_g_num_t key_1 = 0;
        int idx_min_1 = -1;
        PDM_g_num_t min_1 = mesh->n_vtx_abs+1;
        for(int j = 0; j < n_vtx_entity_1; ++j) {
          loc_entity_vtx_1[j] = blk_tot_entity_vtx[beg_1+j];
          key_1 += loc_entity_vtx_1[j];
          if(loc_entity_vtx_1[j] < min_1) {
            min_1 = loc_entity_vtx_1[j];
            idx_min_1 = j;
          };
        }
        PDM_quick_sort_long(loc_entity_vtx_1, 0, n_vtx_entity_1-1);

        for(int i_entity_next = i_entity+1; i_entity_next < n_conflict_entitys; ++i_entity_next) {

          int n_vtx_entity_2 = blk_entity_vtx_n[idx+i_entity_next];

          if(n_vtx_entity_1 == n_vtx_entity_2 ) {

            int beg_2 = blk_entity_vtx_idx[idx+i_entity_next];
            PDM_g_num_t key_2 = 0;
            int idx_min_2 = -1;
            PDM_g_num_t min_2 = mesh->n_vtx_abs+1;
            for(int j = 0; j < n_vtx_entity_1; ++j) {
              loc_entity_vtx_2[j] = blk_tot_entity_vtx[beg_2+j];
              key_2 += loc_entity_vtx_2[j];
              if(loc_entity_vtx_2[j] < min_2) {
                min_2 = loc_entity_vtx_2[j];
                idx_min_2 = j;
              };
            }
            PDM_quick_sort_long(loc_entity_vtx_2, 0, n_vtx_entity_2-1);

            assert(key_1 == key_2);

            int is_same_entity = 1;
            for(int i_vtx = 0; i_vtx < n_vtx_entity_1; ++i_vtx) {
              if(loc_entity_vtx_1[i_vtx] != loc_entity_vtx_2[i_vtx]) {
                is_same_entity = -1;
                break;
              }
            }

            if(is_same_entity == 1 ){
              // printf("idx_min_1 = %i | idx_min_2 = %i \n", idx_min_1, idx_min_2);
              // printf(" test1 :: ");
              // for(int i_vtx = 0; i_vtx < n_vtx_entity_1; ++i_vtx) {
              //   printf(" %i", (int)blk_tot_entity_vtx[beg_1+i_vtx]);
              // }
              // printf(" \n");
              // printf(" test2 :: ");
              // for(int i_vtx = 0; i_vtx < n_vtx_entity_1; ++i_vtx) {
              //   printf(" %i", (int)blk_tot_entity_vtx[beg_2+i_vtx]);
              // }
              // printf(" \n");

              // Determine the sens
              PDM_g_num_t i1 = blk_tot_entity_vtx[beg_1 +  idx_min_1                    ];
              PDM_g_num_t i2 = blk_tot_entity_vtx[beg_1 + (idx_min_1+1) % n_vtx_entity_1];

              PDM_g_num_t j1 = blk_tot_entity_vtx[beg_2 +  idx_min_2                    ];
              PDM_g_num_t j2 = blk_tot_entity_vtx[beg_2 + (idx_min_2+1) % n_vtx_entity_1];

              assert(i1 == j1); // Panic
              if(i2 != j2) {
                sens_entity[idx_next_same_entity] = -1;
                // printf(" idx3 = %i \n", (idx_min_1+n_vtx_entity_1) % n_vtx_entity_1);
                PDM_g_num_t i3 = blk_tot_entity_vtx[beg_1 + (idx_min_1+n_vtx_entity_1-1) % n_vtx_entity_1];
                // printf(" i1 = %i | i2 = %i | i3 = %i | j1 = %i | j2 = %i\n", i1, i2, i3, j1, j2);
                assert(i3 == j2);
              } else {
                sens_entity[idx_next_same_entity] = 1;
                assert(i2 == j2);
              }

              // Check if same sens
              same_entity_idx[idx_next_same_entity++] = i_entity_next;
            }
          } /* End if same number of vertex */
        } /* End loop next entity */

         // printf("[%i] - same_entity_idx::", idx_next_same_entity);
         // for(int ii = 0; ii < idx_next_same_entity; ++ii) {
         //   printf(" %i", same_entity_idx[ii]);
         // }
         // printf("\n");

        _dentity_vtx_idx[i_abs_entity+1] = _dentity_vtx_idx[i_abs_entity] + n_vtx_entity_1;
        for(int i_vtx = 0; i_vtx < n_vtx_entity_1; ++i_vtx) {
          _dentity_vtx[idx_entity_vtx++] = blk_tot_entity_vtx[beg_1+i_vtx];
          // Ecriture à partir du minumun
          // _dentity_vtx[idx_entity_vtx++] = blk_tot_entity_vtx[beg_1+i_vtx];
        }

        _dentity_elmt_idx[i_abs_entity+1] = _dentity_elmt_idx[i_abs_entity];
        for(int i = 0; i < idx_next_same_entity; ++i) {
          int i_same_entity = same_entity_idx[i];
          int sign = 1; // sens_entity[i];
          // Signe à faire
          _dentity_elmt[_dentity_elmt_idx[i_abs_entity+1]++] = sign*blk_elmt_entity_elmt[idx+i_same_entity];
          already_treat[i_same_entity] = 1;
        }
        i_abs_entity++;

      } /* End loop already treated */

    } /* End loop entity in conflict */

    idx += n_conflict_entitys;

  }

  // printf(" realloc dentity_elmt : %i --> %i \n", n_tot_entity_per_key, _dentity_elmt_idx[i_abs_entity]);
  *dentity_elmt = realloc(*dentity_elmt, sizeof(PDM_g_num_t) *  _dentity_elmt_idx[i_abs_entity] );
  _dentity_elmt = *dentity_elmt;

  /*
   * Free all unused structure
   */
  free(loc_entity_vtx_1);
  free(loc_entity_vtx_2);
  free(already_treat);
  free(same_entity_idx);
  free(sens_entity);
  free(delmt_entity);
  free(delmt_entity_vtx);
  free(ln_to_gn);
  free(blk_tot_entity_vtx);
  free(blk_n_entity_per_key);
  free(blk_entity_vtx_n);
  free(blk_elmt_entity_elmt);

  if( 1 == 1 ){
    printf("i_abs_entity::%i \n", i_abs_entity+1);
    PDM_log_trace_array_int(_dentity_vtx_idx, i_abs_entity+1                   , "_dentity_vtx_idx:: " );
    PDM_log_trace_array_long(_dentity_vtx   , _dentity_vtx_idx[i_abs_entity]   , "_dentity_vtx:: "     );
    PDM_log_trace_array_int(_dentity_elmt_idx, i_abs_entity+1                  , "_dentity_elmt_idx:: ");
    PDM_log_trace_array_long(_dentity_elmt    , _dentity_elmt_idx[i_abs_entity], "_dentity_elmt:: "    );
  }

  /*
   * Fill up structure
   */
  *dn_entity = i_abs_entity;
  int _dn_entity = *dn_entity;

  /*
   * Realloc
   */
  *dentity_vtx_idx = (int *        ) realloc(*dentity_vtx_idx, (_dn_entity + 1) * sizeof(int * ) );
  _dentity_vtx_idx = *dentity_vtx_idx;

  *dentity_vtx     = (PDM_g_num_t *) realloc(*dentity_vtx    , _dentity_vtx_idx[_dn_entity] * sizeof(PDM_g_num_t * ));
  _dentity_vtx     = *dentity_vtx;

  /*
   * Generate absolute numerotation of entitys
   */
  *entity_distrib = _make_absolute_entity_numbering(_dn_entity, mesh->pdm_mpi_comm);
  PDM_g_num_t* _entity_distrib = *entity_distrib;

  /*
   * Rebuild elmt entity
   */
  int n_entity_elmt = _dentity_elmt_idx[_dn_entity];

  int*         part_stri_entity_elmt = (int         *) malloc( sizeof(int        ) * n_entity_elmt );
  PDM_g_num_t* ln_to_gn_elem         = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * n_entity_elmt );

  assert( _dn_entity == _entity_distrib[mesh->i_rank+1] - _entity_distrib[mesh->i_rank]);

  for (int i_entity = 0; i_entity < _dn_entity; ++i_entity) {
    for(int i_data = _dentity_elmt_idx[i_entity]; i_data < _dentity_elmt_idx[i_entity+1]; ++i_data) {
      PDM_g_num_t g_num_entity = (PDM_g_num_t) i_entity + _entity_distrib[mesh->i_rank] + 1;
      ln_to_gn_elem        [i_data] = g_num_entity;
      part_stri_entity_elmt[i_data] = 1;
    }
  }

  if(0 == 1 ){
    printf("n_entity_elmt::%i\n", n_entity_elmt);
    PDM_log_trace_array_int(part_stri_entity_elmt, n_entity_elmt, "part_stri_entity_elmt:: ");
    PDM_log_trace_array_long(ln_to_gn_elem       , n_entity_elmt, "ln_to_gn_elem:: ");
    PDM_log_trace_array_long(_dentity_elmt       , n_entity_elmt, "_dentity_elmt:: ");
  }

  /*
   *  Use part_to_block with the elmt numbering
   */
  PDM_part_to_block_t *ptb2 = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       &_dentity_elmt,
                                                       NULL,
                                                       &n_entity_elmt,
                                                       1,
                                                       mesh->pdm_mpi_comm);

  int         *blk_elmt_entity_n = NULL;
  PDM_g_num_t *blk_elmt_entity   = NULL;

  int blk_elmt_entity_size = PDM_part_to_block_exch(          ptb2,
                                                            sizeof(PDM_g_num_t),
                                                            PDM_STRIDE_VAR,
                                                            -1,
                                                            &part_stri_entity_elmt,
                                                  (void **) &ln_to_gn_elem,
                                                            &blk_elmt_entity_n,
                                                  (void **) &blk_elmt_entity);
  PDM_UNUSED(blk_elmt_entity_size);

  /*
   *  Get the size of the current process bloc
   */
  int delmt_tot = PDM_part_to_block_n_elt_block_get(ptb2);
  mesh->dn_elmt = delmt_tot;

  /*
   * Free
   */
  PDM_part_to_block_free(ptb2);
  free(part_stri_entity_elmt);
  free(ln_to_gn_elem      );

  /*
   * Allcoate
   */
  assert(*delmt_entity_out_idx == NULL);
  *delmt_entity_out_idx = (int * ) malloc( (delmt_tot + 1) * sizeof(int) );
  int* _delmt_entity_out_idx = *delmt_entity_out_idx;

  _delmt_entity_out_idx[0] = 0;
  for(int i = 0; i < delmt_tot; i++){
    _delmt_entity_out_idx[i+1] = _delmt_entity_out_idx[i] + blk_elmt_entity_n[i];
  }

  // PDM_log_trace_array_int (blk_elmt_entity_n, delmt_tot         , "blk_elmt_entity_n:: ");
  // PDM_log_trace_array_long(blk_elmt_entity  , blk_elmt_entity_size, "blk_elmt_entity:: ");

  *delmt_entity_out = blk_elmt_entity;
  PDM_g_num_t *_delmt_entity_out = *delmt_entity_out;

  /* Compress connectivity in place */
  PDM_para_graph_compress_connectivity(mesh->dn_elmt,
                                       _delmt_entity_out_idx,
                                       blk_elmt_entity_n,
                                       _delmt_entity_out);


  if( 1 == 1 ){
    printf("mesh->dn_elmt ::%i\n", mesh->dn_elmt );
    PDM_log_trace_array_int(_delmt_entity_out_idx, mesh->dn_elmt+1                     , "_delmt_entity_out_idx:: ");
    PDM_log_trace_array_long(_delmt_entity_out   , _delmt_entity_out_idx[mesh->dn_elmt], "delmt_entity:: ");
  }

  /*
   *  Realloc
   */
  *delmt_entity_out = (PDM_g_num_t *) realloc( *delmt_entity_out, _delmt_entity_out_idx[mesh->dn_elmt] * sizeof(PDM_g_num_t));

  free(blk_elmt_entity_n);
}


static
PDM_dmesh_t*
_generate_faces_from_dmesh_nodal
(
  PDM_dmesh_nodal_t *dmesh_nodal
)
{

  PDM_dmesh_nodal_t* dmn = (PDM_dmesh_nodal_t *) dmesh_nodal;

  int n_face_elt_tot     = 0;
  int n_sum_vtx_face_tot = 0;

  PDM_dmesh_nodal_decompose_faces_get_size(dmn, &n_face_elt_tot, &n_sum_vtx_face_tot);

  PDM_g_num_t* delmt_face         = (PDM_g_num_t*) malloc(  n_face_elt_tot     * sizeof(PDM_g_num_t));
  int*         delmt_face_vtx_idx = (int        *) malloc( (n_face_elt_tot +1) * sizeof(int        ));
  PDM_g_num_t* delmt_face_vtx     = (PDM_g_num_t*) malloc(  n_sum_vtx_face_tot * sizeof(PDM_g_num_t));

  delmt_face_vtx_idx[0] = 0;
  PDM_sections_decompose_faces(dmn,
                               delmt_face_vtx_idx,
                               delmt_face_vtx,
                               delmt_face,
                               NULL, NULL);

  /*
   *  Create empty dmesh
   */
  PDM_dmesh_t* dm = PDM_dmesh_create(-1, -1, -1, -1, -1);

  int               *dentity_elmt_idx;
  _generate_entitiy_connectivity(dmesh_nodal,
                                 n_face_elt_tot,
                                 delmt_face,
                                 delmt_face_vtx_idx,
                                 delmt_face_vtx,
                                 &dmesh_nodal->dn_face,
                                 &dmesh_nodal->face_distrib,
                                 &dmesh_nodal->_dface_vtx_idx,
                                 &dmesh_nodal->_dface_vtx,
                                 &dmesh_nodal->dcell_face_idx,
                                 &dmesh_nodal->dcell_face,
                                 &dmesh_nodal->_dface_cell,
                                 &dentity_elmt_idx);
  free(dentity_elmt_idx);


  // Count the number of edges
  int n_edge_elt_tot = dmesh_nodal->_dface_vtx_idx[dmesh_nodal->dn_face];

  int*         dface_edge_vtx_idx = (int         *) malloc( ( n_edge_elt_tot + 1) * sizeof(int        ) );
  PDM_g_num_t* dface_edge         = (PDM_g_num_t *) malloc(     n_edge_elt_tot    * sizeof(PDM_g_num_t) );
  PDM_g_num_t* dface_edge_vtx     = (PDM_g_num_t *) malloc( 2 * n_edge_elt_tot    * sizeof(PDM_g_num_t) );
  int idx = 0;
  int i_edge = 0;
  dface_edge_vtx_idx[0] = 0;
  for(int i_face = 0; i_face < dmesh_nodal->dn_face; ++i_face ) {
    int beg        = dmesh_nodal->_dface_vtx_idx[i_face  ];
    int n_vtx_elmt = dmesh_nodal->_dface_vtx_idx[i_face+1] - beg;

    for(int ivtx = 0; ivtx < n_vtx_elmt; ++ivtx ) {
      int inext = (ivtx + 1) % n_vtx_elmt;

      dface_edge_vtx[idx++] = dmesh_nodal->_dface_vtx[beg+ivtx ];
      dface_edge_vtx[idx++] = dmesh_nodal->_dface_vtx[beg+inext];

      dface_edge_vtx_idx[i_edge+1] = dface_edge_vtx_idx[i_edge] + 2;
      dface_edge[i_edge] = (PDM_g_num_t) i_face + dmesh_nodal->face_distrib[dmesh_nodal->i_rank] + 1;
      i_edge++;
    }
  }

  if( 1 == 1 ){
    printf("n_edge_elt_tot ::%i\n", n_edge_elt_tot );
    PDM_log_trace_array_int (dface_edge_vtx_idx, n_edge_elt_tot+1              , "dface_edge_vtx_idx:: ");
    PDM_log_trace_array_long(dface_edge_vtx    , dface_edge_vtx_idx[n_edge_elt_tot], "dface_edge_vtx:: ");
    PDM_log_trace_array_long(dface_edge        , n_edge_elt_tot                , "dface_edge:: ");
  }

  int* dedge_face_idx = NULL;
  _generate_entitiy_connectivity(dmesh_nodal,
                                 n_edge_elt_tot,
                                 dface_edge,
                                 dface_edge_vtx_idx,
                                 dface_edge_vtx,
                                 &dmesh_nodal->dn_edge,
                                 &dmesh_nodal->edge_distrib,
                                 &dmesh_nodal->_dedge_vtx_idx,
                                 &dmesh_nodal->_dedge_vtx,
                                 &dmesh_nodal->dface_edge_idx,
                                 &dmesh_nodal->dface_edge,
                                 &dmesh_nodal->_dedge_face,
                                 &dedge_face_idx);

  if( 1 == 1 ){
    printf("dmesh_nodal->dn_edge ::%i\n", dmesh_nodal->dn_edge );
    PDM_log_trace_array_int (dmesh_nodal->_dedge_vtx_idx, dmesh_nodal->dn_edge+1                           , "dmesh_nodal->_dedge_vtx_idx:: ");
    PDM_log_trace_array_long(dmesh_nodal->_dedge_vtx    , dmesh_nodal->_dedge_vtx_idx[dmesh_nodal->dn_edge], "dmesh_nodal->_dedge_vtx:: ");

    PDM_log_trace_array_int (dmesh_nodal->dface_edge_idx, dmesh_nodal->dn_face+1                           , "dmesh_nodal->dface_edge_idx:: ");
    PDM_log_trace_array_long(dmesh_nodal->dface_edge    , dmesh_nodal->dface_edge_idx[dmesh_nodal->dn_face], "dmesh_nodal->dface_edge:: ");

    PDM_log_trace_array_int (dedge_face_idx, dmesh_nodal->dn_edge+1                           , "dmesh_nodal->dedge_face_idx:: ");
    PDM_log_trace_array_long(dmesh_nodal->_dedge_face    , dedge_face_idx[dmesh_nodal->dn_edge], "dmesh_nodal->_dedge_face:: ");

  }

  // free(dface_edge_vtx_idx);
  // free(dface_edge_vtx);
  // free(dface_edge);
  free(dedge_face_idx);
  // dedge_face + dface_cell -> dedge_cell

  return dm;
}



static
PDM_dmesh_t*
_generate_edges_from_dmesh_nodal
(
  PDM_dmesh_nodal_t *dmesh_nodal
)
{
  PDM_dmesh_nodal_t* dmn = (PDM_dmesh_nodal_t *) dmesh_nodal;

  int n_edge_elt_tot     = 0;
  int n_sum_vtx_edge_tot = 0;

  PDM_dmesh_nodal_decompose_edges_get_size(dmn, &n_edge_elt_tot, &n_sum_vtx_edge_tot);

  PDM_g_num_t* delmt_edge         = (PDM_g_num_t*) malloc(  n_edge_elt_tot     * sizeof(PDM_g_num_t));
  int*         delmt_edge_vtx_idx = (int        *) malloc( (n_edge_elt_tot +1) * sizeof(int        ));
  PDM_g_num_t* delmt_edge_vtx     = (PDM_g_num_t*) malloc(  n_sum_vtx_edge_tot * sizeof(PDM_g_num_t));

  delmt_edge_vtx_idx[0] = 0;
  PDM_sections_decompose_edges(dmn,
                               delmt_edge_vtx_idx,
                               delmt_edge_vtx,
                               delmt_edge,
                               NULL, NULL);

  /*
   *  Create empty dmesh
   */
  PDM_dmesh_t* dm = PDM_dmesh_create(-1, -1, -1, -1, -1);

  int               *dentity_elmt_idx;
  _generate_entitiy_connectivity(dmesh_nodal,
                                 n_edge_elt_tot,
                                 delmt_edge,
                                 delmt_edge_vtx_idx,
                                 delmt_edge_vtx,
                                 &dmesh_nodal->dn_face,
                                 &dmesh_nodal->face_distrib,
                                 &dmesh_nodal->_dface_vtx_idx,
                                 &dmesh_nodal->_dface_vtx,
                                 &dmesh_nodal->dcell_face_idx,
                                 &dmesh_nodal->dcell_face,
                                 &dmesh_nodal->_dface_cell,
                                 &dentity_elmt_idx);
  free(dentity_elmt_idx);

  return dm;
}

static
void
_translate_element_group_to_faces
(
  PDM_dmesh_nodal_t *dmesh_nodal,
  PDM_dmesh_t       *dm
)
{

  printf("_translate_element_group_to_faces \n");
  assert(dmesh_nodal != NULL);
  assert(dm          != NULL);

  // Ok we have :
  //     -> delmt_entity
  //     -> delmt_entity_idx
  //     -> dgroup_elmt
  //     -> dgroup_elmt_idx
  // We need : dentity_group and dentity_group_idx
  // In order to apply in parallel the table delmt_enitiy we need to have
  // a block_data by elmt

  /*
   * Prepare exchange protocol
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(dmesh_nodal->face_distrib,
                               (const PDM_g_num_t **) &dmesh_nodal->dgroup_elmt,
                                                      &dmesh_nodal->dgroup_elmt_idx[dmesh_nodal->n_group_elmt],
                                                      1,
                                                      dmesh_nodal->pdm_mpi_comm);

  /*
   * Exchange
   */
  int dn_face = dmesh_nodal->face_distrib[dmesh_nodal->i_rank+1] - dmesh_nodal->face_distrib[dmesh_nodal->i_rank];
  int* dcell_face_n = (int *) malloc( dn_face * sizeof(int));
  for(int i = 0; i < dn_face; ++i) {
    dcell_face_n[i] = dmesh_nodal->dcell_face_idx[i+1] - dmesh_nodal->dcell_face_idx[i];
  }

  int**         part_group_stri;
  PDM_g_num_t** part_group_data;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          dcell_face_n,
             (void *  )   dmesh_nodal->dcell_face,
             (int  ***)  &part_group_stri,
             (void ***)  &part_group_data);
  free(dcell_face_n);

  int*         _part_group_stri = part_group_stri[0];
  PDM_g_num_t* _part_group_data = part_group_data[0];

  // int idx_data = 0;
  // for(int i = 0; i < dmesh_nodal->dgroup_elmt_idx[dmesh_nodal->n_group_elmt]; ++i) {
  //   printf("_part_group_stri[%i] = %i \n", i, _part_group_stri[i]);
  //   for(int i_data = 0; i_data < _part_group_stri[i]; ++i_data) {
  //     printf("  -> _part_group_data[%i] = "PDM_FMT_G_NUM" \n", i, _part_group_data[idx_data++]);
  //   }
  // }

  assert(dm->_dface_bound     == NULL);
  assert(dm->_dface_bound_idx == NULL);

  // dm->n_bnd = dmesh_nodal->n_group_elmt;
  // dm->_dface_bound_idx = (int * ) malloc( (dm->n_bnd+1) * sizeof(int) );

  // int idx_stri = 0;
  // dm->_dface_bound_idx[0] = 0;
  // for(int i_group = 0; i_group < dm->n_bnd; ++i_group) {
  //   dm->_dface_bound_idx[i_group+1] = dm->_dface_bound_idx[i_group];
  //   for(int ielmt = dmesh_nodal->dgroup_elmt_idx[i_group]; ielmt < dmesh_nodal->dgroup_elmt_idx[i_group+1]; ++ielmt) {
  //     dm->_dface_bound_idx[i_group+1] += _part_group_stri[idx_stri++];
  //   }
  // }

  printf("_translate_element_group_to_faces is done but not transfer to dmesh = Leaks or no results !!! \n");
  // dm->_dface_bound = _part_group_data;

  // if(1 == 1) {
  //   PDM_log_trace_array_int (dm->_dface_bound_idx, dm->n_bnd+1                    , "dm->_dface_bound_idx:: ");
  //   PDM_log_trace_array_long(dm->_dface_bound    , dm->_dface_bound_idx[dm->n_bnd], "dm->_dface_bound:: ");
  // }

  free(part_group_stri);
  free(_part_group_stri);
  free(part_group_data);
  free(_part_group_data); // TO Remove when dmesh is OK

  PDM_block_to_part_free(btp);
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
  PDM_dmesh_nodal_to_dmesh_t *dmn_to_dm = (PDM_dmesh_nodal_to_dmesh_t *) malloc(sizeof(PDM_dmesh_nodal_to_dmesh_t));

  dmn_to_dm->comm              = comm;
  dmn_to_dm->owner             = owner;
  dmn_to_dm->results_is_getted = PDM_FALSE;
  dmn_to_dm->n_mesh            = n_mesh;

  dmn_to_dm->dmesh_nodal     = (PDM_dmesh_nodal_t **) malloc(sizeof(PDM_dmesh_nodal_t *));
  dmn_to_dm->dmesh           = (PDM_dmesh_t       **) malloc(sizeof(PDM_dmesh_t       *));

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
        PDM_dmesh_nodal_t          *dmn
)
{
  dmn_to_dm->dmesh_nodal[i_mesh] = dmn;
}


void
PDM_dmesh_nodal_to_dmesh_get_dmesh
(
        PDM_dmesh_nodal_to_dmesh_t  *dmn_to_dm,
  const int                          i_mesh,
        PDM_dmesh_t                **dm
)
{
  *dm = dmn_to_dm->dmesh[i_mesh];
  dmn_to_dm->results_is_getted = PDM_TRUE;
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

  if(( dmn_to_dm->owner == PDM_OWNERSHIP_KEEP ) ||
     ( dmn_to_dm->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !dmn_to_dm->results_is_getted)){

    for(int i_mesh = 0; i_mesh < dmn_to_dm->n_mesh; ++i_mesh) {
      PDM_dmesh_free( (PDM_dmesh_t*)dmn_to_dm->dmesh[i_mesh]);
    }
  }

  free(dmn_to_dm->dmesh_nodal);
  free(dmn_to_dm->dmesh);
}



void
PDM_dmesh_nodal_to_dmesh_compute
(
        PDM_dmesh_nodal_to_dmesh_t                 *dmn_to_dm,
  const PDM_dmesh_nodal_to_dmesh_transform_t        transform_kind,
  const PDM_dmesh_nodal_to_dmesh_translate_group_t  transform_group_kind
)
{


  for(int i_mesh = 0; i_mesh < dmn_to_dm->n_mesh; ++i_mesh) {

    if(dmn_to_dm->dmesh_nodal[i_mesh]->mesh_dimension == 2) {
      assert(transform_kind != PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE);
    }

    switch (transform_kind) {

      case PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE:
        {
          dmn_to_dm->dmesh[i_mesh] = _generate_faces_from_dmesh_nodal(dmn_to_dm->dmesh_nodal[i_mesh]);
        }
        break;

      case PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE:
        {
          dmn_to_dm->dmesh[i_mesh] = _generate_edges_from_dmesh_nodal(dmn_to_dm->dmesh_nodal[i_mesh]);
        }
        break;
    }
  }

  // Boundary management
  for(int i_mesh = 0; i_mesh < dmn_to_dm->n_mesh; ++i_mesh) {
    if(dmn_to_dm->dmesh_nodal[i_mesh]->dgroup_elmt != NULL) {
      switch (transform_group_kind) {
        case PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_NONE:
          break;
        case PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE:
          {
            assert(transform_kind == PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE);
            _translate_element_group_to_faces(dmn_to_dm->dmesh_nodal[i_mesh], dmn_to_dm->dmesh[i_mesh]);
          }
          break;
        case PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE:
          {
            PDM_error (__FILE__, __LINE__, 0, "PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE not implemented \n");
          }
          break;
        case PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_VTX:
          {
            PDM_error (__FILE__, __LINE__, 0, "PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE not implemented \n");
          }
          break;
      }
    }
  }

  // Join management

  // _check_mesh();
  //    --> Sur option
  //    --> Missing boundary faces connected to elements
  //    --> Missing join
  //    --> Other idea
  //
  // Here we have all connectivity for elmt :
  // translate_for_3d_elements = cells so replace face_elmt by face_cell, elmt_face by cell_face ...
  // translate_for_2d_elements = faces etc ...
  // Comment faire avec les const du dmesh ???
  // Quid de l'orientation pour les edges ? edges_cell est orienté ???

  //            |
  //         1  |   2
  //     ------ . -----
  //         3  |   4
  //            |
  //  dedge_cell = (1 2 3 4) =! (1 2 4 3)
  //  dcell_face + dface_edge + dedge_cell

  // face puis edge ---> face --> ???? --> edge_face
  //
  // Face d'abord

  // Il faut eviter de skipper
  // cell_face + face_edge + edge_vtx ---> Il faut aussi déduire le face_vtx -> PDM_deduce_combine_connectivity

  // Si on veut faire un volume de controle autour des edges il faut qu'il soit dans le bon sens non ??
  // En fait il y a l'orientation ET le sens qui importe
  // Comme pour le face_vtx le sens du edge_cell me semble important


  // 1er etape --> Tous les connectivité descendatnes elmnt_***
  // 2eme etape --> dmesh_  --> delmnt_face --> delmt_cell
  //

  // Post-treatment
  // PDM_dmesh_nodal_to_dmesh_transform_to_coherent_dmesh();


}


void
PDM_dmesh_nodal_to_dmesh_transform_to_coherent_dmesh
(
        PDM_dmesh_nodal_to_dmesh_t *dmn_to_dm,
  const int                         extract_dim
)
{
  assert(extract_dim >= 2);

  for(int i_mesh = 0; i_mesh < dmn_to_dm->n_mesh; ++i_mesh) {

    // PDM_dmesh_nodal_t* dmnl = dmn_to_dm->dmesh_nodal[i_mesh];

    // 1. delmt_face --> dcell_face
    // PDM_g_num_t* dcell_face = (PDM_g_num_t * ) malloc( dmnl->dcell_face_idx[dmnl->dn_elmt] ) * sizeof(PDM_g_num_t));

    // Find out in dmesh_nodal where the last element of current proc
    // PDM_g_num_t first_elmt = -1;
    // PDM_g_num_t last_elmt  = -1;
    // int first_section = PDM_binary_search_gap_long(first_elmt, dmnl->section_distribution, dmnl->n_rank+1);
    // int last_section  = PDM_binary_search_gap_long(last_elmt , dmnl->section_distribution, dmnl->n_rank+1);

    // // Ok anyway for std is you pack them, it's simple, but for poly3d the pb is the same !
    // // for(int i_section = 0; i_section < dmnl->n_section_tot; ++i_section) {
    // int idx   = 0;
    // int icell = 0;
    // dcell_face_idx[0] = 0;
    // for(int i_section = first_section; i_section < last_section; ++i_section) {

    //   PDM_section_type_t section_type = dmnl->section_type[i_section];
    //   PDM_g_num_t *section_distrib = NULL; // Faire une foncion du type d'emlt (std/poly2d/poly3d) --> PDM_DMesh_nodal_distrib_section_get

    //   // idem on cherche le premier elmt dans la distribtuion de rank
    //   //
    //   if(section_type == PDM_SECTION_TYPE_STD3D ||
    //      section_type == PDM_SECTION_TYPE_POLY3D )
    //   {
    //     int beg = PDM_MAX(dmnl->section_distribution[i_section  ] - dmnl->elmt_distrib[dmnl->i_rank  ], 0);
    //     int end = PDM_MIN(dmnl->section_distribution[i_section+1], dmnl->elmt_distrib[dmnl->i_rank+1]) - dmnl->section_distribution[i_section  ]; // Impossible d'aller plus loin que la section
    //     // int end = PDM_MIN(dmnl->section_distribution[i_section+1], dmnl->elmt_distrib[dmnl->i_rank+1]) - dmnl->elmt_distrib[i_section  ]; // Impossible d'aller plus loin que la section

    //     // Attention il faut bien l'indice du block elemt
    //     for( int ielmt = beg; ielmt < end; ++ielmt ) {
    //       dcell_face_idx[icell+1] += dcell_face_idx[icell];
    //       for(int iface = delmt_face_idx[ielmt]; iface < delmt_face_idx[ielmt+1]; ++iface ){
    //         dcell_face[idx++] = delmt_face[iface];
    //       }
    //       icell += 1;
    //     }
    //   }


      // for(int ielmt = 0; ielmt < dmnl->dn_elmt; ++ielmt) {

      // }

      // ---------------------------------------
      // 1        10     20       30
      // |  HEXA |  TRI  |  TETRA |

      // ---------------------------------------
      // 1       10       20
      // |  HEXA |  TETRA |

      // face_elmt = [ 1 30 ] --> face_cell = [1 20 ]

      // Si la section est 2D on shift
      // parent_gnum = elmt_to_cell

      // Si 3D on garde ...
      // On construit les cellules suprimée ??? dans un tableau qu'on tri
      // On peut supprimer par range
      // Tout les elements entre x et y doivent être supprimé --> binary_search_gap ?
    // }

    // Optimisation 1 : si tout les elements 2D sont au bout on realloc juste le tableau

    // 2. dface_elmt --> dface_cell
    // Caution : if face_cell is in the bad sens we need to adapt face_vtx (see partitioning with flip)
    // binary_search ???? Pour enlever les cellules pas ???

  }


}
