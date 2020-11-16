
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

  PDM_g_num_t n_face_proc = dn_entity;
  PDM_g_num_t beg_num_abs;

  PDM_MPI_Scan(&n_face_proc, &beg_num_abs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  beg_num_abs -= n_face_proc;

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
_pdm_dmesh_nodal_t *mesh,
int                 n_entity_elt_tot,
PDM_g_num_t        *delmt_entity,
int                *delmt_entity_vtx_idx,
PDM_g_num_t        *delmt_entity_vtx,
int                *dn_entity,
PDM_g_num_t       **entity_distrib,
int               **dentity_vtx_idx,
PDM_g_num_t       **dentity_vtx,
int               **delmt_entity_out_idx,
PDM_g_num_t       **delmt_entity_out,
PDM_g_num_t       **dentity_elmt
)
{
  /*
   * We are now all information flatten - we only need to compute hash_keys for each entitys
   */
  PDM_g_num_t* ln_to_gn = (PDM_g_num_t*) malloc( n_entity_elt_tot * sizeof(PDM_g_num_t));
  PDM_g_num_t key_mod = 4 * mesh->n_vtx_abs;
  // printf("key_mod::%i \n", key_mod);

  // Option des clés
  _compute_keys(n_entity_elt_tot,
                delmt_entity_vtx_idx,
                delmt_entity_vtx,
                ln_to_gn,
                key_mod);

  if(0 == 1) {
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
  for(int i_entity = 0; i_entity < blk_size; ++i_entity) {
    n_max_entity_per_key = PDM_MAX(n_max_entity_per_key, blk_n_entity_per_key[i_entity]);
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
  PDM_g_num_t* loc_entity_vtx_1 = (PDM_g_num_t *) malloc( n_max_vtx            * sizeof(PDM_g_num_t) );
  PDM_g_num_t* loc_entity_vtx_2 = (PDM_g_num_t *) malloc( n_max_vtx            * sizeof(PDM_g_num_t) );
  int*         already_treat    = (int         *) malloc( n_max_entity_per_key * sizeof(int        ) );

  /*
   * Allocate Memory - entity_vtx - entity_elmt
   */
  *dentity_vtx     = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  blk_tot_entity_vtx_size  );
  *dentity_vtx_idx = (int         *) malloc( sizeof(int        ) * (blk_entity_vtx_n_size+1 ));
  *dentity_elmt    = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  blk_entity_elmt_size * 2 );

  PDM_g_num_t *_dentity_vtx     = *dentity_vtx;
  int         *_dentity_vtx_idx = *dentity_vtx_idx;
  PDM_g_num_t *_dentity_elmt    = *dentity_elmt;


  printf("blk_entity_elmt_size::%i\n", blk_entity_elmt_size);

  /*
   * Init global numbering
   */
  int i_abs_entity = 0;
  _dentity_vtx_idx[0] = 0;

  int idx = 0;
  int idx_entity_vtx = 0;
  for(int i_key = 0; i_key < blk_size; ++i_key) {
    // printf("Number of conflicting keys :: %i \n", blk_n_entity_per_key[i_key]);

    int n_conflict_entitys = blk_n_entity_per_key[i_key];

    /* Reset */
    for(int j = 0; j < n_conflict_entitys; ++j) {
      already_treat[j] = -1;
    }

    /* Loop over all entitys in conflict */
    for(int i_entity = 0; i_entity < n_conflict_entitys; ++i_entity) {
      // printf("Number of vtx on entitys %i :: %i with index [%i] \n", i_entity, blk_entity_vtx_n[idx+i_entity], idx+i_entity);

      int n_vtx_entity_1 = blk_entity_vtx_n[idx+i_entity];
      int beg_1        = blk_entity_vtx_idx[idx+i_entity];
      if(already_treat[i_entity] != 1) {

        PDM_g_num_t key_1 = 0;
        for(int j = 0; j < n_vtx_entity_1; ++j) {
          loc_entity_vtx_1[j] = blk_tot_entity_vtx[beg_1+j];
          key_1 += loc_entity_vtx_1[j];
        }
        PDM_quick_sort_long(loc_entity_vtx_1, 0, n_vtx_entity_1-1);

        for(int i_entity_next = i_entity+1; i_entity_next < n_conflict_entitys; ++i_entity_next) {

          int n_vtx_entity_2 = blk_entity_vtx_n[idx+i_entity_next];

          if(n_vtx_entity_1 == n_vtx_entity_2 ) {

            int beg_2 = blk_entity_vtx_idx[idx+i_entity_next];
            PDM_g_num_t key_2 = 0;
            for(int j = 0; j < n_vtx_entity_1; ++j) {
              loc_entity_vtx_2[j] = blk_tot_entity_vtx[beg_2+j];
              key_2 += loc_entity_vtx_2[j];
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
              // printf(" It's a match ! \n");

              _dentity_elmt[2*i_abs_entity  ] = blk_elmt_entity_elmt[idx+i_entity     ];
              _dentity_elmt[2*i_abs_entity+1] = blk_elmt_entity_elmt[idx+i_entity_next];

              for(int i_vtx = 0; i_vtx < n_vtx_entity_1; ++i_vtx) {
                _dentity_vtx[idx_entity_vtx++] = blk_tot_entity_vtx[beg_1+i_vtx];
              }
              _dentity_vtx_idx[i_abs_entity+1] = _dentity_vtx_idx[i_abs_entity] + n_vtx_entity_1;
              i_abs_entity++;

              already_treat[i_entity]      = 1;
              already_treat[i_entity_next] = 1;
            }


          } /* End if same number of vertex */
        } /* End loop next entity */
      } /* End loop already treated */

      /* If a entity is not found a match inside the pool, it's a boundary condition */
      if(already_treat[i_entity] != 1) {

        _dentity_elmt[2*i_abs_entity  ] = blk_elmt_entity_elmt[idx+i_entity];
        _dentity_elmt[2*i_abs_entity+1] = 0;

        for(int i_vtx = 0; i_vtx < n_vtx_entity_1; ++i_vtx) {
          _dentity_vtx[idx_entity_vtx++] = blk_tot_entity_vtx[beg_1+i_vtx];
        }
        _dentity_vtx_idx[i_abs_entity+1] = _dentity_vtx_idx[i_abs_entity] + n_vtx_entity_1;
        i_abs_entity++;

        already_treat[i_entity]      = 1;

      }

    } /* End loop entity in conflict */

    idx += n_conflict_entitys;

  }

  /*
   * Free all unused structure
   */
  free(loc_entity_vtx_1);
  free(loc_entity_vtx_2);
  free(already_treat);
  free(delmt_entity);
  free(delmt_entity_vtx);
  free(ln_to_gn);
  free(blk_tot_entity_vtx);
  free(blk_n_entity_per_key);
  free(blk_entity_vtx_n);
  free(blk_elmt_entity_elmt);

  if( 0 == 1 ){
    printf("i_abs_entity::%i \n", i_abs_entity);
    PDM_log_trace_array_int(_dentity_vtx_idx, i_abs_entity+1                    , "_dentity_vtx_idx:: ");
    PDM_log_trace_array_long(_dentity_vtx   , _dentity_vtx_idx[i_abs_entity], "_dentity_vtx:: ");
    PDM_log_trace_array_long(_dentity_elmt  , 2*i_abs_entity                    , "_dentity_elmt:: ");
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

  *dentity_elmt    = (PDM_g_num_t *) realloc(*dentity_elmt   , _dn_entity * 2 * sizeof(PDM_g_num_t * ));
  _dentity_elmt    = *dentity_elmt;

  /*
   * Generate absolute numerotation of entitys
   */
  *entity_distrib = _make_absolute_entity_numbering(_dn_entity, mesh->pdm_mpi_comm);
  PDM_g_num_t* _entity_distrib = *entity_distrib;

  /*
   * Rebuild elmt entity
   */
  int n_entity_elmt = 2*_dn_entity;

  int*         part_stri_entity_elmt = (int         *) malloc( sizeof(int        ) * n_entity_elmt );
  PDM_g_num_t* dentity_elmt_tmp      = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * n_entity_elmt );
  PDM_g_num_t* ln_to_gn_elem       = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * n_entity_elmt );

  int idx_g = 0;
  assert( _dn_entity == _entity_distrib[mesh->i_rank+1] - _entity_distrib[mesh->i_rank]);

  for (int i_entity = 0; i_entity < _dn_entity; ++i_entity) {

    PDM_g_num_t g_num_entity = (PDM_g_num_t) i_entity + _entity_distrib[mesh->i_rank] + 1;

    dentity_elmt_tmp[idx_g] = _dentity_elmt[2*i_entity];
    ln_to_gn_elem [idx_g] = g_num_entity;
    part_stri_entity_elmt[idx_g] = 1;
    idx_g++;
    if(_dentity_elmt[2*i_entity+1] != 0){
      dentity_elmt_tmp[idx_g] = _dentity_elmt[2*i_entity+1];
      ln_to_gn_elem [idx_g] = g_num_entity;
      part_stri_entity_elmt[idx_g] = 1;
      idx_g++;
    } else {
      dentity_elmt_tmp[idx_g] = _dentity_elmt[2*i_entity];
      ln_to_gn_elem [idx_g] = g_num_entity;
      part_stri_entity_elmt[idx_g] = 0;
      idx_g++;
    }
  }
  n_entity_elmt = idx_g; // Adapt size

  if(0 == 1 ){
    printf("n_entity_elmt::%i\n", n_entity_elmt);
    PDM_log_trace_array_int(part_stri_entity_elmt, n_entity_elmt, "part_stri_entity_elmt:: ");
    PDM_log_trace_array_long(ln_to_gn_elem     , n_entity_elmt, "ln_to_gn_elem:: ");
    PDM_log_trace_array_long(dentity_elmt_tmp    , n_entity_elmt, "dentity_elmt_tmp:: ");
  }

  /*
   *  Use part_to_block with the elmt numbering
   */
  PDM_part_to_block_t *ptb2 = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       &dentity_elmt_tmp,
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
  free(dentity_elmt_tmp     );
  free(ln_to_gn_elem      );

  /*
   * Allcoate
   */
  int* _delmt_entity_out_idx = *delmt_entity_out_idx;
  assert(_delmt_entity_out_idx == NULL);
  _delmt_entity_out_idx = (int * ) malloc( (delmt_tot + 1) * sizeof(int) );

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
_pdm_dmesh_t*
_generate_faces_from_dmesh_nodal
(
  _pdm_dmesh_nodal_t *dmesh_nodal
)
{

  PDM_dmesh_nodal_t* dmn = (PDM_dmesh_nodal_t *) dmesh_nodal;
  PDM_UNUSED(dmesh_nodal);

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

// PDM_g_num_t       **dentity_vtx,
// int               **dentity_vtx_idx,
// PDM_g_num_t       **dentity_elmt
  printf("_generate_query_entitiy \n");
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
                                 &dmesh_nodal->_dface_cell);
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
                                 &dmesh_nodal->_dface_cell);

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
        PDM_dmesh_nodal_t          *dmn
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

      case PDM_dmesh_nodal_tO_DMESH_TRANSFORM_TO_FACE:
        {
          _dmn_to_dm->dmesh[i_mesh] = _generate_faces_from_dmesh_nodal(_dmn_to_dm->dmesh_nodal[i_mesh]);
        }
        break;

      case PDM_dmesh_nodal_tO_DMESH_TRANSFORM_TO_EDGE:
        {
          _dmn_to_dm->dmesh[i_mesh] = _generate_edges_from_dmesh_nodal(_dmn_to_dm->dmesh_nodal[i_mesh]);
        }
        break;
    }
  }

  // Boundary management

  // Join management


}
