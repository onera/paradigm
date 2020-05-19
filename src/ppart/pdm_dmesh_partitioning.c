
/*============================================================================
 * Parallel partitioning
 *============================================================================*/

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
#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_timer.h"
#include "pdm_mpi.h"
#include "pdm_mpi_ext_dependencies.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"
#include "pdm_handles.h"

#include "pdm_dmesh_partitioning.h"
#include "pdm_para_graph_dual.h"

/*----------------------------------------------------------------------------
 *  Optional headers
 *----------------------------------------------------------------------------*/

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

/**
 * \def _PDM_part_MIN(a,b)
 * Computes the minimum of \a x and \a y.
 *
 */

#define _PDM_part_MIN(a,b) ((a) > (b) ? (b) : (a))

/**
 * \def _PDM_part_MAX(a,b)
 * Computes the maximum of \a x and \a y.
 *
 */

#define _PDM_part_MAX(a,b) ((a) < (b) ? (b) : (a))

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _dmesh_partitioning_t
 * \brief  Define a global numberring
 *
 */
typedef struct {


} _dmesh_partitioning_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_dmps   = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Return _dmesh_partitioning_t object from it identifier
 *
 * \param [in]   dmpartitioning_id        ppart identifier
 *
 */
static _dmesh_partitioning_t *
_get_from_id
(
 int  dmpartitioning_id
)
{
  _dmesh_partitioning_t *dmesh_partitioning = (_dmesh_partitioning_t *) PDM_Handles_get (_dmps, dmpartitioning_id);

  if (dmesh_partitioning == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "pdm_dmesh_partitioning error : Bad _dmesh_partitioning identifier\n");
    exit(1);
  }

  return dmesh_partitioning;
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *  \brief Setup cell_ln_to_gn
 */
int
PDM_generate_part_cell_ln_to_gn
(
 const PDM_MPI_Comm    comm,
 PDM_g_num_t          *part_distribution,
 PDM_g_num_t          *cell_distribution,
 int                  *cell_part,
 int                 **n_elmts,
 PDM_g_num_t        ***pcell_ln_to_gn
)
{
  printf("PDM_generate_part_cell_ln_to_gn\n");
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_part = part_distribution[i_rank+1] -  part_distribution[i_rank];
  int dn_cell = cell_distribution[i_rank+1] -  cell_distribution[i_rank];

  /*
   * On recréer un tableau pourtant dans scotch et metis le part est en int64 ...
   *     Il faudrait changer le ext_dependancies ..
   */
  PDM_g_num_t* dpart_ln_to_gn = (PDM_g_num_t * ) malloc( sizeof(PDM_g_num_t) * dn_cell );
  for(int i = 0; i < dn_cell; ++i){
    dpart_ln_to_gn[i] = (PDM_g_num_t) cell_part[i] + 1;
  }

  PDM_g_num_t* part_distribution_ptb = (PDM_g_num_t * ) malloc( sizeof(PDM_g_num_t) * (n_rank+1) );
  for(int i = 0; i < n_rank+1; ++i){
    part_distribution_ptb[i] = part_distribution[i] - 1;
  }

  /*
   * The tricks is to use part_to_block with the part_distribution to have n_part stride of data
   */
  PDM_part_to_block_t *ptb_partition =
   PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                             PDM_PART_TO_BLOCK_POST_MERGE,
                             1.,
                             &dpart_ln_to_gn,
                              part_distribution_ptb,
                             &dn_cell,
                             1,
                             comm);

  const int n_part_block = PDM_part_to_block_n_elt_block_get (ptb_partition);

  /*
   * Generate cell_ln_to_gn
   */
  int* dcell_stri = (int *) malloc( sizeof(int) * dn_cell );
  PDM_g_num_t* dcell_ln_to_gn = (PDM_g_num_t * ) malloc( sizeof(PDM_g_num_t) * dn_cell );

  PDM_g_num_t shift_g = cell_distribution[i_rank];
  for(int i = 0; i < dn_cell; ++i){
    dcell_stri[i] = 1;
    dcell_ln_to_gn[i] = (PDM_g_num_t) shift_g + i;
  }

  /*
   * Exchange
   */
  int*         pcell_stri         = NULL;
  PDM_g_num_t* pcell_ln_to_gn_tmp = NULL;

  PDM_part_to_block_exch (ptb_partition,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          1,
                          &dcell_stri,
                (void **) &dcell_ln_to_gn,
                          &pcell_stri,
                (void **) &pcell_ln_to_gn_tmp);

  /*
   * Free
   */
  free(dcell_stri);
  free(dcell_ln_to_gn);

  /*
   *  Post-traitement
   */
  *n_elmts        = (int *         ) malloc( sizeof(int          ) * n_part_block);
  *pcell_ln_to_gn = (PDM_g_num_t **) malloc( sizeof(PDM_g_num_t *) * n_part_block);

  /*
   *  Shortcut
   */
  int*          _n_elmts       = (int *) *n_elmts;
  PDM_g_num_t** _cell_ln_to_gn = (PDM_g_num_t ** ) *pcell_ln_to_gn;

  int idx_part = 0;
  for(int i_part = 0; i_part < n_part_block; ++i_part){

    /* Suralloc */
    _cell_ln_to_gn[i_part] = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * pcell_stri[i_part]);
    PDM_g_num_t* _part_ln_to_gn = (PDM_g_num_t *) _cell_ln_to_gn[i_part];

    PDM_sort_long(&pcell_ln_to_gn_tmp[idx_part], NULL, pcell_stri[i_part]);

    /* Compress */
    int idx_new_cell = 0;
    PDM_g_num_t last_value = 0; // 0 est une bonne valeur car les numero absolu son [N, -1] - |1, N] :p !
    for(int i_cell = 0; i_cell < pcell_stri[i_part]; ++i_cell){
      if(last_value != pcell_ln_to_gn_tmp[idx_part+i_cell]){
        last_value = pcell_ln_to_gn_tmp[idx_part+i_cell];
        _part_ln_to_gn[idx_new_cell++] = last_value;
      }
    }

    /* For cell the size is the same - No compression */
    assert(idx_new_cell ==  pcell_stri[i_part]);
    _n_elmts[i_part] = idx_new_cell;

    /* Realloc */
    idx_part += pcell_stri[i_part];

    /*
     * Panic verbose
     */
    if(0 == 1){
      printf(" _part_ln_to_gn = ");
      for(int i_data = 0; i_data < idx_new_cell; ++i_data){
        printf("%d ", _part_ln_to_gn[i_data]);
      }
      printf("\n");
    }

  }

  PDM_part_to_block_free (ptb_partition);

  free(pcell_stri);
  free(pcell_ln_to_gn_tmp);
  free(dpart_ln_to_gn);
  free(part_distribution_ptb);

  return n_part_block;

}


/**
 *  \brief Setup cell_ln_to_gn
 */
void
PDM_generate_part_face_group_ln_to_gn
(
 const PDM_MPI_Comm    comm,
 PDM_g_num_t          *face_distribution,
 int                  *dface_group_idx,
 PDM_g_num_t          *dface_group,
 int                   n_part,
 int                   n_face_group,
 int                  *n_faces,
 PDM_g_num_t         **pface_ln_to_gn,
 PDM_g_num_t        ***pface_group_ln_to_gn,
 int                ***pface_group,
 int                ***pface_group_idx
)
{
  printf("PDM_generate_part_face_group_ln_to_gn\n");
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_face = face_distribution[i_rank+1] -  face_distribution[i_rank];
  printf("dn_face::%d\n", dn_face);

  PDM_g_num_t* face_distribution_ptb = (PDM_g_num_t * ) malloc( sizeof(PDM_g_num_t) * (n_rank+1) );
  for(int i = 0; i < n_rank+1; ++i){
    face_distribution_ptb[i] = face_distribution[i] - 1;
  }

  /*
   * Compute the dstribution_face_group - Mandatory to have the link between rank and group
   */
  PDM_g_num_t** face_group_distribution = (PDM_g_num_t **) malloc( n_face_group * sizeof(PDM_g_num_t *));

  for(int i_group = 0; i_group < n_face_group; ++i_group) {
    int dn_face_group = dface_group_idx[i_group+1] - dface_group_idx[i_group];
    face_group_distribution[i_group] = PDM_compute_entity_distribution_long(comm, dn_face_group);
  }

  /*
   * Create exchange protocol
   *    - dface_group contains the absolute number of faces that contains the boundary conditions -> A kind of ln_to_gn
   */
  // int n_face_group = -1;
  int n_total_face_group = dface_group_idx[n_face_group];
  PDM_part_to_block_t *ptb_fg =
   PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                             PDM_PART_TO_BLOCK_POST_MERGE,
                             1.,
                             &dface_group,
                              face_distribution_ptb,
                             &n_total_face_group,
                             1,
                             comm);

  /*
   * Prepare send
   *    Reminder :face_group have double indirection
   *       - For each group we have a list of face
   *       - All rank have a distribution of each group
   * part_data should contains : ( position in distributed array + group_id )
   */
  int*         part_stri = (int         *) malloc( sizeof(int        )     * n_total_face_group );
  PDM_g_num_t* part_data = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * 2 * n_total_face_group );

  int idx_data = 0;
  int idx_stri = 0;
  for(int i_group = 0; i_group < n_face_group; ++i_group) {
    for(int i_face = dface_group_idx[i_group]; i_face < dface_group_idx[i_group+1]; ++i_face ) {
      part_data[idx_data++] = i_group;
      part_data[idx_data++] = i_face + face_group_distribution[i_group][i_rank]; /* Numero dans le tableau distribue */
      part_stri[idx_stri++] = 2;
    }
  }

  /*
   * Exchange
   */
  const int n_face_block = PDM_part_to_block_n_elt_block_get(ptb_fg);
  PDM_g_num_t* blk_gnum  = PDM_part_to_block_block_gnum_get(ptb_fg);

  int*         blk_stri = NULL;
  PDM_g_num_t* blk_data = NULL;
  PDM_part_to_block_exch (ptb_fg,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          1,
                          &part_stri,
                (void **) &part_data,
                          &blk_stri,
                (void **) &blk_data);

  /*
   * Post-treatement
   */
  int idx_face = 0;
  for(int i_face = 0; i_face < n_face_block; ++i_face) {
    printf("[%d] - %d | %d --> ", i_face, blk_stri[i_face], blk_gnum[i_face]);
    for(int i_data = 0; i_data < blk_stri[i_face]; ++i_data) {
      printf("%d ", blk_data[idx_face++]);
    }
    printf("\n");
  }

  /*
   *  Maintenant on peut faire un block_to_part sur le face_ln_to_gn (le vrai) pour recuperer pour chaque partition
   *    le face_group
   */




  free(face_distribution_ptb);
  free(part_data);
  free(part_stri);
  free(blk_stri);
  free(blk_data);
  PDM_part_to_block_free (ptb_fg);
  for(int i_group = 0; i_group < n_face_group; ++i_group) {
    free(face_group_distribution[i_group]);
  }
  free(face_group_distribution);
}





/**
 *  \brief Setup cell_ln_to_gn
 */
void
PDM_generate_part_entity_ln_to_gn
(
 const PDM_MPI_Comm    comm,
 PDM_g_num_t          *part_distribution,
 PDM_g_num_t          *cell_distribution,
 int                  *dcell_face_idx,
 PDM_g_num_t          *dcell_face,
 int                   n_part,
 int                  *n_elmts,
 PDM_g_num_t         **pcell_ln_to_gn,
 int                 **n_faces,
 PDM_g_num_t        ***pface_ln_to_gn,
 int                ***pcell_face_idx,
 int                ***pcell_face
)
{
  printf("PDM_generate_part_entity_ln_to_gn\n");
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_part = part_distribution[i_rank+1] -  part_distribution[i_rank];
  int dn_cell = cell_distribution[i_rank+1] -  cell_distribution[i_rank];

  PDM_g_num_t* cell_distribution_ptb = (PDM_g_num_t * ) malloc( sizeof(PDM_g_num_t) * (n_rank+1) );
  for(int i = 0; i < n_rank+1; ++i){
    cell_distribution_ptb[i] = cell_distribution[i] - 1;
  }

  /*
   * Prepare exchange protocol
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(cell_distribution_ptb,
                               (const PDM_g_num_t **) pcell_ln_to_gn,
                                                      n_elmts,
                                                      n_part,
                                                      comm);

  /*
   * Prepare data
   */
  int* blk_stri = (int *) malloc( sizeof(int) * dn_cell);
  for(int i_cell = 0; i_cell < dn_cell; ++i_cell){
    blk_stri[i_cell] = dcell_face_idx[i_cell+1] - dcell_face_idx[i_cell];
  }

  /*
   * Exchange
   */
  int** cell_stri;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          blk_stri,
             (void *  )   dcell_face,
             (int *** )  &cell_stri,
             (void ***) &*pcell_face);

  free(blk_stri);

  /*
   * Panic verbose
   */
  if(0 == 1){
    for(int i_part = 0; i_part < n_part; ++i_part){
      int idx_data = 0;
      printf("[%d] cell_face:: \n", i_part);
      for(int i_cell = 0; i_cell < n_elmts[i_part]; ++i_cell) {
        printf("[%d] --> ", cell_stri[i_part][i_cell]);
        for(int i_data = 0; i_data < cell_stri[i_part][i_cell]; ++i_data ){
          printf("%d ", (*pcell_face)[i_part][idx_data++] );
        }
        printf("\n");
      }
      printf("\n");
    }
  }


  /*
   * Post-treatment - Caution the recv connectivity can be negative
   */
  *pface_ln_to_gn = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *) );
  *pcell_face_idx = (int         **) malloc( n_part * sizeof(int         *) );
  *n_faces        = (int          *) malloc( n_part * sizeof(int          ) );

  /* Shot cut */
  PDM_g_num_t** _face_ln_to_gn = *pface_ln_to_gn;
  int** _cell_face_idx         = *pcell_face_idx;
  int*  _n_faces               = *n_faces;

  for(int i_part = 0; i_part < n_part; ++i_part){

    /*
     *  First loop to count and setup part_strid_idx
     */
    _cell_face_idx[i_part] = (int *) malloc( (n_elmts[i_part] + 1) * sizeof(int) );
    _cell_face_idx[i_part][0] = 0;
    for(int i_cell = 0; i_cell < n_elmts[i_part]; ++i_cell) {
      _cell_face_idx[i_part][i_cell+1] = _cell_face_idx[i_part][i_cell] + cell_stri[i_part][i_cell];
    }

    /*
     * Save array
     */
    _face_ln_to_gn[i_part] = (PDM_g_num_t *) malloc( _cell_face_idx[i_part][n_elmts[i_part]] * sizeof(PDM_g_num_t));
    PDM_g_num_t* _pface_ln_to_gn = (PDM_g_num_t *) _face_ln_to_gn[i_part];

    int idx_data = 0;
    for(int i_cell = 0; i_cell < n_elmts[i_part]; ++i_cell) {
      for(int i_data = 0; i_data < cell_stri[i_part][i_cell]; ++i_data ){
        _pface_ln_to_gn[idx_data++] = PDM_ABS((*pcell_face)[i_part][idx_data]);
      }
    }

    /*
     * Deduce ln_to_gn
     */
    printf("Sort data between : 0 and %d \n", idx_data);
    int n_elmt_sort = PDM_inpace_unique_long(_pface_ln_to_gn , 0, idx_data);
    _n_faces[i_part] = n_elmt_sort;

    printf("n_elmt_sort::%d\n", n_elmt_sort);
    printf("_pface_ln_to_gn::");
    for(int i = 0; i < idx_data; ++i){
      printf("%d ", _pface_ln_to_gn[i]);
    }
    printf("\n");

    /*
     * Realloc
     */
    _face_ln_to_gn[i_part] = (PDM_g_num_t *) realloc(_face_ln_to_gn[i_part], n_elmt_sort * sizeof(PDM_g_num_t) );
    _pface_ln_to_gn = _face_ln_to_gn[i_part];

    /*
     *  We need to regenerate the connectivity and pass it in local numbering
     */
    idx_data = 0;
    for(int i_cell = 0; i_cell < n_elmts[i_part]; ++i_cell) {
      for(int i_data = 0; i_data < cell_stri[i_part][i_cell]; ++i_data ){
        int         g_sgn  = 1;
        PDM_g_num_t g_elmt = PDM_ABS((*pcell_face)[i_part][idx_data]);
        int l_elmt         = PDM_binary_search_long(g_elmt, _pface_ln_to_gn, n_elmt_sort); /* In [0, n_elmt_sort-1] */

        /* Overwrite the pcell_face with local numbering and reput sign on it */
        (*pcell_face)[i_part][idx_data++] = (l_elmt + 1) * g_sgn ;
      }
    }
  }


  /*
   * Panic verbose
   */
  if(0 == 1){
    for(int i_part = 0; i_part < n_part; ++i_part){
      int idx_data = 0;
      printf("[%d] cell_face:: \n", i_part);
      for(int i_cell = 0; i_cell < n_elmts[i_part]; ++i_cell) {
        printf("[%d] --> ", cell_stri[i_part][i_cell]);
        for(int i_data = 0; i_data < cell_stri[i_part][i_cell]; ++i_data ){
          printf("%d ", (*pcell_face)[i_part][idx_data++] );
        }
        printf("\n");
      }
      printf("\n");
    }
  }

  /*
   * Free
   */
  PDM_block_to_part_free(btp);
  free(cell_distribution_ptb);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(cell_stri[i_part]);
  }
  free(cell_stri);

}



/**
 *  \brief Setup cell_ln_to_gn
 */
// void
// PDM_generate_part_entity_ln_to_gn_old
// (
//  const PDM_MPI_Comm    comm,
//  PDM_g_num_t          *part_distribution,
//  PDM_g_num_t          *cell_distribution,
//  int                  *dcell_face_idx,
//  PDM_g_num_t          *dcell_face,
//  PDM_g_num_t         **pcell_ln_to_gn,
//  int                ***pcell_face_idx,
//  int                ***pcell_face
// )
// {
//   printf("PDM_generate_part_entity_ln_to_gn\n");
//   int i_rank;
//   int n_rank;

//   PDM_MPI_Comm_rank(comm, &i_rank);
//   PDM_MPI_Comm_size(comm, &n_rank);

//   int dn_part = part_distribution[i_rank+1] -  part_distribution[i_rank];
//   int dn_cell = cell_distribution[i_rank+1] -  cell_distribution[i_rank];

//   /*
//    * On recréer un tableau pourtant dans scotch et metis le part est en int64 ...
//    *     Il faudrait changer le ext_dependancies ..
//    */
//   PDM_g_num_t* dpart_ln_to_gn = (PDM_g_num_t * ) malloc( sizeof(PDM_g_num_t) * dn_cell );
//   for(int i = 0; i < dn_cell; ++i){
//     dpart_ln_to_gn[i] = (PDM_g_num_t) cell_part[i] + 1;
//   }

//   PDM_g_num_t* part_distribution_ptb = (PDM_g_num_t * ) malloc( sizeof(PDM_g_num_t) * (n_rank+1) );
//   for(int i = 0; i < n_rank+1; ++i){
//     part_distribution_ptb[i] = part_distribution[i] - 1;
//   }

//   /*
//    * The tricks is to use part_to_block with the part_distribution to have n_part stride of data
//    */
//   PDM_part_to_block_t *ptb_partition =
//    PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
//                              PDM_PART_TO_BLOCK_POST_MERGE,
//                              1.,
//                              &dpart_ln_to_gn,
//                               part_distribution_ptb,
//                              &dn_cell,
//                              1,
//                              comm);

//   const int n_part_block = PDM_part_to_block_n_elt_block_get (ptb_partition);

//   /*
//    * cell_ln_to_gn
//    */
//   int* dcell_stri = (int *) malloc( sizeof(int) * dn_cell );
//   PDM_g_num_t* dcell_ln_to_gn = (PDM_g_num_t * ) malloc( sizeof(PDM_g_num_t) * dn_cell );

//   PDM_g_num_t shift_g = cell_distribution[i_rank];
//   for(int i = 0; i < dn_cell; ++i){
//     dcell_stri[i] = 1;
//     dcell_ln_to_gn[i] = (PDM_g_num_t) shift_g + i;
//   }
//   int* pcell_stri = NULL;
//   PDM_g_num_t* pcell_ln_to_gn_tmp = NULL;
//   PDM_part_to_block_exch (ptb_partition,
//                           sizeof(PDM_g_num_t),
//                           PDM_STRIDE_VAR,
//                           1,
//                           &dcell_stri,
//                 (void **) &dcell_ln_to_gn,
//                           &pcell_stri,
//                 (void **) &pcell_ln_to_gn_tmp);


//   free(dcell_stri);

//   printf("n_part_block::%d\n",n_part_block );

//   /*
//    * Panic verbose
//    */
//   if(0 == 1){
//     int idx_block = 0;
//     for(int i = 0; i < n_part_block; ++i){
//       printf(" pcell_stri = %d ---> ", pcell_stri[i]);
//       for(int i_data = 0; i_data < pcell_stri[i]; ++i_data){
//         printf("%d ", pcell_ln_to_gn_tmp[idx_block]);
//         idx_block++;
//       }
//       printf("\n");
//     }
//   }


//    /*
//     * En changeant un poil le part_to_block on eviterai un tableau temporaire
//     */
//   int* dcell_face_n = (int * ) malloc( dn_cell * sizeof(int));
//   for(int i = 0; i < dn_cell; ++i){
//     dcell_face_n[i] = dcell_face_idx[i+1] - dcell_face_idx[i];
//   }

//   int* cell_face_n = NULL;
//   PDM_g_num_t* cell_face   = NULL;
//   printf("PDM_part_to_block_exch \n");
//   PDM_part_to_block_exch (ptb_partition,
//                           sizeof(PDM_g_num_t),
//                           PDM_STRIDE_VAR,
//                           1,
//                           &dcell_face_n,
//                 (void **) &dcell_face,
//                           &cell_face_n,
//                 (void **) &cell_face);
//   printf("PDM_part_to_block_exch end \n");

//   /*
//    * Panic verbose
//    */
//   // if(1 == 1){
//   //   int idx_block = 0;
//   //   for(int i = 0; i < n_part_block; ++i){
//   //     printf(" pcell_stri = %d ---> ", cell_face_n[i]);
//   //     for(int i_data = 0; i_data < cell_face_n[i]; ++i_data){
//   //       printf("%d ", cell_face[idx_block]);
//   //       idx_block++;
//   //     }
//   //     printf("\n");
//   //   }
//   // }

//   /*
//    * Sort receving array
//    */
//   int idx_block = 0;
//   for(int i = 0; i < n_part_block; ++i){
//     printf(" Tri %d \n", cell_face_n[i]);
//     PDM_sort_long(&cell_face[idx_block] , NULL, cell_face_n[i]);
//     idx_block += cell_face_n[i];
//   }

//   /*
//    * Panic verbose
//    */
//   if(1 == 1){
//     int idx_block = 0;
//     for(int i = 0; i < n_part_block; ++i){
//       printf(" pcell_stri = %d ---> ", cell_face_n[i]);
//       for(int i_data = 0; i_data < cell_face_n[i]; ++i_data){
//         printf("%d ", cell_face[idx_block]);
//         idx_block++;
//       }
//       printf("\n");
//     }
//   }

//   /*
//    *  In the sort cell_face we have multiple same global numbering face
//    *    We need to unique this list and create the face numbering for each part
//    *    Attention au signe si cell_face signé
//    *    La fonction pourrait en meme temps repasser en numero local et en inplace
//    */
//   PDM_g_num_t** part_face_ln_to_gn = (PDM_g_num_t **) malloc( sizeof(PDM_g_num_t *) * n_part_block);
//   idx_block = 0;
//   for(int i_part = 0; i_part < n_part_block; ++i_part){

//     /* Suralloc */
//     part_face_ln_to_gn[i_part] = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * cell_face_n[i_part]);

//     PDM_sort_long(&cell_face[idx_block], NULL, cell_face_n[i_part]);

//     /* Compress */
//     int idx_new_face = 0;
//     // abort();
//     PDM_g_num_t last_value = 0; // 0 est une bonne valeur car les numero absolu son [N, -1] - |1, N] :p !
//     for(int i_face = 0; i_face < cell_face_n[i_part]; ++i_face){
//       if(last_value != cell_face[idx_block+i_face]){
//         last_value = cell_face[idx_block+i_face];
//         part_face_ln_to_gn[i_part][idx_new_face++] = last_value;
//       }
//     }

//     printf(" part_face_ln_to_gn = ");
//     for(int i_data = 0; i_data < idx_new_face; ++i_data){
//       printf("%d ", part_face_ln_to_gn[i_part][i_data]);
//     }
//     printf("\n");

//     /* Realloc */
//     idx_block += cell_face_n[i_part];

//   }


//   /*
//    * Si besoin du cell_face on doit faire un echange pour connaitre le nombre d'element par faces
//    *    --> Pas besoin du cell_face tout le temps non ?
//    */
//   // PDM_sort_long();


//   PDM_part_to_block_free (ptb_partition);

//   free(dpart_ln_to_gn);
//   free(dcell_face_n);
//   free(part_distribution_ptb);
//   for(int i_part = 0; i_part < n_part_block; ++i_part){
//     free(part_face_ln_to_gn[i_part]);
//   }
//   free(part_face_ln_to_gn);
// }


/**
 *
 * \brief Return _dmesh_partitioning_t object from it identifier
 *
 * \param [in]   dmpartitioning_id        ppart identifier
 *
 */
int
PDM_dmesh_partitioning_create
(
 const PDM_MPI_Comm              comm,
 const PDM_partitioning_method_t split_method
)
{
  /*
   * Search a gnum_from_hash_values free id
   */
  if (_dmps == NULL) {
    _dmps = PDM_Handles_create (4);
  }

  _dmesh_partitioning_t *_dmp = (_dmesh_partitioning_t *) malloc(sizeof(_dmesh_partitioning_t));
  int id = PDM_Handles_store (_dmps, _dmp);

  return id;
}


/**
 *
 * \brief Return _dmesh_partitioning_t object from it identifier
 *
 * \param [in]   dmpartitioning_id        ppart identifier
 *
 */
void
PDM_dmesh_partitioning_set_from_dmesh
(
 const int dmesh_partitioning_id,
 const int dmesh_id
)
{
  // _dmesh_partitioning_t* _dmp = _get_from_id(dmesh_partitioning_id);
  // ---> Depuis l'interface du dmesh
  // Recopie
  // PDM_dmesh_dims_get(dmesh_id, &)
}

/**
 *
 * \brief Free
 *
 * \param [in]   id           Identifier
 *
 */
void
PDM_dmesh_partitioning_free
(
 const int id
)
{
  _dmesh_partitioning_t* _dmp = _get_from_id(id);

  free (_dmp);

  PDM_Handles_handle_free (_dmps, id, PDM_FALSE);

  const int n_dmpartitioning = PDM_Handles_n_get (_dmps);

  if (n_dmpartitioning == 0) {
    _dmps = PDM_Handles_free (_dmps);
  }
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
