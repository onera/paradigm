
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
#include "pdm_unique.h"
#include "pdm_binary_search.h"
#include "pdm_handles.h"

#include "pdm_partitioning_algorithm.h"
#include "pdm_distrib.h"
#include "pdm_order.h"
// #include "pdm_para_graph_dual.h"

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

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/


/*============================================================================
 * Private function definitions
 *============================================================================*/

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

  // int dn_part = part_distribution[i_rank+1] -  part_distribution[i_rank];
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
    face_group_distribution[i_group] = PDM_compute_entity_distribution(comm, dn_face_group);
  }

  /*
   * Create exchange protocol
   *    - dface_group contains the absolute number of faces that contains the boundary conditions -> A kind of ln_to_gn
   */
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
  if(1 == 0){
    int idx_face_d = 0;
    for(int i_face = 0; i_face < n_face_block; ++i_face) {
      printf("[%d] - %d | %d --> ", i_face, blk_stri[i_face], blk_gnum[i_face]);
      for(int i_data = 0; i_data < blk_stri[i_face]; ++i_data) {
        printf("%d ", blk_data[idx_face_d++]);
      }
      printf("\n");
    }
    printf("n_face_block::%d\n", n_face_block);
    printf("idx_face_d    ::%d\n", idx_face_d/2);
  }

  /*
   * No choice we need to rebuild a proper blk_stri (but not blk_data)
   */
  int* blk_stri_full = (int *) malloc( dn_face * sizeof(int));
  for(int i = 0; i < dn_face; ++i){
    blk_stri_full[i] = 0;
  }

  /* On remet la stri au bonne endroit */
  // int idx_face = 0;
  for(int i_face = 0; i_face < n_face_block; ++i_face) {
    int l_elmt = blk_gnum[i_face] - face_distribution[i_rank];
    // printf("[%d] --> %d \n", i_face, l_elmt);
    blk_stri_full[l_elmt] = blk_stri[i_face];
  }
  free(blk_stri);



  /*
   *  Maintenant on peut faire un block_to_part sur le face_ln_to_gn (le vrai) pour recuperer pour chaque partition
   *    le face_group
   */
  if( 0 == 1){
    for(int i_part = 0; i_part < n_part; ++i_part) {
      printf("[%d] pface_ln_to_gn:: ", i_part);
      for(int i_face = 0; i_face < n_faces[i_part]; ++i_face){
        printf("%d ", pface_ln_to_gn[i_part][i_face]);
      }
      printf("\n");
    }
  }

  /*
   * Prepare exchange protocol
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(face_distribution_ptb,
                               (const PDM_g_num_t **) pface_ln_to_gn,
                                                      n_faces,
                                                      n_part,
                                                      comm);

  /*
   * Exchange
   */
  int**         part_face_group_stri;
  PDM_g_num_t** part_face_group_data;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          blk_stri_full,
             (void *  )   blk_data,
             (int  ***)  &part_face_group_stri,
             (void ***)  &part_face_group_data);

  /*
   * Post-Treatment
   */
  if(0 == 1){
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int idx_data_fg = 0;
      printf("[%d] part_face_group_data --> ", i_part);
      for(int i_face = 0; i_face < n_faces[i_part]; ++i_face){
        printf(" \t [%d] ", i_face);
        for(int i_data = 0; i_data < part_face_group_stri[i_part][i_face]; ++i_data) {
          printf("%d ", part_face_group_data[i_part][idx_data_fg++]);
        }
        printf("\n");
      }
    }
  }

  /*
   * Donc sur chaque partition on a dans la numerotation des faces la liste des frontières
   * On doit maintenant recréer les tableaux identiques à pdm_part
   */
  *pface_group_ln_to_gn = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *) );
  *pface_group          = (int         **) malloc( n_part * sizeof(int         *) );
  *pface_group_idx      = (int         **) malloc( n_part * sizeof(int         *) );

  /* Shortcut */
  PDM_g_num_t** _face_group_ln_to_gn = *pface_group_ln_to_gn;
  int**         _face_group          = *pface_group;
  int**         _face_group_idx      = *pface_group_idx;

  int* count_fg = (int *) malloc( (n_face_group + 1) * sizeof(int));
  for(int i_part = 0; i_part < n_part; ++i_part) {

    /*
     * All partition allocate is own array of size n_face_group
     */
    _face_group_idx[i_part] = (int *) malloc( (n_face_group + 1) * sizeof(int) );
    for(int i_group = 0; i_group < n_face_group+1; ++i_group){
      _face_group_idx[i_part][i_group] = 0;
      count_fg[i_group] = 0;
    }

    /*
     * First step to count
     */
    int idx_face_bnd = 0;
    for(int i_face = 0; i_face < n_faces[i_part]; ++i_face){
      // TODO : multiple group per face
      if(part_face_group_stri[i_part][i_face] > 0 ){
        PDM_g_num_t i_group = part_face_group_data[i_part][idx_face_bnd];
        _face_group_idx[i_part][i_group+1] += 1;
        idx_face_bnd += 2;
      }
    }

    /*
     *  Deduce idx and allocate
     */
    for(int i_group = 0; i_group < n_face_group; ++i_group){
      _face_group_idx[i_part][i_group+1] += _face_group_idx[i_part][i_group];
    }
    int n_face_bnd = _face_group_idx[i_part][n_face_group];

    _face_group_ln_to_gn[i_part] = (PDM_g_num_t * ) malloc( n_face_bnd * sizeof(PDM_g_num_t));
    _face_group[i_part]          = (int         * ) malloc( n_face_bnd * sizeof(int        ));

    /*
     * Panic verbose
     */
    if( 0 == 1){
      printf("n_face_bnd::%d\n", n_face_bnd );
      printf("[%d] _face_group_idx::\n", i_part );
      for(int i_group = 0; i_group < n_face_group+1; ++i_group){
        printf("%d ", _face_group_idx[i_part][i_group]);
      }
      printf("\n");
    }

    /*
     * Fill data
     */
    idx_face_bnd = 0;
    for(int i_face = 0; i_face < n_faces[i_part]; ++i_face){
      if(part_face_group_stri[i_part][i_face] > 0 ){
        PDM_g_num_t i_group = part_face_group_data[i_part][idx_face_bnd  ];
        PDM_g_num_t g_face  = part_face_group_data[i_part][idx_face_bnd+1];

        int idx = _face_group_idx[i_part][i_group] + count_fg[i_group];

        _face_group         [i_part][idx] = i_face+1;
        _face_group_ln_to_gn[i_part][idx] = g_face;

        count_fg[i_group] += 1;
        idx_face_bnd += 2;
      }
    }

  }


  /*
   * Panic verbose
   */
  if(0 == 1 ){
    for(int i_part = 0; i_part < n_part; ++i_part) {
      printf("[%d] Boundary \n", i_part);
      for(int i_group = 0; i_group < n_face_group; ++i_group){
        printf("\t [%d]_face_group::", i_group);
        for(int idx = _face_group_idx[i_part][i_group]; idx < _face_group_idx[i_part][i_group+1]; ++idx){
          printf("%d ", _face_group[i_part][idx]);
        }
        printf("\n");
        printf("\t [%d]_face_group_ln_to_gn::", i_group);
        for(int idx = _face_group_idx[i_part][i_group]; idx < _face_group_idx[i_part][i_group+1]; ++idx){
          printf("%d ", _face_group_ln_to_gn[i_part][idx]);
        }
        printf("\n");
      }
    }
  }

  free(count_fg);
  free(blk_stri_full);
  free(face_distribution_ptb);
  free(part_data);
  free(part_stri);
  free(blk_data);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(part_face_group_stri[i_part]);
    free(part_face_group_data[i_part]);
  }
  free(part_face_group_stri);
  free(part_face_group_data);
  PDM_part_to_block_free (ptb_fg);
  PDM_block_to_part_free(btp);
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

  // int dn_part = part_distribution[i_rank+1] -  part_distribution[i_rank];
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
  int**         cell_stri;
  PDM_g_num_t** pcell_face_tmp; /* We keep it in double precision because it contains global numbering */
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          blk_stri,
             (void *  )   dcell_face,
             (int  ***)  &cell_stri,
             (void ***)  &pcell_face_tmp);

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
          printf("%d ", pcell_face_tmp[i_part][idx_data++] );
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
  *pcell_face     = (int         **) malloc( n_part * sizeof(int         *) );
  *pcell_face_idx = (int         **) malloc( n_part * sizeof(int         *) );
  *n_faces        = (int          *) malloc( n_part * sizeof(int          ) );

  /* Shortcut */
  PDM_g_num_t** _face_ln_to_gn = *pface_ln_to_gn;
  int** _cell_face             = *pcell_face;
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
    _cell_face[i_part]     = (int *        ) malloc( _cell_face_idx[i_part][n_elmts[i_part]] * sizeof(int        ));

    PDM_g_num_t* _pface_ln_to_gn = _face_ln_to_gn[i_part];
    int*         _pcell_face     = _cell_face[i_part];

    int idx_data = 0;
    for(int i_cell = 0; i_cell < n_elmts[i_part]; ++i_cell) {
      for(int i_data = 0; i_data < cell_stri[i_part][i_cell]; ++i_data ){
        _pface_ln_to_gn[idx_data++] = PDM_ABS(pcell_face_tmp[i_part][idx_data]);
      }
    }

    /*
     * Deduce ln_to_gn
     */
    // printf("Sort data between : 0 and %d \n", idx_data);
    int* unique_order = (int *) malloc( idx_data * sizeof(int));
    int n_elmt_sort = PDM_inplace_unique_long2(_pface_ln_to_gn, unique_order, 0, idx_data-1);
    _n_faces[i_part] = n_elmt_sort;

    if(0 == 1){
      printf("n_elmt_sort::%d\n", n_elmt_sort);
      printf("_pface_ln_to_gn::");
      for(int i = 0; i < n_elmt_sort; ++i){
        printf("%d ", _pface_ln_to_gn[i]);
      }
      printf("\n");
    }

    /*
     * Realloc
     */
    _face_ln_to_gn[i_part] = (PDM_g_num_t *) realloc(_face_ln_to_gn[i_part], n_elmt_sort * sizeof(PDM_g_num_t) );
    _pface_ln_to_gn = _face_ln_to_gn[i_part];

    /*
     *  We need to regenerate the connectivity and pass it in local numbering
     *  This one is vectorisable if we remove PDM_binary_seach and idx_data++
     */
    // idx_data = 0;
    // for(int i_cell = 0; i_cell < n_elmts[i_part]; ++i_cell) {
    //   for(int i_data = 0; i_data < cell_stri[i_part][i_cell]; ++i_data ){
    //     int         g_sgn  = PDM_SIGN(pcell_face_tmp[i_part][idx_data]);
    //     PDM_g_num_t g_elmt = PDM_ABS (pcell_face_tmp[i_part][idx_data]);
    //     int l_elmt         = unique_order[idx_data];
    //     int l_elmt_old     = PDM_binary_search_long(g_elmt, _pface_ln_to_gn, n_elmt_sort); /* In [0, n_elmt_sort-1] */

    //     // printf("[%d] - Search [%d] --> %d | %d \n", idx_data, g_elmt, l_elmt, l_elmt_old);

    //     /* Overwrite the pcell_face with local numbering and reput sign on it */
    //     _pcell_face[idx_data++] = (l_elmt + 1) * g_sgn ;
    //   }
    // }

    /*
     * Overpowered loop
     */
    #pragma simd
    for(int idx = 0; idx < _cell_face_idx[i_part][n_elmts[i_part]]; ++idx) {
      int g_sgn  = PDM_SIGN(pcell_face_tmp[i_part][idx]);
      int l_elmt = unique_order[idx];
      _pcell_face[idx] = (l_elmt + 1) * g_sgn ;
    }

    free(unique_order);
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
    free(pcell_face_tmp[i_part]);
  }
  free(cell_stri);
  free(pcell_face_tmp);

}

/**
 *  \brief
 */
void
PDM_generate_entity_graph_comm
(
 const PDM_MPI_Comm    comm,
 PDM_g_num_t          *part_distribution,
 PDM_g_num_t          *entity_distribution,
 int                   n_part,
 int                  *n_entities,
 PDM_g_num_t         **pentity_ln_to_gn,
 int                ***pproc_bound_idx,
 int                ***ppart_bound_idx,
 int                ***pentity_bound_idx
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_g_num_t* entity_distribution_ptb = (PDM_g_num_t * ) malloc( sizeof(PDM_g_num_t) * (n_rank+1) );
  for(int i = 0; i < n_rank+1; ++i){
    entity_distribution_ptb[i] = entity_distribution[i] - 1;
  }

  /*
   * First part : we put in data the triplet (iRank, iPart, iEntityLoc)
   */
  int** part_stri = (int ** ) malloc( n_part * sizeof(int *));
  int** part_data = (int ** ) malloc( n_part * sizeof(int *));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    part_stri[i_part] = (int *) malloc(     n_entities[i_part] * sizeof(int));
    part_data[i_part] = (int *) malloc( 3 * n_entities[i_part] * sizeof(int));
    for(int i_entity = 0; i_entity < n_entities[i_part]; ++i_entity) {
      part_data[i_part][3*i_entity  ] = i_rank;
      part_data[i_part][3*i_entity+1] = i_part;
      part_data[i_part][3*i_entity+2] = i_entity;

      part_stri[i_part][i_entity] = 3;
    }
  }

  /*
   * Setup protocol exchange
   */
  PDM_part_to_block_t *ptb =
   PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                             PDM_PART_TO_BLOCK_POST_MERGE,
                             1.,
                             pentity_ln_to_gn,
                             entity_distribution_ptb,
                             n_entities,
                             n_part,
                             comm);
   const int n_entity_block = PDM_part_to_block_n_elt_block_get(ptb);

   /*
    * Exchange
    */
  int*         blk_stri = NULL;
  PDM_g_num_t* blk_data = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(int),
                          PDM_STRIDE_VAR,
                          1,
                          part_stri,
                (void **) part_data,
                          &blk_stri,
                (void **) &blk_data);

  /*
   * Free
   */
  PDM_part_to_block_free(ptb);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(part_stri[i_part]);
    free(part_data[i_part]);
  }
  free(part_stri);
  free(part_data);

  /*
   * Panic verbose
   */
  if(0 == 1){
    int idx_debug_data = 0;
    for(int i_block = 0; i_block < n_entity_block; ++i_block){
      printf("[%d]-blk_stri[%d]:: ", i_block, blk_stri[i_block]);
      for(int i_data = 0; i_data < blk_stri[i_block]; ++i_data){
        printf("%d ", blk_data[idx_debug_data++]);
      }
      printf("\n");
    }
  }

  /*
   * Post-treatment : For all non shared data we have blk_stri == 3
   *                  And for others we have n_shared * 3 data
   *                  In blk view we can easily remove all non shared data
   *                  We reexchange with block_to_part only the shared data
   */
  int idx_comp = 0;     /* Compressed index use to fill the buffer */
  int idx_data = 0;     /* Index in the block to post-treat        */

  for(int i_block = 0; i_block < n_entity_block; ++i_block){

    /* Non shared data --> Compression */
    if(blk_stri[i_block] == 3){
      blk_stri[i_block] = 0;
      idx_data += 3;          /* Mv directly to the next */
    } else {
      for(int i_data = 0; i_data < blk_stri[i_block]; ++i_data){
        blk_data[idx_comp++] = blk_data[idx_data++];
      }
    }
  }

  /*
   * Compress data
   */
  blk_data = (int *) realloc(blk_data, idx_comp * sizeof(int));

  /*
   * Panic verbose
   */
  if(0 == 1){
    int idx_debug_data = 0;
    for(int i_block = 0; i_block < n_entity_block; ++i_block){
      printf("[%d]-blk_stri[%d]:: ", i_block, blk_stri[i_block]);
      for(int i_data = 0; i_data < blk_stri[i_block]; ++i_data){
        printf("%d ", blk_data[idx_debug_data++]);
      }
      printf("\n");
    }
  }

  /*
   * All data is now sort we cen resend to partition
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(entity_distribution_ptb,
                               (const PDM_g_num_t **) pentity_ln_to_gn,
                                                      n_entities,
                                                      n_part,
                                                      comm);

  PDM_block_to_part_exch2(btp,
                          sizeof(int),
                          PDM_STRIDE_VAR,
                          blk_stri,
             (void *  )   blk_data,
             (int  ***)  &part_stri,
             (void ***)  &part_data);

  /*
   * Free
   */
  free(blk_data);
  free(blk_stri);


  /*
   * Interface for pdm_part is  :
   *    - Face local number
   *    - Connected process
   *    - Connected partition on the connected process
   *    - Connected face local number in the connected partition
   */

  /* Allocate */
  *pproc_bound_idx   = (int **) malloc( ( n_part ) * sizeof(int * ));
  *ppart_bound_idx   = (int **) malloc( ( n_part ) * sizeof(int * ));
  *pentity_bound_idx = (int **) malloc( ( n_part ) * sizeof(int * ));

  /* Shortcut */
  int** _pproc_bound_idx   = *pproc_bound_idx;
  int** _ppart_bound_idx   = *ppart_bound_idx;
  int** _pentity_bound_idx = *pentity_bound_idx;


  /*
   * Post-treatment in partition
   */
  for(int i_part = 0; i_part < n_part; ++i_part){

    /* Shortcut */
    int* _part_stri = part_stri[i_part];
    int* _part_data = part_data[i_part];

    /* First pass to count */
    int idx_part  = 0;
    int n_connect = 0;
    for(int i_entity = 0; i_entity < n_entities[i_part]; ++i_entity) {
      // printf("[%d] part_stri::[%d] -> ", i_entity, _part_stri[i_entity]);
      for(int i_data = 0; i_data < _part_stri[i_entity]; ++i_data) {
        // printf("%d ", _part_data[idx_part++]);
        n_connect++;
      }
      // printf("\n");
    }

    /* Rebuild in a specific array all the information to sort properly */
    n_connect = n_connect/3;
    PDM_g_num_t* connect_entity = (PDM_g_num_t * ) malloc( n_connect * 3 * sizeof(PDM_g_num_t ));
    int*         connect_info   = (int         * ) malloc( n_connect * 2 * sizeof(int         ));
    // printf("n_connect::%d\n", n_connect);

    idx_part  = 0;
    n_connect = 0;
    for(int i_entity = 0; i_entity < n_entities[i_part]; ++i_entity) {
      for(int i_data = 0; i_data < _part_stri[i_entity]/3; ++i_data) {
        // printf(" idx_part = %d\n", idx_part);
        int opp_rank = _part_data[3*idx_part  ];
        int opp_part = _part_data[3*idx_part+1];

        /* On enleve de part_data les informations inutiles (car connu sur cette partition ) */
        int is_distant_bound = 1;
        if( opp_rank == i_rank ){
          if( opp_part == i_part ) {
            is_distant_bound = 0;
          }
        }

        if( is_distant_bound ) {
          connect_entity[3*n_connect  ] = opp_rank;
          connect_entity[3*n_connect+1] = opp_part;
          connect_entity[3*n_connect+2] = pentity_ln_to_gn[i_part][i_entity];

          connect_info[2*n_connect  ] = i_entity;
          connect_info[2*n_connect+1] = _part_data[3*idx_part+2];

          n_connect += 1;
        }
        idx_part += 1;
      }
    }
    printf("n_connect_sort::%d\n", n_connect);

    /*
     * Sort
     */
    int* order = (int *) malloc( n_connect * sizeof(int) );
    PDM_order_gnum_s(connect_entity, 3, order, n_connect);

    /*
     * Panic verbose
     */
    if( 0 == 1){
      printf("order::\n");
      for(int i = 0; i < n_connect; ++i){
        printf("%d ", order[i]);
      }
      printf("\n");
    }

    /* We need to recompute for each opposite part */
    _pproc_bound_idx[i_part]   = (int *) malloc( ( n_rank + 1                ) * sizeof(int));
    _ppart_bound_idx[i_part]   = (int *) malloc( ( part_distribution[n_rank] ) * sizeof(int));
    _pentity_bound_idx[i_part] = (int *) malloc( ( 4 * n_connect             ) * sizeof(int));

    for(int i = 0; i < n_rank+1; ++i){
      _pproc_bound_idx[i_part][i] = 0;
    }
    for(int i = 0; i < part_distribution[n_rank]; ++i){
      _ppart_bound_idx[i_part][i] = 0;
    }

    /* Rebuild */
    for(int i = 0; i < n_connect; ++i){
      int idx = order[i];

      int opp_proc   = connect_entity[3*idx  ];
      int opp_part   = connect_entity[3*idx+1];

      // printf(" i::%d | idx::%d \n", i, idx);
      // printf(" opp_proc::%d | opp_part::%d \n", opp_proc, opp_part);
      int opp_part_g = opp_part + part_distribution[opp_proc] - 1;

      _pproc_bound_idx[i_part][opp_proc  +1] += 1;
      _ppart_bound_idx[i_part][opp_part_g+1] += 1;

      _pentity_bound_idx[i_part][4*i  ] = connect_info[2*idx  ];
      _pentity_bound_idx[i_part][4*i+3] = connect_info[2*idx+1];

      _pentity_bound_idx[i_part][4*i+1] = opp_proc;
      _pentity_bound_idx[i_part][4*i+2] = opp_part;

    }

    /* Transform in idex */
    for(int i = 0; i < n_rank; ++i){
      _pproc_bound_idx[i_part][i+1] += _pproc_bound_idx[i_part][i];
    }
    for(int i = 0; i < part_distribution[n_rank]-1; ++i){
      _ppart_bound_idx[i_part][i+1] +=  _ppart_bound_idx[i_part][i];
    }

    free(order);
    free(connect_entity);
    free(connect_info);

  }

  /*
   * Panic Verbose
   */
  if(0 == 1){
    for(int i_part = 0; i_part < n_part; ++i_part){

      printf("[%d] _pproc_bound_idx::", i_part);
      for(int i = 0; i < n_rank+1; ++i){
        printf("%d ", _pproc_bound_idx[i_part][i]);
      }
      printf("\n");

      printf("[%d] _ppart_bound_idx::", i_part);
      for(int i = 0; i < part_distribution[n_rank]; ++i){
        printf("%d ", _ppart_bound_idx[i_part][i]);
      }
      printf("\n");
    }

  }



  /*
   * Free
   */
  for(int i_part = 0; i_part < n_part; ++i_part){
    free(part_stri[i_part]);
    free(part_data[i_part]);
  }
  free(part_stri);
  free(part_data);
  free(entity_distribution_ptb);
  PDM_block_to_part_free(btp);



}


#ifdef __cplusplus
}
#endif /* __cplusplus */
