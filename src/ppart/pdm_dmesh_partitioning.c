
/*============================================================================
 * Parallel partitioning
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

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

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"
#include "pdm_handles.h"

#include "pdm_dmesh_partitioning.h"

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
void
PDM_generate_part_cell_ln_to_gn
(
 const PDM_MPI_Comm    comm,
 PDM_g_num_t          *part_distribution,
 PDM_g_num_t          *cell_distribution,
 int                  *dcell_face_idx,
 PDM_g_num_t          *dcell_face,
 int                  *cell_part,
 int                ***pcell_face_idx,
 int                ***pcell_face,
 int                ***pcell_ln_to_gn
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

  // printf("part_distribution::");
  // for(int i = 0; i < n_rank+1; ++i){
  //   printf("%d ", part_distribution[i]);
  // }
  // printf("\n");

  // printf("cell_distribution::");
  // for(int i = 0; i < n_rank+1; ++i){
  //   printf("%d ", cell_distribution[i]);
  // }
  // printf("\n");

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
   * cell_ln_to_gn
   */
  int* dcell_stri = (int * ) malloc( sizeof(int) * dn_cell );
  PDM_g_num_t* dcell_ln_to_gn = (PDM_g_num_t * ) malloc( sizeof(PDM_g_num_t) * dn_cell );

  PDM_g_num_t shift_g = cell_distribution[i_rank];
  for(int i = 0; i < dn_cell; ++i){
    dcell_stri[i] = 1;
    dcell_ln_to_gn[i] = (PDM_g_num_t) shift_g + i;
  }
  int* pcell_stri = NULL;
  PDM_g_num_t* pcell_ln_to_gn_tmp = NULL;
  PDM_part_to_block_exch (ptb_partition,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          1,
                          &dcell_stri,
                (void **) &dcell_ln_to_gn,
                          &pcell_stri,
                (void **) &pcell_ln_to_gn_tmp);


  free(dcell_stri);

  printf("n_part_block::%d\n",n_part_block );
  if(0 == 1){
    int idx_block = 0;
    for(int i = 0; i < n_part_block; ++i){
      printf(" pcell_stri = %d ---> ", pcell_stri[i]);
      for(int i_data = 0; i_data < pcell_stri[i]; ++i_data){
        printf("%d ", pcell_ln_to_gn_tmp[idx_block]);
        idx_block++;
      }
      printf("\n");
    }
  }


   /*
    * En changeant un poil le part_to_block on eviterai un tableau temporaire
    */
  int* dcell_face_n = (int * ) malloc( dn_cell * sizeof(int));
  for(int i = 0; i < dn_cell; ++i){
    dcell_face_n[i] = dcell_face_idx[i+1] - dcell_face_idx[i];
  }

  int* cell_face_n = NULL;
  PDM_g_num_t* cell_face   = NULL;
  printf("PDM_part_to_block_exch \n");
  PDM_part_to_block_exch (ptb_partition,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          1,
                          &dcell_face_n,
                (void **) &dcell_face,
                          &cell_face_n,
                (void **) &cell_face);
  printf("PDM_part_to_block_exch end \n");


  int idx_block = 0;
  for(int i = 0; i < n_part_block; ++i){
    printf(" pcell_stri = %d ---> ", cell_face_n[i]);
    for(int i_data = 0; i_data < cell_face_n[i]; ++i_data){
      printf("%d ", cell_face[idx_block]);
      idx_block++;
    }
    printf("\n");
  }

  /*
   * Si besoin du cell_face on doit faire un echange pour connaitre le nombre d'element par faces
   *    --> Pas besoin du cell_face tout le temps non ?
   */





  PDM_part_to_block_free (ptb_partition);

  free(dpart_ln_to_gn);
  free(dcell_face_n);
  free(part_distribution_ptb);
}


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
  _dmesh_partitioning_t* _dmp = _get_from_id(dmesh_partitioning_id);
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
