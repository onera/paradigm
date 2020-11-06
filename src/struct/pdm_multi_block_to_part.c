/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_multi_block_to_part.h"
#include "pdm_multi_block_to_part_priv.h"
#include "pdm_binary_search.h"
#include "pdm_priv.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a block to partitions redistribution
 *
 * \param [in]   multi_distrib_idx Multiple block distribution (size : \ref size of \ref nblock + 1)
 * \param [in]   block_distrib_idx Block distribution (size : \ref size of \ref comm + 1)
 * \param [in]   gnum_elt          Element global number (size : \ref n_part)
 * \param [in]   n_elt             Local number of elements (size : \ref n_part)
 * \param [in]   n_part            Number of partition
 * \param [in]   comm              MPI communicator
 *
 * \return   Initialized \ref PDM_block_to_part instance
 *
 */

PDM_multi_block_to_part_t *
PDM_multi_block_to_part_create
(
 const PDM_g_num_t   *multi_distrib_idx,
 const int            n_block,
 const PDM_g_num_t  **block_distrib_idx,
 const PDM_g_num_t  **gnum_elt,
 const int           *n_elt,
 const int            n_part,
 const PDM_MPI_Comm   comm
)
{
  printf(" ola que tal PDM_multi_block_to_part_t \n");

  _pdm_multi_block_to_part_t *mbtp =
    (_pdm_multi_block_to_part_t *) malloc (sizeof(_pdm_multi_block_to_part_t));

  mbtp->comm       = comm;
  mbtp->pttopt_comm = 0;

  PDM_MPI_Comm_size (comm, &mbtp->n_rank);
  PDM_MPI_Comm_rank (comm, &mbtp->i_rank);

  mbtp->n_block = n_block;

  mbtp->multi_distrib_idx = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mbtp->n_block + 1));
  for(int i_block = 0; i_block < mbtp->n_block+1; ++i_block){
    mbtp->multi_distrib_idx[i_block] = multi_distrib_idx[i_block];
  }

  mbtp->block_distrib_idx = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t*) * (mbtp->n_block + 1));

  for(int i_block = 0; i_block < mbtp->n_block; ++i_block) {
    mbtp->block_distrib_idx[i_block] = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mbtp->n_rank + 1));
    for(int i = 0; i < mbtp->n_rank + 1; ++i) {
      mbtp->block_distrib_idx[i_block][i] = block_distrib_idx[i_block][i];
    }
  }

  mbtp->n_part  = n_part;

  mbtp->requested_data_idx = malloc (sizeof(int) * (mbtp->n_rank + 1));
  mbtp->requested_data_n   = malloc (sizeof(int) * mbtp->n_rank);
  for (int i = 0; i < mbtp->n_rank; i++) {
    mbtp->requested_data_idx[i] = 0;
    mbtp->requested_data_n  [i] = 0;
  }

  mbtp->n_elt = malloc (sizeof(int) * n_part);
  mbtp->ind   = malloc (sizeof(int *) * n_part);

  for (int i = 0; i < n_part; i++) {

    mbtp->n_elt[i] = n_elt[i];
    mbtp->ind[i]   = malloc (sizeof(int) * n_elt[i]);

    const PDM_g_num_t *_gnum_elt = gnum_elt[i];

    for (int j = 0; j < n_elt[i]; j++) {

      int i_block = PDM_binary_search_gap_long (_gnum_elt[j] - 1,
                                                multi_distrib_idx,
                                                mbtp->n_block + 1);
      int ind     = PDM_binary_search_gap_long (_gnum_elt[j] - 1,
                                                mbtp->block_distrib_idx[i_block],
                                                mbtp->n_rank + 1);

      mbtp->requested_data_n[ind]++;

    }
  }

  for (int i = 0; i < mbtp->n_rank; i++) {
    mbtp->requested_data_idx[i+1] = mbtp->requested_data_idx[i] +
                                    mbtp->requested_data_n[i];
  }

  int s_requested_data = mbtp->requested_data_idx[mbtp->n_rank - 1]
                       + mbtp->requested_data_n[mbtp->n_rank - 1];

  printf("s_requested_data = %i \n", s_requested_data);

  int *requested_data = malloc (sizeof(int) *  s_requested_data);

  for (int i = 0; i < mbtp->n_rank; i++) {
    mbtp->requested_data_n[i] = 0;
  }

  for (int i = 0; i < n_part; i++) {

    const PDM_g_num_t *_gnum_elt = gnum_elt[i];

    for (int j = 0; j < n_elt[i]; j++) {

      int i_block = PDM_binary_search_gap_long (_gnum_elt[j] - 1,
                                                multi_distrib_idx,
                                                mbtp->n_block + 1);
      int ind     = PDM_binary_search_gap_long (_gnum_elt[j] - 1,
                                                mbtp->block_distrib_idx[i_block],
                                                mbtp->n_rank + 1);
      int idx = mbtp->requested_data_idx[ind] + mbtp->requested_data_n[ind]++;

      mbtp->ind[i][j] = idx;

      // Very important here to shift twice !!!
      PDM_g_num_t _mshift         = mbtp->block_distrib_idx[i_block][ind] - mbtp->multi_distrib_idx[i_block];
      PDM_g_num_t _requested_data = _gnum_elt[j] - 1 - _mshift;

      requested_data[idx] = (int) _requested_data;
      // requested_data[i_block][idx] = (int) _requested_data;

    }
  }


  mbtp->distributed_data_n = malloc (sizeof(int) * mbtp->n_rank);

  PDM_MPI_Alltoall (mbtp->requested_data_n,   1, PDM_MPI_INT,
                    mbtp->distributed_data_n, 1, PDM_MPI_INT,
                    comm);

  mbtp->distributed_data_idx = malloc (sizeof(int) * (mbtp->n_rank + 1));
  mbtp->distributed_data_idx[0] = 0;

  for (int i = 0; i < mbtp->n_rank; i++) {
    mbtp->distributed_data_idx[i+1] = mbtp->distributed_data_n[i] +
                                     mbtp->distributed_data_idx[i];
  }

  mbtp->distributed_data = malloc (sizeof(int) *
                                   mbtp->distributed_data_idx[mbtp->n_rank]);


  PDM_MPI_Alltoallv (requested_data,
                     mbtp->requested_data_n,
                     mbtp->requested_data_idx,
                     PDM_MPI_INT,
                     mbtp->distributed_data,
                     mbtp->distributed_data_n,
                     mbtp->distributed_data_idx,
                     PDM_MPI_INT,
                     comm);

  free (requested_data);

  return (PDM_multi_block_to_part_t *) mbtp;
}

/**
 *
 * \brief Initialize an exchange
 * (part_stride and part_data are allocated in function)
 *
 * \param [in]   btp          Block to part structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   block_stride Stride for each block element for \ref PDM_STRIDE_VAR
 *                            Constant stride for \ref PDM_STRIDE_VAR
 * \param [in]   block_data   Block data
 * \param [out]  part_stride  Partition stride or NULL
 * \param [out]  part_data    Partition data
 *
 */

void
PDM_multi_block_to_part_exch2
(
 PDM_multi_block_to_part_t   *mbtp,
 size_t                       s_data,
 PDM_stride_t                 t_stride,
 int                        **block_stride,
 void                       **block_data,
 int                       ***part_stride,
 void                      ***part_data
)
{
  _pdm_multi_block_to_part_t *_mbtp = (_pdm_multi_block_to_part_t *) mbtp;




}


/**
 *
 * \brief Free a block to part structure
 *
 * \param [inout] btp  Block to part structure
 *
 * \return       NULL
 */

PDM_multi_block_to_part_t *
PDM_multi_block_to_part_free
(
 PDM_multi_block_to_part_t *mbtp
)
{
  _pdm_multi_block_to_part_t *_mbtp = (_pdm_multi_block_to_part_t *) mbtp;

  for(int i_block = 0; i_block < _mbtp->n_block; ++i_block) {
    free (_mbtp->block_distrib_idx[i_block]);
  }
  free(_mbtp->block_distrib_idx);

  for (int i = 0; i < _mbtp->n_part; i++) {
    free (_mbtp->ind[i]);
  }
  free (_mbtp->ind);

  free (_mbtp->n_elt);
  free (_mbtp->distributed_data);
  free (_mbtp->distributed_data_idx);
  free (_mbtp->distributed_data_n);
  free (_mbtp->requested_data_idx);
  free (_mbtp->requested_data_n);


  free (_mbtp);

  return NULL;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
