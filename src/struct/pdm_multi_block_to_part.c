/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_multi_block_to_part.h"
#include "pdm_multi_block_to_part_priv.h"
#include "pdm_binary_search.h"
#include "pdm_priv.h"
#include "pdm_logging.h"

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

  mbtp->block_distrib_idx = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t*) * (mbtp->n_block));

  PDM_g_num_t shift = 0;
  for(int i_block = 0; i_block < mbtp->n_block; ++i_block) {
    mbtp->block_distrib_idx[i_block] = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mbtp->n_rank + 1));
    for(int i = 0; i < mbtp->n_rank + 1; ++i) {
      mbtp->block_distrib_idx[i_block][i] = shift + block_distrib_idx[i_block][i];
    }
    shift += mbtp->block_distrib_idx[i_block][mbtp->n_rank];
    PDM_log_trace_array_long(mbtp->block_distrib_idx[i_block], mbtp->n_rank + 1, "mbtp->block_distrib_idx:: ");
  }

  mbtp->n_part  = n_part;

  mbtp->n_data_block = mbtp->n_rank * mbtp->n_block; /* All ranks have the same number of blocks */

  mbtp->requested_block_idx = malloc (sizeof(int) * (mbtp->n_data_block + 1));
  mbtp->requested_block_n   = malloc (sizeof(int) * (mbtp->n_data_block    ));
  for (int i = 0; i < mbtp->n_data_block; i++) {
    mbtp->requested_block_idx[i] = 0;
    mbtp->requested_block_n  [i] = 0;
  }

  mbtp->n_elt = malloc (sizeof(int  ) * n_part);
  mbtp->ind   = malloc (sizeof(int *) * n_part);

  for (int i = 0; i < n_part; i++) {

    mbtp->n_elt[i] = n_elt[i];
    mbtp->ind[i]   = malloc (sizeof(int) * n_elt[i]);

    const PDM_g_num_t *_gnum_elt = gnum_elt[i];

    for (int j = 0; j < n_elt[i]; j++) {

      int idx_block = PDM_binary_search_gap_long(_gnum_elt[j] - 1,
                                                 multi_distrib_idx,
                                                 mbtp->n_block + 1);

      int idx_rank     = PDM_binary_search_gap_long(_gnum_elt[j] - 1,
                                                    mbtp->block_distrib_idx[idx_block],
                                                    mbtp->n_rank + 1);

      // printf("[i_part:%i | ielm : %i -> i_block : %i | ind : % i\n", i, j, i_block, ind);
      int idx_data_block = idx_block + idx_rank*mbtp->n_block;
      mbtp->requested_block_n[idx_data_block]++;

    }
  }

  for (int i = 0; i < mbtp->n_data_block; i++) {
    mbtp->requested_block_idx[i+1] = mbtp->requested_block_idx[i] +
                                    mbtp->requested_block_n[i];
  }

  int s_requested_data = mbtp->requested_block_idx[mbtp->n_data_block - 1]
                       + mbtp->requested_block_n[mbtp->n_data_block - 1];

  printf("s_requested_data = %i \n", s_requested_data);

  int *requested_data = malloc( sizeof(int) * s_requested_data);

  for (int i = 0; i < mbtp->n_data_block; i++) {
    mbtp->requested_block_n[i] = 0;
  }

  for (int i = 0; i < n_part; i++) {

    const PDM_g_num_t *_gnum_elt = gnum_elt[i];

    for (int j = 0; j < n_elt[i]; j++) {

      int idx_block = PDM_binary_search_gap_long(_gnum_elt[j] - 1,
                                                 multi_distrib_idx,
                                                 mbtp->n_block + 1);

      int idx_rank     = PDM_binary_search_gap_long(_gnum_elt[j] - 1,
                                                    mbtp->block_distrib_idx[idx_block],
                                                    mbtp->n_rank + 1);

      int idx_data_block = idx_block + idx_rank*mbtp->n_block;

      int idx = mbtp->requested_block_idx[idx_data_block] + mbtp->requested_block_n[idx_data_block]++;

      // On peut ici faire une double indirection donc l'idx depend du numero de block et du i,j
      mbtp->ind[i][j] = idx;

      // Il faut se faire un double index
      // En fait pour la recuperation on doit pouvoir utiliser
      // les requested_block + idx de deuxieme niveau !!! Car ce sera l'ordre d'arrvié !!!

      // Very important here to shift twice !!!
      PDM_g_num_t _mshift         = mbtp->block_distrib_idx[idx_block][idx_rank] - mbtp->multi_distrib_idx[idx_block];
      PDM_g_num_t _requested_data = _gnum_elt[j] - 1 - _mshift;

      // Dans le requested data il faut soit rajouter les indices de demandes de chaque bloc
      // Soit pour chaque demande le numero de bloc
      requested_data[idx] = (int) _requested_data;
      // requested_block[i_block][idx] = (int) _requested_data;

    }
  }

  mbtp->distributed_block_n = malloc (sizeof(int) * mbtp->n_data_block);

  PDM_MPI_Alltoall (mbtp->requested_block_n,   mbtp->n_block, PDM_MPI_INT,
                    mbtp->distributed_block_n, mbtp->n_block, PDM_MPI_INT,
                    comm);


  mbtp->distributed_block_idx = malloc (sizeof(int) * (mbtp->n_data_block + 1));
  mbtp->distributed_block_idx[0] = 0;

  for (int i = 0; i < mbtp->n_data_block; i++) {
    mbtp->distributed_block_idx[i+1] = mbtp->distributed_block_n[i] +
                                       mbtp->distributed_block_idx[i];
  }

  PDM_log_trace_array_int(mbtp->requested_block_n    , mbtp->n_data_block  , "mbtp->requested_block_n    :: ");
  PDM_log_trace_array_int(mbtp->requested_block_idx  , mbtp->n_data_block+1, "mbtp->requested_block_idx  :: ");
  PDM_log_trace_array_int(mbtp->distributed_block_n  , mbtp->n_data_block  , "mbtp->distributed_block_n  :: ");
  PDM_log_trace_array_int(mbtp->distributed_block_idx, mbtp->n_data_block+1, "mbtp->distributed_block_idx:: ");


  mbtp->distributed_data = malloc(sizeof(int) * mbtp->distributed_data_idx[mbtp->n_data_block]);

  // Les data se deduisent des blocks
  mbtp->requested_data_idx   = malloc (sizeof(int) * (mbtp->n_rank + 1));
  mbtp->requested_data_n     = malloc (sizeof(int) * (mbtp->n_rank    ));
  mbtp->distributed_data_idx = malloc (sizeof(int) * (mbtp->n_rank + 1));
  mbtp->distributed_data_n   = malloc (sizeof(int) * (mbtp->n_rank    ));

  for(int i = 0; i < mbtp->n_rank+1; ++i) {
    int ind = i*mbtp->n_block;
    mbtp->requested_data_idx  [i] = mbtp->requested_block_idx  [ind];
    mbtp->distributed_data_idx[i] = mbtp->distributed_block_idx[ind];
    mbtp->requested_data_n    [i] = 0;
    mbtp->distributed_data_n  [i] = 0;
    for(int i_block = 0; i_block < mbtp->n_block; ++i_block) {
      mbtp->requested_data_n    [i] += mbtp->requested_block_n    [i_block+ind];
      mbtp->distributed_data_n  [i] += mbtp->distributed_block_n  [i_block+ind];
    }
  }
  PDM_log_trace_array_int(mbtp->requested_data_n    , mbtp->n_rank  , "mbtp->requested_data_n    :: ");
  PDM_log_trace_array_int(mbtp->requested_data_idx  , mbtp->n_rank+1, "mbtp->requested_data_idx  :: ");
  PDM_log_trace_array_int(mbtp->distributed_data_n  , mbtp->n_rank  , "mbtp->distributed_data_n  :: ");
  PDM_log_trace_array_int(mbtp->distributed_data_idx, mbtp->n_rank+1, "mbtp->distributed_data_idx:: ");

  PDM_MPI_Alltoallv (requested_data,
                     mbtp->requested_data_n,
                     mbtp->requested_data_idx,
                     PDM_MPI_INT,
                     mbtp->distributed_data,
                     mbtp->distributed_data_n,
                     mbtp->distributed_data_idx,
                     PDM_MPI_INT,
                     comm);

  PDM_log_trace_array_int(mbtp->distributed_data    , mbtp->distributed_data_idx[mbtp->n_rank], "mbtp->distributed_data:: ");

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

  unsigned char **_part_data;

  printf("PDM_multi_block_to_part_exch2\n");

  size_t *i_send_buffer = (size_t *) malloc (sizeof(size_t) * _mbtp->n_rank);
  size_t *i_recv_buffer = (size_t *) malloc (sizeof(size_t) * _mbtp->n_rank);
  int *n_send_buffer = (int *) malloc (sizeof(int) * _mbtp->n_rank);
  int *n_recv_buffer = (int *) malloc (sizeof(int) * _mbtp->n_rank);

  for (int i = 0; i < _mbtp->n_rank; i++) {
    n_send_buffer[i] = 0;
    n_recv_buffer[i] = 0;
    i_send_buffer[i] = 0;
    i_recv_buffer[i] = 0;
  }

  unsigned char *send_buffer = NULL;
  unsigned char *recv_buffer = NULL;

  size_t s_send_buffer = 0;
  size_t s_recv_buffer = 0;

  int n_rank1 = _mbtp->n_rank - 1;

  // We have to treat all block
  // for(int i_block = 0; i_block < _mbtp->n_block; ++i_block) {
  // }


  if (t_stride == PDM_STRIDE_VAR) {
    abort();
  } else if (t_stride == PDM_STRIDE_CST) {
    int cst_stride = -1;
    for(int i_block = 0; i_block < _mbtp->n_block; ++i_block) {
      if(cst_stride == -1){
        cst_stride = block_stride[i_block][0];
        printf("cst_stride   = %i \n", cst_stride  );
      } else {
        printf("block_stride[%i] = %i \n", i_block,  block_stride[i_block][0] );

        assert( cst_stride == block_stride[i_block][0]);
      }
    }

    int s_block_unit = cst_stride * (int) s_data;
    printf("cst_stride   = %i \n", cst_stride  );
    printf("s_block_unit = %i \n", s_block_unit);

    for (int i = 0; i < _mbtp->n_rank; i++) {

      i_send_buffer[i] = _mbtp->distributed_data_idx[i] * s_block_unit;
      i_recv_buffer[i] = _mbtp->requested_data_idx  [i] * s_block_unit;

      n_send_buffer[i] = _mbtp->distributed_data_n[i] * s_block_unit;
      n_recv_buffer[i] = _mbtp->requested_data_n  [i] * s_block_unit;
    }

    s_send_buffer = i_send_buffer[n_rank1] + n_send_buffer[n_rank1];
    s_recv_buffer = i_recv_buffer[n_rank1] + n_recv_buffer[n_rank1];

    send_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_send_buffer);
    recv_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_buffer);

    // int idx1 = 0;
    // for (int i = 0; i < s_distributed_data; i++) {
    //   int ind = _btp->distributed_data[i]; /* List des indices que le rang opposé souhaité recupérer */
    //   unsigned char *_block_data_deb = _block_data + ind * cst_stride * (int) s_data;
    //   for (int k = 0; k < s_block_unit; k++) {
    //     send_buffer[idx1++] = _block_data_deb[k];
    //   }
    // }


    // for(int i_block = 0; i_block < _mbtp->n_block; ++i_block) {
    //   unsigned char* _block_data = (unsigned char *) block_data[i_block];
    //   int* ind_block = &distribute_data[distributed_data_idx[i_block]];
    //   for(int i = 0; i < s_distributed_data[i_block]; ++i) {
    //      int ind = ind_block[i];
    //
    //   }
    // }


  }

  /*
   * Data exchange
   */

  PDM_MPI_Alltoallv_l(send_buffer,
                      n_send_buffer,
                      i_send_buffer,
                      PDM_MPI_BYTE,
                      recv_buffer,
                      n_recv_buffer,
                      i_recv_buffer,
                      PDM_MPI_BYTE,
                      _mbtp->comm);
  free(send_buffer);
  free(n_send_buffer);
  free(i_send_buffer);
  free(n_recv_buffer);
  free(i_recv_buffer);

  if (t_stride == PDM_STRIDE_VAR) {
    abort();
  } else if (t_stride == PDM_STRIDE_CST) {

    // const int cst_stride = *block_stride;
    // const int s_block_unit = cst_stride * (int) s_data;

    // for (int i = 0; i < _btp->n_part; i++) {

    //   _part_data[i] = malloc(sizeof(unsigned char) * s_block_unit * _btp->n_elt[i]);

    //   for (int j = 0; j < _mbtp->n_elt[i]; j++) {

    //     int idx1  = j * s_block_unit;
    //     int idx2 = _btp->ind[i][j] * s_block_unit;

    //     for (int k = 0; k < s_block_unit; k++) {
    //        _part_data[i][idx1+k] = recv_buffer[idx2+k];
    //     }
    //   }
    // }
  }

  free(recv_buffer);

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
  free (_mbtp->distributed_block_idx);
  free (_mbtp->distributed_block_n);
  free (_mbtp->requested_data_idx);
  free (_mbtp->requested_data_n);
  free (_mbtp->requested_block_idx);
  free (_mbtp->requested_block_n);


  free (_mbtp);

  return NULL;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
