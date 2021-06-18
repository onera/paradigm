/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_partgnum1_to_partgnum2.h"
#include "pdm_partgnum1_to_partgnum2_priv.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_gnum_location.h"
#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_timer.h"
#include "pdm_priv.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_printf.h"
#include "pdm_error.h"

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
 * \brief Create a partitions to partitions redistribution
 *
 * \param [in]   gnum_elt1          Element global number (size : \ref n_part1)
 * \param [in]   n_elt1             Local number of elements (size : \ref n_part1)
 * \param [in]   n_part1            Number of partition
 * \param [in]   gnum_elt2          Element global number (size : \ref n_part2)
 * \param [in]   n_elt2             Local number of elements (size : \ref n_part2)
 * \param [in]   n_part2            Number of partition
 * \param [in]   send_to_gnum2_idx  Index of data to send to gnum2 from gnum1 
 *                                  (for each part size : \ref n_elt1+1) 
 * \param [in]   send_to_gnum2      Data to send to gnum2 from gnum1 
 * \param [in]   comm               MPI communicator
 *
 * \return   Initialized \ref PDM_partgnum1_to_partgnum2 instance
 *
 */

PDM_partgnum1_to_partgnum2_t *
PDM_partgnum1_to_partgnum2_create
(
 const PDM_g_num_t   **gnum_elt1,
 const int            *n_elt1,
 const int             n_part1,
 const PDM_g_num_t   **gnum_elt2,
 const int            *n_elt2,
 const int             n_part2,
 const int           **gnum1_to_gnum2_idx,
 const PDM_g_num_t   **gnum1_to_gnum2,
 const PDM_MPI_Comm    comm
)
{
  PDM_partgnum1_to_partgnum2_t *ptp =
    (PDM_partgnum1_to_partgnum2_t *) malloc (sizeof(PDM_partgnum1_to_partgnum2_t));

  /* Init */

  ptp->n_part1                  = n_part1;    
  ptp->gnum_elt1                = gnum_elt1;
  ptp->n_elt1                   = n_elt1;

  ptp->n_part2                  = n_part2;                
  ptp->gnum_elt2                = gnum_elt2;              
  ptp->n_elt2                   = n_elt2;                

  ptp->gnum1_to_gnum2_idx       = gnum1_to_gnum2_idx;    
  ptp->gnum1_to_gnum2           = gnum1_to_gnum2;        
  ptp->comm                     = comm;

  PDM_MPI_Comm_size (comm, &(ptp->n_rank));
  PDM_MPI_Comm_rank (comm, &(ptp->i_rank));
  int i_rank = ptp->i_rank;
  int n_rank = ptp->n_rank;


  ptp->n_ref_gnum2              = NULL;         
  ptp->ref_gnum2                = NULL;          
  ptp->n_unref_gnum2            = NULL;      
  ptp->unref_gnum2              = NULL;        
  ptp->gnum1_come_from_idx      = NULL;
  ptp->gnum1_come_from          = NULL;    

  ptp->gnum1_to_send_buffer     = NULL;  
  ptp->recv_buffer_to_ref_gnum2 = NULL;  

  ptp->default_n_send_buffer    = malloc (sizeof(int) * n_rank);
  ptp->default_i_send_buffer    = malloc (sizeof(int) * (n_rank + 1));
  ptp->default_i_send_buffer[0] = 0;
  ptp->default_n_recv_buffer    = malloc (sizeof(int) * n_rank);
  ptp->default_i_recv_buffer    = malloc (sizeof(int) * (n_rank + 1));
  ptp->default_i_recv_buffer[0] = 0;

  for (int i = 0; i < n_rank; i++) {
    ptp->default_n_send_buffer[i] = 0;
    ptp->default_n_recv_buffer[i] = 0;
  }

  ptp->async_n_exch             = 0;           
  ptp->async_l_array            = 0;          
  ptp->async_s_data             = NULL;          
  ptp->async_cst_stride         = NULL;      
  ptp->async_send_request       = NULL;    
  ptp->async_recv_request       = NULL;    
  ptp->async_send_buffer        = NULL;     
  ptp->async_recv_buffer        = NULL;     
  ptp->async_n_send_buffer      = NULL;   
  ptp->async_i_send_buffer      = NULL;   
  ptp->async_n_recv_buffer      = NULL;   
  ptp->async_i_recv_buffer      = NULL; 

  /* 1 - gnum_location in 2 1D array gnum1_to_gnum2_rank gnum1_to_gnum2_part   gnum1_to_gnum2_part elt*/

  PDM_gnum_location_t *gl = PDM_gnum_location_create (n_part2, n_part1, comm);

  for (int i = 0; i < n_part2; i++) {
    PDM_gnum_location_elements_set (gl, i, n_elt2[i], gnum_elt2[i]);
  }

  for (int i = 0; i < n_part1; i++) {
    PDM_gnum_location_requested_elements_set (gl, i, gnum1_to_gnum2_idx[i][n_elt1[i]], gnum1_to_gnum2[i]);
  }

  PDM_gnum_location_compute(gl);

  int n_total_elt = 0;
  for (int i = 0; i < n_part1; i++) {
    int *location_gnum1_to_gnum2_idx;
    int *location_gnum1_to_gnum2;
    PDM_gnum_location_get (gl,
                           i,
                           &location_gnum1_to_gnum2_idx,
                           &location_gnum1_to_gnum2);
    n_total_elt += location_gnum1_to_gnum2_idx[gnum1_to_gnum2_idx[i][n_elt1[i]]];
  }

  PDM_gnum_location_free (gl, 0);

  int *merge_gnum1_to_gnum2_rank2 = (int *) malloc (sizeof(int) * n_total_elt);
  int *merge_gnum1_to_gnum2_part2 = (int *) malloc (sizeof(int) * n_total_elt);
  int *merge_gnum1_to_gnum2_lnum2 = (int *) malloc (sizeof(int) * n_total_elt);
  int *merge_gnum1_to_gnum2_rank1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *merge_gnum1_to_gnum2_part1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *merge_gnum1_to_gnum2_lnum1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *order                      = (int *) malloc (sizeof(int) * n_total_elt);
  int *i_send_buffer              = (int *) malloc (sizeof(int) * n_total_elt);
  int *n_gnum1_to_gnum2_rank      = (int *) malloc (sizeof(int) * n_rank);
  int *idx_gnum1_to_gnum2_rank    = (int *) malloc (sizeof(int) * (n_rank + 1));

  for (int i = 0; i < n_rank; i++) {
    n_gnum1_to_gnum2_rank[i] = 0;
  }
  idx_gnum1_to_gnum2_rank[0] = 0;

  n_total_elt = 0;
  for (int i = 0; i < n_part1; i++) {
    int *location_gnum1_to_gnum2_idx;
    int *location_gnum1_to_gnum2;
    PDM_gnum_location_get (gl,
                           i,
                           &location_gnum1_to_gnum2_idx,
                           &location_gnum1_to_gnum2);
  
    for (int j = 0; j < n_elt1[i]; j++) {

      for (int k = location_gnum1_to_gnum2_idx[j]; k < location_gnum1_to_gnum2_idx[j+1]; k++) {
        int i_rank2 = location_gnum1_to_gnum2[3*k];
        n_gnum1_to_gnum2_rank[i_rank2]++;
        merge_gnum1_to_gnum2_rank2[n_total_elt] = i_rank2;
        merge_gnum1_to_gnum2_part2[n_total_elt] = location_gnum1_to_gnum2[3*k+1];
        merge_gnum1_to_gnum2_lnum2[n_total_elt] = location_gnum1_to_gnum2[3*k+2];
        merge_gnum1_to_gnum2_rank1[n_total_elt] = i_rank;
        merge_gnum1_to_gnum2_part1[n_total_elt] = i;
        merge_gnum1_to_gnum2_lnum1[n_total_elt] = j; 
        order[n_total_elt]                      = n_total_elt;
        n_total_elt++;
      }
    }
  }

  for (int i = 0; i < n_rank; i++) {
    idx_gnum1_to_gnum2_rank[i+1] = n_gnum1_to_gnum2_rank[i] + 
                                   idx_gnum1_to_gnum2_rank[i];
  }

  free (merge_gnum1_to_gnum2_rank2);
  free (merge_gnum1_to_gnum2_part2);
  free (merge_gnum1_to_gnum2_lnum2);
  free (merge_gnum1_to_gnum2_rank1);
  free (merge_gnum1_to_gnum2_part1);
  free (merge_gnum1_to_gnum2_lnum1);
  free (order                     );

  /* 2 - Sort gnum1_to_gnum2_rank gnum1_to_gnum2_part according successively rank, part and elemt  */

  PDM_sort_int (merge_gnum1_to_gnum2_rank2, order, n_total_elt);

  int *_merge_gnum1_to_gnum2_rank2 = (int *) malloc (sizeof(int) * n_total_elt);
  int *_merge_gnum1_to_gnum2_part2 = (int *) malloc (sizeof(int) * n_total_elt);
  int *_merge_gnum1_to_gnum2_lnum2 = (int *) malloc (sizeof(int) * n_total_elt);
  int *_merge_gnum1_to_gnum2_rank1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *_merge_gnum1_to_gnum2_part1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *_merge_gnum1_to_gnum2_lnum1 = (int *) malloc (sizeof(int) * n_total_elt);

  for (int i = 0; i < n_total_elt; i++) {
    i_send_buffer[i] = -2;
    _merge_gnum1_to_gnum2_part2[i] = merge_gnum1_to_gnum2_part2[order[i]];
    _merge_gnum1_to_gnum2_lnum2[i] = merge_gnum1_to_gnum2_lnum2[order[i]];
    _merge_gnum1_to_gnum2_rank1[i] = merge_gnum1_to_gnum2_rank1[order[i]];
    _merge_gnum1_to_gnum2_part1[i] = merge_gnum1_to_gnum2_part1[order[i]];
    _merge_gnum1_to_gnum2_lnum1[i] = merge_gnum1_to_gnum2_lnum1[order[i]]; 
  }

  int *_tmp_merge_gnum1_to_gnum2_rank2 = merge_gnum1_to_gnum2_rank2;
  int *_tmp_merge_gnum1_to_gnum2_part2 = merge_gnum1_to_gnum2_part2;
  int *_tmp_merge_gnum1_to_gnum2_lnum2 = merge_gnum1_to_gnum2_lnum2;
  int *_tmp_merge_gnum1_to_gnum2_rank1 = merge_gnum1_to_gnum2_rank1;
  int *_tmp_merge_gnum1_to_gnum2_part1 = merge_gnum1_to_gnum2_part1;
  int *_tmp_merge_gnum1_to_gnum2_lnum1 = merge_gnum1_to_gnum2_lnum1;

  merge_gnum1_to_gnum2_rank2 = _merge_gnum1_to_gnum2_rank2;
  merge_gnum1_to_gnum2_part2 = _merge_gnum1_to_gnum2_part2;
  merge_gnum1_to_gnum2_lnum2 = _merge_gnum1_to_gnum2_lnum2;
  merge_gnum1_to_gnum2_rank1 = _merge_gnum1_to_gnum2_rank1;
  merge_gnum1_to_gnum2_part1 = _merge_gnum1_to_gnum2_part1;
  merge_gnum1_to_gnum2_lnum1 = _merge_gnum1_to_gnum2_lnum1;

  _merge_gnum1_to_gnum2_rank2 = _tmp_merge_gnum1_to_gnum2_rank2;
  _merge_gnum1_to_gnum2_part2 = _tmp_merge_gnum1_to_gnum2_part2;
  _merge_gnum1_to_gnum2_lnum2 = _tmp_merge_gnum1_to_gnum2_lnum2;
  _merge_gnum1_to_gnum2_rank1 = _tmp_merge_gnum1_to_gnum2_rank1;
  _merge_gnum1_to_gnum2_part1 = _tmp_merge_gnum1_to_gnum2_part1;
  _merge_gnum1_to_gnum2_lnum1 = _tmp_merge_gnum1_to_gnum2_lnum1;

  int n_part2_max = 0;
  PDM_MPI_Allreduce ((int* )&n_part2, &n_part2_max, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);

  int *n_elt_part = malloc (sizeof(int) * n_part2_max);
  int *idx_elt_part = malloc (sizeof(int) * (n_part2_max + 1));
  idx_elt_part[0] = 0;

  int cpt_buff = -1;
  for (int i = 0; i < n_rank; i++) {
    int n_elt_rank = n_gnum1_to_gnum2_rank[i];
    for (int j = 0; j < n_elt_rank; j++) {
      order[j] = j;
    }

    for (int j = 0; j < n_part2_max; j++) {
      n_elt_part[j] = 0;
    }

    int *rank_i_send_buffer              = i_send_buffer              + idx_gnum1_to_gnum2_rank[i];  
    int *rank_merge_gnum1_to_gnum2_part2 = merge_gnum1_to_gnum2_part2 + idx_gnum1_to_gnum2_rank[i];
    int *rank_merge_gnum1_to_gnum2_lnum2 = merge_gnum1_to_gnum2_lnum2 + idx_gnum1_to_gnum2_rank[i];
    int *rank_merge_gnum1_to_gnum2_rank1 = merge_gnum1_to_gnum2_rank1 + idx_gnum1_to_gnum2_rank[i];
    int *rank_merge_gnum1_to_gnum2_part1 = merge_gnum1_to_gnum2_part1 + idx_gnum1_to_gnum2_rank[i];
    int *rank_merge_gnum1_to_gnum2_lnum1 = merge_gnum1_to_gnum2_lnum1 + idx_gnum1_to_gnum2_rank[i];

    int *_rank_merge_gnum1_to_gnum2_part2 = _merge_gnum1_to_gnum2_part2 + idx_gnum1_to_gnum2_rank[i];
    int *_rank_merge_gnum1_to_gnum2_lnum2 = _merge_gnum1_to_gnum2_lnum2 + idx_gnum1_to_gnum2_rank[i];
    int *_rank_merge_gnum1_to_gnum2_rank1 = _merge_gnum1_to_gnum2_rank1 + idx_gnum1_to_gnum2_rank[i];
    int *_rank_merge_gnum1_to_gnum2_part1 = _merge_gnum1_to_gnum2_part1 + idx_gnum1_to_gnum2_rank[i];
    int *_rank_merge_gnum1_to_gnum2_lnum1 = _merge_gnum1_to_gnum2_lnum1 + idx_gnum1_to_gnum2_rank[i];

    PDM_sort_int (_rank_merge_gnum1_to_gnum2_part2, order, n_elt_rank);

    int _max_part = 0;
    for (int k = 0; k < n_elt_rank; k++) {
      int i_part = rank_merge_gnum1_to_gnum2_part2[order[k]];
      n_elt_part[i_part]++;
      _max_part = PDM_MAX (_max_part, i_part);
      _rank_merge_gnum1_to_gnum2_part2[k] = i_part;
      _rank_merge_gnum1_to_gnum2_lnum2[k] = rank_merge_gnum1_to_gnum2_lnum2[order[k]];
      _rank_merge_gnum1_to_gnum2_rank1[k] = rank_merge_gnum1_to_gnum2_rank1[order[k]];
      _rank_merge_gnum1_to_gnum2_part1[k] = rank_merge_gnum1_to_gnum2_part1[order[k]];
      _rank_merge_gnum1_to_gnum2_lnum1[k] = rank_merge_gnum1_to_gnum2_lnum1[order[k]]; 
    }

    for (int k = 0; k < _max_part; k++) {
      idx_elt_part[k+1] = idx_elt_part[k] + n_elt_part[k]; 
    }

    for (int k1 = 0; k1 < _max_part; k1++) {

      int _n_elt_part = n_elt_part[k1];

      for (int j = 0; j < _n_elt_part; j++) {
        order[j] = j;
      }

      int *_part_rank_merge_gnum1_to_gnum2_part2 = _rank_merge_gnum1_to_gnum2_part2 + idx_elt_part[k1];
      int *_part_rank_merge_gnum1_to_gnum2_lnum2 = _rank_merge_gnum1_to_gnum2_lnum2 + idx_elt_part[k1];
      int *_part_rank_merge_gnum1_to_gnum2_rank1 = _rank_merge_gnum1_to_gnum2_rank1 + idx_elt_part[k1];
      int *_part_rank_merge_gnum1_to_gnum2_part1 = _rank_merge_gnum1_to_gnum2_part1 + idx_elt_part[k1];
      int *_part_rank_merge_gnum1_to_gnum2_lnum1 = _rank_merge_gnum1_to_gnum2_lnum1 + idx_elt_part[k1];

      int *part_rank_merge_gnum1_to_gnum2_part2 = rank_merge_gnum1_to_gnum2_part2 + idx_elt_part[k1];
      int *part_rank_merge_gnum1_to_gnum2_lnum2 = rank_merge_gnum1_to_gnum2_lnum2 + idx_elt_part[k1];
      int *part_rank_merge_gnum1_to_gnum2_rank1 = rank_merge_gnum1_to_gnum2_rank1 + idx_elt_part[k1];
      int *part_rank_merge_gnum1_to_gnum2_part1 = rank_merge_gnum1_to_gnum2_part1 + idx_elt_part[k1];
      int *part_rank_merge_gnum1_to_gnum2_lnum1 = rank_merge_gnum1_to_gnum2_lnum1 + idx_elt_part[k1];
      int *part_rank_i_send_buffer              = rank_i_send_buffer              + idx_elt_part[k1];  

      PDM_sort_int (_part_rank_merge_gnum1_to_gnum2_lnum2, order, _n_elt_part);

      int pre_val = -1;
      for (int k2 = 0; k2 < _n_elt_part; k2++) {
        part_rank_merge_gnum1_to_gnum2_part2[k2] = _part_rank_merge_gnum1_to_gnum2_part2[k2];
        part_rank_merge_gnum1_to_gnum2_lnum2[k2] = _part_rank_merge_gnum1_to_gnum2_lnum2[k2];
        part_rank_merge_gnum1_to_gnum2_rank1[k2] = _part_rank_merge_gnum1_to_gnum2_rank1[order[k2]];
        part_rank_merge_gnum1_to_gnum2_part1[k2] = _part_rank_merge_gnum1_to_gnum2_part1[order[k2]];
        part_rank_merge_gnum1_to_gnum2_lnum1[k2] = _part_rank_merge_gnum1_to_gnum2_lnum1[order[k2]]; 

        if (pre_val != part_rank_merge_gnum1_to_gnum2_lnum2[k2]) {
          cpt_buff++;   
          ptp->default_n_send_buffer[i]++;
        }
        part_rank_i_send_buffer[k2] = cpt_buff;

      }

    }

  }

  /* 3 - Define Default_n_send_buffer and  Default_i_send_buffer*/

  for (int i = 0; i < n_rank; i++) {
    ptp->default_i_send_buffer[i+1] = ptp->default_i_send_buffer[i] + ptp->default_n_send_buffer[i];
  }

  /* 4 - Define gnum1_to_send_buffer */

  ptp->gnum1_to_send_buffer = malloc (sizeof (int*) * n_part1);

  for (int i = 0; i < n_part1; i++) {
    ptp->gnum1_to_send_buffer[i] = malloc (sizeof (int) * n_elt1[i]);    
    for (int k = 0; k < n_elt1[i]; k++) {
      ptp->gnum1_to_send_buffer[i][k] = -1;
    }
  }

  for (int i = 0; i < n_total_elt; i++) {
    int ipart1 = rank_merge_gnum1_to_gnum2_part1[i];
    int ielt1  = rank_merge_gnum1_to_gnum2_lnum1[i];
    int idx    = i_send_buffer[i];
    assert (idx >! 0);
    ptp->gnum1_to_send_buffer[ipart1][ielt1] = idx;
  }

  for (int i = 0; i < n_part1; i++) {
    for (int k = 0; k < n_elt1[i]; k++) {
      assert (ptp->gnum1_to_send_buffer[i][k] != -1);
    }
  }

  free (n_elt_part);
  free (order);

  free (_merge_gnum1_to_gnum2_rank2);
  free (_merge_gnum1_to_gnum2_part2);
  free (_merge_gnum1_to_gnum2_lnum2);
  free (_merge_gnum1_to_gnum2_rank1);
  free (_merge_gnum1_to_gnum2_part1);
  free (_merge_gnum1_to_gnum2_lnum1);

  free (idx_gnum1_to_gnum2_rank);
  free (n_gnum1_to_gnum2_rank);

  free (i_send_buffer);

  /* 5 - Define Default_n_recv_buffer and  Default_i_recv_buffer*/

  


  /* 6 - Alltoall on gnum1_to_gnum2_part gnum1_to_gnum2_elt orig_gnum1 */

  /* 7 - gnum location on orig_gnum1 */

  /* 8 - Define ref_gnum2 and unref_gnum2 */

  /* 9 - Define gnum1_com_from */

  /* 10 - Define recv_buffer_to_ref_gnum2 */


  return ptp;

}


PDM_partgnum1_to_partgnum2_t *
PDM_partgnum1_to_partgnum2_create_cf
(
 const PDM_g_num_t    **gnum_elt1,
 const int            *n_elt1,
 const int             n_part1,
 const PDM_g_num_t    **gnum_elt2,
 const int            *n_elt2,
 const int             n_part2,
 const int           **gnum1_to_gnum2_idx,
 const PDM_g_num_t   **gnum1_to_gnum2,
 const PDM_MPI_Fint    fcomm
)
{
  const PDM_MPI_Comm _comm        = PDM_MPI_Comm_f2c(fcomm);
  return PDM_partgnum1_to_partgnum2_create (gnum_elt1, n_elt1, n_part1,
                                            gnum_elt2, n_elt2, n_part2,
                                            gnum1_to_gnum2_idx, gnum1_to_gnum2, _comm);
}


/**
 *
 * \brief Initialize an exchange
 *
 * \param [in]   ptp           Block to part structure
 * \param [in]   s_data        Data size
 * \param [in]   t_stride      Stride type
 * \param [in]   part1_stride  Partition 1 stride
 * \param [in]   part1_data    Partition 1 data
 * \param [out]  part2_stride  Partition 2 stride
 * \param [out]  part2_data    Partition 2 data
 *
 */

void
PDM_partgnum1_to_partgnum2_exch
(
 PDM_partgnum1_to_partgnum2_t *ptp,
 const size_t                  s_data,
 const PDM_stride_t            t_stride,
 int                         **part1_stride,
 void                        **part1_data,
 int                         **part2_stride,
 void                        **part2_data
)
{
  PDM_UNUSED (ptp);
  PDM_UNUSED (s_data);
  PDM_UNUSED (t_stride);
  PDM_UNUSED (part1_stride);
  PDM_UNUSED (part1_data);
  PDM_UNUSED (part2_stride);
  PDM_UNUSED (part2_data);

  PDM_error(__FILE__, __LINE__, 0,
            "Error PDM_partgnum1_to_partgnum2_exch not yet implemente\n");
}



/**
 *
 * \brief Initialize an exchange with allocation of result arrays
 *
 * \param [in]   ptp           Block to part structure
 * \param [in]   s_data        Data size
 * \param [in]   t_stride      Stride type
 * \param [in]   part1_stride  Partition 1 stride
 * \param [in]   part1_data    Partition 1 data
 * \param [out]  part2_stride  Partition 2 stride
 * \param [out]  part2_data    Partition 2 data
 *
 */

void
PDM_partgnum1_to_partgnum2_exch_with_alloc
(
 PDM_partgnum1_to_partgnum2_t *ptp,
 const size_t               s_data,
 const PDM_stride_t         t_stride,
 int                      **part1_stride,
 void                     **part1_data,
 int                     ***part2_stride,
 void                    ***part2_data
)
{
  PDM_UNUSED (ptp);
  PDM_UNUSED (s_data);
  PDM_UNUSED (t_stride);
  PDM_UNUSED (part1_stride);
  PDM_UNUSED (part1_data);
  PDM_UNUSED (part2_stride);
  PDM_UNUSED (part2_data);

  PDM_error(__FILE__, __LINE__, 0,
            "Error PDM_partgnum1_to_partgnum2_exch not yet implemente\n");
}


/**
 *
 * \brief Get referenced gnum2 elements
 *
 * \param [in]   ptp           Block to part structure
 * \param [out]  n_ref_gnum2   Number of referenced gnum2
 * \param [out]  ref_gnum2     Referenced gnum2
 *
 */

void
PDM_partgnum1_to_partgnum2_ref_gnum2_get
(
 PDM_partgnum1_to_partgnum2_t *ptp,
 int                         **n_ref_gnum2,
 int                        ***ref_gnum2
)
{
  PDM_UNUSED (ptp);
  PDM_UNUSED (n_ref_gnum2);
  PDM_UNUSED (ref_gnum2);

  PDM_error(__FILE__, __LINE__, 0,
            "Error PDM_partgnum1_to_partgnum2_ref_gnum2_get not yet implemente\n");

}


/**
 *
 * \brief Get unreferenced gnum2 elements
 *
 * \param [in]   ptp           Block to part structure
 * \param [out]  n_unref_gnum2   Number of referenced gnum2
 * \param [out]  unref_gnum2     Referenced gnum2
 *
 */

void
PDM_partgnum1_to_partgnum2_unref_gnum2_get
(
 PDM_partgnum1_to_partgnum2_t *ptp,
 int                         **n_unref_gnum2,
 int                        ***unref_gnum2
)
{
  PDM_UNUSED (ptp);
  PDM_UNUSED (n_unref_gnum2);
  PDM_UNUSED (unref_gnum2);

  PDM_error(__FILE__, __LINE__, 0,
            "Error PDM_partgnum1_to_partgnum2_unref_gnum2_get not yet implemente\n");
  
}


/**
 *
 * \brief Get gnum come from gnum1 for each referenced gnum2
 *
 * \param [in]   ptp                 Block to part structure
 * \param [out]  gnum1_come_from_idx Index for gnum1_come_from array (size = \ref n_part2)  
 * \param [out]  gnum1_come_from     Gnum come from gnum1 for each referenced gnum2
 *
 */

void
PDM_partgnum1_to_partgnum2_gnum1_come_from_get
(
 PDM_partgnum1_to_partgnum2_t *ptp,
 int                        ***gnum1_come_from_idx,
 PDM_g_num_t                ***gnum1_come_from
)
{
  PDM_UNUSED (ptp);
  PDM_UNUSED (gnum1_come_from_idx);
  PDM_UNUSED (gnum1_come_from);

  PDM_error(__FILE__, __LINE__, 0,
            "Error PDM_partgnum1_to_partgnum2_gnum1_come_from_get not yet implemente\n");
}


/**
 *
 * \brief Initialize a asynchronus issend
 *
 * \param [in]   ptp           Block to part structure
 * \param [in]   s_data        Data size
 * \param [in]   cst_stride    Constant stride
 * \param [in]   part1_data    Partition 1 data
 * \param [in]   tag           Tag of the exchange 
 * \param [out]  request       Request
 *
 */

void
PDM_partgnum1_to_partgnum2_issend
(
 PDM_partgnum1_to_partgnum2_t *ptp,
 const size_t                  s_data,
 const int                     cst_stride,
 void                        **part1_data,
 int                           tag,
 int                          *request
)
{
  PDM_UNUSED (ptp);
  PDM_UNUSED (s_data);
  PDM_UNUSED (cst_stride);
  PDM_UNUSED (part1_data);
  PDM_UNUSED (tag);
  PDM_UNUSED (request);

  PDM_error(__FILE__, __LINE__, 0,
            "Error PDM_partgnum1_to_partgnum2_issend not yet implemente\n");
  
}


/**
 *
 * \brief Wait a asynchronus issend
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  tag           Tag of the exchange 
 * \param [in]  request       Request
 *
 */

void
PDM_partgnum1_to_partgnum2_issend_wait
(
 PDM_partgnum1_to_partgnum2_t *ptp,
 int                           tag,
 int                        request
)
{
  PDM_UNUSED (ptp);
  PDM_UNUSED (tag);
  PDM_UNUSED (request);

  PDM_error(__FILE__, __LINE__, 0,
            "Error PDM_partgnum1_to_partgnum2_issend_wait not yet implemente\n");
  
}


/**
 *
 * \brief Initialize a asynchronus irecv
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  s_data        Data size
 * \param [in]  cst_stride    Constant stride
 * \param [in]  part1_data    Partition 2 data
 * \param [in]  tag           Tag of the exchange 
 * \param [out] request       Request
 *
 */

void
PDM_partgnum1_to_partgnum2_irecv
(
 PDM_partgnum1_to_partgnum2_t *ptp,
 const size_t                  s_data,
 const int                     cst_stride,
 void                        **part2_data,
 int                           tag,
 int                          *request
)
{
  PDM_UNUSED (ptp);
  PDM_UNUSED (s_data);
  PDM_UNUSED (cst_stride);
  PDM_UNUSED (part2_data);
  PDM_UNUSED (tag);
  PDM_UNUSED (request);

  PDM_error(__FILE__, __LINE__, 0,
            "Error PDM_partgnum1_to_partgnum2_irecv not yet implemente\n");
  
}


/**
 *
 * \brief Initialize a asynchronus irecv
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  tag           Tag of the exchange 
 * \param [in]  request       Request
 *
 */

void
PDM_partgnum1_to_partgnum2_irecv_wait
(
 PDM_partgnum1_to_partgnum2_t *ptp,
 int                           tag,
 int                           request
)
{
  PDM_UNUSED (ptp);
  PDM_UNUSED (tag);
  PDM_UNUSED (request);

  PDM_error(__FILE__, __LINE__, 0,
            "Error PDM_partgnum1_to_partgnum2_irecv_wait not yet implemente\n");
  
}


/**
 *
 * \brief Free a block to part structure
 *
 * \param [inout] ptp  Block to part structure
 *
 * \return       NULL
 */

PDM_partgnum1_to_partgnum2_t *
PDM_partgnum1_to_partgnum2_free
(
 PDM_partgnum1_to_partgnum2_t *ptp
)
{
  if (ptp == NULL) {
    return NULL;
  }

  if (ptp->gnum1_to_send_buffer != NULL) {
    for (int i = 0; i < ptp->n_part1; i++) {    
      free (ptp->gnum1_to_send_buffer[i]);
    }
    free (ptp->gnum1_to_send_buffer);
  }

  if (ptp->recv_buffer_to_ref_gnum2 != NULL) {  
    for (int i = 0; i < ptp->n_part2; i++) {
      free (ptp->ref_gnum2[i]);    
      free (ptp->unref_gnum2[i]);    
      free (ptp->gnum1_come_from_idx[i]);    
      free (ptp->gnum1_come_from[i]);    
      free (ptp->recv_buffer_to_ref_gnum2[i]);
    }
    free (ptp->recv_buffer_to_ref_gnum2);
    free (ptp->ref_gnum2);    
    free (ptp->unref_gnum2);    
    free (ptp->n_ref_gnum2);    
    free (ptp->n_unref_gnum2);    
    free (ptp->gnum1_come_from_idx);    
    free (ptp->gnum1_come_from);    
  }  
  
  if (ptp->async_l_array != 0) {
    for (int i = 0; i < ptp->async_l_array; i++) {    
      free (ptp->async_send_buffer[i]);      
      free (ptp->async_recv_buffer[i]);      
      free (ptp->async_n_send_buffer[i]);  
      free (ptp->async_i_send_buffer[i]);  
      free (ptp->async_n_recv_buffer[i]);  
      free (ptp->async_i_recv_buffer[i]);  
    }
    free (ptp->async_s_data); 
    free (ptp->async_cst_stride);        
    free (ptp->async_send_request);    
    free (ptp->async_recv_request);    
    free (ptp->async_send_buffer);      
    free (ptp->async_recv_buffer);      
    free (ptp->async_n_send_buffer);  
    free (ptp->async_i_send_buffer);  
    free (ptp->async_n_recv_buffer);  
    free (ptp->async_i_recv_buffer);  
  }

  free (ptp->default_n_send_buffer);
  free (ptp->default_i_send_buffer);
  free (ptp->default_n_recv_buffer);
  free (ptp->default_i_recv_buffer);

  ptp->async_n_exch            = 0;           
  ptp->async_l_array           = 0;          

  free(ptp);
  return NULL;
  
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
